#!/usr/bin/env python3
"""Sample a validated surface-release file into paired thermal boundary packets.

This is an upstream boundary artifact, not a gap, Boltzmann or radiation solve.
State masses and geometry conversion are explicit configuration inputs. Example:

  python tools/run_surface_release.py --demo --output-dir runs/release_demo
  python tools/sample_surface_boundary.py \
      --source runs/release_demo/surface_release.h5 \
      --masses configs/manufactured_boundary_masses.json --interval 0 \
      --packets 1000 --seed 42 --output-dir runs/boundary_demo

The example mass and heating parameters are manufactured verification inputs;
they do not constitute a water mechanism or a physical material-property fit.
"""

from __future__ import annotations

import argparse
from dataclasses import fields
import hashlib
import json
import math
from numbers import Real
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.surface_release_io import SCHEMA as SOURCE_SCHEMA, load_release_history
from mhdlab.thermal_boundary import SCHEMA, sample_emission_packets


MASS_CONFIG_SCHEMA = "mhdlab-thermal-boundary-config-v1"
_CODE_PATHS = (
    "mhdlab/thermal_boundary.py", "mhdlab/thermal_release.py",
    "mhdlab/surface_release_history.py", "mhdlab/surface_release_io.py",
    "tools/sample_surface_boundary.py",
)
_ARRAY_UNITS = {
    "time_offset_s": "s", "time_s": "s", "resolved_x_m": "m", "uniform_x_m": "m",
    "reservoir_index": "1", "source_width_index": "1", "temperature_k": "K",
    "mass_kg": "kg", "velocity_m_s": "m s^-1", "momentum_kg_m_s": "kg m s^-1",
    "particle_weight": "physical particles per packet", "internal_energy_ev": "eV per particle",
    "desorption_energy_ev": "eV per particle", "component_reservoir_index": "1",
    "component_width_index": "1", "component_released_particles": "physical particles",
    "component_probabilities": "1", "component_expected_translational_energy_j": "J",
    "component_sampled_translational_energy_j": "J", "component_translational_energy_sampling_residual_j": "J",
    "packet_counts": "packets per source component",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _positive_real(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a positive finite real number")
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be a positive finite real number")
    return result


def load_mass_config(path: Path) -> dict:
    """Read only explicit state masses, status/provenance and face conversion."""
    config = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("schema") != MASS_CONFIG_SCHEMA:
        raise ValueError(f"mass configuration must declare schema {MASS_CONFIG_SCHEMA!r}")
    if config.get("mass_status") not in {"manufactured", "sourced", "unvalidated"}:
        raise ValueError("mass_status must be manufactured, sourced, or unvalidated")
    if not isinstance(config.get("mass_provenance"), str) or not config["mass_provenance"].strip():
        raise ValueError("mass_provenance must describe the supplied exact-state masses")
    masses = config.get("state_masses_kg")
    if not isinstance(masses, dict) or not masses:
        raise ValueError("state_masses_kg must be a nonempty map of exact state IDs to masses")
    validated = {}
    for state, mass in masses.items():
        if not isinstance(state, str) or not state or state != state.strip():
            raise ValueError("mass-map keys must be exact nonempty state IDs without surrounding whitespace")
        validated[state] = _positive_real(mass, f"mass for {state!r}")
    length = _positive_real(config.get("current_direction_length_m"), "current_direction_length_m")
    normal = config.get("outward_normal_sign")
    if isinstance(normal, bool) or not isinstance(normal, int) or normal not in (-1, 1):
        raise ValueError("outward_normal_sign must be explicitly -1 or +1")
    return {"schema": MASS_CONFIG_SCHEMA, "mass_status": config["mass_status"],
            "mass_provenance": config["mass_provenance"], "state_masses_kg": validated,
            "current_direction_length_m": length, "outward_normal_sign": normal}


def _json_value(value):
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def _open_uniforms(seed: int, number: int) -> tuple[np.ndarray, int, str]:
    rng = np.random.default_rng(seed)
    draws = rng.random((number, 6))
    redraws = 0
    zero = draws == 0
    while np.any(zero):
        count = int(np.count_nonzero(zero))
        redraws += count
        draws[zero] = rng.random(count)
        zero = draws == 0
    return draws, redraws, type(rng.bit_generator).__name__


def write_packets(output: Path, packets, *, source_path: Path, masses_path: Path,
                  config: dict, seed: int, redraws: int, bit_generator: str) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    array_payload = {field.name: getattr(packets, field.name) for field in fields(packets)
                     if isinstance(getattr(packets, field.name), np.ndarray)}
    # Explicit Unicode arrays allow a downstream reader to use allow_pickle=False.
    array_payload.update({
        "state_ids": np.asarray(packets.state_ids, dtype=str),
        "reservoir_names": np.asarray(packets.reservoir_names, dtype=str),
        "schema": np.asarray(packets.schema),
        "time_start_s": np.asarray(packets.time_start_s),
        "time_end_s": np.asarray(packets.time_end_s),
        "interval_index": np.asarray(packets.interval_index, dtype=np.int64),
        "current_length_m": np.asarray(packets.current_length_m),
        "outward_normal_sign": np.asarray(packets.outward_normal_sign, dtype=np.int64),
        "surface_attributes_json": np.asarray(json.dumps(dict(packets.surface_attributes), sort_keys=True, allow_nan=False)),
        "parameter_provenance": np.asarray(packets.parameter_provenance),
    })
    archive = output / "boundary_packets.npz"
    np.savez_compressed(archive, **array_payload)
    component_residuals = []
    for component, (pool, width) in enumerate(zip(packets.component_reservoir_index, packets.component_width_index)):
        selected = (packets.reservoir_index == pool) & (packets.source_width_index == width)
        emitted = math.fsum(float(weight) for weight in packets.particle_weight[selected])
        component_residuals.append(math.fsum((emitted, -float(packets.component_released_particles[component]))))
    expected_energy = math.fsum(float(value) for value in packets.component_expected_translational_energy_j)
    sampled_energy = math.fsum(float(value) for value in packets.component_sampled_translational_energy_j)
    energy_residual = math.fsum((sampled_energy, -expected_energy))
    summary = {
        "schema": SCHEMA,
        "scope": "upstream thermal emission packets with paired resolved/uniform width positions; no downstream solver invoked",
        "source_schema": SOURCE_SCHEMA,
        "surface_attributes": dict(packets.surface_attributes),
        "parameter_provenance": packets.parameter_provenance,
        "mass_config": config,
        "input_files": {
            "source": str(source_path.resolve()), "source_sha256": _sha256(source_path),
            "masses": str(masses_path.resolve()), "masses_sha256": _sha256(masses_path),
        },
        "sampling": {"seed": seed, "bit_generator": bit_generator, "numpy_version": np.__version__,
                     "uniforms": "strictly open (0,1); endpoint zero is redrawn", "zero_redraw_count": redraws,
                     "interval_index": packets.interval_index, "packets": len(packets.particle_weight),
                     "required_minimum_n": packets.required_minimum_n,
                     "allocation_model": packets.allocation_model, "time_model": packets.time_model,
                     "velocity_model": packets.velocity_model,
                     "current_direction_length_m": packets.current_length_m,
                     "outward_normal_sign": packets.outward_normal_sign,
                     "unequal_macro_weights": packets.unequal_macro_weights},
        "time_start_s": packets.time_start_s, "time_end_s": packets.time_end_s,
        "state_ids": sorted(set(packets.state_ids)),
        "component_reservoir_index": packets.component_reservoir_index,
        "component_width_index": packets.component_width_index,
        "component_released_particles": packets.component_released_particles,
        "component_probabilities": packets.component_probabilities,
        "packet_counts": packets.packet_counts,
        "component_particle_balance_residual": component_residuals,
        "total_released_particles": packets.total_released_particles,
        "total_packet_particles": math.fsum(float(weight) for weight in packets.particle_weight),
        "thermal_energy": {
            "expected_translational_energy_j": expected_energy,
            "sampled_translational_energy_j": sampled_energy,
            "sampling_residual_j": energy_residual,
            "relative_sampling_residual": energy_residual / expected_energy if expected_energy > 0 else None,
            "component_expected_translational_energy_j": packets.component_expected_translational_energy_j,
            "component_sampled_translational_energy_j": packets.component_sampled_translational_energy_j,
            "component_sampling_residual_j": packets.component_translational_energy_sampling_residual_j,
            "interpretation": "finite velocity samples estimate the 2*kB*T source moment; residual is sampling error, not extra physical wall heating",
        },
        "comparison": "Resolved and uniform hypotheses use these identical time/state/momentum/weight packets; select resolved_x_m or uniform_x_m without resampling",
        "limits": [
            "Unequal packet weights require downstream support; no adapter to an equal-weight collision engine is supplied.",
            "Absolute float64 time_s can lose local offsets; retain time_start_s and time_offset_s separately.",
            "Temperature mixtures, specified internal states and component inventories are retained; no water mechanism or spectrum is inferred.",
        ],
        "array_units": _ARRAY_UNITS,
        "packet_file": archive.name, "packet_file_sha256": _sha256(archive),
        "source_code_sha256": {name: _sha256(ROOT / name) for name in _CODE_PATHS},
    }
    (output / "summary.json").write_text(json.dumps(_json_value(summary), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, required=True, help="Validated surface_release.h5")
    parser.add_argument("--masses", type=Path, required=True, help="Explicit exact-state mass/geometry JSON configuration")
    parser.add_argument("--interval", type=int, required=True)
    parser.add_argument("--packets", type=int, required=True, help="At least one packet per positive source component")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--current-length-m", type=float, help="Explicit override of config current-direction length")
    parser.add_argument("--normal-sign", type=int, choices=(-1, 1), help="Explicit override of config outward normal")
    parser.add_argument("--output-dir", type=Path, default=Path("runs/thermal_boundary_packets"))
    args = parser.parse_args()
    if args.packets < 1 or args.seed < 0:
        parser.error("--packets must be positive and --seed nonnegative")
    try:
        config = load_mass_config(args.masses)
        if args.current_length_m is not None:
            config["current_direction_length_m"] = _positive_real(args.current_length_m, "--current-length-m")
        if args.normal_sign is not None:
            config["outward_normal_sign"] = args.normal_sign
        history = load_release_history(args.source)
        uniform_draws, redraws, generator = _open_uniforms(args.seed, args.packets)
        packets = sample_emission_packets(
            history, config["state_masses_kg"], args.interval, uniform_draws,
            outward_normal_sign=config["outward_normal_sign"],
            current_length_m=config["current_direction_length_m"],
        )
        write_packets(args.output_dir, packets, source_path=args.source, masses_path=args.masses,
                      config=config, seed=args.seed, redraws=redraws, bit_generator=generator)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_dir": str(args.output_dir), "packets": len(packets.particle_weight),
                      "components": packets.required_minimum_n, "mass_status": config["mass_status"],
                      "total_released_particles": packets.total_released_particles}))


if __name__ == "__main__":
    main()
