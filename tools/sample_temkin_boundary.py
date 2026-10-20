#!/usr/bin/env python3
"""Export paired thermal boundary packets with Temkin-specific emission times.

Only emission is represented. No gap, collision, radiation or field solve runs.
Exact state masses and geometry are supplied separately from the source law.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping
from dataclasses import fields
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.temkin_boundary import sample_temkin_emission_packets
from mhdlab.temkin_io import SCHEMA as SOURCE_SCHEMA, load_temkin_history
from mhdlab.temkin_emission_time import TIME_IDENTITY
from tools.sample_surface_boundary import load_mass_config, _open_uniforms, _ARRAY_UNITS


CODE_PATHS = (
    "mhdlab/temkin_boundary.py", "mhdlab/temkin_emission_time.py", "mhdlab/temkin_release.py",
    "mhdlab/temkin_history.py", "mhdlab/temkin_io.py", "mhdlab/thermal_release.py",
    "mhdlab/thermal_boundary.py", "mhdlab/flash_surface.py", "mhdlab/hdf_integrity.py",
    "mhdlab/constants.py", "mhdlab/surface_sources.py", "tools/sample_surface_boundary.py",
    "tools/sample_temkin_boundary.py",
)


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _plain(value):
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_packets(output, packets, *, source_path, mass_path, source_hash, mass_hash, config, seed, redraws, generator,
                  source_run_record=None, additional_code_paths=()):
    arrays = {field.name: getattr(packets, field.name) for field in fields(packets)
              if isinstance(getattr(packets, field.name), np.ndarray)}
    units = dict(_ARRAY_UNITS)
    units.update({"resolved_position_m": "m", "uniform_width_position_m": "m",
                  "source_current_offset_m": "m", "x_edges_m": "m", "z_edges_m": "m",
                  "source_current_index": "1", "component_current_index": "1",
                  "minimum_packets_per_current_row": "packets", "packets_per_current_row": "packets",
                  "released_particles_per_current_row": "physical particles"})
    if set(arrays) - set(units):
        raise ValueError(f"packet units are unspecified for {sorted(set(arrays) - set(units))}")
    metadata = {field.name: getattr(packets, field.name) for field in fields(packets)
                if field.name not in {"state_ids", "reservoir_names", "surface_attributes"}
                and not isinstance(getattr(packets, field.name), np.ndarray)}
    if packets.time_model != TIME_IDENTITY:
        raise ValueError("packets must identify the qualified Temkin conditional-time model")
    metadata["source_parameter_status"] = None if source_run_record is None else source_run_record["parameter_status"]
    payload = dict(arrays)
    payload.update({"state_ids": np.asarray(packets.state_ids, dtype=str),
                    "reservoir_names": np.asarray(packets.reservoir_names, dtype=str),
                    "metadata_json": np.asarray(json.dumps(_plain(metadata), sort_keys=True, allow_nan=False)),
                    "surface_attributes_json": np.asarray(json.dumps(_plain(packets.surface_attributes), sort_keys=True, allow_nan=False))})
    np.savez_compressed(output / "boundary_packets.npz", **payload)
    residuals = []
    for component, (pool, width, count) in enumerate(zip(packets.component_reservoir_index, packets.component_width_index,
                                                        packets.component_released_particles)):
        selected = (packets.reservoir_index == pool) & (packets.source_width_index == width)
        if hasattr(packets, "source_current_index"):
            selected &= packets.source_current_index == packets.component_current_index[component]
        residuals.append(math.fsum((-float(count), *map(float, packets.particle_weight[selected]))))
    expected = math.fsum(map(float, packets.component_expected_translational_energy_j))
    sampled = math.fsum(map(float, packets.component_sampled_translational_energy_j))
    summary = {
        **metadata, "source_schema": packets.source_history_schema,
        "scope": "upstream weighted thermal emission with law-specific event times and paired width positions",
        "surface_attributes": _plain(packets.surface_attributes),
        "input_source": {"path": str(source_path.resolve()), "sha256": source_hash},
        "input_masses": {"path": str(mass_path.resolve()), "sha256": mass_hash},
        "source_run_record": source_run_record,
        "mass_configuration": config,
        "sampling": {"seed": seed, "bit_generator": generator, "numpy_version": np.__version__,
                     "zero_redraw_count": redraws, "uniforms": "strictly inside (0,1)"},
        "packets": len(packets.particle_weight),
        "component_particle_balance_residual": residuals,
        "total_released_particles": packets.total_released_particles,
        "total_packet_particles": math.fsum(map(float, packets.particle_weight)),
        "thermal_energy": {"expected_translational_energy_j": expected,
                           "sampled_translational_energy_j": sampled,
                           "sampling_residual_j": math.fsum((sampled, -expected)),
                           "relative_sampling_residual": (sampled - expected) / expected if expected else None},
        "array_units": {key: units[key] for key in arrays},
        "packet_file_sha256": _sha(output / "boundary_packets.npz"),
        "source_code_sha256": {name: _sha(ROOT / name) for name in dict.fromkeys((*CODE_PATHS, *additional_code_paths))},
        "limits": ["Heating source_kind and mass labels are retained. Source parameter_status requires its hash-associated run summary; absent that record it remains unreported.",
                   "Exact-state populations and masses are supplied; no internal-state chemistry is inferred.",
                   "Downstream modules must support unequal packet weights or provide a verified adapter.",
                   "Resolved and uniform hypotheses share time, state, velocity, mass and weight; choose only their position arrays.",
                   "A coupled wall must balance actual weighted launch energy; sampling error is not an extra heat source.",
                   "Prescribed surface temperatures do not receive automatic wall-energy feedback.",
                   "No gap, collision or radiative evolution is represented."]}
    (output / "summary.json").write_text(json.dumps(_plain(summary), indent=2, allow_nan=False) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--masses", type=Path, required=True)
    parser.add_argument("--source-run-summary", type=Path,
                        help="Optional source-run summary whose output hash binds parameter_status to this exact source file")
    parser.add_argument("--interval", type=int, required=True)
    parser.add_argument("--packets", type=int, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, required=True, help="A new directory; existing paths are refused")
    args = parser.parse_args()
    if args.seed < 0 or args.packets < 1:
        parser.error("--seed must be nonnegative and --packets positive")
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        source_hash, mass_hash = _sha(args.source), _sha(args.masses)
        config = load_mass_config(args.masses)
        history = load_temkin_history(args.source)
        source_run_record = None
        if args.source_run_summary is not None:
            summary_hash = _sha(args.source_run_summary)
            source_summary = json.loads(args.source_run_summary.read_bytes())
            if (not isinstance(source_summary, dict)
                    or not isinstance(source_summary.get("output_sha256"), dict)
                    or source_summary.get("schema") != "mhdlab-temkin-run-v1"
                    or source_summary.get("output_sha256", {}).get("temkin_release.h5") != source_hash
                    or source_summary.get("parameter_status") not in {"manufactured", "sourced", "unvalidated"}
                    or source_summary.get("parameter_provenance") != history.parameter_provenance):
                raise ValueError("source-run summary must bind its parameter status to this exact Temkin source hash and provenance")
            source_run_record = {"path": str(args.source_run_summary.resolve()), "sha256": summary_hash,
                                 "parameter_status": source_summary["parameter_status"],
                                 "parameter_provenance": source_summary["parameter_provenance"]}
        draws, redraws, generator = _open_uniforms(args.seed, args.packets)
        packets = sample_temkin_emission_packets(history, config["state_masses_kg"], args.interval, draws,
                                                 outward_normal_sign=config["outward_normal_sign"],
                                                 current_length_m=config["current_direction_length_m"])
        if _sha(args.source) != source_hash or _sha(args.masses) != mass_hash:
            raise ValueError("an input file changed during sampling")
        if source_run_record is not None and _sha(args.source_run_summary) != source_run_record["sha256"]:
            raise ValueError("source-run summary changed during sampling")
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".temkin-boundary-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            summary = write_packets(output, packets, source_path=args.source, mass_path=args.masses,
                                    source_hash=source_hash, mass_hash=mass_hash, config=config, seed=args.seed,
                                    redraws=redraws, generator=generator, source_run_record=source_run_record)
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during sampling; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_dir": str(args.output_dir), "packets": summary["packets"],
                      "source_components": packets.required_minimum_n,
                      "total_packet_particles": summary["total_packet_particles"]}, allow_nan=False))


if __name__ == "__main__":
    main()
