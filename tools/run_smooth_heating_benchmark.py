"""Export a constant-property linear-ramp magnetic-diffusion reference.

All inputs are explicit SI values. The default fixture is manufactured. No
solver is run; finite-volume field references, infinite-half-space energy and
finite wall/point adiabatic heating are exported for later host verification.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
import math
from numbers import Real
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.heating_comparison import HeatingProfile, load_heating_profile_csv, write_heating_profile_csv
from mhdlab.smooth_heating_benchmark import (
    KERNEL_LOG_TOLERANCE, QUADRATURE_RELATIVE_TOLERANCE, LinearRampConductivityHalfSpace, _log_kernels,
)


CONFIG_SCHEMA = "mhdlab-smooth-heating-benchmark-config-v1"
CODE_PATHS = (
    "mhdlab/smooth_heating_benchmark.py", "mhdlab/heating_benchmarks.py",
    "mhdlab/heating_comparison.py", "tools/run_smooth_heating_benchmark.py",
)
ENERGY_UNITS = {
    "time_s": "s", "magnetic_energy_j_m2": "J m^-2", "joule_heat_j_m2": "J m^-2",
    "poynting_input_j_m2": "J m^-2", "surface_power_w_m2": "W m^-2", "energy_residual_j_m2": "J m^-2",
}
TEMPERATURE_UNITS = {"time_s": "s", "point_depth_m": "m", "point_kind": "1", "adiabatic_temperature_rise_k": "K"}


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate configuration key: {key}")
        result[key] = value
    return result


def _nonfinite(value):
    raise ValueError(f"nonfinite configuration constant: {value}")


def load_configuration(path):
    raw = path.read_bytes()
    config = json.loads(raw, object_pairs_hook=_pairs, parse_constant=_nonfinite)
    keys = {"schema", "parameter_status", "parameter_provenance", "benchmark", "depth_stop_m", "depth_cells", "times_s"}
    if not isinstance(config, dict) or set(config) != keys or config["schema"] != CONFIG_SCHEMA:
        raise ValueError(f"configuration must contain exactly the {CONFIG_SCHEMA} fields")
    if config["parameter_status"] not in {"manufactured", "sourced", "unvalidated"}:
        raise ValueError("parameter_status must be manufactured, sourced or unvalidated")
    if not isinstance(config["parameter_provenance"], str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must describe the supplied constants and drive")
    count = config["depth_cells"]
    if isinstance(count, bool) or not isinstance(count, int) or not 1 <= count <= 100_000:
        raise ValueError("depth_cells must be an integer in [1,100000]")
    stop = config["depth_stop_m"]
    if isinstance(stop, bool) or not isinstance(stop, Real) or not math.isfinite(float(stop)) or stop <= 0:
        raise ValueError("depth_stop_m must be a finite positive SI length")
    times = np.asarray(config["times_s"], dtype=object)
    if (times.ndim != 1 or times.size == 0 or times.size * count > 1_000_000
            or any(isinstance(value, bool) or not isinstance(value, Real)
                   or not math.isfinite(float(value)) or value <= 0 for value in times)):
        raise ValueError("times_s must be positive finite times within the reference memory limit")
    times = np.asarray(times, dtype=float)
    if np.any(np.diff(times) <= 0):
        raise ValueError("times_s must be strictly increasing; canonical field CSV times are positive")
    if not isinstance(config["benchmark"], dict):
        raise ValueError("benchmark must contain all explicit constant-property linear-ramp parameters")
    benchmark = LinearRampConductivityHalfSpace(**config["benchmark"])
    edges = np.linspace(0., float(stop), count + 1)
    if np.any(np.diff(edges) <= 0):
        raise ValueError("depth cells are not representable at the supplied SI length")
    return config, benchmark, edges, times, hashlib.sha256(raw).hexdigest(), raw


def save_plot(output, edges, times, references, budgets):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 2, figsize=(11, 8))
    centers = (edges[:-1] + .5 * np.diff(edges)) * 1e3
    for time, values in zip(times, references):
        label = f"{time * 1e6:g} µs"
        for axis, field in ((axes[0, 0], "magnetic_field_t"), (axes[0, 1], "current_density_a_m2"),
                            (axes[1, 0], "joule_power_w_m3")):
            axis.plot(centers, getattr(values, field), label=label)
    for axis, label in ((axes[0, 0], "Cell-average Bx [T]"), (axes[0, 1], "Cell-average Jz [A/m²]"),
                        (axes[1, 0], "Cell-average Joule power [W/m³]")):
        axis.set_xlabel("Inward material depth [mm]")
        axis.set_ylabel(label)
        axis.legend()
    for key, label in (("magnetic_energy_j_m2", "Magnetic storage"), ("joule_heat_j_m2", "Joule heat"),
                       ("poynting_input_j_m2", "Incoming Poynting energy")):
        axes[1, 1].plot(times * 1e6, [getattr(value, key) for value in budgets], "o-", label=label)
    axes[1, 1].set_xlabel("Time [µs]")
    axes[1, 1].set_ylabel("Infinite-half-space energy [J/m²]")
    axes[1, 1].legend()
    figure.suptitle("Linear-ramp magnetic-diffusion reference")
    figure.tight_layout()
    figure.savefig(output, dpi=150)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/smooth_heating_manufactured.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        source_hashes = {name: _sha(ROOT / name) for name in CODE_PATHS}
        config, benchmark, edges, times, config_hash, config_payload = load_configuration(args.config)
        references = [benchmark.cell_averages(edges, float(time)) for time in times]
        budgets = [benchmark.energy_budget(float(time)) for time in times]
        depths = np.concatenate(([0.], edges[:-1] + .5 * np.diff(edges)))
        temperatures = [benchmark.adiabatic_temperature_rise(depths, float(time)) for time in times]
        wall_temperature = [float(values[0]) for values in temperatures]
        outer_fraction = []
        for time in times:
            if benchmark.boundary_field_rate_t_s == 0:
                outer_fraction.append(None)
            else:
                # Evaluate the dimensionless b(u) directly. Multiplying rate*t
                # or dividing two physical fields can lose an otherwise finite
                # fraction through overflow/underflow of their dimensional scale.
                scale, _ = benchmark._scale(float(time))
                with np.errstate(over="ignore", under="ignore"):
                    coordinate = float(.5 * np.divide(float(edges[-1]), scale))
                outer_fraction.append(math.exp(_log_kernels(coordinate)[1]))
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".smooth-heating-", dir=args.output_dir.parent) as staging:
            output = Path(staging) / "result"
            output.mkdir()
            profiles = output / "profiles"
            profiles.mkdir()
            (output / "input_config.json").write_bytes(config_payload)
            exported = []
            for index, (time, values) in enumerate(zip(times, references)):
                profile = HeatingProfile(float(time), edges, values.magnetic_field_t, values.current_density_a_m2,
                                         values.joule_power_w_m3, "synthetic", "analytic-linear-ramp-reference",
                                         f"linear-ramp B_s(t)=rate*t; {config['parameter_provenance']}; config-sha256={config_hash}")
                csv_path, metadata_path = write_heating_profile_csv(profiles / f"time_{index:03d}.csv", profile)
                restored = load_heating_profile_csv(csv_path, metadata_path=metadata_path)
                for name in ("depth_edges_m", "magnetic_field_t", "current_density_a_m2", "joule_power_w_m3"):
                    if not np.array_equal(getattr(restored, name), getattr(profile, name)):
                        raise ValueError("canonical heating profile did not round-trip exactly")
                exported.append({"time_s": float(time), "csv": str(csv_path.relative_to(output)),
                                 "metadata": str(metadata_path.relative_to(output)),
                                 "csv_sha256": _sha(csv_path), "metadata_sha256": _sha(metadata_path)})
            with (output / "energy_budget.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(ENERGY_UNITS), lineterminator="\n")
                writer.writeheader()
                writer.writerows({"time_s": float(time), **asdict(budget)} for time, budget in zip(times, budgets))
            with (output / "adiabatic_point_temperature.csv").open("w", newline="") as stream:
                writer = csv.writer(stream, lineterminator="\n")
                writer.writerow(TEMPERATURE_UNITS)
                for time, values in zip(times, temperatures):
                    writer.writerows((float(time), float(depth), "wall" if index == 0 else "cell_center", float(value))
                                     for index, (depth, value) in enumerate(zip(depths, values)))
            if args.plot:
                save_plot(output / "overview.png", edges, times, references, budgets)
            if _sha(args.config) != config_hash or any(_sha(ROOT / name) != digest for name, digest in source_hashes.items()):
                raise ValueError("an input or source file changed during reference generation")
            summary = {
                "schema": "mhdlab-smooth-heating-benchmark-run-v1", "model": "constant_property_linear_ramp_half_space",
                "is_flash_result": False, "is_measured_shot": False,
                "parameter_status": config["parameter_status"], "parameter_provenance": config["parameter_provenance"],
                "parameters": asdict(benchmark), "magnetic_diffusivity_m2_s": benchmark.magnetic_diffusivity_m2_s,
                "drive": "B_s(t)=boundary_field_rate_t_s*t for t>=0; ramp continues without a stop",
                "config_sha256": config_hash, "source_sha256": source_hashes,
                "depth_cells": len(edges) - 1, "depth_stop_m": float(edges[-1]), "times_s": times.tolist(),
                "profiles": exported, "energy_budget_units": ENERGY_UNITS, "temperature_units": TEMPERATURE_UNITS,
                "energy_budgets": [asdict(value) for value in budgets],
                "wall_adiabatic_temperature_rise_k": wall_temperature,
                "finite_depth_boundary_B_fraction": outer_fraction,
                "numerical_qualifications": {"kernel_log_tolerance": KERNEL_LOG_TOLERANCE,
                                             "quadrature_relative_tolerance": QUADRATURE_RELATIVE_TOLERANCE},
                "output_sha256": {str(path.relative_to(output)): _sha(path) for path in sorted(output.rglob("*")) if path.is_file()},
                "limits": ["Analytic manufactured reference generation and exact CSV readback do not verify a FLASH run.",
                           "Stationary constant-property conductor; no flow, thermal conduction, Hall term, displacement current, radiation, phase change or EOS closure.",
                           "Each finite-volume Joule mean integrates J squared; squaring a cell-average current gives a different quantity.",
                           "Energy budgets cover the infinite half-space; finite-depth host boundaries need separate placement and refinement checks.",
                           "Adiabatic temperatures are local point rises, including the finite wall; no initial absolute material temperature is supplied.",
                           "The indefinitely continued linear ramp is a verification drive, with no measured pulse or foil-material calibration.",
                           "Canonical CSV times are strictly positive; the underlying reference also has a zero initial state at t=0."]}
            (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during generation; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as error:
        parser.error(str(error))
    print(json.dumps({"output_dir": str(args.output_dir), "profiles": len(times), "parameter_status": config["parameter_status"],
                      "maximum_energy_residual_j_m2": max(abs(value.energy_residual_j_m2) for value in budgets)}, allow_nan=False))


if __name__ == "__main__":
    main()
