"""Export exact periodic resistive-induction and conductive-heating references.

All coefficients are explicit SI inputs.  The default configuration is
manufactured, with no foil-material or measured-shot calibration.  This
command generates a reference; it does not execute FLASH or another solver.
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

from mhdlab.conductive_heating_benchmark import (
    ASSUMPTIONS, MAX_ABSOLUTE_PHASE_RAD, PeriodicConductiveHeating,
)


CONFIG_SCHEMA = "mhdlab-periodic-conductive-heating-config-v1"
RUN_SCHEMA = "mhdlab-periodic-conductive-heating-run-v1"
CODE_PATHS = ("mhdlab/conductive_heating_benchmark.py", "tools/run_conductive_heating_benchmark.py")
CELL_UNITS = {
    "time_s": "s", "cell_index": "1", "left_position_m": "m", "right_position_m": "m",
    "magnetic_field_t": "T", "current_density_a_m2": "A m^-2",
    "joule_power_w_m3": "W m^-3", "temperature_rise_k": "K",
}
ENERGY_UNITS = {
    "time_s": "s", "initial_mean_magnetic_energy_j_m3": "J m^-3",
    "mean_magnetic_energy_j_m3": "J m^-3", "mean_thermal_energy_j_m3": "J m^-3",
    "mean_total_energy_j_m3": "J m^-3", "mean_joule_power_w_m3": "W m^-3",
    "energy_residual_j_m3": "J m^-3",
}


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
    raise ValueError(f"nonfinite JSON constant: {value}")


def _finite_real(value, name):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(float(value)):
        raise ValueError(f"{name} must be a finite real SI value")
    return float(value)


def load_configuration(path):
    payload = path.read_bytes()
    config = json.loads(payload, object_pairs_hook=_pairs, parse_constant=_nonfinite)
    keys = {"schema", "parameter_status", "parameter_provenance", "benchmark",
            "position_origin_m", "periods", "cells_per_period", "times_s"}
    if not isinstance(config, dict) or set(config) != keys or config["schema"] != CONFIG_SCHEMA:
        raise ValueError(f"configuration must contain exactly the {CONFIG_SCHEMA} fields")
    if not isinstance(config["parameter_status"], str) or config["parameter_status"] not in {"manufactured", "sourced", "unvalidated"}:
        raise ValueError("parameter_status must be manufactured, sourced or unvalidated")
    if not isinstance(config["parameter_provenance"], str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must identify the supplied coefficients and initial state")
    for name in ("periods", "cells_per_period"):
        value = config[name]
        if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 100_000:
            raise ValueError(f"{name} must be an integer in [1,100000]")
    count = config["periods"] * config["cells_per_period"]
    if count > 100_000:
        raise ValueError("total periodic cells cannot exceed 100000")
    origin = _finite_real(config["position_origin_m"], "position_origin_m")
    raw_times = np.asarray(config["times_s"], dtype=object)
    if raw_times.ndim != 1 or not raw_times.size or raw_times.size * count > 1_000_000:
        raise ValueError("times_s must be nonempty one-dimensional times within the reference memory limit")
    times = np.asarray([_finite_real(value, "times_s") for value in raw_times], dtype=float)
    if np.any(times < 0) or any(right <= left for left, right in zip(times[:-1], times[1:])):
        raise ValueError("times_s must be nonnegative and strictly increasing")
    if not isinstance(config["benchmark"], dict):
        raise ValueError("benchmark must supply all explicit constant coefficients")
    benchmark = PeriodicConductiveHeating(**config["benchmark"])
    # A finite representable period is needed by this CSV geometry contract.
    period = benchmark.period_m
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        edges = origin + period * (np.arange(count + 1, dtype=float) / config["cells_per_period"])
    if not np.all(np.isfinite(edges)) or any(right <= left for left, right in zip(edges[:-1], edges[1:])):
        raise ValueError("periodic cell edges must be finite and strictly increasing in supplied SI coordinates")
    return config, benchmark, edges, times, hashlib.sha256(payload).hexdigest(), payload


def save_plot(path, edges, times, profiles, budgets):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 3, figsize=(13, 7))
    centers = edges[:-1] / 2 + edges[1:] / 2
    for time, profile in zip(times, profiles):
        for axis, field, label in (
            (axes[0, 0], "magnetic_field_t", "Cell-average Bx [T]"),
            (axes[0, 1], "current_density_a_m2", "Cell-average Jz [A/m²]"),
            (axes[0, 2], "joule_power_w_m3", "Cell-average Joule power [W/m³]"),
            (axes[1, 0], "temperature_rise_k", "Cell-average temperature rise [K]"),
        ):
            axis.plot(centers, getattr(profile, field), label=f"{time:g} s")
            axis.set_xlabel("Periodic material position [m]")
            axis.set_ylabel(label)
    for field, label in (("mean_magnetic_energy_j_m3", "Magnetic"), ("mean_thermal_energy_j_m3", "Thermal"),
                         ("mean_total_energy_j_m3", "Combined")):
        axes[1, 1].plot(times, [getattr(budget, field) for budget in budgets], "o-", label=label)
    axes[1, 1].set_xlabel("Time [s]")
    axes[1, 1].set_ylabel("Spatial mean energy [J/m³]")
    axes[1, 2].plot(times, [budget.energy_residual_j_m3 for budget in budgets], "o-")
    axes[1, 2].set_xlabel("Time [s]")
    axes[1, 2].set_ylabel("Combined energy residual [J/m³]")
    for axis in axes.flat[:-1]:
        axis.legend()
    figure.suptitle("Manufactured periodic resistive and conductive reference")
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/conductive_heating_manufactured.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        source_hashes = {name: _sha(ROOT / name) for name in CODE_PATHS}
        config, benchmark, edges, times, config_hash, payload = load_configuration(args.config)
        profiles = [benchmark.cell_averages(edges, time) for time in times]
        budgets = [benchmark.energy_budget(time) for time in times]
        coefficients = {
            "magnetic_diffusivity_m2_s": benchmark.magnetic_diffusivity_m2_s,
            "thermal_diffusivity_m2_s": benchmark.thermal_diffusivity_m2_s,
            "resistive_energy_decay_rate_s": benchmark.resistive_energy_decay_rate_s,
            "thermal_mode_decay_rate_s": benchmark.thermal_mode_decay_rate_s,
            "period_m": benchmark.period_m,
        }
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".conductive-heating-", dir=args.output_dir.parent) as staging:
            output = Path(staging) / "result"
            output.mkdir()
            (output / "input_config.json").write_bytes(payload)
            with (output / "cell_averages.csv").open("w", newline="") as stream:
                writer = csv.writer(stream, lineterminator="\n")
                writer.writerow(CELL_UNITS)
                for time, profile in zip(times, profiles):
                    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
                        writer.writerow((float(time), index, float(left), float(right),
                                         float(profile.magnetic_field_t[index]), float(profile.current_density_a_m2[index]),
                                         float(profile.joule_power_w_m3[index]), float(profile.temperature_rise_k[index])))
            with (output / "energy_budget.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(ENERGY_UNITS), lineterminator="\n")
                writer.writeheader()
                writer.writerows({"time_s": float(time), **asdict(budget)} for time, budget in zip(times, budgets))
            if args.plot:
                save_plot(output / "overview.png", edges, times, profiles, budgets)
            if _sha(args.config) != config_hash or any(_sha(ROOT / name) != digest for name, digest in source_hashes.items()):
                raise ValueError("an input or source file changed during reference generation")
            summary = {
                "schema": RUN_SCHEMA, "model": "constant_property_periodic_resistive_induction_thermal_conduction",
                "is_flash_result": False, "is_measured_shot": False,
                "parameter_status": config["parameter_status"], "parameter_provenance": config["parameter_provenance"],
                "parameters": asdict(benchmark), "derived_coefficients": coefficients,
                "initial_state": "Bx(y,0)=initial_field_t*cos(wavenumber_rad_m*y); temperature rise zero",
                "current_orientation": "Jz=-dBx/dy/permeability_h_m",
                "position_origin_m": float(edges[0]), "position_end_m": float(edges[-1]),
                "periods": config["periods"], "cells_per_period": config["cells_per_period"],
                "times_s": times.tolist(), "cell_units": CELL_UNITS, "energy_budget_units": ENERGY_UNITS,
                "energy_budgets": [asdict(value) for value in budgets],
                "config_sha256": config_hash, "source_sha256": source_hashes,
                "output_sha256": {str(path.relative_to(output)): _sha(path) for path in sorted(output.rglob("*")) if path.is_file()},
                "assumptions": list(ASSUMPTIONS),
                "numerical_qualifications": {
                    "maximum_absolute_phase_rad": MAX_ABSOLUTE_PHASE_RAD,
                    "final_results": "finite float64; unrepresentable final values raise; sub-float final values round to zero",
                    "small_temperature_difference": "positive bounded 16-point Gauss-Legendre triangle integral when both dimensionless times are <=1",
                    "cli_derived_coefficients": "positive coefficients and period must be representable as float64; zero thermal diffusivity is allowed",
                },
                "limits": [
                    "Analytic reference generation is not a FLASH calculation or measured-shot prediction.",
                    "All property values, including thermal conductivity and volumetric heat capacity, are supplied constants.",
                    "No equation of state, mass-specific Cp-to-Cv conversion, Hall term, flow, radiation or phase change is included.",
                    "Finite-cell Joule power integrates Jz squared independently; it is not the square of mean Jz.",
                    "Temperature is a rise from the supplied zero-rise initial state, not an absolute material temperature.",
                    "The periodic combined magnetic and thermal mean is conserved; conduction redistributes heat without changing the mean.",
                    "A reference's algebraic energy residual does not establish a numerical solver's convergence or conservation.",
                    "A new output directory and one coordinated writer are required; staging and existence checks are not a concurrent-writer lock.",
                ],
            }
            (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during generation; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as error:
        parser.error(str(error))
    print(json.dumps({"output_dir": str(args.output_dir), "times": len(times), "cells_per_time": len(edges) - 1,
                      "parameter_status": config["parameter_status"],
                      "maximum_energy_residual_j_m3": max(abs(value.energy_residual_j_m3) for value in budgets)}, allow_nan=False))


if __name__ == "__main__":
    main()
