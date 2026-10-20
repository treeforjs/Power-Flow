"""Export analytic constant-property diffusion references for FLASH checks.

The default configuration is manufactured. This program runs no MHD solver;
its cell averages and half-space energy integrals are comparison references.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from numbers import Real
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.heating_benchmarks import ConstantConductivityHalfSpace
from mhdlab.heating_comparison import HeatingProfile, load_heating_profile_csv, write_heating_profile_csv


def load_configuration(path: Path):
    raw = path.read_bytes()
    config = json.loads(raw)
    if config.get("schema") != "mhdlab-heating-benchmark-config-v1":
        raise ValueError("unsupported heating benchmark configuration")
    if config.get("parameter_status") not in {"manufactured", "configured"}:
        raise ValueError("parameter_status must be manufactured or configured")
    if not isinstance(config.get("parameter_provenance"), str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance is required")
    count = config["depth_cells"]
    if isinstance(count, bool) or not isinstance(count, int) or not 1 <= count <= 1_000_000:
        raise ValueError("depth_cells must be an integer from 1 to 1000000")
    stop = config["depth_stop_m"]
    if isinstance(stop, bool) or not isinstance(stop, Real) or not np.isfinite(stop) or stop <= 0:
        raise ValueError("depth_stop_m must be finite and positive")
    times = np.asarray(config["times_s"], dtype=object)
    if (times.ndim != 1 or times.size == 0 or any(isinstance(v, bool) or not isinstance(v, Real)
            or not np.isfinite(v) or v <= 0 for v in times)):
        raise ValueError("times_s must contain positive finite real times")
    times = np.asarray(times, dtype=float)
    if np.any(np.diff(times) <= 0):
        raise ValueError("times_s must be strictly increasing")
    benchmark = ConstantConductivityHalfSpace(**config["benchmark"])
    return config, benchmark, np.linspace(0, float(stop), count + 1), times, hashlib.sha256(raw).hexdigest()


def save_plot(output: Path, edges, times, references, budgets):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    center_um = (edges[:-1] + 0.5 * np.diff(edges)) * 1e6
    for time, values in zip(times, references):
        label = f"{time * 1e9:g} ns"
        axes[0, 0].plot(center_um, values.magnetic_field_t, label=label)
        axes[0, 1].plot(center_um, values.current_density_a_m2, label=label)
        axes[1, 0].plot(center_um, values.joule_power_w_m3, label=label)
    for ax, units in [(axes[0, 0], "Cell-average Bx [T]"), (axes[0, 1], "Cell-average Jz [A/m²]"),
                      (axes[1, 0], "Cell-average Joule power [W/m³]")]:
        ax.set_xlabel("Inward material depth [µm]")
        ax.set_ylabel(units)
        ax.legend()
    for key, label in [("magnetic_energy_j_m2", "Magnetic storage"), ("joule_heat_j_m2", "Joule heat"),
                       ("poynting_input_j_m2", "Incoming Poynting energy")]:
        axes[1, 1].plot(times * 1e9, [getattr(value, key) for value in budgets], "o-", label=label)
    axes[1, 1].set_xlabel("Time [ns]")
    axes[1, 1].set_ylabel("Half-space energy [J/m²]")
    axes[1, 1].legend()
    fig.suptitle("Analytic constant-property diffusion reference; no FLASH run")
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/heating_diffusion_manufactured.json")
    parser.add_argument("--output-dir", type=Path, default=Path("runs/heating_benchmark_reference"))
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    try:
        config, benchmark, edges, times, config_hash = load_configuration(args.config)
        references = [benchmark.cell_averages(edges, float(time)) for time in times]
        budgets = [benchmark.energy_budget(float(time)) for time in times]
        output = args.output_dir
        profiles = output / "profiles"
        profiles.mkdir(parents=True, exist_ok=True)
        exported = []
        for index, (time, values) in enumerate(zip(times, references)):
            profile = HeatingProfile(float(time), edges, values.magnetic_field_t, values.current_density_a_m2,
                                     values.joule_power_w_m3, "synthetic", "analytic-half-space-reference",
                                     f"{config['parameter_provenance']}; config-sha256={config_hash}")
            csv_path, metadata_path = write_heating_profile_csv(profiles / f"time_{index:03d}.csv", profile)
            restored = load_heating_profile_csv(csv_path, metadata_path=metadata_path)
            for name in ("depth_edges_m", "magnetic_field_t", "current_density_a_m2", "joule_power_w_m3"):
                np.testing.assert_array_equal(getattr(restored, name), getattr(profile, name))
            exported.append({"time_s": float(time), "csv": str(csv_path.relative_to(output)),
                             "metadata": str(metadata_path.relative_to(output)),
                             "csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                             "metadata_sha256": hashlib.sha256(metadata_path.read_bytes()).hexdigest()})
        with (output / "energy_budget.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["time_s", *asdict(budgets[0])], lineterminator="\n")
            writer.writeheader()
            writer.writerows({"time_s": float(time), **asdict(budget)} for time, budget in zip(times, budgets))
        center = edges[:-1] + 0.5 * np.diff(edges)
        with (output / "adiabatic_point_temperature.csv").open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(["time_s", "cell_center_depth_m", "point_adiabatic_temperature_rise_k"])
            for time in times:
                writer.writerows((float(time), float(depth), float(rise)) for depth, rise in
                                 zip(center, benchmark.adiabatic_temperature_rise(center, float(time))))
        summary = {
            "scope": "analytic constant-conductivity half-space verification references",
            "is_flash_result": False, "is_measured_shot": False,
            "parameter_status": config["parameter_status"], "parameter_provenance": config["parameter_provenance"],
            "parameters": asdict(benchmark), "magnetic_diffusivity_m2_s": benchmark.magnetic_diffusivity_m2_s,
            "config_sha256": config_hash, "depth_cells": len(edges) - 1, "depth_stop_m": float(edges[-1]),
            "profiles": exported, "energy_budgets": [asdict(value) for value in budgets],
            "finite_depth_boundary_B_fraction": [float(benchmark.fields(float(edges[-1]), float(time)).magnetic_field_t / benchmark.boundary_field_t)
                                                  if benchmark.boundary_field_t != 0 else None for time in times],
            "source_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in
                              ["mhdlab/heating_benchmarks.py", "mhdlab/heating_comparison.py", "tools/run_heating_benchmark.py"]},
            "limits": ["constant coefficients, stationary material, imposed field step, no thermal conduction or phase change",
                       "global energy integrals cover the infinite half-space; finite-depth fields need outer-boundary verification",
                       "temperature output is a point diagnostic at cell centres; no finite surface temperature exists for this ideal step",
                       "reference export/readback is not a FLASH solver verification"],
        }
        (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
        if args.plot:
            save_plot(output / "overview.png", edges, times, references, budgets)
    except (OSError, ValueError, TypeError, KeyError) as exc:
        parser.exit(2, f"{exc}\n")
    print(json.dumps({"output_dir": str(output), "profiles": len(times), "parameter_status": config["parameter_status"],
                      "maximum_energy_residual_j_m2": max(abs(value.energy_residual_j_m2) for value in budgets)}))


if __name__ == "__main__":
    main()
