#!/usr/bin/env python3
"""Export the periodic Hall-induction reference for independent verification.

All configured values are explicit manufactured verification inputs. This
program generates exact references; it does not run FLASH or another solver.
Ions are prescribed stationary, so it qualifies only the induction subsystem.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
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

from mhdlab.hall_induction_benchmark import PeriodicHallInduction


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path):
    raw = path.read_bytes()
    config = json.loads(raw)
    required = {"schema", "parameter_status", "parameter_provenance", "parameters", "cells_per_period", "time_s"}
    if not isinstance(config, dict) or set(config) != required or config["schema"] != "mhdlab-hall-induction-config-v1":
        raise ValueError("configuration must contain exactly the mhdlab-hall-induction-config-v1 fields")
    if config["parameter_status"] != "manufactured":
        raise ValueError("this induction reference accepts explicitly manufactured inputs only")
    if not isinstance(config["parameter_provenance"], str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must describe the supplied reference inputs")
    cells = config["cells_per_period"]
    if isinstance(cells, bool) or not isinstance(cells, int) or not 1 <= cells <= 100_000:
        raise ValueError("cells_per_period must be an integer from 1 to 100000")
    times = config["time_s"]
    if not isinstance(times, list) or not times or len(times) > 1000:
        raise ValueError("time_s must be a nonempty list with at most1000 times")
    model = PeriodicHallInduction(**config["parameters"])
    for time in times:
        model.energy_budget(time)
    if any(hi <= lo for lo, hi in zip(times[:-1], times[1:])):
        raise ValueError("time_s must be strictly increasing")
    if cells * len(times) > 2_000_000:
        raise ValueError("requested Hall reference exceeds the artifact size limit")
    return config, model, hashlib.sha256(raw).hexdigest()


def _plot(model, times, edges, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    z = edges[:-1] / 2 + edges[1:] / 2
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for time in times:
        fields = model.cell_averages(edges, time)
        ax[0].plot(z * 1e3, fields["magnetic_field_t"][:, 0], label=f"{time * 1e9:g} ns")
    ax[0].set(xlabel="Periodic coordinate z [mm]", ylabel="Cell-mean transverse Bx [T]")
    ax[0].legend()
    dense_t = np.linspace(times[0], times[-1], 201)
    budgets = [model.energy_budget(time) for time in dense_t]
    ax[1].plot(dense_t * 1e9, [b["magnetic_perturbation_energy_j_m3"] for b in budgets], label="Perturbation magnetic energy")
    ax[1].plot(dense_t * 1e9, [b["joule_heat_j_m3"] for b in budgets], label="Deposited Joule heat")
    ax[1].set(xlabel="Time [ns]", ylabel="Energy density [J/m³]")
    ax[1].legend()
    fig.suptitle("Manufactured prescribed-ion Hall-induction reference")
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/hall_induction_manufactured.json")
    parser.add_argument("--output-dir", type=Path, required=True, help="A new directory; existing paths are refused")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        config, model, config_hash = _load(args.config)
        times = config["time_s"]
        edges = np.linspace(0.0, model.wavelength_m, config["cells_per_period"] + 1)
        if np.any(np.diff(edges) <= 0):
            raise ValueError("configured cell partition is not representable")
        summaries = []
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".hall-reference-", dir=args.output_dir.parent) as stage:
            out = Path(stage) / "result"
            out.mkdir()
            units = {"z_left_m": "m", "z_right_m": "m"}
            for time_index, time in enumerate(times):
                fields = model.cell_averages(edges, time)
                values = [edges[:-1], edges[1:]]
                columns = ["z_left_m", "z_right_m"]
                for name, unit in (("magnetic_field_t", "T"), ("current_density_a_m2", "A m^-2"),
                                   ("resistive_electric_field_v_m", "V m^-1"), ("hall_electric_field_v_m", "V m^-1"),
                                   ("electric_field_v_m", "V m^-1"), ("electron_velocity_m_s", "m s^-1")):
                    for component, axis in enumerate("xyz"):
                        column = f"{name}_{axis}"
                        columns.append(column)
                        units[column] = unit
                        values.append(fields[name][:, component])
                for name in ("joule_power_w_m3", "hall_power_w_m3"):
                    columns.append(name)
                    units[name] = "W m^-3"
                    values.append(fields[name])
                path = out / f"fields_{time_index:03d}.csv"
                with path.open("w", newline="") as stream:
                    writer = csv.writer(stream, lineterminator="\n")
                    writer.writerow(columns)
                    writer.writerows(zip(*values))
                record = {"schema": "mhdlab-hall-induction-reference-v1", "time_s": time,
                          "parameter_status": "manufactured", "parameters": asdict(model),
                          "csv_sha256": _sha(path), "units": units.copy(),
                          "representation": "exact cell means; independently averaged Joule power",
                          "energy_budget": model.energy_budget(time)}
                path.with_suffix(".json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
                summaries.append({"time_s": time, "csv_file": path.name, "csv_sha256": _sha(path),
                                  "energy_budget": record["energy_budget"]})
            summary = {"schema": "mhdlab-hall-induction-run-v1", "parameter_status": "manufactured",
                       "parameter_provenance": config["parameter_provenance"], "parameters": asdict(model),
                       "wavenumber_m_inv": model.wavenumber_m_inv,
                       "decay_rate_s": model.decay_rate_s, "angular_frequency_s": model.angular_frequency_s,
                       "hall_period_s": math.tau / abs(model.angular_frequency_s) if model.angular_frequency_s else None,
                       "snapshots": summaries, "config_sha256": config_hash,
                       "source_sha256": {name: _sha(ROOT / name) for name in (
                           "mhdlab/hall_induction_benchmark.py", "mhdlab/thermal_release.py",
                           "mhdlab/constants.py", "tools/run_hall_induction_benchmark.py")},
                       "limits": ["No FLASH or other numerical solver run is claimed.",
                                  "Ions are prescribed stationary; Lorentz-force acceleration is excluded.",
                                  "Constant scalar resistivity, permeability and electron density; no electron pressure or inertia.",
                                  "No material-property calibration, Hall-MHD fluid qualification or kinetic boundary solution.",
                                  "Cell-mean current or field products do not equal cell-mean physical energy or Joule power."]}
            if args.plot:
                _plot(model, times, edges, out / "overview.png")
            summary["output_sha256"] = {path.name: _sha(path) for path in out.iterdir() if path.is_file()}
            (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
            if _sha(args.config) != config_hash:
                raise ValueError("configuration changed during export")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during export; refusing replacement")
            os.rename(out, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_dir": str(args.output_dir), "snapshots": len(times),
                      "angular_frequency_s": model.angular_frequency_s,
                      "max_energy_balance_residual_j_m3": max(abs(item["energy_budget"]["balance_residual_j_m3"]) for item in summaries)}))


if __name__ == "__main__":
    main()
