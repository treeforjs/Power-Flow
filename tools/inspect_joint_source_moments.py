#!/usr/bin/env python3
"""Export deterministic joint moments for resolved and uniform-width sources."""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.joint_source_moments import joint_source_moments
from mhdlab.planar_release_io import write_planar_surface_history
from mhdlab.source_convergence import manufactured_planar_surface
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import load_temkin_planar_release_history, write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir
from mhdlab.constants import QE

CODE_PATHS = (
    "tools/inspect_joint_source_moments.py", "mhdlab/joint_source_moments.py",
    "mhdlab/emission_time_moments.py", "mhdlab/temkin_emission_time.py",
    "mhdlab/temkin_planar.py", "mhdlab/temkin_planar_io.py", "mhdlab/temkin_history.py",
    "mhdlab/temkin_io.py", "mhdlab/temkin_release.py", "mhdlab/thermal_release.py",
    "mhdlab/planar_release.py", "mhdlab/planar_release_io.py", "mhdlab/source_convergence.py",
    "mhdlab/surface_release_history.py", "mhdlab/surface_sources.py",
    "mhdlab/flash_surface.py", "mhdlab/hdf_integrity.py", "mhdlab/constants.py",
)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path):
    def nonfinite(value):
        raise ValueError(f"nonfinite JSON constant: {value}")
    return json.loads(path.read_bytes(), object_pairs_hook=_pairs, parse_constant=nonfinite)


def demo_history(path):
    config = read_json(path)
    if (not isinstance(config, dict) or set(config) != {"schema", "parameter_status", "parameter_provenance", "geometry", "subdivisions_per_interval", "reservoirs"}
            or config["schema"] != "mhdlab-joint-source-demo-config-v1" or config["parameter_status"] != "manufactured"):
        raise ValueError("demo configuration requires exactly the manufactured mhdlab-joint-source-demo-config-v1 fields")
    if not isinstance(config["parameter_provenance"], str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must be nonempty")
    geometry = config["geometry"]
    if not isinstance(geometry, dict) or set(geometry) != {"width_bins", "current_bins", "snapshot_intervals", "duration_s", "width_m", "current_length_m"}:
        raise ValueError("geometry must declare the manufactured space/time fixture")
    for name in ("width_bins", "current_bins", "snapshot_intervals"):
        if isinstance(geometry[name], bool) or not isinstance(geometry[name], int) or geometry[name] < 1:
            raise ValueError(f"{name} must be a positive integer")
    subdivisions = config["subdivisions_per_interval"]
    if isinstance(subdivisions, bool) or not isinstance(subdivisions, int) or subdivisions < 1:
        raise ValueError("subdivisions_per_interval must be a positive integer")
    records = config["reservoirs"]
    if not isinstance(records, list) or not records:
        raise ValueError("reservoirs must be a nonempty list")
    if geometry["width_bins"] * geometry["current_bins"] * geometry["snapshot_intervals"] * subdivisions * len(records) > 200000:
        raise ValueError("demo exceeds the interactive component budget")
    surface = manufactured_planar_surface(**geometry)
    pools, initial = [], []
    for record in records:
        parameters = dict(record)
        coverage = parameters.pop("initial_coverage_m2")
        if (isinstance(coverage, bool) or not isinstance(coverage, (int, float))
                or not math.isfinite(coverage) or coverage < 0):
            raise ValueError("initial_coverage_m2 must be a finite nonnegative scalar")
        pools.append(TemkinReservoir(**parameters))
        initial.append(coverage)
    coverages = np.broadcast_to(np.asarray(initial)[:, None, None],
                               (len(pools), geometry["current_bins"], geometry["width_bins"]))
    history = evolve_temkin_planar_release(surface, coverages, pools,
                                          parameter_provenance=config["parameter_provenance"],
                                          subdivisions_per_interval=subdivisions)
    return config, surface, history


def associated_status(path, source_hash, history):
    if path is None:
        return "unreported"
    record = read_json(path)
    if (not isinstance(record, dict) or record.get("schema") not in {"mhdlab-temkin-planar-run-v1", "mhdlab-joint-source-run-v1"}
            or record.get("parameter_status") not in {"manufactured", "sourced", "unvalidated"}
            or not isinstance(record.get("output_sha256"), dict)
            or source_hash not in record["output_sha256"].values()
            or record.get("parameter_provenance") != history.parameter_provenance):
        raise ValueError("source-run summary must hash-associate the source, its provenance and a supported parameter status")
    if "reservoirs" in record and record["reservoirs"] != [asdict(pool) for pool in history.reservoirs]:
        raise ValueError("source-run summary reservoir definitions disagree")
    return record["parameter_status"]


def write_conditionals(path, moments, kind):
    first = moments["reservoirs"][0][kind][0]
    stats = list(first["resolved"])
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["reservoir", "state_id", "hypothesis", "bin_index", "bounds_json", *stats])
        for pool in moments["reservoirs"]:
            for index, entry in enumerate(pool[kind]):
                bounds = {key: value for key, value in entry.items() if key not in {"resolved", "uniform_width"}}
                for hypothesis in ("resolved", "uniform_width"):
                    writer.writerow([pool["name"], pool["state_id"], hypothesis, index,
                                     json.dumps(bounds, sort_keys=True, allow_nan=False),
                                     *(entry[hypothesis][key] for key in stats)])


def plot_moments(moments, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for pool in moments["reservoirs"]:
        entries = pool["conditional_width"]
        x = [(.5 * e["x_left_m"] + .5 * e["x_right_m"]) * 1e3 for e in entries]
        for hypothesis, style in (("resolved", "-o"), ("uniform_width", "--s")):
            label = f"{pool['name']}: {hypothesis}"
            for axis, key, factor in ((axes[0], "released_number", 1), (axes[1], "mean_time_s", 1e9),
                                      (axes[2], "mean_expected_translational_energy_j", 1 / QE)):
                values = [None if e[hypothesis][key] is None else factor * e[hypothesis][key] for e in entries]
                axis.plot(x, values, style, label=label)
    for axis, ylabel in zip(axes, ("Emitted physical number per width bin", "Conditional mean emission time [ns]", "Mean local expected translation [eV]")):
        axis.set(xlabel="Width coordinate [mm]", ylabel=ylabel)
        axis.legend(fontsize=6)
        axis.grid(alpha=.25)
    figure.suptitle("Prescribed-temperature source moments; no transport or spectrum calculation")
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--demo", action="store_true")
    source.add_argument("--source", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--source-run-summary", type=Path)
    parser.add_argument("--quadrature-order", type=int, default=8)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    if args.source is not None and args.config is not None:
        parser.error("--config applies only to --demo")
    if args.demo and args.source_run_summary is not None:
        parser.error("--source-run-summary applies only to --source")
    try:
        code_hashes = {name: sha(ROOT / name) for name in CODE_PATHS}
        config_path = (args.config or ROOT / "configs/joint_source_manufactured.json") if args.demo else None
        config_hash = sha(config_path) if config_path else None
        source_hash = sha(args.source) if args.source else None
        summary_hash = sha(args.source_run_summary) if args.source_run_summary else None
        if args.demo:
            config, surface, history = demo_history(config_path)
            parameter_status = config["parameter_status"]
        else:
            surface = None
            history = load_temkin_planar_release_history(args.source)
            parameter_status = associated_status(args.source_run_summary, source_hash, history)
        moments = joint_source_moments(history, quadrature_order=args.quadrature_order)
        report = {"schema": "mhdlab-joint-source-run-v1", "parameter_status": parameter_status,
                  "parameter_provenance": history.parameter_provenance,
                  "reservoirs": [asdict(pool) for pool in history.reservoirs],
                  "joint_moments": moments, "input_config_sha256": config_hash,
                  "input_source_sha256": source_hash, "input_summary_sha256": summary_hash,
                  "source_sha256": code_hashes,
                  "limits": ["Prescribed stationary planar source; no gap evolution, chemistry solve, spectra or applied wall cooling.",
                             "Expected translation 2*kB*T is a component conditional mean; its mixture variance excludes within-component kinetic-energy fluctuations.",
                             "Conditional bins are original source patches/intervals; no truncation at arbitrary diagnostic gate edges.",
                             "Time quadrature is finite; its refinement changes are not certified error bounds.",
                             "Staged publication assumes a coordinated single writer; ordinary pre-existing paths are refused."]}
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".joint-source-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            if args.demo:
                (output / "input_config.json").write_bytes(config_path.read_bytes())
                write_planar_surface_history(output / "surface_input.h5", surface)
                write_temkin_planar_release_history(output / "source_input.h5", history)
            else:
                shutil.copyfile(args.source, output / "source_input.h5")
                if sha(output / "source_input.h5") != source_hash:
                    raise ValueError("source changed while copying")
            write_conditionals(output / "conditional_width.csv", moments, "conditional_width")
            write_conditionals(output / "conditional_interval.csv", moments, "conditional_interval")
            if args.plot:
                plot_moments(moments, output / "joint_moments.png")
            if ((config_path and sha(config_path) != config_hash) or (args.source and sha(args.source) != source_hash)
                    or (args.source_run_summary and sha(args.source_run_summary) != summary_hash)
                    or any(sha(ROOT / name) != value for name, value in code_hashes.items())):
                raise ValueError("input or implementation changed during inspection")
            report["output_sha256"] = {p.name: sha(p) for p in output.iterdir() if p.is_file()}
            (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during inspection; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_dir": str(args.output_dir), "parameter_status": parameter_status,
                      "reservoirs": len(moments["reservoirs"])}))


if __name__ == "__main__":
    main()
