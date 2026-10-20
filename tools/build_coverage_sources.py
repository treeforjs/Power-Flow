#!/usr/bin/env python3
"""Build upstream desorption sources from one explicitly defined parent inventory."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np

from mhdlab.coverage_inventory import (
    InventoryPoolAllocation, MonolayerDefinition, ParentSurfaceInventory,
    build_coverage_source_sweep,
)
from mhdlab.planar_release_io import load_planar_surface_history, write_planar_surface_history
from mhdlab.source_convergence import manufactured_planar_surface
from mhdlab.temkin_planar_io import write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir

CODE_PATHS = (
    "tools/build_coverage_sources.py", "mhdlab/coverage_inventory.py",
    "mhdlab/temkin_planar.py", "mhdlab/temkin_planar_io.py", "mhdlab/temkin_release.py",
    "mhdlab/temkin_history.py", "mhdlab/temkin_io.py", "mhdlab/thermal_release.py",
    "mhdlab/planar_release.py", "mhdlab/planar_release_io.py", "mhdlab/source_convergence.py",
    "mhdlab/surface_release_history.py", "mhdlab/surface_sources.py", "mhdlab/flash_surface.py",
    "mhdlab/hdf_integrity.py", "mhdlab/constants.py",
)
CONFIG_KEYS = {
    "schema", "parameter_status", "parameter_provenance", "geometry",
    "subdivisions_per_interval", "parent", "coverage_ml", "allocations",
}
PARENT_KEYS = {
    "parent_id", "monolayer", "patch_multiplier", "profile_provenance", "status", "provenance",
}
STATUSES = {"manufactured", "reconstructed", "sourced", "unvalidated"}
SUMMARY_FIELDS = (
    "coverage_index", "coverage_ml", "reservoir", "state_id", "initial_parent_number",
    "initial_pool_number", "released_pool_number", "remaining_pool_number",
    "inventory_balance_residual_number", "total_wall_energy_debit_j",
)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_config(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def nonfinite(value):
        raise ValueError(f"nonfinite JSON constant: {value}")

    raw = path.read_bytes()
    config = json.loads(raw, object_pairs_hook=pairs, parse_constant=nonfinite)
    if (not isinstance(config, dict) or set(config) != CONFIG_KEYS
            or config["schema"] != "mhdlab-coverage-source-config-v1"):
        raise ValueError("configuration requires exactly the mhdlab-coverage-source-config-v1 fields")
    if config["parameter_status"] not in STATUSES:
        raise ValueError("parameter_status requires an explicit supported input label")
    if not isinstance(config["parameter_provenance"], str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must be nonempty")
    parent = config["parent"]
    if not isinstance(parent, dict) or set(parent) != PARENT_KEYS:
        raise ValueError("parent requires exactly the explicit monolayer, profile and provenance fields")
    records = config["allocations"]
    if (not isinstance(records, list) or not records
            or any(not isinstance(row, dict) or set(row) != {"reservoir", "fraction", "provenance", "status"}
                   for row in records)):
        raise ValueError("allocations require explicit reservoir, fraction, provenance and status records")
    allocations = tuple(InventoryPoolAllocation(
        reservoir=TemkinReservoir(**row["reservoir"]), fraction=row["fraction"],
        provenance=row["provenance"], status=row["status"],
    ) for row in records)
    parameters = dict(parent)
    parameters["monolayer"] = MonolayerDefinition(**parameters["monolayer"])
    parameters["allocations"] = allocations
    return raw, config, ParentSurfaceInventory(**parameters)


def load_surface(args, config, parent):
    if args.demo:
        if (config["parameter_status"] != "manufactured" or parent.status != "manufactured"
                or parent.monolayer.status != "manufactured"
                or any(row.status != "manufactured" for row in parent.allocations)):
            raise ValueError("--demo requires manufactured parameter, parent, monolayer and allocation labels")
        geometry = config["geometry"]
        fields = {"width_bins", "current_bins", "snapshot_intervals", "duration_s", "width_m", "current_length_m"}
        if not isinstance(geometry, dict) or set(geometry) != fields:
            raise ValueError("--demo requires the complete manufactured geometry fixture")
        return manufactured_planar_surface(**geometry)
    if config["geometry"] is not None:
        raise ValueError("external --surface requires geometry:null; input file owns its geometry")
    return load_planar_surface_history(args.surface)


def integrate_density(value, area):
    # Sum before float64 rounding so small patch contributions can accumulate.
    total = np.sum(np.asarray(value, dtype=np.longdouble)
                   * np.asarray(area, dtype=np.longdouble), dtype=np.longdouble)
    if not np.isfinite(total) or abs(total) > np.finfo(float).max:
        raise ValueError("source summary exceeds the finite numeric range")
    result = float(total)
    if total != 0 and result == 0:
        raise ValueError("nonzero source summary is below the representable numeric range")
    return result


def source_rows(sweep):
    records, rows = [], []
    for index, (coverage, history) in enumerate(zip(sweep.coverage_ml, sweep.histories)):
        area = history.patch_areas_m2
        initial = [integrate_density(history.coverage_m2[0, p], area) for p in range(len(history.reservoirs))]
        released = [integrate_density(history.released_m2[:, p], area) for p in range(len(initial))]
        remaining = [integrate_density(history.coverage_m2[-1, p], area) for p in range(len(initial))]
        residual = [integrate_density(history.inventory_balance_residual_m2[:, p], area) for p in range(len(initial))]
        debit = [integrate_density(history.total_wall_energy_debit_j_m2[:, p], area) for p in range(len(initial))]
        record = dict(index=index, coverage_ml=float(coverage), source_file=f"source_{index:04d}.h5",
                      initial_parent_number=float(sweep.plan.parent_number[index]),
                      initial_pool_number=initial, released_pool_number=released,
                      remaining_pool_number=remaining, total_wall_energy_debit_j=debit)
        records.append(record)
        for p, pool in enumerate(history.reservoirs):
            rows.append(dict(zip(SUMMARY_FIELDS, (
                index, float(coverage), pool.name, pool.state_id, record["initial_parent_number"],
                initial[p], released[p], remaining[p], residual[p], debit[p],
            ))))
    return records, rows


def plot_sources(records, path, parameter_status):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    coverage = [row["coverage_ml"] for row in records]
    axes[0].plot(coverage, [row["initial_parent_number"] for row in records], "o-", label="Initial parent inventory")
    axes[0].plot(coverage, [math.fsum(row["released_pool_number"]) for row in records], "s-", label="Released parent equivalents")
    axes[1].plot(coverage, [math.fsum(row["total_wall_energy_debit_j"]) for row in records], "o-")
    axes[0].set(ylabel="Physical molecule number")
    axes[0].legend(fontsize=8)
    axes[1].set(ylabel="Diagnostic wall energy debit [J]")
    for axis in axes:
        axis.set(xlabel="Declared monolayer equivalents")
        axis.grid(alpha=.25)
    figure.suptitle(f"Parameter status: {parameter_status}\nPrescribed-temperature parent release; no spectral calculation")
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--demo", action="store_true")
    group.add_argument("--surface", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    if args.surface is not None and args.config is None:
        parser.error("--surface requires an explicit --config")
    config_path = args.config or ROOT / "configs/coverage_source_manufactured.json"
    try:
        code_hashes = {name: sha(ROOT / name) for name in CODE_PATHS}
        raw_config, config, parent = read_config(config_path)
        config_hash = hashlib.sha256(raw_config).hexdigest()
        input_surface_hash = sha(args.surface) if args.surface is not None else None
        surface = load_surface(args, config, parent)
        # This budget bounds aggregate output/integration work for this interactive export.
        subdivisions = config["subdivisions_per_interval"]
        knots = config["coverage_ml"]
        if not isinstance(knots, list) or not knots:
            raise ValueError("coverage_ml must be a nonempty list of coverage knots")
        if isinstance(subdivisions, bool) or not isinstance(subdivisions, int) or subdivisions < 1:
            raise ValueError("subdivisions_per_interval must be a positive integer")
        work = (surface.time_s.size - 1) * subdivisions * surface.patch_areas_m2.size * len(parent.allocations) * len(knots)
        if work > 200_000:
            raise ValueError("coverage sweep exceeds the interactive component budget")
        sweep = build_coverage_source_sweep(
            surface, parent, knots, parameter_provenance=config["parameter_provenance"],
            subdivisions_per_interval=subdivisions,
        )
        records, rows = source_rows(sweep)
        report = dict(
            schema="mhdlab-coverage-source-run-v1", parameter_status=config["parameter_status"],
            parameter_provenance=config["parameter_provenance"], input_config_sha256=config_hash,
            input_surface_sha256=input_surface_hash, source_sha256=code_hashes,
            runtime={"python": sys.version, **{name: importlib.metadata.version(name)
                     for name in ("numpy", "scipy", "h5py")}},
            thermal_input_attributes=sweep.histories[0].surface_attributes,
            coverage_source_sweep=sweep.manifest(), sources=records,
            limits=[
                "One physical parent inventory is split among independent pools and supplied outgoing states; these fractions are not dissociation yields or shared-pool branching.",
                "Declared monolayers and kinetic reference coverage are independent; no actual water surface convention or spectrum is inferred.",
                "Prescribed wall temperatures are unchanged by the diagnostic energy debit; no coupled heating/desorption feedback is solved.",
                "Stored matched-width datasets preserve number/energy moments. Use original component temperatures, local emission laws and existing paired boundary interfaces to retain the full time/state/velocity mixture.",
                "No gap evolution, collision/radiation solve, FLASH run or host implementation is performed.",
                "Input status labels are caller declarations; separate thermal, monolayer, allocation and parameter labels are retained.",
                "Staged publication assumes a coordinated single writer; ordinary pre-existing paths are refused.",
            ],
        )
        from mhdlab.planar_release_io import plain_metadata
        report["thermal_input_attributes"] = plain_metadata(report["thermal_input_attributes"])
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".coverage-sources-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            (output / "input_config.json").write_bytes(raw_config)
            if args.demo:
                write_planar_surface_history(output / "surface_input.h5", surface)
            else:
                shutil.copyfile(args.surface, output / "surface_input.h5")
                if sha(output / "surface_input.h5") != input_surface_hash:
                    raise ValueError("input surface changed while copying")
            for record, history in zip(records, sweep.histories):
                target = output / record["source_file"]
                write_temkin_planar_release_history(target, history)
                record["source_sha256"] = sha(target)
            with (output / "summary.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=SUMMARY_FIELDS, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            if args.plot:
                plot_sources(records, output / "coverage_sources.png", config["parameter_status"])
            if (sha(config_path) != config_hash
                    or (args.surface is not None and sha(args.surface) != input_surface_hash)
                    or any(sha(ROOT / name) != digest for name, digest in code_hashes.items())):
                raise ValueError("input or implementation changed during export")
            report["output_sha256"] = {p.name: sha(p) for p in output.iterdir() if p.is_file()}
            (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during export; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as exc:
        parser.error(str(exc))
    print(json.dumps(dict(output_dir=str(args.output_dir), parameter_status=config["parameter_status"],
                          coverage_knots=len(records), parent_id=parent.parent_id)))


if __name__ == "__main__":
    main()
