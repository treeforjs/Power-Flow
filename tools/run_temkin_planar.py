#!/usr/bin/env python3
"""Run local coverage-dependent release on resolved width/current patches.

--demo uses manufactured temperatures and parameters. External input requires a
canonical planar surface plus explicit Temkin parameters. No gap solver runs.
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

from mhdlab.flash_planar import as_planar_boundary_temperature_proxy
from mhdlab.planar_release_io import load_planar_surface_history, write_planar_surface_history, plain_metadata
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import write_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir, LAW_IDENTITY
from mhdlab.thermal_release import SurfaceInventory
from tools.run_planar_release import manufactured_surface, save_plot


CODE_PATHS = (
    "mhdlab/temkin_planar.py", "mhdlab/temkin_planar_io.py", "mhdlab/temkin_release.py",
    "mhdlab/temkin_history.py", "mhdlab/temkin_io.py", "mhdlab/thermal_release.py",
    "mhdlab/planar_release.py", "mhdlab/planar_release_io.py", "mhdlab/flash_planar.py",
    "mhdlab/flash_surface.py", "mhdlab/hdf_integrity.py", "mhdlab/surface_sources.py",
    "mhdlab/constants.py", "tools/run_planar_release.py", "tools/run_temkin_planar.py",
)


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_parameters(path, surface):
    raw = path.read_bytes()
    config = json.loads(raw)
    required = {"schema", "parameter_status", "parameter_provenance", "current_direction_length_m",
                "subdivisions_per_interval", "reservoirs"}
    if not isinstance(config, dict) or set(config) != required or config["schema"] != "mhdlab-temkin-release-config-v1":
        raise ValueError("configuration must contain exactly the mhdlab-temkin-release-config-v1 fields")
    if config["parameter_status"] not in {"manufactured", "sourced", "unvalidated"}:
        raise ValueError("parameter_status must be manufactured, sourced or unvalidated")
    if not isinstance(config["parameter_provenance"], str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must describe the supplied parameters")
    length = config["current_direction_length_m"]
    actual = float(surface.z_edges_m[-1] - surface.z_edges_m[0])
    if (isinstance(length, bool) or not isinstance(length, (int, float)) or not math.isfinite(length)
            or length <= 0 or not math.isclose(length, actual, rel_tol=8 * np.finfo(float).eps, abs_tol=0)):
        raise ValueError("configured current_direction_length_m conflicts with the actual planar z extent")
    shape = (len(surface.z_edges_m) - 1, len(surface.x_edges_m) - 1)
    records = config["reservoirs"]
    if not isinstance(records, list) or not records:
        raise ValueError("reservoirs must be a nonempty list")
    pools, coverage = [], []
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("each reservoir must be a parameter object")
        supplied = dict(record)
        initial = supplied.pop("initial_coverage_m2")
        if np.isscalar(initial):
            initial = np.full(shape, SurfaceInventory([[initial]]).coverage_m2[0, 0])
        elif np.asarray(initial, dtype=object).shape != shape:
            raise ValueError("initial_coverage_m2 must be a scalar or exact (nz,nx) array")
        pools.append(TemkinReservoir(**supplied))
        coverage.append(initial)
    validated = SurfaceInventory(np.asarray(coverage, dtype=object).reshape(len(pools), -1))
    return config, tuple(pools), validated.coverage_m2.reshape((len(pools), *shape)), hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--demo", action="store_true")
    source.add_argument("--surface", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--subdivisions", type=int)
    parser.add_argument("--boundary-proxy-reference")
    parser.add_argument("--include-full-face-uniform", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.surface is not None and args.config is None:
        parser.error("--surface requires explicit --config")
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    path = args.config or ROOT / "configs/temkin_manufactured.json"
    try:
        surface_hash = None if args.demo else _sha(args.surface)
        surface = manufactured_surface() if args.demo else load_planar_surface_history(args.surface)
        if args.boundary_proxy_reference is not None:
            surface = as_planar_boundary_temperature_proxy(surface, reference=args.boundary_proxy_reference)
        config, pools, coverage, config_hash = load_parameters(path, surface)
        subdivisions = args.subdivisions if args.subdivisions is not None else config["subdivisions_per_interval"]
        history = evolve_temkin_planar_release(surface, coverage, pools,
                                               parameter_provenance=config["parameter_provenance"],
                                               subdivisions_per_interval=subdivisions)
        quantities = {name: history.extensive_totals(name) for name in (
            "released_m2", "translational_energy_j_m2", "internal_energy_j_m2",
            "binding_energy_j_m2", "total_wall_energy_debit_j_m2")}
        initial = history.extensive_totals("coverage_m2")[0]
        remaining = history.extensive_totals("coverage_m2")[-1]
        released = quantities["released_m2"].sum(axis=0)
        summary = {
            "schema": "mhdlab-temkin-planar-run-v1", "law_identity": LAW_IDENTITY,
            "scope": "prescribed-temperature local planar depletion; no transport or thermal feedback solve",
            "parameter_status": config["parameter_status"], "parameter_provenance": config["parameter_provenance"],
            "reservoirs": [asdict(pool) for pool in pools],
            "surface_attributes": plain_metadata(history.surface_attributes),
            "patch_shape_z_x": list(history.patch_areas_m2.shape), "face_area_m2": float(history.patch_areas_m2.sum()),
            "initial_particles": initial.tolist(), "remaining_particles": remaining.tolist(),
            "released_particles": released.tolist(),
            "wall_energy_debit_j": quantities["total_wall_energy_debit_j_m2"].sum(axis=0).tolist(),
            "particle_balance_residual": [math.fsum((float(a), -float(b), -float(c))) for a, b, c in zip(initial, remaining, released)],
            "intervals": len(history.time_edges_s) - 1, "subdivisions_per_input_interval": subdivisions,
            "input_config_sha256": config_hash, "input_surface_sha256": surface_hash,
            "uniform_controls": {"default": "width_only_preserving_each_current_patch",
                                 "full_face_included": args.include_full_face_uniform},
            "source_sha256": {name: _sha(ROOT / name) for name in CODE_PATHS},
            "limits": ["Caller-supplied parameters retain their status; the example is manufactured, not H2O/SS304 data.",
                       "No averaging of current-direction temperature precedes nonlinear release.",
                       "Actual patch areas determine all extensive quantities; no extra current length factor is applied.",
                       "Moment controls do not by themselves assign molecular states or launch distributions.",
                       "A planar Temkin particle adapter is separate from this moment-source output.",
                       "Each independent pool uses its own coverage; shared surface occupancy is not inferred."]}
        if _sha(path) != config_hash or (args.surface is not None and _sha(args.surface) != surface_hash):
            raise ValueError("an input file changed during release")
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".temkin-planar-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            write_planar_surface_history(output / "surface_input.h5", surface)
            write_temkin_planar_release_history(output / "temkin_planar_release.h5", history,
                                                include_full_face_uniform=args.include_full_face_uniform)
            with (output / "integrated_source.csv").open("w", newline="") as stream:
                writer = csv.writer(stream, lineterminator="\n")
                writer.writerow(["time_start_s", "time_stop_s", "reservoir", "state_id", "released_particles",
                                 "translational_energy_j", "internal_energy_j", "binding_energy_j", "wall_energy_debit_j"])
                for interval in range(len(history.time_edges_s) - 1):
                    for index, pool in enumerate(pools):
                        writer.writerow([history.time_edges_s[interval], history.time_edges_s[interval + 1], pool.name, pool.state_id,
                                         *(values[interval, index] for values in quantities.values())])
            if args.plot:
                save_plot(history, output / "overview.png")
            summary["output_sha256"] = {p.name: _sha(p) for p in output.iterdir() if p.is_file()}
            (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during release; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_dir": str(args.output_dir), "released_particles": released.tolist(),
                      "parameter_status": config["parameter_status"]}, allow_nan=False))


if __name__ == "__main__":
    main()
