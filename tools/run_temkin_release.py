#!/usr/bin/env python3
"""Run a separately qualified coverage-dependent boundary-source calculation.

Temperatures are prescribed; no wall cooling, gap evolution or spectra are
computed. --demo uses a manufactured width-dependent heating fixture. Every
kinetic/energy parameter and its provenance is supplied in the configuration.
The existing constant-barrier packet-time sampler cannot consume this output.
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

from mhdlab.flash_surface import (
    as_boundary_temperature_proxy, as_current_uniform_temperature_proxy,
    load_surface_history, write_surface_history,
)
from mhdlab.temkin_release import LAW_IDENTITY, TemkinReservoir, TemkinInversionError
from mhdlab.temkin_history import evolve_temkin_release
from mhdlab.temkin_io import write_temkin_history
from mhdlab.thermal_release import BindingReservoir, SurfaceInventory
from mhdlab.surface_release_history import evolve_surface_release
from mhdlab.surface_release_io import write_release_history
from tools.run_surface_release import manufactured_surface, json_value, save_plot


SCHEMA = "mhdlab-temkin-release-config-v1"
CODE_PATHS = (
    "mhdlab/temkin_release.py", "mhdlab/temkin_history.py", "mhdlab/temkin_io.py",
    "mhdlab/thermal_release.py", "mhdlab/surface_release_history.py",
    "mhdlab/surface_release_io.py", "mhdlab/flash_surface.py",
    "mhdlab/hdf_integrity.py", "mhdlab/surface_sources.py", "mhdlab/constants.py",
    "tools/run_surface_release.py", "tools/run_temkin_release.py",
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_parameters(path: Path, width_bins: int):
    raw = path.read_bytes()
    config = json.loads(raw)
    required = {"schema", "parameter_status", "parameter_provenance",
                "current_direction_length_m", "subdivisions_per_interval", "reservoirs"}
    if not isinstance(config, dict) or set(config) != required or config["schema"] != SCHEMA:
        raise ValueError(f"configuration must contain exactly the {SCHEMA} fields")
    if config["parameter_status"] not in {"manufactured", "sourced", "unvalidated"}:
        raise ValueError("parameter_status must be manufactured, sourced, or unvalidated")
    if not isinstance(config["parameter_provenance"], str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must be a nonempty source description")
    length = config["current_direction_length_m"]
    if isinstance(length, bool) or not isinstance(length, Real) or not math.isfinite(length) or length <= 0:
        raise ValueError("current_direction_length_m must be a finite positive number")
    pools = config["reservoirs"]
    if not isinstance(pools, list) or not pools:
        raise ValueError("reservoirs must be a nonempty list")
    reservoirs, coverages = [], []
    for item in pools:
        if not isinstance(item, dict):
            raise ValueError("each reservoir must be a configuration object")
        supplied = dict(item)
        initial = supplied.pop("initial_coverage_m2")
        if np.isscalar(initial):
            initial = np.full(width_bins, SurfaceInventory([[initial]]).coverage_m2[0, 0])
        elif np.asarray(initial, dtype=object).shape != (width_bins,):
            raise ValueError("initial_coverage_m2 must be a scalar or one value per width bin")
        reservoirs.append(TemkinReservoir(**supplied))
        coverages.append(initial)
    return config, tuple(reservoirs), SurfaceInventory(coverages), hashlib.sha256(raw).hexdigest()


def _integrated(history, length: float) -> dict:
    widths = np.diff(history.x_edges_m)
    initial = np.sum(history.coverage_m2[0] * widths, axis=-1) * length
    remaining = np.sum(history.coverage_m2[-1] * widths, axis=-1) * length
    released = np.sum(history.width_integral(), axis=0) * length
    return {
        "initial_particles": initial,
        "remaining_particles": remaining,
        "released_particles": released,
        "wall_energy_debit_j": np.sum(history.width_integral("total_wall_energy_debit_j_m2"), axis=0) * length,
        "particle_balance_residual": [math.fsum((float(a), -float(b), -float(c)))
                                      for a, b, c in zip(initial, remaining, released)],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--demo", action="store_true")
    source.add_argument("--surface", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--subdivisions", type=int)
    parser.add_argument("--boundary-proxy-reference")
    parser.add_argument("--current-uniformity-reference")
    parser.add_argument("--constant-barrier-control", action="store_true",
                        help="Also run the separate alpha=0 kinetic sensitivity, retaining the same wall-energy parameters")
    parser.add_argument("--output-dir", type=Path, required=True, help="A new directory; existing paths are refused")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.surface is not None and args.config is None:
        parser.error("--surface requires an explicit --config")
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path; existing outputs are never overwritten")
    config_path = args.config or ROOT / "configs/temkin_manufactured.json"
    try:
        surface_hash = _sha(args.surface) if args.surface is not None else None
        surface = manufactured_surface(width_bins=24, time_intervals=16) if args.demo else load_surface_history(args.surface)
        if args.boundary_proxy_reference is not None:
            surface = as_boundary_temperature_proxy(surface, reference=args.boundary_proxy_reference)
        if args.current_uniformity_reference is not None:
            surface = as_current_uniform_temperature_proxy(surface, reference=args.current_uniformity_reference)
        config, pools, inventory, config_hash = load_parameters(config_path, len(surface.x_edges_m) - 1)
        subdivisions = args.subdivisions if args.subdivisions is not None else config["subdivisions_per_interval"]
        history = evolve_temkin_release(surface, inventory, pools,
                                       parameter_provenance=config["parameter_provenance"],
                                       subdivisions_per_interval=subdivisions)
        length = config["current_direction_length_m"]
        totals = _integrated(history, length)
        control = None
        if args.constant_barrier_control:
            baseline_pools = tuple(BindingReservoir(**{name: getattr(pool, name) for name in (
                "name", "state_id", "attempt_frequency_s", "activation_energy_ev",
                "desorption_energy_ev", "internal_energy_ev")}) for pool in pools)
            control = evolve_surface_release(surface, inventory, baseline_pools,
                                              parameter_provenance="Constant-barrier sensitivity of " + config["parameter_provenance"],
                                              subdivisions_per_interval=subdivisions)
        if _sha(config_path) != config_hash or (args.surface is not None and _sha(args.surface) != surface_hash):
            raise ValueError("an input file changed during the calculation")
        summary = {
            "schema": "mhdlab-temkin-run-v1", "law_identity": LAW_IDENTITY,
            "scope": "prescribed-temperature coverage-dependent first-order release; no solved wall cooling",
            "parameter_status": config["parameter_status"], "parameter_provenance": config["parameter_provenance"],
            "reservoirs": [asdict(pool) for pool in pools],
            "surface_attributes": dict(history.surface_attributes),
            "config_sha256": config_hash, "surface_input_sha256": surface_hash,
            "current_direction_length_m": length,
            "subdivisions_per_input_interval": subdivisions, "intervals": len(history.time_edges_s) - 1,
            "resolved_source": totals, "source_contrast": history.release_contrast(),
            "spatial_uniform_control": "Redistribute the resolved interval/state number and energy moments; retain the local temperature mixture for velocities",
            "constant_barrier_sensitivity": _integrated(control, length) if control is not None else None,
            "source_sha256": {name: _sha(ROOT / name) for name in CODE_PATHS},
            "limits": [
                "Manufactured examples establish numerical behavior, not H2O/SS304 calibration.",
                "All reservoirs here use constant supplied attempt frequencies; no Eyring prefactor is inferred.",
                "Rate barriers and constant wall-energy costs have independent provenance.",
                "This version has no law-specific packet-time adapter; the first-order packet sampler refuses its schema.",
                "This width-only source requires current-uniform temperature or an explicit justification of the current average.",
                "Varying temperature uses midpoint freezing; refine integration, input snapshot cadence and surface bins separately.",
            ],
        }
        # Stage a complete result beside its destination; failures leave no
        # apparently successful output directory or partially populated input.
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".temkin-", dir=args.output_dir.parent) as stage:
            out = Path(stage) / "result"
            out.mkdir()
            write_surface_history(out / "surface_input.h5", surface)
            write_temkin_history(out / "temkin_release.h5", history)
            if control is not None:
                write_release_history(out / "constant_barrier_control.h5", control)
            with (out / "integrated_source.csv").open("w", newline="") as stream:
                writer = csv.writer(stream, lineterminator="\n")
                writer.writerow(["time_start_s", "time_stop_s", "reservoir", "state_id", "released_particles", "wall_energy_debit_j"])
                counts = history.width_integral() * length
                energies = history.width_integral("total_wall_energy_debit_j_m2") * length
                for interval in range(len(history.time_edges_s) - 1):
                    for index, pool in enumerate(pools):
                        writer.writerow([history.time_edges_s[interval], history.time_edges_s[interval + 1],
                                         pool.name, pool.state_id, counts[interval, index], energies[interval, index]])
            if args.plot:
                save_plot(history, out / "overview.png")
            summary["output_sha256"] = {path.name: _sha(path) for path in out.iterdir() if path.is_file()}
            (out / "summary.json").write_text(json.dumps(json_value(summary), indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during calculation; refusing replacement")
            os.rename(out, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, TemkinInversionError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_dir": str(args.output_dir), "law_identity": LAW_IDENTITY,
                      "released_particles": totals["released_particles"].tolist(),
                      "parameter_status": config["parameter_status"]}, allow_nan=False))


if __name__ == "__main__":
    main()
