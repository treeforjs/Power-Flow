"""Run local planar surface release with an explicit width-only comparison.

--demo uses manufactured patch heating and declared verification parameters.
External heating requires --surface and --config; material-plane temperatures
require an explicit boundary approximation reference. No external solver runs.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.flash_planar import as_planar_boundary_temperature_proxy
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory, evolve_planar_release
from mhdlab.planar_release_io import (
    load_planar_surface_history, plain_metadata, write_planar_release_history,
    write_planar_surface_history,
)
from mhdlab.thermal_release import BindingReservoir, SurfaceInventory


def manufactured_surface() -> PlanarSurfaceHistory:
    """Analytic x/z patch averages with a quadratic time ramp; no shot claim."""
    x = np.linspace(-.0025, .0025, 25)
    z = np.linspace(0.0, .025, 7)
    time = np.linspace(0.0, 120e-9, 13)
    width = x[-1] - x[0]
    width_shape = (4.0 / width**2) * np.diff(x**3) / (3.0 * np.diff(x))
    current_shape = .6 + .8 * (z[:-1] + .5 * np.diff(z)) / z[-1]
    ramp = (time / time[-1])**2
    temperature = 300.0 + ramp[:, None, None] * (
        250.0 + 1200.0 * current_shape[None, :, None] * width_shape[None, None, :]
    )
    base = FlashSurfaceHistory(time, x, temperature[:, 0, :], "synthetic",
                              "manufactured-planar-patch-reference",
                              "tools/run_planar_release.py:manufactured_surface analytic patch averages",
                              "manufactured stationary gap-facing face")
    attributes = base.attributes
    attributes["manufactured_field"] = {
        "formula": "300 + (t/t_end)^2 * (250 + 1200 * mean[(2*x/width)^2] * mean[0.6+0.8*z/length])",
        "claim": "analytical verification field; no experimental or FLASH heating result",
    }
    return PlanarSurfaceHistory(time, x, z, temperature, attributes)


def load_parameters(path: Path, *, nz: int, nx: int, current_length_m: float):
    content = path.read_bytes()
    config = json.loads(content)
    if not isinstance(config, dict) or config.get("schema") != "mhdlab-surface-release-config-v1":
        raise ValueError("unsupported release configuration schema")
    if config.get("parameter_status") not in {"manufactured", "sourced", "unvalidated"}:
        raise ValueError("parameter_status must be manufactured, sourced, or unvalidated")
    if not isinstance(config.get("parameter_provenance"), str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must describe the supplied release parameters")
    if "current_direction_length_m" in config:
        length = config["current_direction_length_m"]
        if (isinstance(length, bool) or not isinstance(length, (int, float)) or not np.isfinite(length)
                or length <= 0.0 or not math.isclose(length, current_length_m, rel_tol=8 * np.finfo(float).eps, abs_tol=0.0)):
            raise ValueError("configured current_direction_length_m conflicts with the actual planar z extent")
    records = config.get("reservoirs")
    if not isinstance(records, list) or not records:
        raise ValueError("reservoirs must be a nonempty list")
    pools, coverage = [], []
    for record in records:
        if not isinstance(record, dict) or "initial_coverage_m2" not in record:
            raise ValueError("each reservoir requires explicit initial_coverage_m2")
        supplied = dict(record)
        initial = supplied.pop("initial_coverage_m2")
        if np.isscalar(initial):
            scalar = SurfaceInventory([[initial]]).coverage_m2[0, 0]
            initial = np.full((nz, nx), scalar)
        elif np.asarray(initial, dtype=object).shape != (nz, nx):
            raise ValueError("initial_coverage_m2 must be a scalar or exact (nz,nx) array per reservoir")
        pools.append(BindingReservoir(**supplied))
        coverage.append(initial)
    inventory = SurfaceInventory(np.asarray(coverage, dtype=object).reshape(len(pools), nz * nx))
    subdivisions = config.get("subdivisions_per_interval")
    if isinstance(subdivisions, bool) or not isinstance(subdivisions, int) or subdivisions < 1:
        raise ValueError("subdivisions_per_interval must be an explicit positive integer")
    return config, tuple(pools), inventory.coverage_m2.reshape(len(pools), nz, nx), hashlib.sha256(content).hexdigest()


def save_plot(history, output: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cumulative = history.released_m2.sum(axis=(0, 1))
    uniform = history.matched_uniform_width()["released_m2"].sum(axis=(0, 1))
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, values, label in (
        (axes[0, 0], history.temperature_k[-1], "Final interval temperature [K]"),
        (axes[0, 1], cumulative, "Cumulative local release [m⁻²]"),
        (axes[1, 0], uniform, "Width-only uniform release [m⁻²]"),
    ):
        mesh = ax.pcolormesh(history.x_edges_m * 1e3, history.z_edges_m * 1e3, values, shading="flat", rasterized=True)
        ax.set_xlabel("Foil width x [mm]")
        ax.set_ylabel("Current direction z [mm]")
        fig.colorbar(mesh, ax=ax, label=label)
    for index, pool in enumerate(history.reservoirs):
        axes[1, 1].plot(history.time_edges_s[1:] * 1e9,
                        np.cumsum(history.extensive_totals()[:, index]), label=pool.name)
    axes[1, 1].set_xlabel("Time [ns]")
    axes[1, 1].set_ylabel("Cumulative released particles")
    axes[1, 1].legend()
    fig.suptitle(f"Prescribed local planar release; source={history.surface_attributes['source_kind']}")
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(args) -> dict:
    config_path = args.config or ROOT / "configs/surface_release_manufactured.json"
    surface = manufactured_surface() if args.demo else load_planar_surface_history(args.surface)
    original_surface_hash = None if args.demo else _sha256(args.surface)
    if args.boundary_proxy_reference is not None:
        surface = as_planar_boundary_temperature_proxy(surface, reference=args.boundary_proxy_reference)
    config, pools, coverage, config_hash = load_parameters(
        config_path, nz=surface.z_edges_m.size - 1, nx=surface.x_edges_m.size - 1,
        current_length_m=float(surface.z_edges_m[-1] - surface.z_edges_m[0]),
    )
    subdivisions = args.subdivisions if args.subdivisions is not None else config["subdivisions_per_interval"]
    history = evolve_planar_release(surface, coverage, pools, parameter_provenance=config["parameter_provenance"],
                                    subdivisions_per_interval=subdivisions)
    output = args.output_dir
    outputs = [output / name for name in ("surface_input.h5", "planar_release.h5", "summary.json", "integrated_source.csv")]
    if args.plot:
        outputs.append(output / "overview.png")
    inputs = {config_path.resolve()}
    if args.surface is not None:
        inputs.add(args.surface.resolve())
    if any(path.resolve() in inputs for path in outputs):
        raise ValueError("output paths must not overwrite the supplied surface or parameter file")
    output.mkdir(parents=True, exist_ok=True)
    write_planar_surface_history(output / "surface_input.h5", surface)
    write_planar_release_history(output / "planar_release.h5", history,
                                 include_full_face_uniform=args.include_full_face_uniform)
    quantities = {name: history.extensive_totals(name) for name in (
        "released_m2", "translational_energy_j_m2", "internal_energy_j_m2",
        "binding_energy_j_m2", "total_wall_energy_debit_j_m2",
    )}
    with (output / "integrated_source.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["time_start_s", "time_stop_s", "reservoir", "state_id", "released_particles",
                         "translational_energy_j", "internal_energy_j", "binding_energy_j", "wall_energy_debit_j"])
        for interval in range(history.time_edges_s.size - 1):
            for index, pool in enumerate(pools):
                writer.writerow([history.time_edges_s[interval], history.time_edges_s[interval + 1], pool.name, pool.state_id,
                                 *(values[interval, index] for values in quantities.values())])
    if args.plot:
        save_plot(history, output / "overview.png")
    initial = history.extensive_totals("coverage_m2")[0]
    remaining = history.extensive_totals("coverage_m2")[-1]
    released = quantities["released_m2"].sum(axis=0)
    debit = quantities["total_wall_energy_debit_j_m2"].sum(axis=0)
    artifacts = ["surface_input.h5", "planar_release.h5", "integrated_source.csv"] + (["overview.png"] if args.plot else [])
    summary = {
        "schema": "mhdlab-planar-release-summary-v1",
        "scope": "prescribed-temperature local planar depletion and delivered number/energy moments",
        "surface_attributes": plain_metadata(history.surface_attributes),
        "parameter_status": config["parameter_status"], "parameter_provenance": config["parameter_provenance"],
        "reservoirs": [asdict(pool) for pool in pools],
        "patch_shape_z_x": list(history.patch_areas_m2.shape), "face_area_m2": float(history.patch_areas_m2.sum()),
        "initial_particles": initial.tolist(), "remaining_particles": remaining.tolist(),
        "released_particles": released.tolist(), "wall_energy_debit_j": debit.tolist(),
        "particle_balance_residual": [math.fsum((float(a), -float(b), -float(c))) for a, b, c in zip(initial, remaining, released)],
        "intervals": history.time_edges_s.size - 1, "subdivisions_per_input_interval": subdivisions,
        "uniform_controls": {
            "default": "width_only",
            "width_only": "preserves each interval/reservoir/z-patch number and energy; current profile retained",
            "full_face_included": bool(args.include_full_face_uniform),
            "distribution_limit": "moment matching does not assign a joint particle velocity/internal-state distribution",
        },
        "input_config": {"path": str(config_path.resolve()), "sha256": config_hash},
        "input_surface": None if args.demo else {"path": str(args.surface.resolve()), "sha256": original_surface_hash},
        "artifacts_sha256": {name: _sha256(output / name) for name in artifacts},
        "source_sha256": {name: _sha256(ROOT / name) for name in (
            "mhdlab/planar_release.py", "mhdlab/planar_release_io.py", "mhdlab/thermal_release.py",
            "mhdlab/surface_sources.py", "mhdlab/flash_surface.py", "mhdlab/flash_planar.py",
            "mhdlab/hdf_integrity.py", "tools/run_planar_release.py",
        )},
        "limitations": ["prescribed heating; wall-energy debit has no temperature feedback",
                        "x/z/time convergence required before interpreting physical release",
                        "no gap transport or radiation calculation is invoked"],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return {"output_dir": str(output), "intervals": summary["intervals"], "released_particles": released.tolist(),
            "source_kind": surface.source_attributes["source_kind"], "parameter_status": config["parameter_status"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--demo", action="store_true")
    source.add_argument("--surface", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("runs/planar_release_reference"))
    parser.add_argument("--subdivisions", type=int)
    parser.add_argument("--boundary-proxy-reference")
    parser.add_argument("--include-full-face-uniform", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.surface is not None and args.config is None:
        parser.error("--surface requires an explicit --config with supplied release parameters")
    try:
        report = run(args)
    except (ValueError, OSError, KeyError, TypeError, OverflowError) as exc:
        parser.error(str(exc))
    print(json.dumps(report, allow_nan=False))


if __name__ == "__main__":
    main()
