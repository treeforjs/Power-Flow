"""Convert a validated surface-temperature history to conservative release.

Run --demo for explicitly manufactured heating and binding parameters, or
supply --surface and --config together for an externally prepared history.
This workflow invokes no gap, Boltzmann, CR or FLASH executable.
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

from mhdlab.flash_surface import (
    FlashSurfaceHistory, as_boundary_temperature_proxy, as_current_uniform_temperature_proxy,
    load_surface_history, write_surface_history,
)
from mhdlab.surface_release_history import evolve_surface_release
from mhdlab.surface_release_io import write_release_history
from mhdlab.thermal_release import BindingReservoir, SurfaceInventory


def manufactured_surface(*, width_bins: int = 80, time_intervals: int = 40) -> FlashSurfaceHistory:
    """Analytic bin-averaged hot-edge temperature fixture, with no shot claim."""
    edges = np.linspace(-0.0025, 0.0025, width_bins + 1)
    time = np.linspace(0, 120e-9, time_intervals + 1)
    scale = 0.00025
    erf = np.vectorize(math.erf, otypes=[float])
    left = 0.5 * scale * math.sqrt(math.pi) * np.diff(erf((edges - edges[0]) / scale)) / np.diff(edges)
    right = -0.5 * scale * math.sqrt(math.pi) * np.diff(erf((edges[-1] - edges) / scale)) / np.diff(edges)
    ramp = (time / time[-1]) ** 2
    temperature = 300 + 300 * ramp[:, None] + 1400 * ramp[:, None] * (left + right)
    return FlashSurfaceHistory(time, edges, temperature, "synthetic", "manufactured-hot-edge-reference",
                               "tools/run_surface_release.py:manufactured_surface analytic bin averages",
                               "manufactured stationary gap-facing face")


def load_parameters(path: Path, width_bins: int):
    config = json.loads(path.read_text())
    if config.get("schema") != "mhdlab-surface-release-config-v1":
        raise ValueError("unsupported release configuration schema")
    if config.get("parameter_status") not in {"manufactured", "sourced", "unvalidated"}:
        raise ValueError("parameter_status must be manufactured, sourced, or unvalidated")
    if not isinstance(config.get("parameter_provenance"), str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must describe the supplied desorption parameters")
    length = config.get("current_direction_length_m")
    if isinstance(length, bool) or not isinstance(length, (int, float)) or not np.isfinite(length) or length <= 0:
        raise ValueError("current_direction_length_m must be finite and positive")
    pools, coverage = [], []
    for supplied in config["reservoirs"]:
        supplied = dict(supplied)
        initial = supplied.pop("initial_coverage_m2")
        if isinstance(initial, bool):
            raise ValueError("initial_coverage_m2 must be numeric")
        if np.isscalar(initial):
            # SurfaceInventory validates the supplied scalar before broadcasting.
            scalar = SurfaceInventory([[initial]]).coverage_m2[0, 0]
            initial = np.full(width_bins, scalar)
        elif np.asarray(initial, dtype=object).shape != (width_bins,):
            raise ValueError("each initial_coverage_m2 must be a scalar or one value per width bin")
        pools.append(BindingReservoir(**supplied))
        coverage.append(initial)
    return config, tuple(pools), SurfaceInventory(coverage)


def json_value(value):
    if isinstance(value, dict):
        return {name: json_value(item) for name, item in value.items()}
    if isinstance(value, (np.ndarray, np.ma.MaskedArray)):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def save_plot(history, output: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    centers = 0.5 * (history.x_edges_m[:-1] + history.x_edges_m[1:]) * 1e3
    midpoint_ns = (history.time_edges_s[:-1] + 0.5 * np.diff(history.time_edges_s)) * 1e9
    flux = history.average_flux_m2_s.sum(axis=1)
    cumulative = history.released_m2.sum(axis=(0, 1))
    uniform = history.matched_uniform()["released_m2"].sum(axis=(0, 1))
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for index in sorted({0, len(midpoint_ns) // 2, len(midpoint_ns) - 1}):
        axes[0, 0].plot(centers, history.temperature_k[index], label=f"{midpoint_ns[index]:.1f} ns")
    axes[0, 0].set_ylabel("Prescribed temperature [K]")
    axes[0, 0].legend()
    mesh = axes[0, 1].pcolormesh(history.x_edges_m * 1e3, history.time_edges_s * 1e9,
                               flux, shading="flat", rasterized=True)
    axes[0, 1].set_ylabel("Time [ns]")
    fig.colorbar(mesh, ax=axes[0, 1], label="Total release flux [m⁻² s⁻¹]")
    axes[1, 0].plot(centers, cumulative, label="Width-resolved release")
    axes[1, 0].plot(centers, uniform, "--", label="Matched uniform source")
    axes[1, 0].set_ylabel("Cumulative release [m⁻²]")
    axes[1, 0].legend()
    initial = history.coverage_m2[0].sum(axis=0)
    remaining = history.coverage_m2[-1].sum(axis=0)
    axes[1, 1].plot(centers, remaining / initial, label="Remaining inventory fraction")
    axes[1, 1].set_ylim(0, 1.05)
    axes[1, 1].set_ylabel("Remaining / initial")
    for ax in axes.flat:
        ax.set_xlabel("Foil-width position [mm]")
    fig.suptitle(f"Prescribed-temperature release; source={history.surface_attributes['source_kind']}")
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--demo", action="store_true")
    source.add_argument("--surface", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--subdivisions", type=int)
    parser.add_argument("--boundary-proxy-reference", help="Explicit evidence reference for approximating boundary T with a material-plane sample")
    parser.add_argument("--current-uniformity-reference", help="Explicit evidence reference for applying a current-averaged T to a nonlinear source")
    parser.add_argument("--output-dir", type=Path, default=Path("runs/surface_release_reference"))
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.surface and args.config is None:
        parser.error("--surface requires an explicit --config with reservoir parameters")
    config_path = args.config or ROOT / "configs/surface_release_manufactured.json"
    surface = manufactured_surface() if args.demo else load_surface_history(args.surface)
    if args.boundary_proxy_reference is not None:
        surface = as_boundary_temperature_proxy(surface, reference=args.boundary_proxy_reference)
    if args.current_uniformity_reference is not None:
        surface = as_current_uniform_temperature_proxy(surface, reference=args.current_uniformity_reference)
    config, pools, inventory = load_parameters(config_path, len(surface.x_edges_m) - 1)
    subdivisions = args.subdivisions if args.subdivisions is not None else config["subdivisions_per_interval"]
    history = evolve_surface_release(surface, inventory, pools,
                                     parameter_provenance=config["parameter_provenance"],
                                     subdivisions_per_interval=subdivisions)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    write_surface_history(output / "input_surface.h5", surface)
    write_release_history(output / "surface_release.h5", history)
    length = config["current_direction_length_m"]
    widths = np.diff(history.x_edges_m)
    initial = np.sum(history.coverage_m2[0] * widths, axis=-1) * length
    remaining = np.sum(history.coverage_m2[-1] * widths, axis=-1) * length
    released = history.width_integral().sum(axis=0) * length
    debit = history.width_integral("total_wall_energy_debit_j_m2").sum(axis=0) * length
    contrast = history.release_contrast()
    summary = {
        "scope": "prescribed-temperature finite-inventory source and matched spatial-redistribution control",
        "surface_attributes": dict(history.surface_attributes),
        "parameter_status": config["parameter_status"],
        "parameter_provenance": config["parameter_provenance"],
        "current_direction_length_m": length,
        "reservoirs": [asdict(pool) for pool in pools],
        "initial_particles": initial, "remaining_particles": remaining,
        "released_particles": released, "wall_energy_debit_j": debit,
        "particle_balance_residual": [math.fsum((float(a), -float(b), -float(c)))
                                      for a, b, c in zip(initial, remaining, released)],
        "subdivisions_per_input_interval": subdivisions,
        "intervals": len(history.time_edges_s) - 1,
        "source_contrast": contrast,
        "uniform_control": "Each interval/reservoir preserves width-integrated number and energy moments; retain the source-temperature mixture for full velocity-distribution preservation",
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "surface_input_sha256": hashlib.sha256(args.surface.read_bytes()).hexdigest() if args.surface else None,
        "source_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in
                          ["mhdlab/flash_surface.py", "mhdlab/thermal_release.py", "mhdlab/surface_sources.py",
                           "mhdlab/surface_release_history.py", "mhdlab/surface_release_io.py", "tools/run_surface_release.py"]},
        "remaining_work": ["validated physical desorption parameters", "coupled wall-energy feedback",
                           "gap and field evolution", "radiative and detector projection"],
    }
    (output / "summary.json").write_text(json.dumps(json_value(summary), indent=2, allow_nan=False) + "\n")
    with (output / "integrated_source.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["time_start_s", "time_stop_s", "reservoir", "state_id", "released_particles", "wall_energy_debit_j"])
        counts = history.width_integral() * length
        energies = history.width_integral("total_wall_energy_debit_j_m2") * length
        for interval in range(len(history.time_edges_s) - 1):
            for index, pool in enumerate(pools):
                writer.writerow([history.time_edges_s[interval], history.time_edges_s[interval + 1],
                                 pool.name, pool.state_id, counts[interval, index], energies[interval, index]])
    if args.plot:
        save_plot(history, output / "overview.png")
    print(json.dumps({"output_dir": str(output), "intervals": len(history.time_edges_s) - 1,
                      "released_particles": released.tolist(), "parameter_status": config["parameter_status"]}))


if __name__ == "__main__":
    main()
