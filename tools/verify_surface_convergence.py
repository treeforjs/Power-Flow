#!/usr/bin/env python3
"""Measure time and width convergence for the manufactured surface source.

Temporal subdivision converges release on the fixed, linearly interpolated
temperature snapshots. It does not refine the snapshots themselves. The width
study refines analytic bin-averaged temperatures at a fixed temporal resolution.
All comparisons use independent finite-resolution references, not exact
solutions or experimental calibration.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.run_surface_release import load_parameters, manufactured_surface
from mhdlab.surface_release_history import evolve_surface_release


REPORT_SCHEMA = "mhdlab-manufactured-surface-convergence-v1"
_QUANTITIES = {
    "released_particles": ("released_m2", "1"),
    "translational_energy_j": ("translational_energy_j_m2", "J"),
    "internal_energy_j": ("internal_energy_j_m2", "J"),
    "binding_energy_j": ("binding_energy_j_m2", "J"),
    "total_wall_energy_debit_j": ("total_wall_energy_debit_j_m2", "J"),
}
_SOURCE_FILES = (
    "tools/verify_surface_convergence.py", "tools/run_surface_release.py",
    "mhdlab/constants.py", "mhdlab/flash_surface.py", "mhdlab/thermal_release.py",
    "mhdlab/surface_sources.py", "mhdlab/surface_release_history.py",
)
_ROUND = 128 * np.finfo(float).eps


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _levels(values, reference, name):
    values = tuple(_positive_integer(value, name) for value in values)
    _positive_integer(reference, f"{name} reference")
    if len(values) < 2 or any(right <= left for left, right in zip(values, values[1:])):
        raise ValueError(f"{name} requires at least two strictly increasing resolutions")
    if reference <= values[-1]:
        raise ValueError(f"{name} reference must be finer than every compared resolution")
    return values


def _surface_array_digest(surface):
    digest = hashlib.sha256()
    for name in ("time_s", "x_edges_m", "temperature_k"):
        array = np.asarray(getattr(surface, name), dtype="<f8", order="C")
        digest.update(json.dumps({"name": name, "dtype": "<f8", "shape": list(array.shape)},
                                 sort_keys=True, separators=(",", ":")).encode("utf-8"))
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _case(config_path, width_bins, input_time_intervals, subdivisions):
    surface = manufactured_surface(width_bins=width_bins, time_intervals=input_time_intervals)
    config, pools, inventory = load_parameters(config_path, width_bins)
    history = evolve_surface_release(surface, inventory, pools,
                                     parameter_provenance=config["parameter_provenance"],
                                     subdivisions_per_interval=subdivisions)
    length = config["current_direction_length_m"]
    totals = []
    for pool_index, pool in enumerate(pools):
        item = {"reservoir": pool.name, "state_id": pool.state_id}
        for name, (field, _) in _QUANTITIES.items():
            item[name] = math.fsum(float(value) * length for value in history.width_integral(field)[:, pool_index])
        totals.append(item)
    return {
        "parameters": {"width_bins": width_bins, "input_time_intervals": input_time_intervals,
                       "subdivisions_per_input_interval": subdivisions,
                       "release_intervals": input_time_intervals * subdivisions,
                       "current_direction_length_m": length,
                       "foil_width_m": float(surface.x_edges_m[-1] - surface.x_edges_m[0]),
                       "time_start_s": float(surface.time_s[0]), "time_stop_s": float(surface.time_s[-1]),
                       "release_time_step_s": float((surface.time_s[-1] - surface.time_s[0]) /
                                                    (input_time_intervals * subdivisions)),
                       "width_bin_size_m": float((surface.x_edges_m[-1] - surface.x_edges_m[0]) / width_bins)},
        "surface_attributes": dict(surface.attributes),
        "surface_arrays_sha256": _surface_array_digest(surface),
        "surface_array_hash_encoding": "time_s,x_edges_m,temperature_k; sorted compact JSON header(name,dtype,shape) followed by contiguous little-endian float64 bytes for each array",
        "totals_per_reservoir": totals,
    }


def _finite_reference_errors(case, reference):
    result = []
    for values, targets in zip(case["totals_per_reservoir"], reference["totals_per_reservoir"]):
        errors = {}
        for name in _QUANTITIES:
            delta = math.fsum((values[name], -targets[name]))
            errors[name] = {"signed_error": delta, "absolute_error": abs(delta),
                            "relative_error": abs(delta) / abs(targets[name]) if targets[name] != 0 else None,
                            "reference_value": targets[name],
                            "relative_error_status": "defined" if targets[name] != 0 else "zero_reference_undefined"}
        result.append({"reservoir": values["reservoir"], "quantity_errors": errors})
    return result


def _self_orders(cases, resolutions):
    result = []
    for start in range(len(cases) - 2):
        triple = resolutions[start:start + 3]
        ratio1, ratio2 = triple[1] / triple[0], triple[2] / triple[1]
        equal_ratio = math.isclose(ratio1, ratio2, rel_tol=1e-12, abs_tol=0)
        per_pool = []
        for pool in range(len(cases[start]["totals_per_reservoir"])):
            records = [case["totals_per_reservoir"][pool] for case in cases[start:start + 3]]
            orders = {}
            for name in _QUANTITIES:
                values = [record[name] for record in records]
                first = math.fsum((values[0], -values[1]))
                second = math.fsum((values[1], -values[2]))
                allowance = _ROUND * max(abs(value) for value in values)
                order, status = None, "defined"
                if not equal_ratio:
                    status = "unequal_refinement_ratios"
                elif all(value == 0 for value in values):
                    status = "identically_zero_at_sampled_resolutions"
                elif abs(first) <= allowance or abs(second) <= allowance:
                    status = "roundoff_or_identical_totals"
                elif math.copysign(1, first) != math.copysign(1, second):
                    status = "nonmonotonic_solution_differences"
                else:
                    order = (math.log(abs(first)) - math.log(abs(second))) / math.log(ratio1)
                orders[name] = {"order": order, "status": status,
                                "coarse_minus_medium": first, "medium_minus_fine": second}
            per_pool.append({"reservoir": records[0]["reservoir"], "quantity_orders": orders})
        result.append({"coarse_resolution": triple[0], "medium_resolution": triple[1],
                       "fine_resolution": triple[2], "orders_per_reservoir": per_pool})
    return result


def _diagnostics(cases, orders, label):
    warnings = []
    observations = []
    for pool, record in enumerate(cases[0]["totals_per_reservoir"]):
        for quantity in _QUANTITIES:
            if all(case["totals_per_reservoir"][pool][quantity] == 0 for case in cases):
                observations.append({"reservoir": record["reservoir"], "quantity": quantity,
                                     "status": "identically_zero_at_sampled_resolutions"})
                continue
            sampled = [item["orders_per_reservoir"][pool]["quantity_orders"][quantity] for item in orders]
            available = [item["order"] for item in sampled if item["order"] is not None]
            status = "consistent_with_second_order_at_sampled_levels"
            if any(item["status"] == "nonmonotonic_solution_differences" for item in sampled):
                status = "nonmonotonic_solution_differences"
            elif len(available) < 2:
                status = "insufficient_resolved_orders_to_assess_asymptotic_regime"
            elif abs(available[-1] - available[-2]) > .5 or not 1.5 <= available[-1] <= 2.5:
                status = "second_order_asymptotic_regime_not_established"
            if status != "consistent_with_second_order_at_sampled_levels":
                warnings.append(f"{label}: {record['reservoir']} / {quantity}: {status}")
            observations.append({"reservoir": record["reservoir"], "quantity": quantity,
                                 "status": status, "last_defined_order": available[-1] if available else None})
    return {"observations": observations, "warnings": warnings,
            "interpretation": "Heuristic diagnostics from sampled solution differences; no physical validation or proof of asymptotic convergence."}


def _study(label, levels, reference_level, fixed, config_path, input_time_intervals):
    def evaluate(resolution):
        width = fixed if label == "temporal_subdivision" else resolution
        subdivisions = resolution if label == "temporal_subdivision" else fixed
        case = _case(config_path, width, input_time_intervals, subdivisions)
        case["resolution"] = resolution
        return case
    cases = [evaluate(level) for level in levels]
    reference = evaluate(reference_level)
    for case in cases:
        case["errors_per_reservoir"] = _finite_reference_errors(case, reference)
    all_cases, all_levels = [*cases, reference], [*levels, reference_level]
    orders = _self_orders(all_cases, all_levels)
    return {"reference_kind": "independent_finite_resolution_not_exact",
            "reference": reference, "cases": cases, "self_convergence_orders": orders,
            "order_definition": "log(abs(Q_coarse-Q_medium)/abs(Q_medium-Q_fine))/log(refinement_ratio); only equal-ratio triples with resolved same-sign differences",
            "diagnostics": _diagnostics(all_cases, orders, label)}


def convergence_report(config_path, *, time_subdivisions=(1, 2, 4, 8), time_reference_subdivisions=16,
                       time_width_bins=80, input_time_intervals=40, width_bins=(20, 40, 80, 160),
                       width_reference_bins=320, width_time_subdivisions=8):
    """Compute separate time and width studies on one manufactured recipe."""
    time_levels = _levels(time_subdivisions, time_reference_subdivisions, "time subdivisions")
    width_levels = _levels(width_bins, width_reference_bins, "width bins")
    for name, value in (("time width bins", time_width_bins), ("input time intervals", input_time_intervals),
                        ("width time subdivisions", width_time_subdivisions)):
        _positive_integer(value, name)
    path = Path(config_path).resolve()
    config_bytes = path.read_bytes()
    config = json.loads(config_bytes)
    if not isinstance(config, dict):
        raise ValueError("configuration must be a JSON object")
    if config.get("parameter_status") != "manufactured":
        raise ValueError("convergence utility requires parameter_status='manufactured'; it cannot qualify physical desorption parameters")
    reservoirs = config.get("reservoirs")
    if not isinstance(reservoirs, list) or not reservoirs or any(not isinstance(pool, dict) for pool in reservoirs):
        raise ValueError("configuration reservoirs must be a nonempty list of records")
    if any(not isinstance(pool.get("initial_coverage_m2"), (int, float))
                                           or isinstance(pool.get("initial_coverage_m2"), bool)
                                           for pool in reservoirs):
        raise ValueError("manufactured width refinement requires scalar initial_coverage_m2 for every reservoir")
    # Use the exact shared parameter loader; retain both supplied and validated
    # values so the report can be reconstructed without assumed defaults.
    all_grids = [(time_width_bins, level) for level in (*time_levels, time_reference_subdivisions)]
    all_grids += [(level, width_time_subdivisions) for level in (*width_levels, width_reference_bins)]
    if any(input_time_intervals * subdivisions > 100_000 or
           input_time_intervals * subdivisions * len(reservoirs) * width > 20_000_000
           for width, subdivisions in all_grids):
        raise ValueError("requested sweep exceeds the reference release memory limit")
    _, pools, _ = load_parameters(path, time_width_bins)
    source_hashes = {name: _sha256(ROOT / name) for name in _SOURCE_FILES}
    temporal = _study("temporal_subdivision", time_levels, time_reference_subdivisions,
                      time_width_bins, path, input_time_intervals)
    spatial = _study("width_refinement", width_levels, width_reference_bins,
                     width_time_subdivisions, path, input_time_intervals)
    if hashlib.sha256(path.read_bytes()).hexdigest() != hashlib.sha256(config_bytes).hexdigest():
        raise ValueError("configuration changed during convergence study")
    if any(_sha256(ROOT / name) != digest for name, digest in source_hashes.items()):
        raise ValueError("source code changed during convergence study")
    return {
        "schema": REPORT_SCHEMA, "source_kind": "synthetic", "parameter_status": "manufactured",
        "scope": "Numerical refinement of a manufactured prescribed-temperature, finite-inventory surface source",
        "configuration": {"path": str(path), "sha256": hashlib.sha256(config_bytes).hexdigest(), "values": config},
        "validated_reservoirs": [asdict(pool) for pool in pools],
        "source_files_sha256": source_hashes,
        "software": {"python": platform.python_version(), "numpy": np.__version__},
        "quantities": {name: {"units": unit, "width_and_current_integrated": True}
                       for name, (_, unit) in _QUANTITIES.items()},
        "temporal_subdivision": temporal, "width_refinement": spatial,
        "warnings": [*temporal["diagnostics"]["warnings"], *spatial["diagnostics"]["warnings"]],
        "limitations": [
            "Finite references retain discretization errors; reported relative errors do not estimate errors against an exact solution.",
            "Temporal subdivision holds the supplied temperature snapshots and their linear interpolation fixed. Refining independent heating snapshots is a separate required study.",
            "Width refinement applies nonlinear release to analytic bin-mean temperatures at fixed temporal resolution; it does not establish convergence of raw heating calculations or current-direction structure.",
            "Observed orders use three-level solution differences, not biased ratios of errors against a nearby finite reference. Unequal-ratio, zero, nonmonotonic, or roundoff-limited triples have no defined order.",
            "All heating and reservoir parameters are manufactured; this report is neither an experimental comparison nor a calibration of physical desorption.",
        ],
    }


def save_plot(report, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 2, figsize=(11, 8))
    for column, (study_name, grid_name, scale, x_label) in enumerate((
        ("temporal_subdivision", "release_time_step_s", 1e9, "Release step [ns]"),
        ("width_refinement", "width_bin_size_m", 1e6, "Width bin [µm]"),
    )):
        study = report[study_name]
        for row, quantity in enumerate(("released_particles", "total_wall_energy_debit_j")):
            axis = axes[row, column]
            for pool, reservoir in enumerate(study["reference"]["totals_per_reservoir"]):
                pairs = [(case["parameters"][grid_name] * scale,
                          case["errors_per_reservoir"][pool]["quantity_errors"][quantity]["relative_error"])
                         for case in study["cases"]]
                pairs = [(spacing, error) for spacing, error in pairs if error is not None and error > 0]
                if pairs:
                    x, y = zip(*pairs)
                    axis.loglog(x, y, "o-", label=reservoir["reservoir"])
            axis.set_xlabel(x_label)
            axis.set_ylabel(f"Relative error: {'released number' if row == 0 else 'wall energy'}")
            axis.grid(True, which="both", alpha=.25)
            if axis.lines:
                axis.legend()
            else:
                axis.text(.5, .5, "No positive defined errors", transform=axis.transAxes, ha="center")
            axis.set_title("Time subdivision" if column == 0 else "Width refinement")
    figure.suptitle("Manufactured surface source: errors against finite references")
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/surface_release_manufactured.json")
    parser.add_argument("--time-subdivisions", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--time-reference-subdivisions", type=int, default=16)
    parser.add_argument("--time-width-bins", type=int, default=80)
    parser.add_argument("--input-time-intervals", type=int, default=40)
    parser.add_argument("--width-bins", type=int, nargs="+", default=[20, 40, 80, 160])
    parser.add_argument("--width-reference-bins", type=int, default=320)
    parser.add_argument("--width-time-subdivisions", type=int, default=8)
    parser.add_argument("--output-dir", type=Path, default=Path("runs/surface_convergence"))
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args(argv)
    try:
        output = args.output_dir.resolve()
        if args.config.resolve() in {output / "report.json", output / "convergence.png"}:
            raise ValueError("output destinations must not overwrite the source configuration")
        report = convergence_report(args.config, time_subdivisions=args.time_subdivisions,
                                    time_reference_subdivisions=args.time_reference_subdivisions,
                                    time_width_bins=args.time_width_bins, input_time_intervals=args.input_time_intervals,
                                    width_bins=args.width_bins, width_reference_bins=args.width_reference_bins,
                                    width_time_subdivisions=args.width_time_subdivisions)
        output.mkdir(parents=True, exist_ok=True)
        (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        if args.plot:
            save_plot(report, output / "convergence.png")
    except (OSError, ValueError, TypeError, KeyError, OverflowError) as error:
        parser.error(str(error))
    print(json.dumps({"report_path": str(output / "report.json"), "source_kind": "synthetic",
                      "warning_count": len(report["warnings"]), "warnings": report["warnings"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
