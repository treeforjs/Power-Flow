#!/usr/bin/env python3
"""Separate manufactured planar source mesh, snapshot and integration changes.

The finest result in each family is finite. Three-level self differences, not
ratios of errors to that finite reference, estimate an observed order. No random
packets, gap evolution or physical parameter fit is used.
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

from mhdlab.emission_time_moments import (
    conditional_time_moments, compare_conditional_time_moments, QUADRATURE_IDENTITY,
)
from mhdlab.source_convergence import (
    FIXTURE_IDENTITY, manufactured_planar_surface, source_observables,
    compare_observables, observed_self_order,
)
from mhdlab.temkin_emission_time import FrozenTemkinEmission
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_release import TemkinReservoir, LAW_IDENTITY


FAMILIES = {
    "width": "width_bins", "current": "current_bins",
    "snapshot": "snapshot_intervals", "integration": "subdivisions_per_interval",
}
CODE_PATHS = (
    "tools/verify_temkin_planar_convergence.py", "mhdlab/source_convergence.py",
    "mhdlab/emission_time_moments.py", "mhdlab/temkin_emission_time.py",
    "mhdlab/temkin_planar.py", "mhdlab/temkin_release.py", "mhdlab/temkin_history.py",
    "mhdlab/thermal_release.py", "mhdlab/planar_release.py", "mhdlab/flash_surface.py",
    "mhdlab/surface_release_history.py", "mhdlab/surface_sources.py", "mhdlab/constants.py",
)
LIMITS = [
    "Manufactured temperatures, independent coverage pools, rates, energies and exact-state masses; no H2O/304 calibration.",
    "Each family changes one resolution while the others remain finite; no continuum or experimental truth is supplied.",
    "Snapshot interpolation and frozen-midpoint integration are distinct approximations.",
    "Conditional event-time estimates use finite-order quadrature; changes are not certified error bounds or physical uncertainties.",
    "No particle sampling noise, transport evolution, thermal feedback, moving wall or actual FLASH execution is included.",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def array_record(value):
    """Hash explicit shape/dtype and canonical little-endian C-order bytes."""
    value = np.asarray(value, dtype="<f8", order="C")
    header = {"shape": list(value.shape), "dtype": "<f8", "order": "C"}
    data = json.dumps(header, sort_keys=True, separators=(",", ":")).encode() + b"\n" + value.tobytes()
    return dict(header, sha256=hashlib.sha256(data).hexdigest())


def _integer(value, name, minimum=1, maximum=100000):
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum},{maximum}]")
    return value


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive number")
    return float(value)


def _unique_object(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"duplicate JSON key: {key}")
        value[key] = item
    return value


def load_config(path):
    config = json.loads(path.read_bytes(), object_pairs_hook=_unique_object)
    expected = {"schema", "parameter_status", "parameter_provenance", "duration_s", "width_m",
                "current_length_m", "base_resolution", "refinements", "selected_metrics", "reservoirs",
                "state_masses_kg", "include_time_moments", "time_quadrature_order", "event_time_quadrature"}
    if (not isinstance(config, dict) or set(config) != expected
            or config["schema"] != "mhdlab-temkin-planar-convergence-config-v1"):
        raise ValueError("configuration must contain exactly the mhdlab-temkin-planar-convergence-config-v1 fields")
    if config["parameter_status"] != "manufactured":
        raise ValueError("this fixed analytic fixture requires parameter_status=manufactured")
    if not isinstance(config["parameter_provenance"], str) or not config["parameter_provenance"].strip():
        raise ValueError("parameter_provenance must be nonempty")
    for name in ("duration_s", "width_m", "current_length_m"):
        _positive(config[name], name)
    for name in ("base_resolution", "refinements"):
        if not isinstance(config[name], dict) or set(config[name]) != set(FAMILIES.values()):
            raise ValueError(f"{name} must declare all four independent resolutions")
    for key in FAMILIES.values():
        _integer(config["base_resolution"][key], key)
        study = config["refinements"][key]
        if not isinstance(study, dict) or set(study) != {"levels", "reference"}:
            raise ValueError("each refinement requires levels and reference")
        levels = study["levels"]
        if not isinstance(levels, list) or len(levels) < 3:
            raise ValueError("each refinement needs at least three increasing geometric levels")
        for level in levels + [study["reference"]]:
            _integer(level, key)
        ratios = [b / a for a, b in zip(levels[:-1], levels[1:])]
        if ratios[0] <= 1 or any(not math.isclose(ratio, ratios[0], rel_tol=1e-12) for ratio in ratios):
            raise ValueError("refinement levels must have one increasing geometric ratio")
        if study["reference"] <= levels[-1]:
            raise ValueError("finite reference must be finer than all study levels")
    selected = config["selected_metrics"]
    if not isinstance(selected, list) or not selected or any(not isinstance(v, str) for v in selected) or len(set(selected)) != len(selected):
        raise ValueError("selected_metrics must be a nonempty list of unique metric names")
    if not isinstance(config["include_time_moments"], bool):
        raise ValueError("include_time_moments must be boolean")
    _integer(config["time_quadrature_order"], "time_quadrature_order", 2, 4096)
    records = config["reservoirs"]
    if not isinstance(records, list) or not records:
        raise ValueError("reservoirs must be a nonempty list")
    pools, initial = [], []
    for record in records:
        supplied = dict(record)
        value = supplied.pop("initial_coverage_m2")
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError("initial_coverage_m2 must be a finite nonnegative scalar; mesh-dependent arrays need an explicit remapping contract")
        pool = TemkinReservoir(**supplied)
        if value > pool.max_coverage_m2:
            raise ValueError("initial coverage exceeds the declared law domain")
        pools.append(pool)
        initial.append(float(value))
    if len({pool.name for pool in pools}) != len(pools):
        raise ValueError("reservoir names must be unique")
    if not isinstance(config["state_masses_kg"], dict) or set(config["state_masses_kg"]) != {p.state_id for p in pools}:
        raise ValueError("state_masses_kg must supply exactly the declared state IDs")
    for value in config["state_masses_kg"].values():
        _positive(value, "exact-state mass")
    event = config["event_time_quadrature"]
    if not isinstance(event, dict) or set(event) != {"orders", "reference_order", "component"}:
        raise ValueError("event_time_quadrature requires orders, reference_order and component")
    orders = event["orders"]
    if not isinstance(orders, list) or len(orders) < 2:
        raise ValueError("event quadrature needs at least two orders")
    for order in orders + [event["reference_order"]]:
        _integer(order, "quadrature order", 2, 4096)
    if any(b <= a for a, b in zip(orders, orders[1:])) or event["reference_order"] <= orders[-1]:
        raise ValueError("event quadrature orders and finite reference must be strictly increasing")
    component = event["component"]
    if not isinstance(component, dict) or set(component) != {"reservoir_index", "initial_coverage_m2", "temperature_k", "dt_s"}:
        raise ValueError("event component fields are invalid")
    _integer(component["reservoir_index"], "reservoir_index", 0, len(pools) - 1)
    for name in ("initial_coverage_m2", "temperature_k", "dt_s"):
        _positive(component[name], name)
    if component["initial_coverage_m2"] > pools[component["reservoir_index"]].max_coverage_m2:
        raise ValueError("event component coverage exceeds the declared domain")
    return config, tuple(pools), initial


def _selected(observables, selected):
    if not set(selected) <= set(observables["metrics"]):
        raise ValueError("selected_metrics contains an unavailable source observable")
    return {"metrics": {key: observables["metrics"][key] for key in selected},
            "units": {key: observables["units"][key] for key in selected}}


def run_study(config, pools, initial, families):
    cases, cache, studies = {}, {}, {}

    def source_case(resolution):
        cache_key = tuple(resolution[key] for key in FAMILIES.values())
        if cache_key in cache:
            return cache[cache_key]
        case_id = "x%d-z%d-s%d-i%d" % cache_key
        print(f"source case {case_id}", file=sys.stderr, flush=True)
        surface = manufactured_planar_surface(**{k: resolution[k] for k in ("width_bins", "current_bins", "snapshot_intervals")},
                                             **{k: config[k] for k in ("duration_s", "width_m", "current_length_m")})
        coverage = np.broadcast_to(np.asarray(initial)[:, None, None],
                                   (len(pools), resolution["current_bins"], resolution["width_bins"]))
        history = evolve_temkin_planar_release(surface, coverage, pools,
                                               parameter_provenance=config["parameter_provenance"],
                                               subdivisions_per_interval=resolution["subdivisions_per_interval"])
        observables = source_observables(history, state_masses_kg=config["state_masses_kg"],
                                         include_time_moments=config["include_time_moments"],
                                         time_quadrature_order=config["time_quadrature_order"])
        _selected(observables, config["selected_metrics"])
        source_names = ("time_edges_s", "x_edges_m", "z_edges_m", "temperature_k", "coverage_m2", "released_m2",
                        "translational_energy_j_m2", "internal_energy_j_m2", "binding_energy_j_m2", "total_wall_energy_debit_j_m2")
        metrics = observables["metrics"]
        balance = [math.fsum((a, -b, -c)) for a, b, c in zip(metrics["initial_number"], metrics["released_number"], metrics["remaining_number"])]
        cases[case_id] = {
            "resolution": dict(resolution), "observables": observables,
            "input_arrays": {name: array_record(getattr(surface, name)) for name in ("time_s", "x_edges_m", "z_edges_m", "temperature_k")},
            "source_arrays": {name: array_record(getattr(history, name)) for name in source_names},
            "max_abs_inversion_log_time_residual": float(np.max(np.abs(history.inversion_log_time_residual))),
            "global_inventory_balance_residual_number": balance,
            "global_inventory_balance_relative": [value / a if a else None for value, a in zip(balance, metrics["initial_number"])],
        }
        cache[cache_key] = case_id
        return case_id

    for family in families:
        if family == "event_time":
            event = config["event_time_quadrature"]
            component = dict(event["component"])
            emission = FrozenTemkinEmission(pools[component.pop("reservoir_index")], **component)
            moments = {str(order): conditional_time_moments(emission, quadrature_order=order)
                       for order in event["orders"] + [event["reference_order"]]}
            estimates = {key: asdict(value) for key, value in moments.items()}
            names = ("mean_offset_s", "second_raw_moment_s2", "variance_offset_s2")
            target = estimates[str(event["reference_order"])]
            comparisons = {}
            for order in event["orders"]:
                value = estimates[str(order)]
                comparisons[str(order)] = compare_observables({k: [value[k]] for k in names}, {k: [target[k]] for k in names})
            studies[family] = {"component": event["component"], "orders": event["orders"],
                               "finite_reference_order": event["reference_order"], "estimates": estimates,
                               "comparisons": comparisons, "quadrature_identity": QUADRATURE_IDENTITY,
                               "log_moment_comparisons": {
                                   str(order): {key: asdict(value) for key, value in compare_conditional_time_moments(
                                       moments[str(order)], moments[str(event["reference_order"]) ]).items()}
                                   for order in event["orders"]},
                               "raw_comparison_qualification": "Raw SI zeros do not establish agreement when underflow flags are set; use log_moment_comparisons, whose relative changes retain the positive log moments.",
                               "qualification": "Finite-order estimates only; no quadrature convergence order or certified bound is assumed. Log moments and underflow flags are retained."}
            continue
        key = FAMILIES[family]
        refinement = config["refinements"][key]
        ids = []
        for level in refinement["levels"] + [refinement["reference"]]:
            resolution = dict(config["base_resolution"], **{key: level})
            ids.append(source_case(resolution))
        selected = [_selected(cases[case_id]["observables"], config["selected_metrics"]) for case_id in ids]
        comparisons = {case_id: compare_observables(obs, selected[-1]) for case_id, obs in zip(ids[:-1], selected[:-1])}
        orders = []
        levels = refinement["levels"]
        for index in range(len(levels) - 2):
            order_metrics = {}
            for name in config["selected_metrics"]:
                values = [selected[j]["metrics"][name] for j in range(index, index + 3)]
                if any(any(v is None for v in value) for value in values):
                    order_metrics[name] = {"observed_order": None, "warnings": ["normalized metric unavailable for an empty released population"]}
                else:
                    order_metrics[name] = observed_self_order(*values, refinement_ratio=levels[index + 1] / levels[index])
            orders.append({"case_ids": ids[index:index + 3], "metrics": order_metrics})
        studies[family] = {"varied_resolution": key, "fixed_resolution": {k: v for k, v in config["base_resolution"].items() if k != key},
                           "levels": levels, "case_ids": ids[:-1], "finite_reference_resolution": refinement["reference"],
                           "finite_reference_case_id": ids[-1], "comparisons": comparisons, "self_orders": orders}
    return cases, studies


def save_plot(report, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for ax, (family, study) in zip(axes.flat, report["studies"].items()):
        event = family == "event_time"
        levels = study["orders"] if event else study["levels"]
        ids = [str(value) for value in levels] if event else study["case_ids"]
        names = list(next(iter(study["comparisons"].values()))["metrics"])
        for name in names:
            values = ([study["log_moment_comparisons"][case_id][name]["relative_change"] for case_id in ids]
                      if event else [study["comparisons"][case_id]["metrics"][name]["relative_l2"] for case_id in ids])
            valid = [(x, y) for x, y in zip(levels, values) if y is not None and y > 0]
            if valid:
                ax.loglog(*zip(*valid), marker="o", label=name)
        ax.set(title=family, xlabel="quadrature order" if event else study["varied_resolution"], ylabel="relative difference to finite reference")
        ax.grid(True, which="both", alpha=.25)
        ax.legend(fontsize=6)
    for ax in list(axes.flat)[len(report["studies"]):]:
        ax.set_visible(False)
    fig.suptitle("Manufactured planar release: independent refinement families")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/temkin_convergence_manufactured.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--families", nargs="+", choices=[*FAMILIES, "event_time"], default=[*FAMILIES, "event_time"])
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        input_hash = sha(args.config)
        code_hashes = {name: sha(ROOT / name) for name in CODE_PATHS}
        config, pools, initial = load_config(args.config)
        families = list(dict.fromkeys(args.families))
        cases, studies = run_study(config, pools, initial, families)
        report = {"schema": "mhdlab-temkin-planar-convergence-report-v1", "parameter_status": "manufactured",
                  "fixture_identity": FIXTURE_IDENTITY, "law_identity": LAW_IDENTITY,
                  "reference_kind": "finite_finest_resolution_reference", "config": config,
                  "input_config_sha256": input_hash, "source_sha256": code_hashes,
                  "unique_source_cases": len(cases), "cases": cases, "studies": studies, "limitations": LIMITS}
        if sha(args.config) != input_hash or any(sha(ROOT / name) != value for name, value in code_hashes.items()):
            raise ValueError("input or implementation changed during the study")
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".source-convergence-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            with (output / "metrics.csv").open("w", newline="") as stream:
                writer = csv.writer(stream, lineterminator="\n")
                writer.writerow(["case_id", *FAMILIES.values(), "metric", "reservoir", "unit", "value"])
                for case_id, case in cases.items():
                    for name, values in case["observables"]["metrics"].items():
                        for pool, value in zip(pools, values):
                            writer.writerow([case_id, *(case["resolution"][k] for k in FAMILIES.values()), name,
                                             pool.name, case["observables"]["units"][name], value])
            if args.plot:
                save_plot(report, output / "convergence.png")
            report["output_sha256"] = {p.name: sha(p) for p in output.iterdir() if p.is_file()}
            (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during the study; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as exc:
        parser.error(str(exc))
    print(json.dumps({"output_dir": str(args.output_dir), "unique_source_cases": len(cases),
                      "families": families, "parameter_status": "manufactured"}, allow_nan=False))


if __name__ == "__main__":
    main()
