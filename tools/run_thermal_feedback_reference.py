#!/usr/bin/env python3
"""Run explicit manufactured zero-power single-patch feedback references.

Each configured case/grid is evolved once. All coefficients are supplied in
the configuration; no material, water, state-chemistry or host model is chosen.
The maintained-temperature control belongs to the independent reference check.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CONFIG_SCHEMA = "mhdlab-thermal-feedback-manufactured-config-v1"
REPORT_SCHEMA = "mhdlab-thermal-feedback-reference-run-v1"
MAX_STEPS = 16384
MAX_TRAJECTORIES = 24
_ID = re.compile(r"[a-z][a-z0-9_]{0,63}\Z")
_METRICS = (
    "coverage_m2", "temperature_k", "released_m2", "binding_energy_j_m2",
    "internal_energy_j_m2", "translational_energy_j_m2",
    "total_wall_energy_debit_j_m2", "gas_energy_j_m2",
    "wall_thermal_change_j_m2", "adsorbate_reference_change_j_m2",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reject_constant(value):
    raise ValueError("nonfinite JSON constant: " + value)


def _json_float(token):
    """Reject numerical lexemes lost before the core can qualify them."""
    try:
        exact = Decimal(token)
        value = float(exact)
    except (InvalidOperation, OverflowError, ValueError) as exc:
        raise ValueError("unsupported JSON numeric representation") from exc
    if not exact.is_finite() or not math.isfinite(value):
        raise ValueError("JSON number overflows its finite binary64 representation")
    if exact != 0 and value == 0:
        raise ValueError("nonzero JSON number underflows its binary64 representation")
    return value


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate configuration key: " + key)
        result[key] = value
    return result


def _positive(value, name):
    if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
        raise ValueError(name + " must be an explicitly supplied positive finite number")
    return float(value)


def load_config(path: Path) -> dict:
    config = json.loads(path.read_text(), object_pairs_hook=_object,
                        parse_float=_json_float, parse_constant=_reject_constant)
    required = {"schema", "parameter_status", "provenance", "external_power_w_m2", "scope",
                "maintained_temperature_control", "refinement_case_id", "cases"}
    if type(config) is not dict or set(config) != required or config["schema"] != CONFIG_SCHEMA:
        raise ValueError("unsupported manufactured feedback configuration")
    if config["parameter_status"] != "manufactured":
        raise ValueError("this reference CLI accepts manufactured parameters only")
    for key in ("provenance", "scope", "maintained_temperature_control"):
        if type(config[key]) is not str or not config[key].strip():
            raise ValueError(key + " must be explicit nonempty text")
    power = config["external_power_w_m2"]
    if type(power) not in (int, float) or not math.isfinite(power) or power != 0:
        raise ValueError("only explicitly zero external power is supported")
    cases = config["cases"]
    if type(cases) is not list or not 1 <= len(cases) <= 6:
        raise ValueError("one to six explicit cases are supported")
    ids, run_count = set(), 0
    for case in cases:
        if type(case) is not dict or set(case) != {"case_id", "purpose", "duration_s", "steps", "problem"}:
            raise ValueError("each case must supply its complete declared input")
        identifier = case["case_id"]
        if type(identifier) is not str or not _ID.fullmatch(identifier) or identifier in ids:
            raise ValueError("case identifiers must be unique bounded lowercase identifiers")
        ids.add(identifier)
        if type(case["purpose"]) is not str or not case["purpose"].strip():
            raise ValueError("each manufactured case needs an explicit purpose")
        _positive(case["duration_s"], "duration_s")
        steps = case["steps"]
        if (type(steps) is not list or not steps or any(type(n) is not int or not 1 <= n <= MAX_STEPS for n in steps)
                or steps != sorted(set(steps))):
            raise ValueError("steps must be distinct increasing positive bounded integers")
        if type(case["problem"]) is not dict or case["problem"].get("parameter_status") != "manufactured":
            raise ValueError("each problem must explicitly declare manufactured status")
        run_count += len(steps)
    if run_count > MAX_TRAJECTORIES:
        raise ValueError("configured reference trajectory work limit exceeded")
    if config["refinement_case_id"] not in ids:
        raise ValueError("refinement case must identify a configured case")
    return config


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _flat_sample(sample):
    record = sample.to_dict()
    extensive = record.pop("extensive")
    if not isinstance(extensive, dict) or set(record) & set(extensive):
        raise ValueError("unsupported feedback sample/extensive mapping")
    flat = {**record, **extensive}
    if any(type(value) not in (int, float) or not math.isfinite(value) for value in flat.values()):
        raise ValueError("a feedback sample has an unsupported numerical representation")
    return flat


def _units(name):
    if name == "time_s":
        return "s"
    if name == "temperature_k":
        return "K"
    if name.endswith("_j_m2"):
        return "J m^-2"
    if name.endswith("_j"):
        return "J"
    if name.endswith("_number"):
        return "parent"
    if name.endswith("_m2"):
        return "parent m^-2"
    return "dimensionless"


def _write_csv(path, rows):
    fields = ("sample_index", "case_id", "steps", "time_s", "coverage_m2", "temperature_k", "released_m2")
    fields += tuple(key for key in rows[0] if key not in fields)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format(value, ".17g") if type(value) is float else value for key, value in row.items()})
    return {key: _units(key) for key in fields if key not in ("case_id", "sample_index", "steps")}


def _refinement(records, case, path):
    family = sorted((record for record in records if record["case_id"] == case["case_id"]), key=lambda row: row["steps"])
    result, previous = [], {}
    for coarse, fine in zip(family[:-1], family[1:]):
        a, b = coarse["final_sample"], fine["final_sample"]
        for metric in _METRICS:
            change = b[metric] - a[metric]
            if metric in ("coverage_m2", "released_m2"):
                budget, convention = b["released_m2"], "fine emitted areal parent inventory"
            elif metric == "temperature_k":
                budget = abs(b[metric] - case["problem"]["initial_temperature_k"])
                convention = "fine temperature change, not initial temperature"
            else:
                budget, convention = abs(b[metric]), "absolute fine channel change/debit"
            old = previous.get(metric)
            order = None
            if old is not None and old > 0 and abs(change) > 0 and fine["steps"] == 2 * coarse["steps"]:
                order = (math.log(old) - math.log(abs(change))) / math.log(2.)
            row = {"case_id": case["case_id"], "coarse_steps": coarse["steps"], "fine_steps": fine["steps"],
                   "quantity": metric, "units": _units(metric), "signed_final_change": change,
                   "absolute_final_change": abs(change), "normalizing_budget": budget,
                   "normalizing_budget_definition": convention,
                   "change_over_budget": None if budget == 0 else abs(change) / budget,
                   "zero_budget_exact_zero_change": budget == 0 and change == 0,
                   "observed_self_order": order,
                   "interpretation": "finite grid changes only; no independent error bound or accuracy certificate"}
            result.append(row)
            previous[metric] = abs(change)
    fields = tuple(result[0]) if result else ("case_id", "coarse_steps", "fine_steps", "quantity")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in result:
            writer.writerow({key: "" if value is None else format(value, ".17g") if type(value) is float else value for key, value in row.items()})
    return result


def _plots(config, saved, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(3, 2, figsize=(12, 10), layout="constrained")
    for ax, case in zip(axes.flat, config["cases"]):
        entries = [entry for entry in saved if entry["case_id"] == case["case_id"]]
        entry = max(entries, key=lambda item: item["steps"])
        rows = entry["rows"]
        ax.plot([row["time_s"] * 1e9 for row in rows], [row["temperature_k"] for row in rows], color="firebrick")
        ax.set(xlabel="Time [ns]", ylabel="Wall temperature [K]", title=case["case_id"])
        ax.grid(alpha=.2)
        twin = ax.twinx()
        n0 = case["problem"]["initial_coverage_m2"]
        if n0 > 0:
            twin.plot([row["time_s"] * 1e9 for row in rows], [row["released_m2"] / n0 for row in rows], "--", color="navy")
            twin.set_ylabel("Emitted / supplied initial inventory")
        else:
            twin.text(.45, .6, "Exactly empty inventory", transform=twin.transAxes)
            twin.set_yticks([])
    for ax in list(axes.flat)[len(config["cases"]):]:
        ax.set_visible(False)
    figure.suptitle("Manufactured zero-power feedback; thermal reservoir and emitted number", fontsize=13)
    figure.savefig(output / "trajectories.png", dpi=160)
    figure.savefig(output / "trajectories.pdf")
    plt.close(figure)
    family = [entry for entry in saved if entry["case_id"] == config["refinement_case_id"]]
    rows = max(family, key=lambda item: item["steps"])["rows"]
    figure, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    for key, label in (("binding_energy_j_m2", "Removal / adsorbate-reference increase"),
                       ("internal_energy_j_m2", "Gas internal"),
                       ("translational_energy_j_m2", "Expected gas translation"),
                       ("total_wall_energy_debit_j_m2", "Total wall debit")):
        ax.plot([row["time_s"] * 1e9 for row in rows], [row[key] for row in rows], label=label)
    ax.plot([row["time_s"] * 1e9 for row in rows], [-row["wall_thermal_change_j_m2"] for row in rows], "--", label="Thermal decrease")
    ax.set(xlabel="Time [ns]", ylabel="Energy per area [J m$^{-2}$]",
           title="Manufactured nonlinear case: separate energy channels")
    ax.legend(fontsize=8)
    ax.grid(alpha=.2)
    figure.savefig(output / "energy_channels.png", dpi=160)
    figure.savefig(output / "energy_channels.pdf")
    plt.close(figure)


def run(config_path: Path, output: Path, *, make_plots: bool):
    config = load_config(config_path)
    if output.exists():
        raise ValueError("output must be a fresh directory; no reference output is overwritten")
    sys.path.insert(0, str(ROOT))
    from mhdlab.thermal_feedback_reference import FeedbackReservoir, FeedbackProblem, evolve_feedback

    problems = []
    for case in config["cases"]:
        supplied = dict(case["problem"])
        pool = FeedbackReservoir(**supplied.pop("reservoir"))
        problems.append((case, FeedbackProblem(reservoir=pool, **supplied)))
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=output.name + ".staging-", dir=output.parent))
    records, saved = [], []
    try:
        shutil.copyfile(config_path, stage / "input_config.json")
        for case, problem in problems:
            for steps in case["steps"]:
                history = evolve_feedback(problem, duration_s=case["duration_s"], steps=steps)
                rows = [{"sample_index": index, "case_id": case["case_id"], "steps": steps, **_flat_sample(sample)}
                        for index, sample in enumerate(history.samples)]
                stem = case["case_id"] + "_steps" + format(steps, "04d")
                csv_path, history_path = stage / (stem + ".csv"), stage / (stem + ".json")
                units = _write_csv(csv_path, rows)
                _write_json(history_path, history.to_dict())
                records.append({"case_id": case["case_id"], "purpose": case["purpose"], "steps": steps,
                                "duration_s": case["duration_s"], "sample_count": len(rows),
                                "csv_file": csv_path.name, "csv_sha256": sha(csv_path),
                                "history_file": history_path.name, "history_sha256": sha(history_path),
                                "input_parameters": case["problem"], "units": units,
                                "final_sample": rows[-1],
                                "max_absolute_residuals": {key: max(abs(row[key]) for row in rows) for key in rows[0] if "residual" in key}})
                saved.append({"case_id": case["case_id"], "steps": steps, "rows": rows})
        refinement_case = next(case for case in config["cases"] if case["case_id"] == config["refinement_case_id"])
        refinement = _refinement(records, refinement_case, stage / "temporal_refinement.csv")
        if make_plots:
            _plots(config, saved, stage)
        code = {"tools/run_thermal_feedback_reference.py": sha(Path(__file__)),
                "mhdlab/thermal_feedback_reference.py": sha(ROOT / "mhdlab/thermal_feedback_reference.py")}
        artifacts = {path.name: sha(path) for path in stage.iterdir() if path.is_file()}
        report = {"schema": REPORT_SCHEMA, "status": "manufactured_reference_only",
                  "config_original_path": str(config_path.resolve()), "config_sha256": sha(config_path),
                  "code_sha256": code, "external_power_w_m2": config["external_power_w_m2"],
                  "parameter_status": "manufactured", "scope": config["scope"],
                  "maintained_temperature_control": config["maintained_temperature_control"],
                  "case_count": len(config["cases"]), "numerical_trajectory_count": len(records),
                  "trajectories": records, "temporal_refinement": refinement,
                  "artifact_path_base": "directory containing report.json",
                  "artifact_sha256": artifacts,
                  "qualifications": {"zero_external_power": True, "explicit_inputs_only": True,
                                     "fixed_kinetic_reference": True, "expected_flux_translation_2kBT": True,
                                     "independent_reference_error_computed_here": False,
                                     "old_source_or_recovery_producer_called": False,
                                     "physical_parameters_or_water_inferred": False,
                                     "actual_FLASH_or_host_coupling": False,
                                     "floor_clamp_or_unsupported_domain_fallback": False},
                  "limits": ["Direct fixed-step manufactured feedback trajectories only; no physical material or parent parameter validation.",
                             "Separate thermal, adsorbate reference, gas internal and continuum expected translational channels remain explicit; signed residuals are not repaired.",
                             "Finite temporal changes and observed self-orders are descriptive, not certified error bounds.",
                             "The independent maintained-T control requires external heating and is not an isolated-wall energy balance.",
                             "No conduction, radiation, phase transition, motion, prescribed input power, state chemistry, gap evolution, spectrum or host transaction is modeled."]}
        _write_json(stage / "report.json", report)
        stage.rename(output)
    except Exception as exc:
        _write_json(stage / "failure.json", {"status": "incomplete_reference_not_admitted", "error_type": type(exc).__name__,
                                             "error": str(exc), "completed_distinct_trajectories": len(records),
                                             "no_fallback_or_retry_performed": True})
        raise ValueError("reference failed closed; incomplete evidence retained at " + str(stage) + ": " + str(exc)) from exc
    print(json.dumps({"report": str(output / "report.json"), "report_sha256": sha(output / "report.json"),
                      "case_count": len(config["cases"]), "numerical_trajectory_count": len(records),
                      "temporal_refinement_rows": len(refinement), "status": "manufactured_reference_only"}, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/thermal_feedback_manufactured.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    try:
        run(args.config, args.output_dir, make_plots=args.plot)
    except (ValueError, TypeError, OverflowError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
