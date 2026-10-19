#!/usr/bin/env python3
"""Save manufactured positive-barrier single-patch forced-feedback references.

Only the separately qualified strictly positive base-barrier model is used. Every
physical coefficient and every finite schedule segment is supplied explicitly.
Raw JSON decimal schedule endpoints are checked before binary64 projection.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CONFIG_SCHEMA = "mhdlab-thermal-forced-feedback-manufactured-config-v1"
REPORT_SCHEMA = "mhdlab-thermal-forced-feedback-reference-run-v1"
MAX_STEPS = 16384
MAX_TRAJECTORIES = 24
_ID = re.compile(r"[a-z][a-z0-9_]{0,63}\Z")
_METRICS = ("coverage_m2", "temperature_k", "released_m2", "input_energy_j_m2",
            "binding_energy_j_m2", "internal_energy_j_m2", "translational_energy_j_m2",
            "total_wall_energy_debit_j_m2", "gas_energy_j_m2", "wall_thermal_change_j_m2")
_POOL_TEXT = {"name", "state_id", "kinetic_provenance", "energy_provenance"}
_POOL_POS = {"attempt_frequency_s", "reference_coverage_m2", "max_coverage_m2"}
_POOL_NONNEG = {"activation_energy_ev", "desorption_energy_ev", "internal_energy_ev", "coverage_coefficient"}
_PROBLEM_TEXT = {"parameter_status", "parameter_provenance", "energy_reference"}
_PROBLEM_POS = {"initial_temperature_k", "areal_heat_capacity_j_m2_k", "patch_area_m2", "temperature_min_k", "temperature_max_k"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reject_constant(token):
    raise ValueError("nonfinite JSON constant: " + token)


def _decimal(token):
    try:
        value = Decimal(token)
    except InvalidOperation as exc:
        raise ValueError("unsupported JSON numeric lexeme") from exc
    if not value.is_finite():
        raise ValueError("nonfinite JSON numeric lexeme")
    projected = float(value)
    if not math.isfinite(projected):
        raise ValueError("JSON numeric lexeme overflows finite binary64 before exact domain arithmetic")
    if value != 0 and projected == 0:
        raise ValueError("nonzero JSON numeric lexeme underflows binary64 before exact domain arithmetic")
    return value


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate configuration key: " + key)
        result[key] = value
    return result


def _number(value, name, *, positive=False):
    if type(value) not in (int, Decimal) or value < 0 or (positive and value == 0):
        raise ValueError(name + " must be an explicitly supplied " + ("positive" if positive else "nonnegative") + " finite number")
    return value


def _text(value, name):
    if type(value) is not str or not value.strip():
        raise ValueError(name + " must be explicit nonempty text")


def _project(value):
    if type(value) is Decimal:
        result = float(value)
        if not math.isfinite(result):
            raise ValueError("JSON number overflows finite binary64")
        if value != 0 and result == 0:
            raise ValueError("nonzero JSON number underflows binary64")
        return result
    if type(value) is int:
        try:
            result = float(value)
        except OverflowError as exc:
            raise ValueError("JSON integer overflows finite binary64") from exc
        if not math.isfinite(result):
            raise ValueError("JSON integer overflows finite binary64")
        return value
    if type(value) is list:
        return [_project(item) for item in value]
    if type(value) is dict:
        return {key: _project(item) for key, item in value.items()}
    return value


def load_config(path: Path) -> dict:
    """Validate original decimal support, then project supplied inputs to SI floats."""
    raw = json.loads(path.read_text(), parse_float=_decimal, parse_constant=_reject_constant, object_pairs_hook=_object)
    required = {"schema", "parameter_status", "provenance", "scope", "schedule_time_convention", "refinement_case_id", "cases"}
    if type(raw) is not dict or set(raw) != required or raw["schema"] != CONFIG_SCHEMA:
        raise ValueError("unsupported manufactured forced-feedback configuration")
    if raw["parameter_status"] != "manufactured":
        raise ValueError("only explicit manufactured parameters are accepted")
    for name in ("provenance", "scope", "schedule_time_convention"):
        _text(raw[name], name)
    if type(raw["cases"]) is not list or not 1 <= len(raw["cases"]) <= 6:
        raise ValueError("one to six explicit cases are supported")
    ids, run_count, schedule_checks = set(), 0, []
    for case in raw["cases"]:
        if type(case) is not dict or set(case) != {"case_id", "purpose", "duration_s", "segments", "step_grids", "problem"}:
            raise ValueError("each case must supply all declared inputs")
        identifier = case["case_id"]
        if type(identifier) is not str or not _ID.fullmatch(identifier) or identifier in ids:
            raise ValueError("case identifiers must be unique bounded lowercase identifiers")
        ids.add(identifier)
        _text(case["purpose"], "purpose")
        duration = _number(case["duration_s"], "duration_s", positive=True)
        segments = case["segments"]
        if type(segments) is not list or not 1 <= len(segments) <= 64:
            raise ValueError("each case needs one to 64 explicit finite segments")
        previous, checks = 0, []
        for index, segment in enumerate(segments):
            if type(segment) is not dict or set(segment) != {"start_s", "stop_s", "power_w_m2", "power_provenance"}:
                raise ValueError("each segment must explicitly declare endpoints, power and provenance")
            start = _number(segment["start_s"], "segment start_s")
            stop = _number(segment["stop_s"], "segment stop_s", positive=True)
            _number(segment["power_w_m2"], "power_w_m2")
            _text(segment["power_provenance"], "power_provenance")
            if start != previous or stop <= start:
                raise ValueError("raw decimal schedule must have positive segments and exact complete adjacency; no gap or overlap")
            a, b = float(start), float(stop)
            if not math.isfinite(a) or not math.isfinite(b) or b <= a:
                raise ValueError("positive raw segment duration collapses or overflows its binary64 endpoints")
            checks.append({"segment_index": index, "raw_start_decimal": str(start), "raw_stop_decimal": str(stop),
                           "projected_start_hex": a.hex(), "projected_stop_hex": b.hex()})
            previous = stop
        if previous != duration:
            raise ValueError("raw decimal schedule must exactly cover duration_s; no implicit final hold")
        grids = case["step_grids"]
        if type(grids) is not list or not grids:
            raise ValueError("step_grids must explicitly supply per-segment positive counts")
        prior_total, seen = 0, set()
        for grid in grids:
            if (type(grid) is not list or len(grid) != len(segments)
                    or any(type(n) is not int or not 1 <= n <= MAX_STEPS for n in grid)
                    or sum(grid) > MAX_STEPS or tuple(grid) in seen or sum(grid) <= prior_total):
                raise ValueError("step grids must be distinct, bounded and increasing in total work with one integer per segment")
            prior_total = sum(grid)
            seen.add(tuple(grid))
        problem = case["problem"]
        if type(problem) is not dict or set(problem) != _PROBLEM_TEXT | _PROBLEM_POS | {"reservoir", "initial_coverage_m2"}:
            raise ValueError("each problem must supply every physical coefficient and provenance field")
        for name in _PROBLEM_TEXT:
            _text(problem[name], name)
        if problem["parameter_status"] != "manufactured":
            raise ValueError("each problem must declare manufactured status")
        for name in _PROBLEM_POS:
            _number(problem[name], name, positive=True)
        _number(problem["initial_coverage_m2"], "initial_coverage_m2")
        pool = problem["reservoir"]
        if type(pool) is not dict or set(pool) != _POOL_TEXT | _POOL_POS | _POOL_NONNEG:
            raise ValueError("each reservoir must supply every physical coefficient and provenance field")
        for name in _POOL_TEXT:
            _text(pool[name], name)
        for name in _POOL_POS:
            _number(pool[name], name, positive=True)
        for name in _POOL_NONNEG:
            _number(pool[name], name)
        if pool["activation_energy_ev"] <= 0:
            raise ValueError("this forced-feedback reference requires strictly positive activation energy")
        if Fraction(pool["coverage_coefficient"])*Fraction(pool["max_coverage_m2"]) > Fraction(pool["reference_coverage_m2"]):
            raise ValueError("raw declared effective barrier becomes negative within reservoir support")
        if problem["initial_coverage_m2"] > pool["max_coverage_m2"]:
            raise ValueError("initial coverage exceeds declared reservoir support")
        if not problem["temperature_min_k"] <= problem["initial_temperature_k"] <= problem["temperature_max_k"] or problem["temperature_min_k"] >= problem["temperature_max_k"]:
            raise ValueError("initial and supported temperatures are inconsistent")
        run_count += len(grids)
        schedule_checks.append({"case_id": identifier, "segments": checks, "raw_adjacency_and_support_checked": True})
    if run_count > MAX_TRAJECTORIES or raw["refinement_case_id"] not in ids:
        raise ValueError("invalid refinement case or configured work limit")
    result = _project(raw)
    result["raw_schedule_validation"] = schedule_checks
    return result


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _flat_sample(sample):
    record = sample.to_dict()
    extensive = record.pop("extensive")
    if type(extensive) is not dict or set(record) & set(extensive):
        raise ValueError("unsupported sample/extensive mapping")
    flat = {**record, **extensive}
    if any(type(value) not in (int, float) or not math.isfinite(value) for value in flat.values()):
        raise ValueError("unsupported numerical sample representation")
    return flat


def _units(name):
    if name == "time_s":
        return "s"
    if name == "temperature_k":
        return "K"
    if name.endswith("_w_m2"):
        return "W m^-2"
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
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format(value, ".17g") if type(value) is float else value for key, value in row.items()})
    return {key: _units(key) for key in rows[0] if key != "case_id"}


def _refinement(records, case, path):
    family = [record for record in records if record["case_id"] == case["case_id"]]
    result, previous = [], {}
    for coarse, fine in zip(family[:-1], family[1:]):
        for metric in _METRICS:
            a, b = coarse["final_sample"][metric], fine["final_sample"][metric]
            change = b - a
            final = fine["final_sample"]
            if metric in ("coverage_m2", "released_m2"):
                budget, convention = final["released_m2"], "fine emitted areal parent inventory"
            elif metric == "temperature_k":
                budget = max(abs(final["input_energy_j_m2"]), abs(final["total_wall_energy_debit_j_m2"])) / case["problem"]["areal_heat_capacity_j_m2_k"]
                convention = "max gross input or wall debit divided by supplied areal heat capacity"
            elif metric == "wall_thermal_change_j_m2":
                budget = max(abs(final["input_energy_j_m2"]), abs(final["total_wall_energy_debit_j_m2"]))
                convention = "max gross input or wall debit; signed net thermal change retained separately"
            else:
                budget, convention = abs(b), "absolute fine channel budget"
            prior = previous.get(metric)
            doubled = all(y == 2*x for x, y in zip(coarse["segment_steps"], fine["segment_steps"]))
            order = (math.log(prior)-math.log(abs(change)))/math.log(2) if prior and change and doubled else None
            result.append({"case_id": case["case_id"], "coarse_total_steps": coarse["total_steps"],
                           "fine_total_steps": fine["total_steps"], "quantity": metric, "units": _units(metric),
                           "signed_final_change": change, "absolute_final_change": abs(change),
                           "normalizing_budget": budget, "normalizing_budget_definition": convention,
                           "change_over_budget": None if budget == 0 else abs(change)/budget,
                           "zero_budget_exact_zero_change": budget == 0 and change == 0,
                           "observed_self_order": order,
                           "interpretation": "finite-grid changes only; no independent error bound or accuracy certificate"})
            previous[metric] = abs(change)
    with path.open("w", newline="") as stream:
        fields = tuple(result[0]) if result else ("case_id", "quantity")
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in result:
            writer.writerow({key: "" if value is None else format(value, ".17g") if type(value) is float else value for key, value in row.items()})
    return result


def _plots(config, saved, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for ax, case in zip(axes.flat, config["cases"]):
        entry = max((entry for entry in saved if entry["case_id"] == case["case_id"]), key=lambda item: item["total_steps"])
        rows = entry["rows"]
        ax.plot([row["time_s"]*1e9 for row in rows], [row["temperature_k"] for row in rows], color="firebrick")
        ax.set(xlabel="Time [ns]", ylabel="Wall temperature [K]", title=case["case_id"])
        ax.grid(alpha=.2)
        twin = ax.twinx()
        n0 = case["problem"]["initial_coverage_m2"]
        if n0:
            twin.plot([row["time_s"]*1e9 for row in rows], [row["released_m2"]/n0 for row in rows], "--", color="navy")
            twin.set_ylabel("Emitted / supplied initial inventory")
        else:
            twin.text(.4, .6, "Exactly empty inventory", transform=twin.transAxes)
            twin.set_yticks([])
    figure.suptitle("Manufactured positive-barrier forced-feedback reference")
    figure.savefig(output/"trajectories.png", dpi=160)
    figure.savefig(output/"trajectories.pdf")
    plt.close(figure)
    family = [entry for entry in saved if entry["case_id"] == config["refinement_case_id"]]
    rows = max(family, key=lambda item: item["total_steps"])["rows"]
    figure, ax = plt.subplots(figsize=(10, 6), layout="constrained")
    channels = (("input_energy_j_m2", "Gross prescribed input"), ("binding_energy_j_m2", "Removal / adsorbate-reference increase"),
                ("internal_energy_j_m2", "Gas internal"), ("translational_energy_j_m2", "Expected gas translation"),
                ("total_wall_energy_debit_j_m2", "Gross wall debit"), ("wall_thermal_change_j_m2", "Signed thermal change"))
    for key, label in channels:
        ax.plot([row["time_s"]*1e9 for row in rows], [row[key] for row in rows], label=label)
    ax.set(xlabel="Time [ns]", ylabel="Energy per area [J m$^{-2}$]", title="Manufactured segmented input: gross transfers and signed net change")
    ax.legend(fontsize=8)
    ax.grid(alpha=.2)
    figure.savefig(output/"energy_channels.png", dpi=160)
    figure.savefig(output/"energy_channels.pdf")
    plt.close(figure)


def run(config_path: Path, output: Path, *, make_plots: bool):
    config = load_config(config_path)
    if output.exists():
        raise ValueError("output must be a fresh directory; reference evidence is never overwritten")
    sys.path.insert(0, str(ROOT))
    from mhdlab.thermal_feedback_reference import FeedbackReservoir, FeedbackProblem
    from mhdlab.thermal_power_reference import PowerSegment, PowerSchedule
    from mhdlab.thermal_forced_feedback_reference import evolve_forced_feedback
    plans = []
    for case in config["cases"]:
        supplied = dict(case["problem"])
        pool = FeedbackReservoir(**supplied.pop("reservoir"))
        problem = FeedbackProblem(reservoir=pool, **supplied)
        for grid_index, grid in enumerate(case["step_grids"]):
            segments = tuple(PowerSegment(**segment, steps=steps) for segment, steps in zip(case["segments"], grid))
            schedule = PowerSchedule(duration_s=case["duration_s"], segments=segments, parameter_status="manufactured", parameter_provenance=config["provenance"])
            plans.append((case, grid_index, problem, schedule, grid))
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=output.name+".staging-", dir=output.parent))
    records, saved = [], []
    try:
        shutil.copyfile(config_path, stage/"input_config.json")
        for case, grid_index, problem, schedule, grid in plans:
            history = evolve_forced_feedback(problem, schedule)
            rows = [{"sample_index": index, "case_id": case["case_id"], "grid_index": grid_index, "total_steps": sum(grid), **_flat_sample(sample)} for index, sample in enumerate(history.samples)]
            stem = case["case_id"]+"_grid"+str(grid_index)+"_steps"+format(sum(grid), "04d")
            csv_path, history_path = stage/(stem+".csv"), stage/(stem+".json")
            units = _write_csv(csv_path, rows)
            _write_json(history_path, history.to_dict())
            records.append({"case_id": case["case_id"], "purpose": case["purpose"], "grid_index": grid_index,
                            "segment_steps": grid, "total_steps": sum(grid), "duration_s": case["duration_s"],
                            "sample_count": len(rows), "csv_file": csv_path.name, "csv_sha256": sha(csv_path),
                            "history_file": history_path.name, "history_sha256": sha(history_path),
                            "input_parameters": case["problem"], "input_segments": case["segments"], "units": units,
                            "final_sample": rows[-1],
                            "max_absolute_residuals": {key: max(abs(row[key]) for row in rows) for key in rows[0] if "residual" in key}})
            saved.append({"case_id": case["case_id"], "total_steps": sum(grid), "rows": rows})
        refinement_case = next(case for case in config["cases"] if case["case_id"] == config["refinement_case_id"])
        refinement = _refinement(records, refinement_case, stage/"temporal_refinement.csv")
        if make_plots:
            _plots(config, saved, stage)
        code = {"tools/run_thermal_forced_feedback_reference.py": sha(Path(__file__)),
                "mhdlab/thermal_forced_feedback_reference.py": sha(ROOT/"mhdlab/thermal_forced_feedback_reference.py"),
                "mhdlab/thermal_power_reference.py": sha(ROOT/"mhdlab/thermal_power_reference.py"),
                "mhdlab/thermal_feedback_reference.py": sha(ROOT/"mhdlab/thermal_feedback_reference.py")}
        artifacts = {path.name: sha(path) for path in stage.iterdir() if path.is_file()}
        report = {"schema": REPORT_SCHEMA, "status": "manufactured_reference_only",
                  "config_original_path": str(config_path.resolve()), "config_sha256": sha(config_path), "code_sha256": code,
                  "parameter_status": "manufactured", "scope": config["scope"],
                  "case_count": len(config["cases"]), "numerical_trajectory_count": len(records),
                  "sample_count": sum(record["sample_count"] for record in records),
                  "trajectories": records, "temporal_refinement": refinement,
                  "raw_schedule_validation": config["raw_schedule_validation"],
                  "schedule_time_convention": config["schedule_time_convention"],
                  "sample_time_convention": "segment_index and segment_step identify the exact logical grid; time_s is a rounded binary64 stamp; one shared state is stored at each boundary, attributed to the completed segment",
                  "numeric_input_convention": "JSON decimal lexemes retained through raw support checks, then explicitly projected to finite binary64 without nonzero-to-zero conversion; integer step counts remain integers",
                  "artifact_path_base": "directory containing report.json", "artifact_sha256": artifacts,
                  "qualifications": {"explicit_finite_prescribed_power": True, "strictly_positive_base_activation_energy": True, "temperature_dependent_coverage_barrier": True,
                                     "gross_input_and_debit_separate": True, "signed_net_thermal_change_preserved": True,
                                     "expected_flux_translation_2kBT": True, "independent_reference_error_computed_here": False,
                                     "old_source_power_or_feedback_evolution_called": False, "physical_parameters_or_water_inferred": False,
                                     "actual_FLASH_or_host_coupling": False, "floor_clamp_or_unsupported_fallback": False},
                  "limits": ["Manufactured single stationary patch with explicit constant areal heat capacity and finite nonnegative power schedule; no physical material validation.",
                             "Strictly positive base activation energy is required; the supported coverage-dependent barrier is nonnegative throughout the declared inventory domain and Nref stays fixed.",
                             "Removal, gas internal, continuum expected translation, gross external input and signed thermal change remain separate; signed residuals are never repaired.",
                             "Finite-grid changes and observed orders are descriptive, not independent errors or certified bounds.",
                             "No conduction, radiation, phase change, motion, source history, gap, chemistry, spectrum, actual MHD or host transaction is modeled."]}
        _write_json(stage/"report.json", report)
        stage.rename(output)
    except Exception as exc:
        _write_json(stage/"failure.json", {"status": "incomplete_reference_not_admitted", "error_type": type(exc).__name__, "error": str(exc), "completed_distinct_trajectories": len(records), "no_fallback_or_retry_performed": True})
        raise ValueError("reference failed closed; incomplete evidence retained at "+str(stage)+": "+str(exc)) from exc
    print(json.dumps({"report": str(output/"report.json"), "report_sha256": sha(output/"report.json"), "case_count": len(config["cases"]), "numerical_trajectory_count": len(records), "sample_count": report["sample_count"], "temporal_refinement_rows": len(refinement), "status": report["status"]}, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"configs/thermal_forced_feedback_manufactured.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    try:
        run(args.config, args.output_dir, make_plots=args.plot)
    except (ValueError, TypeError, OverflowError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
