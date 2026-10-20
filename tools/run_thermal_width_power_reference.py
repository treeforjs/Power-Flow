#!/usr/bin/env python3
"""Save three new manufactured scalar histories and an exact-area width comparison.

Identical width patches reuse immutable histories. The adapter performs no
integration; local temperatures and source mixtures remain separately labeled.
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
CONFIG_SCHEMA = "mhdlab-thermal-width-power-manufactured-config-v1"
REPORT_SCHEMA = "mhdlab-thermal-width-power-reference-run-v1"
_ID = re.compile(r"[a-z][a-z0-9_]{0,63}\Z")
_INT = re.compile(r"-?(?:0|[1-9][0-9]*)\Z")
_POOL_TEXT = {"name", "state_id", "kinetic_provenance", "energy_provenance"}
_POOL_POS = {"attempt_frequency_s", "reference_coverage_m2", "max_coverage_m2"}
_POOL_NONNEG = {"activation_energy_ev", "desorption_energy_ev", "internal_energy_ev", "coverage_coefficient"}
_PROBLEM_TEXT = {"parameter_status", "parameter_provenance", "energy_reference"}
_PROBLEM_POS = {"initial_temperature_k", "areal_heat_capacity_j_m2_k", "patch_area_m2", "temperature_min_k", "temperature_max_k"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _constant(token):
    raise ValueError("nonfinite JSON constant: " + token)


def _decimal(token):
    try:
        value = Decimal(token)
    except InvalidOperation as exc:
        raise ValueError("unsupported numeric lexeme") from exc
    projected = float(value)
    if not value.is_finite() or not math.isfinite(projected):
        raise ValueError("numeric lexeme exceeds finite binary64")
    if value != 0 and projected == 0:
        raise ValueError("nonzero numeric lexeme underflows binary64")
    return value


def _object(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate configuration key: " + key)
        value[key] = item
    return value


def _text(value, name):
    if type(value) is not str or not value.strip():
        raise ValueError(name + " must be nonempty supplied text")


def _number(value, name, *, positive=False):
    if type(value) not in (int, Decimal) or value < 0 or (positive and value == 0):
        raise ValueError(name + " must be a supplied " + ("positive" if positive else "nonnegative") + " finite number")
    return value


def _project(value):
    if type(value) is Decimal:
        return float(value)  # finite/nonzero qualification already happened in the raw parser
    if type(value) is int:
        try:
            if not math.isfinite(float(value)):
                raise ValueError("integer exceeds finite binary64")
        except OverflowError as exc:
            raise ValueError("integer exceeds finite binary64") from exc
        return value
    if type(value) is list:
        return [_project(x) for x in value]
    if type(value) is dict:
        return {k: _project(x) for k, x in value.items()}
    return value


def _rational(record, name):
    if type(record) is not dict or set(record) != {"numerator", "denominator", "value"}:
        raise ValueError(name + " needs an explicit canonical rational and numeric value")
    for key in ("numerator", "denominator"):
        token = record[key]
        if type(token) is not str or len(token) > 128 or not _INT.fullmatch(token) or str(int(token)) != token:
            raise ValueError(name + " rational integers must be bounded canonical strings")
    numerator, denominator = int(record["numerator"]), int(record["denominator"])
    if denominator <= 0:
        raise ValueError(name + " rational denominator must be positive")
    exact = Fraction(numerator, denominator)
    if exact.numerator != numerator or exact.denominator != denominator:
        raise ValueError(name + " rational must be reduced")
    supplied = record["value"]
    if type(supplied) not in (int, Decimal) or Fraction(supplied) != exact:
        raise ValueError(name + " raw numeric value differs from its exact rational")
    try:
        projected = float(exact)
    except OverflowError as exc:
        raise ValueError(name + " geometry exceeds binary64") from exc
    if not math.isfinite(projected) or Fraction.from_float(projected) != exact:
        raise ValueError(name + " geometry must be exactly binary64-representable")
    return exact


def _fraction_record(value):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator), "value": float(value)}


def load_config(path):
    raw = json.loads(Path(path).read_text(), parse_float=_decimal, parse_constant=_constant, object_pairs_hook=_object)
    keys = {"schema", "parameter_status", "provenance", "scope", "geometry", "duration_s", "steps", "common_scalar_problem", "scalar_cases", "hypotheses", "time_convention", "geometry_convention"}
    if type(raw) is not dict or set(raw) != keys or raw["schema"] != CONFIG_SCHEMA:
        raise ValueError("unsupported manufactured width-power configuration")
    if raw["parameter_status"] != "manufactured":
        raise ValueError("only explicitly manufactured inputs are supported")
    for key in ("provenance", "scope", "time_convention", "geometry_convention"):
        _text(raw[key], key)
    duration = _number(raw["duration_s"], "duration_s", positive=True)
    if float(duration) <= 0:
        raise ValueError("finite support collapses during projection")
    if type(raw["steps"]) is not int or raw["steps"] != 64:
        raise ValueError("this bounded reference requires exactly64 explicitly supplied steps per scalar history")
    geometry = raw["geometry"]
    if type(geometry) is not dict or set(geometry) != {"width_edges_m", "common_length_m"}:
        raise ValueError("geometry must explicitly supply width edges and common length")
    if type(geometry["width_edges_m"]) is not list or len(geometry["width_edges_m"]) != 4:
        raise ValueError("exactly three manufactured width patches are supported")
    edges = [_rational(x, "width edge") for x in geometry["width_edges_m"]]
    length = _rational(geometry["common_length_m"], "common length")
    if length <= 0 or any(b <= a for a, b in zip(edges[:-1], edges[1:])):
        raise ValueError("geometry must have positive length and strictly increasing edges")
    areas = [(b-a)*length for a, b in zip(edges[:-1], edges[1:])]
    if any(not math.isfinite(float(x)) or float(x) <= 0 for x in areas):
        raise ValueError("positive patch area is not representable")
    problem = raw["common_scalar_problem"]
    if type(problem) is not dict or set(problem) != _PROBLEM_TEXT | _PROBLEM_POS | {"reservoir", "initial_coverage_m2"}:
        raise ValueError("complete scalar problem parameters are required")
    for key in _PROBLEM_TEXT:
        _text(problem[key], key)
    if problem["parameter_status"] != "manufactured":
        raise ValueError("scalar problem must explicitly be manufactured")
    for key in _PROBLEM_POS:
        _number(problem[key], key, positive=True)
    _number(problem["initial_coverage_m2"], "initial_coverage_m2", positive=True)
    if problem["patch_area_m2"] != 1:
        raise ValueError("scalar reference area must explicitly equal one square metre")
    pool = problem["reservoir"]
    if type(pool) is not dict or set(pool) != _POOL_TEXT | _POOL_POS | _POOL_NONNEG:
        raise ValueError("complete reservoir parameters are required")
    for key in _POOL_TEXT:
        _text(pool[key], key)
    for key in _POOL_POS:
        _number(pool[key], key, positive=True)
    for key in _POOL_NONNEG:
        _number(pool[key], key)
    if pool["activation_energy_ev"] <= 0:
        raise ValueError("strictly positive base activation energy is required")
    if Fraction(pool["coverage_coefficient"])*Fraction(pool["max_coverage_m2"]) > Fraction(pool["reference_coverage_m2"]):
        raise ValueError("effective barrier is negative inside the declared support")
    if problem["initial_coverage_m2"] > pool["max_coverage_m2"]:
        raise ValueError("initial inventory exceeds declared support")
    if not problem["temperature_min_k"] <= problem["initial_temperature_k"] <= problem["temperature_max_k"] or problem["temperature_min_k"] >= problem["temperature_max_k"]:
        raise ValueError("declared temperature support is inconsistent")
    cases = raw["scalar_cases"]
    if type(cases) is not list or len(cases) != 3:
        raise ValueError("exactly three distinct scalar cases are required")
    ids, powers = set(), {}
    for case in cases:
        if type(case) is not dict or set(case) != {"case_id", "power_w_m2", "power_provenance"}:
            raise ValueError("each scalar forcing must be explicitly declared")
        identifier = case["case_id"]
        if type(identifier) is not str or not _ID.fullmatch(identifier) or identifier in ids:
            raise ValueError("scalar case ids must be unique bounded identifiers")
        ids.add(identifier)
        _text(case["power_provenance"], "power_provenance")
        powers[identifier] = Fraction(_number(case["power_w_m2"], "power_w_m2", positive=True))
    if ids != {"uniform_power", "edge_power", "middle_power"} or len(set(powers.values())) != 3:
        raise ValueError("three explicit distinct uniform/edge/middle histories are required")
    expected = {"uniform": ["uniform_power"]*3, "nonuniform": ["edge_power", "middle_power", "edge_power"]}
    if type(raw["hypotheses"]) is not dict or raw["hypotheses"] != expected:
        raise ValueError("the declared bounded reuse map must be explicit and complete")
    totals = {name: sum((area*powers[key] for area, key in zip(areas, assigned)), Fraction()) for name, assigned in expected.items()}
    if totals["uniform"] != totals["nonuniform"]:
        raise ValueError("total prescribed input powers must match exactly without normalization")
    result = _project(raw)
    # Recheck equality in the actual represented model coordinates and powers.
    represented = {case["case_id"]: Fraction.from_float(float(case["power_w_m2"])) for case in result["scalar_cases"]}
    projected_totals = {name: sum((area*represented[key] for area, key in zip(areas, assigned)), Fraction()) for name, assigned in expected.items()}
    if projected_totals["uniform"] != projected_totals["nonuniform"]:
        raise ValueError("projected total powers do not match exactly")
    result["exact_input_matching"] = {"raw_decimal_power_equality_checked": True, "represented_power_equality_checked": True,
        "patch_areas_m2": [_fraction_record(x) for x in areas], "total_area_m2": _fraction_record(sum(areas)),
        "total_input_power_w": _fraction_record(projected_totals["uniform"]),
        "total_input_energy_j": _fraction_record(projected_totals["uniform"]*Fraction.from_float(float(result["duration_s"])))}
    return result


def _write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def _flat_sample(sample):
    value = sample.to_dict()
    extensive = value.pop("extensive")
    if type(extensive) is not dict or set(value) & set(extensive):
        raise ValueError("invalid scalar sample contract")
    flat = {**value, **extensive}
    if any(type(x) not in (int, float) or not math.isfinite(x) for x in flat.values()):
        raise ValueError("unsupported scalar sample representation")
    return flat


def _write_csv(path, rows):
    if not rows:
        raise ValueError("required output table is empty")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format(value, ".17g") if type(value) is float else json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) if type(value) in (dict, list) else value for key, value in row.items()})


def _record_columns(prefix, record):
    # Preserve every signed projection diagnostic, not only its rounded value.
    return {prefix: record["value"], prefix+"__exact_numerator": record["exact"]["numerator"],
            prefix+"__exact_denominator": record["exact"]["denominator"], prefix+"__projection_json": record["projection"]}


def _assembly_csvs(document, output):
    patch_rows, total_rows, difference_rows = [], [], []
    for name, hypothesis in document["hypotheses"].items():
        for patch in hypothesis["patches"]:
            for index, sample in enumerate(patch["samples"]):
                row = {"hypothesis": name, "patch_index": patch["patch_index"], "sample_index": index,
                       "width_left_m": patch["width_left_m"], "width_right_m": patch["width_right_m"],
                       "area_m2": patch["area_m2"], "area_exact_json": patch["area_exact_m2"],
                       "state_id": patch["state_id"], "scalar_history_id": patch["scalar_history_id"],
                       "history_sha256": patch["history_sha256"], "segment_index": sample["segment_index"],
                       "segment_step": sample["segment_step"], "time_s": sample["time_s"], "logical_time_json": sample["logical_time_s"]}
                row.update({"areal__"+key: value for key, value in sample["areal"].items()})
                for key, record in sample["extensive"].items():
                    row.update(_record_columns("extensive__"+key, record))
                patch_rows.append(row)
        for index, sample in enumerate(hypothesis["totals"]):
            row = {"hypothesis": name, "sample_index": index, "segment_index": sample["segment_index"], "segment_step": sample["segment_step"],
                   "time_s": sample["time_s"], "logical_time_json": sample["logical_time_s"], "state_id": sample["state_id"]}
            for section in ("extensive", "recomputed_residuals"):
                for key, record in sample[section].items():
                    row.update(_record_columns(section+"__"+key, record))
            total_rows.append(row)
    for index, sample in enumerate(document["comparison"]["samples"]):
        row = {"sample_index": index, "time_s": sample["time_s"], "logical_time_json": sample["logical_time_s"]}
        for key, record in sample["difference_nonuniform_minus_uniform"].items():
            row.update(_record_columns("difference__"+key, record))
        for key in ("exact_prescribed_input_energy_j", "integrated_input_difference_j", "prescribed_input_difference_j"):
            row.update(_record_columns(key, sample[key]))
        difference_rows.append(row)
    for name, rows in (("patch_samples.csv", patch_rows), ("aggregate_samples.csv", total_rows), ("width_comparison.csv", difference_rows)):
        _write_csv(output/name, rows)
    return {"expanded_patch_samples": len(patch_rows), "aggregate_samples": len(total_rows), "comparison_samples": len(difference_rows)}


def _plots(config, scalar_rows, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    edges = [float(x["value"]) for x in config["geometry"]["width_edges_m"]]
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.8), layout="constrained")
    fields = (("temperature_k", "Final patch temperature [K]"), ("released_m2", "Final emitted parents [m$^{-2}$]"),
              ("translational_energy_j_m2", "Final expected translation [J m$^{-2}$]"))
    for ax, (field, label) in zip(axes, fields):
        for name, assignments in config["hypotheses"].items():
            values = [scalar_rows[key][-1][field] for key in assignments]
            ax.stairs(values, [x*1e3 for x in edges], label=name, linewidth=2, baseline=None)
        ax.set(xlabel="Manufactured width coordinate [mm]", ylabel=label)
        ax.grid(alpha=.2)
    axes[-1].legend()
    figure.suptitle("Manufactured equal-total-input width comparison; independent patches")
    figure.savefig(output/"width_profiles.png", dpi=160)
    figure.savefig(output/"width_profiles.pdf")
    plt.close(figure)
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
    for key, rows in scalar_rows.items():
        axes[0].plot([r["time_s"]*1e9 for r in rows], [r["temperature_k"] for r in rows], label=key)
        axes[1].plot([r["time_s"]*1e9 for r in rows], [r["released_m2"]/config["common_scalar_problem"]["initial_coverage_m2"] for r in rows], label=key)
    axes[0].set(xlabel="Time [ns]", ylabel="Temperature [K]")
    axes[1].set(xlabel="Time [ns]", ylabel="Emitted / supplied initial inventory")
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    figure.suptitle("Three distinct new manufactured scalar histories; reused without re-evolution")
    figure.savefig(output/"scalar_histories.png", dpi=160)
    figure.savefig(output/"scalar_histories.pdf")
    plt.close(figure)


def run(config_path, output, *, make_plots):
    config = load_config(config_path)
    if output.exists():
        raise ValueError("output must be a fresh directory")
    sys.path.insert(0, str(ROOT))
    from mhdlab.thermal_feedback_reference import FeedbackProblem, FeedbackReservoir
    from mhdlab.thermal_power_reference import PowerSegment, PowerSchedule
    from mhdlab.thermal_forced_feedback_reference import evolve_forced_feedback
    from mhdlab.thermal_width_power_reference import WidthPowerGeometry, assemble_width_power_reference, scalar_history_sha256
    supplied = dict(config["common_scalar_problem"])
    problem = FeedbackProblem(reservoir=FeedbackReservoir(**supplied.pop("reservoir")), **supplied)
    geometry = WidthPowerGeometry(width_edges_m=tuple(float(x["value"]) for x in config["geometry"]["width_edges_m"]),
        common_length_m=float(config["geometry"]["common_length_m"]["value"]), parameter_status="manufactured", parameter_provenance=config["provenance"])
    plans = []
    for case in config["scalar_cases"]:
        schedule = PowerSchedule(duration_s=config["duration_s"], segments=(PowerSegment(start_s=0.0, stop_s=config["duration_s"],
            power_w_m2=case["power_w_m2"], steps=config["steps"], power_provenance=case["power_provenance"]),),
            parameter_status="manufactured", parameter_provenance=config["provenance"])
        plans.append((case, schedule))
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=output.name+".staging-", dir=output.parent))
    histories, records, scalar_rows = {}, [], {}
    try:
        shutil.copyfile(config_path, stage/"input_config.json")
        # The only evolution loop: one call per distinct configured scalar input.
        for case, schedule in plans:
            key = case["case_id"]
            history = evolve_forced_feedback(problem, schedule)
            histories[key] = history
            rows = [{"sample_index": i, "case_id": key, **_flat_sample(sample)} for i, sample in enumerate(history.samples)]
            scalar_rows[key] = rows
            history_path, csv_path = stage/(key+".json"), stage/(key+".csv")
            _write_json(history_path, history.to_dict())
            _write_csv(csv_path, rows)
            records.append({"case_id": key, "power_w_m2": case["power_w_m2"], "steps": config["steps"], "sample_count": len(rows),
                "history_file": history_path.name, "history_file_sha256": sha(history_path), "canonical_history_sha256": scalar_history_sha256(history),
                "csv_file": csv_path.name, "csv_sha256": sha(csv_path), "final_sample": rows[-1],
                "max_absolute_residuals": {k: max(abs(r[k]) for r in rows) for k in rows[0] if "residual" in k}})
        digest = {x["case_id"]: x["canonical_history_sha256"] for x in records}
        assembly = assemble_width_power_reference(histories, geometry=geometry, hypotheses=config["hypotheses"], history_sha256=digest,
            history_provenance={case["case_id"]: case["power_provenance"] for case in config["scalar_cases"]}, parameter_status="manufactured", parameter_provenance=config["provenance"])
        document = assembly.to_dict()
        # Export the adapter's cached saved-state RHS, with no extra evaluation.
        for record in records:
            key = record["case_id"]
            cache = document["scalar_histories"][key]
            if len(cache["rhs"]) != len(scalar_rows[key]) or len(cache["release_flux_m2_s"]) != len(scalar_rows[key]):
                raise ValueError("cached scalar RHS support differs from the saved history")
            for row, rhs, flux in zip(scalar_rows[key], cache["rhs"], cache["release_flux_m2_s"]):
                row["state_id"] = config["common_scalar_problem"]["reservoir"]["state_id"]
                row["release_flux_m2_s"] = flux
                row.update({"instantaneous__"+name: value for name, value in rhs.items()})
            _write_csv(stage/record["csv_file"], scalar_rows[key])
            record["csv_sha256"] = sha(stage/record["csv_file"])
        _write_json(stage/"width_assembly.json", document)
        counts = _assembly_csvs(document, stage)
        _write_json(stage/"reuse_map.json", {"hypotheses": config["hypotheses"], "history_sha256": digest, "distinct_scalar_evolution_calls": len(histories),
            "unique_scalar_samples": sum(len(h.samples) for h in histories.values()), "exact_input_matching": config["exact_input_matching"],
            "no_re_evolution_for_reused_patches": True, "native_cache_diagnostics": document["cache"]})
        if make_plots:
            _plots(config, scalar_rows, stage)
        code_paths = ("tools/run_thermal_width_power_reference.py", "mhdlab/thermal_width_power_reference.py", "mhdlab/thermal_forced_feedback_reference.py", "mhdlab/thermal_power_reference.py", "mhdlab/thermal_feedback_reference.py")
        report = {"schema": REPORT_SCHEMA, "status": "manufactured_reference_only", "config_original_path": str(config_path.resolve()), "config_sha256": sha(config_path),
            "code_sha256": {p: sha(ROOT/p) for p in code_paths}, "parameter_status": "manufactured", "scope": config["scope"],
            "distinct_scalar_evolution_calls": len(histories), "distinct_scalar_history_count": len(histories), "unique_scalar_samples": sum(len(h.samples) for h in histories.values()),
            "scalar_histories": records, "expanded_record_counts": counts, "reuse_map": config["hypotheses"], "exact_input_matching": config["exact_input_matching"],
            "assembly_file": "width_assembly.json", "assembly_sha256": sha(stage/"width_assembly.json"),
            "time_convention": config["time_convention"], "geometry_convention": config["geometry_convention"],
            "sample_units":"SI suffixes in scalar fields; areal quantities remain per square metre; extensive quantities use each patch area once. Exact ratios and signed projection diagnostics are preserved next to every extensive CSV value and in the full assembly JSON.",
            "csv_projection_contract":"Columns ending __projection_json contain complete JSON diagnostics; __exact_numerator/denominator are canonical integer strings. Logical time and area JSON retain exact ratios.",
            "artifact_path_base": "directory containing report.json", "artifact_sha256": {p.name: sha(p) for p in stage.iterdir() if p.is_file()},
            "qualifications": {"three_new_distinct_scalar_histories_only": True, "no_re_evolution_for_reused_patch": True, "exact_matched_total_input_without_rescaling": True,
                "full_time_state_temperature_mixture_retained": True, "instantaneous_flux_from_each_unique_saved_coupled_state": True,
                "aggregate_or_difference_sign_assumed": False, "independent_reference_error_computed_here": False, "old_case_or_old_producer_called": False,
                "physical_water_monolayer_material_or_spectral_inference": False, "actual_MHD_Hall_FLASH_host_gap_or_owner_coupling": False},
            "limits": ["Manufactured independent stationary patches, supplied constant areal capacity and positive-base-barrier kinetics only; no lateral conduction or field-derived heating.",
                       "Equal total input is exact for represented geometry and powers, not a physical experimental equivalence.",
                       "Instantaneous source flux and expected translation retain local temperatures; no mean-temperature replacement, particle launch, species mass or velocity distribution is supplied.",
                       "Actual hypothesis differences are saved with signed projections. Independent trajectory, flux and channel error qualification is separate; numerical differences alone do not establish physical observability."]}
        _write_json(stage/"report.json", report)
        stage.rename(output)
    except Exception as exc:
        _write_json(stage/"failure.json", {"status": "incomplete_reference_not_admitted", "error_type": type(exc).__name__, "error": str(exc),
            "completed_distinct_scalar_histories": len(histories), "no_fallback_or_retry_performed": True})
        raise ValueError("reference failed closed; incomplete evidence retained at "+str(stage)+": "+str(exc)) from exc
    print(json.dumps({"report": str(output/"report.json"), "report_sha256": sha(output/"report.json"), "distinct_scalar_histories": len(histories),
        "unique_scalar_samples": report["unique_scalar_samples"], "expanded_record_counts": counts, "status": report["status"]}, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"configs/thermal_width_power_manufactured.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    try:
        run(args.config, args.output_dir, make_plots=args.plot)
    except (ValueError, TypeError, OverflowError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
