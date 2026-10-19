#!/usr/bin/env python3
"""Export quantity-qualified source ledgers from a preserved window study."""
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

from mhdlab.source_window_plan import plan_source_windows
from mhdlab.temkin_windowed_source import condition_temkin_planar_source
from mhdlab.window_source_ledger import window_source_ledger
from tools.compare_windowed_source_resolution import (
    CODE_PATHS as STUDY_CODE_PATHS, configuration, digest, load_case, read_json, sha, shared_inputs,
)


CODE_PATHS = (
    "tools/export_window_source_ledger.py", "mhdlab/window_source_ledger.py",
    "tools/compare_windowed_source_resolution.py", "mhdlab/source_window_plan.py",
    "mhdlab/window_time_moments.py", "mhdlab/windowed_joint_source_moments.py",
    "mhdlab/emission_time_moments.py", "mhdlab/temkin_emission_time.py",
    "mhdlab/temkin_windowed_source.py", "mhdlab/temkin_source_window.py",
    "mhdlab/temkin_emission_window.py", "mhdlab/temkin_planar.py",
    "mhdlab/temkin_planar_io.py", "mhdlab/temkin_history.py", "mhdlab/temkin_release.py",
    "mhdlab/temkin_io.py", "mhdlab/thermal_release.py", "mhdlab/planar_release.py",
    "mhdlab/planar_release_io.py", "mhdlab/surface_release_history.py",
    "mhdlab/flash_surface.py", "mhdlab/surface_sources.py", "mhdlab/hdf_integrity.py",
    "mhdlab/constants.py",
)
PROJECTION_FLAGS = ("exact_zero", "positive", "negative", "finite", "underflowed", "overflowed", "rounded", "subnormal")
GROUP_KEYS = ("width_index", "current_index", "state_id", "reservoir_index", "reservoir_name", "selection_index")
QUANTITY_UNITS = {
    "released_number": "parent", "translational_energy_j": "J", "internal_energy_j": "J",
    "binding_energy_j": "J", "total_wall_energy_debit_j": "J",
}
CSV_FIELDS = (
    "case_id", "window_index", "window_start_s", "window_stop_s", "hypothesis", "group_type", "group_index",
    *GROUP_KEYS, "quantity", "unit", "per_area_unit", "area_m2",
    *("area_" + name for name in PROJECTION_FLAGS),
    "available", "availability_reason", "value", "projection_reason", *PROJECTION_FLAGS,
    "unqualified_count", "unqualified_source_group_ref", "qualified_subtotal_value",
    *("qualified_subtotal_" + name for name in PROJECTION_FLAGS),
    "per_area_available", "per_area_availability_reason", "per_area_value", "per_area_projection_reason",
    *("per_area_" + name for name in PROJECTION_FLAGS),
    "per_area_qualified_subtotal_value", *("per_area_qualified_subtotal_" + name for name in PROJECTION_FLAGS),
)


def regular_file(path):
    """Refuse indirect paths before reading the hash-associated local bytes."""
    path = Path(path)
    if any(part.is_symlink() for part in (path, *path.parents)) or not path.is_file():
        raise ValueError(f"input must be a regular file without symbolic links: {path}")
    return path


def same_json(a, b):
    """Compare JSON associations without treating booleans as numeric counts."""
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return set(a) == set(b) and all(same_json(value, b[key]) for key, value in a.items())
    if isinstance(a, list):
        return len(a) == len(b) and all(same_json(x, y) for x, y in zip(a, b))
    return a == b


def validate_study(directory):
    """Read exact study inputs and associations without replaying prior moments."""
    directory = Path(directory).absolute()
    if any(part.is_symlink() for part in (directory, *directory.parents)) or not directory.is_dir():
        raise ValueError("--study must be a local directory without symbolic links")
    report_path = regular_file(directory / "report.json")
    report_hash = sha(report_path)
    report = read_json(report_path)
    if (not isinstance(report, dict) or report.get("schema") != "mhdlab-windowed-source-resolution-report-v1"
            or report.get("source_evolution_repeated") is not False or report.get("host_budget_accepted") is not False):
        raise ValueError("study must be a preserved source-resolution report with no accepted host transaction")
    config_path = regular_file(directory / "input_config.json")
    raw, config, windows = configuration(config_path)
    expected_files = {"input_config.json"}
    for case in config["cases"]:
        name = case["case_id"]
        expected_files.update((f"source_{name}.h5", f"evidence_{name}.json", f"case_{name}.json"))
    hashes = report.get("output_sha256")
    if (not isinstance(hashes, dict) or set(hashes) != expected_files
            or any(not digest(value) for value in hashes.values())):
        raise ValueError("study output hashes must bind exactly the expected local input, source, evidence and case files")
    for name, expected in hashes.items():
        if sha(regular_file(directory / name)) != expected:
            raise ValueError(f"study file is changed: {name}")
    if (report.get("input_config_sha256") != hashlib.sha256(raw).hexdigest()
            or report.get("parameter_status") != config["parameter_status"]
            or report.get("conditioning_provenance") != config["conditioning_provenance"]
            or report.get("comparison_provenance") != config["comparison_provenance"]
            or not same_json(report.get("orders"), config["orders"])
            or not same_json(report.get("windows"), config["windows"])
            or not same_json(report.get("cases"), config["cases"])):
        raise ValueError("study report disagrees with its exact configuration, status or provenance")
    current_study_code = {name: sha(ROOT / name) for name in STUDY_CODE_PATHS}
    if report.get("source_sha256") != current_study_code:
        raise ValueError("study implementation differs from the pinned current source contract")
    comparisons = report.get("source_comparisons")
    if not isinstance(comparisons, list) or len(comparisons) != len(config["comparisons"]):
        raise ValueError("study comparison records disagree with the configuration")
    expected_pairs = [(index, order) for index in range(len(windows)) for order in config["orders"]]
    for saved, pair in zip(comparisons, config["comparisons"]):
        if (not isinstance(saved, dict) or any(saved.get(name) != value for name, value in pair.items())
                or not isinstance(saved.get("differences"), list)
                or any(not isinstance(item, dict) or type(item.get("window_index")) is not int
                       or type(item.get("quadrature_order")) is not int for item in saved["differences"])
                or [(item["window_index"], item["quadrature_order"]) for item in saved["differences"]] != expected_pairs):
            raise ValueError("study comparison cases, windows or orders disagree")
    loaded = {}
    for case in config["cases"]:
        name = case["case_id"]
        if (hashes[f"source_{name}.h5"] != case["source_sha256"]
                or hashes[f"evidence_{name}.json"] != case["evidence_sha256"]):
            raise ValueError("study source or evidence hash disagrees with its case declaration")
        local_case = {**case, "source_file": str(directory / f"source_{name}.h5"),
                      "evidence_file": str(directory / f"evidence_{name}.json")}
        evidence = read_json(directory / f"evidence_{name}.json")
        if not isinstance(evidence, dict) or not same_json(evidence.get("resolution"), case["resolution"]):
            raise ValueError("study evidence resolution types or counts disagree with the case declaration")
        loaded[name] = load_case(config_path, local_case, config["parameter_status"])[0]
    for pair in config["comparisons"]:
        shared_inputs(loaded[pair["estimate_case_id"]], loaded[pair["reference_case_id"]], pair["family"])
    return directory, report, config, windows, loaded, {**hashes, "report.json": report_hash}


def inspect_case(directory, case, history, config, windows):
    """Reproduce each saved source association, then form a new quantity ledger."""
    name = case["case_id"]
    saved = read_json(directory / f"case_{name}.json")
    plan = plan_source_windows(history.time_edges_s, windows, relative_time_tolerance=config["relative_time_tolerance"])
    if (not isinstance(saved, dict) or not same_json(saved.get("plan"), plan.manifest())
            or not isinstance(saved.get("windows"), list) or len(saved["windows"]) != len(windows)):
        raise ValueError(f"saved source-window plan differs from the canonical history: {name}")
    results = []
    for index, (window, previous) in enumerate(zip(windows, saved["windows"])):
        window_record = {"start_s": window.start_s, "stop_s": window.stop_s}
        if (not isinstance(previous, dict) or type(previous.get("window_index")) is not int
                or previous["window_index"] != index or not same_json(previous.get("window"), window_record)):
            raise ValueError("saved source-window indices or global endpoints differ")
        estimates = previous.get("estimates")
        if (not isinstance(estimates, list) or len(estimates) != len(config["orders"])
                or any(not isinstance(estimate, dict)
                       or estimate.get("schema") != "mhdlab-windowed-joint-source-moments-v1"
                       or type(estimate.get("quadrature_order")) is not int or estimate["quadrature_order"] != order
                       or estimate.get("status") != config["parameter_status"]
                       or estimate.get("provenance") != config["conditioning_provenance"]
                       or estimate.get("parameter_provenance") != history.parameter_provenance
                       or estimate.get("host_budget_accepted") is not False
                       for estimate, order in zip(estimates, config["orders"]))):
            raise ValueError("saved moment metadata differs from its source status, provenance or orders")
        selections = tuple(selection for selection, owner in zip(plan.selections, plan.window_indices) if owner == index)
        source = condition_temkin_planar_source(history, selections, log_mass_tolerance=config["log_mass_tolerance"],
                                                status=config["parameter_status"], provenance=config["conditioning_provenance"])
        if not same_json(previous.get("conditioned_source"), source.manifest()):
            raise ValueError("saved conditioned source differs from its canonical original history windows")
        ledger = window_source_ledger(source, include_components=True)
        results.append({"window_index": index, "window": window_record, "ledger": ledger})
    return {"schema": "mhdlab-window-source-ledger-case-v1", "case_id": name,
            "parameter_status": config["parameter_status"], "parameter_provenance": history.parameter_provenance,
            "conditioning_provenance": config["conditioning_provenance"], "source_sha256": case["source_sha256"],
            "evidence_sha256": case["evidence_sha256"], "resolution": case["resolution"],
            "plan": plan.manifest(), "windows": results}


def numeric(value):
    if value is not None and (type(value) not in (int, float) or not math.isfinite(value)):
        raise ValueError("ledger CSV values must be finite real scalars or unavailable")
    return "" if value is None else value


def projection_columns(projection, prefix=""):
    if projection is None:
        return {prefix + name: "" for name in PROJECTION_FLAGS}
    if not isinstance(projection, dict) or any(type(projection.get(name)) is not bool for name in PROJECTION_FLAGS):
        raise ValueError("ledger projection flags are incomplete")
    return {prefix + name: int(projection[name]) for name in PROJECTION_FLAGS}


def projection_reason(projection):
    if projection is None:
        return "not_available"
    if projection["overflowed"]:
        return "si_projection_overflow"
    if projection["underflowed"]:
        return "si_projection_underflow"
    return "represented"


def quantity_columns(quantity, prefix=""):
    available = quantity["available"]
    if type(available) is not bool:
        raise ValueError("ledger quantity availability must be explicit")
    subtotal = quantity["qualified_subtotal"]
    return {prefix + "available": int(available),
            prefix + "availability_reason": "complete_qualified_aggregation" if available else "unqualified_contributors",
            prefix + "value": numeric(quantity["value"]),
            prefix + "projection_reason": projection_reason(quantity["projection"]),
            **projection_columns(quantity["projection"], prefix),
            prefix + "qualified_subtotal_value": numeric(subtotal["value"]),
            **projection_columns(subtotal["projection"], prefix + "qualified_subtotal_")}


def ledger_rows(records):
    for record in records:
        ledger = record["ledger"]
        window = record.get("window")
        context = {"case_id": record["case_id"], "window_index": record.get("window_index"),
                   "window_start_s": window["start_s"] if window is not None else "",
                   "window_stop_s": window["stop_s"] if window is not None else ""}
        for hypothesis in ("resolved", "matched_uniform_width"):
            for group_type, collection in ledger[hypothesis].items():
                groups = [collection] if group_type == "total" else collection
                if not isinstance(groups, list):
                    raise ValueError("ledger CSV groups must be a total record or lists of records")
                for group_index, group in enumerate(groups):
                    if not isinstance(group["key"], dict) or set(group["key"]) - set(GROUP_KEYS):
                        raise ValueError("ledger group contains an unsupported CSV identifier")
                    if set(group["quantities"]) != set(QUANTITY_UNITS):
                        raise ValueError("ledger quantity names differ from the CSV contract")
                    area = group["area"]
                    for name, quantity in group["quantities"].items():
                        unqualified = quantity["unqualified"]
                        count = unqualified["count"]
                        if type(count) is not int or count < 0:
                            raise ValueError("ledger unqualified contributor count must be a nonnegative integer")
                        yield {**context, "hypothesis": hypothesis, "group_type": group_type,
                               "group_index": group_index, **group["key"], "quantity": name,
                               "unit": QUANTITY_UNITS[name], "per_area_unit": QUANTITY_UNITS[name] + " m^-2",
                               "area_m2": numeric(area["value_m2"]), **projection_columns(area["projection"], "area_"),
                               **quantity_columns(quantity), **quantity_columns(quantity["per_area"], "per_area_"),
                               "unqualified_count": count,
                               "unqualified_source_group_ref": unqualified.get("source_group_ref", "")}


def write_ledger_csv(path, ledger_records):
    """Write group quantities; unknown totals remain distinct from SI underflow.

    Each input record has ``case_id``, optional ``window_index`` and ``window``
    (a start_s/stop_s mapping or None), and the core ``ledger``. Projection flags
    and availability use numeric 0/1 columns; absent projections remain blank.
    Exact rational quantities and contributor catalogues remain in ledger JSON.
    """
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(ledger_rows(ledger_records))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        code_hashes = {name: sha(ROOT / name) for name in CODE_PATHS}
        directory, study_report, config, windows, histories, input_hashes = validate_study(args.study)
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".window-source-ledger-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            copies = output / "study"
            copies.mkdir()
            for name in input_hashes:
                shutil.copyfile(directory / name, copies / name)

            def records():
                for case in config["cases"]:
                    name = case["case_id"]
                    result = inspect_case(directory, case, histories[name], config, windows)
                    (output / f"case_{name}.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
                    for item in result["windows"]:
                        yield {"case_id": name, **item}

            write_ledger_csv(output / "ledger.csv", records())
            report = {"schema": "mhdlab-window-source-ledger-export-v1", "parameter_status": config["parameter_status"],
                      "conditioning_provenance": config["conditioning_provenance"],
                      "comparison_provenance": config["comparison_provenance"],
                      "study_report_sha256": input_hashes["report.json"], "input_study_sha256": input_hashes,
                      "study_source_sha256": study_report["source_sha256"], "source_sha256": code_hashes,
                      "cases": config["cases"], "case_ledger_files": {case["case_id"]: f"case_{case['case_id']}.json" for case in config["cases"]},
                      "source_evolution_repeated": False, "source_moment_quadrature_repeated": False,
                      "packet_sampling_performed": False, "host_budget_accepted": False,
                      "runtime": {"python": sys.version, **{name: importlib.metadata.version(name) for name in ("numpy", "scipy", "h5py")}},
                      "limits": [
                          "Only preserved canonical histories and original SOURCE emission windows are used; photon recording gates remain separate.",
                          "Prior stored moment values are copied and hash-associated; their quadrature is not recomputed.",
                          "Every saved source plan and conditioned source manifest is replayed before a new quantity ledger is formed.",
                          "Quantity availability describes complete qualified aggregation; finite SI projection and positive underflow are separate flags.",
                          "CSV contains finite scalar projections and explicit flags; exact fractions and unqualified contributor evidence remain in JSON.",
                          "Supplied source status and provenance do not establish physical water parameters, coverage, calibration or desorption energy.",
                          "No source evolution, packet generation, host transaction, wall feedback, gap transport, collision, radiation or actual FLASH run is performed.",
                          "Staged output assumes a coordinated single writer; it does not provide a concurrent transaction protocol.",
                      ]}
            for name, expected in input_hashes.items():
                if sha(regular_file(directory / name)) != expected or sha(copies / name) != expected:
                    raise ValueError("study input changed during ledger export")
            if any(sha(ROOT / name) != expected for name, expected in code_hashes.items()):
                raise ValueError("implementation changed during ledger export")
            report["output_sha256"] = {path.relative_to(output).as_posix(): sha(path)
                                       for path in output.rglob("*") if path.is_file()}
            (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during export; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as error:
        parser.error(str(error))
    print(json.dumps({"output_dir": str(args.output_dir), "cases": len(config["cases"]),
                      "parameter_status": config["parameter_status"], "host_budget_accepted": False}))


if __name__ == "__main__":
    main()
