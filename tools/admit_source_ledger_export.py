#!/usr/bin/env python3
"""Check pinned source-ledger artifacts and prepare an owner handoff."""
from __future__ import annotations

import argparse
import ast
import csv
from fractions import Fraction
import hashlib
import itertools
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.window_source_ledger_validation import validate_window_source_ledger


CODE_PATHS = (
    "tools/admit_source_ledger_export.py", "mhdlab/window_source_ledger_validation.py",
    "tools/export_window_source_ledger.py", "tools/export_packet_source_ledger.py",
    "tools/compare_windowed_source_resolution.py", "tools/inspect_windowed_source_moments.py",
    "tools/sample_temkin_source_windows.py",
)
STATUSES = {"manufactured", "reconstructed", "sourced", "unvalidated", "owner_supplied"}
MAX_JSON_BYTES = 64 * 1024 * 1024
PROJECTION_FLAGS = ("exact_zero", "positive", "negative", "finite", "underflowed", "overflowed", "rounded", "subnormal")
GROUP_KEYS = ("width_index", "current_index", "state_id", "reservoir_index", "reservoir_name", "selection_index")
QUANTITY_UNITS = {"released_number": "parent", "translational_energy_j": "J", "internal_energy_j": "J",
                  "binding_energy_j": "J", "total_wall_energy_debit_j": "J"}
CSV_FIELDS = (
    "case_id", "window_index", "window_start_s", "window_stop_s", "hypothesis", "group_type", "group_index",
    *GROUP_KEYS, "quantity", "unit", "per_area_unit", "area_m2", *("area_" + n for n in PROJECTION_FLAGS),
    "available", "availability_reason", "value", "projection_reason", *PROJECTION_FLAGS,
    "unqualified_count", "unqualified_source_group_ref", "qualified_subtotal_value",
    *("qualified_subtotal_" + n for n in PROJECTION_FLAGS),
    "per_area_available", "per_area_availability_reason", "per_area_value", "per_area_projection_reason",
    *("per_area_" + n for n in PROJECTION_FLAGS), "per_area_qualified_subtotal_value",
    *("per_area_qualified_subtotal_" + n for n in PROJECTION_FLAGS),
)
RESOLUTION_KEYS = {"width_bins", "current_bins", "snapshot_intervals", "subdivisions_per_interval"}
STUDY_CONFIG_KEYS = {"schema", "parameter_status", "conditioning_provenance", "comparison_provenance",
                     "relative_time_tolerance", "log_mass_tolerance", "windows", "orders", "cases", "comparisons"}
CASE_KEYS = {"case_id", "source_file", "source_sha256", "evidence_file", "evidence_sha256", "resolution"}
PACKET_CONFIG_KEYS = {"schema", "parameter_status", "conditioning_provenance", "mass_status", "mass_provenance",
                     "state_masses_kg", "selections", "log_mass_tolerance", "face_normal_coordinate_m",
                     "outward_normal_sign", "packet_count", "seed"}
REQUIRED_OWNER_CONTEXT = {name: "pending_owner_supplied" for name in (
    "physical_emission_kernel", "electromagnetic_fields", "hall_feedback", "loss_and_boundary_collection",
    "state_kinetics_and_response", "owner_precision_and_acceptance_policy")}


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def digest(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def safe_name(value):
    require(isinstance(value, str) and value and "\\" not in value, "artifact paths must be relative POSIX names")
    path = PurePosixPath(value)
    require(not path.is_absolute() and path.as_posix() == value and all(p not in (".", "..") for p in path.parts),
            "artifact path traversal or noncanonical name is unsupported")
    return value


def regular_file(path):
    path = Path(path)
    require(not any(p.is_symlink() for p in (path, *path.parents)) and path.is_file(),
            f"artifact must be a regular file without symbolic links: {path}")
    return path


def sha(path):
    result = hashlib.sha256()
    with regular_file(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def nonfinite(value):
        raise ValueError(f"nonfinite JSON constant: {value}")

    def finite_float(value):
        result = float(value)
        require(math.isfinite(result), "nonfinite JSON floating-point literal")
        return result
    with regular_file(path).open("rb") as stream:
        raw = stream.read(MAX_JSON_BYTES+1)
    require(len(raw) <= MAX_JSON_BYTES, "JSON artifact exceeds the admission byte limit")
    try:
        return json.loads(raw, object_pairs_hook=pairs, parse_constant=nonfinite, parse_float=finite_float)
    except RecursionError as error:
        raise ValueError("JSON artifact exceeds the supported nesting depth") from error


def same_json(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return set(a) == set(b) and all(same_json(value, b[key]) for key, value in a.items())
    if isinstance(a, list):
        return len(a) == len(b) and all(same_json(x, y) for x, y in zip(a, b))
    if isinstance(a, float):
        return math.isfinite(a) and math.isfinite(b) and a.hex() == b.hex()
    return a == b


def real(value, name, *, positive=False):
    require(type(value) in (int, float), f"{name} must be a finite real scalar")
    result = float(value)
    require(math.isfinite(result) and (not positive or result > 0), f"{name} is outside its finite domain")
    return result


def text(value, name):
    require(isinstance(value, str) and bool(value.strip()), f"{name} must be nonempty")
    return value


def hash_map(value, expected=None):
    require(isinstance(value, dict) and value and all(digest(v) for v in value.values()), "invalid artifact hash map")
    for name in value:
        safe_name(name)
    if expected is not None:
        require(set(value) == set(expected), "artifact hashes do not bind exactly the expected local files")
    return value


def literal_paths(path):
    """Read a fixed public producer declaration; never execute its module."""
    tree = ast.parse(regular_file(ROOT / path).read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id == "CODE_PATHS" for t in n.targets)]
    require(len(nodes) == 1, "known producer CODE_PATHS declaration is ambiguous")
    values = ast.literal_eval(nodes[0].value)
    require(isinstance(values, tuple) and values and len(set(values)) == len(values), "invalid known producer code paths")
    for name in values:
        safe_name(name)
        require(PurePosixPath(name).parts[0] in ("mhdlab", "tools") and name.endswith(".py"),
                "known producer code must remain inside the public Python contract")
    return values


def producer_pins(kind):
    if kind == "study":
        outer_paths = set(literal_paths("tools/export_window_source_ledger.py"))
        nested_paths = set(literal_paths("tools/compare_windowed_source_resolution.py"))
    else:
        nested_paths = set(literal_paths("tools/sample_temkin_source_windows.py"))
        outer_paths = (set(literal_paths("tools/export_packet_source_ledger.py")) | nested_paths
                       | set(literal_paths("tools/export_window_source_ledger.py"))
                       | set(literal_paths("tools/inspect_windowed_source_moments.py")))
    outer = {name: sha(ROOT / name) for name in sorted(outer_paths)}
    nested = {name: sha(ROOT / name) for name in sorted(nested_paths)}
    admission = {name: sha(ROOT / name) for name in CODE_PATHS}
    return outer, nested, {**outer, **nested, **admission}


def false_flags(document, names):
    require(all(document.get(name) is False for name in names), "artifact escalates an acceptance or recomputation flag")


def validate_resolution(value):
    require(isinstance(value, dict) and set(value) == RESOLUTION_KEYS
            and all(type(v) is int and 1 <= v <= 100_000 for v in value.values()), "invalid declared source resolution")


def study_configuration(value):
    require(isinstance(value, dict) and set(value) == STUDY_CONFIG_KEYS
            and value["schema"] == "mhdlab-windowed-source-resolution-config-v1", "invalid preserved study configuration schema")
    require(value["parameter_status"] in STATUSES, "invalid preserved source status")
    for name in ("conditioning_provenance", "comparison_provenance"):
        text(value[name], name)
    for name in ("relative_time_tolerance", "log_mass_tolerance"):
        require(0 < real(value[name], name) < 1, "source tolerances must be strictly between zero and one")
    windows = value["windows"]
    require(isinstance(windows, list) and 1 <= len(windows) <= 16, "invalid preserved source windows")
    normalized = []
    for window in windows:
        require(isinstance(window, dict) and set(window) == {"start_s", "stop_s"}, "invalid global source-window record")
        start, stop = real(window["start_s"], "start_s"), real(window["stop_s"], "stop_s")
        require(start < stop, "global source windows must have positive duration")
        normalized.append({"start_s": start, "stop_s": stop})
    ordered = sorted(normalized, key=lambda w: w["start_s"])
    require(all(a["stop_s"] <= b["start_s"] for a, b in zip(ordered, ordered[1:])), "global source windows overlap")
    orders = value["orders"]
    require(isinstance(orders, list) and 2 <= len(orders) <= 5
            and all(type(v) is int and 2 <= v <= 4096 for v in orders)
            and all(a < b for a, b in zip(orders, orders[1:])), "invalid preserved quadrature orders")
    cases = value["cases"]
    require(isinstance(cases, list) and 1 <= len(cases) <= 16, "invalid preserved case count")
    identifiers = set()
    for case in cases:
        require(isinstance(case, dict) and set(case) == CASE_KEYS, "invalid preserved case fields")
        name = case["case_id"]
        require(isinstance(name, str) and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", name) is not None
                and name not in identifiers, "invalid or duplicate case identifier")
        identifiers.add(name)
        text(case["source_file"], "original source path")
        text(case["evidence_file"], "original evidence path")
        require(digest(case["source_sha256"]) and digest(case["evidence_sha256"]), "invalid preserved case digest")
        validate_resolution(case["resolution"])
    require(isinstance(value["comparisons"], list) and len(value["comparisons"]) <= 64, "invalid comparison metadata")
    families = {"width": "width_bins", "current": "current_bins", "snapshot": "snapshot_intervals",
                "integration": "subdivisions_per_interval", "snapshot_common_intervals": "snapshot_intervals"}
    by_id = {c["case_id"]: c for c in cases}
    seen = set()
    for pair in value["comparisons"]:
        require(isinstance(pair, dict) and set(pair) == {"estimate_case_id", "reference_case_id", "family"},
                "invalid declared comparison fields")
        a, b, family = pair["estimate_case_id"], pair["reference_case_id"], pair["family"]
        require(isinstance(a, str) and isinstance(b, str) and a in by_id and b in by_id and a != b
                and isinstance(family, str) and family in families and (a, b, family) not in seen,
                "invalid declared source comparison")
        seen.add((a, b, family))
        ca, cb = by_id[a]["resolution"], by_id[b]["resolution"]
        varied = families[family]
        require(cb[varied] > ca[varied], "reference resolution must be finer")
        if family == "snapshot_common_intervals":
            require(all(ca[k] == cb[k] for k in ("width_bins", "current_bins"))
                    and ca[varied]*ca["subdivisions_per_interval"] == cb[varied]*cb["subdivisions_per_interval"],
                    "common-interval snapshots must retain their integration grid size")
        else:
            require(all(ca[k] == cb[k] for k in RESOLUTION_KEYS - {varied}), "comparison changes more than one resolution")
    return normalized


def opaque_digests(value):
    require(isinstance(value, dict) and value and all(isinstance(k, str) and k and digest(v) for k, v in value.items()),
            "invalid preserved provenance digests")


def bound_files(directory, report, expected):
    hashes = hash_map(report.get("output_sha256"), expected)
    for name, expected_hash in hashes.items():
        require(sha(directory / name) == expected_hash, f"bound artifact changed: {name}")
    return hashes


def source_labels(source, *, status, provenance, parameter_provenance):
    require(isinstance(source, dict) and source.get("schema") == "mhdlab-temkin-planar-source-windows-v1"
            and source.get("status") == status and source.get("provenance") == provenance
            and source.get("parameter_provenance") == parameter_provenance
            and source.get("host_budget_accepted") is False, "original conditioned source labels or provenance disagree")


def plan_binding(plan, source, windows, tolerance, index):
    require(isinstance(plan, dict) and plan.get("schema") == "mhdlab-source-window-plan-v1"
            and same_json(plan.get("windows"), windows) and same_json(plan.get("time_edges_s"), source.get("time_edges_s"))
            and same_json(plan.get("relative_time_tolerance"), float(tolerance)) and plan.get("units") == "s",
            "original plan, geometry or global windows disagree")
    false_flags(plan, ("source_emission_windows_are_photon_recording_gates", "physical_timing_uncertainty_estimated"))
    selections, indices = plan.get("selections"), plan.get("window_indices")
    require(isinstance(selections, list) and isinstance(indices, list) and len(selections) == len(indices)
            and all(type(i) is int and 0 <= i < len(windows) for i in indices), "invalid original plan selection ownership")
    require(same_json(source.get("selections"), [selection for selection, owner in zip(selections, indices) if owner == index]),
            "conditioned source selects different original plan intervals")
    edges = source["time_edges_s"]
    require(isinstance(edges, list) and len(edges) >= 2 and all(type(t) is float and math.isfinite(t) for t in edges)
            and all(a < b for a, b in zip(edges, edges[1:])), "invalid original plan time geometry")
    expected = []
    for owner, window in enumerate(windows):
        start, stop = Fraction(window["start_s"]), Fraction(window["stop_s"])
        require(Fraction(edges[0]) <= start < stop <= Fraction(edges[-1]), "global SOURCE window lacks complete support")
        for interval, (left, right) in enumerate(zip(edges, edges[1:])):
            a, b = max(start, Fraction(left)), min(stop, Fraction(right))
            if a < b:
                expected.append((owner, interval, a-Fraction(left), b-Fraction(left), b-a))
    require(len(expected) == len(selections), "original plan omits or adds a positive interval intersection")
    for selection, owner, (expected_owner, interval, start, stop, duration) in zip(selections, indices, expected):
        require(isinstance(selection, dict) and set(selection) == {"interval_index", "start_s", "stop_s"}
                and type(selection["interval_index"]) is int and selection["interval_index"] == interval
                and owner == expected_owner, "original plan interval ownership differs")
        a, b = real(selection["start_s"], "selection start"), real(selection["stop_s"], "selection stop")
        require(a >= 0 and a < b and abs(Fraction(a)-start) <= Fraction(float(tolerance))*duration
                and abs(Fraction(b)-stop) <= Fraction(float(tolerance))*duration,
                "original plan offset projection exceeds its declared geometric tolerance")


def evidence_binding(evidence, case, status, parameter_provenance):
    require(isinstance(evidence, dict) and evidence.get("schema") == "mhdlab-preserved-source-association-v1"
            and evidence.get("parameter_status") == status and evidence.get("parameter_provenance") == parameter_provenance
            and evidence.get("source_sha256") == case["source_sha256"]
            and same_json(evidence.get("resolution"), case["resolution"])
            and evidence.get("source_evolution_repeated") is False, "preserved source evidence association differs")
    opaque_digests(evidence.get("preserved_inputs_sha256"))


def study_inputs(directory, outer, hashes, producer, nested_producer):
    original = read_json(directory / "study/report.json")
    config = read_json(directory / "study/input_config.json")
    windows = study_configuration(config)
    require(isinstance(original, dict) and original.get("schema") == "mhdlab-windowed-source-resolution-report-v1",
            "invalid original study report schema")
    false_flags(original, ("source_evolution_repeated", "host_budget_accepted"))
    expected = {"input_config.json"}
    for case in config["cases"]:
        name = case["case_id"]
        expected.update((f"source_{name}.h5", f"evidence_{name}.json", f"case_{name}.json"))
    nested_hashes = hash_map(original.get("output_sha256"), expected)
    inputs = hash_map(outer.get("input_study_sha256"), expected | {"report.json"})
    require(set(hashes) == {"ledger.csv"} | {f"case_{c['case_id']}.json" for c in config["cases"]}
            | {"study/"+name for name in inputs}, "unexpected study-ledger export files")
    for name, value in inputs.items():
        require(hashes["study/"+name] == value and (name == "report.json" or nested_hashes[name] == value),
                "copied study hashes disagree across report levels")
    require(outer.get("study_report_sha256") == inputs["report.json"]
            and original.get("input_config_sha256") == inputs["input_config.json"], "original study/config pin differs")
    require(same_json(outer.get("source_sha256"), producer)
            and same_json(outer.get("study_source_sha256"), nested_producer)
            and same_json(original.get("source_sha256"), nested_producer), "study producer code pins differ from known current contracts")
    for document in (outer, original):
        for name in ("parameter_status", "conditioning_provenance", "comparison_provenance", "cases"):
            require(same_json(document.get(name), config[name]), f"study {name} association differs")
    require(same_json(original.get("windows"), config["windows"]) and same_json(original.get("orders"), config["orders"]),
            "original study windows or quadrature labels differ")
    require(same_json(outer.get("case_ledger_files"), {c["case_id"]: f"case_{c['case_id']}.json" for c in config["cases"]}),
            "ledger case identifiers or filenames differ")
    false_flags(outer, ("source_evolution_repeated", "source_moment_quadrature_repeated", "packet_sampling_performed", "host_budget_accepted"))
    return config, windows


def packet_configuration(value):
    require(isinstance(value, dict) and set(value) == PACKET_CONFIG_KEYS
            and value["schema"] == "mhdlab-temkin-windowed-packet-config-v1", "invalid original packet configuration")
    for name in ("parameter_status", "mass_status"):
        require(isinstance(value[name], str) and value[name] in STATUSES, "invalid original source or mass status")
    for name in ("conditioning_provenance", "mass_provenance"):
        text(value[name], name)
    require(type(value["packet_count"]) is int and 1 <= value["packet_count"] <= 1_000_000
            and type(value["seed"]) is int and 0 <= value["seed"] < 2**64, "invalid original packet count or seed")
    require(0 < real(value["log_mass_tolerance"], "log_mass_tolerance") < 1, "invalid original source tolerance")
    real(value["face_normal_coordinate_m"], "face_normal_coordinate_m")
    require(type(value["outward_normal_sign"]) is int and value["outward_normal_sign"] in (-1, 1), "invalid original face normal sign")
    masses = value["state_masses_kg"]
    require(isinstance(masses, dict) and masses, "original exact-state masses are absent")
    for state, mass in masses.items():
        text(state, "exact state ID")
        real(mass, "exact-state mass", positive=True)
    selections = value["selections"]
    require(isinstance(selections, list) and 1 <= len(selections) <= 10_000, "original packet selections are absent")
    result = []
    for selection in selections:
        require(isinstance(selection, dict) and set(selection) == {"interval_index", "start_s", "stop_s"}
                and type(selection["interval_index"]) is int and selection["interval_index"] >= 0,
                "invalid original packet selection")
        a, b = real(selection["start_s"], "selection start"), real(selection["stop_s"], "selection stop")
        require(0 <= a < b, "invalid original source selection offsets")
        result.append({"interval_index": selection["interval_index"], "start_s": a, "stop_s": b})
    return result


def packet_inputs(directory, outer, hashes, producer, nested_producer):
    original = read_json(directory / "packet_export/report.json")
    config = read_json(directory / "packet_export/input_config.json")
    selections = packet_configuration(config)
    require(isinstance(original, dict) and original.get("schema") == "mhdlab-temkin-windowed-packet-export-v1",
            "invalid original packet report schema")
    false_flags(original, ("host_budget_accepted",))
    required = {"input_config.json", "source_input.h5", "boundary_packets.npz", "uniforms.npy"}
    nested = original.get("output_sha256")
    require(isinstance(nested, dict) and set(nested) in (required, required | {"source_run_summary.json"}),
            "original packet files differ from the fixed contract")
    hash_map(nested)
    inputs = hash_map(outer.get("input_export_sha256"), set(nested) | {"report.json"})
    require(set(hashes) == {"ledger.json", "ledger.csv"} | {"packet_export/"+name for name in inputs},
            "unexpected packet-ledger export files")
    for name, value in inputs.items():
        require(hashes["packet_export/"+name] == value and (name == "report.json" or nested[name] == value),
                "copied packet hashes disagree across report levels")
    require(original.get("input_config_sha256") == inputs["input_config.json"]
            and original.get("input_source_sha256") == inputs["source_input.h5"], "original packet source/config pin differs")
    require(same_json(outer.get("source_sha256"), producer) and same_json(original.get("source_sha256"), nested_producer),
            "packet producer code pins differ from known current contracts")
    source = original.get("conditioned_source")
    require(isinstance(source, dict), "original packet conditioned source is absent")
    for document in (outer, original):
        require(document.get("parameter_status") == config["parameter_status"]
                and document.get("conditioning_provenance") == config["conditioning_provenance"]
                and document.get("source_parameter_provenance") == source.get("parameter_provenance"),
                "packet source status or provenance differs")
    source_labels(source, status=config["parameter_status"], provenance=config["conditioning_provenance"],
                  parameter_provenance=original["source_parameter_provenance"])
    require(same_json(source.get("selections"), selections), "original packet SOURCE selections differ")
    for outer_name, original_name in (("original_mass_status", "mass_status"), ("original_mass_provenance", "mass_provenance")):
        require(outer.get(outer_name) == original.get(original_name) == config[original_name], "saved mass provenance differs")
    metadata = original.get("packets")
    require(isinstance(metadata, dict) and same_json(metadata.get("source"), source)
            and metadata.get("status") == source["status"] and metadata.get("provenance") == source["provenance"]
            and metadata.get("parameter_provenance") == source["parameter_provenance"]
            and type(metadata.get("packets")) is int and metadata["packets"] == config["packet_count"]
            and same_json(metadata.get("face_normal_coordinate_m"), float(config["face_normal_coordinate_m"]))
            and same_json(metadata.get("outward_normal_sign"), config["outward_normal_sign"]),
            "original packet source or normal-geometry metadata differs")
    false_flags(metadata, ("host_budget_accepted",))
    active_states = {c["source_window"]["reservoir"]["state_id"] for c in source.get("components", [])}
    require(active_states <= set(config["state_masses_kg"]), "an original positive source state has no declared mass")
    require(same_json(metadata.get("state_masses_kg"), {s: float(config["state_masses_kg"][s]) for s in active_states}),
            "original active-state mass metadata differs")
    summary_hash = inputs.get("source_run_summary.json")
    require(original.get("input_summary_sha256") == summary_hash, "original source summary hash differs")
    if summary_hash is None:
        require(config["parameter_status"] == "unvalidated" and original.get("source_run_association") is None,
                "source without a hash-associated run summary must remain unvalidated")
    else:
        summary = read_json(directory / "packet_export/source_run_summary.json")
        require(isinstance(summary, dict) and summary.get("schema") in {
            "mhdlab-temkin-planar-run-v1", "mhdlab-joint-source-run-v1", "mhdlab-coverage-source-run-v1"}
            and summary.get("parameter_status") == config["parameter_status"]
            and summary.get("parameter_provenance") == source["parameter_provenance"], "original source summary provenance differs")
        opaque_digests(summary.get("output_sha256"))
        require(inputs["source_input.h5"] in summary["output_sha256"].values(), "original source summary does not bind its source bytes")
        require(same_json(original.get("source_run_association"), {"sha256": summary_hash,
                "parameter_status": config["parameter_status"], "parameter_provenance": source["parameter_provenance"],
                "schema": summary["schema"]}), "original source-run association differs")
    false_flags(outer, ("host_budget_accepted", "source_evolution_repeated", "packet_sampling_repeated", "source_moments_recomputed"))
    return config, source


def projection_columns(projection, prefix=""):
    return {prefix+name: "" if projection is None else int(projection[name]) for name in PROJECTION_FLAGS}


def projection_reason(projection):
    if projection is None:
        return "not_available"
    if projection["overflowed"]:
        return "si_projection_overflow"
    if projection["underflowed"]:
        return "si_projection_underflow"
    return "represented"


def quantity_columns(quantity, prefix=""):
    return {prefix+"available": int(quantity["available"]),
            prefix+"availability_reason": "complete_qualified_aggregation" if quantity["available"] else "unqualified_contributors",
            prefix+"value": quantity["value"], prefix+"projection_reason": projection_reason(quantity["projection"]),
            **projection_columns(quantity["projection"], prefix),
            prefix+"qualified_subtotal_value": quantity["qualified_subtotal"]["value"],
            **projection_columns(quantity["qualified_subtotal"]["projection"], prefix+"qualified_subtotal_")}


def csv_rows(ledger, context):
    """Independently reconstruct the producer's scalar CSV view of JSON."""
    window = context["window"]
    base = {"case_id": context["case_id"], "window_index": context["window_index"],
            "window_start_s": window["start_s"] if window is not None else "",
            "window_stop_s": window["stop_s"] if window is not None else ""}
    for hypothesis in ("resolved", "matched_uniform_width"):
        for group_type, collection in ledger[hypothesis].items():
            for group_index, group in enumerate([collection] if group_type == "total" else collection):
                area = group["area"]
                for name, quantity in group["quantities"].items():
                    row = {**base, "hypothesis": hypothesis, "group_type": group_type, "group_index": group_index,
                           **group["key"], "quantity": name, "unit": QUANTITY_UNITS[name],
                           "per_area_unit": QUANTITY_UNITS[name]+" m^-2", "area_m2": area["value_m2"],
                           **projection_columns(area["projection"], "area_"), **quantity_columns(quantity),
                           **quantity_columns(quantity["per_area"], "per_area_"),
                           "unqualified_count": quantity["unqualified"]["count"],
                           "unqualified_source_group_ref": quantity["unqualified"].get("source_group_ref", "")}
                    yield {field: "" if row.get(field) is None else str(row.get(field, "")) for field in CSV_FIELDS}


def reference(path, hashes, pointer=""):
    return {"path": "export/"+path, "sha256": hashes[path], "json_pointer": pointer}


def descriptor(ledger, source, summary, context, hashes, *, ledger_path, ledger_pointer,
               source_path, source_pointer, original_source_path, evidence_path, original_report_path, packet_context):
    require(summary.get("algebraically_consistent") is True and summary.get("source_association_bound") is True
            and summary.get("host_budget_accepted") is False, "ledger algebra or original-source binding was not admitted")
    return {**context, "label_origin": ("ledger_export_context" if context["window_index"] is None
                                       else "preserved_study_case_and_global_window"),
            "ledger": reference(ledger_path, hashes, ledger_pointer),
            "original_conditioned_source": reference(source_path, hashes, source_pointer),
            "original_source": reference(original_source_path, hashes),
            "original_evidence": reference(evidence_path, hashes) if evidence_path is not None else None,
            "original_report": reference(original_report_path, hashes),
            "component_joint_source": {
                "ledger_components_pointer": ledger_pointer+"/components",
                "source_components_pointer": source_pointer+"/components",
                "component_count": summary["component_count"], "component_records_preserved": True,
                "state_temperature_timing_are_joint": True,
                "interpretation": "Absolute JSON pointers within the two referenced files retain original interval origins, local windows, exact state/pool identities, temperature and scalar source-law declarations."},
            "source_association": ledger["source_association"], "geometry": ledger["geometry"],
            "resolved_total": ledger["resolved"]["total"],
            "matched_uniform_width_total": ledger["matched_uniform_width"]["total"],
            "width_control": "Only width position is redistributed; original current/state/pool/time/temperature and independent energy-cost mixtures remain joint in the component references.",
            "quantity_availability_meaning": "Complete qualified source aggregation, independent of finite SI projection or host readiness.",
            "unqualified_contributors_reference": reference(ledger_path, hashes, ledger_pointer+"/unqualified_contributors"),
            "packet_context": packet_context, "validation_summary": summary,
            "source_quantity_qualification_recomputed": False, "physical_truth_verified": False,
            "host_budget_accepted": False, "owner_implementation_verified": False,
            "required_owner_context": dict(REQUIRED_OWNER_CONTEXT)}


def study_ledgers(directory, outer, hashes, producer, nested_producer, descriptors):
    config, windows = study_inputs(directory, outer, hashes, producer, nested_producer)
    for case in config["cases"]:
        name = case["case_id"]
        ledger_path, source_path = f"case_{name}.json", f"study/case_{name}.json"
        source_file, evidence_file = f"study/source_{name}.h5", f"study/evidence_{name}.json"
        require(hashes[source_file] == case["source_sha256"] and hashes[evidence_file] == case["evidence_sha256"],
                "source/evidence byte hashes differ from the selected case")
        old = read_json(directory / source_path)
        document = read_json(directory / ledger_path)
        require(isinstance(old, dict) and isinstance(document, dict)
                and document.get("schema") == "mhdlab-window-source-ledger-case-v1"
                and document.get("case_id") == name and document.get("parameter_status") == config["parameter_status"]
                and document.get("conditioning_provenance") == config["conditioning_provenance"]
                and document.get("source_sha256") == case["source_sha256"]
                and document.get("evidence_sha256") == case["evidence_sha256"]
                and same_json(document.get("resolution"), case["resolution"])
                and same_json(document.get("plan"), old.get("plan")), "case ledger labels, source bytes, resolution or plan differ")
        previous, current = old.get("windows"), document.get("windows")
        require(isinstance(previous, list) and isinstance(current, list) and len(previous) == len(current) == len(windows),
                "case ledger omits or adds a global window")
        evidence = read_json(directory / evidence_file)
        evidence_binding(evidence, case, config["parameter_status"], document.get("parameter_provenance"))
        for index, (saved, item, window) in enumerate(zip(previous, current, windows)):
            require(isinstance(saved, dict) and isinstance(item, dict)
                    and type(saved.get("window_index")) is int and type(item.get("window_index")) is int
                    and saved["window_index"] == item["window_index"] == index
                    and same_json(saved.get("window"), window) and same_json(item.get("window"), window),
                    "selected case/global window identity differs")
            source = saved.get("conditioned_source")
            source_labels(source, status=config["parameter_status"], provenance=config["conditioning_provenance"],
                          parameter_provenance=document["parameter_provenance"])
            ledger = item.get("ledger")
            summary = validate_window_source_ledger(ledger, expected_source=source)
            plan_binding(old["plan"], source, windows, config["relative_time_tolerance"], index)
            counts = case["resolution"]
            require(len(source["x_edges_m"])-1 == counts["width_bins"]
                    and len(source["z_edges_m"])-1 == counts["current_bins"]
                    and len(source["time_edges_s"])-1 == counts["snapshot_intervals"]*counts["subdivisions_per_interval"],
                    "JSON source geometry differs from declared resolution")
            context = {"case_id": name, "window_index": index, "window": window}
            descriptors.append(descriptor(ledger, source, summary, context, hashes,
                ledger_path=ledger_path, ledger_pointer=f"/windows/{index}/ledger",
                source_path=source_path, source_pointer=f"/windows/{index}/conditioned_source",
                original_source_path=source_file, evidence_path=evidence_file, original_report_path="study/report.json",
                packet_context={"available": False, "reason": "no_packet_mass_or_normal_geometry_in_study_export",
                                "owner_supplied_required": True}))
            yield from csv_rows(ledger, context)


def packet_ledger(directory, outer, hashes, producer, nested_producer, descriptors):
    config, source = packet_inputs(directory, outer, hashes, producer, nested_producer)
    ledger = read_json(directory / "ledger.json")
    summary = validate_window_source_ledger(ledger, expected_source=source)
    context = {"case_id": "packet_source", "window_index": None, "window": None}
    packet_context = {"available": True, "config_reference": reference("packet_export/input_config.json", hashes),
                      "arrays": [reference("packet_export/"+name, hashes) for name in ("boundary_packets.npz", "uniforms.npy")],
                      "mass_status": config["mass_status"], "mass_provenance": config["mass_provenance"],
                      "state_masses_kg": config["state_masses_kg"],
                      "face_normal_coordinate_m": config["face_normal_coordinate_m"],
                      "outward_normal_sign": config["outward_normal_sign"], "realization_recomputed": False,
                      "array_contents_parsed": False, "normal_geometry_is_supplied_context": True}
    evidence = "packet_export/source_run_summary.json" if "packet_export/source_run_summary.json" in hashes else None
    descriptors.append(descriptor(ledger, source, summary, context, hashes,
        ledger_path="ledger.json", ledger_pointer="", source_path="packet_export/report.json",
        source_pointer="/conditioned_source", original_source_path="packet_export/source_input.h5",
        evidence_path=evidence, original_report_path="packet_export/report.json", packet_context=packet_context))
    yield from csv_rows(ledger, context)


def validate_csv(path, expected_rows):
    count = 0
    with regular_file(path).open(newline="") as stream:
        rows = csv.DictReader(stream, strict=True)
        require(rows.fieldnames == list(CSV_FIELDS), "ledger CSV header differs from the scalar contract")
        for actual, expected in itertools.zip_longest(rows, expected_rows):
            require(actual is not None and expected is not None and actual == expected,
                    "ledger CSV scalar rows differ from typed JSON quantities or projection flags")
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", dest="export_dir", type=Path, required=True)
    parser.add_argument("--expected-report-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not digest(args.expected_report_sha256):
        parser.error("--expected-report-sha256 must be an explicit trusted lowercase 64-character SHA256")
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        directory = args.export_dir.absolute()
        require(directory.is_dir() and not any(p.is_symlink() for p in (directory, *directory.parents)),
                "--export must be a direct local directory")
        require(sha(directory / "report.json") == args.expected_report_sha256, "trusted expected report SHA256 mismatch")
        outer = read_json(directory / "report.json")
        require(isinstance(outer, dict), "export report must be a JSON object")
        schemas = {"mhdlab-window-source-ledger-export-v1": "study", "mhdlab-packet-source-ledger-export-v1": "packet"}
        require(isinstance(outer.get("schema"), str) and outer["schema"] in schemas, "unsupported ledger export schema")
        kind = schemas[outer["schema"]]
        hashes = hash_map(outer.get("output_sha256"))
        require("report.json" not in hashes, "outer report must not self-declare its trusted pin")
        for name, expected in hashes.items():
            require(sha(directory / name) == expected, f"bound export artifact changed: {name}")
        input_hashes = {**hashes, "report.json": args.expected_report_sha256}
        producer, nested_producer, code_hashes = producer_pins(kind)
        descriptors = []
        rows = (study_ledgers if kind == "study" else packet_ledger)(directory, outer, hashes, producer, nested_producer, descriptors)
        csv_count = validate_csv(directory / "ledger.csv", rows)
        require(descriptors, "ledger export contains no selected source descriptors")
        common = {"expected_report_sha256": args.expected_report_sha256, "input_export_schema": outer["schema"],
                  "parameter_status": outer["parameter_status"], "conditioning_provenance": outer["conditioning_provenance"],
                  "host_budget_accepted": False, "owner_implementation_verified": False, "physical_truth_verified": False,
                  "source_quantity_qualification_recomputed": False}
        handoff = {"schema": "mhdlab-source-ledger-owner-handoff-v1", **common, "sources": descriptors,
                   "required_owner_context": dict(REQUIRED_OWNER_CONTEXT),
                   "limits": [
                       "The approved external report pin anchors byte integrity and source association; it is not an authenticity signature or physical qualification.",
                       "Source quantity qualification declarations are bound but not mathematically requalified; only saved JSON algebra and CSV scalar views are checked.",
                       "Quantity availability denotes source completeness, independent of SI projection range, owner precision policy or transaction readiness.",
                       "Original joint state/pool/temperature/current/time source components remain authoritative references; marginals are not an emission kernel.",
                       "The width-only control retains original joint mixtures. Expected 2*kB*T translation and supplied independent internal/removal channels remain separate.",
                       "HDF, NPZ and uniform arrays are hash-bound opaque bytes; array contents, physical source laws and saved packet realizations are not parsed or recomputed.",
                       "The owner must supply and verify physical kernels, fields, Hall feedback, losses/collection, state kinetics/response and numerical acceptance policy.",
                       "No owner implementation, host transaction, source evolution, molecular transport, radiation or photon-gate response is admitted.",
                   ]}
        admission = {"schema": "mhdlab-source-ledger-export-admission-v1", **common,
                     "input_export_sha256": input_hashes, "source_sha256": code_hashes,
                     "producer_code_sha256": producer, "nested_producer_code_sha256": nested_producer,
                     "algebraically_consistent": True, "csv_scalar_rows_verified": csv_count,
                     "sources": [{"case_id": s["case_id"], "window_index": s["window_index"],
                                  "validation_summary": s["validation_summary"]} for s in descriptors],
                     "runtime": {"python": sys.version, "implementation": sys.implementation.name,
                                 "third_party_imports_required": False},
                     "limits": handoff["limits"] + ["Staged publication assumes a coordinated single writer; no concurrent transaction protocol is provided."]}
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".source-ledger-admission-", dir=args.output_dir.parent) as stage:
            output = Path(stage) / "result"
            output.mkdir()
            for name in input_hashes:
                target = output / "export" / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(directory / name, target)
            (output / "owner_handoff.json").write_text(json.dumps(handoff, indent=2, allow_nan=False)+"\n")
            for name, expected in input_hashes.items():
                require(sha(directory / name) == expected and sha(output / "export" / name) == expected,
                        "original or copied export changed during admission")
            require(all(sha(ROOT / name) == expected for name, expected in code_hashes.items()),
                    "implementation or known producer contract changed during admission")
            admission["output_sha256"] = {p.relative_to(output).as_posix(): sha(p) for p in output.rglob("*") if p.is_file()}
            (output / "admission.json").write_text(json.dumps(admission, indent=2, allow_nan=False)+"\n")
            require(not args.output_dir.exists() and not args.output_dir.is_symlink(), "output appeared during admission; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, IndexError, OverflowError, RuntimeError, csv.Error, SyntaxError) as error:
        parser.error(str(error))
    print(json.dumps({"output_dir": str(args.output_dir), "sources": len(descriptors),
                      "host_budget_accepted": False, "owner_implementation_verified": False}))


if __name__ == "__main__":
    main()
