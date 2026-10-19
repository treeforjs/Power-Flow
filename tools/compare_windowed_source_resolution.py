#!/usr/bin/env python3
"""Compare fixed SOURCE windows in preserved canonical release histories."""
from __future__ import annotations

import argparse
from collections.abc import Mapping
from fractions import Fraction
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mhdlab.source_window_plan import SourceTimeWindow, plan_source_windows
from mhdlab.temkin_planar_io import load_temkin_planar_release_history
from mhdlab.temkin_source_window import STATUSES
from mhdlab.temkin_windowed_source import condition_temkin_planar_source
from mhdlab.windowed_joint_source_moments import windowed_joint_source_moments

CODE_PATHS = (
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
RESOLUTION_KEYS = {"width_bins", "current_bins", "snapshot_intervals", "subdivisions_per_interval"}
FAMILIES = {"width": "width_bins", "current": "current_bins",
            "snapshot": "snapshot_intervals", "integration": "subdivisions_per_interval",
            "snapshot_common_intervals": "snapshot_intervals"}
CONFIG_KEYS = {"schema", "parameter_status", "conditioning_provenance", "comparison_provenance",
               "relative_time_tolerance", "log_mass_tolerance", "windows", "orders", "cases", "comparisons"}
CASE_KEYS = {"case_id", "source_file", "source_sha256", "evidence_file", "evidence_sha256", "resolution"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def nonfinite(value):
        raise ValueError(f"nonfinite JSON constant: {value}")
    return json.loads(Path(path).read_bytes(), object_pairs_hook=pairs, parse_constant=nonfinite)


def digest(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def resolution(value):
    if (not isinstance(value, dict) or set(value) != RESOLUTION_KEYS
            or any(type(item) is not int or not 1 <= item <= 100_000 for item in value.values())):
        raise ValueError("resolution requires all four explicit positive integer counts in [1,100000]")
    return value


def configuration(path):
    raw = path.read_bytes()
    value = read_json(path)
    if (not isinstance(value, dict) or set(value) != CONFIG_KEYS
            or value["schema"] != "mhdlab-windowed-source-resolution-config-v1"):
        raise ValueError("configuration requires exactly the windowed-source-resolution-config-v1 fields")
    if not isinstance(value["parameter_status"], str) or value["parameter_status"] not in STATUSES:
        raise ValueError("parameter_status must retain a supported explicit label")
    for name in ("conditioning_provenance", "comparison_provenance"):
        if not isinstance(value[name], str) or not value[name].strip():
            raise ValueError(f"{name} must be nonempty")
    for name in ("relative_time_tolerance", "log_mass_tolerance"):
        item = value[name]
        if type(item) not in (int, float) or not math.isfinite(item) or not 0 < item < 1:
            raise ValueError(f"{name} must be finite and strictly between zero and one")
    windows = value["windows"]
    if (not isinstance(windows, list) or not 1 <= len(windows) <= 16
            or any(not isinstance(row, dict) or set(row) != {"start_s", "stop_s"} for row in windows)):
        raise ValueError("windows requires 1 to 16 explicit global SOURCE start/stop records")
    windows = tuple(SourceTimeWindow(**row) for row in windows)
    orders = value["orders"]
    if (not isinstance(orders, list) or not 2 <= len(orders) <= 5
            or any(type(order) is not int or not 2 <= order <= 4096 for order in orders)
            or any(a >= b for a, b in zip(orders, orders[1:]))):
        raise ValueError("orders requires 2 to 5 strictly increasing integers in [2,4096]")
    cases = value["cases"]
    if not isinstance(cases, list) or not 1 <= len(cases) <= 16:
        raise ValueError("cases requires 1 to 16 preserved canonical histories")
    ids = set()
    for case in cases:
        if not isinstance(case, dict) or set(case) != CASE_KEYS:
            raise ValueError("case requires exact source, evidence, hash and resolution fields")
        name = case["case_id"]
        if not isinstance(name, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", name) is None or name in ids:
            raise ValueError("case_id must be unique and contain only safe identifier characters")
        ids.add(name)
        resolution(case["resolution"])
        for field in ("source_file", "evidence_file"):
            if not isinstance(case[field], str) or not case[field].strip():
                raise ValueError(f"{field} must be a nonempty path")
        if not digest(case["source_sha256"]) or not digest(case["evidence_sha256"]):
            raise ValueError("case source and evidence hashes must be explicit lowercase SHA256")
    comparisons = value["comparisons"]
    if not isinstance(comparisons, list) or len(comparisons) > 64:
        raise ValueError("comparisons must contain at most 64 declared pairs")
    seen = set()
    by_id = {case["case_id"]: case for case in cases}
    for pair in comparisons:
        if (not isinstance(pair, dict) or set(pair) != {"estimate_case_id", "reference_case_id", "family"}
                or pair["family"] not in FAMILIES):
            raise ValueError("comparison requires estimate_case_id, reference_case_id and one supported family")
        a, b = pair["estimate_case_id"], pair["reference_case_id"]
        if not isinstance(a, str) or not isinstance(b, str) or a not in ids or b not in ids or a == b:
            raise ValueError("comparison requires distinct existing case IDs")
        key = (a, b, pair["family"])
        if key in seen:
            raise ValueError("duplicate comparison")
        seen.add(key)
        varied = FAMILIES[pair["family"]]
        ca, cb = by_id[a]["resolution"], by_id[b]["resolution"]
        if pair["family"] == "snapshot_common_intervals":
            if (cb[varied] <= ca[varied] or any(ca[k] != cb[k] for k in ("width_bins", "current_bins"))
                    or ca[varied]*ca["subdivisions_per_interval"] != cb[varied]*cb["subdivisions_per_interval"]):
                raise ValueError("snapshot_common_intervals requires finer snapshots at a fixed total source interval count and spatial grid")
        elif cb[varied] <= ca[varied] or any(ca[k] != cb[k] for k in RESOLUTION_KEYS - {varied}):
            raise ValueError("comparison must refine exactly its declared resolution count")
    return raw, value, windows


def local_input(config_path, value, expected_hash):
    path = Path(value)
    if not path.is_absolute():
        path = config_path.parent / path
    if path.is_symlink() or not path.is_file() or sha(path) != expected_hash:
        raise ValueError(f"input is missing, indirect or changed: {value}")
    return path


def load_case(path, case, status):
    source_path = local_input(path, case["source_file"], case["source_sha256"])
    evidence_path = local_input(path, case["evidence_file"], case["evidence_sha256"])
    evidence = read_json(evidence_path)
    history = load_temkin_planar_release_history(source_path)
    if (not isinstance(evidence, dict) or evidence.get("schema") != "mhdlab-preserved-source-association-v1"
            or evidence.get("parameter_status") != status
            or evidence.get("parameter_provenance") != history.parameter_provenance
            or evidence.get("source_sha256") != case["source_sha256"]
            or evidence.get("resolution") != case["resolution"]
            or evidence.get("source_evolution_repeated") is not False
            or not isinstance(evidence.get("preserved_inputs_sha256"), dict)
            or not evidence["preserved_inputs_sha256"]
            or any(not isinstance(key, str) or not key or not digest(item)
                   for key, item in evidence["preserved_inputs_sha256"].items())):
        raise ValueError("preserved source evidence disagrees with the exact source, status, provenance or resolution")
    counts = case["resolution"]
    if (len(history.x_edges_m)-1 != counts["width_bins"]
            or len(history.z_edges_m)-1 != counts["current_bins"]
            or len(history.time_edges_s)-1 != counts["snapshot_intervals"]*counts["subdivisions_per_interval"]):
        raise ValueError("source dimensions disagree with declared resolution counts")
    remap = history.surface_attributes.get("followup_common_grid_resampling")
    if remap is not None:
        required = {"original_snapshot_intervals", "equivalent_subdivisions_per_original_interval",
                    "actual_integrator_snapshot_intervals", "actual_integrator_subdivisions",
                    "additional_recovered_information", "temporal_smoothing", "original_snapshot_attributes"}
        expected = {"original_snapshot_intervals": counts["snapshot_intervals"],
                    "equivalent_subdivisions_per_original_interval": counts["subdivisions_per_interval"],
                    "actual_integrator_snapshot_intervals": len(history.time_edges_s)-1,
                    "actual_integrator_subdivisions": 1}
        if (not isinstance(remap, Mapping) or set(remap) != required
                or any(type(remap[k]) is not int or remap[k] != v for k, v in expected.items())
                or remap["additional_recovered_information"] is not False or remap["temporal_smoothing"] is not False
                or not isinstance(remap["original_snapshot_attributes"], Mapping)):
            raise ValueError("common-grid resampling metadata disagrees with original/equivalent counts or actual integration")
    return history, source_path, evidence_path


def shared_inputs(a, b, family):
    # This deliberately supports uniform original areal pool inventories only.
    # A general remapped spatial inventory needs its own qualification.
    def metadata(history):
        result = dict(history.surface_attributes)
        if family == "current":
            # Canonical source validation already ties this derived list to z.
            # It changes when the declared current partition is refined.
            result.pop("current_patch_intervals_m", None)
        if family == "snapshot_common_intervals" and "followup_common_grid_resampling" in result:
            remap = dict(result["followup_common_grid_resampling"])
            # load_case already cross-checks these two counts against each case.
            # Every other transformation detail and original attribute remains.
            remap.pop("original_snapshot_intervals")
            remap.pop("equivalent_subdivisions_per_original_interval")
            result["followup_common_grid_resampling"] = remap
        return result
    if (a.reservoirs != b.reservoirs or a.parameter_provenance != b.parameter_provenance
            or metadata(a) != metadata(b)):
        raise ValueError("compared sources must retain identical pool parameters, provenance and surface metadata")
    for name in ("time_edges_s", "x_edges_m", "z_edges_m"):
        va, vb = getattr(a, name), getattr(b, name)
        if va[0] != vb[0] or va[-1] != vb[-1]:
            raise ValueError("compared sources require identical physical time and face domains")
        if (name == "x_edges_m" and family != "width") or (name == "z_edges_m" and family != "current"):
            if not np.array_equal(va, vb):
                raise ValueError("unvaried spatial coordinates must be identical")
        if name == "time_edges_s" and family in {"width", "current", "snapshot_common_intervals"} and not np.array_equal(va, vb):
            raise ValueError("spatial or common-interval snapshot comparisons must retain identical source time intervals")
    for pool in range(len(a.reservoirs)):
        va, vb = a.coverage_m2[0, pool], b.coverage_m2[0, pool]
        if not np.all(va == va.flat[0]) or not np.all(vb == va.flat[0]) or va.flat[0] != vb.flat[0]:
            raise ValueError("resolution comparison currently requires identical uniform initial areal inventory for every pool")


def rational_record(value):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def finite_change(value):
    try:
        projected = float(value)
    except OverflowError:
        projected = None
    finite = projected is not None and math.isfinite(projected)
    return {"value": projected if finite else None, "exact_represented_input_fraction": rational_record(value),
            "underflowed": bool(value and finite and projected == 0.), "overflowed": not finite}


def represented_mean_time(result, hypothesis):
    name = "mean_time_anchor" if hypothesis == "resolved" else "matched_uniform_width_mean_time_anchor"
    anchor = result[name]
    if anchor is None:
        return None
    return sum((Fraction(anchor[key]) for key in ("interval_origin_s", "window_start_offset_s", "mean_time_offset_s")), Fraction(0))


def changes(estimate, reference):
    """Estimate minus reference; time uses retained anchors across histories."""
    result = {}
    for hypothesis in ("resolved", "matched_uniform_width"):
        values = {}
        for name, unit in estimate["units"].items():
            if name == "mean_time_offset_s":
                continue  # Its anchor can differ between histories/orders.
            a, b = estimate[hypothesis][name], reference[hypothesis][name]
            if a is None or b is None:
                values[name] = {"unit": unit, "available": False, "reason": "empty_source_population"}
                continue
            if name == "mean_time_s":
                a, b = represented_mean_time(estimate, hypothesis), represented_mean_time(reference, hypothesis)
            else:
                a, b = Fraction(a), Fraction(b)
            difference = a-b
            values[name] = {"unit": unit, "available": True,
                            "signed_change": finite_change(difference),
                            "absolute_change": finite_change(abs(difference)),
                            "relative_to_reference": finite_change(difference/abs(b)) if b else None,
                            "relative_reason": "represented_reference" if b else "zero_reference"}
        result[hypothesis] = values
    return result


def inspect_case(history, config, windows):
    plan = plan_source_windows(history.time_edges_s, windows, relative_time_tolerance=config["relative_time_tolerance"])
    results = []
    for window_index, window in enumerate(windows):
        selections = tuple(selection for selection, index in zip(plan.selections, plan.window_indices) if index == window_index)
        source = condition_temkin_planar_source(history, selections, log_mass_tolerance=config["log_mass_tolerance"],
                                                status=config["parameter_status"], provenance=config["conditioning_provenance"])
        estimates = []
        for order in config["orders"]:
            estimate = windowed_joint_source_moments(source, quadrature_order=order, include_component_records=False)
            estimate.pop("source")  # Retained once below, rather than at every order.
            estimates.append(estimate)
        results.append({"window_index": window_index, "window": {"start_s": window.start_s, "stop_s": window.stop_s},
                        "conditioned_source": source.manifest(), "estimates": estimates,
                        "quadrature_changes": [{"coarse_order": a["quadrature_order"], "refined_order": b["quadrature_order"],
                                                "refined_minus_coarse": changes(b, a)} for a, b in zip(estimates, estimates[1:])]})
    return {"plan": plan.manifest(), "windows": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        parser.error("--output-dir must be a new path")
    try:
        code_hashes = {name: sha(ROOT/name) for name in CODE_PATHS}
        raw, config, windows = configuration(args.config)
        loaded = {case["case_id"]: load_case(args.config, case, config["parameter_status"]) for case in config["cases"]}
        for pair in config["comparisons"]:
            shared_inputs(loaded[pair["estimate_case_id"]][0], loaded[pair["reference_case_id"]][0], pair["family"])
        estimates = {name: inspect_case(value[0], config, windows) for name, value in loaded.items()}
        comparisons = []
        for pair in config["comparisons"]:
            ca, cb = (estimates[pair[name]] for name in ("estimate_case_id", "reference_case_id"))
            differences = [{"window_index": index, "quadrature_order": order, "estimate_minus_reference": changes(a, b)}
                           for index, (wa, wb) in enumerate(zip(ca["windows"], cb["windows"]))
                           for order, a, b in zip(config["orders"], wa["estimates"], wb["estimates"])]
            comparisons.append({**pair, "differences": differences})
        report = {"schema": "mhdlab-windowed-source-resolution-report-v1", "parameter_status": config["parameter_status"],
                  "conditioning_provenance": config["conditioning_provenance"], "comparison_provenance": config["comparison_provenance"],
                  "input_config_sha256": hashlib.sha256(raw).hexdigest(), "source_sha256": code_hashes,
                  "orders": config["orders"], "windows": config["windows"], "cases": config["cases"],
                  "source_comparisons": comparisons, "source_evolution_repeated": False, "host_budget_accepted": False,
                  "runtime": {"python": sys.version, **{name: importlib.metadata.version(name) for name in ("numpy", "scipy", "h5py")}},
                  "limits": [
                      "Fixed SOURCE emission windows are separate from photon recording gates and diagnostic timing origins.",
                      "Only preserved canonical histories are read and validated; no source evolution or packet generation is performed.",
                      "Evidence binds exact inputs and their declared labels; it does not independently certify their physical origin or all prior transformations.",
                      "One declared count varies per ordinary source comparison. snapshot_common_intervals compensates subdivision count and requires identical source time edges. Common physical drive association remains supplied provenance; no convergence order or certified error is inferred.",
                      "Finite quadrature changes are reported separately from source discretization differences at each fixed order; neither is experimental uncertainty.",
                      "Time differences use retained origins, window starts and projected local means before subtraction; arbitrary local anchors are never compared directly.",
                      "Compared original areal pool inventories must be uniform and identical; general spatial inventory remapping is unsupported.",
                      "No physical water coverage, absolute spectral calibration, gap, chemistry, radiation, host transaction or actual FLASH run is established.",
                      "Staged output assumes a coordinated single writer and refuses ordinary existing destinations.",
                  ]}
        args.output_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".windowed-source-resolution-", dir=args.output_dir.parent) as stage:
            output = Path(stage)/"result"
            output.mkdir()
            (output/"input_config.json").write_bytes(raw)
            for case in config["cases"]:
                name = case["case_id"]
                _, source_path, evidence_path = loaded[name]
                shutil.copyfile(source_path, output/f"source_{name}.h5")
                shutil.copyfile(evidence_path, output/f"evidence_{name}.json")
                (output/f"case_{name}.json").write_text(json.dumps(estimates[name], indent=2, allow_nan=False)+"\n")
                if (sha(source_path) != case["source_sha256"] or sha(output/f"source_{name}.h5") != case["source_sha256"]
                        or sha(evidence_path) != case["evidence_sha256"] or sha(output/f"evidence_{name}.json") != case["evidence_sha256"]):
                    raise ValueError("case input changed during inspection")
            if sha(args.config) != report["input_config_sha256"] or any(sha(ROOT/name) != value for name, value in code_hashes.items()):
                raise ValueError("configuration or implementation changed during inspection")
            for case in config["cases"]:
                _, source_path, evidence_path = loaded[case["case_id"]]
                if sha(source_path) != case["source_sha256"] or sha(evidence_path) != case["evidence_sha256"]:
                    raise ValueError("case input changed during inspection")
            report["output_sha256"] = {p.name: sha(p) for p in output.iterdir() if p.is_file()}
            (output/"report.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
            if args.output_dir.exists() or args.output_dir.is_symlink():
                raise ValueError("output path appeared during inspection; refusing replacement")
            os.rename(output, args.output_dir)
    except (OSError, ValueError, TypeError, KeyError, OverflowError, RuntimeError) as error:
        parser.error(str(error))
    print(json.dumps({"output_dir": str(args.output_dir), "cases": len(loaded), "source_comparisons": len(comparisons),
                      "parameter_status": config["parameter_status"], "source_evolution_repeated": False}))


if __name__ == "__main__":
    main()
