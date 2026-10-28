"""Export literal pilot readouts and guard their use in a physical ensemble.

No physical-to-published reduction is invented. Without qualified manifest
evidence, this exports an owner request packet and withholds comparison.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def local_tool(name):
    spec = importlib.util.spec_from_file_location("_pilot_" + name, ROOT / "tools" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


VALIDATOR = local_tool("validate_pilot_comparison")


def prepare_packet(manifest, observation_id, *, root=ROOT):
    readiness = VALIDATOR.validate_manifest(manifest, root)
    selected = [item for item in manifest["observations"] if item["id"] == observation_id]
    if len(selected) != 1:
        raise ValueError("select one existing observation ID")
    item = selected[0]
    if item.get("scale_role") != "published_derived_unscaled" or "center_column" not in item:
        raise ValueError("this bridge supports only literal bounded derived-display curves")
    path = root / item["artifact"]
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    series_column = "series" if "series" in rows[0] else "species"
    curves, indices = [], {}
    for series in item["series"]:
        selected_rows = [(index, row) for index, row in enumerate(rows) if row[series_column] == series]
        selected_rows.sort(key=lambda pair: float(pair[1][item["coordinate_column"]]))
        if not selected_rows:
            raise ValueError("declared observation series has no retained rows")
        curve_id = item["id"] + ":" + series
        indices[curve_id] = [index for index, _ in selected_rows]
        retained = [row for _, row in selected_rows]
        segment_column = item.get("support_segment_column")
        # With no segment declaration, every point stands alone: no continuity is inferred.
        segments = [int(row[segment_column]) for row in retained] if segment_column else list(range(len(retained)))
        curves.append({"curve_id": curve_id,
            "coordinate": [float(row[item["coordinate_column"]]) for row in retained],
            "center": [float(row[item["center_column"]]) for row in retained],
            "lower": [float(row[item["lower_column"]]) for row in retained],
            "upper": [float(row[item["upper_column"]]) for row in retained],
            "support_segment_ids": segments, "quantity": item["quantity"],
            "signal_units": item["signal_units"], "coordinate_units": item["coordinate_units"],
            "scale_role": item["scale_role"], "provenance": item["sidecar"]})
    selected_readiness = next(row for row in readiness["observations"] if row["id"] == observation_id)
    blockers = list(selected_readiness["blocking_evidence"])
    coverage_blockers = manifest["readiness"].get("coverage_only_blocking_evidence", [])
    blockers.extend(label for label in readiness["readiness"]["blocking_evidence"] if label not in coverage_blockers)
    return {"schema": "mhdlab-pilot-readout-packet-v1", "observation_id": observation_id,
            "shot": item["shot"], "curves": curves, "source_row_indices_zero_based": indices,
            "literal_observation": item, "allowed_C_abs": [1],
            "readiness": selected_readiness, "comparison_blockers": list(dict.fromkeys(blockers)),
            "coverage_only_blockers": coverage_blockers,
            "physical_comparison_ready": not blockers,
            "source_artifact_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "interpretation": "Literal figure readouts for owner prediction requests. No radiance conversion, absolute correction, clock registration or neutral inference is supplied."}


def bind_literal_observation(packet, observation):
    """Reject substituted units, ordinates, support, readouts or unresolved gates."""
    expected = {curve["curve_id"]: curve for curve in packet["curves"]}
    supplied = {curve.curve_id: curve for curve in observation.curves}
    if supplied.keys() != expected.keys():
        raise ValueError("ensemble must use exactly the selected pilot curve IDs")
    for key, curve in supplied.items():
        for field in ("quantity", "signal_units", "coordinate_units", "scale_role"):
            if getattr(curve, field) != expected[key][field]:
                raise ValueError(f"pilot observation has substituted {field}")
        for field in ("coordinate", "center", "lower", "upper", "support_segment_ids"):
            if getattr(curve, field).tolist() != expected[key][field]:
                raise ValueError(f"pilot observation has substituted {field}")
    if not packet["physical_comparison_ready"]:
        raise ValueError("pilot comparison withheld: " + "; ".join(packet["comparison_blockers"]))
    item = packet["literal_observation"]
    if observation.status != "reconstructed" or observation.mapping_status != item["mapping_status"]:
        raise ValueError("ensemble status or mapping differs from the qualified pilot manifest")
    for field in ("gate_start_s", "gate_end_s", "time_registration"):
        if getattr(observation, field) != item[field]:
            raise ValueError(f"ensemble {field} differs from registered pilot evidence")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=ROOT / "data/dissertation/pilot_comparison_manifest.json")
    parser.add_argument("--observation-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ensemble-config", type=Path, help="compare only after pilot physical readiness and literal binding")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("use a fresh output file")
    raw = args.manifest.read_bytes()
    packet = prepare_packet(json.loads(raw), args.observation_id)
    packet["manifest_sha256"] = hashlib.sha256(raw).hexdigest()
    packet["bridge_source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    failed = None
    if args.ensemble_config:
        try:
            engine = local_tool("compare_neutral_ensemble")
            config, observation, cases = engine.load_config(args.ensemble_config)
            bind_literal_observation(packet, observation)
            packet["comparison"] = engine.INFERENCE.compare_neutral_ensemble(observation, cases, config.get("calibration_scales"))
            packet["ensemble_config_sha256"] = hashlib.sha256(args.ensemble_config.read_bytes()).hexdigest()
        except (ValueError, TypeError, KeyError) as exc:
            failed = str(exc)
            packet["comparison_withheld_reason"] = failed
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(packet, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"packet": str(args.output), "physical_comparison_ready": packet["physical_comparison_ready"],
                      "comparison_blockers": packet["comparison_blockers"]}))
    if failed:
        parser.exit(2, failed + "\n")


if __name__ == "__main__":
    main()
