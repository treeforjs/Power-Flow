#!/usr/bin/env python3
"""Validate dissertation evidence identity and report physical-comparison readiness."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data/dissertation/pilot_comparison_manifest.json"


def validate_manifest(manifest: dict, root: Path = ROOT) -> dict:
    if manifest.get("schema") != "mhdlab-dissertation-pilot-comparison-v1" or manifest.get("revision") not in {1, 2}:
        raise ValueError("Unsupported pilot comparison schema/revision")
    if manifest["revision"] == 2:
        backend = manifest.get("heating_backend", {})
        if backend.get("implementation") != "mfem_joule" or backend.get("entrypoint") != "miniapps/electromagnetics/joule.cpp":
            raise ValueError("Pilot revision 2 requires the selected MFEM Joule heating backend")
        if backend.get("flash_acquisition_required") is not False:
            raise ValueError("FLASH acquisition is no longer a comparison requirement")
    scan = manifest["calibration"]["scale_scan"]
    expected = [10 ** (-2 + i / 10) for i in range(41)]
    values = scan["values"]
    if len(values) != 41 or any(not math.isfinite(v) or not math.isclose(v, e, rel_tol=1e-13) for v, e in zip(values, expected)):
        raise ValueError("C_abs must have 41 logarithmic sensitivity points from 0.01 to 100")
    if scan["definition"] != "scenario_corrected_radiance = C_abs * nominal_radiance":
        raise ValueError("C_abs direction must multiply nominal radiance")
    if scan["range_role"] != "exploratory_sensitivity_not_confidence" or scan["apply_to_published_derived_density"] is not False:
        raise ValueError("Exploratory radiance scale cannot become confidence or scale derived density")
    if manifest["calibration"]["lamp_etendue_is_shot_etendue"] is not False:
        raise ValueError("Lamp/slit etendue cannot be assigned to shot relay")
    for family in ("unit_scenarios", "reduction_scenarios", "wavelength_response_scenarios"):
        scenarios = manifest["calibration"][family]
        ids = [item["id"] for item in scenarios]
        if len(set(ids)) != len(ids):
            raise ValueError("Reconstruction scenario IDs must be unique within a family")
        for scenario in scenarios:
            if scenario["status"] in {"unverified_reconstruction_assumption", "requires_sourced_data"} and (
                scenario.get("radiance_admissible") or scenario.get("absolute_inference_admissible")
            ):
                raise ValueError("Unverified reconstruction scenarios cannot admit absolute inference")
    policy = manifest["comparison_policy"]
    if any(policy[key] is not False for key in ("statistical_likelihood", "glyph_bounds_are_noise", "dense_digitization_creates_statistical_weight")):
        raise ValueError("Correlated glyph readout envelopes cannot become a statistical likelihood")
    if policy["heating_benchmark_is_pilot_drive"] is not False:
        raise ValueError("Shot 17133 heating benchmark cannot substitute for pilot shot drive")

    root = root.resolve()
    referenced = {}
    for artifact in manifest["artifacts"]:
        name = artifact["path"]
        path = (root / name).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Artifact path escapes repository")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != artifact["sha256"]:
            raise ValueError(f"Artifact SHA256 mismatch: {name}")
        if name in referenced:
            raise ValueError(f"Duplicate artifact reference: {name}")
        referenced[name] = path

    observations = manifest["observations"]
    if len({item["id"] for item in observations}) != len(observations):
        raise ValueError("Observation IDs must be unique")
    for item in observations:
        if item["quantity"] == "plotted_transition_density":
            if item.get("scale_role") != "published_derived_unscaled" or item.get("allowed_C_abs") != [1]:
                raise ValueError("Literal transition-density displays must remain unscaled")
        if item.get("table_time_is_gate_endpoint"):
            raise ValueError("Table before-breakdown time is not a registered gate endpoint")
        if item.get("time_registration") == "unavailable" and (
            item.get("gate_start_s") is not None or item.get("gate_end_s") is not None
        ):
            raise ValueError("Unregistered photon gate cannot receive invented endpoints")
        if item.get("absolute_vuv_constraint") is True and item["quantity"] == "arbitrary_intensity":
            raise ValueError("Visible arbitrary intensity cannot supply absolute VUV amplitude")
        if "artifact" not in item:
            continue
        if item["artifact"] not in referenced or item["sidecar"] not in referenced:
            raise ValueError("Observation artifacts need versioned hash references")
        with referenced[item["artifact"]].open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != item["rows"]:
            raise ValueError(f"Observation row count mismatch: {item['id']}")
        if "center_column" in item:
            for row in rows:
                center, low, high = [float(row[item[key]]) for key in ("center_column", "lower_column", "upper_column")]
                if not all(math.isfinite(v) for v in (center, low, high)) or not low <= center <= high:
                    raise ValueError("Readout envelope must contain its finite center")
    readiness = manifest["readiness"]
    coverage_only = readiness.get("coverage_only_blocking_evidence", [])
    if not isinstance(coverage_only, list) or set(coverage_only) - {"sourced monolayer convention for coverage"}:
        raise ValueError("Only the monolayer convention can be a coverage-only blocker")
    if readiness["physical_inference_ready"] and readiness["blocking_evidence"]:
        raise ValueError("Physical inference cannot be ready with unresolved blocking evidence")
    def observation_readiness(item: dict) -> dict:
        mapping = item["mapping_status"] == "qualified"
        clock = item.get("time_registration") not in {None, "unavailable"} and all(
            isinstance(item.get(key), (int, float)) and math.isfinite(item[key])
            for key in ("gate_start_s", "gate_end_s")
        )
        collection = item.get("collection_kernel_status") == "qualified"
        return {"id": item["id"], "shot": item["shot"], "status": item["status"],
                "quantity": item["quantity"], "mapping_status": item["mapping_status"],
                "display_support_ready": "artifact" in item,
                "reduction_mapping_ready": mapping, "signed_photon_clock_ready": clock,
                "collection_kernel_ready": collection,
                "physical_observation_admissible": mapping and clock and collection,
                "blocking_evidence": [label for ready, label in (
                    (mapping, "qualified physical-to-published reduction mapping"),
                    (clock, "registered signed photon gate"),
                    (collection, "qualified collection kernel and simulated ray geometry"),
                ) if not ready]}

    return {
        "schema": "mhdlab-dissertation-pilot-readiness-v1", "validation": "passed",
        "manifest_revision": manifest["revision"], "heating_backend": manifest.get("heating_backend"),
        "validated_artifact_count": len(referenced),
        "observations": [observation_readiness(item) for item in observations],
        "calibration_sensitivity_points": len(values), "readiness": readiness,
        "physical_parameters_estimated": False, "absolute_correction_applied": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--require-physical-ready", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    report = validate_manifest(manifest)
    report["manifest_sha256"] = hashlib.sha256(args.manifest.read_bytes()).hexdigest()
    rendered = json.dumps(report, indent=2) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    if args.require_physical_ready and not report["readiness"]["physical_inference_ready"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
