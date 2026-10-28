"""Scientific admission guards and reproducible published VUV readouts."""
import copy
import csv
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def load_tool(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


recovery = load_tool("recover_dissertation_vuv_lineouts")
validator = load_tool("validate_pilot_comparison")


@pytest.fixture
def manifest():
    return json.loads((ROOT / "data/dissertation/pilot_comparison_manifest.json").read_text(encoding="utf-8"))


def test_evidence_manifest_validates_hashes_and_withholds_physical_claims(manifest):
    report = validator.validate_manifest(manifest)
    assert report["validation"] == "passed"
    assert report["validated_artifact_count"] == 14
    assert not report["physical_parameters_estimated"]
    assert not report["readiness"]["physical_inference_ready"]
    assert report["calibration_sensitivity_points"] == 41
    assert all(not o["physical_observation_admissible"] for o in report["observations"])
    assert all(not o["signed_photon_clock_ready"] for o in report["observations"])
    assert all(not o["collection_kernel_ready"] for o in report["observations"])


def test_mfem_backend_selection_removes_flash_acquisition_without_admitting_physics(manifest):
    report = validator.validate_manifest(manifest)
    assert report["manifest_revision"] == 2
    assert report["heating_backend"]["implementation"] == "mfem_joule"
    assert report["heating_backend"]["flash_acquisition_required"] is False
    assert report["heating_backend"]["physical_case_qualified"] is False
    assert "qualified MFEM Joule heating and conservative source provenance" in report["readiness"]["blocking_evidence"]
    changed = copy.deepcopy(manifest)
    changed["heating_backend"]["flash_acquisition_required"] = True
    with pytest.raises(ValueError, match="no longer"):
        validator.validate_manifest(changed)


@pytest.mark.parametrize("mutation,match", [
    (lambda m: m["artifacts"][0].update(sha256="0" * 64), "SHA256"),
    (lambda m: m["artifacts"][0].update(path="../outside.csv"), "escapes"),
    (lambda m: m["calibration"]["scale_scan"]["values"].pop(), "41 logarithmic"),
    (lambda m: m["calibration"]["scale_scan"].update(apply_to_published_derived_density=True), "scale derived"),
    (lambda m: m["calibration"]["scale_scan"].update(definition="nominal_radiance = C_abs * scenario_corrected_radiance"), "direction"),
    (lambda m: m["calibration"].update(lamp_etendue_is_shot_etendue=True), "Lamp/slit"),
    (lambda m: m["calibration"]["unit_scenarios"][1].update(radiance_admissible=True), "Unverified"),
    (lambda m: m["comparison_policy"].update(glyph_bounds_are_noise=True), "glyph"),
    (lambda m: m["comparison_policy"].update(heating_benchmark_is_pilot_drive=True), "17133"),
    (lambda m: m["observations"][0].update(allowed_C_abs=[0.1, 1, 10]), "unscaled"),
    (lambda m: m["observations"][0].update(table_time_is_gate_endpoint=True), "endpoint"),
    (lambda m: m["observations"][0].update(gate_start_s=-30e-9), "invented"),
    (lambda m: m["observations"][-1].update(absolute_vuv_constraint=True), "Visible"),
    (lambda m: m["readiness"].update(physical_inference_ready=True), "blocking"),
    (lambda m: m["readiness"].update(coverage_only_blocking_evidence=["Raikou physical boundary source acceptance"]), "coverage-only"),
])
def test_manifest_rejects_unsafe_scientific_promotions(manifest, mutation, match):
    changed = copy.deepcopy(manifest)
    mutation(changed)
    with pytest.raises(ValueError, match=match):
        validator.validate_manifest(changed)


def test_literal_vuv_artifact_has_partial_segmented_support_and_no_invented_noise():
    metadata = json.loads((ROOT / "data/dissertation/figure_4_20_wavelength_lineouts.json").read_text())
    with (ROOT / "data/dissertation/figure_4_20_wavelength_lineouts.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 157
    assert metadata["row_count_by_series"] == {"top": 40, "middle": 38, "bottom": 79}
    assert metadata["observable"]["quantity"] == "plotted_transition_density"
    assert not metadata["observable"]["absolute_radiance"]
    assert not metadata["observable"]["integrals_supported"]
    for series in recovery.SERIES:
        selected = [r for r in rows if r["series"] == series]
        for left, right in zip(selected, selected[1:]):
            same_segment = left["support_segment_id"] == right["support_segment_id"]
            assert same_segment == (int(right["pixel_x"]) - int(left["pixel_x"]) == 4)


def test_native_pixel_selection_excludes_legend_and_preserves_missing_runs():
    image = np.full((692, 1040, 3), 255, dtype=np.uint8)
    rgb = recovery.SERIES["top"]
    image[600:603, 514] = rgb
    image[500:503, 518] = rgb  # within legend: omitted
    image[600:602, 522] = rgb
    image[600:602, 526] = rgb
    image[590:592, 530] = rgb
    image[600:602, 530] = rgb  # separated runs: omitted
    image[600:612, 534] = rgb  # broad stroke: omitted
    image[600:602, 538] = (127, 128, 0)  # antialias color: omitted
    rows, details = recovery.recover(image)
    assert [r["pixel_x"] for r in rows] == [514, 522, 526]
    assert [r["support_segment_id"] for r in rows] == [0, 1, 1]
    assert rows[0]["transition_density_lower_cm3"] < rows[0]["nominal_transition_density_cm3"] < rows[0]["transition_density_upper_cm3"]
    omitted = {(r["series"], r["pixel_x"]): r["reason"] for r in details["omitted_columns"]}
    assert omitted["top", 518] == "no_unobscured_exact_color"
    assert omitted["top", 530] == "separated_runs"
    assert omitted["top", 534] == "overly_broad_stroke"


def test_wrong_pdf_fails_identity_guard_before_decoder(tmp_path):
    source = tmp_path / "different.pdf"
    source.write_bytes(b"different source")
    with pytest.raises(ValueError, match="hash mismatch"):
        recovery.load_image(source)


@pytest.mark.parametrize("figure", ["4.20", "4.21"])
def test_exact_pdf_reproduces_committed_lineout_values(figure):
    source = os.environ.get("DISSERTATION_PDF")
    if not source:
        pytest.skip("Set DISSERTATION_PDF to external exact source for integration check")
    rows, details = recovery.recover(recovery.load_image(Path(source), figure), figure)
    filename = "figure_" + figure.replace(".", "_") + "_wavelength_lineouts.csv"
    with (ROOT / "data/dissertation" / filename).open(newline="") as handle:
        committed = list(csv.DictReader(handle))
    assert len(rows) == len(committed)
    for actual, expected in zip(rows, committed):
        for key, value in actual.items():
            if isinstance(value, str):
                assert value == expected[key]
            else:
                assert value == float(expected[key])
    assert details["omitted_columns"]
