"""Native-marker identification, skew calibration, units and source guards."""

import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/recover_wilkes_304_resistivity.py"
CSV = ROOT / "data/materials/wilkes_1997_304_resistivity.csv"
METADATA = ROOT / "data/materials/wilkes_1997_304_resistivity.json"
spec = importlib.util.spec_from_file_location("recover_wilkes_304_resistivity", TOOL)
recovery = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recovery)


def diamond(open_face):
    y, x = np.indices((33, 33))
    distance = abs(x-16) + abs(y-16)
    image = np.full((33, 33), 255, dtype=np.uint8)
    image[distance <= 10] = 0
    if open_face:
        image[distance <= 5] = 255
    return image


def test_open_and_filled_markers_require_the_actual_legend_shape():
    x, y, info = recovery.marker_centre(diamond(True), "cooling", 16, 16)
    assert (x, y) == (16, 16)
    assert info["white_hole_pixels"] == 61
    assert info["centre_method"] == "centroid_of_enclosed_white_face"
    x, y, info = recovery.marker_centre(diamond(False), "heating", 16, 16)
    assert (x, y) == (16, 16)
    assert info["black_taxicab_depth_px"] == 11
    assert info["max_depth_plateau_pixels"] == 1
    with pytest.raises(ValueError):
        recovery.marker_centre(diamond(True), "heating", 16, 16)
    with pytest.raises(ValueError):
        recovery.marker_centre(diamond(False), "cooling", 16, 16)


def test_reference_line_and_unenclosed_face_are_not_measurement_markers():
    image = np.full((33, 33), 255, dtype=np.uint8)
    image[15:18] = 0
    with pytest.raises(ValueError, match="filled diamond"):
        recovery.marker_centre(image, "heating", 16, 16)
    image = diamond(True)
    image[:17, 16] = 255  # Break the enclosing border into the outside background.
    with pytest.raises(ValueError, match="enclosed white diamond"):
        recovery.marker_centre(image, "cooling", 16, 16)
    with pytest.raises(ValueError, match="branch"):
        recovery.marker_centre(diamond(False), "unknown", 16, 16)


def test_affine_mapping_recovers_known_skew_and_has_no_axis_unit_swap():
    xy = np.array([[10, 20], [110, 20], [10, 220], [210, 120]], dtype=float)
    expected = np.array([[.5, .01], [-.02, -.1], [10, 100]])
    physical = np.column_stack((xy, np.ones(4))) @ expected
    coefficients, residuals = recovery.affine_fit(np.column_stack((xy, physical)))
    assert coefficients == pytest.approx(expected, abs=1e-12)
    assert residuals == pytest.approx(np.zeros((4, 2)), abs=1e-12)
    assert np.array([60, 70, 1]) @ coefficients == pytest.approx([38.6, 93.6])


@pytest.mark.parametrize("anchors", [
    np.zeros((2, 4)), np.zeros((4, 3)), np.array([[0, 0, 0, 0], [1, 1, 1, 1], [2, 2, 2, 2]]),
    np.array([[0, 0, 0, 0], [1, 0, 1, 0], [0, 1, np.nan, 1]]),
])
def test_invalid_or_rank_deficient_axis_anchors_are_rejected(anchors):
    with pytest.raises(ValueError):
        recovery.affine_fit(anchors)


def test_published_subset_units_source_identity_and_nonphysical_envelopes():
    metadata = json.loads(METADATA.read_text())
    rows = list(csv.DictReader(CSV.open(newline="")))
    assert len(rows) == metadata["output"]["rows"] == 10
    assert metadata["output"]["csv_sha256"] == hashlib.sha256(CSV.read_bytes()).hexdigest()
    assert sum(r["branch"] == "heating" for r in rows) == 6
    assert sum(r["branch"] == "cooling" for r in rows) == 4
    assert {r["specimen_id"] for r in rows} == {"2"}
    assert {r["nominal_foil_thickness_mm"] for r in rows} == {"0.203"}
    assert {r["material_grade"] for r in rows} == {"AISI Type 304"}
    assert rows[0]["point_id"] == "C1" and rows[-1]["point_id"] == "H6"
    assert float(rows[0]["temperature_source_degC"]) == pytest.approx(25.67, abs=.01)
    assert float(rows[-1]["temperature_source_degC"]) == pytest.approx(505.01, abs=.01)
    assert float(rows[0]["resistivity_source_microohm_cm"]) == pytest.approx(73.13, abs=.01)
    assert float(rows[-1]["resistivity_source_microohm_cm"]) == pytest.approx(104.70, abs=.01)
    for row in rows:
        t = float(row["temperature_source_degC"])
        rho = float(row["resistivity_ohm_m"])
        assert float(row["temperature_K"]) == pytest.approx(t + 273.15, abs=1e-8)
        assert rho == pytest.approx(float(row["resistivity_source_microohm_cm"])*1e-8, rel=1e-10)
        assert float(row["temperature_readout_lower_degC"]) < t < float(row["temperature_readout_upper_degC"])
        assert float(row["resistivity_readout_lower_ohm_m"]) < rho < float(row["resistivity_readout_upper_ohm_m"])
        for mask in recovery.EXCLUDED_OVERLAP_ROIS:
            x0,y0,x1,y1=mask["xyxy"]
            assert not (x0 <= float(row["pixel_x"]) <= x1 and y0 <= float(row["pixel_y"]) <= y1)
    assert metadata["source"]["pdf_sha256"] == "5738cdc4407bf6d09da8a8f2379105ff7a32c6d77c6fbb495d5391be67d154e7"
    assert metadata["observable"]["quantity"] == "electrical_resistivity"
    assert metadata["observable"]["SI_per_original_unit"] == 1e-8
    assert metadata["observable"]["coordinate_bounds_are_continuous_curve_support"] is False
    assert metadata["selection"]["original_marker_count_recovered"] is False
    assert metadata["selection"]["interpolation"] is False
    assert metadata["material"]["other_304_specimen"]["absolute_resistivity_recovery_status"] == "not_recovered"
    assert metadata["figure_2_status"]["heating_cooling_legend_available"] is False
    assert metadata["recovery"]["readout"]["physical_measurement_uncertainty"] is False
    assert metadata["source_error_context"]["confidence_level"] is None
    assert metadata["output"]["adopted_as_runtime_defaults"] is False


def test_changed_pdf_fails_before_extraction_and_creates_no_dataset(tmp_path):
    pdf=tmp_path/"changed.pdf";pdf.write_bytes(b"%PDF-1.4\nNot the inspected report")
    output=tmp_path/"table.csv";metadata=tmp_path/"table.json"
    result=subprocess.run([sys.executable,str(TOOL),"--pdf",str(pdf),"--csv",str(output),
        "--metadata",str(metadata)],capture_output=True,text=True,timeout=30)
    assert result.returncode == 2
    assert "PDF differs from the inspected source" in result.stderr
    assert not output.exists() and not metadata.exists()


def test_optional_external_pdf_recovery_matches_published_subset():
    pdf=Path(os.environ.get("WILKES_304_SOURCE_PDF",
        "/workspace/research_sources/ss304_properties/wilkes_1997_ornl_cp96813.pdf"))
    if not pdf.is_file():
        pytest.skip("Validated external source PDF is unavailable")
    pytest.importorskip("pymupdf",reason="Optional PDF extraction dependency is unavailable")
    pages,evidence=recovery.load_pages(pdf)
    rows,details=recovery.recover_markers(pages[6])
    assert recovery.csv_bytes(rows) == CSV.read_bytes()
    assert [e["pdf_page_1based"] for e in evidence] == [4,6,7]
    assert recovery.tick_anchors(pages[6]) == pytest.approx(recovery.INSPECTED_ANCHORS)
    assert details["calibration"]["maximum_absolute_residuals"] == pytest.approx([.464124828102,.092659479564])
    changed=pages[6].copy();changed[1850:1875,673:711]=255
    with pytest.raises(ValueError):
        recovery.tick_anchors(changed)
