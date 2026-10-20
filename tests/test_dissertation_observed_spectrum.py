"""Observed visible-lineout identity, native-pixel selection and source guards."""

import builtins
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/recover_dissertation_observed_spectrum.py"
CSV = ROOT / "data/dissertation/observed_spectrum_4_15.csv"
METADATA = CSV.with_suffix(".json")
SPEC = importlib.util.spec_from_file_location("observed_spectrum_recovery", TOOL)
recovery = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(recovery)


def test_axis_calibration_recovers_known_wavelength_and_downward_intensity_axes():
    wavelength, xres = recovery.axis_fit([10, 30, 70], [400, 420, 460])
    intensity, yres = recovery.axis_fit([5, 25, 45], [4, 2, 0])
    np.testing.assert_allclose(wavelength, [1, 390], rtol=0, atol=1e-13)
    np.testing.assert_allclose(intensity, [-.1, 4.5], rtol=0, atol=1e-14)
    np.testing.assert_allclose(xres, 0, rtol=0, atol=1e-13)
    np.testing.assert_allclose(yres, 0, rtol=0, atol=1e-14)
    assert np.array([30, 1]) @ wavelength == pytest.approx(420)
    assert np.array([25, 1]) @ intensity == pytest.approx(2)


@pytest.mark.parametrize("pixels,values", [
    ([1], [400]), ([1, 1], [400, 500]), ([1, 2], [400]),
    ([1, np.nan], [400, 500]), ([1, 2], [400, np.inf]),
])
def test_axis_calibration_rejects_undefined_fits(pixels, values):
    with pytest.raises(ValueError):
        recovery.axis_fit(pixels, values)


def test_green_mask_requires_literal_series_rgb_and_inclusive_exclusion_bounds():
    image = np.full((12, 14, 3), 255, dtype=np.uint8)
    image[2:6, 3:8] = (60, 179, 113)
    image[8, 1:7] = [(61, 179, 113), (60, 178, 113), (60, 179, 114),
                      (0, 0, 255), (255, 0, 0), (100, 100, 100)]
    expected = np.zeros((12, 14), dtype=bool)
    expected[2:6, 3:8] = True
    assert recovery.GREEN_RGB == (60, 179, 113)
    np.testing.assert_array_equal(recovery.green_mask(image), expected)
    expected[3:5, 4:7] = False
    masked = recovery.green_mask(image, exclusions=((4, 3, 6, 4),))
    assert masked.dtype == np.bool_
    np.testing.assert_array_equal(masked, expected)


def test_column_centre_uses_observed_run_median_without_peak_or_smoothing():
    mask = np.zeros((20, 20), dtype=bool)
    mask[7:12, 9] = True
    centre = recovery.column_center(mask, 9, (2, 17))
    assert centre["pixel_y"] == 9
    assert centre["matched_pixels"] == 5
    assert centre["y_first"] == 7 and centre["y_last"] == 11
    mask[12, 9] = True
    centre = recovery.column_center(mask, 9, (2, 17))
    assert centre["pixel_y"] == 9.5
    assert centre["matched_pixels"] == 6


@pytest.mark.parametrize("rows,min_pixels,max_pixels", [
    ([], 3, 25), ([8, 9], 3, 25), (list(range(5, 11)), 3, 5),
    ([5, 6, 7, 12, 13, 14], 3, 25), ([2, 3, 4], 3, 25),
    ([15, 16, 17], 3, 25),
])
def test_column_rejects_absent_small_broad_separated_and_boundary_runs(rows, min_pixels, max_pixels):
    mask = np.zeros((20, 20), dtype=bool)
    mask[rows, 9] = True
    with pytest.raises(ValueError):
        recovery.column_center(mask, 9, (2, 17), min_pixels=min_pixels, max_pixels=max_pixels)


def test_wrong_pdf_fails_before_optional_pdf_import_and_produces_no_artifact(tmp_path, monkeypatch):
    wrong = tmp_path / "wrong.pdf"
    wrong.write_bytes(b"not the inspected dissertation")
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name in {"pymupdf", "fitz"}:
            pytest.fail("wrong source reached the optional PDF decoder")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    with pytest.raises(ValueError, match="PDF"):
        recovery.load_source_image(wrong)
    assert list(tmp_path.iterdir()) == [wrong]


def test_changed_native_rgb_pixels_are_rejected_before_trace_recovery():
    image = np.full((1320, 1040, 3), 255, dtype=np.uint8)
    with pytest.raises(ValueError):
        recovery.recover_spectrum(image)


def published_artifact():
    with CSV.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return rows, json.loads(METADATA.read_text())


def test_published_artifact_keeps_observed_visible_arbitrary_intensity_and_separate_timing():
    rows, metadata = published_artifact()
    assert metadata["schema"] == "mhdlab-published-observed-visible-spectrum-v1"
    assert metadata["parameter_status"] == "reconstructed_published_observed_visible_trace"
    source = metadata["source"]
    assert source["pdf_sha256"] == "1d7b1f9c191a37e94f1cc3634587750f1af35b07eb13e0b94c60f3da47495ed7"
    assert source["native_rgb_sha256"] == "27a9beba2c5aa9edb0de720ff0f50866157eb845b12c1318864279f1b16d5276"
    assert (source["figure"], source["printed_page"], source["pdf_page_1based"], source["native_image_xref"]) == ("4.15", 75, 100, 2554)
    assert source["native_image_size_px"] == [1040, 1320]
    observable = metadata["observable"]
    assert observable["diagnostic"] == "SVS visible streak spectroscopy"
    assert observable["series_literal"] == "SVS Lineout 3: 330 ns"
    assert observable["axis_label_literal"] == "Intensity [arb.]"
    assert observable["intensity_unit"] == "arbitrary"
    assert observable["wavelength_unit_literal"] == "nm"
    for claim in ("raw_detector_counts", "absolute_radiance", "transition_density", "inferred_covariance"):
        assert observable[claim] is False
    context = metadata["experimental_context"]
    assert (context["shot"], context["foil_material_literal"], context["foil_thickness_um_caption"]) == (15765, "SS304", 400)
    assert context["visible_lineout_time_ns_literal"] == 330
    assert context["visible_exposure_or_averaging_width_ns"] is None
    assert context["source_emission_gate_inferred"] is False
    vuv = context["same_shot_vuv_context_only"]
    assert vuv == {"time_before_breakdown_ns_table_4_1": 29, "diagnostic_gate_ns": 3,
                   "is_this_visible_acquisition": False}
    assert context["collection_region"] == "SVS side-on gap view; distinct from VUV collection"
    geometry = context["collection_geometry_reported"]
    assert (geometry["source_printed_page"], geometry["source_pdf_page_1based"],
            geometry["setup_figure"], geometry["setup_figure_printed_page"]) == (71, 96, "4.13", 72)
    assert geometry["setup_figure_shot_range"] == [15750, 15802]
    assert geometry["focal_point_distance_from_cathode_um"] == 450
    assert geometry["focal_spot_radius_um"] == 175
    assert geometry["nearest_reported_collection_distance_um"] == 275
    assert geometry["primary_lens_focal_length_mm"] == 70.3
    assert geometry["object_distance_from_primary_optic_mm"] == 70.3
    assert geometry["reported_second_focusing_optic_distance_mm"] == 40.2
    assert geometry["reported_demagnification"] == .57
    assert geometry["spatial_kernel_or_etendue_inferred"] is False
    assert "correlated display samples" in metadata["limitations"][0]
    interpretation = metadata["recovery"]["readout_envelopes"]["interpretation"]
    assert "not physical error bars, independent noise, confidence intervals or covariance" in interpretation
    assert metadata["source_local_qa"]["contains_source_pixels"] is True
    assert metadata["source_local_qa"]["redistribution_in_delivery"] is False
    assert len(rows) == metadata["artifact"]["rows"] == 136
    assert list(rows[0]) == metadata["artifact"]["columns"] == list(recovery.FIELDS)
    assert hashlib.sha256(CSV.read_bytes()).hexdigest() == metadata["artifact"]["csv_sha256"] == "53b59e6d0645b718eaa677440097be51469b421d710bcc01327ae9de3e4825c5"


def test_samples_preserve_columns_missing_segments_axis_units_and_readout_envelopes():
    rows, metadata = published_artifact()
    sampling = metadata["recovery"]["sampling"]
    assert sampling["candidate_columns"] == 156
    assert sampling["retained_samples"] == 136
    assert sampling["interpolation"] is False and sampling["smoothing"] is False
    skipped = sampling["skipped_columns"]
    expected_skipped = [264, 268, 272, 276, 280, 300, 304, 308, 316, 320,
                        324, 328, 332, 360, 372, 376, 380, 384, 388, 392]
    assert [record["pixel_x"] for record in skipped] == expected_skipped
    assert all(record["reason"] == "broad_projected_green_run" and record["matched_green_pixels"] > 25 for record in skipped)
    x = np.array([int(row["pixel_x"]) for row in rows])
    wavelength = np.array([float(row["wavelength_nm"]) for row in rows])
    np.testing.assert_array_equal(x, sorted(set(range(220, 841, 4)) - set(expected_skipped)))
    assert np.all(np.diff(wavelength) > 0) and np.any(np.diff(x) > 4)
    np.testing.assert_allclose(wavelength, 350 + (x - 163.5) * 50 / 93, rtol=0, atol=2e-13)
    assert wavelength[0] == pytest.approx(380.3763440860215)
    assert wavelength[-1] == pytest.approx(713.7096774193549)
    for row in rows:
        y = float(row["pixel_y"])
        n = int(row["matched_green_pixels"])
        value = float(row["observed_intensity_arb"])
        assert 3 <= n <= 25 and 516 < y < 788
        # A contiguous pixel run has an integer/half-integer median and span n-1.
        assert 2 * y == int(2 * y)
        assert value == pytest.approx(480.50039179195846 - .6060531857583623 * y, abs=1e-12)
        assert float(row["wavelength_readout_lower_nm"]) < float(row["wavelength_nm"]) < float(row["wavelength_readout_upper_nm"])
        lo, hi = float(row["intensity_readout_lower_arb"]), float(row["intensity_readout_upper_arb"])
        assert lo < value < hi
        assert hi - lo == pytest.approx(.6060531857583623 * (n + 1) + 2 * .20201772858612077, abs=1e-12)
        assert not (626 <= int(row["pixel_x"]) <= 926 and 510 <= y <= 611)


def test_optional_external_pdf_reproduces_artifact_and_rejects_one_changed_native_pixel():
    source = Path(os.environ.get("DISSERTATION_SOURCE_PDF", "/workspace/research_sources/smtrevor_1.pdf"))
    if not source.is_file():
        pytest.skip("The exact dissertation PDF is an external cache artifact")
    try:
        import pymupdf  # noqa: F401
    except ImportError:
        pytest.importorskip("fitz", reason="optional PDF decoder unavailable")
    rows, metadata = recovery.recover(source)
    _, published = published_artifact()
    assert recovery.serialize_csv(rows) == CSV.read_bytes()
    assert metadata["source"] == published["source"]
    assert metadata["observable"] == published["observable"]
    assert metadata["experimental_context"] == published["experimental_context"]
    assert json.loads(json.dumps(metadata["recovery"])) == published["recovery"]
    image, details = recovery.load_source_image(source)
    assert details == published["source"]
    image[0, 0, 0] ^= np.uint8(1)
    with pytest.raises(ValueError, match="SHA256"):
        recovery.recover_spectrum(image)
