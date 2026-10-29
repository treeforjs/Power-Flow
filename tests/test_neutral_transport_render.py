"""Saved-map orientation, exact-zero and encoded-GIF contracts."""
import json

import numpy as np
from PIL import Image
import pytest

from mhdlab.neutral_transport_render import load_transport_maps, render_transport


def fixture(tmp_path):
    density = np.zeros((3, 4, 2, 3))
    density[0, 0, 0, 0] = 1e15
    density[1, 0, 1, 1] = 2e15
    density[2, 0, 1, 2] = 3e15
    density[:, 2, 0, 2] = 5e14
    density[:, 3, 1, 0] = 8e14
    path = tmp_path / "maps.npz"
    np.savez(path, time_s=[0., 1e-8, 2e-8], x_edges_m=[-.003, -.001, 0., .003],
             y_edges_m=[-.002, 0., .002], species=np.array(["H", "H2", "H2O", "OH"]),
             number_density_m3=density, foil_bounds_m=[-.001, .001, -.0001, .0001])
    metadata = tmp_path / "metadata.json"
    metadata.write_text(json.dumps({"status": "manufactured", "time_registration": "authored clock",
        "invariant_current_length_m": .01, "model_description": "orientation and encoding fixture",
        "assumptions": ["Authored densities, no transport or physical reaction model"]}))
    return path, metadata, density


def test_encoding_retains_snapshots_orientation_and_exact_zero_species(tmp_path):
    path, metadata, original = fixture(tmp_path)
    output = tmp_path / "render"
    report = render_transport(path, metadata, output, fps=10, dpi=60)
    assert report["frame_count"] == 3
    assert report["frame_duration_ms"] == [100, 100, 100]
    assert report["display_density_units"] == "cm^-3"
    assert report["display_color_scales"][0]["positive_min_cm3"] == 1e9
    assert report["display_color_scales"][0]["positive_max_cm3"] == 3e9
    assert report["display_color_scales"][1]["log10_color_limits"] is None
    assert report["zeros_replaced_by_floor"] is False
    # H occupies different unequal-width cells. This tests the y,x convention.
    np.testing.assert_allclose(np.array(report["integrated_plot_particle_numbers"])[:, 0],
                               [1e15*.002*.002*.01, 2e15*.001*.002*.01, 3e15*.003*.002*.01])
    arrays, _ = load_transport_maps(path, metadata)
    np.testing.assert_array_equal(arrays["number_density_m3"], original)
    with Image.open(output / "neutral_transport.gif") as image:
        assert image.n_frames == 3
        assert image.info["loop"] == 0
    assert (output / "preview_peak.png").is_file()
    assert report["physical_model_qualified_by_renderer"] is False


def test_overwrite_rejects_without_changing_prior_output(tmp_path):
    path, metadata, _ = fixture(tmp_path)
    output = tmp_path / "occupied"
    output.mkdir()
    marker = output / "existing.txt"
    marker.write_text("preserve")
    with pytest.raises(ValueError, match="new render output"):
        render_transport(path, metadata, output, dpi=60)
    assert marker.read_text() == "preserve"


def test_declared_normal_magnification_changes_display_only(tmp_path):
    path, metadata, original = fixture(tmp_path)
    record = json.loads(metadata.read_text())
    record["normal_axis_magnification"] = 4.0
    metadata.write_text(json.dumps(record))
    report = render_transport(path, metadata, tmp_path / "magnified", dpi=60)
    assert report["geometry_aspect"] == "declared_normal_axis_magnification"
    assert report["normal_axis_magnification"] == 4
    assert report["density_values_changed"] is False
    arrays, _ = load_transport_maps(path, metadata)
    np.testing.assert_array_equal(arrays["number_density_m3"], original)
    record["normal_axis_magnification"] = 0
    metadata.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="magnification"):
        render_transport(path, metadata, tmp_path / "invalid", dpi=60)


def test_color_span_retains_small_positive_values_and_particle_totals(tmp_path):
    path, metadata, original = fixture(tmp_path)
    with np.load(path) as source:
        arrays = {name: source[name] for name in source.files}
    arrays["number_density_m3"][0, 0, 0, 1] = 1e-200
    np.savez(path, **arrays)
    report = render_transport(path, metadata, tmp_path / "span", dpi=60, color_span_decades=8)
    scale = report["display_color_scales"][0]
    assert scale["positive_min_cm3"] == 1e-206
    assert scale["log10_color_limits"] == [np.log10(3e15)-8-6, np.log10(3e15)-6]
    assert scale["positive_bins_below_color_range"] == 1
    assert report["density_values_changed"] is False
    assert report["zeros_replaced_by_floor"] is False
    restored, _ = load_transport_maps(path, metadata)
    np.testing.assert_array_equal(restored["number_density_m3"], arrays["number_density_m3"])
    assert report["integrated_plot_particle_numbers"][0][0] == 1e15*.002*.002*.01


def test_viewport_crops_axes_only_and_keeps_global_colors_times_and_totals(tmp_path):
    path, metadata, original = fixture(tmp_path)
    full = render_transport(path, metadata, tmp_path / "full", dpi=60, color_span_decades=8)
    bounds = [0., .002, -.001, .001]
    crop = render_transport(path, metadata, tmp_path / "crop", dpi=60,
                            color_span_decades=8, view_bounds_m=bounds)
    assert crop["display_view_bounds_m"] == bounds
    assert crop["full_map_bounds_m"] == [-.003, .003, -.002, .002]
    assert crop["view_crops_axes_only"] is True
    assert crop["geometry_aspect"] == "equal_millimetre_units"
    for key in ("input_sha256", "physical_time_s", "display_color_scales", "integrated_plot_particle_numbers"):
        assert crop[key] == full[key]
    assert "whole saved map" in crop["integrated_particle_number_scope"]
    arrays, _ = load_transport_maps(path, metadata)
    np.testing.assert_array_equal(arrays["number_density_m3"], original)


@pytest.mark.parametrize("bounds", [
    [-.004, .001, -.001, .001], [0., .001, -.001, .003],
    [.001, 0., -.001, .001], [0., .001, float("nan"), .001], [0., .001, 0.],
])
def test_viewport_requires_a_finite_ordered_rectangle_inside_saved_domain(tmp_path, bounds):
    path, metadata, _ = fixture(tmp_path)
    with pytest.raises(ValueError, match="view_bounds_m"):
        render_transport(path, metadata, tmp_path / "invalid_view", dpi=60, view_bounds_m=bounds)


@pytest.mark.parametrize("mutation", [
    lambda a: a.update(number_density_m3=-a["number_density_m3"]),
    lambda a: a.update(time_s=np.array([0., 0., 2e-8])),
    lambda a: a.update(number_density_m3=np.transpose(a["number_density_m3"], (0, 1, 3, 2))),
])
def test_invalid_physical_values_or_axis_contract_reject(tmp_path, mutation):
    path, metadata, _ = fixture(tmp_path)
    with np.load(path) as source:
        arrays = {name: source[name] for name in source.files}
    mutation(arrays)
    np.savez(path, **arrays)
    with pytest.raises(ValueError):
        load_transport_maps(path, metadata)
