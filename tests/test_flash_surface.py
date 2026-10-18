"""Interchange and analytical width-integration checks, using synthetic fixtures."""

import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pytest

from mhdlab.flash_surface import (
    COORDINATE_SYSTEM,
    SCHEMA,
    FlashSurfaceHistory,
    load_surface_history,
    thermal_width_report,
    write_surface_history,
)


def history(**overrides):
    values = dict(
        time_s=[0.0, 2.0],
        x_edges_m=[0.0, 0.05, 0.3, 0.55, 0.9, 1.0],
        temperature_k=[[100.0, 200.0, 300.0, 400.0, 500.0]] * 2,
        source_kind="synthetic",
        run_id="unit-test-only",
        source_reference="analytical synthetic fixture; no experimental or FLASH data",
        surface_region="cathode gap-facing surface",
    )
    values.update(overrides)
    return FlashSurfaceHistory(**values)


def test_round_trip_retains_si_values_and_provenance(tmp_path):
    original = history(source_kind="flash", flash_version="4.8", source_reference="unit-test metadata only; not an actual FLASH export")
    target = tmp_path / "surface.h5"
    write_surface_history(target, original)
    loaded = load_surface_history(target)
    for name in ("time_s", "x_edges_m", "temperature_k"):
        np.testing.assert_array_equal(getattr(loaded, name), getattr(original, name))
    assert loaded.attributes == original.attributes
    with h5py.File(target) as handle:
        assert handle.attrs["schema"] == SCHEMA
        assert handle.attrs["coordinate_system"] == COORDINATE_SYSTEM
        assert handle["time_s"].attrs["units"] == "s"
        assert handle["x_edges_m"].attrs["units"] == "m"
        assert handle["temperature_k"].attrs["units"] == "K"


def test_uniform_null_has_no_edge_contrast_despite_unequal_bins():
    sample = history(temperature_k=np.full((2, 5), 700.0))
    metrics = sample.thermal_contrast()
    for name in ("left_edge_mean_temperature_k", "right_edge_mean_temperature_k", "center_mean_temperature_k"):
        np.testing.assert_allclose(metrics[name], 700.0)
    np.testing.assert_allclose(metrics["edge_minus_center_k"], 0.0)
    np.testing.assert_allclose(metrics["edge_to_center_ratio"], 1.0)


def test_unequal_bins_use_exact_interval_overlap():
    # Left 0-.1 overlaps 100 K by .05 and 200 K by .05 => 150 K.
    # Center .4-.6 overlaps 300 K by .15 and 400 K by .05 => 325 K.
    # Right .9-1 is 500 K. Their combined edge mean is 325 K.
    metrics = history().thermal_contrast()
    np.testing.assert_allclose(metrics["left_edge_mean_temperature_k"], 150.0)
    np.testing.assert_allclose(metrics["right_edge_mean_temperature_k"], 500.0)
    np.testing.assert_allclose(metrics["center_mean_temperature_k"], 325.0)
    np.testing.assert_allclose(metrics["edge_minus_center_k"], 0.0, atol=1e-12)


def test_single_spatial_bin_and_touching_regions_are_supported():
    sample = history(x_edges_m=[-0.5, 0.5], temperature_k=[[300.0], [600.0]])
    metrics = sample.thermal_contrast(edge_fraction=0.25, center_fraction=0.5)
    np.testing.assert_allclose(metrics["center_mean_temperature_k"], [300.0, 600.0])
    np.testing.assert_allclose(metrics["edge_to_center_ratio"], 1.0)


def test_numerically_unresolvable_region_raises_value_error():
    with pytest.raises(ValueError, match="too small to resolve"):
        history().thermal_contrast(edge_fraction=1e-300)


def test_interpolation_is_piecewise_linear_and_endpoints_are_included():
    sample = history(time_s=[0.0, 1.0, 3.0], temperature_k=np.array([[100.0] * 5, [300.0] * 5, [700.0] * 5]))
    np.testing.assert_allclose(sample.temperature_at(2.0), 500.0)
    np.testing.assert_allclose(sample.temperature_at([0.0, 0.5, 3.0])[:, 0], [100.0, 200.0, 700.0])
    report = thermal_width_report(sample, time_s=2.0)
    assert report["source_kind"] == "synthetic"
    assert report["metrics"]["center_mean_temperature_k"] == pytest.approx(500.0)
    assert "temperature contrast only" in report["interpretation"]


@pytest.mark.parametrize("query", [-1.0, 2.01, np.nan, np.inf, [], [[0.0]], "invalid"])
def test_interpolation_rejects_extrapolation_and_invalid_queries(query):
    with pytest.raises(ValueError):
        history().temperature_at(query)


@pytest.mark.parametrize("overrides", [
    {"time_s": [0.0]}, {"time_s": [0.0, 0.0]}, {"time_s": [1.0, 0.0]},
    {"time_s": [-1.0, 1.0]}, {"time_s": [0.0, np.inf]}, {"time_s": [[0.0, 1.0]]},
    {"x_edges_m": [0.0]}, {"x_edges_m": [0.0, 0.0]}, {"x_edges_m": [1.0, 0.0]},
    {"x_edges_m": [0.0, np.nan]}, {"temperature_k": np.ones((2, 4))},
    {"temperature_k": np.ones(5)}, {"temperature_k": np.zeros((2, 5))},
    {"temperature_k": np.full((2, 5), -1.0)}, {"temperature_k": np.full((2, 5), np.nan)},
    {"source_kind": "unlabeled"}, {"source_kind": "flash", "flash_version": ""},
    {"run_id": " "}, {"source_reference": ""}, {"surface_region": ""},
])
def test_invalid_history_is_rejected(overrides):
    with pytest.raises(ValueError):
        history(**overrides)


@pytest.mark.parametrize("edge,center", [(0.0, 0.2), (-0.1, 0.2), (0.1, 0.0), (0.3, 0.5), (np.nan, 0.2), (0.1, np.inf)])
def test_invalid_or_overlapping_regions_are_rejected(edge, center):
    with pytest.raises(ValueError):
        history().thermal_contrast(edge_fraction=edge, center_fraction=center)


@pytest.mark.parametrize("mutation", ["units", "missing_units", "schema", "coordinate", "motion", "missing_region", "missing_dataset", "shape"])
def test_malformed_hdf5_rejected_before_use(tmp_path, mutation):
    target = tmp_path / "surface.h5"
    write_surface_history(target, history())
    with h5py.File(target, "a") as handle:
        if mutation == "units":
            handle["x_edges_m"].attrs["units"] = "cm"
        elif mutation == "missing_units":
            del handle["temperature_k"].attrs["units"]
        elif mutation == "schema":
            handle.attrs["schema"] = "flash raw checkpoint"
        elif mutation == "coordinate":
            handle.attrs["coordinate_system"] = "cylindrical"
        elif mutation == "motion":
            handle.attrs["surface_model"] = "moving"
        elif mutation == "missing_region":
            del handle.attrs["surface_region"]
        elif mutation == "missing_dataset":
            del handle["time_s"]
        elif mutation == "shape":
            del handle["temperature_k"]
            handle.create_dataset("temperature_k", data=np.ones((2, 4))).attrs["units"] = "K"
    with pytest.raises(ValueError):
        load_surface_history(target)


def test_missing_or_non_hdf5_file_raises_value_error(tmp_path):
    with pytest.raises(ValueError):
        load_surface_history(tmp_path / "absent.h5")
    target = tmp_path / "text.h5"
    target.write_text("not HDF5")
    with pytest.raises(ValueError):
        load_surface_history(target)


def test_cli_emits_explicitly_labeled_json(tmp_path):
    target = tmp_path / "synthetic.h5"
    write_surface_history(target, history())
    script = Path(__file__).resolve().parents[1] / "tools" / "inspect_flash_surface.py"
    result = subprocess.run([sys.executable, str(script), str(target), "--time-s", "1.0"], check=True, text=True, capture_output=True)
    report = json.loads(result.stdout)
    assert report["source_kind"] == "synthetic"
    assert report["metrics"]["time_s"] == 1.0
    assert report["metrics"]["left_edge_mean_temperature_k"] == pytest.approx(150.0)
