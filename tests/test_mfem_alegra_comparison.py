"""Conditional readout alignment must not become fitted/shifted validation."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("_mfem_alegra_comparison_test", ROOT / "tools/compare_mfem_alegra.py")
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


@pytest.fixture
def snapshots(tmp_path):
    path = tmp_path / "foil_snapshots.csv"
    x = [-.0025, -.001, 0., .001, .0025]
    depth = [0., 10e-6, 25e-6, .0002]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(comparison.COLUMNS)
        for time in (40e-9, 60e-9):
            for iz in range(3):
                for ix in range(4):
                    # Authored synthetic comparison fixture, not a model run.
                    temperature = 300 + 100*iz + (400 if ix in (0, 3) else 0) + time*1e9
                    writer.writerow([time, x[ix], x[ix+1], depth[iz], depth[iz+1], temperature, 1e10, 2e12])
    metadata = {"schema": comparison.SCHEMA, "sampling": "cell_average", "status": "manufactured",
        "units": comparison.UNITS, "csv_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "geometry": {"foil": "upper", "depth_reference": "inward_from_upper_gap_face",
                     "foil_width_m": .005, "foil_thickness_m": .0002, "gap_m": .005},
        "provenance": "Authored test table, no MFEM calculation or ALEGRA reproduction"}
    path.with_suffix(".json").write_text(json.dumps(metadata), encoding="utf-8")
    return path


def test_cell_average_map_and_literal_half_open_sampling(snapshots):
    groups, _ = comparison.load_snapshots(snapshots)
    _, grid = comparison.exact_snapshot(groups, 40e-9)
    np.testing.assert_array_equal(comparison.p0_indices(grid["depth_edges_m"], [0., 10e-6, .0002]), [0, 1, 2])
    center, edge, band = comparison.lateral_profile_columns(grid, .1)
    assert center.tolist() == [1, 2]
    assert edge.tolist() == [0, 3]
    assert band.tolist() == [0, 3]
    assert grid["temperature_k"].shape == (3, 4)


def test_missing_snapshot_and_outside_depth_do_not_interpolate(snapshots):
    groups, _ = comparison.load_snapshots(snapshots)
    with pytest.raises(ValueError, match="temporal interpolation"):
        comparison.exact_snapshot(groups, 50e-9)
    with pytest.raises(ValueError, match="extrapolation"):
        comparison.p0_indices(groups[40e-9]["depth_edges_m"], [-1e-9])


@pytest.mark.parametrize("change,match", [
    (lambda m: m.update(csv_sha256="0"*64), "SHA256"),
    (lambda m: m.update(sampling="point"), "cell-average"),
    (lambda m: m["units"].update(temperature_k="eV"), "SI units"),
    (lambda m: m["geometry"].update(depth_reference="from_back_face"), "inward-depth"),
    (lambda m: m["geometry"].update(foil_thickness_m=.0004), "geometry"),
])
def test_mismatched_quantity_or_geometry_rejects(snapshots, change, match):
    sidecar = snapshots.with_suffix(".json")
    metadata = json.loads(sidecar.read_text())
    change(metadata)
    sidecar.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match=match):
        comparison.load_snapshots(snapshots)


def test_missing_cell_and_duplicate_cell_cannot_hide_in_reconstruction(snapshots):
    lines = snapshots.read_text().splitlines()
    lines[-1] = lines[-2]
    snapshots.write_text("\n".join(lines)+"\n")
    sidecar = snapshots.with_suffix(".json")
    metadata = json.loads(sidecar.read_text())
    metadata["csv_sha256"] = comparison.sha(snapshots)
    sidecar.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="duplicate|complete"):
        comparison.load_snapshots(snapshots)


def test_declared_edge_sensitivity_does_not_change_middle_prediction(snapshots):
    groups, _ = comparison.load_snapshots(snapshots)
    grid = groups[40e-9]
    data = comparison.rows(comparison.DATA / "figure_4_28_temperature.csv")
    result, sampling = comparison.depth_comparison(grid, data, .1)
    first_edge = next(r for r in result if r["profile"] == "foil_edge")
    first_middle = next(r for r in result if r["profile"] == "foil_middle")
    assert first_edge["mfem_p0_temperature_k"] == 740
    assert first_middle["mfem_p0_temperature_k"] == 340
    assert sampling["center_column_indices"] == [1, 2]
    assert not comparison.summarize(result)["statistical_likelihood"]


def test_actual_uniform_joule_field_is_not_replaced_by_temperature_contrast(snapshots):
    groups, _ = comparison.load_snapshots(snapshots)
    contrast = comparison.heating_contrast(groups[40e-9])
    assert contrast["joule_power_w_m3"]["edge_over_center"] == 1.
    assert contrast["temperature_k"]["edge_minus_center"] == 400.
    groups[40e-9]["joule_power_w_m3"][:, 1:3] = 0.
    assert comparison.heating_contrast(groups[40e-9])["joule_power_w_m3"]["edge_over_center"] is None


def test_complete_comparison_is_explicitly_conditional_not_a_physical_fit(snapshots, tmp_path):
    output = tmp_path / "comparison"
    result = comparison.compare(snapshots, output)
    assert result["physical_alegra_reproduction_established"] is False
    assert result["status"] == "provisional_conditional_comparison"
    assert len(result["width_comparison"]["readouts"]) == 137
    assert all(item["summary"]["retained_points"] == 76 for item in result["depth_comparisons"])
    assert "k_B*T" in result["limits"][4]
    assert (output / "mfem_alegra_comparison.png").is_file()
    assert (output / "comparison.json").is_file()
