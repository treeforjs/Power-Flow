import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from mhdlab.thermal_expansion_source import load_matched_history


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value))


def reseal(root):
    """Seal a manufactured fixture; never applied to science source inputs."""
    for name in ("foil_snapshots.json",):
        path = root / name
        meta = json.loads(path.read_text())
        meta["csv_sha256"] = digest(root / "foil_snapshots.csv")
        write_json(path, meta)
    path = root / "source_run.json"
    run = json.loads(path.read_text()) if path.exists() else {
        "schema": "mhdlab-mfem-edge-current-run-v1", "passed": True}
    run["configuration"] = json.loads((root / "configuration.json").read_text())
    run["outputs_sha256"] = {p.name: digest(p) for p in root.iterdir() if p.is_file() and p != path}
    write_json(path, run)
    return digest(path)


@pytest.fixture
def native(tmp_path):
    root = tmp_path / "manufactured_native"
    root.mkdir()
    times = [0.0, 1e-9]
    config = {"width_cells": 2, "thickness_cells": 2, "foil_width_m": .005,
              "foil_thickness_m": .0002, "gap_m": .005, "current_length_m": .025,
              "initial_temperature_k": 300.0}
    write_json(root / "configuration.json", config)
    bulk = ["time_s,x_left_m,x_right_m,depth_left_m,depth_right_m,temperature_k,current_density_a_m2,joule_power_w_m3"]
    for it, t in enumerate(times):
        for ix in (1, 0):  # Deliberately not row-major order.
            for iy in (1, 0):
                bulk.append(f"{t},{-.0025 + ix * .0025},{-.0025 + (ix+1)*.0025},{iy*.0001},{(iy+1)*.0001},{300+it*(10*iy+ix+1)},1,1")
    (root / "foil_snapshots.csv").write_text("\n".join(bulk) + "\n")
    write_json(root / "foil_snapshots.json", {
        "schema": "mhdlab-mfem-foil-snapshots-v1", "sampling": "cell_average",
        "units": {"temperature_k": "K"},
        "geometry": {"foil": "upper", "depth_reference": "inward_from_upper_gap_face"},
        "field_snapshot_times_s": times})
    (root / "em_mesh.mesh").write_text("manufactured mesh\n")
    (root / "em_axes.csv").write_text("manufactured axes\n")
    fields = {"schema": "mhdlab-mfem-Az-field-history-v1", "potential_units": "T m", "rate_units": "V/m",
              "mesh_file": "em_mesh.mesh", "axes_file": "em_axes.csv",
              "mesh_sha256": digest(root / "em_mesh.mesh"), "axes_sha256": digest(root / "em_axes.csv"), "snapshots": []}
    for it, t in enumerate(times):
        entry = {"time_s": t}
        for key in ("potential", "rate"):
            name = f"{key}_{it}.gf"
            (root / name).write_text(f"manufactured {key} {it}\n")
            entry[key + "_file"] = name
            entry[key + "_sha256"] = digest(root / name)
        fields["snapshots"].append(entry)
    write_json(root / "em_field_snapshots.json", fields)
    # Distinct surface fields test reflection orientation as well as offset.
    top = [[0, 0], [1, 2]]
    bottom = [[0, 0], [3, 4]]
    side = [[0, 0], [5, 6]]
    for name, values in {"top": top, "bottom": bottom, "left": side, "right": side,
                         "anode_top": bottom, "anode_bottom": top,
                         "anode_left": np.asarray(side)[:, ::-1], "anode_right": np.asarray(side)[:, ::-1]}.items():
        text = "time_s,z_index,x_index,temperature_value\n"
        for it, t in enumerate(times):
            for ix in range(2):
                text += f"{t},0,{ix},{values[it][ix]}\n"
        (root / (name + ".csv")).write_text(text)
        write_json(root / (name + "_surface.json"), {
            "schema": "mhdlab-mfem-surface-samples-v1", "time_s": times,
            "temperature": {"value_basis": "temperature_increment_kelvin", "initial_temperature_k": 300.0},
            "table": {"path": name + ".csv", "sha256": digest(root / (name + ".csv"))}})
    return root, reseal(root)


def test_native_order_absolute_bulk_and_missing_lower(native):
    root, stamp = native
    history = load_matched_history(root, expected_source_sha256=stamp)
    np.testing.assert_array_equal(history.upper_temperature_k[0], 300)
    np.testing.assert_array_equal(history.upper_temperature_k[1], [[301, 302], [311, 312]])
    np.testing.assert_array_equal(history.face_temperature_k["top"][1], [301, 302])
    assert history.lower_temperature_k is None
    assert len(history.native_fields) == 2
    assert history.manifest["native_fields_verified"]
    assert history.manifest["interpolation"].startswith("none")
    with pytest.raises(ValueError):
        history.upper_temperature_k[0, 0, 0] = 900


def test_explicit_lower_reconstruction_and_y_coordinates(native):
    root, stamp = native
    history = load_matched_history(root, expected_source_sha256=stamp, lower_volume="mirror_with_surface_validation")
    np.testing.assert_array_equal(history.lower_temperature_k[1], [[311, 312], [301, 302]])
    np.testing.assert_allclose(history.upper_y_edges_m, [.0025, .0026, .0027])
    np.testing.assert_allclose(history.lower_y_edges_m, [-.0027, -.0026, -.0025])
    assert "unverified" in history.manifest["lower_volume_basis"]
    assert max(history.manifest["surface_symmetry_max_abs_k"].values()) == 0


def test_admissibility_reports_and_preserves_hot_cells(native):
    root, stamp = native
    history = load_matched_history(root, expected_source_sha256=stamp)
    report = history.admissibility(300, 310)
    assert report["first_invalid_bulk_time_s"] == 1e-9
    assert report["invalid_upper_cell_counts"] == [0, 2]
    assert history.upper_temperature_k.max() == 312
    assert not report["all_bulk_frames_admissible"]
    assert report["first_invalid_surface_time_s"] is None
    assert report["last_common_safe_bulk_time_s"] == 0
    assert not report["all_exported_temperatures_admissible"]
    assert history.admissibility(300, 305)["first_invalid_surface_time_s"] == 1e-9
    assert history.admissibility(300, 312)["all_bulk_frames_admissible"]


def test_surface_violation_blocks_later_cooler_history(native):
    root, stamp = native
    history = load_matched_history(root, expected_source_sha256=stamp)
    faces = dict(history.face_temperature_k)
    face = faces["top"].copy()
    face[0, 0] = 400
    faces["top"] = face
    report = replace(history, face_temperature_k=faces).admissibility(293.15, 350)
    assert report["all_bulk_frames_admissible"]
    assert report["first_invalid_surface_time_s"] == 0
    assert report["last_common_safe_bulk_time_s"] is None
    assert not report["all_exported_temperatures_admissible"]


@pytest.mark.parametrize("which", ["potential_1.gf", "rate_1.gf", "em_mesh.mesh", "foil_snapshots.csv", "anode_top.csv"])
def test_hash_tamper_rejected(native, which):
    root, stamp = native
    with (root / which).open("a") as stream:
        stream.write("tamper\n")
    with pytest.raises(ValueError, match="SHA256"):
        load_matched_history(root, expected_source_sha256=stamp)


def test_old_or_other_source_rejected_by_default(native):
    with pytest.raises(ValueError, match="requested matched run"):
        load_matched_history(native[0])


def test_field_time_mismatch_rejected(native):
    root, _ = native
    path = root / "em_field_snapshots.json"
    fields = json.loads(path.read_text())
    fields["snapshots"][1]["time_s"] = 2e-9
    write_json(path, fields)
    with pytest.raises(ValueError, match="different times"):
        load_matched_history(root, expected_source_sha256=reseal(root))


def test_missing_bulk_cell_rejected(native):
    root, _ = native
    path = root / "foil_snapshots.csv"
    path.write_text("\n".join(path.read_text().splitlines()[:-1]) + "\n")
    with pytest.raises(ValueError, match="missing native bulk"):
        load_matched_history(root, expected_source_sha256=reseal(root))


def test_duplicate_bulk_cell_rejected(native):
    root, _ = native
    path = root / "foil_snapshots.csv"
    rows = path.read_text().splitlines()
    rows[-1] = rows[-2]
    path.write_text("\n".join(rows) + "\n")
    with pytest.raises(ValueError, match="duplicate or missing"):
        load_matched_history(root, expected_source_sha256=reseal(root))


def test_surface_asymmetry_only_blocks_reconstruction(native):
    root, _ = native
    path = root / "anode_top.csv"
    path.write_text(path.read_text().replace("1e-09,0,0,3", "1e-09,0,0,9"))
    meta_path = root / "anode_top_surface.json"
    meta = json.loads(meta_path.read_text())
    meta["table"]["sha256"] = digest(path)
    write_json(meta_path, meta)
    stamp = reseal(root)
    history = load_matched_history(root, expected_source_sha256=stamp)
    assert history.manifest["surface_symmetry_max_abs_k"]["upper_bottom_lower_top_k"] == 6
    with pytest.raises(ValueError, match="mirror tolerance"):
        load_matched_history(root, expected_source_sha256=stamp, lower_volume="mirror_with_surface_validation")


@pytest.mark.parametrize("mode", ["mirror", "native", "", None])
def test_implicit_or_unknown_lower_mode_rejected(native, mode):
    with pytest.raises(ValueError, match="explicitly select"):
        load_matched_history(native[0], expected_source_sha256=native[1], lower_volume=mode)


def test_wrong_geometry_rejected(native):
    root, _ = native
    path = root / "configuration.json"
    config = json.loads(path.read_text())
    config["current_length_m"] = .01
    write_json(path, config)
    with pytest.raises(ValueError, match="geometry"):
        load_matched_history(root, expected_source_sha256=reseal(root))
