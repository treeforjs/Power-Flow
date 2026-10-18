"""Analytic structured-block fixtures, explicitly labeled synthetic."""

from dataclasses import replace
import hashlib
import json

import h5py
import numpy as np
import pytest

from mhdlab.flash_extract import EXTRACTION_SCHEMA, FlashLayout, FlashSnapshot, PlanarFace, extract_surface_history
from mhdlab.flash_surface import SCHEMA, load_surface_history, write_surface_history


def layout(**kw):
    values = dict(bounding_box_dataset="bounding box", node_type_dataset="node type",
                  temperature_dataset="lattice", material_fraction_dataset="foil",
                  bounding_box_axes=("x", "y"), scalar_axes=("block", "y", "x"),
                  physical_axes={"x": "width", "y": "normal"}, leaf_node_value=1,
                  nonleaf_node_values=(2,), coordinate_unit="m", coordinate_scale_to_m=1,
                  temperature_unit="K", temperature_scale_to_k=1, temperature_offset_k=0,
                  min_material_fraction=1 - 1e-8)
    values.update(kw)
    return FlashLayout(**values)


def face(**kw):
    values = dict(x_edges_m=(0, .25, .75, 1), face_normal_coordinate_m=1,
                  outward_sign=1, material_sample_depth_m=.25, current_interval_m=(0, 2))
    values.update(kw)
    return PlanarFace(**values)


def write(path, boxes, temperature, nodes=None, material=None):
    with h5py.File(path, "w") as f:
        f["bounding box"] = np.asarray(boxes, dtype=float)
        f["node type"] = np.ones(len(boxes), dtype=int) if nodes is None else nodes
        f["lattice"] = np.asarray(temperature, dtype=float)
        f["foil"] = np.ones_like(temperature) if material is None else material


def run(tmp_path, boxes, temperature, *, raw_layout=None, sample_face=None, nodes=None, material=None):
    snapshots = []
    for i in range(2):
        path = tmp_path / f"snapshot_{i}.h5"
        write(path, boxes, temperature, nodes, material)
        snapshots.append(FlashSnapshot(path, i, "ns", 1e-9))
    return extract_surface_history(snapshots, raw_layout or layout(), sample_face or face(),
                                   source_kind="synthetic", run_id="analytic-block-fixture",
                                   source_reference="manufactured structured data", surface_region="test face")


def test_unequal_amr_leaf_widths_project_by_overlap_and_parent_is_ignored(tmp_path):
    # Two left cells cover [0,.5], two right cells cover [.5,1]. A parent with
    # an invalid scalar is deliberately not used in the leaf projection.
    boxes = [[[0, .5], [0, 1]], [[.5, 1], [0, 1]], [[0, 1], [0, 1]]]
    temperature = [[[100, 200]], [[300, 500]], [[np.nan, np.nan]]]
    history, manifest = run(tmp_path, boxes, temperature, nodes=[1, 1, 2])
    np.testing.assert_allclose(history.temperature_k, [[100, 250, 500]] * 2)
    assert manifest["source_kind"] == "synthetic"
    assert manifest["schema"] == EXTRACTION_SCHEMA
    assert manifest["surface_history_schema"] == SCHEMA
    assert history.temperature_location == "material_sample_plane"
    assert history.material_sample_depth_m == .25
    assert history.boundary_temperature_reference == ""
    output = tmp_path / "plane.h5"
    write_surface_history(output, history)
    reread = load_surface_history(output)
    assert reread.temperature_location == "material_sample_plane"
    assert reread.material_sample_depth_m == .25
    assert reread.boundary_temperature_reference == ""
    assert reread.current_averaging_interval_m is None
    assert reread.current_uniformity_reference == ""
    assert manifest["snapshots"][0]["selected_leaf_blocks"] == [0, 1]
    assert manifest["snapshots"][0]["sha256"] == hashlib.sha256((tmp_path / "snapshot_0.h5").read_bytes()).hexdigest()
    assert manifest["sampling"]["sample_normal_coordinate_m"] == .75
    assert "not_boundary" in manifest["diagnostic_kind"]
    assert json.loads(json.dumps(manifest))["layout"]["temperature_dataset"] == "lattice"


def test_normal_boundary_is_half_open_and_never_averaged_or_nearest_guessed(tmp_path):
    boxes = [[[0, 1], [0, .5]], [[0, 1], [.5, 1]]]
    history, _ = run(tmp_path, boxes, [[[100]], [[800]]], sample_face=face(material_sample_depth_m=.5))
    np.testing.assert_allclose(history.temperature_k, 800)


def test_axis_permutation_units_temperature_offset_and_signed_coordinates(tmp_path):
    # Raw x is normal, raw y is reversed width in centimetres. Raw storage is
    # (raw_y, block, raw_x), with two raw_x cells. Plane is raw_x=75cm.
    raw_layout = layout(bounding_box_axes=("y", "x"), scalar_axes=("y", "block", "x"),
                        physical_axes={"y": "width", "x": "normal"}, axis_signs={"y": -1, "x": 1},
                        axis_offsets_m={"y": 1, "x": 0}, coordinate_unit="cm", coordinate_scale_to_m=.01,
                        temperature_unit="degC", temperature_offset_k=273.15)
    history, manifest = run(tmp_path, [[[0, 100], [0, 100]]], [[[1, 10]], [[2, 20]]],
                            raw_layout=raw_layout)
    np.testing.assert_allclose(history.temperature_k, [[293.15, 288.15, 283.15]] * 2)
    assert manifest["layout"]["coordinate_scale_to_m"] == .01


def test_3d_projection_averages_resolved_current_extent(tmp_path):
    raw_layout = layout(bounding_box_axes=("x", "y", "z"), scalar_axes=("block", "z", "y", "x"),
                        physical_axes={"x": "width", "y": "normal", "z": "current"})
    boxes = [[[0, 1], [0, 1], [0, 2]]]
    # Two current cells, each spanning half the requested [0.5,1.5] interval.
    history, manifest = run(tmp_path, boxes, [[[[100, 200]], [[300, 600]]]],
                            raw_layout=raw_layout, sample_face=face(current_interval_m=(.5, 1.5)))
    np.testing.assert_allclose(history.temperature_k, [[200, 300, 400]] * 2)
    assert manifest["sampling"]["current_model"] == "resolved_current_interval_average_3d"
    assert history.current_averaging_interval_m == (.5, 1.5)
    assert history.current_uniformity_reference == ""
    path = tmp_path / "three_dimensional_average.h5"
    write_surface_history(path, history)
    assert load_surface_history(path).current_averaging_interval_m == (.5, 1.5)


def test_reversed_normal_axis_and_negative_outward_direction_sample_material_side(tmp_path):
    raw_layout = layout(axis_signs={"x": 1, "y": -1}, axis_offsets_m={"x": 0, "y": 0})
    sample_face = face(face_normal_coordinate_m=0, outward_sign=-1, material_sample_depth_m=.25)
    history, _ = run(tmp_path, [[[0, 1], [-1, 0]]], [[[100], [800]]],
                     raw_layout=raw_layout, sample_face=sample_face)
    np.testing.assert_allclose(history.temperature_k, 800)


def test_2d_supports_only_explicit_singleton_inactive_storage_axis(tmp_path):
    raw_layout = layout(bounding_box_axes=("x", "y", "z"), scalar_axes=("block", "z", "y", "x"))
    boxes = [[[0, 1], [0, 1], [0, 0]]]
    history, _ = run(tmp_path, boxes, [[[[100, 300]]]], raw_layout=raw_layout)
    np.testing.assert_allclose(history.temperature_k, [[100, 200, 300]] * 2)
    with pytest.raises(ValueError, match="inactive storage"):
        run(tmp_path, boxes, [[[[100, 300]], [[200, 400]]]], raw_layout=raw_layout)


@pytest.mark.parametrize("boxes,nodes,match", [
    ([[[0, .4], [0, 1]], [[.5, 1], [0, 1]]], [1, 1], "hole|incomplete"),
    ([[[0, .6], [0, 1]], [[.5, 1], [0, 1]]], [1, 1], "overlap"),
    ([[[0, .5], [0, 1]], [[.5, 1], [0, 1]]], [1, 2], "hole|incomplete"),
    ([[[0, 1], [0, 1]], [[0, 1], [0, 1]]], [3, 2], "unconfigured"),
])
def test_missing_ambiguous_or_nonleaf_coverage_rejected(tmp_path, boxes, nodes, match):
    with pytest.raises(ValueError, match=match):
        run(tmp_path, boxes, [[[100]], [[200]]], nodes=nodes)


def test_equal_total_area_cannot_hide_a_3d_hole_and_overlap(tmp_path):
    raw_layout = layout(bounding_box_axes=("x", "y", "z"), scalar_axes=("block", "z", "y", "x"),
                        physical_axes={"x": "width", "y": "normal", "z": "current"})
    boxes = [[[0, .5], [0, 1], [0, 1]], [[0, .5], [0, 1], [0, 1]],
             [[.5, 1], [0, 1], [0, 2]]]
    with pytest.raises(ValueError, match="overlap"):
        run(tmp_path, boxes, [[[[100]]]] * 3, raw_layout=raw_layout)


@pytest.mark.parametrize("direction", ["width", "current"])
def test_geometry_tolerance_cannot_hide_alternating_holes_and_duplicate_slabs(tmp_path, direction):
    # Every even slab is covered twice and every odd slab is missing. Total
    # area is correct, and each slab is smaller than the declared tolerance.
    # Topology must still reject the duplicated and absent leaf coverage.
    count = 256
    boxes = []
    for i in range(0, count, 2):
        slab = [i / count, (i + 1) / count]
        box = [slab, [0, 1]] if direction == "width" else [[0, 1], [0, 1], slab]
        boxes.extend([box, box])
    sample = face(x_edges_m=(0, .25, .75, 1), current_interval_m=(0, 1),
                  geometry_atol_m=1 / 128)
    if direction == "width":
        raw_layout, temperature = layout(), [[[100]]] * count
    else:
        raw_layout = layout(bounding_box_axes=("x", "y", "z"),
                            scalar_axes=("block", "z", "y", "x"),
                            physical_axes={"x": "width", "y": "normal", "z": "current"})
        temperature = [[[[100]]]] * count
    with pytest.raises(ValueError, match="overlap|hole|incomplete"):
        run(tmp_path, boxes, temperature, raw_layout=raw_layout, sample_face=sample)


@pytest.mark.parametrize("temperature,material", [
    (0, 1), (-1, 1), (np.inf, 1), (np.nan, 1),
    (100, 0), (100, .9), (100, np.nan), (100, 1.1),
])
def test_selected_cells_must_be_positive_temperature_and_pure_material(tmp_path, temperature, material):
    with pytest.raises(ValueError):
        run(tmp_path, [[[0, 1], [0, 1]]], [[[temperature]]], material=[[[material]]])


def test_vacuum_outside_requested_plane_is_not_sampled(tmp_path):
    history, _ = run(tmp_path, [[[0, 1], [0, 1]]], [[[np.nan], [500]]], material=[[[0], [1]]])
    np.testing.assert_allclose(history.temperature_k, 500)


def test_conflicting_unit_metadata_and_external_links_rejected(tmp_path):
    for i in range(2):
        write(tmp_path / f"s{i}.h5", [[[0, 1], [0, 1]]], [[[100]]])
    inputs = [FlashSnapshot(tmp_path / f"s{i}.h5", i, "s", 1) for i in range(2)]
    options = dict(source_kind="synthetic", run_id="test", source_reference="fixture", surface_region="face")
    with h5py.File(inputs[0].path, "a") as f:
        f["lattice"].attrs["units"] = "eV"
    with pytest.raises(ValueError, match="unit attribute"):
        extract_surface_history(inputs, layout(), face(), **options)
    with h5py.File(inputs[0].path, "a") as f:
        del f["lattice"]
        f["lattice"] = h5py.ExternalLink(str(inputs[1].path), "/lattice")
    with pytest.raises(ValueError, match="external links"):
        extract_surface_history(inputs, layout(), face(), **options)


@pytest.mark.parametrize("mutation", ["node_rank", "node_float", "scalar_shape", "boolean_temperature", "virtual", "soft_link"])
def test_unsupported_raw_layouts_and_unhashed_dependencies_rejected(tmp_path, mutation):
    inputs = []
    for i in range(2):
        path = tmp_path / f"s{i}.h5"
        write(path, [[[0, 1], [0, 1]]], [[[100]]])
        inputs.append(FlashSnapshot(path, i, "s", 1))
    with h5py.File(inputs[0].path, "a") as f:
        if mutation.startswith("node"):
            del f["node type"]
            f["node type"] = [[1]] if mutation == "node_rank" else [1.0]
        elif mutation == "scalar_shape":
            del f["foil"]
            f["foil"] = np.ones((1, 2, 1))
        else:
            del f["lattice"]
            if mutation == "boolean_temperature":
                f["lattice"] = np.ones((1, 1, 1), dtype=bool)
            elif mutation == "soft_link":
                f["alias"] = [[[100.0]]]
                f["lattice"] = h5py.SoftLink("/alias")
            elif mutation == "virtual":
                virtual = h5py.VirtualLayout((1, 1, 1), dtype=float)
                virtual[:] = h5py.VirtualSource(str(inputs[1].path), "lattice", (1, 1, 1))
                f.create_virtual_dataset("lattice", virtual)
    with pytest.raises(ValueError):
        extract_surface_history(inputs, layout(), face(), source_kind="synthetic", run_id="fixture",
                                source_reference="fixture", surface_region="face")


@pytest.mark.parametrize("change", [
    {"geometry": "cylindrical"}, {"data_location": "node_centered"},
    {"physical_axes": {"x": "normal", "y": "current"}},
    {"scalar_axes": ("block", "x", "x")}, {"leaf_node_value": 2},
    {"coordinate_scale_to_m": 0}, {"temperature_scale_to_k": -1},
    {"min_material_fraction": .5}, {"axis_signs": {"x": 0, "y": 1}},
])
def test_unsupported_or_ambiguous_layout_configuration_rejected(change):
    with pytest.raises(ValueError):
        layout(**change)


@pytest.mark.parametrize("change", [
    {"material_sample_depth_m": 0}, {"material_sample_depth_m": -1},
    {"outward_sign": 0}, {"current_interval_m": (1, 1)},
    {"x_edges_m": (0, .5, .5, 1)}, {"geometry_atol_m": .2},
])
def test_invalid_face_configuration_rejected(change):
    with pytest.raises(ValueError):
        face(**change)


def test_time_conversion_order_source_classification_and_duplicate_paths_are_explicit(tmp_path):
    path = tmp_path / "a.h5"
    write(path, [[[0, 1], [0, 1]]], [[[100]]])
    inputs = [FlashSnapshot(path, 1, "ns", 1e-9), FlashSnapshot(path, 2, "ns", 1e-9)]
    options = dict(source_kind="synthetic", run_id="fixture", source_reference="fixture", surface_region="face")
    with pytest.raises(ValueError, match="distinct"):
        extract_surface_history(inputs, layout(), face(), **options)
    with pytest.raises(ValueError, match="strictly increasing"):
        extract_surface_history(inputs[::-1], layout(), face(), **options)
    with pytest.raises(ValueError, match="flash_version"):
        extract_surface_history(inputs, layout(), face(), **{**options, "source_kind": "flash"})
    with pytest.raises(ValueError):
        FlashSnapshot(path, -1, "s", 1)
    assert replace(inputs[0], time_value=-1, time_offset_s=2e-9).time_s == pytest.approx(1e-9)
