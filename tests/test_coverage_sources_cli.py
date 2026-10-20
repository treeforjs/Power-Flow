"""Independent inventory, source-file and provenance checks for coverage sweeps."""

import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np
import pytest

from mhdlab.temkin_planar_io import load_temkin_planar_release_history
from mhdlab.planar_release_io import load_planar_surface_history, plain_metadata


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/build_coverage_sources.py"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def config_file(directory):
    config = json.loads((ROOT / "configs/coverage_source_manufactured.json").read_bytes())
    config["geometry"].update(width_bins=3, current_bins=2, snapshot_intervals=2)
    config["subdivisions_per_interval"] = 1
    # Retain the supplied nonzero allocation convention while exercising zero.
    config["coverage_ml"] = [0., .25, .5]
    path = directory / "parameters.json"
    path.write_text(json.dumps(config) + "\n")
    return path, config


def execute(output, *arguments):
    return subprocess.run([sys.executable, str(TOOL), *map(str, arguments),
                           "--output-dir", str(output)], capture_output=True, text=True)


@pytest.fixture(scope="module")
def demo_sweep(tmp_path_factory):
    directory = tmp_path_factory.mktemp("coverage-source-cli")
    path, config = config_file(directory)
    output = directory / "sources"
    process = execute(output, "--demo", "--config", path, "--plot")
    assert process.returncode == 0, process.stderr
    return path, config, output, json.loads((output / "report.json").read_bytes())


def test_exported_files_keep_exact_inputs_hashes_and_manufactured_labels(demo_sweep):
    path, config, output, report = demo_sweep
    assert report["schema"] == "mhdlab-coverage-source-run-v1"
    assert report["parameter_status"] == "manufactured"
    assert report["parameter_provenance"] == config["parameter_provenance"]
    assert report["input_config_sha256"] == sha(path) == sha(output / "input_config.json")
    assert report["input_surface_sha256"] is None
    assert (output / "coverage_sources.png").is_file()
    assert report["coverage_source_sweep"]
    assert len(report["sources"]) == len(config["coverage_ml"])
    assert set(report["output_sha256"]) == {item.name for item in output.iterdir() if item.name != "report.json"}
    for name, digest in report["output_sha256"].items():
        assert sha(output / name) == digest
    for name, digest in report["source_sha256"].items():
        assert sha(ROOT / name) == digest
    for index, source in enumerate(report["sources"]):
        assert source["index"] == index
        assert source["coverage_ml"] == config["coverage_ml"][index]
        assert source["source_file"] == f"source_{index:04d}.h5"
        assert source["source_sha256"] == sha(output / source["source_file"])
    assert not list(output.parent.glob(".coverage-*"))


def test_one_parent_inventory_is_partitioned_without_duplicate_layers(demo_sweep):
    _, config, output, report = demo_sweep
    thermal_history = None
    for source in report["sources"]:
        history = load_temkin_planar_release_history(output / source["source_file"])
        shape = history.patch_areas_m2.shape
        profile = np.broadcast_to(np.asarray(config["parent"]["patch_multiplier"]), shape)
        parent = source["coverage_ml"] * config["parent"]["monolayer"]["areal_density_m2"] * profile
        expected_parent = np.sum(parent * history.patch_areas_m2)
        assert source["initial_parent_number"] == pytest.approx(expected_parent, rel=3e-14, abs=0.)
        fractions = np.array([np.broadcast_to(np.asarray(record["fraction"]), shape)
                              for record in config["allocations"]])
        np.testing.assert_allclose(fractions.sum(axis=0), 1., rtol=0., atol=1e-15)
        expected_pool = parent[None, :, :] * fractions
        np.testing.assert_allclose(history.coverage_m2[0], expected_pool, rtol=3e-14, atol=0.)
        np.testing.assert_allclose(history.coverage_m2[0].sum(axis=0), parent, rtol=3e-14, atol=0.)
        initial = np.sum(expected_pool * history.patch_areas_m2[None, :, :], axis=(1, 2))
        released = np.sum(history.released_m2 * history.patch_areas_m2[None, None, :, :], axis=(0, 2, 3))
        remaining = np.sum(history.coverage_m2[-1] * history.patch_areas_m2[None, :, :], axis=(1, 2))
        debit = np.sum(history.total_wall_energy_debit_j_m2 * history.patch_areas_m2[None, None, :, :], axis=(0, 2, 3))
        for name, expected in (("initial_pool_number", initial), ("released_pool_number", released),
                               ("remaining_pool_number", remaining), ("total_wall_energy_debit_j", debit)):
            np.testing.assert_allclose(source[name], expected, rtol=3e-14, atol=0.)
        np.testing.assert_allclose(released + remaining, initial, rtol=3e-14, atol=0.)
        np.testing.assert_allclose(np.sum(initial), expected_parent, rtol=3e-14, atol=0.)
        for pool, allocation in zip(history.reservoirs, config["allocations"]):
            assert asdict(pool) == allocation["reservoir"]
            assert pool.reference_coverage_m2 == allocation["reservoir"]["reference_coverage_m2"]
        if thermal_history is None:
            thermal_history = history
        else:
            for name in ("temperature_k", "time_edges_s", "x_edges_m", "z_edges_m"):
                np.testing.assert_array_equal(getattr(history, name), getattr(thermal_history, name))
        if source["coverage_ml"] == 0:
            for name in ("coverage_m2", "released_m2", "translational_energy_j_m2",
                         "internal_energy_j_m2", "binding_energy_j_m2", "total_wall_energy_debit_j_m2"):
                np.testing.assert_array_equal(getattr(history, name), 0.)


def test_matched_width_controls_preserve_each_current_time_state_and_energy_mixture(demo_sweep):
    _, _, output, report = demo_sweep
    boltzmann = 1.380649e-23
    for source in report["sources"]:
        path = output / source["source_file"]
        history = load_temkin_planar_release_history(path)
        width = np.diff(history.x_edges_m)
        # The original local temperature mixture remains explicit in the source.
        assert np.ptp(history.temperature_k[-1], axis=-1).max() > 0
        np.testing.assert_allclose(history.translational_energy_j_m2,
                                   history.released_m2 * (2 * boltzmann * history.temperature_k[:, None]),
                                   rtol=3e-14, atol=0.)
        with h5py.File(path, "r") as stream:
            for name in ("released_m2", "translational_energy_j_m2", "internal_energy_j_m2",
                         "binding_energy_j_m2", "total_wall_energy_debit_j_m2"):
                resolved = getattr(history, name)
                uniform = stream[f"matched_uniform_width/{name}"][...]
                np.testing.assert_allclose(uniform, np.broadcast_to(uniform[..., :1], uniform.shape),
                                           rtol=3e-14, atol=0.)
                np.testing.assert_allclose(np.sum(resolved * width, axis=-1),
                                           np.sum(uniform * width, axis=-1), rtol=3e-14, atol=0.)


def test_summary_csv_reproduces_report_and_reopened_hdf_counts(demo_sweep):
    _, config, output, report = demo_sweep
    with (output / "summary.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len(config["coverage_ml"]) * len(config["allocations"])
    for source in report["sources"]:
        for pool_index, allocation in enumerate(config["allocations"]):
            matches = [row for row in rows if int(row["coverage_index"]) == source["index"]
                       and row["reservoir"] == allocation["reservoir"]["name"]]
            assert len(matches) == 1
            row = matches[0]
            assert row["state_id"] == allocation["reservoir"]["state_id"]
            assert float(row["coverage_ml"]) == source["coverage_ml"]
            assert float(row["initial_parent_number"]) == source["initial_parent_number"]
            for name in ("initial_pool_number", "released_pool_number", "remaining_pool_number", "total_wall_energy_debit_j"):
                assert float(row[name]) == source[name][pool_index]
            initial = float(row["initial_pool_number"])
            assert abs(float(row["inventory_balance_residual_number"])) <= 3e-14 * initial


def test_literal_spatial_profile_and_patch_allocations_are_not_normalized(tmp_path):
    path, config = config_file(tmp_path)
    config["coverage_ml"] = [.5]
    profile = np.array([[0., .5, 2.], [1.5, .3, .8]])
    first_fraction = np.array([[.8, .4, .2], [.3, .6, .9]])
    config["parent"]["patch_multiplier"] = profile.tolist()
    config["allocations"][0]["fraction"] = first_fraction.tolist()
    config["allocations"][1]["fraction"] = (1 - first_fraction).tolist()
    path.write_text(json.dumps(config))
    output = tmp_path / "sources"
    process = execute(output, "--demo", "--config", path)
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "report.json").read_bytes())
    history = load_temkin_planar_release_history(output / "source_0000.h5")
    expected_parent = .5 * config["parent"]["monolayer"]["areal_density_m2"] * profile
    expected_pools = expected_parent[None] * np.array([first_fraction, 1 - first_fraction])
    np.testing.assert_allclose(history.coverage_m2[0], expected_pools, rtol=3e-14, atol=0.)
    np.testing.assert_array_equal(history.coverage_m2[:, :, 0, 0], 0.)
    np.testing.assert_array_equal(history.released_m2[:, :, 0, 0], 0.)
    inventory = report["coverage_source_sweep"]["inventory"]
    weighted_mean = np.sum(profile * history.patch_areas_m2) / np.sum(history.patch_areas_m2)
    assert abs(weighted_mean - 1) > .1
    assert inventory["area_weighted_multiplier_mean"] == pytest.approx(weighted_mean, rel=3e-14)
    assert inventory["area_mean_coverage_ml"] == pytest.approx([.5 * weighted_mean], rel=3e-14)
    expected_number = np.sum(expected_parent * history.patch_areas_m2)
    assert report["sources"][0]["initial_parent_number"] == pytest.approx(expected_number, rel=3e-14)


def test_external_surface_copy_and_labels_are_not_promoted_to_physical_heating(tmp_path, demo_sweep):
    _, _, demo, _ = demo_sweep
    path, config = config_file(tmp_path)
    config["geometry"] = None
    path.write_text(json.dumps(config))
    source = tmp_path / "supplied_surface.h5"
    shutil.copyfile(demo / "surface_input.h5", source)
    with h5py.File(source, "r+") as stream:
        attributes = json.loads(stream.attrs["source_attributes_json"])
        attributes["thermal_input_status"] = "owner_supplied"
        attributes["thermal_input_provenance"] = "Manufactured external-input boundary fixture; no actual heating run"
        stream.attrs["source_attributes_json"] = json.dumps(attributes)
    output = tmp_path / "external"
    process = execute(output, "--surface", source, "--config", path)
    assert process.returncode == 0, process.stderr
    report = json.loads((output / "report.json").read_bytes())
    assert report["input_surface_sha256"] == sha(source) == sha(output / "surface_input.h5")
    supplied = load_planar_surface_history(source)
    assert report["thermal_input_attributes"] == plain_metadata(supplied.attributes)
    assert report["thermal_input_attributes"]["source_kind"] == "synthetic"
    assert report["thermal_input_attributes"]["thermal_input_status"] == "owner_supplied"
    assert report["coverage_source_sweep"]["inventory"]["thermal_input_status"] == "owner_supplied"
    for entry in report["sources"]:
        history = load_temkin_planar_release_history(output / entry["source_file"])
        assert plain_metadata(history.surface_attributes) == plain_metadata(supplied.attributes)
    assert report["parameter_status"] == "manufactured"


@pytest.mark.parametrize("existing", ["directory", "file", "symlink"])
def test_existing_destination_is_preserved(tmp_path, existing):
    output = tmp_path / "sources"
    sentinel = b"Preserve completed source evidence"
    if existing == "directory":
        output.mkdir()
        (output / "sentinel").write_bytes(sentinel)
    elif existing == "file":
        output.write_bytes(sentinel)
    else:
        output.symlink_to(tmp_path / "missing-target", target_is_directory=True)
    process = execute(output, "--demo")
    assert process.returncode != 0
    assert "new path" in process.stderr
    if existing == "directory":
        assert list(output.iterdir()) == [output / "sentinel"]
        assert (output / "sentinel").read_bytes() == sentinel
    elif existing == "file":
        assert output.read_bytes() == sentinel
    else:
        assert output.is_symlink()
        assert not (tmp_path / "missing-target").exists()


@pytest.mark.parametrize("arguments", [[], ["--demo", "--surface", "missing.h5"], ["--surface", "missing.h5"]])
def test_source_selection_requires_exactly_one_source_and_external_configuration(tmp_path, arguments):
    output = tmp_path / "sources"
    process = execute(output, *arguments)
    assert process.returncode != 0
    assert not output.exists()


@pytest.mark.parametrize("kind", ["schema", "duplicate", "nonfinite", "unknown_top", "unknown_parent",
                                  "boolean_coverage", "boolean_monolayer", "boolean_fraction", "boolean_geometry",
                                  "boolean_subdivisions", "negative_profile", "wrong_profile_shape",
                                  "nonconserving_fractions", "duplicate_pool", "missing_layer_definition",
                                  "missing_kinetic_reference", "unsupported_coverage", "coupling_range",
                                  "unsorted_knots", "duplicate_knots", "empty_knots", "work_budget",
                                  "parameter_label", "parent_label", "monolayer_label", "allocation_label"])
def test_invalid_or_ambiguous_config_never_creates_partial_sources(tmp_path, kind):
    path, config = config_file(tmp_path)
    if kind == "schema":
        config["schema"] = "unsupported-coverage-source"
    elif kind == "nonfinite":
        config["parent"]["monolayer"]["areal_density_m2"] = float("nan")
    elif kind == "unknown_top":
        config["temperature_k"] = 500.
    elif kind == "unknown_parent":
        config["parent"]["extra_parent_layers"] = 1.
    elif kind == "boolean_coverage":
        config["coverage_ml"] = [0., True]
    elif kind == "boolean_monolayer":
        config["parent"]["monolayer"]["areal_density_m2"] = True
    elif kind == "boolean_fraction":
        config["allocations"][0]["fraction"] = True
    elif kind == "boolean_geometry":
        config["geometry"]["duration_s"] = True
    elif kind == "boolean_subdivisions":
        config["subdivisions_per_interval"] = True
    elif kind == "negative_profile":
        config["parent"]["patch_multiplier"] = -1.
    elif kind == "wrong_profile_shape":
        config["parent"]["patch_multiplier"] = [[1., 1.]]
    elif kind == "nonconserving_fractions":
        config["allocations"][1]["fraction"] = .5
    elif kind == "duplicate_pool":
        config["allocations"][1]["reservoir"]["name"] = config["allocations"][0]["reservoir"]["name"]
    elif kind == "missing_layer_definition":
        del config["parent"]["monolayer"]["definition"]
    elif kind == "missing_kinetic_reference":
        del config["allocations"][0]["reservoir"]["reference_coverage_m2"]
    elif kind == "unsupported_coverage":
        config["coverage_ml"] = [0., 100.]
    elif kind == "coupling_range":
        config["allocations"][0]["reservoir"]["activation_energy_ev"] = 1000.
    elif kind == "unsorted_knots":
        config["coverage_ml"] = [0., .5, .25]
    elif kind == "duplicate_knots":
        config["coverage_ml"] = [0., .25, .25]
    elif kind == "empty_knots":
        config["coverage_ml"] = []
    elif kind == "work_budget":
        config["subdivisions_per_interval"] = 100000
    elif kind == "parameter_label":
        config["parameter_status"] = "sourced"
    elif kind == "parent_label":
        config["parent"]["status"] = "sourced"
    elif kind == "monolayer_label":
        config["parent"]["monolayer"]["status"] = "sourced"
    elif kind == "allocation_label":
        config["allocations"][0]["status"] = "sourced"
    text = json.dumps(config)
    if kind == "duplicate":
        text = text.replace('"schema":', f'"schema": {json.dumps(config["schema"])}, "schema":', 1)
    path.write_text(text)
    output = tmp_path / "sources"
    process = execute(output, "--demo", "--config", path)
    assert process.returncode != 0
    assert not output.exists()
    assert not list(tmp_path.glob(".coverage-*"))


@pytest.mark.parametrize("kind", ["geometry", "external_link", "material_plane"])
def test_external_surface_contract_cannot_be_silently_repaired(tmp_path, demo_sweep, kind):
    _, _, demo, _ = demo_sweep
    path, config = config_file(tmp_path)
    source = tmp_path / "supplied.h5"
    shutil.copyfile(demo / "surface_input.h5", source)
    if kind != "geometry":
        config["geometry"] = None
    if kind == "external_link":
        dependency = tmp_path / "dependency.h5"
        with h5py.File(source, "r+") as stream, h5py.File(dependency, "w") as target:
            stream.copy("temperature_k", target)
            del stream["temperature_k"]
            stream["temperature_k"] = h5py.ExternalLink(str(dependency), "temperature_k")
    elif kind == "material_plane":
        with h5py.File(source, "r+") as stream:
            attributes = json.loads(stream.attrs["source_attributes_json"])
            attributes.update(temperature_location="material_sample_plane", material_sample_depth_m=1e-5,
                              boundary_temperature_reference="")
            stream.attrs["source_attributes_json"] = json.dumps(attributes)
    path.write_text(json.dumps(config))
    output = tmp_path / "sources"
    process = execute(output, "--surface", source, "--config", path)
    assert process.returncode != 0
    assert not output.exists()
    assert not list(tmp_path.glob(".coverage-*"))
