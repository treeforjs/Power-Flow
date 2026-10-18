"""Exercise the external-file CLI with explicitly manufactured inputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pytest

from mhdlab.flash_surface import SCHEMA, load_surface_history


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "export_flash_surface.py"


def invoke(*args, cwd=None):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], cwd=cwd,
                          capture_output=True, text=True, check=False)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_demo_uses_configured_export_and_preserves_plane_labels_units_and_hashes(tmp_path):
    output = tmp_path / "demo"
    result = invoke("--demo", "--output-dir", output, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    summary = json.loads(result.stdout)
    manifest_path = output / "extraction_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    configuration_path = output / "demo_config.json"
    configuration = json.loads(configuration_path.read_text())
    history = load_surface_history(output / "surface.h5")
    assert summary["source_kind"] == manifest["source_kind"] == history.source_kind == "synthetic"
    assert history.temperature_location == "material_sample_plane"
    assert history.material_sample_depth_m == 5e-6
    assert history.boundary_temperature_reference == ""
    assert history.temperature_k.shape == (3, 6)
    np.testing.assert_allclose(history.time_s, [0, 60e-9, 120e-9])
    np.testing.assert_allclose(history.temperature_k[0], 300)
    assert history.temperature_k[-1, 0] > history.temperature_k[-1, 2]
    assert configuration["source_metadata"]["source_kind"] == "synthetic"
    assert configuration["configuration_status"] == "manufactured"
    assert manifest["configuration"]["sha256"] == digest(configuration_path)
    assert manifest["surface_history"]["sha256"] == digest(output / "surface.h5")
    assert manifest["surface_history_schema"] == SCHEMA
    for snapshot in manifest["snapshots"]:
        assert snapshot["sha256"] == digest(snapshot["path"])
        assert Path(snapshot["path"]).is_absolute()
        assert snapshot["time_unit"] == "ns"
        assert snapshot["time_scale_to_s"] == 1e-9
    with h5py.File(output / "surface.h5", "r") as handle:
        assert handle.attrs["extraction_config_path"] == str(configuration_path.resolve())
        assert handle.attrs["extraction_config_sha256"] == digest(configuration_path)
        assert handle.attrs["extraction_manifest_path"] == str(manifest_path.resolve())
        assert handle.attrs["temperature_location"] == "material_sample_plane"
        assert handle.attrs["material_sample_depth_m"] == 5e-6
        assert handle["time_s"].attrs["units"] == "s"
        assert handle["x_edges_m"].attrs["units"] == "m"
        assert handle["temperature_k"].attrs["units"] == "K"
    # Both paths must feed exactly the same loader/projector and temperatures.
    second = tmp_path / "configured-copy"
    result = invoke("--config", configuration_path, "--output-dir", second, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    reread = load_surface_history(second / "surface.h5")
    np.testing.assert_array_equal(reread.temperature_k, history.temperature_k)
    assert json.loads((second / "extraction_manifest.json").read_text())["snapshots"] == manifest["snapshots"]


def external_configuration(tmp_path):
    directory = tmp_path / "external-config"
    raw = directory / "raw"
    raw.mkdir(parents=True)
    for index in range(2):
        with h5py.File(raw / f"snapshot_{index}.h5", "w") as handle:
            handle.create_dataset("boxes_cm", data=[[[0, 1], [-.02, 0]]]).attrs["units"] = "cm"
            handle["node_code"] = np.asarray([7], dtype=np.int32)
            handle.create_dataset("material_C", data=[[[20 + index, 40 + index]]]).attrs["units"] = "degC"
            handle.create_dataset("fraction", data=[[[1., 1.]]]).attrs["units"] = "1"
    config = {
        "schema": "mhdlab-flash-extract-config-v1", "configuration_status": "configured",
        "layout": {"bounding_box_dataset": "boxes_cm", "node_type_dataset": "node_code",
                   "temperature_dataset": "material_C", "material_fraction_dataset": "fraction",
                   "bounding_box_axes": ["x", "y"], "scalar_axes": ["block", "y", "x"],
                   "physical_axes": {"x": "width", "y": "normal"}, "leaf_node_value": 7,
                   "nonleaf_node_values": [3], "coordinate_unit": "cm", "coordinate_scale_to_m": .01,
                   "temperature_unit": "degC", "temperature_scale_to_k": 1, "temperature_offset_k": 273.15,
                   "min_material_fraction": .99999999},
        "face": {"x_edges_m": [0, .005, .01], "face_normal_coordinate_m": 0,
                 "outward_sign": 1, "material_sample_depth_m": 50e-6,
                 "current_interval_m": [0, .02]},
        "snapshots": [{"path": f"raw/snapshot_{i}.h5", "time_value": 10 + i,
                       "time_unit": "ns", "time_scale_to_s": 1e-9} for i in range(2)],
        "source_metadata": {"source_kind": "synthetic", "run_id": "external-manufactured-fixture",
                            "source_reference": "Analytic Celsius/centimetre fixture, not a simulation run",
                            "surface_region": "declared test plane"},
    }
    path = directory / "config.json"
    path.write_text(json.dumps(config), encoding="utf-8")
    return path, config


def test_external_relative_paths_resolve_against_config_and_explicit_units_convert(tmp_path):
    path, config = external_configuration(tmp_path)
    output = tmp_path / "elsewhere" / "export"
    result = invoke("--config", path, "--output-dir", output, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    history = load_surface_history(output / "surface.h5")
    np.testing.assert_allclose(history.temperature_k, [[293.15, 313.15], [294.15, 314.15]])
    np.testing.assert_allclose(history.time_s, [10e-9, 11e-9])
    manifest = json.loads((output / "extraction_manifest.json").read_text())
    assert manifest["layout"]["temperature_offset_k"] == 273.15
    assert manifest["layout"]["coordinate_scale_to_m"] == .01
    for snapshot in manifest["snapshots"]:
        assert Path(snapshot["path"]).parent == path.parent / "raw"
    assert config["snapshots"][0]["path"] == "raw/snapshot_0.h5"


@pytest.mark.parametrize("mutation,message", [
    (lambda c: c.update(unknown_key=True), "unknown keys"),
    (lambda c: c["layout"].update(temperature_field="ignored"), "invalid layout"),
    (lambda c: c["face"].update(material_sample_depth_m=0), "positive"),
    (lambda c: c["layout"].pop("min_material_fraction"), "min_material_fraction"),
    (lambda c: c["snapshots"][0].update(unknown_time=0), "invalid snapshot"),
    (lambda c: c["snapshots"][1].update(time_value=10), "strictly increasing"),
    (lambda c: c["source_metadata"].update(inferred_source="flash"), "unknown keys"),
    (lambda c: c.update(configuration_status="template"), "template configuration cannot execute"),
    (lambda c: c.update(configuration_status="manufactured") or c["source_metadata"].update(source_kind="flash"), "must declare"),
    (lambda c: c.update(schema="unknown-v1"), "schema must equal"),
])
def test_config_mistakes_reject_without_canonical_outputs(tmp_path, mutation, message):
    path, config = external_configuration(tmp_path)
    mutation(config)
    path.write_text(json.dumps(config), encoding="utf-8")
    output = tmp_path / "rejected"
    result = invoke("--config", path, "--output-dir", output)
    assert result.returncode != 0
    assert message in result.stderr
    assert not (output / "surface.h5").exists()
    assert not (output / "extraction_manifest.json").exists()


@pytest.mark.parametrize("text,message", [
    ('{"schema": "one", "schema": "two"}', "duplicate JSON key"),
    ('{"value": NaN}', "nonfinite JSON constant"),
])
def test_ambiguous_or_nonfinite_json_rejects(tmp_path, text, message):
    config = tmp_path / "invalid.json"
    config.write_text(text, encoding="utf-8")
    result = invoke("--config", config, "--output-dir", tmp_path / "out")
    assert result.returncode != 0
    assert message in result.stderr


def test_shipped_template_is_explicitly_non_executable(tmp_path):
    result = invoke("--config", ROOT / "configs" / "flash_extract_example.json", "--output-dir", tmp_path)
    assert result.returncode != 0
    assert "template configuration cannot execute" in result.stderr


def test_cli_requires_exactly_one_input_mode(tmp_path):
    assert invoke("--output-dir", tmp_path).returncode != 0
    assert invoke("--demo", "--config", "unused.json", "--output-dir", tmp_path).returncode != 0


def test_default_output_directory_is_relative_to_calling_directory(tmp_path):
    result = invoke("--demo", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "runs" / "flash_extract_reference" / "surface.h5").is_file()


def test_output_cannot_replace_a_snapshot(tmp_path):
    path, config = external_configuration(tmp_path)
    raw_dir = path.parent / "raw"
    original = raw_dir / "snapshot_0.h5"
    protected = raw_dir / "surface.h5"
    original.rename(protected)
    config["snapshots"][0]["path"] = "raw/surface.h5"
    path.write_text(json.dumps(config), encoding="utf-8")
    before = digest(protected)
    result = invoke("--config", path, "--output-dir", raw_dir)
    assert result.returncode != 0
    assert "must not overwrite" in result.stderr
    assert digest(protected) == before
