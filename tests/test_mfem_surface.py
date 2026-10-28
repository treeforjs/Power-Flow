"""Provenance-safe MFEM table/planar-release integration and unit guards."""
from copy import deepcopy
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pytest

from mhdlab.mfem_surface import export_mfem_surface, load_mfem_surface_samples
from mhdlab.planar_release import MFEM_SURFACE_SCHEMA, PlanarSurfaceHistory, evolve_planar_release
from mhdlab.planar_release_io import load_planar_surface_history
from mhdlab.thermal_release import BindingReservoir
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_io import write_temkin_planar_release_history, load_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("_mfem_surface_export_test", ROOT / "tools/export_mfem_surface.py")
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)


@pytest.fixture
def example(tmp_path):
    return tool.write_demo(tmp_path / "inputs")


def mutate(path, change):
    config = json.loads(path.read_bytes())
    change(config)
    path.write_text(json.dumps(config), encoding="utf-8")


def alter_rows(path, change):
    config = json.loads(path.read_bytes())
    table = path.parent / config["table"]["path"]
    with table.open(newline="") as stream:
        reader = csv.DictReader(stream)
        columns, rows = reader.fieldnames, list(reader)
    change(rows)
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    mutate(path, lambda c: c["table"].update(sha256=hashlib.sha256(table.read_bytes()).hexdigest()))


def test_increment_conversion_and_local_patches_roundtrip_without_flash_label(example, tmp_path):
    # Authored delta-T sample: at time 1, (z=1,x=1) is140K above300K.
    history, manifest = load_mfem_surface_samples(example)
    assert history.temperature_k[1, 1, 1] == 440
    assert history.temperature_k[0, 0, 0] == 300
    np.testing.assert_allclose(history.patch_areas_m2, [[5e-6, 15e-6], [10e-6, 30e-6]])
    assert history.attributes["source_kind"] == "synthetic"
    assert history.attributes["original_source_schema"] == MFEM_SURFACE_SCHEMA
    assert "flash_version" not in history.attributes
    result = export_mfem_surface(example, tmp_path / "export")
    loaded = load_planar_surface_history(result["surface_history"]["path"])
    np.testing.assert_array_equal(loaded.temperature_k, history.temperature_k)
    assert loaded.attributes == history.attributes
    assert result["temperature_offset_applied_once_k"] == 300
    with h5py.File(result["surface_history"]["path"]) as stream:
        assert stream["temperature_k"].attrs["units"] == "K"
    assert manifest["native_mfem_gridfunction_parsed"] is False


def test_configured_mfem_provenance_survives_existing_finite_inventory_release(example):
    mutate(example, lambda c: (c.update(configuration_status="configured"),
        c["source"].update(source_kind="mfem", mfem_version="4.10", mfem_source_revision="a" * 40)))
    history, _ = load_mfem_surface_samples(example)
    pool = BindingReservoir("specified", "H2O:test", 1e6, 0., 0.2)
    coverage = np.full((1, 2, 2), 1e18)
    release = evolve_planar_release(history, coverage, (pool,), parameter_provenance="manufactured constant hazard")
    assert release.surface_attributes["source_kind"] == "mfem"
    assert release.surface_attributes["mfem_version"] == "4.10"
    initial = release.extensive_totals("coverage_m2")[0]
    remaining = release.extensive_totals("coverage_m2")[-1]
    released = release.extensive_totals().sum(axis=0)
    np.testing.assert_allclose(released, initial * (1 - np.exp(-2.)), rtol=3e-15)
    np.testing.assert_allclose(initial - remaining, released, rtol=3e-15)
    assert release.surface_attributes["mfem_solver_scope"]["hall_induction"] is False


def test_configured_mfem_version_cannot_remain_a_placeholder(example):
    mutate(example, lambda c: (c.update(configuration_status="configured"),
        c["source"].update(source_kind="mfem", mfem_version="unknown", mfem_source_revision="a" * 40)))
    with pytest.raises(ValueError, match="placeholder"):
        load_mfem_surface_samples(example)


def test_material_plane_is_diagnostic_until_explicit_boundary_proxy(example):
    mutate(example, lambda c: c["temperature"].update(location="material_sample_plane",
        evaluation="resolved_material_sample_plane_patch_average", material_sample_depth_m=1e-6))
    history, _ = load_mfem_surface_samples(example)
    pool = BindingReservoir("specified", "H2O:test", 1., 0., .2)
    with pytest.raises(ValueError, match="boundary"):
        evolve_planar_release(history, np.ones((1, 2, 2)), (pool,), parameter_provenance="manufactured")


def test_mfem_temkin_planar_release_retains_backend_through_internal_batch_and_io(example, tmp_path):
    mutate(example, lambda c: (c.update(configuration_status="configured"),
        c["source"].update(source_kind="mfem", mfem_version="4.10", mfem_source_revision="b" * 40)))
    history, _ = load_mfem_surface_samples(example)
    raw = json.loads((ROOT / "configs/temkin_manufactured.json").read_bytes())["reservoirs"][0]
    initial = raw.pop("initial_coverage_m2")
    pool = TemkinReservoir(**raw)
    result = evolve_temkin_planar_release(history, np.full((1, 2, 2), initial), (pool,),
                                         parameter_provenance="Manufactured MFEM metadata compatibility fixture")
    target = tmp_path / "temkin.h5"
    write_temkin_planar_release_history(target, result)
    restored = load_temkin_planar_release_history(target)
    assert restored.surface_attributes["source_kind"] == "mfem"
    assert restored.surface_attributes["original_source_schema"] == MFEM_SURFACE_SCHEMA
    assert "flash_version" not in restored.surface_attributes
    np.testing.assert_array_equal(restored.released_m2, result.released_m2)
    np.testing.assert_allclose(restored.extensive_totals("coverage_m2")[0]
        - restored.extensive_totals("coverage_m2")[-1], restored.extensive_totals().sum(axis=0), rtol=3e-12)


@pytest.mark.parametrize("change,match", [
    (lambda c: c.update(configuration_status="template"), "template"),
    (lambda c: c["source"].update(source_kind="mfem"), "separate"),
    (lambda c: c["source"]["si_unit_mapping"].update(time="unknown"), "SI units"),
    (lambda c: c["temperature"].update(initial_temperature_k=0), "positive initial"),
    (lambda c: c["temperature"].update(value_basis="absolute_kelvin"), "offset"),
    (lambda c: c["temperature"].update(evaluation="volume_cell_mean"), "one-sided"),
    (lambda c: c["temperature"].update(material_sample_depth_m=1e-6), "zero sampling"),
    (lambda c: c["face"].update(current_direction_representation="whole_current_averaged"), "inflated"),
    (lambda c: c["face"].update(outward_normal_sign=0), "normal_sign"),
    (lambda c: c["source"]["solver_scope"].update(material_motion=True), "moving material"),
    (lambda c: c["source"].update(unit_mapping_reference="REPLACE_WITH_MAPPING"), "placeholder"),
    (lambda c: c["temperature"].update(zero_convention_reference="unknown"), "placeholder"),
    (lambda c: c["source"].update(heating_backend="flash"), "Joule backend"),
])
def test_ambiguous_units_temperature_and_geometry_reject(example, change, match):
    mutate(example, change)
    with pytest.raises(ValueError, match=match):
        load_mfem_surface_samples(example)


@pytest.mark.parametrize("change,match", [
    (lambda rows: rows.append(dict(rows[0])), "duplicate"),
    (lambda rows: rows.pop(), "cover every"),
    (lambda rows: rows[0].update(time_s="-1"), "timestamp"),
    (lambda rows: rows[0].update(x_index="2"), "outside"),
    (lambda rows: rows[0].update(x_index="0.0"), "integers"),
    (lambda rows: rows[0].update(temperature_value="NaN"), "finite and positive"),
])
def test_table_completeness_and_finiteness_reject(example, change, match):
    alter_rows(example, change)
    with pytest.raises(ValueError, match=match):
        load_mfem_surface_samples(example)


def test_changed_sample_bytes_cannot_pass_original_hash(example):
    config = json.loads(example.read_bytes())
    table = example.parent / config["table"]["path"]
    table.write_bytes(table.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA256"):
        load_mfem_surface_samples(example)


def test_absolute_kelvin_table_receives_no_double_temperature_offset(example):
    alter_rows(example, lambda rows: [row.update(temperature_value=str(float(row["temperature_value"]) + 300)) for row in rows])
    mutate(example, lambda c: c["temperature"].update(value_basis="absolute_kelvin", initial_temperature_k=0))
    history, manifest = load_mfem_surface_samples(example)
    assert history.temperature_k[1, 1, 1] == 440
    assert manifest["temperature_offset_applied_once_k"] == 0


def test_mfem_metadata_cannot_be_recast_into_flash_schema(example):
    history, _ = load_mfem_surface_samples(example)
    attrs = deepcopy(history.attributes)
    attrs.update(source_kind="mfem", original_source_schema="mhdlab-flash-surface-v2", flash_version="4.8")
    with pytest.raises(ValueError, match="MFEM data"):
        PlanarSurfaceHistory(history.time_s, history.x_edges_m, history.z_edges_m, history.temperature_k, attrs)


def test_cli_demo_is_explicitly_synthetic_and_rejects_overwrite(tmp_path):
    args = [sys.executable, str(ROOT / "tools/export_mfem_surface.py"), "--demo", "--output-dir", str(tmp_path / "demo")]
    completed = subprocess.run(args, text=True, capture_output=True, check=True)
    report = json.loads(completed.stdout)
    assert report["source_kind"] == "synthetic"
    assert report["new_heating_simulation_performed"] is False
    assert Path(report["surface_history"]["path"]).is_file()
    assert subprocess.run(args, text=True, capture_output=True).returncode == 2
