"""MFEM boundary semantics must survive direct construction and saved reload."""
import importlib.util
import json
from pathlib import Path

import h5py
import pytest

from mhdlab.mfem_surface import load_mfem_surface_samples
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.planar_release_io import write_planar_surface_history, load_planar_surface_history
from mhdlab.temkin_planar import _batched_metadata
from mhdlab.temkin_history import _qualified_metadata

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("_mfem_metadata_demo", ROOT / "tools/export_mfem_surface.py")
demo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(demo)


@pytest.fixture
def surface(tmp_path):
    return load_mfem_surface_samples(demo.write_demo(tmp_path / "samples"))[0]


@pytest.mark.parametrize("mutate", [
    lambda a: a.update(temperature_evaluation="volume_cell_mean"),
    lambda a: a.update(source_clock_reference="unknown"),
    lambda a: a.update(initial_temperature_offset_k=-1),
    lambda a: a.update(temperature_value_basis="absolute_kelvin"),
    lambda a: a["mfem_solver_scope"].update(material_motion=True),
    lambda a: a["mfem_face"].update(outward_normal_sign=True),
    lambda a: a.update(source_kind="mfem", mfem_version="4.10", mfem_source_revision="pending"),
    lambda a: a.update(sample_table_sha256=""),
])
def test_direct_constructor_cannot_bypass_mfem_semantics(surface, mutate):
    attributes = surface.attributes
    mutate(attributes)
    with pytest.raises(ValueError):
        PlanarSurfaceHistory(surface.time_s, surface.x_edges_m, surface.z_edges_m,
                             surface.temperature_k, attributes)


def test_reloaded_file_cannot_turn_a_volume_mean_into_material_trace(surface, tmp_path):
    path = tmp_path / "surface.h5"
    write_planar_surface_history(path, surface)
    with h5py.File(path, "r+") as handle:
        attributes = json.loads(handle.attrs["source_attributes_json"])
        attributes["temperature_evaluation"] = "volume_cell_mean"
        handle.attrs["source_attributes_json"] = json.dumps(attributes)
    with pytest.raises(ValueError, match="one-sided trace"):
        load_planar_surface_history(path)


def test_internal_temkin_batch_keeps_and_rechecks_mfem_semantics(surface):
    attributes = _batched_metadata(surface.attributes)
    _qualified_metadata(attributes)
    attributes["temperature_evaluation"] = "volume_cell_mean"
    with pytest.raises(ValueError, match="one-sided trace"):
        _qualified_metadata(attributes)
