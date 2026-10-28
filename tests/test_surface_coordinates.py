"""Face-local frames retain provenance, geometry and thermal launch energies."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from mhdlab.mfem_surface import load_mfem_surface_samples
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.planar_release_io import write_planar_surface_history, load_planar_surface_history
from mhdlab.surface_coordinates import (
    MFEM_SIDE_COORDINATE_SYSTEM, to_physical_coordinates, validate_mfem_face_axes,
)
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_boundary import sample_temkin_planar_emission_packets
from mhdlab.temkin_planar_io import write_temkin_planar_release_history, load_temkin_planar_release_history
from mhdlab.temkin_release import TemkinReservoir

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("_mfem_side_demo", ROOT / "tools/export_mfem_surface.py")
demo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(demo)


@pytest.fixture
def side(tmp_path):
    config_path = demo.write_demo(tmp_path / "side")
    config = json.loads(config_path.read_bytes())
    config["face"].update(coordinate_system=MFEM_SIDE_COORDINATE_SYSTEM,
        physical_in_plane_axis="y", physical_normal_axis="x", normal_position_m=-.0025,
        outward_normal_sign=-1, surface_region="manufactured left thickness face")
    config["x_edges_m"] = [.0025, .0026, .0027]
    config_path.write_text(json.dumps(config))
    return config_path


def test_side_material_trace_reload_to_launch_preserves_physical_axes_and_energy(side, tmp_path):
    surface, _ = load_mfem_surface_samples(side)
    surface_path = tmp_path / "side.h5"
    write_planar_surface_history(surface_path, surface)
    restored = load_planar_surface_history(surface_path)
    assert restored.attributes == surface.attributes
    assert restored.attributes["coordinate_system"] == MFEM_SIDE_COORDINATE_SYSTEM
    # 200 um of thickness times 30 mm current length; no foil-width extrusion.
    np.testing.assert_allclose(restored.patch_areas_m2.sum(), .0002*.03, rtol=2e-15)
    pool = TemkinReservoir(name="authored", state_id="H2O:authored", attempt_frequency_s=1e6,
        activation_energy_ev=0., desorption_energy_ev=.2, internal_energy_ev=0.,
        coverage_coefficient=0., reference_coverage_m2=1e18, max_coverage_m2=1e18,
        kinetic_provenance="Authored constant barrier", energy_provenance="Authored passive removal")
    release = evolve_temkin_planar_release(restored, np.full((1,2,2), 1e17), (pool,),
        parameter_provenance="Authored side-face integration fixture")
    release_path = tmp_path / "side_release.h5"
    write_temkin_planar_release_history(release_path, release)
    source = load_temkin_planar_release_history(release_path)
    draws = np.random.default_rng(20261028).uniform(.0001, .9999, size=(256,7))
    packets = sample_temkin_planar_emission_packets(source, {"H2O:authored":3e-26}, 0, draws,
        face_normal_coordinate_m=-.0025, outward_normal_sign=-1)
    face = source.surface_attributes["mfem_face"]
    position = to_physical_coordinates(packets.resolved_position_m, face)
    velocity = to_physical_coordinates(packets.velocity_m_s, face)
    assert packets.coordinate_system == MFEM_SIDE_COORDINATE_SYSTEM
    assert np.all(position[:,0] == -.0025)
    assert np.all((position[:,1] >= .0025) & (position[:,1] <= .0027))
    assert np.all((position[:,2] >= 0) & (position[:,2] <= .03))
    assert np.all(velocity[:,0] < 0)
    assert np.any(velocity[:,1] < 0) and np.any(velocity[:,1] > 0)
    np.testing.assert_allclose(np.sum(velocity**2, axis=1),
                               np.sum(packets.velocity_m_s**2, axis=1), rtol=2e-15)
    expected_number = 1e17*.0002*.03*(-np.expm1(-1.))
    np.testing.assert_allclose(packets.particle_weight.sum(), expected_number, rtol=2e-15)
    assert np.sum(.5*packets.mass_kg*np.sum(velocity**2, axis=1)*packets.particle_weight) == pytest.approx(
        packets.sampled_translational_energy_j, rel=2e-15)
    # Orthogonal component permutation has not mutated the local packet arrays.
    assert np.all(packets.resolved_position_m[:,1] == -.0025)
    with h5py.File(release_path, "r+") as handle:
        attrs = json.loads(handle.attrs["surface_attributes_json"])
        attrs["mfem_face"]["physical_normal_axis"] = "y"
        handle.attrs["surface_attributes_json"] = json.dumps(attrs)
    with pytest.raises(ValueError, match="explicit physical"):
        load_temkin_planar_release_history(release_path)


@pytest.mark.parametrize("mutation", [
    lambda f: f.pop("physical_normal_axis"),
    lambda f: f.pop("physical_in_plane_axis"),
    lambda f: f.update(physical_normal_axis="y"),
    lambda f: f.update(physical_in_plane_axis="z"),
    lambda f: f.update(coordinate_system="arbitrary rotated face"),
])
def test_side_axis_metadata_cannot_be_omitted_or_conflicted(side, mutation):
    config = json.loads(side.read_bytes())
    mutation(config["face"])
    side.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="axes|coordinate system"):
        load_mfem_surface_samples(side)


def test_side_axes_checked_by_direct_constructor_and_saved_reload(side, tmp_path):
    surface, _ = load_mfem_surface_samples(side)
    attrs = surface.attributes
    attrs["mfem_face"].pop("physical_normal_axis")
    with pytest.raises(ValueError, match="explicit physical"):
        PlanarSurfaceHistory(surface.time_s, surface.x_edges_m, surface.z_edges_m,
                             surface.temperature_k, attrs)
    path = tmp_path / "surface.h5"
    write_planar_surface_history(path, surface)
    with h5py.File(path, "r+") as handle:
        attrs = json.loads(handle.attrs["source_attributes_json"])
        attrs["mfem_face"].pop("physical_in_plane_axis")
        handle.attrs["source_attributes_json"] = json.dumps(attrs)
    with pytest.raises(ValueError, match="explicit physical"):
        load_planar_surface_history(path)


def test_side_coordinate_permission_does_not_relax_other_source_schemas(side):
    surface, _ = load_mfem_surface_samples(side)
    attrs = surface.attributes
    attrs.update(original_source_schema="mhdlab-flash-surface-v2", source_kind="synthetic", flash_version="")
    with pytest.raises(ValueError, match="coordinate_system"):
        PlanarSurfaceHistory(surface.time_s, surface.x_edges_m, surface.z_edges_m,
                             surface.temperature_k, attrs)


def test_broad_faces_keep_identity_but_reject_conflicting_axes(tmp_path):
    surface, _ = load_mfem_surface_samples(demo.write_demo(tmp_path / "broad"))
    face = surface.attributes["mfem_face"]
    values = np.array([[1.,2.,3.],[-5.,0.,8.]])
    assert validate_mfem_face_axes(face) == (0,1,2)
    converted = to_physical_coordinates(values, face)
    np.testing.assert_array_equal(converted, values)
    assert not np.shares_memory(converted, values)
    conflict = deepcopy(face)
    conflict["physical_normal_axis"] = "x"
    with pytest.raises(ValueError, match="conflict"):
        to_physical_coordinates(values, conflict)
    with pytest.raises(ValueError, match="finite"):
        to_physical_coordinates([[float("nan"),0.,0.]], face)
