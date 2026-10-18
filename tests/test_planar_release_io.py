"""SI round trips and independent corruption checks for planar artifacts."""

import json

import h5py
import numpy as np
import pytest

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory, evolve_planar_release
from mhdlab.planar_release_io import (
    RELEASE_UNITS, load_planar_release_history, load_planar_surface_history, plain_metadata,
    write_planar_release_history, write_planar_surface_history,
)
from mhdlab.thermal_release import BindingReservoir


def sample_surface(**overrides):
    base = FlashSurfaceHistory([0., 1.], [0., 1.], [[300.], [400.]], "synthetic", "io-test",
                               "manufactured patch persistence check", "test face", **overrides)
    metadata = base.attributes
    metadata["extraction"] = {"patches": [{"interval_m": [0., .01], "hash": "test-source-label"}]}
    return PlanarSurfaceHistory([0., .5, 1.], [0., .2, 1.], [0., .01, .025],
                                np.array([[[300., 600.], [900., 1200.]]] * 3), metadata)


def sample_release():
    return evolve_planar_release(sample_surface(), np.full((1, 2, 2), 1e18),
                                  (BindingReservoir("specified-pool", "M:state", 1., .1, .2, .3),),
                                  parameter_provenance="manufactured source law")


def test_surface_round_trip_preserves_planes_unknown_metadata_and_patch_geometry(tmp_path):
    source = sample_surface(temperature_location="material_sample_plane", material_sample_depth_m=12.5e-6)
    target = tmp_path / "surface.h5"
    write_planar_surface_history(target, source)
    loaded = load_planar_surface_history(target)
    for name in ("time_s", "x_edges_m", "z_edges_m", "temperature_k", "patch_areas_m2"):
        np.testing.assert_array_equal(getattr(loaded, name), getattr(source, name))
    assert loaded.attributes == source.attributes
    assert loaded.source_attributes["temperature_location"] == "material_sample_plane"
    assert loaded.source_attributes["material_sample_depth_m"] == 12.5e-6


@pytest.mark.parametrize("full_face", [False, True])
def test_release_round_trip_revalidates_local_ledgers_and_distinct_controls(tmp_path, full_face):
    original = sample_release()
    target = tmp_path / "release.h5"
    write_planar_release_history(target, original, include_full_face_uniform=full_face)
    loaded = load_planar_release_history(target)
    for name in RELEASE_UNITS:
        np.testing.assert_array_equal(getattr(loaded, name), getattr(original, name))
    assert loaded.reservoirs == original.reservoirs
    assert plain_metadata(loaded.surface_attributes) == plain_metadata(original.surface_attributes)
    with h5py.File(target) as handle:
        assert ("matched_uniform_full_face" in handle) == full_face
        assert "matched_uniform_width/released_m2" in handle
        assert handle["patch_areas_m2"].attrs["units"] == "m^2"


@pytest.mark.parametrize("kind", ["units", "area", "temperature_kind", "current_average", "schema"])
def test_surface_metadata_or_geometry_corruption_is_rejected(tmp_path, kind):
    target = tmp_path / "surface.h5"
    write_planar_surface_history(target, sample_surface())
    with h5py.File(target, "a") as handle:
        if kind == "units":
            handle["z_edges_m"].attrs["units"] = "cm"
        elif kind == "area":
            handle["patch_areas_m2"][0, 0] *= 2
        elif kind == "schema":
            handle.attrs["schema"] = "unknown"
        else:
            metadata = json.loads(handle.attrs["source_attributes_json"])
            if kind == "temperature_kind":
                del metadata["temperature_location"]
            else:
                metadata["current_averaging_interval_m_json"] = "[0, 0.025]"
            handle.attrs["source_attributes_json"] = json.dumps(metadata)
    with pytest.raises(ValueError):
        load_planar_surface_history(target)


@pytest.mark.parametrize("kind", ["energy", "balance", "width_control", "full_control", "reservoir", "undeclared_full"])
def test_release_core_and_derived_corruption_is_rejected(tmp_path, kind):
    target = tmp_path / "release.h5"
    write_planar_release_history(target, sample_release(), include_full_face_uniform=True)
    with h5py.File(target, "a") as handle:
        if kind == "energy":
            handle["translational_energy_j_m2"][0, 0, 0, 0] *= 2
        elif kind == "balance":
            handle["coverage_m2"][1, 0, 0, 0] = 0
        elif kind == "width_control":
            handle["matched_uniform_width/released_m2"][0, 0, 0, 0] *= 2
        elif kind == "full_control":
            handle["matched_uniform_full_face/internal_energy_j_m2"][0, 0, 0, 0] *= 2
        elif kind == "reservoir":
            records = json.loads(handle.attrs["reservoirs_json"])
            del records[0]["internal_energy_ev"]
            handle.attrs["reservoirs_json"] = json.dumps(records)
        else:
            handle.attrs["full_face_uniform_included"] = False
    with pytest.raises(ValueError):
        load_planar_release_history(target)


@pytest.mark.parametrize("storage", ["external_link", "virtual"])
@pytest.mark.parametrize("target_kind", ["surface_temperature", "release_derived"])
def test_arrays_and_derived_controls_cannot_escape_parent_file_hash(tmp_path, storage, target_kind):
    target = tmp_path / "target.h5"
    donor = tmp_path / "donor.h5"
    if target_kind == "surface_temperature":
        write_planar_surface_history(target, sample_surface())
        dataset, units = "temperature_k", "K"
        loader = load_planar_surface_history
    else:
        write_planar_release_history(target, sample_release())
        dataset, units = "matched_uniform_width/released_m2", "m^-2"
        loader = load_planar_release_history
    with h5py.File(target, "r") as handle:
        values = handle[dataset][...]
    with h5py.File(donor, "w") as handle:
        handle.create_dataset("values", data=values).attrs["units"] = units
    with h5py.File(target, "a") as handle:
        del handle[dataset]
        if storage == "external_link":
            handle[dataset] = h5py.ExternalLink(str(donor), "/values")
        else:
            virtual = h5py.VirtualLayout(values.shape, dtype=values.dtype)
            virtual[:] = h5py.VirtualSource(str(donor), "values", values.shape)
            handle.create_virtual_dataset(dataset, virtual).attrs["units"] = units
    with pytest.raises(ValueError, match="hard links|self-contained|nonvirtual"):
        loader(target)


def test_missing_file_and_nonfinite_json_are_rejected(tmp_path):
    with pytest.raises(ValueError):
        load_planar_surface_history(tmp_path / "missing.h5")
    target = tmp_path / "surface.h5"
    write_planar_surface_history(target, sample_surface())
    with h5py.File(target, "a") as handle:
        handle.attrs["source_attributes_json"] = '{"value": NaN}'
    with pytest.raises(ValueError, match="nonfinite"):
        load_planar_surface_history(target)
