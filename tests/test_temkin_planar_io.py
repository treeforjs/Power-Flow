"""Law/energy replay, schema corruption and self-contained source persistence."""

import json

import h5py
import numpy as np
import pytest

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory, evolve_planar_release
from mhdlab.planar_release_io import plain_metadata, write_planar_release_history
from mhdlab.temkin_planar import COVERAGE_MODEL, SCHEMA, evolve_temkin_planar_release
from mhdlab.temkin_planar_io import (
    IDENTITIES, PLANAR_NUMERICAL_QUALIFICATIONS, RELEASE_UNITS,
    load_temkin_planar_release_history, write_temkin_planar_release_history,
)
from mhdlab.temkin_release import TemkinReservoir
from mhdlab.thermal_release import BindingReservoir


def sample_surface():
    source = FlashSurfaceHistory([0., 1.], [0., 1.], [[300.], [300.]], "synthetic", "planar-io-test",
                                "manufactured local source", "test boundary",
                                temperature_location="boundary_temperature_proxy", material_sample_depth_m=1e-5,
                                boundary_temperature_reference="manufactured normal-gradient fixture")
    attrs = source.attributes
    attrs["resolved_patch_extraction"] = {"patch_extraction_manifests": [
        {"sampling": {"current_interval_m": interval, "x_edges_m": [0., .2, 1.], "material_sample_depth_m": 1e-5},
         "snapshots": [{"path": "manufactured-local.h5", "sha256": "synthetic-source-label"}]}
        for interval in ([0., .01], [.01, .025])], "other_source_evidence": {"retained": True}}
    return PlanarSurfaceHistory([0., .5e-7, 1e-7], [0., .2, 1.], [0., .01, .025],
                                np.array([[[450., 650.], [550., 900.]]] * 3), attrs)


def sample_release():
    reservoir = TemkinReservoir("manufactured sites", "M:v0", 1e8, .3, .2, .05, .4, 1e18, 1e18,
                                "manufactured independent-pool barrier", "manufactured energy convention")
    return evolve_temkin_planar_release(sample_surface(), np.array([[[1e18, 5e17], [2e17, 8e17]]]),
                                        [reservoir], parameter_provenance="manufactured source test", subdivisions_per_interval=2)


@pytest.mark.parametrize("full_face", [False, True])
def test_round_trip_retains_law_diagnostics_nested_source_and_both_control_choices(tmp_path, full_face):
    original = sample_release()
    target = tmp_path / "source.h5"
    write_temkin_planar_release_history(target, original, include_full_face_uniform=full_face)
    loaded = load_temkin_planar_release_history(target)
    for name in RELEASE_UNITS:
        np.testing.assert_array_equal(getattr(loaded, name), getattr(original, name))
    assert loaded.reservoirs == original.reservoirs
    assert plain_metadata(loaded.surface_attributes) == plain_metadata(original.surface_attributes)
    assert loaded.schema == SCHEMA
    assert loaded.coverage_model == COVERAGE_MODEL
    assert loaded.parameter_provenance == original.parameter_provenance
    with h5py.File(target) as handle:
        assert ("matched_uniform_full_face" in handle) == full_face
        assert handle["inversion_iterations"].dtype.kind == "i"
        assert handle["inversion_regime"].dtype.kind == "i"
        assert handle["patch_areas_m2"].attrs["units"] == "m^2"
    np.testing.assert_array_equal(loaded.extensive_totals("coverage_m2"), original.extensive_totals("coverage_m2"))


@pytest.mark.parametrize("name", list(IDENTITIES) + list(PLANAR_NUMERICAL_QUALIFICATIONS) + [
    "inversion_regime_definitions_json", "reservoirs_json", "surface_attributes_json", "parameter_provenance",
    "full_face_uniform_included",
])
def test_required_schema_law_time_and_numerical_metadata_cannot_be_omitted(tmp_path, name):
    target = tmp_path / "source.h5"
    write_temkin_planar_release_history(target, sample_release())
    with h5py.File(target, "a") as handle:
        del handle.attrs[name]
    with pytest.raises(ValueError):
        load_temkin_planar_release_history(target)


@pytest.mark.parametrize("kind", [
    "law", "time", "coverage_model", "units", "area", "energy", "balance", "width_control", "full_control",
    "reservoir_missing", "reservoir_extra", "kinetics", "domain", "diagnostic", "regime", "numerical_qualification",
    "undeclared_full", "declared_missing_full", "nonbool_full", "nested_depth", "current_average", "patch_interval",
])
def test_tampered_core_source_geometry_kinetics_and_controls_are_rejected(tmp_path, kind):
    target = tmp_path / "source.h5"
    write_temkin_planar_release_history(target, sample_release(), include_full_face_uniform=True)
    with h5py.File(target, "a") as handle:
        if kind in ("law", "time", "coverage_model"):
            handle.attrs[{"law": "law_identity", "time": "time_model", "coverage_model": "coverage_model"}[kind]] = "unknown"
        elif kind == "units":
            handle["z_edges_m"].attrs["units"] = "cm"
        elif kind == "area":
            handle["patch_areas_m2"][0, 0] *= 2
        elif kind == "energy":
            handle["translational_energy_j_m2"][0, 0, 0, 0] *= 2
        elif kind == "balance":
            handle["coverage_m2"][1, 0, 0, 0] = 0
        elif kind in ("width_control", "full_control"):
            group = "matched_uniform_width" if kind == "width_control" else "matched_uniform_full_face"
            handle[f"{group}/released_m2"][0, 0, 0, 0] *= 2
        elif kind.startswith("reservoir") or kind in ("kinetics", "domain"):
            records = json.loads(handle.attrs["reservoirs_json"])
            if kind == "reservoir_missing":
                del records[0]["energy_provenance"]
            elif kind == "reservoir_extra":
                records[0]["guess"] = 1
            elif kind == "domain":
                records[0]["max_coverage_m2"] = 9e17
            else:
                records[0]["attempt_frequency_s"] *= 2
            handle.attrs["reservoirs_json"] = json.dumps(records)
        elif kind in ("diagnostic", "regime"):
            handle["inversion_iterations" if kind == "diagnostic" else "inversion_regime"][0, 0, 0, 0] = 129 if kind == "diagnostic" else 0
        elif kind == "numerical_qualification":
            handle.attrs["inversion_log_residual_tolerance"] = 1.
        elif kind == "undeclared_full":
            handle.attrs["full_face_uniform_included"] = False
        elif kind == "declared_missing_full":
            del handle["matched_uniform_full_face"]
        elif kind == "nonbool_full":
            handle.attrs["full_face_uniform_included"] = 1
        else:
            metadata = json.loads(handle.attrs["surface_attributes_json"])
            if kind == "nested_depth":
                metadata["resolved_patch_extraction"]["patch_extraction_manifests"][0]["sampling"]["material_sample_depth_m"] = 2e-5
            elif kind == "current_average":
                metadata["current_averaging_interval_m_json"] = "[0, 0.025]"
            else:
                metadata["current_patch_intervals_m"] = [[0., .025]]
            handle.attrs["surface_attributes_json"] = json.dumps(metadata)
    with pytest.raises(ValueError):
        load_temkin_planar_release_history(target)


@pytest.mark.parametrize("storage", ["external_link", "soft_link", "virtual", "external_storage"])
@pytest.mark.parametrize("dataset", ["temperature_k", "matched_uniform_width/released_m2", "patch_areas_m2"])
def test_no_core_or_derived_dataset_can_escape_the_parent_file_hash(tmp_path, storage, dataset):
    target = tmp_path / "source.h5"
    donor = tmp_path / "donor.h5"
    write_temkin_planar_release_history(target, sample_release())
    with h5py.File(target, "r") as handle:
        values, units = handle[dataset][...], handle[dataset].attrs["units"]
    with h5py.File(donor, "w") as handle:
        handle.create_dataset("values", data=values).attrs["units"] = units
    with h5py.File(target, "a") as handle:
        del handle[dataset]
        if storage == "external_link":
            handle[dataset] = h5py.ExternalLink(str(donor), "/values")
        elif storage == "soft_link":
            handle.create_dataset("alternative_values", data=values).attrs["units"] = units
            handle[dataset] = h5py.SoftLink("/alternative_values")
        elif storage == "virtual":
            layout = h5py.VirtualLayout(values.shape, dtype=values.dtype)
            layout[:] = h5py.VirtualSource(str(donor), "values", values.shape)
            handle.create_virtual_dataset(dataset, layout).attrs["units"] = units
        else:
            item = handle.create_dataset(dataset, values.shape, dtype=values.dtype,
                                         external=[(str(tmp_path / "external.bin"), 0, h5py.h5f.UNLIMITED)])
            item[...] = values
            item.attrs["units"] = units
    with pytest.raises(ValueError, match="hard links|self-contained|nonvirtual|external storage"):
        load_temkin_planar_release_history(target)


@pytest.mark.parametrize("metadata", ['{"key": NaN}', '{"schema": "x", "schema": "y"}'])
def test_json_metadata_must_be_finite_and_unambiguous(tmp_path, metadata):
    target = tmp_path / "source.h5"
    write_temkin_planar_release_history(target, sample_release())
    with h5py.File(target, "a") as handle:
        handle.attrs["surface_attributes_json"] = metadata
    with pytest.raises(ValueError, match="nonfinite|duplicate"):
        load_temkin_planar_release_history(target)


@pytest.mark.parametrize("dataset", ["time_edges_s", "coverage_m2", "temperature_k", "matched_uniform_width/released_m2", "patch_areas_m2"])
def test_inconsistent_huge_dataset_shape_is_rejected_before_materialization(tmp_path, monkeypatch, dataset):
    target = tmp_path / "source.h5"
    write_temkin_planar_release_history(target, sample_release())
    with h5py.File(target, "a") as handle:
        original_shape = handle[dataset].shape
        units = handle[dataset].attrs["units"]
        shape = (*original_shape[:-1], 1_000_000_000)
        del handle[dataset]
        handle.create_dataset(dataset, shape=shape, chunks=tuple(1 for _ in shape), dtype="f8").attrs["units"] = units
    original_getitem = h5py.Dataset.__getitem__
    def guarded_getitem(self, key):
        if self.name == f"/{dataset}":
            raise AssertionError("huge malformed dataset must not be read")
        return original_getitem(self, key)
    monkeypatch.setattr(h5py.Dataset, "__getitem__", guarded_getitem)
    with pytest.raises(ValueError, match="memory limit|shape disagrees"):
        load_temkin_planar_release_history(target)


def test_baseline_schema_and_wrong_records_are_not_accepted_as_temkin(tmp_path):
    baseline = evolve_planar_release(sample_surface(), np.ones((1, 2, 2)),
                                     [BindingReservoir("baseline", "M:v0", 1., .1, .2)],
                                     parameter_provenance="manufactured constant-barrier test")
    target = tmp_path / "baseline.h5"
    write_planar_release_history(target, baseline)
    with pytest.raises(ValueError, match="schema"):
        load_temkin_planar_release_history(target)
    with pytest.raises(ValueError, match="TemkinPlanarReleaseHistory"):
        write_temkin_planar_release_history(target, baseline)
    with pytest.raises(ValueError, match="boolean"):
        write_temkin_planar_release_history(tmp_path / "unused.h5", sample_release(), include_full_face_uniform=1)
    with pytest.raises(ValueError):
        load_temkin_planar_release_history(tmp_path / "missing.h5")


def test_resolved_only_archive_retains_local_law_when_uniform_control_is_unrepresentable(tmp_path):
    base = FlashSurfaceHistory([0., 1.], [0., 1.], [[300.], [300.]], 'synthetic', 'local-archive-range',
                               'manufactured floating-point range check', 'test boundary')
    surface = PlanarSurfaceHistory([0., 1.], [0., 1e-200, 1e200], [0., 1.],
                                   np.full((2,1,2),300.), base.attributes)
    pool = TemkinReservoir('local sites','M:v0',1.,0.,0.,0.,0.,1e18,1e18,
                           'manufactured independent pool','manufactured zero binding energy')
    original = evolve_temkin_planar_release(surface, [[[1.,0.]]], [pool], parameter_provenance='manufactured range check')
    with pytest.raises(ValueError, match='below the qualified normal range'):
        write_temkin_planar_release_history(tmp_path/'default.h5',original)
    target=tmp_path/'resolved.h5'
    write_temkin_planar_release_history(target,original,include_width_uniform=False)
    loaded=load_temkin_planar_release_history(target)
    for name in RELEASE_UNITS:
        np.testing.assert_array_equal(getattr(loaded,name),getattr(original,name))
    with pytest.raises(ValueError,match='below the qualified normal range'):
        loaded.matched_uniform_width()
    with h5py.File(target,'a') as f:
        assert not f.attrs['width_uniform_included']
        assert 'matched_uniform_width' not in f
        f['released_m2'][0,0,0,0] *= .5
    with pytest.raises(ValueError):
        load_temkin_planar_release_history(target)


def test_resolved_only_archive_rejects_hidden_control_and_nonboolean_declaration(tmp_path):
    target=tmp_path/'source.h5';write_temkin_planar_release_history(target,sample_release())
    with h5py.File(target,'a') as f: f.attrs['width_uniform_included']=False
    with pytest.raises(ValueError,match='undeclared width-uniform'):
        load_temkin_planar_release_history(target)
    with h5py.File(target,'a') as f: f.attrs['width_uniform_included']='false'
    with pytest.raises(ValueError,match='must be a boolean'):
        load_temkin_planar_release_history(target)
