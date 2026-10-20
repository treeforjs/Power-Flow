"""Temkin law identity, canonical storage and tamper-resistant interval checks."""

import json

import h5py
import numpy as np
import pytest

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.surface_release_history import evolve_surface_release
from mhdlab.surface_release_io import load_release_history, write_release_history
from mhdlab.temkin_history import TemkinReleaseHistory, evolve_temkin_release
from mhdlab.temkin_io import SCHEMA, NUMERICAL_QUALIFICATIONS, load_temkin_history, write_temkin_history
from mhdlab.temkin_release import LAW_IDENTITY, TemkinReservoir
from mhdlab.thermal_release import BindingReservoir, SurfaceInventory


def source():
    surface = FlashSurfaceHistory([0, 1e-8, 3e-8], [0, .2, 1], [[450, 650], [650, 750], [800, 950]],
                                   "synthetic", "io-test", "manufactured Temkin temperature", "test face")
    reservoir = TemkinReservoir("manufactured", "M:v0", 1e8, .3, .2, .05, .4, 1e18, 1e18,
                               "manufactured kinetic parameters", "manufactured constant energy cost")
    return evolve_temkin_release(surface, SurfaceInventory([[1e18, 5e17]]), [reservoir],
                                 parameter_provenance="manufactured IO fixture", subdivisions_per_interval=2)


def test_round_trip_preserves_law_provenance_diagnostics_and_moments(tmp_path):
    path = tmp_path / "temkin.h5"
    original = source()
    write_temkin_history(path, original)
    actual = load_temkin_history(path)
    assert isinstance(actual, TemkinReleaseHistory)
    assert actual.law_identity == original.law_identity == LAW_IDENTITY
    assert actual.reservoirs == original.reservoirs
    assert actual.surface_attributes == original.surface_attributes
    assert actual.parameter_provenance == original.parameter_provenance
    for name in ("time_edges_s", "x_edges_m", "temperature_k", "coverage_m2", "released_m2",
                 "total_wall_energy_debit_j_m2", "inversion_iterations", "inversion_regime", "inversion_log_time_residual"):
        np.testing.assert_array_equal(getattr(actual, name), getattr(original, name))
    with h5py.File(path, "r") as handle:
        assert handle.attrs["schema"] == SCHEMA
        assert handle.attrs["law_identity"] == LAW_IDENTITY
        assert "not an uncertainty" in handle.attrs["inversion_log_time_residual_interpretation"]
        assert handle["inversion_iterations"].dtype == np.dtype("int64")
        assert handle["inversion_regime"].dtype == np.dtype("int8")
        assert handle["inversion_log_time_residual"].attrs["units"] == "1"
        for name, value in NUMERICAL_QUALIFICATIONS.items():
            assert handle.attrs[name] == value
    with pytest.raises(ValueError):
        actual.released_m2.flags.writeable = True


def test_temkin_and_constant_hazard_interchanges_cannot_be_mixed(tmp_path):
    path = tmp_path / "temkin.h5"
    result = source()
    write_temkin_history(path, result)
    with pytest.raises(ValueError):
        load_release_history(path)
    with pytest.raises(ValueError):
        write_release_history(tmp_path / "mislabel.h5", result)
    surface = FlashSurfaceHistory([0, 1e-8], [0, 1], [[600], [600]], "synthetic", "legacy", "manufactured", "face")
    first_order = evolve_surface_release(surface, SurfaceInventory([[1e18]]), [BindingReservoir("sites", "M", 1e8, .3, .2)],
                                        parameter_provenance="manufactured")
    legacy_path = tmp_path / "first_order.h5"
    write_release_history(legacy_path, first_order)
    with pytest.raises(ValueError):
        load_temkin_history(legacy_path)
    with pytest.raises(ValueError):
        write_temkin_history(tmp_path / "mislabel2.h5", first_order)


@pytest.mark.parametrize("mutation", [
    "law_identity", "time_model", "schema", "missing_law", "units", "missing_array", "kinetics",
    "energy", "diagnostic", "uniform", "metadata", "codebook", "qualification", "nonfinite_json", "duplicate_json",
])
def test_tampered_temkin_contract_rejects(tmp_path, mutation):
    path = tmp_path / "bad.h5"
    write_temkin_history(path, source())
    with h5py.File(path, "a") as handle:
        if mutation in {"law_identity", "time_model", "schema"}:
            handle.attrs[mutation] = "constant_hazard"
        elif mutation == "missing_law":
            del handle.attrs["law_identity"]
        elif mutation == "units":
            handle["coverage_m2"].attrs["units"] = "cm^-2"
        elif mutation == "missing_array":
            del handle["inversion_iterations"]
        elif mutation == "kinetics":
            pools = json.loads(handle.attrs["reservoirs_json"])
            pools[0]["attempt_frequency_s"] *= 2
            handle.attrs["reservoirs_json"] = json.dumps(pools)
        elif mutation == "energy":
            handle["total_wall_energy_debit_j_m2"][0, 0, 0] *= 1.1
        elif mutation == "diagnostic":
            handle["inversion_regime"][0, 0, 0] = 7
        elif mutation == "uniform":
            handle["matched_uniform/released_m2"][0, 0, 0] *= 1.1
        elif mutation == "metadata":
            attrs = json.loads(handle.attrs["surface_attributes_json"])
            attrs.pop("current_averaging_interval_m_json")
            handle.attrs["surface_attributes_json"] = json.dumps(attrs)
        elif mutation == "codebook":
            handle.attrs["inversion_regime_definitions_json"] = "{}"
        elif mutation == "qualification":
            handle.attrs["inversion_log_residual_tolerance"] = 1e-3
        elif mutation == "nonfinite_json":
            handle.attrs["surface_attributes_json"] = '{"bad":NaN}'
        else:
            handle.attrs["surface_attributes_json"] = '{"schema":"one","schema":"two"}'
    with pytest.raises(ValueError):
        load_temkin_history(path)


def test_equivalent_inversion_evidence_does_not_require_the_same_library_iteration_count(tmp_path):
    path = tmp_path / "portable.h5"
    write_temkin_history(path, source())
    with h5py.File(path, "a") as handle:
        root = handle["inversion_regime"][...] == 3
        assert np.any(root)
        iterations = handle["inversion_iterations"][...]
        iterations[root] += 1
        handle["inversion_iterations"][...] = iterations
        residual = handle["inversion_log_time_residual"][...]
        residual[root] = 1e-12
        handle["inversion_log_time_residual"][...] = residual
    actual = load_temkin_history(path)
    np.testing.assert_array_equal(actual.inversion_iterations, iterations)
    np.testing.assert_array_equal(actual.inversion_log_time_residual, residual)


@pytest.mark.parametrize("dependency", ["external_link", "soft_link", "virtual", "external_storage", "external_group"])
def test_unhashed_hdf_dependencies_reject(tmp_path, dependency):
    path = tmp_path / "dependent.h5"
    write_temkin_history(path, source())
    external = tmp_path / "external.h5"
    with h5py.File(path, "a") as handle:
        values = handle["temperature_k"][...]
        if dependency == "external_group":
            with h5py.File(external, "w") as target:
                handle.copy("matched_uniform", target)
            del handle["matched_uniform"]
            handle["matched_uniform"] = h5py.ExternalLink(str(external), "/matched_uniform")
        elif dependency == "external_storage":
            del handle["temperature_k"]
            dataset = handle.create_dataset("temperature_k", shape=values.shape, dtype=float,
                                            external=[(str(tmp_path / "temperature.raw"), 0, h5py.h5f.UNLIMITED)])
            dataset[...] = values
            dataset.attrs["units"] = "K"
        else:
            with h5py.File(external, "w") as target:
                target.create_dataset("temperature", data=values).attrs["units"] = "K"
            del handle["temperature_k"]
            if dependency == "external_link":
                handle["temperature_k"] = h5py.ExternalLink(str(external), "/temperature")
            elif dependency == "soft_link":
                handle.create_dataset("local_copy", data=values).attrs["units"] = "K"
                handle["temperature_k"] = h5py.SoftLink("/local_copy")
            else:
                layout = h5py.VirtualLayout(shape=values.shape, dtype=float)
                layout[...] = h5py.VirtualSource(str(external), "temperature", shape=values.shape)
                handle.create_virtual_dataset("temperature_k", layout).attrs["units"] = "K"
    with pytest.raises(ValueError, match="hard links|self-contained"):
        load_temkin_history(path)
