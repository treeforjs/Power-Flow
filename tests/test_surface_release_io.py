import json

import h5py
import numpy as np
import pytest

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.surface_release_history import evolve_surface_release
from mhdlab.surface_release_io import load_release_history, write_release_history
from mhdlab.thermal_release import BindingReservoir, SurfaceInventory


def history():
    surface = FlashSurfaceHistory([0, 1e-8, 3e-8], [-0.0025, -0.001, 0.0025],
                                 [[300, 400], [1000, 600], [1100, 800]], "synthetic",
                                 "test", "manufactured temperature", "cathode face")
    pool = BindingReservoir("weak sites", "H2O:explicit-test-state", 1e12, 0.6, 0.4, 0.02)
    return evolve_surface_release(surface, SurfaceInventory([[1e18, 2e18]]), [pool],
                                  parameter_provenance="manufactured source verification", subdivisions_per_interval=3)


def test_round_trip_preserves_source_and_uniform_moment_null(tmp_path):
    path = tmp_path / "source.h5"
    original = history()
    write_release_history(path, original)
    actual = load_release_history(path)
    assert actual.reservoirs == original.reservoirs
    assert actual.parameter_provenance == original.parameter_provenance
    assert dict(actual.surface_attributes) == dict(original.surface_attributes)
    for name in ("released_m2", "coverage_m2", "temperature_k", "total_wall_energy_debit_j_m2"):
        np.testing.assert_array_equal(getattr(actual, name), getattr(original, name))
    for name in actual.matched_uniform():
        np.testing.assert_array_equal(actual.matched_uniform()[name], original.matched_uniform()[name])
    with pytest.raises(ValueError):
        actual.coverage_m2.flags.writeable = True


@pytest.mark.parametrize("mutation", ["units", "missing", "schema", "time_model", "inventory", "energy", "uniform", "source_kind", "parameters"])
def test_malformed_contract_rejected(tmp_path, mutation):
    path = tmp_path / "bad.h5"
    write_release_history(path, history())
    with h5py.File(path, "a") as f:
        if mutation == "units":
            f["released_m2"].attrs["units"] = "cm^-2"
        elif mutation == "missing":
            del f["released_m2"]
        elif mutation in {"schema", "time_model"}:
            f.attrs[mutation] = "unsupported"
        elif mutation == "inventory":
            f["coverage_m2"][1, 0, 0] *= 0.9
        elif mutation == "energy":
            f["translational_energy_j_m2"][0, 0, 0] *= 1.1
        elif mutation == "uniform":
            f["matched_uniform/released_m2"][0, 0, 0] *= 1.1
        elif mutation == "source_kind":
            attrs = json.loads(f.attrs["surface_attributes_json"])
            attrs["source_kind"] = "unclassified"
            f.attrs["surface_attributes_json"] = json.dumps(attrs)
        elif mutation == "parameters":
            f.attrs["reservoirs_json"] = "not JSON"
    with pytest.raises(ValueError):
        load_release_history(path)


def test_missing_file_rejected(tmp_path):
    with pytest.raises(ValueError):
        load_release_history(tmp_path / "absent.h5")


@pytest.mark.parametrize("mutation", ["missing_current", "missing_depth", "bad_proxy", "boolean_interval", "missing_version"])
def test_semantic_metadata_cannot_bypass_the_surface_contract(tmp_path, mutation):
    path = tmp_path / "semantic.h5"
    write_release_history(path, history())
    with h5py.File(path, "a") as f:
        attrs = json.loads(f.attrs["surface_attributes_json"])
        if mutation == "missing_current":
            del attrs["current_averaging_interval_m_json"]
        elif mutation == "missing_depth":
            del attrs["material_sample_depth_m"]
        elif mutation == "missing_version":
            del attrs["flash_version"]
        elif mutation == "bad_proxy":
            attrs.update(temperature_location="boundary_temperature_proxy", material_sample_depth_m=-1,
                         material_sample_depth_units="cm", boundary_temperature_reference=" ")
        elif mutation == "boolean_interval":
            attrs.update(current_averaging_interval_m_json="false", current_uniformity_reference=True)
        f.attrs["surface_attributes_json"] = json.dumps(attrs)
    with pytest.raises(ValueError):
        load_release_history(path)
