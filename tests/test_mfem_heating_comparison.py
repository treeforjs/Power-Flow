"""MFEM provenance is explicit while analytic comparison stays unchanged."""
import json
import math

import numpy as np
import pytest

from mhdlab.heating_benchmarks import ConstantConductivityHalfSpace
from mhdlab.heating_comparison import (HeatingProfile, compare_heating_profile,
    write_heating_profile_csv, load_heating_profile_csv)


def profile(**changes):
    values = dict(time_s=1e-6, depth_edges_m=[0., 1e-4, 3e-4],
                  magnetic_field_t=[1., 2.], current_density_a_m2=[3., 4.],
                  joule_power_w_m3=[5., 6.], source_kind="mfem", run_id="schema-fixture",
                  source_reference="Authored metadata test, not a physical MFEM run",
                  mfem_version="4.10", source_revision="a" * 40)
    values.update(changes)
    return HeatingProfile(**values)


@pytest.mark.parametrize("kind", ["mfem", "synthetic"])
def test_mfem_v2_roundtrip_keeps_backend_separate_from_flash(kind, tmp_path):
    sample = profile(source_kind=kind)
    path, sidecar = write_heating_profile_csv(tmp_path / "profile.csv", sample)
    loaded = load_heating_profile_csv(path)
    assert loaded.metadata == sample.metadata
    assert loaded.metadata["schema"] == "mhdlab-heating-profile-v2"
    assert loaded.metadata["heating_backend"] == "mfem_joule"
    assert loaded.flash_version == ""
    np.testing.assert_array_equal(loaded.joule_power_w_m3, sample.joule_power_w_m3)


def test_mfem_declaration_cannot_be_downgraded_to_v1(tmp_path):
    path, sidecar = write_heating_profile_csv(tmp_path / "profile.csv", profile())
    metadata = json.loads(sidecar.read_text())
    metadata["schema"] = "mhdlab-heating-profile-v1"
    for key in ("heating_backend", "mfem_version", "source_revision"):
        metadata.pop(key)
    sidecar.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="v1 downgrade"):
        load_heating_profile_csv(path)


@pytest.mark.parametrize("changes", [
    {"mfem_version": ""}, {"mfem_version": "pending"}, {"source_revision": ""},
    {"source_revision": "tag-v4.10"}, {"flash_version": "FLASH 4.8"},
    {"source_kind": "flash", "flash_version": "FLASH 4.8"},
])
def test_missing_or_cross_backend_identity_rejects(changes):
    with pytest.raises(ValueError):
        profile(**changes)


def test_mfem_profile_uses_identical_analytic_cell_average_comparison():
    reference = ConstantConductivityHalfSpace(conductivity_s_m=3e6,
        permeability_h_m=4e-7 * math.pi, volumetric_heat_capacity_j_m3_k=2.8e6,
        boundary_field_t=5.)
    edges = [0., 1e-4, 3e-4]
    exact = reference.cell_averages(edges, 1e-6)
    sample = profile(source_kind="synthetic", magnetic_field_t=exact.magnetic_field_t,
                     current_density_a_m2=exact.current_density_a_m2,
                     joule_power_w_m3=exact.joule_power_w_m3)
    result = compare_heating_profile(sample, reference)
    assert result["profile"]["heating_backend"] == "mfem_joule"
    for row in result["errors"].values():
        assert row["absolute"]["linf"] == 0.
