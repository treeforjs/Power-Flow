"""Native lumped-energy projection, identity, sign and depth-once controls."""
from dataclasses import replace
import hashlib
import json

import numpy as np
import pytest

from mhdlab.solid_foil_mechanics import structured_hex_mesh, mesh_geometry_sha256
from mhdlab.solid_foil_native_mapping import (
    make_native_heat_stage_identity, native_lumped_cell_energy_J,
    project_native_heat_to_solid, signed_cell_exchange_to_native,
)


class Caloric:
    receipt = {"source_seal_sha256": "c" * 64}

    def energy_J_m3(self, temperature):
        d = np.asarray(temperature) - 300.
        return 4e6 * d + 1e4 * d ** 2

    def capacity_J_m3_K(self, temperature):
        return 4e6 + 2e4 * (np.asarray(temperature) - 300.)

    def temperature_K(self, energy):
        # Stable quadratic inverse, valid on this manufactured monotone segment.
        energy = np.asarray(energy)
        return 300. + 2 * energy / (4e6 + np.sqrt(16e12 + 4e4 * energy))


def fixture(stage="trial"):
    mesh = structured_hex_mesh([-.0025, -.001, .0025], [.0025, .00257, .0027], [-.0125, -.01, 0., .0125])
    identity = make_native_heat_stage_identity(mesh, native_x_edges_m=mesh.x_edges_m,
        native_y_edges_m=mesh.y_edges_m, native_height_m=.025,
        native_mesh_revision="native-grid-0", solid_mesh_revision="solid-grid-0", state_revision=7,
        interval_s=(1e-9, 2e-9), stage_time_s=1e-9 if stage == "accepted" else 2e-9,
        stage=stage, material_sha256="a" * 64, source_sha256="b" * 64, caloric_sha256="c" * 64)
    return mesh, identity


def project(mesh, identity, temperature, **changes):
    caloric = Caloric()
    native = native_lumped_cell_energy_J(temperature, x_edges_m=mesh.x_edges_m,
                                       y_edges_m=mesh.y_edges_m, height_m=.025, caloric=caloric)
    arguments = dict(native_x_edges_m=mesh.x_edges_m, native_y_edges_m=mesh.y_edges_m,
                     native_height_m=.025, caloric=caloric, identity=identity, expected_identity=identity,
                     native_cell_base_energy_J=native, native_stored_energy_J=float(native.sum()),
                     admissible_temperature_range_K=[298.2, 773.15])
    arguments.update(changes)
    return project_native_heat_to_solid(mesh, temperature, **arguments)


@pytest.mark.parametrize("temperature", [298.5, 300., 330.])
def test_uniform_native_temperature_maps_exactly_without_depth_duplication(temperature):
    mesh, identity = fixture()
    result = project(mesh, identity, np.full((3, 3), temperature))
    np.testing.assert_allclose(result.cell_temperature_K, temperature, rtol=0, atol=6e-14)
    assert result.cell_base_energy_J.shape == mesh.cell_shape
    volume = mesh.reference_volumes_m3.reshape(mesh.cell_shape)
    np.testing.assert_allclose(result.cell_base_energy_J, volume * Caloric().energy_J_m3(temperature), rtol=5e-14, atol=0.)
    assert result.receipt()["native_stored_energy_J"] == pytest.approx(float(result.cell_base_energy_J.sum()), rel=1e-14)


def test_nonlinear_caloric_average_preserves_native_lumped_energy_not_mean_temperature():
    mesh, identity = fixture()
    temperature = np.array([[300., 305., 310.], [310., 320., 330.], [330., 335., 340.]])
    result = project(mesh, identity, temperature)
    mean_temperature = (temperature[:-1, :-1] + temperature[:-1, 1:] + temperature[1:, 1:] + temperature[1:, :-1]) / 4
    assert np.all(result.cell_temperature_K[0] > mean_temperature)
    # Independent nodal dual-volume sum includes boundary half/quarter ownership.
    xw = np.r_[np.diff(mesh.x_edges_m), 0.] / 2 + np.r_[0., np.diff(mesh.x_edges_m)] / 2
    yw = np.r_[np.diff(mesh.y_edges_m), 0.] / 2 + np.r_[0., np.diff(mesh.y_edges_m)] / 2
    native_total = np.sum(Caloric().energy_J_m3(temperature) * np.outer(yw, xw)) * .025
    assert result.cell_base_energy_J.sum() == pytest.approx(native_total, rel=1e-14)
    reconstructed = np.sum(Caloric().energy_J_m3(result.cell_temperature_K) * mesh.reference_volumes_m3.reshape(mesh.cell_shape))
    assert reconstructed == pytest.approx(native_total, rel=1e-13)
    assert result.receipt()["physical_EOS_ready"] is False


def test_projection_is_deterministic_immutable_and_binds_inputs():
    mesh, identity = fixture()
    temperature = np.arange(9.).reshape(3, 3) + 300.
    first, second = project(mesh, identity, temperature), project(mesh, identity, temperature)
    assert first.receipt_bytes == second.receipt_bytes
    with pytest.raises(ValueError):
        first.cell_temperature_K.flat[0] = 0.
    changed = project(mesh, identity, temperature + 1.)
    assert changed.receipt()["input_sha256"] != first.receipt()["input_sha256"]


def test_native_cell_energy_mismatch_is_caught_even_if_global_total_matches():
    mesh, identity = fixture()
    temperature = np.full((3, 3), 310.)
    native = np.array(native_lumped_cell_energy_J(temperature, x_edges_m=mesh.x_edges_m,
                         y_edges_m=mesh.y_edges_m, height_m=.025, caloric=Caloric()))
    native[0, 0] += 1e-5
    native[0, 1] -= 1e-5
    with pytest.raises(ValueError, match="exact lumped"):
        project(mesh, identity, temperature, native_cell_base_energy_J=native)
    with pytest.raises(ValueError, match="stored-energy"):
        project(mesh, identity, temperature, native_stored_energy_J=42.)


@pytest.mark.parametrize("field,value", [("source_sha256", "d" * 64), ("material_sha256", "d" * 64),
    ("caloric_sha256", "d" * 64), ("state_revision", 8), ("native_mesh_revision", "other"), ("solid_mesh_revision", "other")])
def test_stage_and_source_identities_cannot_be_mixed(field, value):
    mesh, identity = fixture()
    with pytest.raises(ValueError, match="identity mismatch"):
        project(mesh, identity, np.full((3, 3), 310.), expected_identity=replace(identity, **{field: value}))


def test_invalid_stage_clock_geometry_temperature_and_caloric_are_rejected():
    mesh, identity = fixture()
    with pytest.raises(ValueError, match="stage"):
        replace(identity, stage_time_s=1.5e-9)
    with pytest.raises(ValueError, match="mesh"):
        project(mesh, identity, np.full((3, 3), 310.), native_height_m=.026)
    with pytest.raises(ValueError, match="outside solid"):
        project(mesh, identity, np.array([[800., 300., 300.], [300., 300., 300.], [300., 300., 300.]]))
    class WrongInverse(Caloric):
        def temperature_K(self, energy):
            return super().temperature_K(energy) + .1
    with pytest.raises(ValueError, match="roundtrip"):
        project(mesh, identity, np.full((3, 3), 310.), caloric=WrongInverse())


def exchange(mesh, identity):
    packet = dict(schema="solid-thermoelastic-cell-exchange-v2", state_revision=identity.state_revision,
        mesh_revision=identity.solid_mesh_revision, mesh_geometry_sha256=mesh_geometry_sha256(mesh),
        material_sha256=identity.material_sha256, source_sha256=identity.source_sha256,
        caloric_model={"source_seal_sha256": identity.caloric_sha256}, interval_s=list(identity.interval_s),
        stage_time_s=identity.stage_time_s, conditional_thermodynamic_exchange_ready=True,
        species_mutation="none", species_event_ids=[], species_parent_ids=[], transaction_id="work-7",
        cell_heat_to_thermal_J=(np.arange(len(mesh.cells)) - 5.).tolist())
    packet["payload_sha256"] = hashlib.sha256(json.dumps(packet, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    return packet


def native_exchange(mesh, identity, packet):
    return signed_cell_exchange_to_native(mesh, packet, native_x_edges_m=mesh.x_edges_m,
        native_y_edges_m=mesh.y_edges_m, native_height_m=.025, identity=identity, expected_identity=identity)


def test_signed_mechanical_exchange_preserves_local_values_sign_and_total():
    mesh, identity = fixture()
    packet = exchange(mesh, identity)
    result = native_exchange(mesh, identity, packet)
    expected = np.asarray(packet["cell_heat_to_thermal_J"]).reshape(mesh.cell_shape).sum(axis=0)
    np.testing.assert_array_equal(result["cell_heat_J"], expected)
    assert result["source_total_J"] == result["target_total_J"]
    assert result["signed_cell_load_class"] == "reversible_thermoelastic_internal_exchange"
    assert result["native_state_mutated"] is False


@pytest.mark.parametrize("field,value", [("source_sha256", "d" * 64), ("material_sha256", "d" * 64),
    ("stage_time_s", 1.5e-9), ("interval_s", [0., 2e-9]), ("species_event_ids", ["particle"]),
    ("state_revision", 8), ("mesh_revision", "other")])
def test_signed_exchange_rejects_mismatched_or_stale_work(field, value):
    mesh, identity = fixture()
    packet = exchange(mesh, identity)
    packet[field] = value
    with pytest.raises(ValueError):
        native_exchange(mesh, identity, packet)


def test_signed_exchange_rejects_accepted_stage_and_tampered_payload():
    mesh, identity = fixture(stage="accepted")
    with pytest.raises(ValueError, match="trial"):
        native_exchange(mesh, identity, exchange(mesh, identity))
    mesh, identity = fixture()
    packet = exchange(mesh, identity)
    packet["cell_heat_to_thermal_J"][0] = 42.
    with pytest.raises(ValueError, match="altered exchange"):
        native_exchange(mesh, identity, packet)
