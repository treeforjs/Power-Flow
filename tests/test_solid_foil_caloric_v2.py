"""Run centrally in frozen userspace; this module never builds the owner library."""
import json
import os
from pathlib import Path

import numpy as np
import pytest

import mhdlab.solid_foil_caloric_v2 as adapter
from mhdlab.thermal_expansion import load_ss304_material


SEAL = Path(os.environ.get("SS304_V2_SOURCE_SEAL", "/home/treefor/Documents/Codex/2026-10-28/we/outputs/ss304-thermal-properties-v2-20261029/source-manifest.json"))


def _construct(**updates):
    args = dict(source_seal_path=SEAL, expected_source_seal_sha256=adapter.SOURCE_SEAL_SHA256,
                reference_density_kg_m3=7900., qualification_mode="nominal_1800")
    args.update(updates)
    return adapter.SealedNominalCaloricV2(**args)


@pytest.fixture(scope="module")
def model():
    if not SEAL.is_file():
        pytest.skip("sealed external thermal-v2 runtime is unavailable; no owner build is attempted")
    return _construct()


def _solid_u(t):
    return (t-300.)*(472.+.068*(t+300.)-2820000./(t*300.))


def test_cold_reference_is_negative_and_not_clamped(model):
    low = model.evaluate_T(298.2)
    assert low["u_model_j_kg"] == pytest.approx(-866.0792373440698, abs=2e-10)
    assert low["energy_j_m3_reference"] == pytest.approx(7900.*low["u_model_j_kg"])
    assert low["phase"] == 0
    assert model.energy_J_m3(300.) == 0.
    assert model.temperature_K(model.energy_J_m3(298.2)) == 298.2
    assert model.temperature_K(model.energy_J_m3(299.1)) == pytest.approx(299.1, abs=1e-12)
    receipt = model.receipt
    assert receipt["negative_cold_reference_energy_valid"] is True
    assert receipt["mechanics_range_extended"] is False


def test_direct_abi_matches_independent_solid_integral_capacity_and_cold_transport(model):
    temperatures = np.array([[298.2, 299., 300.], [350., 773.15, 1500.]])
    np.testing.assert_allclose(model.energy_J_m3(temperatures), 7900.*_solid_u(temperatures), rtol=2e-15, atol=1e-8)
    cp = 472.+.136*temperatures-2820000./temperatures**2
    np.testing.assert_allclose(model.capacity_J_m3_K(temperatures),7900.*cp,rtol=2e-15)
    t = 299.
    f = (t-293.)/7.
    assert model.conductivity_W_m_K(t) == pytest.approx(14.76+f*(14.89-14.76))
    assert model.electrical_conductivity_S_m(t) == pytest.approx(1/((71.3+f*(71.9-71.3))*1e-8))
    h = 1e-3
    derivative = (model.energy_J_m3(t+h)-model.energy_J_m3(t-h))/(2*h)
    assert derivative == pytest.approx(model.capacity_J_m3_K(t),rel=3e-8)


def test_q1_phase_energy_and_inverse_remain_explicit(model):
    solidus, liquidus = 1672.,1727.
    cp_s = 472.+.136*solidus-2820000./solidus**2
    expected_u = _solid_u(solidus)+.5*(cp_s+800.)*(liquidus-solidus)+160.*liquidus+800.*(1800.-liquidus)
    assert model.evaluate_T(1800.)['u_model_j_kg'] == pytest.approx(expected_u,rel=2e-15)
    t = np.array([298.2, 299.7, 300., 500., 1672., 1700., 1727., 1800.])
    np.testing.assert_allclose(model.temperature_K(model.energy_J_m3(t)),t,rtol=0,atol=5e-12)
    assert model.evaluate_T(1700.)['phase'] == 1
    assert model.evaluate_T(1800.)['phase'] == 2
    # Thermal support does not remove the independent solid-mechanics guard.
    assert model.energy_J_m3(800.) > 0
    with pytest.raises(ValueError,match='supports only'):
        load_ss304_material().check_temperature(800.)


def test_density_mass_and_volume_are_applied_once(model):
    other = _construct(reference_density_kg_m3=3950.)
    t = 299.
    assert other.energy_J_m3(t) == .5*model.energy_J_m3(t)
    assert other.capacity_J_m3_K(t) == .5*model.capacity_J_m3_K(t)
    volume = 2.5e-8
    mass = model.rho*volume
    total = model.energy_J_m3(t)*volume
    assert total == pytest.approx(mass*model.evaluate_T(t)['u_model_j_kg'])
    support = model.assess_energy(current_mass_kg=mass,total_energy_J=total)
    assert support['range_status'] == 0
    assert support['cold_margin_j'] > 0
    assert total < 0  # Negative reference energy is not a negative metal inventory.
    unsupported = model.assess_energy(current_mass_kg=mass,total_energy_J=support['minimum_j']-1e-7)
    assert unsupported['range_status'] == -1
    assert unsupported['cold_margin_j'] == pytest.approx(-1e-7,abs=1e-16)
    receipt = model.receipt
    receipt['reference_density_kg_m3'] = 1.
    assert model.receipt['reference_density_kg_m3'] == 7900.
    with pytest.raises(AttributeError):
        model.rho = 1.


@pytest.mark.parametrize('temperature',[298.1,1800.1,float('nan'),float('inf')])
def test_temperature_domain_rejection(model,temperature):
    with pytest.raises(ValueError):
        model.energy_J_m3(temperature)


def test_negative_energy_below_cold_bound_rejected(model):
    with pytest.raises(ValueError):
        model.temperature_K(model.energy_J_m3(298.2)-1.)
    with pytest.raises(ValueError):
        model.assess_energy(current_mass_kg=-1,total_energy_J=0.)


@pytest.mark.parametrize('mode',['evaluated_solid_1500','source_conditioned_liquid_3000',1,None])
def test_wrong_qualification_rejected_before_loading(mode):
    with pytest.raises(ValueError,match='q1 nominal_1800'):
        _construct(qualification_mode=mode)


def test_wrong_revision_rejected_before_loading():
    with pytest.raises(ValueError,match='reviewed thermal-v2'):
        _construct(expected_source_seal_sha256='0'*64)


def test_seal_runtime_and_library_tampering_rejected(model,tmp_path,monkeypatch):
    bad = tmp_path/'source-manifest.json'
    bad.write_text('{}')
    with pytest.raises(ValueError,match='source seal SHA256'):
        _construct(source_seal_path=bad)
    runtime = tmp_path/'runtime-bindings.json'
    runtime.write_text('{}')
    with pytest.raises(ValueError,match='runtime bindings SHA256'):
        _construct(runtime_bindings_path=runtime)
    library = Path(json.loads(SEAL.read_text())['compiled_library']['path'])
    real_sha = adapter._sha
    monkeypatch.setattr(adapter,'_sha',lambda p:'0'*64 if Path(p)==library else real_sha(p))
    with pytest.raises(ValueError,match='compiled library differs'):
        _construct()
