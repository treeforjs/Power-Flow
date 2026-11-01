import numpy as np
import pytest
from scipy.optimize import brentq

from mhdlab.thermal_expansion import load_ss304_material
from mhdlab.thermoelastic_energy import (
    algorithmic_exchange, isotropic_stiffness, thermodynamic_state,
    ss304_thermodynamic_increment,
)


def test_exact_discrete_exchange_variable_anisotropic_stiffness():
    rng = np.random.default_rng(207)
    for _ in range(20):
        a, b = rng.normal(size=(2,6,6))
        c0, c1 = (a.T@a+np.eye(6))*1e10, (b.T@b+np.eye(6))*1e10
        eps0, eps1, th0, th1 = rng.normal(scale=1e-3, size=(4,6))
        v = 2.5e-8
        out = algorithmic_exchange(eps0, eps1, th0, th1, c0, c1, reference_volume_m3=v)
        e0, e1 = eps0-th0, eps1-th1
        independently_computed = v*(.5*e1@c1@e1-.5*e0@c0@e0)
        assert out['delta_phi_J'] == pytest.approx(independently_computed, abs=1e-16)
        assert out['algorithmic_heat_debit_J'] == pytest.approx(
            out['thermal_eigenstrain_term_J']+out['modulus_term_J'], abs=2e-16)
        assert abs(out['closure_residual_J']) < 2e-16
        # Reversing a step reverses every exchange, including the modulus term.
        reverse = algorithmic_exchange(eps1, eps0, th1, th0, c1, c0, reference_volume_m3=v)
        assert reverse['algorithmic_heat_debit_J'] == pytest.approx(-out['algorithmic_heat_debit_J'])


def test_same_temperature_pure_mechanical_step_has_no_algorithmic_heat():
    c = isotropic_stiffness(200e9,.29)
    eps0 = np.zeros(6)
    eps1 = np.array([.001, 0, 0, .002, 0, 0])
    out = algorithmic_exchange(eps0,eps1,eps0,eps0,c,c,reference_volume_m3=1e-8)
    assert out['algorithmic_heat_debit_J'] == pytest.approx(0, abs=1e-18)
    # Engineering shear gamma contributes mu*gamma^2/2, not 2*mu*gamma^2.
    expected = 1e-8*(c[0,0]*.001**2+c[3,3]*.002**2)/2
    assert out['delta_phi_J'] == pytest.approx(expected)


def _nonlinear_fixture(t, eps):
    d = t-300
    base = isotropic_stiffness(200e9,.29)
    factor, fp, fpp = 1-2e-4*d+4e-8*d*d, -2e-4+8e-8*d, 8e-8
    normal = np.array([1,1,1,0,0,0])
    th = (16e-6*d+2e-9*d*d)*normal
    th_t = (16e-6+4e-9*d)*normal
    th_tt = 4e-9*normal
    return thermodynamic_state(temperature_K=t,strain=eps,thermal_strain=th,
        dthermal_strain_dT=th_t,d2thermal_strain_dT2=th_tt,
        stiffness_Pa=factor*base,dstiffness_dT_Pa_K=fp*base,d2stiffness_dT2_Pa_K2=fpp*base,
        base_energy_J_m3=3.95e6*d,base_capacity_J_m3_K=3.95e6)


def test_thermodynamic_derivatives_against_independent_finite_differences():
    eps = np.array([.001,-.0003,.0002,.0001,-.0002,.0004])
    t, dt = 341.0, .01
    state, plus, minus = (_nonlinear_fixture(x,eps) for x in (t,t+dt,t-dt))
    first = (plus['phi_J_m3']-minus['phi_J_m3'])/(2*dt)
    second = (plus['phi_J_m3']-2*state['phi_J_m3']+minus['phi_J_m3'])/dt**2
    assert state['phi_temperature_derivative_J_m3_K'] == pytest.approx(first, rel=2e-8)
    assert state['phi_second_temperature_derivative_J_m3_K2'] == pytest.approx(second, rel=1e-7)
    assert state['stress_temperature_derivative_Pa_K'] == pytest.approx(
        (plus['stress_Pa']-minus['stress_Pa'])/(2*dt),rel=2e-8)
    assert state['fixed_strain_capacity_J_m3_K'] == pytest.approx(
        (plus['total_internal_energy_J_m3']-minus['total_internal_energy_J_m3'])/(2*dt),rel=2e-8)
    for i in range(6):
        de = np.zeros(6); de[i] = 1e-7
        up = _nonlinear_fixture(t,eps+de)['total_internal_energy_J_m3']
        um = _nonlinear_fixture(t,eps-de)['total_internal_energy_J_m3']
        assert (up-um)/(2e-7) == pytest.approx(
            state['stress_Pa'][i]-t*state['stress_temperature_derivative_Pa_K'][i],rel=1e-7)


def _constant_state(t,eps, tref=320.):
    return thermodynamic_state(temperature_K=t,strain=eps,thermal_strain=16e-6*(t-tref),
        dthermal_strain_dT=16e-6,d2thermal_strain_dT2=0.,stiffness_Pa=200e9,
        dstiffness_dT_Pa_K=0.,d2stiffness_dT2_Pa_K2=0.,
        base_energy_J_m3=3.95e6*(t-tref),base_capacity_J_m3_K=3.95e6)


def test_adiabatic_expansion_cools_and_refines_to_independent_entropy_limit():
    tref, eps_end, ea = 320., .001, 200e9*16e-6
    # Independent reversible endpoint: constant total entropy, not a work-port identity.
    exact = brentq(lambda t: 3.95e6*np.log(t/tref)+ea*(eps_end-16e-6*(t-tref)),319.,320.,xtol=1e-13)
    errors=[]
    for steps in (1,2,4):
        t, eps = tref, 0.
        for end_eps in np.linspace(0,eps_end,steps+1)[1:]:
            old = _constant_state(t,eps)
            def residual(end_t):
                new = _constant_state(end_t,end_eps)
                work = .5*(old['stress_Pa'][0]+new['stress_Pa'][0])*(end_eps-eps)
                return new['total_internal_energy_J_m3']-old['total_internal_energy_J_m3']-work
            new_t = brentq(residual,319.,320.,xtol=1e-13)
            assert new_t < t
            assert abs(residual(new_t)) < 3e-7
            t, eps = new_t, end_eps
        errors.append(abs(t-exact))
    assert errors[0]/errors[1] > 3.8
    assert errors[1]/errors[2] > 3.8


def test_physical_and_algorithmic_heat_are_different_even_isothermally():
    material=load_ss304_material()
    out=ss304_thermodynamic_increment(np.zeros(6),np.array([.001,0,0,0,0,0]),320.,320.,
        material=material,stress_free_temperature_K=320.,poisson_ratio=.29,
        reference_volume_m3=1e-8,base_energy_J_m3=lambda t:3.95e6*(t-300),
        base_capacity_J_m3_K=lambda t:3.95e6)
    assert out['algorithmic_heat_debit_J'] == pytest.approx(0,abs=1e-18)
    assert out['conditional_heat_debit_J'] > 0
    assert out['energy_balance_residual_J'] == out['conditional_heat_debit_J']
    assert out['physical_EOS_ready'] is False


def test_ss304_segment_increment_matches_direct_internal_energy_and_mechanical_balance():
    material=load_ss304_material()
    out=ss304_thermodynamic_increment(np.zeros(6),np.full(6,.0001),320.,325.,
        material=material,stress_free_temperature_K=300.,poisson_ratio=.29,
        reference_volume_m3=2e-8,base_energy_J_m3=lambda t:3.95e6*(t-300),
        base_capacity_J_m3_K=lambda t:3.95e6,supplied_heat_J=.2)
    expected=2e-8*(out['state1']['total_internal_energy_J_m3']-out['state0']['total_internal_energy_J_m3'])-out['work_stress_J']-.2
    assert out['energy_balance_residual_J'] == pytest.approx(expected,abs=1e-16)
    wext=.002
    delta_k=-out['work_stress_J']+wext
    assert out['conditional_heat_debit_J'] == pytest.approx(delta_k+out['delta_elastic_internal_energy_J']-wext)


@pytest.mark.parametrize('t0,t1',[(373.15,374.),(372.,374.),(300.,293.15),(773.15,770.),(300.,800.)])
def test_physical_mode_rejects_knots_crossing_and_unsupported_temperature(t0,t1):
    with pytest.raises(ValueError):
        ss304_thermodynamic_increment(np.zeros(6),np.zeros(6),t0,t1,
            material=load_ss304_material(),stress_free_temperature_K=300.,poisson_ratio=.29,
            reference_volume_m3=1.,base_energy_J_m3=lambda t:3.95e6*(t-300),
            base_capacity_J_m3_K=lambda t:3.95e6)


def test_piecewise_law_has_artificial_internal_energy_jump_at_knot():
    states=[]
    for t in (373.15-1e-7,373.15+1e-7):
        out=ss304_thermodynamic_increment(np.zeros(6),np.zeros(6),t,t,
            material=load_ss304_material(),stress_free_temperature_K=300.,poisson_ratio=.29,
            reference_volume_m3=1.,base_energy_J_m3=lambda t:3.95e6*(t-300),
            base_capacity_J_m3_K=lambda t:3.95e6)
        states.append(out['state0'])
    assert abs(states[1]['phi_J_m3']-states[0]['phi_J_m3']) < .02
    assert abs(states[1]['elastic_internal_energy_J_m3']-states[0]['elastic_internal_energy_J_m3']) > 1e5


@pytest.mark.parametrize('bad',[-1.,float('nan'),float('inf'),0.])
def test_rejects_nonpositive_or_nonfinite_reference_volume(bad):
    with pytest.raises(ValueError):
        algorithmic_exchange(0.,.001,0.,0.,200e9,200e9,reference_volume_m3=bad)


def test_rejects_nonpositive_fixed_strain_capacity_and_nonsymmetric_stiffness():
    with pytest.raises(ValueError,match='capacity'):
        thermodynamic_state(temperature_K=300.,strain=0.,thermal_strain=0.,
            dthermal_strain_dT=16e-6,d2thermal_strain_dT2=0.,stiffness_Pa=200e9,
            dstiffness_dT_Pa_K=0.,d2stiffness_dT2_Pa_K2=0.,
            base_energy_J_m3=0.,base_capacity_J_m3_K=1.)
    with pytest.raises(ValueError,match='symmetric'):
        algorithmic_exchange([0,0],[0,0],[0,0],[0,0],[[1,1],[0,1]],np.eye(2),reference_volume_m3=1.)
