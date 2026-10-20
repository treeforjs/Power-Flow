"""Manufactured physical controls for the local thermal/inventory reference."""
from dataclasses import FrozenInstanceError,replace
from fractions import Fraction
import math

import pytest

from mhdlab.thermal_feedback_reference import (
    FeedbackProblem,
    FeedbackReservoir,
    UnsupportedFeedbackRangeError,
    evolve_feedback,
    feedback_rhs,
)

KB=1.380649e-23
QE=1.602176634e-19


def reservoir(**changes):
    values=dict(name='manufactured_parent',state_id='manufactured_fixed_state',
        attempt_frequency_s=1e7,activation_energy_ev=0.0,desorption_energy_ev=0.12,
        internal_energy_ev=0.03,coverage_coefficient=0.0,
        reference_coverage_m2=4e20,max_coverage_m2=4e20,
        kinetic_provenance='Manufactured kinetics for the local reference only.',
        energy_provenance='Manufactured removal and gas internal energy in one declared reference.')
    values.update(changes)
    return FeedbackReservoir(**values)


def problem(pool=None,**changes):
    values=dict(reservoir=reservoir() if pool is None else pool,
        initial_coverage_m2=2e20,initial_temperature_k=900.0,
        areal_heat_capacity_j_m2_k=0.05,patch_area_m2=2.0**-20,
        temperature_min_k=200.0,temperature_max_k=2000.0,
        parameter_status='manufactured',
        parameter_provenance='Manufactured one-patch reference; no material or water inference.',
        energy_reference='Gas internal energy and parent removal use a common declared manufactured reference.')
    values.update(changes)
    return FeedbackProblem(**values)


@pytest.fixture(scope='module')
def zero_barrier_history():
    return evolve_feedback(problem(),duration_s=1e-7,steps=128)


def test_zero_barrier_joint_solution_and_separate_energy_channels(zero_barrier_history):
    """N remains exponential while wall cooling follows the independent T(N) law."""
    history=zero_barrier_history
    p=history.problem
    pool=p.reservoir
    ca=p.areal_heat_capacity_j_m2_k
    d=pool.desorption_energy_ev*QE
    internal=pool.internal_energy_ev*QE
    a=(d+internal)/(2*KB)
    assert history.steps==128 and isinstance(history.samples,tuple)
    assert len(history.samples)==129
    assert history.samples[0].time_s==0 and history.samples[-1].time_s==1e-7
    for sample in history.samples:
        exponent=-pool.attempt_frequency_s*sample.time_s
        n=p.initial_coverage_m2*math.exp(exponent)
        released=-p.initial_coverage_m2*math.expm1(exponent)
        temperature=p.initial_temperature_k+(p.initial_temperature_k+a)*math.expm1(
            2*KB*(n-p.initial_coverage_m2)/ca)
        translation=ca*(p.initial_temperature_k-temperature)-(d+internal)*released
        assert sample.coverage_m2==pytest.approx(n,rel=8e-10)
        assert sample.released_m2==pytest.approx(released,rel=8e-10,abs=0)
        assert sample.temperature_k==pytest.approx(temperature,rel=8e-10)
        assert sample.binding_energy_j_m2==pytest.approx(d*released,rel=8e-10,abs=2e-14)
        assert sample.internal_energy_j_m2==pytest.approx(internal*released,rel=8e-10,abs=2e-14)
        assert sample.translational_energy_j_m2==pytest.approx(translation,rel=1e-9,abs=2e-14)
        assert sample.total_wall_energy_debit_j_m2==pytest.approx(ca*(p.initial_temperature_k-temperature),rel=1e-9,abs=2e-14)
        assert sample.gas_energy_j_m2==pytest.approx(translation+internal*released,rel=1e-9,abs=2e-14)
        assert sample.adsorbate_reference_energy_j_m2==pytest.approx(-d*sample.coverage_m2,rel=5e-14)
        assert sample.adsorbate_reference_change_j_m2==pytest.approx(d*sample.released_m2,rel=5e-14,abs=0)
        assert sample.adsorbate_reference_inventory_change_j_m2==pytest.approx(d*(p.initial_coverage_m2-sample.coverage_m2),rel=2e-12,abs=1e-14)
        assert sample.wall_thermal_change_j_m2==pytest.approx(ca*(sample.temperature_k-p.initial_temperature_k),rel=5e-14,abs=0)
        energy_scale=max(sample.total_wall_energy_debit_j_m2,1e-20)
        for key in ['total_energy_residual_j_m2','adsorbate_reference_residual_j_m2',
            'wall_debit_residual_j_m2','wall_channel_residual_j_m2',
            'binding_allocation_residual_j_m2','internal_allocation_residual_j_m2']:
            assert math.isfinite(getattr(sample,key))
            assert abs(getattr(sample,key))<=2e-10*energy_scale
        assert abs(math.fsum([sample.wall_thermal_change_j_m2,
            sample.adsorbate_reference_change_j_m2,sample.gas_energy_j_m2]))<=2e-10*energy_scale
    assert history.samples[-1].temperature_k<p.initial_temperature_k


def test_fixed_reference_and_barrier_independent_removal_cost():
    """The rate sees fixed Nref; kinetic Ea does not become a removal energy."""
    pool=reservoir(attempt_frequency_s=1e8,activation_energy_ev=0.35,coverage_coefficient=0.6)
    p=problem(pool)
    n=1e20;t=900.0
    rate=feedback_rhs(p,coverage_m2=n,temperature_k=t)
    expected=pool.attempt_frequency_s*n*math.exp(
        -pool.activation_energy_ev*QE*(1-pool.coverage_coefficient*n/pool.reference_coverage_m2)/(KB*t))
    assert rate['release_flux_m2_s']==pytest.approx(expected,rel=3e-14)
    assert rate['coverage_rate_m2_s']==pytest.approx(-expected,rel=3e-14)
    other_initial=replace(p,initial_coverage_m2=3e20)
    assert feedback_rhs(other_initial,coverage_m2=n,temperature_k=t)==rate
    assert other_initial.reservoir.reference_coverage_m2==p.reservoir.reference_coverage_m2==4e20
    zero=feedback_rhs(problem(reservoir(attempt_frequency_s=1e8)),coverage_m2=n,temperature_k=t)
    assert zero['release_flux_m2_s']>rate['release_flux_m2_s']
    for rhs in [rate,zero]:
        flux=rhs['release_flux_m2_s']
        assert rhs['binding_power_j_m2_s']/flux==pytest.approx(0.12*QE,rel=3e-14)
        assert rhs['internal_power_j_m2_s']/flux==pytest.approx(0.03*QE,rel=3e-14)
        assert rhs['translational_power_j_m2_s']/flux==pytest.approx(2*KB*t,rel=3e-14)
        assert rhs['total_wall_power_j_m2_s']/flux==pytest.approx(0.15*QE+2*KB*t,rel=3e-14)
        assert rhs['temperature_rate_k_s']==pytest.approx(-rhs['total_wall_power_j_m2_s']/p.areal_heat_capacity_j_m2_k,rel=3e-14)
    translation_only=feedback_rhs(problem(reservoir(desorption_energy_ev=0,internal_energy_ev=0)),coverage_m2=n,temperature_k=t)
    assert translation_only['binding_power_j_m2_s']==translation_only['internal_power_j_m2_s']==0
    assert translation_only['total_wall_power_j_m2_s']==translation_only['translational_power_j_m2_s']


def test_area_applied_once_and_detached_immutable_results(zero_barrier_history):
    original=zero_barrier_history
    larger=evolve_feedback(replace(original.problem,patch_area_m2=4*original.problem.patch_area_m2),
        duration_s=original.duration_s,steps=original.steps)
    for first,second in zip(original.samples,larger.samples):
        for key in ['coverage_m2','temperature_k','released_m2','binding_energy_j_m2',
            'internal_energy_j_m2','translational_energy_j_m2','total_wall_energy_debit_j_m2']:
            assert getattr(first,key)==getattr(second,key)
        assert first.extensive['released_number']==first.released_m2*original.problem.patch_area_m2
        for key,value in first.extensive.items():
            assert second.extensive[key]==4*value
        for key in ['binding_energy','internal_energy','translational_energy','total_wall_energy_debit']:
            assert first.extensive[key+'_j']==getattr(first,key+'_j_m2')*original.problem.patch_area_m2
    sample=original.samples[-1]
    with pytest.raises((FrozenInstanceError,AttributeError)):
        sample.temperature_k=100.0
    with pytest.raises(TypeError):
        sample.extensive['released_number']=0
    exported=sample.to_dict()
    exported['extensive']['released_number']=0
    assert sample.extensive['released_number']>0


def test_empty_inventory_remains_stationary_with_exact_zero_channels():
    p=problem(initial_coverage_m2=0)
    history=evolve_feedback(p,duration_s=1e-7,steps=8)
    for sample in history.samples:
        assert sample.coverage_m2==sample.released_m2==0
        assert sample.temperature_k==p.initial_temperature_k
        for key in ['binding_energy_j_m2','internal_energy_j_m2','translational_energy_j_m2',
            'total_wall_energy_debit_j_m2','gas_energy_j_m2','wall_thermal_change_j_m2',
            'adsorbate_reference_energy_j_m2','inventory_residual_m2','total_energy_residual_j_m2']:
            assert getattr(sample,key)==0
        assert all(value==0 for value in sample.extensive.values())


@pytest.mark.parametrize('changes',[
    {'initial_coverage_m2':-1.0},{'initial_coverage_m2':5e20},
    {'initial_temperature_k':True},{'initial_temperature_k':float('nan')},
    {'initial_temperature_k':100.0},{'areal_heat_capacity_j_m2_k':0.0},
    {'areal_heat_capacity_j_m2_k':float('inf')},{'patch_area_m2':0.0},
    {'temperature_min_k':0.0},{'temperature_max_k':100.0},
    {'parameter_status':'sourced'},{'parameter_provenance':''},{'energy_reference':''},
])
def test_invalid_problem_declarations_rejected(changes):
    with pytest.raises(ValueError):
        problem(**changes)


@pytest.mark.parametrize('changes',[
    {'attempt_frequency_s':True},{'attempt_frequency_s':0.0},
    {'activation_energy_ev':-0.1},{'desorption_energy_ev':-0.1},
    {'internal_energy_ev':float('nan')},{'reference_coverage_m2':0.0},
    {'kinetic_provenance':''},
    {'activation_energy_ev':0.35,'coverage_coefficient':1.1},
])
def test_invalid_reservoir_declarations_rejected(changes):
    with pytest.raises(ValueError):
        reservoir(**changes)


@pytest.mark.parametrize('duration,steps',[(0.0,8),(True,8),(float('nan'),8),(1e-7,True),(1e-7,0),(1e-7,1.5),(1e-7,1000001)])
def test_invalid_time_grid_rejected(duration,steps):
    with pytest.raises(ValueError):
        evolve_feedback(problem(),duration_s=duration,steps=steps)


def test_oversized_step_and_insufficient_thermal_domain_reject_without_clipping():
    p=problem()
    with pytest.raises(ValueError):
        evolve_feedback(p,duration_s=1e-6,steps=1)
    with pytest.raises(ValueError):
        evolve_feedback(replace(p,temperature_min_k=899.0),duration_s=1e-7,steps=128)
    assert p.initial_coverage_m2==2e20 and p.initial_temperature_k==900.0


def test_unrepresented_positive_rate_or_feedback_rejected():
    with pytest.raises(UnsupportedFeedbackRangeError):
        feedback_rhs(problem(reservoir(activation_energy_ev=100.0)),coverage_m2=2e20,temperature_k=900.0)
    with pytest.raises(UnsupportedFeedbackRangeError):
        p=problem(initial_coverage_m2=2e-300)
        evolve_feedback(p,duration_s=1e-7,steps=8)


@pytest.mark.parametrize('field',[
    'initial_coverage_m2','activation_energy_ev','desorption_energy_ev',
    'internal_energy_ev','coverage_coefficient',
])
def test_negative_exact_real_cannot_round_to_accepted_zero(field):
    """Nonnegative physical declarations must reject the original negative Real."""
    negative=Fraction(-1,10**1000)
    assert negative<0 and float(negative)==0.0
    with pytest.raises(ValueError):
        if field=='initial_coverage_m2':
            problem(**{field:negative})
        else:
            reservoir(**{field:negative})
