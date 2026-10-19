"""New manufactured controls for finite prescribed-power thermal feedback."""
from dataclasses import FrozenInstanceError,replace
from fractions import Fraction
import math

import pytest

from mhdlab.thermal_power_reference import (
    FeedbackProblem,FeedbackReservoir,PowerSchedule,PowerSegment,
    UnsupportedPowerRangeError,evolve_power,power_rhs,
)

KB=1.380649e-23
QE=1.602176634e-19
DURATION=2e-7


def reservoir(**changes):
    values=dict(name='manufactured_power_parent',state_id='manufactured_fixed_state',
        attempt_frequency_s=5e6,activation_energy_ev=0.0,desorption_energy_ev=0.12,
        internal_energy_ev=0.03,coverage_coefficient=0.6,
        reference_coverage_m2=4e20,max_coverage_m2=4e20,
        kinetic_provenance='Manufactured zero-barrier kinetics; no physical parameter adoption.',
        energy_provenance='Manufactured independent removal and fixed internal reference.')
    values.update(changes)
    return FeedbackReservoir(**values)


def problem(pool=None,**changes):
    values=dict(reservoir=reservoir() if pool is None else pool,
        initial_coverage_m2=2e20,initial_temperature_k=900.0,
        areal_heat_capacity_j_m2_k=0.05,patch_area_m2=2e-4,
        temperature_min_k=600.0,temperature_max_k=1300.0,
        parameter_status='manufactured',parameter_provenance='Manufactured prescribed-power local reference only.',
        energy_reference='Common declared manufactured reference: bound parent -D*N; gas internal I and continuum translation remain separate.')
    values.update(changes)
    return FeedbackProblem(**values)


def schedule(parts=None,*,duration=DURATION,steps=32,**changes):
    if parts is None:
        parts=[(0.0,duration,1e7)]
    segments=tuple(PowerSegment(start,stop,power,steps,
        'Manufactured finite constant segment; no implicit continuation.') for start,stop,power in parts)
    values=dict(duration_s=duration,segments=segments,parameter_status='manufactured',
        parameter_provenance='Manufactured complete finite schedule with explicit breakpoints.')
    values.update(changes)
    return PowerSchedule(**values)


def logical_time(sched,sample):
    segment=sched.segments[sample.segment_index]
    start=Fraction.from_float(segment.start_s)
    width=Fraction.from_float(segment.stop_s)-start
    return start+width*Fraction(sample.segment_step,segment.steps)


def expected_input(sched,sample):
    """Exact integral on the declared canonical binary time coordinates."""
    stop=logical_time(sched,sample)
    total=Fraction(0)
    for segment in sched.segments:
        left=Fraction.from_float(segment.start_s)
        right=min(stop,Fraction.from_float(segment.stop_s))
        if right>left:
            total+=Fraction.from_float(segment.power_w_m2)*(right-left)
    return total


def forced_temperature(p,sched,time):
    """Independent convergent series of the scaled integrating-factor kernel.

    exp(beta*exp(-nu*t))*integral exp(-beta*exp(-nu*s))ds is
    expanded analytically on each constant-power segment. This is independent
    of the joint RK stages and uses no production or oracle helper.
    """
    nu=p.reservoir.attempt_frequency_s
    beta=2*KB*p.initial_coverage_m2/p.areal_heat_capacity_j_m2_k
    removed_fraction=-math.expm1(-nu*time)
    a=(p.reservoir.desorption_energy_ev+p.reservoir.internal_energy_ev)*QE/(2*KB)
    value=p.initial_temperature_k+(p.initial_temperature_k+a)*math.expm1(-beta*removed_fraction)
    terms=[]
    for segment in sched.segments:
        left=segment.start_s;right=min(time,segment.stop_s)
        if right<=left or not segment.power_w_m2:
            continue
        width=right-left
        series=[width]
        coefficient=1.0
        for k in range(1,41):
            coefficient*=-beta/k
            series.append(coefficient*math.exp(-k*nu*left)*(-math.expm1(-k*nu*width))/(k*nu))
        integral=math.exp(beta*math.exp(-nu*time))*math.fsum(series)
        terms.append(segment.power_w_m2/p.areal_heat_capacity_j_m2_k*integral)
    return value+math.fsum(terms)


@pytest.fixture(scope='module')
def empty_history():
    return evolve_power(problem(initial_coverage_m2=0),schedule([(0.0,DURATION,2e7)]))


@pytest.fixture(scope='module')
def constant_histories():
    return {power:evolve_power(problem(),schedule([(0.0,DURATION,power)],steps=128))
            for power in [1e7,6e7]}


@pytest.fixture(scope='module')
def scheduled_history():
    sched=schedule([(0.0,5e-8,6e7),(5e-8,1e-7,0.0),(1e-7,DURATION,3e7)],steps=32)
    return evolve_power(problem(),sched)


def test_empty_inventory_has_exact_prescribed_input_and_heating(empty_history):
    history=empty_history;p=history.problem
    for sample in history.samples:
        energy=float(expected_input(history.schedule,sample))
        assert sample.time_s==float(logical_time(history.schedule,sample))
        assert sample.prescribed_input_energy_j_m2==energy
        assert sample.input_energy_j_m2==pytest.approx(energy,rel=3e-14,abs=0)
        assert sample.temperature_k==pytest.approx(p.initial_temperature_k+energy/p.areal_heat_capacity_j_m2_k,rel=3e-14)
        assert sample.coverage_m2==sample.released_m2==0
        for key in ['binding_energy_j_m2','internal_energy_j_m2','translational_energy_j_m2',
            'total_wall_energy_debit_j_m2','gas_energy_j_m2','adsorbate_reference_energy_j_m2']:
            assert getattr(sample,key)==0
    assert history.samples[-1].input_energy_j_m2>0


def test_constant_power_joint_inventory_temperature_and_channels(constant_histories):
    for power,history in constant_histories.items():
        p=history.problem;nu=p.reservoir.attempt_frequency_s
        d=p.reservoir.desorption_energy_ev*QE;internal=p.reservoir.internal_energy_ev*QE
        rhs=power_rhs(p,coverage_m2=p.initial_coverage_m2,temperature_k=p.initial_temperature_k,power_w_m2=power)
        flux=nu*p.initial_coverage_m2
        assert rhs['release_flux_m2_s']==pytest.approx(flux,rel=3e-14)
        assert rhs['binding_power_j_m2_s']==pytest.approx(d*flux,rel=3e-14)
        assert rhs['internal_power_j_m2_s']==pytest.approx(internal*flux,rel=3e-14)
        assert rhs['translational_power_j_m2_s']==pytest.approx(2*KB*p.initial_temperature_k*flux,rel=3e-14)
        assert rhs['temperature_rate_k_s']==pytest.approx(
            (power-flux*(d+internal+2*KB*p.initial_temperature_k))/p.areal_heat_capacity_j_m2_k,rel=3e-14)
        for sample in history.samples:
            n=p.initial_coverage_m2*math.exp(-nu*sample.time_s)
            emitted=-p.initial_coverage_m2*math.expm1(-nu*sample.time_s)
            temperature=forced_temperature(p,history.schedule,sample.time_s)
            ein=float(expected_input(history.schedule,sample))
            wall=ein-p.areal_heat_capacity_j_m2_k*(temperature-p.initial_temperature_k)
            translation=wall-(d+internal)*emitted
            assert sample.coverage_m2==pytest.approx(n,rel=8e-10)
            assert sample.released_m2==pytest.approx(emitted,rel=8e-10,abs=0)
            assert sample.temperature_k==pytest.approx(temperature,rel=8e-10)
            assert sample.prescribed_input_energy_j_m2==ein
            assert sample.input_energy_j_m2==pytest.approx(ein,rel=3e-14,abs=0)
            assert sample.binding_energy_j_m2==pytest.approx(d*emitted,rel=8e-10,abs=2e-14)
            assert sample.internal_energy_j_m2==pytest.approx(internal*emitted,rel=8e-10,abs=2e-14)
            assert sample.translational_energy_j_m2==pytest.approx(translation,rel=2e-9,abs=3e-14)
            assert sample.total_wall_energy_debit_j_m2==pytest.approx(wall,rel=2e-9,abs=3e-14)
            scale=max(ein,sample.total_wall_energy_debit_j_m2,1e-20)
            assert abs(sample.total_energy_residual_j_m2)<=3e-10*scale
            assert abs(math.fsum([sample.wall_thermal_change_j_m2,
                sample.adsorbate_reference_inventory_change_j_m2,sample.gas_energy_j_m2,-sample.input_energy_j_m2]))<=3e-10*scale
        assert (history.samples[-1].temperature_k>p.initial_temperature_k)==(power==6e7)


def test_discontinuous_schedule_keeps_one_state_and_segment_sided_stages(scheduled_history):
    history=scheduled_history;p=history.problem
    assert len(history.samples)==97
    assert all(right.time_s>left.time_s for left,right in zip(history.samples,history.samples[1:]))
    for boundary,index in [(5e-8,0),(1e-7,1)]:
        matches=[(i,s) for i,s in enumerate(history.samples) if s.time_s==boundary]
        assert len(matches)==1
        i,sample=matches[0]
        assert sample.segment_index==index and sample.segment_step==32
        assert history.samples[i+1].segment_index==index+1 and history.samples[i+1].segment_step==1
    for sample in history.samples:
        assert sample.temperature_k==pytest.approx(forced_temperature(p,history.schedule,sample.time_s),rel=2e-9)
        assert sample.coverage_m2==pytest.approx(p.initial_coverage_m2*math.exp(-p.reservoir.attempt_frequency_s*sample.time_s),rel=1e-9)
        assert sample.prescribed_input_energy_j_m2==float(expected_input(history.schedule,sample))
    assert history.samples[-1].prescribed_input_energy_j_m2==pytest.approx(6.0,rel=2e-15)


def test_exact_initial_balance_retains_positive_gross_transfers():
    """This is an instantaneous balance, not a constant-temperature trajectory."""
    n=5**15;nu=5**14
    p=problem(reservoir(attempt_frequency_s=nu,desorption_energy_ev=0,internal_energy_ev=0),
        initial_coverage_m2=n,initial_temperature_k=1.0,
        areal_heat_capacity_j_m2_k=1.0,temperature_min_k=0.5,temperature_max_k=2.0)
    power=1380649/2**28
    rhs=power_rhs(p,coverage_m2=n,temperature_k=1.0,power_w_m2=power)
    assert rhs['input_power_j_m2_s']==rhs['total_wall_power_j_m2_s']==rhs['translational_power_j_m2_s']==power
    assert rhs['temperature_rate_k_s']==0
    assert rhs['release_flux_m2_s']>0
    assert rhs['binding_power_j_m2_s']==rhs['internal_power_j_m2_s']==0


def test_area_once_and_schedule_immutability(constant_histories):
    history=constant_histories[6e7]
    doubled=evolve_power(replace(history.problem,patch_area_m2=2*history.problem.patch_area_m2),history.schedule)
    for first,second in zip(history.samples,doubled.samples):
        for key in ['coverage_m2','temperature_k','released_m2','input_energy_j_m2','prescribed_input_energy_j_m2',
            'binding_energy_j_m2','internal_energy_j_m2','translational_energy_j_m2','total_wall_energy_debit_j_m2']:
            assert getattr(first,key)==getattr(second,key)
        assert first.extensive['input_energy_j']==first.input_energy_j_m2*history.problem.patch_area_m2
        assert first.extensive['released_number']==first.released_m2*history.problem.patch_area_m2
        assert second.extensive=={key:2*value for key,value in first.extensive.items()}
    assert isinstance(history.schedule.segments,tuple) and isinstance(history.samples,tuple)
    with pytest.raises((FrozenInstanceError,AttributeError)):
        history.schedule.segments[0].power_w_m2=0.0
    with pytest.raises(TypeError):
        history.samples[-1].extensive['input_energy_j']=0.0


def test_positive_barrier_is_explicitly_outside_the_qualified_operator():
    p=problem(reservoir(activation_energy_ev=0.1))
    with pytest.raises(ValueError,match='zero|activation|barrier'):
        power_rhs(p,coverage_m2=2e20,temperature_k=900.0,power_w_m2=1e7)
    with pytest.raises(ValueError,match='zero|activation|barrier'):
        evolve_power(p,schedule())


@pytest.mark.parametrize('field,value',[
    ('power_w_m2',-1.0),('power_w_m2',Fraction(-1,10**1000)),
    ('power_w_m2',Fraction(1,10**1000)),('power_w_m2',True),
    ('power_w_m2',float('nan')),('power_w_m2',float('inf')),
    ('start_s',True),('stop_s',float('inf')),('stop_s',Fraction(1,10)),
    ('steps',True),('steps',0),('steps',1.5),('power_provenance',''),
])
def test_invalid_new_power_segment_declarations_rejected(field,value):
    fields=dict(start_s=0.0,stop_s=DURATION,power_w_m2=1e7,steps=32,
        power_provenance='Manufactured segment declaration.')
    fields[field]=value
    with pytest.raises(ValueError):
        PowerSegment(**fields)


@pytest.mark.parametrize('parts',[
    [],[(1e-8,DURATION,1e7)],[(0.0,1e-7,1e7)],
    [(0.0,1e-7,1e7),(math.nextafter(1e-7,math.inf),DURATION,1e7)],
    [(0.0,1e-7,1e7),(math.nextafter(1e-7,0.0),DURATION,1e7)],
    [(1e-7,DURATION,1e7),(0.0,1e-7,1e7)],
])
def test_complete_nonoverlapping_schedule_required(parts):
    with pytest.raises(ValueError):
        schedule(parts)


def test_exact_boundary_gap_cannot_be_rounded_away():
    exact=Fraction.from_float(1e-7)+Fraction(1,10**1000)
    assert float(exact)==1e-7
    with pytest.raises(ValueError):
        schedule([(0.0,1e-7,1e7),(exact,DURATION,1e7)])


@pytest.mark.parametrize('changes',[
    {'duration_s':True},{'duration_s':0.0},{'duration_s':float('nan')},
    {'parameter_status':'sourced'},{'parameter_provenance':''},{'segments':'not a segment sequence'},
])
def test_invalid_new_schedule_metadata_rejected(changes):
    with pytest.raises(ValueError):
        schedule(**changes)


@pytest.mark.parametrize('power',[True,-1.0,Fraction(-1,10**1000),float('nan')])
def test_rhs_rejects_invalid_new_power_scalar(power):
    with pytest.raises(ValueError):
        power_rhs(problem(),coverage_m2=2e20,temperature_k=900.0,power_w_m2=power)


def test_unsupported_stages_and_positive_gross_input_underflow_reject():
    with pytest.raises(UnsupportedPowerRangeError):
        evolve_power(problem(),schedule([(0.0,1e-6,1e7)],duration=1e-6,steps=1))
    with pytest.raises(UnsupportedPowerRangeError):
        evolve_power(problem(temperature_max_k=901.0),schedule([(0.0,DURATION,6e7)]))
    with pytest.raises(UnsupportedPowerRangeError):
        evolve_power(problem(initial_coverage_m2=0),schedule([(0.0,DURATION,math.ulp(0.0))]))
