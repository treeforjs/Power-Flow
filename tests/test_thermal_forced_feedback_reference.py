"""Independent structural controls for the manufactured forced Arrhenius branch."""
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
import math

import pytest

from mhdlab.thermal_forced_feedback_reference import (
    FeedbackProblem,
    FeedbackReservoir,
    PowerSchedule,
    PowerSegment,
    UnsupportedForcedFeedbackRangeError,
    evolve_forced_feedback,
    forced_feedback_rhs,
)


KB = Fraction("1.380649e-23")
QE = Fraction("1.602176634e-19")
TINY_NEGATIVE = Fraction(-1, 10**1000)
TINY_POSITIVE = -TINY_NEGATIVE


def pool(**changes):
    parameters = dict(
        name="manufactured_pool", state_id="manufactured_parent",
        attempt_frequency_s=1e8, activation_energy_ev=0.35,
        desorption_energy_ev=0.12, internal_energy_ev=0.03,
        coverage_coefficient=0.6, reference_coverage_m2=4e20, max_coverage_m2=4e20,
        kinetic_provenance="Explicit manufactured rate parameters for new controls.",
        energy_provenance="Independent manufactured removal and internal costs.",
    )
    parameters.update(changes)
    return FeedbackReservoir(**parameters)


def problem(reservoir=None, **changes):
    parameters = dict(
        reservoir=pool() if reservoir is None else reservoir,
        initial_coverage_m2=2e20, initial_temperature_k=900.0,
        areal_heat_capacity_j_m2_k=0.05, patch_area_m2=2e-4,
        temperature_min_k=600.0, temperature_max_k=1300.0,
        parameter_status="manufactured",
        parameter_provenance="Manufactured single-patch verification only.",
        energy_reference="Thermal change from initial T; adsorbate reference -D*N.",
    )
    parameters.update(changes)
    return FeedbackProblem(**parameters)


def schedule(powers=(6e7, 0.0, 8e7), *, duration_per_segment=2**-28, steps=(4, 4, 4)):
    segments = tuple(
        PowerSegment(index*duration_per_segment, (index+1)*duration_per_segment,
                     power, count, "Manufactured constant power on this exact support.")
        for index, (power, count) in enumerate(zip(powers, steps))
    )
    return PowerSchedule(len(segments)*duration_per_segment, segments, "manufactured",
                         "Explicit complete manufactured schedule.")


def rhs(value=None, **changes):
    parameters = dict(coverage_m2=2e20, temperature_k=900.0, power_w_m2=6e7)
    parameters.update(changes)
    return forced_feedback_rhs(problem() if value is None else value, **parameters)


@pytest.fixture(scope="module")
def short_history():
    return evolve_forced_feedback(problem(), schedule())


def test_temperature_changes_the_arrhenius_rate_at_fixed_state():
    cold = rhs(temperature_k=700.0)
    hot = rhs(temperature_k=1000.0)
    barrier = 0.35*(1-0.6*2e20/4e20)
    expected_ratio = math.exp(barrier*float(QE/KB)*(1/700-1/1000))
    assert hot["release_flux_m2_s"] > 3*cold["release_flux_m2_s"] > 0
    assert hot["release_flux_m2_s"]/cold["release_flux_m2_s"] == pytest.approx(expected_ratio, rel=5e-13)


def test_removal_and_internal_costs_change_cooling_without_changing_kinetics():
    low = rhs(problem(pool(desorption_energy_ev=0.0, internal_energy_ev=0.0)))
    high = rhs(problem(pool(desorption_energy_ev=0.4, internal_energy_ev=0.2)))
    assert low["release_flux_m2_s"] == high["release_flux_m2_s"]
    assert low["coverage_rate_m2_s"] == high["coverage_rate_m2_s"]
    assert high["temperature_rate_k_s"] < low["temperature_rate_k_s"]
    difference = high["total_wall_power_j_m2_s"]-low["total_wall_power_j_m2_s"]
    assert difference == pytest.approx(low["release_flux_m2_s"]*0.6*float(QE), rel=3e-14)


def test_positive_kinetic_barrier_does_not_create_a_removal_energy_debit():
    values = rhs(problem(pool(desorption_energy_ev=0.0, internal_energy_ev=0.0)))
    assert values["release_flux_m2_s"] > 0
    assert values["binding_power_j_m2_s"] == values["internal_power_j_m2_s"] == 0.0
    assert values["total_wall_power_j_m2_s"] == values["translational_power_j_m2_s"] > 0


def test_reference_coverage_is_not_reset_to_initial_inventory():
    first = problem(initial_coverage_m2=1e20)
    second = problem(initial_coverage_m2=3e20)
    assert first.reservoir.reference_coverage_m2 == second.reservoir.reference_coverage_m2 == 4e20
    assert rhs(first, coverage_m2=5e19) == rhs(second, coverage_m2=5e19)


def test_zero_effective_barrier_at_the_declared_domain_endpoint_is_supported():
    value = problem(pool(coverage_coefficient=1.0), initial_coverage_m2=4e20)
    values = rhs(value, coverage_m2=4e20)
    assert values["release_flux_m2_s"] == pytest.approx(1e8*4e20, rel=2e-15)


def test_positive_gross_channels_allow_both_heating_and_cooling():
    warming = rhs(power_w_m2=6e7)
    cooling = rhs(power_w_m2=1e6)
    assert warming["temperature_rate_k_s"] > 0 > cooling["temperature_rate_k_s"]
    assert warming["release_flux_m2_s"] == cooling["release_flux_m2_s"]
    for values in (warming, cooling):
        assert values["input_power_j_m2_s"] > 0
        assert values["total_wall_power_j_m2_s"] > 0
        assert values["translational_power_j_m2_s"] > 0


def test_population_scaled_rate_survives_an_underflowed_binary_exponential():
    activation = float(Fraction(800)*KB*1000/QE)
    value = problem(pool(activation_energy_ev=activation, attempt_frequency_s=1e308,
                         desorption_energy_ev=0.0, internal_energy_ev=0.0,
                         coverage_coefficient=0.0, reference_coverage_m2=1.0, max_coverage_m2=1.0),
                    initial_coverage_m2=1.0, initial_temperature_k=1000.0)
    beta = float(Fraction.from_float(activation)*QE/(KB*1000))
    assert math.exp(-beta) == 0.0
    expected = 1e308*math.exp(-beta/2)*math.exp(-beta/2)
    values = rhs(value, coverage_m2=1.0, temperature_k=1000.0, power_w_m2=0.0)
    assert values["release_flux_m2_s"] == pytest.approx(expected, rel=1e-12)
    assert values["release_flux_m2_s"] > 0 and values["temperature_rate_k_s"] < 0
    assert values["translational_power_j_m2_s"] > 0


def test_empty_inventory_has_exact_logical_input_and_analytic_heating():
    value = problem(initial_coverage_m2=0.0, areal_heat_capacity_j_m2_k=0.125)
    control = schedule((2**20, 0.0, 2**19), duration_per_segment=2**-20, steps=(3, 2, 4))
    history = evolve_forced_feedback(value, control)
    for sample in history.samples:
        segment = control.segments[sample.segment_index]
        local = (Fraction.from_float(segment.stop_s)-Fraction.from_float(segment.start_s))*sample.segment_step/segment.steps
        emitted_input = sum(Fraction.from_float(old.power_w_m2)*(Fraction.from_float(old.stop_s)-Fraction.from_float(old.start_s))
                            for old in control.segments[:sample.segment_index])
        emitted_input += Fraction.from_float(segment.power_w_m2)*local
        expected_temperature = Fraction(900)+emitted_input/Fraction(1, 8)
        assert sample.temperature_k == pytest.approx(float(expected_temperature), rel=0, abs=5e-13)
        assert sample.input_energy_j_m2 == float(emitted_input)
        assert sample.prescribed_input_energy_j_m2 == float(emitted_input)
        assert sample.coverage_m2 == sample.released_m2 == 0.0
        assert sample.binding_energy_j_m2 == sample.internal_energy_j_m2 == 0.0
        assert sample.translational_energy_j_m2 == sample.total_wall_energy_debit_j_m2 == 0.0
    assert history.samples[-1].temperature_k == pytest.approx(912.0, rel=0, abs=5e-13)


def test_power_jumps_keep_one_boundary_state_and_independent_gross_channels(short_history):
    history = short_history
    assert len(history.samples) == 13
    for segment_index in (0, 1):
        boundary = history.samples[(segment_index+1)*4]
        following = history.samples[(segment_index+1)*4+1]
        assert boundary.segment_index == segment_index and boundary.segment_step == 4
        assert following.segment_index == segment_index+1 and following.segment_step == 1
        assert sum(sample.time_s == boundary.time_s for sample in history.samples) == 1
    assert history.samples[4].temperature_k > history.samples[0].temperature_k
    assert history.samples[8].temperature_k < history.samples[4].temperature_k
    assert history.samples[12].temperature_k > history.samples[8].temperature_k
    assert history.samples[8].input_energy_j_m2 == history.samples[4].input_energy_j_m2
    for left, right in zip(history.samples, history.samples[1:]):
        assert 0 < right.coverage_m2 < left.coverage_m2
        assert right.released_m2 > left.released_m2
        assert right.total_wall_energy_debit_j_m2 > left.total_wall_energy_debit_j_m2
    for sample in history.samples[1:]:
        assert sample.binding_energy_j_m2 == pytest.approx(0.12*float(QE)*sample.released_m2, rel=3e-14)
        assert sample.internal_energy_j_m2 == pytest.approx(0.03*float(QE)*sample.released_m2, rel=3e-14)
        assert sample.translational_energy_j_m2 > 0
    offset = float(Fraction.from_float(0.12+0.03)*QE/(2*KB))
    off_start, off_stop = history.samples[4], history.samples[8]
    off_invariant = math.log((off_stop.temperature_k+offset)/(off_start.temperature_k+offset))
    off_invariant -= 2*float(KB)*(off_stop.coverage_m2-off_start.coverage_m2)/0.05
    assert abs(off_invariant) < 2e-9
    on_start, on_stop = history.samples[0], history.samples[4]
    on_invariant = math.log((on_stop.temperature_k+offset)/(on_start.temperature_k+offset))
    on_invariant -= 2*float(KB)*(on_stop.coverage_m2-on_start.coverage_m2)/0.05
    assert abs(on_invariant) > 1e-4


def test_area_is_applied_once_and_returned_history_is_immutable(short_history):
    scaled = evolve_forced_feedback(problem(patch_area_m2=6e-4), schedule())
    for base, larger in zip(short_history.samples, scaled.samples):
        base_dict, larger_dict = base.to_dict(), larger.to_dict()
        base_dict.pop("extensive"); larger_dict.pop("extensive")
        assert base_dict == larger_dict
        for name, number in base.extensive.items():
            assert larger.extensive[name] == pytest.approx(3*number, rel=4e-15, abs=0)
    with pytest.raises(FrozenInstanceError):
        short_history.samples[0].temperature_k = 1.0
    with pytest.raises(TypeError):
        short_history.samples[-1].extensive["released_number"] = 0
    detached = short_history.to_dict()
    detached["samples"][-1]["extensive"]["released_number"] = -1
    assert short_history.samples[-1].extensive["released_number"] > 0


@pytest.mark.parametrize("changes", [
    {"coverage_m2": True}, {"coverage_m2": TINY_NEGATIVE}, {"coverage_m2": 4.1e20},
    {"temperature_k": True}, {"temperature_k": math.nan}, {"temperature_k": 599.0},
    {"power_w_m2": True}, {"power_w_m2": TINY_NEGATIVE}, {"power_w_m2": math.inf},
    {"power_w_m2": TINY_POSITIVE},
])
def test_new_rhs_rejects_invalid_or_unrepresented_inputs(changes):
    with pytest.raises(ValueError):
        rhs(**changes)


@pytest.mark.parametrize("kind", ["zero_base_barrier", "negative_effective_barrier", "initial_inventory", "status"])
def test_new_entry_revalidates_forged_parameter_data(kind):
    value = problem()
    if kind == "zero_base_barrier":
        value = replace(value, reservoir=pool(activation_energy_ev=0.0))
    elif kind == "negative_effective_barrier":
        object.__setattr__(value.reservoir, "coverage_coefficient", 2.0)
    elif kind == "initial_inventory":
        object.__setattr__(value, "initial_coverage_m2", 8e20)
    else:
        object.__setattr__(value, "parameter_status", "sourced")
    with pytest.raises(ValueError):
        evolve_forced_feedback(value, schedule())


@pytest.mark.parametrize("kind", ["negative_exact_power", "boolean_steps", "incomplete_support"])
def test_new_entry_revalidates_forged_schedule_data(kind):
    control = schedule()
    if kind == "negative_exact_power":
        object.__setattr__(control.segments[0], "power_w_m2", TINY_NEGATIVE)
    elif kind == "boolean_steps":
        object.__setattr__(control.segments[0], "steps", True)
    else:
        object.__setattr__(control, "duration_s", control.duration_s*2)
    with pytest.raises(ValueError):
        evolve_forced_feedback(problem(), control)


@pytest.mark.parametrize("activation_energy_ev", [1000.0, 1e-300])
def test_unrepresentable_arrhenius_hazard_is_rejected(activation_energy_ev):
    with pytest.raises(UnsupportedForcedFeedbackRangeError):
        rhs(problem(pool(activation_energy_ev=activation_energy_ev)))


def test_positive_transfers_that_cannot_change_the_state_are_rejected():
    value = problem(pool(activation_energy_ev=0.01, attempt_frequency_s=1.0))
    control = schedule((1e7,), duration_per_segment=math.ulp(0.0), steps=(1,))
    with pytest.raises(UnsupportedForcedFeedbackRangeError):
        evolve_forced_feedback(value, control)


@pytest.mark.parametrize("kind", ["temperature_stage", "inventory_stage"])
def test_oversized_stages_are_rejected_without_clipping_or_subdivision(kind):
    value = problem(pool(attempt_frequency_s=1e-12))
    power = 1e20
    if kind == "inventory_stage":
        value = problem(pool(activation_energy_ev=0.001, attempt_frequency_s=1e14))
        power = 1.0
    control = schedule((power,), duration_per_segment=1.0, steps=(1,))
    with pytest.raises(UnsupportedForcedFeedbackRangeError):
        evolve_forced_feedback(value, control)
