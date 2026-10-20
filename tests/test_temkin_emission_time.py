"""Independent conditional emission-time verification for the specified law.

These are manufactured rate cases.  Differential-equation oracles neither call
the depletion inverter nor derive parameters for a physical surface species.
"""

from dataclasses import FrozenInstanceError
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from mhdlab.constants import KB, QE
from mhdlab.temkin_emission_time import (
    FrozenTemkinEmission, NoTemkinEmissionError, TemkinEmissionTimeRangeError,
    TIME_IDENTITY,
)
from mhdlab.temkin_release import LAW_IDENTITY, TemkinReservoir, UnsupportedTemkinRangeError
from mhdlab.thermal_release import BindingReservoir


KB_EV_K = KB / QE


def reservoir(**overrides):
    parameters = dict(
        name="manufactured-timing-pool", state_id="manufactured:ground",
        attempt_frequency_s=8e3, activation_energy_ev=.2,
        desorption_energy_ev=.13, internal_energy_ev=.07,
        coverage_coefficient=.6, reference_coverage_m2=1e18,
        max_coverage_m2=1e18,
        kinetic_provenance="manufactured rate for emission-time verification",
        energy_provenance="manufactured independent energy cost",
    )
    parameters.update(overrides)
    return TemkinReservoir(**parameters)


def oracle(pool, initial, temperature, duration):
    """Independent normalized log-survivor ODE and emitted-population CDF."""
    activation = pool.activation_energy_ev * QE / (KB * temperature)
    coefficient = activation * pool.coverage_coefficient * initial / pool.reference_coverage_m2
    initial_rate = pool.attempt_frequency_s * math.exp(-activation + coefficient)
    final_clock = initial_rate * duration
    solution = solve_ivp(
        lambda _, value: [math.exp(coefficient * math.expm1(-float(value[0])))],
        (0., final_clock), [0.], method="DOP853", dense_output=True,
        rtol=2e-13, atol=2e-15,
    )
    assert solution.success
    total_fraction = -math.expm1(-float(solution.y[0, -1]))
    assert total_fraction > 0

    def cdf(offset):
        value = float(solution.sol(initial_rate * offset)[0])
        return -math.expm1(-value) / total_fraction

    return cdf, total_fraction, initial_rate, coefficient


def oracle_quantile(probability, total_fraction, initial_rate, coefficient):
    """Integrate the time ODE at a prescribed population, without inversion.

    The target log depletion follows directly from the conditional released
    population.  Integrating time along that known depletion is an independent
    DOP853 initial-value problem on a fixed [0,1] coordinate, including tiny p.
    """
    target = -math.log1p(-probability * total_fraction)
    solution = solve_ivp(
        lambda location, _: [math.exp(coefficient * (-math.expm1(-target * location)))],
        (0., 1.), [0.], method="DOP853", rtol=2e-13, atol=2e-15,
    )
    assert solution.success
    return target * float(solution.y[0, -1]) / initial_rate


@pytest.mark.parametrize("coefficient", [.5, 3., 20.])
def test_nonlinear_conditional_cdf_matches_independent_time_ode(coefficient):
    temperature, initial = 700., 7.5e17
    alpha = .8
    activation = coefficient / (alpha * initial / 1e18)
    pool = reservoir(activation_energy_ev=activation * KB_EV_K * temperature,
                     coverage_coefficient=alpha)
    initial_rate = pool.attempt_frequency_s * math.exp(-activation + coefficient)
    duration = 12. / initial_rate
    emission = FrozenTemkinEmission(pool, initial, temperature, duration)
    expected, _, _, _ = oracle(pool, initial, temperature, duration)
    offsets = duration * np.array([0., 1e-8, .01, .1, .6, .99, 1.])
    values = emission.cdf(offsets)
    np.testing.assert_allclose(values, [expected(float(offset)) for offset in offsets], rtol=2e-10, atol=0.)
    assert values[0] == 0. and values[-1] == 1.
    assert np.all(np.diff(values) > 0)


@pytest.mark.parametrize("coefficient", [.5, 3., 20.])
def test_nonlinear_quantiles_match_independent_population_time_ode(coefficient):
    temperature, initial = 700., 7.5e17
    activation = coefficient / (.8 * initial / 1e18)
    pool = reservoir(activation_energy_ev=activation * KB_EV_K * temperature,
                     coverage_coefficient=.8)
    duration = 12. / (pool.attempt_frequency_s * math.exp(-activation + coefficient))
    emission = FrozenTemkinEmission(pool, initial, temperature, duration)
    _, total, initial_rate, feedback = oracle(pool, initial, temperature, duration)
    probabilities = np.array([1e-12, .001, .1, .5, .9, .999, 1 - 1e-12])
    expected = np.array([oracle_quantile(float(value), total, initial_rate, feedback) for value in probabilities])
    actual = emission.quantile(probabilities)
    np.testing.assert_allclose(actual / duration, expected / duration, rtol=3e-10, atol=0.)
    np.testing.assert_allclose(emission.cdf(actual), probabilities, rtol=3e-10, atol=0.)
    assert np.all(np.diff(actual) > 0)
    assert np.all(actual <= probabilities * duration * (1 + 2e-12))


@pytest.mark.parametrize("hazard", [1e-8, .2, 4., 80.])
@pytest.mark.parametrize("zero_energy", [False, True])
def test_constant_rate_limits_match_conditional_truncated_exponential(hazard, zero_energy):
    temperature, duration = 600., .003
    activation = 0. if zero_energy else 2.
    pool = reservoir(activation_energy_ev=activation * KB_EV_K * temperature,
                     attempt_frequency_s=hazard * math.exp(activation) / duration,
                     coverage_coefficient=4. if zero_energy else 0.)
    emission = FrozenTemkinEmission(pool, 7e17, temperature, duration)
    offsets = duration * np.array([0., 1e-10, .1, .5, .9, 1.])
    total_fraction = -math.expm1(-hazard)
    expected = -np.expm1(-hazard * offsets / duration) / total_fraction
    np.testing.assert_allclose(emission.cdf(offsets), expected, rtol=2e-11, atol=0.)
    probabilities = np.array([1e-12, .01, .5, .99, 1 - 1e-12])
    expected_times = -np.log1p(-probabilities * total_fraction) / (hazard / duration)
    np.testing.assert_allclose(emission.quantile(probabilities), expected_times, rtol=2e-10, atol=0.)


def test_cdf_is_concave_for_monotonically_decreasing_release_rate():
    pool = reservoir()
    duration = .005
    emission = FrozenTemkinEmission(pool, 9e17, 800., duration)
    fractions = np.linspace(0., 1., 31)
    cdf = emission.cdf(fractions * duration)
    assert np.all(cdf >= fractions - 2e-13)
    assert np.all(np.diff(cdf, n=2) <= 2e-12)
    assert np.min(np.diff(cdf, n=2)) < -1e-4


def test_exact_supported_feedback_endpoint_matches_independent_cdf_and_quantile():
    # These float inputs give Ed*QE/(KB*T)=1024 exactly.  Alpha=.5 and
    # Ninitial=Nref therefore place c at the numerical boundary 512 exactly.
    pool = reservoir(activation_energy_ev=1024 * KB, coverage_coefficient=.5,
                     attempt_frequency_s=1e5)
    initial, temperature = 1e18, QE
    duration = 20. / (pool.attempt_frequency_s * math.exp(-512.))
    emission = FrozenTemkinEmission(pool, initial, temperature, duration)
    expected_cdf, total, rate, feedback = oracle(pool, initial, temperature, duration)
    offsets = duration * np.array([.001, .1, .7])
    np.testing.assert_allclose(emission.cdf(offsets), [expected_cdf(float(value)) for value in offsets],
                               rtol=5e-10, atol=0.)
    probabilities = np.array([.01, .5, .99])
    expected_times = [oracle_quantile(float(value), total, rate, feedback) for value in probabilities]
    np.testing.assert_allclose(emission.quantile(probabilities) / duration,
                               np.array(expected_times) / duration, rtol=5e-10, atol=0.)


def test_feedback_outside_core_qualified_range_is_not_silently_changed():
    pool = reservoir(activation_energy_ev=513 * KB_EV_K, coverage_coefficient=1.,
                     reference_coverage_m2=1., max_coverage_m2=1., attempt_frequency_s=1.)
    with pytest.raises(UnsupportedTemkinRangeError):
        FrozenTemkinEmission(pool, 1., 1., 1.)


def test_fully_depleted_stiff_interval_retains_early_nonlinear_event_times():
    pool = reservoir(activation_energy_ev=KB_EV_K, coverage_coefficient=.5,
                     attempt_frequency_s=1e250, reference_coverage_m2=1e20,
                     max_coverage_m2=1e20)
    emission = FrozenTemkinEmission(pool, 1e20, 1., 1e250)
    assert emission.remaining_coverage_m2 == 0.
    assert emission.released_m2 == 1e20
    probabilities = np.array([1e-10, .1, .5, .99, 1 - 1e-12])
    rate = pool.attempt_frequency_s * math.exp(-.5)
    expected = np.array([oracle_quantile(float(value), 1., rate, .5) for value in probabilities])
    times = emission.quantile(probabilities)
    np.testing.assert_allclose(times, expected, rtol=2e-10, atol=0.)
    assert np.all(times / emission.dt_s == 0.)
    np.testing.assert_allclose(emission.cdf(times), probabilities, rtol=2e-10, atol=0.)


def test_underflowed_hazard_with_representable_release_retains_uniform_time_limit():
    pool = reservoir(activation_energy_ev=1000 * KB_EV_K, coverage_coefficient=.001,
                     attempt_frequency_s=1., reference_coverage_m2=1e300,
                     max_coverage_m2=1e300, desorption_energy_ev=0., internal_energy_ev=0.)
    emission = FrozenTemkinEmission(pool, 1e300, 1., 1.)
    assert emission.released_m2 > 0
    assert emission.remaining_coverage_m2 == 1e300
    probabilities = np.array([1e-12, .01, .5, .99, 1 - 1e-12])
    np.testing.assert_allclose(emission.cdf(probabilities), probabilities, rtol=3e-11, atol=0.)
    np.testing.assert_allclose(emission.quantile(probabilities), probabilities, rtol=3e-11, atol=0.)


def test_saturated_supported_endpoint_retains_tiny_and_late_open_quantiles():
    pool = reservoir(activation_energy_ev=1024 * KB, coverage_coefficient=.5,
                     attempt_frequency_s=math.exp(512.))
    emission = FrozenTemkinEmission(pool, 1e18, QE, 1e300)
    assert emission.remaining_coverage_m2 == 0.
    probability = 1e-300
    initial_rate = pool.attempt_frequency_s * math.exp(-512.)
    assert emission.quantile(probability) == pytest.approx(probability / initial_rate, rel=3e-11, abs=0.)
    # Here N0-p*N0 rounds to N0.  Timing must retain the small released
    # probability instead of recovering it by subtracting populations.
    assert 1e18 - probability * 1e18 == 1e18
    late_probability = np.nextafter(1., 0.)
    expected = oracle_quantile(late_probability, 1., initial_rate, 512.)
    actual = emission.quantile(late_probability)
    assert actual == pytest.approx(expected, rel=8e-10, abs=0.)
    assert actual < 1e-60 * emission.dt_s


@pytest.mark.parametrize("scale", [1e-200, 1e200])
def test_rescaling_clock_preserves_cdf_and_scales_quantiles(scale):
    pool = reservoir()
    duration = .003
    first = FrozenTemkinEmission(pool, 7e17, 700., duration)
    scaled_pool = reservoir(attempt_frequency_s=pool.attempt_frequency_s / scale)
    second = FrozenTemkinEmission(scaled_pool, 7e17, 700., duration * scale)
    offsets = duration * np.array([1e-7, .01, .4, .9])
    np.testing.assert_allclose(second.cdf(offsets * scale), first.cdf(offsets), rtol=3e-11, atol=0.)
    probabilities = np.array([1e-10, .01, .5, .99, 1 - 1e-10])
    np.testing.assert_allclose(second.quantile(probabilities) / scale, first.quantile(probabilities), rtol=3e-11, atol=0.)


def test_subfloat_quantile_has_explicit_error_and_finite_log_time():
    pool = reservoir(activation_energy_ev=KB_EV_K, coverage_coefficient=.5,
                     attempt_frequency_s=1e300, reference_coverage_m2=1., max_coverage_m2=1.)
    emission = FrozenTemkinEmission(pool, 1., 1., 1.)
    probability = np.nextafter(0., 1.)
    expected_log_time = math.log(probability) - math.log(pool.attempt_frequency_s) + .5
    assert emission.log_quantile(probability) == pytest.approx(expected_log_time, rel=0., abs=3e-11)
    with pytest.raises(TemkinEmissionTimeRangeError):
        emission.quantile(probability)


def test_law_and_time_identities_cached_trial_and_results_are_immutable():
    pool = reservoir()
    emission = FrozenTemkinEmission(pool, 6e17, 700., .003)
    assert emission.law_identity == LAW_IDENTITY
    assert emission.time_identity == TIME_IDENTITY == "temkin_conditional_emitted_population_frozen_temperature_v1"
    assert emission.trial.law_identity == LAW_IDENTITY
    assert emission.released_m2 == emission.trial.released_m2[0, 0]
    assert emission.remaining_coverage_m2 == emission.trial.trial_inventory.coverage_m2[0, 0]
    offsets = np.array([.0001, .001, .002])
    probabilities = np.array([.1, .5, .9])
    saved_offsets, saved_probabilities = offsets.copy(), probabilities.copy()
    for method, inputs in ((emission.cdf, offsets), (emission.quantile, probabilities),
                           (emission.log_quantile, probabilities)):
        first = method(inputs)
        np.testing.assert_array_equal(first, method(inputs))
        with pytest.raises(ValueError):
            first.setflags(write=True)
        assert isinstance(method(float(inputs[0])), float)
    np.testing.assert_array_equal(offsets, saved_offsets)
    np.testing.assert_array_equal(probabilities, saved_probabilities)
    with pytest.raises(FrozenInstanceError):
        emission.dt_s = 1.


@pytest.mark.parametrize("initial,temperature,duration", [
    (-1., 700., .003), (1.1e18, 700., .003), (np.nan, 700., .003),
    (np.inf, 700., .003), (True, 700., .003), ("1e17", 700., .003),
    (1e17 + 0j, 700., .003), (1e17, 0., .003), (1e17, True, .003),
    (1e17, "700", .003), (1e17, 700. + 0j, .003),
    (1e17, 700., 0.), (1e17, 700., -.003), (1e17, 700., True),
    (1e17, 700., ".003"), (1e17, 700., np.inf),
])
def test_constructor_rejects_invalid_real_domain(initial, temperature, duration):
    with pytest.raises(ValueError):
        FrozenTemkinEmission(reservoir(), initial, temperature, duration)


def test_constructor_rejects_constant_law_records_and_absent_representable_emission():
    constant = BindingReservoir("constant", "manufactured:ground", 1., 0., 0.)
    with pytest.raises(ValueError):
        FrozenTemkinEmission(constant, 1., 1., 1.)
    with pytest.raises(NoTemkinEmissionError):
        FrozenTemkinEmission(reservoir(), 0., 700., .003)
    pool = reservoir(activation_energy_ev=2000 * KB_EV_K, coverage_coefficient=.001,
                     attempt_frequency_s=1., reference_coverage_m2=1., max_coverage_m2=1.)
    with pytest.raises(NoTemkinEmissionError):
        FrozenTemkinEmission(pool, 1., 1., 1.)


@pytest.mark.parametrize("value", [
    -.001, .004, np.nan, np.inf, True, "0.001", .001 + 0j,
    [], [[.001]], [False, .001], [".001"], [np.nan],
    np.ma.array([.001], mask=[True]),
])
def test_cdf_rejects_outside_times_and_ambiguous_types(value):
    emission = FrozenTemkinEmission(reservoir(), 6e17, 700., .003)
    with pytest.raises(ValueError):
        emission.cdf(value)


@pytest.mark.parametrize("value", [
    0., 1., -.1, 1.1, np.nan, np.inf, True, ".5", .5 + 0j,
    [], [[.5]], [True, .5], [".5"], [0., .5], [.5, 1.],
    np.ma.array([.5], mask=[True]),
])
def test_quantile_methods_require_strict_open_real_probabilities(value):
    emission = FrozenTemkinEmission(reservoir(), 6e17, 700., .003)
    for method in (emission.quantile, emission.log_quantile):
        with pytest.raises(ValueError):
            method(value)
