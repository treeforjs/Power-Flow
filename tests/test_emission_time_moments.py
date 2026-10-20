"""Independent conditional-time integrals and extreme actual-second scales."""

from dataclasses import FrozenInstanceError, replace
import math

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from mhdlab.constants import KB, QE
from mhdlab.emission_time_moments import (
    MAX_QUADRATURE_ORDER, QUADRATURE_IDENTITY, TimeMomentRangeError,
    compare_conditional_time_moments, conditional_time_moments,
    gauss_legendre_probability_rule, refine_conditional_time_moments,
)
from mhdlab.temkin_emission_time import FrozenTemkinEmission, TIME_IDENTITY
from mhdlab.temkin_release import LAW_IDENTITY, TemkinReservoir


def emission(*, rate=2., duration=.75, coefficient=0., barrier=0., coverage=1e18,
             reference=1e18, maximum=1e18, temperature=600.):
    reservoir = TemkinReservoir("manufactured", "M", rate, barrier, .5, .1, coefficient,
                               reference, maximum, "manufactured kinetics", "manufactured energy convention")
    return FrozenTemkinEmission(reservoir, coverage, temperature, duration)


@pytest.mark.parametrize("order", [2, 4, 16])
def test_uniform_conditional_limit_has_analytic_moments_and_immutable_identity(order):
    distribution = emission(rate=1e-100, duration=3.)
    assert distribution.uniform_limit
    result = conditional_time_moments(distribution, quadrature_order=order)
    assert result.mean_offset_s == pytest.approx(1.5, rel=2e-14)
    assert result.second_raw_moment_s2 == pytest.approx(3., rel=2e-14)
    assert result.variance_offset_s2 == pytest.approx(.75, rel=2e-14)
    assert result.log_mean_offset_s == pytest.approx(math.log(1.5), abs=2e-14)
    assert result.log_second_raw_moment_s2 == pytest.approx(math.log(3), abs=2e-14)
    assert result.log_variance_offset_s2 == pytest.approx(math.log(.75), abs=2e-14)
    assert result.law_identity == LAW_IDENTITY
    assert result.time_identity == TIME_IDENTITY
    assert result.quadrature_identity == QUADRATURE_IDENTITY
    assert result.quadrature_order == order
    with pytest.raises(FrozenInstanceError):
        result.mean_offset_s = 0


def test_constant_hazard_truncated_exponential_matches_independent_density_quadrature():
    rate, duration = 2., .75
    distribution = emission(rate=rate, duration=duration)
    normalization = -math.expm1(-rate * duration)
    density = lambda time: rate * math.exp(-rate * time) / normalization
    mean = quad(lambda time: time * density(time), 0, duration, epsabs=1e-14, epsrel=1e-13)[0]
    second = quad(lambda time: time ** 2 * density(time), 0, duration, epsabs=1e-14, epsrel=1e-13)[0]
    result = conditional_time_moments(distribution, quadrature_order=64)
    assert result.mean_offset_s == pytest.approx(mean, rel=3e-14)
    assert result.second_raw_moment_s2 == pytest.approx(second, rel=3e-14)
    assert result.variance_offset_s2 == pytest.approx(second - mean ** 2, rel=5e-14)


def test_nonlinear_time_moments_match_independent_augmented_ode_integrals():
    rate, duration, alpha, temperature, barrier = 1e8, 1e-7, .8, 600., .3
    distribution = emission(rate=rate, duration=duration, coefficient=alpha, barrier=barrier)
    beta = barrier * QE / (KB * temperature)
    # Integrate in normalized time to give inventory and both moment integrals
    # comparable scales; the moment implementation never uses this variable.
    def derivative(time, state):
        released_rate = rate * duration * state[0] * math.exp(-beta * (1 - alpha * state[0]))
        return [-released_rate, time * released_rate, time ** 2 * released_rate]
    solution = solve_ivp(derivative, [0, 1], [1., 0., 0.], rtol=2e-13, atol=1e-15)
    assert solution.success
    emitted = 1 - solution.y[0, -1]
    mean = duration * solution.y[1, -1] / emitted
    second = duration ** 2 * solution.y[2, -1] / emitted
    result = conditional_time_moments(distribution, quadrature_order=96)
    assert result.mean_offset_s == pytest.approx(mean, rel=3e-11)
    assert result.second_raw_moment_s2 == pytest.approx(second, rel=3e-11)
    assert result.variance_offset_s2 == pytest.approx(second - mean ** 2, rel=5e-11)
    initial_rate = rate * math.exp(-beta * (1 - alpha))
    incorrect = 1 / initial_rate - duration / math.expm1(initial_rate * duration)
    assert abs(result.mean_offset_s - incorrect) > .002 * duration


def test_stiff_interval_keeps_finite_actual_mean_when_all_normalized_times_underflow():
    rate = 1e200
    distribution = emission(rate=rate, duration=1e300)
    result = conditional_time_moments(distribution, quadrature_order=256)
    assert result.mean_offset_s > 0
    assert result.mean_offset_s * rate == pytest.approx(1, rel=2e-4)
    assert result.mean_offset_s / distribution.dt_s == 0
    assert not result.mean_underflowed
    assert result.second_raw_moment_s2 == result.variance_offset_s2 == 0
    assert result.second_raw_moment_underflowed and result.variance_underflowed
    assert result.log_second_raw_moment_s2 == pytest.approx(math.log(2) - 2 * math.log(rate), abs=8e-4)
    assert result.log_variance_offset_s2 == pytest.approx(-2 * math.log(rate), abs=2e-3)


def test_nonlinear_stiff_mean_uses_full_coverage_law_not_initial_rate():
    coupling, rate, temperature = 8., 1e150, 600.
    barrier = coupling * KB / QE * temperature
    distribution = emission(rate=rate, duration=1e250, coefficient=1., barrier=barrier)
    result = conditional_time_moments(distribution, quadrature_order=512)
    # For complete emission, swapping the time/coverage integrals gives
    # E[t]=(exp(c)-1)/(c*r_initial). Here initial barrier is exactly zero.
    expected = math.expm1(coupling) / (coupling * rate)
    assert result.mean_offset_s == pytest.approx(expected, rel=8e-5)
    assert result.mean_offset_s / distribution.dt_s == 0
    assert result.mean_offset_s > 100 / rate
    assert result.second_raw_moment_s2 > 0


def test_tiny_emitted_fraction_does_not_erase_conditional_moments():
    distribution = emission(rate=1e-200, duration=1e-200, coverage=1e300,
                            maximum=1e300, reference=1e300)
    assert distribution.released_m2 > 0
    assert math.exp(distribution.log_released_fraction) == 0
    result = conditional_time_moments(distribution, quadrature_order=8)
    assert result.mean_offset_s == pytest.approx(5e-201, rel=8e-14, abs=0)
    assert result.second_raw_moment_s2 == result.variance_offset_s2 == 0
    assert result.log_second_raw_moment_s2 == pytest.approx(2 * math.log(1e-200) - math.log(3), abs=3e-13)


def test_probability_nodes_and_weights_are_cached_open_positive_and_immutable():
    nodes, weights = gauss_legendre_probability_rule(4)
    cached_nodes, cached_weights = gauss_legendre_probability_rule(np.int64(4))
    assert nodes is cached_nodes and weights is cached_weights
    assert np.all((nodes > 0) & (nodes < 1)) and np.all(weights > 0)
    for power in range(8):
        assert math.fsum(weights * nodes ** power) == pytest.approx(1 / (power + 1), rel=4e-15)
    for array in (nodes, weights):
        with pytest.raises(ValueError):
            array.flags.writeable = True


def test_refinement_reports_changes_only_and_retains_underflowed_log_changes():
    refinement = refine_conditional_time_moments(emission(), quadrature_order=4)
    assert refinement.coarse.quadrature_order == 4
    assert refinement.refined.quadrature_order == 8
    assert "not certified" in refinement.interpretation
    for name, change in refinement.changes.items():
        direct = getattr(refinement.refined, name) - getattr(refinement.coarse, name)
        assert change.signed_change == pytest.approx(direct, rel=3e-9, abs=1e-16)
        assert change.absolute_change >= 0
    with pytest.raises(TypeError):
        refinement.changes["mean_offset_s"] = None
    stiff = refine_conditional_time_moments(emission(rate=1e200, duration=1e300), quadrature_order=8)
    change = stiff.changes["second_raw_moment_s2"]
    assert stiff.coarse.second_raw_moment_s2 == stiff.refined.second_raw_moment_s2 == 0
    assert change.change_underflowed
    assert change.log_absolute_change is not None and math.isfinite(change.log_absolute_change)
    assert change.relative_change > 0


def test_arbitrary_order_comparison_retains_underflowed_positive_moment_differences():
    distribution = emission(rate=1e200, duration=1e300)
    estimate = conditional_time_moments(distribution, quadrature_order=8)
    reference = conditional_time_moments(distribution, quadrature_order=64)
    comparisons = compare_conditional_time_moments(estimate, reference)
    change = comparisons["second_raw_moment_s2"]
    assert estimate.second_raw_moment_s2 == reference.second_raw_moment_s2 == 0
    assert estimate.second_raw_moment_underflowed and reference.second_raw_moment_underflowed
    assert change.absolute_change == 0 and change.change_underflowed
    assert change.log_absolute_change is not None and math.isfinite(change.log_absolute_change)
    assert math.copysign(1, change.signed_change) == -1
    expected_relative = abs(math.expm1(estimate.log_second_raw_moment_s2 - reference.log_second_raw_moment_s2))
    assert change.relative_change == pytest.approx(expected_relative, rel=3e-13)
    with pytest.raises(TypeError):
        comparisons["second_raw_moment_s2"] = None
    with pytest.raises(ValueError):
        compare_conditional_time_moments(estimate, None)
    with pytest.raises(ValueError):
        compare_conditional_time_moments(estimate, replace(reference, time_identity="different_law"))


@pytest.mark.parametrize("order", [0, 1, -1, True, 2.5, MAX_QUADRATURE_ORDER + 1])
def test_invalid_orders_reject(order):
    with pytest.raises(ValueError):
        conditional_time_moments(emission(), quadrature_order=order)


def test_wrong_emission_type_and_finite_overflow_fail_explicitly():
    with pytest.raises(ValueError, match="FrozenTemkinEmission"):
        conditional_time_moments(None, quadrature_order=8)
    with pytest.raises(TimeMomentRangeError, match="finite"):
        conditional_time_moments(emission(rate=1e-250, duration=1e200), quadrature_order=8)
