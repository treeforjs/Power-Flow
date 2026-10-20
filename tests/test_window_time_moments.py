"""New manufactured local-window moment checks; no saved oracle integration."""
from copy import copy
from dataclasses import FrozenInstanceError, asdict, replace
import json
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from mhdlab.emission_time_moments import MAX_QUADRATURE_ORDER, MomentChange
from mhdlab.temkin_emission_window import ReducedTemkinEmissionWindow, TemkinWindowRangeError
from mhdlab.temkin_release import LAW_IDENTITY
from mhdlab.window_time_moments import (
    INTERPRETATION, QUADRATURE_IDENTITY, TIME_IDENTITY, WindowTimeMoments,
    WindowTimeMomentRangeError, compare_window_time_moments,
    conditional_window_time_moments, refine_window_time_moments,
)


def window(*, coupling=0., rate=2.5, parent=2., start=.2, stop=1.1, tolerance=1e-7):
    return ReducedTemkinEmissionWindow(coupling, math.log(rate), parent, start, stop,
                                       log_mass_tolerance=tolerance)


def constant_moments(rate, duration):
    denominator = math.expm1(rate*duration)
    mean = 1/rate-duration/denominator
    second = 2/rate**2-(duration**2+2*duration/rate)/denominator
    variance = 1/rate**2-duration**2*math.exp(rate*duration)/denominator**2
    return mean, second, variance


def test_truncated_exponential_matches_three_independent_analytic_moments():
    law = window()
    result = conditional_window_time_moments(law, quadrature_order=64)
    for name, expected in zip(("mean_local_delay_s", "second_raw_moment_s2", "variance_local_delay_s2"),
                              constant_moments(2.5, .9)):
        assert getattr(result, name) == pytest.approx(expected, rel=5e-14)
    assert result.mean_local_delay_s < law.start_s+result.mean_local_delay_s
    assert result.numerical_log_tolerance == law.log_mass_tolerance
    assert result.law_identity == LAW_IDENTITY
    assert result.time_identity == TIME_IDENTITY
    assert result.quadrature_identity == QUADRATURE_IDENTITY
    assert result.interpretation == INTERPRETATION
    json.dumps(asdict(result), allow_nan=False)
    with pytest.raises(FrozenInstanceError):
        result.mean_local_delay_s = 0.


@pytest.mark.parametrize("order", [2, 4, 16])
def test_nearly_uniform_limit_resolves_exact_local_polynomial_moments(order):
    result = conditional_window_time_moments(window(rate=1e-100, parent=7., start=2., stop=5.),
                                             quadrature_order=order)
    assert result.mean_local_delay_s == pytest.approx(1.5, rel=3e-14)
    assert result.second_raw_moment_s2 == pytest.approx(3., rel=3e-14)
    assert result.variance_local_delay_s2 == pytest.approx(.75, rel=4e-14)


@pytest.mark.parametrize("coupling,rate,start,stop", [(3., 4., .3, 1.1), (12., 20., .05, .6)])
def test_nonlinear_window_matches_independent_original_population_ode(coupling, rate, start, stop):
    # The reference evolves the original surviving population, not the runtime
    # integral inversion or rebased window. Initial r is the supplied rate.
    def loss(population):
        return rate*population*math.exp(-coupling*(1-population))
    preceding = solve_ivp(lambda t, y: [-loss(y[0])], [0., start], [1.],
                          method="DOP853", rtol=2e-13, atol=2e-15)
    assert preceding.success
    initial = preceding.y[0, -1]
    def derivative(delay, state):
        release = loss(state[0])
        return [-release, delay*release, delay**2*release]
    reference = solve_ivp(derivative, [0., stop-start], [initial, 0., 0.],
                          method="DOP853", rtol=2e-13, atol=2e-15)
    assert reference.success
    emitted = initial-reference.y[0, -1]
    mean, second = reference.y[1:, -1]/emitted
    variance = second-mean**2
    result = conditional_window_time_moments(window(coupling=coupling, rate=rate, start=start, stop=stop),
                                             quadrature_order=64)
    assert result.mean_local_delay_s == pytest.approx(mean, rel=4e-11)
    assert result.second_raw_moment_s2 == pytest.approx(second, rel=5e-11)
    assert result.variance_local_delay_s2 == pytest.approx(variance, rel=6e-11)


def test_constant_late_window_has_no_full_cdf_subtraction_or_absolute_time_variance():
    early = conditional_window_time_moments(window(rate=1., start=0., stop=1., parent=2.), quadrature_order=64)
    late_law = window(rate=1., start=100_000., stop=100_001., parent=100_002., tolerance=.01)
    assert late_law.start_depletion > 99_000
    assert math.exp(late_law.log_source_fraction) == 0.
    late = conditional_window_time_moments(late_law, quadrature_order=64)
    for name in ("mean_local_delay_s", "second_raw_moment_s2", "variance_local_delay_s2"):
        assert getattr(late, name) == getattr(early, name)


def test_large_start_one_ulp_window_keeps_local_mean_and_centered_variance():
    start = 1e12
    stop = np.nextafter(start, math.inf)
    duration = stop-start
    law = window(rate=1e-100, parent=start+1., start=start, stop=stop)
    result = conditional_window_time_moments(law, quadrature_order=8)
    assert result.mean_local_delay_s == pytest.approx(duration/2, rel=6e-14)
    assert result.variance_local_delay_s2 == pytest.approx(duration**2/12, rel=8e-14)
    assert start+result.mean_local_delay_s in (start, stop)
    assert result.mean_local_delay_s > 0.
    assert result.variance_local_delay_s2 > 0.


def test_tiny_local_seconds_keep_logs_of_underflowed_square_moments():
    duration = 1e-200
    result = conditional_window_time_moments(window(rate=1e-200, parent=1., start=0., stop=duration),
                                             quadrature_order=8)
    assert result.mean_local_delay_s == pytest.approx(duration/2, rel=1e-13, abs=0.)
    assert not result.mean_underflowed
    assert result.second_raw_moment_s2 == result.variance_local_delay_s2 == 0.
    assert result.second_raw_moment_underflowed and result.variance_underflowed
    assert result.log_second_raw_moment_s2 == pytest.approx(2*math.log(duration)-math.log(3), abs=3e-13)
    assert result.log_variance_local_delay_s2 == pytest.approx(2*math.log(duration)-math.log(12), abs=4e-13)


def test_subfloat_mean_is_retained_in_log_seconds_without_linear_quantile_projection():
    duration = float(np.nextafter(0., 1.))
    law = window(rate=1e-300, parent=1e-320, start=0., stop=duration)
    with pytest.raises(TemkinWindowRangeError):
        law.quantile(.1)
    result = conditional_window_time_moments(law, quadrature_order=8)
    assert result.mean_local_delay_s == 0. and result.mean_underflowed
    assert result.second_raw_moment_underflowed and result.variance_underflowed
    assert result.log_mean_local_delay_s == pytest.approx(math.log(duration)-math.log(2), abs=3e-13)


def test_arbitrary_order_comparison_sign_reference_and_underflowed_log_difference():
    law = window(rate=1e200, start=0., stop=1e-190, parent=1e-190)
    estimate = conditional_window_time_moments(law, quadrature_order=8)
    reference = conditional_window_time_moments(law, quadrature_order=64)
    changes = compare_window_time_moments(estimate, reference)
    for name, log_name in (("mean_local_delay_s", "log_mean_local_delay_s"),
                           ("second_raw_moment_s2", "log_second_raw_moment_s2"),
                           ("variance_local_delay_s2", "log_variance_local_delay_s2")):
        change = changes[name]
        assert isinstance(change, MomentChange)
        assert change.relative_change == pytest.approx(abs(math.expm1(getattr(estimate, log_name)-getattr(reference, log_name))), rel=4e-12)
        assert math.copysign(1., change.signed_change) == math.copysign(1., getattr(estimate, log_name)-getattr(reference, log_name))
    assert estimate.second_raw_moment_s2 == reference.second_raw_moment_s2 == 0.
    assert changes["second_raw_moment_s2"].change_underflowed
    assert math.isfinite(changes["second_raw_moment_s2"].log_absolute_change)
    with pytest.raises(TypeError):
        changes["mean_local_delay_s"] = None
    equal = compare_window_time_moments(reference, reference)
    assert all(change == MomentChange(0., 0., 0., None, False) for change in equal.values())


def test_refinement_is_refined_minus_coarse_with_refined_denominator_and_immutable_changes():
    result = refine_window_time_moments(window(), quadrature_order=4)
    assert result.coarse.quadrature_order == 4 and result.refined.quadrature_order == 8
    assert "not certified" in result.interpretation
    for name, change in result.changes.items():
        expected = getattr(result.refined, name)-getattr(result.coarse, name)
        assert change.signed_change == pytest.approx(expected, abs=2e-16, rel=2e-9)
        assert change.relative_change == pytest.approx(abs(expected)/getattr(result.refined, name), abs=2e-15, rel=3e-8)
    with pytest.raises(TypeError):
        result.changes["variance_local_delay_s2"] = None
    with pytest.raises(FrozenInstanceError):
        result.coarse = result.refined


@pytest.mark.parametrize("field,value", [
    ("_local_fraction_log", -.1), ("local_duration_s", .7), ("_query_base_error", 0.),
    ("_diagnostics_json", "{}"), ("log_local_rate_s", 0.),
])
def test_forged_scalar_window_derived_fields_reject_before_quadrature(field, value):
    forged = copy(window())
    assert getattr(forged, field) != value
    object.__setattr__(forged, field, value)
    with pytest.raises(ValueError, match="constructor replay"):
        conditional_window_time_moments(forged, quadrature_order=4)


@pytest.mark.parametrize("order", [True, np.bool_(False), 1, 0, -2, 2.5, "4", None,
                                    np.ma.array(4, mask=False), MAX_QUADRATURE_ORDER+1])
def test_invalid_orders_reject(order):
    with pytest.raises(ValueError, match="quadrature_order"):
        conditional_window_time_moments(window(), quadrature_order=order)


def test_refinement_order_cap_rejects_before_expensive_rule(monkeypatch):
    import mhdlab.window_time_moments as module
    monkeypatch.setattr(module, "conditional_window_time_moments", lambda *a, **k: pytest.fail("cap must precede rule"))
    with pytest.raises(ValueError, match="quadrature_order"):
        module.refine_window_time_moments(window(), quadrature_order=MAX_QUADRATURE_ORDER)


@pytest.mark.parametrize("value", [None, {}, np.ma.array(1), True])
def test_wrong_window_record_rejects(value):
    with pytest.raises(ValueError, match="ReducedTemkinEmissionWindow"):
        conditional_window_time_moments(value, quadrature_order=4)


def test_scalar_query_qualification_failures_propagate(monkeypatch):
    def reject(*a, **k):
        raise TemkinWindowRangeError("manufactured unresolved query")
    monkeypatch.setattr(ReducedTemkinEmissionWindow, "log_quantile", reject)
    with pytest.raises(TemkinWindowRangeError, match="unresolved query"):
        conditional_window_time_moments(window(), quadrature_order=4)


@pytest.mark.parametrize("logs", [[0., -1., 0., 1.], [0., 0., 0., 0.], [0., math.nan, 1., 2.], [0.]])
def test_unresolved_or_nonmonotone_query_arrays_reject(monkeypatch, logs):
    monkeypatch.setattr(ReducedTemkinEmissionWindow, "log_quantile", lambda *a: np.array(logs))
    with pytest.raises(WindowTimeMomentRangeError):
        conditional_window_time_moments(window(), quadrature_order=4)


@pytest.mark.parametrize("field,value", [("mean_local_delay_s", True), ("log_mean_local_delay_s", math.nan),
                                        ("log_second_raw_moment_s2", np.ma.array(-1.)),
                                        ("variance_underflowed", 1), ("quadrature_order", True),
                                        ("numerical_log_tolerance", 0.), ("time_identity", "global_time"),
                                        ("law_identity", "other"), ("quadrature_identity", "certified")])
def test_direct_moment_record_invalid_fields_reject(field, value):
    result = conditional_window_time_moments(window(), quadrature_order=4)
    with pytest.raises(ValueError):
        replace(result, **{field: value})


def test_direct_record_enforces_centered_second_moment_identity_and_compare_revalidates():
    with pytest.raises(ValueError, match="mean squared"):
        WindowTimeMoments(1., 1., 1., 0., 0., 0., False, False, False, 2, 1e-7)
    result = conditional_window_time_moments(window(), quadrature_order=4)
    forged = copy(result)
    object.__setattr__(forged, "mean_underflowed", True)
    with pytest.raises(ValueError):
        compare_window_time_moments(forged, result)
    with pytest.raises(ValueError):
        compare_window_time_moments(None, result)


def test_unrepresentable_second_moment_rejects_without_infinite_output():
    with pytest.raises(WindowTimeMomentRangeError, match="finite"):
        conditional_window_time_moments(window(rate=1e-250, parent=1e200, start=0., stop=1e200),
                                         quadrature_order=4)
