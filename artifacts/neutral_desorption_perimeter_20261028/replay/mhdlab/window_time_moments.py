"""Deterministic local-delay moments conditional on a Temkin source window.

Open probability-space Gauss-Legendre quadrature integrates the qualified
window's log quantiles in actual local seconds. Neither the source-window start
nor an absolute event-time origin is added. Centered scaled differences compute
variance without subtracting two raw moments; logarithms retain positive
estimates whose SI projections underflow.

Finite-order quadrature and n-versus-2n changes are numerical estimates, not
certified error bounds or physical uncertainties. The scalar window's separate
estimated query qualification is retained and can reject a requested rule.
"""
from __future__ import annotations

from dataclasses import dataclass, fields
import math
from numbers import Integral, Real
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .emission_time_moments import (
    MAX_QUADRATURE_ORDER, MomentChange, REFINEMENT_IDENTITY,
    TimeMomentRangeError, gauss_legendre_probability_rule,
)
from .temkin_emission_window import ReducedTemkinEmissionWindow
from .temkin_release import LAW_IDENTITY


TIME_IDENTITY = "temkin_emitted_population_conditional_window_local_delay_v1"
QUADRATURE_IDENTITY = "gauss_legendre_open_probability_window_log_quantiles_actual_local_seconds_v1"
INTERPRETATION = "finite-order conditional source-window local-delay quadrature estimate; not a certified error bound"
_REFINEMENT_INTERPRETATION = "n-versus-2n quadrature changes only; not certified error bounds or physical uncertainties"
_LOG_MAX = math.log(float(np.finfo(float).max))
_EPS = float(np.finfo(float).eps)
_MOMENTS = (("mean_local_delay_s", "log_mean_local_delay_s", "mean_underflowed"),
            ("second_raw_moment_s2", "log_second_raw_moment_s2", "second_raw_moment_underflowed"),
            ("variance_local_delay_s2", "log_variance_local_delay_s2", "variance_underflowed"))


class WindowTimeMomentRangeError(TimeMomentRangeError):
    """A local-delay moment lacks a finite supported numerical estimate."""


def _real(value, name):
    if (np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_))
            or not isinstance(value, Real)):
        raise ValueError(f"{name} must be a finite real scalar, not boolean or masked")
    try:
        result = float(value)
    except (ValueError, OverflowError, TypeError) as error:
        raise ValueError(f"{name} must be finite and representable") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _order(value):
    if (np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_))
            or not isinstance(value, Integral) or not 2 <= value <= MAX_QUADRATURE_ORDER):
        raise ValueError(f"quadrature_order must be an integer in [2,{MAX_QUADRATURE_ORDER}]")
    return int(value)


def _represented(log_value):
    if not math.isfinite(log_value) or log_value > _LOG_MAX:
        raise WindowTimeMomentRangeError("window time moment exceeds the supported finite float64 range")
    try:
        value = math.exp(log_value)
    except OverflowError as error:
        raise WindowTimeMomentRangeError("window time moment exceeds the finite float64 range") from error
    if not math.isfinite(value):
        raise WindowTimeMomentRangeError("window time moment exceeds the finite float64 range")
    return value, value == 0.


def _validated_window(window):
    if not isinstance(window, ReducedTemkinEmissionWindow):
        raise ValueError("window must be a qualified ReducedTemkinEmissionWindow")
    rebuilt = ReducedTemkinEmissionWindow(**{
        item.name: getattr(window, item.name)
        for item in fields(ReducedTemkinEmissionWindow) if item.init
    })
    # Frozen records are still forgeable through object.__setattr__. Check the
    # constructor-derived caches too, then compute only with the rebuilt law.
    for item in fields(ReducedTemkinEmissionWindow):
        supplied, expected = getattr(window, item.name), getattr(rebuilt, item.name)
        if type(supplied) is not type(expected) or supplied != expected:
            raise ValueError(f"window {item.name} disagrees with its qualified constructor replay")
    if window.manifest() != rebuilt.manifest():
        raise ValueError("window manifest disagrees with its qualified constructor replay")
    return rebuilt


@dataclass(frozen=True, slots=True)
class WindowTimeMoments:
    mean_local_delay_s: float
    second_raw_moment_s2: float
    variance_local_delay_s2: float
    log_mean_local_delay_s: float
    log_second_raw_moment_s2: float
    log_variance_local_delay_s2: float
    mean_underflowed: bool
    second_raw_moment_underflowed: bool
    variance_underflowed: bool
    quadrature_order: int
    numerical_log_tolerance: float
    law_identity: str = LAW_IDENTITY
    time_identity: str = TIME_IDENTITY
    quadrature_identity: str = QUADRATURE_IDENTITY
    interpretation: str = INTERPRETATION

    def __post_init__(self):
        object.__setattr__(self, "quadrature_order", _order(self.quadrature_order))
        tolerance = _real(self.numerical_log_tolerance, "numerical_log_tolerance")
        if tolerance <= 0:
            raise ValueError("numerical_log_tolerance must be positive")
        object.__setattr__(self, "numerical_log_tolerance", tolerance)
        for name, expected in (("law_identity", LAW_IDENTITY), ("time_identity", TIME_IDENTITY),
                               ("quadrature_identity", QUADRATURE_IDENTITY), ("interpretation", INTERPRETATION)):
            if getattr(self, name) != expected:
                raise ValueError(f"unsupported window moment {name}")
        for value_name, log_name, flag_name in _MOMENTS:
            log_value = _real(getattr(self, log_name), log_name)
            value = _real(getattr(self, value_name), value_name)
            expected, underflowed = _represented(log_value)
            if type(getattr(self, flag_name)) is not bool or getattr(self, flag_name) != underflowed:
                raise ValueError(f"{flag_name} must match its positive log-moment projection")
            if value != expected:
                raise ValueError(f"{value_name} must match its positive log-moment projection")
            object.__setattr__(self, log_name, log_value)
            object.__setattr__(self, value_name, value)
        allowance = 64*_EPS*math.fsum((1., abs(self.log_second_raw_moment_s2),
                                        abs(2*self.log_mean_local_delay_s), abs(self.log_variance_local_delay_s2)))
        if (self.log_second_raw_moment_s2 < 2*self.log_mean_local_delay_s-allowance
                or self.log_variance_local_delay_s2 > self.log_second_raw_moment_s2+allowance):
            raise ValueError("window moment logs violate positive population moment inequalities")
        larger = max(2*self.log_mean_local_delay_s, self.log_variance_local_delay_s2)
        smaller = min(2*self.log_mean_local_delay_s, self.log_variance_local_delay_s2)
        reconstructed_second = math.fsum((larger, math.log1p(math.exp(smaller-larger))))
        if abs(reconstructed_second-self.log_second_raw_moment_s2) > allowance:
            raise ValueError("window moment logs violate second moment = mean squared + variance")


def conditional_window_time_moments(window: ReducedTemkinEmissionWindow, *, quadrature_order: int) -> WindowTimeMoments:
    """Estimate E[d], E[d²], Var[d] for local delay from the window start.

    Log quantiles, not projected linear delays, feed quadrature. A positive
    subnormal SI value is a rounded projection, while each underflow flag marks
    a positive logarithmic moment projected to zero. The window's own estimated
    log-time query qualification applies at every open node and is separate
    from this finite-order quadrature resolution.
    """
    order = _order(quadrature_order)
    law = _validated_window(window)
    nodes, weights = gauss_legendre_probability_rule(order)
    logs = np.asarray(law.log_quantile(nodes), dtype=np.float64)
    if (logs.shape != (order,) or not np.all(np.isfinite(logs))
            or np.any(logs[1:] < logs[:-1])):
        raise WindowTimeMomentRangeError("window log quantiles must be finite and monotonically increasing")
    anchor = float(logs[-1])
    with np.errstate(under="ignore"):
        scaled = np.exp(logs-anchor)
    if np.any(scaled <= 0):
        raise WindowTimeMomentRangeError("relative local-delay quantiles are numerically unresolved")
    normalization = math.fsum(float(weight) for weight in weights)
    mean_scaled = math.fsum(float(w)*float(value) for w, value in zip(weights, scaled))/normalization
    second_scaled = math.fsum(float(w)*float(value)**2 for w, value in zip(weights, scaled))/normalization
    variance_scaled = math.fsum(float(w)*(float(value)-mean_scaled)**2
                               for w, value in zip(weights, scaled))/normalization
    if min(mean_scaled, second_scaled, variance_scaled) <= 0 or not all(
            math.isfinite(value) for value in (mean_scaled, second_scaled, variance_scaled)):
        raise WindowTimeMomentRangeError("conditional local-delay variance or scaled moment is unresolved")
    log_mean = math.fsum((anchor, math.log(mean_scaled)))
    log_second = math.fsum((2*anchor, math.log(second_scaled)))
    log_variance = math.fsum((2*anchor, math.log(variance_scaled)))
    mean, mean_underflow = _represented(log_mean)
    second, second_underflow = _represented(log_second)
    variance, variance_underflow = _represented(log_variance)
    return WindowTimeMoments(mean, second, variance, log_mean, log_second, log_variance,
                             mean_underflow, second_underflow, variance_underflow,
                             order, law.log_mass_tolerance)


def _validated_moments(record):
    if not isinstance(record, WindowTimeMoments):
        raise ValueError("estimate and reference must be WindowTimeMoments records")
    return WindowTimeMoments(**{item.name: getattr(record, item.name) for item in fields(WindowTimeMoments)})


def _change(estimate_log, reference_log, *, relative_reference_log=None):
    if estimate_log == reference_log:
        return MomentChange(0., 0., 0., None, False)
    upper, lower = max(estimate_log, reference_log), min(estimate_log, reference_log)
    log_difference = math.fsum((upper, math.log(-math.expm1(lower-upper))))
    absolute, underflowed = _represented(log_difference)
    denominator = reference_log if relative_reference_log is None else relative_reference_log
    relative, _ = _represented(log_difference-denominator)
    return MomentChange(math.copysign(absolute, estimate_log-reference_log), absolute,
                        relative, log_difference, underflowed)


def compare_window_time_moments(estimate: WindowTimeMoments, reference: WindowTimeMoments) -> Mapping[str, MomentChange]:
    """Return estimate-minus-reference changes, relative to the reference.

    Both records are finite quadrature estimates. Comparisons use their logs
    even if both SI projections have rounded to zero; callers must associate
    the records with the intended common source before interpreting changes.
    """
    estimate, reference = _validated_moments(estimate), _validated_moments(reference)
    return MappingProxyType({name: _change(getattr(estimate, log_name), getattr(reference, log_name))
                             for name, log_name, _ in _MOMENTS})


@dataclass(frozen=True, slots=True)
class WindowTimeMomentRefinement:
    coarse: WindowTimeMoments
    refined: WindowTimeMoments
    changes: Mapping[str, MomentChange]
    refinement_identity: str = REFINEMENT_IDENTITY
    interpretation: str = _REFINEMENT_INTERPRETATION

    def __post_init__(self):
        coarse, refined = _validated_moments(self.coarse), _validated_moments(self.refined)
        if refined.quadrature_order != 2*coarse.quadrature_order:
            raise ValueError("refinement requires n and 2n quadrature orders")
        if self.refinement_identity != REFINEMENT_IDENTITY or self.interpretation != _REFINEMENT_INTERPRETATION:
            raise ValueError("unsupported window moment refinement identity or interpretation")
        # Refinement relative changes use the refined estimate as denominator.
        expected = {name: _change(getattr(refined, log_name), getattr(coarse, log_name),
                                  relative_reference_log=getattr(refined, log_name))
                    for name, log_name, _ in _MOMENTS}
        if not isinstance(self.changes, Mapping) or dict(self.changes) != expected:
            raise ValueError("refinement changes must match the n and 2n estimates")
        object.__setattr__(self, "coarse", coarse)
        object.__setattr__(self, "refined", refined)
        object.__setattr__(self, "changes", MappingProxyType(expected))


def refine_window_time_moments(window: ReducedTemkinEmissionWindow, *, quadrature_order: int) -> WindowTimeMomentRefinement:
    """Return n/2n estimates and refined-minus-coarse changes, not error bounds."""
    order = _order(quadrature_order)
    _order(2*order)
    coarse = conditional_window_time_moments(window, quadrature_order=order)
    refined = conditional_window_time_moments(window, quadrature_order=2*order)
    changes = {}
    for name, log_name, _ in _MOMENTS:
        changes[name] = _change(getattr(refined, log_name), getattr(coarse, log_name),
                                relative_reference_log=getattr(refined, log_name))
    return WindowTimeMomentRefinement(coarse, refined, MappingProxyType(changes))
