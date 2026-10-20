"""Deterministic conditional Temkin event-time moments in actual seconds.

Open probability-space Gauss-Legendre nodes integrate the law's log quantiles.
Anchoring actual elapsed times avoids dividing early events by an enormous
interval. Centered scaled differences compute variance without subtracting raw
moments. Log moments survive float64 underflow; rounded zero is flagged.

These are finite-order quadrature estimates. Increasing the order measures a
change in estimates, not a certified error bound or physical uncertainty.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
from numbers import Integral
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.special import roots_legendre

from .temkin_emission_time import FrozenTemkinEmission, TIME_IDENTITY
from .temkin_release import LAW_IDENTITY


QUADRATURE_IDENTITY = "gauss_legendre_open_probability_log_quantiles_actual_seconds_v1"
REFINEMENT_IDENTITY = "successive_n_and_2n_quadrature_changes_not_error_bounds_v1"
MAX_QUADRATURE_ORDER = 4096
_LOG_MAX = math.log(float(np.finfo(float).max))


class TimeMomentRangeError(ValueError):
    """A requested moment or quadrature rule lacks a finite supported result."""


def _order(value):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or not 2 <= value <= MAX_QUADRATURE_ORDER):
        raise ValueError(f"quadrature_order must be an integer in [2,{MAX_QUADRATURE_ORDER}]")
    return int(value)


def _immutable(value):
    value = np.asarray(value, dtype=np.float64)
    return np.frombuffer(np.ascontiguousarray(value).tobytes(), dtype=np.float64).reshape(value.shape)


@lru_cache(maxsize=16)
def _cached_rule(order):
    roots, weights = roots_legendre(order)
    nodes, weights = .5 * (roots + 1), .5 * weights
    weights /= math.fsum(float(value) for value in weights)
    if (not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(weights))
            or np.any(nodes <= 0) or np.any(nodes >= 1) or np.any(weights <= 0)):
        raise TimeMomentRangeError("Gauss-Legendre rule does not resolve positive open probability nodes")
    return _immutable(nodes), _immutable(weights)


def gauss_legendre_probability_rule(quadrature_order: int):
    """Return cached immutable (open probability nodes, positive weights).

    The order cap is a computational bound, not an accuracy qualification.
    Callers cannot make the returned byte-backed arrays writable.
    """
    return _cached_rule(_order(quadrature_order))


def _represented(log_value):
    if not math.isfinite(log_value) or log_value > _LOG_MAX:
        raise TimeMomentRangeError("conditional time moment exceeds supported finite float64 range")
    try:
        value = math.exp(log_value)
    except OverflowError as error:
        raise TimeMomentRangeError("conditional time moment exceeds finite float64 range") from error
    if not math.isfinite(value):
        raise TimeMomentRangeError("conditional time moment exceeds finite float64 range")
    return value, value == 0


@dataclass(frozen=True, slots=True)
class TimeMoments:
    mean_offset_s: float
    second_raw_moment_s2: float
    variance_offset_s2: float
    log_mean_offset_s: float
    log_second_raw_moment_s2: float
    log_variance_offset_s2: float
    mean_underflowed: bool
    second_raw_moment_underflowed: bool
    variance_underflowed: bool
    quadrature_order: int
    numerical_log_tolerance: float
    law_identity: str = LAW_IDENTITY
    time_identity: str = TIME_IDENTITY
    quadrature_identity: str = QUADRATURE_IDENTITY
    interpretation: str = "finite-order conditional-time quadrature estimate; not a certified error bound"


def conditional_time_moments(emission: FrozenTemkinEmission, *, quadrature_order: int) -> TimeMoments:
    """Estimate E[t], E[t²], and Var[t] of the emitted-particle conditional law.

    Logs are of elapsed seconds or squared seconds relative to the stated SI
    unit. Flags indicate a positive log moment rounded to zero. Positive
    subnormals remain rounded approximations; log fields retain the estimate
    before that representation step. No initial-rate approximation is used for
    nonlinear emission. The timing kernel's own numerical qualification is
    retained separately from quadrature resolution.
    """
    if not isinstance(emission, FrozenTemkinEmission):
        raise ValueError("emission must be a qualified FrozenTemkinEmission")
    if emission.law_identity != LAW_IDENTITY or emission.time_identity != TIME_IDENTITY:
        raise ValueError("unsupported conditional emission law or time identity")
    order = _order(quadrature_order)
    nodes, weights = gauss_legendre_probability_rule(order)
    logs = np.asarray(emission.log_quantile(nodes), dtype=np.float64)
    if logs.shape != (order,) or not np.all(np.isfinite(logs)) or np.any(logs[1:] < logs[:-1]):
        raise TimeMomentRangeError("conditional log quantiles must be finite and monotonically increasing")
    anchor = float(logs[-1])
    with np.errstate(under="ignore"):
        scaled = np.exp(logs - anchor)
    total_weight = math.fsum(float(value) for value in weights)
    mean_scaled = math.fsum(float(weight) * float(value) for weight, value in zip(weights, scaled)) / total_weight
    second_scaled = math.fsum(float(weight) * float(value) ** 2 for weight, value in zip(weights, scaled)) / total_weight
    variance_scaled = math.fsum(float(weight) * (float(value) - mean_scaled) ** 2
                               for weight, value in zip(weights, scaled)) / total_weight
    if min(mean_scaled, second_scaled, variance_scaled) <= 0 or not all(
            math.isfinite(value) for value in (mean_scaled, second_scaled, variance_scaled)):
        raise TimeMomentRangeError("conditional variance or scaled moment is numerically unresolved")
    log_mean = math.fsum((anchor, math.log(mean_scaled)))
    log_second = math.fsum((2 * anchor, math.log(second_scaled)))
    log_variance = math.fsum((2 * anchor, math.log(variance_scaled)))
    mean, mean_underflow = _represented(log_mean)
    second, second_underflow = _represented(log_second)
    variance, variance_underflow = _represented(log_variance)
    return TimeMoments(mean, second, variance, log_mean, log_second, log_variance,
                       mean_underflow, second_underflow, variance_underflow, order,
                       float(emission.numerical_log_tolerance))


@dataclass(frozen=True, slots=True)
class MomentChange:
    signed_change: float
    absolute_change: float
    relative_change: float
    log_absolute_change: float | None
    change_underflowed: bool


@dataclass(frozen=True, slots=True)
class TimeMomentRefinement:
    coarse: TimeMoments
    refined: TimeMoments
    changes: Mapping[str, MomentChange]
    refinement_identity: str = REFINEMENT_IDENTITY
    interpretation: str = "n-versus-2n quadrature changes only; not certified error bounds or physical uncertainties"


def _change(coarse_log, fine_log, *, relative_reference_log=None):
    if coarse_log == fine_log:
        return MomentChange(0., 0., 0., None, False)
    upper, lower = max(coarse_log, fine_log), min(coarse_log, fine_log)
    log_difference = math.fsum((upper, math.log(-math.expm1(lower - upper))))
    absolute, underflowed = _represented(log_difference)
    denominator = fine_log if relative_reference_log is None else relative_reference_log
    relative, _ = _represented(log_difference - denominator)
    return MomentChange(math.copysign(absolute, fine_log - coarse_log), absolute,
                        relative, log_difference, underflowed)


def compare_conditional_time_moments(estimate: TimeMoments, reference: TimeMoments) -> Mapping[str, MomentChange]:
    """Compare arbitrary orders using log moments, with sign estimate-reference.

    Relative changes are normalized to the supplied reference estimate. That
    estimate is finite quadrature, not an exact solution. Positive moments are
    compared even when both corresponding SI values have rounded to zero; the
    input records retain their original log values and underflow flags.
    """
    if not isinstance(estimate, TimeMoments) or not isinstance(reference, TimeMoments):
        raise ValueError("estimate and reference must be TimeMoments records")
    for name, identity in (("law_identity", LAW_IDENTITY), ("time_identity", TIME_IDENTITY),
                           ("quadrature_identity", QUADRATURE_IDENTITY)):
        if getattr(estimate, name) != identity or getattr(reference, name) != identity:
            raise ValueError(f"time-moment comparisons require the same supported {name}")
    comparisons = {}
    for name, log_name in (("mean_offset_s", "log_mean_offset_s"),
                           ("second_raw_moment_s2", "log_second_raw_moment_s2"),
                           ("variance_offset_s2", "log_variance_offset_s2")):
        estimate_log, reference_log = getattr(estimate, log_name), getattr(reference, log_name)
        if not math.isfinite(estimate_log) or not math.isfinite(reference_log):
            raise ValueError("time-moment logarithms must retain finite positive-moment estimates")
        comparisons[name] = _change(reference_log, estimate_log, relative_reference_log=reference_log)
    return MappingProxyType(comparisons)


def refine_conditional_time_moments(emission: FrozenTemkinEmission, *, quadrature_order: int) -> TimeMomentRefinement:
    """Return n and 2n estimates and their signed/absolute/relative changes.

    Relative change uses the refined log moment, so it remains meaningful when
    both SI moment values have underflowed to zero. An exactly equal pair of
    computed log moments has zero change and no log-absolute-change value.
    """
    order = _order(quadrature_order)
    _order(2 * order)
    coarse = conditional_time_moments(emission, quadrature_order=order)
    refined = conditional_time_moments(emission, quadrature_order=2 * order)
    changes = {name: _change(getattr(coarse, log_name), getattr(refined, log_name))
               for name, log_name in (("mean_offset_s", "log_mean_offset_s"),
                                      ("second_raw_moment_s2", "log_second_raw_moment_s2"),
                                      ("variance_offset_s2", "log_variance_offset_s2"))}
    return TimeMomentRefinement(coarse, refined, MappingProxyType(changes))
