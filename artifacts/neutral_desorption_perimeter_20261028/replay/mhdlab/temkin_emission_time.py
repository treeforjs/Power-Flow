"""Conditional emission times for the distinct frozen-temperature Temkin law.

This scalar-component interface draws no random numbers and advances no gap
particles.  It conditions on particles emitted during one frozen-temperature
interval.  It does not use a truncated exponential with the initial rate for
nonlinear coverage-dependent release.

For s=log(Ninitial/N), c=Ed/(kBT)*alpha*Ninitial/Nref and initial rate r0,
elapsed time is I(s)/r0, where I(s)=integral_0^s exp[c(1-exp(-u))]du.
The cumulative emitted fraction is 1-exp(-s).  The conditional CDF divides this
fraction by its value at dt.  Quantiles target that emitted population, then
evaluate log(I(s))-log(r0); they never normalize by dt before inversion.  Thus
a fully depleted very long interval can still have finite, very early events.

The validated core trial is cached for number/energy accounting.  Conditional
fractions are evaluated in log form from the same separated law, rather than
by dividing rounded extensive coverages; a representable release can have an
unrepresentably small released *fraction*.  The core's c<=512 qualification
and constant-attempt-frequency convention apply unchanged.

The numerical log tolerance is exported below.  Endpoint projections smaller
than this tolerance accommodate equivalent qualified scalar roots; larger
inconsistencies raise an error.  A CDF probability below float64 range rounds
to zero.  A positive quantile time below float64 range raises explicitly;
``log_quantile`` retains its logarithm.  Coarsely quantized positive subnormal
offsets also raise if their log-time error exceeds the qualification.
Returned offsets are local elapsed
times, not absolute timestamps.  No precision at an absolute time origin is
implied.  A certified tiny-hazard limit is conditionally uniform; it is not a
general replacement of the coverage-dependent law by a constant-rate model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math

import numpy as np
from scipy import optimize

from .temkin_release import (
    INVERSION_LOG_RESIDUAL_TOLERANCE, LAW_IDENTITY, TemkinInversionError,
    TemkinReleaseTrial, TemkinReservoir, _log_integral, _thermal_parameters,
    release_step,
)
from .thermal_release import SurfaceInventory, _immutable, _log_hazard, _ratio, _real


TIME_IDENTITY = "temkin_conditional_emitted_population_frozen_temperature_v1"
NUMERICAL_LOG_TOLERANCE = 8 * INVERSION_LOG_RESIDUAL_TOLERANCE
_EPS = float(np.finfo(float).eps)
_LOG_EPS = math.log(_EPS)
_LOG_NEGLIGIBLE_SURVIVOR_DEPLETION = math.log(80.)


class NoTemkinEmissionError(ValueError):
    """The core trial has no positive representable emitted population."""


class TemkinEmissionTimeRangeError(ValueError):
    """Requested timing cannot be represented within the declared qualification."""


def _values(value, name: str) -> np.ndarray:
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked values")
    try:
        raw = np.asarray(value, dtype=object)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a real scalar or one-dimensional array") from exc
    if raw.ndim > 1 or raw.size == 0:
        raise ValueError(f"{name} must be a real scalar or nonempty one-dimensional array")
    for item in raw.flat:
        _real(item, name)
    return np.asarray(raw, dtype=float)


def _output(value: np.ndarray):
    return float(value) if value.ndim == 0 else _immutable(value)


def _law_log_depletion(log_hazard: float, coupling: float) -> float | None:
    """Return log(s), or None when survival is negligible for float64 u<1.

    s>=80 gives survival <1.81e-35.  The closest float64 probability below one
    has 1-u>=1.11e-16, so neglecting that survivor changes a target log-depletion
    by less than 2e-19.  This is a timing representation bound, not a rate floor.
    """
    if coupling == 0:
        return None if log_hazard >= _LOG_NEGLIGIBLE_SURVIVOR_DEPLETION else log_hazard
    if log_hazard + math.log(max(1., coupling)) < _LOG_EPS - math.log(4.):
        return log_hazard
    if log_hazard >= _log_integral(_LOG_NEGLIGIBLE_SURVIVOR_DEPLETION, coupling):
        return None
    lower = log_hazard - coupling
    upper = min(log_hazard, _LOG_NEGLIGIBLE_SURVIVOR_DEPLETION)

    def residual(value):
        return _log_integral(value, coupling) - log_hazard

    low_value, high_value = residual(lower), residual(upper)
    margin = 32 * _EPS * max(1., abs(lower), abs(upper))
    if low_value > 0:
        lower -= margin
        low_value = residual(lower)
    if high_value < 0:
        upper += margin
        high_value = residual(upper)
    if low_value > 0 or high_value < 0:
        raise TemkinInversionError("conditional Temkin time failed to bracket the monotone integral")
    if abs(high_value) <= 4 * _EPS:
        solution = upper
    elif abs(low_value) <= 4 * _EPS:
        solution = lower
    else:
        try:
            solution, result = optimize.brentq(residual, lower, upper, xtol=1e-14,
                                               rtol=4 * _EPS, maxiter=128,
                                               full_output=True, disp=False)
        except (ValueError, RuntimeError) as exc:
            raise TemkinInversionError("conditional Temkin time inversion failed") from exc
        if not result.converged:
            raise TemkinInversionError("conditional Temkin time inversion did not converge")
    if abs(residual(solution)) > INVERSION_LOG_RESIDUAL_TOLERANCE:
        raise TemkinInversionError("conditional Temkin time inversion exceeded its log-time tolerance")
    return solution


def _log_fraction(log_depletion: float | None) -> float:
    if log_depletion is None:
        return 0.
    if log_depletion < _LOG_EPS - math.log(4.):
        return log_depletion
    return math.log(-math.expm1(-math.exp(log_depletion)))


@dataclass(frozen=True, slots=True)
class FrozenTemkinEmission:
    """One explicit reservoir/bin/temperature/interval conditional-time law.

    Constructor arguments are all required. ``trial`` is the cached immutable
    scalar core trial. ``released_m2`` remains its extensive float64 budget,
    while ``log_released_fraction`` describes the corresponding mathematical
    conditional normalization without an extensive-quantity division.
    """

    reservoir: TemkinReservoir
    initial_coverage_m2: float
    temperature_k: float
    dt_s: float
    trial: TemkinReleaseTrial = field(init=False, repr=False)
    _coupling: float = field(init=False, repr=False)
    _log_rate: float = field(init=False, repr=False)
    _log_full_depletion: float | None = field(init=False, repr=False)
    _log_full_fraction: float = field(init=False, repr=False)
    _uniform_limit: bool = field(init=False, repr=False)

    def __post_init__(self):
        if not isinstance(self.reservoir, TemkinReservoir):
            raise ValueError("reservoir must be a TemkinReservoir")
        for name in ("initial_coverage_m2", "temperature_k", "dt_s"):
            object.__setattr__(self, name, _real(getattr(self, name), name,
                                               positive=name != "initial_coverage_m2"))
        trial = release_step(SurfaceInventory([[self.initial_coverage_m2]]), [self.reservoir],
                             [self.temperature_k], self.dt_s)
        object.__setattr__(self, "trial", trial)
        if trial.released_m2[0, 0] == 0:
            raise NoTemkinEmissionError("interval has no positive representable Temkin emission")
        if self.reservoir.coverage_coefficient == 0 or self.reservoir.activation_energy_ev == 0:
            coupling = 0.
            log_rate = _log_hazard(self.reservoir._energy_record(), self.temperature_k, 1.)
        else:
            coupling, initial_barrier = _thermal_parameters(self.reservoir, self.initial_coverage_m2,
                                                           self.temperature_k)
            log_rate = math.fsum((math.log(self.reservoir.attempt_frequency_s), -initial_barrier))
        if not math.isfinite(log_rate):
            raise TemkinEmissionTimeRangeError("positive emitted population has no finite log initial rate")
        log_hazard = math.fsum((log_rate, math.log(self.dt_s)))
        uniform = log_hazard + math.log1p(coupling) < _LOG_EPS - math.log(4.)
        full_depletion = _law_log_depletion(log_hazard, coupling)
        object.__setattr__(self, "_coupling", coupling)
        object.__setattr__(self, "_log_rate", log_rate)
        object.__setattr__(self, "_log_full_depletion", full_depletion)
        object.__setattr__(self, "_log_full_fraction", _log_fraction(full_depletion))
        object.__setattr__(self, "_uniform_limit", uniform)

    @property
    def law_identity(self):
        return LAW_IDENTITY

    @property
    def time_identity(self):
        return TIME_IDENTITY

    @property
    def numerical_log_tolerance(self):
        return NUMERICAL_LOG_TOLERANCE

    @property
    def released_m2(self):
        return float(self.trial.released_m2[0, 0])

    @property
    def remaining_coverage_m2(self):
        return float(self.trial.trial_inventory.coverage_m2[0, 0])

    @property
    def log_released_fraction(self):
        return self._log_full_fraction

    @property
    def uniform_limit(self):
        return self._uniform_limit

    def cdf(self, time_offset_s):
        """Conditional emitted-population CDF, for offsets in [0,dt]."""
        values = _values(time_offset_s, "time_offset_s")
        if np.any(values > self.dt_s):
            raise ValueError("time_offset_s must lie in [0,dt_s]")
        result = np.empty_like(values)
        for location in np.ndindex(values.shape):
            time = float(values[location])
            if time == 0:
                probability = 0.
            elif time == self.dt_s:
                probability = 1.
            elif self._uniform_limit:
                probability = _ratio(time, self.dt_s)
            else:
                log_hazard = math.fsum((self._log_rate, math.log(time)))
                log_fraction = _log_fraction(_law_log_depletion(log_hazard, self._coupling))
                log_probability = math.fsum((log_fraction, -self._log_full_fraction))
                if log_probability > NUMERICAL_LOG_TOLERANCE:
                    raise TemkinEmissionTimeRangeError("conditional CDF exceeded its qualified interval")
                probability = math.exp(min(0., log_probability))
            result[location] = probability
        return _output(result)

    def _target_log_depletion(self, probability):
        log_u = math.log(probability)
        log_target_fraction = math.fsum((log_u, self._log_full_fraction))
        if log_target_fraction < math.log(1e-6):
            if log_target_fraction < _LOG_EPS - math.log(4.):
                return log_target_fraction
            return math.log(-math.log1p(-math.exp(log_target_fraction)))
        # This complementary form keeps tiny survivors when u is near one,
        # even if the released fraction itself has rounded to one.
        log_survivor = (-math.inf if self._log_full_depletion is None
                        else -math.exp(self._log_full_depletion))
        target_log_survival = float(np.logaddexp(math.log1p(-probability), log_u + log_survivor))
        return math.log(-target_log_survival)

    def log_quantile(self, probability):
        """Return log(elapsed seconds), retaining times below float64 range."""
        values = _values(probability, "probability")
        if np.any(values <= 0) or np.any(values >= 1):
            raise ValueError("probability must lie strictly inside (0,1)")
        log_dt = math.log(self.dt_s)
        result = np.empty_like(values)
        for location in np.ndindex(values.shape):
            draw = float(values[location])
            if self._uniform_limit:
                log_time = math.fsum((log_dt, math.log(draw)))
            else:
                log_depletion = self._target_log_depletion(draw)
                log_integral = (log_depletion if self._coupling == 0
                                else _log_integral(log_depletion, self._coupling))
                log_time = math.fsum((log_integral, -self._log_rate))
            if not math.isfinite(log_time) or log_time > log_dt + NUMERICAL_LOG_TOLERANCE:
                raise TemkinEmissionTimeRangeError("conditional quantile exceeded its qualified time interval")
            result[location] = min(log_time, log_dt)
        return _output(result)

    def quantile(self, probability):
        """Return positive elapsed seconds; use log_quantile below float range."""
        logs = np.asarray(self.log_quantile(probability), dtype=float)
        with np.errstate(under="ignore"):
            result = np.exp(logs)
        if np.any(result == 0) or not np.all(np.isfinite(result)):
            raise TemkinEmissionTimeRangeError("positive quantile time is below float64 range; use log_quantile")
        # exp(log(dt)) may round one spacing beyond dt; this projection is
        # within the declared representational endpoint tolerance above.
        result = np.minimum(result, self.dt_s)
        if np.any(np.abs(np.log(result) - logs) > NUMERICAL_LOG_TOLERANCE):
            raise TemkinEmissionTimeRangeError(
                "quantile offset is too coarsely resolved in float64; use log_quantile")
        return _output(result)
