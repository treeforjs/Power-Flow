"""Numerically qualified partial windows of a frozen Temkin source law.

The qualification uses monotonic inversion with estimated integral envelopes:
sixteen times the SciPy quadrature error estimate plus explicit float rounding
allowances. These envelopes are numerical assumptions, not certified bounds.
Large-depletion tails use a positive-integral correction so an unresolved early
boundary layer cannot disappear from quadrature. Source windows are not photon
observation gates. This module performs no conservative host transaction.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
import json
import math
from numbers import Real
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad


SCHEMA = "mhdlab-reduced-temkin-emission-window-v1"
QUALIFICATION_IDENTITY = "estimated_monotonic_log_depletion_envelopes_with_positive_tail_correction_v1"
QUADRATURE_ASSUMPTION = "16 times SciPy absolute quadrature error estimate plus explicit float64 rounding allowances; estimated envelopes, not certified bounds"
EPS = float(np.finfo(float).eps)
MIN_POSITIVE = float(np.nextafter(0., 1.))
LOG_MAX = math.log(float(np.finfo(float).max))
MAX_COUPLING = 512.
TAIL_SWITCH = 32.
_LD = np.longdouble


class TemkinWindowRangeError(ValueError):
    """A mathematical window lacks the requested numerical qualification."""


def _real(value, name):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar, not boolean or masked")
    try:
        result = float(value)
    except (OverflowError, ValueError, TypeError) as error:
        raise ValueError(f"{name} must be finite and representable") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _values(value, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot be masked")
    try:
        objects = np.asarray(value, dtype=object)
    except (ValueError, TypeError) as error:
        raise ValueError(f"{name} must contain finite real values") from error
    if objects.ndim > 1 or objects.size == 0:
        raise ValueError(f"{name} must be a scalar or nonempty one-dimensional array")
    for item in objects.flat:
        _real(item, name)
    return np.asarray(objects, dtype=float)


def _output(value):
    if value.ndim == 0:
        return float(value)
    return np.frombuffer(np.ascontiguousarray(value).tobytes(), dtype=float).reshape(value.shape)


def _rounding(*values, multiplier=64.):
    return multiplier*EPS*math.fsum((1., *(abs(value) for value in values)))


def _quadrature(function):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", IntegrationWarning)
            value, error = quad(function, 0., 1., epsabs=0., epsrel=2.e-13, limit=200)
    except (IntegrationWarning, ValueError, OverflowError) as error:
        raise TemkinWindowRangeError("positive integral failed its quadrature estimate") from error
    allowance = 16*error+128*EPS*abs(value)
    if not math.isfinite(value) or not math.isfinite(allowance) or value <= allowance:
        raise TemkinWindowRangeError("positive integral has no resolved estimated lower bound")
    return value, allowance


@lru_cache(maxsize=256)
def _tail_correction(coupling):
    if coupling < 1.e-6:
        # Alternating Ein(c)=c-c²/4+c³/18-...; omitted positive
        # c³/18 bounds this two-term approximation at c<1e-6.
        value = coupling-coupling*coupling/4
        return value, coupling**3/18+128*EPS*coupling
    return _quadrature(lambda t: -math.expm1(-coupling*t)/t if t else coupling)


def _log_integral_estimate(log_s, coupling, coupling_error=0.):
    if not math.isfinite(log_s) or log_s > LOG_MAX:
        raise TemkinWindowRangeError("depletion exceeds the supported finite float64 range")
    s = math.exp(log_s)
    rounding = _rounding(log_s, coupling)
    if coupling == 0 or s == 0:
        # If s underflows, c*s is smaller than rounding allowance; log(s)
        # remains the authoritative positive depletion coordinate.
        return log_s, rounding+coupling_error, "constant_or_subfloat_depletion"
    if s <= TAIL_SWITCH:
        endpoint = coupling*(-math.expm1(-s))
        value, error = _quadrature(lambda v: math.exp(coupling*math.exp(-s*v)*math.expm1(-s*(1-v))))
        estimate = math.fsum((log_s, endpoint, math.log(value)))
        log_error = max(-math.log1p(-error/value), math.log1p(error/value))
        return estimate, log_error+rounding+coupling_error, "endpoint_scaled_positive_quadrature"
    correction, error = _tail_correction(coupling)
    base = s-correction
    if base <= error:
        raise TemkinWindowRangeError("tail integral has no positive estimated denominator")
    tail = coupling*math.exp(-s)
    allowance = max(-math.log1p(-error/base), math.log1p((error+tail)/base))
    return math.fsum((coupling, math.log(base))), allowance+rounding+coupling_error, "linear_tail_with_bounded_positive_correction"


def _solve(coupling, log_rate, duration, *, rate_error=0., coupling_error=0.):
    if duration == 0:
        return {"log_depletion": -math.inf, "depletion": 0., "log_error": 0., "absolute_error": 0., "iterations": 0, "integral_method": "zero_origin"}
    log_duration = math.log(duration)
    log_h = math.fsum((log_rate, log_duration))
    target_error = _rounding(log_rate, log_duration, multiplier=16.)+rate_error
    if not math.isfinite(log_h) or log_h-coupling > LOG_MAX:
        raise TemkinWindowRangeError("log hazard implies unrepresentable depletion")
    left, right = log_h-coupling-target_error-coupling_error, min(LOG_MAX, log_h+target_error)
    if right < left:
        raise TemkinWindowRangeError("no representable estimated depletion bracket")
    if coupling == 0:
        x, radius, iterations, method = log_h, target_error, 0, "constant_barrier_analytic_log_hazard"
    else:
        # Small-s exponential estimate gives a useful stiff-coupling start;
        # the analytic monotonic bounds remain authoritative.
        if log_h+math.log(coupling) < 700:
            argument = math.log1p(coupling*math.exp(log_h))/coupling
            guess = math.log(argument) if argument > 0 else log_h
        else:
            guess = log_h-coupling
        x = min(right, max(left, guess))
        previous_width = math.inf
        for iterations in range(1, 161):
            estimate, error, method = _log_integral_estimate(x, coupling, coupling_error)
            residual = math.fsum((estimate, -log_h))
            radius = abs(residual)+error+target_error+_rounding(x, multiplier=8.)
            new_left, new_right = max(left, x-radius), min(right, x+radius)
            if new_right < new_left:
                raise TemkinWindowRangeError("estimated integral envelopes give inconsistent monotonic brackets")
            left, right = new_left, new_right
            width = right-left
            if width <= 4*(error+target_error+_rounding(x, multiplier=8.)) or width >= previous_width*.999:
                radius = max(x-left, right-x)+_rounding(x, multiplier=8.)
                break
            previous_width = width
            s = math.exp(x)
            derivative = math.exp(min(LOG_MAX, math.fsum((x, coupling*(-math.expm1(-s)), -estimate))))
            proposed = x-residual/derivative
            x = proposed if left < proposed < right else (left+right)/2
        else:
            raise TemkinWindowRangeError("monotonic depletion inversion exceeded its bounded work limit")
    if x > LOG_MAX or not math.isfinite(x):
        raise TemkinWindowRangeError("depletion lacks a finite supported representation")
    value = math.exp(x)
    try:
        absolute_error = value*math.expm1(radius)+16*EPS*value
    except OverflowError as error:
        raise TemkinWindowRangeError("absolute survival uncertainty is unrepresentable") from error
    if value == 0:
        absolute_error = MIN_POSITIVE
    if not math.isfinite(absolute_error):
        raise TemkinWindowRangeError("absolute survival uncertainty is unrepresentable")
    return {"log_depletion": x, "depletion": value, "log_error": radius, "absolute_error": absolute_error,
            "iterations": iterations, "integral_method": method}


def _log_fraction(log_s):
    if log_s < -18:
        s = math.exp(log_s)
        return math.fsum((log_s, math.log1p(-s/2+s*s/6)))
    s = math.exp(log_s)
    return math.log(-math.expm1(-s)) if s <= 1 else math.log1p(-math.exp(-s))


def _fraction_sensitivity(log_s):
    if log_s < -18:
        return 1.
    s = math.exp(log_s)
    if s > 700:
        return 0.
    return s/math.expm1(s)


@dataclass(frozen=True, slots=True, eq=False)
class ReducedTemkinEmissionWindow:
    coupling: float
    log_initial_rate_s: float
    parent_duration_s: float
    start_s: float
    stop_s: float
    log_mass_tolerance: float = field(kw_only=True)
    local_duration_s: float = field(init=False)
    start_depletion: float = field(init=False)
    local_depletion: float = field(init=False)
    parent_depletion: float = field(init=False)
    log_start_depletion: float = field(init=False)
    log_local_depletion: float = field(init=False)
    log_parent_depletion: float = field(init=False)
    log_local_coupling: float = field(init=False)
    local_coupling: float = field(init=False)
    log_local_rate_s: float = field(init=False)
    log_source_fraction: float = field(init=False)
    log_parent_conditional_probability: float = field(init=False)
    _local_fraction_log: float = field(init=False, repr=False)
    _rate_error: float = field(init=False, repr=False)
    _coupling_error: float = field(init=False, repr=False)
    _query_base_error: float = field(init=False, repr=False)
    _diagnostics_json: str = field(init=False, repr=False)

    def __post_init__(self):
        for name in ("coupling", "log_initial_rate_s", "parent_duration_s", "start_s", "stop_s", "log_mass_tolerance"):
            object.__setattr__(self, name, _real(getattr(self, name), name))
        if not 0 <= self.coupling <= MAX_COUPLING:
            raise ValueError("coupling must lie in the qualified range [0,512]")
        if not 0 <= self.start_s < self.stop_s <= self.parent_duration_s or self.parent_duration_s <= 0:
            raise ValueError("source windows require 0<=start_s<stop_s<=positive parent_duration_s; empty windows are unsupported")
        if self.log_mass_tolerance <= 0:
            raise ValueError("log_mass_tolerance must be explicit and positive")
        duration_ld = _LD(self.stop_s)-_LD(self.start_s)
        duration = float(duration_ld)
        if not math.isfinite(duration) or duration <= 0:
            raise TemkinWindowRangeError("local duration is not positively representable")
        duration_error = float(abs(_LD(duration)-duration_ld)/duration_ld)
        start = _solve(self.coupling, self.log_initial_rate_s, self.start_s)
        parent = _solve(self.coupling, self.log_initial_rate_s, self.parent_duration_s)
        sa = start["depletion"]
        log_c = math.log(self.coupling)-sa if self.coupling else -math.inf
        ca = math.exp(log_c) if self.coupling else 0.
        dropped = self.coupling > 0 and ca == 0
        approximation_bound = MIN_POSITIVE if dropped else 0.
        ca_error = max(MIN_POSITIVE, 64*EPS*ca) if self.coupling else 0.
        ra = math.fsum((self.log_initial_rate_s, -self.coupling*(-math.expm1(-sa))))
        rate_error = _rounding(self.log_initial_rate_s, self.coupling)+duration_error
        local = _solve(ca, ra, duration, rate_error=rate_error, coupling_error=ca_error)
        local_log_fraction = _log_fraction(local["log_depletion"])
        parent_log_fraction = _log_fraction(parent["log_depletion"])
        log_mass = math.fsum((-sa, local_log_fraction))
        log_probability = math.fsum((log_mass, -parent_log_fraction))
        local_fraction_error = local["log_error"]*_fraction_sensitivity(local["log_depletion"]-local["log_error"])+_rounding(local_log_fraction)
        # Coherent rebase sensitivity: |d log(window mass)/d s_a|<=1+c.
        start_effect = (1+self.coupling)*start["absolute_error"]
        mass_error = math.fsum((start_effect, local_fraction_error, approximation_bound, _rounding(sa, local_log_fraction)))
        probability_error = mass_error+parent["log_error"]*_fraction_sensitivity(parent["log_depletion"]-parent["log_error"])+_rounding(parent_log_fraction, log_probability)
        if (not math.isfinite(log_mass) or not math.isfinite(log_probability)
                or max(mass_error, probability_error) > self.log_mass_tolerance):
            raise TemkinWindowRangeError("estimated absolute log-mass/survival uncertainty exceeds log_mass_tolerance")
        if log_probability > probability_error:
            raise TemkinWindowRangeError("estimated window probability exceeds the parent emitted population")
        properties = {
            "local_duration_s": duration, "start_depletion": sa, "local_depletion": local["depletion"], "parent_depletion": parent["depletion"],
            "log_start_depletion": start["log_depletion"], "log_local_depletion": local["log_depletion"], "log_parent_depletion": parent["log_depletion"],
            "log_local_coupling": log_c, "local_coupling": ca, "log_local_rate_s": ra,
            "log_source_fraction": log_mass, "log_parent_conditional_probability": min(0., log_probability),
            "_local_fraction_log": local_log_fraction, "_rate_error": rate_error, "_coupling_error": ca_error,
            "_query_base_error": mass_error,
        }
        for name, value in properties.items():
            object.__setattr__(self, name, value)
        root_records = {name: {key: (None if not math.isfinite(value) else value) if isinstance(value, float) else value for key, value in record.items()}
                        for name, record in (("start", start), ("local", local), ("parent", parent))}
        diagnostics = {
            "qualification_identity": QUALIFICATION_IDENTITY, "quadrature_assumption": QUADRATURE_ASSUMPTION,
            "log_mass_tolerance": self.log_mass_tolerance,
            "estimated_log_source_fraction_error": mass_error, "estimated_log_parent_probability_error": probability_error,
            "start_depletion_error": start["absolute_error"], "local_log_depletion_error": local["log_error"], "parent_log_depletion_error": parent["log_error"],
            "local_coupling_approximated_as_zero": dropped,
            "log_local_coupling_approximation_error_bound": log_c if dropped else None,
            "local_coupling_approximation_error_bound": approximation_bound,
            "start_rebase_log_mass_sensitivity_bound": 1+self.coupling,
            "local_duration_relative_rounding_error": duration_error,
            "local_depletion_underflowed": local["depletion"] == 0,
            "start_depletion_underflowed": self.start_s > 0 and sa == 0,
            "root_evidence": root_records, "tail_switch_depletion": TAIL_SWITCH,
            "certified_error_bound": False, "host_budget_accepted": False,
        }
        object.__setattr__(self, "_diagnostics_json", json.dumps(diagnostics, allow_nan=False))

    @property
    def numerical_diagnostics(self):
        return json.loads(self._diagnostics_json)

    def manifest(self):
        values = {name: getattr(self, name) for name in (
            "coupling", "log_initial_rate_s", "parent_duration_s", "start_s", "stop_s", "local_duration_s",
            "start_depletion", "local_depletion", "parent_depletion", "log_start_depletion", "log_local_depletion",
            "log_parent_depletion", "log_local_coupling", "local_coupling", "log_local_rate_s",
            "log_source_fraction", "log_parent_conditional_probability", "log_mass_tolerance")}
        return {"schema": SCHEMA, **{name: None if not math.isfinite(value) else value for name, value in values.items()},
                "time_coordinate": "local delay from original source-window start; keep start and delay separately",
                "numerical_diagnostics": self.numerical_diagnostics}

    def log_cdf(self, delay_s):
        values = _values(delay_s, "delay_s")
        if np.any(values < 0) or np.any(values > self.local_duration_s):
            raise ValueError("delay_s must lie within the local source window")
        result = np.empty_like(values)
        for location in np.ndindex(values.shape):
            delay = float(values[location])
            if delay == 0:
                result[location] = -math.inf
            elif delay == self.local_duration_s:
                result[location] = 0.
            else:
                trial = _solve(self.local_coupling, self.log_local_rate_s, delay, rate_error=self._rate_error, coupling_error=self._coupling_error)
                log_fraction = _log_fraction(trial["log_depletion"])
                error = 2*self._query_base_error+trial["log_error"]+_rounding(log_fraction, self._local_fraction_log)
                if error > self.log_mass_tolerance:
                    raise TemkinWindowRangeError("conditional CDF lacks the requested estimated log-error qualification")
                probability = math.fsum((log_fraction, -self._local_fraction_log))
                if probability > error:
                    raise TemkinWindowRangeError("conditional CDF exceeds its estimated local envelope")
                result[location] = min(0., probability)
        return _output(result)

    def cdf(self, delay_s):
        logs = np.asarray(self.log_cdf(delay_s), dtype=float)
        with np.errstate(under="ignore"):
            values = np.exp(logs)
        return _output(values)

    def log_quantile(self, probability):
        values = _values(probability, "probability")
        if np.any(values <= 0) or np.any(values >= 1):
            raise ValueError("probability must lie strictly inside (0,1)")
        result = np.empty_like(values)
        for location in np.ndindex(values.shape):
            u = float(values[location])
            log_target_fraction = math.fsum((math.log(u), self._local_fraction_log))
            if log_target_fraction < -18:
                q = math.exp(log_target_fraction)
                log_target = math.fsum((log_target_fraction, math.log1p(q/2+q*q/3)))
            else:
                log_survival = float(np.logaddexp(math.log1p(-u), math.log(u)-self.local_depletion))
                log_target = math.log(-log_survival)
            estimate, integral_error, _ = _log_integral_estimate(log_target, self.local_coupling, self._coupling_error)
            target = math.exp(log_target)
            derivative = math.exp(min(LOG_MAX, math.fsum((log_target, self.local_coupling*(-math.expm1(-target)), -estimate))))
            amplification = math.expm1(target)/target if target else 1.
            error = self._query_base_error*(1+derivative*amplification)+integral_error+self._rate_error+_rounding(estimate, self.log_local_rate_s)
            if not math.isfinite(error) or error > self.log_mass_tolerance:
                raise TemkinWindowRangeError("conditional quantile lacks the requested estimated log-time qualification")
            log_delay = math.fsum((estimate, -self.log_local_rate_s))
            if not math.isfinite(log_delay) or log_delay > math.log(self.local_duration_s)+error:
                raise TemkinWindowRangeError("conditional quantile exceeds its estimated window envelope")
            result[location] = min(math.log(self.local_duration_s), log_delay)
        return _output(result)

    def quantile(self, probability):
        logs = np.asarray(self.log_quantile(probability), dtype=float)
        with np.errstate(under="ignore", over="ignore"):
            values = np.exp(logs)
        if np.any(values == 0) or np.any(~np.isfinite(values)):
            raise TemkinWindowRangeError("local quantile is not positively representable; use log_quantile")
        if np.any(np.abs(np.log(values)-logs) > self.log_mass_tolerance):
            raise TemkinWindowRangeError("local quantile float64 quantization exceeds the requested log-time tolerance")
        values = np.minimum(values, self.local_duration_s)
        return _output(values)
