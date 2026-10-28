"""Exact linear-ramp magnetic-diffusion reference, with finite wall heating.

A stationary uniform conductor fills inward depth y>=0, initially B=0.  Its
tangential boundary field is B_s(t)=R*t for t>=0, with no ramp stop.  The caller
supplies constant electrical conductivity sigma, permeability mu and volumetric
heat capacity Cv, all in SI units.  There is no flow, thermal conduction, Hall
term, displacement current, radiation, phase change or thermodynamic closure.
This reference does not run FLASH or validate material properties.

Write D=1/(mu*sigma), u=y/(2*sqrt(D*t)),
    g(u)=exp(-u²)/sqrt(pi)-u*erfc(u),
    b(u)=(1+2u²)*erfc(u)-2u*exp(-u²)/sqrt(pi).
Then B_x=R*t*b, J_z=2*R*t/(mu*sqrt(D*t))*g and q=4*R²*t/mu*g².
The current sign follows J_z=-dB_x/dy/mu.  Cell q integrates J²/sigma,
independently of the current mean.  Local adiabatic temperature is the time
integral of q/Cv and includes y=0: DeltaT_wall=2*R²*t²/(mu*pi*Cv).

Global energies are per unit half-space face area.  For S=R²/mu*sqrt(D)*t^(5/2),
    U_B=(28-16*sqrt(2))/(15*sqrt(pi))*S,
    Q_J=16*(sqrt(2)-1)/(15*sqrt(pi))*S,
    W_in=4/(5*sqrt(pi))*S.
The instantaneous incoming surface power is 2*R²/mu*sqrt(D/pi)*t^(3/2).

Numerical qualification: near-surface closed forms and positive Gaussian
moment integrals avoid cancellation in b/g tails.  Two positive Gauss-Laguerre
rules must agree within the declared tolerance.  Adaptive positive quadrature
is used for cell Joule averages and adiabatic temperature.  In its scaled tail
coordinate, the omitted integral after 50 is at most exp(-198)/4 (less than
4e-86 relative); this is a quadrature bound, not a physical field or rate cutoff.
Unrepresentable outputs or unverified quadrature raise explicitly; only
mathematically unrepresentable float64 tails round to zero.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import warnings

import numpy as np
from scipy import integrate, special

from .heating_benchmarks import (
    HalfSpaceCellAverages, HalfSpaceEnergyBudget, HalfSpaceFields,
    _depths, _exp, _immutable, _real,
)


KERNEL_LOG_TOLERANCE = 5e-12
QUADRATURE_RELATIVE_TOLERANCE = 2e-11
_LOG_TWO = math.log(2.)
_LOG_PI = math.log(math.pi)
_SQRT_PI = math.sqrt(math.pi)
_LOG_HALF_SUBNORMAL = math.log(float(np.nextafter(0., 1.))) - _LOG_TWO
_LEGENDRE_NODES, _LEGENDRE_WEIGHTS = np.polynomial.legendre.leggauss(16)
_LAGUERRE_RULES = []
for _order in (64, 96):
    _nodes, _weights = special.roots_laguerre(_order)
    _LAGUERRE_RULES.append((_nodes, np.log(_weights), np.log(_nodes)))


class SmoothHeatingNumericalError(ValueError):
    """The requested reference output is outside verified numerical accuracy."""


def _quad(function, lower: float, upper: float) -> float:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", integrate.IntegrationWarning)
            value, error = integrate.quad(function, lower, upper, epsabs=0., epsrel=2e-12, limit=200)
    except (integrate.IntegrationWarning, ValueError, OverflowError) as exc:
        raise SmoothHeatingNumericalError("smooth heating quadrature did not meet its accuracy criteria") from exc
    if not math.isfinite(value) or value <= 0 or error > QUADRATURE_RELATIVE_TOLERANCE * value:
        raise SmoothHeatingNumericalError("smooth heating quadrature did not meet its accuracy criteria")
    return value


def _log_kernels(coordinate: float) -> tuple[float, float, float]:
    """Log(g), log(b), log integral_u^infinity b(v)dv, all positive."""
    if not math.isfinite(coordinate):
        return -math.inf, -math.inf, -math.inf
    if coordinate < 2.:
        scaled_erfc = float(special.erfcx(coordinate))
        square = coordinate * coordinate
        g = 1 / _SQRT_PI - coordinate * scaled_erfc
        b = (1 + 2 * square) * scaled_erfc - 2 * coordinate / _SQRT_PI
        tail_b = (2 * (1 + square) / (3 * _SQRT_PI)
                  - (coordinate + 2 * coordinate**3 / 3) * scaled_erfc)
        if min(g, b, tail_b) <= 0:
            raise SmoothHeatingNumericalError("positive smooth heating kernel lost numerical resolution")
        return tuple(-square + math.log(value) for value in (g, b, tail_b))
    # Positive representation L_n=integral_0^infinity w^n exp(-w²-2uw)dw;
    # substitute x=2uw.  g=2 exp(-u²)L1/sqrt(pi), b=4 exp(-u²)L2/sqrt(pi),
    # tail_b=4 exp(-u²)L3/(3sqrt(pi)). No subtraction of tail-sized terms.
    moments = []
    for nodes, log_weights, log_nodes in _LAGUERRE_RULES:
        correction = -0.25 * (nodes / coordinate)**2
        moments.append(np.array([special.logsumexp(log_weights + power * log_nodes + correction)
                                 for power in (1, 2, 3)]))
    if np.max(np.abs(moments[0] - moments[1])) > KERNEL_LOG_TOLERANCE:
        raise SmoothHeatingNumericalError("positive smooth heating tail rules did not agree")
    square = coordinate * coordinate
    log_u = math.log(coordinate)
    prefactors = (-_LOG_TWO - .5 * _LOG_PI - 2 * log_u,
                  -_LOG_TWO - .5 * _LOG_PI - 3 * log_u,
                  -math.log(12.) - .5 * _LOG_PI - 4 * log_u)
    return tuple(-square + prefactor + moment for prefactor, moment in zip(prefactors, moments[1]))


def _log_difference(first: float, last: float) -> float:
    if last == -math.inf:
        return first
    if last >= first:
        raise SmoothHeatingNumericalError("positive smooth heating cell integral lost numerical resolution")
    return first + math.log(-math.expm1(last - first))


def _validate_edges(value):
    edges = _depths(value)
    if edges.ndim != 1 or edges.size < 2 or np.any(edges[1:] <= edges[:-1]):
        raise ValueError("depth_edges_m must be a strictly increasing 1-D array with at least two edges")
    return edges


@dataclass(frozen=True, slots=True)
class LinearRampConductivityHalfSpace:
    """Caller-supplied constant-property linear-ramp reference in SI units.

    All four physical parameters are required. The signed ramp rate is dB_s/dt
    in T/s. Cv is heat capacity per volume, not specific heat per mass. Time
    and depth are nonnegative; t=0 returns the initial zero field/heating state.
    """

    conductivity_s_m: float
    permeability_h_m: float
    volumetric_heat_capacity_j_m3_k: float
    boundary_field_rate_t_s: float

    def __post_init__(self):
        for name in ("conductivity_s_m", "permeability_h_m", "volumetric_heat_capacity_j_m3_k"):
            object.__setattr__(self, name, _real(getattr(self, name), name, positive=True))
        object.__setattr__(self, "boundary_field_rate_t_s", _real(self.boundary_field_rate_t_s, "boundary_field_rate_t_s"))
        if self.magnetic_diffusivity_m2_s == 0:
            raise ValueError("magnetic diffusivity is below positive float64 range")

    @property
    def _log_diffusivity(self):
        return -math.log(self.permeability_h_m) - math.log(self.conductivity_s_m)

    @property
    def magnetic_diffusivity_m2_s(self):
        return float(_exp(self._log_diffusivity))

    def _time(self, value):
        time = _real(value, "time_s")
        if time < 0:
            raise ValueError("time_s must be nonnegative")
        return time

    def _scale(self, time):
        log_scale = .5 * (self._log_diffusivity + math.log(time))
        scale = float(_exp(log_scale))
        if scale == 0:
            raise ValueError("diffusion length is below positive float64 range")
        return scale, log_scale

    def fields(self, depth_m, time_s) -> HalfSpaceFields:
        """Point B_x [T], J_z [A/m²] and Joule power [W/m³]."""
        depths, time = _depths(depth_m), self._time(time_s)
        magnetic, current, power = (np.zeros_like(depths) for _ in range(3))
        if time > 0 and self.boundary_field_rate_t_s != 0:
            scale, log_scale = self._scale(time)
            log_rate = math.log(abs(self.boundary_field_rate_t_s))
            log_mu, log_t = math.log(self.permeability_h_m), math.log(time)
            sign = math.copysign(1., self.boundary_field_rate_t_s)
            for location in np.ndindex(depths.shape):
                with np.errstate(over="ignore"):
                    coordinate = float(.5 * (depths[location] / scale))
                log_g, log_b, _ = _log_kernels(coordinate)
                magnetic[location] = _exp(log_rate + log_t + log_b, sign=sign)
                current[location] = _exp(_LOG_TWO + log_rate + log_t - log_mu - log_scale + log_g, sign=sign)
                power[location] = _exp(math.log(4.) + 2 * log_rate + log_t - log_mu + 2 * log_g)
        return HalfSpaceFields(*map(_immutable, (depths, magnetic, current, power)))

    def _cell_power_log_mean(self, lower, upper, log_width, log_scale, log_prefactor):
        log_g = _log_kernels(lower)[0]
        log_maximum = log_prefactor + 2 * log_g
        if log_maximum < _LOG_HALF_SUBNORMAL:
            return -math.inf
        log_scaled_width = log_width - _LOG_TWO - log_scale + math.log1p(lower)
        if log_scaled_width <= 0:
            dimensionless_width = math.exp(log_width - _LOG_TWO - log_scale)

            def normalized(location):
                log_value = _log_kernels(lower + dimensionless_width * location)[0]
                return math.exp(2 * (log_value - log_g))

            return log_maximum + math.log(_quad(normalized, 0., 1.))
        maximum = 50. if log_scaled_width >= math.log(50.) else math.exp(log_scaled_width)

        def transformed(location):
            log_value = _log_kernels(lower + location / (1 + lower))[0]
            return math.exp(2 * (log_value - log_g))

        integral = _quad(transformed, 0., maximum)
        return log_maximum + math.log(integral) - log_scaled_width

    def cell_averages(self, depth_edges_m, time_s) -> HalfSpaceCellAverages:
        """Finite-volume B/J/q means; q integrates J² independently of mean J."""
        edges, time = _validate_edges(depth_edges_m), self._time(time_s)
        values = np.zeros((3, edges.size - 1))
        if time > 0 and self.boundary_field_rate_t_s != 0:
            scale, log_scale = self._scale(time)
            log_rate, log_t = math.log(abs(self.boundary_field_rate_t_s)), math.log(time)
            log_mu, sign = math.log(self.permeability_h_m), math.copysign(1., self.boundary_field_rate_t_s)
            log_power_prefactor = math.log(4.) + 2 * log_rate + log_t - log_mu
            for index, (low_y, high_y) in enumerate(zip(edges[:-1], edges[1:])):
                width = float(high_y - low_y)
                log_width = math.log(width)
                with np.errstate(over="ignore"):
                    lower, upper = float(.5 * (low_y / scale)), float(.5 * (high_y / scale))
                _, b_first, tail_first = _log_kernels(lower)
                _, b_last, tail_last = _log_kernels(upper)
                if tail_first == tail_last or tail_last - tail_first > -1e-5:
                    nodes = low_y + .5 * width * (_LEGENDRE_NODES + 1)
                    with np.errstate(over="ignore"):
                        coordinates = .5 * (nodes / scale)
                    kernels = np.array([_log_kernels(float(value))[:2] for value in coordinates])
                    weights = np.log(.5 * _LEGENDRE_WEIGHTS)
                    log_b_mean = float(special.logsumexp(weights + kernels[:, 1]))
                    log_g_mean = float(special.logsumexp(weights + kernels[:, 0]))
                    b_log = log_rate + log_t + log_b_mean
                    j_log = _LOG_TWO + log_rate + log_t - log_mu - log_scale + log_g_mean
                else:
                    b_log = (log_rate + log_t + _LOG_TWO + log_scale - log_width
                             + _log_difference(tail_first, tail_last))
                    j_log = log_rate + log_t - log_mu - log_width + _log_difference(b_first, b_last)
                q_log = self._cell_power_log_mean(lower, upper, log_width, log_scale, log_power_prefactor)
                values[:, index] = [float(_exp(b_log, sign=sign)), float(_exp(j_log, sign=sign)),
                                    float(_exp(q_log))]
        return HalfSpaceCellAverages(*map(_immutable, (edges, *values)))

    def energy_budget(self, time_s) -> HalfSpaceEnergyBudget:
        """Global magnetic/Joule/input energies [J/m²] and surface power [W/m²]."""
        time = self._time(time_s)
        if time == 0 or self.boundary_field_rate_t_s == 0:
            return HalfSpaceEnergyBudget(0., 0., 0., 0., 0.)
        log_common = (2 * math.log(abs(self.boundary_field_rate_t_s)) - math.log(self.permeability_h_m)
                      + .5 * self._log_diffusivity + 2.5 * math.log(time))
        coefficients = ((28 - 16 * math.sqrt(2)) / (15 * _SQRT_PI),
                        16 * (math.sqrt(2) - 1) / (15 * _SQRT_PI), 4 / (5 * _SQRT_PI))
        magnetic, joule, incoming = (float(_exp(log_common + math.log(value))) for value in coefficients)
        power = float(_exp(log_common - math.log(time) + math.log(2 / _SQRT_PI)))
        residual = math.fsum((incoming, -magnetic, -joule))
        return HalfSpaceEnergyBudget(magnetic, joule, incoming, power, residual)

    def adiabatic_temperature_rise(self, depth_m, time_s) -> np.ndarray:
        """Finite local DeltaT [K], integrating q/Cv without thermal conduction."""
        depths, time = _depths(depth_m), self._time(time_s)
        result = np.zeros_like(depths)
        if time == 0 or self.boundary_field_rate_t_s == 0:
            return _immutable(result)
        scale, _ = self._scale(time)
        log_prefactor = (math.log(4.) + 2 * math.log(abs(self.boundary_field_rate_t_s))
                         + 2 * math.log(time) - math.log(self.permeability_h_m)
                         - math.log(self.volumetric_heat_capacity_j_m3_k))
        for location in np.ndindex(depths.shape):
            with np.errstate(over="ignore"):
                coordinate = float(.5 * (depths[location] / scale))
            if coordinate == 0:
                log_h = -math.log(2 * math.pi)
            else:
                log_g = _log_kernels(coordinate)[0]
                if log_prefactor + 2 * log_g - _LOG_TWO < _LOG_HALF_SUBNORMAL:
                    continue
                if coordinate < 1:
                    def integrand(fraction):
                        log_value = _log_kernels(coordinate / math.sqrt(fraction))[0]
                        return fraction * math.exp(2 * (log_value - log_g))

                    integral = _quad(integrand, 0., 1.)
                    log_h = 2 * log_g + math.log(integral)
                else:
                    def integrand(transformed):
                        value = coordinate + transformed / (1 + coordinate)
                        log_value = _log_kernels(value)[0]
                        return math.exp(2 * (log_value - log_g) + 5 * math.log(coordinate / value))

                    integral = _quad(integrand, 0., 50.)
                    log_h = (_LOG_TWO - math.log(coordinate) - math.log1p(coordinate)
                             + 2 * log_g + math.log(integral))
            result[location] = _exp(log_prefactor + log_h)
        return _immutable(result)
