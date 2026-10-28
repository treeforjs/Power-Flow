"""Analytic, constant-property magnetic-diffusion verification in SI units.

This is a benchmark, not a FLASH driver or a heating solver.  A stationary,
uniform conductor fills inward depth y >= 0; its initial field is zero and a
constant tangential field B_s is imposed at y = 0 for t > 0.  Conductivity,
permeability and volumetric heat capacity are constant.  There is no fluid
motion, thermal conduction, displacement current, radiation or phase change.

With magnetic diffusivity D = 1/(mu*sigma),
    B_x = B_s erfc[y/(2 sqrt(D t))],
    J_z = -dB_x/dy / mu,
    q_J = J_z**2/sigma.
The depth direction fixes the sign of J_z.  Signed B_s is allowed.  All material
parameters must be supplied; this module provides no metal-property defaults.

Integrated magnetic energy, accumulated Joule heat and incoming Poynting
energy are per unit surface area of the *infinite* half-space.  Temperature
rise is a local adiabatic diagnostic with constant volumetric heat capacity;
it is not a conductive thermal solution.  Its surface value diverges for the
ideal instantaneous field step, so y = 0 is explicitly excluded for that API.

Cell averages integrate B, J and J**2/sigma separately.  In particular the
Joule average is not the square of the current average.  Closed antiderivatives
are evaluated with logarithmic differences; narrow-bin limiting quadrature
avoids cancellation.  Float64 underflow at unrepresentably small tails is a
representation limit, not a physical field or heating cutoff.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Real

import numpy as np
from scipy import special


_LOG_MAX = math.log(np.finfo(float).max)
_LOG_PI = math.log(math.pi)
_LOG_TWO = math.log(2.0)
_SQRT_PI = math.sqrt(math.pi)
_EULER_GAMMA = 0.5772156649015328606
_GAUSS_NODES, _GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(16)


def _real(value, name: str, *, positive: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number, not bool or complex")
    try:
        value = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(value) or (positive and value <= 0):
        suffix = " and positive" if positive else ""
        raise ValueError(f"{name} must be finite{suffix}")
    return value


def _depths(value, *, positive: bool = False) -> np.ndarray:
    if np.ma.isMaskedArray(value):
        raise ValueError("depth_m cannot contain masked coordinates")
    try:
        raw = np.asarray(value, dtype=object)
    except (ValueError, TypeError) as exc:
        raise ValueError("depth_m must be a nonempty finite real array or scalar") from exc
    if raw.size == 0:
        raise ValueError("depth_m must be nonempty")
    for item in raw.flat:
        number = _real(item, "depth_m", positive=positive)
        if number < 0:
            raise ValueError("depth_m must be nonnegative")
    return np.asarray(raw, dtype=float)


def _immutable(value) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    return np.frombuffer(np.ascontiguousarray(array).tobytes(), dtype=float).reshape(array.shape)


def _exp(log_value, *, sign: float = 1.0):
    logs = np.asarray(log_value, dtype=float)
    if np.any(np.isnan(logs)) or np.any(logs > _LOG_MAX):
        raise ValueError("benchmark result exceeds finite float64 range")
    with np.errstate(under="ignore"):
        result = sign * np.exp(logs)
    return result


def _log_erfc(value):
    # erfc(x) = 2 Phi(-sqrt(2)*x); log_ndtr retains Gaussian tails.
    with np.errstate(over="ignore"):
        return _LOG_TWO + special.log_ndtr(-math.sqrt(2.0) * np.asarray(value))


def _log_difference(log_first: float, log_second: float) -> float:
    if log_second == -math.inf:
        return log_first
    if log_first == log_second:
        return -math.inf
    if log_second > log_first:
        raise ValueError("roundoff prevented evaluation of a positive cell integral")
    return log_first + math.log(-math.expm1(log_second - log_first))


def _log_erfc_tail_integral(value: float) -> float:
    """log integral from value to infinity of erfc(u) du, value >= 0."""
    if not math.isfinite(value):
        return -math.inf
    if value > 1e4:
        # Here even the largest finite B_s times erfc(value) underflows.  This
        # asymptotic expression avoids cancellation in 1-sqrt(pi)*u*erfcx(u).
        return -value * value - math.log(2 * _SQRT_PI) - 2 * math.log(value)
    complement = 1.0 - _SQRT_PI * value * special.erfcx(value)
    return -value * value - 0.5 * _LOG_PI + math.log(complement)


def _log_exp1_from_log_argument(log_argument: np.ndarray) -> np.ndarray:
    """Retain E1 tails when a large dimensional coefficient can restore them."""
    result = np.empty_like(log_argument)
    small = log_argument < -36.0
    # E1(z) = -gamma-log(z)+O(z).  The omitted correction in this branch
    # is below float64 spacing of the retained positive value.
    result[small] = np.log(-_EULER_GAMMA - log_argument[small])
    ordinary = (~small) & (log_argument <= math.log(50.0))
    result[ordinary] = np.log(special.exp1(np.exp(log_argument[ordinary])))
    for index in np.ndindex(log_argument.shape):
        if small[index] or ordinary[index]:
            continue
        if log_argument[index] > _LOG_MAX:
            result[index] = -math.inf
            continue
        argument = math.exp(float(log_argument[index]))
        # E1(z) = exp(-z)/z [1-1/z+2!/z²-...].  For z>50 the
        # optimally truncated asymptotic sum reaches float64 precision well
        # before its terms grow.  Evaluate the exponential only at the end.
        term, series = 1.0, 1.0
        for order in range(1, 128):
            next_term = term * (-order / argument)
            if abs(next_term) >= abs(term):
                break
            series += next_term
            term = next_term
            if abs(term) <= np.finfo(float).eps * abs(series):
                break
        result[index] = -argument - math.log(argument) + math.log(series)
    return result


@dataclass(frozen=True, slots=True, eq=False)
class HalfSpaceFields:
    depth_m: np.ndarray
    magnetic_field_t: np.ndarray
    current_density_a_m2: np.ndarray
    joule_power_w_m3: np.ndarray


@dataclass(frozen=True, slots=True, eq=False)
class HalfSpaceCellAverages:
    depth_edges_m: np.ndarray
    magnetic_field_t: np.ndarray
    current_density_a_m2: np.ndarray
    joule_power_w_m3: np.ndarray


@dataclass(frozen=True, slots=True)
class HalfSpaceEnergyBudget:
    magnetic_energy_j_m2: float
    joule_heat_j_m2: float
    poynting_input_j_m2: float
    surface_power_w_m2: float
    energy_residual_j_m2: float


@dataclass(frozen=True, slots=True)
class ConstantConductivityHalfSpace:
    """Exact magnetic-diffusion reference with caller-supplied SI parameters.

    ``volumetric_heat_capacity_j_m3_k`` is C_v per volume, not specific heat
    per mass.  Time is strictly positive; the imposed step has no finite
    surface-current or surface-power value at t = 0.  Depths can be scalars
    or nonempty arrays; returned arrays preserve their shape and are immutable.
    """

    conductivity_s_m: float
    permeability_h_m: float
    volumetric_heat_capacity_j_m3_k: float
    boundary_field_t: float

    def __post_init__(self) -> None:
        for name in ("conductivity_s_m", "permeability_h_m", "volumetric_heat_capacity_j_m3_k"):
            object.__setattr__(self, name, _real(getattr(self, name), name, positive=True))
        object.__setattr__(self, "boundary_field_t", _real(self.boundary_field_t, "boundary_field_t"))
        # Validate representability without multiplying mu*sigma first.
        diffusivity = float(_exp(self._log_diffusivity))
        if diffusivity == 0:
            raise ValueError("magnetic diffusivity is below positive float64 range")

    @property
    def _log_diffusivity(self) -> float:
        return -math.log(self.permeability_h_m) - math.log(self.conductivity_s_m)

    @property
    def magnetic_diffusivity_m2_s(self) -> float:
        return float(_exp(self._log_diffusivity))

    def _scale(self, time_s) -> tuple[float, float, float]:
        time = _real(time_s, "time_s", positive=True)
        log_scale = 0.5 * (self._log_diffusivity + math.log(time))
        scale = float(_exp(log_scale))
        if scale == 0:
            raise ValueError("diffusion length is below positive float64 range")
        return time, scale, log_scale

    def fields(self, depth_m, time_s) -> HalfSpaceFields:
        """Return point B_x [T], J_z [A/m²] and instantaneous q_J [W/m³]."""
        depths = _depths(depth_m)
        time, scale, log_scale = self._scale(time_s)
        with np.errstate(over="ignore"):
            coordinate = 0.5 * (depths / scale)
            square = coordinate * coordinate
        if self.boundary_field_t == 0:
            magnetic = current = heating = np.zeros_like(depths)
        else:
            log_b = math.log(abs(self.boundary_field_t))
            sign = math.copysign(1.0, self.boundary_field_t)
            magnetic = _exp(log_b + _log_erfc(coordinate), sign=sign)
            current = _exp(log_b - math.log(self.permeability_h_m)
                           - 0.5 * _LOG_PI - log_scale - square, sign=sign)
            heating = _exp(2 * log_b - math.log(self.permeability_h_m)
                           - _LOG_PI - math.log(time) - 2 * square)
        return HalfSpaceFields(*map(_immutable, (depths, magnetic, current, heating)))

    def cell_averages(self, depth_edges_m, time_s) -> HalfSpaceCellAverages:
        """Return finite-volume averages over contiguous increasing depth bins.

        Analytic Gaussian and erfc integrals handle broad bins.  A 16-point
        Gauss rule over an individual bin evaluates only the cancellation-prone
        narrow-bin limit.  No data are sampled or fitted to a FLASH result.
        """
        edges = _depths(depth_edges_m)
        if edges.ndim != 1 or edges.size < 2 or np.any(edges[1:] <= edges[:-1]):
            raise ValueError("depth_edges_m must be a strictly increasing 1-D array with at least two edges")
        time, scale, log_scale = self._scale(time_s)
        values = np.zeros((3, edges.size - 1))
        if self.boundary_field_t != 0:
            log_b = math.log(abs(self.boundary_field_t))
            log_mu = math.log(self.permeability_h_m)
            sign = math.copysign(1.0, self.boundary_field_t)
            for index, (lower, upper) in enumerate(zip(edges[:-1], edges[1:])):
                width = upper - lower
                with np.errstate(over="ignore"):
                    first = float(0.5 * (lower / scale))
                    last = float(0.5 * (upper / scale))
                log_first = _log_erfc_tail_integral(first)
                log_last = _log_erfc_tail_integral(last)
                # The positive antiderivative difference loses precision when
                # its operands almost coincide.  Evaluate that local mean
                # directly on the bin, retaining logarithmic Gaussian tails.
                narrow = log_first == log_last or log_last - log_first > -1e-5
                if narrow:
                    nodes = lower + 0.5 * width * (_GAUSS_NODES + 1)
                    with np.errstate(over="ignore"):
                        scaled = 0.5 * (nodes / scale)
                        square = scaled * scaled
                    log_weights = np.log(0.5 * _GAUSS_WEIGHTS)
                    log_means = [special.logsumexp(log_weights + _log_erfc(scaled)),
                                 special.logsumexp(log_weights - square),
                                 special.logsumexp(log_weights - 2 * square)]
                    logs = [log_b + log_means[0],
                            log_b - log_mu - 0.5 * _LOG_PI - log_scale + log_means[1],
                            2 * log_b - log_mu - _LOG_PI - math.log(time) + log_means[2]]
                else:
                    log_width = math.log(width)
                    log_erfc_difference = _log_difference(float(_log_erfc(first)), float(_log_erfc(last)))
                    log_q_difference = _log_difference(float(_log_erfc(math.sqrt(2) * first)),
                                                       float(_log_erfc(math.sqrt(2) * last)))
                    logs = [log_b + _LOG_TWO + log_scale - log_width
                            + _log_difference(log_first, log_last),
                            log_b - log_mu - log_width + log_erfc_difference,
                            2 * log_b - log_mu - _LOG_PI - math.log(time)
                            + log_scale + 0.5 * (_LOG_PI - _LOG_TWO) - log_width + log_q_difference]
                values[:, index] = [float(_exp(logs[0], sign=sign)),
                                    float(_exp(logs[1], sign=sign)), float(_exp(logs[2]))]
        return HalfSpaceCellAverages(*map(_immutable, (edges, *values)))

    def energy_budget(self, time_s) -> HalfSpaceEnergyBudget:
        """Exact global energies [J/m²] and instantaneous surface power [W/m²]."""
        time, _, log_scale = self._scale(time_s)
        if self.boundary_field_t == 0:
            return HalfSpaceEnergyBudget(0., 0., 0., 0., 0.)
        log_base = 2 * math.log(abs(self.boundary_field_t)) - math.log(self.permeability_h_m)
        base = float(_exp(log_base + log_scale - 0.5 * _LOG_PI))
        magnetic = base * (2 - math.sqrt(2))
        joule = base * math.sqrt(2)
        incoming = float(_exp(math.log(2) + log_base + log_scale - 0.5 * _LOG_PI))
        power = float(_exp(log_base + log_scale - 0.5 * _LOG_PI - math.log(time)))
        return HalfSpaceEnergyBudget(magnetic, joule, incoming, power,
                                    math.fsum((incoming, -magnetic, -joule)))

    def adiabatic_temperature_rise(self, depth_m, time_s) -> np.ndarray:
        """Local ΔT [K] with no conduction; reject the singular surface y = 0."""
        depths = _depths(depth_m, positive=True)
        _, _, log_scale = self._scale(time_s)
        if self.boundary_field_t == 0:
            return _immutable(np.zeros_like(depths))
        log_argument = 2 * (np.log(depths) - log_scale) - _LOG_TWO
        log_e1 = _log_exp1_from_log_argument(log_argument)
        log_coefficient = (2 * math.log(abs(self.boundary_field_t))
                           - math.log(self.permeability_h_m) - _LOG_PI
                           - math.log(self.volumetric_heat_capacity_j_m3_k))
        return _immutable(_exp(log_coefficient + log_e1))
