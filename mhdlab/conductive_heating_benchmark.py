"""Exact periodic resistive induction and conductive heating, in SI units.

A stationary conductor has constant electrical conductivity ``sigma``, magnetic
permeability ``mu``, volumetric heat capacity ``Cv`` and thermal conductivity
``kappa``.  Its initial tangential field is ``B0*cos(k*y)`` and its initial
temperature rise is zero.  There is no boundary energy flux in the periodic
domain, fluid motion, Hall term, radiation, phase change or equation of state.
This is a verification reference, not a FLASH run or a production heating code.

Let eta=1/(mu*sigma), lambda=2*eta*k**2, chi=kappa/Cv, delta=4*chi*k**2.
Then Bx=B0*exp(-lambda*t/2)*cos(k*y), Jz=-dBx/dy/mu, and
q=Jz**2/sigma.  The conductive temperature solution is

    dT = B0**2/(4*mu*Cv) * [1-exp(-x)-H(x,z)*cos(2*k*y)],
    H(x,z) = x*(exp(-x)-exp(-z))/(z-x), x=lambda*t, z=delta*t.

At resonance H(x,x)=x*exp(-x).  At kappa=0, H(x,0)=1-exp(-x),
giving the local adiabatic limit.  Thermal conduction changes the profile but
not the spatially averaged thermal energy.  ``Cv`` is supplied per volume;
no conversion from mass-specific constant-pressure heat capacity is implied.

All coefficients are explicit.  Signed B0 and nonnegative kappa are allowed.
Results are float64; an unrepresentable final result raises, while a final
value below the subnormal range rounds to zero.  Raw rational products and
logarithmic scaling retain representable results when intermediate rates,
dimensionless times, amplitudes or small geometric factors underflow/overflow.
Phases satisfy |k*y| <= MAX_ABSOLUTE_PHASE_RAD: larger coordinates are rejected
rather than claiming accurate range reduction of arbitrary large phases.
Positive quadrature of a bounded smooth kernel evaluates the small-time
temperature at field maxima without subtracting two almost equal terms.
Finite-volume averages integrate q independently of the average current.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Real

import numpy as np


MAX_ABSOLUTE_PHASE_RAD = float(2**20)
ASSUMPTIONS = (
    "periodic one-dimensional stationary conductor",
    "constant electrical conductivity, permeability, volumetric heat capacity and thermal conductivity",
    "initial Bx=B0*cos(k*y) and zero temperature rise",
    "resistive induction and Fourier thermal conduction only; no Hall, flow, radiation, phase change or EOS",
)
_LOG_MAX = math.log(np.finfo(float).max)
_LOG_TWO = math.log(2.0)
_LOG_FOUR = math.log(4.0)
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(16)
_UNIT_NODES, _UNIT_WEIGHTS = (_NODES + 1.0) / 2.0, _WEIGHTS / 2.0


class UnsupportedConductiveBenchmarkRangeError(ValueError):
    """A requested phase or finite result is outside the qualified float range."""


def _real(value, name: str, *, positive: bool = False, nonnegative: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number, not bool, string or complex")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _coordinates(values, name: str, *, edges: bool = False) -> np.ndarray:
    if np.ma.isMaskedArray(values):
        raise ValueError(f"{name} cannot be masked")
    try:
        raw = np.asarray(values, dtype=object)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real scalar or one-dimensional array") from exc
    if raw.ndim > 1 or not raw.size or (edges and (raw.ndim != 1 or raw.size < 2)):
        raise ValueError(f"{name} must be {'a one-dimensional array with at least two edges' if edges else 'a nonempty scalar or one-dimensional array'}")
    result = np.array([_real(item, name) for item in raw.flat], dtype=float).reshape(raw.shape)
    if edges and any(right <= left for left, right in zip(result[:-1], result[1:])):
        raise ValueError(f"{name} must be strictly increasing")
    return result


def _immutable(value) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    return np.frombuffer(np.ascontiguousarray(array).tobytes(), dtype=float).reshape(array.shape)


def _log_fraction(value: Fraction) -> float:
    return math.log(value.numerator) - math.log(value.denominator) if value else -math.inf


def _fraction_float(value: Fraction) -> float:
    try:
        return float(value)
    except OverflowError:
        return math.inf


def _exp(log_value: float, sign: float = 1.0) -> float:
    if log_value == -math.inf or sign == 0:
        return 0.0
    if not math.isfinite(log_value) or log_value > _LOG_MAX:
        raise UnsupportedConductiveBenchmarkRangeError("requested benchmark result exceeds finite float64 range")
    return sign * math.exp(log_value)


def _logadd(first: float, second: float) -> float:
    high, low = max(first, second), min(first, second)
    return high + math.log1p(math.exp(low - high)) if high != -math.inf else -math.inf


def _log_depletion(value: Fraction) -> float:
    """log(1-exp(-value)), retaining a sub-float positive value."""
    if not value:
        return -math.inf
    log_value = _log_fraction(value)
    if log_value < -36.0:
        return log_value  # omitted relative correction is < 1.2e-16
    if value >= 80:
        return 0.0  # complement < 1.9e-35, below rounding of unity
    return math.log(-math.expm1(-float(value)))


def _log_exprel_negative(value: Fraction) -> float:
    """log((1-exp(-value))/value), with value >= 0."""
    return _log_depletion(value) - _log_fraction(value) if value else 0.0


def _thermal_shapes(x: Fraction, z: Fraction) -> tuple[float, float, float]:
    """Return log F, log H and log D, where F=1-e^-x and D=F-H>=0.

    D is the CDF of a sum of two independent exponential waiting times,
    hence symmetric in x,z.  Small rates use the positive triangle integral
    D=x*z*integral_0^1 (1-u)*exp(-x*u)*exprel(-z*(1-u)) du.
    For x,z<=1 the 16-point Gauss rule has a smooth bounded kernel; comparison
    with independent adaptive quadrature is part of the verification suite.
    """
    log_f = _log_depletion(x)
    if not x:
        return log_f, -math.inf, -math.inf
    if not z:
        return log_f, log_f, -math.inf
    smaller, larger = min(x, z), max(x, z)
    difference = larger - smaller
    log_h = _log_fraction(x) - _fraction_float(smaller) + _log_exprel_negative(difference)
    log_small, log_large = _log_fraction(smaller), _log_fraction(larger)
    if larger <= 1:
        a, b = float(smaller), float(larger)
        coefficient = math.fsum(
            float(weight) * (1.0 - float(u)) * math.exp(-a * float(u))
            * (1.0 if b * (1.0 - float(u)) == 0 else -math.expm1(-b * (1.0 - float(u))) / (b * (1.0 - float(u))))
            for u, weight in zip(_UNIT_NODES, _UNIT_WEIGHTS)
        )
        log_d = log_small + log_large + math.log(coefficient)
    elif smaller < Fraction(1, 2):
        fa = math.exp(_log_exprel_negative(smaller))
        fb = math.exp(_log_exprel_negative(larger))
        # Here fa-fb >= exprel(-.5)-exprel(-1) > .154, so no
        # cancellation of a tiny restored temperature is possible.
        ratio = float(smaller / larger)
        log_d = log_small + math.log(fa - fb) - math.log1p(-ratio)
    elif smaller >= 80:
        # 1-D <= (1+smaller)*exp(-smaller) < 1.5e-33.
        log_d = 0.0
    else:
        a = float(smaller)
        f_difference = math.exp(_log_exprel_negative(difference))
        log_survival = -a + math.log1p(a * f_difference)
        log_d = math.log(-math.expm1(log_survival))
    return log_f, log_h, log_d


def _trig(phase: Fraction) -> tuple[float, float, float]:
    """cos(phase), sign(sin(phase)), log(abs(sin(phase)))."""
    if abs(phase) > MAX_ABSOLUTE_PHASE_RAD:
        raise UnsupportedConductiveBenchmarkRangeError(f"|k*y| must not exceed {MAX_ABSOLUTE_PHASE_RAD:g} radians")
    if not phase:
        return 1.0, 0.0, -math.inf
    value = float(phase)
    if abs(phase) < Fraction(1, 100000):
        correction = math.sin(value) / value if value else 1.0
        return math.cos(value), math.copysign(1.0, phase), _log_fraction(abs(phase)) + math.log(correction)
    sine = math.sin(value)
    return math.cos(value), math.copysign(1.0, sine), math.log(abs(sine)) if sine else -math.inf


def _sinc(half_width: Fraction) -> tuple[float, float]:
    if not half_width:
        return 1.0, 0.0
    _, sign, log_sin = _trig(half_width)
    return sign, log_sin - _log_fraction(half_width)


def _log_cell_heating_shape(center: Fraction, half_width: Fraction) -> float:
    """log(1-cos(2*center)*sinc(2*half_width)), including tiny bins."""
    cos_center, _, log_sin_center = _trig(center)
    twice_width = 2 * half_width
    if twice_width <= 1:
        r = float(twice_width)
        # 1-sinc(r) = r²/6*(1-r²/20+r⁴/840-...), evaluated
        # through the r^12 term.  For r<=1 the next term is < 3e-13
        # of this loss; use the direct expression for r>.1 instead.
        if twice_width < Fraction(1, 10):
            r2 = r * r
            correction = 1.0 + r2 * (-1.0 / 20.0 + r2 * (1.0 / 840.0 + r2 * (-1.0 / 60480.0 + r2 * (1.0 / 6652800.0 + r2 * (-1.0 / 1037836800.0)))))
            log_loss = 2 * _log_fraction(twice_width) - math.log(6.0) + math.log(correction)
            log_sinc = math.log1p(-math.exp(log_loss))
        else:
            sinc = math.sin(r) / r
            log_loss, log_sinc = math.log1p(-sinc), math.log(sinc)
        return _logadd(log_loss, log_sinc + _LOG_TWO + 2 * log_sin_center)
    value = float(twice_width)
    sinc = math.sin(value) / value
    cos_twice = 2 * cos_center * cos_center - 1.0
    return math.log(1.0 - cos_twice * sinc)


def _signed_component(log_scale: float, value: float) -> float:
    return _exp(log_scale + math.log(abs(value)), math.copysign(1.0, value)) if value else 0.0


@dataclass(frozen=True, slots=True, eq=False)
class PeriodicConductiveFields:
    position_m: np.ndarray
    magnetic_field_t: np.ndarray
    current_density_a_m2: np.ndarray
    joule_power_w_m3: np.ndarray
    temperature_rise_k: np.ndarray


@dataclass(frozen=True, slots=True, eq=False)
class PeriodicConductiveCellAverages:
    position_edges_m: np.ndarray
    magnetic_field_t: np.ndarray
    current_density_a_m2: np.ndarray
    joule_power_w_m3: np.ndarray
    temperature_rise_k: np.ndarray


@dataclass(frozen=True, slots=True)
class PeriodicConductiveEnergyBudget:
    initial_mean_magnetic_energy_j_m3: float
    mean_magnetic_energy_j_m3: float
    mean_thermal_energy_j_m3: float
    mean_total_energy_j_m3: float
    mean_joule_power_w_m3: float
    energy_residual_j_m3: float


@dataclass(frozen=True, slots=True)
class PeriodicConductiveHeating:
    """Exact single-mode periodic benchmark with caller-supplied SI constants.

    Methods accept finite time >= 0.  Points are scalar or one-dimensional;
    edges are a strictly increasing one-dimensional array.  Arrays returned
    by the methods are immutable.  A zero thermal conductivity explicitly
    selects the local adiabatic limit.  Derived scalar properties raise when
    their positive mathematical value cannot be represented as float64;
    field methods still combine raw factors before rounding final results.
    """

    conductivity_s_m: float
    permeability_h_m: float
    volumetric_heat_capacity_j_m3_k: float
    thermal_conductivity_w_m_k: float
    wavenumber_rad_m: float
    initial_field_t: float

    def __post_init__(self) -> None:
        for name in ("conductivity_s_m", "permeability_h_m", "volumetric_heat_capacity_j_m3_k", "wavenumber_rad_m"):
            object.__setattr__(self, name, _real(getattr(self, name), name, positive=True))
        object.__setattr__(self, "thermal_conductivity_w_m_k", _real(self.thermal_conductivity_w_m_k, "thermal_conductivity_w_m_k", nonnegative=True))
        object.__setattr__(self, "initial_field_t", _real(self.initial_field_t, "initial_field_t"))

    def _constants(self) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction, Fraction]:
        return tuple(Fraction(value) for value in (self.conductivity_s_m, self.permeability_h_m, self.volumetric_heat_capacity_j_m3_k, self.thermal_conductivity_w_m_k, self.wavenumber_rad_m, self.initial_field_t))

    @staticmethod
    def _positive_property(value: Fraction) -> float:
        result = _fraction_float(value)
        if not math.isfinite(result) or (value and result == 0):
            raise UnsupportedConductiveBenchmarkRangeError("derived coefficient is outside positive float64 range")
        return result

    @property
    def magnetic_diffusivity_m2_s(self) -> float:
        sigma, mu, *_ = self._constants()
        return self._positive_property(1 / (mu * sigma))

    @property
    def thermal_diffusivity_m2_s(self) -> float:
        _, _, cv, kappa, *_ = self._constants()
        return self._positive_property(kappa / cv)

    @property
    def resistive_energy_decay_rate_s(self) -> float:
        sigma, mu, _, _, k, _ = self._constants()
        return self._positive_property(2 * k * k / (mu * sigma))

    @property
    def thermal_mode_decay_rate_s(self) -> float:
        _, _, cv, kappa, k, _ = self._constants()
        return self._positive_property(4 * kappa * k * k / cv)

    @property
    def period_m(self) -> float:
        return self._positive_property(Fraction(2 * math.pi) / Fraction(self.wavenumber_rad_m))

    def _scales(self, time_s) -> tuple[Fraction, Fraction, float, float, float, float, float]:
        time = Fraction(_real(time_s, "time_s", nonnegative=True))
        sigma, mu, cv, kappa, k, b = self._constants()
        x, z = 2 * k * k * time / (mu * sigma), 4 * kappa * k * k * time / cv
        log_b = _log_fraction(abs(b))
        log_j = _log_fraction(abs(b) * k / mu)
        log_q = _log_fraction(b * b * k * k / (mu * mu * sigma))
        log_t = _log_fraction(b * b / (4 * mu * cv))
        return x, z, log_b, log_j, log_q, log_t, math.copysign(1.0, self.initial_field_t)

    def fields(self, position_m, time_s) -> PeriodicConductiveFields:
        positions = _coordinates(position_m, "position_m")
        x, z, log_b, log_j, log_q, log_t, sign_b = self._scales(time_s)
        attenuation = _fraction_float(x)
        _, log_h, log_d = _thermal_shapes(x, z)
        arrays = [np.empty_like(positions) for _ in range(4)]
        k = Fraction(self.wavenumber_rad_m)
        for index in np.ndindex(positions.shape):
            cosine, sign_sine, log_sine = _trig(k * Fraction(float(positions[index])))
            arrays[0][index] = _signed_component(log_b - attenuation / 2, sign_b * cosine)
            arrays[1][index] = _exp(log_j - attenuation / 2 + log_sine, sign_b * sign_sine)
            arrays[2][index] = _exp(log_q - attenuation + 2 * log_sine)
            arrays[3][index] = _exp(log_t + _logadd(log_d, log_h + _LOG_TWO + 2 * log_sine))
        return PeriodicConductiveFields(_immutable(positions), *(_immutable(array) for array in arrays))

    def cell_averages(self, position_edges_m, time_s) -> PeriodicConductiveCellAverages:
        edges = _coordinates(position_edges_m, "position_edges_m", edges=True)
        x, z, log_b, log_j, log_q, log_t, sign_b = self._scales(time_s)
        attenuation = _fraction_float(x)
        _, log_h, log_d = _thermal_shapes(x, z)
        arrays = [np.empty(edges.size - 1) for _ in range(4)]
        k = Fraction(self.wavenumber_rad_m)
        for index, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
            left_phase, right_phase = k * Fraction(float(left)), k * Fraction(float(right))
            _trig(left_phase)
            _trig(right_phase)
            center, half_width = (left_phase + right_phase) / 2, (right_phase - left_phase) / 2
            cosine, sign_sine, log_sine = _trig(center)
            sign_sinc, log_sinc = _sinc(half_width)
            arrays[0][index] = _signed_component(log_b - attenuation / 2 + log_sinc, sign_b * sign_sinc * cosine)
            arrays[1][index] = _exp(log_j - attenuation / 2 + log_sinc + log_sine, sign_b * sign_sinc * sign_sine)
            log_shape = _log_cell_heating_shape(center, half_width)
            arrays[2][index] = _exp(log_q - attenuation + log_shape - _LOG_TWO)
            arrays[3][index] = _exp(log_t + _logadd(log_d, log_h + log_shape))
        return PeriodicConductiveCellAverages(_immutable(edges), *(_immutable(array) for array in arrays))

    def energy_budget(self, time_s) -> PeriodicConductiveEnergyBudget:
        time = Fraction(_real(time_s, "time_s", nonnegative=True))
        sigma, mu, _, _, k, b = self._constants()
        rate = 2 * k * k / (mu * sigma)
        x = rate * time
        log_initial = _log_fraction(b * b / (4 * mu))
        attenuation = _fraction_float(x)
        initial = _exp(log_initial)
        magnetic = _exp(log_initial - attenuation)
        thermal = _exp(log_initial + _log_depletion(x))
        total = math.fsum((magnetic, thermal))
        power = _exp(log_initial + _log_fraction(rate) - attenuation)
        residual = math.fsum((magnetic, thermal, -initial))
        return PeriodicConductiveEnergyBudget(initial, magnetic, thermal, total, power, residual)
