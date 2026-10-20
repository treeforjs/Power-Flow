"""Coverage-dependent first-order release at prescribed, frozen temperature.

This separate law is

    dN/dt = -nu*N*exp[-Ed*(1-alpha*N/Nref)/(kB*T)].

``nu`` is a supplied constant attempt frequency, not an Eyring prefactor.
``Nref`` is an independently supplied areal-density convention; it is never
inferred from the initial inventory.  This remaining-population prefactor is
not the dissertation's printed initial-population prefactor.  No water, foil,
monolayer, binding-energy or accommodation parameters are supplied here.

The caller declares a supported coverage interval [0,Nmax].  The kinetic
barrier must be nonnegative throughout it, with no clipping.  Its coverage
dependence does NOT define a surface internal-energy function: the separately
supplied, constant ``desorption_energy_ev`` and gas ``internal_energy_ev`` are
used in the shared energy ledger, alongside mean outgoing translation 2*kB*T.
Prescribed temperature means that wall debit remains a diagnostic.

For positive coupling, solve the monotone exact separated equation in log
depletion s=log(Ninitial/Nfinal):

    integral_0^s exp[c*(1-exp(-u))] du = h_initial,
    c = Ed/(kB*T)*alpha*Ninitial/Nref,
    h_initial = nu*dt*exp[-Ed/(kB*T)*(1-alpha*Ninitial/Nref)].

Endpoint-scaled quadrature and logarithmic bracketing avoid E1 subtraction and
large-barrier cancellation.  This reference implementation qualifies c<=512;
larger nonlinear coupling raises UnsupportedTemkinRangeError, rather than
substituting a constant initial rate.  The limit is numerical, not a physical
cutoff.  Coverage ratios and their deficiencies are formed as exact ratios of
the supplied float64 values before thermal scaling.  Tiny-depletion and
unrepresentable-survivor limits use bounds on the same integral.  Float64
inventory residuals are retained explicitly, and no rate floor is imposed.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Real
from types import MappingProxyType
from typing import Sequence
import warnings

import numpy as np
from scipy import integrate, optimize

from .constants import KB, QE
from .thermal_release import (
    BindingReservoir, SurfaceInventory, _array, _deplete, _immutable, _label,
    _log_hazard, _ratio, _real, energy_budgets as _constant_energy_budgets,
)


LAW_IDENTITY = "temkin_first_order_constant_attempt_frequency_v1"
MAX_DIMENSIONLESS_COVERAGE_EFFECT = 512.0
INVERSION_LOG_RESIDUAL_TOLERANCE = 2e-11
INVERSION_REGIMES = MappingProxyType({
    0: "empty_inventory",
    1: "constant_barrier_analytic",
    2: "tiny_depletion_bounded",
    3: "nonlinear_integral_inversion",
    4: "unrepresentable_survivor_bounded",
})
_LOG_EPS = math.log(np.finfo(float).eps)
_MAX_FLOAT_FRACTION = Fraction(float(np.finfo(float).max))
_LOG_MIN_SUBNORMAL = math.log(float(np.nextafter(0.0, 1.0)))
_LOG_MIN_NORMAL = math.log(float(np.finfo(float).tiny))


class UnsupportedTemkinRangeError(ValueError):
    """The supplied law lies outside the qualified numerical inversion range."""


class TemkinInversionError(RuntimeError):
    """The nonlinear integral did not meet its numerical verification criteria."""


def _positive_real(value, name: str) -> float:
    return _real(value, name, positive=True)


def _nonnegative_real(value, name: str) -> float:
    return _real(value, name)


def _fraction(value: float) -> Fraction:
    return Fraction.from_float(value)


def _immutable_integer(value, dtype) -> np.ndarray:
    array = np.ascontiguousarray(value, dtype=dtype)
    return np.frombuffer(array.tobytes(), dtype=dtype).reshape(array.shape)


@dataclass(frozen=True, slots=True)
class TemkinReservoir:
    """Explicit kinetic/energy hypotheses for one supported binding pool.

    All fields are required.  ``activation_energy_ev`` is the zero-coverage
    kinetic barrier Ed; ``coverage_coefficient`` is dimensionless alpha.
    ``max_coverage_m2`` declares the supported domain, not its initial value.
    At Ed=0 the entire barrier is identically zero for any finite alpha.
    """

    name: str
    state_id: str
    attempt_frequency_s: float
    activation_energy_ev: float
    desorption_energy_ev: float
    internal_energy_ev: float
    coverage_coefficient: float
    reference_coverage_m2: float
    max_coverage_m2: float
    kinetic_provenance: str
    energy_provenance: str

    def __post_init__(self) -> None:
        for name in ("name", "state_id", "kinetic_provenance", "energy_provenance"):
            _label(getattr(self, name), name)
        for name in ("attempt_frequency_s", "reference_coverage_m2", "max_coverage_m2"):
            object.__setattr__(self, name, _positive_real(getattr(self, name), name))
        for name in ("activation_energy_ev", "desorption_energy_ev", "internal_energy_ev", "coverage_coefficient"):
            object.__setattr__(self, name, _nonnegative_real(getattr(self, name), name))
        if self.activation_energy_ev > 0:
            domain_ratio = (_fraction(self.coverage_coefficient) * _fraction(self.max_coverage_m2)
                            / _fraction(self.reference_coverage_m2))
            if domain_ratio > 1:
                raise ValueError("kinetic barrier must be nonnegative throughout [0,max_coverage_m2]; "
                                 "alpha*max_coverage_m2/reference_coverage_m2 exceeds one")

    def _energy_record(self) -> BindingReservoir:
        """Explicit energy-only adapter; it does not transfer this kinetic law."""
        return BindingReservoir(self.name, self.state_id, self.attempt_frequency_s,
                                self.activation_energy_ev, self.desorption_energy_ev,
                                self.internal_energy_ev)


@dataclass(frozen=True, slots=True, eq=False)
class TemkinReleaseTrial:
    trial_inventory: SurfaceInventory
    reservoir_names: tuple[str, ...]
    state_ids: tuple[str, ...]
    temperature_k: np.ndarray
    dt_s: float
    released_m2: np.ndarray
    average_flux_m2_s: np.ndarray
    translational_energy_j_m2: np.ndarray
    internal_energy_j_m2: np.ndarray
    binding_energy_j_m2: np.ndarray
    total_wall_energy_debit_j_m2: np.ndarray
    inventory_balance_residual_m2: np.ndarray
    inversion_iterations: np.ndarray
    inversion_log_time_residual: np.ndarray
    inversion_regime: np.ndarray
    law_identity: str = LAW_IDENTITY


def _pools(reservoirs: Sequence[TemkinReservoir]) -> tuple[TemkinReservoir, ...]:
    if isinstance(reservoirs, (str, bytes)):
        raise ValueError("reservoirs must be a nonempty sequence of TemkinReservoir records")
    try:
        pools = tuple(reservoirs)
    except TypeError as exc:
        raise ValueError("reservoirs must be a nonempty sequence of TemkinReservoir records") from exc
    if not pools or any(not isinstance(pool, TemkinReservoir) for pool in pools):
        raise ValueError("reservoirs must be a nonempty sequence of TemkinReservoir records")
    if len({pool.name for pool in pools}) != len(pools):
        raise ValueError("reservoir names must be unique")
    return pools


def energy_budgets(released_m2, reservoirs: Sequence[TemkinReservoir], temperature_k) -> dict[str, np.ndarray]:
    """The shared constant-cost ledger, with an explicit reservoir conversion.

    This energy-only operation neither evaluates a rate nor infers a binding
    cost from the coverage-dependent activation barrier.  A release amount
    cannot exceed the declared supported inventory of its reservoir.
    """
    pools = _pools(reservoirs)
    if np.ma.isMaskedArray(released_m2) or np.ma.isMaskedArray(temperature_k):
        raise ValueError("release amounts and temperature_k cannot contain masked values")
    released = _array(released_m2, "released_m2", 2)
    if released.shape[0] != len(pools):
        raise ValueError("released_m2 rows must match reservoir count")
    for row, pool in enumerate(pools):
        if np.any(released[row] > pool.max_coverage_m2):
            raise ValueError("released_m2 exceeds the reservoir's max_coverage_m2")
    return _constant_energy_budgets(released, tuple(pool._energy_record() for pool in pools), temperature_k)


def _thermal_parameters(pool: TemkinReservoir, coverage: float, temperature: float) -> tuple[float, float]:
    # Keep the complement exact: beta - beta*f could discard a small positive
    # barrier when beta is enormous and f rounds to one.
    fraction = (_fraction(pool.coverage_coefficient) * _fraction(coverage)
                / _fraction(pool.reference_coverage_m2))
    beta = (_fraction(pool.activation_energy_ev) * _fraction(QE)
            / (_fraction(KB) * _fraction(temperature)))
    coupling_fraction = beta * fraction
    if coupling_fraction > _fraction(MAX_DIMENSIONLESS_COVERAGE_EFFECT):
        raise UnsupportedTemkinRangeError(
            f"dimensionless coverage effect Ed/(kBT)*alpha*N/Nref exceeds "
            f"the qualified limit {MAX_DIMENSIONLESS_COVERAGE_EFFECT:g}")
    coupling = float(coupling_fraction)
    initial_barrier_fraction = beta * (1 - fraction)
    initial_barrier = (math.inf if initial_barrier_fraction > _MAX_FLOAT_FRACTION
                       else float(initial_barrier_fraction))
    return coupling, initial_barrier


def _log_integral(log_depletion: float, coupling: float) -> float:
    depletion = math.exp(log_depletion)
    endpoint_exponent = coupling * (-math.expm1(-depletion))

    def scaled_integrand(location: float) -> float:
        exponent = (coupling * math.exp(-depletion * location)
                    * math.expm1(-depletion * (1 - location)))
        return math.exp(exponent)

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", integrate.IntegrationWarning)
            value, error = integrate.quad(scaled_integrand, 0., 1., epsabs=0., epsrel=2e-13, limit=200)
    except (integrate.IntegrationWarning, ValueError, OverflowError) as exc:
        raise TemkinInversionError("Temkin integral quadrature failed its accuracy criteria") from exc
    if not math.isfinite(value) or value <= 0 or error > 2e-11 * value:
        raise TemkinInversionError("Temkin integral quadrature failed its accuracy criteria")
    return math.fsum((log_depletion, endpoint_exponent, math.log(value)))


def _nonlinear_deplete(coverage: float, coupling: float,
                       log_initial_hazard: float) -> tuple[float, float, int, float, int]:
    # I(s)>=s and I(s)<=s*exp(c).  If c*h_initial is below epsilon,
    # the nonlinear correction is below float64 relative precision.
    if log_initial_hazard + math.log(max(1., coupling)) < _LOG_EPS:
        after, released = _deplete(coverage, log_initial_hazard)
        return after, released, 0, 0., 2
    log_coverage = math.log(coverage)
    # This exceeds the half-subnormal rounding threshold by >1 in s.
    zero_survivor_depletion = math.fsum((log_coverage, -_LOG_MIN_SUBNORMAL, 2.))
    log_zero_survivor = math.log(zero_survivor_depletion)
    if log_initial_hazard >= _log_integral(log_zero_survivor, coupling):
        return 0., coverage, 0, 0., 4
    lower = log_initial_hazard - coupling
    upper = min(log_initial_hazard, log_zero_survivor)

    def residual(value: float) -> float:
        return _log_integral(value, coupling) - log_initial_hazard

    low_value, high_value = residual(lower), residual(upper)
    # Expand only a representationally marginal bracket.  This is a numerical
    # bracket correction, not clipping of the kinetic barrier or inventory.
    margin = 32 * np.finfo(float).eps * max(1., abs(lower), abs(upper))
    if low_value > 0:
        lower -= margin
        low_value = residual(lower)
    if high_value < 0:
        upper += margin
        high_value = residual(upper)
    if low_value > 0 or high_value < 0:
        raise TemkinInversionError("Temkin integral did not bracket the monotone solution")
    if abs(high_value) <= 4 * np.finfo(float).eps:
        solution, iterations = upper, 0
    elif abs(low_value) <= 4 * np.finfo(float).eps:
        solution, iterations = lower, 0
    else:
        try:
            solution, result = optimize.brentq(residual, lower, upper, xtol=1e-14,
                                               rtol=4 * np.finfo(float).eps,
                                               maxiter=128, full_output=True, disp=False)
        except (ValueError, RuntimeError) as exc:
            raise TemkinInversionError("Temkin integral inversion failed") from exc
        if not result.converged:
            raise TemkinInversionError("Temkin integral inversion did not converge")
        iterations = result.iterations
    final_residual = residual(solution)
    if not math.isfinite(final_residual) or abs(final_residual) > INVERSION_LOG_RESIDUAL_TOLERANCE:
        raise TemkinInversionError("Temkin inverse time residual exceeded its qualified tolerance")
    depletion = math.exp(solution)
    released = coverage * (-math.expm1(-depletion))
    remaining = (math.exp(log_coverage - depletion) if -depletion < _LOG_MIN_NORMAL
                 else coverage * math.exp(-depletion))
    return remaining, released, iterations, final_residual, 3


def release_step(inventory: SurfaceInventory, reservoirs: Sequence[TemkinReservoir],
                 temperature_k, dt_s: Real) -> TemkinReleaseTrial:
    """Return a nonmutating frozen-T nonlinear release trial.

    Every diagnostic array is indexed (reservoir,width).  ``inversion_regime``
    uses the exported INVERSION_REGIMES codes.  ``inversion_log_time_residual``
    is log(I(s)/h_initial) for regime 3; zero in other regimes means that no
    nonlinear inversion was required, rather than an estimated physical error.
    Iteration counts likewise concern only the scalar nonlinear root solve.
    """
    if not isinstance(inventory, SurfaceInventory):
        raise ValueError("inventory must be a SurfaceInventory")
    pools = _pools(reservoirs)
    if inventory.coverage_m2.shape[0] != len(pools):
        raise ValueError("reservoir count must match inventory rows")
    if np.ma.isMaskedArray(temperature_k):
        raise ValueError("temperature_k cannot contain masked values")
    temperatures = _array(temperature_k, "temperature_k", 1)
    if temperatures.shape != (inventory.coverage_m2.shape[1],) or np.any(temperatures <= 0):
        raise ValueError("temperature_k must have one strictly positive value per width bin")
    dt = _positive_real(dt_s, "dt_s")
    for row, pool in enumerate(pools):
        if np.any(inventory.coverage_m2[row] > pool.max_coverage_m2):
            raise ValueError("coverage_m2 exceeds the reservoir's max_coverage_m2")
    remaining, released, flux, balance, residuals = (
        np.empty_like(inventory.coverage_m2) for _ in range(5))
    iterations = np.zeros(inventory.coverage_m2.shape, dtype=np.int64)
    regimes = np.zeros(inventory.coverage_m2.shape, dtype=np.int8)
    for row, pool in enumerate(pools):
        constant = pool._energy_record()
        for width, temperature in enumerate(temperatures):
            location = (row, width)
            original = float(inventory.coverage_m2[location])
            count, residual, regime = 0, 0., 0
            if original == 0:
                after, amount = 0., 0.
            elif pool.coverage_coefficient == 0 or pool.activation_energy_ev == 0:
                after, amount = _deplete(original, _log_hazard(constant, float(temperature), dt))
                regime = 1
            else:
                coupling, initial_barrier = _thermal_parameters(pool, original, float(temperature))
                log_hazard = math.fsum((math.log(pool.attempt_frequency_s), math.log(dt), -initial_barrier))
                after, amount, count, residual, regime = _nonlinear_deplete(original, coupling, log_hazard)
            remaining[location], released[location] = after, amount
            flux[location] = _ratio(amount, dt)
            balance[location] = math.fsum((after, -original, amount))
            residuals[location], iterations[location], regimes[location] = residual, count, regime
    if not all(np.all(np.isfinite(array)) for array in (remaining, released, flux, balance, residuals)):
        raise TemkinInversionError("Temkin trial exceeds the finite float64 diagnostic range")
    return TemkinReleaseTrial(
        trial_inventory=SurfaceInventory(remaining),
        reservoir_names=tuple(pool.name for pool in pools), state_ids=tuple(pool.state_id for pool in pools),
        temperature_k=_immutable(temperatures), dt_s=dt,
        released_m2=_immutable(released), average_flux_m2_s=_immutable(flux),
        **energy_budgets(released, pools, temperatures), inventory_balance_residual_m2=_immutable(balance),
        inversion_iterations=_immutable_integer(iterations, np.int64),
        inversion_log_time_residual=_immutable(residuals),
        inversion_regime=_immutable_integer(regimes, np.int8),
    )


temkin_release_step = release_step
