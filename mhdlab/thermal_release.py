"""Finite-inventory thermal release driven by prescribed surface temperatures.

Every reservoir and physical parameter is supplied by the caller. There are no
water properties, site energies, species populations or accommodation defaults.
Each step freezes its given temperature and evaluates the exact first-order law
N_after = N_before * exp[-nu*dt*exp(-Ea/(kB*T))]. No random sampling is performed.

Energy convention: ``desorption_energy_ev`` is the nonnegative energy cost from
the adsorbed reservoir to a specified gas reference with zero internal excitation.
``internal_energy_ev`` is the nonnegative internal energy of ``state_id`` above
that gas reference. Supply each contribution once; do not insert a whole chemical
formation energy again if it is already in the desorption cost. The activation
barrier controls the rate and is NOT an additional wall-energy debit.

The outgoing translational budget is 2*kB*T per released particle, conditional
on a thermally accommodating, flux-weighted Maxwellian emission boundary. This
module reports its mean energy; it does not sample velocities. The energy debit
is a diagnostic for prescribed T, not a coupled wall-temperature update.

Coverage and all outputs use float64. Very small release is retained even when
it is below one spacing of the remaining inventory. The explicit balance
residual records that representational limit; there is no physical rate floor.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Real
from typing import Sequence

import numpy as np

from .constants import KB, QE


_KB_EV_K = KB / QE
_LOG_FLOAT_MAX = math.log(np.finfo(float).max)
_LOG_FLOAT_EPS = math.log(np.finfo(float).eps)
_LOG_FLOAT_TINY = math.log(np.finfo(float).tiny)


def _label(value: str, name: str) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise ValueError(f"{name} must be a nonempty string without surrounding whitespace")
    return value


def _real(value: Real, name: str, *, positive: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number, not bool or complex")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result) or result < 0 or (positive and result == 0):
        qualifier = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must be finite and {qualifier}")
    return result


def _array(value, name: str, ndim: int) -> np.ndarray:
    # The object representation preserves bool/string types in mixed lists;
    # coercing to float first would silently turn True into 1 and '2' into 2.
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} must not contain masked values")
    try:
        raw = np.asarray(value, dtype=object)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a rectangular finite real array") from exc
    if raw.ndim != ndim or any(length == 0 for length in raw.shape):
        raise ValueError(f"{name} must be a nonempty {ndim}-dimensional array")
    for item in raw.flat:
        _real(item, name)
    return np.asarray(raw, dtype=float)


def _immutable(array: np.ndarray) -> np.ndarray:
    # An immutable bytes backing also prevents writeable=True from bypassing
    # the array flag, unlike a writable owned allocation with a read-only view.
    return np.frombuffer(np.ascontiguousarray(array, dtype=float).tobytes(), dtype=float).reshape(array.shape)


def _product(*factors: float) -> float:
    """Multiply nonnegative finite factors without intermediate underflow."""
    if any(factor == 0 for factor in factors):
        return 0.0
    mantissa, exponent = 1.0, 0
    for factor in factors:
        fraction, power = math.frexp(factor)
        mantissa *= fraction
        exponent += power
    try:
        result = math.ldexp(mantissa, exponent)
    except OverflowError as exc:
        raise ValueError("thermal release result exceeds finite float64 range") from exc
    if not math.isfinite(result):
        raise ValueError("thermal release result exceeds finite float64 range")
    return result


def _ratio(numerator: float, denominator: float) -> float:
    if numerator == 0:
        return 0.0
    first, first_power = math.frexp(numerator)
    second, second_power = math.frexp(denominator)
    try:
        result = math.ldexp(first / second, first_power - second_power)
    except OverflowError as exc:
        raise ValueError("average release flux exceeds finite float64 range") from exc
    if not math.isfinite(result):
        raise ValueError("average release flux exceeds finite float64 range")
    return result


@dataclass(frozen=True, slots=True)
class BindingReservoir:
    name: str
    state_id: str
    attempt_frequency_s: float
    activation_energy_ev: float
    desorption_energy_ev: float
    internal_energy_ev: float = 0.0

    def __post_init__(self) -> None:
        _label(self.name, "reservoir name")
        _label(self.state_id, "state_id")
        for name in ("attempt_frequency_s", "activation_energy_ev", "desorption_energy_ev", "internal_energy_ev"):
            object.__setattr__(self, name, _real(getattr(self, name), name, positive=name == "attempt_frequency_s"))


@dataclass(frozen=True, slots=True, eq=False)
class SurfaceInventory:
    """Immutable particle coverage [m^-2], indexed reservoir then width bin.

    Rows must be kept in the same order as the reservoirs passed to release_step.
    Distinct binding pools can emit the same gas state; their names remain unique.
    """

    coverage_m2: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "coverage_m2", _immutable(_array(self.coverage_m2, "coverage_m2", 2)))


@dataclass(frozen=True, slots=True, eq=False)
class ThermalReleaseTrial:
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


def _pools(reservoirs: Sequence[BindingReservoir]) -> tuple[BindingReservoir, ...]:
    if isinstance(reservoirs, (str, bytes)):
        raise ValueError("reservoirs must be a sequence of BindingReservoir records")
    try:
        pools = tuple(reservoirs)
    except TypeError as exc:
        raise ValueError("reservoirs must be a sequence of BindingReservoir records") from exc
    if not pools or any(not isinstance(pool, BindingReservoir) for pool in pools):
        raise ValueError("reservoirs must be a nonempty sequence of BindingReservoir records")
    if len({pool.name for pool in pools}) != len(pools):
        raise ValueError("reservoir names must be unique")
    return pools


def _energy_arrays(released: np.ndarray, pools: tuple[BindingReservoir, ...],
                   temperatures: np.ndarray) -> dict[str, np.ndarray]:
    fields = {name: np.empty_like(released) for name in (
        "translational_energy_j_m2", "internal_energy_j_m2", "binding_energy_j_m2",
        "total_wall_energy_debit_j_m2",
    )}
    for pool_index, pool in enumerate(pools):
        for width_index, temperature in enumerate(temperatures):
            location = (pool_index, width_index)
            amount = float(released[location])
            translation = _product(2.0, KB, float(temperature), amount)
            internal = _product(pool.internal_energy_ev, QE, amount)
            binding = _product(pool.desorption_energy_ev, QE, amount)
            fields["translational_energy_j_m2"][location] = translation
            fields["internal_energy_j_m2"][location] = internal
            fields["binding_energy_j_m2"][location] = binding
            try:
                fields["total_wall_energy_debit_j_m2"][location] = math.fsum((translation, internal, binding))
            except OverflowError as exc:
                raise ValueError("total wall energy debit exceeds finite float64 range") from exc
    if any(not np.all(np.isfinite(value)) for value in fields.values()):
        raise ValueError("thermal release energy budget exceeds finite float64 range")
    return {name: _immutable(value) for name, value in fields.items()}


def energy_budgets(
    released_m2: Sequence[Sequence[float]] | np.ndarray,
    reservoirs: Sequence[BindingReservoir],
    temperature_k: Sequence[float] | np.ndarray,
) -> dict[str, np.ndarray]:
    """Return the four immutable per-reservoir/bin energy arrays for a release.

    The shape of ``released_m2`` is (reservoir,width), with one strictly positive
    temperature per width bin. This evaluates the same scaled products used by
    release_step, including cases where a partial product would under/overflow
    although the final energy is representable. It neither depletes inventory
    nor infers a rate. Activation energy is excluded from these budgets.
    """
    pools = _pools(reservoirs)
    released = _array(released_m2, "released_m2", 2)
    temperatures = _array(temperature_k, "temperature_k", 1)
    if released.shape != (len(pools), len(temperatures)):
        raise ValueError("released_m2 shape must match reservoir count and temperature width bins")
    if np.any(temperatures <= 0):
        raise ValueError("temperature_k must be strictly positive")
    return _energy_arrays(released, pools, temperatures)


def _log_hazard(reservoir: BindingReservoir, temperature: float, dt: float) -> float:
    if reservoir.activation_energy_ev == 0:
        barrier = 0.0
    else:
        log_barrier = math.fsum((math.log(reservoir.activation_energy_ev), -math.log(_KB_EV_K),
                                -math.log(temperature)))
        barrier = math.inf if log_barrier > _LOG_FLOAT_MAX else math.exp(log_barrier)
    return math.fsum((math.log(reservoir.attempt_frequency_s), math.log(dt), -barrier))


def _deplete(coverage: float, log_hazard: float) -> tuple[float, float]:
    if coverage == 0:
        return 0.0, 0.0
    log_coverage = math.log(coverage)
    if log_hazard < _LOG_FLOAT_EPS:
        # For h < epsilon, -expm1(-h) = h to float64 relative precision.
        # Combine logs before exponentiation, so N*h survives even when h by
        # itself is subnormal/unrepresentable. This is not a rate cutoff.
        released = math.exp(math.fsum((log_coverage, log_hazard)))
        return coverage - released, released
    if log_hazard > _LOG_FLOAT_MAX:
        return 0.0, coverage
    hazard = math.exp(log_hazard)
    released = coverage * (-math.expm1(-hazard))
    # Preserve representable N*exp(-h) when exp(-h) alone would underflow.
    if -hazard < _LOG_FLOAT_TINY:
        remaining = math.exp(log_coverage - hazard)
    else:
        remaining = coverage * math.exp(-hazard)
    return remaining, released


def release_step(
    inventory: SurfaceInventory,
    reservoirs: Sequence[BindingReservoir],
    temperature_k: Sequence[float] | np.ndarray,
    dt_s: float,
) -> ThermalReleaseTrial:
    """Return a nonmutating frozen-T depletion trial and particle/energy budgets.

    Every budget array has shape (reservoir, width bin). Total wall energy is
    binding + specified internal + mean thermal translation; the activation
    barrier never appears as a separate debit. Values are per unit face area;
    the caller must supply face areas to form integrated particle/energy totals.
    Temperature and dt must be strictly positive. Rates have no floor.
    """
    if not isinstance(inventory, SurfaceInventory):
        raise ValueError("inventory must be a SurfaceInventory")
    pools = _pools(reservoirs)
    if len(pools) != inventory.coverage_m2.shape[0]:
        raise ValueError("reservoir count must match inventory rows")
    temperatures = _array(temperature_k, "temperature_k", 1)
    if temperatures.shape != (inventory.coverage_m2.shape[1],):
        raise ValueError("temperature_k must have one value for each width bin")
    if np.any(temperatures <= 0):
        raise ValueError("temperature_k must be strictly positive")
    dt = _real(dt_s, "dt_s", positive=True)
    remaining, released, flux, balance = (
        np.empty_like(inventory.coverage_m2) for _ in range(4)
    )
    for pool_index, pool in enumerate(pools):
        for width_index, temperature in enumerate(temperatures):
            location = (pool_index, width_index)
            original = float(inventory.coverage_m2[location])
            after, amount = _deplete(original, _log_hazard(pool, float(temperature), dt))
            remaining[location], released[location] = after, amount
            flux[location] = _ratio(amount, dt)
            balance[location] = math.fsum((after, amount, -original))
    if not all(np.all(np.isfinite(array)) for array in (remaining, released, flux, balance)):
        raise ValueError("thermal release result exceeds finite float64 range")
    return ThermalReleaseTrial(
        trial_inventory=SurfaceInventory(remaining),
        reservoir_names=tuple(pool.name for pool in pools),
        state_ids=tuple(pool.state_id for pool in pools),
        temperature_k=_immutable(temperatures), dt_s=dt,
        released_m2=_immutable(released), average_flux_m2_s=_immutable(flux),
        **_energy_arrays(released, pools, temperatures),
        inventory_balance_residual_m2=_immutable(balance),
    )
