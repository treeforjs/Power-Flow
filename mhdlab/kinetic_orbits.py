"""CPU reference Boris orbits in prescribed electric and magnetic fields.

Positions and mechanical momenta are three dimensional, including the component
outside a two-dimensional spatial mesh. The caller samples E and B at the desired
step midpoint and owns field interpolation, subcycling, boundaries, and collisions.
There is no field solve, charge deposition, space-charge closure, or Hall feedback
here. A large Boris rotation is bounded but does not resolve the physical orbit;
the returned gyro angles let the caller control temporal accuracy without clipping.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


SPEED_OF_LIGHT_M_S = 299_792_458.0


def _numeric_array(value: Any, name: str) -> np.ndarray:
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must contain real numeric values")
    try:
        result = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must contain finite numeric values") from exc
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _vectors(value: Any, name: str, count: int | None = None) -> np.ndarray:
    result = _numeric_array(value, name)
    if result.ndim != 2 or result.shape[1] != 3 or (count is not None and result.shape[0] != count):
        raise ValueError(f"{name} must have shape (n, 3), with the same n as the momenta")
    return result


def _scalars(value: Any, name: str, count: int, *, positive: bool = False) -> np.ndarray:
    result = _numeric_array(value, name)
    if result.shape != (count,):
        raise ValueError(f"{name} must have shape (n,); scalar broadcasting is not supported")
    if positive and np.any(result <= 0.0):
        raise ValueError(f"{name} must be strictly positive")
    return result


def _validate_relativistic(relativistic: bool) -> None:
    if not isinstance(relativistic, (bool, np.bool_)):
        raise ValueError("relativistic must be a boolean")


def _gamma(momentum: np.ndarray, mass: np.ndarray, relativistic: bool) -> np.ndarray:
    if not relativistic:
        return np.ones(mass.shape, dtype=np.float64)
    norm = np.hypot.reduce(momentum, axis=1)
    return np.hypot(1.0, norm / (mass * SPEED_OF_LIGHT_M_S))


def _energy(momentum: np.ndarray, mass: np.ndarray, relativistic: bool) -> np.ndarray:
    norm = np.hypot.reduce(momentum, axis=1)
    if relativistic:
        # Rationalized (gamma-1) m c^2 avoids cancellation for cold particles.
        return norm * (norm / (mass * (_gamma(momentum, mass, True) + 1.0)))
    return 0.5 * norm * (norm / mass)


def _finite_result(*arrays: np.ndarray) -> None:
    if any(not np.all(np.isfinite(array)) for array in arrays):
        raise ValueError("orbit calculation exceeded the finite numeric range")


def kinetic_energy_j(
    momentum_kg_m_s: np.ndarray,
    mass_kg: np.ndarray,
    *,
    relativistic: bool = True,
) -> np.ndarray:
    """Per-particle kinetic energy; input shapes are (n,3) and (n,)."""
    _validate_relativistic(relativistic)
    momentum = _vectors(momentum_kg_m_s, "momentum_kg_m_s")
    mass = _scalars(mass_kg, "mass_kg", momentum.shape[0], positive=True)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            energy = _energy(momentum, mass, relativistic)
    except FloatingPointError as exc:
        raise ValueError("kinetic energy exceeded the finite numeric range") from exc
    _finite_result(energy)
    return energy


def velocity_m_s(
    momentum_kg_m_s: np.ndarray,
    mass_kg: np.ndarray,
    *,
    relativistic: bool = True,
) -> np.ndarray:
    """Mechanical momentum to velocity; nonrelativistic mode uses p/m."""
    _validate_relativistic(relativistic)
    momentum = _vectors(momentum_kg_m_s, "momentum_kg_m_s")
    mass = _scalars(mass_kg, "mass_kg", momentum.shape[0], positive=True)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            velocity = momentum / (mass * _gamma(momentum, mass, relativistic))[:, None]
    except FloatingPointError as exc:
        raise ValueError("velocity exceeded the finite numeric range") from exc
    _finite_result(velocity)
    return velocity


@dataclass(frozen=True)
class OrbitStep:
    """One prescribed-field step and per-particle accuracy diagnostics.

    kinetic_delta_j is the particle energy change. electric_displacement_work_j
    is q E dot (x_new-x_old), using the supplied E and trapezoidal displacement.
    Their difference is reported, since they need not agree at finite timestep
    in relativistic mode. gyro_angle_rad is |q| |B| dt/(m gamma_minus), whereas
    boris_rotation_angle_rad is the actual bounded rotation 2 atan(gyro_angle/2).
    magnetic_rotation_relative_error measures norm preservation during the B
    rotation alone; it is zero for an exactly zero momentum before rotation.
    """

    position_m: np.ndarray
    momentum_kg_m_s: np.ndarray
    kinetic_delta_j: np.ndarray
    electric_displacement_work_j: np.ndarray
    work_residual_j: np.ndarray
    gyro_angle_rad: np.ndarray
    boris_rotation_angle_rad: np.ndarray
    magnetic_rotation_relative_error: np.ndarray


def advance_orbits(
    position_m: np.ndarray,
    momentum_kg_m_s: np.ndarray,
    mass_kg: np.ndarray,
    charge_c: np.ndarray,
    dt_s: float,
    electric_field_v_m: np.ndarray,
    magnetic_field_t: np.ndarray,
    *,
    relativistic: bool = True,
) -> OrbitStep:
    """Advance a batch by an electric half-kick, Boris rotation, half-kick.

    Phase and field arrays must be (n,3); masses and charges must be (n,), with
    finite, positive masses and a finite, positive scalar timestep. Charges may
    have either sign or be zero. Zero charge is ballistic. Relativistic rotation
    uses gamma from momentum after the first electric half-kick (standard
    relativistic Boris). Positions use the mean of old and new velocities.

    Inputs are never mutated. E/B are prescribed samples, not callbacks, and
    must be sampled by the caller for each substep. No orbit-angle clipping or
    automatic accuracy acceptance is performed. Nonrelativistic mode requires
    the caller to remain in its low-speed physical regime.
    """
    _validate_relativistic(relativistic)
    momentum = _vectors(momentum_kg_m_s, "momentum_kg_m_s")
    count = momentum.shape[0]
    position = _vectors(position_m, "position_m", count)
    mass = _scalars(mass_kg, "mass_kg", count, positive=True)
    charge = _scalars(charge_c, "charge_c", count)
    electric = _vectors(electric_field_v_m, "electric_field_v_m", count)
    magnetic = _vectors(magnetic_field_t, "magnetic_field_t", count)
    dt = _numeric_array(dt_s, "dt_s")
    if dt.ndim != 0 or dt <= 0.0:
        raise ValueError("dt_s must be a finite, positive scalar")
    dt = float(dt)

    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            half_kick = (0.5 * dt * charge)[:, None] * electric
            momentum_minus = momentum + half_kick
            gamma_minus = _gamma(momentum_minus, mass, relativistic)
            t = (0.5 * dt * charge / mass / gamma_minus)[:, None] * magnetic
            t_norm = np.hypot.reduce(t, axis=1)
            inverse_norm = 1.0 / np.hypot(1.0, t_norm)
            scaled_t = t * inverse_norm[:, None]
            # Algebraically the usual p' = p- + p- x t; p+ = p- + p' x s,
            # s = 2t/(1+t.t). Scaling t first avoids squaring a large gyro angle.
            cross_once = np.cross(momentum_minus, scaled_t)
            momentum_rotated = momentum_minus + 2.0 * (
                inverse_norm[:, None] * cross_once + np.cross(cross_once, scaled_t)
            )
            momentum_new = momentum_rotated + half_kick
            velocity_old = momentum / (mass * _gamma(momentum, mass, relativistic))[:, None]
            velocity_new = momentum_new / (mass * _gamma(momentum_new, mass, relativistic))[:, None]
            displacement = 0.5 * dt * (velocity_old + velocity_new)
            position_new = position + displacement
            kinetic_delta = _energy(momentum_new, mass, relativistic) - _energy(momentum, mass, relativistic)
            electric_work = charge * np.sum(electric * displacement, axis=1)
            residual = kinetic_delta - electric_work
            before_rotation = np.hypot.reduce(momentum_minus, axis=1)
            after_rotation = np.hypot.reduce(momentum_rotated, axis=1)
            rotation_error = np.divide(
                np.abs(after_rotation - before_rotation), before_rotation,
                out=np.zeros(count), where=before_rotation > 0.0,
            )
            gyro_angle = 2.0 * t_norm
            rotation_angle = 2.0 * np.arctan(t_norm)
    except FloatingPointError as exc:
        raise ValueError("orbit calculation exceeded the finite numeric range; check inputs and timestep") from exc
    _finite_result(position_new, momentum_new, kinetic_delta, electric_work, residual, gyro_angle, rotation_angle, rotation_error)
    return OrbitStep(
        position_new, momentum_new, kinetic_delta, electric_work, residual,
        gyro_angle, rotation_angle, rotation_error,
    )
