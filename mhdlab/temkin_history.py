"""Versioned coverage-dependent release histories at prescribed temperature.

This record is separate from the constant-hazard surface history. Its moments
cannot justify a constant-hazard emission-time distribution. Each interval uses
a frozen midpoint temperature and the qualified implicit Temkin depletion law;
time subdivision converges only the supplied temperature interpolation.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, fields
import json
import math
from types import MappingProxyType
from typing import Any

import numpy as np

from .flash_surface import FlashSurfaceHistory, validate_surface_attributes
from .surface_sources import matched_uniform_moments, source_contrast
from .thermal_release import SurfaceInventory
from .temkin_release import (
    LAW_IDENTITY, INVERSION_REGIMES, INVERSION_LOG_RESIDUAL_TOLERANCE,
    TemkinReservoir, energy_budgets, temkin_release_step,
)


SCHEMA = "mhdlab-temkin-release-v1"
TIME_MODEL = "piecewise_frozen_midpoint_temperature_implicit_coverage_dependent_depletion"
UNIFORM_MODEL = "spatial_redistribution_preserving_interval_reservoir_number_and_energy_moments"
REGIME_DEFINITIONS = {str(code): label for code, label in INVERSION_REGIMES.items()}
RESIDUAL_INTERPRETATION = "log(I(s)/h_initial) in regime 3; zero elsewhere means no inversion required, not an uncertainty estimate"
REPLAY_RELATIVE_TOLERANCE = max(128 * np.finfo(float).eps, 8 * INVERSION_LOG_RESIDUAL_TOLERANCE)
RELEASE_FIELDS = ("released_m2", "translational_energy_j_m2", "internal_energy_j_m2",
                  "binding_energy_j_m2", "total_wall_energy_debit_j_m2")
STEP_FIELDS = (*RELEASE_FIELDS, "inventory_balance_residual_m2")
DIAGNOSTIC_FIELDS = ("inversion_iterations", "inversion_log_time_residual", "inversion_regime")


def _array(value, name, *, integer=False, maximum=None):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked values")
    if any(isinstance(item, (bool, np.bool_)) for item in np.asarray(value, dtype=object).flat):
        raise ValueError(f"{name} cannot contain booleans")
    raw = np.asarray(value)
    if raw.dtype.kind not in ("iu" if integer else "iuf"):
        raise ValueError(f"{name} must contain {'integer' if integer else 'real numeric'} values without booleans")
    if integer:
        maximum = np.iinfo(np.int64).max if maximum is None else maximum
        if raw.size and (int(np.min(raw)) < 0 or int(np.max(raw)) > maximum):
            raise ValueError(f"{name} lies outside its supported integer domain")
    result = np.asarray(raw, dtype=np.int64 if integer else np.float64)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite values")
    return np.frombuffer(np.ascontiguousarray(result).tobytes(), dtype=result.dtype).reshape(result.shape)


def _freeze(value):
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) or not key for key in value):
            raise ValueError("source metadata keys must be nonempty strings")
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (tuple, list, np.ndarray)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, np.generic):
        return _freeze(value.item())
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float) and np.isfinite(value):
        return value
    raise ValueError("source metadata must contain finite JSON-compatible values")


def _plain(value):
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_plain(item) for item in value]
    return value


def _qualified_metadata(value):
    validate_surface_attributes(value)
    if value["temperature_location"] == "material_sample_plane":
        raise ValueError("Temkin release requires a boundary temperature or explicitly justified boundary proxy")
    if (json.loads(value["current_averaging_interval_m_json"]) is not None
            and not value["current_uniformity_reference"].strip()):
        raise ValueError("Temkin release from a current-averaged temperature requires explicit uniformity qualification")
    return _freeze(value)


@dataclass(frozen=True, slots=True, eq=False)
class TemkinReleaseHistory:
    """Immutable interval moments with shape (interval,reservoir,width).

    Coverage includes both initial and final interval boundaries. Width
    integrals are per unit current length. Inversion diagnostics are numerical
    qualifications of the implicit interval solve, not physical error bars.
    Stored intervals are rechecked against the identified law on construction.
    """

    time_edges_s: np.ndarray
    x_edges_m: np.ndarray
    temperature_k: np.ndarray
    coverage_m2: np.ndarray
    released_m2: np.ndarray
    translational_energy_j_m2: np.ndarray
    internal_energy_j_m2: np.ndarray
    binding_energy_j_m2: np.ndarray
    total_wall_energy_debit_j_m2: np.ndarray
    inventory_balance_residual_m2: np.ndarray
    inversion_iterations: np.ndarray
    inversion_log_time_residual: np.ndarray
    inversion_regime: np.ndarray
    reservoirs: tuple[TemkinReservoir, ...]
    surface_attributes: Mapping[str, Any]
    parameter_provenance: str
    law_identity: str = LAW_IDENTITY
    time_model: str = TIME_MODEL
    schema: str = SCHEMA

    def __post_init__(self):
        if self.law_identity != LAW_IDENTITY or self.time_model != TIME_MODEL or self.schema != SCHEMA:
            raise ValueError("unsupported Temkin history law identity, time model or schema")
        for name in ("time_edges_s", "x_edges_m", "temperature_k", "coverage_m2", *STEP_FIELDS,
                     "inversion_log_time_residual"):
            object.__setattr__(self, name, _array(getattr(self, name), name))
        object.__setattr__(self, "inversion_iterations", _array(self.inversion_iterations, "inversion_iterations", integer=True, maximum=128))
        regime = _array(self.inversion_regime, "inversion_regime", integer=True, maximum=4)
        object.__setattr__(self, "inversion_regime", np.frombuffer(regime.astype(np.int8).tobytes(), dtype=np.int8).reshape(regime.shape))
        if (self.time_edges_s.ndim != 1 or self.time_edges_s.size < 2 or np.any(self.time_edges_s < 0)
                or np.any(self.time_edges_s[1:] <= self.time_edges_s[:-1])):
            raise ValueError("time_edges_s must be nonnegative and strictly increasing")
        if (self.x_edges_m.ndim != 1 or self.x_edges_m.size < 2
                or np.any(self.x_edges_m[1:] <= self.x_edges_m[:-1])
                or not np.isfinite(self.x_edges_m[-1] - self.x_edges_m[0])):
            raise ValueError("x_edges_m must increase over a finite width")
        pools = tuple(self.reservoirs)
        if not pools or any(not isinstance(pool, TemkinReservoir) for pool in pools):
            raise ValueError("reservoirs must contain distinct TemkinReservoir records")
        if len({pool.name for pool in pools}) != len(pools):
            raise ValueError("Temkin reservoir names must be unique")
        object.__setattr__(self, "reservoirs", pools)
        object.__setattr__(self, "surface_attributes", _qualified_metadata(self.surface_attributes))
        if not isinstance(self.parameter_provenance, str) or not self.parameter_provenance.strip():
            raise ValueError("parameter_provenance must be nonempty")
        nt, nr, nx = self.time_edges_s.size - 1, len(pools), self.x_edges_m.size - 1
        if self.temperature_k.shape != (nt, nx) or np.any(self.temperature_k <= 0):
            raise ValueError("temperature_k must be positive interval-by-width values")
        if self.coverage_m2.shape != (nt + 1, nr, nx) or np.any(self.coverage_m2 < 0):
            raise ValueError("coverage_m2 must be nonnegative boundary-by-reservoir-by-width values")
        for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS):
            if getattr(self, name).shape != (nt, nr, nx):
                raise ValueError(f"{name} must have shape (interval,reservoir,width)")
        if any(np.any(getattr(self, name) < 0) for name in RELEASE_FIELDS):
            raise ValueError("released number and energy moments must be nonnegative")
        root = self.inversion_regime == 3
        if (np.any(np.abs(self.inversion_log_time_residual[root]) > INVERSION_LOG_RESIDUAL_TOLERANCE)
                or np.any(self.inversion_log_time_residual[~root] != 0)
                or np.any(self.inversion_iterations[~root] != 0)):
            raise ValueError("Temkin inversion diagnostics exceed their declared numerical qualification")
        if np.any((self.inversion_regime == 0) != (self.coverage_m2[:-1] == 0)):
            raise ValueError("empty-inventory diagnostic regime disagrees with initial coverage")
        if np.any((self.inversion_regime == 4) & (self.coverage_m2[1:] != 0)):
            raise ValueError("survivor-underflow regime requires zero stored survivor")
        for pool_index, pool in enumerate(pools):
            positive = self.coverage_m2[:-1, pool_index] > 0
            regime = self.inversion_regime[:, pool_index]
            constant = pool.coverage_coefficient == 0 or pool.activation_energy_ev == 0
            if np.any(positive & ((regime != 1) if constant else (regime == 1))):
                raise ValueError("analytic diagnostic regime disagrees with the declared coverage dependence")
        tolerance = 128 * np.finfo(float).eps
        if np.any(self.coverage_m2[1:] > self.coverage_m2[:-1]):
            raise ValueError("Temkin desorption cannot increase inventory")
        residual = np.array([
            math.fsum((float(after), -float(before), float(released)))
            for after, before, released in zip(self.coverage_m2[1:].flat,
                                              self.coverage_m2[:-1].flat,
                                              self.released_m2.flat)
        ]).reshape(self.released_m2.shape)
        if (not np.all(np.isfinite(residual)) or np.any(np.abs(residual) > tolerance * self.coverage_m2[:-1])
                or np.any(np.abs(self.inventory_balance_residual_m2) > tolerance * self.coverage_m2[:-1])):
            raise ValueError("Temkin coverage and release do not balance")
        if not np.array_equal(residual, self.inventory_balance_residual_m2):
            raise ValueError("reported inventory residual disagrees with the compensated balance")
        for interval, dt in enumerate(np.diff(self.time_edges_s)):
            trial = temkin_release_step(SurfaceInventory(self.coverage_m2[interval]), pools,
                                        self.temperature_k[interval], float(dt))
            if trial.law_identity != self.law_identity:
                raise ValueError("core trial law disagrees with Temkin history")
            if not np.array_equal(self.inversion_regime[interval], trial.inversion_regime):
                raise ValueError("stored inversion regime disagrees with the qualified Temkin interval")
            # For monotone I(s), dlog(I)/dlog(s)>=1. Release has logarithmic
            # sensitivity <=1 to s; survivor absolute sensitivity is bounded
            # on the initial-inventory scale. This permits equivalent qualified
            # roots from different numerical libraries, not arbitrary rates.
            survivor_difference = np.abs(self.coverage_m2[interval + 1] - trial.trial_inventory.coverage_m2)
            if np.any(survivor_difference > REPLAY_RELATIVE_TOLERANCE * self.coverage_m2[interval]):
                raise ValueError("coverage history disagrees with the identified Temkin interval law")
            if not np.allclose(self.released_m2[interval], trial.released_m2,
                               rtol=REPLAY_RELATIVE_TOLERANCE, atol=0):
                raise ValueError("released_m2 disagrees with the identified Temkin interval law")
            budget = energy_budgets(self.released_m2[interval], pools, self.temperature_k[interval])
            for name in RELEASE_FIELDS[1:]:
                if not np.allclose(getattr(self, name)[interval], budget[name], rtol=tolerance, atol=0):
                    raise ValueError(f"{name} disagrees with the declared Temkin energy convention")

    @property
    def average_flux_m2_s(self):
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            result = self.released_m2 / np.diff(self.time_edges_s)[:, None, None]
        return _array(result, "average_flux_m2_s")

    def width_integral(self, quantity="released_m2"):
        if quantity not in RELEASE_FIELDS:
            raise ValueError("unsupported Temkin interval quantity")
        with np.errstate(over="ignore", invalid="ignore"):
            result = np.sum(getattr(self, quantity) * np.diff(self.x_edges_m), axis=-1)
        return _array(result, "width_integral")

    def matched_uniform(self):
        """Redistribute delivered moments, retaining interval and reservoir."""
        return matched_uniform_moments({name: getattr(self, name) for name in RELEASE_FIELDS}, self.x_edges_m)

    def release_contrast(self, *, edge_fraction=.1, center_fraction=.2):
        return source_contrast(self.released_m2, self.x_edges_m,
                               edge_fraction=edge_fraction, center_fraction=center_fraction)


def evolve_temkin_release(surface: FlashSurfaceHistory, inventory: SurfaceInventory,
                          reservoirs: Sequence[TemkinReservoir], *, parameter_provenance: str,
                          subdivisions_per_interval: int = 1) -> TemkinReleaseHistory:
    """Integrate local Temkin depletion with qualified midpoint temperatures."""
    if not isinstance(surface, FlashSurfaceHistory) or not isinstance(inventory, SurfaceInventory):
        raise ValueError("surface and inventory must be their validated record types")
    if (isinstance(subdivisions_per_interval, (bool, np.bool_))
            or not isinstance(subdivisions_per_interval, (int, np.integer)) or subdivisions_per_interval < 1):
        raise ValueError("subdivisions_per_interval must be a positive integer")
    if not isinstance(parameter_provenance, str) or not parameter_provenance.strip():
        raise ValueError("parameter_provenance must be nonempty")
    surface = FlashSurfaceHistory(**{field.name: getattr(surface, field.name) for field in fields(FlashSurfaceHistory)})
    _qualified_metadata(surface.attributes)
    pools = tuple(reservoirs)
    if not pools or any(not isinstance(pool, TemkinReservoir) for pool in pools):
        raise ValueError("reservoirs must contain TemkinReservoir records")
    count, nx = (surface.time_s.size - 1) * int(subdivisions_per_interval), surface.x_edges_m.size - 1
    if count > 100_000 or count * len(pools) * nx > 20_000_000:
        raise ValueError("Temkin history exceeds the reference memory limit")
    partitions = [np.linspace(lo, hi, int(subdivisions_per_interval) + 1)
                  for lo, hi in zip(surface.time_s[:-1], surface.time_s[1:])]
    edges = np.concatenate((surface.time_s[:1], *(partition[1:] for partition in partitions)))
    if np.any(edges[1:] <= edges[:-1]):
        raise ValueError("time subdivision is not representable at this time origin")
    temperatures = np.concatenate([
        (1 - fraction[:, None]) * surface.temperature_k[index] + fraction[:, None] * surface.temperature_k[index + 1]
        for index, partition in enumerate(partitions)
        for fraction in [(partition[:-1] - partition[0]) / (partition[-1] - partition[0])
                         + .5 * np.diff(partition) / (partition[-1] - partition[0])]
    ])
    current = SurfaceInventory(inventory.coverage_m2)
    coverage = [current.coverage_m2]
    values = {name: [] for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS)}
    for temperature, dt in zip(temperatures, np.diff(edges)):
        trial = temkin_release_step(current, pools, temperature, float(dt))
        if trial.law_identity != LAW_IDENTITY:
            raise ValueError("core trial has an unsupported Temkin law identity")
        current = trial.trial_inventory
        coverage.append(current.coverage_m2)
        for name in values:
            values[name].append(getattr(trial, name))
    return TemkinReleaseHistory(edges, surface.x_edges_m, temperatures, np.stack(coverage),
                                **{name: np.stack(items) for name, items in values.items()},
                                reservoirs=pools, surface_attributes=surface.attributes,
                                parameter_provenance=parameter_provenance)
