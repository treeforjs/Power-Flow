"""Prescribed-temperature surface release upstream of gap and radiation solvers.

The integration is exact for each frozen-temperature first-order interval.
Midpoint sampling approximates a varying temperature history; callers must
refine the time partition. Width bins remain distinct through depletion.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from types import MappingProxyType
from typing import Mapping, Sequence

import numpy as np

from .flash_surface import FlashSurfaceHistory, validate_surface_attributes
from .surface_sources import matched_uniform_moments, source_contrast
from .thermal_release import BindingReservoir, SurfaceInventory, energy_budgets, release_step


_RELEASE_FIELDS = (
    "released_m2", "translational_energy_j_m2", "internal_energy_j_m2",
    "binding_energy_j_m2", "total_wall_energy_debit_j_m2",
)
_STEP_FIELDS = (*_RELEASE_FIELDS, "inventory_balance_residual_m2")


def _immutable(value: np.ndarray) -> np.ndarray:
    def contains_bool(item):
        if isinstance(item, (bool, np.bool_)):
            return True
        if isinstance(item, (list, tuple)):
            return any(contains_bool(part) for part in item)
        if isinstance(item, np.ndarray):
            return item.dtype.kind == "b" or (item.dtype.kind == "O" and any(contains_bool(part) for part in item.flat))
        return False
    if contains_bool(value):
        raise ValueError("surface release arrays cannot contain booleans")
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.dtype.kind == "b":
        raise ValueError("surface release arrays must contain real numeric values")
    result = np.array(raw, dtype=float, copy=True)
    if not np.all(np.isfinite(result)):
        raise ValueError("surface release arrays must remain finite")
    return np.frombuffer(result.tobytes(), dtype=float).reshape(result.shape)


@dataclass(frozen=True)
class SurfaceReleaseHistory:
    """Cell-area quantities for one stationary face, per unit current length.

    Releases/energies have shape (interval,reservoir,width); coverage includes
    the initial state and has shape (interval+1,reservoir,width). ``temperature_k``
    has shape (interval,width). Interval quantities are integrated over time,
    not endpoint rates. Number and energy width integrals have units m^-1 and
    J/m respectively, until multiplied by a declared current-direction length.
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
    reservoirs: tuple[BindingReservoir, ...]
    surface_attributes: Mapping[str, str | float]
    parameter_provenance: str

    def __post_init__(self) -> None:
        for name in ("time_edges_s", "x_edges_m", "temperature_k", "coverage_m2", *_STEP_FIELDS):
            object.__setattr__(self, name, _immutable(getattr(self, name)))
        object.__setattr__(self, "reservoirs", tuple(self.reservoirs))
        object.__setattr__(self, "surface_attributes", MappingProxyType(dict(self.surface_attributes)))
        if (self.time_edges_s.ndim != 1 or len(self.time_edges_s) < 2
                or np.any(self.time_edges_s < 0) or np.any(np.diff(self.time_edges_s) <= 0)):
            raise ValueError("time_edges_s must be nonnegative and strictly increasing")
        if (self.x_edges_m.ndim != 1 or len(self.x_edges_m) < 2
                or np.any(np.diff(self.x_edges_m) <= 0)
                or not np.isfinite(self.x_edges_m[-1] - self.x_edges_m[0])):
            raise ValueError("x_edges_m must be strictly increasing over a finite width")
        nt, nx, nr = len(self.time_edges_s) - 1, len(self.x_edges_m) - 1, len(self.reservoirs)
        if nr == 0 or not all(isinstance(pool, BindingReservoir) for pool in self.reservoirs):
            raise ValueError("reservoirs must contain BindingReservoir records")
        if len({pool.name for pool in self.reservoirs}) != nr:
            raise ValueError("reservoir names must be unique")
        if self.temperature_k.shape != (nt, nx) or np.any(self.temperature_k <= 0):
            raise ValueError("temperature_k must have positive interval-by-width values")
        if self.coverage_m2.shape != (nt + 1, nr, nx) or np.any(self.coverage_m2 < 0):
            raise ValueError("coverage_m2 must have nonnegative boundary-by-reservoir-by-width values")
        for name in _RELEASE_FIELDS:
            value = getattr(self, name)
            if value.shape != (nt, nr, nx) or np.any(value < 0):
                raise ValueError(f"{name} must have nonnegative interval-by-reservoir-by-width values")
        if self.inventory_balance_residual_m2.shape != (nt, nr, nx):
            raise ValueError("inventory_balance_residual_m2 must have interval-by-reservoir-by-width shape")
        if not isinstance(self.parameter_provenance, str) or not self.parameter_provenance.strip():
            raise ValueError("parameter_provenance must be a nonempty source description")
        validate_surface_attributes(self.surface_attributes)
        if self.surface_attributes.get("temperature_location") not in {"surface_boundary", "boundary_temperature_proxy"}:
            raise ValueError("release requires a boundary temperature or an explicitly justified boundary proxy")
        if (self.surface_attributes.get("temperature_location") == "boundary_temperature_proxy"
                and not self.surface_attributes.get("boundary_temperature_reference")):
            raise ValueError("a boundary-temperature proxy must retain its justification")
        if (json.loads(self.surface_attributes.get("current_averaging_interval_m_json", "null")) is not None
                and not self.surface_attributes.get("current_uniformity_reference")):
            raise ValueError("nonlinear release from a current-averaged temperature requires an explicit uniformity reference")
        tolerance = 64 * np.finfo(float).eps
        if np.any(self.coverage_m2[1:] > self.coverage_m2[:-1]):
            raise ValueError("a desorption-only history cannot increase surface coverage")
        residual = self.coverage_m2[1:] + self.released_m2 - self.coverage_m2[:-1]
        if np.any(np.abs(residual) > tolerance * self.coverage_m2[:-1]):
            raise ValueError("surface coverage and released inventory do not balance")
        if np.any(np.abs(self.inventory_balance_residual_m2) > tolerance * self.coverage_m2[:-1]):
            raise ValueError("reported inventory residual exceeds the float64 balance allowance")
        expected = [energy_budgets(amount, self.reservoirs, temperature)
                    for amount, temperature in zip(self.released_m2, self.temperature_k)]
        for name in _RELEASE_FIELDS[1:]:
            target = np.stack([budget[name] for budget in expected])
            if not np.all(np.isfinite(target)) or not np.allclose(getattr(self, name), target, rtol=tolerance, atol=0):
                raise ValueError(f"{name} disagrees with the declared surface energy convention")

    @property
    def average_flux_m2_s(self) -> np.ndarray:
        return self.released_m2 / np.diff(self.time_edges_s)[:, None, None]

    def width_integral(self, quantity: str = "released_m2") -> np.ndarray:
        """Integrate one delivered interval quantity over foil width."""
        if quantity not in _RELEASE_FIELDS:
            raise ValueError(f"quantity must be one of {_RELEASE_FIELDS}")
        result = np.sum(getattr(self, quantity) * np.diff(self.x_edges_m), axis=-1)
        if not np.all(np.isfinite(result)):
            raise ValueError("width-integrated source exceeds the finite numeric range")
        return result

    def matched_uniform(self) -> dict[str, np.ndarray]:
        """Spatially uniform number and energy moments at each time/state.

        This is a controlled redistribution of the width-resolved release.
        It does not recompute desorption at an averaged foil temperature.
        Preserving the full emitted velocity distribution also requires retaining
        the original source-temperature mixture in the downstream sampler.
        """
        return matched_uniform_moments(
            {name: getattr(self, name) for name in _RELEASE_FIELDS}, self.x_edges_m,
        )

    def release_contrast(self, *, edge_fraction: float = 0.1, center_fraction: float = 0.2) -> dict:
        return source_contrast(
            self.released_m2, self.x_edges_m,
            edge_fraction=edge_fraction, center_fraction=center_fraction,
        )


def evolve_surface_release(
    surface: FlashSurfaceHistory,
    inventory: SurfaceInventory,
    reservoirs: Sequence[BindingReservoir],
    *,
    parameter_provenance: str,
    subdivisions_per_interval: int = 1,
) -> SurfaceReleaseHistory:
    """Propose a complete release history without modifying its inputs.

    Linear interpolation follows the validated surface history. Each input
    interval is subdivided uniformly and evaluated at midpoint temperature.
    Reservoir depletion is exact at this frozen temperature. No CR, Boltzmann,
    electromagnetic, surface-reaction or temperature-feedback solve is invoked.
    """
    if not isinstance(surface, FlashSurfaceHistory) or not isinstance(inventory, SurfaceInventory):
        raise ValueError("surface and inventory must use their validated record types")
    if (isinstance(subdivisions_per_interval, (bool, np.bool_))
            or not isinstance(subdivisions_per_interval, (int, np.integer))
            or subdivisions_per_interval < 1):
        raise ValueError("subdivisions_per_interval must be a positive integer")
    if not isinstance(parameter_provenance, str) or not parameter_provenance.strip():
        raise ValueError("parameter_provenance must be a nonempty source description")
    # Revalidate the mutable NumPy storage behind the externally supplied record.
    surface = FlashSurfaceHistory(
        surface.time_s, surface.x_edges_m, surface.temperature_k,
        surface.source_kind, surface.run_id, surface.source_reference,
        surface.surface_region, surface.flash_version,
        surface.temperature_location, surface.material_sample_depth_m, surface.boundary_temperature_reference,
        surface.current_averaging_interval_m, surface.current_uniformity_reference,
    )
    if surface.temperature_location == "material_sample_plane":
        raise ValueError("a material sample plane is not a boundary temperature; an explicit justified boundary proxy is required")
    if surface.current_averaging_interval_m is not None and not surface.current_uniformity_reference:
        raise ValueError("a current-averaged temperature requires an explicit uniformity reference or separately resolved slices")
    pools = tuple(reservoirs)
    nintervals = (len(surface.time_s) - 1) * int(subdivisions_per_interval)
    nx = len(surface.x_edges_m) - 1
    if nintervals > 100_000 or nintervals * max(len(pools), 1) * nx > 20_000_000:
        raise ValueError("release history exceeds the configured reference memory limit")
    partitions = [np.linspace(start, stop, int(subdivisions_per_interval) + 1)
                  for start, stop in zip(surface.time_s[:-1], surface.time_s[1:])]
    time_edges = np.concatenate((surface.time_s[:1], *(partition[1:] for partition in partitions)))
    if np.any(np.diff(time_edges) <= 0):
        raise ValueError("requested time subdivision is not representable at this time origin")
    # Evaluate at the mathematical midpoint in each original interpolation
    # interval. Constructing its absolute time can round onto an endpoint when
    # two stored times are adjacent representable floats.
    temperatures = np.concatenate([
        (1 - fraction[:, None]) * surface.temperature_k[index]
        + fraction[:, None] * surface.temperature_k[index + 1]
        for index, partition in enumerate(partitions)
        for fraction in [
            (partition[:-1] - partition[0]) / (partition[-1] - partition[0])
            + 0.5 * (np.diff(partition) / (partition[-1] - partition[0]))
        ]
    ])
    current = SurfaceInventory(inventory.coverage_m2)
    coverage = [current.coverage_m2]
    quantities: dict[str, list[np.ndarray]] = {name: [] for name in _STEP_FIELDS}
    for temperature, dt in zip(temperatures, np.diff(time_edges)):
        update = release_step(current, pools, temperature, float(dt))
        current = update.trial_inventory
        coverage.append(current.coverage_m2)
        for name in _STEP_FIELDS:
            quantities[name].append(getattr(update, name))
    return SurfaceReleaseHistory(
        time_edges, surface.x_edges_m, temperatures, np.stack(coverage),
        **{name: np.stack(values) for name, values in quantities.items()},
        reservoirs=pools, surface_attributes=surface.attributes,
        parameter_provenance=parameter_provenance,
    )
