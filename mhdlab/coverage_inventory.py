"""Explicit physical parent inventory and fixed-kinetics Temkin source sweeps.

The layer coordinate scales one parent's inventory. Pool/state fractions divide
that inventory; they are not dissociation yields, competing shared sites, or a
replacement for separately owned chemistry. Spatial multipliers are literal and
are never normalized. The kinetic reference coverage remains supplied and fixed.
"""

from __future__ import annotations

from dataclasses import KW_ONLY, asdict, dataclass
import math
from numbers import Integral, Real

import numpy as np

from .planar_release import PlanarSurfaceHistory, _plain
from .temkin_planar import (
    TemkinPlanarReleaseHistory, _memory_limit, _planar_attributes,
    evolve_temkin_planar_release,
)
from .temkin_release import (
    MAX_DIMENSIONLESS_COVERAGE_EFFECT, TemkinReservoir, _thermal_parameters,
)


SCHEMA = "mhdlab-coverage-source-inventory-v1"
SWEEP_SCHEMA = "mhdlab-coverage-source-sweep-v1"
ALLOCATION_MODEL = "single_parent_initial_population_partition_no_dissociation_yields_or_shared_site_competition"
PROFILE_MODEL = "literal_local_inventory_multiplier_without_normalization"
STATUSES = frozenset({"manufactured", "reconstructed", "sourced", "unvalidated", "owner_supplied", "observed", "source_evaluated"})
FRACTION_SUM_TOLERANCE = 64*float(np.finfo(float).eps)
MAX_COVERAGE_KNOTS = 1000
MAX_SWEEP_INTERVAL_RESERVOIR_PATCHES = 20_000_000
LD = np.longdouble


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be an explicit nonempty string")
    return value


def _status(value):
    if not isinstance(value, str) or value not in STATUSES:
        raise ValueError(f"status must be one of {sorted(STATUSES)}")


def _real(value, name, *, positive=False):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar, not boolean or masked")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be a finite representable real scalar") from error
    if not math.isfinite(result) or result < 0 or (positive and result == 0):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}")
    return result


def _immutable(value):
    array = np.asarray(value, dtype=float)
    return np.frombuffer(np.ascontiguousarray(array).tobytes(), dtype=float).reshape(array.shape)


def _array(value, name, ndim=None, *, nonnegative=True):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked values")
    try:
        objects = np.asarray(value, dtype=object)
    except (ValueError, TypeError) as error:
        raise ValueError(f"{name} must be a finite real array") from error
    if objects.size == 0 or (ndim is not None and objects.ndim != ndim):
        raise ValueError(f"{name} has unsupported or empty shape")
    for item in objects.flat:
        if nonnegative:
            _real(item, name)
        else:
            if isinstance(item, (bool, np.bool_)) or not isinstance(item, Real) or not math.isfinite(float(item)):
                raise ValueError(f"{name} must contain finite real numbers")
    return _immutable(np.asarray(objects, dtype=float))


def _map(value, name):
    result = _array(value, name)
    if result.ndim not in (0, 2):
        raise ValueError(f"{name} must be a scalar or a (current,width) patch array")
    return result


def _broadcast(value, shape, name):
    if value.ndim == 0:
        return _immutable(np.full(shape, float(value)))
    if value.shape != shape:
        raise ValueError(f"{name} must match the surface (current,width) shape")
    return value


def _represented(value, name):
    value = np.asarray(value, dtype=LD)
    if np.any(~np.isfinite(value)) or np.any(np.abs(value) > np.finfo(float).max):
        raise ValueError(f"{name} exceeds finite float64 range")
    with np.errstate(under="ignore"):
        rounded = np.asarray(value, dtype=float)
    if np.any((value != 0) & (rounded == 0)):
        raise ValueError(f"{name} is nonzero but below representable float64 range")
    return _immutable(rounded)


def _number(density, areas, name):
    return _represented(np.sum(np.asarray(density, dtype=LD)*np.asarray(areas, dtype=LD), axis=(-2, -1), dtype=LD), name)


def _time_partitions(surface, subdivisions):
    partitions = [np.linspace(lo, hi, subdivisions+1) for lo, hi in zip(surface.time_s[:-1], surface.time_s[1:])]
    edges = np.concatenate((surface.time_s[:1], *(partition[1:] for partition in partitions)))
    if np.any(edges[1:] <= edges[:-1]):
        raise ValueError("source time subdivision is not representable at the supplied time origin")
    return partitions, edges


@dataclass(frozen=True, slots=True)
class MonolayerDefinition:
    areal_density_m2: float
    definition: str
    provenance: str
    status: str

    def __post_init__(self):
        object.__setattr__(self, "areal_density_m2", _real(self.areal_density_m2, "monolayer areal_density_m2", positive=True))
        _text(self.definition, "monolayer definition")
        _text(self.provenance, "monolayer provenance")
        _status(self.status)


@dataclass(frozen=True, slots=True, eq=False)
class InventoryPoolAllocation:
    reservoir: TemkinReservoir
    fraction: np.ndarray
    provenance: str
    status: str

    def __post_init__(self):
        if not isinstance(self.reservoir, TemkinReservoir):
            raise ValueError("allocation reservoir must be a supplied TemkinReservoir")
        fraction = _map(self.fraction, "allocation fraction")
        if np.any(fraction > 1):
            raise ValueError("each allocation fraction must lie in [0,1]")
        object.__setattr__(self, "fraction", fraction)
        _text(self.provenance, "allocation provenance")
        _status(self.status)


@dataclass(frozen=True, slots=True, eq=False)
class ParentSurfaceInventory:
    parent_id: str
    monolayer: MonolayerDefinition
    allocations: tuple[InventoryPoolAllocation, ...]
    patch_multiplier: np.ndarray
    _: KW_ONLY
    profile_provenance: str
    status: str
    provenance: str

    def __post_init__(self):
        _text(self.parent_id, "parent_id")
        if not isinstance(self.monolayer, MonolayerDefinition):
            raise ValueError("monolayer must be an explicit MonolayerDefinition")
        try:
            allocations = tuple(self.allocations)
        except TypeError as error:
            raise ValueError("allocations must contain supplied InventoryPoolAllocation records") from error
        if not allocations or any(not isinstance(item, InventoryPoolAllocation) for item in allocations):
            raise ValueError("allocations must contain supplied InventoryPoolAllocation records")
        if len({item.reservoir.name for item in allocations}) != len(allocations):
            raise ValueError("allocation reservoir names must be unique; one pool cannot duplicate a parent")
        object.__setattr__(self, "allocations", allocations)
        object.__setattr__(self, "patch_multiplier", _map(self.patch_multiplier, "patch_multiplier"))
        _text(self.profile_provenance, "profile_provenance")
        _text(self.provenance, "parent provenance")
        _status(self.status)


@dataclass(frozen=True, slots=True, eq=False)
class CoverageInventoryPlan:
    surface: PlanarSurfaceHistory
    parent: ParentSurfaceInventory
    coverage_ml: np.ndarray
    subdivisions_per_interval: int
    patch_multiplier: np.ndarray
    allocation_fractions: np.ndarray
    initial_parent_coverage_m2: np.ndarray
    initial_pool_coverage_m2: np.ndarray
    parent_number: np.ndarray
    pool_number: np.ndarray
    allocation_balance_residual_m2: np.ndarray
    max_dimensionless_coverage_effect: np.ndarray

    def __post_init__(self):
        if not isinstance(self.surface, PlanarSurfaceHistory) or not isinstance(self.parent, ParentSurfaceInventory):
            raise ValueError("inventory plans require supplied validated surface and parent records")
        if isinstance(self.subdivisions_per_interval, (bool, np.bool_)) or not isinstance(self.subdivisions_per_interval, Integral) or self.subdivisions_per_interval < 1:
            raise ValueError("plan subdivisions_per_interval must be a positive integer")
        for name, ndim in (("coverage_ml", 1), ("patch_multiplier", 2), ("allocation_fractions", 3),
                           ("initial_parent_coverage_m2", 3), ("initial_pool_coverage_m2", 4),
                           ("parent_number", 1), ("pool_number", 2),
                           ("allocation_balance_residual_m2", 3), ("max_dimensionless_coverage_effect", 1)):
            object.__setattr__(self, name, _array(getattr(self, name), name, ndim,
                                                nonnegative=name != "allocation_balance_residual_m2"))
        nc, nr = len(self.coverage_ml), len(self.parent.allocations)
        shape = (len(self.surface.z_edges_m)-1, len(self.surface.x_edges_m)-1)
        if len(self.coverage_ml) > MAX_COVERAGE_KNOTS or np.any(self.coverage_ml[1:] <= self.coverage_ml[:-1]):
            raise ValueError("plan coverage knots must be finite, strictly increasing and within the knot limit")
        expected_shapes = {"patch_multiplier": shape, "allocation_fractions": (nr, *shape),
                           "initial_parent_coverage_m2": (nc, *shape), "initial_pool_coverage_m2": (nc, nr, *shape),
                           "parent_number": (nc,), "pool_number": (nc, nr),
                           "allocation_balance_residual_m2": (nc, *shape), "max_dimensionless_coverage_effect": (nr,)}
        if any(getattr(self, name).shape != expected for name, expected in expected_shapes.items()):
            raise ValueError("inventory plan arrays must match its coverage, pool and physical patch geometry")
        expected_multiplier = _broadcast(self.parent.patch_multiplier, shape, "patch_multiplier")
        expected_fractions = np.stack([_broadcast(item.fraction, shape, "allocation fraction") for item in self.parent.allocations])
        if not np.array_equal(self.patch_multiplier, expected_multiplier) or not np.array_equal(self.allocation_fractions, expected_fractions):
            raise ValueError("plan spatial multiplier or allocations disagree with the supplied parent definition")
        if np.any(np.abs(np.sum(np.asarray(self.allocation_fractions, dtype=LD), axis=0, dtype=LD)-1) > LD(FRACTION_SUM_TOLERANCE)):
            raise ValueError("plan pool fractions must partition one physical parent at every patch")
        expected_parent = _represented(np.asarray(self.coverage_ml, dtype=LD)[:, None, None]*LD(self.parent.monolayer.areal_density_m2)*np.asarray(self.patch_multiplier, dtype=LD)[None], "plan parent density")
        expected_pools = _represented(np.asarray(expected_parent, dtype=LD)[:, None]*np.asarray(self.allocation_fractions, dtype=LD)[None], "plan pool density")
        expected_residual = _represented(np.asarray(expected_parent, dtype=LD)-np.sum(np.asarray(expected_pools, dtype=LD), axis=1, dtype=LD), "plan allocation residual")
        areas = self.surface.patch_areas_m2
        _represented(np.sum(np.asarray(areas, dtype=LD)*np.asarray(self.patch_multiplier, dtype=LD), dtype=LD)/np.sum(areas, dtype=LD), "area-weighted multiplier mean")
        _represented(np.asarray(expected_parent, dtype=LD)*np.asarray(areas, dtype=LD), "initial patch parent number")
        _represented(np.asarray(expected_pools, dtype=LD)*np.asarray(areas, dtype=LD), "initial patch pool number")
        for actual, expected in ((self.initial_parent_coverage_m2, expected_parent), (self.initial_pool_coverage_m2, expected_pools),
                                 (self.allocation_balance_residual_m2, expected_residual),
                                 (self.parent_number, _number(expected_parent, areas, "plan parent number")),
                                 (self.pool_number, _number(expected_pools, areas, "plan pool number"))):
            if not np.array_equal(actual, expected):
                raise ValueError("inventory plan density, allocation residual or actual-area number ledger is inconsistent")
        nt = (len(self.surface.time_s)-1)*int(self.subdivisions_per_interval)
        _memory_limit(nt, nr, *shape)
        if nc*nt*nr*math.prod(shape) > MAX_SWEEP_INTERVAL_RESERVOIR_PATCHES:
            raise ValueError("inventory plan exceeds the aggregate source sweep memory limit")
        _planar_attributes(self.surface.attributes, self.surface.x_edges_m, self.surface.z_edges_m)
        _time_partitions(self.surface, int(self.subdivisions_per_interval))
        minimum_t = np.min(self.surface.temperature_k, axis=0)
        effects = []
        for row, allocation in enumerate(self.parent.allocations):
            maximum = np.max(self.initial_pool_coverage_m2[:, row], axis=0)
            if np.any(maximum > allocation.reservoir.max_coverage_m2):
                raise ValueError("plan initial pool density exceeds its supplied maximum")
            effects.append(max(_thermal_parameters(allocation.reservoir, float(value), float(temperature))[0]
                               for value, temperature in zip(maximum.flat, minimum_t.flat)))
        if not np.array_equal(self.max_dimensionless_coverage_effect, effects):
            raise ValueError("plan coupling qualification disagrees with its supplied pools and local temperature bounds")

    def manifest(self):
        areas = self.surface.patch_areas_m2
        multiplier_mean = np.sum(np.asarray(areas, dtype=LD)*np.asarray(self.patch_multiplier, dtype=LD), dtype=LD)/np.sum(areas, dtype=LD)
        return {
            "schema": SCHEMA, "parent_id": self.parent.parent_id,
            "parent_status": self.parent.status, "parent_provenance": self.parent.provenance,
            "monolayer": {"areal_density_m2": self.parent.monolayer.areal_density_m2, "definition": self.parent.monolayer.definition, "provenance": self.parent.monolayer.provenance, "status": self.parent.monolayer.status},
            "coverage_ml": self.coverage_ml.tolist(), "coverage_coordinate": "nominal layer amplitude before literal local multiplier",
            "profile_model": PROFILE_MODEL, "profile_provenance": self.parent.profile_provenance,
            "patch_multiplier": self.patch_multiplier.tolist(), "area_weighted_multiplier_mean": _represented(multiplier_mean, "area-weighted multiplier mean").item(),
            "area_mean_coverage_ml": _represented(np.asarray(self.coverage_ml, dtype=LD)*multiplier_mean, "area-mean coverage").tolist(),
            "allocation_model": ALLOCATION_MODEL, "allocation_fractions": self.allocation_fractions.tolist(),
            "allocation_fraction_sum_tolerance": FRACTION_SUM_TOLERANCE,
            "allocation_balance_residual_m2": self.allocation_balance_residual_m2.tolist(),
            "parent_number": self.parent_number.tolist(), "pool_number": self.pool_number.tolist(),
            "pool_allocation_number_residual": _represented(np.asarray(self.parent_number, dtype=LD)-np.sum(np.asarray(self.pool_number, dtype=LD), axis=1, dtype=LD), "allocation number residual").tolist(),
            "pools": [{"name": item.reservoir.name, "state_id": item.reservoir.state_id, "allocation_provenance": item.provenance, "allocation_status": item.status, "reference_coverage_m2": item.reservoir.reference_coverage_m2, "max_coverage_m2": item.reservoir.max_coverage_m2, "kinetic_provenance": item.reservoir.kinetic_provenance, "energy_provenance": item.reservoir.energy_provenance, "reservoir": asdict(item.reservoir)} for item in self.parent.allocations],
            "kinetic_reference_policy": "supplied reservoir reference_coverage_m2 remains fixed across all coverage knots",
            "max_dimensionless_coverage_effect": self.max_dimensionless_coverage_effect.tolist(),
            "qualified_max_dimensionless_coverage_effect": MAX_DIMENSIONLESS_COVERAGE_EFFECT,
            "coupling_preflight": "maximum initial pool density across knots at each patch's minimum supplied temperature; conservative bound over all later depleted inventories and linearly interpolated temperatures",
            "subdivisions_per_interval": self.subdivisions_per_interval,
            "x_edges_m": self.surface.x_edges_m.tolist(), "z_edges_m": self.surface.z_edges_m.tolist(),
            "patch_areas_m2": areas.tolist(), "surface_area_m2": float(np.sum(areas, dtype=LD)),
            "thermal_input_status": self.surface.attributes.get("thermal_input_status", "unspecified_in_supplied_surface"),
            "memory_limits": {"max_coverage_knots": MAX_COVERAGE_KNOTS, "max_sweep_interval_reservoir_patches": MAX_SWEEP_INTERVAL_RESERVOIR_PATCHES},
            "units": {"coverage_ml": "explicitly defined parent monolayers", "areal_density_m2": "parent m^-2", "parent_number": "parent", "pool_number": "initial parent population assigned to pool", "patch_areas_m2": "m^2"},
        }


@dataclass(frozen=True, slots=True, eq=False)
class CoverageSourceSweep:
    plan: CoverageInventoryPlan
    histories: tuple[TemkinPlanarReleaseHistory, ...]

    def __post_init__(self):
        if not isinstance(self.plan, CoverageInventoryPlan):
            raise ValueError("source sweeps require a validated CoverageInventoryPlan")
        try:
            histories = tuple(self.histories)
        except TypeError as error:
            raise ValueError("source sweeps require validated TemkinPlanarReleaseHistory records") from error
        if len(histories) != len(self.plan.coverage_ml) or any(not isinstance(history, TemkinPlanarReleaseHistory) for history in histories):
            raise ValueError("source sweep must retain one law-specific history per supplied coverage knot")
        pools = tuple(item.reservoir for item in self.plan.parent.allocations)
        partitions, expected_edges = _time_partitions(self.plan.surface, self.plan.subdivisions_per_interval)
        expected_temperature = np.concatenate([
            (1-fraction[:, None, None])*self.plan.surface.temperature_k[index]
            + fraction[:, None, None]*self.plan.surface.temperature_k[index+1]
            for index, partition in enumerate(partitions)
            for fraction in [(partition[:-1]-partition[0])/(partition[-1]-partition[0])
                             + .5*np.diff(partition)/(partition[-1]-partition[0])]
        ])
        attributes = _plain(_planar_attributes(self.plan.surface.attributes, self.plan.surface.x_edges_m, self.plan.surface.z_edges_m))
        for index, history in enumerate(histories):
            if (history.reservoirs != pools or not np.array_equal(history.x_edges_m, self.plan.surface.x_edges_m)
                    or not np.array_equal(history.z_edges_m, self.plan.surface.z_edges_m)
                    or not np.array_equal(history.coverage_m2[0], self.plan.initial_pool_coverage_m2[index])
                    or not np.array_equal(history.time_edges_s, expected_edges)
                    or not np.array_equal(history.temperature_k, expected_temperature)
                    or _plain(history.surface_attributes) != attributes):
                raise ValueError("source sweep history changes its fixed pools, geometry, initial inventory, temperature/time mixture or source metadata")
        object.__setattr__(self, "histories", histories)

    @property
    def coverage_ml(self):
        return self.plan.coverage_ml

    def manifest(self):
        records = []
        areas = np.asarray(self.plan.surface.patch_areas_m2, dtype=LD)
        for index, history in enumerate(self.histories):
            remaining = np.sum(np.asarray(history.coverage_m2[-1], dtype=LD)*areas, axis=(-2, -1), dtype=LD)
            released = np.sum(np.asarray(history.released_m2, dtype=LD)*areas, axis=(0, 2, 3), dtype=LD)
            residual = np.asarray(self.plan.pool_number[index], dtype=LD)-np.asarray(remaining, dtype=LD)-released
            records.append({
                "coverage_ml": float(self.coverage_ml[index]), "initial_parent_number": float(self.plan.parent_number[index]),
                "initial_pool_number": self.plan.pool_number[index].tolist(), "remaining_pool_number": _represented(remaining, "remaining pool number").tolist(),
                "released_pool_number": _represented(released, "released pool number").tolist(),
                "pool_number_balance_residual": _represented(residual, "source inventory balance residual").tolist(),
                "energy_j": {name: _represented(np.sum(np.asarray(getattr(history, name), dtype=LD)*areas, axis=(0, 2, 3), dtype=LD), "source energy").tolist() for name in ("translational_energy_j_m2", "internal_energy_j_m2", "binding_energy_j_m2", "total_wall_energy_debit_j_m2")},
                "source_schema": history.schema, "source_law_identity": history.law_identity,
            })
        return {"schema": SWEEP_SCHEMA, "inventory": self.plan.manifest(), "sources": records,
                "limits": ["One physical parent inventory; allocations are initial independent-pool hypotheses, not reaction branches or daughter yields.", "Prescribed temperature and diagnostic wall debit; no thermal feedback, gap evolution, collisions, radiation or spectrum synthesis.", "Matched width controls redistribute existing interval/state/current source moments; no new surface is independently depleted or temperature-averaged.", "Layer convention, spatial profile, allocation and kinetic/energy parameters are supplied hypotheses; source construction does not calibrate them."]}


def prepare_coverage_inventory(surface: PlanarSurfaceHistory, parent: ParentSurfaceInventory,
                               coverage_ml, *, subdivisions_per_interval=1) -> CoverageInventoryPlan:
    """Validate every knot's parent allocation and conservative kinetic domain.

    The entire sweep is checked before any depletion integration is called.
    Positive inventories or extensive counts rounding to zero are rejected;
    no density floor, normalized multiplier or repaired allocation is inserted.
    """
    if not isinstance(surface, PlanarSurfaceHistory) or not isinstance(parent, ParentSurfaceInventory):
        raise ValueError("surface and parent must be validated PlanarSurfaceHistory and ParentSurfaceInventory records")
    surface = PlanarSurfaceHistory(surface.time_s, surface.x_edges_m, surface.z_edges_m, surface.temperature_k, surface.attributes)
    _planar_attributes(surface.attributes, surface.x_edges_m, surface.z_edges_m)
    if isinstance(subdivisions_per_interval, (bool, np.bool_)) or not isinstance(subdivisions_per_interval, Integral) or subdivisions_per_interval < 1:
        raise ValueError("subdivisions_per_interval must be a positive integer")
    subdivisions = int(subdivisions_per_interval)
    coverage = _array(coverage_ml, "coverage_ml", 1)
    if len(coverage) > MAX_COVERAGE_KNOTS or np.any(coverage[1:] <= coverage[:-1]):
        raise ValueError("coverage_ml must be strictly increasing, nonnegative and within the sweep knot limit")
    nr, nz, nx = len(parent.allocations), surface.z_edges_m.size-1, surface.x_edges_m.size-1
    nt = (surface.time_s.size-1)*subdivisions
    _memory_limit(nt, nr, nz, nx)
    if len(coverage)*nt*nr*nz*nx > MAX_SWEEP_INTERVAL_RESERVOIR_PATCHES:
        raise ValueError("coverage source sweep exceeds the aggregate interval/reservoir/patch memory limit")
    _time_partitions(surface, subdivisions)
    shape = (nz, nx)
    multiplier = _broadcast(parent.patch_multiplier, shape, "patch_multiplier")
    fractions = _immutable(np.stack([_broadcast(item.fraction, shape, "allocation fraction") for item in parent.allocations]))
    if np.any(np.abs(np.sum(np.asarray(fractions, dtype=LD), axis=0, dtype=LD)-1) > LD(FRACTION_SUM_TOLERANCE)):
        raise ValueError("initial allocation fractions must sum to one for the single physical parent at every patch; no normalization is applied")
    initial_parent = _represented(np.asarray(coverage, dtype=LD)[:, None, None]*LD(parent.monolayer.areal_density_m2)*np.asarray(multiplier, dtype=LD)[None], "initial parent density")
    initial_pools = _represented(np.asarray(initial_parent, dtype=LD)[:, None]*np.asarray(fractions, dtype=LD)[None], "initial pool density")
    residual = _represented(np.asarray(initial_parent, dtype=LD)-np.sum(np.asarray(initial_pools, dtype=LD), axis=1, dtype=LD), "initial allocation residual")
    areas = surface.patch_areas_m2
    _represented(np.asarray(initial_parent, dtype=LD)*np.asarray(areas, dtype=LD), "initial patch parent number")
    _represented(np.asarray(initial_pools, dtype=LD)*np.asarray(areas, dtype=LD), "initial patch pool number")
    parent_number = _number(initial_parent, areas, "initial parent number")
    pool_number = _number(initial_pools, areas, "initial pool number")
    minimum_t = np.min(surface.temperature_k, axis=0)
    maximum_coupling = []
    for pool_index, allocation in enumerate(parent.allocations):
        pool = allocation.reservoir
        maximum = np.max(initial_pools[:, pool_index], axis=0)
        if np.any(maximum > pool.max_coverage_m2):
            raise ValueError(f"initial allocation for pool {pool.name!r} exceeds its supplied max_coverage_m2")
        effects = [_thermal_parameters(pool, float(value), float(temperature))[0] for value, temperature in zip(maximum.flat, minimum_t.flat)]
        maximum_coupling.append(max(effects))
    return CoverageInventoryPlan(surface, parent, coverage, subdivisions, multiplier, fractions, initial_parent, initial_pools,
                                 parent_number, pool_number, residual, _immutable(maximum_coupling))


def build_coverage_source_sweep(surface: PlanarSurfaceHistory, parent: ParentSurfaceInventory,
                                coverage_ml, *, parameter_provenance,
                                subdivisions_per_interval=1) -> CoverageSourceSweep:
    """Integrate the preflighted inventory knots with unchanged supplied pools."""
    _text(parameter_provenance, "parameter_provenance")
    plan = prepare_coverage_inventory(surface, parent, coverage_ml, subdivisions_per_interval=subdivisions_per_interval)
    pools = tuple(item.reservoir for item in parent.allocations)
    histories = tuple(evolve_temkin_planar_release(plan.surface, initial, pools, parameter_provenance=parameter_provenance,
                                                  subdivisions_per_interval=plan.subdivisions_per_interval)
                      for initial in plan.initial_pool_coverage_m2)
    return CoverageSourceSweep(plan, histories)
