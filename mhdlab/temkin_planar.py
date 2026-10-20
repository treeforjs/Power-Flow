"""Qualified coverage-dependent release on local stationary (z,x) patches.

Each reservoir uses its own local population and explicitly supplied reference
coverage. This independent-pool hypothesis includes no shared-site competition.
Nonlinear depletion precedes spatial averaging; temperatures are never reduced
over the current extent before evaluating the law. Prescribed wall temperatures
make the reported wall-energy debit a diagnostic rather than thermal feedback.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
from typing import Any

import numpy as np

from .flash_surface import SCHEMA as SOURCE_SCHEMA
from .planar_release import (
    PLANAR_SURFACE_SCHEMA, PlanarSurfaceHistory, _areas, _edges, _immutable,
    _plain, _source_attributes, _text,
)
from .temkin_history import (
    DIAGNOSTIC_FIELDS, RELEASE_FIELDS, STEP_FIELDS, TIME_MODEL,
    TemkinReleaseHistory, _array,
)
from .temkin_release import LAW_IDENTITY, TemkinReservoir, temkin_release_step
from .thermal_release import SurfaceInventory


SCHEMA = "mhdlab-temkin-planar-release-v1"
TEMKIN_PLANAR_RELEASE_SCHEMA = SCHEMA
COVERAGE_MODEL = "independent_reservoir_local_coverage_and_reference_density_without_shared_site_competition"
WIDTH_UNIFORM_MODEL = "width_redistribution_preserving_each_interval_reservoir_current_patch_number_and_energy"
FULL_FACE_UNIFORM_MODEL = "full_face_redistribution_preserving_each_interval_reservoir_number_and_energy"
PATCH_TEMPERATURE_APPROXIMATION = "local x/z patch averages; refine x,z,time before nonlinear release"
MAX_INTERVALS = 100_000
MAX_INTERVAL_RESERVOIR_PATCHES = 20_000_000
MIN_POSITIVE_UNIFORM_MOMENT_DENSITY = float(np.finfo(float).tiny)


def _pools(reservoirs):
    try:
        pools = tuple(reservoirs)
    except TypeError as error:
        raise ValueError("reservoirs must contain TemkinReservoir records") from error
    if not pools or any(not isinstance(pool, TemkinReservoir) for pool in pools):
        raise ValueError("reservoirs must contain TemkinReservoir records")
    if len({pool.name for pool in pools}) != len(pools):
        raise ValueError("Temkin reservoir names must be unique")
    return pools


def _memory_limit(nt, nr, nz, nx):
    if nt > MAX_INTERVALS or nt * nr * nz * nx > MAX_INTERVAL_RESERVOIR_PATCHES:
        raise ValueError("Temkin planar history exceeds the reference memory limit")


def _planar_attributes(value, x_edges, z_edges):
    """Check resolved patch geometry while retaining every nested source field."""
    attributes = _source_attributes(value)
    if attributes["schema"] != PLANAR_SURFACE_SCHEMA:
        raise ValueError("Temkin planar release requires resolved planar source metadata")
    if attributes["temperature_location"] == "material_sample_plane":
        raise ValueError("a material sample plane is not a boundary temperature; explicit justification is required")
    intervals = [[float(lo), float(hi)] for lo, hi in zip(z_edges[:-1], z_edges[1:])]
    if (_plain(attributes.get("current_patch_intervals_m")) != intervals
            or attributes.get("current_direction_representation") != "resolved_piecewise_constant_surface_patches"
            or attributes.get("patch_temperature_approximation") != PATCH_TEMPERATURE_APPROXIMATION):
        raise ValueError("Temkin planar source metadata must preserve every resolved current patch and local temperature approximation")
    extraction = attributes.get("resolved_patch_extraction")
    if extraction is not None:
        if not isinstance(extraction, Mapping):
            raise ValueError("resolved patch extraction metadata must be a mapping")
        manifests = extraction.get("patch_extraction_manifests")
        if not isinstance(manifests, tuple) or len(manifests) != len(intervals):
            raise ValueError("resolved patch extraction must retain one manifest per current patch")
        for manifest, interval in zip(manifests, intervals):
            if not isinstance(manifest, Mapping) or not isinstance(manifest.get("sampling"), Mapping):
                raise ValueError("every resolved patch extraction must retain its sampling metadata")
            sampling = manifest["sampling"]
            if (_plain(sampling.get("current_interval_m")) != interval
                    or _plain(sampling.get("x_edges_m")) != x_edges.tolist()
                    or sampling.get("material_sample_depth_m") != attributes["material_sample_depth_m"]):
                raise ValueError("nested extraction geometry or sampling depth disagrees with the planar source")
    return attributes


def _batched_metadata(attributes):
    """Adapt metadata only for independent-patch kinetic/energy validation.

    The flattened coordinate is an index, never a physical width. No spatial
    integral or uniform-source construction uses this temporary representation.
    Original resolved metadata remains authoritative in the planar record.
    """
    adapted = _plain(attributes)
    adapted["schema"] = SOURCE_SCHEMA
    adapted["planar_batch_coordinate_interpretation"] = "independent_patch_index_for_kinetic_validation_only_no_physical_width"
    return adapted


def _weighted_mean(source, weights, total, axes):
    """Form positive weighted means without underflowing normalized weights.

    Float64 mantissas/exponents keep a small patch's large source contribution
    until after the spatial sum. Terms lost relative to the largest positive
    term are below floating-point significance, rather than lost merely because
    a geometric fraction is below the subnormal range.

    A positive mean must be at least float64's minimum normal value in each
    dataset's declared SI units. Subnormal means are rejected because their
    rounding can substantially change the reconstructed number/energy moment.
    Exact zero moments remain supported. HDF metadata records this restriction.
    """
    source_mantissa, source_exponent = np.frexp(source)
    weight_mantissa, weight_exponent = np.frexp(weights)
    exponent = source_exponent + weight_exponent
    positive = source > 0
    pivot = np.max(np.where(positive, exponent, -4096), axis=axes, keepdims=True)
    any_positive = np.any(positive, axis=axes, keepdims=True)
    pivot = np.where(any_positive, pivot, 0)
    total_mantissa, total_exponent = math.frexp(float(total))
    with np.errstate(over="ignore", invalid="ignore", under="ignore"):
        scaled = np.ldexp(source_mantissa * weight_mantissa, exponent - pivot)
        mean = np.ldexp(np.sum(scaled, axis=axes, keepdims=True) / total_mantissa,
                        pivot - total_exponent)
    mean = _array(mean, "weighted mean")
    if np.any(any_positive & (mean < MIN_POSITIVE_UNIFORM_MOMENT_DENSITY)):
        raise ValueError("uniform source density is below the qualified normal range; subnormal rounding cannot preserve relative source moments")
    return mean


@dataclass(frozen=True, slots=True, eq=False)
class TemkinPlanarReleaseHistory:
    """Law-specific local source moments, shape (interval,reservoir,z,x).

    Coverage includes both time boundaries. Construction replays each interval
    through the qualified Temkin record, including kinetic-domain, inventory,
    energy and numerical-diagnostic checks. Uniform controls redistribute source
    moments; they are not separately depleted surfaces or velocity distributions.
    """

    time_edges_s: np.ndarray
    x_edges_m: np.ndarray
    z_edges_m: np.ndarray
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
    coverage_model: str = COVERAGE_MODEL
    schema: str = SCHEMA

    def __post_init__(self):
        if (self.law_identity != LAW_IDENTITY or self.time_model != TIME_MODEL
                or self.coverage_model != COVERAGE_MODEL or self.schema != SCHEMA):
            raise ValueError("unsupported Temkin planar schema, law, time or independent-pool coverage model")
        for name in ("time_edges_s", "x_edges_m", "z_edges_m"):
            object.__setattr__(self, name, _edges(getattr(self, name), name, time=name == "time_edges_s"))
        pools = _pools(self.reservoirs)
        object.__setattr__(self, "reservoirs", pools)
        nt, nr = self.time_edges_s.size - 1, len(pools)
        nz, nx = self.z_edges_m.size - 1, self.x_edges_m.size - 1
        _memory_limit(nt, nr, nz, nx)
        for name in ("temperature_k", "coverage_m2", *STEP_FIELDS, "inversion_log_time_residual"):
            object.__setattr__(self, name, _array(getattr(self, name), name))
        object.__setattr__(self, "inversion_iterations", _array(self.inversion_iterations, "inversion_iterations", integer=True, maximum=128))
        object.__setattr__(self, "inversion_regime", _array(self.inversion_regime, "inversion_regime", integer=True, maximum=4))
        if self.temperature_k.shape != (nt, nz, nx) or np.any(self.temperature_k <= 0):
            raise ValueError("temperature_k must be positive with shape (interval,current,width)")
        if self.coverage_m2.shape != (nt + 1, nr, nz, nx):
            raise ValueError("coverage_m2 must have shape (boundary,reservoir,current,width)")
        if any(getattr(self, name).shape != (nt, nr, nz, nx) for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS)):
            raise ValueError("Temkin planar source and diagnostic arrays must have shape (interval,reservoir,current,width)")
        attributes = _planar_attributes(self.surface_attributes, self.x_edges_m, self.z_edges_m)
        object.__setattr__(self, "surface_attributes", attributes)
        _text(self.parameter_provenance, "parameter_provenance")
        _areas(self.x_edges_m, self.z_edges_m)
        # A law-specific validation record batches independent patches. Its
        # temporary index edges are not exposed as physical source geometry.
        batch = TemkinReleaseHistory(
            time_edges_s=self.time_edges_s, x_edges_m=np.arange(nz * nx + 1, dtype=float),
            temperature_k=self.temperature_k.reshape(nt, nz * nx),
            coverage_m2=self.coverage_m2.reshape(nt + 1, nr, nz * nx),
            **{name: getattr(self, name).reshape(nt, nr, nz * nx) for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS)},
            reservoirs=pools, surface_attributes=_batched_metadata(attributes),
            parameter_provenance=self.parameter_provenance, law_identity=self.law_identity, time_model=self.time_model,
        )
        # Retain canonical, strictly immutable integer diagnostic dtypes.
        for name in DIAGNOSTIC_FIELDS:
            object.__setattr__(self, name, getattr(batch, name).reshape(nt, nr, nz, nx))

    @property
    def patch_areas_m2(self):
        return _areas(self.x_edges_m, self.z_edges_m)

    @property
    def average_flux_m2_s(self):
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            return _array(self.released_m2 / np.diff(self.time_edges_s)[:, None, None, None], "average_flux_m2_s")

    def patch_extensive(self, quantity="released_m2"):
        """Use actual x/z area once; coverage and balance residuals are allowed."""
        if quantity not in (*STEP_FIELDS, "coverage_m2"):
            raise ValueError("unsupported Temkin planar extensive quantity")
        with np.errstate(over="ignore", invalid="ignore"):
            return _array(getattr(self, quantity) * self.patch_areas_m2, "patch_extensive")

    def extensive_totals(self, quantity="released_m2"):
        """Integrate actual face area, retaining interval and reservoir."""
        with np.errstate(over="ignore", invalid="ignore"):
            return _array(np.sum(self.patch_extensive(quantity), axis=(-2, -1)), "extensive_totals")

    def matched_uniform_width(self):
        """Preserve each time/reservoir/current-patch number and energy moment."""
        widths = np.diff(self.x_edges_m)
        total = self.x_edges_m[-1] - self.x_edges_m[0]
        result = {}
        for name in RELEASE_FIELDS:
            source = getattr(self, name)
            mean = _weighted_mean(source, widths, total, (-1,))
            result[name] = _array(np.broadcast_to(mean, source.shape), name)
        return result

    def matched_uniform_full_face(self):
        """Explicit alternate control that also erases current-direction profiles."""
        areas = self.patch_areas_m2
        result = {}
        for name in RELEASE_FIELDS:
            source = getattr(self, name)
            mean = _weighted_mean(source, areas, np.sum(areas), (-2, -1))
            result[name] = _array(np.broadcast_to(mean, source.shape), name)
        return result


def evolve_temkin_planar_release(
    surface: PlanarSurfaceHistory,
    coverage_m2: np.ndarray,
    reservoirs: Sequence[TemkinReservoir],
    *,
    parameter_provenance: str,
    subdivisions_per_interval: int = 1,
) -> TemkinPlanarReleaseHistory:
    """Deplete independent local pools using qualified frozen-midpoint solves.

    Flattening only batches the (z,x) patches through the same nonlinear law;
    temperature, coverage, reservoir and current patches are never averaged.
    Refine width/current patches and temperature-history cadence separately.
    """
    if not isinstance(surface, PlanarSurfaceHistory):
        raise ValueError("surface must be a PlanarSurfaceHistory")
    _text(parameter_provenance, "parameter_provenance")
    if (isinstance(subdivisions_per_interval, (bool, np.bool_))
            or not isinstance(subdivisions_per_interval, (int, np.integer)) or subdivisions_per_interval < 1):
        raise ValueError("subdivisions_per_interval must be a positive integer")
    # Revalidate geometry and metadata even if the caller bypassed immutability.
    surface = PlanarSurfaceHistory(surface.time_s, surface.x_edges_m, surface.z_edges_m,
                                  surface.temperature_k, surface.attributes)
    attributes = _planar_attributes(surface.attributes, surface.x_edges_m, surface.z_edges_m)
    pools = _pools(reservoirs)
    nr, nz, nx = len(pools), surface.z_edges_m.size - 1, surface.x_edges_m.size - 1
    count = (surface.time_s.size - 1) * int(subdivisions_per_interval)
    _memory_limit(count, nr, nz, nx)
    coverage = _immutable(coverage_m2, "coverage_m2", 3, nonnegative=True)
    if coverage.shape != (nr, nz, nx):
        raise ValueError("coverage_m2 must have shape (reservoir,current,width)")
    partitions = [np.linspace(lo, hi, int(subdivisions_per_interval) + 1)
                  for lo, hi in zip(surface.time_s[:-1], surface.time_s[1:])]
    edges = np.concatenate((surface.time_s[:1], *(partition[1:] for partition in partitions)))
    if np.any(edges[1:] <= edges[:-1]):
        raise ValueError("time subdivision is not representable at this time origin")
    temperatures = np.concatenate([
        (1 - fraction[:, None, None]) * surface.temperature_k[index]
        + fraction[:, None, None] * surface.temperature_k[index + 1]
        for index, partition in enumerate(partitions)
        for fraction in [(partition[:-1] - partition[0]) / (partition[-1] - partition[0])
                         + .5 * np.diff(partition) / (partition[-1] - partition[0])]
    ])
    current = SurfaceInventory(coverage.reshape(nr, nz * nx))
    coverage_edges = [coverage]
    values = {name: [] for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS)}
    for temperature, dt in zip(temperatures, np.diff(edges)):
        trial = temkin_release_step(current, pools, temperature.reshape(nz * nx), float(dt))
        if trial.law_identity != LAW_IDENTITY:
            raise ValueError("core trial has an unsupported Temkin law identity")
        current = trial.trial_inventory
        coverage_edges.append(current.coverage_m2.reshape(nr, nz, nx))
        for name in values:
            values[name].append(getattr(trial, name).reshape(nr, nz, nx))
    return TemkinPlanarReleaseHistory(
        edges, surface.x_edges_m, surface.z_edges_m, temperatures, np.stack(coverage_edges),
        **{name: np.stack(items) for name, items in values.items()}, reservoirs=pools,
        surface_attributes=attributes, parameter_provenance=parameter_provenance,
    )
