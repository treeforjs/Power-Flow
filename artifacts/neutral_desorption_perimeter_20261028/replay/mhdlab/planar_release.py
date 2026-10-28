"""Local release on a stationary planar surface resolved in width and current.

Temperature and reservoir coverage remain correlated on every (z,x) patch.
Nonlinear release and depletion precede area integration. Patch temperatures are
piecewise constant spatial approximations; refine x, z, and time before using
their source totals. No temperature conversion, radiation or transport is inferred.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import json
import re
from numbers import Real
from types import MappingProxyType
from typing import Any

import numpy as np

from .flash_surface import COORDINATE_SYSTEM, SCHEMA as SOURCE_SCHEMA, SURFACE_MODEL
from .surface_sources import matched_uniform_moments
from .thermal_release import BindingReservoir, SurfaceInventory, energy_budgets, release_step
from .surface_coordinates import MFEM_SIDE_COORDINATE_SYSTEM, validate_mfem_face_axes


PLANAR_SURFACE_SCHEMA = "mhdlab-planar-surface-v1"
MFEM_SURFACE_SCHEMA = "mhdlab-mfem-surface-v1"
PLANAR_RELEASE_SCHEMA = "mhdlab-planar-release-v1"
_RELEASE_FIELDS = (
    "released_m2", "translational_energy_j_m2", "internal_energy_j_m2",
    "binding_energy_j_m2", "total_wall_energy_debit_j_m2",
)
_STEP_FIELDS = (*_RELEASE_FIELDS, "inventory_balance_residual_m2")


def _immutable(value: Any, name: str, rank: int, *, nonnegative: bool = False) -> np.ndarray:
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked values")
    try:
        raw = np.asarray(value, dtype=object)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a rectangular real array") from exc
    if raw.ndim != rank or any(size == 0 for size in raw.shape):
        raise ValueError(f"{name} must be a nonempty rank-{rank} array")
    for item in raw.flat:
        if isinstance(item, (bool, np.bool_)) or not isinstance(item, Real):
            raise ValueError(f"{name} must contain real numbers, not booleans or complex values")
    try:
        with np.errstate(over="raise", invalid="raise"):
            result = np.asarray(raw, dtype=np.float64)
    except (TypeError, ValueError, OverflowError, FloatingPointError) as exc:
        raise ValueError(f"{name} must contain finite real numbers") from exc
    if not np.all(np.isfinite(result)) or (nonnegative and np.any(result < 0.0)):
        raise ValueError(f"{name} must contain finite {'nonnegative ' if nonnegative else ''}values")
    return np.frombuffer(result.tobytes(), dtype=np.float64).reshape(result.shape)


def _edges(value: Any, name: str, *, time: bool = False) -> np.ndarray:
    edges = _immutable(value, name, 1, nonnegative=time)
    if edges.size < 2 or np.any(edges[1:] <= edges[:-1]):
        raise ValueError(f"{name} must contain at least two strictly increasing values")
    with np.errstate(over="ignore"):
        span = edges[-1] - edges[0]
    if not np.isfinite(span):
        raise ValueError(f"{name} must span a finite interval")
    return edges


def _text(value: Any, name: str, *, empty: bool = False) -> str:
    if not isinstance(value, str) or (not empty and not value.strip()):
        raise ValueError(f"{name} must be an explicit {'string' if empty else 'nonempty string'}")
    return value


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) or not key for key in value):
            raise ValueError("source metadata keys must be nonempty strings")
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, np.ndarray):
        return _freeze(value.tolist())
    if isinstance(value, (tuple, list)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, np.generic):
        return _freeze(value.item())
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float) and np.isfinite(value):
        return value
    raise ValueError("source metadata must contain finite JSON-compatible values")


def _plain(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_plain(item) for item in value]
    return value


def _mfem_attributes(attributes: Mapping[str, Any], source_kind: str) -> None:
    """Apply MFEM semantics to direct construction, reload and release alike."""
    def evidence(value, name):
        value = _text(value, name)
        normalized = value.strip().lower()
        if normalized.startswith("replace_") or normalized in {
                "pending", "unknown", "unavailable", "none", "n/a", "tbd", "missing"}:
            raise ValueError(f"{name} cannot be placeholder evidence")
        return value

    for key in ("run_id", "source_reference", "surface_region", "unit_mapping_reference",
                "source_clock_reference", "temperature_zero_convention_reference"):
        evidence(attributes.get(key), key)
    if source_kind == "mfem":
        evidence(attributes.get("mfem_version"), "mfem_version")
        revision = evidence(attributes.get("mfem_source_revision"), "mfem_source_revision")
        if re.fullmatch(r"[0-9a-f]{40}", revision) is None:
            raise ValueError("MFEM source revision requires a full source commit")
    else:
        for key in ("mfem_version", "mfem_source_revision"):
            _text(attributes.get(key), key, empty=True)
    location = attributes.get("temperature_location")
    expected = ("one_sided_material_boundary_trace_patch_average" if location == "surface_boundary"
                else "resolved_material_sample_plane_patch_average")
    if attributes.get("temperature_evaluation") != expected:
        raise ValueError("MFEM temperature evaluation must preserve its one-sided trace or sampled-plane association")
    basis = attributes.get("temperature_value_basis")
    offset = attributes.get("initial_temperature_offset_k")
    if isinstance(offset, bool) or not isinstance(offset, Real) or not np.isfinite(offset):
        raise ValueError("MFEM initial temperature offset must be finite real data")
    if (basis == "absolute_kelvin" and offset != 0 or
            basis == "temperature_increment_kelvin" and offset <= 0 or
            basis not in {"absolute_kelvin", "temperature_increment_kelvin"}):
        raise ValueError("MFEM temperature basis and initial offset are inconsistent")
    scope = attributes.get("mfem_solver_scope")
    if not isinstance(scope, Mapping) or scope.get("stationary_material") is not True:
        raise ValueError("MFEM surface requires stationary solver scope")
    for key in ("conductivity_temperature_feedback", "hall_induction", "material_motion"):
        if type(scope.get(key)) is not bool:
            raise ValueError("MFEM solver scope must explicitly declare supported physics")
    if scope["material_motion"]:
        raise ValueError("moving MFEM material requires a different surface contract")
    evidence(scope.get("qualification_reference"), "MFEM solver qualification")
    face = attributes.get("mfem_face")
    if not isinstance(face, Mapping):
        raise ValueError("MFEM source requires its physical face metadata")
    validate_mfem_face_axes(face)
    for key in ("coordinate_system", "surface_model", "surface_region"):
        if face.get(key) != attributes.get(key):
            raise ValueError("MFEM face metadata conflicts with the surface association")
    if face.get("current_direction_representation") != "resolved_patches":
        raise ValueError("MFEM face must retain resolved current patches")
    if type(face.get("outward_normal_sign")) is not int or face["outward_normal_sign"] not in {-1, 1}:
        raise ValueError("MFEM face requires a signed local normal")
    if type(face.get("material_boundary_attribute")) is not int or face["material_boundary_attribute"] <= 0:
        raise ValueError("MFEM face requires its positive material boundary attribute")
    position = face.get("normal_position_m")
    if isinstance(position, bool) or not isinstance(position, Real) or not np.isfinite(position):
        raise ValueError("MFEM face position must be finite real data")
    for key in ("material_selection_reference", "full_width_coverage_reference"):
        evidence(face.get(key), key)
    for key in ("sample_table_sha256", "sample_configuration_sha256"):
        value = attributes.get(key)
        if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
            raise ValueError("MFEM surface must preserve sample/configuration SHA256 identities")


def _source_attributes(value: Mapping[str, Any]) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError("source_attributes must be a metadata mapping")
    attributes = _freeze(value)
    for key, expected in (("surface_model", SURFACE_MODEL), ("material_sample_depth_units", "m")):
        if attributes.get(key) != expected:
            raise ValueError(f"source_attributes must preserve {key}={expected!r}")
    schema = _text(attributes.get("schema"), "source schema")
    if schema not in {SOURCE_SCHEMA, MFEM_SURFACE_SCHEMA, PLANAR_SURFACE_SCHEMA}:
        raise ValueError("source_attributes must identify the supported surface metadata schema")
    original_schema = attributes.get("original_source_schema", schema)
    if schema == PLANAR_SURFACE_SCHEMA and original_schema not in {SOURCE_SCHEMA, MFEM_SURFACE_SCHEMA}:
        raise ValueError("planar metadata must retain its original supported source schema")
    allowed_coordinates = ({COORDINATE_SYSTEM, MFEM_SIDE_COORDINATE_SYSTEM}
                           if original_schema == MFEM_SURFACE_SCHEMA else {COORDINATE_SYSTEM})
    if attributes.get("coordinate_system") not in allowed_coordinates:
        raise ValueError("source_attributes must preserve its supported coordinate_system")
    for key in ("run_id", "source_reference", "surface_region"):
        _text(attributes.get(key), key)
    source_kind = _text(attributes.get("source_kind"), "source_kind")
    if source_kind not in {"flash", "mfem", "synthetic"}:
        raise ValueError("source_attributes must explicitly identify flash, mfem or synthetic data")
    if original_schema == MFEM_SURFACE_SCHEMA:
        if source_kind not in {"mfem", "synthetic"} or attributes.get("heating_backend") != "mfem_joule":
            raise ValueError("MFEM surface metadata must identify the MFEM Joule backend")
        _mfem_attributes(attributes, source_kind)
    else:
        if source_kind == "mfem":
            raise ValueError("MFEM data must retain its own supported source schema")
        _text(attributes.get("flash_version"), "flash_version", empty=source_kind != "flash")
    location = _text(attributes.get("temperature_location"), "temperature_location")
    if location not in {"surface_boundary", "material_sample_plane", "boundary_temperature_proxy"}:
        raise ValueError("source_attributes must preserve the temperature association")
    depth = attributes.get("material_sample_depth_m")
    if isinstance(depth, bool) or not isinstance(depth, (float, int)):
        raise ValueError("material_sample_depth_m must be finite and nonnegative")
    try:
        depth = float(depth)
    except OverflowError as exc:
        raise ValueError("material_sample_depth_m must be finite and nonnegative") from exc
    if not np.isfinite(depth) or depth < 0:
        raise ValueError("material_sample_depth_m must be finite and nonnegative")
    reference = _text(attributes.get("boundary_temperature_reference"), "boundary_temperature_reference", empty=True)
    if location == "surface_boundary" and depth != 0:
        raise ValueError("a boundary temperature must retain zero sampling depth")
    if location != "surface_boundary" and depth == 0:
        raise ValueError("material-plane and proxy temperatures require positive sampling depth")
    if location == "material_sample_plane" and reference:
        raise ValueError("a material sample plane cannot claim a boundary conversion reference")
    if location == "boundary_temperature_proxy" and not reference.strip():
        raise ValueError("a boundary proxy requires its explicit justification reference")
    averaging = attributes.get("current_averaging_interval_m_json")
    if not isinstance(averaging, str):
        raise ValueError("source_attributes must explicitly preserve current averaging metadata")
    try:
        interval = json.loads(averaging)
    except (TypeError, ValueError) as exc:
        raise ValueError("invalid current averaging metadata") from exc
    if interval is not None or attributes.get("current_averaging_interval_m") is not None:
        raise ValueError("a current-averaged trace cannot be inflated into resolved planar patches")
    _text(attributes.get("current_uniformity_reference"), "current_uniformity_reference", empty=True)
    return attributes


def _areas(x_edges: np.ndarray, z_edges: np.ndarray) -> np.ndarray:
    try:
        with np.errstate(over="raise", invalid="raise"):
            area = np.diff(z_edges)[:, None] * np.diff(x_edges)[None, :]
            total = np.sum(area)
    except FloatingPointError as exc:
        raise ValueError("surface areas exceed the finite numeric range") from exc
    if not np.all(np.isfinite(area)) or np.any(area <= 0.0) or not np.isfinite(total):
        raise ValueError("surface patch areas must be positive and numerically resolved")
    return _immutable(area, "patch areas", 2, nonnegative=True)


@dataclass(frozen=True)
class PlanarSurfaceHistory:
    """Stationary (z,x) patch temperatures, shape (time,current,width).

    Unknown JSON-compatible source fields are retained immutably. Material-side
    planes can be represented for diagnostics, but are not boundary temperatures
    and cannot drive release. A declared proxy keeps its original positive depth
    and justification. No whole-current averaged input can become resolved data.
    """

    time_s: np.ndarray
    x_edges_m: np.ndarray
    z_edges_m: np.ndarray
    temperature_k: np.ndarray
    source_attributes: Mapping[str, Any]

    def __post_init__(self):
        time = _edges(self.time_s, "time_s", time=True)
        x = _edges(self.x_edges_m, "x_edges_m")
        z = _edges(self.z_edges_m, "z_edges_m")
        temperature = _immutable(self.temperature_k, "temperature_k", 3)
        if temperature.shape != (time.size, z.size - 1, x.size - 1) or np.any(temperature <= 0.0):
            raise ValueError("temperature_k must be positive with shape (nt,nz,nx)")
        attributes = _source_attributes(self.source_attributes)
        for name, value in (("time_s", time), ("x_edges_m", x), ("z_edges_m", z),
                            ("temperature_k", temperature), ("source_attributes", attributes)):
            object.__setattr__(self, name, value)
        _areas(x, z)
        # Refuse inconsistent preexisting patch geometry instead of overwriting it.
        self.attributes

    @property
    def patch_areas_m2(self) -> np.ndarray:
        return _areas(self.x_edges_m, self.z_edges_m)

    @property
    def attributes(self) -> dict[str, Any]:
        attributes = _plain(self.source_attributes)
        additions = {
            "current_direction_representation": "resolved_piecewise_constant_surface_patches",
            "current_patch_intervals_m": [[float(lo), float(hi)] for lo, hi in zip(self.z_edges_m[:-1], self.z_edges_m[1:])],
            "patch_temperature_approximation": "local x/z patch averages; refine x,z,time before nonlinear release",
        }
        for key, value in additions.items():
            if key in attributes and attributes[key] != value:
                raise ValueError(f"source metadata conflicts with resolved planar {key}")
        attributes.update(additions)
        attributes["original_source_schema"] = attributes.get("original_source_schema", attributes["schema"])
        attributes["schema"] = PLANAR_SURFACE_SCHEMA
        return attributes


@dataclass(frozen=True)
class PlanarReleaseHistory:
    """Time-integrated number/energy densities on local (reservoir,z,x) patches."""

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
    reservoirs: tuple[BindingReservoir, ...]
    surface_attributes: Mapping[str, Any]
    parameter_provenance: str

    def __post_init__(self):
        for name in ("time_edges_s", "x_edges_m", "z_edges_m"):
            object.__setattr__(self, name, _edges(getattr(self, name), name, time=name == "time_edges_s"))
        object.__setattr__(self, "temperature_k", _immutable(self.temperature_k, "temperature_k", 3))
        for name in ("coverage_m2", *_STEP_FIELDS):
            object.__setattr__(self, name, _immutable(getattr(self, name), name, 4,
                                                     nonnegative=name != "inventory_balance_residual_m2"))
        pools = tuple(self.reservoirs)
        if not pools or any(not isinstance(pool, BindingReservoir) for pool in pools):
            raise ValueError("reservoirs must contain BindingReservoir records")
        if len({pool.name for pool in pools}) != len(pools):
            raise ValueError("reservoir names must be unique")
        object.__setattr__(self, "reservoirs", pools)
        attributes = _source_attributes(self.surface_attributes)
        if attributes["schema"] != PLANAR_SURFACE_SCHEMA:
            raise ValueError("release history must preserve resolved planar surface metadata")
        if attributes["temperature_location"] == "material_sample_plane":
            raise ValueError("release requires a boundary temperature or an explicitly justified boundary proxy")
        object.__setattr__(self, "surface_attributes", attributes)
        _text(self.parameter_provenance, "parameter_provenance")
        nt, nr, nz, nx = self.time_edges_s.size - 1, len(pools), self.z_edges_m.size - 1, self.x_edges_m.size - 1
        if self.temperature_k.shape != (nt, nz, nx) or np.any(self.temperature_k <= 0.0):
            raise ValueError("temperature_k must be positive with shape (nsteps,nz,nx)")
        if self.coverage_m2.shape != (nt + 1, nr, nz, nx):
            raise ValueError("coverage_m2 must have shape (nsteps+1,nres,nz,nx)")
        if any(getattr(self, name).shape != (nt, nr, nz, nx) for name in _STEP_FIELDS):
            raise ValueError("release and residual arrays must have shape (nsteps,nres,nz,nx)")
        _areas(self.x_edges_m, self.z_edges_m)
        expected_intervals = [[float(lo), float(hi)] for lo, hi in zip(self.z_edges_m[:-1], self.z_edges_m[1:])]
        if (_plain(attributes.get("current_patch_intervals_m")) != expected_intervals
                or attributes.get("current_direction_representation") != "resolved_piecewise_constant_surface_patches"):
            raise ValueError("release metadata must preserve every resolved current patch interval")
        tolerance = 64.0 * np.finfo(float).eps
        if np.any(self.coverage_m2[1:] > self.coverage_m2[:-1]):
            raise ValueError("desorption-only local coverage cannot increase")
        with np.errstate(over="ignore", invalid="ignore"):
            residual = self.coverage_m2[1:] + self.released_m2 - self.coverage_m2[:-1]
        allowance = tolerance * self.coverage_m2[:-1]
        if (not np.all(np.isfinite(residual)) or np.any(np.abs(residual) > allowance)
                or np.any(np.abs(self.inventory_balance_residual_m2) > allowance)):
            raise ValueError("patch-local release and coverage do not balance")
        for interval in range(nt):
            budget = energy_budgets(self.released_m2[interval].reshape(nr, nz * nx), pools,
                                    self.temperature_k[interval].reshape(nz * nx))
            for name in _RELEASE_FIELDS[1:]:
                expected = budget[name].reshape(nr, nz, nx)
                if not np.allclose(getattr(self, name)[interval], expected, rtol=tolerance, atol=0.0):
                    raise ValueError(f"{name} disagrees with the patch-local energy convention")

    @property
    def patch_areas_m2(self) -> np.ndarray:
        return _areas(self.x_edges_m, self.z_edges_m)

    @property
    def average_flux_m2_s(self) -> np.ndarray:
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                return _immutable(self.released_m2 / np.diff(self.time_edges_s)[:, None, None, None],
                                  "average flux", 4, nonnegative=True)
        except FloatingPointError as exc:
            raise ValueError("average source flux exceeds the finite numeric range") from exc

    def patch_extensive(self, quantity: str = "released_m2") -> np.ndarray:
        """Multiply a delivered density or coverage by each actual patch area."""
        if quantity not in (*_STEP_FIELDS, "coverage_m2"):
            raise ValueError("unsupported surface quantity")
        try:
            with np.errstate(over="raise", invalid="raise"):
                return _immutable(getattr(self, quantity) * self.patch_areas_m2,
                                  "patch extensive quantity", 4)
        except FloatingPointError as exc:
            raise ValueError("patch extensive source exceeds the finite numeric range") from exc

    def extensive_totals(self, quantity: str = "released_m2") -> np.ndarray:
        """Integrate over x and z, retaining interval and reservoir dimensions."""
        try:
            with np.errstate(over="raise", invalid="raise"):
                return _immutable(np.sum(self.patch_extensive(quantity), axis=(-2, -1)),
                                  "extensive totals", 2)
        except FloatingPointError as exc:
            raise ValueError("integrated source exceeds the finite numeric range") from exc

    def matched_uniform_width(self) -> dict[str, np.ndarray]:
        """Redistribute delivered number/energy across x at each time/reservoir/z.

        The original current-direction source profile and each reservoir remain
        separate. This changes delivered moments, not stored coverage or source
        temperatures, and does not assign a joint particle velocity distribution.
        """
        return matched_uniform_moments({name: getattr(self, name) for name in _RELEASE_FIELDS}, self.x_edges_m)

    def matched_uniform_full_face(self) -> dict[str, np.ndarray]:
        """Explicit alternate control that erases both width and current profiles."""
        area = self.patch_areas_m2
        weights = area / np.sum(area)
        result = {}
        for name in _RELEASE_FIELDS:
            source = getattr(self, name)
            try:
                with np.errstate(over="raise", invalid="raise"):
                    mean = np.sum(source * weights, axis=(-2, -1), keepdims=True)
                result[name] = _immutable(np.broadcast_to(mean, source.shape), name, 4, nonnegative=True)
            except FloatingPointError as exc:
                raise ValueError("uniform full-face source exceeds the finite numeric range") from exc
        return result


def evolve_planar_release(
    surface: PlanarSurfaceHistory,
    coverage_m2: np.ndarray,
    reservoirs: Sequence[BindingReservoir],
    *,
    parameter_provenance: str,
    subdivisions_per_interval: int = 1,
) -> PlanarReleaseHistory:
    """Deplete every local reservoir at fractional midpoint temperatures.

    Flattening here only batches independent patches through release_step; no
    temperature, coverage, or reservoir averaging takes place. The frozen-T
    depletion is exact in each interval; varying-T accuracy requires refinement.
    """
    if not isinstance(surface, PlanarSurfaceHistory):
        raise ValueError("surface must be a PlanarSurfaceHistory")
    _text(parameter_provenance, "parameter_provenance")
    if (isinstance(subdivisions_per_interval, (bool, np.bool_))
            or not isinstance(subdivisions_per_interval, (int, np.integer)) or subdivisions_per_interval < 1):
        raise ValueError("subdivisions_per_interval must be a positive integer")
    if surface.source_attributes["temperature_location"] == "material_sample_plane":
        raise ValueError("a material sample plane is not a boundary temperature; explicit justification is required")
    pools = tuple(reservoirs)
    if not pools or any(not isinstance(pool, BindingReservoir) for pool in pools):
        raise ValueError("reservoirs must contain BindingReservoir records")
    if len({pool.name for pool in pools}) != len(pools):
        raise ValueError("reservoir names must be unique")
    nr, nz, nx = len(pools), surface.z_edges_m.size - 1, surface.x_edges_m.size - 1
    coverage = _immutable(coverage_m2, "coverage_m2", 3, nonnegative=True)
    if coverage.shape != (nr, nz, nx):
        raise ValueError("coverage_m2 must have shape (nres,nz,nx)")
    count = (surface.time_s.size - 1) * int(subdivisions_per_interval)
    if count > 100_000 or count * nr * nz * nx > 20_000_000:
        raise ValueError("planar release exceeds the reference memory limit")
    partitions = [np.linspace(start, stop, int(subdivisions_per_interval) + 1)
                  for start, stop in zip(surface.time_s[:-1], surface.time_s[1:])]
    time_edges = np.concatenate((surface.time_s[:1], *(partition[1:] for partition in partitions)))
    if np.any(time_edges[1:] <= time_edges[:-1]):
        raise ValueError("time subdivision is not representable at the supplied time origin")
    temperatures = np.concatenate([
        (1.0 - fraction[:, None, None]) * surface.temperature_k[index]
        + fraction[:, None, None] * surface.temperature_k[index + 1]
        for index, partition in enumerate(partitions)
        for fraction in [(partition[:-1] - partition[0]) / (partition[-1] - partition[0])
                         + 0.5 * (np.diff(partition) / (partition[-1] - partition[0]))]
    ])
    inventory = SurfaceInventory(coverage.reshape(nr, nz * nx))
    coverage_edges = [coverage]
    delivered = {name: [] for name in _STEP_FIELDS}
    for temperature, dt in zip(temperatures, np.diff(time_edges)):
        trial = release_step(inventory, pools, temperature.reshape(nz * nx), float(dt))
        inventory = trial.trial_inventory
        coverage_edges.append(inventory.coverage_m2.reshape(nr, nz, nx))
        for name in _STEP_FIELDS:
            delivered[name].append(getattr(trial, name).reshape(nr, nz, nx))
    return PlanarReleaseHistory(
        time_edges, surface.x_edges_m, surface.z_edges_m, temperatures, np.stack(coverage_edges),
        **{name: np.stack(values) for name, values in delivered.items()},
        reservoirs=pools, surface_attributes=surface.attributes, parameter_provenance=parameter_provenance,
    )
