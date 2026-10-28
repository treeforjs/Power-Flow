"""Paired thermal launch packets from resolved stationary surface patches.

Each current row is passed to the existing thermal boundary sampler. The
comparison source redistributes width positions only: current position, normal
position, time, internal state, weight and momentum are identical in each pair.
The current-row allocation is deterministic largest remainder after seeding
every positive reservoir/current/width component. No field or particle advance
is performed, and unequal packet weights require downstream support.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
import math
from numbers import Integral, Real
from types import MappingProxyType
from typing import Any

import numpy as np

from .planar_release import PLANAR_RELEASE_SCHEMA, PlanarReleaseHistory
from .surface_release_history import SurfaceReleaseHistory
from .thermal_boundary import NoThermalEmissionError, sample_emission_packets
from .thermal_release import _array, _product, _real


SCHEMA = "mhdlab-planar-boundary-packets-v1"
_ROW_FIELDS = (
    "released_m2", "translational_energy_j_m2", "internal_energy_j_m2",
    "binding_energy_j_m2", "total_wall_energy_debit_j_m2", "inventory_balance_residual_m2",
)


def _immutable(value, dtype=None):
    result = np.ascontiguousarray(value, dtype=dtype)
    return np.frombuffer(result.tobytes(), dtype=result.dtype).reshape(result.shape)


def _coordinate(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be an explicit finite real coordinate")
    try:
        value = float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be finite") from error
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _normal_sign(value):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or value not in (-1, 1)):
        raise ValueError("outward_normal_sign must be explicitly -1 or +1")
    return int(value)


def _validate_extracted_face(history, coordinate, sign):
    """Check geometry retained by the raw patch exporter when it is present."""
    if "mfem_face" in history.surface_attributes:
        face = history.surface_attributes["mfem_face"]
        if (not isinstance(face, Mapping)
                or _coordinate(face.get("normal_position_m"), "MFEM face coordinate") != coordinate
                or _normal_sign(face.get("outward_normal_sign")) != sign):
            raise ValueError("launch face conflicts with the preserved MFEM material face")
    if "resolved_patch_extraction" not in history.surface_attributes:
        return
    extraction = history.surface_attributes["resolved_patch_extraction"]
    if not isinstance(extraction, Mapping):
        raise ValueError("resolved_patch_extraction must preserve its metadata mapping")
    manifests = extraction.get("patch_extraction_manifests")
    if (not isinstance(manifests, Sequence) or isinstance(manifests, (str, bytes))
            or len(manifests) != history.z_edges_m.size - 1):
        raise ValueError("raw patch metadata must contain one manifest per current interval")
    for index, manifest in enumerate(manifests):
        if not isinstance(manifest, Mapping) or not isinstance(manifest.get("sampling"), Mapping):
            raise ValueError("raw patch manifest must retain its sampling geometry")
        sampling = manifest["sampling"]
        raw_coordinate = _coordinate(sampling.get("face_normal_coordinate_m"), "preserved face coordinate")
        raw_sign = _normal_sign(sampling.get("outward_sign"))
        if raw_coordinate != coordinate or raw_sign != sign:
            raise ValueError("launch face coordinate/sign conflicts with preserved raw extraction geometry")
        interval = sampling.get("current_interval_m")
        if (not isinstance(interval, Sequence) or isinstance(interval, (str, bytes)) or len(interval) != 2
                or tuple(_coordinate(value, "preserved current interval") for value in interval)
                != tuple(history.z_edges_m[index:index + 2])):
            raise ValueError("raw patch geometry conflicts with the resolved current interval")
        depth = sampling.get("material_sample_depth_m")
        if _coordinate(depth, "preserved material depth") != history.surface_attributes["material_sample_depth_m"]:
            raise ValueError("raw patch geometry conflicts with the preserved sampling depth")


def _row_history(history, current_index):
    """A width-only sampling adapter for one already resolved current patch."""
    lo, hi = (float(value) for value in history.z_edges_m[current_index:current_index + 2])
    attributes = dict(history.surface_attributes)
    attributes["schema"] = attributes["original_source_schema"]
    attributes["planar_row_adapter"] = MappingProxyType({
        "parent_history_schema": PLANAR_RELEASE_SCHEMA,
        "current_index": current_index, "current_interval_m": (lo, hi),
        "representation": "one already resolved piecewise-constant current row; no whole-current trace reconstruction",
        "parent_surface_attributes": history.surface_attributes,
    })
    return SurfaceReleaseHistory(
        time_edges_s=history.time_edges_s, x_edges_m=history.x_edges_m,
        temperature_k=history.temperature_k[:, current_index, :],
        coverage_m2=history.coverage_m2[:, :, current_index, :],
        **{name: getattr(history, name)[:, :, current_index, :] for name in _ROW_FIELDS},
        reservoirs=history.reservoirs, surface_attributes=attributes,
        parameter_provenance=history.parameter_provenance,
    )


def _allocate_rows(minimum, populations, packet_number):
    """Allocate extras by exact rational largest remainder, ties by z index."""
    result = np.array(minimum, dtype=np.int64, copy=True)
    remaining = packet_number - int(np.sum(result))
    if remaining < 0:
        raise ValueError(f"at least {int(np.sum(result))} packets are required to represent every positive source component")
    if remaining == 0:
        return result
    active = [index for index, value in enumerate(populations) if value > 0]
    weights = {index: Fraction.from_float(float(populations[index])) for index in active}
    normalizer = sum(weights.values(), Fraction())
    quotas = {index: remaining * weight / normalizer for index, weight in weights.items()}
    extras = {index: quota.numerator // quota.denominator for index, quota in quotas.items()}
    leftovers = remaining - sum(extras.values())
    order = sorted(active, key=lambda index: (-(quotas[index] - extras[index]), index))
    for index in order[:leftovers]:
        extras[index] += 1
    for index, count in extras.items():
        result[index] += count
    return result


@dataclass(frozen=True, slots=True, eq=False)
class PlanarEmissionPackets:
    interval_index: int
    time_start_s: float
    time_end_s: float
    time_offset_s: np.ndarray
    time_s: np.ndarray
    resolved_position_m: np.ndarray
    uniform_width_position_m: np.ndarray
    source_current_offset_m: np.ndarray
    reservoir_index: np.ndarray
    source_current_index: np.ndarray
    source_width_index: np.ndarray
    state_ids: tuple[str, ...]
    reservoir_names: tuple[str, ...]
    temperature_k: np.ndarray
    mass_kg: np.ndarray
    velocity_m_s: np.ndarray
    momentum_kg_m_s: np.ndarray
    particle_weight: np.ndarray
    internal_energy_ev: np.ndarray
    desorption_energy_ev: np.ndarray
    component_reservoir_index: np.ndarray
    component_current_index: np.ndarray
    component_width_index: np.ndarray
    component_released_particles: np.ndarray
    component_probabilities: np.ndarray
    component_expected_translational_energy_j: np.ndarray
    component_sampled_translational_energy_j: np.ndarray
    component_translational_energy_sampling_residual_j: np.ndarray
    packet_counts: np.ndarray
    minimum_packets_per_current_row: np.ndarray
    packets_per_current_row: np.ndarray
    released_particles_per_current_row: np.ndarray
    required_minimum_n: int
    total_released_particles: float
    expected_translational_energy_j: float
    sampled_translational_energy_j: float
    translational_energy_sampling_residual_j: float
    face_normal_coordinate_m: float
    outward_normal_sign: int
    unequal_macro_weights: bool
    surface_attributes: Mapping[str, Any]
    row_surface_attributes: tuple[Mapping[str, Any], ...]
    row_current_indices: tuple[int, ...]
    parameter_provenance: str
    schema: str = SCHEMA
    coordinate_system: str = "x=foil_width;y=gap_normal;z=current"
    comparison_model: str = "width_only_redistribution_preserving_each_packet_current_position_time_state_weight_and_momentum"
    allocation_model: str = "one_per_positive_component_then_largest_remainder_current_rows_and_row_inverse_cdf_refinement"
    time_model: str = "component_conditional_truncated_exponential_at_frozen_temperature"
    velocity_model: str = "nonrelativistic_outward_flux_weighted_maxwellian"


def sample_planar_emission_packets(
    history: PlanarReleaseHistory,
    state_masses_kg: Mapping[str, float],
    interval_index: int,
    uniforms,
    *,
    face_normal_coordinate_m: float,
    outward_normal_sign: int,
) -> PlanarEmissionPackets:
    """Return paired stationary-face launches, with no downstream integration.

    All n-by-7 variates must lie strictly inside (0,1). Rows of variates are
    assigned in increasing current-index order. The first six columns pass
    unchanged to ``sample_emission_packets``; column seven places the packet
    uniformly within that resolved current interval. Every positive local
    reservoir/current/width component receives at least one packet.

    Positions are best float64 coordinates and can round onto bin edges. The
    current offset and emission-time offset are retained separately. Thermal
    energy is sampled, not forced to its expected 2*kB*T moment; the component
    energy diagnostics are transferred without rescaling velocities. Weights
    generally differ and require a downstream weighted-particle interface.
    Component probabilities are diagnostics: a tiny probability can round to
    zero while its mandatory packet still retains the positive source weight.
    """
    if not isinstance(history, PlanarReleaseHistory):
        raise ValueError("history must be a validated PlanarReleaseHistory")
    if (isinstance(interval_index, (bool, np.bool_)) or not isinstance(interval_index, Integral)
            or not 0 <= interval_index < history.time_edges_s.size - 1):
        raise ValueError("interval_index must identify one source interval")
    index = int(interval_index)
    coordinate, sign = _coordinate(face_normal_coordinate_m, "face_normal_coordinate_m"), _normal_sign(outward_normal_sign)
    _validate_extracted_face(history, coordinate, sign)
    draws = _array(uniforms, "uniforms", 2)
    if draws.shape[1] != 7 or np.any(draws <= 0) or np.any(draws >= 1):
        raise ValueError("uniforms must have shape (n,7), strictly inside (0,1)")
    if not isinstance(state_masses_kg, Mapping):
        raise ValueError("state_masses_kg must map exact source states to positive masses")
    masses = {}
    for pool in history.reservoirs:
        if pool.state_id not in state_masses_kg:
            raise ValueError(f"missing exact-state mass for {pool.state_id!r}; species fallback is forbidden")
        masses[pool.state_id] = _real(state_masses_kg[pool.state_id], "state mass_kg", positive=True)
    source = history.released_m2[index]
    widths, current_widths = np.diff(history.x_edges_m), np.diff(history.z_edges_m)
    minimum = np.count_nonzero(source > 0, axis=(0, 2))
    if not np.any(minimum):
        raise NoThermalEmissionError("source interval has no positive release; no packets may be generated")
    populations = []
    for current in range(len(current_widths)):
        particles = [_product(float(source[pool, current, width]), float(widths[width]), float(current_widths[current]))
                     for pool, width in zip(*np.nonzero(source[:, current, :] > 0))]
        if any(value == 0 for value in particles):
            raise ValueError("positive local component is below the representable particle-weight range")
        try:
            populations.append(math.fsum(particles))
        except OverflowError as error:
            raise ValueError("current-row source exceeds finite particle-count range") from error
    counts = _allocate_rows(minimum, populations, len(draws))
    rows, row_indices, row_attributes, positions, uniform_positions, current_offsets = [], [], [], [], [], []
    start = 0
    for current, count in enumerate(counts):
        if count == 0:
            continue
        stop = start + int(count)
        row = _row_history(history, current)
        packets = sample_emission_packets(row, masses, index, draws[start:stop, :6],
                                          outward_normal_sign=sign, current_length_m=float(current_widths[current]))
        offset = current_widths[current] * draws[start:stop, 6]
        z = history.z_edges_m[current] + offset
        positions.append(np.column_stack((packets.resolved_x_m, np.full(int(count), coordinate), z)))
        uniform_positions.append(np.column_stack((packets.uniform_x_m, np.full(int(count), coordinate), z)))
        current_offsets.append(offset)
        rows.append(packets)
        row_indices.append(current)
        row_attributes.append(packets.surface_attributes)
        start = stop
    def concatenate(name):
        return _immutable(np.concatenate([getattr(row, name) for row in rows], axis=0))
    all_positions, all_uniform = np.concatenate(positions), np.concatenate(uniform_positions)
    if not np.all(np.isfinite(all_positions)) or not np.all(np.isfinite(all_uniform)):
        raise ValueError("launch positions exceed the finite coordinate range")
    component_particles = concatenate("component_released_particles")
    expected_energy = concatenate("component_expected_translational_energy_j")
    sampled_energy = concatenate("component_sampled_translational_energy_j")
    try:
        total = math.fsum(float(value) for value in component_particles)
        expected_total = math.fsum(float(value) for value in expected_energy)
        sampled_total = math.fsum(float(value) for value in sampled_energy)
    except OverflowError as error:
        raise ValueError("planar packet number/energy totals exceed finite numeric range") from error
    return PlanarEmissionPackets(
        interval_index=index, time_start_s=rows[0].time_start_s, time_end_s=rows[0].time_end_s,
        **{name: concatenate(name) for name in ("time_offset_s", "time_s", "reservoir_index", "source_width_index",
             "temperature_k", "mass_kg", "velocity_m_s", "momentum_kg_m_s", "particle_weight",
             "internal_energy_ev", "desorption_energy_ev", "component_reservoir_index", "component_width_index",
             "component_translational_energy_sampling_residual_j", "packet_counts")},
        resolved_position_m=_immutable(all_positions), uniform_width_position_m=_immutable(all_uniform),
        source_current_offset_m=_immutable(np.concatenate(current_offsets)),
        source_current_index=_immutable(np.concatenate([np.full(len(row.particle_weight), current, dtype=np.int64)
                                                       for row, current in zip(rows, row_indices)])),
        state_ids=tuple(state for row in rows for state in row.state_ids),
        reservoir_names=tuple(name for row in rows for name in row.reservoir_names),
        component_current_index=_immutable(np.concatenate([np.full(len(row.packet_counts), current, dtype=np.int64)
                                                          for row, current in zip(rows, row_indices)])),
        component_released_particles=component_particles,
        component_probabilities=_immutable([float(value) / total for value in component_particles]),
        component_expected_translational_energy_j=expected_energy,
        component_sampled_translational_energy_j=sampled_energy,
        minimum_packets_per_current_row=_immutable(minimum), packets_per_current_row=_immutable(counts),
        released_particles_per_current_row=_immutable(populations), required_minimum_n=int(np.sum(minimum)),
        total_released_particles=total, expected_translational_energy_j=expected_total,
        sampled_translational_energy_j=sampled_total,
        translational_energy_sampling_residual_j=math.fsum((sampled_total, -expected_total)),
        face_normal_coordinate_m=coordinate, outward_normal_sign=sign,
        unequal_macro_weights=any(float(value) != float(rows[0].particle_weight[0]) for row in rows for value in row.particle_weight),
        surface_attributes=history.surface_attributes, row_surface_attributes=tuple(row_attributes),
        row_current_indices=tuple(row_indices), parameter_provenance=history.parameter_provenance,
    )
