"""Temkin launches from resolved stationary current/width surface patches.

Local row samplers retain the coverage-dependent conditional time law and the
source-temperature mixture. The comparison changes width position only. Real
patch areas determine number and energy weights; no extra extrusion length or
whole-face temperature average is introduced.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
import math
from numbers import Integral
from types import MappingProxyType
from typing import Any

import numpy as np

from .flash_surface import validate_surface_attributes
from .planar_release import MFEM_SURFACE_SCHEMA, _source_attributes
from .planar_boundary import _allocate_rows, _coordinate, _immutable, _normal_sign, _validate_extracted_face
from .temkin_boundary import sample_temkin_emission_packets
from .temkin_emission_time import NoTemkinEmissionError, TIME_IDENTITY
from .temkin_history import TemkinReleaseHistory, STEP_FIELDS, DIAGNOSTIC_FIELDS
from .temkin_planar import TemkinPlanarReleaseHistory, SCHEMA as SOURCE_SCHEMA
from .temkin_release import LAW_IDENTITY
from .thermal_release import _array, _product, _real


SCHEMA = "mhdlab-temkin-planar-boundary-packets-v1"
ROW_ADAPTER_IDENTITY = "temkin_resolved_current_row_sampling_v1"


def _selected_history(history, index):
    arguments = {field.name: getattr(history, field.name) for field in fields(TemkinPlanarReleaseHistory)}
    arguments.update(time_edges_s=history.time_edges_s[index:index + 2],
                     temperature_k=history.temperature_k[index:index + 1],
                     coverage_m2=history.coverage_m2[index:index + 2])
    arguments.update({name: getattr(history, name)[index:index + 1] for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS)})
    return TemkinPlanarReleaseHistory(**arguments)


def _row_history(history, current, source_interval):
    attributes = dict(history.surface_attributes)
    attributes["schema"] = attributes["original_source_schema"]
    attributes = (dict(_source_attributes(attributes)) if attributes["schema"] == MFEM_SURFACE_SCHEMA
                  else validate_surface_attributes(attributes))
    attributes["temkin_planar_row_adapter"] = MappingProxyType({
        "identity": ROW_ADAPTER_IDENTITY, "parent_history_schema": SOURCE_SCHEMA,
        "current_index": current, "source_interval_index": source_interval,
        "current_interval_m": tuple(float(value) for value in history.z_edges_m[current:current + 2]),
        "qualification": "physical width and local temperatures of one resolved piecewise-constant current row; no whole-face uniformity or reconstruction of a current-averaged trace",
        "parent_surface_attributes_reference": "TemkinPlanarEmissionPackets.surface_attributes",
        "parent_source_schema": history.surface_attributes["schema"],
        "parent_run_id": history.surface_attributes["run_id"],
        "parent_source_reference": history.surface_attributes["source_reference"],
    })
    return TemkinReleaseHistory(
        time_edges_s=history.time_edges_s, x_edges_m=history.x_edges_m,
        temperature_k=history.temperature_k[:, current], coverage_m2=history.coverage_m2[:, :, current],
        **{name: getattr(history, name)[:, :, current] for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS)},
        reservoirs=history.reservoirs, surface_attributes=attributes,
        parameter_provenance=history.parameter_provenance, law_identity=history.law_identity,
        time_model=history.time_model,
    )


@dataclass(frozen=True, slots=True, eq=False)
class TemkinPlanarEmissionPackets:
    interval_index: int
    time_start_s: float
    time_end_s: float
    time_offset_s: np.ndarray
    time_s: np.ndarray
    resolved_x_m: np.ndarray
    uniform_x_m: np.ndarray
    resolved_position_m: np.ndarray
    uniform_width_position_m: np.ndarray
    source_current_offset_m: np.ndarray
    x_edges_m: np.ndarray
    z_edges_m: np.ndarray
    current_interval_m: tuple[float, float]
    current_length_m: float
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
    source_time_model: str
    coverage_model: str
    schema: str = SCHEMA
    law_identity: str = LAW_IDENTITY
    time_identity: str = TIME_IDENTITY
    time_model: str = TIME_IDENTITY
    source_history_schema: str = SOURCE_SCHEMA
    coordinate_system: str = "x=foil_width;y=gap_normal;z=current"
    comparison_model: str = "width_only_redistribution_preserving_each_packet_current_position_time_state_weight_and_momentum"
    allocation_model: str = "one_per_positive_component_then_largest_remainder_current_rows_and_row_inverse_cdf_refinement"
    velocity_model: str = "nonrelativistic_outward_flux_weighted_maxwellian"


def sample_temkin_planar_emission_packets(
    history: TemkinPlanarReleaseHistory, state_masses_kg: Mapping[str, float], interval_index: int,
    uniforms, *, face_normal_coordinate_m: float, outward_normal_sign: int,
) -> TemkinPlanarEmissionPackets:
    """Create paired 3D launches using explicit n-by-7 open variates.

    The first six columns feed the law-specific row sampler unchanged; column
    seven samples position within its current patch. One packet is seeded per
    positive reservoir/z/x component. Remaining row counts use exact-rational
    largest remainder, with current-index tie breaking. Unequal particle weights
    preserve component inventories and require a weighted downstream interface.

    Positions can round onto cell edges; current and time offsets are retained
    independently. Reported expected and sampled translational energies remain
    separate. No wall motion, temperature averaging or velocity adjustment is
    inferred from the supplied face coordinate.
    """
    if not isinstance(history, TemkinPlanarReleaseHistory):
        raise ValueError("history must be a distinct validated TemkinPlanarReleaseHistory")
    if (isinstance(interval_index, (bool, np.bool_)) or not isinstance(interval_index, Integral)
            or not 0 <= interval_index < history.time_edges_s.size - 1):
        raise ValueError("interval_index must identify one Temkin planar interval")
    index = int(interval_index)
    selected = _selected_history(history, index)
    coordinate, sign = _coordinate(face_normal_coordinate_m, "face_normal_coordinate_m"), _normal_sign(outward_normal_sign)
    _validate_extracted_face(selected, coordinate, sign)
    if np.ma.isMaskedArray(uniforms):
        raise ValueError("uniforms cannot contain masked variates")
    draws = _array(uniforms, "uniforms", 2)
    if draws.shape[1] != 7 or np.any(draws <= 0) or np.any(draws >= 1):
        raise ValueError("uniforms must have shape (n,7), strictly inside (0,1)")
    if not isinstance(state_masses_kg, Mapping):
        raise ValueError("state_masses_kg must map exact source states to positive masses")
    masses = {}
    for pool in selected.reservoirs:
        if pool.state_id not in state_masses_kg:
            raise ValueError(f"missing exact-state mass for {pool.state_id!r}; species fallback is forbidden")
        masses[pool.state_id] = _real(state_masses_kg[pool.state_id], "state mass_kg", positive=True)
    source = selected.released_m2[0]
    widths, current_widths = np.diff(selected.x_edges_m), np.diff(selected.z_edges_m)
    minimum = np.count_nonzero(source > 0, axis=(0, 2))
    if not np.any(minimum):
        raise NoTemkinEmissionError("selected planar interval has no positive Temkin release")
    populations = []
    for current in range(len(current_widths)):
        particles = [_product(float(source[pool, current, width]), float(widths[width]), float(current_widths[current]))
                     for pool, width in zip(*np.nonzero(source[:, current] > 0))]
        if any(value == 0 for value in particles):
            raise ValueError("positive local Temkin component is below representable particle-weight range")
        try:
            populations.append(math.fsum(particles))
        except OverflowError as error:
            raise ValueError("Temkin current-row population exceeds finite range") from error
    counts = _allocate_rows(minimum, populations, len(draws))
    rows, row_indices, row_attributes = [], [], []
    positions, uniform_positions, current_offsets = [], [], []
    start = 0
    for current, count in enumerate(counts):
        if count == 0:
            continue
        stop = start + int(count)
        row = _row_history(selected, current, index)
        packets = sample_temkin_emission_packets(row, masses, 0, draws[start:stop, :6],
                                                 outward_normal_sign=sign, current_length_m=float(current_widths[current]))
        offset = current_widths[current] * draws[start:stop, 6]
        z = selected.z_edges_m[current] + offset
        positions.append(np.column_stack((packets.resolved_x_m, np.full(int(count), coordinate), z)))
        uniform_positions.append(np.column_stack((packets.uniform_x_m, np.full(int(count), coordinate), z)))
        current_offsets.append(offset)
        rows.append(packets)
        row_indices.append(current)
        row_attributes.append(packets.surface_attributes)
        start = stop
    def concatenate(name):
        return _immutable(np.concatenate([getattr(row, name) for row in rows], axis=0))
    resolved, uniform = np.concatenate(positions), np.concatenate(uniform_positions)
    if not np.all(np.isfinite(resolved)) or not np.all(np.isfinite(uniform)):
        raise ValueError("Temkin planar launch coordinates exceed finite range")
    component_particles = concatenate("component_released_particles")
    expected_energy, sampled_energy = (concatenate(name) for name in
                                       ("component_expected_translational_energy_j", "component_sampled_translational_energy_j"))
    try:
        total, expected_total, sampled_total = (math.fsum(float(value) for value in values)
                                               for values in (component_particles, expected_energy, sampled_energy))
    except OverflowError as error:
        raise ValueError("Temkin planar packet totals exceed finite range") from error
    return TemkinPlanarEmissionPackets(
        interval_index=index, time_start_s=rows[0].time_start_s, time_end_s=rows[0].time_end_s,
        **{name: concatenate(name) for name in ("time_offset_s", "time_s", "resolved_x_m", "uniform_x_m",
             "reservoir_index", "source_width_index", "temperature_k", "mass_kg", "velocity_m_s",
             "momentum_kg_m_s", "particle_weight", "internal_energy_ev", "desorption_energy_ev",
             "component_reservoir_index", "component_width_index", "component_translational_energy_sampling_residual_j", "packet_counts")},
        resolved_position_m=_immutable(resolved), uniform_width_position_m=_immutable(uniform),
        source_current_offset_m=_immutable(np.concatenate(current_offsets)),
        x_edges_m=selected.x_edges_m, z_edges_m=selected.z_edges_m,
        current_interval_m=tuple(float(value) for value in selected.z_edges_m[[0, -1]]),
        current_length_m=float(selected.z_edges_m[-1] - selected.z_edges_m[0]),
        source_current_index=_immutable(np.concatenate([np.full(len(row.particle_weight), current, dtype=np.int64)
                                                       for row, current in zip(rows, row_indices)])),
        component_current_index=_immutable(np.concatenate([np.full(len(row.packet_counts), current, dtype=np.int64)
                                                          for row, current in zip(rows, row_indices)])),
        state_ids=tuple(state for row in rows for state in row.state_ids),
        reservoir_names=tuple(name for row in rows for name in row.reservoir_names),
        component_released_particles=component_particles, component_probabilities=_immutable(component_particles / total),
        component_expected_translational_energy_j=expected_energy, component_sampled_translational_energy_j=sampled_energy,
        minimum_packets_per_current_row=_immutable(minimum), packets_per_current_row=_immutable(counts),
        released_particles_per_current_row=_immutable(populations), required_minimum_n=int(np.sum(minimum)),
        total_released_particles=total, expected_translational_energy_j=expected_total,
        sampled_translational_energy_j=sampled_total,
        translational_energy_sampling_residual_j=math.fsum((sampled_total, -expected_total)),
        face_normal_coordinate_m=coordinate, outward_normal_sign=sign,
        unequal_macro_weights=any(float(value) != float(rows[0].particle_weight[0]) for row in rows for value in row.particle_weight),
        surface_attributes=selected.surface_attributes, row_surface_attributes=tuple(row_attributes),
        row_current_indices=tuple(row_indices), parameter_provenance=selected.parameter_provenance,
        source_time_model=selected.time_model, coverage_model=selected.coverage_model,
    )
