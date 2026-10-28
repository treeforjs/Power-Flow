"""Weighted launch packets with the conditional Temkin emission-time law.

The two source hypotheses share every packet's state, time, thermal velocity,
momentum and weight. Only width position is redistributed. Local temperature
mixtures remain intact. This distinct record invokes no constant-hazard packet
adapter, field advance, collision model or radiative solver.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
import math
from numbers import Integral
from typing import Any

import numpy as np

from .constants import KB
from .temkin_history import TemkinReleaseHistory, SCHEMA as SOURCE_SCHEMA, STEP_FIELDS, DIAGNOSTIC_FIELDS
from .temkin_release import LAW_IDENTITY
from .temkin_emission_time import FrozenTemkinEmission, NoTemkinEmissionError, TIME_IDENTITY
from .thermal_release import _array, _product, _ratio, _real


SCHEMA = "mhdlab-temkin-boundary-packets-v1"
_C_M_S = 299_792_458.0


def _immutable(value, dtype=None):
    result = np.ascontiguousarray(value, dtype=dtype)
    return np.frombuffer(result.tobytes(), dtype=result.dtype).reshape(result.shape)


def _selected_interval(history, index):
    """Replay only the selected interval, retaining the actual Temkin type."""
    arguments = {field.name: getattr(history, field.name) for field in fields(TemkinReleaseHistory)}
    arguments.update(time_edges_s=history.time_edges_s[index:index + 2],
                     temperature_k=history.temperature_k[index:index + 1],
                     coverage_m2=history.coverage_m2[index:index + 2])
    arguments.update({name: getattr(history, name)[index:index + 1] for name in (*STEP_FIELDS, *DIAGNOSTIC_FIELDS)})
    return TemkinReleaseHistory(**arguments)


@dataclass(frozen=True, slots=True, eq=False)
class TemkinEmissionPackets:
    interval_index: int
    time_start_s: float
    time_end_s: float
    time_offset_s: np.ndarray
    time_s: np.ndarray
    resolved_x_m: np.ndarray
    uniform_x_m: np.ndarray
    reservoir_index: np.ndarray
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
    component_width_index: np.ndarray
    component_released_particles: np.ndarray
    component_probabilities: np.ndarray
    component_expected_translational_energy_j: np.ndarray
    component_sampled_translational_energy_j: np.ndarray
    component_translational_energy_sampling_residual_j: np.ndarray
    packet_counts: np.ndarray
    required_minimum_n: int
    total_released_particles: float
    current_length_m: float
    outward_normal_sign: int
    unequal_macro_weights: bool
    surface_attributes: Mapping[str, Any]
    parameter_provenance: str
    source_time_model: str
    schema: str = SCHEMA
    law_identity: str = LAW_IDENTITY
    time_identity: str = TIME_IDENTITY
    time_model: str = TIME_IDENTITY
    source_history_schema: str = SOURCE_SCHEMA
    allocation_model: str = "one_per_positive_component_then_inverse_cdf_refinement"
    velocity_model: str = "nonrelativistic_outward_flux_weighted_maxwellian"


def sample_temkin_emission_packets(
    history: TemkinReleaseHistory,
    state_masses_kg: Mapping[str, float],
    interval_index: int,
    uniforms,
    *,
    outward_normal_sign: int,
    current_length_m: float,
) -> TemkinEmissionPackets:
    """Sample a law-qualified frozen-temperature Temkin source interval.

    The n-by-6 strictly open variates select component refinement, conditional
    emission time, width position, outward normal speed, tangential Gaussian
    radius, and tangential angle. Every positive reservoir/width component gets
    one mandatory packet before component-CDF refinement. Unequal weights
    conserve every component independently and require downstream support.

    Conditional times use the integrated coverage-dependent release law, not
    a fixed initial hazard. Local offsets are retained separately from rounded
    absolute times. Sampled thermal energy is compared with the expected 2kT
    source moment; no velocity rescaling forces those finite sums to agree.
    Tiny component probabilities may round to zero while mandatory packets
    retain their positive inventories. No random-number generator is used.
    """
    if not isinstance(history, TemkinReleaseHistory):
        raise ValueError("history must be a distinct validated TemkinReleaseHistory")
    if (isinstance(interval_index, (bool, np.bool_)) or not isinstance(interval_index, Integral)
            or not 0 <= interval_index < len(history.time_edges_s) - 1):
        raise ValueError("interval_index must identify one Temkin source interval")
    index = int(interval_index)
    if (isinstance(outward_normal_sign, (bool, np.bool_)) or not isinstance(outward_normal_sign, Integral)
            or outward_normal_sign not in (-1, 1)):
        raise ValueError("outward_normal_sign must be explicitly -1 or +1")
    sign = int(outward_normal_sign)
    length = _real(current_length_m, "current_length_m", positive=True)
    if np.ma.isMaskedArray(uniforms):
        raise ValueError("uniforms cannot contain masked variates")
    draws = _array(uniforms, "uniforms", 2)
    if draws.shape[1] != 6 or np.any(draws <= 0) or np.any(draws >= 1):
        raise ValueError("uniforms must have shape (n,6), strictly inside (0,1)")
    selected_history = _selected_interval(history, index)
    if not isinstance(state_masses_kg, Mapping):
        raise ValueError("state_masses_kg must map exact source states to positive masses")
    masses_by_state = {}
    for pool in selected_history.reservoirs:
        if pool.state_id not in state_masses_kg:
            raise ValueError(f"missing exact-state mass for {pool.state_id!r}; species fallback is forbidden")
        masses_by_state[pool.state_id] = _real(state_masses_kg[pool.state_id], "state mass_kg", positive=True)
    start, stop = (float(value) for value in selected_history.time_edges_s)
    dt = stop - start
    source, temperature = selected_history.released_m2[0], selected_history.temperature_k[0]
    widths = np.diff(selected_history.x_edges_m)
    components, particles = [], []
    for pool, width in zip(*np.nonzero(source > 0)):
        count = _product(float(source[pool, width]), float(widths[width]), length)
        if count == 0:
            raise ValueError("positive Temkin component is below representable particle-weight range")
        components.append((int(pool), int(width)))
        particles.append(count)
    if not components:
        raise NoTemkinEmissionError("selected Temkin interval has no positive release")
    required, number = len(components), len(draws)
    if number < required:
        raise ValueError(f"at least {required} packets are required to represent every positive Temkin component")
    try:
        total = math.fsum(particles)
    except OverflowError as error:
        raise ValueError("Temkin released particle total exceeds finite numeric range") from error
    probabilities = np.array([_ratio(value, total) for value in particles])
    cumulative = np.minimum(np.cumsum(probabilities), 1.)
    cumulative[-1] = 1
    assigned = np.empty(number, dtype=np.int64)
    assigned[:required] = np.arange(required)
    assigned[required:] = np.searchsorted(cumulative, draws[required:, 0], side="right")
    counts = np.bincount(assigned, minlength=required)
    component_weights = np.array([_ratio(value, float(count)) for value, count in zip(particles, counts)])
    if np.any(component_weights <= 0):
        raise ValueError("positive Temkin packet weight is not representable")
    reservoir_index = np.array([components[component][0] for component in assigned], dtype=np.int64)
    width_index = np.array([components[component][1] for component in assigned], dtype=np.int64)
    state_ids = tuple(selected_history.reservoirs[pool].state_id for pool in reservoir_index)
    names = tuple(selected_history.reservoirs[pool].name for pool in reservoir_index)
    masses = np.array([masses_by_state[state] for state in state_ids])
    weights, temperatures = component_weights[assigned], temperature[width_index]
    offsets = np.empty(number)
    for component, (pool_index, width) in enumerate(components):
        emission = FrozenTemkinEmission(selected_history.reservoirs[pool_index],
                                        float(selected_history.coverage_m2[0, pool_index, width]),
                                        float(temperature[width]), dt)
        if emission.law_identity != LAW_IDENTITY or emission.time_identity != TIME_IDENTITY:
            raise ValueError("conditional emission law/time identity does not match the Temkin source")
        packets = np.flatnonzero(assigned == component)
        offsets[packets] = emission.quantile(draws[packets, 1])
    velocities = np.empty((number, 3))
    for packet in range(number):
        sigma = math.exp(.5 * math.fsum((math.log(KB), math.log(float(temperatures[packet])), -math.log(float(masses[packet])))))
        normal = sigma * math.sqrt(-2 * math.log1p(-float(draws[packet, 3])))
        radius = sigma * math.sqrt(-2 * math.log(float(draws[packet, 4])))
        angle = 2 * math.pi * float(draws[packet, 5])
        velocities[packet] = (radius * math.cos(angle), sign * normal, radius * math.sin(angle))
        if not np.all(np.isfinite(velocities[packet])) or math.hypot(*velocities[packet]) > .05 * _C_M_S:
            raise ValueError("Temkin thermal sample exceeds the nonrelativistic speed limit 0.05c")
    momenta = masses[:, None] * velocities
    expected_energy = np.array([_product(float(selected_history.translational_energy_j_m2[0, pool, width]),
                                          float(widths[width]), length) for pool, width in components])
    sampled_energy = np.empty(required)
    for component in range(required):
        try:
            sampled_energy[component] = math.fsum(
                _product(.5, float(masses[packet]), math.hypot(*velocities[packet]),
                         math.hypot(*velocities[packet]), float(weights[packet]))
                for packet in np.flatnonzero(assigned == component))
        except OverflowError as error:
            raise ValueError("sampled Temkin translational energy exceeds finite numeric range") from error
    x0, x1 = (float(value) for value in selected_history.x_edges_m[[0, -1]])
    resolved_x = selected_history.x_edges_m[width_index] + widths[width_index] * draws[:, 2]
    uniform_x = x0 + (x1 - x0) * draws[:, 2]
    if not all(np.all(np.isfinite(value)) for value in (momenta, offsets, resolved_x, uniform_x)):
        raise ValueError("Temkin packet output exceeds finite numeric range")
    return TemkinEmissionPackets(
        interval_index=index, time_start_s=start, time_end_s=stop,
        time_offset_s=_immutable(offsets), time_s=_immutable(start + offsets),
        resolved_x_m=_immutable(resolved_x), uniform_x_m=_immutable(uniform_x),
        reservoir_index=_immutable(reservoir_index), source_width_index=_immutable(width_index),
        state_ids=state_ids, reservoir_names=names, temperature_k=_immutable(temperatures),
        mass_kg=_immutable(masses), velocity_m_s=_immutable(velocities), momentum_kg_m_s=_immutable(momenta),
        particle_weight=_immutable(weights),
        internal_energy_ev=_immutable([selected_history.reservoirs[pool].internal_energy_ev for pool in reservoir_index]),
        desorption_energy_ev=_immutable([selected_history.reservoirs[pool].desorption_energy_ev for pool in reservoir_index]),
        component_reservoir_index=_immutable([pool for pool, _ in components], dtype=np.int64),
        component_width_index=_immutable([width for _, width in components], dtype=np.int64),
        component_released_particles=_immutable(particles), component_probabilities=_immutable(probabilities),
        component_expected_translational_energy_j=_immutable(expected_energy),
        component_sampled_translational_energy_j=_immutable(sampled_energy),
        component_translational_energy_sampling_residual_j=_immutable(sampled_energy - expected_energy),
        packet_counts=_immutable(counts), required_minimum_n=required, total_released_particles=total,
        current_length_m=length, outward_normal_sign=sign, unequal_macro_weights=bool(np.any(weights != weights[0])),
        surface_attributes=selected_history.surface_attributes, parameter_provenance=selected_history.parameter_provenance,
        source_time_model=selected_history.time_model,
    )
