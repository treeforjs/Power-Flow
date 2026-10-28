"""Weighted thermal boundary packets for separately owned downstream solvers.

Only emission is represented. No particle is advanced, collided or assigned a
new internal state. All masses, states, temperatures, reservoir parameters and
uniform variates are explicit inputs; no random-number generator is used.

The resolved and matched-uniform sources share exactly the same packets,
weights, emission times, internal states and momenta. Only their width positions
differ. A uniform source therefore retains the original source-temperature
mixture instead of being rethermalized at an averaged wall temperature.

Allocation is weighted stratification: first allocate one packet to every
positive (reservoir,width) component, then allocate remaining packets using
inverse-CDF component draws. Each component's released inventory is divided by
its allocated packet count. Weights can differ; this is NOT an equal-weight
collision-solver input. The downstream owner must support the weighted contract.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math
from numbers import Integral
from types import MappingProxyType

import numpy as np

from .constants import KB
from .surface_release_history import SurfaceReleaseHistory
from .thermal_release import SurfaceInventory, _array, _log_hazard, _product, _ratio, _real, release_step


_LOG_EPS = math.log(np.finfo(float).eps)
_LOG_MAX = math.log(np.finfo(float).max)
_C_M_S = 299_792_458.0
SCHEMA = "mhdlab-thermal-boundary-packets-v1"


class NoThermalEmissionError(ValueError):
    """The requested interval contains no representable positive source."""


def _immutable(value, dtype=None) -> np.ndarray:
    result = np.ascontiguousarray(value, dtype=dtype)
    return np.frombuffer(result.tobytes(), dtype=result.dtype).reshape(result.shape)


@dataclass(frozen=True, slots=True, eq=False)
class ThermalEmissionPackets:
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
    surface_attributes: Mapping[str, str]
    parameter_provenance: str
    schema: str = SCHEMA
    allocation_model: str = "one_per_positive_component_then_inverse_cdf_refinement"
    time_model: str = "component_conditional_truncated_exponential_at_frozen_temperature"
    velocity_model: str = "nonrelativistic_outward_flux_weighted_maxwellian"


def _conditional_time_offset(log_hazard: float, dt: float, quantile: float) -> float:
    """Invert the released-particle CDF, with no premature small-rate cutoff."""
    if log_hazard < _LOG_EPS:
        # The normalized conditional law differs from uniform by O(h), below
        # float64 relative precision here. The number release remains exact in
        # the source core and is never discarded by this conditional limit.
        return dt * quantile
    log_rate = log_hazard - math.log(dt)
    if log_hazard > _LOG_MAX:
        log_released_fraction = 0.0
    else:
        hazard = math.exp(log_hazard)
        log_released_fraction = math.log(-math.expm1(-hazard))
    log_cdf_argument = math.log(quantile) + log_released_fraction
    if log_cdf_argument < _LOG_EPS:
        log_elapsed_hazard = log_cdf_argument
    else:
        # The open quantile guarantees the argument is below one. Computing
        # its product directly is more accurate near one than exp(log(u)).
        fraction = math.exp(log_released_fraction)
        elapsed_hazard = -math.log1p(-quantile * fraction)
        log_elapsed_hazard = math.log(elapsed_hazard)
    elapsed = math.exp(log_elapsed_hazard - log_rate)
    if elapsed > dt:
        if elapsed > dt * (1.0 + 64 * np.finfo(float).eps):
            raise ValueError("emission offset is not representable inside its source interval")
        elapsed = dt  # only the documented last-bit endpoint rounding
    return elapsed


def sample_emission_packets(
    history: SurfaceReleaseHistory,
    state_masses_kg: Mapping[str, float],
    interval_index: int,
    uniforms,
    *,
    outward_normal_sign: int,
    current_length_m: float,
) -> ThermalEmissionPackets:
    """Construct paired resolved/uniform packets for one validated interval.

    ``uniforms`` is an n-by-6 array, strictly inside (0,1) in every column:
    component refinement, emission-time CDF, width position, normal-speed CDF,
    tangential Gaussian radius and tangential Gaussian angle. Its first
    ``required_minimum_n`` component draws are unused because every positive
    source component receives one mandatory packet, in reservoir/width order.

    Components have physical particle counts released_m2 * width * declared
    current_length_m. Dividing by their individual packet counts conserves each
    component. These unequal weights are not suitable for a downstream module
    that assumes equal-weight macroparticles without a validated adapter.

    Coordinates follow x=width,y=gap normal,z=current. The normal sign is
    required explicitly; velocities contain all three components. Absolute time
    and positions are best float64 values and may round to an interval/bin edge.
    Retain time_start_s + time_offset_s separately when sub-time precision is
    important. No additional precision is implied by ``time_s``.
    Sampled thermal energy is a weighted Monte Carlo/quadrature estimate of
    the expected 2*kB*T moment. Both source hypotheses share that same estimate;
    the returned component diagnostics retain its difference from the source
    budget. Velocities are never rescaled to force an exact thermal-energy sum.
    """
    if not isinstance(history, SurfaceReleaseHistory):
        raise ValueError("history must be a validated SurfaceReleaseHistory")
    if (isinstance(interval_index, (bool, np.bool_)) or not isinstance(interval_index, Integral)
            or not 0 <= interval_index < len(history.time_edges_s) - 1):
        raise ValueError("interval_index must identify one source interval")
    index = int(interval_index)
    if (isinstance(outward_normal_sign, (bool, np.bool_)) or not isinstance(outward_normal_sign, Integral)
            or outward_normal_sign not in (-1, 1)):
        raise ValueError("outward_normal_sign must be explicitly -1 or +1")
    sign = int(outward_normal_sign)
    length = _real(current_length_m, "current_length_m", positive=True)
    draws = _array(uniforms, "uniforms", 2)
    if draws.shape[1] != 6 or np.any(draws <= 0) or np.any(draws >= 1):
        raise ValueError("uniforms must have shape (n,6), strictly inside (0,1)")
    if not isinstance(state_masses_kg, Mapping):
        raise ValueError("state_masses_kg must map exact source state IDs to positive masses")
    masses_by_state = {}
    for pool in history.reservoirs:
        if pool.state_id not in state_masses_kg:
            raise ValueError(f"missing exact-state mass for {pool.state_id!r}; species fallback is forbidden")
        masses_by_state[pool.state_id] = _real(state_masses_kg[pool.state_id], "state mass_kg", positive=True)
    start, stop = (float(item) for item in history.time_edges_s[index:index + 2])
    dt = stop - start
    source = history.released_m2[index]
    temperatures = history.temperature_k[index]
    # The stored interval number/energy moments alone cannot justify a temporal
    # law. Confirm their declared frozen-temperature first-order provenance.
    expected = release_step(SurfaceInventory(history.coverage_m2[index]), history.reservoirs, temperatures, dt)
    tolerance = 128 * np.finfo(float).eps
    if (not np.allclose(source, expected.released_m2, rtol=tolerance, atol=0)
            or not np.allclose(history.coverage_m2[index + 1], expected.trial_inventory.coverage_m2,
                               rtol=tolerance, atol=0)):
        raise ValueError("source interval is inconsistent with the frozen-temperature first-order emission law")
    widths = np.diff(history.x_edges_m)
    components = []
    component_particles = []
    for pool_index, width_index in zip(*np.nonzero(source > 0)):
        particles = _product(float(source[pool_index, width_index]), float(widths[width_index]), length)
        if particles == 0:
            raise ValueError("positive component inventory is below the representable particle-weight range")
        components.append((int(pool_index), int(width_index)))
        component_particles.append(particles)
    if not components:
        raise NoThermalEmissionError("source interval has no positive release; no packets may be generated")
    required = len(components)
    packet_number = len(draws)
    if packet_number < required:
        raise ValueError(f"at least {required} packets are required to represent every positive source component")
    try:
        total_particles = math.fsum(component_particles)
    except OverflowError as exc:
        raise ValueError("total released particle count exceeds finite float64 range") from exc
    probabilities = np.array([_ratio(value, total_particles) for value in component_particles])
    cdf = np.cumsum(probabilities)
    cdf[-1] = 1.0  # normalize only accumulated floating-point endpoint roundoff
    assigned = np.empty(packet_number, dtype=np.int64)
    assigned[:required] = np.arange(required)
    assigned[required:] = np.searchsorted(cdf, draws[required:, 0], side="right")
    counts = np.bincount(assigned, minlength=required)
    component_weights = np.array([_ratio(value, float(count))
                                  for value, count in zip(component_particles, counts)])
    if np.any(component_weights <= 0):
        raise ValueError("positive component packet weight is not representable")
    reservoir_index = np.array([components[component][0] for component in assigned], dtype=np.int64)
    width_index = np.array([components[component][1] for component in assigned], dtype=np.int64)
    packet_temperatures = temperatures[width_index]
    state_ids = tuple(history.reservoirs[pool].state_id for pool in reservoir_index)
    names = tuple(history.reservoirs[pool].name for pool in reservoir_index)
    masses = np.array([masses_by_state[state] for state in state_ids])
    weights = component_weights[assigned]
    offsets = np.empty(packet_number)
    velocities = np.empty((packet_number, 3))
    for packet, (pool_index, source_width) in enumerate(zip(reservoir_index, width_index)):
        pool = history.reservoirs[pool_index]
        temperature = float(temperatures[source_width])
        offsets[packet] = _conditional_time_offset(_log_hazard(pool, temperature, dt), dt, float(draws[packet, 1]))
        sigma = math.exp(0.5 * math.fsum((math.log(KB), math.log(temperature), -math.log(masses[packet]))))
        normal = sigma * math.sqrt(-2.0 * math.log1p(-float(draws[packet, 3])))
        radius = sigma * math.sqrt(-2.0 * math.log(float(draws[packet, 4])))
        angle = 2.0 * math.pi * float(draws[packet, 5])
        velocities[packet] = (radius * math.cos(angle), sign * normal, radius * math.sin(angle))
        if not np.all(np.isfinite(velocities[packet])) or math.hypot(*velocities[packet]) > 0.05 * _C_M_S:
            raise ValueError("thermal boundary sample exceeds the supported nonrelativistic speed limit 0.05c")
    momenta = velocities * masses[:, None]
    expected_energy = np.array([
        _product(float(history.translational_energy_j_m2[index, pool, width]), float(widths[width]), length)
        for pool, width in components
    ])
    sampled_energy = np.empty(required)
    for component in range(required):
        selected = np.flatnonzero(assigned == component)
        try:
            sampled_energy[component] = math.fsum(
                _product(0.5, float(masses[packet]), math.hypot(*velocities[packet]),
                         math.hypot(*velocities[packet]), float(weights[packet]))
                for packet in selected
            )
        except OverflowError as exc:
            raise ValueError("sampled translational energy exceeds finite float64 range") from exc
    x0, x1 = (float(item) for item in history.x_edges_m[[0, -1]])
    resolved_x = history.x_edges_m[width_index] + widths[width_index] * draws[:, 2]
    uniform_x = x0 + (x1 - x0) * draws[:, 2]
    if not all(np.all(np.isfinite(value)) for value in (momenta, resolved_x, uniform_x, offsets)):
        raise ValueError("thermal boundary output exceeds finite float64 range")
    return ThermalEmissionPackets(
        interval_index=index, time_start_s=start, time_end_s=stop,
        time_offset_s=_immutable(offsets), time_s=_immutable(start + offsets),
        resolved_x_m=_immutable(resolved_x), uniform_x_m=_immutable(uniform_x),
        reservoir_index=_immutable(reservoir_index), source_width_index=_immutable(width_index),
        state_ids=state_ids, reservoir_names=names,
        temperature_k=_immutable(packet_temperatures), mass_kg=_immutable(masses),
        velocity_m_s=_immutable(velocities), momentum_kg_m_s=_immutable(momenta),
        particle_weight=_immutable(weights),
        internal_energy_ev=_immutable([history.reservoirs[pool].internal_energy_ev for pool in reservoir_index]),
        desorption_energy_ev=_immutable([history.reservoirs[pool].desorption_energy_ev for pool in reservoir_index]),
        component_reservoir_index=_immutable([pair[0] for pair in components], dtype=np.int64),
        component_width_index=_immutable([pair[1] for pair in components], dtype=np.int64),
        component_released_particles=_immutable(component_particles), component_probabilities=_immutable(probabilities),
        component_expected_translational_energy_j=_immutable(expected_energy),
        component_sampled_translational_energy_j=_immutable(sampled_energy),
        component_translational_energy_sampling_residual_j=_immutable(sampled_energy - expected_energy),
        packet_counts=_immutable(counts), required_minimum_n=required,
        total_released_particles=total_particles, current_length_m=length, outward_normal_sign=sign,
        unequal_macro_weights=bool(np.any(weights != weights[0])),
        surface_attributes=MappingProxyType(dict(history.surface_attributes)),
        parameter_provenance=history.parameter_provenance,
    )
