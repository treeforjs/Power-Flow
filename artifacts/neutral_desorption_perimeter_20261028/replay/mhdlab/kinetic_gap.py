"""Coupled prescribed-field, homogeneous-cell kinetic reference operator.

This small reference couples actual particle orbits to state-to-state collisions.
It is a verification target for a native solver port, not a self-consistent PIC
solver. All equal-weight particles share one collision volume; spatially varying
collision cells, walls, sources and charge-conserving field deposition are not
implemented. A reaction occurs at the reactants' mass-weighted position. This
cell collision approximation must not be used as subcell PIC source deposition.

Streaming then collisions is a first-order operator split. Rates are recomputed
after every Gillespie event in the frozen collision cell. Every channel requires
supplied joint energy-sharing/direction data through an explicit kernel callback.
There is no built-in water mechanism, isotropic-fragment fallback or Hall force.
"""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from functools import wraps
import math
from types import MappingProxyType
from typing import Callable, Mapping

import numpy as np

from .kinetic_orbits import advance_orbits, kinetic_energy_j
from .reaction_kinematics import react
from .state_kinetics import ELEMENTARY_CHARGE_C, ReactionChannel, StateCatalog


def _vector(value, name: str) -> tuple[float, float, float]:
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be real")
    array = np.asarray(value, dtype=float)
    if array.shape != (3,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite three-component vector")
    return tuple(float(item) for item in array)


@dataclass(frozen=True)
class Particle:
    state_id: str
    position_m: tuple[float, float, float]
    momentum_kg_m_s: tuple[float, float, float]
    weight: float = 1.0

    def __post_init__(self):
        if not isinstance(self.state_id, str) or not self.state_id.strip():
            raise ValueError("particle requires an explicit state_id")
        object.__setattr__(self, "position_m", _vector(self.position_m, "position_m"))
        object.__setattr__(self, "momentum_kg_m_s", _vector(self.momentum_kg_m_s, "momentum_kg_m_s"))
        if isinstance(self.weight, bool) or not np.isfinite(self.weight) or self.weight <= 0:
            raise ValueError("weight must be positive and finite")


@dataclass(frozen=True)
class GapState:
    particles: tuple[Particle, ...]
    time_s: float = 0.0

    def __post_init__(self):
        particles = tuple(self.particles)
        if any(not isinstance(p, Particle) for p in particles):
            raise ValueError("particles must contain Particle records")
        if not np.isfinite(self.time_s) or self.time_s < 0:
            raise ValueError("time_s must be finite and nonnegative")
        object.__setattr__(self, "particles", particles)


@dataclass(frozen=True)
class GapStep:
    state: GapState
    reaction_counts: Mapping[str, int]
    orbit_substeps: int
    field_kinetic_gain_j: float
    electric_displacement_work_j: float
    field_work_quadrature_residual_j: float
    reaction_energy_residual_j: float
    total_energy_residual_j: float
    reaction_momentum_residual_kg_m_s: tuple[float, float, float]


def total_energy_j(state: GapState, catalog: StateCatalog) -> float:
    """Nonrelativistic translational plus common-reference internal energy."""
    terms = []
    for p in state.particles:
        kinetic = _particle_kinetic_j(p, catalog)
        terms.extend((p.weight * kinetic, p.weight * catalog.energy_ev(p.state_id) * ELEMENTARY_CHARGE_C))
    if not all(math.isfinite(v) for v in terms):
        raise ValueError("weighted particle energy exceeds finite numeric range")
    try:
        return math.fsum(terms)
    except OverflowError as exc:
        raise ValueError("total particle energy exceeds finite numeric range") from exc


def _particle_kinetic_j(p: Particle, catalog: StateCatalog) -> float:
    return float(kinetic_energy_j(np.array([p.momentum_kg_m_s]),
                                 np.array([catalog.mass_kg(p.state_id)]), relativistic=False)[0])


def sample_thermal_emission(
    catalog: StateCatalog, state_ids: tuple[str, ...], position_m, outward_unit_normals,
    surface_temperature_k, *, rng: np.random.Generator, weight: float = 1.0,
    wall_velocity_m_s=None,
) -> tuple[Particle, ...]:
    """Sample flux-weighted thermal launches with explicitly supplied states.

    Each outgoing normal velocity is Rayleigh distributed; the two tangential
    components are Gaussian. Mean wall-frame translational energy is 2 kB T.
    Internal states are inputs and are not assigned a Boltzmann distribution.
    This samples velocities only: an upstream surface model must debit inventory
    and the consistent emitted/binding/internal-energy flux. Shape is (n,3) for
    positions/normals/wall velocities and (n,) for temperatures; no unit guessing.
    """
    if not isinstance(rng, np.random.Generator):
        raise ValueError("an explicit numpy random Generator is required")
    ids = tuple(state_ids)
    n = len(ids)
    if any(np.iscomplexobj(a) for a in (position_m, outward_unit_normals, surface_temperature_k)) or np.iscomplexobj(wall_velocity_m_s):
        raise ValueError("thermal launch arrays must be real")
    positions = np.asarray(position_m, dtype=float)
    normals = np.asarray(outward_unit_normals, dtype=float)
    temperatures = np.asarray(surface_temperature_k, dtype=float)
    wall = np.zeros((n, 3)) if wall_velocity_m_s is None else np.asarray(wall_velocity_m_s, dtype=float)
    if any(a.shape != (n, 3) or not np.all(np.isfinite(a)) for a in (positions, normals, wall)):
        raise ValueError("position, normal and wall velocity arrays must be finite (n,3)")
    if temperatures.shape != (n,) or not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0):
        raise ValueError("surface_temperature_k must be a positive finite (n,) array")
    if np.any(np.abs(np.linalg.norm(normals, axis=1) - 1.0) > 1e-12):
        raise ValueError("outward normals must be unit vectors")
    masses = np.array([catalog.mass_kg(s) for s in ids])
    sigma = np.sqrt(1.380649e-23 * temperatures / masses)
    gaussian = rng.normal(size=(n, 3)) * sigma[:, None]
    tangent = gaussian - np.sum(gaussian * normals, axis=1)[:, None] * normals
    normal_speed = rng.rayleigh(scale=sigma)
    velocities = wall + tangent + normal_speed[:, None] * normals
    return tuple(Particle(s, tuple(x), tuple(m * v), weight)
                 for s, x, m, v in zip(ids, positions, masses, velocities))


def cell_moments(state: GapState, catalog: StateCatalog, volume_m3: float) -> dict:
    """State densities/current diagnostics; no field or deposition solve implied."""
    if not np.isfinite(volume_m3) or volume_m3 <= 0:
        raise ValueError("volume_m3 must be finite and positive")
    density: dict[str, float] = {}
    current = np.zeros(3)
    charge = 0.0
    for p in state.particles:
        species = catalog.species_for_state(p.state_id)
        density[p.state_id] = density.get(p.state_id, 0.0) + p.weight / volume_m3
        charge += species.charge_c * p.weight / volume_m3
        current += species.charge_c * p.weight * np.asarray(p.momentum_kg_m_s) / (species.mass_kg * volume_m3)
    return {"state_density_m3": density, "charge_density_c_m3": charge,
            "current_density_a_m2": current.tolist()}


def _restore_random_state_on_failure(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        rng = kwargs.get("rng")
        checkpoint = deepcopy(rng.bit_generator.state) if isinstance(rng, np.random.Generator) else None
        try:
            return function(*args, **kwargs)
        except Exception:
            if checkpoint is not None:
                rng.bit_generator.state = checkpoint
            raise
    return wrapped


@_restore_random_state_on_failure
def advance_gap(
    state: GapState,
    catalog: StateCatalog,
    dt_s: float,
    *,
    collision_volume_m3: float,
    electric_field_v_m,
    magnetic_field_t,
    channels: tuple[ReactionChannel, ...] = (),
    differential_kernels: Mapping[str, Callable] | None = None,
    rng: np.random.Generator,
    max_events: int = 10_000,
    max_particles: int = 10_000,
    max_gyro_angle_rad: float = 0.2,
    max_orbit_substeps: int = 10_000,
    max_nonrel_beta: float = 0.05,
) -> GapStep:
    """Orbit and collide equal-weight particles with explicitly supplied inputs.

    A kernel receives `(channel, catalog, incoming_momenta, rng)` and returns
    `(jacobi_energy_fractions, jacobi_directions)`. Product order is channel order.
    Rates for each unordered pair are `weight * sigma(g) * g / volume`.
    Querying missing cross-section support raises, rather than assuming zero.
    Unsupported steps raise without changing input particles or the supplied
    RNG state. Kernel callbacks must not modify external state. A broader host
    proposal must still checkpoint its fields, caches and integration history.

    This coupled operator is nonrelativistic and enforces a speed-validity cap;
    the separate orbit module supports relativistic motion. Large gyro steps
    are subcycled with a hard work limit, not silently accepted or clipped.
    Uniform prescribed E/B are the only field representation here.
    """
    for name, value in (("dt_s", dt_s), ("collision_volume_m3", collision_volume_m3),
                        ("max_gyro_angle_rad", max_gyro_angle_rad)):
        if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    final_time = state.time_s + dt_s
    if not np.isfinite(final_time) or final_time <= state.time_s:
        raise ValueError("dt_s must produce finite, strictly increasing simulation time")
    for name, value in (("max_events", max_events), ("max_particles", max_particles),
                        ("max_orbit_substeps", max_orbit_substeps)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value <= 0:
            raise ValueError(f"{name} must be a positive integer")
    if not isinstance(rng, np.random.Generator):
        raise ValueError("an explicit numpy random Generator is required")
    if not np.isfinite(max_nonrel_beta) or not 0 < max_nonrel_beta <= 0.05:
        raise ValueError("max_nonrel_beta must lie in (0, 0.05]")
    electric = np.asarray(_vector(electric_field_v_m, "electric_field_v_m"))
    magnetic = np.asarray(_vector(magnetic_field_t, "magnetic_field_t"))
    kernels = differential_kernels or {}
    channels = tuple(channels)
    names = [c.name for c in channels]
    if len(names) != len(set(names)):
        raise ValueError("channel names must be unique")
    for channel in channels:
        channel.validate(catalog, production=True)
        if not channel.enabled or len(channel.reactants) != 2 or len(channel.products) < 2:
            raise ValueError("enabled two-reactant channels with at least two products are required")
        if channel.kinematics_model != "jacobi_energy_partition":
            raise ValueError("coupled operator requires jacobi_energy_partition channels")
        if channel.name not in kernels or not callable(kernels[channel.name]):
            raise ValueError(f"missing differential kernel for {channel.name}")
    particles = list(state.particles)
    if len(particles) > max_particles:
        raise ValueError("particle count exceeds max_particles")
    if particles and any(p.weight != particles[0].weight for p in particles):
        raise ValueError("equal macro-particle weights are required; splitting is not implemented")
    weight = particles[0].weight if particles else 1.0
    total_energy_j(state, catalog)  # Reject a nonrepresentable weighted ledger.
    masses = np.array([catalog.mass_kg(p.state_id) for p in particles])
    charges = np.array([catalog.species_for_state(p.state_id).charge_c for p in particles])
    x = np.array([p.position_m for p in particles], dtype=float).reshape(-1, 3)
    momentum = np.array([p.momentum_kg_m_s for p in particles], dtype=float).reshape(-1, 3)

    def check_speed(p, mass):
        if len(mass) and np.max(np.linalg.norm(p / mass[:, None], axis=1)) > max_nonrel_beta * 299_792_458.0:
            raise ValueError("nonrelativistic reaction/operator speed validity exceeded")

    check_speed(momentum, masses)
    angle = float(np.max(np.abs(charges) / masses, initial=0.0)) * np.linalg.norm(magnetic) * dt_s
    if not np.isfinite(angle) or angle / max_gyro_angle_rad > max_orbit_substeps:
        raise ValueError("required orbit subcycling exceeds max_orbit_substeps")
    substeps = max(1, math.ceil(angle / max_gyro_angle_rad))
    field_gain = 0.0
    displacement_work = 0.0
    for _ in range(substeps):
        step = advance_orbits(x, momentum, masses, charges, dt_s / substeps,
                              np.broadcast_to(electric, x.shape), np.broadcast_to(magnetic, x.shape),
                              relativistic=False)
        x, momentum = step.position_m, step.momentum_kg_m_s
        check_speed(momentum, masses)
        field_gain = math.fsum((field_gain, weight * math.fsum(step.kinetic_delta_j)))
        displacement_work = math.fsum((displacement_work, weight * math.fsum(step.electric_displacement_work_j)))
    particles = [Particle(p.state_id, tuple(x[i]), tuple(momentum[i]), weight)
                 for i, p in enumerate(particles)]

    elapsed = 0.0
    counts = dict.fromkeys(names, 0)
    reaction_residual = 0.0
    momentum_residual = np.zeros(3)
    while channels:
        candidates = []
        for i in range(len(particles)):
            for j in range(i + 1, len(particles)):
                for channel in channels:
                    ids = (particles[i].state_id, particles[j].state_id)
                    if ids == channel.reactants:
                        indices = (i, j)
                    elif ids[::-1] == channel.reactants:
                        indices = (j, i)
                    else:
                        continue
                    pair = [particles[k] for k in indices]
                    m1, m2 = (catalog.mass_kg(p.state_id) for p in pair)
                    relative = np.asarray(pair[0].momentum_kg_m_s) / m1 - np.asarray(pair[1].momentum_kg_m_s) / m2
                    speed = float(np.linalg.norm(relative))
                    if speed == 0:
                        continue
                    reduced_mass = m1 * (m2 / (m1 + m2))
                    energy_ev = 0.5 * reduced_mass * speed**2 / ELEMENTARY_CHARGE_C
                    rate = weight * float(channel.sigma_at(energy_ev, catalog)) * speed / collision_volume_m3
                    if not np.isfinite(rate):
                        raise ValueError("collision rate exceeds finite numeric range")
                    if rate > 0:
                        candidates.append((rate, indices, channel))
        total_rate = math.fsum(item[0] for item in candidates)
        if total_rate == 0:
            break
        if not np.isfinite(total_rate):
            raise ValueError("total collision rate exceeds finite numeric range")
        wait = float(rng.exponential(1.0 / total_rate))
        if elapsed + wait >= dt_s:
            break
        if sum(counts.values()) >= max_events:
            raise ValueError("collision work exceeds max_events; shorten the step")
        elapsed += wait
        pick = float(rng.random()) * total_rate
        cumulative = 0.0
        selected = candidates[-1]
        for candidate in candidates:
            cumulative += candidate[0]
            if pick < cumulative:
                selected = candidate
                break
        _, indices, channel = selected
        pair = [particles[k] for k in indices]
        incoming = np.array([p.momentum_kg_m_s for p in pair])
        fractions, directions = kernels[channel.name](channel, catalog, incoming.copy(), rng)
        outcome = react(channel, catalog, incoming, jacobi_energy_fractions=fractions,
                        jacobi_directions=directions, max_nonrel_beta=max_nonrel_beta)
        if len(particles) - 2 + len(channel.products) > max_particles:
            raise ValueError("reaction products would exceed max_particles")
        pair_mass = np.array([catalog.mass_kg(p.state_id) for p in pair])
        location = np.average([p.position_m for p in pair], axis=0, weights=pair_mass)
        out_momentum = outcome.outgoing_momenta_kg_m_s
        reaction_residual = math.fsum((reaction_residual, weight * outcome.energy_residual_j))
        momentum_residual += weight * (out_momentum.sum(axis=0) - incoming.sum(axis=0))
        particles = [p for k, p in enumerate(particles) if k not in indices] + [
            Particle(s, tuple(location), tuple(p), weight) for s, p in zip(channel.products, out_momentum)]
        counts[channel.name] += 1
    final = GapState(tuple(particles), final_time)
    total_energy_j(final, catalog)
    before_counts = Counter(p.state_id for p in state.particles)
    after_counts = Counter(p.state_id for p in final.particles)
    # First form exact integer population changes, then sum common-reference
    # state energies in eV. Do not subtract large rounded absolute energies.
    internal_change = weight * ELEMENTARY_CHARGE_C * math.fsum(
        (after_counts[s] - before_counts[s]) * catalog.energy_ev(s)
        for s in before_counts.keys() | after_counts.keys())
    kinetic_changes = [p.weight * _particle_kinetic_j(p, catalog) for p in final.particles]
    kinetic_changes.extend(-p.weight * _particle_kinetic_j(p, catalog) for p in state.particles)
    balance_residual = math.fsum((*kinetic_changes, internal_change, -field_gain))
    if not np.all(np.isfinite([field_gain, displacement_work, reaction_residual, balance_residual, *momentum_residual])):
        raise ValueError("weighted energy or momentum ledger exceeds finite numeric range")
    return GapStep(final, MappingProxyType(counts), substeps, field_gain,
                   displacement_work, field_gain - displacement_work,
                   reaction_residual, balance_residual,
                   tuple(float(v) for v in momentum_residual))
