"""Exact collisionless trajectories and conservative width/normal observations.

Inputs are already debited, weighted source packets with their real birth times.
No chemical yields, launch law, electric/magnetic force or heating is invented.
The current direction is invariant/periodic; all three velocity components are
retained for energy, while density is a projection over its declared length.
"""
from __future__ import annotations

import math
import numpy as np

SPECIES = ("H", "H2", "H2O", "OH")


def _array(value, name, ndim):
    raw = np.asarray(value)
    if raw.ndim != ndim or raw.dtype.kind not in "iuf" or not raw.size:
        raise ValueError(f"{name} must be a nonempty real rank-{ndim} array")
    result = np.asarray(raw, dtype=float)
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must be finite")
    return result


def packet_loss_times(position_m, velocity_m_s, domain_bounds_m, foil_bounds_m):
    """Return time since birth of first box exit/solid entry; tangencies pass.

    Box edges absorb at first crossing. The solid is a convex rectangle with
    open interior; a point/ray merely touching its boundary is not absorbed.
    This exact segment test does not depend on animation frame spacing.
    """
    positions = _array(position_m, "position_m", 2)
    velocities = _array(velocity_m_s, "velocity_m_s", 2)
    if positions.shape != velocities.shape or positions.shape[1] != 3:
        raise ValueError("position and velocity arrays must share shape (particle,3)")
    domain = _array(domain_bounds_m, "domain_bounds_m", 1)
    foil = _array(foil_bounds_m, "foil_bounds_m", 1)
    if domain.shape != (4,) or foil.shape != (4,):
        raise ValueError("bounds must be [xmin,xmax,ymin,ymax]")
    if any(bounds[0] >= bounds[1] or bounds[2] >= bounds[3] for bounds in (domain, foil)):
        raise ValueError("rectangles must have positive extents")
    if not (domain[0] < foil[0] < foil[1] < domain[1] and domain[2] < foil[2] < foil[3] < domain[3]):
        raise ValueError("the observation box must strictly contain the foil")
    inside = ((positions[:, 0] > foil[0]) & (positions[:, 0] < foil[1])
              & (positions[:, 1] > foil[2]) & (positions[:, 1] < foil[3]))
    if inside.any():
        raise ValueError("source particles cannot start in the solid interior")
    count = len(positions)
    exits = np.full(count, np.inf)
    enter, leave = np.full(count, -np.inf), np.full(count, np.inf)
    possible = np.ones(count, dtype=bool)
    for axis in (0, 1):
        p, v = positions[:, axis], velocities[:, axis]
        lower, upper = domain[2*axis:2*axis+2]
        if np.any((p < lower) | (p > upper)):
            raise ValueError("source particles must start within the observation box")
        positive, negative = v > 0, v < 0
        times = np.full(count, np.inf)
        times[positive] = (upper - p[positive]) / v[positive]
        times[negative] = (lower - p[negative]) / v[negative]
        exits = np.minimum(exits, times)
        lower, upper = foil[2*axis:2*axis+2]
        moving = v != 0
        possible &= moving | ((p > lower) & (p < upper))
        left, right = np.full(count, -np.inf), np.full(count, np.inf)
        a = (lower - p[moving]) / v[moving]
        b = (upper - p[moving]) / v[moving]
        left[moving], right[moving] = np.minimum(a, b), np.maximum(a, b)
        enter, leave = np.maximum(enter, left), np.minimum(leave, right)
    hits = np.where(possible & (leave > np.maximum(enter, 0)), np.maximum(enter, 0), np.inf)
    return exits, hits


def observe_vacuum_packets(*, birth_time_s, position_m, velocity_m_s, particle_weight,
                           species_index, mass_kg, time_s, x_edges_m, y_edges_m,
                           foil_bounds_m, invariant_current_length_m, initial_surface_particles):
    """Observe independent weighted packets at exact requested frame times."""
    birth = _array(birth_time_s, "birth_time_s", 1)
    positions = _array(position_m, "position_m", 2)
    velocity = _array(velocity_m_s, "velocity_m_s", 2)
    weights = _array(particle_weight, "particle_weight", 1)
    masses = _array(mass_kg, "mass_kg", 1)
    species = np.asarray(species_index)
    n = len(birth)
    if (positions.shape != (n, 3) or velocity.shape != positions.shape or weights.shape != (n,)
            or masses.shape != (n,) or species.shape != (n,) or species.dtype.kind not in "iu"
            or np.any((species < 0) | (species >= len(SPECIES)))
            or np.any(weights <= 0) or np.any(masses <= 0)):
        raise ValueError("packet arrays require consistent sizes, positive masses/weights and canonical species indices")
    times = _array(time_s, "time_s", 1)
    x = _array(x_edges_m, "x_edges_m", 1)
    y = _array(y_edges_m, "y_edges_m", 1)
    foil = _array(foil_bounds_m, "foil_bounds_m", 1)
    for values in (times, x, y):
        if len(values) < 2 or np.any(np.diff(values) <= 0):
            raise ValueError("frame times and grid edges must strictly increase")
    if not np.isfinite(invariant_current_length_m) or invariant_current_length_m <= 0:
        raise ValueError("invariant current length must be positive and finite")
    initial = _array(initial_surface_particles, "initial_surface_particles", 1)
    if initial.shape != (4,) or np.any(initial < 0):
        raise ValueError("four nonnegative initial reservoir counts are required")
    # Cell volumes and the solid mask are exact: no cut-cell volume guesses.
    for boundary, edges in ((foil[0], x), (foil[1], x), (foil[2], y), (foil[3], y)):
        if not np.any(edges == boundary):
            raise ValueError("grid edges must align exactly with each foil boundary")
    exits, hits = packet_loss_times(positions, velocity, [x[0], x[-1], y[0], y[-1]], foil)
    loss = np.minimum(exits, hits)
    loss_kind = np.where(exits <= hits, 0, 1)  # box wins an exact tie
    volumes = np.diff(y)[:, None] * np.diff(x)[None, :] * invariant_current_length_m
    solid = ((.5*(y[:-1]+y[1:]))[:, None] > foil[2]) & ((.5*(y[:-1]+y[1:]))[:, None] < foil[3])
    solid = solid & ((.5*(x[:-1]+x[1:]))[None, :] > foil[0]) & ((.5*(x[:-1]+x[1:]))[None, :] < foil[1])
    density = np.zeros((len(times), 4, len(y)-1, len(x)-1))
    ledgers = []
    energy = .5 * masses * np.sum(velocity * velocity, axis=1) * weights
    if not np.isfinite(energy).all():
        raise ValueError("sampled launch energy exceeds the finite domain")
    for frame, stamp in enumerate(times):
        age = stamp - birth
        born, lost = age >= 0, (age >= 0) & (age >= loss)
        live = born & ~lost
        points = positions[live, :2] + age[live, None] * velocity[live, :2]
        if not np.isfinite(points).all():
            raise ValueError("ballistic positions exceed the finite domain")
        # Boundary points occupy the gas side for histogram ownership only.
        # Exact birth positions/trajectories are not advanced or displaced.
        bin_points = points.copy()
        for axis in (0, 1):
            other = 1-axis
            face_span = ((points[:, other] >= foil[2*other])
                         & (points[:, other] <= foil[2*other+1]))
            for endpoint, direction in ((foil[2*axis], -np.inf), (foil[2*axis+1], np.inf)):
                at_face = face_span & (points[:, axis] == endpoint)
                bin_points[at_face, axis] = np.nextafter(endpoint, direction)
        live_species, live_weights = species[live], weights[live]
        rows = []
        for s, name in enumerate(SPECIES):
            select = live_species == s
            numbers, _, _ = np.histogram2d(bin_points[select, 1], bin_points[select, 0], bins=(y, x), weights=live_weights[select])
            if np.any(numbers[solid] != 0):
                raise ValueError("a live particle entered a solid cell")
            density[frame, s] = numbers / volumes
            this = species == s
            launched = math.fsum(map(float, weights[born & this]))
            inside = math.fsum(map(float, weights[live & this]))
            escaped = math.fsum(map(float, weights[lost & this & (loss_kind == 0)]))
            absorbed = math.fsum(map(float, weights[lost & this & (loss_kind == 1)]))
            remaining = math.fsum((float(initial[s]), -launched))
            if remaining < -1e-12 * max(float(initial[s]), launched):
                raise ValueError("source packets exceed the initial surface inventory")
            histogram = float(np.sum(density[frame, s] * volumes))
            if not math.isclose(histogram, inside, rel_tol=2e-12, abs_tol=0):
                raise ValueError("histogram volume integral lost live particles")
            rows.append({"species": name, "initial_surface_particles": float(initial[s]),
                         "remaining_surface_particles": remaining, "launched_particles": launched,
                         "inside_particles": inside, "escaped_particles": escaped, "absorbed_particles": absorbed,
                         "particle_balance_residual": math.fsum((launched, -inside, -escaped, -absorbed)),
                         "histogram_balance_residual": histogram-inside,
                         "launched_translational_energy_j": math.fsum(map(float, energy[born & this]))})
        ledgers.append({"time_s": float(stamp), "species": rows})
    return {"time_s": times, "x_edges_m": x, "y_edges_m": y, "species": np.asarray(SPECIES),
            "number_density_m3": density, "foil_bounds_m": foil}, ledgers
