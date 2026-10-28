"""Explicit numerical representation of exhausted source tails in packets.

The physical release history is never changed. A sampling-only interval may
omit components below a declared particle-count threshold, retaining an exact
decimal account and a strict per-component upper bound. Qualified depletion
and launch laws remain unchanged for every represented independent component.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
import math
from numbers import Integral, Real

import numpy as np

from .planar_release import PlanarSurfaceHistory
from .temkin_planar import TemkinPlanarReleaseHistory, evolve_temkin_planar_release


def prepare_sampling_interval(source, index, minimum_component_particles=1e-250):
    """Return a validated sampling history/index and explicit omitted-tail report.

    The cutoff is a numerical packet-representation choice, not a desorption
    rate, inventory floor or modification of the saved physical source. Exact
    Decimal products avoid underflow even on Windows where longdouble=float64.
    """
    if not isinstance(source, TemkinPlanarReleaseHistory):
        raise ValueError("source must be a TemkinPlanarReleaseHistory")
    if (isinstance(index, (bool, np.bool_)) or not isinstance(index, Integral)
            or not 0 <= index < len(source.time_edges_s)-1):
        raise ValueError("index must identify one source interval")
    if (isinstance(minimum_component_particles, (bool, np.bool_))
            or not isinstance(minimum_component_particles, Real)
            or not math.isfinite(minimum_component_particles) or minimum_component_particles <= 0):
        raise ValueError("minimum_component_particles must be finite and positive")
    index = int(index)
    cutoff = float(minimum_component_particles)
    released = source.released_m2[index]
    areas = source.patch_areas_m2
    with np.errstate(under="ignore"):
        products = released.astype(np.longdouble) * areas.astype(np.longdouble)[None]
    candidates = (released > 0) & (products < np.longdouble(cutoff))
    mask = np.zeros_like(candidates)
    omitted = [Decimal(0) for _ in source.reservoirs]
    bounds = [Decimal(0) for _ in source.reservoirs]
    with localcontext() as context:
        # Exact products/sums of finite binary64 inputs require at most ~2200
        # decimal digits, including subnormal operands and cancellation-free sums.
        context.prec = 4096
        threshold = Decimal.from_float(cutoff)
        for pool, current, width in zip(*np.nonzero(candidates)):
            count = (Decimal.from_float(float(released[pool, current, width]))
                     * Decimal.from_float(float(areas[current, width])))
            if count < threshold:
                mask[pool, current, width] = True
                omitted[pool] += count
                bounds[pool] += threshold
        total = sum(omitted, Decimal(0))
        total_bound = sum(bounds, Decimal(0))
    report = {
        "minimum_component_particles": cutoff,
        "omitted_positive_component_count": int(np.count_nonzero(mask)),
        "omitted_component_count": int(np.count_nonzero(mask)),
        "time_interval_s": source.time_edges_s[index:index+2].tolist(),
        "omitted_particles": float(total),
        "omitted_particles_exact_decimal": str(total),
        "omitted_particles_upper_bound": float(total_bound),
        "omitted_particles_upper_bound_exact_decimal": str(total_bound),
        "per_pool": [{"pool": pool.name, "state_id": pool.state_id,
            "omitted_particles": float(count), "omitted_particles_exact_decimal": str(count),
            "omitted_particles_upper_bound": float(bound),
            "omitted_particles_upper_bound_exact_decimal": str(bound)}
            for pool, count, bound in zip(source.reservoirs, omitted, bounds)],
        "physical_source_history_changed": False,
        "sampling_representation": "original_history" if not mask.any() else "single_interval_with_declared_unrepresented_tails",
    }
    if not mask.any():
        return source, index, report
    if not np.any((released > 0) & ~mask):
        return None, 0, report
    coverage = source.coverage_m2[index].copy()
    coverage[mask] = 0
    surface = PlanarSurfaceHistory(source.time_edges_s[index:index+2], source.x_edges_m, source.z_edges_m,
        np.stack((source.temperature_k[index], source.temperature_k[index])), source.surface_attributes)
    sampling = evolve_temkin_planar_release(surface, coverage, source.reservoirs,
        parameter_provenance=source.parameter_provenance)
    # Independent pools/patches must be bit-identical outside the omitted tail.
    for name in ("released_m2", "translational_energy_j_m2", "internal_energy_j_m2",
                 "binding_energy_j_m2", "total_wall_energy_debit_j_m2"):
        if not np.array_equal(getattr(sampling, name)[0][~mask], getattr(source, name)[index][~mask]):
            raise ValueError("sampling view changed a represented source component")
    return sampling, 0, report
