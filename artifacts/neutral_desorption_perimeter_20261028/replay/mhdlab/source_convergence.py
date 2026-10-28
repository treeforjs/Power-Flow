"""Manufactured planar source observables and finite-resolution comparisons.

The temperature fixture has exact spatial bin averages and explicit snapshot
knots. Source statistics integrate the locally depleted piecewise-constant
patches, including within-patch position variance. They never infer a spectrum,
advance gap particles or interpret a finite finest mesh as an exact solution.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import fields
import math
from numbers import Real

import numpy as np

from .constants import KB
from .flash_surface import FlashSurfaceHistory
from .planar_release import PlanarSurfaceHistory
from .temkin_planar import TemkinPlanarReleaseHistory


SCHEMA = "mhdlab-planar-source-observables-v1"
FIXTURE_IDENTITY = "quadratic_time_quadratic_width_linear_current_temperature_v1"
_EPS = float(np.finfo(float).eps)


def _positive(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite positive real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} must be a finite positive real number") from error
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be a finite positive real number")
    return result


def _count(value, name):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
            or value < 1 or value > 100_000):
        raise ValueError(f"{name} must be an integer in [1,100000]")
    return int(value)


def manufactured_planar_surface(*, width_bins, current_bins, snapshot_intervals,
                               duration_s, width_m, current_length_m):
    """Return exact spatial bin means of an explicitly manufactured field.

    Coordinates span x=[-W/2,W/2], z=[0,L], and t=[0,duration]. The point field is
    300+(t/duration)^2*[250+1200*(2*x/W)^2*(0.6+0.8*z/L)] kelvin. Snapshot
    temperatures are values at their knots, not temporal averages. A downstream
    source uses piecewise-linear knot interpolation; refining its subdivisions
    alone does not remove the snapshot interpolation error.
    """
    nx, nz, ns = (_count(value, name) for value, name in (
        (width_bins, "width_bins"), (current_bins, "current_bins"),
        (snapshot_intervals, "snapshot_intervals")))
    duration = _positive(duration_s, "duration_s")
    width = _positive(width_m, "width_m")
    length = _positive(current_length_m, "current_length_m")
    if (ns + 1) * nx * nz > 20_000_000:
        raise ValueError("manufactured surface exceeds the reference memory limit")
    x = np.linspace(-.5 * width, .5 * width, nx + 1)
    z = np.linspace(0., length, nz + 1)
    time = np.linspace(0., duration, ns + 1)
    # Work in dimensionless coordinates to avoid squaring dimensional widths.
    ux = x / width
    uz = z / length
    x_square_mean = 4 * (ux[:-1] ** 2 + ux[:-1] * ux[1:] + ux[1:] ** 2) / 3
    current_factor = .6 + .8 * (.5 * uz[:-1] + .5 * uz[1:])
    amplitude = 250 + 1200 * current_factor[:, None] * x_square_mean[None, :]
    temperature = 300 + (time / duration)[:, None, None] ** 2 * amplitude
    source = FlashSurfaceHistory([0., 1.], [0., 1.], [[300.], [300.]], "synthetic",
                                "manufactured-source-convergence", "analytic polynomial patch-average fixture",
                                "stationary manufactured single face")
    attributes = source.attributes
    attributes.update(manufactured_fixture_identity=FIXTURE_IDENTITY, thermal_input_status="manufactured",
                      manufactured_point_temperature_formula="300+(t/duration_s)^2*(250+1200*(2*x/width_m)^2*(0.6+0.8*z/current_length_m)) K",
                      manufactured_snapshot_interpolation="exact knot values; downstream linear interpolation requires independent snapshot refinement")
    return PlanarSurfaceHistory(time, x, z, temperature, attributes)


def _sum(value):
    result = math.fsum(float(item) for item in np.asarray(value).flat)
    if not math.isfinite(result):
        raise ValueError("source observable exceeds the finite numeric range")
    return result


def _finite(value, name):
    if not math.isfinite(value):
        raise ValueError(f"{name} exceeds the finite numeric range")
    return float(value)


def _thermal_speed(temperature, mass):
    """Evaluate sqrt(kB*T/m) without losing kB*T before mass division."""
    temperature_mantissa, temperature_exponent = np.frexp(temperature)
    kb_mantissa, kb_exponent = math.frexp(KB)
    mass_mantissa, mass_exponent = math.frexp(mass)
    exponent = temperature_exponent + kb_exponent - mass_exponent
    mantissa = temperature_mantissa * kb_mantissa / mass_mantissa
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        result = np.ldexp(np.sqrt(mantissa * np.where(exponent % 2 == 0, 1., 2.)), exponent // 2)
    if not np.all(np.isfinite(result)) or np.any(result <= 0):
        raise ValueError("positive thermal speed is outside the representable source-observable range")
    return result


def source_observables(history: TemkinPlanarReleaseHistory, *, state_masses_kg,
                       include_time_moments=False, time_quadrature_order=16):
    """Integrate local source moments over time and actual stationary face area.

    Each returned metric is a per-reservoir list. Normalized moments are null
    for a reservoir with zero emitted number. Width redistribution preserves
    each current patch's time/state/thermal mixture, so its temporal and thermal
    moments equal the resolved source's. Velocity means refer to outward-normal
    magnitude in an outgoing flux Maxwellian; no global signed y is assigned.
    Chemistry parameter status is not inferred from free-text provenance.
    """
    if not isinstance(history, TemkinPlanarReleaseHistory):
        raise ValueError("history must be a TemkinPlanarReleaseHistory")
    if not isinstance(include_time_moments, (bool, np.bool_)):
        raise ValueError("include_time_moments must be a boolean")
    if (isinstance(time_quadrature_order, (bool, np.bool_))
            or not isinstance(time_quadrature_order, (int, np.integer)) or time_quadrature_order < 2):
        raise ValueError("time_quadrature_order must be an integer at least two")
    if not isinstance(state_masses_kg, Mapping):
        raise ValueError("state_masses_kg must map exact state IDs to positive masses")
    masses = []
    for pool in history.reservoirs:
        if pool.state_id not in state_masses_kg:
            raise ValueError(f"missing exact-state mass for {pool.state_id!r}; no species fallback is allowed")
        masses.append(_positive(state_masses_kg[pool.state_id], f"mass for {pool.state_id}"))
    history = TemkinPlanarReleaseHistory(**{field.name: getattr(history, field.name)
                                          for field in fields(TemkinPlanarReleaseHistory)})
    energy_fields = {
        "translation_energy_j": "translational_energy_j_m2", "internal_energy_j": "internal_energy_j_m2",
        "binding_energy_j": "binding_energy_j_m2", "total_wall_energy_debit_j": "total_wall_energy_debit_j_m2",
    }
    units = {name: "1" for name in ("initial_number", "released_number", "remaining_number", "inventory_balance_residual_number")}
    units.update({name: "J" for name in energy_fields})
    for prefix in ("resolved", "uniform_width"):
        units.update({f"{prefix}_mean_x_m": "m", f"{prefix}_mean_z_m": "m",
                      f"{prefix}_second_raw_x_m2": "m^2", f"{prefix}_second_raw_z_m2": "m^2",
                      f"{prefix}_variance_x_m2": "m^2", f"{prefix}_variance_z_m2": "m^2",
                      f"{prefix}_covariance_xz_m2": "m^2"})
    units.update(mean_velocity_normal_m_s="m s^-1", mean_speed_squared_m2_s2="m^2 s^-2")
    if include_time_moments:
        from .temkin_emission_time import FrozenTemkinEmission
        from .emission_time_moments import conditional_time_moments, QUADRATURE_IDENTITY
        units.update(emission_mean_time_s="s", emission_variance_time_s2="s^2")
    metrics = {name: [] for name in units}
    area = history.patch_areas_m2
    x_origin, z_origin = float(history.x_edges_m[0]), float(history.z_edges_m[0])
    x_center_offset = .5 * (history.x_edges_m[:-1] - x_origin) + .5 * (history.x_edges_m[1:] - x_origin)
    z_center_offset = .5 * (history.z_edges_m[:-1] - z_origin) + .5 * (history.z_edges_m[1:] - z_origin)
    x_variance = np.diff(history.x_edges_m) ** 2 / 12
    z_variance = np.diff(history.z_edges_m) ** 2 / 12
    uniform_x_mean = .5 * history.x_edges_m[0] + .5 * history.x_edges_m[-1]
    uniform_x_variance = (history.x_edges_m[-1] - history.x_edges_m[0]) ** 2 / 12
    empty_reservoirs = []
    time_underflow = []
    for p, pool in enumerate(history.reservoirs):
        number = history.released_m2[:, p] * area
        total = _sum(number)
        metrics["initial_number"].append(_sum(history.coverage_m2[0, p] * area))
        metrics["released_number"].append(total)
        metrics["remaining_number"].append(_sum(history.coverage_m2[-1, p] * area))
        metrics["inventory_balance_residual_number"].append(_sum(history.inventory_balance_residual_m2[:, p] * area))
        for name, field in energy_fields.items():
            metrics[name].append(_sum(getattr(history, field)[:, p] * area))
        if total == 0:
            empty_reservoirs.append(pool.name)
            for name in metrics:
                if len(metrics[name]) <= p:
                    metrics[name].append(None)
            continue
        weights = number / total
        if np.any((number > 0) & (weights < np.finfo(float).tiny)):
            raise ValueError("positive source probability is below the qualified normal range for observables")
        mean_x_offset = _sum(weights * x_center_offset[None, None, :])
        mean_z_offset = _sum(weights * z_center_offset[None, :, None])
        mean_x, mean_z = x_origin + mean_x_offset, z_origin + mean_z_offset
        variance_x = _sum(weights * ((x_center_offset[None, None, :] - mean_x_offset) ** 2 + x_variance[None, None, :]))
        variance_z = _sum(weights * ((z_center_offset[None, :, None] - mean_z_offset) ** 2 + z_variance[None, :, None]))
        covariance = _sum(weights * (x_center_offset[None, None, :] - mean_x_offset) * (z_center_offset[None, :, None] - mean_z_offset))
        resolved = dict(mean_x_m=mean_x, mean_z_m=mean_z,
                        second_raw_x_m2=mean_x ** 2 + variance_x,
                        second_raw_z_m2=mean_z ** 2 + variance_z,
                        variance_x_m2=variance_x, variance_z_m2=variance_z, covariance_xz_m2=covariance)
        uniform = dict(resolved, mean_x_m=uniform_x_mean, variance_x_m2=uniform_x_variance,
                       second_raw_x_m2=uniform_x_mean ** 2 + uniform_x_variance, covariance_xz_m2=0.)
        for prefix, moments in (("resolved", resolved), ("uniform_width", uniform)):
            for key, value in moments.items():
                metrics[f"{prefix}_{key}"].append(_finite(value, key))
        thermal_speed = np.zeros(number.shape)
        active = number > 0
        thermal_speed[active] = _thermal_speed(history.temperature_k[active], masses[p])
        with np.errstate(over="ignore", invalid="ignore", under="ignore"):
            speed_squared = 4 * thermal_speed ** 2
        if not np.all(np.isfinite(speed_squared[active])) or np.any(speed_squared[active] <= 0):
            raise ValueError("positive thermal speed squared is outside the representable source-observable range")
        metrics["mean_velocity_normal_m_s"].append(_sum(weights * (math.sqrt(math.pi / 2) * thermal_speed)))
        metrics["mean_speed_squared_m2_s2"].append(_sum(weights * speed_squared))
        if include_time_moments:
            local_times = np.zeros(number.shape)
            local_variances = np.zeros(number.shape)
            origin = float(history.time_edges_s[0])
            for interval, z, x in np.argwhere(number > 0):
                emission = FrozenTemkinEmission(pool, float(history.coverage_m2[interval, p, z, x]),
                                                float(history.temperature_k[interval, z, x]),
                                                float(history.time_edges_s[interval + 1] - history.time_edges_s[interval]))
                timing = conditional_time_moments(emission, quadrature_order=int(time_quadrature_order))
                local_times[interval, z, x] = float(history.time_edges_s[interval] - origin) + timing.mean_offset_s
                local_variances[interval, z, x] = timing.variance_offset_s2
                if timing.mean_underflowed or timing.variance_underflowed:
                    time_underflow.append({"reservoir": pool.name, "interval": int(interval), "current_patch": int(z), "width_patch": int(x),
                                           "mean_underflowed": timing.mean_underflowed, "variance_underflowed": timing.variance_underflowed})
            mean_time_offset = _sum(weights * local_times)
            metrics["emission_mean_time_s"].append(_finite(origin + mean_time_offset, "emission mean time"))
            metrics["emission_variance_time_s2"].append(_sum(weights * (local_variances + (local_times - mean_time_offset) ** 2)))
    return {
        "schema": SCHEMA, "units": units, "metrics": metrics,
        "reservoir_names": [pool.name for pool in history.reservoirs],
        "state_ids": [pool.state_id for pool in history.reservoirs], "state_masses_kg": masses,
        "source_kind": history.surface_attributes["source_kind"],
        "thermal_input_status": history.surface_attributes.get("thermal_input_status", "not_inferred_from_history"),
        "parameter_provenance": history.parameter_provenance,
        "law_identity": history.law_identity, "source_time_model": history.time_model,
        "source_history_schema": history.schema, "coverage_model": history.coverage_model,
        "qualifications": {
            "spatial_model": "independent uniform rectangles within each source patch; exact within-patch second moments",
            "uniform_width_model": "width positions redistributed independently while retaining each current/time/state/thermal mixture",
            "velocity_model": "classical outgoing flux Maxwellian; normal mean is outward magnitude; no signed global y or transport advance",
            "empty_release_reservoirs": empty_reservoirs,
            "normalized_zero_release_metrics": "null; no emitted-population distribution exists",
            "include_time_moments": bool(include_time_moments),
            "time_quadrature_order": int(time_quadrature_order) if include_time_moments else None,
            "time_quadrature_identity": QUADRATURE_IDENTITY if include_time_moments else None,
            "time_moment_underflow_components": time_underflow,
            "time_reference": "absolute mean time; centered variance assembled relative to first supplied time edge",
            "finite_reference": "no exact continuum source or experimental validation is inferred",
            "minimum_positive_source_probability": float(np.finfo(float).tiny),
        },
    }


def _numeric(value, name, *, allow_none=False):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked values")
    raw = np.asarray(value, dtype=object)
    if raw.size == 0:
        raise ValueError(f"{name} must be nonempty")
    result = np.empty(raw.shape, dtype=float)
    missing = np.zeros(raw.shape, dtype=bool)
    for index in np.ndindex(raw.shape):
        item = raw[index]
        if item is None and allow_none:
            missing[index] = True
            result[index] = 0.
        elif isinstance(item, (bool, np.bool_)) or not isinstance(item, Real):
            raise ValueError(f"{name} must contain finite real numbers")
        else:
            try:
                converted = float(item)
            except (ValueError, OverflowError) as error:
                raise ValueError(f"{name} must contain finite real numbers") from error
            if not math.isfinite(converted):
                raise ValueError(f"{name} must contain finite real numbers")
            result[index] = converted
    return result, missing


def _json_array(value, missing=None):
    result = np.asarray(value, dtype=object)
    if missing is not None:
        result = result.copy()
        result[missing] = None
    return result.tolist()


def _norm(value):
    return _finite(math.hypot(*(float(item) for item in np.asarray(value).flat)), "comparison norm")


def compare_observables(observables, reference):
    """Compare numeric metrics against an explicitly finite finest reference."""
    if not isinstance(observables, Mapping) or not isinstance(reference, Mapping):
        raise ValueError("observables and reference must be metric mappings")
    left, right = observables.get("metrics", observables), reference.get("metrics", reference)
    if not isinstance(left, Mapping) or not isinstance(right, Mapping) or not left or set(left) != set(right):
        raise ValueError("observables must contain the same nonempty metric key set")
    for name in ("units", "reservoir_names", "state_ids", "state_masses_kg"):
        if name in observables and name in reference and observables[name] != reference[name]:
            raise ValueError(f"observable comparison metadata disagree: {name}")
    metrics = {}
    for name in left:
        value, missing_left = _numeric(left[name], name, allow_none=True)
        target, missing_right = _numeric(right[name], name, allow_none=True)
        if value.shape != target.shape:
            raise ValueError(f"metric shapes disagree: {name}")
        missing = missing_left | missing_right
        with np.errstate(over="ignore", invalid="ignore"):
            difference = np.asarray(np.abs(value - target))
        if not np.all(np.isfinite(difference)):
            raise ValueError(f"metric difference exceeds finite range: {name}")
        difference[missing] = 0.
        zero_reference = target == 0
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            relative = np.divide(difference, np.abs(target), out=np.zeros_like(difference), where=~(missing | zero_reference))
        relative_unresolved = ~np.isfinite(relative)
        warnings = []
        if np.any(missing):
            warnings.append("undefined normalized metric; null entries excluded from norms")
        if np.any(zero_reference & ~missing):
            warnings.append("zero reference component; relative difference is undefined")
        if np.any(relative_unresolved):
            warnings.append("relative component exceeds finite range; represented as null")
        delta_norm = _norm(difference) if np.any(~missing) else None
        ref_norm = _norm(target[~missing]) if np.any(~missing) else None
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            relative_norm = float(np.divide(delta_norm, ref_norm)) if ref_norm is not None and ref_norm > 0 else None
        if relative_norm is not None and not math.isfinite(relative_norm):
            relative_norm = None
            warnings.append("relative norm exceeds finite range; represented as null")
        metrics[name] = {
            "absolute_difference": _json_array(difference, missing),
            "relative_difference": _json_array(relative, missing | zero_reference | relative_unresolved),
            "difference_l2": delta_norm, "reference_l2": ref_norm,
            "relative_l2": relative_norm,
            "warnings": warnings,
        }
    return {"schema": "mhdlab-source-observable-comparison-v1",
            "reference_kind": "finite_finest_resolution_reference", "metrics": metrics}


def observed_self_order(coarse, medium, fine, *, refinement_ratio=2.):
    """Estimate order from three actual solution differences, without a proxy truth.

    Roundoff-scale differences, sign reversal/cancellation and nondecreasing
    differences produce explicit warnings. No finite-reference error ratio is
    substituted for the three-level estimate.
    """
    ratio = _positive(refinement_ratio, "refinement_ratio")
    if ratio <= 1:
        raise ValueError("refinement_ratio must exceed one")
    values = [_numeric(value, name)[0] for value, name in ((coarse, "coarse"), (medium, "medium"), (fine, "fine"))]
    if any(value.shape != values[0].shape for value in values):
        raise ValueError("three refinement levels must have identical metric shapes")
    with np.errstate(over="ignore", invalid="ignore"):
        first, second = values[0] - values[1], values[1] - values[2]
    if not np.all(np.isfinite(first)) or not np.all(np.isfinite(second)):
        raise ValueError("solution differences exceed the finite numeric range")
    first_norm, second_norm = _norm(first), _norm(second)
    scale = max(_norm(value) for value in values)
    warnings, order = [], None
    if first_norm <= 64 * _EPS * scale or second_norm <= 64 * _EPS * scale:
        warnings.append("zero or roundoff-scale solution difference; observed order is unresolved")
    else:
        alignment = _sum((first / first_norm) * (second / second_norm))
        if alignment <= 0:
            warnings.append("difference reversal or cancellation; refinement is not demonstrably asymptotic")
        else:
            if first_norm <= second_norm:
                warnings.append("nondecreasing solution difference; refinement is not demonstrably asymptotic")
            else:
                order = (math.log(first_norm) - math.log(second_norm)) / math.log(ratio)
    return {"schema": "mhdlab-source-self-order-v1", "refinement_ratio": ratio,
            "coarse_difference_l2": first_norm, "fine_difference_l2": second_norm,
            "observed_order": order, "warnings": warnings}
