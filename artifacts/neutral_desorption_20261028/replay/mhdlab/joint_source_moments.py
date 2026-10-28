"""Deterministic joint statistics of a stationary planar Temkin boundary source.

An emitted component has its original interval, reservoir, current/width patch,
frozen temperature and law-specific conditional emission time. Position is
uniform within that rectangle and independent of time inside the component.
Spatial second/fourth moments are integrated analytically. The comparison
redistributes width position only and retains the original time/z/state/thermal
mixture. No particle sampling, gap evolution or radiative calculation occurs.

Energy here is the component's expected outgoing translation, 2*kB*T per
particle. Its reported variance/covariances describe variation of that local
conditional mean. They exclude kinetic-energy scatter within a Maxwellian,
binding/internal costs and sampled packet energy.
"""

from __future__ import annotations

from dataclasses import asdict, fields
import math
from numbers import Integral

import numpy as np

from .constants import KB
from .emission_time_moments import (
    MAX_QUADRATURE_ORDER, QUADRATURE_IDENTITY, conditional_time_moments,
)
from .planar_release_io import plain_metadata
from .temkin_emission_time import FrozenTemkinEmission, TIME_IDENTITY
from .temkin_planar import TemkinPlanarReleaseHistory
from .thermal_release import _product


SCHEMA = "mhdlab-joint-source-moments-v1"
MIN_POSITIVE_SOURCE_PROBABILITY = float(np.finfo(float).tiny)
UNITS = {
    "released_number": "1", "mean_time_s": "s", "variance_time_s2": "s^2",
    "mean_x_m": "m", "variance_x_m2": "m^2", "mean_z_m": "m", "variance_z_m2": "m^2",
    "mean_width_centered_square_m2": "m^2", "variance_width_centered_square_m4": "m^4",
    "mean_expected_translational_energy_j": "J", "variance_expected_translational_energy_j2": "J^2",
    "cov_time_x_s_m": "s m", "cov_time_width_centered_square_s_m2": "s m^2",
    "cov_expected_translational_energy_width_centered_square_j_m2": "J m^2",
    "cov_time_z_s_m": "s m", "cov_time_expected_translational_energy_s_j": "s J",
}


class JointSourceMomentRangeError(ValueError):
    """A source statistic is outside the explicitly supported float64 range."""


def _sum(values):
    try:
        result = math.fsum(float(value) for value in values)
    except (ValueError, OverflowError) as error:
        raise JointSourceMomentRangeError("joint source statistic exceeds the finite numeric range") from error
    if not math.isfinite(result):
        raise JointSourceMomentRangeError("joint source statistic exceeds the finite numeric range")
    return result


def _signed_product(*factors):
    sign = -1. if sum(value < 0 for value in factors) % 2 else 1.
    magnitude = _product(*(abs(float(value)) for value in factors))
    if magnitude == 0 and all(value != 0 for value in factors):
        raise JointSourceMomentRangeError("nonzero joint-statistic product is below the representable range")
    return sign * magnitude


def _weighted(weights, *values):
    return _sum(_signed_product(weight, *items) for weight, *items in zip(weights, *values))


def _centered(weights, values):
    """Anchor at the dominant component, retaining exact constant-value limits."""
    anchor = float(values[int(np.argmax(weights))])
    offsets = values - anchor
    mean_offset = _weighted(weights, offsets)
    return _sum((anchor, mean_offset)), offsets - mean_offset


def _two_part_centered(weights, bases, offsets, dominant):
    """Keep base/local information separate until population centering.

    A rare early component must not set the numerical origin for a dominant
    late cluster. Covariances retain both centered parts separately as well.
    Reported absolute means still have their ordinary float64 representation.
    """
    origin = float(bases[dominant])
    relative = np.asarray([_sum((base, -origin)) for base in bases])
    mean_base, centered_base = _centered(weights, relative)
    mean_local, centered_local = _centered(weights, offsets)
    centered = np.asarray([_sum((base, local)) for base, local in zip(centered_base, centered_local)])
    mean = _sum((origin, mean_base, mean_local))
    return mean, (centered_base, centered_local), centered


def _covariance(weights, left_parts, right_parts):
    return _sum(_weighted(weights, left, right) for left in left_parts for right in right_parts)


def _count_fraction(number, width, full_width):
    """Compute N*dx/W before rounding small normalized width fractions."""
    if number == 0:
        return 0.
    n_m, n_e = math.frexp(number)
    w_m, w_e = math.frexp(width)
    total_m, total_e = math.frexp(full_width)
    try:
        result = math.ldexp(n_m * w_m / total_m, n_e + w_e - total_e)
    except OverflowError as error:
        raise JointSourceMomentRangeError("uniform-width conditional number exceeds finite range") from error
    if not math.isfinite(result) or result <= 0:
        raise JointSourceMomentRangeError("positive uniform-width conditional number is not representable")
    return result


def _x_geometry(center, width):
    variance = _signed_product(width, width, 1 / 12)
    square_mean = _sum((_signed_product(center, center), variance))
    square_variance = _sum((_signed_product(center, center, width, width, 1 / 3),
                            _signed_product(width, width, width, width, 1 / 180)))
    return center, variance, square_mean, square_variance


def _empty():
    return {name: 0. if name == "released_number" else None for name in UNITS}


def _statistics(batch, selection, *, x_face_center, x_face_left, x_face_right, z_origin, time_origin,
                uniform_x_geometry=None, emitted_number=None):
    """Use centered mixture formulas; optional x law is independent of mixture."""
    if not np.any(selection):
        return _empty()
    number = batch["number"][selection]
    mixture_total = _sum(number)
    if mixture_total == 0:
        return _empty()
    weights = number / mixture_total
    if np.any(weights < MIN_POSITIVE_SOURCE_PROBABILITY):
        raise JointSourceMomentRangeError("positive source probability is below the qualified normal range")
    count = mixture_total if emitted_number is None else emitted_number
    # Keep the component's interval start/local elapsed mean and patch
    # endpoint halves separate until this population is anchored. Adding
    # a tiny local offset to a distant global origin would lose correlations.
    dominant = int(np.argmax(number))
    starts, offsets = batch["time_start"][selection], batch["time_offset"][selection]
    mean_time_absolute, time_parts, centered_time = _two_part_centered(weights, starts, offsets, dominant)
    energy = batch["energy"][selection]
    z_first, z_last = batch["z_first_half"][selection], batch["z_last_half"][selection]
    mean_z_absolute, z_parts, centered_z = _two_part_centered(weights, z_first, z_last, dominant)
    mean_energy, centered_energy = _centered(weights, energy)
    time_variance = _sum((_weighted(weights, batch["time_variance"][selection]),
                          _weighted(weights, centered_time, centered_time)))
    z_variance = _sum((_weighted(weights, batch["z_variance"][selection]),
                       _weighted(weights, centered_z, centered_z)))
    energy_variance = _weighted(weights, centered_energy, centered_energy)
    if uniform_x_geometry is None:
        x_first, x_last = batch["x_first_half"][selection], batch["x_last_half"][selection]
        mean_x_absolute, x_parts, centered_x = _two_part_centered(weights, x_first, x_last, dominant)
        # Keep the face-center displacement separate from the local centered
        # spread. Cov(A,(x-facecenter)^2)=2*m*Cov(A,x)+Cov(A,(x-meanx)^2).
        # This retains a tiny cluster's correlation far from the face center.
        first_anchor, last_anchor = float(x_first[dominant]), float(x_last[dominant])
        relative_first = np.asarray([_sum((first, -first_anchor)) for first in x_first])
        relative_last = np.asarray([_sum((last, -last_anchor)) for last in x_last])
        mean_first = _weighted(weights, relative_first)
        mean_last = _weighted(weights, relative_last)
        mean_x_face_displacement = _sum((first_anchor, last_anchor, -.5 * x_face_left, -.5 * x_face_right, mean_first, mean_last))
        x_local_variance = batch["x_variance"][selection]
        square_local_variance = batch["square_variance"][selection]
        x_variance = _sum((_weighted(weights, x_local_variance), _weighted(weights, centered_x, centered_x)))
        centered_square_mean = np.asarray([_sum((_signed_product(center, center), variance))
                                           for center, variance in zip(centered_x, x_local_variance)])
        _, centered_square_mean_delta = _centered(weights, centered_square_mean)
        shifted_square = np.asarray([_sum((_signed_product(2., mean_x_face_displacement, center), square))
                                     for center, square in zip(centered_x, centered_square_mean)])
        _, shifted_square_delta = _centered(weights, shifted_square)
        mean_square = _sum((_signed_product(mean_x_face_displacement, mean_x_face_displacement), x_variance))
        square_variance = _sum((_weighted(weights, square_local_variance),
                                _weighted(weights, shifted_square_delta, shifted_square_delta)))
        time_x = _covariance(weights, time_parts, x_parts)
        time_square = _sum((_signed_product(2., mean_x_face_displacement, time_x),
                            _covariance(weights, time_parts, (centered_square_mean_delta,))))
        energy_x = _covariance(weights, (centered_energy,), x_parts)
        energy_square = _sum((_signed_product(2., mean_x_face_displacement, energy_x),
                              _covariance(weights, (centered_energy,), (centered_square_mean_delta,))))
    else:
        _, x_variance, mean_square, square_variance, mean_x_absolute = uniform_x_geometry
        # The randomized x position is independent of the retained joint
        # time/z/state/temperature mixture, including within conditional bins.
        time_x, time_square, energy_square = 0., 0., 0.
    values = {
        "released_number": count, "mean_time_s": mean_time_absolute,
        "variance_time_s2": time_variance,
        "mean_x_m": mean_x_absolute, "variance_x_m2": x_variance,
        "mean_z_m": mean_z_absolute, "variance_z_m2": z_variance,
        "mean_width_centered_square_m2": mean_square, "variance_width_centered_square_m4": square_variance,
        "mean_expected_translational_energy_j": mean_energy,
        "variance_expected_translational_energy_j2": energy_variance,
        "cov_time_x_s_m": time_x, "cov_time_width_centered_square_s_m2": time_square,
        "cov_expected_translational_energy_width_centered_square_j_m2": energy_square,
        "cov_time_z_s_m": _covariance(weights, time_parts, z_parts),
        "cov_time_expected_translational_energy_s_j": _covariance(weights, time_parts, (centered_energy,)),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise JointSourceMomentRangeError("joint source statistic exceeds the finite numeric range")
    return {name: float(value) for name, value in values.items()}


def joint_source_moments(history: TemkinPlanarReleaseHistory, *, quadrature_order: int,
                         include_component_records: bool = False):
    """Return joint emitted-population moments and original-bin conditionals.

    The conditional time law is integrated at finite supplied probability
    quadrature order. Its representation flags and log-moment qualifications
    are retained; these are not physical uncertainty or certified error bounds.
    Empty populations have zero emitted number and null normalized statistics.
    All covariances use centered variables in local coordinates, never a
    subtraction of large raw moments. The width-square reference is the center
    of the entire original face for every conditional and source hypothesis.
    """
    if not isinstance(history, TemkinPlanarReleaseHistory):
        raise ValueError("history must be a distinct TemkinPlanarReleaseHistory")
    if (isinstance(quadrature_order, (bool, np.bool_)) or not isinstance(quadrature_order, Integral)
            or not 2 <= quadrature_order <= MAX_QUADRATURE_ORDER):
        raise ValueError(f"quadrature_order must be an integer in [2,{MAX_QUADRATURE_ORDER}]")
    if not isinstance(include_component_records, (bool, np.bool_)):
        raise ValueError("include_component_records must be a boolean")
    order = int(quadrature_order)
    history = TemkinPlanarReleaseHistory(**{field.name: getattr(history, field.name)
                                          for field in fields(TemkinPlanarReleaseHistory)})
    x0, z0 = float(history.x_edges_m[0]), float(history.z_edges_m[0])
    full_width = float(history.x_edges_m[-1] - history.x_edges_m[0])
    x_face_center = _sum((.5 * x0, .5 * float(history.x_edges_m[-1])))
    time_origin = float(history.time_edges_s[0])
    dx, dz = np.diff(history.x_edges_m), np.diff(history.z_edges_m)
    x_centers = np.asarray([_sum((.5 * float(left), .5 * float(right), -.5 * x0,
                                 -.5 * float(history.x_edges_m[-1])))
                            for left, right in zip(history.x_edges_m[:-1], history.x_edges_m[1:])])
    x_geometry = [(*_x_geometry(float(center), float(width)), _sum((.5 * float(left), .5 * float(right))))
                  for center, width, left, right in zip(x_centers, dx, history.x_edges_m[:-1], history.x_edges_m[1:])]
    full_geometry = (*_x_geometry(0., full_width), x_face_center)
    area = history.patch_areas_m2
    result, component_records, underflows, subnormal_energies = [], [], [], []
    positive_components, max_time_tolerance = 0, 0.
    for p, pool in enumerate(history.reservoirs):
        keys = ("number", "interval", "width", "time_start", "time_offset", "time_variance", "energy", "x_first_half", "x_last_half",
                "x_variance", "square_mean", "square_variance", "z_first_half", "z_last_half", "z_variance")
        lists = {key: [] for key in keys}
        for interval, current, width in np.argwhere(history.released_m2[:, p] > 0):
            interval, current, width = int(interval), int(current), int(width)
            count = _product(float(history.released_m2[interval, p, current, width]), float(area[current, width]))
            if count == 0:
                raise JointSourceMomentRangeError("positive emitted patch population is below representable range")
            temperature = float(history.temperature_k[interval, current, width])
            expected_energy = _product(2., KB, temperature)
            if expected_energy == 0:
                raise JointSourceMomentRangeError("positive per-particle expected thermal energy is below representable range")
            emission = FrozenTemkinEmission(pool, float(history.coverage_m2[interval, p, current, width]), temperature,
                                            float(history.time_edges_s[interval + 1] - history.time_edges_s[interval]))
            timing = conditional_time_moments(emission, quadrature_order=order)
            timing_record = asdict(timing)
            indices = dict(reservoir_index=p, interval_index=interval, current_index=current, width_index=width)
            if timing.mean_underflowed or timing.second_raw_moment_underflowed or timing.variance_underflowed:
                underflows.append({**indices, "timing": timing_record})
            if expected_energy < np.finfo(float).tiny:
                subnormal_energies.append({**indices, "expected_translational_energy_j": expected_energy,
                                          "log_expected_translational_energy_j": math.fsum((math.log(2.), math.log(KB), math.log(temperature)))})
            if include_component_records:
                component_records.append({**indices, "released_number": count, "temperature_k": temperature,
                                          "expected_translational_energy_j": expected_energy, "timing": timing_record})
            positive_components += 1
            max_time_tolerance = max(max_time_tolerance, timing.numerical_log_tolerance)
            center, variance, square, square_variance, _ = x_geometry[width]
            values = (count, interval, width, float(history.time_edges_s[interval]), timing.mean_offset_s,
                      timing.variance_offset_s2, expected_energy, .5 * float(history.x_edges_m[width]), .5 * float(history.x_edges_m[width + 1]), variance, square, square_variance,
                      .5 * float(history.z_edges_m[current]), .5 * float(history.z_edges_m[current + 1]),
                      _signed_product(float(dz[current]), float(dz[current]), 1 / 12))
            for key, value in zip(keys, values):
                lists[key].append(value)
        batch = {key: np.asarray(values, dtype=int if key in {"interval", "width"} else float)
                 for key, values in lists.items()}
        all_components = np.ones(len(batch["number"]), dtype=bool)
        common = dict(x_face_center=x_face_center, x_face_left=x0, x_face_right=float(history.x_edges_m[-1]),
                      z_origin=z0, time_origin=time_origin)
        resolved = _statistics(batch, all_components, **common)
        uniform = _statistics(batch, all_components, uniform_x_geometry=full_geometry, **common)
        conditional_width = []
        for width in range(len(dx)):
            resolved_bin = _statistics(batch, batch["width"] == width, **common)
            uniform_bin = _statistics(batch, all_components, uniform_x_geometry=x_geometry[width],
                                      emitted_number=_count_fraction(resolved["released_number"], float(dx[width]), full_width), **common)
            conditional_width.append({"width_index": width, "x_left_m": float(history.x_edges_m[width]),
                                      "x_right_m": float(history.x_edges_m[width + 1]),
                                      "resolved": resolved_bin, "uniform_width": uniform_bin})
        conditional_interval = []
        for interval in range(len(history.time_edges_s) - 1):
            selected = batch["interval"] == interval
            conditional_interval.append({"interval_index": interval, "time_start_s": float(history.time_edges_s[interval]),
                                         "time_stop_s": float(history.time_edges_s[interval + 1]),
                                         "resolved": _statistics(batch, selected, **common),
                                         "uniform_width": _statistics(batch, selected, uniform_x_geometry=full_geometry, **common)})
        result.append({"name": pool.name, "state_id": pool.state_id, "released_number": resolved["released_number"],
                       "resolved": resolved, "uniform_width": uniform,
                       "conditional_width": conditional_width, "conditional_interval": conditional_interval})
    report = {
        "schema": SCHEMA, "units": dict(UNITS), "reservoirs": result,
        "coordinate_reference": {"x_face_center_m": x_face_center,
                                 "x_face_center_endpoint_parts_m": [.5 * x0, .5 * float(history.x_edges_m[-1])],
                                 "time_origin_s": time_origin, "z_origin_m": z0},
        "source_history_schema": history.schema, "law_identity": history.law_identity,
        "time_identity": TIME_IDENTITY, "quadrature_identity": QUADRATURE_IDENTITY, "quadrature_order": order,
        "coverage_model": history.coverage_model, "source_time_model": history.time_model,
        "surface_attributes": plain_metadata(history.surface_attributes), "parameter_provenance": history.parameter_provenance,
        "diagnostics": {"positive_components": positive_components, "component_time_underflows": underflows,
                        "maximum_component_numerical_log_tolerance": max_time_tolerance,
                        "subnormal_expected_translational_energy_components": subnormal_energies},
        "qualifications": {
            "population": "emitted population only; reservoirs remain independent and are never pooled into shared coverage",
            "source_position_model": "stationary independent uniform rectangle within each local source patch",
            "width_square_reference": "squared distance from the center of the entire original face, also within conditional bins",
            "width_null_model": "independent width redistribution preserving each original time/current/state/temperature mixture",
            "conditional_width_null": "each original width bin receives its width fraction of the full original mixture",
            "expected_energy_model": "2*kB*T per particle: conditional mean outgoing translation; variance/covariance of these local means only",
            "excluded_energy_fluctuations": "no within-Maxwellian kinetic-energy scatter, sampled packet energy, binding/internal energy or total wall debit",
            "time_model": "law-specific frozen-temperature conditional emitted-population time quadrature; finite order, not a certified bound",
            "representation": "time log moments and underflow flags retained; rounded SI zeros do not establish exact zero or agreement",
            "minimum_positive_source_probability": MIN_POSITIVE_SOURCE_PROBABILITY,
            "statistic_range": "nonzero geometry or weighted-statistic products rounding to zero raise; no physical floor is imposed",
            "spatial_statistics": "analytic within-patch second and fourth moments; centered covariances assembled in relative coordinates",
            "centering": "component start/local-time and spatial endpoint halves kept separate until each population is anchored; covariances retain centered time parts; reported absolute means may round at their origin",
            "scope": "no gap evolution, spectral inference, feedback, wall motion or actual FLASH execution",
        },
    }
    if include_component_records:
        report["component_records"] = component_records
    return report
