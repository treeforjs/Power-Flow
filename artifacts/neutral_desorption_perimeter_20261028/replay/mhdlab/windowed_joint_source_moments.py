"""Count-weighted conditional moments of selected planar emission windows.

Exact rational arithmetic combines the represented component counts, local
timing estimates, rectangle geometry and thermal moments. Absolute origins are
never added before centering. The matched hypothesis replaces only the width
law with a uniform law over the original face. Quadrature resolution remains a
separate numerical estimate, rather than an uncertainty or certified bound.
"""
from __future__ import annotations

from dataclasses import asdict, fields
from fractions import Fraction
import math
from numbers import Integral

import numpy as np

from .constants import KB
from .emission_time_moments import MAX_QUADRATURE_ORDER
from .temkin_release import TemkinReservoir
from .temkin_source_window import FrozenTemkinSourceWindow
from .temkin_windowed_source import TemkinPlanarSourceWindows
from .window_time_moments import conditional_window_time_moments


SCHEMA = "mhdlab-windowed-joint-source-moments-v1"
_ZERO = Fraction(0)
_EPS = float(np.finfo(float).eps)
UNITS = {
    "released_number": "parent",
    "mean_time_s": "s", "mean_time_offset_s": "s", "variance_time_s2": "s^2",
    "mean_x_m": "m", "variance_x_m2": "m^2",
    "mean_z_m": "m", "variance_z_m2": "m^2",
    "mean_width_centered_square_m2": "m^2", "variance_width_centered_square_m4": "m^4",
    "mean_expected_translational_energy_j": "J",
    "variance_expected_translational_energy_j2": "J^2",
    "mean_conditional_kinetic_energy_variance_j2": "J^2",
    "variance_kinetic_energy_j2": "J^2", "second_raw_kinetic_energy_j2": "J^2",
    "cov_time_x_s_m": "s m", "cov_time_z_s_m": "s m",
    "cov_time_expected_translational_energy_s_j": "s J",
    "cov_expected_translational_energy_x_j_m": "J m",
    "cov_expected_translational_energy_z_j_m": "J m",
    "cov_time_width_centered_square_s_m2": "s m^2",
    "cov_expected_translational_energy_width_centered_square_j_m2": "J m^2",
}


class WindowedJointMomentRangeError(ValueError):
    """A positive component or requested SI moment is unrepresentable."""


def _project(value, name, tolerance):
    """Project a signed exact rational without hiding underflow/quantization."""
    try:
        result = float(value)
    except OverflowError as error:
        raise WindowedJointMomentRangeError(f"{name} exceeds finite float64 range") from error
    if not math.isfinite(result) or (value and result == 0.):
        raise WindowedJointMomentRangeError(f"{name} is nonzero outside supported finite float64 range")
    if value and abs(Fraction(result)-value) > abs(value)*Fraction(max(tolerance, 256*_EPS)):
        raise WindowedJointMomentRangeError(f"{name} projection has unresolved float64 quantization")
    return result


def _validated_source(source):
    if not isinstance(source, TemkinPlanarSourceWindows):
        raise ValueError("source must be a TemkinPlanarSourceWindows record")
    try:
        source = source.validated_copy()
        for component in source.components:
            supplied = component.source_window
            parameters = {item.name: getattr(supplied, item.name)
                          for item in fields(FrozenTemkinSourceWindow) if item.init}
            parameters["reservoir"] = TemkinReservoir(**asdict(supplied.reservoir))
            replay = FrozenTemkinSourceWindow(**parameters)
            if replay.manifest() != supplied.manifest():
                raise ValueError("source-window quantities or reduced timing were modified")
        return source
    except (AttributeError, TypeError) as error:
        raise ValueError("source-window association is incomplete or modified") from error


def _rectangle(bounds, face_center):
    left, right = map(Fraction, bounds)
    mean = (left+right)/2
    width = right-left
    variance = width**2/12
    centered = mean-face_center
    return mean, variance, centered**2+variance, centered**2*width**2/3+width**4/180


def _diagnostic_positive(value):
    """Retain an optional higher moment without restricting required moments."""
    try:
        projected = float(value)
    except OverflowError:
        projected = None
    return {"value": projected, "log_value": math.log(value.numerator)-math.log(value.denominator),
            "underflowed": projected == 0., "overflowed": projected is None}


def _time_anchor(records, means, total, tolerance):
    # Count selection is independent of quadrature order. The base contains no
    # local quadrature mean, so order refinement compares the same offset.
    dominant = max(range(len(records)), key=lambda index: records[index]["number"])
    record = records[dominant]
    origin = record["component"].interval_origin_s
    start = record["component"].source_window.start_s
    base = Fraction(origin)+Fraction(start)
    mean = sum((item["number"]*value for item, value in zip(records, means)), _ZERO)/total
    offset = mean-base
    projected = _project(mean, "mean_time_s", tolerance)
    projected_offset = _project(offset, "mean_time_offset_s", tolerance)
    error = Fraction(projected)-mean
    qualified = not error or (offset != 0 and abs(error) <= abs(offset)*Fraction(tolerance))
    absolute_error = float(abs(error))
    return projected, projected_offset, {
        "component_index": record["index"],
        "interval_origin_s": origin, "window_start_offset_s": start,
        "mean_time_offset_s": projected_offset,
        "best_effort_mean_time_s": projected,
        "global_mean_time_precision_qualified": bool(qualified),
        "global_mean_time_absolute_rounding_error_s": absolute_error,
        "global_mean_time_rounding_error_underflowed": bool(error and absolute_error == 0.),
        "offset_definition": "mean time minus exact sum of the retained interval origin and window start",
    }


def _summarize(records, face_bounds, tolerance, *, matched):
    if not records:
        return {name: (0. if name == "released_number" else None) for name in UNITS}, None
    total = sum((record["number"] for record in records), _ZERO)
    center = sum(map(Fraction, face_bounds), _ZERO)/2
    xmoments = [_rectangle(face_bounds if matched else record["component"].x_bounds_m, center)
                for record in records]
    zmoments = [_rectangle(record["component"].z_bounds_m, _ZERO) for record in records]
    columns = {
        "time": [record["time"] for record in records],
        "x": [moment[0] for moment in xmoments],
        "z": [moment[0] for moment in zmoments],
        "width_square": [moment[2] for moment in xmoments],
        "energy": [record["energy"] for record in records],
    }
    means = {name: sum((record["number"]*value for record, value in zip(records, values)), _ZERO)/total
             for name, values in columns.items()}
    def average(values):
        return sum((record["number"]*value for record, value in zip(records, values)), _ZERO)/total
    def covariance(left, right):
        return average([(a-means[left])*(b-means[right])
                        for a, b in zip(columns[left], columns[right])])
    within_energy = average([record["energy_variance"] for record in records])
    between_energy = covariance("energy", "energy")
    mean_time, mean_offset, anchor = _time_anchor(records, columns["time"], total, tolerance)
    exact = {
        "released_number": total,
        "variance_time_s2": covariance("time", "time")+average([record["time_variance"] for record in records]),
        "mean_x_m": means["x"], "variance_x_m2": covariance("x", "x")+average([m[1] for m in xmoments]),
        "mean_z_m": means["z"], "variance_z_m2": covariance("z", "z")+average([m[1] for m in zmoments]),
        "mean_width_centered_square_m2": means["width_square"],
        "variance_width_centered_square_m4": covariance("width_square", "width_square")+average([m[3] for m in xmoments]),
        "mean_expected_translational_energy_j": means["energy"],
        "variance_expected_translational_energy_j2": between_energy,
        "mean_conditional_kinetic_energy_variance_j2": within_energy,
        "variance_kinetic_energy_j2": within_energy+between_energy,
        "second_raw_kinetic_energy_j2": within_energy+between_energy+means["energy"]**2,
        "cov_time_x_s_m": covariance("time", "x"),
        "cov_time_z_s_m": covariance("time", "z"),
        "cov_time_expected_translational_energy_s_j": covariance("time", "energy"),
        "cov_expected_translational_energy_x_j_m": covariance("energy", "x"),
        "cov_expected_translational_energy_z_j_m": covariance("energy", "z"),
        "cov_time_width_centered_square_s_m2": covariance("time", "width_square"),
        "cov_expected_translational_energy_width_centered_square_j_m2": covariance("energy", "width_square"),
    }
    values = {name: _project(value, name, tolerance) for name, value in exact.items()}
    values.update(mean_time_s=mean_time, mean_time_offset_s=mean_offset)
    return {name: values[name] for name in UNITS}, anchor


def windowed_joint_source_moments(source, *, quadrature_order, include_component_records=True):
    """Return detached resolved and release-matched uniform-width statistics.

    Positive extensive counts are mandatory. Primitive conditional local timing
    moments and normalized SI output must be representable; unsupported ranges
    reject explicitly, while the timing API retains its independent log fields.
    Areal underflow alone does not discard a represented extensive component.
    """
    if (np.ma.isMaskedArray(quadrature_order) or isinstance(quadrature_order, (bool, np.bool_))
            or not isinstance(quadrature_order, Integral) or not 2 <= quadrature_order <= MAX_QUADRATURE_ORDER):
        raise ValueError(f"quadrature_order must be an integer in [2,{MAX_QUADRATURE_ORDER}]")
    if type(include_component_records) is not bool:
        raise ValueError("include_component_records must be boolean")
    source = _validated_source(source)
    tolerance = source.log_mass_tolerance
    records = []
    for index, component in enumerate(source.components):
        window = component.source_window
        quantity = window.manifest()["quantities"]["released_number"]
        if not quantity["qualified"] or quantity["value"] is None or quantity["value"] <= 0 or quantity["exact_zero"]:
            raise WindowedJointMomentRangeError("every positive component requires a qualified extensive released_number")
        timing = conditional_window_time_moments(window.time_law, quadrature_order=int(quadrature_order))
        if timing.mean_underflowed or timing.variance_underflowed or timing.second_raw_moment_underflowed:
            raise WindowedJointMomentRangeError("positive primitive local timing moments are unrepresented; inspect their independent log moments")
        # Check coarse SI projections against their retained log estimates too.
        for name, log_name in (("mean_local_delay_s", "log_mean_local_delay_s"),
                               ("variance_local_delay_s2", "log_variance_local_delay_s2"),
                               ("second_raw_moment_s2", "log_second_raw_moment_s2")):
            if abs(math.log(getattr(timing, name))-getattr(timing, log_name)) > tolerance:
                raise WindowedJointMomentRangeError("primitive local timing float64 projection is unresolved")
        kbt = Fraction(KB)*Fraction(window.temperature_k)
        records.append({"index": index, "component": component, "number": Fraction(quantity["value"]),
                        "timing": timing, "time": Fraction(component.interval_origin_s)+Fraction(window.start_s)+Fraction(timing.mean_local_delay_s),
                        "time_variance": Fraction(timing.variance_local_delay_s2),
                        "energy": 2*kbt, "energy_variance": 2*kbt**2,
                        "energy_second": 6*kbt**2, "energy_fourth": 120*kbt**4})
    bounds = (float(source.x_edges_m[0]), float(source.x_edges_m[-1]))
    def group(selected):
        resolved, anchor = _summarize(selected, bounds, tolerance, matched=False)
        uniform, uniform_anchor = _summarize(selected, bounds, tolerance, matched=True)
        return {"resolved": resolved, "matched_uniform_width": uniform,
                "mean_time_anchor": anchor, "matched_uniform_width_mean_time_anchor": uniform_anchor}
    result = {"schema": SCHEMA, **group(records), "units": dict(UNITS),
              "quadrature_order": int(quadrature_order), "status": source.status,
              "provenance": source.provenance, "parameter_provenance": source.parameter_provenance,
              "source": source.manifest(), "host_budget_accepted": False,
              "per_selection": [{"selection_index": index, **selection.manifest(),
                                 **group([record for record in records if record["component"].selection_index == index])}
                                for index, selection in enumerate(source.selections)],
              "per_reservoir": [{"reservoir_index": index, "reservoir_name": reservoir.name, "state_id": reservoir.state_id,
                                 **group([record for record in records if record["component"].reservoir_index == index])}
                                for index, reservoir in enumerate(source.source_reservoirs)],
              "numerical_qualification": {
                  "mixture_arithmetic": "exact rational arithmetic on represented counts, geometry, thermal constants and conditional timing estimates; no probability cutoff",
                  "output_relative_projection_tolerance": max(tolerance, 256*_EPS),
                  "quadrature_interpretation": "finite-order local-delay quadrature estimates; not certified bounds or physical uncertainty",
                  "time_definition": "original interval origin plus source-window start plus conditional local delay; centered before SI projection",
                  "kinetic_energy_definition": "outward flux-Maxwellian Gamma(shape=2, scale=k_B T); conditional mean 2 k_B T and variance 2 (k_B T)^2",
              },
              "limits": ["Source-emission windows are distinct from photon recording gates.",
                         "Matched uniform width changes x only and preserves the complete conditional time/state/temperature/current/energy mixture.",
                         "No density, count, energy or moment floor; no positive component is omitted.",
                         "Required local timing and requested normalized SI moments must be represented; requesting component records additionally requires their standalone second moments to be represented. Optional fourth kinetic moments retain diagnostic logs outside the SI range.",
                         "Represented counts and timing estimates are numerical inputs, not exact physical source parameters.",
                         "No gap transport, molecular collisions, radiation, charged trajectories or wall-budget acceptance is evaluated."]}
    if include_component_records:
        total = sum((record["number"] for record in records), _ZERO)
        components = []
        for record in records:
            component, timing = record["component"], record["timing"]
            window = component.source_window
            probability = record["number"]/total
            projected_probability = float(probability)
            entry = {"component_index": record["index"],
                     **{name: getattr(component, name) for name in ("selection_index", "interval_index", "reservoir_index", "current_index", "width_index", "interval_origin_s")},
                     "window_start_offset_s": window.start_s, "window_stop_offset_s": window.stop_s,
                     "released_number": float(record["number"]), "temperature_k": window.temperature_k,
                     "reservoir_name": window.reservoir.name, "state_id": window.reservoir.state_id,
                     "x_bounds_m": list(component.x_bounds_m), "z_bounds_m": list(component.z_bounds_m),
                     "timing": asdict(timing), "probability": projected_probability,
                     "log_probability": math.log(probability.numerator)-math.log(probability.denominator),
                     "probability_underflowed": projected_probability == 0.,
                     "probability_projection_qualified": projected_probability > 0. and abs(Fraction(projected_probability)-probability) <= probability*Fraction(max(tolerance, 256*_EPS)),
                     **{name: _project(record[key], name, tolerance) for name, key in (
                         ("mean_expected_translational_energy_j", "energy"), ("variance_kinetic_energy_j2", "energy_variance"),
                         ("second_raw_kinetic_energy_j2", "energy_second"))},
                     "fourth_raw_kinetic_energy_j4": _diagnostic_positive(record["energy_fourth"]),
                     **group([record])}
            components.append(entry)
        result["components"] = components
    return result
