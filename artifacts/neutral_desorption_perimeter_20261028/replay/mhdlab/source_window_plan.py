"""Map represented physical source-time windows to original interval offsets.

All intersections and projection errors use exact Fractions of the supplied
binary64 inputs. This is a coordinate/representation qualification, not an
estimate of physical timing uncertainty. Windows describe source emission,
not photon recording gates. Exact intersections are formed without truncation
or dropping support. Qualified nearest-float endpoint projection can move cuts
and close small exact gaps; every motion is recorded. No clamping or manual
repair is used to make a local selection fit, and projected overlaps reject.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections.abc import Sequence
from dataclasses import dataclass, field, fields
from fractions import Fraction
import json
import math
from numbers import Real

import numpy as np

from .temkin_windowed_source import MAX_SELECTIONS, SourceWindowSelection


SCHEMA = "mhdlab-source-window-plan-v1"
QUALIFICATION_IDENTITY = "exact_represented_fraction_intersections_relative_local_projection_v1"
MAX_TIME_EDGES = 100_001
MAX_WINDOWS = 10_000


class SourceWindowPlanRangeError(ValueError):
    """An exact positive intersection lacks an accepted local representation."""


def _real(value, name):
    if (np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_))
            or not isinstance(value, Real)):
        raise ValueError(f"{name} must be a finite real scalar, not boolean or masked")
    try:
        result = float(value)
    except (ValueError, OverflowError, TypeError) as error:
        raise ValueError(f"{name} must be finite and representable") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _sequence(value, name, limit, minimum=1):
    if (np.ma.isMaskedArray(value) or isinstance(value, (str, bytes))
            or not (isinstance(value, Sequence) or isinstance(value, np.ndarray))):
        raise ValueError(f"{name} must be a bounded finite sequence, not a generator or mask")
    if isinstance(value, np.ndarray) and value.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    count = len(value)
    if not minimum <= count <= limit:
        raise ValueError(f"{name} requires {minimum}..{limit} entries within its work limit")
    return tuple(value[index] for index in range(count))


def _fraction(value):
    return Fraction.from_float(value)


def _rational(value):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _error(value, duration=None):
    try:
        projected = float(value)
    except OverflowError as error:
        raise SourceWindowPlanRangeError("signed projection error exceeds finite diagnostic range") from error
    result = {"exact_s": _rational(value), "value_s": projected,
              "value_underflowed": value != 0 and projected == 0.}
    if duration is not None:
        relative = value/duration
        try:
            represented = float(relative)
        except OverflowError as error:
            raise SourceWindowPlanRangeError("relative projection error exceeds finite diagnostic range") from error
        result.update(relative_to_intersection=_rational(relative),
                      relative_value=represented,
                      relative_value_underflowed=relative != 0 and represented == 0.)
    return result


@dataclass(frozen=True, slots=True)
class SourceTimeWindow:
    start_s: float
    stop_s: float

    def __post_init__(self):
        object.__setattr__(self, "start_s", _real(self.start_s, "start_s"))
        object.__setattr__(self, "stop_s", _real(self.stop_s, "stop_s"))
        if not 0 <= self.start_s < self.stop_s:
            raise ValueError("source windows require finite 0<=start_s<stop_s")

    def manifest(self):
        # Replay public input records too, instead of blessing forged fields.
        rebuilt = SourceTimeWindow(self.start_s, self.stop_s)
        return {"start_s": rebuilt.start_s, "stop_s": rebuilt.stop_s}


def _inputs(time_edges_s, windows, relative_time_tolerance):
    edges = tuple(_real(value, "time_edges_s") for value in
                  _sequence(time_edges_s, "time_edges_s", MAX_TIME_EDGES, minimum=2))
    if edges[0] < 0 or any(right <= left for left, right in zip(edges[:-1], edges[1:])):
        raise ValueError("time_edges_s must be nonnegative and strictly increasing")
    supplied = _sequence(windows, "windows", MAX_WINDOWS)
    if any(not isinstance(value, SourceTimeWindow) for value in supplied):
        raise ValueError("windows must contain SourceTimeWindow records")
    windows = tuple(SourceTimeWindow(value.start_s, value.stop_s) for value in supplied)
    ordered = sorted(windows, key=lambda value: (value.start_s, value.stop_s))
    if any(right.start_s < left.stop_s for left, right in zip(ordered[:-1], ordered[1:])):
        raise ValueError("source windows overlap or duplicate; exact adjacency is allowed")
    if any(value.start_s < edges[0] or value.stop_s > edges[-1] for value in windows):
        raise ValueError("every source window must be fully supported; truncation is forbidden")
    tolerance = _real(relative_time_tolerance, "relative_time_tolerance")
    if tolerance <= 0:
        raise ValueError("relative_time_tolerance must be explicit and positive")
    return edges, windows, tolerance


def _same(left, right):
    if type(left) is not type(right):
        return False
    if isinstance(left, tuple):
        return len(left) == len(right) and all(_same(a, b) for a, b in zip(left, right))
    if isinstance(left, (SourceTimeWindow, SourceWindowSelection)):
        return all(_same(getattr(left, item.name), getattr(right, item.name)) for item in fields(type(right)))
    return left == right


def _build(edges, windows, tolerance):
    plans = []
    count = 0
    # Bound expansion before producing rational records or any selections.
    for index, window in enumerate(windows):
        first = bisect_right(edges, window.start_s)-1
        end = bisect_left(edges, window.stop_s)
        count += end-first
        if count > MAX_SELECTIONS:
            raise ValueError(f"source-window intersection count exceeds work limit {MAX_SELECTIONS}")
        plans.append((index, window, first, end))
    selections, window_indices, records, window_records = [], [], [], []
    tolerance_q = _fraction(tolerance)
    by_interval = {}
    for window_index, window, first, end in plans:
        window_start, window_stop = _fraction(window.start_s), _fraction(window.stop_s)
        exact_total = Fraction(0)
        projected_total = Fraction(0)
        canonical_total = Fraction(0)
        first_selection = len(selections)
        for interval in range(first, end):
            origin = _fraction(edges[interval])
            parent_stop = _fraction(edges[interval+1])
            a, b = max(window_start, origin), min(window_stop, parent_stop)
            duration = b-a
            if duration <= 0:
                raise SourceWindowPlanRangeError("interval intersection is not strictly positive")
            local_a, local_b = a-origin, b-origin
            start, stop = float(local_a), float(local_b)
            parent_duration = float(edges[interval+1]-edges[interval])
            if (not all(math.isfinite(value) for value in (start, stop, parent_duration))
                    or not 0 <= start < stop <= parent_duration):
                raise SourceWindowPlanRangeError("positive intersection collapses or lies outside canonical parent duration after projection")
            projected_a, projected_b = _fraction(start), _fraction(stop)
            projected_duration = projected_b-projected_a
            local_duration = float(stop-start)
            if not math.isfinite(local_duration) or local_duration <= 0:
                raise SourceWindowPlanRangeError("positive local duration is not representable")
            canonical_duration = _fraction(local_duration)
            parent_error = _fraction(parent_duration)-(parent_stop-origin)
            touches_parent_stop = b == parent_stop
            errors = {
                "local_start_projection_error": projected_a-local_a,
                "local_stop_projection_error": projected_b-local_b,
                "projected_offset_duration_error": projected_duration-duration,
                "local_duration_total_error": canonical_duration-duration,
            }
            if touches_parent_stop:
                errors["parent_duration_boundary_error"] = parent_error
            if any(abs(error) > tolerance_q*duration for error in errors.values()):
                raise SourceWindowPlanRangeError("local coordinate error exceeds relative_time_tolerance against positive intersection duration")
            selection = SourceWindowSelection(interval, start, stop)
            by_interval.setdefault(interval, []).append(selection)
            selections.append(selection)
            window_indices.append(window_index)
            global_start, global_stop = edges[interval]+start, edges[interval]+stop
            if not math.isfinite(global_start) or not math.isfinite(global_stop):
                raise SourceWindowPlanRangeError("global diagnostic projection is not finite")
            global_start_q, global_stop_q = _fraction(global_start), _fraction(global_stop)
            records.append({
                "selection_index": len(selections)-1, "window_index": window_index,
                "interval_index": interval, "interval_origin_s": edges[interval],
                "selection": selection.manifest(),
                "exact_global_start_s": _rational(a), "exact_global_stop_s": _rational(b),
                "exact_local_start_s": _rational(local_a), "exact_local_stop_s": _rational(local_b),
                "exact_intersection_duration_s": _rational(duration),
                "exact_projected_offset_duration_s": _rational(projected_duration),
                "canonical_local_duration_s": local_duration,
                "exact_parent_duration_s": _rational(parent_stop-origin),
                "canonical_parent_duration_s": parent_duration,
                "touches_parent_stop": touches_parent_stop,
                "parent_duration_error_qualified_against_intersection": touches_parent_stop,
                **{name: _error(error, duration) for name, error in errors.items()},
                "local_duration_subtraction_rounding_error": _error(canonical_duration-projected_duration, duration),
                "canonical_parent_duration_error": _error(parent_error, duration if touches_parent_stop else None),
                "global_projected_start_s": global_start, "global_projected_stop_s": global_stop,
                "global_start_projection_error": _error(global_start_q-a),
                "global_stop_projection_error": _error(global_stop_q-b),
                "global_start_addition_rounding_error": _error(global_start_q-(origin+projected_a)),
                "global_stop_addition_rounding_error": _error(global_stop_q-(origin+projected_b)),
                "precision_flags": {
                    "local_start_rounded": projected_a != local_a,
                    "local_stop_rounded": projected_b != local_b,
                    "canonical_local_duration_rounded": canonical_duration != projected_duration,
                    "canonical_parent_duration_rounded": parent_error != 0,
                    "global_start_addition_rounded": global_start_q != origin+projected_a,
                    "global_stop_addition_rounded": global_stop_q != origin+projected_b,
                    "global_start_differs_from_requested": global_start_q != a,
                    "global_stop_differs_from_requested": global_stop_q != b,
                },
            })
            exact_total += duration
            projected_total += projected_duration
            canonical_total += canonical_duration
        if exact_total != window_stop-window_start:
            raise SourceWindowPlanRangeError("exact intersections do not partition the supplied window")
        window_records.append({"window_index": window_index,
                               "selection_indices": list(range(first_selection, len(selections))),
                               "exact_window_duration_s": _rational(exact_total),
                               "exact_projected_offset_duration_sum_s": _rational(projected_total),
                               "exact_canonical_local_duration_sum_s": _rational(canonical_total),
                               "projected_offset_duration_sum_error": _error(projected_total-exact_total, exact_total),
                               "canonical_local_duration_sum_error": _error(canonical_total-exact_total, exact_total)})
    # Check the operational selections too. Rounding must not introduce an
    # overlap/duplication even if the original global windows were disjoint.
    for values in by_interval.values():
        ordered = sorted(values, key=lambda value: (value.start_s, value.stop_s))
        if any(right.start_s < left.stop_s for left, right in zip(ordered[:-1], ordered[1:])):
            raise SourceWindowPlanRangeError("projected local selections overlap or duplicate")
    manifest = {
        "schema": SCHEMA, "qualification_identity": QUALIFICATION_IDENTITY,
        "time_edges_s": list(edges), "windows": [value.manifest() for value in windows],
        "relative_time_tolerance": tolerance,
        "selections": [value.manifest() for value in selections],
        "window_indices": window_indices, "selection_records": records,
        "window_records": window_records, "units": "s",
        "input_interpretation": "exact represented binary64 endpoints; rational numerator/denominator strings",
        "global_projection_interpretation": "diagnostic projections only; preserve original interval origin and local offsets separately",
        "projection_limitations": "qualified nearest-float endpoints may move cuts or close exact gaps; no clipping, dropping, clamping or manual repair; geometric tolerance is not a source-mass or physical timing uncertainty bound",
        "source_emission_windows_are_photon_recording_gates": False,
        "physical_timing_uncertainty_estimated": False,
        "limits": {"max_time_edges": MAX_TIME_EDGES, "max_windows": MAX_WINDOWS,
                   "max_selections": MAX_SELECTIONS},
    }
    return tuple(selections), tuple(window_indices), json.dumps(manifest, allow_nan=False)


@dataclass(frozen=True, slots=True, eq=False)
class SourceWindowPlan:
    time_edges_s: tuple[float, ...]
    windows: tuple[SourceTimeWindow, ...]
    relative_time_tolerance: float
    selections: tuple[SourceWindowSelection, ...] = field(init=False)
    window_indices: tuple[int, ...] = field(init=False)
    _manifest_json: str = field(init=False, repr=False)

    def __post_init__(self):
        edges, windows, tolerance = _inputs(self.time_edges_s, self.windows, self.relative_time_tolerance)
        selections, indices, manifest = _build(edges, windows, tolerance)
        for name, value in (("time_edges_s", edges), ("windows", windows),
                            ("relative_time_tolerance", tolerance), ("selections", selections),
                            ("window_indices", indices), ("_manifest_json", manifest)):
            object.__setattr__(self, name, value)

    def validated_copy(self):
        rebuilt = SourceWindowPlan(self.time_edges_s, self.windows, self.relative_time_tolerance)
        if any(not _same(getattr(self, item.name), getattr(rebuilt, item.name)) for item in fields(SourceWindowPlan)):
            raise ValueError("source-window plan disagrees with exact constructor replay")
        return rebuilt

    def manifest(self):
        return json.loads(self.validated_copy()._manifest_json)


def plan_source_windows(time_edges_s, windows, *, relative_time_tolerance) -> SourceWindowPlan:
    """Intersect disjoint global source windows without clipping or repair.

    Errors are relative to each positive intersection duration, rather than a
    potentially enormous global coordinate. Canonical parent-duration rounding
    is always reported; it enters that intersection's tolerance test when the
    selection actually reaches the parent's stop. Interior selections are not
    rejected merely because a distant parent stop has a rounding discrepancy.
    """
    return SourceWindowPlan(time_edges_s, windows, relative_time_tolerance)
