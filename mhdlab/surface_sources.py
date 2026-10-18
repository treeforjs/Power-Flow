"""Conservative one-dimensional surface source comparisons in SI coordinates.

Values are bin averages per unit surface area, including time-integrated release
in m^-2. Width integrals are per unit length in the invariant surface direction.
Every other array index (time, reservoir, component) is handled independently.
These operations neither advance desorption nor assign particle internal states,
velocities, collisions, or radiation. No source is inferred from temperature.
"""

from __future__ import annotations

from collections.abc import Collection, Mapping
import operator
from typing import Any

import numpy as np


def _contains_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return True
    if isinstance(value, np.ndarray):
        return value.dtype.kind == "b" or (
            value.dtype.kind == "O" and any(_contains_bool(item) for item in value.flat)
        )
    if isinstance(value, (list, tuple)):
        return any(_contains_bool(item) for item in value)
    return False


def _real_array(value: Any, name: str, *, nonnegative: bool = False) -> np.ndarray:
    if _contains_bool(value) or np.iscomplexobj(value):
        raise ValueError(f"{name} must contain real numbers, without booleans or complex values")
    try:
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf":
            raise ValueError(f"{name} must contain numeric real values")
        with np.errstate(over="raise", invalid="raise"):
            result = np.array(raw, dtype=np.float64, copy=True)
    except (TypeError, ValueError, OverflowError, FloatingPointError) as exc:
        raise ValueError(f"{name} must contain finite real numbers") from exc
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    if nonnegative and np.any(result < 0.0):
        raise ValueError(f"{name} must be nonnegative")
    return result


def _edges(value: Any, name: str) -> tuple[np.ndarray, np.ndarray]:
    edges = _real_array(value, name)
    if edges.ndim != 1 or edges.size < 2 or np.any(edges[1:] <= edges[:-1]):
        raise ValueError(f"{name} must contain at least two strictly increasing edges")
    try:
        with np.errstate(over="raise", invalid="raise"):
            widths = np.diff(edges)
            span = edges[-1] - edges[0]
    except FloatingPointError as exc:
        raise ValueError(f"{name} must span a finite width") from exc
    if not np.isfinite(span) or not np.all(np.isfinite(widths)):
        raise ValueError(f"{name} must span a finite width")
    return edges, widths


def _source(value: Any, widths: np.ndarray, axis: int, *, allow_signed: bool = False) -> tuple[np.ndarray, int]:
    if not isinstance(allow_signed, (bool, np.bool_)):
        raise ValueError("allow_signed must be a boolean")
    array = _real_array(value, "values_per_area", nonnegative=not allow_signed)
    if isinstance(axis, (bool, np.bool_)):
        raise ValueError("axis must be an integer")
    try:
        selected_axis = operator.index(axis)
    except TypeError as exc:
        raise ValueError("axis must be an integer") from exc
    if array.ndim == 0 or not -array.ndim <= selected_axis < array.ndim:
        raise ValueError("axis must identify an existing spatial bin dimension")
    selected_axis %= array.ndim
    if array.shape[selected_axis] != widths.size:
        raise ValueError("spatial bin count must equal len(x_edges_m)-1")
    return np.moveaxis(array, selected_axis, -1), selected_axis


def _finite_result(value: np.ndarray) -> np.ndarray:
    if not np.all(np.isfinite(value)):
        raise ValueError("source operation exceeded the finite numeric range")
    return value


def conservative_remap(
    values_per_area: Any,
    source_x_edges_m: Any,
    target_x_edges_m: Any,
    *,
    axis: int = -1,
    allow_signed: bool = False,
) -> np.ndarray:
    """Remap piecewise-constant bin averages by exact interval overlap.

    Arbitrary leading and trailing dimensions are supported via ``axis``. Each
    target average is sum(source_average * overlap_width)/target_width, so the
    width integral is conserved separately for every other-index combination.
    Endpoints must match exactly; there is no extrapolation or domain trimming.
    Signed fields require explicit ``allow_signed=True``; physical number and
    energy source fields remain nonnegative by default.
    """
    source_edges, source_widths = _edges(source_x_edges_m, "source_x_edges_m")
    target_edges, target_widths = _edges(target_x_edges_m, "target_x_edges_m")
    if source_edges[0] != target_edges[0] or source_edges[-1] != target_edges[-1]:
        raise ValueError("source and target must have exactly the same width-domain endpoints")
    source, selected_axis = _source(values_per_area, source_widths, axis, allow_signed=allow_signed)
    target = np.zeros((*source.shape[:-1], target_widths.size), dtype=np.float64)
    old_bin = new_bin = 0
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            while old_bin < source_widths.size and new_bin < target_widths.size:
                overlap = min(source_edges[old_bin + 1], target_edges[new_bin + 1]) - max(
                    source_edges[old_bin], target_edges[new_bin]
                )
                if overlap > 0.0:
                    target[..., new_bin] += source[..., old_bin] * (overlap / target_widths[new_bin])
                old_right, new_right = source_edges[old_bin + 1], target_edges[new_bin + 1]
                if old_right <= new_right:
                    old_bin += 1
                if new_right <= old_right:
                    new_bin += 1
    except FloatingPointError as exc:
        raise ValueError("source remap exceeded the finite numeric range") from exc
    return np.moveaxis(_finite_result(target), -1, selected_axis).copy()


def _uniform(values: Any, edges: np.ndarray, widths: np.ndarray, axis: int, *, allow_signed: bool) -> np.ndarray:
    source, selected_axis = _source(values, widths, axis, allow_signed=allow_signed)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            mean = np.sum(source * (widths / (edges[-1] - edges[0])), axis=-1, keepdims=True)
    except FloatingPointError as exc:
        raise ValueError("uniform source exceeded the finite numeric range") from exc
    result = np.broadcast_to(_finite_result(mean), source.shape).copy()
    return np.moveaxis(result, -1, selected_axis).copy()


def matched_uniform_source(released_m2: Any, x_edges_m: Any, *, axis: int = -1) -> np.ndarray:
    """Return a number-only uniform null preserving each width-integrated release.

    ``released_m2`` is time-integrated release per unit area in each spatial bin,
    not the number already integrated over that bin's area. The output has the
    same shape and a constant spatial value for each time/reservoir combination.
    Only supplied number integrals are constrained; energy and velocity/internal
    state moments are not implicitly preserved by a number-only comparison.
    """
    edges, widths = _edges(x_edges_m, "x_edges_m")
    return _uniform(released_m2, edges, widths, axis, allow_signed=False)


def matched_uniform_moments(
    moments_per_area: Mapping[str, Any],
    x_edges_m: Any,
    *,
    axis: int = -1,
    allow_signed: Collection[str] = (),
) -> dict[str, np.ndarray]:
    """Uniformize explicitly supplied moments, conserving each separate integral.

    Number, energy, and named momentum components can be supplied independently.
    Only names explicitly included in ``allow_signed`` may contain negative
    entries; unknown names are errors. Each array uses the same specified axis
    and may have its own nonspatial dimensions. This moment-constrained null is
    distinct from a number-only null. It constructs moment arrays, not a physical
    particle randomization: joint-distribution realizability and actual sampling
    remain the caller's responsibility.
    """
    if not isinstance(moments_per_area, Mapping) or not moments_per_area:
        raise ValueError("moments_per_area must be a nonempty mapping of named moments")
    if any(not isinstance(name, str) or not name.strip() for name in moments_per_area):
        raise ValueError("moment names must be nonempty strings")
    if isinstance(allow_signed, (str, bytes)) or not isinstance(allow_signed, Collection):
        raise ValueError("allow_signed must be a collection of supplied moment names")
    if any(not isinstance(name, str) or name not in moments_per_area for name in allow_signed):
        raise ValueError("allow_signed contains an unknown moment name")
    signed_names = set(allow_signed)
    edges, widths = _edges(x_edges_m, "x_edges_m")
    return {
        name: _uniform(values, edges, widths, axis, allow_signed=name in signed_names)
        for name, values in moments_per_area.items()
    }


def source_contrast(
    values_per_area: Any,
    x_edges_m: Any,
    *,
    axis: int = -1,
    edge_fraction: float = 0.1,
    center_fraction: float = 0.2,
) -> dict[str, np.ndarray | np.ma.MaskedArray]:
    """Overlap-weighted nonnegative source means over edges and the center.

    Each outer interval occupies ``edge_fraction`` of the full width; the
    centered interval occupies ``center_fraction``. Intervals may touch but not
    overlap. Results omit the spatial axis. ``edge_to_center_ratio`` is a masked
    array with an explicit mask wherever the center mean is zero, including the
    all-zero source limit. No infinity or assumed zero ratio is returned.
    """
    edges, widths = _edges(x_edges_m, "x_edges_m")
    source, _ = _source(values_per_area, widths, axis)
    edge = _real_array(edge_fraction, "edge_fraction")
    center = _real_array(center_fraction, "center_fraction")
    if edge.ndim != 0 or center.ndim != 0 or edge <= 0.0 or center <= 0.0:
        raise ValueError("region fractions must be positive finite real scalars")
    edge_value, center_value = float(edge), float(center)
    if 2.0 * edge_value + center_value > 1.0:
        raise ValueError("edge and center intervals must not overlap or exceed the width")
    start, stop = edges[[0, -1]]
    width = stop - start
    midpoint = start + 0.5 * width

    def mean(lo: float, hi: float) -> np.ndarray:
        overlap = np.maximum(0.0, np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo))
        covered_width = float(np.sum(overlap))
        if covered_width <= 0.0:
            raise ValueError("region is too small to resolve numerically on this x axis")
        return np.sum(source * (overlap / covered_width), axis=-1)

    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            left = mean(start, start + edge_value * width)
            right = mean(stop - edge_value * width, stop)
            middle = mean(midpoint - 0.5 * center_value * width, midpoint + 0.5 * center_value * width)
            combined = 0.5 * left + 0.5 * right
            undefined = middle == 0.0
            ratio = np.divide(combined, middle, out=np.zeros_like(middle), where=~undefined)
            result = {
                "left_edge_mean": left,
                "right_edge_mean": right,
                "combined_edge_mean": combined,
                "center_mean": middle,
                "edge_minus_center": combined - middle,
                "edge_to_center_ratio": np.ma.array(ratio, mask=undefined),
            }
    except FloatingPointError as exc:
        raise ValueError("source contrast exceeded the finite numeric range") from exc
    for value in result.values():
        _finite_result(np.ma.getdata(value))
    return result
