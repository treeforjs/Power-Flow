"""Publication input contract for an explicitly sourced electric-current history.

The trace is piecewise linear in aligned SI time. Its segment slopes describe
that interpolation, not an independent B-dot measurement or a differentiable raw
signal. This module does not infer voltage, fit a circuit, smooth data, or extend
the supplied time support. Digitized figures must be labeled reconstructed.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
import hashlib
import io
import math
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np


def _contains_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return True
    if isinstance(value, np.ndarray):
        return value.dtype.kind == "b" or (
            value.dtype.kind == "O" and any(_contains_bool(item) for item in value.flat)
        )
    if isinstance(value, (tuple, list)):
        return any(_contains_bool(item) for item in value)
    return False


def _real(value: Any, name: str) -> np.ndarray:
    if _contains_bool(value) or np.iscomplexobj(value) or np.ma.isMaskedArray(value):
        raise ValueError(f"{name} must contain unmasked real numbers, not booleans or complex values")
    try:
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf":
            raise ValueError(f"{name} must be numeric")
        with np.errstate(over="raise", invalid="raise"):
            result = np.array(raw, dtype=np.float64, copy=True)
    except (TypeError, ValueError, OverflowError, FloatingPointError) as exc:
        raise ValueError(f"{name} must contain finite real numbers") from exc
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _scalar(value: Any, name: str, *, positive: bool = False) -> float:
    result = _real(value, name)
    if result.ndim != 0 or (positive and result <= 0.0):
        raise ValueError(f"{name} must be a finite {'positive ' if positive else ''}real scalar")
    return float(result)


def _text(value: Any, name: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str) or (not allow_empty and not value.strip()):
        raise ValueError(f"{name} must be a {'string' if allow_empty else 'nonempty string'}")
    return value


@dataclass(frozen=True)
class DriveTrace:
    """One signed current waveform in amperes on strictly increasing SI times.

    ``time_s`` is already aligned; construction never applies a second offset.
    ``diagnostic_time_offset_s`` records the explicit offset used to align the
    original diagnostic time. Negative pretrigger times and signed currents are
    valid. Source kind is measured, reconstructed, or synthetic; reconstructed
    traces additionally require the original figure label in figure_reference.
    """

    time_s: np.ndarray
    current_a: np.ndarray
    source_kind: str
    shot_id: str
    source_reference: str
    diagnostic_time_offset_s: float
    figure_reference: str = ""
    csv_import_metadata: Mapping[str, str | float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        time = _real(self.time_s, "time_s")
        current = _real(self.current_a, "current_a")
        if time.ndim != 1 or time.size < 2 or np.any(time[1:] <= time[:-1]):
            raise ValueError("time_s must contain at least two strictly increasing times")
        if current.shape != time.shape:
            raise ValueError("current_a must have the same one-dimensional shape as time_s")
        try:
            with np.errstate(over="raise", invalid="raise"):
                span = time[-1] - time[0]
        except FloatingPointError as exc:
            raise ValueError("time_s must span a finite interval") from exc
        if not np.isfinite(span):
            raise ValueError("time_s must span a finite interval")
        for name in ("source_kind", "shot_id", "source_reference"):
            _text(getattr(self, name), name)
        if self.source_kind not in {"measured", "reconstructed", "synthetic"}:
            raise ValueError("source_kind must be measured, reconstructed, or synthetic")
        _text(self.figure_reference, "figure_reference", allow_empty=self.source_kind != "reconstructed")
        offset = _scalar(self.diagnostic_time_offset_s, "diagnostic_time_offset_s")
        if not isinstance(self.csv_import_metadata, Mapping):
            raise ValueError("csv_import_metadata must be a mapping")
        metadata = {}
        for key, value in self.csv_import_metadata.items():
            _text(key, "CSV metadata key")
            metadata[key] = value if isinstance(value, str) else _scalar(value, f"CSV metadata {key}")
        time.setflags(write=False)
        current.setflags(write=False)
        object.__setattr__(self, "time_s", time)
        object.__setattr__(self, "current_a", current)
        object.__setattr__(self, "diagnostic_time_offset_s", offset)
        object.__setattr__(self, "csv_import_metadata", MappingProxyType(metadata))

    def _query(self, time_s: Any) -> np.ndarray:
        query = _real(time_s, "query time_s")
        if query.ndim > 1 or query.size == 0:
            raise ValueError("query time_s must be a scalar or nonempty one-dimensional array")
        if np.any(query < self.time_s[0]) or np.any(query > self.time_s[-1]):
            raise ValueError("query lies outside trace support; extrapolation is forbidden")
        return query

    def current_at(self, time_s: Any) -> float | np.ndarray:
        """Interpolate signed current within the supplied time support only."""
        query = self._query(time_s)
        index = np.clip(np.searchsorted(self.time_s, query, side="right") - 1, 0, self.time_s.size - 2)
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                fraction = (query - self.time_s[index]) / (self.time_s[index + 1] - self.time_s[index])
                value = (1.0 - fraction) * self.current_a[index] + fraction * self.current_a[index + 1]
        except FloatingPointError as exc:
            raise ValueError("current interpolation exceeded the finite numeric range") from exc
        if not np.all(np.isfinite(value)):
            raise ValueError("current interpolation exceeded the finite numeric range")
        return value

    def slope_a_s(self, time_s: Any, *, knot_policy: str) -> float | np.ndarray:
        """Return the piecewise-constant interpolation slope, not raw dI/dt.

        At interior knots choose ``left`` or ``right`` explicitly, or ``raise``
        to reject those queries. At either endpoint the sole inward segment is
        used for all three policies. No differentiability across a knot is implied.
        """
        if not isinstance(knot_policy, str) or knot_policy not in {"left", "right", "raise"}:
            raise ValueError("knot_policy must be left, right, or raise")
        query = self._query(time_s)
        if knot_policy == "raise" and np.any(np.isin(query, self.time_s[1:-1])):
            raise ValueError("slope is undefined at an interior knot under knot_policy='raise'")
        side = "right" if knot_policy == "right" else "left"
        index = np.clip(np.searchsorted(self.time_s, query, side=side) - 1, 0, self.time_s.size - 2)
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                slopes = np.diff(self.current_a) / np.diff(self.time_s)
        except FloatingPointError as exc:
            raise ValueError("segment slope exceeded the finite numeric range") from exc
        return slopes[index]

    def _segments(self, start_s: Any, end_s: Any) -> tuple[np.ndarray, np.ndarray]:
        start = _scalar(start_s, "start_s")
        end = _scalar(end_s, "end_s")
        self._query([start, end])
        if end < start:
            raise ValueError("end_s must be greater than or equal to start_s")
        interior = self.time_s[(self.time_s > start) & (self.time_s < end)]
        times = np.concatenate(([start], interior, [end]))
        return np.diff(times), np.asarray(self.current_at(times))

    def integrate_current_a_s(self, start_s: Any, end_s: Any) -> float:
        """Exact integral of the piecewise-linear I(t), including partial bins."""
        duration, current = self._segments(start_s, end_s)
        try:
            with np.errstate(over="raise", invalid="raise"):
                terms = duration * (0.5 * current[:-1] + 0.5 * current[1:])
            result = math.fsum(float(value) for value in terms)
        except (FloatingPointError, OverflowError) as exc:
            raise ValueError("current integral exceeded the finite numeric range") from exc
        if not math.isfinite(result):
            raise ValueError("current integral exceeded the finite numeric range")
        return result

    def integrate_current_squared_a2_s(self, start_s: Any, end_s: Any) -> float:
        """Exact integral of I(t)^2, not a trapezoidal integral of squared samples."""
        duration, current = self._segments(start_s, end_s)
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                magnitude = np.maximum(np.abs(current[:-1]), np.abs(current[1:]))
                left = np.divide(current[:-1], magnitude, out=np.zeros_like(magnitude), where=magnitude > 0.0)
                right = np.divide(current[1:], magnitude, out=np.zeros_like(magnitude), where=magnitude > 0.0)
                shape = (left * left + left * right + right * right) / 3.0
                terms = ((duration * magnitude) * magnitude) * shape
            result = math.fsum(float(value) for value in terms)
        except (FloatingPointError, OverflowError) as exc:
            raise ValueError("squared-current integral exceeded the finite numeric range") from exc
        if not math.isfinite(result):
            raise ValueError("squared-current integral exceeded the finite numeric range")
        return result


def load_current_csv(
    path: str | Path,
    *,
    time_column: str,
    current_column: str,
    time_scale_to_s: float,
    current_scale_to_a: float,
    diagnostic_time_offset_s: float,
    source_kind: str,
    shot_id: str,
    source_reference: str,
    figure_reference: str = "",
) -> DriveTrace:
    """Load explicitly selected columns with declared SI conversion/alignment.

    aligned_time_s = CSV_time * time_scale_to_s + diagnostic_time_offset_s;
    current_a = CSV_current * current_scale_to_a. Unit scales are strictly positive;
    current polarity is represented by the data. No header-based unit guessing,
    time sorting, smoothing, baseline subtraction, or out-of-support padding occurs.
    The returned metadata records the exact imported CSV bytes' SHA-256 digest.
    """
    _text(time_column, "time_column")
    _text(current_column, "current_column")
    if time_column == current_column:
        raise ValueError("time_column and current_column must be different")
    time_scale = _scalar(time_scale_to_s, "time_scale_to_s", positive=True)
    current_scale = _scalar(current_scale_to_a, "current_scale_to_a", positive=True)
    offset = _scalar(diagnostic_time_offset_s, "diagnostic_time_offset_s")
    try:
        source_path = Path(path)
        content = source_path.read_bytes()
        text = content.decode("utf-8-sig")
    except (OSError, TypeError, UnicodeError) as exc:
        raise ValueError(f"cannot read current CSV: {exc}") from exc
    try:
        reader = csv.DictReader(io.StringIO(text), strict=True)
        names = reader.fieldnames
        if (names is None or any(not name.strip() for name in names)
                or len(names) != len(set(names))):
            raise ValueError("CSV requires nonempty, unique column headers")
        if time_column not in names or current_column not in names:
            raise ValueError("explicit time/current column is missing from CSV")
        times, currents = [], []
        for row in reader:
            if None in row:
                raise ValueError(f"CSV row {reader.line_num} contains extra columns")
            try:
                times.append(float(row[time_column]))
                currents.append(float(row[current_column]))
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(f"CSV row {reader.line_num} has an invalid selected value") from exc
        with np.errstate(over="raise", invalid="raise"):
            aligned_time = np.asarray(times, dtype=float) * time_scale + offset
            current = np.asarray(currents, dtype=float) * current_scale
    except (csv.Error, FloatingPointError) as exc:
        raise ValueError(f"invalid current CSV: {exc}") from exc
    return DriveTrace(
        aligned_time, current, source_kind, shot_id, source_reference, offset,
        figure_reference=figure_reference,
        csv_import_metadata={
            "path": str(source_path.resolve()),
            "sha256": hashlib.sha256(content).hexdigest(),
            "time_column": time_column,
            "current_column": current_column,
            "time_scale_to_s": time_scale,
            "current_scale_to_a": current_scale,
            "alignment": "aligned_time_s = CSV_time * time_scale_to_s + diagnostic_time_offset_s",
        },
    )
