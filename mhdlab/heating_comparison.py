"""Cell-average magnetic-heating verification and a fixed-SI lineout contract.

This module compares one depth lineout against the constant-property half-space
reference; it neither reads raw FLASH meshes nor runs a solver.  B and J can be
signed.  Joule power is an independent cell average of J**2/sigma, so it is
never reconstructed by squaring a cell-averaged current.  Temperature is not
compared: the ideal boundary field step has a singular adiabatic surface
temperature, and a cell average is not a centroid sample.

CSV files have a required JSON sidecar containing units, sampling semantics,
provenance and a digest of the exact CSV bytes.  The digest detects changes to
the pair; it is not an authentication mechanism for its scientific provenance.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import io
import json
import math
from numbers import Real
from pathlib import Path

import numpy as np

from .heating_benchmarks import ConstantConductivityHalfSpace


HEATING_PROFILE_SCHEMA = "mhdlab-heating-profile-v1"
HEATING_COMPARISON_SCHEMA = "mhdlab-heating-comparison-v1"
CSV_COLUMNS = ("depth_left_m", "depth_right_m", "bx_t", "jz_a_m2", "joule_w_m3")
PROFILE_UNITS = {
    "time_s": "s", "depth_left_m": "m", "depth_right_m": "m",
    "bx_t": "T", "jz_a_m2": "A m^-2", "joule_w_m3": "W m^-3",
}
COORDINATE_SYSTEM = "y=conductor_inward_depth;Bx=tangential;Jz=-dBx/dy/mu"
_FIELDS = (
    ("magnetic_field_t", "T"), ("current_density_a_m2", "A m^-2"),
    ("joule_power_w_m3", "W m^-3"),
)
_LOG_MAX = math.log(np.finfo(float).max)
_ASSUMPTIONS = (
    "stationary uniform conductor on the half-space y>=0",
    "constant electrical conductivity and magnetic permeability",
    "constant volumetric heat capacity; B/J/q comparisons do not use temperature",
    "initial magnetic field zero; tangential boundary field steps to a constant at t>0",
    "no fluid motion, displacement current, radiation or phase change",
    "no thermal conduction; temperature is not a comparison observable",
    "finite-volume B, J and Joule-power averages are integrated separately",
    "norms cover only the supplied depth interval, not the infinite half-space energy",
)


def _real(value, name: str, *, positive: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be finite real, not bool, complex or string")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite real") from exc
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return result


def _array(value, name: str) -> np.ndarray:
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot be a masked array")
    try:
        raw = np.asarray(value, dtype=object)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{name} must be a finite real 1-D array") from exc
    if raw.ndim != 1 or raw.size == 0:
        raise ValueError(f"{name} must be a nonempty 1-D array")
    array = np.array([_real(item, name) for item in raw], dtype=float)
    return np.frombuffer(array.tobytes(), dtype=float)


def _text(value, name: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str) or (not allow_empty and not value.strip()):
        raise ValueError(f"{name} must be a nonempty string" if not allow_empty else f"{name} must be a string")
    if any(ord(character) < 32 for character in value):
        raise ValueError(f"{name} cannot contain control characters")
    return value


@dataclass(frozen=True, slots=True, eq=False)
class HeatingProfile:
    """One finite-volume depth lineout in SI, with explicit origin labels.

    ``depth_edges_m`` is nonnegative inward conductor depth; the lineout may
    start below the surface.  All three fields have one value per cell and
    refer to the same strictly positive time after the boundary step.  A
    FLASH label requires a supplied version; labels do not establish that a
    real case satisfies the benchmark assumptions.
    """

    time_s: float
    depth_edges_m: np.ndarray
    magnetic_field_t: np.ndarray
    current_density_a_m2: np.ndarray
    joule_power_w_m3: np.ndarray
    source_kind: str
    run_id: str
    source_reference: str
    flash_version: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "time_s", _real(self.time_s, "time_s", positive=True))
        edges = _array(self.depth_edges_m, "depth_edges_m")
        if edges.size < 2 or edges[0] < 0 or np.any(edges[1:] <= edges[:-1]):
            raise ValueError("depth_edges_m must have at least two nonnegative strictly increasing edges")
        object.__setattr__(self, "depth_edges_m", edges)
        for name, _ in _FIELDS:
            value = _array(getattr(self, name), name)
            if value.size != edges.size - 1:
                raise ValueError(f"{name} must have one cell average per depth interval")
            if name == "joule_power_w_m3" and np.any(value < 0):
                raise ValueError("joule_power_w_m3 must be nonnegative")
            object.__setattr__(self, name, value)
        if self.source_kind not in ("flash", "synthetic") or not isinstance(self.source_kind, str):
            raise ValueError("source_kind must be 'flash' or 'synthetic'")
        _text(self.run_id, "run_id")
        _text(self.source_reference, "source_reference")
        _text(self.flash_version, "flash_version", allow_empty=self.source_kind != "flash")

    @property
    def metadata(self) -> dict:
        """Return an independent JSON-compatible provenance and units record."""
        return {
            "schema": HEATING_PROFILE_SCHEMA,
            "time_s": self.time_s, "source_kind": self.source_kind,
            "run_id": self.run_id, "source_reference": self.source_reference,
            "flash_version": self.flash_version,
            "coordinate_system": COORDINATE_SYSTEM,
            "cell_sampling": "cell_average", "units": dict(PROFILE_UNITS),
        }


def _log_sum(terms: list[float]) -> float:
    if not terms:
        return -math.inf
    scale = max(terms)
    return scale + math.log(math.fsum(math.exp(term - scale) for term in terms))


def _log_norms(log_values: list[float], widths: np.ndarray, span: float) -> dict[str, float]:
    # Keep width and magnitude in log space until after summation: normalizing
    # a narrow cell's width first can underflow even when width*value is finite.
    active = [(math.log(float(width)), value) for width, value in zip(widths, log_values)
              if value != -math.inf]
    if not active:
        return {"l1": -math.inf, "l2": -math.inf, "linf": -math.inf}
    upper = max(log_values)
    # These are normalized positive moments, hence both are <= Linf.  Enforce
    # that analytic bound after log subtraction so cancellation at maxfloat
    # cannot turn a representable constant profile into a spurious overflow.
    return {
        "l1": min(upper, _log_sum([width + value for width, value in active]) - math.log(span)),
        "l2": min(upper, 0.5 * (_log_sum([width + 2 * value for width, value in active]) - math.log(span))),
        "linf": upper,
    }


def _log_abs_difference(first: float, second: float) -> float:
    if first == second:
        return -math.inf
    if (first >= 0) == (second >= 0):
        return math.log(abs(first - second))
    # Opposite signs can overflow a direct subtraction even when a weighted
    # L1/L2 norm is representable.  Sum their magnitudes before exponentiating.
    high, low = max(abs(first), abs(second)), min(abs(first), abs(second))
    return math.log(high) + math.log1p(low / high)


def _finite_exp(value: float, label: str) -> float:
    if value == -math.inf:
        return 0.0
    if not math.isfinite(value) or value > _LOG_MAX:
        raise ValueError(f"{label} exceeds finite float64 range")
    try:
        return math.exp(value)
    except OverflowError as exc:
        raise ValueError(f"{label} exceeds finite float64 range") from exc


def compare_heating_profile(profile: HeatingProfile, benchmark: ConstantConductivityHalfSpace) -> dict:
    """Return depth-weighted absolute and relative B/J/q errors as JSON data.

    L1 is mean absolute error and L2 is root mean square error over the
    supplied interval, using exact cell widths.  Each relative norm divides
    by that norm of the analytic reference; an identically zero reference has
    ``None`` ratios.  Zero here means all evaluated float64 reference cells
    are zero; sufficiently small analytic tails can underflow in the reference.
    Log moments retain very wide numeric ranges.  A positive
    norm below float64 range rounds to zero, while ``reference_is_zero`` still
    distinguishes nonzero represented reference cells.  Unrepresentable large
    errors or ratios raise rather than emitting invalid JSON infinity.
    """
    if not isinstance(profile, HeatingProfile):
        raise ValueError("profile must be a HeatingProfile")
    if not isinstance(benchmark, ConstantConductivityHalfSpace):
        raise ValueError("benchmark must be a ConstantConductivityHalfSpace")
    exact = benchmark.cell_averages(depth_edges_m=profile.depth_edges_m, time_s=profile.time_s)
    widths = np.diff(profile.depth_edges_m)
    span = float(profile.depth_edges_m[-1] - profile.depth_edges_m[0])
    errors = {}
    for name, units in _FIELDS:
        observed, expected = getattr(profile, name), getattr(exact, name)
        ref_logs = [math.log(abs(float(value))) if value != 0 else -math.inf for value in expected]
        diff_logs = [_log_abs_difference(float(value), float(reference))
                     for value, reference in zip(observed, expected)]
        ref_norms = _log_norms(ref_logs, widths, span)
        diff_norms = _log_norms(diff_logs, widths, span)
        is_zero = all(value == -math.inf for value in ref_logs)
        errors[name] = {
            "units": units,
            "relative_units": "1",
            "absolute": {key: _finite_exp(value, f"{name} absolute {key}") for key, value in diff_norms.items()},
            "relative": {key: None if is_zero else _finite_exp(diff_norms[key] - value, f"{name} relative {key}")
                         for key, value in ref_norms.items()},
            "reference_norms": {key: _finite_exp(value, f"{name} reference {key}") for key, value in ref_norms.items()},
            "reference_is_zero": is_zero,
        }
    return {
        "schema": HEATING_COMPARISON_SCHEMA,
        "profile": profile.metadata,
        "reference": {
            "model": "constant_conductivity_half_space_boundary_step",
            "conductivity_s_m": benchmark.conductivity_s_m,
            "permeability_h_m": benchmark.permeability_h_m,
            "volumetric_heat_capacity_j_m3_k": benchmark.volumetric_heat_capacity_j_m3_k,
            "boundary_field_t": benchmark.boundary_field_t,
            "magnetic_diffusivity_m2_s": benchmark.magnetic_diffusivity_m2_s,
            "assumptions": list(_ASSUMPTIONS),
        },
        "grid": {"cells": int(widths.size), "depth_start_m": float(profile.depth_edges_m[0]),
                 "depth_end_m": float(profile.depth_edges_m[-1]), "depth_span_m": span,
                 "weighting": "cell_width_normalized_by_depth_span", "sampling": "cell_average"},
        "zero_reference_policy": "null_relative_norms_when_all_evaluated_reference_cells_are_zero",
        "errors": errors,
    }


def _paths(path, metadata_path) -> tuple[Path, Path]:
    csv_path = Path(path)
    sidecar_path = Path(str(csv_path) + ".json") if metadata_path is None else Path(metadata_path)
    if csv_path.resolve() == sidecar_path.resolve():
        raise ValueError("CSV and metadata paths must be distinct")
    return csv_path, sidecar_path


def write_heating_profile_csv(path, profile: HeatingProfile, *, metadata_path=None) -> tuple[Path, Path]:
    """Write canonical SI CSV and its required hash-linked JSON sidecar."""
    if not isinstance(profile, HeatingProfile):
        raise ValueError("profile must be a HeatingProfile")
    csv_path, sidecar_path = _paths(path, metadata_path)
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(CSV_COLUMNS)
    for index in range(profile.depth_edges_m.size - 1):
        writer.writerow(format(float(value), ".17g") for value in (
            profile.depth_edges_m[index], profile.depth_edges_m[index + 1],
            profile.magnetic_field_t[index], profile.current_density_a_m2[index],
            profile.joule_power_w_m3[index],
        ))
    payload = stream.getvalue().encode("utf-8")
    metadata = profile.metadata
    metadata["csv_sha256"] = hashlib.sha256(payload).hexdigest()
    csv_path.write_bytes(payload)
    sidecar_path.write_text(json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return csv_path, sidecar_path


def _json_pairs(pairs):
    record = {}
    for key, value in pairs:
        if key in record:
            raise ValueError(f"duplicate metadata key: {key}")
        record[key] = value
    return record


def _invalid_json_constant(value):
    raise ValueError(f"nonfinite JSON constant: {value}")


def load_heating_profile_csv(path, *, metadata_path=None) -> HeatingProfile:
    """Load fixed-SI cell averages; reject omitted, inconsistent or changed metadata."""
    csv_path, sidecar_path = _paths(path, metadata_path)
    try:
        payload = csv_path.read_bytes()
        metadata = json.loads(sidecar_path.read_text(encoding="utf-8"), object_pairs_hook=_json_pairs,
                              parse_constant=_invalid_json_constant)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read heating CSV and required metadata: {exc}") from exc
    expected_keys = {"schema", "time_s", "source_kind", "run_id", "source_reference", "flash_version",
                     "coordinate_system", "cell_sampling", "units", "csv_sha256"}
    if not isinstance(metadata, dict) or set(metadata) != expected_keys:
        raise ValueError("heating metadata must contain exactly the canonical schema fields")
    if metadata["schema"] != HEATING_PROFILE_SCHEMA or metadata["units"] != PROFILE_UNITS:
        raise ValueError("heating metadata schema or fixed SI units are inconsistent")
    if metadata["coordinate_system"] != COORDINATE_SYSTEM or metadata["cell_sampling"] != "cell_average":
        raise ValueError("heating metadata must describe the canonical cell-average coordinates")
    if metadata["csv_sha256"] != hashlib.sha256(payload).hexdigest():
        raise ValueError("CSV SHA256 does not match its metadata")
    try:
        rows = list(csv.reader(io.StringIO(payload.decode("utf-8"), newline=""), strict=True))
    except (UnicodeError, csv.Error) as exc:
        raise ValueError("CSV must be a valid UTF-8 table") from exc
    if len(rows) < 2 or tuple(rows[0]) != CSV_COLUMNS:
        raise ValueError("CSV requires canonical fixed-SI columns and at least one cell")
    values = []
    for row in rows[1:]:
        if len(row) != len(CSV_COLUMNS):
            raise ValueError("CSV must contain exactly five values per cell")
        try:
            values.append([float(item) for item in row])
        except (ValueError, OverflowError) as exc:
            raise ValueError("CSV cells must be finite real numeric values") from exc
    table = np.asarray(values, dtype=float)
    if np.any(~np.isfinite(table)) or np.any(table[:, 1] <= table[:, 0]):
        raise ValueError("CSV cells must be finite with positive depth widths")
    if np.any(table[1:, 0] != table[:-1, 1]):
        raise ValueError("CSV depth intervals must be exactly contiguous without gaps or overlaps")
    return HeatingProfile(
        time_s=metadata["time_s"], depth_edges_m=np.concatenate((table[:, 0], table[-1:, 1])),
        magnetic_field_t=table[:, 2], current_density_a_m2=table[:, 3], joule_power_w_m3=table[:, 4],
        source_kind=metadata["source_kind"], run_id=metadata["run_id"],
        source_reference=metadata["source_reference"], flash_version=metadata["flash_version"],
    )
