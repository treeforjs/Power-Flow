"""Explicit, bounded solid-material transport inputs in SI units.

Electrical resistivity and thermal conductivity are different observables.
Their piecewise-linear interpolation in temperature is a declared modeling
choice.  Electrical conductivity is derived as the reciprocal of interpolated
resistivity; it is not independently interpolated.  No value outside a supplied
table is inferred, and this module installs no material or solver defaults.
"""

from __future__ import annotations

from collections.abc import Mapping
import csv
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import hashlib
import io
import json
import math
from numbers import Integral, Real
from pathlib import Path
import re
from types import MappingProxyType
from typing import Any, ClassVar

import numpy as np


REPRESENTATION = "piecewise_linear_in_temperature"
TABLE_SCHEMA = "mhdlab-bounded-solid-transport-v1"
FLYNN_MATERIAL_GRADE = "AISI Type 304"
FLYNN_CSV_SHA256 = "ed63e23758d874f2d149ab04c76d0044f6c3cd9f5203d0da29996bc9fd91dbfb"
FLYNN_PDF_SHA256 = "ff6afb078d512aba2a07b40228a3997857f8bb5448eba465dca9c3d06c74edfa"
_FLYNN_COLUMNS = (
    "source_row", "temperature_source_degC", "temperature_K", "resistivity_source_mantissa",
    "resistivity_ohm_m", "thermal_conductivity_W_m_K",
)
_UNITS = {"temperature_k": "K", "resistivity_ohm_m": "ohm m",
          "thermal_conductivity_w_m_k": "W m^-1 K^-1", "electrical_conductivity_s_m": "S m^-1"}


def _real(value, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite positive real value, not bool, string or complex")
    try:
        number = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite and positive") from exc
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return number


def _values(value, name: str, *, scalar: bool = False) -> np.ndarray:
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked values")
    try:
        raw = np.asarray(value, dtype=object)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{name} must be a finite real array") from exc
    if raw.size == 0 or (raw.ndim not in (0, 1) if scalar else raw.ndim != 1):
        raise ValueError(f"{name} must be " + ("a scalar or nonempty 1-D array" if scalar else "a nonempty 1-D array"))
    return np.array([_real(item, name) for item in raw.flat], dtype=float).reshape(raw.shape)


def _immutable(value):
    array = np.ascontiguousarray(value, dtype=float)
    return np.frombuffer(array.tobytes(), dtype=float).reshape(np.asarray(value).shape)


def _output(value: np.ndarray):
    return float(value) if value.ndim == 0 else _immutable(value)


def _label(value, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"{name} must be a nonempty string without surrounding whitespace")
    if any(ord(character) < 32 for character in value):
        raise ValueError(f"{name} cannot contain control characters")
    return value


def _freeze(value):
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise ValueError("source metadata keys must be strings")
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Real) and math.isfinite(float(value)):
        return float(value)
    raise ValueError("source metadata must contain finite JSON-compatible values")


def _plain(value):
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_plain(item) for item in value]
    return value


@dataclass(frozen=True, slots=True, eq=False)
class TransportValues:
    """Immutable interpolated values, with electrical sigma derived from rho."""

    temperature_k: Any
    resistivity_ohm_m: Any
    thermal_conductivity_w_m_k: Any
    electrical_conductivity_s_m: Any = field(init=False)
    representation: ClassVar[str] = REPRESENTATION

    def __post_init__(self):
        arrays = {name: _values(getattr(self, name), name, scalar=True) for name in
                  ("temperature_k", "resistivity_ohm_m", "thermal_conductivity_w_m_k")}
        shape = arrays["temperature_k"].shape
        if any(value.shape != shape for value in arrays.values()):
            raise ValueError("transport result arrays must share a scalar or 1-D shape")
        with np.errstate(over="ignore", divide="ignore"):
            sigma = 1. / arrays["resistivity_ohm_m"]
        if np.any(~np.isfinite(sigma)) or np.any(sigma <= 0):
            raise ValueError("derived electrical conductivity exceeds finite positive float64 range")
        for name, value in arrays.items():
            object.__setattr__(self, name, _output(value))
        object.__setattr__(self, "electrical_conductivity_s_m", _output(sigma))


@dataclass(frozen=True, slots=True, eq=False)
class BoundedTransportTable:
    """One supplied solid phase and specimen context, without extrapolation.

    All physical inputs and provenance labels are required.  The status is
    either an explicitly manufactured case or the supported recovered published
    smoothed table.  A future source needs its own reviewed provenance adapter.
    Unknown sample conditioning must be stated as unknown, not filled by a
    grade-based guess.  No density, heat capacity, EOS or phase transition is
    supplied or inferred by this record.
    """

    temperature_k: np.ndarray
    resistivity_ohm_m: np.ndarray
    thermal_conductivity_w_m_k: np.ndarray
    material_grade: str
    sample_condition: str
    phase: str
    parameter_status: str
    source_reference: str
    source_data_sha256: str
    source_metadata: Mapping[str, Any] = field(default_factory=dict)
    representation: ClassVar[str] = REPRESENTATION

    def __post_init__(self):
        arrays = {name: _values(getattr(self, name), name) for name in
                  ("temperature_k", "resistivity_ohm_m", "thermal_conductivity_w_m_k")}
        temperatures = arrays["temperature_k"]
        if temperatures.size < 2 or np.any(temperatures[1:] <= temperatures[:-1]):
            raise ValueError("temperature_k must have at least two strictly increasing positive values")
        if any(value.shape != temperatures.shape for value in arrays.values()):
            raise ValueError("temperature, resistivity and thermal conductivity shapes must match")
        for name, value in arrays.items():
            object.__setattr__(self, name, _immutable(value))
        for name in ("material_grade", "sample_condition", "phase", "parameter_status", "source_reference"):
            _label(getattr(self, name), name)
        if self.phase != "solid":
            raise ValueError("only a single explicitly solid phase is supported")
        if self.parameter_status not in ("manufactured", "recovered_published_smoothed_table"):
            raise ValueError("parameter_status must explicitly label manufactured or recovered published smoothed data")
        if not isinstance(self.source_data_sha256, str) or re.fullmatch("[0-9a-f]{64}", self.source_data_sha256) is None:
            raise ValueError("source_data_sha256 must be a lowercase SHA256 hexadecimal digest")
        if not isinstance(self.source_metadata, Mapping):
            raise ValueError("source_metadata must be a JSON-compatible mapping")
        object.__setattr__(self, "source_metadata", _freeze(self.source_metadata))

    @property
    def supported_temperature_range_k(self) -> tuple[float, float]:
        return float(self.temperature_k[0]), float(self.temperature_k[-1])

    @property
    def metadata(self) -> dict:
        return {
            "schema": TABLE_SCHEMA, "material_grade": self.material_grade,
            "sample_condition": self.sample_condition, "phase": self.phase,
            "parameter_status": self.parameter_status, "source_reference": self.source_reference,
            "source_data_sha256": self.source_data_sha256,
            "supported_temperature_range_k": list(self.supported_temperature_range_k),
            "representation": self.representation, "units": dict(_UNITS),
            "electrical_conductivity_representation": "reciprocal_of_interpolated_resistivity",
            "source_metadata": _plain(self.source_metadata),
        }

    def evaluate(self, temperature_k) -> TransportValues:
        """Interpolate rho/kappa inside the table and derive sigma=1/rho.

        Domain endpoints are supported; every out-of-range query raises.
        Ratios and products use exact values of the supplied float inputs before
        rounding each final property to float64. This retains a dimensional
        contribution even when its interpolation weight alone would underflow.
        No interpolation across phase boundaries occurs.
        """
        query = _values(temperature_k, "temperature_k", scalar=True)
        if np.any(query < self.temperature_k[0]) or np.any(query > self.temperature_k[-1]):
            raise ValueError("temperature lies outside the supported table interval; extrapolation is not permitted")
        flat = query.reshape(-1)
        index = np.minimum(np.searchsorted(self.temperature_k, flat, side="right") - 1,
                           self.temperature_k.size - 2)
        result = [np.empty_like(flat), np.empty_like(flat)]
        for location, (temperature, bin_index) in enumerate(zip(flat, index)):
            lower = Fraction.from_float(float(self.temperature_k[bin_index]))
            upper = Fraction.from_float(float(self.temperature_k[bin_index + 1]))
            weight = (Fraction.from_float(float(temperature)) - lower) / (upper - lower)
            for output, values in zip(result, (self.resistivity_ohm_m, self.thermal_conductivity_w_m_k)):
                first = Fraction.from_float(float(values[bin_index]))
                last = Fraction.from_float(float(values[bin_index + 1]))
                output[location] = float((1 - weight) * first + weight * last)
        result = [values.reshape(query.shape) for values in result]
        return TransportValues(query, *result)


def _json_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate metadata key: {key}")
        result[key] = value
    return result


def _invalid_constant(value):
    raise ValueError(f"nonfinite metadata constant: {value}")


def _at(metadata, path):
    value = metadata
    for key in path.split("."):
        if not isinstance(value, dict) or key not in value:
            raise ValueError(f"missing Flynn source metadata field: {path}")
        value = value[key]
    return value


def load_flynn_304_transport(csv_path, metadata_path, *, expected_material_grade: str) -> BoundedTransportTable:
    """Load the audited Flynn table, requiring explicit intended grade matching.

    This named adapter pins the recovered CSV and source-PDF digests.  Recovered
    rows are also checked in decimal arithmetic before conversion to float64.
    A different table/recovery needs review and a new adapter/version; changing
    CSV values and recomputing its sidecar hash cannot inherit this citation.
    """
    if expected_material_grade != FLYNN_MATERIAL_GRADE or not isinstance(expected_material_grade, str):
        raise ValueError("expected_material_grade must explicitly be 'AISI Type 304'; no grade substitution is supported")
    try:
        payload = Path(csv_path).read_bytes()
        metadata = json.loads(Path(metadata_path).read_text(encoding="utf-8"), object_pairs_hook=_json_pairs,
                              parse_constant=_invalid_constant)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read Flynn CSV and source metadata: {exc}") from exc
    expected = {
        "schema": "mhdlab-published-smoothed-transport-table-v1",
        "status": "recovered_published_smoothed_table",
        "source.author": "D. R. Flynn", "source.report_date": "1968-06-25",
        "source.report": "NBS Report 9869",
        "source.title": "Thermal conductivity and electrical resistivity of a sample of AISI type 304 stainless steel",
        "source.doi": "10.6028/NBS.RPT.9869",
        "source.url": "https://nvlpubs.nist.gov/nistpubs/Legacy/RPT/nbsreport9869.pdf",
        "source.pdf_sha256": FLYNN_PDF_SHA256,
        "source.table_text_sha256": "249300bc3e3ccd2509b2f3c55983f080c3a061aa0bf0f03013952e1e11544167",
        "source.printed_page": 12, "source.pdf_page_1based": 16,
        "recovery.rows": 10, "recovery.missing_rows_invented": False,
        "recovery.graph_digitization": False, "recovery.interpolation": False,
        "recovery.extrapolation": False, "recovery.fit": False,
        "observable.material": "AISI Type 304, one supplied solid plate sample",
        "observable.value_kind": "Published smoothed results from measurements, not original acquisition samples",
        "observable.sample_composition": None, "observable.sample_thermal_history": None,
        "observable.sample_context_evidence.printed_page": 1,
        "observable.sample_context_evidence.pdf_page_1based": 5,
        "observable.actual_dissertation_foil_match": "Unknown; source sample is not the dissertation foil",
        "observable.temperature_source_unit": "degC", "observable.temperature_si_unit": "K",
        "observable.supported_temperature_range_degC": [-20, 150],
        "observable.supported_temperature_range_K": [253.15, 423.15],
        "observable.resistivity_unit": "ohm m", "observable.conductivity_unit": "W m^-1 K^-1",
        "conversion.resistivity_scale_ohm_m": "1e-7",
        "source_error_statements.coverage_factor": None,
        "source_error_statements.confidence_level": None,
        "source_error_statements.standard_deviation_interpretation": False,
        "source_error_statements.between_heat_or_foil_variability_quantified": False,
        "source_error_statements.printed_rounding_is_separate_from_source_error": True,
        "source_error_statements.evidence.printed_page": 12,
        "source_error_statements.evidence.pdf_page_1based": 16,
        "source_error_statements.resistivity.plus_minus_percent": .5,
        "source_error_statements.resistivity.interpretation": "Author believes error no greater than this value",
        "source_error_statements.conductivity.percent_range": [2, 3],
        "source_error_statements.conductivity.interpretation": "Author estimates uncertainty no greater than two to three percent; chiefly uncorrected insulation heat loss",
        "output.csv_filename": "flynn_1968_304_transport.csv",
        "output.csv_sha256": FLYNN_CSV_SHA256, "output.rows": 10,
        "output.columns": list(_FLYNN_COLUMNS), "output.adopted_as_runtime_defaults": False,
    }
    for path, required in expected.items():
        actual = _at(metadata, path)
        if type(actual) is not type(required) or actual != required:
            raise ValueError(f"Flynn source metadata inconsistent with the reviewed table: {path}")
    digest = hashlib.sha256(payload).hexdigest()
    if digest != FLYNN_CSV_SHA256:
        raise ValueError("Flynn CSV SHA256 does not match the reviewed recovered table")
    try:
        rows = list(csv.reader(io.StringIO(payload.decode("utf-8"), newline=""), strict=True))
    except (UnicodeError, csv.Error) as exc:
        raise ValueError("Flynn CSV must be a valid UTF-8 table") from exc
    if len(rows) != 11 or tuple(rows[0]) != _FLYNN_COLUMNS:
        raise ValueError("Flynn CSV must have the reviewed six-column header and ten source rows")
    temperatures, resistivities, conductivities = [], [], []
    for number, row in enumerate(rows[1:], start=1):
        if len(row) != 6:
            raise ValueError("Flynn CSV source rows must have six columns")
        try:
            source_row, celsius, kelvin, mantissa, rho, kappa = map(Decimal, row)
        except InvalidOperation as exc:
            raise ValueError("Flynn CSV source cells must be finite decimal numbers") from exc
        if not all(value.is_finite() for value in (source_row, celsius, kelvin, mantissa, rho, kappa)):
            raise ValueError("Flynn CSV source cells must be finite")
        if source_row != number or celsius != celsius.to_integral_value():
            raise ValueError("Flynn source row numbers and integer Celsius temperatures are inconsistent")
        if kelvin != celsius + Decimal("273.15") or rho != mantissa * Decimal("1e-7"):
            raise ValueError("Flynn Kelvin or resistivity scale conversion is inconsistent")
        temperatures.append(float(kelvin))
        resistivities.append(float(rho))
        conductivities.append(float(kappa))
    if (temperatures[0], temperatures[-1], resistivities[0], resistivities[-1],
            conductivities[0], conductivities[-1]) != (253.15, 423.15, 6.82e-7, 8.19e-7, 13.9, 17.9):
        raise ValueError("Flynn recovered temperature and transport anchors are inconsistent")
    return BoundedTransportTable(
        temperatures, resistivities, conductivities, FLYNN_MATERIAL_GRADE,
        "Supplied solid plate sample; composition and thermal history unknown",
        "solid", metadata["status"],
        "D. R. Flynn, NBS Report 9869 (1968), DOI 10.6028/NBS.RPT.9869, printed p.12 (PDF p.16)",
        digest, metadata,
    )
