"""Explicit, bounded Mills (2004) austenitic-solid mass-specific Cp and H.

The reviewed evaluated equations are not a foil-specific EOS.  A caller must
select a solid interval, material/phase qualifications and enthalpy convention.
No density, Cv, internal energy, fusion, liquid branch or transport splice is
inferred.  This module contains no installed solver or material default.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from numbers import Integral, Real
from pathlib import Path
from types import MappingProxyType
from typing import Any, ClassVar

import numpy as np


MILLS_MANIFEST_SHA256 = "ff64f70a2c605cf1024cc8ecfff9ac113d325981feb4efb9cc69dc5340205a8a"
MILLS_PDF_SHA256 = "bf120dcbdf20d5f9b00494efd5e05643793edc0838deee8ba2c52988f9074836"
MILLS_MATERIAL_FAMILY = "austenitic stainless steel"
MILLS_PHASE = "solid"
QUALIFIED_SOURCE_TEMPERATURE_RANGE_K = (300.0, 1700.0)
REFERENCE_CONVENTIONS = ("printed_eq11", "rezeroed_eq11")
MODEL_SCHEMA = "mhdlab-mills2004-bounded-solid-enthalpy-v1"


def _real(value, name: str) -> float:
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number, not bool, masked, string or complex")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


def _values(value, name: str) -> np.ndarray:
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked values")
    try:
        raw = np.asarray(value, dtype=object)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{name} must be a nonempty finite real scalar or array") from exc
    if raw.size == 0:
        raise ValueError(f"{name} cannot be empty")
    return np.array([_real(item, name) for item in raw.flat], dtype=float).reshape(raw.shape)


def _output(value):
    array = np.asarray(value, dtype=float)
    if array.ndim == 0:
        return float(array.item())
    # A bytes owner prevents callers from re-enabling the NumPy writeable flag.
    return np.frombuffer(array.tobytes(order="C"), dtype=float).reshape(array.shape)


def _label(value, name: str) -> str:
    if (not isinstance(value, str) or not value.strip() or value != value.strip()
            or any(ord(character) < 32 for character in value)):
        raise ValueError(f"{name} must be a nonempty string without surrounding whitespace or control characters")
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


def _fraction(temperature):
    return Fraction.from_float(float(temperature))


def _map_source(function, temperature):
    query = np.asarray(temperature, dtype=float)
    values = np.array([float(function(_fraction(t))) for t in query.flat]).reshape(query.shape)
    return float(values.item()) if values.ndim == 0 else values


def _cp_fraction(temperature):
    return 472 + Fraction(17, 125) * temperature - 2820000 / temperature**2


def _printed_h_fraction(temperature):
    return 472 * temperature + Fraction(17, 250) * temperature**2 + 2820000 / temperature - 156000


def _delta_h_fraction(temperature, reference):
    # Integral of Eq. 10, factored before evaluation to retain small differences.
    return ((temperature - reference)
            * (472 + Fraction(17, 250) * (temperature + reference)
               - 2820000 / (temperature * reference)))


def _cp(temperature):
    return _map_source(_cp_fraction, temperature)


def _printed_h(temperature):
    return _map_source(_printed_h_fraction, temperature)


def _delta_h(temperature, reference):
    ref = _fraction(reference)
    return _map_source(lambda t: _delta_h_fraction(t, ref), temperature)


@dataclass(frozen=True, slots=True, eq=False)
class SolidEnthalpyValues:
    """Immutable scalar or same-shaped array evaluation, all mass-specific SI."""

    temperature_k: Any
    cp_j_kg_k: Any
    printed_enthalpy_j_kg: Any
    delta_enthalpy_from_reference_j_kg: Any
    enthalpy_j_kg: Any

    def __post_init__(self):
        names = ("temperature_k", "cp_j_kg_k", "printed_enthalpy_j_kg",
                 "delta_enthalpy_from_reference_j_kg", "enthalpy_j_kg")
        arrays = {name: _values(getattr(self, name), name) for name in names}
        shape = arrays["temperature_k"].shape
        if any(value.shape != shape for value in arrays.values()):
            raise ValueError("enthalpy result arrays must have the same scalar or array shape")
        if np.any(arrays["temperature_k"] <= 0) or np.any(arrays["cp_j_kg_k"] <= 0):
            raise ValueError("temperature and isobaric heat capacity must be positive")
        for name, value in arrays.items():
            object.__setattr__(self, name, _output(value))


@dataclass(frozen=True, slots=True, eq=False, init=False)
class Mills2004AusteniticSolidEnthalpy:
    """Reviewed evaluated solid branch with explicit support and reference.

    Construct with :func:`load_mills_2004_solid_enthalpy`.  Source qualification
    pins the assessment bytes and source-PDF metadata; no external PDF or network
    access is needed to evaluate the reviewed equations.
    """

    temperature_min_k: float
    temperature_max_k: float
    reference_temperature_k: float
    reference_convention: str
    material_family: str
    phase: str
    material_qualification: str
    phase_qualification: str
    source_manifest_sha256: str
    source_pdf_sha256: str
    source_metadata: Mapping[str, Any]
    parameter_status: ClassVar[str] = "evaluated_source_correlation"

    def __init__(self, *args, **kwargs):
        raise TypeError("use load_mills_2004_solid_enthalpy with explicit source and qualifications")

    @property
    def supported_temperature_range_k(self) -> tuple[float, float]:
        return self.temperature_min_k, self.temperature_max_k

    @property
    def printed_enthalpy_at_reference_j_kg(self) -> float:
        return float(_printed_h(self.reference_temperature_k))

    @property
    def supported_enthalpy_range_j_kg(self) -> tuple[float, float]:
        return (self.enthalpy_j_kg(self.temperature_min_k),
                self.enthalpy_j_kg(self.temperature_max_k))

    @property
    def metadata(self) -> dict:
        """Detached JSON-compatible metadata; modifying it cannot alter the model."""
        return {
            "schema": MODEL_SCHEMA,
            "material_family": self.material_family,
            "phase": self.phase,
            "material_qualification": self.material_qualification,
            "phase_qualification": self.phase_qualification,
            "parameter_status": self.parameter_status,
            "source_manifest_sha256": self.source_manifest_sha256,
            "source_pdf_sha256": self.source_pdf_sha256,
            "source_reference": "Mills et al. (2004), DOI 10.2355/isijinternational.44.1661, p.1663, Eqs.10-11",
            "supported_temperature_range_k": list(self.supported_temperature_range_k),
            "qualified_source_temperature_range_k": list(QUALIFIED_SOURCE_TEMPERATURE_RANGE_K),
            "supported_enthalpy_range_j_kg": list(self.supported_enthalpy_range_j_kg),
            "reference_temperature_k": self.reference_temperature_k,
            "reference_convention": self.reference_convention,
            "printed_enthalpy_at_reference_j_kg": self.printed_enthalpy_at_reference_j_kg,
            "declared_enthalpy_at_reference_j_kg": self.enthalpy_j_kg(self.reference_temperature_k),
            "printed_eq11_nominal_source_reference_k": 298.0,
            "printed_eq11_reference_rounding": _plain(self.source_metadata["reference_sanity"]),
            "units": {"temperature_k": "K", "cp_j_kg_k": "J kg^-1 K^-1",
                      "printed_enthalpy_j_kg": "J kg^-1",
                      "delta_enthalpy_from_reference_j_kg": "J kg^-1", "enthalpy_j_kg": "J kg^-1"},
            "cp_is_cv": False,
            "enthalpy_is_internal_energy": False,
            "actual_dissertation_foil_validated": False,
            "density_or_fusion_supplied": False,
            "extrapolation_or_clamping": False,
            "adopted_as_solver_default": False,
            "source_metadata": _plain(self.source_metadata),
        }

    def _temperatures(self, temperature_k):
        query = _values(temperature_k, "temperature_k")
        if np.any(query < self.temperature_min_k) or np.any(query > self.temperature_max_k):
            raise ValueError("temperature outside the explicitly supported solid interval; no clamping or extrapolation")
        return query

    def cp_j_kg_k(self, temperature_k):
        """Eq. 10 mass-specific isobaric heat capacity, not Cv."""
        return _output(_cp(self._temperatures(temperature_k)))

    def printed_enthalpy_j_kg(self, temperature_k):
        """Eq. 11 exactly as printed, retaining its rounded integration constant."""
        return _output(_printed_h(self._temperatures(temperature_k)))

    def delta_enthalpy_from_reference_j_kg(self, temperature_k):
        """Stable integral of Cp from the explicitly declared in-domain Tref."""
        return _output(_delta_h(self._temperatures(temperature_k), self.reference_temperature_k))

    def enthalpy_j_kg(self, temperature_k):
        """Enthalpy in the selected printed or explicitly re-zeroed convention."""
        query = self._temperatures(temperature_k)
        if self.reference_convention == "printed_eq11":
            return _output(_printed_h(query))
        return _output(_delta_h(query, self.reference_temperature_k))

    def evaluate(self, temperature_k) -> SolidEnthalpyValues:
        query = self._temperatures(temperature_k)
        printed = _printed_h(query)
        delta = _delta_h(query, self.reference_temperature_k)
        declared = printed if self.reference_convention == "printed_eq11" else delta
        return SolidEnthalpyValues(query, _cp(query), printed, delta, declared)

    def _inverse_one(self, enthalpy: float) -> float:
        lower, upper = self.supported_temperature_range_k
        bounds = self.supported_enthalpy_range_j_kg
        if enthalpy == bounds[0]:
            return lower
        if enthalpy == bounds[1]:
            return upper
        reference = self.reference_temperature_k
        reference_exact = _fraction(reference)
        target = _fraction(enthalpy)
        if self.reference_convention == "printed_eq11":
            target -= _printed_h_fraction(reference_exact)
        trial = reference + float(target) / _cp(reference)
        if not lower < trial < upper:
            trial = lower + (upper - lower) / 2.0
        # Safeguarded Newton iterations use the stable integral as their residual.
        # The final adjacent-float comparison requires no arbitrary Kelvin tolerance.
        for _ in range(256):
            residual = _delta_h_fraction(_fraction(trial), reference_exact) - target
            if residual == 0:
                return float(trial)
            if residual < 0:
                lower = trial
            else:
                upper = trial
            if np.nextafter(lower, upper) >= upper:
                break
            candidate = trial - float(residual) / _cp(trial)
            if not lower < candidate < upper:
                candidate = lower + (upper - lower) / 2.0
            if candidate == trial:
                break
            trial = candidate
        candidates = (lower, upper, trial, float(np.nextafter(trial, lower)),
                      float(np.nextafter(trial, upper)))
        return float(min(candidates, key=lambda t:
                         abs(_delta_h_fraction(_fraction(t), reference_exact) - target)))

    def temperature_from_enthalpy_k(self, enthalpy_j_kg):
        """Invert the selected monotone branch inside its exact declared bounds.

        Endpoints are included. Every input outside the declared enthalpy interval
        raises; values are never clipped. Float64 H resolution limits temperature
        resolution in the printed-offset convention. A re-zeroed near-reference
        query retains the small integral before its conversion back to temperature.
        """
        query = _values(enthalpy_j_kg, "enthalpy_j_kg")
        lower, upper = self.supported_enthalpy_range_j_kg
        if np.any(query < lower) or np.any(query > upper):
            raise ValueError("enthalpy outside the explicitly supported solid interval; no clamping or extrapolation")
        result = np.array([self._inverse_one(float(value)) for value in query.flat]).reshape(query.shape)
        return _output(result)


def _json_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate source manifest key: {key}")
        result[key] = value
    return result


def _invalid_constant(value):
    raise ValueError(f"nonfinite source manifest constant: {value}")


def _reviewed_manifest(manifest_path):
    try:
        payload = Path(manifest_path).read_bytes()
    except (OSError, TypeError) as exc:
        raise ValueError(f"cannot read reviewed Mills source manifest: {exc}") from exc
    if hashlib.sha256(payload).hexdigest() != MILLS_MANIFEST_SHA256:
        raise ValueError("Mills source manifest SHA256 does not match the reviewed source assessment")
    try:
        metadata = json.loads(payload.decode("utf-8"), object_pairs_hook=_json_pairs,
                              parse_constant=_invalid_constant)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError("Mills source manifest must be valid UTF-8 JSON") from exc
    expected = {
        "schema": "mhdlab-mills2004-source-assessment-v1",
        "status": "official_full_text_inspected_evaluated_correlations",
        "source.doi": "10.2355/isijinternational.44.1661",
        "source.pdf_sha256": MILLS_PDF_SHA256,
        "source.pdf_pages": 8,
        "source.new_characterized_dissertation_foil_measurement": False,
        "material.actual_dissertation_foil_composition": None,
        "material.actual_dissertation_foil_condition": None,
        "material.actual_dissertation_foil_phase_fractions": None,
        "material.solidus_or_mushy_fraction_law_supplied": False,
        "reference_sanity.d_enthalpy_dT_matches_cp": True,
        "reference_sanity.reference_correction_installed": False,
    }
    for path, required in expected.items():
        actual = metadata
        for key in path.split("."):
            if not isinstance(actual, dict) or key not in actual:
                raise ValueError(f"missing reviewed Mills source metadata: {path}")
            actual = actual[key]
        if type(actual) is not type(required) or actual != required:
            raise ValueError(f"inconsistent reviewed Mills source identity or phase qualification: {path}")
    if "3-series metastable austenitic" not in metadata["material"]["families"]:
        raise ValueError("reviewed Mills metadata does not identify the austenitic material family")
    properties = {item["name"]: item for item in metadata["properties"]}
    for name, equation, formula in (
        ("mass_specific_isobaric_heat_capacity_solid", 10, "cp(T)=472+0.136*T-2.82e6/T^2"),
        ("mass_specific_enthalpy_relative_to_H298_solid", 11,
         "H(T)-H298=472*T+0.068*T^2+2.82e6/T-156000"),
    ):
        item = properties.get(name, {})
        if (item.get("equation") != equation or item.get("formula") != formula
                or item.get("evidence") != {"pdf_page_1based": 3, "printed_page": 1663}):
            raise ValueError("reviewed Mills solid equations or source pages do not match")
    return _freeze(metadata)


def load_mills_2004_solid_enthalpy(
    manifest_path, *, temperature_min_k, temperature_max_k, reference_temperature_k,
    reference_convention: str, expected_material_family: str, expected_phase: str,
    material_qualification: str, phase_qualification: str,
) -> Mills2004AusteniticSolidEnthalpy:
    """Select the reviewed evaluated austenitic solid with no physical defaults.

    All selections are required. The caller's interval must lie inside the
    conservative 300--1700 K qualification and its Tref must lie in that interval.
    Qualifying a generic evaluated branch does not validate the actual foil.
    A changed manifest or source needs review and a separately versioned adapter.
    """
    if not isinstance(expected_material_family, str) or expected_material_family != MILLS_MATERIAL_FAMILY:
        raise ValueError("expected_material_family must explicitly be 'austenitic stainless steel'")
    if not isinstance(expected_phase, str) or expected_phase != MILLS_PHASE:
        raise ValueError("expected_phase must explicitly be 'solid'; no liquid or mixed branch is supplied")
    if not isinstance(reference_convention, str) or reference_convention not in REFERENCE_CONVENTIONS:
        raise ValueError("reference_convention must explicitly be 'printed_eq11' or 'rezeroed_eq11'")
    material_qualification = _label(material_qualification, "material_qualification")
    phase_qualification = _label(phase_qualification, "phase_qualification")
    lower = _real(temperature_min_k, "temperature_min_k")
    upper = _real(temperature_max_k, "temperature_max_k")
    reference = _real(reference_temperature_k, "reference_temperature_k")
    source_lower, source_upper = QUALIFIED_SOURCE_TEMPERATURE_RANGE_K
    if not source_lower <= lower < upper <= source_upper:
        raise ValueError("supported solid interval must satisfy 300 <= Tmin < Tmax <= 1700 K")
    if not lower <= reference <= upper:
        raise ValueError("declared reference_temperature_k must lie inside the supported solid interval")
    metadata = _reviewed_manifest(manifest_path)
    model = object.__new__(Mills2004AusteniticSolidEnthalpy)
    fields = {
        "temperature_min_k": lower, "temperature_max_k": upper,
        "reference_temperature_k": reference, "reference_convention": reference_convention,
        "material_family": expected_material_family, "phase": expected_phase,
        "material_qualification": material_qualification, "phase_qualification": phase_qualification,
        "source_manifest_sha256": MILLS_MANIFEST_SHA256, "source_pdf_sha256": MILLS_PDF_SHA256,
        "source_metadata": metadata,
    }
    for name, value in fields.items():
        object.__setattr__(model, name, value)
    # Cp'=.136+5.64e6/T^3 > 0 here; positive lower-end Cp proves monotone H.
    if _cp(lower) <= 0:
        raise ValueError("selected source interval does not have positive isobaric heat capacity")
    return model
