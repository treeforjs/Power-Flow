"""Set-valued comparison of complete neutral forward cases with curve readouts.

This module runs no kinetics or transport, rescales no population history and
interpolates no forward case. A curve envelope describes figure-readout support,
not independent noise samples. Its worst normalized departure is therefore a
compatibility score, never a likelihood, chi-square or confidence interval.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Real

import numpy as np


SCHEMA = "mhdlab-neutral-ensemble-comparison-v1"
RADIANCE_ROLES = frozenset({"nominal_energy_radiance", "nominal_photon_radiance"})
SCALE_ROLES = RADIANCE_ROLES | {"published_derived_unscaled"}
PHYSICAL_STATUSES = frozenset({"observed", "reconstructed", "owner_supplied", "source_evaluated"})
_MISSING_EVIDENCE = frozenset({
    "unknown", "unavailable", "pending", "missing", "none", "null", "n/a", "na",
    "tbd", "todo", "unspecified", "unregistered", "unverified", "blocked",
    "not available", "not provided", "not supplied", "not established", "not qualified",
    "to be determined", "to be supplied", "to be provided", "not applicable",
})


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be an explicit nonempty string")
    return value


def _evidence_text(value, name, status):
    """Reject missing-evidence sentinels, without authenticating free-text claims."""
    _text(value, name)
    normalized = " ".join(value.strip().casefold().replace("_", " ").replace("-", " ").split())
    normalized = normalized.strip(" .,:;!?<>[](){}")
    if status in PHYSICAL_STATUSES and normalized in _MISSING_EVIDENCE:
        raise ValueError(f"physical {name} requires concrete evidence; missing-evidence sentinel {value!r} is inadmissible")
    return value


def _number(value, name, *, positive=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not np.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{name} must be finite{' and positive' if positive else ''}")
    return result


def _array(value, name):
    if np.ma.isMaskedArray(value):
        raise ValueError(f"{name} cannot contain masked values")
    objects = np.asarray(value, dtype=object)
    if objects.ndim != 1 or objects.size == 0:
        raise ValueError(f"{name} must be a nonempty real vector")
    values = np.array([_number(item, name) for item in objects], dtype=float)
    return np.frombuffer(values.tobytes(), dtype=float)


def _segments(value, length):
    values = np.asarray(value)
    booleans = any(isinstance(item, (bool, np.bool_)) for item in np.asarray(value, dtype=object).flat)
    if values.shape != (length,) or values.dtype.kind not in "iu" or booleans or np.any(values < 0):
        raise ValueError("support_segment_ids must be nonnegative integers, one per retained coordinate")
    if np.any(values > np.iinfo(np.int64).max):
        raise ValueError("support segment IDs must fit in signed 64-bit integers")
    # A segment may occupy only one contiguous run. No missing intervals are filled.
    runs = values[np.r_[True, values[1:] != values[:-1]]]
    if len(np.unique(runs)) != len(runs):
        raise ValueError("a support segment ID cannot reappear after another segment")
    return np.frombuffer(np.asarray(values, dtype=np.int64).tobytes(), dtype=np.int64)


def _curve_fields(curve):
    coordinate = _array(curve.coordinate, "coordinate")
    if np.any(np.diff(coordinate) <= 0):
        raise ValueError("retained coordinates must be strictly increasing")
    object.__setattr__(curve, "coordinate", coordinate)
    object.__setattr__(curve, "support_segment_ids", _segments(curve.support_segment_ids, len(coordinate)))
    for name in ("curve_id", "quantity", "signal_units", "coordinate_units", "provenance"):
        _text(getattr(curve, name), name)
    return coordinate


@dataclass(frozen=True, eq=False)
class ReadoutEnvelope:
    curve_id: str
    coordinate: np.ndarray
    center: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    support_segment_ids: np.ndarray
    quantity: str
    signal_units: str
    coordinate_units: str
    scale_role: str
    provenance: str

    def __post_init__(self):
        coordinate = _curve_fields(self)
        for name in ("center", "lower", "upper"):
            values = _array(getattr(self, name), name)
            if values.shape != coordinate.shape:
                raise ValueError(f"{name} must match retained coordinates")
            object.__setattr__(self, name, values)
        if np.any(self.lower > self.center) or np.any(self.upper < self.center):
            raise ValueError("readout lower <= center <= upper is required")
        if self.scale_role not in SCALE_ROLES:
            raise ValueError(f"scale_role must be one of {sorted(SCALE_ROLES)}")
        if self.scale_role == "nominal_energy_radiance" and self.quantity not in {
            "spectral_energy_radiance", "integrated_energy_radiance",
        }:
            raise ValueError("nominal_energy_radiance role requires a declared energy-radiance quantity")
        if self.scale_role == "nominal_photon_radiance" and self.quantity not in {
            "spectral_photon_radiance", "integrated_photon_radiance",
        }:
            raise ValueError("nominal_photon_radiance role requires a declared photon-radiance quantity")


@dataclass(frozen=True, eq=False)
class ProjectedCurve:
    curve_id: str
    coordinate: np.ndarray
    values: np.ndarray
    support_segment_ids: np.ndarray
    quantity: str
    signal_units: str
    coordinate_units: str
    provenance: str

    def __post_init__(self):
        coordinate = _curve_fields(self)
        values = _array(self.values, "values")
        if values.shape != coordinate.shape:
            raise ValueError("prediction values must match retained coordinates")
        object.__setattr__(self, "values", values)


def _observation_fields(item):
    _text(item.reconstruction_id, "reconstruction_id")
    for name in ("time_registration", "observation_contract", "excitation_contract", "provenance"):
        _evidence_text(getattr(item, name), name, item.status)
    start = _number(item.gate_start_s, "gate_start_s")
    end = _number(item.gate_end_s, "gate_end_s")
    if end <= start:
        raise ValueError("registered gate must have positive duration")
    object.__setattr__(item, "gate_start_s", start)
    object.__setattr__(item, "gate_end_s", end)


@dataclass(frozen=True, eq=False)
class EnsembleObservation:
    reconstruction_id: str
    curves: tuple[ReadoutEnvelope, ...]
    status: str
    mapping_status: str
    time_registration: str
    gate_start_s: float
    gate_end_s: float
    observation_contract: str
    excitation_contract: str
    calibration_provenance: str
    provenance: str

    def __post_init__(self):
        _observation_fields(self)
        _text(self.calibration_provenance, "calibration_provenance")
        if self.status not in {"manufactured", "reconstructed", "observed"}:
            raise ValueError("observation status must be manufactured, reconstructed or observed")
        if self.mapping_status not in {"qualified", "scenario_assumption", "unavailable"}:
            raise ValueError("mapping_status must be qualified, scenario_assumption or unavailable")
        curves = tuple(self.curves)
        if not curves or any(not isinstance(c, ReadoutEnvelope) for c in curves):
            raise ValueError("observation requires explicit ReadoutEnvelope values")
        if len({c.curve_id for c in curves}) != len(curves):
            raise ValueError("observation curve IDs must be unique")
        for curve in curves:
            _evidence_text(curve.provenance, f"curve {curve.curve_id} provenance", self.status)
            if curve.scale_role in RADIANCE_ROLES:
                _evidence_text(self.calibration_provenance, "calibration_provenance", self.status)
        object.__setattr__(self, "curves", curves)


def _diagnostic_map(value, name):
    if not isinstance(value, dict) or not value:
        raise ValueError(f"{name} requires explicit keyed neutral species diagnostics")
    result = {}
    for key, number in value.items():
        _text(key, f"{name} species ID")
        result[key] = _number(number, name)
        if result[key] < 0:
            raise ValueError(f"{name} must be nonnegative")
    return result


@dataclass(frozen=True, eq=False)
class NeutralForwardCase:
    case_id: str
    hypothesis: str
    reconstruction_id: str
    parent_water_areal_inventory_m2: float
    thermal_parameter_K: float
    thermal_coordinate: str
    thermal_coordinate_definition: str
    curves: tuple[ProjectedCurve, ...]
    status: str
    time_registration: str
    gate_start_s: float
    gate_end_s: float
    observation_contract: str
    excitation_contract: str
    source_receipt: str
    population_provenance: str
    transport_provenance: str
    conservation_status: str
    conservation_provenance: str
    optical_admission: str
    numerical_qualification: str
    neutral_columns_m2: dict[str, float]
    weighted_neutral_densities_m3: dict[str, float]
    neutral_diagnostic_definition: str
    provenance: str

    def __post_init__(self):
        _observation_fields(self)
        for name in ("case_id", "hypothesis", "thermal_coordinate"):
            _text(getattr(self, name), name)
        for name in ("thermal_coordinate_definition", "source_receipt", "population_provenance", "transport_provenance",
                     "conservation_provenance", "optical_admission", "numerical_qualification",
                     "neutral_diagnostic_definition"):
            _evidence_text(getattr(self, name), name, self.status)
        if self.status not in {"manufactured", "owner_supplied", "source_evaluated"}:
            raise ValueError("case status must be manufactured, owner_supplied or source_evaluated")
        if self.conservation_status != "passed":
            raise ValueError("forward cases require passed number/energy accounting and its provenance")
        inventory = _number(self.parent_water_areal_inventory_m2, "parent_water_areal_inventory_m2")
        if inventory < 0:
            raise ValueError("parent water inventory must be nonnegative")
        object.__setattr__(self, "parent_water_areal_inventory_m2", inventory)
        object.__setattr__(self, "thermal_parameter_K", _number(self.thermal_parameter_K, "thermal_parameter_K", positive=True))
        curves = tuple(self.curves)
        if not curves or any(not isinstance(c, ProjectedCurve) for c in curves):
            raise ValueError("case requires explicit ProjectedCurve values")
        if len({c.curve_id for c in curves}) != len(curves):
            raise ValueError("prediction curve IDs must be unique")
        for curve in curves:
            _evidence_text(curve.provenance, f"curve {curve.curve_id} provenance", self.status)
        object.__setattr__(self, "curves", curves)
        for name in ("neutral_columns_m2", "weighted_neutral_densities_m3"):
            object.__setattr__(self, name, _diagnostic_map(getattr(self, name), name))
        if self.neutral_columns_m2.keys() != self.weighted_neutral_densities_m3.keys():
            raise ValueError("column and density diagnostic species IDs must agree")


def exploratory_calibration_scales():
    """The agreed 41-point sensitivity range, never an uncertainty distribution."""
    return np.power(10.0, np.linspace(-2.0, 2.0, 41))


def _score(values, envelope, scale):
    # Long-double intermediates avoid losing small authored readout widths when
    # applying the global scale. No unreported tolerance/noise floor is added.
    target = np.asarray(envelope.center, dtype=np.longdouble) * np.longdouble(scale)
    lower = np.asarray(envelope.lower, dtype=np.longdouble) * np.longdouble(scale)
    upper = np.asarray(envelope.upper, dtype=np.longdouble) * np.longdouble(scale)
    delta = np.asarray(values, dtype=np.longdouble) - target
    width = np.where(delta >= 0, upper-target, target-lower)
    if np.any((width == 0) & (delta != 0)):
        return float("inf")
    normalized = np.divide(np.abs(delta), width, out=np.zeros_like(delta), where=width != 0)
    score = np.max(normalized)
    if not np.isfinite(score) or score > np.finfo(float).max:
        return float("inf")
    return float(score)


def _case_metadata(case):
    return {
        "case_id": case.case_id, "hypothesis": case.hypothesis,
        "parent_water_areal_inventory_m2": case.parent_water_areal_inventory_m2,
        "thermal_parameter_K": case.thermal_parameter_K,
        "neutral_columns_m2": dict(case.neutral_columns_m2),
        "weighted_neutral_densities_m3": dict(case.weighted_neutral_densities_m3),
        "neutral_diagnostic_definition": case.neutral_diagnostic_definition,
        "status": case.status, "source_receipt": case.source_receipt,
        "population_provenance": case.population_provenance,
        "transport_provenance": case.transport_provenance,
        "conservation_provenance": case.conservation_provenance,
        "optical_admission": case.optical_admission,
        "numerical_qualification": case.numerical_qualification,
        "provenance": case.provenance,
    }


def compare_neutral_ensemble(observation, cases, calibration_scales=None):
    """Retain every supplied case/scale alternative, without interpolation.

    Each complete forward prediction is left untouched. ``C_abs`` multiplies
    nominal radiance readout centers and bounds only. Derived published curves
    cannot acquire an absolute interpretation through a multiplicative scan.
    """
    if not isinstance(observation, EnsembleObservation):
        raise ValueError("observation must be an EnsembleObservation")
    if observation.mapping_status == "unavailable":
        raise ValueError("published-to-forward mapping unavailable; withhold inference until a named reconstruction is admissible")
    cases = tuple(cases)
    if not cases or any(not isinstance(case, NeutralForwardCase) for case in cases):
        raise ValueError("a nonempty set of complete NeutralForwardCase values is required")
    if len({case.case_id for case in cases}) != len(cases):
        raise ValueError("forward case IDs must be unique")
    radiance = {c.curve_id for c in observation.curves if c.scale_role in RADIANCE_ROLES}
    scales = _array(exploratory_calibration_scales() if calibration_scales is None and radiance
                    else [1.0] if calibration_scales is None else calibration_scales, "calibration_scales")
    if np.any(scales <= 0) or np.any(np.diff(scales) <= 0):
        raise ValueError("calibration scales must be positive and strictly increasing")
    if not radiance and (len(scales) != 1 or scales[0] != 1):
        raise ValueError("a derived-display-only comparison admits only C_abs=1; it has no radiance scale mapping")
    observations = {curve.curve_id: curve for curve in observation.curves}
    first = cases[0]
    for case in cases:
        if (observation.status == "manufactured") != (case.status == "manufactured"):
            raise ValueError("manufactured and physical observations/cases cannot be mixed")
        for name in ("reconstruction_id", "time_registration", "gate_start_s", "gate_end_s",
                     "observation_contract", "excitation_contract"):
            if getattr(case, name) != getattr(observation, name):
                raise ValueError(f"case {case.case_id} has mismatched {name}")
        for name in ("thermal_coordinate", "thermal_coordinate_definition"):
            if getattr(case, name) != getattr(first, name):
                raise ValueError("all cases must use the same explicitly defined thermal coordinate")
        predicted = {curve.curve_id: curve for curve in case.curves}
        if predicted.keys() != observations.keys():
            raise ValueError("every case must cover exactly the retained observation curve IDs")
        for curve_id, envelope in observations.items():
            curve = predicted[curve_id]
            for name in ("quantity", "signal_units", "coordinate_units"):
                if getattr(curve, name) != getattr(envelope, name):
                    raise ValueError(f"case {case.case_id}, curve {curve_id}: mismatched {name}")
            for name in ("coordinate", "support_segment_ids"):
                if not np.array_equal(getattr(curve, name), getattr(envelope, name)):
                    raise ValueError(f"case {case.case_id}, curve {curve_id}: exact retained {name} required; no interpolation")

    results, groups = [], []
    hypotheses = sorted({case.hypothesis for case in cases})
    compatible_any = set()
    for scale in scales:
        scored = []
        for case in cases:
            curve_scores = {curve.curve_id: _score(curve.values, observations[curve.curve_id],
                                                  scale if curve.curve_id in radiance else 1.0)
                            for curve in case.curves}
            score = max(curve_scores.values())
            compatible = score <= 1
            if compatible:
                compatible_any.add(case.case_id)
            row = {"case_id": case.case_id, "hypothesis": case.hypothesis, "C_abs": float(scale),
                   "max_normalized_departure": score if np.isfinite(score) else None,
                   "unbounded_departure": not np.isfinite(score),
                   "readout_compatible": compatible,
                   "curve_scores": {key: value if np.isfinite(value) else None for key, value in curve_scores.items()}}
            results.append(row)
            scored.append((case, score))
        for hypothesis in hypotheses:
            selected = [(case, score) for case, score in scored if case.hypothesis == hypothesis]
            best = min(score for _, score in selected)
            compatible = [case for case, score in selected if score <= 1]
            groups.append({"hypothesis": hypothesis, "C_abs": float(scale),
                           "minimum_normalized_departure": best if np.isfinite(best) else None,
                           "all_departures_unbounded": not np.isfinite(best),
                           "minimum_case_ids": [case.case_id for case, score in selected if score == best],
                           "compatible_case_ids": [case.case_id for case in compatible],
                           "compatible_inventory_values_m2": sorted({case.parent_water_areal_inventory_m2 for case in compatible}),
                           "compatible_thermal_values_K": sorted({case.thermal_parameter_K for case in compatible})})
    compatible_cases = [case for case in cases if case.case_id in compatible_any]
    return {
        "schema": SCHEMA, "reconstruction_id": observation.reconstruction_id,
        "mapping_status": observation.mapping_status, "observation_status": observation.status,
        "score_definition": "maximum over all retained curves/points of absolute displacement from center divided by the lower or upper readout half-width on that side; <=1 means every point is inside its envelope",
        "zero_half_width_policy": "exact equality has zero departure; any displacement toward a zero half-width is unbounded and serialized as null with an explicit flag",
        "calibration": {"definition": "C_abs multiplies nominal-to-scenario-corrected radiance observations, never the forward population or detector responsivity",
                        "scales": scales.tolist(), "scaled_curve_ids": sorted(radiance),
                        "provenance": observation.calibration_provenance,
                        "interpretation": "exploratory sensitivity scenarios; not a confidence interval or probability distribution"},
        "registered_observation": {name: getattr(observation, name) for name in
                                   ("time_registration", "gate_start_s", "gate_end_s", "observation_contract", "excitation_contract", "provenance")},
        "thermal_coordinate": first.thermal_coordinate,
        "thermal_coordinate_definition": first.thermal_coordinate_definition,
        "case_metadata": [_case_metadata(case) for case in cases],
        "case_scale_comparisons": results, "hypothesis_scale_summaries": groups,
        "diagnostics": {
            "compatible_at_any_scanned_scale_case_ids": [case.case_id for case in compatible_cases],
            "compatible_inventory_values_m2": sorted({case.parent_water_areal_inventory_m2 for case in compatible_cases}),
            "compatible_thermal_values_K": sorted({case.thermal_parameter_K for case in compatible_cases}),
            "multiple_compatible_inventory_values": len({case.parent_water_areal_inventory_m2 for case in compatible_cases}) > 1,
            "multiple_compatible_thermal_values": len({case.thermal_parameter_K for case in compatible_cases}) > 1,
            "no_compatible_supplied_case": not compatible_cases,
            "forward_case_interpolation": False, "extrapolation": False,
            "populations_rescaled": False, "confidence_intervals_calculated": False,
        },
        "limits": [
            "Readout envelopes are reconstruction support, not stochastic errors; point counts do not increase statistical weight.",
            "All minima and compatible values refer only to supplied complete cases and scanned scales; unsampled alternatives remain unresolved.",
            "Matching contract labels and owner attestations are admission requirements, not independent validation of kinetics, transport, conservation, opacity or convergence.",
            "Source inventory is the single fitted amplitude coordinate; no independent density or monolayer amplitude is fitted.",
            "Neutral columns/densities and the declared thermal coordinate are supplied forward diagnostics; internal populations do not automatically define translational gas temperature.",
            "A categorical reconstruction assumption or exploratory calibration scale does not establish a verified absolute calibration or physical estimate.",
        ],
    }
