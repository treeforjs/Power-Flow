"""Manufactured matched-input width assembly of supplied coupled histories.

This module performs no time integration. Three independently supplied unit-area
histories are admitted as data, their instantaneous fluxes are evaluated once per
source/sample identity, and their areal channels are assigned to exact geometric
areas. Exact sums below mean sums of represented inputs, not exact IVP solutions.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, fields
from fractions import Fraction
from hashlib import sha256
import json
import math
from numbers import Real
from types import MappingProxyType
from typing import Any

from .thermal_forced_feedback_reference import (
    FeedbackProblem, FeedbackReservoir, ForcedFeedbackHistory, ForcedFeedbackSample,
    PowerSchedule, PowerSegment, forced_feedback_rhs,
)


MODEL_IDENTITY = "manufactured-three-patch-matched-input-width-assembly-v1"
HISTORY_POWERS_W_M2 = MappingProxyType({
    "uniform_power": 60_000_000.0,
    "edge_power": 90_000_000.0,
    "middle_power": 30_000_000.0,
})
LEDGER_ROUNDING_FACTOR = 16
_KB = Fraction("1.380649e-23")
_QE = Fraction("1.602176634e-19")
_BASE = (
    "coverage_m2", "temperature_k", "released_m2", "input_energy_j_m2",
    "binding_energy_j_m2", "internal_energy_j_m2", "translational_energy_j_m2",
    "total_wall_energy_debit_j_m2",
)
_CUMULATIVE = _BASE[2:]
_SIGNED = {
    "input_energy_residual_j_m2", "wall_thermal_change_j_m2",
    "adsorbate_reference_energy_j_m2", "inventory_residual_m2",
    "adsorbate_reference_residual_j_m2", "wall_debit_residual_j_m2",
    "wall_channel_residual_j_m2", "binding_allocation_residual_j_m2",
    "internal_allocation_residual_j_m2", "total_energy_residual_j_m2",
}


class UnsupportedWidthPowerRangeError(ValueError):
    """An association, finite support, ledger or projection is unsupported."""


def _label(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _number(value: Any, name: str, *, signed: bool = False, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar, not a boolean")
    if not signed and value < 0:
        raise ValueError(f"{name} must be nonnegative before conversion")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise UnsupportedWidthPowerRangeError(f"{name} is outside finite float range") from exc
    if not math.isfinite(result) or (result == 0 and math.copysign(1.0, result) < 0):
        raise ValueError(f"{name} must be finite without negative zero")
    if value != 0 and result == 0:
        raise UnsupportedWidthPowerRangeError(f"nonzero {name} underflows")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _f(value: float) -> Fraction:
    return Fraction.from_float(value)


def _ratio(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _project(value: Fraction, name: str) -> float:
    if not value:
        return 0.0
    try:
        result = float(value)
    except OverflowError as exc:
        raise UnsupportedWidthPowerRangeError(f"{name} overflows") from exc
    if not math.isfinite(result) or result == 0:
        raise UnsupportedWidthPowerRangeError(f"nonzero {name} is not representable")
    return result


def _record(value: Fraction, name: str) -> dict[str, Any]:
    projected = _project(value, name)
    error = _f(projected) - value
    return {
        "exact": _ratio(value), "value": projected,
        "projection": {
            "exact_zero": value == 0, "positive": value > 0, "negative": value < 0,
            "finite": True, "underflowed": False, "overflowed": False,
            "rounded": error != 0,
            "subnormal": projected != 0 and abs(projected) < float.fromhex("0x1p-1022"),
            "exact_rounding_error": _ratio(error),
        },
    }


def _same_float(actual: float, expected: Fraction, name: str) -> None:
    if actual.hex() != _project(expected, name).hex():
        raise ValueError(f"scalar history {name} does not match its represented-input definition")


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def _thaw(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _thaw(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw(item) for item in value]
    return value


@dataclass(frozen=True, slots=True)
class WidthPowerGeometry:
    width_edges_m: Sequence[float]
    common_length_m: float
    parameter_status: str
    parameter_provenance: str

    def __post_init__(self) -> None:
        if self.parameter_status != "manufactured":
            raise ValueError("width geometry requires parameter_status='manufactured'")
        _label(self.parameter_provenance, "geometry parameter_provenance")
        if isinstance(self.width_edges_m, (str, bytes)) or not isinstance(self.width_edges_m, Sequence):
            raise ValueError("width_edges_m must be four canonical coordinates")
        if len(self.width_edges_m) != 4:
            raise ValueError("exactly three width patches require four edges")
        edges = tuple(_number(item, "width edge", signed=True) for item in self.width_edges_m)
        if any(raw != represented for raw, represented in zip(self.width_edges_m, edges)):
            raise UnsupportedWidthPowerRangeError("width coordinates must be exactly represented")
        if any(left >= right for left, right in zip(edges, edges[1:])):
            raise ValueError("width edges must be strictly increasing")
        length = _number(self.common_length_m, "common_length_m", positive=True)
        if self.common_length_m != length:
            raise UnsupportedWidthPowerRangeError("length must be exactly represented")
        areas = [(_f(right)-_f(left))*_f(length) for left, right in zip(edges, edges[1:])]
        total = sum(areas, Fraction())
        if tuple(area/total for area in areas) != (Fraction(1, 4), Fraction(1, 2), Fraction(1, 4)):
            raise ValueError("manufactured width area weights must be exactly 1/4,1/2,1/4")
        for area in (*areas, total):
            _project(area, "geometric area")
        object.__setattr__(self, "width_edges_m", edges)
        object.__setattr__(self, "common_length_m", length)


@dataclass(frozen=True, slots=True)
class WidthPowerAssembly:
    document: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "document", _freeze(_thaw(self.document)))

    def to_dict(self) -> dict[str, Any]:
        return _thaw(self.document)


def scalar_history_sha256(history: ForcedFeedbackHistory) -> str:
    """Canonical JSON association of a typed scalar history; no model evolution."""
    if type(history) is not ForcedFeedbackHistory:
        raise ValueError("each scalar history must be a ForcedFeedbackHistory")
    try:
        payload = ForcedFeedbackHistory.to_dict(history)
        encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    except (AttributeError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError("scalar history must have finite JSON data") from exc
    return sha256(encoded).hexdigest()


def _problem_snapshot(original: FeedbackProblem) -> FeedbackProblem:
    if type(original) is not FeedbackProblem or type(original.reservoir) is not FeedbackReservoir:
        raise ValueError("history requires explicit FeedbackProblem and FeedbackReservoir records")
    pool = FeedbackReservoir(**{field.name: getattr(original.reservoir, field.name) for field in fields(FeedbackReservoir)})
    problem = FeedbackProblem(**{
        field.name: pool if field.name == "reservoir" else getattr(original, field.name)
        for field in fields(FeedbackProblem)
    })
    for source in (original, original.reservoir):
        for field in fields(source):
            value = getattr(source, field.name)
            if isinstance(value, Real):
                _number(value, field.name)
    if problem.reservoir.activation_energy_ev <= 0:
        raise ValueError("width reference requires strictly positive base activation energy")
    if problem.patch_area_m2 != 1.0:
        raise ValueError("scalar histories must use unit reference area; scale areal fields only")
    return problem


def _schedule_snapshot(original: PowerSchedule, power: float) -> PowerSchedule:
    if type(original) is not PowerSchedule or not isinstance(original.segments, Sequence) or len(original.segments) != 1:
        raise ValueError("each width history requires one finite constant-power segment")
    if type(original.segments[0]) is not PowerSegment:
        raise ValueError("history requires a PowerSegment record")
    segment = PowerSegment(**asdict(original.segments[0]))
    schedule = PowerSchedule(original.duration_s, (segment,), original.parameter_status, original.parameter_provenance)
    if segment.power_w_m2 != power:
        raise ValueError("scalar history power does not match its assigned source identity")
    return schedule


def _extensive_key(key: str) -> str | None:
    if key.endswith("_j_m2"):
        return key[:-3]
    if key in ("coverage_m2", "released_m2", "inventory_residual_m2"):
        return key[:-3] + "_number"
    if key == "release_flux_m2_s":
        return "release_rate_s"
    return None


def _derived(problem: FeedbackProblem, sample: dict[str, float], prescribed: Fraction) -> dict[str, Fraction]:
    n, t, m, ein, binding, internal, translation, wall = (_f(sample[key]) for key in _BASE)
    d = _f(problem.reservoir.desorption_energy_ev)*_QE
    i = _f(problem.reservoir.internal_energy_ev)*_QE
    thermal = _f(problem.areal_heat_capacity_j_m2_k)*(t-_f(problem.initial_temperature_k))
    ads_change = d*m
    ads_inventory = d*(_f(problem.initial_coverage_m2)-n)
    gas = internal+translation
    return {
        "prescribed_input_energy_j_m2": prescribed,
        "input_energy_residual_j_m2": ein-prescribed,
        "gas_energy_j_m2": gas,
        "wall_thermal_change_j_m2": thermal,
        "adsorbate_reference_energy_j_m2": -d*n,
        "adsorbate_reference_change_j_m2": ads_change,
        "adsorbate_reference_inventory_change_j_m2": ads_inventory,
        "inventory_residual_m2": n+m-_f(problem.initial_coverage_m2),
        "adsorbate_reference_residual_j_m2": ads_inventory-ads_change,
        "wall_debit_residual_j_m2": thermal-ein+wall,
        "wall_channel_residual_j_m2": wall-binding-gas,
        "binding_allocation_residual_j_m2": binding-d*m,
        "internal_allocation_residual_j_m2": internal-i*m,
        "total_energy_residual_j_m2": thermal+ads_inventory+gas-ein,
    }


def _rounding_budget(values: Sequence[Fraction], steps: int) -> Fraction:
    # This admits accumulated represented-state arithmetic, not IVP truncation
    # error. Its input-relative size and purpose are recorded in the document.
    return LEDGER_ROUNDING_FACTOR*(steps+1)*sum(
        (_f(math.ulp(_project(abs(value), "ledger budget scale"))) for value in values), Fraction()
    )


def _history_snapshot(history: ForcedFeedbackHistory, power: float) -> tuple[ForcedFeedbackHistory, list[Fraction]]:
    problem = _problem_snapshot(history.problem)
    schedule = _schedule_snapshot(history.schedule, power)
    segment = schedule.segments[0]
    if type(history.samples) is not tuple or len(history.samples) != segment.steps+1:
        raise ValueError("history sample count must match the complete source step grid")
    times = [_f(schedule.duration_s)*step/segment.steps for step in range(segment.steps+1)]
    copied: list[ForcedFeedbackSample] = []
    previous: dict[str, float] | None = None
    for step, (source, time) in enumerate(zip(history.samples, times)):
        if type(source) is not ForcedFeedbackSample or type(source.segment_index) is not int or type(source.segment_step) is not int:
            raise ValueError("samples require typed records and strict integer source indices")
        if (source.segment_index, source.segment_step) != (0, step):
            raise ValueError("scalar sample indices do not match the complete logical grid")
        if type(source.time_s) is not float:
            raise ValueError("sample time must be a canonical float")
        _number(source.time_s, "sample time")
        _same_float(source.time_s, time, "sample time")
        values: dict[str, float] = {}
        for field in fields(ForcedFeedbackSample):
            if field.name in ("segment_index", "segment_step", "time_s", "extensive"):
                continue
            value = getattr(source, field.name)
            if type(value) is not float:
                raise ValueError(f"scalar sample {field.name} must be a finite float")
            values[field.name] = _number(value, field.name, signed=field.name in _SIGNED)
        n, t = values["coverage_m2"], values["temperature_k"]
        if not 0 <= n <= problem.reservoir.max_coverage_m2 or not problem.temperature_min_k <= t <= problem.temperature_max_k:
            raise ValueError("scalar state leaves its supplied coverage or temperature support")
        if values["released_m2"] > problem.initial_coverage_m2:
            raise ValueError("cumulative release exceeds initial inventory")
        if step == 0:
            if n != problem.initial_coverage_m2 or t != problem.initial_temperature_k or any(values[key] != 0 for key in _CUMULATIVE):
                raise ValueError("initial scalar snapshot does not match its supplied problem")
        elif previous is not None:
            if n > previous["coverage_m2"] or any(values[key] < previous[key] for key in _CUMULATIVE):
                raise ValueError("scalar coverage or gross channels violate monotonicity")
        prescribed = _f(power)*time
        derived = _derived(problem, values, prescribed)
        for name, exact in derived.items():
            _same_float(values[name], exact, name)
        number_budget = _rounding_budget([_f(problem.initial_coverage_m2), _f(n), _f(values["released_m2"])], segment.steps)
        energy_scales = [_f(values[key]) for key in _CUMULATIVE if key.endswith("_j_m2")]
        energy_scales.extend([
            _f(problem.areal_heat_capacity_j_m2_k)*_f(problem.initial_temperature_k),
            _f(problem.areal_heat_capacity_j_m2_k)*_f(t),
            _f(problem.reservoir.desorption_energy_ev)*_QE*_f(problem.initial_coverage_m2),
            _f(problem.reservoir.internal_energy_ev)*_QE*_f(problem.initial_coverage_m2),
        ])
        energy_budget = _rounding_budget(energy_scales, segment.steps)
        if abs(derived["inventory_residual_m2"]) > number_budget:
            raise ValueError("scalar inventory residual exceeds its arithmetic admission budget")
        for name, exact in derived.items():
            if name.endswith("_residual_j_m2") and abs(exact) > energy_budget:
                raise ValueError(f"scalar {name} exceeds its arithmetic admission budget")
        expected_extensive = {_extensive_key(key): value for key, value in values.items() if _extensive_key(key) is not None}
        if not isinstance(source.extensive, Mapping) or set(source.extensive) != set(expected_extensive):
            raise ValueError("unit-area scalar extensive fields are incomplete or unexpected")
        for key, expected in expected_extensive.items():
            actual = source.extensive[key]
            if type(actual) is not float or actual.hex() != expected.hex():
                raise ValueError("unit-area scalar extensive value does not match its areal source")
        copied.append(ForcedFeedbackSample(0, step, source.time_s, **values, extensive=dict(source.extensive)))
        previous = values
    return ForcedFeedbackHistory(problem, schedule, tuple(copied)), times


def _aggregate_residuals(problem: FeedbackProblem, values: dict[str, Fraction], area: Fraction) -> dict[str, Fraction]:
    d = _f(problem.reservoir.desorption_energy_ev)*_QE
    i = _f(problem.reservoir.internal_energy_ev)*_QE
    n, m = values["coverage_number"], values["released_number"]
    binding, internal = values["binding_energy_j"], values["internal_energy_j"]
    translation, wall = values["translational_energy_j"], values["total_wall_energy_debit_j"]
    thermal, ein = values["wall_thermal_change_j"], values["input_energy_j"]
    ads_inventory = d*(_f(problem.initial_coverage_m2)*area-n)
    return {
        "inventory_residual_number": n+m-_f(problem.initial_coverage_m2)*area,
        "input_energy_residual_j": ein-values["prescribed_input_energy_j"],
        "adsorbate_reference_residual_j": ads_inventory-d*m,
        "wall_debit_residual_j": thermal-ein+wall,
        "wall_channel_residual_j": wall-binding-internal-translation,
        "binding_allocation_residual_j": binding-d*m,
        "internal_allocation_residual_j": internal-i*m,
        "total_energy_residual_j": thermal+ads_inventory+internal+translation-ein,
    }


def assemble_width_power_reference(
    histories: Mapping[str, ForcedFeedbackHistory], *, geometry: WidthPowerGeometry,
    hypotheses: Mapping[str, Sequence[str]], history_sha256: Mapping[str, str],
    history_provenance: Mapping[str, str], parameter_status: str, parameter_provenance: str,
) -> WidthPowerAssembly:
    """Admit and assemble three distinct supplied histories without re-evolution."""
    if parameter_status != "manufactured":
        raise ValueError("width assembly requires parameter_status='manufactured'")
    _label(parameter_provenance, "assembly parameter_provenance")
    expected_keys = set(HISTORY_POWERS_W_M2)
    for name, mapping in (("histories", histories), ("history_sha256", history_sha256), ("history_provenance", history_provenance)):
        if not isinstance(mapping, Mapping) or set(mapping) != expected_keys:
            raise ValueError(f"{name} must contain exactly the three declared scalar source identities")
    originals = dict(histories)
    hashes = dict(history_sha256)
    provenance = {key: _label(value, "history provenance") for key, value in history_provenance.items()}
    expected_hypotheses = {
        "uniform": ("uniform_power",)*3,
        "nonuniform": ("edge_power", "middle_power", "edge_power"),
    }
    if not isinstance(hypotheses, Mapping) or set(hypotheses) != set(expected_hypotheses):
        raise ValueError("both explicit uniform and nonuniform hypotheses are required")
    assignments = {}
    for name, expected in expected_hypotheses.items():
        raw = hypotheses[name]
        if isinstance(raw, (str, bytes)) or not isinstance(raw, Sequence) or tuple(raw) != expected:
            raise ValueError("width hypothesis source assignments do not match the declared experiment")
        assignments[name] = tuple(raw)
    if type(geometry) is not WidthPowerGeometry:
        raise ValueError("geometry must be a WidthPowerGeometry record")
    geometry = WidthPowerGeometry(geometry.width_edges_m, geometry.common_length_m, geometry.parameter_status, geometry.parameter_provenance)
    edges = geometry.width_edges_m
    areas = [(_f(right)-_f(left))*_f(geometry.common_length_m) for left, right in zip(edges, edges[1:])]
    total_area = sum(areas, Fraction())
    admitted: dict[str, ForcedFeedbackHistory] = {}
    times_by_key: dict[str, list[Fraction]] = {}
    for key, power in HISTORY_POWERS_W_M2.items():
        digest = hashes[key]
        if type(digest) is not str or len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            raise ValueError("history_sha256 must supply canonical lower-case SHA256 strings")
        if scalar_history_sha256(originals[key]) != digest:
            raise ValueError("scalar history source hash association mismatch")
        admitted[key], times_by_key[key] = _history_snapshot(originals[key], power)
    common = admitted["uniform_power"]
    common_problem = asdict(common.problem)
    common_times = times_by_key["uniform_power"]
    for key, history in admitted.items():
        if asdict(history.problem) != common_problem or times_by_key[key] != common_times:
            raise ValueError("histories must share exact problem/state/cost/provenance and complete logical grid")
        schedule = history.schedule
        if schedule.parameter_status != common.schedule.parameter_status or schedule.parameter_provenance != common.schedule.parameter_provenance:
            raise ValueError("power schedule status and common provenance must match")
        if scalar_history_sha256(history) != hashes[key]:
            raise ValueError("admitted snapshot changed its scalar source association")
    power_sums = {
        name: sum((area*_f(HISTORY_POWERS_W_M2[key]) for area, key in zip(areas, keys)), Fraction())
        for name, keys in assignments.items()
    }
    if power_sums["uniform"] != power_sums["nonuniform"]:
        raise ValueError("prescribed total input powers must match exactly without normalization")

    fluxes: dict[str, list[float]] = {}
    rhs_records: dict[str, list[dict[str, float]]] = {}
    rhs_evaluations = 0
    for key, history in admitted.items():
        rates = []
        cached_rhs = []
        for sample in history.samples:
            rhs = forced_feedback_rhs(history.problem, coverage_m2=sample.coverage_m2,
                                      temperature_k=sample.temperature_k,
                                      power_w_m2=HISTORY_POWERS_W_M2[key])
            expected_rhs = {
                "release_flux_m2_s", "coverage_rate_m2_s", "temperature_rate_k_s",
                "input_power_j_m2_s", "binding_power_j_m2_s", "internal_power_j_m2_s",
                "translational_power_j_m2_s", "total_wall_power_j_m2_s",
            }
            if not isinstance(rhs, Mapping) or set(rhs) != expected_rhs:
                raise ValueError("complete instantaneous RHS rate records are required")
            for quantity, value in rhs.items():
                if type(value) is not float:
                    raise ValueError("RHS rates must be finite floats")
                _number(value, quantity, signed=quantity in ("coverage_rate_m2_s", "temperature_rate_k_s"))
            rate = rhs["release_flux_m2_s"]
            if (sample.coverage_m2 > 0) != (rate > 0):
                raise UnsupportedWidthPowerRangeError("instantaneous flux loses positive inventory support")
            rates.append(rate)
            cached_rhs.append(dict(rhs))
            rhs_evaluations += 1
        fluxes[key] = rates
        rhs_records[key] = cached_rhs

    hypothesis_docs = {}
    totals_exact: dict[str, list[dict[str, Fraction]]] = {}
    for name, keys in assignments.items():
        patches = []
        accumulated: list[dict[str, Fraction]] = [{} for _ in common_times]
        for patch_index, (key, area) in enumerate(zip(keys, areas)):
            history = admitted[key]
            patch_samples = []
            for sample_index, sample in enumerate(history.samples):
                areal = {field.name: getattr(sample, field.name) for field in fields(ForcedFeedbackSample)
                         if field.name not in ("segment_index", "segment_step", "time_s", "extensive")}
                areal["release_flux_m2_s"] = fluxes[key][sample_index]
                exact_values = {_extensive_key(quantity): _f(value)*area for quantity, value in areal.items()
                                if _extensive_key(quantity) is not None}
                for quantity, value in exact_values.items():
                    accumulated[sample_index][quantity] = accumulated[sample_index].get(quantity, Fraction())+value
                patch_samples.append({
                    "segment_index": sample.segment_index, "segment_step": sample.segment_step,
                    "time_s": sample.time_s, "logical_time_s": _ratio(common_times[sample_index]),
                    "state_id": history.problem.reservoir.state_id, "scalar_history_id": key,
                    "areal": areal,
                    "extensive": {quantity: _record(value, f"patch {quantity}") for quantity, value in exact_values.items()},
                })
            patches.append({
                "patch_index": patch_index, "width_left_m": edges[patch_index], "width_right_m": edges[patch_index+1],
                "area_m2": _project(area, "patch area"), "area_exact_m2": _ratio(area),
                "state_id": history.problem.reservoir.state_id, "scalar_history_id": key,
                "history_sha256": hashes[key], "samples": patch_samples,
            })
        totals = []
        for index, (values, time) in enumerate(zip(accumulated, common_times)):
            sample = common.samples[index]
            residuals = _aggregate_residuals(common.problem, values, total_area)
            totals.append({
                "segment_index": sample.segment_index, "segment_step": sample.segment_step,
                "time_s": sample.time_s, "logical_time_s": _ratio(time),
                "state_id": common.problem.reservoir.state_id,
                "extensive": {key: _record(value, f"total {key}") for key, value in values.items()},
                "recomputed_residuals": {key: _record(value, f"aggregate {key}") for key, value in residuals.items()},
            })
        totals_exact[name] = accumulated
        hypothesis_docs[name] = {"patches": patches, "totals": totals}
    comparison_samples = []
    for index, time in enumerate(common_times):
        uniform = totals_exact["uniform"][index]
        nonuniform = totals_exact["nonuniform"][index]
        comparison_samples.append({
            "time_s": common.samples[index].time_s, "logical_time_s": _ratio(time),
            "difference_nonuniform_minus_uniform": {
                key: _record(nonuniform[key]-uniform[key], f"hypothesis difference {key}") for key in uniform
            },
            "exact_prescribed_input_energy_j": _record(power_sums["uniform"]*time, "exact common prescribed input"),
            "integrated_input_difference_j": _record(nonuniform["input_energy_j"]-uniform["input_energy_j"], "integrated input difference"),
            "prescribed_input_difference_j": _record((power_sums["nonuniform"]-power_sums["uniform"])*time, "exact prescribed input difference"),
        })
    for key, original in originals.items():
        if scalar_history_sha256(original) != hashes[key]:
            raise ValueError("original scalar history changed during width assembly")
    document = {
        "schema": "mhdlab-thermal-width-power-reference-v1", "model_identity": MODEL_IDENTITY,
        "parameter_status": parameter_status, "parameter_provenance": parameter_provenance,
        "geometry": {
            "width_edges_m": list(edges), "width_edges_exact_m": [_ratio(_f(edge)) for edge in edges],
            "common_length_m": geometry.common_length_m, "common_length_exact_m": _ratio(_f(geometry.common_length_m)),
            "patch_areas_m2": [_project(area, "patch area") for area in areas],
            "patch_areas_exact_m2": [_ratio(area) for area in areas],
            "total_area_m2": _project(total_area, "total area"), "total_area_exact_m2": _ratio(total_area),
            "area_weights": [_ratio(area/total_area) for area in areas],
            "parameter_status": geometry.parameter_status, "parameter_provenance": geometry.parameter_provenance,
        },
        "scalar_histories": {
            key: {"history_sha256": hashes[key], "parameter_provenance": provenance[key],
                  "history": history.to_dict(), "logical_time_s": [_ratio(time) for time in times_by_key[key]],
                  "release_flux_m2_s": fluxes[key], "rhs": rhs_records[key]}
            for key, history in admitted.items()
        },
        "hypotheses": hypothesis_docs,
        "comparison": {"exact_total_input_power_w": _record(power_sums["uniform"], "common input power"),
                       "samples": comparison_samples},
        "cache": {"distinct_scalar_histories": 3, "scalar_samples": sum(len(history.samples) for history in admitted.values()),
                  "rhs_evaluations": rhs_evaluations, "patch_history_references": 6, "scalar_evolutions_by_assembly": 0,
                  "key": "scalar_history_id plus saved sample index; reused patches share the same evaluated flux"},
        "numerical_contract": {
            "exact_arithmetic": "Fractions of represented areal values and canonical geometry; projected once for each patch/total/difference record.",
            "arithmetic_admission": f"Derived scalar fields must match exactly; inventory/energy residuals must fit {LEDGER_ROUNDING_FACTOR}*(steps+1) times summed relevant represented-scale ulps. This is not an IVP truncation-error or physical-uncertainty certificate.",
            "flux": "Existing positive-barrier RHS numerical contract retained; no new rate approximation or time evolution.",
            "source_admission": "Canonical SHA256 associations, common supplied data and finite-grid/ledger consistency; independent numerical references must establish trajectory accuracy.",
            "projection": "Unsupported nonzero underflow/overflow rejects; no clipping, renormalization or repair.",
        },
        "host_budget_accepted": False, "owner_implementation_verified": False,
        "limits": [
            "Manufactured independent stationary width patches only; no lateral conduction, phase change, moving boundaries, radiation or field-derived power.",
            "Prescribed input equality is exact. Separately integrated input channels and their rounding differences remain visible.",
            "Full per-patch temperature/time/state mixtures are retained; no mean-temperature or single-velocity substitution is made.",
            "Signs and magnitudes of nonlinear release and energy differences are outcomes, not assumptions.",
            "No physical water inventory, actual foil geometry/material kinetics, outgoing-state chemistry, velocity launch, spectra, FLASH, Hall, gap or owner acceptance is supplied.",
        ],
    }
    return WidthPowerAssembly(document)
