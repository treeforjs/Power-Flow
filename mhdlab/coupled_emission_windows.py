"""Positive polynomial emission windows from explicitly admitted saved data.

No history evolution, rate evaluation, assembly or transport module is imported.
Exact arithmetic describes a declared interpolant of represented source values;
it does not make that interpolant an exact kinetic trajectory.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, fields
from fractions import Fraction
from hashlib import sha256
import json
import math
import re
from types import MappingProxyType
from typing import Any


MODEL_IDENTITY = "manufactured-constant-power-positive-coupled-emission-windows-v1"
MAX_JSON_BYTES = 64*1024*1024
MAX_TOTAL_JSON_BYTES = 256*1024*1024
MAX_INTEGER_DIGITS = 4096
MAX_FRACTION_BITS = 65536
MAX_CELLS = 64
_KB = Fraction("1.380649e-23")
_QE = Fraction("1.602176634e-19")
_POWERS = {"uniform_power": 60_000_000.0, "edge_power": 90_000_000.0, "middle_power": 30_000_000.0}
_ASSIGNMENTS = {"uniform": ["uniform_power"]*3, "nonuniform": ["edge_power", "middle_power", "edge_power"]}
_ALIASES = {"report", "assembly", "input_config", "reuse_map", *_POWERS,
            "independent_completion", "independent_report", "scientific_review"}
_RAW_MOMENTS = (
    "released_m2", "local_time_first_moment_s_m2", "local_time_second_moment_s2_m2",
    "temperature_first_moment_k_m2", "temperature_second_moment_k2_m2",
    "local_time_temperature_moment_k_s_m2",
)
_RHS_FIELDS = {"release_flux_m2_s", "coverage_rate_m2_s", "temperature_rate_k_s",
               "input_power_j_m2_s", "binding_power_j_m2_s", "internal_power_j_m2_s",
               "translational_power_j_m2_s", "total_wall_power_j_m2_s"}
_SOURCE_CODE_SHA256 = {
    "tools/run_thermal_width_power_reference.py": "57e4717b48920d1657f3e4ad156bad7a72066c96137b971d47552570475ff72b",
    "mhdlab/thermal_width_power_reference.py": "ed796e54e9a46a9981459aabc5a09d4da28479f776ef4ce601f85e108cba930c",
    "mhdlab/thermal_forced_feedback_reference.py": "16033899c949896ca346eb6cfb6beb9cdd9d9253c1a110fa9f057a9b8bdeac4d",
    "mhdlab/thermal_power_reference.py": "741ab61b45cc0440b2f133732df387e34bd7630d6b6d5eeb48be7f6e471f1e00",
    "mhdlab/thermal_feedback_reference.py": "3fbded5bc801fac4e847195ecd5ad369142d3ab56e70848c4501ea1425d61cf8",
}


class UnsupportedCoupledWindowRangeError(ValueError):
    """Unsupported source, positive interpolant, finite support or projection."""


def _require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError(reason)


def _label(value: Any, name: str) -> str:
    _require(isinstance(value, str) and bool(value.strip()), f"{name} must be a nonempty string")
    return value


def _bounded(value: Fraction) -> Fraction:
    if max(value.numerator.bit_length(), value.denominator.bit_length()) > MAX_FRACTION_BITS:
        raise UnsupportedCoupledWindowRangeError("exact rational work exceeds the declared bit limit")
    return value


def _exact(value: Any, name: str, *, signed: bool = False, positive: bool = False) -> Fraction:
    if type(value) is int or type(value) is Fraction:
        result = Fraction(value)
    elif type(value) is float:
        _require(math.isfinite(value) and not (value == 0 and math.copysign(1.0, value) < 0), f"{name} must be finite without negative zero")
        result = Fraction.from_float(value)
    else:
        raise ValueError(f"{name} requires a Fraction, integer or finite represented float; booleans reject")
    _bounded(result)
    _require(signed or result >= 0, f"{name} must be nonnegative")
    _require(not positive or result > 0, f"{name} must be positive")
    return result


def _ratio(value: Fraction) -> dict[str, str]:
    _bounded(value)
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _from_ratio(value: Any) -> Fraction:
    _require(isinstance(value, Mapping) and set(value) == {"numerator", "denominator"}, "canonical rational record required")
    n, d = value["numerator"], value["denominator"]
    _require(type(n) is str and type(d) is str and len(n) <= MAX_INTEGER_DIGITS and len(d) <= MAX_INTEGER_DIGITS, "bounded rational integer strings required")
    _require(re.fullmatch(r"0|-?[1-9][0-9]*", n) is not None and re.fullmatch(r"[1-9][0-9]*", d) is not None, "noncanonical rational integers")
    result = _bounded(Fraction(int(n), int(d)))
    _require(_ratio(result) == dict(value), "rational must be reduced")
    return result


def _project(value: Fraction, name: str) -> float:
    _bounded(value)
    if not value:
        return 0.0
    try:
        result = float(value)
    except OverflowError as exc:
        raise UnsupportedCoupledWindowRangeError(f"{name} overflows") from exc
    if not math.isfinite(result) or result == 0:
        raise UnsupportedCoupledWindowRangeError(f"nonzero {name} is not representable")
    return result


def _record(value: Fraction, name: str) -> dict[str, Any]:
    result = _project(value, name)
    error = Fraction.from_float(result)-value
    return {"exact": _ratio(value), "value": result, "projection": {
        "exact_zero": value == 0, "positive": value > 0, "negative": value < 0,
        "finite": True, "underflowed": False, "overflowed": False,
        "rounded": error != 0, "subnormal": result != 0 and abs(result) < float.fromhex("0x1p-1022"),
        "exact_rounding_error": _ratio(error)}}


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _thaw(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _thaw(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw(item) for item in value]
    return value


def _same(actual: Any, expected: Any, name: str) -> None:
    if type(actual) is not type(expected):
        raise ValueError(f"{name} has a type/association mismatch")
    if type(actual) is float:
        equal = actual.hex() == expected.hex()
    elif isinstance(actual, dict):
        equal = set(actual) == set(expected)
        if equal:
            for key in actual:
                _same(actual[key], expected[key], name)
    elif isinstance(actual, list):
        equal = len(actual) == len(expected)
        if equal:
            for left, right in zip(actual, expected):
                _same(left, right, name)
    else:
        equal = actual == expected
    _require(equal, f"{name} association mismatch")


def _json_sha(value: Any) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _sha(value: Any) -> str:
    _require(type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value) is not None, "expected SHA256 must be lower-case hexadecimal")
    return value


def _strict_json(raw: bytes) -> dict[str, Any]:
    def pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
        result = {}
        for key, value in items:
            _require(key not in result, "duplicate JSON object key")
            result[key] = value
        return result

    def integer(value: str) -> int:
        _require(len(value.lstrip("-")) <= MAX_INTEGER_DIGITS, "JSON integer exceeds digit limit")
        _require(value != "-0", "negative-zero JSON integer is not canonical")
        return int(value)

    def floating(value: str) -> float:
        _require(len(value) <= MAX_INTEGER_DIGITS, "JSON float token exceeds digit limit")
        result = float(value)
        _require(math.isfinite(result), "nonfinite JSON numeric literal")
        mantissa = value.lower().split("e", 1)[0]
        raw_nonzero = any(character in "123456789" for character in mantissa)
        _require(not (raw_nonzero and result == 0), "nonzero JSON float underflows its representation")
        return result

    def bad_constant(value: str) -> None:
        raise ValueError(f"nonfinite JSON constant {value}")

    def walk(value: Any, depth: int = 0) -> None:
        _require(depth <= 128, "JSON nesting exceeds limit")
        if type(value) is float:
            _require(math.isfinite(value) and not (value == 0 and math.copysign(1.0, value) < 0), "nonfinite or negative-zero JSON float")
        elif isinstance(value, dict):
            for item in value.values():
                walk(item, depth+1)
        elif isinstance(value, list):
            for item in value:
                walk(item, depth+1)

    try:
        result = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs, parse_int=integer,
                            parse_float=floating, parse_constant=bad_constant)
        _require(type(result) is dict, "each saved file must contain a JSON object")
        walk(result)
    except (UnicodeError, RecursionError, OverflowError) as exc:
        raise ValueError("unsupported JSON encoding, nesting or numeric range") from exc
    return result


@dataclass(frozen=True, slots=True)
class CoupledEmissionCell:
    start_s: Any
    stop_s: Any
    emitted_controls_m2: Sequence[Any]
    temperature_controls_k: Sequence[Any]
    temperature_min_k: Any
    temperature_max_k: Any
    initial_coverage_m2: Any
    state_id: str
    source_history_id: str
    parameter_status: str
    parameter_provenance: str

    def __post_init__(self) -> None:
        start = _exact(self.start_s, "cell start")
        stop = _exact(self.stop_s, "cell stop", positive=True)
        _require(start < stop, "cell requires positive finite support")
        n0 = _exact(self.initial_coverage_m2, "initial coverage")
        tmin = _exact(self.temperature_min_k, "minimum temperature", positive=True)
        tmax = _exact(self.temperature_max_k, "maximum temperature", positive=True)
        _require(tmin < tmax, "temperature support must have positive width")
        controls = []
        for name, raw in (("emitted_controls_m2", self.emitted_controls_m2), ("temperature_controls_k", self.temperature_controls_k)):
            _require(isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)) and len(raw) == 4, "four cubic Bernstein controls required")
            controls.append(tuple(_exact(item, name) for item in raw))
        mass, temperature = controls
        _require(all(left <= right for left, right in zip(mass, mass[1:])) and 0 <= mass[0] <= mass[-1] <= n0, "emission Bernstein controls must be nondecreasing within finite inventory")
        _require(all(tmin <= item <= tmax for item in temperature), "temperature Bernstein controls leave supplied support")
        _require(self.parameter_status == "manufactured", "cell requires parameter_status='manufactured'")
        for key in ("state_id", "source_history_id", "parameter_provenance"):
            _label(getattr(self, key), key)
        for name, value in (("start_s", start), ("stop_s", stop), ("initial_coverage_m2", n0),
                            ("temperature_min_k", tmin), ("temperature_max_k", tmax),
                            ("emitted_controls_m2", mass), ("temperature_controls_k", temperature)):
            object.__setattr__(self, name, value)


@dataclass(frozen=True, slots=True)
class WindowQuery:
    window_id: str
    start_fraction: Any
    stop_fraction: Any
    parameter_provenance: str

    def __post_init__(self) -> None:
        _label(self.window_id, "window_id")
        _label(self.parameter_provenance, "query provenance")
        start = _exact(self.start_fraction, "query start fraction")
        stop = _exact(self.stop_fraction, "query stop fraction", positive=True)
        _require(0 <= start < stop <= 1, "query must lie inside finite canonical source duration")
        object.__setattr__(self, "start_fraction", start)
        object.__setattr__(self, "stop_fraction", stop)


@dataclass(frozen=True, slots=True)
class AdmittedCoupledSource:
    assembly: Mapping[str, Any]
    histories: Mapping[str, Any]
    admission: Mapping[str, Any]

    def __post_init__(self) -> None:
        for key in ("assembly", "histories", "admission"):
            object.__setattr__(self, key, _freeze(getattr(self, key)))

    def to_dict(self) -> dict[str, Any]:
        return _thaw(self.admission)


@dataclass(frozen=True, slots=True)
class CoupledEmissionWindows:
    document: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "document", _freeze(self.document))

    def to_dict(self) -> dict[str, Any]:
        return _thaw(self.document)


def _split(controls: Sequence[Fraction], point: Fraction) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    rows = [tuple(controls)]
    while len(rows[-1]) > 1:
        rows.append(tuple(_bounded((1-point)*left+point*right) for left, right in zip(rows[-1], rows[-1][1:])))
    return tuple(row[0] for row in rows), tuple(row[-1] for row in reversed(rows))


def _restrict(controls: Sequence[Fraction], start: Fraction, stop: Fraction) -> tuple[Fraction, ...]:
    _require(0 <= start < stop <= 1, "invalid polynomial restriction")
    left, _ = _split(controls, stop)
    _, result = _split(left, start/stop)
    return result


def _product(left: Sequence[Fraction], right: Sequence[Fraction]) -> tuple[Fraction, ...]:
    m, n = len(left)-1, len(right)-1
    _require(m+n <= 8, "polynomial degree exceeds bounded moment scope")
    result = [Fraction(0)]*(m+n+1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] = _bounded(result[i+j]+a*b*Fraction(math.comb(m, i)*math.comb(n, j), math.comb(m+n, i+j)))
    return tuple(result)


def _integral(*polynomials: Sequence[Fraction]) -> Fraction:
    coefficients = (Fraction(1),)
    for polynomial in polynomials:
        coefficients = _product(coefficients, polynomial)
    return _bounded(sum(coefficients, Fraction())/len(coefficients))


def _statistics(raw: Mapping[str, Fraction]) -> dict[str, Any]:
    mass = raw["released_m2"]
    keys = ("mean_local_time_s", "variance_local_time_s2", "mean_temperature_k", "variance_temperature_k2",
            "mean_local_time_temperature_k_s", "covariance_local_time_temperature_k_s")
    if mass == 0:
        return {key: {"available": False, "reason": "exact_zero_emitted_mass", "record": None} for key in keys}
    _require(mass > 0, "negative emitted-number measure")
    mt = raw["local_time_first_moment_s_m2"]/mass
    mtemp = raw["temperature_first_moment_k_m2"]/mass
    vt = raw["local_time_second_moment_s2_m2"]/mass-mt*mt
    vtemp = raw["temperature_second_moment_k2_m2"]/mass-mtemp*mtemp
    joint = raw["local_time_temperature_moment_k_s_m2"]/mass
    covariance = joint-mt*mtemp
    _require(vt >= 0 and vtemp >= 0 and covariance*covariance <= vt*vtemp, "exact moment covariance is not positive semidefinite")
    values = (mt, vt, mtemp, vtemp, joint, covariance)
    return {key: {"available": True, "reason": None, "record": _record(value, key)} for key, value in zip(keys, values)}


def integrate_emission_window(
    cells: Sequence[CoupledEmissionCell], *, start_s: Any, stop_s: Any,
    desorption_energy_ev: Any, internal_energy_ev: Any, areal_heat_capacity_j_m2_k: Any,
    power_w_m2: Any, query_id: str, parameter_provenance: str,
) -> dict[str, Any]:
    """Integrate an admitted positive piecewise polynomial measure exactly."""
    _label(query_id, "query_id")
    _label(parameter_provenance, "window provenance")
    _require(isinstance(cells, Sequence) and not isinstance(cells, (str, bytes)) and 1 <= len(cells) <= MAX_CELLS, "one to64 source cells required")
    copied = []
    for cell in cells:
        _require(type(cell) is CoupledEmissionCell, "typed CoupledEmissionCell required")
        copied.append(CoupledEmissionCell(**{field.name: getattr(cell, field.name) for field in fields(CoupledEmissionCell)}))
    first = copied[0]
    for cell in copied:
        for name in ("state_id", "source_history_id", "temperature_min_k", "temperature_max_k", "initial_coverage_m2", "parameter_status", "parameter_provenance"):
            _require(getattr(cell, name) == getattr(first, name), "source cells must share identity/domain/provenance")
    for left, right in zip(copied, copied[1:]):
        _require(left.stop_s == right.start_s and left.emitted_controls_m2[-1] == right.emitted_controls_m2[0]
                 and left.temperature_controls_k[-1] == right.temperature_controls_k[0], "source cell support or endpoints are discontinuous")
        for name in ("emitted_controls_m2", "temperature_controls_k"):
            a, b = getattr(left, name), getattr(right, name)
            _require((a[-1]-a[-2])/(left.stop_s-left.start_s) == (b[1]-b[0])/(right.stop_s-right.start_s), "constant-power endpoint derivative association is discontinuous")
    start = _exact(start_s, "window start")
    stop = _exact(stop_s, "window stop", positive=True)
    _require(first.start_s <= start < stop <= copied[-1].stop_s, "window may not extrapolate beyond finite source support")
    d = _exact(desorption_energy_ev, "removal cost")*_QE
    internal = _exact(internal_energy_ev, "internal cost")*_QE
    capacity = _exact(areal_heat_capacity_j_m2_k, "areal heat capacity", positive=True)
    power = _exact(power_w_m2, "constant power")
    raw = {key: Fraction(0) for key in _RAW_MOMENTS}
    pieces = []
    temperature_start = temperature_stop = None
    for index, cell in enumerate(copied):
        a, b = max(start, cell.start_s), min(stop, cell.stop_s)
        if a >= b:
            continue
        h = cell.stop_s-cell.start_s
        ua, ub = (a-cell.start_s)/h, (b-cell.start_s)/h
        m = _restrict(cell.emitted_controls_m2, ua, ub)
        t = _restrict(cell.temperature_controls_k, ua, ub)
        dm = tuple(3*(right-left) for left, right in zip(m, m[1:]))
        _require(all(value >= 0 for value in dm), "restriction lost the positive emission measure")
        local = (a-start, b-start)
        moments = (_integral(dm), _integral(local, dm), _integral(local, local, dm),
                   _integral(t, dm), _integral(t, t, dm), _integral(local, t, dm))
        _require(moments[0] == m[-1]-m[0], "exact polynomial mass identity failed")
        for key, value in zip(_RAW_MOMENTS, moments):
            raw[key] = _bounded(raw[key]+value)
        if temperature_start is None:
            temperature_start = t[0]
        temperature_stop = t[-1]
        pieces.append({
            "source_cell_index": index, "source_start_exact_s": _ratio(cell.start_s), "source_stop_exact_s": _ratio(cell.stop_s),
            "start_exact_s": _ratio(a), "stop_exact_s": _ratio(b), "source_u_start": _ratio(ua), "source_u_stop": _ratio(ub),
            "source_emitted_controls_m2": [_ratio(value) for value in cell.emitted_controls_m2],
            "source_temperature_controls_k": [_ratio(value) for value in cell.temperature_controls_k],
            "emitted_controls_m2": [_ratio(value) for value in m], "temperature_controls_k": [_ratio(value) for value in t],
            "emitted_derivative_controls_m2": [_ratio(value) for value in dm],
        })
    _require(temperature_start is not None and temperature_stop is not None, "window has no admitted pieces")
    mass = raw["released_m2"]
    binding, gas_internal = d*mass, internal*mass
    translation = 2*_KB*raw["temperature_first_moment_k_m2"]
    thermal = capacity*(temperature_stop-temperature_start)
    ein = power*(stop-start)
    raw.update({"binding_energy_j_m2": binding, "internal_energy_j_m2": gas_internal,
                "translational_energy_j_m2": translation, "total_wall_energy_debit_j_m2": binding+gas_internal+translation,
                "input_energy_j_m2": ein, "wall_thermal_change_j_m2": thermal,
                "adsorbate_reference_change_j_m2": binding, "gas_energy_j_m2": gas_internal+translation,
                "total_energy_residual_j_m2": thermal+binding+gas_internal+translation-ein})
    return {"query_id": query_id, "source_history_id": first.source_history_id, "state_id": first.state_id,
            "window_start_exact_s": _ratio(start), "window_stop_exact_s": _ratio(stop), "local_time_anchor_exact_s": _ratio(start),
            "pieces": pieces, "areal": {key: _record(value, key) for key, value in raw.items()},
            "statistics": _statistics(raw), "parameter_provenance": parameter_provenance}


def _strict_int(value: Any, expected: int, name: str) -> None:
    _require(type(value) is int and value == expected, f"{name} strict integer association mismatch")


def _source_float(value: Any, name: str, *, nonnegative: bool = False) -> Fraction:
    _require(type(value) is float, f"saved {name} must be a represented float")
    return _exact(value, name, signed=not nonnegative)


def _parameter_same(actual: Any, expected: Any, name: str) -> None:
    # Raw configuration integers may project to canonical scalar floats. This
    # allowance is confined to the numerical parameter comparison, never IDs.
    if type(expected) is float:
        _require(type(actual) in (int, float), f"{name} numeric parameter type mismatch")
        _require(_exact(actual, name, signed=True) == _exact(expected, name, signed=True), f"{name} parameter mismatch")
    elif type(expected) is dict:
        _require(type(actual) is dict and set(actual) == set(expected), f"{name} parameter shape mismatch")
        for key in expected:
            _parameter_same(actual[key], expected[key], name)
    else:
        _same(actual, expected, name)


def _source_extensive_key(key: str) -> str | None:
    if key.endswith("_j_m2"):
        return key[:-3]
    if key in ("coverage_m2", "released_m2", "inventory_residual_m2"):
        return key[:-3]+"_number"
    return None


def admit_saved_width_source(
    files: Mapping[str, bytes], *, expected_sha256: Mapping[str, str], parameter_provenance: str,
) -> AdmittedCoupledSource:
    """Admit explicitly hash-pinned source/acceptance bytes, without old loaders.

    The caller supplies trusted expected hashes and owns filesystem path/copy
    safety. All raw hashes are checked before any JSON graph is interpreted.
    """
    _label(parameter_provenance, "admission provenance")
    _require(isinstance(files, Mapping) and set(files) == _ALIASES, "all ten saved-source file aliases are required")
    _require(isinstance(expected_sha256, Mapping) and set(expected_sha256) == _ALIASES, "all ten expected file hashes are required")
    captured = dict(files)
    hashes = {key: _sha(value) for key, value in expected_sha256.items()}
    total_bytes = 0
    for key, raw in captured.items():
        _require(type(raw) is bytes and 0 < len(raw) <= MAX_JSON_BYTES, "saved JSON bytes exceed per-file support")
        total_bytes += len(raw)
        _require(sha256(raw).hexdigest() == hashes[key], f"raw {key} SHA256 mismatch")
    _require(total_bytes <= MAX_TOTAL_JSON_BYTES, "saved source bundle exceeds total byte limit")
    docs = {key: _strict_json(raw) for key, raw in captured.items()}
    report, assembly = docs["report"], docs["assembly"]
    config, reuse = docs["input_config"], docs["reuse_map"]
    independent, reference, review = (docs[key] for key in ("independent_completion", "independent_report", "scientific_review"))
    _require(report.get("schema") == "mhdlab-thermal-width-power-reference-run-v1" and report.get("status") == "manufactured_reference_only", "unsupported source report schema/status")
    _same(report.get("code_sha256"), _SOURCE_CODE_SHA256, "qualified source numerical-code version vector")
    _require(report.get("parameter_status") == "manufactured" and config.get("schema") == "mhdlab-thermal-width-power-manufactured-config-v1" and config.get("parameter_status") == "manufactured", "source configuration is not qualified manufactured data")
    _require(assembly.get("schema") == "mhdlab-thermal-width-power-reference-v1" and assembly.get("model_identity") == "manufactured-three-patch-matched-input-width-assembly-v1", "unsupported saved assembly identity")
    _require(assembly.get("parameter_status") == "manufactured" and assembly.get("host_budget_accepted") is False and assembly.get("owner_implementation_verified") is False, "unsupported physical/owner source acceptance")
    _label(assembly.get("parameter_provenance"), "source provenance")
    _require(independent.get("schema") == "mitl-independent-width-power-completion-v1" and independent.get("status") == "passed" and independent.get("all_checks_passed") is True, "actual independent source acceptance is required")
    _require(reference.get("status") == "passed" and reference.get("all_checks_passed") is True, "independent source report not passed")
    _require(independent.get("report_sha256") == hashes["independent_report"] and independent.get("producer_report_sha256") == hashes["report"], "independent source/report hash association mismatch")
    _require(review.get("schema") == "mitl-width-power-static-oracle-review-v1" and review.get("status") == "accepted" and review.get("further_changes_requested") is False, "final scientific source review is required")
    _require(review.get("owner_report_sha256") == hashes["independent_report"] and review.get("completed_report_bound_without_numerical_replay") is True, "scientific review does not bind actual independent source report")
    _require(type(review.get("findings")) is list and all(type(item) is dict and item.get("status") == "repaired_before_first_reference" for item in review["findings"]), "unresolved scientific source findings")
    for key, count in (("new_scalar_cases", 3), ("unique_scalar_reference_samples", 195), ("unique_scalar_flux_references", 195),
                       ("expanded_patch_samples_checked", 390), ("aggregate_samples_checked", 130), ("difference_samples_checked", 65),
                       ("post_source_cache_records_verified", 195), ("independent_RHS_field_comparisons", 1560),
                       ("new_independent_oracle_runs", 1), ("old_case_or_oracle_runs", 0)):
        _strict_int(independent.get(key), count, key)
    _same(independent.get("producer_imports_or_calls"), [], "independent producer-use declaration")
    _require(independent.get("all_eight_RHS_fields_checked") is True and independent.get("exact_area_and_power_identity_verified") is True
             and independent.get("full_patch_mixtures_preserved") is True, "source cache/area/mixture acceptance missing")
    _require(independent.get("dense_Taylor_coefficients_or_off_knot_queries_retained") is False, "unsupported independent checkpoint scope")
    input_pins = reference.get("input_sha256")
    _require(type(input_pins) is dict and all(type(key) is str and type(value) is str for key, value in input_pins.items()), "independent input hash catalogue required")
    _require(all(hashes[key] in input_pins.values() for key in ("report", "assembly", "input_config")), "independent source input association mismatch")
    artifact_pins = report.get("artifact_sha256")
    _require(type(artifact_pins) is dict and report.get("assembly_file") == "width_assembly.json", "source artifact catalogue required")
    names = {"assembly": "width_assembly.json", "input_config": "input_config.json", "reuse_map": "reuse_map.json",
             **{key: key+".json" for key in _POWERS}}
    for key, filename in names.items():
        _require(artifact_pins.get(filename) == hashes[key], f"source file/report association mismatch: {key}")
    _require(report.get("assembly_sha256") == hashes["assembly"] and report.get("config_sha256") == hashes["input_config"], "source report/config/assembly mismatch")
    for key, count in (("distinct_scalar_evolution_calls", 3), ("distinct_scalar_history_count", 3), ("unique_scalar_samples", 195)):
        _strict_int(report.get(key), count, key)
    _same(report.get("reuse_map"), _ASSIGNMENTS, "source report reuse map")
    _same(reuse.get("hypotheses"), _ASSIGNMENTS, "saved reuse map")
    _same(config.get("hypotheses"), _ASSIGNMENTS, "source configuration reuse map")
    _require(reuse.get("no_re_evolution_for_reused_patches") is True, "source reuse qualification missing")
    cache = assembly.get("cache")
    _require(type(cache) is dict, "source cache catalogue required")
    for key, count in (("distinct_scalar_histories", 3), ("scalar_samples", 195), ("rhs_evaluations", 195), ("patch_history_references", 6), ("scalar_evolutions_by_assembly", 0)):
        _strict_int(cache.get(key), count, key)
    _same(reuse.get("native_cache_diagnostics"), cache, "saved cache reuse association")

    histories = {key: docs[key] for key in _POWERS}
    source_records = assembly.get("scalar_histories")
    _require(type(source_records) is dict and set(source_records) == set(_POWERS), "three source history records required")
    rows = report.get("scalar_histories")
    _require(type(rows) is list and len(rows) == 3 and all(type(row) is dict for row in rows), "three source report history associations required")
    _require({row.get("case_id") for row in rows} == set(_POWERS), "source report case identities mismatch")
    report_rows = {row["case_id"]: row for row in rows}
    common_problem = common_grid = None
    duration = None
    for key, power in _POWERS.items():
        history, record, row = histories[key], source_records[key], report_rows[key]
        _require(history.get("schema") == "mhdlab-thermal-forced-feedback-reference-v1" and history.get("model_identity") == "manufactured-positive-barrier-piecewise-power-joint-rk4-v1", "unsupported scalar history model")
        canonical = _json_sha(history)
        _require(type(record) is dict and record.get("history_sha256") == canonical and row.get("canonical_history_sha256") == canonical, "canonical history association mismatch")
        _require(row.get("history_file") == key+".json" and row.get("history_file_sha256") == hashes[key], "raw scalar file association mismatch")
        _require(reuse.get("history_sha256", {}).get(key) == canonical, "source reuse canonical association mismatch")
        _same(record.get("history"), history, "inline scalar history")
        problem, schedule = history.get("problem"), history.get("schedule")
        _require(type(problem) is dict and type(schedule) is dict and problem.get("parameter_status") == "manufactured" and schedule.get("parameter_status") == "manufactured", "explicit manufactured scalar data required")
        _label(problem.get("parameter_provenance"), "problem provenance")
        _label(problem.get("energy_reference"), "energy reference")
        _require(type(problem.get("patch_area_m2")) is float and problem["patch_area_m2"] == 1.0, "unit reference area required")
        pool = problem.get("reservoir")
        _require(type(pool) is dict, "explicit reservoir data required")
        for name in ("name", "state_id", "kinetic_provenance", "energy_provenance"):
            _label(pool.get(name), name)
        for name in ("attempt_frequency_s", "reference_coverage_m2", "max_coverage_m2", "activation_energy_ev"):
            _require(_source_float(pool.get(name), name, nonnegative=True) > 0, f"positive {name} required")
        for name in ("desorption_energy_ev", "internal_energy_ev", "coverage_coefficient"):
            _source_float(pool.get(name), name, nonnegative=True)
        _require(_exact(pool["coverage_coefficient"], "alpha")*_exact(pool["max_coverage_m2"], "Nmax") <= _exact(pool["reference_coverage_m2"], "Nref"), "negative effective barrier domain")
        for name in ("areal_heat_capacity_j_m2_k", "initial_temperature_k", "temperature_min_k", "temperature_max_k"):
            _require(_source_float(problem.get(name), name, nonnegative=True) > 0, f"positive {name} required")
        _require(problem["temperature_min_k"] < problem["temperature_max_k"], "source temperature support must have positive width")
        n0 = _source_float(problem.get("initial_coverage_m2"), "initial coverage", nonnegative=True)
        _require(n0 <= _exact(pool["max_coverage_m2"], "Nmax"), "initial inventory outside source support")
        if common_problem is None:
            common_problem = problem
        else:
            _same(problem, common_problem, "common problem/state/cost/energy reference/provenance")
        _parameter_same(config.get("common_scalar_problem"), problem, "configuration scalar problem")
        segments = schedule.get("segments")
        _require(type(segments) is list and len(segments) == 1 and type(segments[0]) is dict, "constant-power one-segment source required")
        segment = segments[0]
        _strict_int(segment.get("steps"), 64, "source steps")
        _strict_int(row.get("steps"), 64, "source report steps")
        _strict_int(row.get("sample_count"), 65, "source report samples")
        _require(_source_float(segment.get("start_s"), "segment start") == 0, "source must start at canonical zero")
        local_duration = _source_float(schedule.get("duration_s"), "source duration", nonnegative=True)
        _require(local_duration > 0 and _source_float(segment.get("stop_s"), "segment stop") == local_duration, "source support mismatch")
        _require(_exact(segment.get("power_w_m2"), "source power") == _exact(power, "role power") and _exact(row.get("power_w_m2"), "reported power") == _exact(power, "role power"), "constant power/source identity mismatch")
        _label(segment.get("power_provenance"), "power provenance")
        _label(schedule.get("parameter_provenance"), "schedule provenance")
        if duration is None:
            duration = local_duration
        _require(local_duration == duration, "source durations must agree exactly")
        grid = record.get("logical_time_s")
        samples, rhs, fluxes = history.get("samples"), record.get("rhs"), record.get("release_flux_m2_s")
        _require(all(type(items) is list and len(items) == 65 for items in (grid, samples, rhs, fluxes)), "complete65-node source/cache grid required")
        expected_grid = [duration*index/64 for index in range(65)]
        _require([_from_ratio(item) for item in grid] == expected_grid, "exact logical grid mismatch")
        if common_grid is None:
            common_grid = grid
        _same(grid, common_grid, "common logical grid")
        previous = None
        for index, (sample, rates, flux, time) in enumerate(zip(samples, rhs, fluxes, expected_grid)):
            _require(type(sample) is dict and type(rates) is dict and set(rates) == _RHS_FIELDS, "complete scalar/RHS record required")
            _strict_int(sample.get("segment_index"), 0, "source segment index")
            _strict_int(sample.get("segment_step"), index, "source segment step")
            _same(sample.get("time_s"), _project(time, "source time"), "represented source timestamp")
            for name, value in sample.items():
                if name not in ("segment_index", "segment_step", "extensive"):
                    _source_float(value, name)
            gross_names = ("released_m2", "input_energy_j_m2", "binding_energy_j_m2", "internal_energy_j_m2",
                           "translational_energy_j_m2", "total_wall_energy_debit_j_m2")
            for name in (*gross_names, "gas_energy_j_m2", "adsorbate_reference_change_j_m2", "adsorbate_reference_inventory_change_j_m2"):
                _source_float(sample.get(name), name, nonnegative=True)
            n = _exact(sample["coverage_m2"], "coverage")
            m = _exact(sample["released_m2"], "released")
            temp = _exact(sample["temperature_k"], "temperature")
            _require(0 <= n <= _exact(pool["max_coverage_m2"], "Nmax") and 0 <= m <= n0, "source inventory support violation")
            _require(_exact(problem["temperature_min_k"], "Tmin") <= temp <= _exact(problem["temperature_max_k"], "Tmax"), "source temperature support violation")
            if index == 0:
                _same(sample["coverage_m2"], problem["initial_coverage_m2"], "initial coverage")
                _same(sample["temperature_k"], problem["initial_temperature_k"], "initial temperature")
                _require(m == 0, "source initial emitted count must vanish")
                _require(all(sample[name] == 0 for name in gross_names), "source initial gross channels must vanish")
            if previous is not None:
                _require(n <= _exact(previous["coverage_m2"], "previous N") and m >= _exact(previous["released_m2"], "previous M"), "source inventory monotonicity violation")
                _require(all(sample[name] >= previous[name] for name in gross_names), "source cumulative gross-channel monotonicity violation")
            for name, value in rates.items():
                _source_float(value, name, nonnegative=name not in ("coverage_rate_m2_s", "temperature_rate_k_s"))
            gamma = _source_float(flux, "saved flux", nonnegative=True)
            _require((n > 0) == (gamma > 0), "cached flux loses positive source support")
            _same(rates["release_flux_m2_s"], flux, "cached flux/RHS association")
            _require(_exact(rates["coverage_rate_m2_s"], "coverage rate", signed=True) == -gamma, "cached coverage-rate sign mismatch")
            _require(_exact(rates["input_power_j_m2_s"], "input power") == _exact(power, "role power"), "cached input power mismatch")
            extensive = sample.get("extensive")
            expected_extensive = {_source_extensive_key(name): value for name, value in sample.items() if _source_extensive_key(name) is not None}
            _same(extensive, expected_extensive, "unit-area extensive/areal association")
            previous = sample

    geometry = assembly.get("geometry")
    _require(type(geometry) is dict and geometry.get("parameter_status") == "manufactured", "manufactured exact geometry required")
    edges = [Fraction(-1, 512), Fraction(-1, 1024), Fraction(1, 1024), Fraction(1, 512)]
    length = Fraction(1, 128)
    areas = [(right-left)*length for left, right in zip(edges, edges[1:])]
    _require([_from_ratio(item) for item in geometry.get("width_edges_exact_m", [])] == edges and _from_ratio(geometry.get("common_length_exact_m")) == length, "unsupported exact source geometry")
    _same(geometry.get("width_edges_m"), [float(item) for item in edges], "represented geometry edges")
    _same(geometry.get("common_length_m"), float(length), "represented length")
    _require([_from_ratio(item) for item in geometry.get("patch_areas_exact_m2", [])] == areas and _from_ratio(geometry.get("total_area_exact_m2")) == sum(areas), "source area association mismatch")
    _same(geometry.get("patch_areas_m2"), [float(item) for item in areas], "represented patch areas")
    _same(geometry.get("total_area_m2"), float(sum(areas)), "represented total area")
    _require([_from_ratio(item) for item in geometry.get("area_weights", [])] == [Fraction(1, 4), Fraction(1, 2), Fraction(1, 4)], "source area weights mismatch")
    hypothesis_docs = assembly.get("hypotheses")
    _require(type(hypothesis_docs) is dict and set(hypothesis_docs) == set(_ASSIGNMENTS), "two source hypotheses required")
    for name, assignments in _ASSIGNMENTS.items():
        patches = hypothesis_docs[name].get("patches")
        _require(type(patches) is list and len(patches) == 3, "three source patch references required")
        for index, (patch, key, area) in enumerate(zip(patches, assignments, areas)):
            _strict_int(patch.get("patch_index"), index, "patch index")
            _require(patch.get("scalar_history_id") == key and patch.get("state_id") == common_problem["reservoir"]["state_id"], "source patch/state assignment mismatch")
            _require(patch.get("history_sha256") == source_records[key]["history_sha256"] and _from_ratio(patch.get("area_exact_m2")) == area, "source patch history/area mismatch")
            _same(patch.get("width_left_m"), float(edges[index]), "patch left edge")
            _same(patch.get("width_right_m"), float(edges[index+1]), "patch right edge")
            _same(patch.get("area_m2"), float(area), "patch represented area")
            patch_samples = patch.get("samples")
            _require(type(patch_samples) is list and len(patch_samples) == 65, "source patch complete time mixture required")
            for node, source_sample in enumerate(histories[key]["samples"]):
                saved = patch_samples[node]
                _strict_int(saved.get("segment_index"), 0, "patch segment index")
                _strict_int(saved.get("segment_step"), node, "patch sample index")
                _same(saved.get("time_s"), source_sample["time_s"], "patch time association")
                _same(saved.get("logical_time_s"), common_grid[node], "patch logical time")
                _require(saved.get("state_id") == patch["state_id"] and saved.get("scalar_history_id") == key, "patch node source/state mismatch")
                areal = {field: value for field, value in source_sample.items() if field not in ("segment_index", "segment_step", "time_s", "extensive")}
                areal["release_flux_m2_s"] = source_records[key]["release_flux_m2_s"][node]
                _same(saved.get("areal"), areal, "patch full areal-mixture association")
    _require(sum(area*_exact(_POWERS[key], "power") for area, key in zip(areas, _ASSIGNMENTS["uniform"])) ==
             sum(area*_exact(_POWERS[key], "power") for area, key in zip(areas, _ASSIGNMENTS["nonuniform"])), "exact prescribed source input does not match")
    admission = {"schema": "mhdlab-coupled-emission-saved-source-admission-v1", "status": "accepted_saved_manufactured_source",
                 "parameter_provenance": parameter_provenance, "file_sha256": hashes,
                 "canonical_history_sha256": {key: source_records[key]["history_sha256"] for key in _POWERS},
                 "assembly_canonical_sha256": _json_sha(assembly),
                 "numerical_code_sha256": report.get("code_sha256"), "accepted_independent_scope": {
                     key: independent[key] for key in ("new_scalar_cases", "unique_scalar_reference_samples", "unique_scalar_flux_references", "all_eight_RHS_fields_checked", "dense_Taylor_coefficients_or_off_knot_queries_retained")},
                 "source_calls": {"history_evolution": 0, "RHS": 0, "assembly": 0},
                 "limits": ["Trusted expected hashes are explicit caller inputs; no self-authenticating provenance claim.",
                            "Prior accepted source/cache/geometry associations are consumed as data, without recomputing Arrhenius or old assembly totals.",
                            "Saved RHS powers were projected independently; they are not reconstructed from rounded Gamma.",
                            "Prior acceptance supplies knot reference evidence only; new window/interpolation/reference error is separate."],
                 "host_budget_accepted": False, "owner_implementation_verified": False}
    return AdmittedCoupledSource(assembly, histories, admission)


def _extensive_key(key: str) -> str:
    if key in ("released_m2", "coverage_m2", "inventory_residual_m2"):
        return key[:-3]+"_number"
    _require(key.endswith("_m2"), "extensive conversion requires an areal quantity")
    return key[:-3]


def _record_value(record: Mapping[str, Any]) -> Fraction:
    return _from_ratio(record["exact"])


def _unavailable_reference() -> dict[str, Any]:
    return {"available": False, "reason": "new_window_IVP_query_error_is_verified_separately", "record": None}


def _source_comparison(history: dict[str, Any], grid: list[Fraction], result: dict[str, Any]) -> dict[str, Any]:
    start = _from_ratio(result["window_start_exact_s"])
    stop = _from_ratio(result["window_stop_exact_s"])
    moment_unavailable = {key: {"available": False, "reason": "source_does_not_store_cumulative_time_or_temperature_moments", "record": None}
                          for key in _RAW_MOMENTS[1:]}
    comparison = {"available": False, "reason": "window_bounds_not_both_saved_source_knots", "bounds_source_indices": None,
                  "stored_channel_deltas": {}, "polynomial_minus_source": {}, "source_residual_deltas": {},
                  "energy_defect_decomposition": {}, "wall_debit_difference_decomposition": {},
                  "reconstructed_source_total_energy_residual_j_m2": None, "energy_defect_difference_j_m2": None,
                  "moment_comparisons": moment_unavailable}
    if start not in grid or stop not in grid:
        return comparison
    ia, ib = grid.index(start), grid.index(stop)
    a, b = history["samples"][ia], history["samples"][ib]
    deltas = {key: _exact(b[key], key, signed=True)-_exact(a[key], key, signed=True)
              for key in a if key.endswith("_j_m2") or key in ("coverage_m2", "released_m2", "inventory_residual_m2")}
    polynomial = {key: _record_value(value) for key, value in result["areal"].items()}
    problem, pool = history["problem"], history["problem"]["reservoir"]
    d = _exact(pool["desorption_energy_ev"], "D")*_QE
    i = _exact(pool["internal_energy_ev"], "I")*_QE
    delta_n = _exact(a["coverage_m2"], "Na")-_exact(b["coverage_m2"], "Nb")
    thermal = _exact(problem["areal_heat_capacity_j_m2_k"], "capacity")*(
        _exact(b["temperature_k"], "Tb")-_exact(a["temperature_k"], "Ta"))
    source_residual = thermal+d*delta_n+deltas["internal_energy_j_m2"]+deltas["translational_energy_j_m2"]-deltas["input_energy_j_m2"]
    mass = polynomial["released_m2"]
    decomposition = {
        "source_inventory_removal_difference_j_m2": d*(mass-delta_n),
        "source_internal_allocation_difference_j_m2": i*mass-deltas["internal_energy_j_m2"],
        "translation_interpolant_minus_source_j_m2": polynomial["translational_energy_j_m2"]-deltas["translational_energy_j_m2"],
        "source_input_accumulation_difference_j_m2": deltas["input_energy_j_m2"]-polynomial["input_energy_j_m2"],
    }
    defect_difference = polynomial["total_energy_residual_j_m2"]-source_residual
    _require(sum(decomposition.values(), Fraction()) == defect_difference, "exact energy-defect decomposition failed")
    wall_decomposition = {
        "source_binding_allocation_difference_j_m2": d*mass-deltas["binding_energy_j_m2"],
        "source_internal_allocation_difference_j_m2": i*mass-deltas["internal_energy_j_m2"],
        "translation_interpolant_minus_source_j_m2": polynomial["translational_energy_j_m2"]-deltas["translational_energy_j_m2"],
        "negative_source_wall_channel_difference_j_m2": -(deltas["total_wall_energy_debit_j_m2"]-deltas["binding_energy_j_m2"]-deltas["internal_energy_j_m2"]-deltas["translational_energy_j_m2"]),
    }
    _require(sum(wall_decomposition.values(), Fraction()) == polynomial["total_wall_energy_debit_j_m2"]-deltas["total_wall_energy_debit_j_m2"], "exact wall-debit decomposition failed")
    comparison.update({
        "available": True, "reason": None, "bounds_source_indices": [ia, ib],
        "stored_channel_deltas": {key: _record(value, f"source delta {key}") for key, value in deltas.items()},
        "polynomial_minus_source": {key: _record(value-deltas[key], f"polynomial/source {key}") for key, value in polynomial.items() if key in deltas},
        "source_residual_deltas": {key: _record(value, f"source residual delta {key}") for key, value in deltas.items()
                                   if "residual" in key},
        "energy_defect_decomposition": {key: _record(value, key) for key, value in decomposition.items()},
        "wall_debit_difference_decomposition": {key: _record(value, key) for key, value in wall_decomposition.items()},
        "reconstructed_source_total_energy_residual_j_m2": _record(source_residual, "unprojected source-state energy residual"),
        "energy_defect_difference_j_m2": _record(defect_difference, "polynomial minus reconstructed source defect"),
    })
    return comparison


_COMPARISON_MAPS = ("stored_channel_deltas", "polynomial_minus_source", "source_residual_deltas",
                    "energy_defect_decomposition", "wall_debit_difference_decomposition")
_COMPARISON_SCALARS = ("reconstructed_source_total_energy_residual_j_m2", "energy_defect_difference_j_m2")


def _scale_comparison(comparison: dict[str, Any], area: Fraction, source_key: str) -> dict[str, Any]:
    result = {"available": comparison["available"], "reason": comparison["reason"], "source_window_key": source_key, "extensive": {}}
    if not comparison["available"]:
        return result
    for category in _COMPARISON_MAPS:
        result["extensive"][category] = {_extensive_key(key): _record(_record_value(record)*area, f"patch source comparison {key}")
                                            for key, record in comparison[category].items()}
    for key in _COMPARISON_SCALARS:
        result["extensive"][_extensive_key(key)] = _record(_record_value(comparison[key])*area, f"patch {key}")
    return result


def _sum_comparisons(comparisons: list[dict[str, Any]]) -> dict[str, Any]:
    result = {"available": all(item["available"] for item in comparisons), "reason": None,
              "source_window_keys": [item["source_window_key"] for item in comparisons], "extensive": {}}
    if not result["available"]:
        result["reason"] = "window_bounds_not_both_saved_source_knots"
        return result
    for category in _COMPARISON_MAPS:
        keys = comparisons[0]["extensive"][category]
        result["extensive"][category] = {key: _record(sum((_record_value(item["extensive"][category][key]) for item in comparisons), Fraction()), f"aggregate source {key}") for key in keys}
    for key in _COMPARISON_SCALARS:
        extensive_key = _extensive_key(key)
        result["extensive"][extensive_key] = _record(sum((_record_value(item["extensive"][extensive_key]) for item in comparisons), Fraction()), f"aggregate {key}")
    return result


def build_coupled_emission_windows(
    source: AdmittedCoupledSource, queries: Sequence[WindowQuery], *,
    parameter_status: str, parameter_provenance: str,
) -> CoupledEmissionWindows:
    """Construct nine unique windows and reuse their full joint piece mixtures."""
    _require(parameter_status == "manufactured", "window export requires parameter_status='manufactured'")
    _label(parameter_provenance, "export provenance")
    _require(type(source) is AdmittedCoupledSource, "explicitly admitted saved coupled source required")
    admission = source.to_dict()
    _require(admission.get("status") == "accepted_saved_manufactured_source", "saved source admission missing")
    _require(isinstance(queries, Sequence) and not isinstance(queries, (str, bytes)) and len(queries) == 3, "exactly three declared queries required")
    copied_queries = []
    for query in queries:
        _require(type(query) is WindowQuery, "typed WindowQuery required")
        copied_queries.append(WindowQuery(**{field.name: getattr(query, field.name) for field in fields(WindowQuery)}))
    _require(len({query.window_id for query in copied_queries}) == 3, "query identities must be unique")
    _require([(query.start_fraction, query.stop_fraction) for query in copied_queries] == [
        (Fraction(0), Fraction(1)), (Fraction(1, 16), Fraction(5, 16)), (Fraction(23, 128), Fraction(71, 128))], "unsupported bounded query set or order")
    assembly, histories = _thaw(source.assembly), _thaw(source.histories)
    _require(_json_sha(assembly) == admission.get("assembly_canonical_sha256"), "admitted saved assembly changed")
    for key in _POWERS:
        _require(_json_sha(histories[key]) == admission["canonical_history_sha256"][key], "admitted canonical history changed")
    first_history = histories["uniform_power"]
    duration = _exact(first_history["schedule"]["duration_s"], "duration")
    windows = []
    lookup: dict[tuple[str, str], dict[str, Any]] = {}
    for key in _POWERS:
        history = histories[key]
        problem = history["problem"]
        record = assembly["scalar_histories"][key]
        grid = [_from_ratio(item) for item in record["logical_time_s"]]
        cells = []
        for index, (start, stop) in enumerate(zip(grid, grid[1:])):
            a, b = history["samples"][index:index+2]
            h = stop-start
            gamma_a = _exact(record["rhs"][index]["release_flux_m2_s"], "Gamma0")
            gamma_b = _exact(record["rhs"][index+1]["release_flux_m2_s"], "Gamma1")
            tdot_a = _exact(record["rhs"][index]["temperature_rate_k_s"], "Tdot0", signed=True)
            tdot_b = _exact(record["rhs"][index+1]["temperature_rate_k_s"], "Tdot1", signed=True)
            ma, mb = _exact(a["released_m2"], "M0"), _exact(b["released_m2"], "M1")
            ta, tb = _exact(a["temperature_k"], "T0"), _exact(b["temperature_k"], "T1")
            cells.append(CoupledEmissionCell(start, stop, (ma, ma+h*gamma_a/3, mb-h*gamma_b/3, mb),
                                            (ta, ta+h*tdot_a/3, tb-h*tdot_b/3, tb),
                                            problem["temperature_min_k"], problem["temperature_max_k"], problem["initial_coverage_m2"],
                                            problem["reservoir"]["state_id"], key, "manufactured", admission["parameter_provenance"]))
        for query in copied_queries:
            result = integrate_emission_window(cells, start_s=duration*query.start_fraction, stop_s=duration*query.stop_fraction,
                                               desorption_energy_ev=problem["reservoir"]["desorption_energy_ev"],
                                               internal_energy_ev=problem["reservoir"]["internal_energy_ev"],
                                               areal_heat_capacity_j_m2_k=problem["areal_heat_capacity_j_m2_k"],
                                               power_w_m2=_POWERS[key], query_id=query.window_id, parameter_provenance=query.parameter_provenance)
            result.update({"source_window_key": key+"/"+query.window_id,
                           "canonical_history_sha256": admission["canonical_history_sha256"][key],
                           "raw_history_file_sha256": admission["file_sha256"][key], "constant_power_w_m2": _POWERS[key],
                           "source_RK_comparison": _source_comparison(history, grid, result), "independent_window_reference": _unavailable_reference()})
            for piece in result["pieces"]:
                index = piece["source_cell_index"]
                piece["source_endpoint_references"] = [
                    {"sample_index": node, "coverage_m2": history["samples"][node]["coverage_m2"],
                     "released_m2": history["samples"][node]["released_m2"],
                     "release_flux_m2_s": record["rhs"][node]["release_flux_m2_s"],
                     "temperature_rate_k_s": record["rhs"][node]["temperature_rate_k_s"],
                     "inventory_residual_m2": history["samples"][node]["inventory_residual_m2"],
                     "scalar_residuals": {name: value for name, value in history["samples"][node].items() if name.endswith("_residual_j_m2")}}
                    for node in (index, index+1)]
            windows.append(result)
            lookup[key, query.window_id] = result
    patches = []
    aggregates = []
    aggregate_lookup = {}
    areas = [_from_ratio(item) for item in assembly["geometry"]["patch_areas_exact_m2"]]
    for hypothesis, assignments in _ASSIGNMENTS.items():
        for query in copied_queries:
            accumulated = {}
            patch_indices = []
            comparisons = []
            for index, (key, area) in enumerate(zip(assignments, areas)):
                result = lookup[key, query.window_id]
                extensive = {_extensive_key(name): _record_value(value)*area for name, value in result["areal"].items()}
                for name, value in extensive.items():
                    accumulated[name] = _bounded(accumulated.get(name, Fraction())+value)
                comparison = _scale_comparison(result["source_RK_comparison"], area, result["source_window_key"])
                comparisons.append(comparison)
                patch_indices.append(len(patches))
                patches.append({"hypothesis": hypothesis, "patch_index": index, "query_id": query.window_id,
                                "source_history_id": key, "state_id": result["state_id"], "area_exact_m2": _ratio(area),
                                "source_window_key": result["source_window_key"],
                                "extensive": {name: _record(value, f"patch {name}") for name, value in extensive.items()},
                                "statistics": _thaw(_freeze(result["statistics"])), "source_RK_comparison": comparison})
            representative = lookup[assignments[0], query.window_id]
            raw_moments = {name: accumulated[_extensive_key(name)] for name in _RAW_MOMENTS}
            aggregate = {"hypothesis": hypothesis, "query_id": query.window_id, "state_id": representative["state_id"],
                         "window_start_exact_s": representative["window_start_exact_s"], "window_stop_exact_s": representative["window_stop_exact_s"],
                         "local_time_anchor_exact_s": representative["local_time_anchor_exact_s"], "patch_references": patch_indices,
                         "extensive": {name: _record(value, f"aggregate {name}") for name, value in accumulated.items()},
                         "statistics": _statistics(raw_moments), "source_RK_comparison": _sum_comparisons(comparisons)}
            aggregates.append(aggregate)
            aggregate_lookup[hypothesis, query.window_id] = aggregate
    contrasts = []
    for query in copied_queries:
        uniform = aggregate_lookup["uniform", query.window_id]
        nonuniform = aggregate_lookup["nonuniform", query.window_id]
        statistics = {}
        for name in uniform["statistics"]:
            a, b = uniform["statistics"][name], nonuniform["statistics"][name]
            statistics[name] = ({"available": True, "reason": None, "record": _record(_record_value(b["record"])-_record_value(a["record"]), f"statistic difference {name}")}
                                if a["available"] and b["available"] else {"available": False, "reason": "one_or_both_emission_measures_have_zero_mass", "record": None})
        contrasts.append({"query_id": query.window_id,
                          "quantity_difference_nonuniform_minus_uniform": {name: _record(_record_value(nonuniform["extensive"][name])-_record_value(uniform["extensive"][name]), f"window contrast {name}") for name in uniform["extensive"]},
                          "statistics_difference_nonuniform_minus_uniform": statistics,
                          "independent_window_reference": _unavailable_reference()})
    areal_units = {"released_m2": "parents m^-2", "local_time_first_moment_s_m2": "parent s m^-2",
                   "local_time_second_moment_s2_m2": "parent s^2 m^-2", "temperature_first_moment_k_m2": "parent K m^-2",
                   "temperature_second_moment_k2_m2": "parent K^2 m^-2", "local_time_temperature_moment_k_s_m2": "parent K s m^-2"}
    areal_units.update({key: "J m^-2" for key in windows[0]["areal"] if key.endswith("_j_m2")})
    document = {"schema": "mhdlab-coupled-emission-windows-v1", "model_identity": MODEL_IDENTITY,
                "parameter_status": parameter_status, "parameter_provenance": parameter_provenance,
                "source_admission": admission,
                "queries": [{"window_id": query.window_id, "start_fraction": _ratio(query.start_fraction), "stop_fraction": _ratio(query.stop_fraction), "parameter_provenance": query.parameter_provenance} for query in copied_queries],
                "source_windows": windows, "patch_windows": patches, "aggregate_windows": aggregates, "width_window_contrasts": contrasts,
                "cache": {"distinct_source_windows": 9, "patch_references": 18, "aggregate_windows": 6, "contrasts": 3,
                          "source_calls": {"history_evolution": 0, "RHS": 0, "assembly": 0}},
                "quantity_units": {"areal": areal_units,
                                   "extensive": {_extensive_key(key): ("J" if key.endswith("_j_m2") else unit.replace(" m^-2", "")) for key, unit in areal_units.items()},
                                   "statistics": {"mean_local_time_s": "s", "variance_local_time_s2": "s^2", "mean_temperature_k": "K", "variance_temperature_k2": "K^2", "mean_local_time_temperature_k_s": "K s", "covariance_local_time_temperature_k_s": "K s"}},
                "limits": ["Exact rational algebra qualifies the declared positive interpolant, not an exact continuous kinetic trajectory.",
                           "Only the three admitted manufactured constant-power sources and finite declared windows are supported; scheduled jumps need a separate one-sided derivative contract.",
                           "All temperatures, times, states, source pieces and area associations remain explicit; normalized moments are statistics, not a replacement thermal launch model.",
                           "Polynomial-versus-source channel discrepancies, source scalar residuals and new independent IVP/query error are separate layers; none is repaired by fitted temperature or normalization.",
                           "No source evolution/RHS/assembly is repeated. No packet sampler, physical velocity, species mass, molecular state data, CR, kinetic gap, FLASH, Hall, spectra or owner coupling is supplied."],
                "host_budget_accepted": False, "owner_implementation_verified": False}
    return CoupledEmissionWindows(document)
