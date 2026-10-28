"""Strict algebraic admission of full-component source-window ledger JSON.

This module uses only standard-library JSON/type and Fraction arithmetic. It
does not import or evaluate a source law, history, ledger producer, moment rule
or sampler. Component quantity qualifications are declared input evidence.
Optional binding compares the complete declarations with a separately supplied
conditioned-source manifest; neither operation proves kinetic or physical truth
or accepts a host transaction.

Reduced coupling caches are separately rounded scalar declarations. An upward
log-to-exp round trip is admitted only when it exactly matches the producer's
declared scalar transform, rather than an arbitrary closeness envelope. That
check uses the current platform's math library and is not a law inversion.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
import re
import sys


SCHEMA = "mhdlab-window-source-ledger-admission-v1"
LEDGER_SCHEMA = "mhdlab-window-source-ledger-v1"
VALIDATION_IDENTITY = "strict_json_full_component_exact_fraction_ledger_algebra_v1"
MAX_COMPONENTS = 200_000
MAX_GROUPS = 250_000
MAX_JSON_NODES = 5_000_000
MAX_JSON_DEPTH = 64
MAX_STRING_LENGTH = 1_000_000
MAX_TOTAL_STRING_LENGTH = 128_000_000
MAX_RATIONAL_DIGITS = 4096
QUANTITIES = ("released_number", "translational_energy_j", "internal_energy_j",
              "binding_energy_j", "total_wall_energy_debit_j")
_KB = Fraction.from_float(1.380649e-23)
_ZERO = Fraction(0)
_ONE = Fraction(1)
_INTEGER = re.compile(r"(?:0|-[1-9][0-9]*|[1-9][0-9]*)\Z")
_LAW_IDENTITY = "temkin_first_order_constant_attempt_frequency_v1"
_POOL_TEXT = ("name", "state_id", "kinetic_provenance", "energy_provenance")
_POOL_POSITIVE = ("attempt_frequency_s", "reference_coverage_m2", "max_coverage_m2")
_POOL_NONNEGATIVE = ("activation_energy_ev", "desorption_energy_ev", "internal_energy_ev", "coverage_coefficient")
_ASSOCIATION_KEYS = {"component_index", "selection_index", "interval_index", "reservoir_index", "reservoir_name",
                     "state_id", "current_index", "width_index", "interval_origin_s", "x_bounds_m", "z_bounds_m"}
_ORIGINAL_COMPONENT_KEYS = _ASSOCIATION_KEYS-{"component_index", "reservoir_name", "state_id"}
_SOURCE_QUANTITIES = ("initial_coverage_m2", "initial_number", "remaining_start_m2", "remaining_start_number",
                      "remaining_stop_m2", "remaining_stop_number", "released_m2", "released_number",
                      "translational_energy_j_m2", "translational_energy_j", "internal_energy_j_m2", "internal_energy_j",
                      "binding_energy_j_m2", "binding_energy_j", "total_wall_energy_debit_j_m2", "total_wall_energy_debit_j")
_REPRESENTATION_KEYS = {"value", "log_value", "exact_zero", "qualified", "reason", "log_rounding_error", "log_rounding_allowance"}
_WINDOW_SCALARS = ("initial_coverage_m2", "temperature_k", "parent_duration_s", "start_s", "stop_s", "patch_area_m2", "log_mass_tolerance")
_WINDOW_KEYS = {"schema", "status", "provenance", "reservoir", "law_identity", *_WINDOW_SCALARS,
                "time_window", "quantities", "quantity_units", "all_budget_values_qualified", "inventory_balance_residual_m2",
                "inventory_balance_residual_number", "energy_balance_residual_j_m2", "energy_balance_residual_j",
                "host_budget_accepted", "numerical_allowance_allocation", "limits"}
_TIME_SCALARS = ("coupling", "log_initial_rate_s", "parent_duration_s", "start_s", "stop_s", "local_duration_s",
                 "start_depletion", "local_depletion", "parent_depletion", "log_start_depletion", "log_local_depletion",
                 "log_parent_depletion", "log_local_coupling", "local_coupling", "log_local_rate_s", "log_source_fraction",
                 "log_parent_conditional_probability", "log_mass_tolerance")
_SOURCE_LIMITS = [
    "Partial source emission in one stationary frozen-temperature patch; no photons, transport or chemistry are evolved.",
    "Density, actual patch number and diagnostic energy each retain their own representation flags; a finite log is not an accepted budget.",
    "Initial independent pool/state inventory remains attached to the supplied reservoir; no new parent, state branch or daughter yield is created.",
    "Signed inventory rounding residuals are retained and require the owner ledger's explicit treatment; no automatic residual repair or transaction is performed.",
    "Translation is expected outgoing flux energy 2*kB*T. The supplied internal/removal costs are separate from the kinetic barrier; no wall cooling feedback or individual velocity sample is calculated.",
    "Source start and local delay remain separate from any global timestamp or recorded photon gate.",
]
_LEDGER_LIMITS = [
    "Unqualified positive contributor amounts remain unavailable; qualified subtotals are not complete totals.",
    "Exact-zero source costs and empty declared bins are distinct from unrepresented positive amounts.",
    "SI projections can underflow, overflow or round coarsely; exact rational amounts and signed projection errors remain authoritative diagnostics, not accepted budgets.",
    "Density divides the unique exact geometric area once; canonical float patch areas used by source quantities are retained through signed area residuals.",
    "Matched uniform width changes only width position; current/state/time/temperature and supplied independent energy-cost mixtures remain those of the original components.",
    "Translation is expected outgoing flux energy 2*kB*T, not sampled particle energy; binding/removal cost is independent of the kinetic barrier.",
    "This ledger accepts no host transaction, changes no source history, and evolves no gap, molecular chemistry, radiation, photon gate or thermal feedback.",
]
_CONDITIONED_SOURCE_LIMITS = [
    "Positive original initial inventories are retained regardless of individual partial-quantity representation flags.",
    "Exact-zero canonical initial densities are skipped explicitly; no cross-interval mathematical survival is reconstructed from a rounded zero.",
    "No source history is evolved, no inventory is changed and no particle or energy transaction is accepted.",
    "Status and provenance are caller declarations; original surface attributes, pool records and parameter provenance remain separate.",
    "Source-emission restrictions are not photon observation gates; no transport, collisions, chemistry or spectra are evaluated.",
]


class LedgerAdmissionError(ValueError):
    """A JSON declaration fails strict structural, algebraic or binding checks."""


def _fail(path, message):
    raise LedgerAdmissionError(f"{path}: {message}")


def _keys(value, expected, path):
    if type(value) is not dict or set(value) != set(expected):
        _fail(path, "object keys differ from the supported schema")


def _same(actual, expected, path):
    if type(actual) is not type(expected):
        _fail(path, "JSON type differs from the required declaration")
    if type(expected) is float:
        if actual.hex() != expected.hex():
            _fail(path, "float representation differs (including signed zero)")
    elif type(expected) is dict:
        _keys(actual, expected, path)
        for key, value in expected.items():
            _same(actual[key], value, path+"/"+key)
    elif type(expected) is list:
        if len(actual) != len(expected):
            _fail(path, "list length differs")
        for index, (left, right) in enumerate(zip(actual, expected)):
            _same(left, right, path+f"/{index}")
    elif actual != expected:
        _fail(path, "value differs from the required declaration or exact algebra")


def _canonical_fraction(value, path):
    _keys(value, {"numerator", "denominator"}, path)
    for name in ("numerator", "denominator"):
        text = value[name]
        if type(text) is not str or len(text.lstrip("-")) > MAX_RATIONAL_DIGITS:
            _fail(path, "rational integer string exceeds its digit/type limit")
        if not _INTEGER.fullmatch(text):
            _fail(path, "rational integer strings must be canonical decimal integers")
    numerator, denominator = int(value["numerator"]), int(value["denominator"])
    if denominator <= 0 or math.gcd(numerator, denominator) != 1:
        _fail(path, "rational must be reduced with positive denominator")
    return Fraction(numerator, denominator)


def _json_domain(document, path):
    stack, nodes, text_size = [(document, 0, path)], 0, 0
    while stack:
        value, depth, current = stack.pop()
        nodes += 1
        if nodes > MAX_JSON_NODES or depth > MAX_JSON_DEPTH:
            _fail(current, "JSON node/depth work limit exceeded")
        kind = type(value)
        if value is None or kind is bool:
            continue
        if kind is float:
            if not math.isfinite(value):
                _fail(current, "nonfinite JSON floats are forbidden")
        elif kind is int:
            if value.bit_length() > 4096:
                _fail(current, "JSON integer work limit exceeded")
        elif kind is str:
            text_size += len(value)
            if len(value) > MAX_STRING_LENGTH or text_size > MAX_TOTAL_STRING_LENGTH:
                _fail(current, "JSON string work limit exceeded")
        elif kind is list:
            if len(value) > MAX_JSON_NODES-nodes:
                _fail(current, "JSON collection work limit exceeded")
            stack.extend((item, depth+1, current+f"/{index}") for index, item in enumerate(value))
        elif kind is dict:
            if any(type(key) is not str for key in value):
                _fail(current, "JSON object keys must be plain strings")
            if set(value) == {"numerator", "denominator"}:
                _canonical_fraction(value, current)
            for key, item in value.items():
                text_size += len(key)
                if len(key) > MAX_STRING_LENGTH or text_size > MAX_TOTAL_STRING_LENGTH:
                    _fail(current, "JSON object-key string work limit exceeded")
                stack.append((item, depth+1, current+"/"+key))
        else:
            _fail(current, "only plain strict JSON types are admitted")


def _float(value, path, *, positive=False, nonnegative=False):
    if type(value) is not float or not math.isfinite(value):
        _fail(path, "a finite JSON float is required")
    if positive and value <= 0 or nonnegative and value < 0:
        _fail(path, "float violates its positive/nonnegative domain")
    return value


def _index(value, limit, path):
    if type(value) is not int or not 0 <= value < limit:
        _fail(path, "integer index is outside its declared support")
    return value


def _text(value, path):
    if type(value) is not str or not value.strip():
        _fail(path, "nonempty explicit string required")


def _edges(value, path, *, time=False):
    if type(value) is not list or not 2 <= len(value) <= MAX_COMPONENTS+1:
        _fail(path, "edge list is empty or exceeds its work limit")
    for item in value:
        _float(item, path)
    if (time and value[0] < 0) or any(right <= left for left, right in zip(value[:-1], value[1:])):
        _fail(path, "edges must be strictly increasing with nonnegative source time")
    return value


def _pool(value, path):
    _keys(value, {*_POOL_TEXT, *_POOL_POSITIVE, *_POOL_NONNEGATIVE}, path)
    for name in _POOL_TEXT:
        _text(value[name], path+"/"+name)
    for name in _POOL_POSITIVE:
        _float(value[name], path+"/"+name, positive=True)
    for name in _POOL_NONNEGATIVE:
        _float(value[name], path+"/"+name, nonnegative=True)
    if value["activation_energy_ev"] > 0 and Fraction(value["coverage_coefficient"])*Fraction(value["max_coverage_m2"]) > Fraction(value["reference_coverage_m2"]):
        _fail(path, "declared reservoir coverage domain permits a negative barrier")


def _representation(record, zero, half_tolerance, path):
    _keys(record, _REPRESENTATION_KEYS, path)
    if type(record["qualified"]) is not bool or type(record["exact_zero"]) is not bool:
        _fail(path, "qualification/exact-zero flags must be booleans")
    _same(record["exact_zero"], zero, path+"/exact_zero")
    _float(record["log_rounding_allowance"], path+"/log_rounding_allowance", nonnegative=True)
    if zero:
        _same(record, {"value": 0., "log_value": None, "exact_zero": True, "qualified": True,
                       "reason": "exact_zero_supplied_energy_cost", "log_rounding_error": 0., "log_rounding_allowance": 0.}, path)
        return
    if record["log_value"] is not None:
        _float(record["log_value"], path+"/log_value")
    if record["value"] is not None:
        _float(record["value"], path+"/value", nonnegative=True)
    if record["log_rounding_error"] is not None:
        _float(record["log_rounding_error"], path+"/log_rounding_error", nonnegative=True)
    reason = record["reason"]
    if record["qualified"] or reason == "log_rounding_allowance_exceeded":
        if record["value"] is None or record["value"] <= 0 or record["log_value"] is None or record["log_rounding_error"] is None:
            _fail(path, "positive represented quantity lacks required value/log/rounding evidence")
        qualifies = record["log_rounding_error"]+record["log_rounding_allowance"] <= half_tolerance
        _same(record["qualified"], qualifies, path+"/qualified")
        _same(reason, "represented_within_declared_log_allowance" if qualifies else "log_rounding_allowance_exceeded", path+"/reason")
    elif reason == "positive_amount_underflow":
        _same(record["value"], 0., path+"/value")
        if record["log_value"] is None or record["log_rounding_error"] is not None:
            _fail(path, "positive underflow requires a log and no projected log-rounding error")
    elif reason in ("positive_amount_overflow", "log_amount_outside_float64_range"):
        if record["value"] is not None or record["log_rounding_error"] is not None:
            _fail(path, "overflow/out-of-range amount must not claim a represented value")
        if (record["log_value"] is None) != (reason == "log_amount_outside_float64_range"):
            _fail(path, "overflow log evidence disagrees with its declared reason")
    else:
        _fail(path, "unsupported positive representation reason")


def _time_window(value, physical, path):
    _keys(value, {"schema", *_TIME_SCALARS, "time_coordinate", "numerical_diagnostics"}, path)
    _same(value["schema"], "mhdlab-reduced-temkin-emission-window-v1", path+"/schema")
    nullable = {"log_start_depletion", "log_local_coupling"}
    for name in _TIME_SCALARS:
        if value[name] is not None or name not in nullable:
            _float(value[name], path+"/"+name)
    _float(value["log_mass_tolerance"], path+"/log_mass_tolerance", positive=True)
    for name in ("parent_duration_s", "start_s", "stop_s"):
        _same(value[name], physical[name], path+"/"+name)
    _same(value["local_duration_s"], float(physical["stop_s"]-physical["start_s"]), path+"/local_duration_s")
    _same(value["log_mass_tolerance"], physical["log_mass_tolerance"]/2, path+"/log_mass_tolerance")
    if (not 0 <= value["coupling"] <= 512. or value["local_coupling"] < 0
            or min(value["start_depletion"], value["local_depletion"], value["parent_depletion"]) < 0
            or value["parent_depletion"] < value["start_depletion"]
            or value["log_source_fraction"] > 0 or value["log_parent_conditional_probability"] > 0):
        _fail(path, "reduced-window declarations violate their structural domains")
    if value["local_coupling"] > value["coupling"]:
        # The mathematical c_a<=c does not imply the same exact ordering for
        # the independent binary64 caches c and exp(log(c)-s_a). Permit only
        # that documented producer round trip, with strict scalar identities.
        if value["coupling"] <= 0 or value["log_local_coupling"] is None:
            _fail(path, "upward coupling cache lacks its declared log projection")
        _same(value["log_local_coupling"], math.log(value["coupling"])-value["start_depletion"], path+"/log_local_coupling_roundtrip")
        _same(value["local_coupling"], math.exp(value["log_local_coupling"]), path+"/local_coupling_roundtrip")
    _same(value["time_coordinate"], "local delay from original source-window start; keep start and delay separately", path+"/time_coordinate")
    diagnostic = value["numerical_diagnostics"]
    diagnostic_keys = {"qualification_identity", "quadrature_assumption", "log_mass_tolerance", "estimated_log_source_fraction_error",
                       "estimated_log_parent_probability_error", "start_depletion_error", "local_log_depletion_error", "parent_log_depletion_error",
                       "local_coupling_approximated_as_zero", "log_local_coupling_approximation_error_bound", "local_coupling_approximation_error_bound",
                       "start_rebase_log_mass_sensitivity_bound", "local_duration_relative_rounding_error", "local_depletion_underflowed",
                       "start_depletion_underflowed", "root_evidence", "tail_switch_depletion", "certified_error_bound", "host_budget_accepted"}
    _keys(diagnostic, diagnostic_keys, path+"/numerical_diagnostics")
    _same(diagnostic["qualification_identity"], "estimated_monotonic_log_depletion_envelopes_with_positive_tail_correction_v1", path+"/qualification_identity")
    _same(diagnostic["quadrature_assumption"], "16 times SciPy absolute quadrature error estimate plus explicit float64 rounding allowances; estimated envelopes, not certified bounds", path+"/quadrature_assumption")
    _same(diagnostic["log_mass_tolerance"], value["log_mass_tolerance"], path+"/diagnostic_tolerance")
    for name in ("estimated_log_source_fraction_error", "estimated_log_parent_probability_error", "start_depletion_error",
                 "local_log_depletion_error", "parent_log_depletion_error", "local_coupling_approximation_error_bound",
                 "start_rebase_log_mass_sensitivity_bound", "local_duration_relative_rounding_error", "tail_switch_depletion"):
        _float(diagnostic[name], path+"/"+name, nonnegative=True)
    if max(diagnostic["estimated_log_source_fraction_error"], diagnostic["estimated_log_parent_probability_error"]) > value["log_mass_tolerance"]:
        _fail(path, "declared estimated errors exceed the declared log-mass tolerance")
    for name in ("certified_error_bound", "host_budget_accepted"):
        _same(diagnostic[name], False, path+"/"+name)
    for name, expected in (("local_depletion_underflowed", value["local_depletion"] == 0.),
                           ("start_depletion_underflowed", value["start_s"] > 0 and value["start_depletion"] == 0.)):
        _same(diagnostic[name], expected, path+"/"+name)
    if type(diagnostic["local_coupling_approximated_as_zero"]) is not bool:
        _fail(path, "coupling approximation flag must be boolean")
    approximation_log = diagnostic["log_local_coupling_approximation_error_bound"]
    if approximation_log is not None:
        _float(approximation_log, path+"/log_local_coupling_approximation_error_bound")
    if diagnostic["local_coupling_approximated_as_zero"] != (value["coupling"] > 0 and value["local_coupling"] == 0.):
        _fail(path, "coupling approximation declaration disagrees with retained values")
    _same(diagnostic["start_rebase_log_mass_sensitivity_bound"], 1.+value["coupling"], path+"/sensitivity")
    _same(diagnostic["tail_switch_depletion"], 32., path+"/tail_switch")
    roots = diagnostic["root_evidence"]
    _keys(roots, {"start", "local", "parent"}, path+"/root_evidence")
    for name in ("start", "local", "parent"):
        record = roots[name]
        _keys(record, {"log_depletion", "depletion", "log_error", "absolute_error", "iterations", "integral_method"}, path+"/root_evidence/"+name)
        for field in ("depletion", "log_error", "absolute_error"):
            _float(record[field], path+"/root_evidence/"+field, nonnegative=True)
        if record["log_depletion"] is not None:
            _float(record["log_depletion"], path+"/root_evidence/log_depletion")
        _index(record["iterations"], 161, path+"/root_evidence/iterations")
        _text(record["integral_method"], path+"/root_evidence/integral_method")
        _same(record["depletion"], value[name+"_depletion"], path+"/root_evidence/"+name+"/depletion")
        _same(record["log_depletion"], value["log_"+name+"_depletion"], path+"/root_evidence/"+name+"/log_depletion")


def _physical_window(value, pool, selection, origin, parent_duration, canonical_area, association, path):
    _keys(value, _WINDOW_KEYS, path)
    _same(value["schema"], "mhdlab-temkin-source-window-v1", path+"/schema")
    _same(value["law_identity"], _LAW_IDENTITY, path+"/law_identity")
    _same(value["reservoir"], pool, path+"/reservoir")
    for name in _WINDOW_SCALARS:
        _float(value[name], path+"/"+name, positive=name not in ("start_s",), nonnegative=True)
    if value["initial_coverage_m2"] > pool["max_coverage_m2"]:
        _fail(path, "initial inventory exceeds the declared reservoir domain")
    for name in ("status", "provenance", "log_mass_tolerance"):
        _same(value[name], association[name], path+"/"+name)
    _same(value["parent_duration_s"], parent_duration, path+"/parent_duration_s")
    _same(value["patch_area_m2"], canonical_area, path+"/patch_area_m2")
    _same(value["start_s"], selection["start_s"], path+"/start_s")
    _same(value["stop_s"], selection["stop_s"], path+"/stop_s")
    _same(value["host_budget_accepted"], False, path+"/host_budget_accepted")
    _same(value["limits"], _SOURCE_LIMITS, path+"/limits")
    _time_window(value["time_window"], value, path+"/time_window")
    records = value["quantities"]
    _keys(records, _SOURCE_QUANTITIES, path+"/quantities")
    for name in _SOURCE_QUANTITIES:
        zero = name.startswith("internal_energy_") and pool["internal_energy_ev"] == 0 or name.startswith("binding_energy_") and pool["desorption_energy_ev"] == 0
        _representation(records[name], zero, value["log_mass_tolerance"]/2, path+"/quantities/"+name)
    for name in ("initial_coverage_m2", "initial_number"):
        _same(records[name]["qualified"], True, path+"/quantities/"+name+"/qualified")
    _same(records["initial_coverage_m2"]["value"], value["initial_coverage_m2"], path+"/initial_coverage_value")
    _same(value["all_budget_values_qualified"], all(record["qualified"] for record in records.values()), path+"/all_budget_values_qualified")
    units = {name: "J m^-2" if name.endswith("_j_m2") else "J" if name.endswith("_j") else "parent" if name.endswith("_number") else "parent m^-2" for name in records}
    _same(value["quantity_units"], units, path+"/quantity_units")
    for suffix in ("m2", "number"):
        operands = [records[name+suffix]["value"] for name in ("remaining_stop_", "remaining_start_", "released_")]
        residual = None if any(item is None for item in operands) else math.fsum((operands[0], -operands[1], operands[2]))
        _same(value["inventory_balance_residual_"+suffix], residual, path+"/inventory_balance_residual_"+suffix)
    for suffix in ("j_m2", "j"):
        operands = [records[name+suffix]["value"] for name in ("total_wall_energy_debit_", "translational_energy_", "internal_energy_", "binding_energy_")]
        residual = None if any(item is None for item in operands) else math.fsum((operands[0], -operands[1], -operands[2], -operands[3]))
        _same(value["energy_balance_residual_"+suffix], residual, path+"/energy_balance_residual_"+suffix)
    _same(value["numerical_allowance_allocation"], {"time_law_log_mass": value["log_mass_tolerance"]/2,
          "representation_log_amount": value["log_mass_tolerance"]/2,
          "assumption": "Extended intermediates plus conservative float64 elementary-operation allowance; not a certified interval enclosure."}, path+"/numerical_allowance_allocation")


def _rat(value):
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _number(value):
    try:
        projected = float(value)
    except OverflowError:
        projected = None
    if projected is not None and not math.isfinite(projected):
        projected = None
    error = None if projected is None else Fraction(projected)-value
    return {"exact_sum": _rat(value), "value": projected, "projection": {
        "exact_zero": value == 0, "positive": value > 0, "negative": value < 0, "finite": projected is not None,
        "underflowed": value != 0 and projected == 0., "overflowed": projected is None,
        "rounded": error is None or error != 0,
        "subnormal": projected is not None and 0 < abs(projected) < sys.float_info.min,
        "exact_rounding_error": None if error is None else _rat(error)}}


class _Totals:
    def __init__(self):
        self.sums = {name: _ZERO for name in QUANTITIES}
        self.missing = {name: [] for name in QUANTITIES}
        self.thermal = _ZERO
        self.count = 0

    def add(self, values, missing, temperature):
        for name in QUANTITIES:
            if values[name] is None:
                self.missing[name].append(missing[name])
            else:
                self.sums[name] += values[name]
        if values["released_number"] is not None:
            self.thermal += values["released_number"]*2*_KB*Fraction(temperature)
        self.count += 1


def _amount(value, missing, area, factor, reference):
    value *= factor
    absent = {"exact_sum": None, "value": None, "projection": None}
    evidence = {"count": len(missing)}
    evidence["catalog_indices" if reference is None else "source_group_ref"] = list(missing) if reference is None else reference
    return {"available": not missing, **(_number(value) if not missing else absent),
            "qualified_subtotal": _number(value), "per_area": {"available": not missing,
            **(_number(value/area) if not missing else absent), "qualified_subtotal": _number(value/area)}, "unqualified": evidence}


def _group(totals, key, area, canonical, definition, *, factor=_ONE, reference=None):
    quantities = {name: _amount(totals.sums[name], totals.missing[name], area, factor,
                                None if reference is None else reference+"/quantities/"+name) for name in QUANTITIES}
    def residual(value, names):
        missing = [name for name in names if totals.missing[name]]
        return {"available": not missing, **(_number(value*factor) if not missing else {"exact_sum": None, "value": None, "projection": None}),
                "unavailable_quantities": missing}
    number = _number(area)
    result = {"key": key, "area": {"exact_m2": number["exact_sum"], "value_m2": number["value"], "projection": number["projection"],
                "canonical_patch_area_sum": _number(canonical), "canonical_minus_exact_residual_m2": _number(canonical-area), "definition": definition},
              "component_count": totals.count, "quantities": quantities, "residuals": {
                  "wall_minus_parts_j": residual(totals.sums[QUANTITIES[-1]]-sum((totals.sums[name] for name in QUANTITIES[1:-1]), _ZERO), QUANTITIES[1:]),
                  "translation_minus_count_2kbt_j": residual(totals.sums["translational_energy_j"]-totals.thermal, ("translational_energy_j", "released_number")),
                  "count_times_2kbt_j": _amount(totals.thermal, totals.missing["released_number"], area, factor,
                                                None if reference is None else reference+"/residuals/count_times_2kbt_j")}, "host_budget_accepted": False}
    if reference is not None:
        result["redistribution"] = {"source_group_ref": reference, "fraction": _rat(factor),
                                    "interpretation": "change width only; preserve original current/state/time/temperature/cost mixture"}
    return result


def _bind(document, expected):
    keys = {"schema", "selections", "components", "skipped_zero_initials", "source_reservoirs", "x_edges_m", "z_edges_m", "time_edges_s",
            "surface_attributes", "parameter_provenance", "status", "provenance", "log_mass_tolerance", "host_budget_accepted", "time_coordinate", "limits"}
    _keys(expected, keys, "expected_source")
    _same(expected["schema"], "mhdlab-temkin-planar-source-windows-v1", "expected_source/schema")
    _same(expected["host_budget_accepted"], False, "expected_source/host_budget_accepted")
    _same(expected["time_coordinate"], "interval_origin_s remains separate; selection start/stop are original interval offsets; emission quantiles are local delays from selection start", "expected_source/time_coordinate")
    _same(expected["limits"], _CONDITIONED_SOURCE_LIMITS, "expected_source/limits")
    for name in ("x_edges_m", "z_edges_m", "time_edges_s"):
        _same(document["geometry"][name], expected[name], "binding/"+name)
    for name, value in document["source_association"].items():
        _same(value, expected[name], "binding/"+name)
    originals = [{**{name: component[name] for name in _ORIGINAL_COMPONENT_KEYS}, "source_window": component["source_window"]}
                 for component in document["components"]]
    _same(originals, expected["components"], "binding/components")
    _same(document["skipped_zero_initials"], expected["skipped_zero_initials"], "binding/skipped_zero_initials")


def validate_window_source_ledger(document, *, expected_source=None):
    """Admit full-component JSON algebra, optionally bound to original declarations.

    No source equations or quantity logs are independently solved. The result
    distinguishes algebraic consistency, external declaration binding and host
    acceptance. Only the first two can be established here; host acceptance is
    always false. Component-free budget summaries are intentionally inadmissible.
    """
    _json_domain(document, "ledger")
    if expected_source is not None:
        _json_domain(expected_source, "expected_source")
    _keys(document, {"schema", "resolved", "matched_uniform_width", "uniform_width_bases", "unqualified_contributors",
                     "skipped_zero_initials", "source_association", "geometry", "quantity_units", "density_units", "component_count",
                     "components_included", "host_budget_accepted", "interpretation", "limits", "work_limits", "components"}, "ledger")
    _same(document["schema"], LEDGER_SCHEMA, "ledger/schema")
    _same(document["components_included"], True, "ledger/components_included")
    _same(document["host_budget_accepted"], False, "ledger/host_budget_accepted")
    _same(document["interpretation"], "exact aggregation of qualified represented component amounts; SI projection flags are independent of source qualification", "ledger/interpretation")
    _same(document["limits"], _LEDGER_LIMITS, "ledger/limits")
    _same(document["work_limits"], {"max_ledger_groups": MAX_GROUPS, "max_components": MAX_COMPONENTS}, "ledger/work_limits")
    _same(document["quantity_units"], {name: "parent" if name == "released_number" else "J" for name in QUANTITIES}, "ledger/quantity_units")
    _same(document["density_units"], {name: "parent m^-2" if name == "released_number" else "J m^-2" for name in QUANTITIES}, "ledger/density_units")
    association = document["source_association"]
    _keys(association, {"status", "provenance", "parameter_provenance", "log_mass_tolerance", "surface_attributes", "source_reservoirs", "selections"}, "ledger/source_association")
    _text(association["status"], "ledger/source_association/status")
    if association["status"] not in {"manufactured", "reconstructed", "sourced", "unvalidated", "owner_supplied"}:
        _fail("ledger/source_association/status", "unsupported declared status")
    for name in ("provenance", "parameter_provenance"):
        _text(association[name], "ledger/source_association/"+name)
    _float(association["log_mass_tolerance"], "ledger/source_association/log_mass_tolerance", positive=True)
    if type(association["surface_attributes"]) is not dict:
        _fail("ledger/source_association/surface_attributes", "plain object required")
    geometry = document["geometry"]
    _keys(geometry, {"x_edges_m", "z_edges_m", "time_edges_s", "width_m", "current_length_m"}, "ledger/geometry")
    x, z, times = (_edges(geometry[name], "ledger/geometry/"+name, time=name == "time_edges_s") for name in ("x_edges_m", "z_edges_m", "time_edges_s"))
    dx, dz = ([Fraction(b)-Fraction(a) for a, b in zip(values[:-1], values[1:])] for values in (x, z))
    width, length = Fraction(x[-1])-Fraction(x[0]), Fraction(z[-1])-Fraction(z[0])
    _same(geometry["width_m"], _number(width), "ledger/geometry/width_m")
    _same(geometry["current_length_m"], _number(length), "ledger/geometry/current_length_m")
    pools, selections = association["source_reservoirs"], association["selections"]
    if type(pools) is not list or not 0 < len(pools) <= MAX_COMPONENTS:
        _fail("ledger/source_association/source_reservoirs", "nonempty bounded pool list required")
    for index, pool in enumerate(pools):
        _pool(pool, f"ledger/source_association/source_reservoirs/{index}")
    if len({pool["name"] for pool in pools}) != len(pools):
        _fail("ledger/source_association/source_reservoirs", "duplicate reservoir name")
    if type(selections) is not list or not 0 < len(selections) <= 10_000:
        _fail("ledger/source_association/selections", "nonempty bounded selection list required")
    interval_selections = {}
    for index, selection in enumerate(selections):
        path = f"ledger/source_association/selections/{index}"
        _keys(selection, {"interval_index", "start_s", "stop_s"}, path)
        interval = _index(selection["interval_index"], len(times)-1, path+"/interval_index")
        _float(selection["start_s"], path+"/start_s", nonnegative=True)
        _float(selection["stop_s"], path+"/stop_s", positive=True)
        if not selection["start_s"] < selection["stop_s"] <= float(times[interval+1]-times[interval]):
            _fail(path, "selection is empty or outside its canonical parent interval")
        interval_selections.setdefault(interval, []).append(selection)
    for values in interval_selections.values():
        ordered = sorted(values, key=lambda item: (item["start_s"], item["stop_s"]))
        if any(b["start_s"] < a["stop_s"] for a, b in zip(ordered[:-1], ordered[1:])):
            _fail("ledger/source_association/selections", "overlap or duplicate")
    nx, nz, nr, ns = len(dx), len(dz), len(pools), len(selections)
    states = list(dict.fromkeys(pool["state_id"] for pool in pools))
    state_indices = {state: index for index, state in enumerate(states)}
    group_count = 2*(1+nx+nz+len(states)+nr+ns+nz*len(states)*nx)+nz*len(states)
    slots = nx*nz*nr*ns
    if slots > MAX_COMPONENTS or group_count > MAX_GROUPS:
        _fail("ledger", "component slot/group work limit exceeded")
    canonical = [[float((z[j+1]-z[j])*(x[i+1]-x[i])) for i in range(nx)] for j in range(nz)]
    if any(not math.isfinite(value) or value <= 0 for row in canonical for value in row):
        _fail("ledger/geometry", "canonical patch area must be finite and positive")
    cq = [[Fraction(value) for value in row] for row in canonical]
    full_canonical = sum((sum(row, _ZERO) for row in cq), _ZERO)
    width_canonical = [sum((row[i] for row in cq), _ZERO) for i in range(nx)]
    current_canonical = [sum(row, _ZERO) for row in cq]
    total = _Totals()
    widths, currents, reservoirs, selected = ([_Totals() for _ in range(count)] for count in (nx, nz, nr, ns))
    state_totals = {state: _Totals() for state in states}
    patches = {(j, state, i): _Totals() for j in range(nz) for state in states for i in range(nx)}
    bases = {(j, state): _Totals() for j in range(nz) for state in states}
    components = document["components"]
    if type(components) is not list or len(components) > slots:
        _fail("ledger/components", "bounded full component list required")
    _same(document["component_count"], len(components), "ledger/component_count")
    keys, catalog = [], []
    temperatures, original_inventory = {}, {}
    for index, component in enumerate(components):
        path = f"ledger/components/{index}"
        _keys(component, _ASSOCIATION_KEYS|{"source_window", "geometric_patch_area_m2", "canonical_minus_geometric_patch_area_m2"}, path)
        _same(component["component_index"], index, path+"/component_index")
        s = _index(component["selection_index"], ns, path+"/selection_index")
        r = _index(component["reservoir_index"], nr, path+"/reservoir_index")
        j = _index(component["current_index"], nz, path+"/current_index")
        i = _index(component["width_index"], nx, path+"/width_index")
        selection, pool = selections[s], pools[r]
        interval = selection["interval_index"]
        for name, expected in (("interval_index", interval), ("interval_origin_s", times[interval]), ("reservoir_name", pool["name"]),
                               ("state_id", pool["state_id"]), ("x_bounds_m", x[i:i+2]), ("z_bounds_m", z[j:j+2])):
            _same(component[name], expected, path+"/"+name)
        window = component["source_window"]
        _physical_window(window, pool, selection, times[interval], float(times[interval+1]-times[interval]), canonical[j][i], association, path+"/source_window")
        key = (s, r, j, i)
        keys.append(key)
        temperature_key = (interval, j, i)
        if temperature_key in temperatures:
            _same(window["temperature_k"], temperatures[temperature_key], path+"/source_window/temperature_k")
        else:
            temperatures[temperature_key] = window["temperature_k"]
        original_key = (interval, r, j, i)
        original_declaration = {"initial_coverage_m2": window["initial_coverage_m2"],
                                "initial_number": window["quantities"]["initial_number"],
                                "coupling": window["time_window"]["coupling"],
                                "log_initial_rate_s": window["time_window"]["log_initial_rate_s"],
                                "parent_depletion": window["time_window"]["parent_depletion"],
                                "log_parent_depletion": window["time_window"]["log_parent_depletion"]}
        if original_key in original_inventory:
            _same(original_declaration, original_inventory[original_key], path+"/original_interval_declaration")
        else:
            original_inventory[original_key] = original_declaration
        _same(component["geometric_patch_area_m2"], _number(dx[i]*dz[j]), path+"/geometric_patch_area_m2")
        _same(component["canonical_minus_geometric_patch_area_m2"], _number(cq[j][i]-dx[i]*dz[j]), path+"/canonical_minus_geometric_patch_area_m2")
        association_record = {name: component[name] for name in _ASSOCIATION_KEYS}
        values, missing = {}, {}
        for name in QUANTITIES:
            record = window["quantities"][name]
            if record["qualified"]:
                values[name] = Fraction(record["value"])
            else:
                values[name] = None
                missing[name] = len(catalog)
                catalog.append({"catalog_index": len(catalog), **association_record, "quantity_name": name,
                                "positive": not record["exact_zero"], "representation": record})
        for accumulator in (total, widths[i], currents[j], state_totals[pool["state_id"]], reservoirs[r], selected[s],
                            patches[j, pool["state_id"], i], bases[j, pool["state_id"]]):
            accumulator.add(values, missing, window["temperature_k"])
    if keys != sorted(set(keys)):
        _fail("ledger/components", "component associations must be unique and canonically ordered")
    skips, skip_keys = document["skipped_zero_initials"], []
    if type(skips) is not list or len(skips) > slots:
        _fail("ledger/skipped_zero_initials", "bounded skip list required")
    for index, skip in enumerate(skips):
        path = f"ledger/skipped_zero_initials/{index}"
        _keys(skip, {"selection_index", "interval_index", "reservoir_index", "current_index", "width_index", "reason"}, path)
        s = _index(skip["selection_index"], ns, path+"/selection_index")
        r = _index(skip["reservoir_index"], nr, path+"/reservoir_index")
        j = _index(skip["current_index"], nz, path+"/current_index")
        i = _index(skip["width_index"], nx, path+"/width_index")
        _same(skip["interval_index"], selections[s]["interval_index"], path+"/interval_index")
        _same(skip["reason"], "exact_zero_canonical_initial_inventory", path+"/reason")
        if (selections[s]["interval_index"], r, j, i) in original_inventory:
            _fail(path, "the same original initial inventory cannot be both positive and an exact-zero skip")
        skip_keys.append((s, r, j, i))
    if skip_keys != sorted(set(skip_keys)) or set(keys)&set(skip_keys) or len(keys)+len(skip_keys) != slots:
        _fail("ledger", "components and declared zero skips must partition all original slots exactly once")
    _same(document["unqualified_contributors"], catalog, "ledger/unqualified_contributors")
    area = width*length
    specs = {
        "total": [(total, {}, area, full_canonical, "full physical face once")],
        "by_width": [(widths[i], {"width_index": i}, dx[i]*length, width_canonical[i], "width bin times full current length once") for i in range(nx)],
        "by_current": [(currents[j], {"current_index": j}, width*dz[j], current_canonical[j], "full width times current bin once") for j in range(nz)],
        "by_state": [(state_totals[state], {"state_id": state}, area, full_canonical, "full physical face once per declared exact state ID") for state in states],
        "by_reservoir": [(reservoirs[r], {"reservoir_index": r, "reservoir_name": pool["name"], "state_id": pool["state_id"]}, area, full_canonical,
                           "full physical face once per declared reservoir") for r, pool in enumerate(pools)],
        "by_selection": [(selected[s], {"selection_index": s}, area, full_canonical, "full physical face once per original selection") for s in range(ns)],
        "patch_state": [(patches[j, state, i], {"current_index": j, "state_id": state, "width_index": i}, dx[i]*dz[j], cq[j][i],
                          "one physical patch once per declared exact state ID") for j in range(nz) for state in states for i in range(nx)],
    }
    for hypothesis in ("resolved", "matched_uniform_width"):
        _keys(document[hypothesis], specs, "ledger/"+hypothesis)
        for family, entries in specs.items():
            records = [document[hypothesis][family]] if family == "total" else document[hypothesis][family]
            if type(records) is not list or len(records) != len(entries):
                _fail("ledger/"+hypothesis+"/"+family, "missing or extra geometry/state group")
            for index, entry in enumerate(entries):
                accumulator, key, exact_area, canonical_area, definition = entry
                factor, reference = _ONE, None
                if hypothesis == "matched_uniform_width" and family == "by_width":
                    accumulator, factor, reference = total, dx[index]/width, "resolved/total"
                elif hypothesis == "matched_uniform_width" and family == "patch_state":
                    j, state, i = key["current_index"], key["state_id"], key["width_index"]
                    base_index = j*len(states)+state_indices[state]
                    accumulator, factor, reference = bases[j, state], dx[i]/width, f"uniform_width_bases/{base_index}"
                _same(records[index], _group(accumulator, key, exact_area, canonical_area, definition, factor=factor, reference=reference),
                      "ledger/"+hypothesis+"/"+family+f"/{index}")
    base_specs = [(bases[j, state], {"current_index": j, "state_id": state}, width*dz[j], current_canonical[j],
                   "full width times one current bin, retaining original component mixture") for j in range(nz) for state in states]
    if type(document["uniform_width_bases"]) is not list or len(document["uniform_width_bases"]) != len(base_specs):
        _fail("ledger/uniform_width_bases", "missing or extra mixture base")
    for index, entry in enumerate(base_specs):
        _same(document["uniform_width_bases"][index], _group(*entry), f"ledger/uniform_width_bases/{index}")
    if expected_source is not None:
        _bind(document, expected_source)
    digest = hashlib.sha256(json.dumps(document, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()
    return {"schema": SCHEMA, "ledger_schema": LEDGER_SCHEMA, "validation_identity": VALIDATION_IDENTITY,
            "algebraically_consistent": True, "source_association_bound": expected_source is not None,
            "component_count": len(components), "skipped_zero_initial_count": len(skips), "group_count": group_count,
            "unqualified_contributor_count": len(catalog), "total_quantity_availability": {name: not total.missing[name] for name in QUANTITIES},
            "ledger_sha256": digest, "host_budget_accepted": False,
            "limits": ["Source qualifications, logs and zero-initial skips are supplied declarations, not independently solved kinetic or physical truth.",
                       "Source binding is an exact typed comparison with the separately supplied conditioned-source manifest; its external provenance remains the caller's responsibility.",
                       "Algebraic admission does not qualify a rounded SI projection or accept a host transaction, wall debit, gap evolution, molecular state data or photon response."],
            "work_limits": {"max_components": MAX_COMPONENTS, "max_groups": MAX_GROUPS, "max_json_nodes": MAX_JSON_NODES,
                            "max_json_depth": MAX_JSON_DEPTH, "max_rational_digits": MAX_RATIONAL_DIGITS}}
