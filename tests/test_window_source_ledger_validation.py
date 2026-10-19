"""New algebraic admission checks with fresh manufactured component evidence."""
import ast
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_release import TemkinReservoir
from mhdlab.temkin_windowed_source import SourceWindowSelection, condition_temkin_planar_source
from mhdlab.window_source_ledger import window_source_ledger
from mhdlab.window_source_ledger_validation import (
    LedgerAdmissionError, MAX_RATIONAL_DIGITS, SCHEMA, validate_window_source_ledger,
)


def fixture(*, scale=1., empty=False, late=False, overflow=False):
    rate, maximum = (1., 1e100) if late else ((100., 1e308) if overflow else (.7, 1e6))
    states = ["M:A"] if late or overflow else ["M:A", "M:A", "M:empty"]
    pools = [TemkinReservoir(f"pool-{i}", state, rate, 0., 0. if late or overflow else .4,
                             0. if late or overflow else .2+i*.1, 0., maximum, maximum,
                             "Manufactured admission kinetics", "Manufactured independent costs") for i, state in enumerate(states)]
    initial = ([[[1e-100, 1e100]]] if late else [[[1e308, 1e308]]] if overflow
               else [[[20.*scale, 0.]], [[0., 10.*scale]], [[0., 0.]]])
    if empty:
        initial = np.zeros((len(pools), 1, 2))
    times = [0., 1000.] if late else [0., 1.]
    x, z = [0., 1., 2.] if late or overflow else [-1., 0., 2.], [0., 1.] if overflow else [0., 3.]
    temperature = np.array([[[300., 600.]], [[400., 700.]]])
    base = FlashSurfaceHistory([0., 1.], [0., 1.], [[300.], [300.]], "synthetic", "admission-test",
                               "Manufactured admission temperatures", "Manufactured stationary patch")
    attrs = base.attributes
    attrs["admission_evidence"] = {"label": "manufactured", "nested": [1, True, 0.]}
    surface = PlanarSurfaceHistory(times, x, z, temperature, attrs)
    history = evolve_temkin_planar_release(surface, initial, pools,
                                           parameter_provenance="Manufactured admission component inputs")
    selections = ([SourceWindowSelection(0, 800., 801.)] if late else
                  [SourceWindowSelection(0, 0., 1.)] if overflow else
                  [SourceWindowSelection(0, 0., .25), SourceWindowSelection(0, .5, .9)])
    source = condition_temkin_planar_source(history, selections, log_mass_tolerance=1e-7,
                                            status="manufactured", provenance="Manufactured admission source windows")
    return window_source_ledger(source), source.manifest()


@pytest.fixture(scope="module")
def evidence():
    return fixture()


def test_bound_and_unbound_admission_replay_exact_algebra_with_distinct_trust_flags(evidence):
    ledger, source = evidence
    unbound = validate_window_source_ledger(ledger)
    bound = validate_window_source_ledger(ledger, expected_source=source)
    assert bound["schema"] == SCHEMA and bound["algebraically_consistent"] is True
    assert bound["source_association_bound"] is True and unbound["source_association_bound"] is False
    assert bound["host_budget_accepted"] is False
    assert bound["component_count"] == len(ledger["components"]) == 4
    assert bound["skipped_zero_initial_count"] == 8
    assert bound["group_count"] == 2*(1+2+1+2+3+2+1*2*2)+1*2
    assert bound["total_quantity_availability"] == {name: True for name in ledger["resolved"]["total"]["quantities"]}
    canonical = json.dumps(ledger, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    assert bound["ledger_sha256"] == hashlib.sha256(canonical).hexdigest()
    json.dumps(bound, allow_nan=False)


def test_two_distinct_self_consistent_ledgers_need_external_source_binding(evidence):
    _, source_a = evidence
    ledger_b, source_b = fixture(scale=2.)
    assert validate_window_source_ledger(ledger_b)["algebraically_consistent"]
    assert validate_window_source_ledger(ledger_b, expected_source=source_b)["source_association_bound"]
    with pytest.raises(LedgerAdmissionError, match="binding"):
        validate_window_source_ledger(ledger_b, expected_source=source_a)


def test_coherent_declared_zero_source_is_not_a_proof_of_original_positive_inventory(evidence):
    ledger, source = fixture(empty=True)
    summary = validate_window_source_ledger(ledger, expected_source=source)
    assert summary["component_count"] == 0 and summary["skipped_zero_initial_count"] == 12
    assert summary["total_quantity_availability"]["released_number"]
    with pytest.raises(LedgerAdmissionError, match="binding"):
        validate_window_source_ledger(ledger, expected_source=evidence[1])


def rename(value, old, new):
    if type(value) is dict:
        return {key: rename(item, old, new) for key, item in value.items()}
    if type(value) is list:
        return [rename(item, old, new) for item in value]
    return new if type(value) is str and value == old else value


def test_consistent_exact_state_rename_is_algebraic_but_does_not_match_original_identity(evidence):
    ledger = rename(evidence[0], "M:A", "M:renamed")
    assert validate_window_source_ledger(ledger)["algebraically_consistent"]
    with pytest.raises(LedgerAdmissionError, match="binding"):
        validate_window_source_ledger(ledger, expected_source=evidence[1])


def test_unqualified_positive_fields_keep_partial_subtotals_and_catalog(evidence):
    ledger, source = fixture(late=True)
    summary = validate_window_source_ledger(ledger, expected_source=source)
    assert not summary["total_quantity_availability"]["released_number"]
    assert not summary["total_quantity_availability"]["translational_energy_j"]
    assert summary["total_quantity_availability"]["internal_energy_j"]
    assert summary["unqualified_contributor_count"] > 0
    tampered = deepcopy(ledger)
    value = tampered["resolved"]["total"]["quantities"]["released_number"]
    value.update(available=True, **value["qualified_subtotal"])
    value["per_area"].update(available=True, **value["per_area"]["qualified_subtotal"])
    value["unqualified"] = {"count": 0, "catalog_indices": []}
    tampered["unqualified_contributors"] = []
    with pytest.raises(LedgerAdmissionError):
        validate_window_source_ledger(tampered)


def test_complete_available_sum_with_overflowing_display_is_valid_not_zero():
    ledger, source = fixture(overflow=True)
    quantity = ledger["resolved"]["total"]["quantities"]["released_number"]
    assert quantity["available"] and quantity["value"] is None and quantity["projection"]["overflowed"]
    result = validate_window_source_ledger(ledger, expected_source=source)
    assert result["total_quantity_availability"]["released_number"]


def test_summary_only_ledger_is_not_admissible(evidence):
    ledger = deepcopy(evidence[0])
    del ledger["components"]
    ledger["components_included"] = False
    with pytest.raises(LedgerAdmissionError):
        validate_window_source_ledger(ledger)


def test_admission_is_stdlib_only_and_never_calls_scientific_producers(evidence, monkeypatch):
    import mhdlab.temkin_planar as planar
    import mhdlab.temkin_source_window as physical
    import mhdlab.window_source_ledger as producer
    import mhdlab.window_source_ledger_validation as validator
    imports = ast.parse(Path(validator.__file__).read_text())
    names = {node.module for node in ast.walk(imports) if isinstance(node, ast.ImportFrom)}
    names |= {alias.name for node in ast.walk(imports) if isinstance(node, ast.Import) for alias in node.names}
    assert names <= {"__future__", "fractions", "hashlib", "json", "math", "re", "sys"}
    reject = lambda *a, **k: pytest.fail("admission must not replay a source/law/history/producer")
    monkeypatch.setattr(planar, "evolve_temkin_planar_release", reject)
    monkeypatch.setattr(physical.FrozenTemkinSourceWindow, "__post_init__", reject)
    monkeypatch.setattr(producer, "window_source_ledger", reject)
    assert validator.validate_window_source_ledger(evidence[0], expected_source=evidence[1])["source_association_bound"]


def path_get(value, path):
    for key in path:
        value = value[key]
    return value


@pytest.mark.parametrize("path,value", [
    (("host_budget_accepted",), 0), (("host_budget_accepted",), True),
    (("components_included",), 1), (("component_count",), 4.),
    (("components", 0, "component_index"), False),
    (("resolved", "total", "quantities", "released_number", "available"), 1),
    (("resolved", "total", "quantities", "released_number", "value"), 0.),
    (("resolved", "total", "quantities", "binding_energy_j", "projection", "exact_zero"), 0),
    (("resolved", "total", "residuals", "count_times_2kbt_j", "unqualified", "count"), False),
    (("matched_uniform_width", "by_width", 0, "quantities", "released_number", "unqualified", "source_group_ref"), "resolved/total/quantities/translational_energy_j"),
    (("matched_uniform_width", "by_width", 0, "redistribution", "source_group_ref"), "uniform_width_bases/0"),
    (("components", 0, "interval_origin_s"), -0.),
    (("components", 0, "source_window", "temperature_k"), 301.),
    (("components", 0, "source_window", "start_s"), .1),
    (("components", 0, "source_window", "reservoir", "state_id"), "M:wrong"),
    (("components", 0, "source_window", "quantities", "released_number", "qualified"), 1),
    (("components", 0, "source_window", "quantities", "released_number", "exact_zero"), True),
    (("components", 0, "source_window", "time_window", "numerical_diagnostics", "host_budget_accepted"), True),
])
def test_typed_availability_projection_time_and_state_tampers_reject(evidence, path, value):
    ledger = deepcopy(evidence[0])
    target = path_get(ledger, path[:-1])
    target[path[-1]] = value
    with pytest.raises(LedgerAdmissionError):
        validate_window_source_ledger(ledger)


@pytest.mark.parametrize("numerator,denominator", [("0", "2"), ("-0", "1"), ("+1", "1"), ("01", "1"),
                                                    ("1", "01"), ("1", "-1"), ("1", "0"), ("2", "2"),
                                                    (1, "1"), ("9"*(MAX_RATIONAL_DIGITS+1), "1")])
def test_noncanonical_or_oversized_rational_records_reject(evidence, numerator, denominator):
    ledger = deepcopy(evidence[0])
    ledger["resolved"]["total"]["quantities"]["released_number"]["exact_sum"] = {"numerator": numerator, "denominator": denominator}
    with pytest.raises(LedgerAdmissionError):
        validate_window_source_ledger(ledger)


def test_positive_zero_projection_sign_and_empty_geometry_state_group_are_checked(evidence):
    ledger = deepcopy(evidence[0])
    empty = next(row for row in ledger["resolved"]["by_state"] if row["key"]["state_id"] == "M:empty")
    assert empty["quantities"]["released_number"]["value"].hex() == 0.0.hex()
    empty["quantities"]["released_number"]["value"] = -0.
    with pytest.raises(LedgerAdmissionError, match="signed zero"):
        validate_window_source_ledger(ledger)
    ledger = deepcopy(evidence[0])
    ledger["resolved"]["by_state"] = ledger["resolved"]["by_state"][:-1]
    with pytest.raises(LedgerAdmissionError, match="group"):
        validate_window_source_ledger(ledger)


@pytest.mark.parametrize("kind", ["extra_top", "extra_projection", "duplicate_component", "missing_component", "missing_skip", "duplicate_skip", "extra_catalog", "missing_base", "duplicate_width", "wrong_area"])
def test_exact_keys_and_complete_group_component_catalog_partition_required(evidence, kind):
    ledger = deepcopy(evidence[0])
    if kind == "extra_top":
        ledger["extra"] = None
    elif kind == "extra_projection":
        ledger["resolved"]["total"]["quantities"]["released_number"]["projection"]["extra"] = False
    elif kind == "duplicate_component":
        ledger["components"][1] = deepcopy(ledger["components"][0])
        ledger["components"][1]["component_index"] = 1
    elif kind == "missing_component":
        ledger["components"].pop()
        ledger["component_count"] -= 1
    elif kind == "missing_skip":
        ledger["skipped_zero_initials"].pop()
    elif kind == "duplicate_skip":
        ledger["skipped_zero_initials"][1] = deepcopy(ledger["skipped_zero_initials"][0])
    elif kind == "extra_catalog":
        ledger["unqualified_contributors"].append({})
    elif kind == "missing_base":
        ledger["uniform_width_bases"].pop()
    elif kind == "duplicate_width":
        ledger["resolved"]["by_width"][1] = deepcopy(ledger["resolved"]["by_width"][0])
    elif kind == "wrong_area":
        ledger["resolved"]["by_state"][0]["area"]["exact_m2"] = {"numerator": "18", "denominator": "1"}
    with pytest.raises(LedgerAdmissionError):
        validate_window_source_ledger(ledger)


@pytest.mark.parametrize("value", [np.float64(0.), np.int64(1), np.ma.array(1), (1, 2), {1: "bad"}, float("nan"), float("inf")])
def test_non_json_or_nonfinite_nested_metadata_rejects(evidence, value):
    ledger = deepcopy(evidence[0])
    ledger["source_association"]["surface_attributes"]["invalid"] = value
    with pytest.raises(LedgerAdmissionError):
        validate_window_source_ledger(ledger)


def test_repeated_original_inventory_and_positive_versus_skip_conflicts_reject(evidence):
    ledger = deepcopy(evidence[0])
    # Repeated selections refer to the same original frozen interval/patch.
    repeated = ledger["components"][2]
    repeated["source_window"]["initial_coverage_m2"] *= 2
    repeated["source_window"]["quantities"]["initial_coverage_m2"]["value"] *= 2
    with pytest.raises(LedgerAdmissionError, match="original_interval_declaration"):
        validate_window_source_ledger(ledger)
    ledger = deepcopy(evidence[0])
    ledger["skipped_zero_initials"][4] = {"selection_index": 1, "interval_index": 0, "reservoir_index": 0,
                                         "current_index": 0, "width_index": 0, "reason": "exact_zero_canonical_initial_inventory"}
    with pytest.raises(LedgerAdmissionError, match="positive and an exact-zero"):
        validate_window_source_ledger(ledger)


def test_expected_source_requires_original_schema_and_strict_typed_components(evidence):
    source = deepcopy(evidence[1])
    source["components"][0]["interval_origin_s"] = -0.
    with pytest.raises(LedgerAdmissionError, match="binding"):
        validate_window_source_ledger(evidence[0], expected_source=source)
    source = deepcopy(evidence[1])
    source["host_budget_accepted"] = 0
    with pytest.raises(LedgerAdmissionError):
        validate_window_source_ledger(evidence[0], expected_source=source)


@pytest.mark.parametrize("field,value", [("time_coordinate", "global photon gate"), ("limits", []), ("limits", {"wrong": True})])
def test_expected_source_time_convention_and_limit_declarations_are_bound(evidence, field, value):
    source = deepcopy(evidence[1])
    source[field] = value
    with pytest.raises(LedgerAdmissionError):
        validate_window_source_ledger(evidence[0], expected_source=source)


@pytest.mark.parametrize("status", [[], {}, 1, True])
def test_status_type_rejects_as_admission_error(evidence, status):
    ledger = deepcopy(evidence[0])
    ledger["source_association"]["status"] = status
    with pytest.raises(LedgerAdmissionError):
        validate_window_source_ledger(ledger)


def test_underflowed_canonical_patch_area_is_not_admissible_even_for_zero_slots(evidence):
    from mhdlab.window_source_ledger_validation import _number
    ledger = deepcopy(evidence[0])
    ledger["geometry"].update(x_edges_m=[0., 1e-200], z_edges_m=[0., 1e-200],
                               width_m=_number(Fraction(1e-200)), current_length_m=_number(Fraction(1e-200)))
    with pytest.raises(LedgerAdmissionError, match="canonical patch area"):
        validate_window_source_ledger(ledger)


def test_reduced_window_tolerance_must_remain_positive_after_half_allocation(evidence):
    from mhdlab.window_source_ledger_validation import _time_window
    physical = deepcopy(evidence[0]["components"][0]["source_window"])
    physical["log_mass_tolerance"] = float(np.nextafter(0., 1.))
    physical["time_window"]["log_mass_tolerance"] = 0.
    with pytest.raises(LedgerAdmissionError, match="log_mass_tolerance"):
        _time_window(physical["time_window"], physical, "window")


def cached_coupling_declaration():
    # Literal reduced-window declaration from preserved snapshot8, window0,
    # component112. Its caches come from independent log/exp rounding.
    duration = 2.3437499999999776e-10
    value = {
        "schema": "mhdlab-reduced-temkin-emission-window-v1", "coupling": 2.785835717475691,
        "log_initial_rate_s": 14.398722238710507, "parent_duration_s": duration,
        "start_s": 0., "stop_s": duration, "local_duration_s": duration,
        "start_depletion": 0., "local_depletion": 0.00041970391116034184,
        "parent_depletion": 0.00041970391116034184, "log_start_depletion": None,
        "log_local_depletion": -7.775961068725234, "log_parent_depletion": -7.775961068725234,
        "log_local_coupling": 1.0245479064530423, "local_coupling": 2.7858357174756914,
        "log_local_rate_s": 14.398722238710507, "log_source_fraction": -7.776170913341174,
        "log_parent_conditional_probability": 0., "log_mass_tolerance": 5e-9,
        "time_coordinate": "local delay from original source-window start; keep start and delay separately",
        "numerical_diagnostics": {
            "qualification_identity": "estimated_monotonic_log_depletion_envelopes_with_positive_tail_correction_v1",
            "quadrature_assumption": "16 times SciPy absolute quadrature error estimate plus explicit float64 rounding allowances; estimated envelopes, not certified bounds",
            "log_mass_tolerance": 5e-9, "estimated_log_source_fraction_error": 1.0828604219396013e-12,
            "estimated_log_parent_probability_error": 1.7435266184094407e-12, "start_depletion_error": 0.,
            "local_log_depletion_error": 8.336015630105557e-13, "parent_log_depletion_error": 5.360617924110138e-13,
            "local_coupling_approximated_as_zero": False, "log_local_coupling_approximation_error_bound": None,
            "local_coupling_approximation_error_bound": 0., "start_rebase_log_mass_sensitivity_bound": 3.785835717475691,
            "local_duration_relative_rounding_error": 0., "local_depletion_underflowed": False,
            "start_depletion_underflowed": False, "tail_switch_depletion": 32.,
            "certified_error_bound": False, "host_budget_accepted": False,
            "root_evidence": {
                "start": {"log_depletion": None, "depletion": 0., "log_error": 0., "absolute_error": 0., "iterations": 0, "integral_method": "zero_origin"},
                "local": {"log_depletion": -7.775961068725234, "depletion": 0.00041970391116034184,
                          "log_error": 8.336015630105557e-13, "absolute_error": 3.513569241712756e-16,
                          "iterations": 2, "integral_method": "endpoint_scaled_positive_quadrature"},
                "parent": {"log_depletion": -7.775961068725234, "depletion": 0.00041970391116034184,
                           "log_error": 5.360617924110138e-13, "absolute_error": 2.264783187248114e-16,
                           "iterations": 2, "integral_method": "endpoint_scaled_positive_quadrature"},
            },
        },
    }
    physical = {"parent_duration_s": duration, "start_s": 0., "stop_s": duration, "log_mass_tolerance": 1e-8}
    return value, physical


def test_saved_scalar_log_exp_roundtrip_can_exceed_original_coupling_by_one_ulp():
    from mhdlab.window_source_ledger_validation import _time_window
    value, physical = cached_coupling_declaration()
    assert value["local_coupling"] == math.nextafter(value["coupling"], math.inf)
    assert value["log_local_coupling"].hex() == math.log(value["coupling"]).hex()
    assert value["local_coupling"].hex() == math.exp(value["log_local_coupling"]).hex()
    _time_window(value, physical, "saved_window")


@pytest.mark.parametrize("field,change", [("local_coupling", "nextafter"), ("local_coupling", "large"),
                                          ("log_local_coupling", "nextafter"), ("coupling", "out_of_range")])
def test_upward_coupling_exception_requires_exact_declared_roundtrip(field, change):
    from mhdlab.window_source_ledger_validation import _time_window
    value, physical = cached_coupling_declaration()
    value[field] = (math.nextafter(value[field], math.inf) if change == "nextafter"
                    else value[field]*1.001 if change == "large" else 512.0001)
    with pytest.raises(LedgerAdmissionError):
        _time_window(value, physical, "saved_window")


def test_admission_summary_is_detached_and_work_limits_precede_unbounded_processing(evidence, monkeypatch):
    import mhdlab.window_source_ledger_validation as validator
    summary = validate_window_source_ledger(evidence[0], expected_source=evidence[1])
    summary["total_quantity_availability"]["released_number"] = False
    assert validate_window_source_ledger(evidence[0])["total_quantity_availability"]["released_number"]
    monkeypatch.setattr(validator, "MAX_JSON_NODES", 2)
    with pytest.raises(LedgerAdmissionError, match="work limit"):
        validator.validate_window_source_ledger(evidence[0])
