"""New exact-represented inventory/energy and unique-area ledger checks."""
from copy import copy
from dataclasses import replace
from fractions import Fraction
import json
import math
from types import MappingProxyType

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_release import TemkinReservoir
from mhdlab.temkin_windowed_source import SourceWindowSelection, condition_temkin_planar_source
from mhdlab.window_source_ledger import (
    QUANTITIES, SCHEMA, WindowSourceLedgerRangeError, window_source_ledger,
)


def q(value):
    return Fraction.from_float(float(value))


def exact(record):
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def pool(name="a", state="M:v0", **changes):
    values = dict(name=name, state_id=state, attempt_frequency_s=.7, activation_energy_ev=0.,
                  desorption_energy_ev=.4, internal_energy_ev=.2, coverage_coefficient=0.,
                  reference_coverage_m2=1e6, max_coverage_m2=1e6,
                  kinetic_provenance="Manufactured ledger kinetics", energy_provenance="Manufactured independent ledger costs")
    values.update(changes)
    return TemkinReservoir(**values)


def make_source(initial, pools, *, x=(0., 1., 3.), z=(0., 2., 5.), times=(0., 1., 2.), selections=None, temperature=300.):
    initial = np.asarray(initial, dtype=float)
    temperatures = np.full((len(times), len(z)-1, len(x)-1), temperature, dtype=float)
    if not np.isscalar(temperature):
        temperatures = np.asarray(temperature, dtype=float)
    base = FlashSurfaceHistory([0., 1.], [0., 1.], [[300.], [300.]], "synthetic", "ledger-test",
                               "Manufactured ledger temperatures", "Manufactured stationary source")
    attributes = base.attributes
    attributes["ledger_evidence"] = {"label": "manufactured", "retained": [1, 2]}
    surface = PlanarSurfaceHistory(times, x, z, temperatures, attributes)
    history = evolve_temkin_planar_release(surface, initial, pools, parameter_provenance="Manufactured ledger independent pool inputs")
    if selections is None:
        selections = [SourceWindowSelection(0, .1, .4), SourceWindowSelection(1, .2, .8)]
    return condition_temkin_planar_source(history, selections, log_mass_tolerance=1e-7,
                                          status="manufactured", provenance="Manufactured ledger source restrictions")


@pytest.fixture(scope="module")
def source():
    pools = [pool(), pool("b", "M:v0", attempt_frequency_s=1.1, internal_energy_ev=.6), pool("c", "M:v1")]
    initial = [[[100., 0.], [200., 300.]], [[0., 400.], [500., 600.]], [[0., 0.], [0., 0.]]]
    temperatures = [[[300., 500.], [700., 900.]], [[400., 600.], [800., 1000.]], [[500., 700.], [900., 1100.]]]
    return make_source(initial, pools, temperature=temperatures)


def quantity_sum(components, name):
    return sum((q(component.source_window.manifest()["quantities"][name]["value"]) for component in components), Fraction())


def matches(component, key):
    return all((component.source_window.reservoir.state_id if name == "state_id" else
                component.source_window.reservoir.name if name == "reservoir_name" else getattr(component, name)) == value
               for name, value in key.items())


def test_exact_resolved_groups_match_independent_component_partition_and_analytic_release(source):
    result = window_source_ledger(source)
    assert result["schema"] == SCHEMA
    assert result["host_budget_accepted"] is False
    assert result["source_association"]["status"] == "manufactured"
    assert result["source_association"]["parameter_provenance"] == source.parameter_provenance
    for family, rows in result["resolved"].items():
        rows = [rows] if family == "total" else rows
        for row in rows:
            selected = [component for component in source.components if matches(component, row["key"])]
            for name in QUANTITIES:
                value = row["quantities"][name]
                expected = quantity_sum(selected, name)
                assert value["available"]
                assert exact(value["exact_sum"]) == expected
                assert exact(value["qualified_subtotal"]["exact_sum"]) == expected
                assert exact(value["per_area"]["exact_sum"]) == expected/exact(row["area"]["exact_m2"])
    for component in source.components:
        window = component.source_window
        expected = (window.initial_coverage_m2*window.patch_area_m2*math.exp(-window.reservoir.attempt_frequency_s*window.start_s)
                    *(-math.expm1(-window.reservoir.attempt_frequency_s*(window.stop_s-window.start_s))))
        assert window.manifest()["quantities"]["released_number"]["value"] == pytest.approx(expected, rel=3e-14)


def test_unique_group_areas_include_zero_slots_and_combine_same_state_pools_once(source):
    result = window_source_ledger(source)
    resolved = result["resolved"]
    assert exact(resolved["total"]["area"]["exact_m2"]) == 15
    assert [exact(row["area"]["exact_m2"]) for row in resolved["by_width"]] == [5, 10]
    assert [exact(row["area"]["exact_m2"]) for row in resolved["by_current"]] == [6, 9]
    for family in ("by_state", "by_reservoir", "by_selection"):
        assert all(exact(row["area"]["exact_m2"]) == 15 for row in resolved[family])
    assert len(resolved["by_state"]) == 2 and len(resolved["by_reservoir"]) == 3
    zero = next(row for row in resolved["by_state"] if row["key"]["state_id"] == "M:v1")
    assert zero["component_count"] == 0
    for value in zero["quantities"].values():
        assert value["available"] and value["value"] == 0.
        assert value["projection"]["exact_zero"] and not value["projection"]["positive"]
    assert len(resolved["patch_state"]) == 8
    for row in resolved["patch_state"]:
        expected_area = [[2, 4], [3, 6]][row["key"]["current_index"]][row["key"]["width_index"]]
        assert exact(row["area"]["exact_m2"]) == expected_area
    assert len(result["skipped_zero_initials"]) == 12


def test_matched_width_equals_naive_exact_component_by_target_width_control(source):
    result = window_source_ledger(source)
    uniform = result["matched_uniform_width"]
    for name in ("total", "by_current", "by_state", "by_reservoir", "by_selection"):
        assert uniform[name] == result["resolved"][name]
    for family in ("by_width", "patch_state"):
        for row in uniform[family]:
            target = row["key"]["width_index"]
            factor = Fraction((1, 2)[target], 3)
            for name in QUANTITIES:
                expected = Fraction()
                for component in source.components:
                    if family == "patch_state" and (component.current_index != row["key"]["current_index"]
                            or component.source_window.reservoir.state_id != row["key"]["state_id"]):
                        continue
                    expected += q(component.source_window.manifest()["quantities"][name]["value"])*factor
                assert exact(row["quantities"][name]["exact_sum"]) == expected
            assert exact(row["redistribution"]["fraction"]) == factor
            assert row["quantities"]["released_number"]["unqualified"]["source_group_ref"]
    for hypothesis in (result["resolved"], uniform):
        total = exact(hypothesis["total"]["quantities"]["released_number"]["exact_sum"])
        for family in ("by_width", "by_current", "by_state", "by_reservoir", "by_selection", "patch_state"):
            assert sum((exact(row["quantities"]["released_number"]["exact_sum"]) for row in hypothesis[family]), Fraction()) == total


def test_signed_energy_and_count_temperature_residuals_are_exact_and_unrepaired(source):
    result = window_source_ledger(source)
    expected_count_energy = sum((q(c.source_window.manifest()["quantities"]["released_number"]["value"])*2*q(KB)*q(c.source_window.temperature_k)
                                 for c in source.components), Fraction())
    row = result["resolved"]["total"]
    quantities, residuals = row["quantities"], row["residuals"]
    expected_wall = exact(quantities["total_wall_energy_debit_j"]["exact_sum"])-sum(
        (exact(quantities[name]["exact_sum"]) for name in ("translational_energy_j", "internal_energy_j", "binding_energy_j")), Fraction())
    assert exact(residuals["wall_minus_parts_j"]["exact_sum"]) == expected_wall
    assert exact(residuals["count_times_2kbt_j"]["exact_sum"]) == expected_count_energy
    assert exact(residuals["translation_minus_count_2kbt_j"]["exact_sum"]) == exact(quantities["translational_energy_j"]["exact_sum"])-expected_count_energy
    assert expected_wall != 0


def test_canonical_patch_area_and_exact_geometry_remain_distinct_without_density_repair():
    source = make_source([[[20., 30.]]], [pool()], x=(0., .1, .3), z=(0., .7), times=(0., 1.),
                         selections=[SourceWindowSelection(0, .1, .5)])
    result = window_source_ledger(source)
    row = result["resolved"]["total"]
    geometric = q(.3)*q(.7)
    canonical = sum((q(c.source_window.patch_area_m2) for c in source.components), Fraction())
    assert exact(row["area"]["exact_m2"]) == geometric
    assert exact(row["area"]["canonical_patch_area_sum"]["exact_sum"]) == canonical
    assert exact(row["area"]["canonical_minus_exact_residual_m2"]["exact_sum"]) == canonical-geometric != 0
    assert exact(row["quantities"]["released_number"]["per_area"]["exact_sum"]) == quantity_sum(source.components, "released_number")/geometric


def test_independent_field_underflow_preserves_qualified_subtotal_and_exact_zero_cost():
    source = make_source([[[1e80, 1e100]]], [pool(attempt_frequency_s=1., reference_coverage_m2=1e100,
                         max_coverage_m2=1e100, internal_energy_ev=0., desorption_energy_ev=1e20)],
                         x=(0., 1., 2.), z=(0., 1.), times=(0., 1000.),
                         selections=[SourceWindowSelection(0, 900., 901.)])
    result = window_source_ledger(source, include_components=False)
    total = result["resolved"]["total"]
    translation = total["quantities"]["translational_energy_j"]
    assert not translation["available"] and translation["exact_sum"] is None and translation["value"] is None
    expected_subtotal = q(source.components[1].source_window.manifest()["quantities"]["translational_energy_j"]["value"])
    assert exact(translation["qualified_subtotal"]["exact_sum"]) == expected_subtotal > 0
    assert len(result["unqualified_contributors"]) == 1
    failure = result["unqualified_contributors"][0]
    assert failure["positive"] and failure["quantity_name"] == "translational_energy_j"
    assert failure["representation"]["reason"] == "positive_amount_underflow"
    assert math.isfinite(failure["representation"]["log_value"])
    assert not failure["representation"]["exact_zero"]
    for name in ("released_number", "internal_energy_j", "binding_energy_j", "total_wall_energy_debit_j"):
        assert total["quantities"][name]["available"]
    assert total["quantities"]["internal_energy_j"]["projection"]["exact_zero"]
    assert not total["residuals"]["wall_minus_parts_j"]["available"]
    assert not total["residuals"]["translation_minus_count_2kbt_j"]["available"]
    assert total["residuals"]["count_times_2kbt_j"]["available"]
    assert "components" not in result
    for row in result["matched_uniform_width"]["by_width"]:
        assert not row["quantities"]["translational_energy_j"]["available"]
        assert row["quantities"]["translational_energy_j"]["unqualified"]["count"] == 1


def test_unqualified_positive_number_is_not_reported_as_empty_zero():
    source = make_source([[[1e-100]]], [pool(attempt_frequency_s=1., reference_coverage_m2=1e-100,
                         max_coverage_m2=1e-100, internal_energy_ev=0., desorption_energy_ev=0.)],
                         x=(0., 1.), z=(0., 1.), times=(0., 1000.), selections=[SourceWindowSelection(0, 800., 801.)])
    result = window_source_ledger(source)
    value = result["resolved"]["total"]["quantities"]["released_number"]
    assert not value["available"] and value["value"] is None and value["exact_sum"] is None
    assert value["qualified_subtotal"]["value"] == 0.
    assert result["component_count"] == 1 and len(result["components"]) == 1
    assert result["unqualified_contributors"][0]["positive"]
    assert not result["resolved"]["total"]["residuals"]["count_times_2kbt_j"]["available"]


def test_unqualified_areal_fields_do_not_block_qualified_physical_amounts():
    source = make_source([[[1e-100]]], [pool(attempt_frequency_s=1., reference_coverage_m2=1e-100,
                         max_coverage_m2=1e-100, internal_energy_ev=0., desorption_energy_ev=0.)],
                         x=(0., 1e150), z=(0., 1e150), times=(0., 1000.), selections=[SourceWindowSelection(0, 800., 801.)])
    assert not source.components[0].source_window.all_budget_values_qualified
    result = window_source_ledger(source)
    assert all(row["available"] for row in result["resolved"]["total"]["quantities"].values())
    assert not result["unqualified_contributors"]
    assert result["resolved"]["total"]["quantities"]["released_number"]["value"] > 0


def test_all_zero_inventory_has_complete_zero_groups_for_entire_declared_geometry():
    source = make_source(np.zeros((1, 2, 2)), [pool()])
    result = window_source_ledger(source)
    assert not result["components"] and not result["unqualified_contributors"]
    assert len(result["skipped_zero_initials"]) == 8
    for hypothesis in ("resolved", "matched_uniform_width"):
        for family, rows in result[hypothesis].items():
            for row in ([rows] if family == "total" else rows):
                assert exact(row["area"]["exact_m2"]) > 0
                assert all(value["available"] and value["projection"]["exact_zero"] for value in row["quantities"].values())


def test_exact_uniform_width_factor_does_not_underflow_before_restoring_representable_amount():
    source = make_source([[[1e308, 0.]]], [pool(attempt_frequency_s=1., reference_coverage_m2=1e308,
                         max_coverage_m2=1e308, internal_energy_ev=0., desorption_energy_ev=0.)],
                         x=(0., 1e-180, 1e180), z=(0., 1.), times=(0., 1.), selections=[SourceWindowSelection(0, 0., 1.)])
    result = window_source_ledger(source)
    tiny = result["matched_uniform_width"]["by_width"][0]
    factor = exact(tiny["redistribution"]["fraction"])
    assert float(factor) == 0.
    for name in ("released_number", "translational_energy_j"):
        expected = exact(result["resolved"]["total"]["quantities"][name]["exact_sum"])*factor
        assert exact(tiny["quantities"][name]["exact_sum"]) == expected > 0
        assert tiny["quantities"][name]["value"] > 0
    assert exact(result["resolved"]["total"]["area"]["canonical_minus_exact_residual_m2"]["exact_sum"]) > 0


def test_complete_exact_aggregate_can_overflow_display_without_becoming_unqualified_source():
    source = make_source([[[1e308, 1e308]]], [pool(attempt_frequency_s=100., reference_coverage_m2=1e308,
                         max_coverage_m2=1e308, internal_energy_ev=0., desorption_energy_ev=0.)],
                         x=(0., 1., 2.), z=(0., 1.), times=(0., 1.), selections=[SourceWindowSelection(0, 0., 1.)])
    result = window_source_ledger(source)
    value = result["resolved"]["total"]["quantities"]["released_number"]
    assert value["available"] and value["exact_sum"] is not None
    assert value["value"] is None and value["projection"]["overflowed"]
    assert value["projection"]["positive"] and not result["unqualified_contributors"]
    assert value["per_area"]["value"] > 0


def test_projection_flags_distinguish_exact_zero_underflow_and_coarse_subnormal_rounding():
    from mhdlab.window_source_ledger import _number
    tiny = q(float(np.nextafter(0., 1.)))
    zero, small, rounded = _number(Fraction()), _number(tiny/4), _number(3*tiny/4)
    assert zero["projection"]["exact_zero"] and not zero["projection"]["positive"]
    assert small["value"] == 0. and small["projection"]["underflowed"] and small["projection"]["positive"]
    assert not small["projection"]["exact_zero"]
    assert rounded["value"] == float(tiny) and rounded["projection"]["rounded"] and rounded["projection"]["subnormal"]
    assert exact(rounded["projection"]["exact_rounding_error"]) == tiny/4


def test_ledger_is_detached_json_and_preserves_original_component_time_cost_associations(source):
    result = window_source_ledger(source)
    json.dumps(result, allow_nan=False)
    for row, component in zip(result["components"], source.components):
        assert row["interval_origin_s"] == component.interval_origin_s
        assert row["selection_index"] == component.selection_index
        assert row["source_window"] == component.source_window.manifest()
    result["source_association"]["surface_attributes"]["ledger_evidence"]["retained"][0] = 99
    result["components"][0]["source_window"]["quantities"]["released_number"]["value"] = -1.
    assert source.surface_attributes["ledger_evidence"]["retained"][0] == 1
    assert source.components[0].source_window.manifest()["quantities"]["released_number"]["value"] > 0


def test_single_canonical_validation_and_factored_linear_component_work_no_evolution(source, monkeypatch):
    import mhdlab.temkin_planar as planar
    import mhdlab.window_source_ledger as ledger
    from mhdlab.temkin_windowed_source import TemkinPlanarSourceWindows
    validated = TemkinPlanarSourceWindows.validated_copy
    add = ledger._Accumulator.add
    counters = {"validated": 0, "add": 0}
    def validate_once(self):
        counters["validated"] += 1
        return validated(self)
    def add_once(self, *args):
        counters["add"] += 1
        return add(self, *args)
    monkeypatch.setattr(TemkinPlanarSourceWindows, "validated_copy", validate_once)
    monkeypatch.setattr(ledger._Accumulator, "add", add_once)
    monkeypatch.setattr(planar, "evolve_temkin_planar_release", lambda *a, **k: pytest.fail("no source evolution"))
    ledger.window_source_ledger(source, include_components=False)
    assert counters == {"validated": 1, "add": 8*len(source.components)}


@pytest.mark.parametrize("field", ["value", "qualified", "exact_zero", "reason", "log_value"])
def test_forged_physical_quantity_records_reject(source, field):
    forged = copy(source)
    component = copy(source.components[0])
    window = copy(component.source_window)
    records = dict(window._quantities)
    record = copy(records["released_number"])
    changes = {"value": 0., "qualified": 1, "exact_zero": True, "reason": "forged", "log_value": 0.}
    object.__setattr__(record, field, changes[field])
    records["released_number"] = record
    object.__setattr__(window, "_quantities", MappingProxyType(records))
    object.__setattr__(component, "source_window", window)
    object.__setattr__(forged, "components", (component, *source.components[1:]))
    with pytest.raises(ValueError, match="constructor replay"):
        window_source_ledger(forged)


def test_forged_private_time_cache_and_canonical_associations_reject(source):
    forged = copy(source)
    component = copy(source.components[0])
    window = copy(component.source_window)
    law = copy(window.time_law)
    object.__setattr__(law, "_query_base_error", 0.)
    object.__setattr__(window, "time_law", law)
    object.__setattr__(component, "source_window", window)
    object.__setattr__(forged, "components", (component, *source.components[1:]))
    with pytest.raises(ValueError, match="timing cache"):
        window_source_ledger(forged)
    forged = copy(source)
    object.__setattr__(forged, "parameter_provenance", "forged")
    with pytest.raises(ValueError):
        window_source_ledger(forged)


@pytest.mark.parametrize("value", [None, True, {}, np.ma.array(1.)])
def test_noncanonical_source_rejects(value):
    with pytest.raises(ValueError, match="TemkinPlanarSourceWindows"):
        window_source_ledger(value)


@pytest.mark.parametrize("value", [1, 0, None, "true", np.bool_(True), np.ma.array(True)])
def test_include_components_requires_explicit_boolean(source, value):
    with pytest.raises(ValueError, match="include_components"):
        window_source_ledger(source, include_components=value)


def test_group_work_cap_precedes_scalar_replay(source, monkeypatch):
    import mhdlab.window_source_ledger as ledger
    monkeypatch.setattr(ledger, "MAX_LEDGER_GROUPS", 1)
    monkeypatch.setattr(ledger, "_replayed_window", lambda *a: pytest.fail("work cap precedes scalar replay"))
    with pytest.raises(WindowSourceLedgerRangeError, match="work limit"):
        ledger.window_source_ledger(source)
