"""Independent physical counts and energy for a law-specific source window."""

from dataclasses import FrozenInstanceError, asdict, replace
from decimal import Decimal, localcontext
import json
import math

import numpy as np
import pytest

from mhdlab.temkin_release import TemkinReservoir
from mhdlab.temkin_source_window import FrozenTemkinSourceWindow
import mhdlab.temkin_source_window as source_window_module


KB = 1.380649e-23
EV = 1.602176634e-19


def reservoir(**changes):
    parameters = dict(name="manufactured_window_pool", state_id="M:manufactured-window-state",
                      attempt_frequency_s=2.5, activation_energy_ev=0., desorption_energy_ev=.6,
                      internal_energy_ev=.2, coverage_coefficient=.35, reference_coverage_m2=1e19,
                      max_coverage_m2=1e20, kinetic_provenance="Manufactured constant-hazard verification",
                      energy_provenance="Manufactured removal and internal references, independent of zero rate barrier")
    parameters.update(changes)
    return TemkinReservoir(**parameters)


def window(**changes):
    parameters = dict(reservoir=reservoir(), initial_coverage_m2=3e18, temperature_k=700.,
                      parent_duration_s=2., start_s=.3, stop_s=.9, patch_area_m2=2e-6,
                      log_mass_tolerance=1e-8, status="manufactured", provenance="Manufactured physical-window fixture")
    parameters.update(changes)
    return FrozenTemkinSourceWindow(**parameters)


def constant_hazard_reference(source):
    with localcontext() as context:
        context.prec = 80
        n0, area, attempt, left, right = map(Decimal.from_float,
            (source.initial_coverage_m2, source.patch_area_m2, source.reservoir.attempt_frequency_s,
             source.start_s, source.stop_s))
        barrier = (Decimal.from_float(source.reservoir.activation_energy_ev) * Decimal.from_float(EV)
                   / (Decimal.from_float(KB) * Decimal.from_float(source.temperature_k)))
        rate = attempt * (-barrier).exp()
        start, stop = n0 * (-rate * left).exp(), n0 * (-rate * right).exp()
        released = start * (1 - (-rate * (right-left)).exp())
        densities = {"initial_coverage_m2": n0, "remaining_start_m2": start,
                     "remaining_stop_m2": stop, "released_m2": released}
        values = {name: float(value) for name, value in densities.items()}
        values.update({name: float(value * area) for name, value in
                       (("initial_number", n0), ("remaining_start_number", start),
                        ("remaining_stop_number", stop), ("released_number", released))})
        per_particle = {"translational_energy": Decimal(2) * Decimal.from_float(KB) * Decimal.from_float(source.temperature_k),
                        "internal_energy": Decimal.from_float(EV) * Decimal.from_float(source.reservoir.internal_energy_ev),
                        "binding_energy": Decimal.from_float(EV) * Decimal.from_float(source.reservoir.desorption_energy_ev)}
        per_particle["total_wall_energy_debit"] = sum(per_particle.values())
        for name, energy in per_particle.items():
            values[name + "_j_m2"] = float(energy * released)
            values[name + "_j"] = float(energy * released * area)
        return values


def test_closed_hazard_density_counts_and_energy_with_original_inventory_and_fixed_reference():
    source = window()
    manifest = source.manifest()
    expected = constant_hazard_reference(source)
    assert manifest["schema"] == "mhdlab-temkin-source-window-v1"
    assert manifest["status"] == "manufactured"
    assert manifest["provenance"] == source.provenance
    assert manifest["reservoir"] == asdict(source.reservoir)
    assert source.reservoir.reference_coverage_m2 == 1e19
    assert source.time_law is source.time_law
    assert manifest["host_budget_accepted"] is False
    assert manifest["all_budget_values_qualified"] is True
    assert set(manifest["quantities"]) == set(expected)
    for name, reference in expected.items():
        record = manifest["quantities"][name]
        assert record["qualified"] is True
        assert record["exact_zero"] is False
        assert record["value"] == pytest.approx(reference, rel=3e-12, abs=0.)
        assert record["log_value"] == pytest.approx(math.log(reference), abs=3e-12)
        assert record["log_rounding_error"] + record["log_rounding_allowance"] <= source.log_mass_tolerance / 2
    assert expected["remaining_start_number"] < expected["initial_number"]
    for suffix in ("m2", "number"):
        remaining_start = manifest["quantities"][f"remaining_start_{suffix}"]["value"]
        remaining_stop = manifest["quantities"][f"remaining_stop_{suffix}"]["value"]
        released = manifest["quantities"][f"released_{suffix}"]["value"]
        residual = manifest[f"inventory_balance_residual_{suffix}"]
        assert residual == math.fsum((remaining_stop, -remaining_start, released))
    budget = source.represented_budget()
    assert budget["schema"] == "mhdlab-temkin-window-represented-budget-v1"
    assert budget["host_budget_accepted"] is False
    assert budget["values"] == {name: record["value"] for name, record in manifest["quantities"].items()}
    json.dumps(manifest, allow_nan=False)
    json.dumps(budget, allow_nan=False)


def test_patch_area_enters_once_without_rescaling_density_or_time_law():
    first = window()
    second = window(patch_area_m2=4 * first.patch_area_m2)
    a, b = first.manifest(), second.manifest()
    assert a["time_window"] == b["time_window"]
    for name, record in a["quantities"].items():
        factor = 1. if name.endswith("m2") else 4.
        assert b["quantities"][name]["value"] == pytest.approx(factor * record["value"], rel=3e-12)
    assert b["reservoir"] == a["reservoir"]


def test_arrhenius_rate_barrier_is_independent_of_declared_wall_and_internal_energy():
    pool = reservoir(attempt_frequency_s=2.5e6, activation_energy_ev=.3, coverage_coefficient=0.)
    source = window(reservoir=pool, parent_duration_s=4e-4, start_s=4e-5, stop_s=1e-4)
    expected = constant_hazard_reference(source)
    manifest = source.manifest()
    for name, reference in expected.items():
        assert manifest["quantities"][name]["value"] == pytest.approx(reference, rel=3e-12, abs=0.)
    changed = window(reservoir=replace(pool, desorption_energy_ev=1.2, internal_energy_ev=.8),
                     parent_duration_s=4e-4, start_s=4e-5, stop_s=1e-4).manifest()
    assert changed["time_window"] == manifest["time_window"]
    for name in ("released_number", "released_m2", "remaining_start_m2", "remaining_stop_m2"):
        assert changed["quantities"][name] == manifest["quantities"][name]
    assert changed["quantities"]["binding_energy_j"]["value"] == pytest.approx(
        2 * manifest["quantities"]["binding_energy_j"]["value"], rel=3e-12)
    assert changed["quantities"]["internal_energy_j"]["value"] == pytest.approx(
        4 * manifest["quantities"]["internal_energy_j"]["value"], rel=3e-12)


def test_disjoint_windows_partition_the_same_parent_release_and_energy():
    whole = window(start_s=0., stop_s=2.).represented_budget()["values"]
    pieces = [window(start_s=start, stop_s=stop).represented_budget()["values"]
              for start, stop in ((0., .3), (.3, .9), (.9, 2.))]
    for name in ("released_m2", "released_number", "translational_energy_j_m2", "internal_energy_j_m2",
                 "binding_energy_j_m2", "total_wall_energy_debit_j_m2", "translational_energy_j",
                 "internal_energy_j", "binding_energy_j", "total_wall_energy_debit_j"):
        assert math.fsum(piece[name] for piece in pieces) == pytest.approx(whole[name], rel=3e-12)
    assert all(piece["initial_number"] == whole["initial_number"] for piece in pieces)
    for left, right in zip(pieces[:-1], pieces[1:]):
        assert left["remaining_stop_number"] == pytest.approx(right["remaining_start_number"], rel=3e-12)


def test_underflowed_window_fraction_is_rescued_by_original_inventory_before_rounding():
    pool = reservoir(attempt_frequency_s=1., reference_coverage_m2=1e300, max_coverage_m2=1e300)
    source = window(reservoir=pool, initial_coverage_m2=1e300, parent_duration_s=1000.,
                    start_s=800., stop_s=801.)
    log_fraction = -800 + math.log(-math.expm1(-1.))
    assert math.exp(log_fraction) == 0
    manifest = source.manifest()
    assert manifest["all_budget_values_qualified"] is True
    record = manifest["quantities"]["released_number"]
    expected_log = math.log(1e300) + math.log(source.patch_area_m2) + log_fraction
    assert record["log_value"] == pytest.approx(expected_log, abs=3e-10)
    assert record["value"] > 0
    assert record["value"] == pytest.approx(math.exp(expected_log), rel=3e-10)
    assert source.represented_budget()["values"]["released_number"] == record["value"]


def test_positive_count_cannot_hide_underflowed_outgoing_energy():
    pool = reservoir(attempt_frequency_s=1., desorption_energy_ev=0., internal_energy_ev=0.)
    source = window(reservoir=pool, initial_coverage_m2=1e-300, temperature_k=1e-100,
                    parent_duration_s=1., start_s=0., stop_s=1., patch_area_m2=1.)
    manifest = source.manifest()
    number = manifest["quantities"]["released_number"]
    assert number["value"] > 0
    assert number["qualified"] is True
    energy = manifest["quantities"]["translational_energy_j"]
    assert energy["value"] == 0
    assert energy["exact_zero"] is False
    assert math.isfinite(energy["log_value"])
    assert energy["qualified"] is False
    assert manifest["all_budget_values_qualified"] is False
    with pytest.raises(ValueError):
        source.represented_budget()


def test_area_can_rescue_physical_number_without_claiming_a_represented_density_budget():
    source = window(reservoir=reservoir(attempt_frequency_s=1.), initial_coverage_m2=1.,
                    parent_duration_s=1000., start_s=800., stop_s=801., patch_area_m2=1e300)
    manifest = source.manifest()
    density = manifest["quantities"]["released_m2"]
    number = manifest["quantities"]["released_number"]
    assert density["value"] == 0
    assert density["qualified"] is False
    assert number["value"] > 0
    assert number["qualified"] is True
    expected_log = math.log(1e300) - 800 + math.log(-math.expm1(-1.))
    assert number["log_value"] == pytest.approx(expected_log, abs=3e-10)
    assert manifest["all_budget_values_qualified"] is False
    with pytest.raises(ValueError):
        source.represented_budget()


def test_overflowed_original_inventory_cannot_be_hidden_by_finite_late_release():
    pool = reservoir(attempt_frequency_s=1., reference_coverage_m2=1e300, max_coverage_m2=1e300)
    with pytest.raises(ValueError, match="initial parent"):
        window(reservoir=pool, initial_coverage_m2=1e300, patch_area_m2=1e300,
               parent_duration_s=1000., start_s=800., stop_s=801.)


def test_overflowed_outgoing_energy_retains_a_finite_log_and_never_returns_infinity():
    pool = reservoir(attempt_frequency_s=1., reference_coverage_m2=1e300, max_coverage_m2=1e300)
    source = window(reservoir=pool, initial_coverage_m2=1e300, temperature_k=1e300,
                    parent_duration_s=1., start_s=0., stop_s=1., patch_area_m2=1.)
    manifest = source.manifest()
    energy = manifest["quantities"]["translational_energy_j"]
    assert energy["value"] is None
    assert energy["qualified"] is False
    assert math.isfinite(energy["log_value"])
    assert energy["log_value"] > math.log(np.finfo(float).max)
    assert manifest["quantities"]["released_number"]["value"] > 0
    assert manifest["all_budget_values_qualified"] is False
    json.dumps(manifest, allow_nan=False)
    with pytest.raises(ValueError):
        source.represented_budget()


def test_positive_rounded_subnormal_count_requires_its_own_log_qualification():
    pool = reservoir(attempt_frequency_s=1.)
    source = window(reservoir=pool, initial_coverage_m2=float(np.nextafter(0., 1.)),
                    parent_duration_s=1., start_s=0., stop_s=1., patch_area_m2=1.)
    record = source.manifest()["quantities"]["released_number"]
    assert record["value"] > 0
    assert record["exact_zero"] is False
    assert record["log_rounding_error"] > record["log_rounding_allowance"]
    assert record["qualified"] is False
    with pytest.raises(ValueError):
        source.represented_budget()


def test_exact_zero_energy_is_distinct_from_lost_positive_energy():
    source = window(reservoir=reservoir(desorption_energy_ev=0., internal_energy_ev=0.))
    for name in ("binding_energy_j", "binding_energy_j_m2", "internal_energy_j", "internal_energy_j_m2"):
        record = source.manifest()["quantities"][name]
        assert record["value"] == 0
        assert record["log_value"] is None
        assert record["exact_zero"] is True
        assert record["qualified"] is True
    assert source.represented_budget()["values"]["total_wall_energy_debit_j"] > 0


def test_records_and_detached_outputs_cannot_change_the_source_or_accept_a_host_budget():
    source = window(status="sourced", provenance="Explicit supplied parameter declaration for a numerical fixture")
    manifest = source.manifest()
    original = json.loads(json.dumps(manifest))
    manifest["quantities"]["released_number"]["value"] = -1.
    manifest["reservoir"]["reference_coverage_m2"] = -1.
    manifest["time_window"].clear()
    manifest["status"] = "manufactured"
    assert source.manifest() == original
    budget = source.represented_budget()
    budget["values"]["released_number"] = -1.
    budget["host_budget_accepted"] = True
    assert source.represented_budget()["values"]["released_number"] > 0
    assert source.represented_budget()["host_budget_accepted"] is False
    with pytest.raises((FrozenInstanceError, AttributeError)):
        source.patch_area_m2 = 5.


def test_quantity_units_are_explicit_and_detached_in_manifest_and_budget():
    source = window()
    manifest, budget = source.manifest(), source.represented_budget()
    units = manifest["quantity_units"]
    assert set(units) == set(manifest["quantities"])
    assert units == budget["quantity_units"]
    for name in ("initial_coverage_m2", "remaining_start_m2", "remaining_stop_m2", "released_m2"):
        assert units[name] == "parent m^-2"
    for name in ("initial_number", "remaining_start_number", "remaining_stop_number", "released_number"):
        assert units[name] == "parent"
    for stem in ("translational_energy", "internal_energy", "binding_energy", "total_wall_energy_debit"):
        assert units[stem + "_j_m2"] == "J m^-2"
        assert units[stem + "_j"] == "J"
    units.clear()
    assert source.manifest()["quantity_units"] == budget["quantity_units"]


def test_energy_residuals_report_signed_represented_ledgers_or_unavailable_overflow():
    source = window()
    manifest, budget = source.manifest(), source.represented_budget()
    for suffix in ("j_m2", "j"):
        values = [manifest["quantities"][name + "_" + suffix]["value"] for name in
                  ("total_wall_energy_debit", "translational_energy", "internal_energy", "binding_energy")]
        expected = math.fsum((values[0], -values[1], -values[2], -values[3]))
        assert manifest["energy_balance_residual_" + suffix] == expected
        assert budget["energy_balance_residual_" + suffix] == expected
    pool = reservoir(attempt_frequency_s=1., reference_coverage_m2=1e300, max_coverage_m2=1e300)
    overflow = window(reservoir=pool, initial_coverage_m2=1e300, temperature_k=1e300,
                      parent_duration_s=1., start_s=0., stop_s=1., patch_area_m2=1.).manifest()
    assert overflow["energy_balance_residual_j"] is None
    assert overflow["energy_balance_residual_j_m2"] is None


@pytest.mark.parametrize("cost_case", ["translation", "supplied_costs"])
def test_positive_energy_rescue_does_not_depend_on_extended_exponent_range(monkeypatch, cost_case):
    # Some platforms give longdouble the same exponent range as float64.
    # Positive per-particle costs can underflow there while their source-total
    # energy remains representable; compute the independent total in Decimal.
    monkeypatch.setattr(source_window_module, "LD", np.float64)
    tiny = float(np.nextafter(0., 1.))
    pool = reservoir(attempt_frequency_s=1., reference_coverage_m2=1e300, max_coverage_m2=1e300,
                     desorption_energy_ev=0. if cost_case == "translation" else tiny,
                     internal_energy_ev=0. if cost_case == "translation" else tiny)
    source = window(reservoir=pool, initial_coverage_m2=1e300,
                    temperature_k=tiny if cost_case == "translation" else 700.,
                    parent_duration_s=1., start_s=0., stop_s=1., patch_area_m2=1.)
    expected = constant_hazard_reference(source)
    manifest = source.manifest()
    names = (("translational_energy_j", "translational_energy_j_m2") if cost_case == "translation" else
             ("binding_energy_j", "binding_energy_j_m2", "internal_energy_j", "internal_energy_j_m2"))
    for name in names:
        record = manifest["quantities"][name]
        assert expected[name] > 0
        assert record["value"] == pytest.approx(expected[name], rel=5e-12, abs=0.)
        assert record["exact_zero"] is False
        assert record["qualified"] is True
    assert source.represented_budget()["values"][names[0]] > 0


@pytest.mark.parametrize("field", ["initial_coverage_m2", "temperature_k", "parent_duration_s", "patch_area_m2", "log_mass_tolerance"])
@pytest.mark.parametrize("value", [0., -1., True, math.inf, math.nan, np.ma.masked])
def test_positive_constructor_inputs_are_explicit_finite_unmasked_numbers(field, value):
    with pytest.raises(ValueError):
        window(**{field: value})


@pytest.mark.parametrize("changes", [dict(start_s=-.1), dict(start_s=.9, stop_s=.3), dict(start_s=.3, stop_s=.3),
                                     dict(stop_s=3.), dict(start_s=True), dict(stop_s=math.nan),
                                     dict(start_s=np.ma.masked), dict(reservoir=None), dict(status="unsupported"),
                                     dict(provenance=""), dict(initial_coverage_m2=1e25)])
def test_invalid_window_and_association_are_rejected(changes):
    with pytest.raises(ValueError):
        window(**changes)
