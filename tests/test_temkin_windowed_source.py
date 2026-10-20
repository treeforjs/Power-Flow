"""Original-interval conditioning, local physics and source association."""
from copy import copy
from collections.abc import Sequence
from dataclasses import FrozenInstanceError, replace
import json
import math

import numpy as np
import pytest

from mhdlab.constants import KB
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_release import TemkinReservoir
from mhdlab.temkin_windowed_source import (
    SCHEMA, SKIP_REASON, SourceWindowSelection, TemkinPlanarSourceWindows,
    TemkinPlanarWindowComponent, condition_temkin_planar_source,
)


def pool(name="pool-a", state="M:v0", **changes):
    parameters = dict(name=name, state_id=state, attempt_frequency_s=.4,
                      activation_energy_ev=0., desorption_energy_ev=.2,
                      internal_energy_ev=.05, coverage_coefficient=0.,
                      reference_coverage_m2=1e9, max_coverage_m2=1e9,
                      kinetic_provenance="Manufactured independent constant hazard",
                      energy_provenance="Manufactured costs separate from kinetics")
    parameters.update(changes)
    return TemkinReservoir(**parameters)


def history(*, origin=0., initial=None, times=None, x=None, z=None, reservoirs=None):
    times = [origin, origin+1., origin+3.] if times is None else times
    x = [-2., -1., 2.] if x is None else x
    z = [0., 2., 5.] if z is None else z
    if reservoirs is None:
        reservoirs = (pool(), pool("pool-b", "M:v1", attempt_frequency_s=.9))
    shape = (len(times), len(z)-1, len(x)-1)
    temperatures = np.arange(np.prod(shape)).reshape(shape)*20.+400.
    base = FlashSurfaceHistory([0., 1.], [0., 1.], [[400.], [400.]], "synthetic",
                               "windowed-source-test", "Manufactured local patch temperatures",
                               "Manufactured stationary boundary")
    attributes = base.attributes
    attributes["window_test_evidence"] = {"original": ["retained", {"value": 7}]}
    surface = PlanarSurfaceHistory(times, x, z, temperatures, attributes)
    if initial is None:
        initial = np.array([[[1e8, 0.], [2e8, 3e8]], [[0., 4e8], [5e8, 6e8]]])
    return evolve_temkin_planar_release(surface, initial, reservoirs,
                                        parameter_provenance="Original manufactured source parameters")


def condition(sample, selections=None, **changes):
    values = dict(log_mass_tolerance=1e-7, status="manufactured",
                  provenance="Explicit manufactured original-component source restriction")
    values.update(changes)
    return condition_temkin_planar_source(sample, [SourceWindowSelection(0, .1, .7)] if selections is None else selections, **values)


def test_selection_order_original_interval_initials_and_local_temperatures_are_preserved():
    sample = history()
    choices = [SourceWindowSelection(1, .5, 1.5), SourceWindowSelection(0, 0., .4), SourceWindowSelection(0, .4, 1.)]
    result = condition(sample, choices)
    assert isinstance(result, TemkinPlanarSourceWindows)
    assert result.selections == tuple(choices)
    expected = [(s, r, z, x) for s, selection in enumerate(choices) for r in range(2)
                for z in range(2) for x in range(2) if sample.coverage_m2[selection.interval_index, r, z, x] > 0]
    assert [(c.selection_index, c.reservoir_index, c.current_index, c.width_index) for c in result.components] == expected
    for component in result.components:
        s, r, z, x = component.selection_index, component.reservoir_index, component.current_index, component.width_index
        selection = choices[s]
        t = selection.interval_index
        window = component.source_window
        assert component.interval_index == t
        assert component.interval_origin_s == sample.time_edges_s[t]
        assert window.initial_coverage_m2 == sample.coverage_m2[t, r, z, x]
        assert window.temperature_k == sample.temperature_k[t, z, x]
        assert window.reservoir is sample.reservoirs[r]
        assert window.reservoir.reference_coverage_m2 == 1e9
        assert window.parent_duration_s == np.diff(sample.time_edges_s)[t]
        assert (window.start_s, window.stop_s) == (selection.start_s, selection.stop_s)
        assert component.x_bounds_m == tuple(sample.x_edges_m[x:x+2])
        assert component.z_bounds_m == tuple(sample.z_edges_m[z:z+2])
        assert window.patch_area_m2 == sample.patch_areas_m2[z, x]
        expected_density = window.initial_coverage_m2*math.exp(-window.reservoir.attempt_frequency_s*selection.start_s)*(-math.expm1(-window.reservoir.attempt_frequency_s*(selection.stop_s-selection.start_s)))
        quantities = window.manifest()["quantities"]
        assert quantities["released_m2"]["value"] == pytest.approx(expected_density, rel=3e-14)
        assert quantities["released_number"]["value"] == pytest.approx(expected_density*sample.patch_areas_m2[z, x], rel=3e-14)
        assert quantities["translational_energy_j"]["value"] == pytest.approx(expected_density*sample.patch_areas_m2[z, x]*2*KB*sample.temperature_k[t, z, x], rel=4e-14)


def test_adjacent_disjoint_windows_partition_original_interval_without_restarting_inventory():
    sample = history()
    result = condition(sample, [SourceWindowSelection(0, .5, 1.), SourceWindowSelection(0, 0., .5)])
    for r, z, x in np.ndindex(2, 2, 2):
        components = [c for c in result.components if (c.reservoir_index, c.current_index, c.width_index) == (r, z, x)]
        actual = math.fsum(c.source_window.manifest()["quantities"]["released_number"]["value"] for c in components)
        expected = sample.released_m2[0, r, z, x]*sample.patch_areas_m2[z, x]
        assert actual == pytest.approx(expected, rel=3e-14)


def test_interval_origin_is_not_added_to_local_window_or_quantile_coordinates():
    sample = history(origin=1e12)
    selection = SourceWindowSelection(1, 1e-5, 1e-4)
    result = condition(sample, [selection])
    component = result.components[0]
    assert component.interval_origin_s == 1e12+1
    assert component.source_window.start_s == 1e-5
    delays = component.source_window.time_law.quantile([.1, .5, .9])
    assert np.all(np.diff(delays) > 0)
    assert np.all(delays < 9e-5)
    assert len(set(component.interval_origin_s+selection.start_s+delays)) < 3


def test_zero_initials_are_explicitly_skipped_but_never_inferred_from_release():
    sample = history()
    result = condition(sample)
    assert len(result.components) == 6
    assert len(result.skipped_zero_initials) == 2
    assert {(row["reservoir_index"], row["current_index"], row["width_index"]) for row in result.skipped_zero_initials} == {(0, 0, 1), (1, 0, 0)}
    assert all(row["reason"] == SKIP_REASON for row in result.skipped_zero_initials)
    empty = condition(history(initial=np.zeros((2, 2, 2))))
    assert empty.components == ()
    assert len(empty.skipped_zero_initials) == 8
    assert empty.manifest()["host_budget_accepted"] is False


def test_positive_initials_with_underflowed_density_retain_qualified_extensive_source_fields():
    reservoir = pool(reference_coverage_m2=1e-100, max_coverage_m2=1e-100, attempt_frequency_s=1.,
                     desorption_energy_ev=0., internal_energy_ev=0.)
    sample = history(times=[0., 1000.], x=[0., 1e150], z=[0., 1e150],
                     reservoirs=[reservoir], initial=np.array([[[1e-100]]]))
    result = condition(sample, [SourceWindowSelection(0, 800., 801.)])
    assert len(result.components) == 1 and not result.skipped_zero_initials
    window = result.components[0].source_window
    quantities = window.manifest()["quantities"]
    assert quantities["released_m2"]["value"] == 0
    assert quantities["released_m2"]["qualified"] is False
    assert quantities["released_number"]["qualified"] is True
    expected_log = math.log(1e-100)+math.log(1e300)-800+math.log(-math.expm1(-1.))
    assert quantities["released_number"]["log_value"] == pytest.approx(expected_log, abs=3e-12)
    assert quantities["released_number"]["value"] > 0
    assert quantities["translational_energy_j"]["qualified"] is True
    assert window.all_budget_values_qualified is False


def test_original_parameter_provenance_and_nested_metadata_are_immutable_and_manifest_is_detached():
    sample = history()
    result = condition(sample, status="unvalidated", provenance="Owner supplied restriction, no inferred source status")
    assert result.parameter_provenance == sample.parameter_provenance
    assert result.status == "unvalidated"
    assert result.components[0].source_window.status == "unvalidated"
    assert result.components[0].source_window.provenance == result.provenance
    assert result.surface_attributes["window_test_evidence"]["original"][1]["value"] == 7
    with pytest.raises(TypeError):
        result.surface_attributes["window_test_evidence"]["original"][1]["value"] = 9
    with pytest.raises(TypeError):
        result.skipped_zero_initials[0]["reason"] = "changed"
    for name in ("x_edges_m", "z_edges_m", "time_edges_s"):
        with pytest.raises(ValueError):
            getattr(result, name).setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        result.status = "sourced"
    manifest = result.manifest()
    assert manifest["schema"] == SCHEMA
    assert manifest["host_budget_accepted"] is False
    assert len(manifest["source_reservoirs"]) == 2
    json.dumps(manifest, allow_nan=False)
    manifest["surface_attributes"]["window_test_evidence"]["original"][1]["value"] = 99
    manifest["components"][0]["source_window"]["reservoir"]["name"] = "changed"
    manifest["skipped_zero_initials"].clear()
    assert result.manifest()["surface_attributes"]["window_test_evidence"]["original"][1]["value"] == 7
    assert result.components[0].source_window.reservoir.name == "pool-a"
    assert len(result.skipped_zero_initials) == 2


def test_conditioning_never_calls_the_source_evolution_driver_and_does_not_mutate_input(monkeypatch):
    import mhdlab.temkin_planar as source_module
    sample = history()
    before = sample.coverage_m2.tobytes(), sample.released_m2.tobytes()
    monkeypatch.setattr(source_module, "evolve_temkin_planar_release", lambda *a, **k: pytest.fail("new source evolution is forbidden"))
    result = condition(sample)
    assert result.components
    assert before == (sample.coverage_m2.tobytes(), sample.released_m2.tobytes())


@pytest.mark.parametrize("index,start,stop", [
    (True, 0., 1.), (1., 0., 1.), (-1, 0., 1.), (np.ma.array(0), 0., 1.),
    (0, True, 1.), (0, 0., False), (0, float("nan"), 1.), (0, 0., float("inf")),
    (0, -.1, 1.), (0, 1., 1.), (0, 2., 1.), (0, np.ma.array(0.), 1.),
])
def test_selection_rejects_invalid_indices_offsets_masks_and_empty_ranges(index, start, stop):
    with pytest.raises(ValueError):
        SourceWindowSelection(index, start, stop)


@pytest.mark.parametrize("selections", [
    [], [None], [dict(interval_index=0, start_s=0, stop_s=1)], np.ma.array([0]),
    [SourceWindowSelection(2, 0., .1)], [SourceWindowSelection(0, 0., 1.1)],
    [SourceWindowSelection(0, 0., .6), SourceWindowSelection(0, .5, 1.)],
    [SourceWindowSelection(0, .2, .4), SourceWindowSelection(0, .2, .4)],
    [SourceWindowSelection(0, .3, .5), SourceWindowSelection(0, .1, .7)],
])
def test_invalid_or_overlapping_selections_fail_before_any_scalar_windows(monkeypatch, selections):
    import mhdlab.temkin_windowed_source as module
    sample = history()
    monkeypatch.setattr(module, "FrozenTemkinSourceWindow", lambda *a, **k: pytest.fail("invalid selection reached window construction"))
    with pytest.raises(ValueError):
        condition(sample, selections)


@pytest.mark.parametrize("changes", [
    {"log_mass_tolerance": 0.}, {"log_mass_tolerance": -1.}, {"log_mass_tolerance": True},
    {"log_mass_tolerance": float("nan")}, {"log_mass_tolerance": np.ma.array(1e-7)},
    {"status": "invented"}, {"status": []}, {"provenance": ""}, {"provenance": None},
])
def test_window_qualification_labels_are_validated_even_when_every_original_inventory_is_zero(changes):
    sample = history(initial=np.zeros((2, 2, 2)))
    with pytest.raises(ValueError):
        condition(sample, **changes)


@pytest.mark.parametrize("field", ["schema", "law_identity", "coverage_model", "parameter_provenance",
                                   "temperature_k", "released_m2", "total_wall_energy_debit_j_m2", "inversion_regime", "reservoirs"])
def test_forged_critical_canonical_source_fields_are_rejected(field):
    original = history()
    forged = copy(original)
    if field in {"schema", "law_identity", "coverage_model"}:
        value = "unsupported"
    elif field == "parameter_provenance":
        value = ""
    elif field == "reservoirs":
        value = (replace(original.reservoirs[0], attempt_frequency_s=.5), original.reservoirs[1])
    else:
        value = np.array(getattr(original, field), copy=True)
        value.flat[0] += 1 if field != "inversion_regime" else 2
    object.__setattr__(forged, field, value)
    with pytest.raises(ValueError):
        condition(forged)


def test_forged_selection_is_revalidated_and_noncanonical_input_is_rejected():
    sample = history()
    forged = SourceWindowSelection(0, 0., .1)
    object.__setattr__(forged, "interval_index", True)
    with pytest.raises(ValueError):
        condition(sample, [forged])
    with pytest.raises(ValueError):
        condition(object())


@pytest.mark.parametrize("changes", [
    {"components": ()}, {"components": (object(),)}, {"skipped_zero_initials": ()},
    {"parameter_provenance": ""}, {"status": "unknown"}, {"log_mass_tolerance": False},
])
def test_direct_result_reconstruction_rejects_inconsistent_component_or_source_association(changes):
    with pytest.raises(ValueError):
        replace(condition(history()), **changes)


def test_component_and_direct_result_reject_changed_area_origin_or_source_offsets():
    result = condition(history())
    component = result.components[0]
    with pytest.raises(ValueError):
        replace(component, x_bounds_m=(-2., 0.))
    with pytest.raises(ValueError):
        replace(result, components=(replace(component, interval_origin_s=1.), *result.components[1:]))
    with pytest.raises(ValueError):
        replace(result, selections=(SourceWindowSelection(0, .2, .7),))


def test_selected_component_limit_is_checked_before_scalar_conditioning(monkeypatch):
    import mhdlab.temkin_windowed_source as module
    sample = history()
    monkeypatch.setattr(module, "MAX_SELECTED_COMPONENTS", 1)
    with pytest.raises(ValueError, match="selected-component"):
        condition(sample)


@pytest.mark.parametrize("selections", ["not selections", iter([SourceWindowSelection(0, 0., .1)])])
def test_selection_inputs_require_finite_sequences_instead_of_arbitrary_iterators(selections):
    with pytest.raises(ValueError, match="finite Sequence"):
        condition(history(), selections)


def test_selection_length_limit_precedes_any_indexed_expansion():
    import mhdlab.temkin_windowed_source as module
    class Oversized(Sequence):
        def __len__(self):
            return module.MAX_SELECTIONS+1
        def __getitem__(self, index):
            pytest.fail("oversized selection list must not be expanded")
    with pytest.raises(ValueError, match="selection-count"):
        condition(history(), Oversized())


@pytest.mark.parametrize("attribute", ["initial_coverage_m2", "temperature_k"])
def test_direct_result_cannot_replace_a_component_with_another_valid_but_noncanonical_window(attribute):
    result = condition(history())
    component = result.components[0]
    old = component.source_window
    changed = replace(old, **{attribute: getattr(old, attribute)*1.01})
    altered = replace(component, source_window=changed)
    with pytest.raises(ValueError, match="association"):
        replace(result, components=(altered, *result.components[1:]))


def test_direct_result_cannot_relabel_a_positive_component_as_an_exact_zero_skip():
    result = condition(history())
    component = result.components[0]
    skip = {name: getattr(component, name) for name in ("selection_index", "interval_index", "reservoir_index", "current_index", "width_index")}
    skip["reason"] = SKIP_REASON
    skipped = tuple(sorted((*result.skipped_zero_initials, skip), key=lambda row: (row["selection_index"], row["reservoir_index"], row["current_index"], row["width_index"])))
    with pytest.raises(ValueError, match="positive canonical"):
        replace(result, components=result.components[1:], skipped_zero_initials=skipped)


def test_validated_copy_rechecks_the_authoritative_source_and_forged_component_fields():
    result = condition(history())
    checked = result.validated_copy()
    assert checked.manifest() == result.manifest()
    forged_history = copy(result._source_history)
    altered_release = np.array(forged_history.released_m2)
    altered_release.flat[0] += 1.
    object.__setattr__(forged_history, "released_m2", altered_release)
    forged = copy(result)
    object.__setattr__(forged, "_source_history", forged_history)
    with pytest.raises(ValueError):
        forged.validated_copy()
    forged_component = copy(result.components[0])
    object.__setattr__(forged_component, "width_index", True)
    forged = copy(result)
    object.__setattr__(forged, "components", (forged_component, *result.components[1:]))
    with pytest.raises(ValueError):
        forged.validated_copy()


def test_direct_result_rejects_other_source_provenance_and_missing_authoritative_history():
    result = condition(history())
    with pytest.raises(ValueError, match="canonical"):
        replace(result, parameter_provenance="A different supplied source declaration")
    with pytest.raises(ValueError, match="authoritative"):
        replace(result, _source_history=None)


def test_finite_sequence_materialization_uses_bounded_indexing_without_arbitrary_iterators():
    class IndexedOnly(Sequence):
        def __init__(self, values):
            self.values = values
        def __len__(self):
            return len(self.values)
        def __getitem__(self, index):
            return self.values[index]
        def __iter__(self):
            pytest.fail("a custom iterator must not control finite materialization")
    result = condition(history(), IndexedOnly([SourceWindowSelection(0, .1, .7)]))
    checked = replace(result, components=IndexedOnly(result.components), source_reservoirs=IndexedOnly(result.source_reservoirs),
                      skipped_zero_initials=IndexedOnly(result.skipped_zero_initials))
    assert checked.manifest() == result.manifest()


@pytest.mark.parametrize("name,value", [("reference_coverage_m2", 0.), ("attempt_frequency_s", -1.), ("kinetic_provenance", "")])
def test_forged_nested_pool_declarations_are_rejected_before_source_replay(name, value):
    sample = history()
    forged_pool = copy(sample.reservoirs[0])
    object.__setattr__(forged_pool, name, value)
    forged = copy(sample)
    object.__setattr__(forged, "reservoirs", (forged_pool, sample.reservoirs[1]))
    with pytest.raises(ValueError):
        condition(forged)
