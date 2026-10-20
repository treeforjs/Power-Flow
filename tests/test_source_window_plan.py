"""Independent exact-coordinate source-window planning checks."""
from collections.abc import Sequence
from copy import copy
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
import json
import math

import numpy as np
import pytest

from mhdlab.source_window_plan import (
    MAX_TIME_EDGES, MAX_WINDOWS, QUALIFICATION_IDENTITY, SCHEMA,
    SourceTimeWindow, SourceWindowPlan, SourceWindowPlanRangeError, plan_source_windows,
)
from mhdlab.temkin_windowed_source import SourceWindowSelection


def q(value):
    return Fraction.from_float(float(value))


def unpack(value):
    return Fraction(int(value["numerator"]), int(value["denominator"]))


def plan(edges=(0., 1., 2., 4.), windows=None, tolerance=1e-12):
    if windows is None:
        windows = [SourceTimeWindow(.25, 3.5)]
    return plan_source_windows(edges, windows, relative_time_tolerance=tolerance)


def test_exact_intersections_preserve_input_window_order_and_local_interval_association():
    windows = [SourceTimeWindow(2.5, 3.5), SourceTimeWindow(.25, 1.5)]
    result = plan(windows=windows)
    assert result.window_indices == (0, 1, 1)
    assert result.selections == (SourceWindowSelection(2, .5, 1.5),
                                 SourceWindowSelection(0, .25, 1.),
                                 SourceWindowSelection(1, 0., .5))
    assert result.windows == tuple(windows)
    report = result.manifest()
    assert report["schema"] == SCHEMA
    assert report["qualification_identity"] == QUALIFICATION_IDENTITY
    assert report["source_emission_windows_are_photon_recording_gates"] is False
    assert report["physical_timing_uncertainty_estimated"] is False
    for record, selection in zip(report["selection_records"], result.selections):
        origin = q(result.time_edges_s[selection.interval_index])
        assert unpack(record["exact_local_start_s"])+origin == unpack(record["exact_global_start_s"])
        assert unpack(record["exact_local_stop_s"])+origin == unpack(record["exact_global_stop_s"])
        assert unpack(record["exact_intersection_duration_s"]) > 0
        assert record["selection"] == selection.manifest()
    for record, window in zip(report["window_records"], windows):
        assert unpack(record["exact_window_duration_s"]) == q(window.stop_s)-q(window.start_s)


def test_binary64_decimal_adjacent_windows_tile_without_zero_boundary_pieces():
    result = plan((0., .1, .3), [SourceTimeWindow(.05, .2), SourceTimeWindow(.2, .3)])
    assert len(result.selections) == 3
    assert result.window_indices == (0, 0, 1)
    assert tuple(value.interval_index for value in result.selections) == (0, 1, 1)
    assert result.selections[1].stop_s == result.selections[2].start_s
    report = result.manifest()
    for index, window in enumerate(result.windows):
        total = sum((unpack(item["exact_intersection_duration_s"]) for item in report["selection_records"]
                     if item["window_index"] == index), Fraction())
        assert total == q(window.stop_s)-q(window.start_s)
    with pytest.raises(ValueError, match="overlap"):
        plan((0., .1, .3), [SourceTimeWindow(.05, .2), SourceTimeWindow(np.nextafter(.2, 0.), .3)])


def test_full_and_edge_aligned_windows_do_not_create_empty_intersections():
    full = plan(windows=[SourceTimeWindow(0., 4.)])
    assert full.selections == (SourceWindowSelection(0, 0., 1.), SourceWindowSelection(1, 0., 1.),
                               SourceWindowSelection(2, 0., 2.))
    aligned = plan(windows=[SourceTimeWindow(1., 2.)])
    assert aligned.selections == (SourceWindowSelection(1, 0., 1.),)


@pytest.mark.parametrize("origin,parent_stop,window_start,window_stop,sign", [
    (1., 1e16, 1e16-4., 1e16, 1),
    (3., float(2**54), float(2**54-4), float(2**54), -1),
])
def test_signed_projection_errors_qualified_against_small_intersection_not_global_time(origin, parent_stop, window_start, window_stop, sign):
    windows = [SourceTimeWindow(window_start, window_stop)]
    with pytest.raises(SourceWindowPlanRangeError, match="relative_time_tolerance"):
        plan((origin, parent_stop), windows, tolerance=.24)
    result = plan((origin, parent_stop), windows, tolerance=.25)
    record = result.manifest()["selection_records"][0]
    assert unpack(record["exact_intersection_duration_s"]) == 4
    for name in ("local_start_projection_error", "local_stop_projection_error", "canonical_parent_duration_error"):
        assert unpack(record[name]["exact_s"]) == sign
        assert record[name]["relative_value"] == sign*.25
    assert record["touches_parent_stop"]
    assert record["parent_duration_error_qualified_against_intersection"]
    assert record["precision_flags"]["local_start_rounded"]
    assert record["precision_flags"]["global_start_addition_rounded"]
    assert record["global_start_projection_error"]["value_s"] == 0.


@pytest.mark.parametrize("edges,window", [
    ((1., 2e16), SourceTimeWindow(1e16, np.nextafter(1e16, math.inf))),
    ((3., float(2**54)), SourceTimeWindow(np.nextafter(float(2**54), -math.inf), float(2**54))),
])
def test_distinct_global_endpoints_that_collapse_locally_are_rejected(edges, window):
    assert window.stop_s > window.start_s
    with pytest.raises(SourceWindowPlanRangeError, match="collapses"):
        plan(edges, [window], tolerance=1.)


def test_irrelevant_parent_stop_rounding_does_not_reject_exact_interior_window():
    report = plan((1., 1e16), [SourceTimeWindow(2., 3.)], tolerance=1e-30).manifest()
    record = report["selection_records"][0]
    assert record["selection"] == {"interval_index": 0, "start_s": 1., "stop_s": 2.}
    assert unpack(record["canonical_parent_duration_error"]["exact_s"]) == 1
    assert "relative_value" not in record["canonical_parent_duration_error"]
    assert not record["parent_duration_error_qualified_against_intersection"]


def test_subnormal_interior_window_keeps_true_positive_support():
    tiny = float(np.nextafter(0., 1.))
    result = plan((tiny, 1.), [SourceTimeWindow(tiny, 2*tiny)], tolerance=1e-30)
    assert result.selections == (SourceWindowSelection(0, 0., tiny),)
    record = result.manifest()["selection_records"][0]
    assert unpack(record["exact_intersection_duration_s"]) == q(tiny)
    assert unpack(record["canonical_parent_duration_error"]["exact_s"]) == q(tiny)
    assert not record["touches_parent_stop"]


def test_unrepresentably_small_relative_error_retains_exact_nonzero_rational_and_flag():
    result = plan((0., 1e-300, 1e300), [SourceTimeWindow(0., 1e300)])
    record = result.manifest()["selection_records"][1]
    error = record["canonical_parent_duration_error"]
    assert unpack(error["exact_s"]) == q(1e-300)
    assert unpack(error["relative_to_intersection"]) > 0
    assert error["relative_value"] == 0.
    assert error["relative_value_underflowed"]
    assert record["precision_flags"]["canonical_parent_duration_rounded"]


def test_large_global_origin_keeps_source_offsets_independent_of_global_summation():
    origin = 1e12
    ulp = np.nextafter(origin, math.inf)-origin
    result = plan((origin, origin+16*ulp), [SourceTimeWindow(origin+3*ulp, origin+5*ulp)], tolerance=1e-30)
    assert result.selections == (SourceWindowSelection(0, 3*ulp, 5*ulp),)
    record = result.manifest()["selection_records"][0]
    assert record["interval_origin_s"] == origin
    assert unpack(record["exact_intersection_duration_s"]) == q(2*ulp)


def test_unrepresentable_global_diagnostic_projection_rejects_explicitly():
    origin, stop = float(3*2**970), float(np.finfo(float).max)
    with pytest.raises(SourceWindowPlanRangeError, match="global diagnostic projection"):
        plan((origin, stop), [SourceTimeWindow(origin, stop)])


def test_qualified_endpoint_motion_can_close_small_global_gap_and_is_reported():
    edge = float(2**54)
    result = plan((3., float(2**55)), [SourceTimeWindow(3., edge-2.),
                                      SourceTimeWindow(edge, float(2**55))])
    assert result.windows[1].start_s-result.windows[0].stop_s == 2.
    assert result.selections[0].stop_s == result.selections[1].start_s
    report = result.manifest()
    assert "close exact gaps" in report["projection_limitations"]
    first, second = report["selection_records"]
    assert unpack(first["local_stop_projection_error"]["exact_s"]) == 1
    assert unpack(second["local_start_projection_error"]["exact_s"]) == -1
    assert first["precision_flags"]["local_stop_rounded"]
    assert second["precision_flags"]["local_start_rounded"]


def test_plan_and_nested_records_are_immutable_and_manifest_is_detached():
    result = plan()
    assert isinstance(result.time_edges_s, tuple) and isinstance(result.windows, tuple)
    with pytest.raises(FrozenInstanceError):
        result.selections = ()
    with pytest.raises(TypeError):
        result.window_indices[0] = 9
    with pytest.raises(FrozenInstanceError):
        result.windows[0].start_s = 0.
    report = result.manifest()
    json.dumps(report, allow_nan=False)
    report["selection_records"][0]["exact_local_start_s"]["numerator"] = "wrong"
    report["selections"][0]["interval_index"] = 99
    assert result.manifest()["selection_records"][0]["exact_local_start_s"]["numerator"] != "wrong"
    assert result.selections[0].interval_index == 0
    assert result.validated_copy().manifest() == result.manifest()


@pytest.mark.parametrize("field,value", [("selections", ()), ("window_indices", (True, 0, 0)),
                                        ("_manifest_json", "{}"), ("relative_time_tolerance", .5),
                                        ("time_edges_s", (0., 1., 2., 5.))])
def test_forged_derived_plan_fields_rejected_by_validated_copy_and_manifest(field, value):
    result = copy(plan())
    object.__setattr__(result, field, value)
    with pytest.raises(ValueError, match="constructor replay"):
        result.validated_copy()
    with pytest.raises(ValueError, match="constructor replay"):
        result.manifest()


def test_forged_nested_selection_boolean_index_is_not_equal_to_valid_integer():
    result = copy(plan())
    selection = copy(result.selections[1])
    object.__setattr__(selection, "interval_index", True)
    object.__setattr__(result, "selections", (result.selections[0], selection, result.selections[2]))
    with pytest.raises(ValueError, match="constructor replay"):
        result.validated_copy()


@pytest.mark.parametrize("start,stop", [(True, 1.), (0., False), (np.ma.array(0.), 1.),
                                      (0., math.inf), (math.nan, 1.), (-1., 1.), (1., 1.), (2., 1.), ("0", 1.)])
def test_invalid_global_window_scalars_reject(start, stop):
    with pytest.raises(ValueError):
        SourceTimeWindow(start, stop)


@pytest.mark.parametrize("edges", [[], [0.], [0., 0.], [1., 0.], [-1., 1.], [0., math.inf],
                                  [0., math.nan], [False, 1.], [[0., 1.]], np.ma.array([0., 1.], mask=False),
                                  (x for x in [0., 1.]), "0,1"])
def test_invalid_time_edge_sequences_reject(edges):
    with pytest.raises(ValueError):
        plan(edges, [SourceTimeWindow(.25, .5)])


@pytest.mark.parametrize("windows", [[], [None], [(.25, .5)], True, "windows", np.ma.array([1.]),
                                    (SourceTimeWindow(.25, .5) for _ in range(1))])
def test_invalid_window_collections_reject(windows):
    with pytest.raises(ValueError):
        plan(windows=windows)


@pytest.mark.parametrize("windows", [[SourceTimeWindow(0., 1.), SourceTimeWindow(0., 1.)],
                                    [SourceTimeWindow(.5, 2.), SourceTimeWindow(1., 3.)],
                                    [SourceTimeWindow(.5, 3.), SourceTimeWindow(1., 2.)]])
def test_overlap_duplicate_or_nesting_reject_without_repair(windows):
    with pytest.raises(ValueError, match="overlap"):
        plan(windows=windows)


@pytest.mark.parametrize("window", [SourceTimeWindow(0., 2.), SourceTimeWindow(2., 5.), SourceTimeWindow(0., 5.)])
def test_out_of_support_windows_reject_instead_of_truncation(window):
    with pytest.raises(ValueError, match="fully supported"):
        plan((1., 2., 4.), [window])


@pytest.mark.parametrize("tolerance", [True, 0., -1., math.inf, math.nan, None, np.ma.array(1e-3)])
def test_invalid_explicit_relative_tolerance_rejects(tolerance):
    with pytest.raises(ValueError, match="relative_time_tolerance"):
        plan(tolerance=tolerance)


class IndexedOnly(Sequence):
    def __init__(self, values):
        self.values = values
    def __len__(self):
        return len(self.values)
    def __getitem__(self, index):
        return self.values[index]
    def __iter__(self):
        raise AssertionError("bounded indexed reads are required")


class Oversized(Sequence):
    def __init__(self, count):
        self.count = count
    def __len__(self):
        return self.count
    def __getitem__(self, index):
        raise AssertionError("length cap must precede indexed reads")


def test_finite_indexed_sequences_and_numpy_edges_supported_without_arbitrary_iteration():
    one = plan(IndexedOnly([0., 1.]), IndexedOnly([SourceTimeWindow(.25, .75)]))
    two = plan(np.array([0., 1.]), [SourceTimeWindow(.25, .75)])
    assert one.manifest() == two.manifest()


def test_input_count_limits_precede_materialization():
    with pytest.raises(ValueError, match="work limit"):
        plan(Oversized(MAX_TIME_EDGES+1), [SourceTimeWindow(0., 1.)])
    with pytest.raises(ValueError, match="work limit"):
        plan(windows=Oversized(MAX_WINDOWS+1))


def test_intersection_expansion_is_capped_before_selection_construction(monkeypatch):
    import mhdlab.source_window_plan as module
    monkeypatch.setattr(module, "MAX_SELECTIONS", 1)
    monkeypatch.setattr(module, "SourceWindowSelection", lambda *a: pytest.fail("expansion cap must precede construction"))
    with pytest.raises(ValueError, match="intersection count"):
        module.plan_source_windows((0., 1., 2.), [SourceTimeWindow(0., 2.)], relative_time_tolerance=1e-12)


def test_public_plan_constructor_rebuilds_immutable_inputs_and_cannot_accept_derived_overrides():
    result = SourceWindowPlan([0., 1.], [SourceTimeWindow(.25, .5)], 1e-12)
    assert result.time_edges_s == (0., 1.)
    with pytest.raises(TypeError):
        SourceWindowPlan([0., 1.], [SourceTimeWindow(.25, .5)], 1e-12, selections=())
    malformed = copy(SourceTimeWindow(.25, .5))
    object.__setattr__(malformed, "start_s", True)
    with pytest.raises(ValueError):
        plan((0., 1.), [malformed])
