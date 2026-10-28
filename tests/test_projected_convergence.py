import pytest

from mhdlab.projected_convergence import AXES, qualify_projected_convergence


def request():
    hypotheses = {}
    for name, signal in (("localized", 100.), ("uniform", 90.)):
        counts = {key: 32 for key in AXES}
        case = {"comparison_contract": "same source, clock, optics, calibration and excitation",
                "finest": {"run_id": name + "_fine", "provenance": "manufactured",
                           "resolution": counts, "values": {"line": signal}},
                "coarser_by_axis": {}}
        for key in AXES:
            coarse_counts = dict(counts, **{key: 16})
            case["coarser_by_axis"][key] = {"run_id": name + "_coarse_" + key,
                "provenance": "manufactured", "resolution": coarse_counts, "values": {"line": signal + .01}}
        hypotheses[name] = case
    return {"schema": "mhdlab-projected-convergence-v1", "status": "manufactured",
        "comparison_contract": "same source, clock, optics, calibration and excitation",
        "features": [{"id": "line", "units": "W m^-2 sr^-1"}],
        "claim_feature_ids": ["line"], "hypotheses": hypotheses}


def with_bounds(value):
    for row in value["hypotheses"].values():
        row["numerical_uncertainty"] = {"bounds": {"line": .1}, "provenance": "manufactured bound",
                                       "method": "synthetic exact-solution error"}
    return value


def test_small_resolution_changes_do_not_substitute_for_error_bounds():
    report = qualify_projected_convergence(request())
    assert report["resolution_criteria_met"] is True
    assert report["contrast_criteria_met"] is False
    assert report["contrasts"]["line"]["combined_owner_uncertainty_bound"] is None


def test_explicit_bound_and_fixed_independent_axes_qualify_synthetic_contrast():
    report = qualify_projected_convergence(with_bounds(request()))
    assert report["contrast_criteria_met"] is True
    assert report["physical_discrimination_established"] is False


def test_zero_contrast_does_not_establish_discrimination_even_with_zero_uncertainty():
    value = with_bounds(request())
    for row in value["hypotheses"].values():
        row["finest"]["values"]["line"] = 100.
        row["numerical_uncertainty"]["bounds"]["line"] = 0.
        for coarse in row["coarser_by_axis"].values():
            coarse["values"]["line"] = 100.
    report = qualify_projected_convergence(value)
    assert report["resolution_criteria_met"] is True
    assert report["contrast_criteria_met"] is False


def test_missing_state_refinement_remains_incomplete():
    value = with_bounds(request())
    del value["hypotheses"]["localized"]["coarser_by_axis"]["state"]
    report = qualify_projected_convergence(value)
    assert report["hypotheses"]["localized"]["missing_axes"] == ["state"]
    assert report["contrast_criteria_met"] is False


def test_opposite_changes_cannot_cancel_in_resolution_budget():
    value = request()
    case = value["hypotheses"]["localized"]
    case["coarser_by_axis"]["mesh"]["values"]["line"] = 100.5
    case["coarser_by_axis"]["time"]["values"]["line"] = 99.5
    report = qualify_projected_convergence(value)
    assert report["hypotheses"]["localized"]["sum_absolute_resolution_changes"]["line"] == pytest.approx(1.02)


def test_two_axes_cannot_be_refined_in_one_independence_check():
    value = request()
    value["hypotheses"]["localized"]["coarser_by_axis"]["mesh"]["resolution"]["time"] = 16
    with pytest.raises(ValueError, match="exactly the declared axis"):
        qualify_projected_convergence(value)


def test_zero_feature_has_undefined_relative_change_and_cannot_hide_with_scale_floor():
    value = request()
    case = value["hypotheses"]["localized"]
    case["finest"]["values"]["line"] = 0.
    report = qualify_projected_convergence(value)
    assert report["hypotheses"]["localized"]["axes"]["mesh"]["relative_change"]["line"] is None
    value["features"][0]["relative_scale_floor"] = 1e6
    with pytest.raises(ValueError, match="arbitrary normalization floor"):
        qualify_projected_convergence(value)


def test_report_preserves_run_identity_units_and_uncertainty_method():
    value = with_bounds(request())
    report = qualify_projected_convergence(value)
    assert report["input_evidence"] == value
    value["hypotheses"]["localized"]["finest"]["values"]["line"] = -100
    assert report["input_evidence"]["hypotheses"]["localized"]["finest"]["values"]["line"] == 100


def test_run_identity_cannot_describe_conflicting_resolutions():
    value = request()
    rows = value["hypotheses"]["localized"]["coarser_by_axis"]
    rows["time"]["run_id"] = rows["mesh"]["run_id"]
    with pytest.raises(ValueError, match="distinct refinement evidence"):
        qualify_projected_convergence(value)


def test_owner_uncertainty_method_cannot_be_placeholder_evidence():
    value = with_bounds(request())
    value["status"] = "owner_supplied"
    value["hypotheses"]["localized"]["numerical_uncertainty"]["method"] = "pending"
    with pytest.raises(ValueError, match="placeholder evidence"):
        qualify_projected_convergence(value)


def test_exact_one_percent_and_ten_percent_do_not_meet_strict_thresholds():
    value = with_bounds(request())
    value["hypotheses"]["localized"]["coarser_by_axis"]["state"]["values"]["line"] = 101.
    assert qualify_projected_convergence(value)["resolution_criteria_met"] is False
    value = with_bounds(request())
    for case in value["hypotheses"].values():
        case["numerical_uncertainty"]["bounds"]["line"] = .5
    assert qualify_projected_convergence(value)["contrast_criteria_met"] is False


@pytest.mark.parametrize("mutation", [
    lambda v: v["hypotheses"]["localized"].update(comparison_contract="different gate"),
    lambda v: v["features"][0].update(relative_scale_floor=1e100),
    lambda v: v["hypotheses"]["localized"]["finest"]["values"].update(line=True),
    lambda v: v.update(claim_feature_ids=["unknown"]),
    lambda v: v["hypotheses"]["localized"]["numerical_uncertainty"].update(bounds={"line": -.1}),
])
def test_incompatible_or_invalid_evidence_rejects(mutation):
    value = with_bounds(request())
    mutation(value)
    with pytest.raises(ValueError):
        qualify_projected_convergence(value)
