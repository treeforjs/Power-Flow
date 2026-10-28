"""Identifiability and readout semantics, independent of any CR implementation."""

from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest

from mhdlab.neutral_spectral_inference import (
    EnsembleObservation, NeutralForwardCase, ProjectedCurve, ReadoutEnvelope,
    compare_neutral_ensemble, exploratory_calibration_scales,
)


ROOT = Path(__file__).resolve().parents[1]


def fixture():
    config = json.loads((ROOT / "configs/neutral_ensemble_manufactured.json").read_bytes())
    observation = dict(config["observation"])
    observation["curves"] = tuple(ReadoutEnvelope(**curve) for curve in observation["curves"])
    cases = []
    for values in config["cases"]:
        case = dict(values)
        case["curves"] = tuple(ProjectedCurve(**curve) for curve in case["curves"])
        cases.append(NeutralForwardCase(**case))
    return EnsembleObservation(**observation), cases


def row(report, case_id, scale=1.):
    return next(item for item in report["case_scale_comparisons"] if item["case_id"] == case_id and item["C_abs"] == scale)


def test_manufactured_recovery_reports_all_cases_and_separate_scenario_scale():
    observed, cases = fixture()
    report = compare_neutral_ensemble(observed, cases)
    scales = report["calibration"]["scales"]
    np.testing.assert_allclose(scales, 10**np.linspace(-2, 2, 41))
    assert len(report["case_scale_comparisons"]) == 41*len(cases)
    target = row(report, "localized_source_N2_T600")
    assert target["max_normalized_departure"] == 0
    assert target["readout_compatible"]
    assert not row(report, "uniform_source_N2_T600")["readout_compatible"]
    assert report["diagnostics"]["multiple_compatible_inventory_values"]
    assert report["diagnostics"]["compatible_thermal_values_K"] == [600.]
    assert not report["diagnostics"]["populations_rescaled"]
    assert not report["diagnostics"]["confidence_intervals_calculated"]


def test_calibration_multiplies_observation_not_prediction_or_inverse_response():
    observed, cases = fixture()
    original = cases[3]
    values = original.curves[0].values.copy()
    report = compare_neutral_ensemble(observed, [replace(original, curves=(replace(original.curves[0], values=2*values),))], [.5, 1., 2.])
    assert row(report, original.case_id, 2.)["max_normalized_departure"] == 0
    assert not row(report, original.case_id, .5)["readout_compatible"]
    np.testing.assert_array_equal(original.curves[0].values, values)


def test_asymmetric_envelope_uses_half_width_on_prediction_side():
    observed, cases = fixture()
    curve = observed.curves[0]
    # Upper first-point half-width is .2; lower half-width is .1.
    case = cases[3]
    above = replace(case, curves=(replace(case.curves[0], values=[2.1, 4., 2.]),))
    below = replace(case, case_id="below", curves=(replace(case.curves[0], values=[1.9, 4., 2.]),))
    report = compare_neutral_ensemble(observed, [above, below], [1.])
    assert row(report, above.case_id)["max_normalized_departure"] == pytest.approx(.5)
    assert row(report, below.case_id)["max_normalized_departure"] == 1


def test_denser_correlated_digitization_does_not_increase_score():
    observed, cases = fixture()
    envelope, prediction = observed.curves[0], cases[3].curves[0]
    baseline = compare_neutral_ensemble(observed, [replace(cases[3], curves=(replace(prediction, values=[2.1, 4., 2.]),))], [1.])
    # Repeat the worst display sample at a new retained abscissa. Max, not sum.
    dense_observed = replace(observed, curves=(replace(envelope, coordinate=[119., 120., 130., 150.],
                              center=[2., 2., 4., 2.], lower=[1.9, 1.9, 3.8, 1.7],
                              upper=[2.2, 2.2, 4.3, 2.1], support_segment_ids=[0, 0, 0, 1]),))
    dense_case = replace(cases[3], curves=(replace(prediction, coordinate=[119., 120., 130., 150.],
                         values=[2.1, 2.1, 4., 2.], support_segment_ids=[0, 0, 0, 1]),))
    dense = compare_neutral_ensemble(dense_observed, [dense_case], [1.])
    assert row(dense, dense_case.case_id)["max_normalized_departure"] == row(baseline, dense_case.case_id)["max_normalized_departure"]


def test_thermal_and_hypothesis_degeneracies_retain_alternative_cases():
    observed, cases = fixture()
    case = cases[3]
    duplicate = replace(case, case_id="alternative_temperature", thermal_parameter_K=900.)
    null_case = replace(case, case_id="linear_matched_null", hypothesis="uniform_source")
    report = compare_neutral_ensemble(observed, [case, duplicate, null_case], [1.])
    assert report["diagnostics"]["compatible_thermal_values_K"] == [600., 900.]
    assert report["diagnostics"]["multiple_compatible_thermal_values"]
    assert len(report["diagnostics"]["compatible_at_any_scanned_scale_case_ids"]) == 3
    summaries = report["hypothesis_scale_summaries"]
    assert all(summary["minimum_normalized_departure"] == 0 for summary in summaries)
    assert summaries[0]["minimum_case_ids"] == [case.case_id, duplicate.case_id]


def test_independent_inventory_cases_calibration_degeneracy_is_not_hidden():
    observed, cases = fixture()
    cases = [case for case in cases if case.hypothesis == "localized_source" and case.thermal_parameter_K == 600]
    report = compare_neutral_ensemble(observed, cases, [.5, 1.])
    assert row(report, cases[0].case_id, .5)["max_normalized_departure"] == 0
    assert row(report, cases[1].case_id, 1.)["max_normalized_departure"] == 0
    assert report["diagnostics"]["compatible_inventory_values_m2"] == [1e18, 2e18]
    assert len(report["case_metadata"]) == 2
    assert report["diagnostics"]["forward_case_interpolation"] is False


def test_excited_spectrum_alone_cannot_distinguish_supplied_ground_neutral_density_alternatives():
    observed, cases = fixture()
    case = cases[3]
    # The same excited-state spectrum may accompany distinct ground-neutral
    # inventories in complete external cases. The inverse layer adds no closure.
    alternative = replace(case, case_id="same_emission_other_neutral_inventory",
                          parent_water_areal_inventory_m2=4e18,
                          neutral_columns_m2={"H2O[X]": 8e16},
                          weighted_neutral_densities_m3={"H2O[X]": 1.6e19})
    report = compare_neutral_ensemble(observed, [case, alternative], [1.])
    assert report["diagnostics"]["compatible_inventory_values_m2"] == [2e18, 4e18]
    assert len(report["diagnostics"]["compatible_at_any_scanned_scale_case_ids"]) == 2
    assert all(item["max_normalized_departure"] == 0 for item in report["case_scale_comparisons"])


def test_fixed_linear_full_width_observation_cannot_distinguish_matched_localized_source():
    from mhdlab.resolved_emission_observation import ResolvedTransition, project_resolved_emission
    observed, cases = fixture()
    settings = dict(time_s=[-3e-9, 0.], path_weights_m=[.001, .001],
                    transitions=[ResolvedTransition("synthetic_branch", "X1", 130., 2.)],
                    bin_edges_nm=[129., 130., 131.], gate_start_s=-3e-9, gate_end_s=0.,
                    fwhm_nm=.2, time_registration=observed.time_registration,
                    population_provenance="Manufactured matched joint-source null",
                    geometry_provenance="Equal full-width cell collection weights, no transport or excitation feedback",
                    status="manufactured", optical_model="optically_thin")
    localized = project_resolved_emission(upper_state_density_m3={"X1": [[1e12, 3e12], [2e12, 4e12]]}, **settings)
    uniform = project_resolved_emission(upper_state_density_m3={"X1": [[2e12, 2e12], [3e12, 3e12]]}, **settings)
    y_local = np.asarray(localized["bin_mean_radiance_w_m2_sr_nm"])
    y_uniform = np.asarray(uniform["bin_mean_radiance_w_m2_sr_nm"])
    np.testing.assert_allclose(y_local, y_uniform, rtol=1e-15, atol=0)
    envelope = replace(observed.curves[0], coordinate=[129.5, 130.5], center=y_local,
                       lower=.999*y_local, upper=1.001*y_local, support_segment_ids=[0, 0])
    prediction = replace(cases[3].curves[0], coordinate=[129.5, 130.5], values=y_local, support_segment_ids=[0, 0])
    local_case = replace(cases[3], curves=(prediction,))
    uniform_case = replace(local_case, case_id="uniform_matched_source", hypothesis="uniform_source",
                           curves=(replace(prediction, values=y_uniform),))
    report = compare_neutral_ensemble(replace(observed, curves=(envelope,)), [local_case, uniform_case], [1.])
    assert all(item["max_normalized_departure"] == 0 for item in report["case_scale_comparisons"])


def test_derived_displays_do_not_admit_absolute_scale_scan():
    observed, cases = fixture()
    envelope = replace(observed.curves[0], quantity="transition_density", signal_units="cm^-3", scale_role="published_derived_unscaled")
    observed = replace(observed, curves=(envelope,), mapping_status="scenario_assumption")
    case = replace(cases[3], curves=(replace(cases[3].curves[0], quantity="transition_density", signal_units="cm^-3"),))
    report = compare_neutral_ensemble(observed, [case])
    assert report["calibration"]["scales"] == [1.]
    assert report["calibration"]["scaled_curve_ids"] == []
    assert report["mapping_status"] == "scenario_assumption"
    with pytest.raises(ValueError, match="derived-display-only"):
        compare_neutral_ensemble(observed, [case], exploratory_calibration_scales())
    with pytest.raises(ValueError, match="radiance quantity"):
        replace(envelope, scale_role="nominal_energy_radiance")


def test_unavailable_mapping_withholds_inference():
    observed, cases = fixture()
    with pytest.raises(ValueError, match="withhold inference"):
        compare_neutral_ensemble(replace(observed, mapping_status="unavailable"), cases)


def test_exact_bounds_are_exact_constraints_and_unbounded_mismatches_stay_json_safe():
    observed, cases = fixture()
    envelope = observed.curves[0]
    observed = replace(observed, curves=(replace(envelope, lower=envelope.center, upper=envelope.center),))
    exact, different = cases[3], cases[2]
    report = compare_neutral_ensemble(observed, [exact, different], [1.])
    assert row(report, exact.case_id)["max_normalized_departure"] == 0
    assert row(report, different.case_id)["max_normalized_departure"] is None
    assert row(report, different.case_id)["unbounded_departure"]
    assert not row(report, different.case_id)["readout_compatible"]
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("change", [
    {"reconstruction_id": "another_reconstruction"}, {"gate_start_s": -4e-9},
    {"gate_end_s": 1e-9}, {"time_registration": "another_origin"},
    {"observation_contract": "another_LSF"}, {"excitation_contract": "different_electrons"},
])
def test_incomparable_case_contracts_reject(change):
    observed, cases = fixture()
    with pytest.raises(ValueError, match="mismatched"):
        compare_neutral_ensemble(observed, [replace(cases[3], **change)], [1.])


@pytest.mark.parametrize("change", [
    {"quantity": "unmapped_transition_density"}, {"signal_units": "arbitrary"},
    {"coordinate_units": "pixel"}, {"coordinate": [119., 130., 150.]},
    {"support_segment_ids": [0, 0, 0]},
])
def test_units_quantity_and_retained_support_require_exact_admission(change):
    observed, cases = fixture()
    case = replace(cases[3], curves=(replace(cases[3].curves[0], **change),))
    with pytest.raises(ValueError):
        compare_neutral_ensemble(observed, [case], [1.])


@pytest.mark.parametrize("change", [
    {"conservation_status": "unknown"}, {"population_provenance": ""},
    {"transport_provenance": ""}, {"source_receipt": ""},
    {"optical_admission": ""}, {"numerical_qualification": ""},
    {"thermal_parameter_K": 0}, {"parent_water_areal_inventory_m2": -1},
    {"weighted_neutral_densities_m3": {"H2O[X]": -1}},
    {"neutral_columns_m2": {"H[X]": 1}},
])
def test_missing_physical_qualification_or_invalid_coordinates_reject(change):
    _, cases = fixture()
    with pytest.raises(ValueError):
        replace(cases[3], **change)


def test_no_thermal_descriptor_substitution_or_manufactured_physical_mix():
    observed, cases = fixture()
    with pytest.raises(ValueError, match="thermal coordinate"):
        compare_neutral_ensemble(observed, [cases[3], replace(cases[3], case_id="rotational", thermal_coordinate="rotational_temperature_K")], [1.])
    with pytest.raises(ValueError, match="cannot be mixed"):
        compare_neutral_ensemble(observed, [replace(cases[3], status="owner_supplied")], [1.])


@pytest.mark.parametrize("scales", [[0], [True], [.1, True], [1., .1], [1., 1.], [np.nan], []])
def test_invalid_or_duplicate_scales_reject(scales):
    observed, cases = fixture()
    with pytest.raises(ValueError):
        compare_neutral_ensemble(observed, cases, scales)


def test_arrays_are_snapshots_and_segment_gaps_cannot_reappear():
    observed, cases = fixture()
    with pytest.raises(ValueError):
        observed.curves[0].center[0] = 42
    with pytest.raises(ValueError, match="reappear"):
        replace(observed.curves[0], support_segment_ids=[0, 1, 0])
    with pytest.raises(ValueError, match="readout lower"):
        replace(observed.curves[0], lower=[3., 3.8, 1.7])
    with pytest.raises(ValueError, match="unique"):
        compare_neutral_ensemble(observed, [cases[3], cases[3]])
    with pytest.raises(ValueError, match="finite real"):
        replace(observed.curves[0], center=[True, 4., 2.])
    with pytest.raises(ValueError, match="nonnegative integers"):
        replace(observed.curves[0], support_segment_ids=[False, 0, 1])


@pytest.mark.parametrize("field", [
    "time_registration", "observation_contract", "excitation_contract", "provenance",
    "calibration_provenance",
])
def test_physical_observation_rejects_missing_evidence_even_in_conditional_reconstruction(field):
    observed, _ = fixture()
    with pytest.raises(ValueError, match="missing-evidence sentinel"):
        replace(observed, status="reconstructed", mapping_status="scenario_assumption", **{field: "unavailable"})


@pytest.mark.parametrize("field", [
    "time_registration", "observation_contract", "excitation_contract", "source_receipt",
    "population_provenance", "transport_provenance", "conservation_provenance",
    "optical_admission", "numerical_qualification", "provenance",
    "thermal_coordinate_definition", "neutral_diagnostic_definition",
])
def test_physical_forward_case_rejects_missing_evidence(field):
    _, cases = fixture()
    with pytest.raises(ValueError, match="missing-evidence sentinel"):
        replace(cases[3], status="owner_supplied", **{field: "pending"})


@pytest.mark.parametrize("sentinel", [
    " UNKNOWN ", "unavailable", "<pending>", "missing", "None", "null", "N/A", "TBD.",
    "not_available", "not-provided", "unregistered", "unverified", "not applicable",
])
def test_physical_clock_sentinel_normalization(sentinel):
    observed, _ = fixture()
    with pytest.raises(ValueError, match="time_registration.*missing-evidence sentinel"):
        replace(observed, status="observed", time_registration=sentinel)


def test_physical_nested_curve_provenance_cannot_be_missing():
    observed, cases = fixture()
    with pytest.raises(ValueError, match="curve.*provenance.*missing-evidence sentinel"):
        replace(observed, status="reconstructed", curves=(replace(observed.curves[0], provenance="unknown"),))
    with pytest.raises(ValueError, match="curve.*provenance.*missing-evidence sentinel"):
        replace(cases[3], status="source_evaluated", curves=(replace(cases[3].curves[0], provenance="not supplied"),))


def test_descriptive_manufactured_records_remain_valid_but_do_not_authenticate_physical_evidence():
    observed, cases = fixture()
    # Synthetic fixtures may explain that a physical channel is unavailable;
    # descriptive text is not itself a certified clock/transport measurement.
    case = replace(cases[3], transport_provenance="Manufactured fixture: physical transport output unavailable.")
    report = compare_neutral_ensemble(observed, [case], [1.])
    assert report["case_metadata"][0]["status"] == "manufactured"
    assert report["case_scale_comparisons"][0]["readout_compatible"]
