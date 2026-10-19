"""Width assembly controls using explicitly synthetic ledger records only.

The fixtures satisfy the typed ledger algebra but are not trajectories of the
Arrhenius IVP. No scalar evolution is executed or represented as production truth.
"""
from dataclasses import fields, replace
from fractions import Fraction
import math

import pytest

import mhdlab.thermal_forced_feedback_reference as scalar_data
import mhdlab.thermal_width_power_reference as width_assembly
from mhdlab.thermal_forced_feedback_reference import (
    FeedbackProblem, FeedbackReservoir, ForcedFeedbackHistory, ForcedFeedbackSample,
    PowerSchedule, PowerSegment,
)
from mhdlab.thermal_width_power_reference import (
    WidthPowerGeometry, assemble_width_power_reference, scalar_history_sha256,
)


KB = Fraction("1.380649e-23")
POWER = {"uniform_power": 6e7, "edge_power": 9e7, "middle_power": 3e7}
RELEASE = {"uniform_power": 2**70, "edge_power": 2**71, "middle_power": 2**69}
HYPOTHESES = {"uniform": ["uniform_power"]*3,
              "nonuniform": ["edge_power", "middle_power", "edge_power"]}
LABEL = "Manufactured assembly-only synthetic ledger; not an ODE solution."


@pytest.fixture(autouse=True)
def forbid_scalar_evolution(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("adapter tests must not execute scalar evolution")
    monkeypatch.setattr(scalar_data, "evolve_forced_feedback", forbidden)
    calls = []
    def synthetic_rhs(problem, *, coverage_m2, temperature_k, power_w_m2):
        calls.append((coverage_m2, temperature_k, power_w_m2))
        translation = 2*KB*Fraction.from_float(temperature_k)*Fraction.from_float(coverage_m2)
        # The prescribed fixture rate is1/s, not an Arrhenius result.
        return {
            "release_flux_m2_s": float(coverage_m2), "coverage_rate_m2_s": -float(coverage_m2) if coverage_m2 else 0.0,
            "temperature_rate_k_s": float((Fraction.from_float(power_w_m2)-translation)/Fraction.from_float(problem.areal_heat_capacity_j_m2_k)),
            "input_power_j_m2_s": float(power_w_m2), "binding_power_j_m2_s": 0.0,
            "internal_power_j_m2_s": 0.0, "translational_power_j_m2_s": float(translation),
            "total_wall_power_j_m2_s": float(translation),
        }
    monkeypatch.setattr(width_assembly, "forced_feedback_rhs", synthetic_rhs)
    return calls


def geometry(**changes):
    values = dict(width_edges_m=(-1.0, -0.5, 0.5, 1.0), common_length_m=2.0,
                  parameter_status="manufactured", parameter_provenance=LABEL)
    values.update(changes)
    return WidthPowerGeometry(**values)


def synthetic_history(case, *, steps=2, state_id="synthetic_parent", area=1.0, empty=False):
    """Construct self-consistent records; the release sequence is supplied data."""
    reservoir = FeedbackReservoir(
        name="synthetic_pool", state_id=state_id, attempt_frequency_s=1e8,
        activation_energy_ev=0.35, desorption_energy_ev=0.0, internal_energy_ev=0.0,
        coverage_coefficient=0.6, reference_coverage_m2=float(2**75), max_coverage_m2=float(2**75),
        kinetic_provenance=LABEL, energy_provenance=LABEL,
    )
    initial = 0.0 if empty else float(2**74)
    problem = FeedbackProblem(
        reservoir=reservoir, initial_coverage_m2=initial, initial_temperature_k=900.0,
        areal_heat_capacity_j_m2_k=1.0, patch_area_m2=area,
        temperature_min_k=600.0, temperature_max_k=1300.0,
        parameter_status="manufactured", parameter_provenance=LABEL, energy_reference=LABEL,
    )
    duration = 2**-20
    segment = PowerSegment(0.0, duration, POWER[case], steps, LABEL+case)
    schedule = PowerSchedule(duration, (segment,), "manufactured", LABEL)
    samples = []
    for index in range(steps+1):
        time = Fraction.from_float(duration)*index/steps
        released = 0.0 if empty else float(Fraction(RELEASE[case])*index*2/steps)
        remaining = initial-released
        energy_input = float(Fraction.from_float(POWER[case])*time)
        amount = Fraction.from_float(released)
        # A declared linear-T emission interpolation closes this fixture's
        # thermal/translation ledger; it does not supply a kinetic time law.
        temperature = float((Fraction(900)+Fraction.from_float(energy_input)-KB*amount*900)/(1+KB*amount))
        translation = float(KB*amount*(900+Fraction.from_float(temperature)))
        thermal = Fraction.from_float(temperature)-900
        residual = float(thermal+Fraction.from_float(translation)-Fraction.from_float(energy_input))
        quantities = {field.name: 0.0 for field in fields(ForcedFeedbackSample)
                      if field.name not in {"segment_index", "segment_step", "time_s", "extensive"}}
        quantities.update(
            coverage_m2=remaining, temperature_k=temperature, released_m2=released,
            input_energy_j_m2=energy_input, prescribed_input_energy_j_m2=energy_input,
            translational_energy_j_m2=translation, total_wall_energy_debit_j_m2=translation,
            gas_energy_j_m2=translation, wall_thermal_change_j_m2=float(thermal),
            wall_debit_residual_j_m2=residual, total_energy_residual_j_m2=residual,
        )
        extensive = {}
        for name, value in quantities.items():
            if name.endswith("_j_m2"):
                extensive[name[:-5]+"_j"] = float(Fraction.from_float(value)*Fraction.from_float(area))
            elif name in {"coverage_m2", "released_m2", "inventory_residual_m2"}:
                extensive[name[:-3]+"_number"] = float(Fraction.from_float(value)*Fraction.from_float(area))
        samples.append(ForcedFeedbackSample(segment_index=0, segment_step=index, time_s=float(time),
                                            **quantities, extensive=extensive))
    return ForcedFeedbackHistory(problem, schedule, tuple(samples))


def synthetic_histories(**changes):
    return {case: synthetic_history(case, **changes) for case in POWER}


def assemble(histories=None, *, shape=None, hypotheses=None, hashes=None, provenance=None):
    source = synthetic_histories() if histories is None else histories
    return assemble_width_power_reference(
        source, geometry=geometry() if shape is None else shape,
        hypotheses={name: list(value) for name, value in HYPOTHESES.items()} if hypotheses is None else hypotheses,
        history_sha256={name: scalar_history_sha256(value) for name, value in source.items()} if hashes is None else hashes,
        history_provenance={name: LABEL for name in source} if provenance is None else provenance,
        parameter_status="manufactured", parameter_provenance=LABEL,
    )


def rational(record):
    return Fraction(int(record["numerator"]), int(record["denominator"]))


def assert_quantity(record, expected):
    assert rational(record["exact"]) == expected
    assert record["value"] == float(expected)
    error = Fraction.from_float(record["value"])-expected
    flags = record["projection"]
    assert rational(flags["exact_rounding_error"]) == error
    assert flags["rounded"] is bool(error)
    assert flags["exact_zero"] is (expected == 0)
    assert flags["underflowed"] is flags["overflowed"] is False


def test_exact_geometry_retains_rounding_of_nonbinary_area_products():
    shape = geometry(width_edges_m=(-0.1, -0.05, 0.05, 0.1), common_length_m=0.3)
    source = synthetic_histories()
    report = assemble(source, shape=shape).to_dict()
    expected = (Fraction.from_float(0.1)-Fraction.from_float(0.05))*Fraction.from_float(0.3)
    assert rational(report["geometry"]["patch_areas_exact_m2"][0]) == expected
    assert expected != Fraction.from_float(float(expected))
    for index, weight in enumerate((Fraction(1, 4), Fraction(1, 2), Fraction(1, 4))):
        assert rational(report["geometry"]["area_weights"][index]) == weight
    record = report["hypotheses"]["uniform"]["patches"][0]["samples"][1]["extensive"]["input_energy_j"]
    assert_quantity(record, expected*Fraction.from_float(source["uniform_power"].samples[1].input_energy_j_m2))
    assert record["projection"]["rounded"] is True


def test_all_extensives_scale_once_before_exact_patch_and_total_sums():
    source = synthetic_histories()
    report = assemble(source).to_dict()
    for hypothesis, assignments in HYPOTHESES.items():
        totals = report["hypotheses"][hypothesis]["totals"]
        for sample_index in range(3):
            expected_totals = {}
            for patch_index, (case, area) in enumerate(zip(assignments, (1, 2, 1))):
                original = source[case].samples[sample_index]
                expected = {key: Fraction.from_float(value)*area for key, value in original.extensive.items()}
                expected["release_rate_s"] = Fraction.from_float(original.coverage_m2)*area
                actual = report["hypotheses"][hypothesis]["patches"][patch_index]["samples"][sample_index]["extensive"]
                assert set(actual) == set(expected)
                for key, value in expected.items():
                    assert_quantity(actual[key], value)
                    expected_totals[key] = expected_totals.get(key, Fraction())+value
            for key, value in expected_totals.items():
                assert_quantity(totals[sample_index]["extensive"][key], value)


def test_matched_power_inventory_and_difference_are_exact_without_repair():
    report = assemble().to_dict()
    assert_quantity(report["comparison"]["exact_total_input_power_w"], Fraction(240_000_000))
    uniform = report["hypotheses"]["uniform"]["totals"]
    nonuniform = report["hypotheses"]["nonuniform"]["totals"]
    for name in ("uniform", "nonuniform"):
        assert_quantity(report["hypotheses"][name]["totals"][0]["extensive"]["coverage_number"], Fraction(4*2**74))
    for index, comparison in enumerate(report["comparison"]["samples"]):
        time = rational(comparison["logical_time_s"])
        assert_quantity(comparison["exact_prescribed_input_energy_j"], 240_000_000*time)
        assert_quantity(comparison["prescribed_input_difference_j"], Fraction())
        assert_quantity(comparison["integrated_input_difference_j"], Fraction())
        for key, record in comparison["difference_nonuniform_minus_uniform"].items():
            expected = rational(nonuniform[index]["extensive"][key]["exact"])-rational(uniform[index]["extensive"][key]["exact"])
            assert_quantity(record, expected)
    assert rational(nonuniform[-1]["extensive"]["released_number"]["exact"]) != rational(uniform[-1]["extensive"]["released_number"]["exact"])


def test_three_source_histories_are_reused_with_the_complete_patch_state(forbid_scalar_evolution):
    source = synthetic_histories()
    report = assemble(source).to_dict()
    assert set(report["scalar_histories"]) == set(POWER)
    assert len(forbid_scalar_evolution) == 9
    assert report["cache"]["rhs_evaluations"] == report["cache"]["scalar_samples"] == 9
    assert report["cache"]["distinct_scalar_histories"] == 3
    assert report["cache"]["patch_history_references"] == 6
    assert report["cache"]["scalar_evolutions_by_assembly"] == 0
    rhs_fields = {"release_flux_m2_s", "coverage_rate_m2_s", "temperature_rate_k_s",
                  "input_power_j_m2_s", "binding_power_j_m2_s", "internal_power_j_m2_s",
                  "translational_power_j_m2_s", "total_wall_power_j_m2_s"}
    for case, record in report["scalar_histories"].items():
        assert len(record["rhs"]) == 3
        for index, values in enumerate(record["rhs"]):
            assert set(values) == rhs_fields
            assert values["release_flux_m2_s"] == record["release_flux_m2_s"][index]
            assert values["coverage_rate_m2_s"] == -source[case].samples[index].coverage_m2
    for hypothesis, assignments in HYPOTHESES.items():
        for patch, case in zip(report["hypotheses"][hypothesis]["patches"], assignments):
            assert patch["scalar_history_id"] == case
            assert patch["history_sha256"] == scalar_history_sha256(source[case])
            for row, original in zip(patch["samples"], source[case].samples):
                assert row["state_id"] == source[case].problem.reservoir.state_id
                assert row["areal"]["temperature_k"] == original.temperature_k
                assert row["areal"]["released_m2"] == original.released_m2
                assert row["logical_time_s"] == report["scalar_histories"][case]["logical_time_s"][original.segment_step]


def test_empty_parent_inventory_keeps_true_zero_release_and_energy():
    report = assemble(synthetic_histories(empty=True)).to_dict()
    for hypothesis in ("uniform", "nonuniform"):
        for row in report["hypotheses"][hypothesis]["totals"]:
            for key in ("coverage_number", "released_number", "translational_energy_j", "total_wall_energy_debit_j", "release_rate_s"):
                assert_quantity(row["extensive"][key], Fraction())
    assert report["host_budget_accepted"] is False
    assert report["owner_implementation_verified"] is False


def test_assembly_snapshot_and_detached_documents_do_not_alias_inputs():
    source = synthetic_histories()
    result = assemble(source)
    before = result.to_dict()
    detached = result.to_dict()
    detached["hypotheses"]["uniform"]["patches"][0]["samples"][1]["areal"]["temperature_k"] = -1
    object.__setattr__(source["uniform_power"].samples[1], "temperature_k", 0.0)
    assert result.to_dict() == before
    with pytest.raises(TypeError):
        result.document["parameter_status"] = "sourced"


@pytest.mark.parametrize("changes", [
    {"width_edges_m": (-1.0, -0.5, -0.5, 1.0)},
    {"width_edges_m": (1.0, 0.5, -0.5, -1.0)},
    {"width_edges_m": (-1.0, False, 0.5, 1.0)},
    {"width_edges_m": (-1.0, math.nan, 0.5, 1.0)},
    {"common_length_m": True}, {"common_length_m": math.inf},
    {"common_length_m": Fraction(-1, 10**1000)},
    {"parameter_status": "sourced"},
])
def test_invalid_geometry_and_truth_labels_are_rejected(changes):
    with pytest.raises(ValueError):
        geometry(**changes)


def test_missing_history_is_rejected():
    source = synthetic_histories()
    source.pop("middle_power")
    with pytest.raises(ValueError):
        assemble(source)


def test_extra_history_is_rejected():
    source = synthetic_histories()
    source["unrequested_case"] = source["middle_power"]
    with pytest.raises(ValueError):
        assemble(source)


def test_source_hash_is_bound_to_the_supplied_complete_history():
    source = synthetic_histories()
    hashes = {name: scalar_history_sha256(value) for name, value in source.items()}
    hashes["middle_power"] = "0"*64
    with pytest.raises(ValueError):
        assemble(source, hashes=hashes)


@pytest.mark.parametrize("provenance", [
    {"uniform_power": LABEL, "edge_power": LABEL},
    {"uniform_power": LABEL, "edge_power": LABEL, "middle_power": ""},
])
def test_complete_nonempty_source_provenance_is_required(provenance):
    with pytest.raises(ValueError):
        assemble(provenance=provenance)


@pytest.mark.parametrize("hypotheses", [
    {"uniform": ["uniform_power"]*3},
    {"uniform": ["uniform_power"]*3, "nonuniform": ["edge_power"]*3},
    {"uniform": ["uniform_power"]*3, "nonuniform": ["edge_power", "unknown", "edge_power"]},
])
def test_hypotheses_cannot_be_silently_filled_or_renormalized(hypotheses):
    with pytest.raises(ValueError):
        assemble(hypotheses=hypotheses)


@pytest.mark.parametrize("kind", ["unit_area", "state_identity", "logical_grid"])
def test_incompatible_history_context_is_rejected(kind):
    source = synthetic_histories()
    if kind == "unit_area":
        source["middle_power"] = synthetic_history("middle_power", area=2.0)
    elif kind == "state_identity":
        source["middle_power"] = synthetic_history("middle_power", state_id="other_synthetic_parent")
    else:
        source["middle_power"] = synthetic_history("middle_power", steps=1)
    with pytest.raises(ValueError):
        assemble(source)


@pytest.mark.parametrize("kind", ["time_support", "boolean_primary", "extensive", "derived_residual"])
def test_altered_sample_records_are_rejected_even_with_updated_hashes(kind):
    source = synthetic_histories()
    history = source["middle_power"]
    samples = list(history.samples)
    if kind == "time_support":
        samples[-1] = replace(samples[-1], time_s=math.nextafter(samples[-1].time_s, math.inf))
    elif kind == "boolean_primary":
        samples[-1] = replace(samples[-1], temperature_k=True)
    elif kind == "extensive":
        amounts = dict(samples[-1].extensive)
        amounts["released_number"] *= 2
        samples[-1] = replace(samples[-1], extensive=amounts)
    else:
        samples[-1] = replace(samples[-1], total_energy_residual_j_m2=1.0)
    source["middle_power"] = replace(history, samples=tuple(samples))
    with pytest.raises(ValueError):
        assemble(source)


def test_tiny_exact_power_mismatch_is_rejected_without_tolerance_or_scaling():
    # The visually equal0.2 widths differ as exact represented inputs.
    with pytest.raises(ValueError):
        shape = geometry(width_edges_m=(0.0, 0.1, 0.3, 0.4), common_length_m=1.0)
        assemble(shape=shape)


@pytest.mark.parametrize("changes", [
    {"width_edges_m": (-2*math.ulp(0.0), -math.ulp(0.0), math.ulp(0.0), 2*math.ulp(0.0)), "common_length_m": math.ulp(0.0)},
    {"width_edges_m": (-1e308, -5e307, 5e307, 1e308), "common_length_m": 1e308},
])
def test_nonzero_unrepresentable_geometry_is_rejected_without_floors(changes):
    with pytest.raises(ValueError):
        geometry(**changes)


def test_a_history_cannot_be_relabeled_as_another_prescribed_power():
    source = synthetic_histories()
    history = source["middle_power"]
    segment = replace(history.schedule.segments[0], power_w_m2=3e7+1.0)
    source["middle_power"] = replace(history, schedule=replace(history.schedule, segments=(segment,)))
    with pytest.raises(ValueError):
        assemble(source)


def test_initial_inventory_snapshot_is_revalidated():
    source = synthetic_histories()
    history = source["middle_power"]
    samples = (replace(history.samples[0], coverage_m2=0.0), *history.samples[1:])
    source["middle_power"] = replace(history, samples=samples)
    with pytest.raises(ValueError):
        assemble(source)
