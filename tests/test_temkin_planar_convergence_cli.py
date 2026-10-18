"""Exercise the complete refinement report without a FLASH installation."""

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from tools.verify_temkin_planar_convergence import array_record, load_config


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/verify_temkin_planar_convergence.py"


def config():
    value = json.loads((ROOT / "configs/temkin_convergence_manufactured.json").read_text())
    value["base_resolution"] = dict(width_bins=4, current_bins=2, snapshot_intervals=4, subdivisions_per_interval=2)
    value["refinements"] = {
        "width_bins": {"levels": [2, 4, 8], "reference": 16},
        "current_bins": {"levels": [1, 2, 4], "reference": 8},
        "snapshot_intervals": {"levels": [2, 4, 8], "reference": 16},
        "subdivisions_per_interval": {"levels": [1, 2, 4], "reference": 8},
    }
    value["event_time_quadrature"].update(orders=[2, 4], reference_order=8)
    return value


def invoke(tmp_path, value, *extra):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(value))
    out = tmp_path / "result"
    run = subprocess.run([sys.executable, str(TOOL), "--config", str(path), "--output-dir", str(out), *extra],
                         capture_output=True, text=True)
    return run, out, path


@pytest.fixture(scope="module")
def report(tmp_path_factory):
    path = tmp_path_factory.mktemp("convergence")
    run, out, source = invoke(path, config())
    assert run.returncode == 0, run.stderr
    return json.loads((out / "report.json").read_text()), out, source


def test_separate_axes_and_reuse_identical_base_case(report):
    data, _, _ = report
    axes = {"width": "width_bins", "current": "current_bins", "snapshot": "snapshot_intervals", "integration": "subdivisions_per_interval"}
    base = data["config"]["base_resolution"]
    for name, key in axes.items():
        family = data["studies"][name]
        for case_id in [*family["case_ids"], family["finite_reference_case_id"]]:
            resolution = data["cases"][case_id]["resolution"]
            assert all(resolution[other] == base[other] for other in base if other != key)
        assert family["finite_reference_resolution"] > max(family["levels"])
    assert data["unique_source_cases"] == 13  # Four sets of four, sharing one base case.
    assert len(data["cases"]) == 13
    assert data["parameter_status"] == "manufactured"
    assert data["reference_kind"] == "finite_finest_resolution_reference"


def test_orders_use_three_solutions_instead_of_finite_reference_errors(report):
    import math
    data, _, _ = report
    checked = 0
    for family in data["studies"].values():
        for record in family.get("self_orders", []):
            for name, result in record["metrics"].items():
                order = result["observed_order"]
                if order is None:
                    assert result["warnings"]
                    continue
                values = [data["cases"][case_id]["observables"]["metrics"][name][0] for case_id in record["case_ids"]]
                expected = math.log(abs((values[0] - values[1]) / (values[1] - values[2]))) / math.log(2)
                assert order == pytest.approx(expected, abs=2e-13)
                checked += 1
    assert checked > 0


def test_report_has_reproducible_input_code_array_and_output_hashes(report):
    data, out, source = report
    assert data["input_config_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    for name, expected in data["source_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected
    for name, expected in data["output_sha256"].items():
        assert hashlib.sha256((out / name).read_bytes()).hexdigest() == expected
    for case in data["cases"].values():
        assert all(abs(value) < 1e-14 for value in case["global_inventory_balance_relative"])
        assert len(case["source_arrays"]["released_m2"]["shape"]) == 4
        for group in ("input_arrays", "source_arrays"):
            assert all(len(value["sha256"]) == 64 for value in case[group].values())
    # Shape is part of identity even when flattened bytes match.
    assert array_record([[1., 2.]])["sha256"] != array_record([1., 2.])["sha256"]


def test_event_family_retains_conditional_law_and_finite_quadrature_qualification(report):
    data, _, _ = report
    family = data["studies"]["event_time"]
    assert family["finite_reference_order"] == 8
    assert family["estimates"]["8"]["mean_offset_s"] > 0
    assert family["estimates"]["8"]["law_identity"] == data["law_identity"]
    assert family["estimates"]["8"]["time_identity"].startswith("temkin_conditional")
    assert "certified" in family["qualification"]


def test_existing_result_is_never_replaced(report):
    data, out, source = report
    before = (out / "report.json").read_bytes()
    result = subprocess.run([sys.executable, str(TOOL), "--config", str(source), "--output-dir", str(out)], capture_output=True, text=True)
    assert result.returncode == 2
    assert "new path" in result.stderr
    assert (out / "report.json").read_bytes() == before


@pytest.mark.parametrize("mutation", [
    lambda d: d.update(parameter_status="sourced"),
    lambda d: d["base_resolution"].update(width_bins=True),
    lambda d: d["refinements"]["width_bins"].update(levels=[2, 4, 9]),
    lambda d: d["refinements"]["width_bins"].update(reference=8),
    lambda d: d["reservoirs"][0].update(initial_coverage_m2=[[1e18]]),
    lambda d: d["state_masses_kg"].update({"M:manufactured-ground": 0}),
    lambda d: d["event_time_quadrature"].update(orders=[4, 2]),
    lambda d: d["event_time_quadrature"]["component"].update(reservoir_index=True),
    lambda d: d.update(include_time_moments="false"),
])
def test_invalid_refinement_inputs_rejected(tmp_path, mutation):
    value = copy.deepcopy(config())
    mutation(value)
    path = tmp_path / "invalid.json"
    path.write_text(json.dumps(value))
    with pytest.raises((ValueError, TypeError)):
        load_config(path)


def test_invalid_metric_does_not_leave_partial_artifacts(tmp_path):
    value = config()
    value["selected_metrics"] = ["not_a_source_observable"]
    run, out, _ = invoke(tmp_path, value, "--families", "width")
    assert run.returncode == 2
    assert "unavailable" in run.stderr
    assert not out.exists()


def test_event_only_study_needs_no_planar_cases(tmp_path):
    run, out, _ = invoke(tmp_path, config(), "--families", "event_time")
    assert run.returncode == 0, run.stderr
    data = json.loads((out / "report.json").read_text())
    assert data["unique_source_cases"] == 0
    assert list(data["studies"]) == ["event_time"]
