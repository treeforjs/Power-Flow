"""Check the manufactured convergence report and finite-reference distinction."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from tools.verify_surface_convergence import _finite_reference_errors, _self_orders


ROOT = Path(__file__).resolve().parents[1]


def test_tiny_cli_reports_separate_reservoir_energies_hashes_and_finite_references(tmp_path):
    result = subprocess.run([
        sys.executable, str(ROOT / "tools/verify_surface_convergence.py"),
        "--time-subdivisions", "1", "2", "--time-reference-subdivisions", "4",
        "--time-width-bins", "4", "--input-time-intervals", "4",
        "--width-bins", "4", "8", "--width-reference-bins", "16",
        "--width-time-subdivisions", "2", "--output-dir", str(tmp_path), "--plot",
    ], cwd=tmp_path, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    report = json.loads((tmp_path / "report.json").read_text())
    assert report["source_kind"] == "synthetic"
    assert report["parameter_status"] == "manufactured"
    config_path = Path(report["configuration"]["path"])
    assert report["configuration"]["sha256"] == hashlib.sha256(config_path.read_bytes()).hexdigest()
    assert report["configuration"]["values"] == json.loads(config_path.read_text())
    for name, digest in report["source_files_sha256"].items():
        assert digest == hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
    assert (tmp_path / "convergence.png").stat().st_size > 1000
    for study_name in ("temporal_subdivision", "width_refinement"):
        study = report[study_name]
        assert study["reference_kind"] == "independent_finite_resolution_not_exact"
        assert len(study["cases"]) == 2
        assert len(study["self_convergence_orders"]) == 1
        reference = study["reference"]
        for case in [*study["cases"], reference]:
            assert len(case["totals_per_reservoir"]) == 2
            assert case["surface_attributes"]["source_kind"] == "synthetic"
            assert len(case["surface_arrays_sha256"]) == 64
            for pool in case["totals_per_reservoir"]:
                assert pool["released_particles"] > 0
                assert pool["total_wall_energy_debit_j"] == pytest.approx(
                    pool["translational_energy_j"] + pool["binding_energy_j"] + pool["internal_energy_j"], rel=3e-14)
                assert pool["internal_energy_j"] == 0
        for case in study["cases"]:
            for pool, errors in enumerate(case["errors_per_reservoir"]):
                error = errors["quantity_errors"]["released_particles"]
                actual = case["totals_per_reservoir"][pool]["released_particles"]
                target = reference["totals_per_reservoir"][pool]["released_particles"]
                assert error["relative_error"] == pytest.approx(abs(actual - target) / target)
                assert errors["quantity_errors"]["internal_energy_j"]["relative_error"] is None
        for pool in study["self_convergence_orders"][0]["orders_per_reservoir"]:
            assert pool["quantity_orders"]["internal_energy_j"]["order"] is None
    assert any("snapshots" in limitation for limitation in report["limitations"])
    assert any("manufactured" in limitation for limitation in report["limitations"])
    assert json.loads(result.stdout)["warning_count"] == len(report["warnings"])


def test_second_order_self_differences_avoid_finite_reference_order_bias():
    # Q(n)=Q_exact+C/n^2 has exact self-convergence order two, whereas errors
    # against the finite n=16 result give ratio five between n=4 and n=8.
    levels = [4, 8, 16]
    cases = []
    for level in levels:
        value = 10 + 64 / level ** 2
        cases.append({"totals_per_reservoir": [{
            "reservoir": "manufactured", "state_id": "M",
            "released_particles": value, "translational_energy_j": value,
            "internal_energy_j": 0, "binding_energy_j": value, "total_wall_energy_debit_j": 2 * value,
        }]})
    coarse = _finite_reference_errors(cases[0], cases[-1])[0]["quantity_errors"]["released_particles"]["absolute_error"]
    medium = _finite_reference_errors(cases[1], cases[-1])[0]["quantity_errors"]["released_particles"]["absolute_error"]
    assert coarse / medium == 5
    orders = _self_orders(cases, levels)[0]["orders_per_reservoir"][0]["quantity_orders"]
    assert orders["released_particles"]["order"] == pytest.approx(2)
    assert orders["total_wall_energy_debit_j"]["order"] == pytest.approx(2)
    assert orders["internal_energy_j"]["order"] is None
