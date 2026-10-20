"""Independent CLI checks of log-moment changes and strict report provenance."""

import json
import math
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "tools/verify_temkin_planar_convergence.py"
CONFIG = ROOT / "configs/temkin_convergence_manufactured.json"


def small_event_config():
    config = json.loads(CONFIG.read_text())
    config["reservoirs"][0]["attempt_frequency_s"] = 1e250
    config["event_time_quadrature"]["orders"] = [2, 4]
    config["event_time_quadrature"]["reference_order"] = 8
    return config


def execute(path, output):
    return subprocess.run(
        [sys.executable, str(CLI), "--config", str(path), "--output-dir", str(output),
         "--families", "event_time"], cwd=ROOT, text=True, capture_output=True,
        check=False, timeout=60,
    )


def test_underflowed_si_moments_retain_nonzero_log_relative_quadrature_changes(tmp_path):
    path, output = tmp_path / "underflow.json", tmp_path / "output"
    path.write_text(json.dumps(small_event_config()))
    result = execute(path, output)
    assert result.returncode == 0, result.stderr
    report = json.loads((output / "report.json").read_text())
    assert report["unique_source_cases"] == 0
    study = report["studies"]["event_time"]
    qualification = study["raw_comparison_qualification"].lower()
    assert "underflow" in qualification and "agreement" in qualification and "log" in qualification
    names = {
        "second_raw_moment_s2": ("log_second_raw_moment_s2", "second_raw_moment_underflowed"),
        "variance_offset_s2": ("log_variance_offset_s2", "variance_underflowed"),
    }
    reference = study["estimates"]["8"]
    for order in ("2", "4"):
        estimate = study["estimates"][order]
        for name, (log_name, underflow_flag) in names.items():
            assert estimate[name] == reference[name] == 0.
            assert estimate[underflow_flag] is True and reference[underflow_flag] is True
            log_ratio = estimate[log_name] - reference[log_name]
            assert log_ratio < -.1  # a material change hidden by the SI zeros
            expected_relative = -math.expm1(log_ratio)
            expected_log_change = reference[log_name] + math.log(expected_relative)
            change = study["log_moment_comparisons"][order][name]
            assert change["relative_change"] == pytest.approx(expected_relative, rel=3e-13, abs=0.)
            assert change["relative_change"] > .1
            assert change["log_absolute_change"] == pytest.approx(expected_log_change, rel=0., abs=3e-12)
            assert change["change_underflowed"] is True
            assert change["absolute_change"] == 0.
            assert math.copysign(1., change["signed_change"]) == -1.  # estimate-reference
    json.dumps(report, allow_nan=False)


def test_duplicate_configuration_key_cannot_override_provenance(tmp_path):
    path, output = tmp_path / "duplicate.json", tmp_path / "output"
    body = json.dumps(small_event_config())
    # This is valid JSON syntax with a repeated key, not a malformed-object test.
    path.write_text('{"parameter_status":"sourced",' + body[1:])
    result = execute(path, output)
    assert result.returncode != 0
    assert "duplicate" in result.stderr.lower()
    assert not output.exists()


def test_duplicate_reservoir_names_reject_even_for_an_event_only_study(tmp_path):
    path, output = tmp_path / "duplicate-pool.json", tmp_path / "output"
    config = small_event_config()
    config["reservoirs"].append(dict(config["reservoirs"][0]))
    path.write_text(json.dumps(config))
    result = execute(path, output)
    assert result.returncode != 0
    assert not output.exists()
