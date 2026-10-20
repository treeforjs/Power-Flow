"""Fresh reference exports retain SI, finite wall heating and provenance."""

import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from mhdlab.heating_comparison import load_heating_profile_csv


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/run_smooth_heating_benchmark.py"


def configuration(tmp_path, **changes):
    config = json.loads((ROOT / "configs/smooth_heating_manufactured.json").read_text())
    config.update(depth_cells=6, times_s=[5e-7, 2e-6])
    config.update(changes)
    path = tmp_path / "parameters.json"
    path.write_text(json.dumps(config) + "\n")
    return path, config


def execute(path, output, *extra):
    return subprocess.run([sys.executable, str(TOOL), "--config", str(path), "--output-dir", str(output), *extra],
                          capture_output=True, text=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_export_round_trip_hashes_energy_and_finite_wall_temperature(tmp_path):
    path, config = configuration(tmp_path)
    output = tmp_path / "reference"
    result = execute(path, output, "--plot")
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["parameter_status"] == "manufactured"
    summary = json.loads((output / "summary.json").read_text())
    assert summary["is_flash_result"] is False
    assert summary["is_measured_shot"] is False
    assert summary["config_sha256"] == sha(path) == sha(output / "input_config.json")
    assert summary["parameter_provenance"] == config["parameter_provenance"]
    assert summary["parameters"] == config["benchmark"]
    assert summary["numerical_qualifications"]["quadrature_relative_tolerance"] > 0
    for name, digest in summary["source_sha256"].items():
        assert sha(ROOT / name) == digest
    for name, digest in summary["output_sha256"].items():
        assert sha(output / name) == digest
    assert (output / "overview.png").is_file()
    for entry in summary["profiles"]:
        profile = load_heating_profile_csv(output / entry["csv"])
        assert profile.source_kind == "synthetic"
        assert "linear-ramp" in profile.source_reference
        assert profile.time_s == entry["time_s"]
        assert len(profile.magnetic_field_t) == 6
        assert sha(output / entry["csv"]) == entry["csv_sha256"]
        assert np.all(profile.joule_power_w_m3 >= 0)
    with (output / "energy_budget.csv").open() as stream:
        budgets = list(csv.DictReader(stream))
    assert len(budgets) == 2
    for row in budgets:
        incoming, magnetic, heat = (float(row[name]) for name in
                                    ("poynting_input_j_m2", "magnetic_energy_j_m2", "joule_heat_j_m2"))
        assert abs(incoming - magnetic - heat) <= 3e-14 * incoming
    with (output / "adiabatic_point_temperature.csv").open() as stream:
        temperatures = list(csv.DictReader(stream))
    assert len(temperatures) == 2 * (6 + 1)
    wall = [row for row in temperatures if row["point_kind"] == "wall"]
    assert len(wall) == 2
    parameters = config["benchmark"]
    for row, reported in zip(wall, summary["wall_adiabatic_temperature_rise_k"]):
        assert float(row["point_depth_m"]) == 0
        time = float(row["time_s"])
        expected = (2 * parameters["boundary_field_rate_t_s"] ** 2 * time ** 2
                    / (parameters["permeability_h_m"] * math.pi * parameters["volumetric_heat_capacity_j_m3_k"]))
        assert math.isfinite(reported)
        np.testing.assert_allclose(reported, expected, rtol=2e-14)
        assert float(row["adiabatic_temperature_rise_k"]) == reported
    assert not list(tmp_path.glob(".smooth-heating-*"))


@pytest.mark.parametrize("rate", [0., -2e6])
def test_signed_and_zero_drive_remain_explicit_in_profiles_and_energy(tmp_path, rate):
    path, config = configuration(tmp_path, depth_cells=3, times_s=[1e-6])
    config["benchmark"]["boundary_field_rate_t_s"] = rate
    path.write_text(json.dumps(config))
    output = tmp_path / "reference"
    result = execute(path, output)
    assert result.returncode == 0, result.stderr
    summary = json.loads((output / "summary.json").read_text())
    profile = load_heating_profile_csv(output / summary["profiles"][0]["csv"])
    if rate == 0:
        np.testing.assert_array_equal(profile.magnetic_field_t, 0)
        np.testing.assert_array_equal(profile.current_density_a_m2, 0)
        np.testing.assert_array_equal(profile.joule_power_w_m3, 0)
        assert summary["finite_depth_boundary_B_fraction"] == [None]
        assert summary["wall_adiabatic_temperature_rise_k"] == [0.]
    else:
        assert np.all(profile.magnetic_field_t < 0)
        assert np.all(profile.current_density_a_m2 < 0)
        assert np.all(profile.joule_power_w_m3 > 0)
        assert summary["wall_adiabatic_temperature_rise_k"][0] > 0


@pytest.mark.parametrize("parameters,time,depth", [
    ({"conductivity_s_m": 1e15, "permeability_h_m": 1e308, "volumetric_heat_capacity_j_m3_k": 1e308,
      "boundary_field_rate_t_s": 1e308}, 10., 1e-160),
    ({"conductivity_s_m": 1e308, "permeability_h_m": 1e-308, "volumetric_heat_capacity_j_m3_k": 1e-308,
      "boundary_field_rate_t_s": 1e-308}, 1e-20, 1e-10),
])
def test_dimensionless_boundary_fraction_survives_unrepresentable_drive_scale(tmp_path, parameters, time, depth):
    from scipy.integrate import quad
    path, _ = configuration(tmp_path, benchmark=parameters, times_s=[time], depth_stop_m=depth, depth_cells=1)
    output = tmp_path / "reference"
    result = execute(path, output)
    assert result.returncode == 0, result.stderr
    summary = json.loads((output / "summary.json").read_text())
    fraction = summary["finite_depth_boundary_B_fraction"][0]
    coordinate = math.exp(math.log(depth) - math.log(2.)
                          + .5 * (math.log(parameters["permeability_h_m"]) + math.log(parameters["conductivity_s_m"]) - math.log(time)))
    expected = 4 / math.sqrt(math.pi) * quad(lambda w: w*w * math.exp(-(coordinate+w)**2),
                                           0., np.inf, epsabs=0., epsrel=2e-12)[0]
    assert fraction is not None and fraction > 0
    np.testing.assert_allclose(fraction, expected, rtol=1e-10)


@pytest.mark.parametrize("existing", ["directory", "file", "symlink"])
def test_existing_destination_is_preserved_without_partial_generation(tmp_path, existing):
    path, _ = configuration(tmp_path)
    output = tmp_path / "reference"
    sentinel = b"retain existing work"
    if existing == "directory":
        output.mkdir()
        (output / "sentinel").write_bytes(sentinel)
    elif existing == "file":
        output.write_bytes(sentinel)
    else:
        output.symlink_to(tmp_path / "missing-target", target_is_directory=True)
    result = execute(path, output)
    assert result.returncode != 0
    assert "must be a new path" in result.stderr
    if existing == "directory":
        assert list(output.iterdir()) == [output / "sentinel"]
        assert (output / "sentinel").read_bytes() == sentinel
    elif existing == "file":
        assert output.read_bytes() == sentinel
    else:
        assert output.is_symlink()
        assert not (tmp_path / "missing-target").exists()
    assert not list(tmp_path.glob(".smooth-heating-*"))


@pytest.mark.parametrize("kind", ["zero_time", "missing_constant", "boolean_length", "status", "duplicate_json", "nonfinite_json"])
def test_invalid_or_ambiguous_inputs_leave_no_result(tmp_path, kind):
    path, config = configuration(tmp_path)
    if kind == "zero_time":
        config["times_s"] = [0.]
    elif kind == "missing_constant":
        del config["benchmark"]["volumetric_heat_capacity_j_m3_k"]
    elif kind == "boolean_length":
        config["depth_stop_m"] = True
    elif kind == "status":
        config["parameter_status"] = "measured-shot"
    path.write_text(json.dumps(config))
    if kind == "duplicate_json":
        path.write_text('{"schema":"x","schema":"y"}')
    elif kind == "nonfinite_json":
        path.write_text('{"value":NaN}')
    output = tmp_path / "reference"
    result = execute(path, output)
    assert result.returncode != 0
    assert not output.exists()
    assert not list(tmp_path.glob(".smooth-heating-*"))
