"""Periodic exports preserve finite-cell physics, explicit inputs and hashes."""

import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from mhdlab.conductive_heating_benchmark import PeriodicConductiveHeating


ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools/run_conductive_heating_benchmark.py"


def configuration(tmp_path, **changes):
    config = json.loads((ROOT / "configs/conductive_heating_manufactured.json").read_text())
    config.update(cells_per_period=8, times_s=[0.0, 2e-7, 1e-6])
    config.update(changes)
    path = tmp_path / "parameters.json"
    path.write_text(json.dumps(config) + "\n")
    return path, config


def execute(path, output, *extra):
    return subprocess.run([sys.executable, str(TOOL), "--config", str(path), "--output-dir", str(output), *extra],
                          capture_output=True, text=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_export_finite_cell_physics_energy_units_and_hashes(tmp_path):
    path, config = configuration(tmp_path, periods=2, position_origin_m=-0.001)
    output = tmp_path / "reference"
    process = execute(path, output, "--plot")
    assert process.returncode == 0, process.stderr
    assert json.loads(process.stdout)["parameter_status"] == "manufactured"
    summary = json.loads((output / "summary.json").read_text())
    assert summary["schema"] == "mhdlab-periodic-conductive-heating-run-v1"
    assert summary["is_flash_result"] is False
    assert summary["is_measured_shot"] is False
    assert summary["config_sha256"] == sha(path) == sha(output / "input_config.json")
    assert summary["parameters"] == config["benchmark"]
    assert summary["parameter_provenance"] == config["parameter_provenance"]
    assert summary["cell_units"]["temperature_rise_k"] == "K"
    assert summary["cell_units"]["joule_power_w_m3"] == "W m^-3"
    assert summary["energy_budget_units"]["mean_thermal_energy_j_m3"] == "J m^-3"
    assert summary["numerical_qualifications"]["maximum_absolute_phase_rad"] > 0
    for name, digest in summary["source_sha256"].items():
        assert sha(ROOT / name) == digest
    for name, digest in summary["output_sha256"].items():
        assert sha(output / name) == digest
    assert (output / "overview.png").is_file()
    model = PeriodicConductiveHeating(**config["benchmark"])
    with (output / "cell_averages.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    count = config["periods"] * config["cells_per_period"]
    assert len(rows) == count * len(config["times_s"])
    edges = np.array([float(row["left_position_m"]) for row in rows[:count]] + [float(rows[count - 1]["right_position_m"])])
    for interval, time in enumerate(config["times_s"]):
        block = rows[interval * count:(interval + 1) * count]
        expected = model.cell_averages(edges, time)
        assert [int(row["cell_index"]) for row in block] == list(range(count))
        assert {float(row["time_s"]) for row in block} == {time}
        for name in ("magnetic_field_t", "current_density_a_m2", "joule_power_w_m3", "temperature_rise_k"):
            np.testing.assert_array_equal([float(row[name]) for row in block], getattr(expected, name))
        if time == 0:
            assert all(float(row["temperature_rise_k"]) == 0 for row in block)
    with (output / "energy_budget.csv").open() as stream:
        energy_rows = list(csv.DictReader(stream))
    for row, recorded in zip(energy_rows, summary["energy_budgets"]):
        initial = float(row["initial_mean_magnetic_energy_j_m3"])
        magnetic, thermal = float(row["mean_magnetic_energy_j_m3"]), float(row["mean_thermal_energy_j_m3"])
        assert abs(initial - magnetic - thermal) < 5e-14 * initial
        assert float(row["mean_total_energy_j_m3"]) == recorded["mean_total_energy_j_m3"]
    assert not list(tmp_path.glob(".conductive-heating-*"))


@pytest.mark.parametrize("field,kappa", [(0.0, 0.0), (-2.0, 0.0), (-2.0, 200000.0)])
def test_explicit_zero_and_signed_field_and_zero_conduction(tmp_path, field, kappa):
    path, config = configuration(tmp_path, times_s=[1e-6])
    config["benchmark"].update(initial_field_t=field, thermal_conductivity_w_m_k=kappa)
    path.write_text(json.dumps(config))
    output = tmp_path / "reference"
    process = execute(path, output)
    assert process.returncode == 0, process.stderr
    summary = json.loads((output / "summary.json").read_text())
    assert summary["parameters"]["thermal_conductivity_w_m_k"] == kappa
    assert summary["derived_coefficients"]["thermal_diffusivity_m2_s"] == kappa / config["benchmark"]["volumetric_heat_capacity_j_m3_k"]
    with (output / "cell_averages.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if field == 0:
        for row in rows:
            assert all(float(row[key]) == 0 for key in ("magnetic_field_t", "current_density_a_m2", "joule_power_w_m3", "temperature_rise_k"))
    else:
        assert float(rows[0]["magnetic_field_t"]) < 0
        assert float(rows[0]["current_density_a_m2"]) < 0
        assert all(float(row["joule_power_w_m3"]) > 0 for row in rows)
        assert all(float(row["temperature_rise_k"]) > 0 for row in rows)


@pytest.mark.parametrize("existing", ["directory", "file", "symlink"])
def test_existing_outputs_are_preserved_without_staging(tmp_path, existing):
    path, _ = configuration(tmp_path)
    output = tmp_path / "reference"
    if existing == "directory":
        output.mkdir()
        (output / "sentinel").write_text("preserve")
    elif existing == "file":
        output.write_text("preserve")
    else:
        output.symlink_to(tmp_path / "missing-target", target_is_directory=True)
    process = execute(path, output)
    assert process.returncode != 0
    assert "new path" in process.stderr
    if existing == "directory":
        assert (output / "sentinel").read_text() == "preserve"
        assert len(list(output.iterdir())) == 1
    elif existing == "file":
        assert output.read_text() == "preserve"
    else:
        assert output.is_symlink()
        assert not (tmp_path / "missing-target").exists()
    assert not list(tmp_path.glob(".conductive-heating-*"))


@pytest.mark.parametrize("kind", [
    "negative_time", "duplicate_time", "missing_cv", "missing_kappa", "unknown_constant",
    "boolean_origin", "boolean_count", "status", "nonstring_status", "empty_provenance",
    "excess_memory", "large_phase", "duplicate_json", "nonfinite_json",
])
def test_bad_configuration_leaves_no_partial_output(tmp_path, kind):
    path, config = configuration(tmp_path)
    if kind == "negative_time":
        config["times_s"] = [-1]
    elif kind == "duplicate_time":
        config["times_s"] = [1, 1]
    elif kind == "missing_cv":
        del config["benchmark"]["volumetric_heat_capacity_j_m3_k"]
    elif kind == "missing_kappa":
        del config["benchmark"]["thermal_conductivity_w_m_k"]
    elif kind == "unknown_constant":
        config["benchmark"]["mass_specific_cp_j_kg_k"] = 500
    elif kind == "boolean_origin":
        config["position_origin_m"] = True
    elif kind == "boolean_count":
        config["cells_per_period"] = True
    elif kind == "status":
        config["parameter_status"] = "measured-shot"
    elif kind == "nonstring_status":
        config["parameter_status"] = ["manufactured"]
    elif kind == "empty_provenance":
        config["parameter_provenance"] = " "
    elif kind == "excess_memory":
        config.update(periods=100000, cells_per_period=100000)
    elif kind == "large_phase":
        config["position_origin_m"] = 1e8
    path.write_text(json.dumps(config))
    if kind == "duplicate_json":
        path.write_text('{"schema":"first","schema":"second"}')
    elif kind == "nonfinite_json":
        path.write_text('{"value":NaN}')
    output = tmp_path / "reference"
    process = execute(path, output)
    assert process.returncode != 0
    assert not output.exists()
    assert not list(tmp_path.glob(".conductive-heating-*"))
