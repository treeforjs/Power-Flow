"""Independent checks of a prescribed-ion periodic Hall induction reference.

All coefficients are manufactured.  These checks qualify the induction and
energy identities, rather than a coupled fluid or thermal solution.
"""

import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy import integrate

from mhdlab.hall_induction_benchmark import PeriodicHallInduction


ELEMENTARY_CHARGE_C = 1.602176634e-19
VECTOR_FIELDS = (
    "magnetic_field_t", "current_density_a_m2",
    "resistive_electric_field_v_m", "hall_electric_field_v_m",
    "electric_field_v_m", "electron_velocity_m_s",
)
SCALAR_FIELDS = ("joule_power_w_m3", "hall_power_w_m3")


def reference(**overrides):
    parameters = dict(
        resistivity_ohm_m=2e-7, permeability_h_m=1.25e-6,
        electron_density_m3=1e25, background_field_t=3.,
        perturbation_field_t=0.7, wavelength_m=0.004,
        helicity=1, phase_rad=0.23,
    )
    parameters.update(overrides)
    return PeriodicHallInduction(**parameters)


def curl_from_z_derivative(derivative):
    """Independent Cartesian curl for fields depending only on z."""
    result = np.zeros_like(derivative)
    result[:, 0] = -derivative[:, 1]
    result[:, 1] = derivative[:, 0]
    return result


@pytest.mark.parametrize("helicity", [-1, 1])
@pytest.mark.parametrize("background", [-3., 3.])
def test_ampere_faraday_and_signed_hall_rotation(helicity, background):
    model = reference(helicity=helicity, background_field_t=background)
    time = 2.7e-7
    z = model.wavelength_m * np.array([-0.39, -0.07, 0.13, 0.34])
    dz = model.wavelength_m * 2e-6
    center = model.fields(z, time)
    below, above = model.fields(z - dz, time), model.fields(z + dz, time)
    derivative_b = (above["magnetic_field_t"] - below["magnetic_field_t"]) / (2 * dz)
    np.testing.assert_allclose(
        curl_from_z_derivative(derivative_b) / model.permeability_h_m,
        center["current_density_a_m2"], rtol=5e-10, atol=1e-8,
    )
    derivative_e = (above["electric_field_v_m"] - below["electric_field_v_m"]) / (2 * dz)
    dt = 2e-6 / max(model.decay_rate_s, abs(model.angular_frequency_s))
    derivative_time_b = (
        model.fields(z, time + dt)["magnetic_field_t"]
        - model.fields(z, time - dt)["magnetic_field_t"]
    ) / (2 * dt)
    np.testing.assert_allclose(
        derivative_time_b, -curl_from_z_derivative(derivative_e),
        rtol=3e-9, atol=1e-4,
    )
    # For the explicit +sin convention, positive B0 rotates +x toward +y.
    at_origin = reference(helicity=helicity, background_field_t=background,
                          phase_rad=0., resistivity_ohm_m=0.)
    rotated = at_origin.fields([0.], 0.01 / abs(at_origin.angular_frequency_s))
    assert np.sign(rotated["magnetic_field_t"][0, 1]) == np.sign(background)


def test_ohm_law_electron_velocity_and_no_hall_work():
    model = reference(helicity=-1)
    fields = model.fields(np.linspace(-0.005, 0.007, 19), 3e-7)
    magnetic, current = fields["magnetic_field_t"], fields["current_density_a_m2"]
    denominator = model.electron_density_m3 * ELEMENTARY_CHARGE_C
    expected_hall = np.cross(current, magnetic) / denominator
    np.testing.assert_allclose(fields["hall_electric_field_v_m"], expected_hall,
                               rtol=3e-14, atol=1e-12)
    np.testing.assert_allclose(fields["electron_velocity_m_s"], -current / denominator,
                               rtol=3e-14, atol=0.)
    np.testing.assert_allclose(fields["hall_electric_field_v_m"],
                               -np.cross(fields["electron_velocity_m_s"], magnetic),
                               rtol=3e-14, atol=1e-12)
    np.testing.assert_allclose(fields["resistive_electric_field_v_m"],
                               model.resistivity_ohm_m * current, rtol=3e-14, atol=0.)
    np.testing.assert_allclose(fields["electric_field_v_m"],
                               fields["hall_electric_field_v_m"]
                               + fields["resistive_electric_field_v_m"], rtol=3e-14)
    joule = model.resistivity_ohm_m * np.einsum("ij,ij->i", current, current)
    np.testing.assert_allclose(fields["joule_power_w_m3"], joule, rtol=3e-14)
    work_scale = np.linalg.norm(current, axis=1) * np.linalg.norm(expected_hall, axis=1)
    hall_dot = np.einsum("ij,ij->i", current, fields["hall_electric_field_v_m"])
    assert np.all(np.abs(hall_dot) <= 8 * np.finfo(float).eps * work_scale)
    np.testing.assert_array_equal(fields["hall_power_w_m3"], np.zeros(len(current)))
    assert np.all(fields["current_density_a_m2"][:, 2] == 0.)


def test_ideal_hall_rotates_without_damping_or_deposited_heat():
    model = reference(resistivity_ohm_m=0., phase_rad=0.)
    period = 2 * math.pi / abs(model.angular_frequency_s)
    at_quarter = model.fields([0.], period / 4)
    np.testing.assert_allclose(at_quarter["magnetic_field_t"][0],
                               [0., model.perturbation_field_t, model.background_field_t],
                               rtol=0., atol=2e-15)
    z = np.linspace(0., model.wavelength_m, 13)
    initial = model.fields(z, 0.)
    returned = model.fields(z, period)
    np.testing.assert_allclose(returned["magnetic_field_t"], initial["magnetic_field_t"],
                               rtol=0., atol=3e-15)
    np.testing.assert_array_equal(returned["joule_power_w_m3"], np.zeros(len(z)))
    energy = model.energy_budget(100 * period)
    assert energy["joule_heat_j_m3"] == 0.
    assert energy["magnetic_perturbation_energy_j_m3"] == pytest.approx(
        energy["initial_perturbation_energy_j_m3"], rel=2e-15)


def test_resistive_only_decay_preserves_phase_and_background():
    model = reference(background_field_t=0.)
    z, time = np.linspace(-0.006, 0.003, 17), 4e-7
    initial, later = model.fields(z, 0.), model.fields(z, time)
    np.testing.assert_allclose(later["magnetic_field_t"],
                               initial["magnetic_field_t"] * math.exp(-model.decay_rate_s * time),
                               rtol=3e-14, atol=0.)
    np.testing.assert_array_equal(later["hall_electric_field_v_m"], np.zeros((len(z), 3)))
    assert model.angular_frequency_s == 0.


def test_background_reversal_changes_rotation_but_not_heating():
    positive, negative = reference(phase_rad=0.), reference(background_field_t=-3., phase_rad=0.)
    time = 3e-7
    first, second = positive.fields([0.], time), negative.fields([0.], time)
    np.testing.assert_allclose(second["magnetic_field_t"][0],
                               first["magnetic_field_t"][0] * [1., -1., -1.], rtol=2e-14)
    np.testing.assert_array_equal(second["joule_power_w_m3"], first["joule_power_w_m3"])
    assert negative.angular_frequency_s == -positive.angular_frequency_s
    assert negative.energy_budget(time) == positive.energy_budget(time)


def test_cell_averages_match_independent_quadrature_and_preserve_background():
    model = reference(helicity=-1)
    time = 3.1e-7
    edges = np.array([-0.0037, -0.00365, -0.0013, 0.0009, 0.0048])
    result = model.cell_averages(edges, time)
    wave = model.helicity * 2 * math.pi / model.wavelength_m
    decay = model.resistivity_ohm_m * wave**2 / model.permeability_h_m
    omega = model.background_field_t * wave**2 / (
        model.permeability_h_m * model.electron_density_m3 * ELEMENTARY_CHARGE_C)
    amplitude = model.perturbation_field_t * math.exp(-decay * time)

    def physical_fields(position):
        angle = wave * position + omega * time + model.phase_rad
        magnetic = np.array([amplitude * math.cos(angle), amplitude * math.sin(angle),
                             model.background_field_t])
        current = np.array([-wave * magnetic[0], -wave * magnetic[1], 0.]) / model.permeability_h_m
        resistive = model.resistivity_ohm_m * current
        hall = np.cross(current, magnetic) / (model.electron_density_m3 * ELEMENTARY_CHARGE_C)
        return {
            "magnetic_field_t": magnetic,
            "current_density_a_m2": current,
            "resistive_electric_field_v_m": resistive,
            "hall_electric_field_v_m": hall,
            "electric_field_v_m": resistive + hall,
            "electron_velocity_m_s": -current / (model.electron_density_m3 * ELEMENTARY_CHARGE_C),
            "joule_power_w_m3": model.resistivity_ohm_m * float(current @ current),
            "hall_power_w_m3": 0.,
        }

    for key in VECTOR_FIELDS + SCALAR_FIELDS:
        expected = []
        for lower, upper in zip(edges[:-1], edges[1:]):
            integral = integrate.quad_vec(lambda position: physical_fields(position)[key], lower, upper,
                                          epsabs=1e-12, epsrel=5e-13)[0]
            expected.append(integral / (upper - lower))
        np.testing.assert_allclose(result[key], expected, rtol=3e-12, atol=1e-9)
    np.testing.assert_array_equal(result["magnetic_field_t"][:, 2],
                                  np.full(len(edges) - 1, model.background_field_t))


def test_joule_cell_average_is_not_square_of_cell_average_current():
    model = reference()
    result = model.cell_averages([0., model.wavelength_m], 2e-7)
    assert np.linalg.norm(result["current_density_a_m2"][0]) < 1e-6
    assert result["joule_power_w_m3"][0] > 1e10
    point = model.fields([0.], 2e-7)
    assert result["joule_power_w_m3"][0] == pytest.approx(point["joule_power_w_m3"][0], rel=3e-14)
    average_squared = model.resistivity_ohm_m * float(
        result["current_density_a_m2"][0] @ result["current_density_a_m2"][0])
    assert result["joule_power_w_m3"][0] > 1e20 * average_squared


def test_narrow_cell_limits_remain_finite_and_agree_with_points():
    model = reference()
    for lower in (0., 0.0013, -0.0037):
        upper = np.nextafter(lower, np.inf) if lower else 1e-18
        cell = model.cell_averages([lower, upper], 1e-7)
        point = model.fields([lower], 1e-7)
        for key in VECTOR_FIELDS + SCALAR_FIELDS:
            assert np.all(np.isfinite(cell[key]))
            np.testing.assert_allclose(cell[key], point[key], rtol=3e-13, atol=1e-9)


def test_energy_balance_and_time_integral_of_joule_power():
    model, time = reference(), 8e-7
    energy = model.energy_budget(time)
    initial = model.perturbation_field_t**2 / (2 * model.permeability_h_m)
    assert energy["initial_perturbation_energy_j_m3"] == pytest.approx(initial, rel=3e-15)
    final_fields = model.fields([0.00037], time)
    transverse = final_fields["magnetic_field_t"][0, :2]
    assert energy["magnetic_perturbation_energy_j_m3"] == pytest.approx(
        float(transverse @ transverse) / (2 * model.permeability_h_m), rel=3e-14)
    deposited = integrate.quad(lambda tau: model.fields([0.00037], tau)["joule_power_w_m3"][0],
                               0., time, epsabs=1e-9, epsrel=2e-13)[0]
    assert energy["joule_heat_j_m3"] == pytest.approx(deposited, rel=3e-13)
    assert energy["magnetic_perturbation_energy_j_m3"] + energy["joule_heat_j_m3"] == pytest.approx(
        initial, rel=4e-15)
    assert abs(energy["balance_residual_j_m3"]) <= 4 * np.finfo(float).eps * initial
    dt = time * 1e-5
    derivative_heat = (model.energy_budget(time + dt)["joule_heat_j_m3"]
                       - model.energy_budget(time - dt)["joule_heat_j_m3"]) / (2 * dt)
    assert derivative_heat == pytest.approx(final_fields["joule_power_w_m3"][0], rel=2e-10)


def test_energy_initial_tiny_time_depletion_and_zero_perturbation_limits():
    model = reference()
    initial = model.energy_budget(0.)
    assert initial["joule_heat_j_m3"] == 0.
    assert initial["magnetic_perturbation_energy_j_m3"] == initial["initial_perturbation_energy_j_m3"]
    tiny_time = 1e-20 / model.decay_rate_s
    tiny = model.energy_budget(tiny_time)
    # expm1 is needed: the heat is positive even when exp(-2 gamma t)
    # rounds to one in the retained magnetic energy.
    assert tiny["joule_heat_j_m3"] > 0.
    expected = model.fields([0.], 0.)["joule_power_w_m3"][0] * tiny_time
    assert tiny["joule_heat_j_m3"] == pytest.approx(expected, rel=3e-14, abs=0.)
    depleted = model.energy_budget(1000. / model.decay_rate_s)
    assert depleted["magnetic_perturbation_energy_j_m3"] == 0.
    assert depleted["joule_heat_j_m3"] == depleted["initial_perturbation_energy_j_m3"]
    zero = reference(perturbation_field_t=0.)
    zero_energy = zero.energy_budget(1e-6)
    for key in ("initial_perturbation_energy_j_m3", "magnetic_perturbation_energy_j_m3",
                "joule_heat_j_m3", "balance_residual_j_m3"):
        assert zero_energy[key] == 0.
    assert zero_energy["background_magnetic_energy_j_m3"] == pytest.approx(
        zero.background_field_t**2 / (2 * zero.permeability_h_m), rel=3e-15)
    fields = zero.fields([-0.1, 0., 0.1], 1e-6)
    for key in VECTOR_FIELDS:
        if key != "magnetic_field_t":
            np.testing.assert_array_equal(fields[key], np.zeros((3, 3)))
    np.testing.assert_array_equal(fields["magnetic_field_t"],
                                  [[0., 0., zero.background_field_t]] * 3)


def test_scaled_decay_retains_representable_fields_and_energy_after_exponential_underflow():
    # Dimensionless gamma=1 s^-1 and k=1 m^-1 isolate decay evaluation.
    # None of the expected requested outputs exceeds float64 range.
    model = reference(resistivity_ohm_m=1., permeability_h_m=1., electron_density_m3=1e19,
                      background_field_t=0., perturbation_field_t=1e100,
                      wavelength_m=math.tau, phase_rad=0.)
    fields = model.fields([0.], 800.)
    expected_field = math.exp(math.log(model.perturbation_field_t) - 800.)
    assert expected_field > 0.
    assert fields["magnetic_field_t"][0, 0] == pytest.approx(expected_field, rel=3e-13, abs=0.)
    assert fields["current_density_a_m2"][0, 0] == pytest.approx(-expected_field, rel=3e-13, abs=0.)
    energy = model.energy_budget(400.)
    expected_energy = math.exp(math.log(0.5) + 2 * math.log(model.perturbation_field_t) - 800.)
    assert expected_energy > 0.
    assert energy["magnetic_perturbation_energy_j_m3"] == pytest.approx(
        expected_energy, rel=3e-13, abs=0.)


def test_scaled_heat_survives_unrepresentably_small_frozen_hazard():
    model = reference(resistivity_ohm_m=1e-200, permeability_h_m=1.,
                      electron_density_m3=1e30, background_field_t=0.,
                      perturbation_field_t=1e150, wavelength_m=math.tau, phase_rad=0.)
    time = 1e-200
    # gamma*t=1e-400 underflows, but q(0)*dt=1e-100 is representable.
    power = model.fields([0.], 0.)["joule_power_w_m3"][0]
    expected = power * time
    assert expected > 0.
    heat = model.energy_budget(time)["joule_heat_j_m3"]
    assert heat == pytest.approx(expected, rel=3e-13, abs=0.)


def test_maximum_finite_energy_ledger_avoids_intermediate_sum_overflow():
    maximum = np.finfo(float).max
    permeability = maximum / 2
    model = reference(resistivity_ohm_m=permeability / 4,
                      permeability_h_m=permeability, electron_density_m3=maximum,
                      background_field_t=0., perturbation_field_t=maximum,
                      wavelength_m=math.tau, phase_rad=0.)
    energy = model.energy_budget(3.0408086612830174)
    assert energy["initial_perturbation_energy_j_m3"] == maximum
    assert all(math.isfinite(value) for value in energy.values())
    assert 0. < energy["magnetic_perturbation_energy_j_m3"] < maximum
    assert 0. < energy["joule_heat_j_m3"] < maximum
    # The ledger's mathematical residual is finite even when adding its two
    # positive terms before subtracting the initial energy would overflow.
    assert abs(energy["balance_residual_j_m3"]) <= 8 * np.finfo(float).eps * maximum


def test_scaled_rate_is_combined_with_time_before_underflow():
    model = reference(resistivity_ohm_m=1e-200, permeability_h_m=1.,
                      electron_density_m3=1e30, background_field_t=0.,
                      perturbation_field_t=1e150, wavelength_m=math.tau * 1e100,
                      phase_rad=0.)
    time = 1e300
    # The dimensional gamma=1e-400 is below float64 range.  q and its
    # accumulated heat remain representable after the other scales enter.
    power = model.fields([0.], 0.)["joule_power_w_m3"][0]
    expected = power * time
    assert expected == pytest.approx(1e200, rel=3e-14)
    heat = model.energy_budget(time)["joule_heat_j_m3"]
    assert heat == pytest.approx(expected, rel=3e-13, abs=0.)


def test_underflowed_hall_frequency_retains_representable_rotation_at_large_time():
    model = reference(resistivity_ohm_m=0., permeability_h_m=1., electron_density_m3=1e30,
                      background_field_t=1e-200, perturbation_field_t=1e150,
                      wavelength_m=math.tau * 1e100, phase_rad=0.)
    time = 1e300
    # omega is below float64 range; omega*t and A*sin(omega*t) are not.
    # An independently grouped scale avoids the unrepresentable rate.
    phase = (model.background_field_t * time) * 1e-200 / (
        model.permeability_h_m * model.electron_density_m3 * ELEMENTARY_CHARGE_C)
    expected_y = model.perturbation_field_t * math.sin(phase)
    assert expected_y > 0.
    fields = model.fields([0.], time)
    assert fields["magnetic_field_t"][0, 1] == pytest.approx(expected_y, rel=3e-13, abs=0.)


def test_manufactured_export_workflow_preserves_units_hashes_energy_and_existing_output(tmp_path):
    root = Path(__file__).resolve().parents[1]
    config_path = root / "configs/hall_induction_manufactured.json"
    config = json.loads(config_path.read_text())
    output = tmp_path / "hall-reference"
    command = [sys.executable, str(root / "tools/run_hall_induction_benchmark.py"),
               "--output-dir", str(output)]
    completed = subprocess.run(command, cwd=root, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    reported = json.loads(completed.stdout)
    assert reported["snapshots"] == len(config["time_s"])
    summary = json.loads((output / "summary.json").read_text())
    assert summary["schema"] == "mhdlab-hall-induction-run-v1"
    assert summary["parameter_status"] == "manufactured"
    assert summary["parameters"] == config["parameters"]
    assert summary["parameter_provenance"] == config["parameter_provenance"]
    assert summary["config_sha256"] == hashlib.sha256(config_path.read_bytes()).hexdigest()
    assert any("No FLASH or other numerical solver run is claimed" in item for item in summary["limits"])
    assert any("prescribed stationary" in item for item in summary["limits"])
    assert set(summary["source_sha256"]) == {
        "mhdlab/hall_induction_benchmark.py", "mhdlab/thermal_release.py",
        "mhdlab/constants.py", "tools/run_hall_induction_benchmark.py",
    }
    for name, expected in summary["source_sha256"].items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected
    for name, expected in summary["output_sha256"].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == expected
    for time, snapshot in zip(config["time_s"], summary["snapshots"], strict=True):
        path = output / snapshot["csv_file"]
        sidecar = json.loads(path.with_suffix(".json").read_text())
        assert sidecar["schema"] == "mhdlab-hall-induction-reference-v1"
        assert sidecar["time_s"] == time == snapshot["time_s"]
        assert sidecar["parameter_status"] == "manufactured"
        assert sidecar["parameters"] == config["parameters"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert sidecar["csv_sha256"] == digest == snapshot["csv_sha256"]
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            assert set(reader.fieldnames) == set(sidecar["units"])
            rows = list(reader)
        assert len(rows) == config["cells_per_period"]
        assert sidecar["units"]["z_left_m"] == sidecar["units"]["z_right_m"] == "m"
        for name, unit in (("magnetic_field_t", "T"), ("current_density_a_m2", "A m^-2"),
                           ("electric_field_v_m", "V m^-1"), ("electron_velocity_m_s", "m s^-1")):
            assert all(sidecar["units"][f"{name}_{axis}"] == unit for axis in "xyz")
        assert sidecar["units"]["joule_power_w_m3"] == sidecar["units"]["hall_power_w_m3"] == "W m^-3"
        widths = np.array([float(row["z_right_m"]) - float(row["z_left_m"]) for row in rows])
        joule = np.array([float(row["joule_power_w_m3"]) for row in rows])
        assert np.all(widths > 0.)
        assert sum(widths) == pytest.approx(config["parameters"]["wavelength_m"], rel=3e-15)
        energy = sidecar["energy_budget"]
        assert energy == snapshot["energy_budget"]
        expected_power = 2 * summary["decay_rate_s"] * energy["magnetic_perturbation_energy_j_m3"]
        assert float(widths @ joule / sum(widths)) == pytest.approx(expected_power, rel=3e-13)
        assert energy["joule_heat_j_m3"] + energy["magnetic_perturbation_energy_j_m3"] == pytest.approx(
            energy["initial_perturbation_energy_j_m3"], rel=3e-15)
        assert all(float(row["hall_power_w_m3"]) == 0. for row in rows)
    saved = {path.name: path.read_bytes() for path in output.iterdir()}
    rejected = subprocess.run(command, cwd=root, text=True, capture_output=True)
    assert rejected.returncode == 2
    assert "must be a new path" in rejected.stderr
    assert {path.name: path.read_bytes() for path in output.iterdir()} == saved


@pytest.mark.parametrize("name", [
    "resistivity_ohm_m", "permeability_h_m", "electron_density_m3",
    "background_field_t", "perturbation_field_t", "wavelength_m", "phase_rad",
])
@pytest.mark.parametrize("value", [True, "1.0", complex(1., 0.), np.nan, np.inf])
def test_parameters_reject_nonphysical_numeric_coercions(name, value):
    with pytest.raises(ValueError):
        reference(**{name: value})


@pytest.mark.parametrize("name,value", [
    ("resistivity_ohm_m", -1.), ("permeability_h_m", 0.), ("permeability_h_m", -1.),
    ("electron_density_m3", 0.), ("electron_density_m3", -1.),
    ("perturbation_field_t", -1.), ("wavelength_m", 0.), ("wavelength_m", -1.),
    ("helicity", 0), ("helicity", 2), ("helicity", -2), ("helicity", True), ("helicity", "1"),
])
def test_parameters_reject_outside_declared_domain(name, value):
    with pytest.raises(ValueError):
        reference(**{name: value})


@pytest.mark.parametrize("time", [-1., True, "1", np.nan, np.inf, 1 + 0j, [0.]])
def test_time_requires_nonnegative_finite_real_scalar(time):
    model = reference()
    for method, args in ((model.fields, ([0.], time)),
                          (model.cell_averages, ([0., 0.1], time)),
                          (model.energy_budget, (time,))):
        with pytest.raises(ValueError):
            method(*args)


@pytest.mark.parametrize("coordinates", [
    [], [[0., 1.]], [True], ["0"], [1 + 0j], [np.nan], [np.inf],
    np.ma.array([0., 1.], mask=[False, True]),
])
def test_coordinate_samples_require_nonempty_finite_real_vector(coordinates):
    with pytest.raises(ValueError):
        reference().fields(coordinates, 0.)


@pytest.mark.parametrize("edges", [
    [], [0.], [0., 0.], [1., 0.], [0., 1., 0.5], [True, False],
    ["0", "1"], [0j, 1j], [0., np.nan], [0., np.inf], [[0., 1.]],
    np.ma.array([0., 1.], mask=[False, True]),
])
def test_cell_edges_require_strictly_increasing_finite_real_vector(edges):
    with pytest.raises(ValueError):
        reference().cell_averages(edges, 0.)
