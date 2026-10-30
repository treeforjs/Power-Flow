"""Independent elasticity, conservation, and refinement limits."""
import numpy as np
import pytest
import json
from pathlib import Path

from mhdlab.thermal_expansion import (
    SolidExpansionMaterial, load_ss304_material, isotropic_elastic_stress,
    static_bar, transient_bar,
)


def manufactured():
    return SolidExpansionMaterial(np.array([300., 800.]), np.array([0., .01]),
                                  np.array([300., 800.]), np.array([200e9, 200e9]),
                                  8000., 300., "manufactured constant E, linear cumulative expansion")


def test_frozen_grade_data_and_integrated_mean_expansion():
    mat = load_ss304_material()
    np.testing.assert_allclose(mat.strain([373.15, 473.15, 773.15], stress_free_temperature_K=293.15),
                               [80*16e-6, 180*16.5e-6, 480*18e-6], rtol=1e-14)
    # Mean alpha=18e-6 at 500C is NOT a local tangent applied to every interval.
    assert mat.strain(773.15, stress_free_temperature_K=673.15) == pytest.approx(.00199/1.00665)
    assert mat.modulus(523.15) == pytest.approx(182.5e9)
    a = mat.integrated_alpha(773.15, stress_free_temperature_K=373.15)
    b = mat.integrated_alpha(573.15, stress_free_temperature_K=373.15)
    c = mat.integrated_alpha(773.15, stress_free_temperature_K=573.15)
    assert a == pytest.approx(b+c, abs=1e-16)


@pytest.mark.parametrize("temperature", [293.14, 773.16, 1700., np.nan, np.inf])
def test_no_extrapolation_or_nonfinite_material(temperature):
    with pytest.raises(ValueError):
        load_ss304_material().strain(temperature, stress_free_temperature_K=300.)


@pytest.mark.parametrize("boundary", ["free", "fixed_left", "fixed_both"])
def test_zero_change(boundary):
    mat = load_ss304_material()
    x = np.linspace(0, .025, 17)
    result = static_bar(x, np.full(16, 300.), material=mat,
                        stress_free_temperature_K=300., boundary=boundary)
    np.testing.assert_array_equal(result["displacement_m"], 0.)
    np.testing.assert_array_equal(result["stress_Pa"], 0.)
    dynamic = transient_bar(x, [0., 1e-9, 2e-9], np.full((3, 16), 300.), material=mat,
                            stress_free_temperature_K=300., boundary=boundary)
    np.testing.assert_array_equal(dynamic["displacement_m"], 0.)


def test_uniform_free_and_fixed_stress_and_mass():
    mat = manufactured(); x = np.linspace(0, .025, 11)
    free = static_bar(x, np.full(10, 400.), material=mat,
                      stress_free_temperature_K=300., boundary="free", area_m2=1e-6)
    np.testing.assert_allclose(free["displacement_m"], .002*(x-.0125), atol=1e-19)
    np.testing.assert_allclose(free["stress_Pa"], 0., atol=1e-5)
    assert free["reference_mass_kg"] == pytest.approx(free["current_mass_kg"], rel=2e-15)
    fixed = static_bar(x, np.full(10, 400.), material=mat,
                       stress_free_temperature_K=300., boundary="fixed_both", area_m2=1e-6)
    np.testing.assert_allclose(fixed["displacement_m"], 0., atol=1e-19)
    np.testing.assert_allclose(fixed["stress_Pa"], -400e6, rtol=1e-14)
    np.testing.assert_allclose(fixed["reaction_N"], [400., -400.])


def test_three_dimensional_free_vs_fully_constrained_trial_stress():
    mat = manufactured(); eps = .002
    free = isotropic_elastic_stress(np.eye(3)*eps, 400., material=mat,
                                   stress_free_temperature_K=300., poisson_ratio=.3)
    fixed = isotropic_elastic_stress(np.zeros((3, 3)), 400., material=mat,
                                    stress_free_temperature_K=300., poisson_ratio=.3)
    np.testing.assert_allclose(free, 0., atol=1e-5)
    np.testing.assert_allclose(fixed, -1e9*np.eye(3), rtol=1e-14)


def test_nonuniform_equilibrium_variable_modulus_independent_compliance():
    mat = load_ss304_material(); x = np.array([0., .001, .003, .007, .025])
    t = np.array([310., 400., 500., 600.])
    result = static_bar(x, t, material=mat, stress_free_temperature_K=300., boundary="fixed_both")
    expected = -sum(mat.strain(t, stress_free_temperature_K=300.)*np.diff(x))/sum(np.diff(x)/mat.modulus(t))
    np.testing.assert_allclose(result["stress_Pa"], expected, rtol=2e-14)
    assert abs(result["displacement_m"][-1]) < 1e-18


def test_nonuniform_heating_spatial_refinement_against_integral():
    mat = manufactured(); errors = []
    # T=300+100 sin(pi x/L); exact free extension=alpha*100*2L/pi.
    for n in (8, 16, 32):
        x = np.linspace(0, .025, n+1); centers = (x[:-1]+x[1:])/2
        result = static_bar(x, 300+100*np.sin(np.pi*centers/.025), material=mat,
                            stress_free_temperature_K=300., boundary="fixed_left")
        errors.append(abs(result["displacement_m"][-1] - .002*2*.025/np.pi))
    assert min(np.array(errors[:-1])/errors[1:]) > 3.98


def _mode(n, nt, duration_factor=.37):
    mat = manufactured(); length = .0002; speed = np.sqrt(200e9/8000.)
    x = np.linspace(0, length, n+1); times = np.linspace(0, duration_factor*length/speed, nt)
    result = transient_bar(x, times, np.full((nt, n), 300.), material=mat,
                           stress_free_temperature_K=300., boundary="free", area_m2=.005*.025,
                           initial_displacement_m=1e-8*np.cos(np.pi*x/length))
    return x, times, result, speed


def test_newmark_time_refinement_against_exact_discrete_eigenmode():
    errors=[]
    for nt in (11, 21, 41):
        x, times, result, speed = _mode(16, nt)
        omega = 2*speed/(x[1]-x[0])*np.sin(np.pi/(2*16))
        exact = 1e-8*np.cos(omega*times[-1])*np.cos(np.pi*x/x[-1])
        errors.append(np.max(np.abs(result["displacement_m"][-1]-exact)))
    assert min(np.array(errors[:-1])/errors[1:]) > 3.95


def test_dynamic_mesh_refinement_against_continuum_wave():
    errors=[]
    for n in (8, 16, 32):
        x, times, result, speed = _mode(n, 401)
        exact = 1e-8*np.cos(np.pi*speed*times[-1]/x[-1])*np.cos(np.pi*x/x[-1])
        errors.append(np.max(np.abs(result["displacement_m"][-1]-exact)))
    assert min(np.array(errors[:-1])/errors[1:]) > 3.8


def test_unforced_energy_mass_momentum_and_rigid_motion():
    _, _, result, _ = _mode(24, 101)
    energy = result["kinetic_energy_J"]+result["elastic_energy_J"]
    np.testing.assert_allclose(energy, energy[0], rtol=2e-11)
    np.testing.assert_allclose(result["current_mass_kg"], result["reference_mass_kg"], rtol=2e-15)
    momentum = result["velocity_m_s"]@result["nodal_reference_mass_kg"]
    momentum_scale = np.max(np.abs(result["velocity_m_s"])@result["nodal_reference_mass_kg"])
    assert np.max(np.abs(momentum))/momentum_scale < 5e-12
    mat = manufactured(); x=np.linspace(0, .001, 11); times=np.linspace(0, 1e-9, 21)
    rigid=transient_bar(x, times, np.full((21, 10), 300.), material=mat,
                       stress_free_temperature_K=300., boundary="free", initial_velocity_m_s=np.full(11, 3.))
    np.testing.assert_allclose(rigid["displacement_m"], np.broadcast_to(3*times[:, None], (21, 11)), atol=2e-20)


def test_transient_nonuniform_heating_has_inertia_and_constraint_reactions():
    mat=manufactured(); x=np.linspace(0,.025,17); times=np.array([0.,1e-12])
    temperatures=np.stack([np.full(16,300.), np.linspace(301.,400.,16)])
    dynamic=transient_bar(x,times,temperatures,material=mat,stress_free_temperature_K=300.,boundary="free")
    static=static_bar(x,temperatures[-1],material=mat,stress_free_temperature_K=300.,boundary="free")
    assert np.max(abs(dynamic["displacement_m"][-1])) < 1e-5*np.max(abs(static["displacement_m"]))
    fixed=transient_bar([0.,.025],times,[[300.],[400.]],material=mat,stress_free_temperature_K=300.,boundary="fixed_both",area_m2=1e-6)
    np.testing.assert_allclose(fixed["stress_Pa"][-1],[-400e6])
    np.testing.assert_allclose(fixed["reaction_N"][-1],[400.,-400.])


def test_invalid_hot_cell_rejected_before_any_average():
    with pytest.raises(ValueError,match="no clipping"):
        static_bar([0., .001, .002], [300., 3000.], material=load_ss304_material(),
                   stress_free_temperature_K=300., boundary="free")


def test_reject_silent_304l_substitution(tmp_path):
    record=json.loads((Path(__file__).parents[1]/"data/materials/ss304-v1.json").read_text())
    record["material"]["UNS"]="S30403"
    candidate=tmp_path/"substituted.json";candidate.write_text(json.dumps(record))
    with pytest.raises(ValueError,match="Type 304"):
        load_ss304_material(candidate)


def test_guard_large_initial_deformation():
    with pytest.raises(ValueError,match="small-strain"):
        transient_bar([0.,.001],[0.,1e-9],[[300.],[300.]],material=manufactured(),
                      stress_free_temperature_K=300.,boundary="free",initial_displacement_m=[0.,.0001])


def test_material_table_is_immutable():
    material=load_ss304_material()
    with pytest.raises(ValueError):
        material.relative_length_change.setflags(write=True)


def test_underflowing_time_step_cannot_return_nonfinite_mechanics():
    with pytest.raises(ValueError,match="time-step square"):
        transient_bar([0.,.001],[0.,1e-200],[[300.],[300.]],material=manufactured(),
                      stress_free_temperature_K=300.,boundary="free")
