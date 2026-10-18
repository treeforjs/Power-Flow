"""Independent physical checks of the constant-property half-space reference.

Parameters here are manufactured numerical fixtures, not measured metal data.
"""

from dataclasses import FrozenInstanceError
import inspect
import math

import numpy as np
import pytest
from scipy import integrate, special

from mhdlab.heating_benchmarks import ConstantConductivityHalfSpace


def reference(**overrides):
    parameters = dict(conductivity_s_m=3e6, permeability_h_m=4e-7 * math.pi,
                      volumetric_heat_capacity_j_m3_k=2.8e6, boundary_field_t=5.0)
    parameters.update(overrides)
    return ConstantConductivityHalfSpace(**parameters)


def test_diffusion_solution_boundary_sign_and_local_ohmic_identity():
    model = reference()
    time = 2e-6
    diffusivity = 1 / (model.permeability_h_m * model.conductivity_s_m)
    assert model.magnetic_diffusivity_m2_s == pytest.approx(diffusivity, rel=2e-15)
    scale = math.sqrt(diffusivity * time)
    depths = np.array([0, 0.2, 1, 3, 10]) * scale
    result = model.fields(depths, time)
    np.testing.assert_allclose(result.magnetic_field_t,
                               model.boundary_field_t * special.erfc(depths / (2 * scale)), rtol=3e-14)
    np.testing.assert_allclose(result.joule_power_w_m3,
                               result.current_density_a_m2**2 / model.conductivity_s_m, rtol=2e-14)
    assert result.magnetic_field_t[0] == pytest.approx(model.boundary_field_t)
    assert np.all(result.current_density_a_m2 > 0)
    negative = reference(boundary_field_t=-5).fields(depths, time)
    np.testing.assert_array_equal(negative.magnetic_field_t, -result.magnetic_field_t)
    np.testing.assert_array_equal(negative.current_density_a_m2, -result.current_density_a_m2)
    np.testing.assert_array_equal(negative.joule_power_w_m3, result.joule_power_w_m3)


def test_current_obeys_differential_ampere_law_and_magnetic_diffusion():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    y = np.array([0.2, 1, 2, 4]) * scale
    step = 2e-4 * scale
    minus, middle, plus = [model.fields(y + offset, time) for offset in (-step, 0, step)]
    gradient = (plus.magnetic_field_t - minus.magnetic_field_t) / (2 * step)
    np.testing.assert_allclose(-gradient / model.permeability_h_m,
                               middle.current_density_a_m2, rtol=3e-8)
    time_step = time * 2e-4
    time_derivative = (model.fields(y, time + time_step).magnetic_field_t
                       - model.fields(y, time - time_step).magnetic_field_t) / (2 * time_step)
    curvature = (plus.magnetic_field_t - 2 * middle.magnetic_field_t
                 + minus.magnetic_field_t) / step**2
    np.testing.assert_allclose(time_derivative, model.magnetic_diffusivity_m2_s * curvature, rtol=2e-7)


def test_analytic_global_energy_and_surface_poynting_closure():
    model, time = reference(), 2e-6
    energy = model.energy_budget(time)
    prefactor = (model.boundary_field_t**2 / model.permeability_h_m
                 * math.sqrt(model.magnetic_diffusivity_m2_s * time / math.pi))
    assert energy.magnetic_energy_j_m2 == pytest.approx(prefactor * (2 - math.sqrt(2)), rel=5e-15)
    assert energy.joule_heat_j_m2 == pytest.approx(prefactor * math.sqrt(2), rel=5e-15)
    assert energy.poynting_input_j_m2 == pytest.approx(2 * prefactor, rel=5e-15)
    assert abs(energy.energy_residual_j_m2) < 5e-15 * energy.poynting_input_j_m2
    surface = model.fields(0., time)
    # Ohm's E_z and tangential B_x give inward S_y = E_z*B_x/mu.
    inward_flux = (float(surface.current_density_a_m2) / model.conductivity_s_m
                   * model.boundary_field_t / model.permeability_h_m)
    assert energy.surface_power_w_m2 == pytest.approx(inward_flux, rel=1e-14)
    dt = time * 1e-5
    deposited_power = (model.energy_budget(time + dt).poynting_input_j_m2
                       - model.energy_budget(time - dt).poynting_input_j_m2) / (2 * dt)
    assert deposited_power == pytest.approx(inward_flux, rel=1e-9)


def test_independent_space_time_quadrature_and_truncated_domain_convergence():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    energy = model.energy_budget(time)
    # Integrate dimensionless independent formulas.  The accumulated heat
    # integral uses E1, not the analytic global Q formula implemented above.
    magnetic_integral = integrate.quad(lambda u: math.erfc(u / 2)**2,
                                       0, np.inf, epsabs=2e-12, epsrel=2e-12)[0]
    magnetic = model.boundary_field_t**2 / (2 * model.permeability_h_m) * scale * magnetic_integral
    heat_integral = integrate.quad(lambda u: special.exp1(u*u / 2),
                                   0, np.inf, epsabs=2e-11, epsrel=2e-11)[0]
    heat = model.boundary_field_t**2 / (model.permeability_h_m * math.pi) * scale * heat_integral
    assert magnetic == pytest.approx(energy.magnetic_energy_j_m2, rel=3e-12)
    assert heat == pytest.approx(energy.joule_heat_j_m2, rel=3e-12)
    truncated_u, truncated_q = [], []
    for upper in (0.5, 1, 2, 4, 8):
        integral_b = integrate.quad(lambda u: math.erfc(u / 2)**2, 0, upper)[0]
        integral_q = integrate.quad(lambda u: special.exp1(u*u / 2), 0, upper)[0]
        truncated_u.append(model.boundary_field_t**2 / (2 * model.permeability_h_m) * scale * integral_b)
        truncated_q.append(model.boundary_field_t**2 / (model.permeability_h_m * math.pi) * scale * integral_q)
    assert np.all(np.diff(truncated_u) > 0)
    assert np.all(np.diff(truncated_q) > 0)
    assert truncated_u[-1] == pytest.approx(energy.magnetic_energy_j_m2, rel=1e-13)
    assert truncated_q[-1] == pytest.approx(energy.joule_heat_j_m2, rel=1e-13)
    # Independent time integral of instantaneous integrated Joule power.
    power = lambda tau: model.boundary_field_t**2 / model.permeability_h_m * math.sqrt(
        model.magnetic_diffusivity_m2_s / (2 * math.pi * tau))
    integrated_power = integrate.quad(power, 0, time, epsabs=1e-6, epsrel=2e-12)[0]
    assert integrated_power == pytest.approx(energy.joule_heat_j_m2, rel=3e-13)


def test_cell_averages_ampere_and_independent_quadrature_not_square_mean_current():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    edges = np.array([0, 0.01, 0.3, 1, 3, 8]) * scale
    cells = model.cell_averages(edges, time)
    boundary_b = model.fields(edges, time).magnetic_field_t
    np.testing.assert_allclose(cells.current_density_a_m2 * np.diff(edges),
                               (boundary_b[:-1] - boundary_b[1:]) / model.permeability_h_m, rtol=1e-12)
    expected = np.zeros((3, len(edges) - 1))
    for index, (lower, upper) in enumerate(zip(edges[:-1] / scale, edges[1:] / scale)):
        for row, integrand in enumerate((lambda u: model.boundary_field_t * math.erfc(u / 2),
                                        lambda u: model.boundary_field_t / (model.permeability_h_m * math.sqrt(math.pi) * scale)
                                        * math.exp(-u*u/4),
                                        lambda u: model.boundary_field_t**2 / (model.permeability_h_m * math.pi * time)
                                        * math.exp(-u*u/2))):
            expected[row, index] = integrate.quad(integrand, lower, upper, epsabs=0., epsrel=2e-13)[0] / (upper-lower)
    np.testing.assert_allclose(cells.magnetic_field_t, expected[0], rtol=2e-12)
    np.testing.assert_allclose(cells.current_density_a_m2, expected[1], rtol=2e-12)
    np.testing.assert_allclose(cells.joule_power_w_m3, expected[2], rtol=2e-12)
    assert np.all(cells.joule_power_w_m3 >= cells.current_density_a_m2**2 / model.conductivity_s_m)
    assert cells.joule_power_w_m3[-1] > 2 * cells.current_density_a_m2[-1]**2 / model.conductivity_s_m


def test_cell_means_are_stable_for_adjacent_float_edges_and_far_tails():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    for y in (0., scale, 8 * scale, 50 * scale):
        upper = np.nextafter(y, np.inf) if y else scale * 1e-14
        cells = model.cell_averages([y, upper], time)
        point = model.fields(y, time)
        for average, expected in zip((cells.magnetic_field_t, cells.current_density_a_m2, cells.joule_power_w_m3),
                                     (point.magnetic_field_t, point.current_density_a_m2, point.joule_power_w_m3)):
            assert average[0] == pytest.approx(float(expected), rel=2e-12, abs=0.)
            assert np.isfinite(average[0])


def test_similarity_scaling_and_temperature_time_derivative():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    depths = np.array([0.05, 0.5, 2, 5]) * scale
    fields = model.fields(depths, time)
    scaled = model.fields(3 * depths, 9 * time)
    np.testing.assert_allclose(scaled.magnetic_field_t, fields.magnetic_field_t, rtol=1e-14)
    np.testing.assert_allclose(3 * scaled.current_density_a_m2, fields.current_density_a_m2, rtol=1e-14)
    np.testing.assert_allclose(9 * scaled.joule_power_w_m3, fields.joule_power_w_m3, rtol=2e-14)
    np.testing.assert_allclose(model.adiabatic_temperature_rise(3 * depths, 9 * time),
                               model.adiabatic_temperature_rise(depths, time), rtol=1e-14)
    assert model.energy_budget(9 * time).joule_heat_j_m2 == pytest.approx(3 * model.energy_budget(time).joule_heat_j_m2)
    delta = time * 1e-5
    derivative = (model.adiabatic_temperature_rise(depths, time + delta)
                  - model.adiabatic_temperature_rise(depths, time - delta)) / (2 * delta)
    np.testing.assert_allclose(derivative, fields.joule_power_w_m3 / model.volumetric_heat_capacity_j_m3_k,
                               rtol=2e-9)


def test_temperature_pointwise_time_integration_and_integrable_surface_singularity():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    for depth in (0.1 * scale, scale, 3 * scale):
        integrated = integrate.quad(lambda fraction: float(model.fields(depth, time * fraction).joule_power_w_m3)
                                    * time / model.volumetric_heat_capacity_j_m3_k,
                                    0, 1, epsabs=1e-11, epsrel=1e-11)[0]
        assert float(model.adiabatic_temperature_rise(depth, time)) == pytest.approx(integrated, rel=2e-11)
    approaching = model.adiabatic_temperature_rise(scale * np.array([1e-3, 1e-6, 1e-12, 1e-200]), time)
    assert np.all(np.isfinite(approaching)) and np.all(np.diff(approaching) > 0)
    for depths in (0., [0, scale]):
        with pytest.raises(ValueError, match="positive"):
            model.adiabatic_temperature_rise(depths, time)


def test_temperature_keeps_representable_scaled_e1_tail():
    # The unscaled E1(1000) underflows, while this manufactured dimensional
    # product is finite.  Independent integral: exp(z)*E1(z)=int exp(-u)/(z+u)du.
    model = reference(conductivity_s_m=1., permeability_h_m=1.,
                      volumetric_heat_capacity_j_m3_k=1., boundary_field_t=1e155)
    argument = 1000.
    scaled_integral = integrate.quad(lambda u: math.exp(-u)/(argument + u),
                                     0, np.inf, epsabs=1e-15, epsrel=1e-12)[0]
    expected = math.exp(2 * math.log(model.boundary_field_t) - math.log(math.pi) - argument) * scaled_integral
    actual = float(model.adiabatic_temperature_rise(math.sqrt(2 * argument), 1.))
    assert actual > 0 and math.isfinite(actual)
    assert actual == pytest.approx(expected, rel=2e-12, abs=0.)


def test_zero_drive_outputs_and_immutable_shape_preserving_arrays():
    model = reference(boundary_field_t=0)
    fields = model.fields([[0, 1], [2, 3]], 1)
    cells = model.cell_averages([0, 1, 2], 1)
    assert fields.magnetic_field_t.shape == (2, 2)
    for array in (fields.depth_m, fields.magnetic_field_t, fields.current_density_a_m2,
                  fields.joule_power_w_m3, cells.depth_edges_m, cells.joule_power_w_m3,
                  model.adiabatic_temperature_rise([1, 2], 1)):
        if array is not fields.depth_m and array is not cells.depth_edges_m:
            assert np.all(array == 0)
        with pytest.raises(ValueError):
            array.flags.writeable = True
    assert model.energy_budget(1).poynting_input_j_m2 == 0
    with pytest.raises(FrozenInstanceError):
        model.boundary_field_t = 1
    with pytest.raises(ValueError, match="positive"):
        model.adiabatic_temperature_rise(0, 1)


@pytest.mark.parametrize("name", ["conductivity_s_m", "permeability_h_m", "volumetric_heat_capacity_j_m3_k"])
@pytest.mark.parametrize("value", [0, -1, True, 1+0j, "1", np.nan, np.inf])
def test_material_parameter_guards(name, value):
    with pytest.raises(ValueError):
        reference(**{name: value})


@pytest.mark.parametrize("value", [True, 1j, "1", np.inf, np.nan])
def test_boundary_field_guards(value):
    with pytest.raises(ValueError):
        reference(boundary_field_t=value)


@pytest.mark.parametrize("time", [0, -1, True, 1j, "1", np.inf, np.nan])
def test_time_guards_for_every_api(time):
    model = reference()
    for function, args in ((model.fields, ([0, 1], time)),
                           (model.cell_averages, ([0, 1], time)),
                           (model.energy_budget, (time,)),
                           (model.adiabatic_temperature_rise, ([1, 2], time))):
        with pytest.raises(ValueError):
            function(*args)


@pytest.mark.parametrize("depths", [[], [-1], [True], [1, False], [1j], ["1"], [np.inf], [np.nan],
                                   np.ma.array([1.], mask=[True])])
def test_depth_guards(depths):
    with pytest.raises(ValueError):
        reference().fields(depths, 1)


@pytest.mark.parametrize("edges", [[0], [[0, 1]], [0, 0], [1, 0], [0, 1, 1]])
def test_cell_geometry_guards(edges):
    with pytest.raises(ValueError):
        reference().cell_averages(edges, 1)


def test_no_physical_defaults_and_explicit_float_range_failure():
    signature = inspect.signature(ConstantConductivityHalfSpace)
    assert all(parameter.default is inspect.Parameter.empty for parameter in signature.parameters.values())
    with pytest.raises(TypeError):
        ConstantConductivityHalfSpace()
    with pytest.raises(ValueError, match="float64"):
        reference(conductivity_s_m=1e-300, permeability_h_m=1e-300)
    with pytest.raises(ValueError, match="float64"):
        reference(conductivity_s_m=1e300, permeability_h_m=1e300)
    with pytest.raises(ValueError, match="float64"):
        reference(boundary_field_t=1e308).energy_budget(1)
