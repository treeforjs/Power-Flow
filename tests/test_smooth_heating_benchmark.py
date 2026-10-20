"""Independent linear-ramp diffusion/energy checks using manufactured inputs."""

from dataclasses import FrozenInstanceError
import inspect
import math

import numpy as np
import pytest
from scipy import integrate

from mhdlab.smooth_heating_benchmark import LinearRampConductivityHalfSpace


def reference(**overrides):
    values = dict(conductivity_s_m=3e6, permeability_h_m=4e-7 * math.pi,
                  volumetric_heat_capacity_j_m3_k=2.8e6, boundary_field_rate_t_s=2e6)
    values.update(overrides)
    return LinearRampConductivityHalfSpace(**values)


def direct_b(u):
    return (1 + 2*u*u) * math.erfc(u) - 2*u * math.exp(-u*u) / math.sqrt(math.pi)


def direct_g(u):
    return math.exp(-u*u) / math.sqrt(math.pi) - u * math.erfc(u)


def test_linear_ramp_fields_match_independent_duhamel_convolution():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    coordinates = np.array([0., .1, .7, 1.9, 2., 4.])
    result = model.fields(2 * scale * coordinates, time)
    for index, u in enumerate(coordinates):
        b = integrate.quad(lambda v: math.erfc(u / math.sqrt(v)), 0, 1,
                           epsabs=1e-14, epsrel=1e-12)[0]
        j = integrate.quad(lambda v: math.exp(-u*u/v) / math.sqrt(v), 0, 1,
                           epsabs=1e-14, epsrel=1e-12)[0]
        expected_b = model.boundary_field_rate_t_s * time * b
        expected_j = (model.boundary_field_rate_t_s / model.permeability_h_m
                      * math.sqrt(time / (math.pi * model.magnetic_diffusivity_m2_s)) * j)
        assert result.magnetic_field_t[index] == pytest.approx(expected_b, rel=3e-11, abs=0.)
        assert result.current_density_a_m2[index] == pytest.approx(expected_j, rel=3e-11, abs=0.)
    np.testing.assert_allclose(result.joule_power_w_m3,
                               result.current_density_a_m2**2 / model.conductivity_s_m, rtol=3e-14)
    assert result.magnetic_field_t[0] == pytest.approx(model.boundary_field_rate_t_s * time)
    assert np.all(result.current_density_a_m2 > 0)


def test_differential_ampere_and_diffusion_pde():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    depths = np.array([.15, .6, 1., 2.5]) * scale
    step = 2e-4 * scale
    below, middle, above = [model.fields(depths + offset, time) for offset in (-step, 0, step)]
    gradient = (above.magnetic_field_t - below.magnetic_field_t) / (2 * step)
    np.testing.assert_allclose(-gradient / model.permeability_h_m,
                               middle.current_density_a_m2, rtol=2e-8)
    curvature = (above.magnetic_field_t - 2 * middle.magnetic_field_t + below.magnetic_field_t) / step**2
    dt = time * 1e-4
    derivative = (model.fields(depths, time + dt).magnetic_field_t
                  - model.fields(depths, time - dt).magnetic_field_t) / (2 * dt)
    np.testing.assert_allclose(derivative, model.magnetic_diffusivity_m2_s * curvature, rtol=3e-7)


def test_cells_match_independent_integrals_and_integral_ampere_not_square_mean_j():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    u_edges = np.array([0., .01, .4, 1., 2., 4., 20.])
    edges = 2 * scale * u_edges
    cells = model.cell_averages(edges, time)
    points = model.fields(edges, time)
    np.testing.assert_allclose(cells.current_density_a_m2 * np.diff(edges),
                               (points.magnetic_field_t[:-1] - points.magnetic_field_t[1:]) / model.permeability_h_m,
                               rtol=2e-11, atol=0.)
    expected = np.zeros((3, len(edges) - 1))
    for index, (lower, upper) in enumerate(zip(u_edges[:-1], u_edges[1:])):
        b = integrate.quad(direct_b, lower, upper, epsabs=1e-16, epsrel=2e-12)[0] / (upper-lower)
        g = integrate.quad(direct_g, lower, upper, epsabs=1e-16, epsrel=2e-12)[0] / (upper-lower)
        g2 = integrate.quad(lambda u: direct_g(u)**2, lower, upper,
                            epsabs=1e-24, epsrel=2e-12)[0] / (upper-lower)
        expected[:, index] = [model.boundary_field_rate_t_s * time * b,
                              2 * model.boundary_field_rate_t_s * time / (model.permeability_h_m * scale) * g,
                              4 * model.boundary_field_rate_t_s**2 * time / model.permeability_h_m * g2]
    for array, correct in zip((cells.magnetic_field_t, cells.current_density_a_m2, cells.joule_power_w_m3), expected):
        np.testing.assert_allclose(array, correct, rtol=5e-10, atol=0.)
    assert np.all(cells.joule_power_w_m3 >= cells.current_density_a_m2**2 / model.conductivity_s_m)
    assert cells.joule_power_w_m3[-1] > 10 * cells.current_density_a_m2[-1]**2 / model.conductivity_s_m


def test_narrow_cells_match_point_limits():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    for depth in (0., .2*scale, 4*scale, 20*scale):
        upper = np.nextafter(depth, np.inf) if depth else scale * 1e-14
        cells = model.cell_averages([depth, upper], time)
        points = model.fields(depth, time)
        for name in ("magnetic_field_t", "current_density_a_m2", "joule_power_w_m3"):
            assert getattr(cells, name)[0] == pytest.approx(float(getattr(points, name)), rel=3e-11, abs=0.)


def test_broad_finite_volume_keeps_integrated_skin_layer_in_extreme_bin():
    model = reference(conductivity_s_m=1., permeability_h_m=1.,
                      volumetric_heat_capacity_j_m3_k=1., boundary_field_rate_t_s=1.)
    cells = model.cell_averages([0., 1e300], 1.)
    assert cells.magnetic_field_t[0] == pytest.approx(4/(3*math.sqrt(math.pi))*1e-300, rel=3e-12, abs=0.)
    assert cells.current_density_a_m2[0] == pytest.approx(1e-300, rel=3e-12, abs=0.)
    assert cells.joule_power_w_m3[0] == pytest.approx(
        8*(math.sqrt(2)-1)/(3*math.sqrt(math.pi))*1e-300, rel=3e-12, abs=0.)


def test_far_tail_asymptotics_and_truncated_energy_convergence():
    model = reference(conductivity_s_m=1., permeability_h_m=1.,
                      volumetric_heat_capacity_j_m3_k=1., boundary_field_rate_t_s=1.)
    coordinates = np.array([10., 20.])
    fields = model.fields(2*coordinates, 1.)
    b_leading_ratio = fields.magnetic_field_t * np.exp(coordinates**2) * math.sqrt(math.pi) * coordinates**3
    j_leading_ratio = fields.current_density_a_m2 * np.exp(coordinates**2) * math.sqrt(math.pi) * coordinates**2
    for ratio in (b_leading_ratio, j_leading_ratio):
        assert np.all(ratio < 1)
        assert 3.7 < (1-ratio[0])/(1-ratio[1]) < 4.1
    energy = model.energy_budget(1.)
    magnetic, heat = [], []
    for upper in (.5, 1., 2., 4., 10.):
        magnetic.append(integrate.quad(lambda u: direct_b(u)**2, 0, upper)[0])
        heat.append(16/5 * integrate.quad(lambda u: direct_g(u)**2, 0, upper)[0])
    assert np.all(np.diff(magnetic[:4]) > 0) and np.all(np.diff(heat[:4]) > 0)
    # The remaining positive tails are below spacing of the integrated totals.
    assert magnetic[-1] >= magnetic[-2] - 4*np.finfo(float).eps*magnetic[-1]
    assert heat[-1] >= heat[-2] - 4*np.finfo(float).eps*heat[-1]
    assert magnetic[-1] == pytest.approx(energy.magnetic_energy_j_m2, rel=2e-12)
    assert heat[-1] == pytest.approx(energy.joule_heat_j_m2, rel=2e-12)


def test_global_energy_constants_from_independent_space_and_time_integrals():
    model, time = reference(), 2e-6
    energy = model.energy_budget(time)
    ib = integrate.quad(lambda u: direct_b(u)**2, 0, 10, epsabs=1e-13, epsrel=1e-12)[0]
    ig = integrate.quad(lambda u: direct_g(u)**2, 0, 10, epsabs=1e-13, epsrel=1e-12)[0]
    common = (model.boundary_field_rate_t_s**2 / model.permeability_h_m
              * math.sqrt(model.magnetic_diffusivity_m2_s) * time**2.5)
    assert energy.magnetic_energy_j_m2 == pytest.approx(common * ib, rel=2e-12)
    time_integral = integrate.quad(lambda fraction: fraction**1.5, 0, 1,
                                   epsabs=1e-13, epsrel=1e-12)[0]
    assert energy.joule_heat_j_m2 == pytest.approx(8 * common * ig * time_integral, rel=2e-12)
    surface = model.fields(0., time)
    power = (float(surface.current_density_a_m2) / model.conductivity_s_m
             * model.boundary_field_rate_t_s * time / model.permeability_h_m)
    assert energy.surface_power_w_m2 == pytest.approx(power, rel=3e-14)
    assert energy.poynting_input_j_m2 == pytest.approx(power * time * time_integral, rel=2e-12)
    assert abs(energy.energy_residual_j_m2) <= 5e-14 * energy.poynting_input_j_m2
    assert energy.magnetic_energy_j_m2 + energy.joule_heat_j_m2 == pytest.approx(
        energy.poynting_input_j_m2, rel=5e-14)
    dt = time * 1e-5
    earlier, later = model.energy_budget(time - dt), model.energy_budget(time + dt)
    assert (later.poynting_input_j_m2 - earlier.poynting_input_j_m2) / (2*dt) == pytest.approx(power, rel=3e-10)
    instantaneous_joule = (8 * model.boundary_field_rate_t_s**2 / model.permeability_h_m
                           * math.sqrt(model.magnetic_diffusivity_m2_s) * time**1.5 * ig)
    assert (later.joule_heat_j_m2 - earlier.joule_heat_j_m2) / (2*dt) == pytest.approx(instantaneous_joule, rel=3e-10)


def test_finite_surface_temperature_and_independent_local_time_deposition():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    depths = 2 * scale * np.array([0., .03, .7, 1., 3.])
    temperature = model.adiabatic_temperature_rise(depths, time)
    wall = (2 * model.boundary_field_rate_t_s**2 * time**2
            / (model.permeability_h_m * math.pi * model.volumetric_heat_capacity_j_m3_k))
    assert temperature[0] == pytest.approx(wall, rel=1e-14)
    for index, u in enumerate(depths / (2*scale)):
        integral = integrate.quad(lambda fraction: fraction * direct_g(u / math.sqrt(fraction))**2,
                                  0, 1, epsabs=1e-18, epsrel=2e-11)[0]
        expected = (4 * model.boundary_field_rate_t_s**2 * time**2
                    / (model.permeability_h_m * model.volumetric_heat_capacity_j_m3_k) * integral)
        assert temperature[index] == pytest.approx(expected, rel=5e-10, abs=0.)
    dt = time * 1e-5
    derivative = (model.adiabatic_temperature_rise(depths, time + dt)
                  - model.adiabatic_temperature_rise(depths, time - dt)) / (2*dt)
    np.testing.assert_allclose(derivative,
                               model.fields(depths, time).joule_power_w_m3 / model.volumetric_heat_capacity_j_m3_k,
                               rtol=8e-9)


def test_similarity_and_drive_sign_scaling():
    model, time = reference(), 2e-6
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    depths = np.array([0., .5, 2., 7.]) * scale
    base, scaled = model.fields(depths, time), model.fields(3*depths, 9*time)
    np.testing.assert_allclose(scaled.magnetic_field_t, 9*base.magnetic_field_t, rtol=3e-13)
    np.testing.assert_allclose(scaled.current_density_a_m2, 3*base.current_density_a_m2, rtol=3e-13)
    np.testing.assert_allclose(scaled.joule_power_w_m3, 9*base.joule_power_w_m3, rtol=3e-13)
    np.testing.assert_allclose(model.adiabatic_temperature_rise(3*depths, 9*time),
                               81 * model.adiabatic_temperature_rise(depths, time), rtol=5e-12)
    assert model.energy_budget(9*time).joule_heat_j_m2 == pytest.approx(243*model.energy_budget(time).joule_heat_j_m2, rel=3e-13)
    negative = reference(boundary_field_rate_t_s=-model.boundary_field_rate_t_s)
    flipped = negative.fields(depths, time)
    np.testing.assert_array_equal(flipped.magnetic_field_t, -base.magnetic_field_t)
    np.testing.assert_array_equal(flipped.current_density_a_m2, -base.current_density_a_m2)
    np.testing.assert_array_equal(flipped.joule_power_w_m3, base.joule_power_w_m3)
    np.testing.assert_array_equal(negative.adiabatic_temperature_rise(depths, time),
                                  model.adiabatic_temperature_rise(depths, time))


def test_tail_fields_and_temperature_survive_unscaled_gaussian_underflow():
    model = reference(conductivity_s_m=1., permeability_h_m=1e-300,
                      volumetric_heat_capacity_j_m3_k=1., boundary_field_rate_t_s=1e300)
    time = 1e-200
    scale = math.sqrt(model.magnetic_diffusivity_m2_s * time)
    u = 30.
    assert math.exp(-u*u) == 0
    result = model.fields(2*scale*u, time)
    assert all(float(value) > 0 and math.isfinite(float(value))
               for value in (result.magnetic_field_t, result.current_density_a_m2, result.joule_power_w_m3))
    temperature = float(model.adiabatic_temperature_rise(2*scale*u, time))
    assert temperature > 0 and math.isfinite(temperature)
    # Independent positive, unscaled Gaussian-moment quadrature in w.
    moment = integrate.quad(lambda w: w * math.exp(-2*u*w-w*w), 0, np.inf,
                            epsabs=1e-16, epsrel=1e-12)[0]
    log_g = -u*u + math.log(2 / math.sqrt(math.pi)) + math.log(moment)
    expected_log_q = (math.log(4.) + 2*math.log(model.boundary_field_rate_t_s)
                      + math.log(time) - math.log(model.permeability_h_m) + 2*log_g)
    assert math.log(float(result.joule_power_w_m3)) == pytest.approx(expected_log_q, abs=1e-10)


def test_initial_state_and_zero_drive_are_finite_immutable_and_shape_preserving():
    model = reference()
    for driven, time in ((model, 0.), (reference(boundary_field_rate_t_s=0.), 1.)):
        values = driven.fields([[0., 1.], [2., 3.]], time)
        assert values.magnetic_field_t.shape == (2, 2)
        for array in (values.magnetic_field_t, values.current_density_a_m2, values.joule_power_w_m3,
                      driven.adiabatic_temperature_rise([0., 1.], time),
                      driven.cell_averages([0., 1., 2.], time).joule_power_w_m3):
            assert np.all(array == 0)
            with pytest.raises(ValueError):
                array.flags.writeable = True
        assert driven.energy_budget(time).poynting_input_j_m2 == 0
    with pytest.raises(FrozenInstanceError):
        model.boundary_field_rate_t_s = 1.


def test_small_time_compensated_drive_retains_finite_wall_temperature_and_energy():
    model = reference(conductivity_s_m=1., permeability_h_m=1.,
                      volumetric_heat_capacity_j_m3_k=1., boundary_field_rate_t_s=1e300)
    time = 1e-300
    assert time*time == 0
    assert float(model.adiabatic_temperature_rise(0., time)) == pytest.approx(2/math.pi, rel=3e-13)
    assert model.energy_budget(time).poynting_input_j_m2 > 0


@pytest.mark.parametrize("name", ["conductivity_s_m", "permeability_h_m", "volumetric_heat_capacity_j_m3_k"])
@pytest.mark.parametrize("value", [0, -1, True, 1+0j, "1", np.nan, np.inf])
def test_material_guards(name, value):
    with pytest.raises(ValueError):
        reference(**{name: value})


@pytest.mark.parametrize("value", [True, 1j, "1", np.inf, np.nan])
def test_drive_guards(value):
    with pytest.raises(ValueError):
        reference(boundary_field_rate_t_s=value)


@pytest.mark.parametrize("time", [-1, True, 1j, "1", np.inf, np.nan])
def test_time_guards(time):
    model = reference()
    for function, args in ((model.fields, ([0, 1], time)), (model.cell_averages, ([0, 1], time)),
                           (model.adiabatic_temperature_rise, ([0, 1], time)), (model.energy_budget, (time,))):
        with pytest.raises(ValueError):
            function(*args)


@pytest.mark.parametrize("depth", [[], [-1], [True], [1j], ["1"], [np.inf], [np.nan],
                                  np.ma.array([1.], mask=True)])
def test_depth_guards(depth):
    with pytest.raises(ValueError):
        reference().fields(depth, 1.)


@pytest.mark.parametrize("edges", [[0], [[0, 1]], [0, 0], [1, 0], [0, 1, 1]])
def test_cell_geometry_guards(edges):
    with pytest.raises(ValueError):
        reference().cell_averages(edges, 1.)


def test_parameters_have_no_defaults_and_unrepresentable_outputs_raise():
    assert all(value.default is inspect.Parameter.empty for value in inspect.signature(LinearRampConductivityHalfSpace).parameters.values())
    with pytest.raises(TypeError):
        LinearRampConductivityHalfSpace()
    with pytest.raises(ValueError):
        reference(conductivity_s_m=1e300, permeability_h_m=1e300)
    with pytest.raises(ValueError):
        reference(boundary_field_rate_t_s=1e308).energy_budget(1e308)
