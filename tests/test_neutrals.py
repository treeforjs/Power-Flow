import numpy as np

from mhdlab.constants import SPECIES_MASS_KG
from mhdlab.geometry import Geometry, boundary_mask
from mhdlab.neutrals import VelocityGrid, half_range_flux_moment, shift_fractional_no_wrap, thermal_weights
from mhdlab.neutrals import KineticNeutralSolver


def test_half_range_flux_moment_positive_and_scales_with_temperature():
    mass = SPECIES_MASS_KG["H2O"]
    cold = half_range_flux_moment(300.0, mass)
    hot = half_range_flux_moment(1200.0, mass)
    assert cold > 0.0
    assert np.isclose(hot / cold, 2.0, rtol=1e-12)


def test_thermal_velocity_weights_are_normalized():
    grid = VelocityGrid.polar(max_speed_m_s=5000.0, n_speed=4, n_angle=8)
    weights = thermal_weights(grid, SPECIES_MASS_KG["H2O"], 600.0)
    assert np.isclose(weights.sum(), 1.0)
    assert np.all(weights >= 0.0)


def test_surface_emission_deposits_neutrals_into_neighboring_vacuum():
    geom = Geometry.from_json("examples/geometry/mykonos_foil_gap_5mm.json")
    raster = geom.rasterize(nx=40, ny=48)
    grid = VelocityGrid.polar(max_speed_m_s=5000.0, n_speed=3, n_angle=6)
    solver = KineticNeutralSolver(raster, ["H2O"], grid, {"wall_sticking": 0.0})
    state = solver.initial_state()
    material = raster.mask_by_kind("cathode") | raster.mask_by_kind("anode")
    surface = boundary_mask(material, raster.mask_by_kind("vacuum") | ~material)
    source = np.zeros(raster.shape)
    source[surface] = 1.0e20
    temperature = np.full(raster.shape, 600.0)

    next_state = solver.step(
        state,
        dt_s=1.0e-9,
        surface_source_m2_s=source,
        surface_temperature_k=temperature,
        source_species="H2O",
        surface_mask=surface,
    )

    density = next_state.density("H2O")
    assert density.sum() > 0.0
    assert np.all(density[~solver.vacuum_mask] == 0.0)


def test_fractional_shift_moves_density_without_rounding_to_zero():
    arr = np.zeros((5, 5))
    arr[2, 2] = 1.0

    shifted = shift_fractional_no_wrap(arr, sx=0.5, sy=0.0)

    assert np.isclose(shifted.sum(), 1.0)
    assert np.isclose(shifted[2, 2], 0.5)
    assert np.isclose(shifted[2, 3], 0.5)
