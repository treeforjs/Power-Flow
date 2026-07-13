import numpy as np

from mhdlab.constants import SPECIES_MASS_KG
from mhdlab.neutrals import VelocityGrid, half_range_flux_moment, thermal_weights


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
