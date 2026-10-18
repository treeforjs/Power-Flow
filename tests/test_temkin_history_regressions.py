"""Conservative history validation at float64 inventory limits."""

from dataclasses import replace
import math

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.flash_surface import FlashSurfaceHistory
from mhdlab.temkin_history import evolve_temkin_release
from mhdlab.temkin_release import TemkinReservoir
from mhdlab.thermal_release import SurfaceInventory


def make_history(initial, *, attempt=1., barrier=KB, alpha=.5, temperature=QE, dt=.5):
    pool = TemkinReservoir("manufactured", "M:ground", attempt, barrier, 0., 0.,
                           alpha, initial, initial, "manufactured rate", "manufactured energy")
    surface = FlashSurfaceHistory([0., dt], [0., 1.], [[temperature], [temperature]],
                                  "synthetic", "regression", "manufactured", "stationary face")
    return evolve_temkin_release(surface, SurfaceInventory([[initial]]), [pool],
                                 parameter_provenance="numerical limit verification")


def test_compensated_balance_is_retained_below_the_remaining_inventory_spacing():
    history = make_history(1e300, attempt=1e-200, barrier=0., alpha=0.,
                           temperature=300., dt=1e-200)
    assert history.coverage_m2[-1, 0, 0] == 1e300
    assert history.released_m2[0, 0, 0] == pytest.approx(1e-100, rel=3e-13, abs=0.)
    assert history.inventory_balance_residual_m2[0, 0, 0] == history.released_m2[0, 0, 0]
    with pytest.raises(ValueError, match="compensated"):
        replace(history, inventory_balance_residual_m2=np.zeros_like(history.released_m2))


def test_maximum_inventory_balance_does_not_overflow_an_intermediate_sum():
    maximum = float(np.finfo(float).max)
    history = make_history(maximum)
    after = float(history.coverage_m2[-1, 0, 0])
    released = float(history.released_m2[0, 0, 0])
    expected = math.fsum((after, -maximum, released))
    assert math.isfinite(expected)
    assert history.inventory_balance_residual_m2[0, 0, 0] == expected


def test_survivor_underflow_regime_cannot_retain_a_positive_survivor():
    history = make_history(1., attempt=1e300, dt=1.)
    assert history.inversion_regime[0, 0, 0] == 4
    assert history.coverage_m2[-1, 0, 0] == 0.
    coverage = history.coverage_m2.copy()
    smallest = float(np.nextafter(0., 1.))
    coverage[-1, 0, 0] = smallest
    with pytest.raises(ValueError, match="survivor"):
        replace(history, coverage_m2=coverage,
                inventory_balance_residual_m2=np.full_like(history.released_m2, smallest))


def test_a_nonlinear_root_cannot_be_relabelled_as_the_tiny_depletion_limit():
    history = make_history(1.)
    assert history.inversion_regime[0, 0, 0] == 3
    with pytest.raises(ValueError, match="regime"):
        replace(history, inversion_regime=np.full_like(history.inversion_regime, 2),
                inversion_iterations=np.zeros_like(history.inversion_iterations),
                inversion_log_time_residual=np.zeros_like(history.inversion_log_time_residual))
