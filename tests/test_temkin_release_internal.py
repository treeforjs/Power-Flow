"""Numerical edge contracts beyond the independent differential-equation tests.

Every parameter is manufactured. No species/material interpretation is implied.
"""

import inspect
import math

import numpy as np
import pytest

from mhdlab.constants import KB, QE
from mhdlab.temkin_release import (
    INVERSION_LOG_RESIDUAL_TOLERANCE, INVERSION_REGIMES, LAW_IDENTITY,
    TemkinReservoir, UnsupportedTemkinRangeError, energy_budgets, release_step,
    temkin_release_step,
)
from mhdlab.thermal_release import SurfaceInventory


def pool(**overrides):
    values = dict(name="manufactured", state_id="M:ground", attempt_frequency_s=1.,
                  activation_energy_ev=KB, desorption_energy_ev=0., internal_energy_ev=0.,
                  coverage_coefficient=.5, reference_coverage_m2=1., max_coverage_m2=1.,
                  kinetic_provenance="manufactured mathematical law",
                  energy_provenance="manufactured independent zero binding/internal costs")
    values.update(overrides)
    return TemkinReservoir(**values)


def test_no_defaults_and_identity_alias_and_diagnostic_dtypes():
    assert all(parameter.default is inspect.Parameter.empty
               for parameter in inspect.signature(TemkinReservoir).parameters.values())
    assert temkin_release_step is release_step
    trial = release_step(SurfaceInventory([[0., 1.]]), [pool()], [QE, QE], 2.)
    assert trial.law_identity == LAW_IDENTITY
    assert trial.inversion_regime.dtype == np.int8
    assert trial.inversion_iterations.dtype == np.int64
    assert trial.inversion_log_time_residual.dtype == np.float64
    assert trial.inversion_regime[0, 0] == 0
    assert trial.inversion_regime[0, 1] == 3
    assert abs(trial.inversion_log_time_residual[0, 1]) <= INVERSION_LOG_RESIDUAL_TOLERANCE
    for array in (trial.inversion_regime, trial.inversion_iterations, trial.inversion_log_time_residual):
        with pytest.raises(ValueError):
            array.flags.writeable = True
    with pytest.raises(TypeError):
        INVERSION_REGIMES[3] = "different model"


def test_huge_beta_with_huge_coupling_fails_instead_of_cancelling_barrier():
    reservoir = pool(activation_energy_ev=1e300, coverage_coefficient=1.)
    # The initial barrier is exactly zero, but the barrier immediately rises.
    # Freezing the initial rate would give an entirely different finite step.
    with pytest.raises(UnsupportedTemkinRangeError):
        release_step(SurfaceInventory([[1.]]), [reservoir], [1.], 1.)


def test_beta_beyond_float_range_is_safe_when_coverage_coupling_is_supported():
    reservoir = pool(activation_energy_ev=1e300, coverage_coefficient=1e-304)
    # beta is around 1e604 and alpha*N/Nref around 1e-604, so c is moderate.
    # The initial barrier remains enormous and the true release rounds to zero.
    trial = release_step(SurfaceInventory([[1e-300]]), [reservoir], [1e-300], 1e300)
    assert trial.trial_inventory.coverage_m2[0, 0] == 1e-300
    assert trial.released_m2[0, 0] == 0
    assert trial.inversion_regime[0, 0] == 2


def test_tiny_depletion_bound_at_largest_supported_coupling_keeps_release():
    reservoir = pool(activation_energy_ev=512 * KB, coverage_coefficient=1.,
                     reference_coverage_m2=1e300, max_coverage_m2=1e300)
    trial = release_step(SurfaceInventory([[1e300]]), [reservoir], [QE], 1e-20)
    # c*h_initial=512e-20 < epsilon. The difference from the initial rate
    # is below float64 relative precision, while N*h itself is representable.
    assert trial.released_m2[0, 0] == pytest.approx(1e280, rel=2e-13, abs=0.)
    assert trial.trial_inventory.coverage_m2[0, 0] == 1e300
    assert trial.inventory_balance_residual_m2[0, 0] == trial.released_m2[0, 0]
    assert trial.inversion_regime[0, 0] == 2


def test_float_max_inventory_residual_never_overflows_intermediate_sum():
    maximum = float(np.finfo(float).max)
    reservoir = pool(reference_coverage_m2=maximum, max_coverage_m2=maximum)
    trial = release_step(SurfaceInventory([[maximum]]), [reservoir], [QE], .5)
    assert math.isfinite(trial.inventory_balance_residual_m2[0, 0])
    assert 0 < trial.released_m2[0, 0] < maximum
    assert 0 < trial.trial_inventory.coverage_m2[0, 0] < maximum
    assert abs(trial.inventory_balance_residual_m2[0, 0]) <= 2 * np.spacing(maximum / 2)


def test_subnormal_remaining_population_and_strict_energy_domain():
    smallest = float(np.nextafter(0., 1.))
    reservoir = pool(reference_coverage_m2=smallest, max_coverage_m2=smallest)
    trial = release_step(SurfaceInventory([[smallest]]), [reservoir], [QE], 10.)
    assert trial.trial_inventory.coverage_m2[0, 0] == 0
    assert trial.released_m2[0, 0] == smallest
    assert trial.inventory_balance_residual_m2[0, 0] == 0
    with pytest.raises(ValueError, match="max_coverage"):
        energy_budgets([[2 * smallest]], [reservoir], [QE])
    with pytest.raises(ValueError, match="masked"):
        energy_budgets(np.ma.array([[smallest]], mask=True), [reservoir], [QE])
