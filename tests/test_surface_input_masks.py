"""Missing surface data must not turn into physical source quantities."""

import numpy as np
import pytest

from mhdlab.thermal_release import BindingReservoir, SurfaceInventory, energy_budgets, release_step


@pytest.mark.parametrize("mask", [True, False])
def test_surface_inventory_requires_explicit_unmasked_values(mask):
    values = np.ma.array([[1e18, 9e30]], mask=[[False, mask]])
    with pytest.raises(ValueError, match="masked"):
        SurfaceInventory(values)


@pytest.mark.parametrize("operation", ["temperature", "released_energy"])
def test_masked_source_inputs_cannot_use_hidden_backing_values(operation):
    pool = BindingReservoir("manufactured", "M:ground", 1e12, .8, .6)
    with pytest.raises(ValueError, match="masked"):
        if operation == "temperature":
            release_step(SurfaceInventory([[1e18]]), [pool], np.ma.array([900.], mask=True), 1e-8)
        else:
            energy_budgets(np.ma.array([[1e18]], mask=True), [pool], [900.])
