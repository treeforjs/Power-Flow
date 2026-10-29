"""Unit conversion, Richardson slopes and the absence of a hidden current cap."""
import math
import numpy as np
import pytest
from mhdlab.constants import KB,QE
from mhdlab.thermionic_supply import IDEAL_RICHARDSON_A_M2_K2 as A0,log_thermionic_supply,thermionic_supply


def test_ideal_constant_is120_ampere_per_square_centimetre_per_kelvin_squared():
    assert A0/1e4==pytest.approx(120.173,rel=1e-5)


def test_work_function_sets_richardson_plot_slope_and_prefactor_sets_scale():
    t=np.array([1000.,1500.,2000.,2500.]);j=thermionic_supply(t,work_function_ev=4.3,richardson_a_m2_k2=A0)
    slope=np.polyfit(1/t,np.log(j/t**2),1)[0]
    assert slope==pytest.approx(-4.3*QE/KB,rel=1e-13)
    np.testing.assert_allclose(thermionic_supply(t,work_function_ev=4.3,richardson_a_m2_k2=.1*A0),.1*j,rtol=1e-13)


def test_supply_is_monotone_and_does_not_saturate_at_an_assumed_space_charge_limit():
    t=np.array([1500.,2000.,3000.,5000.,10000.]);j=thermionic_supply(t,work_function_ev=4.5,richardson_a_m2_k2=A0)
    assert np.all(np.diff(j)>0) and j[-1]>1e10
    assert j[1]==pytest.approx(22.0,rel=5e-4)


@pytest.mark.parametrize('temperature',[0.,-1.,np.nan,np.inf,True,'1000'])
def test_invalid_surface_temperatures_are_rejected(temperature):
    with pytest.raises(ValueError):thermionic_supply(temperature,work_function_ev=4.3,richardson_a_m2_k2=A0)
