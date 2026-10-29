from fractions import Fraction
import numpy as np
import pytest
from mhdlab.target_optics import RelayGeometry,focused_relay_conjugate

def test_known_symmetric_relay_conjugate():
    r=focused_relay_conjugate(.2,.4,.1,.1)
    assert r['image_type']=='real'
    assert r['second_to_slit_m']==pytest.approx(.2)
    assert r['magnification_abs']==pytest.approx(1)
    assert r['object_plane_slit_width_m']==pytest.approx(50e-6)

def test_focused_scenario_reproduces_existing_ray_matrix():
    u,d,f=.4826,.2794,.1;r=focused_relay_conjugate(u,d,f,f)
    g=RelayGeometry(u,d,r['second_to_slit_m'],f,f)
    np.testing.assert_allclose(r['source_to_slit_abcd'],g.transfer_matrix(),atol=1e-13)
    assert abs(r['focus_coefficient_B_m'])<1e-14
    assert r['determinant']==pytest.approx(1)
    uf,df,ff=Fraction(4826,10000),Fraction(2794,10000),Fraction(1,10)
    q1=1/(1/ff-1/uf);p2=df-q1;q2=1/(1/ff-1/p2)
    assert r['magnification_abs']==pytest.approx(float(abs(q1*q2/(uf*p2))),rel=1e-13)

def test_approximate_eight_inch_case_is_not_a_physical_slit_image():
    r=focused_relay_conjugate(.4826,.2032,.1,.1)
    assert r['image_type']=='virtual'
    assert r['second_to_slit_m']<0
    assert r['magnification_abs'] is None
    assert r['object_plane_slit_width_m'] is None

def test_first_mirror_collimates_without_division_by_zero():
    r=focused_relay_conjugate(.1,.25,.1,.1)
    assert r['image_type']=='real'
    assert r['second_to_slit_m']==pytest.approx(.1)

@pytest.mark.parametrize('index,value',[(0,0),(1,-1),(2,float('nan')),(3,True)])
def test_invalid_measurement_rejects(index,value):
    args=[.4826,.2794,.1,.1];args[index]=value
    with pytest.raises(ValueError):focused_relay_conjugate(*args)
