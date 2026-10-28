"""Independent geometry, solid-angle, clipping and photon-unit checks."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from mhdlab.target_optics import (RelayGeometry, build_collection_kernel,
    segment_intersects_foil, _solid_angle_z_integral, _intersect_affine)


def geometry(**overrides):
    values=dict(target_to_first_m=.2,mirror_separation_m=.4,second_to_slit_m=.2,
        focal_length_first_m=.1,focal_length_second_m=.1,
        diameter_first_m=.02,diameter_second_m=2.,slit_width_m=2.,slit_height_m=2.,
        spectrometer_f_number=.1)
    values.update(overrides)
    return RelayGeometry(**values)


def test_focused_transfer_matrix_and_first_pupil_fraction():
    g=geometry()
    matrix=g.transfer_matrix()
    np.testing.assert_allclose(matrix,[[1.,0.],[20.,1.]],atol=3e-15)
    assert np.linalg.det(matrix)==pytest.approx(1.)
    report=g.diagnostics()
    expected=2*math.pi*(1-.2/math.sqrt(.2**2+.01**2))
    assert report['on_axis_first_pupil_solid_angle_sr']==pytest.approx(expected,rel=1e-13)
    assert report['on_axis_first_pupil_isotropic_fraction']==pytest.approx(expected/(4*math.pi))


def test_z_primitive_matches_independent_numerical_integral():
    d,dy,pz,lo,hi=.4826,.009,-.012,-.005,.004
    expected=quad(lambda z:d/(d*d+dy*dy+(pz-z)**2)**1.5,lo,hi,epsabs=0,epsrel=1e-12)[0]
    assert _solid_angle_z_integral(d,dy,pz,lo,hi)==pytest.approx(expected,rel=3e-15)


def test_open_downstream_train_matches_exact_on_axis_solid_angle():
    g=geometry()
    dx=dy=2e-6
    zspan=1e-10
    kernel,report=build_collection_kernel([-dx/2,dx/2],[-dy/2,dy/2],[-1.,1.],
        [.01,.02,.01,.02],[-zspan/2,zspan/2],g,nradial=12,nazimuth=32)
    omega=2*math.pi*(1-.2/math.sqrt(.2**2+.01**2))
    expected=dx*dy*zspan*omega/(4*math.pi)
    assert kernel.sum()==pytest.approx(expected,rel=2e-12)
    assert kernel.sum()<=dx*dy*zspan/2
    assert report['kernel_units']=='m3'
    assert report['spectral_throughput_applied'] is False


def test_affine_interval_positive_negative_and_constant_constraints():
    a=np.array([2.,-2.,0.,0.])
    b=np.array([.1,.1,0.,1.])
    lo,hi=_intersect_affine(np.full(4,-1.),np.full(4,1.),a,b,.05)
    np.testing.assert_allclose(lo[:3],[-.075,.025,-1.])
    np.testing.assert_allclose(hi[:3],[-.025,.075,1.])
    assert lo[3]>hi[3]


def test_opaque_foil_blocks_far_side_but_not_near_side_or_tangencies():
    foil=[-.0025,.0025,.0025,.0027]
    source=[[-.003,.0026],[.003,.0026],[-.003,.0027],[.0025,.0026]]
    dest=[[.48,.0026],[.48,.0026],[.48,.0027],[.48,.0026]]
    np.testing.assert_array_equal(segment_intersects_foil(source,dest,foil),[True,False,False,False])
    u=.4826;f=.1;q=f*u/(u-f)
    g=geometry(target_to_first_m=u,mirror_separation_m=2*q,second_to_slit_m=u,
        focal_length_first_m=f,focal_length_second_m=f,diameter_first_m=.001,
        diameter_second_m=.0254,slit_width_m=50e-6,slit_height_m=.0254,
        spectrometer_f_number=4.5)
    args=([-.003,-.0026,.0026,.003],[.00255,.00265],[-.01,.01],foil,[-.005,.005],g)
    blocked,_=build_collection_kernel(*args)
    open_kernel,_=build_collection_kernel(*args,opaque_foil=False)
    assert blocked[0,0,0]==0
    assert open_kernel[0,0,0]>0
    assert blocked[0,0,2]>0
    assert np.all(blocked[:,:,1]==0)  # metal does not emit even in control


def test_narrow_slit_reduces_collection_and_bins_conserve_accepted_photons():
    args=([-1e-5,1e-5],[-1e-5,1e-5],[-.1,0.,.1],[.01,.02,.01,.02],[-.005,.005])
    wide,_=build_collection_kernel(*args,geometry())
    narrow,_=build_collection_kernel(*args,geometry(slit_width_m=50e-6))
    assert 0<narrow.sum()<wide.sum()
    one,_=build_collection_kernel(args[0],args[1],[-.1,.1],*args[3:],geometry(slit_width_m=50e-6))
    assert one.sum()==pytest.approx(narrow.sum(),rel=2e-15)


@pytest.mark.parametrize('mutation',[{'target_to_first_m':0.},{'focal_length_first_m':np.nan},
    {'spectrometer_f_number':True},{'detector_magnification':0.}])
def test_invalid_geometry_rejects(mutation):
    with pytest.raises(ValueError): geometry(**mutation)


def test_invalid_domain_and_nonpositive_source_to_mirror_reject():
    args=([-.001,.001],[-.001,.001],[-.01,.01],[.01,.02,.01,.02],[-.005,.005],geometry())
    with pytest.raises(ValueError,match='exactly two'):
        build_collection_kernel(*args[:4],[-.005,0.,.005],args[-1])
    with pytest.raises(ValueError,match='before'):
        build_collection_kernel([.19,.21],*args[1:])
    with pytest.raises(ValueError,match='increasing'):
        build_collection_kernel(*args[:2],[0.,0.],*args[3:])
