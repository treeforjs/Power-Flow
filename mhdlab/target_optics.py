"""Finite-pupil, conditional paraxial collection of isotropic volume emission.

The two reflective elements are unfolded into positive thin optical elements.
The emitting gas and opaque rectangular foil are invariant in the current
direction z. Source emissivity is integrated analytically along an explicitly
bounded z interval. No coating, grating efficiency, quantum efficiency, gain,
diffraction or measured camera calibration is inferred by this geometry model.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from numbers import Real, Integral
import math

import numpy as np


def _real(value, name, *, positive=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} requires finite real data")
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} has invalid finite range")
    return value


def _edges(value, name):
    array = np.asarray(value)
    if (np.ma.isMaskedArray(value) or array.ndim != 1 or len(array) < 2
            or array.dtype.kind not in 'iuf'):
        raise ValueError(f"{name} requires real increasing edges")
    array = np.asarray(array, dtype=float)
    if not np.isfinite(array).all() or np.any(np.diff(array) <= 0):
        raise ValueError(f"{name} requires finite increasing edges")
    return array


@dataclass(frozen=True)
class RelayGeometry:
    target_to_first_m: float
    mirror_separation_m: float
    second_to_slit_m: float
    focal_length_first_m: float
    focal_length_second_m: float
    diameter_first_m: float = .0254
    diameter_second_m: float = .0254
    slit_width_m: float = 50e-6
    slit_height_m: float = .0254
    spectrometer_f_number: float = 4.5
    optical_axis_y_m: float = 0.
    detector_magnification: float = 1.

    def __post_init__(self):
        for name, value in asdict(self).items():
            object.__setattr__(self, name, _real(value, name, positive=name != 'optical_axis_y_m'))

    def transfer_matrix(self, source_x_m=0.):
        """ABCD matrix from one physical source-x plane to the slit plane."""
        distance = self.target_to_first_m-_real(source_x_m, 'source_x_m')
        if distance <= 0:
            raise ValueError('source plane must be before the first mirror')
        def p(d): return np.array([[1., d], [0., 1.]])
        def lens(f): return np.array([[1., 0.], [-1/f, 1.]])
        return (p(self.second_to_slit_m) @ lens(self.focal_length_second_m)
                @ p(self.mirror_separation_m) @ lens(self.focal_length_first_m) @ p(distance))

    def diagnostics(self, source_x_m=0.):
        matrix = self.transfer_matrix(source_x_m)
        d = self.target_to_first_m-float(source_x_m)
        radius = self.diameter_first_m/2
        omega = 2*math.pi*(1-d/math.hypot(d, radius))
        return {'geometry': asdict(self), 'source_plane_x_m': float(source_x_m),
                'source_to_slit_abcd': matrix.tolist(),
                'source_plane_signed_magnification_A': float(matrix[0,0]),
                'source_plane_focus_coefficient_B_m': float(matrix[0,1]),
                'source_plane_is_paraxially_focused': bool(abs(matrix[0,1]) <= 1e-10),
                'magnification_interpretation': 'ABCD A is image magnification only at B=0; otherwise it is a transfer coefficient',
                'source_plane_matrix_determinant': float(np.linalg.det(matrix)),
                'on_axis_first_pupil_solid_angle_sr': omega,
                'on_axis_first_pupil_isotropic_fraction': omega/(4*math.pi),
                'on_axis_first_pupil_geometric_semiangle_rad': math.atan2(radius,d)}


def focused_relay_conjugate(target_to_first_m, mirror_separation_m,
                            focal_length_first_m, focal_length_second_m,
                            slit_width_m=50e-6):
    """Infer a slit conjugate for an unfolded, near-normal two-concave relay.

    Distances follow the ray encounters, not mechanical mirror names. A
    virtual image is returned as a diagnostic and is never assigned a physical
    object-plane slit width. Focus observation and uncertain distances must be
    reconciled independently; this routine does not fit measurements.
    """
    values=[_real(v,n,positive=True) for v,n in zip(
        [target_to_first_m,mirror_separation_m,focal_length_first_m,
         focal_length_second_m,slit_width_m],
        ['target_to_first_m','mirror_separation_m','focal_length_first_m',
         'focal_length_second_m','slit_width_m'])]
    u,d,f1,f2,width=values
    p=lambda length:np.array([[1.,length],[0.,1.]])
    power=lambda focal:np.array([[1.,0.],[-1/focal,1.]])
    pre=power(f2)@p(d)@power(f1)@p(u)
    if abs(pre[1,1])<64*np.finfo(float).eps:
        return dict(image_type='at_infinity',second_to_slit_m=None,
                    magnification_abs=None,object_plane_slit_width_m=None)
    final=-pre[0,1]/pre[1,1]
    system=p(final)@pre
    real=final>0
    return dict(image_type='real' if real else 'virtual',
        second_to_slit_m=float(final),formal_signed_magnification=float(system[0,0]),
        magnification_abs=float(abs(system[0,0])) if real else None,
        object_plane_slit_width_m=float(width/abs(system[0,0])) if real else None,
        source_to_slit_abcd=system.tolist(),focus_coefficient_B_m=float(system[0,1]),
        determinant=float(np.linalg.det(system)),
        qualification='Near-normal paraxial scenario; not an independently measured optical prescription')


def _foil_bounds(value):
    bounds = np.asarray(value)
    if (bounds.shape != (4,) or bounds.dtype.kind not in 'iuf'
            or not np.isfinite(bounds).all() or bounds[0] >= bounds[1] or bounds[2] >= bounds[3]):
        raise ValueError('foil_bounds_m must be finite [xmin,xmax,ymin,ymax]')
    return np.asarray(bounds, dtype=float)


def _blocked(sx, sy, ex, ey, foil):
    """Broadcast segment intersections with the open rectangle interior."""
    sx, sy, ex, ey = np.broadcast_arrays(sx, sy, ex, ey)
    enter, leave = np.full(sx.shape,-np.inf), np.full(sx.shape,np.inf)
    possible = np.ones(sx.shape,dtype=bool)
    for origin, end, lo, hi in ((sx,ex,foil[0],foil[1]),(sy,ey,foil[2],foil[3])):
        delta = end-origin
        moving = delta != 0
        possible &= moving | ((origin > lo) & (origin < hi))
        safe = np.where(moving,delta,1.)
        a,b = (lo-origin)/safe, (hi-origin)/safe
        enter = np.maximum(enter,np.where(moving,np.minimum(a,b),-np.inf))
        leave = np.minimum(leave,np.where(moving,np.maximum(a,b),np.inf))
    return possible & (np.minimum(leave,1.) > np.maximum(enter,0.))


def segment_intersects_foil(source_xy_m, destination_xy_m, foil_bounds_m):
    """Return open-interior foil interception, retaining exact tangencies."""
    a,b = np.asarray(source_xy_m), np.asarray(destination_xy_m)
    for values in (a,b):
        if (values.ndim < 1 or values.shape[-1] != 2 or values.dtype.kind not in 'iuf'
                or not np.isfinite(values).all()):
            raise ValueError('segment coordinates must be finite real (...,2) arrays')
    return _blocked(a[...,0],a[...,1],b[...,0],b[...,1],_foil_bounds(foil_bounds_m))


def _intersect_affine(lo, hi, a, b, half_width):
    """Intersect an interval in z with |a*z+b| <= half_width."""
    zero = a == 0
    safe = np.where(zero,1.,a)
    left,right = (-half_width-b)/safe, (half_width-b)/safe
    lower = np.where(zero,-np.inf,np.minimum(left,right))
    upper = np.where(zero,np.inf,np.maximum(left,right))
    invalid = zero & (np.abs(b) > half_width)
    return np.where(invalid,np.inf,np.maximum(lo,lower)), np.minimum(hi,upper)


def _solid_angle_z_integral(d, dy, pz, lo, hi):
    base = d*d+dy*dy
    def primitive(z):
        t = z-pz
        return d*t/(base*np.sqrt(base+t*t))
    return np.maximum(0.,primitive(hi)-primitive(lo))


def build_collection_kernel(x_edges_m, y_edges_m, detector_edges_m, foil_bounds_m,
                            source_z_interval_m, geometry, *, nradial=8, nazimuth=24,
                            opaque_foil=True):
    """Return (kernel, diagnostics) for a fixed conditional optical geometry.

    Kernel shape is (detector_row,y_cell,x_cell), with units m³. Multiplying it
    by an isotropic, ALL-DIRECTION volume photon emissivity [photons/(m³ s)]
    and summing source cells gives photons/s at the selected image rows, before
    spectral throughput. Do not apply another source-volume, z-length, etendue
    or 1/(4*pi) factor. Density alone is not a photon emissivity.

    Source-cell geometry is evaluated at its x/y center. The z integral is
    analytic; finite-pupil integration uses Gauss radial and uniform azimuth
    quadrature. Detector coordinates are slit-image coordinates times an
    explicitly supplied magnification, not a recovered iCMOS pixel mapping.
    The opaque foil is infinitely extruded in z, consistent with this invariant
    source model; end effects and actual three-dimensional alignment are absent.
    """
    x,y,rows = (_edges(value,name) for value,name in ((x_edges_m,'x_edges_m'),
        (y_edges_m,'y_edges_m'),(detector_edges_m,'detector_edges_m')))
    foil = _foil_bounds(foil_bounds_m)
    z = _edges(source_z_interval_m,'source_z_interval_m')
    if len(z) != 2:
        raise ValueError('source_z_interval_m requires exactly two endpoints')
    if not isinstance(geometry,RelayGeometry):
        raise ValueError('geometry must be RelayGeometry')
    for value,name,minimum in ((nradial,'nradial',1),(nazimuth,'nazimuth',4)):
        if isinstance(value,(bool,np.bool_)) or not isinstance(value,Integral) or value < minimum:
            raise ValueError(f'{name} requires an integer >= {minimum}')
    if type(opaque_foil) is not bool:
        raise ValueError('opaque_foil must be explicit bool')
    if x[-1] >= geometry.target_to_first_m:
        raise ValueError('source map must lie before the first mirror')
    xc,yc = .5*(x[:-1]+x[1:]), .5*(y[:-1]+y[1:])
    sx,sy = np.meshgrid(xc,yc)
    sx,sy = sx.ravel(),sy.ravel()
    area = (np.diff(y)[:,None]*np.diff(x)[None,:]).ravel()
    solid = ((sx > foil[0]) & (sx < foil[1]) & (sy > foil[2]) & (sy < foil[3]))
    radial,weight = np.polynomial.legendre.leggauss(int(nradial))
    radial,weight = .5*(radial+1), .5*weight
    theta = (np.arange(int(nazimuth))+.5)*2*np.pi/int(nazimuth)
    radius = geometry.diameter_first_m/2
    py = (radius*np.sqrt(radial)[:,None]*np.cos(theta)[None,:]).ravel()[None,:]
    pz = (radius*np.sqrt(radial)[:,None]*np.sin(theta)[None,:]).ravel()[None,:]
    pupil_weight = np.repeat(np.pi*radius*radius*weight/int(nazimuth),int(nazimuth))[None,:]
    kernel = np.zeros((len(rows)-1,len(sx)))
    first_volume = np.zeros(len(sx))
    blocked_samples = 0
    second_radius2 = (geometry.diameter_second_m/2)**2
    cone2 = (1/(2*geometry.spectrometer_f_number))**2
    for start in range(0,len(sx),256):
        stop = min(start+256,len(sx))
        xx, yy = sx[start:stop,None], sy[start:stop,None]
        d = geometry.target_to_first_m-xx
        dy = geometry.optical_axis_y_m+py-yy
        theta1y = dy/d-py/geometry.focal_length_first_m
        y2 = py+geometry.mirror_separation_m*theta1y
        theta2y = theta1y-y2/geometry.focal_length_second_m
        slit_y = y2+geometry.second_to_slit_m*theta2y
        image_y = geometry.detector_magnification*slit_y
        bin_index = np.searchsorted(rows,image_y,side='right')-1
        bin_index = np.minimum(bin_index,len(rows)-2)
        valid = ((np.abs(slit_y) <= geometry.slit_height_m/2) & (image_y >= rows[0])
                 & (image_y <= rows[-1]) & (y2*y2 <= second_radius2) & (theta2y*theta2y <= cone2))
        valid &= ~solid[start:stop,None]
        if opaque_foil:
            blocked = _blocked(xx,yy,geometry.target_to_first_m,geometry.optical_axis_y_m+py,foil)
            blocked_samples += int(np.count_nonzero(blocked & ~solid[start:stop,None]))
            valid &= ~blocked
        a1 = -1/d
        b1 = pz/d-pz/geometry.focal_length_first_m
        a2 = geometry.mirror_separation_m*a1
        b2 = pz+geometry.mirror_separation_m*b1
        at = a1-a2/geometry.focal_length_second_m
        bt = b1-b2/geometry.focal_length_second_m
        az = a2+geometry.second_to_slit_m*at
        bz = b2+geometry.second_to_slit_m*bt
        lo,hi = np.full(valid.shape,z[0]), np.full(valid.shape,z[1])
        lo,hi = _intersect_affine(lo,hi,az,bz,geometry.slit_width_m/2)
        lo,hi = _intersect_affine(lo,hi,a2,b2,np.sqrt(np.maximum(0.,second_radius2-y2*y2)))
        lo,hi = _intersect_affine(lo,hi,at,bt,np.sqrt(np.maximum(0.,cone2-theta2y*theta2y)))
        valid &= hi > lo
        # Invalid intervals never enter the primitive, avoiding inf/inf.
        low,high = np.where(valid,lo,0.), np.where(valid,hi,0.)
        weights = (_solid_angle_z_integral(d,dy,pz,low,high)*pupil_weight
                   *area[start:stop,None]/(4*np.pi))
        weights = np.where(valid,weights,0.)
        source_index = np.broadcast_to(np.arange(start,stop)[:,None],valid.shape)
        np.add.at(kernel,(bin_index[valid],source_index[valid]),weights[valid])
        unfiltered = (_solid_angle_z_integral(d,dy,pz,z[0],z[1])*pupil_weight
                      *area[start:stop,None]/(4*np.pi))
        first_volume[start:stop] = np.where(solid[start:stop],0.,np.sum(unfiltered,axis=1))
    volumes = area*(z[1]-z[0])
    total = kernel.sum(axis=0)
    if (not np.isfinite(kernel).all() or np.any(kernel < 0)
            or np.any(total > first_volume*(1+2e-12)+np.finfo(float).tiny)):
        raise ValueError('collection kernel violates finite nonnegative first-pupil support')
    diagnostics = geometry.diagnostics()
    diagnostics.update({'schema':'mhdlab-finite-pupil-collection-kernel-v1',
        'status':'conditional_geometry','kernel_units':'m3',
        'emissivity_units':'photons/(m3 s), isotropic all-direction total',
        'output_units_after_emissivity':'photons/s before spectral throughput',
        'kernel_shape':[len(rows)-1,len(y)-1,len(x)-1],
        'source_z_interval_m':z.tolist(),'detector_edges_m':rows.tolist(),
        'opaque_foil':opaque_foil,'foil_occlusion_model':'open rectangular x/y interior, invariant z',
        'metal_cells_excluded':True,'source_cell_geometry':'cell-center x/y; exact stated z integral',
        'pupil_quadrature':{'radial_gauss_order':int(nradial),'azimuth_midpoint_count':int(nazimuth),
                            'points_per_source_cell':int(nradial*nazimuth)},
        'blocked_gas_pupil_sample_count':blocked_samples,
        'first_pupil_volume_kernel_sum_m3':float(first_volume.sum()),
        'accepted_volume_kernel_sum_m3':float(total.sum()),
        'maximum_cell_geometric_collection_fraction':float(np.max(total/volumes)),
        'spectral_throughput_applied':False,'quantum_efficiency_applied':False,
        'camera_gain_applied':False,'physical_detector_prediction_qualified':False,
        'quadrature_convergence_established':False,
        'detector_coordinate_reference':'slit image coordinate times authored detector_magnification; no recovered iCMOS pixel mapping',
        'limitations':['Unfolded positive thin elements; paraxial slopes with exact first-pupil solid-angle measure.',
                      'No off-axis mirror astigmatism, diffraction, aberration or measured line-spread function.',
                      'Unspecified window clear aperture and camera active-area clipping are not modeled.',
                      'No calibrated mirror/window/grating/detector spectral transmission is inserted.',
                      'Source density requires a separately supplied CR photon-emission rate.']})
    return kernel.reshape(len(rows)-1,len(y)-1,len(x)-1),diagnostics
