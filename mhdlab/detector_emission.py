"""Project supplied CR upper-state populations into finite-pupil photon counts.

This evolves no populations and supplies no excitation rates, throughput or
line widths. Geometry already contains accepted volume and isotropic solid
angle; no additional optical depth, volume, etendue or 4pi factor is inserted.
"""
from __future__ import annotations
from collections.abc import Mapping
from numbers import Real
import numpy as np
from scipy.special import ndtr
from .resolved_emission_observation import ResolvedTransition


def _array(value,name,rank):
    a=np.asarray(value)
    if a.ndim!=rank or a.dtype.kind not in 'iuf' or not a.size:
        raise ValueError(f'{name} requires real rank{rank} data')
    a=np.array(a,dtype=float,copy=True)
    if not np.isfinite(a).all(): raise ValueError(f'{name} must be finite')
    return a


def project_detector_photons(time_s, upper_state_density_m3, collection_kernel_m3,
        transitions, bin_edges_nm, gate_start_s, gate_end_s, fwhm_nm,
        throughput_by_transition, *, time_registration, population_provenance,
        geometry_provenance, throughput_provenance, profile_separability_reference,
        response_reference_plane, status):
    """Return expected photons per spatial row and wavelength bin in one gate.

    Kernel is (row,cell); each upper population is (time,cell). Transition
    throughput is an explicit fraction in [0,1] into the declared detector
    reference plane. The caller justifies separable spectral redistribution
    after spatial collection; a measured joint response needs a joint kernel.
    ADU, intensifier amplification and measurement covariance are separate.
    """
    for name,value in locals().copy().items():
        if name in {'time_registration','population_provenance','geometry_provenance',
                'throughput_provenance','profile_separability_reference','response_reference_plane'}:
            if not isinstance(value,str) or not value.strip(): raise ValueError(f'{name} requires evidence')
    if status not in {'manufactured','source_evaluated','owner_supplied'}: raise ValueError('explicit population status required')
    t=_array(time_s,'time_s',1);k=_array(collection_kernel_m3,'collection_kernel_m3',2)
    edges=_array(bin_edges_nm,'bin_edges_nm',1)
    if len(t)<2 or np.any(np.diff(t)<=0) or np.any(k<0) or len(edges)<2 or np.any(edges<=0) or np.any(np.diff(edges)<=0):
        raise ValueError('increasing time/wavelength axes and nonnegative kernel required')
    if (any(isinstance(v,(bool,np.bool_)) or not isinstance(v,Real) for v in (gate_start_s,gate_end_s))
            or not np.isfinite([gate_start_s,gate_end_s]).all()
            or not t[0]<=gate_start_s<gate_end_s<=t[-1]):
        raise ValueError('gate must be contained in supplied populations; no extrapolation')
    if (not isinstance(transitions,(list,tuple)) or not transitions
            or any(not isinstance(v,ResolvedTransition) for v in transitions)
            or len({v.transition_id for v in transitions})!=len(transitions)):
        raise ValueError('distinct explicitly resolved transitions required')
    if not isinstance(upper_state_density_m3,Mapping) or not isinstance(throughput_by_transition,Mapping):
        raise ValueError('state populations and transition throughput require mappings')
    if set(throughput_by_transition)!={v.transition_id for v in transitions}:
        raise ValueError('throughput must identify every exact transition, without missing or extra factors')
    widths=np.asarray(fwhm_nm)
    if widths.ndim==0: widths=np.full(len(transitions),widths.item())
    widths=_array(widths,'fwhm_nm',1)
    if len(widths)!=len(transitions) or np.any(widths<=0): raise ValueError('positive supplied width per transition required')
    knots=np.r_[gate_start_s,t[(t>gate_start_s)&(t<gate_end_s)],gate_end_s]
    indices=np.clip(np.searchsorted(t,knots,side='right')-1,0,len(t)-2)
    fractions=(knots-t[indices])/(t[indices+1]-t[indices])
    duration=float(gate_end_s-gate_start_s);means={}
    for state in {v.upper_state_id for v in transitions}:
        if state not in upper_state_density_m3: raise ValueError('missing exact upper-state population: '+state)
        n=_array(upper_state_density_m3[state],state,2)
        if n.shape!=(len(t),k.shape[1]) or np.any(n<0): raise ValueError('population requires nonnegative(time,cell) densities')
        samples=(1-fractions[:,None])*n[indices]+fractions[:,None]*n[indices+1]
        means[state]=np.einsum('t,tc->c',np.diff(knots)/duration,.5*(samples[:-1]+samples[1:]))
    photons=np.zeros((len(k),len(edges)-1));lines=[]
    with np.errstate(over='raise',invalid='raise'):
        for transition,width in zip(transitions,widths):
            raw=throughput_by_transition[transition.transition_id]
            if isinstance(raw,(bool,np.bool_)) or not isinstance(raw,Real) or not np.isfinite(raw) or not 0<=raw<=1:
                raise ValueError('optical throughput must be an explicit photon fraction in[0,1]')
            sigma=width/np.sqrt(8*np.log(2));lo=(edges[:-1]-transition.vacuum_wavelength_nm)/sigma;hi=(edges[1:]-transition.vacuum_wavelength_nm)/sigma
            fraction=np.where(lo>0,ndtr(-lo)-ndtr(-hi),ndtr(hi)-ndtr(lo))
            all_wavelength=duration*transition.einstein_a_per_s*float(raw)*(k@means[transition.upper_state_id])
            photons+=all_wavelength[:,None]*fraction[None]
            lines.append({'transition_id':transition.transition_id,'upper_state_id':transition.upper_state_id,
                'vacuum_wavelength_nm':transition.vacuum_wavelength_nm,'einstein_a_per_s':transition.einstein_a_per_s,
                'photon_throughput':float(raw),'fraction_in_wavelength_window':float(fraction.sum()),
                'all_wavelength_photons_by_row':all_wavelength.tolist()})
    if not np.isfinite(photons).all(): raise ValueError('photon counts exceed finite numeric range')
    return {'schema':'mhdlab-state-resolved-detector-photons-v1','status':status,
        'reference_plane':response_reference_plane,'expected_photons_by_row_and_bin':photons.tolist(),
        'gate_s':[float(gate_start_s),float(gate_end_s)],'bin_edges_nm':edges.tolist(),'lines':lines,
        'factor_ownership':'Collection kernel already contains accepted volume and1/(4pi); A, gate duration, photon throughput and spectral-bin fraction applied exactly once.',
        'population_evolution_performed':False,'ADU_conversion_performed':False,'physical_fit_performed':False,
        'time_registration':time_registration,'population_provenance':population_provenance,
        'geometry_provenance':geometry_provenance,'throughput_provenance':throughput_provenance,
        'profile_separability_reference':profile_separability_reference,
        'spectral_profile':'unit-area photon Gaussian; finite-window losses not renormalized',
        'gate_interpolation':'piecewise linear supplied populations; exact integration of that interpolant'}
