import numpy as np
import pytest
from mhdlab.detector_emission import project_detector_photons
from mhdlab.resolved_emission_observation import ResolvedTransition

def request():
    return dict(time_s=[0.,1.,3.],upper_state_density_m3={'u':[[2.],[4.],[8.]]},
        collection_kernel_m3=[[3e-9],[6e-9]],transitions=[ResolvedTransition('line','u',120.,5.)],
        bin_edges_nm=[100.,120.,140.],gate_start_s=.5,gate_end_s=2.5,fwhm_nm=.1,
        throughput_by_transition={'line':.2},time_registration='authored clock',
        population_provenance='authored linear upper state',geometry_provenance='authored accepted volumes',
        throughput_provenance='authored fraction',profile_separability_reference='authored separable profile',
        response_reference_plane='authored intensifier entrance',status='manufactured')

def test_absolute_photon_budget_gate_geometry_and_throughput_enter_once():
    r=project_detector_photons(**request())
    # n(t)=2+2t: mean5 over[.5,2.5], duration2, A5, T.2; Gaussian half per bin.
    np.testing.assert_allclose(r['expected_photons_by_row_and_bin'],[[15e-9,15e-9],[30e-9,30e-9]],rtol=3e-15)
    assert r['ADU_conversion_performed'] is False

def test_shared_upper_state_branches_and_finite_spectral_window_are_not_renormalized():
    p=request();p['transitions'].append(ResolvedTransition('outside','u',180.,7.));p['throughput_by_transition']['outside']=.8
    r=project_detector_photons(**p)
    assert r['lines'][1]['fraction_in_wavelength_window']==0
    np.testing.assert_allclose(r['expected_photons_by_row_and_bin'],project_detector_photons(**request())['expected_photons_by_row_and_bin'])

@pytest.mark.parametrize('mutation',[
    lambda p:p.update(gate_start_s=-.1),
    lambda p:p.update(upper_state_density_m3={'wrong':[[2.],[4.],[8.]]}),
    lambda p:p.update(throughput_by_transition={}),
    lambda p:p.update(throughput_by_transition={'line':1.1}),
])
def test_extrapolation_missing_upper_states_and_ambiguous_response_reject(mutation):
    p=request();mutation(p)
    with pytest.raises(ValueError): project_detector_photons(**p)
