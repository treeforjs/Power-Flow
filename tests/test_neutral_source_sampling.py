"""Unrepresentable packet tails are accounted without changing source physics."""
from decimal import Decimal
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from mhdlab.mfem_surface import load_mfem_surface_samples
from mhdlab.neutral_source_sampling import prepare_sampling_interval
from mhdlab.planar_release import PlanarSurfaceHistory
from mhdlab.temkin_planar import evolve_temkin_planar_release
from mhdlab.temkin_planar_boundary import sample_temkin_planar_emission_packets
from mhdlab.temkin_release import TemkinReservoir

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("_sampling_tail_demo", ROOT/'tools/export_mfem_surface.py')
demo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(demo)


@pytest.fixture
def source(tmp_path):
    surface, _ = load_mfem_surface_samples(demo.write_demo(tmp_path/'demo'))
    hot = PlanarSurfaceHistory([0., 1e-9], surface.x_edges_m, surface.z_edges_m,
        np.full((2,2,2), 1000.), surface.attributes)
    pools = tuple(TemkinReservoir(name=name, state_id=state, attempt_frequency_s=1e12,
        activation_energy_ev=0., desorption_energy_ev=.8, internal_energy_ev=0.,
        coverage_coefficient=0., reference_coverage_m2=1e16, max_coverage_m2=1e17,
        kinetic_provenance="Authored exhausted-tail test", energy_provenance="Authored independent energy debit")
        for name, state in [('ordinary','H2O:test'), ('tail','H:test')])
    coverage = np.stack((np.full((2,2),1e16), np.full((2,2),1e-320)))
    return evolve_temkin_planar_release(hot, coverage, pools, parameter_provenance="Authored packet underflow fixture")


def test_tail_representation_preserves_original_and_other_pool(source):
    before = source.coverage_m2.copy(), source.released_m2.copy()
    draws = np.random.default_rng(3).uniform(.001,.999,(128,7))
    face = source.surface_attributes['mfem_face']
    keywords = dict(face_normal_coordinate_m=face['normal_position_m'], outward_normal_sign=face['outward_normal_sign'])
    masses = {'H2O:test':3e-26, 'H:test':1.6e-27}
    with pytest.raises(ValueError, match="below representable particle-weight"):
        sample_temkin_planar_emission_packets(source, masses, 0, draws, **keywords)
    view, index, report = prepare_sampling_interval(source, 0)
    assert index == 0 and report['omitted_component_count'] == 4
    assert report['physical_source_history_changed'] is False
    assert Decimal(report['omitted_particles_exact_decimal']) > 0
    assert Decimal(report['omitted_particles_exact_decimal']) < Decimal(report['omitted_particles_upper_bound_exact_decimal'])
    assert report['omitted_particles_upper_bound'] == 4e-250
    np.testing.assert_array_equal(source.coverage_m2, before[0])
    np.testing.assert_array_equal(source.released_m2, before[1])
    np.testing.assert_array_equal(view.released_m2[0,0], source.released_m2[0,0])
    np.testing.assert_array_equal(view.translational_energy_j_m2[0,0], source.translational_energy_j_m2[0,0])
    packets = sample_temkin_planar_emission_packets(view, masses, index, draws, **keywords)
    assert np.isfinite(packets.particle_weight).all() and (packets.particle_weight > 0).all()
    assert set(packets.state_ids) == {'H2O:test'}


def test_all_filtered_returns_no_sampling_interval(source):
    view, index, report = prepare_sampling_interval(source, 0, minimum_component_particles=1e20)
    assert view is None and index == 0
    assert report['omitted_component_count'] == 8


def test_no_filtered_component_keeps_original_object(source):
    # An ordinary-only law record contains no exhausted tail.
    surface = PlanarSurfaceHistory(source.time_edges_s, source.x_edges_m, source.z_edges_m,
        np.stack((source.temperature_k[0],source.temperature_k[0])), source.surface_attributes)
    ordinary = evolve_temkin_planar_release(surface, source.coverage_m2[0,:1], source.reservoirs[:1],
        parameter_provenance=source.parameter_provenance)
    view, index, report = prepare_sampling_interval(ordinary, 0)
    assert view is ordinary and index == 0
    assert report['omitted_component_count'] == 0
    assert report['omitted_particles_upper_bound'] == 0


@pytest.mark.parametrize('cutoff',[0.,-1.,float('nan'),float('inf'),True])
def test_invalid_cutoff_rejects(source, cutoff):
    with pytest.raises(ValueError, match='finite and positive'):
        prepare_sampling_interval(source, 0, cutoff)
