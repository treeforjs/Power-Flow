import numpy as np
import pytest

from mhdlab.neutral_vacuum_map import packet_loss_times, observe_vacuum_packets


def request():
    return dict(birth_time_s=[0., .5], position_m=[[0., .1, 0.], [0., -.1, 0.]],
        velocity_m_s=[[0., 1., 30.], [0., -1., 50.]], particle_weight=[2., 3.],
        species_index=[0, 2], mass_kg=[1., 18.], time_s=[0., .5, 1., 2.],
        x_edges_m=[-2., -.5, .5, 2.], y_edges_m=[-2., -.1, .1, 2.],
        foil_bounds_m=[-.5, .5, -.1, .1], invariant_current_length_m=2.,
        initial_surface_particles=[2., 0., 3., 0.])


def test_exact_birth_times_and_volume_weighting():
    arrays, ledgers = observe_vacuum_packets(**request())
    assert ledgers[0]["species"][0]["inside_particles"] == 2
    assert ledgers[0]["species"][2]["launched_particles"] == 0
    assert ledgers[1]["species"][2]["inside_particles"] == 3
    assert ledgers[-1]["species"][0]["escaped_particles"] == 2
    assert np.all(arrays["number_density_m3"][:, :, 1, 1] == 0)


def test_frame_partition_does_not_change_trajectories_or_losses():
    data = request()
    _, coarse = observe_vacuum_packets(**data)
    data["time_s"] = np.arange(0, 2.01, .25)
    _, fine = observe_vacuum_packets(**data)
    assert coarse == [fine[i] for i in (0, 2, 4, 8)]


def test_solid_crossing_between_frames_is_absorbed_but_tangent_passes():
    exits, hits = packet_loss_times([[-1, 0, 0], [-1, .1, 0], [0, .1, 0]],
                                    [[10, 0, 0], [10, 0, 0], [0, 1, 0]],
                                    [-2, 2, -2, 2], [-.5, .5, -.1, .1])
    assert hits[0] == pytest.approx(.05)
    assert np.isinf(hits[1:]).all()


def test_current_velocity_is_retained_for_energy_without_false_z_escape():
    _, ledgers = observe_vacuum_packets(**request())
    assert ledgers[0]["species"][0]["launched_translational_energy_j"] == 901
    assert ledgers[0]["species"][0]["inside_particles"] == 2


def test_unaligned_solid_bins_and_overdrawn_inventory_reject():
    data = request()
    data["y_edges_m"] = [-2, 0, 2]
    with pytest.raises(ValueError, match="align"):
        observe_vacuum_packets(**data)
    data = request()
    data["initial_surface_particles"][0] = 1
    with pytest.raises(ValueError, match="exceed"):
        observe_vacuum_packets(**data)


def test_all_four_faces_birth_into_gas_bins_and_leave_without_false_absorption():
    data = request()
    data.update(birth_time_s=[0., 0., 0., 0.],
        position_m=[[-.5, 0., 0.], [.5, 0., 0.], [0., -.1, 0.], [0., .1, 0.]],
        velocity_m_s=[[-1., 0., 2.], [1., 0., -2.], [0., -1., 3.], [0., 1., -3.]],
        particle_weight=[2., 3., 4., 5.], species_index=[0, 1, 2, 3],
        mass_kg=[1., 2., 18., 17.], time_s=[0., .1], initial_surface_particles=[2., 3., 4., 5.])
    arrays, ledgers = observe_vacuum_packets(**data)
    assert arrays['number_density_m3'][0, 0, 1, 0] > 0
    assert arrays['number_density_m3'][0, 1, 1, 2] > 0
    assert arrays['number_density_m3'][0, 2, 0, 1] > 0
    assert arrays['number_density_m3'][0, 3, 2, 1] > 0
    assert np.all(arrays['number_density_m3'][:, :, 1, 1] == 0)
    for frame in ledgers:
        for species, weight in zip(frame['species'], [2., 3., 4., 5.]):
            assert species['inside_particles'] == weight
            assert species['absorbed_particles'] == 0
            assert species['particle_balance_residual'] == 0
