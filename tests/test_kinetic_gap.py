"""Analytic and conservation checks with explicitly manufactured channel data."""

import math

import numpy as np
import pytest

from mhdlab.kinetic_gap import GapState, Particle, advance_gap, cell_moments, sample_thermal_emission, total_energy_j
from mhdlab.state_kinetics import (
    ELEMENTARY_CHARGE_C as Q,
    CrossSectionDomainError,
    InternalState,
    ReactionChannel,
    Species,
    StateCatalog,
)


ME = 9.1093837139e-31
MP = 3.0e-26


def fixture():
    catalog = StateCatalog(
        [Species("e", ME, -1, {}), Species("P", MP, 0, {"H": 2, "O": 1}),
         Species("P+", MP - ME, 1, {"H": 2, "O": 1})],
        [InternalState("e0", "e", 0, {}),
         InternalState("p0", "P", 0, {"v": [0, 0, 0], "J": 0}),
         InternalState("p1", "P", 1, {"v": [1, 0, 0], "J": 0}),
         InternalState("pi", "P+", 2, {"v": [0, 0, 0], "N": 0})],
    )
    speed = math.sqrt(2 * 5 * Q / ME)
    state = GapState((Particle("e0", (0, 0, 0), (ME * speed, 0, 0)),
                      Particle("p0", (0, 0, 0), (0, 0, 0))))
    channel = ReactionChannel(
        "analytic_excitation", ("e0", "p0"), ("e0", "p1"),
        (0, 1, 2, 100), (0, 0, 1e-18, 1e-18), "jacobi_energy_partition",
        "Manufactured integral cross section for an analytic verification problem",
        metadata={"kinematics_provenance": "Explicit deterministic direction for this verification problem"},
    )
    return catalog, state, channel, speed


def kernel(channel, catalog, incoming, rng):
    return [1.0], [[0.0, 1.0, 0.0]]


def advance(state, catalog, **kw):
    options = dict(collision_volume_m3=1.0, electric_field_v_m=[0, 0, 0],
                   magnetic_field_t=[0, 0, 0], rng=np.random.default_rng(3))
    options.update(kw)
    return advance_gap(state, catalog, options.pop("dt_s", 1e-10), **options)


def test_ballistic_streaming_and_exact_state_identity():
    catalog, state, _, speed = fixture()
    result = advance(state, catalog)
    assert result.state.time_s == 1e-10
    assert result.state.particles[0].position_m[0] == pytest.approx(speed * 1e-10)
    assert result.state.particles[0].state_id == "e0"
    assert result.total_energy_residual_j == 0
    assert state.particles[0].position_m == (0, 0, 0)


def test_electric_work_and_state_weighted_current():
    catalog, _, _, _ = fixture()
    state = GapState((Particle("e0", (0, 0, 0), (0, 0, 0), 10),))
    dt, field = 1e-12, 1000.0
    result = advance(state, catalog, dt_s=dt, electric_field_v_m=[field, 0, 0])
    expected_p = -Q * field * dt
    expected_gain = 10 * expected_p**2 / (2 * ME)
    assert result.state.particles[0].momentum_kg_m_s[0] == pytest.approx(expected_p, rel=1e-12, abs=0)
    assert result.field_kinetic_gain_j == pytest.approx(expected_gain, rel=1e-12, abs=0)
    assert abs(result.total_energy_residual_j) < 1e-14 * expected_gain
    assert abs(result.field_work_quadrature_residual_j) < 1e-14 * expected_gain
    moments = cell_moments(result.state, catalog, 2.0)
    assert moments["charge_density_c_m3"] == pytest.approx(-5 * Q, rel=1e-12, abs=0)
    assert moments["current_density_a_m2"][0] == pytest.approx(-5 * Q * expected_p / ME, rel=1e-12, abs=0)


def test_magnetic_subcycling_preserves_energy():
    catalog, state, _, _ = fixture()
    result = advance(state, catalog, dt_s=1e-10, magnetic_field_t=[0, 0, 1])
    assert result.orbit_substeps > 1
    assert abs(result.field_kinetic_gain_j) < 2e-13 * total_energy_j(state, catalog)
    assert abs(result.total_energy_residual_j) < 2e-13 * total_energy_j(state, catalog)
    assert result.electric_displacement_work_j == 0


def test_excitation_changes_internal_state_and_recoil_together():
    catalog, state, channel, _ = fixture()
    result = advance(state, catalog, collision_volume_m3=1e-28,
                     channels=(channel,), differential_kernels={channel.name: kernel})
    assert result.reaction_counts[channel.name] == 1
    assert [p.state_id for p in result.state.particles] == ["e0", "p1"]
    assert result.state.particles[0].momentum_kg_m_s[1] > 0
    assert result.state.particles[1].momentum_kg_m_s[1] < 0
    assert abs(result.total_energy_residual_j) < 1e-14 * total_energy_j(state, catalog)
    assert abs(result.reaction_energy_residual_j) < 1e-14 * Q
    np.testing.assert_allclose(result.reaction_momentum_residual_kg_m_s, 0, atol=1e-39)


def test_gillespie_first_event_matches_analytic_survival_probability():
    catalog, state, channel, speed = fixture()
    dt = 1e-10
    # Set lambda*dt=0.8 for one constant-hazard pair. The first event removes
    # the only p0 donor, so the exact number of events is Bernoulli.
    volume = 1e-18 * speed * dt / 0.8
    rng = np.random.default_rng(718)
    events = sum(advance(state, catalog, dt_s=dt, collision_volume_m3=volume,
                         channels=(channel,), differential_kernels={channel.name: kernel},
                         rng=rng).reaction_counts[channel.name] for _ in range(600))
    probability = 1 - math.exp(-0.8)
    standard_deviation = math.sqrt(600 * probability * (1 - probability))
    assert abs(events - 600 * probability) < 5 * standard_deviation


def test_ionization_creates_charged_products_without_charge_or_energy_loss():
    catalog, state, _, _ = fixture()
    ionization = ReactionChannel(
        "analytic_ionization", ("e0", "p0"), ("e0", "e0", "pi"),
        (0, 2, 3, 100), (0, 0, 1e-18, 1e-18), "jacobi_energy_partition",
        "Manufactured ionization verification channel",
        metadata={"kinematics_provenance": "Specified Jacobi fractions and orthogonal directions"},
    )
    def ion_kernel(channel, catalog, incoming, rng):
        return [.3, .7], [[1, 0, 0], [0, 1, 0]]
    result = advance(state, catalog, collision_volume_m3=1e-28,
                     channels=(ionization,), differential_kernels={ionization.name: ion_kernel})
    assert len(result.state.particles) == 3
    assert result.reaction_counts[ionization.name] == 1
    assert cell_moments(result.state, catalog, 1)["charge_density_c_m3"] == pytest.approx(-Q, rel=1e-12, abs=0)
    assert abs(result.total_energy_residual_j) < 1e-14 * total_energy_j(state, catalog)


def test_missing_cross_section_support_is_rejected():
    catalog, state, channel, _ = fixture()
    fast = Particle("e0", (0, 0, 0), (ME * math.sqrt(2 * 200 * Q / ME), 0, 0))
    with pytest.raises(CrossSectionDomainError):
        advance(GapState((fast, state.particles[1])), catalog,
                channels=(channel,), differential_kernels={channel.name: kernel})


def test_missing_differential_data_are_rejected():
    catalog, state, channel, _ = fixture()
    with pytest.raises(ValueError, match="missing differential kernel"):
        advance(state, catalog, channels=(channel,))


def test_unequal_weights_require_explicit_splitting():
    catalog, state, _, _ = fixture()
    unequal = GapState((state.particles[0], Particle("p0", (0, 0, 0), (0, 0, 0), 2)))
    with pytest.raises(ValueError, match="equal macro-particle"):
        advance(unequal, catalog)


def test_subcycling_limit_refuses_unresolved_orbit():
    catalog, state, _, _ = fixture()
    with pytest.raises(ValueError, match="subcycling"):
        advance(state, catalog, magnetic_field_t=[0, 0, 1000], max_orbit_substeps=2)


def test_event_failure_preserves_input_particle_state():
    catalog, state, channel, _ = fixture()
    def invalid_kernel(*args):
        return [1.0], [[0, 0, 0]]
    before = total_energy_j(state, catalog)
    rng = np.random.default_rng(137)
    random_state = rng.bit_generator.state
    with pytest.raises(ValueError):
        advance(state, catalog, collision_volume_m3=1e-28,
                channels=(channel,), differential_kernels={channel.name: invalid_kernel}, rng=rng)
    assert total_energy_j(state, catalog) == before
    assert state.time_s == 0
    assert state.particles[1].state_id == "p0"
    assert rng.bit_generator.state == random_state


def test_empty_cell_is_well_defined():
    catalog, _, _, _ = fixture()
    result = advance(GapState(()), catalog)
    assert result.state.particles == ()
    assert result.total_energy_residual_j == 0
    assert cell_moments(result.state, catalog, 1)["state_density_m3"] == {}


def test_thermal_launch_has_outward_flux_energy_and_supplied_internal_states():
    catalog, _, _, _ = fixture()
    n, temperature = 20_000, 1000.0
    normals = np.tile([0.6, 0.8, 0.0], (n, 1))
    particles = sample_thermal_emission(catalog, ("p1",) * n, np.zeros((n, 3)), normals,
                                      np.full(n, temperature), rng=np.random.default_rng(21))
    velocity = np.array([p.momentum_kg_m_s for p in particles]) / MP
    vn = np.sum(velocity * normals, axis=1)
    assert np.all(vn > 0)
    energy = .5 * MP * np.sum(velocity**2, axis=1)
    assert np.mean(energy) == pytest.approx(2 * 1.380649e-23 * temperature, rel=.02, abs=0)
    assert np.mean(vn) == pytest.approx(math.sqrt(math.pi * 1.380649e-23 * temperature / (2 * MP)), rel=.02)
    assert {p.state_id for p in particles} == {"p1"}


def test_thermal_launch_rejects_nonunit_normals():
    catalog, _, _, _ = fixture()
    with pytest.raises(ValueError, match="unit vectors"):
        sample_thermal_emission(catalog, ("p0",), [[0, 0, 0]], [[0, 2, 0]], [300],
                                rng=np.random.default_rng(2))


def test_energy_residual_is_invariant_to_large_conserved_chemical_reference():
    catalog, state, channel, _ = fixture()
    shifted = StateCatalog(catalog.species.values(), [
        InternalState(s.state_id, s.species, s.energy_ev + (1e15 if s.species != "e" else 0), s.quantum_numbers)
        for s in catalog.states.values()])
    options = dict(collision_volume_m3=1e-28, channels=(channel,),
                   differential_kernels={channel.name: kernel})
    first = advance(state, catalog, **options)
    second = advance(state, shifted, **options)
    assert first.reaction_counts == second.reaction_counts
    assert second.total_energy_residual_j == first.total_energy_residual_j
    assert second.reaction_energy_residual_j == first.reaction_energy_residual_j


def test_weight_and_volume_scaling_preserves_events_and_trajectories():
    catalog, state, channel, _ = fixture()
    scaled = GapState(tuple(Particle(p.state_id, p.position_m, p.momentum_kg_m_s, 8) for p in state.particles))
    options = dict(channels=(channel,), differential_kernels={channel.name: kernel}, electric_field_v_m=[0, 1, 0])
    first = advance(state, catalog, collision_volume_m3=1e-28, **options)
    second = advance(scaled, catalog, collision_volume_m3=8e-28, **options)
    assert first.reaction_counts == second.reaction_counts
    assert first.field_kinetic_gain_j * 8 == second.field_kinetic_gain_j
    for a, b in zip(first.state.particles, second.state.particles):
        assert a.state_id == b.state_id
        assert a.position_m == b.position_m
        assert a.momentum_kg_m_s == b.momentum_kg_m_s


def test_identical_reactants_count_each_unordered_pair_once():
    catalog, _, _, _ = fixture()
    speed = math.sqrt(5 * Q / (.75 * MP))
    velocities = [(speed, 0, 0), (-speed / 2, math.sqrt(3) * speed / 2, 0),
                  (-speed / 2, -math.sqrt(3) * speed / 2, 0)]
    state = GapState(tuple(Particle("p0", (0, 0, 0), tuple(MP * np.array(v))) for v in velocities))
    channel = ReactionChannel("pair_excitation", ("p0", "p0"), ("p1", "p1"),
                              (0, 2, 3, 100), (0, 0, 1e-18, 1e-18), "jacobi_energy_partition",
                              "Manufactured three-particle hazard check",
                              metadata={"kinematics_provenance": "Specified direction for verification"})
    dt, probability = 1e-10, 1 - math.exp(-.8)
    volume = 3 * 1e-18 * math.sqrt(3) * speed * dt / .8
    rng = np.random.default_rng(952)
    events = sum(advance(state, catalog, dt_s=dt, collision_volume_m3=volume,
                         channels=(channel,), differential_kernels={channel.name: kernel}, rng=rng)
                 .reaction_counts[channel.name] for _ in range(600))
    assert abs(events - 600 * probability) < 5 * math.sqrt(600 * probability * (1 - probability))


def test_relative_energy_below_threshold_does_not_react():
    catalog, state, channel, _ = fixture()
    # Projectile lab energy exceeds the level defect, but COM relative energy
    # is below it because a finite-mass stationary target must recoil.
    energy_lab_ev = 1 + .5 * ME / MP
    particle = Particle("e0", (0, 0, 0), (math.sqrt(2 * ME * energy_lab_ev * Q), 0, 0))
    result = advance(GapState((particle, state.particles[1])), catalog, collision_volume_m3=1e-40,
                     channels=(channel,), differential_kernels={channel.name: kernel})
    assert result.reaction_counts[channel.name] == 0
    assert result.state.particles[1].state_id == "p0"


def test_thermal_launch_transforms_with_wall_velocity_and_mass():
    catalog, _, _, _ = fixture()
    args = (np.zeros((4, 3)), np.tile([0, 1, 0], (4, 1)), np.full(4, 500))
    heavy = sample_thermal_emission(catalog, ("p0",) * 4, *args, rng=np.random.default_rng(88))
    light = sample_thermal_emission(catalog, ("e0",) * 4, *args, rng=np.random.default_rng(88))
    wall = np.tile([100, -200, 300], (4, 1))
    moving = sample_thermal_emission(catalog, ("p0",) * 4, *args, wall_velocity_m_s=wall,
                                    rng=np.random.default_rng(88))
    vh = np.array([p.momentum_kg_m_s for p in heavy]) / MP
    vl = np.array([p.momentum_kg_m_s for p in light]) / ME
    vm = np.array([p.momentum_kg_m_s for p in moving]) / MP
    np.testing.assert_allclose(vl / vh, math.sqrt(MP / ME), rtol=1e-13, atol=0)
    np.testing.assert_allclose(vm - vh, wall, rtol=1e-13, atol=1e-10)


def test_timestep_must_advance_representable_time():
    catalog, state, _, _ = fixture()
    with pytest.raises(ValueError, match="strictly increasing"):
        advance(GapState(state.particles, 1e12), catalog, dt_s=1e-10)


def test_complex_coordinates_and_temperatures_are_rejected():
    with pytest.raises(ValueError, match="real"):
        Particle("p0", np.array([0, 1j, 0]), [0, 0, 0])
    catalog, _, _, _ = fixture()
    with pytest.raises(ValueError, match="real"):
        sample_thermal_emission(catalog, ("p0",), [[0, 0, 0]], [[0, 1, 0]], np.array([300 + 1j]),
                                rng=np.random.default_rng(1))


def test_nonfinite_weighted_energy_ledger_is_rejected():
    catalog = StateCatalog([Species("X", 1e5, 0, {"X": 1})], [InternalState("x", "X", 0, {})])
    state = GapState((Particle("x", [0, 0, 0], [1e9, 0, 0], 1e308),))
    with pytest.raises(ValueError, match="finite numeric"):
        advance(state, catalog)


def test_competing_state_channels_follow_total_hazard_and_branch_ratio():
    catalog, state, first, speed = fixture()
    catalog = StateCatalog(catalog.species.values(), [*catalog.states.values(),
        InternalState("p_other", "P", 1, {"v": [0, 1, 0], "J": 0})])
    second = ReactionChannel("other_excitation", first.reactants, ("e0", "p_other"),
                             first.energy_ev, tuple(2 * x for x in first.sigma_m2),
                             first.kinematics_model, first.provenance, first.metadata)
    dt, rng, n = 1e-10, np.random.default_rng(615), 900
    volume = 3e-18 * speed * dt / .8
    counts = {first.name: 0, second.name: 0}
    for _ in range(n):
        result = advance(state, catalog, dt_s=dt, collision_volume_m3=volume,
                         channels=(first, second), differential_kernels={first.name: kernel, second.name: kernel}, rng=rng)
        for name in counts:
            counts[name] += result.reaction_counts[name]
    collision_probability = 1 - math.exp(-.8)
    for name, fraction in [(first.name, 1 / 3), (second.name, 2 / 3)]:
        probability = collision_probability * fraction
        assert abs(counts[name] - n * probability) < 5 * math.sqrt(n * probability * (1 - probability))
