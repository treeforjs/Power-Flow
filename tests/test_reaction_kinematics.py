"""Conservation/analytic checks using synthetic inputs, not a water mechanism."""

from dataclasses import replace

import numpy as np
import pytest

from mhdlab.reaction_kinematics import (
    InaccessibleReactionError,
    NonrelativisticLimitError,
    SPEED_OF_LIGHT_M_S,
    react,
)
from mhdlab.state_kinetics import (
    ELEMENTARY_CHARGE_C as EV_J,
    InternalState,
    ReactionChannel,
    Species,
    StateCatalog,
)


def reference_catalog():
    # Artificial, exactly additive masses and common-reference energies make
    # the expected kinematics independent of a physical molecular dataset.
    m_e, m_a, m_b = 1e-30, 2e-27, 6e-27
    species = [Species("e", m_e, -1, {}),
               Species("A", m_a, 0, {"A": 1}),
               Species("B", m_b, 0, {"B": 1}),
               Species("AB", m_a + m_b, 0, {"A": 1, "B": 1}),
               Species("AB+", m_a + m_b - m_e, 1, {"A": 1, "B": 1}),
               Species("A2B", 2 * m_a + m_b, 0, {"A": 2, "B": 1})]
    states = [InternalState("e", "e", 0, {}),
              InternalState("A:g", "A", 0, {"term": "ground"}),
              InternalState("A:x", "A", 2, {"term": "excited"}),
              InternalState("B:g", "B", 0, {"term": "ground"}),
              InternalState("AB:g", "AB", -5, {"v": 0, "J": 0}),
              InternalState("AB+:g", "AB+", 3, {"v": 0, "J": 0}),
              InternalState("A2B:g", "A2B", -8, {"modes": (0, 0, 0), "J": 0})]
    return StateCatalog(species, states)


def channel(reactants=("A:g", "B:g"), products=("A:g", "B:g"), **kwargs):
    values = dict(name="analytic_conservation_reference", reactants=reactants,
                  products=products, energy_ev=(0, 10_000), sigma_m2=(0, 0),
                  kinematics_model="jacobi_energy_partition",
                  provenance="Analytic kinematics fixture, no physical cross sections",
                  metadata={"kinematics_provenance": "Supplied deterministic Jacobi modes for analytic checks"})
    values.update(kwargs)
    return ReactionChannel(**values)


def incoming_com_momenta(catalog, states, relative_energy_ev, boost=(0, 0, 0)):
    masses = np.array([catalog.mass_kg(state) for state in states])
    reduced = masses[0] * masses[1] / masses.sum()
    magnitude = np.sqrt(2 * reduced * relative_energy_ev * EV_J)
    return np.array([[magnitude, 0, 0], [-magnitude, 0, 0]]) + masses[:, None] * np.array(boost)


def assert_conservative(outcome, tolerance=3e-13):
    scale = max(outcome.kinetic_energy_before_j, outcome.kinetic_energy_after_j,
                abs(outcome.energy_defect_ev) * EV_J)
    assert abs(outcome.energy_residual_j) <= tolerance * scale
    momentum_scale = max(np.linalg.norm(outcome.outgoing_momenta_kg_m_s, axis=1).max(),
                         np.linalg.norm(outcome.momentum_before_kg_m_s))
    assert np.linalg.norm(outcome.momentum_residual_kg_m_s) <= tolerance * momentum_scale
    assert outcome.outgoing_com_kinetic_energy_j == pytest.approx(outcome.available_energy_j, rel=tolerance, abs=0)


def test_two_body_closed_form_unequal_mass_recoil_and_readonly_ledger():
    catalog, reaction = reference_catalog(), channel()
    incoming = incoming_com_momenta(catalog, reaction.reactants, 4)
    outcome = react(reaction, catalog, incoming, jacobi_energy_fractions=[1],
                    jacobi_directions=[[0, 1, 0]])
    masses = np.array([catalog.mass_kg(state) for state in reaction.products])
    kinetic = np.sum(outcome.outgoing_momenta_kg_m_s ** 2, axis=1) / (2 * masses)
    np.testing.assert_allclose(kinetic / EV_J, [3, 1], rtol=2e-15)
    np.testing.assert_allclose(outcome.outgoing_momenta_kg_m_s[0],
                               -outcome.outgoing_momenta_kg_m_s[1], atol=1e-38)
    assert outcome.product_state_ids == reaction.products
    assert outcome.available_energy_j / EV_J == pytest.approx(4)
    assert_conservative(outcome)
    with pytest.raises(ValueError):
        outcome.outgoing_momenta_kg_m_s[0, 0] = 0
    with pytest.raises(TypeError):
        outcome.metadata["physical_default"] = True


@pytest.mark.parametrize("initial,final,defect,available", [
    ("A:g", "A:x", 2, 3), ("A:x", "A:g", -2, 7),
])
def test_state_changes_transfer_internal_energy_and_preserve_moving_com(initial, final, defect, available):
    catalog = reference_catalog()
    reaction = channel(("e", initial), ("e", final))
    incoming = incoming_com_momenta(catalog, reaction.reactants, 5, boost=(320, -120, 81))
    outcome = react(reaction, catalog, incoming, jacobi_energy_fractions=[1],
                    jacobi_directions=[[0, 0, 1]])
    assert outcome.energy_defect_ev == defect
    assert outcome.available_energy_j / EV_J == pytest.approx(available)
    assert outcome.internal_energy_after_j - outcome.internal_energy_before_j == pytest.approx(defect * EV_J, rel=1e-12, abs=0)
    np.testing.assert_allclose(outcome.center_of_mass_velocity_m_s, [320, -120, 81], rtol=3e-15)
    assert_conservative(outcome)


def test_three_body_fragmentation_has_correct_electron_and_heavy_pair_recoil():
    catalog = reference_catalog()
    reaction = channel(("e", "AB:g"), ("e", "A:g", "B:g"))
    incoming = incoming_com_momenta(catalog, reaction.reactants, 15)
    outcome = react(reaction, catalog, incoming, jacobi_energy_fractions=[0.6, 0.4],
                    jacobi_directions=[[1, 0, 0], [0, 1, 0]])
    assert outcome.available_energy_j / EV_J == pytest.approx(10)
    masses = np.array([catalog.mass_kg(state) for state in reaction.products])
    velocities = outcome.outgoing_momenta_kg_m_s / masses[:, None]
    electron_energy = np.dot(outcome.outgoing_momenta_kg_m_s[0],
                             outcome.outgoing_momenta_kg_m_s[0]) / (2 * masses[0])
    assert electron_energy / EV_J == pytest.approx(6 * masses[1:].sum() / masses.sum())
    pair_relative_energy = 0.5 * masses[1] * masses[2] / masses[1:].sum() * np.sum((velocities[1] - velocities[2]) ** 2)
    assert pair_relative_energy / EV_J == pytest.approx(4)
    assert_conservative(outcome)


def test_ionization_secondary_electron_momenta_and_charge_are_conserved():
    catalog = reference_catalog()
    reaction = channel(("e", "AB:g"), ("e", "e", "AB+:g"))
    outcome = react(reaction, catalog, incoming_com_momenta(catalog, reaction.reactants, 20, (20, 30, 40)),
                    jacobi_energy_fractions=[0.8, 0.2],
                    jacobi_directions=[[0, 1, 0], [1, 0, 0]])
    assert outcome.product_state_ids == ("e", "e", "AB+:g")
    assert outcome.energy_defect_ev == 8
    assert outcome.available_energy_j / EV_J == pytest.approx(12)
    assert sum(catalog.charge_number(state) for state in reaction.reactants) == sum(
        catalog.charge_number(state) for state in outcome.product_state_ids)
    assert_conservative(outcome)


def test_general_four_product_jacobi_modes_conserve_energy_and_momentum():
    catalog = reference_catalog()
    reaction = channel(("e", "A2B:g"), ("e", "A:g", "A:g", "B:g"))
    outcome = react(reaction, catalog, incoming_com_momenta(catalog, reaction.reactants, 20, (33, -16, 80)),
                    jacobi_energy_fractions=[0.5, 0.2, 0.3],
                    jacobi_directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert outcome.available_energy_j / EV_J == pytest.approx(12)
    assert_conservative(outcome)


def test_product_order_is_part_of_the_kernel_and_not_silently_sorted():
    catalog = reference_catalog()
    first = channel(("e", "AB:g"), ("e", "A:g", "B:g"))
    permuted = replace(first, products=("A:g", "e", "B:g"))
    incoming = incoming_com_momenta(catalog, first.reactants, 15)
    kwargs = dict(jacobi_energy_fractions=[0.6, 0.4], jacobi_directions=[[1, 0, 0], [0, 1, 0]])
    original, changed = react(first, catalog, incoming, **kwargs), react(permuted, catalog, incoming, **kwargs)
    assert changed.product_state_ids == permuted.products
    assert not np.allclose(changed.outgoing_momenta_kg_m_s[1], original.outgoing_momenta_kg_m_s[0], atol=1e-30)
    assert_conservative(original)
    assert_conservative(changed)


def test_rotation_and_galilean_boost_commute_with_constructor():
    catalog = reference_catalog()
    reaction = channel(("e", "AB:g"), ("e", "A:g", "B:g"))
    incoming = incoming_com_momenta(catalog, reaction.reactants, 15)
    directions = np.array([[1, 0, 0], [0, 1, 0]])
    base = react(reaction, catalog, incoming, jacobi_energy_fractions=[0.6, 0.4], jacobi_directions=directions)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    rotated = react(reaction, catalog, incoming @ rotation.T, jacobi_energy_fractions=[0.6, 0.4],
                    jacobi_directions=directions @ rotation.T)
    np.testing.assert_allclose(rotated.outgoing_momenta_kg_m_s, base.outgoing_momenta_kg_m_s @ rotation.T,
                               rtol=2e-15, atol=1e-38)
    boost = np.array([400, -200, 80])
    masses_in = np.array([catalog.mass_kg(state) for state in reaction.reactants])
    masses_out = np.array([catalog.mass_kg(state) for state in reaction.products])
    boosted = react(reaction, catalog, incoming + masses_in[:, None] * boost,
                    jacobi_energy_fractions=[0.6, 0.4], jacobi_directions=directions)
    np.testing.assert_allclose(boosted.outgoing_momenta_kg_m_s,
                               base.outgoing_momenta_kg_m_s + masses_out[:, None] * boost,
                               rtol=3e-15, atol=1e-38)
    assert_conservative(boosted)


def test_endothermic_threshold_is_not_replaced_by_an_energy_floor():
    catalog = reference_catalog()
    reaction = channel(("e", "AB:g"), ("e", "A:g", "B:g"))
    kwargs = dict(jacobi_energy_fractions=[0.6, 0.4], jacobi_directions=[[1, 0, 0], [0, 1, 0]])
    with pytest.raises(InaccessibleReactionError):
        react(reaction, catalog, incoming_com_momenta(catalog, reaction.reactants, 5 - 1e-8), **kwargs)
    threshold = react(reaction, catalog, incoming_com_momenta(catalog, reaction.reactants, 5), **kwargs)
    assert threshold.available_energy_j <= 64 * np.finfo(float).eps * 5 * EV_J
    just_above = react(reaction, catalog, incoming_com_momenta(catalog, reaction.reactants, 5 + 1e-8), **kwargs)
    assert just_above.available_energy_j / EV_J == pytest.approx(1e-8, rel=1e-6)


def test_incoming_and_exothermic_outgoing_speed_caps_are_enforced():
    catalog = reference_catalog()
    reaction = channel()
    incoming = np.array([[catalog.mass_kg("A:g") * 0.06 * SPEED_OF_LIGHT_M_S, 0, 0], [0, 0, 0]])
    with pytest.raises(NonrelativisticLimitError, match="incoming"):
        react(reaction, catalog, incoming, jacobi_energy_fractions=[1], jacobi_directions=[[1, 0, 0]])
    reverse = channel(("e", "A:x"), ("e", "A:g"))
    with pytest.raises(NonrelativisticLimitError, match="outgoing"):
        react(reverse, catalog, np.zeros((2, 3)), jacobi_energy_fractions=[1],
              jacobi_directions=[[1, 0, 0]], max_nonrel_beta=0.001)


@pytest.mark.parametrize("overrides,match", [
    ({"products": ("A:g", "A:g")}, "nuclei"),
    ({"products": ("e", "AB:g")}, "charge"),
    ({"kinematics_model": "isotropic_com_two_body"}, "unsupported kinematics_model"),
    ({"enabled": False}, "disabled"),
    ({"metadata": {}}, "kinematics_provenance"),
    ({"provisional": True}, "provisional"),
])
def test_channel_validation_is_not_bypassed(overrides, match):
    catalog = reference_catalog()
    reaction = channel(**overrides)
    with pytest.raises(ValueError, match=match):
        react(reaction, catalog, np.zeros((2, 3)), jacobi_energy_fractions=[1], jacobi_directions=[[1, 0, 0]])


@pytest.mark.parametrize("fractions,directions,match", [
    ([0.5], [[1, 0, 0]], "sum to one"),
    ([-1], [[1, 0, 0]], "nonnegative"),
    ([float("nan")], [[1, 0, 0]], "finite real array"),
    ([1, 0], [[1, 0, 0]], "shape"),
    ([1], [[0, 0, 0]], "unit vectors"),
    ([1], [[2, 0, 0]], "unit vectors"),
    ([1], [[1, float("inf"), 0]], "finite real array"),
    ([1], [[1, 0]], "shape"),
])
def test_explicit_partition_validation(fractions, directions, match):
    with pytest.raises(ValueError, match=match):
        react(channel(), reference_catalog(), np.zeros((2, 3)),
              jacobi_energy_fractions=fractions, jacobi_directions=directions)


@pytest.mark.parametrize("momenta", [np.zeros((3, 3)), [[0, 0, np.nan], [0, 0, 0]],
                                    [[1j, 0, 0], [0, 0, 0]]])
def test_invalid_incoming_momenta_are_rejected(momenta):
    with pytest.raises(ValueError, match="incoming_momenta_kg_m_s"):
        react(channel(), reference_catalog(), momenta,
              jacobi_energy_fractions=[1], jacobi_directions=[[1, 0, 0]])


def test_inconsistent_nonrelativistic_mass_convention_is_rejected():
    catalog = reference_catalog()
    species = dict(catalog.species)
    species["AB+"] = replace(species["AB+"], mass_kg=catalog.mass_kg("AB:g"))
    invalid_catalog = StateCatalog(species, catalog.states)
    reaction = channel(("e", "AB:g"), ("e", "e", "AB+:g"))
    with pytest.raises(ValueError, match="mass"):
        react(reaction, invalid_catalog, np.zeros((2, 3)),
              jacobi_energy_fractions=[0.5, 0.5],
              jacobi_directions=[[1, 0, 0], [0, 1, 0]])


@pytest.mark.parametrize("cap", [0, -0.1, 1, np.nan, True, "not-a-speed"])
def test_invalid_nonrelativistic_caps_are_rejected(cap):
    with pytest.raises(ValueError, match="max_nonrel_beta"):
        react(channel(), reference_catalog(), np.zeros((2, 3)),
              jacobi_energy_fractions=[1], jacobi_directions=[[1, 0, 0]], max_nonrel_beta=cap)


def test_only_two_reactants_and_multiple_massive_products_are_supported():
    catalog = reference_catalog()
    for reaction in [channel(("AB:g",), ("A:g", "B:g")),
                     channel(("A:g", "B:g"), ("AB:g",))]:
        with pytest.raises(ValueError, match="exactly two reactants"):
            react(reaction, catalog, np.zeros((2, 3)),
                  jacobi_energy_fractions=[1], jacobi_directions=[[1, 0, 0]])
