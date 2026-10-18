"""Analytic reference contracts, deliberately not a physical water mechanism."""

from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from mhdlab.state_kinetics import (
    CrossSectionDomainError,
    InternalState,
    RadiativeTransition,
    ReactionChannel,
    Species,
    StateCatalog,
)


def catalog():
    # Artificial masses/energies make the conservation checks independently exact.
    electron_mass, atom_mass = 1e-30, 10e-27
    species = [Species("e", electron_mass, -1, {}),
               Species("A", atom_mass, 0, {"A": 1}),
               Species("A+", atom_mass - electron_mass, 1, {"A": 1}),
               Species("A2", 2 * atom_mass, 0, {"A": 2})]
    states = [InternalState("electron", "e", 0.0, {}),
              InternalState("A:g", "A", -1.0, {"n": 1}),
              InternalState("A:x", "A", 1.0, {"n": 2}),
              InternalState("A+:g", "A+", 3.0, {"term": "ground"}),
              InternalState("A2:000", "A2", -5.0, {"mode_tuple": (0, 0, 0), "J": 0})]
    return StateCatalog(species, states)


def excitation(**overrides):
    values = dict(name="analytic_excitation", reactants=("electron", "A:g"),
                  products=("electron", "A:x"), energy_ev=(2.0, 4.0, 6.0),
                  sigma_m2=(0.0, 2e-20, 4e-20), kinematics_model="isotropic_com_two_body",
                  provenance="Analytically specified reference fixture; no molecular dataset",
                  metadata={"kinematics_provenance": "Explicit isotropic COM reference model"})
    values.update(overrides)
    return ReactionChannel(**values)


def test_explicit_polyatomic_assignments_and_recursive_immutability():
    assignments = {"modes": [1, 0, 2], "J": 3, "Ka": 1, "Kc": 2,
                   "symmetry": {"rovib": "A1", "parity": "+"}}
    nuclei = {"H": 2, "O": 1}
    water = Species("H2O", 3e-26, 0, nuclei)
    state = InternalState("source:water:17", "H2O", -10.0, assignments)
    sample = StateCatalog([water], [state])
    assignments["modes"][0] = 99
    nuclei["H"] = 99
    assert sample.state("source:water:17").quantum_numbers["modes"] == (1, 0, 2)
    assert sample.species["H2O"].nuclei["H"] == 2
    with pytest.raises(TypeError):
        sample.states["new"] = state
    with pytest.raises(TypeError):
        state.quantum_numbers["symmetry"]["parity"] = "-"
    with pytest.raises(FrozenInstanceError):
        state.energy_ev = 1


def test_lookup_never_falls_back_to_species_or_unlisted_quantum_state():
    sample = catalog()
    with pytest.raises(KeyError, match="fallback is forbidden"):
        sample.state("A")
    with pytest.raises(KeyError):
        excitation(products=("electron", "A:missing")).validate(sample)
    with pytest.raises(ValueError, match="required species"):
        sample.require_species(["H2O"])


def test_fingerprint_tracks_energy_and_quantum_labels_independent_of_order():
    sample = catalog()
    reordered = StateCatalog(reversed(list(sample.species.values())), reversed(list(sample.states.values())))
    assert reordered.fingerprint == sample.fingerprint
    states = dict(sample.states)
    states["A:x"] = replace(states["A:x"], energy_ev=2)
    changed = StateCatalog(sample.species, states)
    assert changed.fingerprint != sample.fingerprint
    with pytest.raises(ValueError, match="fingerprint"):
        excitation(catalog_fingerprint=sample.fingerprint).validate(changed)


def test_threshold_uses_formation_and_excitation_energy_and_linear_sigma():
    sample, channel = catalog(), excitation()
    channel.validate(sample)
    assert channel.energy_defect_ev(sample) == 2
    assert channel.threshold_ev(sample) == 2
    np.testing.assert_allclose(channel.sigma_at([0, 1, 2, 3, 4, 6], sample),
                               [0, 0, 0, 1e-20, 2e-20, 4e-20], atol=0)


def test_supported_threshold_and_unspecified_domain_are_not_conflated():
    sample = catalog()
    incomplete = excitation(energy_ev=(3.0, 6.0), sigma_m2=(1e-20, 4e-20))
    assert incomplete.sigma_at(1.99, sample) == 0
    for query in (2, 2.5, 6.1, [1.0, 7.0]):
        with pytest.raises(CrossSectionDomainError):
            incomplete.sigma_at(query, sample)
    with pytest.raises(ValueError, match="below.*threshold"):
        excitation(energy_ev=(1.0, 3.0), sigma_m2=(1e-20, 2e-20)).validate(sample)


def test_exothermic_channel_has_zero_threshold_and_no_zero_below_grid():
    sample = catalog()
    reverse = excitation(reactants=("electron", "A:x"), products=("electron", "A:g"),
                         energy_ev=(0.0, 5.0), sigma_m2=(3e-20, 1e-20))
    assert reverse.energy_defect_ev(sample) == -2
    assert reverse.threshold_ev(sample) == 0
    assert reverse.sigma_at(0, sample) == 3e-20
    truncated = replace(reverse, energy_ev=(1.0, 5.0))
    with pytest.raises(CrossSectionDomainError):
        truncated.sigma_at(0.5, sample)


def test_ionization_and_dissociation_stoichiometry_use_all_state_energies():
    sample = catalog()
    ionization = excitation(products=("electron", "electron", "A+:g"),
                            energy_ev=(4, 6), sigma_m2=(0, 1e-20),
                            kinematics_model="electron_scattered_three_body")
    ionization.validate(sample)
    assert ionization.threshold_ev(sample) == 4
    dissociation = excitation(reactants=("electron", "A2:000"),
                             products=("electron", "A:g", "A:x"),
                             energy_ev=(5, 6), sigma_m2=(0, 1e-20),
                             kinematics_model="electron_scattered_three_body")
    dissociation.validate(sample)
    assert dissociation.threshold_ev(sample) == 5


@pytest.mark.parametrize("products,match", [
    (("electron", "A:g", "A:g"), "nuclei"),
    (("electron", "A+:g"), "charge"),
])
def test_unconserved_reaction_rejected(products, match):
    with pytest.raises(ValueError, match=match):
        excitation(products=products).validate(catalog())


def test_inconsistent_species_mass_convention_rejected():
    sample = catalog()
    species = dict(sample.species)
    species["A+"] = replace(species["A+"], mass_kg=species["A"].mass_kg)
    inconsistent = StateCatalog(species, sample.states)
    channel = excitation(products=("electron", "electron", "A+:g"),
                         energy_ev=(4, 6), sigma_m2=(0, 1e-20))
    with pytest.raises(ValueError, match="mass"):
        channel.validate(inconsistent)


@pytest.mark.parametrize("overrides,match", [
    ({"provisional": True}, "provisional"),
    ({"provenance": "Placeholder pending measurements"}, "placeholder"),
    ({"kinematics_model": "provisional_isotropic"}, "provisional"),
    ({"metadata": {}}, "kinematics_provenance"),
    ({"metadata": {"exclusive_reaction": False}}, "nonexclusive"),
    ({"metadata": {"data_kind": "electron_impact_partial_ion_yield"}}, "yield/total"),
    ({"metadata": {"data_kind": "electron_impact_total_ionization"}}, "yield/total"),
])
def test_production_requires_ready_exclusive_channel_and_kinematics_provenance(overrides, match):
    channel = excitation(**overrides)
    with pytest.raises(ValueError, match=match):
        channel.validate(catalog(), production=True)
    channel.validate(catalog(), production=False)


def test_disabled_provisional_channels_cannot_be_evaluated():
    channel = excitation(enabled=False, provisional=True)
    channel.validate(catalog(), production=True)
    with pytest.raises(ValueError, match="disabled"):
        channel.sigma_at(3, catalog())


@pytest.mark.parametrize("overrides", [
    {"energy_ev": [1, 1, 2]}, {"energy_ev": [3, 2, 1]},
    {"energy_ev": [1, 2]}, {"sigma_m2": [0, -1, 1]},
    {"sigma_m2": [0, np.nan, 1]}, {"energy_ev": [0, np.inf, 2]},
    {"energy_axis": "incident_lab_energy_ev"}, {"reactants": "A:g"},
    {"products": ()}, {"provenance": ""}, {"metadata": {"uncertainty": np.inf}},
])
def test_malformed_channel_inputs_rejected(overrides):
    with pytest.raises(ValueError):
        excitation(**overrides)


@pytest.mark.parametrize("query", [-1, np.nan, np.inf, [], [[3]], "bad", "3", True, ["3"]])
def test_invalid_query_rejected(query):
    with pytest.raises(ValueError):
        excitation().sigma_at(query, catalog())


def test_catalog_rejects_duplicate_ids_unknown_species_and_unassigned_molecule():
    sample = catalog()
    with pytest.raises(ValueError, match="duplicate"):
        StateCatalog(sample.species, list(sample.states.values()) * 2)
    with pytest.raises(ValueError, match="unknown species"):
        StateCatalog(sample.species, [InternalState("x", "absent", 0, {})])
    with pytest.raises(ValueError, match="quantum assignments"):
        StateCatalog([sample.species["A2"]], [InternalState("A2:g", "A2", 0, {})])


def test_radiative_transition_uses_explicit_upper_population_state():
    sample = catalog()
    line = RadiativeTransition("A:x", "A:g", 1e8, "Analytic reference transition")
    assert line.photon_energy_ev(sample) == 2
    with pytest.raises(ValueError, match="positive"):
        replace(line, upper_state="A:g", lower_state="A:x").validate(sample)
    with pytest.raises(ValueError, match="same species"):
        replace(line, lower_state="A+:g").validate(sample)
    with pytest.raises(ValueError):
        replace(line, a_per_s=0)
