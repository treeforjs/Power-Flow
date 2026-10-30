"""Rollback, identity, clocks, work ownership and stale/replay protections."""
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from mhdlab.thermal_expansion_mapping import affine_foil_geometry, conservative_remap_1d
from mhdlab.thermal_expansion_port import GeometryPort


MATERIAL, SOURCE = "a" * 64, "b" * 64


def test_parent_heat_and_material_bindings_match_frozen_contracts():
    root=Path(__file__).resolve().parents[1]
    contract=json.loads((root/"docs/interfaces/ss304-geometry-port-v1.json").read_text())
    for key,path in (("parent_requirements_sha256", "docs/interfaces/parent-coupling-port-requirements-20261029.json"),
                     ("heat_feedback_interface_sha256", "docs/interfaces/heat-feedback-port-20261029.json"),
                     ("mechanics_material_sha256", "data/materials/ss304-v1.json")):
        assert contract[key] == hashlib.sha256((root/path).read_bytes()).hexdigest()


def handoff(stretch=1.):
    return dict(
        geometries=[affine_foil_geometry(np.eye(3) * stretch, reference_density_kg_m3=7900.,
                                        center_reference_m=(0., y, 0.)) for y in (-.0026, .0026)],
        temperatures_k=[300., 300.], admissible_temperature_range_k=[293.15, 800.],
        transfers=[conservative_remap_1d([0, 1], [0, .5, 1], [1.])],
        mechanics_scope="uniform_affine_free_benchmark", support_description="free benchmark; mounting unknown",
    )


def port():
    return GeometryPort(mesh_revision="mesh-0", accepted_time_s=0., material_sha256=MATERIAL,
                        source_sha256=SOURCE, **handoff())


def prepare(owner, **overrides):
    arguments = dict(transaction_id="geometry-1", expected_revision=0, expected_mesh_revision="mesh-0",
                     next_mesh_revision="mesh-1", interval_s=[0., 1e-9], stage_time_s=1e-9,
                     material_sha256=MATERIAL, source_sha256=SOURCE, **handoff(1.001))
    arguments.update(overrides)
    return owner.prepare(**arguments)


def test_prepare_is_pure_for_accepted_state_and_commit_is_atomic():
    owner = port()
    initial = owner.state
    proposal = prepare(owner)
    assert owner.state is initial
    committed = owner.commit(proposal, expected_revision=0, expected_mesh_revision="mesh-0")
    assert owner.state is committed
    assert committed.revision == 1 and committed.mesh_revision == "mesh-1"
    assert committed.accepted_time_s == 1e-9
    assert committed.receipt()["geometries"][0]["mass_kg"] == initial.receipt()["geometries"][0]["mass_kg"]


def test_rollback_restores_exact_state_and_rejects_replay():
    owner = port()
    initial = owner.state
    proposal = prepare(owner)
    assert owner.rollback(proposal) is initial
    with pytest.raises(ValueError, match="replayed"):
        prepare(owner)
    with pytest.raises(ValueError, match="consumed"):
        owner.commit(proposal, expected_revision=0, expected_mesh_revision="mesh-0")


def test_stale_parallel_proposals_and_duplicate_commit_cannot_mutate():
    owner = port()
    first = prepare(owner)
    second = prepare(owner, transaction_id="geometry-2")
    accepted = owner.commit(first, expected_revision=0, expected_mesh_revision="mesh-0")
    for proposal in (first, second):
        with pytest.raises(ValueError):
            owner.commit(proposal, expected_revision=1, expected_mesh_revision="mesh-1")
        assert owner.state is accepted
    assert owner.rollback(second) is accepted


@pytest.mark.parametrize("override", [
    {"expected_revision": 1}, {"expected_mesh_revision": "wrong"}, {"source_sha256": "c" * 64},
    {"material_sha256": "c" * 64}, {"interval_s": [1e-9, 2e-9]}, {"interval_s": [0., 0.]},
    {"stage_time_s": 2e-9}, {"stage_time_s": .5e-9}, {"stage_time_s": 0.},
    {"stage_time_s": float("nan")}, {"next_mesh_revision": "mesh-0"},
    {"kinetic_energy_J": -1.}, {"elastic_trial_energy_J": float("inf")},
    {"species_event_ids": ["e1"]}, {"species_parent_ids": ["p1"]},
    {"temperatures_k": [900.]}, {"mechanics_scope": "nonuniform_transient"},
    {"application_mode": "physical_heat_field_pilot"},
])
def test_rejected_proposals_leave_state_untouched(override):
    owner = port()
    initial = owner.state
    with pytest.raises(ValueError):
        prepare(owner, **override)
    assert owner.state is initial
    # Rejection never consumes the transaction ID.
    prepare(owner)


def test_receipts_are_copies_and_forged_proposal_cannot_commit():
    owner = port()
    proposal = prepare(owner)
    view = proposal.receipt()
    view["material_sha256"] = "changed"
    assert proposal.receipt()["material_sha256"] == MATERIAL
    forged = replace(proposal, candidate=replace(proposal.candidate, accepted_time_s=2e-9))
    with pytest.raises(ValueError, match="altered"):
        owner.commit(forged, expected_revision=0, expected_mesh_revision="mesh-0")
    owner.commit(proposal, expected_revision=0, expected_mesh_revision="mesh-0")


def test_energy_diagnostics_never_claim_heat_debit_or_physical_closure():
    owner = port()
    receipt = prepare(owner, kinetic_energy_J=1., elastic_trial_energy_J=2.).receipt()
    assert receipt["kinetic_energy_J"] == 1. and receipt["elastic_trial_energy_J"] == 2.
    assert receipt["heat_debit_J"] is None and receipt["mechanical_work_J"] is None
    assert receipt["physical_coupled_ready"] is False
    assert receipt["species_event_ids"] == receipt["species_parent_ids"] == []
    assert receipt["species_mutation"] == "none"
    assert any("caloric EOS" in reason for reason in receipt["blockers"])


def test_material_mass_and_reference_coordinates_cannot_change_in_transaction():
    owner = port()
    bad = handoff()
    bad["geometries"] = [affine_foil_geometry(np.eye(3), reference_density_kg_m3=8000.,
                                            center_reference_m=(0., y, 0.)) for y in (-.0026, .0026)]
    with pytest.raises(ValueError, match="reference geometry/mass"):
        prepare(owner, **bad)
    with pytest.raises(AttributeError):
        owner.material_sha256 = "c" * 64


def test_old_mesh_revision_cannot_be_reused_after_commit():
    owner = port()
    proposal = prepare(owner)
    owner.commit(proposal, expected_revision=0, expected_mesh_revision="mesh-0")
    with pytest.raises(ValueError, match="mesh revision cannot be reused"):
        prepare(owner, transaction_id="geometry-2", expected_revision=1, expected_mesh_revision="mesh-1",
                next_mesh_revision="mesh-0", interval_s=[1e-9, 2e-9], stage_time_s=2e-9, **handoff(1.002))
