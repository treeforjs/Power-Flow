"""Local atomic geometry transactions for qualified affine benchmarks only.

The caller owns the global coupled commit/rollback coordinator. No physical heat
or field pilot may consume this port: a caloric EOS and qualified mechanical /
thermal work exchange are missing. Mechanical energies are diagnostics only.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import re

from .thermal_expansion_mapping import qualify_field_handoff


BLOCKERS = (
    "Nonuniform transient foil/support mechanics and field mesh transfer are not qualified.",
    "No qualified mechanical/thermal work exchange or caloric EOS; no physical heat debit is defined.",
    "Local transactions require an external atomic coupled coordinator and qualified physical material state.",
)


def _label(value, name):
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError(f"{name} must be a nonempty identifier")
    return value


def _digest(value, name):
    if not isinstance(value, str) or re.fullmatch("[0-9a-f]{64}", value) is None:
        raise ValueError(f"{name} must be a lowercase SHA256")
    return value


def _number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return float(value)


def _freeze(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


@dataclass(frozen=True)
class GeometryState:
    revision: int
    mesh_revision: str
    accepted_time_s: float
    payload: bytes

    def receipt(self):
        return json.loads(self.payload)


@dataclass(frozen=True)
class GeometryProposal:
    transaction_id: str
    base_revision: int
    base_mesh_revision: str
    candidate: GeometryState
    receipt_bytes: bytes

    def receipt(self):
        return json.loads(self.receipt_bytes)


class GeometryPort:
    """Single-process local transaction owner; external coordinator handles all ports."""

    def __init__(self, *, mesh_revision, accepted_time_s, material_sha256, source_sha256, **handoff):
        self._material_sha256 = _digest(material_sha256, "material_sha256")
        self._source_sha256 = _digest(source_sha256, "source_sha256")
        clock = _number(accepted_time_s, "accepted_time_s")
        if clock < 0:
            raise ValueError("accepted clock must be nonnegative")
        qualified = qualify_field_handoff(**handoff)
        self._state = GeometryState(0, _label(mesh_revision, "mesh_revision"), clock, _freeze(qualified))
        self._pending, self._consumed = {}, set()
        self._mesh_history = {mesh_revision}

    @property
    def material_sha256(self):
        return self._material_sha256

    @property
    def source_sha256(self):
        return self._source_sha256

    @property
    def state(self):
        return self._state

    def readiness(self):
        return {"schema": "solid-geometry-transaction-v1", "physical_coupled_ready": False,
                "benchmark_transactions_ready": True, "blockers": list(BLOCKERS),
                "temperature_stage": "trial_interval_end; source emission uses previous accepted start temperature",
                "heat_debit_J": None, "mechanical_work_J": None,
                "replay_protection_scope": "current port instance; checkpoint persistence not implemented",
                "rollback_semantics": "discard prepared proposal; accepted commits are final",
                "energy_ownership": "KE and elastic trial energy are diagnostics, not heat debit or a closure ledger"}

    def _expected(self, revision, mesh_revision):
        if type(revision) is not int or revision != self._state.revision or mesh_revision != self._state.mesh_revision:
            raise ValueError("stale state/mesh revision")

    def prepare(self, *, transaction_id, expected_revision, expected_mesh_revision, next_mesh_revision,
                interval_s, stage_time_s, material_sha256, source_sha256, kinetic_energy_J=0.,
                elastic_trial_energy_J=0., species_event_ids=(), species_parent_ids=(),
                application_mode="affine_benchmark", **handoff):
        if application_mode != "affine_benchmark":
            raise ValueError("physical coupled geometry blocked: " + " ".join(BLOCKERS))
        self._expected(expected_revision, expected_mesh_revision)
        transaction_id = _label(transaction_id, "transaction_id")
        if transaction_id in self._pending or transaction_id in self._consumed:
            raise ValueError("duplicate or replayed geometry transaction")
        if material_sha256 != self.material_sha256 or source_sha256 != self.source_sha256:
            raise ValueError("material/source SHA256 binding mismatch")
        if species_event_ids or species_parent_ids:
            raise ValueError("geometry port performs no species mutation; species event/parent IDs must be empty")
        if len(interval_s) != 2:
            raise ValueError("accepted interval must contain [t, t+dt]")
        start, end = (_number(t, "interval time") for t in interval_s)
        stage = _number(stage_time_s, "stage_time_s")
        if start != self._state.accepted_time_s or end <= start or stage != end:
            raise ValueError("interval/stage clock does not match accepted state")
        kinetic, elastic = (_number(kinetic_energy_J, "kinetic_energy_J"),
                            _number(elastic_trial_energy_J, "elastic_trial_energy_J"))
        if min(kinetic, elastic) < 0:
            raise ValueError("mechanical diagnostic energies must be nonnegative")
        qualified = qualify_field_handoff(**handoff)
        next_mesh_revision = _label(next_mesh_revision, "next_mesh_revision")
        old_geometry = self._state.receipt()["geometries"]
        if len(qualified["geometries"]) != len(old_geometry):
            raise ValueError("material topology changes require a separate conservative remesh port")
        for old, new in zip(old_geometry, qualified["geometries"]):
            for key in ("reference_dimensions_m", "center_reference_m", "mass_kg"):
                if old[key] != new[key]:
                    raise ValueError("material reference geometry/mass must remain invariant")
        if qualified["geometries"] != old_geometry and next_mesh_revision == expected_mesh_revision:
            raise ValueError("changed geometry requires a new mesh revision")
        if next_mesh_revision != expected_mesh_revision and next_mesh_revision in self._mesh_history:
            raise ValueError("previously accepted mesh revision cannot be reused")
        candidate = GeometryState(expected_revision + 1, next_mesh_revision, end, _freeze(qualified))
        receipt = dict(self.readiness(), transaction_id=transaction_id, state_revision=expected_revision,
                       mesh_revision=expected_mesh_revision, next_mesh_revision=next_mesh_revision,
                       interval_s=[start, end], stage_time_s=stage, material_sha256=self.material_sha256,
                       source_sha256=self.source_sha256, candidate_geometry_sha256=hashlib.sha256(candidate.payload).hexdigest(),
                       species_event_ids=[], species_parent_ids=[], species_mutation="none",
                       kinetic_energy_J=kinetic, elastic_trial_energy_J=elastic, geometry=qualified)
        proposal = GeometryProposal(transaction_id, expected_revision, expected_mesh_revision, candidate, _freeze(receipt))
        self._pending[transaction_id] = proposal
        return proposal

    def commit(self, proposal, *, expected_revision, expected_mesh_revision):
        self._expected(expected_revision, expected_mesh_revision)
        if not isinstance(proposal, GeometryProposal) or self._pending.get(proposal.transaction_id) is not proposal:
            raise ValueError("unknown, consumed, or altered proposal")
        if proposal.base_revision != expected_revision or proposal.base_mesh_revision != expected_mesh_revision:
            raise ValueError("stale proposal")
        self._state = proposal.candidate
        self._mesh_history.add(self._state.mesh_revision)
        del self._pending[proposal.transaction_id]
        self._consumed.add(proposal.transaction_id)
        return self._state

    def rollback(self, proposal):
        if not isinstance(proposal, GeometryProposal) or self._pending.get(proposal.transaction_id) is not proposal:
            raise ValueError("unknown, consumed, or altered proposal")
        del self._pending[proposal.transaction_id]
        self._consumed.add(proposal.transaction_id)
        return self._state
