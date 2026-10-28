"""Validate an external source receipt and resolved population observation bundle.

This module runs no source, transport or CR model. Hashes bind supplied evidence;
they do not independently establish the origin or physical validity of that evidence.
"""

from __future__ import annotations

import json
from hashlib import sha256
from numbers import Real
from typing import Mapping

import numpy as np


ARTIFACTS = frozenset({"source", "receipt", "states", "rates", "populations", "mesh",
                       "collection", "electrons", "opacity", "runtime_output"})
MOMENTS = ("released_number", "local_time_first_moment_s", "local_time_second_moment_s2",
           "temperature_first_moment_k", "temperature_second_moment_k2",
           "local_time_temperature_moment_k_s", "binding_energy_j", "internal_energy_j",
           "translational_energy_j")
CHANNELS = frozenset({"excitation", "radiative_branching", "cascades", "quenching",
                      "predissociation", "ionization", "recombination", "fragmentation"})


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be nonempty text")
    return value


def _number(value, name, *, nonnegative=False, positive=False):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be finite numeric data")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite numeric data") from exc
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite numeric data")
    if (nonnegative and value < 0) or (positive and value <= 0):
        raise ValueError(f"{name} has an invalid sign")
    return value


def _array(value, name, shape=None, *, nonnegative=False):
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.size == 0:
        raise ValueError(f"{name} must be nonempty numeric data")
    array = np.asarray(raw, dtype=float)
    if (shape is not None and array.shape != shape) or not np.isfinite(array).all():
        raise ValueError(f"{name} has invalid shape or nonfinite values")
    if nonnegative and np.any(array < 0):
        raise ValueError(f"{name} must be nonnegative")
    return array


def _unique(values, name):
    if not isinstance(values, list) or not values:
        raise ValueError(f"{name} must be a nonempty list")
    for value in values:
        _text(value, name)
    if len(set(values)) != len(values):
        raise ValueError(f"{name} contains duplicate identifiers")
    return values


def _digest(value, name):
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"{name} must be a lowercase SHA256")
    return value


def strict_json(raw: bytes):
    """Read finite JSON without duplicate keys; source bytes remain the hash identity."""
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def constant(value):
        raise ValueError(f"nonfinite JSON value: {value}")

    def floating(value):
        result = float(value)
        if not np.isfinite(result):
            raise ValueError("nonfinite JSON number")
        return result

    try:
        result = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs,
                            parse_constant=constant, parse_float=floating)
    except (UnicodeError, json.JSONDecodeError, RecursionError) as exc:
        raise ValueError("invalid JSON artifact") from exc
    return result


def collected_path_weights(cell_bounds_m, ray_start_m, ray_end_m, collection_weights):
    """Intersect finite rays with disjoint axis-aligned simulation cells.

    Weights are dimensionless radiance averages and must sum to one. Cells own
    [lower, upper) on parallel faces, except the global outer maximum is closed.
    This gives a ray on an internal face exactly one owner. No etendue, gate,
    material thickness or population factor is included in the returned metres.
    """
    bounds = _array(cell_bounds_m, "cell_bounds_m")
    if bounds.ndim != 3 or bounds.shape[1:] != (2, 3) or np.any(bounds[:, 1] <= bounds[:, 0]):
        raise ValueError("cell_bounds_m must have shape (cell,2,3) with positive extents")
    starts = _array(ray_start_m, "ray_start_m")
    if starts.ndim != 2 or starts.shape[1] != 3:
        raise ValueError("ray_start_m must have shape (ray,3)")
    stops = _array(ray_end_m, "ray_end_m", starts.shape)
    weights = _array(collection_weights, "collection_weights", (len(starts),), nonnegative=True)
    if not np.isclose(np.sum(weights), 1.0, rtol=0, atol=1e-12):
        raise ValueError("collection weights must already sum to one; no normalization is applied")
    for index, box in enumerate(bounds):
        if np.any(np.all(np.minimum(box[1], bounds[index + 1:, 1])
                         > np.maximum(box[0], bounds[index + 1:, 0]), axis=1)):
            raise ValueError("simulation cells overlap in their interiors")
    directions = stops - starts
    lengths = np.linalg.norm(directions, axis=1)
    if np.any(lengths <= 0) or not np.isfinite(lengths).all():
        raise ValueError("rays must have finite positive lengths")
    outer_max = np.max(bounds[:, 1], axis=0)
    intersections = np.zeros((len(starts), len(bounds)))
    for ray, (start, direction, length) in enumerate(zip(starts, directions, lengths)):
        for cell, box in enumerate(bounds):
            low, high = 0.0, 1.0
            for axis in range(3):
                if direction[axis] == 0:
                    if (start[axis] < box[0, axis] or start[axis] > box[1, axis]
                            or (start[axis] == box[1, axis] and box[1, axis] != outer_max[axis])):
                        high = low
                        break
                else:
                    a = (box[0, axis] - start[axis]) / direction[axis]
                    b = (box[1, axis] - start[axis]) / direction[axis]
                    low, high = max(low, min(a, b)), min(high, max(a, b))
                    if high <= low:
                        break
            if high > low:
                intersections[ray, cell] = (high - low) * length
    path_weights = weights @ intersections
    if not np.isfinite(path_weights).all():
        raise ValueError("ray projection exceeds finite arithmetic")
    return {"path_weights_m": path_weights.tolist(), "ray_cell_lengths_m": intersections.tolist(),
            "collection_weights": weights.tolist(), "geometry_applied_once": True,
            "boundary_ownership": "half_open_global_outer_maximum_closed",
            "detector_etendue_applied": False}


def _states(document):
    species = document["species"]
    state_rows = document["states"]
    species_ids = _unique([x["species_id"] for x in species], "species IDs")
    state_ids = _unique([x["state_id"] for x in state_rows], "state IDs")
    species_map = dict(zip(species_ids, species))
    state_map = dict(zip(state_ids, state_rows))
    for row in species:
        _text(row["formula"], "formula")
        if row["kind"] not in {"atom", "molecule"} or type(row["charge_e"]) is not int:
            raise ValueError("species must declare atomic/molecular kind and integer charge")
    represented = set()
    hydrogen_n2 = set()
    for row in state_rows:
        parent = species_map.get(row["species_id"])
        if parent is None:
            raise ValueError("state refers to an unknown canonical species")
        represented.add(row["species_id"])
        _text(row["electronic_state"], "electronic_state")
        if parent["kind"] == "molecule":
            count = parent["vibrational_mode_count"]
            vibrational = row["vibrational_quantum_numbers"]
            if (type(count) is not int or count <= 0 or not isinstance(vibrational, list)
                    or len(vibrational) != count or any(type(x) is not int or x < 0 for x in vibrational)):
                raise ValueError("every retained molecule requires explicit resolved vibrational quantum numbers")
            rotor = parent["rotor_type"]
            keys = {"diatomic": {"J"}, "linear": {"J"}, "symmetric": {"J", "K"},
                    "asymmetric": {"J", "Ka", "Kc"}}.get(rotor)
            rotational = row["rotational_quantum_numbers"]
            if keys is None or not isinstance(rotational, dict) or not keys <= rotational.keys():
                raise ValueError("molecular rotational quantum numbers do not resolve the declared rotor")
            for key in keys:
                value = _number(rotational[key], "rotational quantum number", nonnegative=True)
                if not (2 * value).is_integer():
                    raise ValueError("rotational quantum numbers must be integer or half-integer")
            _text(row["nuclear_spin_symmetry"], "nuclear_spin_symmetry")
        elif parent["formula"] == "H" and parent["charge_e"] == 0:
            n, orbital = row["principal_n"], row["orbital_l"]
            if type(n) is not int or type(orbital) is not int or n <= 0 or not 0 <= orbital < n:
                raise ValueError("neutral H requires resolved principal_n/orbital_l")
            if n == 2:
                hydrogen_n2.add(orbital)
        else:
            _text(row["level_id"], "atomic level_id")
    if represented != set(species_ids):
        raise ValueError("every retained species requires a state mapping")
    if hydrogen_n2 and hydrogen_n2 != {0, 1}:
        raise ValueError("neutral H n=2 mapping must distinguish both 2s and 2p")
    return state_map


def _source_moments(source, state_map):
    patches = source["patches"]
    patch_ids = _unique([x["patch_id"] for x in patches], "source patch IDs")
    for patch in patches:
        _number(patch["area_m2"], "source patch area", positive=True)
        normal = _array(patch["outward_normal"], "outward_normal", (3,))
        if not np.isclose(np.linalg.norm(normal), 1, rtol=0, atol=1e-12):
            raise ValueError("patch outward normal must be unit length")
        if patch["wall_model"] != "stationary_supported":
            raise ValueError("this receipt validator supports only stationary supported walls")
    groups = source["groups"]
    _unique([x["group_id"] for x in groups], "source group IDs")
    summaries = {}
    for group in groups:
        if group["patch_id"] not in patch_ids or group["state_id"] not in state_map:
            raise ValueError("source group uses unmapped patch/state")
        for name in ("reservoir_id", "source_history_id", "source_query_id"):
            _text(group[name], name)
        anchor = _number(group["local_time_anchor_s"], "source local time anchor")
        nodes = group["nodes"]
        _unique([x["node_id"] for x in nodes], "source node IDs")
        values = []
        for node in nodes:
            number = _number(node["physical_particle_count"], "physical_particle_count", positive=True)
            tau = _number(node["time_s"], "node time") - anchor
            temp = _number(node["temperature_k"], "source temperature", positive=True)
            energies = [_number(node[key], key, nonnegative=True) for key in (
                "removal_energy_per_parent_j", "internal_energy_per_parent_j",
                "expected_translational_energy_per_parent_j")]
            # Area and interval scale have already entered physical_particle_count.
            values.append(number * np.array([1, tau, tau * tau, temp, temp * temp, tau * temp, *energies]))
        total = np.sum(values, axis=0)
        if not np.isfinite(total).all():
            raise ValueError("source moments exceed finite arithmetic")
        summaries[group["group_id"]] = {"identity": {key: group[key] for key in (
            "patch_id", "state_id", "reservoir_id", "source_history_id", "source_query_id", "local_time_anchor_s")},
            "moments": dict(zip(MOMENTS, total.tolist()))}
    return summaries


def _source_placement(patches, readbacks, cell_ids, cell_bounds):
    """Check stationary rectangular source supports against receiving gas faces."""
    if any("placement" not in patch for patch in patches) or readbacks is None:
        return False
    ids = _unique([x["patch_id"] for x in readbacks], "receiver placement patch IDs")
    if set(ids) != {x["patch_id"] for x in patches}:
        raise ValueError("receiver placement must retain every source patch")
    receiving = dict(zip(ids, readbacks))
    bounds = np.asarray(cell_bounds, dtype=float)
    for patch in patches:
        placement = patch["placement"]
        origin = _array(placement["origin_m"], "source patch origin", (3,))
        u = _array(placement["span_u_m"], "source patch span u", (3,))
        v = _array(placement["span_v_m"], "source patch span v", (3,))
        u_axes, v_axes = np.flatnonzero(u), np.flatnonzero(v)
        if len(u_axes) != 1 or len(v_axes) != 1 or u_axes[0] == v_axes[0]:
            raise ValueError("source placement adapter requires axis-aligned rectangular patch spans")
        area = float(np.linalg.norm(np.cross(u, v)))
        if not np.isclose(area, patch["area_m2"], rtol=1e-12, atol=0):
            raise ValueError("source patch placement area disagrees with its physical area")
        normal = np.asarray(patch["outward_normal"], dtype=float)
        if not np.allclose(np.cross(u, v) / area, normal, rtol=0, atol=1e-12):
            raise ValueError("source patch span orientation disagrees with its outward normal")
        row = receiving[patch["patch_id"]]
        if row["placement"] != placement:
            raise ValueError("receiver changed the declared source patch placement")
        _text(row["insertion_geometry_evidence"], "receiver insertion geometry evidence")
        receiving_ids = _unique(row["receiving_cell_ids"], "source receiving cell IDs")
        if not set(receiving_ids) <= set(cell_ids):
            raise ValueError("source placement refers to unknown receiving gas cells")
        axis = int(({0, 1, 2} - {int(u_axes[0]), int(v_axes[0])}).pop())
        tangent = [int(u_axes[0]), int(v_axes[0])]
        vertices = np.stack([origin, origin + u, origin + v, origin + u + v])
        patch_lo, patch_hi = vertices.min(axis=0), vertices.max(axis=0)
        covered_area = 0.0
        for cell_id in receiving_ids:
            box = bounds[cell_ids.index(cell_id)]
            face = box[0 if normal[axis] > 0 else 1, axis]
            if face != origin[axis]:
                raise ValueError("source placement does not touch the outward receiving cell face")
            extent = np.minimum(patch_hi[tangent], box[1, tangent]) - np.maximum(patch_lo[tangent], box[0, tangent])
            if np.any(extent <= 0):
                raise ValueError("source placement has no area on a declared receiving gas cell")
            covered_area += float(np.prod(extent))
        if not np.isclose(covered_area, area, rtol=1e-12, atol=0):
            raise ValueError("receiving gas faces do not cover the complete source patch")
    return True


def validate_raikou_handoff(document: dict, artifact_bytes: Mapping[str, bytes]) -> dict:
    """Check hash-bound owner exports; return admission limits without running Raikou.

    A passing manufactured fixture is never a physical receiver receipt. An
    owner-supplied receipt remains attributed to its supplied execution evidence.
    Numerical tolerances are per-quantity absolute tolerances plus one declared
    relative tolerance; no unit-sized floor suppresses small physical channels.
    """
    try:
        return _validate(document, artifact_bytes)
    except KeyError as exc:
        raise ValueError(f"required handoff field missing: {exc.args[0]}") from exc


def prepare_resolved_projection(document: dict, artifact_bytes: Mapping[str, bytes], instrument: dict,
                                *, allow_manufactured_control: bool = False) -> dict:
    """Prepare explicit projector inputs after admission, without evolving any state.

    Effective Gaussian widths are supplied by the instrument/profile owner. The
    checked velocity closure is not automatically converted into a line profile.
    The return contains transition dictionaries, to be constructed as the exact
    projector module's ``ResolvedTransition`` class by the caller.
    """
    report = validate_raikou_handoff(document, artifact_bytes)
    try:
        if instrument["schema"] != "mhdlab-resolved-gaussian-instrument-v1":
            raise ValueError("unsupported instrument/profile schema")
        if instrument["status"] != document["status"]:
            raise ValueError("instrument and population evidence status must match")
        for key in ("effective_profile_justification", "velocity_treatment", "profile_qualification_reference"):
            _text(instrument[key], key)
        if instrument["population_history_sha256"] != report["evidence_sha256"]["populations"]:
            raise ValueError("effective profile is not qualified for this population history")
        if "resolved_radiative_transfer_required" in report["blockers"]:
            raise ValueError("thin projector cannot replace required resolved radiative transfer")
        if document["status"] == "manufactured":
            if allow_manufactured_control is not True:
                raise ValueError("manufactured projection requires explicit allow_manufactured_control=True")
            if instrument["effective_profile_status"] != "manufactured":
                raise ValueError("manufactured effective profiles must retain their status")
        elif report["blockers"]:
            raise ValueError("physical projection is blocked: " + ", ".join(report["blockers"]))
        elif instrument["effective_profile_status"] != "owner_qualified_for_this_gate_and_geometry":
            raise ValueError("physical projection requires owner qualification of effective Gaussian profiles")
        populations = strict_json(artifact_bytes["populations"])
        rates = strict_json(artifact_bytes["rates"])
        edges = _array(instrument["bin_edges_nm"], "instrument bin_edges_nm")
        if edges.ndim != 1 or len(edges) < 2 or np.any(edges <= 0) or np.any(np.diff(edges) <= 0):
            raise ValueError("instrument bins must be positive and strictly increasing")
        widths = _array(np.atleast_1d(instrument["effective_gaussian_fwhm_nm"]), "effective Gaussian FWHM")
        if widths.ndim != 1 or len(widths) not in {1, len(rates["transitions"])} or np.any(widths <= 0):
            raise ValueError("one positive Gaussian FWHM or one per transition is required")
        used_states = {x["upper_state_id"] for x in rates["transitions"]}
        return {"admission": report, "profile_quality_independently_verified": False,
                "velocity_projection_performed_by_adapter": False,
                "instrument_declaration": instrument,
                "projection_kwargs": {
                    "time_s": report["registered_photon_time_s"],
                    "upper_state_density_m3": {key: value for key, value in populations["density_m3"].items()
                                                if key in used_states},
                    "path_weights_m": report["geometry"]["path_weights_m"],
                    "transitions": [{key: row[key] for key in ("transition_id", "upper_state_id",
                                      "vacuum_wavelength_nm", "einstein_a_per_s")} for row in rates["transitions"]],
                    "bin_edges_nm": edges.tolist(), "gate_start_s": report["photon_gate_s"][0],
                    "gate_end_s": report["photon_gate_s"][1],
                    "fwhm_nm": float(widths[0]) if len(widths) == 1 else widths.tolist(),
                    "time_registration": json.dumps(populations["clock_registration"], sort_keys=True),
                    "population_provenance": "hash-bound owner bundle: " + report["evidence_sha256"]["populations"],
                    "geometry_provenance": "mesh/collection SHA256: " + report["evidence_sha256"]["mesh"]
                                           + "/" + report["evidence_sha256"]["collection"],
                    "status": document["status"], "optical_model": "optically_thin"}}
    except KeyError as exc:
        raise ValueError(f"required instrument/profile field missing: {exc.args[0]}") from exc


def _validate(document, artifact_bytes):
    if document["schema"] != "mhdlab-raikou-handoff-v1":
        raise ValueError("unsupported handoff schema")
    status = document["status"]
    if status not in {"manufactured", "owner_supplied"}:
        raise ValueError("status must be manufactured or owner_supplied")
    _text(document["case_id"], "case_id")
    if set(document["artifacts"]) != ARTIFACTS or set(artifact_bytes) != ARTIFACTS:
        raise ValueError("all named evidence artifacts must be supplied exactly once")
    hashes, data = {}, {}
    for name in sorted(ARTIFACTS):
        raw = artifact_bytes[name]
        if not isinstance(raw, bytes):
            raise ValueError("artifact content must be immutable bytes")
        expected = _digest(document["artifacts"][name]["sha256"], name)
        hashes[name] = sha256(raw).hexdigest()
        if hashes[name] != expected:
            raise ValueError(f"{name} evidence hash differs")
        if name != "runtime_output":
            data[name] = strict_json(raw)
    source, receipt, populations = data["source"], data["receipt"], data["populations"]
    blockers = []
    if any(value["status"] != status for value in data.values()):
        raise ValueError("manufactured/owner status cannot change within a bundle")
    if not artifact_bytes["runtime_output"]:
        raise ValueError("receiver runtime-output evidence must not be empty")
    if source["representation"] != "positive_joint_boundary_quadrature":
        raise ValueError("unsupported source representation")
    state_map = _states(data["states"])
    summaries = _source_moments(source, state_map)
    if receipt["source_sha256"] != hashes["source"] or receipt["runtime_output_sha256"] != hashes["runtime_output"]:
        raise ValueError("receiver receipt is not bound to source and runtime-output bytes")
    if receipt["accepted"] is not True:
        raise ValueError("receiver has not accepted this source transaction")
    _text(receipt["transaction_id"], "transaction_id")
    for key in ("duplicate_transaction_policy", "rollback_evidence", "restart_evidence"):
        _text(receipt[key], key)
    execution = receipt["execution"]
    if type(execution["exit_code"]) is not int or execution["exit_code"] != 0:
        raise ValueError("receiver execution did not succeed")
    if status == "owner_supplied":
        if execution["kind"] != "runtime_receipt":
            raise ValueError("owner-supplied acceptance requires runtime receipt evidence")
        _digest(execution["executable_sha256"], "receiver executable")
        _text(execution["source_revision"], "receiver source revision")
        _text(execution["command"], "receiver command")
    elif execution["kind"] != "manufactured_control":
        raise ValueError("manufactured status requires manufactured execution labeling")

    tolerance = document["conservation_tolerance"]
    rtol = _number(tolerance["relative"], "relative conservation tolerance", nonnegative=True)
    if rtol > 1e-3:
        raise ValueError("relative tolerance is too loose for source admission")
    absolute = tolerance["absolute"]
    used_residuals = []

    def close(actual, expected, quantity):
        actual = _number(actual, quantity)
        expected = _number(expected, quantity)
        atol = _number(absolute[quantity], f"absolute tolerance {quantity}", nonnegative=True)
        if quantity == "released_number" and atol != 0:
            raise ValueError("released_number requires zero absolute tolerance to retain rare positive pools")
        residual = actual - expected
        allowed = atol + rtol * max(abs(actual), abs(expected))
        if abs(residual) > allowed:
            raise ValueError(f"receiver {quantity} does not conserve its bound source/ledger")
        used_residuals.append({"quantity": quantity, "residual": residual, "allowance": allowed})

    received = receipt["groups"]
    _unique([x["group_id"] for x in received], "received group IDs")
    if {x["group_id"] for x in received} != summaries.keys():
        raise ValueError("receiver must retain every source component, including rare positive components")
    for row in received:
        expected = summaries[row["group_id"]]
        if row["identity"] != expected["identity"]:
            raise ValueError("receiver changed the source component identity or time anchor")
        if set(row["moments"]) != set(MOMENTS):
            raise ValueError("receiver must retain all joint number/time/temperature/energy moments")
        for key in MOMENTS:
            close(row["moments"][key], expected["moments"][key], key)
        number = _number(row["moments"]["released_number"], "received physical count", positive=True)
        variance_t = (row["moments"]["local_time_second_moment_s2"] / number
                      - (row["moments"]["local_time_first_moment_s"] / number) ** 2)
        variance_T = (row["moments"]["temperature_second_moment_k2"] / number
                      - (row["moments"]["temperature_first_moment_k"] / number) ** 2)
        covariance = (row["moments"]["local_time_temperature_moment_k_s"] / number
                      - row["moments"]["local_time_first_moment_s"]
                      * row["moments"]["temperature_first_moment_k"] / number ** 2)
        close(row["covariance_local_time_temperature_k_s"], covariance, "covariance_local_time_temperature_k_s")
        # Preserve expected source translation and actual applied gas energy separately.
        ledger = row["stationary_energy_ledger"]
        for key in ("gas_added_translational_j", "gas_added_internal_j", "adsorption_reference_change_j"):
            _number(ledger[key], key, nonnegative=True)
        close(ledger["gas_added_internal_j"], row["moments"]["internal_energy_j"], "internal_energy_j")
        close(ledger["adsorption_reference_change_j"], row["moments"]["binding_energy_j"], "binding_energy_j")
        close(ledger["wall_total_change_j"], ledger["wall_thermal_change_j"]
              + ledger["adsorption_reference_change_j"], "energy_ledger_j")
        close(ledger["wall_total_change_j"] + ledger["gas_added_translational_j"]
              + ledger["gas_added_internal_j"], 0, "energy_ledger_j")
        before = _number(ledger["wall_thermal_before_j"], "wall thermal account before", nonnegative=True)
        after = _number(ledger["wall_thermal_after_j"], "wall thermal account after", nonnegative=True)
        close(after - before, ledger["wall_thermal_change_j"], "energy_ledger_j")
        gas_p = _array(row["gas_added_momentum_kg_m_per_s"], "gas momentum", (3,))
        support_p = _array(row["external_support_impulse_kg_m_per_s"], "support impulse", (3,))
        for value in gas_p - support_p:
            close(value, 0, "momentum_ledger_kg_m_per_s")
        row_variances = (variance_t, variance_T)
        if not all(np.isfinite(row_variances)):
            raise ValueError("receiver covariance exceeds finite arithmetic")

    expected_pools = {}
    for summary in summaries.values():
        identity = summary["identity"]
        key = (identity["reservoir_id"], identity["patch_id"])
        expected_pools[key] = expected_pools.get(key, 0.0) + summary["moments"]["released_number"]
    inventory = receipt["inventory_ledgers"]
    actual_pools = [(x["reservoir_id"], x["patch_id"]) for x in inventory]
    if len(set(actual_pools)) != len(actual_pools) or set(actual_pools) != set(expected_pools):
        raise ValueError("each owning reservoir/patch inventory must be debited exactly once")
    for pool, key in zip(inventory, actual_pools):
        before = _number(pool["before_parent_count"], "inventory before", nonnegative=True)
        after = _number(pool["after_parent_count"], "inventory after", nonnegative=True)
        added = _number(pool["gas_inserted_parent_count"], "gas inserted parent count", positive=True)
        close(before - after, expected_pools[key], "released_number")
        close(added, expected_pools[key], "released_number")

    mesh, collection = data["mesh"], data["collection"]
    if mesh["units"] != "m" or mesh["fixed_over_gate"] is not True:
        raise ValueError("mesh must use metres and be fixed over the photon gate")
    _text(mesh["simulation_provenance"], "simulation geometry provenance")
    cells = _unique(mesh["cell_ids"], "mesh cell IDs")
    _array(mesh["cell_bounds_m"], "cell bounds", (len(cells), 2, 3))
    if collection["mesh_sha256"] != hashes["mesh"]:
        raise ValueError("collection refers to a different simulation mesh")
    _text(collection["optics_provenance"], "optics provenance")
    geometry = collected_path_weights(mesh["cell_bounds_m"], collection["ray_start_m"],
                                      collection["ray_end_m"], collection["radiance_collection_weights"])
    placement_checked = _source_placement(source["patches"], receipt.get("patch_placements"), cells,
                                          mesh["cell_bounds_m"])
    if not placement_checked:
        blockers.append("physical_source_placement_unverified")
    if (populations["source_receipt_sha256"] != hashes["receipt"]
            or populations["state_mapping_sha256"] != hashes["states"]
            or populations["mesh_sha256"] != hashes["mesh"]
            or populations["rates_sha256"] != hashes["rates"]):
        raise ValueError("population history uses a different receipt, mesh, state map or rates")
    times = _array(populations["source_time_s"], "source population time")
    if times.ndim != 1 or len(times) < 2 or np.any(np.diff(times) <= 0):
        raise ValueError("population times must be strictly increasing")
    registration = populations["clock_registration"]
    for key in ("source_clock_id", "photon_clock_id", "provenance"):
        _text(registration[key], key)
    offset = _number(registration["photon_minus_source_time_s"], "clock offset")
    photon_times = times + offset
    if not np.isfinite(photon_times).all() or np.any(np.diff(photon_times) <= 0):
        raise ValueError("registered photon clock is not finite/increasing")
    gate = _array(populations["photon_gate_s"], "photon gate", (2,))
    if gate[1] <= gate[0] or gate[0] < photon_times[0] or gate[1] > photon_times[-1]:
        raise ValueError("photon gate must be fully contained in the registered population history")
    if populations["cell_ids"] != cells or set(populations["density_m3"]) != set(state_map):
        raise ValueError("population state/cell identities must be complete and preserve mesh ordering")
    shape = (len(times), len(cells))
    for state, values in populations["density_m3"].items():
        _array(values, f"density {state}", shape, nonnegative=True)
    rates = data["rates"]
    if rates["state_mapping_sha256"] != hashes["states"]:
        raise ValueError("rates use a different canonical state map")
    _text(rates["data_revision"], "rate data revision")
    coverage = rates["channel_coverage"]
    if set(coverage) != CHANNELS:
        raise ValueError("rate coverage must explicitly address every required channel family")
    for channel, evidence in coverage.items():
        if evidence["status"] not in {"evaluated", "not_applicable", "missing"}:
            raise ValueError("invalid rate-channel coverage status")
        _text(evidence["evidence"], "rate-channel evidence")
        if evidence["status"] == "missing":
            blockers.append(f"missing_rate_family:{channel}")
    transitions = rates["transitions"]
    transition_ids = _unique([x["transition_id"] for x in transitions], "transition IDs")
    for line in transitions:
        if line["upper_state_id"] not in state_map or line["lower_state_id"] not in state_map:
            raise ValueError("radiative transition refers to an unmapped resolved state")
        _number(line["einstein_a_per_s"], "Einstein A", positive=True)
        _number(line["vacuum_wavelength_nm"], "vacuum wavelength", positive=True)
    velocity = populations["velocity_closures"]
    if set(velocity) != {x["upper_state_id"] for x in transitions}:
        raise ValueError("every emitting state requires its own velocity distribution or declared closure")
    for state, closure in velocity.items():
        _text(closure["justification"], "velocity closure justification")
        if closure["model"] == "isotropic_maxwellian":
            temperature = _array(closure["translational_temperature_k"], "neutral translational temperature", shape)
            if np.any(temperature <= 0):
                raise ValueError("velocity closure temperature must be positive")
            _array(closure["mean_velocity_m_per_s"], "mean velocity", (*shape, 3))
        elif closure["model"] == "resolved_velocity_grid":
            nodes = _array(closure["velocity_nodes_m_per_s"], "velocity nodes")
            if nodes.ndim != 2 or nodes.shape[1] != 3:
                raise ValueError("velocity nodes must have shape (velocity_node,3)")
            probability = _array(closure["number_fraction"], "velocity number fraction",
                                 (*shape, len(nodes)), nonnegative=True)
            if not np.allclose(np.sum(probability, axis=-1), 1, rtol=0, atol=1e-12):
                raise ValueError("resolved velocity number fractions must already sum to one")
        else:
            raise ValueError("unsupported emitting-state velocity closure")
    electrons = data["electrons"]
    if electrons["source_time_s"] != populations["source_time_s"] or electrons["cell_ids"] != cells:
        raise ValueError("electron excitation history must share the resolved time/cell coordinates")
    _text(electrons["eedf_model"], "electron EEDF model")
    _text(electrons["provenance"], "electron excitation provenance")
    _array(electrons["density_m3"], "electron density", shape, nonnegative=True)
    _array(electrons["mean_energy_ev"], "electron mean energy", shape, nonnegative=True)
    if populations["electron_history_sha256"] != hashes["electrons"]:
        raise ValueError("population history is not bound to electron excitation history")
    opacity = data["opacity"]
    if opacity["populations_sha256"] != hashes["populations"] or opacity["collection_sha256"] != hashes["collection"]:
        raise ValueError("optical-depth assessment uses different populations or collection geometry")
    _text(opacity["method"], "optical depth method")
    if set(opacity["transition_ids"]) != set(transition_ids) or len(opacity["transition_ids"]) != len(transition_ids):
        raise ValueError("optical-depth assessment must cover every retained transition")
    tau = _number(opacity["maximum_optical_depth"], "maximum optical depth", nonnegative=True)
    allowed_tau = _number(opacity["admission_limit"], "optical depth admission limit", positive=True)
    if allowed_tau >= 1:
        raise ValueError("thin-emission optical-depth admission limit must be less than one")
    thin_budget = _number(opacity["relative_radiance_error_budget"], "thin radiance error budget", positive=True)
    if thin_budget > 0.01:
        raise ValueError("thin-emission relative radiance error budget must be at most 1 percent")
    attenuation_bound = float(-np.expm1(-tau))
    if (opacity["optical_model"] != "optically_thin" or tau > min(allowed_tau, 0.01)
            or attenuation_bound > thin_budget):
        blockers.append("resolved_radiative_transfer_required")
    if status == "manufactured":
        blockers.append("manufactured_evidence_only")
    return {"schema": "mhdlab-raikou-handoff-validation-v1", "case_id": document["case_id"],
            "status": status, "evidence_sha256": hashes, "receiver_receipt_structurally_valid": True,
            "physical_receiver_receipt_present": status == "owner_supplied",
            "independent_runtime_execution_verified": False,
            "thin_projection_admissible_under_owner_declarations": not blockers,
            "source_placement_checked_against_receiving_faces": placement_checked,
            "thin_emission_assessment": {"maximum_optical_depth": tau,
                "hard_optical_depth_limit": 0.01, "declared_relative_radiance_error_budget": thin_budget,
                "fixed_emissivity_absorption_attenuation_bound": attenuation_bound},
            "blockers": blockers, "geometry": geometry, "registered_photon_time_s": photon_times.tolist(),
            "photon_gate_s": gate.tolist(), "source_group_count": len(summaries),
            "resolved_state_count": len(state_map), "transition_count": len(transitions),
            "conservation_checks": used_residuals,
            "limits": ["No external model was executed by this validation.",
                       "Hashes bind evidence; its physical origin and modeling closure remain owner declarations.",
                       "No physical rate completeness, absolute calibration or neutral parameter inference is established.",
                       "Source expected translational energy and actual receiver gas energy remain distinct.",
                       "Fixed axis-aligned cells, finite rays and stationary supported walls only."]}
