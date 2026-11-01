"""Conservative reference-energy mapping for the dedicated native heat pilot.

The native scheme owns mass-lumped *nodal* h(T), not quadrature h(Q1 T).
Canonical nodal/cell arrays are (y,x), x fastest. Solid arrays are (z,y,x).
This initial adapter requires identical XY grids and an arbitrary axial split.
It never commits state, evolves species or authorizes deformed field geometry.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import re

import numpy as np

from .solid_foil_coupling import project_exchange_to_thermal_xy
from .solid_foil_geometry import extrude_extensive_xy
from .solid_foil_mechanics import mesh_geometry_sha256


RTOL = 5e-13


def _finite(value, name):
    raw = np.asarray(value)
    if np.ma.isMaskedArray(value) or raw.dtype.kind not in "fiu" or not raw.size:
        raise ValueError(f"{name} must be a nonempty real array")
    result = np.asarray(value, dtype=float)
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must be finite")
    return result


def _edges(value, name):
    result = _finite(value, name)
    if result.ndim != 1 or len(result) < 2 or np.any(np.diff(result) <= 0):
        raise ValueError(f"{name} must be increasing cell edges")
    return result


def _immutable(value):
    array = np.asarray(value)
    return np.frombuffer(np.ascontiguousarray(array).tobytes(), dtype=array.dtype).reshape(array.shape)


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _sha(value, name):
    if not isinstance(value, str) or re.fullmatch("[0-9a-f]{64}", value) is None:
        raise ValueError(f"{name} must be a lowercase SHA256")


def native_reference_mesh_sha256(x_edges_m, y_edges_m, height_m):
    x, y = _edges(x_edges_m, "native x edges"), _edges(y_edges_m, "native y edges")
    height = _finite(height_m, "native height")
    if height.ndim != 0 or height <= 0:
        raise ValueError("positive scalar native axial extent required")
    digest = hashlib.sha256(b"native-Q1-reference-xy-full-height-v1")
    for value in (x, y, height):
        digest.update(_json(list(value.shape)))
        digest.update(np.asarray(value, dtype="<f8").tobytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class NativeHeatStageIdentity:
    """Both stages use the accepted base revision; stage selects t or t+dt."""

    native_mesh_revision: str
    solid_mesh_revision: str
    state_revision: int
    interval_s: tuple[float, float]
    stage_time_s: float
    stage: str
    material_sha256: str
    source_sha256: str
    caloric_sha256: str
    native_reference_mesh_sha256: str
    solid_reference_mesh_sha256: str

    def __post_init__(self):
        for name in ("native_mesh_revision", "solid_mesh_revision"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip() or value != value.strip():
                raise ValueError(f"explicit {name} required")
        if type(self.state_revision) is not int or self.state_revision < 0:
            raise ValueError("nonnegative accepted base revision required")
        interval = _finite(self.interval_s, "accepted interval")
        time = _finite(self.stage_time_s, "stage time")
        if interval.shape != (2,) or interval[0] < 0 or interval[1] <= interval[0] or time.ndim != 0:
            raise ValueError("positive accepted interval and scalar stage time required")
        if self.stage not in ("accepted", "trial") or time != interval[int(self.stage == "trial")]:
            raise ValueError("accepted stage must be t; trial stage must be t+dt")
        object.__setattr__(self, "interval_s", tuple(float(x) for x in interval))
        object.__setattr__(self, "stage_time_s", float(time))
        for name in ("material_sha256", "source_sha256", "caloric_sha256", "native_reference_mesh_sha256", "solid_reference_mesh_sha256"):
            _sha(getattr(self, name), name)


def make_native_heat_stage_identity(mesh, *, native_x_edges_m, native_y_edges_m, native_height_m, **identities):
    return NativeHeatStageIdentity(
        native_reference_mesh_sha256=native_reference_mesh_sha256(native_x_edges_m, native_y_edges_m, native_height_m),
        solid_reference_mesh_sha256=mesh_geometry_sha256(mesh), **identities)


def _check_identity(mesh, x, y, height, identity, expected):
    if not isinstance(identity, NativeHeatStageIdentity) or identity != expected:
        raise ValueError("native/solid source, material, mesh, revision or stage identity mismatch")
    if (identity.native_reference_mesh_sha256 != native_reference_mesh_sha256(x, y, height)
            or identity.solid_reference_mesh_sha256 != mesh_geometry_sha256(mesh)):
        raise ValueError("reference mesh differs from immutable geometry identity")
    if not np.array_equal(mesh.x_edges_m, x) or not np.array_equal(mesh.y_edges_m, y):
        raise ValueError("initial native mapping requires identical complete XY grids; no extrapolation")
    if float(mesh.z_edges_m[-1] - mesh.z_edges_m[0]) != float(height):
        raise ValueError("native and solid full axial extents must match exactly")
    expected_volumes = (np.diff(mesh.z_edges_m)[:, None, None]
                        * np.diff(y)[None, :, None] * np.diff(x)[None, None, :])
    if (tuple(mesh.cell_shape) != expected_volumes.shape or len(mesh.cells) != expected_volumes.size
            or not np.allclose(mesh.reference_volumes_m3, expected_volumes.ravel(), rtol=RTOL, atol=0)):
        raise ValueError("solid reference volumes or cell ordering inconsistent with tensor grid")


def native_lumped_cell_energy_J(nodal_temperature_K, *, x_edges_m, y_edges_m, height_m, caloric):
    """Independent exact-discretization reconstruction, not h(mean T).

    Each rectangular cell owns one quarter of each corner's reference energy.
    Summing cells gives native dual-cell nodal weights, including half/quarter
    boundary contributions. Nonlinear h is evaluated at the nodes themselves.
    """
    x, y = _edges(x_edges_m, "x edges"), _edges(y_edges_m, "y edges")
    native_reference_mesh_sha256(x, y, height_m)
    temperature = _finite(nodal_temperature_K, "native nodal temperature")
    if temperature.shape != (len(y), len(x)):
        raise ValueError("native nodal temperature must be canonical (ny+1,nx+1), x fastest")
    energy = _finite(caloric.energy_J_m3(temperature), "nodal reference caloric energy")
    if energy.shape != temperature.shape:
        raise ValueError("caloric callback changed nodal array shape")
    mean = (energy[:-1, :-1] + energy[:-1, 1:] + energy[1:, 1:] + energy[1:, :-1]) / 4
    volume = np.diff(y)[:, None] * np.diff(x)[None, :] * float(height_m)
    return _immutable(mean * volume)


@dataclass(frozen=True)
class NativeSolidHeatProjection:
    identity: NativeHeatStageIdentity
    cell_temperature_K: np.ndarray
    cell_base_energy_J: np.ndarray
    native_cell_base_energy_J: np.ndarray
    receipt_bytes: bytes

    def receipt(self):
        return json.loads(self.receipt_bytes)


def project_native_heat_to_solid(mesh, native_temperature_K, *, native_x_edges_m, native_y_edges_m,
                                native_height_m, caloric, identity, expected_identity,
                                native_cell_base_energy_J, native_stored_energy_J,
                                admissible_temperature_range_K):
    """Conserve native lumped caloric energy, then invert cell-average h for T.

    The provided native cell-energy getter is authoritative; independent corner
    h reconstruction rejects discretization, callback and ordering mismatches.
    Native full-height cell J are distributed by dz/height, exactly once. Native
    energy remains the authoritative ledger; the inverse's representable-T
    residual is checked against an explicit capacity/temperature roundoff bound.
    This projection reduces subcell variation and is not a nonlinear Q1 energy
    quadrature or a mechanical convergence claim.
    """
    x, y = _edges(native_x_edges_m, "native x edges"), _edges(native_y_edges_m, "native y edges")
    _check_identity(mesh, x, y, native_height_m, identity, expected_identity)
    if caloric.receipt.get("source_seal_sha256") != identity.caloric_sha256:
        raise ValueError("caloric callback source seal differs from stage identity")
    temperature = _finite(native_temperature_K, "native nodal temperature")
    bounds = _finite(admissible_temperature_range_K, "admissible temperature range")
    if bounds.shape != (2,) or not 0 < bounds[0] < bounds[1]:
        raise ValueError("positive ordered solid temperature bounds required")
    if np.any(temperature < bounds[0]) or np.any(temperature > bounds[1]):
        raise ValueError("native nodal temperature outside solid admission; averaging cannot hide hot nodes")
    reconstructed = native_lumped_cell_energy_J(temperature, x_edges_m=x, y_edges_m=y,
                                                height_m=native_height_m, caloric=caloric)
    native_cells = _finite(native_cell_base_energy_J, "native cell base energy")
    total = _finite(native_stored_energy_J, "native stored energy")
    if native_cells.shape != reconstructed.shape or total.ndim != 0:
        raise ValueError("native cell energy must be canonical (ny,nx) with scalar stored energy")
    nodal_energy = np.abs(_finite(caloric.energy_J_m3(temperature), "nodal energy"))
    volume_xy = np.diff(y)[:, None] * np.diff(x)[None, :] * float(native_height_m)
    local_scale = (nodal_energy[:-1, :-1] + nodal_energy[:-1, 1:] + nodal_energy[1:, 1:] + nodal_energy[1:, :-1]) / 4 * volume_xy
    native_local_error = np.abs(native_cells - reconstructed)
    if np.any(native_local_error > RTOL * local_scale):
        raise ValueError("native cell energy differs from exact lumped nodal caloric reconstruction")
    total_bound = RTOL * float(local_scale.sum())
    if abs(float(native_cells.sum()) - float(total)) > total_bound:
        raise ValueError("native cell energy sum differs from stored-energy ledger")
    energy = extrude_extensive_xy(native_cells, mesh.z_edges_m)
    volumes = np.asarray(mesh.reference_volumes_m3).reshape(mesh.cell_shape)
    cell_temperature = _finite(caloric.temperature_K(energy / volumes), "inverted cell temperature")
    if cell_temperature.shape != energy.shape or np.any(cell_temperature < bounds[0]) or np.any(cell_temperature > bounds[1]):
        raise ValueError("caloric inverse violates solid range or cell shape")
    roundtrip_density = _finite(caloric.energy_J_m3(cell_temperature), "cell caloric roundtrip")
    if roundtrip_density.shape != energy.shape:
        raise ValueError("caloric roundtrip changed cell array shape")
    roundtrip = roundtrip_density * volumes
    capacity = _finite(caloric.capacity_J_m3_K(cell_temperature), "cell caloric capacity")
    if capacity.shape != energy.shape or np.any(capacity <= 0):
        raise ValueError("positive shape-preserving caloric capacity required")
    inverse_bound = (256 * np.finfo(float).eps * volumes * capacity * np.abs(cell_temperature)
                     + RTOL * np.abs(energy))
    if np.any(np.abs(roundtrip - energy) > inverse_bound):
        raise ValueError("caloric inverse roundtrip violates energy conservation bound")
    input_digest = hashlib.sha256(_json(asdict(identity)) + np.asarray(temperature, dtype="<f8").tobytes()
                                 + np.asarray(native_cells, dtype="<f8").tobytes()).hexdigest()
    receipt = {
        "schema": "native-lumped-heat-to-solid-v1", "identity": asdict(identity), "input_sha256": input_digest,
        "caloric_callback_receipt_sha256": hashlib.sha256(_json(dict(caloric.receipt))).hexdigest(),
        "caloric_callback": dict(caloric.receipt),
        "native_stored_energy_J": float(total), "solid_cell_base_energy_J": float(energy.sum()),
        "solid_caloric_roundtrip_energy_J": float(roundtrip.sum()),
        "native_reconstruction_max_cell_error_J": float(native_local_error.max()),
        "native_total_roundoff_bound_J": total_bound,
        "inverse_roundtrip_max_cell_error_J": float(np.abs(roundtrip - energy).max()),
        "inverse_roundtrip_total_error_J": float(roundtrip.sum() - energy.sum()),
        "inverse_roundtrip_total_bound_J": float(inverse_bound.sum()),
        "temperature_range_K": [float(cell_temperature.min()), float(cell_temperature.max())],
        "discretization": "native lumped corner h(T), full-height cell J, axial fractions, caloric cell inverse",
        "cell_order": "z,y,x; x fastest", "depth_convention": "native cell J already includes full height; split once by dz/height",
        "accepted_state_mutated": False, "physical_EOS_ready": False, "field_handoff_authorized": False,
    }
    return NativeSolidHeatProjection(identity, _immutable(cell_temperature), _immutable(energy),
                                     _immutable(native_cells), _json(receipt))


def signed_cell_exchange_to_native(mesh, packet, *, native_x_edges_m, native_y_edges_m, native_height_m,
                                   identity, expected_identity):
    """Verify a completed mechanical trial and emit dedicated signed native cell J.

    No classification as authored/Joule heat, no absolute value, no second depth
    multiplier, no particle/source mutation and no native state commit occur.
    """
    x, y = _edges(native_x_edges_m, "native x edges"), _edges(native_y_edges_m, "native y edges")
    _check_identity(mesh, x, y, native_height_m, identity, expected_identity)
    if identity.stage != "trial":
        raise ValueError("mechanical exchange belongs to trial interval-end stage")
    if (packet.get("schema") != "solid-thermoelastic-cell-exchange-v2"
            or packet.get("material_sha256") != identity.material_sha256
            or packet.get("source_sha256") != identity.source_sha256
            or packet.get("caloric_model", {}).get("source_seal_sha256") != identity.caloric_sha256
            or tuple(packet.get("interval_s", ())) != identity.interval_s
            or packet.get("stage_time_s") != identity.stage_time_s
            or packet.get("conditional_thermodynamic_exchange_ready") is not True
            or packet.get("species_mutation") != "none"
            or packet.get("species_event_ids") != [] or packet.get("species_parent_ids") != []):
        raise ValueError("mechanical exchange material/source/stage/caloric/species identity mismatch")
    mapped = project_exchange_to_thermal_xy(mesh, packet, target_x_edges_m=x, target_y_edges_m=y,
                                          expected_revision=identity.state_revision,
                                          expected_mesh_revision=identity.solid_mesh_revision)
    full = _finite(packet["cell_heat_to_thermal_J"], "signed mechanical cell J")
    collapsed = full.reshape(mesh.cell_shape).sum(axis=0)
    scale = np.abs(full).reshape(mesh.cell_shape).sum(axis=0)
    if np.any(np.abs(mapped["cell_heat_J"] - collapsed) > RTOL * scale):
        raise ValueError("same-grid native signed cell map changed local energy")
    if abs(float(collapsed.sum()) - float(full.sum())) > RTOL * float(np.abs(full).sum()):
        raise ValueError("native signed cell map changed global signed energy")
    # On the mandated identical grid retain the direct axial sum, avoiding
    # unnecessary divide/multiply rounding in an identity XY remap.
    return dict(mapped, cell_heat_J=_immutable(collapsed), source_total_J=float(full.sum()),
                target_total_J=float(collapsed.sum()), identity=asdict(identity),
                native_reference_mesh_sha256=identity.native_reference_mesh_sha256,
                native_mesh_revision=identity.native_mesh_revision, material_sha256=identity.material_sha256,
                source_sha256=identity.source_sha256, caloric_sha256=identity.caloric_sha256,
                signed_cell_load_class="reversible_thermoelastic_internal_exchange",
                native_state_mutated=False, physical_field_handoff_authorized=False)
