"""Sealed q1 thermal-v2 callbacks, including the reviewed 298.2-300 K interval.

This is a separate adapter revision. It never relaxes v1 guards or changes the
owner's source. Reference energy remains zero at 300 K; negative cold energy is
valid. Full thermal support ends at 1800 K; mechanics independently stops at
773.15 K and its smooth-segment energy guards still apply. No vapor, mass loss,
mixture transport, or compressible physical EOS is supplied by this adapter.
"""
from __future__ import annotations

import ctypes
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


SOURCE_SEAL_SHA256 = "1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d"
ADAPTER_REVISION = "ss304-caloric-v2-q1-cold298p2-reference300-v1"


class _ValuesV2(ctypes.Structure):
    # ss304v2::Values is an alias of the unchanged ss304::Values ABI.
    _fields_ = [(name, ctypes.c_double) for name in (
        "temperature_k", "h_j_kg", "u_model_j_kg", "cp_sensible_j_kg_k",
        "apparent_heat_capacity_j_kg_k", "k_w_m_k", "resistivity_ohm_m",
        "sigma_s_m", "liquid_fraction", "energy_j_m3_reference", "capacity_j_m3_reference_k",
    )] + [("phase", ctypes.c_int), ("synthesized_high_temperature", ctypes.c_int)]


class _EnergySupportV2(ctypes.Structure):
    _fields_ = [(name, ctypes.c_double) for name in (
        "minimum_j", "maximum_j", "cold_margin_j", "hot_margin_j",
    )] + [("range_status", ctypes.c_int)]


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _positive(value, name):
    a = np.asarray(value, dtype=float)
    if a.ndim != 0 or not np.isfinite(a) or a <= 0:
        raise ValueError(f"{name} must be a positive finite scalar")
    return float(a)


class SealedNominalCaloricV2:
    """Read-only C ABI binding to explicit q1 nominal_1800.

    Run inside the frozen userspace. Besides source/library hashes, construction
    checks the owner's runtime-binding files as seen by the executing process.
    rho_ref is accepted metal mass/reference volume, not current density. This
    fixed-inventory adapter is invalid after unhandled mass transfer.
    """

    def __init__(self, source_seal_path, *, expected_source_seal_sha256,
                 reference_density_kg_m3, qualification_mode, runtime_bindings_path=None):
        if qualification_mode != "nominal_1800":
            raise ValueError("v2 adapter requires explicit q1 nominal_1800")
        if expected_source_seal_sha256 != SOURCE_SEAL_SHA256:
            raise ValueError("adapter revision requires the reviewed thermal-v2 source seal")
        self._rho = _positive(reference_density_kg_m3, "reference density")
        seal_path = Path(source_seal_path)
        if _sha(seal_path) != expected_source_seal_sha256:
            raise ValueError("thermal-v2 source seal SHA256 mismatch")
        seal = json.loads(seal_path.read_text())
        if seal.get("schema") != "ss304-thermal-v2-source-seal":
            raise ValueError("unexpected thermal-v2 source seal schema")
        root = Path(seal["source_root"])
        required = {"include/ss304_thermal_v2.hpp", "include/sealed1800.hpp",
                    "c_api.cpp", "data/scenario-v2.json"}
        if not required.issubset(seal["files"]):
            raise ValueError("thermal-v2 transitive implementation/scenario seal is incomplete")
        for name, digest in seal["files"].items():
            relative = Path(name)
            if relative.is_absolute() or ".." in relative.parts or _sha(root/relative) != digest:
                raise ValueError(f"thermal-v2 sealed source differs: {name}")
        library = Path(seal["compiled_library"]["path"])
        if _sha(library) != seal["compiled_library"]["sha256"]:
            raise ValueError("thermal-v2 compiled library differs from sealed source")
        runtime_path = (Path(runtime_bindings_path) if runtime_bindings_path is not None
                        else seal_path.parent/"runtime-bindings.json")
        if _sha(runtime_path) != seal["runtime_bindings_sha256"]:
            raise ValueError("thermal-v2 runtime bindings SHA256 mismatch")
        runtime = json.loads(runtime_path.read_text())
        if runtime.get("schema") != "ss304-v2-runtime-bindings":
            raise ValueError("unexpected thermal-v2 runtime bindings schema")
        for item in (runtime["compiler"], runtime["python_verifier"]):
            if _sha(item["path"]) != item["sha256"]:
                raise ValueError("thermal-v2 frozen compiler/Python runtime differs")
        for path, digest in runtime["shared_libraries_sha256"].items():
            if _sha(path) != digest:
                raise ValueError(f"thermal-v2 frozen library runtime differs: {path}")
        if _sha(sys.executable) != runtime["python_verifier"]["sha256"]:
            raise ValueError("execute thermal-v2 adapter with the frozen Python interpreter")

        self._lib = ctypes.CDLL(str(library))
        self._lib.ss304v2_evaluate.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double,
                                             ctypes.POINTER(_ValuesV2)]
        self._lib.ss304v2_evaluate.restype = ctypes.c_int
        self._lib.ss304v2_inverse.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double,
                                            ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
        self._lib.ss304v2_inverse.restype = ctypes.c_int
        self._lib.ss304v2_energy_support.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double,
                                                    ctypes.c_double, ctypes.POINTER(_EnergySupportV2)]
        self._lib.ss304v2_energy_support.restype = ctypes.c_int
        self._lib.ss304v2_error.argtypes = []
        self._lib.ss304v2_error.restype = ctypes.c_char_p
        # Constructor/domain qualification runs through the owner, not duplicated fits.
        cold, zero, hot = (self._evaluate(t) for t in (298.2, 300., 1800.))
        if cold.u_model_j_kg >= 0 or zero.u_model_j_kg != 0 or hot.phase != 2:
            raise ValueError("thermal-v2 q1 ABI/reference-energy contract mismatch")
        receipt = {
            "adapter_revision": ADAPTER_REVISION,
            "source_seal_sha256": expected_source_seal_sha256,
            "source_root": str(root),
            "library_sha256": seal["compiled_library"]["sha256"],
            "scenario_sha256": seal["files"]["data/scenario-v2.json"],
            "runtime_bindings_sha256": seal["runtime_bindings_sha256"],
            "qualification_mode": "nominal_1800",
            "qualification_id": 1,
            "temperature_range_K": [298.2, 1800.],
            "reference_density_kg_m3": self._rho,
            "energy_zero_temperature_K": 300.,
            "minimum_specific_energy_J_kg": cold.u_model_j_kg,
            "energy_basis": "rho_ref*u_model(T), fixed-reference-volume Cp-as-Cv approximation",
            "negative_cold_reference_energy_valid": True,
            "inventory_basis": "fixed accepted metal mass/reference volume; no mass-transfer binding",
            "physical_EOS": False,
            "mechanics_range_extended": False,
            "mechanical_upper_temperature_K": 773.15,
            "phase_or_vapor_mechanics_supplied": False,
        }
        self._receipt_bytes = json.dumps(receipt, sort_keys=True, allow_nan=False).encode()

    @property
    def rho(self):
        return self._rho

    @property
    def receipt(self):
        return json.loads(self._receipt_bytes)

    def _error(self):
        message = self._lib.ss304v2_error()
        return message.decode() if message else "thermal-v2 owner rejected input"

    def _evaluate(self, t):
        values = _ValuesV2()
        if self._lib.ss304v2_evaluate(1, self._rho, float(t), ctypes.byref(values)):
            raise ValueError(self._error())
        if any(not np.isfinite(getattr(values, name)) for name, typ in values._fields_
               if typ is ctypes.c_double):
            raise ValueError("nonfinite thermal-v2 owner result")
        return values

    def evaluate_T(self, temperature_K):
        """Scalar intrinsic owner values; phase flags remain explicit."""
        t = np.asarray(temperature_K, dtype=float)
        if t.ndim != 0 or not np.isfinite(t):
            raise ValueError("evaluate_T requires one finite scalar temperature")
        values = self._evaluate(float(t))
        return {name: getattr(values, name) for name, _ in values._fields_}

    def _vector(self, function, value):
        array = np.asarray(value, dtype=float)
        if array.size == 0 or not np.isfinite(array).all():
            raise ValueError("finite nonempty caloric argument required")
        unique, index = np.unique(array, return_inverse=True)
        result = np.array([function(float(x)) for x in unique])[index].reshape(array.shape)
        return float(result) if result.ndim == 0 else result

    def energy_J_m3(self, temperature_K):
        return self._vector(lambda t: self._evaluate(t).energy_j_m3_reference, temperature_K)

    def capacity_J_m3_K(self, temperature_K):
        return self._vector(lambda t: self._evaluate(t).capacity_j_m3_reference_k, temperature_K)

    def conductivity_W_m_K(self, temperature_K):
        return self._vector(lambda t: self._evaluate(t).k_w_m_k, temperature_K)

    def electrical_conductivity_S_m(self, temperature_K):
        return self._vector(lambda t: self._evaluate(t).sigma_s_m, temperature_K)

    def temperature_K(self, energy_J_m3):
        def inverse(energy):
            result = ctypes.c_double()
            if self._lib.ss304v2_inverse(1, self._rho, energy, 1, ctypes.byref(result)):
                raise ValueError(self._error())
            return result.value
        return self._vector(inverse, energy_J_m3)

    def assess_energy(self, *, current_mass_kg, total_energy_J):
        """Signed support margins in joules; no mutation or endpoint clamping."""
        mass = _positive(current_mass_kg, "current mass")
        energy = np.asarray(total_energy_J, dtype=float)
        if energy.ndim != 0 or not np.isfinite(energy):
            raise ValueError("total energy must be a finite scalar")
        out = _EnergySupportV2()
        if self._lib.ss304v2_energy_support(1, self._rho, mass, float(energy), ctypes.byref(out)):
            raise ValueError(self._error())
        return {name: getattr(out, name) for name, _ in out._fields_}
