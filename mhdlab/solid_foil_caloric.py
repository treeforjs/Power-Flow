"""Read-only binding to the separately owned, sealed nominal thermal C ABI.

This adapter does not add or refit cp/k/phase properties. Its selected model is
the owner's fixed-reference-volume Cp-as-Cv approximation, not a physical EOS.
Density is supplied from the mechanics reference geometry, never default 8000.
"""
from __future__ import annotations

import ctypes
import hashlib
import json
from pathlib import Path

import numpy as np


class _Values(ctypes.Structure):
    _fields_ = [(name, ctypes.c_double) for name in (
        "temperature_k", "h_j_kg", "u_model_j_kg", "cp_sensible_j_kg_k",
        "apparent_heat_capacity_j_kg_k", "k_w_m_k", "resistivity_ohm_m",
        "sigma_s_m", "liquid_fraction", "energy_j_m3_reference", "capacity_j_m3_reference_k"
    )] + [("phase", ctypes.c_int), ("synthesized_high_temperature", ctypes.c_int)]


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class SealedNominalCaloric:
    """Scalar C ABI callbacks, with optional vectorization and strict source seal."""

    def __init__(self, source_seal_path, *, expected_source_seal_sha256,
                 reference_density_kg_m3, qualification_mode):
        if qualification_mode != "evaluated_solid_1500":
            raise ValueError("mechanics adapter requires explicit evaluated_solid_1500 mode")
        if not np.isfinite(reference_density_kg_m3) or reference_density_kg_m3 <= 0:
            raise ValueError("positive reference mass density required")
        seal_path=Path(source_seal_path)
        if _sha(seal_path) != expected_source_seal_sha256:
            raise ValueError("thermal source seal SHA256 mismatch")
        seal=json.loads(seal_path.read_text())
        root=Path(seal["source_root"])
        for name in ("include/ss304_thermal.hpp", "ss304_c_api.cpp", "data/scenario.json"):
            if _sha(root/name) != seal["files"][name]:
                raise ValueError("thermal implementation/scenario differs from sealed source")
        library=Path(seal["compiled_library"]["path"])
        if _sha(library) != seal["compiled_library"]["sha256"]:
            raise ValueError("thermal compiled library differs from sealed source")
        self.rho=float(reference_density_kg_m3)
        self._lib=ctypes.CDLL(str(library))
        self._lib.ss304_evaluate.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.POINTER(_Values)]
        self._lib.ss304_evaluate.restype=ctypes.c_int
        self._lib.ss304_inverse.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.POINTER(ctypes.c_double)]
        self._lib.ss304_inverse.restype=ctypes.c_int
        self._lib.ss304_error.restype=ctypes.c_char_p
        self.receipt={"source_seal_sha256":expected_source_seal_sha256,
            "source_root":str(root),"library_sha256":seal["compiled_library"]["sha256"],
            "scenario_sha256":seal["scenario_sha256"],"qualification_mode":qualification_mode,
            "reference_density_kg_m3":self.rho,
            "energy_basis":"rho_ref*u_model(T), fixed-reference-volume Cp-as-Cv approximation; zero at300 K",
            "physical_EOS":False,"mechanics_range_extended":False}

    def _evaluate(self, t):
        values=_Values()
        if self._lib.ss304_evaluate(0,self.rho,float(t),ctypes.byref(values)):
            raise ValueError(self._lib.ss304_error().decode())
        if values.phase != 0:
            raise ValueError("solid mechanics cannot consume a nonsolid caloric state")
        return values

    def _vector(self, function, value):
        array=np.asarray(value,dtype=float)
        if array.size == 0 or not np.isfinite(array).all():
            raise ValueError("finite nonempty caloric argument required")
        # Repeated native temperatures share one scalar lookup, preserving exact ABI.
        unique,index=np.unique(array,return_inverse=True)
        result=np.array([function(float(x)) for x in unique])[index].reshape(array.shape)
        return float(result) if result.ndim == 0 else result

    def energy_J_m3(self, temperature_K):
        return self._vector(lambda t:self._evaluate(t).energy_j_m3_reference,temperature_K)

    def capacity_J_m3_K(self, temperature_K):
        return self._vector(lambda t:self._evaluate(t).capacity_j_m3_reference_k,temperature_K)

    def temperature_K(self, energy_J_m3):
        def inverse(e):
            result=ctypes.c_double()
            if self._lib.ss304_inverse(0,self.rho,e,1,ctypes.byref(result)):
                raise ValueError(self._lib.ss304_error().decode())
            return result.value
        return self._vector(inverse,energy_J_m3)
