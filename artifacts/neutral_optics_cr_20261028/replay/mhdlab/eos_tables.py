"""Bounded SI pressure/internal-energy input for explicitly supplied EOS tables.

This is an offline table handoff, not an EOS generator or material default.
Interpolation preserves the supplied fields; it cannot establish thermodynamic
consistency, phase equilibrium, or the validity of an upstream extrapolation.
"""
from __future__ import annotations

from collections.abc import Mapping
import csv
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import io
import json
import math
from numbers import Real
from pathlib import Path
from types import MappingProxyType

import numpy as np

SCHEMA = "mhdlab-bounded-eos-table-v1"
COLUMNS = ("density_kg_m3", "temperature_k", "pressure_pa", "internal_energy_j_kg")
UNITS = {"density": "kg/m3", "temperature": "K", "pressure": "Pa", "internal_energy": "J/kg"}


def _scalar(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real scalar")
    try:
        value=float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} must be a finite real scalar") from error
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real scalar")
    return value


def _array(values, name, rank):
    if np.ma.isMaskedArray(values):
        raise ValueError(f"{name} cannot be masked")
    raw = np.asarray(values, dtype=object)
    if raw.ndim != rank or not raw.size:
        raise ValueError(f"{name} must be nonempty rank-{rank} data")
    result = np.asarray([_scalar(v, name) for v in raw.flat]).reshape(raw.shape)
    return np.frombuffer(result.tobytes(), dtype=float).reshape(raw.shape)


def _freeze(value):
    if isinstance(value, Mapping):
        return MappingProxyType({str(k): _freeze(v) for k, v in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(v) for v in value)
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    if isinstance(value, float) and math.isfinite(value):
        return value
    raise ValueError("EOS metadata must contain finite JSON-compatible data")


def _evidence(value, name):
    if (not isinstance(value, str) or not value.strip()
            or value.strip().lower().startswith("replace_")
            or value.strip().lower() in {"unknown", "pending", "none", "tbd", "unavailable"}):
        raise ValueError(f"{name} requires explicit source/selection evidence")


def _fraction(value):
    return Fraction.from_float(float(value))


def _bracket(axis, value, name):
    value = _scalar(value, name)
    if not axis[0] <= value <= axis[-1]:
        raise ValueError(f"{name} is outside the EOS table; extrapolation is disabled")
    i = min(int(np.searchsorted(axis, value, side="right"))-1, len(axis)-2)
    weight = (_fraction(value)-_fraction(axis[i]))/(_fraction(axis[i+1])-_fraction(axis[i]))
    return value, i, weight


@dataclass(frozen=True, slots=True)
class EOSValues:
    density_kg_m3: float
    temperature_k: float
    pressure_pa: float
    internal_energy_j_kg: float


@dataclass(frozen=True, slots=True, eq=False)
class BoundedEOSTable:
    density_kg_m3: np.ndarray
    temperature_k: np.ndarray
    pressure_pa: np.ndarray
    internal_energy_j_kg: np.ndarray
    metadata: Mapping

    def __post_init__(self):
        for name, rank in (("density_kg_m3",1), ("temperature_k",1), ("pressure_pa",2), ("internal_energy_j_kg",2)):
            object.__setattr__(self, name, _array(getattr(self,name),name,rank))
        for axis in (self.density_kg_m3, self.temperature_k):
            if len(axis)<2 or np.any(axis<=0) or np.any(np.diff(axis)<=0) or not np.isfinite(axis[-1]-axis[0]):
                raise ValueError("EOS density/temperature axes require finite positive increasing support")
        shape=(len(self.density_kg_m3),len(self.temperature_k))
        if self.pressure_pa.shape!=shape or self.internal_energy_j_kg.shape!=shape:
            raise ValueError("EOS field shape must be (density,temperature)")
        if np.any(self.internal_energy_j_kg[:,1:]<=self.internal_energy_j_kg[:,:-1]):
            raise ValueError("internal energy must strictly increase with temperature for a bounded inverse")
        m=self.metadata
        if not isinstance(m,Mapping) or m.get("schema")!=SCHEMA or m.get("units")!=UNITS:
            raise ValueError("EOS metadata must preserve schema and explicit SI units")
        if m.get("parameter_status") not in {"source_table","manufactured"}:
            raise ValueError("EOS table status must identify source_table or manufactured")
        if m.get("phase_policy") not in {"single_phase","source_equilibrium_multiphase"}:
            raise ValueError("EOS phase policy must be explicit")
        for name in ("material_grade","table_version","source_reference","energy_reference","phase_reference","thermodynamic_consistency_reference"):
            _evidence(m.get(name),name)
        if m.get("interpolation")!="bilinear_in_density_and_temperature" or m.get("extrapolation")!="error":
            raise ValueError("EOS handoff requires its declared interpolation and no extrapolation")
        digest=m.get("csv_sha256")
        if not isinstance(digest,str) or len(digest)!=64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("EOS metadata requires the source CSV SHA256")
        object.__setattr__(self,"metadata",_freeze(m))

    def _interpolate(self, field, ir, wr, it, wt):
        exact=sum((a*b*_fraction(field[r,t])
                   for r,a in ((ir,1-wr),(ir+1,wr))
                   for t,b in ((it,1-wt),(it+1,wt))),Fraction())
        try:
            value=float(exact)
        except OverflowError as error:
            raise ValueError("interpolated EOS value is not representable") from error
        if not math.isfinite(value):
            raise ValueError("interpolated EOS value is not finite")
        return value

    def evaluate(self,density_kg_m3,temperature_k):
        density,ir,wr=_bracket(self.density_kg_m3,density_kg_m3,"density_kg_m3")
        temperature,it,wt=_bracket(self.temperature_k,temperature_k,"temperature_k")
        return EOSValues(density,temperature,self._interpolate(self.pressure_pa,ir,wr,it,wt),
                         self._interpolate(self.internal_energy_j_kg,ir,wr,it,wt))

    def temperature_from_internal_energy_k(self,density_kg_m3,internal_energy_j_kg):
        _,ir,wr=_bracket(self.density_kg_m3,density_kg_m3,"density_kg_m3")
        energy=_fraction(_scalar(internal_energy_j_kg,"internal_energy_j_kg"))
        curve=[(1-wr)*_fraction(self.internal_energy_j_kg[ir,i])+wr*_fraction(self.internal_energy_j_kg[ir+1,i])
               for i in range(len(self.temperature_k))]
        # The forward API rounds a node once. Include those exact rounded limits.
        if not float(curve[0])<=float(energy)<=float(curve[-1]):
            raise ValueError("internal_energy_j_kg is outside the isochoric EOS support")
        if energy<=curve[0]: return float(self.temperature_k[0])
        if energy>=curve[-1]: return float(self.temperature_k[-1])
        i=next(j for j in range(len(curve)-1) if curve[j]<=energy<=curve[j+1])
        fraction=(energy-curve[i])/(curve[i+1]-curve[i])
        return float(_fraction(self.temperature_k[i])+fraction*(_fraction(self.temperature_k[i+1])-_fraction(self.temperature_k[i])))


def load_eos_table(csv_path,metadata_path,*,expected_material_grade):
    """Load a complete hash-bound SI rectangular table; install no model default."""
    raw=Path(csv_path).read_bytes()
    metadata=json.loads(Path(metadata_path).read_bytes())
    if metadata.get("csv_sha256")!=hashlib.sha256(raw).hexdigest():
        raise ValueError("EOS CSV SHA256 mismatch")
    if metadata.get("material_grade")!=expected_material_grade:
        raise ValueError("EOS material grade does not match the caller's explicit selection")
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8"),newline=""))
    if tuple(reader.fieldnames or ())!=COLUMNS:
        raise ValueError(f"EOS columns must be exactly {COLUMNS}")
    rows=[]
    for row in reader:
        if None in row or any(v is None for v in row.values()):
            raise ValueError("malformed EOS row")
        rows.append(tuple(float(row[k]) for k in COLUMNS))
    if not rows or not np.isfinite(rows).all():
        raise ValueError("EOS table must contain finite rows")
    density=np.unique([r[0] for r in rows]);temperature=np.unique([r[1] for r in rows])
    p=np.full((len(density),len(temperature)),np.nan);u=p.copy()
    for rho,t,pressure,energy in rows:
        ir,it=int(np.searchsorted(density,rho)),int(np.searchsorted(temperature,t))
        if not np.isnan(p[ir,it]): raise ValueError("duplicate EOS density/temperature row")
        p[ir,it],u[ir,it]=pressure,energy
    if np.isnan(p).any(): raise ValueError("EOS rows must cover a complete rectangular grid")
    return BoundedEOSTable(density,temperature,p,u,metadata)
