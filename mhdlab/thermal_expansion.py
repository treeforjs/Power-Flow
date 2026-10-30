"""Bounded solid SS304 eigenstrain and explicitly reduced axial mechanics.

All quantities are SI. These are uniaxial, small-strain elastic benchmarks,
not a fixture model for the three-dimensional electrodes. Magnetic forces,
plasticity, phase changes and thermal feedback are deliberately separate.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded

MAX_BENCHMARK_STRAIN = 0.02  # implementation guard, not a material yield criterion


def _small_strain(strain):
    if np.max(np.abs(_finite(strain, "computed strain"))) > MAX_BENCHMARK_STRAIN:
        raise ValueError("2 percent small-strain benchmark guard exceeded; not a yield criterion")


def _finite(value, name):
    a = np.asarray(value, dtype=float)
    if a.size == 0 or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be nonempty and finite")
    return a


@dataclass(frozen=True)
class SolidExpansionMaterial:
    """Cumulative length law s(T) = L(T)/L(table_reference) - 1.

    Linear interpolation of *cumulative* expansion integrates the tabulated
    mean expansion coefficients without mistaking them for tangent alpha.
    Between knots the true tangent alpha is s'(T)/(1+s(T)); its integral is
    log((1+s(T))/(1+s(Tref))). No extrapolation is permitted.
    """
    temperature_K: np.ndarray
    relative_length_change: np.ndarray
    modulus_temperature_K: np.ndarray
    young_modulus_Pa: np.ndarray
    density_reference_kg_m3: float
    density_reference_temperature_K: float
    provenance: str

    def __post_init__(self):
        for name in ("temperature_K", "relative_length_change", "modulus_temperature_K", "young_modulus_Pa"):
            a = _finite(getattr(self, name), name).copy()
            if a.ndim != 1 or len(a) < 2:
                raise ValueError(f"{name} must be a 1-D table with at least two values")
            a = np.frombuffer(a.tobytes(), dtype=float)
            object.__setattr__(self, name, a)
        for t, v in ((self.temperature_K, self.relative_length_change),
                     (self.modulus_temperature_K, self.young_modulus_Pa)):
            if t.shape != v.shape or np.any(np.diff(t) <= 0) or np.any(t <= 0):
                raise ValueError("table temperatures must increase and match values")
        if np.any(self.relative_length_change <= -1) or np.any(self.young_modulus_Pa <= 0):
            raise ValueError("positive stretches and elastic moduli required")
        if not np.isfinite(self.density_reference_kg_m3) or self.density_reference_kg_m3 <= 0:
            raise ValueError("positive finite reference density required")
        self.check_temperature(self.density_reference_temperature_K)
        if not self.provenance:
            raise ValueError("material provenance required")

    @property
    def admissible_temperature_K(self):
        return (max(self.temperature_K[0], self.modulus_temperature_K[0]),
                min(self.temperature_K[-1], self.modulus_temperature_K[-1]))

    def check_temperature(self, temperature_K):
        t = _finite(temperature_K, "temperature")
        lo, hi = self.admissible_temperature_K
        if np.any((t < lo) | (t > hi)):
            raise ValueError(f"solid elastic material supports only [{lo}, {hi}] K; no clipping or extrapolation")
        return t

    def strain(self, temperature_K, *, stress_free_temperature_K):
        t = self.check_temperature(temperature_K)
        ref = self.check_temperature(stress_free_temperature_K)
        if ref.ndim != 0:
            raise ValueError("one explicit scalar stress-free temperature required")
        s = np.interp(t, self.temperature_K, self.relative_length_change)
        sr = np.interp(ref, self.temperature_K, self.relative_length_change)
        return (s - sr) / (1.0 + sr)

    def integrated_alpha(self, temperature_K, *, stress_free_temperature_K):
        return np.log1p(self.strain(temperature_K, stress_free_temperature_K=stress_free_temperature_K))

    def modulus(self, temperature_K):
        t = self.check_temperature(temperature_K)
        return np.interp(t, self.modulus_temperature_K, self.young_modulus_Pa)

    def density_at_stress_free_temperature(self, stress_free_temperature_K):
        stretch = 1 + self.strain(stress_free_temperature_K,
                                 stress_free_temperature_K=self.density_reference_temperature_K)
        return float(self.density_reference_kg_m3 / stretch**3)


def load_ss304_material(path=None):
    path = Path(path) if path is not None else Path(__file__).resolve().parents[1] / "data/materials/ss304-v1.json"
    raw = path.read_bytes()
    data = json.loads(raw)
    if (data.get("record_id") != "ss304-solid-expansion-v1" or data.get("schema_version") != 1
            or data.get("material", {}).get("UNS") != "S30400"
            or data.get("material", {}).get("EN") != "1.4301"):
        raise ValueError("expected versioned Type 304 / S30400 / 1.4301 material record")
    expansion, modulus, density = (data[k] for k in ("thermal_expansion", "young_modulus", "density_reference"))
    material = SolidExpansionMaterial(
        np.array(expansion["temperature_K"]), np.array(expansion["relative_length_change"]),
        np.array(modulus["temperature_K"]), np.array(modulus["value_Pa"]),
        density["value_kg_m3"], density["temperature_K"],
        f"{path.name};sha256={hashlib.sha256(raw).hexdigest()}")
    if tuple(data["valid_temperature_K"]) != material.admissible_temperature_K:
        raise ValueError("material admissible domain differs from its property tables")
    return material


def isotropic_elastic_stress(strain_tensor, temperature_K, *, material,
                             stress_free_temperature_K, poisson_ratio):
    """3-D Hooke law; caller must supply/document nu (not a hidden SS304 datum).

    Free uniform expansion has strain=eps_th*I and zero stress. Fully fixed
    strain=0 gives -E*eps_th/(1-2nu)*I, unlike a uniaxial restrained bar.
    This returns an elastic trial stress, not a yield/plasticity prediction.
    """
    strain = _finite(strain_tensor, "strain tensor")
    if strain.shape != (3, 3) or not np.allclose(strain, strain.T, atol=1e-15, rtol=0):
        raise ValueError("symmetric 3x3 strain tensor required")
    _small_strain(np.linalg.eigvalsh(strain))
    nu = float(poisson_ratio)
    if not np.isfinite(nu) or not -1 < nu < .5:
        raise ValueError("Poisson ratio must lie in (-1, .5)")
    e = float(material.modulus(temperature_K))
    elastic = strain - float(material.strain(temperature_K, stress_free_temperature_K=stress_free_temperature_K)) * np.eye(3)
    return e/(1+nu)*elastic + e*nu/((1+nu)*(1-2*nu))*np.trace(elastic)*np.eye(3)


def _bar_inputs(nodes_m, cell_temperature_K, material, stress_free_temperature_K, area_m2):
    x = _finite(nodes_m, "nodes")
    temp = material.check_temperature(cell_temperature_K)
    if x.ndim != 1 or len(x) < 2 or np.any(np.diff(x) <= 0):
        raise ValueError("strictly increasing 1-D nodes required")
    if temp.shape[-1:] != (len(x)-1,):
        raise ValueError("one temperature per reference element required")
    if not np.isfinite(area_m2) or area_m2 <= 0:
        raise ValueError("finite positive reference area required")
    h = _finite(np.diff(x), "element lengths")
    eps = material.strain(temp, stress_free_temperature_K=stress_free_temperature_K)
    e = material.modulus(temp)
    rho = material.density_at_stress_free_temperature(stress_free_temperature_K)
    return x, temp, h, eps, e, rho


def static_bar(nodes_m, cell_temperature_K, *, material, stress_free_temperature_K,
               boundary, area_m2=1.0):
    """Exact equilibrium for piecewise-constant 1-D E and eigenstrain.

    free: both ends traction free, mass-center translation set to zero only
    as a rigid-motion gauge. fixed_left: u(left)=0, right traction free.
    fixed_both: u(left)=u(right)=0. These are *benchmark* boundary cases.
    """
    x, temp, h, eps, e, rho = _bar_inputs(nodes_m, cell_temperature_K, material, stress_free_temperature_K, area_m2)
    if temp.ndim != 1 or boundary not in ("free", "fixed_left", "fixed_both"):
        raise ValueError("1-D temperatures and an explicit benchmark boundary required")
    force = -np.sum(eps*h)/np.sum(h/(e*area_m2)) if boundary == "fixed_both" else 0.0
    strain = eps + force/(e*area_m2)
    _small_strain(strain)
    u = np.r_[0.0, np.cumsum(h*strain)]
    if boundary == "free":
        u -= np.sum(h*(u[:-1]+u[1:])/2)/np.sum(h)
    if np.any(1+strain <= 0):
        raise ValueError("nonpositive rod Jacobian")
    return dict(displacement_m=u, strain=strain, thermal_strain=eps,
                stress_Pa=e*(strain-eps), axial_force_N=float(force),
                reaction_N=np.array([-force, force]),
                density_kg_m3=rho/(1+strain), jacobian=1+strain,
                reference_mass_kg=float(rho*area_m2*np.sum(h)),
                current_mass_kg=float(np.sum(rho/(1+strain)*area_m2*h*(1+strain))),
                boundary=boundary, scope="uniaxial_constant_area_small_strain_benchmark")


def _operators(h, e, eps, area):
    stiffness = area*e/h
    diagonal = np.r_[stiffness, 0.] + np.r_[0., stiffness]
    thermal_force = np.r_[-area*e*eps, 0.] + np.r_[0., area*e*eps]
    return diagonal, -stiffness, thermal_force


def _apply(diagonal, off, u):
    result = diagonal*u
    result[:-1] += off*u[1:]
    result[1:] += off*u[:-1]
    return result


def transient_bar(nodes_m, times_s, cell_temperature_K, *, material,
                  stress_free_temperature_K, boundary, area_m2=1.0,
                  initial_displacement_m=None, initial_velocity_m_s=None):
    """Lagrangian P1 bar, lumped reference mass, Newmark average acceleration.

    Temperatures are supplied at every time; variable E and integrated thermal
    strain are evaluated at each endpoint. No thermal-history extrapolation.
    Initial conditions default to zero displacement/velocity, NOT equilibrium
    at a hot first frame. The initial acceleration follows the initial forces.
    No damping or applied magnetic pressure is silently introduced.
    """
    x, temp, h, eps, e, rho = _bar_inputs(nodes_m, cell_temperature_K, material, stress_free_temperature_K, area_m2)
    times = _finite(times_s, "times")
    if times.ndim != 1 or len(times) < 2 or np.any(np.diff(times) <= 0) or temp.shape != (len(times), len(h)):
        raise ValueError("increasing times and (time, element) temperatures required")
    if boundary not in ("free", "fixed_left", "fixed_both"):
        raise ValueError("explicit benchmark boundary required")
    n = len(x)
    constrained = [] if boundary == "free" else ([0] if boundary == "fixed_left" else [0, n-1])
    active = np.arange(n)[slice(0 if boundary == "free" else 1, n-1 if boundary == "fixed_both" else n)]
    mass = rho*area_m2*(np.r_[h, 0.]+np.r_[0., h])/2
    u = np.zeros((len(times), n)); v = np.zeros_like(u); a = np.zeros_like(u)
    for target, value, name in ((u, initial_displacement_m, "initial displacement"), (v, initial_velocity_m_s, "initial velocity")):
        if value is not None:
            value = _finite(value, name)
            if value.shape != (n,) or np.any(value[constrained] != 0):
                raise ValueError("initial state must match nodes and homogeneous constraints")
            target[0] = value
    d, off, f = _operators(h, e[0], eps[0], area_m2)
    a[0, active] = ((f-_apply(d, off, u[0]))/mass)[active]
    reactions = np.zeros_like(u)
    reactions[0, constrained] = (_apply(d, off, u[0])-f)[constrained]
    for k, dt in enumerate(np.diff(times), 1):
        dt2 = float(dt)*float(dt)
        if not np.isfinite(dt2) or dt2 <= 0:
            raise ValueError("time-step square is outside finite floating-point support")
        d, off, f = _operators(h, e[k], eps[k], area_m2)
        predictor = u[k-1]+dt*v[k-1]+.25*dt2*a[k-1]
        rhs = f+4*mass*predictor/dt2
        band = np.zeros((3, len(active)))
        band[1] = (d+4*mass/dt2)[active]
        if len(active) > 1:
            band[0, 1:] = off[active[:-1]]
            band[2, :-1] = off[active[:-1]]
        if len(active):
            u[k, active] = solve_banded((1, 1), band, rhs[active], check_finite=True)
        a[k] = 4*(u[k]-predictor)/dt2
        v[k] = v[k-1]+.5*dt*(a[k-1]+a[k])
        a[k, constrained] = 0.; v[k, constrained] = 0.
        reactions[k, constrained] = (_apply(d, off, u[k])-f)[constrained]
    strain = np.diff(u, axis=1)/h
    _small_strain(strain)
    jacobian = 1+strain
    if np.any(jacobian <= 0):
        raise ValueError("nonpositive rod Jacobian")
    kinetic = .5*np.sum(mass*v*v, axis=1)
    elastic = .5*np.sum(e*(strain-eps)**2*h*area_m2, axis=1)
    for value in (u, v, a, reactions, kinetic, elastic):
        _finite(value, "computed mechanical state")
    return dict(times_s=times, displacement_m=u, velocity_m_s=v, acceleration_m_s2=a,
                strain=strain, thermal_strain=eps, stress_Pa=e*(strain-eps),
                reaction_N=reactions, nodal_reference_mass_kg=mass,
                density_kg_m3=rho/jacobian, jacobian=jacobian,
                reference_mass_kg=float(np.sum(mass)),
                current_mass_kg=np.sum(rho/jacobian*jacobian*h*area_m2, axis=1),
                kinetic_energy_J=kinetic, elastic_energy_J=elastic,
                boundary=boundary, scope="uniaxial_constant_area_small_strain_transient_benchmark")
