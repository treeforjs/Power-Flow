"""Reference-volume thermoelastic energy identities and conditional heat ports.

The algorithmic port conserves kinetic + elastic free + reservoir energy. The
conditional Helmholtz port instead uses elastic INTERNAL energy Phi-T*Phi_T.
Neither certifies the temperature-dependent material inputs as a physical EOS.
Engineering Voigt order is xx, yy, zz, gamma_xy, gamma_yz, gamma_xz; its conjugate
stress order is xx, yy, zz, tau_xy, tau_yz, tau_xz. Dot products give work directly.
All densities are per FIXED reference volume; multiply by that volume once.
"""
from __future__ import annotations

import numpy as np


STRAIN_BASIS = "engineering_voigt_xx_yy_zz_xy_yz_xz"


def _scalar(value, name, *, positive=False):
    a = np.asarray(value, dtype=float)
    if a.ndim != 0 or not np.isfinite(a) or (positive and a <= 0):
        raise ValueError(f"{name} must be a finite {'positive ' if positive else ''}scalar")
    return float(a)


def _vector(value, name):
    a = np.atleast_1d(np.asarray(value, dtype=float))
    if a.ndim != 1 or a.size == 0 or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be a nonempty finite vector")
    return a


def _matrix(value, n, name, *, stiffness=False):
    a = np.asarray(value, dtype=float)
    if a.ndim == 0 and n == 1:
        a = a.reshape(1, 1)
    if a.shape != (n, n) or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be finite with shape {(n, n)}")
    if not np.allclose(a, a.T, rtol=1e-13, atol=0):
        raise ValueError(f"{name} must be symmetric in energy-conjugate coordinates")
    if stiffness:
        try:
            np.linalg.cholesky(a)
        except np.linalg.LinAlgError as exc:
            raise ValueError("stiffness must be positive definite") from exc
    return a


def _matched_vectors(*values):
    arrays = [_vector(value, f"vector {i}") for i, value in enumerate(values)]
    if any(a.shape != arrays[0].shape for a in arrays[1:]):
        raise ValueError("strain and thermal-strain vectors must have equal shape")
    return arrays


def algorithmic_exchange(strain0, strain1, thermal_strain0, thermal_strain1,
                         stiffness0_Pa, stiffness1_Pa, *, reference_volume_m3):
    """Exact endpoint identity; positive debit removes reservoir energy.

    sigma_bar = (C0 e0 + C1 e1)/2, Wstress = V*sigma_bar.dot(eps1-eps0).
    DeltaPhi-Wstress = -V*sigma_bar.dot(eth1-eth0)+V*e1.dot((C1-C0)e0)/2.
    With midpoint-compatible momentum, DeltaK=-Wstress+Wexternal, so the
    reported debit also equals Delta(K+Phi)-Wexternal. It is NOT physical heat.
    """
    eps0, eps1, th0, th1 = _matched_vectors(strain0, strain1, thermal_strain0, thermal_strain1)
    c0 = _matrix(stiffness0_Pa, eps0.size, "C0", stiffness=True)
    c1 = _matrix(stiffness1_Pa, eps0.size, "C1", stiffness=True)
    volume = _scalar(reference_volume_m3, "reference volume", positive=True)
    e0, e1 = eps0-th0, eps1-th1
    sig0, sig1 = c0@e0, c1@e1
    sigbar = (sig0+sig1)/2
    phi0, phi1 = .5*e0@sig0, .5*e1@sig1
    work = volume*(sigbar@(eps1-eps0))
    delta = volume*(phi1-phi0)
    eigen = -volume*(sigbar@(th1-th0))
    modulus = .5*volume*(e1@((c1-c0)@e0))
    result = {
        "basis": STRAIN_BASIS if eps0.size == 6 else "energy_conjugate_vector",
        "reference_volume_m3": volume,
        "work_stress_J": float(work),
        "delta_phi_J": float(delta),
        "algorithmic_heat_debit_J": float(delta-work),
        "thermal_eigenstrain_term_J": float(eigen),
        "modulus_term_J": float(modulus),
        "closure_residual_J": float(delta-work-eigen-modulus),
        "physical_EOS_ready": False,
        "energy_model": "algorithmic_free_energy_plus_thermal_reservoir",
    }
    if not all(np.isfinite(v) for v in result.values() if isinstance(v, float)):
        raise ValueError("nonfinite energy calculation")
    return result


def thermodynamic_state(*, temperature_K, strain, thermal_strain,
                        dthermal_strain_dT, d2thermal_strain_dT2,
                        stiffness_Pa, dstiffness_dT_Pa_K,
                        d2stiffness_dT2_Pa_K2, base_energy_J_m3,
                        base_capacity_J_m3_K):
    """Evaluate a smooth, explicitly supplied Helmholtz constitutive model.

    psi=psi0(T)+Phi, Phi=.5*(eps-eth):C:(eps-eth), U0=psi0-T*psi0_T.
    C0heat=dU0/dT is the caloric capacity along the stress-free path. A Cp-as-Cv
    proxy can be passed only as a declared conditional model, not physical EOS.
    Returned entropy is elastic excess only; base entropy needs its own datum.
    The caller owns differentiability, material validity, and derivative units.
    """
    t = _scalar(temperature_K, "temperature", positive=True)
    eps, th, th_t, th_tt = _matched_vectors(strain, thermal_strain, dthermal_strain_dT,
                                            d2thermal_strain_dT2)
    c = _matrix(stiffness_Pa, eps.size, "C", stiffness=True)
    ct = _matrix(dstiffness_dT_Pa_K, eps.size, "C_T")
    ctt = _matrix(d2stiffness_dT2_Pa_K2, eps.size, "C_TT")
    u0 = _scalar(base_energy_J_m3, "base energy")
    capacity0 = _scalar(base_capacity_J_m3_K, "base capacity", positive=True)
    e = eps-th
    stress = c@e
    stress_t = ct@e-c@th_t
    phi = .5*e@stress
    phi_t = .5*e@(ct@e)-stress@th_t
    phi_tt = .5*e@(ctt@e)-2*th_t@(ct@e)+th_t@(c@th_t)-stress@th_tt
    uelastic = phi-t*phi_t
    capacity_eps = capacity0-t*phi_tt
    if not np.isfinite(capacity_eps) or capacity_eps <= 0:
        raise ValueError("conditional fixed-strain capacity must remain positive")
    out = {
        "temperature_K": t,
        "stress_Pa": stress,
        "stress_temperature_derivative_Pa_K": stress_t,
        "phi_J_m3": float(phi),
        "phi_temperature_derivative_J_m3_K": float(phi_t),
        "phi_second_temperature_derivative_J_m3_K2": float(phi_tt),
        "elastic_entropy_J_m3_K": float(-phi_t),
        "elastic_internal_energy_J_m3": float(uelastic),
        "base_energy_J_m3": u0,
        "total_internal_energy_J_m3": float(u0+uelastic),
        "base_capacity_J_m3_K": capacity0,
        "fixed_strain_capacity_J_m3_K": float(capacity_eps),
        "physical_EOS_ready": False,
        "basis": STRAIN_BASIS if eps.size == 6 else "energy_conjugate_vector",
    }
    if any(not np.all(np.isfinite(v)) for v in out.values()
           if isinstance(v, (float, np.ndarray))):
        raise ValueError("nonfinite thermodynamic state")
    return out


def isotropic_stiffness(young_modulus_Pa, poisson_ratio):
    """Engineering Voigt stiffness; nu is an explicit caller-owned assumption."""
    young = _scalar(young_modulus_Pa, "Young modulus", positive=True)
    nu = _scalar(poisson_ratio, "Poisson ratio")
    if not -1 < nu < .5:
        raise ValueError("Poisson ratio must lie in (-1,.5)")
    mu = young/(2*(1+nu))
    lam = young*nu/((1+nu)*(1-2*nu))
    c = np.diag([2*mu]*3+[mu]*3)
    c[:3, :3] += lam
    return c


def _open_segment(t, knots):
    if t <= knots[0] or t >= knots[-1] or np.any(t == knots):
        raise ValueError("thermodynamic mode requires temperatures inside open smooth coefficient segments; knots are not differentiable")
    return int(np.searchsorted(knots, t)-1)


def _ss304_segment_values(material, t, reference, nu):
    material.check_temperature(t)
    material.check_temperature(reference)
    i = _open_segment(t, material.temperature_K)
    j = _open_segment(t, material.modulus_temperature_K)
    sref = np.interp(reference, material.temperature_K, material.relative_length_change)
    eth = float(material.strain(t, stress_free_temperature_K=reference))
    eth_t = (np.diff(material.relative_length_change)[i]/np.diff(material.temperature_K)[i])/(1+sref)
    young = float(material.modulus(t))
    young_t = np.diff(material.young_modulus_Pa)[j]/np.diff(material.modulus_temperature_K)[j]
    unit_c = isotropic_stiffness(1.0, nu)
    normal = np.array([1., 1., 1., 0., 0., 0.])
    return (i, j), eth*normal, eth_t*normal, young*unit_c, young_t*unit_c


def ss304_thermodynamic_increment(strain0, strain1, temperature0_K, temperature1_K, *,
                                  material, stress_free_temperature_K, poisson_ratio,
                                  reference_volume_m3, base_energy_J_m3,
                                  base_capacity_J_m3_K, supplied_heat_J=0.0):
    """Single Gauss-point conditional internal-energy increment in joules.

    The two base arguments are scalar callbacks in temperature, evaluated at
    each endpoint. They must use the SAME accepted reference mass/volume as
    mechanics. Calls to the separate thermal owner retain its Cp-as-Cv proxy
    status and must enforce its own domain, currently T>=300 K.

    No crossing or landing on a coefficient knot is allowed in this physical-
    form mode. The original piecewise-linear laws are not a globally smooth
    Helmholtz EOS. Algorithmic_exchange remains available across such knots.

    Positive conditional_heat_debit_J is removed from base thermal energy:
    DeltaU0 = supplied_heat - (DeltaUelastic-Wstress).
    This debit also equals Delta(K+Uelastic)-Wexternal ONLY for momentum updates
    satisfying DeltaK=-Wstress+Wexternal. Caller must check that ledger itself.
    """
    eps0, eps1 = _matched_vectors(strain0, strain1)
    if eps0.size != 6:
        raise ValueError("SS304 isotropic helper requires six engineering Voigt strains")
    t0 = _scalar(temperature0_K, "T0", positive=True)
    t1 = _scalar(temperature1_K, "T1", positive=True)
    ref = _scalar(stress_free_temperature_K, "stress-free temperature", positive=True)
    volume = _scalar(reference_volume_m3, "reference volume", positive=True)
    heat = _scalar(supplied_heat_J, "supplied heat")
    values0 = _ss304_segment_values(material, t0, ref, poisson_ratio)
    values1 = _ss304_segment_values(material, t1, ref, poisson_ratio)
    if values0[0] != values1[0]:
        raise ValueError("thermodynamic increment cannot cross expansion or stiffness coefficient knots")
    states = []
    for eps, t, values in ((eps0, t0, values0), (eps1, t1, values1)):
        _, eth, eth_t, c, ct = values
        states.append(thermodynamic_state(
            temperature_K=t, strain=eps, thermal_strain=eth,
            dthermal_strain_dT=eth_t, d2thermal_strain_dT2=np.zeros(6),
            stiffness_Pa=c, dstiffness_dT_Pa_K=ct, d2stiffness_dT2_Pa_K2=np.zeros((6,6)),
            base_energy_J_m3=base_energy_J_m3(t), base_capacity_J_m3_K=base_capacity_J_m3_K(t)))
    s0, s1 = states
    out = algorithmic_exchange(eps0, eps1, values0[1], values1[1], values0[3], values1[3],
                                reference_volume_m3=volume)
    delta_uel = volume*(s1["elastic_internal_energy_J_m3"]-s0["elastic_internal_energy_J_m3"])
    delta_base = volume*(s1["base_energy_J_m3"]-s0["base_energy_J_m3"])
    debit = delta_uel-out["work_stress_J"]
    out.update({
        "energy_model": "conditional_smooth_segment_Helmholtz_internal_energy",
        "delta_elastic_internal_energy_J": float(delta_uel),
        "conditional_heat_debit_J": float(debit),
        "delta_base_energy_J": float(delta_base),
        "supplied_heat_J": heat,
        "energy_balance_residual_J": float(delta_base+debit-heat),
        "state0": s0,
        "state1": s1,
        "coefficient_segment_indices": values0[0],
        "caloric_qualification": "caller_owned; physical_EOS_ready remains false",
    })
    if any(not np.isfinite(v) for v in out.values() if isinstance(v, float)):
        raise ValueError("nonfinite conditional energy increment")
    return out
