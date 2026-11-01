// chem_kinematics.h -- relativistic collision kinematics of chemistry v3.
//
// Units: energies, masses (rest energies) and momenta in eV (c = 1). Carrier
// states enter and leave as proper velocities u = gamma v in m/s (the native
// arena convention); p = m u / c. Every function uses only +, -, *, /, sqrt
// and det_math.h, in a fixed order, so host and NVRTC builds agree bit for bit.
//
// Frames. A collision is evaluated in the rest frame of the sampled target
// (target frame, "TF"); the pair's centre-of-momentum frame ("CM") moves in
// TF with proper velocity p'/sqrt(s). Boost of a 4-momentum (E, p) from a
// frame moving with proper velocity w (gamma_w = sqrt(1 + w^2)) to its parent:
//   p_parent = p + w (w.p)/(gamma_w + 1) + w E,   E_parent = gamma_w E + w.p;
// and from the parent into the frame:
//   p_frame  = p + w (w.p)/(gamma_w + 1) - w E.
// Kinetic energies are always formed from momenta, K = p^2/(sqrt(p^2+m^2)+m),
// never as differences of total energies. The CM kinetic energy of a
// projectile with TF kinetic energy T on a target of mass M is
//   K_cm = 2 M T / (sqrt((m+M)^2 + 2 M T) + m + M)
// and a two-body final state (masses m3, m4) with available CM kinetic energy
// K has momentum
//   p*^2 = K (K + 2 m3)(K + 2 m4)(K + 2 m3 + 2 m4) / (4 (K + m3 + m4)^2).
// Reaction energies enter as reference-energy differences Q (eV), never as
// differences of rest energies, so no cancellation occurs.
#pragma once
#include "det_math.h"

namespace pfc {

#define PFC_C_M_S 299792458.0

struct V3 {
    double x, y, z;
};
PFC_FN V3 v3(double x, double y, double z) {
    V3 r;
    r.x = x;
    r.y = y;
    r.z = z;
    return r;
}
PFC_FN V3 add(V3 a, V3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
PFC_FN V3 sub(V3 a, V3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
PFC_FN V3 scale(V3 a, double s) { return v3(a.x * s, a.y * s, a.z * s); }
PFC_FN double dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
PFC_FN double norm2(V3 a) { return a.x * a.x + a.y * a.y + a.z * a.z; }
PFC_FN V3 cross(V3 a, V3 b) { return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }

// Canonical kinetic energy (eV) of a carrier of rest energy m (eV) stored with
// proper velocity u (m/s). This is the ledger definition of carrier energy.
PFC_FN double ke_from_u(double m, V3 u) {
    const double c = PFC_C_M_S;
    const double w2 = (u.x / c) * (u.x / c) + (u.y / c) * (u.y / c) + (u.z / c) * (u.z / c);
    return m * gamma_minus_one(w2);
}
PFC_FN V3 p_from_u(double m, V3 u) {
    const double c = PFC_C_M_S;
    return v3(m * (u.x / c), m * (u.y / c), m * (u.z / c));
}
PFC_FN V3 u_from_p(double m, V3 p) {
    const double c = PFC_C_M_S;
    return v3((p.x / m) * c, (p.y / m) * c, (p.z / m) * c);
}
// Dimensionless proper velocity of a frame from m/s.
PFC_FN V3 w_from_u(V3 u) {
    const double c = PFC_C_M_S;
    return v3(u.x / c, u.y / c, u.z / c);
}
PFC_FN double energy_of(double m, V3 p) { return m + kinetic_from_p2(m, norm2(p)); }

// Boost from a frame with proper velocity w (dimensionless) to its parent.
PFC_FN V3 boost_out(V3 w, double gw, V3 p, double E) {
    const double a = dot(w, p) / (gw + 1.0) + E;
    return add(p, scale(w, a));
}
// Boost from the parent into the frame with proper velocity w.
PFC_FN V3 boost_in(V3 w, double gw, V3 p, double E) {
    const double a = dot(w, p) / (gw + 1.0) - E;
    return add(p, scale(w, a));
}

// Unit vector, or (0,0,1) for the zero vector.
PFC_FN V3 unit_or_z(V3 a) {
    const double n2 = norm2(a);
    if (!(n2 > 0.0)) return v3(0.0, 0.0, 1.0);
    const double n = pfc_sqrt(n2);
    return v3(a.x / n, a.y / n, a.z / n);
}
// Isotropic unit vector from two open uniforms.
PFC_FN V3 iso_dir(double u1, double u2) {
    const double mu = 2.0 * u1 - 1.0;
    const double s = pfc_sqrt((1.0 - mu) * (1.0 + mu));
    double c, sn;
    sincos_2pi(u2, &c, &sn);
    return v3(s * c, s * sn, mu);
}
// Rotate unit vector d by polar angle (cos mu) and azimuth (u2) about itself.
PFC_FN V3 deflect(V3 d, double mu, double u2) {
    const double st = pfc_sqrt((1.0 - mu) * (1.0 + mu));
    double cp, sp;
    sincos_2pi(u2, &cp, &sp);
    // Orthonormal basis (e1, e2, d): e1 from the coordinate axis least aligned with d.
    V3 a;
    const double ax = d.x < 0 ? -d.x : d.x, ay = d.y < 0 ? -d.y : d.y, az = d.z < 0 ? -d.z : d.z;
    if (ax <= ay && ax <= az) a = v3(1.0, 0.0, 0.0);
    else if (ay <= az) a = v3(0.0, 1.0, 0.0);
    else a = v3(0.0, 0.0, 1.0);
    const V3 e1 = unit_or_z(cross(d, a));
    const V3 e2 = cross(d, e1);
    return add(add(scale(d, mu), scale(e1, st * cp)), scale(e2, st * sp));
}
// Screened (Okhrimovskyy-type) deflection law f(mu) = (1-xi^2)/(2(1-xi mu)^2)
// on [-1,1]; <mu> = 1/xi - (1-xi^2)/(2 xi^2) ln((1+xi)/(1-xi)); xi = 0 is
// isotropic. Inverse CDF: mu = (2u - 1 + xi)/(1 + xi(2u - 1)).
PFC_FN double screened_mu(double xi, double u) {
    const double a = 2.0 * u - 1.0;
    return (a + xi) / (1.0 + xi * a);
}

// CM kinetic energy of a projectile (mass m, TF kinetic energy T) on a target
// at rest (mass M).
PFC_FN double kcm_from_T(double m, double M, double T) {
    const double S = m + M;
    const double x = 2.0 * M * T;
    return x / (pfc_sqrt(S * S + x) + S);
}
// Two-body momentum for available CM kinetic energy K >= 0 and masses m3, m4.
PFC_FN double pstar(double K, double m3, double m4) {
    if (!(K > 0.0)) return 0.0;
    const double num = K * (K + 2.0 * m3) * (K + 2.0 * m4) * (K + 2.0 * m3 + 2.0 * m4);
    const double den = 2.0 * (K + m3 + m4);
    return pfc_sqrt(num) / den;
}

// One two-body final state in a frame F (given in TF by proper velocity w,
// gamma gw): particle masses m3 (along +d) and m4 (along -d) with momentum ps.
// Results are TF momenta.
struct TwoBody {
    V3 p3, p4;
};
PFC_FN TwoBody two_body_tf(V3 w, double gw, double m3, double m4, double ps, V3 d) {
    const V3 a = scale(d, ps), b = scale(d, -ps);
    const double E3 = energy_of(m3, a), E4 = energy_of(m4, b);
    TwoBody r;
    r.p3 = boost_out(w, gw, a, E3);
    r.p4 = boost_out(w, gw, b, E4);
    return r;
}

}  // namespace pfc
