// phase_law.h -- SS304 phase law v1 (298.2-5000 K) and a constant-property benchmark law, evaluated with
// identical arithmetic on the host and in NVRTC kernels. Provenance of every number: DATA-SOURCES.json.
//
// Caloric law (J/kg, energy zero at 300 K solid; fixed-reference-volume convention of the sealed heat-v2 law):
//   298.2 <= T < 1672 K  solid    h = (T-300)(472 + 0.068(T+300) - 2.82e6/(300 T))           Mills 2004 Eq. 10/11
//   1672 <= T < 1727 K   mushy    sensible part with Cp blended linearly to 800, plus L_f (T-1672)/55,
//                                  L_f = 160 J/(K kg) * 1727 K = 276320 J/kg                    Mills Eq. 12/13,
//                                  solidus/liquidus = producer melting range 1399-1454 C (Nippon Steel TP 304)
//   1727 <= T <= 5000 K  liquid   h = h(1727) + 800 (T-1727)                                   Mills Eq. 12
// These are the sealed ss304::Model / ss304v2::Model expressions, written in the same operation order, so h, Cp,
// k, resistivity and liquid fraction are bit-identical to the sealed law wherever the sealed law is defined
// (298.2-1800 K); above 1800 K the liquid branch is the same formula continued.
// Transport: Ho & Chu 1977 Table 2 (PWL, sealed) to 1800 K; above 1800 K k follows the Mills 2004 Table 3(b)
// liquid slope (29 -> 30 W/m/K from 1800 to 1900 K) continued linearly, and the resistivity follows the constant
// Wiedemann-Franz ratio anchored at the Ho & Chu 1800 K row (the anchored ratio is 2.443e-8 W Ohm/K^2).
// Vapor: pure-liquid vapor pressures of Fe, Cr, Ni from NIST-JANAF (gas minus liquid Gibbs energy at the table
// rows; liquid continued with its tabulated constant Cp beyond the table end), ln p interpolated linearly in 1/T
// between rows; Raoult activity a_i = gamma_i x_i (gamma_i = 1 by default; Mills 2004: Fe-Cr and Fe-Ni nearly
// ideal); Hertz-Knudsen-Langmuir flux j_i = (1-beta_r) alpha_e a_i p_i / sqrt(2 pi m_i k T).
#pragma once
#include "det_math.h"

#define PT_MAX_KNOTS 72
#define PT_NEL 3   // Fe, Cr, Ni

// Status codes written by the per-cell arithmetic (0 = ok).
#define PT_OK 0
#define PT_ERR_T_LOW 1          // specific energy below h(T_min)
#define PT_ERR_T_HIGH 2         // specific energy above h(T_max)
#define PT_ERR_NONFINITE 3
#define PT_ERR_FIXED_RANGE 4    // fixed-point conversion out of range / overflow
#define PT_ERR_DEPLETED 5       // a cell lost more than the admitted fraction of its reference mass in one step
#define PT_ERR_NEGATIVE_INV 6   // an element inventory would become negative

struct pt_law {
    int kind;                 // 0: SS304 phase law v1, 1: constant-property benchmark law
    double T_min, T_max, h_min, h_max;
    // constant-property benchmark law (kind 1): h = 0 at c_T0
    double c_T0, c_cps, c_cpl, c_Tm, c_dTm, c_L, c_ks, c_kl, c_rhoe;
    // vapor law (kind 0 only)
    int vapor_enabled;
    int n_knots[PT_NEL];
    double knot_T[PT_NEL][PT_MAX_KNOTS], knot_lnp[PT_NEL][PT_MAX_KNOTS];   // K, ln(p/Pa)
    double atom_mass[PT_NEL], eps0[PT_NEL], gamma[PT_NEL];                 // kg, J per atom, 1
    double alpha_e, beta_r;
};

#define PT_KB 1.380649e-23
#define PT_PI 3.14159265358979323846

// --------------------------------------------------------------------------------------------------------
// SS304 phase law v1 (sealed expressions).
PT_HD double pt_ss_cp_solid(double t) { return 472.0 + 0.136 * t - 2820000.0 / (t * t); }
PT_HD double pt_ss_h_solid(double t) { return (t - 300.0) * (472.0 + 0.068 * (t + 300.0) - 2820000.0 / (t * 300.0)); }
PT_HD double pt_ss_h_mushy_sensible(double dt) {
    const double cs = pt_ss_cp_solid(1672.0);
    return pt_ss_h_solid(1672.0) + cs * dt + 0.5 * (800.0 - cs) * dt * dt / (1727.0 - 1672.0);
}
#define PT_SS_LATENT (160.0 * 1727.0)
PT_HD double pt_ss_h(double t) {
    const double width = 1727.0 - 1672.0;
    if (t < 1672.0) return pt_ss_h_solid(t);
    if (t < 1727.0) { const double dt = t - 1672.0; return pt_ss_h_mushy_sensible(dt) + PT_SS_LATENT * (dt / width); }
    return pt_ss_h_mushy_sensible(width) + PT_SS_LATENT + 800.0 * (t - 1727.0);
}
// Apparent capacity dh/dT (right-sided at the phase boundaries), as the sealed evaluate().
PT_HD double pt_ss_cap(double t) {
    if (t >= 1727.0) return 800.0;
    if (t >= 1672.0) {
        const double f = (t - 1672.0) / (1727.0 - 1672.0);
        const double cp = pt_ss_cp_solid(1672.0) + (800.0 - pt_ss_cp_solid(1672.0)) * f;
        return cp + PT_SS_LATENT / (1727.0 - 1672.0);
    }
    return pt_ss_cp_solid(t);
}
PT_HD double pt_ss_liquid_fraction(double t) {
    if (t >= 1727.0) return 1.0;
    if (t >= 1672.0) return (t - 1672.0) / (1727.0 - 1672.0);
    return 0.0;
}
// Resistivity (Ohm m) and thermal conductivity (W/m/K).
PT_HD void pt_ss_transport(double t, double* rho, double* k) {
    if (t < 300.0) {
        // Ho & Chu Table 2 rows 293 K and 300 K (sealed ss304v2 cold branch).
        const double f = (t - 293) / 7;
        *rho = (71.3 + f * (71.9 - 71.3)) * 1e-8;
        *k = 14.76 + f * (14.89 - 14.76);
        return;
    }
    if (t <= 1800.0) {
        const double ts[18] = {300,350,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1672,1727,1800};
        const double rs[18] = {71.9,76.0,79.8,86.8,93.3,99.2,104.3,108.6,112.5,115.8,118.7,121.4,124.0,126.4,128.7,130.4,150.0,151.7};
        const double ks[18] = {14.89,15.79,16.61,18.28,19.77,21.21,22.59,23.99,25.33,26.58,27.81,29.18,30.34,31.55,32.70,33.53,28.15,28.99};
        int upper = 18;
        for (int i = 0; i < 18; ++i) if (ts[i] > t) { upper = i; break; }
        int i = upper - 1; if (i > 16) i = 16;
        const double f = (t - ts[i]) / (ts[i + 1] - ts[i]);
        *rho = (rs[i] + f * (rs[i + 1] - rs[i])) * 1e-8;
        *k = ks[i] + f * (ks[i + 1] - ks[i]);
        return;
    }
    // Liquid continuation: Mills 2004 Table 3(b) slope (1 W/m/K per 100 K), constant anchored Lorenz ratio.
    const double kk = 28.99 + 0.01 * (t - 1800.0);
    const double lorenz = (28.99 * 151.7e-8) / 1800.0;
    *k = kk;
    *rho = lorenz * t / kk;
}
// Density (kg/m3), diagnostic only (the owner keeps a fixed reference volume): Mills 2004 Table 3 / Eq. 5, 6.
PT_HD double pt_ss_density(double t) {
    if (t < 1672.0) return 7921.0 - 0.5 * (t - 298.2);
    if (t >= 1727.0) return 6908.0 - 0.86 * (t - 1723.0);
    const double f = (t - 1672.0) / (1727.0 - 1672.0);
    const double vs = 1.0 / (7921.0 - 0.5 * (1672.0 - 298.2)), vl = 1.0 / (6908.0 - 0.86 * (1727.0 - 1723.0));
    return 1.0 / ((1.0 - f) * vs + f * vl);
}

// --------------------------------------------------------------------------------------------------------
// Constant-property benchmark law (kind 1).
PT_HD double pt_c_h(const pt_law* L, double t) {
    const double t1 = L->c_Tm - 0.5 * L->c_dTm, t2 = L->c_Tm + 0.5 * L->c_dTm;
    const double h1 = L->c_cps * (t1 - L->c_T0);
    if (t < t1) return L->c_cps * (t - L->c_T0);
    const double cm = 0.5 * (L->c_cps + L->c_cpl);
    if (t < t2) return h1 + (cm + L->c_L / L->c_dTm) * (t - t1);
    const double h2 = h1 + cm * L->c_dTm + L->c_L;
    return h2 + L->c_cpl * (t - t2);
}
PT_HD double pt_c_cap(const pt_law* L, double t) {
    const double t1 = L->c_Tm - 0.5 * L->c_dTm, t2 = L->c_Tm + 0.5 * L->c_dTm;
    if (t < t1) return L->c_cps;
    if (t < t2) return 0.5 * (L->c_cps + L->c_cpl) + L->c_L / L->c_dTm;
    return L->c_cpl;
}
PT_HD double pt_c_fl(const pt_law* L, double t) {
    const double t1 = L->c_Tm - 0.5 * L->c_dTm, t2 = L->c_Tm + 0.5 * L->c_dTm;
    if (t < t1) return 0.0;
    if (t < t2) return (t - t1) / L->c_dTm;
    return 1.0;
}
PT_HD double pt_c_k(const pt_law* L, double t) {
    const double f = pt_c_fl(L, t);
    return L->c_ks + f * (L->c_kl - L->c_ks);
}

// --------------------------------------------------------------------------------------------------------
// Law dispatch.
PT_HD double pt_h(const pt_law* L, double t) { return L->kind == 0 ? pt_ss_h(t) : pt_c_h(L, t); }
PT_HD double pt_cap(const pt_law* L, double t) { return L->kind == 0 ? pt_ss_cap(t) : pt_c_cap(L, t); }
PT_HD double pt_fl(const pt_law* L, double t) { return L->kind == 0 ? pt_ss_liquid_fraction(t) : pt_c_fl(L, t); }
PT_HD double pt_k(const pt_law* L, double t) {
    if (L->kind != 0) return pt_c_k(L, t);
    double r, k; pt_ss_transport(t, &r, &k); return k;
}
PT_HD double pt_rhoe(const pt_law* L, double t) {
    if (L->kind != 0) return L->c_rhoe;
    double r, k; pt_ss_transport(t, &r, &k); return r;
}

// Choose among T and its neighbours the representable temperature whose forward enthalpy is closest to h
// (at most 4 steps each way), inside [T_min, T_max].
PT_HD double pt_polish(const pt_law* L, double h, double t) {
    for (int it = 0; it < 4; ++it) {
        const double e = pt_fabs(pt_h(L, t) - h);
        const double up = pt_next_up(t), dn = pt_next_down(t);
        double best = t, be = e;
        if (up <= L->T_max) { const double eu = pt_fabs(pt_h(L, up) - h); if (eu < be) { best = up; be = eu; } }
        if (dn >= L->T_min) { const double ed = pt_fabs(pt_h(L, dn) - h); if (ed < be) { best = dn; be = ed; } }
        if (best == t) break;
        t = best;
    }
    return t;
}
// Temperature from specific energy (J/kg). *status receives PT_ERR_T_LOW/HIGH/NONFINITE outside the domain.
PT_HD double pt_T_of_h(const pt_law* L, double h, int* status) {
    if (!pt_finite(h)) { *status = PT_ERR_NONFINITE; return 0.0; }
    if (h < L->h_min) { *status = PT_ERR_T_LOW; return L->T_min; }
    if (h > L->h_max) { *status = PT_ERR_T_HIGH; return L->T_max; }
    if (h == L->h_min) return L->T_min;
    if (h == L->h_max) return L->T_max;
    double t;
    if (L->kind == 0) {
        const double hs = pt_ss_h(1672.0), hl = pt_ss_h(1727.0);
        if (h >= hl) {
            t = 1727.0 + (h - hl) / 800.0;
        } else if (h >= hs) {
            const double cs = pt_ss_cp_solid(1672.0), w = 1727.0 - 1672.0;
            const double B = cs + PT_SS_LATENT / w, C = 0.5 * (800.0 - cs) / w, a = h - hs;
            t = 1672.0 + 2.0 * a / (B + pt_sqrt(B * B + 4.0 * C * a));
        } else {
            const double lo = L->T_min, hi = 1672.0;
            t = lo + (h - L->h_min) * ((hi - lo) / (hs - L->h_min));
            for (int it = 0; it < 8; ++it) {
                t = t - (pt_ss_h_solid(t) - h) / pt_ss_cp_solid(t);
                if (t < lo) t = lo;
                if (t > hi) t = hi;
            }
        }
    } else {
        const double t1 = L->c_Tm - 0.5 * L->c_dTm, t2 = L->c_Tm + 0.5 * L->c_dTm;
        const double h1 = pt_c_h(L, t1), h2 = pt_c_h(L, t2);
        if (h < h1) t = L->c_T0 + h / L->c_cps;
        else if (h < h2) t = t1 + (h - h1) / (0.5 * (L->c_cps + L->c_cpl) + L->c_L / L->c_dTm);
        else t = t2 + (h - h2) / L->c_cpl;
    }
    if (t < L->T_min) t = L->T_min;
    if (t > L->T_max) t = L->T_max;
    return pt_polish(L, h, t);
}

// --------------------------------------------------------------------------------------------------------
// Vapor law. ln(p_i / Pa) of the pure liquid element at T (linear in 1/T between JANAF rows).
PT_HD double pt_lnp(const pt_law* L, int el, double t, int* status) {
    const int n = L->n_knots[el];
    const double* T = L->knot_T[el];
    if (!(t >= T[0]) || !(t <= T[n - 1])) { *status = PT_ERR_NONFINITE; return -1e300; }
    int k = 0;
    while (k + 2 < n && T[k + 1] <= t) ++k;
    const double w = (1.0 / t - 1.0 / T[k]) / (1.0 / T[k + 1] - 1.0 / T[k]);
    return L->knot_lnp[el][k] + w * (L->knot_lnp[el][k + 1] - L->knot_lnp[el][k]);
}
// One-sided HKL number flux (atoms/m2/s) of element el from a surface at T with mole fraction x.
PT_HD double pt_flux(const pt_law* L, int el, double t, double x, int* status) {
    if (!L->vapor_enabled || x <= 0.0) return 0.0;
    const double p = pt_exp(pt_lnp(L, el, t, status));
    const double a = L->gamma[el] * x;
    return (1.0 - L->beta_r) * L->alpha_e * a * p / pt_sqrt(2.0 * PT_PI * L->atom_mass[el] * PT_KB * t);
}
