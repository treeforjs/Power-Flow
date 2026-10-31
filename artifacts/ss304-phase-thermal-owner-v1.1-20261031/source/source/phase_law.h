// phase_law.h -- SS304 phase law v1.1, evaluated with identical arithmetic on the host and in NVRTC kernels.
// Provenance of every number: DATA-SOURCES.json. Every v1.1 extension is selected by a mode field of pq_law;
// with all modes 0 the arithmetic is the v1 law operation for operation (bit-identical, tests/test_regression).
//
// v1 law (modes 0): see the v1 header text, unchanged:
//   solid 298.2-1672 K Mills 2004 Eq. 10/11 (sealed expressions); mushy 1672-1727 K with L_f = 160 J/(K kg) x 1727 K;
//   liquid Cp 800; transport Ho & Chu 1977 to 1800 K, then Mills slope and constant Lorenz ratio; vapor of Fe, Cr, Ni
//   from JANAF; Raoult with constant gamma; HKL flux with fixed recondensation (1 - beta_r).
// v1.1 extensions:
//   caloric = 1  liquid (T >= 1727 K) from the pulse-heating measurements of Pichler, Simonds, Sowards, Pottlacher,
//                J. Mater. Sci. 55 (2020) 4081 (NIST SRM 1155a, Cr-Ni-Mo austenitic steel, 1708-2900 K):
//                H(T) - H(298.15 K) = -335 + 0.847 T kJ/kg (2.9 %, k = 2); rho_e,corr = 1.154 + 1.893e-4 T uOhm m
//                (4.5 %); D = 8065 - 0.661 T kg/m3 (3 %); thermal conductivity from the Wiedemann-Franz law with
//                the Sommerfeld Lorenz number L0 = (pi^2/3)(k_B/e)^2; linear in T across the mushy zone; solid
//                unchanged (Mills / Ho & Chu, 304); beyond 2900 K the linear fits are extrapolated.
//   activity = 1 regular solution over Fe, Cr, Ni, Mn, Si, C with binary interaction energies Omega_ij (J/mol):
//                ln gamma_k = [sum_j Omega_kj x_j - sum_{i<j} Omega_ij x_i x_j] / (R T)
//   gas = 1      the energy removed per evaporated atom is eps0_i + E_int,i(T) + kinetic, with E_int tabulated at
//                the JANAF rows so that eps0 + E_int + 5/2 kT - m h(T) equals the JANAF vaporization enthalpy
//                (gas electronic excitation and the element-versus-alloy caloric difference are carried by E_int)
//   knudsen = 1  Knight (1979) jump conditions with back pressure: net flux fraction and recoil pressure from the
//                Mach number at which the Knudsen-layer pressure equals the back pressure (sonic for vacuum)
//   radiation=1  grey-body loss eps sigma (T^4 - T_env^4) on exposed faces
#pragma once
#include "det_math.h"

#define PQ_MAX_KNOTS 72
#define PQ_NEL 6   // Fe, Cr, Ni, Mn, Si, C
#define PQ_KL_N 1025

// Status codes written by the per-cell arithmetic (0 = ok).
#define PQ_OK 0
#define PQ_ERR_T_LOW 1          // specific energy below h(T_min)
#define PQ_ERR_T_HIGH 2         // specific energy above h(T_max)
#define PQ_ERR_NONFINITE 3
#define PQ_ERR_FIXED_RANGE 4    // fixed-point conversion out of range / overflow
#define PQ_ERR_DEPLETED 5       // a cell lost more than the admitted fraction of its reference mass in one step
#define PQ_ERR_NEGATIVE_INV 6   // an element inventory would become negative
#define PQ_ERR_FLOW 7           // melt flow would move more than the admitted fraction of a donor cell in one sub-step
#define PQ_ERR_PHASE_EXPLOSION 8 // a cell reached the phase-explosion temperature (0.9 T_c)

struct pq_law {
    int kind;                 // 0: SS304 phase law, 1: constant-property benchmark law
    double T_min, T_max, h_min, h_max;
    // constant-property benchmark law (kind 1): h = 0 at c_T0
    double c_T0, c_cps, c_cpl, c_Tm, c_dTm, c_L, c_ks, c_kl, c_rhoe;
    double c_rho, c_mu, c_sigma, c_dsdT;               // benchmark melt properties (kind 1)
    // v1.1 modes
    int caloric, activity, gas, knudsen, radiation;
    double ph_hA, ph_cp, ph_hs, ph_hl, ph_lorenz;      // pulse-heating liquid: h = ph_hA + ph_cp T; mushy ph_hs -> ph_hl
    double emissivity, T_env, T_explosion;             // radiation; phase-explosion temperature (0.9 T_c)
    // vapor law
    int vapor_enabled;
    int n_knots[PQ_NEL];
    double knot_T[PQ_NEL][PQ_MAX_KNOTS], knot_lnp[PQ_NEL][PQ_MAX_KNOTS], knot_eint[PQ_NEL][PQ_MAX_KNOTS];
    double atom_mass[PQ_NEL], eps0[PQ_NEL], gamma[PQ_NEL];
    double omega[PQ_NEL][PQ_NEL];                      // J/mol, symmetric, zero diagonal
    double alpha_e, beta_r;
    // Knudsen-layer table over the Mach number grid M_k = k/(PQ_KL_N-1): p_v/p_s (decreasing), net/one-sided flux,
    // recoil/p_s
    double kl_pv[PQ_KL_N], kl_fnet[PQ_KL_N], kl_prec[PQ_KL_N];
    // melt (Mills 2004 Eq. 26-28)
    double sulfur_wt;
};

#define PQ_KB 1.380649e-23
#define PQ_PI 3.14159265358979323846
#define PQ_SIGMA_SB 5.670374419e-8
#define PQ_R_GAS (1.380649e-23 * 6.02214076e23)

// --------------------------------------------------------------------------------------------------------
// SS304 phase law v1 (sealed expressions, mode 0).
PQ_HD double pq_ss_cp_solid(double t) { return 472.0 + 0.136 * t - 2820000.0 / (t * t); }
PQ_HD double pq_ss_h_solid(double t) { return (t - 300.0) * (472.0 + 0.068 * (t + 300.0) - 2820000.0 / (t * 300.0)); }
PQ_HD double pq_ss_h_mushy_sensible(double dt) {
    const double cs = pq_ss_cp_solid(1672.0);
    return pq_ss_h_solid(1672.0) + cs * dt + 0.5 * (800.0 - cs) * dt * dt / (1727.0 - 1672.0);
}
#define PQ_SS_LATENT (160.0 * 1727.0)
PQ_HD double pq_ss_h(double t) {
    const double width = 1727.0 - 1672.0;
    if (t < 1672.0) return pq_ss_h_solid(t);
    if (t < 1727.0) { const double dt = t - 1672.0; return pq_ss_h_mushy_sensible(dt) + PQ_SS_LATENT * (dt / width); }
    return pq_ss_h_mushy_sensible(width) + PQ_SS_LATENT + 800.0 * (t - 1727.0);
}
PQ_HD double pq_ss_cap(double t) {
    if (t >= 1727.0) return 800.0;
    if (t >= 1672.0) {
        const double f = (t - 1672.0) / (1727.0 - 1672.0);
        const double cp = pq_ss_cp_solid(1672.0) + (800.0 - pq_ss_cp_solid(1672.0)) * f;
        return cp + PQ_SS_LATENT / (1727.0 - 1672.0);
    }
    return pq_ss_cp_solid(t);
}
PQ_HD double pq_ss_liquid_fraction(double t) {
    if (t >= 1727.0) return 1.0;
    if (t >= 1672.0) return (t - 1672.0) / (1727.0 - 1672.0);
    return 0.0;
}
PQ_HD void pq_ss_transport(double t, double* rho, double* k) {
    if (t < 300.0) {
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
    const double kk = 28.99 + 0.01 * (t - 1800.0);
    const double lorenz = (28.99 * 151.7e-8) / 1800.0;
    *k = kk;
    *rho = lorenz * t / kk;
}
PQ_HD double pq_ss_density(double t) {
    if (t < 1672.0) return 7921.0 - 0.5 * (t - 298.2);
    if (t >= 1727.0) return 6908.0 - 0.86 * (t - 1723.0);
    const double f = (t - 1672.0) / (1727.0 - 1672.0);
    const double vs = 1.0 / (7921.0 - 0.5 * (1672.0 - 298.2)), vl = 1.0 / (6908.0 - 0.86 * (1727.0 - 1723.0));
    return 1.0 / ((1.0 - f) * vs + f * vl);
}

// --------------------------------------------------------------------------------------------------------
// Pulse-heating liquid (caloric = 1).
PQ_HD double pq_ph_h(const pq_law* L, double t) {
    if (t < 1672.0) return pq_ss_h_solid(t);
    if (t < 1727.0) return L->ph_hs + (L->ph_hl - L->ph_hs) * ((t - 1672.0) / (1727.0 - 1672.0));
    return L->ph_hA + L->ph_cp * t;
}
PQ_HD double pq_ph_cap(const pq_law* L, double t) {
    if (t >= 1727.0) return L->ph_cp;
    if (t >= 1672.0) return (L->ph_hl - L->ph_hs) / (1727.0 - 1672.0);
    return pq_ss_cp_solid(t);
}
PQ_HD double pq_ph_rho_liquid(double t) { return (1.154 + 1.893e-4 * t) * 1e-6; }
PQ_HD void pq_ph_transport(const pq_law* L, double t, double* rho, double* k) {
    if (t < 1672.0) { pq_ss_transport(t, rho, k); return; }
    const double rl = pq_ph_rho_liquid(t < 1727.0 ? 1727.0 : t);
    const double kl = L->ph_lorenz * (t < 1727.0 ? 1727.0 : t) / rl;
    if (t >= 1727.0) { *rho = rl; *k = kl; return; }
    double r0, k0; pq_ss_transport(1672.0, &r0, &k0);
    const double f = (t - 1672.0) / (1727.0 - 1672.0);
    *rho = r0 + f * (rl - r0);
    *k = k0 + f * (kl - k0);
}
PQ_HD double pq_ph_density(double t) {
    if (t < 1672.0) return 7921.0 - 0.5 * (t - 298.2);
    if (t >= 1727.0) return 8065.0 - 0.661 * t;
    const double f = (t - 1672.0) / (1727.0 - 1672.0);
    const double vs = 1.0 / (7921.0 - 0.5 * (1672.0 - 298.2)), vl = 1.0 / (8065.0 - 0.661 * 1727.0);
    return 1.0 / ((1.0 - f) * vs + f * vl);
}

// --------------------------------------------------------------------------------------------------------
// Constant-property benchmark law (kind 1).
PQ_HD double pq_c_h(const pq_law* L, double t) {
    const double t1 = L->c_Tm - 0.5 * L->c_dTm, t2 = L->c_Tm + 0.5 * L->c_dTm;
    const double h1 = L->c_cps * (t1 - L->c_T0);
    if (t < t1) return L->c_cps * (t - L->c_T0);
    const double cm = 0.5 * (L->c_cps + L->c_cpl);
    if (t < t2) return h1 + (cm + L->c_L / L->c_dTm) * (t - t1);
    const double h2 = h1 + cm * L->c_dTm + L->c_L;
    return h2 + L->c_cpl * (t - t2);
}
PQ_HD double pq_c_cap(const pq_law* L, double t) {
    const double t1 = L->c_Tm - 0.5 * L->c_dTm, t2 = L->c_Tm + 0.5 * L->c_dTm;
    if (t < t1) return L->c_cps;
    if (t < t2) return 0.5 * (L->c_cps + L->c_cpl) + L->c_L / L->c_dTm;
    return L->c_cpl;
}
PQ_HD double pq_c_fl(const pq_law* L, double t) {
    const double t1 = L->c_Tm - 0.5 * L->c_dTm, t2 = L->c_Tm + 0.5 * L->c_dTm;
    if (t < t1) return 0.0;
    if (t < t2) return (t - t1) / L->c_dTm;
    return 1.0;
}
PQ_HD double pq_c_k(const pq_law* L, double t) {
    const double f = pq_c_fl(L, t);
    return L->c_ks + f * (L->c_kl - L->c_ks);
}

// --------------------------------------------------------------------------------------------------------
// Law dispatch.
PQ_HD double pq_h(const pq_law* L, double t) { return L->kind == 0 ? (L->caloric == 0 ? pq_ss_h(t) : pq_ph_h(L, t)) : pq_c_h(L, t); }
PQ_HD double pq_cap(const pq_law* L, double t) { return L->kind == 0 ? (L->caloric == 0 ? pq_ss_cap(t) : pq_ph_cap(L, t)) : pq_c_cap(L, t); }
PQ_HD double pq_fl(const pq_law* L, double t) { return L->kind == 0 ? pq_ss_liquid_fraction(t) : pq_c_fl(L, t); }
PQ_HD double pq_k(const pq_law* L, double t) {
    if (L->kind != 0) return pq_c_k(L, t);
    double r, k;
    if (L->caloric == 0) pq_ss_transport(t, &r, &k); else pq_ph_transport(L, t, &r, &k);
    return k;
}
PQ_HD double pq_rhoe(const pq_law* L, double t) {
    if (L->kind != 0) return L->c_rhoe;
    double r, k;
    if (L->caloric == 0) pq_ss_transport(t, &r, &k); else pq_ph_transport(L, t, &r, &k);
    return r;
}
PQ_HD double pq_density(const pq_law* L, double t) {
    if (L->kind != 0) return 0.0;
    return L->caloric == 0 ? pq_ss_density(t) : pq_ph_density(t);
}

PQ_HD double pq_polish(const pq_law* L, double h, double t) {
    for (int it = 0; it < 4; ++it) {
        const double e = pq_fabs(pq_h(L, t) - h);
        const double up = pq_next_up(t), dn = pq_next_down(t);
        double best = t, be = e;
        if (up <= L->T_max) { const double eu = pq_fabs(pq_h(L, up) - h); if (eu < be) { best = up; be = eu; } }
        if (dn >= L->T_min) { const double ed = pq_fabs(pq_h(L, dn) - h); if (ed < be) { best = dn; be = ed; } }
        if (best == t) break;
        t = best;
    }
    return t;
}
PQ_HD double pq_T_of_h(const pq_law* L, double h, int* status) {
    if (!pq_finite(h)) { *status = PQ_ERR_NONFINITE; return 0.0; }
    if (h < L->h_min) { *status = PQ_ERR_T_LOW; return L->T_min; }
    if (h > L->h_max) { *status = PQ_ERR_T_HIGH; return L->T_max; }
    if (h == L->h_min) return L->T_min;
    if (h == L->h_max) return L->T_max;
    double t;
    if (L->kind == 0) {
        const double hs = pq_ss_h(1672.0);
        const double hl = L->caloric == 0 ? pq_ss_h(1727.0) : L->ph_hl;
        if (h >= hl) {
            t = L->caloric == 0 ? 1727.0 + (h - hl) / 800.0 : (h - L->ph_hA) / L->ph_cp;
        } else if (h >= hs) {
            if (L->caloric == 0) {
                const double cs = pq_ss_cp_solid(1672.0), w = 1727.0 - 1672.0;
                const double B = cs + PQ_SS_LATENT / w, C = 0.5 * (800.0 - cs) / w, a = h - hs;
                t = 1672.0 + 2.0 * a / (B + pq_sqrt(B * B + 4.0 * C * a));
            } else {
                t = 1672.0 + (1727.0 - 1672.0) * ((h - L->ph_hs) / (L->ph_hl - L->ph_hs));
            }
        } else {
            const double lo = L->T_min, hi = 1672.0;
            t = lo + (h - L->h_min) * ((hi - lo) / (hs - L->h_min));
            for (int it = 0; it < 8; ++it) {
                t = t - (pq_ss_h_solid(t) - h) / pq_ss_cp_solid(t);
                if (t < lo) t = lo;
                if (t > hi) t = hi;
            }
        }
    } else {
        const double t1 = L->c_Tm - 0.5 * L->c_dTm, t2 = L->c_Tm + 0.5 * L->c_dTm;
        const double h1 = pq_c_h(L, t1), h2 = pq_c_h(L, t2);
        if (h < h1) t = L->c_T0 + h / L->c_cps;
        else if (h < h2) t = t1 + (h - h1) / (0.5 * (L->c_cps + L->c_cpl) + L->c_L / L->c_dTm);
        else t = t2 + (h - h2) / L->c_cpl;
    }
    if (t < L->T_min) t = L->T_min;
    if (t > L->T_max) t = L->T_max;
    return pq_polish(L, h, t);
}

// --------------------------------------------------------------------------------------------------------
// Vapor law. ln(p_i / Pa) of the pure element in its standard state at T (linear in 1/T between JANAF rows).
PQ_HD double pq_lnp(const pq_law* L, int el, double t, int* status) {
    const int n = L->n_knots[el];
    const double* T = L->knot_T[el];
    if (!(t >= T[0]) || !(t <= T[n - 1])) { *status = PQ_ERR_NONFINITE; return -1e300; }
    int k = 0;
    while (k + 2 < n && T[k + 1] <= t) ++k;
    const double w = (1.0 / t - 1.0 / T[k]) / (1.0 / T[k + 1] - 1.0 / T[k]);
    return L->knot_lnp[el][k] + w * (L->knot_lnp[el][k + 1] - L->knot_lnp[el][k]);
}
// Internal energy (J per atom) carried by an atom evaporated at T (gas = 1), linear in T between JANAF rows.
PQ_HD double pq_eint(const pq_law* L, int el, double t) {
    if (!L->gas) return 0.0;
    const int n = L->n_knots[el];
    const double* T = L->knot_T[el];
    int k = 0;
    while (k + 2 < n && T[k + 1] <= t) ++k;
    const double w = (t - T[k]) / (T[k + 1] - T[k]);
    return L->knot_eint[el][k] + w * (L->knot_eint[el][k + 1] - L->knot_eint[el][k]);
}
// v1 one-sided HKL number flux with fixed recondensation (mode 0 arithmetic).
PQ_HD double pq_flux(const pq_law* L, int el, double t, double x, int* status) {
    if (!L->vapor_enabled || x <= 0.0) return 0.0;
    const double p = pq_exp(pq_lnp(L, el, t, status));
    const double a = L->gamma[el] * x;
    return (1.0 - L->beta_r) * L->alpha_e * a * p / pq_sqrt(2.0 * PQ_PI * L->atom_mass[el] * PQ_KB * t);
}
// Activity coefficient of element el in the regular-solution model at composition x (activity = 1).
PQ_HD double pq_gamma(const pq_law* L, int el, double t, const double* x) {
    if (!L->activity) return L->gamma[el];
    double s1 = 0.0, s2 = 0.0;
    for (int j = 0; j < PQ_NEL; ++j) s1 += L->omega[el][j] * x[j];
    for (int i = 0; i < PQ_NEL; ++i) for (int j = i + 1; j < PQ_NEL; ++j) s2 += L->omega[i][j] * x[i] * x[j];
    return pq_exp((s1 - s2) / (PQ_R_GAS * t));
}
// Knudsen layer: net-flux fraction and recoil pressure ratio for back-pressure ratio r = p_back / p_s.
PQ_HD void pq_knudsen(const pq_law* L, double r, double* fnet, double* prec) {
    const double* pv = L->kl_pv;
    if (!(r < 1.0)) { *fnet = 0.0; *prec = 1.0; return; }
    if (r <= pv[PQ_KL_N - 1]) { *fnet = L->kl_fnet[PQ_KL_N - 1]; *prec = L->kl_prec[PQ_KL_N - 1]; return; }
    int lo = 0, hi = PQ_KL_N - 1;                         // pv[lo] > r >= pv[hi]
    while (hi - lo > 1) { const int mid = (lo + hi) / 2; if (pv[mid] > r) lo = mid; else hi = mid; }
    const double w = (pv[lo] - r) / (pv[lo] - pv[hi]);
    *fnet = L->kl_fnet[lo] + w * (L->kl_fnet[hi] - L->kl_fnet[lo]);
    *prec = L->kl_prec[lo] + w * (L->kl_prec[hi] - L->kl_prec[lo]);
}
// Total equilibrium vapor pressure over a surface of composition x at T (Pa).
PQ_HD double pq_psat_total(const pq_law* L, double t, const double* x, int* status) {
    double ps = 0.0;
    for (int el = 0; el < PQ_NEL; ++el) if (x[el] > 0.0) ps += pq_gamma(L, el, t, x) * x[el] * pq_exp(pq_lnp(L, el, t, status));
    return ps;
}
// Net-flux fraction of the one-sided HKL flux: 1 - beta_r (knudsen = 0, the v1 value) or the Knight net fraction for
// the back pressure pb (Pa).
PQ_HD double pq_flux_fraction(const pq_law* L, double t, const double* x, double pb, int* status) {
    if (!L->knudsen) return 1.0 - L->beta_r;
    const double ps = pq_psat_total(L, t, x, status);
    double f, prec;
    pq_knudsen(L, ps > 0.0 ? pb / ps : 1.0, &f, &prec);
    return f;
}
// Net HKL number flux (atoms/m2/s) of element el for a surface of composition x at T, with net fraction f.
// With activity = 0 and f = 1 - beta_r this is the v1 expression operation for operation.
PQ_HD double pq_flux_f(const pq_law* L, int el, double t, const double* x, double f, int* status) {
    if (!L->vapor_enabled || x[el] <= 0.0) return 0.0;
    const double p = pq_exp(pq_lnp(L, el, t, status));
    const double a = pq_gamma(L, el, t, x) * x[el];
    return f * L->alpha_e * a * p / pq_sqrt(2.0 * PQ_PI * L->atom_mass[el] * PQ_KB * t);
}
// Recoil pressure on the surface (Pa); without the Knudsen-layer model the sonic value is used.
PQ_HD double pq_recoil(const pq_law* L, double t, const double* x, double pb, int* status) {
    if (!L->vapor_enabled) return 0.0;
    const double ps = pq_psat_total(L, t, x, status);
    if (!(ps > 0.0)) return 0.0;
    double f, prec;
    pq_knudsen(L, L->knudsen ? pb / ps : 0.0, &f, &prec);
    return prec * ps;
}

// --------------------------------------------------------------------------------------------------------
// Liquid surface tension (N/m) and its temperature derivative (Mills 2004 Eq. 26, 27; S in mass %), viscosity
// (Pa s, Mills 2004 Eq. 28).
PQ_HD void pq_surface_tension(const pq_law* L, double t, double* g, double* dgdT) {
    const double S = L->sulfur_wt;
    if (!(S > 0.0)) { *g = (1840.0 - 0.4 * (t - 1823.0)) * 1e-3; *dgdT = -0.4e-3; return; }
    const double B = 1.0 + pq_exp(28798.0 / t - 8.5647) * (0.68 * S);
    *g = (1840.0 - 0.4 * (t - 1823.0) - 0.056 * t * pq_log(B)) * 1e-3;
    *dgdT = (-0.4 - 0.056 * (28798.0 * (1.0 - B) / (B * t) + pq_log(B))) * 1e-3;
}
PQ_HD double pq_viscosity(double t) { return 1e-3 * pq_exp(-2.396 + 7950.0 / t); }
