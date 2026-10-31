// phase_law_host.hpp -- host construction of the pq_law record (v1.1): JANAF vapor knots for Fe, Cr, Ni, Mn, Si, C,
// per-atom reference and internal energies, regular-solution interaction energies, Knudsen-layer table,
// pulse-heating liquid constants, enthalpy bounds and the law identity digest.
#pragma once
#include "phase_kernels.h"
#include "janaf_data.h"
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <openssl/sha.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ss304_phase_v1_1 {

inline constexpr double avogadro = 6.02214076e23;          // SI exact
inline constexpr double boltzmann = 1.380649e-23;          // SI exact
inline constexpr double elementary_charge = 1.602176634e-19;   // SI exact
inline constexpr double gas_constant = boltzmann * avogadro;
// CIAAW standard atomic weights (abridged): Fe, Cr, Ni, Mn, Si, C.
inline constexpr std::array<double, PQ_NEL> molar_mass_g = {55.845, 51.9961, 58.6934, 54.938043, 28.085, 12.011};
inline constexpr std::array<const char*, PQ_NEL> element_name = {"Fe", "Cr", "Ni", "Mn", "Si", "C"};
// v1 composition: Mills 2004 Table 1 typical 304 with Mn, Si, C lumped into Fe.
inline constexpr std::array<double, PQ_NEL> v1_mass_fraction = {0.717, 0.190, 0.093, 0.0, 0.0, 0.0};
// v1.1 composition: Mills 2004 Table 1 typical 304 (19 Cr, 9.3 Ni, 1.0 Mn, 0.5 Si, 0.04 C mass %, Fe balance).
inline constexpr std::array<double, PQ_NEL> nominal_mass_fraction = {0.7016, 0.190, 0.093, 0.010, 0.005, 0.0004};

struct LawOptions {
    int kind = 0;
    double T_min = 298.2, T_max = 5000.0;
    double c_T0 = 300, c_cps = 500, c_cpl = 500, c_Tm = 1000, c_dTm = 1, c_L = 0, c_ks = 20, c_kl = 20, c_rhoe = 1e-6;
    double c_rho = 0, c_mu = 0, c_sigma = 0, c_dsdT = 0;
    bool vapor = true;
    double alpha_e = 1.0, beta_r = 0.18;
    std::array<double, PQ_NEL> gamma = {1, 1, 1, 1, 1, 1};
    double eps0_anchor_K = 3000.0;
    // v1.1 modes (0 = v1 behaviour)
    int caloric = 0, activity = 0, gas = 0, knudsen = 0, radiation = 0;
    double emissivity = 0.37, T_env = 300.0;
    double T_critical = 0.0, explosion_fraction = 0.9;       // phase-explosion temperature = fraction x T_c (0: off)
    std::array<std::array<double, PQ_NEL>, PQ_NEL> omega{};  // regular-solution interaction energies (J/mol)
    double sulfur_wt = 0.0;                                   // mass % S for the surface tension (Mills Eq. 26)

    static LawOptions v1() { return LawOptions{}; }
    // All v1.1 extensions on (see DATA-SOURCES.json).
    static LawOptions v1_1() {
        LawOptions o;
        o.T_max = 6000.0;
        o.caloric = 1; o.activity = 1; o.gas = 1; o.knudsen = 1; o.radiation = 1;
        // Fe-Ni: activity coefficient of Ni at infinite dilution in liquid Fe 0.66 at 1600 C (Zellars et al. 1960);
        // regular-solution energy Omega = R T ln(0.66) at 1873 K. Fe-Cr, Fe-Mn and the remaining pairs: 0 (Gilby &
        // St. Pierre 1969: negligible interactions in iron-rich Fe-Cr-Ni; Roy & Hultgren 1965: Fe-Mn nearly ideal).
        const double w = gas_constant * 1873.0 * std::log(0.66);
        o.omega[0][2] = o.omega[2][0] = w;
        o.T_critical = 8000.0;   // lowest published estimate cited (SESAME 4272, 2000 K below Mattsson & Desjarlais 2004)
        return o;
    }
};

struct VaporKnot { double T, lnp_Pa, dHvap_J_mol; bool continued; };

inline std::vector<VaporKnot> janaf_knots(const janaf::Table& gas, const janaf::Table& cond, double T_lo, double T_end) {
    std::vector<VaporKnot> out;
    const auto& last = cond.rows[cond.n - 1];
    const double dH0 = gas.dfH298_kJ_mol - cond.dfH298_kJ_mol;
    for (int i = 0; i < gas.n; ++i) {
        const auto& g = gas.rows[i];
        if (g.T < T_lo) continue;
        double S = 0, H = 0; bool cont = false, found = false;
        if (g.T <= last.T) {
            for (int j = 0; j < cond.n; ++j) if (cond.rows[j].T == g.T) { S = cond.rows[j].S; H = cond.rows[j].HmH298; found = true; break; }
            if (!found) continue;
        } else {
            S = last.S + last.Cp * pq_log(g.T / last.T);
            H = last.HmH298 + last.Cp * (g.T - last.T) / 1000.0;
            cont = true;
        }
        const double dH = dH0 + (g.HmH298 - H);
        const double dG = 1000.0 * dH - g.T * (g.S - S);
        const double lnp_bar = -dG / (gas_constant * g.T);
        out.push_back({g.T, lnp_bar + pq_log(1e5), 1000.0 * dH, cont});
        if (g.T >= T_end) break;
    }
    return out;
}
inline const janaf::Table& gas_table(int el) {
    switch (el) { case 0: return janaf::fe_gas; case 1: return janaf::cr_gas; case 2: return janaf::ni_gas; case 3: return janaf::mn_gas; case 4: return janaf::si_gas; default: return janaf::c_gas; }
}
inline const janaf::Table& condensed_table(int el) {
    switch (el) { case 0: return janaf::fe_liquid; case 1: return janaf::cr_liquid; case 2: return janaf::ni_liquid; case 3: return janaf::mn_liquid; case 4: return janaf::si_liquid; default: return janaf::c_graphite; }
}
inline double knot_dHvap(const std::vector<VaporKnot>& k, double T) {
    for (const auto& x : k) if (x.T == T) return x.dHvap_J_mol;
    throw std::invalid_argument("vapor reference anchor is not a JANAF row");
}

// erf by its Taylor series (|x| <= 1.5; basic operations only, deterministic on every host).
inline double erf_series(double x) {
    double term = x, sum = x;
    for (int n = 1; n < 80; ++n) { term *= -x * x / n; sum += term / (2 * n + 1); }
    return 2.0 / std::sqrt(PQ_PI) * sum;
}
// Knight (1979) jump conditions for a monatomic vapor (gamma = 5/3) at Mach number M (formulas as reproduced in
// Schreter-Fleischhacker et al., arXiv:2501.18781, App. D, Eq. D.1-D.3).
struct KnudsenPoint { double M, Tr, rhor, pr, fnet, prec; };
inline KnudsenPoint knight(double M) {
    const double g = 5.0 / 3.0, m = std::sqrt(g / 2.0) * M, sp = std::sqrt(PQ_PI);
    const double a = sp * (g - 1.0) / (g + 1.0) * m / 2.0;
    const double s = std::sqrt(1.0 + a * a) - a;
    const double Tr = s * s;
    const double em2 = pq_exp(m * m), erfcm = 1.0 - erf_series(m);
    const double rhor = std::sqrt(1.0 / Tr) * ((m * m + 0.5) * em2 * erfcm - m / sp) + 0.5 / Tr * (1.0 - sp * m * em2 * erfcm);
    const double pr = rhor * Tr;
    const double fnet = rhor * M * std::sqrt(g * Tr) * std::sqrt(2.0 * PQ_PI);
    return {M, Tr, rhor, pr, fnet, pr * (1.0 + g * M * M)};
}

inline pq_law make_law(const LawOptions& o) {
    pq_law L;
    std::memset(&L, 0, sizeof L);
    L.kind = o.kind;
    for (int el = 0; el < PQ_NEL; ++el) L.atom_mass[el] = molar_mass_g[el] * 1e-3 / avogadro;
    L.caloric = o.caloric; L.activity = o.activity; L.gas = o.gas; L.knudsen = o.knudsen; L.radiation = o.radiation;
    if (!(o.caloric >= 0 && o.caloric <= 1 && o.activity >= 0 && o.activity <= 1 && o.gas >= 0 && o.gas <= 1 && o.knudsen >= 0 && o.knudsen <= 1 && o.radiation >= 0 && o.radiation <= 1))
        throw std::invalid_argument("law mode selection");
    L.emissivity = o.emissivity; L.T_env = o.T_env;
    if (o.radiation && !(o.emissivity > 0 && o.emissivity <= 1 && o.T_env > 0)) throw std::invalid_argument("emissivity / environment temperature");
    L.T_explosion = o.T_critical > 0 ? o.explosion_fraction * o.T_critical : 0.0;
    L.sulfur_wt = o.sulfur_wt;
    // pulse-heating liquid (Pichler et al. 2020): H - H(298.15 K) = -335 + 0.847 T kJ/kg, zero of this law at 300 K
    L.ph_cp = 847.0;
    L.ph_hA = -335000.0 + pq_ss_h_solid(298.15);
    L.ph_hs = pq_ss_h_solid(1672.0);
    L.ph_hl = L.ph_hA + L.ph_cp * 1727.0;
    L.ph_lorenz = PQ_PI * PQ_PI / 3.0 * (boltzmann / elementary_charge) * (boltzmann / elementary_charge);
    if (o.kind == 1) {
        if (!(o.c_cps > 0 && o.c_cpl > 0 && o.c_dTm > 0 && o.c_L >= 0 && o.c_ks > 0 && o.c_kl > 0 && o.c_rhoe > 0 && o.T_max > o.T_min))
            throw std::invalid_argument("constant benchmark law parameters");
        L.c_T0 = o.c_T0; L.c_cps = o.c_cps; L.c_cpl = o.c_cpl; L.c_Tm = o.c_Tm; L.c_dTm = o.c_dTm; L.c_L = o.c_L;
        L.c_ks = o.c_ks; L.c_kl = o.c_kl; L.c_rhoe = o.c_rhoe;
        L.c_rho = o.c_rho; L.c_mu = o.c_mu; L.c_sigma = o.c_sigma; L.c_dsdT = o.c_dsdT;
        L.T_min = o.T_min; L.T_max = o.T_max;
        L.vapor_enabled = 0;
        L.caloric = 0; L.activity = 0; L.gas = 0; L.knudsen = 0;
    } else {
        const double Tvap_end = o.T_max > 5000.0 ? 6000.0 : 5000.0;
        if (!(o.T_min >= 298.2 && o.T_max <= 6000.0 && o.T_max > o.T_min)) throw std::invalid_argument("SS304 phase law v1.1 domain is 298.2-6000 K");
        if (o.caloric == 0 && o.T_max > 5000.0) throw std::invalid_argument("the v1 caloric law is limited to 5000 K");
        L.T_min = o.T_min; L.T_max = o.T_max;
        L.vapor_enabled = o.vapor ? 1 : 0;
        if (!(o.alpha_e > 0 && o.alpha_e <= 1 && o.beta_r >= 0 && o.beta_r < 1)) throw std::invalid_argument("evaporation/recondensation coefficients");
        L.alpha_e = o.alpha_e; L.beta_r = o.beta_r;
        for (int i = 0; i < PQ_NEL; ++i) for (int j = 0; j < PQ_NEL; ++j) {
            if (!std::isfinite(o.omega[i][j]) || o.omega[i][j] != o.omega[j][i] || (i == j && o.omega[i][j] != 0)) throw std::invalid_argument("interaction energies must be finite, symmetric, zero on the diagonal");
            L.omega[i][j] = o.omega[i][j];
        }
        for (int el = 0; el < PQ_NEL; ++el) {
            if (!(o.gamma[el] > 0)) throw std::invalid_argument("activity coefficient");
            L.gamma[el] = o.gamma[el];
            const auto k = janaf_knots(gas_table(el), condensed_table(el), 298.15, Tvap_end);
            if (k.size() < 2 || k.size() > PQ_MAX_KNOTS || k.front().T > o.T_min || k.back().T < o.T_max) throw std::logic_error("vapor knot coverage");
            L.n_knots[el] = int(k.size());
            for (std::size_t i = 0; i < k.size(); ++i) { L.knot_T[el][i] = k[i].T; L.knot_lnp[el][i] = k[i].lnp_Pa; }
            const double TA = o.eps0_anchor_K;
            // v1 reference (gas = 0). With gas = 1, eps0 is the minimum of m h + dHvap/N_A - 5/2 kT over the domain
            // and E_int carries the rest at every JANAF row.
            if (!o.gas) {
                L.eps0[el] = L.atom_mass[el] * pq_ss_h(TA) + knot_dHvap(k, TA) / avogadro - 2.5 * boltzmann * TA;
            } else {
                double mn = INFINITY;
                std::vector<double> u(k.size());
                for (std::size_t i = 0; i < k.size(); ++i) {
                    const double T = k[i].T < o.T_min ? o.T_min : k[i].T;
                    u[i] = L.atom_mass[el] * pq_h(&L, T) + k[i].dHvap_J_mol / avogadro - 2.5 * boltzmann * k[i].T;
                    if (k[i].T >= o.T_min && k[i].T <= o.T_max) mn = std::min(mn, u[i]);
                }
                L.eps0[el] = mn;
                for (std::size_t i = 0; i < k.size(); ++i) L.knot_eint[el][i] = u[i] - mn;
            }
        }
        // Knudsen-layer table (Mach 0..1)
        for (int i = 0; i < PQ_KL_N; ++i) {
            const auto kp = knight(double(i) / (PQ_KL_N - 1));
            L.kl_pv[i] = kp.pr; L.kl_fnet[i] = kp.fnet; L.kl_prec[i] = kp.prec;
        }
    }
    L.h_min = pq_h(&L, L.T_min);
    L.h_max = pq_h(&L, L.T_max);
    return L;
}

inline std::string hex(const unsigned char* d, std::size_t n) {
    std::ostringstream o;
    for (std::size_t i = 0; i < n; ++i) o << std::hex << std::setw(2) << std::setfill('0') << unsigned(d[i]);
    return o.str();
}
inline std::string sha256(const void* p, std::size_t n) { unsigned char d[32]; SHA256(static_cast<const unsigned char*>(p), n, d); return hex(d, 32); }

struct Encoder {
    std::vector<unsigned char> bytes;
    void u64(std::uint64_t v) { for (int i = 0; i < 8; ++i) bytes.push_back(static_cast<unsigned char>(v >> (8 * i))); }
    void i64(std::int64_t v) { u64(static_cast<std::uint64_t>(v)); }
    void f64(double v) { std::uint64_t u; std::memcpy(&u, &v, 8); u64(u); }
    void text(const std::string& s) { u64(s.size()); bytes.insert(bytes.end(), s.begin(), s.end()); }
    void w(pq_i128 v) { u64(v.lo); i64(v.hi); }
    std::string digest() const { return sha256(bytes.data(), bytes.size()); }
};

inline std::string law_sha256(const pq_law& L) {
    Encoder e;
    e.text("ss304-phase-law-v1.1");
    // the whole record, field by field (no padding bytes)
    e.i64(L.kind);
    for (double v : {L.T_min, L.T_max, L.h_min, L.h_max, L.c_T0, L.c_cps, L.c_cpl, L.c_Tm, L.c_dTm, L.c_L, L.c_ks, L.c_kl, L.c_rhoe,
                     L.c_rho, L.c_mu, L.c_sigma, L.c_dsdT, L.ph_hA, L.ph_cp, L.ph_hs, L.ph_hl, L.ph_lorenz, L.emissivity, L.T_env,
                     L.T_explosion, L.alpha_e, L.beta_r, L.sulfur_wt}) e.f64(v);
    for (int v : {L.caloric, L.activity, L.gas, L.knudsen, L.radiation, L.vapor_enabled}) e.i64(v);
    for (int el = 0; el < PQ_NEL; ++el) {
        e.i64(L.n_knots[el]);
        for (int i = 0; i < L.n_knots[el]; ++i) { e.f64(L.knot_T[el][i]); e.f64(L.knot_lnp[el][i]); e.f64(L.knot_eint[el][i]); }
        e.f64(L.atom_mass[el]); e.f64(L.eps0[el]); e.f64(L.gamma[el]);
        for (int j = 0; j < PQ_NEL; ++j) e.f64(L.omega[el][j]);
    }
    for (int i = 0; i < PQ_KL_N; ++i) { e.f64(L.kl_pv[i]); e.f64(L.kl_fnet[i]); e.f64(L.kl_prec[i]); }
    for (double T : {298.2, 300.0, 1000.0, 1500.0, 1672.0, 1700.0, 1727.0, 1800.0, 3000.0, 5000.0}) {
        e.f64(pq_h(&L, T)); e.f64(pq_cap(&L, T)); e.f64(pq_rhoe(&L, T)); e.f64(pq_k(&L, T));
    }
    return e.digest();
}

}  // namespace ss304_phase_v1_1
