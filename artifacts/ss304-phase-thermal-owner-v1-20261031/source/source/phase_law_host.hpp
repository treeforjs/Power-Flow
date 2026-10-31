// phase_law_host.hpp -- host construction of the pt_law record: JANAF vapor knots, per-atom reference energies,
// enthalpy bounds and the law identity digest.
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

namespace ss304_phase_v1 {

inline constexpr double avogadro = 6.02214076e23;          // SI exact
inline constexpr double boltzmann = 1.380649e-23;          // SI exact
inline constexpr double gas_constant = boltzmann * avogadro;
inline constexpr double atomic_mass_unit = 1.66053906660e-27;   // CODATA 2018
// CIAAW standard atomic weights (abridged values used by the v2 vapor model).
inline constexpr std::array<double, 3> molar_mass_g = {55.845, 51.9961, 58.6934};
inline constexpr std::array<const char*, 3> element_name = {"Fe", "Cr", "Ni"};
// Nominal 304 composition: Mills et al. 2004 Table 1 (typical 304: 19 Cr, 9.3 Ni, 1.0 Mn, 0.5 Si, 0.04 C mass %);
// Fe carries the balance including Mn, Si and C (minor elements are not tracked separately in v1).
inline constexpr std::array<double, 3> nominal_mass_fraction = {0.717, 0.190, 0.093};

struct LawOptions {
    int kind = 0;                              // 0 SS304 phase law v1, 1 constant-property benchmark law
    double T_min = 298.2, T_max = 5000.0;      // kind 0 domain (refusals only outside it)
    // kind 1 parameters
    double c_T0 = 300, c_cps = 500, c_cpl = 500, c_Tm = 1000, c_dTm = 1, c_L = 0, c_ks = 20, c_kl = 20, c_rhoe = 1e-6;
    // vapor (kind 0)
    bool vapor = true;
    double alpha_e = 1.0;                      // evaporation coefficient (ideal Langmuir upper bound)
    double beta_r = 0.18;                      // recondensation fraction of the emitted flux (Anisimov 1968 sonic Knudsen layer)
    std::array<double, 3> gamma = {1.0, 1.0, 1.0};   // Raoult activity coefficients
    double eps0_anchor_K = 3000.0;             // anchor of the constant per-atom vapor reference energies
};

struct VaporKnot { double T, lnp_Pa, dHvap_J_mol; bool liquid_continued; };

// ln(p/Pa) and the vaporization enthalpy of the pure liquid element at the JANAF gas-table rows (<= T_end).
inline std::vector<VaporKnot> janaf_knots(const janaf::Table& gas, const janaf::Table& liq, double T_lo, double T_end) {
    std::vector<VaporKnot> out;
    const auto& last = liq.rows[liq.n - 1];
    const double dH0 = gas.dfH298_kJ_mol - liq.dfH298_kJ_mol;   // kJ/mol, g - l at 298.15 K
    for (int i = 0; i < gas.n; ++i) {
        const auto& g = gas.rows[i];
        if (g.T < T_lo) continue;
        double S = 0, H = 0; bool cont = false, found = false;
        if (g.T <= last.T) {
            for (int j = 0; j < liq.n; ++j) if (liq.rows[j].T == g.T) { S = liq.rows[j].S; H = liq.rows[j].HmH298; found = true; break; }
            if (!found) continue;
        } else {
            // constant-Cp continuation of the JANAF liquid beyond its last row
            S = last.S + last.Cp * pt_log(g.T / last.T);
            H = last.HmH298 + last.Cp * (g.T - last.T) / 1000.0;
            cont = true;
        }
        const double dH = dH0 + (g.HmH298 - H);               // kJ/mol
        const double dG = 1000.0 * dH - g.T * (g.S - S);       // J/mol
        const double lnp_bar = -dG / (gas_constant * g.T);
        out.push_back({g.T, lnp_bar + pt_log(1e5), 1000.0 * dH, cont});
        if (g.T >= T_end) break;
    }
    return out;
}

inline const janaf::Table& gas_table(int el) { return el == 0 ? janaf::fe_gas : el == 1 ? janaf::cr_gas : janaf::ni_gas; }
inline const janaf::Table& liquid_table(int el) { return el == 0 ? janaf::fe_liquid : el == 1 ? janaf::cr_liquid : janaf::ni_liquid; }

inline double knot_dHvap(const std::vector<VaporKnot>& k, double T) {
    for (const auto& x : k) if (x.T == T) return x.dHvap_J_mol;
    throw std::invalid_argument("vapor reference anchor is not a JANAF row");
}

inline pt_law make_law(const LawOptions& o) {
    pt_law L;
    std::memset(&L, 0, sizeof L);
    L.kind = o.kind;
    for (int el = 0; el < 3; ++el) L.atom_mass[el] = molar_mass_g[el] * 1e-3 / avogadro;   // kg per atom (both kinds)
    if (o.kind == 1) {
        if (!(o.c_cps > 0 && o.c_cpl > 0 && o.c_dTm > 0 && o.c_L >= 0 && o.c_ks > 0 && o.c_kl > 0 && o.c_rhoe > 0 && o.T_max > o.T_min))
            throw std::invalid_argument("constant benchmark law parameters");
        L.c_T0 = o.c_T0; L.c_cps = o.c_cps; L.c_cpl = o.c_cpl; L.c_Tm = o.c_Tm; L.c_dTm = o.c_dTm; L.c_L = o.c_L;
        L.c_ks = o.c_ks; L.c_kl = o.c_kl; L.c_rhoe = o.c_rhoe;
        L.T_min = o.T_min; L.T_max = o.T_max;
        L.vapor_enabled = 0;
    } else {
        if (!(o.T_min >= 298.2 && o.T_max <= 5000.0 && o.T_max > o.T_min)) throw std::invalid_argument("SS304 phase law v1 domain is 298.2-5000 K");
        L.T_min = o.T_min; L.T_max = o.T_max;
        L.vapor_enabled = o.vapor ? 1 : 0;
        if (!(o.alpha_e > 0 && o.alpha_e <= 1 && o.beta_r >= 0 && o.beta_r < 1)) throw std::invalid_argument("evaporation/recondensation coefficients");
        L.alpha_e = o.alpha_e; L.beta_r = o.beta_r;
        for (int el = 0; el < 3; ++el) {
            if (!(o.gamma[el] > 0)) throw std::invalid_argument("activity coefficient");
            L.gamma[el] = o.gamma[el];
            const auto k = janaf_knots(gas_table(el), liquid_table(el), 298.15, 5000.0);
            if (k.size() < 2 || k.size() > PT_MAX_KNOTS || k.front().T > o.T_min || k.back().T < o.T_max) throw std::logic_error("vapor knot coverage");
            L.n_knots[el] = int(k.size());
            for (std::size_t i = 0; i < k.size(); ++i) { L.knot_T[el][i] = k[i].T; L.knot_lnp[el][i] = k[i].lnp_Pa; }
            // Constant per-atom vapor reference: eps0 = m h(T_A) + dHvap(T_A)/N_A - 5/2 k T_A, so that the energy
            // removed per evaporated atom (eps0 + kinetic) equals the JANAF vaporization enthalpy at T_A.
            const double TA = o.eps0_anchor_K;
            L.eps0[el] = L.atom_mass[el] * pt_ss_h(TA) + knot_dHvap(k, TA) / avogadro - 2.5 * boltzmann * TA;
        }
    }
    L.h_min = pt_h(&L, L.T_min);
    L.h_max = pt_h(&L, L.T_max);
    return L;
}

inline std::string hex(const unsigned char* d, std::size_t n) {
    std::ostringstream o;
    for (std::size_t i = 0; i < n; ++i) o << std::hex << std::setw(2) << std::setfill('0') << unsigned(d[i]);
    return o.str();
}
inline std::string sha256(const void* p, std::size_t n) { unsigned char d[32]; SHA256(static_cast<const unsigned char*>(p), n, d); return hex(d, 32); }

// Canonical byte encoder for identities (little-endian words, length-prefixed text).
struct Encoder {
    std::vector<unsigned char> bytes;
    void u64(std::uint64_t v) { for (int i = 0; i < 8; ++i) bytes.push_back(static_cast<unsigned char>(v >> (8 * i))); }
    void i64(std::int64_t v) { u64(static_cast<std::uint64_t>(v)); }
    void f64(double v) { std::uint64_t u; std::memcpy(&u, &v, 8); u64(u); }
    void text(const std::string& s) { u64(s.size()); bytes.insert(bytes.end(), s.begin(), s.end()); }
    void w(pt_i128 v) { u64(v.lo); i64(v.hi); }
    std::string digest() const { return sha256(bytes.data(), bytes.size()); }
};

inline std::string law_sha256(const pt_law& L) {
    Encoder e;
    e.text("ss304-phase-law-v1");
    e.i64(L.kind);
    for (double v : {L.T_min, L.T_max, L.h_min, L.h_max, L.c_T0, L.c_cps, L.c_cpl, L.c_Tm, L.c_dTm, L.c_L, L.c_ks, L.c_kl, L.c_rhoe, L.alpha_e, L.beta_r}) e.f64(v);
    e.i64(L.vapor_enabled);
    for (int el = 0; el < 3; ++el) {
        e.i64(L.n_knots[el]);
        for (int i = 0; i < L.n_knots[el]; ++i) { e.f64(L.knot_T[el][i]); e.f64(L.knot_lnp[el][i]); }
        e.f64(L.atom_mass[el]); e.f64(L.eps0[el]); e.f64(L.gamma[el]);
    }
    // sealed-expression probes bind the arithmetic of the condensed law
    for (double T : {298.2, 300.0, 1000.0, 1500.0, 1672.0, 1700.0, 1727.0, 1800.0, 3000.0, 5000.0}) {
        e.f64(pt_ss_h(T)); e.f64(pt_ss_cap(T)); double r, k; pt_ss_transport(T, &r, &k); e.f64(r); e.f64(k);
    }
    return e.digest();
}

}  // namespace ss304_phase_v1
