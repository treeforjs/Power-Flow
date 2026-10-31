// test_common.hpp -- deterministic test populations and ledger digests shared by
// the CPU and GPU qualification programs.
#pragma once
#include "chem_host.hpp"

namespace pfc_test {
using namespace pfc;
inline std::uint64_t splitmix(std::uint64_t& s) {
    std::uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}
inline double u01(std::uint64_t& s) { return (double(splitmix(s) >> 11) + 0.5) * 0x1p-53; }
// ------------------------------------------------------------- populations
struct Population {
    std::vector<CarrierIn> carriers;
    std::vector<CellTargets> cells;
};
inline Population make_population(const Blob& b, std::size_t n_carriers, std::size_t n_cells, std::uint64_t seed) {
    Population P;
    std::uint64_t s = seed;
    // projectile species that have groups
    std::vector<pfc_u32> proj;
    for (pfc_u32 g = 0; g < b.n_groups; ++g)
        if (std::find(proj.begin(), proj.end(), b.view.groups[g].projectile) == proj.end()) proj.push_back(b.view.groups[g].projectile);
    std::vector<pfc_u32> targets;
    for (pfc_u32 g = 0; g < b.n_groups; ++g)
        if (std::find(targets.begin(), targets.end(), b.view.groups[g].target) == targets.end()) targets.push_back(b.view.groups[g].target);
    const double kB = 8.617333262e-5;   // eV/K
    for (std::size_t c = 0; c < n_cells; ++c) {
        CellTargets ct{};
        ct.n_targets = pfc_u32(std::min<std::size_t>(targets.size(), PFC_MAX_TARGETS));
        for (pfc_u32 k = 0; k < ct.n_targets; ++k) {
            const pfc_u32 sp = targets[k];
            ct.species[k] = sp;
            const double m = b.view.species[sp].mass_eV;
            const double T = 300.0 + 3000.0 * u01(s);
            ct.s[k] = std::sqrt(kB * T / m) * PFC_C_M_S;
            ct.n[k] = (b.view.species[sp].kind == PFC_KIND_NEUTRAL ? 1e22 : 1e18) * std::exp(4 * u01(s));
            for (int a = 0; a < 3; ++a) ct.drift[k][a] = (u01(s) - 0.5) * 2000.0;
        }
        for (pfc_u32 k = ct.n_targets; k < PFC_MAX_TARGETS; ++k) ct.species[k] = 0xFFFFFFFFu;
        P.cells.push_back(ct);
    }
    for (std::size_t i = 0; i < n_carriers; ++i) {
        CarrierIn c{};
        c.id = 1 + i * 7 + (splitmix(s) % 5) + (i % 3 == 0 ? (1ull << 40) : 0);
        c.species = proj[splitmix(s) % proj.size()];
        c.cell = pfc_u32(splitmix(s) % n_cells);
        c.weight = (i % 4 == 0) ? 1e6 : std::ldexp(1.0 + u01(s), 10 + int(splitmix(s) % 20));
        const double m = b.view.species[c.species].mass_eV;
        // energies log-uniform 1e-3 .. 1e6 eV for electrons, 1e-3 .. 1e3 eV for heavy
        const double emax = c.species == 0 ? 1e6 : 1e3;
        const double E = 1e-3 * std::exp(std::log(emax / 1e-3) * u01(s));
        const double p = std::sqrt(E * (E + 2 * m));
        const V3 d = iso_dir(u01(s), u01(s));
        const V3 u = u_from_p(m, scale(d, p));
        c.u[0] = u.x;
        c.u[1] = u.y;
        c.u[2] = u.z;
        P.carriers.push_back(c);
    }
    // ids must be unique
    std::sort(P.carriers.begin(), P.carriers.end(), [](const CarrierIn& a, const CarrierIn& b) { return a.id < b.id; });
    for (std::size_t i = 1; i < P.carriers.size(); ++i)
        if (P.carriers[i].id == P.carriers[i - 1].id) P.carriers[i].id = P.carriers[i - 1].id + 1;
    // shuffle to a non-canonical input order
    for (std::size_t i = P.carriers.size(); i > 1; --i) std::swap(P.carriers[i - 1], P.carriers[splitmix(s) % i]);
    return P;
}

inline Ledger ledger_of(const Blob& b, const StepResult& r, pfc_u32 max_events) {
    Ledger L;
    L.init(b);
    for (std::size_t i = 0; i < r.outs.size(); ++i) {
        ledger_add_carrier(L, r.outs[i]);
        for (pfc_u32 k = 0; k < r.outs[i].n_events; ++k) ledger_add_event(b, L, r.events[i * max_events + k]);
    }
    ledger_finish(L);
    return L;
}
inline std::string ledger_digest(const Ledger& L) {
    Sha256 h;
    h.text("pfc-chem-v3-ledger-v1");
    for (auto* x : {&L.ke_in_projectile, &L.ke_in_target, &L.ke_out_survivor, &L.ke_out_birth, &L.ke_out_reservoir, &L.ke_out_gas, &L.q_reference, &L.defect})
        h.text(x->hex());
    for (auto& x : L.internal) h.text(x.hex());
    for (auto& x : L.created) h.text(x.hex());
    for (auto& x : L.destroyed) h.text(x.hex());
    h.text(L.charge_created.hex());
    h.text(L.charge_destroyed.hex());
    for (auto& x : L.region_exposure) h.text(x.hex());
    for (auto& x : L.region_bound) h.text(x.hex());
    for (auto& p : L.process) {
        h.pod(p.events);
        h.text(p.weight.hex());
    }
    return h.hex();
}

}  // namespace pfc_test
