// test_regression.cpp -- v1.1 with every new option off against the committed v1 owner (read in place from
// runtimes/ss304-phase-thermal-owner-v1-20261031/source): every energy word, the Fe/Cr/Ni atom words, the void
// flags and the ledger words must be bit-identical after every step; the Mn/Si/C words stay zero.
//   test_regression OUT_DIR
#include "check.hpp"
#include "/home/treefor/projects/runtimes/ss304-phase-thermal-owner-v1-20261031/source/phase_owner.hpp"
#include "phase_owner.hpp"
#include <functional>
#include <random>

namespace A = ss304_phase_v1;     // v1 (committed)
namespace B = ss304_phase_v1_1;   // v1.1, options off

struct Pair {
    A::Owner a; B::Owner b; long steps = 0, mismatches = 0;
    Pair(const A::Config& ca, const B::Config& cb) : a(ca), b(cb) {}
    bool same_state() {
        const auto sa = a.state(); const auto sb = b.state();
        if (sa.E.size() != sb.E.size()) return false;
        for (std::size_t c = 0; c < sa.E.size(); ++c) {
            if (sa.E[c].lo != sb.E[c].lo || sa.E[c].hi != sb.E[c].hi || sa.flags[c] != sb.flags[c]) return false;
            for (int el = 0; el < 3; ++el) { const auto x = sa.N[3 * c + el]; const auto y = sb.N[PQ_NEL * c + el]; if (x.lo != y.lo || x.hi != y.hi) return false; }
            for (int el = 3; el < PQ_NEL; ++el) if (!pq_is_zero(sb.N[PQ_NEL * c + el])) return false;
        }
        const auto& la = a.ledger(); const auto& lb = b.ledger();
        auto eq = [](pt_i128 x, pq_i128 y) { return x.lo == y.lo && x.hi == y.hi; };
        bool ok = eq(la.E0, lb.E0) && eq(la.joule, lb.joule) && eq(la.loads, lb.loads) && eq(la.deposit_E, lb.deposit_E) && eq(la.boundary, lb.boundary) &&
                  eq(la.evap_ref, lb.evap_ref) && eq(la.evap_kin, lb.evap_kin) && eq(la.release_E, lb.release_E) && pq_is_zero(lb.evap_int) && pq_is_zero(lb.radiation);
        for (int el = 0; el < 3; ++el) ok = ok && eq(la.N0[el], lb.N0[el]) && eq(la.evap_N[el], lb.evap_N[el]) && eq(la.dep_N[el], lb.dep_N[el]) && eq(la.rel_N[el], lb.rel_N[el]);
        return ok;
    }
};

template <class C> static void geometry(C& c, int S, double W, double H, int nx_em, int rx, std::vector<double> rows) {
    for (int s = 0; s < S; ++s) c.slices.push_back({s % 2, s / 2, -0.5 * W, s % 2 ? 0.0025 : -0.0027});
    c.width = W; c.thickness = H; c.height = 2e-4; c.nx_em = nx_em; c.ny_em = 1; c.rx = rx; c.y_edges = rows;
}
static std::vector<double> uniform_rows(double H, int n) { std::vector<double> e(n + 1); for (int j = 0; j <= n; ++j) e[j] = H * j / n; e.back() = H; return e; }

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : ".";
    // 1. evaporation, loads, deposits, hot profile, graded rows, 2 slices
    {
        const double W = 0.005, H = 0.0002;
        A::Config ca; B::Config cb;
        geometry(ca, 2, W, H, 5, 2, A::graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
        geometry(cb, 2, W, H, 5, 2, B::graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
        Pair p(ca, cb);
        auto prof = [&](int s, double x, double) { return 3500.0 + 2.8e5 * (x + 0.5 * W) + 50.0 * s; };
        p.a.initialize_temperature(prof); p.b.initialize_temperature(prof);
        std::uint64_t seq = 0;
        for (int k = 0; k < 600; ++k) {
            A::StepInput ia; B::StepInput ib; ia.dt = ib.dt = 1e-11;
            ia.em_power_W.resize(p.a.dims().n_em); for (int e = 0; e < p.a.dims().n_em; ++e) ia.em_power_W[e] = 6e6 * (1 + 0.1 * (e % 5));
            ib.em_power_W = ia.em_power_W;
            for (int s = 0; s < 2; ++s) {
                ++seq; const double y0 = ca.slices[s].y0;
                ia.loads.push_back({s, -0.5 * W + 0.0011, y0 + H, 2e-9, "t/wall", "wall", seq}); ib.loads.push_back({s, -0.5 * W + 0.0011, y0 + H, 2e-9, "t/wall", "wall", seq});
                ia.loads.push_back({s, -0.5 * W + 0.0031, y0, -1e-9, "t/em", "emission", seq}); ib.loads.push_back({s, -0.5 * W + 0.0031, y0, -1e-9, "t/em", "emission", seq});
                ia.deposits.push_back({s, -0.5 * W + 0.0021, y0 + H, 1, 1e7, 1e-12, "t/cr", seq}); ib.deposits.push_back({s, -0.5 * W + 0.0021, y0 + H, 1, 1e7, 1e-12, "t/cr", seq});
            }
            std::string ea, eb;
            try { p.a.advance(ia); } catch (const A::PhaseRefusal& e) { ea = e.kind; }
            try { p.b.advance(ib); } catch (const B::PhaseRefusal& e) { eb = e.kind; }
            ++p.steps;
            if (ea != eb || !p.same_state()) ++p.mismatches;
            if (!ea.empty()) break;
        }
        chk::check(p.mismatches == 0, "v1.1 (options off) = v1 bit for bit: evaporation, loads, deposits");
        chk::put("regression_evaporation_steps", double(p.steps));
    }
    // 2. quantized releases
    {
        const double W = 0.002, H = 0.0002;
        A::Config ca; B::Config cb;
        geometry(ca, 1, W, H, 2, 2, A::graded_rows(H, 0.5e-6, 1.4, 20e-6, 1)); geometry(cb, 1, W, H, 2, 2, B::graded_rows(H, 0.5e-6, 1.4, 20e-6, 1));
        ca.evap_mode = cb.evap_mode = 2;
        ca.patches = {{0, 1, 0.0, 0.5}, {0, 1, 0.5, 1.0}}; cb.patches = {{0, 1, 0.0, 0.5}, {0, 1, 0.5, 1.0}};
        Pair p(ca, cb);
        p.a.initialize_temperature([](int, double, double) { return 3800.0; }); p.b.initialize_temperature([](int, double, double) { return 3800.0; });
        std::mt19937_64 rng(5); std::uniform_real_distribution<double> U(0, 1);
        for (int k = 0; k < 200; ++k) {
            const auto qa = p.a.quote(1e-12); const auto qb = p.b.quote(1e-12);
            A::StepInput ia; B::StepInput ib; ia.dt = ib.dt = 1e-12;
            ia.em_power_W.assign(p.a.dims().n_em, 4e5); ib.em_power_W = ia.em_power_W;
            ia.quote_sha256 = qa.sha256; ib.quote_sha256 = qb.sha256;
            for (int pp = 0; pp < p.a.n_patches(); ++pp) for (int el = 0; el < 3; ++el) {
                if (qa.expected_atoms[pp][el] != qb.expected_atoms[pp][el]) ++p.mismatches;
                const double f = qa.expected_atoms[pp][el] / 1e6, fl = std::floor(f);
                const std::uint64_t n = std::uint64_t(fl) + (U(rng) < f - fl ? 1 : 0);
                if (!n) continue;
                const double K = double(n) * 1e6 * 2 * A::boltzmann * qa.mean_T[pp];
                ia.releases.push_back({pp, el, n, 1e6, K}); ib.releases.push_back({pp, el, n, 1e6, K});
            }
            p.a.advance(ia); p.b.advance(ib); ++p.steps;
            if (!p.same_state()) ++p.mismatches;
        }
        chk::check(p.mismatches == 0, "v1.1 (options off) = v1 bit for bit: quantized releases and quotes");
    }
    // 3. Stefan problem: constant law, prescribed temperature, many sub-steps
    {
        A::Config ca; B::Config cb;
        geometry(ca, 1, 1e-4, 6e-4, 1, 1, uniform_rows(6e-4, 150)); geometry(cb, 1, 1e-4, 6e-4, 1, 1, uniform_rows(6e-4, 150));
        ca.law.kind = 1; ca.law.c_cps = ca.law.c_cpl = 500; ca.law.c_ks = ca.law.c_kl = 20; ca.law.c_Tm = 1000; ca.law.c_dTm = 0.05; ca.law.c_L = 2.5e5; ca.law.T_min = 200; ca.law.T_max = 2e5;
        cb.law.kind = 1; cb.law.c_cps = cb.law.c_cpl = 500; cb.law.c_ks = cb.law.c_kl = 20; cb.law.c_Tm = 1000; cb.law.c_dTm = 0.05; cb.law.c_L = 2.5e5; cb.law.T_min = 200; cb.law.T_max = 2e5;
        ca.evap_mode = cb.evap_mode = 0; ca.bc_kind = cb.bc_kind = {1, 0, 0, 0}; ca.bc_T = cb.bc_T = {1200, 0, 0, 0}; ca.T0 = cb.T0 = 900;
        Pair p(ca, cb);
        for (int k = 0; k < 60; ++k) { A::StepInput ia; B::StepInput ib; ia.dt = ib.dt = 1e-5; ia.em_power_W = {0.0}; ib.em_power_W = {0.0}; p.a.advance(ia); p.b.advance(ib); ++p.steps; if (!p.same_state()) ++p.mismatches; }
        chk::check(p.mismatches == 0, "v1.1 (options off) = v1 bit for bit: Stefan problem with sub-steps");
    }
    // 4. surface-cell merge
    {
        const double H = 0.0002;
        std::vector<double> rows{0.0, 20e-9, 40e-9, 80e-9, 160e-9, 320e-9, 640e-9, 1.28e-6, 2.56e-6, 5.12e-6, 10.24e-6, 20.48e-6, 40.96e-6, 81.92e-6, H};
        A::Config ca; B::Config cb;
        geometry(ca, 1, 1e-4, H, 1, 1, rows); geometry(cb, 1, 1e-4, H, 1, 1, rows);
        ca.evap_face = cb.evap_face = {1, 0, 0, 0};
        Pair p(ca, cb);
        p.a.initialize_temperature([](int, double, double) { return 4990.0; }); p.b.initialize_temperature([](int, double, double) { return 4990.0; });
        int merged = 0;
        for (int k = 0; k < 150000 && !merged; ++k) {
            A::StepInput ia; B::StepInput ib; ia.dt = ib.dt = 1e-12; ia.em_power_W = {0.0}; ib.em_power_W = {0.0};
            const auto ra = p.a.advance(ia); p.b.advance(ib); ++p.steps;
            merged = ra.merges;
            if (merged || k % 5000 == 0) if (!p.same_state()) ++p.mismatches;
        }
        chk::check(merged && p.mismatches == 0, "v1.1 (options off) = v1 bit for bit through a surface-cell merge");
        chk::put("regression_merge_steps", double(p.steps));
    }
    return chk::finish("test_regression", out);
}
