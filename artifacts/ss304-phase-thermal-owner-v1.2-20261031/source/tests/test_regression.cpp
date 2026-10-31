// test_regression.cpp -- v1.2 with every new option off against the committed v1.1 owner (read in place from
// runtimes/ss304-phase-thermal-owner-v1.1-20261031/source): the energy, atom, flag and face-mass words, the ledger
// words, the step records, the cell and face-column Merkle roots and the patch reports must be bit-identical after
// every step, with the v1.1 options off (v1 behaviour) and on (Config::enable_v1_1 in both owners).
//   test_regression OUT_DIR [suite ...]   suites: v1 quantized v11 merge
#include "check.hpp"
#include "/home/treefor/projects/runtimes/ss304-phase-thermal-owner-v1.1-20261031/source/phase_owner.hpp"
#include "phase_owner.hpp"
#include <cstring>
#include <functional>
#include <random>
#include <set>

namespace A = ss304_phase_v1_1;   // v1.1 (committed)
namespace B = ss304_phase_v1_2;   // v1.2, new options off

static bool eqw(pq_i128 x, px_i128 y) { return x.lo == y.lo && x.hi == y.hi; }

struct Pair {
    A::Owner a; B::Owner b; long steps = 0, mismatches = 0, record_mismatches = 0;
    std::string stop_a, stop_b;
    Pair(const A::Config& ca, const B::Config& cb) : a(ca), b(cb) {}
    bool same_state() {
        const auto sa = a.state(); const auto sb = b.state();
        if (sa.E.size() != sb.E.size() || sa.fcm.size() != sb.fcm.size()) return false;
        for (std::size_t c = 0; c < sa.E.size(); ++c) {
            if (!eqw(sa.E[c], sb.E[c]) || sa.flags[c] != sb.flags[c]) return false;
            for (int el = 0; el < PX_NEL; ++el) if (!eqw(sa.N[PQ_NEL * c + el], sb.N[PX_NEL * c + el])) return false;
        }
        for (std::size_t q = 0; q < sa.fcm.size(); ++q) if (!eqw(sa.fcm[q], sb.fcm[q])) return false;
        for (double v : sb.fcV) if (v != 0.0) return false;
        for (auto w : sb.fcK) if (!px_is_zero(w)) return false;
        const auto& la = a.ledger(); const auto& lb = b.ledger();
        bool ok = eqw(la.E0, lb.E0) && eqw(la.joule, lb.joule) && eqw(la.loads, lb.loads) && eqw(la.deposit_E, lb.deposit_E) && eqw(la.boundary, lb.boundary) &&
                  eqw(la.evap_ref, lb.evap_ref) && eqw(la.evap_kin, lb.evap_kin) && eqw(la.evap_int, lb.evap_int) && eqw(la.radiation, lb.radiation) &&
                  eqw(la.release_E, lb.release_E);
        for (int el = 0; el < PX_NEL; ++el) ok = ok && eqw(la.N0[el], lb.N0[el]) && eqw(la.evap_N[el], lb.evap_N[el]) && eqw(la.dep_N[el], lb.dep_N[el]) && eqw(la.rel_N[el], lb.rel_N[el]);
        for (auto w : {lb.K0, lb.lorentz, lb.surface, lb.dissipation, lb.eject_E, lb.eject_K}) ok = ok && px_is_zero(w);
        const auto ra = a.state_roots(); const auto rb = b.state_roots();
        ok = ok && std::memcmp(ra.cells.data(), rb.cells.data(), 32) == 0 && std::memcmp(ra.fcm.data(), rb.fcm.data(), 32) == 0 && !rb.has_flow;
        return ok;
    }
    bool same_reports() {
        const auto pa = a.patch_reports(); const auto pb = b.patch_reports();
        if (pa.size() != pb.size()) return false;
        for (std::size_t i = 0; i < pa.size(); ++i) {
            const auto& x = pa[i]; const auto& y = pb[i];
            const double va[] = {x.area_m2, x.mean_T, x.max_T, x.liquid_fraction, x.melt_depth_m, x.surface_density, x.recession_m, x.displacement_mean_m,
                                 x.displacement_max_abs_m, x.recoil_pressure_Pa, x.melt_speed_m_s, x.released_internal_J, x.reynolds_max, x.cumulative_evaporated_kg};
            const double vb[] = {y.area_m2, y.mean_T, y.max_T, y.liquid_fraction, y.melt_depth_m, y.surface_density, y.recession_m, y.displacement_mean_m,
                                 y.displacement_max_abs_m, y.recoil_pressure_Pa, y.melt_speed_m_s, y.released_internal_J, y.reynolds_max, y.cumulative_evaporated_kg};
            if (std::memcmp(va, vb, sizeof va) != 0 || x.id != y.id) return false;
            for (int el = 0; el < PX_NEL; ++el)
                if (x.flux_atoms_m2_s[el] != y.flux_atoms_m2_s[el] || x.expected_atoms[el] != y.expected_atoms[el] || x.surface_atom_fraction[el] != y.surface_atom_fraction[el]) return false;
        }
        const auto sa = a.sigma_em(); const auto sb = b.sigma_em();
        return sa.size() == sb.size() && std::memcmp(sa.data(), sb.data(), 8 * sa.size()) == 0;
    }
    static bool same_record(const A::StepRecord& x, const B::StepRecord& y) {
        bool ok = x.n_sub == y.n_sub && x.merges == y.merges && x.max_T == y.max_T && x.max_displacement_fraction == y.max_displacement_fraction &&
                  x.max_reynolds == y.max_reynolds && x.phase_explosion_margin == y.phase_explosion_margin && x.tau == y.tau;
        for (auto [p, q] : {std::pair{x.dE, y.dE}, {x.joule, y.joule}, {x.loads, y.loads}, {x.deposit_E, y.deposit_E}, {x.boundary, y.boundary},
                            {x.evap_ref, y.evap_ref}, {x.evap_kin, y.evap_kin}, {x.evap_int, y.evap_int}, {x.radiation, y.radiation}, {x.release_E, y.release_E}}) ok = ok && eqw(p, q);
        for (int el = 0; el < PX_NEL; ++el) ok = ok && eqw(x.dN[el], y.dN[el]) && eqw(x.evap_N[el], y.evap_N[el]) && eqw(x.rel_N[el], y.rel_N[el]);
        return ok && px_is_zero(y.dK) && px_is_zero(y.lorentz) && px_is_zero(y.eject_E);
    }
    // one step in both owners; refusals must agree in kind
    template <class FA, class FB> bool step(FA&& ia, FB&& ib) {
        std::string ea, eb; A::StepRecord ra; B::StepRecord rb;
        try { ra = a.advance(ia); } catch (const A::PhaseRefusal& e) { ea = e.kind; stop_a = e.what(); }
        try { rb = b.advance(ib); } catch (const B::PhaseRefusal& e) { eb = e.kind; stop_b = e.what(); }
        ++steps;
        if (ea != eb) ++mismatches;
        else if (ea.empty() && !same_record(ra, rb)) ++record_mismatches;
        if (!same_state()) ++mismatches;
        return ea.empty();
    }
};

template <class C> static void geometry(C& c, int S, double W, double H, int nx_em, int rx, std::vector<double> rows) {
    for (int s = 0; s < S; ++s) c.slices.push_back({s % 2, s / 2, -0.5 * W, s % 2 ? 0.0025 : -0.0027});
    c.width = W; c.thickness = H; c.height = 2e-4; c.nx_em = nx_em; c.ny_em = 1; c.rx = rx; c.y_edges = rows;
}

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : ".";
    std::set<std::string> want; for (int i = 2; i < argc; ++i) want.insert(argv[i]);
    auto on = [&](const char* n) { return want.empty() || want.count(n); };
    long total_steps = 0;
    // 1. options off: evaporation, loads, deposits, hot profile, graded rows, 2 slices
    if (on("v1")) {
        const double W = 0.005, H = 0.0002;
        A::Config ca; B::Config cb;
        geometry(ca, 2, W, H, 5, 2, A::graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
        geometry(cb, 2, W, H, 5, 2, B::graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
        Pair p(ca, cb);
        auto prof = [&](int s, double x, double) { return 3500.0 + 2.8e5 * (x + 0.5 * W) + 50.0 * s; };
        p.a.initialize_temperature(prof); p.b.initialize_temperature(prof);
        std::uint64_t seq = 0;
        for (int k = 0; k < 400; ++k) {
            A::StepInput ia; B::StepInput ib; ia.dt = ib.dt = 1e-11;
            ia.em_power_W.resize(p.a.dims().n_em); for (int e = 0; e < p.a.dims().n_em; ++e) ia.em_power_W[e] = 6e6 * (1 + 0.1 * (e % 5));
            ib.em_power_W = ia.em_power_W;
            for (int s = 0; s < 2; ++s) {
                ++seq; const double y0 = ca.slices[s].y0;
                ia.loads.push_back({s, -0.5 * W + 0.0011, y0 + H, 2e-9, "t/wall", "wall", seq}); ib.loads.push_back({s, -0.5 * W + 0.0011, y0 + H, 2e-9, "t/wall", "wall", seq});
                ia.loads.push_back({s, -0.5 * W + 0.0031, y0, -1e-9, "t/em", "emission", seq}); ib.loads.push_back({s, -0.5 * W + 0.0031, y0, -1e-9, "t/em", "emission", seq});
                ia.deposits.push_back({s, -0.5 * W + 0.0021, y0 + H, 1, 1e7, 1e-12, "t/cr", seq}); ib.deposits.push_back({s, -0.5 * W + 0.0021, y0 + H, 1, 1e7, 1e-12, "t/cr", seq});
            }
            if (!p.step(ia, ib)) break;
        }
        const bool rep = p.same_reports();
        chk::check(p.mismatches == 0 && p.record_mismatches == 0 && rep, "v1.2 (options off) = v1.1 bit for bit: evaporation, loads, deposits (state, ledger, records, roots, reports)");
        chk::put("regression_v1_steps", double(p.steps)); total_steps += p.steps;
        chk::put_text("regression_v1_stop_v1_1", p.stop_a); chk::put_text("regression_v1_stop_v1_2", p.stop_b);
    }
    // 2. options off: quantized releases and quotes
    if (on("quantized")) {
        const double W = 0.002, H = 0.0002;
        A::Config ca; B::Config cb;
        geometry(ca, 1, W, H, 2, 2, A::graded_rows(H, 0.5e-6, 1.4, 20e-6, 1)); geometry(cb, 1, W, H, 2, 2, B::graded_rows(H, 0.5e-6, 1.4, 20e-6, 1));
        ca.evap_mode = cb.evap_mode = 2;
        ca.patches = {{0, 1, 0.0, 0.5}, {0, 1, 0.5, 1.0}}; cb.patches = {{0, 1, 0.0, 0.5}, {0, 1, 0.5, 1.0}};
        Pair p(ca, cb);
        p.a.initialize_temperature([](int, double, double) { return 3800.0; }); p.b.initialize_temperature([](int, double, double) { return 3800.0; });
        std::mt19937_64 rng(5); std::uniform_real_distribution<double> U(0, 1);
        long quote_mismatch = 0;
        for (int k = 0; k < 200; ++k) {
            const auto qa = p.a.quote(1e-12); const auto qb = p.b.quote(1e-12);
            A::StepInput ia; B::StepInput ib; ia.dt = ib.dt = 1e-12;
            ia.em_power_W.assign(p.a.dims().n_em, 4e5); ib.em_power_W = ia.em_power_W;
            ia.quote_sha256 = qa.sha256; ib.quote_sha256 = qb.sha256;
            for (int pp = 0; pp < p.a.n_patches(); ++pp) for (int el = 0; el < 3; ++el) {
                if (qa.expected_atoms[pp][el] != qb.expected_atoms[pp][el]) ++quote_mismatch;
                const double f = qa.expected_atoms[pp][el] / 1e6, fl = std::floor(f);
                const std::uint64_t n = std::uint64_t(fl) + (U(rng) < f - fl ? 1 : 0);
                if (!n) continue;
                const double K = double(n) * 1e6 * 2 * A::boltzmann * qa.mean_T[pp];
                ia.releases.push_back({pp, el, n, 1e6, K}); ib.releases.push_back({pp, el, n, 1e6, K});
            }
            if (!p.step(ia, ib)) break;
        }
        chk::check(p.mismatches == 0 && p.record_mismatches == 0 && quote_mismatch == 0 && p.same_reports(), "v1.2 (options off) = v1.1 bit for bit: quantized releases and quotes");
        total_steps += p.steps;
    }
    // 3. v1.1 options on (both owners): melt flow (lubrication) driven by J x B, Mn, activities, Knudsen layer with back
    //    pressures, radiation, internal energy, loads and deposits
    if (on("v11")) {
        const double W = 0.005, H = 0.0002;
        A::Config ca; B::Config cb;
        geometry(ca, 2, W, H, 5, 2, A::graded_rows(H, 0.5e-6, 1.35, 10e-6, 1)); geometry(cb, 2, W, H, 5, 2, B::graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
        A::Config::enable_v1_1(ca); B::Config::enable_v1_1(cb);
        Pair p(ca, cb);
        auto prof = [&](int s, double x, double) { return 2400.0 + 2.0e5 * (x + 0.5 * W) + 30.0 * s; };
        p.a.initialize_temperature(prof); p.b.initialize_temperature(prof);
        const int nem = p.a.dims().n_em;
        std::vector<double> JB(6 * nem);
        for (int e = 0; e < nem; ++e) { const double j = 5e8 * (1 + 0.02 * (e % 5)); JB[6 * e + 2] = j; JB[6 * e + 3] = 20.0; JB[6 * e + 4] = -0.05; }
        std::vector<double> pb(p.a.n_patches(), 0.0);
        for (int q = 0; q < p.a.n_patches(); ++q) pb[q] = (q % 3) * 2e3;
        std::uint64_t seq = 0;
        for (int k = 0; k < 600; ++k) {
            A::StepInput ia; B::StepInput ib; ia.dt = ib.dt = 1e-11;
            ia.em_power_W.assign(nem, 3e6); ib.em_power_W = ia.em_power_W; ia.em_JB = ib.em_JB = JB; ia.back_pressure_Pa = ib.back_pressure_Pa = pb;
            for (int s = 0; s < 2; ++s) {
                ++seq; const double y1 = ca.slices[s].y0 + H;
                ia.loads.push_back({s, -0.5 * W + 0.0011, y1, 2e-9, "t/wall", "wall", seq}); ib.loads.push_back({s, -0.5 * W + 0.0011, y1, 2e-9, "t/wall", "wall", seq});
                ia.deposits.push_back({s, -0.5 * W + 0.0021, y1, 3, 1e6, 1e-13, "t/mn-return", seq}); ib.deposits.push_back({s, -0.5 * W + 0.0021, y1, 3, 1e6, 1e-13, "t/mn-return", seq});
            }
            if (!p.step(ia, ib)) break;
            if (k % 100 == 99 && !p.same_reports()) ++p.mismatches;
        }
        chk::check(p.mismatches == 0 && p.record_mismatches == 0 && p.steps >= 300, "v1.2 with enable_v1_1 = v1.1 with enable_v1_1 bit for bit: lubrication flow, Mn, Knudsen, radiation");
        chk::put("regression_v11_options_steps", double(p.steps)); total_steps += p.steps;
        chk::put("regression_v11_mismatches", double(p.mismatches)); chk::put("regression_v11_record_mismatches", double(p.record_mismatches));
        chk::put_text("regression_v11_stop_v1_1", p.stop_a); chk::put_text("regression_v11_stop_v1_2", p.stop_b);
        const auto sa = p.a.state_sha256(); const auto sb = p.b.state_sha256();
        chk::put_text("regression_v11_state_sha256_v1_1", sa); chk::put_text("regression_v11_state_sha256_v1_2", sb);
    }
    // 4. options off: surface-cell merge
    if (on("merge")) {
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
            const auto ra = p.a.advance(ia); const auto rb = p.b.advance(ib); ++p.steps;
            if (!Pair::same_record(ra, rb)) ++p.record_mismatches;
            merged = ra.merges;
            if (merged || k % 5000 == 0) if (!p.same_state()) ++p.mismatches;
        }
        chk::check(merged && p.mismatches == 0 && p.record_mismatches == 0 && p.same_reports(), "v1.2 (options off) = v1.1 bit for bit through a surface-cell merge");
        chk::put("regression_merge_steps", double(p.steps)); total_steps += p.steps;
    }
    chk::put("regression_total_steps", double(total_steps));
    std::string name = "test_regression"; for (const auto& w : want) name += "-" + w;
    return chk::finish(name, out);
}
