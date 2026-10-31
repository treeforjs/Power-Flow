// test_gpu.cpp -- v1.2: CPU reference versus NVRTC GPU engine, bit identity of every accepted state (block-Merkle
// digest computed on the host for the CPU engine and on the device for the GPU engine, with the layer-flow section),
// step record, patch report, quote and sigma vector; device Merkle leaves against a host recomputation;
// caller-provided session; sub-cycling with host and device accumulation; timing per field step on the Hall foil.
//   test_gpu OUT_DIR NVRTC_LIBRARY GPU_UUID [suite ...]
//   suites: v1 v11 quantized leveling merkle borrowed hall v12 dam eject accdev timing
#include "check.hpp"
#include "phase_owner.hpp"
#include "gpu_engine.hpp"
#include "common/hall_gpu_common.hpp"
#include <chrono>
#include <cstring>
#include <functional>
#include <random>
#include <set>

using namespace ss304_phase_v1_2;
using Clock = std::chrono::steady_clock;
static GpuOptions gopt;

static bool same_doubles(const std::vector<double>& a, const std::vector<double>& b) { return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * 8) == 0; }
static bool same_reports(const std::vector<PatchReport>& a, const std::vector<PatchReport>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        auto pack = [](const PatchReport& p) {
            std::vector<double> v{p.area_m2, p.mean_T, p.max_T, p.liquid_fraction, p.melt_depth_m, p.surface_density, p.recession_m, p.displacement_mean_m,
                                  p.displacement_max_abs_m, p.recoil_pressure_Pa, p.melt_speed_m_s, p.released_internal_J, p.reynolds_max,
                                  p.crater_depth_m, p.layer_ux_m_s, p.layer_uz_m_s, p.layer_thickness_m, p.cumulative_ejected_J, p.cumulative_ejected_kinetic_J,
                                  p.cumulative_lorentz_work_J};
            for (int el = 0; el < PX_NEL; ++el) { v.push_back(p.flux_atoms_m2_s[el]); v.push_back(p.expected_atoms[el]); v.push_back(p.surface_atom_fraction[el]); v.push_back(p.cumulative_evaporated_atoms[el]); v.push_back(p.cumulative_ejected_atoms[el]); }
            return v;
        };
        if (!same_doubles(pack(a[i]), pack(b[i]))) return false;
    }
    return true;
}
static bool same_record(const StepRecord& a, const StepRecord& b) {
    auto eq = [](px_i128 x, px_i128 y) { return px_eq(x, y); };
    bool ok = a.n_sub == b.n_sub && a.merges == b.merges && std::memcmp(&a.max_T, &b.max_T, 8) == 0 && a.inputs_sha256 == b.inputs_sha256 &&
              std::memcmp(&a.max_displacement_fraction, &b.max_displacement_fraction, 8) == 0 && std::memcmp(&a.max_reynolds, &b.max_reynolds, 8) == 0;
    for (auto [x, y] : {std::pair{a.dE, b.dE}, {a.joule, b.joule}, {a.loads, b.loads}, {a.deposit_E, b.deposit_E}, {a.boundary, b.boundary}, {a.evap_ref, b.evap_ref},
                        {a.evap_kin, b.evap_kin}, {a.evap_int, b.evap_int}, {a.radiation, b.radiation}, {a.release_E, b.release_E}}) ok = ok && eq(x, y);
    for (int el = 0; el < PX_NEL; ++el) ok = ok && eq(a.dN[el], b.dN[el]) && eq(a.evap_N[el], b.evap_N[el]) && eq(a.rel_N[el], b.rel_N[el]) && eq(a.dep_N[el], b.dep_N[el]);
    // v1.2
    for (auto [x, y] : {std::pair{a.dK, b.dK}, {a.lorentz, b.lorentz}, {a.surface, b.surface}, {a.dissipation, b.dissipation}, {a.eject_E, b.eject_E}, {a.eject_K, b.eject_K}}) ok = ok && eq(x, y);
    for (int el = 0; el < PX_NEL; ++el) ok = ok && eq(a.eject_N[el], b.eject_N[el]);
    ok = ok && a.cfl_retries == b.cfl_retries && a.window_steps == b.window_steps && std::memcmp(&a.sigma_change_max, &b.sigma_change_max, 8) == 0 &&
         std::memcmp(&a.window_heating_K, &b.window_heating_K, 8) == 0 && std::memcmp(&a.max_signal_speed, &b.max_signal_speed, 8) == 0 &&
         std::memcmp(&a.max_crater_fraction, &b.max_crater_fraction, 8) == 0 && std::memcmp(&a.joule_requested_J, &b.joule_requested_J, 8) == 0;
    for (std::size_t p = 0; p < a.patch_lorentz.size() && p < b.patch_lorentz.size(); ++p) ok = ok && eq(a.patch_lorentz[p], b.patch_lorentz[p]) && eq(a.patch_eject_E[p], b.patch_eject_E[p]);
    return ok;
}
struct Pair {
    Owner cpu, gpu;
    long steps = 0, mismatches = 0;
    Pair(const Config& c, const GpuOptions& g) : cpu(c), gpu(c, std::make_unique<GpuEngine>(g)) {}
    template <class F> void init(F f) { cpu.initialize_temperature(f); gpu.initialize_temperature(f); }
    bool compare(bool full) {
        bool ok = cpu.state_sha256() == gpu.state_sha256();
        if (full) ok = ok && same_reports(cpu.patch_reports(), gpu.patch_reports()) && same_doubles(cpu.sigma_em(), gpu.sigma_em()) &&
                       same_doubles(cpu.face_displacements(), gpu.face_displacements());
        return ok;
    }
    void accumulate(const std::vector<double>& P, const std::vector<double>& JB, double dt) { cpu.accumulate(P, JB, dt); gpu.accumulate(P, JB, dt); }
    void step(const StepInput& in, bool full) {
        StepRecord a, b; std::string ea, eb;
        try { a = cpu.advance(in); } catch (const PhaseRefusal& e) { ea = e.kind; }
        try { b = gpu.advance(in); } catch (const PhaseRefusal& e) { eb = e.kind; }
        ++steps;
        if (ea != eb || (ea.empty() && !same_record(a, b)) || !compare(full)) ++mismatches;
        if (!ea.empty()) throw PhaseRefusal(ea, "pair", "pair", 0, 0, "both refused");
    }
};
static Config cfg(int S, double W, double H, int nx_em, int rx, std::vector<double> rows) {
    Config c;
    for (int s = 0; s < S; ++s) c.slices.push_back({s % 2, s / 2, -0.5 * W, s % 2 ? 0.0025 : -0.0027});
    c.width = W; c.thickness = H; c.height = 2e-4; c.nx_em = nx_em; c.ny_em = 1; c.rx = rx; c.y_edges = std::move(rows);
    return c;
}

// v1 scenario in v1 mode (as v1's parity_evaporation)
static void v1_evaporation() {
    const double W = 0.005, H = 0.0002;
    Config c = cfg(2, W, H, 5, 2, graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
    Pair p(c, gopt);
    p.init([&](int s, double x, double) { return 3500.0 + 2.8e5 * (x + 0.5 * W) + 50.0 * s; });
    std::uint64_t seq = 0;
    for (int k = 0; k < 200; ++k) {
        StepInput in; in.dt = 1e-11; in.em_power_W.resize(p.cpu.dims().n_em);
        for (int e = 0; e < p.cpu.dims().n_em; ++e) in.em_power_W[e] = 6e6 * (1 + 0.1 * (e % 5));
        for (int s = 0; s < 2; ++s) {
            ++seq; const double y0 = c.slices[s].y0;
            in.loads.push_back({s, -0.5 * W + 0.0011, y0 + H, 2e-9, "t/wall", "wall", seq});
            in.deposits.push_back({s, -0.5 * W + 0.0021, y0 + H, 1, 1e7, 1e-12, "t/cr-return", seq});
        }
        try { p.step(in, k % 25 == 0); } catch (const PhaseRefusal&) { break; }
    }
    chk::check(p.mismatches == 0, "GPU = CPU bit for bit, v1 mode");
    chk::put("parity_v1_steps", double(p.steps));
}
// all v1.1 options: Mn, activities, Knudsen with back pressures, radiation, internal energy, melt flow from J x B
static void v11_ledger() {
    const double W = 0.005, H = 0.0002;
    Config c = cfg(2, W, H, 5, 2, graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
    Config::enable_v1_1(c);
    Pair p(c, gopt);
    p.init([&](int s, double x, double) { return 2400.0 + 2.0e5 * (x + 0.5 * W) + 30.0 * s; });
    const int nem = p.cpu.dims().n_em;
    std::vector<double> JB(6 * nem);
    for (int e = 0; e < nem; ++e) { JB[6 * e + 2] = 5e8 * (1 + 0.02 * (e % 5)); JB[6 * e + 3] = 20.0; JB[6 * e + 4] = -0.05; }
    std::vector<double> pb(p.cpu.n_patches()); for (int q = 0; q < p.cpu.n_patches(); ++q) pb[q] = (q % 3) * 2e3;
    std::uint64_t seq = 0;
    for (int k = 0; k < 300; ++k) {
        StepInput in; in.dt = 1e-11; in.em_power_W.assign(nem, 3e6); in.em_JB = JB; in.back_pressure_Pa = pb;
        for (int s = 0; s < 2; ++s) { ++seq; const double y0 = c.slices[s].y0; in.deposits.push_back({s, -0.5 * W + 0.0021, y0 + H, 3, 1e6, 1e-13, "t/mn", seq}); }
        try { p.step(in, k % 20 == 0); } catch (const PhaseRefusal& e) { chk::put_text("parity_v11_stop", e.what()); break; }
    }
    chk::check(p.mismatches == 0 && p.steps == 300, "GPU = CPU bit for bit with every v1.1 option (melt flow, Mn, Knudsen, radiation, internal energy)");
    chk::check(p.gpu.ledger_closes(), "exact ledgers on the GPU engine (v1.1)");
    chk::put("parity_v11_steps", double(p.steps));
}
static void v11_quantized() {
    const double W = 0.002, H = 0.0002;
    Config c = cfg(1, W, H, 2, 2, graded_rows(H, 0.5e-6, 1.4, 20e-6, 1));
    Config::enable_v1_1(c); c.flow = false; c.evap_mode = 2;
    c.patches = {{0, 1, 0.0, 0.5}, {0, 1, 0.5, 1.0}, {0, 0, 0.0, 0.3}, {0, 0, 0.3, 1.0}};
    Pair p(c, gopt);
    p.init([](int, double, double) { return 3600.0; });
    std::mt19937_64 rng(9); std::uniform_real_distribution<double> U(0, 1);
    for (int k = 0; k < 120; ++k) {
        const Quote qc = p.cpu.quote(1e-12), qg = p.gpu.quote(1e-12);
        if (qc.sha256 != qg.sha256) ++p.mismatches;
        StepInput in; in.dt = 1e-12; in.em_power_W.assign(p.cpu.dims().n_em, 4e5); in.quote_sha256 = qc.sha256;
        for (int pp = 0; pp < p.cpu.n_patches(); ++pp) for (int el = 0; el < PX_NEL; ++el) {
            const double f = qc.expected_atoms[pp][el] / 1e6, fl = std::floor(f);
            const std::uint64_t n = std::uint64_t(fl) + (U(rng) < f - fl ? 1 : 0);
            if (n) in.releases.push_back({pp, el, n, 1e6, double(n) * 1e6 * 2 * boltzmann * qc.mean_T[pp]});
        }
        p.step(in, k % 10 == 0);
    }
    chk::check(p.mismatches == 0, "GPU = CPU bit for bit: quantized releases with internal energy and activities");
}
static void v11_leveling() {
    const int nx = 20; const double W = 2e-3, H = 1e-4, delta = 2e-5;
    std::vector<double> rows; for (int j = 0; j <= 4; ++j) rows.push_back((H - delta) * j / 4); for (int j = 1; j <= 4; ++j) rows.push_back(H - delta + delta * j / 4); rows.back() = H;
    Config c = cfg(1, W, H, nx, 1, rows);
    LawOptions o; o.kind = 1; o.c_cps = o.c_cpl = 500; o.c_ks = o.c_kl = 1e-12; o.c_Tm = 1000; o.c_dTm = 1; o.T_min = 200; o.T_max = 5000; o.c_rho = 7000; o.c_mu = 6e-3; o.c_sigma = 1.8;
    c.law = o; c.rho_ref = 7000; c.evap_mode = 0; c.flow = true; c.T0 = 800; c.slices[0].y0 = 0.0;
    Pair p(c, gopt);
    p.init([&](int, double, double y) { return y > H - delta ? 1200.0 : 800.0; });
    auto disp = [&](int, int f, double x) { return f == 1 ? 1e-8 * std::cos(M_PI * (x + 0.5 * W) / W) : 0.0; };
    p.cpu.initialize_displacement(disp); p.gpu.initialize_displacement(disp);
    for (int k = 0; k < 10; ++k) { StepInput in; in.dt = 1e-4; in.em_power_W.assign(nx, 0.0); p.step(in, true); }
    chk::check(p.mismatches == 0, "GPU = CPU bit for bit: capillary melt flow with sub-steps");
    chk::put("parity_leveling_substeps_per_step", double(p.cpu.substeps_for(1e-4)));
}
static void merkle() {
    const double W = 0.005, H = 0.0002;
    Config c = cfg(4, W, H, 25, 1, graded_rows(H, 0.25e-6, 1.3, 10e-6, 1));
    Config::enable_v1_1(c);
    Owner g(c, std::make_unique<GpuEngine>(gopt)), h(c);
    StepInput in; in.dt = 1e-14; in.em_power_W.assign(g.dims().n_em, 1e5);
    g.advance(in); h.advance(in);
    chk::check(g.verify_state_leaves(), "device Merkle leaves equal a host recomputation from the downloaded state");
    chk::check(g.state_sha256() == h.state_sha256(), "device Merkle digest = host digest of the CPU engine");
    chk::put("merkle_cell_leaves", double(g.state_roots().cell_leaves));
}
static void borrowed() {
    hall_gpu::SessionOptions so; so.require_uuid = gopt.require_uuid; so.nvrtc_library = gopt.nvrtc_library;
    hall_gpu::Session session(so);
    GpuOptions b = gopt; b.session = &session;
    const double W = 0.002, H = 0.0002;
    Config c = cfg(2, W, H, 4, 1, graded_rows(H, 0.5e-6, 1.4, 20e-6, 1));
    Config::enable_v1_1(c);
    bool flag = false;
    {
        auto eng = std::make_unique<GpuEngine>(b); flag = eng->borrowed_session();
        Owner g(c, std::move(eng)), h(c);
        g.initialize_temperature([](int, double, double) { return 3000.0; }); h.initialize_temperature([](int, double, double) { return 3000.0; });
        long mism = 0;
        for (int k = 0; k < 30; ++k) { StepInput in; in.dt = 1e-11; in.em_power_W.assign(g.dims().n_em, 2e6); g.advance(in); h.advance(in); if (g.state_sha256() != h.state_sha256()) ++mism; }
        chk::check(mism == 0, "engine on a caller-provided session: GPU = CPU bit for bit");
    }
    session.synchronize("after borrowed engine");
    chk::check(flag && !session.failed(), "caller-provided session used and still healthy after the engine is destroyed");
    const auto rep = session.close();
    chk::check(rep.released && !rep.quarantined, "caller's session closes cleanly");
}
static void hall_timing() {
    Config c;
    for (int k = 0; k < 125; ++k) for (int body = 0; body < 2; ++body) c.slices.push_back({body, k, -0.0025, body ? 0.0025 : -0.0027});
    c.width = 0.005; c.thickness = 0.0002; c.height = 0.0002; c.nx_em = 25; c.ny_em = 1; c.rx = 1;
    c.y_edges = graded_rows(0.0002, 0.25e-6, 1.3, 10e-6, 1);
    Config::enable_v1_1(c);
    Owner g(c, std::make_unique<GpuEngine>(gopt)), h(c);
    const int nem = g.dims().n_em;
    std::vector<double> JB(6 * nem); for (int e = 0; e < nem; ++e) { JB[6 * e + 2] = 1e8; JB[6 * e + 3] = 5.0; JB[6 * e + 4] = -0.01; }
    StepInput in; in.dt = 1e-14; in.em_power_W.assign(nem, 2e5); in.em_JB = JB;
    for (int k = 0; k < 3; ++k) { g.advance(in); h.advance(in); }
    chk::check(g.state_sha256() == h.state_sha256(), "full Hall geometry, v1.1 options: GPU = CPU after 3 steps");
    auto t0 = Clock::now();
    for (int k = 0; k < 3; ++k) h.advance(in);
    const double cpu_ms = 1e3 * std::chrono::duration<double>(Clock::now() - t0).count() / 3;
    t0 = Clock::now();
    const int n = 200;
    for (int k = 0; k < n; ++k) g.advance(in);
    const double gpu_ms = 1e3 * std::chrono::duration<double>(Clock::now() - t0).count() / n;
    t0 = Clock::now();
    const std::string dg = g.state_sha256();
    const double dig_ms = 1e3 * std::chrono::duration<double>(Clock::now() - t0).count();
    chk::put("hall_full_cells", double(g.dims().n_cells)); chk::put("hall_v11_cpu_ms_per_step", cpu_ms);
    chk::put("hall_v11_gpu_ms_per_step", gpu_ms); chk::put("hall_v11_device_digest_ms", dig_ms);
}

// ------------------------------------------------------------------------------------------------ v1.2
// every v1.2 option (inertial layer along x and z, ejection, diagnostics) with sub-cycling K = 4 (host accumulation),
// four slices (two bodies, z neighbours), J x B, back pressures, loads and deposits
static void v12_parity() {
    const double W = 0.005, H = 0.0002;
    Config c = cfg(4, W, H, 5, 2, graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
    Config::enable_v1_2(c);
    Pair p(c, gopt);
    p.init([&](int s, double x, double) { return 2600.0 + 2.0e5 * (x + 0.5 * W) + 30.0 * s; });
    const int nem = p.cpu.dims().n_em;
    std::vector<double> JB(6 * nem);
    for (int e = 0; e < nem; ++e) { JB[6 * e + 2] = 5e8 * (1 + 0.02 * (e % 5)); JB[6 * e + 3] = 20.0; JB[6 * e + 4] = -0.05; }
    std::vector<double> pb(p.cpu.n_patches()); for (int q = 0; q < p.cpu.n_patches(); ++q) pb[q] = (q % 3) * 2e3;
    std::uint64_t seq = 0;
    for (int k = 0; k < 200; ++k) {
        for (int q = 0; q < 4; ++q) { std::vector<double> P(nem); for (int e = 0; e < nem; ++e) P[e] = 3e6 * (1 + 0.2 * std::sin(0.3 * (4 * k + q) + e)); p.accumulate(P, JB, 2.5e-12); }
        StepInput in; in.use_accumulated = true; in.back_pressure_Pa = pb;
        for (int s = 0; s < 4; ++s) { ++seq; const double y0 = c.slices[s].y0; in.deposits.push_back({s, -0.5 * W + 0.0021, y0 + H, 3, 1e6, 1e-13, "t/mn", seq}); }
        try { p.step(in, k % 20 == 0); } catch (const PhaseRefusal& e) { chk::put_text("parity_v12_stop", e.what()); break; }
    }
    chk::check(p.mismatches == 0 && p.steps == 200, "GPU = CPU bit for bit with every v1.2 option (inertial layer x/z, ejection, sub-cycling K = 4)");
    chk::check(p.gpu.ledger_closes(), "exact ledgers on the GPU engine (v1.2)");
    chk::check(p.gpu.verify_state_leaves(), "device Merkle leaves with the layer-flow section = host recomputation");
    chk::put("parity_v12_steps", double(p.steps));
    chk::put("parity_v12_kinetic_J", px_to_double(p.gpu.kinetic_word(), PX_E_INV));
}
// shallow-water dam break (inviscid layer, CFL sub-steps with halvings) and edge plus departure ejection
static Config sw_cfg(int nx, double W, double H) {
    Config c = cfg(1, W, H, nx, 1, {0.0, 0.5 * H, H});
    LawOptions o; o.kind = 1; o.c_cps = o.c_cpl = 500; o.c_ks = o.c_kl = 1e-12; o.c_Tm = 1000; o.c_dTm = 1; o.T_min = 200; o.T_max = 5000;
    o.c_rho = 7000; o.c_mu = 1e-9; o.vapor = false;
    c.law = o; c.rho_ref = 7000; c.evap_mode = 0; c.flow = true; c.flow_model = 2; c.T0 = 1200; c.height = 2.5e-4; c.slices[0].y0 = 0.0; c.max_substeps = 100000000;
    return c;
}
static void dam_parity() {
    Config c = sw_cfg(200, 0.2, 0.002);
    Pair p(c, gopt);
    auto disp = [&](int, int f, double x) { return f == 1 ? (x < 0 ? 0.4e-3 : -0.4e-3) : 0.0; };
    p.cpu.initialize_displacement(disp); p.gpu.initialize_displacement(disp);
    std::vector<double> JB(6 * p.cpu.dims().n_em, 0.0);
    for (int e = 0; e < p.cpu.dims().n_em; ++e) { JB[6 * e + 2] = 1e6; JB[6 * e + 3] = -7.0; }
    long retries = 0;
    for (int k = 0; k < 20; ++k) { StepInput in; in.dt = 1e-3; in.em_power_W.assign(p.cpu.dims().n_em, 0.0); in.em_JB = JB; p.step(in, true); retries += p.gpu.ledger_closes() ? 0 : 1000; }
    chk::check(p.mismatches == 0 && retries == 0, "GPU = CPU bit for bit: dam break with CFL halvings");
}
static void eject_parity() {
    const int nx = 20; const double W = 2e-3, H = 1e-4, delta = 2e-5;
    std::vector<double> rows; for (int j = 0; j <= 4; ++j) rows.push_back((H - delta) * j / 4); for (int j = 1; j <= 4; ++j) rows.push_back(H - delta + delta * j / 4); rows.back() = H;
    for (int open : {1, 0}) {
        Config c = cfg(1, W, H, nx, 1, rows);
        LawOptions o; o.kind = 1; o.c_cps = o.c_cpl = 500; o.c_ks = o.c_kl = 1e-12; o.c_Tm = 1000; o.c_dTm = 1; o.T_min = 200; o.T_max = 5000; o.c_rho = 7000; o.c_mu = 6e-3; o.c_sigma = 1.8; o.vapor = false;
        c.law = o; c.rho_ref = 7000; c.evap_mode = 0; c.flow = true; c.flow_model = 2; c.eject = true; c.open_x = open; c.eject_fraction = 0.02; c.T0 = 800; c.slices[0].y0 = 0.0;
        Pair p(c, gopt);
        p.init([&](int, double, double y) { return y > H - delta ? 1200.0 : 800.0; });
        std::vector<double> JB(6 * nx, 0.0); for (int e = 0; e < nx; ++e) { JB[6 * e + 2] = 1e7; JB[6 * e + 4] = -1.0; }
        for (int k = 0; k < 60; ++k) { StepInput in; in.dt = 2e-5; in.em_power_W.assign(nx, 0.0); in.em_JB = JB; p.step(in, k % 10 == 0); }
        chk::check(p.mismatches == 0 && !px_is_zero(p.gpu.ledger().eject_E), open ? "GPU = CPU bit for bit: edge outflow into the ejected-melt reservoir" : "GPU = CPU bit for bit: departure ejection");
    }
}
// device accumulation (field data resident on the caller's session) = host accumulation, bit for bit
static void accdev() {
    hall_gpu::SessionOptions so; so.require_uuid = gopt.require_uuid; so.nvrtc_library = gopt.nvrtc_library;
    hall_gpu::Session session(so);
    GpuOptions b = gopt; b.session = &session;
    const double W = 0.005, H = 0.0002;
    Config c = cfg(4, W, H, 5, 2, graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
    Config::enable_v1_2(c);
    bool same = true;
    {
        Owner g(c, std::make_unique<GpuEngine>(b)), h(c);
        auto prof = [&](int s, double x, double) { return 2600.0 + 2.0e5 * (x + 0.5 * W) + 30.0 * s; };
        g.initialize_temperature(prof); h.initialize_temperature(prof);
        const int nem = g.dims().n_em;
        std::vector<double> JB(6 * nem); for (int e = 0; e < nem; ++e) { JB[6 * e + 2] = 5e8; JB[6 * e + 3] = 20.0; JB[6 * e + 4] = -0.05; }
        auto dP = session.alloc(nem * sizeof(double), "test power"), dJB = session.alloc(6 * nem * sizeof(double), "test JB");
        auto stage = session.alloc_pinned(7 * nem * sizeof(double), "test stage");
        for (int k = 0; k < 40; ++k) {
            for (int q = 0; q < 8; ++q) {
                std::vector<double> P(nem); for (int e = 0; e < nem; ++e) P[e] = 3e6 * (1 + 0.2 * std::sin(0.3 * (8 * k + q) + e));
                session.synchronize("test stage");
                std::memcpy(stage.pointer, P.data(), nem * sizeof(double)); std::memcpy(static_cast<char*>(stage.pointer) + nem * sizeof(double), JB.data(), 6 * nem * sizeof(double));
                session.upload(dP, 0, stage, 0, nem * sizeof(double)); session.upload(dJB, 0, stage, nem * sizeof(double), 6 * nem * sizeof(double));
                g.accumulate_device(reinterpret_cast<const double*>(dP.address), reinterpret_cast<const double*>(dJB.address), 1.25e-12);
                h.accumulate(P, JB, 1.25e-12);
            }
            StepInput in; in.use_accumulated = true;
            g.advance(in); h.advance(in);
            same = same && g.state_sha256() == h.state_sha256();
        }
        session.synchronize("test end");
        session.release(dP); session.release(dJB); session.release(stage);
    }
    chk::check(same, "device accumulation on the caller's session = host accumulation on the CPU engine, bit for bit (40 windows of 8 field steps)");
    session.close();
}
// Hall foil (125 k x 2 bodies, 25 columns, 262,500 cells), every v1.2 option: cost per field step with sub-cycling
// (device accumulation from resident power and J, B; owner step every K field steps), device time per kernel group.
static void timing() {
    hall_gpu::SessionOptions so; so.require_uuid = gopt.require_uuid; so.nvrtc_library = gopt.nvrtc_library;
    hall_gpu::Session session(so);
    GpuOptions b = gopt; b.session = &session;
    Config c;
    for (int k = 0; k < 125; ++k) for (int body = 0; body < 2; ++body) c.slices.push_back({body, k, -0.0025, body ? 0.0025 : -0.0027});
    c.width = 0.005; c.thickness = 0.0002; c.height = 0.0002; c.nx_em = 25; c.ny_em = 1; c.rx = 1;
    c.y_edges = graded_rows(0.0002, 0.25e-6, 1.3, 10e-6, 1);
    Config::enable_v1_2(c);
    {
        GpuOptions bp = b; bp.profile = true;
        auto eng = std::make_unique<GpuEngine>(bp); GpuEngine* ge = eng.get();
        Owner g(c, std::move(eng));
        g.initialize_temperature([](int, double, double) { return 2500.0; });   // molten surface: the layer pass is active
        const int nem = g.dims().n_em;
        std::vector<double> JB(6 * nem); for (int e = 0; e < nem; ++e) { JB[6 * e + 2] = 1e8; JB[6 * e + 3] = 5.0; JB[6 * e + 4] = -0.01; }
        std::vector<double> P(nem, 2e5);
        auto dP = session.alloc(nem * sizeof(double), "power"), dJB = session.alloc(6 * nem * sizeof(double), "JB");
        auto stage = session.alloc_pinned(7 * nem * sizeof(double), "stage");
        std::memcpy(stage.pointer, P.data(), nem * sizeof(double)); std::memcpy(static_cast<char*>(stage.pointer) + nem * sizeof(double), JB.data(), 6 * nem * sizeof(double));
        session.upload(dP, 0, stage, 0, nem * sizeof(double)); session.upload(dJB, 0, stage, nem * sizeof(double), 6 * nem * sizeof(double));
        session.synchronize("resident field data");
        const double dtf = 1e-14;
        const auto* pP = reinterpret_cast<const double*>(dP.address); const auto* pJ = reinterpret_cast<const double*>(dJB.address);
        // warm-up (compile caches, first CFL estimate)
        for (int k = 0; k < 3; ++k) { g.accumulate_device(pP, pJ, dtf); StepInput in; in.use_accumulated = true; g.advance(in); }
        chk::put("timing_cells", double(g.dims().n_cells)); chk::put("timing_stability_tau_s", g.stability_tau()); chk::put("timing_field_dt_s", dtf);
        for (int K : {1, 4, 8, 16, 32}) {
            const int windows = K >= 16 ? 12 : 24, n = K * windows;
            session.synchronize("timing start");
            double owner_ms = 0, engine_ms = 0, acc_ms = 0;
            std::int64_t nsub = 0;
            const auto t0 = Clock::now();
            for (int f = 0; f < n; ++f) {
                const auto ta = Clock::now();
                g.accumulate_device(pP, pJ, dtf);
                acc_ms += 1e3 * std::chrono::duration<double>(Clock::now() - ta).count();
                if ((f + 1) % K == 0) {
                    const auto tb = Clock::now();
                    StepInput in; in.use_accumulated = true;
                    nsub += g.advance(in).n_sub;
                    owner_ms += 1e3 * std::chrono::duration<double>(Clock::now() - tb).count();
                    engine_ms += ge->last_step_ms();
                }
            }
            session.synchronize("timing end");
            const double total = 1e3 * std::chrono::duration<double>(Clock::now() - t0).count();
            const std::string T = "timing_v12_K" + std::to_string(K);
            chk::put(T + "_ms_per_field_step", total / n);
            chk::put(T + "_owner_ms_per_window", owner_ms / windows); chk::put(T + "_engine_ms_per_window", engine_ms / windows);
            chk::put(T + "_accumulate_host_ms_per_field_step", acc_ms / n);
            chk::put(T + "_substeps_per_window", double(nsub) / windows);
        }
        chk::put_raw("timing_v12_profile", ge->profile_json());
        chk::check(g.ledger_closes(), "exact ledgers after the timing run (full Hall foil, v1.2)");
        session.release(dP); session.release(dJB); session.release(stage);
    }
    // v1.1 options and v1.2 options without sub-cycling (one owner step per field step), for comparison
    for (int v : {11, 12}) {
        Config cc = c; if (v == 11) { cc.flow_model = 1; cc.eject = false; cc.diagnostics = false; }
        Owner g(cc, std::make_unique<GpuEngine>(b));
        g.initialize_temperature([](int, double, double) { return 2500.0; });
        const int nem = g.dims().n_em;
        std::vector<double> JB(6 * nem); for (int e = 0; e < nem; ++e) { JB[6 * e + 2] = 1e8; JB[6 * e + 3] = 5.0; JB[6 * e + 4] = -0.01; }
        StepInput in; in.dt = 1e-14; in.em_power_W.assign(nem, 2e5); in.em_JB = JB;
        for (int k = 0; k < 3; ++k) g.advance(in);
        const auto t0 = Clock::now();
        for (int k = 0; k < 50; ++k) g.advance(in);
        chk::put("timing_v" + std::to_string(v) + "_single_step_ms", 1e3 * std::chrono::duration<double>(Clock::now() - t0).count() / 50);
    }
    session.synchronize("timing");
    const auto rep = session.close();
    chk::check(rep.released && !rep.quarantined, "timing session closes cleanly");
}

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: test_gpu OUT NVRTC UUID [suite...]\n"); return 2; }
    const std::string out = argv[1];
    gopt.nvrtc_library = argv[2]; gopt.require_uuid = argv[3];
    std::set<std::string> want; for (int i = 4; i < argc; ++i) want.insert(argv[i]);
    { GpuEngine probe(gopt); chk::put_raw("gpu_record", probe.record_json()); }
    auto timed = [&](const char* s, const std::function<void()>& f) {
        if (!want.empty() && !want.count(s)) return;
        const auto t = Clock::now(); f();
        chk::put(std::string("seconds_") + s, std::chrono::duration<double>(Clock::now() - t).count());
        std::printf("suite %s done\n", s); std::fflush(stdout);
    };
    timed("v1", v1_evaporation);
    timed("v11", v11_ledger);
    timed("quantized", v11_quantized);
    timed("leveling", v11_leveling);
    timed("merkle", merkle);
    timed("borrowed", borrowed);
    timed("hall", hall_timing);
    timed("v12", v12_parity);
    timed("dam", dam_parity);
    timed("eject", eject_parity);
    timed("accdev", accdev);
    timed("timing", timing);
    std::string name = "test_gpu"; for (const auto& w : want) name += "-" + w;
    return chk::finish(name, out);
}
