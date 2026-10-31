// test_gpu.cpp -- CPU reference versus NVRTC GPU engine: bit identity of every accepted state, step record,
// patch report, quote and sigma vector; timing on the full Hall foil geometry.
//   test_gpu OUT_DIR NVRTC_LIBRARY GPU_UUID [suite ...]
#include "check.hpp"
#include "phase_owner.hpp"
#include "gpu_engine.hpp"
#include <chrono>
#include <cstring>
#include <functional>
#include <random>
#include <set>

using namespace ss304_phase_v1;
using Clock = std::chrono::steady_clock;
static GpuOptions gopt;

static bool same_doubles(const std::vector<double>& a, const std::vector<double>& b) { return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * 8) == 0; }
static bool same_reports(const std::vector<PatchReport>& a, const std::vector<PatchReport>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double x[] = {a[i].area_m2, a[i].mean_T, a[i].max_T, a[i].liquid_fraction, a[i].melt_depth_m, a[i].surface_density, a[i].flux_atoms_m2_s[0], a[i].flux_atoms_m2_s[1], a[i].flux_atoms_m2_s[2], a[i].expected_atoms[0], a[i].expected_atoms[1], a[i].expected_atoms[2], a[i].cumulative_evaporated_atoms[0], a[i].recession_m};
        const double y[] = {b[i].area_m2, b[i].mean_T, b[i].max_T, b[i].liquid_fraction, b[i].melt_depth_m, b[i].surface_density, b[i].flux_atoms_m2_s[0], b[i].flux_atoms_m2_s[1], b[i].flux_atoms_m2_s[2], b[i].expected_atoms[0], b[i].expected_atoms[1], b[i].expected_atoms[2], b[i].cumulative_evaporated_atoms[0], b[i].recession_m};
        if (std::memcmp(x, y, sizeof x)) return false;
    }
    return true;
}
static bool same_record(const StepRecord& a, const StepRecord& b) {
    auto eq = [](pt_i128 x, pt_i128 y) { return pt_eq(x, y); };
    bool ok = a.n_sub == b.n_sub && a.merges == b.merges && std::memcmp(&a.max_T, &b.max_T, 8) == 0 && a.inputs_sha256 == b.inputs_sha256;
    for (auto [x, y] : {std::pair{a.dE, b.dE}, {a.joule, b.joule}, {a.loads, b.loads}, {a.deposit_E, b.deposit_E}, {a.boundary, b.boundary}, {a.evap_ref, b.evap_ref}, {a.evap_kin, b.evap_kin}, {a.release_E, b.release_E}}) ok = ok && eq(x, y);
    for (int el = 0; el < 3; ++el) ok = ok && eq(a.dN[el], b.dN[el]) && eq(a.evap_N[el], b.evap_N[el]) && eq(a.rel_N[el], b.rel_N[el]) && eq(a.dep_N[el], b.dep_N[el]);
    return ok;
}

// Runs the same input sequence through a CPU and a GPU owner and compares after every step.
struct Pair {
    Owner cpu, gpu;
    long steps = 0, mismatches = 0;
    explicit Pair(const Config& c) : cpu(c), gpu(c, std::make_unique<GpuEngine>(gopt)) {}
    template <class F> void init(F f) { cpu.initialize_temperature(f); gpu.initialize_temperature(f); }
    bool compare(bool full) {
        bool ok = cpu.state_sha256() == gpu.state_sha256();
        if (full) ok = ok && same_reports(cpu.patch_reports(), gpu.patch_reports()) && same_doubles(cpu.sigma_em(), gpu.sigma_em());
        return ok;
    }
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

static void parity_evaporation() {
    const double W = 0.005, H = 0.0002;
    Config c = cfg(2, W, H, 5, 2, graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
    Pair p(c);
    p.init([&](int s, double x, double y) { (void)y; return 3500.0 + 2.8e5 * (x + 0.5 * W) + 50.0 * s; });
    std::uint64_t seq = 0;
    for (int k = 0; k < 300; ++k) {
        StepInput in; in.dt = 1e-11; in.em_power_W.resize(p.cpu.dims().n_em);
        for (int e = 0; e < p.cpu.dims().n_em; ++e) in.em_power_W[e] = 6e6 * (1 + 0.1 * (e % 5));
        for (int s = 0; s < 2; ++s) {
            ++seq;
            const double y0 = c.slices[s].y0;
            in.loads.push_back({s, -0.5 * W + 0.0011, y0 + H, 2e-9, "t/wall", "wall", seq});
            in.loads.push_back({s, -0.5 * W + 0.0031, y0, -1e-9, "t/emission", "emission", seq});
            in.deposits.push_back({s, -0.5 * W + 0.0021, y0 + H, 1, 1e7, 1e-12, "t/cr-return", seq});
        }
        try { p.step(in, k % 25 == 0); } catch (const PhaseRefusal&) { break; }
    }
    chk::check(p.mismatches == 0, "GPU equals CPU bit for bit: evaporation, loads, deposits");
    chk::check(p.cpu.ledger_closes() && p.gpu.ledger_closes(), "exact ledgers on both engines");
    chk::put("parity_evaporation_steps", double(p.steps)); chk::put("parity_evaporation_mismatches", double(p.mismatches));
    chk::put("parity_evaporation_max_T", p.gpu.max_temperature_K());
}

static void parity_quantized() {
    const double W = 0.002, H = 0.0002;
    Config c = cfg(1, W, H, 2, 2, graded_rows(H, 0.5e-6, 1.4, 20e-6, 1));
    c.evap_mode = 2;
    c.patches = {{0, 1, 0.0, 0.5}, {0, 1, 0.5, 1.0}, {0, 0, 0.0, 0.3}, {0, 0, 0.3, 1.0}};
    Pair p(c);
    p.init([](int, double, double) { return 3900.0; });
    std::mt19937_64 rng(9); std::uniform_real_distribution<double> U(0, 1);
    for (int k = 0; k < 150; ++k) {
        const Quote qc = p.cpu.quote(1e-12), qg = p.gpu.quote(1e-12);
        if (qc.sha256 != qg.sha256) ++p.mismatches;
        StepInput in; in.dt = 1e-12; in.em_power_W.assign(p.cpu.dims().n_em, 4e5); in.quote_sha256 = qc.sha256;
        for (int pp = 0; pp < p.cpu.n_patches(); ++pp) for (int el = 0; el < 3; ++el) {
            const double f = qc.expected_atoms[pp][el] / 1e6, fl = std::floor(f);
            const std::uint64_t n = std::uint64_t(fl) + (U(rng) < f - fl ? 1 : 0);
            if (n) in.releases.push_back({pp, el, n, 1e6, double(n) * 1e6 * 2 * boltzmann * qc.mean_T[pp]});
        }
        p.step(in, k % 10 == 0);
    }
    chk::check(p.mismatches == 0, "GPU equals CPU bit for bit: quantized releases and quotes");
    chk::put("parity_quantized_steps", double(p.steps));
}

static void parity_stefan() {
    const double H = 6e-4;
    Config c = cfg(1, 1e-4, H, 1, 1, [] { std::vector<double> e; for (int j = 0; j <= 150; ++j) e.push_back(6e-4 * j / 150); e.back() = 6e-4; return e; }());
    LawOptions o; o.kind = 1; o.c_T0 = 300; o.c_cps = o.c_cpl = 500; o.c_ks = o.c_kl = 20; o.c_Tm = 1000; o.c_dTm = 0.05; o.c_L = 2.5e5; o.T_min = 200; o.T_max = 2e5;
    c.law = o; c.evap_mode = 0; c.bc_kind = {1, 0, 0, 0}; c.bc_T = {1200, 0, 0, 0}; c.T0 = 900;
    Pair p(c);
    for (int k = 0; k < 60; ++k) { StepInput in; in.dt = 1e-5; in.em_power_W = {0.0}; p.step(in, k % 10 == 0); }
    chk::check(p.mismatches == 0, "GPU equals CPU bit for bit: Stefan problem with prescribed temperature and sub-steps");
    chk::put("parity_stefan_substeps_per_step", double(p.cpu.substeps_for(1e-5)));
}

static void parity_merge() {
    const double H = 0.0002;
    Config c = cfg(1, 1e-4, H, 1, 1, {0.0, 20e-9, 40e-9, 80e-9, 160e-9, 320e-9, 640e-9, 1.28e-6, 2.56e-6, 5.12e-6, 10.24e-6, 20.48e-6, 40.96e-6, 81.92e-6, H});
    c.evap_face = {1, 0, 0, 0};
    Pair p(c);
    p.init([](int, double, double) { return 4990.0; });
    int merges = 0;
    for (int k = 0; k < 150000 && merges < 1; ++k) {
        StepInput in; in.dt = 1e-12; in.em_power_W = {0.0};
        p.step(in, false);
        merges = int(std::lround(double(p.gpu.state().flags[0] & 1)));
        if (k % 1000 == 0 && p.mismatches) break;
    }
    chk::check(merges == 1 && p.mismatches == 0, "GPU equals CPU bit for bit through a surface-cell merge");
    chk::put("parity_merge_steps", double(p.steps));
}

static void hall_timing() {
    Config c;
    for (int k = 0; k < 125; ++k) for (int body = 0; body < 2; ++body) c.slices.push_back({body, k, -0.0025, body ? 0.0025 : -0.0027});
    c.width = 0.005; c.thickness = 0.0002; c.height = 0.0002; c.nx_em = 25; c.ny_em = 1; c.rx = 1;
    c.y_edges = graded_rows(0.0002, 0.25e-6, 1.3, 10e-6, 1);
    Owner g(c, std::make_unique<GpuEngine>(gopt)), h(c);
    std::vector<double> P(g.dims().n_em, 2e5);
    StepInput in; in.dt = 1e-14; in.em_power_W = P;
    for (int k = 0; k < 3; ++k) { g.advance(in); h.advance(in); }
    chk::check(g.state_sha256() == h.state_sha256(), "full Hall geometry: GPU equals CPU after 3 steps");
    auto t0 = Clock::now();
    for (int k = 0; k < 3; ++k) h.advance(in);
    const double cpu_ms = 1e3 * std::chrono::duration<double>(Clock::now() - t0).count() / 3;
    t0 = Clock::now();
    const int n = 300;
    for (int k = 0; k < n; ++k) g.advance(in);
    const double gpu_ms = 1e3 * std::chrono::duration<double>(Clock::now() - t0).count() / n;
    chk::put("hall_full_cells", double(g.dims().n_cells)); chk::put("hall_full_slices", double(g.dims().S));
    chk::put("hall_full_cpu_ms_per_step", cpu_ms); chk::put("hall_full_gpu_ms_per_step_incl_owner", gpu_ms);
}

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: test_gpu OUT NVRTC UUID [suite...]\n"); return 2; }
    const std::string out = argv[1];
    gopt.nvrtc_library = argv[2]; gopt.require_uuid = argv[3];
    std::set<std::string> want; for (int i = 4; i < argc; ++i) want.insert(argv[i]);
    auto on = [&](const char* s) { return want.empty() || want.count(s); };
    {
        GpuEngine probe(gopt);
        chk::put_raw("gpu_record", probe.record_json());
    }
    auto timed = [&](const char* s, const std::function<void()>& f) {
        if (!on(s)) return; const auto t = Clock::now(); f();
        chk::put(std::string("seconds_") + s, std::chrono::duration<double>(Clock::now() - t).count());
        std::printf("suite %s done\n", s); std::fflush(stdout);
    };
    timed("evaporation", parity_evaporation);
    timed("quantized", parity_quantized);
    timed("stefan", parity_stefan);
    timed("merge", parity_merge);
    timed("hall", hall_timing);
    std::string name = "test_gpu"; for (const auto& w : want) name += "-" + w;
    return chk::finish(name, out);
}
