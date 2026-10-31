// test_owner.cpp -- benchmarks and ledger tests of the SS304 phase thermal owner v1 (CPU reference).
//   test_owner OUT_DIR [suite ...]      suites: uniform conduction mode stefan ledger quantized merge bennett
//                                                checkpoint refusals hall
#include "check.hpp"
#include "phase_owner.hpp"
#include <chrono>
#include <cstring>
#include <functional>
#include <random>
#include <set>

using namespace ss304_phase_v1;
using Clock = std::chrono::steady_clock;

static std::vector<double> uniform_rows(double H, int n) { std::vector<double> e(n + 1); for (int j = 0; j <= n; ++j) e[j] = H * j / n; e.back() = H; return e; }
static Config make_cfg(int S, double W, double H, double dz, int nx_em, int ny_em, int rx, std::vector<double> rows) {
    Config c;
    for (int s = 0; s < S; ++s) c.slices.push_back({0, s, -0.5 * W, 0.0});
    c.width = W; c.thickness = H; c.height = dz; c.nx_em = nx_em; c.ny_em = ny_em; c.rx = rx; c.y_edges = std::move(rows);
    return c;
}
static LawOptions constant_law(double cp, double k, double Tm, double dTm, double L, double Tmin, double Tmax) {
    LawOptions o; o.kind = 1; o.c_T0 = 300; o.c_cps = o.c_cpl = cp; o.c_ks = o.c_kl = k; o.c_Tm = Tm; o.c_dTm = dTm; o.c_L = L; o.T_min = Tmin; o.T_max = Tmax; o.c_rhoe = 1e-6;
    return o;
}
static double ierfc(double x) { return std::exp(-x * x) / std::sqrt(M_PI) - x * std::erfc(x); }
static std::string jarr(const std::vector<double>& v) { std::string s = "["; for (std::size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + chk::num(v[i]); return s + "]"; }
static double secs(Clock::time_point a) { return std::chrono::duration<double>(Clock::now() - a).count(); }

// ------------------------------------------------------------------------------------------------ uniform
// Adiabatic slice with uniform Joule power density: exact solution h(t) = h0 + p t / rho_ref through the
// solid, the mushy zone, the liquid and up to the domain limit; refusal just above 5000 K.
static void uniform_heating() {
    Config c = make_cfg(1, 0.004, 0.0002, 0.0002, 2, 1, 2, uniform_rows(0.0002, 4));
    c.evap_mode = 0;
    Owner o(c);
    const double p = 5e17;   // W/m3
    std::vector<double> P(o.dims().n_em);
    for (int e = 0; e < o.dims().n_em; ++e) P[e] = p * (0.002 * 0.0002 * 0.0002);
    const double dt = 2e-11, rho = c.rho_ref;
    double max_rel = 0, t = 0; bool passed1500 = false, melted = false, closes = true; long steps = 0;
    double h0 = pt_ss_h(300.0);
    for (;;) {
        StepInput in; in.dt = dt; in.em_power_W = P;
        const double h_next = h0 + p * (t + dt) / rho;
        if (h_next > pt_ss_h(5000.0)) {
            const std::string before = o.state_sha256();
            bool refused = false; std::string kind;
            try { o.prepare(in); } catch (const PhaseRefusal& r) { refused = true; kind = r.kind; }
            chk::check(refused && kind == "temperature-above-domain", "typed refusal just above 5000 K");
            chk::check(o.state_sha256() == before, "refused step leaves the accepted state unchanged");
            break;
        }
        o.advance(in); t += dt; ++steps;
        int st = 0; const double T_exact = pt_T_of_h(&o.law(), h_next, &st);
        const auto T = o.temperatures();
        for (double x : T) max_rel = std::max(max_rel, std::abs(x - T_exact) / T_exact);
        if (T[0] > 1500) passed1500 = true;
        const auto rep = o.patch_reports();
        if (rep[0].liquid_fraction > 0 && rep[0].liquid_fraction < 1) melted = true;
        if (steps % 50 == 0) closes = closes && o.ledger_closes();
    }
    chk::check(passed1500 && melted, "no refusal at 1500 K or in the mushy zone");
    chk::check(max_rel < 1e-12, "uniform heating follows the exact enthalpy solution");
    chk::check(closes && o.ledger_closes(), "exact ledger closure (uniform heating)");
    chk::put("uniform_steps", double(steps)); chk::put("uniform_max_rel_T_error", max_rel); chk::put("uniform_final_T", o.max_temperature_K());
}

// ------------------------------------------------------------------------------------------------ conduction
// 1-D semi-infinite solid, constant surface flux q: T(y,t) = T0 + (2q/k) sqrt(alpha t) ierfc(y / (2 sqrt(alpha t))).
static double conduction_case(int rows, double* surf_err) {
    const double H = 0.001, k = 20, cp = 500, rho = 8000, q = 1e7, tend = 1e-3, alpha = k / (rho * cp);
    Config c = make_cfg(1, 0.001, H, 0.001, 1, 1, 1, uniform_rows(H, rows));
    c.law = constant_law(cp, k, 1e5, 1, 0, 200, 2e5); c.evap_mode = 0;
    Owner o(c);
    const int nsteps = 400; const double dt = tend / nsteps;
    for (int s = 0; s < nsteps; ++s) {
        StepInput in; in.dt = dt; in.em_power_W = {0.0};
        in.loads.push_back({0, 0.0, 0.0, q * 0.001 * 0.001 * dt, "bench-flux", "constant-flux", std::uint64_t(s + 1)});
        o.advance(in);
    }
    const auto T = o.temperatures();
    double err = 0, scale = 0;
    for (int j = 0; j < rows; ++j) {
        const double y = 0.5 * (c.y_edges[j] + c.y_edges[j + 1]);
        const double Te = 300 + 2 * q / k * std::sqrt(alpha * tend) * ierfc(y / (2 * std::sqrt(alpha * tend)));
        err = std::max(err, std::abs(T[j] - Te)); scale = std::max(scale, Te - 300);
    }
    const double Ts = 300 + 2 * q * std::sqrt(tend / M_PI) / std::sqrt(k * rho * cp);
    *surf_err = std::abs(T[0] - Ts) / (Ts - 300);
    chk::check(o.ledger_closes(), "exact ledger closure (conduction)");
    return err / scale;
}
static void conduction() {
    double s1, s2, s3;
    const double e1 = conduction_case(250, &s1), e2 = conduction_case(500, &s2), e3 = conduction_case(1000, &s3);
    chk::put("flux_benchmark_profile_err_250_500_1000", 0); chk::put_raw("flux_benchmark_profile_rel_err", jarr({e1, e2, e3}));
    chk::put_raw("flux_benchmark_surface_rel_err", jarr({s1, s2, s3}));
    chk::check(e3 < e2 && e2 < e1 && e3 < 0.01, "constant-flux conduction converges to the erfc solution");
    chk::put("flux_benchmark_order_last", std::log2(e2 / e3));
}

// ------------------------------------------------------------------------------------------------ mode decay
// 2-D adiabatic slice, T = T0 + A cos(pi x/W) cos(pi y/H): explicit FV amplitude factor (1 - alpha tau lambda_h)
// per sub-step (exact discrete solution) and exp(-alpha pi^2 (1/W^2 + 1/H^2) t) (continuum).
static void mode_decay() {
    const double W = 4e-4, H = 2e-4, k = 20, cp = 500, rho = 8000, alpha = k / (rho * cp);
    const int nx = 40, ny = 20;
    Config c = make_cfg(1, W, H, 1e-4, 10, 1, 4, uniform_rows(H, ny));
    c.law = constant_law(cp, k, 1e5, 1, 0, 200, 2e5); c.evap_mode = 0;
    Owner o(c);
    const double A = 50, x0 = -0.5 * W;
    o.initialize_temperature([&](int, double x, double y) { return 1000 + A * std::cos(M_PI * (x - x0) / W) * std::cos(M_PI * y / H); });
    const double dx = W / nx, dy = H / ny;
    const double lam = (2 - 2 * std::cos(M_PI / nx)) / (dx * dx) + (2 - 2 * std::cos(M_PI / ny)) / (dy * dy);
    const double tend = 0.05 * H * H / alpha; const int steps = 50; const double dt = tend / steps;
    std::int64_t nsub_total = 0;
    for (int s = 0; s < steps; ++s) { StepInput in; in.dt = dt; in.em_power_W.assign(o.dims().n_em, 0.0); nsub_total += o.advance(in).n_sub; }
    const double tau = dt / double(o.substeps_for(dt));
    const double disc = A * std::pow(1 - alpha * tau * lam, double(nsub_total)), cont = A * std::exp(-alpha * M_PI * M_PI * (1 / (W * W) + 1 / (H * H)) * tend);
    const auto T = o.temperatures();
    double ed = 0, ec = 0;
    for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const double sh = std::cos(M_PI * (i + 0.5) / nx) * std::cos(M_PI * (j + 0.5) / ny);
        ed = std::max(ed, std::abs(T[j * nx + i] - 1000 - disc * sh)); ec = std::max(ec, std::abs(T[j * nx + i] - 1000 - cont * sh));
    }
    chk::put("mode_amplitude_discrete", disc); chk::put("mode_amplitude_continuum", cont);
    chk::put("mode_max_err_vs_discrete_K", ed); chk::put("mode_max_err_vs_continuum_K", ec);
    chk::check(ed < 1e-9 * A, "2-D mode decays with the exact discrete eigenvalue");
    chk::check(ec < 0.01 * A, "2-D mode within 1% of the continuum decay");
    chk::check(o.ledger_closes(), "exact ledger closure (mode decay)");
}

// ------------------------------------------------------------------------------------------------ Stefan
// Two-phase Neumann problem: solid at Ti < Tm, face held at Ts > Tm from t = 0 (equal phase properties).
static double neumann_lambda(double St_l, double ratio) {   // ratio = (Tm-Ti)/(Ts-Tm)
    auto f = [&](double l) { return std::exp(-l * l) / std::erf(l) - ratio * std::exp(-l * l) / std::erfc(l) - l * std::sqrt(M_PI) / St_l; };
    double a = 1e-6, b = 3;
    for (int i = 0; i < 200; ++i) { const double m = 0.5 * (a + b); (f(a) * f(m) <= 0 ? b : a) = m; }
    return 0.5 * (a + b);
}
static void stefan() {
    const double k = 20, cp = 500, rho = 8000, L = 2.5e5, Tm = 1000, Ti = 900, Ts = 1200, alpha = k / (rho * cp), H = 6e-4, tend = 3e-3;
    const double lam = neumann_lambda(cp * (Ts - Tm) / L, (Tm - Ti) / (Ts - Tm));
    const double s_exact = 2 * lam * std::sqrt(alpha * tend);
    std::vector<double> errs, terrs;
    for (int rows : {150, 300, 600}) {
        Config c = make_cfg(1, 1e-4, H, 1e-4, 1, 1, 1, uniform_rows(H, rows));
        c.law = constant_law(cp, k, Tm, 0.05, L, 200, 2e5); c.evap_mode = 0;
        c.bc_kind = {1, 0, 0, 0}; c.bc_T = {Ts, 0, 0, 0};
        c.T0 = Ti;
        Owner o(c);
        const int steps = 300; const double dt = tend / steps;
        for (int s = 0; s < steps; ++s) { StepInput in; in.dt = dt; in.em_power_W = {0.0}; o.advance(in); }
        const auto rep = o.patch_reports();
        const double s_num = rep[0].melt_depth_m;   // LowY patch: contiguous melt layer
        errs.push_back(std::abs(s_num - s_exact) / s_exact);
        const auto T = o.temperatures();
        double te = 0;
        for (int j = 0; j < rows; ++j) {
            const double y = (j + 0.5) * H / rows, eta = y / (2 * std::sqrt(alpha * tend));
            const double Te = y < s_exact ? Ts - (Ts - Tm) * std::erf(eta) / std::erf(lam) : Ti + (Tm - Ti) * std::erfc(eta) / std::erfc(lam);
            if (std::abs(y - s_exact) > 3 * H / rows) te = std::max(te, std::abs(T[j] - Te));
        }
        terrs.push_back(te);
        chk::check(o.ledger_closes(), "exact ledger closure (Stefan, boundary energy included)");
    }
    chk::put("stefan_lambda", lam); chk::put("stefan_front_exact_m", s_exact);
    chk::put_raw("stefan_front_rel_err_150_300_600", jarr(errs)); chk::put_raw("stefan_T_max_err_K_150_300_600", jarr(terrs));
    chk::check(errs[2] < 0.01 && errs[2] < errs[0], "Stefan front within 1% and converging");
    chk::check(terrs[2] < 1.0, "Stefan temperature profile within 1 K away from the front");
}

// ------------------------------------------------------------------------------------------------ ledger
// SS304 law through melting into strong vaporization: nonuniform Joule, wall heat, emission cooling, Fe
// deposits, expected HKL evaporation on all faces. Exact closure every step; patch cumulative atoms equal the
// ledger.
static void ledger_run() {
    const double W = 0.005, H = 0.0002, dz = 0.0002;
    Config c = make_cfg(2, W, H, dz, 5, 1, 2, graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
    c.evap_mode = 1;
    Owner o(c);
    const int nem = o.dims().n_em;
    std::vector<double> base(nem);
    for (int e = 0; e < nem; ++e) base[e] = 4.5e17 * (1.0 + 0.3 * ((e % 5) - 2) / 2.0) * (W / 5 * H * dz);
    const double dt = 1e-11;
    long steps = 0; double maxT = 0; bool exact = true; std::uint64_t seq = 0;
    const auto t0 = Clock::now();
    double evap_kg = 0;
    while (steps < 20000) {
        StepInput in; in.dt = dt; in.em_power_W = base;
        const double tt = o.time_s();
        for (int s = 0; s < 2; ++s) {
            ++seq;
            in.loads.push_back({s, -0.5 * W + 0.0011, H, 2e-9, "test/wall-heat", "collected-wall-heat", seq});
            in.loads.push_back({s, -0.5 * W + 0.0031, 0.0, -1e-9, "test/emission-cooling", "emission-cooling", seq});
            in.deposits.push_back({s, -0.5 * W + 0.0021, H, 0, 1e7, 1e7 * 1e-19, "test/fe-return", seq});
        }
        (void)tt;
        StepRecord r;
        try { r = o.advance(in); } catch (const PhaseRefusal& e) { chk::put_text("ledger_run_stop", e.what()); break; }
        ++steps; maxT = std::max(maxT, r.max_T);
        if (steps % 200 == 0) exact = exact && o.ledger_closes();
        if (maxT > 4500) break;
    }
    const double elapsed = secs(t0);
    exact = exact && o.ledger_closes();
    const auto rep = o.patch_reports();
    for (std::size_t p = 0; p < rep.size(); ++p) evap_kg += rep[p].cumulative_evaporated_kg;
    // patch cumulative evaporated atoms (reported as doubles) versus the exact ledger words
    double patch_atoms = 0, ledger_atoms = 0;
    for (const auto& p : rep) for (int el = 0; el < 3; ++el) patch_atoms += p.cumulative_evaporated_atoms[el];
    for (int el = 0; el < 3; ++el) ledger_atoms += pt_to_double(o.ledger().evap_N[el], PT_N_INV);
    chk::check(exact, "exact energy and atom ledgers through melting and vaporization (every 200 steps and final)");
    chk::check(maxT > 4000, "run reached strong vaporization (> 4000 K)");
    chk::check(std::abs(patch_atoms - ledger_atoms) <= 1e-12 * ledger_atoms, "patch evaporation totals equal the ledger");
    const auto& L = o.ledger();
    chk::put("ledger_steps", double(steps)); chk::put("ledger_max_T", maxT); chk::put("ledger_seconds", elapsed);
    chk::put("ledger_joule_J", pt_to_double(L.joule, PT_E_INV)); chk::put("ledger_loads_J", pt_to_double(L.loads, PT_E_INV));
    chk::put("ledger_deposit_J", pt_to_double(L.deposit_E, PT_E_INV)); chk::put("ledger_evap_reference_J", pt_to_double(L.evap_ref, PT_E_INV));
    chk::put("ledger_evap_kinetic_J", pt_to_double(L.evap_kin, PT_E_INV));
    for (int el = 0; el < 3; ++el) chk::put(std::string("ledger_evaporated_atoms_") + element_name[el], pt_to_double(L.evap_N[el], PT_N_INV));
    chk::put("ledger_evaporated_kg_net", evap_kg);
    double maxdepth = 0; for (const auto& p : rep) maxdepth = std::max(maxdepth, p.melt_depth_m);
    chk::put("ledger_max_melt_depth_m", maxdepth);
    chk::put_text("ledger_state_sha256", o.state_sha256());
}

// ------------------------------------------------------------------------------------------------ quantized
static void quantized() {
    const double W = 0.002, H = 0.0002, dz = 0.0002;
    Config c = make_cfg(1, W, H, dz, 2, 1, 2, graded_rows(H, 0.5e-6, 1.4, 20e-6, 1));
    c.evap_mode = 2;
    c.patches = {{0, 1, 0.0, 0.5}, {0, 1, 0.5, 1.0}};
    Owner o(c);
    o.initialize_temperature([](int, double, double) { return 3800.0; });
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> U(0, 1);
    const double dt = 1e-12, w0 = 1e6;
    long released_macros = 0; bool exact = true;
    for (int s = 0; s < 400; ++s) {
        const Quote q = o.quote(dt);
        StepInput in; in.dt = dt; in.em_power_W.assign(o.dims().n_em, 4e5); in.quote_sha256 = q.sha256;
        for (int p = 0; p < o.n_patches(); ++p) for (int el = 0; el < 3; ++el) {
            const double D = q.expected_atoms[p][el], f = D / w0, fl = std::floor(f);
            const std::uint64_t n = std::uint64_t(fl) + (U(rng) < f - fl ? 1 : 0);
            if (n == 0) continue;
            in.releases.push_back({p, el, n, w0, double(n) * w0 * 2 * boltzmann * q.mean_T[p]});
            released_macros += long(n);
        }
        if (s == 3) {   // binding checks on a copy
            StepInput bad = in; bad.quote_sha256 = std::string(64, '0');
            bool r1 = false; try { o.prepare(bad); } catch (const PhaseRefusal& e) { r1 = e.kind == "release-binding"; }
            chk::check(r1, "release with a foreign quote digest is refused");
            if (!in.releases.empty()) {
                StepInput b2 = in; b2.releases[0].macro_count += 2;
                bool r2 = false; try { o.prepare(b2); } catch (const PhaseRefusal& e) { r2 = e.kind == "release-binding"; }
                chk::check(r2, "macro count outside {floor, ceil}(D/w0) is refused");
                StepInput b3 = in; b3.releases.push_back(b3.releases[0]);
                bool r3 = false; try { o.prepare(b3); } catch (const PhaseRefusal& e) { r3 = e.kind == "release-binding"; }
                chk::check(r3, "duplicate release is refused");
            }
        }
        o.advance(in);
        if (s % 50 == 0) exact = exact && o.ledger_closes();
    }
    exact = exact && o.ledger_closes();
    double rel = 0; for (int el = 0; el < 3; ++el) rel += pt_to_double(o.ledger().rel_N[el], PT_N_INV);
    chk::check(exact, "exact ledgers with quantized releases");
    chk::check(rel == double(released_macros) * w0, "released atoms equal n * w0 exactly");
    chk::put("quantized_released_macros", double(released_macros)); chk::put("quantized_release_energy_J", pt_to_double(o.ledger().release_E, PT_E_INV));
}

// ------------------------------------------------------------------------------------------------ merge
static void merge_run() {
    const double H = 0.0002;
    std::vector<double> rows{0.0, 20e-9, 40e-9, 80e-9, 160e-9, 320e-9, 640e-9, 1.28e-6, 2.56e-6, 5.12e-6, 10.24e-6, 20.48e-6, 40.96e-6, 81.92e-6, H};
    Config c = make_cfg(1, 1e-4, H, 1e-4, 1, 1, 1, rows);
    c.evap_face = {1, 0, 0, 0};   // evaporation from the LowY face only
    Owner o(c);
    o.initialize_temperature([](int, double, double) { return 4990.0; });   // evaporative cooling only
    const double dt = 1e-12; int merges = 0; long steps = 0;
    while (steps < 150000 && merges < 2) {
        StepInput in; in.dt = dt; in.em_power_W = {0.0};
        StepRecord r;
        try { r = o.advance(in); } catch (const PhaseRefusal& e) { chk::put_text("merge_stop", e.what()); break; }
        merges += r.merges; ++steps;
    }
    const auto rep = o.patch_reports();
    chk::check(merges >= 1, "depleted surface cell merged into its inward neighbour");
    chk::check(o.ledger_closes(), "exact ledgers across merges");
    const auto st = o.state();
    chk::check(st.flags[0] & 1, "merged cell is void and exposure moved inward");
    chk::put("merge_steps", double(steps)); chk::put("merge_count", merges); chk::put("merge_recession_m", rep[0].recession_m);
    chk::put("merge_surface_T", rep[0].mean_T);
}

// ------------------------------------------------------------------------------------------------ Bennett
// Bennett et al., Phys. Rev. Accel. Beams 26, 040401 (2023), Eq. (1): for a linearly rising surface field on a
// conductor of constant sigma with no heat conduction, the surface temperature rise is 1.273 mu0 H^2/(2 c_v);
// with c_v = 3.9 J/cm3/K stainless steel "melts near 1700 K, or 104 T". The owner is driven by the analytic
// magnetic-diffusion Joule power (H = H_s(t) 4 i2erfc(xi)) row by row.
static void bennett() {
    const double mu0 = 1.25663706212e-6, sigma = 1.3e6, trise = 100e-9, Bend = 150.0, Hdot = Bend / mu0 / trise;
    const double Hslab = 1.5e-3; const int rows = 1500;
    auto jcur = [&](double y, double t) { if (t <= 0) return 0.0; const double d = std::sqrt(t / (mu0 * sigma)); return 2 * Hdot * t / d * ierfc(y / (2 * d)); };
    auto run = [&](bool constant, std::vector<std::pair<double, double>>* trace) {
        Config c = make_cfg(1, 1e-4, Hslab, 1e-4, 1, rows, 1, uniform_rows(Hslab, rows));
        c.evap_mode = 0;
        if (constant) c.law = constant_law(3.9e6 / 8000, 1e-9, 1.9e4, 1, 0, 200, 2e4);
        Owner o(c);
        const int steps = 2000; const double dt = trise / steps;
        for (int s = 0; s < steps; ++s) {
            StepInput in; in.dt = dt; in.em_power_W.resize(rows);
            const double tm = (s + 0.5) * dt;
            for (int j = 0; j < rows; ++j) { const double y = (j + 0.5) * Hslab / rows, jj = jcur(y, tm); in.em_power_W[j] = jj * jj / sigma * (1e-4 * Hslab / rows * 1e-4); }
            try { o.advance(in); } catch (const PhaseRefusal& e) { chk::put_text("bennett_stop", e.what()); break; }
            if (trace) trace->push_back({mu0 * Hdot * (s + 1) * dt, o.temperatures()[0]});
        }
        return o.temperatures();
    };
    // constant c_v, no conduction: compare with Eq. (1) and the exact cell-centre integral
    const auto Tc = run(true, nullptr);
    const double y0 = 0.5 * Hslab / rows;
    double integral = 0; const int nq = 200000;
    for (int i = 0; i < nq; ++i) { const double t = (i + 0.5) * trise / nq, jj = jcur(y0, t); integral += jj * jj / sigma * trise / nq; }
    const double dT_cell = integral / 3.9e6, dT_eq1 = 1.273 * mu0 * (Hdot * trise) * (Hdot * trise) / (2 * 3.9e6);
    chk::put("bennett_eq1_surface_dT_K", dT_eq1); chk::put("bennett_cell_centre_exact_dT_K", dT_cell); chk::put("bennett_owner_dT_K", Tc[0] - 300);
    chk::check(std::abs(Tc[0] - 300 - dT_cell) / dT_cell < 1e-3, "owner reproduces the analytic cell-centre Joule heating (constant c_v)");
    chk::check(std::abs(Tc[0] - 300 - dT_eq1) / dT_eq1 < 0.01, "owner reproduces Bennett 2023 Eq. (1) at the surface cell");
    const double B_melt_eq1 = std::sqrt((1700 - 300) * 2 * 3.9e6 * mu0 / 1.273);
    chk::put("bennett_eq1_B_at_1700K_T", B_melt_eq1);
    // SS304 phase law (T-dependent Cp, latent heat, conduction), same constant-sigma drive
    std::vector<std::pair<double, double>> tr;
    run(false, &tr);
    double B_sol = 0, B_liq = 0, B_3000 = 0, B_eq1_sol = std::sqrt((1672 - 300) * 2 * 3.9e6 * mu0 / 1.273);
    for (const auto& [B, T] : tr) { if (!B_sol && T >= 1672) B_sol = B; if (!B_liq && T >= 1727) B_liq = B; if (!B_3000 && T >= 3000) B_3000 = B; }
    chk::put("bennett_ss304_B_at_solidus_T", B_sol); chk::put("bennett_ss304_B_at_liquidus_T", B_liq); chk::put("bennett_ss304_B_at_3000K_T", B_3000);
    chk::put("bennett_eq1_B_at_solidus_T", B_eq1_sol);
    chk::check(B_sol > 0 && std::abs(B_sol - 104) / 104 < 0.15, "SS304 surface reaches the solidus within 15% of the published 104 T");
}

// ------------------------------------------------------------------------------------------------ checkpoint
static void checkpoint(const std::string& out) {
    const double W = 0.002, H = 0.0002;
    Config c = make_cfg(1, W, H, 2e-4, 2, 1, 1, graded_rows(H, 1e-6, 1.5, 20e-6, 1));
    Owner a(c), b(c);
    a.initialize_temperature([](int, double x, double) { return 2500 + 1e5 * x; });
    b.initialize_temperature([](int, double x, double) { return 2500 + 1e5 * x; });
    auto in = [&](int s) { StepInput i; i.dt = 1e-11; i.em_power_W = {3e6, 2e6}; i.loads.push_back({0, -0.0005, H, 1e-8, "ck", "wall", std::uint64_t(s + 1)}); return i; };
    for (int s = 0; s < 30; ++s) { a.advance(in(s)); b.advance(in(s)); }
    const std::string path = out + "/checkpoint-test.bin";
    std::filesystem::remove(path);
    a.save(path);
    Owner r(c);
    r.load(path);
    chk::check(r.state_sha256() == a.state_sha256(), "restored state digest equals the saved state");
    for (int s = 30; s < 60; ++s) { a.advance(in(s)); r.advance(in(s)); }
    chk::check(r.state_sha256() == a.state_sha256(), "restored owner continues bit-identically");
    bool stale = false; try { auto x = in(5); r.prepare(x); } catch (const PhaseRefusal& e) { stale = e.kind == "load-identity"; }
    chk::check(stale, "restored producer watermarks refuse a replayed sequence");
    std::filesystem::remove(path);
}

// ------------------------------------------------------------------------------------------------ refusals
static void refusals() {
    const double W = 0.002, H = 0.0002;
    Config c = make_cfg(1, W, H, 2e-4, 2, 1, 1, uniform_rows(H, 4));
    Owner o(c);
    const std::string s0 = o.state_sha256();
    auto expect = [&](const StepInput& in, const std::string& kind, const std::string& what) {
        bool ok = false; std::string got;
        try { o.prepare(in); } catch (const PhaseRefusal& e) { got = e.kind; ok = e.kind == kind; }
        chk::check(ok && o.state_sha256() == s0, what + " (got " + got + ")");
    };
    StepInput a; a.dt = 0; a.em_power_W = {0, 0}; expect(a, "invalid-step", "zero dt refused");
    StepInput b; b.dt = 1e-12; b.em_power_W = {NAN, 0}; expect(b, "invalid-input", "nonfinite Joule refused");
    StepInput d; d.dt = 1e-12; d.em_power_W = {0, 0}; d.loads = {{0, 0.0, 0.0001, 1e-9, "p", "ch", 1}}; expect(d, "off-surface", "interior load point refused");
    StepInput e; e.dt = 1e-12; e.em_power_W = {0, 0}; e.loads = {{0, 0.0, 0.0, 1e-9, "p", "ch", 1}, {0, 0.0, 0.0, 1e-9, "p", "ch", 1}}; expect(e, "load-identity", "duplicate producer sequence refused");
    StepInput f; f.dt = 1e-12; f.em_power_W = {0, 0}; f.loads = {{0, 0.0, 0.0, -1e3, "p", "ch", 1}}; expect(f, "temperature-below-domain", "cooling below 298.2 K refused");
    StepInput g; g.dt = 1.0; g.em_power_W = {0, 0}; expect(g, "substep-budget", "step beyond the sub-step budget refused");
    StepInput h; h.dt = 1e-12; h.em_power_W = {0, 0}; h.releases = {{0, 0, 1, 1e6, 0}}; expect(h, "release-binding", "release outside quantized mode refused");
    // a prepared trial is invalidated by the commit of another trial
    StepInput ok; ok.dt = 1e-12; ok.em_power_W = {1e3, 1e3};
    Trial t1 = o.prepare(ok), t2 = o.prepare(ok);
    chk::check(o.commit(std::move(t2)), "second trial commits");
    chk::check(!o.preflight(t1) && !o.commit(std::move(t1)), "stale trial fails preflight after another commit");
}

// ------------------------------------------------------------------------------------------------ super-cycling
// The same Joule energy delivered as 100 steps of power or as 10 steps of accumulated energy: both ledgers exact,
// equal admitted Joule words, temperatures agree to the time-discretization level.
static void supercycle() {
    const double W = 0.002, H = 0.0002;
    Config c = make_cfg(1, W, H, 2e-4, 2, 1, 1, graded_rows(H, 1e-6, 1.5, 20e-6, 1));
    c.evap_mode = 0;
    Owner a(c), b(c);
    const std::vector<double> P = {3e6, 1e6};
    const double dt = 1e-11;
    for (int s = 0; s < 100; ++s) { StepInput in; in.dt = dt; in.em_power_W = P; a.advance(in); }
    for (int s = 0; s < 10; ++s) {
        StepInput in; in.dt = 10 * dt; in.em_energy_J.assign(2, 0.0);
        for (int e = 0; e < 2; ++e) for (int k = 0; k < 10; ++k) in.em_energy_J[e] += P[e] * dt;
        b.advance(in);
    }
    const auto Ta = a.temperatures(), Tb = b.temperatures();
    double d = 0; for (std::size_t i = 0; i < Ta.size(); ++i) d = std::max(d, std::abs(Ta[i] - Tb[i]));
    chk::check(a.ledger_closes() && b.ledger_closes(), "exact ledgers with power and with accumulated-energy input");
    chk::check(std::abs(pt_to_double(a.ledger().joule, PT_E_INV) / pt_to_double(b.ledger().joule, PT_E_INV) - 1) < 1e-14, "same admitted Joule energy");
    chk::put("supercycle_max_dT_K", d);
    chk::check(d < 0.5, "super-cycled temperatures within 0.5 K");
    StepInput both; both.dt = dt; both.em_power_W = P; both.em_energy_J = P;
    bool r = false; try { a.prepare(both); } catch (const PhaseRefusal& e) { r = e.kind == "invalid-input"; }
    chk::check(r, "power and energy vectors together are refused");
}

// ------------------------------------------------------------------------------------------------ hall geometry
static void hall_geometry() {
    Config c;
    for (int body = 0; body < 2; ++body) for (int k = 0; k < 4; ++k) c.slices.push_back({body, k, -0.0025, body ? 0.0025 : -0.0027});
    c.width = 0.005; c.thickness = 0.0002; c.height = 0.0002; c.nx_em = 25; c.ny_em = 1; c.rx = 1;
    c.y_edges = graded_rows(0.0002, 0.25e-6, 1.3, 10e-6, 1);
    Owner o(c);
    chk::put("hall_rows", double(o.dims().ny)); chk::put("hall_cells", double(o.dims().n_cells)); chk::put("hall_tau_max_s", o.stability_tau());
    chk::put("hall_substeps_at_10fs", double(o.substeps_for(1e-14)));
    const auto sig = o.sigma_em();
    double rho, k; pt_ss_transport(300.0, &rho, &k);
    double dev = 0; for (double s : sig) dev = std::max(dev, std::abs(s * rho - 1));
    chk::check(dev < 1e-15, "uniform 300 K: sigma per EM cell equals 1/rho(300 K)");
    std::vector<double> P(o.dims().n_em, 1e5);
    const auto t0 = Clock::now();
    for (int s = 0; s < 200; ++s) { StepInput in; in.dt = 1e-14; in.em_power_W = P; o.advance(in); }
    chk::put("hall_cpu_ms_per_step", 1e3 * secs(t0) / 200);
    chk::check(o.ledger_closes(), "exact ledger closure (Hall geometry)");
}

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : ".";
    std::set<std::string> want;
    for (int i = 2; i < argc; ++i) want.insert(argv[i]);
    auto on = [&](const char* s) { return want.empty() || want.count(s); };
    auto timed = [&](const char* s, const std::function<void()>& f) { if (!on(s)) return; const auto t = Clock::now(); f(); chk::put(std::string("seconds_") + s, secs(t)); std::printf("suite %s done (%.1f s)\n", s, secs(t)); std::fflush(stdout); };
    timed("uniform", uniform_heating);
    timed("conduction", conduction);
    timed("mode", mode_decay);
    timed("stefan", stefan);
    timed("refusals", refusals);
    timed("checkpoint", [&] { checkpoint(out); });
    timed("quantized", quantized);
    timed("supercycle", supercycle);
    timed("hall", hall_geometry);
    timed("bennett", bennett);
    timed("merge", merge_run);
    timed("ledger", ledger_run);
    std::string name = "test_owner";
    if (!want.empty()) for (const auto& w : want) name += "-" + w;
    return chk::finish(name, out);
}
