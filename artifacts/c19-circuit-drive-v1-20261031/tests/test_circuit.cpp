// Circuit-only tests of c19-circuit-drive-v1 (no Hall dependency).
//  T1 ledger closure, full pulse with the reduced Hall vacuum load (dt = 100 ps and 10 ps) and a short dt = 10 fs window
//  T2 implicit midpoint accuracy: series RLC into a resistor vs the analytic solution, second-order convergence
//  T3 port coupling on the reduced load (affine port): evaluations per step and final |F|
//  T4 open switch before t_trigger: no current, stored energy constant, ledger exact
//  T5 determinism: two identical runs give identical state bytes
// Usage: test_circuit PARAMS_FILE LAW_TABLE
#include "../src/c19_circuit.hpp"
#include "../src/reduced_hall_load.hpp"
#include <cstdio>
#include <cstring>
#include <map>

using namespace c19_circuit_v1;

static int failures = 0;
static void check(bool ok, const char* name, const std::string& detail) {
    std::printf("%s %s %s\n", ok ? "PASS" : "FAIL", name, detail.c_str());
    if (!ok) ++failures;
}
static std::string fmt(const char* f, double a, double b = 0, double c = 0, double d = 0) {
    char buf[512]; std::snprintf(buf, sizeof buf, f, a, b, c, d); return buf;
}

static Parameters read_params(const char* path, double& L_H) {
    std::map<std::string, double> kv; std::ifstream in(path); std::string line;
    while (std::getline(in, line)) { if (line.empty() || line[0] == '#') continue; std::istringstream r(line); std::string k; double v; if (r >> k >> v) kv[k] = v; }
    Parameters p; auto g = [&](const char* k, double d) { auto it = kv.find(k); return it == kv.end() ? d : it->second; };
    p.C_F = g("C_F", p.C_F); p.V0_V = g("V0_V", p.V0_V); p.L_ltd_H = g("L_ltd_H", p.L_ltd_H); p.L_feed_H = g("L_feed_H", p.L_feed_H);
    p.R_inf_ohm = g("R_inf_ohm", 0); p.R_sw_ohm = g("R_sw_ohm", 0); p.tau_sw_s = g("tau_sw_s", 1e-9);
    p.inv_L_core_per_H = g("inv_L_core_per_H", 0); p.G_core_S = g("G_core_S", 0); p.R_feed_ohm = g("R_feed_ohm", 0);
    L_H = g("L_H", 1.1045115428330026e-08);
    return p;
}

struct RunResult { double closure_over_E0, max_res_over_bound, load_closure_rel; int max_eval; double max_F_rel; std::vector<double> history; };

static RunResult run_reduced(Parameters p, double L_H, const LawTable& law, double dt, double t_end, bool keep_history) {
    LumpedLtdCircuit c(p, 0.); ReducedHallLoad load; load.law = &law; load.L_H = L_H;
    PortCoupler coupler; RunResult rr{}; rr.max_eval = 0;
    while (c.state().time_s < t_end * (1 - 1e-12)) {
        const double tm = c.state().time_s + .5 * dt;
        auto r = coupler.step(c, dt, [&](double I) { return load.voltage(dt, tm, I); });
        rr.max_eval = std::max(rr.max_eval, r.evaluations);
        rr.max_F_rel = std::max(rr.max_F_rel, std::abs(r.F) / std::max(1., std::abs(r.I)));
        auto sl = c.accept(dt, r.I, r.V_T, dt * r.I * r.V_T);
        rr.max_res_over_bound = std::max(rr.max_res_over_bound, std::abs(sl.residual_J) / sl.rounding_bound_J);
        load.accept(dt, tm, r.I, r.V_T);
        if (keep_history) { const auto& s = c.state(); for (double v : {s.time_s, s.v_C, s.i_s, s.i_m, s.i_f}) rr.history.push_back(v); }
    }
    const auto& L = c.ledger(); const double E0 = L.initial_stored_J, E1 = c.stored_energy(c.state());
    rr.closure_over_E0 = (E1 - E0 + L.dissipated_J.value() + L.delivered_J.value()) / E0;
    rr.load_closure_rel = (load.received_J.value() - load.joule_J.value() - load.magnetic_J.value()) / load.received_J.value();
    return rr;
}

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: test_circuit PARAMS LAW\n"); return 2; }
    double L_H; const Parameters p = read_params(argv[1], L_H); LawTable law(argv[2]);

    // T1
    for (double dt : {1e-10, 1e-11}) {
        auto r = run_reduced(p, L_H, law, dt, 260e-9, false);
        // The circuit ledger must close at the rounding level of each step (residual/bound <= 1); the cumulative closure is
        // bounded by steps * 16 eps relative to the stored energy, the reduced-load fixture by the same per-step rule.
        const double steps = 260e-9 / dt, cum = steps * 16 * eps;
        check(std::abs(r.closure_over_E0) < cum && r.max_res_over_bound <= 1 && std::abs(r.load_closure_rel) < cum,
              "T1_ledger_full_pulse", fmt("dt=%.0e closure/E0=%.3e max_step_residual/bound=%.3f load_closure_rel=%.3e", dt, r.closure_over_E0,
                                          r.max_res_over_bound, r.load_closure_rel) + fmt(" (cumulative bound %.2e)", cum));
    }
    {   // dt = 10 fs (the Hall step) over 20 ns from closure
        auto r = run_reduced(p, L_H, law, 1e-14, 20e-9, false);
        check(std::abs(r.closure_over_E0) < 2e6 * 16 * eps && r.max_res_over_bound <= 1, "T1_ledger_hall_step_dt",
              fmt("dt=1e-14 steps=2e6 closure/E0=%.3e max_step_residual/bound=%.3f max_evaluations=%.0f", r.closure_over_E0, r.max_res_over_bound, r.max_eval));
    }
    // T2: C=144 nF, V0=800 kV, L=36 nH, R=0.2 ohm series, port resistor 0.5 ohm; feed 1 pH. Analytic underdamped series RLC.
    {
        Parameters q; q.L_feed_H = 1e-12; q.R_inf_ohm = .2; q.R_sw_ohm = 0; q.G_core_S = 0; q.inv_L_core_per_H = 0;
        const double R = .2 + .5, L = q.L_ltd_H + q.L_feed_H, C = q.C_F, a = R / (2 * L), w0 = 1 / std::sqrt(L * C), wd = std::sqrt(w0 * w0 - a * a);
        double err[3]; const double dts[3] = {1e-10, 5e-11, 2.5e-11};
        for (int k = 0; k < 3; ++k) {
            LumpedLtdCircuit c(q, 0.); PortCoupler cp; double maxe = 0; const double dt = dts[k];
            while (c.state().time_s < 300e-9 * (1 - 1e-12)) {
                auto r = cp.step(c, dt, [&](double I) { return .5 * I; });
                c.accept(dt, r.I, r.V_T, dt * r.I * r.V_T);
                const double t = c.state().time_s, exact = q.V0_V / (wd * L) * std::exp(-a * t) * std::sin(wd * t);
                maxe = std::max(maxe, std::abs(c.state().i_f - exact));
            }
            err[k] = maxe;
        }
        const double o1 = std::log2(err[0] / err[1]), o2 = std::log2(err[1] / err[2]);
        check(o1 > 1.9 && o2 > 1.9 && err[0] < 2e3, "T2_midpoint_order_RLC",
              fmt("max|I-I_exact| = %.3e, %.3e, %.3e A (dt 100/50/25 ps)", err[0], err[1], err[2]) + fmt("; observed order %.3f / %.3f", o1, o2));
    }
    // T3
    {
        auto r = run_reduced(p, L_H, law, 1e-11, 260e-9, false);
        check(r.max_eval <= 3 && r.max_F_rel < 1e-13, "T3_port_coupling_reduced_load", fmt("max evaluations per step=%.0f max|F|/|I|=%.3e", r.max_eval, r.max_F_rel));
    }
    // T4
    {
        Parameters q = p; q.t_trigger_s = 5e-9; LumpedLtdCircuit c(q, 0.); PortCoupler cp; double maxI = 0;
        for (int n = 0; n < 400; ++n) { auto r = cp.step(c, 1e-11, [&](double I) { return 1e-3 * I; }); c.accept(1e-11, r.I, r.V_T, 1e-11 * r.I * r.V_T); maxI = std::max(maxI, std::abs(c.state().i_f)); }
        const double E = c.stored_energy(c.state());
        check(maxI == 0 && E == c.ledger().initial_stored_J && c.ledger().dissipated_J.value() == 0, "T4_open_switch",
              fmt("max|I| before trigger=%.3e A stored change=%.3e J", maxI, E - c.ledger().initial_stored_J));
    }
    // T5
    {
        auto a = run_reduced(p, L_H, law, 1e-10, 260e-9, true), b = run_reduced(p, L_H, law, 1e-10, 260e-9, true);
        check(a.history.size() == b.history.size() && std::memcmp(a.history.data(), b.history.data(), a.history.size() * sizeof(double)) == 0,
              "T5_determinism", fmt("%.0f state words compared", double(a.history.size())));
    }
    std::printf("%s %d failure(s)\n", failures ? "FAILED" : "ALL PASSED", failures);
    return failures ? 1 : 0;
}
