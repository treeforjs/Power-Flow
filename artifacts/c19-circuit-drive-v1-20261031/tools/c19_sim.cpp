// c19_sim: the lumped Mykonos LTD circuit (src/c19_circuit.hpp) coupled through PortCoupler to the reduced Hall vacuum
// load (src/reduced_hall_load.hpp). Used by the calibration (calibration/calibrate.py) and the scenario runs, so the
// calibrated model is the same code that runs inside the coupled step.
// Usage: c19_sim PARAMS_FILE LAW_TABLE_CSV OUT_CSV [every]
// PARAMS_FILE: "key value" lines; keys = Parameters fields, plus dt_s, t_end_s, L_H, t_breakdown_s, R_plasma_ohm,
// resistive_load_ohm (if > 0 the port load is a pure resistor instead of the Hall load: rating check).
#include "../src/c19_circuit.hpp"
#include "../src/reduced_hall_load.hpp"
#include <cstdio>
#include <map>

using namespace c19_circuit_v1;

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: c19_sim PARAMS LAW OUT [every]\n"); return 2; }
    std::map<std::string, double> kv;
    { std::ifstream in(argv[1]); std::string k; double v; while (in >> k >> v) kv[k] = v; }
    auto get = [&](const char* k, double d) { auto it = kv.find(k); return it == kv.end() ? d : it->second; };
    Parameters p;
    p.C_F = get("C_F", p.C_F); p.V0_V = get("V0_V", p.V0_V); p.L_ltd_H = get("L_ltd_H", p.L_ltd_H); p.L_feed_H = get("L_feed_H", p.L_feed_H);
    p.R_inf_ohm = get("R_inf_ohm", 0); p.R_sw_ohm = get("R_sw_ohm", 0); p.tau_sw_s = get("tau_sw_s", 1e-9);
    p.inv_L_core_per_H = get("inv_L_core_per_H", 0); p.G_core_S = get("G_core_S", 0); p.R_feed_ohm = get("R_feed_ohm", 0);
    p.t_trigger_s = 0;
    const double dt = get("dt_s", 1e-10), t_end = get("t_end_s", 260e-9), Rres = get("resistive_load_ohm", 0);
    const long every = argc > 4 ? std::atol(argv[4]) : 1;
    LawTable law(argv[2]);
    ReducedHallLoad load; load.law = &law; load.L_H = get("L_H", load.L_H);
    load.t_breakdown_s = get("t_breakdown_s", std::numeric_limits<double>::infinity()); load.R_plasma_ohm = get("R_plasma_ohm", 0);
    LumpedLtdCircuit c(p, 0.);
    CouplerOptions o; o.max_evaluations = 30; PortCoupler coupler(o);
    std::FILE* f = std::fopen(argv[3], "w");
    std::fprintf(f, "t_s,I_port_A,V_port_V,v_A_V,v_C_V,i_s_A,i_m_A,T_foil_K,R_load_ohm,evaluations\n");
    std::fprintf(f, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n", 0., 0., 0., 0., c.state().v_C, 0., 0., load.T, Rres > 0 ? Rres : load.R_now(0), 0);
    long n = 0; int max_eval = 0; double max_res_over_bound = 0;
    while (c.state().time_s < t_end * (1 - 1e-12)) {
        const double tm = c.state().time_s + .5 * dt;
        auto port = [&](double I) { return Rres > 0 ? Rres * I : load.voltage(dt, tm, I); };
        auto r = coupler.step(c, dt, port);
        if (!r.converged) {
            std::fprintf(stderr, "port coupling did not converge: t=%.6e I=%.17g F=%.3e a=%.17g Y=%.6e Z=%.6e evaluations=%d\n",
                         tm, r.I, r.F, r.a, r.Y, r.Z, r.evaluations);
            return 3;
        }
        max_eval = std::max(max_eval, r.evaluations);
        const double work = dt * r.I * r.V_T;
        auto sl = c.accept(dt, r.I, r.V_T, work);
        max_res_over_bound = std::max(max_res_over_bound, std::abs(sl.residual_J) / sl.rounding_bound_J);
        // node A voltage at the midpoint (water-line output / dI/dt monitor)
        if (Rres <= 0) load.accept(dt, tm, r.I, r.V_T);
        ++n;
        if (n % every == 0) {
            const auto& s = c.state();
            // midpoint quantities: port current and voltage of the step; state quantities at the end of the step.
            const double vA = r.V_T + p.L_feed_H * (s.i_f - (2 * r.I - s.i_f)) / dt + p.R_feed_ohm * r.I;
            std::fprintf(f, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n", tm, r.I, r.V_T, vA, s.v_C, s.i_s, s.i_m, load.T,
                         Rres > 0 ? Rres : load.R_now(tm), r.evaluations);
        }
    }
    std::fclose(f);
    const auto& L = c.ledger();
    const double E0 = L.initial_stored_J, E1 = c.stored_energy(c.state());
    std::printf("steps %ld max_evaluations %d\n", n, max_eval);
    std::printf("ledger initial_stored_J %.17g final_stored_J %.17g dissipated_J %.17g delivered_J %.17g\n", E0, E1, L.dissipated_J.value(), L.delivered_J.value());
    std::printf("ledger closure_J %.6e closure_over_E0 %.3e increment_sum_minus_difference_J %.6e max_step_residual_over_bound %.3f\n",
                E1 - E0 + L.dissipated_J.value() + L.delivered_J.value(), (E1 - E0 + L.dissipated_J.value() + L.delivered_J.value()) / E0,
                L.dE_J.value() - (E1 - E0), max_res_over_bound);
    std::printf("join port_iteration_defect_J %.6e hall_dot_rounding_J %.6e\n", L.port_iteration_defect_J.value(), L.hall_dot_rounding_J.value());
    if (Rres <= 0)
        std::printf("load joule_J %.17g foil_joule_J %.17g magnetic_J %.17g received_J %.17g load_closure_J %.3e T_final_K %.3f\n", load.joule_J.value(),
                    load.foil_joule_J.value(), load.magnetic_J.value(), load.received_J.value(),
                    load.received_J.value() - load.joule_J.value() - load.magnetic_J.value(), load.T);
    if (Rres <= 0) std::printf("load t_law_exceeded_s %.6e (sigma held at 1500 K beyond)\n", load.t_law_exceeded_s);
    return 0;
}
