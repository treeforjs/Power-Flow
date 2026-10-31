// Hall-port tests of c19-circuit-drive-v1. Read-only includes: Hall corrective v3.1 field_contract.hpp (the C19
// CurrentDrive of record) and the sealed Hall axial field step (axial::Fields: the Hall operator with an injectable
// impressed current; dt <= 10 fs as in Hall).
//  R1 regression: PrescribedCurrentDrive == hall_corrective_v3::CurrentDrive bit for bit (knots, 1e6 random times,
//     1e6 Hall step midpoints generated with Hall's own time recurrence).
//  R2 regression: N Hall field steps driven by the Hall CurrentDrive vs by the circuit module's prescribed mode:
//     every receipt scalar and the final E, B, E_low, B_low arrays bit for bit; prescribed ledger = sum circuit_work.
//  H1 coupled circuit + Hall field step (vacuum): per-step Newton/secant iteration, evaluations, |F|, joined ledger.
//  H2 coupled with an iterate-dependent gap conduction current switched on mid-run (breakdown mimic): joint fixed
//     point of the port current and the conduction current, joined ledger including plasma work.
// Usage: test_hall_port TRACE_CSV PARAMS LAW_TABLE N_R2 N_H1 N_H2 [OUT_CSV_PREFIX]
#include "field_contract.hpp"
#include "../src/c19_circuit.hpp"
#include "../src/reduced_hall_load.hpp"
#include <cstdio>
#include <cstring>
#include <map>
#include <random>

using namespace c19_circuit_v1;

static int failures = 0;
static void check(bool ok, const char* name, const std::string& detail) {
    std::printf("%s %s %s\n", ok ? "PASS" : "FAIL", name, detail.c_str()); std::fflush(stdout);
    if (!ok) ++failures;
}
static std::string fmt(const char* f, double a, double b = 0, double c = 0, double d = 0, double e = 0) {
    char buf[1024]; std::snprintf(buf, sizeof buf, f, a, b, c, d, e); return buf;
}
static bool same(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }
static bool same(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
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

struct HallSums { Sum2 field_delta, plasma_work, joule, outer_loss, outer_supply, circuit_work, balance; double max_balance_rel = 0; };
static void add(HallSums& h, const axial::Receipt& r) {
    h.field_delta.add(r.field_delta_J); h.plasma_work.add(r.plasma_work_J); h.joule.add(r.joule_J); h.outer_loss.add(r.outer_loss_J);
    h.outer_supply.add(r.outer_supply_J); h.circuit_work.add(r.circuit_work_J); h.balance.add(r.balance_J);
    const double act = std::abs(r.field_delta_J) + std::abs(r.plasma_work_J) + std::abs(r.joule_J) + std::abs(r.outer_loss_J) + std::abs(r.outer_supply_J) + std::abs(r.circuit_work_J);
    h.max_balance_rel = std::max(h.max_balance_rel, std::abs(r.balance_J) / std::max(act, 1e-300));
}

int main(int argc, char** argv) {
    if (argc < 7) { std::fprintf(stderr, "usage\n"); return 2; }
    const std::filesystem::path trace = argv[1];
    double L_H; const Parameters p = read_params(argv[2], L_H); LawTable law(argv[3]);
    const long N2 = std::atol(argv[4]), NH1 = std::atol(argv[5]), NH2 = std::atol(argv[6]);
    const std::string prefix = argc > 7 ? argv[7] : "";
    const double dt = 1e-14;

    // R1 -------------------------------------------------------------------------------------------------------
    const hall_corrective_v3::CurrentDrive hall(trace);   // verifies the sha256 of the bytes of record
    const PrescribedCurrentDrive mine(trace);
    {
        long n = 0, bad = 0;
        for (const auto& k : mine.knots()) { ++n; bad += !same(hall.current(k.first), mine.current(k.first)); }
        std::mt19937_64 rng(20261031); std::uniform_real_distribution<double> u(hall.start(), hall.end());
        for (int i = 0; i < 1000000; ++i) { const double t = u(rng); ++n; bad += !same(hall.current(t), mine.current(t)); }
        double t = hall.start();
        for (long i = 0; i < 1000000; ++i) { ++n; bad += !same(hall.current(t + .5 * dt), mine.step_current(t, dt)); t = t + dt; }
        ++n; bad += !same(hall.current(hall.end()), mine.current(mine.end()));
        check(bad == 0 && same(hall.start(), mine.start()) && same(hall.end(), mine.end()), "R1_prescribed_drive_bitwise",
              fmt("%.0f queries, %.0f differ", double(n), double(bad)));
    }

    // Hall-like field step on the Hall mesh.
    axial::Mesh g(45, 4, .009, 1);
    const double sigma300 = law.sigma_at(300.);
    axial::Vec s300(g.cells.size()); for (std::size_t c = 0; c < s300.size(); ++c) if (g.cells[c].body >= 0) s300[c] = sigma300;
    const axial::Vec zero(g.ne);

    // R2 -------------------------------------------------------------------------------------------------------
    {
        axial::Fields f(g, s300, hall.current(hall.start()));
        auto a = f.start(hall.start()), b = a;
        PrescribedLedger led; HallSums ha; long mismatched = 0;
        for (long i = 0; i < N2; ++i) {
            const double Ia = hall.current(a.time_s + .5 * dt);          // literally the Hall Model expression
            const double Ib = mine.step_current(b.time_s, dt);           // the circuit module, prescribed mode
            auto [na, ra] = f.prepare(a, zero, s300, dt, Ia);
            auto [nb, rb] = f.prepare(b, zero, s300, dt, Ib);
            mismatched += !(same(Ia, Ib) && same(ra.terminal_V, rb.terminal_V) && same(ra.circuit_work_J, rb.circuit_work_J) &&
                            same(ra.field_delta_J, rb.field_delta_J) && same(ra.joule_J, rb.joule_J) && same(ra.balance_J, rb.balance_J));
            led.accept(dt, Ib, rb.circuit_work_J); add(ha, ra);
            a = std::move(na); b = std::move(nb);
        }
        const bool state_same = same(a.E, b.E) && same(a.B, b.B) && same(a.E_low, b.E_low) && same(a.B_low, b.B_low) && same(a.time_s, b.time_s);
        check(mismatched == 0 && state_same && same(led.source_supplied_J.value(), ha.circuit_work.value()), "R2_prescribed_mode_hall_steps_bitwise",
              fmt("%.0f Hall steps (dt 10 fs) from the first C19 knot: receipt mismatches %.0f, final E/B/E_low/B_low identical %.0f, "
                  "source ledger %.17g J = sum circuit_work %.17g J", double(N2), double(mismatched), double(state_same),
                  led.source_supplied_J.value(), ha.circuit_work.value()));
    }

    // H1/H2 setup: advance the circuit with the reduced Hall load to t_s, then hand over to the Hall field step.
    const double t_s = 40e-9;
    LumpedLtdCircuit circ(p, 0.);
    ReducedHallLoad red; red.law = &law; red.L_H = L_H;
    {
        PortCoupler cp; const double h = 1e-11;
        while (circ.state().time_s < t_s * (1 - 1e-12)) {
            const double tm = circ.state().time_s + .5 * h;
            auto r = cp.step(circ, h, [&](double I) { return red.voltage(h, tm, I); });
            circ.accept(h, r.I, r.V_T, h * r.I * r.V_T); red.accept(h, tm, r.I, r.V_T);
        }
    }
    const double sigmaT = law.sigma_at(std::min(red.T, 1500.));
    axial::Vec sT(g.cells.size()); for (std::size_t c = 0; c < sT.size(); ++c) if (g.cells[c].body >= 0) sT[c] = sigmaT;
    std::printf("INFO handover t=%.6e s I=%.6e A v_C=%.6e V foil T=%.1f K sigma=%.6e S/m\n", circ.state().time_s, circ.state().i_f, circ.state().v_C, red.T, sigmaT);

    auto coupled = [&](const char* name, long N, bool plasma_on, long plasma_start, double sigma_p, LumpedLtdCircuit c, PortCoupler& cp) {
        axial::Fields f(g, sT, c.state().i_f);
        auto st = f.start(c.state().time_s);
        // ledger snapshot at the handover: the coupled segment is measured by differences
        const double E0 = c.stored_energy(c.state()), D0 = c.ledger().dissipated_J.value(), P0 = c.ledger().delivered_J.value();
        const double PD0 = c.ledger().port_iteration_defect_J.value(), HR0 = c.ledger().hall_dot_rounding_J.value();
        const double dEsum0 = c.ledger().dE_J.value();
        // gap conduction edges: y-edges strictly between the foils, |x| <= 2.5 mm, within the foil length
        std::vector<int> gap;
        for (int e = 0; e < g.ne; ++e) {
            const auto& ed = g.edges[e]; if (ed.axis != 1) continue;
            const auto pa = g.node_position(ed.a), pb = g.node_position(ed.b);
            if (std::abs(pa[0]) <= .0025 + 1e-12 && pa[1] > -.0025 + 1e-12 && pb[1] < .0025 - 1e-12 && pa[2] >= -.0125 - 1e-12 && pa[2] <= .0125 + 1e-12) gap.push_back(e);
        }
        axial::Vec plasma(g.ne), next_plasma(g.ne);
        HallSums hs; Sum2 ec_res, segE, segD, segP, segDefect, segRound, segBound, hallAllowed; long total_eval = 0; int max_eval = 0; double maxF = 0, max_dp = 0; long unconverged = 0;
        std::FILE* out = prefix.empty() ? nullptr : std::fopen((prefix + name + ".csv").c_str(), "w");
        if (out) std::fprintf(out, "step,t_mid_s,I_port_A,V_T_V,evaluations,F_A,Z_ohm,plasma_work_J\n");
        double dp_rel = 0;
        for (long i = 0; i < N; ++i) {
            const bool on = plasma_on && i >= plasma_start;
            std::pair<axial::State, axial::Receipt> cand;
            auto evaluate = [&](double I) {
                cand = f.prepare(st, plasma, sT, dt, I);
                // next iterate of the conduction current from this trial field (fixed point, as the coordinator does with motion)
                if (on) {
                    double num = 0, den = 0;
                    for (int e : gap) { next_plasma[e] = sigma_p * g.VE[e] * cand.second.mid_E[e]; num = std::max(num, std::abs(next_plasma[e] - plasma[e])); den = std::max(den, std::abs(next_plasma[e])); }
                    dp_rel = den > 0 ? num / den : 0; std::swap(plasma, next_plasma);
                } else dp_rel = 0;
                return cand.second.terminal_V;
            };
            // the conduction current used by the last prepared candidate is the one before the swap; converged when
            // the change is below 1e-13 relative (then the swapped vector equals the used one to that level)
            auto others = [&] { return dp_rel <= 1e-13; };
            auto r = on ? cp.step(c, dt, evaluate, others) : cp.step(c, dt, evaluate);
            if (!r.converged) ++unconverged;
            total_eval += r.evaluations; max_eval = std::max(max_eval, r.evaluations); maxF = std::max(maxF, std::abs(r.F) / std::max(1., std::abs(r.I)));
            max_dp = std::max(max_dp, dp_rel);
            auto sl = c.accept(dt, r.I, r.V_T, cand.second.circuit_work_J);
            ec_res.add(sl.residual_J); segE.add(sl.dE_stored_J); segD.add(sl.dissipated_J); segP.add(sl.delivered_J);
            segDefect.add(sl.delivered_J - dt * r.I * r.V_T); segRound.add(dt * r.I * r.V_T - cand.second.circuit_work_J);
            segBound.add(sl.rounding_bound_J);
            {   const auto& q = cand.second;   // the Hall step's own admitted balance bound: 1e-20 + 1e-10 activity
                hallAllowed.add(1e-20 + 1e-10 * (std::abs(q.field_delta_J) + std::abs(q.plasma_work_J) + std::abs(q.joule_J) + std::abs(q.outer_loss_J) +
                                                 std::abs(q.outer_supply_J) + std::abs(q.circuit_work_J))); }
            add(hs, cand.second);
            if (out && (i % 10 == 0 || i == N - 1)) std::fprintf(out, "%ld,%.17g,%.17g,%.17g,%d,%.6e,%.6e,%.6e\n", i, st.time_s + .5 * dt, r.I, r.V_T, r.evaluations, r.F, r.Z, cand.second.plasma_work_J);
            st = std::move(cand.first);
        }
        if (out) std::fclose(out);
        const auto& L = c.ledger();
        const double E1 = c.stored_energy(c.state());
        const double dEdiff = E1 - E0, Dc = segD.value(), Pc = segP.value(), W = hs.circuit_work.value();
        (void)D0; (void)P0;
        // stored-energy change in increment form (sum of per-step C dv v_mid + L di i_mid terms); the difference of the
        // absolute stored energies (~4e4 J) is reported beside it and agrees to its own rounding (~eps E0 per evaluation)
        const double dEc = segE.value(); (void)dEsum0;
        const double circuit_closure = dEc + Dc + Pc;
        const double join = Pc - W;          // = port iteration defect + Hall dot rounding
        const double hall_closure = hs.field_delta.value() + hs.plasma_work.value() + hs.joule.value() + hs.outer_loss.value() - hs.outer_supply.value() - W;
        const double total = dEc + Dc + hs.field_delta.value() + hs.plasma_work.value() + hs.joule.value() + hs.outer_loss.value() - hs.outer_supply.value();
        const double activity = std::abs(Dc) + std::abs(Pc) + std::abs(hs.field_delta.value()) + std::abs(hs.joule.value()) + std::abs(hs.plasma_work.value()) + std::abs(hs.outer_loss.value());
        std::printf("INFO %s steps=%ld t_end=%.9e I_end=%.9e mean_evaluations=%.3f max_evaluations=%d max|F|/|I|=%.3e max_conduction_change=%.3e Z_est=%.6e unconverged=%ld\n",
                    name, N, st.time_s, c.state().i_f, double(total_eval) / N, max_eval, maxF, max_dp, cp.impedance_estimate(), unconverged);
        std::printf("INFO %s circuit: dE=%.12e D=%.12e P=%.12e closure=%.3e (sum of step residuals %.3e)\n", name, dEc, Dc, Pc, circuit_closure, ec_res.value());
        std::printf("INFO %s join: P_circuit - W_hall = %.3e J (port iteration defect %.3e, Hall dot rounding %.3e)\n", name, join,
                    segDefect.value(), segRound.value()); (void)PD0; (void)HR0;
        std::printf("INFO %s hall: field_delta=%.12e joule=%.12e plasma_work=%.12e outer_loss=%.12e outer_supply=%.12e circuit_work=%.12e closure=%.3e max step balance/activity=%.3e\n",
                    name, hs.field_delta.value(), hs.joule.value(), hs.plasma_work.value(), hs.outer_loss.value(), hs.outer_supply.value(), W, hall_closure, hs.max_balance_rel);
        std::printf("INFO %s total system: dE_circuit + D_circuit + field_delta + joule + plasma_work + outer_loss - outer_supply = %.3e J (activity %.3e J, relative %.3e)\n",
                    name, total, activity, total / activity);
        std::printf("INFO %s stored-energy change: increment form %.12e J, difference of absolute stored energies %.12e J (rounding of a %.3e J store ~%.1e J)\n",
                    name, dEc, dEdiff, E0, eps * E0);
        // circuit closure at the rounding level of the stored energy (it is the difference of two 4e4 J numbers)
        // Acceptance: the circuit closes within the sum of its per-step rounding bounds (dominated by the representation
        // of a ~4.5e4 J store in doubles), the join within 1e-12 of the delivered energy, Hall within its own admitted
        // per-step bounds, and the total system within the sum of these.
        std::printf("INFO %s bounds: circuit sum of step bounds %.3e J, Hall admitted %.3e J\n", name, segBound.value(), hallAllowed.value());
        const bool ok = unconverged == 0 && std::abs(circuit_closure) <= segBound.value() && std::abs(join) <= 1e-12 * std::abs(Pc) + 1e-18 &&
                        std::abs(hall_closure) <= hallAllowed.value() &&
                        std::abs(total) <= segBound.value() + hallAllowed.value() + std::abs(join);
        return std::make_tuple(ok, double(total_eval) / N, max_eval, maxF, total / activity, join, c);
    };
    CouplerOptions opt; opt.max_evaluations = 40; PortCoupler cp1(opt);
    auto [ok1, mean1, max1, F1, tot1, join1, c1] = coupled("H1_vacuum", NH1, false, 0, 0., circ, cp1);
    check(ok1 && max1 <= 4, "H1_coupled_iteration_hall_vacuum",
          fmt("mean evaluations %.3f, max %.0f, max|F|/|I| %.3e, total-system closure/activity %.3e, join %.3e J", mean1, max1, F1, tot1, join1));
    PortCoupler cp2 = cp1;   // carries the port-impedance estimate
    auto [ok2, mean2, max2, F2, tot2, join2, c2] = coupled("H2_gap_conduction", NH2, true, NH2 / 4, 20., c1, cp2);
    check(ok2, "H2_coupled_iteration_gap_conduction",
          fmt("mean evaluations %.3f, max %.0f, max|F|/|I| %.3e, total-system closure/activity %.3e, join %.3e J", mean2, max2, F2, tot2, join2));
    std::printf("%s %d failure(s)\n", failures ? "FAILED" : "ALL PASSED", failures);
    return failures ? 1 : 0;
}
