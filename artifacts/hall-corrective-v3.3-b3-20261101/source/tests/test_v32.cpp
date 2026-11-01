// Hall corrective v3.2: qualification of the admitted Gauss ledger (absolute Gauss admission relative to
// the ledger) and of the resident-owner checkpoint.
//   test_v32 a09              TRACE OUT NVRTC UUID
//       A09 counterexamples (inconsistent charge -e x 1e9 at the audit point, 0.01 T single-face divB) on the
//       host Model/Owner and the resident owner (bootstrap, device admission, resealed host checkpoints),
//       nonzero-ledger bootstrap, host checkpoint ledger tamper.
//   test_v32 ledger-run       TRACE OUT STEPS EVERY CONTINUE NVRTC UUID
//       field-only run at 10 fs from the supported start past 1.24 ns with the ledger ratios every EVERY steps;
//       at the end: injected charge at every node (host) and at sampled nodes (device), the admissible window,
//       the A09 counterexamples on the late state, ledger tamper (device slot, checkpoint with and without
//       reseal), checkpoint save/restore with CONTINUE steps of bit-identical continuation (in process, and
//       restore into a Model with another dt upper bound).
//   test_v32 restore-continue TRACE CHECKPOINT REFERENCE_JSON OUT NVRTC UUID
//       fresh-process restore of the ledger-run checkpoint and continuation compared with the reference chain.
//   test_v32 dense            TRACE OUT WEIGHT EV DT STEPS cpu|gpu NVRTC UUID [v6|v7 [sinks]]   (v7: native v7 motion forms;
//                             sinks: 20 extra electron/ion pairs at the lower foil face, electrons collected in the first steps)
//   test_v32 prerun           TRACE OUT DT T_END_S SAMPLE_S CHECKPOINT_S DT_MAX NVRTC UUID
//       field-only (vacuum) run of the resident owner from the supported start (DC seed at the first C19 knot) to
//       T_END_S on an exact landing grid (multiples of SAMPLE_S from 0; each interval split into equal steps <= DT);
//       a time-series row at every landing time and a checkpoint at every multiple of CHECKPOINT_S and at T_END_S
//       (CHECKPOINT_S = 0: only at T_END_S). Model: resident backend, dt upper bound DT_MAX, native v7 motion forms.
//   test_v32 compare          TRACE CHECKPOINT_A CHECKPOINT_B OUT
//       energy-norm differences of E and B and the terminal voltage of two resident checkpoints at the same time.
//       coupled closed steps of the coordinator v7 dense blob (400 co-located electron/H2O+ pairs of weight WEIGHT in
//       native x 2.6-3.4 mm, y 5.1-5.9 mm, z -3..3 mm, Maxwellian electrons at EV, ions at 0.025 eV, generator seed
//       20261101; without the qual-fixture sources) on the resident owner at fixed DT (a step refused by the midpoint
//       fixed point is retried as two half steps); per step the work identity, its v3.1 and v3.2 tolerances and the
//       absolute Gauss values (v3.1 form and v3.2 guard); at the step of smallest net exchange the work guard is
//       exercised end to end through MotionAssembly::assemble with inconsistent kinetic changes.
#include "resident_owner.hpp"
#include "gpu/gpu_backends.hpp"
#include "hall-finite-foil-em-20261029/axial_particles.hpp"
namespace hc = hall_corrective_v3;
#include "motion_fixtures.hpp"
#include <algorithm>
#include <numeric>
#include <malloc.h>
#include <fstream>
#include <iostream>

namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
namespace gp = hall_corrective_v3::gpu;
using J = hc::Json;
using steady = std::chrono::steady_clock;

namespace {
constexpr long double eps_ld = std::numeric_limits<double>::epsilon();
double ms_since(steady::time_point t) { return std::chrono::duration<double, std::milli>(steady::now() - t).count(); }
struct Results {
    J cases = J::array();
    int failed = 0;
    void record(const std::string& group, const std::string& name, bool ok, J evidence = J::object()) {
        cases.push_back({{"group", group}, {"name", name}, {"passed", ok}, {"evidence", std::move(evidence)}});
        if (!ok) ++failed;
        std::cout << group << " | " << name << ": " << (ok ? "pass" : "FAIL") << std::endl;
        if (!ok) std::cout << "  evidence: " << cases.back()["evidence"].dump().substr(0, 3000) << std::endl;
    }
};
struct Boundary { std::string field_sha256, heat_sha256, source_sha256, mesh_sha256, material_sha256, common_clock_sha256; double time_s = 0; std::uint64_t accepted_step = 0; };
template <class F> J attempt(F&& f) {
    try { f(); return J{{"accepted", true}}; }
    catch (const std::exception& e) { return J{{"accepted", false}, {"refusal", e.what()}}; }
}
bool refused_with(const J& a, const std::string& text) { return !a.at("accepted").get<bool>() && a.at("refusal").get<std::string>().find(text) != std::string::npos; }
hc::ChargeLedger declared(const hc::Model& m, const axial::State& s, hc::FreeCharge particles) {
    const auto& g = m.mesh(); hc::ChargeLedger q; q.particles = std::move(particles);
    q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
    q.conductor_declaration = hc::sha256("parity: declared DC conductor/PEC induced charge");
    q.exterior_declaration = hc::sha256("parity: declared DC exterior termination charge");
    const auto field = hc::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = field[k]; else if (g.outer_node[k]) q.exterior_C[k] = field[k]; }
    return q;
}
struct Setup {
    std::shared_ptr<gp::Device> dev;
    std::shared_ptr<const hc::Model> cpu, resident, resident_dt13;
    std::string descriptor;
};
Setup setup(const std::string& trace, const std::string& nvrtc, const std::string& uuid) {
    Setup s;
    s.dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    s.cpu = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    s.descriptor = hc::resident_descriptor_sha256(s.dev, *s.cpu);
    s.resident = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<hc::ResidentBackend>(s.descriptor));
    s.resident_dt13 = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<hc::ResidentBackend>(s.descriptor), hc::CarrierLimits{},
                                                        hc::ThermalDeclaration{}, 1e-13);
    return s;
}
// Values-only digests: physical state (field words, clock, free/conductor/exterior charge, heat; the v3.1
// definition, so v3.1 and v3.2 values compare) and the admitted Gauss ledger (empty = zero).
std::string physical_sha256(const axial::State& s, const hc::ChargeLedger& q, const hc::ThermalV3& th) {
    hc::Digest d("physical-state-v1");
    d.f64s(s.E); d.f64s(s.B); d.f64s(s.E_low); d.f64s(s.B_low); d.f64(s.time_s); d.u64(s.revision);
    hc::encode(d, q.particles.net_C); d.f64s(q.particles.gross_C); hc::encode(d, q.conductor_C); hc::encode(d, q.exterior_C);
    for (const auto& h : th.states()) {
        d.f64s(h.T); d.u64(h.inventory.size()); for (const auto& c : h.inventory) { d.f64(c.hi); d.f64(c.lo); }
        d.u64(h.watermarks.size()); for (const auto& [k, v] : h.watermarks) { d.text(k); d.u64(v); }
        for (double x : h.cumulative) d.f64(x);
        d.f64(h.time); d.u64(h.revision);
    }
    return d.hex();
}
std::string ledger_sha256(const hc::ChargeLedger& q, int nn) {
    hc::Digest d("admitted-gauss-ledger-values-v1");
    if (q.gauss_admitted_C.empty()) { d.u64(nn); for (int i = 0; i < nn; ++i) { d.f64(0.); d.f64(0.); } }
    else hc::encode(d, q.gauss_admitted_C);
    return d.hex();
}
// Independent per-node recomputation of the v3.2 absolute Gauss predicate (long double, as admit_constraints):
// x_n = D_n - L_n, bound_n, |D_n|, L_n.
struct Scan {
    std::vector<long double> x, bound, D, L;
    long double max_ratio = 0, max_v31_ratio = 0, max_implied = 0, max_L = 0; int node_ratio = -1, node_L = -1;
};
Scan gauss_scan(const hc::Model& m, const axial::State& s, const hc::ChargeLedger& q) {
    const auto& g = m.mesh(); const auto E = hc::words(s.E, s.E_low);
    std::vector<long double> act(g.nn);
    for (int e = 0; e < g.ne; ++e) { const auto& ed = g.edges[e]; const long double a = hc::magnitude(E[e]) * g.M[e] / ed.length; act[ed.a] += a; act[ed.b] += a; }
    const auto Q = hc::field_charge(g, s);
    Scan r; r.x.resize(g.nn); r.bound.resize(g.nn); r.D.resize(g.nn); r.L.resize(g.nn);
    for (int k = 0; k < g.nn; ++k) {
        const hc::Words L = q.gauss_admitted_C.empty() ? hc::Words() : q.gauss_admitted_C[k];
        const auto D = Q[k] - q.particles.net_C[k] - q.conductor_C[k] - q.exterior_C[k];
        const auto x = D - L;
        const long double b = (long double)m.initial_gauss_solve_bound_C() +
                              256 * eps_ld * (act[k] + q.particles.gross_C[k] + hc::magnitude(q.conductor_C[k]) + hc::magnitude(q.exterior_C[k]));
        r.x[k] = x.wide(); r.bound[k] = b; r.D[k] = D.wide(); r.L[k] = L.wide();
        const long double ratio = std::abs(r.x[k]) / b;
        if (ratio > r.max_ratio) { r.max_ratio = ratio; r.node_ratio = k; }
        r.max_v31_ratio = std::max(r.max_v31_ratio, std::abs(r.D[k]) / b);
        r.max_implied = std::max(r.max_implied, std::abs(r.D[k]) / (b + std::abs(r.L[k])));
        if (std::abs(r.L[k]) > r.max_L) { r.max_L = std::abs(r.L[k]); r.node_L = k; }
    }
    return r;
}
J node_json(const Scan& sc, int k) {
    if (k < 0) return nullptr;
    return {{"node", k}, {"D_C", double(sc.D[k])}, {"L_C", double(sc.L[k])}, {"D_minus_L_C", double(sc.x[k])}, {"bound_C", double(sc.bound[k])},
            {"guard_ratio", double(std::abs(sc.x[k]) / sc.bound[k])}, {"v31_form_ratio", double(std::abs(sc.D[k]) / sc.bound[k])},
            {"D_over_L", sc.L[k] != 0 ? double(sc.D[k] / sc.L[k]) : 0.}};
}
J scan_json(const Scan& sc) {
    return {{"max_guard_ratio_abs_D_minus_L_over_bound", double(sc.max_ratio)}, {"node_of_max_guard_ratio", node_json(sc, sc.node_ratio)},
            {"max_v31_form_ratio_abs_D_over_bound", double(sc.max_v31_ratio)}, {"max_implied_ratio_abs_D_over_bound_plus_abs_L", double(sc.max_implied)},
            {"max_abs_L_C", double(sc.max_L)}, {"node_of_max_abs_L", node_json(sc, sc.node_L)}};
}
// divB per cell as admit_constraints (long double); returns the cell with the largest |divB|/bound.
J divb_scan(const axial::Mesh& g, const axial::State& s) {
    const auto B = hc::words(s.B, s.B_low);
    long double worst = 0, wd = 0, wb = 0; int cell = -1; long violations = 0;
    for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
        hc::Words d; long double scale = 0;
        for (auto [a, b, h] : std::array<std::tuple<int, int, double>, 3>{{{g.bx(i + 1, j, k), g.bx(i, j, k), g.dx}, {g.by(i, j + 1, k), g.by(i, j, k), g.dy}, {g.bz(i, j, k + 1), g.bz(i, j, k), g.dz}}}) {
            d = d + (B[a] - B[b]) / h; scale += (hc::magnitude(B[a]) + hc::magnitude(B[b])) / h; }
        const long double bound = 128 * eps_ld * scale, v = std::abs(d.wide());
        if (!(v <= bound)) ++violations;
        const long double r = bound > 0 ? v / bound : (v > 0 ? INFINITY : 0);
        if (r > worst) { worst = r; wd = v; wb = bound; cell = (k * g.ny + j) * g.nx + i; }
    }
    return {{"max_divB_over_bound", double(worst)}, {"divB_T_per_m", double(wd)}, {"bound_T_per_m", double(wb)}, {"cell", cell}, {"cells_violating", violations}};
}
std::vector<n::ArenaRateEntry> a09_census(const axial::Mesh& g) {
    return {mf::row(1, "A09-manufactured-electron", "electron", {.00021 + g.side / 2, -.00017 + g.side / 2, .001}, {0, 0, 0}, -n::pic_e, mf::me, 1e9)};
}
// Free charge plus an injected charge (net and gross) at node k.
hc::ChargeLedger inject(hc::ChargeLedger q, int k, double charge) {
    q.particles.net_C[k] = q.particles.net_C[k] + hc::Words(charge);
    q.particles.gross_C[k] += std::abs(charge);
    return q;
}
hc::ChargeLedger add_free(hc::ChargeLedger q, const hc::FreeCharge& add) {
    for (std::size_t k = 0; k < q.particles.net_C.size(); ++k) { q.particles.net_C[k] = q.particles.net_C[k] + add.net_C[k]; q.particles.gross_C[k] += add.gross_C[k]; }
    return q;
}
// Device admission of a modified state in the candidate slot of a scratch resident owner.
J device_admit(hc::ResidentOwner& scratch, const axial::State& s, const hc::ChargeLedger& q) {
    return attempt([&] { scratch.device_field().upload(1, s, q); auto m = scratch.accepted_meta(); m.residue_present = !q.residue_C.empty();
                         (void)scratch.device_field().admit(1, m); });
}
J read_json(const std::filesystem::path& p) { std::ifstream in(p); J j; in >> j; return j; }
void write_json(const std::filesystem::path& p, const J& j) { std::ofstream(p, std::ios::binary) << j.dump() << '\n'; }
// Resident digest of a modified state under the configuration (scratch owner candidate slot).
std::string resident_digest(hc::ResidentOwner& scratch, const axial::State& s, const hc::ChargeLedger& q) {
    scratch.device_field().upload(1, s, q); auto m = scratch.accepted_meta(); m.residue_present = !q.residue_C.empty();
    m.time_s = s.time_s; m.revision = s.revision;
    return scratch.device_field().state_digest(1, hc::api, scratch.configuration_sha256(), m);
}
struct Chain { std::vector<std::string> field, heat; std::string physical, ledger; };
template <class OwnerT> Chain continue_steps(OwnerT& res, long steps, double dt) {
    Chain c;
    for (long k = 0; k < steps; ++k) {
        auto p = res.prepare_vacuum(dt);
        Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
        hc::need(res.preflight_vacuum_final(p, b), "continuation preflight");
        c.field.push_back(p.candidate_field_sha256()); c.heat.push_back(p.candidate_heat_sha256());
        res.publish_vacuum_prevalidated(std::move(p));
    }
    axial::State s; hc::ChargeLedger q; res.accepted_state(s, q);
    c.physical = physical_sha256(s, q, res.thermal()); c.ledger = ledger_sha256(q, int(q.particles.net_C.size()));
    return c;
}
J chain_json(const Chain& c) { return {{"field_sha256", c.field}, {"heat_sha256", c.heat}, {"physical_sha256", c.physical}, {"ledger_sha256", c.ledger}}; }

// ---------------------------------------------------------------- A09 (bootstrap / host checkpoints)
int a09(const std::string& trace, const std::filesystem::path& out, const std::string& nvrtc, const std::string& uuid) {
    Results R; auto S = setup(trace, nvrtc, uuid);
    std::filesystem::create_directories(out);
    const auto& g = S.cpu->mesh(); const auto o = mf::geometry(g); const auto seed = S.cpu->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    const auto control = declared(*S.cpu, seed, hc::charge_of_census(g, o, none));
    const auto census = a09_census(g);
    const auto cc = hc::charge_of_census(g, o, census);
    const auto honest = declared(*S.cpu, seed, cc);
    hc::ResidentOwner scratch(S.resident, S.dev, seed, control, o, none);
    // ---- A09 Gauss
    {
        const auto sc = gauss_scan(*S.cpu, seed, honest);
        const auto sc0 = gauss_scan(*S.cpu, seed, control);
        J ev{{"scan_with_A09_charge", scan_json(sc)}, {"scan_control", scan_json(sc0)}, {"initial_gauss_solve_bound_C", S.cpu->initial_gauss_solve_bound_C()},
             {"census_charge_total_C", double([&] { long double t = 0; for (const auto& w : cc.net_C) t += w.wide(); return t; }())}};
        ev["Model_admit"] = attempt([&] { S.cpu->admit(seed, honest); });
        ev["Owner_bootstrap"] = attempt([&] { hc::Owner ow(S.cpu, seed, honest, o, census); });
        ev["ResidentOwner_bootstrap"] = attempt([&] { hc::ResidentOwner ow(S.resident, S.dev, seed, honest, o, census); });
        ev["device_admission"] = device_admit(scratch, seed, honest);
        ev["control_Model_admit"] = attempt([&] { S.cpu->admit(seed, control); });
        ev["control_device_admission"] = device_admit(scratch, seed, control);
        const std::string t = "absolute Gauss constraint disagrees with actual declared charge";
        const bool ok = refused_with(ev["Model_admit"], t) && refused_with(ev["Owner_bootstrap"], t) && refused_with(ev["ResidentOwner_bootstrap"], t) &&
                        refused_with(ev["device_admission"], t) && ev["control_Model_admit"]["accepted"].get<bool>() &&
                        ev["control_device_admission"]["accepted"].get<bool>() && sc.max_ratio > 1e9;
        ev["margin_worst_node_abs_D_minus_L_over_bound"] = double(sc.max_ratio);
        R.record("A09-gauss", "4.81e-11 C inconsistent charge refused at bootstrap (host Model/Owner, resident bootstrap, device admission); "
                 "margin |D-L|/bound at the worst node", ok, ev);
    }
    // ---- a nonzero admitted Gauss ledger cannot be supplied at bootstrap (the A09 mismatch "explained" by a forged ledger)
    {
        const auto sc = gauss_scan(*S.cpu, seed, honest);
        auto forged = honest; forged.gauss_admitted_C.assign(g.nn, hc::Words());
        for (int k = 0; k < g.nn; ++k) if (sc.D[k] != 0) forged.gauss_admitted_C[k] = hc::Words(double(sc.D[k]));
        J ev;
        ev["Owner_bootstrap"] = attempt([&] { hc::Owner ow(S.cpu, seed, forged, o, census); });
        ev["ResidentOwner_bootstrap"] = attempt([&] { hc::ResidentOwner ow(S.resident, S.dev, seed, forged, o, census); });
        const std::string t = "bootstrap admitted Gauss ledger must be zero";
        R.record("A09-gauss", "a bootstrap ledger that would cover the A09 mismatch is refused (host and resident)",
                 refused_with(ev["Owner_bootstrap"], t) && refused_with(ev["ResidentOwner_bootstrap"], t), ev);
    }
    // ---- A09 divB
    {
        const int face = g.bx(22, 22, g.kshort + 1);
        auto bad = seed; bad.B[face] += .01;
        J ev{{"face", face}, {"divB_scan_perturbed", divb_scan(g, bad)}, {"divB_scan_control", divb_scan(g, seed)}};
        ev["Model_admit"] = attempt([&] { S.cpu->admit(bad, control); });
        ev["Owner_bootstrap"] = attempt([&] { hc::Owner ow(S.cpu, bad, control, o, none); });
        ev["ResidentOwner_bootstrap"] = attempt([&] { hc::ResidentOwner ow(S.resident, S.dev, bad, control, o, none); });
        ev["device_admission"] = device_admit(scratch, bad, control);
        const std::string t = "absolute divB constraint failed";
        R.record("A09-divB", "50 T/m single-face divB refused (host Model/Owner, resident bootstrap, device admission); margin divB/bound",
                 refused_with(ev["Model_admit"], t) && refused_with(ev["Owner_bootstrap"], t) && refused_with(ev["ResidentOwner_bootstrap"], t) &&
                     refused_with(ev["device_admission"], t) && ev["divB_scan_perturbed"]["max_divB_over_bound"].get<double>() > 1e9,
                 ev);
    }
    // ---- host Owner checkpoints after 3 steps (nonzero ledger): restore, resealed divB and A09 charge, ledger tamper
    {
        hc::Owner live(S.cpu, seed, control, o, none);
        for (int k = 0; k < 3; ++k) {
            auto p = live.prepare_vacuum(1e-14);
            Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
            hc::need(live.preflight_vacuum_final(p, b), "a09 host preflight");
            live.publish_vacuum_prevalidated(std::move(p));
        }
        const auto dir = out / "host-checkpoint"; std::filesystem::remove_all(dir); live.save(dir);
        auto variant = [&](const std::string& name, auto&& edit, bool reseal) {
            const auto d2 = out / ("host-checkpoint-" + name); std::filesystem::remove_all(d2);
            std::filesystem::copy(dir, d2, std::filesystem::copy_options::recursive);
            J j = read_json(dir / "checkpoint.json"); edit(j);
            if (reseal) {
                auto st = hc::checkpoint_detail::field(j.at("field"), g); auto q = hc::checkpoint_detail::charge(j.at("charge"), g);
                for (const auto& x : j.at("bank_residue_C")) q.residue_C.push_back(x.get<double>());
                q.gauss_admitted_C = hc::checkpoint_detail::ledger_words(j.at("admitted_gauss_C"), g.nn);
                j["field_sha256"] = S.cpu->state_sha256(st, q);
            }
            write_json(d2 / "checkpoint.json", j);
            return d2;
        };
        J ev;
        ev["ledger_nonzero_after_3_steps"] = !hc::zero_words(live.charge().gauss_admitted_C);
        ev["restore_unmodified"] = attempt([&] { auto b = hc::Owner::restore(S.cpu, o, none, dir);
                                                 hc::need(b->field_sha256() == live.field_sha256(), "restored digest differs"); });
        const int face = g.bx(22, 22, g.kshort + 1);
        ev["restore_resealed_divB"] = attempt([&] { (void)hc::Owner::restore(S.cpu, o, none, variant("divB", [&](J& j) {
                                                       j["field"]["B"][face] = j["field"]["B"][face].get<double>() + .01; }, true)); });
        ev["restore_resealed_A09_charge"] = attempt([&] { (void)hc::Owner::restore(S.cpu, o, census, variant("a09", [&](J& j) {
                                                             j["charge"]["free_C"] = hc::record(cc.net_C); j["charge"]["gross_free_C"] = cc.gross_C; }, true)); });
        int node = 0; { const auto& L = live.charge().gauss_admitted_C; for (int k = 0; k < g.nn; ++k) if (std::abs(L[k].hi) > std::abs(L[node].hi)) node = k; }
        ev["ledger_tamper_node"] = node;
        ev["restore_ledger_one_ulp_not_resealed"] = attempt([&] { (void)hc::Owner::restore(S.cpu, o, none, variant("ledger-ulp", [&](J& j) {
            auto& w = j["admitted_gauss_C"][node]; w[0] = std::nextafter(w[0].get<double>(), INFINITY); }, false)); });
        ev["restore_ledger_plus_1e-20_resealed"] = attempt([&] { (void)hc::Owner::restore(S.cpu, o, none, variant("ledger-1e-20", [&](J& j) {
            auto& w = j["admitted_gauss_C"][node]; w[0] = w[0].get<double>() + 1e-20; w[1] = 0.; }, true)); });
        R.record("checkpoint-host", "host Owner checkpoint carries the ledger (nonzero after 3 steps) and restores to the same digest",
                 ev["ledger_nonzero_after_3_steps"].get<bool>() && ev["restore_unmodified"]["accepted"].get<bool>(), ev);
        R.record("A09-divB", "resealed host checkpoint with the 0.01 T face perturbation refused", refused_with(ev["restore_resealed_divB"], "absolute divB constraint failed"));
        R.record("A09-gauss", "resealed host checkpoint with the A09 free charge (and its census) refused by the absolute Gauss guard",
                 refused_with(ev["restore_resealed_A09_charge"], "absolute Gauss constraint disagrees"));
        R.record("ledger-tamper", "host checkpoint: one-ulp ledger change without reseal refused (seal)",
                 refused_with(ev["restore_ledger_one_ulp_not_resealed"], "checkpoint field/charge seal mismatch"));
        R.record("ledger-tamper", "host checkpoint: resealed ledger change of 1e-20 C refused by the absolute Gauss guard",
                 refused_with(ev["restore_ledger_plus_1e-20_resealed"], "absolute Gauss constraint disagrees"));
    }
    J summary{{"schema", "hall-v3.2-a09-v1"}, {"failed", R.failed}, {"cases", R.cases}, {"configuration_sha256_cpu", S.cpu->configuration_sha256()},
              {"configuration_sha256_resident", S.resident->configuration_sha256()}, {"device", J::parse(S.dev->device_json())}};
    std::ofstream(out / "a09.json") << summary.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// ---------------------------------------------------------------- vacuum run past 1.24 ns and the late-state tests
int ledger_run(const std::string& trace, const std::filesystem::path& out, long steps, long every, long cont, const std::string& nvrtc, const std::string& uuid) {
    Results R; auto S = setup(trace, nvrtc, uuid);
    std::filesystem::create_directories(out);
    const auto& g = S.resident->mesh(); const auto o = mf::geometry(g); const auto seed = S.resident->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    const double dt = 1e-14;
    const auto q0 = declared(*S.resident, seed, hc::charge_of_census(g, o, none));
    const auto scan0 = gauss_scan(*S.resident, seed, q0);   // initial (bootstrap) defects D0
    hc::ResidentOwner res(S.resident, S.dev, seed, q0, o, none);
    std::ofstream trace_out(out / "ledger-trace.jsonl");
    hc::Vec cen, bnd, raw, led;
    double max_step_guard_ratio = 0, max_step_implied = 0; long first_v31_over_1 = -1; double first_v31_time = 0;
    long refused_at = -1; std::string refusal;
    const int watch = 3747;
    const auto t0 = steady::now();
    for (long step = 1; step <= steps; ++step) {
        try {
            auto p = res.prepare_vacuum(dt);
            Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
            hc::need(res.preflight_vacuum_final(p, b), "ledger-run preflight");
            const auto& c = p.candidate_constraints();
            max_step_guard_ratio = std::max(max_step_guard_ratio, c.max_gauss_node_ratio);
            max_step_implied = std::max(max_step_implied, c.max_gauss_defect_over_bound_plus_ledger);
            const double V = p.field_receipt().terminal_V, I = p.field_receipt().current_A;
            res.publish_vacuum_prevalidated(std::move(p));
            if (step == 1 || step % every == 0 || step == steps) {
                const auto cs = res.device_field().gauss_profile(0, res.accepted_meta(), cen, bnd, &raw, &led);
                double v31 = 0; int v31_node = -1; double top_led = 0; int top = -1;
                for (int k = 0; k < g.nn; ++k) {
                    const double r = raw[k] / bnd[k]; if (r > v31) { v31 = r; v31_node = k; }
                    if (std::abs(led[k]) > top_led) { top_led = std::abs(led[k]); top = k; }
                }
                if (v31 > 1 && first_v31_over_1 < 0) { first_v31_over_1 = step; first_v31_time = res.time_s(); }
                // ledger explanation at the 16 nodes of largest |L|
                std::vector<int> idx(g.nn); std::iota(idx.begin(), idx.end(), 0);
                std::partial_sort(idx.begin(), idx.begin() + 16, idx.end(), [&](int a, int b2) { return std::abs(led[a]) > std::abs(led[b2]); });
                double rmin = INFINITY, rmax = 0;
                for (int i = 0; i < 16; ++i) { const int k = idx[i]; if (led[k] != 0) { const double r = raw[k] / std::abs(led[k]); rmin = std::min(rmin, r); rmax = std::max(rmax, r); } }
                trace_out << J{{"step", step}, {"time_s", res.time_s()}, {"terminal_V", V}, {"current_A", I},
                               {"guard_ratio_max_abs_D_minus_L_over_bound", cs.max_gauss_node_ratio}, {"max_abs_D_minus_L_C", cs.max_gauss_C},
                               {"v31_form_ratio_max_abs_D_over_bound", v31}, {"v31_form_worst_node", v31_node},
                               {"implied_ratio_max_abs_D_over_bound_plus_abs_L", cs.max_gauss_defect_over_bound_plus_ledger},
                               {"max_abs_L_C", cs.max_gauss_ledger_C}, {"max_abs_D_C", cs.max_gauss_defect_C}, {"node_of_max_abs_L", top},
                               {"abs_D_over_abs_L_top16_min", rmin}, {"abs_D_over_abs_L_top16_max", rmax},
                               {"node3747", {{"abs_D_C", raw[watch]}, {"L_C", led[watch]}, {"abs_D_minus_L_C", cen[watch]}, {"bound_C", bnd[watch]}}},
                               {"max_T_K", res.thermal().max_temperature_K()}}.dump() << '\n';
                trace_out.flush();
            }
        } catch (const std::exception& e) { refused_at = step; refusal = e.what(); break; }
    }
    const double wall_s = ms_since(t0) / 1e3;
    J run{{"steps_requested", steps}, {"refused_at_step", refused_at}, {"refusal", refusal}, {"end_time_s", res.time_s()}, {"revision", res.revision()},
          {"wall_s", wall_s}, {"ms_per_step", 1e3 * wall_s / std::max<long>(1, res.revision())}, {"max_step_guard_ratio", max_step_guard_ratio},
          {"max_step_implied_ratio", max_step_implied}, {"first_sample_with_v31_form_ratio_over_1", {{"step", first_v31_over_1}, {"time_s", first_v31_time}}},
          {"v31_refusal_of_record", {{"step", 48650}, {"time_s", 1.2363e-9}, {"node", 3747}}}};
    R.record("ledger-run", "field-only run at 10 fs from the supported start passes 1.24 ns (v3.1 refused at step 48,650, 1.2363 ns) with every guard",
             refused_at < 0 && res.time_s() > 1.24e-9 && max_step_guard_ratio <= 1, run);
    axial::State s; hc::ChargeLedger q; res.accepted_state(s, q);
    const auto sc = gauss_scan(*S.resident, s, q);
    {
        J ev = scan_json(sc); ev["node3747"] = node_json(sc, watch);
        // ledger explanation: the drift D_n - D0_n against L_n at the nodes with |L_n| >= 1e-25 C, and the largest
        // |D_n - L_n - D0_n| over all nodes (rounding of the word operations)
        long considered = 0; long double rmin = INFINITY, rmax = 0, resid = 0;
        for (int k = 0; k < g.nn; ++k) {
            resid = std::max(resid, std::abs(sc.x[k] - scan0.x[k]));
            if (std::abs(sc.L[k]) >= 1e-25L) { ++considered; const long double r = (sc.D[k] - scan0.D[k]) / sc.L[k]; rmin = std::min(rmin, r); rmax = std::max(rmax, r); }
        }
        ev["drift_over_L_min"] = double(rmin); ev["drift_over_L_max"] = double(rmax); ev["nodes_with_abs_L_ge_1e-25"] = considered;
        ev["max_abs_D_minus_L_minus_D0_C"] = double(resid); ev["initial_scan"] = scan_json(scan0);
        ev["host_admission"] = attempt([&] { S.resident->admit(s, q); });
        const bool past = s.time_s > 1.24e-9;
        R.record("ledger-run", "end state: host recomputation of the guard (|D-L|/bound <= 1); drift D-D0 equals L at every node with |L| >= 1e-25 C; "
                 "the v3.1 form |D|/bound exceeds 1 (past 1.24 ns)",
                 ev["host_admission"]["accepted"].get<bool>() && sc.max_ratio <= 1 && (!past || sc.max_v31_ratio > 1) &&
                     (considered == 0 || (rmin > 1 - 1e-9 && rmax < 1 + 1e-9)), ev);
    }
    // ---- checkpoint of the end state, then continuation of the original owner
    const auto ck = out / "checkpoint-end"; std::filesystem::remove_all(ck);
    const std::string saved_field = res.field_sha256(), saved_heat = res.heat_sha256();
    res.save(ck);
    const std::string saved_physical = physical_sha256(s, q, res.thermal());
    const auto ref = continue_steps(res, cont, dt);
    {
        J refj = chain_json(ref); refj["saved_field_sha256"] = saved_field; refj["saved_heat_sha256"] = saved_heat; refj["dt_s"] = dt; refj["steps"] = cont;
        write_json(out / "continuation-reference.json", refj);
    }
    {
        auto r2 = hc::ResidentOwner::restore(S.resident, S.dev, o, none, ck);
        axial::State s2; hc::ChargeLedger q2; r2->accepted_state(s2, q2);
        const bool same_start = r2->field_sha256() == saved_field && r2->heat_sha256() == saved_heat && physical_sha256(s2, q2, r2->thermal()) == saved_physical &&
                                ledger_sha256(q2, g.nn) == ledger_sha256(q, g.nn);
        const auto c2 = continue_steps(*r2, cont, dt);
        const bool same = c2.field == ref.field && c2.heat == ref.heat && c2.physical == ref.physical && c2.ledger == ref.ledger;
        R.record("checkpoint-resident", "resident save/restore (same process, new device state): restored digests equal the saved ones and " + std::to_string(cont) +
                 " continued steps are bit-identical (every field and heat digest, final physical state and ledger)", same_start && same,
                 {{"same_start", same_start}, {"same_chain", same}, {"final_physical_sha256", c2.physical}, {"final_ledger_sha256", c2.ledger}});
        // transfer into a Model whose dt upper bound is 1e-13 (other configuration identity), same steps at 10 fs
        auto r3 = hc::ResidentOwner::restore(S.resident_dt13, S.dev, o, none, ck);
        const auto c3 = continue_steps(*r3, cont, dt);
        R.record("checkpoint-resident", "restore into a Model with dt upper bound 1e-13 (configuration recomputed from the recorded operands): "
                 "same physical state and ledger after the same steps, different identities",
                 c3.physical == ref.physical && c3.ledger == ref.ledger && c3.field.back() != ref.field.back() &&
                     S.resident_dt13->configuration_sha256() != S.resident->configuration_sha256(),
                 {{"configuration_default", S.resident->configuration_sha256()}, {"configuration_dt13", S.resident_dt13->configuration_sha256()}});
        // a different configuration that is not a dt / backend variant is refused
        auto other = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<hc::ResidentBackend>(S.descriptor), hc::CarrierLimits{50000, 1ull << 30});
        R.record("checkpoint-resident", "restore into a Model with other carrier limits refused",
                 refused_with(attempt([&] { (void)hc::ResidentOwner::restore(other, S.dev, o, none, ck); }), "checkpoint field/forcing/mesh configuration changed"));
        // ---- injected charge at every node (host arithmetic, end state of the run)
        {
            const double e = n::pic_e;
            long refused = 0, wrong_node = 0, tested = 0; J first_bad = nullptr;
            const auto t1 = steady::now();
            for (int k = 0; k < g.nn; ++k) for (double sgn : {1., -1.}) {
                ++tested;
                const auto a = attempt([&] { S.resident->admit(s, inject(q, k, sgn * e)); });
                if (refused_with(a, "absolute Gauss constraint disagrees")) {
                    ++refused;
                    if (a["refusal"].get<std::string>().find("node=" + std::to_string(k) + " ") == std::string::npos) { ++wrong_node; if (first_bad.is_null()) first_bad = a; }
                } else if (first_bad.is_null()) first_bad = {{"node", k}, {"sign", sgn}, {"result", a}};
            }
            long double qmax = 0; int kmax = -1;
            for (int k = 0; k < g.nn; ++k) { const long double w = sc.bound[k] + std::abs(sc.x[k]); if (w > qmax) { qmax = w; kmax = k; } }
            R.record("injection", "an injected elementary charge (+e and -e, free charge) at every node of the late state is refused by the absolute Gauss guard at that node",
                     refused == tested && wrong_node == 0,
                     {{"tested", tested}, {"refused", refused}, {"refused_at_other_node", wrong_node}, {"first_unexpected", first_bad}, {"host_s", ms_since(t1) / 1e3},
                      {"largest_admissible_injection_anywhere_C_max_n(bound_n+|D_n-L_n|)", double(qmax)}, {"at_node", kmax},
                      {"elementary_charge_over_largest_admissible_injection", double(e / qmax)}, {"max_abs_L_C", double(sc.max_L)}});
        }
        // ---- admissible window at selected nodes (host and device): |q| > bound + |D-L| refused for both signs; q = D-L admitted
        {
            std::vector<int> nodes{watch, sc.node_L, sc.node_ratio};
            { int kb = 0; for (int k = 0; k < g.nn; ++k) if (sc.bound[k] > sc.bound[kb]) kb = k; nodes.push_back(kb); }
            std::mt19937_64 rng(17); std::uniform_int_distribution<int> pick(0, g.nn - 1);
            for (int i = 0; i < 24; ++i) nodes.push_back(pick(rng));
            bool ok = true; J rows = J::array();
            for (int k : nodes) {
                const double w = double((sc.bound[k] + std::abs(sc.x[k])) * (1 + 1e-6L));
                // free charge enters D with a minus sign: injecting +q moves D - L by -q
                const auto hp = attempt([&] { S.resident->admit(s, inject(q, k, w)); }), hm = attempt([&] { S.resident->admit(s, inject(q, k, -w)); });
                const auto dp = device_admit(*r2, s, inject(q, k, w)), dm = device_admit(*r2, s, inject(q, k, -w));
                const auto h0 = attempt([&] { S.resident->admit(s, inject(q, k, double(sc.x[k]))); });
                const auto d0 = device_admit(*r2, s, inject(q, k, double(sc.x[k])));
                const bool row_ok = !hp["accepted"].get<bool>() && !hm["accepted"].get<bool>() && !dp["accepted"].get<bool>() && !dm["accepted"].get<bool>() &&
                                    h0["accepted"].get<bool>() && d0["accepted"].get<bool>();
                ok = ok && row_ok;
                rows.push_back({{"node", k}, {"bound_C", double(sc.bound[k])}, {"D_minus_L_C", double(sc.x[k])}, {"L_C", double(sc.L[k])}, {"threshold_C", w},
                                {"refused_both_signs_host", !hp["accepted"].get<bool>() && !hm["accepted"].get<bool>()},
                                {"refused_both_signs_device", !dp["accepted"].get<bool>() && !dm["accepted"].get<bool>()},
                                {"q_equal_D_minus_L_admitted_host_device", h0["accepted"].get<bool>() && d0["accepted"].get<bool>()}, {"ok", row_ok}});
            }
            R.record("injection", "admissible window per node is [D-L-bound, D-L+bound] (seed level, independent of |L|): |q| = (bound+|D-L|)(1+1e-6) refused "
                     "for both signs on host and device at 28 nodes (3747, max |L|, max ratio, max bound, 24 random); q = D-L admitted", ok, {{"rows", rows}});
        }
        // ---- the literal form |D| <= bound + |L| would admit q = 2L at the node of largest |L|; the v3.2 guard admits it only
        //      when |D - L - 2L| <= bound (|L| small against the bound); once |L| exceeds the bound it is refused
        {
            const int k = sc.node_L; const double qv = double(2 * sc.L[k]);
            const auto qq = inject(q, k, qv);
            const auto sc2 = gauss_scan(*S.resident, s, qq);
            const bool literal_admits = std::abs(sc2.D[k]) <= sc2.bound[k] + std::abs(sc2.L[k]);
            const bool predicted_refusal = std::abs(sc2.x[k]) > sc2.bound[k];
            const auto h = attempt([&] { S.resident->admit(s, qq); });
            const auto d = device_admit(*r2, s, qq);
            const bool host_refused = refused_with(h, "absolute Gauss"), dev_refused = refused_with(d, "absolute Gauss");
            R.record("injection", "charge 2L at the node of largest |L| (admissible under |D| <= bound + |L|): v3.2 decision (host and device) "
                     "= predicted |D-L| > bound; refused once |L| > bound",
                     literal_admits && host_refused == predicted_refusal && dev_refused == predicted_refusal && (predicted_refusal || sc.max_L <= sc.bound[k]),
                     {{"node", k}, {"injected_C", qv}, {"literal_form_admits", literal_admits}, {"refused", host_refused && dev_refused},
                      {"abs_L_over_bound", double(sc.max_L / sc.bound[k])}, {"host", h}, {"device", d}, {"after", node_json(sc2, k)}});
        }
        // ---- A09 counterexamples on the late state
        {
            const auto cc = hc::charge_of_census(g, o, a09_census(g));
            const auto qq = add_free(q, cc);
            const auto sc2 = gauss_scan(*S.resident, s, qq);
            const auto h = attempt([&] { S.resident->admit(s, qq); });
            const auto d = device_admit(*r2, s, qq);
            const double literal_margin = [&] { long double m = 0; for (int k = 0; k < g.nn; ++k) m = std::max(m, std::abs(sc2.D[k]) / (sc2.bound[k] + std::abs(sc2.L[k]))); return double(m); }();
            R.record("A09-gauss", "A09 charge injected into the late state (t > 1.24 ns, ledger accumulated) refused on host and device; margins stated",
                     refused_with(h, "absolute Gauss constraint disagrees") && refused_with(d, "absolute Gauss constraint disagrees"),
                     {{"margin_abs_D_minus_L_over_bound", double(sc2.max_ratio)}, {"margin_literal_form_abs_D_over_bound_plus_abs_L", literal_margin},
                      {"worst_node", node_json(sc2, sc2.node_ratio)}, {"time_s", s.time_s}, {"host", h}, {"device", d}});
            auto bad = s; const int face = g.bx(22, 22, g.kshort + 1); bad.B[face] += .01;
            const auto hb = attempt([&] { S.resident->admit(bad, q); });
            const auto db = device_admit(*r2, bad, q);
            R.record("A09-divB", "0.01 T face perturbation of the late state refused on host and device; margin stated",
                     refused_with(hb, "absolute divB constraint failed") && refused_with(db, "absolute divB constraint failed"),
                     {{"divB_scan", divb_scan(g, bad)}, {"divB_scan_unperturbed", divb_scan(g, s)}, {"host", hb}, {"device", db}});
        }
        // ---- ledger tamper: accepted device slot (preflight), checkpoint without reseal (seal), resealed checkpoint (guard)
        {
            const int k = sc.node_L;
            J ev{{"node", k}};
            {   // one ulp of one accepted ledger word on the device
                auto p = r3->prepare_vacuum(dt);
                axial::State s3; hc::ChargeLedger q3; r3->accepted_state(s3, q3);
                q3.gauss_admitted_C[k].lo = std::nextafter(q3.gauss_admitted_C[k].lo, INFINITY);
                r3->device_field().upload(0, s3, q3);
                Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
                ev["preflight_after_device_ledger_ulp"] = r3->preflight_vacuum_final(p, b);
            }
            auto variant = [&](const std::string& name, auto&& edit, bool reseal) {
                const auto d2 = out / ("checkpoint-" + name); std::filesystem::remove_all(d2);
                std::filesystem::copy(ck, d2, std::filesystem::copy_options::recursive);
                J j = read_json(ck / "checkpoint.json"); edit(j);
                if (reseal) {
                    auto st = hc::checkpoint_detail::field(j.at("field"), g); auto qq = hc::checkpoint_detail::charge(j.at("charge"), g);
                    qq.gauss_admitted_C = hc::checkpoint_detail::ledger_words(j.at("admitted_gauss_C"), g.nn);
                    j["field_sha256"] = resident_digest(*r2, st, qq);
                }
                write_json(d2 / "checkpoint.json", j);
                return d2;
            };
            ev["restore_ledger_ulp_not_resealed"] = attempt([&] { (void)hc::ResidentOwner::restore(S.resident, S.dev, o, none, variant("ledger-ulp", [&](J& j) {
                auto& w = j["admitted_gauss_C"][k]; w[1] = std::nextafter(w[1].get<double>(), INFINITY); }, false)); });
            ev["restore_ledger_plus_1e-20_resealed"] = attempt([&] { (void)hc::ResidentOwner::restore(S.resident, S.dev, o, none, variant("ledger-1e-20", [&](J& j) {
                auto& w = j["admitted_gauss_C"][k]; const hc::Words x = hc::Words(w[0].get<double>(), w[1].get<double>()) + hc::Words(1e-20);
                w[0] = x.hi; w[1] = x.lo; }, true)); });
            ev["restore_resealed_control"] = attempt([&] { (void)hc::ResidentOwner::restore(S.resident, S.dev, o, none, variant("resealed-control", [&](J&) {}, true)); });
            R.record("ledger-tamper", "one-ulp change of one accepted ledger word on the device refused at preflight (leaves recomputed)",
                     !ev["preflight_after_device_ledger_ulp"].get<bool>(), ev);
            R.record("ledger-tamper", "resident checkpoint: one-ulp ledger change without reseal refused (seal)",
                     refused_with(ev["restore_ledger_ulp_not_resealed"], "checkpoint field/charge seal mismatch"));
            R.record("ledger-tamper", "resident checkpoint: resealed ledger change of 1e-20 C refused by the absolute Gauss guard; resealed unmodified control restores",
                     refused_with(ev["restore_ledger_plus_1e-20_resealed"], "absolute Gauss constraint disagrees") && ev["restore_resealed_control"]["accepted"].get<bool>());
        }
    }
    J summary{{"schema", "hall-v3.2-ledger-run-v1"}, {"failed", R.failed}, {"cases", R.cases}, {"run", run}, {"end_scan", scan_json(sc)},
              {"configuration_sha256", S.resident->configuration_sha256()}, {"checkpoint", ck.string()}, {"device", J::parse(S.dev->device_json())}};
    std::ofstream(out / "ledger-run.json") << summary.dump(1) << '\n';
    return R.failed ? 1 : 0;
}


// ---------------------------------------------------------------- dense blob (work identity, absolute Gauss)
std::vector<n::ArenaRateEntry> blob_census(double weight, double eV, long pairs) {
    // the coordinator v7 dense generator (tests/driver_v7.cpp dense_packets), without the qual fixture
    std::mt19937_64 rng(20261101);
    std::uniform_real_distribution<double> u(0, 1);
    std::normal_distribution<double> nd(0, 1);
    const double water_mass_kg = 2.9915e-26;
    const double se = std::sqrt(eV * n::pic_e / mf::me), si = std::sqrt(.025 * n::pic_e / water_mass_kg);
    std::vector<n::ArenaRateEntry> out;
    for (long k = 0; k < pairs; ++k) {
        const n::PicVector x{.0026 + 8e-4 * u(rng), .0051 + 8e-4 * u(rng), -.003 + .006 * u(rng)};
        const n::PicVector ve{se * nd(rng), se * nd(rng), se * nd(rng)}, vi{si * nd(rng), si * nd(rng), si * nd(rng)};
        out.push_back(mf::row(std::uint64_t(2 * k + 1), "be" + std::to_string(k), "electron", x, ve, -n::pic_e, mf::me, weight));
        out.push_back(mf::row(std::uint64_t(2 * k + 2), "bi" + std::to_string(k), "H2O+", x, vi, n::pic_e, mf::ion_mass(), weight));
    }
    return out;
}
// Optional sink carriers: 20 electron/H2O+ pairs 20-400 nm above the lower foil inner face (native y = 2.0 mm), electrons
// at 100 eV moving into the foil (first-hit sinks in the first steps), ions thermal.
void add_sink_pairs(std::vector<n::ArenaRateEntry>& out, double weight) {
    std::uint64_t id = out.empty() ? 1 : out.back().native.id + 1;
    for (int k = 0; k < 20; ++k) {
        const n::PicVector x{.0035 + 1e-4 * k, .002 + 2e-8 * (k + 1), -.002 + 2e-4 * k};
        auto u = mf::electron_u(100); const n::PicVector ue{.1 * u[0], -u[0], 0};
        out.push_back(mf::row(id++, "se" + std::to_string(k), "electron", x, ue, -n::pic_e, mf::me, weight));
        out.push_back(mf::row(id++, "si" + std::to_string(k), "H2O+", x, {0, 300, 0}, n::pic_e, mf::ion_mass(), weight));
    }
}
hall_driver::MotionBinding dense_binding(const hc::ResidentOwner& res, const hc::Model& m, long step, double dt) {
    hall_driver::MotionBinding b; b.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = res.revision(); b.staged_arena_revision = res.revision() + 1; b.begin_s = res.time_s(); b.dt_s = dt; b.end_s = b.begin_s + dt;
    b.staged_arena_sha256 = hc::sha256("dense/staged/" + std::to_string(step)); b.mesh_sha256 = res.motion_mesh_sha256();
    b.material_sha256 = m.material_sha256(); b.species_sha256 = hc::sha256("dense/species");
    b.common_clock_sha256 = hc::sha256("dense/clock/" + std::to_string(step)); b.kinetic_policy_sha256 = res.motion_context().policy_sha256;
    return b;
}
J scales_json(const hc::WorkScales& w) {
    return {{"G_K_J", w.kinetic_gross_J}, {"C_K_J", w.kinetic_content_J}, {"G_W_J", w.field_work_gross_J}, {"Q_C", w.charge_abs_C}, {"E_hat_V_m", w.field_bound_V_m},
            {"tolerance_J", w.tolerance_J}, {"v31_tolerance_J", w.v31_tolerance_J},
            {"terms_J", {{"relative_gross", 1e-10 * (w.kinetic_gross_J + w.field_work_gross_J)}, {"kinetic_rounding", 32 * 2.220446049250313e-16 * w.kinetic_content_J},
                         {"fixed_point", 8e-18 * w.field_bound_V_m * w.charge_abs_C}}}};
}
int dense(const std::string& trace, const std::filesystem::path& out, double weight, double eV, double dt, long steps, bool gpu_push,
          const std::string& nvrtc, const std::string& uuid, bool forms_v7, bool sinks) {
    Results R; auto S = setup(trace, nvrtc, uuid);
    if (forms_v7) S.resident = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<hc::ResidentBackend>(S.descriptor), hc::CarrierLimits{},
                                                              hc::ThermalDeclaration{}, 1e-14, hc::MotionForms{true, 1e-14});
    std::filesystem::create_directories(out);
    const auto& g = S.resident->mesh(); const auto o = mf::geometry(g); const auto seed = S.resident->dc_seed();
    auto census = blob_census(weight, eV, 400);
    if (sinks) add_sink_pairs(census, weight);
    hc::ResidentOwner res(S.resident, S.dev, seed, declared(*S.resident, seed, hc::charge_of_census(g, o, census)), o, census);
    std::unique_ptr<hc::MotionProvider> provider;
    if (gpu_push) provider = std::make_unique<gp::GpuMotionProvider>(S.dev, res.motion_context(), census.size());
    else provider = std::make_unique<hc::CpuMotionProvider>();
    std::ofstream led(out / "dense-steps.jsonl");
    const std::string src = hc::sha256("dense/source");
    double phase = 0, max_ratio = 0, max_joint = 0, max_combined = 0, max_gauss_guard = 0, max_gauss_v31 = 0, min_net_rel = INFINITY;
    long first_v31_work = -1, first_v31_gauss = -1, retries = 0, refused_at = -1, accepted = 0, census_size_end = 0;
    double first_v31_work_phase = 0, first_v31_gauss_time = 0; std::string refusal;
    // state at the step of the smallest |K| / G_K (for the end-to-end work-guard tests)
    std::vector<n::ArenaRateEntry> zc_census; hc::Vec zcE, zcB; double zc_dt = dt; long zc_step = -1; std::uint64_t zc_rev = 0; double zc_time = 0;
    hc::Vec cen, bnd, raw, ledv;
    const auto t0 = steady::now();
    long step = 0;
    auto one = [&](double h) {   // one accepted closed step of length h
        auto b = dense_binding(res, *S.resident, step, h);
        hc::Vec E0, B0; res.device_field().accepted_rounded(E0, B0);
        auto p = hc::prepare_closed_proposal(res, o, census, b, src, *provider, hc::StepInputs{});
        const auto& bb = p.motion().binding();
        Boundary bd{p.candidate_field_sha256(), p.candidate_heat_sha256(), src, bb.mesh_sha256, bb.material_sha256, bb.common_clock_sha256, p.end_s(), p.endpoint_revision()};
        hc::need(res.preflight_final(p, census, bd), "dense preflight");
        const auto& m = p.motion(); const auto& bu = m.budget(); const auto& w = m.work_scales(); const auto& r = p.field_heat().field_receipt();
        // A1 admission inputs: the trial field digest and the current digest bind the spans handed to native
        hc::need(hall_driver::total_field_sha256(p.trial_E(), p.trial_B()) == m.binding().trial_field_sha256 && &p.current_A_m() == &m.current_A_m(),
                 "trial E/B spans do not reproduce trial_field_sha256");
        const double defect = bu.work_defect_J, ratio = std::abs(defect) / w.tolerance_J, v31 = std::abs(defect) / w.v31_tolerance_J;
        const double joint = std::abs(bu.kinetic_change_J - r.plasma_work_J) / w.tolerance_J;
        const auto& en = p.energy();
        const double comb = std::abs(p.combined_balance_J());
        phase += bu.plasma_dt;
        max_ratio = std::max(max_ratio, ratio); max_joint = std::max(max_joint, joint);
        if (v31 > 1 && first_v31_work < 0) { first_v31_work = step; first_v31_work_phase = phase; }
        const double net_rel = std::abs(bu.kinetic_change_J) / std::max(w.kinetic_gross_J, 1e-300);
        if (net_rel < min_net_rel && step > 2) { min_net_rel = net_rel; zc_census = census; zcE = E0; zcB = B0; zc_dt = h; zc_step = step; zc_rev = res.revision(); zc_time = res.time_s(); }
        const auto c = p.field_heat().candidate_constraints();
        max_gauss_guard = std::max(max_gauss_guard, c.max_gauss_node_ratio);
        auto next = hc::next_census(census, m);
        res.publish_prevalidated(std::move(p));
        ++accepted;
        // v3.1 form of the absolute Gauss check on the accepted state (per node |D| / bound)
        (void)res.device_field().gauss_profile(0, res.accepted_meta(), cen, bnd, &raw, &ledv);
        double v31g = 0; for (int k = 0; k < g.nn; ++k) v31g = std::max(v31g, raw[k] / bnd[k]);
        max_gauss_v31 = std::max(max_gauss_v31, v31g);
        if (v31g > 1 && first_v31_gauss < 0) { first_v31_gauss = step; first_v31_gauss_time = res.time_s(); }
        max_combined = std::max(max_combined, comb / (2e-20 + 1e-10 * en.scale_J + w.tolerance_J));
        led << J{{"step", step}, {"time_s", res.time_s()}, {"dt_s", h}, {"omega_p_dt", bu.plasma_dt}, {"omega_p_t", phase}, {"carriers", census.size()},
                 {"K_J", bu.kinetic_change_J}, {"W_J", bu.field_work_J}, {"defect_J", defect}, {"plasma_work_J", r.plasma_work_J},
                 {"defect_over_v32_tolerance", ratio}, {"defect_over_v31_tolerance", v31}, {"joint_over_v32_tolerance", joint},
                 {"net_over_gross", net_rel}, {"scales", scales_json(w)}, {"combined_balance_J", p.combined_balance_J()},
                 {"gauss_v32_guard_ratio", c.max_gauss_node_ratio}, {"gauss_v31_form_ratio", v31g}, {"max_abs_L_C", c.max_gauss_ledger_C},
                 {"max_abs_D_C", c.max_gauss_defect_C}, {"iterations", p.iterations()}, {"forms_v7", m.forms_v7()},
                 {"census_sha256", m.census_sha256()}, {"record_v7_sha256", m.record_v7_sha256()}, {"record_v3_sha256", m.record_v3_sha256()},
                 {"sinks", m.sinks().size()}}.dump() << '\n';
        census = std::move(next);
    };
    for (step = 1; step <= steps; ++step) {
        try { one(dt); }
        catch (const std::exception& e) {
            const std::string why = e.what();
            if (why.find("did not reach the identical authoritative midpoint") == std::string::npos) { refused_at = step; refusal = why; break; }
            ++retries;
            try { one(.5 * dt); one(.5 * dt); }
            catch (const std::exception& e2) { refused_at = step; refusal = std::string("after retry: ") + e2.what(); break; }
        }
    }
    census_size_end = long(census.size());
    J run{{"weight", weight}, {"electron_eV", eV}, {"dt_s", dt}, {"push", provider->name()}, {"steps_requested", steps}, {"accepted_steps", accepted},
          {"retries_two_half_steps", retries}, {"refused_at_step", refused_at}, {"refusal", refusal}, {"end_time_s", res.time_s()},
          {"elapsed_from_start_s", res.time_s() - seed.time_s}, {"omega_p_t_rad", phase}, {"carriers_end", census_size_end},
          {"max_defect_over_v32_tolerance", max_ratio}, {"max_joint_over_v32_tolerance", max_joint}, {"max_combined_over_bound", max_combined},
          {"first_step_v31_work_tolerance_exceeded", first_v31_work}, {"omega_p_t_at_first_v31_work_excess", first_v31_work_phase},
          {"min_abs_K_over_G_K", min_net_rel}, {"max_gauss_v32_guard_ratio", max_gauss_guard}, {"max_gauss_v31_form_ratio", max_gauss_v31},
          {"first_step_v31_gauss_form_exceeded", first_v31_gauss}, {"time_at_first_v31_gauss_excess_s", first_v31_gauss_time},
          {"wall_s", ms_since(t0) / 1e3}};
    R.record("dense", "dense blob runs the requested steps without refusal (v3.2 work tolerance and admitted Gauss ledger)", refused_at < 0, run);
    // ---- end-to-end work guard through MotionAssembly::assemble at the step of the smallest net exchange
    if (zc_step > 0) {
        auto b = dense_binding(res, *S.resident, zc_step, zc_dt); b.accepted_step = zc_rev; b.staged_arena_revision = zc_rev + 1; b.begin_s = zc_time; b.end_s = zc_time + zc_dt;
        b.trial_field_sha256 = hall_driver::total_field_sha256(zcE, zcB);
        const auto& ctx = res.motion_context();
        const auto r0 = hc::MotionAssembly::cpu(ctx, zc_census, b, {}, zcE, zcB);
        const auto w = r0.work_scales();
        auto arrays = [&] {
            hc::MotionArrays m;
            m.endpoints.assign(r0.endpoints().begin(), r0.endpoints().end()); m.sinks.assign(r0.sinks().begin(), r0.sinks().end());
            m.births.assign(r0.births().begin(), r0.births().end()); m.paths.assign(r0.paths().begin(), r0.paths().end());
            m.current_A_m = r0.current_A_m(); m.initial_charge_C = r0.initial_charge_C(); m.endpoint_charge_C = r0.endpoint_charge_C();
            m.budget = r0.budget(); m.work = r0.work_scales(); m.birth_scope = r0.birth_scope();
            return m;
        };
        auto with_K = [&](double dK) {
            return attempt([&] { auto m = arrays(); m.budget.kinetic_change_J += dK; m.budget.work_defect_J = m.budget.kinetic_change_J - m.budget.field_work_J;
                                 (void)hc::MotionAssembly::assemble(ctx, hc::MotionAssembly::admit(ctx, zc_census, b, {}, zcE, zcB), b, std::move(m)); });
        };
        // inconsistencies: one electron's kinetic change off by 1e-3 of its mean thermal energy (1.5 kT), the field work of a
        // field 1e-6 different (relative) from the push field, one carrier's full exchange counted with the wrong sign
        const double kT = 1.5 * eV * n::pic_e * weight;
        const double one_carrier = w.field_work_gross_J / double(std::max<std::size_t>(zc_census.size(), 1));
        struct Case { std::string name; double dK; };
        std::vector<Case> cases{{"one electron kinetic change off by 1e-3 of 1.5 kT w", 1e-3 * kT}, {"push field 1e-6 (relative) away from the work field", 1e-6 * w.field_work_gross_J},
                                {"one carrier's mean gross exchange with the wrong sign", 2 * one_carrier}};
        J rows = J::array(); bool ok = true;
        const auto control = with_K(0.);
        for (const auto& c : cases) {
            const auto a = with_K(c.dK);
            const double margin = (std::abs(r0.budget().work_defect_J + c.dK)) / w.tolerance_J;
            const bool refused = refused_with(a, "device motion budget outside the sealed guards");
            ok = ok && refused == (margin > 1);
            rows.push_back({{"case", c.name}, {"injected_J", c.dK}, {"defect_over_tolerance", margin}, {"refused", refused}, {"result", a}});
        }
        const double ulp = std::nextafter(r0.budget().kinetic_change_J, INFINITY) - r0.budget().kinetic_change_J;
        const auto aulp = with_K(ulp);
        R.record("work-guard", "end-to-end at the step of the smallest |K|/G_K: unmodified motion admitted; inconsistent kinetic changes refused exactly when "
                 "they exceed the v3.2 tolerance (margins stated)",
                 control["accepted"].get<bool>() && ok,
                 {{"step", zc_step}, {"K_J", r0.budget().kinetic_change_J}, {"W_J", r0.budget().field_work_J}, {"defect_J", r0.budget().work_defect_J},
                  {"scales", scales_json(w)}, {"control", control}, {"cases", rows},
                  {"one_ulp_of_K", {{"ulp_J", ulp}, {"result", aulp}, {"note", "below every tolerance (v3.1 1e-10|K| included); exactness of K is the native x87 recomputation and the record digests"}}}});
    }
    J summary{{"schema", "hall-v3.2-dense-v1"}, {"failed", R.failed}, {"cases", R.cases}, {"run", run}, {"configuration_sha256", S.resident->configuration_sha256()}};
    provider.reset();
    std::ofstream(out / "dense.json") << summary.dump(1) << '\n';
    return R.failed ? 1 : 0;
}


// ---------------------------------------------------------------- pre-run (field-only) and checkpoint comparison
struct FaceSpec { std::string name; int j; };   // node plane of the face; the vacuum-side normal edges are ey(i, j_edge, k)
J field_row(const hc::Model& m, hc::ResidentOwner& res, const axial::State& s, const hc::ChargeLedger& q) {
    const auto& g = m.mesh(); J r = J::object();
    long double we = 0, wb = 0;
    for (int e = 0; e < g.ne; ++e) { const long double x = (long double)s.E[e] + s.E_low[e]; we += .5L * g.M[e] * x * x; }
    for (int b = 0; b < g.nb; ++b) { const long double x = (long double)s.B[b] + s.B_low[b]; wb += .5L * g.H[b] * x * x; }
    r["W_E_J"] = double(we); r["W_B_J"] = double(wb);
    // surface field: max |E_y| on the vacuum-side y edges of each foil face over the foil extent (|x| <= 2.5 mm, z between
    // the short and feed planes)
    auto jn = [&](double y) { return int(std::lround((y + g.side / 2) / g.dy)); };
    const int i0 = int(std::lround((-.0025 + g.side / 2) / g.dx)), i1 = int(std::lround((.0025 + g.side / 2) / g.dx));
    struct F { const char* name; int jedge; };
    const F faces[4] = {{"lower_inner", jn(-.0025)}, {"lower_outer", jn(-.0027) - 1}, {"upper_inner", jn(.0025) - 1}, {"upper_outer", jn(.0027)}};
    J sf = J::object();
    for (const auto& f : faces) {
        double mx = 0; int mi = -1, mk = -1;
        for (int k = g.kshort; k <= g.kfeed; ++k) for (int i = i0; i <= i1; ++i) {
            const int e = g.ey(i, f.jedge, k); const double v = std::abs(s.E[e] + s.E_low[e]);
            if (v > mx) { mx = v; mi = i; mk = k; }
        }
        sf[f.name] = {{"max_abs_E_normal_V_m", mx}, {"i", mi}, {"k", mk}};
    }
    r["surface_field"] = sf;
    // foil temperatures (8 slices: index 2 (k - kshort) + body)
    const auto st = res.thermal().states();
    for (int body = 0; body < 2; ++body) {
        double mx = 0, sum = 0; long cnt = 0;
        for (std::size_t i = 0; i < st.size(); ++i) if (int(i % 2) == body) for (double T : st[i].T) { mx = std::max(mx, T); sum += T; ++cnt; }
        r[body ? "upper_foil" : "lower_foil"] = {{"max_T_K", mx}, {"mean_T_K", cnt ? sum / cnt : 0.}};
    }
    (void)q;
    return r;
}
int prerun(const std::string& trace, const std::filesystem::path& out, double dt, double t_end, double sample_s, double checkpoint_s, double dt_max,
           const std::string& nvrtc, const std::string& uuid) {
    auto dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    auto probe = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto desc = hc::resident_descriptor_sha256(dev, *probe);
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<hc::ResidentBackend>(desc), hc::CarrierLimits{}, hc::ThermalDeclaration{},
                                                   dt_max, hc::MotionForms{true, 1e-14});
    const auto& g = model->mesh(); const auto o = mf::geometry(g); const auto seed = model->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    hc::need(dt > 0 && dt <= dt_max && t_end > seed.time_s && t_end <= model->drive().end() && sample_s > 0, "pre-run arguments");
    hc::ResidentOwner res(model, dev, seed, declared(*model, seed, hc::charge_of_census(g, o, none)), o, none);
    std::filesystem::create_directories(out);
    std::ofstream series(out / "series.jsonl");
    const long double ratio_sample = (long double)sample_s, ratio_ck = (long double)checkpoint_s;
    // landing grid: m / (1/sample_s) as a correctly rounded quotient (exact nearest double of m * sample_s for integral 1/sample_s)
    const double per = std::round(1.0 / sample_s);
    hc::need(std::abs(per * sample_s - 1) < 1e-12, "SAMPLE_S must be the reciprocal of an integer");
    const long long m_end = std::llround(t_end * per);
    hc::need(double(m_end) / per == t_end, "T_END_S must lie on the landing grid");
    long long m = (long long)std::floor(seed.time_s * per) + 1;
    const long long ck_every = checkpoint_s > 0 ? std::llround(checkpoint_s * per) : 0;
    (void)ratio_sample; (void)ratio_ck;
    long steps = 0; double max_guard = 0, max_v31 = 0, max_divb = 0, max_energy = 0, max_resid = 0; J checkpoints = J::array();
    hc::Vec cen, bnd, raw, led, dvb, dbd;
    double max_cell_divb = 0, min_bound_seen = INFINITY; J worst_divb_row = nullptr;
    const auto t0 = steady::now();
    double V_last = 0, I_last = 0;
    auto row = [&](long long mm) {
        axial::State s; hc::ChargeLedger q; res.accepted_state(s, q);
        J r = field_row(*model, res, s, q);
        const auto cs = res.device_field().gauss_profile(0, res.accepted_meta(), cen, bnd, &raw, &led, &dvb, &dbd);
        {   // v3.3: per-cell divB / bound, the smallest bound (its floor at field zero crossings) and the ratio there
            std::size_t qm = 0, qb = 0; double rm = 0;
            for (std::size_t c = 0; c < dvb.size(); ++c) { const double rr = dbd[c] > 0 ? dvb[c] / dbd[c] : (dvb[c] > 0 ? INFINITY : 0);
                if (rr > rm) { rm = rr; qm = c; } if (dbd[c] < dbd[qb]) qb = c; }
            auto cj = [&](std::size_t c) { return J{{"cell", c}, {"i", int(c % g.nx)}, {"j", int((c / g.nx) % g.ny)}, {"k", int(c / (std::size_t(g.nx) * g.ny))},
                                                   {"abs_divB_T_m", dvb[c]}, {"bound_T_m", dbd[c]}}; };
            r["divB_cell"] = {{"max_ratio", rm}, {"at", cj(qm)}, {"min_bound", cj(qb)}, {"ratio_at_min_bound", dbd[qb] > 0 ? dvb[qb] / dbd[qb] : 0.}};
            if (rm > max_cell_divb) { max_cell_divb = rm; worst_divb_row = {{"time_s", res.time_s()}, {"max_ratio", rm}, {"at", cj(qm)}}; }
            min_bound_seen = std::min(min_bound_seen, dbd[qb]);
        }
        double v31 = 0; for (int k = 0; k < g.nn; ++k) v31 = std::max(v31, raw[k] / bnd[k]);
        max_guard = std::max(max_guard, cs.max_gauss_node_ratio); max_v31 = std::max(max_v31, v31);
        const double divb = cs.max_divB_T_m / std::max(cs.max_divB_bound_T_m, 1e-300); max_divb = std::max(max_divb, divb);
        r["m"] = mm; r["time_s"] = res.time_s(); r["steps"] = steps; r["terminal_V"] = V_last; r["drive_current_A"] = I_last;
        r["drive_current_at_t_A"] = model->drive().current(res.time_s());
        r["gauss"] = {{"guard_ratio_max_abs_D_minus_L_over_bound", cs.max_gauss_node_ratio}, {"v31_form_ratio_max_abs_D_over_bound", v31},
                      {"implied_ratio_max_abs_D_over_bound_plus_abs_L", cs.max_gauss_defect_over_bound_plus_ledger}, {"max_abs_L_C", cs.max_gauss_ledger_C},
                      {"max_abs_D_C", cs.max_gauss_defect_C}};
        r["divB_over_bound_max_max"] = divb; r["wall_s"] = ms_since(t0) / 1e3;
        series << r.dump() << '\n'; series.flush();
    };
    row(m - 1);
    for (; m <= m_end; ++m) {
        const double target = double(m) / per;
        while (res.time_s() < target) {
            const double remaining = target - res.time_s();
            long long nleft = (long long)std::ceil(remaining / dt);
            double h = remaining / double(std::max<long long>(nleft, 1));
            if (h > dt) h = remaining / double(nleft + 1);
            if (nleft <= 1 && remaining <= dt) h = remaining;   // last step lands exactly (Sterbenz: target - t exact)
            auto p = res.prepare_vacuum(h);
            Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
            hc::need(res.preflight_vacuum_final(p, b), "pre-run preflight");
            const auto& rc = p.field_receipt(); const auto& sr = p.solve_record();
            V_last = rc.terminal_V; I_last = rc.current_A;
            const double activity = std::abs(rc.field_delta_J) + std::abs(rc.plasma_work_J) + std::abs(rc.joule_J) + std::abs(rc.outer_loss_J) +
                                    std::abs(rc.outer_supply_J) + std::abs(rc.circuit_work_J);
            max_energy = std::max(max_energy, std::abs(rc.balance_J) / (1e-20 + 1e-10 * activity));
            max_resid = std::max(max_resid, sr.independent_residual / sr.target);
            res.publish_vacuum_prevalidated(std::move(p));
            ++steps;
        }
        hc::need(res.time_s() == target, "landing time not reached exactly");
        row(m);
        if ((ck_every > 0 && m % ck_every == 0) || m == m_end) {
            char name[64]; std::snprintf(name, sizeof name, "checkpoint-%.6fns", target * 1e9);
            const auto dir = out / name; std::filesystem::remove_all(dir); res.save(dir);
            checkpoints.push_back({{"time_s", target}, {"dir", dir.string()}, {"field_sha256", res.field_sha256()}, {"heat_sha256", res.heat_sha256()},
                                   {"revision", res.revision()}});
            std::ofstream(out / "checkpoints.json") << checkpoints.dump(1) << '\n';
        }
    }
    const double wall = ms_since(t0) / 1e3;
    axial::State s; hc::ChargeLedger q; res.accepted_state(s, q);
    J summary{{"schema", "hall-v3.2-prerun-v1"}, {"dt_s", dt}, {"dt_max_s", dt_max}, {"start_s", seed.time_s}, {"end_time_s", res.time_s()}, {"steps", steps},
              {"wall_s", wall}, {"ms_per_step", 1e3 * wall / std::max(1L, steps)}, {"configuration_sha256", model->configuration_sha256()},
              {"final_field_sha256", res.field_sha256()}, {"final_heat_sha256", res.heat_sha256()},
              {"final_physical_sha256", physical_sha256(s, q, res.thermal())}, {"final_ledger_sha256", ledger_sha256(q, g.nn)},
              {"max_gauss_guard_ratio", max_guard}, {"max_v31_form_gauss_ratio", max_v31}, {"max_divB_over_bound", max_divb},
              {"max_energy_over_bound", max_energy}, {"max_independent_residual_over_target", max_resid}, {"checkpoints", checkpoints},
              {"max_cell_divB_over_bound", max_cell_divb}, {"worst_cell_divB", worst_divb_row}, {"min_cell_divB_bound_T_m", min_bound_seen},
              {"seed", {{"coulomb_cycles", model->seed_coulomb_cycles()}, {"coulomb_residual", model->seed_coulomb_residual()}}},
              {"final", field_row(*model, res, s, q)}, {"device", J::parse(dev->device_json())}};
    std::ofstream(out / "prerun.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return 0;
}
int compare(const std::string& trace, const std::filesystem::path& a, const std::filesystem::path& b, const std::filesystem::path& out) {
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto& g = model->mesh();
    auto load = [&](const std::filesystem::path& d) { J j = read_json(d / "checkpoint.json"); return hc::checkpoint_detail::field(j.at("field"), g); };
    const auto A = load(a), B = load(b);
    long double sE = 0, nE = 0, sB = 0, nB = 0;
    for (int e = 0; e < g.ne; ++e) { const long double x = (long double)A.E[e] + A.E_low[e], y = (long double)B.E[e] + B.E_low[e]; sE += g.M[e] * (x - y) * (x - y); nE += g.M[e] * y * y; }
    for (int f = 0; f < g.nb; ++f) { const long double x = (long double)A.B[f] + A.B_low[f], y = (long double)B.B[f] + B.B_low[f]; sB += g.H[f] * (x - y) * (x - y); nB += g.H[f] * y * y; }
    const auto& feed = model->initial_fields().feed;
    long double va = 0, vb = 0; for (int e = 0; e < g.ne; ++e) { va -= feed[e] * ((long double)A.E[e] + A.E_low[e]); vb -= feed[e] * ((long double)B.E[e] + B.E_low[e]); }
    // ulp distances of the rounded field values (b as reference): max and the count above 4 ulp; energies
    auto ulps = [](const hc::Vec& xh, const hc::Vec& xl, const hc::Vec& yh, const hc::Vec& yl, double& mx, long& above4, double& at) {
        mx = 0; above4 = 0; at = 0;
        for (std::size_t i = 0; i < xh.size(); ++i) { const double x = xh[i] + xl[i], y = yh[i] + yl[i]; if (x == y) continue;
            const double u = std::abs(std::nextafter(y, INFINITY) - y); const double d = u > 0 ? std::abs(x - y) / u : INFINITY;
            if (d > mx) { mx = d; at = y; } if (d > 4) ++above4; } };
    double uE, uB, atE, atB; long aE, aB; ulps(A.E, A.E_low, B.E, B.E_low, uE, aE, atE); ulps(A.B, A.B_low, B.B, B.B_low, uB, aB, atB);
    long double wEa = 0, wEb = 0, wBa = 0, wBb = 0;
    for (int e = 0; e < g.ne; ++e) { const long double x = (long double)A.E[e] + A.E_low[e], y = (long double)B.E[e] + B.E_low[e]; wEa += .5L * g.M[e] * x * x; wEb += .5L * g.M[e] * y * y; }
    for (int f = 0; f < g.nb; ++f) { const long double x = (long double)A.B[f] + A.B_low[f], y = (long double)B.B[f] + B.B_low[f]; wBa += .5L * g.H[f] * x * x; wBb += .5L * g.H[f] * y * y; }
    J r{{"schema", "hall-v3.3-checkpoint-compare-v1"}, {"a", a.string()}, {"b", b.string()}, {"time_a_s", A.time_s}, {"time_b_s", B.time_s},
        {"max_ulp_E", uE}, {"edges_above_4_ulp", aE}, {"value_at_max_ulp_E", atE}, {"max_ulp_B", uB}, {"faces_above_4_ulp", aB}, {"value_at_max_ulp_B", atB},
        {"W_E_rel_diff", double(std::abs(wEa - wEb) / wEb)}, {"W_B_rel_diff", double(std::abs(wBa - wBb) / wBb)},
        {"E_energy_norm_rel_diff", double(std::sqrt(sE / nE))}, {"B_energy_norm_rel_diff", double(std::sqrt(sB / nB))},
        {"endpoint_terminal_V_a", double(va)}, {"endpoint_terminal_V_b", double(vb)}, {"endpoint_terminal_V_rel_diff", double(std::abs(va - vb) / std::max(std::abs(vb), 1e-300L))}};
    std::filesystem::create_directories(out);
    std::ofstream(out / "compare.json") << r.dump(1) << '\n';
    std::cout << r.dump(1) << std::endl;
    return 0;
}

// v3.3: DC seed construction on a mesh (CPU): Coulomb solve cycles and double-word residual, seed divergence statistics.
int seed_check(const std::string& trace, int nt, int axial_cells, int air, const std::filesystem::path& out) {
    const auto t0 = steady::now();
    J r{{"schema", "hall-v3.3-seed-check-v1"}, {"mesh", {{"transverse", nt}, {"axial", axial_cells}, {"air", air}}}};
    try {
        hc::Model m(hc::MeshSpec{nt, axial_cells, air, .009}, trace);
        const auto& g = m.mesh(); const auto s = m.dc_seed();
        long double worst = 0; double wd = 0, wb = 0; const auto B = hc::words(s.B, s.B_low);
        for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
            hc::Words d; long double scale = 0;
            for (auto [a, b, h] : std::array<std::tuple<int, int, double>, 3>{{{g.bx(i + 1, j, k), g.bx(i, j, k), g.dx}, {g.by(i, j + 1, k), g.by(i, j, k), g.dy}, {g.bz(i, j, k + 1), g.bz(i, j, k), g.dz}}}) {
                d = d + (B[a] - B[b]) / h; scale += (hc::magnitude(B[a]) + hc::magnitude(B[b])) / h; }
            const long double bound = 128 * eps_ld * scale, v = std::abs(d.wide());
            if (bound > 0 && v / bound > worst) { worst = v / bound; wd = double(v); wb = double(bound); }
        }
        r["admitted"] = true; r["coulomb_cycles"] = m.seed_coulomb_cycles(); r["coulomb_residual"] = m.seed_coulomb_residual();
        r["edges"] = g.ne; r["seed_max_divB_over_bound"] = double(worst); r["at_divB_T_m"] = wd; r["at_bound_T_m"] = wb;
        r["configuration_sha256"] = m.configuration_sha256();
    } catch (const std::exception& e) { r["admitted"] = false; r["refusal"] = e.what(); }
    r["wall_s"] = ms_since(t0) / 1e3;
    std::filesystem::create_directories(out);
    std::ofstream(out / "seed-check.json") << r.dump(1) << '\n';
    std::cout << r.dump(1) << std::endl;
    return r["admitted"].get<bool>() ? 0 : 1;
}

int restore_continue(const std::string& trace, const std::filesystem::path& ck, const std::filesystem::path& reference, const std::filesystem::path& out,
                     const std::string& nvrtc, const std::string& uuid) {
    Results R; auto S = setup(trace, nvrtc, uuid);
    std::filesystem::create_directories(out);
    const auto& g = S.resident->mesh(); const auto o = mf::geometry(g); std::vector<n::ArenaRateEntry> none;
    const J ref = read_json(reference);
    auto r = hc::ResidentOwner::restore(S.resident, S.dev, o, none, ck);
    const bool start = r->field_sha256() == ref.at("saved_field_sha256").get<std::string>() && r->heat_sha256() == ref.at("saved_heat_sha256").get<std::string>();
    const auto c = continue_steps(*r, ref.at("steps").get<long>(), ref.at("dt_s").get<double>());
    const bool same = J(c.field) == ref.at("field_sha256") && J(c.heat) == ref.at("heat_sha256") && c.physical == ref.at("physical_sha256").get<std::string>() &&
                      c.ledger == ref.at("ledger_sha256").get<std::string>();
    R.record("checkpoint-resident", "fresh-process restore of the resident checkpoint: digests equal the saved ones and the continuation is bit-identical to the "
             "uninterrupted run", start && same, {{"steps", c.field.size()}, {"final_physical_sha256", c.physical}, {"final_ledger_sha256", c.ledger}});
    J summary{{"schema", "hall-v3.2-restore-continue-v1"}, {"failed", R.failed}, {"cases", R.cases}, {"checkpoint", ck.string()}};
    std::ofstream(out / "restore-continue.json") << summary.dump(1) << '\n';
    return R.failed ? 1 : 0;
}
}  // namespace

int main(int argc, char** argv) {
    mallopt(M_MMAP_THRESHOLD, 256 << 20); mallopt(M_TRIM_THRESHOLD, 512 << 20);
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "a09" && argc == 6) return a09(argv[2], argv[3], argv[4], argv[5]);
        if (mode == "ledger-run" && argc == 9) return ledger_run(argv[2], argv[3], std::stol(argv[4]), std::stol(argv[5]), std::stol(argv[6]), argv[7], argv[8]);
        if (mode == "restore-continue" && argc == 8) return restore_continue(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
        if (mode == "dense" && argc >= 11 && argc <= 13)
            return dense(argv[2], argv[3], std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stol(argv[7]), std::string(argv[8]) == "gpu", argv[9], argv[10],
                         argc >= 12 && std::string(argv[11]) == "v7", argc == 13 && std::string(argv[12]) == "sinks");
        if (mode == "prerun" && argc == 11)
            return prerun(argv[2], argv[3], std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), argv[9], argv[10]);
        if (mode == "compare" && argc == 6) return compare(argv[2], argv[3], argv[4], argv[5]);
        if (mode == "seed-check" && argc == 7) return seed_check(argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), argv[6]);
        std::cerr << "usage: see the header of test_v32.cpp\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
}
