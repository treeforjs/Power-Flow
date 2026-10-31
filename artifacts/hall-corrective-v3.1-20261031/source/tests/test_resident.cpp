// Hall corrective v3.1 device-resident field step: parity, negative controls and timing.
//   test_resident units          TRACE OUT NVRTC UUID
//   test_resident parity-vacuum  TRACE OUT STEPS NVRTC UUID
//   test_resident parity-coupled TRACE OUT STEPS CARRIERS NVRTC UUID
//   test_resident negative       TRACE OUT NVRTC UUID
//   test_resident vacuum-run     TRACE OUT STEPS NVRTC UUID [verify-all] [LEDGER_EVERY]
//   test_resident coupled-run    TRACE OUT STEPS CARRIERS NVRTC UUID [cpu|gpu]
//   test_resident gauss-trace    TRACE OUT STEPS EVERY NVRTC UUID [DT]   (diagnostic)
//   test_resident dt-study       TRACE OUT NVRTC UUID WINDOW_S DT1 DT2 ... (finest first; diagnostic)
// Parity: the reference is the v3.1 host Owner (Hall v3 arithmetic; cpu-pcg or gpu1 backend),
// run in lockstep in the same process; physical state compared word for word every step.
#include "resident_owner.hpp"
#include "gpu/gpu_backends.hpp"
#include "hall-finite-foil-em-20261029/axial_particles.hpp"
namespace hc = hall_corrective_v3;
#include "motion_fixtures.hpp"
#include <algorithm>
#include <malloc.h>
#include <fstream>
#include <iostream>

namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
namespace gp = hall_corrective_v3::gpu;
using J = hc::Json;
using steady = std::chrono::steady_clock;

namespace {
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
struct Stat {
    std::vector<double> v;
    J json() const {
        if (v.empty()) return J::object();
        auto s = v; std::sort(s.begin(), s.end()); double sum = 0; for (double x : s) sum += x;
        return {{"n", s.size()}, {"median_ms", s[s.size() / 2]}, {"mean_ms", sum / s.size()}, {"p10_ms", s[s.size() / 10]}, {"p90_ms", s[s.size() * 9 / 10]},
                {"min_ms", s.front()}, {"max_ms", s.back()}};
    }
};
struct Stats {
    std::map<std::string, Stat> m;
    void add(const std::string& k, double x) { m[k].v.push_back(x); }
    J json() const { J j = J::object(); for (const auto& [k, s] : m) j[k] = s.json(); return j; }
};
struct Boundary { std::string field_sha256, heat_sha256, source_sha256, mesh_sha256, material_sha256, common_clock_sha256; double time_s = 0; std::uint64_t accepted_step = 0; };
hc::ChargeLedger declared(const hc::Model& m, const axial::State& s, hc::FreeCharge particles) {
    const auto& g = m.mesh(); hc::ChargeLedger q; q.particles = std::move(particles);
    q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
    q.conductor_declaration = hc::sha256("parity: declared DC conductor/PEC induced charge");
    q.exterior_declaration = hc::sha256("parity: declared DC exterior termination charge");
    const auto field = hc::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = field[k]; else if (g.outer_node[k]) q.exterior_C[k] = field[k]; }
    return q;
}
hall_driver::MotionBinding binding_for(double t, std::uint64_t rev, const std::string& mesh, const std::string& material, const std::string& policy, long step,
                                       double dt = 1e-14) {
    hall_driver::MotionBinding b; b.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = rev; b.staged_arena_revision = rev + 1; b.begin_s = t; b.dt_s = dt; b.end_s = b.begin_s + b.dt_s;
    b.staged_arena_sha256 = hc::sha256("resident/staged/" + std::to_string(step)); b.mesh_sha256 = mesh;
    b.material_sha256 = material; b.species_sha256 = hc::sha256("resident/species");
    b.common_clock_sha256 = hc::sha256("resident/clock/" + std::to_string(step)); b.kinetic_policy_sha256 = policy;
    return b;
}
bool same_words(const hc::WordVec& a, const hc::WordVec& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (a[i].hi != b[i].hi || a[i].lo != b[i].lo) return false;
    return true;
}
bool same_bits(const hc::Vec& a, const hc::Vec& b) { return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * 8) == 0); }
// Physical state comparison (words, clock, complete ledger).
J compare_state(const axial::State& a, const hc::ChargeLedger& qa, const axial::State& b, const hc::ChargeLedger& qb) {
    J j = J::object();
    j["E"] = same_bits(a.E, b.E); j["B"] = same_bits(a.B, b.B); j["E_low"] = same_bits(a.E_low, b.E_low); j["B_low"] = same_bits(a.B_low, b.B_low);
    j["time"] = a.time_s == b.time_s; j["revision"] = a.revision == b.revision;
    j["free_net"] = same_words(qa.particles.net_C, qb.particles.net_C); j["free_gross"] = same_bits(qa.particles.gross_C, qb.particles.gross_C);
    j["conductor"] = same_words(qa.conductor_C, qb.conductor_C); j["exterior"] = same_words(qa.exterior_C, qb.exterior_C);
    j["residue"] = same_bits(qa.residue_C, qb.residue_C);
    bool all = true; for (auto& [k, v] : j.items()) all = all && v.get<bool>();
    j["all"] = all;
    return j;
}
double rel(double a, double b) { const double s = std::max(std::abs(a), std::abs(b)); return s == 0 ? 0 : std::abs(a - b) / s; }
// Guard sums: the host value is a sequential binary64 sum (error <= (n-1) eps sum|t|), the device
// value a double-word fixed-tree sum of the same rounded terms. Returns |device - host| / ((n-1) eps sum|t|)
// for each sum (<= 1 means the difference is within the host's own rounding bound).
J guard_sum_check(const hc::Model& m, const axial::State& old, const axial::State& next, const axial::Receipt& host, const axial::Receipt& dev, const hc::Vec& plasma) {
    const auto& g = m.mesh(); const auto& in = m.initial_fields(); const double dt = next.time_s - old.time_s;
    const double eps = std::numeric_limits<double>::epsilon();
    long double a_pw = 0, a_j = 0, a_c = 0, a_ol = 0, a_os = 0, a_t = 0, a_fd = 0;
    const auto& mid = host.mid_E;
    for (int e = 0; e < g.ne; ++e) {
        a_pw += std::abs(dt * mid[e] * plasma[e]); a_j += std::abs(dt * mid[e] * host.solid_current[e]); a_c += std::abs(dt * host.current_A * in.feed[e] * mid[e]);
        a_ol += std::abs(dt * g.Y[e] * mid[e] * mid[e]); a_os += std::abs(dt * in.background_force[e] * mid[e]); a_t += std::abs((long double)in.feed[e] * mid[e]);
        const auto o = hc::Words(old.E[e], old.E_low[e]), nx = hc::Words(next.E[e], next.E_low[e]);
        a_fd += std::abs((long double)g.M[e] * ((nx + o) * .5).wide() * (nx - o).wide());
    }
    for (int b = 0; b < g.nb; ++b) {
        const auto o = hc::Words(old.B[b], old.B_low[b]), nx = hc::Words(next.B[b], next.B_low[b]);
        a_fd += std::abs((long double)g.H[b] * ((nx + o) * .5).wide() * (nx - o).wide());
    }
    auto ratio = [&](double h, double d, long double abs_sum, int n) { const long double bound = (n - 1) * eps * abs_sum + 1e-300L; return double(std::abs((long double)h - d) / bound); };
    J j = J::object();
    j["plasma_work_J"] = ratio(host.plasma_work_J, dev.plasma_work_J, a_pw, g.ne); j["joule_J"] = ratio(host.joule_J, dev.joule_J, a_j, g.ne);
    j["circuit_work_J"] = ratio(host.circuit_work_J, dev.circuit_work_J, a_c, g.ne); j["outer_loss_J"] = ratio(host.outer_loss_J, dev.outer_loss_J, a_ol, g.ne);
    j["outer_supply_J"] = ratio(host.outer_supply_J, dev.outer_supply_J, a_os, g.ne); j["terminal_V"] = ratio(host.terminal_V, dev.terminal_V, a_t, g.ne);
    j["field_delta_J"] = ratio(host.field_delta_J, dev.field_delta_J, a_fd, g.ne + g.nb);
    double worst = 0; for (auto& [k, v] : j.items()) worst = std::max(worst, v.get<double>());
    j["worst"] = worst;
    return j;
}
J receipt_diff(const axial::Receipt& a, const axial::Receipt& b) {
    J j = J::object();
    j["field_delta_J"] = rel(a.field_delta_J, b.field_delta_J); j["plasma_work_J"] = rel(a.plasma_work_J, b.plasma_work_J);
    j["joule_J"] = rel(a.joule_J, b.joule_J); j["circuit_work_J"] = rel(a.circuit_work_J, b.circuit_work_J);
    j["outer_loss_J"] = rel(a.outer_loss_J, b.outer_loss_J); j["outer_supply_J"] = rel(a.outer_supply_J, b.outer_supply_J);
    j["terminal_V"] = rel(a.terminal_V, b.terminal_V); j["balance_J_host"] = a.balance_J; j["balance_J_device"] = b.balance_J;
    j["gauss_increment_C"] = rel(a.gauss_increment_C, b.gauss_increment_C);
    j["kcl_residual_A_bitwise"] = a.kcl_residual_A == b.kcl_residual_A; j["ampere_residual_bitwise"] = a.ampere_residual_A_m == b.ampere_residual_A_m;
    j["region_conduction_A"] = rel(a.region_conduction_A, b.region_conduction_A); j["region_displacement_A"] = rel(a.region_displacement_A, b.region_displacement_A);
    j["cell_joule_W_bitwise"] = same_bits(a.cell_joule_W, b.cell_joule_W);
    j["current_A_bitwise"] = a.current_A == b.current_A;
    return j;
}
struct Setup {
    std::shared_ptr<gp::Device> dev;
    std::shared_ptr<const hc::Model> cpu, gpu1, resident, resident_cpu;
};
Setup setup(const std::string& trace, const std::string& nvrtc, const std::string& uuid, hc::CarrierLimits limits = {}, bool with_gpu1 = true, double dt_max = 1e-14) {
    Setup s;
    s.dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    s.cpu = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, hc::make_cpu_pcg(), limits, hc::ThermalDeclaration{}, dt_max);
    if (with_gpu1) s.gpu1 = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, gp::make_gpu_field_solver(s.dev), limits, hc::ThermalDeclaration{}, dt_max);
    const auto desc = hc::resident_descriptor_sha256(s.dev, *s.cpu);
    s.resident = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<hc::ResidentBackend>(desc), limits, hc::ThermalDeclaration{}, dt_max);
    return s;
}
std::string state_values_sha256(const axial::State& s, const hc::ChargeLedger& q, const hc::ThermalV3& th) {
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
J control_json(const gp::ControlReceipt& c) {
    return {{"omega_p_dt_max", c.omega_p_dt_max}, {"B_max_T", c.B_max_T}, {"omega_c_dt_max", c.omega_c_dt_max}, {"drive_current_A", c.drive_current_A},
            {"drive_rel_change", c.drive_rel_change}, {"drive_knot_interval_s", c.drive_knot_interval_s}, {"increment_rel", c.increment_rel},
            {"third_difference_rel", c.third_difference_rel}, {"history", c.history}};
}

// ---------------------------------------------------------------- units
int units(const std::string& trace, const std::filesystem::path& out, const std::string& nvrtc, const std::string& uuid) {
    Results R;
    auto S = setup(trace, nvrtc, uuid, {}, false);
    const auto& g = S.cpu->mesh();
    // HostApply equals the sealed apply lambda bit for bit (random vectors, heated conductivity).
    {
        std::mt19937_64 rng(5); std::normal_distribution<double> nd;
        hc::ThermalV3 th(g, S.cpu->drive().start());
        auto sigma = th.conductivity(); for (std::size_t c = 0; c < sigma.size(); ++c) if (g.cells[c].body >= 0) sigma[c] *= 1 + 0.01 * nd(rng);
        const auto S_ref = g.conduction(sigma);
        std::vector<int> metal; for (std::size_t c = 0; c < g.cells.size(); ++c) if (g.cells[c].body >= 0) metal.push_back(int(c));
        hc::Vec S_fast; gp::conduction_metal(g, metal, sigma, S_fast);
        R.record("units", "conduction from metal cells equals Mesh::conduction bit for bit", same_bits(S_ref, S_fast));
        gp::HostApply ha(g); bool all = true; const double dt = 1e-14;
        hc::Vec sqrtM(g.ne); for (int e = 0; e < g.ne; ++e) sqrtM[e] = std::sqrt(g.M[e]);
        auto apply = [&](const hc::Vec& v) {
            hc::Vec E(g.ne);
            for (int e = 0; e < g.ne; ++e) if (!g.fixed_E[e]) E[e] = v[e] / sqrtM[e];
            auto K = g.curl_adjoint(g.curl(E));
            hc::Vec result(g.ne);
            for (int e = 0; e < g.ne; ++e) result[e] = g.fixed_E[e] ? v[e] : v[e] + (.5 * dt * (S_ref[e] + g.Y[e]) * E[e] + .25 * dt * dt * K[e]) / sqrtM[e];
            return result;
        };
        for (int k = 0; k < 5; ++k) {
            hc::Vec v(g.ne); for (double& x : v) x = nd(rng) * std::pow(10., k - 2);
            hc::Vec out; ha.apply(S_fast, dt, v, out);
            all = all && same_bits(out, apply(v));
        }
        R.record("units", "host apply of the independent residual equals the sealed apply bit for bit (5 random vectors)", all);
    }
    // Device leaves equal host SHA-256 of every block (full verification on a loaded state).
    {
        auto seed = S.cpu->dc_seed(); std::vector<n::ArenaRateEntry> none; const auto o = mf::geometry(g);
        auto q = declared(*S.cpu, seed, hc::charge_of_census(g, o, none));
        gp::ResidentOptions opt; opt.verify_all = true;
        hc::ResidentOwner owner(S.resident, S.dev, seed, q, o, none, opt);
        R.record("units", "bootstrap digest with every leaf re-hashed on the host", hc::digest(owner.field_sha256()),
                 {{"field_sha256", owner.field_sha256()}, {"stats_leaves_host", owner.device_field().digest_stats().leaves_hashed_on_host}});
        axial::State s2; hc::ChargeLedger q2; owner.accepted_state(s2, q2);
        R.record("units", "upload/download round trip of the accepted state is exact", compare_state(seed, q, s2, q2)["all"].get<bool>());
        // admission refusals on the device follow the sealed texts
        auto bad = seed; bad.E_low[100] = bad.E[100] == 0 ? 1e-300 : bad.E[100];
        std::string why;
        try { hc::ResidentOwner o2(S.resident, S.dev, bad, q, o, none); } catch (const std::exception& e) { why = e.what(); }
        R.record("units", "noncanonical E words refused (host bootstrap admission)", why == "noncanonical field words", {{"refusal", why}});
    }
    J summary{{"schema", "hall-v3.1-units-v1"}, {"failed", R.failed}, {"cases", R.cases}, {"device", J::parse(S.dev->device_json())}};
    std::filesystem::create_directories(out);
    std::ofstream(out / "units.json") << summary.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// ---------------------------------------------------------------- vacuum parity
int parity_vacuum(const std::string& trace, const std::filesystem::path& out, long steps, const std::string& nvrtc, const std::string& uuid) {
    Results R;
    auto S = setup(trace, nvrtc, uuid);
    const auto& g = S.cpu->mesh(); const auto o = mf::geometry(g); const auto seed = S.cpu->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    const auto q0 = declared(*S.cpu, seed, hc::charge_of_census(g, o, none));
    J rows = J::array();
    for (int pass = 0; pass < 2; ++pass) {
        const bool device_cg = pass == 1;
        auto ref_model = device_cg ? S.gpu1 : S.cpu;
        hc::Owner ref(ref_model, seed, q0, o, none);
        gp::ResidentOptions opt; opt.solve = device_cg ? gp::ResidentOptions::Solve::device_cg : gp::ResidentOptions::Solve::host_cpu_pcg;
        opt.verify_all = true;
        hc::ResidentOwner res(S.resident, S.dev, seed, q0, o, none, opt);
        bool all_state = true, all_heat = true, all_mid = true, all_joule = true, all_solve = true, all_bitwise_guards = true, ctl_ok = true;
        double worst_receipt = 0, worst_incr = 0; J worst = J::object(); J first_bad = nullptr; J ctl_last = nullptr;
        for (long step = 1; step <= steps; ++step) {
            auto f = ref.prepare_vacuum(1e-14);
            Boundary bf{f.candidate_field_sha256(), f.candidate_heat_sha256(), "", "", "", "", f.end_s(), f.endpoint_revision()};
            hc::need(ref.preflight_vacuum_final(f, bf), "reference vacuum preflight");
            auto p = res.prepare_vacuum(1e-14);
            Boundary bp{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
            const bool pre = res.preflight_vacuum_final(p, bp);
            // candidate comparison
            axial::State cs; hc::ChargeLedger cq; res.candidate_state(cs, cq, p);
            const auto cmp = compare_state(f.candidate_field(), f.candidate_charge(), cs, cq);
            hc::Vec mE, mB; res.device_field().midpoint_download(mE, mB);
            const bool mid = same_bits(mE, f.midpoint().E()) && same_bits(mB, f.midpoint().B());
            const auto rd = receipt_diff(f.field_receipt(), p.field_receipt());
            const bool joule = rd["cell_joule_W_bitwise"].get<bool>();
            const bool solve_same = f.solve_record().rhs_norm == p.solve_record().rhs_norm && f.solve_record().independent_residual == p.solve_record().independent_residual &&
                                    f.solve_record().iterations == p.solve_record().iterations;
            const bool bguards = rd["kcl_residual_A_bitwise"].get<bool>() && rd["ampere_residual_bitwise"].get<bool>();
            {
                const auto gs = guard_sum_check(*S.cpu, ref.field(), f.candidate_field(), f.field_receipt(), p.field_receipt(), hc::Vec(g.ne));
                const double v = gs["worst"].get<double>();
                if (v >= worst_receipt) { worst_receipt = v; worst = {{"step", step}, {"ratio_to_host_rounding_bound", gs}, {"relative", rd}}; }
            }
            {   // controller receipt: device values against the host form on the reference proposal
                const auto hcr = hc::host_control_receipt(g, S.cpu->drive(), ref.field(), f.candidate_field(), f.candidate_charge().particles, f.midpoint().B(), 1e-14);
                const auto& dcr = p.control();
                ctl_ok = ctl_ok && hcr.omega_p_dt_max == dcr.omega_p_dt_max && hcr.B_max_T == dcr.B_max_T && hcr.omega_c_dt_max == dcr.omega_c_dt_max &&
                         hcr.drive_rel_change == dcr.drive_rel_change && hcr.drive_knot_interval_s == dcr.drive_knot_interval_s && hcr.drive_current_A == dcr.drive_current_A;
                worst_incr = std::max(worst_incr, rel(hcr.increment_rel, dcr.increment_rel));
                ctl_last = {{"host", control_json(hcr)}, {"device", control_json(dcr)}};
            }
            ref.publish_vacuum_prevalidated(std::move(f));
            hc::need(pre, "resident vacuum preflight");
            res.publish_vacuum_prevalidated(std::move(p));
            const bool heat = hc::heat_sha256_of(S.cpu->material_sha256(), ref.thermal()) == res.heat_sha256();
            all_state = all_state && cmp["all"].get<bool>(); all_heat = all_heat && heat; all_mid = all_mid && mid; all_joule = all_joule && joule;
            all_solve = all_solve && (device_cg ? solve_same : f.solve_record().independent_residual == p.solve_record().independent_residual);
            all_bitwise_guards = all_bitwise_guards && bguards;
            if (first_bad.is_null() && !(cmp["all"].get<bool>() && heat && mid && joule)) first_bad = {{"step", step}, {"state", cmp}, {"heat", heat}, {"mid", mid}, {"receipt", rd}};
            if (step == 1 || step == steps || step % std::max(1L, steps / 5) == 0)
                rows.push_back({{"pass", device_cg ? "gpu1-vs-resident-device-cg" : "cpu-pcg-vs-resident-host-pcg"}, {"step", step}, {"state", cmp}, {"heat_equal", heat},
                                {"midpoint_equal", mid}, {"receipt", rd}, {"solve_ref", {{"iterations", ref.revision()}}},
                                {"control", control_json(gp::ControlReceipt{})}});
        }
        axial::State fs; hc::ChargeLedger fq; res.accepted_state(fs, fq);
        const auto ev = J{{"steps", steps}, {"first_mismatch", first_bad}, {"worst_receipt_relative", worst},
                          {"reference_physical_sha256", state_values_sha256(ref.field(), ref.charge(), ref.thermal())},
                          {"resident_physical_sha256", state_values_sha256(fs, fq, res.thermal())}, {"resident_field_sha256", res.field_sha256()},
                          {"reference_field_sha256", ref.field_sha256()}, {"heat_sha256", res.heat_sha256()},
                          {"leaves_hashed_on_host", res.device_field().digest_stats().leaves_hashed_on_host}};
        const std::string tag = device_cg ? "gpu1 (host arithmetic) vs resident (device arithmetic, device CG)"
                                          : "cpu-pcg (host arithmetic) vs resident (device arithmetic, sealed CPU PCG on the device rhs)";
        R.record("parity-vacuum", tag + ": E, B, E_low, B_low, clock and complete ledger bit-identical every step", all_state, ev);
        R.record("parity-vacuum", tag + ": heat state digest (unchanged v3 record) identical every step", all_heat);
        R.record("parity-vacuum", tag + ": authoritative midpoint bit-identical every step", all_mid);
        R.record("parity-vacuum", tag + ": cell Joule power bit-identical every step", all_joule);
        R.record("parity-vacuum", tag + ": solve record (rhs norm, independent residual, iterations) identical", all_solve);
        R.record("parity-vacuum", tag + ": KCL and Ampere residual maxima bit-identical", all_bitwise_guards);
        R.record("parity-vacuum", tag + ": guard sums differ from the host sequential sums by less than the host rounding bound (n-1) eps sum|t|",
                 worst_receipt <= 1, worst);
        R.record("parity-vacuum", tag + ": controller receipt (omega_p dt, |B| max, omega_c dt, drive change) equal to the host form; increment norm within 1e-9",
                 ctl_ok && worst_incr <= 1e-9, {{"increment_rel_relative_difference", worst_incr}, {"last", ctl_last}});
    }
    J summary{{"schema", "hall-v3.1-parity-vacuum-v1"}, {"failed", R.failed}, {"cases", R.cases}, {"rows", rows}, {"device", J::parse(S.dev->device_json())},
              {"teardown", J::parse(S.dev->close_json())}};
    std::filesystem::create_directories(out);
    std::ofstream(out / "parity-vacuum.json") << summary.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// ---------------------------------------------------------------- coupled parity
int parity_coupled(const std::string& trace, const std::filesystem::path& out, long steps, std::size_t carriers, const std::string& nvrtc, const std::string& uuid) {
    Results R;
    hc::CarrierLimits lim{std::max<std::size_t>(100000, carriers), 1ull << 34};
    auto S = setup(trace, nvrtc, uuid, lim);
    const auto& g = S.cpu->mesh(); const auto o = mf::geometry(g); const auto seed = S.cpu->dc_seed();
    const auto census0 = carriers == 2 ? mf::pair({g.side / 2, g.side / 2, 0.}, 1.) : mf::pairs(carriers, 7, 1e3);
    const auto q0 = declared(*S.cpu, seed, hc::charge_of_census(g, o, census0, lim.maximum_carriers));
    J rows = J::array();
    for (int pass = 0; pass < 2; ++pass) {
        const bool device_cg = pass == 1;
        hc::Owner ref(device_cg ? S.gpu1 : S.cpu, seed, q0, o, census0);
        gp::ResidentOptions opt; opt.solve = device_cg ? gp::ResidentOptions::Solve::device_cg : gp::ResidentOptions::Solve::host_cpu_pcg;
        opt.verify_all = true;
        hc::ResidentOwner res(S.resident, S.dev, seed, q0, o, census0, opt);
        std::unique_ptr<hc::MotionProvider> pr_ref, pr_res;
        if (device_cg) { pr_ref = std::make_unique<gp::GpuMotionProvider>(S.dev, ref.motion_context(), carriers); pr_res = std::make_unique<gp::GpuMotionProvider>(S.dev, res.motion_context(), carriers); }
        else { pr_ref = std::make_unique<hc::CpuMotionProvider>(); pr_res = std::make_unique<hc::CpuMotionProvider>(); }
        auto census = census0;
        bool all_state = true, all_heat = true, all_mid = true, all_motion = true, all_iter = true;
        J first_bad = nullptr;
        const std::string src = hc::sha256("resident/source");
        for (long step = 1; step <= steps; ++step) {
            auto bref = binding_for(ref.time_s(), ref.revision(), ref.motion_mesh_sha256(), S.cpu->material_sha256(), ref.motion_context().policy_sha256, step);
            auto p1 = hc::prepare_closed_proposal(ref, o, census, bref, src, *pr_ref, hc::StepInputs{});
            const auto& b1 = p1.motion().binding();
            Boundary bd1{p1.candidate_field_sha256(), p1.candidate_heat_sha256(), src, b1.mesh_sha256, b1.material_sha256, b1.common_clock_sha256, p1.end_s(), p1.endpoint_revision()};
            hc::need(ref.preflight_final(p1, census, bd1), "reference closed preflight");
            auto p2 = hc::prepare_closed_proposal(res, o, census, bref, src, *pr_res, hc::StepInputs{});
            const auto& b2 = p2.motion().binding();
            Boundary bd2{p2.candidate_field_sha256(), p2.candidate_heat_sha256(), src, b2.mesh_sha256, b2.material_sha256, b2.common_clock_sha256, p2.end_s(), p2.endpoint_revision()};
            const bool pre = res.preflight_final(p2, census, bd2);
            axial::State cs; hc::ChargeLedger cq; res.candidate_state(cs, cq, p2.field_heat());
            const auto cmp = compare_state(p1.field_heat().candidate_field(), p1.field_heat().candidate_charge(), cs, cq);
            const bool mid = same_bits(p1.midpoint().E(), p2.midpoint_E()) && same_bits(p1.midpoint().B(), p2.midpoint_B());
            const bool motion = p1.motion().record_v3_sha256() == p2.motion().record_v3_sha256() && p1.motion().seal_sha256() == p2.motion().seal_sha256() &&
                                p1.motion().census_sha256() == p2.motion().census_sha256() && b1.trial_field_sha256 == b2.trial_field_sha256;
            const bool iters = p1.iterations() == p2.iterations();
            const auto rd = receipt_diff(p1.field_heat().field_receipt(), p2.field_heat().field_receipt());
            auto next = hc::next_census(census, p1.motion());
            ref.publish_prevalidated(std::move(p1));
            hc::need(pre, "resident closed preflight");
            res.publish_prevalidated(std::move(p2));
            const bool heat = hc::heat_sha256_of(S.cpu->material_sha256(), ref.thermal()) == res.heat_sha256();
            all_state = all_state && cmp["all"].get<bool>(); all_heat = all_heat && heat; all_mid = all_mid && mid; all_motion = all_motion && motion; all_iter = all_iter && iters;
            if (first_bad.is_null() && !(cmp["all"].get<bool>() && heat && mid && motion && iters))
                first_bad = {{"step", step}, {"state", cmp}, {"heat", heat}, {"mid", mid}, {"motion", motion}, {"iterations", iters}};
            if (step == 1 || step == steps) rows.push_back({{"pass", pass}, {"step", step}, {"state", cmp}, {"receipt", rd}});
            census = std::move(next);
        }
        const std::string tag = std::string(device_cg ? "gpu1 + GPU push (host field arithmetic) vs resident + GPU push" : "cpu-pcg + CPU push vs resident (sealed CPU PCG) + CPU push") +
                                ", " + std::to_string(carriers) + " carriers";
        R.record("parity-coupled", tag + ": accepted field words and ledger bit-identical every step", all_state, {{"first_mismatch", first_bad}, {"steps", steps}});
        R.record("parity-coupled", tag + ": heat state identical every step", all_heat);
        R.record("parity-coupled", tag + ": fixed-point midpoint and iteration counts identical", all_mid && all_iter);
        R.record("parity-coupled", tag + ": motion records (v3 digest, v2 seal, census, trial field digest) identical", all_motion);
    }
    J summary{{"schema", "hall-v3.1-parity-coupled-v1"}, {"failed", R.failed}, {"cases", R.cases}, {"rows", rows}, {"teardown", J::parse(S.dev->close_json())}};
    std::filesystem::create_directories(out);
    std::ofstream(out / ("parity-coupled-" + std::to_string(carriers) + ".json")) << summary.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// ---------------------------------------------------------------- negative controls
int negative(const std::string& trace, const std::filesystem::path& out, const std::string& nvrtc, const std::string& uuid) {
    Results R;
    auto S = setup(trace, nvrtc, uuid, {}, false);
    const auto& g = S.cpu->mesh(); const auto o = mf::geometry(g); const auto seed = S.cpu->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    const auto q0 = declared(*S.cpu, seed, hc::charge_of_census(g, o, none));
    gp::ResidentOptions opt; opt.verify_samples = 16;
    hc::ResidentOwner res(S.resident, S.dev, seed, q0, o, none, opt);
    for (int i = 0; i < 3; ++i) {
        auto p = res.prepare_vacuum(1e-14);
        Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
        hc::need(res.preflight_vacuum_final(p, b), "negative setup preflight");
        res.publish_vacuum_prevalidated(std::move(p));
    }
    auto tamper = [&](int which) {   // 0: one ulp of one candidate E word, 1: one accepted ledger word
        auto p = res.prepare_vacuum(1e-14);
        axial::State s; hc::ChargeLedger q;
        if (which == 0) { res.candidate_state(s, q, p); s.E[20000] = std::nextafter(s.E[20000], INFINITY); res.device_field().upload(1, s, q); }
        else { res.accepted_state(s, q); q.conductor_C[500].lo = std::nextafter(q.conductor_C[500].lo, INFINITY); res.device_field().upload(0, s, q); }
        Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
        return res.preflight_vacuum_final(p, b);
    };
    R.record("negative", "one-ulp change of one candidate E word after sealing is refused at preflight (leaves recomputed)", !tamper(0));
    R.record("negative", "one-ulp change of one accepted ledger word is refused at preflight", !tamper(1));
    // stale proposal: a second prepare invalidates the first
    {
        axial::State s; hc::ChargeLedger q;
        hc::ResidentOwner r2(S.resident, S.dev, seed, q0, o, none, opt);
        auto a = r2.prepare_vacuum(1e-14);
        auto b2 = r2.prepare_vacuum(1e-14);
        Boundary ba{a.candidate_field_sha256(), a.candidate_heat_sha256(), "", "", "", "", a.end_s(), a.endpoint_revision()};
        R.record("negative", "a proposal whose candidate slot was overwritten by a later prepare is refused", !r2.preflight_vacuum_final(a, ba));
        (void)b2;
    }
    // dt above the configured bound is refused with the sealed text
    {
        std::string why;
        try { (void)res.prepare_vacuum(1.5e-14); } catch (const std::exception& e) { why = e.what(); }
        R.record("negative", "dt above the configured 1e-14 bound refused", why == "authoritative field interval", {{"refusal", why}});
    }
    J summary{{"schema", "hall-v3.1-negative-v1"}, {"failed", R.failed}, {"cases", R.cases}, {"teardown", J::parse(S.dev->close_json())}};
    std::filesystem::create_directories(out);
    std::ofstream(out / "negative.json") << summary.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// ---------------------------------------------------------------- timing runs
int vacuum_run(const std::string& trace, const std::filesystem::path& out, long steps, const std::string& nvrtc, const std::string& uuid, bool verify_all, long every) {
    auto S = setup(trace, nvrtc, uuid, {}, false);
    const auto& g = S.cpu->mesh(); const auto o = mf::geometry(g); const auto seed = S.cpu->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    gp::ResidentOptions opt; opt.verify_all = verify_all; opt.sample_seed = 0;
    hc::ResidentOwner res(S.resident, S.dev, seed, declared(*S.cpu, seed, hc::charge_of_census(g, o, none)), o, none, opt);
    std::filesystem::create_directories(out);
    std::ofstream ledger(out / "vacuum-ledger.jsonl");
    Stats T; hc::Digest chain("hall-v3.1-vacuum-run-chain-v1");
    double max_gauss_ratio = 0, max_divB_ratio = 0, max_energy_ratio = 0, max_resid_ratio = 0, max_heat_balance = 0;
    const auto t_all = steady::now();
    for (long step = 1; step <= steps; ++step) {
        auto t = steady::now();
        auto p = res.prepare_vacuum(1e-14);
        const double prep = ms_since(t); t = steady::now();
        Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
        const bool ok = res.preflight_vacuum_final(p, b);
        const double pre = ms_since(t);
        hc::need(ok, "vacuum preflight");
        const auto& r = p.field_receipt(); const auto& c = p.candidate_constraints(); const auto& sr = p.solve_record();
        max_gauss_ratio = std::max(max_gauss_ratio, c.max_gauss_C / std::max(c.max_gauss_bound_C, 1e-300));
        max_divB_ratio = std::max(max_divB_ratio, c.max_divB_T_m / std::max(c.max_divB_bound_T_m, 1e-300));
        const double activity = std::abs(r.field_delta_J) + std::abs(r.plasma_work_J) + std::abs(r.joule_J) + std::abs(r.outer_loss_J) + std::abs(r.outer_supply_J) + std::abs(r.circuit_work_J);
        max_energy_ratio = std::max(max_energy_ratio, std::abs(r.balance_J) / (1e-20 + 1e-10 * activity));
        max_resid_ratio = std::max(max_resid_ratio, sr.independent_residual / sr.target);
        max_heat_balance = std::max(max_heat_balance, std::abs(p.field_heat_balance_J()));
        for (const auto& [k, v] : p.step().host_ms) T.add("prepare." + k, v);
        T.add("step.prepare_vacuum", prep); T.add("step.preflight_vacuum_final", pre);
        T.add("solve.iterations", sr.iterations);
        chain.text(p.candidate_field_sha256()); chain.text(p.candidate_heat_sha256());
        if (step == 1 || step % every == 0 || step == steps)
            ledger << J{{"step", step}, {"time_s", p.end_s()}, {"field_sha256", p.candidate_field_sha256()}, {"heat_sha256", p.candidate_heat_sha256()},
                        {"current_A", r.current_A}, {"terminal_V", r.terminal_V}, {"joule_J", r.joule_J}, {"field_delta_J", r.field_delta_J},
                        {"balance_J", r.balance_J}, {"max_gauss_C", c.max_gauss_C}, {"max_divB_T_m", c.max_divB_T_m}, {"residual_ratio", sr.independent_residual / sr.target},
                        {"iterations", sr.iterations}, {"max_T_K", res.thermal().max_temperature_K()}, {"control", control_json(p.control())}}.dump() << '\n';
        t = steady::now();
        res.publish_vacuum_prevalidated(std::move(p));
        T.add("step.publish", ms_since(t));
        T.add("step.total", prep + pre + T.m["step.publish"].v.back());
    }
    const double wall = ms_since(t_all);
    J summary{{"schema", "hall-v3.1-vacuum-run-v1"}, {"steps", steps}, {"wall_ms", wall}, {"ms_per_step", wall / steps}, {"start_time_s", seed.time_s},
              {"end_time_s", res.time_s()}, {"verify_all", verify_all}, {"final_field_sha256", res.field_sha256()}, {"final_heat_sha256", res.heat_sha256()},
              {"chain_sha256", chain.hex()}, {"max_T_K", res.thermal().max_temperature_K()},
              {"guards", {{"max_gauss_over_bound", max_gauss_ratio}, {"max_divB_over_bound", max_divB_ratio}, {"max_energy_over_bound", max_energy_ratio},
                          {"max_independent_residual_over_target", max_resid_ratio}, {"max_abs_field_heat_balance_J", max_heat_balance}}},
              {"digests", {{"computed", res.device_field().digest_stats().digests}, {"leaves_hashed_on_host", res.device_field().digest_stats().leaves_hashed_on_host},
                           {"samples", res.device_field().digest_stats().samples}}},
              {"phases", T.json()}, {"device", J::parse(S.dev->device_json())}, {"memory", J::parse(S.dev->memory_json())}};
    {
        axial::State fs; hc::ChargeLedger fq; res.accepted_state(fs, fq);
        summary["final_physical_sha256"] = state_values_sha256(fs, fq, res.thermal());
    }
    summary["teardown"] = J::parse(S.dev->close_json());
    std::ofstream(out / "vacuum-run.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return 0;
}

int coupled_run(const std::string& trace, const std::filesystem::path& out, long steps, std::size_t carriers, const std::string& nvrtc, const std::string& uuid,
                bool gpu_motion) {
    hc::CarrierLimits lim{std::max<std::size_t>(100000, carriers), 1ull << 34};
    auto S = setup(trace, nvrtc, uuid, lim, false);
    const auto& g = S.cpu->mesh(); const auto o = mf::geometry(g); const auto seed = S.cpu->dc_seed();
    auto census = carriers == 2 ? mf::pair({g.side / 2, g.side / 2, 0.}, 1.) : mf::pairs(carriers, 7, 1e3);
    gp::ResidentOptions opt;
    hc::ResidentOwner res(S.resident, S.dev, seed, declared(*S.cpu, seed, hc::charge_of_census(g, o, census, lim.maximum_carriers)), o, census, opt);
    std::unique_ptr<hc::MotionProvider> provider;
    if (gpu_motion) provider = std::make_unique<gp::GpuMotionProvider>(S.dev, res.motion_context(), carriers);
    else provider = std::make_unique<hc::CpuMotionProvider>();
    Stats T; hc::Digest chain("hall-v3.1-coupled-run-chain-v1");
    std::filesystem::create_directories(out);
    std::ofstream ledger(out / "coupled-ledger.jsonl");
    const std::string src = hc::sha256("resident/source");
    const auto t_all = steady::now();
    for (long step = 1; step <= steps; ++step) {
        const auto h0 = hc::motion_host_timers();
        auto t = steady::now();
        auto b = binding_for(res.time_s(), res.revision(), res.motion_mesh_sha256(), S.resident->material_sha256(), res.motion_context().policy_sha256, step);
        auto p = hc::prepare_closed_proposal(res, o, census, b, src, *provider, hc::StepInputs{});
        const double prep = ms_since(t); t = steady::now();
        const auto& bb = p.motion().binding();
        Boundary bd{p.candidate_field_sha256(), p.candidate_heat_sha256(), src, bb.mesh_sha256, bb.material_sha256, bb.common_clock_sha256, p.end_s(), p.endpoint_revision()};
        const bool ok = res.preflight_final(p, census, bd);
        const double pre = ms_since(t);
        hc::need(ok, "coupled preflight");
        const auto& TT = p.timings();
        T.add("step.prepare_closed_proposal", prep); T.add("step.preflight_final", pre);
        T.add("proposal.admit_census", 1e3 * TT.admit_census_s); T.add("proposal.motion", 1e3 * TT.motion_s); T.add("proposal.support", 1e3 * TT.support_s);
        T.add("proposal.field_prepare(device)", 1e3 * TT.field_prepare_s); T.add("proposal.motion_v2_seal", 1e3 * TT.motion_v2_seal_s); T.add("proposal.heat+digests", 1e3 * TT.heat_s);
        T.add("proposal.iterations", TT.iterations);
        const auto h1 = hc::motion_host_timers();
        T.add("motion_host.census_hash", 1e3 * (h1.census_hash_s - h0.census_hash_s)); T.add("motion_host.v3_digest", 1e3 * (h1.record_v3_s - h0.record_v3_s));
        T.add("motion_host.v2_seal", 1e3 * (h1.record_seal_s - h0.record_seal_s));
        {   // probe: deposit support and bank transfer of the accepted motion, timed in isolation
            auto t2 = steady::now(); (void)hc::DepositSupportBuilder::build(g, o, census, p.motion()); T.add("probe.deposit_support", ms_since(t2));
            t2 = steady::now(); (void)hc::DepositSupportBuilder::build(g, o, census, p.motion()); T.add("probe.deposit_support_repeat", ms_since(t2));
            t2 = steady::now(); (void)hc::DepositSupportBuilder::build(g, o, census, p.motion()); T.add("probe.deposit_support_repeat2", ms_since(t2));
            t2 = steady::now(); (void)hc::bank_transfer(g, o, p.motion(), S.resident->exterior_support()); T.add("probe.bank_transfer", ms_since(t2));
        }
        chain.text(p.candidate_field_sha256()); chain.text(p.candidate_heat_sha256()); chain.text(p.motion().record_v3_sha256());
        if (step == 1 || step % 100 == 0 || step == steps)
            ledger << J{{"step", step}, {"time_s", p.end_s()}, {"field_sha256", p.candidate_field_sha256()}, {"motion_v3", p.motion().record_v3_sha256()},
                        {"iterations", p.iterations()}, {"combined_balance_J", p.combined_balance_J()}, {"control", control_json(p.field_heat().control())},
                        {"gyro_max", p.motion().budget().maximum_gyro}, {"displacement_cells_max", p.motion().budget().maximum_displacement_cells},
                        {"plasma_dt", p.motion().budget().plasma_dt}}.dump() << '\n';
        auto next = hc::next_census(census, p.motion());
        t = steady::now();
        res.publish_prevalidated(std::move(p));
        T.add("step.publish", ms_since(t));
        T.add("step.total", prep + pre + T.m["step.publish"].v.back());
        census = std::move(next);
    }
    const double wall = ms_since(t_all);
    J summary{{"schema", "hall-v3.1-coupled-run-v1"}, {"steps", steps}, {"carriers", carriers}, {"motion_provider", provider->name()}, {"wall_ms", wall},
              {"ms_per_step", wall / steps}, {"chain_sha256", chain.hex()}, {"final_field_sha256", res.field_sha256()}, {"phases", T.json()},
              {"device", J::parse(S.dev->device_json())}};
    provider.reset();
    summary["teardown"] = J::parse(S.dev->close_json());
    std::ofstream(out / ("coupled-run-" + std::to_string(carriers) + ".json")) << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return 0;
}
// Absolute Gauss defect trace of the vacuum step (diagnostic): per-node |defect| and bound of the accepted
// state every EVERY steps, with the accumulated per-node bound of the solves' independent residuals
// sum_k sum_{e at n} sqrt(M_e) |rho_k,e| / len_e (rho = b - A x, the vector the independent check measures).
int gauss_trace(const std::string& trace, const std::filesystem::path& out, long steps, long every, const std::string& nvrtc, const std::string& uuid, double dt) {
    auto S = setup(trace, nvrtc, uuid, {}, false, std::max(dt, 1e-14));   // dt > 1e-14: diagnostic configuration only
    const auto& g = S.cpu->mesh(); const auto o = mf::geometry(g); const auto seed = S.cpu->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    gp::ResidentOptions opt; opt.trace_residual_charge = true;
    hc::ResidentOwner res(S.resident, S.dev, seed, declared(*S.cpu, seed, hc::charge_of_census(g, o, none)), o, none, opt);
    std::filesystem::create_directories(out);
    std::ofstream ledger(out / "gauss-trace.jsonl");
    hc::Vec d0, b0, d, b;
    res.device_field().gauss_profile(0, res.accepted_meta(), d0, b0);
    double sum_rhs = 0, max_ratio = 0, max_corrected = 0, max_explained = 0, max_gx = 0, max_gc = 0; long refused_at = -1; std::string refusal;
    const auto t0 = steady::now();
    for (long step = 1; step <= steps; ++step) {
        try {
            auto p = res.prepare_vacuum(dt);
            Boundary bd{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
            hc::need(res.preflight_vacuum_final(p, bd), "vacuum preflight");
            const auto sr = p.solve_record(); const auto r = p.field_receipt(); const auto st = p.step();
            sum_rhs += sr.rhs_norm;
            res.publish_vacuum_prevalidated(std::move(p));
            if (step % every == 0 || step == 1 || step == steps) {
                res.device_field().gauss_profile(0, res.accepted_meta(), d, b);
                const auto& cum = res.device_field().residual_charge_cumulative();
                hc::Vec gacc; res.device_field().gauss_increment_cumulative(gacc);
                double mr = 0, mc = 0, me = 0, md = 0, mb = 0, mcum = 0, drift = 0, mg = 0, mgx = 0, mgc = 0; int node = -1;
                for (int k = 0; k < g.nn; ++k) {
                    const double ratio = d[k] / b[k];
                    if (ratio > mr) { mr = ratio; node = k; }
                    mc = std::max(mc, d[k] / (b[k] + cum[k]));
                    if (cum[k] > 0) me = std::max(me, (d[k] - d0[k]) / cum[k]);
                    md = std::max(md, d[k]); mb = std::max(mb, b[k]); mcum = std::max(mcum, cum[k]); drift = std::max(drift, std::abs(d[k] - d0[k]));
                    mg = std::max(mg, gacc[k]); if (gacc[k] > 0) mgx = std::max(mgx, std::abs(d[k] - d0[k]) / gacc[k]); mgc = std::max(mgc, d[k] / (b[k] + gacc[k]));
                }
                max_ratio = std::max(max_ratio, mr); max_corrected = std::max(max_corrected, mc); max_explained = std::max(max_explained, me);
                ledger << J{{"step", step}, {"time_s", res.time_s()}, {"terminal_V", r.terminal_V}, {"current_A", r.current_A}, {"rhs_norm", sr.rhs_norm},
                            {"sum_rhs_norm", sum_rhs}, {"independent_residual", sr.independent_residual}, {"target", sr.target},
                            {"max_defect_C", md}, {"max_bound_C", mb}, {"max_node_ratio", mr}, {"worst_node", node},
                            {"max_drift_C", drift}, {"max_cumulative_residual_charge_C", mcum}, {"max_drift_over_cumulative", me},
                            {"max_ratio_with_cumulative_allowance", mc}, {"step_residual_charge_max_C", st.residual_charge_step_max_C},
                            {"increment_rel", st.control.increment_rel}, {"max_T_K", res.thermal().max_temperature_K()},
                            {"max_cumulative_gauss_increment_C", mg}, {"max_drift_over_cumulative_gauss_increment", mgx},
                            {"max_ratio_with_cumulative_gauss_increment_allowance", mgc}, {"worst_node_defect_C", node >= 0 ? d[node] : 0.},
                            {"worst_node_bound_C", node >= 0 ? b[node] : 0.}, {"worst_node_cumulative_gauss_increment_C", node >= 0 ? gacc[node] : 0.}}.dump() << '\n';
                max_gx = std::max(max_gx, mgx); max_gc = std::max(max_gc, mgc);
                ledger.flush();
            }
        } catch (const std::exception& e) { refused_at = step; refusal = e.what(); break; }
    }
    J summary{{"schema", "hall-v3.1-gauss-trace-v1"}, {"diagnostic_only", true}, {"dt_s", dt}, {"steps_requested", steps}, {"refused_at_step", refused_at}, {"refusal", refusal},
              {"end_time_s", res.time_s()}, {"wall_s", ms_since(t0) / 1e3}, {"max_node_ratio", max_ratio}, {"max_ratio_with_cumulative_allowance", max_corrected},
              {"max_drift_over_cumulative_residual_charge", max_explained}, {"max_drift_over_cumulative_gauss_increment", max_gx},
              {"max_ratio_with_cumulative_gauss_increment_allowance", max_gc}, {"initial_max_defect_C", *std::max_element(d0.begin(), d0.end())},
              {"initial_max_bound_C", *std::max_element(b0.begin(), b0.end())}, {"initial_gauss_solve_bound_C", S.resident->initial_gauss_solve_bound_C()}};
    std::ofstream(out / "gauss-trace.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    summary["teardown"] = J::parse(S.dev->close_json());
    return 0;
}
// Diagnostic step-size study of the particle-free (vacuum) field step: the same window from the
// supported start is run with several dt on a Model whose dt upper bound is raised for this study
// only (non-default configuration identity, not a qualification of a larger step). Differences of
// the final state are measured against the finest step.
int dt_study(const std::string& trace, const std::filesystem::path& out, const std::string& nvrtc, const std::string& uuid, double window_s, std::vector<double> dts) {
    auto dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    const double cap = *std::max_element(dts.begin(), dts.end());
    auto probe = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto desc = hc::resident_descriptor_sha256(dev, *probe);
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<hc::ResidentBackend>(desc), hc::CarrierLimits{}, hc::ThermalDeclaration{},
                                                   std::max(cap, 1e-14));
    const auto& g = model->mesh(); const auto o = mf::geometry(g); const auto seed = model->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    const auto q0 = declared(*model, seed, hc::charge_of_census(g, o, none));
    struct Run { double dt; long steps; axial::State s; double joule_J = 0, thermal_J = 0, maxT = 0, ms_per_step = 0; long iterations = 0; std::vector<double> V; double max_third = 0, max_incr = 0, max_omega_c = 0; };
    std::vector<Run> runs;
    const double sample_every_s = 1e-11;
    for (double dt : dts) {
        hc::ResidentOwner res(model, dev, seed, q0, o, none, gp::ResidentOptions{});
        Run r; r.dt = dt; r.steps = std::lround(window_s / dt);
        const long every = std::max(1L, std::lround(sample_every_s / dt));
        const auto t0 = steady::now();
        for (long k = 1; k <= r.steps; ++k) {
            auto p = res.prepare_vacuum(dt);
            Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
            hc::need(res.preflight_vacuum_final(p, b), "dt study preflight");
            r.joule_J += p.field_receipt().joule_J; r.iterations += p.solve_record().iterations;
            r.max_third = std::max(r.max_third, p.control().third_difference_rel); r.max_incr = std::max(r.max_incr, p.control().increment_rel);
            r.max_omega_c = std::max(r.max_omega_c, p.control().omega_c_dt_max);
            if (k % every == 0) r.V.push_back(p.field_receipt().terminal_V);
            res.publish_vacuum_prevalidated(std::move(p));
        }
        r.ms_per_step = ms_since(t0) / r.steps;
        hc::ChargeLedger q; res.accepted_state(r.s, q);
        r.thermal_J = res.thermal().energy_J(); r.maxT = res.thermal().max_temperature_K();
        runs.push_back(std::move(r));
        std::cout << "dt " << dt << " steps " << runs.back().steps << " ms/step " << runs.back().ms_per_step << std::endl;
    }
    const auto& ref = runs.front();   // finest step first
    auto norms = [&](const axial::State& a, const axial::State& b, double& dE, double& nE, double& dB, double& nB) {
        long double sE = 0, sEn = 0, sB = 0, sBn = 0;
        for (int e = 0; e < g.ne; ++e) { const long double x = (long double)a.E[e] + a.E_low[e], y = (long double)b.E[e] + b.E_low[e]; sE += g.M[e] * (x - y) * (x - y); sEn += g.M[e] * y * y; }
        for (int f = 0; f < g.nb; ++f) { const long double x = (long double)a.B[f] + a.B_low[f], y = (long double)b.B[f] + b.B_low[f]; sB += g.H[f] * (x - y) * (x - y); sBn += g.H[f] * y * y; }
        dE = std::sqrt(double(sE)); nE = std::sqrt(double(sEn)); dB = std::sqrt(double(sB)); nB = std::sqrt(double(sBn));
    };
    J rows = J::array();
    for (const auto& r : runs) {
        double dE, nE, dB, nB; norms(r.s, ref.s, dE, nE, dB, nB);
        double dV = 0, nV = 0;
        const std::size_t m = std::min(r.V.size(), ref.V.size());
        for (std::size_t i = 0; i < m; ++i) { dV = std::max(dV, std::abs(r.V[i] - ref.V[i])); nV = std::max(nV, std::abs(ref.V[i])); }
        rows.push_back({{"dt_s", r.dt}, {"steps", r.steps}, {"end_time_s", r.s.time_s}, {"ms_per_step", r.ms_per_step}, {"mean_cg_iterations", double(r.iterations) / r.steps},
                        {"E_energy_norm_rel_diff", dE / nE}, {"B_energy_norm_rel_diff", dB / nB}, {"terminal_V_max_abs_diff", dV}, {"terminal_V_max_abs", nV},
                        {"joule_J", r.joule_J}, {"joule_rel_diff", rel(r.joule_J, ref.joule_J)}, {"thermal_energy_J", r.thermal_J},
                        {"thermal_rel_diff", rel(r.thermal_J, ref.thermal_J)}, {"max_T_K", r.maxT}, {"max_third_difference_rel", r.max_third},
                        {"max_increment_rel", r.max_incr}, {"max_omega_c_dt", r.max_omega_c}});
    }
    J summary{{"schema", "hall-v3.1-dt-study-v1"}, {"diagnostic_only", true},
              {"note", "non-default dt upper bound used only in this study; the default and qualified bound stays 1e-14 s"},
              {"window_s", window_s}, {"start_s", seed.time_s}, {"configuration_sha256", model->configuration_sha256()}, {"rows", rows},
              {"sigma_initial_S_m", *std::max_element(model->initial_conductivity().begin(), model->initial_conductivity().end())},
              {"mesh", {{"dx_m", g.dx}, {"dy_m", g.dy}, {"dz_m", g.dz}, {"nx", g.nx}, {"ny", g.ny}, {"nz", g.nz}, {"side_m", g.side}}},
              {"device", J::parse(dev->device_json())}};
    std::filesystem::create_directories(out);
    std::ofstream(out / "dt-study.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    dev->close_json();
    return 0;
}
}  // namespace

int main(int argc, char** argv) {
    // Large host work arrays come from the heap, not from fresh mmap/munmap pairs per step (with the
    // device driver's threads in the process, every munmap is a cross-core TLB shootdown).
    mallopt(M_MMAP_THRESHOLD, 256 << 20); mallopt(M_TRIM_THRESHOLD, 512 << 20);
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "units" && argc == 6) return units(argv[2], argv[3], argv[4], argv[5]);
        if (mode == "parity-vacuum" && argc == 7) return parity_vacuum(argv[2], argv[3], std::stol(argv[4]), argv[5], argv[6]);
        if (mode == "parity-coupled" && argc == 8) return parity_coupled(argv[2], argv[3], std::stol(argv[4]), std::stoul(argv[5]), argv[6], argv[7]);
        if (mode == "negative" && argc == 6) return negative(argv[2], argv[3], argv[4], argv[5]);
        if (mode == "vacuum-run" && argc >= 7)
            return vacuum_run(argv[2], argv[3], std::stol(argv[4]), argv[5], argv[6], argc > 7 && std::string(argv[7]) == "verify-all", argc > 8 ? std::stol(argv[8]) : 100);
        if (mode == "coupled-run" && argc >= 8)
            return coupled_run(argv[2], argv[3], std::stol(argv[4]), std::stoul(argv[5]), argv[6], argv[7], !(argc > 8 && std::string(argv[8]) == "cpu"));
        if (mode == "gauss-trace" && argc >= 8) return gauss_trace(argv[2], argv[3], std::stol(argv[4]), std::stol(argv[5]), argv[6], argv[7], argc > 8 ? std::stod(argv[8]) : 1e-14);
        if (mode == "dt-study" && argc >= 8) {
            std::vector<double> dts; for (int i = 7; i < argc; ++i) dts.push_back(std::stod(argv[i]));
            return dt_study(argv[2], argv[3], argv[4], argv[5], std::stod(argv[6]), dts);
        }
        std::cerr << "usage: see the header of test_resident.cpp\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
}
