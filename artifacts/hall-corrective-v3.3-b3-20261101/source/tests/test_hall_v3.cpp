// Hall corrective v3 component qualification (CPU backends).
//   test_hall_v3 components TRACE OUTDIR     sinks, banks, Gauss/ledger/energy closure, partial-path current,
//                                            node planes, grazing, ties, births, loads, typed refusals
//   test_hall_v3 parity-dump TRACE OUTFILE   v2/v3 physical-state parity scenarios (values only)
//   test_hall_v3 motion-dump TRACE OUTFILE   v2 fixture set record digests (compared with the v2 dumps)
//   test_hall_v3 sink-dump TRACE OUTFILE     sink fixture records/arrays (compared with the GPU provider)
//   test_hall_v3 multistep TRACE OUTDIR STEPS BIRTH_EVERY BIRTHS   emitted electrons crossing the gap
#include "root_api.hpp"
#include "hall-finite-foil-em-20261029/axial_particles.hpp"
#include "sink_core.hpp"
namespace hc = hall_corrective_v3;
#include "parity_scenarios.hpp"
#include "sink_fixtures.hpp"
#include "sink_dump.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <optional>
#include <bit>
#include <map>

namespace hd = hall_driver;
namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
namespace sf = sink_fixtures;
namespace eh = electrode_heat_v2;
using J = hc::Json;
using steady = std::chrono::steady_clock;

namespace hall_corrective_v3 {
struct ModelTestAccess {
    static FieldCandidate prepare(const Model& m, const axial::State& old, const ChargeLedger& before, const Vec& plasma,
                                  const FreeCharge& after, const DepositSupport& support, double dt, const BankTransfer* bank) {
        const auto sigma = m.initial_conductivity();
        return m.prepare(old, before, m.state_sha256(old, before), plasma, after, support, sigma, Model::conductivity_sha256(sigma), dt, bank);
    }
};
}  // namespace hall_corrective_v3

namespace {
constexpr double dt = 1e-14;
double since(steady::time_point t) { return std::chrono::duration<double>(steady::now() - t).count(); }
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
template <class F> std::string refusal(F&& f) {
    try { f(); } catch (const std::exception& e) { return e.what(); }
    return "";
}
using parity::declared;
hd::MotionBinding binding_for(const hc::Owner& o, const std::string& tag) {
    hd::MotionBinding b; b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = o.revision(); b.staged_arena_revision = o.revision() + 1; b.begin_s = o.time_s(); b.dt_s = dt; b.end_s = b.begin_s + dt;
    b.staged_arena_sha256 = hc::sha256(tag + "/staged"); b.mesh_sha256 = o.motion_mesh_sha256(); b.material_sha256 = o.model().material_sha256();
    b.species_sha256 = hc::sha256(tag + "/species"); b.common_clock_sha256 = hc::sha256(tag + "/clock"); b.kinetic_policy_sha256 = o.motion_context().policy_sha256;
    return b;
}
parity::Boundary boundary_of(const hc::JointProposal& p, const std::string& source) {
    const auto& b = p.motion().binding();
    return {p.candidate_field_sha256(), p.candidate_heat_sha256(), source, b.mesh_sha256, b.material_sha256, b.common_clock_sha256, p.end_s(), p.endpoint_revision()};
}
long double total_charge(const hc::ChargeLedger& q) {
    long double s = 0;
    for (std::size_t k = 0; k < q.particles.net_C.size(); ++k) s += q.particles.net_C[k].wide() + q.conductor_C[k].wide() + q.exterior_C[k].wide();
    return s;
}
long double charge_magnitude(const hc::ChargeLedger& q) {
    long double s = 0;
    for (std::size_t k = 0; k < q.particles.net_C.size(); ++k) s += q.particles.gross_C[k] + hc::magnitude(q.conductor_C[k]) + hc::magnitude(q.exterior_C[k]);
    return s;
}
J energy_json(const hc::StepEnergy& e) {
    return {{"field_delta_J", e.field_delta_J}, {"plasma_work_J", e.plasma_work_J}, {"joule_J", e.joule_J}, {"outer_loss_J", e.outer_loss_J},
            {"outer_supply_J", e.outer_supply_J}, {"circuit_work_J", e.circuit_work_J}, {"kinetic_change_J", e.kinetic_change_J},
            {"plasma_kinetic_store_change_J", e.plasma_kinetic_store_change_J},
            {"thermal_delta_J", e.thermal_delta_J}, {"wall_heat_J", e.wall_heat_J}, {"external_J", e.external_J},
            {"collected_kinetic_J", e.collected_kinetic_J}, {"short_bridge_kinetic_J", e.short_bridge_kinetic_J},
            {"escaped_kinetic_J", e.escaped_kinetic_J}, {"birth_kinetic_J", e.birth_kinetic_J}, {"birth_supplied_J", e.birth_supplied_J},
            {"combined_balance_J", e.combined_balance_J}, {"residual_J", e.residual_J}, {"scale_J", e.scale_J}};
}
bool energy_closes(const hc::StepEnergy& e, double roundoff) { return std::isfinite(e.residual_J) && std::abs(e.residual_J) <= 4e-20 + 2e-10 * e.scale_J + roundoff; }

// Independent re-derivation of the native v6 admission rules for the Hall records
// (API-006.txt sections 2 and 4): first hit of (start, segment) with the sealed Hall
// helper (native's port is bit-identical), hit point, flight F f, discarded F - flight,
// path row, neutral straight drift, kind/surface, and the flows record digest with a
// test-side encoder of full_time_motion_record_v3_sha256.
J native_conformance(const hall_driver::MotionReceipt& m, std::span<const n::ArenaRateEntry> census) {
    J ev = J::object(); bool ok = true; long checked = 0;
    std::map<std::uint64_t, const n::ArenaRateEntry*> byid; for (const auto& r : census) byid[r.native.id] = &r;
    std::map<std::uint64_t, double> flight; for (const auto& b : m.births()) flight[b.native_id] = b.flight_s;
    std::map<std::uint64_t, const hall_driver::CarrierPath*> path; for (const auto& p : m.paths()) path[p.native_id] = &p;
    const auto& pol = m.kinetic_policy(); const auto shift = m.frame().native_to_field_m;
    for (const auto& k : m.sinks()) {
        const auto& old = byid.at(k.native_id)->native; const double F = flight.count(k.native_id) ? flight[k.native_id] : m.binding().dt_s;
        bool r = (k.kind == 1) == (k.boundary == 0 || k.boundary == 1 || k.boundary == 3) && (k.kind == 1 || k.kind == 2) && k.axis >= 0 && k.axis <= 2;
        r = r && k.fraction >= 0 && k.fraction <= 1 && k.flight_s == F * k.fraction && k.flight_s <= F && k.discarded_flight_s == F - k.flight_s;
        hall_driver::KineticHit hit;
        if (flight.count(k.native_id)) { hgp::KineticHit h; const hgp::Policy P{pol.width_m, pol.gap_m, pol.short_z_m, pol.feed_z_m,
            {pol.lower.left_m, pol.lower.right_m, pol.lower.bottom_m, pol.lower.top_m}, {pol.upper.left_m, pol.upper.right_m, pol.upper.bottom_m, pol.upper.top_m}};
            (void)hgp::first_kinetic_hit_core(P, {{old.position_m[0], old.position_m[1], old.position_m[2]}}, {{k.segment_m[0], k.segment_m[1], k.segment_m[2]}}, h);
            hit.sink = hall_driver::KineticSink(h.sink); hit.fraction = h.fraction; hit.axis = h.axis; hit.plane_m = h.plane; }
        else hit = hall_driver::first_kinetic_hit(pol, old.position_m, k.segment_m);
        r = r && int(hit.sink) == k.boundary && hit.fraction == k.fraction && hit.axis == k.axis && hit.plane_m == k.plane_m;
        n::PicVector point{}; for (int d = 0; d < 3; ++d) point[d] = old.position_m[d] + k.fraction * k.segment_m[d]; point[k.axis] = k.plane_m;
        r = r && point == k.hit_point_m;
        if (old.charge_C == 0) { double u2 = 0; for (double v : old.proper_velocity_m_s) u2 += v * v; const double gamma = std::sqrt(1 + u2 / (n::pic_c * n::pic_c));
            for (int d = 0; d < 3; ++d) r = r && k.segment_m[d] == F * old.proper_velocity_m_s[d] / gamma;
            r = r && k.incident_proper_velocity_m_s == old.proper_velocity_m_s; }
        const auto* pr = path.at(k.native_id);
        r = r && pr->flight_s == k.flight_s;
        for (int d = 0; d < 3; ++d) r = r && pr->field_displacement_m[d] == (k.hit_point_m[d] + shift[d]) - pr->field_start_m[d];
        n::GapParticle inc = old; inc.proper_velocity_m_s = k.incident_proper_velocity_m_s;
        const double energy = old.weight * n::pic_kinetic_energy_J(inc);
        r = r && std::isfinite(energy) && energy >= 0 && energy == k.weighted_kinetic_J;
        ok = ok && r; ++checked;
    }
    // test-side flows digest encoder (native full_time_motion_record_v3_sha256, flows schema)
    if (!m.v2_identical()) {
        hc::Digest unused("unused"); (void)unused;
        std::string bytes;
        auto u64 = [&](std::uint64_t v) { for (int i = 0; i < 8; ++i) bytes.push_back(char((v >> (8 * i)) & 255)); };
        auto f64 = [&](double v) { u64(std::bit_cast<std::uint64_t>(v)); };
        auto text = [&](const std::string& t) { u64(t.size()); bytes += t; };
        auto arr = [&](const std::vector<double>& v) { std::string keep; keep.swap(bytes); text("motion-array-binary64-v1"); u64(v.size()); for (double x : v) f64(x);
            const auto h = hc::bytes_sha256(bytes); bytes.swap(keep); return h; };
        const auto& b = m.binding();
        const auto cur = arr(m.current_A_m()), endc = arr(m.endpoint_charge_C()), ini = arr(m.initial_charge_C());
        text("canonical-full-time-total-field-motion-v3-flows-binary"); text("canonical-full-time-total-field-motion-v2"); text("binding");
        text("full_time_position_and_proper_velocity"); u64(b.accepted_step); u64(b.staged_arena_revision); u64(b.trial_stage); f64(b.begin_s); f64(b.end_s); f64(b.dt_s);
        for (const auto* x : {&b.staged_arena_sha256, &b.trial_field_sha256, &b.mesh_sha256, &b.material_sha256, &b.species_sha256, &b.common_clock_sha256, &b.kinetic_policy_sha256}) text(*x);
        const auto& u = m.budget(); text("budget");
        for (double x : {u.kinetic_change_J, u.field_work_J, u.work_defect_J, u.continuity_C, u.plasma_dt, u.maximum_gyro, u.maximum_displacement_cells}) f64(x);
        text("census_sha256"); text(m.census_sha256()); text("coverage"); text("all_native_carriers_births_and_first_hit_sinks");
        text("current_A_m_sha256"); text(cur); text("endpoint_charge_C_sha256"); text(endc);
        text("endpoints"); u64(m.endpoints().size());
        for (const auto& e : m.endpoints()) { u64(e.id); for (double v : e.position_m) f64(v); for (double v : e.proper_velocity_m_s) f64(v); f64(e.charge_C); f64(e.mass_kg); f64(e.weight); }
        text("initial_charge_C_sha256"); text(ini);
        text("kinetic_policy"); text("native-finite-target-kinetic-policy-v1");
        for (double x : {pol.width_m, pol.gap_m, pol.short_z_m, pol.feed_z_m, pol.lower.left_m, pol.lower.right_m, pol.lower.bottom_m, pol.lower.top_m,
                         pol.upper.left_m, pol.upper.right_m, pol.upper.bottom_m, pol.upper.top_m}) f64(x);
        for (const char* x : {"lower foil x span and inner y gap only", "absorbing escape", "axial escape", "x y lower upper axial; later exact tie wins"}) text(x);
        text("native_DKD_called"); u64(0); text("native_Poisson_called"); u64(0); text("native_to_field_m"); for (double v : shift) f64(v);
        text("neutral_motion_owner"); text("same_full_time_q0_drift");
        text("paths"); u64(m.paths().size());
        for (const auto& p : m.paths()) { u64(p.native_id); for (double v : p.field_start_m) f64(v); for (double v : p.field_displacement_m) f64(v); f64(p.flight_s); }
        text("ordinary_total_field");
        text("sinks"); u64(m.sinks().size());
        for (const auto& k : m.sinks()) { u64(k.native_id); u64(std::uint64_t(k.kind)); u64(std::uint64_t(k.boundary)); u64(std::uint64_t(k.axis)); f64(k.plane_m); f64(k.fraction);
            for (double v : k.segment_m) f64(v); for (double v : k.hit_point_m) f64(v); for (double v : k.incident_proper_velocity_m_s) f64(v); f64(k.flight_s); f64(k.discarded_flight_s); }
        text("births"); u64(m.births().size());
        for (const auto& nb : m.births()) { u64(nb.native_id); f64(nb.birth_time_s); f64(nb.flight_s); }
        const auto digest = hc::bytes_sha256(bytes);
        ev["flows_digest_independent"] = digest; ev["flows_digest_receipt"] = m.record_v3_sha256();
        ok = ok && digest == m.record_v3_sha256();
    }
    ev["sink_records_checked"] = checked; ev["ok"] = ok;
    return ev;
}
// One closed step with checks of every closure; publishes when preflight passes.
struct StepCheck {
    J ev = J::object(); bool ok = false; std::vector<n::ArenaRateEntry> next;
};
StepCheck closed_step(hc::Owner& owner, const n::GapPicOptions& o, const std::vector<n::ArenaRateEntry>& census, const std::string& tag,
                      const hc::StepInputs& in = {}, hc::MotionProvider* provider = nullptr) {
    StepCheck s;
    const auto before = owner.charge(); const long double q0 = total_charge(before);
    const auto& g = owner.mesh();
    hc::CpuMotionProvider cpu;
    auto p = hc::prepare_closed_proposal(owner, o, census, binding_for(owner, tag), hc::sha256(tag + "/source"), provider ? *provider : cpu, in);
    const auto& fh = p.field_heat(); const auto& c = fh.candidate_constraints(); const auto& q1 = fh.candidate_charge();
    const auto bd = boundary_of(p, hc::sha256(tag + "/source"));
    // Bank changes of this step = candidate banks - before banks - metal/outer deltas: verified through absolute Gauss
    // (admission) and the total-charge identity: free + conductor + exterior is conserved by every transfer.
    const long double q1t = total_charge(q1);
    long double free0 = 0, free1 = 0;
    for (int k = 0; k < g.nn; ++k) { free0 += before.particles.net_C[k].wide(); free1 += q1.particles.net_C[k].wide(); }
    long double sink_q = 0, birth_q = 0;
    for (const auto& k : p.sinks()) sink_q += (long double)(hc::Words(k.weight) * k.charge_C).wide();
    for (const auto& b : p.births()) birth_q += (long double)(hc::Words(b.weight) * b.charge_C).wide();
    // native v6 sink-record conformance (API-006.txt section 2) on every record of the step
    s.ev["native_v6_conformance"] = native_conformance(p.motion(), census);
    long double gross = 0; for (double x : before.particles.gross_C) gross += x;
    const auto bank = hc::bank_transfer(g, o, p.motion(), owner.model().exterior_support());
    long double bank_q = 0; if (!bank.empty) for (int k = 0; k < g.nn; ++k) bank_q += bank.conductor_C[k].wide() + bank.exterior_C[k].wide();
    s.ev["sinks"] = p.sinks().size(); s.ev["births"] = p.births().size(); s.ev["iterations"] = p.iterations();
    s.ev["max_gauss_C"] = c.max_gauss_C; s.ev["max_gauss_bound_C"] = c.max_gauss_bound_C; s.ev["max_divB_T_m"] = c.max_divB_T_m;
    s.ev["continuity_ratio"] = p.timings().continuity_ratio;
    s.ev["total_charge_before_C"] = double(q0); s.ev["total_charge_change_C"] = double(q1t - q0);
    s.ev["free_charge_change_C"] = double(free1 - free0); s.ev["sink_macro_charge_C"] = double(sink_q); s.ev["birth_macro_charge_C"] = double(birth_q);
    s.ev["bank_transfer_total_C"] = double(bank_q);
    s.ev["energy"] = energy_json(p.energy()); s.ev["combined_balance_J"] = p.combined_balance_J();
    s.ev["heat"] = {{"thermal_delta_J", fh.heat().thermal_delta_J}, {"joule_J", fh.heat().admitted_joule_J}, {"wall_heat_J", fh.heat().wall_heat_J},
                    {"external_J", fh.heat().external_J}, {"loads", fh.heat().loads}, {"roundoff_bound_J", fh.heat().roundoff_bound_J}};
    J sinks = J::array();
    for (const auto& k : p.sinks()) sinks.push_back({{"id", k.native_id}, {"boundary", hc::sink_boundary_names[k.boundary]}, {"axis", k.axis},
        {"fraction", k.fraction}, {"flight_s", k.flight_s}, {"discarded_flight_s", k.discarded_flight_s},
        {"hit_point_m", {k.hit_point_m[0], k.hit_point_m[1], k.hit_point_m[2]}}, {"plane_m", k.plane_m}, {"weighted_kinetic_J", k.weighted_kinetic_J}});
    s.ev["sink_records"] = sinks;
    const double qscale = double(gross) + std::abs(double(sink_q)) + std::abs(double(birth_q)) + 1e-30;
    // Total charge (free + conductor + exterior) is conserved by every transfer of the step up to rounding of the banks.
    const double qabs = double(charge_magnitude(before) + charge_magnitude(q1));
    s.ev["charge_magnitude_C"] = qabs;
    const bool charge_ok = std::abs(double(q1t - q0)) <= 1e-12 * qabs;
    // free charge lost = sink macro-charge - birth macro-charge; bank gained the same.
    const bool partition_ok = std::abs(double((free0 - free1) - (sink_q - birth_q))) <= 1e-12 * qscale &&
                              std::abs(double(bank_q - (sink_q - birth_q))) <= 1e-12 * qscale;
    const bool energy_ok = energy_closes(p.energy(), fh.heat().roundoff_bound_J);
    const bool pre = owner.preflight_final(p, census, bd);
    s.ev["charge_conserved"] = charge_ok; s.ev["free_bank_partition"] = partition_ok; s.ev["energy_closed"] = energy_ok; s.ev["preflight"] = pre;
    s.next = hc::next_census(census, p.motion());
    if (pre) owner.publish_prevalidated(std::move(p));
    s.ok = charge_ok && partition_ok && energy_ok && pre && c.max_gauss_C <= c.max_gauss_bound_C && s.ev["native_v6_conformance"].value("ok", false);
    return s;
}

// ------------------------------------------------------------------ components
int components(const std::string& trace, const std::filesystem::path& out) {
    Results R;
    J report = J::object();
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto& g = model->mesh();
    const auto o = mf::geometry(g);
    const auto seed = model->dc_seed();
    report["configuration_sha256"] = model->configuration_sha256();
    report["material_sha256"] = model->material_sha256();
    const auto ctx = hc::make_motion_context(g, o, model->mesh_json());
    report["sink_support"] = ctx.sink_support_text;
    R.record("support", "kinetic boundary planes are aligned with field node planes and every sink kind has an admitted bank support",
             ctx.sink_support_aligned && ctx.sink_supported[0] && ctx.sink_supported[1] && ctx.sink_supported[2] && ctx.sink_supported[3] && ctx.sink_supported[4],
             {{"support", ctx.sink_support_text}});

    // ---- 1. every fixture: expected classification, closures, publication
    for (const auto& f : sf::fixture_set()) {
        if (f.field != 0) continue;   // E x 100 fixtures are motion-only (sink-dump / GPU parity)
        J ev{{"fixture", f.name}};
        bool ok = false;
        try {
            std::vector<hc::Newborn> newborns;
            for (const auto& [id, frac] : f.newborns) newborns.push_back({id, model->drive().start() + frac * dt});
            const auto residents = hc::residents(f.census, newborns);
            hc::Owner owner(model, seed, declared(*model, seed, hc::charge_of_census(g, o, residents)), o, residents);
            auto s = closed_step(owner, o, f.census, f.name, hc::StepInputs{newborns, {}});
            ev["step"] = s.ev;
            bool classes = true;
            std::vector<int> got;
            for (const auto& c : f.census) {
                int kind = -1;
                for (const auto& k : s.ev["sink_records"]) if (k["id"].get<std::uint64_t>() == c.native.id) {
                    const auto name = k["boundary"].get<std::string>(); for (int b = 0; b < 5; ++b) if (name == hc::sink_boundary_names[b]) kind = b; }
                got.push_back(kind);
            }
            ev["classification"] = got;
            if (!f.expected_sinks.empty()) classes = got == f.expected_sinks;
            ok = s.ok && classes;
            // a second step from the published state (survivors only) also closes
            if (ok && !s.next.empty()) { auto s2 = closed_step(owner, o, s.next, f.name + "/2"); ev["second_step"] = s2.ev; ok = s2.ok; }
        } catch (const std::exception& e) { ev["refusal"] = e.what(); ok = false; }
        R.record("sinks", f.name + ": classification, absolute Gauss, charge-bank partition, energy closure incl. wall heat, preflight/publish", ok, ev);
    }

    // ---- 2. partial-path current: the deposit covers only the path start -> hit point
    {
        const auto census = sf::fixture_set()[0].census;   // lower-foil-top-face
        const auto mid = model->at(seed);
        auto b = parity::binding(*std::make_unique<hc::Owner>(model, seed, declared(*model, seed, hc::charge_of_census(g, o, census)), o, census), 1);
        b.trial_field_sha256 = hd::total_field_sha256(mid.E(), mid.B());
        auto m = hc::MotionAssembly::cpu(ctx, census, b, {}, mid.E(), mid.B());
        J ev;
        const auto& k = m.sinks()[0]; const auto& path = m.paths()[0]; const auto& partner = m.paths()[1];
        hc::Vec J1(g.ne); n::PicVector D{};
        for (int d = 0; d < 3; ++d) D[d] = k.field_point_m[d] - path.field_start_m[d];
        axial::Path(g, path.field_start_m, path.field_displacement_m).deposit(J1, census[0].native.weight * census[0].native.charge_C, path.field_displacement_m, dt);
        axial::Path(g, partner.field_start_m, partner.field_displacement_m).deposit(J1, census[1].native.weight * census[1].native.charge_C, partner.field_displacement_m, dt);
        hc::Vec Jfull(g.ne); n::PicVector Df{};
        {const auto& u = census[0].native.proper_velocity_m_s; const double gam = std::sqrt(1 + (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) / (n::pic_c * n::pic_c));
         for (int d = 0; d < 3; ++d) Df[d] = dt * u[d] / gam;   // the untruncated full-interval straight segment
         axial::Path(g, path.field_start_m, Df).deposit(Jfull, census[0].native.weight * census[0].native.charge_C, Df, dt);
         axial::Path(g, partner.field_start_m, partner.field_displacement_m).deposit(Jfull, census[1].native.weight * census[1].native.charge_C, partner.field_displacement_m, dt);}
        double diff_full = 0; for (int e = 0; e < g.ne; ++e) diff_full = std::max(diff_full, std::abs(Jfull[e] - m.current_A_m()[e]));
        ev["receipt_current_equals_truncated_path_deposit"] = J1 == m.current_A_m();
        ev["path_displacement_equals_hit_point_minus_start"] = D == path.field_displacement_m;
        ev["path_flight_equals_sink_flight"] = path.flight_s == k.flight_s && k.flight_s < dt && k.discarded_flight_s == dt - k.flight_s && k.flight_s == dt * k.fraction;
        ev["max_difference_from_untruncated_deposit_A_m"] = diff_full;
        ev["hit_y_native_on_plane"] = k.hit_point_m[1] == k.plane_m || std::abs(k.hit_point_m[1] - k.plane_m) <= 1e-18;
        ev["budget"] = {{"continuity_C", m.budget().continuity_C}, {"work_defect_J", m.budget().work_defect_J}, {"kinetic_change_J", m.budget().kinetic_change_J}};
        R.record("partial-path", "the current deposit is exactly the deposit of the path from the start to the exact hit point (flight < dt)",
                 ev["receipt_current_equals_truncated_path_deposit"].get<bool>() && ev["path_displacement_equals_hit_point_minus_start"].get<bool>() &&
                     ev["path_flight_equals_sink_flight"].get<bool>() && diff_full > 0, ev);
    }

    // ---- 3. negative controls: the bank transfer and the free-charge removal are both required
    {
        const auto census = sf::fixture_set()[0].census;
        const auto before = declared(*model, seed, hc::charge_of_census(g, o, census));
        hc::Owner owner(model, seed, before, o, census);
        const auto mid = model->at(seed);
        auto b = binding_for(owner, "negative"); b.trial_field_sha256 = hd::total_field_sha256(mid.E(), mid.B());
        auto m = hc::MotionAssembly::cpu(owner.motion_context(), census, b, {}, mid.E(), mid.B());
        const auto after = hc::charge_of_endpoints(g, o, m.endpoints());
        const auto bank = hc::bank_transfer(g, o, m, model->exterior_support());
        const auto support = hc::DepositSupportBuilder::build(g, o, census, m);
        J ev;
        ev["control"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, seed, before, m.current_A_m(), after, support, dt, &bank); });
        ev["bank_transfer_omitted"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, seed, before, m.current_A_m(), after, support, dt, nullptr); });
        auto wrong = bank; for (auto& w : wrong.conductor_C) w = w * .5;
        ev["half_bank_transfer"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, seed, before, m.current_A_m(), after, support, dt, &wrong); });
        auto moved = bank; std::swap(moved.conductor_C, moved.exterior_C);
        ev["collected_charge_booked_as_exterior"] = refusal([&] { (void)hc::ModelTestAccess::prepare(*model, seed, before, m.current_A_m(), after, support, dt, &moved); });
        const std::string text = "actual particle endpoint charge/current mismatch";
        R.record("negative", "omitted or halved bank transfer refuses (endpoint guard); collected charge booked on the exterior bank refuses (support admission)",
                 ev["control"] == "" && ev["bank_transfer_omitted"] == text && ev["half_bank_transfer"] == text &&
                     ev["collected_charge_booked_as_exterior"].get<std::string>().find("exterior charge outside") != std::string::npos, ev);
    }

    // ---- 4. newborns (native v6 births over [tau, end]) on the lower-foil face and a collected electron in one step
    {
        auto census = sf::with_partners({sf::electron(1, {.0045, .007 - 3e-8, .001}, 100, {0, 1, 0}, 1e6)});
        hc::Owner owner(model, seed, declared(*model, seed, hc::charge_of_census(g, o, census)), o, census);
        std::vector<hc::Newborn> newborns;
        auto with_births = census;
        for (int k = 0; k < 3; ++k) {
            auto e = sf::electron(10 + k, {.003 + .001 * k, sf::lower_top, -.002 + .002 * k}, 1e6, {0, 1, 0}, 1e6);   // on the lower-foil top face, leaving it
            with_births.push_back(e); newborns.push_back({e.native.id, owner.time_s() + (.25 + .25 * k) * dt});
        }
        {   auto h = mf::row(20, "h2o-newborn", "H2O", {.0045, sf::lower_top, .003}, {0, 650, 0}, 0., mf::mn, 1e9);   // neutral packet on the face (source v5 desorption)
            with_births.push_back(h); newborns.push_back({20, owner.time_s()}); }
        StepCheck s; J ev;
        try { s = closed_step(owner, o, with_births, "newborns", hc::StepInputs{newborns, {}}); ev = s.ev; } catch (const std::exception& e) { ev["refusal"] = e.what(); }
        const bool ok = s.ok && s.ev.value("births", 0) == 4 && s.ev.value("sinks", 0) == 1 && s.next.size() == 5;
        ev["next_census"] = s.next.size();
        bool second = false;
        if (ok) { auto s2 = closed_step(owner, o, s.next, "newborns/2"); ev["second_step"] = s2.ev; second = s2.ok; }
        ev["ledgers"] = {{"sinks", owner.kinetic_ledgers().sinks.json()}, {"births", owner.kinetic_ledgers().births.json()}};
        R.record("newborns", "three electrons and one H2O packet born on the lower-foil face over [tau, end] and one electron collected at the upper "
                             "foil: conductor -> free -> conductor charge, flows record (native v6 rules), closures, second step", ok && second, ev);
        J bad;
        auto try_births = [&](std::vector<n::ArenaRateEntry> c, std::uint64_t id, double frac) {
            return refusal([&] { hc::Owner x(model, seed, declared(*model, seed, hc::charge_of_census(g, o, census)), o, census); hc::CpuMotionProvider cpu;
                std::vector<hc::Newborn> nb{{id, x.time_s() + frac * dt}};
                (void)hc::prepare_closed_proposal(x, o, c, binding_for(x, "bad"), hc::sha256("bad"), cpu, hc::StepInputs{nb, {}}); }); };
        {auto c = census; c.push_back(sf::electron(10, {.003, .0019, 0.}, 1e6, {0, 1, 0}, 1e6)); bad["start_inside_metal"] = try_births(c, 10, .5);}
        {auto c = census; c.push_back(sf::electron(10, {.003, sf::lower_top, 0.}, 1e6, {0, -1, 0}, 1e6)); bad["velocity_into_metal"] = try_births(c, 10, .5);}
        {auto c = census; c.push_back(sf::electron(10, {sf::lower_left, sf::lower_top, 0.}, 1e6, {0, 1, 0}, 1e6)); bad["start_on_foil_corner"] = try_births(c, 10, .5);}
        {auto c = census; c.push_back(sf::electron(10, {.003, .003, 0.}, 1e6, {0, 1, 0}, 1e6)); bad["charged_interior_newborn"] = try_births(c, 10, .5);}
        {auto c = census; c.push_back(sf::electron(10, {.003, sf::lower_top, 0.}, 1e6, {0, 1, 0}, 1e6)); bad["birth_time_at_end"] = try_births(c, 10, 1.);}
        {bad["unknown_newborn_id"] = try_births(census, 99, .5);}
        const std::string face = "newborn start neither a strictly interior vacuum point nor a foil-face point leaving the metal";
        bool all = bad["start_inside_metal"] == face && bad["velocity_into_metal"] == face && bad["start_on_foil_corner"] == face &&
                   bad["charged_interior_newborn"].get<std::string>().find("charged interior newborn") != std::string::npos &&
                   bad["birth_time_at_end"].get<std::string>().find("birth time outside the step") != std::string::npos &&
                   bad["unknown_newborn_id"] == "newborn native ID absent from the census";
        R.record("newborns", "newborn refusals: start inside metal, velocity into the metal, foil corner, charged interior newborn (no source path), "
                             "no remaining flight, unknown ID", all, bad);
    }

    // ---- 5. external surface loads: exactly-once accounting, double-Joule detection, typed refusals
    {
        auto load = [&](const std::string& producer, std::uint64_t seq, double energy, int body, int k, const std::string& channel, double t) {
            hc::SurfaceLoad l; l.body = body; l.slice = k; l.channel = channel; auto& x = l.transfer;
            x.producer = producer; x.energy_owner = producer + "/energy"; x.sequence = x.owner_sequence = seq; x.state_id = "electrode";
            x.parent_id = "electrode"; x.channel_id = channel; x.event_time_s = t; x.surface = true; x.x = .001; x.y = body ? .0025 : -.0025;
            x.authored_heat_J = energy; return l; };
        J ev;
        hc::Owner owner(model, seed, declared(*model, seed, {hc::WordVec(g.nn), hc::Vec(g.nn)}), o, {});
        const double t = owner.time_s() + .5 * dt;
        std::vector<hc::SurfaceLoad> L{load("test-emission", 1, -2e-9, 0, g.kshort + 1, "emission-cooling", t), load("test-desorption", 1, -5e-10, 1, g.kshort, "desorption-binding", t),
                                       load("test-ion-impact", 1, 3e-9, 1, g.kshort + 2, "ion-impact", t)};
        auto f = owner.prepare_vacuum(dt, L);
        const auto& h = f.heat();
        ev["thermal_delta_J"] = h.thermal_delta_J; ev["joule_J"] = h.admitted_joule_J; ev["external_J"] = h.external_J; ev["roundoff_bound_J"] = h.roundoff_bound_J;
        ev["closure_residual_J"] = h.thermal_delta_J - h.admitted_joule_J - h.external_J;
        ev["field_heat_balance_J"] = f.field_heat_balance_J();
        const bool closes = std::abs(h.thermal_delta_J - h.admitted_joule_J - h.external_J) <= h.roundoff_bound_J + 2e-12 * (h.external_heating_J + h.external_cooling_J);
        parity::Boundary bd{f.candidate_field_sha256(), f.candidate_heat_sha256(), "", "", "", "", f.end_s(), f.endpoint_revision()};
        const bool pre = owner.preflight_vacuum_final(f, bd);
        owner.publish_vacuum_prevalidated(std::move(f));
        ev["ledger_after"] = {{"loads", owner.kinetic_ledgers().loads}, {"external_J", owner.kinetic_ledgers().external_J.rounded()}};
        const double t2 = owner.time_s() + .5 * dt;
        ev["replayed_sequence_next_step"] = refusal([&] { std::vector<hc::SurfaceLoad> x{load("test-emission", 1, -2e-9, 0, g.kshort + 1, "emission-cooling", t2)}; (void)owner.prepare_vacuum(dt, x); });
        ev["duplicate_in_one_step"] = refusal([&] { std::vector<hc::SurfaceLoad> x{load("test-emission", 2, -1e-9, 0, g.kshort + 1, "emission-cooling", t2),
                                                                                    load("test-emission", 2, -1e-9, 0, g.kshort + 1, "emission-cooling", t2)}; (void)owner.prepare_vacuum(dt, x); });
        ev["energy_owner_reused_by_another_producer"] = refusal([&] { auto a = load("test-a", 5, 1e-10, 0, g.kshort, "ion-impact", t2); auto b2 = load("test-b", 5, 1e-10, 0, g.kshort, "ion-impact", t2);
            b2.transfer.energy_owner = a.transfer.energy_owner; std::vector<hc::SurfaceLoad> x{a, b2}; (void)owner.prepare_vacuum(dt, x); });
        ev["external_joule_channel"] = refusal([&] { std::vector<hc::SurfaceLoad> x{load("source-private-heat", 1, 1e-9, 0, g.kshort, "joule", t2)}; (void)owner.prepare_vacuum(dt, x); });
        ev["external_joule_owner"] = refusal([&] { auto l = load("source-v5", 1, 1e-9, 0, g.kshort, "private-heat", t2); l.transfer.energy_owner = "source-v5/Joule";
            std::vector<hc::SurfaceLoad> x{l}; (void)owner.prepare_vacuum(dt, x); });
        ev["reserved_hall_producer"] = refusal([&] { auto l = load("hall-v3/kinetic-sink-wall-heat", 99999999999, 1e-9, 0, g.kshort, "ion-impact", t2); std::vector<hc::SurfaceLoad> x{l}; (void)owner.prepare_vacuum(dt, x); });
        ev["accepted_load_once_more_after_refusals"] = refusal([&] { std::vector<hc::SurfaceLoad> x{load("test-emission", 3, -1e-9, 0, g.kshort + 1, "emission-cooling", t2)}; (void)owner.prepare_vacuum(dt, x); });
        const std::string wm = "electrode heat: duplicate/stale event or canonical energy owner; missing species/parent/channel identity";
        R.record("loads", "external loads: thermal delta = Joule + loads (exactly once), ledger updated; replayed sequence, duplicate and shared energy owner refused",
                 closes && pre && ev["replayed_sequence_next_step"].get<std::string>().find("duplicate/stale") != std::string::npos &&
                     ev["duplicate_in_one_step"].get<std::string>().find("duplicate/stale") != std::string::npos &&
                     ev["energy_owner_reused_by_another_producer"].get<std::string>().find("duplicate/stale") != std::string::npos &&
                     ev["accepted_load_once_more_after_refusals"] == "", ev);
        R.record("double-joule", "an external Joule load (channel or energy owner) and a Hall-reserved producer are refused typed: Joule is counted once",
                 ev["external_joule_channel"].get<std::string>().find("[external-joule]") != std::string::npos &&
                     ev["external_joule_owner"].get<std::string>().find("[external-joule]") != std::string::npos &&
                     ev["reserved_hall_producer"].get<std::string>().find("[reserved-producer]") != std::string::npos, ev);
        J tr;
        const double t3 = owner.time_s() + .5 * dt;
        auto typed = [&](std::vector<hc::SurfaceLoad> x) -> J { try { (void)owner.prepare_vacuum(dt, x); return {{"admitted", true}}; }
            catch (const hc::ThermalRefusal& r) { return {{"kind", r.kind}, {"quantity", r.quantity}, {"identity", r.identity}, {"observed", r.observed}, {"bound", r.bound}, {"what", r.what()}}; }
            catch (const std::exception& e) { return {{"untyped", e.what()}}; } };
        const auto cap = owner.thermal().capacity(owner.thermal().index(0, g.kshort));
        tr["slice_capacity"] = {{"heating_room_J", cap.heating_room_J}, {"cooling_room_J", cap.cooling_room_J}, {"full_range_J", cap.full_range_J}};
        tr["stock_1e200_scale_load"] = typed({load("test-huge", 10, -1e181, 0, g.kshort, "desorption-binding", t3)});
        tr["release_5e21_scale_load"] = typed({load("test-5e21", 10, -5e21 * 1.3e-19, 0, g.kshort, "desorption-binding", t3)});
        tr["net_heating_above_room"] = typed({load("test-h1", 10, .6 * cap.full_range_J, 0, g.kshort, "ion-impact", t3), load("test-h2", 10, .6 * cap.full_range_J, 0, g.kshort, "ion-impact", t3)});
        tr["nonfinite_load"] = typed({load("test-nan", 10, std::nan(""), 0, g.kshort, "ion-impact", t3)});
        tr["concentrated_load_within_capacity"] = typed({load("test-point", 10, .9 * cap.heating_room_J, 0, g.kshort, "ion-impact", t3)});
        tr["owner_unchanged"] = owner.revision() == 1;
        R.record("typed-refusals", "heat-engine extremes are typed scoped refusals before the engine (1e181 J, 5e21-molecule-scale load, net heating above room, "
                                   "nonfinite) or typed engine refusals (concentrated load); no abort, owner unchanged",
                 tr["stock_1e200_scale_load"].value("kind", "") == "load-exceeds-material-range" && tr["release_5e21_scale_load"].value("kind", "") == "load-exceeds-material-range" &&
                     tr["net_heating_above_room"].value("kind", "") == "net-heating-exceeds-capacity" && tr["nonfinite_load"].value("kind", "") == "nonfinite-load" &&
                     (tr["concentrated_load_within_capacity"].value("kind", "") == "thermal-engine" || tr["concentrated_load_within_capacity"].value("admitted", false)) &&
                     tr["owner_unchanged"].get<bool>(), tr);
    }

    // ---- 6. heat-only source patch declaration leaves the thermal evolution bit-identical
    {
        hc::ThermalDeclaration decl;
        eh::Patch p; p.id = "lower-face"; p.face = eh::Face::HighY; p.source_hash = "source-v6-scenario";
        for (int k = g.kshort; k < g.kfeed; ++k) decl.slices[{0, k}] = {"source-v6-scenario", {p}};
        auto declared_model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, hc::make_cpu_pcg(), hc::CarrierLimits{}, decl);
        hc::Owner a(model, seed, declared(*model, seed, {hc::WordVec(g.nn), hc::Vec(g.nn)}), o, {});
        hc::Owner b(declared_model, seed, declared(*declared_model, seed, {hc::WordVec(g.nn), hc::Vec(g.nn)}), o, {});
        for (int s = 0; s < 2; ++s) {
            auto fa = a.prepare_vacuum(dt); auto fb = b.prepare_vacuum(dt);
            parity::Boundary ba{fa.candidate_field_sha256(), fa.candidate_heat_sha256(), "", "", "", "", fa.end_s(), fa.endpoint_revision()};
            parity::Boundary bb{fb.candidate_field_sha256(), fb.candidate_heat_sha256(), "", "", "", "", fb.end_s(), fb.endpoint_revision()};
            (void)a.preflight_vacuum_final(fa, ba); (void)b.preflight_vacuum_final(fb, bb);
            a.publish_vacuum_prevalidated(std::move(fa)); b.publish_vacuum_prevalidated(std::move(fb)); }
        bool same = true; const auto sa = a.thermal().states(), sb = b.thermal().states();
        for (std::size_t i = 0; i < sa.size(); ++i) same = same && sa[i].T == sb[i].T && sa[i].cumulative == sb[i].cumulative;
        R.record("declaration", "declared heat-only source patches leave temperatures bit-identical; the material/configuration identities change",
                 same && declared_model->material_sha256() != model->material_sha256() && declared_model->configuration_sha256() != model->configuration_sha256(),
                 {{"patch_surface_T", b.thermal().slice(0, g.kshort).surface_temperatures()}, {"material_plain", model->material_sha256()},
                  {"material_declared", declared_model->material_sha256()}});
    }

    // ---- 6b. empty birth scope (arrival committed every step, no births this step)
    {
        const auto census = mf::fixture_set(g)[1].census;   // pair-z0-w1: no hits
        hc::Owner a(model, seed, declared(*model, seed, hc::charge_of_census(g, o, census)), o, census);
        hc::Owner b(model, seed, declared(*model, seed, hc::charge_of_census(g, o, census)), o, census);
        hc::CpuMotionProvider cpu;
        auto pa = hc::prepare_closed_proposal(a, o, census, binding_for(a, "scope"), hc::sha256("scope/source"), cpu, hc::StepInputs{});
        hc::StepInputs in; in.birth_scope = true;
        auto pb = hc::prepare_closed_proposal(b, o, census, binding_for(b, "scope"), hc::sha256("scope/source"), cpu, in);
        const auto conf = native_conformance(pb.motion(), census);
        const bool bank_empty = hc::bank_transfer(g, o, pb.motion(), model->exterior_support()).empty;
        const bool pre = b.preflight_final(pb, census, boundary_of(pb, hc::sha256("scope/source")));
        J ev{{"plain_record_v3", pa.motion().record_v3_sha256()}, {"flows_record_v3", pb.motion().record_v3_sha256()}, {"conformance", conf},
             {"bank_transfer_empty", bank_empty}, {"preflight", pre}};
        R.record("birth-scope", "an empty birth scope gives the flows record with empty sink/birth sections (native v6 digest), no bank transfer, "
                                "and the candidate field/heat of the plain step bit for bit",
                 pa.motion().v2_identical() && !pb.motion().v2_identical() && pb.motion().birth_scope() && pb.motion().sinks().empty() &&
                     pb.motion().births().empty() && pa.candidate_field_sha256() == pb.candidate_field_sha256() &&
                     pa.candidate_heat_sha256() == pb.candidate_heat_sha256() && pa.motion().record_v3_sha256() != pb.motion().record_v3_sha256() &&
                     conf.value("ok", false) && conf.contains("flows_digest_independent") && bank_empty && pre,
                 ev);
    }
    // ---- 7. checkpoint round trip keeps the kinetic ledgers
    {
        const auto census = sf::fixture_set()[0].census;
        hc::Owner owner(model, seed, declared(*model, seed, hc::charge_of_census(g, o, census)), o, census);
        auto s = closed_step(owner, o, census, "checkpoint");
        const auto dir = out / "checkpoint-after-sink"; std::filesystem::remove_all(dir); owner.save(dir);
        std::string why; std::unique_ptr<hc::Owner> back;
        try { back = hc::Owner::restore(model, o, s.next, dir); } catch (const std::exception& e) { why = e.what(); }
        R.record("checkpoint", "save/restore after a collection step restores the banks and the kinetic ledgers",
                 back && back->kinetic_ledgers_sha256() == owner.kinetic_ledgers_sha256() && back->field_sha256() == owner.field_sha256() &&
                     owner.kinetic_ledgers().sinks.count() == 1, {{"refusal", why}, {"ledgers", owner.kinetic_ledgers().sinks.json()}});
    }
    report["cases"] = R.cases;
    report["failed"] = R.failed;
    std::ofstream(out / "components.json") << report.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// ------------------------------------------------------------------ dumps
int motion_dump(const std::string& trace, const std::string& file) {
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto& g = model->mesh(); const auto o = mf::geometry(g);
    const auto ctx = hc::make_motion_context(g, o, model->mesh_json());
    const auto mid = model->at(model->dc_seed());
    hc::Vec E100 = mid.E(); for (double& v : E100) v *= 100;
    J out = J::array();
    for (const auto& f : mf::fixture_set(g)) {
        const hc::Vec& E = f.field ? E100 : mid.E();
        hd::MotionBinding b; b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
        b.accepted_step = 0; b.staged_arena_revision = 1; b.begin_s = model->drive().start(); b.dt_s = dt; b.end_s = b.begin_s + dt;
        b.staged_arena_sha256 = hc::sha256("dump/staged"); b.trial_field_sha256 = hd::total_field_sha256(E, mid.B());
        b.mesh_sha256 = ctx.mesh_sha256; b.material_sha256 = model->material_sha256(); b.species_sha256 = hc::sha256("dump/species");
        b.common_clock_sha256 = hc::sha256("dump/clock"); b.kinetic_policy_sha256 = ctx.policy_sha256;
        for (int api = 0; api < 2; ++api) {
            J r{{"fixture", f.name}, {"carriers", f.census.size()}, {"entry", api ? "prepare_canonical_motion" : "MotionAssembly::cpu"}};
            const auto t = steady::now();
            try {
                auto m = api ? hd::prepare_canonical_motion(g, o, f.census, b, {}, E, mid.B()) : hc::MotionAssembly::cpu(ctx, f.census, b, {}, E, mid.B());
                hc::MotionAssembly::seal_v2(m);
                J rec = hd::motion_record(m); const auto text = rec.dump(); hc::drain(rec);
                r["record_sha256"] = hc::sha256(text); r["seal_sha256"] = m.seal_sha256(); r["record_v3_sha256"] = m.record_v3_sha256();
                try { hc::Digest d("test-after-charge"); auto q = hc::charge_of_endpoints(g, o, m.endpoints()); hc::encode(d, q.net_C); d.f64s(q.gross_C); r["after_sha256"] = d.hex(); }
                catch (const std::exception& e) { r["after_refusal"] = e.what(); }
            } catch (const std::exception& e) { r["refusal"] = e.what(); }
            r["wall_s"] = since(t);
            out.push_back(r);
        }
    }
    std::ofstream(file) << out.dump(1) << '\n';
    return 0;
}

// ------------------------------------------------------------------ multistep: emitted electrons cross the gap
int multistep(const std::string& trace, const std::filesystem::path& out, long steps, long birth_every, long births_total) {
    std::filesystem::create_directories(out);
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto& g = model->mesh(); const auto o = mf::geometry(g); const auto seed = model->dc_seed();
    std::vector<n::ArenaRateEntry> census;
    hc::Owner owner(model, seed, declared(*model, seed, {hc::WordVec(g.nn), hc::Vec(g.nn)}), o, census);
    std::ofstream ledger(out / "step-ledger.jsonl");
    std::mt19937_64 rng(20261030); std::uniform_real_distribution<double> ux(.0025, .0065), uz(-.003, .003), ue(0, 1);
    std::uint64_t next_id = 1; long born = 0; std::string stop = "steps"; long done = 0;
    J summary{{"schema", "hall-v3-multistep-summary-v1"}, {"scope", "Hall v3 closed loop with newborn cathode electrons (native v6 birth model: "
                                                                     "start on the lower-foil face, flight end - tau, charge from the conductor bank; injected, "
                                                                     "not sampled by the source owner) and first-hit collection; no native arena or chemistry"},
              {"configuration_sha256", model->configuration_sha256()}, {"steps_requested", steps}, {"birth_every", birth_every}, {"births_requested", births_total}};
    long double max_gauss_ratio = 0, max_energy_rel = 0; double t_total = 0;
    for (long step = 1; step <= steps; ++step) {
        std::vector<hc::Newborn> births;
        if (born < births_total && (step - 1) % birth_every == 0) {
            const double x = ux(rng), z = uz(rng), spread = .15 * (2 * ue(rng) - 1);
            auto e = sf::electron(1000000 + (next_id++), {x, sf::lower_top, z}, 2e6, {spread, 1, .1 * spread}, 1e6);   // on the cathode face, leaving it
            census.push_back(e); births.push_back({e.native.id, owner.time_s() + .999 * ue(rng) * dt}); ++born;
        }
        J row{{"step", step}};
        const auto t0 = steady::now();
        try {
            const auto before_q = total_charge(owner.charge());
            auto s = closed_step(owner, o, census, "multistep/" + std::to_string(step), hc::StepInputs{births, {}});
            row["closures"] = s.ev; row["ok"] = s.ok;
            const auto& L = owner.kinetic_ledgers();
            row["census_after"] = s.next.size(); row["total_charge_C"] = double(total_charge(owner.charge())); row["total_charge_change_C"] = double(total_charge(owner.charge()) - before_q);
            row["cumulative"] = {{"collected", L.sinks.boundary[0].count + L.sinks.boundary[1].count + L.sinks.boundary[3].count},
                                 {"collected_upper_foil", L.sinks.boundary[1].count}, {"collected_lower_foil", L.sinks.boundary[0].count},
                                 {"escaped", L.sinks.boundary[2].count + L.sinks.boundary[4].count}, {"births", L.births.count()},
                                 {"collected_charge_C", L.sinks.collected_charge_C().rounded()}, {"escaped_charge_C", L.sinks.escaped_charge_C().rounded()},
                                 {"wall_heat_J", L.wall_heat_J.rounded()}, {"escaped_J", L.escaped_J.rounded()}, {"birth_kinetic_J", L.births.kinetic_J().rounded()}};
            max_gauss_ratio = std::max<long double>(max_gauss_ratio, s.ev["max_gauss_C"].get<double>() / std::max(1e-300, s.ev["max_gauss_bound_C"].get<double>()));
            const auto& e = s.ev["energy"];
            max_energy_rel = std::max<long double>(max_energy_rel, std::abs(e["residual_J"].get<double>()) / std::max(1e-300, e["scale_J"].get<double>()));
            census = std::move(s.next);
            if (!s.ok) { stop = "closure check failed at step " + std::to_string(step); ledger << row.dump() << '\n'; break; }
        } catch (const std::exception& e) { row["refusal"] = e.what(); stop = std::string("refused at step ") + std::to_string(step) + ": " + e.what(); ledger << row.dump() << '\n'; break; }
        row["step_ms"] = 1e3 * since(t0); t_total += since(t0);
        ledger << row.dump() << '\n'; ++done;
        if (step % 100 == 0) { ledger.flush(); std::cout << "step " << step << " census " << census.size() << " collected " << row["cumulative"]["collected"] << std::endl; }
        if (born == births_total && census.empty()) { stop = "all born electrons left the domain"; break; }
    }
    const auto& L = owner.kinetic_ledgers();
    summary["steps_completed"] = done; summary["stop"] = stop; summary["births"] = born; summary["census_final"] = census.size();
    summary["sink_ledger"] = L.sinks.json(); summary["birth_ledger"] = L.births.json();
    summary["wall_heat_J"] = L.wall_heat_J.rounded(); summary["escaped_J"] = L.escaped_J.rounded(); summary["short_bridge_J"] = L.short_bridge_J.rounded();
    summary["max_gauss_over_bound"] = double(max_gauss_ratio); summary["max_energy_residual_over_scale"] = double(max_energy_rel);
    summary["mean_step_ms"] = done ? 1e3 * t_total / done : 0; summary["field_sha256"] = owner.field_sha256(); summary["kinetic_ledgers_sha256"] = owner.kinetic_ledgers_sha256();
    summary["max_temperature_K"] = owner.thermal().max_temperature_K();
    std::ofstream(out / "multistep-summary.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return stop.rfind("refused", 0) == 0 || stop.rfind("closure", 0) == 0 ? 1 : 0;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 4) throw std::runtime_error("usage: test_hall_v3 components|parity-dump|motion-dump|sink-dump|multistep TRACE OUT [...]");
        const std::string mode = argv[1];
        if (mode == "components") { std::filesystem::create_directories(argv[3]); return components(argv[2], argv[3]); }
        if (mode == "parity-dump") { std::ofstream(argv[3]) << parity::run(argv[2]).dump(1) << '\n'; return 0; }
        if (mode == "motion-dump") return motion_dump(argv[2], argv[3]);
        if (mode == "sink-dump") { auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[2]); hc::CpuMotionProvider cpu;
            std::ofstream(argv[3]) << sink_dump::run(model, cpu).dump(1) << '\n'; return 0; }
        if (mode == "multistep" && argc == 7) return multistep(argv[2], argv[3], std::stol(argv[4]), std::stol(argv[5]), std::stol(argv[6]));
        throw std::runtime_error("unknown mode");
    } catch (const std::exception& e) { std::cerr << "fatal: " << e.what() << '\n'; return 2; }
}
