// Native v6 admission of Hall corrective v3 motion (genuine Hall v3 CPU provider,
// identically zero total field, the native v6 storage-control fixture on the Hall
// geometry). Every step is opened, sealed and published by the native v6 arena
// (native-full-time-candidate006, compiled here from its sealed sources); the motion
// is the Hall v3 receipt bound through hall_v3_flows_adapter.hpp. Native recomputes
// the record digest, the first hit of every sink, the hit point, the truncated
// path, the newborn flights and the kinetic budget, and refuses any inconsistency.
//   1. plain v3 record of a v2-identical Hall v3 receipt (no newborns, no hits);
//   2. first-hit sinks of every surface (native sink population): admission, native
//      sink records equal the Hall records field by field, native per-surface sink
//      totals equal the Hall SinkLedger exactly (count, weight, kinetic energy, +/- e);
//   3. births on the upper-foil face moving into the gap, sink-prone births and
//      interior births (arrival v1 plans through the native birth scope);
//   4. births and sinks in the same step;
//   5. tampered Hall records are refused typed by native;
//   6. multi-step run: electrons born on the lower-foil face each step for the first
//      steps cross the 5 mm gap and are collected at the upper foil; every step's
//      ledgers compared (native vs Hall);
//   7. arrival policy v1.1 (quantized neutral release) through the native v6 birth
//      scope: neutral leaves of exactly w0 born on the upper-foil face, native stream
//      debit = n w0 exactly, Hall motion admitted every step.
// usage: test_native_v6_admission OUTPUT_DIRECTORY [STEPS]
#include "v6_fixture.hpp"
#include "arrival_v1_births.hpp"
#include "pf_arrival_v11/native_v6_births.hpp"
#include "hall_v3_flows_adapter.hpp"
#include "sink_ledger.hpp"
#include "hall-finite-foil-em-20261029/axial_mesh.hpp"
#include <bit>
#include <functional>

using namespace v6;
namespace pa = pf_arrival;
namespace av = m::arrival_v1;
namespace hc = hall_corrective_v3;
namespace hn = hall_v3_native_v6;
namespace pa11 = pf_arrival_v11;
namespace av11 = m::arrival_v11;
using FK = m::FullTimeFlowRefusal::Kind;

namespace {
J measurements = J::object();
const axial::Mesh& mesh() { static const axial::Mesh g; return g; }
const std::string& mesh_json() {
    static const std::string s = [] { m::ArenaJson j = mesh().json(); auto t = j.dump(); hc::drain(j); return t; }();
    return s;
}
const hc::MotionContext& context() { static const auto c = hc::make_motion_context(mesh(), options(), mesh_json()); return c; }
m::FullTimeGenesisBinding hall_genesis() { auto g = genesis(); g.mesh_sha256 = hall_driver::motion_mesh_sha256(mesh(), options()); return g; }
std::unique_ptr<m::PicArena> fresh(std::vector<m::FullTimeAuthoredPacket> pop = base_population()) {
    return m::PicArena::create_full_time(options(), states(), channels(), std::move(pop), hall_genesis(), {});
}
// Copy of the native v6 sink population (test_v6_flows.cpp): carriers 2e-7 m from
// each kinetic surface moving into it, a co-located two-carrier packet, one survivor.
std::vector<m::FullTimeAuthoredPacket> sink_population() {
    auto pop = base_population();
    const double d = 2e-7;
    pop.push_back({"to-lower", "e", "manufactured/sink", {3, 0}, {.004, .002 + d, .001}, {0, -1e8, 0}});
    pop.push_back({"to-upper", "i", "manufactured/sink", {5, 0}, {.005, .007 - 1e-12, .002}, {0, 6e3, 0}});
    pop.push_back({"to-outer-x", "e", "manufactured/sink", {7, 0}, {.009 - d, .0045, .003}, {1e8, 0, 0}});
    pop.push_back({"to-outer-y", "n", "manufactured/sink", {1.5, 0}, {.0010, 1e-12, .004}, {0, -650, 0}});
    pop.push_back({"to-bridge", "e", "manufactured/sink", {2, 0}, {.0045, .0047, -.0125 + d}, {0, 0, -1e8}});
    pop.push_back({"to-feed", "e", "manufactured/sink", {6, 0}, {.0045, .0043, .0125 - d}, {0, 3e5, 1e8}});
    pop.push_back({"to-axial-short", "e", "manufactured/sink", {9, 0}, {.0010, .0045, -.0125 + d}, {0, 0, -1e8}});
    pop.push_back({"pair-to-x0", "e", "manufactured/sink", {1.0, std::ldexp(1.0, -60)}, {d, .0050, .005}, {-1e8, 0, 0}});
    pop.push_back({"stays", "e", "manufactured/sink", {4, 0}, {.0060, .0045, .006}, {1e6, 0, 0}});
    return pop;
}
bool same(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }
bool same(const m::PicVector& a, const m::PicVector& b) { return same(a[0], b[0]) && same(a[1], b[1]) && same(a[2], b[2]); }

// ---- the Hall v3 motion hook
struct HallStep {
    std::vector<hc::SinkRecord> sinks;
    std::vector<hc::NewbornRecord> births;
    std::shared_ptr<hn::FlowsRecord> record;
    std::string hall_record_sha256;
    std::size_t census = 0;
    std::map<std::uint64_t, m::PicVector> start;   // native census positions
};
hall_driver::MotionReceipt hall_motion(m::PrivateTransportTrial& trial, const m::FullTimeStepBinding& b, std::size_t cap, HallStep* keep, bool scope) {
    const auto census = trial.native_census(cap);
    if (keep) for (const auto& row : census) keep->start[row.native.id] = row.native.position_m;
    static const std::vector<double> E(mesh().ne, 0.0), B(mesh().nb, 0.0);
    hall_driver::MotionBinding mb;
    mb.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
    mb.accepted_step = b.accepted_step;
    mb.staged_arena_revision = trial.staged_arena_revision();
    mb.begin_s = b.begin_s;
    mb.end_s = b.end_s;
    mb.dt_s = b.dt_s;
    mb.staged_arena_sha256 = trial.state_prefix();
    mb.trial_field_sha256 = hall_driver::total_field_sha256(E, B);
    mb.mesh_sha256 = b.mesh_sha256;
    mb.material_sha256 = b.material_sha256;
    mb.species_sha256 = sha("hall-v3-native-v6-species");
    mb.common_clock_sha256 = b.common_clock_sha256;
    mb.kinetic_policy_sha256 = context().policy_sha256;
    std::vector<hc::Newborn> newborns;
    for (const auto& f : trial.birth_flights()) newborns.push_back({f.native_id, f.birth_time_s});
    return hc::MotionAssembly::cpu(context(), census, mb, {}, E, B, newborns, scope);
}
// scope: the native trial carries a birth scope (arrival committed every step), so the
// Hall record is the flows record even without newborns or sinks.
Hooks::Motion hall_hook(HallStep* keep, std::size_t cap = 4096, std::function<void(hn::FlowsRecord&)> mutate = {}, bool reseal = true, bool scope = false) {
    return [=](m::PrivateTransportTrial& trial, m::FullTimeMotionBinding& binding) {
        auto r = hall_motion(trial, binding.common, cap, keep, scope);
        auto x = hn::flows_record(r);
        if (keep) {
            keep->sinks.assign(r.sinks().begin(), r.sinks().end());
            keep->births.assign(r.births().begin(), r.births().end());
            keep->hall_record_sha256 = r.record_v3_sha256();
            keep->census = r.paths().size();
        }
        if (mutate) {
            mutate(*x);
            x->bind();
            if (reseal) x->record.record_sha256 = m::full_time_motion_record_v3_sha256(x->record);
        }
        if (keep) keep->record = x;
        binding.transport_trial_sha256 = x->record.record_sha256;
        return trial.prepare_total_field_motion(binding, m::bind_total_field_motion_v3(x->record, x));
    };
}
// Native sink ledger of the sealed step vs the Hall records and the Hall SinkLedger.
J compare_ledgers(const m::FullTimeSinkLedger* native, const HallStep& h, bool& ok) {
    J ev = J::object();
    const auto L = hc::SinkLedger::of(h.sinks);
    if (!native) {
        ok = h.sinks.empty();
        ev["native_ledger"] = "none";
        ev["hall_sinks"] = h.sinks.size();
        return ev;
    }
    bool records = native->records.size() == h.sinks.size() && native->carriers == h.sinks.size();
    for (std::size_t i = 0; records && i < h.sinks.size(); ++i) {
        const auto& a = native->records[i];
        const auto& k = h.sinks[i];
        records = a.native_id == k.native_id && int(a.kind) == k.kind && int(a.surface) == k.boundary && int(a.axis) == k.axis &&
                  same(a.plane_m, k.plane_m) && same(a.fraction, k.fraction) && same(a.segment_m, k.segment_m) &&
                  same(a.hit_point_m, k.hit_point_m) && same(a.incident_proper_velocity_m_s, k.incident_proper_velocity_m_s) &&
                  same(a.flight_s, k.flight_s) && same(a.discarded_flight_s, k.discarded_flight_s);
    }
    bool totals = native->surfaces.size() == 5;
    J per = J::object();
    for (std::size_t s = 0; totals && s < 5; ++s) {
        const auto& n = native->surfaces[s];
        const auto& t = L.boundary[int(n.surface)];
        const bool eq = n.carriers == t.count && n.weight == t.weight && n.kinetic_energy_J == t.kinetic_exact &&
                        n.positive_charge_e == t.positive_e && n.negative_charge_e == t.negative_e && t.non_elementary == 0;
        per[hc::sink_boundary_names[int(n.surface)]] = {{"carriers", n.carriers}, {"weight", n.weight.zero() ? 0. : n.weight.rounded_for_rate()},
                                                         {"kinetic_J", n.kinetic_energy_J.zero() ? 0. : n.kinetic_energy_J.rounded_for_rate()},
                                                         {"equal_to_hall", eq}};
        totals = totals && eq;
    }
    ev["records_equal"] = records;
    ev["surface_totals_equal"] = totals;
    ev["surfaces"] = per;
    ok = records && totals;
    return ev;
}

// ---- open step (zero-event first exposure) for refusal controls
struct OpenStep {
    m::FullTimeStepBinding b;
    m::PrivateFullTimeStep outer;
    std::optional<m::CollisionScope> first;
};
OpenStep open_step(m::PicArena& a, StepState& st, std::size_t cap = 64, std::size_t records = 64, double dt = 5e-15) {
    const auto k = a.full_time_boundary()->accepted_step;
    auto b = binding_of(a, dt, sha(step_text("v6-common-clock", k)));
    OpenStep s{b, a.begin_full_time_step(b, cap, records), std::nullopt};
    m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s, .5 * b.dt_s, st.sequence, b.accepted_field_sha256};
    s.first.emplace(s.outer.first_exposure(e1));
    auto t = s.first->prepare_compound_with_completion(1, [&](m::CollisionScope&) { return completion(e1, st.rng, sha(step_text("v6-rng-first", k))); });
    if (!s.first->try_publish(t)) throw std::runtime_error("first exposure");
    return s;
}
m::FullTimeMotionBinding motion_binding(const m::FullTimeStepBinding& b, const std::string& source_after) {
    return {b, "", b.accepted_field_sha256, sha(step_text("v6-field-after", b.accepted_step)), sha(step_text("v6-heat-after", b.accepted_step)), source_after};
}

// ---- arrival births (arrival owner v1 through the native birth scope)
std::shared_ptr<const pa::Policy> arrival_policy(double we = 4, double wn = 2) {
    pa::Configuration c;
    c.run_seed = 20261030;
    c.streams = {{0, pa::Kind::electron, 0, "e", we, 4096}, {1, pa::Kind::neutral, 0, "n", wn, 4096}};
    return std::make_shared<const pa::Policy>(c);
}
pa::StepInput arrival_input(std::uint64_t step, double begin, double dt, double electrons, double neutral_debit) {
    pa::StepInput in;
    in.step = step;
    if (electrons > 0) in.electrons = {{0, electrons, std::log(electrons), true, begin, dt, begin + dt}};
    if (neutral_debit > 0) in.neutrals = {{1, pa::Count(neutral_debit), std::log(neutral_debit), pa::Count(1e20), begin, begin + dt}};
    return in;
}
// Placement 0: upper-foil lower face (y = 7 mm) into the gap; 1: sink-prone (electrons
// 1e-7 m from the outer x face moving +x, neutrals 1e-12 m above the lower foil moving
// -y); 2: interior gap points; 3: lower-foil upper face (y = 2 mm) moving +y at 2e8 m/s
// (electrons; crossing the gap).
std::vector<av::BirthKinematics> placement(const pa::Staged& s, int kind) {
    std::vector<av::BirthKinematics> out;
    for (const auto& b : s.births()) {
        av::BirthKinematics k;
        const double dx = 1e-6 * double(b.index % 1000), dz = 1e-6 * double(b.index / 1000);
        const bool e = b.kind == pa::Kind::electron;
        if (kind == 0) { k.position_m = {.0030 + dx, .007, -.001 + dz + (e ? 0 : .002)}; k.proper_velocity_m_s = {0, e ? -2e6 : -650, 0}; }
        else if (kind == 1) {
            if (e) { k.position_m = {.009 - 1e-7, .0045 + dx, dz}; k.proper_velocity_m_s = {1e8, 0, 0}; }
            else { k.position_m = {.004 + dx, .002 + 1e-12, dz}; k.proper_velocity_m_s = {0, -650, 0}; }
        } else if (kind == 2) { k.position_m = {.004 + dx, .0045, .001 + dz + (e ? 0 : .001)}; k.proper_velocity_m_s = {e ? 1e6 : 300., 0, 0}; }
        else { k.position_m = {.0025 + 1e-4 * double(b.step % 40) + dx, .002, -.004 + 2e-4 * double(b.step % 40) + dz}; k.proper_velocity_m_s = {1e6, 2e8, -3e5}; }
        k.origin = std::string(e ? "source/electron-emission/event-" : "source/neutral-release/event-") + std::to_string(b.ordinal + 1);
        out.push_back(k);
    }
    return out;
}
void bind_neutral_times(pa::Staged& s) {
    for (const auto& d : s.plan().decisions())
        if (d.kind == pa::Kind::neutral && d.births) {
            std::vector<double> t;
            for (const auto& b : s.births())
                if (b.stream_id == d.stream_id) t.push_back(d.begin_s + b.marks[5] * (d.end_s - d.begin_s));
            s.bind_source_times(d.stream_id, t);
        }
}
struct ArrivalRun {
    std::shared_ptr<const pa::Policy> P;
    std::unique_ptr<pa::Owner> owner;
    int place = 0;
    double electrons = 6, neutral_debit = 3;
    std::optional<pa::Staged> staged;
    std::uint64_t source_index = 0;
    double scope_begin = 0, scope_dt = 0;
    explicit ArrivalRun(std::shared_ptr<const pa::Policy> p) : P(std::move(p)), owner(std::make_unique<pa::Owner>(P, 0)) {}
    Hooks hooks(Hooks::Motion motion, bool births_now) {
        Hooks h;
        h.births = [this, births_now](m::BirthScope& scope) -> std::optional<m::FullTimeBirthCredit> {
            const auto step = scope.accepted_step();
            auto plan = owner->plan(births_now ? arrival_input(step, scope_begin, scope_dt, electrons, neutral_debit) : arrival_input(step, scope_begin, scope_dt, 0, 0));
            staged.emplace(owner->stage(std::move(plan), scope.base_native_id()));
            bind_neutral_times(*staged);
            const auto kin = placement(*staged, place);
            return av::credit(*P, *staged, kin, av::plan_debits(staged->plan()), scope.source_before_sha256(),
                              sha(step_text("source", ++source_index)));
        };
        h.motion = std::move(motion);
        return h;
    }
    StepResult step(m::PicArena& a, StepState& st, Hooks::Motion motion, bool births_now = true, std::size_t cap = 4096) {
        StepPlan plan;
        plan.maximum_carriers = cap;
        plan.maximum_records = 64;
        scope_begin = a.full_time_boundary()->time_s;
        scope_dt = plan.dt_s;
        auto r = full_step(a, plan, st, hooks(std::move(motion), births_now));
        if (!r.flows || !r.flows->births) throw std::runtime_error("arrival step without native birth receipt");
        (void)owner->commit(std::move(*staged), av::receipt(*r.flows->births));
        staged.reset();
        return r;
    }
};
// Arrival owner v1.1 through the native birth scope (upper-foil face births).
struct ArrivalRun11 {
    std::shared_ptr<const pa11::Policy> P;
    std::unique_ptr<pa11::Owner> owner;
    std::optional<pa11::Staged> staged;
    double electrons = 0, neutral_debit = 0;
    std::uint64_t source_index = 0;
    double scope_begin = 0, scope_dt = 0;
    explicit ArrivalRun11(std::shared_ptr<const pa11::Policy> p) : P(std::move(p)), owner(std::make_unique<pa11::Owner>(P, 0)) {}
    StepResult step(m::PicArena& a, StepState& st, Hooks::Motion motion) {
        StepPlan plan;
        plan.maximum_carriers = 4096;
        plan.maximum_records = 64;
        scope_begin = a.full_time_boundary()->time_s;
        scope_dt = plan.dt_s;
        Hooks h;
        h.births = [this](m::BirthScope& scope) -> std::optional<m::FullTimeBirthCredit> {
            pa11::StepInput in;
            in.step = scope.accepted_step();
            if (electrons > 0) in.electrons = {{0, electrons, std::log(electrons), true, scope_begin, scope_dt, scope_begin + scope_dt}};
            if (neutral_debit > 0)
                in.neutrals = {{1, pa11::Count(neutral_debit), std::log(neutral_debit), pa11::Count(1e20), scope_begin, scope_begin + scope_dt}};
            staged.emplace(owner->stage(owner->plan(in), scope.base_native_id()));
            for (const auto& d : staged->plan().decisions())
                if (d.kind == pa11::Kind::neutral && d.births) {
                    std::vector<double> t;
                    for (const auto& b : staged->births())
                        if (b.stream_id == d.stream_id) t.push_back(d.begin_s + b.marks[5] * (d.end_s - d.begin_s));
                    staged->bind_source_times(d.stream_id, t);
                }
            std::vector<av11::BirthKinematics> kin;
            for (const auto& b : staged->births()) {
                av11::BirthKinematics k;
                const double dx = 1e-6 * double(b.index % 1000), dz = 1e-5 * double(b.step % 400);
                const bool e = b.kind == pa11::Kind::electron;
                k.position_m = {.0030 + dx, .007, -.002 + dz + (e ? 0 : .001)};
                k.proper_velocity_m_s = {0, e ? -2e6 : -650, 0};
                k.origin = std::string(e ? "source/electron-emission/event-" : "source/neutral-release/event-") + std::to_string(b.ordinal + 1);
                kin.push_back(k);
            }
            return av11::credit(*P, *staged, kin, av11::plan_debits(staged->plan()), scope.source_before_sha256(),
                                sha(step_text("source-v1.1", ++source_index)));
        };
        h.motion = std::move(motion);
        auto r = full_step(a, plan, st, h);
        if (!r.flows || !r.flows->births) throw std::runtime_error("arrival v1.1 step without native birth receipt");
        (void)owner->commit(std::move(*staged), av11::receipt(*r.flows->births));
        staged.reset();
        return r;
    }
};
bool births_equal(const m::FullTimeBirthReceipt* rec, const HallStep& h) {
    if (!rec) return h.births.empty();
    if (rec->entries.size() != h.births.size()) return false;
    for (std::size_t i = 0; i < h.births.size(); ++i)
        if (rec->entries[i].native_id != h.births[i].native_id || !same(rec->entries[i].weight, h.births[i].weight)) return false;
    return true;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2 && argc != 3) throw std::runtime_error("usage: test_native_v6_admission OUTPUT_DIRECTORY [STEPS]");
        const std::string dir = argv[1];
        const long steps = argc == 3 ? std::stol(argv[2]) : 12000;
        const auto t0 = std::chrono::steady_clock::now();
        // 1. plain v3 record (v2-identical Hall v3 receipt)
        {
            auto a = fresh();
            StepState st;
            HallStep h;
            Hooks hk;
            hk.motion = hall_hook(&h);
            auto r = full_step(*a, StepPlan{}, st, hk);
            check(h.record && h.record->record.schema == std::string(m::motion_v3_schema) && h.sinks.empty() && h.births.empty() && (!r.flows || !r.flows->sinks),
                  "plain: a Hall v3 receipt without newborns and hits is the plain binary v3 record, admitted and published by native v6");
        }
        // 2. sinks of every surface
        {
            auto a = fresh(sink_population());
            StepState st;
            HallStep h;
            Hooks hk;
            hk.motion = hall_hook(&h);
            auto r = full_step(*a, StepPlan{}, st, hk);
            bool ok = false;
            auto ev = compare_ledgers(r.flows ? r.flows->sinks.get() : nullptr, h, ok);
            std::set<int> kinds;
            for (const auto& k : h.sinks) kinds.insert(k.boundary);
            ev["surfaces_hit"] = J(std::vector<int>(kinds.begin(), kinds.end()));
            ev["hall_record_sha256"] = h.hall_record_sha256;
            measurements["sinks"] = ev;
            check(h.record->record.schema == std::string(m::motion_v3_flows_schema) && h.sinks.size() >= 9 && kinds.size() == 5,
                  "sinks: Hall v3 truncates carriers at all five kinetic surfaces (flows record)");
            check(ok, "sinks: native admits the Hall flows record; native sink records equal the Hall records and native per-surface totals "
                      "equal the Hall SinkLedger exactly (carriers, weight, kinetic energy, +/- elementary charge)");
            bool port = true;
            for (const auto& k : h.record->sinks) {
                // re-derive from the census start and the Hall segment
                const auto hit = m::full_time_first_kinetic_hit(options(), h.start.at(k.native_id), k.segment_m);
                port = port && hit.surface == int(k.surface) && same(hit.fraction, k.fraction) && hit.axis == int(k.axis) && same(hit.plane_m, k.plane_m);
            }
            check(port, "sinks: native full_time_first_kinetic_hit of every Hall sink segment equals the Hall surface, fraction, axis and plane");
        }
        // 3. births: upper-foil face, sink-prone, interior
        for (int place : {0, 1, 2}) {
            auto P = arrival_policy();
            auto a = fresh();
            StepState st;
            ArrivalRun run(P);
            run.place = place;
            run.electrons = 4 * 3.3;
            run.neutral_debit = 5;
            HallStep h;
            auto r = run.step(*a, st, hall_hook(&h, 4096, {}, true, true));
            bool ok = false;
            auto ev = compare_ledgers(r.flows->sinks.get(), h, ok);
            ev["births"] = h.births.size();
            std::size_t face = 0;
            for (const auto& b : h.births) face += b.origin >= 0;
            ev["face_births"] = face;
            bool partial = true;
            for (const auto& k : h.sinks) {
                const auto it = std::find_if(h.births.begin(), h.births.end(), [&](const auto& b) { return b.native_id == k.native_id; });
                if (it != h.births.end()) partial = partial && same(k.flight_s + k.discarded_flight_s, it->flight_s) && it->flight_s < 5e-15;
            }
            ev["newborn_sinks_within_flight"] = partial;
            measurements["births-place-" + std::to_string(place)] = ev;
            const char* what[] = {"upper-foil face births moving into the gap", "sink-prone births (sinks within end - tau)", "interior births"};
            check(!h.births.empty() && births_equal(r.flows->births.get(), h) && ok && partial && (place != 0 || face == h.births.size()) &&
                      (place != 1 || !h.sinks.empty()),
                  std::string("births: ") + what[place] + ": native admits the Hall flows record (births with flight end - tau)");
        }
        // 4. births and sinks in one step
        {
            auto P = arrival_policy();
            auto a = fresh(sink_population());
            StepState st;
            ArrivalRun run(P);
            run.place = 1;
            run.electrons = 4 * 3.3;
            run.neutral_debit = 5;
            HallStep h;
            auto r = run.step(*a, st, hall_hook(&h, 4096, {}, true, true));
            bool ok = false;
            auto ev = compare_ledgers(r.flows->sinks.get(), h, ok);
            std::size_t newborn_sunk = 0;
            for (const auto& k : h.sinks) newborn_sunk += k.native_id >= r.flows->births->base;
            ev["newborn_sinks"] = newborn_sunk;
            measurements["same-step"] = ev;
            check(ok && newborn_sunk > 0 && h.sinks.size() > newborn_sunk && births_equal(r.flows->births.get(), h),
                  "same step: resident and newborn sinks with births, native admits; ledgers equal");
        }
        // 5. tampered Hall records are refused (typed)
        {
            // one fresh arena and open step per control (native bounds the forks of one step)
            auto run = [&](std::function<void(hn::FlowsRecord&)> mutate, bool reseal) {
                auto a = fresh(sink_population());
                StepState st;
                auto o = open_step(*a, st);
                auto trial = o.outer.fork_transport_trial();
                auto mb = motion_binding(o.b, a->full_time_boundary()->source_sha256);
                (void)hall_hook(nullptr, 64, mutate, reseal)(trial, mb);
            };
            const double dt = 5e-15;   // open_step interval
            check(flow_refuses([&] { run([](hn::FlowsRecord& x) { x.sinks[0].fraction = std::nextafter(x.sinks[0].fraction, 2.); }, true); }, FK::sink_outside_step),
                  "tamper: one-ulp change of a Hall sink fraction alone (resealed) is refused (sink_outside_step: flight != F * fraction)");
            check(flow_refuses([&] { run([dt](hn::FlowsRecord& x) { auto& k = x.sinks[0]; k.fraction = std::nextafter(k.fraction, 0.);
                                                                   k.flight_s = dt * k.fraction; k.discarded_flight_s = dt - k.flight_s; }, true); }, FK::sink_not_crossing),
                  "tamper: one-ulp change of a Hall sink fraction with consistent flight (resealed) is refused (sink_not_crossing)");
            check(flow_refuses([&] { run([](hn::FlowsRecord& x) { x.sinks[0].hit_point_m[0] = std::nextafter(x.sinks[0].hit_point_m[0], 1.); }, true); }, FK::sink_record),
                  "tamper: one-ulp change of a Hall hit point (resealed) is refused (sink_record)");
            check(flow_refuses([&] { run([](hn::FlowsRecord& x) { x.sinks[0].flight_s *= 1 + 0x1p-50; }, false); }, FK::motion_record),
                  "tamper: a Hall record changed without resealing is refused (motion_record digest)");
            check(flow_refuses([&] { run([](hn::FlowsRecord& x) { x.record.coverage = "all_native_carriers_no_sources_no_boundaries"; }, true); }, FK::motion_record),
                  "tamper: flows record with the plain coverage is refused (motion_record)");
        }
        // 6. multi-step: electrons born on the lower-foil face cross the gap and are collected at the upper foil
        {
            auto P = arrival_policy(4, 2);
            auto a = fresh({});
            StepState st;
            ArrivalRun run(P);
            run.place = 3;
            run.electrons = 4 * 2.5;
            run.neutral_debit = 0;
            const long birth_steps = 40;
            std::uint64_t born = 0, collected_upper = 0, other = 0, admitted = 0, ledger_equal = 0;
            hc::SinkLedger cumulative;
            std::ofstream rows(dir + "/native-v6-multistep.jsonl");
            const auto t1 = std::chrono::steady_clock::now();
            long done = 0;
            for (long s = 0; s < steps; ++s) {
                HallStep h;
                auto r = run.step(*a, st, hall_hook(&h, 4096, {}, true, true), s < birth_steps);
                ++admitted;
                bool ok = false;
                auto ev = compare_ledgers(r.flows->sinks.get(), h, ok);
                ledger_equal += ok;
                born += h.births.size();
                for (const auto& k : h.sinks) (k.boundary == 1 ? collected_upper : other) += 1;
                cumulative.add(hc::SinkLedger::of(h.sinks));
                ++done;
                if (!h.births.empty() || !h.sinks.empty())
                    rows << J({{"step", s}, {"census", h.census}, {"births", h.births.size()}, {"sinks", h.sinks.size()}, {"ledger_equal", ok},
                               {"record_schema", h.record->record.schema}})
                                .dump()
                         << '\n';
                if (s >= birth_steps && h.census == 0) break;
            }
            const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
            const auto totals = a->full_time_flow_totals();
            measurements["multistep"] = {{"steps", done}, {"births", born}, {"collected_upper_foil", collected_upper}, {"other_sinks", other},
                                         {"native_births", totals.births}, {"native_sunk_carriers", totals.sunk_carriers},
                                         {"steps_with_equal_ledgers", ledger_equal}, {"cumulative_hall_sink_ledger", cumulative.json()},
                                         {"wall_s", wall}, {"ms_per_step", 1e3 * wall / std::max<long>(1, done)}};
            a->validate_full_time_history();
            check(admitted == std::uint64_t(done) && ledger_equal == std::uint64_t(done) && born > 0 && collected_upper == born && other == 0 &&
                      totals.births == born && totals.sunk_carriers == born,
                  "multistep: every step admitted by native v6; electrons born on the lower-foil face cross the gap and are all collected at the upper "
                  "foil; per-step native ledgers equal the Hall ledgers; native cumulative births = sunk carriers");
        }
        // 7. arrival policy v1.1: quantized neutral leaves through the native v6 birth scope
        {
            const double w0 = 2;
            pa11::Configuration c;
            c.run_seed = 20261030;
            c.streams = {{0, pa11::Kind::electron, 0, "e", 4, 4096}, {1, pa11::Kind::neutral, 0, "n", w0, 4096}};
            auto P = std::make_shared<const pa11::Policy>(c);
            auto a = fresh();
            StepState st;
            ArrivalRun11 run(P);
            run.electrons = 4 * 1.5;
            run.neutral_debit = 5;
            const std::uint64_t n_steps = 40;
            std::uint64_t admitted = 0, neutral = 0, equal = 0;
            bool all_w0 = true, debit_exact = true;
            for (std::uint64_t s = 0; s < n_steps; ++s) {
                HallStep h;
                auto r = run.step(*a, st, hall_hook(&h, 4096, {}, true, true));
                ++admitted;
                const auto& R = *r.flows->births;
                std::uint64_t nn = 0;
                for (const auto& e : R.entries)
                    if (e.state == "n") { ++nn; all_w0 = all_w0 && e.weight == w0; }
                for (const auto& x : R.streams)
                    if (x.stream_id == 1) debit_exact = debit_exact && (nn ? x.source_debit == pa11::Count(w0 * double(nn)) : x.source_debit.zero());
                neutral += nn;
                equal += births_equal(&R, h);
            }
            const double expected = double(n_steps) * 5 / w0;
            measurements["arrival-v1.1"] = {{"steps", admitted}, {"neutral_leaves", neutral}, {"expected_leaves", expected}, {"all_leaves_w0", all_w0},
                                            {"stream_debit_equals_n_w0", debit_exact}, {"hall_births_equal_native", equal},
                                            {"native_births_total", a->full_time_flow_totals().births}};
            a->validate_full_time_history();
            check(admitted == n_steps && all_w0 && debit_exact && equal == n_steps && neutral > 0 && std::abs(double(neutral) - expected) <= 5 * std::sqrt(.25 * n_steps),
                  "arrival v1.1: quantized neutral leaves of exactly w0 on the upper-foil face through the native v6 birth scope; native stream debit "
                  "= n w0 exactly; Hall flows records admitted every step; leaf count consistent with D/w0");
        }
        std::ofstream(dir + "/native-v6-admission.json")
            << J({{"status", "passed"}, {"check_count", checks}, {"checks", results}, {"measurements", measurements},
                  {"seconds", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()}})
                   .dump(1)
            << '\n';
        std::cout << checks << " native v6 admission checks passed\n";
        return 0;
    } catch (const std::exception& x) {
        std::cerr << "FAILED after " << checks << " checks: " << x.what() << '\n';
        if (argc >= 2)
            std::ofstream(std::string(argv[1]) + "/native-v6-admission.json")
                << J({{"status", "failed"}, {"check_count", checks}, {"checks", results}, {"failure", x.what()}, {"measurements", measurements}}).dump(1) << '\n';
        return 1;
    }
}
