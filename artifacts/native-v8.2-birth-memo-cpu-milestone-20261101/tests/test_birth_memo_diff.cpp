// Bounded host birth/sink differential fixture against two actual release archives.
// Exact physical bytes cover population, motion records, birth/source credits,
// sink exposures, exact ledgers and RNG history. Provenance-bound hashes remain
// separate. Invalid, moved/stale-scope and changed-cap refusals must preserve
// accepted state and agree between binaries. No production timing forecast.
#include "v6_fixture.hpp"
#include "arrival_v1_births.hpp"
#include "forms_gpu.hpp"
#include "hall_gpu_common.hpp"
#include "raikou/mfem_pic_full_time_device.h"
#include <filesystem>
#include <fstream>
#include <iostream>

using namespace v6;
namespace pa = pf_arrival;
namespace av = m::arrival_v1;
namespace {
std::size_t carrier_count = 10000, cap = 12000;
double birth_seconds = 0;
int n_checks = 0, n_fail = 0;
void need(bool ok, const std::string &what, const J &d = J::object()) {
  ++n_checks;
  if (!ok) {
    ++n_fail;
    std::cout << "FAIL " << what << " " << d.dump() << "\n";
  }
}
std::shared_ptr<const pa::Policy> arrival_policy() {
  pa::Configuration c;
  c.run_seed = 20261031;
  c.streams = {{0, pa::Kind::electron, 0, "e", 4, 4096}, {1, pa::Kind::neutral, 0, "n", 2, 4096}};
  return std::make_shared<const pa::Policy>(c);
}
pa::StepInput arrival_input(std::uint64_t step, double begin, double dt, double electrons, double neutral_debit) {
  pa::StepInput in;
  in.step = step;
  in.electrons = {{0, electrons, std::log(electrons), true, begin, dt, begin + dt}};
  in.neutrals = {{1, pa::Count(neutral_debit), std::log(neutral_debit), pa::Count(1e20), begin, begin + dt}};
  return in;
}
// Placement 0: emitting-face births (upper foil lower face, into the gap); 1: sink-prone births
// (electrons at the outer x face moving +x, sinking in their first flight; neutrals just above the
// lower foil moving -y); 2: interior births.
std::vector<av::BirthKinematics> placement(const pa::Staged &s, int kind) {
  std::vector<av::BirthKinematics> out;
  for (const auto &b : s.births()) {
    av::BirthKinematics k;
    const double dx = 1e-6 * double(b.index % 1000), dz = 1e-6 * double(b.index / 1000);
    const bool e = b.kind == pa::Kind::electron;
    if (kind == 0) {
      k.position_m = {.0030 + dx, .007, -.001 + dz + (e ? 0 : .002)};
      k.proper_velocity_m_s = {0, e ? -2e6 : -650, 0};
    } else if (kind == 1) {
      if (e) {
        k.position_m = {.009 - 1e-7, .0045 + dx, dz};
        k.proper_velocity_m_s = {1e8, 0, 0};
      } else {
        k.position_m = {.004 + dx, .002 + 1e-12, dz};
        k.proper_velocity_m_s = {0, -650, 0};
      }
    } else {
      k.position_m = {.004 + dx, .0045, .001 + dz + (e ? 0 : .001)};
      k.proper_velocity_m_s = {e ? 1e6 : 300., 0, 0};
    }
    k.origin = std::string(e ? "source/electron-emission/event-" : "source/neutral-release/event-") +
               std::to_string(b.ordinal + 1);
    out.push_back(k);
  }
  return out;
}
void bind_neutral_times(pa::Staged &s) {
  for (const auto &d : s.plan().decisions())
    if (d.kind == pa::Kind::neutral && d.births) {
      std::vector<double> t;
      for (const auto &b : s.births())
        if (b.stream_id == d.stream_id)
          t.push_back(d.begin_s + b.marks[5] * (d.end_s - d.begin_s));
      s.bind_source_times(d.stream_id, t);
    }
}
std::vector<m::FullTimeAuthoredPacket> population() {
  auto p = base_population();
  // co-located two-carrier packets: one leaving through the outer x face (sunk whole), one staying
  p.push_back({"pair-out", "e", "manufactured/source", {1, 1e-20}, {.009 - 2e-7, .0050, .002}, {1e8, 0, 0}});
  p.push_back({"pair-in", "n", "manufactured/source", {2, 1e-19}, {.0050, .0050, -.002}, {0, 0, 300}});
  for (std::size_t k = 0; k < carrier_count; ++k)
    p.push_back({"bulk-" + std::to_string(k), k % 2 ? "e" : "i", "manufactured/source", {1.5 + k, 0},
                 {.0025 + 1e-5 * (k % 300), .0045 + 1e-6 * (k % 100), -.003 + 1e-5 * (k % 600)},
                 {k % 3 ? 2e5 : -1e5, 1e4 * (k % 5), 3e4}});
  return p;
}
StepResult timed_full_step(m::PicArena &a, const StepPlan &opt, StepState &st,
                            const Hooks &hooks) {
  const auto before = a.full_time_boundary();
  const std::uint64_t k = before->accepted_step;
  auto b = binding_of(a, opt.dt_s, sha(step_text("v6-common-clock", k)));
  auto outer = a.begin_full_time_step(b, opt.maximum_carriers, opt.maximum_records);
  m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s,
                        .5 * b.dt_s, st.sequence, b.accepted_field_sha256};
  const auto rng_a = sha(step_text("v6-rng-first", k)),
             rng_b = sha(step_text("v6-rng-second", k));
  const auto budget = std::size_t(std::max({1, opt.first_events, opt.second_events}));
  auto s1 = outer.first_exposure(e1);
  auto t1 = s1.prepare_compound_with_completion(budget, [&](m::CollisionScope &) {
    return completion(e1, st.rng, rng_a);
  });
  if (!s1.try_publish(t1))
    throw std::runtime_error("first exposure publication");
  StepResult result;
  std::optional<m::PreparedBirths> births;
  std::string source_after = before->source_sha256;
  if (hooks.births) {
    auto scope = outer.birth_scope();
    result.birth_base = scope.base_native_id();
    if (auto credit = hooks.births(scope)) {
      auto timer = std::chrono::steady_clock::now();
      births.emplace(scope.prepare(*credit));
      birth_seconds += since(timer);
      source_after = credit->source_after_sha256;
    }
  }
  auto trial = births ? outer.fork_transport_trial(*births) : outer.fork_transport_trial();
  const auto field_after = sha(step_text("v6-field-after", k));
  const auto heat_after = sha(step_text("v6-heat-after", k));
  m::FullTimeMotionBinding mb{b, "", b.accepted_field_sha256, field_after,
                              heat_after, source_after};
  auto moved = hooks.motion(trial, mb);
  outer.accept_transport_trial(std::move(moved));
  m::ExposureBinding e2 = e1;
  e2.slot = m::FullTimeSlot::after_transport;
  e2.frozen_physical_stage_s = b.end_s;
  e2.owner_stage_sequence = st.sequence + 1;
  e2.frozen_field_sha256 = field_after;
  auto s2 = outer.second_exposure(e2);
  auto t2 = s2.prepare_compound_with_completion(budget, [&](m::CollisionScope &) {
    return completion(e2, rng_a, rng_b);
  });
  if (!s2.try_publish(t2))
    throw std::runtime_error("second exposure publication");
  auto token = outer.seal({field_after, heat_after, source_after, rng_b});
  result.flows = token.candidate_full_time_flows();
  if (hooks.before_publish)
    hooks.before_publish(token);
  if (!a.can_publish(token) || !a.try_publish(token))
    throw std::runtime_error("outer publication refused");
  st.sequence += 2;
  st.rng = rng_b;
  return result;
}

struct Harness {
  std::shared_ptr<const pa::Policy> P = arrival_policy();
  std::unique_ptr<pa::Owner> owner = std::make_unique<pa::Owner>(P, 0);
  std::optional<pa::Staged> staged;
  std::uint64_t source_index = 0;
  int place = 0;
  std::function<void(m::FullTimeBirthCredit &)> mutate;
  double begin = 0, dt = 0;
  Hooks hooks(Hooks::Motion motion) {
    Hooks h;
    h.births = [this](m::BirthScope &scope) -> std::optional<m::FullTimeBirthCredit> {
      auto plan = owner->plan(arrival_input(scope.accepted_step(), begin, dt, 6, 1.5));
      staged.emplace(owner->stage(std::move(plan), scope.base_native_id()));
      bind_neutral_times(*staged);
      auto c = av::credit(*P, *staged, placement(*staged, place), av::plan_debits(staged->plan()),
                        scope.source_before_sha256(), sha("source/" + std::to_string(++source_index)));
      if (mutate) mutate(c);
      return c;
    };
    h.motion = std::move(motion);
    return h;
  }
  StepResult step(m::PicArena &a, StepState &st, Hooks::Motion motion = reference_hook(cap)) {
    StepPlan plan;
    plan.dt_s = 5e-15;
    plan.maximum_carriers = cap;
    begin = a.full_time_boundary()->time_s;
    dt = plan.dt_s;
    auto r = timed_full_step(a, plan, st, hooks(std::move(motion)));
    (void)owner->commit(std::move(*staged), av::receipt(*r.flows->births));
    staged.reset();
    return r;
  }
};
std::unique_ptr<m::PicArena> fresh_v8() {
  m::FullTimeConfiguration c;
  c.identity = m::FullTimeIdentitySchema::v8;
  c.limits.maximum_carriers = cap;
  c.limits.memory_budget_bytes = cap * m::full_time_bytes_per_carrier + m::full_time_fixed_bytes;
  c.history.anchor_interval_steps = 16;
  return m::PicArena::create_full_time(options(), states(), channels(), population(), genesis(), c);
}
void count(Bytes &b, const m::ArenaCount &x) {
  const auto n = x.copy_units();
  char *hex = BN_bn2hex(n.get());
  if (!hex) throw std::bad_alloc();
  b.text(hex);
  OPENSSL_free(hex);
}
void particle(Bytes &b, const m::GapParticle &p) {
  b.u64(p.id); b.vec(p.position_m); b.vec(p.proper_velocity_m_s);
  b.f64(p.charge_C); b.f64(p.mass_kg); b.f64(p.weight);
}
void sink(Bytes &b, const m::FullTimeSinkRecord &s) {
  b.u64(s.native_id); b.u64(unsigned(s.kind)); b.u64(unsigned(s.surface)); b.u64(s.axis);
  b.f64(s.plane_m); b.f64(s.fraction); b.vec(s.segment_m); b.vec(s.hit_point_m);
  b.vec(s.incident_proper_velocity_m_s); b.f64(s.flight_s); b.f64(s.discarded_flight_s);
}
void surface(Bytes &b, const m::FullTimeSinkSurfaceTotals &s) {
  b.u64(unsigned(s.surface)); b.u64(unsigned(s.kind)); b.u64(s.carriers);
  count(b, s.weight); count(b, s.kinetic_energy_J); count(b, s.positive_charge_e); count(b, s.negative_charge_e);
}
// Every binary64 is emitted as its exact little-endian bits, never decimal.
// Counts use the exact integer in 2^-1074 units. No object padding is copied.
Bytes physical(const m::PicArena &a, const StepResult &r, const ReferenceMotion &motion, const StepState &st) {
  Bytes b;
  b.text("birth-memo-physical-differential-v1");
  b.u64(a.full_time_boundary()->accepted_step); b.f64(a.full_time_boundary()->time_s);
  b.text(st.rng); b.u64(st.sequence); b.u64(st.events);
  const auto v = a.full_time_population_view();
  b.u64(v.rows.size()); b.u64(v.states.size());
  for (const auto &s : v.states) b.text(s);
  for (std::size_t k = 0; k < v.rows.size(); ++k) {
    particle(b, v.rows[k]); b.u64(v.packet_keys[k]); b.u64(v.state_indices[k]);
  }
  for (const auto &s : v.totals) { b.u64(s.carriers); b.f64(s.weight); b.f64(s.charge_C); }
  b.u64(motion.census.size()); for (const auto &p : motion.census) particle(b, p);
  b.u64(motion.endpoints.size()); for (const auto &p : motion.endpoints) particle(b, p);
  b.u64(motion.paths.size());
  for (const auto &p : motion.paths) { b.u64(p.native_id); b.vec(p.field_start_m); b.vec(p.field_displacement_m); b.f64(p.flight_s); }
  b.u64(motion.births.size());
  for (const auto &p : motion.births) { b.u64(p.native_id); b.f64(p.birth_time_s); b.f64(p.flight_s); }
  for (double x : motion.record.budget) b.f64(x);
  for (double x : motion.record.kinetic_policy) b.f64(x);
  for (auto array : {motion.record.trial_E, motion.record.trial_B, motion.record.current_A_m}) {
    b.u64(array.size()); for (double x : array) b.f64(x);
  }
  const auto &birth = *r.flows->births;
  b.u64(birth.accepted_step); b.u64(birth.base); b.u64(birth.next_id_after); b.u64(birth.arrival_next_id_after);
  b.text(birth.plan_sha256); b.text(birth.staged_token_sha256);
  b.text(birth.source_before_sha256); b.text(birth.source_after_sha256);
  b.u64(birth.entries.size());
  for (const auto &e : birth.entries) { b.text(e.packet_id); b.text(e.state); b.u64(e.native_id); b.f64(e.weight); b.u64(e.returned); }
  b.u64(birth.streams.size());
  for (const auto &s : birth.streams) { b.u64(s.stream_id); b.text(s.state); count(b, s.source_debit); }
  b.text(birth.internal_energy_J_digits); b.f64(birth.internal_energy_J);
  b.u64(bool(r.flows->sinks));
  if (r.flows->sinks) {
    const auto &s = *r.flows->sinks;
    b.u64(s.accepted_step); b.u64(s.carriers); b.u64(s.records.size());
    for (const auto &x : s.records) sink(b, x);
    b.u64(s.removed_packets.size()); for (const auto &x : s.removed_packets) b.text(x);
    b.u64(s.tallies.size());
    for (const auto &x : s.tallies) {
      b.u64(unsigned(x.surface)); b.u64(unsigned(x.kind)); b.text(x.state); b.u64(x.charge_e); b.u64(x.carriers);
      count(b, x.weight); count(b, x.kinetic_energy_J);
    }
    b.u64(s.surfaces.size()); for (const auto &x : s.surfaces) surface(b, x);
    b.u64(s.exposed.size());
    for (const auto &x : s.exposed) {
      sink(b, x.record); b.text(x.packet_id); b.text(x.state); b.f64(x.weight);
      b.f64(x.charge_C); b.f64(x.mass_kg); b.vec(x.start_position_m); b.f64(x.internal_energy_J);
    }
  }
  const auto t = a.full_time_flow_totals();
  for (auto n : {t.accepted_step, t.next_native_id, t.births, t.birth_rows, t.sunk_carriers, t.sink_rows,
                 std::uint64_t(t.has_last_birth), t.last_birth_step, t.last_birth_base, t.last_birth_count}) b.u64(n);
  b.text(t.last_birth_plan_sha256); b.text(t.last_birth_token_sha256);
  b.u64(t.surfaces.size()); for (const auto &x : t.surfaces) surface(b, x);
  return b;
}
template <class F> J refusal(F f) {
  try { f(); }
  catch (const m::FullTimeFlowRefusal &e) { return {{"type", "flow"}, {"kind", int(e.kind())}, {"text", e.what()}}; }
  catch (const m::FullTimeResourceRefusal &e) { return {{"type", "resource"}, {"kind", int(e.kind())}, {"text", e.what()}}; }
  catch (const std::exception &e) { return {{"type", "other"}, {"kind", -1}, {"text", e.what()}}; }
  return {{"type", "none"}};
}
J refused_controls() {
  J rows = J::array();
  using C = m::FullTimeBirthCredit;
  const std::vector<std::pair<std::string, std::function<void(C &)>>> bad = {
      {"source binding", [](C &c) { c.source_before_sha256 = sha("wrong source"); }},
      {"step binding", [](C &c) { ++c.arrival_step; }},
      {"native id", [](C &c) { ++c.births.at(0).native_id; }},
      {"packet id", [](C &c) { c.births.at(0).packet_id = "electron"; }},
      {"state", [](C &c) { c.births.at(0).state = "missing"; }},
      {"zero weight", [](C &c) { c.births.at(0).weight = 0; }},
      {"source debit", [](C &c) { c.streams.at(0).source_debit = c.streams.at(0).source_debit.credited(1); }},
      {"nonfinite velocity", [](C &c) { c.births.at(0).proper_velocity_m_s[0] = INFINITY; }},
      {"conductor endpoint", [](C &c) { c.births.at(0).position_m = {.004, .0071, 0}; }},
      {"birth time", [](C &c) { c.births.at(0).birth_time_s -= 1e-10; }},
      {"internal energy", [](C &c) { c.births.at(0).internal_energy_J = -1; }},
      {"duplicate birth", [](C &c) { c.births.insert(c.births.begin() + 1, c.births.at(0)); }}};
  for (const auto &[name, mutate] : bad) {
    auto a = fresh_v8(); Harness h; StepState st; h.mutate = mutate;
    const auto before = a->state_prefix();
    auto row = refusal([&] { (void)h.step(*a, st); });
    row["name"] = name; row["unchanged"] = a->state_prefix() == before && st.sequence == 0 && st.rng == genesis().rng_sha256;
    const bool expected = name == "nonfinite velocity"
        ? row["type"] == "other" && row["text"] == "nonfinite native product velocity"
        : row["type"] == "flow";
    need(expected && row["unchanged"].get<bool>(), "native refusal: " + name, row);
    rows.push_back(row);
  }
  // Moved/stale scope and cap mismatch exercise the memo's owner/scope boundary.
  auto a = fresh_v8(); Harness h; const auto before = a->state_prefix();
  std::optional<m::BirthScope> expired;
  std::optional<C> good;
  {
    auto b = binding_of(*a, 5e-15, sha("scope-control"));
    auto outer = a->begin_full_time_step(b, cap, 32);
    m::ExposureBinding e{b, m::FullTimeSlot::before_transport, b.begin_s, .5*b.dt_s, 0, b.accepted_field_sha256};
    auto first = outer.first_exposure(e);
    auto token = first.prepare_compound_with_completion(1, [&](auto &) { return completion(e, genesis().rng_sha256, sha("rng-a")); });
    if (!first.try_publish(token)) throw std::runtime_error("scope first exposure");
    auto scope = outer.birth_scope(); h.begin=b.begin_s; h.dt=b.dt_s;
    good = h.hooks(reference_hook(cap)).births(scope);
    (void)scope.prepare(*good); // valid before destruction, and primes validated immutable base
    expired.emplace(std::move(scope));
    auto row = refusal([&] { (void)scope.prepare(*good); });
    row["name"]="moved scope"; row["unchanged"]=a->state_prefix()==before;
    need(row["type"]!="none" && row["unchanged"].get<bool>(), "moved scope",row); rows.push_back(row);
    row = refusal([&] { (void)a->begin_full_time_step(b, cap, 32); });
    row["name"]="retained scope blocks second outer"; row["unchanged"]=a->state_prefix()==before;
    need(row["type"]=="other" && row["text"]=="full-time outer requires accepted full-time boundary" &&
         row["unchanged"].get<bool>(), "retained scope blocks second outer",row); rows.push_back(row);
    // Retaining a BirthScope retains its outer frame; a second outer must not
    // start. Advance the phase of this same outer to test a stale birth scope.
    auto births = expired->prepare(*good);
    auto trial = outer.fork_transport_trial(births);
    m::FullTimeMotionBinding mb{b, "", b.accepted_field_sha256,
        sha("scope-field-after"), sha("scope-heat-after"), good->source_after_sha256};
    auto moved = reference_hook(cap)(trial, mb);
    outer.accept_transport_trial(std::move(moved));
    row = refusal([&] { (void)expired->prepare(*good); });
    row["name"]="stale scope after transport phase"; row["unchanged"]=a->state_prefix()==before;
    need(row["type"]=="flow" && row["unchanged"].get<bool>(), "stale phase scope",row); rows.push_back(row);
  }
  // The scope owns the outer frame, but must not access a destroyed arena.
  a.reset();
  auto row = refusal([&] { (void)expired->prepare(*good); });
  row["name"]="stale scope after owner destruction";
  need(row["type"]=="other" && row["text"]=="stale full-time outer owner", "stale owner scope",row); rows.push_back(row);
  expired.reset();
  const auto old_cap=cap;
  auto capped=fresh_v8(); const auto capped_prefix=capped->state_prefix();
  cap=capped->full_time_population_view().rows.size();
  StepState st; Harness hc;
  row=refusal([&] { (void)hc.step(*capped,st); });
  row["name"]="changed exact carrier cap"; row["unchanged"]=capped->state_prefix()==capped_prefix;
  need(row["type"]=="resource" && row["unchanged"].get<bool>(), "changed cap",row); rows.push_back(row);
  cap=old_cap;
  return rows;
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 4) throw std::runtime_error("usage: test_birth_memo_diff OUT_DIR BULK_CARRIERS STEPS");
    const std::filesystem::path out=argv[1]; std::filesystem::create_directories(out);
    carrier_count=std::stoull(argv[2]); const int steps=std::stoi(argv[3]);
    if (carrier_count>50000 || steps<1 || steps>64) throw std::runtime_error("bounded fixture size");
    cap=carrier_count+2048;
    auto a=fresh_v8(); Harness h; StepState st;
    if (std::string(m::PicArena::source_identity())!=PIC_REACTION_NATIVE_SOURCE_SHA256)
      throw std::runtime_error("linked archive differs from fixture source identity");
    J report={{"source_identity",m::PicArena::source_identity()},{"bulk_carriers",carrier_count},{"steps",steps}};
    std::ofstream bytes(out/"physical.bin",std::ios::binary);
    J seals=J::array(), physical_hashes=J::array();
    double step_seconds=0; std::uint64_t births=0,sinks=0,newborn_sink_steps=0;
    for (int k=0;k<steps;++k) {
      h.place=k%3; std::shared_ptr<ReferenceMotion> motion;
      auto timer=std::chrono::steady_clock::now();
      const auto r=h.step(*a,st,reference_hook(cap,{},true,true,&motion));
      step_seconds+=since(timer);
      auto b=physical(*a,r,*motion,st);
      bytes.write(b.b.data(),std::streamsize(b.b.size())); physical_hashes.push_back(b.sha());
      births+=r.flows->births->entries.size();
      bool newborn=false;
      if (r.flows->sinks) { sinks+=r.flows->sinks->carriers; for (const auto &s:r.flows->sinks->records) newborn|=s.native_id>=r.birth_base; }
      newborn_sink_steps+=newborn;
      const auto history = a->full_time_history_status();
      seals.push_back({{"accepted",a->full_time_boundary()->accepted_arena_sha256},{"prefix",a->state_prefix()},
          {"anchor_root",history.anchor_root_sha256},{"anchor_head",history.anchor_head_sha256},
          {"anchor_chain",sha(a->export_full_time_anchor_chain())},
          {"flows",a->full_time_flow_totals().flows_sha256},{"motion",motion->record.record_sha256},
          {"birth_row",r.flows->births->row_sha256},{"birth_root",r.flows->births->births_root_sha256},
          {"sink_ledger",r.flows->sinks ? r.flows->sinks->ledger_sha256 : ""},
          {"sink_records",r.flows->sinks ? r.flows->sinks->records_root_sha256 : ""}});
    }
    bytes.close(); a->validate_full_time_history();
    report["birth_prepare_seconds"]=birth_seconds; report["full_step_seconds"]=step_seconds;
    report["births"]=births;report["sunk"]=sinks;report["newborn_sink_steps"]=newborn_sink_steps;
    report["final_carriers"]=a->full_time_population_view().rows.size();
    need(births>std::uint64_t(steps) && sinks>0 && newborn_sink_steps>0,"populated births and sinks");
    report["physical_step_sha256"]=physical_hashes; report["release_provenance_seals"]=seals;
    // Refusals use the same populated base, after valid memo proof re-use above.
    report["refusals"]=refused_controls();
    report["checks"]=n_checks;report["failures"]=n_fail;
    report["notes"]="Manufactured zero-field motion/zero-event chemistry fixture; no production runtime forecast. Exact physical bytes are separate from provenance-bound release seals.";
    std::ofstream(out/"report.json")<<report.dump(2)<<"\n";
    std::cout<<n_checks<<" differential controls, "<<n_fail<<" failures\n";
    return n_fail ? 1 : 0;
  } catch (const std::exception &e) { std::cerr<<"error: "<<e.what()<<"\n"; return 2; }
}
