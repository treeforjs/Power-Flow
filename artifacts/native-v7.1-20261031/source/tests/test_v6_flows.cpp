// Native v6 controls ported to native v7: births (arrival owner v1 plans,
// contract (a)-(f)), first-hit sinks with exact sink ledgers, births and sinks
// in one step, motion record v7 admission (genuine Hall corrective v2 receipts
// through the transitional v7 adapter, and the zero-field reference flows
// record), the native kinetic first-hit port against the sealed Hall helper,
// bounded history with births and sinks over >= 1e4 steps, checkpoint/arrival
// restart alignment and the v6 common predicate. v7 identity changes:
// API-007.txt section 5.
// usage: test_v6_flows OUTPUT_DIRECTORY [reduced]
#include "device_mode.hpp"
#include "v6_fixture.hpp"
#include "../source/src/full_time_contract_validation.hpp"
#include "arrival_v1_births.hpp"
#include "common_commit_v6.hpp"
#include "hall-finite-foil-em-20261029/axial_mesh.hpp"
#include "motion_receipt.hpp" // Hall corrective v2
#include "hall_receipt_v7.hpp"
#include <random>

using namespace v6;
namespace pa = pf_arrival;
namespace av = m::arrival_v1;
using FK = m::FullTimeFlowRefusal::Kind;
using RK = m::FullTimeResourceRefusal::Kind;
int forbidden_calls = 0;
extern "C" void
__wrap__ZN6raikou12mfem_backend6GapPic12solve_fieldsEv(m::GapPic *) {
  ++forbidden_calls;
  throw std::runtime_error("forbidden native solver trap");
}
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic7advanceEd(m::GapPic *,
                                                                 double) {
  ++forbidden_calls;
  throw std::runtime_error("forbidden native advance trap");
}
J measurements = J::object();
bool reduced = false;
std::uint64_t long_steps = 0; // "long=N": only the long run, N steps

const axial::Mesh &mesh() {
  static const axial::Mesh g;
  return g;
}
m::FullTimeGenesisBinding hall_genesis() {
  auto g = genesis();
  g.mesh_sha256 = hall_driver::motion_mesh_sha256(mesh(), options());
  return g;
}
std::unique_ptr<m::PicArena>
fresh(std::vector<m::FullTimeAuthoredPacket> pop = base_population(),
      const m::FullTimeConfiguration &c = {}) {
  return m::PicArena::create_full_time(options(), states(), channels(),
                                       std::move(pop), hall_genesis(), c);
}
std::uint64_t next_id(const m::PicArena &a) {
  return a.full_time_flow_totals().next_native_id;
}
// ---------------------------------------------------------------------------
// Genuine Hall corrective v2 provider at identically zero total field.
// v7: both former routes (v2 JSON seal, binary v3) enter through the v7
// adapter (the record is the motion record v7).
enum class HallRoute { v2_json, v3_binary };
const std::vector<double> &zero_E() {
  static const std::vector<double> E(mesh().ne, 0.0);
  return E;
}
const std::vector<double> &zero_B() {
  static const std::vector<double> B(mesh().nb, 0.0);
  return B;
}
hall_driver::MotionReceipt hall_receipt(m::PrivateTransportTrial &trial,
                                        const m::FullTimeStepBinding &b,
                                        std::size_t cap) {
  const auto census = trial.native_census(cap);
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
  mb.species_sha256 = sha("v6-species");
  mb.common_clock_sha256 = b.common_clock_sha256;
  mb.kinetic_policy_sha256 =
      hall_driver::kinetic_policy_sha256(hall_driver::declared_kinetic_policy(options()));
  return hall_driver::prepare_canonical_motion(mesh(), options(), census, mb, {}, E, B);
}
Hooks::Motion hall_hook(HallRoute route, std::size_t cap) {
  return [=](m::PrivateTransportTrial &trial, m::FullTimeMotionBinding &binding) {
    (void)route;
    auto receipt = std::make_shared<hall_driver::MotionReceipt>(hall_receipt(trial, binding.common, cap));
    auto v7 = m::v7_adapter::from_hall_receipt(*receipt, trial, zero_E(), zero_B());
    binding.transport_trial_sha256 = v7->record.record_sha256;
    auto storage = std::make_shared<std::pair<std::shared_ptr<hall_driver::MotionReceipt>,
                                              std::shared_ptr<m::v7_adapter::MotionV7>>>(receipt, v7);
    return trial.prepare_total_field_motion(binding, m::bind_total_field_motion_v7(v7->record, storage));
  };
}
// Test copy of the adapter mapping (for mutation controls).
struct OwnedRecord {
  std::vector<m::GapParticle> endpoints;
  std::vector<m::FullTimeMotionPath> paths;
  std::vector<m::FullTimeSinkRecord> sinks;
  std::vector<m::FullTimeBirthFlight> births;
  m::FullTimeMotionRecordV7 record;
  void bind() {
    record.endpoints = endpoints;
    record.paths = paths;
    record.sinks = sinks;
    record.births = births;
  }
};
std::shared_ptr<OwnedRecord> record_of(const hall_driver::MotionReceipt &r,
                                      const m::PrivateTransportTrial &trial) {
  auto o = std::make_shared<OwnedRecord>();
  const auto v7 = m::v7_adapter::from_hall_receipt(r, trial, zero_E(), zero_B());
  o->endpoints.assign(r.endpoints().begin(), r.endpoints().end());
  o->paths = v7->paths;
  o->sinks = v7->sinks;
  o->births = v7->births;
  o->record = v7->record;
  o->bind();
  return o;
}
// Open a step up to the transport fork (zero-event first exposure).
struct OpenStep {
  m::FullTimeStepBinding b;
  m::PrivateFullTimeStep outer;
  std::optional<m::CollisionScope> first;
};
OpenStep open_step(m::PicArena &a, StepState &st, std::size_t cap = 64,
                   std::size_t records = 32, double dt = 5e-15) {
  const auto k = a.full_time_boundary()->accepted_step;
  auto b = binding_of(a, dt, sha(step_text("v6-common-clock", k)));
  OpenStep s{b, a.begin_full_time_step(b, cap, records), std::nullopt};
  m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s,
                        .5 * b.dt_s, st.sequence, b.accepted_field_sha256};
  s.first.emplace(s.outer.first_exposure(e1));
  auto t = s.first->prepare_compound_with_completion(1, [&](m::CollisionScope &) {
    return completion(e1, st.rng, sha(step_text("v6-rng-first", k)));
  });
  if (!s.first->try_publish(t))
    throw std::runtime_error("first exposure");
  return s;
}
m::FullTimeMotionBinding motion_binding(const m::FullTimeStepBinding &b,
                                        const std::string &source_after) {
  return {b, "", b.accepted_field_sha256,
          sha(step_text("v6-field-after", b.accepted_step)),
          sha(step_text("v6-heat-after", b.accepted_step)), source_after};
}
// ---------------------------------------------------------------------------
// 1. Native kinetic first-hit port vs the sealed Hall helper.
void kinetic_port_controls() {
  const auto o = options();
  const auto policy = hall_driver::declared_kinetic_policy(o);
  std::mt19937_64 g(20261030);
  std::uniform_real_distribution<double> ux(1e-9, o.width_m - 1e-9),
      uy(1e-9, o.gap_m - 1e-9), uz(o.short_z_m + 1e-9, o.feed_z_m - 1e-9),
      ud(-1, 1), scale(-12, -1.5);
  std::size_t compared = 0, mismatched = 0;
  std::array<std::size_t, 6> hits{};
  const std::size_t n = reduced ? 40000 : 200000;
  for (std::size_t i = 0; i < n; ++i) {
    m::PicVector x{ux(g), uy(g), uz(g)};
    try {
      hall_driver::require_kinetic_interior(policy, x);
    } catch (const std::exception &) {
      continue;
    }
    const double s = std::pow(10.0, scale(g));
    m::PicVector D{s * ud(g), s * ud(g), s * ud(g)};
    if (i % 7 == 0)
      D[i % 3] = 0; // axis-parallel segments
    if (i % 11 == 0) { // land exactly on a plane (fraction 1 ties)
      const int axis = int(i % 2);
      D[axis] = (axis ? o.gap_m : o.width_m) - x[axis];
    }
    const auto h = hall_driver::first_kinetic_hit(policy, x, D);
    const auto n6 = m::full_time_first_kinetic_hit(o, x, D);
    ++compared;
    const bool same = int(h.sink) == n6.surface && h.fraction == n6.fraction &&
                      h.axis == n6.axis && h.plane_m == n6.plane_m;
    mismatched += !same;
    ++hits[std::size_t(int(h.sink) + 1)];
  }
  measurements["kinetic_port"] = {{"segments_compared", compared},
                                  {"mismatches", mismatched},
                                  {"none", hits[0]}, {"lower_foil", hits[1]},
                                  {"upper_foil", hits[2]}, {"outer_xy", hits[3]},
                                  {"short_bridge", hits[4]}, {"axial", hits[5]}};
  check(mismatched == 0 && compared > n / 2,
        "native kinetic first-hit port equals the sealed Hall first_kinetic_hit "
        "bitwise (surface, fraction, axis, plane)");
  check(hits[1] && hits[2] && hits[3] && hits[4] && hits[5] && hits[0],
        "port comparison covers no-hit and all five kinetic surfaces");
  // Foil-face start leaving into vacuum reports no hit on its own face.
  const auto leave = m::full_time_first_kinetic_hit(o, {.004, .007, 0}, {0, -1e-9, 0});
  const auto enter = m::full_time_first_kinetic_hit(o, {.004, .007 - 1e-9, 0}, {0, 2e-9, 0});
  check(leave.surface == -1 && enter.surface == 1 && enter.fraction > 0 &&
            enter.fraction < 1 && enter.plane_m == .007 && enter.axis == 1,
        "foil-face start leaving into vacuum has no hit; entering it hits the foil");
}
// ---------------------------------------------------------------------------
// 2. Motion record v7 (genuine Hall corrective v2 receipts through the v7
//    adapter): encoders, census v2, exact budget and tampering.
void plain_v3_controls() {
  auto a = fresh(), b = fresh();
  StepState sa, sb;
  StepPlan plan;
  full_step(*a, plan, sa, {nullptr, hall_hook(HallRoute::v2_json, 64), nullptr});
  full_step(*b, plan, sb, {nullptr, hall_hook(HallRoute::v3_binary, 64), nullptr});
  const auto pa_ = a->rate_population(64), pb = b->rate_population(64);
  bool same_rows = pa_.size() == pb.size();
  for (std::size_t i = 0; same_rows && i < pa_.size(); ++i)
    same_rows = pa_[i].packet_id == pb[i].packet_id &&
                std::memcmp(&pa_[i].native, &pb[i].native, sizeof(m::GapParticle)) == 0;
  check(same_rows && a->state_prefix() == b->state_prefix(),
        "v7: the same Hall motion gives identical native rows and identities");
  a->validate_full_time_history();
  b->validate_full_time_history();
  check(true, "histories with v7 motion rows revalidate");
  // Three independent encoders of the v7 record agree (native, test encoder
  // from API-007.txt, pf_v7 library header + sections).
  auto o = open_step(*b, sb);
  auto trial = o.outer.fork_transport_trial();
  auto receipt = hall_receipt(trial, o.b, 64);
  auto rec = record_of(receipt, trial);
  {
    pf_v7::MotionHeaderV7 h;
    const auto &x = rec->record;
    h.accepted_step = x.accepted_step;
    h.staged_arena_revision = x.staged_arena_revision;
    h.trial_stage = x.trial_stage;
    h.begin_s = x.begin_s;
    h.end_s = x.end_s;
    h.dt_s = x.dt_s;
    h.staged_arena_sha256 = x.staged_arena_sha256;
    h.trial_field_sha256 = x.trial_field_sha256;
    h.mesh_sha256 = x.mesh_sha256;
    h.material_sha256 = x.material_sha256;
    h.species_sha256 = x.species_sha256;
    h.common_clock_sha256 = x.common_clock_sha256;
    h.kinetic_policy_sha256 = x.kinetic_policy_sha256;
    h.budget = x.budget;
    h.census_sha256 = x.census_sha256;
    h.current_A_m_sha256 = x.current_A_m_sha256;
    h.endpoint_charge_C_sha256 = x.endpoint_charge_C_sha256;
    h.initial_charge_C_sha256 = x.initial_charge_C_sha256;
    h.kinetic_policy = x.kinetic_policy;
    h.native_to_field_m = x.native_to_field_m;
    const auto sec = pf_v7::motion_sections_v7(
        {reinterpret_cast<const pf_v7::Carrier *>(x.endpoints.data()), x.endpoints.size()},
        {reinterpret_cast<const pf_v7::Path *>(x.paths.data()), x.paths.size()}, {}, {});
    check(m::full_time_motion_record_v7_sha256(x) == reference_record_sha(x) &&
              reference_record_sha(x) == pf_v7::motion_record_v7_sha256(h, sec) &&
              x.record_sha256 == reference_record_sha(x) && x.record_sha256 != receipt.record_v3_sha256() &&
              x.record_sha256 != receipt.seal_sha256(),
          "native, test (API-007.txt) and library encoders agree on the v7 record digest");
  }
  const auto census = trial.native_census(64);
  const auto &c1 = trial.census_sha256();
  const auto &c2 = trial.census_sha256();
  {
    std::vector<std::string> states(trial.census_v7().states.begin(), trial.census_v7().states.end());
    std::vector<pf_v7::HostCensusEntry> entries;
    for (const auto &row : census)
      entries.push_back({row.packet_id, row.state, reinterpret_cast<const pf_v7::Carrier *>(&row.native)});
    check(&c1 == &c2 && c1 == trial.census_v7().sha256 &&
              c1 == pf_v7::census_v2_sha256(pf_v7::census_records(entries, states), states) &&
              c1 != hall_driver::canonical_census_sha256(census) && c1 == rec->record.census_sha256,
          "trial census v2 equals the rebuild from host texts, is cached, differs from v1");
  }
  const auto base = b->state_prefix();
  auto attempt = [&](std::function<void(OwnedRecord &)> mutate, bool reseal,
                     std::function<void(m::FullTimeMotionBinding &)> bind = {}) {
    auto x = record_of(receipt, trial);
    mutate(*x);
    x->bind();
    if (reseal)
      x->record.record_sha256 = m::full_time_motion_record_v7_sha256(x->record);
    auto mb = motion_binding(o.b, b->full_time_boundary()->source_sha256);
    mb.transport_trial_sha256 = x->record.record_sha256;
    if (bind)
      bind(mb);
    return [&trial, x, mb]() {
      (void)trial.prepare_total_field_motion(mb, m::bind_total_field_motion_v7(x->record, x));
    };
  };
  auto ulp = [](double v) { return std::nextafter(v, 1.0); };
  check(flow_refuses(attempt([&](auto &x) { x.endpoints[0].position_m[0] = ulp(x.endpoints[0].position_m[0]); }, false), FK::motion_record),
        "v7: one endpoint (one leaf) changed after sealing refuses (digest)");
  check(flow_refuses(attempt([&](auto &x) { x.paths[1].field_displacement_m[1] = ulp(x.paths[1].field_displacement_m[1]); }, false), FK::motion_record),
        "v7: one path (one leaf) changed after sealing refuses (digest)");
  check(flow_refuses(attempt([](auto &) {}, true, [](auto &mb) { mb.transport_trial_sha256 = sha("other"); }), FK::motion_record),
        "v7: binding transport identity differs from the record digest");
  check(flow_refuses(attempt([](auto &x) { x.record.census_sha256 = sha("foreign-census"); }, true), FK::motion_record),
        "v7: census digest differs from the trial census");
  check(flow_refuses(attempt([&](auto &x) { x.record.census_sha256 = hall_driver::canonical_census_sha256(census); }, true),
                     FK::motion_record),
        "v7: census digest in the old canonical-motion-census-v1 form refuses");
  check(flow_refuses(attempt([](auto &x) { x.record.native_to_field_m[0] *= 2; }, true), FK::motion_record),
        "v7: frame shift differs");
  check(flow_refuses(attempt([&](auto &x) { x.record.kinetic_policy[4] = ulp(x.record.kinetic_policy[4]); }, true), FK::motion_record),
        "v7: kinetic policy operand differs");
  check(flow_refuses(attempt([](auto &x) { std::swap(x.endpoints[0], x.endpoints[1]); }, true), FK::sink_duplicate),
        "v7: endpoint rows out of census order");
  check(flow_refuses(attempt([](auto &x) { x.endpoints[1].id = 999; }, true), FK::sink_duplicate),
        "v7: endpoint of a carrier outside the census");
  check(flow_refuses(attempt([](auto &x) { x.paths.pop_back(); }, true), FK::motion_record),
        "v7: missing path row (endpoint/path bijection)");
  check(flow_refuses(attempt([&](auto &x) { x.paths[0].field_start_m[2] = ulp(x.paths[0].field_start_m[2]); }, true), FK::motion_record),
        "v7: path start differs from the census position");
  check(flow_refuses(attempt([](auto &x) { x.paths[1].flight_s *= .5; }, true), FK::motion_record),
        "v7: resident carrier flight differs from dt");
  check(refuses(attempt([](auto &x) {
          x.endpoints[0].position_m[0] = options().width_m;
          x.paths[0].field_displacement_m[0] = (options().width_m - .0045) - x.paths[0].field_start_m[0];
        }, true)),
        "v7: endpoint on the outer plane without a sink record (strict interior)");
  check(refuses(attempt([](auto &x) { x.endpoints[1].weight *= 2; }, true)),
        "v7: endpoint altered carrier weight");
  check(flow_refuses(attempt([](auto &x) { x.record.budget[0] = 1e-30; }, true), FK::motion_record),
        "v7: kinetic budget differs from the native recomputation");
  {
    // One ulp of the kinetic sum (defect and work adjusted so only the
    // rounding form differs) refuses; so does the v6 x87 form when it differs.
    auto x = record_of(receipt, trial);
    const double K = x->record.budget[0];
    check(flow_refuses(attempt([&](auto &y) {
            y.record.budget[0] = std::nextafter(K, 1.0);
            y.record.budget[2] = y.record.budget[0] - y.record.budget[1];
          }, true), FK::motion_record),
          "v7: one-ulp change of kinetic_change_J (guard-consistent) refuses");
    measurements["v7_kinetic_vs_receipt_x87"] = {{"exact", K}, {"x87", receipt.budget().kinetic_change_J},
                                                 {"identical", K == receipt.budget().kinetic_change_J}};
  }
  check(flow_refuses(attempt([](auto &x) { x.record.staged_arena_sha256 = sha("other-prefix"); }, true), FK::motion_record),
        "v7: staged arena prefix differs");
  check(flow_refuses(attempt([](auto &x) { x.record.coverage = "all_native_carriers_no_sources_no_boundaries"; }, true), FK::motion_record),
        "v7: coverage text differs (the v6 plain coverage refuses)");
  check(flow_refuses(attempt([](auto &x) { x.record.schema = "canonical-full-time-total-field-motion-v3-binary"; }, false),
                     FK::motion_record),
        "v7: a v6 schema text refuses");
  check(flow_refuses(attempt([](auto &x) { x.record.neutral_stream_dt_max_s = 2e-14; }, true), FK::motion_record),
        "v7: neutral stream bound differs from the native configuration");
  check(flow_refuses(attempt([](auto &x) { x.record.native_DKD_called = true; }, true), FK::motion_record),
        "v7: native D-K-D flag set");
  check(refuses(attempt([](auto &x) {
          m::FullTimeSinkRecord k;
          k.native_id = x.endpoints[0].id;
          x.sinks.push_back(k);
          x.endpoints.erase(x.endpoints.begin());
        }, false)),
        "v7: unsealed sink insertion refuses");
  check(b->state_prefix() == base, "every refused v7 admission left the accepted state unchanged");
  auto mb = motion_binding(o.b, b->full_time_boundary()->source_sha256);
  auto x = record_of(receipt, trial);
  mb.transport_trial_sha256 = x->record.record_sha256;
  auto ok = trial.prepare_total_field_motion(mb, m::bind_total_field_motion_v7(x->record, x));
  check(ok.audit().at("record_schema") == m::motion_v7_schema &&
            ok.audit().at("record_sha256") == x->record.record_sha256,
        "the same trial admits the unmutated genuine v7 record after refusals");
}
// ---------------------------------------------------------------------------
// Arrival fixture.
std::shared_ptr<const pa::Policy> arrival_policy(double we = 4, double wn = 2,
                                                 std::uint64_t seed = 20261030) {
  pa::Configuration c;
  c.run_seed = seed;
  c.streams = {{0, pa::Kind::electron, 0, "e", we, 4096}, {1, pa::Kind::neutral, 0, "n", wn, 4096}};
  return std::make_shared<const pa::Policy>(c);
}
pa::StepInput arrival_input(std::uint64_t step, const m::FullTimeStepBinding &b,
                            double electrons, double neutral_debit) {
  pa::StepInput in;
  in.step = step;
  if (electrons > 0)
    in.electrons = {{0, electrons, std::log(electrons), true, b.begin_s, b.dt_s, b.end_s}};
  if (neutral_debit > 0)
    in.neutrals = {{1, pa::Count(neutral_debit), std::log(neutral_debit), pa::Count(1e20), b.begin_s, b.end_s}};
  return in;
}
// Birth placement: kind 0 = emitting-face births (electrons and neutrals on
// the upper foil lower face, moving into the gap), 1 = sink-prone births
// (electrons near the outer x face moving +x, neutrals just above the lower
// foil moving -y), 2 = interior births in the gap.
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
std::string source_seal(std::uint64_t k) { return sha(step_text("source", k)); }
// One step with arrival births through the native birth scope (committing
// the arrival owner with the native receipt after publication).
struct ArrivalRun {
  std::shared_ptr<const pa::Policy> P;
  std::unique_ptr<pa::Owner> owner;
  int place = 1;
  double electrons = 6, neutral_debit = 3;
  std::optional<pa::Staged> staged;
  std::uint64_t committed = 0; // records are verified, not retained
  std::uint64_t births = 0, source_index = 0;
  explicit ArrivalRun(std::shared_ptr<const pa::Policy> p, std::uint64_t first = 0)
      : P(std::move(p)), owner(std::make_unique<pa::Owner>(P, first)) {}
  Hooks::Motion motion = reference_hook(256);
  Hooks hooks() {
    Hooks h;
    h.births = [this](m::BirthScope &scope) -> std::optional<m::FullTimeBirthCredit> {
      const auto step = scope.accepted_step();
      m::FullTimeStepBinding b;
      b.begin_s = scope_begin;
      b.dt_s = scope_dt;
      b.end_s = scope_begin + scope_dt;
      auto plan = owner->plan(arrival_input(step, b, electrons, neutral_debit));
      staged.emplace(owner->stage(std::move(plan), scope.base_native_id()));
      bind_neutral_times(*staged);
      const auto kin = placement(*staged, place);
      return av::credit(*P, *staged, kin, av::plan_debits(staged->plan()),
                        scope.source_before_sha256(), source_seal(++source_index));
    };
    h.motion = motion;
    return h;
  }
  double scope_begin = 0, scope_dt = 0;
  StepResult step(m::PicArena &a, StepState &st, StepPlan plan = {}) {
    plan.maximum_carriers = std::max<std::size_t>(plan.maximum_carriers, 256);
    scope_begin = a.full_time_boundary()->time_s;
    scope_dt = plan.dt_s;
    auto r = full_step(a, plan, st, hooks());
    if (!r.flows || !r.flows->births)
      throw std::runtime_error("arrival step without native birth receipt");
    const auto record = owner->commit(std::move(*staged), av::receipt(*r.flows->births));
    if (record.chain_head != owner->ledger_head())
      throw std::runtime_error("arrival ledger head");
    ++committed;
    staged.reset();
    births += r.flows->births->entries.size();
    return r;
  }
};
// ---------------------------------------------------------------------------
// 3. Births: canonical credit, every refusal in (b), (a) B_s, (c), (d), (e).
void birth_controls(const std::string &dir) {
  auto P = arrival_policy();
  auto a = fresh();
  StepState st;
  const auto base_prefix = a->state_prefix();
  std::optional<m::FullTimeBirthCredit> replayed;
  {
    auto o = open_step(*a, st);
    auto scope = o.outer.birth_scope();
    const auto B = scope.base_native_id();
    check(B == next_id(*a) && B == 3 && scope.accepted_step() == 0,
          "(a) B_s is native's next ID after the first exposure, read before source preparation");
    check(flow_refuses([&] { (void)o.outer.birth_scope(); }, FK::birth_scope),
          "birth scope is reserved once per outer step");
    pa::Owner owner(P, 0);
    auto staged = owner.stage(owner.plan(arrival_input(0, o.b, 9.5, 7)), B);
    bind_neutral_times(staged);
    const auto n = staged.births().size();
    check(n >= 4, "arrival plan has electron and neutral births");
    const auto kin = placement(staged, 0);
    const auto good = av::credit(*P, staged, kin, av::plan_debits(staged.plan()),
                                 scope.source_before_sha256(), source_seal(1));
    auto p1 = scope.prepare(good), p2 = scope.prepare(good);
    check(p1.receipt().row_sha256 == p2.receipt().row_sha256 &&
              p1.receipt().base == B && p2.receipt().base == B &&
              p1.receipt().entries.size() == n,
          "(a) re-evaluated credits from the same base give identical B_s and receipts");
    bool canonical = p1.receipt().next_id_after == B + n;
    for (std::size_t r = 0; r < n; ++r) {
      const auto &e = p1.receipt().entries[r];
      const auto &b = staged.births()[r];
      canonical = canonical && e.native_id == B + r && e.native_id == b.native_id &&
                  e.packet_id == b.packet_id() &&
                  std::bit_cast<std::uint64_t>(e.weight) == std::bit_cast<std::uint64_t>(b.weight);
    }
    check(canonical, "(b) one packet per leaf, native IDs B_s + rank, canonical arrival IDs");
    bool exact = true;
    for (const auto &s : p1.receipt().streams)
      for (const auto &d : staged.plan().decisions())
        if (d.stream_id == s.stream_id)
          exact = exact && s.source_debit == d.debit;
    check(exact && p1.receipt().streams.size() == 2,
          "(d) one exact ExactParticleCount equality per stream: native credit = source debit");
    // Every refusal of the credit (each leaves the scope usable).
    auto bad = [&](std::function<void(m::FullTimeBirthCredit &)> f) {
      auto c = good;
      f(c);
      return [&scope, c]() { (void)scope.prepare(c); };
    };
    check(flow_refuses(bad([](auto &c) { std::swap(c.births[0], c.births[1]); }), FK::birth_order),
          "(b) refuses births out of canonical order");
    check(flow_refuses(bad([](auto &c) { std::reverse(c.births.begin(), c.births.end()); }), FK::birth_order),
          "(b) refuses reversed order");
    check(flow_refuses(bad([](auto &c) { c.births.insert(c.births.begin() + 1, c.births[0]); }), FK::birth_order),
          "(b) refuses a duplicated birth");
    check(flow_refuses(bad([](auto &c) { c.births.erase(c.births.begin()); }), FK::birth_order),
          "(b) refuses a missing birth (indices 0..n-1)");
    check(flow_refuses(bad([](auto &c) { c.births[0].packet_id = "arrival-v1/s0/p0/k00"; }), FK::birth_identity),
          "(b) refuses a non-canonical (foreign) packet ID spelling");
    check(flow_refuses(bad([](auto &c) { c.births[0].packet_id = "electron"; }), FK::birth_identity),
          "(b) refuses a foreign packet ID outside the namespace (existing packet)");
    check(flow_refuses(bad([](auto &c) { c.births[1].native_id += 1; }), FK::birth_identity),
          "(b) refuses native ID != B_s + rank");
    check(flow_refuses(bad([](auto &c) {
            for (auto &b : c.births) ++b.native_id;
            ++c.base_native_id;
          }), FK::birth_binding),
          "(a) refuses a block staged against another base");
    check(flow_refuses(bad([](auto &c) {
            for (auto &b : c.births) {
              b.step = 1;
              b.packet_id = m::full_time_birth_packet_id(1, b.stream_id, b.index);
            }
            c.arrival_step = 1;
          }), FK::birth_binding),
          "(b) refuses births of another accepted step (replayed or future IDs)");
    check(flow_refuses(bad([](auto &c) { c.births[0].origin.clear(); }), FK::birth_identity),
          "refuses a birth without its source event origin");
    check(flow_refuses(bad([](auto &c) { c.streams[0].source_debit = c.streams[0].source_debit.credited(std::ldexp(1.0, -1074)); }), FK::birth_count),
          "(d) refuses native credit != source debit by one unit 2^-1074");
    check(flow_refuses(bad([](auto &c) { c.births[0].weight = -c.births[0].weight; }), FK::birth_count),
          "refuses a nonpositive birth weight");
    check(flow_refuses(bad([](auto &c) { c.streams.erase(c.streams.begin()); }), FK::birth_binding),
          "refuses births of an undeclared stream");
    check(flow_refuses(bad([](auto &c) { c.streams[1].state = "e"; }), FK::birth_binding),
          "refuses a stream/birth species mismatch");
    check(flow_refuses(bad([](auto &c) { c.plan_sha256 = "x"; }), FK::birth_binding),
          "refuses a missing plan digest");
    check(flow_refuses(bad([](auto &c) { c.source_before_sha256 = sha("other-source"); }), FK::birth_binding),
          "refuses a source transition not starting at the accepted source seal");
    check(flow_refuses(bad([](auto &c) { c.births[0].position_m = {.004, .0071, 0}; }), FK::birth_kinematics),
          "(c) refuses a birth inside metal");
    check(flow_refuses(bad([](auto &c) { c.births[0].proper_velocity_m_s = {0, 2e6, 0}; }), FK::birth_kinematics),
          "(c) refuses a foil-face birth moving into the metal");
    check(flow_refuses(bad([](auto &c) { c.births[0].position_m = {.002, .007, 0}; }), FK::birth_kinematics),
          "(c) refuses a birth on a foil corner/edge");
    check(flow_refuses(bad([](auto &c) { c.births[0].birth_time_s -= 1e-15; }), FK::birth_kinematics),
          "(c) refuses remaining flight != end - tau");
    check(flow_refuses(bad([](auto &c) {
            c.births[0].birth_time_s = c.births[0].birth_time_s + 1e-14;
            c.births[0].remaining_flight_s = 0;
          }), FK::birth_kinematics),
          "(c) refuses a birth time outside the step");
    check(flow_refuses([] { m::full_time_detail::birth_id_block_guard(~0ull - 5, 10); }, FK::birth_exhaustion) &&
              flow_refuses([] { m::full_time_detail::birth_id_block_guard(0, 1); }, FK::birth_exhaustion) &&
              !refuses([] { m::full_time_detail::birth_id_block_guard(~0ull - 10, 10); }),
          "(b) native ID exhaustion refuses (B_s + n must stay below 2^64)");
    check(a->state_prefix() == base_prefix, "every refused credit left the accepted state unchanged");
    // (c) birth time and remaining flight: newborns move only over [tau, end].
    auto trial = o.outer.fork_transport_trial(p1);
    const auto flights = trial.birth_flights();
    bool flights_ok = flights.size() == n;
    for (std::size_t r = 0; flights_ok && r < n; ++r)
      flights_ok = flights[r].native_id == B + r &&
                   flights[r].birth_time_s == staged.births()[r].time_s &&
                   flights[r].flight_s == staged.births()[r].remaining_s &&
                   flights[r].flight_s == o.b.end_s - flights[r].birth_time_s;
    check(flights_ok && trial.source_after_sha256() == source_seal(1) &&
              trial.native_census(64).size() == 2 + n,
          "(c) trial census includes newborns with flight end - tau and the post-birth source seal");
    auto mb = motion_binding(o.b, source_seal(1));
    auto try_motion = [&](std::function<void(ReferenceMotion &)> mutate, bool reseal = true) {
      auto hook = reference_hook(64, mutate, reseal);
      auto copy = mb;
      return [hook, copy, &trial]() mutable { (void)hook(trial, copy); };
    };
    check(flow_refuses(try_motion([](auto &r) { r.paths.back().flight_s = r.record.dt_s; }), FK::motion_record),
          "(c) newborn moved over the full step refuses (flight != end - tau)");
    check(flow_refuses(try_motion([](auto &r) { r.births.clear(); }), FK::motion_record),
          "record without the trial's birth flights refuses");
    check(flow_refuses(try_motion([](auto &r) { r.paths.back().field_start_m[0] += 1e-9; }), FK::motion_record),
          "newborn path not starting at its birth position refuses");
    check(flow_refuses([&] {
            auto r = reference_motion(trial, o.b, options(), 64, trial.census_sha256(), false);
            r->births.clear();
            r->seal();
            auto copy = mb;
            copy.transport_trial_sha256 = r->record.record_sha256;
            (void)trial.prepare_total_field_motion(copy, admit(r));
          }, FK::motion_record),
          "sealed record without the newborn flights cannot admit a births trial");
    check(refuses([&] {
            auto copy = mb;
            copy.source_after_sha256 = a->full_time_boundary()->source_sha256;
            (void)reference_hook(64)(trial, copy);
          }),
          "motion binding must carry the post-birth source seal");
    {
      std::shared_ptr<ReferenceMotion> kept;
      auto hook = reference_hook(64, {}, true, true, &kept);
      auto copy = mb;
      auto moved = hook(trial, copy);
      bool partial = kept->births.size() == n;
      for (std::size_t r = 0; partial && r < kept->paths.size(); ++r)
        if (kept->paths[r].native_id >= B)
          partial = kept->paths[r].flight_s == flights[kept->paths[r].native_id - B].flight_s;
      check(partial && kept->transcribed_hits == n,
            "(c) admitted newborn paths cover exactly [tau, end] from the emitting face");
      o.outer.accept_transport_trial(std::move(moved));
    }
    const auto k = o.b.accepted_step;
    m::ExposureBinding e2{o.b, m::FullTimeSlot::after_transport, o.b.end_s, .5 * o.b.dt_s,
                          st.sequence + 1, sha(step_text("v6-field-after", k))};
    auto s2 = o.outer.second_exposure(e2);
    auto t2 = s2.prepare_compound_with_completion(1, [&](m::CollisionScope &) {
      return completion(e2, sha(step_text("v6-rng-first", k)), sha(step_text("v6-rng-second", k)));
    });
    check(s2.try_publish(t2), "second exposure after births publishes");
    auto token = o.outer.seal({sha(step_text("v6-field-after", k)), sha(step_text("v6-heat-after", k)),
                               source_seal(1), sha(step_text("v6-rng-second", k))});
    const auto flows = token.candidate_full_time_flows();
    check(flows && flows->births && !flows->sinks &&
              flows->births->row_sha256 == p1.receipt().row_sha256 &&
              flows->births->staged_token_sha256 == staged.token() &&
              flows->births->plan_sha256 == staged.plan().sha256() &&
              token.candidate_full_time_boundary()->source_sha256 == source_seal(1),
          "(e) sealed step carries the birth receipt bound to the staged token and plan digest");
    check(a->try_publish(token), "outer step with births publishes");
    st.sequence += 2;
    st.rng = sha(step_text("v6-rng-second", k));
    const auto rec = owner.commit(std::move(staged), av::receipt(*flows->births));
    check(rec.chain_head == owner.ledger_head() && owner.next_ordinal() == n,
          "(e) arrival Owner::commit accepts the native receipt");
    const auto totals = a->full_time_flow_totals();
    check(totals.births == n && totals.birth_rows == 1 && totals.has_last_birth &&
              totals.last_birth_base == B && totals.last_birth_count == n &&
              totals.last_birth_token_sha256 == flows->births->staged_token_sha256 &&
              totals.next_native_id == B + n && a->snapshot()->ever_packet_ids.size() == 2,
          "accepted flows: last birth block, next ID; arrival IDs are not retained in the ever-ID set");
    a->save(dir + "/births-step1.json");
    {
      const auto doc = J::parse(std::ifstream(dir + "/births-step1.json"));
      std::size_t births_rows = 0;
      std::string row_sha;
      for (const auto &r : doc["payload"]["metadata"]["history"]["window"])
        if (r["kind"] == "full_time_birth") {
          ++births_rows;
          row_sha = r["sha256"];
        }
      check(births_rows == 1 && row_sha == flows->births->row_sha256 &&
                doc["payload"]["metadata"].contains("flows"),
            "(e) journal window holds one birth row whose digest is the receipt row digest");
    }
    a->validate_full_time_history();
    check(true, "history with a birth row revalidates (grammar, chain, source transition)");
    replayed = good;
  }
  // Replay of step-0 births in step 1 refuses; the namespace is reserved.
  {
    auto o2 = open_step(*a, st);
    auto sc = o2.outer.birth_scope();
    auto replay = *replayed;
    replay.base_native_id = sc.base_native_id();
    replay.source_before_sha256 = sc.source_before_sha256();
    check(flow_refuses([&] { (void)sc.prepare(replay); }, FK::birth_binding),
          "(b) replayed births of an earlier step refuse in the next step");
  }
  check(flow_refuses([&] {
          auto pop = base_population();
          pop.push_back({"arrival-v1/s0/p0/k0", "e", "manufactured", {4, 0}, {.004, .005, 0}, {0, 0, 0}});
          (void)fresh(pop);
        }, FK::birth_namespace),
        "(b) genesis cannot author a packet in the arrival-v1/ namespace");
  for (int which = 0; which < 2; ++which) {
    // A fresh outer: first exposure with one forced event whose product or
    // untouched remainder uses the reserved namespace.
    auto y = fresh();
    StepState sy;
    auto by = binding_of(*y, 5e-15, sha("ns-clock"));
    auto outer = y->begin_full_time_step(by, 64, 32);
    m::ExposureBinding ey{by, m::FullTimeSlot::before_transport, by.begin_s, .5 * by.dt_s, 0,
                          by.accepted_field_sha256};
    auto scope_y = outer.first_exposure(ey);
    const bool typed = flow_refuses([&] {
      (void)scope_y.prepare_compound_with_completion(2, [&](m::CollisionScope &s) {
        m::ArenaEvent e;
        e.id = "manufactured/ns";
        e.channel = channels()[0].id;
        e.primary = "electron";
        e.target = "neutral";
        e.event_weight = which ? 2 : 4;
        e.phase = s.phase();
        const auto c = channels()[0];
        e.provenance = {m::ArenaProvenance::Kind::ConditionalPhysicalData, c.dataset_sha256,
                        c.kernel_manifest_sha256, c.kernel_binding_file_sha256, c.model,
                        "forced_kernel_integration_only", s.state_prefix(), 0};
        if (which == 0) {
          e.products = {{"arrival-v1/s0/p9/k0", "e", {0, 0, 1e6}, m::ArenaProduct::At::Primary},
                        {"ns/neutral", "n", {0, 0, 0}, m::ArenaProduct::At::Target}};
        } else {
          e.untouched_primary_id = "arrival-v1/s0/p9/k1";
          e.products = {{"ns/e", "e", {0, 0, 1e6}, m::ArenaProduct::At::Primary},
                        {"ns/neutral", "n", {0, 0, 0}, m::ArenaProduct::At::Target}};
        }
        auto child = s.prepare(e);
        (void)s.try_publish(child);
        return completion(ey, genesis().rng_sha256, sha("x"));
      });
    }, FK::birth_namespace);
    check(typed, which ? "(b) a reaction cannot name an untouched remainder in the arrival-v1/ namespace"
                       : "(b) a reaction product cannot use the arrival-v1/ namespace");
  }
  // Resource refusals before allocation.
  {
    auto x = fresh();
    StepState sx;
    auto ox = open_step(*x, sx, 3, 32);
    auto sc = ox.outer.birth_scope();
    pa::Owner ow(P, 0);
    auto sg = ow.stage(ow.plan(arrival_input(0, ox.b, 9.5, 7)), sc.base_native_id());
    bind_neutral_times(sg);
    auto c = av::credit(*P, sg, placement(sg, 0), av::plan_debits(sg.plan()), sc.source_before_sha256(), source_seal(1));
    check(resource_refuses([&] { (void)sc.prepare(c); }, RK::carrier_cap),
          "births above the outer carrier budget refuse (typed carrier_cap)");
  }
  {
    auto x = fresh();
    StepState sx;
    auto ox = open_step(*x, sx, 64, 7);
    auto sc = ox.outer.birth_scope();
    pa::Owner ow(P, 0);
    auto sg = ow.stage(ow.plan(arrival_input(0, ox.b, 9.5, 7)), sc.base_native_id());
    bind_neutral_times(sg);
    auto c = av::credit(*P, sg, placement(sg, 0), av::plan_debits(sg.plan()), sc.source_before_sha256(), source_seal(1));
    check(resource_refuses([&] { (void)sc.prepare(c); }, RK::step_records),
          "birth row beyond the declared outer record budget refuses (typed step_records)");
  }
  {
    auto x = fresh();
    StepState sx;
    auto ox = open_step(*x, sx);
    auto sc = ox.outer.birth_scope();
    auto y = fresh();
    StepState sy;
    auto oy = open_step(*y, sy);
    auto scy = oy.outer.birth_scope();
    pa::Owner ow(P, 0);
    auto sg = ow.stage(ow.plan(arrival_input(0, oy.b, 9.5, 7)), scy.base_native_id());
    bind_neutral_times(sg);
    auto py = scy.prepare(av::credit(*P, sg, placement(sg, 0), av::plan_debits(sg.plan()), scy.source_before_sha256(), source_seal(1)));
    check(flow_refuses([&] { (void)ox.outer.fork_transport_trial(py); }, FK::birth_scope),
          "prepared births of another outer cannot be forked");
    auto z = fresh();
    StepState sz;
    auto bz = binding_of(*z, 5e-15, sha("early"));
    auto outer = z->begin_full_time_step(bz, 64, 32);
    check(flow_refuses([&] { (void)outer.birth_scope(); }, FK::birth_scope),
          "birth scope before the first exposure refuses");
  }
  // Births with zero leaves still bind the plan (every step is committed).
  {
    auto x = fresh();
    StepState sx;
    ArrivalRun run(P);
    run.electrons = 0;
    run.neutral_debit = 0;
    auto r = run.step(*x, sx);
    check(r.flows && r.flows->births && r.flows->births->entries.empty() &&
              run.owner->next_step() == 1 && !x->full_time_flow_totals().has_last_birth &&
              x->full_time_flow_totals().birth_rows == 1,
          "an empty birth credit binds the plan and commits the arrival step");
  }
  // (newborn reaction) A second-slot reaction may use a newborn packet.
  {
    auto x = fresh();
    StepState sx;
    auto ox = open_step(*x, sx);
    auto sc = ox.outer.birth_scope();
    pa::Owner ow(P, 0);
    auto sg = ow.stage(ow.plan(arrival_input(0, ox.b, 1.2 * 4, 0)), sc.base_native_id());
    auto pb = sc.prepare(av::credit(*P, sg, placement(sg, 2), av::plan_debits(sg.plan()),
                                    sc.source_before_sha256(), source_seal(1)));
    const auto newborn = pb.receipt().entries.front().packet_id;
    auto trial = ox.outer.fork_transport_trial(pb);
    auto mbx = motion_binding(ox.b, source_seal(1));
    // Interior newborns: the genuine Hall v2 provider produces a full-step
    // record (flight dt for every carrier); native refuses it for births.
    check(flow_refuses([&] {
            auto copy = mbx;
            (void)hall_hook(HallRoute::v2_json, 64)(trial, copy);
          }, FK::motion_record),
          "Hall v2 motion (adapted) of a births trial refuses (no partial flights)");
    check(flow_refuses([&] {
            auto copy = mbx;
            (void)hall_hook(HallRoute::v3_binary, 64)(trial, copy);
          }, FK::motion_record),
          "Hall corrective v2 motion (adapted) of a births trial refuses (no newborn flights)");
    ox.outer.accept_transport_trial(reference_hook(64)(trial, mbx));
    const auto k = ox.b.accepted_step;
    m::ExposureBinding e2{ox.b, m::FullTimeSlot::after_transport, ox.b.end_s, .5 * ox.b.dt_s,
                          sx.sequence + 1, sha(step_text("v6-field-after", k))};
    auto s2 = ox.outer.second_exposure(e2);
    auto t2 = s2.prepare_compound_with_completion(1, [&](m::CollisionScope &s) {
      m::ArenaEvent e;
      e.id = "manufactured/newborn-elastic";
      e.channel = channels()[0].id;
      e.primary = newborn;
      e.target = "neutral";
      e.event_weight = 4;
      e.scattered_primary_id = newborn;
      e.phase = s.phase();
      const auto c = channels()[0];
      e.provenance = {m::ArenaProvenance::Kind::ConditionalPhysicalData, c.dataset_sha256,
                      c.kernel_manifest_sha256, c.kernel_binding_file_sha256, c.model,
                      "forced_kernel_integration_only", s.state_prefix(), 0};
      m::PicVector u{};
      for (const auto &row : s.rate_population(64))
        if (row.packet_id == newborn)
          u = row.native.proper_velocity_m_s;
      e.products = {{newborn, "e", u, m::ArenaProduct::At::Primary},
                    {"newborn-elastic/neutral", "n", {0, 0, 0}, m::ArenaProduct::At::Target}};
      auto child = s.prepare(e);
      if (!s.try_publish(child))
        throw std::runtime_error("newborn reaction");
      return completion(e2, sha(step_text("v6-rng-first", k)), sha(step_text("v6-rng-second", k)));
    });
    check(s2.try_publish(t2), "second-slot reaction with a newborn packet (scattered primary keeps its arrival ID)");
    auto tok = ox.outer.seal({sha(step_text("v6-field-after", k)), sha(step_text("v6-heat-after", k)),
                              source_seal(1), sha(step_text("v6-rng-second", k))});
    check(x->try_publish(tok) && x->snapshot()->packets.count(newborn),
          "step with births and a newborn reaction publishes");
    x->validate_full_time_history();
  }
}
// ---------------------------------------------------------------------------
// 4. Sinks of each kind with exact ledgers, and refusals.
std::vector<m::FullTimeAuthoredPacket> sink_population() {
  auto pop = base_population();
  const double d = 2e-7;
  // electrons u = 1e8 m/s (D ~ 4.7e-7 m per 5 fs) at d = 2e-7 m from a plane
  pop.push_back({"to-lower", "e", "manufactured/sink", {3, 0}, {.004, .002 + d, .001}, {0, -1e8, 0}});
  pop.push_back({"to-upper", "i", "manufactured/sink", {5, 0}, {.005, .007 - 1e-12, .002}, {0, 6e3, 0}});
  pop.push_back({"to-outer-x", "e", "manufactured/sink", {7, 0}, {.009 - d, .0045, .003}, {1e8, 0, 0}});
  pop.push_back({"to-outer-y", "n", "manufactured/sink", {1.5, 0}, {.0010, 1e-12, .004}, {0, -650, 0}});
  pop.push_back({"to-bridge", "e", "manufactured/sink", {2, 0}, {.0045, .0047, -.0125 + d}, {0, 0, -1e8}});
  pop.push_back({"to-feed", "e", "manufactured/sink", {6, 0}, {.0045, .0043, .0125 - d}, {0, 3e5, 1e8}});
  pop.push_back({"to-axial-short", "e", "manufactured/sink", {9, 0}, {.0010, .0045, -.0125 + d}, {0, 0, -1e8}});
  // two co-located carriers (count 1 + 2^-60), sinking together
  pop.push_back({"pair-to-x0", "e", "manufactured/sink", {1.0, std::ldexp(1.0, -60)}, {d, .0050, .005}, {-1e8, 0, 0}});
  pop.push_back({"stays", "e", "manufactured/sink", {4, 0}, {.0060, .0045, .006}, {1e6, 0, 0}});
  return pop;
}
void sink_controls() {
  auto a = fresh(sink_population());
  StepState st;
  const auto before = a->state_prefix();
  const auto population = a->rate_population(64);
  std::map<std::string, std::pair<std::string, double>> by_packet; // id -> state, weight
  // Refusals on one open step (fresh trials; nothing is consumed).
  auto o = open_step(*a, st);
  auto find_sink = [](ReferenceMotion &r, m::FullTimeSinkSurface s) -> m::FullTimeSinkRecord & {
    for (auto &k : r.sinks)
      if (k.surface == s)
        return k;
    throw std::runtime_error("fixture sink missing");
  };
  // fork budget is 4 per outer: run refusal batches on separate outers.
  auto batch = [&](std::vector<std::pair<std::function<void(ReferenceMotion &)>, std::pair<FK, std::string>>> cases) {
    for (auto &[mutate, expect] : cases) {
      auto x = fresh(sink_population());
      StepState sx;
      auto ox = open_step(*x, sx);
      auto trial = ox.outer.fork_transport_trial();
      auto mb = motion_binding(ox.b, x->full_time_boundary()->source_sha256);
      auto hook = reference_hook(64, mutate, true);
      check(flow_refuses([&] { (void)hook(trial, mb); }, expect.first), expect.second);
    }
  };
  using S = m::FullTimeSinkSurface;
  batch({
      {[&](ReferenceMotion &r) {
         auto &k = find_sink(r, S::lower_foil);
         k.fraction = 1.5;
         k.flight_s = r.record.dt_s * 1.5;
         k.discarded_flight_s = r.record.dt_s - k.flight_s;
       },
       {FK::sink_outside_step, "sink hit fraction beyond the carrier flight (hit outside the step) refuses"}},
      {[&](ReferenceMotion &r) {
         auto &k = find_sink(r, S::upper_foil);
         k.flight_s = 2 * r.record.dt_s;
       },
       {FK::sink_outside_step, "sink flight longer than the step refuses"}},
      {[&](ReferenceMotion &r) {
         // carrier "stays" (a survivor) reported as sunk on its own segment
         auto it = std::find_if(r.endpoints.begin(), r.endpoints.end(), [](const auto &e) {
           return e.weight == 4 && e.position_m[0] > .0059 && e.position_m[0] < .0061;
         });
         const auto &c = *std::find_if(r.census.begin(), r.census.end(),
                                       [&](const auto &x) { return x.id == it->id; });
         auto &path = *std::find_if(r.paths.begin(), r.paths.end(),
                                    [&](const auto &p) { return p.native_id == it->id; });
         m::FullTimeSinkRecord k;
         k.native_id = it->id;
         k.kind = m::FullTimeSinkKind::escaped_exterior;
         k.surface = S::outer_xy;
         k.axis = 0;
         k.plane_m = .009;
         k.fraction = 1;
         for (int d = 0; d < 3; ++d)
           k.segment_m[d] = it->position_m[d] - c.position_m[d];
         k.hit_point_m = it->position_m;
         k.hit_point_m[0] = .009;
         k.incident_proper_velocity_m_s = it->proper_velocity_m_s;
         k.flight_s = r.record.dt_s;
         k.discarded_flight_s = 0;
         path.field_displacement_m[0] = (k.hit_point_m[0] + r.record.native_to_field_m[0]) - path.field_start_m[0];
         r.sinks.push_back(k);
         std::sort(r.sinks.begin(), r.sinks.end(), [](auto &x, auto &y) { return x.native_id < y.native_id; });
         r.endpoints.erase(it);
       },
       {FK::sink_not_crossing, "sink record for a carrier whose segment crosses no boundary refuses"}},
      {[&](ReferenceMotion &r) {
         auto &k = find_sink(r, S::lower_foil);
         k.surface = S::upper_foil;
       },
       {FK::sink_not_crossing, "sink surface differing from the native first hit refuses"}},
  });
  batch({
      {[&](ReferenceMotion &r) { find_sink(r, S::outer_xy).kind = m::FullTimeSinkKind::collected_at_conductor; },
       {FK::sink_record, "sink kind inconsistent with its surface refuses"}},
      {[&](ReferenceMotion &r) {
         auto &k = find_sink(r, S::short_bridge);
         k.hit_point_m[0] = std::nextafter(k.hit_point_m[0], 1.0);
       },
       {FK::sink_record, "sink hit point one ulp off refuses"}},
      {[&](ReferenceMotion &r) {
         // duplicate removal: the same sink twice (one survivor dropped to keep counts)
         r.sinks.insert(r.sinks.begin() + 1, r.sinks.front());
         r.endpoints.pop_back();
       },
       {FK::sink_duplicate, "duplicate removal of one carrier refuses"}},
      {[&](ReferenceMotion &r) {
         auto k = r.sinks.back();
         k.native_id = 99999;
         r.sinks.push_back(k);
         r.endpoints.pop_back();
       },
       {FK::sink_duplicate, "sink of a carrier outside the census refuses"}},
  });
  batch({
      {[&](ReferenceMotion &r) {
         // pair-to-x0 has two co-located carriers: keep the second one as a
         // survivor at its census position (zero displacement).
         std::size_t at = r.sinks.size();
         for (std::size_t i = 0; i < r.sinks.size() && at == r.sinks.size(); ++i)
           if (r.sinks[i].surface == S::outer_xy && r.sinks[i].segment_m[0] < 0)
             at = i;
         const auto second = r.sinks.at(at + 1);
         r.sinks.erase(r.sinks.begin() + std::ptrdiff_t(at) + 1);
         const auto &c = *std::find_if(r.census.begin(), r.census.end(),
                                       [&](const auto &x) { return x.id == second.native_id; });
         for (auto &row : r.paths)
           if (row.native_id == second.native_id) {
             for (int d = 0; d < 3; ++d)
               row.field_displacement_m[d] = (c.position_m[d] + r.record.native_to_field_m[d]) - row.field_start_m[d];
             row.flight_s = r.record.dt_s;
           }
         r.endpoints.push_back(c);
         std::sort(r.endpoints.begin(), r.endpoints.end(), [](auto &x, auto &y) { return x.id < y.id; });
       },
       {FK::sink_partial_packet, "sink records removing only one carrier of a two-carrier packet refuse"}},
      {[&](ReferenceMotion &r) {
         for (auto &s : r.sinks)
           if (s.surface == S::outer_xy && s.segment_m[1] < 0)
             s.segment_m[1] *= 1.0000001; // the neutral "to-outer-y"
       },
       {FK::sink_not_crossing, "neutral sink segment differing from its straight drift refuses"}},
      {[&](ReferenceMotion &r) {
         auto &k = find_sink(r, S::upper_foil);
         for (auto &p : r.paths)
           if (p.native_id == k.native_id)
             p.flight_s = r.record.dt_s;
       },
       {FK::sink_record, "sink path flight differing from the hit flight refuses"}},
  });
  {
    auto x = fresh(sink_population());
    StepState sx;
    auto ox = open_step(*x, sx);
    auto trial = ox.outer.fork_transport_trial();
    auto mb = motion_binding(ox.b, x->full_time_boundary()->source_sha256);
    auto hook = reference_hook(64, [&](ReferenceMotion &r) {
      // drop the axial (feed) sink and report the carrier at its hit point
      auto it = std::find_if(r.sinks.begin(), r.sinks.end(), [](const auto &k) {
        return k.surface == S::axial && k.plane_m > 0;
      });
      const auto k = *it;
      auto e = *std::find_if(r.census.begin(), r.census.end(),
                             [&](const auto &c) { return c.id == k.native_id; });
      e.position_m = k.hit_point_m;
      r.sinks.erase(it);
      r.endpoints.push_back(e);
      std::sort(r.endpoints.begin(), r.endpoints.end(), [](auto &a, auto &b) { return a.id < b.id; });
    });
    const auto msg = message_of([&] { (void)hook(trial, mb); });
    check(msg.find("axial endpoint") != std::string::npos,
          "strict-interior rule stays: a carrier on a boundary plane without a sink record refuses");
  }

  check(a->state_prefix() == before, "sink refusals left the accepted state unchanged");
  // Positive step with every kind.
  std::shared_ptr<ReferenceMotion> kept;
  auto trial = o.outer.fork_transport_trial();
  auto mb = motion_binding(o.b, a->full_time_boundary()->source_sha256);
  auto moved = reference_hook(64, {}, true, true, &kept)(trial, mb);
  o.outer.accept_transport_trial(std::move(moved));
  const auto k = o.b.accepted_step;
  m::ExposureBinding e2{o.b, m::FullTimeSlot::after_transport, o.b.end_s, .5 * o.b.dt_s,
                        st.sequence + 1, sha(step_text("v6-field-after", k))};
  auto s2 = o.outer.second_exposure(e2);
  auto t2 = s2.prepare_compound_with_completion(1, [&](m::CollisionScope &) {
    return completion(e2, sha(step_text("v6-rng-first", k)), sha(step_text("v6-rng-second", k)));
  });
  check(s2.try_publish(t2), "second exposure after sinks");
  auto token = o.outer.seal({sha(step_text("v6-field-after", k)), sha(step_text("v6-heat-after", k)),
                             a->full_time_boundary()->source_sha256, sha(step_text("v6-rng-second", k))});
  const auto flows = token.candidate_full_time_flows();
  check(flows && flows->sinks && !flows->births, "sealed step carries the sink ledger");
  const auto &ledger = *flows->sinks;
  std::array<std::uint64_t, 5> carriers{};
  for (const auto &r : kept->sinks)
    ++carriers[std::size_t(r.surface)];
  check(carriers[0] == 1 && carriers[1] == 1 && carriers[2] == 4 && carriers[3] == 1 && carriers[4] == 2,
        "sinks of every surface: lower foil, upper foil, outer x/y, short bridge, axial (feed and non-bridge short)");
  // Independent exact ledger recomputation from the census.
  std::map<std::uint64_t, m::ArenaRateEntry> census;
  for (const auto &r : population)
    census.emplace(r.native.id, r);
  std::array<m::ArenaCount, 5> w, ke, pos, neg;
  for (const auto &r : kept->sinks) {
    const auto &c = census.at(r.native_id);
    auto at = c.native;
    at.proper_velocity_m_s = r.incident_proper_velocity_m_s;
    const auto s = std::size_t(r.surface);
    w[s] = w[s].credited(c.native.weight);
    ke[s] = ke[s].credited(c.native.weight * m::pic_kinetic_energy_J(at));
    const int q = c.state == "e" ? -1 : c.state == "i" ? 1 : 0;
    for (int j = 0; j < std::abs(q); ++j)
      (q > 0 ? pos[s] : neg[s]) = (q > 0 ? pos[s] : neg[s]).credited(c.native.weight);
  }
  bool exact = ledger.surfaces.size() == 5 && ledger.carriers == kept->sinks.size();
  for (std::size_t s = 0; exact && s < 5; ++s) {
    const auto &t = ledger.surfaces[s];
    exact = t.carriers == carriers[s] && t.weight == w[s] && t.kinetic_energy_J == ke[s] &&
            t.positive_charge_e == pos[s] && t.negative_charge_e == neg[s] &&
            t.kind == ((s == 0 || s == 1 || s == 3) ? m::FullTimeSinkKind::collected_at_conductor
                                                    : m::FullTimeSinkKind::escaped_exterior);
  }
  check(exact, "exact sink ledger per surface: carriers, weights, incident kinetic energy, "
               "positive/negative elementary charge equal an independent exact recomputation");
  check(ledger.removed_packets.size() == 7 + 1 &&
            std::count(ledger.removed_packets.begin(), ledger.removed_packets.end(), "pair-to-x0") == 1,
        "removed logical packets listed once each (two-carrier packet removed whole)");
  bool tallies = !ledger.tallies.empty();
  for (const auto &t : ledger.tallies)
    tallies = tallies && (t.state == "e" ? t.charge_e == -1 : t.state == "i" ? t.charge_e == 1 : t.charge_e == 0);
  check(tallies && m::canonical::is_digest(ledger.ledger_sha256) &&
            m::canonical::is_digest(ledger.records_root_sha256),
        "per (surface, species) tallies with charge numbers and ledger/record digests");
  check(a->try_publish(token), "step with sinks publishes");
  st.sequence += 2;
  st.rng = sha(step_text("v6-rng-second", k));
  const auto snap = a->snapshot();
  bool gone = true;
  for (const auto &id : ledger.removed_packets)
    gone = gone && !snap->packets.count(id) && snap->ever_packet_ids.count(id);
  check(gone && a->rate_population(64).size() == population.size() - kept->sinks.size(),
        "removed carriers leave the arena; their packet IDs remain retired (never reused)");
  const auto totals = a->full_time_flow_totals();
  bool cumulative = totals.sunk_carriers == kept->sinks.size() && totals.sink_rows == 1;
  for (std::size_t s = 0; s < 5; ++s)
    cumulative = cumulative && totals.surfaces[s].weight == w[s];
  check(cumulative, "accepted cumulative sink totals equal the step ledger");
  a->validate_full_time_history();
  check(true, "history with a sink row revalidates");
  measurements["sink_ledger_example"] = {{"carriers", ledger.carriers},
                                         {"removed_packets", ledger.removed_packets},
                                         {"ledger_sha256", ledger.ledger_sha256}};
}
// ---------------------------------------------------------------------------
// 5. Births and sinks in the same step.
void same_step_controls() {
  auto P = arrival_policy();
  auto a = fresh(sink_population());
  StepState st;
  ArrivalRun run(P);
  run.place = 1; // sink-prone births
  run.electrons = 4 * 3.3;
  run.neutral_debit = 5;
  std::uint64_t newborn_sunk = 0;
  auto r = run.step(*a, st);
  check(r.flows && r.flows->births && r.flows->sinks,
        "one step carries both a birth receipt and a sink ledger");
  const auto B = r.flows->births->base;
  for (const auto &k : r.flows->sinks->records)
    if (k.native_id >= B)
      ++newborn_sunk;
  bool partial = true;
  for (const auto &k : r.flows->sinks->records)
    if (k.native_id >= B)
      partial = partial && k.flight_s + k.discarded_flight_s < 5e-15 &&
                k.flight_s <= k.flight_s + k.discarded_flight_s;
  check(newborn_sunk > 0 && partial,
        "newborns sink within their remaining flight (flight + discarded = end - tau < dt)");
  a->validate_full_time_history();
  check(run.owner->next_ordinal() == a->full_time_flow_totals().births,
        "arrival ordinal equals native cumulative births");
}
// ---------------------------------------------------------------------------
// 6. Long run: >= 1e4 steps with births and sinks, bounded history/memory.
void long_run(const std::string &dir) {
  const std::uint64_t steps = long_steps ? long_steps : reduced ? 2000 : 10000;
  auto P = arrival_policy(4, 2, 99);
  m::FullTimeConfiguration c;
  c.history.anchor_interval_steps = 1024;
  auto a = fresh(base_population(), c);
  StepState st;
  ArrivalRun run(P);
  run.place = 1;
  run.electrons = 4 * 0.75;
  run.neutral_debit = 1.5;
  J samples = J::array();
  std::uint64_t births = 0, sunk = 0, max_population = 0, steps_with_both = 0;
  std::array<m::ArenaCount, 5> weights;
  const auto t0 = std::chrono::steady_clock::now();
  long rss_after_warmup = 0;
  for (std::uint64_t k = 0; k < steps; ++k) {
    auto r = run.step(*a, st);
    births += r.flows->births->entries.size();
    if (r.flows->sinks) {
      sunk += r.flows->sinks->carriers;
      for (std::size_t s = 0; s < 5; ++s)
        weights[s] = [&] {
          m::ArenaCount x = weights[s];
          for (double v : r.flows->sinks->surfaces[s].weight.positive_carriers())
            x = x.credited(v);
          return x;
        }();
      if (!r.flows->births->entries.empty())
        ++steps_with_both;
    }
    max_population = std::max<std::uint64_t>(max_population, a->snapshot()->packets.size());
    if (k == 1000)
      rss_after_warmup = status_kib("VmRSS:");
    if (k % (steps >= 100000 ? 10000 : 1000) == 999) {
      const auto h = a->full_time_history_status();
      samples.push_back({{"step", k + 1}, {"VmRSS_KiB", status_kib("VmRSS:")},
                         {"packets", a->snapshot()->packets.size()},
                         {"ever_packet_ids", a->snapshot()->ever_packet_ids.size()},
                         {"anchors", h.anchors}, {"window_rows", h.window_rows},
                         {"window_bytes", h.window_bytes}, {"births", births}, {"sunk", sunk}});
    }
  }
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  const auto h = a->full_time_history_status();
  const auto totals = a->full_time_flow_totals();
  a->validate_full_time_history();
  const auto chain = a->export_full_time_anchor_chain();
  const auto anchors = m::verify_full_time_anchor_chain(chain, h.anchor_root_sha256, h.anchor_head_sha256);
  std::uint64_t anchored_births = 0, anchored_sunk = 0, birth_rows = 0;
  for (const auto &x : anchors) {
    anchored_births += x.window_births;
    anchored_sunk += x.window_sunk_carriers;
    birth_rows += x.window_birth_rows;
  }
  bool weights_ok = true;
  for (std::size_t s = 0; s < 5; ++s)
    weights_ok = weights_ok && totals.surfaces[s].weight == weights[s];
  const long rss_end = status_kib("VmRSS:");
  measurements["long_run"] = {{"steps", steps}, {"seconds", seconds},
                              {"births", births}, {"sunk_carriers", sunk},
                              {"steps_with_births_and_sinks", steps_with_both},
                              {"maximum_packets", max_population},
                              {"anchors", h.anchors}, {"anchor_bytes", h.anchor_bytes},
                              {"window_rows", h.window_rows}, {"window_bytes", h.window_bytes},
                              {"cumulative_rows", h.cumulative_rows},
                              {"VmRSS_KiB_after_1000_steps", rss_after_warmup},
                              {"VmRSS_KiB_end", rss_end}, {"VmHWM_KiB", status_kib("VmHWM:")},
                              {"ever_packet_ids", a->snapshot()->ever_packet_ids.size()},
                              {"samples", samples}};
  check(totals.births == births && totals.sunk_carriers == sunk && births > steps / 2 &&
            sunk > steps / 2 && steps_with_both > 0 && weights_ok,
        "long run: cumulative native births/sinks equal the per-step receipts and ledgers");
  check(run.owner->next_ordinal() == births && run.owner->next_step() == steps,
        "long run: arrival owner committed every step with the native receipts");
  check(h.anchors == steps / 1024 && h.window_rows <= 7 * 1024 &&
            anchors.size() == h.anchors && anchored_births <= births && anchored_sunk <= sunk &&
            birth_rows == h.anchors * 1024,
        "long run: anchors every K steps verify; one birth row per anchored step");
  check(max_population < 64 && a->snapshot()->ever_packet_ids.size() == 2,
        "long run: population bounded (births balanced by sinks); births add no identity entries");
#ifndef __SANITIZE_ADDRESS__
  // Normal build only: the AddressSanitizer quarantine holds freed memory.
  check(rss_end - rss_after_warmup < 8 * 1024,
        "long run: resident memory bounded after warm-up (< 8 MiB growth)");
#endif
  (void)dir;
}
// ---------------------------------------------------------------------------
// 7. Checkpoint: native next ID + arrival checkpoint consistency (f).
void checkpoint_controls(const std::string &dir) {
  auto P = arrival_policy(4, 2, 7);
  // Uninterrupted reference: 20 steps from genesis.
  auto run_from = [&](m::PicArena &a, ArrivalRun &run, StepState &st, int n) {
    std::vector<std::string> prefixes;
    for (int i = 0; i < n; ++i) {
      run.step(a, st);
      prefixes.push_back(a.state_prefix());
    }
    return prefixes;
  };
  auto ref = fresh(sink_population());
  ArrivalRun rref(P);
  StepState sref;
  const auto expected = run_from(*ref, rref, sref, 20);
  // Genesis checkpoint of both owners, restored, then the same 20 steps.
  auto g = fresh(sink_population());
  ArrivalRun rg(P);
  g->save(dir + "/genesis-v6.json");
  const auto arrival_bytes = rg.owner->save();
  auto restored = m::PicArena::restore_full_time(dir + "/genesis-v6.json", options(), states(),
                                                 channels(), *g->full_time_boundary());
  ArrivalRun rr(P);
  rr.owner = pa::Owner::restore(P, arrival_bytes);
  av::require_restart_alignment(restored->full_time_flow_totals(), *rr.owner);
  StepState srr;
  const auto replay = run_from(*restored, rr, srr, 20);
  check(replay == expected && rr.owner->ledger_head() == rref.owner->ledger_head(),
        "(f) genesis checkpoint of native + arrival restored: 20 steps with births/sinks reproduce "
        "every accepted prefix and the arrival ledger head");
  // After 20 steps: native checkpoint carries next ID and last birth block.
  ref->save(dir + "/step20-v6.json");
  const auto from_file = m::PicArena::read_full_time_flow_totals(dir + "/step20-v6.json");
  const auto live = ref->full_time_flow_totals();
  check(from_file.next_native_id == live.next_native_id && from_file.births == live.births &&
            from_file.last_birth_base == live.last_birth_base &&
            from_file.last_birth_count == live.last_birth_count &&
            from_file.flows_sha256 == live.flows_sha256 && from_file.accepted_step == 20,
        "(f) native checkpoint records next native ID, last birth block and flows digest");
  const auto a20 = rref.owner->save();
  auto arrival20 = pa::Owner::restore(P, a20);
  av::require_restart_alignment(from_file, *arrival20);
  check(true, "(f) arrival checkpoint restored alongside aligns with the native checkpoint");
  auto stale = pa::Owner::restore(P, arrival_bytes);
  check(refuses([&] { av::require_restart_alignment(from_file, *stale); }),
        "(f) an arrival checkpoint of another step does not align (refused)");
  // Next B_s is read from native.
  StepState s21 = sref;
  auto o = open_step(*ref, s21, 256);
  auto scope = o.outer.birth_scope();
  check(scope.base_native_id() == from_file.next_native_id,
        "(f) B_s of the next step equals the checkpointed native next ID");
  // Nonempty native restore validates the window (birth/sink rows), then defers.
  bool deferred = false;
  try {
    (void)m::PicArena::restore_full_time(dir + "/step20-v6.json", options(), states(), channels(),
                                         *ref->full_time_boundary());
  } catch (const m::FullTimeResourceRefusal &e) {
    deferred = e.kind() == RK::restore_deferred;
  }
  check(deferred, "nonempty restore validates birth/sink rows in the exported window, then refuses (deferred)");
}
// ---------------------------------------------------------------------------
// 8. Anchor interval independence and tamper detection with births/sinks.
void anchor_controls() {
  auto P = arrival_policy(4, 2, 11);
  std::vector<std::string> reference;
  std::string chain3, root3, head3;
  for (std::uint64_t K : {1ull, 3ull, 1024ull}) {
    m::FullTimeConfiguration c;
    c.history.anchor_interval_steps = K;
    auto a = fresh(sink_population(), c);
    ArrivalRun run(P);
    StepState st;
    std::vector<std::string> prefixes;
    for (int i = 0; i < 30; ++i) {
      run.step(*a, st);
      prefixes.push_back(a->state_prefix());
    }
    if (K == 1)
      reference = prefixes;
    check(prefixes == reference, "accepted states with births/sinks identical for every anchor interval K");
    a->validate_full_time_history();
    const auto h = a->full_time_history_status();
    const auto chain = a->export_full_time_anchor_chain();
    const auto fields = m::verify_full_time_anchor_chain(chain, h.anchor_root_sha256, h.anchor_head_sha256);
    std::uint64_t births = 0, rows = 0;
    for (const auto &f : fields) {
      births += f.window_births;
      rows += f.window_birth_rows;
    }
    check(fields.size() == 30 / K && rows == fields.size() * K &&
              (K == 1024 || births > 0) &&
              (fields.empty() || m::canonical::is_digest(fields.back().flows_sha256)),
          "v6 anchors commit birth/sink window counters and the flows digest");
    if (K == 3) {
      chain3 = chain;
      root3 = h.anchor_root_sha256;
      head3 = h.anchor_head_sha256;
    }
  }
  std::size_t positions = 0, refused = 0;
  for (std::size_t at = 0; at < chain3.size(); at += 53, ++positions) {
    auto x = chain3;
    x[at] = char(x[at] ^ 1);
    refused += refuses([&] { m::verify_full_time_anchor_chain(x, root3, head3); });
  }
  check(positions > 50 && refused == positions, "every sampled bit flip in a v6 anchor chain with births refuses");
}
// ---------------------------------------------------------------------------
// 9. The v6 common predicate (coordinator-link notes).
void predicate_controls() {
  namespace cc = raikou::native_water::common_commit;
  namespace c6 = raikou::native_water::common_commit_v6;
  auto P = arrival_policy();
  auto a = fresh();
  StepState st;
  const auto before = *a->full_time_boundary();
  ArrivalRun run(P);
  std::shared_ptr<const m::FullTimeStepFlows> flows;
  std::shared_ptr<const m::FullTimeAcceptedBoundary> after;
  std::string token, plan;
  std::uint64_t base = 0, count = 0;
  auto hooks = run.hooks();
  hooks.before_publish = [&](m::PreparedArena &t) {
    flows = t.candidate_full_time_flows();
    after = t.candidate_full_time_boundary();
    token = run.staged->token();
    plan = run.staged->plan().sha256();
    base = run.staged->native_base();
    count = run.staged->births().size();
  };
  run.scope_begin = before.time_s;
  run.scope_dt = 5e-15;
  StepPlan sp;
  sp.maximum_carriers = 256;
  (void)full_step(*a, sp, st, hooks);
  auto observed = [](const m::FullTimeAcceptedBoundary &b) {
    cc::ObservedOwners o;
    o.native = b;
    o.actual_native_prefix = b.accepted_arena_sha256;
    o.actual_native_revision = b.accepted_arena_revision;
    o.actual_native_step = b.accepted_step;
    o.field = b.field_sha256;
    o.heat = b.heat_sha256;
    o.source = b.source_sha256;
    o.rng = b.rng_sha256;
    o.mesh = b.mesh_sha256;
    o.material = b.material_sha256;
    o.frame = b.frame_sha256;
    o.configuration = sha("coordinator-configuration");
    o.field_revision = o.heat_revision = b.accepted_step;
    o.field_time_s = o.heat_time_s = b.time_s;
    return o;
  };
  cc::ExpectedStep e{before, 5e-15, before.time_s + 5e-15, after->common_clock_sha256,
                     sha("coordinator-configuration")};
  cc::ProposedOwners p{observed(*after), sha("midpoint"), sha("midpoint")};
  c6::ArrivalEvidence ev{true, plan, token, base, count};
  check(!cc::accepts_closed_step(e, observed(before), p),
        "the unchanged V3 predicate refuses a source-changing (birth) step");
  check(c6::accepts_closed_or_birth_step(e, observed(before), p, flows.get(), ev),
        "the v6 predicate accepts the source change bound by the native birth row and the arrival token");
  auto wrong = ev;
  wrong.token_sha256 = sha("other-token");
  check(!c6::accepts_closed_or_birth_step(e, observed(before), p, flows.get(), wrong),
        "v6 predicate refuses a foreign arrival token");
  check(!c6::accepts_closed_or_birth_step(e, observed(before), p, nullptr, ev),
        "v6 predicate refuses a source change without a native birth receipt");
  auto q = p;
  q.after.source = q.after.native.source_sha256 = sha("unbound-source");
  check(!c6::accepts_closed_or_birth_step(e, observed(before), q, flows.get(), ev),
        "v6 predicate refuses a source seal other than the birth row's source_after");
  auto unchanged = p;
  unchanged.after.source = unchanged.after.native.source_sha256 = before.source_sha256;
  check(!c6::accepts_closed_or_birth_step(e, observed(before), unchanged, flows.get(), ev),
        "v6 predicate refuses a birth row whose source transition the candidate does not show");
}
int main(int argc, char **argv) {
  try {
    if (argc != 2 && !(argc == 3 && (std::string(argv[2]) == "reduced" ||
                                     std::string(argv[2]).rfind("long=", 0) == 0)))
      throw std::runtime_error("usage: test_v6_flows OUTPUT_DIRECTORY [reduced|long=N]");
    const std::string dir = argv[1];
    v71::install_from_environment();
    struct ExecutionReport {
      std::string path;
      ~ExecutionReport() { v71::write_report(path); }
    } execution_report{dir + "/execution.json"};
    reduced = argc == 3 && std::string(argv[2]) == "reduced";
    if (argc == 3 && !reduced)
      long_steps = std::strtoull(argv[2] + 5, nullptr, 10);
    const auto t0 = std::chrono::steady_clock::now();
    if (long_steps) {
      long_run(dir);
      std::ofstream(dir + "/v6-flows-long-run.json")
          << J({{"status", "passed"}, {"check_count", checks}, {"checks", results},
                {"measurements", measurements}}).dump(2) << '\n';
      std::cout << checks << " long-run checks passed\n";
      return 0;
    }
    kinetic_port_controls();
    plain_v3_controls();
    birth_controls(dir);
    sink_controls();
    same_step_controls();
    checkpoint_controls(dir);
    anchor_controls();
    predicate_controls();
    long_run(dir);
    check(forbidden_calls == 0, "no native solver/advance call in any v6 control");
    J output = {{"status", "passed"},
                {"check_count", checks},
                {"checks", results},
                {"measurements", measurements},
                {"forbidden_native_solver_or_advance_calls", forbidden_calls},
                {"reduced", reduced},
                {"seconds", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()},
                {"motion_providers",
                 "genuine Hall corrective v2 provider at zero total field through the v7 adapter; "
                 "zero-field reference flows provider (births/sinks) with the sealed Hall first-hit helper"}};
    std::ofstream(dir + "/v6-flows-tests.json") << output.dump(2) << '\n';
    std::cout << checks << " native v6 flows controls passed\n";
    return 0;
  } catch (const std::exception &x) {
    std::cerr << "FAILED after " << checks << " checks: " << x.what() << '\n';
    J output = {{"status", "failed"}, {"check_count", checks}, {"checks", results},
                {"failure", x.what()}, {"measurements", measurements}};
    if (argc >= 2)
      std::ofstream(std::string(argv[1]) + "/v6-flows-tests.json") << output.dump(2) << '\n';
    return 1;
  }
}
