#pragma once
// Native v6 test fixture: storage-control species/channel/genesis of the v5
// fixture on the Hall geometry (z in [-depth/2, +depth/2]), a complete closed
// full-time step driver with optional births and a pluggable motion source,
// and a zero-field reference motion provider for binary flows records.
//
// The reference provider is a kinematic test oracle, not a field solver: at
// identically zero total field every carrier (charged or neutral) drifts
// straight, D = F u / gamma (the sealed Hall neutral-stream arithmetic), with
// F = dt for resident carriers and F = end - tau for newborns. Its first hit is
// the sealed Hall first_kinetic_hit (kinetic_boundary.cpp at 8dddec0) for an
// interior start and a transcription of the same arithmetic for a foil-face
// birth (Hall's helper requires an interior start). The record digest is
// computed by an independent byte encoder written from API-006.txt, so native
// admission recomputing it is a format check.
#include "raikou/mfem_pic_full_time.h"
#include "raikou/mfem_pic_canonical.h"
#include "hall-canonical-driver-20261029/kinetic_boundary.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <openssl/sha.h>
#include <sstream>

namespace m = raikou::mfem_backend;
using J = m::ArenaJson;
namespace v6 {
inline int checks = 0;
inline J results = J::array();
inline void check(bool x, const std::string &name) {
  ++checks;
  results.push_back({{"name", name}, {"passed", x}});
  if (!x)
    throw std::runtime_error("check failed: " + name);
}
template <class F> std::string message_of(F f) {
  try {
    f();
  } catch (const std::exception &e) {
    return e.what();
  }
  return {};
}
template <class F> bool refuses(F f) { return !message_of(f).empty(); }
template <class F>
bool flow_refuses(F f, m::FullTimeFlowRefusal::Kind kind) {
  try {
    f();
  } catch (const m::FullTimeFlowRefusal &e) {
    if (e.kind() != kind)
      std::cerr << "  refusal kind " << int(e.kind()) << " (expected "
                << int(kind) << "): " << e.what() << "\n";
    return e.kind() == kind;
  } catch (const std::exception &e) {
    std::cerr << "  untyped refusal: " << e.what() << "\n";
    return false;
  }
  std::cerr << "  no refusal\n";
  return false;
}
template <class F>
bool resource_refuses(F f, m::FullTimeResourceRefusal::Kind kind) {
  try {
    f();
  } catch (const m::FullTimeResourceRefusal &e) {
    return e.kind() == kind;
  } catch (const std::exception &e) {
    std::cerr << "  other refusal: " << e.what() << "\n";
    return false;
  }
  return false;
}
inline std::string sha(const std::string &s) {
  return m::canonical::sha256_hex(s);
}
inline m::GapPicOptions options() {
  m::GapPicOptions o;
  o.nx = o.ny = 45;
  o.width_m = o.gap_m = .009;
  o.depth_m = .025;
  o.finite_foils = true;
  o.finite_axial_extent = true;
  o.short_z_m = -.0125;
  o.feed_z_m = .0125;
  o.lower_foil = {.002, .007, .0018, .002};
  o.upper_foil = {.002, .007, .007, .0072};
  return o;
}
inline std::vector<m::ArenaStateSpec> states() {
  return {{"e", 0, 0, -1, m::pic_me, 0, 0},
          {"n", 2, 1, 0, 2.9915e-26, 0, 0},
          {"i", 2, 1, 1, 2.9915e-26 - m::pic_me, 0, 0}};
}
inline std::vector<m::ArenaChannel> channels() {
  return {{"conditional-storage-control",
           "e",
           "n",
           {"e", "n"},
           sha("data"),
           sha("kernel"),
           sha("binding"),
           "manufactured-storage-control-only"}};
}
inline m::FullTimeGenesisBinding genesis() {
  m::FullTimeGenesisBinding g;
  g.time_s = 26.31715e-9;
  g.origin_provenance = "manufactured/full-time-v6-flows-control-no-simulation";
  g.field_sha256 = sha("field");
  g.mesh_sha256 = sha("mesh");
  g.material_sha256 = sha("material");
  g.frame_sha256 = sha("frame");
  g.common_clock_sha256 = sha("initial-clock");
  g.heat_sha256 = sha("heat");
  g.source_sha256 = sha("source/0");
  g.rng_sha256 = sha("chemistry-RNG-and-config");
  return g;
}
inline std::vector<m::FullTimeAuthoredPacket> base_population() {
  return {{"electron", "e", "manufactured/source", {4, 0}, {.004, .004, 0},
           {0, 0, 1e6}},
          {"neutral", "n", "manufactured/source", {8, 0}, {.00403, .004, 0},
           {0, 0, 0}}};
}
inline m::ArenaStochasticCompletion completion(const m::ExposureBinding &b,
                                               const std::string &before,
                                               const std::string &after) {
  const auto c = channels()[0];
  m::ArenaStochasticCompletion x;
  x.integration_owner_id = "storage-control-owner";
  x.before_rng_clock_sha256 = before;
  x.after_rng_clock_sha256 = after;
  x.common_clock_sha256 = b.common.common_clock_sha256;
  x.channel_id = c.id;
  x.dataset_sha256 = c.dataset_sha256;
  x.kernel_manifest_sha256 = c.kernel_manifest_sha256;
  x.kernel_binding_file_sha256 = c.kernel_binding_file_sha256;
  x.model = c.model;
  x.stage_sequence = b.owner_stage_sequence;
  x.frozen_physical_stage_s = b.frozen_physical_stage_s;
  x.common_interval_begin_s = b.common.begin_s;
  x.common_interval_end_s = b.common.end_s;
  return x;
}
// ---------------------------------------------------------------------------
// Independent byte encoder (API-006.txt section 4) for the reference seal.
struct Bytes {
  std::string b;
  void u64(std::uint64_t v) {
    for (int i = 0; i < 8; ++i)
      b.push_back(char((v >> (8 * i)) & 255));
  }
  void f64(double x) {
    std::uint64_t u;
    std::memcpy(&u, &x, 8);
    u64(u);
  }
  void text(const std::string &s) {
    u64(s.size());
    b += s;
  }
  void vec(const m::PicVector &v) {
    for (double x : v)
      f64(x);
  }
  std::string sha() const { return m::canonical::sha256_hex(b); }
};
inline std::string array_digest(const std::vector<double> &v) {
  Bytes e;
  e.text("motion-array-binary64-v1");
  e.u64(v.size());
  for (double x : v)
    e.f64(x);
  return e.sha();
}
inline std::string reference_record_sha(const m::FullTimeMotionRecordV3 &r) {
  Bytes e;
  const bool flows = r.schema == m::motion_v3_flows_schema;
  e.text(r.schema);
  e.text("canonical-full-time-total-field-motion-v2");
  e.text("binding");
  e.text("full_time_position_and_proper_velocity");
  e.u64(r.accepted_step);
  e.u64(r.staged_arena_revision);
  e.u64(r.trial_stage);
  e.f64(r.begin_s);
  e.f64(r.end_s);
  e.f64(r.dt_s);
  for (const auto *x : {&r.staged_arena_sha256, &r.trial_field_sha256,
                        &r.mesh_sha256, &r.material_sha256, &r.species_sha256,
                        &r.common_clock_sha256, &r.kinetic_policy_sha256})
    e.text(*x);
  e.text("budget");
  for (double x : r.budget)
    e.f64(x);
  e.text("census_sha256");
  e.text(r.census_sha256);
  e.text("coverage");
  e.text(r.coverage);
  e.text("current_A_m_sha256");
  e.text(r.current_A_m_sha256);
  e.text("endpoint_charge_C_sha256");
  e.text(r.endpoint_charge_C_sha256);
  e.text("endpoints");
  e.u64(r.endpoints.size());
  for (const auto &p : r.endpoints) {
    e.u64(p.id);
    e.vec(p.position_m);
    e.vec(p.proper_velocity_m_s);
    e.f64(p.charge_C);
    e.f64(p.mass_kg);
    e.f64(p.weight);
  }
  e.text("initial_charge_C_sha256");
  e.text(r.initial_charge_C_sha256);
  e.text("kinetic_policy");
  e.text("native-finite-target-kinetic-policy-v1");
  for (double x : r.kinetic_policy)
    e.f64(x);
  for (const char *x : {"lower foil x span and inner y gap only",
                        "absorbing escape", "axial escape",
                        "x y lower upper axial; later exact tie wins"})
    e.text(x);
  e.text("native_DKD_called");
  e.u64(r.native_DKD_called);
  e.text("native_Poisson_called");
  e.u64(r.native_Poisson_called);
  e.text("native_to_field_m");
  e.vec(r.native_to_field_m);
  e.text("neutral_motion_owner");
  e.text(r.neutral_motion_owner);
  e.text("paths");
  e.u64(r.paths.size());
  for (const auto &p : r.paths) {
    e.u64(p.native_id);
    e.vec(p.field_start_m);
    e.vec(p.field_displacement_m);
    e.f64(p.flight_s);
  }
  e.text("ordinary_total_field");
  if (flows) {
    e.text("sinks");
    e.u64(r.sinks.size());
    for (const auto &k : r.sinks) {
      e.u64(k.native_id);
      e.u64(std::uint64_t(k.kind));
      e.u64(std::uint64_t(k.surface));
      e.u64(k.axis);
      e.f64(k.plane_m);
      e.f64(k.fraction);
      e.vec(k.segment_m);
      e.vec(k.hit_point_m);
      e.vec(k.incident_proper_velocity_m_s);
      e.f64(k.flight_s);
      e.f64(k.discarded_flight_s);
    }
    e.text("births");
    e.u64(r.births.size());
    for (const auto &x : r.births) {
      e.u64(x.native_id);
      e.f64(x.birth_time_s);
      e.f64(x.flight_s);
    }
  }
  return e.sha();
}
// Transcribed first-hit arithmetic (for a foil-face start; Hall's helper
// refuses non-interior starts). Identical operations to first_kinetic_hit.
inline hall_driver::KineticHit transcribed_first_hit(const m::GapPicOptions &o,
                                                     const m::PicVector &x,
                                                     const m::PicVector &D) {
  hall_driver::KineticHit hit;
  auto take = [&](double f, hall_driver::KineticSink sink, int axis,
                  double plane) {
    if (f < 0 || f > hit.fraction)
      return;
    hit = {sink, f, axis, plane, {}};
  };
  for (int axis = 0; axis < 2; ++axis)
    if (D[axis] != 0) {
      double plane = D[axis] > 0 ? (axis ? o.gap_m : o.width_m) : 0;
      take((plane - x[axis]) / D[axis], hall_driver::KineticSink::outer_xy_escape,
           axis, plane);
    }
  int body = 0;
  for (const auto &foil : {o.lower_foil, o.upper_foil}) {
    double enter = -std::numeric_limits<double>::infinity(),
           leave = std::numeric_limits<double>::infinity(), plane = 0;
    int enter_axis = -1;
    bool possible = true;
    for (int axis = 0; axis < 2; ++axis) {
      double low = axis ? foil.bottom_m : foil.left_m,
             high = axis ? foil.top_m : foil.right_m;
      if (D[axis] == 0) {
        if (x[axis] < low || x[axis] > high)
          possible = false;
        continue;
      }
      double a = (low - x[axis]) / D[axis], b = (high - x[axis]) / D[axis],
             entry_plane = low;
      if (a > b) {
        std::swap(a, b);
        entry_plane = high;
      }
      if (a >= enter) {
        enter = a;
        enter_axis = axis;
        plane = entry_plane;
      }
      leave = std::min(leave, b);
    }
    if (possible && leave > 0 && enter >= 0 && enter <= leave && enter_axis >= 0)
      take(enter,
           body ? hall_driver::KineticSink::upper_foil
                : hall_driver::KineticSink::lower_foil,
           enter_axis, plane);
    ++body;
  }
  if (D[2] != 0) {
    double plane = D[2] > 0 ? o.feed_z_m : o.short_z_m,
           f = (plane - x[2]) / D[2];
    const double xx = x[0] + f * D[0], yy = x[1] + f * D[1];
    bool bridge = D[2] < 0 && xx >= o.lower_foil.left_m &&
                  xx <= o.lower_foil.right_m && yy >= o.lower_foil.top_m &&
                  yy <= o.upper_foil.bottom_m;
    take(f,
         bridge ? hall_driver::KineticSink::short_bridge
                : hall_driver::KineticSink::axial_escape,
         2, plane);
  }
  return hit;
}
// Reference flows motion for one trial (zero field, straight drift).
struct ReferenceMotion {
  std::vector<m::GapParticle> endpoints;
  std::vector<m::FullTimeMotionPath> paths;
  std::vector<m::FullTimeSinkRecord> sinks;
  std::vector<m::FullTimeBirthFlight> births;
  std::vector<m::GapParticle> census; // the trial census (ascending ID)
  m::FullTimeMotionRecordV3 record;
  std::uint64_t hall_hits = 0, transcribed_hits = 0;
  void seal() {
    record.endpoints = endpoints;
    record.paths = paths;
    record.sinks = sinks;
    record.births = births;
    record.record_sha256 = reference_record_sha(record);
  }
};
inline std::shared_ptr<ReferenceMotion>
reference_motion(const m::PrivateTransportTrial &trial,
                 const m::FullTimeStepBinding &b, const m::GapPicOptions &o,
                 std::size_t cap, const std::string &census_sha256,
                 bool flows = true) {
  auto r = std::make_shared<ReferenceMotion>();
  auto census = trial.native_census(cap);
  std::sort(census.begin(), census.end(),
            [](const auto &x, const auto &y) { return x.native.id < y.native.id; });
  r->births = trial.birth_flights();
  const auto policy = hall_driver::declared_kinetic_policy(o);
  const m::PicVector shift{-.5 * o.width_m, -.5 * o.gap_m, 0};
  std::size_t nb = 0;
  long double kinetic = 0;
  double cells = 0;
  for (const auto &row : census) {
    const auto &p = row.native;
    r->census.push_back(p);
    double F = b.dt_s;
    if (nb < r->births.size() && r->births[nb].native_id == p.id)
      F = r->births[nb++].flight_s;
    double u2 = 0;
    for (double v : p.proper_velocity_m_s)
      u2 += v * v;
    const double gamma = std::sqrt(1 + u2 / (m::pic_c * m::pic_c));
    m::PicVector D{};
    for (int d = 0; d < 3; ++d)
      D[d] = F * p.proper_velocity_m_s[d] / gamma;
    hall_driver::KineticHit hit;
    bool interior = true;
    try {
      hall_driver::require_kinetic_interior(policy, p.position_m);
    } catch (const std::exception &) {
      interior = false;
    }
    if (interior) {
      hit = hall_driver::first_kinetic_hit(policy, p.position_m, D);
      ++r->hall_hits;
    } else {
      hit = transcribed_first_hit(o, p.position_m, D);
      ++r->transcribed_hits;
    }
    m::FullTimeMotionPath path;
    path.native_id = p.id;
    for (int d = 0; d < 3; ++d)
      path.field_start_m[d] = p.position_m[d] + shift[d];
    cells = std::max(cells, std::sqrt(D[0] * D[0] + D[1] * D[1] + D[2] * D[2]) /
                                (o.width_m / o.nx));
    if (hit.sink != hall_driver::KineticSink::none) {
      m::FullTimeSinkRecord k;
      k.native_id = p.id;
      k.surface = static_cast<m::FullTimeSinkSurface>(int(hit.sink));
      const int s = int(hit.sink);
      k.kind = (s == 0 || s == 1 || s == 3)
                   ? m::FullTimeSinkKind::collected_at_conductor
                   : m::FullTimeSinkKind::escaped_exterior;
      k.axis = std::uint32_t(hit.axis);
      k.plane_m = hit.plane_m;
      k.fraction = hit.fraction;
      k.segment_m = D;
      for (int d = 0; d < 3; ++d)
        k.hit_point_m[d] = p.position_m[d] + hit.fraction * D[d];
      k.hit_point_m[hit.axis] = hit.plane_m;
      k.incident_proper_velocity_m_s = p.proper_velocity_m_s;
      k.flight_s = F * hit.fraction;
      k.discarded_flight_s = F - k.flight_s;
      for (int d = 0; d < 3; ++d)
        path.field_displacement_m[d] =
            (k.hit_point_m[d] + shift[d]) - path.field_start_m[d];
      path.flight_s = k.flight_s;
      r->sinks.push_back(k);
      m::GapParticle at = p;
      at.proper_velocity_m_s = k.incident_proper_velocity_m_s;
      kinetic += static_cast<long double>(p.weight) *
                 (m::pic_kinetic_energy_J(at) - m::pic_kinetic_energy_J(p));
    } else {
      auto e = p;
      for (int d = 0; d < 3; ++d)
        e.position_m[d] = p.position_m[d] + D[d];
      for (int d = 0; d < 3; ++d)
        path.field_displacement_m[d] =
            (e.position_m[d] + shift[d]) - path.field_start_m[d];
      path.flight_s = F;
      r->endpoints.push_back(e);
      kinetic += static_cast<long double>(p.weight) *
                 (m::pic_kinetic_energy_J(e) - m::pic_kinetic_energy_J(p));
    }
    r->paths.push_back(path);
  }
  auto &x = r->record;
  x.schema = flows ? m::motion_v3_flows_schema : m::motion_v3_schema;
  x.accepted_step = b.accepted_step;
  x.staged_arena_revision = trial.staged_arena_revision();
  x.trial_stage = 0;
  x.begin_s = b.begin_s;
  x.end_s = b.end_s;
  x.dt_s = b.dt_s;
  x.staged_arena_sha256 = trial.state_prefix();
  x.trial_field_sha256 = sha("zero-total-field");
  x.mesh_sha256 = b.mesh_sha256;
  x.material_sha256 = b.material_sha256;
  x.species_sha256 = sha("v6-fixture-species");
  x.common_clock_sha256 = b.common_clock_sha256;
  x.kinetic_policy_sha256 = hall_driver::kinetic_policy_sha256(policy);
  x.budget = {double(kinetic), 0, double(kinetic) - 0, 0, 0, 0, cells};
  x.census_sha256 = census_sha256;
  x.coverage = flows ? "all_native_carriers_births_and_first_hit_sinks"
                     : "all_native_carriers_no_sources_no_boundaries";
  x.current_A_m_sha256 = array_digest(std::vector<double>(8, 0.0));
  x.initial_charge_C_sha256 = array_digest(std::vector<double>(4, 0.0));
  x.endpoint_charge_C_sha256 = array_digest(std::vector<double>(4, 0.0));
  x.kinetic_policy = {o.width_m,          o.gap_m,
                      o.short_z_m,        o.feed_z_m,
                      o.lower_foil.left_m, o.lower_foil.right_m,
                      o.lower_foil.bottom_m, o.lower_foil.top_m,
                      o.upper_foil.left_m, o.upper_foil.right_m,
                      o.upper_foil.bottom_m, o.upper_foil.top_m};
  x.native_to_field_m = shift;
  x.neutral_motion_owner = "same_full_time_q0_drift";
  r->seal();
  return r;
}
inline m::FullTimeMotionAdmission admit(const std::shared_ptr<ReferenceMotion> &r) {
  return m::bind_total_field_motion_v3(r->record, r);
}
// ---------------------------------------------------------------------------
// Step driver.
struct StepPlan {
  double dt_s = 5e-15;
  std::size_t maximum_carriers = 64;
  std::size_t maximum_records = 32;
  int first_events = 0, second_events = 0;
  std::string source_after; // empty: unchanged accepted source seal
};
struct StepState {
  std::uint64_t sequence = 0;
  std::string rng = genesis().rng_sha256;
  std::uint64_t events = 0;
};
inline m::FullTimeStepBinding binding_of(const m::PicArena &a, double dt,
                                         const std::string &clock) {
  const auto b = a.full_time_boundary();
  return {b->accepted_step,  b->accepted_arena_revision, b->time_s,
          b->time_s + dt,    dt,
          b->accepted_arena_sha256, b->field_sha256, b->mesh_sha256,
          b->material_sha256, clock};
}
inline std::string step_text(const char *what, std::uint64_t k) {
  return std::string(what) + "/" + std::to_string(k);
}
inline double since(std::chrono::steady_clock::time_point &t) {
  const auto now = std::chrono::steady_clock::now();
  const double d = std::chrono::duration<double>(now - t).count();
  t = now;
  return d;
}
inline long status_kib(const char *key) {
  FILE *f = std::fopen("/proc/self/status", "r");
  if (!f)
    return -1;
  char line[256];
  long value = -1;
  const auto n = std::strlen(key);
  while (std::fgets(line, sizeof line, f))
    if (!std::strncmp(line, key, n)) {
      value = std::strtol(line + n + 1, nullptr, 10);
      break;
    }
  std::fclose(f);
  return value;
}
// Hooks of one step: births (given the scope) and the motion admission.
struct Hooks {
  using Motion = std::function<m::PreparedPrivateMotion(
      m::PrivateTransportTrial &, m::FullTimeMotionBinding &)>;
  std::function<std::optional<m::FullTimeBirthCredit>(m::BirthScope &)> births;
  Motion motion;
  std::function<void(m::PreparedArena &)> before_publish;
};
struct StepResult {
  std::shared_ptr<const m::FullTimeStepFlows> flows;
  std::uint64_t birth_base = 0;
};
inline StepResult full_step(m::PicArena &a, const StepPlan &opt, StepState &st,
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
      births.emplace(scope.prepare(*credit));
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
// Motion hook: reference flows record (zero field) bound as typed binary v3.
// `mutate` (optional) edits the record before the admission; `reseal`
// recomputes the reference digest afterwards (semantic refusal tests).
inline Hooks::Motion reference_hook(
    std::size_t cap, std::function<void(ReferenceMotion &)> mutate = {},
    bool reseal = true, bool flows = true,
    std::shared_ptr<ReferenceMotion> *keep = nullptr) {
  return [=](m::PrivateTransportTrial &trial, m::FullTimeMotionBinding &binding) {
    auto r = reference_motion(trial, binding.common, options(), cap,
                              trial.census_sha256(), flows);
    if (mutate) {
      mutate(*r);
      r->record.endpoints = r->endpoints;
      r->record.paths = r->paths;
      r->record.sinks = r->sinks;
      r->record.births = r->births;
      if (reseal)
        r->record.record_sha256 = reference_record_sha(r->record);
    }
    if (keep)
      *keep = r;
    binding.transport_trial_sha256 = r->record.record_sha256;
    return trial.prepare_total_field_motion(binding, admit(r));
  };
}
} // namespace v6
