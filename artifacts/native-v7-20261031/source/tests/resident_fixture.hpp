#pragma once
// Resident full-time step fixture shared by the long-run driver and the v5
// history tests. It reuses the candidate003 storage-control species, channel,
// genesis seals and exposure/completion recipes (tests/full_time_fixture.hpp)
// on the Hall-compatible geometry (z in [-depth/2, +depth/2]); motion comes
// from the genuine Hall v2 provider (prepare_canonical_motion) evaluated with
// identically zero total-field arrays. No field solve is performed; field and
// heat seals after each step are synthetic digests.
// Native v7: the receipt enters through the transitional v7 adapter
// (source/adapters/hall_receipt_v7.hpp: exact sums, Merkle record, census v2).
// Compiled with PARITY_NATIVE_V6 against the sealed v6 headers it uses the v6
// JSON route (parity trace).
#include "full_time_fixture.hpp"
#include "hall-motion-v2/motion_receipt.hpp"
#ifndef PARITY_NATIVE_V6
#include "hall_receipt_v7.hpp"
#endif
#include "hall-finite-foil-em-20261029/axial_mesh.hpp"
#include <chrono>
#include <cstdio>
#include <sys/resource.h>

namespace resident {
inline m::GapPicOptions hall_options() {
  auto o = options();
  o.short_z_m = -.0125;
  o.feed_z_m = .0125;
  return o;
}
inline const axial::Mesh &mesh() {
  static const axial::Mesh g;
  return g;
}
inline m::FullTimeGenesisBinding hall_genesis() {
  auto g = genesis();
  g.mesh_sha256 = hall_driver::motion_mesh_sha256(mesh(), hall_options());
  return g;
}
inline std::vector<m::FullTimeAuthoredPacket> population() {
  return {{"electron",
           "e",
           "manufactured/source",
           {4, 0},
           {.004, .004, 0},
           {0, 0, 1e6}},
          {"neutral",
           "n",
           "manufactured/source",
           {8, 0},
           {.00403, .004, 0},
           {0, 0, 0}}};
}
// n single-carrier packets (alternating electrons moving along +z and neutrals
// at rest) on a lattice strictly inside the inter-foil gap.
inline std::vector<m::FullTimeAuthoredPacket> lattice(std::size_t n) {
  std::vector<m::FullTimeAuthoredPacket> out;
  out.reserve(n);
  const m::ArenaCount four(4, 0);
  for (std::size_t i = 0; i < n; ++i) {
    const double x = .0021 + 4.6e-3 * double(i % 997) / 997.0;
    const double y = .0021 + 4.6e-3 * double((i / 997) % 1009) / 1009.0;
    const double z = -.01 + 1e-8 * double(i / (997 * 1009));
    char id[24];
    std::snprintf(id, sizeof id, "p%07zu", i);
    const bool electron = i % 2 == 0;
    out.push_back({id, electron ? "e" : "n", "manufactured/lattice", four,
                   {x, y, z}, {0, 0, electron ? 1e6 : 0}});
  }
  return out;
}
inline std::string step_text(const char *what, std::uint64_t k) {
  return std::string(what) + "/" + std::to_string(k);
}
// One complete closed full-time step: zero-event (or forced elastic) exposures
// around one genuine Hall v2 motion. Returns false only if publication refuses.
struct StepOptions {
  double dt_s = 5e-15;
  std::size_t maximum_carriers = 8;
  std::size_t maximum_records = 32;
  int first_events = 0;  // forced manufactured elastic events, first slot
  int second_events = 0; // second slot
};
struct StepState {
  std::uint64_t sequence = 0; // chemistry stage sequence (two per step)
  std::string rng = genesis().rng_sha256;
  std::uint64_t events = 0;
  std::string target = "neutral";
  // Phase timers (s): native owner calls vs the Hall provider and its record.
  double begin_s = 0, first_s = 0, census_s = 0, hall_s = 0, admit_s = 0,
         second_s = 0, seal_s = 0, publish_s = 0;
};
inline double since(std::chrono::steady_clock::time_point &t) {
  const auto now = std::chrono::steady_clock::now();
  const double d = std::chrono::duration<double>(now - t).count();
  t = now;
  return d;
}
inline m::ArenaEvent forced_elastic(m::CollisionScope &s, const std::string &id,
                                    const std::string &target) {
  m::ArenaEvent e;
  e.id = id;
  e.channel = channels()[0].id;
  e.primary = "electron";
  e.target = target;
  e.event_weight = 4;
  e.scattered_primary_id = e.primary;
  e.phase = s.phase();
  const auto c = channels()[0];
  e.provenance = {m::ArenaProvenance::Kind::ConditionalPhysicalData,
                  c.dataset_sha256,
                  c.kernel_manifest_sha256,
                  c.kernel_binding_file_sha256,
                  c.model,
                  "forced_kernel_integration_only",
                  s.state_prefix(),
                  0};
  // Zero total field: the electron keeps its genesis proper velocity, so the
  // elastic products close energy and momentum exactly.
  e.products = {{e.primary, "e", {0, 0, 1e6}, m::ArenaProduct::At::Primary},
                {id + "/neutral", "n", {0, 0, 0}, m::ArenaProduct::At::Target}};
  return e;
}
inline void exposure_events(m::CollisionScope &scope, int events,
                            StepState &state, std::uint64_t step) {
  for (int i = 0; i < events; ++i) {
    const auto id = "manufactured/forced/" + std::to_string(step) + "/" +
                    std::to_string(state.events++);
    auto child = scope.prepare(forced_elastic(scope, id, state.target));
    if (!scope.try_publish(child))
      throw std::runtime_error("forced private event publication");
    state.target = id + "/neutral";
  }
}
inline void full_step(m::PicArena &a, const StepOptions &opt, StepState &st) {
  namespace h = hall_driver;
  const auto before = a.full_time_boundary();
  const std::uint64_t k = before->accepted_step;
  m::FullTimeStepBinding b{before->accepted_step,
                           before->accepted_arena_revision,
                           before->time_s,
                           before->time_s + opt.dt_s,
                           opt.dt_s,
                           before->accepted_arena_sha256,
                           before->field_sha256,
                           before->mesh_sha256,
                           before->material_sha256,
                           sha(step_text("common-clock", k))};
  auto clock0 = std::chrono::steady_clock::now();
  auto outer =
      a.begin_full_time_step(b, opt.maximum_carriers, opt.maximum_records);
  st.begin_s += since(clock0);
  m::ExposureBinding e1{b,         m::FullTimeSlot::before_transport,
                        b.begin_s, .5 * b.dt_s,
                        st.sequence, b.accepted_field_sha256};
  auto completion_for = [&](const m::ExposureBinding &e,
                            const std::string &before_rng,
                            const std::string &after_rng) {
    auto c = completion(e);
    c.before_rng_clock_sha256 = before_rng;
    c.after_rng_clock_sha256 = after_rng;
    return c;
  };
  const auto rng_a = sha(step_text("chemistry-rng-first", k));
  const auto rng_b = sha(step_text("chemistry-rng-second", k));
  auto s1 = outer.first_exposure(e1);
  // One declared reaction budget M for both slots (second must be <= first).
  const auto budget = std::size_t(std::max({1, opt.first_events, opt.second_events}));
  auto t1 = s1.prepare_compound_with_completion(
      budget, [&](m::CollisionScope &scope) {
        exposure_events(scope, opt.first_events, st, k);
        return completion_for(e1, st.rng, rng_a);
      });
  if (!s1.try_publish(t1))
    throw std::runtime_error("first exposure publication");
  st.first_s += since(clock0);
  auto trial = outer.fork_transport_trial();
  const auto census = trial.native_census(opt.maximum_carriers);
  const auto staged_revision = trial.staged_arena_revision();
  const auto staged_prefix = trial.state_prefix();
  st.census_s += since(clock0);
  static const std::vector<double> E(mesh().ne, 0.0), B(mesh().nb, 0.0);
  static const std::string field_mid = h::total_field_sha256(E, B);
  h::MotionBinding mb;
  mb.phase = h::ParticlePhase::full_time_position_and_proper_velocity;
  mb.accepted_step = b.accepted_step;
  mb.staged_arena_revision = staged_revision;
  mb.begin_s = b.begin_s;
  mb.end_s = b.end_s;
  mb.dt_s = b.dt_s;
  mb.staged_arena_sha256 = staged_prefix;
  mb.trial_field_sha256 = field_mid;
  mb.mesh_sha256 = b.mesh_sha256;
  mb.material_sha256 = b.material_sha256;
  mb.species_sha256 = sha("resident-fixture-species");
  mb.common_clock_sha256 = b.common_clock_sha256;
  mb.kinetic_policy_sha256 =
      h::kinetic_policy_sha256(h::declared_kinetic_policy(a.native_options()));
  auto receipt = h::prepare_canonical_motion(mesh(), a.native_options(), census,
                                             mb, {}, E, B);
  const auto field_after = sha(step_text("field-after", k));
  const auto heat_after = sha(step_text("heat-after", k));
#ifdef PARITY_NATIVE_V6
  m::FullTimeMotionBinding binding{b,           receipt.seal_sha256(),
                                   b.accepted_field_sha256, field_after,
                                   heat_after,  before->source_sha256};
  auto admission = m::bind_hall_total_field_motion(receipt);
#else
  auto v7 = m::v7_adapter::from_hall_receipt(
      receipt, trial, E, B, a.full_time_motion_policy().neutral_stream_dt_max_s);
  m::FullTimeMotionBinding binding{b,           v7->record.record_sha256,
                                   b.accepted_field_sha256, field_after,
                                   heat_after,  before->source_sha256};
  auto admission = m::v7_adapter::admit(v7);
#endif
  st.hall_s += since(clock0);
  auto moved = trial.prepare_total_field_motion(binding, std::move(admission));
  outer.accept_transport_trial(std::move(moved));
  st.admit_s += since(clock0);
  m::ExposureBinding e2 = e1;
  e2.slot = m::FullTimeSlot::after_transport;
  e2.frozen_physical_stage_s = b.end_s;
  e2.owner_stage_sequence = st.sequence + 1;
  e2.frozen_field_sha256 = field_after;
  auto s2 = outer.second_exposure(e2);
  auto t2 = s2.prepare_compound_with_completion(
      budget, [&](m::CollisionScope &scope) {
        exposure_events(scope, opt.second_events, st, k);
        return completion_for(e2, rng_a, rng_b);
      });
  if (!s2.try_publish(t2))
    throw std::runtime_error("second exposure publication");
  st.second_s += since(clock0);
  auto token = outer.seal({field_after, heat_after, before->source_sha256, rng_b});
  st.seal_s += since(clock0);
  if (!a.can_publish(token) || !a.try_publish(token))
    throw std::runtime_error("outer publication refused");
  st.publish_s += since(clock0);
  st.sequence += 2;
  st.rng = rng_b;
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
} // namespace resident
