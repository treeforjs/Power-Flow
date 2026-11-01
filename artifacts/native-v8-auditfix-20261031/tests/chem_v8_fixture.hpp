#pragma once
// Native v8 chemistry fixture: a full-time arena with identity schema v8 and
// chemistry v3 (the production blob, read in place), a load population
// (electrons, H2O+ and H2O carriers), and complete outer steps: begin, stage 1
// (prepare_chemistry_v3), fork, census, a straight-drift motion record at zero
// field built from the trial census, accept, stage 2 with the release, seal and
// publish. Every native call is timed (thread CPU and wall).
#include "raikou/mfem_pic_full_time.h"
#include "raikou/mfem_pic_canonical.h"
#include "pf_v7/forms.hpp"
#include "hall-canonical-driver-20261029/kinetic_boundary.hpp"
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace chemfix {
namespace m = raikou::mfem_backend;
// Tables: chemistry v3 (default) or v3.1 (PF_V8_CHEM=v3.1); native v8 compiles the v3.1 headers, which give
// the v3 records for the v3 tables.
inline bool v31() {
  const char *e = std::getenv("PF_V8_CHEM");
  return e && std::string(e) == "v3.1";
}
inline const char *blob_path() {
  return v31() ? "/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031/data/tables/chemistry-v3.1.blob"
               : "/home/treefor/projects/runtimes/pic-water-chemistry-v3-20261031/data/tables/chemistry-v3.blob";
}
inline const char *blob_sha() {
  return v31() ? "70fcc2adc2fdaa62abed7ed02d68ffc09b1b18554eda8fa99c5451120ae454b8"
               : "3ca3ea495630e6b6e125c05ebf45488fd89a899544a2aa138cefb9b1dd415b0a";
}
inline std::string sha(const std::string &s) { return m::canonical::sha256_hex(s); }
inline double cpu_now() {
  timespec t;
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
  return double(t.tv_sec) + 1e-9 * double(t.tv_nsec);
}
inline double wall_now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
struct Phase {
  double cpu = 0, wall = 0;
};
struct Timer {
  double c = 0, w = 0;
  void start() {
    c = cpu_now();
    w = wall_now();
  }
  void stop(Phase &p) {
    p.cpu += cpu_now() - c;
    p.wall += wall_now() - w;
  }
};
// Hall axial mesh transverse cells: 45 (200 um) or 180 (50 um).
inline m::GapPicOptions options(int n = 45) {
  m::GapPicOptions o;
  o.nx = o.ny = n;
  o.width_m = o.gap_m = .009;
  o.depth_m = .025;
  o.finite_foils = o.finite_axial_extent = true;
  o.short_z_m = -.0125;
  o.feed_z_m = .0125;
  o.lower_foil = {.002, .007, .0018, .002};
  o.upper_foil = {.002, .007, .007, .0072};
  return o;
}
struct Population {
  std::size_t carriers = 1000;
  double electron_fraction = .49, ion_fraction = .49; // rest: H2O
  double water_weight = 1e6, charged_weight = 1;
  double electron_eV_min = 1, electron_eV_max = 200;
  std::uint64_t seed = 7;
  bool clustered = false; // place all carriers in a few collision cells
};
inline std::vector<m::FullTimeAuthoredPacket> population(const Population &P) {
  std::vector<m::FullTimeAuthoredPacket> out;
  out.reserve(P.carriers);
  std::mt19937_64 rng(P.seed);
  std::uniform_real_distribution<double> U(0, 1);
  const double me_c2_eV = 510998.95069, c = 299792458.0;
  for (std::size_t i = 0; i < P.carriers; ++i) {
    double x, y, z;
    if (P.clustered) {
      x = .0021 + 1e-4 * U(rng);
      y = .0021 + 1e-4 * U(rng);
      z = -.01 + 1e-3 * U(rng);
    } else {
      x = .0021 + 4.6e-3 * U(rng);
      y = .0021 + 4.6e-3 * U(rng);
      z = -.0115 + .023 * U(rng);
    }
    char id[32];
    std::snprintf(id, sizeof id, "load/%07zu", i);
    const double f = double(i % 100) / 100.0;
    std::string state;
    double w = P.charged_weight;
    m::PicVector u{0, 0, 0};
    if (f < P.electron_fraction) {
      state = "pfc-v3/e";
      const double E = P.electron_eV_min * std::pow(P.electron_eV_max / P.electron_eV_min, U(rng));
      const double gamma = 1 + E / me_c2_eV;
      const double un = c * std::sqrt(gamma * gamma - 1);
      const double ct = 2 * U(rng) - 1, ph = 6.283185307179586 * U(rng), st = std::sqrt(1 - ct * ct);
      u = {un * st * std::cos(ph), un * st * std::sin(ph), un * ct};
    } else if (f < P.electron_fraction + P.ion_fraction) {
      state = "pfc-v3/H2O+";
      for (auto &v : u)
        v = 400 * (2 * U(rng) - 1);
    } else {
      state = "pfc-v3/H2O";
      w = P.water_weight;
      for (auto &v : u)
        v = 400 * (2 * U(rng) - 1);
    }
    out.push_back({id, state, "manufactured/chemistry-load", m::ArenaCount(w, 0), {x, y, z}, u});
  }
  return out;
}
inline m::FullTimeChemistryConfiguration chemistry_configuration(const m::GapPicOptions &o, unsigned partitions = 1) {
  m::FullTimeChemistryConfiguration c;
  c.enabled = true;
  c.diagnostic_field_target_approximation = true; // component qualification, not full-system closure
  c.blob_path = blob_path();
  c.blob_sha256 = blob_sha();
  c.run_seed = 0x5EED5EED20261031ull;
  c.release_weight = 1e6;
  c.grid_counts = {std::uint32_t(o.nx), std::uint32_t(o.ny), 5};
  c.host_partitions = partitions;
  c.retain_stage_inputs = false;
  // Admitted: recommended list (placeholders, 1e7 eV decade) and the regions
  // below the tables (target-frame energies under the first grid point:
  // carriers co-moving with their targets).
  const auto reg = m::full_time_chemistry_registry([&] {
    auto x = c;
    return x;
  }(), o);
  // Fixture policy: every unbounded region admitted (exposure recorded), bounded
  // regions reported only; the recommended list is exercised by the refusal
  // controls.
  for (std::size_t i = 0; i < reg.region_names.size(); ++i)
    if (!reg.region_bounded[i])
      c.policy.admitted_regions.push_back(reg.region_names[i]);
  c.policy.event_limit = 1e30; // bounded regions: report only (fixture)
  return c;
}
struct Arena {
  m::GapPicOptions o;
  m::FullTimeConfiguration cfg;
  m::FullTimeChemistryRegistry reg;
  std::vector<m::ArenaStateSpec> states;
  std::vector<m::ArenaChannel> channels;
  m::FullTimeGenesisBinding g;
  std::unique_ptr<m::PicArena> arena;
};
inline Arena make_arena(const m::GapPicOptions &o, const m::FullTimeChemistryConfiguration &chem,
                        const std::vector<m::FullTimeAuthoredPacket> &packets, std::size_t cap,
                        std::uint64_t anchor_interval = 1024) {
  Arena a;
  a.cfg.history.anchor_interval_steps = anchor_interval;
  a.o = o;
  a.cfg.identity = m::FullTimeIdentitySchema::v8;
  a.cfg.chemistry = chem;
  a.cfg.motion.neutral_stream_dt_max_s = 1e-11; // steps up to 10 ps (chemistry controls)
  a.cfg.limits.maximum_carriers = std::max<std::size_t>(cap, 1000);
  a.cfg.limits.memory_budget_bytes =
      a.cfg.limits.maximum_carriers * m::full_time_bytes_per_carrier + m::full_time_fixed_bytes;
  a.cfg.limits.maximum_identity_entries = std::max<std::size_t>(2 * cap, 1u << 22);
  a.reg = m::full_time_chemistry_registry(chem, o);
  a.states = m::full_time_chemistry_states(a.reg);
  a.channels = {a.reg.channel};
  a.g.time_s = 26.31715e-9;
  a.g.origin_provenance = "manufactured/chemistry-v8-fixture";
  for (auto *x : {&a.g.field_sha256, &a.g.mesh_sha256, &a.g.material_sha256, &a.g.frame_sha256,
                  &a.g.common_clock_sha256, &a.g.heat_sha256, &a.g.source_sha256})
    *x = sha("chem-fixture/" + std::to_string(std::uintptr_t(x) - std::uintptr_t(&a.g)));
  a.g.rng_sha256 = a.reg.initial_rng_sha256;
  a.arena = m::PicArena::create_full_time(o, a.states, a.channels, packets, a.g, a.cfg);
  return a;
}
struct StepReport {
  std::shared_ptr<const m::FullTimeChemistryStageReceipt> stage1, stage2;
  std::map<std::string, Phase> phases;
  std::string boundary_sha256, census_sha256;
  std::size_t census = 0;
};
// One complete outer step of length dt with straight drift at zero field.
// before_motion / before_seal: optional hooks (tests).
inline StepReport step(Arena &A, double dt, int k,
                       const std::function<void(m::PrivateFullTimeStep &)> &before_motion = {},
                       const std::function<void(m::PrivateFullTimeStep &)> &before_seal = {}) {
  StepReport R;
  Timer T;
  auto &arena = *A.arena;
  const auto &o = A.o;
  const auto bnd = arena.full_time_boundary();
  m::FullTimeStepBinding b{bnd->accepted_step, bnd->accepted_arena_revision, bnd->time_s, bnd->time_s + dt, dt,
                           bnd->accepted_arena_sha256, bnd->field_sha256, bnd->mesh_sha256, bnd->material_sha256,
                           sha("clock/" + std::to_string(k))};
  const auto stages = arena.full_time_chemistry().stages;
  T.start();
  auto outer = arena.begin_full_time_step(b, A.cfg.limits.maximum_carriers, 64);
  T.stop(R.phases["begin"]);
  m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s, .5 * dt, stages, b.accepted_field_sha256};
  T.start();
  auto s1 = outer.first_exposure(e1);
  auto t1 = s1.prepare_chemistry_v3();
  R.stage1 = t1.candidate_chemistry_stage();
  if (!s1.try_publish(t1))
    throw std::runtime_error("first chemistry exposure publication");
  T.stop(R.phases["stage1"]);
  if (before_motion)
    before_motion(outer);
  T.start();
  auto trial = outer.fork_transport_trial();
  R.census_sha256 = trial.census_sha256();
  T.stop(R.phases["fork_census"]);
  // ---- producer (untimed): straight drift of the trial census
  const auto &census = trial.census_v7();
  std::vector<m::GapParticle> rows(census.carriers.begin(), census.carriers.end());
  R.census = rows.size();
  const m::PicVector shift{-.5 * o.width_m, -.5 * o.gap_m, 0};
  std::vector<m::GapParticle> ends = rows;
  std::vector<m::FullTimeMotionPath> paths(rows.size());
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const auto &p = rows[i];
    double u2 = 0;
    for (double v : p.proper_velocity_m_s)
      u2 += v * v;
    const double gamma = std::sqrt(1 + u2 / (m::pic_c * m::pic_c));
    for (int d = 0; d < 3; ++d)
      ends[i].position_m[d] = p.position_m[d] + dt * p.proper_velocity_m_s[d] / gamma;
    auto &pa = paths[i];
    pa.native_id = p.id;
    pa.flight_s = dt;
    for (int d = 0; d < 3; ++d) {
      pa.field_start_m[d] = p.position_m[d] + shift[d];
      pa.field_displacement_m[d] = (ends[i].position_m[d] + shift[d]) - pa.field_start_m[d];
    }
  }
  const pf_v7::HallAxialMesh hm{std::uint64_t(o.nx), 6};
  static std::map<std::uint64_t, std::vector<double>> Es, Bs, Js;
  auto &E = Es[hm.n], &B = Bs[hm.n], &Jc = Js[hm.n];
  if (E.empty()) {
    E.assign(hm.ne(), 1e3);
    B.assign(hm.nb(), 0.0);
    Jc.assign(hm.ne(), 0.0);
  }
  m::FullTimeMotionRecordV7 x;
  x.accepted_step = b.accepted_step;
  x.staged_arena_revision = trial.staged_arena_revision();
  x.begin_s = b.begin_s;
  x.end_s = b.end_s;
  x.dt_s = dt;
  x.staged_arena_sha256 = trial.state_prefix();
  x.trial_field_sha256 = pf_v7::trial_field_blocks_sha256(E, B);
  x.mesh_sha256 = b.mesh_sha256;
  x.material_sha256 = b.material_sha256;
  x.species_sha256 = sha("species");
  x.common_clock_sha256 = b.common_clock_sha256;
  x.kinetic_policy_sha256 = hall_driver::kinetic_policy_sha256(hall_driver::declared_kinetic_policy(o));
  std::vector<m::GapParticle> start_rows(ends.size());
  for (std::size_t i = 0; i < ends.size(); ++i) {
    start_rows[i] = ends[i];
    for (int d = 0; d < 3; ++d)
      start_rows[i].position_m[d] = paths[i].field_start_m[d] - shift[d];
  }
  const auto K = pf_v7::kinetic_change_J({reinterpret_cast<const pf_v7::Carrier *>(start_rows.data()), start_rows.size()},
                                         {reinterpret_cast<const pf_v7::Carrier *>(ends.data()), ends.size()}, {});
  const auto W = pf_v7::field_work_J(dt, E, Jc);
  x.budget = {K.value, W.value, K.value - W.value, 0, 0, 0, 0};
  x.census_sha256 = R.census_sha256;
  x.current_A_m_sha256 = pf_v7::current_blocks_sha256(Jc);
  x.endpoint_charge_C_sha256 = sha("endpoint-charge");
  x.initial_charge_C_sha256 = sha("initial-charge");
  x.kinetic_policy = {o.width_m, o.gap_m, o.short_z_m, o.feed_z_m, o.lower_foil.left_m, o.lower_foil.right_m,
                      o.lower_foil.bottom_m, o.lower_foil.top_m, o.upper_foil.left_m, o.upper_foil.right_m,
                      o.upper_foil.bottom_m, o.upper_foil.top_m};
  x.native_to_field_m = shift;
  x.neutral_stream_dt_max_s = A.cfg.motion.neutral_stream_dt_max_s;
  x.endpoints = ends;
  x.paths = paths;
  x.record_sha256 = m::full_time_motion_record_v7_sha256(x);
  x.trial_E = E;
  x.trial_B = B;
  x.current_A_m = Jc;
  m::FullTimeMotionBinding mb{b, x.record_sha256, b.accepted_field_sha256, sha("f/" + std::to_string(k)),
                              sha("h/" + std::to_string(k)), bnd->source_sha256};
  T.start();
  auto moved = trial.prepare_total_field_motion(mb, m::bind_total_field_motion_v7(x));
  outer.accept_transport_trial(std::move(moved));
  T.stop(R.phases["admission_accept"]);
  m::ExposureBinding e2 = e1;
  e2.slot = m::FullTimeSlot::after_transport;
  e2.frozen_physical_stage_s = b.end_s;
  e2.owner_stage_sequence = stages + 1;
  e2.frozen_field_sha256 = mb.field_after_sha256;
  T.start();
  auto s2 = outer.second_exposure(e2);
  auto t2 = s2.prepare_chemistry_v3();
  R.stage2 = t2.candidate_chemistry_stage();
  if (!s2.try_publish(t2))
    throw std::runtime_error("second chemistry exposure publication");
  T.stop(R.phases["stage2"]);
  if (before_seal)
    before_seal(outer);
  T.start();
  auto token = outer.seal({mb.field_after_sha256, mb.heat_after_sha256, bnd->source_sha256, R.stage2->state_after_sha256});
  if (!arena.can_publish(token) || !arena.try_publish(token))
    throw std::runtime_error("publish");
  T.stop(R.phases["seal_publish"]);
  R.boundary_sha256 = arena.full_time_boundary()->accepted_arena_sha256;
  return R;
}
} // namespace chemfix
