// Diagnostic fixture derived from pinned native chem_v8_fixture.hpp.
// SHA256 4bcd7663d241cf12f4524952dabff8915e9210bce58ba33882f9de8895e6127c
// Only changes: live stage inspection hooks and abandon instead of outer seal/publication.
#pragma once
#include "chem_v8_fixture.hpp"
namespace adapter_fixture {
using namespace chemfix;
inline StepReport staged_step(Arena &A, double dt, int k,
  const std::function<void(m::CollisionScope &,m::PreparedArena &)> &inspect) {
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
  inspect(s1,t1);
  if (!s1.try_publish(t1))
    throw std::runtime_error("first chemistry exposure publication");
  T.stop(R.phases["stage1"]);
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
  inspect(s2,t2);
  if (!s2.try_publish(t2))
    throw std::runtime_error("second chemistry exposure publication");
  T.stop(R.phases["stage2"]);
  // No outer publication: all staged native population and RNG changes roll back.
  return R;
}
}
