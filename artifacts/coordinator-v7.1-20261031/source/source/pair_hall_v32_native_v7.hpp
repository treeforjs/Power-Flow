#pragma once
// Coordinator v7.1 adapter of the owner pair Hall corrective v3.2 (native v7 motion forms) +
// native candidate007 (native v7): the native motion record v7 of an accepted Hall proposal
// (API-007.txt 1.2 with amendment A1; HALL-LINK.txt 3e).
//   header: binding fields, budget (Hall v3.2 forms: exact kinetic change and field work),
//     census v2 digest (Hall's, computed with the forms library from the census of the native
//     trial and native's state table), array digests (motion-array-binary64-v1), kinetic
//     policy, frame, neutral-stream bound, record digest = Hall's record_v7_sha256;
//   spans: endpoints (the receipt's GapParticle array), paths (the receipt's CarrierPath array,
//     the same 64-byte layout as FullTimeMotionPath, checked below), sinks and births (converted,
//     O(sinks + newborns));
//   A1 admission inputs: the proposal's trial midpoint E and B and the plasma current.
// Native recomputes the exact sums, the work tolerance, the four section roots, the record digest
// and the census digest and refuses any difference. The spans reference the proposal, which
// lives in the trial until publication; the converted arrays live in the admission storage.
#include "pf_v7/forms.hpp"
#include "raikou/mfem_pic_full_time.h"
#include "root_api.hpp"
#include <cstddef>
#include <memory>
#include <span>
#include <vector>

namespace raikou::native_water::v7_1::pair_v32_v7 {
namespace hc = ::hall_corrective_v3;
namespace nm = ::raikou::mfem_backend;

static_assert(sizeof(::hall_driver::CarrierPath) == sizeof(nm::FullTimeMotionPath) && sizeof(nm::FullTimeMotionPath) == 64,
              "Hall CarrierPath and native FullTimeMotionPath must have the canonical 64-byte layout");
static_assert(offsetof(::hall_driver::CarrierPath, native_id) == offsetof(nm::FullTimeMotionPath, native_id) &&
                  offsetof(::hall_driver::CarrierPath, field_start_m) == offsetof(nm::FullTimeMotionPath, field_start_m) &&
                  offsetof(::hall_driver::CarrierPath, field_displacement_m) == offsetof(nm::FullTimeMotionPath, field_displacement_m) &&
                  offsetof(::hall_driver::CarrierPath, flight_s) == offsetof(nm::FullTimeMotionPath, flight_s),
              "Hall CarrierPath and native FullTimeMotionPath field offsets differ");

struct Storage {
  std::vector<nm::FullTimeSinkRecord> sinks;
  std::vector<nm::FullTimeBirthFlight> births;
  nm::FullTimeMotionRecordV7 record;
};

inline nm::FullTimeMotionAdmission admission(const hc::JointProposal& p) {
  const auto& r = p.motion();
  if (!r.forms_v7() || r.record_v7_sha256().empty())
    throw std::runtime_error("coordinator v7.1: Hall proposal without the native v7 motion forms (Model MotionForms::native_v7)");
  auto s = std::make_shared<Storage>();
  s->sinks.reserve(r.sinks().size());
  for (const auto& k : r.sinks()) {
    nm::FullTimeSinkRecord x;
    x.native_id = k.native_id;
    x.kind = nm::FullTimeSinkKind(k.kind);
    x.surface = nm::FullTimeSinkSurface(k.boundary);
    x.axis = std::uint32_t(k.axis);
    x.plane_m = k.plane_m;
    x.fraction = k.fraction;
    x.segment_m = k.segment_m;
    x.hit_point_m = k.hit_point_m;
    x.incident_proper_velocity_m_s = k.incident_proper_velocity_m_s;
    x.flight_s = k.flight_s;
    x.discarded_flight_s = k.discarded_flight_s;
    s->sinks.push_back(x);
  }
  s->births.reserve(r.births().size());
  for (const auto& b : r.births()) s->births.push_back({b.native_id, b.birth_time_s, b.flight_s});
  auto& x = s->record;
  const auto& b = r.binding();
  x.accepted_step = b.accepted_step;
  x.staged_arena_revision = b.staged_arena_revision;
  x.trial_stage = b.trial_stage;
  x.begin_s = b.begin_s;
  x.end_s = b.end_s;
  x.dt_s = b.dt_s;
  x.staged_arena_sha256 = b.staged_arena_sha256;
  x.trial_field_sha256 = b.trial_field_sha256;
  x.mesh_sha256 = b.mesh_sha256;
  x.material_sha256 = b.material_sha256;
  x.species_sha256 = b.species_sha256;
  x.common_clock_sha256 = b.common_clock_sha256;
  x.kinetic_policy_sha256 = b.kinetic_policy_sha256;
  const auto& u = r.budget();
  x.budget = {u.kinetic_change_J, u.field_work_J, u.work_defect_J, u.continuity_C, u.plasma_dt, u.maximum_gyro, u.maximum_displacement_cells};
  x.census_sha256 = r.census_sha256();
  x.current_A_m_sha256 = pf_v7::array_sha256(r.current_A_m());
  x.endpoint_charge_C_sha256 = pf_v7::array_sha256(r.endpoint_charge_C());
  x.initial_charge_C_sha256 = pf_v7::array_sha256(r.initial_charge_C());
  x.endpoints = r.endpoints();
  const auto& k = r.kinetic_policy();
  x.kinetic_policy = {k.width_m,      k.gap_m,          k.short_z_m,      k.feed_z_m,      k.lower.left_m,  k.lower.right_m,
                      k.lower.bottom_m, k.lower.top_m, k.upper.left_m, k.upper.right_m, k.upper.bottom_m, k.upper.top_m};
  x.native_to_field_m = r.frame().native_to_field_m;
  x.neutral_stream_dt_max_s = r.neutral_stream_dt_max_s();
  const auto paths = r.paths();
  x.paths = {reinterpret_cast<const nm::FullTimeMotionPath*>(paths.data()), paths.size()};
  x.sinks = s->sinks;
  x.births = s->births;
  x.trial_E = p.trial_E();
  x.trial_B = p.trial_B();
  x.current_A_m = p.current_A_m();
  x.record_sha256 = r.record_v7_sha256();
  return nm::bind_total_field_motion_v7(s->record, s);
}
}  // namespace raikou::native_water::v7_1::pair_v32_v7
