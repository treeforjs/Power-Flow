#pragma once
// Copy: adapter of the owner pair Hall corrective v3.4 (resident owner, device coupled step, ResidentOptions of the
// device lane: no host record arrays, block-form field/current digests) + native candidate007.1 (native v7.1 device
// lane, HALL-LINK-7.1.txt): the native motion record v7 of an accepted Hall proposal.
//   header: binding fields (trial_field_sha256 in the block form total-midpoint-E-B-blocks-v1), budget, census v2
//     digest (Hall's, from native's device census), the sealed array digests of the record (current in the block form
//     motion-array-blocks-v1; endpoint and initial charge motion-array-binary64-v1), kinetic policy, frame,
//     neutral-stream bound, record digest = Hall's record_v7_sha256;
//   device: the addresses and counts of the record arrays (endpoints 80 B, paths 64 B, sinks 136 B, births 24 B) and of
//     the trial midpoint E, B and the current J on the shared session (DeviceStep::record_view of the sealed iterate,
//     unchanged until the next prepare of the field owner);
//   host spans: sinks and births (O(sinks + newborns), converted); no endpoints, paths, E, B or J on the host (native's
//     host path, used for refusal texts, verified mode and steps with births or sinks, reads the device arrays back).
// Native recomputes the exact sums, the work tolerance, the section roots, the record digest, the census digest and the
// block-form field/current digests on the device and refuses any difference.
#include "pf_v7/forms.hpp"
#include "raikou/mfem_pic_full_time.h"
#include "resident_owner.hpp"
#include <cstddef>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

namespace raikou::native_water::v7_1::pair_v34_v71 {
namespace hc = ::hall_corrective_v3;
namespace nm = ::raikou::mfem_backend;

static_assert(sizeof(::hall_driver::CarrierPath) == sizeof(nm::FullTimeMotionPath) && sizeof(nm::FullTimeMotionPath) == 64,
              "Hall CarrierPath and native FullTimeMotionPath must have the canonical 64-byte layout");

struct Storage {
  std::vector<nm::FullTimeSinkRecord> sinks;
  std::vector<nm::FullTimeBirthFlight> births;
  nm::FullTimeMotionRecordV7 record;
};

inline nm::FullTimeMotionAdmission admission(const hc::ResidentJointProposal& p, hc::gpu::DeviceStep& step) {
  const auto& r = p.motion();
  if (!r.forms_v7() || r.record_v7_sha256().empty() || !r.device_arrays() || !p.device_path())
    throw std::runtime_error("coordinator v7.1 copy: native v7.1 lane needs a Hall v3.4 device-path proposal with the v7 forms");
  const auto& v = step.record_view();
  if (v.seal_generation != p.device_seal_generation())
    throw std::runtime_error("coordinator v7.1 copy: the field owner's device record arrays belong to a later proposal");
  const auto& digests = r.record_v7_array_digests();
  if (digests[0].empty() || r.binding().trial_field_sha256 != v.trial_field_blocks_sha256 || digests[0] != v.current_blocks_sha256)
    throw std::runtime_error("coordinator v7.1 copy: Hall record without the block-form field/current digests of its device arrays");
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
  x.current_A_m_sha256 = digests[0];
  x.endpoint_charge_C_sha256 = digests[1];
  x.initial_charge_C_sha256 = digests[2];
  const auto& k = r.kinetic_policy();
  x.kinetic_policy = {k.width_m,      k.gap_m,          k.short_z_m,      k.feed_z_m,      k.lower.left_m,  k.lower.right_m,
                      k.lower.bottom_m, k.lower.top_m, k.upper.left_m, k.upper.right_m, k.upper.bottom_m, k.upper.top_m};
  x.native_to_field_m = r.frame().native_to_field_m;
  x.neutral_stream_dt_max_s = r.neutral_stream_dt_max_s();
  x.sinks = s->sinks;
  x.births = s->births;
  x.record_sha256 = r.record_v7_sha256();
  auto& d = x.device;
  d.endpoints = v.endpoints; d.paths = v.paths; d.sinks = v.sinks; d.births = v.births;
  d.n_endpoints = v.n_endpoints; d.n_paths = v.n_paths; d.n_sinks = v.n_sinks; d.n_births = v.n_births;
  d.E = v.E; d.B = v.B; d.J = v.J; d.n_E = v.n_E; d.n_B = v.n_B;
  return nm::bind_total_field_motion_v7(s->record, s);
}
}  // namespace raikou::native_water::v7_1::pair_v34_v71
