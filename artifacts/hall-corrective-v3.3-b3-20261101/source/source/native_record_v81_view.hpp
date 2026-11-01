#pragma once
#include "motion_receipt.hpp"
#include "pf_v7/forms.hpp"
#include <memory>
#include <cstddef>
#include <type_traits>
namespace hall_corrective_v3 {
namespace nm=raikou::mfem_backend;
static_assert(std::is_standard_layout_v<hall_driver::CarrierPath> && std::is_standard_layout_v<nm::FullTimeMotionPath>);
static_assert(sizeof(hall_driver::CarrierPath)==sizeof(nm::FullTimeMotionPath) && alignof(hall_driver::CarrierPath)==alignof(nm::FullTimeMotionPath));
static_assert(offsetof(hall_driver::CarrierPath,native_id)==offsetof(nm::FullTimeMotionPath,native_id) && offsetof(hall_driver::CarrierPath,field_start_m)==offsetof(nm::FullTimeMotionPath,field_start_m) && offsetof(hall_driver::CarrierPath,field_displacement_m)==offsetof(nm::FullTimeMotionPath,field_displacement_m) && offsetof(hall_driver::CarrierPath,flight_s)==offsetof(nm::FullTimeMotionPath,flight_s));
struct NativeMotionRecordV81 {
 std::vector<nm::FullTimeSinkRecord>sinks;std::vector<nm::FullTimeBirthFlight>births;
 nm::FullTimeMotionRecordV7 record;
 // Optional retention only for an admission that can outlive its proposal.
 std::vector<double> retained_E,retained_B,retained_current,retained_extra;
 std::vector<nm::GapParticle> retained_endpoints;
 std::vector<nm::FullTimeMotionPath> retained_paths;
 std::vector<nm::FullTimeBirthBridge> retained_bridges;
 NativeMotionRecordV81(const NativeMotionRecordV81&)=delete;
 NativeMotionRecordV81&operator=(const NativeMotionRecordV81&)=delete;
 explicit NativeMotionRecordV81(const hall_driver::MotionReceipt&r,std::span<const double>E={},std::span<const double>B={},bool retain=false){
  if(!r.forms_v7())throw std::runtime_error("native v8.1 record requires Hall native-v7 forms");
  sinks.reserve(r.sinks().size());
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
    sinks.push_back(x);
  }
  births.reserve(r.births().size());
  for (const auto& b : r.births()) births.push_back({b.native_id, b.birth_time_s, b.flight_s});
  auto& x = record;
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
  x.paths = paths;
  x.sinks = sinks;
  x.births = births;
  x.trial_E = E;
  x.trial_B = B;
  x.current_A_m = r.current_A_m();
  x.schema=r.bridges().empty()?nm::motion_v7_schema:nm::motion_v81_schema;x.bridges=r.bridges();x.extra_current_A_m=r.extra_current_A_m();
  if(!x.bridges.empty())x.extra_current_A_m_sha256=pf_v7::array_sha256(x.extra_current_A_m);
  x.extra_work_J=r.extra_work_J();x.record_sha256=r.bridges().empty()?nm::full_time_motion_record_v7_sha256(x):nm::full_time_motion_record_v81_sha256(x);
  if(retain){
   retained_E.assign(E.begin(),E.end());retained_B.assign(B.begin(),B.end());retained_current.assign(x.current_A_m.begin(),x.current_A_m.end());retained_extra.assign(x.extra_current_A_m.begin(),x.extra_current_A_m.end());retained_endpoints.assign(x.endpoints.begin(),x.endpoints.end());retained_paths.assign(x.paths.begin(),x.paths.end());retained_bridges.assign(x.bridges.begin(),x.bridges.end());
   x.trial_E=retained_E;x.trial_B=retained_B;x.current_A_m=retained_current;x.extra_current_A_m=retained_extra;x.endpoints=retained_endpoints;x.paths=retained_paths;x.bridges=retained_bridges;
  }
 }
};
template<class Proposal>nm::FullTimeMotionAdmission matching_motion_admission(const Proposal&p){
 auto stored=std::make_shared<NativeMotionRecordV81>(p.motion(),p.trial_E(),p.trial_B(),true);
 if(stored->record.record_sha256!=p.motion().record_v7_sha256())throw std::runtime_error("Hall/native v8.1 motion digest mismatch");
 return nm::bind_total_field_motion_v7(stored->record,stored);
}
}
