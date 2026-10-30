#pragma once
#include "native_parent_subflow.hpp"
#include "axial_joint_v2.hpp"
#include "common_commit_contract.hpp"

// Integration fixture for resident particles and the conditional parent channel.
// Surface emission/desorption/collection and broader chemistry are not supplied
// by this port. Their absence must never be relabeled as the requested full run.
namespace raikou::native_water {
namespace h = ::hall_driver;
struct ClosedParentInputs {
  std::string species_sha256;
  std::string unchanged_source_sha256;
  std::size_t maximum_native_records=4096;
};
class ClosedParentTrial {
  friend ClosedParentTrial prepare_closed_parent_joint(
    m::PicArena&,ParentSubflow&,const axial::Mesh&,const axial::Fields&,
    const axial::State&,h::Thermal&,const h::CurrentTrace&,
    const ClosedParentInputs&,double,std::size_t,int);
  friend bool can_publish_closed_parent_joint(const m::PicArena&,
    const ParentSubflow&,const axial::State&,const h::Thermal&,
    const std::string&,const ClosedParentTrial&) noexcept;
  friend bool try_publish_closed_parent_joint(m::PicArena&,ParentSubflow&,
    axial::State&,h::Thermal&,const std::string&,ClosedParentTrial&) noexcept;
  std::unique_ptr<m::PrivateFullTimeStep> private_step;
  std::unique_ptr<m::PreparedArena> native;
  PreparedParentState random;
  std::unique_ptr<h::AxialJointTrial> field_heat;
  std::vector<m::ArenaRateEntry> post_first_census;
  std::string source_sha256;
  Json audit_;
  common_commit::ExpectedStep expected;
  bool consumed=false;
  ClosedParentTrial()=default;
public:
  ClosedParentTrial(ClosedParentTrial&&) noexcept=default;
  ClosedParentTrial& operator=(ClosedParentTrial&&) noexcept=default;
  ClosedParentTrial(const ClosedParentTrial&)=delete;
  ClosedParentTrial& operator=(const ClosedParentTrial&)=delete;
  // Diagnostic data only. No candidate owner, token, mutable state or vector
  // escapes this move-only outer holder.
  const Json& audit()const noexcept{return audit_;}
};

// Failure injection is a test hook: every location is still before common
// publication. Only the accepted arena owns particles/counts/IDs permanently.
inline ClosedParentTrial prepare_closed_parent_joint(
    m::PicArena& arena,ParentSubflow& chemistry,const axial::Mesh& mesh,
    const axial::Fields& fields,const axial::State& field,h::Thermal& heat,
    const h::CurrentTrace& current,const ClosedParentInputs& inputs,
    double dt_s,std::size_t maximum_carriers,int failure_point=0){
  pc::require(pc::digest(inputs.species_sha256)&&pc::digest(inputs.unchanged_source_sha256),
    "closed parent integration needs explicit species/source seals");
  pc::require(failure_point>=0&&failure_point<=5,"unknown manufactured failure point");
  const auto fail=[&](int point){if(failure_point==point)throw std::runtime_error("manufactured common-step preparation failure");};
  const auto boundary=arena.full_time_boundary();
  pc::require(boundary&&boundary->field_sha256==h::field_state_sha256(field)&&
    boundary->heat_sha256==h::thermal_state_sha256(heat)&&
    boundary->source_sha256==inputs.unchanged_source_sha256&&
    boundary->mesh_sha256==h::motion_mesh_sha256(mesh,arena.native_options())&&
    boundary->material_sha256==h::material_binding_sha256()&&
    boundary->time_s==field.time_s&&boundary->time_s==heat.time_s()&&
    boundary->accepted_step==field.revision&&boundary->accepted_step==heat.revision(),
    "accepted native/field/thermal/source boundaries differ");
  pc::require(std::isfinite(dt_s)&&dt_s>0&&dt_s<=1e-14,
    "closed parent integration requires the unchanged native full-time dt guard");
  const Json random_before=chemistry.checkpoint();
  const auto first_sequence=random_before.at("next_interval").get<std::uint64_t>();
  pc::require(first_sequence<=UINT64_MAX-3,"two restorable collision exposure counters required");
  m::FullTimeStepBinding common;
  common.accepted_step=boundary->accepted_step;
  common.accepted_arena_revision=boundary->accepted_arena_revision;
  common.begin_s=boundary->time_s;common.dt_s=dt_s;common.end_s=common.begin_s+dt_s;
  common.accepted_arena_sha256=boundary->accepted_arena_sha256;
  common.accepted_field_sha256=boundary->field_sha256;
  common.mesh_sha256=boundary->mesh_sha256;common.material_sha256=boundary->material_sha256;
  common.common_clock_sha256=sha(Json{{"schema","closed-parent-common-clock-v1"},
    {"step",common.accepted_step},{"revision",common.accepted_arena_revision},
    {"begin",common.begin_s},{"end",common.end_s},{"dt",common.dt_s},
    {"arena",common.accepted_arena_sha256},{"field",common.accepted_field_sha256},
    {"heat",boundary->heat_sha256},{"source",boundary->source_sha256},
    {"rng",boundary->rng_sha256},{"mesh",common.mesh_sha256},
    {"material",common.material_sha256},{"frame",boundary->frame_sha256},
    {"species",inputs.species_sha256}}.dump());
  ClosedParentTrial out;out.source_sha256=inputs.unchanged_source_sha256;
  out.expected.before=*boundary;out.expected.dt_s=common.dt_s;
  out.expected.end_s=common.end_s;out.expected.common_clock=common.common_clock_sha256;
  // The corrective Hall Owner must supply its actual immutable configuration
  // identity. This legacy-signature staging copy is not an admitted adapter.
  out.private_step=std::make_unique<m::PrivateFullTimeStep>(
    arena.begin_full_time_step(common,maximum_carriers,inputs.maximum_native_records));
  m::ExposureBinding exposure{common,m::FullTimeSlot::before_transport,
    common.begin_s,.5*common.dt_s,first_sequence,boundary->field_sha256};
  auto first_scope=out.private_step->first_exposure(exposure);
  auto first=chemistry.prepare_on(first_scope);fail(1);
  auto motion_scope=out.private_step->fork_transport_trial();
  out.post_first_census=motion_scope.native_census(maximum_carriers);
  h::MotionBinding motion;
  motion.phase=h::ParticlePhase::full_time_position_and_proper_velocity;
  motion.accepted_step=common.accepted_step;
  motion.staged_arena_revision=motion_scope.staged_arena_revision();
  motion.begin_s=common.begin_s;motion.end_s=common.end_s;motion.dt_s=common.dt_s;
  motion.staged_arena_sha256=motion_scope.state_prefix();
  motion.mesh_sha256=common.mesh_sha256;motion.material_sha256=common.material_sha256;
  motion.species_sha256=inputs.species_sha256;motion.common_clock_sha256=common.common_clock_sha256;
  // Filled by the sealed kinetic-boundary API when that version is installed.
  motion.kinetic_policy_sha256=h::kinetic_policy_sha256(h::declared_kinetic_policy(arena.native_options()));
  out.field_heat=std::make_unique<h::AxialJointTrial>(h::prepare_axial_joint(
    mesh,fields,field,arena.native_options(),out.post_first_census,motion,heat,current));
  fail(2);
  const auto& proposal=*out.field_heat;
  const auto field_after=h::field_state_sha256(proposal.field);
  const auto heat_after=h::thermal_trial_state_sha256(proposal.heat);
  m::FullTimeMotionBinding native_motion{common,proposal.motion->seal_sha256(),
    boundary->field_sha256,field_after,heat_after,inputs.unchanged_source_sha256};
  Json envelope{{"schema","arena_total_field_motion_receipt_v1"},
    {"transport_trial_sha256",native_motion.transport_trial_sha256},
    {"common_clock_sha256",common.common_clock_sha256},{"mesh_sha256",common.mesh_sha256},
    {"field_before_sha256",native_motion.field_before_sha256},{"field_after_sha256",field_after},
    {"begin_s",common.begin_s},{"end_s",common.end_s},{"dt_s",common.dt_s},
    {"owner_audit",h::motion_record(*proposal.motion)}};
  const auto span=proposal.motion->endpoints();
  std::vector<m::GapParticle> endpoints(span.begin(),span.end());
  auto moved=motion_scope.prepare_total_field_motion(native_motion,endpoints,envelope);
  out.private_step->accept_transport_trial(std::move(moved));fail(3);
  exposure.slot=m::FullTimeSlot::after_transport;exposure.frozen_physical_stage_s=common.end_s;
  exposure.owner_stage_sequence=first_sequence+1;exposure.frozen_field_sha256=field_after;
  auto second_scope=out.private_step->second_exposure(exposure);
  auto second=chemistry.prepare_on(second_scope,&first);fail(4);
  out.native=std::make_unique<m::PreparedArena>(out.private_step->seal(
    {field_after,heat_after,inputs.unchanged_source_sha256,chemistry.staged_rng_sha256(second)}));
  out.random=chemistry.seal(first,second,*out.native);fail(5);
  out.audit_={{"schema","closed-parent-field-heat-step-v2-staged"},
    {"begin_s",common.begin_s},{"end_s",common.end_s},{"authoritative_dt_s",common.dt_s},
    {"common_clock_sha256",common.common_clock_sha256},{"first_exposure",first.audit()},
    {"second_exposure",second.audit()},{"transport",h::motion_record(*proposal.motion)},
    {"field_heat_balance_J",proposal.combined_balance_J},{"native",out.native->audit()},
    {"surface_source_and_collection",false},{"conditional_parent_only",true},{"full_plasma_run",false}};
  return out;
}
inline bool can_publish_closed_parent_joint(const m::PicArena& arena,const ParentSubflow& chemistry,
    const axial::State& field,const h::Thermal& heat,const std::string& accepted_source,
    const ClosedParentTrial& trial) noexcept {
  try {
    return !trial.consumed&&trial.native&&trial.field_heat&&accepted_source==trial.source_sha256&&
      arena.can_publish(*trial.native)&&chemistry.can_publish(trial.random,*trial.native)&&
      h::preflight_axial_joint(field,heat,trial.post_first_census,*trial.field_heat);
  }catch(...){return false;}
}
inline bool try_publish_closed_parent_joint(m::PicArena& arena,ParentSubflow& chemistry,
    axial::State& field,h::Thermal& heat,const std::string& accepted_source,
    ClosedParentTrial& trial) noexcept {
  if(!can_publish_closed_parent_joint(arena,chemistry,field,heat,accepted_source,trial))return false;
  if(!arena.try_publish(*trial.native))return false;
  // No allocation, throwing computation, callbacks or I/O between owner swaps.
  chemistry.publish_prevalidated(trial.random);
  h::publish_field_heat_prevalidated(field,heat,std::move(*trial.field_heat));
  trial.consumed=true;
  return true;
}
} // namespace raikou::native_water
