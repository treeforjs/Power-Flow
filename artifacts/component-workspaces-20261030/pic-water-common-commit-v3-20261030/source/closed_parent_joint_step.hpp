#pragma once
#include "native_parent_subflow.hpp"
#include "root_api.hpp"
#include "common_commit_contract.hpp"

// Closed resident-particle component only. Surface sources, collecting sinks,
// evolving layers and general water chemistry require a different contract.
namespace raikou::native_water {
namespace h = ::hall_driver;
namespace hc = ::hall_corrective;

struct ClosedParentInputs {
  std::string species_sha256, frame_sha256;
  std::size_t maximum_native_records=4096, maximum_carriers=100000;
};

inline std::string closed_parent_configuration_sha256(
    const std::string& chemistry_configuration,const std::string& field_configuration,
    const ClosedParentInputs& in){
  pc::require(pc::digest(chemistry_configuration)&&pc::digest(field_configuration)&&
    pc::digest(in.species_sha256)&&pc::digest(in.frame_sha256)&&
    in.maximum_native_records>0&&in.maximum_carriers>0,"closed component requires explicit identities and finite resource caps");
  return sha(Json{{"schema","closed-parent-corrective-common-configuration-v3"},
    {"field_configuration",field_configuration},
    {"chemistry_configuration",chemistry_configuration},
    {"species",in.species_sha256},{"frame",in.frame_sha256},
    {"maximum_native_records",in.maximum_native_records},{"maximum_carriers",in.maximum_carriers},
    {"surface_source",false},{"collection",false},{"active_layers",false}}.dump());
}

// The string overload also permits an explicit fresh genesis to bind its
// configuration before constructing the native arena. Runtime preflight always
// obtains these values from the actual owners through this overload.
inline std::string closed_parent_configuration_sha256(
    const ParentSubflow& chemistry,const hc::Owner& owner,const ClosedParentInputs& in){
  return closed_parent_configuration_sha256(chemistry.configuration_sha256(),
    owner.configuration_sha256(),in);
}

// This is the identity of a declared inactive source contract, not an active
// source inventory. Binding the entire closed configuration here keeps it in
// the accepted native genesis/restart boundary across successive steps.
inline std::string closed_parent_inactive_source_sha256(const std::string& config){
  pc::require(pc::digest(config),"inactive source needs a closed configuration identity");
  return sha(Json{{"schema","closed-parent-inactive-source-v3"},
    {"configuration",config},{"births",false},{"sinks",false},{"layer_changes",false}}.dump());
}

inline common_commit::ObservedOwners observe_closed_accepted(
    const m::PicArena& arena,const ParentSubflow& chemistry,const hc::Owner& owner,
    const ClosedParentInputs& in){
  const auto boundary=arena.full_time_boundary();
  pc::require(bool(boundary),"closed component requires an actual full-time native boundary");
  common_commit::ObservedOwners out;out.native=*boundary;
  out.actual_native_prefix=arena.state_prefix();
  const auto snapshot=arena.snapshot();
  out.actual_native_revision=snapshot->revision;out.actual_native_step=snapshot->accepted_steps;
  out.field=owner.field_sha256();out.heat=owner.heat_sha256();
  out.rng=chemistry.accepted_rng_sha256();
  out.mesh=h::motion_mesh_sha256(owner.mesh(),arena.native_options());
  out.material=h::material_binding_sha256();out.frame=in.frame_sha256;
  out.configuration=closed_parent_configuration_sha256(chemistry,owner,in);
  out.source=closed_parent_inactive_source_sha256(out.configuration);
  // Owner admits and privately owns both field and thermal clocks. Its final
  // preflight checks both against this single authoritative aggregate clock.
  out.field_time_s=out.heat_time_s=owner.time_s();
  out.field_revision=out.heat_revision=owner.revision();
  return out;
}

class ClosedParentTrial {
  friend ClosedParentTrial prepare_closed_parent_joint(m::PicArena&,ParentSubflow&,
    hc::Owner&,const ClosedParentInputs&,double,int);
  friend bool can_publish_closed_parent_joint(const m::PicArena&,const ParentSubflow&,
    const hc::Owner&,const ClosedParentTrial&) noexcept;
  friend bool try_publish_closed_parent_joint(m::PicArena&,ParentSubflow&,
    hc::Owner&,ClosedParentTrial&) noexcept;
  std::unique_ptr<m::PrivateFullTimeStep> private_step_;
  std::unique_ptr<m::PreparedArena> native_;
  PreparedParentState random_;
  std::unique_ptr<hc::JointProposal> field_heat_;
  std::vector<m::ArenaRateEntry> post_first_census_;
  ClosedParentInputs inputs_;
  common_commit::ExpectedStep expected_;
  Json audit_;
  bool consumed_=false;
  ClosedParentTrial()=default;
public:
  ClosedParentTrial(ClosedParentTrial&&) noexcept=default;
  ClosedParentTrial& operator=(ClosedParentTrial&&) noexcept=default;
  ClosedParentTrial(const ClosedParentTrial&)=delete;
  ClosedParentTrial& operator=(const ClosedParentTrial&)=delete;
  const Json& audit()const noexcept{return audit_;}
};

inline ClosedParentTrial prepare_closed_parent_joint(
    m::PicArena& arena,ParentSubflow& chemistry,hc::Owner& owner,
    const ClosedParentInputs& inputs,double dt_s,int failure_point=0){
  pc::require(failure_point>=0&&failure_point<=5,"unknown manufactured failure point");
  const auto fail=[&](int point){if(failure_point==point)throw std::runtime_error("manufactured common-step preparation failure");};
  const auto before=observe_closed_accepted(arena,chemistry,owner,inputs);
  pc::require(common_commit::owner_values_match(before),"accepted native/field/heat/RNG/closed-source boundaries differ");
  pc::require(std::isfinite(dt_s)&&dt_s>0&&dt_s<=1e-14&&.5*dt_s>0,
    "closed component requires the unchanged native full-time dt guard");
  const auto& boundary=before.native;
  const auto first_sequence=chemistry.checkpoint().at("next_interval").get<std::uint64_t>();
  pc::require(first_sequence<=UINT64_MAX-3,"two restorable collision exposure counters required");
  m::FullTimeStepBinding common;
  common.accepted_step=boundary.accepted_step;common.accepted_arena_revision=boundary.accepted_arena_revision;
  common.begin_s=boundary.time_s;common.dt_s=dt_s;common.end_s=common.begin_s+dt_s;
  pc::require(std::isfinite(common.end_s)&&common.end_s>common.begin_s,"common endpoint is not representable");
  common.accepted_arena_sha256=before.actual_native_prefix;common.accepted_field_sha256=before.field;
  common.mesh_sha256=before.mesh;common.material_sha256=before.material;
  common.common_clock_sha256=sha(Json{{"schema","closed-parent-common-clock-v3"},
    {"step",common.accepted_step},{"revision",common.accepted_arena_revision},
    {"begin",common.begin_s},{"end",common.end_s},{"dt",common.dt_s},
    {"arena",before.actual_native_prefix},{"field",before.field},{"heat",before.heat},
    {"source",before.source},{"rng",before.rng},{"mesh",before.mesh},
    {"material",before.material},{"frame",before.frame},{"configuration",before.configuration}}.dump());
  ClosedParentTrial out;out.inputs_=inputs;
  out.expected_={boundary,common.dt_s,common.end_s,common.common_clock_sha256,before.configuration};
  out.private_step_=std::make_unique<m::PrivateFullTimeStep>(
    arena.begin_full_time_step(common,inputs.maximum_carriers,inputs.maximum_native_records));
  m::ExposureBinding exposure{common,m::FullTimeSlot::before_transport,
    common.begin_s,.5*common.dt_s,first_sequence,before.field};
  auto first_scope=out.private_step_->first_exposure(exposure);
  auto first=chemistry.prepare_on(first_scope);fail(1);
  auto motion_scope=out.private_step_->fork_transport_trial();
  out.post_first_census_=motion_scope.native_census(inputs.maximum_carriers);
  h::MotionBinding motion;
  motion.phase=h::ParticlePhase::full_time_position_and_proper_velocity;
  motion.accepted_step=common.accepted_step;motion.staged_arena_revision=motion_scope.staged_arena_revision();
  motion.begin_s=common.begin_s;motion.end_s=common.end_s;motion.dt_s=common.dt_s;
  motion.staged_arena_sha256=motion_scope.state_prefix();motion.mesh_sha256=before.mesh;
  motion.material_sha256=before.material;motion.species_sha256=inputs.species_sha256;
  motion.common_clock_sha256=common.common_clock_sha256;
  motion.kinetic_policy_sha256=h::kinetic_policy_sha256(h::declared_kinetic_policy(arena.native_options()));
  out.field_heat_=std::make_unique<hc::JointProposal>(hc::prepare_closed_proposal(
    owner,arena.native_options(),out.post_first_census_,motion,before.source));fail(2);
  const auto& proposal=*out.field_heat_;
  const auto& field_after=proposal.candidate_field_sha256();
  const auto& heat_after=proposal.candidate_heat_sha256();
  m::FullTimeMotionBinding native_motion{common,proposal.motion().seal_sha256(),
    before.field,field_after,heat_after,before.source};
  // Live admission comes only from the genuine immutable Hall object. Audit
  // JSON and copied endpoint vectors cannot authorize private native motion.
  auto moved=motion_scope.prepare_total_field_motion(native_motion,
    m::bind_hall_total_field_motion(proposal.motion()));
  out.private_step_->accept_transport_trial(std::move(moved));fail(3);
  exposure.slot=m::FullTimeSlot::after_transport;exposure.frozen_physical_stage_s=common.end_s;
  exposure.owner_stage_sequence=first_sequence+1;exposure.frozen_field_sha256=field_after;
  auto second_scope=out.private_step_->second_exposure(exposure);
  auto second=chemistry.prepare_on(second_scope,&first);fail(4);
  out.native_=std::make_unique<m::PreparedArena>(out.private_step_->seal(
    {field_after,heat_after,before.source,chemistry.staged_rng_sha256(second)}));
  out.random_=chemistry.seal(first,second,*out.native_);fail(5);
  out.audit_={{"schema","closed-parent-corrective-common-step-v3"},
    {"begin_s",common.begin_s},{"end_s",common.end_s},{"authoritative_dt_s",common.dt_s},
    {"common_clock_sha256",common.common_clock_sha256},{"configuration_sha256",before.configuration},
    {"first_exposure",first.audit()},{"second_exposure",second.audit()},
    {"transport",h::motion_record(proposal.motion())},{"midpoint_role_sha256",proposal.midpoint().seal()},
    {"field_heat_balance_J",proposal.combined_balance_J()},{"native",out.native_->audit()},
    {"surface_source_and_collection",false},{"conditional_parent_only",true},{"full_plasma_run",false}};
  return out;
}

inline bool can_publish_closed_parent_joint(const m::PicArena& arena,const ParentSubflow& chemistry,
    const hc::Owner& owner,const ClosedParentTrial& trial) noexcept{
  try{
    if(trial.consumed_||!trial.native_||!trial.field_heat_||
       !arena.can_publish(*trial.native_)||!chemistry.can_publish(trial.random_,*trial.native_))return false;
    const auto boundary=trial.native_->candidate_full_time_boundary();
    const auto exposures=trial.native_->candidate_full_time_exposures();
    const auto& proposal=*trial.field_heat_;
    if(!boundary||!exposures||exposures->transport_trial_sha256!=proposal.motion().seal_sha256())return false;
    const auto current=observe_closed_accepted(arena,chemistry,owner,trial.inputs_);
    common_commit::ProposedOwners candidate;auto& after=candidate.after;
    after.native=*boundary;after.actual_native_prefix=trial.native_->candidate_state_prefix();
    // These counters are the immutable native candidate's own values, validated
    // by native can_publish; there is no second counter owner or JSON parser.
    after.actual_native_revision=boundary->accepted_arena_revision;
    after.actual_native_step=boundary->accepted_step;
    after.field=proposal.candidate_field_sha256();after.heat=proposal.candidate_heat_sha256();
    after.source=current.source;after.rng=chemistry.candidate_rng_sha256(trial.random_,*trial.native_);
    after.mesh=current.mesh;after.material=current.material;after.frame=current.frame;
    after.configuration=current.configuration;
    after.field_time_s=after.heat_time_s=proposal.end_s();
    after.field_revision=after.heat_revision=proposal.endpoint_revision();
    candidate.motion_midpoint=proposal.motion().binding().trial_field_sha256;
    candidate.field_midpoint=h::total_field_sha256(proposal.midpoint().E(),proposal.midpoint().B());
    if(proposal.configuration_sha256()!=owner.configuration_sha256()||
       !common_commit::accepts_closed_step(trial.expected_,current,candidate))return false;
    // Final Hall approval is last: all hashing, allocation and fallible checks
    // complete before the first accepted native/RNG/field/heat publication.
    return owner.preflight_final(proposal,trial.post_first_census_,*boundary);
  }catch(...){return false;}
}

inline bool try_publish_closed_parent_joint(m::PicArena& arena,ParentSubflow& chemistry,
    hc::Owner& owner,ClosedParentTrial& trial) noexcept{
  if(!can_publish_closed_parent_joint(arena,chemistry,owner,trial))return false;
  if(!arena.try_publish(*trial.native_))return false;
  chemistry.publish_prevalidated(trial.random_);
  owner.publish_prevalidated(std::move(*trial.field_heat_));
  trial.consumed_=true;
  return true;
}
} // namespace raikou::native_water
