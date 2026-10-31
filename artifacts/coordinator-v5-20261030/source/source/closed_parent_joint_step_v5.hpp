#pragma once
#include "identity_v5.hpp"
#include "common_commit_contract.hpp"
#include <chrono>

// Closed resident-particle common step, version 5 = the version 4 flow linked
// against Hall corrective v2 (binary identities, explicit field-solve backend
// bound into the configuration identity, explicit motion provider, node-plane
// deposit-rounding tolerance) and still native candidate003. Surface sources,
// collecting sinks, evolving layers and general water chemistry need a
// different contract. Test instrumentation is reachable only through
// TrialTestAccess, which is declared here and never defined in source.
namespace raikou::native_water::v5 {
struct TrialTestAccess;

struct ClosedParentInputsV5 {
  std::vector<m::ArenaStateSpec> states;
  std::vector<m::ArenaChannel> channels;
  std::string species_sha256, frame_sha256;
  std::size_t maximum_native_records=4096, maximum_carriers=100000;
  // Diagnostic only (not an identity): include the complete Hall motion record
  // in the trial audit. Multi-step timing runs keep only its seal.
  bool audit_transport_record=true;
};

inline std::string closed_parent_configuration_sha256(
    const std::string& chemistry_configuration,const std::string& field_configuration,
    const ClosedParentInputsV5& in){
  pc::require(pc::digest(chemistry_configuration)&&pc::digest(field_configuration)&&
    pc::digest(in.species_sha256)&&pc::digest(in.frame_sha256)&&
    in.maximum_native_records>0&&in.maximum_carriers>0,"closed component requires explicit identities and finite resource caps");
  pc::require(species_sha256(in.states,in.channels)==in.species_sha256,"declared species identity differs from its catalog");
  Encoding e;e.text("closed-parent-corrective-common-configuration-v4");
  e.text(field_configuration);e.text(chemistry_configuration);e.text(in.species_sha256);e.text(in.frame_sha256);
  e.u64(in.maximum_native_records);e.u64(in.maximum_carriers);
  e.u64(0);e.u64(0);e.u64(0); // surface source, collection, active layers: absent
  return e.finish();
}

inline std::string closed_parent_configuration_sha256(
    const ParentSubflow& chemistry,const hc::Owner& owner,const ClosedParentInputsV5& in){
  return closed_parent_configuration_sha256(chemistry.configuration_sha256(),owner.configuration_sha256(),in);
}

// Declared inactive source contract (no births, sinks or layer changes).
inline std::string closed_parent_inactive_source_sha256(const std::string& config){
  pc::require(pc::digest(config),"inactive source needs a closed configuration identity");
  Encoding e;e.text("closed-parent-inactive-source-v4");e.text(config);
  e.u64(0);e.u64(0);e.u64(0); // births, sinks, layer changes: none
  return e.finish();
}

inline common_commit::ObservedOwners observe_closed_accepted(
    const m::PicArena& arena,const ParentSubflow& chemistry,const hc::Owner& owner,
    const ClosedParentInputsV5& in){
  const auto boundary=arena.full_time_boundary();
  pc::require(bool(boundary),"closed component requires an actual full-time native boundary");
  common_commit::ObservedOwners out;out.native=*boundary;
  out.actual_native_prefix=arena.state_prefix();
  const auto snapshot=arena.snapshot();
  out.actual_native_revision=snapshot->revision;out.actual_native_step=snapshot->accepted_steps;
  out.field=owner.field_sha256();out.heat=owner.heat_sha256();
  out.rng=chemistry.accepted_rng_sha256();
  out.mesh=owner.model().motion_mesh_sha256(arena.native_options());
  pc::require(out.mesh==owner.motion_mesh_sha256(),"Hall owner geometry differs from the live native geometry");
  out.material=owner.model().material_sha256();
  require_arena_species(arena,in.states,in.channels,in.maximum_carriers);
  out.frame=frame_sha256(arena.native_options(),owner.model());
  pc::require(out.frame==in.frame_sha256,"declared frame differs from the live native/Hall translation frame");
  out.configuration=closed_parent_configuration_sha256(chemistry,owner,in);
  out.source=closed_parent_inactive_source_sha256(out.configuration);
  out.field_time_s=out.heat_time_s=owner.time_s();
  out.field_revision=out.heat_revision=owner.revision();
  return out;
}

struct PhaseTimesV5 {
  double observe_s=0,begin_s=0,first_exposure_s=0,census_s=0,hall_proposal_s=0,
         binding_checks_s=0,motion_admission_s=0,second_exposure_s=0,seal_s=0,audit_s=0;
  unsigned hall_iterations=0;
  hc::ProposalTimings hall;   // split of hall_proposal_s (motion, field, solve, heat)
};
inline constexpr int prepare_fault_points_v5=12;

// Clears nested JSON bottom-up so that the JSON destructor, which otherwise
// allocates a flattening stack, releases the audit without allocating.
inline void drain(Json& j) noexcept{
  if(j.is_object()){auto& o=j.get_ref<Json::object_t&>();for(auto& item:o)drain(item.second);o.clear();}
  else if(j.is_array()){auto& a=j.get_ref<Json::array_t&>();for(auto& item:a)drain(item);a.clear();}
}

class ClosedParentTrialV5 {
  friend ClosedParentTrialV5 prepare_closed_parent_joint_v5(m::PicArena&,ParentSubflow&,
    hc::Owner&,hc::MotionProvider&,const ClosedParentInputsV5&,double,int);
  friend bool can_publish_closed_parent_joint_v5(const m::PicArena&,const ParentSubflow&,
    const hc::Owner&,const ClosedParentTrialV5&) noexcept;
  friend bool try_publish_closed_parent_joint_v5(m::PicArena&,ParentSubflow&,
    hc::Owner&,ClosedParentTrialV5&) noexcept;
  friend struct TrialTestAccess;
  std::unique_ptr<m::PrivateFullTimeStep> private_step_;
  std::unique_ptr<m::PreparedArena> native_;
  PreparedParentState random_;
  std::unique_ptr<hc::JointProposal> field_heat_;
  std::vector<m::ArenaRateEntry> post_first_census_;
  ClosedParentInputsV5 inputs_;
  common_commit::ExpectedStep expected_;
  Json audit_;
  PhaseTimesV5 times_{};
  bool consumed_=false;
  // Test instrumentation. Default values leave behavior unchanged.
  int can_publish_fault_=0;
  mutable int refusal_point_=0;
  void (*publish_probe_)(int)=nullptr;
  std::unique_ptr<m::PreparedArena>* retain_native_=nullptr;
  std::unique_ptr<hc::JointProposal>* retain_field_heat_=nullptr;
  PreparedParentState* retain_random_=nullptr;
  ClosedParentTrialV5()=default;
public:
  ClosedParentTrialV5(ClosedParentTrialV5&&) noexcept=default;
  ClosedParentTrialV5& operator=(ClosedParentTrialV5&&) noexcept=default;
  ClosedParentTrialV5(const ClosedParentTrialV5&)=delete;
  ClosedParentTrialV5& operator=(const ClosedParentTrialV5&)=delete;
  ~ClosedParentTrialV5(){drain(audit_);}
  const Json& audit()const noexcept{return audit_;}
  const PhaseTimesV5& times()const noexcept{return times_;}
};

inline ClosedParentTrialV5 prepare_closed_parent_joint_v5(
    m::PicArena& arena,ParentSubflow& chemistry,hc::Owner& owner,hc::MotionProvider& provider,
    const ClosedParentInputsV5& inputs,double dt_s,int failure_point=0){
  using clock=std::chrono::steady_clock;
  pc::require(failure_point>=0&&failure_point<=prepare_fault_points_v5,"unknown manufactured preparation fault point");
  const auto fail=[&](int point){if(failure_point==point)throw std::runtime_error("manufactured v5 preparation fault");};
  const auto seconds=[](clock::time_point a){return std::chrono::duration<double>(clock::now()-a).count();};
  ClosedParentTrialV5 out;
  auto t=clock::now();
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
  {
    Encoding e;e.text("closed-parent-common-clock-v4");
    e.u64(common.accepted_step);e.u64(common.accepted_arena_revision);
    e.f64(common.begin_s);e.f64(common.end_s);e.f64(common.dt_s);
    for(const std::string* x:std::initializer_list<const std::string*>{&before.actual_native_prefix,&before.field,
        &before.heat,&before.source,&before.rng,&before.mesh,&before.material,&before.frame,&before.configuration,
        &inputs.species_sha256})e.text(*x);
    common.common_clock_sha256=e.finish();
  }
  out.inputs_=inputs;
  out.expected_={boundary,common.dt_s,common.end_s,common.common_clock_sha256,before.configuration};
  out.times_.observe_s=seconds(t);fail(1);
  t=clock::now();
  out.private_step_=std::make_unique<m::PrivateFullTimeStep>(
    arena.begin_full_time_step(common,inputs.maximum_carriers,inputs.maximum_native_records));
  out.times_.begin_s=seconds(t);fail(2);
  t=clock::now();
  m::ExposureBinding exposure{common,m::FullTimeSlot::before_transport,
    common.begin_s,.5*common.dt_s,first_sequence,before.field};
  auto first_scope=out.private_step_->first_exposure(exposure);
  auto first=chemistry.prepare_on(first_scope);
  out.times_.first_exposure_s=seconds(t);fail(3);
  t=clock::now();
  auto motion_scope=out.private_step_->fork_transport_trial();
  out.post_first_census_=motion_scope.native_census(inputs.maximum_carriers);
  require_census_species(out.post_first_census_,inputs.states);
  h::MotionBinding motion;
  motion.phase=h::ParticlePhase::full_time_position_and_proper_velocity;
  motion.accepted_step=common.accepted_step;motion.staged_arena_revision=motion_scope.staged_arena_revision();
  motion.begin_s=common.begin_s;motion.end_s=common.end_s;motion.dt_s=common.dt_s;
  motion.staged_arena_sha256=motion_scope.state_prefix();motion.mesh_sha256=before.mesh;
  motion.material_sha256=before.material;motion.species_sha256=inputs.species_sha256;
  motion.common_clock_sha256=common.common_clock_sha256;
  motion.kinetic_policy_sha256=h::kinetic_policy_sha256(h::declared_kinetic_policy(arena.native_options()));
  out.times_.census_s=seconds(t);fail(4);
  t=clock::now();
  out.field_heat_=std::make_unique<hc::JointProposal>(hc::prepare_closed_proposal(
    owner,arena.native_options(),out.post_first_census_,motion,before.source,provider));
  out.times_.hall_proposal_s=seconds(t);out.times_.hall_iterations=out.field_heat_->iterations();
  out.times_.hall=out.field_heat_->timings();fail(5);
  t=clock::now();
  const auto& proposal=*out.field_heat_;
  // The translation actually applied by the Hall motion provider must be the
  // derived frame, and motion must be bound to the authoritative midpoint
  // digest (the binary64 E/B digest, never the midpoint role seal).
  pc::require(proposal.motion().frame().native_to_field_m==native_to_field_shift(arena.native_options()),
    "Hall motion translation differs from the derived native/Hall frame");
  const auto midpoint_digest=h::total_field_sha256(proposal.midpoint().E(),proposal.midpoint().B());
  pc::require(proposal.motion().binding().trial_field_sha256==midpoint_digest&&
    midpoint_digest!=proposal.midpoint().seal(),"motion is not bound to the actual authoritative midpoint digest");
  const auto& field_after=proposal.candidate_field_sha256();
  const auto& heat_after=proposal.candidate_heat_sha256();
  out.times_.binding_checks_s=seconds(t);fail(6);
  t=clock::now();
  m::FullTimeMotionBinding native_motion{common,proposal.motion().seal_sha256(),
    before.field,field_after,heat_after,before.source};
  // Live admission comes only from the genuine immutable Hall object.
  auto moved=motion_scope.prepare_total_field_motion(native_motion,
    m::bind_hall_total_field_motion(proposal.motion()));
  fail(7);
  out.private_step_->accept_transport_trial(std::move(moved));
  out.times_.motion_admission_s=seconds(t);fail(8);
  t=clock::now();
  exposure.slot=m::FullTimeSlot::after_transport;exposure.frozen_physical_stage_s=common.end_s;
  exposure.owner_stage_sequence=first_sequence+1;exposure.frozen_field_sha256=field_after;
  auto second_scope=out.private_step_->second_exposure(exposure);
  auto second=chemistry.prepare_on(second_scope,&first);
  out.times_.second_exposure_s=seconds(t);fail(9);
  t=clock::now();
  out.native_=std::make_unique<m::PreparedArena>(out.private_step_->seal(
    {field_after,heat_after,before.source,chemistry.staged_rng_sha256(second)}));
  fail(10);
  out.random_=chemistry.seal(first,second,*out.native_);
  out.times_.seal_s=seconds(t);fail(11);
  t=clock::now();
  {
    Json& a=out.audit_;
    a=Json::object();
    a["schema"]="closed-parent-corrective-common-step-v5";
    a["begin_s"]=common.begin_s;a["end_s"]=common.end_s;a["authoritative_dt_s"]=common.dt_s;
    a["common_clock_sha256"]=common.common_clock_sha256;a["configuration_sha256"]=before.configuration;
    a["species_sha256"]=inputs.species_sha256;a["frame_sha256"]=before.frame;
    a["first_exposure"]=first.audit();a["second_exposure"]=second.audit();
    if(inputs.audit_transport_record)a["transport"]=h::motion_record(proposal.motion());
    else a["transport_seal_sha256"]=proposal.motion().seal_sha256();
    a["transport_record_v3_sha256"]=proposal.motion().record_v3_sha256();
    a["field_solve_backend"]=owner.model().solver().name();a["motion_provider"]=provider.name();
    a["motion_midpoint_sha256"]=midpoint_digest;a["midpoint_role_sha256"]=proposal.midpoint().seal();
    a["hall_iterations"]=proposal.iterations();a["field_heat_balance_J"]=proposal.combined_balance_J();
    a["native"]=out.native_->audit();
    a["candidate_field_sha256"]=field_after;a["candidate_heat_sha256"]=heat_after;
    a["surface_source_and_collection"]=false;a["conditional_parent_only"]=true;a["full_plasma_run"]=false;
  }
  out.times_.audit_s=seconds(t);fail(12);
  return out;
}

inline bool can_publish_closed_parent_joint_v5(const m::PicArena& arena,const ParentSubflow& chemistry,
    const hc::Owner& owner,const ClosedParentTrialV5& trial) noexcept{
  int point=0;
  const auto fault=[&]{++point;if(trial.can_publish_fault_==point)throw std::runtime_error("manufactured v5 can_publish fault");};
  const auto refuse=[&](int code){trial.refusal_point_=code;return false;};
  trial.refusal_point_=0;
  try{
    if(trial.consumed_||!trial.native_||!trial.field_heat_)return refuse(1);
    fault();
    if(!arena.can_publish(*trial.native_))return refuse(2);
    fault();
    if(!chemistry.can_publish(trial.random_,*trial.native_))return refuse(3);
    fault();
    const auto boundary=trial.native_->candidate_full_time_boundary();
    const auto exposures=trial.native_->candidate_full_time_exposures();
    const auto& proposal=*trial.field_heat_;
    if(!boundary||!exposures||exposures->transport_trial_sha256!=proposal.motion().seal_sha256())return refuse(4);
    fault();
    // Fresh observations: field/heat rehash, RNG payload recomputed, species
    // and frame rederived from the live arena and Hall owner.
    const auto current=observe_closed_accepted(arena,chemistry,owner,trial.inputs_);
    fault();
    common_commit::ProposedOwners candidate;auto& after=candidate.after;
    after.native=*boundary;after.actual_native_prefix=trial.native_->candidate_state_prefix();
    after.actual_native_revision=boundary->accepted_arena_revision;
    after.actual_native_step=boundary->accepted_step;
    after.field=proposal.candidate_field_sha256();after.heat=proposal.candidate_heat_sha256();
    after.source=current.source;after.rng=chemistry.candidate_rng_sha256(trial.random_,*trial.native_);
    after.mesh=current.mesh;after.material=current.material;after.frame=current.frame;
    after.configuration=current.configuration;
    after.field_time_s=after.heat_time_s=proposal.end_s();
    after.field_revision=after.heat_revision=proposal.endpoint_revision();
    fault();
    candidate.motion_midpoint=proposal.motion().binding().trial_field_sha256;
    candidate.field_midpoint=h::total_field_sha256(proposal.midpoint().E(),proposal.midpoint().B());
    if(candidate.field_midpoint==proposal.midpoint().seal())return refuse(5);
    fault();
    const auto& motion=proposal.motion().binding();
    if(proposal.motion().frame().native_to_field_m!=native_to_field_shift(arena.native_options())||
       motion.species_sha256!=trial.inputs_.species_sha256||motion.mesh_sha256!=current.mesh||
       motion.kinetic_policy_sha256!=h::kinetic_policy_sha256(h::declared_kinetic_policy(arena.native_options())))
      return refuse(6);
    fault();
    if(proposal.configuration_sha256()!=owner.configuration_sha256()||
       !common_commit::accepts_closed_step(trial.expected_,current,candidate))return refuse(7);
    fault();
    // Final Hall approval is last: all hashing, allocation and fallible checks
    // complete before the first accepted native/RNG/field/heat publication.
    if(!owner.preflight_final(proposal,trial.post_first_census_,*boundary))return refuse(8);
    return true;
  }catch(...){trial.refusal_point_=100+point;return false;}
}

// Publication order native -> chemistry RNG -> Hall field/heat, with no I/O or
// callbacks between the first and last swap. The published native token and
// the private step retain the outer reservation; both are released afterwards
// so that the next step can begin.
inline bool try_publish_closed_parent_joint_v5(m::PicArena& arena,ParentSubflow& chemistry,
    hc::Owner& owner,ClosedParentTrialV5& trial) noexcept{
  if(!can_publish_closed_parent_joint_v5(arena,chemistry,owner,trial))return false;
  if(trial.publish_probe_)trial.publish_probe_(0);
  if(!arena.try_publish(*trial.native_)){trial.refusal_point_=90;return false;}
  chemistry.publish_prevalidated(trial.random_);
  owner.publish_prevalidated(std::move(*trial.field_heat_));
  if(trial.publish_probe_)trial.publish_probe_(1);
  trial.consumed_=true;
  if(trial.retain_native_)*trial.retain_native_=std::move(trial.native_);else trial.native_.reset();
  if(trial.retain_field_heat_)*trial.retain_field_heat_=std::move(trial.field_heat_);else trial.field_heat_.reset();
  if(trial.retain_random_)*trial.retain_random_=std::move(trial.random_);
  trial.private_step_.reset();
  return true;
}
} // namespace raikou::native_water::v5
