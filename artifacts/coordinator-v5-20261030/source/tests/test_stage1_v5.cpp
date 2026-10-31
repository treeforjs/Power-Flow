// Stage 1: first real all-owner step (serial, one process, sealed sources
// unchanged). Modes:
//   t1     genesis + one published resident-only joint step (P1 and P2),
//          deterministic receipt, phase timers.
//   matrix T2-T7 (plus the in-process T1 reference used for T5).
#include "fixture_v5.hpp"
#include <algorithm>
#include <bit>
#include <random>
#include <optional>
#include <map>
#include <set>

namespace stage1 {
int forbidden_solve_calls = 0, forbidden_advance_calls = 0;
}
using namespace stage1;

// ---------------------------------------------------------------- forbidden legacy native calls
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic12solve_fieldsEv(m::GapPic*) {
  ++forbidden_solve_calls;
  throw std::runtime_error("forbidden full-time native solver trap");
}
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic7advanceEd(m::GapPic*, double) {
  ++forbidden_advance_calls;
  throw std::runtime_error("forbidden full-time native advance trap");
}
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic13advance_emptyEd(m::GapPic*, double) {
  ++forbidden_advance_calls;
  throw std::runtime_error("forbidden full-time native empty advance trap");
}
extern "C" void
__wrap__ZN6raikou12mfem_backend6GapPic16advance_observedEdRKSt8functionIFvRKNS0_11GapParticleEidEE(
    m::GapPic*, double, const m::GapPic::SinkObserver&) {
  ++forbidden_advance_calls;
  throw std::runtime_error("forbidden full-time observed advance trap");
}

// ---------------------------------------------------------------- publication probe
namespace stage1 {
long long probe_allocations[2]={-1,-1};
steady::time_point probe_times[2];
int probe_calls=0;
void publish_probe(int k){probe_allocations[k]=stage1_allocation_count();probe_times[k]=steady::now();++probe_calls;}
double seconds_between(steady::time_point a,steady::time_point b){return std::chrono::duration<double>(b-a).count();}

J publish_measured(v5::JointOwnersV5& f,v5::ClosedParentTrialV5& t){
  probe_allocations[0]=probe_allocations[1]=-1;probe_calls=0;
  TA::set_publish_probe(t,&publish_probe);
  const long long a0=stage1_allocation_count();const auto t0=steady::now();
  const bool ok=publish(f,t);
  const auto t1=steady::now();const long long a1=stage1_allocation_count();
  J r={{"published",ok},{"refusal_point",TA::refusal(t)},{"total_s",seconds_between(t0,t1)},{"probe_calls",probe_calls}};
  if(ok&&probe_calls==2){
    r["can_publish_s"]=seconds_between(t0,probe_times[0]);
    r["swap_window_s"]=seconds_between(probe_times[0],probe_times[1]);
    r["release_s"]=seconds_between(probe_times[1],t1);
    r["allocations_can_publish"]=probe_allocations[0]-a0;
    r["allocations_between_first_and_last_swap"]=probe_allocations[1]-probe_allocations[0];
    r["allocations_release"]=a1-probe_allocations[1];
  }
  return r;
}

J observed_json(const nw::common_commit::ObservedOwners& o){
  return {{"native_boundary",boundary_json(o.native)},{"actual_native_prefix",o.actual_native_prefix},
    {"field",o.field},{"heat",o.heat},{"source",o.source},{"rng",o.rng},{"mesh",o.mesh},
    {"material",o.material},{"frame",o.frame},{"configuration",o.configuration},
    {"native_revision",o.actual_native_revision},{"native_step",o.actual_native_step},
    {"field_revision",o.field_revision},{"field_time_s",o.field_time_s}};
}
std::size_t count_events(const J& audit){
  return audit.at("first_exposure").at("events").size()+audit.at("second_exposure").at("events").size();
}
J hall_json(const hc::ProposalTimings& h){
  return {{"admit_census_s",h.admit_census_s},{"motion_s",h.motion_s},{"deposit_support_s",h.support_s},
    {"field_prepare_s",h.field_prepare_s},{"field_solve_s",h.solve_s},{"independent_residual_check_s",h.residual_check_s},
    {"field_seal_s",h.seal_s},{"heat_s",h.heat_s},{"total_s",h.total_s},{"motion_v2_json_seal_s",h.motion_v2_seal_s},{"iterations",h.iterations},
    {"continuity_defect_over_tolerance",h.continuity_ratio},{"motion_provider",h.motion_provider},
    {"last_solve",{{"backend",h.last_solve.backend},{"iterations",h.last_solve.iterations},{"corrections",h.last_solve.corrections},
      {"independent_residual",h.last_solve.independent_residual},{"target",h.last_solve.target},{"backend_residual",h.last_solve.backend_residual}}}};
}
J times_json(const v5::PhaseTimesV5& t){
  return {{"observe_s",t.observe_s},{"begin_s",t.begin_s},{"first_exposure_s",t.first_exposure_s},
    {"census_s",t.census_s},{"hall_proposal_s",t.hall_proposal_s},{"binding_checks_s",t.binding_checks_s},
    {"motion_admission_s",t.motion_admission_s},{"second_exposure_s",t.second_exposure_s},
    {"seal_s",t.seal_s},{"audit_s",t.audit_s},{"hall_iterations",t.hall_iterations},{"hall",hall_json(t.hall)}};
}

// Timed replays of individual identity/validation operations on the accepted
// post-step state (the Hall proposal itself is split by its own timers).
J phase_probes(v5::JointOwnersV5& F,int repeats){
  const auto& model=*F.model;const auto o=F.arena->native_options();
  const auto census=F.arena->rate_population(F.inputs.maximum_carriers);
  const auto mid=model.at(F.owner->field());
  auto median=[&](auto&& fn){std::vector<double> v;for(int i=0;i<repeats;++i){const auto t=steady::now();fn();v.push_back(since(t));}
    std::sort(v.begin(),v.end());return v[v.size()/2];};
  J p;
  p["repeats"]=repeats;
  p["model_state_sha256_binary_s"]=median([&]{(void)model.state_sha256(F.owner->field(),F.owner->charge());});
  p["model_admit_constraints_s"]=median([&]{(void)model.admit(F.owner->field(),F.owner->charge());});
  p["heat_sha256_binary_s"]=median([&]{(void)hc::heat_sha256_of(model.material_sha256(),F.owner->thermal());});
  p["total_field_sha256_s"]=median([&]{(void)h::total_field_sha256(mid.E(),mid.B());});
  p["canonical_census_sha256_s"]=median([&]{(void)h::canonical_census_sha256(census);});
  p["motion_mesh_sha256_cached_s"]=median([&]{(void)model.motion_mesh_sha256(o);});
  p["observe_closed_accepted_s"]=median([&]{(void)v5::observe_closed_accepted(*F.arena,*F.chemistry,*F.owner,F.inputs);});
  p["chemistry_accepted_rng_sha256_s"]=median([&]{(void)F.chemistry->accepted_rng_sha256();});
  p["census_rows"]=census.size();
  return p;
}

// ---------------------------------------------------------------- T1
J run_t1_population(const Inputs& in,const std::string& label,std::vector<m::FullTimeAuthoredPacket> packets,
                    const fs::path& out,const fs::path& scratch,Results& R,J& timings,bool probes){
  const std::string T="T1";
  J rec={{"population",label}};J time;
  auto t=steady::now();
  auto F=genesis(in,std::move(packets));
  time["genesis_s"]=since(t);
  rec["genesis"]=F.genesis;
  R.record(T,label+": joint genesis admitted (absolute Gauss and divB against the actual census)",true,F.genesis.at("admitted_constraints"));
  const auto before=v5::observe_closed_accepted(*F.arena,*F.chemistry,*F.owner,F.inputs);
  rec["observed_before"]=observed_json(before);
  t=steady::now();
  rec["before"]=snapshot(F,scratch,label+"-before");
  time["owner_snapshot_s"]=since(t);
  const long long a0=stage1_allocation_count();
  t=steady::now();
  std::optional<v5::ClosedParentTrialV5> trial;
  try{trial.emplace(prepare(F));}
  catch(const std::exception& e){
    R.record(T,label+": prepare joint step",false,{{"error",e.what()}});
    timings[label]=time;return rec;
  }
  time["prepare_s"]=since(t);time["prepare_allocations"]=stage1_allocation_count()-a0;
  time["prepare_phases"]=times_json(trial->times());
  const auto& proposal=TA::proposal(*trial);
  const auto& audit=trial->audit();
  write_text(out/(label+"-trial-audit.json"),audit.dump(1)+"\n");
  rec["trial_audit_sha256"]=nw::sha(audit.dump());
  rec["hall_iterations"]=proposal.iterations();
  rec["combined_field_motion_heat_balance_J"]=proposal.combined_balance_J();
  const auto& budget=proposal.motion().budget();
  rec["motion_budget"]={{"kinetic_change_J",budget.kinetic_change_J},{"field_work_J",budget.field_work_J},
    {"work_defect_J",budget.work_defect_J},{"continuity_C",budget.continuity_C},{"plasma_dt",budget.plasma_dt},
    {"maximum_gyro",budget.maximum_gyro},{"maximum_displacement_cells",budget.maximum_displacement_cells}};
  const auto events=count_events(audit);
  rec["parent_events"]={{"first_exposure",audit.at("first_exposure").at("events").size()},
                        {"second_exposure",audit.at("second_exposure").at("events").size()}};
  const auto motion_digest=proposal.motion().binding().trial_field_sha256;
  const auto midpoint_digest=h::total_field_sha256(proposal.midpoint().E(),proposal.midpoint().B());
  rec["motion_midpoint_sha256"]=motion_digest;rec["midpoint_role_seal"]=proposal.midpoint().seal();
  R.record(T,label+": A10 motion digest equals the actual authoritative midpoint digest and differs from the role seal",
    motion_digest==midpoint_digest&&motion_digest!=proposal.midpoint().seal(),
    {{"motion",motion_digest},{"midpoint",midpoint_digest},{"role_seal",proposal.midpoint().seal()}});
  const axial::State old_field=F.owner->field();
  const std::vector<double> midE=proposal.midpoint().E(),midB=proposal.midpoint().B();
  const auto candidate_field=proposal.candidate_field_sha256(),candidate_heat=proposal.candidate_heat_sha256();
  rec["candidate_field_sha256"]=candidate_field;rec["candidate_heat_sha256"]=candidate_heat;
  const bool can=can_publish(F,*trial);
  R.record(T,label+": can_publish accepts the genuine sealed trial",can,{{"refusal_point",TA::refusal(*trial)}});
  const auto pub=publish_measured(F,*trial);
  time["publish"]=pub;
  J pub_det=pub;for(const char* k:{"total_s","can_publish_s","swap_window_s","release_s"})pub_det.erase(k);
  rec["publish"]=pub_det;
  R.record(T,label+": real positive resident-only joint step published",pub.at("published").get<bool>(),pub_det);
  if(!pub.at("published").get<bool>()){timings[label]=time;return rec;}
  R.record(T,label+": zero allocations between the first and last swap",
    pub.at("allocations_between_first_and_last_swap").get<long long>()==0,pub_det);
  const auto after=v5::observe_closed_accepted(*F.arena,*F.chemistry,*F.owner,F.inputs);
  rec["observed_after"]=observed_json(after);
  const double t1_expected=before.native.time_s+dt;
  R.record(T,label+": all owners advanced to t0+dt and step 1",
    after.native.time_s==t1_expected&&after.native.accepted_step==1&&F.owner->time_s()==t1_expected&&
    F.owner->revision()==1&&nw::common_commit::owner_values_match(after)&&!after.native.authored_full_time_genesis,
    {{"t0_s",before.native.time_s},{"t1_s",after.native.time_s},{"dt_s",dt}});
  R.record(T,label+": field, heat and RNG changed; source, mesh, material, frame and configuration unchanged",
    after.field!=before.field&&after.heat!=before.heat&&after.rng!=before.rng&&after.source==before.source&&
    after.mesh==before.mesh&&after.material==before.material&&after.frame==before.frame&&
    after.configuration==before.configuration&&after.field==candidate_field&&after.heat==candidate_heat,
    {{"field",{before.field,after.field}},{"heat",{before.heat,after.heat}},{"rng",{before.rng,after.rng}},{"source",after.source}});
  // Compensated midpoint versus the high-word average (A10 nonzero compensation).
  const auto& new_field=F.owner->field();
  long e_high=0,b_high=0,e_words=0,b_words=0,e_low_nonzero=0,b_low_nonzero=0;
  for(std::size_t e=0;e<midE.size();++e){
    if(midE[e]!=.5*(old_field.E[e]+new_field.E[e]))++e_high;
    if(midE[e]!=((hc::Words(old_field.E[e],old_field.E_low[e])+hc::Words(new_field.E[e],new_field.E_low[e]))*.5).rounded())++e_words;
    if(new_field.E_low[e]!=0)++e_low_nonzero;
  }
  for(std::size_t b=0;b<midB.size();++b){
    if(midB[b]!=.5*(old_field.B[b]+new_field.B[b]))++b_high;
    if(midB[b]!=((hc::Words(old_field.B[b],old_field.B_low[b])+hc::Words(new_field.B[b],new_field.B_low[b]))*.5).rounded())++b_words;
    if(new_field.B_low[b]!=0)++b_low_nonzero;
  }
  rec["midpoint_compensation"]={{"E_edges",midE.size()},{"B_faces",midB.size()},
    {"E_edges_compensated_differs_from_high_word_average",e_high},{"B_faces_compensated_differs_from_high_word_average",b_high},
    {"E_edges_not_equal_to_recomputed_compensated_average",e_words},{"B_faces_not_equal_to_recomputed_compensated_average",b_words},
    {"E_low_nonzero_after",e_low_nonzero},{"B_low_nonzero_after",b_low_nonzero}};
  R.record(T,label+": authoritative midpoint equals the recomputed compensated average; differing high-word edges recorded",
    e_words==0&&b_words==0,rec["midpoint_compensation"]);
  R.record(T,label+": zero calls to wrapped GapPic::solve_fields/advance/advance_empty/advance_observed",
    forbidden_solve_calls==0&&forbidden_advance_calls==0,{{"solve_fields",forbidden_solve_calls},{"advance_family",forbidden_advance_calls}});
  t=steady::now();
  const fs::path keep=out/(label+"-after-owner-bytes");fs::create_directories(keep);
  rec["after"]=snapshot(F,scratch,label+"-after",keep);
  time["owner_snapshot_after_s"]=since(t);
  R.record(T,label+": native save after a completed step",!rec["after"].contains("native_save_error"),
    {{"native_save_error",rec["after"].value("native_save_error","")}});
  // Outer reservation released: a second step can be prepared and dropped.
  t=steady::now();
  try{
    auto second=prepare(F);
    time["second_prepare_s"]=since(t);time["second_prepare_phases"]=times_json(second.times());
    const bool second_ok=can_publish(F,second);
    rec["second_step_audit_sha256"]=nw::sha(second.audit().dump());
    R.record(T,label+": outer reservation released; a second step prepares and preflights",second_ok,
      {{"refusal_point",TA::refusal(second)},{"second_step_begin_s",second.audit().at("begin_s")}});
  }catch(const std::exception& e){
    R.record(T,label+": outer reservation released; a second step prepares and preflights",false,{{"error",e.what()}});
  }
  const auto again=snapshot(F,scratch,label+"-after-drop");
  R.record(T,label+": dropping the prepared second step leaves every owner byte-identical",
    again.at("digest")==rec["after"].at("digest"),{{"after",rec["after"].at("digest")},{"after_drop",again.at("digest")}});
  rec["parent_event_count"]=events;
  if(probes){
    t=steady::now();time["phase_probes"]=phase_probes(F,5);time["phase_probes_wall_s"]=since(t);
  }
  timings[label]=time;
  return rec;
}

int run_t1_mode(const Inputs& in,const fs::path& out,double water_weight){
  Results R;J receipt,timings;
  const fs::path scratch=out/"scratch";fs::create_directories(scratch);
  // Plan position z = 0 (Hall node plane) first; then the Stage 1 positions
  // (z = 1 mm) for the bit-for-bit comparison with the Stage 1 receipts.
  receipt["P1"]=run_t1_population(in,"P1",packets_p1(centre_plan),out,scratch,R,timings,true);
  receipt["P1"]["position_m"]=centre_plan;
  receipt["P2"]=run_t1_population(in,"P2",packets_p2(water_weight,centre_plan),out,scratch,R,timings,false);
  receipt["P2"]["water_weight"]=water_weight;receipt["P2"]["position_m"]=centre_plan;
  receipt["P1-stage1-z1mm"]=run_t1_population(in,"P1-stage1-z1mm",packets_p1(centre_stage1),out,scratch,R,timings,false);
  receipt["P1-stage1-z1mm"]["position_m"]=centre_stage1;
  receipt["P2-stage1-z1mm"]=run_t1_population(in,"P2-stage1-z1mm",packets_p2(water_weight,centre_stage1),out,scratch,R,timings,false);
  receipt["P2-stage1-z1mm"]["water_weight"]=water_weight;receipt["P2-stage1-z1mm"]["position_m"]=centre_stage1;
  const bool p2_events=receipt["P2"].contains("parent_event_count")&&receipt["P2"]["parent_event_count"].get<std::size_t>()>=1;
  R.record("T1","P2 (z=0): at least one parent ionization event inside the published step",p2_events,
    receipt["P2"].value("parent_events",J::object()));
  R.record("T1","no forbidden native solver/advance call in the whole process",
    forbidden_solve_calls==0&&forbidden_advance_calls==0,{{"solve_fields",forbidden_solve_calls},{"advance_family",forbidden_advance_calls}});
  receipt["schema"]="stage1-v5-t1-deterministic-receipt-v1";
  receipt["backends"]={{"field",backends.field},{"push",backends.push},{"configuration_sha256",in.model->configuration_sha256()}};
  write_text(out/"t1-receipt.json",receipt.dump(1)+"\n");
  write_text(out/"t1-timings.json",timings.dump(1)+"\n");
  J result={{"schema","stage1-v5-t1-result-v1"},{"cases",R.cases},{"summary",R.summary()},
            {"receipt_sha256",nw::sha(receipt.dump(1)+"\n")}};
  write_text(out/"t1-result.json",result.dump(1)+"\n");
  fs::remove_all(scratch);
  for(const auto& c:R.cases)if(!c.at("passed").get<bool>())return 1;
  return 0;
}

// ---------------------------------------------------------------- matrix
J output_of(const J& child){
  if(child.contains("output")&&child.at("output").is_object())return child.at("output");
  return J::object();
}
struct Matrix {
  const Inputs& in;fs::path out,scratch;Results R;
  v5::JointOwnersV5 F;J S0;std::string S0d,Rd;J reference;
  long fresh_runs=0,fresh_matches=0,restored_runs=0,restored_matches=0;
  explicit Matrix(const Inputs& i,const fs::path& o):in(i),out(o),scratch(o/"scratch"),F(genesis(i,packets_p1())){
    fs::create_directories(scratch);
  }
  std::string digest(){return snapshot(F,scratch,"m").at("digest").get<std::string>();}
  bool unchanged(){return digest()==S0d;}
  // T5: in a child process, drop any outstanding trial, prepare a fresh one,
  // publish it and compare the post-step owner bytes with the T1 reference.
  J fresh(v5::ClosedParentTrialV5* outstanding){
    const auto c=in_child([&]{
      if(outstanding){auto dead=std::move(*outstanding);}
      auto t=prepare(F);const bool ok=publish(F,t);
      return J{{"published",ok},{"refusal_point",TA::refusal(t)},
               {"matches_reference",ok&&snapshot(F,scratch,"fresh").at("digest")==Rd}}.dump();
    });
    ++fresh_runs;J j=child_json(c);
    if(child_ok(j))++fresh_matches;
    return j;
  }
  // The restored (unmutated) trial still publishes, in a child process.
  J restored(v5::ClosedParentTrialV5& t){
    const auto c=in_child([&]{
      const bool ok=publish(F,t);
      return J{{"published",ok},{"refusal_point",TA::refusal(t)},
               {"matches_reference",ok&&snapshot(F,scratch,"restored").at("digest")==Rd}}.dump();
    });
    ++restored_runs;J j=child_json(c);
    if(child_ok(j))++restored_matches;
    return j;
  }
  static bool child_ok(const J& j){return j.value("exit_code",-1)==0&&output_of(j).value("matches_reference",false);}
};

void t6_genesis(Matrix& X){
  const std::string T="T6";
  X.R.record(T,"control: neutral P1 genesis admits with actual-census bootstrap",true,X.F.genesis.at("admitted_constraints"));
  auto attempt=[&](const std::string& name,v5::GenesisRequestV5 rq,v5::GenesisTestAccess::PacketHook ph,
                   v5::GenesisTestAccess::LedgerHook lh,int expected_stage,const std::string& expected_text){
    int reached=0;
    const auto msg=refusal_message([&]{auto owners=v5::GenesisTestAccess::admit(rq,*X.in.bridge,ph,lh,&reached);});
    X.R.record(T,name,!msg.empty()&&reached==expected_stage&&msg.find(expected_text)!=std::string::npos,
      {{"refusal",msg},{"stage_reached",reached},{"expected_stage",expected_stage},{"expected_guard",expected_text}});
  };
  {
    std::vector<m::FullTimeAuthoredPacket> lone{packets_p1().front()};
    attempt("lone electron with the DC seed refuses (absolute Gauss)",request(X.in,lone),nullptr,nullptr,4,
            "absolute Gauss constraint disagrees");
  }
  {
    auto rq=request(X.in,packets_p1());const auto& g=X.in.model->mesh();
    const int face=g.bx(22,22,g.kshort+1);rq.seed.B[face]+=.01;
    attempt("about 50 T/m divB seed refuses (single face +0.01 T over 200 um)",rq,nullptr,nullptr,4,"absolute divB constraint failed");
  }
  attempt("ledger gross-charge mismatch refuses",request(X.in,packets_p1()),nullptr,
    [](hc::ChargeLedger& q){for(auto& g:q.particles.gross_C)if(g>0){g*=2;break;}},4,
    "initial/restored net charge differs from actual canonical carriers");
  attempt("1-ulp ion move after the provisional hash refuses",request(X.in,packets_p1()),
    [](std::vector<m::FullTimeAuthoredPacket>& p){for(auto& x:p)if(x.id=="i0")x.position_m[0]=std::nextafter(x.position_m[0],1.);},
    nullptr,8,"actual native genesis census differs from the provisional Hall prediction");
  {
    auto rq=request(X.in,packets_p1());rq.seed.time_s=std::nextafter(rq.seed.time_s,INFINITY);
    attempt("genesis time different from the C19 drive start refuses",rq,nullptr,nullptr,1,
            "genesis time differs from the first supported C19 knot");
  }
}

m::FullTimeStepBinding raw_common(const v5::JointOwnersV5& f){
  const auto b=f.arena->full_time_boundary();
  return {b->accepted_step,b->accepted_arena_revision,b->time_s,b->time_s+dt,dt,b->accepted_arena_sha256,
          b->field_sha256,b->mesh_sha256,b->material_sha256,nw::sha("stage1/raw-native-reentry-probe-clock")};
}

void t7_entry(Matrix& X){
  const std::string T="T7";auto& F=X.F;
  int callbacks=0;
  const int solve0=forbidden_solve_calls,advance0=forbidden_advance_calls;
  J legacy;
  legacy["prepare_compound"]=refusal_message([&]{F.arena->prepare_compound(1,[&](m::PicArena&){++callbacks;});});
  legacy["prepare_compound_with_completion"]=refusal_message([&]{F.arena->prepare_compound_with_completion(1,
    [&](m::PicArena&){++callbacks;return m::ArenaStochasticCompletion{};});});
  legacy["prepare_transport"]=refusal_message([&]{F.arena->prepare_transport(dt);});
  legacy["prepare_event"]=refusal_message([&]{F.arena->prepare(m::ArenaEvent{});});
  {
    nw::Interval iv;iv.sequence=F.chemistry->checkpoint().at("next_interval").get<std::uint64_t>();
    iv.end_s=dt;iv.common_clock_sha256=nw::sha("stage1/legacy-parent-interval");
    iv.phase.state_prefix=F.arena->state_prefix();iv.phase.mesh_prefix=F.chemistry->grid_identity();
    iv.phase.clock_prefix=iv.common_clock_sha256;
    legacy["parent_legacy_prepare"]=refusal_message([&]{(void)F.chemistry->prepare(iv);});
  }
  bool all=true;for(const auto& [k,v]:legacy.items()){(void)k;all=all&&!v.get<std::string>().empty();}
  X.R.record(T,"A01: every legacy entry point refuses in full-time mode",all,legacy);
  X.R.record(T,"A01: legacy refusals ran zero callbacks and zero native solver/advance calls",
    callbacks==0&&forbidden_solve_calls==solve0&&forbidden_advance_calls==advance0,
    {{"callbacks",callbacks},{"solve_fields",forbidden_solve_calls-solve0},{"advance_family",forbidden_advance_calls-advance0}});
  X.R.record(T,"A01: legacy refusals leave every owner byte-identical",X.unchanged());
  {
    auto TF=prepare(F);J ev;
    ev["second_begin_full_time_step"]=refusal_message([&]{auto o=F.arena->begin_full_time_step(raw_common(F),8,32);});
    ev["sibling_v5_prepare"]=refusal_message([&]{auto o=prepare(F);});
    const auto common=raw_common(F);
    m::ExposureBinding e{common,m::FullTimeSlot::before_transport,common.begin_s,.5*dt,0,common.accepted_field_sha256};
    ev["sibling_first_exposure"]=refusal_message([&]{auto s=TA::private_step(TF)->first_exposure(e);});
    e.slot=m::FullTimeSlot::after_transport;e.frozen_physical_stage_s=common.end_s;
    ev["sibling_second_exposure"]=refusal_message([&]{auto s=TA::private_step(TF)->second_exposure(e);});
    ev["late_transport_fork"]=refusal_message([&]{auto s=TA::private_step(TF)->fork_transport_trial();});
    ev["late_reseal"]=refusal_message([&]{auto s=TA::private_step(TF)->seal({});});
    bool ok=true;for(const auto& [k,v]:ev.items()){(void)k;ok=ok&&!v.get<std::string>().empty();}
    X.R.record(T,"A03: second begin_full_time_step, sibling prepare, sibling exposure slots, late fork and reseal refuse",ok,ev);
    X.R.record(T,"A03: sibling refusals leave every owner byte-identical",X.unchanged());
    const auto pub=X.restored(TF);
    X.R.record(T,"A03: the outstanding trial still publishes after the sibling refusals (child process)",Matrix::child_ok(pub),pub);
  }
  {
    // Reentry from inside a private native collision builder.
    const auto common=raw_common(F);
    m::ExposureBinding e{common,m::FullTimeSlot::before_transport,common.begin_s,.5*dt,0,common.accepted_field_sha256};
    J ev;int inner=0;
    {
      auto outer=F.arena->begin_full_time_step(common,8,32);
      auto scope=outer.first_exposure(e);
      ev["builder"]=refusal_message([&]{(void)scope.prepare_compound_with_completion(1,[&](m::CollisionScope& s)->m::ArenaStochasticCompletion{
        ++inner;
        ev["v4_prepare_from_builder"]=refusal_message([&]{auto t=prepare(F);});
        ev["nested_scope_builder"]=refusal_message([&]{(void)s.prepare_compound_with_completion(1,[&](m::CollisionScope&){++callbacks;return m::ArenaStochasticCompletion{};});});
        ev["parent_prepare_on_from_builder"]=refusal_message([&]{(void)F.chemistry->prepare_on(s);});
        ev["outer_begin_from_builder"]=refusal_message([&]{auto o=F.arena->begin_full_time_step(common,8,32);});
        throw std::runtime_error("stage1 reentry probe ends the builder");});});
    }
    bool ok=inner==1;
    for(const char* k:{"v4_prepare_from_builder","nested_scope_builder","parent_prepare_on_from_builder","outer_begin_from_builder"})
      ok=ok&&ev.contains(k)&&!ev[k].get<std::string>().empty();
    ev["builder_callbacks"]=inner;ev["nested_callbacks"]=callbacks;
    X.R.record(T,"A03: reentry from inside a private collision builder refuses (v4 prepare, nested builder, parent prepare_on, outer begin)",
      ok&&callbacks==0,ev);
    X.R.record(T,"A03: reentry refusals leave every owner byte-identical",X.unchanged());
    const auto fresh=X.fresh(nullptr);
    X.R.record("T5","after T7 reentry refusals a fresh trial reproduces the T1 receipt (child process)",Matrix::child_ok(fresh),fresh);
  }
}

void t2_mutations(Matrix& X){
  const std::string T="T2";auto& F=X.F;
  auto Tw=genesis(X.in,packets_p1());
  auto geometry=request(X.in,packets_p1());
  geometry.options.lower_foil.left_m=std::nextafter(geometry.options.lower_foil.left_m,1.);
  geometry.options.upper_foil.left_m=geometry.options.lower_foil.left_m;
  std::optional<v5::JointOwnersV5> Gx;
  const auto gx_msg=refusal_message([&]{Gx.emplace(v5::JointGenesisV5::admit(geometry,*X.in.bridge));});
  X.R.record(T,"geometry twin (nextafter foil x coordinate) admits with a different mesh/frame identity",
    gx_msg.empty()&&Gx&&Gx->genesis.at("mesh_sha256")!=F.genesis.at("mesh_sha256")&&Gx->genesis.at("frame_sha256")!=F.genesis.at("frame_sha256"),
    {{"refusal",gx_msg},{"twin_mesh",Gx?Gx->genesis.at("mesh_sha256"):J()},{"twin_frame",Gx?Gx->genesis.at("frame_sha256"):J()}});
  X.R.record(T,"value-identical twin fixture has identical accepted bytes",snapshot(Tw,X.scratch,"twin").at("digest")==X.S0d);
  auto TT=prepare(Tw);
  std::optional<v5::ClosedParentTrialV5> GT;
  if(Gx)GT.emplace(prepare(*Gx));
  auto run_case=[&](const std::string& name,const std::function<void(v5::ClosedParentTrialV5&)>& mutate,
                    const std::function<void(v5::ClosedParentTrialV5&)>& restore,int expected_refusal){
    try{
      auto TF=prepare(F);
      mutate(TF);
      const bool published=publish(F,TF);
      const int refusal=TA::refusal(TF);
      const auto now=X.digest();
      J ev={{"published",published},{"refusal_point",refusal},{"owners_digest_before",X.S0d},{"owners_digest_after",now}};
      if(expected_refusal)ev["expected_refusal_point"]=expected_refusal;
      X.R.record(T,name+": refuses with all owners byte-identical",
        !published&&now==X.S0d&&(!expected_refusal||refusal==expected_refusal),ev);
      if(restore){
        restore(TF);
        const auto pub=X.restored(TF);
        X.R.record(T,name+": the unmutated trial still publishes afterwards (child process, T1 bytes)",Matrix::child_ok(pub),pub);
      }
      const auto fresh=X.fresh(&TF);
      X.R.record("T5","after T2 "+name+" a fresh trial reproduces the T1 receipt (child process)",Matrix::child_ok(fresh),fresh);
    }catch(const std::exception& e){X.R.record(T,name+": blocked",false,{{"error",e.what()}});}
  };
  using Trial=v5::ClosedParentTrialV5;
  run_case("native token from a value-identical twin arena",
    [&](Trial& t){std::swap(TA::native(t),TA::native(TT));},[&](Trial& t){std::swap(TA::native(t),TA::native(TT));},2);
  run_case("field/heat proposal from a value-identical twin Hall owner",
    [&](Trial& t){std::swap(TA::field_heat(t),TA::field_heat(TT));},[&](Trial& t){std::swap(TA::field_heat(t),TA::field_heat(TT));},0);
  run_case("RNG token from a twin chemistry owner",
    [&](Trial& t){std::swap(TA::random(t),TA::random(TT));},[&](Trial& t){std::swap(TA::random(t),TA::random(TT));},3);
  if(GT){
    run_case("geometry: native token from the nextafter-foil twin",
      [&](Trial& t){std::swap(TA::native(t),TA::native(*GT));},[&](Trial& t){std::swap(TA::native(t),TA::native(*GT));},2);
    run_case("geometry: Hall proposal from the nextafter-foil twin",
      [&](Trial& t){std::swap(TA::field_heat(t),TA::field_heat(*GT));},[&](Trial& t){std::swap(TA::field_heat(t),TA::field_heat(*GT));},0);
  }else X.R.record(T,"geometry twin swaps",false,{{"blocked","geometry twin genesis refused"}});
  std::unique_ptr<hc::JointProposal> held;
  run_case("heat: same-owner proposal for dt/2",
    [&](Trial& t){auto b=TA::proposal(t).motion().binding();b.dt_s=.5*dt;b.end_s=b.begin_s+b.dt_s;
      held=std::make_unique<hc::JointProposal>(hc::prepare_closed_proposal(*F.owner,F.arena->native_options(),TA::census(t),b,
        TA::expected(t).before.source_sha256));std::swap(TA::field_heat(t),held);},
    [&](Trial& t){std::swap(TA::field_heat(t),held);held.reset();},0);
  run_case("source: same-owner proposal bound to a foreign source digest (refused by preflight_final)",
    [&](Trial& t){const auto b=TA::proposal(t).motion().binding();
      held=std::make_unique<hc::JointProposal>(hc::prepare_closed_proposal(*F.owner,F.arena->native_options(),TA::census(t),b,
        nw::sha("stage1/foreign-source")));std::swap(TA::field_heat(t),held);},
    [&](Trial& t){std::swap(TA::field_heat(t),held);held.reset();},8);
  std::string saved;
  run_case("config: species identity",
    [&](Trial& t){saved=TA::inputs(t).species_sha256;TA::inputs(t).species_sha256=nw::sha("stage1/other-species");},
    [&](Trial& t){TA::inputs(t).species_sha256=saved;},0);
  run_case("config: frame identity",
    [&](Trial& t){saved=TA::inputs(t).frame_sha256;TA::inputs(t).frame_sha256=nw::sha("stage1/other-frame");},
    [&](Trial& t){TA::inputs(t).frame_sha256=saved;},0);
  run_case("config: carrier cap",
    [&](Trial& t){--TA::inputs(t).maximum_carriers;},[&](Trial& t){++TA::inputs(t).maximum_carriers;},0);
  run_case("config: native record cap",
    [&](Trial& t){--TA::inputs(t).maximum_native_records;},[&](Trial& t){++TA::inputs(t).maximum_native_records;},0);
  // Advance the twin by one step, then use its step-2 proposal (other clock).
  {
    const bool twin=publish(Tw,TT);
    X.R.record(T,"twin fixture publishes its own step (control)",twin,{{"refusal_point",TA::refusal(TT)}});
    std::optional<Trial> TT2;
    if(twin)TT2.emplace(prepare(Tw));
    if(TT2)run_case("field: proposal from another clock (twin at step 1)",
      [&](Trial& t){std::swap(TA::field_heat(t),TA::field_heat(*TT2));},[&](Trial& t){std::swap(TA::field_heat(t),TA::field_heat(*TT2));},0);
  }
  // Last: chemistry restore between prepare and publish (stale RNG token).
  run_case("RNG: chemistry restore() between prepare and publish",
    [&](Trial&){F.chemistry->restore(F.chemistry->checkpoint());},nullptr,3);
}

// Classifies a terminated sweep child from its symbolized backtrace: whether
// the injected failure was raised inside a nlohmann::basic_json destructor
// (its flattening stack allocates under noexcept) and which component's frame
// is the innermost non-library caller.
// Owner of the innermost non-library frame of a symbolized stack (same table as crashes).
std::string attribute_frames(const std::string& text);
J attribute_crash(const Child& c){
  J a={{"exit_code",c.code},{"signal",c.signal}};
  const auto at=c.out.find("@@terminate-backtrace@@");
  if(at==std::string::npos){a["kind"]="no-backtrace";a["component"]="unknown";return a;}
  std::vector<std::string> frames;std::string line;
  for(char ch:c.out.substr(at)){if(ch=='\n'){if(!line.empty())frames.push_back(line);line.clear();}else line+=ch;}
  if(!line.empty())frames.push_back(line);
  bool json_destructor=false;
  for(const auto& f:frames)if(f.find("basic_json")!=std::string::npos&&(f.find("destroy")!=std::string::npos||f.find("D1Ev")!=std::string::npos||f.find("D2Ev")!=std::string::npos))json_destructor=true;
  a["kind"]=json_destructor?"json-destructor-allocation-under-noexcept":"other-noexcept-allocation";
  // Most specific first. Hall v2 code: namespace hall_corrective_v2 and the v2
  // motion provider in hall_driver (motion_provider.cpp); the sealed kinetic
  // boundary and thermal functions in hall_driver are listed before it.
  static const std::pair<const char*,const char*> owners[]={{"12native_water2v5","v5-coordinator"},
    {"18hall_corrective_v2","hall-v2"},{"11hall_driver7Thermal","hall-sealed-thermal"},
    {"11hall_driver17first_kinetic_hit","hall-sealed-kinetic-boundary"},{"11hall_driver24require_kinetic_interior","hall-sealed-kinetic-boundary"},
    {"11hall_driver23declared_kinetic_policy","hall-sealed-kinetic-boundary"},{"11hall_driver21kinetic_policy_sha256","hall-sealed-kinetic-boundary"},
    {"11hall_driver22prepare_neutral_stream","hall-sealed-kinetic-boundary"},{"11hall_driver","hall-v2-motion"},
    {"5axial","hall-sealed-field-headers"},{"4foil","hall-sealed-field-headers"},{"17electrode_heat_v2","hall-sealed-thermal"},
    {"8hall_gpu","gpu-common-push"},{"12mfem_backend","native"},{"12native_water","v3-parent-subflow"},{"14pic_collisions","v3-collision-scheduler"},
    {"9pic_water","water-bridge"},{"6stage1","stage1-test"}};
  std::string component="unattributed";J innermost=J::array();
  for(std::size_t i=1;i<frames.size()&&innermost.size()<8;++i)innermost.push_back(frames[i].substr(0,220));
  for(std::size_t i=1;i<frames.size()&&component=="unattributed";++i){
    if(frames[i].find("report_terminate")!=std::string::npos||frames[i].find("basic_json")!=std::string::npos)continue;
    for(const auto& [key,name]:owners)if(frames[i].find(key)!=std::string::npos){component=name;break;}
  }
  a["component"]=component;a["innermost_frames"]=innermost;
  return a;
}
std::string attribute_frames(const std::string& text){
  Child c;c.out="@@terminate-backtrace@@\n"+text;
  return attribute_crash(c).at("component").get<std::string>();
}

void t3_faults(Matrix& X,double sweep_budget_s){
  const std::string T="T3";auto& F=X.F;
  for(int k=1;k<=v5::prepare_fault_points_v5;++k){
    const auto msg=refusal_message([&]{auto t=prepare(F,k);});
    const bool same=X.unchanged();
    X.R.record(T,"prepare throw point "+std::to_string(k)+" refuses with all owners byte-identical",
      msg.find("manufactured v5 preparation fault")!=std::string::npos&&same,{{"refusal",msg}});
    const auto fresh=X.fresh(nullptr);
    X.R.record("T5","after prepare throw point "+std::to_string(k)+" a fresh trial reproduces the T1 receipt (child process)",
      Matrix::child_ok(fresh),fresh);
  }
  for(int k=1;k<=9;++k){
    auto TF=prepare(F);TA::set_can_publish_fault(TF,k);
    const bool published=publish(F,TF);const int refusal=TA::refusal(TF);
    X.R.record(T,"can_publish throw point "+std::to_string(k)+" refuses with all owners byte-identical",
      !published&&refusal==100+k&&X.unchanged(),{{"refusal_point",refusal}});
    TA::set_can_publish_fault(TF,0);
    const auto pub=X.restored(TF);
    X.R.record(T,"after can_publish throw point "+std::to_string(k)+" the same trial publishes (child process, T1 bytes)",Matrix::child_ok(pub),pub);
    const auto fresh=X.fresh(&TF);
    X.R.record("T5","after can_publish throw point "+std::to_string(k)+" a fresh trial reproduces the T1 receipt (child process)",
      Matrix::child_ok(fresh),fresh);
  }
  {
    // Fault index 10 is past the last throw point: the trial must publish.
    auto TF=prepare(F);TA::set_can_publish_fault(TF,10);
    const auto pub=X.restored(TF);
    X.R.record(T,"can_publish has exactly nine throw points before the final Hall preflight (index 10 publishes, child process)",
      Matrix::child_ok(pub),pub);
  }
  // Allocation census of one complete step, in a child.
  const auto census=in_child([&]{
    const long long a0=stage1_allocation_count();auto t=prepare(F);const long long a1=stage1_allocation_count();
    J p=publish_measured(F,t);p["prepare_allocations"]=a1-a0;
    p["matches_reference"]=p.value("published",false)&&snapshot(F,X.scratch,"alloc-census").at("digest")==X.Rd;
    return p.dump();});
  J cj=child_json(census);const J co=output_of(cj);
  X.R.record(T,"allocation census step publishes the T1 bytes with zero allocations between the first and last swap",
    census.code==0&&co.value("matches_reference",false)&&co.value("allocations_between_first_and_last_swap",-1LL)==0,cj);
  if(census.code!=0||!co.value("published",false))return;
  const long long Ap=co.at("prepare_allocations").get<long long>();
  const long long Ac=co.at("allocations_can_publish").get<long long>();
  const long long As=co.at("allocations_between_first_and_last_swap").get<long long>();
  const long long Ar=co.at("allocations_release").get<long long>();
  const long long total=Ap+Ac+As+Ar;
  std::vector<long long> order;std::set<long long> seen;
  auto add=[&](long long k){if(k>=0&&k<=total+1&&seen.insert(k).second)order.push_back(k);};
  for(long long k=0;k<64;++k)add(k);
  for(long long k=Ap-16;k<Ap+16;++k)add(k);
  for(long long k=Ap+Ac-16;k<=total+1;++k)add(k);
  for(long long stride=std::max(1LL,total/64);;stride/=2){for(long long k=0;k<total;k+=stride)add(k);if(stride<=1)break;}
  const auto start=steady::now();
  J points=J::array();long executed=0,clean=0,crashes=0,published_points=0,reproduced=0,reproduced_ok=0;
  std::map<std::string,long> crash_kinds,silent;
  for(long long k:order){
    if(since(start)>sweep_budget_s)break;
    const bool reproduce=executed<10||executed%25==0;
    const auto c=in_child([&]{
      J r={{"k",k}};int phase=1;bool published=false;std::string error;
      stage1_injection_report_fd=child_report_fd;
      stage1_set_allocation_countdown(k);
      try{auto t=prepare(F);phase=2;published=publish(F,t);phase=3;}
      catch(const std::bad_alloc&){error="std::bad_alloc";}
      catch(const std::exception& e){error=e.what();}
      stage1_set_allocation_countdown(-1);stage1_injection_report_fd=-1;
      r["phase_reached"]=phase;r["published"]=published;r["error"]=error;
      r["injected"]=stage1_injected_allocation_failures()>0;
      if(published)r["matches_reference"]=snapshot(F,X.scratch,"sweep").at("digest")==X.Rd;
      else{
        r["owners_identical"]=snapshot(F,X.scratch,"sweep").at("digest")==X.S0d;
        if(reproduce){auto t2=prepare(F);const bool ok=publish(F,t2);
          r["fresh_matches_reference"]=ok&&snapshot(F,X.scratch,"sweep-fresh").at("digest")==X.Rd;}
      }
      return r.dump();});
    ++executed;J j=child_json(c);j["k"]=k;
    j["region"]=k<Ap?"prepare":k<Ap+Ac?"can_publish":k<Ap+Ac+As?"between_swaps":k<total?"release":"beyond_last_allocation";
    bool ok=false;
    if(j.contains("injection_backtrace")){j["injection_component"]=attribute_frames(j.at("injection_backtrace").get<std::string>());}
    if(c.code==0&&j.contains("output")&&j.at("output").is_object()){
      const J o=j.at("output");
      if(o.value("published",false)){++published_points;ok=o.value("matches_reference",false)&&!o.value("injected",true);
        if(!ok&&o.value("injected",false)){++silent[j.value("injection_component",std::string("unattributed"))];}}
      else ok=o.value("injected",false)&&o.value("owners_identical",false);
      if(o.contains("fresh_matches_reference")){++reproduced;if(o["fresh_matches_reference"].get<bool>())++reproduced_ok;}
    }else{
      ++crashes;
      const auto attribution=attribute_crash(c);
      j["crash"]=attribution;
      ++crash_kinds[attribution.at("kind").get<std::string>()+" | "+attribution.at("component").get<std::string>()];
      j.erase("output_text");
    }
    if(ok)++clean;
    j["clean"]=ok;
    if(!ok||executed<=20||executed%50==0)points.push_back(j);
  }
  long v2_crashes=0,silent_total=0,silent_v2=0;
  for(const auto& [k,v]:crash_kinds)if(k.find("| hall-v2")!=std::string::npos)v2_crashes+=v;
  for(const auto& [k,v]:silent){silent_total+=v;if(k.rfind("hall-v2",0)==0)silent_v2+=v;}
  J sweep={{"crashes_in_hall_v2_code",v2_crashes},{"silent_publications_after_injected_failure",silent_total},
    {"silent_publications_by_injection_site",silent},{"silent_publications_injected_in_hall_v2_code",silent_v2},{"allocations_prepare",Ap},{"allocations_can_publish",Ac},{"allocations_between_swaps",As},
    {"allocations_release",Ar},{"total_allocation_sites_in_one_step",total},{"schedule_points",order.size()},
    {"executed_points",executed},{"clean_points",clean},{"crashed_points",crashes},
    {"points_beyond_last_allocation_published",published_points},{"fresh_reproductions",reproduced},
    {"fresh_reproductions_matching",reproduced_ok},{"budget_s",sweep_budget_s},{"elapsed_s",since(start)},
    {"crash_attribution",crash_kinds},
    {"recorded_points",points}};
  write_text(X.out/"allocation-sweep.json",sweep.dump(1)+"\n");
  X.R.record(T,"global operator new countdown sweep: no std::terminate attributed to Hall v2 code (hall_corrective_v2, v2 motion provider)",
    executed>0&&v2_crashes==0,{{"executed",executed},{"crashed",crashes},{"crashes_in_hall_v2_code",v2_crashes},{"crash_attribution",crash_kinds}});
  X.R.record(T,"global operator new countdown sweep: no injected failure inside Hall v2 code is swallowed into a publication",
    executed>0&&silent_v2==0,{{"silent_publications",silent_total},{"by_injection_site",silent}});
  long sealed_crashes=0;
  for(const auto& [k,v]:crash_kinds)for(const char* c:{"| native","| v3-parent-subflow","| v3-collision-scheduler","| water-bridge",
      "| hall-sealed-thermal","| hall-sealed-kinetic-boundary","| hall-sealed-field-headers"})if(k.find(c)!=std::string::npos)sealed_crashes+=v;
  long sealed_silent=0;
  for(const auto& [k,v]:silent)for(const char* c:{"native","v3-parent-subflow","v3-collision-scheduler","water-bridge",
      "hall-sealed-thermal","hall-sealed-kinetic-boundary","hall-sealed-field-headers"})if(k==c)sealed_silent+=v;
  X.R.record(T,"global operator new countdown sweep: every injected failure refuses with all owners byte-identical, except crashes and "
    "silent publications whose injection site is in a sealed component (open sealed-code defects: Stage 1 defect 2; swallowed failures)",
    executed>0&&clean+crashes+silent_total==executed&&sealed_crashes==crashes&&sealed_silent==silent_total,{{"executed",executed},{"clean",clean},{"crashed",crashes},{"total_allocation_sites",total},
      {"crash_attribution",crash_kinds},
      {"coverage_fraction",total>0?double(executed)/double(total+2):1.0},{"detail","allocation-sweep.json"}});
  X.R.record("T5","sweep subset: fresh trials after injected allocation failures reproduce the T1 receipt",
    reproduced>0&&reproduced==reproduced_ok,{{"reproduced",reproduced},{"matching",reproduced_ok}});
}

void t4_single_use(Matrix& X){
  const std::string T="T4";auto& F=X.F;using Trial=v5::ClosedParentTrialV5;
  auto A=prepare(F);
  Trial moved=std::move(A);
  {const bool refused=!publish(F,A);const int point=TA::refusal(A);
   X.R.record(T,"moved-from trial refuses",refused&&point==1&&X.unchanged(),{{"refusal_point",point}});}
  std::unique_ptr<m::PreparedArena> held_native;std::unique_ptr<hc::JointProposal> held_hall;nw::PreparedParentState held_random;
  TA::retain(moved,&held_native,&held_hall,&held_random);
  const bool first=publish(F,moved);
  const auto S1=snapshot(F,X.scratch,"t4-step1");
  X.R.record(T,"the moved trial publishes step 1 with the T1 bytes",first&&S1.at("digest")==X.Rd,
    {{"digest",S1.at("digest")},{"reference",X.Rd}});
  const std::string S1d=S1.at("digest").get<std::string>();
  auto same=[&]{return snapshot(F,X.scratch,"t4").at("digest")==S1d;};
  {const bool refused=!publish(F,moved);const int point=TA::refusal(moved);
   X.R.record(T,"second try_publish of the published trial refuses",refused&&point==1&&same(),{{"refusal_point",point}});}
  // While the published native token is retained the outer reservation is
  // still held: the next step cannot begin. v4 releases it after publication.
  const auto blocked=refusal_message([&]{auto t=prepare(F);});
  X.R.record(T,"a retained published native token keeps the outer reservation (next begin refuses)",!blocked.empty(),{{"refusal",blocked}});
  {
    Trial shell=TA::shell();TA::copy_contract(shell,moved);
    TA::native(shell)=std::move(held_native);TA::field_heat(shell)=std::move(held_hall);TA::random(shell)=std::move(held_random);
    const bool arena_accepts=F.arena->can_publish(*TA::native(shell));
    const bool reused=publish(F,shell);const int point=TA::refusal(shell);
    X.R.record(T,"reused (already published) native token refuses",!arena_accepts&&!reused&&point==2&&same(),
      {{"native_can_publish",arena_accepts},{"refusal_point",point}});
    held_native=std::move(TA::native(shell));held_hall=std::move(TA::field_heat(shell));held_random=std::move(TA::random(shell));
  }
  held_native.reset();
  std::optional<Trial> B;
  const auto b_msg=refusal_message([&]{B.emplace(prepare(F));});
  X.R.record(T,"after releasing the token, step 2 prepares",b_msg.empty()&&B.has_value(),{{"refusal",b_msg}});
  if(!B)return;
  const auto census=TA::census(*B);
  X.R.record(T,"consumed Hall proposal fails Hall preflight directly",!F.owner->preflight(*held_hall,census));
  std::swap(TA::field_heat(*B),held_hall);
  const bool hall=publish(F,*B);const int hall_refusal=TA::refusal(*B);
  std::swap(TA::field_heat(*B),held_hall);
  X.R.record(T,"consumed Hall proposal refuses inside a step-2 trial",!hall&&same(),{{"refusal_point",hall_refusal}});
  std::swap(TA::random(*B),held_random);
  const bool rng=publish(F,*B);const int rng_refusal=TA::refusal(*B);
  std::swap(TA::random(*B),held_random);
  X.R.record(T,"stale chemistry generation (consumed step-1 RNG token) refuses",!rng&&rng_refusal==3&&same(),{{"refusal_point",rng_refusal}});
  const bool second=publish(F,*B);
  X.R.record(T,"the unmutated step-2 trial publishes (two consecutive published steps)",
    second&&F.arena->full_time_boundary()->accepted_step==2&&F.owner->revision()==2,
    {{"accepted_step",F.arena->full_time_boundary()->accepted_step},{"time_s",F.owner->time_s()}});
}

// Diagnostic: allocation-failure points K0..K1 (one child process each, fresh
// genesis); reports every point where an injected failure still published and
// which owner's bytes then differ from an uninjected reference step.
int run_sweep_point_mode(const Inputs& in,const fs::path& out,long long k0,long long k1){
  const fs::path scratch=out/"scratch";fs::create_directories(scratch);
  J ref;{auto F=genesis(in,packets_p1());auto t=prepare(F);if(!publish(F,t))throw std::runtime_error("reference step refused");ref=snapshot(F,scratch,"ref");}
  J points=J::array();long anomalies=0;
  for(long long k=k0;k<=k1;++k){
    const auto c=in_child([&]{
      auto F=genesis(in,packets_p1());
      J r={{"k",k}};std::string error;bool published=false;
      stage1_set_allocation_countdown(k);
      try{auto t=prepare(F);published=publish(F,t);}catch(const std::exception& e){error=e.what();}
      stage1_set_allocation_countdown(-1);
      r["injected"]=stage1_injected_allocation_failures()>0;r["published"]=published;r["error"]=error;
      const auto now=snapshot(F,scratch,"point");
      J diff=J::object();for(const auto& [key,v]:ref.items())if(key!="digest"&&(!now.contains(key)||now.at(key)!=v))diff[key]={{"reference",v},{"point",now.value(key,J())}};
      r["differing_snapshot_fields"]=diff;return r.dump();});
    J j=child_json(c);const J o=output_of(j);
    const bool anomaly=o.value("injected",false)&&o.value("published",false)&&!o.value("differing_snapshot_fields",J::object()).empty();
    if(anomaly){++anomalies;j["anomaly"]=true;}
    if(anomaly||c.code!=0)points.push_back(j);
  }
  J r={{"k0",k0},{"k1",k1},{"anomalies",anomalies},{"points",points}};
  write_text(out/"sweep-points.json",r.dump(1)+"\n");std::cout<<r.dump(1).substr(0,20000)<<std::endl;fs::remove_all(scratch);return 0;
}

int run_matrix_mode(const Inputs& in,const fs::path& out,double sweep_budget_s){
  auto start=steady::now();
  Matrix X(in,out);J stage;
  X.S0=snapshot(X.F,X.scratch,"S0");X.S0d=X.S0.at("digest").get<std::string>();
  {
    const auto c=in_child([&]{auto t=prepare(X.F);J p=publish_measured(X.F,t);
      p["snapshot"]=snapshot(X.F,X.scratch,"reference");return p.dump();});
    const J j=child_json(c);const J o=output_of(j);
    const bool ok=c.code==0&&o.value("published",false)&&o.contains("snapshot");
    if(ok){X.reference=o.at("snapshot");X.Rd=X.reference.at("digest").get<std::string>();}
    X.R.record("T1","matrix reference: fresh positive joint step publishes (child process)",ok,
      {{"exit_code",c.code},{"signal",c.signal},{"reference_digest",X.Rd}});
    if(!ok){write_text(out/"matrix-result.json",J{{"cases",X.R.cases},{"fatal","reference step did not publish"}}.dump(1)+"\n");return 1;}
  }
  auto group=[&](const std::string& name,const std::function<void()>& fn){
    const auto t=steady::now();
    try{fn();}catch(const std::exception& e){X.R.record(name,"group blocked by an exception",false,{{"error",e.what()}});}
    stage[name+"_s"]=since(t);
  };
  group("T6",[&]{t6_genesis(X);});
  group("T7",[&]{t7_entry(X);});
  group("T2",[&]{t2_mutations(X);});
  group("T3",[&]{t3_faults(X,sweep_budget_s);});
  X.R.record("T5","failure atomicity tally: all fresh trials after refusals reproduce the T1 receipt",
    X.fresh_runs>0&&X.fresh_runs==X.fresh_matches,{{"fresh_runs",X.fresh_runs},{"matching",X.fresh_matches},
      {"restored_runs",X.restored_runs},{"restored_matching",X.restored_matches}});
  group("T4",[&]{t4_single_use(X);});
  X.R.record("T7","no forbidden native solver/advance call in the whole matrix process",
    forbidden_solve_calls==0&&forbidden_advance_calls==0,{{"solve_fields",forbidden_solve_calls},{"advance_family",forbidden_advance_calls}});
  J result={{"schema","stage1-matrix-result-v1"},{"cases",X.R.cases},{"summary",X.R.summary()},
    {"S0",X.S0},{"reference_after_step",X.reference},{"elapsed_s",since(start)},{"group_wall_s",stage}};
  write_text(out/"matrix-result.json",result.dump(1)+"\n");
  std::error_code ignored;fs::remove_all(X.scratch,ignored);
  for(const auto& c:X.R.cases)if(!c.at("passed").get<bool>())return 1;
  return 0;
}
// ================================================================ v5 modes: GPU step, multi-step, Hall-only loop
// Runs the device provider and the CPU provider on the same census and midpoint
// for every fixed-point iterate, compares record bytes, seal and endpoint free
// charge bit for bit, and returns the device result.
class ParityMotionProvider final:public hc::MotionProvider {
  std::shared_ptr<hc::MotionProvider> device_;hc::CpuMotionProvider cpu_;std::string name_;
public:
  long compared=0,identical=0;J details=J::array();
  explicit ParityMotionProvider(std::shared_ptr<hc::MotionProvider> d):device_(std::move(d)),name_(device_->name()+"+cpu-parity"){}
  const std::string& name()const noexcept override{return name_;}
  void begin(const hc::MotionContext& c,std::span<const m::ArenaRateEntry> census)override{device_->begin(c,census);cpu_.begin(c,census);}
  h::MotionReceipt motion(const h::MotionBinding& b,const hc::Vec& E,const hc::Vec& B,hc::FreeCharge& after)override{
    auto t=steady::now();auto g=device_->motion(b,E,B,after);const double tg=since(t);
    t=steady::now();hc::FreeCharge ca;auto c=cpu_.motion(b,E,B,ca);const double tc=since(t);
    // Binary motion v3 record digest covers every record field (binding, census, policy,
    // frame, endpoints, paths, array digests, budget); the arrays are compared directly.
    bool words=ca.net_C.size()==after.net_C.size()&&ca.gross_C==after.gross_C;
    for(std::size_t k=0;words&&k<ca.net_C.size();++k)
      words=std::bit_cast<std::uint64_t>(ca.net_C[k].hi)==std::bit_cast<std::uint64_t>(after.net_C[k].hi)&&
            std::bit_cast<std::uint64_t>(ca.net_C[k].lo)==std::bit_cast<std::uint64_t>(after.net_C[k].lo);
    const bool v3=g.record_v3_sha256()==c.record_v3_sha256();
    const bool arrays=g.current_A_m()==c.current_A_m()&&g.initial_charge_C()==c.initial_charge_C()&&g.endpoint_charge_C()==c.endpoint_charge_C();
    const bool same=v3&&words&&arrays;
    ++compared;if(same)++identical;
    if(details.size()<64)details.push_back({{"carriers",g.endpoints().size()},{"record_v3_identical",v3},
      {"endpoint_free_charge_words_identical",words},{"current_and_charge_arrays_identical",arrays},{"record_v3_sha256",g.record_v3_sha256()},
      {"device_s",tg},{"cpu_s",tc}});
    return g;
  }
  void end()noexcept override{device_->end();cpu_.end();}
};

J guard_receipt(const hc::JointProposal& p){
  const auto& fh=p.field_heat();const auto& r=fh.field_receipt();const auto& c=fh.candidate_constraints();const auto& s=fh.solve_record();
  return {{"backend",s.backend},{"independent_residual",s.independent_residual},{"residual_target",s.target},
    {"residual_over_target",s.target>0?s.independent_residual/s.target:0.},{"backend_residual",s.backend_residual},{"solver_iterations",s.iterations},
    {"max_gauss_C",c.max_gauss_C},{"max_gauss_bound_C",c.max_gauss_bound_C},{"max_divB_T_m",c.max_divB_T_m},{"max_divB_bound_T_m",c.max_divB_bound_T_m},
    {"gauss_increment_C",r.gauss_increment_C},{"field_energy_balance_J",r.balance_J},{"field_heat_balance_J",fh.field_heat_balance_J()},
    {"combined_field_motion_heat_balance_J",p.combined_balance_J()},{"fixed_point_iterations",p.iterations()},
    {"candidate_field_sha256",p.candidate_field_sha256()},{"candidate_heat_sha256",p.candidate_heat_sha256()},
    {"configuration_sha256",p.configuration_sha256()}};
}
long double sum_words(const hc::WordVec& v){long double s=0;for(const auto& w:v)s+=w.wide();return s;}
long double sum_vec(const hc::Vec& v){long double s=0;for(double x:v)s+=x;return s;}
// Per-step ledger values available before publication (proposal) and after (owners).
J step_ledger(const hc::JointProposal& p){
  const auto& fh=p.field_heat();const auto& r=fh.field_receipt();const auto& c=fh.candidate_constraints();const auto& q=fh.candidate_charge();
  const auto& b=p.motion().budget();
  return {{"charge",{{"free_net_C",double(sum_words(q.particles.net_C))},{"free_gross_C",double(sum_vec(q.particles.gross_C))},
      {"conductor_bank_C",double(sum_words(q.conductor_C))},{"exterior_bank_C",double(sum_words(q.exterior_C))},
      {"max_gauss_defect_C",c.max_gauss_C},{"max_gauss_bound_C",c.max_gauss_bound_C},{"max_divB_T_m",c.max_divB_T_m},
      {"gauss_increment_C",r.gauss_increment_C},{"motion_continuity_C",b.continuity_C}}},
    {"energy",{{"field_delta_J",r.field_delta_J},{"plasma_work_J",r.plasma_work_J},{"joule_J",r.joule_J},{"outer_loss_J",r.outer_loss_J},
      {"outer_supply_J",r.outer_supply_J},{"circuit_work_J",r.circuit_work_J},{"kinetic_change_J",b.kinetic_change_J},
      {"thermal_delta_J",fh.heat().thermal_delta_J},{"field_balance_J",r.balance_J},{"combined_balance_J",p.combined_balance_J()}}},
    {"motion",{{"plasma_dt",b.plasma_dt},{"maximum_gyro",b.maximum_gyro},{"maximum_displacement_cells",b.maximum_displacement_cells}}},
    {"circuit",{{"current_A",r.current_A},{"terminal_V",r.terminal_V}}},
    {"field_solve",{{"backend",fh.solve_record().backend},{"iterations",fh.solve_record().iterations},
      {"residual_over_target",fh.solve_record().target>0?fh.solve_record().independent_residual/fh.solve_record().target:0.}}}};
}

// (b) the same step with the GPU push and the GPU field solve.
int run_gpu_step_mode(const Inputs& in,const fs::path& out,double water_weight){
  Results R;J receipt;
  if(backends.field!="gpu"||backends.push!="gpu")throw std::runtime_error("gpu-step requires --field gpu --push gpu");
  const fs::path scratch=out/"scratch";fs::create_directories(scratch);
  receipt["device"]=J::parse(backends.gpu()->device_json());
  receipt["configuration_sha256"]=in.model->configuration_sha256();
  receipt["field_solve_backend"]={{"name",in.model->solver().name()},{"descriptor_sha256",in.model->solver().descriptor_sha256()},
    {"compile",J::parse(hg::gpu_field_solver_compile_json(in.model->solver()))}};
  for(const auto& [label,packets]:std::vector<std::pair<std::string,std::vector<m::FullTimeAuthoredPacket>>>{
        {"P1",packets_p1(centre_plan)},{"P2",packets_p2(water_weight,centre_plan)}}){
    J rec;
    auto F=genesis(in,packets);
    auto parity=std::make_shared<ParityMotionProvider>(F.motion);F.motion=parity;
    rec["genesis"]=F.genesis;
    const auto t=steady::now();
    std::optional<v5::ClosedParentTrialV5> trial;
    try{trial.emplace(prepare(F));}catch(const std::exception& e){R.record("GPU",label+": prepare with GPU push and GPU field",false,{{"error",e.what()}});continue;}
    rec["prepare_s"]=since(t);rec["prepare_phases"]=times_json(trial->times());
    const auto& proposal=TA::proposal(*trial);
    rec["guard_receipt"]=guard_receipt(proposal);
    rec["motion_seal_sha256"]=proposal.motion().seal_sha256();
    rec["parent_events"]=count_events(trial->audit());
    const auto pub=publish_measured(F,*trial);
    J pd=pub;for(const char* k:{"total_s","can_publish_s","swap_window_s","release_s"})pd.erase(k);
    rec["publish"]=pd;rec["publish_timing"]=pub;
    rec["push_parity"]={{"iterates_compared",parity->compared},{"iterates_identical",parity->identical},{"per_iterate",parity->details}};
    const auto& g=rec["guard_receipt"];
    R.record("GPU",label+": GPU push bit-identical to the CPU provider on every fixed-point iterate (motion v3 record, current, charge arrays, endpoint charge words)",
      parity->compared>0&&parity->compared==parity->identical,rec["push_parity"]);
    R.record("GPU",label+": GPU field accepted through the guard receipt (independent residual, absolute Gauss, divB, Gauss increment, energy)",
      g.at("backend")=="gpu1-factored-cg"&&g.at("residual_over_target").get<double>()<=1&&
      g.at("max_gauss_C").get<double>()<=g.at("max_gauss_bound_C").get<double>()&&g.at("max_divB_T_m").get<double>()<=g.at("max_divB_bound_T_m").get<double>(),g);
    R.record("GPU",label+": joint step with GPU backends published (zero allocations between first and last swap)",
      pub.at("published").get<bool>()&&pub.value("allocations_between_first_and_last_swap",-1LL)==0,pd);
    if(pub.at("published").get<bool>()){
      rec["after"]=snapshot(F,scratch,label+"-gpu-after");
      rec["owner_bytes_digest"]=rec["after"]["digest"];
      R.record("GPU",label+": the backend is recorded in the configuration identity",
        F.owner->configuration_sha256()==in.model->configuration_sha256()&&F.genesis.at("field_solve_backend")=="gpu1-factored-cg");
    }
    receipt[label]=rec;
  }
  receipt["teardown_deferred"]="session closed at process exit";
  receipt["schema"]="stage1-v5-gpu-step-receipt-v1";
  write_text(out/"gpu-step-receipt.json",receipt.dump(1)+"\n");
  J deterministic=J::object();
  for(const char* l:{"P1","P2"})if(receipt.contains(l))deterministic[l]={{"guard_receipt",receipt[l]["guard_receipt"]},
    {"motion_seal_sha256",receipt[l]["motion_seal_sha256"]},{"owner_bytes_digest",receipt[l].value("owner_bytes_digest",J())}};
  write_text(out/"gpu-step-deterministic.json",deterministic.dump(1)+"\n");
  write_text(out/"gpu-step-result.json",J{{"cases",R.cases},{"summary",R.summary()}}.dump(1)+"\n");
  fs::remove_all(scratch);
  for(const auto& c:R.cases)if(!c.at("passed").get<bool>())return 1;
  return 0;
}

// Population of N carriers: N/2 co-located electron/H2O+ packet pairs (counts
// {1,0}) at seeded random positions in the vacuum gap, electrons 1..300 eV.
std::vector<m::FullTimeAuthoredPacket> pair_packets(std::size_t carriers,std::uint64_t seed){
  if(carriers==2)return packets_p1(centre_plan);
  std::mt19937_64 rng(seed);std::uniform_real_distribution<double> ux(.0035,.0055),uz(-.004,.004),ue(0,1);
  std::vector<m::FullTimeAuthoredPacket> out;out.reserve(carriers);
  for(std::size_t k=0;k<carriers/2;++k){
    const m::PicVector x{ux(rng),ux(rng),uz(rng)};
    const double eV=1+299*ue(rng),ct=2*ue(rng)-1,ph=6.283185307179586*ue(rng),st=std::sqrt(std::max(0.,1-ct*ct));
    const double xe=eV*pc::qe/(w::me*pc::c*pc::c),u=pc::c*std::sqrt(xe*(xe+2));
    out.push_back({"e"+std::to_string(k),"electron","multistep/authored-electron",{1,0},x,{u*st*std::cos(ph),u*st*std::sin(ph),u*ct}});
    out.push_back({"i"+std::to_string(k),"H2O+","multistep/authored-thermal-ion",{1,0},x,{372*(2*ue(rng)-1),372*(2*ue(rng)-1),372*(2*ue(rng)-1)}});
  }
  return out;
}
J phase_row(const v5::PhaseTimesV5& t,const J& pub,const hc::MotionHostTimers& host0,const hc::MotionHostTimers& host1){
  const auto& hl=t.hall;
  return {{"coordinator_observe_ms",1e3*t.observe_s},{"native_begin_ms",1e3*t.begin_s},{"chemistry_first_exposure_ms",1e3*t.first_exposure_s},
    {"native_census_ms",1e3*t.census_s},{"hall_admit_census_ms",1e3*hl.admit_census_s},{"hall_motion_ms",1e3*hl.motion_s},
    {"hall_motion_census_hash_ms",1e3*(host1.census_hash_s-host0.census_hash_s)},{"hall_motion_record_v3_ms",1e3*(host1.record_v3_s-host0.record_v3_s)},
    {"hall_motion_v2_json_seal_ms",1e3*hl.motion_v2_seal_s},
    {"hall_deposit_support_ms",1e3*hl.support_s},{"hall_field_prepare_ms",1e3*hl.field_prepare_s},{"hall_field_solve_ms",1e3*hl.solve_s},
    {"hall_independent_residual_ms",1e3*hl.residual_check_s},{"hall_field_seal_ms",1e3*hl.seal_s},{"hall_heat_ms",1e3*hl.heat_s},
    {"hall_proposal_total_ms",1e3*t.hall_proposal_s},{"coordinator_binding_checks_ms",1e3*t.binding_checks_s},
    {"native_motion_admission_ms",1e3*t.motion_admission_s},{"chemistry_second_exposure_ms",1e3*t.second_exposure_s},
    {"native_and_chemistry_seal_ms",1e3*t.seal_s},{"coordinator_audit_ms",1e3*t.audit_s},
    {"can_publish_ms",1e3*pub.value("can_publish_s",0.)},{"swap_window_ms",1e3*pub.value("swap_window_s",0.)},{"release_ms",1e3*pub.value("release_s",0.)},
    {"fixed_point_iterations",t.hall_iterations}};
}
J phase_summary(const std::vector<J>& rows){
  J s=J::object();if(rows.empty())return s;
  for(const auto& [k,v]:rows.front().items()){(void)v;std::vector<double> x;for(const auto& r:rows)x.push_back(r.at(k).get<double>());
    std::sort(x.begin(),x.end());double mean=0;for(double y:x)mean+=y;mean/=double(x.size());
    s[k]={{"median",x[x.size()/2]},{"mean",mean},{"p90",x[std::min(x.size()-1,std::size_t(.9*double(x.size())))]},{"max",x.back()}};}
  return s;
}
// (c) multi-step joint run until candidate003's history cap (or max_steps).
int run_multistep_mode(const Inputs& in,const fs::path& out,std::size_t carriers,long max_steps,long parity_every){
  J summary={{"schema","stage1-v5-multistep-summary-v1"},{"carriers_requested",carriers},{"field_backend",backends.field},{"push",backends.push},
    {"configuration_sha256",in.model->configuration_sha256()},{"max_steps",max_steps},{"parity_every",parity_every}};
  if(backends.device)summary["device"]=J::parse(backends.device->device_json());
  std::ofstream ledger(out/"step-ledger.jsonl");
  const auto t_gen=steady::now();
  std::optional<v5::JointOwnersV5> owners;
  try{owners.emplace(genesis(in,pair_packets(carriers,20261030),std::max<std::size_t>(carriers+64,1024)));}
  catch(const std::exception& e){
    summary["genesis_refusal"]=e.what();summary["steps_completed"]=0;summary["genesis_s"]=since(t_gen);
    write_text(out/"multistep-summary.json",summary.dump(1)+"\n");std::cout<<"genesis refused: "<<e.what()<<std::endl;return 0;}
  auto& F=*owners;F.inputs.audit_transport_record=false;
  summary["genesis_s"]=since(t_gen);summary["census_rows"]=F.arena->rate_population(F.inputs.maximum_carriers).size();
  summary["motion_provider"]=F.motion->name();
  const auto device_motion=F.motion;
  std::vector<J> rows;long steps=0;std::string stop="max_steps";long parity_compared=0,parity_identical=0;
  const auto wall=steady::now();
  for(long step=1;step<=max_steps;++step){
    std::shared_ptr<ParityMotionProvider> parity;
    if(backends.push=="gpu"&&parity_every>0&&(step==1||step%parity_every==0)){parity=std::make_shared<ParityMotionProvider>(device_motion);F.motion=parity;}
    else F.motion=device_motion;
    const auto host0=hc::motion_host_timers();
    const auto t0=steady::now();
    std::optional<v5::ClosedParentTrialV5> trial;
    try{trial.emplace(prepare(F));}catch(const std::exception& e){stop=std::string("prepare refused at step ")+std::to_string(step)+": "+e.what();break;}
    const double prepare_s=since(t0);
    const auto host1=hc::motion_host_timers();
    const auto& proposal=TA::proposal(*trial);
    J row={{"step",step},{"ledger",step_ledger(proposal)},{"prepare_ms",1e3*prepare_s},{"parent_events",count_events(trial->audit())}};
    const auto times=trial->times();
    const auto pub=publish_measured(F,*trial);
    if(!pub.at("published").get<bool>()){stop="publication refused at step "+std::to_string(step)+" (refusal point "+pub.at("refusal_point").dump()+")";break;}
    row["phases_ms"]=phase_row(times,pub,host0,host1);
    row["step_total_ms"]=1e3*(prepare_s+pub.at("total_s").get<double>());
    row["time_s"]=F.owner->time_s();row["field_sha256"]=F.owner->field_sha256();row["heat_sha256"]=F.owner->heat_sha256();
    row["native_prefix"]=F.arena->state_prefix();
#ifdef COORDINATOR_NATIVE_V5
    {const auto hs=F.arena->full_time_history_status();row["journal_rows"]=hs.cumulative_rows;row["journal_bytes"]=hs.cumulative_bytes;
      row["history_anchors"]=hs.anchors;row["history_window_rows"]=hs.window_rows;}
#else
    row["journal_bytes"]=F.arena->snapshot()->journal_bytes;row["journal_rows"]=F.arena->snapshot()->journal.size();
#endif
    if(parity){row["push_parity"]={{"compared",parity->compared},{"identical",parity->identical}};parity_compared+=parity->compared;parity_identical+=parity->identical;}
    ledger<<row.dump()<<'\n';ledger.flush();
    rows.push_back(row["phases_ms"]);++steps;
    if(step%50==0)std::cout<<"step "<<step<<" "<<row["step_total_ms"]<<" ms"<<std::endl;
  }
  summary["steps_completed"]=steps;summary["stop"]=stop;summary["wall_s"]=since(wall);
  summary["final_time_s"]=F.owner->time_s();summary["final_census_rows"]=F.arena->rate_population(F.inputs.maximum_carriers).size();
#ifdef COORDINATOR_NATIVE_V5
  {const auto hs=F.arena->full_time_history_status();summary["final_journal_bytes"]=hs.cumulative_bytes;summary["final_journal_rows"]=hs.cumulative_rows;
    summary["history"]={{"anchor_interval_steps",hs.anchor_interval_steps},{"anchors",hs.anchors},{"anchor_bytes",hs.anchor_bytes},
      {"window_rows",hs.window_rows},{"window_bytes",hs.window_bytes},{"history_policy_sha256",hs.history_policy_sha256}};
    summary["native"]="native v5 (candidate005)";summary["native_carrier_cap"]=F.arena->full_time_carrier_cap();}
#else
  summary["final_journal_bytes"]=F.arena->snapshot()->journal_bytes;summary["final_journal_rows"]=F.arena->snapshot()->journal.size();
  summary["native"]="native candidate003";
#endif
  summary["hall_carrier_cap"]=in.model->limits().maximum_carriers;
  summary["phase_ms"]=phase_summary(rows);
  {std::vector<double> tot;std::ifstream in2(out/"step-ledger.jsonl");std::string line;
    while(std::getline(in2,line)){auto j=J::parse(line);tot.push_back(j.at("step_total_ms").get<double>());}
    if(!tot.empty()){auto sorted=tot;std::sort(sorted.begin(),sorted.end());double mean=0;for(double x:tot)mean+=x;
      summary["step_total_ms"]={{"first",tot.front()},{"median",sorted[sorted.size()/2]},{"mean",mean/double(tot.size())},{"max",sorted.back()}};}}
  if(backends.push=="gpu")summary["push_parity"]={{"iterates_compared",parity_compared},{"iterates_identical",parity_identical}};
  if(backends.field=="gpu"){const auto st=hg::gpu_field_solver_stats(in.model->solver());summary["gpu_field_solves"]={{"solves",st.solves},{"total_s",st.total_s}};}
  if(backends.device)summary["device_memory"]=J::parse(backends.device->memory_json());
  write_text(out/"multistep-summary.json",summary.dump(1)+"\n");
  std::cout<<summary.dump(1)<<std::endl;
  return parity_compared==parity_identical?0:1;
}

// Hall-only closed loop (no native arena or chemistry): used where candidate003's
// journal byte cap blocks the joint step (1e4 after a few steps, 1e5 at genesis).
struct LoopBoundary {std::string field_sha256,heat_sha256,source_sha256,mesh_sha256,material_sha256,common_clock_sha256;double time_s=0;std::uint64_t accepted_step=0;};
int run_hall_loop_mode(const Inputs& in,const fs::path& out,std::size_t carriers,long steps,long parity_every){
  const auto& model=*in.model;const auto& g=model.mesh();const auto o=hall_geometry(g);const auto seed=model.dc_seed();
  std::vector<m::ArenaRateEntry> census;
  {std::uint64_t id=1;for(const auto& p:pair_packets(carriers,20261030)){m::ArenaRateEntry r;r.packet_id=p.id;r.state=p.state;r.native.id=id++;
    r.native.position_m=p.position_m;r.native.proper_velocity_m_s=p.proper_velocity_m_s;
    const bool e=p.state=="electron";r.native.charge_C=e?-m::pic_e:m::pic_e;r.native.mass_kg=e?w::me:states()[2].mass_kg;r.native.weight=1;census.push_back(r);}}
  J summary={{"schema","stage1-v5-hall-loop-summary-v1"},{"scope","Hall v2 closed loop only: no native arena, chemistry or common commit"},
    {"carriers",census.size()},{"field_backend",backends.field},{"push",backends.push},{"configuration_sha256",model.configuration_sha256()}};
  auto ledger=v5::GenesisTestAccess::declared_ledger(model,seed);
  ledger.particles=hc::charge_of_census(g,o,census,model.limits().maximum_carriers);
  hc::Owner owner(in.model,seed,ledger,o,census);
  std::shared_ptr<hc::MotionProvider> device;
  if(backends.push=="gpu")device=std::make_shared<hg::GpuMotionProvider>(backends.gpu(),owner.motion_context(),census.size());
  else device=std::make_shared<hc::CpuMotionProvider>();
  std::ofstream ledger_file(out/"hall-loop-ledger.jsonl");std::vector<J> rows;long done=0,pc=0,pi=0;std::string stop="steps";
  const std::string source=nw::sha("hall-loop/closed-source");
  for(long step=1;step<=steps;++step){
    std::shared_ptr<hc::MotionProvider> provider=device;std::shared_ptr<ParityMotionProvider> parity;
    if(backends.push=="gpu"&&parity_every>0&&(step==1||step%parity_every==0)){parity=std::make_shared<ParityMotionProvider>(device);provider=parity;}
    h::MotionBinding b;b.phase=h::ParticlePhase::full_time_position_and_proper_velocity;b.accepted_step=owner.revision();b.staged_arena_revision=owner.revision()+1;
    b.begin_s=owner.time_s();b.dt_s=dt;b.end_s=b.begin_s+dt;b.staged_arena_sha256=nw::sha("hall-loop/staged/"+std::to_string(step));
    b.mesh_sha256=owner.motion_mesh_sha256();b.material_sha256=model.material_sha256();b.species_sha256=nw::sha("hall-loop/species");
    b.common_clock_sha256=nw::sha("hall-loop/clock/"+std::to_string(step));b.kinetic_policy_sha256=owner.motion_context().policy_sha256;
    const auto host0=hc::motion_host_timers();const auto t0=steady::now();
    std::optional<hc::JointProposal> p;
    try{p.emplace(hc::prepare_closed_proposal(owner,o,census,b,source,*provider));}
    catch(const std::exception& e){stop=std::string("refused at step ")+std::to_string(step)+": "+e.what();break;}
    const double prep=since(t0);const auto host1=hc::motion_host_timers();
    LoopBoundary bd{p->candidate_field_sha256(),p->candidate_heat_sha256(),source,b.mesh_sha256,b.material_sha256,b.common_clock_sha256,p->end_s(),p->endpoint_revision()};
    J row={{"step",step},{"ledger",step_ledger(*p)}};
    auto t=steady::now();const bool ok=owner.preflight_final(*p,census,bd);const double pre=since(t);
    if(!ok){stop="preflight refused at step "+std::to_string(step);break;}
    std::vector<m::GapParticle> ends(p->motion().endpoints().begin(),p->motion().endpoints().end());
    const auto& hl=p->timings();
    t=steady::now();owner.publish_prevalidated(std::move(*p));const double publ=since(t);
    {std::vector<m::ArenaRateEntry*> sorted;for(auto& r:census)sorted.push_back(&r);
      std::sort(sorted.begin(),sorted.end(),[](auto a,auto c){return a->native.id<c->native.id;});
      for(std::size_t i=0;i<sorted.size();++i)sorted[i]->native=ends[i];}
    J ph={{"hall_admit_census_ms",1e3*hl.admit_census_s},{"hall_motion_ms",1e3*hl.motion_s},
      {"hall_motion_census_hash_ms",1e3*(host1.census_hash_s-host0.census_hash_s)},{"hall_motion_record_v3_ms",1e3*(host1.record_v3_s-host0.record_v3_s)},
      {"hall_motion_v2_json_seal_ms",1e3*hl.motion_v2_seal_s},
      {"hall_deposit_support_ms",1e3*hl.support_s},{"hall_field_prepare_ms",1e3*hl.field_prepare_s},{"hall_field_solve_ms",1e3*hl.solve_s},
      {"hall_independent_residual_ms",1e3*hl.residual_check_s},{"hall_field_seal_ms",1e3*hl.seal_s},{"hall_heat_ms",1e3*hl.heat_s},
      {"hall_proposal_total_ms",1e3*prep},{"hall_preflight_final_ms",1e3*pre},{"hall_publish_ms",1e3*publ},{"fixed_point_iterations",hl.iterations}};
    row["phases_ms"]=ph;row["step_total_ms"]=1e3*(prep+pre+publ);row["field_sha256"]=owner.field_sha256();
    if(parity){row["push_parity"]={{"compared",parity->compared},{"identical",parity->identical}};pc+=parity->compared;pi+=parity->identical;}
    ledger_file<<row.dump()<<'\n';ledger_file.flush();rows.push_back(ph);++done;
  }
  summary["steps_completed"]=done;summary["stop"]=stop;summary["phase_ms"]=phase_summary(rows);
  if(backends.push=="gpu")summary["push_parity"]={{"iterates_compared",pc},{"iterates_identical",pi}};
  if(backends.device)summary["device_memory"]=J::parse(backends.device->memory_json());
  write_text(out/"hall-loop-summary.json",summary.dump(1)+"\n");
  std::cout<<summary.dump(1)<<std::endl;
  return pc==pi?0:1;
}
} // namespace stage1

int main(int argc,char** argv){
  try{
    if(argc<7)throw std::runtime_error("usage: MODE TRACE H2O_DATA BINDING OUTPUT_DIRECTORY NATIVE_SOURCE_IDENTITY [ARGS] "
      "[--field cpu|gpu] [--push cpu|gpu] [--uuid UUID] [--nvrtc PATH]");
    std::vector<std::string> args;
    for(int i=7;i<argc;++i){const std::string a=argv[i];
      auto value=[&]{if(i+1>=argc)throw std::runtime_error("missing option value");return std::string(argv[++i]);};
      if(a=="--field")backends.field=value();else if(a=="--push")backends.push=value();
      else if(a=="--uuid")backends.uuid=value();else if(a=="--nvrtc")backends.nvrtc=value();
      else if(a=="--carrier-cap")caps.carriers=std::size_t(std::stod(value()));
      else if(a=="--hall-budget")caps.hall_budget=std::uint64_t(std::stod(value()));
      else if(a=="--native-budget")caps.native_budget=std::uint64_t(std::stod(value()));else args.push_back(a);}
    const std::string mode=argv[1];const fs::path out=argv[5];
    if(std::string(m::PicArena::source_identity())!=argv[6])throw std::runtime_error("linked native source identity differs");
    fs::create_directories(out);
    const auto in=load_inputs(argv[2],argv[3],argv[4]);
    auto arg=[&](std::size_t k,double d){return args.size()>k?std::stod(args[k]):d;};
    int rc=2;
    if(mode=="t1")rc=run_t1_mode(in,out,arg(0,1e18));
    else if(mode=="matrix")rc=run_matrix_mode(in,out,arg(0,900.));
    else if(mode=="gpu-step")rc=run_gpu_step_mode(in,out,arg(0,1e18));
    else if(mode=="multistep")rc=run_multistep_mode(in,out,std::size_t(arg(0,2)),long(arg(1,1000)),long(arg(2,0)));
    else if(mode=="sweep-point")rc=run_sweep_point_mode(in,out,(long long)arg(0,0),(long long)arg(1,0));
    else if(mode=="hall-loop")rc=run_hall_loop_mode(in,out,std::size_t(arg(0,2)),long(arg(1,20)),long(arg(2,0)));
    else throw std::runtime_error("unknown mode");
    if(backends.device){const auto td=backends.device->close_json();std::ofstream(out/"gpu-teardown.json")<<td<<'\n';}
    return rc;
  }catch(const std::exception& e){std::cerr<<"fatal: "<<e.what()<<'\n';return 2;}
}
