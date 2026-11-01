#include "chem_v8_fixture.hpp"
#include "raikou/mfem_pic_full_time_device.h"
#include "native_gas_targets.hpp"
#include "target_step_fixture.hpp"
#include <fstream>
#include <filesystem>
#include <iostream>
namespace n=raikou::mfem_backend;namespace g=pf_gas_v1;namespace a=native_gas_stage;namespace t=native_gas_targets;
using J=nlohmann::json;int checks=0;
void check(bool b,const char *s){++checks;if(!b)throw std::runtime_error(s);}
template<class F>void refuses(F f,const char*s){bool bad=false;try{f();}catch(const std::exception&){bad=true;}check(bad,s);}
template<class F>void typed(F f,n::FullTimeChemistryRefusal::Kind k,const char*s){bool bad=false;try{f();}catch(const n::FullTimeChemistryRefusal&e){bad=e.kind()==k;}check(bad,s);}
g::Configuration gas_config(const n::GapPicOptions&o,const n::FullTimeChemistryRegistry&r,const pfc::Blob&b){
 g::Configuration c;c.preset="diagnostic-host-target-feedback";c.feedback=true;auto grid=n::chem8::make_grid(o,{std::uint32_t(o.nx),std::uint32_t(o.ny),5});
 c.grid={grid.low,grid.high,grid.n,grid.volume};c.chemistry_blob_sha256=b.sha256;c.chemistry_registry_sha256=r.registry_sha256;
 for(const auto&s:r.species)c.species.push_back({s.state,s.mass_kg,s.charge_e,s.target});
 for(std::uint32_t k=0;k<b.n_proc;++k){const auto&p=b.view.procs[k];c.processes.push_back({b.view.groups[p.group].projectile,p.n_out,a::internal_carrier(p.kind,p.n_out),p.kind});}return c;
}
n::FullTimeStepBinding binding(const chemfix::Arena&A,double dt){const auto b=*A.arena->full_time_boundary();return {b.accepted_step,b.accepted_arena_revision,b.time_s,b.time_s+dt,dt,b.accepted_arena_sha256,b.field_sha256,b.mesh_sha256,b.material_sha256,chemfix::sha("target-feedback-control-clock")};}
struct Scenario {
 n::GapPicOptions o; n::FullTimeChemistryConfiguration c;std::shared_ptr<g::Owner> gas;chemfix::Arena A;
 Scenario(const pfc::Blob&blob,n::FullTimeChemistryConfiguration config,const std::vector<n::FullTimeAuthoredPacket>&packets)
 :o(chemfix::options()),c(std::move(config)) {
  auto reg=n::full_time_chemistry_registry(c,o);gas=std::make_shared<g::Owner>(gas_config(o,reg,blob));
  c.host_target_feedback=true;c.host_target_contract_sha256=t::contract(c,o,*gas);
  A=chemfix::make_arena(o,c,packets,120000,2);
 }
 std::shared_ptr<t::Provider> provider(const pfc::Blob&blob){return std::make_shared<t::Provider>(gas,A.reg,blob,c,o);}
};
std::shared_ptr<const n::FullTimeChemistryStageReceipt> first(chemfix::Arena&A,std::shared_ptr<n::FullTimeHostTargetProvider> p,
 const std::function<void(n::CollisionScope&,n::PreparedArena&)>&inspect={}){
 const auto before=A.arena->full_time_boundary()->accepted_arena_sha256;const auto rng=A.arena->full_time_chemistry().state_sha256;
 std::shared_ptr<const n::FullTimeChemistryStageReceipt> result;
 {auto b=binding(A,1e-12);auto outer=A.arena->begin_full_time_step(b,120000,64);
  n::ExposureBinding e{b,n::FullTimeSlot::before_transport,b.begin_s,.5*b.dt_s,A.arena->full_time_chemistry().stages,b.accepted_field_sha256};
  auto scope=outer.first_exposure(e);auto token=p?scope.prepare_chemistry_v3_host_targets(std::move(p)):scope.prepare_chemistry_v3();
  result=token.candidate_chemistry_stage();if(inspect)inspect(scope,token);
 }
 check(before==A.arena->full_time_boundary()->accepted_arena_sha256 && rng==A.arena->full_time_chemistry().state_sha256,"abandoned stage mutated accepted state/RNG");return result;
}
// Deliberately untrusted provider demonstrates native validation, independently
// of the genuine adapter's own stricter result check.
struct Fault: n::FullTimeHostTargetProvider {
 std::shared_ptr<t::Provider> delegate;int mode;int calls=0;std::function<void()> hook,validate_hook;
 Fault(std::shared_ptr<t::Provider>p,int m):delegate(std::move(p)),mode(m){}
 n::FullTimeHostTargetDraft apply(std::shared_ptr<const n::FullTimeHostTargetInputs>x) override {
  ++calls;if(hook)hook();auto y=delegate->apply(x);
  if(mode==1)y.cells.at(0).n[0]*=2;
  if(mode==2)y.cells.at(0).pad=99;
  if(mode==3)y.before_targets_root=chemfix::sha("wrong-root");
  if(mode==4)y.cells.at(0).s[0]=-1;
  if(mode==5)y.cells.at(0).drift[0][0]=std::numeric_limits<double>::infinity();
  if(mode!=3)y.after_targets_root=n::full_time_host_target_cells_sha256(y.cells);
  return y;
 }
 void validate_current(std::shared_ptr<const n::FullTimeHostTargetInputs>,const n::FullTimeHostTargetDraft&)const override{if(validate_hook)validate_hook();}
};
int main(int argc,char**argv){try{
 if(argc!=2)throw std::runtime_error("output directory required");
 std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
 const auto blob=pfc::load_blob(chemfix::blob_path(),chemfix::blob_sha());const auto o=chemfix::options();auto base=chemfix::chemistry_configuration(o);base.retain_stage_inputs=true;
 chemfix::Population pop;pop.carriers=3000;pop.seed=11;pop.clustered=true;pop.electron_fraction=.49;pop.ion_fraction=.49;pop.water_weight=2e13;pop.charged_weight=1e-10;pop.electron_eV_min=5;pop.electron_eV_max=300;
 const auto packets=chemfix::population(pop);Scenario S(blob,base,packets);auto reference=chemfix::make_arena(o,base,packets,120000,2);
 const auto expected=first(reference,{});S.gas->begin_step(0);auto provider=S.provider(blob);
 auto actual=first(S.A,provider,[&](auto&scope,auto&token){
  check(scope.can_publish(token),"actual live native token");token.validate_candidate_host_targets();
  const auto r=token.candidate_chemistry_stage();const auto f=r->host_target_feedback;
  check(f && f->input->native_rows==pop.carriers && f->input->census_sha256.size()==64,"authentic all-species census retained");
  check(f->output.before_targets_root==f->output.after_targets_root && provider->feedback().changed==0,"empty gas correction is exact no-op");
  check(r->admitted_export->host_target_feedback_sha256==f->sha256,"admitted export binds feedback");
  auto copy=std::make_shared<n::FullTimeHostTargetInputs>(*f->input);
  refuses([&]{provider->validate_current(copy,f->output);},"equal-content foreign native input refused");
  refuses([&]{(void)provider->apply(f->input);},"one-shot provider reused");
  auto draft=f->output;draft.input_sha256=chemfix::sha("tampered");refuses([&]{provider->validate_current(f->input,draft);},"changed provider result refused");
  const auto digest=S.gas->candidate_sha256();S.gas->prepare_keys({},{},{});
  check(S.gas->candidate_sha256()==digest,"same-digest mutation control");
  refuses([&]{token.validate_candidate_host_targets();},"same-digest gas mutation invalidates feedback capability");
 });
 check(actual->events>0 && actual->admitted_events>0,"populated eventful diagnostic fixture");
 check(actual->records_root==expected->records_root && actual->records_sha256==expected->records_sha256 && actual->ledger_sha256==expected->ledger_sha256 && actual->cell_bytes==expected->cell_bytes && actual->carrier_bytes==expected->carrier_bytes,"zero-correction physical records/ledgers/inputs identical");
 check(actual->inputs_sha256!=expected->inputs_sha256 && actual->state_after_sha256!=expected->state_after_sha256,"configuration/feedback identities intentionally distinct");
 check(!actual->full_system_energy_closure,"diagnostic stage does not claim closure");S.gas->abandon();
 for(int mode=1;mode<=5;++mode){S.gas->begin_step(0);auto f=std::make_shared<Fault>(S.provider(blob),mode);typed([&]{first(S.A,f);},n::FullTimeChemistryRefusal::Kind::target_fields,"native target-byte fault not refused");S.gas->abandon();}
 typed([&]{first(S.A,{});},n::FullTimeChemistryRefusal::Kind::target_fields,"missing provider silently accepted");
 for(auto mode:{n::FullTimeExecution::device,n::FullTimeExecution::device_verified})typed([&]{S.A.arena->attach_full_time_device({},mode);},n::FullTimeChemistryRefusal::Kind::target_fields,"device target feedback silently fell back");
 S.gas->begin_step(1);refuses([&]{first(S.A,S.provider(blob));},"gas/native step mismatch accepted");S.gas->abandon();
 S.gas->begin_step(0);first(S.A,S.provider(blob),[&](auto&,auto&token){auto q=token.candidate_chemistry_stage();S.gas->abandon();S.gas->begin_step(0);refuses([&]{token.validate_candidate_host_targets();},"same-state retry accepted");check(q->host_target_feedback!=nullptr,"immutable facts survive retry");});S.gas->abandon();
 auto strict=base;strict.diagnostic_field_target_approximation=false;Scenario T(blob,strict,packets);T.gas->begin_step(0);
 typed([&]{first(T.A,T.provider(blob));},n::FullTimeChemistryRefusal::Kind::target_energy_owner,"strict eventful target-owner gate opened");T.gas->abandon();
 // Normal public construction order proves scope -> gas registry -> target
 // feedback contract is acyclic. No hand-authored configuration hashes.
 auto scoped=base;scoped.water_scope_contract=n::full_time_minimum_water_scope_contract(scoped,o);scoped.water_scope_contract_sha256=chemfix::sha(scoped.water_scope_contract);scoped.water_scope_mode=n::FullTimeWaterScopeMode::minimum_water;
 std::vector<n::FullTimeAuthoredPacket> neutrals{{"cold/e","pfc-v3/e","manufactured/cold",n::ArenaCount(1),{.00215,.00215,-.01},{0,0,0}},
  {"cold/water","pfc-v3/H2O","manufactured/cold",n::ArenaCount(1e6),{.00215,.00215,-.01},{0,0,0}}};
 Scenario C(blob,scoped,neutrals);check(C.A.reg.configuration_sha256.size()==64 && C.c.host_target_contract_sha256.size()==64,"normal scope/target construction");
 // Manufactured positive moment correction through the chemistry-stage API.
 // There is no corresponding physical native debit: this isolates target
 // injection and is deliberately abandoned, not a closure qualification.
 C.gas->begin_step(0);std::shared_ptr<const n::FullTimeHostTargetReceipt> stage1_targets;
 auto correction_report=chemfix::step_with_targets(C.A,1e-15,0,[&](n::CollisionScope&scope){
  auto tok=scope.prepare_chemistry_v3_host_targets(C.provider(blob));tok.validate_candidate_host_targets();
  const auto r=tok.candidate_chemistry_stage();const auto f=r->host_target_feedback;
  check(r->events==0,"controlled correction fixture unexpectedly reactive");
  if(r->stage==1){stage1_targets=f;std::vector<g::FlowRecord> heat;
   for(std::size_t i=0;i<f->input->cells.size();++i)for(std::uint32_t k=0;k<f->input->cells[i].n_targets;++k){const auto sp=f->input->cells[i].species[k];if(C.A.reg.species[sp].charge_e!=0)continue;
    for(double sign:{-1.,1.}){g::FlowRecord q{};q.kind=pfg::kFlowDebit;q.negative=1;q.grid_cell=f->input->atlas[i];q.species=sp;q.native_id=1+heat.size();q.weight=1e5;q.u[0]=sign*1e4;heat.push_back(q);}}
   check(!heat.empty(),"manufactured correction has neutral target");C.gas->take(0,1,{},heat);
  }else{
   check(f->output.before_targets_root!=f->output.after_targets_root,"nonzero correction reaches actual MCC target bytes");
   std::vector<g::CellTargetsAbi> independent(f->input->cells.size());std::memcpy(independent.data(),f->input->cells.data(),independent.size()*sizeof(g::CellTargetsAbi));
   const auto ir=C.gas->feedback(0,2,f->input->atlas,independent);
   check(std::memcmp(independent.data(),f->output.cells.data(),independent.size()*sizeof(g::CellTargetsAbi))==0 && t::feedback_digest(ir)==f->output.owner_feedback_sha256,"native targets equal independent gas feedback");
   check(f->input->sha256!=stage1_targets->input->sha256,"corrected stage2 authentic input differs");C.gas->take(0,2,{},{});
  }return tok;
 },{},{},false);
 check(C.A.arena->full_time_boundary()->accepted_step==0,"manufactured correction outer step was published");C.gas->abandon();
 auto run=[&](Scenario&X,int k){X.gas->begin_step(X.A.arena->full_time_boundary()->accepted_step);
  auto report=chemfix::step_with_targets(X.A,1e-15,k,[&](n::CollisionScope&s){auto p=X.provider(blob);auto tok=s.prepare_chemistry_v3_host_targets(p);tok.validate_candidate_host_targets();auto r=tok.candidate_chemistry_stage();check(r->events==0,"cold continuation fixture unexpectedly reactive");X.gas->take(r->accepted_step,r->stage,{},{});return tok;});
  X.gas->commit(X.gas->candidate_sha256());return report;};
 for(int k=0;k<3;++k){auto rr=run(C,k);check(rr.stage1->host_target_feedback->input->sha256!=rr.stage2->host_target_feedback->input->sha256,"stage2 reused stage1 input");
  const auto st=C.A.arena->full_time_chemistry();auto boundary=*C.A.arena->full_time_boundary();auto path=out/("checkpoint-"+std::to_string(k)+".json");C.A.arena->save(path.string());
  C.A.arena=n::PicArena::restore_full_time(path.string(),o,C.A.states,C.A.channels,boundary);auto restored=C.A.arena->full_time_chemistry();check(st.host_target_feedback_sha256.size()==64 && st.host_target_feedback_sha256==restored.host_target_feedback_sha256 && st.host_target_history_sha256==restored.host_target_history_sha256,"feedback roots lost at restart");
  J doc;std::ifstream(path)>>doc;auto payload=doc.at("payload");payload["configuration"]["v8"]["chemistry"]["host_target_contract_sha256"]=chemfix::sha("wrong-target-contract");auto wrong=out/("wrong-"+std::to_string(k)+".json");std::ofstream(wrong)<<J{{"payload",payload},{"sha256",chemfix::sha(payload.dump())}}.dump();refuses([&]{(void)n::PicArena::restore_full_time(wrong.string(),o,C.A.states,C.A.channels,boundary);},"re-signed feedback configuration tamper accepted");
  auto changed=doc.at("payload");auto hex=changed["metadata"]["chemistry"]["state_hex"].get<std::string>();std::string needle;const char*digits="0123456789abcdef";for(unsigned char v:st.host_target_feedback_sha256){needle+=digits[v>>4];needle+=digits[v&15];}auto pos=hex.find(needle);check(pos!=std::string::npos,"feedback root absent in codec");hex.replace(pos,2,hex.substr(pos,2)=="30"?"31":"30");changed["metadata"]["chemistry"]["state_hex"]=hex;auto badroot=out/("badroot-"+std::to_string(k)+".json");std::ofstream(badroot)<<J{{"payload",changed},{"sha256",chemfix::sha(changed.dump())}}.dump();refuses([&]{(void)n::PicArena::restore_full_time(badroot.string(),o,C.A.states,C.A.channels,boundary);},"re-signed feedback state root tamper accepted");
 }
 // Arbitrary callback lifetime/reentry failures must not dereference a dead
 // native owner. These controls do not imply the provider owns publication.
 for(bool late:{false,true}){
  Scenario dead(blob,base,packets);dead.gas->begin_step(0);auto fault=std::make_shared<Fault>(dead.provider(blob),0);
  auto destroy=[&]{dead.A.arena.reset();};if(late)fault->validate_hook=destroy;else fault->hook=destroy;
  refuses([&]{first(dead.A,fault);},"destroyed native owner survived callback validation");
  check(!dead.A.arena && dead.gas->step_open(),"native callback lifetime control");dead.gas->abandon();
 }
 S.gas->begin_step(0);
 {auto b=binding(S.A,1e-12);auto outer=S.A.arena->begin_full_time_step(b,120000,64);
  n::ExposureBinding e{b,n::FullTimeSlot::before_transport,b.begin_s,.5*b.dt_s,S.A.arena->full_time_chemistry().stages,b.accepted_field_sha256};
  auto scope=outer.first_exposure(e);auto fault=std::make_shared<Fault>(S.provider(blob),0);
  fault->hook=[&]{refuses([&]{(void)scope.prepare_chemistry_v3_host_targets(S.provider(blob));},"nested callback preparation accepted");};
  auto token=scope.prepare_chemistry_v3_host_targets(fault);check(scope.can_publish(token),"reentry refusal corrupted outer attempt");
 }S.gas->abandon();
 // Required e+OH encounter refuses before any provider call.
 std::vector<n::FullTimeAuthoredPacket> missing{{"electron","pfc-v3/e","manufactured/scope",n::ArenaCount(1),{.00215,.00215,-.01},{0,0,0}},{"radical","pfc-v3/OH","manufactured/scope",n::ArenaCount(1),{.00215,.00215,-.01},{0,0,0}}};
 Scenario M(blob,scoped,missing);M.gas->begin_step(0);auto counter=std::make_shared<Fault>(M.provider(blob),0);typed([&]{first(M.A,counter);},n::FullTimeChemistryRefusal::Kind::network_scope,"scope gate bypassed");check(counter->calls==0,"feedback ran before scope admission");M.gas->abandon();
 J result={{"passed",true},{"checks",checks},{"eventful_events",actual->events},{"admitted_events",actual->admitted_events},{"zero_correction_raw_records_equal",true},{"manufactured_nonzero_feedback",true},{"production_qualified",false},{"joint_energy_closure",false},{"device_feedback_qualified",false}};std::ofstream(out/"result.json")<<result.dump(2)<<'\n';std::cout<<result.dump(2)<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
