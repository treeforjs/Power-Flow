#include "coupled_step_fixture.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace n=raikou::mfem_backend;namespace g=pf_gas_v1;namespace a=native_gas_stage;namespace t=native_gas_targets;namespace c=native_gas_common;
using J=nlohmann::json;int checks=0;
void check(bool b,const char*s){++checks;if(!b)throw std::runtime_error(s);}
g::Configuration gas_config(const n::GapPicOptions&o,const n::FullTimeChemistryRegistry&r,const pfc::Blob&b){
 g::Configuration x;x.preset="bounded-strict-native-gas-reference";x.feedback=true;x.transport=g::Transport::census_conservative;
 auto grid=n::chem8::make_grid(o,{std::uint32_t(o.nx),std::uint32_t(o.ny),5});x.grid={grid.low,grid.high,grid.n,grid.volume};
 x.chemistry_blob_sha256=b.sha256;x.chemistry_registry_sha256=r.registry_sha256;
 for(const auto&s:r.species)x.species.push_back({s.state,s.mass_kg,s.charge_e,s.target});
 for(std::uint32_t k=0;k<b.n_proc;++k){const auto&p=b.view.procs[k];x.processes.push_back({b.view.groups[p.group].projectile,p.n_out,a::internal_carrier(p.kind,p.n_out),p.kind});}
 return x;
}
struct Scenario {
 n::GapPicOptions o;n::FullTimeChemistryConfiguration cfg;std::shared_ptr<g::Owner>gas;chemfix::Arena A;
 Scenario(const pfc::Blob&blob,const std::vector<n::FullTimeAuthoredPacket>&packets,bool minimum=false)
 :o(chemfix::options()),cfg(chemfix::chemistry_configuration(o)),
 gas([&]{cfg.diagnostic_field_target_approximation=false;if(minimum){cfg.water_scope_contract=n::full_time_minimum_water_scope_contract(cfg,o);cfg.water_scope_contract_sha256=native_gas_targets::sha(cfg.water_scope_contract);cfg.water_scope_mode=n::FullTimeWaterScopeMode::minimum_water;}auto r=n::full_time_chemistry_registry(cfg,o);return std::make_shared<g::Owner>(gas_config(o,r,blob));}()),
 A([&]{cfg.coupled_target_contract_sha256=c::coupled_contract(cfg,o,*gas);cfg.host_target_feedback=true;cfg.host_target_contract_sha256=t::contract(cfg,o,*gas);return chemfix::make_arena(o,cfg,packets,120000,2);}()){}
};
int main(int argc,char**argv){try{
 if(argc!=3)throw std::runtime_error("output path and control required");
 const std::string control=argv[2];std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
 check(control=="motion-empty-flows" || control.starts_with("allocation-motion-") || control=="snapshot-controls" || control=="correction-exits" || control=="correction-exterior" || control=="empty-stream-caps" || control=="motion-flows" || control=="same-digest-motion" || control=="normal" || control=="empty" || control=="abandon" || control=="same-digest-stage1" || control=="same-digest-stage2" || control=="abandon-final" || control=="regular-refusal" || control=="move-final" || control=="cell-crossing" || control=="two-step" || control=="destroyed-native" || control=="allocation-count" || control.starts_with("allocation-seal-"),"unknown control");
 auto blob=pfc::load_blob(chemfix::blob_path(),chemfix::blob_sha());chemfix::Population p;
 p.carriers=3000;p.seed=11;p.clustered=true;p.water_weight=2e13;p.charged_weight=1e-10;p.electron_eV_min=5;p.electron_eV_max=300;
 const bool empty_initial=control=="empty" || control=="empty-stream-caps" || control=="motion-empty-flows" || control.starts_with("allocation-motion-");
 auto packets=empty_initial?std::vector<n::FullTimeAuthoredPacket>{}:chemfix::population(p);
 // Place the bounded reference fixture away from cell edges. The unchanged
 // original placement remains a genuine crossing control.
 if(control!="cell-crossing")for(auto &packet:packets)for(int axis=0;axis<2;++axis)
   packet.position_m[axis]=.00205+.5*(packet.position_m[axis]-.0021);
 if(control=="correction-exits" || control=="correction-exterior")for(auto&packet:packets) {
   const bool water=packet.state=="pfc-v3/H2O";
   if(control=="correction-exits"){packet.position_m={.00405,water?.002+1e-10:.00205,0};if(water)packet.proper_velocity_m_s={0,-650,0};}
   else {packet.position_m={water?.009-1e-10:.00895,.00405,0};if(water)packet.proper_velocity_m_s={650,0,0};}
 }
 Scenario S(blob,packets,empty_initial);const auto native0=S.A.arena->state_prefix(),rng0=S.A.arena->full_time_chemistry().state_sha256,gas0=S.gas->accepted_sha256();
 std::shared_ptr<const c::CoupledTrial::Continuation> continuation;
 J record;bool retry_accepted=false;std::string retry_receipt_sha256;long seal_allocations=-1;bool refused=false,regular_refused=false;std::string error;chemfix::StepReport R;
 try{
  R=chemfix::coupled_step(S.A,1e-12,0,S.gas,
   [&](c::CoupledTrial&,n::CollisionScope&scope,unsigned stage){
    if(control=="regular-refusal" && stage==1){try{(void)scope.prepare_chemistry_v3_host_targets(std::make_shared<t::Provider>(S.gas,S.A.reg,blob,S.cfg,S.o));}
      catch(const n::FullTimeChemistryRefusal&e){regular_refused=e.kind()==n::FullTimeChemistryRefusal::Kind::target_energy_owner;}
      check(regular_refused,"regular public chemistry must not open coupled owner gate");}
    if(control=="same-digest-stage1" && stage==2){const auto prior=S.gas->candidate_sha256();S.gas->prepare_keys({},{},{});check(prior==S.gas->candidate_sha256(),"fixture mutation must preserve content digest");}
   },
   [&](c::CoupledTrial&){if(control.starts_with("allocation-seal-"))allocation_probe::fail_after=std::stol(control.substr(16));if(control=="allocation-count")allocation_probe::fail_after=std::numeric_limits<long>::max();if(control=="same-digest-stage2"){const auto prior=S.gas->candidate_sha256();S.gas->prepare_keys({},{},{});check(prior==S.gas->candidate_sha256(),"stage2 same-digest mutation");}},
   [&](c::CoupledTrial&trial){
    if(control=="allocation-count"){seal_allocations=std::numeric_limits<long>::max()-allocation_probe::fail_after;allocation_probe::fail_after=-1;}
    if(control=="destroyed-native"){S.A.arena.reset();check(!trial.can_publish() && !trial.try_publish(),"destroyed native owner must invalidate final token");}
    if(control=="abandon-final"){S.gas->abandon();check(!trial.can_publish() && !trial.try_publish(),"abandoned gas invalidates common publication");}
    if(control=="move-final"){c::CoupledTrial moved(std::move(trial));check(!trial.can_publish() && moved.can_publish(),"moved common capability");trial=std::move(moved);check(!moved.can_publish() && trial.can_publish(),"moved-back common capability");}
   },record,continuation,control!="abandon",control);
 }catch(const std::exception&e){allocation_probe::fail_after=-1;refused=true;error=e.what();}
 const bool expected_refusal=(control.starts_with("allocation-motion-") && control!="allocation-motion-count") || control=="same-digest-motion" || control=="same-digest-stage1" || control=="same-digest-stage2" || control=="abandon-final" || control=="destroyed-native" || control.starts_with("allocation-seal-");
 if(refused!=expected_refusal)throw std::runtime_error("unexpected common trial refusal outcome: "+error);
 const bool accepted=!refused && control!="abandon";
 if(accepted){check(S.A.arena->state_prefix()!=native0 && S.A.arena->full_time_boundary()->accepted_step==1 && !S.gas->step_open(),"both accepted owners advanced");check(record.at("stages").size()==2,"complete two-stage transcript");
  for(const auto&r:record.at("stages")){
   check(r.at("energy_exact_zero") && r.at("momentum_exact_zero") && r.at("quarantine_verified"),"verified quarantined actual inventory proof");
   for(const auto&v:r.at("momentum_residual_kg_m_s"))check(v.is_object() && v.at("digits")==g::digits_text(g::Exact{}.normalised().data()),"serialized exact momentum proof");
   for(const auto&v:r.at("raw_stored_momentum_balance"))check(v.is_object() && v.contains("digits"),"serialized raw momentum balance");
  }
  check(record.at("motion").at("correction_transport_exact_zero"),"motion correction proof");
  if(control=="motion-flows")check(record.at("motion").at("births")==4 && record.at("motion").at("sinks")==2 && record.at("motion").at("actual_sink_destinations").size()==2,"actual birth/conductor/exterior fixture");
  if(control=="correction-exits" || control=="correction-exterior") {
    const auto key=control=="correction-exits"?"conductor_correction":"exterior_correction";bool nonzero=false;
    for(const auto&v:record.at("motion").at(key))nonzero |= v.at("digits")!=g::digits_text(g::Exact{}.normalised().data());
    check(nonzero && record.at("motion").at("sinks").get<unsigned>()>0,"actual nonzero held correction exit required");
  }
  if(empty_initial)check(R.stage1->admitted_events==0 && R.stage2->admitted_events==0,"empty stage transcript must remain present");
  else check(R.stage1->admitted_events>0 && R.stage2->released_packets>0,"actual populated stages and release required");
 }else check((control=="destroyed-native" || (S.A.arena->state_prefix()==native0 && S.A.arena->full_time_chemistry().state_sha256==rng0)) && S.gas->accepted_sha256()==gas0 && !S.gas->step_open(),"refusal/abandon changed accepted native RNG or gas");
 if(control.starts_with("allocation-seal-") || (control.starts_with("allocation-motion-") && control!="allocation-motion-count")) {
  J retry;(void)chemfix::coupled_step(S.A,1e-12,0,S.gas,{},{},{},retry,continuation,true,control.starts_with("allocation-motion-")?"motion-empty-flows":"");
  check(S.A.arena->full_time_boundary()->accepted_step==1 && !S.gas->step_open() && continuation && continuation->next_step()==1,
        "allocation refusal must permit exact accepted-pair retry");
  retry_accepted=true;retry_receipt_sha256=retry.at("sha256").get<std::string>();
 }
 if(control=="two-step") {
  auto old=continuation;const auto first_native=S.A.arena->state_prefix(),first_gas=S.gas->accepted_sha256();
  std::shared_ptr<const c::CoupledTrial::Continuation> missing;
  bool refused_missing=false;try{J discarded;(void)chemfix::coupled_step(S.A,1e-12,1,S.gas,{},{},{},discarded,missing);}catch(const std::invalid_argument&){refused_missing=true;}
  check(refused_missing && S.A.arena->state_prefix()==first_native && S.gas->accepted_sha256()==first_gas && !S.gas->step_open(),"unattested continuation must refuse without accepted changes");
  J second;(void)chemfix::coupled_step(S.A,1e-12,1,S.gas,{},{},{},second,continuation);
  check(continuation && continuation->next_step()==2 && second.at("parent_common_receipt")==old->common_receipt_sha256() &&
        S.A.arena->full_time_boundary()->accepted_step==2,"actual authenticated second common step");
  const auto second_native=S.A.arena->state_prefix(),second_gas=S.gas->accepted_sha256();
  bool stale=false;try{J discarded;(void)chemfix::coupled_step(S.A,1e-12,2,S.gas,{},{},{},discarded,old);}catch(const std::invalid_argument&){stale=true;}
  check(stale && S.A.arena->state_prefix()==second_native && S.gas->accepted_sha256()==second_gas && !S.gas->step_open(),"successful publication invalidates prior accepted pair");
  Scenario foreign(blob,packets);bool refused_foreign=false;
  try{J discarded;(void)chemfix::coupled_step(foreign.A,1e-12,0,foreign.gas,{},{},{},discarded,continuation);}catch(const std::invalid_argument&){refused_foreign=true;}
  check(refused_foreign && foreign.A.arena->full_time_boundary()->accepted_step==0 && !foreign.gas->step_open(),"foreign equal-content pair refuses");
  record=J{{"first",record},{"second",second},{"continuation",continuation->common_receipt_sha256()}};
 }
 J phases=J::object();for(const auto&[name,value]:R.phases){phases[name]=J{{"wall_s",value.wall}};phases[name]["cpu_s"]=name.starts_with("motion_")?J():J(value.cpu);}
 J result={{"passed",true},{"control",control},{"checks",checks},{"refused",refused},{"refusal",error},{"common_accepted_publication",accepted},{"native_source",n::PicArena::source_identity()},
  {"timings_scope","concurrent-load component diagnostic, no isolated speedup claim"},{"phase_seconds",phases},{"census_carriers",R.census},{"retry_accepted",retry_accepted},{"retry_receipt_sha256",retry_receipt_sha256},{"seal_allocations",seal_allocations},{"motion_allocations",chemfix::measured_motion_allocations},{"strict_target_owner_transaction",true},{"production_qualified",false},{"record",record}};
 std::ofstream(out/"result.json")<<result.dump(2)<<"\n";std::cout<<result.dump()<<"\n";return 0;
 }catch(const std::exception&e){allocation_probe::fail_after=-1;std::cerr<<e.what()<<"\n";return 1;}}
