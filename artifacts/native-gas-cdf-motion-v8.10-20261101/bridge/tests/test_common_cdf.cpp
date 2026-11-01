#include "coupled_step_fixture.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace n=raikou::mfem_backend;namespace g=pf_gas_v1;namespace a=native_gas_stage;namespace t=native_gas_targets;namespace c=native_gas_common;
using J=nlohmann::json;unsigned checks=0;
void check(bool b,const char*s){++checks;if(!b)throw std::runtime_error(s);}
const char*blob_path="/home/treefor/projects/runtimes/chemistry-dynamic-cdf-v1.1-20261101/data/chemistry-cdf.blob";
const char*blob_sha="8435e873b7c68146fb4c429d96c3e8051334c22525fd508729299e8a2d65ea0f";
g::Configuration gas_config(const n::GapPicOptions&o,const n::FullTimeChemistryRegistry&r,const pfc::Blob&b){
 g::Configuration x;x.preset="diagnostic-CDF-coupled-bookkeeping";x.feedback=true;x.transport=g::Transport::census_conservative;
 const auto grid=n::chem8::make_grid(o,{std::uint32_t(o.nx),std::uint32_t(o.ny),5});x.grid={grid.low,grid.high,grid.n,grid.volume};
 x.chemistry_blob_sha256=b.sha256;x.chemistry_registry_sha256=r.registry_sha256;
 for(const auto&s:r.species)x.species.push_back({s.state,s.mass_kg,s.charge_e,s.target});
 for(std::uint32_t k=0;k<b.n_proc;++k){const auto&p=b.view.procs[k];x.processes.push_back({b.view.groups[p.group].projectile,p.n_out,a::internal_carrier(p.kind,p.n_out),p.kind});}return x;
}
std::vector<n::FullTimeAuthoredPacket> population(const pfc::Blob&b,double E,double W,unsigned count){
 auto it=std::find(b.species_names.begin(),b.species_names.end(),"e");check(it!=b.species_names.end(),"electron registry");const double me=b.view.species[it-b.species_names.begin()].mass_eV;
 std::vector<n::FullTimeAuthoredPacket> out;const double u=n::pic_c*std::sqrt(E*(E+2*me))/me;
 for(unsigned i=0;i<count;++i)out.push_back({"cdf/e/"+std::to_string(i),"pfc-v3/e","diagnostic-CDF-bookkeeping",n::ArenaCount(1.,0),{.0041,.0041,.0001},{0,0,u}});
 out.push_back({"cdf/H2","pfc-v3/H2","diagnostic-CDF-bookkeeping",n::ArenaCount(W,0),{.0041,.0041,.0001},{0,0,0}});return out;
}
struct Scenario {
 n::GapPicOptions o=chemfix::options();n::FullTimeChemistryConfiguration cfg;std::shared_ptr<g::Owner>gas;chemfix::Arena arena;
 Scenario(const pfc::Blob&b,const std::vector<n::FullTimeAuthoredPacket>&rows,bool strict_scope=false)
 :cfg(chemfix::chemistry_configuration(o)),gas([&]{cfg.blob_path=blob_path;cfg.blob_sha256=blob_sha;cfg.diagnostic_field_target_approximation=false;cfg.release_weight=1;
   cfg.policy.admitted_regions.clear();const auto r=n::full_time_chemistry_registry(cfg,o);
   // Explicitly diagnostic encounter coverage, independent of the strict energy owner.
   if(!strict_scope)for(unsigned i=0;i<r.region_names.size();++i)if(!r.region_bounded[i])cfg.policy.admitted_regions.push_back(r.region_names[i]);
   if(strict_scope){cfg.water_scope_contract=n::full_time_minimum_water_scope_contract(cfg,o);cfg.water_scope_contract_sha256=t::sha(cfg.water_scope_contract);cfg.water_scope_mode=n::FullTimeWaterScopeMode::minimum_water;}
   return std::make_shared<g::Owner>(gas_config(o,n::full_time_chemistry_registry(cfg,o),b));}()),
 arena([&]{cfg.coupled_target_contract_sha256=c::coupled_contract(cfg,o,*gas);cfg.host_target_feedback=true;cfg.host_target_contract_sha256=t::contract(cfg,o,*gas);return chemfix::make_arena(o,cfg,rows,20000,2);}()){}
};
int main(int argc,char**argv){try {
 check(argc==2,"output required");const std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
 auto blob=pfc::load_blob(blob_path,blob_sha);check(blob.view.blob_version==2 && PFC_RNG_LAYOUT_VERSION==2,"explicit blob/RNG schema");
 Scenario s(blob,population(blob,100,1e18,128));J receipt;std::shared_ptr<const c::CoupledTrial::Continuation> next;
 const auto report=chemfix::coupled_step(s.arena,1e-13,0,s.gas,{},{},{},receipt,next);
 std::uint64_t events=0,releases=0;
 for(const auto&stage:{report.stage1,report.stage2}) {
  events+=stage->admitted_events;releases+=stage->released_packets;
  for(const auto&item:stage->admitted_export->events){pfc::EventRec e{};check(item.event_bytes.size()==sizeof e,"event ABI");std::memcpy(&e,item.event_bytes.data(),sizeof e);
    check(blob.view.procs[e.process].kind==32 && e.n_out==3 && e.e_internal==0 && item.consumes_tracked_target,"actual coupled CDF consumption/no internal sink");}
 }
 check(events>0 && releases>0 && next && next->next_step()==1,"actual populated CDF common publication and release");
 for(const auto&stage:receipt.at("stages"))check(stage.at("energy_exact_zero") && stage.at("momentum_exact_zero"),"actual CDF combined stock accounting");
 std::vector<std::string> refusals;
 for(bool scope:{false,true}) {
  Scenario bad(blob,population(blob,scope?100:7,scope?1e18:1e22,scope?1:3),scope);
  const auto native=bad.arena.arena->state_prefix(),rng=bad.arena.arena->full_time_chemistry().state_sha256,gas=bad.gas->accepted_sha256();
  auto snapshot=bad.arena.arena->snapshot();std::shared_ptr<const c::CoupledTrial::Continuation> none;std::string why;
  try{J discarded;(void)chemfix::coupled_step(bad.arena,1e-13,0,bad.gas,{},{},{},discarded,none);}catch(const std::exception&e){why=e.what();}
  check(!why.empty() && why.find(scope?"water":"CDF_SUPPORT")!=std::string::npos,"strict incomplete group/selected unsupported CDF refusal");
  check(bad.arena.arena->state_prefix()==native && bad.arena.arena->snapshot()==snapshot && bad.arena.arena->full_time_chemistry().state_sha256==rng && bad.gas->accepted_sha256()==gas && !bad.gas->step_open(),"whole native/gas rollback on CDF refusal");
  refusals.push_back(why);
 }
 J result={{"passed",true},{"checks",checks},{"native_source",n::PicArena::source_identity()},{"events",events},{"releases",releases},{"refusals",refusals},{"receipt",receipt},{"coverage_mode","explicit diagnostic missing-channel admission; strict water scope separately refuses"},{"production_qualified",false}};
 std::ofstream(out/"result.json")<<result.dump(2)<<'\n';std::cout<<result.dump()<<'\n';return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
