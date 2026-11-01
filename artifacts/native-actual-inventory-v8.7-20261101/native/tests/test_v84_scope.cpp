// Native-owned minimum scope gate: manufactured populations, no physical validation.
#include "chem_v8_fixture.hpp"
#include "raikou/mfem_pic_full_time_device.h"
#include "raikou/json.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace m=raikou::mfem_backend;
using J=nlohmann::json;
unsigned checks=0;
void check(bool b,const char *why){++checks;if(!b)throw std::runtime_error(why);}
template<class F> void refuses(F f,m::FullTimeChemistryRefusal::Kind kind,const char *why){
  bool hit=false;try{f();}catch(const m::FullTimeChemistryRefusal &e){hit=e.kind()==kind;}
  check(hit,why);
}
template<class F> void any_refusal(F f,const char *why){bool hit=false;try{f();}catch(const std::exception &){hit=true;}check(hit,why);}
m::FullTimeChemistryConfiguration scoped(m::FullTimeChemistryConfiguration c,const m::GapPicOptions &o){
  c.water_scope_contract=m::full_time_minimum_water_scope_contract(c,o);
  c.water_scope_contract_sha256=chemfix::sha(c.water_scope_contract);
  c.water_scope_mode=m::FullTimeWaterScopeMode::minimum_water;return c;
}
std::vector<m::FullTimeAuthoredPacket> pair(const std::string &target,double ux=1e6,double x2=.00215){
  return {{"load/e","pfc-v3/e","manufactured/scope",m::ArenaCount(1),{.00215,.00215,-.01},{ux,0,0}},
          {"load/target","pfc-v3/"+target,"manufactured/scope",m::ArenaCount(2),{x2,.00215,-.01},{0,0,0}}};
}
void refusal_unchanged(chemfix::Arena &A,m::FullTimeChemistryRefusal::Kind kind,const char *why,double dt=1e-15){
  const auto a=A.arena->full_time_boundary()->accepted_arena_sha256;
  const auto b=A.arena->full_time_chemistry().state_sha256;
  refuses([&]{chemfix::step(A,dt,0);},kind,why);
  check(A.arena->full_time_boundary()->accepted_arena_sha256==a && A.arena->full_time_chemistry().state_sha256==b,
        "refusal changed accepted population/history/RNG");
}
int main(int argc,char **argv){try{
  if(argc!=2)throw std::runtime_error("usage: test_v84_scope OUT");
  std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
  const auto o=chemfix::options();const auto base=chemfix::chemistry_configuration(o);const auto c=scoped(base,o);
  using K=m::FullTimeChemistryRefusal::Kind;
  for(const auto *name:{"OH","O","H2"}){
    auto A=chemfix::make_arena(o,c,pair(name),1000,2);refusal_unchanged(A,K::network_scope,"missing required neutral encounter not refused");
    // Diagnostic legacy behavior is preserved: missing groups are not admitted by the new mode accidentally.
    auto L=chemfix::make_arena(o,base,pair(name),1000,2);auto r=chemfix::step(L,1e-15,0);
    check(r.stage1->water_scope_binding_sha256.empty() && r.stage2->water_scope_binding_sha256.empty(),"legacy scope changed");
  }
  auto U=chemfix::make_arena(o,c,pair("OH+",0),1000,2);
  refusal_unchanged(U,K::network_scope,"wholly unknown recombination must refuse even cold-at-rest");
  auto separated=chemfix::make_arena(o,c,pair("OH",0,.00255),1000,2);
  check(chemfix::step(separated,1e-15,0).stage1->water_scope_binding_sha256.size()==64,"different-cell gap should be inactive");
  auto converge=chemfix::make_arena(o,c,pair("OH",4e7,.00255),1000,2);
  const auto converge_before=converge.arena->full_time_boundary()->accepted_arena_sha256;
  const auto converge_rng=converge.arena->full_time_chemistry().state_sha256;bool stage1_passed=false;
  refuses([&]{chemfix::step(converge,1e-11,0,[&](auto &){stage1_passed=true;});},K::network_scope,
          "new stage2 encounter after motion not refused");
  check(stage1_passed && converge.arena->full_time_boundary()->accepted_arena_sha256==converge_before &&
        converge.arena->full_time_chemistry().state_sha256==converge_rng,"stage2 refusal did not preserve outer rollback");
  auto one=pair("OH");one.pop_back();auto absent=chemfix::make_arena(o,c,one,1000,2);
  check(chemfix::step(absent,1e-15,0).stage1->water_scope_census_sha256.size()==64,"absent target handling");
  // Genuine all-row mapping cannot be replaced by the current target flags.
  auto unknown=pair("OH");unknown[1].state="unmapped/metal";
  any_refusal([&]{(void)chemfix::make_arena(o,c,unknown,1000,2);},"unregistered authored native state refused");
  auto unmapped=chemfix::make_arena(o,c,{},1000,2);
  auto extra=*std::find_if(unmapped.states.begin(),unmapped.states.end(),[](const auto &s){return s.id=="pfc-v3/OH";});
  extra.id="unmapped/metal";unmapped.states.push_back(extra);
  unmapped.arena=m::PicArena::create_full_time(o,unmapped.states,unmapped.channels,unknown,unmapped.g,unmapped.cfg);
  refusal_unchanged(unmapped,K::network_scope,"registered native state silently omitted from chemistry census");
  for(auto execution:{m::FullTimeChemistryExecution::device,m::FullTimeChemistryExecution::device_verified}){
    auto d=c;d.execution=execution;auto D=chemfix::make_arena(o,d,pair("H"),1000,2);
    refusal_unchanged(D,K::network_scope,"strict device mode fell back to host");
  }
  for(auto execution:{m::FullTimeExecution::device,m::FullTimeExecution::device_verified}){
    auto D=chemfix::make_arena(o,c,pair("H"),1000,2);
    refuses([&]{D.arena->attach_full_time_device({},execution);},K::network_scope,
            "strict host chemistry accepted a native device lane request");
    check(D.arena->full_time_execution()==m::FullTimeExecution::host,"refused device attach changed native lane");
  }
  // Scope support is not permission to bypass the existing region policy.
  auto strict=base;strict.policy.admitted_regions.clear();strict=scoped(strict,o);
  auto C=chemfix::make_arena(o,strict,pair("H2O",0),1000,2);
  refusal_unchanged(C,K::coverage_policy,"existing cold unknown-region guard relaxed");
  for(int k=0;k<4;++k){auto bad=c;
    if(k==0)bad.water_scope_contract+="weakened\n";
    if(k==1)bad.water_scope_contract_sha256=chemfix::sha("wrong");
    if(k==2)bad.run_seed++;
    if(k==3)bad.water_scope_mode=m::FullTimeWaterScopeMode(99);
    any_refusal([&]{(void)m::full_time_chemistry_registry(bad,o);},"scope/base configuration tamper admitted");
  }
  for(int k=0;k<2;++k){auto empty=chemfix::make_arena(o,base,{},1000,2);auto disabled=empty.cfg;
    disabled.chemistry.enabled=false;
    if(k==0)disabled.chemistry.water_scope_mode=m::FullTimeWaterScopeMode::minimum_water;
    else disabled.chemistry.water_scope_contract="stray contract";
    refuses([&]{(void)m::PicArena::create_full_time(o,empty.states,empty.channels,{},empty.g,disabled);},
            K::configuration,"disabled chemistry silently lost scope request");
  }
  // Two actual native stages and repeated restart; inherited physical fields are synthetic.
  auto A=chemfix::make_arena(o,c,pair("H"),1000,2);auto B=chemfix::make_arena(o,c,pair("H"),1000,2);
  for(int k=0;k<3;++k){
    const auto r=chemfix::step(A,1e-15,k);const auto q=chemfix::step(B,1e-15,k);
    check(r.stage1->water_scope_contract_sha256==c.water_scope_contract_sha256 &&
          r.stage2->water_scope_contract_sha256==c.water_scope_contract_sha256,"stage scope contract lost");
    check(r.stage1->water_scope_binding_sha256.size()==64 && r.stage2->water_scope_binding_sha256.size()==64 &&
          r.stage1->water_scope_binding_sha256!=r.stage2->water_scope_binding_sha256,"stage2 reused stage1 binding");
    check(r.boundary_sha256==q.boundary_sha256 && r.stage2->state_after_sha256==q.stage2->state_after_sha256,
          "restored continuation differs");
    const auto st=B.arena->full_time_chemistry();const auto boundary=*B.arena->full_time_boundary();
    const auto path=out/("scope-"+std::to_string(k)+".json");B.arena->save(path.string());
    B.arena=m::PicArena::restore_full_time(path.string(),o,B.states,B.channels,boundary);
    const auto restored=B.arena->full_time_chemistry();
    check(st.water_scope_history_sha256.size()==64 && st.water_scope_history_sha256==restored.water_scope_history_sha256 &&
          st.water_scope_census_sha256==restored.water_scope_census_sha256 && st.water_scope_binding_sha256==restored.water_scope_binding_sha256,
          "scope roots lost on restore");
    J doc;std::ifstream(path)>>doc;auto payload=doc.at("payload");
    payload["configuration"]["v8"]["chemistry"]["water_scope_contract_sha256"]=chemfix::sha("tampered scope");
    J altered={{"payload",payload},{"sha256",chemfix::sha(payload.dump())}};auto wrong=out/("wrong-"+std::to_string(k)+".json");
    std::ofstream(wrong)<<altered.dump();
    any_refusal([&]{(void)m::PicArena::restore_full_time(wrong.string(),o,B.states,B.channels,boundary);},"re-signed scope checkpoint tamper accepted");
    auto changed=doc.at("payload");auto hex=changed["metadata"]["chemistry"]["state_hex"].get<std::string>();
    static constexpr char digits[]="0123456789abcdef";std::string needle;
    for(unsigned char v:st.water_scope_binding_sha256){needle+=digits[v>>4];needle+=digits[v&15];}
    const auto pos=hex.find(needle);check(pos!=std::string::npos,"checkpoint must contain latest native-held scope binding");
    // Change one ASCII hex digit to another valid digest digit, then re-sign only the envelope.
    hex.replace(pos,2,hex.substr(pos,2)=="30"?"31":"30");changed["metadata"]["chemistry"]["state_hex"]=hex;
    auto wrong_root=out/("wrong-root-"+std::to_string(k)+".json");
    std::ofstream(wrong_root)<<J{{"payload",changed},{"sha256",chemfix::sha(changed.dump())}}.dump();
    any_refusal([&]{(void)m::PicArena::restore_full_time(wrong_root.string(),o,B.states,B.channels,boundary);},
                "re-signed native scope binding tamper accepted");
  }
  J result={{"passed",true},{"checks",checks},{"actual_native_fixture",true},{"scope_contract_sha256",c.water_scope_contract_sha256},
            {"network_complete",false},{"physical_model_validated",false},{"device_census_qualified",false},{"production_qualified",false}};
  std::ofstream(out/"result.json")<<result.dump(2)<<'\n';std::cout<<result.dump(2)<<'\n';return 0;
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
