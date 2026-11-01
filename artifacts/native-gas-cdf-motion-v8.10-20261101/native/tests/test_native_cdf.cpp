#include "population_step_fixture.hpp"
#include "chem_v3_native.hpp"
#include "cdf_legacy_common.hpp"
#include "forms_gpu.hpp"
#include "hall_gpu_common.hpp"
#include "raikou/mfem_pic_full_time_device.h"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace n=raikou::mfem_backend;
using J=nlohmann::json;
unsigned checks=0;
void check(bool x,const char*s){++checks;if(!x)throw std::runtime_error(s);}
const char* CDF="/home/treefor/projects/runtimes/chemistry-dynamic-cdf-v1.1-20261101/data/chemistry-cdf.blob";
const char* SHA="8435e873b7c68146fb4c429d96c3e8051334c22525fd508729299e8a2d65ea0f";
auto config(const n::GapPicOptions&o,bool diagnostic=true){auto c=chemfix::chemistry_configuration(o);c.blob_path=CDF;c.blob_sha256=SHA;c.policy.admitted_regions.clear();c.release_weight=1;
 auto reg=n::full_time_chemistry_registry(c,o);if(diagnostic)for(unsigned i=0;i<reg.region_names.size();++i)if(!reg.region_bounded[i])c.policy.admitted_regions.push_back(reg.region_names[i]);return c;}
std::vector<n::FullTimeAuthoredPacket> population(double E,double W,unsigned count=128,bool mixed=false){std::vector<n::FullTimeAuthoredPacket>out;const double me=510998.95069;
 for(unsigned i=0;i<count;++i){const double energy=mixed && i==0?100:E;const double u=n::pic_c*std::sqrt(energy*(energy+2*me))/me;out.push_back({"cdf/e/"+std::to_string(i),"pfc-v3/e","diagnostic-CDF",n::ArenaCount(1.,0),{.0041,.0041,.0001},{0,0,u}});}
 out.push_back({"cdf/target","pfc-v3/H2","diagnostic-CDF",n::ArenaCount(W,0),{.0041,.0041,.0001},{0,0,0}});return out;}
J legacy(const n::GapPicOptions&o){const auto cfg=chemfix::chemistry_configuration(o);const auto setup=n::chem8::make_setup(cfg,o);const auto&b=*setup->reg->blob;std::uint64_t kinds[7]{};pfc::Sha256 total;
 for(unsigned case_id=0;case_id<3;++case_id){auto P=pfc_test::make_population(b,10000,31,1000+case_id);n::chem8::StageInputs in;in.carriers=P.carriers;in.cells=P.cells;auto&c=in.cfg;c.seed=1234567+case_id;c.step=42+case_id;c.dt_s=case_id==0?4e-11:case_id==1?5e-12:2e-10;c.trunc_sigmas=7;c.purpose=PFC_PURPOSE_MCC;c.max_candidates=4096;c.max_events=32;
  const auto native=n::chem8::run_host(*setup,in,3);pfc::StepResult ref;pfc::cpu_step(b,c,P.carriers,P.cells,ref);
  check(native.outs.size()==ref.outs.size() && std::memcmp(native.outs.data(),ref.outs.data(),ref.outs.size()*sizeof(pfc::CarrierOut))==0,"legacy native/shared carrier parity");
  for(std::size_t i=0;i<in.carriers.size();++i)for(unsigned j=0;j<ref.outs[i].n_events;++j){const auto&ev=ref.events[i*c.max_events+j];check(std::memcmp(&ev,&native.events[native.event_begin[i]+j],sizeof ev)==0,"legacy native/shared event parity");++kinds[b.view.procs[ev.process].kind];}
  total.text(pfc::step_digest(ref.outs,ref.events,c.max_events,P.carriers));total.text(pfc_test::ledger_digest(pfc_test::ledger_of(b,ref,c.max_events)));
 }
 const auto digest=total.hex();check(digest=="73bc2fa6b1c8f8d9e6ae599f52f0b233ee93780bce7a6122f847f31a04a34172","qualified RNG2 legacy record/ledger digest");for(unsigned k=1;k<7;++k)check(kinds[k]>0,"all legacy kinds exercised");return J{{"digest",digest},{"kinds",{kinds[1],kinds[2],kinds[3],kinds[4],kinds[5],kinds[6]}}};}
int main(int argc,char**argv){try{check(argc==2||argc==4,"output [NVRTC UUID]");const std::filesystem::path out=argv[1];std::filesystem::create_directories(out);const auto o=chemfix::options();const auto blob=pfc::load_blob(CDF,SHA);
 check(blob.view.blob_version==2 && blob.view.n_cdf==1 && blob.view.procs[87].kind==32,"actual versioned CDF descriptor");
 std::unique_ptr<hall_gpu::Session>session;std::unique_ptr<pf_v7::gpu::Forms>forms;std::shared_ptr<n::FullTimeDevice>dev;
 if(argc==4){hall_gpu::SessionOptions so;so.require_uuid=argv[3];so.nvrtc_library=argv[2];session=std::make_unique<hall_gpu::Session>(so);forms=std::make_unique<pf_v7::gpu::Forms>(*session);dev=std::make_shared<n::FullTimeDevice>(*session,*forms);}
 auto attach=[&](auto&A,bool gpu){if(gpu)A.arena->attach_full_time_device(dev,n::FullTimeExecution::device);};
 J result={{"legacy",legacy(o)}};std::vector<std::string>trace,mixed_refusals;
 for(bool gpu:{false,true}){if(gpu&&!dev)continue;auto c=config(o);if(gpu)c.execution=n::FullTimeChemistryExecution::device;
  auto A=chemfix::make_arena(o,c,population(100.,1e18),20000);attach(A,gpu);const auto run=chemfix::step(A,1e-13,0);const auto ne=run.stage1->admitted_events+run.stage2->admitted_events;check(ne>0,"actual native supported CDF events");
  for(const auto&s:{run.stage1,run.stage2})for(const auto&a:s->admitted_export->events){pfc::EventRec e{};check(a.event_bytes.size()==sizeof e,"native event ABI");std::memcpy(&e,a.event_bytes.data(),sizeof e);check(blob.view.procs[e.process].kind==32 && e.n_out==3 && e.e_internal==0 && a.consumes_tracked_target,"kind32 consumes H2 with no internal sink");check(blob.species_names[e.out_species[0]]=="e"&&blob.species_names[e.out_species[1]]=="H"&&blob.species_names[e.out_species[2]]=="H","kind32 exact product mapping");}
  const auto key=run.stage1->records_root+run.stage2->records_root+run.stage1->ledger_sha256+run.stage2->ledger_sha256+run.boundary_sha256;trace.push_back(key);
  result[gpu?"device_events":"host_events"]=ne;
  // Actual low-energy native stage must refuse, including raw earlier events.
  auto B=chemfix::make_arena(o,c,population(7.,1e22,3,true),20000);attach(B,gpu);const auto state=B.arena->state_prefix(),rng=B.arena->full_time_chemistry().state_sha256,coverage=B.arena->full_time_chemistry().coverage_json;
  auto snapshot=B.arena->snapshot();std::string why;std::uint64_t earlier=0;unsigned bad=0;
  auto prepare=[&](n::CollisionScope&s){if(!gpu){const auto rows=s.rate_population(20000);std::vector<n::GapParticle> p;std::vector<std::uint32_t>sp;std::vector<std::uint8_t>single;auto setup=n::chem8::make_setup(c,o);for(const auto&r:rows){p.push_back(r.native);sp.push_back(unsigned(std::find(setup->reg->states.begin(),setup->reg->states.end(),r.state)-setup->reg->states.begin()));single.push_back(1);}n::chem8::RowView view{p.data(),p.size(),sp.data(),single.data()};auto in=n::chem8::build_inputs(*setup,view,{0,1,.5e-13,20000});const auto rec=n::chem8::run_host(*setup,in,1);earlier=rec.events.size();for(const auto&v:rec.outs)bad+=v.status==PFC_ST_CDF_SUPPORT;}
   return s.prepare_chemistry_v3();};
  try{(void)chemfix::step_with_population(B,1e-13,0,prepare);}catch(const n::FullTimeChemistryRefusal&e){why=e.what();}
  check(why.find("status CDF_SUPPORT")!=std::string::npos,"selected unsupported CDF typed refusal");check(B.arena->state_prefix()==state&&B.arena->full_time_chemistry().state_sha256==rng&&B.arena->full_time_chemistry().coverage_json==coverage&&B.arena->snapshot()==snapshot,"whole-stage refusal preserves accepted particles/RNG/stocks/coverage");
  if(!gpu){check(earlier>0&&bad>0,"real earlier events coexist with selected unsupported refusal");result["discarded_earlier_raw_events"]=earlier;}
  result[gpu?"device_refusal":"host_refusal"]=why;
  // First bad carrier has CDF_SUPPORT8 while another has EVENT_CAP3.
  auto mc=c;mc.max_events=1;auto packets=population(7.,1e22,2);packets[1].proper_velocity_m_s=population(100.,1.,1)[0].proper_velocity_m_s;
  auto M=chemfix::make_arena(o,mc,packets,20000);attach(M,gpu);const auto mb=M.arena->state_prefix(),mr=M.arena->full_time_chemistry().state_sha256;auto ms=M.arena->snapshot();std::string mw;
  auto mp=[&](n::CollisionScope&scope){if(!gpu){const auto rows=scope.rate_population(20000);std::vector<n::GapParticle> p;std::vector<std::uint32_t>sp;std::vector<std::uint8_t>single;auto setup=n::chem8::make_setup(mc,o);for(const auto&r:rows){p.push_back(r.native);sp.push_back(unsigned(std::find(setup->reg->states.begin(),setup->reg->states.end(),r.state)-setup->reg->states.begin()));single.push_back(1);}n::chem8::RowView view{p.data(),p.size(),sp.data(),single.data()};auto in=n::chem8::build_inputs(*setup,view,{0,1,.5e-13,20000});const auto rec=n::chem8::run_host(*setup,in,1);bool has8=false,has3=false;for(const auto&x:rec.outs){has8|=x.status==PFC_ST_CDF_SUPPORT;has3|=x.status==PFC_ST_EVENT_CAP;}check(has8&&has3&&rec.outs.front().status==PFC_ST_CDF_SUPPORT,"actual mixed-status first carrier fixture");}
   return scope.prepare_chemistry_v3();};
  try{(void)chemfix::step_with_population(M,1e-13,0,mp);}catch(const n::FullTimeChemistryRefusal&e){mw=e.what();}
  check(mw.find("status CDF_SUPPORT")!=std::string::npos,"first bad carrier retains its actual status8");check(M.arena->state_prefix()==mb&&M.arena->full_time_chemistry().state_sha256==mr&&M.arena->snapshot()==ms,"mixed-status refusal leaves accepted state unchanged");mixed_refusals.push_back(mw);
 }
 if(dev){check(trace.size()==2&&trace[0]==trace[1],"actual native CPU/device records ledgers accepted identity");check(mixed_refusals.size()==2&&mixed_refusals[0]==mixed_refusals[1],"mixed-status CPU/device exact carrier/status parity");}
 // Checkpoint/restart carries the new configuration/codec/RNG identity.
 {auto c=config(o);auto A=chemfix::make_arena(o,c,population(100.,1e17,16),20000,2);auto B=chemfix::make_arena(o,c,population(100.,1e17,16),20000,2);
  for(int k=0;k<3;++k){const auto a=chemfix::step(A,1e-13,k),b=chemfix::step(B,1e-13,k);
   check(a.boundary_sha256==b.boundary_sha256 && a.stage2->state_after_sha256==b.stage2->state_after_sha256,"CDF repeated restart equals uninterrupted continuation");
   const auto boundary=*B.arena->full_time_boundary();const auto state=B.arena->full_time_chemistry().state_sha256;const auto path=out/("cdf-restart-"+std::to_string(k)+".json");B.arena->save(path.string());
   B.arena=n::PicArena::restore_full_time(path.string(),o,B.states,B.channels,boundary);check(B.arena->full_time_chemistry().state_sha256==state,"CDF stored RNG/reservoir state survives restore");
   J doc;std::ifstream(path)>>doc;auto payload=doc.at("payload");payload["configuration"]["v8"]["identity"]="v8.7-actual-reservoir-stock-v1";J tamper={{"payload",payload},{"sha256",chemfix::sha(payload.dump())}};
   const auto bad=out/("cdf-old-schema-"+std::to_string(k)+".json");std::ofstream(bad)<<tamper.dump();bool refused=false;try{(void)n::PicArena::restore_full_time(bad.string(),o,B.states,B.channels,boundary);}catch(const std::exception&){refused=true;}check(refused,"old chemistry schema cannot masquerade as CDF/RNG2");
  }
 }
 // Strict incomplete-pair scope remains a hard gate with the new group present.
 {auto c=config(o);const auto reg=n::full_time_chemistry_registry(c,o);auto setup=n::chem8::make_setup(c,o);water_scope::Catalogue cat(std::span(reinterpret_cast<const unsigned char*>(blob.base()),blob.bytes),blob.sha256);c.water_scope_mode=n::FullTimeWaterScopeMode::minimum_water;c.water_scope_contract=water_scope::minimum_contract(cat,reg.registry_sha256,setup->base_configuration_sha256);c.water_scope_contract_sha256=water_scope::sha256(c.water_scope_contract);
  auto A=chemfix::make_arena(o,c,population(100.,1e18,1),20000);const auto before=A.arena->state_prefix();bool refused=false;try{(void)chemfix::step(A,1e-13,0);}catch(const std::exception&e){refused=std::string(e.what()).find("water")!=std::string::npos;}check(refused&&A.arena->state_prefix()==before,"incomplete H2 scope stays refused");}
 result["passed"]=true;result["checks"]=checks;result["production_qualified"]=false;result["network_complete"]=false;result["source"]=n::PicArena::source_identity();std::ofstream(out/"result.json")<<result.dump(2)<<'\n';std::cout<<result.dump()<<'\n';dev.reset();forms.reset();if(session)(void)session->close();return 0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
