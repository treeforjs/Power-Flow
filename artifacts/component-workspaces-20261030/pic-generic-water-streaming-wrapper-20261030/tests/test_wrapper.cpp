#include "streaming_water_bridge.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <chrono>
namespace w=raikou::streaming_water;namespace g=raikou::generic_water;namespace pc=raikou::pic_collisions;namespace e=raikou::chemistry_evidence;using J=e::Json;
std::size_t checks=0,error_cases=0;
void check(bool x,const char*s){++checks;if(!x)throw std::runtime_error(s);}
J read(const std::string&p){std::ifstream f(p);J j;f>>j;if(!f)throw std::runtime_error("input read");return j;}
void write(const std::string&p,const J&j){std::ofstream f(p);f<<j.dump(2)<<'\n';if(!f)throw std::runtime_error("output write");}
std::string rng_state(const std::mt19937_64&r){std::ostringstream s;s<<r;return s.str();}
std::string candidate_bytes(const std::vector<pc::Candidate>&rows){std::string out;for(const auto&a:rows){for(const auto*s:{&a.projectile_id,&a.target_id,&a.channel_id,&a.data_sha256,&a.kernel_sha256})e::word(out,*s);e::u64(out,a.cell);for(double x:{a.event_weight,a.projectile_before,a.target_before,a.projectile_remainder.high,a.projectile_remainder.low,a.target_remainder.high,a.target_remainder.low,a.rate_s,a.energy_eV,a.flux_m_s,a.sigma_m2,a.volume_m3})e::u64(out,std::bit_cast<std::uint64_t>(x));}return out;}
J work(const pc::CandidateWork&x){return {{"input_string_bytes",x.input_string_bytes},{"cell_lookups",x.cell_lookups},{"target_index_entries",x.target_index_entries},{"channel_index_entries",x.channel_index_entries},{"projectiles_visited",x.projectiles_visited},{"pair_visits",x.pair_visits},{"kinematics_calls",x.kinematics_calls},{"channel_group_lookups",x.channel_group_lookups},{"channel_queries",x.channel_queries},{"candidate_rows",x.candidate_rows},{"candidate_string_bytes",x.candidate_string_bytes}};}
J outcome(const g::Outcome&o){J products=J::array();for(const auto&p:o.products)products.push_back({p.state.id,p.state.H,p.state.O,p.state.q,p.state.mass_kg,p.state.internal_J,p.state.reference_J,p.proper_velocity_m_s,p.at_primary,p.weight});return {{"channel",o.channel},{"model",o.model},{"data",o.data_sha256},{"kernel",o.kernel_sha256},{"weight",o.event_weight},{"energy",o.energy_eV},{"sigma",o.sigma_m2},{"products",products},{"capture",o.captures_primary},{"OH_level",o.OH_level},{"energy_defect",o.weighted_energy_defect_J},{"energy_bound",o.weighted_energy_bound_J},{"momentum_defect",o.weighted_momentum_defect},{"momentum_bound",o.weighted_momentum_bound},{"random_variates",o.random_variates},{"projectile_remainder",{o.projectile_remainder.high,o.projectile_remainder.low}},{"target_remainder",{o.target_remainder.high,o.target_remainder.low}}};}
template<class F>std::string error(F f){try{f();return "ADMITTED";}catch(const pc::Error&x){return std::string("pc::Error:")+x.what();}catch(const g::k::CollisionCoverageError&x){return std::string("CollisionCoverageError:")+x.what();}catch(const std::exception&x){return std::string("other:")+x.what();}}
constexpr double amu=1.66053906660e-27,M=18.01056468403*amu,H=1.00782503223*amu;
g::PhysicalPacket electron(double E,std::string id,double weight=.03125,std::uint64_t cell=7){const double x=E*g::qe/(g::me*g::c*g::c),u=g::c*std::sqrt(x*(x+2));return {{std::move(id),"electron",cell,weight,g::me,{.6*u,0,.8*u}},{.004,.004,.0125}};}
g::PhysicalPacket water(std::string id="water",double weight=7,std::uint64_t cell=7){return {{std::move(id),"H2O",cell,weight,M,{}},{.004,.004,.0125}};}
void clocks(const std::vector<pc::Candidate>&a,const std::vector<pc::Candidate>&b){long double sum=0;for(const auto&x:a)sum+=x.rate_s;const double dt=sum>0?static_cast<double>(1/sum):1;pc::Clock stage{0,0,0,dt,std::string(64,'a'),std::string(64,'b'),std::string(64,'c'),0,0};for(int seed=0;seed<32;++seed){pc::Scheduler left(seed+71),right(seed+71);const auto x=left.prepare(stage,a),y=right.prepare(stage,b);check(bool(x.event())==bool(y.event())&&x.event_subflow_elapsed_s()==y.event_subflow_elapsed_s(),"scheduler choice/wait parity");if(x.event())check(candidate_bytes({*x.event()})==candidate_bytes({*y.event()}),"selected candidate parity");check(left.try_publish(x)&&right.try_publish(y)&&left.checkpoint()==right.checkpoint(),"scheduler RNG parity");}}
int main(int argc,char**argv){try{
 if(argc!=5)return 2;
 const auto config=read(argv[2]);const std::string out=argv[3],mode=argv[4];auto data=std::make_shared<const g::k::StrictCollisionDataset>(g::k::StrictCollisionDataset::read(argv[1]));const auto bindings=config.at("stream_bindings");
 g::Bridge legacy(data,bindings.at("dataset_sha256"),bindings.at("kernel_manifest_sha256"),M,H,g::SitePolicy::PreserveReactantSites);
 w::Bridge fork(data,bindings.at("dataset_sha256"),bindings.at("kernel_manifest_sha256"),M,H,g::SitePolicy::PreserveReactantSites);
 w::EvidenceRequest request{bindings,config.at("scope_context_sha256"),config.at("wrapper_source_sha256"),config.at("observed_builder_source_sha256"),{}};
 const std::vector<pc::Cell> cells={{7,1e-12},{8,2e-12}};J cases=J::array(),failures=J::array();std::size_t projections=0,products=0;
 auto compare=[&](std::string name,const std::vector<g::PhysicalPacket>&es,const std::vector<g::PhysicalPacket>&ns,bool old_admits){
  const auto started=std::chrono::steady_clock::now();const auto actual=fork.candidates(es,ns,cells,g::Scope::CompareSupportedSubset,request);const double wrapper_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();std::vector<pc::Packet> pe,pn;for(const auto&p:es)pe.push_back(p.packet);for(const auto&p:ns)pn.push_back(p.packet);pc::CandidateWork physical_work;const auto reference=pc::candidates(pe,pn,cells,legacy.comparison_channels(g::Scope::CompareSupportedSubset),1000000,1000000,100000,pc::CandidateInputLimits{},&physical_work);
  check(candidate_bytes(actual.rows)==candidate_bytes(reference),"full physical candidate bytes/order differ");check(work(actual.work)==work(physical_work),"physical work counters differ");
  if(old_admits){const auto old=legacy.candidates(es,ns,cells,g::Scope::CompareSupportedSubset);check(candidate_bytes(actual.rows)==candidate_bytes(old.rows)&&work(actual.work)==work(old.work),"admitted old entrypoint parity");}
  auto expected_bindings=bindings;expected_bindings["query_input_sha256"]=e::sha256(actual.evidence.at("query_input_descriptor").dump());e::Stream expected(legacy.coverage(0),expected_bindings);std::uint64_t ordinal=0;
  for(const auto&p:es)for(const auto&n:ns)if(p.packet.cell==n.packet.cell){const auto cell=std::find_if(cells.begin(),cells.end(),[&](const auto&c){return c.id==p.packet.cell;});const auto kin=pc::relative(p.packet,n.packet);expected.observe({ordinal++,w::pair_id(p.packet,n.packet,*cell,kin),request.scope_context_sha256,kin.target_frame_energy_eV,p.packet.weight},legacy.coverage(kin.target_frame_energy_eV));++projections;}
  check(expected.snapshot()==actual.evidence.at("stream"),"cached complete evidence differs from frozen provider");clocks(actual.rows,reference);
  struct rusage usage{};getrusage(RUSAGE_SELF,&usage);cases.push_back({{"name",name},{"wrapper_wall_seconds",wrapper_seconds},{"projectiles",es.size()},{"targets",ns.size()},{"candidate_rows",actual.rows.size()},{"candidate_bytes_sha256",e::sha256(candidate_bytes(actual.rows))},{"work",work(actual.work)},{"evidence",actual.evidence},{"peak_RSS_KiB",usage.ru_maxrss},{"old_entrypoint_admitted",old_admits}});
  return actual;
 };
 if(mode=="small"){
  for(double E:{0.,6.5,8.5,20.,100.,1000.,8084.470604111418,12000.}){
   auto p=electron(E,"e");auto n=water();auto result=compare("energy/"+std::to_string(E),{p},{n},true);
   for(const auto&candidate:result.rows){std::mt19937_64 r1(93),r2(93);check(outcome(legacy.sample(candidate,p,n,r1))==outcome(fork.sample(candidate,p,n,r2))&&rng_state(r1)==rng_state(r2),"delegated product/RNG parity");
    const auto&specs=legacy.specs();const auto s=std::find_if(specs.begin(),specs.end(),[&](const auto&spec){return spec.channel==candidate.channel_id;});
    check(outcome(legacy.sample(candidate,p,n,std::vector<double>(s->variates,.25)))==outcome(fork.sample(candidate,p,n,std::vector<double>(s->variates,.25))),"delegated explicit-variates product parity");++products;}
  }
  compare("empty",{}, {water()},true);compare("sparse",{electron(20,"e0",1,7),electron(25,"e1",2,8)},{water("n0",3,7),water("n1",5,8)},true);
  compare("repeated and unequal",{electron(20,"e0",1),electron(20,"e1",2)},{water("n0",3),water("n1",5)},true);
  compare("maximum-length IDs",{electron(20,std::string(4096,'a'))},{water(std::string(4096,'b'))},true);
  // Exercise support edges using the actual relative-energy calculation; the
  // frozen provider remains the pointwise oracle even when rounding moves E.
  for(double edge:{4.5,6.25,6.75,9.,12.5,13.5,977.7,1000.,10000.})for(double E:{std::nextafter(edge,0.),edge,std::nextafter(edge,INFINITY)})
   compare("support-edge/"+std::to_string(E),{electron(E,"edge")},{water()},true);
  auto moving_water=water();moving_water.packet.proper_velocity_m_s={730.,-420.,310.};
  compare("moving-target",{electron(20,"moving",.125)},{moving_water},true);
  const auto es=std::vector<g::PhysicalPacket>{electron(20,"e0"),electron(24,"e1")};const auto ns=std::vector<g::PhysicalPacket>{water()};const auto baseline=fork.candidates(es,ns,cells,g::Scope::CompareSupportedSubset,request);
  for(int kind=0;kind<14;++kind){auto a=es;auto n=ns;auto c=cells;auto scope=g::Scope::CompareSupportedSubset;std::size_t pairs=1000000,queries=1000000,rows=100000;
   if(kind==0)a[0].packet.state="OH_X_v0";
   if(kind==1)a[0].packet.weight=0;
   if(kind==2)n[0].packet.mass_kg=g::me;
   if(kind==3)a[0].position_m[0]=INFINITY;
   if(kind==4)c.clear();
   if(kind==5)c.push_back(c[0]);
   if(kind==6)a[1].packet.id=a[0].packet.id;
   if(kind==7)a[0].packet.proper_velocity_m_s[0]=NAN;
   if(kind==8)c[0].volume_m3=0;
   if(kind==9)pairs=0;
   if(kind==10)queries=7;
   if(kind==11)queries=8;
   if(kind==12)rows=2;
   if(kind==13)scope=g::Scope::RequireCompleteMechanism;
   const auto old_error=error([&]{legacy.candidates(a,n,c,scope,pairs,queries,rows);});const auto new_error=error([&]{fork.candidates(a,n,c,scope,request,pairs,queries,rows);});check(old_error!="ADMITTED"&&new_error==old_error,"physical error parity changed");++error_cases;failures.push_back({{"case",kind},{"old",old_error},{"new",new_error},{"classification","physical-builder parity"}});
   const auto after=fork.candidates(es,ns,cells,g::Scope::CompareSupportedSubset,request);check(after.evidence==baseline.evidence&&candidate_bytes(after.rows)==candidate_bytes(baseline.rows),"failed build leaked private evidence/state");
  }
  // Approved stream-v1 numerical-domain limits are distinct from physical
  // builder errors. Never mislabel these earlier evidence refusals as parity.
  for(double weight:{std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::max()}){
   const auto a=std::vector<g::PhysicalPacket>{electron(20,"numeric-domain",weight)};const auto old_error=error([&]{legacy.candidates(a,ns,cells,g::Scope::CompareSupportedSubset);});const auto new_error=error([&]{fork.candidates(a,ns,cells,g::Scope::CompareSupportedSubset,request);});check(old_error=="ADMITTED"&&new_error!="ADMITTED","evidence numerical-domain difference not visible");failures.push_back({{"weight",weight},{"old",old_error},{"new",new_error},{"classification","new explicit evidence-domain refusal; not physical error parity"}});
  }
  for(int kind=0;kind<6;++kind){auto bad=request;
   if(kind==0)bad.scope_context_sha256="not-a-digest";
   if(kind==1)bad.stream_bindings["dataset_sha256"]=std::string(64,'f');
   if(kind==2)bad.stream_bindings["weight_measure"]="not projectile query weight";
   if(kind==3)bad.limits.queries=1; // Fail after first complete private pair.
   if(kind==4)bad.limits.record_bytes=1;
   if(kind==5)bad.limits.encoded_bytes=1;
   const auto refusal=error([&]{fork.candidates(es,ns,cells,g::Scope::CompareSupportedSubset,bad);});check(refusal!="ADMITTED","evidence contract/cap refusal missing");
   const auto after=fork.candidates(es,ns,cells,g::Scope::CompareSupportedSubset,request);check(after.evidence==baseline.evidence&&candidate_bytes(after.rows)==candidate_bytes(baseline.rows),"failed evidence preparation leaked state");
   failures.push_back({{"case",kind},{"new",refusal},{"classification","new evidence contract/cap refusal; private state discarded"}});
  }
  write(out+"/small-cases.json",cases);write(out+"/error-cases.json",failures);write(out+"/small-results.json",{{"status","passed"},{"assertions",checks},{"physical_error_cases",error_cases},{"full_projection_comparisons",projections},{"actual_product_comparisons",products},{"evidence_domain_differences",2},{"physical_kernel_changed",false},{"native_or_field_evolution",false}});
 }else if(mode=="broad"){
  std::vector<g::PhysicalPacket> es;for(std::size_t i=0;i<2048;++i)es.push_back(electron(5+4*i,"broad/"+std::to_string(i)));
  const auto old_error=error([&]{legacy.candidates(es,{water()},cells,g::Scope::CompareSupportedSubset,1000000,1000000,100000,1024);});check(old_error=="pc::Error:coverage evidence cap","old maximum did not reject broad evidence");
  auto result=compare("2048 distinct energies beyond old maximum",es,{water()},false);check(result.evidence.at("stream").at("queries")==2048&&result.evidence.at("stream").at("retained_category_slots")==280,"broad evidence grew/dropped queries");
  for(auto&p:es)p.packet.proper_velocity_m_s=es[0].packet.proper_velocity_m_s;
  auto repeated=compare("2048 equal energies admitted by legacy",es,{water()},true);check(repeated.evidence.at("last_energy_record_reuses")==2047,"bounded last-record memo not used");
  write(out+"/broad-cases.json",cases);write(out+"/broad-results.json",{{"status","passed"},{"assertions",checks},{"full_projection_comparisons",projections},{"new_entrypoint_broad_queries",2048},{"old_maximum_refusal",old_error},{"physical_reference","unchanged frozen indexed pc::candidates"},{"physical_kernel_changed",false},{"native_or_field_evolution",false}});
 }else throw std::runtime_error("unknown test mode");
 std::cout<<mode<<": "<<checks<<" assertions, "<<projections<<" full provider projections\n";
}catch(const std::exception&x){std::cerr<<x.what()<<'\n';return 1;}}
