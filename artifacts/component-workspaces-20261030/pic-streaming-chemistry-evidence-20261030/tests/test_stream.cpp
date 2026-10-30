#include "streaming_coverage.hpp"
#include "generic_water_bridge.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/resource.h>
namespace e=raikou::chemistry_evidence;namespace g=raikou::generic_water;using J=e::Json;
std::size_t checks=0;
void check(bool x,const char*s){++checks;if(!x)throw std::runtime_error(s);}
J read(const std::string&p){std::ifstream f(p);J j;f>>j;if(!f)throw std::runtime_error("input read");return j;}
void write(const std::string&p,const J&j){std::ofstream f(p);f<<j.dump(2)<<'\n';if(!f)throw std::runtime_error("output write");}
std::string bits(double x){std::ostringstream s;s<<std::hex<<std::bit_cast<std::uint64_t>(x);return s.str();}
std::vector<std::string> split(const std::string&s){std::vector<std::string> out;std::istringstream f(s);std::string x;while(std::getline(f,x,','))out.push_back(x);return out;}
J sample_row(const e::Query&q,const J&coverage){return {{"ordinal",q.ordinal},{"id",q.id},{"context_sha256",q.context_sha256},{"energy_eV",q.energy_eV},{"weight",q.weight},{"energy_bits",bits(q.energy_eV)},{"weight_bits",bits(q.weight)},{"provider_canonical_json",coverage.dump()}};}
template<class F>void refuses(e::Stream&s,F f,const char*why){const auto before=s.snapshot();bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,why);check(s.snapshot()==before,"refusal mutated stream");}
int main(int argc,char**argv){try{
 if(argc!=6)return 2;
 const std::string mode=argv[1],output=argv[5];auto data=std::make_shared<const raikou::kinetics::StrictCollisionDataset>(raikou::kinetics::StrictCollisionDataset::read(argv[2]));const auto run=read(argv[3]);J bindings=run.at("bindings");
 constexpr double amu=1.66053906660e-27,M=18.01056468403*amu,H=1.00782503223*amu;
 g::Bridge bridge(data,bindings.at("dataset_sha256"),bindings.at("kernel_manifest_sha256"),M,H,g::SitePolicy::PreserveReactantSites);const auto reference=bridge.coverage(0);
 if(mode=="small"){
  const auto existing=read(argv[4]);J inputs=J::array();e::Stream stream(reference,bindings),chunked(reference,bindings);std::size_t i=0;
  const std::array<double,7> weights={.03125,1.,2.,0.,.5,1e-20,7.};
  for(const auto&coverage:existing){const double E=coverage.at("energy_eV");check(bridge.coverage(E)==coverage,"old pointwise fixture differs");e::Query q{i,"existing-fixture/"+std::to_string(i),std::string(64,'a'),E,weights.at(i)};stream.observe(q,coverage);inputs.push_back(sample_row(q,coverage));++i;}
  // Chunk boundaries are a caller loop choice, with no aggregation reset.
  for(std::size_t start=0;start<inputs.size();start+=3)for(std::size_t n=start;n<std::min(inputs.size(),start+3);++n){const auto&r=inputs[n];chunked.observe({r.at("ordinal"),r.at("id"),r.at("context_sha256"),r.at("energy_eV"),r.at("weight")},J::parse(r.at("provider_canonical_json").get<std::string>()));}
  check(stream.snapshot()==chunked.snapshot(),"chunking changes evidence");
  // Old admitted candidate batch: preserve every row's scalar operand and the
  // frozen provider coverage, rather than replacing its physics callbacks.
  std::vector<g::PhysicalPacket> electrons;for(std::size_t n=0;n<16;++n){const double E=20+.1*n,x=E*g::qe/(g::me*g::c*g::c),u=g::c*std::sqrt(x*(x+2));electrons.push_back({{"e/"+std::to_string(n),"electron",7,.125*(n+1),g::me,{u,0,0}},{.004,.004,.0125}});}
  const std::vector<g::PhysicalPacket> water={{{"water","H2O",7,7,M,{}},{.004,.004,.0125}}};
  const auto legacy=bridge.candidates(electrons,water,{{7,1e-12}},g::Scope::CompareSupportedSubset);e::Stream actual(reference,bindings);J batch_inputs=J::array();
  for(std::size_t n=0;n<electrons.size();++n){const double E=g::pc::relative(electrons[n].packet,water[0].packet).target_frame_energy_eV;const auto cov=bridge.coverage(E);const auto found=std::find_if(legacy.coverage.begin(),legacy.coverage.end(),[&](const auto&r){return r.at("energy_eV")==E;});check(found!=legacy.coverage.end()&&*found==cov,"legacy batch evidence differs");e::Query q{n,electrons[n].packet.id,std::string(64,'b'),E,electrons[n].packet.weight};actual.observe(q,cov);batch_inputs.push_back(sample_row(q,cov));}
  const auto after=bridge.candidates(electrons,water,{{7,1e-12}},g::Scope::CompareSupportedSubset);check(after.rows.size()==legacy.rows.size(),"candidate count changed");
  for(std::size_t n=0;n<legacy.rows.size();++n){const auto&a=legacy.rows[n];const auto&b=after.rows[n];check(a.channel_id==b.channel_id&&a.projectile_id==b.projectile_id&&a.target_id==b.target_id,"candidate identity changed");check(bits(a.rate_s)==bits(b.rate_s)&&bits(a.sigma_m2)==bits(b.sigma_m2)&&bits(a.energy_eV)==bits(b.energy_eV)&&bits(a.event_weight)==bits(b.event_weight)&&bits(a.flux_m_s)==bits(b.flux_m_s)&&bits(a.volume_m3)==bits(b.volume_m3)&&bits(a.projectile_remainder.high)==bits(b.projectile_remainder.high)&&bits(a.projectile_remainder.low)==bits(b.projectile_remainder.low)&&bits(a.target_remainder.high)==bits(b.target_remainder.high)&&bits(a.target_remainder.low)==bits(b.target_remainder.low),"candidate numeric bytes changed");}
  g::pc::Scheduler clock1(811),clock2(811);g::pc::Clock clock{0,0,0,1e-10,std::string(64,'a'),std::string(64,'b'),std::string(64,'c'),0,0};const auto p1=clock1.prepare(clock,legacy.rows),p2=clock2.prepare(clock,after.rows);check(clock1.try_publish(p1)&&clock2.try_publish(p2)&&clock1.checkpoint()==clock2.checkpoint(),"rate/RNG parity changed");
  auto repeated=electrons;repeated.resize(3);for(auto&x:repeated)x.packet.proper_velocity_m_s=electrons[0].packet.proper_velocity_m_s;
  const auto old_repeated=bridge.candidates(repeated,water,{{7,1e-12}},g::Scope::CompareSupportedSubset);check(old_repeated.coverage.size()==1,"duplicate fixture did not exercise legacy energy deduplication");e::Stream duplicates(reference,bindings);J duplicate_inputs=J::array();
  for(std::size_t n=0;n<repeated.size();++n){const double E=g::pc::relative(repeated[n].packet,water[0].packet).target_frame_energy_eV;const auto cov=bridge.coverage(E);check(cov==old_repeated.coverage[0],"repeated-energy source evidence changed");e::Query q{n,"repeated/"+std::to_string(n),std::string(64,'d'),E,repeated[n].packet.weight};duplicates.observe(q,cov);duplicate_inputs.push_back(sample_row(q,cov));}
  check(duplicates.queries()==3,"repeated queries were silently deduplicated");write(output+"/duplicate-inputs.json",duplicate_inputs);write(output+"/duplicate-receipt.json",duplicates.snapshot());
  // Content commitment includes query order/weight and all provider details.
  auto changed_inputs=inputs;changed_inputs[0]["weight"]=.0625;e::Stream changed(reference,bindings);for(const auto&r:changed_inputs)changed.observe({r.at("ordinal"),r.at("id"),r.at("context_sha256"),r.at("energy_eV"),r.at("weight")},J::parse(r.at("provider_canonical_json").get<std::string>()));check(changed.commitment()!=stream.commitment(),"weight missing from commitment");
  const e::Query q{0,"probe",std::string(64,'c'),20,1};const auto good=bridge.coverage(20);e::Stream invalid(reference,bindings);
  for(int kind=0;kind<10;++kind){auto j=good;auto bad=q;
   if(kind==0)j["missing_rates"].erase(0);
   if(kind==1)j["conditional_models"][0]["source_sigma_m2"]=INFINITY;
   if(kind==2)j["complete_mechanism"]=true;
   if(kind==3)j["conditional_models"][0]["selected_subset"]=false;
   if(kind==4)bad.ordinal=1;
   if(kind==5)bad.energy_eV=21;
   if(kind==6)bad.weight=NAN;
   if(kind==7)j["missing_rates"][1]["sigma_m2"]=0;
   if(kind==8)j["unapplied_observables"][0]["rate_sum_eligible"]=true;
   if(kind==9)j["conditional_models"].push_back(j["conditional_models"][0]);
   refuses(invalid,[&]{invalid.observe(bad,j);},"invalid evidence admitted");}
  {auto bad=q;bad.weight=std::numeric_limits<double>::max();refuses(invalid,[&]{invalid.observe(bad,good);},"moment overflow admitted");}
  {auto bad=q;bad.weight=std::numeric_limits<double>::denorm_min();refuses(invalid,[&]{invalid.observe(bad,good);},"positive moment underflow admitted");}
  {e::Stream limited(reference,bindings,{1,UINT64_MAX,262144});limited.observe(q,good);auto bad=q;bad.ordinal=1;refuses(limited,[&]{limited.observe(bad,good);},"query cap ignored");}
  {e::Stream limited(reference,bindings,{1,1,262144});refuses(limited,[&]{limited.observe(q,good);},"byte cap ignored");}
  std::fesetround(FE_UPWARD);refuses(invalid,[&]{invalid.observe(q,good);},"changed rounding environment admitted");std::fesetround(FE_TONEAREST);
  for(int kind=0;kind<2;++kind){auto ref=reference;auto bind=bindings;if(kind==0)ref["site_policy"]=std::string(513,'x');else bind["query_input_sha256"]=std::string(4097,'a');bool caught=false;try{e::Stream rejected(ref,bind);}catch(const std::exception&){caught=true;}check(caught,"oversized constructor input admitted");}
  check(invalid.queries()==0,"invalid queries leaked count");
  write(output+"/small-inputs.json",inputs);write(output+"/small-receipt.json",stream.snapshot());write(output+"/legacy-batch-inputs.json",batch_inputs);write(output+"/legacy-batch-receipt.json",actual.snapshot());write(output+"/small-results.json",{{"status","passed"},{"assertions",checks},{"old_unique_points",7},{"old_admitted_batch_queries",16},{"legacy_candidate_rows",legacy.rows.size()},{"physics_or_RNG_changed",false}});std::cout<<checks<<" small-fixture and refusal assertions\n";
 }else if(mode=="capture"){
  e::Stream stream(reference,bindings);J checkpoints=J::array();std::ifstream f(argv[4]);std::string line;std::getline(f,line);const auto names=split(line);std::map<std::string,std::size_t> col;for(std::size_t n=0;n<names.size();++n)col[names[n]]=n;g::pc::Packet electron,target;target.mass_kg=M;std::uint64_t count=0;
  while(std::getline(f,line)){const auto row=split(line);check(row.size()==names.size(),"captured CSV field count");electron.mass_kg=std::stod(row.at(col.at("mass_kg")));electron.proper_velocity_m_s={std::stod(row.at(col.at("ux_m_s"))),std::stod(row.at(col.at("uy_m_s"))),std::stod(row.at(col.at("uz_m_s")))};const double E=g::pc::relative(electron,target).target_frame_energy_eV,W=std::stod(row.at(col.at("weight")));stream.observe({count,"captured-pre-kick/electron/"+row.at(col.at("id")),bindings.at("query_input_sha256"),E,W},bridge.coverage(E));++count;
   if(count==256||count==1024||count==4096||count==16384||count==46726){struct rusage usage{};getrusage(RUSAGE_SELF,&usage);checkpoints.push_back({{"queries",count},{"retained_keys",stream.retained_keys()},{"retained_category_slots",stream.retained_keys()*7},{"peak_RSS_KiB",usage.ru_maxrss},{"content_sha256",stream.commitment()}});}}
  check(count==46726&&stream.queries()==count,"captured query inventory mismatch");check(stream.retained_keys()==40,"evidence inventory grew/changed");write(output+"/capture-receipt.json",stream.snapshot());write(output+"/memory-checkpoints.json",checkpoints);write(output+"/capture-results.json",{{"status","passed"},{"queries",count},{"retained_keys",stream.retained_keys()},{"retained_category_slots",280},{"stationary_H2O_target_assumption",true},{"input_velocity_columns","pre-kick ux_m_s,uy_m_s,uz_m_s"},{"new_EEDF_or_evolution",false},{"physical_time_advanced",false}});std::cout<<count<<" captured scalar queries, fixed40-key/280-category storage\n";
 }else throw std::runtime_error("unknown mode");
}catch(const std::exception&x){std::cerr<<x.what()<<'\n';return 1;}}
