#include "ledger_codec.hpp"
#include "common_step.hpp"
#include "request_fixture.hpp"
namespace lc=powerflow::common_ledger_codec_v1;
namespace cs=powerflow::common_step_v9;
int main(int argc,char**argv){
 if(argc!=2)return 2;J result={{"scope","actual neutral common ledger serialization plus labelled arithmetic corruption controls; no restored owner authority or common restart claim"},{"cases",J::array()},{"passed",true}};int failed=0;
 auto check=[&](const char*name,bool ok){result["cases"].push_back({{"name",name},{"passed",ok}});failed+=!ok;};
 auto refuses=[&](const char*name,auto action){bool refused=false;try{action();}catch(const std::exception&){refused=true;}check(name,refused);};
 try{
  auto run=cs::Run::create(request(argv[1]));ss::s7::BirthPolicy policy;policy.subcell_extent_m=20e-6;policy.return_extent_m=0;
  const auto states=n::full_time_chemistry_states(run->owners().registry);
  for(unsigned step=0;step<2;++step){
   auto next=cs::Step::prepare(run,1e-14,policy);check("actual common transition accepts",next.try_publish());const auto&ledger=*run->ledger();
   const auto bytes=lc::encode(ledger);auto restored=lc::decode(bytes,ledger.sha256,ledger.accepted_step,ledger.time_s,states);
   check("actual ledger exact canonical roundtrip",lc::encode(*restored)==bytes&&restored->sha256==ledger.sha256);
   check("actual cumulative endpoint survives exact codec",restored->accepted_energy_J==ledger.accepted_energy_J&&cs::same_inventory(restored->accepted_material,ledger.accepted_material)&&restored->admitted_bound_sum_J==ledger.admitted_bound_sum_J);
   refuses("wrong accepted ledger root",[&]{lc::decode(bytes,std::string(64,'0'),ledger.accepted_step,ledger.time_s,states);});
   refuses("wrong accepted step",[&]{lc::decode(bytes,ledger.sha256,ledger.accepted_step+1,ledger.time_s,states);});
   refuses("wrong accepted time",[&]{lc::decode(bytes,ledger.sha256,ledger.accepted_step,std::nextafter(ledger.time_s,INFINITY),states);});
   refuses("wrong actual species catalogue",[&]{auto wrong=states;wrong.front().id="unknown";lc::decode(bytes,ledger.sha256,ledger.accepted_step,ledger.time_s,wrong);});
   refuses("duplicate actual species catalogue",[&]{auto wrong=states;wrong.push_back(states.front());lc::decode(bytes,ledger.sha256,ledger.accepted_step,ledger.time_s,wrong);});
   refuses("trailing bytes",[&]{lc::decode(bytes+"x",ledger.sha256,ledger.accepted_step,ledger.time_s,states);});
   for(std::size_t cut:{std::size_t(0),std::size_t(7),std::size_t(8),std::size_t(40),bytes.size()/2,bytes.size()-1})refuses("truncated ledger",[&]{lc::decode(std::string_view(bytes).substr(0,cut),ledger.sha256,ledger.accepted_step,ledger.time_s,states);});
   refuses("oversized input before parsing",[&]{lc::decode(std::string(lc::byte_cap+1,'x'),ledger.sha256,ledger.accepted_step,ledger.time_s,states);});
   // First inventory exact digit begins after schema and four 64-bit clocks.
   const std::size_t first_digit=8+std::string("common-ledger-checkpoint-v1").size()+32;
   refuses("noncanonical signed lower digit",[&]{auto bad=bytes;for(unsigned j=0;j<8;++j)bad[first_digit+j]=char(255);lc::decode(bad,ledger.sha256,ledger.accepted_step,ledger.time_s,states);});
   refuses("overflowing signed top digit",[&]{auto bad=bytes;for(unsigned j=0;j<8;++j)bad[first_digit+8*(pfg::kDigits-1)+j]=char(j==7?127:255);lc::decode(bad,ledger.sha256,ledger.accepted_step,ledger.time_s,states);});
   refuses("rehashed inconsistent cumulative energy",[&]{auto bad=ledger;bad.external_net_J.add_product(1.,1.);bad.sha256=cs::ledger_hash(bad);lc::encode(bad);});
   refuses("rehashed inconsistent cumulative nuclei",[&]{auto bad=ledger;bad.accepted_material.hydrogen.add_product(1.,1.);bad.sha256=cs::ledger_hash(bad);lc::encode(bad);});
   refuses("rehashed negative cumulative bound",[&]{auto bad=ledger;bad.admitted_bound_sum_J={};bad.admitted_bound_sum_J.add_product(1.,1.,1);bad.sha256=cs::ledger_hash(bad);lc::encode(bad);});
   result["last_actual"]={{"step",ledger.accepted_step},{"bytes",bytes.size()},{"sha256",ledger.sha256},{"time_s",ledger.time_s}};
  }
  // Deliberately manufactured arithmetic data, not a physical or owner restore.
  for(bool negative:{false,true}){
   cs::CommonLedger tiny;tiny.accepted_step=1;tiny.begin_s=1;tiny.time_s=2;
   for(const auto&s:states){tiny.conductor_species.emplace(s.id,g::Exact{});tiny.exterior_species.emplace(s.id,g::Exact{});}
   const double x=std::numeric_limits<double>::denorm_min();
   if(negative)tiny.external_net_J.add_product(x,x);else tiny.accepted_energy_J.add_product(x,x);
   tiny.raw_residual_J.add_product(x,x,negative);tiny.accounted_residual_J=tiny.raw_residual_J;tiny.admitted_bound_sum_J.add_product(x,1);tiny.sha256=cs::ledger_hash(tiny);
   const auto bytes=lc::encode(tiny);auto decoded=lc::decode(bytes,tiny.sha256,1,2,states);
   check("signed below-binary64 exact energy is retained",!decoded->accounted_residual_J.zero()&&decoded->accounted_residual_J.value()==0&&decoded->accounted_residual_J==tiny.accounted_residual_J&&lc::encode(*decoded)==bytes);
   refuses("sub-binary64 exact imbalance fails zero bound",[&]{tiny.admitted_bound_sum_J={};tiny.sha256=cs::ledger_hash(tiny);lc::encode(tiny);});
  }
 }catch(const std::exception&e){result["error"]=e.what();check("unexpected refusal",false);}
 result["passed"]=failed==0;result["failures"]=failed;std::cout<<result.dump(2)<<'\n';return failed?1:0;
}
