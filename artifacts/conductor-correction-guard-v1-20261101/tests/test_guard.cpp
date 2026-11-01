#include "conductor_correction_guard.hpp"
#include <iostream>
#include <limits>
using namespace pf_gas_v1;
namespace guard=powerflow::conductor_correction_guard_v1;
unsigned checks=0;
void check(bool b,const char*s){++checks;if(!b)throw std::runtime_error(s);}
template<class F>void refuses(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}check(caught,"expected refusal");}
Configuration config(){Configuration c;c.preset="diagnostic-per-exit-conductor-refusal";c.transport=Transport::census_conservative;
 c.species={{"water",2.9915e-26,0,true},{"electron",9.109383713868749e-31,-1,false}};
 c.grid={{0,0,0},{2,1,1},{2,1,1},1};c.chemistry_blob_sha256=std::string(64,'a');c.chemistry_registry_sha256=std::string(64,'b');return c;}
CensusRow row(std::uint64_t id,unsigned cell){CensusRow r{};r.native_id=id;r.weight=1;r.species=0;r.position_m[0]=cell+.5;r.position_m[1]=r.position_m[2]=.5;return r;}
void seed_opposite(Owner&o){FlowRecord f[2]{};for(unsigned i=0;i<2;++i){f[i].kind=pfg::kFlowDebit;f[i].negative=1-i;f[i].grid_cell=i;f[i].weight=1;f[i].native_id=i+1;f[i].u[0]=8000;f[i].u[1]=-2000;f[i].u[2]=3000;}o.take(0,1,{},f);}
bool aggregate_unchanged(const Owner&o){const auto L=o.layout();for(unsigned sp=0;sp<o.configuration().species.size();++sp)for(unsigned q=0;q<5;++q)if(o.global_digits(L.conductor(sp,q),true)!=o.global_digits(L.conductor(sp,q),false))return false;return true;}
int main(){try{
 const auto cfg=config();const CensusRow before[]={row(1,0),row(2,1)};const ExitRecord conductor[]={{1,ExitDestination::conductor},{2,ExitDestination::conductor}};
 Owner o(cfg);const auto pristine=o.encode();const auto accepted=o.accepted_sha256();auto attempt=o.begin_attempt(0);seed_opposite(o);
 const auto r=o.transport_census(0,before,{},conductor);
 check(r.exits.size()==2,"actual per-exit receipt");check(aggregate_unchanged(o),"preimage aggregate conductor guard must pass cancellation");
 for(unsigned q=0;q<5;++q){Exact net;for(const auto&e:r.exits)net.add_digits(e.digits[q].data());check(net.zero(),"actual opposing exit shares sum exactly zero");}
 check(!std::all_of(r.exits[0].digits[0].begin(),r.exits[0].digits[0].end(),[](auto x){return x==0;}),"actual first exit energy nonzero");
 const auto candidate=o.candidate_sha256();const auto binding=o.candidate_binding();
 refuses([&]{guard::require_no_unrouted_conductor_correction(r,0,cfg);});
 check(o.candidate_sha256()==candidate&&o.accepted_sha256()==accepted&&o.matches_candidate(binding),"predicate has no mutation");
 check(o.abandon_attempt(attempt)&&o.encode()==pristine,"refused common attempt rolls back exactly");
 // Restored accepted boundary reproduces actual cancellation receipt and refusal.
 auto restored=Owner::restore(pristine,cfg,accepted);restored->begin_step(0);seed_opposite(*restored);
 const auto replay=restored->transport_census(0,before,{},conductor);
 check(replay.inputs_sha256==r.inputs_sha256&&replay.candidate_after_sha256==r.candidate_after_sha256,"actual checkpoint/replay");
 for(unsigned i=0;i<2;++i)check(replay.exits[i].digits==r.exits[i].digits,"replay exact exits");
 refuses([&]{guard::require_no_unrouted_conductor_correction(replay,0,cfg);});
 // Actual conductor sinks with no correction are admitted, including no exits.
 Owner zero(cfg);zero.begin_step(0);zero.take(0,1,{},{});const auto z=zero.transport_census(0,before,{},conductor);
 auto proof=guard::require_no_unrouted_conductor_correction(z,0,cfg);check(proof.conductor_exits==2&&proof.exterior_exits==0,"zero conductor allowed");
 Owner empty(cfg);empty.begin_step(0);empty.take(0,1,{},{});const auto none=empty.transport_census(0,{}, {},{});
 proof=guard::require_no_unrouted_conductor_correction(none,0,cfg);check(proof.conductor_exits==0&&proof.exterior_exits==0,"empty authentic receipt allowed");
 // Actual nonzero exterior corrections retain their existing exterior ownership.
 Owner ext(cfg);ext.begin_step(0);seed_opposite(ext);const ExitRecord exterior[]={{1,ExitDestination::exterior},{2,ExitDestination::exterior}};
 const auto e=ext.transport_census(0,before,{},exterior);proof=guard::require_no_unrouted_conductor_correction(e,0,cfg);check(proof.exterior_exits==2&&proof.conductor_exits==0,"nonzero exterior allowed");
 // Operand fault controls are copies, not claims of authentic native receipts.
 for(unsigned q=0;q<5;++q){auto bad=z;bad.exits[0].digits[q][0]=1;Exact tiny;tiny.add_digits(bad.exits[0].digits[q].data());check(tiny.value()==0&&!tiny.zero(),"below-binary64 exact correction preimage");refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});
  Exact negative;negative.merge(tiny,1);bad.exits[0].digits[q]=negative.normalised();refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});}
 auto bad=z;bad.exits[0].destination=ExitDestination(99);refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});
 bad=z;bad.exits[0].species=99;refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});
 bad=z;bad.exits[0].species=1;refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});
 bad=z;bad.exits[1].native_id=bad.exits[0].native_id;refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});
 bad=z;bad.exits[0].native_id=0;refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});
 bad=z;std::reverse(bad.exits.begin(),bad.exits.end());refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});
 for(unsigned k=0;k<3;++k){bad=z;(k==0?bad.inputs_sha256:k==1?bad.candidate_before_sha256:bad.candidate_after_sha256).clear();refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});}
 bad=z;bad.inputs_sha256=std::string(64,'z');refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});
 refuses([&]{guard::require_no_unrouted_conductor_correction(z,1,cfg);});
 for(auto digit:{std::int64_t(-1),std::int64_t(0x100000000ull)}){bad=e;bad.exits[0].digits[0][0]=digit;refuses([&]{guard::require_no_unrouted_conductor_correction(bad,0,cfg);});}
 auto small=cfg;small.max_transport_exit_bytes=sizeof(ExactExitShare);refuses([&]{guard::require_no_unrouted_conductor_correction(z,0,small);});
 small=cfg;small.max_transport_carriers=1;refuses([&]{guard::require_no_unrouted_conductor_correction(z,0,small);});
 small=cfg;small.transport=Transport::carrier_shares;refuses([&]{guard::require_no_unrouted_conductor_correction(z,0,small);});
 for(auto mode:{Execution::device,Execution::device_verified}){small=cfg;small.execution=mode;refuses([&]{guard::require_no_unrouted_conductor_correction(z,0,small);});}
 std::cout<<"{\"passed\":true,\"checks\":"<<checks<<",\"actual_aggregate_preimage_passed\":true,\"actual_per_exit_refused\":true,\"actual_zero_and_exterior_allowed\":true,\"native_publication_authority\":false,\"production_qualified\":false}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
