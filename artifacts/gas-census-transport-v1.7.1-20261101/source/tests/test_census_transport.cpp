#include "pf_gas_v1/gas_owner.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
using namespace pf_gas_v1;
unsigned checks=0;
void check(bool x,const char*why){++checks;if(!x)throw std::runtime_error(why);}
template<class F>void refuse(F f){bool yes=false;try{f();}catch(const std::exception&){yes=true;}check(yes,"expected refusal");}
Configuration config(){Configuration c;c.preset="diagnostic-conservative-census";c.transport=Transport::census_conservative;
 c.species={{"water",2.9915e-26,0,true},{"electron",9.109383713868749e-31,-1,false}};
 c.grid={{0,0,0},{3,1,1},{3,1,1},1};c.chemistry_blob_sha256=std::string(64,'a');c.chemistry_registry_sha256=std::string(64,'b');return c;}
CensusRow row(std::uint64_t id,double w,unsigned cell=0,unsigned sp=0){CensusRow r{};r.native_id=id;r.weight=w;r.species=sp;r.position_m[0]=cell+.5;r.position_m[1]=r.position_m[2]=.5;return r;}
void seed(Owner&o,unsigned cell=0,double weight=1){FlowRecord f{};f.kind=pfg::kFlowDebit;f.negative=1;f.grid_cell=cell;f.weight=weight;f.native_id=1000+cell;f.u[0]=8000;f.u[1]=-2000;f.u[2]=3000;f.internal_J=3e-19;
 o.take(0,1,{},std::span(&f,1));}
Exact global(Owner&o,unsigned index){Exact x;auto d=o.global_digits(index,true);x.add_digits(d.data());return x;}
std::array<Exact,5> total(Owner&o){std::array<Exact,5>x;for(unsigned q=0;q<5;++q)x[q]=o.entry_total(q,true);return x;}
void conservation(Owner&o,const std::array<Exact,5>&initial,const TransportReceipt&r){
 for(unsigned q=0;q<5;++q){Exact e=o.entry_total(q,true);auto c=global(o,o.layout().conductor(0,q)),x=global(o,o.layout().exterior(0,q));e.merge(c);e.merge(x);check(e==initial[q],"exact total+resolved exits");
  Exact shares;for(const auto&exit:r.exits)shares.add_digits(exit.digits[q].data());c.merge(x);check(c==global(o,o.layout().wall(0,q)),"resolved equals aggregate");check(c==shares,"per-exit exact shares");}}
bool has_origin(const Owner&o,unsigned cell){for(const auto&e:o.entries(true))if(e.cell==cell)return true;return false;}
void complete(Owner&o){o.take(0,2,{},{});auto p=o.prepare_commit(o.candidate_sha256());check(o.try_publish(p),"publication");}
int main(){try{
 const auto cfg=config();
 // Same weights and deposit that expose the v1.6 preimage, with both exit destinations.
 std::vector<CensusRow>b={row(1,2),row(2,6),row(4,4),row(9,3,0,1)},a={b[0]};a[0].position_m[0]=1.5;
 std::vector<ExitRecord>ex={{2,ExitDestination::conductor},{4,ExitDestination::exterior},{9,ExitDestination::conductor}};
 Owner o(cfg);o.begin_step(0);seed(o);auto initial=total(o);auto old=o.candidate_binding();
 FlowRecord born{};born.kind=pfg::kFlowBirth;born.weight=4;born.native_id=4;o.apply_flows(0,std::span(&born,1));
 auto r=o.transport_census(0,b,a,ex);check(r.drained_cells==1 && r.movers==3 && r.exits.size()==2,"actual census cardinality");check(!o.matches_candidate(old),"mutation invalidates binding");
 check(!has_origin(o,0),"empty source exactly absent");conservation(o,initial,r);
 for(unsigned q=0;q<5;++q){Exact conductor;conductor.add_product(initial[q].value()*.5,1);check(conductor==global(o,o.layout().conductor(0,q)),"conductor ordinary half share");
  Exact expected=initial[q];Exact first;first.add_product(initial[q].value()*(2./12.),1);expected.merge(first,1);expected.merge(conductor,1);
  check(expected==global(o,o.layout().exterior(0,q)),"highest-ID newborn exterior receives exact residue");}
 complete(o);auto encoded=o.encode();const auto accepted=o.accepted_sha256();auto restored=Owner::restore(encoded,cfg,accepted);check(restored->encode()==encoded,"resolved ledger checkpoint");
 auto corrupt=encoded;corrupt[8]=1;refuse([&]{Owner::restore(corrupt,cfg,accepted);});corrupt=encoded;corrupt.back()^=1;refuse([&]{Owner::restore(corrupt,cfg,accepted);});
 // Input permutations preserve candidate, exact exit records and input commitment.
 Owner perm(cfg);perm.begin_step(0);seed(perm);perm.apply_flows(0,std::span(&born,1));
 std::reverse(b.begin(),b.end());std::reverse(a.begin(),a.end());std::reverse(ex.begin(),ex.end());
 auto pr=perm.transport_census(0,b,a,ex);check(pr.inputs_sha256==r.inputs_sha256 && pr.candidate_after_sha256==r.candidate_after_sha256,"canonical permutation commitment");
 for(unsigned i=0;i<r.exits.size();++i)check(pr.exits[i].digits==r.exits[i].digits,"permutation exact exit digits");
 // Non-dyadic unequal weights; highest ID exits but other carrier moves to cell1.
 for(auto weights: {std::array<double,3>{.1,.2,.3},std::array<double,3>{.125,.375,.5}}){
  Owner x(cfg);x.begin_step(0);seed(x);auto t=total(x);std::vector<CensusRow>bb={row(1,weights[0]),row(2,weights[1]),row(3,weights[2])},aa={bb[0],bb[1]};for(auto&v:aa)v.position_m[0]=1.5;
  ExitRecord sink{3,ExitDestination::exterior};auto rr=x.transport_census(0,bb,aa,std::span(&sink,1));check(!has_origin(x,0),"split weights clear emptied source");conservation(x,t,rr);}
 // A tiny real resident is decisive even when RN(before count)==RN(outgoing weight).
 {Owner x(cfg);x.begin_step(0);seed(x,0,.1);auto t=total(x);std::vector<CensusRow>bb={row(1,1),row(2,0x1p-60)},aa=bb;aa[0].position_m[0]=1.5;
  auto rr=x.transport_census(0,bb,aa,{});check(rr.drained_cells==0,"tiny resident prevents empty inference");conservation(x,t,rr);
  bool expected_nonzero=false;for(unsigned q=0;q<5;++q){Exact residual=t[q];residual.add_product(t[q].value(),1,1);expected_nonzero|=!residual.zero();}
  check(expected_nonzero && has_origin(x,0),"tiny resident keeps exact rounded-share residue");}
 // Old residents leave cell1, but an actual incoming survivor retains its correction.
 {Owner x(cfg);x.begin_step(0);FlowRecord f[2]{};for(unsigned i=0;i<2;++i){f[i].kind=pfg::kFlowDebit;f[i].negative=1;f[i].weight=.1;f[i].grid_cell=i;f[i].u[0]=8000;f[i].native_id=100+i;f[i].internal_J=3e-19;}
  x.take(0,1,{},f);auto t=total(x);std::vector<CensusRow>bb={row(1,1,0),row(2,1,1)},aa=bb;aa[0].position_m[0]=1.5;aa[1].position_m[0]=2.5;
  auto rr=x.transport_census(0,bb,aa,{});check(rr.drained_cells==1 && has_origin(x,1),"incoming carrier owns cell1 residue");conservation(x,t,rr);}
 // Strong refusal: candidate, accepted state and read binding survive invalid inputs.
 auto bad=[&](auto change){Owner x(cfg);x.begin_step(0);seed(x);auto bb=std::vector<CensusRow>{row(1,2),row(2,6)},aa=std::vector<CensusRow>{bb[0]};aa[0].position_m[0]=1.5;auto ee=std::vector<ExitRecord>{{2,ExitDestination::conductor}};
  change(bb,aa,ee);const auto c=x.candidate_sha256(),acc=x.accepted_sha256();auto token=x.candidate_binding();refuse([&]{x.transport_census(0,bb,aa,ee);});check(x.candidate_sha256()==c && x.accepted_sha256()==acc && x.matches_candidate(token),"refusal leaves exact state and binding");};
 bad([](auto&bb,auto&,auto&){bb[0].native_id=0;});bad([](auto&bb,auto&,auto&){bb.push_back(bb[0]);});
 bad([](auto&,auto&aa,auto&){aa[0].weight=3;});bad([](auto&,auto&aa,auto&){aa[0].species=1;});
 bad([](auto&,auto&aa,auto&){aa.push_back(row(88,1));});bad([](auto&,auto&aa,auto&){aa.clear();});
 bad([](auto&,auto&,auto&ee){ee[0].native_id=999;});bad([](auto&,auto&,auto&ee){ee.push_back(ee[0]);});
 bad([](auto&,auto&,auto&ee){ee[0].destination=static_cast<ExitDestination>(99);});
 bad([](auto&bb,auto&,auto&){bb[0].position_m[0]=std::numeric_limits<double>::quiet_NaN();});
 bad([](auto&bb,auto&,auto&){bb[0].weight=0;});
 bad([](auto&bb,auto&,auto&){bb[0].proper_velocity_m_s[0]=std::numeric_limits<double>::infinity();});
 {Owner x(cfg);x.begin_step(0);seed(x);const auto c=x.candidate_sha256();refuse([&]{x.transport_census(0,{}, {}, {});});check(x.candidate_sha256()==c,"orphan correction not erased");refuse([&]{x.transport(0,{});});}
 {auto c=cfg;c.max_transport_exit_bytes=sizeof(ExactExitShare);Owner x(c);x.begin_step(0);seed(x);const auto old=x.candidate_sha256();std::vector<CensusRow>bb={row(1,1),row(2,1)};ExitRecord e[2]={{1,ExitDestination::conductor},{2,ExitDestination::exterior}};refuse([&]{x.transport_census(0,bb,{},e);});check(x.candidate_sha256()==old,"payload cap before mutation");}
 {auto c=cfg;c.max_transport_carriers=1;Owner x(c);x.begin_step(0);std::vector<CensusRow>bb={row(1,1),row(2,1)};refuse([&]{x.transport_census(0,bb,bb,{});});}
 for(auto lane:{Execution::device,Execution::device_verified}){auto c=cfg;c.execution=lane;refuse([&]{Owner x(c);});}
 {auto c=cfg;c.transport=Transport::carrier_shares;Owner x(c);x.begin_step(0);refuse([&]{x.transport_census(0,{}, {}, {});});}
 // Abandon successful motion, then restore/retry exactly from accepted boundary.
 auto retry=Owner::restore(encoded,cfg,accepted);auto attempt=retry->begin_attempt(1);std::vector<CensusRow>bb={row(1,2,1)},aa=bb;aa[0].position_m[0]=2.5;
 retry->take(1,1,{},{});const auto rr=retry->transport_census(1,bb,aa,{});check(retry->abandon_attempt(attempt) && retry->encode()==encoded,"rollback exact checkpoint");
 retry->begin_step(1);retry->take(1,1,{},{});const auto rt=retry->transport_census(1,bb,aa,{});check(rt.inputs_sha256==rr.inputs_sha256 && rt.candidate_after_sha256==rr.candidate_after_sha256,"checkpoint restart deterministic");
 std::cout<<"{\"passed\":true,\"checks\":"<<checks<<",\"production_qualified\":false,\"device_qualified\":false}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
