#include "pf_gas_v1/gas_owner.hpp"
#include "allocation_probe.hpp"
#include <iostream>
using namespace pf_gas_v1;
void check(bool x,const char*s){if(!x)throw std::runtime_error(s);}
Configuration config(){Configuration c;c.transport=Transport::census_conservative;c.species={{"gas",2.9915e-26,0,true}};c.grid={{0,0,0},{2,1,1},{2,1,1},1};return c;}
void initialize(Owner&o){o.begin_step(0);FlowRecord f{};f.kind=pfg::kFlowDebit;f.negative=1;f.weight=.1;f.u[0]=8000;f.internal_J=3e-19;o.take(0,1,{},std::span(&f,1));}
int main(){try{const auto c=config();std::vector<CensusRow>b(3),a(1);const double w[]={2,6,4};for(unsigned i=0;i<3;++i){b[i].native_id=i+1;b[i].weight=w[i];b[i].position_m[0]=b[i].position_m[1]=b[i].position_m[2]=.5;}a[0]=b[0];a[0].position_m[0]=1.5;
 ExitRecord e[2]={{2,ExitDestination::conductor},{3,ExitDestination::exterior}};
 Owner normal(c);initialize(normal);allocation_probe::allocations=0;auto r=normal.transport_census(0,b,a,e);const auto count=allocation_probe::allocations;check(count>0&&count<4096,"bounded fault sweep");
 unsigned caught=0;
 for(std::size_t k=0;k<count;++k){Owner o(c);initialize(o);const auto before=o.candidate_sha256(),accepted=o.accepted_sha256();auto binding=o.candidate_binding();bool refused=false;
  allocation_probe::fail_after=long(k);try{(void)o.transport_census(0,b,a,e);}catch(const std::bad_alloc&){refused=true;}catch(...){allocation_probe::fail_after=-1;throw;}allocation_probe::fail_after=-1;
  check(refused,"expected allocation refusal");++caught;check(before==o.candidate_sha256()&&accepted==o.accepted_sha256()&&o.matches_candidate(binding),"allocation changed private candidate/binding");
  auto retry=o.transport_census(0,b,a,e);check(retry.candidate_after_sha256==r.candidate_after_sha256&&retry.inputs_sha256==r.inputs_sha256,"allocation retry changed exact result");
 }
 std::cout<<"{\"passed\":true,\"allocation_refusals\":"<<caught<<",\"candidate_and_binding_unchanged\":true,\"retry_identical\":true}\n";
 }catch(const std::exception&e){allocation_probe::fail_after=-1;std::cerr<<e.what()<<'\n';return 1;}}
