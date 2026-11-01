#include "pf_gas_v1/gas_owner.hpp"
#include <iostream>
using namespace pf_gas_v1;
int main(){try{Configuration c;c.transport=Transport::carrier_shares;c.species={{"water",2.9915e-26,0,true}};c.grid={{0,0,0},{2,1,1},{2,1,1},1};
 Owner o(c);o.begin_step(0);FlowRecord f{};f.kind=pfg::kFlowDebit;f.negative=1;f.weight=1;f.u[0]=8000;f.internal_J=3e-19;o.take(0,1,{},std::span(&f,1));
 std::vector<CensusRow>b(3),a(2);double weights[]={2,6,4};for(unsigned i=0;i<3;++i){b[i].native_id=i+1;b[i].weight=weights[i];b[i].position_m[0]=b[i].position_m[1]=b[i].position_m[2]=.5;}a[0]=b[0];a[1]=b[2];for(auto&r:a)r.position_m[0]=1.5;
 std::uint64_t sink=2;const auto m=movers(c,b,a,std::span(&sink,1));o.transport(0,m);bool residue=false;for(const auto&e:o.entries(true))residue|=e.cell==0;
 if(!residue) { throw std::runtime_error("preimage did not reproduce"); }
 std::cout<<"REPRODUCED exact correction stranded in empty origin under v1.6\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
