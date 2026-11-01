#include "pf_gas_v1/gas_owner.hpp"
#include <iostream>
#include <stdexcept>
using namespace pf_gas_v1;
template<class F> bool allowed(F f) {try{f();return true;}catch(const Refusal&){return false;}}
Configuration cfg(){Configuration c;c.species={{"neutral",1,0,true}};c.grid.low={0,0,0};c.grid.high={1,1,1};c.grid.n={1,1,1};c.grid.volume=1;return c;}
int main(int argc,char**argv){
 if(argc!=2)return 2;bool old=std::string(argv[1])=="baseline";
 double v[5]={-2,0,0,0,0},d[3]={0,0,0},s=1;int clamp=0;
 bool variance=pfg::feedback_slot(1,1,1,v,d,&s,&clamp)>=0;
 Owner a(cfg());a.begin_step(1);bool early=allowed([&]{a.commit(a.candidate_sha256());});a.abandon();
 Owner b(cfg());b.begin_step(1);b.take(1,1,{},{});bool duplicate=allowed([&]{b.take(1,1,{},{});});b.abandon();
 Owner c(cfg());c.begin_step(1);c.take(1,1,{},{});c.take(1,2,{},{});c.commit(c.candidate_sha256());
 bool replay=allowed([&]{c.begin_step(1);});
 std::cout<<"negative_variance_admitted="<<variance<<" early_commit_admitted="<<early<<" duplicate_stage_admitted="<<duplicate<<" repeated_step_admitted="<<replay<<"\n";
 return variance==old && early==old && duplicate==old && replay==old?0:1;
}
