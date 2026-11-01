#include "test_common.hpp"
#include <iostream>
using namespace pfc;
using namespace pfc_test;
int main(int argc,char**argv) {
 require(argc==3,"base blob/SHA");const auto b=load_blob(argv[1],argv[2]);
 std::uint64_t kinds[7]{};Sha256 total;
 for(unsigned case_id=0;case_id<3;++case_id) {
  auto P=make_population(b,10000,31,1000+case_id);
  StepConfig cfg{};cfg.seed=1234567+case_id;cfg.step=42+case_id;cfg.dt_s=case_id==0?4e-11:case_id==1?5e-12:2e-10;cfg.trunc_sigmas=7;cfg.purpose=PFC_PURPOSE_MCC;cfg.max_candidates=4096;cfg.max_events=32;
  StepResult result;cpu_step(b,cfg,P.carriers,P.cells,result);
  total.text(step_digest(result.outs,result.events,cfg.max_events,P.carriers));total.text(ledger_digest(ledger_of(b,result,cfg.max_events)));
  for(std::size_t i=0;i<P.carriers.size();++i)for(unsigned j=0;j<result.outs[i].n_events;++j) {const auto& ev=result.events[i*cfg.max_events+j];const auto kind=b.view.procs[ev.process].kind;require(kind>0 && kind<7,"unexpected base kind");++kinds[kind];}
 }
 for(unsigned k=1;k<=6;++k)require(kinds[k]>0,"legacy process kind unexercised");
 std::cout<<"{\"record_ledger_digest\":\""<<total.hex()<<"\",\"events_by_kind\":[";
 for(unsigned k=1;k<=6;++k){if(k>1)std::cout<<',';std::cout<<kinds[k];}std::cout<<"],\"passed\":true}\n";
}
