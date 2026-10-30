#include "phase_profile.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

int main(int argc,char** argv) {
 using namespace pic_phase;
 {Scope startup(Phase::CudaOther);}
 assert(steps.empty());
 begin_step(37);interval(1.25e-9,1e-14);
 {
  Scope selection(Phase::DriverSelection);
  {Scope diagnostics(Phase::Diagnostics);}
 }
 {
  Scope advance(Phase::Advance);
  {Scope drift(Phase::DriftFirst);}
  {
   Scope fields(Phase::SolveFieldsFirst);
   assert(solve_phase(Phase::CudaFirst,Phase::CudaSecond,Phase::CudaOther)==Phase::CudaFirst);
   {Scope cuda(solve_phase(Phase::CudaFirst,Phase::CudaSecond,Phase::CudaOther));}
   {Scope residual(solve_phase(Phase::ResidualFirst,Phase::ResidualSecond,Phase::ResidualOther));}
  }
  {Scope force(Phase::ForceLoop);force.stop();force.stop();}
  {Scope drift(Phase::DriftSecond);}
  {
   Scope fields(Phase::SolveFieldsSecond);
   assert(solve_phase(Phase::CudaFirst,Phase::CudaSecond,Phase::CudaOther)==Phase::CudaSecond);
  }
 }
 assert(current_scope==nullptr);end_step();
 begin_step(38);interval(1.25001e-9,1e-14);
 {Scope diagnostics(Phase::Diagnostics);}
 {Scope empty(Phase::AdvanceEmpty);}
 end_step();
 begin_step(39);
 {Scope terminal_diagnostics(Phase::Diagnostics);}
 end_terminal_condition();
 assert(current_step==nullptr && current_scope==nullptr);
 assert(steps.size()==2 && steps[0].step==37 && steps[1].step==38);
 assert(totals[static_cast<std::size_t>(Phase::Diagnostics)].calls==2);
 assert(terminal.phases[static_cast<std::size_t>(Phase::Diagnostics)].calls==1);
 assert(totals[static_cast<std::size_t>(Phase::ForceLoop)].calls==1);
 assert(totals[static_cast<std::size_t>(Phase::CudaFirst)].calls==1);
 assert(totals[static_cast<std::size_t>(Phase::CudaOther)].calls==0);
 assert(totals[static_cast<std::size_t>(Phase::ResidualFirst)].calls==1);
 for(const auto& step:steps) {
  assert(std::abs(exclusive_sum(step.phases)-step.wall_s)<1e-12);
  for(const auto& entry:step.phases) {
   assert(entry.inclusive>=entry.exclusive && entry.exclusive>=0);
  }
 }
 assert(std::abs(exclusive_sum(totals)-total_wall_s)<1e-12);
 if(argc>1)write(argv[1]);
 std::cout<<"Nested additive accounting, scope routing, startup exclusion, terminal exclusion, and counts passed.\n";
}
