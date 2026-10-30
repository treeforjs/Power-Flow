// Public values only. No native/Hall/source/RNG authority is constructed.
#define main preserved_original_predicate_test_main
#include "snapshot/tests/test_common_commit_contract.cpp"
#undef main
#include <limits>

int main(){
  const auto accepts=[](const Fixture& x){return c::accepts_closed_step(x.e,x.current,x.proposed);};
  auto change_time=[](Fixture& x,double begin,double dt){
    x.e.before.time_s=begin;x.e.dt_s=dt;x.e.end_s=begin+dt;
    x.current.native=x.e.before;observe(x.current);
    x.proposed.after.native.time_s=x.e.end_s;observe(x.proposed.after);
  };
  auto above=fixture();change_time(above,5e-9,2e-14);
  auto tiny=fixture();change_time(tiny,0,std::numeric_limits<double>::denorm_min());
  auto genesis=fixture();genesis.e.before.authored_full_time_genesis=true;
  genesis.current.native=genesis.e.before;observe(genesis.current);
  auto one=fixture();one.proposed.after.native.accepted_arena_revision=one.e.before.accepted_arena_revision+1;
  observe(one.proposed.after);
  auto maximum=fixture();maximum.e.before.accepted_step=UINT64_MAX-2;
  maximum.e.before.accepted_arena_revision=UINT64_MAX-2;
  maximum.current.native=maximum.e.before;observe(maximum.current);
  maximum.proposed.after.native.accepted_step=UINT64_MAX-1;
  maximum.proposed.after.native.accepted_arena_revision=UINT64_MAX-1;observe(maximum.proposed.after);
  auto exhausted=maximum;exhausted.e.before.accepted_step=UINT64_MAX-1;
  exhausted.current.native=exhausted.e.before;observe(exhausted.current);
  exhausted.proposed.after.native.accepted_step=UINT64_MAX;observe(exhausted.proposed.after);
  std::cout<<std::boolalpha
    <<"{\"status\":\"public_value_observations\",\"above_native_10fs_guard_accepted\":"<<accepts(above)
    <<",\"minimum_subnormal_dt_accepted\":"<<accepts(tiny)
    <<",\"minimum_subnormal_half_is_zero\":"<<(.5*tiny.e.dt_s==0)
    <<",\"nonzero_step_authored_genesis_accepted\":"<<accepts(genesis)
    <<",\"one_revision_per_step_accepted\":"<<accepts(one)
    <<",\"penultimate_counter_boundary_accepted\":"<<accepts(maximum)
    <<",\"exhausting_counter_boundary_accepted\":"<<accepts(exhausted)
    <<",\"native_or_Hall_owner_constructed\":false,\"numerical_operator_executed\":false,\"A02_closed\":false}\n";
  return 0;
}
