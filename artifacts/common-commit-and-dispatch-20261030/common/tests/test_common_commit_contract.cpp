#include "common_commit_contract.hpp"
#include <functional>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace c = raikou::native_water::common_commit;
std::string hash(char x){return std::string(64,x);}
struct Fixture {c::ExpectedStep e;c::ObservedOwners current;c::ProposedOwners proposed;};
void observe(c::ObservedOwners& o){
  const auto& b=o.native;o.actual_native_prefix=b.accepted_arena_sha256;
  o.actual_native_revision=b.accepted_arena_revision;o.actual_native_step=b.accepted_step;
  o.field=b.field_sha256;o.heat=b.heat_sha256;o.source=b.source_sha256;o.rng=b.rng_sha256;
  o.mesh=b.mesh_sha256;o.material=b.material_sha256;o.frame=b.frame_sha256;
  o.configuration=hash('0');
  o.field_revision=o.heat_revision=b.accepted_step;o.field_time_s=o.heat_time_s=b.time_s;
}
Fixture fixture(){
  Fixture f;auto& b=f.e.before;b.accepted_step=7;b.accepted_arena_revision=19;b.time_s=5e-9;
  b.accepted_arena_sha256=hash('1');b.field_sha256=hash('2');b.heat_sha256=hash('3');
  b.source_sha256=hash('4');b.rng_sha256=hash('5');b.mesh_sha256=hash('6');
  b.material_sha256=hash('7');b.frame_sha256=hash('8');b.common_clock_sha256=hash('9');
  f.e.dt_s=1e-14;f.e.end_s=b.time_s+f.e.dt_s;f.e.common_clock=hash('e');f.e.configuration=hash('0');
  f.current.native=b;observe(f.current);
  auto& a=f.proposed.after.native;a=b;a.accepted_step=8;a.accepted_arena_revision=22;
  a.time_s=f.e.end_s;a.accepted_arena_sha256=hash('a');a.field_sha256=hash('b');
  a.heat_sha256=hash('c');a.rng_sha256=hash('d');a.common_clock_sha256=f.e.common_clock;
  observe(f.proposed.after);f.proposed.motion_midpoint=f.proposed.field_midpoint=hash('f');
  return f;
}
int main(){try{
  int count=0;auto require=[&](bool ok,const char* name){++count;if(!ok)throw std::runtime_error(name);};
  const auto accepts=[](const Fixture& f){return c::accepts_closed_step(f.e,f.current,f.proposed);};
  auto f=fixture();require(accepts(f),"consistent value-only fixture refused");
  require(f.e.end_s-f.e.before.time_s!=f.e.dt_s,"fixture must exercise rounded absolute subtraction");
  auto reject=[&](const char* name,const std::function<void(Fixture&)>& change){auto x=fixture();change(x);require(!accepts(x),name);};
  // Direct audit counterexamples: actual candidate values change while native
  // seals remain fixed; these must fail before the first accepted owner swap.
  reject("changed actual proposed field accepted",[](auto& x){x.proposed.after.field=hash('0');});
  reject("valid alternate heat candidate accepted",[](auto& x){x.proposed.after.heat=hash('0');});
  reject("changed candidate time accepted",[](auto& x){x.proposed.after.field_time_s=std::nextafter(x.e.end_s,INFINITY);});
  reject("changed heat revision accepted",[](auto& x){++x.proposed.after.heat_revision;});
  reject("different work midpoint accepted",[](auto& x){x.proposed.motion_midpoint=hash('0');});
  reject("stale accepted field accepted",[](auto& x){x.current.field=hash('0');});
  reject("stale accepted heat accepted",[](auto& x){x.current.heat=hash('0');});
  reject("stale accepted RNG accepted",[](auto& x){x.current.rng=hash('0');});
  reject("live accepted mesh drift accepted",[](auto& x){x.current.mesh=hash('f');});
  reject("candidate material observation drift accepted",[](auto& x){x.proposed.after.material=hash('f');});
  reject("candidate frame observation drift accepted",[](auto& x){x.proposed.after.frame=hash('f');});
  reject("changed current drive configuration accepted",[](auto& x){x.current.configuration=hash('f');});
  reject("changed candidate forcing configuration accepted",[](auto& x){x.proposed.after.configuration=hash('f');});
  reject("missing expected configuration accepted",[](auto& x){x.e.configuration="";});
  reject("source mutation in closed scope accepted",[](auto& x){x.proposed.after.native.source_sha256=hash('0');observe(x.proposed.after);});
  reject("new geometry in closed scope accepted",[](auto& x){x.proposed.after.native.mesh_sha256=hash('0');observe(x.proposed.after);});
  reject("new material in closed scope accepted",[](auto& x){x.proposed.after.native.material_sha256=hash('0');observe(x.proposed.after);});
  reject("new frame in closed scope accepted",[](auto& x){x.proposed.after.native.frame_sha256=hash('0');observe(x.proposed.after);});
  reject("candidate RNG seal mismatch accepted",[](auto& x){x.proposed.after.rng=hash('0');});
  reject("candidate native prefix mismatch accepted",[](auto& x){x.proposed.after.actual_native_prefix=hash('0');});
  reject("candidate native revision mismatch accepted",[](auto& x){++x.proposed.after.actual_native_revision;});
  reject("candidate native step mismatch accepted",[](auto& x){++x.proposed.after.actual_native_step;});
  reject("wrong common clock accepted",[](auto& x){x.proposed.after.native.common_clock_sha256=hash('0');observe(x.proposed.after);});
  reject("restamped genesis accepted",[](auto& x){x.proposed.after.native.authored_full_time_genesis=true;});
  reject("non-incrementing native revision accepted",[](auto& x){x.proposed.after.native.accepted_arena_revision=x.e.before.accepted_arena_revision;observe(x.proposed.after);});
  reject("skipped common step accepted",[](auto& x){++x.proposed.after.native.accepted_step;observe(x.proposed.after);});
  reject("wrong physical end accepted",[](auto& x){x.e.end_s=std::nextafter(x.e.end_s,INFINITY);x.proposed.after.native.time_s=x.e.end_s;observe(x.proposed.after);});
  reject("zero duration accepted",[](auto& x){x.e.dt_s=0;});
  reject("negative duration accepted",[](auto& x){x.e.dt_s=-1e-14;});
  reject("NaN duration accepted",[](auto& x){x.e.dt_s=NAN;});
  reject("duration above native guard accepted",[](auto& x){x.e.dt_s=2e-14;x.e.end_s=x.e.before.time_s+x.e.dt_s;x.proposed.after.native.time_s=x.e.end_s;observe(x.proposed.after);});
  reject("zero representable half-exposure accepted",[](auto& x){x.e.before.time_s=0;x.current.native=x.e.before;observe(x.current);x.e.dt_s=std::numeric_limits<double>::denorm_min();x.e.end_s=x.e.dt_s;x.proposed.after.native.time_s=x.e.end_s;observe(x.proposed.after);});
  reject("infinite endpoint accepted",[](auto& x){x.e.end_s=INFINITY;});
  reject("negative genesis time accepted",[](auto& x){x.e.before.time_s=-1;});
  reject("unrepresentable endpoint advance accepted",[](auto& x){x.e.dt_s=1e-100;x.e.end_s=x.e.before.time_s;x.proposed.after.native.time_s=x.e.end_s;observe(x.proposed.after);});
  reject("exhausted accepted step accepted",[](auto& x){x.e.before.accepted_step=UINT64_MAX;x.current.native=x.e.before;observe(x.current);x.proposed.after.native.accepted_step=0;observe(x.proposed.after);});
  reject("unrestorable candidate revision accepted",[](auto& x){x.proposed.after.native.accepted_arena_revision=UINT64_MAX;observe(x.proposed.after);});
  reject("missing midpoint identity accepted",[](auto& x){x.proposed.field_midpoint=x.proposed.motion_midpoint="";});
  reject("malformed field digest accepted",[](auto& x){x.proposed.after.native.field_sha256=hash('z');observe(x.proposed.after);});
  reject("short common digest accepted",[](auto& x){x.e.common_clock="abc";x.proposed.after.native.common_clock_sha256="abc";});
  reject("changed all accepted observations without new trial accepted",[](auto& x){x.current.native.field_sha256=hash('0');observe(x.current);});
  require(c::same_boundary(f.e.before,f.current.native),"read-only predicate changed accepted input");
  require(accepts(f),"read-only repeat changed result");
  std::cout<<"{\"status\":\"passed\",\"checks\":"<<count
    <<",\"scope\":\"cross-owner value predicate only\",\"native_or_field_solve\":false,\"full_joint_test\":false}\n";
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
