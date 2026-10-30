#include "native_parent_subflow.hpp"
#include <iostream>
#include <fstream>
namespace nw=raikou::native_water;
namespace m=raikou::mfem_backend;
namespace pc=raikou::pic_collisions;
namespace w=raikou::pic_water;
using J=nw::Json;
J checks=J::array(),receipts=J::array();
void check(bool ok,const std::string& name){checks.push_back({{"name",name},{"passed",ok}});if(!ok)throw std::runtime_error(name);}
template<class F>void reject(F f,const std::string& s){bool threw=false;try{f();}catch(const std::exception&){threw=true;}check(threw,s);}
constexpr double mn=2.9915e-26;
m::GapPicOptions options(){
  m::GapPicOptions o;o.nx=o.ny=45;o.width_m=o.gap_m=.009;o.depth_m=.025;
  o.finite_foils=true;o.finite_axial_extent=true;o.short_z_m=0;o.feed_z_m=.025;
  o.lower_foil={.002,.007,.0018,.002};o.upper_foil={.002,.007,.007,.0072};return o;
}
m::PicVector electron_u(double E){double x=E*pc::qe/(w::me*pc::c*pc::c);return {pc::c*std::sqrt(x*(x+2)),0,0};}
std::vector<m::ArenaStateSpec> states(){return {
  {"electron",0,0,-1,w::me,0,0},{"H2O",2,1,0,mn,0,0},
  {"H2O+",2,1,1,mn-w::me+w::threshold_eV*pc::qe/(pc::c*pc::c),0,w::threshold_eV*pc::qe}};}
std::vector<m::ArenaChannel> channels(const nw::Binding& b){return {{w::parent_channel,"electron","H2O",{"electron","electron","H2O+"},b.dataset_sha256,b.kernel_sha256,b.kernel_binding_file_sha256,w::parent_model}};}
struct Case {
  m::GapPic pic{options()};std::unique_ptr<m::PicArena> arena;
  Case(const nw::Binding& b,double ew=1,double nwgt=16,double E=100,bool separate=false){
    m::GapParticle e,n;e.weight=ew;e.mass_kg=w::me;e.position_m={.004,.004,.0125};e.proper_velocity_m_s=electron_u(E);
    n.weight=nwgt;n.mass_kg=mn;n.charge_C=0;n.position_m={separate?.005:.00403,.004,.0125};
    auto ids=pic.add_authored_batch({e,n});pic.solve_fields();
    arena=std::make_unique<m::PicArena>(pic,states(),channels(b),std::vector<m::ArenaPacket>{
      {"e","electron","authored_fixture","",{ew,0},{ids[0]}},
      {"n","H2O","authored_fixture","",{nwgt,0},{ids[1]}}},"native-water-subflow-fixture/v1");
  }
};
nw::Interval interval(m::PicArena& arena,const nw::ParentSubflow& flow,double duration,std::uint64_t sequence=0){
  nw::Interval i;i.sequence=sequence;i.end_s=duration;i.common_clock_sha256=nw::sha("frozen-native-parent-fixture/"+std::to_string(sequence));
  i.phase.state_prefix=arena.state_prefix();i.phase.mesh_prefix=flow.grid_identity();i.phase.clock_prefix=i.common_clock_sha256;return i;
}
long double weight(const m::PicArena& a,const std::string& state){
  long double total=0;for(const auto& [_,p]:a.snapshot()->packets)if(p.state==state)total+=static_cast<long double>(p.count.hi)+p.count.lo;return total;
}
