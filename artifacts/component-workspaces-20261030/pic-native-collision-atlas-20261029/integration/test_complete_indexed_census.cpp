#include "native_parent_subflow.hpp"
#include "complete_indexed_census.hpp"
#include <iostream>
#include <type_traits>
namespace cg=raikou::collision_geometry;namespace m=raikou::mfem_backend;namespace nw=raikou::native_water;
using NativeView=decltype(cg::complete_indexed_census(std::declval<const m::PicArena&>(),std::declval<const cg::CollisionAtlas&>()));
static_assert(std::is_same_v<decltype(std::declval<NativeView>().entries),std::vector<m::ArenaRateEntry>>);
struct FakeSnapshot {struct Packet{std::vector<std::uint64_t>native_ids;};std::map<std::string,Packet>packets;};
struct FakeArena {
 mutable std::shared_ptr<FakeSnapshot>accepted=std::make_shared<FakeSnapshot>();
 std::vector<m::ArenaRateEntry>rows;mutable unsigned reads=0;bool change_prefix=false,change_snapshot=false;mutable bool changed=false;
 auto snapshot()const{return accepted;}
 std::string state_prefix()const{return nw::sha(changed?"changed":"accepted-fixture");}
 auto rate_population(std::size_t cap)const{raikou::pic_collisions::require(rows.size()<=cap,"fixture cap");++reads;if(change_prefix)changed=true;if(change_snapshot)accepted=std::make_shared<FakeSnapshot>(*accepted);return rows;}
};
void need(bool ok,const char*s){if(!ok)throw std::runtime_error(s);}
template<class F>void refuses(F f){bool refused=false;try{f();}catch(const std::exception&){refused=true;}need(refused,"helper refusal missing");}
FakeArena fixture(){FakeArena f;for(auto pair:{std::pair<std::uint64_t,const char*>{3,"ion"},{1,"electron"},{2,"H2O"}}){m::ArenaRateEntry e;e.packet_id="packet/"+std::to_string(pair.first);e.state=pair.second;e.native.id=pair.first;e.native.position_m={.004,.003+pair.first*.001,.0125};e.native.weight=pair.first*.25;e.native.mass_kg=pair.first==1?9.1093837139e-31:2.9915e-26;e.native.charge_C=pair.first==1?-1.602176634e-19:pair.first==3?1.602176634e-19:0;f.rows.push_back(e);f.accepted->packets[e.packet_id].native_ids.push_back(e.native.id);}return f;}
int main(){try{m::GapPicOptions o;o.width_m=o.gap_m=.009;o.depth_m=.025;o.short_z_m=0;o.feed_z_m=.025;o.finite_axial_extent=o.finite_foils=true;o.lower_foil={.002,.007,.0018,.002};o.upper_foil={.002,.007,.007,.0072};auto g=cg::Geometry::from_native_options(o);auto atlas=cg::CollisionAtlas::legacy_uniform(g,cg::FrameBinding(nw::sha("canonical-only-helper-test")),{45,45,5});auto f=fixture();auto view=cg::complete_indexed_census(f,atlas);need(f.reads==1&&view.entries.size()==3&&view.occupied.entry_cell_ids.size()==3,"one complete native read");for(std::size_t i=0;i<f.rows.size();++i)need(view.entries[i].native.id==f.rows[i].native.id&&view.entries[i].state==f.rows[i].state&&view.occupied.entry_cell_ids[i]==atlas.locate({f.rows[i].native.position_m}),"helper filters/reorders/transforms carrier");
 auto missing=fixture();missing.rows.pop_back();refuses([&]{cg::complete_indexed_census(missing,atlas);});auto prefix=fixture();prefix.change_prefix=true;refuses([&]{cg::complete_indexed_census(prefix,atlas);});auto snap=fixture();snap.change_snapshot=true;refuses([&]{cg::complete_indexed_census(snap,atlas);});auto cap=fixture();cg::CensusLimits limits;limits.carriers=1;refuses([&]{cg::complete_indexed_census(cap,atlas,limits);});need(cap.reads==0,"canonical cap checked too late");
 std::cout<<"Complete mixed-state census helper PASS: native API instantiated; one read; original order/coordinates; missing carrier, changed prefix/snapshot and pre-read caps refused; no chemistry/RNG/PIC/GPU\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
