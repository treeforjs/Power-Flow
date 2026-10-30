// Metadata/type smoke only against the observed full-time draft declarations.
// No missing native implementation is bypassed or claimed to be qualified.
#include "native_parent_subflow.hpp"
#include "complete_indexed_census.hpp"
#include <type_traits>
#include <iostream>
namespace m=raikou::mfem_backend;namespace cg=raikou::collision_geometry;namespace nw=raikou::native_water;
using ArenaView=decltype(cg::complete_indexed_census(std::declval<const m::PicArena&>(),std::declval<const cg::CollisionAtlas&>()));
using ScopeView=decltype(cg::complete_indexed_census(std::declval<const m::CollisionScope&>(),std::declval<const cg::CollisionAtlas&>()));
static_assert(std::is_same_v<decltype(std::declval<ArenaView>().entries),std::vector<m::ArenaRateEntry>>);
static_assert(std::is_same_v<ArenaView,ScopeView>);
int main(){m::GapPicOptions o;o.width_m=o.gap_m=.009;o.depth_m=.025;o.short_z_m=0;o.feed_z_m=.025;o.finite_axial_extent=o.finite_foils=true;o.lower_foil={.002,.007,.0018,.002};o.upper_foil={.002,.007,.007,.0072};nw::Grid old(o,{45,45,5},1000000);auto g=cg::Geometry::from_native_options(o);auto atlas=cg::CollisionAtlas::legacy_uniform(g,cg::FrameBinding(nw::sha("full-time-canonical-helper-compile")),{45,45,5});if(atlas.legacy_uniform_identity()!=old.identity)return 1;std::cout<<"Full-time header + selected collision header + PicArena/CollisionScope helper types PASS; metadata only; no native implementation/event/field/evolution qualification\n";}
