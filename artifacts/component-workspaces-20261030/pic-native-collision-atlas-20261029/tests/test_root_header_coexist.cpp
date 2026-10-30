// This uses the ACTUAL immutable root header and its original collision header.
// The atlas's pinned dependency copy stays private to its separate library TU.
#include "native_parent_subflow.hpp"
#include "physical_collision_atlas.hpp"
#include <iostream>
int main(){
 namespace m=raikou::mfem_backend;namespace nw=raikou::native_water;namespace cg=raikou::collision_geometry;
 m::GapPicOptions o;o.width_m=o.gap_m=.009;o.depth_m=.025;o.short_z_m=0;o.feed_z_m=.025;o.finite_axial_extent=o.finite_foils=true;o.lower_foil={.002,.007,.0018,.002};o.upper_foil={.002,.007,.007,.0072};
 nw::Grid old(o,{45,45,5},1000000);auto g=cg::Geometry::from_native_options(o);auto a=cg::CollisionAtlas::legacy_uniform(g,cg::FrameBinding(nw::sha("canonical-only-coexist-test")),{45,45,5});
 if(a.legacy_uniform_identity()!=old.identity||a.locate({{.004,.004,.0125}})!=old.locate({.004,.004,.0125}))return 1;
 std::cout<<"Actual root header + native atlas coexist with one original public collision header; metadata-only PASS\n";
}
