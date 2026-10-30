#pragma once
#include "collision_scheduler.hpp"
#include <string_view>

namespace raikou::collision_geometry {
namespace pc=pic_collisions;
using Vec=std::array<double,3>;
struct CanonicalPointSI {Vec value;};
struct NativeCenteredXYPointSI {Vec value;};
struct FoilXY {double left_m,right_m,bottom_m,top_m;};
struct Geometry {
 Vec low_m{},high_m{};
 double physical_depth_m=0;
 std::vector<FoilXY> foils;
 template<class NativeOptions>static Geometry from_native_options(const NativeOptions&o){
  pc::require(o.finite_axial_extent,"collision atlas requires physical finite axial bounds");
  Geometry g{{0,0,o.short_z_m},{o.width_m,o.gap_m,o.feed_z_m},o.depth_m,{}};
  if(o.finite_foils)for(const auto&f:{o.lower_foil,o.upper_foil})g.foils.push_back({f.left_m,f.right_m,f.bottom_m,f.top_m});
  return g;
 }
};
// Arena positions already use the canonical frame: atlas calls do NOT shift them.
// Raw native-field sampling has a distinct type, preventing a second transform.
class FrameBinding {
 std::string canonical_,native_,identity_;
 std::array<double,2> origin_{};
public:
 explicit FrameBinding(std::string canonical_frame_sha256); // Canonical TOTAL-E/B owner: no native-field map.
 FrameBinding(const Geometry&,std::string canonical_frame_sha256,std::string native_centered_xy_frame_sha256);
 NativeCenteredXYPointSI native_field_point(CanonicalPointSI)const;
 const std::string& identity()const{return identity_;}
 const std::string& canonical_identity()const{return canonical_;}
 const std::array<double,2>& native_xy_origin_m()const{return origin_;}
 bool has_native_centered_xy_map()const{return !native_.empty();}
};
struct AtlasLimits {std::size_t geometric_bins=1000000,xy_nodes=200000,xy_leaves=100000;};
struct CensusLimits {std::size_t carriers=100000,occupied_cells=100000,string_bytes=64*1024*1024;};
struct GradedOptions {
 double normal_target_m=1e-6,distance_growth=.3;
 double broad_tangent_max_m=125e-6,side_tangent_max_m=25e-6,bulk_max_m=500e-6;
 std::uint32_t axial_bins=50,near_axial_bins=0;
 double axial_refinement_band_m=0;
};
struct CanonicalCensusBinding {
 std::string accepted_state_prefix;
 // Derived independently from the accepted canonical snapshot, before filtering.
 std::vector<std::uint64_t> expected_native_ids;
};
struct OccupiedCellView {
 std::vector<pc::Cell> cells;                 // Ascending IDs; no empty cells.
 std::vector<std::uint64_t> entry_cell_ids;   // Original complete census order.
 std::string accepted_state_prefix,geometry_identity,census_identity;
};
class CollisionAtlas {
 struct Box {double x0,x1,y0,y1;};
 struct Node {std::uint32_t left=0,right=0,leaf=0;int axis=-1;double cut=0;};
 struct Root {Box box;std::uint32_t node;};
 struct Leaf {Box box;std::uint32_t z_family;std::uint64_t offset;};
 struct CarrierView {std::uint64_t id;std::string_view packet,state;const Vec*position;const Vec*proper_velocity;double weight,mass,charge;};
 Geometry geometry_;
 FrameBinding frame_;
 AtlasLimits limits_;
 bool uniform_=false;
 std::array<std::uint32_t,3> counts_{};
 std::vector<std::array<std::uint64_t,4>> solid_indices_;
 double uniform_volume_=0;
 std::uint64_t geometric_bins_=0,vacuum_bins_=0;
 std::vector<Node> nodes_;std::vector<Root> roots_;std::vector<Leaf> leaves_;
 std::vector<std::vector<double>> z_families_;
 std::string identity_,legacy_identity_;
 CollisionAtlas(Geometry,FrameBinding,AtlasLimits);
 std::uint32_t split(Box,const GradedOptions&,unsigned);
 OccupiedCellView occupied_impl(const std::vector<CarrierView>&,const CanonicalCensusBinding&,CensusLimits)const;
 void finish_identity(const GradedOptions*);
 void require_point(CanonicalPointSI)const;
public:
 static CollisionAtlas legacy_uniform(Geometry,FrameBinding,std::array<std::uint32_t,3>,AtlasLimits={});
 static CollisionAtlas graded(Geometry,FrameBinding,GradedOptions,AtlasLimits={});
 std::uint64_t locate(CanonicalPointSI)const;
 double physical_volume(std::uint64_t cell_id)const;
 const std::string& identity()const{return identity_;}
 const std::string& legacy_uniform_identity()const{return legacy_identity_;}
 std::uint64_t geometric_bin_count()const{return geometric_bins_;}
 std::uint64_t vacuum_bin_count()const{return vacuum_bins_;}
 std::size_t transverse_leaf_count()const{return leaves_.size();}
 const Geometry& geometry()const{return geometry_;}
 template<class NativeEntry>OccupiedCellView occupied_cells(const std::vector<NativeEntry>&entries,
    const CanonicalCensusBinding&binding,CensusLimits limits={})const{
  pc::require(limits.carriers<=1000000&&limits.occupied_cells<=1000000&&limits.string_bytes<=256*1024*1024&&entries.size()<=limits.carriers&&binding.expected_native_ids.size()<=limits.carriers,"bounded complete collision census required");
  std::size_t bytes=0;
  for(const auto&e:entries)for(const auto*text:{&e.packet_id,&e.state}){
   pc::require(text->size()<=limits.string_bytes-bytes,"collision census string byte cap");bytes+=text->size();
  }
  std::vector<CarrierView>views;views.reserve(entries.size());
  for(const auto&e:entries){const auto&n=e.native;views.push_back({n.id,e.packet_id,e.state,&n.position_m,&n.proper_velocity_m_s,n.weight,n.mass_kg,n.charge_C});}
  return occupied_impl(views,binding,limits);
 }
};
} // namespace raikou::collision_geometry
