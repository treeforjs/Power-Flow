#pragma once
#include "axial_particles.hpp"
#include "surface_path_port.hpp"
namespace axial {
inline foil::SurfacePathReceipt import_surface_paths(const Mesh&g,const foil::CommonInterval&common,
    const foil::CommonInterval&source,V3 patch_origin_m,std::uint64_t canonical_particle_id,
    const std::set<std::uint64_t>&excluded_from_ordinary_push,
    const std::vector<foil::SurfacePathSI>&paths,const Vec&same_stage_E){
    return foil::import_surface_paths_t<Mesh,Path>(g,common,source,patch_origin_m,canonical_particle_id,excluded_from_ordinary_push,paths,same_stage_E);
}
inline foil::SurfacePathReceipt import_global_surface_paths(const Mesh&g,const foil::CommonInterval&common,
    const foil::CommonInterval&source,std::uint64_t canonical_particle_id,
    const std::set<std::uint64_t>&excluded_from_ordinary_push,
    const std::vector<foil::SurfacePathSI>&global_paths,const Vec&same_stage_E){
    return foil::import_surface_paths_t<Mesh,Path>(g,common,source,{0,0,0},canonical_particle_id,excluded_from_ordinary_push,global_paths,same_stage_E);
}
}
