#pragma once
#include "foil_particles.hpp"
namespace foil {
struct CommonInterval {std::uint64_t accepted_revision,mesh_revision,material_revision;double begin_s,end_s,field_stage_s,dt_s;
    std::string field_hash,arena_hash,mesh_hash,material_hash;};
struct SurfacePathSI {V3 from_m,to_m;double start_time_s,end_time_s,charge_C,field_work_J;};
struct SurfacePathReceipt {Vec current_load_C_m_s,endpoint_charge_impulse_C;double field_work_J=0,claimed_field_work_J=0,work_defect_J=0,continuity_C=0;int zero_duration_bridges=0;};
// Source paths are already weighted. from/to are surface-relative GLOBAL vector
// components, as in the sealed micro-layer API. Add the physical patch origin;
// never rotate B or multiply charge/work by packet weight or depth a second time.
// This operator returns a staged current contribution. It does not push the
// canonical active electron or commit the source's thermal/image/reference ports.
template<class CompatibleMesh,class CompatiblePath>
inline SurfacePathReceipt import_surface_paths_t(const CompatibleMesh&g,const CommonInterval&common,
    const CommonInterval&source,V3 patch_origin_m,std::uint64_t canonical_particle_id,
    const std::set<std::uint64_t>&excluded_from_ordinary_push,
    const std::vector<SurfacePathSI>&paths,const Vec&same_stage_E){
    hall::require(common.accepted_revision==source.accepted_revision&&common.mesh_revision==source.mesh_revision&&common.material_revision==source.material_revision&&common.begin_s==source.begin_s&&common.end_s==source.end_s&&common.field_stage_s==source.field_stage_s&&common.dt_s==source.dt_s,"surface current stage binding mismatch");
    hall::require(!common.field_hash.empty()&&!common.arena_hash.empty()&&!common.mesh_hash.empty()&&!common.material_hash.empty()&&
        common.field_hash==source.field_hash&&common.arena_hash==source.arena_hash&&common.mesh_hash==source.mesh_hash&&common.material_hash==source.material_hash,
        "surface field/arena/mesh/material identity mismatch");
    const double dt=common.dt_s;hall::require(std::isfinite(dt)&&dt>0&&dt<=1e-14&&common.end_s==common.begin_s+dt&&common.field_stage_s>=common.begin_s&&common.field_stage_s<=common.end_s,"invalid common surface interval");
    hall::require(same_stage_E.size()==std::size_t(g.ne)&&finite(same_stage_E),"surface total-field binding mismatch");
    hall::require(canonical_particle_id>0&&excluded_from_ordinary_push.contains(canonical_particle_id),"active micro-layer particle must have exactly one motion owner");
    SurfacePathReceipt r{Vec(g.ne),Vec(g.nn)};double gross_charge=0,work_scale=0;
    const SurfacePathSI*previous=nullptr;
    for(const auto&p:paths){
        if(previous)hall::require(p.from_m==previous->to_m&&p.start_time_s>=previous->end_time_s&&p.charge_C==previous->charge_C,"micro-layer path sequence or canonical charge changed");previous=&p;
        hall::require(std::isfinite(p.start_time_s)&&std::isfinite(p.end_time_s)&&p.start_time_s>=common.begin_s&&p.end_time_s<=common.end_s&&p.end_time_s>=p.start_time_s&&std::isfinite(p.charge_C)&&std::isfinite(p.field_work_J),"surface path time/charge/work outside admission");
        V3 from=patch_origin_m,to=patch_origin_m,D{};for(int d=0;d<3;++d){hall::require(std::isfinite(p.from_m[d])&&std::isfinite(p.to_m[d])&&std::isfinite(patch_origin_m[d]),"nonfinite surface path");from[d]+=p.from_m[d];to[d]+=p.to_m[d];D[d]=to[d]-from[d];}
        CompatiblePath path(g,from,D);path.deposit(r.current_load_C_m_s,p.charge_C,D,dt);
        point_deposit(g,r.endpoint_charge_impulse_C,from,-p.charge_C);point_deposit(g,r.endpoint_charge_impulse_C,to,p.charge_C);
        double W=p.charge_C*hall::dot(path.gather_E(same_stage_E),D);r.field_work_J+=W;r.claimed_field_work_J+=p.field_work_J;work_scale+=std::abs(W)+std::abs(p.field_work_J);gross_charge+=std::abs(p.charge_C);if(p.start_time_s==p.end_time_s)++r.zero_duration_bridges;
    }
    auto div=g.gradient_transpose(r.current_load_C_m_s);for(int n=0;n<g.nn;++n)r.continuity_C=std::max(r.continuity_C,std::abs(r.endpoint_charge_impulse_C[n]-dt*div[n]));
    r.work_defect_J=dt*dot(same_stage_E,r.current_load_C_m_s)-r.claimed_field_work_J;
    hall::require(r.continuity_C<=1e-12*std::max(gross_charge,1e-30),"surface path charge continuity");
    hall::require(std::abs(r.work_defect_J)<=1e-10*work_scale,"micro-layer frozen field work differs from compatible field work; re-evaluate coupled stage, do not correct energy");return r;
}
inline SurfacePathReceipt import_surface_paths(const Mesh&g,const CommonInterval&common,
    const CommonInterval&source,V3 patch_origin_m,std::uint64_t canonical_particle_id,
    const std::set<std::uint64_t>&excluded_from_ordinary_push,
    const std::vector<SurfacePathSI>&paths,const Vec&same_stage_E){
    return import_surface_paths_t<Mesh,Path>(g,common,source,patch_origin_m,canonical_particle_id,excluded_from_ordinary_push,paths,same_stage_E);
}
// Heat/source v2 BoundLayerTrial.global_paths are ALREADY translated.
inline SurfacePathReceipt import_global_surface_paths(const Mesh&g,const CommonInterval&common,
    const CommonInterval&source,std::uint64_t canonical_particle_id,
    const std::set<std::uint64_t>&excluded_from_ordinary_push,
    const std::vector<SurfacePathSI>&global_paths,const Vec&same_stage_E){
    return import_surface_paths_t<Mesh,Path>(g,common,source,{0,0,0},canonical_particle_id,excluded_from_ordinary_push,global_paths,same_stage_E);
}
}
