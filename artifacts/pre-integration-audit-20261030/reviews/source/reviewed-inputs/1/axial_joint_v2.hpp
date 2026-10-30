#pragma once
#include "axial_thermal_v2.hpp"
#include "motion_receipt.hpp"
#include <iomanip>
#include <openssl/sha.h>

namespace hall_driver {
inline std::string text_sha256(const std::string&bytes){unsigned char d[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size(),d);std::ostringstream out;
    for(auto b:d)out<<std::hex<<std::setfill('0')<<std::setw(2)<<unsigned(b);return out.str();}
inline std::string material_binding_sha256(){return text_sha256(Thermal::material_binding().dump());}
inline std::string field_state_sha256(const axial::State&s){return text_sha256(Json{
    {"E",s.E},{"B",s.B},{"E_low",s.E_low},{"B_low",s.B_low},{"time_s",s.time_s},{"revision",s.revision}}.dump());}
inline Json thermal_state_record(const electrode_heat_v2::State&s){Json inventory=Json::array();
    for(const auto&count:s.inventory)inventory.push_back({{"hi",count.hi},{"lo",count.lo}});
    return {{"T",s.T},{"time_s",s.time},{"revision",s.revision},{"cumulative",s.cumulative},{"watermarks",s.watermarks},{"inventory",inventory}};}
inline std::string thermal_state_sha256(const Thermal&t){Json state=Json::array();for(const auto&s:t.states())state.push_back(thermal_state_record(s));
    return text_sha256(Json{{"material",Thermal::material_binding()},{"slices",state}}.dump());}
inline std::string thermal_trial_state_sha256(const ThermalTrial&t){Json state=Json::array();for(const auto&body:t.bodies)state.push_back(thermal_state_record(body.state()));
    return text_sha256(Json{{"material",Thermal::material_binding()},{"slices",state}}.dump());}
struct AxialJointTrial {
    std::string base_field_sha256,base_thermal_sha256,base_census_sha256;
    axial::State field;axial::Receipt field_receipt;
    std::optional<MotionReceipt> motion;ThermalTrial heat;
    double nonlinear_relative=0,particle_work_defect_J=0,combined_balance_J=0;
    unsigned iterations=0;
};
// Pure common preparation over the outer arena's PRIVATE, post-first-exposure
// census. All source/boundary population changes are explicitly unsupported by
// this first port. The root coordinator owns arena/RNG/source publication.
inline AxialJointTrial prepare_axial_joint(const Mesh&g,const axial::Fields&fields,
    const axial::State&accepted,const native::GapPicOptions&native_geometry,
    std::span<const native::ArenaRateEntry>private_census,const MotionBinding&common,
    Thermal&thermal,const CurrentTrace&trace,std::span<const std::uint64_t>active_layers={}){
    hall::require(accepted.time_s==thermal.time_s()&&accepted.revision==thermal.revision(),"v2 joint field/heat clock mismatch");
    hall::require(common.phase==ParticlePhase::full_time_position_and_proper_velocity&&
        common.accepted_step==accepted.revision&&common.begin_s==accepted.time_s&&common.end_s==common.begin_s+common.dt_s&&
        std::isfinite(common.dt_s)&&common.dt_s>0&&common.dt_s<=1e-14,"v2 joint full-time common binding");
    hall::require(common.material_sha256==material_binding_sha256()&&
        common.mesh_sha256==motion_mesh_sha256(g,native_geometry),"v2 joint explicit material/mesh binding");
    AxialJointTrial out;out.base_field_sha256=field_state_sha256(accepted);out.base_thermal_sha256=thermal_state_sha256(thermal);
    out.base_census_sha256=canonical_census_sha256(private_census);const auto sigma=thermal.conductivity();
    const double current=trace.current(common.begin_s+common.dt_s/2);axial::State guess=accepted;bool converged=false;
    auto midpoint=[&](const axial::State&s){std::pair<Vec,Vec> result{Vec(g.ne),Vec(g.nb)};
        for(int i=0;i<g.ne;++i)result.first[i]=.5*(accepted.E[i]+s.E[i]);
        for(int i=0;i<g.nb;++i)result.second[i]=.5*(accepted.B[i]+s.B[i]);return result;};
    for(unsigned iteration=0;iteration<16;++iteration){auto [E,B]=midpoint(guess);auto binding=common;
        binding.trial_field_sha256=total_field_sha256(E,B);
        auto motion=prepare_canonical_motion(g,native_geometry,private_census,binding,active_layers,E,B);
        auto field=fields.prepare(accepted,motion.current_A_m(),sigma,common.dt_s,current);
        auto [nextE,nextB]=midpoint(field.first);double error=0,scale=1;
        for(int i=0;i<g.ne;++i){error=std::max(error,std::abs(nextE[i]-E[i]));scale=std::max(scale,std::abs(nextE[i]));}
        for(int i=0;i<g.nb;++i){error=std::max(error,pic_c*std::abs(nextB[i]-B[i]));scale=std::max(scale,pic_c*std::abs(nextB[i]));}
        out.nonlinear_relative=error/scale;out.iterations=iteration+1;
        out.field=std::move(field.first);out.field_receipt=std::move(field.second);out.motion.emplace(std::move(motion));
        // The immutable receipt must bind exactly the rounded total midpoint
        // field used by the accepted Maxwell work, not merely a nearby iterate.
        if(nextE==E&&nextB==B){converged=true;break;}guess=out.field;
    }
    hall::require(converged,"v2 joint did not reach identical rounded midpoint fields within 16 iterations");
    out.particle_work_defect_J=out.motion->budget().kinetic_change_J-out.field_receipt.plasma_work_J;
    hall::require(std::abs(out.particle_work_defect_J)<=1e-20+1e-10*(
        std::abs(out.motion->budget().kinetic_change_J)+std::abs(out.field_receipt.plasma_work_J)),"v2 joint particle/Maxwell work mismatch");
    out.heat=thermal.prepare(out.field_receipt,common.dt_s,field_state_sha256(out.field));
    out.combined_balance_J=out.field_receipt.balance_J+out.particle_work_defect_J+
        out.heat.thermal_delta_J-out.field_receipt.joule_J;
    const double transfers=std::abs(out.field_receipt.field_delta_J)+std::abs(out.motion->budget().kinetic_change_J)+
        std::abs(out.heat.thermal_delta_J)+std::abs(out.field_receipt.outer_loss_J)+std::abs(out.field_receipt.outer_supply_J)+std::abs(out.field_receipt.circuit_work_J);
    hall::require(std::abs(out.combined_balance_J)<=2e-20+1e-10*transfers+out.heat.roundoff_bound_J,"v2 joint complete field/motion/heat ledger");
    return out;
}
inline bool preflight_axial_joint(const axial::State&accepted,const Thermal&thermal,
    std::span<const native::ArenaRateEntry>private_census,const AxialJointTrial&trial){
    return trial.motion&&field_state_sha256(accepted)==trial.base_field_sha256&&
        thermal_state_sha256(thermal)==trial.base_thermal_sha256&&canonical_census_sha256(private_census)==trial.base_census_sha256&&
        thermal.preflight(trial.heat);
}
// Call only AFTER root's common preflight of arena, both RNGs, source and this
// field/heat token. It performs no allocation or I/O, and never publishes arena.
inline void publish_field_heat_prevalidated(axial::State&accepted,Thermal&thermal,AxialJointTrial&&trial)noexcept{
    static_assert(std::is_nothrow_swappable_v<axial::State>);
    if(!thermal.commit(std::move(trial.heat)))std::terminate();using std::swap;swap(accepted,trial.field);
}
} // namespace hall_driver
