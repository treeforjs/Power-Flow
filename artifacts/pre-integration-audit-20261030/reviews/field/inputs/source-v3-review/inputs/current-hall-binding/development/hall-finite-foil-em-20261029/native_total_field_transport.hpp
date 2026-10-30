#pragma once
#include "foil_particles.hpp"
#include <span>

namespace foil {
// This API has no GapPic object and never calls solve_fields, electric_field,
// set_external_fields, or advance. Its input is the one canonical arena's
// immutable accepted view. Momentum uses native pic_boris; transport uses the
// separately qualified work displacement and conjugate path gather/deposition.
enum class FieldMode { legacy_electrostatic_plus_external, maxwell_total };
enum class NativePhase { unbound, legacy_leapfrog_dkd, full_time_work_conjugate_v1 };
struct NativeInterval {
    FieldMode mode=FieldMode::maxwell_total;
    std::uint64_t accepted_revision=0, trial_stage=0;
    double begin_s=0,end_s=0,dt_s=0;
    std::string arena_hash,field_hash,mesh_hash,material_hash,species_hash;
    NativePhase particle_phase=NativePhase::unbound;
    void validate() const {
        hall::require(mode==FieldMode::maxwell_total,"EM transport requires sole total-field ownership");
        hall::require(particle_phase==NativePhase::full_time_work_conjugate_v1,"native EM transport requires qualified full-time x/u phase; legacy D-K-D is not relabeled");
        hall::require(std::isfinite(begin_s)&&std::isfinite(end_s)&&end_s>begin_s&&
                      std::isfinite(dt_s)&&dt_s>0&&dt_s<=1e-14&&end_s==begin_s+dt_s,
                      "invalid native common interval");
        hall::require(!arena_hash.empty()&&!field_hash.empty()&&!mesh_hash.empty()&&!material_hash.empty()&&!species_hash.empty(),
                      "unbound native total-field stage");
    }
    bool operator==(const NativeInterval&)const=default;
};
struct NativeBoundary {
    GapParticle endpoint;
    int surface;
    double kinetic_J,flight_s;
};
struct NativeTransportTrial {
    NativeInterval interval;
    // Proposals only. The caller owns IDs, membership, time, and the sole commit.
    // Excluded layer-owned particles are absent from updates and current.
    std::vector<GapParticle> updates;
    std::vector<NativeBoundary> boundaries;
    Vec current_load_A_m,initial_charge_C,endpoint_charge_C;
    double field_work_J=0,kinetic_change_J=0,work_defect_J=0,continuity_C=0;
    double plasma_dt=0;
};

inline NativeTransportTrial prepare_native_total_field_transport(
    const Mesh& mesh,const NativeInterval& interval,std::span<const GapParticle> arena,
    const std::set<std::uint64_t>& layer_owned_ids,const Vec& total_mid_E_V_m,
    const Vec& total_mid_B_T) {
    interval.validate();
    hall::require(total_mid_E_V_m.size()==std::size_t(mesh.ne)&&total_mid_B_T.size()==std::size_t(mesh.nb)&&
                  finite(total_mid_E_V_m)&&finite(total_mid_B_T),"total-field mesh mismatch");
    const double dt=interval.dt_s;
    NativeTransportTrial out;out.interval=interval;out.current_load_A_m.resize(mesh.ne);
    out.initial_charge_C.resize(mesh.nn);out.endpoint_charge_C.resize(mesh.nn);
    Vec gross(mesh.nn),omega2_load(mesh.nn);std::set<std::uint64_t> ids;
    for(const auto&p:arena) {
        hall::require(p.id!=0&&ids.insert(p.id).second&&std::isfinite(p.weight)&&p.weight>0&&
                      std::isfinite(p.mass_kg)&&p.mass_kg>0&&std::isfinite(p.charge_C)&&p.charge_C!=0,
                      "invalid canonical charged particle");
        for(const auto&v:{p.position_m,p.proper_velocity_m_s})for(double x:v)
            hall::require(std::isfinite(x),"nonfinite canonical particle");
        point_deposit(mesh,omega2_load,p.position_m,p.weight*p.charge_C*p.charge_C/(eps*p.mass_kg));
    }
    for(auto id:layer_owned_ids)hall::require(ids.count(id),"layer exclusion is absent from canonical arena");
    for(int j=0;j<=mesh.ny;++j)for(int i=0;i<=mesh.nx;++i)
        out.plasma_dt=std::max(out.plasma_dt,dt*std::sqrt(omega2_load[mesh.node(i,j)]/
                              (mesh.dual_x(i)*mesh.dual_y(j)*mesh.depth)));
    hall::require(out.plasma_dt<=.1,"native total-field full-arena plasma guard");
    for(const auto&p:arena) {
        if(layer_owned_ids.count(p.id))continue;
        Species s{"canonical",p.charge_C,p.mass_kg};
        Particle input{p.id,0,p.position_m,p.proper_velocity_m_s,p.weight};
        auto trial=push(mesh,input,s,total_mid_E_V_m,total_mid_B_T,dt);
        Path path(mesh,p.position_m,trial.displacement);
        path.deposit(out.current_load_A_m,p.weight*p.charge_C,trial.displacement,dt);
        point_deposit(mesh,out.initial_charge_C,p.position_m,p.weight*p.charge_C);
        point_deposit(mesh,out.endpoint_charge_C,trial.particle.x,p.weight*p.charge_C);
        point_deposit(mesh,gross,p.position_m,p.weight*std::abs(p.charge_C));
        auto q=p;q.position_m=trial.particle.x;q.proper_velocity_m_s=trial.particle.u;
        out.kinetic_change_J+=p.weight*(pic_kinetic_energy_J(q)-pic_kinetic_energy_J(p));
        if(trial.hit.surface<0)out.updates.push_back(q);
        else out.boundaries.push_back({q,trial.hit.surface,p.weight*pic_kinetic_energy_J(q),trial.flight_s});
    }
    out.field_work_J=dt*dot(total_mid_E_V_m,out.current_load_A_m);
    out.work_defect_J=out.kinetic_change_J-out.field_work_J;
    auto flux=mesh.gradient_transpose(out.current_load_A_m);
    for(int n=0;n<mesh.nn;++n)out.continuity_C=std::max(out.continuity_C,
        std::abs(out.endpoint_charge_C[n]-out.initial_charge_C[n]-dt*flux[n]));
    hall::require(out.continuity_C<=1e-12*std::max(maximum(gross),1e-30),"native path continuity failed");
    hall::require(std::abs(out.work_defect_J)<=1e-20+1e-10*(std::abs(out.kinetic_change_J)+std::abs(out.field_work_J)),
                  "native particle total-field work mismatch");
    return out;
}
} // namespace foil
