#pragma once
#include "axial_particles.hpp"
#include "axial_thermal.hpp"
#include <bit>

namespace axial {
inline std::string fingerprint(std::span<const double> values){
    std::uint64_t h=14695981039346656037ULL;for(double x:values){auto bits=std::bit_cast<std::uint64_t>(x);
        for(int j=0;j<8;++j){h^=(bits>>(8*j))&255;h*=1099511628211ULL;}}
    return std::to_string(h);
}
inline std::string arena_fingerprint(std::span<const GapParticle>arena){
    Vec values;for(const auto&p:arena){values.push_back(std::bit_cast<double>(p.id));
        for(double x:p.position_m)values.push_back(x);for(double x:p.proper_velocity_m_s)values.push_back(x);
        values.push_back(p.charge_C);values.push_back(p.mass_kg);values.push_back(p.weight);}
    return fingerprint(values);
}
inline std::string field_fingerprint(const State&s){return fingerprint(s.E)+"/"+fingerprint(s.B)+"/"+fingerprint(s.E_low)+"/"+fingerprint(s.B_low)+"/"+fingerprint(std::span<const double>(&s.time_s,1))+"/"+std::to_string(s.revision);}
struct JointTrial {
    std::string base_field,base_arena;std::uint64_t revision=0;double time_s=0;
    State field;foil::NativeTransportTrial transport;ThermalTrial heat;Receipt receipt;
    double nonlinear_relative=0,particle_work_defect_J=0,total_balance_J=0;int iterations=0;
};
// Reference composition for a Joule-only pilot. The production common driver
// may instead use the prepare ports individually to add its source/event and
// boundary owners. There is no accepted particle inventory in this object.
inline JointTrial prepare_joint(const Mesh&g,const Fields&fields,const State&accepted,
    std::span<const GapParticle>arena,Thermal&thermal,const CurrentTrace&trace,double dt){
    hall::require(accepted.time_s==thermal.time_s()&&accepted.revision==thermal.revision(),"EM/thermal common clock mismatch");
    hall::require(dt>0&&dt<=1e-14&&accepted.time_s+dt>accepted.time_s,"common dt guard");
    const double current=trace.current(accepted.time_s+.5*dt);auto sigma=thermal.conductivity();State guess=accepted;
    JointTrial out;out.base_field=field_fingerprint(accepted);out.base_arena=arena_fingerprint(arena);
    out.revision=accepted.revision;out.time_s=accepted.time_s;bool converged=false;
    for(int iteration=0;iteration<16;++iteration){Vec E(g.ne),B(g.nb);
        for(int e=0;e<g.ne;++e)E[e]=.5*(accepted.E[e]+guess.E[e]);
        for(int b=0;b<g.nb;++b)B[b]=.5*(accepted.B[b]+guess.B[b]);
        foil::NativeInterval interval{foil::FieldMode::maxwell_total,accepted.revision,std::uint64_t(iteration),
            accepted.time_s,accepted.time_s+dt,dt,out.base_arena,fingerprint(E)+"/"+fingerprint(B),
            g.json().dump(),std::string(ss304_source_seal_sha256)+"/"+fingerprint(sigma),out.base_arena,foil::NativePhase::full_time_work_conjugate_v1};
        out.transport=prepare_native_total_field_transport(g,interval,arena,{},E,B);
        hall::require(out.transport.boundaries.empty(),"joint Joule-only pilot requires explicit boundary owner before collection");
        auto next=fields.prepare(accepted,out.transport.current_load_A_m,sigma,dt,current);
        out.field=std::move(next.first);out.receipt=std::move(next.second);
        double difference=0,scale=std::max(1.,maximum(out.field.E));
        for(int e=0;e<g.ne;++e)difference=std::max(difference,std::abs(out.field.E[e]-guess.E[e]));
        for(int b=0;b<g.nb;++b){difference=std::max(difference,pic_c*std::abs(out.field.B[b]-guess.B[b]));scale=std::max(scale,pic_c*std::abs(out.field.B[b]));}
        out.nonlinear_relative=difference/scale;out.iterations=iteration+1;guess=out.field;
        if(out.nonlinear_relative<=1e-13){converged=true;break;}
    }
    hall::require(converged,"axial common field/transport iteration failed");
    out.particle_work_defect_J=out.transport.kinetic_change_J-out.receipt.plasma_work_J;
    hall::require(std::abs(out.particle_work_defect_J)<=1e-20+1e-10*(std::abs(out.transport.kinetic_change_J)+std::abs(out.receipt.plasma_work_J)),"axial accepted total-field particle work failed");
    out.heat=thermal.prepare(out.receipt,dt,field_fingerprint(out.field));
    out.total_balance_J=out.receipt.balance_J+out.particle_work_defect_J+out.heat.thermal_delta_J-out.receipt.joule_J;
    const double transfers=std::abs(out.receipt.field_delta_J)+std::abs(out.transport.kinetic_change_J)+
        std::abs(out.heat.thermal_delta_J)+std::abs(out.receipt.outer_loss_J)+std::abs(out.receipt.outer_supply_J)+std::abs(out.receipt.circuit_work_J);
    hall::require(std::abs(out.total_balance_J)<=2e-20+1e-10*transfers+out.heat.roundoff_bound_J,
                  "axial field particle thermal circuit ledger failed");
    return out;
}
inline bool commit_joint(State&accepted,std::vector<GapParticle>&canonical_arena,Thermal&thermal,JointTrial&&trial){
    // All allocation/checks precede these nonthrowing swaps. No component may
    // advance asynchronously between preflight and this single-thread commit.
    if(accepted.revision!=trial.revision||accepted.time_s!=trial.time_s||
       field_fingerprint(accepted)!=trial.base_field||arena_fingerprint(canonical_arena)!=trial.base_arena||
       !thermal.preflight(trial.heat))return false;
    static_assert(std::is_nothrow_swappable_v<State>);
    if(!thermal.commit(std::move(trial.heat)))std::terminate();
    using std::swap;swap(accepted,trial.field);canonical_arena.swap(trial.transport.updates);return true;
}
} // namespace axial
