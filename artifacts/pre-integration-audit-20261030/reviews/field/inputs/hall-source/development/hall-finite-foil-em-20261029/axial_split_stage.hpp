#pragma once
#include "axial_system.hpp"

namespace axial {
enum class ChemistrySide { before_transport,after_transport };
struct ChemistryStage {
    foil::NativeInterval common;
    ChemistrySide side;
    double physical_stage_s,local_exposure_s;
    bool operator==(const ChemistryStage&)const=default;
};
struct NeutralAnchor {std::string parent_id,owner_state_hash;V3 position_m;};
struct ChemistryPrepared {
    ChemistryStage stage;
    std::vector<GapParticle>charged;
    std::vector<NeutralAnchor>neutral_anchors;
    // Root owns and preflights this staged token, including inventory/RNG/event
    // journal. This field operator cannot commit or duplicate those owners.
    std::string source_token,other_energy_owner;
    double other_owner_energy_delta_J=0;
    V3 other_owner_momentum_delta_kg_m_s{};
};
struct ChemistryReceipt {
    double charged_delta_K_J=0,other_delta_J=0,energy_defect_J=0,charge_impulse_C=0;
    V3 charged_delta_P{},other_delta_P{};
    std::string source_token;
};
inline double charged_kinetic(std::span<const GapParticle>p){long double sum=0;for(const auto&v:p)sum+=v.weight*pic_kinetic_energy_J(v);return double(sum);}
inline V3 charged_momentum(std::span<const GapParticle>p){V3 result{};for(const auto&v:p)for(int d=0;d<3;++d)result[d]+=v.weight*v.mass_kg*v.proper_velocity_m_s[d];return result;}
inline ChemistryReceipt admit_chemistry(const Mesh&g,const ChemistryStage&request,
    std::span<const GapParticle>before,const ChemistryPrepared&trial){
    hall::require(trial.stage==request&&!trial.source_token.empty()&&!trial.other_energy_owner.empty()&&
        trial.other_energy_owner!="field"&&trial.other_energy_owner!="Maxwell"&&trial.other_energy_owner!="circuit",
        "chemistry stage or counterpart owner mismatch");
    hall::require(request.local_exposure_s==.5*request.common.dt_s&&request.physical_stage_s==
        (request.side==ChemistrySide::before_transport?request.common.begin_s:request.common.end_s),"chemistry physical stage/exposure mismatch");
    std::map<std::uint64_t,GapParticle>old;std::vector<V3>anchors;
    Vec q0(g.nn),q1(g.nn),gross(g.nn);std::set<std::uint64_t>ids;
    for(const auto&p:before){old.emplace(p.id,p);anchors.push_back(p.position_m);
        point_deposit(g,q0,p.position_m,p.weight*p.charge_C);point_deposit(g,gross,p.position_m,p.weight*std::abs(p.charge_C));}
    for(const auto&a:trial.neutral_anchors){hall::require(!a.parent_id.empty()&&!a.owner_state_hash.empty(),"unbound neutral reaction anchor");anchors.push_back(a.position_m);}
    for(const auto&p:trial.charged){hall::require(p.id>0&&ids.insert(p.id).second&&std::isfinite(p.weight)&&p.weight>0&&
        std::isfinite(p.mass_kg)&&p.mass_kg>0&&std::isfinite(p.charge_C)&&p.charge_C!=0,"invalid chemistry charged proposal");
        for(const auto&v:{p.position_m,p.proper_velocity_m_s})for(double x:v)hall::require(std::isfinite(x),"nonfinite chemistry proposal");
        if(old.count(p.id))hall::require(p.position_m==old.at(p.id).position_m,"SSA exposure advanced an orbit");
        else hall::require(std::find(anchors.begin(),anchors.end(),p.position_m)!=anchors.end(),"reaction product has no frozen participant anchor");
        point_deposit(g,q1,p.position_m,p.weight*p.charge_C);point_deposit(g,gross,p.position_m,p.weight*std::abs(p.charge_C));
    }
    ChemistryReceipt r;r.source_token=trial.source_token;r.charged_delta_K_J=charged_kinetic(trial.charged)-charged_kinetic(before);
    r.other_delta_J=trial.other_owner_energy_delta_J;r.energy_defect_J=r.charged_delta_K_J+r.other_delta_J;
    auto P0=charged_momentum(before),P1=charged_momentum(trial.charged);double momentum_scale=0,momentum_transfer=0,momentum_error=0;
    for(int d=0;d<3;++d){r.charged_delta_P[d]=P1[d]-P0[d];r.other_delta_P[d]=trial.other_owner_momentum_delta_kg_m_s[d];
        hall::require(std::isfinite(r.other_delta_P[d]),"nonfinite collision momentum owner");
        momentum_scale+=std::abs(P0[d])+std::abs(P1[d])+std::abs(r.other_delta_P[d]);
        momentum_transfer+=std::abs(r.charged_delta_P[d])+std::abs(r.other_delta_P[d]);
        momentum_error=std::max(momentum_error,std::abs(r.charged_delta_P[d]+r.other_delta_P[d]));}
    r.charge_impulse_C=maximum(subtract(q1,q0));
    hall::require(r.charge_impulse_C<=1e-12*std::max(maximum(gross),1e-30),"instantaneous reaction is not nodewise charge balanced; explicit compatible source current required");
    double energy_roundoff=128*std::numeric_limits<double>::epsilon()*(charged_kinetic(before)+charged_kinetic(trial.charged));
    hall::require(std::isfinite(r.other_delta_J)&&std::abs(r.energy_defect_J)<=energy_roundoff+1e-10*(std::abs(r.charged_delta_K_J)+std::abs(r.other_delta_J)),"collision charged/other-owner energy failed");
    hall::require(momentum_error<=128*std::numeric_limits<double>::epsilon()*momentum_scale+1e-10*momentum_transfer,"collision momentum counterpart failed");
    return r;
}
struct SplitTrial {
    JointTrial candidate;
    ChemistryPrepared before,after;
    ChemistryReceipt before_receipt,after_receipt;
    double charged_total_delta_K_J=0,total_balance_J=0;
};
using ChemistryPrepare=std::function<ChemistryPrepared(const ChemistryStage&,std::span<const GapParticle>)>;
// Strang operator composition for FULL-TIME x/u: C(dt/2) at t_n, EM(dt),
// C(dt/2) at t_(n+1). The SSA waiting-time coordinate is local exposure only;
// no trajectory is advanced during either frozen-position collision operator.
// This function prepares only. Root must preflight BOTH source tokens alongside
// candidate heat/field/arena before its one accepted transaction commits.
inline SplitTrial prepare_chemistry_em_chemistry(const Mesh&g,const Fields&fields,const State&accepted,
    std::span<const GapParticle>arena,Thermal&thermal,const CurrentTrace&trace,double dt,
    foil::NativePhase incoming_phase,const std::set<std::uint64_t>&active_layer_ids,const ChemistryPrepare&prepare){
    hall::require(incoming_phase==foil::NativePhase::full_time_work_conjugate_v1,"split stage rejects old native D-K-D phase");
    hall::require(active_layer_ids.empty(),"reference split coordinator is Joule-only; active source layers require lower transport/path ports and root source coordinator");
    foil::NativeInterval common{foil::FieldMode::maxwell_total,accepted.revision,0,accepted.time_s,accepted.time_s+dt,dt,
        arena_fingerprint(arena),field_fingerprint(accepted),g.json().dump(),ss304_source_seal_sha256,arena_fingerprint(arena),incoming_phase};
    common.validate();SplitTrial out;ChemistryStage first{common,ChemistrySide::before_transport,accepted.time_s,.5*dt};
    out.before=prepare(first,arena);out.before_receipt=admit_chemistry(g,first,arena,out.before);
    out.candidate=prepare_joint(g,fields,accepted,out.before.charged,thermal,trace,dt);
    auto end_binding=common;end_binding.trial_stage=1;end_binding.field_hash=field_fingerprint(out.candidate.field);
    end_binding.arena_hash=arena_fingerprint(out.candidate.transport.updates);
    ChemistryStage last{end_binding,ChemistrySide::after_transport,common.end_s,.5*dt};
    out.after=prepare(last,out.candidate.transport.updates);
    out.after_receipt=admit_chemistry(g,last,out.candidate.transport.updates,out.after);
    hall::require(out.before.source_token!=out.after.source_token,"collision stage token reused");
    out.charged_total_delta_K_J=charged_kinetic(out.after.charged)-charged_kinetic(arena);
    out.total_balance_J=out.candidate.total_balance_J+out.before_receipt.energy_defect_J+out.after_receipt.energy_defect_J;
    out.candidate.base_arena=arena_fingerprint(arena);out.candidate.transport.updates=out.after.charged;
    // The transport receipt retains EM-only work; collision energy remains in
    // these explicit receipts. No fictitious current or charge repair is added.
    return out;
}
} // namespace axial
