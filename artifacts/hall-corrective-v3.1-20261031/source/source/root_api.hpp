#pragma once
// Hall corrective owner v3, entry point (namespace hall_corrective_v3).
// Derived from Hall corrective v2 (hall-field-correction-20261030/root_api.hpp,
// commit 8dddec0). Unchanged: the opaque JointProposal, the 16-iteration
// bitwise-stable compensated-midpoint fixed point, the particle/Maxwell work and
// complete energy bounds, preflight/preflight_final and nonallocating publication.
// v2: explicit motion provider, deposit-rounding support for the endpoint guard,
// cached binary identities and per-phase timings of the proposal.
// v3: first-hit kinetic sinks with conductor/exterior bank transfers in the same
// field step, surface births, Hall-owned collected wall heat and external
// surface loads for the one thermal owner, the step energy ledger (StepEnergy)
// with every new transfer counted once. Without sinks, births and loads every
// value is computed by the v2 expressions.
#include "owner.hpp"
#include "checkpoint.hpp"

namespace hall_corrective_v3 {
// Step energy ledger (J). Field: field_delta + plasma_work + joule + outer_loss -
// outer_supply - circuit = 0. Plasma: kinetic_change = plasma_work (newborns over
// [tau, end] included). Kinetic store: kinetic_change - sum(sink incident kinetic)
// + sum(newborn kinetic at birth, supplied by the source).
// Thermal: thermal_delta = joule + wall_heat + external (+ roundoff). Every sink's
// impact kinetic energy is exactly one of wall heat (foils), short-bridge energy
// (no thermal owner; ideal PEC sheet) or escaped energy. residual_J closes all of it.
struct StepEnergy {
    double field_delta_J=0,plasma_work_J=0,joule_J=0,outer_loss_J=0,outer_supply_J=0,circuit_work_J=0;
    double kinetic_change_J=0,plasma_kinetic_store_change_J=0;
    double thermal_delta_J=0,wall_heat_J=0,external_J=0,collected_kinetic_J=0,short_bridge_kinetic_J=0,escaped_kinetic_J=0;
    double birth_kinetic_J=0,birth_supplied_J=0,combined_balance_J=0,residual_J=0,scale_J=0;
    static StepEnergy of(const axial::Receipt&r,const hall_driver::MotionReceipt&m,const ThermalTrialV3&h,double combined,double scale){
        StepEnergy e;e.field_delta_J=r.field_delta_J;e.plasma_work_J=r.plasma_work_J;e.joule_J=r.joule_J;e.outer_loss_J=r.outer_loss_J;
        e.outer_supply_J=r.outer_supply_J;e.circuit_work_J=r.circuit_work_J;e.kinetic_change_J=m.budget().kinetic_change_J;
        const auto sinks=SinkLedger::of(m.sinks());const auto births=BirthLedger::of(m.births());
        e.collected_kinetic_J=sinks.collected_kinetic_J().rounded();e.short_bridge_kinetic_J=sinks.boundary[3].kinetic_J.rounded();
        e.escaped_kinetic_J=sinks.escaped_kinetic_J().rounded();
        const Words bk=births.kinetic_J();e.birth_kinetic_J=bk.rounded();e.birth_supplied_J=e.birth_kinetic_J;
        const Words all_sinks=sinks.collected_kinetic_J()+sinks.escaped_kinetic_J();
        e.plasma_kinetic_store_change_J=(Words(e.kinetic_change_J)-all_sinks+bk).rounded();
        e.thermal_delta_J=h.thermal_delta_J;e.wall_heat_J=h.wall_heat_J;e.external_J=h.external_J;e.combined_balance_J=combined;
        const long double res=(long double)e.field_delta_J+e.plasma_kinetic_store_change_J+e.thermal_delta_J+e.escaped_kinetic_J+e.short_bridge_kinetic_J+
            e.outer_loss_J-e.outer_supply_J-e.circuit_work_J-e.external_J-e.birth_supplied_J;
        e.residual_J=double(res);e.scale_J=scale+std::abs(e.escaped_kinetic_J)+std::abs(e.short_bridge_kinetic_J)+std::abs(e.birth_kinetic_J)+std::abs(e.collected_kinetic_J);
        return e;
    }
};
// Thermal loads of one step: Hall's own wall heat of carriers collected on the
// foils (incoming kinetic energy at the impact point, foil slice of the hit z),
// then the external loads in the caller's order.
// v3.1: generic over the owner type (Owner, ResidentOwner); unchanged body.
template<class OwnerT>inline std::vector<SurfaceLoad> step_loads(const OwnerT&owner,const hall_driver::MotionReceipt&m,std::span<const SurfaceLoad>external,
                                           const hall_driver::MotionBinding&b){
    std::vector<SurfaceLoad> out;if(m.sinks().empty()&&external.empty())return out;
    const auto&g=owner.mesh();need(owner.revision()<(std::uint64_t(1)<<43),"wall heat sequence space");
    const std::uint64_t base=(owner.revision()+1)<<20;std::uint64_t i=0;
    for(const auto&k:m.sinks()){if(!k.collected||k.conductor>1)continue;
        need(i+1<(std::uint64_t(1)<<20),"more than 2^20 collected carriers in one step");
        SurfaceLoad l;l.body=k.conductor;l.slice=std::clamp(int(std::floor((k.field_point_m[2]-g.zmin)/g.dz)),g.kshort,g.kfeed-1);l.channel=wall_heat_channel;
        auto&x=l.transfer;x.producer=wall_heat_producer;x.energy_owner=wall_heat_owner;x.sequence=x.owner_sequence=base+(++i);
        x.state_id=k.state.empty()?"kinetic-sink":k.state;x.parent_id=k.packet_id.empty()?"kinetic-sink":k.packet_id;x.channel_id=wall_heat_channel;
        x.event_time_s=b.begin_s+k.flight_s;x.surface=true;x.x=k.field_point_m[0];x.y=k.field_point_m[1];x.incoming_K_J=k.weighted_kinetic_J;
        out.push_back(std::move(l));}
    out.insert(out.end(),external.begin(),external.end());
    return out;
}
class JointProposal {
    friend class Owner;
    friend JointProposal prepare_closed_proposal(Owner&,const native::GapPicOptions&,
        std::span<const native::ArenaRateEntry>,hall_driver::MotionBinding,const std::string&,MotionProvider&,const StepInputs&);
    FieldHeatProposal field_heat_;hall_driver::MotionReceipt motion_;std::string source_;
    double combined_balance_=0;unsigned iterations_=0;ProposalTimings timings_;StepEnergy energy_;
    JointProposal(FieldHeatProposal&&p,hall_driver::MotionReceipt&&motion,std::string source,unsigned iterations,ProposalTimings timings):
        field_heat_(std::move(p)),motion_(std::move(motion)),source_(std::move(source)),iterations_(iterations),timings_(std::move(timings)){
        const auto&r=field_heat_.field_.receipt();const auto&h=field_heat_.heat_;
        const bool births=!motion_.births().empty();
        const double work=motion_.budget().kinetic_change_J-r.plasma_work_J;
        need(std::isfinite(work)&&std::abs(work)<=1e-20+1e-10*(std::abs(motion_.budget().kinetic_change_J)+std::abs(r.plasma_work_J)),
            "unchanged particle/Maxwell work bound");
        combined_balance_=field_heat_.field_heat_balance_J()+work;
        double transfers=std::abs(r.field_delta_J)+std::abs(motion_.budget().kinetic_change_J)+
            std::abs(h.thermal_delta_J)+std::abs(r.outer_loss_J)+std::abs(r.outer_supply_J)+std::abs(r.circuit_work_J);
        if(births||h.loads)transfers+=std::abs(h.wall_heat_J)+h.external_heating_J+h.external_cooling_J;
        need(std::isfinite(combined_balance_)&&std::abs(combined_balance_)<=2e-20+1e-10*transfers+h.roundoff_bound_J,
            "unchanged complete field/motion/heat energy bound");
        energy_=StepEnergy::of(r,motion_,h,combined_balance_,transfers);
        // Collected kinetic energy on the foils is exactly the Hall wall heat applied to the thermal owner.
        const auto sinks=SinkLedger::of(motion_.sinks());const double foil=(sinks.boundary[0].kinetic_J+sinks.boundary[1].kinetic_J).rounded();
        need(std::abs(foil-h.wall_heat_J)<=1e-12*std::abs(foil)+1e-300,"collected foil kinetic energy differs from the applied wall heat");
    }
public:
    JointProposal(const JointProposal&)=delete;JointProposal&operator=(const JointProposal&)=delete;
    JointProposal(JointProposal&&)noexcept=default;JointProposal&operator=(JointProposal&&)noexcept=default;
    const std::string&candidate_field_sha256()const noexcept{return field_heat_.candidate_field_sha256();}
    const std::string&candidate_heat_sha256()const noexcept{return field_heat_.candidate_heat_sha256();}
    double end_s()const noexcept{return field_heat_.end_s();}
    std::uint64_t endpoint_revision()const noexcept{return field_heat_.endpoint_revision();}
    const Midpoint&midpoint()const noexcept{return field_heat_.midpoint();}
    const hall_driver::MotionReceipt&motion()const noexcept{return motion_;}
    double combined_balance_J()const noexcept{return combined_balance_;}
    const std::string&configuration_sha256()const noexcept{return field_heat_.configuration_sha256();}
    unsigned iterations()const noexcept{return iterations_;}
    const FieldHeatProposal&field_heat()const noexcept{return field_heat_;}
    const ProposalTimings&timings()const noexcept{return timings_;}
    const StepEnergy&energy()const noexcept{return energy_;}
    std::span<const SinkRecord> sinks()const noexcept{return motion_.sinks();}
    std::span<const NewbornRecord> births()const noexcept{return motion_.births();}
};

// Closed, source-free preparation only. Replays immutable actual carriers from
// the same post-first-exposure private census. This never creates an arena,
// allocates an ID, samples a source, runs native Poisson/DKD, or publishes.
inline JointProposal prepare_closed_proposal(Owner&owner,const native::GapPicOptions&o,
    std::span<const native::ArenaRateEntry>census,hall_driver::MotionBinding common,const std::string&source,MotionProvider&provider,
    const StepInputs&inputs){
    using clock=std::chrono::steady_clock;const auto since=[](clock::time_point a){return std::chrono::duration<double>(clock::now()-a).count();};
    const auto start=clock::now();ProposalTimings T;T.motion_provider=provider.name();
    owner.clock_guard();owner.admit_census(residents(census,inputs.newborns));need(digest(source),"closed proposal requires unchanged source identity");
    T.admit_census_s=since(start);
    need(common.phase==hall_driver::ParticlePhase::full_time_position_and_proper_velocity&&
        common.accepted_step==owner.revision()&&common.begin_s==owner.time_s()&&common.end_s==common.begin_s+common.dt_s&&
        std::isfinite(common.dt_s)&&common.dt_s>0&&common.dt_s<=owner.model_->dt_max(),"closed proposal common clock mismatch");
    need(common.mesh_sha256==owner.motion_.mesh_sha256&&owner.model_->motion_mesh_sha256(o)==owner.motion_.mesh_sha256&&
        common.material_sha256==owner.model_->material_sha256(),"closed proposal mesh/material identity mismatch");
    const auto sigma=owner.thermal_->conductivity();const auto sigma_sha=Model::conductivity_sha256(sigma);
    auto midpoint=owner.model_->at(owner.field_);
    provider.set_newborns(inputs.newborns,inputs.birth_scope);provider.begin(owner.motion_,census);
    struct End {MotionProvider&p;~End(){p.end();}} end{provider};
    for(unsigned iteration=0;iteration<16;++iteration){
        auto binding=common;binding.trial_field_sha256=midpoint.total_field_sha256();
        provider.set_verified_field(&midpoint.E(),&midpoint.B(),binding.trial_field_sha256);
        auto t=clock::now();FreeCharge after;
        auto motion=provider.motion(binding,midpoint.E(),midpoint.B(),after);
        need(motion.births().size()==inputs.newborns.size(),"motion provider did not move every newborn");
        T.motion_s+=since(t);t=clock::now();
        const auto bank=bank_transfer(owner.mesh(),owner.geometry_,motion,owner.model_->exterior_support());
        const auto support=DepositSupportBuilder::build(owner.mesh(),owner.geometry_,census,motion);T.support_s+=since(t);t=clock::now();
        auto field=owner.prepare_field(motion.current_A_m(),after,support,sigma,sigma_sha,common.dt_s,&bank);T.field_prepare_s+=since(t);
        T.solve_s+=field.solve().solve_s;T.residual_check_s+=field.solve().check_s;T.last_solve=field.solve();
        T.continuity_ratio=std::max(T.continuity_ratio,field.continuity_ratio());
        if(midpoint.E()==field.midpoint().E()&&midpoint.B()==field.midpoint().B()){
            t=clock::now();MotionAssembly::seal_v2(motion);T.motion_v2_seal_s=since(t);
            t=clock::now();owner.seal_field(field);T.seal_s=since(t);t=clock::now();
            const auto loads=step_loads(owner,motion,inputs.loads,common);
            auto heat=owner.heat_for(std::move(field),loads,&motion);T.heat_s=since(t);T.iterations=iteration+1;T.total_s=since(start);
            return JointProposal(std::move(heat),std::move(motion),source,iteration+1,std::move(T));
        }
        midpoint=field.midpoint();
    }
    throw std::runtime_error("closed proposal did not reach the identical authoritative midpoint within 16 iterations");
}
// v2 signatures: no births and no external loads.
inline JointProposal prepare_closed_proposal(Owner&owner,const native::GapPicOptions&o,
    std::span<const native::ArenaRateEntry>census,hall_driver::MotionBinding common,const std::string&source,MotionProvider&provider){
    return prepare_closed_proposal(owner,o,census,std::move(common),source,provider,StepInputs{});
}
// Default provider: the sealed CPU arithmetic.
inline JointProposal prepare_closed_proposal(Owner&owner,const native::GapPicOptions&o,
    std::span<const native::ArenaRateEntry>census,hall_driver::MotionBinding common,const std::string&source){
    CpuMotionProvider cpu;return prepare_closed_proposal(owner,o,census,std::move(common),source,cpu,StepInputs{});
}
inline bool Owner::preflight(const JointProposal&p,std::span<const native::ArenaRateEntry>census)const{
    try{if(!preflight_field_heat(p.field_heat_))return false;const auto&m=p.motion();const auto&b=m.binding();
        {std::vector<Newborn> nb;for(const auto&x:m.births())nb.push_back({x.native_id,x.birth_time_s});admit_census(residents(census,nb));}
        if(m.census_sha256()!=canonical_census_sha256(census,motion_.maximum_carriers)||!digest(m.seal_sha256())||
            !digest(m.record_v3_sha256())||MotionAssembly::record_v3(m)!=m.record_v3_sha256()||
            b.phase!=hall_driver::ParticlePhase::full_time_position_and_proper_velocity||b.accepted_step!=revision()||
            b.begin_s!=time_s()||b.end_s!=p.end_s()||b.dt_s!=p.field_heat_.field_.dt_s()||
            b.mesh_sha256!=motion_.mesh_sha256||b.material_sha256!=model_->material_sha256()||
            b.kinetic_policy_sha256!=motion_.policy_sha256||
            b.trial_field_sha256!=p.midpoint().total_field_sha256()||
            m.current_A_m()!=p.field_heat_.field_.current())return false;
        auto actual=charge_of_endpoints(mesh(),geometry_,m.endpoints(),motion_.maximum_carriers);
        if(!m.v2_identical()&&(MotionAssembly::sink_birth_digest(m).empty()||
            p.field_heat_.candidate_ledgers().sinks.count()!=kinetic_ledgers().sinks.count()+m.sinks().size()||
            p.field_heat_.candidate_ledgers().births.count()!=kinetic_ledgers().births.count()+m.births().size()))return false;
        if(!same_words(actual.net_C,p.field_heat_.candidate_charge().particles.net_C)||
            actual.gross_C!=p.field_heat_.candidate_charge().particles.gross_C)return false;
        return digest(p.source_)&&std::isfinite(p.combined_balance_J());
    }catch(const std::exception&){return false;}
}
// The caller MUST pass the read-only boundary returned by its actual opaque
// native candidate. This predicate compares values; native token authenticity,
// provider-seal binding and native/parent can_publish remain root-owned.
template<class Boundary>bool Owner::preflight_final(const JointProposal&p,
    std::span<const native::ArenaRateEntry>census,const Boundary&b)const{
    p.field_heat_.final_approved_=false;const auto&m=p.motion().binding();
    const bool ok=preflight(p,census)&&b.field_sha256==p.candidate_field_sha256()&&b.heat_sha256==p.candidate_heat_sha256()&&
        b.source_sha256==p.source_&&b.time_s==p.end_s()&&b.accepted_step==p.endpoint_revision()&&
        b.mesh_sha256==m.mesh_sha256&&b.material_sha256==m.material_sha256&&b.common_clock_sha256==m.common_clock_sha256;
    p.field_heat_.final_approved_=ok;return ok;
}
inline void Owner::publish_prevalidated(JointProposal&&p)noexcept{publish(std::move(p.field_heat_));}
} // namespace hall_corrective_v3
