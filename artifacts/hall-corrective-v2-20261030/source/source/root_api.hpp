#pragma once
// Hall corrective owner v2, entry point (namespace hall_corrective_v2).
// Derived from hall-field-correction-20261030/root_api.hpp (commit 8dddec0).
// Unchanged: the opaque JointProposal, the 16-iteration bitwise-stable
// compensated-midpoint fixed point, the particle/Maxwell work and complete
// energy bounds, preflight/preflight_final and nonallocating publication.
// Changes: explicit motion provider, deposit-rounding support for the endpoint
// guard, cached binary identities and per-phase timings of the proposal.
#include "owner.hpp"
#include "checkpoint.hpp"

namespace hall_corrective_v2 {
class JointProposal {
    friend class Owner;
    friend JointProposal prepare_closed_proposal(Owner&,const native::GapPicOptions&,
        std::span<const native::ArenaRateEntry>,hall_driver::MotionBinding,const std::string&,MotionProvider&);
    FieldHeatProposal field_heat_;hall_driver::MotionReceipt motion_;std::string source_;
    double combined_balance_=0;unsigned iterations_=0;ProposalTimings timings_;
    JointProposal(FieldHeatProposal&&p,hall_driver::MotionReceipt&&motion,std::string source,unsigned iterations,ProposalTimings timings):
        field_heat_(std::move(p)),motion_(std::move(motion)),source_(std::move(source)),iterations_(iterations),timings_(std::move(timings)){
        const auto&r=field_heat_.field_.receipt();const double work=motion_.budget().kinetic_change_J-r.plasma_work_J;
        need(std::isfinite(work)&&std::abs(work)<=1e-20+1e-10*(std::abs(motion_.budget().kinetic_change_J)+std::abs(r.plasma_work_J)),
            "unchanged particle/Maxwell work bound");
        combined_balance_=field_heat_.field_heat_balance_J()+work;
        const double transfers=std::abs(r.field_delta_J)+std::abs(motion_.budget().kinetic_change_J)+
            std::abs(field_heat_.heat_.thermal_delta_J)+std::abs(r.outer_loss_J)+std::abs(r.outer_supply_J)+std::abs(r.circuit_work_J);
        need(std::isfinite(combined_balance_)&&std::abs(combined_balance_)<=2e-20+1e-10*transfers+field_heat_.heat_.roundoff_bound_J,
            "unchanged complete field/motion/heat energy bound");
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
};

// Closed, source-free preparation only. Replays immutable actual carriers from
// the same post-first-exposure private census. This never creates an arena,
// allocates an ID, samples a source, runs native Poisson/DKD, or publishes.
inline JointProposal prepare_closed_proposal(Owner&owner,const native::GapPicOptions&o,
    std::span<const native::ArenaRateEntry>census,hall_driver::MotionBinding common,const std::string&source,MotionProvider&provider){
    using clock=std::chrono::steady_clock;const auto since=[](clock::time_point a){return std::chrono::duration<double>(clock::now()-a).count();};
    const auto start=clock::now();ProposalTimings T;T.motion_provider=provider.name();
    owner.clock_guard();owner.admit_census(census);need(digest(source),"closed proposal requires unchanged source identity");
    T.admit_census_s=since(start);
    need(common.phase==hall_driver::ParticlePhase::full_time_position_and_proper_velocity&&
        common.accepted_step==owner.revision()&&common.begin_s==owner.time_s()&&common.end_s==common.begin_s+common.dt_s&&
        std::isfinite(common.dt_s)&&common.dt_s>0&&common.dt_s<=1e-14,"closed proposal common clock mismatch");
    need(common.mesh_sha256==owner.motion_.mesh_sha256&&owner.model_->motion_mesh_sha256(o)==owner.motion_.mesh_sha256&&
        common.material_sha256==owner.model_->material_sha256(),"closed proposal mesh/material identity mismatch");
    const auto sigma=owner.thermal_->conductivity();const auto sigma_sha=Model::conductivity_sha256(sigma);
    auto midpoint=owner.model_->at(owner.field_);
    provider.begin(owner.motion_,census);
    struct End {MotionProvider&p;~End(){p.end();}} end{provider};
    for(unsigned iteration=0;iteration<16;++iteration){
        auto binding=common;binding.trial_field_sha256=hall_driver::total_field_sha256(midpoint.E(),midpoint.B());
        auto t=clock::now();FreeCharge after;
        auto motion=provider.motion(binding,midpoint.E(),midpoint.B(),after);T.motion_s+=since(t);t=clock::now();
        const auto support=DepositSupportBuilder::build(owner.mesh(),owner.geometry_,census,motion.endpoints());T.support_s+=since(t);t=clock::now();
        auto field=owner.prepare_field(motion.current_A_m(),after,support,sigma,sigma_sha,common.dt_s);T.field_prepare_s+=since(t);
        T.solve_s+=field.solve().solve_s;T.residual_check_s+=field.solve().check_s;T.last_solve=field.solve();
        T.continuity_ratio=std::max(T.continuity_ratio,field.continuity_ratio());
        if(midpoint.E()==field.midpoint().E()&&midpoint.B()==field.midpoint().B()){
            t=clock::now();MotionAssembly::seal_v2(motion);T.motion_v2_seal_s=since(t);
            t=clock::now();owner.seal_field(field);T.seal_s=since(t);t=clock::now();
            auto heat=owner.heat_for(std::move(field));T.heat_s=since(t);T.iterations=iteration+1;T.total_s=since(start);
            return JointProposal(std::move(heat),std::move(motion),source,iteration+1,std::move(T));
        }
        midpoint=field.midpoint();
    }
    throw std::runtime_error("closed proposal did not reach the identical authoritative midpoint within 16 iterations");
}
// Default provider: the sealed CPU arithmetic.
inline JointProposal prepare_closed_proposal(Owner&owner,const native::GapPicOptions&o,
    std::span<const native::ArenaRateEntry>census,hall_driver::MotionBinding common,const std::string&source){
    CpuMotionProvider cpu;return prepare_closed_proposal(owner,o,census,std::move(common),source,cpu);
}
inline bool Owner::preflight(const JointProposal&p,std::span<const native::ArenaRateEntry>census)const{
    try{if(!preflight_field_heat(p.field_heat_))return false;admit_census(census);const auto&m=p.motion();const auto&b=m.binding();
        if(m.census_sha256()!=canonical_census_sha256(census,motion_.maximum_carriers)||!digest(m.seal_sha256())||
            !digest(m.record_v3_sha256())||MotionAssembly::record_v3(m)!=m.record_v3_sha256()||
            b.phase!=hall_driver::ParticlePhase::full_time_position_and_proper_velocity||b.accepted_step!=revision()||
            b.begin_s!=time_s()||b.end_s!=p.end_s()||b.dt_s!=p.field_heat_.field_.dt_s()||
            b.mesh_sha256!=motion_.mesh_sha256||b.material_sha256!=model_->material_sha256()||
            b.kinetic_policy_sha256!=motion_.policy_sha256||
            b.trial_field_sha256!=hall_driver::total_field_sha256(p.midpoint().E(),p.midpoint().B())||
            m.current_A_m()!=p.field_heat_.field_.current())return false;
        const auto actual=charge_of_endpoints(mesh(),geometry_,m.endpoints(),motion_.maximum_carriers);
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
} // namespace hall_corrective_v2
