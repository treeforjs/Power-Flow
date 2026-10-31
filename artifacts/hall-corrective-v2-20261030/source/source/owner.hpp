#pragma once
// Hall corrective owner, v2. Derived from hall-field-correction-20261030/owner.hpp
// (commit 8dddec0); the exclusive-owner protocol (move-only privately built
// proposals, lifetime/generation tickets, final approval, nonallocating
// publication) is unchanged. Changes: binary identities hashed once per
// immutable object and cached (accepted field and heat digests are replaced by
// the candidate's digests at publication, by nothrow swaps); preflight still
// recomputes the accepted and candidate field/heat digests once; motion comes
// from an explicit provider (CPU sealed arithmetic or the device provider).
#include "census_charge.hpp"
#include <chrono>

namespace hall_corrective_v2 {
class Owner;class JointProposal;class MotionProvider;
JointProposal prepare_closed_proposal(Owner&,const native::GapPicOptions&,
    std::span<const native::ArenaRateEntry>,hall_driver::MotionBinding,const std::string&,MotionProvider&);
struct OwnerLifetime final {};

// Canonical binary heat identities. The material binding digest is the
// unchanged 8dddec0 value (SHA-256 of the material binding JSON text).
inline void encode_heat_state(Digest&d,const electrode_heat_v2::State&s){
    d.text("hall-v2-heat-slice-v1");d.f64s(s.T);
    d.u64(s.inventory.size());for(const auto&c:s.inventory){d.f64(c.hi);d.f64(c.lo);}
    d.u64(s.watermarks.size());for(const auto&[key,value]:s.watermarks){d.text(key);d.u64(value);}
    d.u64(s.cumulative.size());for(double x:s.cumulative)d.f64(x);
    d.f64(s.time);d.u64(s.revision);
}
inline std::string heat_sha256_of(const std::string&material,const hall_driver::Thermal&t){
    Digest d("hall-v2-heat-state-v1");d.text(material);const auto states=t.states();d.u64(states.size());
    for(const auto&s:states)encode_heat_state(d,s);return d.hex();
}
inline std::string heat_sha256_of(const std::string&material,const hall_driver::ThermalTrial&t){
    Digest d("hall-v2-heat-state-v1");d.text(material);d.u64(t.bodies.size());
    for(const auto&b:t.bodies)encode_heat_state(d,b.state());return d.hex();
}

// Motion for one fixed-point iterate. Providers produce identical receipts
// (sealed record bytes) and the endpoint free charge equal to charge_of_endpoints.
class MotionProvider {
public:
    virtual ~MotionProvider()=default;
    virtual const std::string& name()const noexcept=0;
    // Once per proposal, with the immutable post-first-exposure census.
    virtual void begin(const MotionContext&,std::span<const native::ArenaRateEntry> census)=0;
    virtual hall_driver::MotionReceipt motion(const hall_driver::MotionBinding&,const Vec&E,const Vec&B,FreeCharge&after)=0;
    virtual void end()noexcept{}
};
class CpuMotionProvider final:public MotionProvider {
    std::string name_="cpu-sealed-arithmetic";const MotionContext*context_=nullptr;std::span<const native::ArenaRateEntry> census_;
public:
    const std::string& name()const noexcept override{return name_;}
    void begin(const MotionContext&c,std::span<const native::ArenaRateEntry> census)override{context_=&c;census_=census;}
    hall_driver::MotionReceipt motion(const hall_driver::MotionBinding&b,const Vec&E,const Vec&B,FreeCharge&after)override{
        need(context_!=nullptr,"motion provider used outside a proposal");
        auto r=MotionAssembly::cpu(*context_,census_,b,{},E,B);
        after=charge_of_endpoints(*context_->mesh,context_->options,r.endpoints(),context_->maximum_carriers);return r;
    }
    void end()noexcept override{context_=nullptr;census_={};}
};

struct ProposalTimings {
    double admit_census_s=0,motion_s=0,support_s=0,field_prepare_s=0,solve_s=0,residual_check_s=0,seal_s=0,heat_s=0,total_s=0;
    double motion_v2_seal_s=0;   // JSON-text seal of the accepted motion (native admission), once per proposal
    unsigned iterations=0;double continuity_ratio=0;SolveRecord last_solve;std::string motion_provider;
};

class FieldHeatProposal {
    friend class Owner;friend class JointProposal;
    FieldCandidate field_;hall_driver::ThermalTrial heat_;
    std::weak_ptr<const OwnerLifetime> owner_;std::uint64_t generation_=0;
    std::string before_heat_,after_heat_;mutable bool final_approved_=false;bool live_=true;
    FieldHeatProposal(FieldCandidate&&f,hall_driver::ThermalTrial&&h,const std::string&material,
        const std::shared_ptr<const OwnerLifetime>&owner,std::uint64_t generation,std::string before_heat):
        field_(std::move(f)),heat_(std::move(h)),owner_(owner),generation_(generation),before_heat_(std::move(before_heat)),
        after_heat_(heat_sha256_of(material,heat_)){
        const auto&r=field_.receipt();const double transfers=std::abs(r.field_delta_J)+std::abs(r.plasma_work_J)+std::abs(heat_.thermal_delta_J)+
            std::abs(r.outer_loss_J)+std::abs(r.outer_supply_J)+std::abs(r.circuit_work_J);
        need(std::isfinite(field_heat_balance_J())&&std::abs(field_heat_balance_J())<=2e-20+1e-10*transfers+heat_.roundoff_bound_J,
            "unchanged field/heat component energy bound");
    }
public:
    FieldHeatProposal(const FieldHeatProposal&)=delete;FieldHeatProposal&operator=(const FieldHeatProposal&)=delete;
    FieldHeatProposal(FieldHeatProposal&&)noexcept=default;FieldHeatProposal&operator=(FieldHeatProposal&&)noexcept=default;
    const std::string&candidate_field_sha256()const noexcept{return field_.after_sha256();}
    const std::string&candidate_heat_sha256()const noexcept{return after_heat_;}
    double end_s()const noexcept{return field_.field().time_s;}
    std::uint64_t endpoint_revision()const noexcept{return field_.field().revision;}
    const Midpoint&midpoint()const noexcept{return field_.midpoint();}
    const std::string&configuration_sha256()const noexcept{return field_.configuration_sha256();}
    const axial::State&candidate_field()const noexcept{return field_.field();}
    const ChargeLedger&candidate_charge()const noexcept{return field_.charge();}
    const axial::Receipt&field_receipt()const noexcept{return field_.receipt();}
    const Constraints&candidate_constraints()const noexcept{return field_.constraints();}
    const SolveRecord&solve_record()const noexcept{return field_.solve();}
    const hall_driver::ThermalTrial&heat()const noexcept{return heat_;}
    double field_heat_balance_J()const noexcept{return field_.receipt().balance_J+heat_.thermal_delta_J-field_.receipt().joule_J;}
};

class Owner {
    friend JointProposal prepare_closed_proposal(Owner&,const native::GapPicOptions&,
        std::span<const native::ArenaRateEntry>,hall_driver::MotionBinding,const std::string&,MotionProvider&);
    std::shared_ptr<const Model> model_;const native::GapPicOptions geometry_;const MotionContext motion_;
    axial::State field_;ChargeLedger charge_;std::unique_ptr<hall_driver::Thermal> thermal_;
    const std::shared_ptr<const OwnerLifetime> lifetime_=std::make_shared<const OwnerLifetime>();
    std::uint64_t generation_=0;bool vacuum_=false;std::string field_sha256_,heat_sha256_;
    struct RestoreTag {bool after_local_reaction=false;};
    static std::shared_ptr<const Model> checked_model(std::shared_ptr<const Model>m){need(bool(m),"missing immutable model");return m;}
    // Numeric equality of every word (the 8dddec0 check compared JSON records).
    static bool same_words(const WordVec&a,const WordVec&b){if(a.size()!=b.size())return false;
        for(std::size_t i=0;i<a.size();++i)if(a[i].hi!=b[i].hi||a[i].lo!=b[i].lo)return false;return true;}
    Owner(std::shared_ptr<const Model>m,axial::State field,ChargeLedger charge,const native::GapPicOptions&o,
        std::span<const native::ArenaRateEntry>census,RestoreTag tag):model_(checked_model(std::move(m))),geometry_(o),
        motion_(make_motion_context(model_->mesh(),geometry_,model_->mesh_json(),model_->limits().maximum_carriers)),field_(std::move(field)),charge_(std::move(charge)),vacuum_(census.empty()){
        const auto actual=charge_of_census(mesh(),geometry_,census,model_->limits().maximum_carriers);
        need(tag.after_local_reaction?same_net_charge(charge_.particles,actual):
            (same_words(actual.net_C,charge_.particles.net_C)&&actual.gross_C==charge_.particles.gross_C),
            "initial/restored net charge differs from actual canonical carriers");
        model_->admit(field_,charge_);admit_census(census);
        thermal_=std::make_unique<hall_driver::Thermal>(mesh(),model_->drive().start());
        field_sha256_=model_->state_sha256(field_,charge_);heat_sha256_=heat_sha256_of(model_->material_sha256(),*thermal_);
    }
    void admit_census(std::span<const native::ArenaRateEntry>census)const{
        auto actual=charge_of_census(mesh(),geometry_,census,model_->limits().maximum_carriers);need(same_net_charge(charge_.particles,actual),"accepted charge does not match actual canonical census");
        auto ledger=charge_;ledger.particles=std::move(actual);model_->admit(field_,ledger);
    }
    std::string conductivity_sha256()const{return Model::conductivity_sha256(thermal_->conductivity());}
    FieldCandidate prepare_field(const Vec&plasma,const FreeCharge&after,const DepositSupport&support,const Vec&sigma,
                                 const std::string&sigma_sha256,double dt)const{
        return model_->prepare(field_,charge_,field_sha256_,plasma,after,support,sigma,sigma_sha256,dt);
    }
    void seal_field(FieldCandidate&f)const{model_->seal(f);}
    FieldHeatProposal heat_for(FieldCandidate&&field){
        if(!field.sealed())model_->seal(field);
        need(field.before_sha256()==field_sha256_&&field.configuration_sha256()==configuration_sha256()&&
            field.conductivity_sha256()==conductivity_sha256(),"field preparation is stale or belongs to another configuration");
        model_->admit(field.field(),field.charge());
        auto heat=thermal_->prepare(field.receipt(),field.dt_s(),field.after_sha256());
        return FieldHeatProposal(std::move(field),std::move(heat),model_->material_sha256(),lifetime_,generation_,heat_sha256_);
    }
    bool quick_preflight(const FieldHeatProposal&p)const noexcept{
        return p.live_&&p.owner_.lock()==lifetime_&&p.generation_==generation_&&thermal_->preflight(p.heat_);
    }
    void publish(FieldHeatProposal&&p)noexcept{
        static_assert(std::is_nothrow_swappable_v<axial::State>&&std::is_nothrow_swappable_v<ChargeLedger>&&
                      std::is_nothrow_swappable_v<std::string>);
        if(!p.final_approved_||!quick_preflight(p)||generation_==UINT64_MAX)std::terminate();
        if(!thermal_->commit(std::move(p.heat_)))std::terminate();using std::swap;
        swap(field_,p.field_.field_);swap(charge_,p.field_.charge_);
        swap(field_sha256_,p.field_.after_);swap(heat_sha256_,p.after_heat_);
        ++generation_;p.live_=false;p.final_approved_=false;
    }
    void clock_guard()const{need(field_.time_s==thermal_->time_s()&&field_.revision==thermal_->revision(),"accepted field/heat clock mismatch");}
public:
    // Bootstrap consumes an actual read-only census. The field component never
    // allocates carriers or IDs and never repairs its declared charge.
    Owner(std::shared_ptr<const Model>m,axial::State field,ChargeLedger charge,const native::GapPicOptions&o,
        std::span<const native::ArenaRateEntry>census):Owner(std::move(m),std::move(field),std::move(charge),o,census,RestoreTag{}){clock_guard();}
    Owner(const Owner&)=delete;Owner&operator=(const Owner&)=delete;Owner(Owner&&)=delete;Owner&operator=(Owner&&)=delete;
    const axial::State&field()const noexcept{return field_;}const ChargeLedger&charge()const noexcept{return charge_;}
    const Mesh&mesh()const noexcept{return model_->mesh();}double time_s()const noexcept{return field_.time_s;}
    const Model&model()const noexcept{return *model_;}
    std::uint64_t revision()const noexcept{return field_.revision;}
    // Cached digests of the accepted (immutable between publications) state.
    const std::string&field_sha256()const noexcept{return field_sha256_;}
    const std::string&heat_sha256()const noexcept{return heat_sha256_;}
    const std::string&configuration_sha256()const noexcept{return model_->configuration_sha256();}
    const std::string&motion_mesh_sha256()const noexcept{return motion_.mesh_sha256;}
    const MotionContext&motion_context()const noexcept{return motion_;}
    const hall_driver::Thermal&thermal()const noexcept{return *thermal_;}
    FieldHeatProposal prepare_vacuum(double dt){need(vacuum_,"field-only publication requires an empty actual census");clock_guard();
        const auto sigma=thermal_->conductivity();
        auto field=prepare_field(Vec(mesh().ne),charge_.particles,DepositSupport::empty(mesh()),sigma,Model::conductivity_sha256(sigma),dt);
        return heat_for(std::move(field));}
    bool preflight_field_heat(const FieldHeatProposal&p)const{
        try{if(!quick_preflight(p)||p.field_.before_sha256()!=field_sha256_||p.before_heat_!=heat_sha256_||
                p.configuration_sha256()!=configuration_sha256()||p.field_.conductivity_sha256()!=conductivity_sha256())return false;
            clock_guard();model_->admit(field_,charge_);model_->admit(p.candidate_field(),p.candidate_charge());
            // One independent recomputation of the accepted and candidate digests.
            if(model_->state_sha256(field_,charge_)!=field_sha256_||heat_sha256_of(model_->material_sha256(),*thermal_)!=heat_sha256_||
                p.candidate_field_sha256()!=model_->state_sha256(p.candidate_field(),p.candidate_charge())||
                p.candidate_heat_sha256()!=heat_sha256_of(model_->material_sha256(),p.heat_)||
                p.end_s()!=time_s()+p.field_.dt_s()||p.endpoint_revision()!=revision()+1||p.midpoint().E()!=p.field_.receipt().mid_E)return false;
            const auto&f=p.candidate_field();Vec E(mesh().ne),B(mesh().nb);
            for(int e=0;e<mesh().ne;++e)E[e]=((Words(field_.E[e],field_.E_low[e])+Words(f.E[e],f.E_low[e]))*.5).rounded();
            for(int b=0;b<mesh().nb;++b)B[b]=((Words(field_.B[b],field_.B_low[b])+Words(f.B[b],f.B_low[b]))*.5).rounded();
            if(E!=p.midpoint().E()||B!=p.midpoint().B()||Midpoint(std::move(E),std::move(B)).seal()!=p.midpoint().seal())return false;
            for(const auto&body:p.heat_.bodies)if(body.state().time!=p.end_s()||body.state().revision!=p.endpoint_revision())return false;
            return true;
        }catch(const std::exception&){return false;}
    }
    template<class Boundary>bool preflight_vacuum_final(const FieldHeatProposal&p,const Boundary&b)const{
        p.final_approved_=false;
        const bool ok=vacuum_&&preflight_field_heat(p)&&b.field_sha256==p.candidate_field_sha256()&&b.heat_sha256==p.candidate_heat_sha256()&&
            b.time_s==p.end_s()&&b.accepted_step==p.endpoint_revision();p.final_approved_=ok;return ok;
    }
    void publish_vacuum_prevalidated(FieldHeatProposal&&p)noexcept{if(!vacuum_)std::terminate();publish(std::move(p));}
    bool preflight(const JointProposal&,std::span<const native::ArenaRateEntry>)const;
    template<class Boundary>bool preflight_final(const JointProposal&,std::span<const native::ArenaRateEntry>,const Boundary&)const;
    void publish_prevalidated(JointProposal&&)noexcept;
    void save(const std::filesystem::path&dir)const;
    static std::unique_ptr<Owner> restore(std::shared_ptr<const Model>,const native::GapPicOptions&,
        std::span<const native::ArenaRateEntry>,const std::filesystem::path&);
};
} // namespace hall_corrective_v2
