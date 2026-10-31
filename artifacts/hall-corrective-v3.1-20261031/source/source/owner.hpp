#pragma once
// Hall corrective owner, v3 (derived from v2, itself derived from
// hall-field-correction-20261030/owner.hpp, commit 8dddec0); the exclusive-owner
// protocol (move-only privately built proposals, lifetime/generation tickets,
// final approval, nonallocating publication) is unchanged.
// v3: the thermal owner is ThermalV3 (one electrode thermal owner with external
// surface loads); the owner keeps cumulative kinetic sink, surface birth and
// surface load ledgers (accepted state, swapped at publication) with their own
// digest; field-only steps may carry external surface loads.
// v2: binary identities hashed once per immutable object and cached; preflight
// recomputes the accepted and candidate field/heat digests once; motion comes
// from an explicit provider (CPU sealed arithmetic or the device provider).
#include "census_charge.hpp"
#include <chrono>

namespace hall_corrective_v3 {
class Owner;class JointProposal;class MotionProvider;
// v3 step inputs beyond the census: the newborns of the step (native v6
// trial.birth_flights(): ID and birth time; they are in the census at their birth
// positions) and external surface loads for the thermal owner.
// birth_scope: the native step carries a birth scope (arrival committed every step,
// possibly with no births): the motion record is the flows record even when empty.
struct StepInputs {std::span<const Newborn> newborns;std::span<const SurfaceLoad> loads;bool birth_scope=false;};
// Census rows that are not newborns of this step (the accepted free charge covers these).
inline std::vector<native::ArenaRateEntry> residents(std::span<const native::ArenaRateEntry> census,std::span<const Newborn> newborns){
    std::vector<native::ArenaRateEntry> out;if(newborns.empty()){out.assign(census.begin(),census.end());return out;}
    std::set<std::uint64_t> ids;for(const auto&b:newborns)ids.insert(b.native_id);out.reserve(census.size());
    for(const auto&r:census)if(!ids.contains(r.native.id))out.push_back(r);
    return out;
}
JointProposal prepare_closed_proposal(Owner&,const native::GapPicOptions&,
    std::span<const native::ArenaRateEntry>,hall_driver::MotionBinding,const std::string&,MotionProvider&,const StepInputs&);
struct OwnerLifetime final {};
// v3.1: motion context carrying the configured field-step upper bound.
inline MotionContext with_dt_max(MotionContext c,double dt_max){c.dt_max=dt_max;return c;}

// Canonical binary heat identities. The material binding digest is the
// unchanged 8dddec0 value (SHA-256 of the material binding JSON text).
inline void encode_heat_state(Digest&d,const electrode_heat_v2::State&s){
    d.text("hall-v2-heat-slice-v1");d.f64s(s.T);
    d.u64(s.inventory.size());for(const auto&c:s.inventory){d.f64(c.hi);d.f64(c.lo);}
    d.u64(s.watermarks.size());for(const auto&[key,value]:s.watermarks){d.text(key);d.u64(value);}
    d.u64(s.cumulative.size());for(double x:s.cumulative)d.f64(x);
    d.f64(s.time);d.u64(s.revision);
}
inline std::string heat_sha256_of(const std::string&material,const ThermalV3&t){
    Digest d("hall-v3-heat-state-v1");d.text(material);const auto states=t.states();d.u64(states.size());
    for(const auto&s:states)encode_heat_state(d,s);return d.hex();
}
inline std::string heat_sha256_of(const std::string&material,const ThermalTrialV3&t){
    Digest d("hall-v3-heat-state-v1");d.text(material);d.u64(t.bodies.size());
    for(const auto&b:t.bodies)encode_heat_state(d,b.state());return d.hex();
}
// Cumulative accepted ledgers of kinetic sinks, surface births and surface loads.
struct KineticLedgers {
    SinkLedger sinks;BirthLedger births;std::map<std::string,Words> load_channel_J;std::uint64_t loads=0;
    Words wall_heat_J,external_J,short_bridge_J,escaped_J;
    std::string sha256()const{Digest d("hall-v3-kinetic-ledgers-v1");sinks.encode(d);births.encode(d);d.u64(loads);d.u64(load_channel_J.size());
        for(const auto&[k,v]:load_channel_J){d.text(k);d.f64(v.hi);d.f64(v.lo);}
        for(const auto*w:{&wall_heat_J,&external_J,&short_bridge_J,&escaped_J}){d.f64(w->hi);d.f64(w->lo);}return d.hex();}
};

// Motion for one fixed-point iterate. Providers produce identical receipts
// (sealed record bytes) and the endpoint free charge equal to charge_of_endpoints.
class MotionProvider {
protected:
    std::span<const Newborn> newborns_;bool birth_scope_=false;
    VerifiedField verified_;   // v3.1: total-field digest of the exact midpoint objects of the next motion() call
public:
    virtual ~MotionProvider()=default;
    // v3: newborns of the next proposal (flight end - tau, possible foil-face start).
    void set_newborns(std::span<const Newborn> n,bool birth_scope=false)noexcept{newborns_=n;birth_scope_=birth_scope||!n.empty();}
    // v3.1: the caller computed binding.trial_field_sha256 from these exact objects (see VerifiedField).
    void set_verified_field(const Vec*E,const Vec*B,const std::string&digest){verified_.E=E;verified_.B=B;verified_.digest=digest;}
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
        const VerifiedField v=std::exchange(verified_,VerifiedField{});   // one call only
        auto r=MotionAssembly::cpu(*context_,census_,b,{},E,B,newborns_,birth_scope_,&v);
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
    FieldCandidate field_;ThermalTrialV3 heat_;
    std::weak_ptr<const OwnerLifetime> owner_;std::uint64_t generation_=0;
    std::string before_heat_,after_heat_;mutable bool final_approved_=false;bool live_=true;
    KineticLedgers ledgers_;std::string ledgers_sha256_;
    FieldHeatProposal(FieldCandidate&&f,ThermalTrialV3&&h,const std::string&material,
        const std::shared_ptr<const OwnerLifetime>&owner,std::uint64_t generation,std::string before_heat,KineticLedgers ledgers):
        field_(std::move(f)),heat_(std::move(h)),owner_(owner),generation_(generation),before_heat_(std::move(before_heat)),
        after_heat_(heat_sha256_of(material,heat_)),ledgers_(std::move(ledgers)),ledgers_sha256_(ledgers_.sha256()){
        const auto&r=field_.receipt();double transfers=std::abs(r.field_delta_J)+std::abs(r.plasma_work_J)+std::abs(heat_.thermal_delta_J)+
            std::abs(r.outer_loss_J)+std::abs(r.outer_supply_J)+std::abs(r.circuit_work_J);
        if(heat_.loads)transfers+=std::abs(heat_.wall_heat_J)+heat_.external_heating_J+heat_.external_cooling_J;
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
    const ThermalTrialV3&heat()const noexcept{return heat_;}
    const KineticLedgers&candidate_ledgers()const noexcept{return ledgers_;}
    const std::string&candidate_ledgers_sha256()const noexcept{return ledgers_sha256_;}
    // v3: the thermal change is Joule + collected wall heat + external surface loads.
    double field_heat_balance_J()const noexcept{
        if(!heat_.loads)return field_.receipt().balance_J+heat_.thermal_delta_J-field_.receipt().joule_J;
        return field_.receipt().balance_J+heat_.thermal_delta_J-field_.receipt().joule_J-heat_.wall_heat_J-heat_.external_J;}
};

class Owner {
    friend JointProposal prepare_closed_proposal(Owner&,const native::GapPicOptions&,
        std::span<const native::ArenaRateEntry>,hall_driver::MotionBinding,const std::string&,MotionProvider&,const StepInputs&);
    std::shared_ptr<const Model> model_;const native::GapPicOptions geometry_;const MotionContext motion_;
    axial::State field_;ChargeLedger charge_;std::unique_ptr<ThermalV3> thermal_;KineticLedgers ledgers_;std::string ledgers_sha256_;
    const std::shared_ptr<const OwnerLifetime> lifetime_=std::make_shared<const OwnerLifetime>();
    std::uint64_t generation_=0;bool vacuum_=false;std::string field_sha256_,heat_sha256_;
    // v3.1: generation whose accepted state was admitted and whose field/heat digests were recomputed (by a
    // preflight, or as the candidate of the published proposal); the accepted state changes only at publication
    // or restore (new owner).
    mutable std::uint64_t accepted_verified_generation_=UINT64_MAX;
    struct RestoreTag {bool after_local_reaction=false;};
    static std::shared_ptr<const Model> checked_model(std::shared_ptr<const Model>m){need(bool(m),"missing immutable model");return m;}
    // Numeric equality of every word (the 8dddec0 check compared JSON records).
    static bool same_words(const WordVec&a,const WordVec&b){if(a.size()!=b.size())return false;
        for(std::size_t i=0;i<a.size();++i)if(a[i].hi!=b[i].hi||a[i].lo!=b[i].lo)return false;return true;}
    Owner(std::shared_ptr<const Model>m,axial::State field,ChargeLedger charge,const native::GapPicOptions&o,
        std::span<const native::ArenaRateEntry>census,RestoreTag tag):model_(checked_model(std::move(m))),geometry_(o),
        motion_(with_dt_max(make_motion_context(model_->mesh(),geometry_,model_->mesh_json(),model_->limits().maximum_carriers),model_->dt_max())),field_(std::move(field)),charge_(std::move(charge)),vacuum_(census.empty()){
        const auto actual=charge_of_census(mesh(),geometry_,census,model_->limits().maximum_carriers);
        need(tag.after_local_reaction?same_net_charge(charge_.particles,actual):
            (same_words(actual.net_C,charge_.particles.net_C)&&actual.gross_C==charge_.particles.gross_C),
            "initial/restored net charge differs from actual canonical carriers");
        model_->admit(field_,charge_);admit_census(census);
        thermal_=std::make_unique<ThermalV3>(mesh(),model_->drive().start(),model_->thermal_declaration(),model_->dt_max());
        field_sha256_=model_->state_sha256(field_,charge_);heat_sha256_=heat_sha256_of(model_->material_sha256(),*thermal_);ledgers_sha256_=ledgers_.sha256();
    }
    void admit_census(std::span<const native::ArenaRateEntry>census)const{
        auto actual=charge_of_census(mesh(),geometry_,census,model_->limits().maximum_carriers);need(same_net_charge(charge_.particles,actual),"accepted charge does not match actual canonical census");
        auto ledger=charge_;ledger.particles=std::move(actual);model_->admit(field_,ledger);
    }
    std::string conductivity_sha256()const{return Model::conductivity_sha256(thermal_->conductivity());}
    FieldCandidate prepare_field(const Vec&plasma,const FreeCharge&after,const DepositSupport&support,const Vec&sigma,
                                 const std::string&sigma_sha256,double dt,const BankTransfer*bank=nullptr)const{
        return model_->prepare(field_,charge_,field_sha256_,plasma,after,support,sigma,sigma_sha256,dt,bank);
    }
    void seal_field(FieldCandidate&f)const{model_->seal(f);}
    FieldHeatProposal heat_for(FieldCandidate&&field,std::span<const SurfaceLoad>loads={},const hall_driver::MotionReceipt*motion=nullptr){
        if(!field.sealed())model_->seal(field);
        need(field.before_sha256()==field_sha256_&&field.configuration_sha256()==configuration_sha256()&&
            field.conductivity_sha256()==conductivity_sha256(),"field preparation is stale or belongs to another configuration");
        // v3.1: no second admission here: Model::prepare admitted this (private, immutable) candidate as its last
        // statement, and the final preflight admits it again.
        auto heat=thermal_->prepare(field.receipt(),field.dt_s(),field.after_sha256(),loads);
        KineticLedgers next=ledgers_;
        if(motion&&!motion->v2_identical()){const auto step=SinkLedger::of(motion->sinks());next.sinks.add(step);next.births.add(BirthLedger::of(motion->births()));
            next.short_bridge_J=next.short_bridge_J+step.boundary[3].kinetic_J;next.escaped_J=next.escaped_J+step.escaped_kinetic_J();}
        if(heat.loads){next.loads+=heat.loads;for(const auto&[k,v]:heat.channel_J){auto&w=next.load_channel_J[k];w=w+v;}
            next.wall_heat_J=next.wall_heat_J+Words(heat.wall_heat_J);next.external_J=next.external_J+Words(heat.external_J);}
        return FieldHeatProposal(std::move(field),std::move(heat),model_->material_sha256(),lifetime_,generation_,heat_sha256_,std::move(next));
    }
    bool quick_preflight(const FieldHeatProposal&p)const noexcept{
        return p.live_&&p.owner_.lock()==lifetime_&&p.generation_==generation_&&thermal_->preflight(p.heat_);
    }
    void publish(FieldHeatProposal&&p)noexcept{
        static_assert(std::is_nothrow_swappable_v<axial::State>&&std::is_nothrow_swappable_v<ChargeLedger>&&
                      std::is_nothrow_swappable_v<std::string>&&std::is_nothrow_swappable_v<KineticLedgers>);
        if(!p.final_approved_||!quick_preflight(p)||generation_==UINT64_MAX)std::terminate();
        if(!thermal_->commit(std::move(p.heat_)))std::terminate();using std::swap;
        swap(field_,p.field_.field_);swap(charge_,p.field_.charge_);
        swap(field_sha256_,p.field_.after_);swap(heat_sha256_,p.after_heat_);
        swap(ledgers_,p.ledgers_);swap(ledgers_sha256_,p.ledgers_sha256_);
        ++generation_;p.live_=false;p.final_approved_=false;
        // The published candidate was admitted and its field/heat digests recomputed by the final preflight of
        // this proposal; the swap does not change its bytes, so the new accepted state starts verified.
        accepted_verified_generation_=generation_;
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
    const ThermalV3&thermal()const noexcept{return *thermal_;}
    const KineticLedgers&kinetic_ledgers()const noexcept{return ledgers_;}
    const std::string&kinetic_ledgers_sha256()const noexcept{return ledgers_sha256_;}
    FieldHeatProposal prepare_vacuum(double dt,std::span<const SurfaceLoad>loads={}){need(vacuum_,"field-only publication requires an empty actual census");clock_guard();
        const auto sigma=thermal_->conductivity();
        auto field=prepare_field(Vec(mesh().ne),charge_.particles,DepositSupport::empty(mesh()),sigma,Model::conductivity_sha256(sigma),dt);
        return heat_for(std::move(field),loads);}
    bool preflight_field_heat(const FieldHeatProposal&p)const{
        try{if(!quick_preflight(p)||p.field_.before_sha256()!=field_sha256_||p.before_heat_!=heat_sha256_||
                p.configuration_sha256()!=configuration_sha256()||p.field_.conductivity_sha256()!=conductivity_sha256())return false;
            clock_guard();
            // One independent recomputation of the accepted digests and admission per owner generation (v3.1: cached
            // until the next publication; the values and the refusal on a mismatch are unchanged), and of the
            // candidate digests on every call.
            if(accepted_verified_generation_!=generation_){model_->admit(field_,charge_);
                if(model_->state_sha256(field_,charge_)!=field_sha256_||heat_sha256_of(model_->material_sha256(),*thermal_)!=heat_sha256_)return false;
                accepted_verified_generation_=generation_;}
            model_->admit(p.candidate_field(),p.candidate_charge());
            if(p.candidate_field_sha256()!=model_->state_sha256(p.candidate_field(),p.candidate_charge())||
                p.candidate_heat_sha256()!=heat_sha256_of(model_->material_sha256(),p.heat_)||
                p.end_s()!=time_s()+p.field_.dt_s()||p.endpoint_revision()!=revision()+1||p.midpoint().E()!=p.field_.receipt().mid_E||
                p.ledgers_sha256_!=p.ledgers_.sha256())return false;
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
} // namespace hall_corrective_v3
