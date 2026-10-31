#pragma once
// Hall corrective field model, v3.2: unchanged v3/v3.1 arithmetic and guards except the
// absolute Gauss admission, which is measured relative to the admitted Gauss ledger carried
// in the charge ledger (field_contract.hpp, admitted_gauss_ledger_rule). prepare accumulates
// the ledger from the per-node values its own Gauss increment and endpoint continuity guards
// admitted; the rule text enters the configuration identity.
// Hall corrective field model, v3 (derived from v2, itself derived from
// hall-field-correction-20261030/field_model.hpp, commit 8dddec0). Unchanged:
// the implicit edge operator, the compensated right-hand side, the update, the
// authoritative compensated midpoint, current/Joule/energy receipts and every
// guard bound of 8dddec0, and the v2 changes below.
// v3 changes:
//  a. kinetic sinks and surface births move charge between the free charge and
//     the conductor/exterior banks inside the same field step (BankTransfer):
//     the endpoint continuity guard closes after + bank - before = div(J) dt and
//     the banks of the candidate ledger receive the transfer, so absolute Gauss
//     and the charge ledger close; the exterior bank may hold escaped charge on
//     the kinetic escape support (exterior_support());
//  b. the thermal owner is ThermalV3 (external surface loads, typed refusals,
//     optional heat-only source patch declaration bound into the material and
//     configuration identities).
// v2 changes:
//  1. identities are canonical binary encodings, each hashed once per immutable
//     object (configuration at construction; candidate state when sealed);
//  2. the field solve goes through an explicit backend (solver_backend.hpp)
//     whose descriptor is bound into the configuration identity, followed by
//     the mandatory independent FP64 residual check;
//  3. the endpoint charge/current guard has an absolute deposit-rounding term
//     derived from the actual carriers (DepositSupport, census_charge.hpp), so
//     charge-conserving motion next to a node plane is admitted;
//  4. prepare is private: it is reached only through the Owner and
//     prepare_closed_proposal, which supply the owner's accepted digest.
#include "field_contract.hpp"
#include "motion_receipt.hpp"   // v3.2: motion_work_rule (configuration identity)
#include "solver_backend.hpp"
#include "thermal_v3.hpp"

namespace hall_corrective_v3 {
class Model;class Owner;class ResidentOwner;class DepositSupportBuilder;struct ModelTestAccess;
// Per-node absolute bound (coulombs) on the rounding of the actual carriers'
// point deposits and path current deposits. Constructible only from actual
// carriers (DepositSupportBuilder) or as the empty-census zero bound.
class DepositSupport {
    friend class DepositSupportBuilder;
    Vec absolute_C_;std::size_t carriers_=0;
    DepositSupport()=default;
public:
    static DepositSupport empty(const Mesh&g){DepositSupport s;s.absolute_C_.assign(g.nn,0.);return s;}
    const Vec& absolute_C()const noexcept{return absolute_C_;}
    std::size_t carriers()const noexcept{return carriers_;}
};
// Text of the endpoint continuity tolerance; part of the configuration identity.
inline constexpr char endpoint_tolerance_rule[]=
    "endpoint-continuity-v2: |after-before-div(J)dt|_n <= 1e-12*(gross_before+gross_after+|div(J)dt|)_n "
    "+ eps*(sum_{p in S(n)} Q_p*(2C+14+(16C+67)R_p) + 24*M_n*sum_{p in S(n)} Q_p*R_p); "
    "C=nx+ny+nz, Q_p=|w q|, R_p=sum_d |dx_d|/h_d, S(n)= carriers whose start/end/path cell range "
    "(widened by eps*(3c_d+8) cells) contains node n, M_n=|S(n)|";
// v3 rules (part of the configuration identity).
inline constexpr char sink_bank_rule[]=
    "kinetic-sink-banks-v2: the represented step of a carrier is the straight segment from its start x to x+D, D the "
    "full sealed work-kick displacement over its flight F; the first hit is first_kinetic_hit of "
    "native-finite-target-kinetic-policy-v1 on that segment (tie order x y lower upper axial, later exact tie wins) with "
    "fraction f; the carrier ends at the hit point x+f*D with the hit axis coordinate set to the plane; flight F*f, "
    "discarded flight F-F*f; the path current is the straight path deposit from the field start to the hit point; "
    "the incident kinetic energy is K0+q*E(path)*(hit-start) with E the path gather of that deposit, the incident "
    "proper velocity has that energy along D; the macro-charge moves from free charge to the conductor bank (surfaces "
    "0 lower foil, 1 upper foil, 3 short bridge) or the exterior bank (2 outer xy, 4 axial) by the Words CIC deposit at "
    "the hit point in the same field step; CIC weight on a node outside the bank support is residue (<= 1e-12 of the "
    "macro-charge per deposit, bounded per node); collected foil kinetic energy is wall heat of the thermal owner, "
    "short-bridge and escaped kinetic energy are recorded; a neutral carrier drifts straight with u/gamma and ends "
    "at its first hit with the same record (no charge transfer)";
inline constexpr char birth_bank_rule[]=
    "newborns-v1: a newborn (native census row with birth time tau in the step) moves only over [tau, end], flight "
    "F=end-tau, by the sealed work-kick (charged) or the straight drift (neutral) and may hit within F; it starts at a "
    "strictly interior vacuum point or exactly on a foil face with velocity into the vacuum; a charged newborn on a "
    "foil face takes its macro-charge from that foil's conductor bank at its start (Words CIC deposit, same residue "
    "rule) in the same field step; a charged interior newborn is refused; its birth kinetic energy fl(w K) is "
    "supplied by the source and recorded per origin";
inline constexpr char exterior_support_rule[]=
    "exterior-support-v1: exterior charge only on EM exterior nodes and on the kinetic target-end node planes k=kshort, k=kfeed";
class FieldCandidate {
    friend class Model;
    friend class Owner;
    axial::State field_;ChargeLedger charge_;Midpoint midpoint_{{},{}};axial::Receipt receipt_;
    Constraints constraints_;std::string before_,after_,configuration_,conductivity_;
    Vec current_;double dt_=0;SolveRecord solve_;double continuity_ratio_=0;bool sealed_=false;bool bank_transfer_=false;
    FieldCandidate()=default;
public:
    FieldCandidate(const FieldCandidate&)=delete;FieldCandidate&operator=(const FieldCandidate&)=delete;
    FieldCandidate(FieldCandidate&&)noexcept=default;FieldCandidate&operator=(FieldCandidate&&)noexcept=default;
    const axial::State&field()const noexcept{return field_;}const ChargeLedger&charge()const noexcept{return charge_;}
    const Midpoint&midpoint()const noexcept{return midpoint_;}const axial::Receipt&receipt()const noexcept{return receipt_;}
    const Constraints&constraints()const noexcept{return constraints_;}const Vec&current()const noexcept{return current_;}
    const std::string&before_sha256()const noexcept{return before_;}const std::string&after_sha256()const noexcept{return after_;}
    const std::string&configuration_sha256()const noexcept{return configuration_;}
    const std::string&conductivity_sha256()const noexcept{return conductivity_;}double dt_s()const noexcept{return dt_;}
    const SolveRecord&solve()const noexcept{return solve_;}
    // Largest |defect|/tolerance of the endpoint continuity guard (diagnostic).
    double continuity_ratio()const noexcept{return continuity_ratio_;}
    bool bank_transfer()const noexcept{return bank_transfer_;}
    bool sealed()const noexcept{return sealed_;}
};
// Configured Hall carrier cap with a host memory budget. The per-carrier model
// covers the census row and index, endpoints, paths, deposit support, device
// staging and the JSON motion record built for native admission (the largest
// term). Bound into the configuration identity.
struct CarrierLimits {
    std::size_t maximum_carriers=100000;
    std::uint64_t memory_budget_bytes=1ull<<30;
    static constexpr std::uint64_t bytes_per_carrier=4096,fixed_bytes=256ull<<20,absolute_cap=std::uint64_t(1)<<27;
};
// v3.2: motion forms of the owner's motion receipts. native_v7: census v2, exact kinetic/work sums and the motion
// record v7 (native v7 API-007 sections 1-3, HALL-LINK.txt) with the resident-neutral interval bound of the native
// v7 configuration; otherwise the v3/v3.1 forms (native candidate003/v5/v6). Bound into the configuration identity.
struct MotionForms {bool native_v7=false;double neutral_stream_dt_max_s=1e-14;};
class Model {
    friend class Owner;friend class ResidentOwner;friend struct ModelTestAccess;
    const Mesh mesh_;const CurrentDrive drive_;Vec initial_sigma_;axial::Initial initial_;std::string configuration_;double initial_gauss_solve_bound_C_=0;
    std::shared_ptr<FieldSolveBackend> solver_;CarrierLimits limits_;std::string mesh_json_,material_sha256_;
    ThermalDeclaration thermal_declaration_;std::vector<bool> exterior_support_;
    // v3.1: configured upper bound of the field step (default 1e-14 s = the v3 literal); bound into
    // the configuration identity only when it differs from the default (default identities unchanged).
    double dt_max_=1e-14;
    MotionForms forms_{};   // v3.2
    static MotionForms checked_forms(MotionForms f){need(std::isfinite(f.neutral_stream_dt_max_s)&&f.neutral_stream_dt_max_s>0,
        "native v7 neutral stream bound must be finite and positive");return f;}
    static double checked_dt_max(double v){need(std::isfinite(v)&&v>0&&v<=1e-9,"Hall dt upper bound outside (0, 1 ns]");return v;}
    static CarrierLimits checked_limits(CarrierLimits l){
        need(l.maximum_carriers>0&&l.maximum_carriers<=CarrierLimits::absolute_cap,"Hall carrier cap outside 1..2^27");
        need((long double)l.maximum_carriers*CarrierLimits::bytes_per_carrier+CarrierLimits::fixed_bytes<=(long double)l.memory_budget_bytes,
             "Hall carrier limits exceed the memory budget");
        return l;}
    static Mesh checked_mesh(MeshSpec s){
        need(s.transverse>=9&&s.transverse<=500&&s.axial>=2&&s.air>=1&&std::isfinite(s.width_m)&&s.width_m>.0054,
            "invalid immutable mesh specification");
        const std::uint64_t n=s.transverse,z=std::uint64_t(s.axial)+2*std::uint64_t(s.air);
        need(2*n*(n+1)*(z+1)+(n+1)*(n+1)*z<=std::uint64_t(INT_MAX),"mesh index product overflow");
        return Mesh(s.transverse,s.axial,s.width_m,s.air);
    }
    static std::shared_ptr<FieldSolveBackend> checked_solver(std::shared_ptr<FieldSolveBackend>s){need(bool(s),"field solve backend required");return s;}
    Vec curl_words(const WordVec&E)const{Vec B(mesh_.nb);for(int b=0;b<mesh_.nb;++b){Words sum;
        for(int k=0;k<4;++k)sum=sum+E[mesh_.C[b].index[k]]*mesh_.C[b].coefficient[k];B[b]=sum.rounded();}return B;}
    Vec ampere_words(const WordVec&B)const{WordVec result(mesh_.ne);for(int b=0;b<mesh_.nb;++b){const auto HB=B[b]*mesh_.H[b];
        for(int k=0;k<4;++k){int e=mesh_.C[b].index[k];result[e]=result[e]+HB*mesh_.C[b].coefficient[k];}}
        Vec out(mesh_.ne);for(int e=0;e<mesh_.ne;++e)out[e]=result[e].rounded();return out;}
    // Private preparation; before_sha256 is the owner's cached accepted digest.
    FieldCandidate prepare(const axial::State&old,const ChargeLedger&before,const std::string&before_sha256,const Vec&plasma,
        const FreeCharge&after,const DepositSupport&support,const Vec&sigma,const std::string&sigma_sha256,double dt,
        const BankTransfer*bank=nullptr)const{
        using clock=std::chrono::steady_clock;
        const auto&g=mesh_;admit(old,before);free_charge_guard(g,after);
        need(std::isfinite(dt)&&dt>0&&dt<=dt_max_&&old.time_s+dt>old.time_s&&old.time_s+dt<=drive_.end()&&old.revision<UINT64_MAX-1,"authoritative field interval");
        need(plasma.size()==std::size_t(g.ne)&&axial::finite(plasma),"field current shape/finiteness");
        need(support.absolute_C().size()==std::size_t(g.nn),"deposit support shape");
        FieldCandidate out;
        const bool transfer=bank&&!bank->empty;
        if(transfer)need(bank->conductor_C.size()==std::size_t(g.nn)&&bank->exterior_C.size()==std::size_t(g.nn)&&bank->gross_C.size()==std::size_t(g.nn),
            "bank transfer shape");
        out.bank_transfer_=transfer;
        auto divergence=incidence(g,plasma);
        WordVec continuity(g.nn);   // v3.2: per-node endpoint continuity defect c_n (admitted below)
        for(int n=0;n<g.nn;++n){const auto defect=transfer?after.net_C[n]+bank->conductor_C[n]+bank->exterior_C[n]-before.particles.net_C[n]-divergence[n]*dt:
                after.net_C[n]-before.particles.net_C[n]-divergence[n]*dt;
            const long double scale=before.particles.gross_C[n]+after.gross_C[n]+magnitude(divergence[n]*dt)+(transfer?(long double)bank->gross_C[n]:0.L);
            const long double tolerance=1e-12L*scale+(long double)support.absolute_C()[n];
            const long double size=std::abs(defect.wide());
            need(size<=tolerance,"actual particle endpoint charge/current mismatch");
            continuity[n]=defect;
            if(size>0)out.continuity_ratio_=std::max(out.continuity_ratio_,double(size/tolerance));}
        const double current=drive_.current(old.time_s+.5*dt);auto S=g.conduction(sigma);Vec diagonal(g.ne,1),sqrtM(g.ne);
        for(int e=0;e<g.ne;++e){sqrtM[e]=std::sqrt(g.M[e]);if(!g.fixed_E[e])diagonal[e]+=.5*dt*(S[e]+g.Y[e])/g.M[e];}
        for(int b=0;b<g.nb;++b)for(int k=0;k<4;++k){int e=g.C[b].index[k];if(!g.fixed_E[e])diagonal[e]+=.25*dt*dt*g.H[b]*std::pow(g.C[b].coefficient[k],2)/g.M[e];}
        auto apply=[&](const Vec&v){Vec E(g.ne);for(int e=0;e<g.ne;++e)if(!g.fixed_E[e])E[e]=v[e]/sqrtM[e];
            auto K=g.curl_adjoint(g.curl(E));Vec result(g.ne);for(int e=0;e<g.ne;++e)result[e]=g.fixed_E[e]?v[e]:v[e]+(.5*dt*(S[e]+g.Y[e])*E[e]+.25*dt*dt*K[e])/sqrtM[e];return result;};
        const auto oldE=words(old.E,old.E_low),oldB=words(old.B,old.B_low);auto curl=curl_words(oldE);auto estimatedB=oldB;
        for(int b=0;b<g.nb;++b)estimatedB[b]=estimatedB[b]-Words(curl[b])*(.5*dt);auto ampere=ampere_words(estimatedB);Vec rhs(g.ne);
        for(int e=0;e<g.ne;++e)if(!g.fixed_E[e])rhs[e]=((Words(ampere[e])-Words(plasma[e])-Words(current)*initial_.feed[e]
            -oldE[e]*(S[e]+g.Y[e])+Words(initial_.background_force[e]))*dt/sqrtM[e]).rounded();
        // Backend solve, then the mandatory independent FP64 residual check with
        // the sealed apply and the sealed long-double norm.
        auto&sr=out.solve_;sr.backend=solver_->name();
        auto t=clock::now();
        const Apply closure=std::cref(apply);
        const auto increment=solver_->solve(closure,rhs,200,SolveContext{g,dt,S,sqrtM,diagonal},sr);
        sr.solve_s=std::chrono::duration<double>(clock::now()-t).count();t=clock::now();
        need(increment.size()==rhs.size()&&axial::finite(increment),"field solve backend returned a misshaped or nonfinite increment");
        {const auto Ax=apply(increment);Vec r(rhs.size());for(std::size_t i=0;i<r.size();++i)r[i]=rhs[i]-Ax[i];
            sr.rhs_norm=foil::norm(rhs);sr.target=std::max(1e-30,2e-13*sr.rhs_norm);sr.independent_residual=foil::norm(r);
            if(!(sr.independent_residual<=sr.target)){std::ostringstream why;why<<std::setprecision(17)
                <<"independent FP64 field residual check failed: backend="<<sr.backend<<" residual="<<sr.independent_residual<<" target="<<sr.target;
                throw std::runtime_error(why.str());}}
        sr.check_s=std::chrono::duration<double>(clock::now()-t).count();
        out.field_=old;out.charge_=before;out.charge_.particles=after;out.dt_=dt;
        out.before_=before_sha256;out.configuration_=configuration_;out.conductivity_=sigma_sha256;out.current_=plasma;
        auto&r=out.receipt_;r.current_A=current;r.mid_E.resize(g.ne);r.solid_current.resize(g.ne);r.short_current.resize(g.ne);
        WordVec nextE(g.ne),nextB(g.nb);Vec midB(g.nb);
        for(int e=0;e<g.ne;++e){nextE[e]=oldE[e]+Words(increment[e]/sqrtM[e]);out.field_.E[e]=nextE[e].hi;out.field_.E_low[e]=nextE[e].lo;
            r.mid_E[e]=((oldE[e]+nextE[e])*.5).rounded();}
        const auto cmid=g.curl(r.mid_E);
        for(int b=0;b<g.nb;++b){nextB[b]=oldB[b]-Words(cmid[b])*dt;out.field_.B[b]=nextB[b].hi;out.field_.B_low[b]=nextB[b].lo;
            midB[b]=((oldB[b]+nextB[b])*.5).rounded();}
        out.midpoint_=Midpoint(r.mid_E,midB);ampere=g.curl_adjoint(midB);
        Vec metal_current(g.ne),exterior_current(g.ne),total(g.ne),displacement(g.ne);
        for(int e=0;e<g.ne;++e){const auto delta=nextE[e]-oldE[e];const auto true_mid=(nextE[e]+oldE[e])*.5;
            displacement[e]=(delta*(g.M[e]/dt)).rounded();r.solid_current[e]=S[e]*r.mid_E[e];
            exterior_current[e]=g.Y[e]*r.mid_E[e]-initial_.background_force[e];
            metal_current[e]=r.solid_current[e]+current*initial_.feed[e];
            total[e]=metal_current[e]+plasma[e]+exterior_current[e]+displacement[e];const double residual=ampere[e]-total[e];
            if(g.fixed_E[e]){r.short_current[e]=residual;metal_current[e]+=residual;total[e]+=residual;}
            else r.ampere_residual_A_m=std::max(r.ampere_residual_A_m,std::abs(residual));
            r.field_delta_J+=double(g.M[e]*true_mid.wide()*delta.wide());
            r.plasma_work_J+=dt*r.mid_E[e]*plasma[e];r.joule_J+=dt*r.mid_E[e]*r.solid_current[e];
            r.circuit_work_J-=dt*current*initial_.feed[e]*r.mid_E[e];r.outer_loss_J+=dt*g.Y[e]*r.mid_E[e]*r.mid_E[e];
            r.outer_supply_J+=dt*initial_.background_force[e]*r.mid_E[e];}
        for(int b=0;b<g.nb;++b)r.field_delta_J+=double(g.H[b]*((nextB[b]+oldB[b])*.5).wide()*(nextB[b]-oldB[b]).wide());
        const auto metal_delta=incidence(g,metal_current),outer_delta=incidence(g,exterior_current);
        // v3: kinetic sink / surface birth bank transfers of this step (skipped
        // entirely without transfers, so v2 ledgers are reproduced bit for bit).
        if(transfer){for(int n=0;n<g.nn;++n){
            if(bank->conductor_C[n].hi!=0||bank->conductor_C[n].lo!=0)out.charge_.conductor_C[n]=out.charge_.conductor_C[n]+bank->conductor_C[n];
            if(bank->exterior_C[n].hi!=0||bank->exterior_C[n].lo!=0)out.charge_.exterior_C[n]=out.charge_.exterior_C[n]+bank->exterior_C[n];}
            bool any=false;for(double x:bank->residue_C)any=any||x!=0;
            if(any){if(out.charge_.residue_C.empty())out.charge_.residue_C.assign(g.nn,0.);
                for(int n=0;n<g.nn;++n)if(bank->residue_C[n]!=0){const long double x=(long double)out.charge_.residue_C[n]+bank->residue_C[n];
                    double v=double(x);if((long double)v<x)v=std::nextafter(v,INFINITY);out.charge_.residue_C[n]=v;}}}
        for(int n=0;n<g.nn;++n){out.charge_.conductor_C[n]=out.charge_.conductor_C[n]+metal_delta[n]*dt;
            out.charge_.exterior_C[n]=out.charge_.exterior_C[n]+outer_delta[n]*dt;}
        r.terminal_V=-axial::dot(initial_.feed,r.mid_E);r.cell_joule_W=g.cell_joule_W(sigma,r.mid_E);
        r.balance_J=r.field_delta_J+r.plasma_work_J+r.joule_J+r.outer_loss_J-r.outer_supply_J-r.circuit_work_J;
        r.kcl_residual_A=axial::maximum(g.gradient_transpose(total));
        const auto Q0=field_charge(g,old),Q1=field_charge(g,out.field_);double activityQ=0;
        // v3.2: admitted Gauss ledger L_new = L_old + (delta - c) per node (compensated words; the same word
        // operations as the device kernel k_gauss_increment, so both ledgers are bitwise equal).
        auto&ledger=out.charge_.gauss_admitted_C;if(ledger.empty())ledger.assign(g.nn,Words());
        for(int n=0;n<g.nn;++n){const auto change=Q1[n]-Q0[n];const auto expected=(divergence[n]+metal_delta[n]+outer_delta[n])*dt;
            const auto delta=change-expected;
            r.gauss_increment_C=std::max(r.gauss_increment_C,double(std::abs(delta.wide())));activityQ=std::max(activityQ,double(magnitude(change)));
            ledger[n]=ledger[n]+(delta-continuity[n]);}
        need(r.gauss_increment_C<=1e-22+2e-10*activityQ,"existing field Gauss increment bound");
        const int cut=g.kshort+g.nz_foil/2-1;r.region_conduction_A=g.upper_region_inflow(r.solid_current,cut);
        r.region_plasma_A=g.upper_region_inflow(plasma,cut);r.region_displacement_A=g.upper_region_inflow(displacement,cut);r.region_boundary_A=g.upper_region_inflow(exterior_current,cut);
        const double activity=std::abs(r.field_delta_J)+std::abs(r.plasma_work_J)+std::abs(r.joule_J)+std::abs(r.outer_loss_J)+std::abs(r.outer_supply_J)+std::abs(r.circuit_work_J);
        need(std::isfinite(r.balance_J)&&std::abs(r.balance_J)<=1e-20+1e-10*activity,"unchanged field/current/Joule energy bound");
        out.field_.time_s=old.time_s+dt;++out.field_.revision;out.constraints_=admit(out.field_,out.charge_);
        return out;
    }
    // Hashes the accepted candidate once; prepare does not hash discarded iterates.
    void seal(FieldCandidate&c)const{need(!c.sealed_&&c.configuration_==configuration_,"field candidate already sealed or foreign");
        c.after_=state_sha256(c.field_,c.charge_);c.sealed_=true;}
public:
    Model(MeshSpec spec,const std::filesystem::path&trace,std::shared_ptr<FieldSolveBackend> solver=make_cpu_pcg(),CarrierLimits limits={},
          ThermalDeclaration thermal={},double dt_max=1e-14,MotionForms forms={})
        :mesh_(checked_mesh(spec)),drive_(trace),solver_(checked_solver(std::move(solver))),limits_(checked_limits(limits)),
         thermal_declaration_(std::move(thermal)),dt_max_(checked_dt_max(dt_max)),forms_(checked_forms(forms)){
        auto law=electrode_heat_v2::conditional_ss304_v2(ss304v2::Qualification::evaluated_solid_1500,
            ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000.,electrode_heat_v2::material_headers_sha256);
        initial_sigma_.resize(mesh_.cells.size());for(std::size_t c=0;c<initial_sigma_.size();++c)if(mesh_.cells[c].body>=0)initial_sigma_[c]=law.sigma(300.);
        initial_=axial::dc_initial(mesh_,initial_sigma_,drive_.current(drive_.start()));
        Vec known(mesh_.nn);for(int n=0;n<mesh_.nn;++n)if(mesh_.node_body[n]>=0||mesh_.short_node[n])known[n]=initial_.potential[n];
        auto flux=mesh_.gradient(known);for(int e=0;e<mesh_.ne;++e)flux[e]*=mesh_.VE[e];
        auto rhs=mesh_.gradient_transpose(flux);for(int n=0;n<mesh_.nn;++n)if(mesh_.node_body[n]>=0||mesh_.short_node[n])rhs[n]=0;
        initial_gauss_solve_bound_C_=axial::eps*std::max(1e-30,2e-13*axial::norm(rhs));
        // Identities computed once: the sealed material binding digest (same value
        // as 8dddec0 hall_driver::material_binding_sha256) and the mesh JSON text
        // used by the unchanged motion mesh identity.
        {Json j=ThermalV3::material_binding(thermal_declaration_);material_sha256_=sha256(j.dump());drain(j);}
        exterior_support_.assign(mesh_.nn,false);
        for(int n=0;n<mesh_.nn;++n){const int k=n/((mesh_.nx+1)*(mesh_.ny+1));exterior_support_[n]=mesh_.outer_node[n]||k==mesh_.kshort||k==mesh_.kfeed;}
        {Json j=mesh_.json();mesh_json_=j.dump();drain(j);}
        configuration_=configuration_digest(dt_max_,solver_->name(),solver_->descriptor_sha256());
    }
    // Configuration identity (computed once at construction). v3.2: the admitted Gauss ledger rule is bound;
    // configuration_with() recomputes it for another field-step upper bound and/or solve backend record
    // (checkpoint transfer between such configurations, ResidentOwner::restore).
    std::string configuration_digest(double dt_max,const std::string&solver_name,const std::string&solver_descriptor)const{
        return configuration_digest(dt_max,solver_name,solver_descriptor,forms_);}
    std::string configuration_digest(double dt_max,const std::string&solver_name,const std::string&solver_descriptor,const MotionForms&forms)const{
        Digest d("hall-v3-field-configuration-v1");
        d.text(api);encode_mesh(d,mesh_);
        d.text(material_sha256_);d.text(c19_sha256);d.text("impressed_total_current_only");
        d.text("literal_C19_support_no_shift_no_extrapolation");d.text("declared_300K_DC_at_first_supported_knot");
        for(const Vec*v:{&initial_sigma_,&initial_.E,&initial_.B,&initial_.potential,&initial_.feed,&initial_.background_force,&initial_.short_current})d.f64s(*v);
        d.text("hall-v2-constraints-v1");d.f64(initial_gauss_solve_bound_C_);d.u64(256);d.u64(128);
        d.text("epsilon*max(1e-30,2e-13*initial_vacuum_rhs_L2)");d.f64(0);d.text("conductor_EM_exterior_and_kinetic_escape_banks");
        d.text(sink_bank_rule);d.text(birth_bank_rule);d.text(exterior_support_rule);
        d.text("hall-v3-thermal-owner: ThermalV3 surface loads, exactly-once watermarks, typed capacity refusals");
        d.text(thermal_declaration_.empty()?std::string("no-source-patch-declaration"):thermal_declaration_.sha256());
        d.text(endpoint_tolerance_rule);
        d.text("independent FP64 residual ||b-Ax||_2<=max(1e-30,2e-13||b||_2) with the sealed apply and norm, every backend");
        d.text(solver_name);d.text(solver_descriptor);
        d.text("hall-v2-carrier-limits-v1");d.u64(limits_.maximum_carriers);d.u64(limits_.memory_budget_bytes);
        d.u64(CarrierLimits::bytes_per_carrier);d.u64(CarrierLimits::fixed_bytes);
        if(dt_max!=1e-14){d.text("hall-v3.1-dt-upper-bound-v1");d.f64(dt_max);}
        d.text(admitted_gauss_ledger_rule);d.text(motion_work_rule);
        if(forms.native_v7){d.text("native-v7-motion-forms-v1: census v2, exact kinetic/work sums, motion record v7 merkle (API-007 1.1-1.3)");
            d.text("neutral_stream_dt_max_s");d.f64(forms.neutral_stream_dt_max_s);}
        else d.text("native-v6-motion-forms: census v1, x87 kinetic/work sums, motion record v3 and v2 JSON seal");
        return d.hex();
    }
    Model(const Model&)=delete;Model&operator=(const Model&)=delete;
    std::string configuration_with(double dt_max,const std::string&solver_name,const std::string&solver_descriptor,const MotionForms&forms)const{
        return configuration_digest(checked_dt_max(dt_max),solver_name,solver_descriptor,checked_forms(forms));}
    const Mesh&mesh()const noexcept{return mesh_;}const CurrentDrive&drive()const noexcept{return drive_;}
    double initial_gauss_solve_bound_C()const noexcept{return initial_gauss_solve_bound_C_;}
    const std::string&configuration_sha256()const noexcept{return configuration_;}const Vec&initial_conductivity()const noexcept{return initial_sigma_;}
    const FieldSolveBackend&solver()const noexcept{return *solver_;}
    const CarrierLimits&limits()const noexcept{return limits_;}
    const std::string&material_sha256()const noexcept{return material_sha256_;}
    const std::string&mesh_json()const noexcept{return mesh_json_;}
    const ThermalDeclaration&thermal_declaration()const noexcept{return thermal_declaration_;}
    // v3.1: configured field-step upper bound and the immutable DC initial operands (feed, background force).
    double dt_max()const noexcept{return dt_max_;}
    const MotionForms& motion_forms()const noexcept{return forms_;}
    const axial::Initial&initial_fields()const noexcept{return initial_;}
    // Nodes where the exterior bank may hold charge (EM exterior and the kinetic target-end planes).
    const std::vector<bool>&exterior_support()const noexcept{return exterior_support_;}
    // Unchanged motion mesh identity (equals hall_driver::motion_mesh_sha256), no JSON per call.
    std::string motion_mesh_sha256(const raikou::mfem_backend::GapPicOptions&o)const;
    axial::State dc_seed()const{return {initial_.E,initial_.B,drive_.start(),0,Vec(mesh_.ne),Vec(mesh_.nb)};}
    std::string state_sha256(const axial::State&s,const ChargeLedger&q)const{Digest d("hall-v2-field-state-v1");d.text(api);
        d.text(configuration_);encode_state(d,s);encode_ledger(d,q);return d.hex();}
    static std::string conductivity_sha256(const Vec&sigma){Digest d("hall-v2-conductivity-v1");d.f64s(sigma);return d.hex();}
    Constraints admit(const axial::State&s,const ChargeLedger&q)const{need(s.time_s>=drive_.start()&&s.time_s<=drive_.end(),"field state outside C19 clock support");return admit_constraints(mesh_,s,q,initial_gauss_solve_bound_C_,&exterior_support_);}
    Midpoint at(const axial::State&s)const{const auto E=words(s.E,s.E_low),B=words(s.B,s.B_low);Vec er(E.size()),br(B.size());
        for(std::size_t i=0;i<E.size();++i)er[i]=E[i].rounded();for(std::size_t i=0;i<B.size();++i)br[i]=B[i].rounded();return Midpoint(std::move(er),std::move(br));}
};
} // namespace hall_corrective_v3
