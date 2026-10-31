#pragma once
// Hall corrective field model, v2. Derived from hall-field-correction-20261030/
// field_model.hpp (commit 8dddec0). Unchanged: the implicit edge operator, the
// compensated right-hand side, the update, the authoritative compensated
// midpoint, current/Joule/energy receipts and every guard bound of 8dddec0.
// Changes:
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
#include "solver_backend.hpp"
#include "hall-canonical-driver-20261029/axial_thermal_v2.hpp"

namespace hall_corrective_v2 {
class Model;class Owner;class DepositSupportBuilder;struct ModelTestAccess;
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
class FieldCandidate {
    friend class Model;
    friend class Owner;
    axial::State field_;ChargeLedger charge_;Midpoint midpoint_{{},{}};axial::Receipt receipt_;
    Constraints constraints_;std::string before_,after_,configuration_,conductivity_;
    Vec current_;double dt_=0;SolveRecord solve_;double continuity_ratio_=0;bool sealed_=false;
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
class Model {
    friend class Owner;friend struct ModelTestAccess;
    const Mesh mesh_;const CurrentDrive drive_;Vec initial_sigma_;axial::Initial initial_;std::string configuration_;double initial_gauss_solve_bound_C_=0;
    std::shared_ptr<FieldSolveBackend> solver_;CarrierLimits limits_;std::string mesh_json_,material_sha256_;
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
        const FreeCharge&after,const DepositSupport&support,const Vec&sigma,const std::string&sigma_sha256,double dt)const{
        using clock=std::chrono::steady_clock;
        const auto&g=mesh_;admit(old,before);free_charge_guard(g,after);
        need(std::isfinite(dt)&&dt>0&&dt<=1e-14&&old.time_s+dt>old.time_s&&old.time_s+dt<=drive_.end()&&old.revision<UINT64_MAX-1,"authoritative field interval");
        need(plasma.size()==std::size_t(g.ne)&&axial::finite(plasma),"field current shape/finiteness");
        need(support.absolute_C().size()==std::size_t(g.nn),"deposit support shape");
        FieldCandidate out;
        auto divergence=incidence(g,plasma);
        for(int n=0;n<g.nn;++n){const auto defect=after.net_C[n]-before.particles.net_C[n]-divergence[n]*dt;
            const long double scale=before.particles.gross_C[n]+after.gross_C[n]+magnitude(divergence[n]*dt);
            const long double tolerance=1e-12L*scale+(long double)support.absolute_C()[n];
            const long double size=std::abs(defect.wide());
            need(size<=tolerance,"actual particle endpoint charge/current mismatch");
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
        for(int n=0;n<g.nn;++n){out.charge_.conductor_C[n]=out.charge_.conductor_C[n]+metal_delta[n]*dt;
            out.charge_.exterior_C[n]=out.charge_.exterior_C[n]+outer_delta[n]*dt;}
        r.terminal_V=-axial::dot(initial_.feed,r.mid_E);r.cell_joule_W=g.cell_joule_W(sigma,r.mid_E);
        r.balance_J=r.field_delta_J+r.plasma_work_J+r.joule_J+r.outer_loss_J-r.outer_supply_J-r.circuit_work_J;
        r.kcl_residual_A=axial::maximum(g.gradient_transpose(total));
        const auto Q0=field_charge(g,old),Q1=field_charge(g,out.field_);double activityQ=0;
        for(int n=0;n<g.nn;++n){const auto change=Q1[n]-Q0[n];const auto expected=(divergence[n]+metal_delta[n]+outer_delta[n])*dt;
            r.gauss_increment_C=std::max(r.gauss_increment_C,double(std::abs((change-expected).wide())));activityQ=std::max(activityQ,double(magnitude(change)));}
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
    Model(MeshSpec spec,const std::filesystem::path&trace,std::shared_ptr<FieldSolveBackend> solver=make_cpu_pcg(),CarrierLimits limits={})
        :mesh_(checked_mesh(spec)),drive_(trace),solver_(checked_solver(std::move(solver))),limits_(checked_limits(limits)){
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
        {Json j=hall_driver::Thermal::material_binding();material_sha256_=sha256(j.dump());drain(j);}
        {Json j=mesh_.json();mesh_json_=j.dump();drain(j);}
        Digest d("hall-v2-field-configuration-v1");d.text(api);encode_mesh(d,mesh_);
        d.text(material_sha256_);d.text(c19_sha256);d.text("impressed_total_current_only");
        d.text("literal_C19_support_no_shift_no_extrapolation");d.text("declared_300K_DC_at_first_supported_knot");
        for(const Vec*v:{&initial_sigma_,&initial_.E,&initial_.B,&initial_.potential,&initial_.feed,&initial_.background_force,&initial_.short_current})d.f64s(*v);
        d.text("hall-v2-constraints-v1");d.f64(initial_gauss_solve_bound_C_);d.u64(256);d.u64(128);
        d.text("epsilon*max(1e-30,2e-13*initial_vacuum_rhs_L2)");d.f64(0);d.text("conductor_and_EM_exterior_no_kinetic_sinks");
        d.text(endpoint_tolerance_rule);
        d.text("independent FP64 residual ||b-Ax||_2<=max(1e-30,2e-13||b||_2) with the sealed apply and norm, every backend");
        d.text(solver_->name());d.text(solver_->descriptor_sha256());
        d.text("hall-v2-carrier-limits-v1");d.u64(limits_.maximum_carriers);d.u64(limits_.memory_budget_bytes);
        d.u64(CarrierLimits::bytes_per_carrier);d.u64(CarrierLimits::fixed_bytes);
        configuration_=d.hex();
    }
    Model(const Model&)=delete;Model&operator=(const Model&)=delete;
    const Mesh&mesh()const noexcept{return mesh_;}const CurrentDrive&drive()const noexcept{return drive_;}
    double initial_gauss_solve_bound_C()const noexcept{return initial_gauss_solve_bound_C_;}
    const std::string&configuration_sha256()const noexcept{return configuration_;}const Vec&initial_conductivity()const noexcept{return initial_sigma_;}
    const FieldSolveBackend&solver()const noexcept{return *solver_;}
    const CarrierLimits&limits()const noexcept{return limits_;}
    const std::string&material_sha256()const noexcept{return material_sha256_;}
    const std::string&mesh_json()const noexcept{return mesh_json_;}
    // Unchanged motion mesh identity (equals hall_driver::motion_mesh_sha256), no JSON per call.
    std::string motion_mesh_sha256(const raikou::mfem_backend::GapPicOptions&o)const;
    axial::State dc_seed()const{return {initial_.E,initial_.B,drive_.start(),0,Vec(mesh_.ne),Vec(mesh_.nb)};}
    std::string state_sha256(const axial::State&s,const ChargeLedger&q)const{Digest d("hall-v2-field-state-v1");d.text(api);
        d.text(configuration_);encode_state(d,s);encode_ledger(d,q);return d.hex();}
    static std::string conductivity_sha256(const Vec&sigma){Digest d("hall-v2-conductivity-v1");d.f64s(sigma);return d.hex();}
    Constraints admit(const axial::State&s,const ChargeLedger&q)const{need(s.time_s>=drive_.start()&&s.time_s<=drive_.end(),"field state outside C19 clock support");return admit_constraints(mesh_,s,q,initial_gauss_solve_bound_C_);}
    Midpoint at(const axial::State&s)const{const auto E=words(s.E,s.E_low),B=words(s.B,s.B_low);Vec er(E.size()),br(B.size());
        for(std::size_t i=0;i<E.size();++i)er[i]=E[i].rounded();for(std::size_t i=0;i<B.size();++i)br[i]=B[i].rounded();return Midpoint(std::move(er),std::move(br));}
};
} // namespace hall_corrective_v2
