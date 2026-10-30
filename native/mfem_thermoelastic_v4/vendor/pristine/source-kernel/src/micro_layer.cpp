#include "surface_emission.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace surface_emission {
namespace {
constexpr double pi=3.1415926535897932384626433832795,hbar=planck/(2*pi),image_c=qe/(16*pi*eps0);
void need(bool okay,const char* message){if(!okay)throw std::invalid_argument(message);}
double dot(const Vector& a,const Vector& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
Vector add(const Vector& a,const Vector& b,double scale=1){return {a[0]+scale*b[0],a[1]+scale*b[1],a[2]+scale*b[2]};}
Vector mul(Vector a,double scale){for(double& v:a)v*=scale;return a;}
Vector cross(const Vector& a,const Vector& b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double norm(const Vector& a){return std::sqrt(dot(a,a));}
bool finite(const Vector& a){return std::all_of(a.begin(),a.end(),[](double x){return std::isfinite(x);});}
double gamma(const Vector& u){return std::sqrt(1+dot(u,u)/(c_light*c_light));}
void validate(const LayerPolicy& p) {
    need(std::isfinite(p.local_linear_vacuum_extent_m) && p.local_linear_vacuum_extent_m>0 &&
         p.image_force_fraction>0 && p.image_force_fraction<=.2 && p.image_energy_fraction>0 && p.image_energy_fraction<=.2 &&
         p.maximum_quantum_magnetic_energy_fraction>0 && p.maximum_quantum_magnetic_energy_fraction<=.3 &&
         p.maximum_gyroangle>0 && p.maximum_gyroangle<=.05 && p.maximum_relative_normal_step>0 && p.maximum_relative_normal_step<=.05 &&
         p.maximum_substeps>=1 && p.maximum_substeps<=16384 && !p.PIC_already_applies_microscopic_single_electron_image_force,
         "Invalid bounded Lorentz/image-layer policy or duplicate microscopic image force");
}
Vector momentum(const Vector& n,double normal_eV,double transverse_eV,double azimuth) {
    const double K=(normal_eV+transverse_eV)*qe;
    need(normal_eV>=-1e-12 && K>0 && std::isfinite(K),"Invalid overlap kinetic energy");
    const double ut=std::sqrt(2*qe*transverse_eV/me),un=std::sqrt(std::max(0.,2*qe*normal_eV/me)+std::pow(K/(me*c_light),2));
    const Vector axis=std::abs(n[0])<.8?Vector{1,0,0}:Vector{0,1,0};auto tangent=cross(n,axis);tangent=mul(tangent,1/norm(tangent));const auto other=cross(n,tangent);
    return add(mul(n,un),mul(add(mul(tangent,std::cos(2*pi*azimuth)),mul(other,std::sin(2*pi*azimuth))),ut));
}
LayerExchange exchange(const Ledger& l){return {l.physical_electron_number,l.particle_charge_C,l.electrode_charge_C,l.wall_heat_J,l.outgoing_K_J-l.incoming_K_J,l.electron_reference_change_J,l.image_reference_change_J,l.resolved_matching_field_work_J};}
void accumulate(LayerExchange& a,const LayerExchange& b){
    a.physical_electron_number_change+=b.physical_electron_number_change;a.particle_charge_change_C+=b.particle_charge_change_C;a.electrode_charge_change_C+=b.electrode_charge_change_C;
    a.wall_heat_J+=b.wall_heat_J;a.particle_K_change_J+=b.particle_K_change_J;a.electron_reference_change_J+=b.electron_reference_change_J;a.image_energy_change_J+=b.image_energy_change_J;a.field_work_J+=b.field_work_J;
}
struct Endpoint {Vector r,u;};
bool implicit_step(const Endpoint& a,const Vector& n,const Stage& field,double h,Endpoint& out) {
    const double x0=dot(a.r,n),g0=gamma(a.u),qm=-qe/me;
    if(!(x0>0))return false;
    Vector u=a.u,r=add(a.r,mul(a.u,h/g0));
    for(int iteration=0;iteration<60;++iteration) {
        const double x1=dot(r,n);if(!(x1>0) || !finite(u) || !finite(r))return false;
        const auto E=add(field.total_E,mul(n,image_c/(x0*x1)));
        const auto k=mul(field.total_B_T,qm*h/(g0+gamma(u)));
        const auto rhs=add(add(a.u,mul(E,qm*h)),cross(a.u,k));
        const auto next=mul(add(add(rhs,cross(rhs,k)),mul(k,dot(rhs,k))),1/(1+dot(k,k)));
        const auto next_r=add(a.r,mul(add(next,a.u),h/(g0+gamma(next))));
        const double uerr=norm(add(next,u,-1)),rerr=norm(add(next_r,r,-1));
        u=next;r=next_r;
        if(uerr<=2e-13*std::max(1.,norm(u)) && rerr<=2e-13*std::max(x0,norm(r))) {out={r,u};return true;}
    }
    return false;
}
}
double LayerExchange::residual_J() const{return wall_heat_J+particle_K_change_J+electron_reference_change_J+image_energy_change_J-field_work_J;}
LayerTrial begin_layer(const PreparedFlux& f,const LayerPolicy& p,double weight,double un,double ut,double ua,std::uint64_t id) {
    validate(p);need(un>0 && un<1 && ut>0 && ut<1 && ua>0 && ua<1 && std::isfinite(weight) && weight>0 && id>0,"Invalid caller-owned micro-layer sample");
    LayerTrial trial;auto& state=trial.state;state.outward_normal=f.stage.outward_normal;state.current_time_s=f.stage.stage_time;
    state.mesh_revision=f.stage.mesh_revision;state.material_revision=f.stage.material_revision;state.work_function_eV=f.spectrum.config.work_function_eV;
    if(f.spectrum.extraction_field_V_m<=0) {
        const MatchingPlane surface{PlaneKind::ThermalVacuumSurfaceStep,0,true,true,true};state.particle=sample(f,surface,weight,un,ut,ua,id);
        state.status=LayerStatus::PICHandedOff;trial.exchange=exchange(state.particle.ledger);return trial;
    }
    need(!f.spectrum.nodes.empty(),"Empty transmitted spectrum");
    const auto it=std::lower_bound(f.spectrum.nodes.begin(),f.spectrum.nodes.end(),un,[](const Node& n,double u){return n.cdf<u;});need(it!=f.spectrum.nodes.end(),"Invalid normal CDF");
    const double E=it->normal_energy_eV,t=transverse_quantile_eV(E,f.spectrum.T_K,ut),F=f.spectrum.extraction_field_V_m,phi=state.work_function_eV;
    state.source_energy_scale_eV=std::max(kb*f.spectrum.T_K/qe,f.spectrum.fermi_decay_width_eV);need(state.source_energy_scale_eV>0,"Missing source energy scale");
    const double overlap=p.image_energy_fraction*state.source_energy_scale_eV,A=phi-E+overlap,lowering=f.spectrum.barrier_lowering_eV;
    const double x=A>lowering?(A+std::sqrt((A-lowering)*(A+lowering)))/(2*F):std::sqrt(image_c/F);
    need(x>0 && x<p.local_linear_vacuum_extent_m,"Quantum overlap exceeds actual local vacuum geometry");
    const auto& n=state.outward_normal;const double Bnormal=dot(f.stage.total_B_T,n);const auto Bt=add(f.stage.total_B_T,mul(n,-Bnormal));
    const double transverse_speed=std::sqrt(2*qe*t/me);
    const double magnetic_shift=norm(Bt)*x*transverse_speed+qe*dot(Bt,Bt)*x*x/(2*me)+.5*hbar*std::abs(Bnormal)/me;
    state.quantum_magnetic_energy_fraction=magnetic_shift/state.source_energy_scale_eV;
    need(state.quantum_magnetic_energy_fraction<=p.maximum_quantum_magnetic_energy_fraction,
         "Unmagnetized MG barrier approximation exceeds explicit magnetic energy criterion; no packet resampling");
    const double image=-image_c/x,normal=E-phi+F*x-image;
    state.position_from_surface_m=mul(n,x);state.return_overlap_distance_m=x;
    state.handoff_distance_m=std::max(std::sqrt(image_c/(p.image_force_fraction*F)),image_c/(p.image_energy_fraction*state.source_energy_scale_eV));
    const auto u=momentum(n,normal,t,ua);const double K=kinetic_energy_J(u),scale=qe*weight;
    Ledger ledger{weight,-scale,scale,-(E+t)*scale,K*weight,0,phi*scale,image*scale,F*x*scale};
    state.particle={id,f.stage.accepted_state,f.stage.patch_id,weight,E,t,u,ledger};trial.exchange=exchange(ledger);
    trial.paths.push_back({{},state.position_from_surface_m,state.current_time_s,state.current_time_s,-scale,F*x*scale});
    return trial;
}
LayerTrial advance_layer(const LayerState& initial,const Stage& field,const LayerPolicy& p,double duration) {
    validate(p);need(initial.status==LayerStatus::Active,"Layer transport already handed off or recollected");
    need(field.accepted_state>0 && field.field_accepted_state==field.accepted_state && field.field_time==field.stage_time &&
         field.mesh_revision==initial.mesh_revision && field.material_revision==initial.material_revision && field.patch_id==initial.particle.patch_id &&
         field.outward_normal==initial.outward_normal && finite(field.total_E) && finite(field.total_B_T),"Local Lorentz field/geometry revision mismatch");
    need(std::isfinite(duration) && duration>0 && std::isfinite(initial.current_time_s+duration) && initial.current_time_s>=field.t &&
         initial.current_time_s+duration<=field.t+field.dt+1e-24,"Micro-layer flight exceeds accepted field interval");
    LayerTrial trial;trial.state=initial;auto& s=trial.state;const auto& n=s.outward_normal;
    const double end=initial.current_time_s+duration,weight=s.particle.weight,charge=-qe*weight,B=norm(field.total_B_T);
    trial.paths.reserve(std::min(p.maximum_substeps,256));
    while(s.current_time_s<end) {
        need(trial.substeps<p.maximum_substeps,"Bounded micro-layer substeps exceeded; rollback/refine accepted interval");
        const Endpoint a{s.position_from_surface_m,s.particle.proper_velocity_m_s};const double x0=dot(a.r,n),g0=gamma(a.u),v=norm(a.u)/g0;
        need(x0>0 && norm(a.r)<=p.local_linear_vacuum_extent_m,"Active layer left its actual locally valid field/geometry region");
        const double acceleration=qe/me*(norm(field.total_E)+image_c/(x0*x0)+v*B);
        double h=std::min(end-s.current_time_s,std::sqrt(p.maximum_relative_normal_step*x0/std::max(acceleration,1.)));
        h=std::min(h,p.maximum_relative_normal_step*x0/std::max(std::abs(dot(a.u,n))/g0,1.));
        if(B>0)h=std::min(h,p.maximum_gyroangle*g0*me/(qe*B));
        need(h>1e-25 && s.current_time_s+h>s.current_time_s,"Nonprogressing bounded micro-layer step");
        Endpoint b;bool success=false;
        for(int retry=0;retry<24;++retry){if(implicit_step(a,n,field,h,b)){success=true;break;}h*=.5;}
        need(success,"Implicit Lorentz/image endpoint did not converge; no partial accepted commit");
        bool return_crossing=dot(b.r,n)<=s.return_overlap_distance_m && dot(b.u,n)<0;
        if(return_crossing) {
            double lo=0,hi=h;
            for(int bisect=0;bisect<48;++bisect) {
                const double mid=(lo+hi)/2;Endpoint candidate;need(implicit_step(a,n,field,mid,candidate),"Return boundary interpolation failed");
                if(dot(candidate.r,n)>s.return_overlap_distance_m)lo=mid;else hi=mid;
            }
            h=(lo+hi)/2;need(implicit_step(a,n,field,h,b),"Return endpoint failed");
        }
        const double x1=dot(b.r,n),K0=kinetic_energy_J(a.u),K1=kinetic_energy_J(b.u),U0=-qe*image_c/x0,U1=-qe*image_c/x1;
        const auto displacement=add(b.r,a.r,-1);const double work=-qe*dot(field.total_E,displacement);
        const double residual=K1-K0+U1-U0-work,energy_scale=std::max({K0,K1,std::abs(U0),std::abs(U1),qe*s.source_energy_scale_eV});
        need(std::isfinite(residual) && std::abs(residual)<=3e-11*energy_scale,"Conservative Lorentz/image discrete work guard failed");
        need(norm(b.r)<=p.local_linear_vacuum_extent_m,"Lorentz/image path leaves admitted local field region; no sheath/gap leap");
        accumulate(trial.exchange,{0,0,0,0,(K1-K0)*weight,0,(U1-U0)*weight,work*weight});
        trial.paths.push_back({a.r,b.r,s.current_time_s,s.current_time_s+h,charge,work*weight});
        s.current_time_s+=h;s.position_from_surface_m=b.r;s.particle.proper_velocity_m_s=b.u;++trial.substeps;
        if(return_crossing) {
            need(p.own_return_unit_accommodation,"Own returned-electron accommodation model not selected");
            const auto surface=add(b.r,mul(n,-x1));const double bridge=-qe*dot(field.total_E,add(surface,b.r,-1));
            const double phi=qe*s.work_function_eV;
            accumulate(trial.exchange,{-weight,-charge,charge,(K1+phi+U1+bridge)*weight,-K1*weight,-phi*weight,-U1*weight,bridge*weight});
            trial.paths.push_back({b.r,surface,s.current_time_s,s.current_time_s,charge,bridge*weight});
            s.position_from_surface_m=surface;s.status=LayerStatus::Recollected;break;
        }
        const double normal=dot(b.u,n),perp2=std::max(0.,dot(b.u,b.u)-normal*normal);
        const double target=K1+U1; // image-tail restoration: K decreases by |U_img|, image reference rises to0
        const double total_u2=target>0?target*(target+2*me*c_light*c_light)/(me*me*c_light*c_light):-1;
        const double Fnormal=std::abs(dot(field.total_E,n));
        if(normal>0 && x1>=s.handoff_distance_m && Fnormal>0 && image_c/(x1*x1)<=p.image_force_fraction*Fnormal &&
           image_c/x1<=p.image_energy_fraction*s.source_energy_scale_eV && total_u2>=perp2) {
            s.particle.proper_velocity_m_s=add(b.u,mul(n,std::sqrt(total_u2-perp2)-normal));
            accumulate(trial.exchange,{0,0,0,0,U1*weight,0,-U1*weight,0});s.status=LayerStatus::PICHandedOff;break;
        }
    }
    const double scale=std::abs(trial.exchange.particle_K_change_J)+std::abs(trial.exchange.image_energy_change_J)+std::abs(trial.exchange.field_work_J)+std::abs(trial.exchange.electron_reference_change_J)+std::abs(trial.exchange.wall_heat_J);
    need(std::isfinite(trial.exchange.residual_J()) && std::abs(trial.exchange.residual_J())<=1e-10*std::max(scale,qe*weight*s.source_energy_scale_eV),"Micro-layer aggregate energy closure failed");
    return trial;
}
}
