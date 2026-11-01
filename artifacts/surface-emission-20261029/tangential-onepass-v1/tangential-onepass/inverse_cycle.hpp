#pragma once
#include "tangential_joint.hpp"
#include "inward_layer_adapter.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace joint_cycle_test {
namespace se=surface_emission;
inline double dot(const se::Vector& a,const se::Vector& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline se::Vector add(const se::Vector& a,const se::Vector& b,double s=1){return {a[0]+s*b[0],a[1]+s*b[1],a[2]+s*b[2]};}
inline se::Vector mul(se::Vector a,double s){for(double& x:a)x*=s;return a;}
inline se::Vector cross(const se::Vector& a,const se::Vector& b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
inline double norm(const se::Vector& a){return std::sqrt(dot(a,a));}
inline double gamma(const se::Vector& u){return std::sqrt(1+dot(u,u)/(se::c_light*se::c_light));}
inline void need(bool b,const char* s){if(!b)throw std::runtime_error(s);}
inline void plus(se::LayerExchange& a,const se::LayerExchange& b){a.physical_electron_number_change+=b.physical_electron_number_change;a.particle_charge_change_C+=b.particle_charge_change_C;a.electrode_charge_change_C+=b.electrode_charge_change_C;a.wall_heat_J+=b.wall_heat_J;a.particle_K_change_J+=b.particle_K_change_J;a.electron_reference_change_J+=b.electron_reference_change_J;a.image_energy_change_J+=b.image_energy_change_J;a.field_work_J+=b.field_work_J;}
struct Endpoint {se::Vector r,u;};
inline Endpoint ordinary_midpoint(const Endpoint& a,const se::Stage& f,double h) {
    auto u=a.u;const double g0=gamma(a.u),qm=-se::qe/se::me;
    for(int i=0;i<60;++i){const auto k=mul(f.total_B_T,qm*h/(g0+gamma(u))),rhs=add(add(a.u,mul(f.total_E,qm*h)),cross(a.u,k));const auto next=mul(add(add(rhs,cross(rhs,k)),mul(k,dot(rhs,k))),1/(1+dot(k,k)));const double error=norm(add(next,u,-1));u=next;if(error<2e-13*std::max(1.,norm(u)))return {add(a.r,mul(add(u,a.u),h/(g0+gamma(u)))),u};}
    throw std::runtime_error("Ordinary midpoint cycle did not converge");
}
struct Receipt {se::LayerExchange total;int handoffs=0,reentries=0,epochs=0,quantum_reflections=0,quantum_absorptions=0;bool completed=false;double transported_scale_J=0,return_time_s=0,max_reconstructed_phase_error_eV=0;};
inline Receipt cycle(se::LayerTrial birth,se::Stage stage,const se::LayerPolicy& policy,double common_dt,const tangential_onepass::Config& config) {
    Receipt receipt;receipt.total=birth.exchange;receipt.transported_scale_J=std::abs(birth.exchange.particle_K_change_J)+std::abs(birth.exchange.wall_heat_J)+std::abs(birth.exchange.electron_reference_change_J)+std::abs(birth.exchange.image_energy_change_J)+std::abs(birth.exchange.field_work_J);
    auto state=birth.state;double accepted_time=stage.t;stage.dt=common_dt;
    while(state.status!=se::LayerStatus::Recollected && accepted_time<4e-12) {
        need(++receipt.epochs<=2000,"Bounded inverse cycle exceeded epochs");
        stage.t=stage.stage_time=stage.field_time=stage.thermal_time=accepted_time;stage.accepted_state=stage.field_accepted_state=stage.thermal_accepted_state=receipt.epochs;
        const double end=accepted_time+common_dt;int transitions=0;
        while(state.status!=se::LayerStatus::Recollected && state.current_time_s<end) {
            need(++transitions<=8,"Bounded inverse cycle interface chatter");
            if(state.status==se::LayerStatus::Active) {
                auto trial=se::advance_layer(state,stage,policy,end-state.current_time_s);
                if(trial.state.status==se::LayerStatus::Recollected) {
                    const double draw=receipt.quantum_reflections==0?.999999999:1e-12;
                    const auto resolved=tangential_onepass::resolve_quantum_return(stage,config,policy,trial,draw);
                    receipt.max_reconstructed_phase_error_eV=std::max({receipt.max_reconstructed_phase_error_eV,std::abs(resolved.reconstructed_normal_eV-birth.state.particle.normal_metal_energy_eV),std::abs(resolved.reconstructed_transverse_eV-birth.state.particle.transverse_energy_eV)});
                    if(resolved.absorbed)++receipt.quantum_absorptions;else ++receipt.quantum_reflections;
                    trial=resolved.trial;
                }
                state=trial.state;plus(receipt.total,trial.exchange);
                receipt.transported_scale_J+=std::abs(trial.exchange.particle_K_change_J)+std::abs(trial.exchange.image_energy_change_J)+std::abs(trial.exchange.field_work_J)+std::abs(trial.exchange.wall_heat_J);
                if(state.status==se::LayerStatus::PICHandedOff)++receipt.handoffs;
            }
            if(state.status==se::LayerStatus::PICHandedOff && state.current_time_s<end) {
                const Endpoint a{state.position_from_surface_m,state.particle.proper_velocity_m_s};double h=end-state.current_time_s;auto b=ordinary_midpoint(a,stage,h);
                const bool inward=dot(b.r,state.outward_normal)<=state.handoff_distance_m && dot(b.u,state.outward_normal)<0;
                if(inward){double lo=0,hi=h;for(int i=0;i<48;++i){const double mid=(lo+hi)/2;const auto candidate=ordinary_midpoint(a,stage,mid);if(dot(candidate.r,state.outward_normal)>state.handoff_distance_m)lo=mid;else hi=mid;}h=(lo+hi)/2;b=ordinary_midpoint(a,stage,h);}
                const double K0=se::kinetic_energy_J(a.u),K1=se::kinetic_energy_J(b.u),work=-se::qe*dot(stage.total_E,add(b.r,a.r,-1)),weight=state.particle.weight;
                need(std::abs(K1-K0-work)<1e-10*std::max({K0,K1,1e-25}),"Ordinary inverse-cycle current/work mismatch");
                plus(receipt.total,{0,0,0,0,weight*(K1-K0),0,0,weight*work});receipt.transported_scale_J+=weight*(std::abs(K1-K0)+std::abs(work));
                state.position_from_surface_m=b.r;state.particle.proper_velocity_m_s=b.u;state.current_time_s+=h;
                if(inward){const auto trial=surface_emission_integration::reenter_image_layer(state,stage,policy);state=trial.state;plus(receipt.total,trial.exchange);++receipt.reentries;}
            }
        }
        if(state.status!=se::LayerStatus::Recollected) {
            const auto& n=state.outward_normal;const double x=dot(state.position_from_surface_m,n),B=norm(stage.total_B_T),F=se::extraction_field(stage.total_E,n);
            const auto direction=mul(cross(n,stage.total_B_T),1/B);auto metal=state.particle.proper_velocity_m_s;const double un=dot(metal,n);
            for(int j=0;j<3;++j)metal[j]+=-un*n[j]+se::qe*B*x/se::me*direction[j];
            const double Et=se::me*dot(metal,metal)/(2*se::qe),image=state.status==se::LayerStatus::Active?-se::qe/(16*3.14159265358979323846*se::eps0*x):0;
            const double E=se::kinetic_energy_J(state.particle.proper_velocity_m_s)/se::qe+state.work_function_eV+image-F*x-Et;
            receipt.max_reconstructed_phase_error_eV=std::max({receipt.max_reconstructed_phase_error_eV,std::abs(E-birth.state.particle.normal_metal_energy_eV),std::abs(Et-birth.state.particle.transverse_energy_eV)});
            need(receipt.max_reconstructed_phase_error_eV<1e-5,"Canonical phase drift across the classical image/PIC interfaces");
        }
        accepted_time=end;
    }
    receipt.completed=state.status==se::LayerStatus::Recollected;
    const double remaining=receipt.completed?0:state.particle.weight;
    need(std::abs(receipt.total.physical_electron_number_change-remaining)<1e-12 && std::abs(receipt.total.particle_charge_change_C+se::qe*remaining)<1e-30 && std::abs(receipt.total.electrode_charge_change_C-se::qe*remaining)<1e-30,"Joint inverse cycle leaked or duplicated its retained/collected particle");
    need(std::abs(receipt.total.residual_J())<1e-10*receipt.transported_scale_J,"Joint inverse-cycle complete reference/work/heat mismatch");
    if(receipt.completed)need(std::abs(receipt.total.wall_heat_J)<1e-8*receipt.transported_scale_J && receipt.quantum_absorptions==1,"Static closed source-return cycle created material heat or duplicate absorption");
    else need(receipt.quantum_absorptions==0,"An absorbed packet remained in the inventory");
    need(receipt.max_reconstructed_phase_error_eV<1e-5,"Reciprocal return failed canonical phase reconstruction");
    receipt.return_time_s=state.current_time_s;return receipt;
}
}
