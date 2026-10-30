#include "tangential_joint.hpp"
#include "inverse_cycle.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sched.h>
#include <sys/resource.h>

namespace se=surface_emission;namespace tj=tangential_onepass;
namespace {
void need(bool b,const char* s){if(!b)throw std::runtime_error(s);}
se::Stage stage(double T,double B) {
    se::Stage s;s.accepted_state=s.field_accepted_state=s.thermal_accepted_state=s.mesh_revision=s.material_revision=s.patch_id=1;
    s.T_K=T;s.total_E={0,-97724679.40458213,0};s.total_B_T={B,0,0};s.outward_normal={0,1,0};s.area_m2=1.25e-4;s.dt=3e-12;
    s.material_temperature_phase_admitted=s.locally_planar_barrier_admitted=s.local_equilibrium_electron_supply_admitted=true;return s;
}
}
int main(int argc,char** argv) {
    try {
        need(argc==5,"Expected OUT T B refinement");cpu_set_t mask;CPU_ZERO(&mask);need(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(7,&mask) && getpriority(PRIO_PROCESS,0)>=15,"CPU7/nice15 required");
        const std::filesystem::path out=argv[1];std::filesystem::create_directories(out);const double T=std::stod(argv[2]),B=std::stod(argv[3]);const int refinement=std::stoi(argv[4]);
        auto s=stage(T,B);tj::Config config{4.3};config.continuum_free_electron_reservoir_selected=config.incoherent_onepass_matching_selected=true;
        config.transverse_nodes=16*refinement;config.azimuth_nodes=16*refinement;config.normal_panels=4*refinement;config.action_nodes=32*refinement;
        se::LayerPolicy policy;policy.local_linear_vacuum_extent_m=25e-6;
        const auto d=tj::evaluate(s,config,policy);need(d.truncation_bounds_within_declared_budget,"Unresolved phase/tail bound too large for this tangential benchmark");
        const auto nominal=se::prepare(s,se::Config{4.3});
        if(B==0){need(std::abs(d.known_component_current_A_m2/nominal.spectrum.current_A_m2-1)<2e-5,"B0 joint integral failed immutable MG current comparison");need(std::abs(d.known_component_mean_metal_energy_eV-nominal.spectrum.mean_metal_energy_eV)<1e-5,"B0 joint Nottingham mean mismatch");}
        long double probability=0,energy=0;for(const auto& node:d.known_component_nodes){probability+=node.probability;energy+=node.probability*(node.normal_eV+node.transverse_eV);}
        need(std::abs(probability-1)<1e-11 && std::abs(energy-d.known_component_mean_metal_energy_eV)<1e-11,"Joint phase-space normalization/heat mismatch");
        double maximum_residual=0;long double drawn_mean=0;int handoffs=0;
        for(int i=0;i<10000;++i){const auto birth=tj::sample_known_component(d,policy,1,(i+.5)/10000,i+1);drawn_mean+=birth.state.particle.normal_metal_energy_eV+birth.state.particle.transverse_energy_eV;maximum_residual=std::max(maximum_residual,std::abs(birth.exchange.residual_J()));}
        need(maximum_residual<1e-30,"Joint sampled birth energy ledger failed");
        need(std::abs(drawn_mean/10000-d.known_component_mean_metal_energy_eV)<.01,"Same-joint-spectrum drawn heat mean mismatch");
        for(double u:{.1,.5,.9}){const auto birth=tj::sample_known_component(d,policy,1,u,1);const auto flight=se::advance_layer(birth.state,s,policy,s.dt);need(flight.state.status==se::LayerStatus::PICHandedOff || flight.state.status==se::LayerStatus::Recollected,"Bounded new-source micro flight did not finish");if(flight.state.status==se::LayerStatus::PICHandedOff)++handoffs;}
        if(B>0){std::ofstream cycles(out/("inverse-cycle-"+std::to_string(int(T))+'-'+std::to_string(int(B))+'-'+std::to_string(refinement)+".csv"));cycles<<std::setprecision(17)<<"joint_uniform,common_dt_s,handoffs,reentries,epochs,return_time_s,net_wall_heat_J,field_work_J,complete_residual_J,transported_energy_scale_J,quantum_reflections,quantum_absorptions,phase_reconstruction_error_eV,completed_material_return\n";
            const tj::Node* mid=nullptr;
            for(const auto& node:d.known_component_nodes)if(node.probability>1e-8 && node.log_transmission>std::log(.3) && node.log_transmission<std::log(.7) && (!mid || node.probability>mid->probability))mid=&node;
            need(mid,"No partial-transmission quantum return test phase");const double partial_uniform=mid->cdf-mid->probability/2;
            for(double u:{.1,.5,.9,partial_uniform})for(double dt:{1e-14,5e-15,2.5e-15}){const auto receipt=joint_cycle_test::cycle(tj::sample_known_component(d,policy,1,u,1),s,policy,dt,config);if(u==partial_uniform)need(receipt.completed && receipt.quantum_reflections>0 && receipt.reentries>0,"Partial-transmission test did not complete its reflected/inward/material cycle");cycles<<u<<','<<dt<<','<<receipt.handoffs<<','<<receipt.reentries<<','<<receipt.epochs<<','<<receipt.return_time_s<<','<<receipt.total.wall_heat_J<<','<<receipt.total.field_work_J<<','<<receipt.total.residual_J()<<','<<receipt.transported_scale_J<<','<<receipt.quantum_reflections<<','<<receipt.quantum_absorptions<<','<<receipt.max_reconstructed_phase_error_eV<<','<<receipt.completed<<'\n';}}
        bool rejected=false;auto oblique=s;oblique.total_B_T[1]=1;try{tj::evaluate(oblique,config,policy);}catch(const std::invalid_argument&){rejected=true;}need(rejected,"Tangential prototype accepted oblique B");
        rejected=false;auto tangentE=s;tangentE.total_E[0]=1e6;try{tj::evaluate(tangentE,config,policy);}catch(const std::invalid_argument&){rejected=true;}need(rejected,"Tangential prototype accepted non-normal E");
        const std::string id=std::to_string(int(T))+'-'+std::to_string(int(B))+'-'+std::to_string(refinement);
        std::ofstream report(out/("joint-"+id+".json"));report<<std::setprecision(17)<<"{\"passed\":true,\"T_K\":"<<T<<",\"B_tangent_T\":"<<B<<",\"refinement\":"<<refinement<<",\"known_component_current_A_m2\":"<<d.known_component_current_A_m2<<",\"nominal_unmagnetized_current_A_m2\":"<<nominal.spectrum.current_A_m2<<",\"known_component_complete_heat_W_m2\":"<<d.known_component_complete_heat_W_m2<<",\"known_component_mean_metal_energy_eV\":"<<d.known_component_mean_metal_energy_eV<<",\"unresolved_current_upper_A_m2\":"<<d.unresolved_current_upper_A_m2<<",\"unresolved_absolute_heat_upper_W_m2\":"<<d.unresolved_absolute_heat_upper_W_m2<<",\"magnetic_reservoir_energy_scale_over_kT\":"<<d.magnetic_reservoir_energy_scale_over_kT<<",\"known_joint_nodes\":"<<d.known_component_nodes.size()<<",\"unresolved_phase_points\":"<<d.unresolved_quadrature_phase_points<<",\"sampled_birth_max_energy_residual_J\":"<<maximum_residual<<",\"three_samples_handed_off\":"<<handoffs<<",\"oblique_B_rejected\":true,\"physical_foil_model_admitted\":false,\"exact_magnetized_solid_supply_solved\":false,\"GPU_executed\":false}\n";
        std::cout<<"Passed joint prototype "<<id<<" Jknown="<<d.known_component_current_A_m2<<" unknown_upper="<<d.unresolved_current_upper_A_m2<<'\n';
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
