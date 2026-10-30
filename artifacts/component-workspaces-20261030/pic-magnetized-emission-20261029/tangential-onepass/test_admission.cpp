#include "tangential_joint.hpp"
#include "inverse_cycle.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sched.h>
#include <sys/resource.h>

namespace se=surface_emission;namespace tj=tangential_onepass;
int main(int argc,char** argv) {
    try {
        auto need=[](bool b,const char* s){if(!b)throw std::runtime_error(s);};
        cpu_set_t mask;CPU_ZERO(&mask);need(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(7,&mask) && getpriority(PRIO_PROCESS,0)>=15,"CPU7/nice15 required");need(argc==2,"Expected output directory");
        se::Stage s;s.accepted_state=s.field_accepted_state=s.thermal_accepted_state=s.mesh_revision=s.material_revision=s.patch_id=1;s.T_K=1490;s.total_E={0,-97724679.40458213,0};s.total_B_T={30,0,0};s.outward_normal={0,1,0};s.area_m2=1.25e-4;s.dt=3e-12;s.material_temperature_phase_admitted=s.locally_planar_barrier_admitted=s.local_equilibrium_electron_supply_admitted=true;
        tj::Config cfg{4.3};cfg.continuum_free_electron_reservoir_selected=cfg.incoherent_onepass_matching_selected=true;
        se::LayerPolicy p;p.local_linear_vacuum_extent_m=25e-6;
        const auto nominal=tj::evaluate(s,cfg,p);p.image_energy_fraction=.01;const auto tighter=tj::evaluate(s,cfg,p);
        const double current_change=std::abs(tighter.known_component_current_A_m2/nominal.known_component_current_A_m2-1),heat_change=std::abs(tighter.known_component_complete_heat_W_m2/nominal.known_component_complete_heat_W_m2-1);
        need(current_change<1e-6 && heat_change<1e-6,"Matching-distance source rate/heat sensitivity exceeded declared benchmark");
        int cycles=0;for(double u:{.113,.527,.889}){const auto receipt=joint_cycle_test::cycle(tj::sample_known_component(tighter,p,1,u,1),s,p,5e-15,cfg);need(receipt.completed,"Tighter matching representative return remained in flight");++cycles;}
        p.image_energy_fraction=.03;auto poor=s;poor.total_E={0,-5e7,0};poor.total_B_T={40,0,0};const auto unresolved=tj::evaluate(poor,cfg,p);
        need(!unresolved.truncation_bounds_within_declared_budget && unresolved.unresolved_current_upper_A_m2>unresolved.known_component_current_A_m2,"Large unresolved contribution was relabeled as zero or admitted");
        bool rejected=false;try{tj::sample_known_component(unresolved,p,1,.4,1);}catch(const std::invalid_argument&){rejected=true;}need(rejected,"Sampler ignored a large unresolved current bound");
        auto undeclared=cfg;undeclared.continuum_free_electron_reservoir_selected=false;rejected=false;try{tj::evaluate(s,undeclared,p);}catch(const std::invalid_argument&){rejected=true;}need(rejected,"Undeclared electronic reservoir was accepted");
        const std::filesystem::path out=argv[1];std::filesystem::create_directories(out);std::ofstream report(out/"admission-and-matching.json");report<<std::setprecision(17)<<"{\"passed\":true,\"eta_energy_nominal\":0.03,\"eta_energy_tighter\":0.01,\"matching_current_relative_change\":"<<current_change<<",\"matching_heat_relative_change\":"<<heat_change<<",\"tighter_matching_completed_reciprocal_cycles\":"<<cycles<<",\"unresolved_case_F_V_m\":50000000,\"unresolved_case_B_tangent_T\":40,\"unresolved_case_known_current_A_m2\":"<<unresolved.known_component_current_A_m2<<",\"unresolved_case_current_upper_A_m2\":"<<unresolved.unresolved_current_upper_A_m2<<",\"unresolved_case_sampler_refused\":true,\"undeclared_reservoir_refused\":true,\"physical_foil_admitted\":false}\n";
        std::cout<<"Admission and inverse matching checks passed\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
