#include "surface_emission.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <sched.h>
#include <sys/resource.h>
namespace se=surface_emission;
namespace {
void check(bool good,const char* text){if(!good)throw std::runtime_error(text);}
se::Stage source(double F,double B) {
    se::Stage s{1,1,1,1,1,1,0,2e-12,0,0,0,{0,-F,0},{0,1,0},1500,1.25e-4,true,true,true};s.total_B_T={B,0,0};return s;
}
void add(se::LayerExchange& a,const se::LayerExchange& b) {
    a.physical_electron_number_change+=b.physical_electron_number_change;a.particle_charge_change_C+=b.particle_charge_change_C;a.electrode_charge_change_C+=b.electrode_charge_change_C;
    a.wall_heat_J+=b.wall_heat_J;a.particle_K_change_J+=b.particle_K_change_J;a.electron_reference_change_J+=b.electron_reference_change_J;a.image_energy_change_J+=b.image_energy_change_J;a.field_work_J+=b.field_work_J;
}
double norm(const se::Vector& r){return std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);}
}
int main(int argc,char** argv) {
    try {
        check(argc==2,"Expected report directory");const std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
        cpu_set_t mask;CPU_ZERO(&mask);check(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(7,&mask) && getpriority(PRIO_PROCESS,0)>=15,"CPU7/nice15 required");
        std::ofstream csv(out/"micro-layer-cases.csv");csv<<std::setprecision(17)<<"F_V_m,B_T,eta_energy,status,substeps,paths,flight_time_s,normal_distance_m,lateral_distance_m,quantum_magnetic_fraction,wall_heat_J,particle_K_change_J,image_change_J,field_work_J,residual_J\n";
        int cases=0,returned=0,handed=0,total_steps=0;double max_residual=0;
        for(const auto& condition:std::vector<std::array<double,4>>{{1e8,0,.03,.1},{1e8,30,.03,.1},{1e7,30,.003,.1},{1e8,30,.01,.5}}) {
            auto s=source(condition[0],condition[1]);const auto flux=se::prepare(s,se::Config{4.3});se::LayerPolicy policy;
            policy.local_linear_vacuum_extent_m=25e-6;policy.image_energy_fraction=condition[2];policy.maximum_quantum_magnetic_energy_fraction=.3;
            const auto birth=se::begin_layer(flux,policy,1,.6,.35,condition[3],cases+1);se::LayerExchange ledger=birth.exchange;
            check(birth.state.status==se::LayerStatus::Active,"Positive-F source did not create an actual active micro-layer packet");
            const auto flight=se::advance_layer(birth.state,s,policy,s.dt);add(ledger,flight.exchange);
            const auto status=flight.state.status;if(status==se::LayerStatus::Recollected)++returned;if(status==se::LayerStatus::PICHandedOff)++handed;
            check(status!=se::LayerStatus::Active,"Bounded manufactured packet did not finish flight/collection");
            // Scale a closed emit/return cycle by transported energy, not its
            // cancellation-to-zero net ledger. Production/work guards unchanged.
            double scale=std::abs(birth.exchange.wall_heat_J)+std::abs(birth.exchange.particle_K_change_J)+
                std::abs(birth.exchange.field_work_J)+std::abs(birth.exchange.electron_reference_change_J)+std::abs(birth.exchange.image_energy_change_J)+
                std::abs(flight.exchange.wall_heat_J)+std::abs(flight.exchange.particle_K_change_J)+std::abs(flight.exchange.field_work_J)+
                std::abs(flight.exchange.electron_reference_change_J)+std::abs(flight.exchange.image_energy_change_J);
            const double residual=std::abs(ledger.residual_J())/std::max(scale,1e-30);max_residual=std::max(max_residual,residual);check(residual<2e-10,"Birth/flight/handoff/recollection energy not closed");
            check(ledger.particle_charge_change_C+ledger.electrode_charge_change_C==0,"Layer particle/electrode charge mismatch");
            for(const auto& path:flight.paths){check(path.end_time_s>=path.start_time_s,"Flight time reversed");check(norm(path.to_m)<=policy.local_linear_vacuum_extent_m,"Layer skipped admitted geometry");}
            if(status==se::LayerStatus::Recollected)check(ledger.physical_electron_number_change==0 && std::abs(ledger.image_energy_change_J)<1e-28,"Returned packet not removed from one inventory/image reference");
            if(condition[1]>0)check(std::abs(flight.state.position_from_surface_m[2])>1e-10,"Strong B flight erased lateral deflection");
            csv<<s.total_E[1]*-1<<','<<condition[1]<<','<<condition[2]<<','<<int(status)<<','<<flight.substeps<<','<<birth.paths.size()+flight.paths.size()<<','<<flight.state.current_time_s<<','<<flight.state.position_from_surface_m[1]<<','<<std::hypot(flight.state.position_from_surface_m[0],flight.state.position_from_surface_m[2])<<','<<birth.state.quantum_magnetic_energy_fraction<<','<<ledger.wall_heat_J<<','<<ledger.particle_K_change_J<<','<<ledger.image_energy_change_J<<','<<ledger.field_work_J<<','<<ledger.residual_J()<<'\n';
            total_steps+=flight.substeps;++cases;
        }
        check(returned>0 && handed>0,"Tests did not exercise both magnetic return and actual PIC handoff");
        auto s=source(1e8,30);const auto f=se::prepare(s,se::Config{4.3});se::ConditionalBoundary straight;
        straight.local_linear_vacuum_extent_m=25e-6;straight.source_log_current_timescale_s=1e-9;
        bool rejected=false;try{se::sample_conditional(f,straight,1,.6,.35,.1,999);}catch(const std::exception&){rejected=true;}check(rejected,"Strong B straight matching was silently admitted");
        se::LayerPolicy policy;policy.local_linear_vacuum_extent_m=25e-6;policy.maximum_quantum_magnetic_energy_fraction=.3;
        const auto birth=se::begin_layer(f,policy,1,.6,.35,.1,1000);const auto copy=birth.state;
        const auto partial=se::advance_layer(birth.state,s,policy,1e-15);check(partial.state.status==se::LayerStatus::Active && partial.state.current_time_s==1e-15,"In-flight particle was lost/fast-forwarded instead of persisted");
        check(birth.state.current_time_s==copy.current_time_s && birth.state.particle.proper_velocity_m_s==copy.particle.proper_velocity_m_s,"Prepare mutated accepted layer state");
        const auto repeat=se::advance_layer(copy,s,policy,1e-15);check(repeat.state.position_from_surface_m==partial.state.position_from_surface_m && repeat.state.particle.proper_velocity_m_s==partial.state.particle.proper_velocity_m_s,"Retry/restart from same immutable layer differs");
        policy.PIC_already_applies_microscopic_single_electron_image_force=true;rejected=false;try{se::begin_layer(f,policy,1,.6,.35,.1,1001);}catch(const std::exception&){rejected=true;}check(rejected,"Microscopic image force double counted");
        policy.PIC_already_applies_microscopic_single_electron_image_force=false;policy.local_linear_vacuum_extent_m=1e-12;rejected=false;try{se::begin_layer(f,policy,1,.6,.35,.1,1002);}catch(const std::exception&){rejected=true;}check(rejected,"Unexplained barrier/geometry jump admitted");
        std::ofstream report(out/"micro-qualification.json");report<<std::setprecision(17)<<"{\"passed\":true,\"CPU\":7,\"GPU_executed\":false,\"physical_foil_model_admitted\":false,\"completed_cases\":"<<cases<<",\"handed_to_PIC\":"<<handed<<",\"recollected\":"<<returned<<",\"microsteps\":"<<total_steps<<",\"max_scaled_complete_energy_residual\":"<<max_residual<<",\"inflight_persisted\":true,\"prepare_pure_repeatable\":true,\"strong_B_straight_projection_rejected\":true}\n";
        std::cout<<"{\"passed\":true,\"micro_layer_cases\":"<<cases<<",\"handed\":"<<handed<<",\"returned\":"<<returned<<"}\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
