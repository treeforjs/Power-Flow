// Prepared benchmark adapter only. Execution is parent-owned and remains on
// HOLD until fresh v2 review. There is no integration into the application.
#include "native_reference.hpp"
#include "cuda_batched_force.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sched.h>
#include <sys/resource.h>

using namespace pic_batch;
namespace io=snapshot_io;
namespace gp=gather_push;
using Clock=std::chrono::steady_clock;
struct Sample {double reference_ms,gpu_api_ms,commit_ms,total_ms;BatchTiming gpu;};
namespace {
void save_cleanup(const std::filesystem::path& output,const CleanupReport& cleanup,const char* original_error=nullptr) {
    std::filesystem::create_directories(output);std::ofstream file(output/"cleanup.json");
    file<<"{\"issue_count\":"<<cleanup.issue_count<<",\"resources_released\":"<<(cleanup.resources_released?"true":"false")
        <<",\"resources_retained\":"<<(cleanup.resources_retained?"true":"false")<<",\"caller_context_restored\":"<<(cleanup.caller_context_restored?"true":"false");
    if(original_error)file<<",\"original_error\":"<<std::quoted(original_error);
    file<<",\"issues\":[";
    for(std::size_t i=0;i<cleanup.issue_count && i<cleanup.issues.size();++i)
        file<<(i?",":"")<<"{\"operation\":"<<std::quoted(cleanup.issues[i].operation)<<",\"code\":"<<cleanup.issues[i].code<<'}';
    file<<"]}\n";io::require(bool(file),"Cannot save cleanup report");
}
}
int main(int argc,char** argv) {
    CleanupReport cleanup;
    try {
        io::require(argc==6 && std::string(argv[4])=="--parent-authorized" && std::string(argv[5])=="--standalone-primary",
                    "Prepared benchmark requires ACTUAL_FIXTURE PTX OUTPUT --parent-authorized --standalone-primary after fresh review");
        cpu_set_t mask;CPU_ZERO(&mask);
        io::require(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(7,&mask) && getpriority(PRIO_PROCESS,0)>=15,"Parent GPU timing requires CPU7/nice15 and external GPU/core7 leases");
        const auto setup_start=Clock::now();const auto snapshot=HostSnapshot::load(argv[1],true);
        native_reference::NativeReference legacy(snapshot,argv[1],"<workspace>/outputs/gap-pulse-midplane110kV.csv");
        std::vector<gp::Outputs> fresh_reference(snapshot.particles.size()),first_result(snapshot.particles.size());
        auto live=snapshot.particles;std::vector<Sample> samples;samples.reserve(20);
        CudaBatchedForce backend(argv[2],snapshot,snapshot.particles.size(),{ContextPolicy::StandalonePrimary,&cleanup});
        const double setup_ms=std::chrono::duration<double,std::milli>(Clock::now()-setup_start).count();
        const auto loop_start=Clock::now();
        for(int round=0;round<23;++round) {
            io::require(std::chrono::duration<double>(Clock::now()-loop_start).count()<15,"Reference-inclusive GPU timing15-second loop guard exceeded");
            std::copy(snapshot.particles.begin(),snapshot.particles.end(),live.begin());
            const auto begin=Clock::now();
            // This cost is mandatory under the CURRENT backend API. Do not
            // replace it with a stored reference or the same shared core.
            legacy.generate(snapshot.particles,fresh_reference);
            const auto generated=Clock::now();
            const auto result=backend.evaluate(snapshot.particles,snapshot.phi,snapshot.parameters,fresh_reference);
            const auto evaluated=Clock::now();backend.commit_velocities(live);const auto committed=Clock::now();
            auto ms=[](auto a,auto b){return std::chrono::duration<double,std::milli>(b-a).count();};
            const Sample sample{ms(begin,generated),ms(generated,evaluated),ms(evaluated,committed),ms(begin,committed),backend.timing()};
            // Golden fixture checks are qualification instrumentation outside
            // the measured API pipeline. API reference gates are INCLUDED above.
            io::require(std::memcmp(fresh_reference.data(),snapshot.reference.data(),fresh_reference.size()*sizeof(gp::Outputs))==0,"Fresh native reference differs from captured E/B/u");
            if(round==0)std::copy(result.begin(),result.end(),first_result.begin());
            else io::require(std::memcmp(first_result.data(),result.data(),first_result.size()*sizeof(gp::Outputs))==0,"Repeated read-only GPU result changed");
            for(std::size_t i=0;i<live.size();++i) {
                auto restored=live[i];restored.proper_velocity_m_s=snapshot.particles[i].proper_velocity_m_s;
                io::require(std::memcmp(&restored,&snapshot.particles[i],sizeof(restored))==0,"Commit changed metadata");
            }
            if(round>=3)samples.push_back(sample);
        }
        const std::filesystem::path output=argv[3];std::filesystem::create_directories(output);
        io::write(output/"device.bin",first_result);
        std::ofstream comparison(output/"comparison.json");batch_compare::print_comparison(backend.comparison(),comparison);comparison.close();
        cleanup=backend.close();
        save_cleanup(output,cleanup);
        io::require(cleanup.issue_count==0 && cleanup.resources_released && !cleanup.resources_retained && cleanup.caller_context_restored,"Cleanup failed; do not qualify this benchmark");
        std::ofstream json(output/"reference-inclusive-timing.json");json<<std::setprecision(17)
            <<"{\"passed\":true,\"particles\":"<<snapshot.particles.size()<<",\"setup_excluded_ms\":"<<setup_ms
            <<",\"fresh_independent_CPU_reference_every_batch\":true,\"reference_generation_charged\":true,\"reference_gates_charged\":true"
            <<",\"commit_charged\":true,\"python_force_loop\":false,\"application_integrated\":false,\"samples_ms\":[";
        for(std::size_t i=0;i<samples.size();++i) {
            const auto& s=samples[i];json<<(i?",":"")<<"{\"native_reference_generation\":"<<s.reference_ms<<",\"GPU_API_including_pack_and_gates\":"<<s.gpu_api_ms
                <<",\"commit\":"<<s.commit_ms<<",\"force_pipeline_total\":"<<s.total_ms<<",\"packing\":"<<s.gpu.pack_ms
                <<",\"H2D\":"<<s.gpu.h2d_ms<<",\"kernel\":"<<s.gpu.kernel_ms<<",\"D2H\":"<<s.gpu.d2h_ms
                <<",\"stage_through_sync\":"<<s.gpu.stage_through_sync_ms<<",\"output_reference_validation\":"<<s.gpu.validation_ms<<'}';
        }
        json<<"],\"scope\":\"Fixed actual midpoint FORCE pipeline, not full timestep. Setup/native refresh and resetting immutable fixture inputs excluded; fresh independent native CPU reference, all backend checks/pack/transfers/sync/gates and transactional commit INCLUDED.\"}\n";
        io::require(bool(json),"Cannot save reference-inclusive timing");
    }catch(const std::exception& error){
        if(argc>=4)try{save_cleanup(argv[3],cleanup,error.what());}catch(const std::exception& report_error){std::cerr<<"Cleanup report write failed: "<<report_error.what()<<'\n';}
        std::cerr<<error.what()<<'\n';return 1;
    }
}
