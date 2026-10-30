#include "cuda_batched_force.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sched.h>
#include <sys/resource.h>

namespace io=snapshot_io;
namespace gp=gather_push;
using namespace pic_batch;
int main(int argc,char** argv) {
 try {
    std::filesystem::path fixture,ptx,out;bool authorized=false,validate_only=false;int warmup=3,repeats=20;
    for(int i=1;i<argc;++i) {
        const std::string arg=argv[i];
        if(arg=="--parent-authorized")authorized=true;
        else if(arg=="--validate-only")validate_only=true;
        else {
            io::require(i+1<argc,"Missing option value");const std::string value=argv[++i];
            if(arg=="--fixture")fixture=value;else if(arg=="--ptx")ptx=value;else if(arg=="--out")out=value;
            else if(arg=="--warmup")warmup=std::stoi(value);else if(arg=="--repeats")repeats=std::stoi(value);
            else throw std::invalid_argument("Unknown benchmark option");
        }
    }
    io::require(!fixture.empty() && !out.empty() && (validate_only || (authorized && !ptx.empty())),"Expected actual --fixture, --out and parent authorization/PTX for GPU work");
    io::require(warmup>=0 && warmup<=5 && repeats>=1 && repeats<=40,"Bounded C++ benchmark iteration guard");
    auto snapshot=HostSnapshot::load(fixture,true);
    std::filesystem::create_directories(out);
    if(validate_only) {
        std::vector<gp::Outputs> result;result.reserve(snapshot.inputs.size());const auto tables=snapshot.tables();
        for(const auto& q:snapshot.inputs)result.push_back(gp::evaluate(tables,q));
        const auto comparison=batch_compare::compare(snapshot.reference,result,snapshot.particles);
        batch_compare::print_comparison(comparison);io::require(comparison.failures==0,"Actual midpoint CPU core differs from original callback reference");
        io::write(out/"shared_cpu.bin",result);return 0;
    }
    cpu_set_t cpus;CPU_ZERO(&cpus);io::require(sched_getaffinity(0,sizeof(cpus),&cpus)==0 && CPU_COUNT(&cpus)==1 && CPU_ISSET(7,&cpus) && getpriority(PRIO_PROCESS,0)>=15,"Parent GPU backend benchmark requires CPU7/nice15");
    const auto setup_start=std::chrono::steady_clock::now();
    CudaBatchedForce backend(ptx,snapshot,snapshot.particles.size());
    const double setup_ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-setup_start).count();
    std::vector<BatchTiming> timings;timings.reserve(repeats);
    const auto loop_start=std::chrono::steady_clock::now();std::span<const gp::Outputs> result;
    for(int i=0;i<warmup+repeats;++i) {
        io::require(std::chrono::duration<double>(std::chrono::steady_clock::now()-loop_start).count()<=5,"C++ steady benchmark five-second guard exceeded");
        result=backend.evaluate(snapshot.particles,snapshot.phi,snapshot.parameters,snapshot.reference);
        if(i>=warmup)timings.push_back(backend.timing());
        if(i+1<warmup+repeats)backend.discard();
    }
    io::write(out/"device.bin",std::vector<gp::Outputs>(result.begin(),result.end()));
    batch_compare::print_comparison(backend.comparison());
    std::ofstream comparison_file(out/"comparison.json");batch_compare::print_comparison(backend.comparison(),comparison_file);
    io::require(bool(comparison_file),"Cannot write independent C++ comparison receipt");
    auto committed=snapshot.particles;
    const auto commit_start=std::chrono::steady_clock::now();backend.commit_velocities(committed);
    const double commit_ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-commit_start).count();
    for(std::size_t i=0;i<committed.size();++i) {
        auto restored=committed[i];restored.proper_velocity_m_s=snapshot.particles[i].proper_velocity_m_s;
        io::require(std::memcmp(&restored,&snapshot.particles[i],sizeof(restored))==0,"Commit changed identity/position/q/m/weight");
        io::require(std::memcmp(committed[i].proper_velocity_m_s.data(),result[i].u.v,24)==0,"Commit changed returned proper velocity");
    }
    std::ofstream json(out/"cpp-timing.json");json<<std::setprecision(17)
        <<"{\"gpu_executed\":true,\"real_midpoint_fixture\":true,\"particles\":"<<snapshot.particles.size()
        <<",\"setup_ms\":"<<setup_ms<<",\"explicit_device_bytes\":"<<backend.explicit_device_bytes()
        <<",\"commit_ms\":"<<commit_ms<<",\"reference_gate_every_batch\":true,\"per_step_allocations\":0"
        <<",\"python_inner_loop\":false,\"full_application_integrated\":false,\"samples\":[";
    for(std::size_t i=0;i<timings.size();++i) {
        const auto& t=timings[i];json<<(i?",":"")<<"{\"pack_ms\":"<<t.pack_ms<<",\"h2d_ms\":"<<t.h2d_ms
            <<",\"kernel_ms\":"<<t.kernel_ms<<",\"d2h_ms\":"<<t.d2h_ms<<",\"stage_through_sync_ms\":"<<t.stage_through_sync_ms
            <<",\"validation_ms\":"<<t.validation_ms<<",\"total_with_validation_ms\":"<<t.total_with_validation_ms<<'}';
    }
    json<<"],\"limitations\":\"Actual midpoint offline reference is already captured; generating a fresh CPU reference and full PIC advance are excluded. Static snapshot refresh costs are setup, not per-step samples.\"}\n";
    io::require(bool(json),"Cannot write C++ timing receipt");
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 2;}
}
