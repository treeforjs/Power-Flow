#include "native_reference.hpp"
#include "batch_compare.hpp"
#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sched.h>
#include <sys/resource.h>

using namespace raikou::mfem_backend;
namespace io=snapshot_io;
namespace gp=gather_push;
using Clock=std::chrono::steady_clock;
namespace {
void commit(std::span<GapParticle> live,std::span<const GapParticle> before,std::span<const gp::Outputs> result) {
    io::require(live.size()==before.size() && result.size()==before.size(),"Commit probe shape mismatch");
    for(std::size_t i=0;i<live.size();++i)
        io::require(std::memcmp(&live[i],&before[i],sizeof(GapParticle))==0,"Commit probe metadata changed");
    for(std::size_t i=0;i<live.size();++i)live[i].proper_velocity_m_s=io::vec(result[i].u);
}
void verify_velocities(std::span<const GapParticle> after,const pic_batch::HostSnapshot& s) {
    for(std::size_t i=0;i<after.size();++i) {
        auto restored=after[i];restored.proper_velocity_m_s=s.particles[i].proper_velocity_m_s;
        io::require(std::memcmp(&restored,&s.particles[i],sizeof(GapParticle))==0,"Force probe changed identity/position/q/m/weight");
        io::require(std::memcmp(after[i].proper_velocity_m_s.data(),s.reference[i].u.v,24)==0,"Native force probe differs from captured velocity bits");
    }
}
}
int main(int argc,char** argv) {
 try {
    io::require(argc==3,"Usage: host_force_timing ACTUAL_FIXTURE OUT_JSON");
    cpu_set_t mask;CPU_ZERO(&mask);
    io::require(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(6,&mask) && getpriority(PRIO_PROCESS,0)>=15,"CPU timing requires CPU6/nice15 and caller-held testing-core.lock");
    const auto setup_begin=Clock::now();const auto s=pic_batch::HostSnapshot::load(argv[1],true);
    native_reference::NativeReference legacy(s,argv[1],"<workspace>/outputs/gap-pulse-midplane110kV.csv");
    const auto table=s.tables();const auto n=s.particles.size();
    std::vector<gp::Outputs> native(n),shared(n),fresh(n);
    std::vector<gp::Inputs> packed(n);std::vector<GapParticle> work=s.particles;
    const double setup_ms=std::chrono::duration<double,std::milli>(Clock::now()-setup_begin).count();
    legacy.generate(s.particles,native);
    io::require(std::memcmp(native.data(),s.reference.data(),n*sizeof(gp::Outputs))==0,"Independent native reference generation changed captured E/B/u bits");
    for(std::size_t i=0;i<n;++i)shared[i]=gp::evaluate(table,s.inputs[i]);
    io::require(std::memcmp(shared.data(),s.reference.data(),n*sizeof(gp::Outputs))==0,"Shared core differs from actual captured E/B/u bits");
    constexpr int warmups=3,repetitions=21,paths=6;
    const std::array<const char*,paths> names={"native_force_update","native_reference_generate_EBu","shared_core_prepacked_EBu","pack_plus_shared_core_EBu","independent_reference_gate","transactional_velocity_commit"};
    std::array<std::vector<double>,paths> times;for(auto& item:times)item.reserve(repetitions);
    const auto loop_begin=Clock::now();
    for(int round=0;round<warmups+repetitions;++round) {
        for(int order=0;order<paths;++order) {
            io::require(std::chrono::duration<double>(Clock::now()-loop_begin).count()<15,"CPU force timing15-second guard exceeded");
            const int offset=round%paths;const int path=(offset+(round%2?paths-1-order:order))%paths;
            if(path==0 || path==5)std::copy(s.particles.begin(),s.particles.end(),work.begin()); // input reset is excluded uniformly
            batch_compare::Comparison gate;
            const auto begin=Clock::now();
            switch(path) {
                case 0:legacy.force_update(work);break;
                case 1:legacy.generate(s.particles,native);break;
                case 2:for(std::size_t i=0;i<n;++i)shared[i]=gp::evaluate(table,s.inputs[i]);break;
                case 3:
                    for(std::size_t i=0;i<n;++i){const auto& p=s.particles[i];packed[i]={p.position_m[0],p.position_m[1],p.position_m[2],io::vec(p.proper_velocity_m_s),p.charge_C/p.mass_kg};}
                    for(std::size_t i=0;i<n;++i)shared[i]=gp::evaluate(table,packed[i]);
                    break;
                case 4:gate=batch_compare::compare(native,shared,s.particles);break;
                case 5:commit(work,s.particles,shared);break;
            }
            const double ms=std::chrono::duration<double,std::milli>(Clock::now()-begin).count();
            if(round>=warmups)times[path].push_back(ms);
            // Preserve the complete qualification check after every measured
            // variant. Its separately measured cost is charged in pipeline totals.
            if(path==0 || path==5)verify_velocities(work,s);
            else if(path==4)io::require(gate.failures==0,"Shared independent reference gate failed");
            else {
                const auto& candidate=path==1?native:shared;
                io::require(std::memcmp(candidate.data(),s.reference.data(),n*sizeof(gp::Outputs))==0,"Timed force variant differs from captured E/B/u bits");
                io::require(batch_compare::compare(s.reference,candidate,s.particles).failures==0,"Timed componentwise/finite/weighted qualification failed");
            }
        }
    }
    // Separately time a complete, still-reference-gated CPU force pipeline:
    // fresh independent native reference + fresh AoS packing/shared arithmetic
    // + unchanged finite/componentwise/weighted gates + velocity-only commit.
    std::vector<double> complete;complete.reserve(repetitions);
    for(int round=0;round<warmups+repetitions;++round) {
        std::copy(s.particles.begin(),s.particles.end(),work.begin());
        const auto begin=Clock::now();legacy.generate(s.particles,fresh);
        for(std::size_t i=0;i<n;++i){const auto& p=s.particles[i];packed[i]={p.position_m[0],p.position_m[1],p.position_m[2],io::vec(p.proper_velocity_m_s),p.charge_C/p.mass_kg};}
        for(std::size_t i=0;i<n;++i)shared[i]=gp::evaluate(table,packed[i]);
        io::require(batch_compare::compare(fresh,shared,s.particles).failures==0,"Full CPU reference-gated pipeline failed");
        commit(work,s.particles,shared);
        const double ms=std::chrono::duration<double,std::milli>(Clock::now()-begin).count();
        if(round>=warmups)complete.push_back(ms);
        verify_velocities(work,s);io::require(std::memcmp(fresh.data(),s.reference.data(),n*sizeof(gp::Outputs))==0,"Fresh independent reference changed");
    }
    std::ofstream json(argv[2]);json<<std::setprecision(17)<<"{\"passed\":true,\"real_CUDA_initialized_or_executed\":false,\"particles\":"<<n
        <<",\"CPU\":6,\"nice\":"<<getpriority(PRIO_PROCESS,0)<<",\"warmups\":"<<warmups<<",\"repetitions\":"<<repetitions
        <<",\"setup_excluded_ms\":"<<setup_ms<<",\"all_EBu_buffers_bitwise_equal\":true,\"original_gates_preserved\":true,\"samples_ms\":{";
    auto array=[&](const std::vector<double>& data){json<<'[';for(std::size_t i=0;i<data.size();++i)json<<(i?",":"")<<data[i];json<<']';};
    for(int path=0;path<paths;++path){json<<(path?",":"")<<std::quoted(names[path])<<':';array(times[path]);}
    json<<",\"complete_reference_gated_CPU_force_pipeline\":";array(complete);
    json<<"},\"method\":\"Warm fixed actual midpoint; alternating/rotated path order; separate self/unit MFEM meshes and original Fields/Voltage/std::function/pic_boris reference. All storage preallocated; input reset and post-variant golden checks excluded, with gate separately timed and included in complete pipeline. No solver/drift/source step or application integration.\"}\n";
    io::require(bool(json),"Cannot write CPU timing receipt");std::cout<<"CPU-only force timing completed; every arithmetic and qualification comparison passed\n";
 }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
