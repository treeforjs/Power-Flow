#include "evolving_force.hpp"
#include "runtime_guards.hpp"
#include <iostream>
#include <limits>
using namespace pic_batch;
namespace io=snapshot_io;namespace gp=gather_push;
template<class F>void rejects(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}io::require(rejected,"Expected runtime guard rejection");}
int main(int argc,char** argv) {
    try {
        io::require(argc==2,"Expected actual midpoint fixture");const auto s=HostSnapshot::load(argv[1],true);
        const auto state=RuntimeSnapshot::from_host(s);state.validate();
        std::vector<gp::Outputs> result;result.reserve(s.particles.size());
        for(const auto& p:s.inputs)result.push_back(gp::evaluate(s.tables(),p));
        validate_runtime_outputs(result,s.particles);
        const auto c=batch_compare::compare(s.reference,result,s.particles);io::require(c.failures==0,"Original actual-midpoint reference gates failed");
        for(int which=0;which<4;++which){auto bad=result;
            if(which==0)bad.back().E.v[0]=std::numeric_limits<double>::quiet_NaN();
            if(which==1)bad.back().B.v[0]=std::numeric_limits<double>::infinity();
            if(which==2)bad.back().u.v[0]=std::numeric_limits<double>::max();
            if(which==3)bad.back().status=1;
            rejects([&]{validate_runtime_outputs(bad,s.particles);});
        }
        rejects([&]{pic_evolving::candidate_environment();}); // CPU6 ownership gate, before CUDA
        rejects([&]{pic_evolving::validate_output("/home/treefor/sim-workspace/powerflow/outputs/pic-profile/full60");});
        rejects([&]{pic_evolving::bound_restart(0,2.6478441131125727e-08,2.6480441131125727e-08,437411);});
        rejects([&]{pic_evolving::bound_restart(0,2.6478441131125727e-08,2.6479441131125727e-08,437412);});
        pic_evolving::validate_output("/home/treefor/projects/runtimes/pic-cuda-force-evolving-v3-20261029/qualification-runs/private-preflight");
        std::cout<<"{\"passed\":true,\"real_CUDA_initialized_or_executed\":false,\"particles\":"<<s.particles.size()
            <<",\"actual_CPU_EBu_reference_gates\":true,\"reference_free_runtime_finite_status_energy_guards\":true,\"CPU6_deployment_gate\":true}\n";
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
