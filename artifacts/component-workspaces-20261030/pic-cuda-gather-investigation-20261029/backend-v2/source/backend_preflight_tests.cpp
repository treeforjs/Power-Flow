#include "cuda_batched_force.hpp"
#include <iostream>

int main(int argc,char** argv) {
 try {
    snapshot_io::require(argc==2,"Expected actual midpoint fixture");
    const auto snapshot=pic_batch::HostSnapshot::load(argv[1],true);
    auto rejected=[&](std::size_t capacity,const std::string& expected) {
        bool pass=false;
        try {pic_batch::CudaBatchedForce forbidden("not-loaded.ptx",snapshot,capacity);}
        catch(const std::exception& e){pass=std::string(e.what()).find(expected)!=std::string::npos;}
        snapshot_io::require(pass,"Pre-CUDA resource check did not reject expected invalid configuration");
    };
    // Both exits are in the public constructor before its Impl, where cuInit
    // lives, can be created. No CUDA API is called by these tests.
    rejected(0,"capacity");rejected(2000000,"before CUDA initialization");
    std::cout<<"{\"passed\":true,\"GPU_initialized_or_executed\":false,\"capacity_and_memory_rejected_before_CUDA\":true}\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
