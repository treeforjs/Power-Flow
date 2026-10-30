#pragma once
#include "cuda_batched_force.hpp" // unchanged v2 context/cleanup/timing types
#include "runtime_snapshot.hpp"
namespace pic_batch {
// NEW, unqualified runtime candidate. The accepted v2 mandatory-reference API
// is untouched. Independent shadow comparisons live outside this executor.
class CudaForceRuntime {
public:
    CudaForceRuntime(const std::filesystem::path&,const RuntimeSnapshot&,std::size_t capacity,BackendOptions={});
    ~CudaForceRuntime();
    CudaForceRuntime(const CudaForceRuntime&)=delete;
    CudaForceRuntime& operator=(const CudaForceRuntime&)=delete;
    void set_snapshot(const RuntimeSnapshot&);
    std::span<const gather_push::Outputs> evaluate(
        std::span<const raikou::mfem_backend::GapParticle>,std::span<const double> phi,const Parameters&);
    void commit_velocities(std::span<raikou::mfem_backend::GapParticle>);
    void discard();
    CleanupReport close() noexcept;
    const BatchTiming& timing() const;
    std::size_t explicit_device_bytes() const;
private:
    struct Impl;
    std::unique_ptr<Impl> p_;
};
}
