#pragma once
#include "batch_snapshot.hpp"
#include "batch_compare.hpp"
#include <memory>

namespace pic_batch {
struct BatchTiming {
    double pack_ms=0,h2d_ms=0,kernel_ms=0,d2h_ms=0,stage_through_sync_ms=0;
    double validation_ms=0,total_with_validation_ms=0;
};
// Optional qualification backend, not connected to GapPic::advance. The CUDA
// kernel and host core share gather_push_core.hpp; CPU references are captured
// independently from the unchanged application's actual callbacks/pic_boris.
// Reference validation is mandatory for every batch in this first C++ version.
// Construction is the first operation that initializes CUDA. One owning host
// thread/context, explicit static snapshot refresh, no per-step allocations.
class CudaBatchedForce {
public:
    CudaBatchedForce(const std::filesystem::path& ptx,const HostSnapshot&,std::size_t capacity);
    ~CudaBatchedForce();
    CudaBatchedForce(const CudaBatchedForce&)=delete;
    CudaBatchedForce& operator=(const CudaBatchedForce&)=delete;
    // Re-upload on ANY layout/unit/native-frame change. A size change requires
    // a new backend. There is no inference from mutable callback addresses.
    void set_snapshot(const HostSnapshot&);
    std::span<const gather_push::Outputs> evaluate(
        std::span<const raikou::mfem_backend::GapParticle>,std::span<const double> phi,
        const Parameters& midpoint,std::span<const gather_push::Outputs> cpu_reference);
    // Transactional host commit: exact input metadata/state must still match.
    // Validate whole batch before changing only proper_velocity_m_s.
    void commit_velocities(std::span<raikou::mfem_backend::GapParticle>);
    void discard();
    const BatchTiming& timing() const;
    const batch_compare::Comparison& comparison() const;
    std::size_t explicit_device_bytes() const;
private:
    struct Impl;
    std::unique_ptr<Impl> p_;
};
}
