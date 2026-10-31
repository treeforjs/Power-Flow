// ms_gpu.hpp -- device execution of the merge/split transaction (NVRTC kernels
// gpu/ms_kernels.cu over the hall-gpu-resident-v1 common session and deterministic
// primitives). Produces the same Result as merge_split::transact, bit for bit.
#pragma once
#include "hall_gpu_common.hpp"
#include "hall_gpu_primitives.hpp"
#include "merge_split.hpp"

#include <memory>

namespace merge_split {

struct GpuTimings {
    double upload_ms = 0;       // host -> device (pinned staging, one copy)
    double sort_ms = 0;         // classify, two stable radix sorts, gather, segment bounds, group counts, scan
    double kinematics_ms = 0;   // group descriptors, kinematics, output/relocation scans
    double emit_ledger_ms = 0;  // emission, removed-ID compaction, exact ledger reduction
    double download_ms = 0;     // device -> host of every result array
    double device_compute_ms() const { return sort_ms + kinematics_ms + emit_ledger_ms; }
    double host_wall_ms = 0;    // whole run() call (allocation, synchronizations, unpacking)
    std::size_t device_bytes = 0, pinned_bytes = 0;
};

class GpuTransaction {
public:
    explicit GpuTransaction(hall_gpu::Session&);
    ~GpuTransaction();
    GpuTransaction(const GpuTransaction&) = delete;
    GpuTransaction& operator=(const GpuTransaction&) = delete;
    // Carriers must be in strictly ascending native ID order (the native census order).
    Result run(const MsParams&, const Carriers&, GpuTimings* timings = nullptr);
    const hall_gpu::Module& module() const;

private:
    struct Impl;
    std::unique_ptr<Impl> p_;
};

}  // namespace merge_split
