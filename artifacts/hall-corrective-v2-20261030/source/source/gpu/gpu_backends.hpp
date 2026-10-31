#pragma once
// GPU backends of Hall corrective v2, sharing ONE device session and ONE ordered
// stream (hall_gpu::Session from hall-gpu-resident-v1-20261030/common):
//  * GpuFieldSolver   deterministic FP64 Jacobi CG on the factored Hall edge
//                     operator (algorithm and kernels of hall-gpu-resident-v1/
//                     field, driver wrapper replaced by the common session);
//  * GpuMotionProvider the push engine of hall-gpu-resident-v1/push, assembled
//                     into the v2 MotionReceipt with the sealed record and seal.
// The header exposes no CUDA types. Device code is compiled at run time by NVRTC
// with --fmad=false --prec-div=true --prec-sqrt=true --ftz=false for the actual
// device architecture.
#include "../owner.hpp"
#include <filesystem>
#include <memory>

namespace hall_gpu {class Session;}
namespace hall_corrective_v2::gpu {
struct DeviceOptions {
    std::string uuid;                       // required device UUID (GPU-xxxxxxxx-...)
    std::filesystem::path nvrtc_library;    // absolute path of libnvrtc.so.12
};
class Device {
    struct Impl;std::unique_ptr<Impl> p_;
public:
    explicit Device(const DeviceOptions&);
    ~Device();
    Device(const Device&)=delete;Device&operator=(const Device&)=delete;
    hall_gpu::Session& session();
    // Device identity, memory and teardown records as JSON text (execution envelope).
    std::string device_json()const;
    std::string memory_json();
    // Drains the stream and releases everything; idempotent.
    std::string close_json()noexcept;
};
struct GpuSolveStats {long solves=0;double total_s=0;std::vector<int> last_iterations;double last_device_residual=0;};
// Field backend "gpu1-factored-cg". Its descriptor binds the algorithm, launch
// geometry, kernel source digest, NVRTC options and the compiled target.
std::shared_ptr<FieldSolveBackend> make_gpu_field_solver(std::shared_ptr<Device>);
std::string gpu_field_solver_compile_json(const FieldSolveBackend&);
GpuSolveStats gpu_field_solver_stats(const FieldSolveBackend&);

struct GpuMotionStats {
    long runs=0;double upload_carriers_s=0,admit_s=0,device_run_s=0,assemble_s=0,total_s=0;
    double last_gpu_total_ms=0,last_host_fold_ms=0;std::uint32_t last_edge_records=0;
};
class GpuMotionProvider final:public MotionProvider {
    struct Impl;std::unique_ptr<Impl> p_;
public:
    GpuMotionProvider(std::shared_ptr<Device>,const MotionContext&,std::size_t capacity);
    ~GpuMotionProvider()override;
    const std::string& name()const noexcept override;
    void begin(const MotionContext&,std::span<const native::ArenaRateEntry> census)override;
    hall_driver::MotionReceipt motion(const hall_driver::MotionBinding&,const Vec&E,const Vec&B,FreeCharge&after)override;
    void end()noexcept override;
    std::string compile_json()const;
    const GpuMotionStats& stats()const noexcept;
};
} // namespace hall_corrective_v2::gpu
