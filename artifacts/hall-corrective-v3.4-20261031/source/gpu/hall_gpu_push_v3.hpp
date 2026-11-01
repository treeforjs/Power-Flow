// hall_gpu_push_v3.hpp -- GPU backend of the Hall v3 motion path (first-hit kinetic
// sinks). Same engine as hall-gpu-resident-v1-20261030/push MotionEngine (residency,
// one ordered stream, sorted segmented reductions, no float atomics, host folds in
// the reference refusal order) with the v3 carrier body and endpoint reduction, and
// the per-carrier sink outputs. Types of the v1 header are reused unchanged.
#pragma once
#include "push/hall_gpu_push.hpp"
#include "work_scales_v32.hpp"

namespace hall_gpu::push_v3 {
using namespace hall_gpu::push;

struct MotionViewV3 : MotionView {
    WorkScalesV32 work;   // v3.2 work scales and tolerance of the fold
    std::span<const int> sink_kind, sink_axis;   // -1 for surviving carriers
    std::span<const double> sink_flight, sink_remaining, sink_fraction, sink_plane, sink_segment_x, sink_segment_y, sink_segment_z, gyro, dcells;
};
// v3 sink geometry (filled by the provider from the Hall motion context).
struct SinkSupport { int lo[5][3]{}, hi[5][3]{}, supported[5]{}; bool aligned = false; };

// Reference text of a status code: the v1 texts, and the v3 sink texts.
const char* code_message_v3(int code);

class MotionEngineV3 {
public:
    // v3.4: module = a compiled push program of this session (shared_module), or nullptr to compile.
    MotionEngineV3(Session&, const MeshTables&, const Geometry&, const SinkSupport&, const EngineOptions&, const Module* module = nullptr);
    // The push program compiled once per device (hall_corrective_v3::gpu::Device::cached).
    template <class D> static std::shared_ptr<const Module> shared_module(D& device) {
        return std::static_pointer_cast<const Module>(device.cached("hall-v3-push-module", [&] {
            return std::static_pointer_cast<void>(std::make_shared<Module>(compile(device.session()))); }));
    }
    static Module compile(Session&);
    ~MotionEngineV3();
    MotionEngineV3(const MotionEngineV3&) = delete;
    MotionEngineV3& operator=(const MotionEngineV3&) = delete;
    const CompileRecord& compile_record() const noexcept;
    void set_carriers(const CarrierView&);
    void set_field(std::span<const double> E, std::span<const double> B);
    // Per-carrier flight (dt; newborn end - tau), newborn face start (-1) and newborn flag,
    // in carrier (ID) order; empty spans = residents only.
    void set_flights(std::span<const double> flight, std::span<const int> face, std::span<const char> newborn);
    // v3.2: native v7 motion forms (exact kinetic change and field work in the host fold); default false.
    void set_motion_forms(bool native_v7);
    MotionViewV3 run(double dt, bool diagnostics);
    std::size_t device_bytes() const noexcept;
    // ---- v3.4 device path (device_step_v34.cpp): carriers and flights written on the device by the caller,
    // the trial field read from device buffers, every output left on the device.
    struct DeviceRun {
        std::size_t n = 0;
        std::uint64_t outs = 0;          // [15][n]: endpoint x y z ux uy uz, start x y z, disp x y z, dke, gyro, dcells
        std::uint64_t nodes = 0;         // [8][nn]: omega, initial charge, gross, endpoint charge, net hi, net lo, free gross, plasma
        std::uint64_t current = 0;       // [ne]
        std::uint64_t sinks_i = 0;       // [2][n]: sink kind, axis
        std::uint64_t sinks_r = 0;       // [7][n]: flight, discarded, fraction, plane, segment x y z
        std::uint64_t status = 0;        // 8 u32: census rank/code, motion rank/code, charge rank/code, edge records, 0
        std::uint64_t maxima = 0;        // 3 f64: plasma dt, gyro, displacement cells
        std::uint32_t edge_records = 0;
        DeviceBuffer outs_buf, nodes_buf, current_buf, sinks_i_buf, sinks_r_buf, status_buf, maxima_buf;   // handles (status at byte 64, maxima at 0)
    };
    std::uint64_t carriers_address() const noexcept;   // [9][capacity]: SoA with stride n
    std::uint64_t flights_address() const noexcept;    // flight [n] f64, face [n] i32, newborn [n] i32
    std::size_t capacity() const noexcept;
    void use_device_carriers(std::size_t n, bool flights);
    void set_field_device(std::uint64_t E, std::uint64_t B);
    DeviceRun run_device(double dt);
    // Stream timings of the last run_device (call after the caller's synchronize).
    std::vector<Phase> device_phases();
private:
    struct Impl;
    Impl* p_;
};
}  // namespace hall_gpu::push_v3
