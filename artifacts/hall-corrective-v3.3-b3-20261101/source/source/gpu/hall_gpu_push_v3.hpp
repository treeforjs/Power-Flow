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
    MotionEngineV3(Session&, const MeshTables&, const Geometry&, const SinkSupport&, const EngineOptions&);
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
private:
    struct Impl;
    Impl* p_;
};
}  // namespace hall_gpu::push_v3
