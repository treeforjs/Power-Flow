// gpu_engine.hpp -- NVRTC backend of the phase owner v1.2 on the sealed hall-gpu-resident-v1 session layer
// (common/hall_gpu_common.hpp: one device, one stream, NVRTC with --fmad=false --prec-div=true --prec-sqrt=true
// --ftz=false). The kernels are the phase_kernels.h / state_digest.h functions; the launch sequence is the
// CpuEngine sequence, so every state word, ledger sum, patch output and Merkle leaf is bit-identical to the CPU
// reference. The engine opens its own session, or runs on a caller-provided session (and therefore on that
// session's stream), e.g. the Hall resident owner's session.
#pragma once
#include "engine.hpp"
#include <filesystem>
#include <string>

namespace hall_gpu { class Session; }

namespace ss304_phase_v1_2 {

struct GpuOptions {
    std::string require_uuid;               // device UUID (GPU-...); empty: ordinal 0 (own session only)
    std::filesystem::path nvrtc_library;    // absolute path of libnvrtc.so.12 (own session only)
    hall_gpu::Session* session = nullptr;   // caller-provided session; the caller keeps it open longer than the engine
    unsigned block = 128;
    bool profile = false;                   // record device time per kernel group of each step (stream markers)
};

class GpuEngine final : public Engine {
public:
    explicit GpuEngine(const GpuOptions& o);
    ~GpuEngine() override;
    std::string name() const override { return "gpu-nvrtc"; }
    void configure(const px_law& law, const px_dims& dims, const HostGeometry& g) override;
    void set_elements(const HostElements& el) override;
    void write_state(int slot, const HostState& s) override;
    HostState read_state(int slot) override;
    void evaluate(int slot, const px_dims& d, const std::vector<double>& patch_pb, EngineResult& r, std::vector<double>* sigma_em) override;
    void step(int from, const px_dims& d, const StepData& in, EngineResult& r) override;
    std::vector<double> temperatures(int slot, const px_dims& d) override;
    std::vector<double> face_displacements(int slot, const px_dims& d) override;
    PxStateRoots state_roots(int slot, const px_dims& d, std::vector<std::array<px_u32, 8>>* cell_leaves) override;
    void accumulate(const std::vector<double>& power, const std::vector<double>& JB, double dt) override;
    void reset_accumulator() override;
    // power [n_em] W and JB [6 n_em] (or null) on the device of the engine's session; read in stream order
    void accumulate_device(const double* power, const double* JB, double dt) override;
    std::string record_json() const;
    // accumulated device milliseconds per kernel group and the number of profiled steps (profile = true)
    std::string profile_json() const;
    double last_step_ms() const { return last_step_ms_; }
    bool borrowed_session() const;
private:
    struct Impl;
    Impl* p_;
    double last_step_ms_ = 0;
};

}  // namespace ss304_phase_v1_2
