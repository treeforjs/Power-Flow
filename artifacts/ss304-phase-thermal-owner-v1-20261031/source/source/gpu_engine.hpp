// gpu_engine.hpp -- NVRTC backend of the phase owner on the sealed hall-gpu-resident-v1 session layer
// (common/hall_gpu_common.hpp: one device, one stream, NVRTC with --fmad=false --prec-div=true
// --prec-sqrt=true --ftz=false). The kernels are the phase_kernels.h functions; the launch sequence is
// CpuEngine::step/evaluate, so every state word, ledger sum and patch output is bit-identical to the CPU
// reference. State lives on the device; per step the host uploads the EM Joule energies, the sparse surface
// inputs and the patch release totals, and downloads the per-slice sums/statuses and the per-patch outputs.
#pragma once
#include "engine.hpp"
#include <filesystem>
#include <string>

namespace hall_gpu { class Session; }

namespace ss304_phase_v1 {

struct GpuOptions {
    std::string require_uuid;               // device UUID (GPU-...); empty: ordinal 0
    std::filesystem::path nvrtc_library;    // absolute path of libnvrtc.so.12
    unsigned block = 128;
};

class GpuEngine final : public Engine {
public:
    explicit GpuEngine(const GpuOptions& o);
    ~GpuEngine() override;
    std::string name() const override { return "gpu-nvrtc"; }
    void configure(const pt_law& law, const pt_dims& dims, const HostGeometry& g) override;
    void set_elements(const HostElements& el) override;
    void write_state(int slot, const HostState& s) override;
    HostState read_state(int slot) override;
    void evaluate(int slot, const pt_dims& d, EngineResult& r, std::vector<double>* sigma_em) override;
    void step(int from, const pt_dims& d, const StepData& in, EngineResult& r) override;
    std::vector<double> temperatures(int slot, const pt_dims& d) override;
    // Device, NVRTC and compile records (JSON text) for receipts.
    std::string record_json() const;
    double last_step_ms() const { return last_step_ms_; }
private:
    struct Impl;
    Impl* p_;
    double last_step_ms_ = 0;
};

}  // namespace ss304_phase_v1
