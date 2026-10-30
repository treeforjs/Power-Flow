// Minimal CUDA driver API + NVRTC wrapper for the Hall field backend.
// Candidate for common/ (device context, NVRTC helper) once that interface
// exists; kept local here so the field backend does not depend on it yet.
//
// Lifetime and ordering rules (closing the two findings on the earlier
// gather backend):
//  * Ordering: every copy, memset and launch is issued on ONE created
//    non-blocking stream; the legacy default stream is never used. Host
//    sources of uploads are pinned buffers owned by the session; static
//    uploads (operator arrays) are staged through a pinned buffer and the
//    stream is synchronized before the staging buffer is reused and before
//    the operator is reported ready. No pageable-memory copy depends on
//    implicit legacy-stream synchronization.
//  * Teardown: one Session owns every device allocation, pinned buffer,
//    module, graph, event and the stream. close() first makes the session's
//    primary context current, drains the stream (falls back to a context
//    synchronize), and only then frees memory, unloads the module, destroys
//    the stream, pops the context and releases its own primary reference. If
//    completion cannot be established (sticky device error), nothing is freed,
//    the session is marked quarantined and the caller must exit the process.
//    Cleanup never throws and never replaces the original exception; issues
//    are recorded in a caller-owned CleanupReport.
//  * Context: the process must not have a current context; the session
//    retains device 0's primary context (optionally with blocking sync) and
//    never resets or destroys it.
#pragma once
#include <cuda.h>
#include <nvrtc.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <thread>
#include <vector>
#include "csr_operator.hpp"

namespace hgf::cu {

inline std::string describe(CUresult r) {
    const char *name = nullptr, *text = nullptr;
    cuGetErrorName(r, &name);
    cuGetErrorString(r, &text);
    return std::string(name ? name : "CUDA_ERROR") + " (" + std::to_string(int(r)) + "): " + (text ? text : "");
}
inline void ok(CUresult r, const char* what) {
    if (r != CUDA_SUCCESS) throw std::runtime_error(std::string(what) + ": " + describe(r));
}
inline void nvrtc_ok(nvrtcResult r, const char* what) {
    if (r != NVRTC_SUCCESS) throw std::runtime_error(std::string(what) + ": " + nvrtcGetErrorString(r));
}

struct CleanupIssue { std::string operation; int code; };
struct CleanupReport {
    std::vector<CleanupIssue> issues;
    bool closed = false, drained = false, quarantined = false, primary_released = false, context_restored = false;
    std::size_t freed_device = 0, freed_pinned = 0;
};

struct DeviceInfo {
    int ordinal = 0, major = 0, minor = 0, sm_count = 0, driver_api_version = 0;
    std::string name, uuid, pci_bus_id;
    std::size_t total_memory = 0;
};

struct Compiled {
    std::vector<char> image;
    bool ptx = false;
    int nvrtc_major = 0, nvrtc_minor = 0;
    std::vector<int> supported_archs;
    std::vector<std::string> options;
    std::string log, image_sha256, source_sha256, arch_decision;
};

// NVRTC compile for the device's own architecture when NVRTC supports it
// (SASS, no driver JIT), otherwise the newest supported virtual architecture
// below it (PTX, JIT by the driver). FP options are fixed.
inline Compiled compile(const std::string& source, const char* name, int major, int minor) {
    Compiled c;
    nvrtc_ok(nvrtcVersion(&c.nvrtc_major, &c.nvrtc_minor), "nvrtcVersion");
    int count = 0;
    nvrtc_ok(nvrtcGetNumSupportedArchs(&count), "nvrtcGetNumSupportedArchs");
    c.supported_archs.resize(count);
    nvrtc_ok(nvrtcGetSupportedArchs(c.supported_archs.data()), "nvrtcGetSupportedArchs");
    const int cc = major * 10 + minor;
    int chosen = -1;
    for (int a : c.supported_archs) if (a == cc) chosen = a;
    std::string arch;
    if (chosen == cc) {
        arch = "--gpu-architecture=sm_" + std::to_string(cc);
        c.arch_decision = "device architecture supported by NVRTC: SASS for sm_" + std::to_string(cc);
    } else {
        for (int a : c.supported_archs) if (a < cc) chosen = std::max(chosen, a);
        need(chosen > 0, "NVRTC supports no architecture at or below the device");
        arch = "--gpu-architecture=compute_" + std::to_string(chosen);
        c.ptx = true;
        c.arch_decision = "device sm_" + std::to_string(cc) + " not in NVRTC list: PTX for compute_" + std::to_string(chosen) + " (driver JIT)";
    }
    c.options = {arch, "--fmad=false", "--prec-div=true", "--prec-sqrt=true", "--ftz=false", "--std=c++17"};
    c.source_sha256 = sha256_bytes(source.data(), source.size());
    nvrtcProgram prog = nullptr;
    nvrtc_ok(nvrtcCreateProgram(&prog, source.c_str(), name, 0, nullptr, nullptr), "nvrtcCreateProgram");
    std::vector<const char*> opts;
    for (auto& o : c.options) opts.push_back(o.c_str());
    const nvrtcResult res = nvrtcCompileProgram(prog, int(opts.size()), opts.data());
    std::size_t size = 0;
    nvrtcGetProgramLogSize(prog, &size);
    c.log.assign(size, '\0');
    if (size) nvrtcGetProgramLog(prog, c.log.data());
    while (!c.log.empty() && c.log.back() == '\0') c.log.pop_back();
    if (res != NVRTC_SUCCESS) {
        nvrtcDestroyProgram(&prog);
        throw std::runtime_error("NVRTC compilation failed: " + c.log);
    }
    if (c.ptx) {
        nvrtcGetPTXSize(prog, &size); c.image.resize(size); nvrtcGetPTX(prog, c.image.data());
    } else {
        nvrtcGetCUBINSize(prog, &size); c.image.resize(size); nvrtcGetCUBIN(prog, c.image.data());
    }
    nvrtcDestroyProgram(&prog);
    c.image_sha256 = sha256_bytes(c.image.data(), c.image.size());
    return c;
}

class Session {
public:
    Session(int ordinal, bool blocking_sync, CleanupReport& report) : report_(report) {
        try {
            ok(cuInit(0), "cuInit");
            CUcontext current = nullptr;
            ok(cuCtxGetCurrent(&current), "cuCtxGetCurrent");
            need(current == nullptr, "a CUDA context is already current; this standalone session does not borrow contexts");
            ok(cuDeviceGet(&device_, ordinal), "cuDeviceGet");
            info_.ordinal = ordinal;
            char name[256] = {};
            ok(cuDeviceGetName(name, sizeof(name), device_), "cuDeviceGetName");
            info_.name = name;
            CUuuid id;
            ok(cuDeviceGetUuid_v2(&id, device_), "cuDeviceGetUuid_v2");
            char u[64];
            const auto* b = reinterpret_cast<const unsigned char*>(id.bytes);
            std::snprintf(u, sizeof(u), "GPU-%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
                          b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);
            info_.uuid = u;
            char pci[64] = {};
            ok(cuDeviceGetPCIBusId(pci, sizeof(pci), device_), "cuDeviceGetPCIBusId");
            info_.pci_bus_id = pci;
            ok(cuDeviceGetAttribute(&info_.major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, device_), "cc major");
            ok(cuDeviceGetAttribute(&info_.minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, device_), "cc minor");
            ok(cuDeviceGetAttribute(&info_.sm_count, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, device_), "sm count");
            ok(cuDeviceTotalMem(&info_.total_memory, device_), "cuDeviceTotalMem");
            ok(cuDriverGetVersion(&info_.driver_api_version), "cuDriverGetVersion");
            if (blocking_sync) {
                unsigned flags = 0; int active = 0;
                ok(cuDevicePrimaryCtxGetState(device_, &flags, &active), "cuDevicePrimaryCtxGetState");
                need(!active, "primary context already active; scheduling flags cannot be set");
                ok(cuDevicePrimaryCtxSetFlags(device_, CU_CTX_SCHED_BLOCKING_SYNC), "cuDevicePrimaryCtxSetFlags");
            }
            ok(cuDevicePrimaryCtxRetain(&ctx_, device_), "cuDevicePrimaryCtxRetain");
            retained_ = true;
            ok(cuCtxPushCurrent(ctx_), "cuCtxPushCurrent");
            pushed_ = true;
            ok(cuStreamCreate(&stream_, CU_STREAM_NON_BLOCKING), "cuStreamCreate");
            staging_bytes_ = std::size_t(32) << 20;
            staging_ = static_cast<char*>(pinned(staging_bytes_));
            thread_ = std::this_thread::get_id();
        } catch (...) {
            close();
            throw;
        }
    }
    Session(const Session&) = delete;
    Session& operator=(const Session&) = delete;
    ~Session() { close(); }

    const DeviceInfo& info() const { return info_; }
    CUstream stream() const { return stream_; }

    CUdeviceptr device(std::size_t bytes) {
        CUdeviceptr p = 0;
        ok(cuMemAlloc(&p, std::max<std::size_t>(bytes, 8)), "cuMemAlloc");
        device_allocs_.push_back(p);
        return p;
    }
    void* pinned(std::size_t bytes) {
        void* p = nullptr;
        ok(cuMemHostAlloc(&p, std::max<std::size_t>(bytes, 8), 0), "cuMemHostAlloc");
        pinned_allocs_.push_back(p);
        return p;
    }
    // Static upload from any host memory, staged through the pinned buffer on
    // the session stream. Returns only after the copy has completed.
    void upload_static(CUdeviceptr dst, const void* src, std::size_t bytes) {
        const char* s = static_cast<const char*>(src);
        for (std::size_t off = 0; off < bytes; off += staging_bytes_) {
            const std::size_t chunk = std::min(staging_bytes_, bytes - off);
            std::memcpy(staging_, s + off, chunk);
            ok(cuMemcpyHtoDAsync(dst + off, staging_, chunk, stream_), "cuMemcpyHtoDAsync (static)");
            sync();
        }
    }
    // Asynchronous copies between session-owned pinned buffers and device memory.
    void h2d(CUdeviceptr dst, const void* pinned_src, std::size_t bytes) { ok(cuMemcpyHtoDAsync(dst, pinned_src, bytes, stream_), "cuMemcpyHtoDAsync"); }
    void d2h(void* pinned_dst, CUdeviceptr src, std::size_t bytes) { ok(cuMemcpyDtoHAsync(pinned_dst, src, bytes, stream_), "cuMemcpyDtoHAsync"); }
    void d2d(CUdeviceptr dst, CUdeviceptr src, std::size_t bytes) { ok(cuMemcpyDtoDAsync(dst, src, bytes, stream_), "cuMemcpyDtoDAsync"); }
    void zero(CUdeviceptr dst, std::size_t bytes) { ok(cuMemsetD8Async(dst, 0, bytes, stream_), "cuMemsetD8Async"); }
    void sync() { ok(cuStreamSynchronize(stream_), "cuStreamSynchronize"); }

    void load(const Compiled& c) {
        need(module_ == nullptr, "module already loaded");
        std::vector<char> log(16384, '\0');
        CUjit_option keys[2] = {CU_JIT_ERROR_LOG_BUFFER, CU_JIT_ERROR_LOG_BUFFER_SIZE_BYTES};
        void* values[2] = {log.data(), reinterpret_cast<void*>(std::uintptr_t(log.size()))};
        const CUresult r = cuModuleLoadDataEx(&module_, c.image.data(), 2, keys, values);
        if (r != CUDA_SUCCESS) throw std::runtime_error("cuModuleLoadDataEx: " + describe(r) + " " + std::string(log.data()));
    }
    CUfunction function(const char* name) {
        CUfunction f = nullptr;
        ok(cuModuleGetFunction(&f, module_, name), name);
        return f;
    }
    void launch(CUfunction f, unsigned grid, unsigned block, void** args) {
        ok(cuLaunchKernel(f, grid, 1, 1, block, 1, 1, 0, stream_, args, nullptr), "cuLaunchKernel");
    }
    CUevent event() {
        CUevent e = nullptr;
        ok(cuEventCreate(&e, CU_EVENT_DEFAULT), "cuEventCreate");
        events_.push_back(e);
        return e;
    }
    void record(CUevent e) { ok(cuEventRecord(e, stream_), "cuEventRecord"); }
    float elapsed_ms(CUevent a, CUevent b) { float ms = 0; ok(cuEventElapsedTime(&ms, a, b), "cuEventElapsedTime"); return ms; }
    // Graph capture of work issued on the session stream.
    void begin_capture() { ok(cuStreamBeginCapture(stream_, CU_STREAM_CAPTURE_MODE_THREAD_LOCAL), "cuStreamBeginCapture"); }
    CUgraphExec end_capture() {
        CUgraph g = nullptr;
        ok(cuStreamEndCapture(stream_, &g), "cuStreamEndCapture");
        graphs_.push_back(g);
        CUgraphExec e = nullptr;
        ok(cuGraphInstantiate(&e, g, 0), "cuGraphInstantiate");
        execs_.push_back(e);
        return e;
    }
    void launch_graph(CUgraphExec e) { ok(cuGraphLaunch(e, stream_), "cuGraphLaunch"); }

    void close() noexcept {
        if (report_.closed) return;
        report_.closed = true;
        if (thread_ != std::thread::id() && std::this_thread::get_id() != thread_) {
            issue("close from a different thread", CUDA_ERROR_INVALID_CONTEXT);
            report_.quarantined = true;
            print();
            return;
        }
        if (retained_) {
            CUcontext now = nullptr;
            const CUresult g = cuCtxGetCurrent(&now);
            issue("cuCtxGetCurrent", g);
            if (g != CUDA_SUCCESS || now != ctx_) {
                const CUresult p = cuCtxPushCurrent(ctx_);
                issue("cuCtxPushCurrent (cleanup)", p);
                if (p != CUDA_SUCCESS) { report_.quarantined = true; print(); return; }
                extra_push_ = true;
            }
        }
        bool drained = true;
        if (stream_) {
            const CUresult s = cuStreamSynchronize(stream_);
            issue("cuStreamSynchronize (cleanup)", s);
            if (s != CUDA_SUCCESS) {
                const CUresult c = cuCtxSynchronize();
                issue("cuCtxSynchronize (cleanup)", c);
                drained = c == CUDA_SUCCESS;
            }
        }
        report_.drained = drained;
        if (!drained) {
            // Sticky error: completion of submitted work cannot be proven, so
            // nothing is released. The process must terminate.
            report_.quarantined = true;
            print();
            return;
        }
        for (auto e : execs_) issue("cuGraphExecDestroy", cuGraphExecDestroy(e));
        for (auto g : graphs_) issue("cuGraphDestroy", cuGraphDestroy(g));
        for (auto e : events_) issue("cuEventDestroy", cuEventDestroy(e));
        for (auto p : device_allocs_) { const CUresult r = cuMemFree(p); issue("cuMemFree", r); report_.freed_device += r == CUDA_SUCCESS; }
        for (auto p : pinned_allocs_) { const CUresult r = cuMemFreeHost(p); issue("cuMemFreeHost", r); report_.freed_pinned += r == CUDA_SUCCESS; }
        if (module_) issue("cuModuleUnload", cuModuleUnload(module_));
        if (stream_) issue("cuStreamDestroy", cuStreamDestroy(stream_));
        execs_.clear(); graphs_.clear(); events_.clear(); device_allocs_.clear(); pinned_allocs_.clear();
        module_ = nullptr; stream_ = nullptr;
        bool restored = true;
        if (extra_push_) { CUcontext popped = nullptr; const CUresult r = cuCtxPopCurrent(&popped); issue("cuCtxPopCurrent (cleanup)", r); restored = r == CUDA_SUCCESS && popped == ctx_; }
        if (pushed_) { CUcontext popped = nullptr; const CUresult r = cuCtxPopCurrent(&popped); issue("cuCtxPopCurrent", r); restored = restored && r == CUDA_SUCCESS && popped == ctx_; }
        report_.context_restored = restored;
        if (retained_ && restored) {
            const CUresult r = cuDevicePrimaryCtxRelease(device_);
            issue("cuDevicePrimaryCtxRelease", r);
            report_.primary_released = r == CUDA_SUCCESS;
        }
        print();
    }

private:
    void issue(const char* op, CUresult r) noexcept {
        if (r == CUDA_SUCCESS) return;
        try { report_.issues.push_back({op, int(r)}); } catch (...) {}
    }
    void print() const noexcept {
        for (const auto& i : report_.issues) std::fprintf(stderr, "CUDA cleanup issue: %s (%d)\n", i.operation.c_str(), i.code);
        if (report_.quarantined) std::fprintf(stderr, "CUDA session quarantined: resources retained; terminate this process.\n");
    }
    CleanupReport& report_;
    CUdevice device_ = 0;
    CUcontext ctx_ = nullptr;
    CUstream stream_ = nullptr;
    CUmodule module_ = nullptr;
    bool retained_ = false, pushed_ = false, extra_push_ = false;
    std::vector<CUdeviceptr> device_allocs_;
    std::vector<void*> pinned_allocs_;
    std::vector<CUevent> events_;
    std::vector<CUgraph> graphs_;
    std::vector<CUgraphExec> execs_;
    char* staging_ = nullptr;
    std::size_t staging_bytes_ = 0;
    DeviceInfo info_;
    std::thread::id thread_{};
};
}  // namespace hgf::cu
