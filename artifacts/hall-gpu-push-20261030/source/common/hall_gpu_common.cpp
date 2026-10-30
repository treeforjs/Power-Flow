// Implementation of hall_gpu_common.hpp. Driver API and NVRTC are resolved
// with dlopen/dlsym; cuda.h and nvrtc.h are used for declarations only.
#include "hall_gpu_common.hpp"

#include <cuda.h>
#include <nvrtc.h>
#include <dlfcn.h>
#include <openssl/evp.h>

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <thread>

#define HGC_STR2(x) #x
#define HGC_STR(x) HGC_STR2(x)

namespace hall_gpu {

CudaError::CudaError(std::string call_name, int code, const std::string& detail)
    : Error(call_name + " failed (status " + std::to_string(code) + "): " + detail), status(code), call(std::move(call_name)) {}

std::string sha256_hex(const void* data, std::size_t bytes) {
    unsigned char out[EVP_MAX_MD_SIZE];
    unsigned int length = 0;
    if (EVP_Digest(data, bytes, out, &length, EVP_sha256(), nullptr) != 1) throw Error("EVP_Digest sha256 failed");
    static const char* hex = "0123456789abcdef";
    std::string s;
    s.reserve(2 * length);
    for (unsigned i = 0; i < length; ++i) { s.push_back(hex[out[i] >> 4]); s.push_back(hex[out[i] & 15]); }
    return s;
}

std::string sha256_file(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw Error("cannot read " + path.string());
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    if (!ctx || EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr) != 1) { EVP_MD_CTX_free(ctx); throw Error("sha256 init failed"); }
    std::vector<char> buffer(1 << 20);
    while (in) {
        in.read(buffer.data(), std::streamsize(buffer.size()));
        if (in.gcount() > 0 && EVP_DigestUpdate(ctx, buffer.data(), std::size_t(in.gcount())) != 1) {
            EVP_MD_CTX_free(ctx); throw Error("sha256 update failed");
        }
    }
    unsigned char out[EVP_MAX_MD_SIZE];
    unsigned int length = 0;
    const int ok = EVP_DigestFinal_ex(ctx, out, &length);
    EVP_MD_CTX_free(ctx);
    if (ok != 1) throw Error("sha256 final failed");
    static const char* hex = "0123456789abcdef";
    std::string s;
    for (unsigned i = 0; i < length; ++i) { s.push_back(hex[out[i] >> 4]); s.push_back(hex[out[i] & 15]); }
    return s;
}

namespace {

std::string mapped_path_of(const void* address) {
    std::ifstream maps("/proc/self/maps");
    std::string line;
    const auto a = reinterpret_cast<std::uintptr_t>(address);
    while (std::getline(maps, line)) {
        unsigned long long lo = 0, hi = 0;
        if (std::sscanf(line.c_str(), "%llx-%llx", &lo, &hi) != 2) continue;
        if (a >= lo && a < hi) {
            const auto slash = line.find('/');
            return slash == std::string::npos ? std::string() : line.substr(slash);
        }
    }
    return {};
}

struct DriverApi {
    void* library = nullptr;
#define HGC_DRIVER_FUNCTIONS(X) \
    X(cuInit) X(cuDriverGetVersion) X(cuDeviceGetCount) X(cuDeviceGet) X(cuDeviceGetName) \
    X(cuDeviceGetUuid_v2) X(cuDeviceGetAttribute) X(cuDeviceTotalMem) X(cuDeviceGetPCIBusId) \
    X(cuDevicePrimaryCtxSetFlags) X(cuDevicePrimaryCtxGetState) X(cuDevicePrimaryCtxRetain) \
    X(cuDevicePrimaryCtxRelease) X(cuCtxPushCurrent) X(cuCtxPopCurrent) X(cuCtxGetCurrent) \
    X(cuCtxSynchronize) X(cuStreamCreate) X(cuStreamDestroy) X(cuStreamSynchronize) \
    X(cuEventCreate) X(cuEventDestroy) X(cuEventRecord) X(cuEventSynchronize) X(cuEventElapsedTime) X(cuMemAlloc) \
    X(cuMemFree) X(cuMemHostAlloc) X(cuMemFreeHost) X(cuMemcpyHtoDAsync) X(cuMemcpyDtoHAsync) \
    X(cuMemsetD8Async) X(cuMemGetInfo) X(cuModuleLoadDataEx) X(cuModuleUnload) \
    X(cuModuleGetFunction) X(cuLaunchKernel) X(cuGetErrorString) X(cuGetErrorName)
#define HGC_DECLARE(name) decltype(&name) name##_ = nullptr;
    HGC_DRIVER_FUNCTIONS(HGC_DECLARE)
#undef HGC_DECLARE
    void load() {
        library = dlopen("libcuda.so.1", RTLD_NOW | RTLD_LOCAL);
        if (!library) throw Error(std::string("dlopen libcuda.so.1 failed: ") + dlerror());
#define HGC_LOAD(name)                                                                    \
        name##_ = reinterpret_cast<decltype(name##_)>(dlsym(library, HGC_STR(name)));     \
        if (!name##_) throw Error(std::string("libcuda.so.1 lacks ") + HGC_STR(name));
        HGC_DRIVER_FUNCTIONS(HGC_LOAD)
#undef HGC_LOAD
    }
    std::string describe(CUresult r) const {
        const char* name = nullptr;
        const char* text = nullptr;
        if (cuGetErrorName_) cuGetErrorName_(r, &name);
        if (cuGetErrorString_) cuGetErrorString_(r, &text);
        return std::string(name ? name : "CUDA_ERROR") + ": " + (text ? text : "no description");
    }
};

struct NvrtcApi {
    void* library = nullptr;
    std::filesystem::path path;
#define HGC_NVRTC_FUNCTIONS(X) \
    X(nvrtcVersion) X(nvrtcGetNumSupportedArchs) X(nvrtcGetSupportedArchs) X(nvrtcCreateProgram) \
    X(nvrtcCompileProgram) X(nvrtcGetProgramLogSize) X(nvrtcGetProgramLog) X(nvrtcGetCUBINSize) \
    X(nvrtcGetCUBIN) X(nvrtcGetPTXSize) X(nvrtcGetPTX) X(nvrtcDestroyProgram) X(nvrtcGetErrorString)
#define HGC_DECLARE(name) decltype(&name) name##_ = nullptr;
    HGC_NVRTC_FUNCTIONS(HGC_DECLARE)
#undef HGC_DECLARE
    void load(const std::filesystem::path& p) {
        if (p.empty() || !p.is_absolute()) throw Error("NVRTC library must be given as an absolute path");
        path = p;
        library = dlopen(p.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (!library) throw Error(std::string("dlopen NVRTC failed: ") + dlerror());
#define HGC_LOAD(name)                                                                    \
        name##_ = reinterpret_cast<decltype(name##_)>(dlsym(library, HGC_STR(name)));     \
        if (!name##_) throw Error(std::string("NVRTC library lacks ") + HGC_STR(name));
        HGC_NVRTC_FUNCTIONS(HGC_LOAD)
#undef HGC_LOAD
    }
    void check(nvrtcResult r, const char* call) const {
        if (r != NVRTC_SUCCESS) throw CudaError(call, int(r), nvrtcGetErrorString_ ? nvrtcGetErrorString_(r) : "NVRTC error");
    }
};

std::string format_uuid(const CUuuid& u) {
    const auto* b = reinterpret_cast<const unsigned char*>(u.bytes);
    char s[64];
    std::snprintf(s, sizeof s, "GPU-%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
                  b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]);
    return s;
}

std::string first_line(const char* path) {
    std::ifstream in(path);
    std::string line;
    if (in) std::getline(in, line);
    return line;
}

bool forbidden_option(const std::string& o) {
    for (const char* prefix : {"--fmad", "-fmad", "--prec-div", "-prec-div", "--prec-sqrt", "-prec-sqrt", "--ftz", "-ftz",
                               "--use_fast_math", "-use_fast_math", "--gpu-architecture", "-arch", "--extra-device-vectorization"})
        if (o.rfind(prefix, 0) == 0) return true;
    return false;
}

}  // namespace

struct Session::Impl {
    struct Slot {
        bool device = true, live = false;
        CUdeviceptr dptr = 0;
        void* hptr = nullptr;
        std::size_t bytes = 0;
        std::string tag;
    };
    SessionOptions opt;
    DriverApi d;
    NvrtcApi n;
    bool nvrtc_loaded = false;
    DeviceInfo info;
    CUdevice dev = 0;
    CUcontext ctx = nullptr;
    bool retained = false, pushed = false;
    CUstream stream = nullptr;
    std::vector<CUevent> events;
    std::vector<std::string> labels;
    std::size_t used_events = 0;
    std::vector<Slot> slots;
    std::vector<CUmodule> modules;
    bool failed = false, closed = false;
    std::string first_error;
    std::size_t in_use = 0, peak = 0, pinned = 0, free_at_open = 0, total = 0, min_free = SIZE_MAX;
    TeardownReport report;
    std::thread::id owner = std::this_thread::get_id();

    explicit Impl(const SessionOptions& o) : opt(o) {}

    [[noreturn]] void fail(const char* call, CUresult r) {
        const std::string detail = d.describe(r);
        if (!failed) { failed = true; first_error = std::string(call) + ": " + detail; }
        throw CudaError(call, int(r), detail);
    }
    void check(CUresult r, const char* call) { if (r != CUDA_SUCCESS) fail(call, r); }
    void ready() const {
        if (closed) throw Error("GPU session is closed");
        if (failed) throw Error("GPU session failed earlier: " + first_error);
        if (std::this_thread::get_id() != owner) throw Error("GPU session used from a thread other than its owner");
    }
    void issue(const char* op, CUresult r) noexcept {
        if (r == CUDA_SUCCESS) return;
        if (report.issue_count < report.issues.size()) report.issues[report.issue_count] = {op, int(r)};
        ++report.issue_count;
    }
    Slot& live_slot(int slot, bool device, std::size_t bytes) {
        if (slot < 0 || std::size_t(slot) >= slots.size() || !slots[slot].live || slots[slot].device != device || slots[slot].bytes != bytes)
            throw Error("stale or foreign buffer handle");
        return slots[slot];
    }
    void observe_free() {
        std::size_t f = 0, t = 0;
        check(d.cuMemGetInfo_(&f, &t), "cuMemGetInfo");
        min_free = std::min(min_free, f);
    }

    void open() {
        d.load();
        check(d.cuInit_(0), "cuInit");
        check(d.cuDriverGetVersion_(&info.driver_api_version), "cuDriverGetVersion");
        info.driver_library = mapped_path_of(reinterpret_cast<const void*>(d.cuInit_));
        if (!info.driver_library.empty()) info.driver_library_sha256 = sha256_file(info.driver_library);
        info.kernel_module_version = first_line("/proc/driver/nvidia/version");
        int count = 0;
        check(d.cuDeviceGetCount_(&count), "cuDeviceGetCount");
        if (count <= 0) throw Error("no CUDA device visible");
        int chosen = -1;
        for (int i = 0; i < count; ++i) {
            CUdevice candidate = 0;
            check(d.cuDeviceGet_(&candidate, i), "cuDeviceGet");
            CUuuid u{};
            check(d.cuDeviceGetUuid_v2_(&u, candidate), "cuDeviceGetUuid_v2");
            const auto text = format_uuid(u);
            if (opt.require_uuid.empty() ? i == opt.ordinal : text == opt.require_uuid) { chosen = i; dev = candidate; info.uuid = text; break; }
        }
        if (chosen < 0) throw Error("requested GPU not visible: " + (opt.require_uuid.empty() ? std::to_string(opt.ordinal) : opt.require_uuid));
        info.ordinal = chosen;
        char name[256] = {};
        check(d.cuDeviceGetName_(name, int(sizeof name) - 1, dev), "cuDeviceGetName");
        info.name = name;
        char bus[64] = {};
        check(d.cuDeviceGetPCIBusId_(bus, int(sizeof bus) - 1, dev), "cuDeviceGetPCIBusId");
        info.pci_bus_id = bus;
        check(d.cuDeviceGetAttribute_(&info.cc_major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, dev), "cuDeviceGetAttribute major");
        check(d.cuDeviceGetAttribute_(&info.cc_minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, dev), "cuDeviceGetAttribute minor");
        check(d.cuDeviceGetAttribute_(&info.multiprocessors, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, dev), "cuDeviceGetAttribute SMs");
        check(d.cuDeviceTotalMem_(&info.total_memory, dev), "cuDeviceTotalMem");
        unsigned flags = 0;
        int active = 0;
        check(d.cuDevicePrimaryCtxGetState_(dev, &flags, &active), "cuDevicePrimaryCtxGetState");
        info.primary_context_was_active = active != 0;
        if (opt.blocking_sync && !active) check(d.cuDevicePrimaryCtxSetFlags_(dev, CU_CTX_SCHED_BLOCKING_SYNC), "cuDevicePrimaryCtxSetFlags");
        check(d.cuDevicePrimaryCtxRetain_(&ctx, dev), "cuDevicePrimaryCtxRetain");
        retained = true;
        check(d.cuCtxPushCurrent_(ctx), "cuCtxPushCurrent");
        pushed = true;
        check(d.cuDevicePrimaryCtxGetState_(dev, &info.primary_context_flags, &active), "cuDevicePrimaryCtxGetState");
        if (opt.blocking_sync && (info.primary_context_flags & CU_CTX_SCHED_MASK) != CU_CTX_SCHED_BLOCKING_SYNC)
            throw Error("primary context is active with a scheduling policy other than blocking sync");
        check(d.cuStreamCreate_(&stream, CU_STREAM_NON_BLOCKING), "cuStreamCreate");
        events.assign(opt.event_pool, nullptr);
        labels.assign(opt.event_pool, std::string());
        for (auto& e : events) check(d.cuEventCreate_(&e, CU_EVENT_DEFAULT), "cuEventCreate");
        check(d.cuMemGetInfo_(&free_at_open, &total), "cuMemGetInfo");
        min_free = free_at_open;
    }

    TeardownReport close() noexcept {
        if (closed) return report;
        closed = true;
        report.attempted = true;
        report.first_error = first_error;
        if (!d.library || !retained) { report.released = true; return report; }
        if (std::this_thread::get_id() != owner) {
            issue("teardown on a thread other than the owner", CUDA_ERROR_INVALID_CONTEXT);
            report.quarantined = true;
            std::fprintf(stderr, "hall_gpu: teardown from a foreign thread; resources retained, exit this process\n");
            return report;
        }
        bool own_push = false;
        CUcontext current = nullptr;
        const auto got = d.cuCtxGetCurrent_(&current);
        issue("teardown cuCtxGetCurrent", got);
        if (got != CUDA_SUCCESS || current != ctx) {
            const auto r = d.cuCtxPushCurrent_(ctx);
            issue("teardown cuCtxPushCurrent", r);
            if (r != CUDA_SUCCESS) {
                report.quarantined = true;
                std::fprintf(stderr, "hall_gpu: cannot make the session context current; resources retained, exit this process\n");
                return report;
            }
            own_push = true;
        }
        bool drained = true;
        if (stream) {
            const auto r = d.cuStreamSynchronize_(stream);
            issue("teardown cuStreamSynchronize", r);
            drained = r == CUDA_SUCCESS;
        }
        if (!drained) {
            report.context_fallback_used = true;
            const auto r = d.cuCtxSynchronize_();
            issue("teardown cuCtxSynchronize", r);
            drained = r == CUDA_SUCCESS;
        }
        if (!drained) {
            // Completion of queued copies and kernels cannot be proven. Pinned
            // host pages may still be DMA targets and device memory may still
            // be in use: free nothing, keep the context reference.
            report.quarantined = true;
            std::fprintf(stderr, "hall_gpu: stream drain failed (%s); %zu device and pinned buffers retained, exit this process\n",
                         first_error.c_str(), slots.size());
            return report;
        }
        report.stream_drained = true;
        const std::size_t issues_before = report.issue_count;
        for (auto& e : events) if (e) { issue("cuEventDestroy", d.cuEventDestroy_(e)); e = nullptr; }
        for (auto& m : modules) if (m) { issue("cuModuleUnload", d.cuModuleUnload_(m)); m = nullptr; }
        for (auto& s : slots) {
            if (!s.live) continue;
            if (s.device) issue("cuMemFree", d.cuMemFree_(s.dptr));
            else issue("cuMemFreeHost", d.cuMemFreeHost_(s.hptr));
            s.live = false;
        }
        if (stream) { issue("cuStreamDestroy", d.cuStreamDestroy_(stream)); stream = nullptr; }
        CUcontext popped = nullptr;
        if (own_push) issue("teardown cuCtxPopCurrent", d.cuCtxPopCurrent_(&popped));
        if (pushed) { issue("cuCtxPopCurrent", d.cuCtxPopCurrent_(&popped)); pushed = false; }
        const auto r = d.cuDevicePrimaryCtxRelease_(dev);
        issue("cuDevicePrimaryCtxRelease", r);
        report.context_released = r == CUDA_SUCCESS;
        retained = false;
        report.released = report.issue_count == issues_before;
        return report;
    }
};

Session::Session(const SessionOptions& o) : p_(std::make_unique<Impl>(o)) {
    try {
        p_->open();
    } catch (...) {
        p_->close();
        throw;
    }
}

Session::~Session() {
    if (p_) {
        const auto r = p_->close();
        if (r.quarantined) std::fprintf(stderr, "hall_gpu: session destroyed in quarantine\n");
    }
}

const DeviceInfo& Session::device() const { return p_->info; }
bool Session::failed() const noexcept { return p_->failed; }
const SessionOptions& Session::options() const noexcept { return p_->opt; }

DeviceBuffer Session::alloc(std::size_t bytes, const char* tag) {
    auto& s = *p_;
    s.ready();
    if (bytes == 0) throw Error("zero-byte device allocation requested");
    CUdeviceptr ptr = 0;
    s.check(s.d.cuMemAlloc_(&ptr, bytes), "cuMemAlloc");
    s.slots.push_back({true, true, ptr, nullptr, bytes, tag ? tag : ""});
    s.in_use += bytes;
    s.peak = std::max(s.peak, s.in_use);
    s.observe_free();
    return {std::uint64_t(ptr), bytes, int(s.slots.size() - 1)};
}

HostBuffer Session::alloc_pinned(std::size_t bytes, const char* tag) {
    auto& s = *p_;
    s.ready();
    if (bytes == 0) throw Error("zero-byte pinned allocation requested");
    void* ptr = nullptr;
    s.check(s.d.cuMemHostAlloc_(&ptr, bytes, 0), "cuMemHostAlloc");
    s.slots.push_back({false, true, 0, ptr, bytes, tag ? tag : ""});
    s.pinned += bytes;
    return {ptr, bytes, int(s.slots.size() - 1)};
}

void Session::release(DeviceBuffer& b) {
    auto& s = *p_;
    if (!b) return;
    s.ready();
    auto& slot = s.live_slot(b.slot, true, b.bytes);
    s.check(s.d.cuStreamSynchronize_(s.stream), "cuStreamSynchronize before cuMemFree");
    s.check(s.d.cuMemFree_(slot.dptr), "cuMemFree");
    slot.live = false;
    s.in_use -= slot.bytes;
    b = {};
}

void Session::release(HostBuffer& b) {
    auto& s = *p_;
    if (!b) return;
    s.ready();
    auto& slot = s.live_slot(b.slot, false, b.bytes);
    s.check(s.d.cuStreamSynchronize_(s.stream), "cuStreamSynchronize before cuMemFreeHost");
    s.check(s.d.cuMemFreeHost_(slot.hptr), "cuMemFreeHost");
    slot.live = false;
    s.pinned -= slot.bytes;
    b = {};
}

void Session::upload(const DeviceBuffer& dst, std::size_t dst_offset, const HostBuffer& src, std::size_t src_offset, std::size_t bytes) {
    auto& s = *p_;
    s.ready();
    s.live_slot(dst.slot, true, dst.bytes);
    s.live_slot(src.slot, false, src.bytes);
    if (bytes == 0) return;
    if (dst_offset + bytes > dst.bytes || src_offset + bytes > src.bytes) throw Error("upload range exceeds a buffer");
    s.check(s.d.cuMemcpyHtoDAsync_(CUdeviceptr(dst.address + dst_offset), static_cast<const char*>(src.pointer) + src_offset, bytes, s.stream),
            "cuMemcpyHtoDAsync");
    if (s.opt.check_each_launch) synchronize("after upload");
}

void Session::download(const HostBuffer& dst, std::size_t dst_offset, const DeviceBuffer& src, std::size_t src_offset, std::size_t bytes) {
    auto& s = *p_;
    s.ready();
    s.live_slot(dst.slot, false, dst.bytes);
    s.live_slot(src.slot, true, src.bytes);
    if (bytes == 0) return;
    if (dst_offset + bytes > dst.bytes || src_offset + bytes > src.bytes) throw Error("download range exceeds a buffer");
    s.check(s.d.cuMemcpyDtoHAsync_(static_cast<char*>(dst.pointer) + dst_offset, CUdeviceptr(src.address + src_offset), bytes, s.stream),
            "cuMemcpyDtoHAsync");
    if (s.opt.check_each_launch) synchronize("after download");
}

void Session::zero(const DeviceBuffer& dst, std::size_t offset, std::size_t bytes) {
    auto& s = *p_;
    s.ready();
    s.live_slot(dst.slot, true, dst.bytes);
    if (bytes == 0) return;
    if (offset + bytes > dst.bytes) throw Error("memset range exceeds a buffer");
    s.check(s.d.cuMemsetD8Async_(CUdeviceptr(dst.address + offset), 0, bytes, s.stream), "cuMemsetD8Async");
    if (s.opt.check_each_launch) synchronize("after memset");
}

void Session::launch(const Kernel& k, Dim3 grid, Dim3 block, unsigned shared_bytes, void** args) {
    auto& s = *p_;
    s.ready();
    if (!k.handle) throw Error("launch of an unresolved kernel");
    if (grid.x == 0 || grid.y == 0 || grid.z == 0) return;   // empty problem: nothing to enqueue
    const auto r = s.d.cuLaunchKernel_(static_cast<CUfunction>(k.handle), grid.x, grid.y, grid.z, block.x, block.y, block.z,
                                       shared_bytes, s.stream, args, nullptr);
    if (r != CUDA_SUCCESS) s.fail(("cuLaunchKernel " + k.name).c_str(), r);
    if (s.opt.check_each_launch) synchronize(k.name.c_str());
}

void Session::synchronize(const char* where) {
    auto& s = *p_;
    s.ready();
    const auto r = s.d.cuStreamSynchronize_(s.stream);
    if (r != CUDA_SUCCESS) s.fail((std::string("cuStreamSynchronize at ") + (where ? where : "?")).c_str(), r);
}

Module Session::compile_and_load(const CompileRequest& req) {
    auto& s = *p_;
    s.ready();
    if (!s.nvrtc_loaded) { s.n.load(s.opt.nvrtc_library); s.nvrtc_loaded = true; }
    Module m;
    auto& rec = m.record;
    rec.nvrtc_library = s.opt.nvrtc_library.string();
    rec.nvrtc_library_sha256 = sha256_file(s.opt.nvrtc_library);
    for (const auto& entry : std::filesystem::directory_iterator(s.opt.nvrtc_library.parent_path())) {
        const auto name = entry.path().filename().string();
        if (name.rfind("libnvrtc-builtins.so.", 0) == 0) rec.nvrtc_builtins_sha256 = name + ":" + sha256_file(entry.path());
    }
    s.n.check(s.n.nvrtcVersion_(&rec.nvrtc_major, &rec.nvrtc_minor), "nvrtcVersion");
    int archs = 0;
    s.n.check(s.n.nvrtcGetNumSupportedArchs_(&archs), "nvrtcGetNumSupportedArchs");
    rec.supported_archs.resize(std::size_t(archs));
    if (archs > 0) s.n.check(s.n.nvrtcGetSupportedArchs_(rec.supported_archs.data()), "nvrtcGetSupportedArchs");
    rec.device_cc = s.info.cc_major * 10 + s.info.cc_minor;
    const bool exact = std::find(rec.supported_archs.begin(), rec.supported_archs.end(), rec.device_cc) != rec.supported_archs.end();
    if (exact) {
        rec.target_cc = rec.device_cc;
        rec.arch_option = "--gpu-architecture=sm_" + std::to_string(rec.target_cc);
        rec.image_kind = "cubin";
    } else {
        int best = -1;
        for (int a : rec.supported_archs) if (a <= rec.device_cc && a > best) best = a;
        if (best < 0) throw Error("NVRTC supports no architecture at or below the device compute capability " + std::to_string(rec.device_cc));
        rec.target_cc = best;
        rec.arch_option = "--gpu-architecture=compute_" + std::to_string(best);
        rec.image_kind = "ptx-jit";
    }
    rec.options = {rec.arch_option, "--fmad=false", "--prec-div=true", "--prec-sqrt=true", "--ftz=false", "--std=c++17"};
    for (const auto& o : req.extra_options) {
        if (forbidden_option(o)) throw Error("compile option overrides the deterministic arithmetic or architecture: " + o);
        rec.options.push_back(o);
    }
    rec.source_sha256.push_back({req.main.name, sha256_hex(req.main.text.data(), req.main.text.size())});
    std::vector<const char*> header_text, header_name;
    for (const auto& h : req.headers) {
        header_text.push_back(h.text.c_str());
        header_name.push_back(h.name.c_str());
        rec.source_sha256.push_back({h.name, sha256_hex(h.text.data(), h.text.size())});
    }
    nvrtcProgram program = nullptr;
    s.n.check(s.n.nvrtcCreateProgram_(&program, req.main.text.c_str(), req.main.name.c_str(), int(header_text.size()),
                                      header_text.data(), header_name.data()), "nvrtcCreateProgram");
    std::vector<const char*> options;
    for (const auto& o : rec.options) options.push_back(o.c_str());
    const auto start = std::chrono::steady_clock::now();
    const auto compiled = s.n.nvrtcCompileProgram_(program, int(options.size()), options.data());
    rec.compile_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    std::size_t log_size = 0;
    if (s.n.nvrtcGetProgramLogSize_(program, &log_size) == NVRTC_SUCCESS && log_size > 1) {
        std::string log(log_size, '\0');
        if (s.n.nvrtcGetProgramLog_(program, log.data()) == NVRTC_SUCCESS) rec.log = log.c_str();
    }
    try {
        s.n.check(compiled, "nvrtcCompileProgram");
        std::size_t size = 0;
        if (rec.image_kind == "cubin") {
            s.n.check(s.n.nvrtcGetCUBINSize_(program, &size), "nvrtcGetCUBINSize");
            m.image.resize(size);
            s.n.check(s.n.nvrtcGetCUBIN_(program, m.image.data()), "nvrtcGetCUBIN");
        } else {
            s.n.check(s.n.nvrtcGetPTXSize_(program, &size), "nvrtcGetPTXSize");
            m.image.resize(size);
            s.n.check(s.n.nvrtcGetPTX_(program, m.image.data()), "nvrtcGetPTX");
        }
    } catch (const CudaError& e) {
        s.n.nvrtcDestroyProgram_(&program);
        throw Error(std::string(e.what()) + "\n" + rec.log);
    }
    s.n.nvrtcDestroyProgram_(&program);
    rec.image_sha256 = sha256_hex(m.image.data(), m.image.size());
    CUmodule module = nullptr;
    s.check(s.d.cuModuleLoadDataEx_(&module, m.image.data(), 0, nullptr, nullptr), "cuModuleLoadDataEx");
    s.modules.push_back(module);
    m.handle = module;
    return m;
}

Kernel Session::kernel(const Module& m, const char* name) {
    auto& s = *p_;
    s.ready();
    CUfunction f = nullptr;
    s.check(s.d.cuModuleGetFunction_(&f, static_cast<CUmodule>(m.handle), name), "cuModuleGetFunction");
    return {f, name};
}

int Session::mark(const char* label) {
    auto& s = *p_;
    s.ready();
    if (s.used_events >= s.events.size()) throw Error("timing marker pool exhausted; clear_markers() after synchronize()");
    s.check(s.d.cuEventRecord_(s.events[s.used_events], s.stream), "cuEventRecord");
    s.labels[s.used_events] = label ? label : "";
    return int(s.used_events++);
}

void Session::wait_marker(int m, const char* where) {
    auto& s = *p_;
    s.ready();
    if (m < 0 || std::size_t(m) >= s.used_events) throw Error("unknown timing marker");
    const auto r = s.d.cuEventSynchronize_(s.events[std::size_t(m)]);
    if (r != CUDA_SUCCESS) s.fail((std::string("cuEventSynchronize at ") + (where ? where : "?")).c_str(), r);
}

double Session::elapsed_ms(int from, int to) const {
    auto& s = *p_;
    s.ready();
    if (from < 0 || to < 0 || std::size_t(from) >= s.used_events || std::size_t(to) >= s.used_events) throw Error("unknown timing marker");
    float ms = 0;
    const auto r = s.d.cuEventElapsedTime_(&ms, s.events[std::size_t(from)], s.events[std::size_t(to)]);
    if (r != CUDA_SUCCESS) s.fail("cuEventElapsedTime", r);
    return double(ms);
}

std::size_t Session::marker_count() const noexcept { return p_->used_events; }
const std::string& Session::marker_label(int i) const { return p_->labels.at(std::size_t(i)); }
void Session::clear_markers() { p_->ready(); p_->used_events = 0; }

MemoryReport Session::memory() {
    auto& s = *p_;
    s.ready();
    MemoryReport r;
    std::size_t f = 0, t = 0;
    s.check(s.d.cuMemGetInfo_(&f, &t), "cuMemGetInfo");
    s.min_free = std::min(s.min_free, f);
    r.own_bytes_in_use = s.in_use;
    r.own_peak_bytes = s.peak;
    r.own_pinned_bytes = s.pinned;
    r.device_free_at_open = s.free_at_open;
    r.device_total = t;
    r.device_min_free_observed = s.min_free;
    r.device_free_now = f;
    return r;
}

TeardownReport Session::close() noexcept { return p_->close(); }

std::string json_string(const std::string& s) {
    std::string o = "\"";
    for (unsigned char c : s) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\r': o += "\\r"; break;
            case '\t': o += "\\t"; break;
            default:
                if (c < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); o += b; }
                else o.push_back(char(c));
        }
    }
    return o + "\"";
}

std::string json_number(double v) {
    if (!std::isfinite(v)) return "null";
    char b[64];
    std::snprintf(b, sizeof b, "%.17g", v);
    return b;
}

std::string device_json(const DeviceInfo& i) {
    std::ostringstream o;
    o << "{\"ordinal\":" << i.ordinal << ",\"name\":" << json_string(i.name) << ",\"uuid\":" << json_string(i.uuid)
      << ",\"pci_bus_id\":" << json_string(i.pci_bus_id) << ",\"compute_capability\":" << json_string(std::to_string(i.cc_major) + "." + std::to_string(i.cc_minor))
      << ",\"multiprocessors\":" << i.multiprocessors << ",\"driver_api_version\":" << i.driver_api_version
      << ",\"kernel_module_version\":" << json_string(i.kernel_module_version) << ",\"total_memory_bytes\":" << i.total_memory
      << ",\"driver_library\":" << json_string(i.driver_library) << ",\"driver_library_sha256\":" << json_string(i.driver_library_sha256)
      << ",\"primary_context_flags\":" << i.primary_context_flags << ",\"primary_context_was_active\":" << (i.primary_context_was_active ? "true" : "false")
      << "}";
    return o.str();
}

std::string compile_json(const CompileRecord& r) {
    std::ostringstream o;
    o << "{\"nvrtc_version\":" << json_string(std::to_string(r.nvrtc_major) + "." + std::to_string(r.nvrtc_minor))
      << ",\"nvrtc_library\":" << json_string(r.nvrtc_library) << ",\"nvrtc_library_sha256\":" << json_string(r.nvrtc_library_sha256)
      << ",\"nvrtc_builtins\":" << json_string(r.nvrtc_builtins_sha256) << ",\"supported_archs\":[";
    for (std::size_t i = 0; i < r.supported_archs.size(); ++i) o << (i ? "," : "") << r.supported_archs[i];
    o << "],\"device_cc\":" << r.device_cc << ",\"target_cc\":" << r.target_cc << ",\"arch_option\":" << json_string(r.arch_option)
      << ",\"image_kind\":" << json_string(r.image_kind) << ",\"options\":[";
    for (std::size_t i = 0; i < r.options.size(); ++i) o << (i ? "," : "") << json_string(r.options[i]);
    o << "],\"sources\":{";
    for (std::size_t i = 0; i < r.source_sha256.size(); ++i)
        o << (i ? "," : "") << json_string(r.source_sha256[i].first) << ":" << json_string(r.source_sha256[i].second);
    o << "},\"image_sha256\":" << json_string(r.image_sha256) << ",\"compile_ms\":" << json_number(r.compile_ms)
      << ",\"log\":" << json_string(r.log) << "}";
    return o.str();
}

std::string memory_json(const MemoryReport& m) {
    std::ostringstream o;
    o << "{\"own_bytes_in_use\":" << m.own_bytes_in_use << ",\"own_peak_device_bytes\":" << m.own_peak_bytes
      << ",\"own_pinned_host_bytes\":" << m.own_pinned_bytes << ",\"device_free_at_open\":" << m.device_free_at_open
      << ",\"device_total\":" << m.device_total << ",\"device_min_free_observed\":" << m.device_min_free_observed
      << ",\"device_free_now\":" << m.device_free_now
      << ",\"device_peak_used_by_all_processes\":" << (m.device_total - m.device_min_free_observed) << "}";
    return o.str();
}

std::string teardown_json(const TeardownReport& r) {
    std::ostringstream o;
    o << "{\"attempted\":" << (r.attempted ? "true" : "false") << ",\"stream_drained\":" << (r.stream_drained ? "true" : "false")
      << ",\"context_fallback_used\":" << (r.context_fallback_used ? "true" : "false") << ",\"released\":" << (r.released ? "true" : "false")
      << ",\"quarantined\":" << (r.quarantined ? "true" : "false") << ",\"context_released\":" << (r.context_released ? "true" : "false")
      << ",\"first_error\":" << json_string(r.first_error) << ",\"issues\":[";
    for (std::size_t i = 0; i < r.issue_count && i < r.issues.size(); ++i)
        o << (i ? "," : "") << "{\"operation\":" << json_string(r.issues[i].operation ? r.issues[i].operation : "?") << ",\"status\":" << r.issues[i].status << "}";
    o << "],\"issue_count\":" << r.issue_count << "}";
    return o.str();
}

}  // namespace hall_gpu
