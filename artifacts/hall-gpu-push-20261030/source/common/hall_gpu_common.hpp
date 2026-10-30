// hall_gpu_common.hpp -- device session, NVRTC compilation and stream-ordered
// memory for the Hall GPU-resident step. Shared by push/ and (later) field/.
//
// Contract summary (details next to each declaration):
//  * One Session = one device (selected by UUID or ordinal), its retained
//    primary context, ONE non-blocking stream and a pool of timing events.
//    Every upload, download, memset and kernel launch is enqueued on that
//    stream, so data written by an upload is visible to every later launch
//    by stream order alone. Nothing uses the legacy default stream or
//    synchronous (pageable) copies.
//  * Host staging memory is page-locked (HostBuffer). The host may write a
//    HostBuffer, or read a downloaded one, only after a Session::synchronize()
//    that follows the last queued operation touching it.
//  * Every driver call is checked synchronously. The first failure marks the
//    session failed and throws CudaError; later calls refuse with Error.
//    Asynchronous kernel faults surface at the next synchronize().
//  * Teardown (close() or the destructor) drains the stream before releasing
//    anything. If the drain fails, a context-wide synchronize is attempted.
//    If completion still cannot be established, nothing is freed (pinned
//    memory may still be a DMA target): the resources are quarantined, the
//    report says so, and the owning process must exit non-zero.
//  * NVRTC: compile() selects --gpu-architecture from the actual device
//    compute capability (sm_XY cubin when NVRTC supports it, otherwise the
//    highest supported compute_XY PTX for driver JIT) and always adds
//    --fmad=false --prec-div=true --prec-sqrt=true --ftz=false. Requests that
//    try to override those options are refused. The record keeps the NVRTC
//    version, library hashes, options, log, source and image hashes.
//
// The header exposes no CUDA types; the implementation loads the driver
// (libcuda.so.1) and the pinned NVRTC library with dlopen at run time.
#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace hall_gpu {

struct Error : std::runtime_error {
    using std::runtime_error::runtime_error;
};
// A failed driver or NVRTC call. status is the CUresult / nvrtcResult value.
struct CudaError : Error {
    int status = 0;
    std::string call;
    CudaError(std::string call_name, int code, const std::string& detail);
};

std::string sha256_hex(const void* data, std::size_t bytes);
std::string sha256_file(const std::filesystem::path& path);

struct DeviceInfo {
    int ordinal = -1;
    std::string name, uuid, pci_bus_id;   // uuid formatted as GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
    int cc_major = 0, cc_minor = 0, multiprocessors = 0;
    int driver_api_version = 0;           // cuDriverGetVersion
    std::string kernel_module_version;    // first line of /proc/driver/nvidia/version, if readable
    std::size_t total_memory = 0;
    std::string driver_library, driver_library_sha256;   // resolved path of the loaded libcuda.so.1
    unsigned primary_context_flags = 0;
    bool primary_context_was_active = false;
};

struct SessionOptions {
    std::string require_uuid;              // if non-empty the device must carry this UUID
    int ordinal = 0;                       // used only when require_uuid is empty
    std::filesystem::path nvrtc_library;   // absolute path of libnvrtc.so.12 (required for compile())
    bool blocking_sync = true;             // CU_CTX_SCHED_BLOCKING_SYNC: host waits yield the core
    bool check_each_launch = false;        // synchronize after every launch (fault localisation)
    std::size_t event_pool = 96;           // timing markers available between clear_markers()
};

struct Dim3 { unsigned x = 1, y = 1, z = 1; };

// Opaque handles. Copies are views; ownership stays with the Session, which
// releases everything at close() unless released earlier with release().
struct DeviceBuffer {
    std::uint64_t address = 0;
    std::size_t bytes = 0;
    int slot = -1;
    template<class T> T* as() const noexcept { return reinterpret_cast<T*>(address); }
    explicit operator bool() const noexcept { return address != 0; }
};
struct HostBuffer {
    void* pointer = nullptr;
    std::size_t bytes = 0;
    int slot = -1;
    template<class T> T* as() const noexcept { return static_cast<T*>(pointer); }
    template<class T> std::span<T> span(std::size_t count, std::size_t offset_elems = 0) const {
        if ((offset_elems + count) * sizeof(T) > bytes) throw Error("HostBuffer span exceeds allocation");
        return {static_cast<T*>(pointer) + offset_elems, count};
    }
    explicit operator bool() const noexcept { return pointer != nullptr; }
};
struct Kernel {
    void* handle = nullptr;
    std::string name;
};

struct ProgramSource {
    std::string name;   // file name as seen by #include
    std::string text;
};
struct CompileRequest {
    ProgramSource main;
    std::vector<ProgramSource> headers;
    std::vector<std::string> extra_options;   // must not touch arch or floating-point options
};
struct CompileRecord {
    int nvrtc_major = 0, nvrtc_minor = 0;
    std::string nvrtc_library, nvrtc_library_sha256, nvrtc_builtins_sha256;
    std::vector<int> supported_archs;
    int device_cc = 0, target_cc = 0;
    std::string arch_option, image_kind;   // image_kind: "cubin" or "ptx-jit"
    std::vector<std::string> options;
    std::string log;
    std::vector<std::pair<std::string, std::string>> source_sha256;   // name -> sha256
    std::string image_sha256;
    double compile_ms = 0;
};
struct Module {
    void* handle = nullptr;
    CompileRecord record;
    std::vector<char> image;
};

struct MemoryReport {
    std::size_t own_bytes_in_use = 0, own_peak_bytes = 0, own_pinned_bytes = 0;
    std::size_t device_free_at_open = 0, device_total = 0, device_min_free_observed = 0;
    std::size_t device_free_now = 0;
};
struct TeardownIssue { const char* operation = nullptr; int status = 0; };
struct TeardownReport {
    bool attempted = false, stream_drained = false, context_fallback_used = false;
    bool released = false, quarantined = false, context_released = false;
    std::string first_error;   // first error that failed the session, if any
    std::array<TeardownIssue, 32> issues{};
    std::size_t issue_count = 0;
};

class Session {
public:
    explicit Session(const SessionOptions&);
    ~Session();
    Session(const Session&) = delete;
    Session& operator=(const Session&) = delete;

    const DeviceInfo& device() const;
    bool failed() const noexcept;
    const SessionOptions& options() const noexcept;

    // Allocation. Device allocations count toward own_bytes_in_use/peak.
    DeviceBuffer alloc(std::size_t bytes, const char* tag);
    HostBuffer alloc_pinned(std::size_t bytes, const char* tag);
    // Drains the stream, then frees. The handle is cleared.
    void release(DeviceBuffer&);
    void release(HostBuffer&);

    // Stream-ordered operations.
    void upload(const DeviceBuffer& dst, std::size_t dst_offset, const HostBuffer& src,
                std::size_t src_offset, std::size_t bytes);
    void download(const HostBuffer& dst, std::size_t dst_offset, const DeviceBuffer& src,
                  std::size_t src_offset, std::size_t bytes);
    void zero(const DeviceBuffer& dst, std::size_t offset, std::size_t bytes);
    void launch(const Kernel&, Dim3 grid, Dim3 block, unsigned shared_bytes, void** args);
    // Blocks until all queued work has completed; surfaces asynchronous faults.
    void synchronize(const char* where);

    // NVRTC compilation for this device and module loading.
    Module compile_and_load(const CompileRequest&);
    Kernel kernel(const Module&, const char* name);

    // Timing markers (events recorded on the stream). Valid after synchronize().
    int mark(const char* label);
    // Blocks until the stream has passed marker m (work queued after it may continue).
    void wait_marker(int m, const char* where);
    double elapsed_ms(int from_marker, int to_marker) const;
    std::size_t marker_count() const noexcept;
    const std::string& marker_label(int) const;
    void clear_markers();   // requires a preceding synchronize()

    MemoryReport memory();
    TeardownReport close() noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> p_;
};

// Minimal JSON text helpers for receipts (no external dependency).
std::string json_string(const std::string&);
std::string json_number(double);   // %.17g, "null" for non-finite
std::string device_json(const DeviceInfo&);
std::string compile_json(const CompileRecord&);
std::string memory_json(const MemoryReport&);
std::string teardown_json(const TeardownReport&);

}  // namespace hall_gpu
