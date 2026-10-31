// nvrtc_archs.cpp -- compile the owner's NVRTC program for several target architectures without a device
// (cubin), with the fixed floating-point options of the session layer, and record the image digests.
//   nvrtc_archs OUT_DIR NVRTC_LIBRARY arch...      e.g. 86 89 120
#include "check.hpp"
#include <dlfcn.h>
#include <nvrtc.h>
#include <openssl/sha.h>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
static const char kKernels[] = {
#embed "../source/phase_kernels.h"
    , 0};
static const char kLaw[] = {
#embed "../source/phase_law.h"
    , 0};
static const char kMath[] = {
#embed "../source/det_math.h"
    , 0};
static const char kDigest[] = {
#embed "../source/state_digest.h"
    , 0};
#pragma GCC diagnostic pop
// identical to the main program of gpu_engine.cpp
static const char kMain[] =
    "#include \"phase_kernels.h\"\n"
    "#include \"state_digest.h\"\n"
    "extern \"C\" __global__ void k_copy_i32(const int* a, int* b, int n) {\n"
    "    const int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) b[i] = a[i];\n"
    "}\n"
    "extern \"C\" __global__ void k_chunks(int section, pq_u64 n, const pq_i128* E, const pq_i128* N, const int* flags, const pq_i128* fcm, pq_u32* out) {\n"
    "    const pq_u64 j = (pq_u64)blockIdx.x * blockDim.x + threadIdx.x; if (j < (n + 31) / 32) pq_chunk_digest(section, j, n, E, N, flags, fcm, out + 8 * j);\n"
    "}\n"
    "extern \"C\" __global__ void k_leaves(int section, pq_u64 n, const pq_u32* chunks, pq_u32* out) {\n"
    "    const pq_u64 k = (pq_u64)blockIdx.x * blockDim.x + threadIdx.x; if (k < (n + 1023) / 1024) pq_leaf_digest(section, k, n, chunks, out + 8 * k);\n"
    "}\n";

static std::string hexdigest(const void* p, std::size_t n) {
    unsigned char d[32]; SHA256(static_cast<const unsigned char*>(p), n, d);
    std::ostringstream o; for (auto c : d) o << std::hex << std::setw(2) << std::setfill('0') << unsigned(c); return o.str();
}

int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: nvrtc_archs OUT NVRTC arch...\n"); return 2; }
    void* h = dlopen(argv[2], RTLD_NOW | RTLD_LOCAL);
    if (!h) { std::fprintf(stderr, "dlopen failed: %s\n", dlerror()); return 2; }
    auto sym = [&](const char* n) { void* f = dlsym(h, n); if (!f) { std::fprintf(stderr, "missing %s\n", n); std::exit(2); } return f; };
    auto create = reinterpret_cast<nvrtcResult (*)(nvrtcProgram*, const char*, const char*, int, const char* const*, const char* const*)>(sym("nvrtcCreateProgram"));
    auto compile = reinterpret_cast<nvrtcResult (*)(nvrtcProgram, int, const char* const*)>(sym("nvrtcCompileProgram"));
    auto logsize = reinterpret_cast<nvrtcResult (*)(nvrtcProgram, size_t*)>(sym("nvrtcGetProgramLogSize"));
    auto getlog = reinterpret_cast<nvrtcResult (*)(nvrtcProgram, char*)>(sym("nvrtcGetProgramLog"));
    auto cubinsize = reinterpret_cast<nvrtcResult (*)(nvrtcProgram, size_t*)>(sym("nvrtcGetCUBINSize"));
    auto getcubin = reinterpret_cast<nvrtcResult (*)(nvrtcProgram, char*)>(sym("nvrtcGetCUBIN"));
    auto destroy = reinterpret_cast<nvrtcResult (*)(nvrtcProgram*)>(sym("nvrtcDestroyProgram"));
    auto version = reinterpret_cast<nvrtcResult (*)(int*, int*)>(sym("nvrtcVersion"));
    int maj = 0, min = 0; version(&maj, &min);
    chk::put_text("nvrtc_version", std::to_string(maj) + "." + std::to_string(min));
    const char* headers[] = {kKernels, kLaw, kMath, kDigest};
    const char* names[] = {"phase_kernels.h", "phase_law.h", "det_math.h", "state_digest.h"};
    std::ostringstream rec; rec << "[";
    for (int a = 3; a < argc; ++a) {
        const std::string arch = std::string("--gpu-architecture=sm_") + argv[a];
        const char* opts[] = {arch.c_str(), "--fmad=false", "--prec-div=true", "--prec-sqrt=true", "--ftz=false", "--std=c++17"};
        nvrtcProgram prog;
        create(&prog, kMain, "phase_owner_kernels.cu", 4, headers, names);
        const nvrtcResult r = compile(prog, 6, opts);
        size_t ls = 0; logsize(prog, &ls); std::string log(ls, '\0'); if (ls) getlog(prog, log.data());
        size_t cs = 0; std::string digest;
        if (r == NVRTC_SUCCESS) { cubinsize(prog, &cs); std::vector<char> img(cs); getcubin(prog, img.data()); digest = hexdigest(img.data(), cs); }
        destroy(&prog);
        chk::check(r == NVRTC_SUCCESS, std::string("NVRTC compiles the owner kernels for sm_") + argv[a]);
        rec << (a > 3 ? "," : "") << "{\"arch\":\"sm_" << argv[a] << "\",\"status\":" << int(r) << ",\"cubin_bytes\":" << cs << ",\"cubin_sha256\":" << chk::str(digest)
            << ",\"log\":" << chk::str(log.substr(0, 2000)) << "}";
    }
    rec << "]";
    chk::put_raw("architectures", rec.str());
    chk::put_text("source_sha256", hexdigest(kMain, sizeof kMain - 1) + "/" + hexdigest(kKernels, sizeof kKernels - 1) + "/" + hexdigest(kLaw, sizeof kLaw - 1) + "/" +
                  hexdigest(kMath, sizeof kMath - 1) + "/" + hexdigest(kDigest, sizeof kDigest - 1));
    return chk::finish("nvrtc_archs", argv[1]);
}
