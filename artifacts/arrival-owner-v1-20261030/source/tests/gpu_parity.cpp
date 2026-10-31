// gpu_parity.cpp -- host/device bit parity of the arrival owner v1 RNG layer.
// The host computes every value with include/pf_arrival/philox_core.h; the device
// computes the same values from the same header bytes (embedded at build time and
// compiled by NVRTC through the shared GPU layer). Per key record: one Philox
// block (4 words), 3 conversions, 4 electron marks, 6 neutral marks, one
// stochastic rounding and one exact Bernoulli draw.
// usage: gpu_parity OUTPUT_DIRECTORY NVRTC_LIBRARY DEVICE_UUID RECORDS REPEATS
// exit: 0 all bits equal, 4 parity mismatch, 2 in-process repeat mismatch, 3 teardown quarantine, 1 error.
#include "common/hall_gpu_common.hpp"
#include "pf_arrival/philox_core.h"

#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kHeader[] = {
#embed "../include/pf_arrival/philox_core.h"
    , 0};
const char kKernel[] = {
#embed "../device/arrival_parity.cu"
    , 0};
#pragma GCC diagnostic pop

struct KeyRecord {
    std::uint64_t seed, step;
    std::uint32_t stream, candidate, purpose, block;
};
static_assert(sizeof(KeyRecord) == 32 && alignof(KeyRecord) == 8);

std::uint64_t splitmix(std::uint64_t& s) {
    std::uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}
double from_bits(std::uint64_t b) { return std::bit_cast<double>(b); }
// Rounding inputs in [0, 2^52) of several classes.
double make_x(std::uint64_t i, std::uint64_t& s) {
    const auto r = splitmix(s);
    switch (i % 8) {
        case 0: return double(r >> 11) * 0x1p-41;                                        // [0, 4096)
        case 1: return from_bits(r & 0x000FFFFFFFFFFFFFull);                            // subnormal
        case 2: return from_bits((r & 0x000FFFFFFFFFFFFFull) | (std::uint64_t(1 + r % 120) << 52));   // tiny normal
        case 3: { const double n = double(1 + (r >> 52)); return (r & 1) ? std::nextafter(n, 0.) : std::nextafter(n, 1e300); }
        case 4: return double(r >> 40);                                                  // integers
        case 5: return double(r >> 11) * 0x1p-53;                                        // [0, 1)
        case 6: return double(r % 1024) * 0x1p-10;                                       // dyadic
        default: return double(r >> 12);                                                 // < 2^52
    }
}
double make_p(std::uint64_t i, std::uint64_t& s) {
    const auto r = splitmix(s);
    switch (i % 4) {
        case 0: return double(r >> 11) * 0x1p-53;
        case 1: return from_bits(r & 0x000FFFFFFFFFFFFFull);
        case 2: return from_bits((r & 0x000FFFFFFFFFFFFFull) | (std::uint64_t(1 + r % 1022) << 52));
        default: return double(r % 8) * 0.125;
    }
}
struct Outputs {
    std::vector<std::uint32_t> raw;
    std::vector<double> uniforms, electron, neutral;
    std::vector<std::uint64_t> rounded;
    std::vector<std::uint32_t> bernoulli;
    void resize(std::size_t n) {
        raw.assign(4 * n, 0);
        uniforms.assign(3 * n, 0);
        electron.assign(4 * n, 0);
        neutral.assign(6 * n, 0);
        rounded.assign(n, 0);
        bernoulli.assign(n, 0);
    }
};
template <class T>
std::string sha(const std::vector<T>& v) { return hall_gpu::sha256_hex(v.data(), v.size() * sizeof(T)); }
template <class T>
std::size_t mismatches(const std::vector<T>& a, const std::vector<T>& b) {
    std::size_t n = 0;
    for (std::size_t i = 0; i < a.size(); ++i) n += std::memcmp(&a[i], &b[i], sizeof(T)) != 0;
    return n;
}
double ms_since(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count(); }
}  // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "usage: gpu_parity OUTPUT_DIRECTORY NVRTC_LIBRARY DEVICE_UUID RECORDS REPEATS\n";
        return 1;
    }
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    const std::size_t n = std::stoull(argv[4]);
    const int repeats = std::stoi(argv[5]);
    std::ostringstream js;
    int rc = 0;
    hall_gpu::TeardownReport teardown;
    try {
        if (n == 0 || n > (1u << 26) || repeats < 1) throw std::runtime_error("records must be 1..2^26, repeats >= 1");
        // Inputs (deterministic), including the domain edges.
        std::vector<KeyRecord> keys(n);
        std::vector<double> x(n), p(n);
        std::uint64_t s = 0xA11CE5EEDull;
        for (std::size_t i = 0; i < n; ++i) {
            keys[i] = {splitmix(s), splitmix(s) & (PFA_STEP_LIMIT - 1), std::uint32_t(splitmix(s) & 0xFFFF), std::uint32_t(splitmix(s)),
                       std::uint32_t(1 + splitmix(s) % 255), std::uint32_t(splitmix(s))};
            x[i] = make_x(i, s);
            p[i] = make_p(i, s);
        }
        keys[0] = {0, 0, 0, 0, 1, 0};
        keys[1 % n] = {~0ull, PFA_STEP_LIMIT - 1, 0xFFFF, 0xFFFFFFFFu, 255, 0xFFFFFFFFu};
        // Host reference.
        const auto th = std::chrono::steady_clock::now();
        Outputs host;
        host.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& k = keys[i];
            const pfa::StreamKey sk{k.seed, k.step, k.stream, k.candidate, k.purpose};
            const auto b = pfa::stream_block(sk, k.block);
            for (int j = 0; j < 4; ++j) host.raw[4 * i + j] = b.v[j];
            host.uniforms[3 * i] = pfa::u01_open(pfa::low_word(b));
            host.uniforms[3 * i + 1] = pfa::u01_open(pfa::high_word(b));
            host.uniforms[3 * i + 2] = pfa::u01_half_open(pfa::low_word(b));
            pfa::electron_marks(sk, &host.electron[4 * i]);
            pfa::neutral_marks(sk, &host.neutral[6 * i]);
            host.rounded[i] = pfa::stochastic_round(sk, x[i]);
            host.bernoulli[i] = std::uint32_t(pfa::bernoulli_exact(sk, p[i], k.block));
        }
        const double host_ms = ms_since(th);
        // Device.
        hall_gpu::SessionOptions so;
        so.require_uuid = argv[3];
        so.nvrtc_library = argv[2];
        hall_gpu::Session session(so);
        hall_gpu::CompileRequest req;
        req.main = {"arrival_parity.cu", kKernel};
        req.headers = {{"philox_core.h", kHeader}};
        const auto module = session.compile_and_load(req);
        const auto kernel = session.kernel(module, "pfa_parity");
        const std::size_t in_bytes = n * (sizeof(KeyRecord) + 16), out_bytes = n * (16 + 24 + 32 + 48 + 8 + 4);
        auto d_in = session.alloc(in_bytes, "inputs");
        auto d_out = session.alloc(out_bytes, "outputs");
        auto h_in = session.alloc_pinned(in_bytes, "input staging");
        auto h_out = session.alloc_pinned(out_bytes, "output staging");
        std::memcpy(h_in.as<char>(), keys.data(), n * sizeof(KeyRecord));
        std::memcpy(h_in.as<char>() + n * sizeof(KeyRecord), x.data(), n * 8);
        std::memcpy(h_in.as<char>() + n * (sizeof(KeyRecord) + 8), p.data(), n * 8);
        std::uint64_t a_keys = d_in.address, a_x = a_keys + n * sizeof(KeyRecord), a_p = a_x + n * 8;
        std::uint64_t a_raw = d_out.address, a_uni = a_raw + 16 * n, a_em = a_uni + 24 * n, a_nm = a_em + 32 * n, a_sr = a_nm + 48 * n, a_b = a_sr + 8 * n;
        unsigned count = unsigned(n);
        std::vector<Outputs> runs;
        std::ostringstream timing;
        for (int r = 0; r < repeats; ++r) {
            const auto t0 = std::chrono::steady_clock::now();
            session.upload(d_in, 0, h_in, 0, in_bytes);
            session.zero(d_out, 0, out_bytes);
            void* args[] = {&a_keys, &a_x, &a_p, &count, &a_raw, &a_uni, &a_em, &a_nm, &a_sr, &a_b};
            const int m0 = session.mark("kernel-begin");
            session.launch(kernel, {unsigned((n + 127) / 128)}, {128}, 0, args);
            const int m1 = session.mark("kernel-end");
            session.download(h_out, 0, d_out, 0, out_bytes);
            session.synchronize("parity download");
            Outputs dev;
            dev.resize(n);
            const char* o = h_out.as<char>();
            std::memcpy(dev.raw.data(), o, 16 * n);
            std::memcpy(dev.uniforms.data(), o + 16 * n, 24 * n);
            std::memcpy(dev.electron.data(), o + 40 * n, 32 * n);
            std::memcpy(dev.neutral.data(), o + 72 * n, 48 * n);
            std::memcpy(dev.rounded.data(), o + 120 * n, 8 * n);
            std::memcpy(dev.bernoulli.data(), o + 128 * n, 4 * n);
            timing << (r ? "," : "") << "{\"repeat\":" << r << ",\"kernel_ms\":" << hall_gpu::json_number(session.elapsed_ms(m0, m1))
                   << ",\"wall_ms\":" << hall_gpu::json_number(ms_since(t0)) << "}";
            session.clear_markers();
            runs.push_back(std::move(dev));
        }
        const auto& dev = runs.front();
        bool repeat_identical = true;
        for (const auto& rr : runs)
            repeat_identical = repeat_identical && rr.raw == dev.raw && sha(rr.uniforms) == sha(dev.uniforms) && sha(rr.electron) == sha(dev.electron) &&
                               sha(rr.neutral) == sha(dev.neutral) && rr.rounded == dev.rounded && rr.bernoulli == dev.bernoulli;
        const std::size_t mm_raw = mismatches(host.raw, dev.raw), mm_uni = mismatches(host.uniforms, dev.uniforms), mm_em = mismatches(host.electron, dev.electron),
                          mm_nm = mismatches(host.neutral, dev.neutral), mm_sr = mismatches(host.rounded, dev.rounded), mm_b = mismatches(host.bernoulli, dev.bernoulli);
        const std::size_t total_mm = mm_raw + mm_uni + mm_em + mm_nm + mm_sr + mm_b;
        std::size_t ones = 0;
        for (auto b : host.bernoulli) ones += b;
        const std::uint64_t values = 4 * n + 3 * n + 4 * n + 6 * n + n + n;
        js << "{\"records\":" << n << ",\"values_compared\":" << values << ",\"mismatches\":{\"raw_words\":" << mm_raw << ",\"conversions\":" << mm_uni << ",\"electron_marks\":" << mm_em
           << ",\"neutral_marks\":" << mm_nm << ",\"stochastic_round\":" << mm_sr << ",\"bernoulli\":" << mm_b << ",\"total\":" << total_mm << "}"
           << ",\"sha256\":{\"host\":{\"raw\":\"" << sha(host.raw) << "\",\"conversions\":\"" << sha(host.uniforms) << "\",\"electron\":\"" << sha(host.electron)
           << "\",\"neutral\":\"" << sha(host.neutral) << "\",\"rounded\":\"" << sha(host.rounded) << "\",\"bernoulli\":\"" << sha(host.bernoulli)
           << "\"},\"device\":{\"raw\":\"" << sha(dev.raw) << "\",\"conversions\":\"" << sha(dev.uniforms) << "\",\"electron\":\"" << sha(dev.electron)
           << "\",\"neutral\":\"" << sha(dev.neutral) << "\",\"rounded\":\"" << sha(dev.rounded) << "\",\"bernoulli\":\"" << sha(dev.bernoulli) << "\"}}"
           << ",\"bernoulli_ones\":" << ones << ",\"repeats_in_process\":" << repeats << ",\"repeat_outputs_identical\":" << (repeat_identical ? "true" : "false")
           << ",\"host_reference_ms\":" << hall_gpu::json_number(host_ms) << ",\"runs\":[" << timing.str() << "]"
           << ",\"header_sha256\":\"" << hall_gpu::sha256_hex(kHeader, std::strlen(kHeader)) << "\",\"kernel_sha256\":\""
           << hall_gpu::sha256_hex(kKernel, std::strlen(kKernel)) << "\",\"device\":" << hall_gpu::device_json(session.device())
           << ",\"compile\":" << hall_gpu::compile_json(module.record) << ",\"memory\":" << hall_gpu::memory_json(session.memory());
        session.release(d_in);
        session.release(d_out);
        session.release(h_in);
        session.release(h_out);
        teardown = session.close();
        js << ",\"teardown\":" << hall_gpu::teardown_json(teardown) << "}";
        if (total_mm) rc = 4;
        else if (!repeat_identical) rc = 2;
        if (teardown.quarantined || !teardown.released) rc = 3;
        std::cout << "gpu_parity records=" << n << " values=" << values << " mismatches=" << total_mm << " repeat_identical=" << repeat_identical << "\n";
    } catch (const std::exception& e) {
        std::cerr << "gpu_parity failed: " << e.what() << "\n";
        js.str("");
        js << "{\"error\":" << hall_gpu::json_string(e.what()) << "}";
        rc = 1;
    }
    std::ofstream(out / "gpu-parity.json") << js.str() << "\n";
    return rc;
}
