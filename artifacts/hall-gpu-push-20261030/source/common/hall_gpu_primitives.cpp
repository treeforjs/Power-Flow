#include "hall_gpu_primitives.hpp"

#include <climits>

namespace hall_gpu {
namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kPrimitivesSource[] = {
#embed "device/deterministic_primitives.cuh"
    , 0};
#pragma GCC diagnostic pop

unsigned u32_count(std::size_t n) {
    if (n > std::size_t(UINT_MAX) / 2) throw Error("deterministic primitive size exceeds the 32-bit index range");
    return unsigned(n);
}
}  // namespace

const ProgramSource& Primitives::device_source() {
    static const ProgramSource source{"deterministic_primitives.cuh", kPrimitivesSource};
    return source;
}

Primitives::Primitives(Session& s, const Module& m)
    : s_(s),
      scan_reduce_(s.kernel(m, "hgd_scan_reduce")),
      scan_partials_(s.kernel(m, "hgd_scan_partials")),
      scan_apply_(s.kernel(m, "hgd_scan_apply")),
      histogram_(s.kernel(m, "hgd_radix_histogram")),
      scatter_(s.kernel(m, "hgd_radix_scatter")),
      iota_(s.kernel(m, "hgd_iota")),
      bounds_(s.kernel(m, "hgd_segment_bounds")),
      max_(s.kernel(m, "hgd_reduce_max_f64")),
      min_(s.kernel(m, "hgd_reduce_min_u32")) {}

void Primitives::exclusive_scan(std::uint64_t in, std::size_t n, std::uint64_t out, std::uint64_t scratch, std::uint64_t total) {
    unsigned count = u32_count(n);
    unsigned t = unsigned(tiles(n));
    if (count == 0) {
        // total = 0 through a one-thread partial scan over zero tiles.
        void* a[] = {&scratch, &t, &total};
        s_.launch(scan_partials_, {1}, {block}, 0, a);
        return;
    }
    void* a1[] = {&in, &count, &scratch};
    s_.launch(scan_reduce_, {t}, {block}, 0, a1);
    void* a2[] = {&scratch, &t, &total};
    s_.launch(scan_partials_, {1}, {block}, 0, a2);
    void* a3[] = {&in, &count, &scratch, &out};
    s_.launch(scan_apply_, {t}, {block}, 0, a3);
}

bool Primitives::radix_sort(std::uint64_t key_a, std::uint64_t pay_a, std::uint64_t key_b, std::uint64_t pay_b,
                            std::size_t n, unsigned key_bits, std::uint64_t scratch) {
    unsigned count = u32_count(n);
    if (count == 0) return false;
    unsigned t = unsigned(tiles(n));
    const std::size_t hist_elems = std::size_t(256) * t;
    std::uint64_t hist = scratch;
    std::uint64_t hist_scan = scratch + hist_elems * 4;
    std::uint64_t hist_total = hist_scan + tiles(hist_elems) * 4;
    bool in_b = false;
    for (unsigned shift = 0; shift < key_bits; shift += 8) {
        std::uint64_t ki = in_b ? key_b : key_a, pi = in_b ? pay_b : pay_a;
        std::uint64_t ko = in_b ? key_a : key_b, po = in_b ? pay_a : pay_b;
        void* h[] = {&ki, &count, &shift, &t, &hist};
        s_.launch(histogram_, {t}, {block}, 0, h);
        exclusive_scan(hist, hist_elems, hist, hist_scan, hist_total);
        void* sc[] = {&ki, &pi, &count, &shift, &hist, &t, &ko, &po};
        s_.launch(scatter_, {t}, {block}, 0, sc);
        in_b = !in_b;
    }
    return in_b;
}

void Primitives::iota(std::uint64_t out, std::size_t n) {
    unsigned count = u32_count(n);
    void* a[] = {&out, &count};
    s_.launch(iota_, {unsigned((n + block - 1) / block)}, {block}, 0, a);
}

void Primitives::segment_bounds(std::uint64_t sorted_keys, std::size_t n, std::uint64_t begin, std::uint64_t end) {
    unsigned count = u32_count(n);
    void* a[] = {&sorted_keys, &count, &begin, &end};
    s_.launch(bounds_, {unsigned((n + block - 1) / block)}, {block}, 0, a);
}

void Primitives::reduce_max_f64(std::uint64_t values, std::size_t n, std::uint64_t scratch, std::uint64_t out) {
    std::uint64_t src = values;
    std::uint64_t ping = scratch, pong = scratch + (tiles(n) + 1) * 8;
    unsigned count = u32_count(n);
    do {
        unsigned t = unsigned(tiles(count ? count : 1));
        std::uint64_t dst = (t == 1) ? out : ping;
        void* a[] = {&src, &count, &dst};
        s_.launch(max_, {t}, {block}, 0, a);
        src = dst;
        count = t;
        std::swap(ping, pong);
    } while (count > 1);
}

void Primitives::reduce_min_u32(std::uint64_t values, std::size_t n, std::uint64_t scratch, std::uint64_t out) {
    std::uint64_t src = values;
    std::uint64_t ping = scratch, pong = scratch + (tiles(n) + 1) * 8;
    unsigned count = u32_count(n);
    do {
        unsigned t = unsigned(tiles(count ? count : 1));
        std::uint64_t dst = (t == 1) ? out : ping;
        void* a[] = {&src, &count, &dst};
        s_.launch(min_, {t}, {block}, 0, a);
        src = dst;
        count = t;
        std::swap(ping, pong);
    } while (count > 1);
}

}  // namespace hall_gpu
