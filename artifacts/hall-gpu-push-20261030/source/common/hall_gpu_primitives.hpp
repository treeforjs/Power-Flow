// hall_gpu_primitives.hpp -- host launchers for common/device/deterministic_primitives.cuh.
// Add Primitives::device_source() to an NVRTC program (as header
// "deterministic_primitives.cuh"), compile it with Session::compile_and_load,
// then construct Primitives on the resulting Module. All work is enqueued on
// the session stream; addresses are raw device addresses (DeviceBuffer::address
// plus a byte offset) so callers can address sub-ranges.
#pragma once
#include "hall_gpu_common.hpp"

namespace hall_gpu {

class Primitives {
public:
    static constexpr unsigned block = 256, items = 8, tile = block * items;
    static const ProgramSource& device_source();   // embedded deterministic_primitives.cuh

    Primitives(Session&, const Module&);

    static std::size_t tiles(std::size_t n) { return (n + tile - 1) / tile; }
    // Scratch (u32 elements) needed by exclusive_scan for n inputs.
    static std::size_t scan_scratch(std::size_t n) { return tiles(n) + 1; }
    // Scratch (u32 elements) needed by radix_sort for n records.
    static std::size_t sort_scratch(std::size_t n) { return 256 * tiles(n) + scan_scratch(256 * tiles(n)); }
    // Scratch (8-byte elements) needed by reduce_max_f64 / reduce_min_u32.
    static std::size_t reduce_scratch(std::size_t n) { return 2 * (tiles(n) + 1); }

    // out[i] = sum in[0..i) (u32, exact). Grand total written to total (one u32).
    void exclusive_scan(std::uint64_t in, std::size_t n, std::uint64_t out, std::uint64_t scratch, std::uint64_t total);
    // Stable LSD sort of (key, payload) by the low key_bits bits of key. Buffers a and b
    // are ping-pong storage of n entries each; returns true when the result is in b.
    bool radix_sort(std::uint64_t key_a, std::uint64_t pay_a, std::uint64_t key_b, std::uint64_t pay_b,
                    std::size_t n, unsigned key_bits, std::uint64_t scratch);
    void iota(std::uint64_t out, std::size_t n);
    // begin/end (u32 per key) must be zeroed beforehand.
    void segment_bounds(std::uint64_t sorted_keys, std::size_t n, std::uint64_t begin, std::uint64_t end);
    // Exact max of non-negative finite doubles / min of u32; result written to out (one element).
    void reduce_max_f64(std::uint64_t values, std::size_t n, std::uint64_t scratch, std::uint64_t out);
    void reduce_min_u32(std::uint64_t values, std::size_t n, std::uint64_t scratch, std::uint64_t out);

private:
    Session& s_;
    Kernel scan_reduce_, scan_partials_, scan_apply_, histogram_, scatter_, iota_, bounds_, max_, min_;
};

}  // namespace hall_gpu
