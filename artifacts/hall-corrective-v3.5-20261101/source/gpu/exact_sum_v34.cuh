// exact_sum_v34.cuh -- device exact sums of Hall v3.4 (NVRTC; include after forms_core.h of the native v7
// forms library, namespace pf7). The same superaccumulator definition as pf_v7::ExactSum (API-007.txt 1.1:
// every product exact, one rounding of the total, status bits); the kernel layout follows the forms
// library's sum_body (per-warp shared accumulators, lanes 0..4 add the five 32-bit chunks of one term,
// per-block partials, one finishing block). Any exact method gives the same bits, independent of the
// order and the partition of the terms. blockDim.x must be XS_BLOCK.
#pragma once

#define XS_BLOCK 256
#define XS_WARPS 8
#define XS_FULL 0xffffffffu

// Device result: the layout of pf_v7::gpu::SumOut.
struct XsOut {
    double value;                         // RN(exact sum)
    double scaled;                        // RN(scale x value)
    pf7_u32 status, parts;
    pf7_u64 terms;
    pf7_i64 digit[pf7::acc_digits];
};

template <class F>
__device__ void xs_body(pf7_u64 n, F term_of, pf7_i64* partial, pf7_u32* status) {
    using namespace pf7;
    __shared__ pf7_i64 acc[XS_WARPS][acc_digits];
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    for (int k = threadIdx.x; k < XS_WARPS * acc_digits; k += XS_BLOCK) (&acc[0][0])[k] = 0;
    __syncthreads();
    pf7_u32 flags = n >= (1ull << 31) ? sum_too_many_terms : 0u;
    const pf7_u64 stride = (pf7_u64)gridDim.x * XS_BLOCK;
    for (pf7_u64 base = (pf7_u64)blockIdx.x * XS_BLOCK + 32u * warp; base < n; base += stride) {
        const pf7_u64 i = base + lane;
        Term t;
        t.zero = 1; t.bad = 0; t.d = 0; t.neg = 0;
        for (int k = 0; k < 5; ++k) t.c[k] = 0;
        if (i < n) t = term_of(i);
        if (t.bad) flags |= sum_nonfinite_input;
        for (int j = 0; j < 32; ++j) {
            const int skip = __shfl_sync(XS_FULL, t.zero | t.bad, j);
            if (skip) continue;
            const int d = __shfl_sync(XS_FULL, t.d, j);
            const int neg = __shfl_sync(XS_FULL, t.neg, j);
            const pf7_u32 c0 = __shfl_sync(XS_FULL, t.c[0], j), c1 = __shfl_sync(XS_FULL, t.c[1], j), c2 = __shfl_sync(XS_FULL, t.c[2], j);
            const pf7_u32 c3 = __shfl_sync(XS_FULL, t.c[3], j), c4 = __shfl_sync(XS_FULL, t.c[4], j);
            if (lane < 5) {
                const pf7_u32 c = lane == 0 ? c0 : lane == 1 ? c1 : lane == 2 ? c2 : lane == 3 ? c3 : c4;
                acc[warp][d + lane] += neg ? -(pf7_i64)c : (pf7_i64)c;
            }
            __syncwarp();
        }
    }
    const int any = __syncthreads_or((int)flags);
    __shared__ pf7_u32 block_flags[XS_BLOCK];
    block_flags[threadIdx.x] = flags;
    __syncthreads();
    for (int k = threadIdx.x; k < acc_digits; k += XS_BLOCK) {
        pf7_i64 s = 0;
        for (int w = 0; w < XS_WARPS; ++w) s += acc[w][k];
        partial[(pf7_u64)blockIdx.x * acc_digits + k] = s;
    }
    if (threadIdx.x == 0) {
        pf7_u32 f = 0;
        if (any) for (int k = 0; k < XS_BLOCK; ++k) f |= block_flags[k];
        status[blockIdx.x] = f;
    }
}

// One block of XS_BLOCK threads: merge the partials in block order, normalise, round once.
__device__ void xs_finish(pf7_u32 parts, const pf7_i64* partial, const pf7_u32* status, pf7_u64 terms, double scale,
                          XsOut* out) {
    using namespace pf7;
    __shared__ pf7_i64 d[acc_digits];
    for (int k = threadIdx.x; k < acc_digits; k += XS_BLOCK) {
        pf7_i64 s = 0;
        for (pf7_u32 p = 0; p < parts; ++p) s += partial[(pf7_u64)p * acc_digits + k];
        d[k] = s;
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        pf7_u32 st = 0;
        for (pf7_u32 p = 0; p < parts; ++p) st |= status[p];
        normalise(d);
        for (int k = 0; k < acc_digits; ++k) out->digit[k] = d[k];
        const double v = round_sum(d, &st);
        out->value = v;
        out->scaled = scale * v;
        if (((bits_of(out->scaled) >> 52) & 0x7ff) == 0x7ff) st |= sum_overflow;
        out->status = st;
        out->parts = parts;
        out->terms = terms;
    }
}
