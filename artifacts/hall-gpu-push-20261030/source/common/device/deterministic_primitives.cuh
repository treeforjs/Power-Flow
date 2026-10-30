// deterministic_primitives.cuh -- NVRTC device code shared by GPU programs.
// Every result is independent of thread scheduling:
//  * exclusive scans of unsigned 32-bit counts (integer, exact);
//  * stable LSD radix-sort passes on 32-bit keys with a 32-bit payload
//    (8-bit digits; per-tile histograms; order-preserving scatter that ranks
//    equal digits by warp match/popc in index order, never by atomics);
//  * segment bounds of a sorted key array;
//  * max of non-negative finite doubles and min of u32 (exact operations).
// Histogram counting uses shared-memory integer atomics; counts do not depend
// on their order. No floating-point value is ever accumulated with atomics.
// All kernels assume blockDim.x == HGD_BLOCK (256).
#pragma once

#define HGD_BLOCK 256
#define HGD_WARPS (HGD_BLOCK / 32)
#define HGD_ITEMS 8
#define HGD_TILE (HGD_BLOCK * HGD_ITEMS)

typedef unsigned int hgd_u32;
typedef unsigned long long hgd_u64;

// Exclusive scan across the block of one value per thread; *total receives the block sum.
__device__ inline hgd_u32 hgd_block_exclusive_scan(hgd_u32 v, hgd_u32* total) {
    __shared__ hgd_u32 warp_sums[HGD_WARPS];
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    hgd_u32 x = v;
    for (int o = 1; o < 32; o <<= 1) {
        const hgd_u32 y = __shfl_up_sync(0xffffffffu, x, o);
        if (lane >= o) x += y;
    }
    if (lane == 31) warp_sums[warp] = x;
    __syncthreads();
    if (warp == 0) {
        hgd_u32 w = lane < HGD_WARPS ? warp_sums[lane] : 0u;
        for (int o = 1; o < 32; o <<= 1) {
            const hgd_u32 y = __shfl_up_sync(0xffffffffu, w, o);
            if (lane >= o) w += y;
        }
        if (lane < HGD_WARPS) warp_sums[lane] = w;
    }
    __syncthreads();
    const hgd_u32 prefix = warp ? warp_sums[warp - 1] : 0u;
    *total = warp_sums[HGD_WARPS - 1];
    __syncthreads();
    return prefix + x - v;
}

// Scan step 1: per-tile sums. grid = ceil(n / HGD_TILE).
extern "C" __global__ void hgd_scan_reduce(const hgd_u32* in, hgd_u32 n, hgd_u32* partial) {
    const hgd_u32 base = blockIdx.x * HGD_TILE + threadIdx.x * HGD_ITEMS;
    hgd_u32 local = 0;
    for (int k = 0; k < HGD_ITEMS; ++k) if (base + k < n) local += in[base + k];
    hgd_u32 total = 0;
    (void)hgd_block_exclusive_scan(local, &total);
    if (threadIdx.x == 0) partial[blockIdx.x] = total;
}

// Scan step 2: exclusive scan of the tile sums in place, by one block. total_out[0] = grand total.
extern "C" __global__ void hgd_scan_partials(hgd_u32* partial, hgd_u32 tiles, hgd_u32* total_out) {
    hgd_u32 carry = 0;
    for (hgd_u32 start = 0; start < tiles; start += HGD_BLOCK) {
        const hgd_u32 i = start + threadIdx.x;
        const hgd_u32 v = i < tiles ? partial[i] : 0u;
        hgd_u32 chunk = 0;
        const hgd_u32 excl = hgd_block_exclusive_scan(v, &chunk);
        if (i < tiles) partial[i] = carry + excl;
        carry += chunk;
    }
    if (threadIdx.x == 0) total_out[0] = carry;
}

// Scan step 3: exclusive scan within each tile plus the tile offset. in may alias out.
extern "C" __global__ void hgd_scan_apply(const hgd_u32* in, hgd_u32 n, const hgd_u32* partial, hgd_u32* out) {
    const hgd_u32 base = blockIdx.x * HGD_TILE + threadIdx.x * HGD_ITEMS;
    hgd_u32 item[HGD_ITEMS];
    hgd_u32 local = 0;
    for (int k = 0; k < HGD_ITEMS; ++k) { item[k] = base + k < n ? in[base + k] : 0u; local += item[k]; }
    hgd_u32 total = 0;
    hgd_u32 running = partial[blockIdx.x] + hgd_block_exclusive_scan(local, &total);
    for (int k = 0; k < HGD_ITEMS; ++k) {
        if (base + k < n) out[base + k] = running;
        running += item[k];
    }
}

// Radix pass step 1: digit histogram of each tile, stored digit-major: hist[d * tiles + tile].
extern "C" __global__ void hgd_radix_histogram(const hgd_u32* keys, hgd_u32 n, hgd_u32 shift, hgd_u32 tiles, hgd_u32* hist) {
    __shared__ hgd_u32 counts[256];
    counts[threadIdx.x] = 0u;
    __syncthreads();
    for (int r = 0; r < HGD_ITEMS; ++r) {
        const hgd_u32 i = blockIdx.x * HGD_TILE + r * HGD_BLOCK + threadIdx.x;
        if (i < n) atomicAdd(&counts[(keys[i] >> shift) & 255u], 1u);
    }
    __syncthreads();
    hist[threadIdx.x * tiles + blockIdx.x] = counts[threadIdx.x];
}

// Radix pass step 2 (after an exclusive scan of hist): stable scatter. Items of a
// tile are processed in rounds of HGD_BLOCK consecutive indices; equal digits are
// ranked by warp, then lane, i.e. by index. Requires sm_70+ (__match_any_sync).
extern "C" __global__ void hgd_radix_scatter(const hgd_u32* key_in, const hgd_u32* pay_in, hgd_u32 n, hgd_u32 shift,
                                             const hgd_u32* offsets, hgd_u32 tiles, hgd_u32* key_out, hgd_u32* pay_out) {
    __shared__ hgd_u32 base[256];
    __shared__ hgd_u32 warp_count[HGD_WARPS][256];
    const hgd_u32 tid = threadIdx.x, lane = tid & 31u, warp = tid >> 5;
    base[tid] = offsets[tid * tiles + blockIdx.x];
    for (int r = 0; r < HGD_ITEMS; ++r) {
        for (int w = 0; w < HGD_WARPS; ++w) warp_count[w][tid] = 0u;
        __syncthreads();
        const hgd_u32 i = blockIdx.x * HGD_TILE + r * HGD_BLOCK + tid;
        const bool valid = i < n;
        hgd_u32 key = 0u, pay = 0u, digit = 256u;
        if (valid) { key = key_in[i]; pay = pay_in[i]; digit = (key >> shift) & 255u; }
        const hgd_u32 peers = __match_any_sync(0xffffffffu, digit);
        const hgd_u32 rank = __popc(peers & ((1u << lane) - 1u));
        if (valid && hgd_u32(__ffs(int(peers)) - 1) == lane) warp_count[warp][digit] = __popc(peers);
        __syncthreads();
        {
            hgd_u32 run = base[tid];
            for (int w = 0; w < HGD_WARPS; ++w) { const hgd_u32 c = warp_count[w][tid]; warp_count[w][tid] = run; run += c; }
            base[tid] = run;
        }
        __syncthreads();
        if (valid) {
            const hgd_u32 pos = warp_count[warp][digit] + rank;
            key_out[pos] = key;
            pay_out[pos] = pay;
        }
        __syncthreads();
    }
}

extern "C" __global__ void hgd_iota(hgd_u32* out, hgd_u32 n) {
    const hgd_u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = i;
}

// begin/end must be zeroed (key range) before the launch; empty keys keep begin == end == 0.
extern "C" __global__ void hgd_segment_bounds(const hgd_u32* sorted_keys, hgd_u32 n, hgd_u32* begin, hgd_u32* end) {
    const hgd_u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const hgd_u32 k = sorted_keys[i];
    if (i == 0 || sorted_keys[i - 1] != k) begin[k] = i;
    if (i == n - 1 || sorted_keys[i + 1] != k) end[k] = i + 1;
}

// Max of non-negative finite doubles: out[block] = max over the block's tile.
extern "C" __global__ void hgd_reduce_max_f64(const double* v, hgd_u32 n, double* out) {
    __shared__ double s[HGD_BLOCK];
    double m = 0.0;
    for (int r = 0; r < HGD_ITEMS; ++r) {
        const hgd_u32 i = blockIdx.x * HGD_TILE + r * HGD_BLOCK + threadIdx.x;
        if (i < n) { const double x = v[i]; m = (m < x) ? x : m; }
    }
    s[threadIdx.x] = m;
    __syncthreads();
    for (int o = HGD_BLOCK / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o) { const double x = s[threadIdx.x + o]; s[threadIdx.x] = (s[threadIdx.x] < x) ? x : s[threadIdx.x]; }
        __syncthreads();
    }
    if (threadIdx.x == 0) out[blockIdx.x] = s[0];
}

extern "C" __global__ void hgd_reduce_min_u32(const hgd_u32* v, hgd_u32 n, hgd_u32* out) {
    __shared__ hgd_u32 s[HGD_BLOCK];
    hgd_u32 m = 0xffffffffu;
    for (int r = 0; r < HGD_ITEMS; ++r) {
        const hgd_u32 i = blockIdx.x * HGD_TILE + r * HGD_BLOCK + threadIdx.x;
        if (i < n) { const hgd_u32 x = v[i]; m = x < m ? x : m; }
    }
    s[threadIdx.x] = m;
    __syncthreads();
    for (int o = HGD_BLOCK / 2; o > 0; o >>= 1) {
        if (threadIdx.x < o) { const hgd_u32 x = s[threadIdx.x + o]; s[threadIdx.x] = x < s[threadIdx.x] ? x : s[threadIdx.x]; }
        __syncthreads();
    }
    if (threadIdx.x == 0) out[blockIdx.x] = s[0];
}
