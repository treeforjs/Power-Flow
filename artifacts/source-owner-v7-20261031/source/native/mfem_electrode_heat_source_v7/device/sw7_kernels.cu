// sw7_kernels.cu -- NVRTC entry points of the v7 wall kernel (compiled at run time by hall_gpu::Session with
// --fmad=false --prec-div=true --prec-sqrt=true --ftz=false, architecture of the device selected at run time).
// One thread evaluates one wall event with the same function as the CPU reference (sw7::wall_event), so the
// outputs are bit-identical. Ledger sums: fixed shared-memory tree per block, then one sequential pass over the
// block partials; integers only (64-bit fields split into 32-bit limbs), no atomics, independent of scheduling.
#include "sw7_wall.h"

#define SW7_NF 11
#define SW7_RB 256

extern "C" __global__ void sw7_sizes(unsigned long long* out) {
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    SW7Out o;
    SW7Params p;
    out[0] = sizeof(SW7Params);
    out[1] = sizeof(SW7Hit);
    out[2] = sizeof(SW7Out);
    out[3] = sizeof(SW7Species);
    out[4] = (unsigned long long)((char*)&o.k_in_q - (char*)&o);
    out[5] = (unsigned long long)((char*)&o.sigma_q - (char*)&o);
    out[6] = (unsigned long long)((char*)&p.pn[0] - (char*)&p);
    out[7] = (unsigned long long)((char*)&p.enable_isd - (char*)&p);
    out[8] = (unsigned long long)((char*)&o.u[0][0] - (char*)&o);
}

extern "C" __global__ void sw7_wall_events(const SW7Params* params, const SW7Species* table, const SW7Hit* hits, unsigned n,
                                           unsigned long long seed, unsigned long long step, SW7Out* out) {
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    sw7::wall_event(params[hits[i].patch], table, hits[i], seed, step, &out[i]);
}

__device__ inline void sw7_fields(const SW7Out& o, long long* f) {
    f[0] = o.e_in_q; f[1] = o.e_out_q; f[2] = o.wall_q; f[3] = o.k_in_q; f[4] = o.ads_q; f[5] = o.sigma_q;
    f[6] = o.n_out; f[7] = o.electrons_out; f[8] = o.stuck; f[9] = o.charge_in; f[10] = o.charge_out;
}

// blockDim.x must be SW7_RB. partial[block * 2 NF + 2k] = sum of high limbs (x >> 32), [.. + 1] = sum of low limbs.
extern "C" __global__ void sw7_reduce_blocks(const SW7Out* out, unsigned n, long long* partial) {
    __shared__ long long s[SW7_RB][2 * SW7_NF];
    const unsigned t = threadIdx.x, i = blockIdx.x * SW7_RB + t;
    long long f[SW7_NF];
    for (int k = 0; k < SW7_NF; ++k) f[k] = 0;
    if (i < n) sw7_fields(out[i], f);
    for (int k = 0; k < SW7_NF; ++k) {
        s[t][2 * k] = f[k] >> 32;
        s[t][2 * k + 1] = (long long)((unsigned long long)f[k] & 0xFFFFFFFFull);
    }
    __syncthreads();
    for (unsigned stride = SW7_RB / 2; stride > 0; stride >>= 1) {
        if (t < stride)
            for (int j = 0; j < 2 * SW7_NF; ++j) s[t][j] += s[t + stride][j];
        __syncthreads();
    }
    if (t == 0)
        for (int j = 0; j < 2 * SW7_NF; ++j) partial[(unsigned long long)blockIdx.x * 2 * SW7_NF + j] = s[0][j];
}

extern "C" __global__ void sw7_reduce_final(const long long* partial, unsigned blocks, long long* total) {
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    for (int j = 0; j < 2 * SW7_NF; ++j) {
        long long acc = 0;
        for (unsigned b = 0; b < blocks; ++b) acc += partial[(unsigned long long)b * 2 * SW7_NF + j];
        total[j] = acc;
    }
}
