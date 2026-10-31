// push_kernels_v3.cu -- NVRTC program of the Hall v3 motion path: the v1 kernels
// (hall-gpu-resident-v1-20261030/push/push_kernels.cu) with the v3 carrier body
// (first-hit truncation) and the v3 endpoint node reduction (sink records excluded
// from the free-charge census). No atomics on floating-point data.
#include "push_bodies.hpp"
#include "push_bodies_v3.hpp"
#include "deterministic_primitives.cuh"

using namespace hgp;

extern "C" __global__ void hgp_carrier_v3(const CarrierArgsV3 a) { carrier_body_v3(a, blockIdx.x * blockDim.x + threadIdx.x); }
extern "C" __global__ void hgp_emit_edges(const EmitArgs a) { emit_body(a, blockIdx.x * blockDim.x + threadIdx.x); }
extern "C" __global__ void hgp_reduce_edges(const NodeReduceArgs a) { reduce_edges_body(a, blockIdx.x * blockDim.x + threadIdx.x); }
extern "C" __global__ void hgp_reduce_start_nodes(const NodeReduceArgs a) { reduce_start_body(a, blockIdx.x * blockDim.x + threadIdx.x); }
extern "C" __global__ void hgp_reduce_end_nodes_v3(const EndNodeArgsV3 a) { reduce_end_body_v3(a, blockIdx.x * blockDim.x + threadIdx.x); }

extern "C" __global__ void hgp_plasma(const double* omega, const double* VN, hgp_u32 nn, double dt, double* out) {
    const hgp_u32 k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nn) return;
    out[k] = dt * HGP_SQRT(omega[k] / VN[k]);
}

extern "C" __global__ void hgp_finalize(const FinalizeArgs a) {
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    const hgp_u32 r1 = a.min_fail[0], r2 = a.min_fail[1], r3 = a.min_fail[2];
    a.status[0] = r1; a.status[1] = r1 != 0xffffffffu ? hgp_u32(a.code1[r1]) : 0u;
    a.status[2] = r2; a.status[3] = r2 != 0xffffffffu ? hgp_u32(a.code2[r2]) : 0u;
    a.status[4] = r3; a.status[5] = r3 != 0xffffffffu ? hgp_u32(a.code3[r3]) : 0u;
    a.status[6] = a.edge_total[0];
    a.status[7] = 0;
}
