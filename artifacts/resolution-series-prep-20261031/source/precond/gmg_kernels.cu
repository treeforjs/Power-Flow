// Deterministic NVRTC kernels of the GMG-preconditioned CG (resolution-series prep). One thread per row
// (grid-stride over NBLK x TPB threads for the reductions, exact launch for the row kernels); no atomics; dot
// products by a fixed per-thread order followed by fixed shared-memory trees (same tree as the Hall gpu1 CG).
#include "gmg_core.h"
#define TPB 256
#define NBLK 512
#define ROWS(i, n) for (int i = blockIdx.x * TPB + threadIdx.x; i < (n); i += NBLK * TPB)
#define ROW(i, n) const int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= (n)) return

__device__ __forceinline__ double block_tree(double* sh, double v) {
    const int t = threadIdx.x;
    sh[t] = v;
    __syncthreads();
    for (int s = TPB / 2; s > 0; s >>= 1) {
        if (t < s) sh[t] = sh[t] + sh[t + s];
        __syncthreads();
    }
    return sh[0];
}
// r = b - K x
extern "C" __global__ void g_residual(int n, const int* __restrict__ rp, const int* __restrict__ ci, const double* __restrict__ v,
                                      const double* __restrict__ x, const double* __restrict__ b, double* __restrict__ r) {
    ROW(i, n);
    r[i] = gmg_row_residual(i, rp, ci, v, x, b);
}
// y = K x
extern "C" __global__ void g_mul(int n, const int* __restrict__ rp, const int* __restrict__ ci, const double* __restrict__ v,
                                 const double* __restrict__ x, double* __restrict__ y) {
    ROW(i, n);
    y[i] = gmg_row_mul(i, rp, ci, v, x);
}
// y = y + K x
extern "C" __global__ void g_mul_add(int n, const int* __restrict__ rp, const int* __restrict__ ci, const double* __restrict__ v,
                                     const double* __restrict__ x, double* __restrict__ y) {
    ROW(i, n);
    y[i] = y[i] + gmg_row_mul(i, rp, ci, v, x);
}
// x = x + d * r
extern "C" __global__ void g_jacobi(int n, const double* __restrict__ d, const double* __restrict__ r, double* __restrict__ x) {
    ROW(i, n);
    x[i] = x[i] + d[i] * r[i];
}
// y = d * x
extern "C" __global__ void g_scale(int n, const double* __restrict__ d, const double* __restrict__ x, double* __restrict__ y) {
    ROW(i, n);
    y[i] = d[i] * x[i];
}
// z = fixed ? r : s * x   (level-0 output of the V-cycle, A = M^-1/2 K M^-1/2)
extern "C" __global__ void g_finish(int n, const unsigned char* __restrict__ fixed, const double* __restrict__ s, const double* __restrict__ r,
                                    const double* __restrict__ x, double* __restrict__ z) {
    ROW(i, n);
    z[i] = fixed[i] ? r[i] : s[i] * x[i];
}
// b0 = fixed ? 0 : s * r
extern "C" __global__ void g_start(int n, const unsigned char* __restrict__ fixed, const double* __restrict__ s, const double* __restrict__ r,
                                   double* __restrict__ b) {
    ROW(i, n);
    b[i] = fixed[i] ? 0.0 : s[i] * r[i];
}
extern "C" __global__ void g_zero(int n, double* __restrict__ x) {
    ROW(i, n);
    x[i] = 0.0;
}
extern "C" __global__ void g_copy(int n, const double* __restrict__ a, double* __restrict__ b) {
    ROW(i, n);
    b[i] = a[i];
}
// x = x + alpha p ; r = r - alpha q   (alpha read from the device scalar slot)
extern "C" __global__ void g_cg_update(int n, const double* __restrict__ sc, int slot, const double* __restrict__ p, const double* __restrict__ q,
                                       double* __restrict__ x, double* __restrict__ r) {
    ROW(i, n);
    const double alpha = sc[slot];
    x[i] = x[i] + alpha * p[i];
    r[i] = r[i] - alpha * q[i];
}
// p = z + beta p
extern "C" __global__ void g_cg_direction(int n, const double* __restrict__ sc, int slot, const double* __restrict__ z, double* __restrict__ p) {
    ROW(i, n);
    p[i] = z[i] + sc[slot] * p[i];
}
// partial sums of a . b (fixed order per thread, fixed block tree)
extern "C" __global__ void g_dot(int n, const double* __restrict__ a, const double* __restrict__ b, double* part) {
    __shared__ double sh[TPB];
    double s = 0.0;
    ROWS(i, n) s = s + a[i] * b[i];
    const double v = block_tree(sh, s);
    if (threadIdx.x == 0) part[blockIdx.x] = v;
}
// final tree over the NBLK partials into sc[slot]
extern "C" __global__ void g_fin(const double* part, double* sc, int slot) {
    __shared__ double sh[NBLK];
    const int t = threadIdx.x;
    sh[t] = part[t];
    __syncthreads();
    for (int s = NBLK / 2; s > 0; s >>= 1) {
        if (t < s) sh[t] = sh[t] + sh[t + s];
        __syncthreads();
    }
    if (t == 0) sc[slot] = sh[0];
}
// scalar step of CG on the device (no host round trip): sc[dst] = sc[a] / sc[b]
extern "C" __global__ void g_div(double* sc, int dst, int a, int b) {
    if (threadIdx.x == 0 && blockIdx.x == 0) sc[dst] = sc[a] / sc[b];
}
// r = b - q
extern "C" __global__ void g_sub(int n, const double* __restrict__ b, const double* __restrict__ q, double* __restrict__ r) {
    ROW(i, n);
    r[i] = b[i] - q[i];
}
// sc[dst] = sc[src]
extern "C" __global__ void g_mov(double* sc, int dst, int src) {
    if (threadIdx.x == 0 && blockIdx.x == 0) sc[dst] = sc[src];
}
