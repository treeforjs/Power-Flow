// Device code for the deterministic FP64 Jacobi CG (compiled at run time by
// NVRTC with --fmad=false --prec-div=true --prec-sqrt=true --ftz=false).
// Determinism: row-per-thread products, fixed grid (NBLK x TPB) grid-stride
// loops for every reduction, shared-memory halving trees, one finishing block
// of NBLK threads, no atomics. Iteration kernels read the control block and
// return immediately once the stopping rule has fired, so launching a fixed
// batch of iterations never changes the result.
#pragma once
#include <cstdint>

namespace hgf {
// Host mirror of the device control block (same declaration order).
struct Ctl {
    double rz, rr, rnorm, alpha, beta, pap, target, bnorm, tnorm, spare;
    std::int32_t it, done, error, maxit, corr, pad0, pad1, pad2;
};
static_assert(sizeof(Ctl) == 112, "control block layout");

inline const char* kernel_source = R"CUDA(
#define TPB 256
#define NBLK 512
struct Ctl {
    double rz, rr, rnorm, alpha, beta, pap, target, bnorm, tnorm, spare;
    int it, done, error, maxit, corr, pad0, pad1, pad2;
};
__device__ __forceinline__ bool stopped(const Ctl* c) { return c->done != 0 || c->error != 0; }
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
__device__ __forceinline__ double final_tree(double* sh, const double* part) {
    const int t = threadIdx.x;
    sh[t] = part[t];
    __syncthreads();
    for (int s = NBLK / 2; s > 0; s >>= 1) {
        if (t < s) sh[t] = sh[t] + sh[t + s];
        __syncthreads();
    }
    return sh[0];
}
#define ROWS(i, n) for (int i = blockIdx.x * TPB + threadIdx.x; i < (n); i += NBLK * TPB)

extern "C" __global__ void k_dot_self(int n, const double* __restrict__ b, double* part) {
    __shared__ double sh[TPB];
    double s = 0.0;
    ROWS(i, n) s += b[i] * b[i];
    const double v = block_tree(sh, s);
    if (threadIdx.x == 0) part[blockIdx.x] = v;
}
extern "C" __global__ void k_fin_target(const double* part, Ctl* c, int maxit) {
    __shared__ double sh[NBLK];
    const double v = final_tree(sh, part);
    if (threadIdx.x == 0) {
        c->bnorm = sqrt(v);
        c->target = fmax(1e-30, 2e-13 * c->bnorm);
        c->done = 0; c->error = 0; c->it = 0; c->maxit = maxit; c->corr = 0;
    }
}
extern "C" __global__ void k_init(int n, const double* __restrict__ r, const double* __restrict__ diag,
                                  double* __restrict__ z, double* __restrict__ p, double* part_rz, double* part_rr) {
    __shared__ double sh[TPB];
    double s1 = 0.0, s2 = 0.0;
    ROWS(i, n) {
        const double ri = r[i];
        const double zi = ri / diag[i];
        z[i] = zi; p[i] = zi;
        s1 += ri * zi; s2 += ri * ri;
    }
    double v = block_tree(sh, s1);
    if (threadIdx.x == 0) part_rz[blockIdx.x] = v;
    __syncthreads();
    v = block_tree(sh, s2);
    if (threadIdx.x == 0) part_rr[blockIdx.x] = v;
}
extern "C" __global__ void k_fin_init(const double* part_rz, const double* part_rr, Ctl* c, int corr, double* hist_rnorm) {
    __shared__ double sh[NBLK];
    const double rz = final_tree(sh, part_rz);
    __syncthreads();
    const double rr = final_tree(sh, part_rr);
    if (threadIdx.x == 0) {
        c->rz = rz; c->rr = rr; c->rnorm = sqrt(rr); c->it = 0; c->corr = corr; c->error = 0;
        c->done = (0 < c->maxit && c->rnorm > c->target * .2) ? 0 : 1;
        hist_rnorm[corr * (c->maxit + 1)] = c->rnorm;
    }
}
// Assembled CSR: Ap = A p, optional partial of p.Ap.
extern "C" __global__ void k_spmv_csr(int n, const int* __restrict__ rp, const int* __restrict__ ci, const double* __restrict__ val,
                                      const double* __restrict__ p, double* __restrict__ Ap, double* part, const Ctl* c, int gated, int with_dot) {
    if (gated && stopped(c)) return;
    __shared__ double sh[TPB];
    double s = 0.0;
    ROWS(i, n) {
        double y = 0.0;
        const int end = rp[i + 1];
        for (int j = rp[i]; j < end; ++j) y += val[j] * p[ci[j]];
        Ap[i] = y;
        s += p[i] * y;
    }
    if (with_dot) {
        const double v = block_tree(sh, s);
        if (threadIdx.x == 0) part[blockIdx.x] = v;
    }
}
// Factored apply, part 1: Bh = H .* (C E), E = v / sqrt(M) off PEC edges.
extern "C" __global__ void k_face(int nb, const int* __restrict__ fe, const double* __restrict__ fc, const double* __restrict__ H,
                                  const double* __restrict__ sqrtM, const unsigned char* __restrict__ fixed,
                                  const double* __restrict__ v, double* __restrict__ Bh, const Ctl* c, int gated) {
    if (gated && stopped(c)) return;
    const int b = blockIdx.x * TPB + threadIdx.x;
    if (b >= nb) return;
    double s = 0.0;
    for (int k = 0; k < 4; ++k) {
        const int e = fe[4 * b + k];
        const double E = fixed[e] ? 0.0 : v[e] / sqrtM[e];
        s += fc[4 * b + k] * E;
    }
    Bh[b] = s * H[b];
}
// Factored apply, part 2: Ap = v + (a E + q C^T Bh) / sqrt(M) off PEC edges.
extern "C" __global__ void k_edge(int n, const int* __restrict__ ctp, const int* __restrict__ ctf, const double* __restrict__ ctc,
                                  const double* __restrict__ Bh, const double* __restrict__ sqrtM, const double* __restrict__ a, double q,
                                  const unsigned char* __restrict__ fixed, const double* __restrict__ v, double* __restrict__ Ap,
                                  double* part, const Ctl* c, int gated, int with_dot) {
    if (gated && stopped(c)) return;
    __shared__ double sh[TPB];
    double s = 0.0;
    ROWS(i, n) {
        double y;
        if (fixed[i]) y = v[i];
        else {
            const double E = v[i] / sqrtM[i];
            double K = 0.0;
            const int end = ctp[i + 1];
            for (int j = ctp[i]; j < end; ++j) K += ctc[j] * Bh[ctf[j]];
            y = v[i] + (a[i] * E + q * K) / sqrtM[i];
        }
        Ap[i] = y;
        s += v[i] * y;
    }
    if (with_dot) {
        const double r = block_tree(sh, s);
        if (threadIdx.x == 0) part[blockIdx.x] = r;
    }
}
extern "C" __global__ void k_fin_alpha(const double* part, Ctl* c, double* hist_alpha) {
    if (stopped(c)) return;
    __shared__ double sh[NBLK];
    const double pap = final_tree(sh, part);
    if (threadIdx.x == 0) {
        c->pap = pap;
        if (!(pap > 0.0 && isfinite(pap))) { c->error = 1; return; }
        c->alpha = c->rz / pap;
        hist_alpha[c->corr * c->maxit + c->it] = c->alpha;
    }
}
extern "C" __global__ void k_update(int n, double* __restrict__ x, double* __restrict__ r, const double* __restrict__ p,
                                    const double* __restrict__ Ap, const double* __restrict__ diag, double* __restrict__ z,
                                    double* part_rz, double* part_rr, const Ctl* c) {
    if (stopped(c)) return;
    __shared__ double sh[TPB];
    const double alpha = c->alpha;
    double s1 = 0.0, s2 = 0.0;
    ROWS(i, n) {
        x[i] = x[i] + alpha * p[i];
        const double ri = r[i] - alpha * Ap[i];
        r[i] = ri;
        const double zi = ri / diag[i];
        z[i] = zi;
        s1 += ri * zi; s2 += ri * ri;
    }
    double v = block_tree(sh, s1);
    if (threadIdx.x == 0) part_rz[blockIdx.x] = v;
    __syncthreads();
    v = block_tree(sh, s2);
    if (threadIdx.x == 0) part_rr[blockIdx.x] = v;
}
extern "C" __global__ void k_fin_beta(const double* part_rz, const double* part_rr, Ctl* c, double* hist_beta, double* hist_rnorm) {
    if (stopped(c)) return;
    __shared__ double sh[NBLK];
    const double nrz = final_tree(sh, part_rz);
    __syncthreads();
    const double rr = final_tree(sh, part_rr);
    if (threadIdx.x == 0) {
        const double beta = nrz / c->rz;
        c->beta = beta; c->rz = nrz; c->rr = rr; c->rnorm = sqrt(rr);
        hist_beta[c->corr * c->maxit + c->it] = beta;
        const int it = c->it + 1;
        c->it = it;
        hist_rnorm[c->corr * (c->maxit + 1) + it] = c->rnorm;
        c->done = (it < c->maxit && c->rnorm > c->target * .2) ? 0 : 1;
    }
}
extern "C" __global__ void k_direction(int n, const double* __restrict__ z, double* __restrict__ p, const Ctl* c) {
    if (stopped(c)) return;
    const double beta = c->beta;
    const int i = blockIdx.x * TPB + threadIdx.x;
    if (i < n) p[i] = z[i] + beta * p[i];
}
extern "C" __global__ void k_residual(int n, const double* __restrict__ b, const double* __restrict__ Ax, double* __restrict__ r, double* part) {
    __shared__ double sh[TPB];
    double s = 0.0;
    ROWS(i, n) {
        const double ri = b[i] - Ax[i];
        r[i] = ri;
        s += ri * ri;
    }
    const double v = block_tree(sh, s);
    if (threadIdx.x == 0) part[blockIdx.x] = v;
}
// sigma(T) update: only rows with conduction change (CSR diagonal, Jacobi
// diagonal and a = (.5 dt)(S+Y)); values are computed on the host.
extern "C" __global__ void k_refresh(int m, const int* __restrict__ rows, const int* __restrict__ pos,
                                     const double* __restrict__ new_diag, const double* __restrict__ new_a,
                                     double* __restrict__ val, double* __restrict__ diag, double* __restrict__ a) {
    const int k = blockIdx.x * TPB + threadIdx.x;
    if (k >= m) return;
    const int e = rows[k];
    val[pos[k]] = new_diag[k];
    diag[e] = new_diag[k];
    a[e] = new_a[k];
}
extern "C" __global__ void k_fin_true(const double* part, Ctl* c) {
    __shared__ double sh[NBLK];
    const double v = final_tree(sh, part);
    if (threadIdx.x == 0) c->tnorm = sqrt(v);
}
)CUDA";
}  // namespace hgf
