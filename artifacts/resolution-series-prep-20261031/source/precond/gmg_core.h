// gmg_core.h -- per-row arithmetic of the deterministic geometric multigrid preconditioner (resolution-series prep).
//
// Compiled twice from this one file: by the host C++ compiler (-ffp-contract=off -fno-fast-math) for the CPU reference,
// and by NVRTC (--fmad=false --prec-div=true --prec-sqrt=true --ftz=false) inside gmg_kernels.cu. Every function is a
// fixed-order loop over one row (ascending column index of a CSR row), using only + - * /, so the host and device builds
// give the same bits for the same inputs. No atomics; reductions use the fixed two-stage tree of gmg_kernels.cu, which
// the CPU reference reproduces (gmg.hpp, tree_dot).
#pragma once
#if defined(__CUDACC_RTC__)
#define GMG_FN __device__ __forceinline__
#else
#define GMG_FN inline
#endif

// y = row i of (K x), ascending column order
GMG_FN double gmg_row_mul(int i, const int* rp, const int* ci, const double* v, const double* x) {
    double s = 0.0;
    const int end = rp[i + 1];
    for (int k = rp[i]; k < end; ++k) s = s + v[k] * x[ci[k]];
    return s;
}
// r_i = b_i - (K x)_i
GMG_FN double gmg_row_residual(int i, const int* rp, const int* ci, const double* v, const double* x, const double* b) {
    return b[i] - gmg_row_mul(i, rp, ci, v, x);
}
