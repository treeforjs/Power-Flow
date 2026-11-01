// gas_kernels.cu -- gas target owner v1 device kernels (NVRTC; gas_core.h and
// forms_core.h compiled from the same bytes as the host reference). Exact
// accumulation by integer digit atomics (order independent); every binary64
// operation is the gas_core.h function of the host path.
#include "gas_core.h"

using pfg::i64;
using pfg::u32;
using pfg::u64;

__device__ inline void add_atomic(i64 *d, const pfg::Chunks &t, int negate) {
  if (t.zero || t.bad)
    return;
  const int neg = t.neg ^ negate;
  for (int k = 0; k < t.n; ++k)
    if (t.c[k])
      atomicAdd((unsigned long long *)(d + t.d + k),
                (unsigned long long)(neg ? -(long long)t.c[k] : (long long)t.c[k]));
}

struct DevSink {
  i64 *e;
  i64 *g;
  unsigned long long *c;
  __device__ void entry(i64 idx, u32 q, const pfg::Chunks &x, int negate) {
    add_atomic(e + ((u64)idx * pfg::kEntryQ + q) * pfg::kDigits, x, negate);
  }
  __device__ void global(u32 idx, const pfg::Chunks &x, int negate) { add_atomic(g + (u64)idx * pfg::kDigits, x, negate); }
  __device__ void status(u64 bits) { atomicOr(c + pfg::C_STATUS, (unsigned long long)bits); }
  __device__ void closure(double relE, double relP) {
    if (!(relE <= 0x1p-40))
      atomicAdd(c + pfg::C_CLOSURE_E, 1ull);
    if (!(relP <= 0x1p-40))
      atomicAdd(c + pfg::C_CLOSURE_P, 1ull);
    atomicMax(c + pfg::C_MAXREL_E, (unsigned long long)pf7::bits_of(relE));
    atomicMax(c + pfg::C_MAXREL_P, (unsigned long long)pf7::bits_of(relP));
  }
};

extern "C" __global__ void gas_apply_records(pfg::Tables T, const pfg::RecoilRecord *r, u64 n, i64 *de, i64 *dg,
                                             unsigned long long *dc) {
  const u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  DevSink S{de, dg, dc};
  pfg::apply_record(S, T, r[i]);
}
extern "C" __global__ void gas_apply_flows(pfg::Tables T, const pfg::FlowRecord *f, u64 n, i64 *de, i64 *dg,
                                           unsigned long long *dc) {
  const u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  DevSink S{de, dg, dc};
  pfg::apply_flow(S, T, f[i]);
}
extern "C" __global__ void gas_apply_movers(pfg::Tables T, const pfg::MoverRecord *m, u64 n, const double *value,
                                            i64 *de, i64 *dg, unsigned long long *dc) {
  const u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  DevSink S{de, dg, dc};
  pfg::apply_mover(S, T, m[i], value);
}
// dst[i] += src[i]
extern "C" __global__ void gas_add_digits(i64 *dst, const i64 *src, u64 n) {
  const u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    dst[i] += src[i];
}
extern "C" __global__ void gas_copy_digits(i64 *dst, const i64 *src, u64 n) {
  const u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n)
    dst[i] = src[i];
}
// New key table: entry block of new key k comes from old block from[k] (-1: zero).
extern "C" __global__ void gas_rekey(i64 *dst, const i64 *src, const long long *from, u64 nkeys) {
  const u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  const u64 per = (u64)pfg::kEntryQ * pfg::kDigits;
  if (i >= nkeys * per)
    return;
  const u64 k = i / per, j = i % per;
  const long long f = from[k];
  dst[i] = f < 0 ? 0 : src[(u64)f * per + j];
}
// Normalise accumulators; entries also refresh their rounded value.
extern "C" __global__ void gas_normalise(i64 *d, double *value, u64 n, unsigned long long *dc) {
  const u64 a = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (a >= n)
    return;
  i64 *x = d + a * pfg::kDigits;
  pf7::normalise(x);
  if (value) {
    i64 y[pfg::kDigits];
    for (int k = 0; k < pfg::kDigits; ++k)
      y[k] = x[k];
    pf7_u32 st = 0;
    value[a] = pf7::round_sum(y, &st);
    if (st)
      atomicOr(dc + pfg::C_STATUS, (unsigned long long)pfg::kStOverflow);
  }
}
// Feedback: one thread per cell. out: eligible, changed, clamped.
extern "C" __global__ void gas_feedback(const u32 *atlas, pfg::CellTargetsAbi *cells, u64 ncells, const u64 *keys,
                                        u64 nkeys, const double *value, const pfg::SpeciesInfo *species, u32 S,
                                        double V, unsigned long long *out) {
  const u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= ncells)
    return;
  pfg::CellTargetsAbi &c = cells[i];
  for (u32 k = 0; k < c.n_targets && k < pfg::kMaxTargets; ++k) {
    const u32 sp = c.species[k];
    if (sp >= S || !(species[sp].flags & pfg::kSpFeedback))
      continue;
    atomicAdd(out + 0, 1ull);
    const i64 e = pfg::find_key(keys, nkeys, pfg::entry_key(atlas[i], sp));
    if (e < 0)
      continue;
    int clamped = 0;
    if (pfg::feedback_slot(c.n[k], V, species[sp].mass_kg, value + (u64)e * pfg::kEntryQ, c.drift[k], &c.s[k],
                           &clamped)) {
      atomicAdd(out + 1, 1ull);
      if (clamped)
        atomicAdd(out + 2, 1ull);
    }
  }
}
// Leaves: 8 words per accumulator block and a nonzero flag.
extern "C" __global__ void gas_entry_leaves(const i64 *d, const u64 *keys, u64 n, unsigned int *leaf,
                                            unsigned int *nonzero) {
  const u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  const i64 *x = d + i * pfg::kEntryQ * pfg::kDigits;
  nonzero[i] = pfg::digits_zero(x, (int)(pfg::kEntryQ * pfg::kDigits)) ? 0u : 1u;
  if (nonzero[i])
    pfg::entry_leaf(keys[i], x, leaf + 8 * i);
}
extern "C" __global__ void gas_global_leaves(const i64 *d, u64 n, unsigned int *leaf, unsigned int *nonzero) {
  const u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  const i64 *x = d + i * pfg::kDigits;
  nonzero[i] = pfg::digits_zero(x, pfg::kDigits) ? 0u : 1u;
  if (nonzero[i])
    pfg::global_leaf(i, x, leaf + 8 * i);
}
