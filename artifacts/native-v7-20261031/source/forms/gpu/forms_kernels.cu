// forms_kernels.cu -- NVRTC program of the native v7 producer forms.
// Compile with the hall-gpu-resident-v1 session (--fmad=false --prec-div=true
// --prec-sqrt=true --ftz=false, architecture from the device) together with
// the headers forms_core.h and deterministic_primitives.cuh. blockDim.x must
// be 256 for every kernel. No atomics: exact sums are integer digit sums in a
// fixed order; failure searches are fixed-tree minima; hashes are pure.
//
// Entry points (API-007.txt section 3):
//  (a) pf7_sum_products, pf7_kinetic_aligned, pf7_kinetic_split -> per-block
//      partial digits; pf7_sum_finish -> SumOut (value, dt-scaled value,
//      status, normalised digits)
//  (b) pf7_chunks -> pf7_blocks -> pf7_tree (section root)
//  (c) pf7_census_pack, pf7_census_check + pf7_census_finish,
//      pf7_census_split_ids, pf7_gather_hi, pf7_census_gather (ordering with
//      the primitives' stable radix sort)
#include "forms_core.h"
#include "deterministic_primitives.cuh"

using namespace pf7;
#define PF7_BLOCK 256
#define PF7_WARPS 8
#define PF7_FULL 0xffffffffu

struct SumOut {
  double value;      // RN(exact sum)
  double scaled;     // RN(scale x value)
  pf7_u32 status;    // pf7::sum_* bits
  pf7_u32 parts;     // partial accumulators merged
  pf7_u64 terms;     // terms of the launch
  pf7_i64 digit[acc_digits]; // normalised exact sum (interchange form)
};

// Per-warp shared accumulators; each lane builds one term, then the warp
// adds the 32 terms one after another (lanes 0..4 add the five chunks).
template <class F>
__device__ void sum_body(pf7_u64 n, F term_of, pf7_i64 *partial, pf7_u32 *status) {
  __shared__ pf7_i64 acc[PF7_WARPS][acc_digits];
  const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
  for (int k = threadIdx.x; k < PF7_WARPS * acc_digits; k += PF7_BLOCK)
    (&acc[0][0])[k] = 0;
  __syncthreads();
  pf7_u32 flags = n >= (1ull << 31) ? sum_too_many_terms : 0u;
  const pf7_u64 stride = (pf7_u64)gridDim.x * PF7_BLOCK;
  for (pf7_u64 base = (pf7_u64)blockIdx.x * PF7_BLOCK + 32u * warp; base < n; base += stride) {
    const pf7_u64 i = base + lane;
    Term t;
    t.zero = 1;
    t.bad = 0;
    t.d = 0;
    t.neg = 0;
    for (int k = 0; k < 5; ++k)
      t.c[k] = 0;
    if (i < n)
      t = term_of(i, flags);
    if (t.bad)
      flags |= sum_nonfinite_input;
    for (int j = 0; j < 32; ++j) {
      const int skip = __shfl_sync(PF7_FULL, t.zero | t.bad, j);
      if (skip)
        continue;
      const int d = __shfl_sync(PF7_FULL, t.d, j);
      const int neg = __shfl_sync(PF7_FULL, t.neg, j);
      const pf7_u32 c0 = __shfl_sync(PF7_FULL, t.c[0], j);
      const pf7_u32 c1 = __shfl_sync(PF7_FULL, t.c[1], j);
      const pf7_u32 c2 = __shfl_sync(PF7_FULL, t.c[2], j);
      const pf7_u32 c3 = __shfl_sync(PF7_FULL, t.c[3], j);
      const pf7_u32 c4 = __shfl_sync(PF7_FULL, t.c[4], j);
      if (lane < 5) {
        const pf7_u32 c = lane == 0 ? c0 : lane == 1 ? c1 : lane == 2 ? c2 : lane == 3 ? c3 : c4;
        acc[warp][d + lane] += neg ? -(pf7_i64)c : (pf7_i64)c;
      }
      __syncwarp();
    }
  }
  const int any = __syncthreads_or((int)flags);
  // Block flags: OR over threads (fixed result).
  __shared__ pf7_u32 block_flags[PF7_BLOCK];
  block_flags[threadIdx.x] = flags;
  __syncthreads();
  for (int k = threadIdx.x; k < acc_digits; k += PF7_BLOCK) {
    pf7_i64 s = 0;
    for (int w = 0; w < PF7_WARPS; ++w)
      s += acc[w][k];
    partial[(pf7_u64)blockIdx.x * acc_digits + k] = s;
  }
  if (threadIdx.x == 0) {
    pf7_u32 f = 0;
    if (any)
      for (int k = 0; k < PF7_BLOCK; ++k)
        f |= block_flags[k];
    status[blockIdx.x] = f;
  }
}

// Exact sum of a_i x b_i (strides in doubles).
extern "C" __global__ void pf7_sum_products(pf7_u64 n, const double *a, pf7_u64 a_stride, const double *b,
                                            pf7_u64 b_stride, pf7_i64 *partial, pf7_u32 *status) {
  sum_body(
      n, [&](pf7_u64 i, pf7_u32 &) { return product_term(a[i * a_stride], b[i * b_stride]); }, partial,
      status);
}
// Kinetic terms of carriers aligned with their end proper velocity:
// w_i x RN(KE(m_i, u_end_i) - KE(m_i, u_start_i)). Strides in doubles; each
// u array holds 3 consecutive components per carrier.
extern "C" __global__ void pf7_kinetic_aligned(pf7_u64 n, const double *w, pf7_u64 w_stride, const double *m,
                                               pf7_u64 m_stride, const double *u0, pf7_u64 u0_stride,
                                               const double *u1, pf7_u64 u1_stride, pf7_i64 *partial,
                                               pf7_u32 *status) {
  sum_body(
      n,
      [&](pf7_u64 i, pf7_u32 &) {
        const double *a = u0 + i * u0_stride, *b = u1 + i * u1_stride;
        double ua[3] = {a[0], a[1], a[2]}, ub[3] = {b[0], b[1], b[2]};
        return product_term(w[i * w_stride], kinetic_change(m[i * m_stride], ua, ub));
      },
      partial, status);
}
// Kinetic terms of survivors (endpoint records, 10 words) and sinks (sink
// records, 17 words; incident proper velocity), each matched by native ID to
// its census carrier (records of census_words words, ascending ID; carrier
// words 0..9 as the endpoint layout). Missing IDs set sum_census_missing.
extern "C" __global__ void pf7_kinetic_split(pf7_u64 ne, const pf7_u64 *endpoints, pf7_u64 ns,
                                             const pf7_u64 *sinks, pf7_u64 nc, const pf7_u64 *census,
                                             pf7_u64 census_words, pf7_i64 *partial, pf7_u32 *status) {
  sum_body(
      ne + ns,
      [&](pf7_u64 i, pf7_u32 &flags) {
        const pf7_u64 *r = i < ne ? endpoints + 10 * i : sinks + 17 * (i - ne);
        const pf7_u64 id = r[0];
        pf7_u64 lo = 0, hi = nc;
        while (lo < hi) {
          const pf7_u64 mid = lo + (hi - lo) / 2;
          if (census[mid * census_words] < id)
            lo = mid + 1;
          else
            hi = mid;
        }
        Term t;
        t.zero = 1;
        t.bad = 0;
        t.d = 0;
        t.neg = 0;
        for (int k = 0; k < 5; ++k)
          t.c[k] = 0;
        if (lo >= nc || census[lo * census_words] != id) {
          flags |= sum_census_missing;
          return t;
        }
        const pf7_u64 *c = census + lo * census_words;
        double ua[3], ub[3];
        const int u_end = i < ne ? 4 : 12;
        for (int k = 0; k < 3; ++k) {
          ua[k] = double_of(c[4 + k]);
          ub[k] = double_of(r[u_end + k]);
        }
        return product_term(double_of(c[9]), kinetic_change(double_of(c[8]), ua, ub));
      },
      partial, status);
}
// Merges the per-block partials in block order (integers: any order gives
// the same digits), normalises and rounds. One block.
extern "C" __global__ void pf7_sum_finish(pf7_u32 parts, const pf7_i64 *partial, const pf7_u32 *status,
                                          pf7_u64 terms, double scale, SumOut *out) {
  __shared__ pf7_i64 d[acc_digits];
  for (int k = threadIdx.x; k < acc_digits; k += PF7_BLOCK) {
    pf7_i64 s = 0;
    for (pf7_u32 p = 0; p < parts; ++p)
      s += partial[(pf7_u64)p * acc_digits + k];
    d[k] = s;
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    pf7_u32 st = 0;
    for (pf7_u32 p = 0; p < parts; ++p)
      st |= status[p];
    normalise(d);
    for (int k = 0; k < acc_digits; ++k)
      out->digit[k] = d[k];
    const double v = round_sum(d, &st);
    out->value = v;
    out->scaled = scale * v;
    const pf7_u64 sb = bits_of(out->scaled);
    if (((sb >> 52) & 0x7ff) == 0x7ff)
      st |= sum_overflow;
    out->status = st;
    out->parts = parts;
    out->terms = terms;
  }
}

// ---------------------------------------------------------------- (b) block lists
// One thread per chunk of 32 records (records as little-endian u64 words).
extern "C" __global__ void pf7_chunks(pf7_u32 section, pf7_u64 n, const pf7_u64 *records, pf7_u32 *chunk_out) {
  const pf7_u64 j = (pf7_u64)blockIdx.x * PF7_BLOCK + threadIdx.x;
  const pf7_u64 chunks = (n + chunk_records - 1) / chunk_records;
  if (j >= chunks)
    return;
  const pf7_u64 m = n - j * chunk_records < chunk_records ? n - j * chunk_records : chunk_records;
  const pf7_u64 rw = section_record_bytes(section) / 8;
  chunk_digest(section, j, m, records + j * chunk_records * rw, chunk_out + 8 * j);
}
// One thread per leaf of 1024 records (32 chunk digests).
extern "C" __global__ void pf7_blocks(pf7_u32 section, pf7_u64 n, const pf7_u32 *chunks, pf7_u32 *leaf_out) {
  const pf7_u64 k = (pf7_u64)blockIdx.x * PF7_BLOCK + threadIdx.x;
  const pf7_u64 leaves = (n + leaf_records - 1) / leaf_records;
  if (k >= leaves)
    return;
  const pf7_u64 nk = n - k * leaf_records < leaf_records ? n - k * leaf_records : leaf_records;
  const pf7_u64 q = (nk + chunk_records - 1) / chunk_records;
  block_digest(section, k, nk, q, chunks + 8 * leaf_chunks * k, leaf_out + 8 * k);
}
// Native v5 tree over the leaves and the section root. One block; a and b
// are ping-pong scratch of ceil(leaves / 2) digests each.
extern "C" __global__ void pf7_tree(pf7_u32 section, pf7_u64 n, pf7_u64 leaves, const pf7_u32 *leaf,
                                    pf7_u32 *a, pf7_u32 *b, pf7_u32 *root_out) {
  __shared__ pf7_u32 tree[8];
  if (leaves == 0) {
    if (threadIdx.x == 0)
      merkle_empty(tree);
  } else if (leaves == 1) {
    if (threadIdx.x < 8)
      tree[threadIdx.x] = leaf[threadIdx.x];
  } else {
    const pf7_u32 *src = leaf;
    pf7_u32 *dst = a;
    pf7_u64 count = leaves;
    while (count > 1) {
      const pf7_u64 up = (count + 1) / 2;
      for (pf7_u64 i = threadIdx.x; i < up; i += PF7_BLOCK) {
        if (2 * i + 1 < count)
          merkle_node(src + 16 * i, src + 16 * i + 8, dst + 8 * i);
        else
          for (int w = 0; w < 8; ++w)
            dst[8 * i + w] = src[16 * i + w];
      }
      __syncthreads();
      src = dst;
      dst = dst == a ? b : a;
      count = up;
    }
    if (threadIdx.x < 8)
      tree[threadIdx.x] = src[threadIdx.x];
  }
  __syncthreads();
  if (threadIdx.x == 0)
    section_root(section, n, leaves, tree, root_out);
}

// ---------------------------------------------------------------- (c) census
// Census records (12 words) from carrier rows (row_words words per row,
// carrier words 0..9 first) and per-row packet keys and state indices.
extern "C" __global__ void pf7_census_pack(pf7_u64 n, const pf7_u64 *rows, pf7_u64 row_words, const pf7_u64 *keys,
                                           const pf7_u64 *states, pf7_u64 *out) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * PF7_BLOCK + threadIdx.x;
  if (i >= n)
    return;
  for (int w = 0; w < 10; ++w)
    out[12 * i + w] = rows[row_words * i + w];
  out[12 * i + 10] = keys[i];
  out[12 * i + 11] = states[i];
}
// Per block: smallest failing record index (or ~0) and its code.
extern "C" __global__ void pf7_census_check(pf7_u64 n, const pf7_u64 *rec, pf7_u64 states, pf7_u64 *block_index,
                                            pf7_u32 *block_code) {
  __shared__ pf7_u64 best[PF7_BLOCK];
  __shared__ pf7_u32 code[PF7_BLOCK];
  const pf7_u64 i = (pf7_u64)blockIdx.x * PF7_BLOCK + threadIdx.x;
  const pf7_u32 c = i < n ? census_check(rec, n, i, states) : 0u;
  best[threadIdx.x] = c ? i : ~0ull;
  code[threadIdx.x] = c;
  __syncthreads();
  for (int s = PF7_BLOCK / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s && best[threadIdx.x + s] < best[threadIdx.x]) {
      best[threadIdx.x] = best[threadIdx.x + s];
      code[threadIdx.x] = code[threadIdx.x + s];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    block_index[blockIdx.x] = best[0];
    block_code[blockIdx.x] = code[0];
  }
}
// First failing record over the blocks: out[0] = index (~0 when none),
// out[1] = code. One thread.
extern "C" __global__ void pf7_census_finish(pf7_u64 blocks, const pf7_u64 *block_index, const pf7_u32 *block_code,
                                             pf7_u64 *out) {
  if (threadIdx.x || blockIdx.x)
    return;
  out[0] = ~0ull;
  out[1] = 0;
  for (pf7_u64 k = 0; k < blocks; ++k)
    if (block_index[k] != ~0ull) {
      out[0] = block_index[k];
      out[1] = block_code[k];
      return;
    }
}
// Ordering: split the native IDs (record word 0, rec_words per record) into
// u32 halves for the two stable radix passes (low half first).
extern "C" __global__ void pf7_census_split_ids(pf7_u64 n, const pf7_u64 *rec, pf7_u64 rec_words, pf7_u32 *lo,
                                                pf7_u32 *hi) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * PF7_BLOCK + threadIdx.x;
  if (i >= n)
    return;
  const pf7_u64 id = rec[rec_words * i];
  lo[i] = (pf7_u32)(id & 0xffffffffull);
  hi[i] = (pf7_u32)(id >> 32);
}
// key_out[i] = hi[perm[i]] (second pass keys in first-pass order).
extern "C" __global__ void pf7_gather_hi(pf7_u64 n, const pf7_u32 *hi, const pf7_u32 *perm, pf7_u32 *key_out) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * PF7_BLOCK + threadIdx.x;
  if (i < n)
    key_out[i] = hi[perm[i]];
}
// out record i = rec record perm[i].
extern "C" __global__ void pf7_census_gather(pf7_u64 n, const pf7_u64 *rec, pf7_u64 rec_words, const pf7_u32 *perm,
                                             pf7_u64 *out) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * PF7_BLOCK + threadIdx.x;
  if (i >= n)
    return;
  for (pf7_u64 w = 0; w < rec_words; ++w)
    out[rec_words * i + w] = rec[rec_words * perm[i] + w];
}
