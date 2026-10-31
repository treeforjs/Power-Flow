// native_kernels.cu -- NVRTC program of the native v7.1 device lane (motion
// admission without births and sinks). One thread per census carrier repeats
// the host path's per-carrier checks of FullTimeAccess::motion_v7 in the same
// arithmetic (closed_position with native_particle_guard, segment_guard, path
// start/displacement/flight, resident-neutral bound, endpoint identity and
// species/count, co-located carriers of one logical packet) and writes the
// terms of the exact sums (kinetic change, motion-work-gross-v1 carrier
// scales). Any failing check marks the carrier; the host then re-runs the host
// path for the refusal text. blockDim.x == 256; no atomics.
#include "forms_core.h"

using namespace pf7;
#define PF7N_BLOCK 256

struct Geometry {
  double width, gap, short_z, feed_z;
  double lower[4], upper[4]; // left, right, bottom, top
  double shift[3];
  double dt, neutral_bound;
  pf7_u32 finite_foils, finite_axial;
};
struct AdmitArgs {
  pf7_u64 n;
  const pf7_u64 *base;  // census rows (10 words), ascending native ID
  const pf7_u64 *ends;  // endpoint rows (10 words), one per census carrier
  const pf7_u64 *paths; // path rows (8 words)
  const pf7_u64 *keys;  // packet key per census carrier
  Geometry g;
  // term arrays (n or 2n doubles)
  double *k_w, *k_d;    // kinetic change: w, d
  double *g_w, *g_d;    // G_K: w (charged) , |d|
  double *c_w, *c_k;    // C_K over 2n: w (charged) twice, K0 then K1
  double *q_a;          // Q: |q| (charged)
  pf7_u64 *block_fail;  // per block: smallest failing index (~0 none)
  pf7_u32 *block_code;
};
PF7_HD double dbl(pf7_u64 w) { return double_of(w); }
// native_particle_guard of a GapParticle given as words (mfem_pic_arena.cpp).
PF7_HD int particle_guard(const pf7_u64 *p) {
  double u2 = 0;
  for (int d = 0; d < 3; ++d)
    if (!finite(dbl(p[1 + d])))
      return 0;
  for (int d = 0; d < 3; ++d) {
    const double u = dbl(p[4 + d]);
    if (!finite(u))
      return 0;
    u2 += u * u;
  }
  const double q = dbl(p[7]), m = dbl(p[8]), w = dbl(p[9]);
  const double ke = kinetic_energy(m, dbl(p[4]), dbl(p[5]), dbl(p[6]));
  if (!(finite(u2) && finite(q) && finite(m) && m > 0 && finite(w) && w > 0 && finite(q * w) && finite(q / m) &&
        finite(ke * w)))
    return 0;
  if (!(q == 0 || q * w != 0))
    return 0;
  if (!(ke == 0 || ke * w != 0))
    return 0;
  return 1;
}
// closed_position (mfem_pic_full_time.inc).
PF7_HD int closed_position(const Geometry &g, const pf7_u64 *p) {
  if (!particle_guard(p))
    return 0;
  const double x = dbl(p[1]), y = dbl(p[2]), z = dbl(p[3]);
  if (!(x > 0 && x < g.width && y > 0 && y < g.gap))
    return 0;
  if (g.finite_axial && !(z > g.short_z && z < g.feed_z))
    return 0;
  if (g.finite_foils)
    for (int b = 0; b < 2; ++b) {
      const double *f = b ? g.upper : g.lower;
      if (x >= f[0] && x <= f[1] && y >= f[2] && y <= f[3])
        return 0;
    }
  return 1;
}
// segment_guard (mfem_pic_full_time.inc): std::max / std::min semantics.
PF7_HD int segment_guard(const Geometry &g, const double *a, const double *b) {
  if (!g.finite_foils)
    return 1;
  for (int k = 0; k < 2; ++k) {
    const double *f = k ? g.upper : g.lower;
    double enter = 0, leave = 1;
    int possible = 1;
    for (int d = 0; d < 2; ++d) {
      const double low = d ? f[2] : f[0], high = d ? f[3] : f[1], D = b[d] - a[d];
      if (D == 0) {
        if (a[d] < low || a[d] > high)
          possible = 0;
        continue;
      }
      double l = (low - a[d]) / D, h = (high - a[d]) / D;
      if (l > h) {
        const double t = l;
        l = h;
        h = t;
      }
      enter = (enter < l) ? l : enter;
      leave = (h < leave) ? h : leave;
    }
    if (!(!possible || enter > leave))
      return 0;
  }
  return 1;
}
PF7_HD pf7_u64 find_id(const pf7_u64 *rows, pf7_u64 words, pf7_u64 n, pf7_u64 id) {
  pf7_u64 lo = 0, hi = n;
  while (lo < hi) {
    const pf7_u64 mid = lo + (hi - lo) / 2;
    if (rows[words * mid] < id)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo;
}
// Failure codes (diagnostic; the host path produces the refusal text).
enum : pf7_u32 {
  fail_none = 0,
  fail_path_id = 1,
  fail_neutral_bound = 2,
  fail_path_start = 3,
  fail_endpoint_id = 4,
  fail_closed = 5,
  fail_species = 6,
  fail_flight = 7,
  fail_displacement = 8,
  fail_segment = 9,
  fail_colocated = 10
};
PF7_HD pf7_u32 check_carrier(const AdmitArgs &a, pf7_u64 i) {
  const pf7_u64 *old = a.base + 10 * i, *path = a.paths + 8 * i, *e = a.ends + 10 * i;
  const Geometry &g = a.g;
  if (path[0] != old[0])
    return fail_path_id;
  const double q = dbl(old[7]);
  if (q == 0 && !(g.dt <= g.neutral_bound))
    return fail_neutral_bound;
  double start[3], pos_old[3], pos_e[3];
  for (int d = 0; d < 3; ++d) {
    pos_old[d] = dbl(old[1 + d]);
    pos_e[d] = dbl(e[1 + d]);
    start[d] = dbl(path[1 + d]);
    if (!(start[d] == pos_old[d] + g.shift[d]))
      return fail_path_start;
  }
  if (e[0] != old[0])
    return fail_endpoint_id;
  if (!closed_position(g, e))
    return fail_closed;
  if (!(dbl(e[7]) == dbl(old[7]) && dbl(e[8]) == dbl(old[8]) && dbl(e[9]) == dbl(old[9])))
    return fail_species;
  if (!(dbl(path[7]) == g.dt))
    return fail_flight;
  for (int d = 0; d < 3; ++d)
    if (!(dbl(path[4 + d]) == (pos_e[d] + g.shift[d]) - start[d]))
      return fail_displacement;
  if (!segment_guard(g, pos_old, pos_e))
    return fail_segment;
  // Co-located carriers of one logical packet stay co-located (the key
  // carrier, the smallest native ID of the packet, is in the census).
  const pf7_u64 key = a.keys[i];
  if (key != old[0]) {
    const pf7_u64 k = find_id(a.base, 10, a.n, key);
    if (k >= a.n || a.base[10 * k] != key)
      return fail_colocated;
    const pf7_u64 *f = a.ends + 10 * k;
    for (int w = 1; w <= 6; ++w)
      if (!(dbl(f[w]) == dbl(e[w])))
        return fail_colocated;
  }
  return fail_none;
}
extern "C" __global__ void pf7n_admit(const AdmitArgs a) {
  __shared__ pf7_u64 best[PF7N_BLOCK];
  __shared__ pf7_u32 code[PF7N_BLOCK];
  const pf7_u64 i = (pf7_u64)blockIdx.x * PF7N_BLOCK + threadIdx.x;
  pf7_u32 c = fail_none;
  if (i < a.n) {
    c = check_carrier(a, i);
    // Terms (written for every carrier; a failing step is decided on the host).
    const pf7_u64 *old = a.base + 10 * i, *e = a.ends + 10 * i;
    const double w = dbl(old[9]), m = dbl(old[8]), q = dbl(old[7]);
    const double K0 = kinetic_energy(m, dbl(old[4]), dbl(old[5]), dbl(old[6]));
    const double d = kinetic_energy(dbl(e[8]), dbl(e[4]), dbl(e[5]), dbl(e[6])) - K0;
    const double K1 = K0 + d;
    const int charged = q != 0;
    a.k_w[i] = w;
    a.k_d[i] = d;
    a.g_w[i] = charged ? w : 0.0;
    a.g_d[i] = d < 0 ? -d : d;
    a.c_w[i] = charged ? w : 0.0;
    a.c_w[a.n + i] = charged ? w : 0.0;
    a.c_k[i] = K0;
    a.c_k[a.n + i] = K1;
    a.q_a[i] = charged ? (q < 0 ? -q : q) : 0.0;
  }
  best[threadIdx.x] = c ? i : ~0ull;
  code[threadIdx.x] = c;
  __syncthreads();
  for (int s = PF7N_BLOCK / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s && best[threadIdx.x + s] < best[threadIdx.x]) {
      best[threadIdx.x] = best[threadIdx.x + s];
      code[threadIdx.x] = code[threadIdx.x + s];
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    a.block_fail[blockIdx.x] = best[0];
    a.block_code[blockIdx.x] = code[0];
  }
}
// First failing carrier over the blocks: out[0] = index (~0 none), out[1] = code.
extern "C" __global__ void pf7n_fail_finish(pf7_u64 blocks, const pf7_u64 *block_fail, const pf7_u32 *block_code,
                                           pf7_u64 *out) {
  if (threadIdx.x || blockIdx.x)
    return;
  out[0] = ~0ull;
  out[1] = 0;
  for (pf7_u64 k = 0; k < blocks; ++k)
    if (block_fail[k] != ~0ull) {
      out[0] = block_fail[k];
      out[1] = block_code[k];
      return;
    }
}
// dst[i] = src[i] for n words.
extern "C" __global__ void pf7n_copy(pf7_u64 n, const pf7_u64 *src, pf7_u64 *dst) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * PF7N_BLOCK + threadIdx.x;
  if (i < n)
    dst[i] = src[i];
}
