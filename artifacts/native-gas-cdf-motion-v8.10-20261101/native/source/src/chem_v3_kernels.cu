// chem_v3_kernels.cu -- NVRTC program of the native v8 chemistry v3 device
// stage (API-008.txt 3.7). Compiled by the shared GPU layer with --fmad=false
// --prec-div=true --prec-sqrt=true --ftz=false from the sealed chemistry v3
// headers (pfc::mcc_carrier, the C3 per-carrier collision stage) and the
// forms core (exact superaccumulator, SHA-256). Every value written equals the
// host path of chem_v3_native.cpp bit for bit:
//   classification (species from the census state index, binary64 grid
//   locate), projectile and cell compaction in ascending native ID / atlas ID,
//   exact per-(cell, species) target sums (integer digits: any order),
//   CellTargets, the C3 stage, per-carrier totals (exact coverage sums by
//   integer digit atomics, counts, max), 256-carrier record leaves, event and
//   active-carrier compaction (canonical order), row edits and compaction.
#include "forms_core.h"
#include "chem_event.h"

using namespace pf7;
#define CH_NONE 0xFFFFFFFFu
#define CH_S 22u
typedef unsigned long long ull;

// ------------------------------------------------------------ scans (u32 -> u64 exclusive)
// Block-local inclusive scan of 256 values; block totals out.
extern "C" __global__ void ch_scan_blocks(const pf7_u32 *in, pf7_u64 n, pf7_u64 *out, pf7_u64 *totals) {
  __shared__ pf7_u64 s[256];
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  s[threadIdx.x] = i < n ? in[i] : 0;
  __syncthreads();
  for (unsigned o = 1; o < 256; o <<= 1) {
    const pf7_u64 v = threadIdx.x >= o ? s[threadIdx.x - o] : 0;
    __syncthreads();
    s[threadIdx.x] += v;
    __syncthreads();
  }
  if (i < n)
    out[i] = s[threadIdx.x] - (i < n ? in[i] : 0); // exclusive within the block
  if (threadIdx.x == 255)
    totals[blockIdx.x] = s[255];
}
// Exclusive scan of the block totals (one thread; blocks <= a few 10^4).
extern "C" __global__ void ch_scan_totals(pf7_u64 *totals, pf7_u64 blocks, pf7_u64 *grand) {
  if (threadIdx.x || blockIdx.x)
    return;
  pf7_u64 acc = 0;
  for (pf7_u64 b = 0; b < blocks; ++b) {
    const pf7_u64 v = totals[b];
    totals[b] = acc;
    acc += v;
  }
  *grand = acc;
}
extern "C" __global__ void ch_scan_add(pf7_u64 *out, pf7_u64 n, const pf7_u64 *totals) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (i < n)
    out[i] += totals[blockIdx.x];
}
extern "C" __global__ void ch_zero_words(pf7_u64 *p, pf7_u64 n) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (i < n)
    p[i] = 0;
}
extern "C" __global__ void ch_fill_u32(pf7_u32 *p, pf7_u64 n, pf7_u32 v) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (i < n)
    p[i] = v;
}

// ------------------------------------------------------------ classification
struct ClassArgs {
  pf7_u64 n;
  const pf7_u64 *rows;   // 10 words per row
  const pf7_u64 *states; // census state index per row
  const pf7_u64 *keys;   // packet key per row
  const pf7_u32 *state_species; // census state index -> species (CH_NONE)
  pf7_u64 n_states;
  const pf7_u32 *proj, *targ;   // per species flags
  double low[3], high[3];
  pf7_u32 cnt[3], pad;
  const pf7_u32 *vacuum;        // per grid cell
  pf7_u32 *species, *cell, *flags, *pflag;
  pf7_u32 *occupied;            // per grid cell
  ull *err;                     // [4]: smallest failing row per code (grid, packet, id, cell targets)
};
__device__ pf7_u32 ch_locate(const ClassArgs &a, const double *x) {
  pf7_u64 idx[3];
  for (int d = 0; d < 3; ++d) {
    if (!(x[d] >= a.low[d] && x[d] <= a.high[d]))
      return CH_NONE;
    const double f = (x[d] - a.low[d]) / (a.high[d] - a.low[d]);
    pf7_u64 i = (pf7_u64)(f * (double)a.cnt[d]);
    if (i >= a.cnt[d])
      i = a.cnt[d] - 1;
    idx[d] = i;
  }
  const pf7_u64 id = idx[0] + (pf7_u64)a.cnt[0] * (idx[1] + (pf7_u64)a.cnt[1] * idx[2]);
  return a.vacuum[id] ? (pf7_u32)id : CH_NONE;
}
extern "C" __global__ void ch_classify(ClassArgs a) {
  const pf7_u64 r = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (r >= a.n)
    return;
  const pf7_u64 *row = a.rows + 10 * r;
  const pf7_u64 st = a.states[r];
  const pf7_u32 sp = st < a.n_states ? a.state_species[st] : CH_NONE;
  pf7_u32 flags = 0, cell = CH_NONE;
  if (sp != CH_NONE) {
    double x[3] = {double_of(row[1]), double_of(row[2]), double_of(row[3])};
    cell = ch_locate(a, x);
    if (cell == CH_NONE)
      atomicMin(&a.err[0], (ull)r);
    else {
      if (a.proj[sp])
        flags |= 1u;
      if (a.targ[sp])
        flags |= 2u;
    }
  }
  a.species[r] = sp;
  a.cell[r] = cell;
  a.flags[r] = flags;
}
// Rows of multi-carrier packets (key != own ID marks the secondary row and
// its packet's first row).
extern "C" __global__ void ch_multi(ClassArgs a) {
  const pf7_u64 r = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (r >= a.n)
    return;
  const pf7_u64 id = a.rows[10 * r], key = a.keys[r];
  if (key == id)
    return;
  atomicOr(&a.flags[r], 4u);
  pf7_u64 lo = 0, hi = a.n;
  while (lo < hi) {
    const pf7_u64 mid = lo + (hi - lo) / 2;
    if (a.rows[10 * mid] < key)
      lo = mid + 1;
    else
      hi = mid;
  }
  if (lo < a.n && a.rows[10 * lo] == key)
    atomicOr(&a.flags[lo], 4u);
}
extern "C" __global__ void ch_mark(ClassArgs a) {
  const pf7_u64 r = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (r >= a.n)
    return;
  const pf7_u32 f = a.flags[r];
  a.pflag[r] = f & 1u;
  if (!(f & 1u))
    return;
  if (f & 4u)
    atomicMin(&a.err[1], (ull)r);
  if (a.rows[10 * r] >= 0x1000000000000ull)
    atomicMin(&a.err[2], (ull)r);
  a.occupied[a.cell[r]] = 1u;
}
extern "C" __global__ void ch_cells(const pf7_u32 *occupied, const pf7_u64 *rank, pf7_u64 n_cells, pf7_u32 *cell_atlas) {
  const pf7_u64 c = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (c < n_cells && occupied[c])
    cell_atlas[rank[c]] = (pf7_u32)c;
}
extern "C" __global__ void ch_carriers(ClassArgs a, const pf7_u64 *prank, const pf7_u64 *crank, pfc::CarrierIn *out,
                                       pf7_u32 *carrier_row) {
  const pf7_u64 r = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (r >= a.n || !(a.flags[r] & 1u))
    return;
  const pf7_u64 *row = a.rows + 10 * r;
  const pf7_u64 k = prank[r];
  pfc::CarrierIn c;
  c.id = row[0];
  c.u[0] = double_of(row[4]);
  c.u[1] = double_of(row[5]);
  c.u[2] = double_of(row[6]);
  c.weight = double_of(row[9]);
  c.species = a.species[r];
  c.cell = (pf7_u32)crank[a.cell[r]];
  out[k] = c;
  carrier_row[k] = (pf7_u32)r;
}
// ------------------------------------------------------------ targets
extern "C" __global__ void ch_seg_count(ClassArgs a, const pf7_u64 *crank, pf7_u32 *seg_count) {
  const pf7_u64 r = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (r >= a.n || !(a.flags[r] & 2u) || !a.occupied[a.cell[r]])
    return;
  atomicAdd(&seg_count[crank[a.cell[r]] * CH_S + a.species[r]], 1u);
}
extern "C" __global__ void ch_seg_fill(ClassArgs a, const pf7_u64 *crank, const pf7_u64 *seg_off, pf7_u32 *seg_fill,
                                       pf7_u32 *seg_rows) {
  const pf7_u64 r = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (r >= a.n || !(a.flags[r] & 2u) || !a.occupied[a.cell[r]])
    return;
  const pf7_u64 seg = crank[a.cell[r]] * CH_S + a.species[r];
  const pf7_u32 slot = atomicAdd(&seg_fill[seg], 1u);
  seg_rows[seg_off[seg] + slot] = (pf7_u32)r;
}
struct SegArgs {
  pf7_u64 n_seg; // n_occ * CH_S
  const pf7_u64 *rows;
  const pf7_u32 *seg_count, *seg_rows;
  const pf7_u64 *seg_off, *seg_rank; // exclusive offsets; rank among present segments
  double volume;
  pfc::CellTargets *cells;
  pf7_i64 *inv;   // present segments x acc_digits: exact sum of weights (normalised digits)
  ull *err;
};
__device__ double ch_round(pf7_i64 *d) {
  pf7_u32 st = 0;
  return round_sum(d, &st);
}
extern "C" __global__ void ch_seg_sum(SegArgs a) {
  const pf7_u64 seg = (pf7_u64)blockIdx.x * 64 + threadIdx.x;
  if (seg >= a.n_seg)
    return;
  const pf7_u32 cnt = a.seg_count[seg];
  if (!cnt)
    return;
  const pf7_u64 c = seg / CH_S;
  const pf7_u32 sp = (pf7_u32)(seg % CH_S);
  pf7_u32 slot = 0;
  for (pf7_u32 s = 0; s < sp; ++s)
    slot += a.seg_count[c * CH_S + s] ? 1u : 0u;
  if (slot >= PFC_MAX_TARGETS) {
    atomicMin(&a.err[3], (ull)c);
    return;
  }
  pf7_i64 sw[acc_digits], su[acc_digits];
  for (int k = 0; k < acc_digits; ++k)
    sw[k] = 0;
  const pf7_u32 *rows = a.seg_rows + a.seg_off[seg];
  for (pf7_u32 j = 0; j < cnt; ++j)
    add_term(sw, product_term(double_of(a.rows[10 * (pf7_u64)rows[j] + 9]), 1.0));
  normalise(sw);
  pf7_i64 *inv = a.inv + a.seg_rank[seg] * acc_digits;
  for (int k = 0; k < acc_digits; ++k)
    inv[k] = sw[k];
  const double W = ch_round(sw);
  double drift[3];
  for (int d = 0; d < 3; ++d) {
    for (int k = 0; k < acc_digits; ++k)
      su[k] = 0;
    for (pf7_u32 j = 0; j < cnt; ++j) {
      const pf7_u64 *r = a.rows + 10 * (pf7_u64)rows[j];
      add_term(su, product_term(double_of(r[9]), double_of(r[4 + d])));
    }
    drift[d] = ch_round(su) / W;
  }
  for (int k = 0; k < acc_digits; ++k)
    su[k] = 0;
  for (pf7_u32 j = 0; j < cnt; ++j) {
    const pf7_u64 *r = a.rows + 10 * (pf7_u64)rows[j];
    const double t0 = double_of(r[4]) - drift[0], t1 = double_of(r[5]) - drift[1], t2 = double_of(r[6]) - drift[2];
    const double q2 = t0 * t0 + t1 * t1 + t2 * t2;
    add_term(su, product_term(double_of(r[9]), q2));
  }
  const double S2 = ch_round(su);
  pfc::CellTargets &ct = a.cells[c];
  ct.species[slot] = sp;
  ct.n[slot] = W / a.volume;
  for (int d = 0; d < 3; ++d)
    ct.drift[slot][d] = drift[d];
  ct.s[slot] = sqrt((S2 / W) / 3.0);
}
extern "C" __global__ void ch_cell_count(const pf7_u32 *seg_count, pf7_u64 n_occ, pfc::CellTargets *cells) {
  const pf7_u64 c = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (c >= n_occ)
    return;
  pf7_u32 n = 0;
  for (pf7_u32 s = 0; s < CH_S; ++s)
    n += seg_count[c * CH_S + s] ? 1u : 0u;
  cells[c].n_targets = n < PFC_MAX_TARGETS ? n : PFC_MAX_TARGETS;
}
// ------------------------------------------------------------ the C3 stage (one chunk)
extern "C" __global__ void ch_mcc(pfc::TableView t, pfc::CoverView cv, pfc::StepConfig cfg, const pfc::CarrierIn *carriers,
                                  const pfc::CellTargets *cells, unsigned n_cells, unsigned begin, unsigned count,
                                  pfc::EventRec *events, pfc::CarrierOut *outs) {
  const unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= count)
    return;
  const pfc::CarrierIn c = carriers[begin + j];
  if (c.cell >= n_cells) {
    pfc::CarrierOut o;
    o.id = c.id;
    o.status = PFC_ST_BAD_INPUT;
    o.n_events = 0;
    o.n_candidates = 0;
    o.species_final = c.species;
    o.u_final[0] = c.u[0];
    o.u_final[1] = c.u[1];
    o.u_final[2] = c.u[2];
    o.nu_dt_max = 0.0;
    o.n_cover = 0;
    o.pad = 0;
    for (int k = 0; k < PFC_MAX_COVER_SLOTS; ++k) {
      o.cover_region[k] = 0xFFFFFFFFu;
      o.cover_exposure[k] = 0.0;
      o.cover_bound[k] = 0.0;
    }
    outs[j] = o;
    return;
  }
  pfc::mcc_carrier(t, cv, cfg, c, cells[c.cell], events + (unsigned long long)j * cfg.max_events, outs + j);
}
struct SumArgs {
  unsigned begin, count;
  const pfc::CarrierOut *outs;
  pf7_u32 *nev, *active;       // per projectile: events; 1 if events
  ull *bad;                    // minimum (uint32 projectile index << 32) | actual bounded status
  ull *candidates, *status;    // status[9] (includes selected unsupported CDF)
  ull *numax;                  // bits of max nu_dt_max
  pf7_i64 *region;             // regions x 2 x acc_digits (exposure, bound)
  ull *support;                // reachable carrier/region pairs; n_proj bounded by uint32
};
__device__ void ch_atomic_term(pf7_i64 *digit, const Term &t) {
  if (t.zero || t.bad)
    return;
  for (int k = 0; k < 5; ++k)
    if (t.c[k])
      atomicAdd((ull *)&digit[t.d + k], (ull)(t.neg ? -(pf7_i64)t.c[k] : (pf7_i64)t.c[k]));
}
extern "C" __global__ void ch_summary(SumArgs a) {
  const unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= a.count)
    return;
  const unsigned i = a.begin + j;
  const pfc::CarrierOut &o = a.outs[j];
  a.nev[i] = o.n_events;
  a.active[i] = o.n_events ? 1u : 0u;
  if (o.status != PFC_ST_OK && o.status != PFC_ST_CONSUMED)
    atomicMin(a.bad, ((ull)i << 32) | (ull)(o.status <= PFC_ST_CDF_SUPPORT ? o.status : PFC_ST_BAD_INPUT));
  atomicAdd(a.candidates, (ull)o.n_candidates);
  atomicAdd(&a.status[o.status <= PFC_ST_CDF_SUPPORT ? o.status : PFC_ST_BAD_INPUT], 1ull);
  atomicMax(a.numax, (ull)bits_of(o.nu_dt_max));
  for (pf7_u32 k = 0; k < o.n_cover && k < PFC_MAX_COVER_SLOTS; ++k) {
    atomicAdd(&a.support[o.cover_region[k]], 1ull);
    pf7_i64 *base = a.region + (pf7_u64)o.cover_region[k] * 2 * acc_digits;
    ch_atomic_term(base, product_term(o.cover_exposure[k], 1.0));
    ch_atomic_term(base + acc_digits, product_term(o.cover_bound[k], 1.0));
  }
}
// One thread per 256-carrier leaf of the chunk (begin is a multiple of 256):
// SHA256(text tag, u64 leaf, u64 count, records of each carrier).
extern "C" __global__ void ch_leaves(unsigned begin, unsigned count, unsigned max_events, const pfc::CarrierOut *outs,
                                     const pfc::EventRec *events, const char *tag, unsigned tag_len, pf7_u32 *leaves) {
  const unsigned l = blockIdx.x * blockDim.x + threadIdx.x;
  const unsigned lo = l * 256;
  if (lo >= count)
    return;
  const unsigned hi = lo + 256 < count ? lo + 256 : count;
  Sha s;
  sha_init(s);
  sha_text(s, tag, tag_len);
  sha_u64(s, (begin + lo) / 256);
  sha_u64(s, hi - lo);
  for (unsigned j = lo; j < hi; ++j) {
    sha_u64_array(s, (const pf7_u64 *)(outs + j), sizeof(pfc::CarrierOut) / 8);
    sha_u64_array(s, (const pf7_u64 *)(events + (unsigned long long)j * max_events),
                  (pf7_u64)outs[j].n_events * (sizeof(pfc::EventRec) / 8));
  }
  pf7_u32 out[8];
  sha_finish(s, out);
  for (int k = 0; k < 8; ++k)
    leaves[(pf7_u64)((begin + lo) / 256) * 8 + k] = out[k];
}
// Event compaction of one chunk: offsets = exclusive scan of nev within the
// chunk; base = running total (device counter, stream ordered).
extern "C" __global__ void ch_events_copy(unsigned begin, unsigned count, unsigned max_events, const pf7_u32 *nev,
                                          const pf7_u64 *offset, const ull *base, ull cap, const pfc::EventRec *events,
                                          pfc::EventRec *compact, ull *overflow) {
  const unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= count)
    return;
  const pf7_u32 n = nev[begin + j];
  if (!n)
    return;
  const ull at = *base + offset[j];
  if (at + n > cap) {
    atomicMax(overflow, 1ull);
    return;
  }
  for (pf7_u32 k = 0; k < n; ++k)
    compact[at + k] = events[(unsigned long long)j * max_events + k];
}
extern "C" __global__ void ch_add_total(ull *base, const pf7_u64 *grand) {
  if (threadIdx.x == 0 && blockIdx.x == 0)
    *base += *grand;
}
// ------------------------------------------------------------ active carriers
struct ActiveRec { // 64 bytes
  pf7_u64 id;
  pf7_u32 index, row, cell, species, nev, status;
  double x[3];
  double weight;
};
extern "C" __global__ void ch_active(pf7_u64 n, const pf7_u32 *active, const pf7_u64 *arank, const pf7_u32 *nev,
                                     const pfc::CarrierIn *carriers, const pf7_u32 *carrier_row, const pf7_u64 *rows,
                                     ActiveRec *out) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (i >= n || !active[i])
    return;
  ActiveRec r;
  const pf7_u32 row = carrier_row[i];
  r.id = carriers[i].id;
  r.index = (pf7_u32)i;
  r.row = row;
  r.cell = carriers[i].cell;
  r.species = carriers[i].species;
  r.nev = nev[i];
  r.status = 0;
  for (int d = 0; d < 3; ++d)
    r.x[d] = double_of(rows[10 * (pf7_u64)row + 1 + d]);
  r.weight = double_of(rows[10 * (pf7_u64)row + 9]);
  out[arank[i]] = r;
}
// Inventory digits of requested present segments.
extern "C" __global__ void ch_gather_inv(pf7_u64 m, const pf7_u64 *which, const pf7_i64 *inv, pf7_i64 *out) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (i >= m)
    return;
  for (int k = 0; k < acc_digits; ++k)
    out[i * acc_digits + k] = inv[which[i] * acc_digits + k];
}
// Rows of one segment: row, ID, key, flags, weight and six exact x/u words.
extern "C" __global__ void ch_gather_seg(pf7_u32 cnt, const pf7_u32 *seg_rows, const pf7_u64 *rows, const pf7_u64 *keys,
                                         const pf7_u32 *flags, pf7_u64 *out) {
  const pf7_u32 j = blockIdx.x * 256 + threadIdx.x;
  if (j >= cnt)
    return;
  const pf7_u32 r = seg_rows[j];
  out[11 * (pf7_u64)j] = r;
  out[11 * (pf7_u64)j + 1] = rows[10 * (pf7_u64)r];
  out[11 * (pf7_u64)j + 2] = keys[r];
  out[11 * (pf7_u64)j + 3] = flags[r];
  out[11 * (pf7_u64)j + 4] = rows[10 * (pf7_u64)r + 9];
  for (int d = 0; d < 6; ++d)
    out[11 * (pf7_u64)j + 5 + d] = rows[10 * (pf7_u64)r + 1 + d];
}
extern "C" __global__ void ch_nonzero(const pf7_u32 *in, pf7_u64 n, pf7_u32 *out) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (i < n)
    out[i] = in[i] ? 1u : 0u;
}
// Generic leaf hash of fixed-size records (8-byte multiples):
// SHA256(text tag, u64 leaf, u64 count, record words).
extern "C" __global__ void ch_record_leaves(pf7_u64 n, pf7_u64 words, pf7_u64 per_leaf, const pf7_u64 *data,
                                            const char *tag, unsigned tag_len, pf7_u32 *leaves) {
  const pf7_u64 l = (pf7_u64)blockIdx.x * blockDim.x + threadIdx.x;
  const pf7_u64 lo = l * per_leaf;
  if (lo >= n)
    return;
  const pf7_u64 hi = lo + per_leaf < n ? lo + per_leaf : n;
  Sha s;
  sha_init(s);
  sha_text(s, tag, tag_len);
  sha_u64(s, l);
  sha_u64(s, hi - lo);
  sha_u64_array(s, data + lo * words, (hi - lo) * words);
  pf7_u32 out[8];
  sha_finish(s, out);
  for (int k = 0; k < 8; ++k)
    leaves[l * 8 + k] = out[k];
}
extern "C" __global__ void ch_widen(pf7_u64 n, const pf7_u32 *in, pf7_u64 *out) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (i < n)
    out[i] = in[i];
}
// ------------------------------------------------------------ row edits
// keep[r] = 0 for dropped rows; velocity / weight edits applied in place on
// the copy; compaction writes kept rows (ascending) then the appended rows.
extern "C" __global__ void ch_edit_velocity(pf7_u64 m, const pf7_u64 *edits, pf7_u64 *rows) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (i >= m)
    return;
  const pf7_u64 r = edits[4 * i];
  for (int d = 0; d < 3; ++d)
    rows[10 * r + 4 + d] = edits[4 * i + 1 + d];
}
extern "C" __global__ void ch_edit_weight(pf7_u64 m, const pf7_u64 *edits, pf7_u64 *rows) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (i >= m)
    return;
  rows[10 * edits[2 * i] + 9] = edits[2 * i + 1];
}
extern "C" __global__ void ch_drop(pf7_u64 m, const pf7_u64 *list, pf7_u32 *keep) {
  const pf7_u64 i = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (i < m)
    keep[list[i]] = 0u;
}
extern "C" __global__ void ch_compact(pf7_u64 n, const pf7_u32 *keep, const pf7_u64 *rank, const pf7_u64 *rows,
                                      const pf7_u64 *keys, const pf7_u64 *states, pf7_u64 *orows, pf7_u64 *okeys,
                                      pf7_u64 *ostates) {
  const pf7_u64 r = (pf7_u64)blockIdx.x * 256 + threadIdx.x;
  if (r >= n || !keep[r])
    return;
  const pf7_u64 k = rank[r];
  for (int w = 0; w < 10; ++w)
    orows[10 * k + w] = rows[10 * r + w];
  okeys[k] = keys[r];
  ostates[k] = states[r];
}
