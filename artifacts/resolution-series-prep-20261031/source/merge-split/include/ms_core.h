// ms_core.h -- conservative macroparticle merge/split transaction, shared core.
//
// Compiled twice from this one file: by the host C++ compiler (-ffp-contract=off
// -fno-fast-math) for the CPU reference, and by NVRTC (--fmad=false
// --prec-div=true --prec-sqrt=true --ftz=false) for the device kernels. Only
// IEEE binary64 + - * / sqrt, floor (exact), comparisons and 32/64-bit integer
// arithmetic are used, so both builds produce the same bits.
//
// Contents
//  * classification of one carrier: Hall mesh cell, kinetic energy (the native
//    pic_kinetic_energy_J formula, operation for operation), velocity-space bin;
//  * error-free transformations (TwoSum, Veltkamp/Dekker TwoProduct) and an exact
//    fixed-point superaccumulator in units of 2^-1074 (the ExactParticleCount unit);
//  * group formation of one (cell, species) segment (merge and split policy);
//  * merge kinematics of one group (two outputs conserving momentum and kinetic
//    energy in exact arithmetic) and the split check;
//  * the Philox stream layout of the transaction (Philox4x32-10 core of the
//    arrival owner v1, philox_core.h, unchanged).
#pragma once
#include "philox_core.h"

#if defined(__CUDACC_RTC__)
#define MS_FN __device__ inline
typedef unsigned int ms_u32;
typedef unsigned long long ms_u64;
typedef long long ms_i64;
MS_FN double ms_from_bits(ms_u64 b) { return __longlong_as_double((long long)b); }
MS_FN double ms_sqrt(double x) { return sqrt(x); }
MS_FN double ms_floor(double x) { return floor(x); }
#else
#include <bit>
#include <cmath>
#include <cstdint>
#define MS_FN inline
typedef std::uint32_t ms_u32;
typedef std::uint64_t ms_u64;
typedef std::int64_t ms_i64;
MS_FN double ms_from_bits(ms_u64 b) { return std::bit_cast<double>(b); }
MS_FN double ms_sqrt(double x) { return std::sqrt(x); }
MS_FN double ms_floor(double x) { return std::floor(x); }
#endif
MS_FN ms_u64 ms_bits(double x) { return pfa_bits(x); }

#define MS_C_LIGHT 299792458.
#define MS_MAX_SPECIES 4
#define MS_MAX_THRESH 63
#define MS_MAX_SPLIT_LEVELS 4      // at most 16 clones per split
#define MS_MAX_GROUP 64            // k_max bound (one thread per group, sequential sums)
#define MS_ACC_DIGITS 68           // 32-bit digits in int64 lanes: bits 0..2175 of units 2^-1074
#define MS_FLAG_MERGE 1u           // carrier may be a merge input
#define MS_FLAG_SPLIT 2u           // carrier may be split
#define MS_KIND_MERGE 1u
#define MS_KIND_SPLIT 2u
// Group status (deterministic skips; the inputs of a skipped group survive unchanged).
#define MS_OK 0u
#define MS_SKIP_NO_REAL_SOLUTION 1u   // E/W below the energy of |P/W| (q^2 < 0) or nonfinite
#define MS_SKIP_OPERAND_RANGE 2u      // an exact ledger product would under/overflow
#define MS_SKIP_RNG_EXHAUSTED 3u      // no accepted disk point in 512 trials
#define MS_STATUS_KINDS 4u
// Group flags
#define MS_GFLAG_CENTROID_FALLBACK 1u  // centroid not strictly inside the cell; first input position used
#define MS_GFLAG_ISOTROPIC 2u          // P/W == 0: isotropic output axis
// Philox purposes (stream layout merge-split-v1)
#define MS_PURPOSE_AZIMUTH 1u
#define MS_PURPOSE_ISOTROPIC 2u
#define MS_KEY_DOMAIN0 0x6d657267u     // "merg": key = run seed XOR domain (separates from arrival v1)
#define MS_KEY_DOMAIN1 0x73706c74u     // "splt"

struct MsSpecies {
    double mass_kg, charge_C, mc2_J, split_min_weight;
    double e_thresh[MS_MAX_THRESH + 1];   // ascending kinetic-energy bin edges (J), first n_thresh used
    ms_u32 n_thresh, merge_above, merge_target, split_below;
    ms_u32 split_target, k_min, k_max, split_levels;
};
struct MsParams {
    ms_u64 seed, step, base_id;
    ms_u32 n_species, nsub;               // nsub: direction sub-cells per cube face edge
    ms_u32 nx, ny, nz, kz_lo, kz_hi, cell_bits;   // kinetic cells: k in [kz_lo, kz_hi)
    ms_u32 sub[3], pad0;                  // spatial sub-bins per cell axis (1 = none); bounds relocation length
    double shift[3];                      // native -> field frame (Hall: (-w/2, -g/2, 0))
    double origin[3], h[3];               // Hall mesh origin and spacing (field frame)
    MsSpecies species[MS_MAX_SPECIES];
};
// One group, in canonical order (cell, species, bin, first input ID). 64 bytes.
struct MsGroup {
    double W;                              // exact sum of the input weights (binary64)
    ms_u32 start, end, k, kind;            // members: eligible carriers of [start, end) (sorted indices)
    ms_u32 cell, species, bin, index;      // index: group number within its (cell, species) segment
    ms_u32 status, n_out, out_offset, reloc_offset;
    ms_u32 gflags, pad;
};
// Merge kinematics result of one group (unused for splits). 72 bytes.
struct MsMergeResult { double c[3], ua[3], ub[3]; };
// Sorted carrier view (structure of arrays, canonical order).
struct MsView {
    const ms_u64* id;
    const ms_u32* species;
    const ms_u32* flags;
    const ms_u32* bin;
    const double* x[3];
    const double* u[3];
    const double* w;
    const double* ke;
};
struct MsAcc { ms_i64 d[MS_ACC_DIGITS]; };

// ---------------------------------------------------------------- arithmetic
// Kinetic energy per physical particle, operation for operation as native
// pic_kinetic_energy_J (mfem_gap_pic.cpp:116, norm2/gamma at :25-26):
//   m*norm2(u)/(gamma(u)+1), norm2 = u0*u0+u1*u1+u2*u2, gamma = sqrt(1+norm2/(c*c)).
MS_FN double ms_kinetic_energy(double m, double u0, double u1, double u2) {
    const double c = MS_C_LIGHT;
    const double n2 = u0 * u0 + u1 * u1 + u2 * u2;
    const double g = ms_sqrt(1 + n2 / (c * c));
    return m * n2 / (g + 1);
}
MS_FN double ms_abs(double x) { return x < 0 ? -x : x; }
MS_FN int ms_finite(double x) { return ((ms_bits(x) >> 52) & 0x7ffu) != 0x7ffu; }
MS_FN int ms_exponent_field(double x) { return (int)((ms_bits(x) >> 52) & 0x7ffu); }
// Knuth TwoSum: s + e == a + b exactly (no overflow).
MS_FN void ms_two_sum(double a, double b, double* s, double* e) {
    const double x = a + b, bv = x - a, av = x - bv, br = b - bv, ar = a - av;
    *s = x;
    *e = ar + br;
}
// Veltkamp split and Dekker TwoProduct (Shewchuk's Two_Product_Tail form):
// p + e == a * b exactly when ms_product_safe(a, b).
MS_FN void ms_split(double a, double* hi, double* lo) {
    const double c = 134217729.0 * a;   // 2^27 + 1
    const double abig = c - a;
    *hi = c - abig;
    *lo = a - *hi;
}
MS_FN void ms_two_product(double a, double b, double* p, double* e) {
    const double x = a * b;
    double ah, al, bh, bl;
    ms_split(a, &ah, &al);
    ms_split(b, &bh, &bl);
    const double err1 = x - ah * bh, err2 = err1 - al * bh, err3 = err2 - ah * bl;
    *p = x;
    *e = al * bl - err3;
}
// Exactness guard of ms_two_product: both operands zero or normal with
// |unbiased exponent| <= 990 and the exponent sum in [-900, 900]; the Dekker
// error term is then exact (needs >= -970) and the split cannot overflow.
MS_FN int ms_product_safe(double a, double b) {
    if (a == 0 || b == 0) return 1;
    const int ea = ms_exponent_field(a), eb = ms_exponent_field(b);
    if (ea == 0 || eb == 0 || ea == 0x7ff || eb == 0x7ff) return 0;
    const int sa = ea - 1023, sb = eb - 1023;
    if (sa > 990 || sb > 990 || sa < -990 || sb < -990) return 0;
    const int s = sa + sb;
    return s >= -900 && s <= 900;
}
// Exact three-factor product a*b*m as four binary64 terms (guards included).
MS_FN int ms_product3_safe(double a, double b, double m) {
    if (!ms_product_safe(a, b)) return 0;
    double p, e;
    ms_two_product(a, b, &p, &e);
    return ms_product_safe(p, m) && ms_product_safe(e, m);
}
// Double-word accumulation used only for the merge construction (P/W, E/W, centroid).
MS_FN void ms_dd_add(double* hi, double* lo, double a, double b) {
    double s, t;
    ms_two_sum(*hi, a, &s, &t);
    t = t + (*lo + b);
    ms_two_sum(s, t, hi, lo);
}
MS_FN void ms_dd_add_product(double* hi, double* lo, double a, double b) {
    double p, e;
    ms_two_product(a, b, &p, &e);
    ms_dd_add(hi, lo, p, e);
}

// ---------------------------------------------------------------- superaccumulator
// value = sum_i d[i] * 2^(32 i) units of 2^-1074. Adding a finite binary64 is exact;
// digits stay below 2^63 in magnitude for < 2^31 additions between normalizations.
MS_FN void ms_acc_zero(MsAcc* a) {
    for (int i = 0; i < MS_ACC_DIGITS; ++i) a->d[i] = 0;
}
MS_FN void ms_acc_add(MsAcc* a, double x) {
    const ms_u64 b = ms_bits(x);
    ms_u32 e = (ms_u32)((b >> 52) & 0x7ffu);
    ms_u64 m = b & 0xFFFFFFFFFFFFFull;
    if (e == 0 && m == 0) return;
    if (e) m |= 1ull << 52;
    else e = 1;
    const ms_u32 p = e - 1;            // units 2^-1074 of mantissa bit 0
    const ms_u32 i = p >> 5, s = p & 31u;
    const ms_u64 lo = (m & 0xffffffffull) << s;
    const ms_u64 hi = (m >> 32) << s;
    const ms_u64 t1 = (lo >> 32) + hi;
    const ms_i64 d0 = (ms_i64)(lo & 0xffffffffull), d1 = (ms_i64)(t1 & 0xffffffffull), d2 = (ms_i64)(t1 >> 32);
    if (b >> 63) { a->d[i] -= d0; a->d[i + 1] -= d1; a->d[i + 2] -= d2; }
    else { a->d[i] += d0; a->d[i + 1] += d1; a->d[i + 2] += d2; }
}
// Exact a*b*m (four binary64 terms); caller checked ms_product3_safe.
MS_FN void ms_acc_add_product3(MsAcc* acc, double a, double b, double m) {
    double p, e, p1, e1, p2, e2;
    ms_two_product(a, b, &p, &e);
    ms_two_product(p, m, &p1, &e1);
    ms_two_product(e, m, &p2, &e2);
    ms_acc_add(acc, p1);
    ms_acc_add(acc, e1);
    ms_acc_add(acc, p2);
    ms_acc_add(acc, e2);
}
MS_FN void ms_acc_add_product2(MsAcc* acc, double a, double b) {
    double p, e;
    ms_two_product(a, b, &p, &e);
    ms_acc_add(acc, p);
    ms_acc_add(acc, e);
}
// Carry propagation: digits 0..66 in [0, 2^32), digit 67 signed.
MS_FN void ms_acc_normalize(MsAcc* a) {
    for (int i = 0; i < MS_ACC_DIGITS - 1; ++i) {
        const ms_i64 r = (ms_i64)((ms_u64)a->d[i] & 0xffffffffull);
        const ms_i64 c = (a->d[i] - r) / 4294967296LL;
        a->d[i] = r;
        a->d[i + 1] += c;
    }
}

// ---------------------------------------------------------------- classification
// Hall point-deposit cell of one axis (census_charge.hpp deposit_charge,
// axial_particles.hpp point_deposit): v = (x - o)/h, clamp to [0, count],
// cell = min(floor(v), count - 1). open: o + cell h < x < o + (cell+1) h in cell units.
// frac: v - cell in [0, 1] (cell-local coordinate).
MS_FN int ms_axis_cell(double xf, double o, double h, ms_u32 count, ms_u32* cell, int* open, double* frac) {
    double v = (xf - o) / h;
    if (!(v >= -1e-10 && v <= (double)count + 1e-10)) return 0;
    if (v < 0) v = 0;
    if (v > (double)count) v = (double)count;
    ms_u32 c = (ms_u32)ms_floor(v);
    if (c > count - 1) c = count - 1;
    *cell = c;
    *open = (v > (double)c) && (v < (double)c + 1);
    *frac = v - (double)c;
    return 1;
}
// Cell id (k*ny + j)*nx + i (Hall Mesh::cells order). Returns 0 outside the
// kinetic cells. open: strictly inside the open cell box on every axis.
// spatial (optional): sub-bin (sz*sub1 + sy)*sub0 + sx of the cell-local coordinates.
MS_FN int ms_cell_of_sub(const MsParams* P, double x0, double x1, double x2, ms_u32* cell, int* open, ms_u32* spatial) {
    const double xf[3] = {x0 + P->shift[0], x1 + P->shift[1], x2 + P->shift[2]};
    const ms_u32 counts[3] = {P->nx, P->ny, P->nz};
    ms_u32 c[3], sb[3];
    int all_open = 1;
    for (int d = 0; d < 3; ++d) {
        int o = 0;
        double f = 0;
        if (!ms_axis_cell(xf[d], P->origin[d], P->h[d], counts[d], &c[d], &o, &f)) return 0;
        all_open = all_open && o;
        const ms_u32 ns = P->sub[d];
        ms_u32 k = (ms_u32)ms_floor(f * (double)ns);
        sb[d] = k >= ns ? ns - 1 : k;
    }
    if (c[2] < P->kz_lo || c[2] >= P->kz_hi) return 0;
    *cell = (c[2] * P->ny + c[1]) * P->nx + c[0];
    *open = all_open;
    if (spatial) *spatial = (sb[2] * P->sub[1] + sb[1]) * P->sub[0] + sb[0];
    return 1;
}
MS_FN int ms_cell_of(const MsParams* P, double x0, double x1, double x2, ms_u32* cell, int* open) {
    return ms_cell_of_sub(P, x0, x1, x2, cell, open, (ms_u32*)0);
}
MS_FN ms_u32 ms_energy_bin(const MsSpecies* s, double ke) {
    ms_u32 lo = 0, hi = s->n_thresh;   // number of edges <= ke
    while (lo < hi) {
        const ms_u32 mid = (lo + hi) >> 1;
        if (s->e_thresh[mid] <= ke) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}
MS_FN ms_u32 ms_sub_cell(double s, ms_u32 nsub) {
    const double v = (s + 1) * (0.5 * (double)nsub);
    ms_u32 i = v <= 0 ? 0u : (ms_u32)ms_floor(v);
    return i >= nsub ? nsub - 1 : i;
}
// Cube-sphere direction cell: dominant axis a (ties: lowest index), face 2a + (u_a < 0),
// nsub x nsub sub-cells of the two other components divided by |u_a|. u == 0 -> 0.
MS_FN ms_u32 ms_direction_bin(ms_u32 nsub, double u0, double u1, double u2) {
    if (u0 == 0 && u1 == 0 && u2 == 0) return 0;
    const double a0 = ms_abs(u0), a1 = ms_abs(u1), a2 = ms_abs(u2);
    int a = 0;
    double am = a0;
    if (a1 > am) { a = 1; am = a1; }
    if (a2 > am) { a = 2; am = a2; }
    const double ua = a == 0 ? u0 : (a == 1 ? u1 : u2);
    const ms_u32 face = 2u * (ms_u32)a + (ua < 0 ? 1u : 0u);
    const double s = (a == 0 ? u1 : u0) / am;
    const double t = (a == 2 ? u1 : u2) / am;
    return (face * nsub + ms_sub_cell(s, nsub)) * nsub + ms_sub_cell(t, nsub);
}
// bin = (spatial * (n_thresh + 1) + energy bin) * n_dir + direction cell,
// n_dir = 6 nsub^2, or 1 with nsub == 0 (no direction binning).
MS_FN ms_u32 ms_bin(const MsParams* P, const MsSpecies* s, ms_u32 spatial, double ke, double u0, double u1, double u2) {
    const ms_u32 e = spatial * (s->n_thresh + 1u) + ms_energy_bin(s, ke);
    if (P->nsub == 0) return e;
    return e * (6u * P->nsub * P->nsub) + ms_direction_bin(P->nsub, u0, u1, u2);
}
// Classification of one carrier; returns 0 for an invalid census row.
MS_FN int ms_classify(const MsParams* P, ms_u32 species, double x0, double x1, double x2, double u0, double u1,
                      double u2, double w, ms_u32* cell, ms_u32* bin, double* ke) {
    if (species >= P->n_species) return 0;
    if (!(ms_finite(x0) && ms_finite(x1) && ms_finite(x2) && ms_finite(u0) && ms_finite(u1) && ms_finite(u2) &&
          ms_finite(w) && w > 0))
        return 0;
    int open = 0;
    ms_u32 spatial = 0;
    if (!ms_cell_of_sub(P, x0, x1, x2, cell, &open, &spatial)) return 0;
    const MsSpecies* s = &P->species[species];
    const double k = ms_kinetic_energy(s->mass_kg, u0, u1, u2);
    if (!ms_finite(k)) return 0;
    *ke = k;
    *bin = ms_bin(P, s, spatial, k, u0, u1, u2);
    return 1;
}

// ---------------------------------------------------------------- group formation
// One (cell, species) segment [begin, end) of the canonical order (bin, ID inside).
// N = end - begin carriers of the species in the cell.
//  merge (N > merge_above): reduce toward merge_target. In each bin, a group starts
//    at an eligible carrier and takes the following eligible carriers of the bin
//    while k < min(k_max, excess - reduction + 2) and the running weight sum stays
//    exactly representable (TwoSum error 0): a group never needs a rounded split.
//    k >= k_min -> group (reduction += k - 2); otherwise the start carrier is left.
//    W/2 must be normal (W >= 2^-1021), else the group is not formed (refusal count).
//  split (0 < N < split_below): eligible carriers with w >= split_min_weight, in
//    canonical order, become 2^m clones of weight w 2^-m (exact: the result must be
//    normal) until the added count reaches split_target - N.
// groups == 0: count only. Returns the number of groups; *refused counts groups not
// formed for an unrepresentable weight.
MS_FN ms_u32 ms_form_segment(const MsParams* P, MsView v, ms_u32 begin, ms_u32 end, ms_u32 cell, ms_u32 species,
                             MsGroup* groups, ms_u32* refused) {
    const MsSpecies* s = &P->species[species];
    const ms_u32 N = end - begin;
    ms_u32 count = 0, rej = 0;
    if (N > s->merge_above) {
        const ms_u32 excess = N - s->merge_target;
        ms_u32 reduction = 0;
        ms_u32 i = begin;
        int stop = 0;
        while (i < end && !stop) {
            const ms_u32 b = v.bin[i];
            ms_u32 j = i;
            while (j < end && v.bin[j] == b) ++j;
            ms_u32 p = i;
            while (p < j) {
                if (reduction >= excess) { stop = 1; break; }
                ms_u32 kmax = excess - reduction + 2;
                if (kmax > s->k_max) kmax = s->k_max;
                if (kmax < s->k_min) { stop = 1; break; }
                if (!(v.flags[p] & MS_FLAG_MERGE)) { ++p; continue; }
                double S = v.w[p];
                ms_u32 k = 1, last = p, q = p + 1;
                while (q < j && k < kmax) {
                    if (!(v.flags[q] & MS_FLAG_MERGE)) { ++q; continue; }
                    double t, e;
                    ms_two_sum(S, v.w[q], &t, &e);
                    if (e != 0 || !ms_finite(t)) break;
                    S = t;
                    ++k;
                    last = q;
                    ++q;
                }
                if (k >= s->k_min) {
                    if (!(S >= 0x1p-1021)) { ++rej; ++p; continue; }
                    if (groups) {
                        MsGroup g;
                        g.W = S; g.start = p; g.end = last + 1; g.k = k; g.kind = MS_KIND_MERGE;
                        g.cell = cell; g.species = species; g.bin = b; g.index = count;
                        g.status = 0; g.n_out = 0; g.out_offset = 0; g.reloc_offset = 0; g.gflags = 0; g.pad = 0;
                        groups[count] = g;
                    }
                    ++count;
                    reduction += k - 2;
                    p = last + 1;
                } else {
                    ++p;
                }
            }
            i = j;
        }
    } else if (N > 0 && N < s->split_below) {
        const ms_u32 deficit = s->split_target > N ? s->split_target - N : 0u;
        const ms_u32 clones = 1u << s->split_levels;
        const double scale = 1.0 / (double)clones;
        ms_u32 added = 0;
        for (ms_u32 i = begin; i < end && added < deficit; ++i) {
            if (!(v.flags[i] & MS_FLAG_SPLIT)) continue;
            const double w = v.w[i];
            if (!(w >= s->split_min_weight)) continue;
            const double wm = w * scale;
            if (ms_exponent_field(wm) == 0) { ++rej; continue; }
            if (groups) {
                MsGroup g;
                g.W = w; g.start = i; g.end = i + 1; g.k = 1; g.kind = MS_KIND_SPLIT;
                g.cell = cell; g.species = species; g.bin = v.bin[i]; g.index = count;
                g.status = 0; g.n_out = 0; g.out_offset = 0; g.reloc_offset = 0; g.gflags = 0; g.pad = 0;
                groups[count] = g;
            }
            ++count;
            added += clones - 1;
        }
    }
    if (refused) *refused = rej;
    return count;
}

// ---------------------------------------------------------------- Philox stream layout merge-split-v1
//  key     k0 = seed bits 0..31 XOR "merg", k1 = seed bits 32..63 XOR "splt"
//  counter c0 = block (bits 0..7) | species (8..15) | purpose (16..23) | step bits 32..39 (24..31)
//          c1 = group index within its (cell, species) segment
//          c2 = cell id
//          c3 = step bits 0..31
//  domain  step < 2^40, species < 256, purpose 1..255, block < 256: the map is injective.
MS_FN pfa::Block ms_philox(const MsParams* P, ms_u32 cell, ms_u32 species, ms_u32 group, ms_u32 purpose, ms_u32 block) {
    pfa::Key2 k;
    k.v[0] = (pfa_u32)P->seed ^ MS_KEY_DOMAIN0;
    k.v[1] = (pfa_u32)(P->seed >> 32) ^ MS_KEY_DOMAIN1;
    pfa::Block c;
    c.v[0] = (block & 0xffu) | ((species & 0xffu) << 8) | ((purpose & 0xffu) << 16) |
             ((pfa_u32)((P->step >> 32) & 0xffu) << 24);
    c.v[1] = group;
    c.v[2] = cell;
    c.v[3] = (pfa_u32)P->step;
    return pfa::philox4x32_10(c, k);
}
// Uniform point of the open unit disk minus the origin by rejection from the
// square (-1,1)^2: trial b uses block b, x = 2 U(R(2b)) - 1, y = 2 U(R(2b+1)) - 1
// with U = pfa::u01_open (both steps exact). Returns 0 if all 256 trials were
// rejected (probability (1 - pi/4)^256 < 1e-170).
MS_FN int ms_disk_point(const MsParams* P, ms_u32 cell, ms_u32 species, ms_u32 group, ms_u32 purpose, double* x,
                        double* y, double* s) {
    for (ms_u32 b = 0; b < 256u; ++b) {
        const pfa::Block r = ms_philox(P, cell, species, group, purpose, b);
        const double px = 2 * pfa::u01_open(pfa::low_word(r)) - 1, py = 2 * pfa::u01_open(pfa::high_word(r)) - 1;
        const double ss = px * px + py * py;
        if (ss > 0 && ss < 1) { *x = px; *y = py; *s = ss; return 1; }
    }
    return 0;
}

// ---------------------------------------------------------------- merge kinematics
// k inputs (same species, cell and bin) -> 2 outputs of weight W/2 at the centroid c:
//   U = P/W (P = sum w u, proper velocity), kappa = E/W (E = sum w K),
//   |u_out|^2 = (kappa/m) (2 + kappa/(m c^2))   [K(u_out) = kappa exactly in real arithmetic],
//   q^2 = |u_out|^2 - |U|^2 (no real solution when negative: skipped),
//   u_a,b = U +/- q n, n perpendicular to U with azimuth from Philox:
//   (cos, sin) = ((x^2 - y^2)/s, 2xy/s) for a uniform disk point (x, y), s = x^2 + y^2.
//   U == 0: n isotropic (Marsaglia: (2x sqrt(1-s), 2y sqrt(1-s), 1-2s)).
// Real-arithmetic conservation: W/2 (u_a + u_b) = W U = P; W/2 (K_a + K_b) = W kappa = E.
// Rounding residuals are booked exactly by the ledger.
MS_FN ms_u32 ms_merge_group(const MsParams* P, MsView v, const MsGroup* g, MsMergeResult* r, ms_u32* gflags) {
    const MsSpecies* sp = &P->species[g->species];
    const double m = sp->mass_kg, W = g->W;
    ms_u32 flags = 0;
    // centroid (double-word sums of w x in canonical order)
    double ch[3] = {0, 0, 0}, cl[3] = {0, 0, 0};
    double ph[3] = {0, 0, 0}, pl[3] = {0, 0, 0}, eh = 0, el = 0;
    ms_u32 first = 0xffffffffu;
    for (ms_u32 i = g->start; i < g->end; ++i) {
        if (!(v.flags[i] & MS_FLAG_MERGE)) continue;
        if (first == 0xffffffffu) first = i;
        const double w = v.w[i];
        for (int d = 0; d < 3; ++d) {
            ms_dd_add_product(&ch[d], &cl[d], w, v.x[d][i]);
            ms_dd_add_product(&ph[d], &pl[d], w, v.u[d][i]);
        }
        ms_dd_add_product(&eh, &el, w, v.ke[i]);
    }
    double c[3];
    for (int d = 0; d < 3; ++d) c[d] = (ch[d] + cl[d]) / W;
    {
        ms_u32 cc = 0;
        int open = 0;
        const int inside = ms_cell_of(P, c[0], c[1], c[2], &cc, &open);
        if (!(inside && open && cc == g->cell)) {
            for (int d = 0; d < 3; ++d) c[d] = v.x[d][first];
            flags |= MS_GFLAG_CENTROID_FALLBACK;
        }
    }
    double U[3];
    for (int d = 0; d < 3; ++d) U[d] = (ph[d] + pl[d]) / W;
    const double kappa = (eh + el) / W;
    const double uout2 = (kappa / m) * (2 + kappa / sp->mc2_J);
    const double U2 = U[0] * U[0] + U[1] * U[1] + U[2] * U[2];
    const double q2 = uout2 - U2;
    *gflags = flags;
    if (!(ms_finite(q2) && q2 >= 0)) return MS_SKIP_NO_REAL_SOLUTION;
    const double q = ms_sqrt(q2);
    double n[3];
    double x, y, s;
    if (U2 > 0) {
        if (!ms_disk_point(P, g->cell, g->species, g->index, MS_PURPOSE_AZIMUTH, &x, &y, &s)) return MS_SKIP_RNG_EXHAUSTED;
        const double cs = (x * x - y * y) / s, sn = (2 * x * y) / s;
        int a = 0;
        double am = ms_abs(U[0]);
        if (ms_abs(U[1]) < am) { a = 1; am = ms_abs(U[1]); }
        if (ms_abs(U[2]) < am) { a = 2; }
        double e1[3];
        if (a == 0) { e1[0] = 0; e1[1] = U[2]; e1[2] = -U[1]; }
        else if (a == 1) { e1[0] = -U[2]; e1[1] = 0; e1[2] = U[0]; }
        else { e1[0] = U[1]; e1[1] = -U[0]; e1[2] = 0; }
        const double n1 = ms_sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
        for (int d = 0; d < 3; ++d) e1[d] = e1[d] / n1;
        double e2[3] = {U[1] * e1[2] - U[2] * e1[1], U[2] * e1[0] - U[0] * e1[2], U[0] * e1[1] - U[1] * e1[0]};
        const double n2 = ms_sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
        for (int d = 0; d < 3; ++d) e2[d] = e2[d] / n2;
        for (int d = 0; d < 3; ++d) n[d] = cs * e1[d] + sn * e2[d];
    } else {
        if (!ms_disk_point(P, g->cell, g->species, g->index, MS_PURPOSE_ISOTROPIC, &x, &y, &s)) return MS_SKIP_RNG_EXHAUSTED;
        const double t = 2 * ms_sqrt(1 - s);
        n[0] = x * t; n[1] = y * t; n[2] = 1 - 2 * s;
        flags |= MS_GFLAG_ISOTROPIC;
        *gflags = flags;
    }
    double ua[3], ub[3];
    for (int d = 0; d < 3; ++d) {
        const double qn = q * n[d];
        ua[d] = U[d] + qn;
        ub[d] = U[d] - qn;
    }
    // exactness guards of every ledger product (inputs and outputs)
    const double Wh = W * 0.5;
    const double ka = ms_kinetic_energy(m, ua[0], ua[1], ua[2]), kb = ms_kinetic_energy(m, ub[0], ub[1], ub[2]);
    int safe = ms_finite(ka) && ms_finite(kb) && ms_product_safe(Wh, ka) && ms_product_safe(Wh, kb);
    for (int d = 0; d < 3 && safe; ++d)
        safe = ms_finite(ua[d]) && ms_finite(ub[d]) && ms_product3_safe(Wh, ua[d], m) && ms_product3_safe(Wh, ub[d], m);
    for (ms_u32 i = g->start; i < g->end && safe; ++i) {
        if (!(v.flags[i] & MS_FLAG_MERGE)) continue;
        const double w = v.w[i];
        safe = ms_product_safe(w, v.ke[i]);
        for (int d = 0; d < 3 && safe; ++d) safe = ms_product3_safe(w, v.u[d][i], m);
    }
    if (!safe) return MS_SKIP_OPERAND_RANGE;
    for (int d = 0; d < 3; ++d) { r->c[d] = c[d]; r->ua[d] = ua[d]; r->ub[d] = ub[d]; }
    return MS_OK;
}
// Split exactness guard (formation already checked that w 2^-m is normal).
MS_FN ms_u32 ms_split_group(const MsParams* P, MsView v, const MsGroup* g) {
    const MsSpecies* sp = &P->species[g->species];
    const double m = sp->mass_kg, w = g->W, wm = w * (1.0 / (double)(1u << sp->split_levels));
    const ms_u32 i = g->start;
    int safe = ms_product_safe(w, v.ke[i]) && ms_product_safe(wm, v.ke[i]);
    for (int d = 0; d < 3 && safe; ++d) safe = ms_product3_safe(w, v.u[d][i], m) && ms_product3_safe(wm, v.u[d][i], m);
    return safe ? MS_OK : MS_SKIP_OPERAND_RANGE;
}
