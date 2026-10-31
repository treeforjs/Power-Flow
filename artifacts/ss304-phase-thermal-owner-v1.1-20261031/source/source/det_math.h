// det_math.h -- arithmetic shared bit for bit by the host reference and the NVRTC kernels.
//
//  * pq_i128: signed 128-bit two's-complement fixed-point word (lo, hi). Energies are held in units of
//    2^-96 J and atom inventories in units of 2^-32 atoms; addition/subtraction is exact and associative,
//    so every ledger sum is exact and independent of summation order.
//  * pq_from_double / pq_to_double: deterministic conversions (integer bit manipulation; round to
//    nearest, ties to even, from a double; two-word conversion to a double).
//  * pq_exp / pq_log: ports of the fdlibm e_exp.c / e_log.c algorithms using only IEEE +, -, *, / and
//    bit operations. Compiled with -ffp-contract=off (host) and --fmad=false (device) they give the same
//    bits on both sides (libm and CUDA exp/log are not used anywhere in the per-cell arithmetic).
// No headers are included so that the file compiles unchanged under NVRTC.
#pragma once

#if defined(__CUDACC_RTC__) || defined(__CUDACC__)
#define PQ_HD __device__ __forceinline__
#define PQ_DEVICE 1
#else
#define PQ_HD inline
#define PQ_DEVICE 0
#endif

typedef unsigned long long pq_u64;
typedef long long pq_i64;
typedef unsigned int pq_u32;

PQ_HD pq_u64 pq_bits(double x) {
#if PQ_DEVICE
    return (pq_u64)__double_as_longlong(x);
#else
    pq_u64 u; __builtin_memcpy(&u, &x, 8); return u;
#endif
}
PQ_HD double pq_from_bits(pq_u64 u) {
#if PQ_DEVICE
    return __longlong_as_double((pq_i64)u);
#else
    double x; __builtin_memcpy(&x, &u, 8); return x;
#endif
}
PQ_HD bool pq_finite(double x) { return ((pq_bits(x) >> 52) & 0x7ffULL) != 0x7ffULL; }
PQ_HD double pq_fabs(double x) { return pq_from_bits(pq_bits(x) & 0x7fffffffffffffffULL); }
PQ_HD double pq_min(double a, double b) { return b < a ? b : a; }
PQ_HD double pq_max(double a, double b) { return b > a ? b : a; }
// Correctly rounded square root on both sides (IEEE sqrt; --prec-sqrt=true on the device).
PQ_HD double pq_sqrt(double x) {
#if PQ_DEVICE
    return sqrt(x);
#else
    return __builtin_sqrt(x);
#endif
}
PQ_HD double pq_floor(double x) {
#if PQ_DEVICE
    return floor(x);
#else
    return __builtin_floor(x);
#endif
}
// Next representable double toward +inf / -inf (bit stepping; finite nonzero-safe).
PQ_HD double pq_next_up(double x) {
    if (!pq_finite(x)) return x;
    if (x == 0.0) return pq_from_bits(1ULL);
    pq_u64 u = pq_bits(x);
    return pq_from_bits((u >> 63) ? u - 1 : u + 1);
}
PQ_HD double pq_next_down(double x) { return -pq_next_up(-x); }

// ---------------------------------------------------------------------------------------------------------
// 128-bit fixed-point words.
struct pq_i128 { pq_u64 lo; pq_i64 hi; };
PQ_HD pq_i128 pq_zero() { pq_i128 r; r.lo = 0; r.hi = 0; return r; }
PQ_HD bool pq_is_zero(pq_i128 a) { return a.lo == 0 && a.hi == 0; }
PQ_HD bool pq_is_neg(pq_i128 a) { return a.hi < 0; }
PQ_HD bool pq_eq(pq_i128 a, pq_i128 b) { return a.lo == b.lo && a.hi == b.hi; }
PQ_HD int pq_cmp(pq_i128 a, pq_i128 b) {
    if (a.hi != b.hi) return a.hi < b.hi ? -1 : 1;
    if (a.lo != b.lo) return a.lo < b.lo ? -1 : 1;
    return 0;
}
// Wrapping add; *overflow is set (never cleared) on signed overflow.
PQ_HD pq_i128 pq_add(pq_i128 a, pq_i128 b, int* overflow) {
    pq_i128 r;
    r.lo = a.lo + b.lo;
    const pq_u64 carry = r.lo < a.lo ? 1ULL : 0ULL;
    r.hi = (pq_i64)((pq_u64)a.hi + (pq_u64)b.hi + carry);
    if (((a.hi < 0) == (b.hi < 0)) && ((r.hi < 0) != (a.hi < 0))) *overflow = 1;
    return r;
}
PQ_HD pq_i128 pq_neg(pq_i128 a, int* overflow) {
    pq_i128 r;
    r.lo = ~a.lo + 1ULL;
    r.hi = (pq_i64)(~(pq_u64)a.hi + (r.lo == 0 ? 1ULL : 0ULL));
    if (a.hi == (pq_i64)0x8000000000000000ULL && a.lo == 0) *overflow = 1;
    return r;
}
PQ_HD pq_i128 pq_sub(pq_i128 a, pq_i128 b, int* overflow) { return pq_add(a, pq_neg(b, overflow), overflow); }

// Unsigned 128-bit magnitude helpers.
struct pq_u128 { pq_u64 lo, hi; };
PQ_HD pq_u128 pq_mag(pq_i128 a, int* overflow) {
    if (a.hi < 0) a = pq_neg(a, overflow);
    pq_u128 m; m.lo = a.lo; m.hi = (pq_u64)a.hi; return m;
}
// Floor division by a positive 32-bit divisor with remainder in [0, d): a = q*d + r.
PQ_HD pq_i128 pq_divmod_u32(pq_i128 a, pq_u32 d, pq_u32* rem, int* overflow) {
    const bool neg = a.hi < 0;
    pq_u128 m = pq_mag(a, overflow);
    // Long division in 32-bit limbs (most significant first).
    pq_u32 limb[4] = {(pq_u32)(m.hi >> 32), (pq_u32)m.hi, (pq_u32)(m.lo >> 32), (pq_u32)m.lo};
    pq_u64 r = 0;
    for (int i = 0; i < 4; ++i) {
        const pq_u64 cur = (r << 32) | limb[i];
        limb[i] = (pq_u32)(cur / d);
        r = cur % d;
    }
    pq_i128 q; q.hi = (pq_i64)(((pq_u64)limb[0] << 32) | limb[1]); q.lo = ((pq_u64)limb[2] << 32) | limb[3];
    if (!neg) { *rem = (pq_u32)r; return q; }
    // -(m) = -(q d + r) = (-q-1) d + (d - r) when r > 0.
    q = pq_neg(q, overflow);
    if (r == 0) { *rem = 0; return q; }
    pq_i128 one; one.lo = 1; one.hi = 0;
    q = pq_sub(q, one, overflow);
    *rem = (pq_u32)(d - r);
    return q;
}
// Round-to-nearest-even of x * 2^scale_pow2 as a 128-bit integer; *bad set if x is not finite or the
// magnitude reaches 2^126.
PQ_HD pq_i128 pq_from_double(double x, int scale_pow2, int* bad) {
    pq_i128 r = pq_zero();
    const pq_u64 u = pq_bits(x);
    const int sign = (int)(u >> 63);
    const int e = (int)((u >> 52) & 0x7ffULL);
    pq_u64 m = u & 0x000fffffffffffffULL;
    if (e == 0x7ff) { *bad = 1; return r; }
    int s;
    if (e == 0) s = -1074 + scale_pow2; else { m |= 0x0010000000000000ULL; s = e - 1075 + scale_pow2; }
    if (m == 0) return r;
    pq_u128 v; v.lo = 0; v.hi = 0;
    if (s >= 0) {
        if (s + 53 > 126) { *bad = 1; return r; }
        if (s >= 64) { v.hi = m << (s - 64); v.lo = 0; }
        else if (s == 0) { v.lo = m; v.hi = 0; }
        else { v.lo = m << s; v.hi = m >> (64 - s); }
    } else {
        const int sh = -s;
        if (sh >= 64) { v.lo = 0; }   // |value| < 2^-10 (m < 2^53): rounds to zero
        else {
            pq_u64 q = m >> sh;
            const pq_u64 rem = m & ((1ULL << sh) - 1ULL);
            const pq_u64 half = 1ULL << (sh - 1);
            if (rem > half || (rem == half && (q & 1ULL))) ++q;
            v.lo = q;
        }
    }
    r.lo = v.lo; r.hi = (pq_i64)v.hi;
    if (sign) { int ov = 0; r = pq_neg(r, &ov); if (ov) *bad = 1; }
    return r;
}
// Deterministic two-word conversion of a * 2^-scale_pow2 (not necessarily correctly rounded).
PQ_HD double pq_to_double(pq_i128 a, double inv_scale) {
    int ov = 0;
    const bool neg = a.hi < 0;
    pq_u128 m = pq_mag(a, &ov);
    const double d = ((double)m.hi * 18446744073709551616.0 + (double)m.lo) * inv_scale;
    return neg ? -d : d;
}

// Units of the fixed-point state.
#define PQ_E_POW2 96          // energy word = J * 2^96
#define PQ_N_POW2 32          // atom word = atoms * 2^32
#define PQ_E_INV 0x1p-96
#define PQ_N_INV 0x1p-32

// ---------------------------------------------------------------------------------------------------------
// fdlibm e_exp.c (Sun Microsystems, 1993; freely redistributable) ported to bit helpers.
PQ_HD double pq_exp(double x) {
    const double o_threshold = 7.09782712893383973096e+02, u_threshold = -7.45133219101941108420e+02;
    const double ln2HI = 6.93147180369123816490e-01, ln2LO = 1.90821492927058770002e-10;
    const double invln2 = 1.44269504088896338700e+00;
    const double P1 = 1.66666666666666019037e-01, P2 = -2.77777777770155933842e-03,
                 P3 = 6.61375632143793436117e-05, P4 = -1.65339022054652515390e-06,
                 P5 = 4.13813679705723846039e-08;
    if (x != x) return x;
    if (x > o_threshold) return pq_from_bits(0x7ff0000000000000ULL);
    if (x < u_threshold) return 0.0;
    const pq_u64 ux = pq_bits(x);
    const pq_u32 hx = (pq_u32)((ux >> 32) & 0x7fffffffU);
    double hi = x, lo = 0.0, r = x;
    int k = 0;
    if (hx > 0x3fd62e42U) {
        const double half = x < 0 ? -0.5 : 0.5;
        k = (int)(invln2 * x + half);
        const double t = (double)k;
        hi = x - t * ln2HI;
        lo = t * ln2LO;
        r = hi - lo;
    } else if (hx < 0x3e300000U) {
        return 1.0 + x;
    }
    const double t = r * r;
    const double c = r - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
    double y;
    if (k == 0) return 1.0 - ((r * c) / (c - 2.0) - r);
    y = 1.0 - ((lo - (r * c) / (2.0 - c)) - hi);
    pq_u64 uy = pq_bits(y);
    // exponent adjustment in unsigned arithmetic (two's-complement wrap of a negative k is intended)
    if (k >= -1021) {
        uy = uy + ((pq_u64)(pq_i64)k << 52);
        return pq_from_bits(uy);
    }
    uy = uy + ((pq_u64)(pq_i64)(k + 1000) << 52);
    return pq_from_bits(uy) * 9.33263618503218878990e-302;
}

// fdlibm e_log.c ported to bit helpers; defined for finite x > 0 (returns NaN otherwise).
PQ_HD double pq_log(double x) {
    const double ln2_hi = 6.93147180369123816490e-01, ln2_lo = 1.90821492927058770002e-10,
                 two54 = 1.80143985094819840000e+16;
    const double Lg1 = 6.666666666666735130e-01, Lg2 = 3.999999999940941908e-01, Lg3 = 2.857142874366239149e-01,
                 Lg4 = 2.222219843214978396e-01, Lg5 = 1.818357216161805012e-01, Lg6 = 1.531383769920937332e-01,
                 Lg7 = 1.479819860511658591e-01;
    if (!(x > 0.0) || !pq_finite(x)) return pq_from_bits(0x7ff8000000000000ULL);
    pq_u64 ux = pq_bits(x);
    int hx = (int)(ux >> 32);
    int k = 0;
    if (hx < 0x00100000) { k -= 54; x *= two54; ux = pq_bits(x); hx = (int)(ux >> 32); }
    k += (hx >> 20) - 1023;
    hx &= 0x000fffff;
    int i = (hx + 0x95f64) & 0x100000;
    ux = ((pq_u64)(pq_u32)(hx | (i ^ 0x3ff00000)) << 32) | (ux & 0xffffffffULL);
    x = pq_from_bits(ux);
    k += (i >> 20);
    const double f = x - 1.0;
    double dk;
    if ((0x000fffff & (2 + hx)) < 3) {
        if (f == 0.0) { if (k == 0) return 0.0; dk = (double)k; return dk * ln2_hi + dk * ln2_lo; }
        const double R = f * f * (0.5 - 0.33333333333333333 * f);
        if (k == 0) return f - R;
        dk = (double)k;
        return dk * ln2_hi - ((R - dk * ln2_lo) - f);
    }
    const double s = f / (2.0 + f);
    dk = (double)k;
    const double z = s * s;
    i = hx - 0x6147a;
    const double w = z * z;
    const int j = 0x6b851 - hx;
    const double t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
    const double t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
    i |= j;
    const double R = t2 + t1;
    if (i > 0) {
        const double hfsq = 0.5 * f * f;
        if (k == 0) return f - (hfsq - s * (hfsq + R));
        return dk * ln2_hi - ((hfsq - (s * (hfsq + R) + dk * ln2_lo)) - f);
    }
    if (k == 0) return f - s * (f - R);
    return dk * ln2_hi - ((s * (f - R) - dk * ln2_lo) - f);
}
