// sw7_math.h -- deterministic elementary functions of surface owner v7.
//
// Compiled twice from this one file: by the host C++ compiler (-ffp-contract=off
// -fno-fast-math, SSE2 binary64) and by NVRTC (--fmad=false --prec-div=true
// --prec-sqrt=true --ftz=false). Only binary64 +, -, *, /, sqrt (each correctly
// rounded to nearest on both), exact integer <-> binary64 conversions and bit
// manipulation are used, in a fixed order, so both builds give the same bits
// (checked by tests/test_wall_gpu_v7.cpp). No libm/libdevice transcendental is
// called. Accuracy (measured against the host libm in tests/test_surface_v7.cpp):
// exp, log within 2 ulp; sincos_turn within 2e-16 absolute; lgamma within
// 1e-14 relative for a in [0.5, 200]; regularized incomplete gamma within 1e-13.
#pragma once

#if defined(__CUDACC_RTC__)
#define SW7_FN __device__ inline
typedef unsigned int sw7_u32;
typedef unsigned long long sw7_u64;
typedef long long sw7_i64;
SW7_FN sw7_u64 sw7_bits(double x) { return (sw7_u64)__double_as_longlong(x); }
SW7_FN double sw7_from_bits(sw7_u64 b) { return __longlong_as_double((long long)b); }
SW7_FN double sw7_sqrt(double x) { return sqrt(x); }
#else
#include <bit>
#include <cmath>
#include <cstdint>
#define SW7_FN inline
typedef std::uint32_t sw7_u32;
typedef std::uint64_t sw7_u64;
typedef std::int64_t sw7_i64;
SW7_FN sw7_u64 sw7_bits(double x) { return std::bit_cast<sw7_u64>(x); }
SW7_FN double sw7_from_bits(sw7_u64 b) { return std::bit_cast<double>(b); }
SW7_FN double sw7_sqrt(double x) { return std::sqrt(x); }
#endif

#define SW7_INF (sw7_from_bits(0x7FF0000000000000ull))
#define SW7_NAN (sw7_from_bits(0x7FF8000000000000ull))

namespace sw7 {

SW7_FN bool finite(double x) { return (sw7_bits(x) & 0x7FF0000000000000ull) != 0x7FF0000000000000ull; }
SW7_FN double fabs_(double x) { return sw7_from_bits(sw7_bits(x) & 0x7FFFFFFFFFFFFFFFull); }
SW7_FN double fmin_(double a, double b) { return a < b ? a : b; }
SW7_FN double fmax_(double a, double b) { return a > b ? a : b; }
// 2^k for integer k in [-1074, 1023] (exact).
SW7_FN double pow2i(int k) {
    if (k >= -1022) return sw7_from_bits((sw7_u64)(k + 1023) << 52);
    return sw7_from_bits(1ull << (k + 1074));
}
// Round to nearest integer (ties away from zero), |x| < 2^52.
SW7_FN double round_(double x) {
    const double a = fabs_(x) + 0.5;
    const double r = (double)(sw7_i64)a;   // truncation (exact)
    return x < 0 ? -r : r;
}
SW7_FN double floor_(double x) {
    const double t = (double)(sw7_i64)x;   // truncation toward zero, |x| < 2^62
    return (t > x) ? t - 1.0 : t;
}

// exp(x): Cody-Waite reduction x = k ln2 + r, |r| <= ln2/2, degree-13 Taylor polynomial.
SW7_FN double exp_(double x) {
    if (!(x == x)) return x;
    if (x > 709.782712893384) return SW7_INF;
    if (x < -745.2) return 0.0;
    const double inv_ln2 = 1.4426950408889634074;
    const double ln2_hi = 6.93147180369123816490e-01;   // 0x3FE62E42FEE00000: trailing zeros, k*ln2_hi exact
    const double ln2_lo = 1.90821492927058770002e-10;
    const double k = round_(x * inv_ln2);
    const double r = (x - k * ln2_hi) - k * ln2_lo;
    double p = 1.0 / 6227020800.0;                      // 1/13!
    p = p * r + 1.0 / 479001600.0;
    p = p * r + 1.0 / 39916800.0;
    p = p * r + 1.0 / 3628800.0;
    p = p * r + 1.0 / 362880.0;
    p = p * r + 1.0 / 40320.0;
    p = p * r + 1.0 / 5040.0;
    p = p * r + 1.0 / 720.0;
    p = p * r + 1.0 / 120.0;
    p = p * r + 1.0 / 24.0;
    p = p * r + 1.0 / 6.0;
    p = p * r + 0.5;
    p = p * r + 1.0;
    p = p * r + 1.0;
    const int ki = (int)k;
    if (ki >= -1021 && ki <= 1023) return p * pow2i(ki);
    if (ki > 1023) return (p * pow2i(1023)) * pow2i(ki - 1023);
    return (p * pow2i(ki + 600)) * pow2i(-600);
}

// log(x), x > 0: x = m 2^e, m in [sqrt(1/2), sqrt(2)); log m = 2 atanh((m-1)/(m+1)).
SW7_FN double log_(double x) {
    if (!(x > 0)) return x == 0 ? -SW7_INF : SW7_NAN;
    if (!finite(x)) return x;
    int e = 0;
    sw7_u64 b = sw7_bits(x);
    if ((b >> 52) == 0) {   // subnormal: scale by 2^54
        x = x * pow2i(54);
        b = sw7_bits(x);
        e = -54;
    }
    e += (int)((b >> 52) & 0x7FF) - 1023;
    double m = sw7_from_bits((b & 0x000FFFFFFFFFFFFFull) | 0x3FF0000000000000ull);   // [1, 2)
    if (m > 1.4142135623730951) { m = m * 0.5; e += 1; }
    const double f = (m - 1.0) / (m + 1.0);
    const double s = f * f;
    double p = 1.0 / 25.0;
    p = p * s + 1.0 / 23.0;
    p = p * s + 1.0 / 21.0;
    p = p * s + 1.0 / 19.0;
    p = p * s + 1.0 / 17.0;
    p = p * s + 1.0 / 15.0;
    p = p * s + 1.0 / 13.0;
    p = p * s + 1.0 / 11.0;
    p = p * s + 1.0 / 9.0;
    p = p * s + 1.0 / 7.0;
    p = p * s + 1.0 / 5.0;
    p = p * s + 1.0 / 3.0;
    const double logm = 2.0 * f + 2.0 * f * (s * p);
    const double ln2_hi = 6.93147180369123816490e-01, ln2_lo = 1.90821492927058770002e-10;
    const double de = (double)e;
    return de * ln2_hi + (de * ln2_lo + logm);
}
SW7_FN double log1p_small(double x) { return log_(1.0 + x); }   // used only where 1 + x is exact enough
SW7_FN double pow_(double x, double y) {
    if (y == 0) return 1.0;
    if (x == 0) return y > 0 ? 0.0 : SW7_INF;
    if (x == 1) return 1.0;
    return exp_(y * log_(x));
}
SW7_FN double cbrt_cube(double x) { return x * x * x; }

// sin and cos of 2*pi*t for t in [0, 1): reduction in turns (exact), Taylor on [0, pi/4].
SW7_FN void sincos_turn(double t, double* s_out, double* c_out) {
    const double q = floor_(t * 4.0);           // quadrant (t*4 exact)
    double r = t - q * 0.25;                    // [0, 1/4), exact
    bool swap = false;
    if (r > 0.125) { r = 0.25 - r; swap = true; }   // exact
    const double x = r * 6.283185307179586477;
    const double x2 = x * x;
    double sp = -1.0 / 355687428096000.0;       // -1/17!
    sp = sp * x2 + 1.0 / 1307674368000.0;
    sp = sp * x2 - 1.0 / 6227020800.0;
    sp = sp * x2 + 1.0 / 39916800.0;
    sp = sp * x2 - 1.0 / 362880.0;
    sp = sp * x2 + 1.0 / 5040.0;
    sp = sp * x2 - 1.0 / 120.0;
    sp = sp * x2 + 1.0 / 6.0;
    const double sn = x - x * (x2 * sp);
    double cp = 1.0 / 20922789888000.0;         // 1/16!
    cp = cp * x2 - 1.0 / 87178291200.0;
    cp = cp * x2 + 1.0 / 479001600.0;
    cp = cp * x2 - 1.0 / 3628800.0;
    cp = cp * x2 + 1.0 / 40320.0;
    cp = cp * x2 - 1.0 / 720.0;
    cp = cp * x2 + 1.0 / 24.0;
    cp = cp * x2 - 0.5;
    const double cs = 1.0 + x2 * cp;
    double s = swap ? cs : sn, c = swap ? sn : cs;
    const int qi = (int)q & 3;
    if (qi == 0) { *s_out = s; *c_out = c; }
    else if (qi == 1) { *s_out = c; *c_out = -s; }
    else if (qi == 2) { *s_out = -s; *c_out = -c; }
    else { *s_out = -c; *c_out = s; }
}

// lgamma(a), a > 0: Lanczos (g = 7, n = 9), recurrence below 0.5.
SW7_FN double lgamma_(double a) {
    double shift = 0;
    while (a < 0.5) { shift -= log_(a); a += 1.0; }
    const double g = 7.0;
    const double c0 = 0.99999999999980993, c1 = 676.5203681218851, c2 = -1259.1392167224028, c3 = 771.32342877765313,
                 c4 = -176.61502916214059, c5 = 12.507343278686905, c6 = -0.13857109526572012, c7 = 9.9843695780195716e-6,
                 c8 = 1.5056327351493116e-7;
    const double z = a - 1.0;
    double x = c0;
    x = x + c1 / (z + 1.0);
    x = x + c2 / (z + 2.0);
    x = x + c3 / (z + 3.0);
    x = x + c4 / (z + 4.0);
    x = x + c5 / (z + 5.0);
    x = x + c6 / (z + 6.0);
    x = x + c7 / (z + 7.0);
    x = x + c8 / (z + 8.0);
    const double t = z + g + 0.5;
    return shift + (0.91893853320467274178 + (z + 0.5) * log_(t) - t + log_(x));
}

// Regularized lower incomplete gamma P(a, x), a > 0, x >= 0 (series / Lentz continued fraction).
SW7_FN double gammainc_P(double a, double x) {
    if (x <= 0) return 0.0;
    const double lfront = -x + a * log_(x) - lgamma_(a);
    if (x < a + 1.0) {
        double ap = a, del = 1.0 / a, sum = del;
        for (int n = 0; n < 500; ++n) {
            ap += 1.0;
            del = del * x / ap;
            sum += del;
            if (fabs_(del) < fabs_(sum) * 1e-17) break;
        }
        const double v = sum * exp_(lfront);
        return v > 1.0 ? 1.0 : v;
    }
    const double tiny = 1e-300;
    double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
    for (int i = 1; i < 500; ++i) {
        const double an = -(double)i * ((double)i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs_(d) < tiny) d = tiny;
        c = b + an / c;
        if (fabs_(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (fabs_(del - 1.0) < 1e-17) break;
    }
    const double q = exp_(lfront) * h;
    const double v = 1.0 - q;
    return v < 0 ? 0.0 : v;
}
// Inverse of P(a, .) on [0, xmax] by 64 bisections of the target P(a, x) = p (deterministic).
SW7_FN double gammainc_P_inverse(double a, double p, double xmax) {
    double lo = 0, hi = xmax;
    for (int i = 0; i < 64; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (mid <= lo || mid >= hi) break;
        if (gammainc_P(a, mid) < p) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}
// erf(x) for x >= 0 from the incomplete gamma: erf(x) = P(1/2, x^2).
SW7_FN double erf_pos(double x) { return x <= 0 ? 0.0 : gammainc_P(0.5, x * x); }

}  // namespace sw7
