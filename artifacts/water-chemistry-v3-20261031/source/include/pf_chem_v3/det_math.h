// det_math.h -- deterministic elementary functions for chemistry v3.
//
// Compiled from these bytes by the host C++ compiler (-ffp-contract=off
// -fno-fast-math, IEEE binary64, round to nearest) and by NVRTC (--fmad=false
// --prec-div=true --prec-sqrt=true --ftz=false). Only +, -, *, /, sqrt
// (correctly rounded on both), comparisons, integer arithmetic and exact
// bit manipulation are used, in one fixed operation order, so host and device
// produce identical bits. No libm or libdevice transcendental is called.
//
//   pfc_log(x)          natural logarithm, x > 0 finite. Algorithm of the
//                       FreeBSD/fdlibm e_log.c (argument reduction to
//                       [sqrt(2)/2, sqrt(2)], s = f/(2+f), degree-14 series in
//                       s with the published Lg1..Lg7 coefficients); error
//                       below 1 ulp (checked against the host libm in
//                       tests/test_core.cpp).
//   pfc_log1p(x)        log(1+x), x > -1, via u = 1+x and log(u) * x/(u-1)
//                       (exact compensation of the rounding of u); a few ulp.
//   pfc_sincos_2pi(u)   cos(2 pi u), sin(2 pi u) for u in [0, 1): quadrant
//                       reduction 4u = k + r (exact), reflection at r = 1/2
//                       (exact), Taylor polynomials on [0, pi/4]; a few ulp.
//   pfc_kinetic_eV(m, p2)  K = p2 / (sqrt(p2 + m^2) + m), the kinetic energy
//                       of rest energy m and squared momentum p2 (eV units,
//                       c = 1); no cancellation.
//
// Permission notice of the fdlibm source whose algorithm and coefficients
// pfc_log follows: "Copyright (C) 1993 by Sun Microsystems, Inc. All rights
// reserved. Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this software is freely
// granted, provided that this notice is preserved."
#pragma once

#if defined(__CUDACC_RTC__)
#define PFC_FN __device__ inline
typedef unsigned int pfc_u32;
typedef unsigned long long pfc_u64;
typedef long long pfc_i64;
PFC_FN pfc_u64 pfc_bits(double x) { return (pfc_u64)__double_as_longlong(x); }
PFC_FN double pfc_from_bits(pfc_u64 b) { return __longlong_as_double((long long)b); }
PFC_FN double pfc_sqrt(double x) { return __dsqrt_rn(x); }
#else
#include <bit>
#include <cmath>
#include <cstdint>
#define PFC_FN inline
typedef std::uint32_t pfc_u32;
typedef std::uint64_t pfc_u64;
typedef std::int64_t pfc_i64;
PFC_FN pfc_u64 pfc_bits(double x) { return std::bit_cast<pfc_u64>(x); }
PFC_FN double pfc_from_bits(pfc_u64 b) { return std::bit_cast<double>(b); }
PFC_FN double pfc_sqrt(double x) { return std::sqrt(x); }   // IEEE correctly rounded
#endif

namespace pfc {

// fdlibm e_log.c constants, written as exact bit patterns.
#define PFC_LN2_HI 0x3FE62E42FEE00000ull   // 6.93147180369123816490e-01
#define PFC_LN2_LO 0x3DEA39EF35793C76ull   // 1.90821492927058770002e-10
#define PFC_TWO54 0x4350000000000000ull    // 2^54
#define PFC_LG1 0x3FE5555555555593ull
#define PFC_LG2 0x3FD999999997FA04ull
#define PFC_LG3 0x3FD2492494229359ull
#define PFC_LG4 0x3FCC71C51D8E78AFull
#define PFC_LG5 0x3FC7466496CB03DEull
#define PFC_LG6 0x3FC39A09D078C69Full
#define PFC_LG7 0x3FC2F112DF3E5244ull

// Natural logarithm of a positive finite x. Returns -inf for 0 and NaN for
// negative or NaN input (callers never pass those; checked by the host).
PFC_FN double log(double x) {
    const double ln2_hi = pfc_from_bits(PFC_LN2_HI), ln2_lo = pfc_from_bits(PFC_LN2_LO);
    const double Lg1 = pfc_from_bits(PFC_LG1), Lg2 = pfc_from_bits(PFC_LG2), Lg3 = pfc_from_bits(PFC_LG3),
                 Lg4 = pfc_from_bits(PFC_LG4), Lg5 = pfc_from_bits(PFC_LG5), Lg6 = pfc_from_bits(PFC_LG6),
                 Lg7 = pfc_from_bits(PFC_LG7);
    pfc_u64 b = pfc_bits(x);
    int hx = (int)(pfc_i64)(b >> 32);   // signed high word
    const pfc_u32 lx = (pfc_u32)b;
    int k = 0;
    if (hx < 0x00100000) {   // x < 2^-1022
        if ((((pfc_u32)hx & 0x7fffffffu) | lx) == 0u) return -pfc_from_bits(0x7FF0000000000000ull);
        if (hx < 0) return pfc_from_bits(0x7FF8000000000000ull);
        k -= 54;
        x = x * pfc_from_bits(PFC_TWO54);
        b = pfc_bits(x);
        hx = (int)(pfc_i64)(b >> 32);
    }
    if (hx >= 0x7ff00000) return x + x;
    k += (hx >> 20) - 1023;
    hx &= 0x000fffff;
    const int i0 = (hx + 0x95f64) & 0x100000;
    // normalize x or x/2
    b = pfc_bits(x);
    b = ((pfc_u64)(pfc_u32)(hx | (i0 ^ 0x3ff00000)) << 32) | (b & 0xFFFFFFFFull);
    x = pfc_from_bits(b);
    k += (i0 >> 20);
    const double f = x - 1.0;
    if ((0x000fffff & (2 + hx)) < 3) {   // -2^-20 <= f < 2^-20
        if (f == 0.0) {
            if (k == 0) return 0.0;
            const double dk = (double)k;
            return dk * ln2_hi + dk * ln2_lo;
        }
        const double R = f * f * (0.5 - 0.33333333333333333 * f);
        if (k == 0) return f - R;
        const double dk = (double)k;
        return dk * ln2_hi - ((R - dk * ln2_lo) - f);
    }
    const double s = f / (2.0 + f);
    const double dk = (double)k;
    const double z = s * s;
    int i = hx - 0x6147a;
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

// log(1 + x) for x > -1. u = 1 + x is rounded; log(u) * (x / (u - 1)) corrects
// the rounding of u to first order (u - 1 is exact by Sterbenz near 0).
PFC_FN double log1p(double x) {
    const double u = 1.0 + x;
    if (u == 1.0) return x;
    return pfc::log(u) * (x / (u - 1.0));
}

// Taylor polynomials of sin and cos on [0, pi/4] in Horner form. Truncation
// below 1e-20 relative on that interval (next terms a^21/21!, a^20/20!).
// cos a = 1 - z q(z), q(z) = 1/2! - z/4! + ... + z^8/18!.
// sin a = a - a z q(z), q(z) = 1/3! - z/5! + z^2/7! - ... - z^7/17! + z^8/19!.
PFC_FN double sin_poly(double a) {
    const double z = a * a;
    double q = 1.0 / 121645100408832000.0;            // 1/19!
    q = 1.0 / 355687428096000.0 - z * q;              // 1/17!
    q = 1.0 / 1307674368000.0 - z * q;                // 1/15!
    q = 1.0 / 6227020800.0 - z * q;                   // 1/13!
    q = 1.0 / 39916800.0 - z * q;                     // 1/11!
    q = 1.0 / 362880.0 - z * q;                       // 1/9!
    q = 1.0 / 5040.0 - z * q;                         // 1/7!
    q = 1.0 / 120.0 - z * q;                          // 1/5!
    q = 1.0 / 6.0 - z * q;                            // 1/3!
    return a - a * (z * q);
}
PFC_FN double cos_poly(double a) {
    const double z = a * a;
    double q = 1.0 / 6402373705728000.0;              // 1/18!
    q = 1.0 / 20922789888000.0 - z * q;               // 1/16!
    q = 1.0 / 87178291200.0 - z * q;                  // 1/14!
    q = 1.0 / 479001600.0 - z * q;                    // 1/12!
    q = 1.0 / 3628800.0 - z * q;                      // 1/10!
    q = 1.0 / 40320.0 - z * q;                        // 1/8!
    q = 1.0 / 720.0 - z * q;                          // 1/6!
    q = 1.0 / 24.0 - z * q;                           // 1/4!
    q = 0.5 - z * q;                                  // 1/2!
    return 1.0 - z * q;
}

// cos(2 pi u), sin(2 pi u) for u in [0, 1).
PFC_FN void sincos_2pi(double u, double* c, double* s) {
    const double half_pi = 1.5707963267948966192313216916397514;   // rounded pi/2
    const double y = 4.0 * u;                 // exact
    const int k = (int)y;                      // quadrant 0..3
    const double r = y - (double)k;            // exact, [0, 1)
    double cr, sr;
    if (r <= 0.5) {
        const double a = half_pi * r;
        cr = cos_poly(a);
        sr = sin_poly(a);
    } else {
        const double a = half_pi * (1.0 - r);  // 1 - r exact (Sterbenz)
        cr = sin_poly(a);
        sr = cos_poly(a);
    }
    switch (k & 3) {
        case 0: *c = cr; *s = sr; break;
        case 1: *c = -sr; *s = cr; break;
        case 2: *c = -cr; *s = -sr; break;
        default: *c = sr; *s = -cr; break;
    }
}

// Kinetic energy from rest energy m and squared momentum p2 (same units).
PFC_FN double kinetic_from_p2(double m, double p2) {
    if (p2 == 0.0) return 0.0;
    return p2 / (pfc_sqrt(p2 + m * m) + m);
}
// gamma - 1 from squared proper velocity u2 (u = gamma beta), no cancellation.
PFC_FN double gamma_minus_one(double u2) {
    return u2 / (pfc_sqrt(1.0 + u2) + 1.0);
}

}  // namespace pfc
