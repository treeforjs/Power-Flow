#pragma once
#include <cfloat>
#include <cmath>

// Exact value of fma(a,b,-x), where x is the rounded product. The shortcut
// uses Dekker's split only for x87's 64-significand-bit extended format,
// extended precision / nearest rounding, and a conservative normal range.
// Other formats, precision/rounding modes and exceptional values retain fma.
// This is value-equivalent, not an exception-flag emulation: splitting may
// set FE_INEXACT even when a*b is exact. No caller uses floating-point flags.
namespace pf_exact {
inline bool product_error_environment() {
#if LDBL_MANT_DIG == 64 && LDBL_MAX_EXP == 16384 && FLT_RADIX == 2 && (defined(__x86_64__) || defined(__i386__))
    unsigned short control;
    __asm__ volatile ("fnstcw %0" : "=m" (control));
    // Precision=11 (64 bits), rounding=00 (nearest), exceptions masked.
    return (control & 0x0f3fu) == 0x033fu;
#else
    return false;
#endif
}
inline bool product_error_fast_domain(long double a,long double b,long double x) {
#if LDBL_MANT_DIG == 64 && LDBL_MAX_EXP == 16384 && FLT_RADIX == 2
    if (!product_error_environment()) return false;
    const long double fa=a<0?-a:a,fb=b<0?-b:b;
    constexpr long double lo=0x1p-4000L,hi=0x1p4000L;
    return fa>=lo&&fa<=hi&&fb>=lo&&fb<=hi&&a*b==x;
#else
    (void)a;(void)b;(void)x;
    return false;
#endif
}
inline long double product_error(long double a,long double b,long double x) {
    if (!product_error_fast_domain(a,b,x)) return std::fma(a,b,-x);
    constexpr long double split=0x1p32L+1.0L;
    long double t=split*a;
    const long double ah=t-(t-a),al=a-ah;
    t=split*b;
    const long double bh=t-(t-b),bl=b-bh;
    return ((ah*bh-x)+ah*bl+al*bh)+al*bl;
}
}
