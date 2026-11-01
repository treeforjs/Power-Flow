// x87_v34.hpp -- exact emulation of the x87 80-bit extended format (64-bit significand, round to nearest
// even, the precision of `long double` on x86-64 Linux) for the long-double arithmetic of the sealed host
// code that Hall v3.4 moves to the device: the deposit support of the endpoint continuity guard
// (census_charge.hpp DepositSupportBuilder) and the free-charge comparison of the census admission
// (same_net_charge). Compiled by NVRTC (device) and by the host compiler (CPU emulation used in the
// qualification against the real x87 arithmetic).
//
// Values are (-1)^s * m * 2^e with a 64-bit significand m whose top bit is set, or zero (m == 0). The
// exponent is an unbounded int: every value of the emulated expressions lies far inside the x87 exponent
// range (|e| < 1200), so x87 overflow, underflow and denormals never occur and are not modelled.
// Operations are correctly rounded (nearest, ties to even) to 64 significand bits, as one x87 FADD/FSUB/
// FMUL with precision control "extended":
//   add/sub: exact 192-bit window (an operand more than 65 binary orders below the other cannot change
//            the rounded result; it is dropped), normalisation, one rounding;
//   mul:     exact 128-bit product, one rounding.
// Conversions: from binary64 (exact, also subnormal), from small integers (exact), to binary64 with
// nearest-even rounding and optionally the sealed "if (long double)double(x) < x then nextafter(+inf)"
// step (rounding toward +inf for positive values), floor to integer (exact).
#pragma once

#if defined(__CUDACC_RTC__)
#define X87_FN __device__ inline
typedef unsigned long long x87_u64;
typedef long long x87_i64;
X87_FN int x87_clz(x87_u64 v) { return __clzll((long long)v); }
X87_FN x87_u64 x87_mulhi(x87_u64 a, x87_u64 b) { return __umul64hi(a, b); }
X87_FN x87_u64 x87_bits(double d) { return (x87_u64)__double_as_longlong(d); }
X87_FN double x87_double(x87_u64 b) { return __longlong_as_double((long long)b); }
#else
#include <bit>
#include <cstdint>
#define X87_FN inline
typedef std::uint64_t x87_u64;
typedef std::int64_t x87_i64;
X87_FN int x87_clz(x87_u64 v) { return v ? std::countl_zero(v) : 64; }
X87_FN x87_u64 x87_mulhi(x87_u64 a, x87_u64 b) { return (x87_u64)(((unsigned __int128)a * b) >> 64); }
X87_FN x87_u64 x87_bits(double d) { return std::bit_cast<x87_u64>(d); }
X87_FN double x87_double(x87_u64 b) { return std::bit_cast<double>(b); }
#endif

namespace x87 {

struct X { x87_u64 m; int e; int s; };   // zero: m == 0 (s == 0)

X87_FN X zero() { X r; r.m = 0; r.e = 0; r.s = 0; return r; }
X87_FN bool is_zero(X a) { return a.m == 0; }
X87_FN X neg(X a) { if (a.m != 0) a.s ^= 1; return a; }
X87_FN X abs(X a) { a.s = 0; return a; }

// Exact conversion of a finite binary64 value.
X87_FN X from_double(double d) {
    const x87_u64 bits = x87_bits(d);
    X r;
    r.s = int(bits >> 63);
    const int exponent = int((bits >> 52) & 0x7ff);
    const x87_u64 fraction = bits & 0x000fffffffffffffull;
    if (exponent == 0) {
        if (fraction == 0) return zero();
        const int shift = x87_clz(fraction);
        r.m = fraction << shift; r.e = -1074 - shift;
        return r;
    }
    r.m = (fraction | 0x0010000000000000ull) << 11;
    r.e = exponent - 1075 - 11;
    return r;
}
// Exact conversion of a non-negative integer.
X87_FN X from_u64(x87_u64 v) {
    if (v == 0) return zero();
    const int z = x87_clz(v);
    X r; r.m = v << z; r.e = -z; r.s = 0;
    return r;
}
// Exact scaling by 2^k.
X87_FN X scale2(X a, int k) { if (a.m != 0) a.e += k; return a; }

// Magnitude comparison |a| >= |b| (normalised significands).
X87_FN bool mag_ge(X a, X b) {
    if (b.m == 0) return true;
    if (a.m == 0) return false;
    if (a.e != b.e) return a.e > b.e;
    return a.m >= b.m;
}

// Rounds the normalised 192-bit value (w2 top bit set) times 2^(e - 128) to 64 bits.
X87_FN X round192(x87_u64 w2, x87_u64 w1, x87_u64 w0, int e, int s, bool sticky_in) {
    X r; r.s = s;
    x87_u64 m = w2;
    const bool guard = (w1 >> 63) != 0;
    const bool rest = (w1 << 1) != 0 || w0 != 0 || sticky_in;
    r.e = e;
    if (guard && (rest || (m & 1))) {
        m += 1;
        if (m == 0) { m = 0x8000000000000000ull; r.e += 1; }
    }
    r.m = m;
    return r;
}

// Correctly rounded sum.
X87_FN X add(X a, X b) {
    if (b.m == 0) return a;
    if (a.m == 0) return b;
    if (!mag_ge(a, b)) { const X t = a; a = b; b = t; }
    const int d = a.e - b.e;   // >= 0
    if (d >= 66) return a;     // |b| < ulp(a)/4: the rounded result is a (sum and difference)
    // a at bits 191..128 of the window, b shifted right by d (exact for d <= 65)
    x87_u64 b2, b1, b0;
    if (d == 0) { b2 = b.m; b1 = 0; b0 = 0; }
    else if (d < 64) { b2 = b.m >> d; b1 = b.m << (64 - d); b0 = 0; }
    else if (d == 64) { b2 = 0; b1 = b.m; b0 = 0; }
    else { b2 = 0; b1 = b.m >> 1; b0 = b.m << 63; }   // d == 65
    if (a.s == b.s) {
        // (a.m, 0, 0) + (b2, b1, b0)
        const x87_u64 r0 = b0, r1 = b1;
        const x87_u64 r2 = a.m + b2;
        if (r2 < a.m) {   // carry out: shift right by one
            const bool sticky = (r0 & 1) != 0;
            const x87_u64 n2 = 0x8000000000000000ull | (r2 >> 1), n1 = (r2 << 63) | (r1 >> 1), n0 = (r1 << 63) | (r0 >> 1);
            return round192(n2, n1, n0, a.e + 1, a.s, sticky);
        }
        return round192(r2, r1, r0, a.e, a.s, false);
    }
    // (a.m, 0, 0) - (b2, b1, b0), non-negative
    const x87_u64 r0 = 0 - b0;
    const x87_u64 borrow0 = b0 != 0 ? 1 : 0;
    const x87_u64 t1 = 0 - b1;
    const x87_u64 r1 = t1 - borrow0;
    const x87_u64 borrow1 = (b1 != 0 || (b1 == 0 && borrow0 != 0)) ? 1 : 0;
    const x87_u64 r2 = a.m - b2 - borrow1;
    if (r2 == 0 && r1 == 0 && r0 == 0) return zero();
    // normalise left
    int shift;
    x87_u64 n2 = r2, n1 = r1, n0 = r0;
    if (n2 == 0) {
        if (n1 == 0) { n2 = n0; n1 = 0; n0 = 0; shift = 128; }
        else { n2 = n1; n1 = n0; n0 = 0; shift = 64; }
    } else shift = 0;
    const int z = x87_clz(n2);
    if (z) { n2 = (n2 << z) | (n1 >> (64 - z)); n1 = (n1 << z) | (n0 >> (64 - z)); n0 <<= z; }
    shift += z;
    return round192(n2, n1, n0, a.e - shift, a.s, false);
}
X87_FN X sub(X a, X b) { return add(a, neg(b)); }

// Correctly rounded product.
X87_FN X mul(X a, X b) {
    if (a.m == 0 || b.m == 0) return zero();
    x87_u64 hi = x87_mulhi(a.m, b.m), lo = a.m * b.m;
    int e = a.e + b.e + 64;
    if ((hi >> 63) == 0) { hi = (hi << 1) | (lo >> 63); lo <<= 1; e -= 1; }
    return round192(hi, lo, 0, e, a.s ^ b.s, false);
}

// a <= b (values).
X87_FN bool le(X a, X b) {
    if (a.m == 0 && b.m == 0) return true;
    if (a.m == 0) return b.s == 0;
    if (b.m == 0) return a.s == 1;
    if (a.s != b.s) return a.s == 1;
    if (a.s == 0) return mag_ge(b, a);
    return mag_ge(a, b);
}
X87_FN bool lt(X a, X b) { return !le(b, a); }

// floor(a) as a signed integer (|a| < 2^62 in every use).
X87_FN x87_i64 floor_int(X a) {
    if (a.m == 0) return 0;
    x87_u64 ip; bool frac;
    if (a.e >= 0) { ip = a.m << a.e; frac = false; }
    else if (a.e <= -64) { ip = 0; frac = true; }
    else { ip = a.m >> (-a.e); frac = (a.m << (64 + a.e)) != 0; }
    if (a.s == 0) return x87_i64(ip);
    return -x87_i64(ip) - (frac ? 1 : 0);
}

// Binary64 value: nearest even; upward_if_below: then the sealed nextafter(+inf) step when the rounded
// magnitude is below the exact one (rounding toward +inf for positive values).
X87_FN double to_double(X x, bool upward_if_below) {
    if (x.m == 0) return 0.0;
    const int top = x.e + 63;
    x87_u64 q, rem, half, high_part;
    if (top >= -1022) {
        q = x.m >> 11; rem = x.m & 0x7ffull; half = 0x400ull;
        high_part = x87_u64(top + 1022) << 52;
    } else {
        const int shift = -(x.e + 1074);
        high_part = 0;
        if (shift > 64) { q = 0; rem = x.m; half = 0; }
        else if (shift == 64) { q = 0; rem = x.m; half = 0x8000000000000000ull; }
        else { q = x.m >> shift; rem = x.m & ((x87_u64(1) << shift) - 1); half = x87_u64(1) << (shift - 1); }
    }
    bool up;
    if (half == 0) up = false;
    else up = rem > half || (rem == half && (q & 1));
    x87_u64 bits = high_part + q + (up ? 1 : 0);
    if (top > 1023 || bits >= 0x7ff0000000000000ull) bits = 0x7ff0000000000000ull;
    if (upward_if_below && x.s == 0 && rem != 0 && !up && bits < 0x7ff0000000000000ull) bits += 1;
    if (x.s) bits |= 0x8000000000000000ull;
    return x87_double(bits);
}

// Interchange form (device arrays, kernel arguments): significand and exponent/sign word.
struct XBits { x87_u64 m; x87_i64 es; };
X87_FN XBits pack(X a) { XBits b; b.m = a.m; b.es = x87_i64(a.e) * 2 + x87_i64(a.s); return b; }
X87_FN X unpack(XBits b) {
    X a; a.m = b.m; a.s = int(b.es & 1); a.e = int((b.es - x87_i64(a.s)) / 2);
    if (a.m == 0) { a.e = 0; a.s = 0; }
    return a;
}

}  // namespace x87
