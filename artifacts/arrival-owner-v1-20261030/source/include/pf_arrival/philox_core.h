// philox_core.h -- counter-based random numbers of the arrival owner v1.
//
// Compiled twice from this one file: by the host C++ compiler
// (-ffp-contract=off -fno-fast-math) and by NVRTC (--fmad=false
// --prec-div=true --prec-sqrt=true --ftz=false). Only 32/64-bit unsigned
// integer arithmetic, one exact integer-to-binary64 conversion and one
// multiplication by a power of two are used, so the two builds produce the
// same bits by construction (checked by tests/gpu_parity.cpp).
//
// Generator: Philox4x32-10 (J. K. Salmon, M. A. Moraes, R. O. Dror,
// D. E. Shaw, "Parallel random numbers: as easy as 1, 2, 3", SC11), the
// exact round function, multipliers and Weyl key increments of Random123.
// Known-answer vectors: tests/test_core.cpp (Random123 kat_vectors).
//
// Stream layout v1 (normative text in API-v1.txt):
//   key     k0 = seed bits 0..31, k1 = seed bits 32..63
//   counter c0 = block index b (32 bits)
//           c1 = candidate index (32 bits)
//           c2 = stream id (bits 0..15) | purpose tag (bits 16..23)
//                | step bits 32..39 (bits 24..31)
//           c3 = step bits 0..31
//   domain  step < 2^40, stream id < 2^16, purpose tag 1..255,
//           candidate < 2^32, block < 2^32; the map is injective on it.
//   words   block b = Philox4x32-10(counter, key) = (w0, w1, w2, w3) yields
//           R(2b) = (w1<<32 | w0) and R(2b+1) = (w3<<32 | w2).
// No state is carried between draws: every value is a function of the key
// tuple and the word index only.
#pragma once

#if defined(__CUDACC_RTC__)
#define PFA_FN __device__ inline
typedef unsigned int pfa_u32;
typedef unsigned long long pfa_u64;
PFA_FN pfa_u64 pfa_bits(double x) { return (pfa_u64)__double_as_longlong(x); }
#else
#include <bit>
#include <cstdint>
#define PFA_FN inline
typedef std::uint32_t pfa_u32;
typedef std::uint64_t pfa_u64;
PFA_FN pfa_u64 pfa_bits(double x) { return std::bit_cast<pfa_u64>(x); }
static_assert(sizeof(double) == 8 && sizeof(pfa_u64) == 8 && sizeof(pfa_u32) == 4);
#endif

namespace pfa {

struct Block { pfa_u32 v[4]; };
struct Key2 { pfa_u32 v[2]; };

// Philox constants (Random123 philox.h).
#define PFA_PHILOX_M0 0xD2511F53u
#define PFA_PHILOX_M1 0xCD9E8D57u
#define PFA_PHILOX_W0 0x9E3779B9u
#define PFA_PHILOX_W1 0xBB67AE85u

PFA_FN void mulhilo32(pfa_u32 a, pfa_u32 b, pfa_u32& hi, pfa_u32& lo) {
    const pfa_u64 p = (pfa_u64)a * (pfa_u64)b;
    hi = (pfa_u32)(p >> 32);
    lo = (pfa_u32)p;
}

PFA_FN Block philox4x32_round(Block c, Key2 k) {
    pfa_u32 hi0, lo0, hi1, lo1;
    mulhilo32(PFA_PHILOX_M0, c.v[0], hi0, lo0);
    mulhilo32(PFA_PHILOX_M1, c.v[2], hi1, lo1);
    Block r;
    r.v[0] = hi1 ^ c.v[1] ^ k.v[0];
    r.v[1] = lo1;
    r.v[2] = hi0 ^ c.v[3] ^ k.v[1];
    r.v[3] = lo0;
    return r;
}

// Ten rounds; the key is bumped by the Weyl increments before rounds 2..10.
PFA_FN Block philox4x32_10(Block c, Key2 k) {
    c = philox4x32_round(c, k);
    for (int i = 1; i < 10; ++i) {
        k.v[0] += PFA_PHILOX_W0;
        k.v[1] += PFA_PHILOX_W1;
        c = philox4x32_round(c, k);
    }
    return c;
}

// Purpose tags of stream layout v1. 0 is invalid. 0xF0..0xFF are reserved for
// qualification fixtures and never used by the production policy.
#define PFA_PURPOSE_COUNT 1u
#define PFA_PURPOSE_ELECTRON_MARKS 2u
#define PFA_PURPOSE_NEUTRAL_MARKS 3u

#define PFA_STEP_LIMIT 0x10000000000ull   // 2^40
#define PFA_STREAM_LIMIT 0x10000u         // 2^16

struct StreamKey {
    pfa_u64 seed;
    pfa_u64 step;        // accepted step index, < 2^40
    pfa_u32 stream;      // source stream ("patch") id, < 2^16
    pfa_u32 candidate;   // candidate index within (step, stream, purpose)
    pfa_u32 purpose;     // purpose tag 1..255
};

// Caller guarantees the domain (checked on the host by the policy before any draw).
PFA_FN Block stream_block(const StreamKey& s, pfa_u32 block) {
    Key2 k;
    k.v[0] = (pfa_u32)s.seed;
    k.v[1] = (pfa_u32)(s.seed >> 32);
    Block c;
    c.v[0] = block;
    c.v[1] = s.candidate;
    c.v[2] = (s.stream & 0xFFFFu) | ((s.purpose & 0xFFu) << 16) | ((pfa_u32)((s.step >> 32) & 0xFFu) << 24);
    c.v[3] = (pfa_u32)s.step;
    return philox4x32_10(c, k);
}
PFA_FN pfa_u64 low_word(const Block& b) { return ((pfa_u64)b.v[1] << 32) | (pfa_u64)b.v[0]; }
PFA_FN pfa_u64 high_word(const Block& b) { return ((pfa_u64)b.v[3] << 32) | (pfa_u64)b.v[2]; }
// R(j): the j-th 64-bit word of a stream.
PFA_FN pfa_u64 stream_word(const StreamKey& s, pfa_u32 j) {
    const Block b = stream_block(s, j >> 1);
    return (j & 1u) ? high_word(b) : low_word(b);
}

// Exact integer-to-binary64 conversions (no rounding occurs in either):
//  u01_open(R)      = (2*floor(R/2^12) + 1) * 2^-53, open interval (0,1):
//                     2^52 equally likely values 2^-53, 3*2^-53, ..., 1-2^-53,
//                     symmetric under u -> 1-u, mean exactly 1/2.
//  u01_half_open(R) = floor(R/2^11) * 2^-53, half-open [0, 1-2^-53]:
//                     2^53 equally likely values (not used by policy v1).
// The integer operand has at most 53 significant bits, so the conversion to
// binary64 is exact; the product with 2^-53 is an exponent change only.
PFA_FN double u01_open(pfa_u64 r) { return (double)(((r >> 12) << 1) | 1ull) * 0x1p-53; }
PFA_FN double u01_half_open(pfa_u64 r) { return (double)(r >> 11) * 0x1p-53; }

// Digit j (j = 0..16, most significant first) of floor(p * 2^1088) in base 2^64
// for a finite binary64 p in [0,1). Every bit of such p lies in 2^-1..2^-1074,
// so the 17 digits hold p exactly.
#define PFA_FRACTION_DIGITS 17
PFA_FN pfa_u64 fraction_digit(double p, int j) {
    const pfa_u64 b = pfa_bits(p);
    const int ex = (int)((b >> 52) & 0x7FFu);
    pfa_u64 m = b & 0xFFFFFFFFFFFFFull;
    int s;   // p = m * 2^(s - 1088)
    if (ex) {
        m |= 1ull << 52;
        s = ex - 1075 + 1088;
    } else {
        s = 14;   // -1074 + 1088
    }
    if (m == 0) return 0;
    const int lo = 1024 - 64 * j;   // digit j covers bit positions [lo, lo + 64)
    if (s + 53 <= lo || s >= lo + 64) return 0;
    return s >= lo ? (m << (s - lo)) : (m >> (lo - s));
}

// Exact Bernoulli(p) for binary64 p in [0,1): returns 1 iff U < p where
// U = 0.R(first) R(first+1) ... (base 2^64) is uniform on [0,1). Digits are
// compared most significant first; the first unequal digit decides. If all 17
// digits of p are matched, U >= p (p has no further bits). Hence
// P(1) = p exactly. p == 0 returns 0 without consuming words.
// The digit comparison is written over any word source w(j) so that tests can
// drive it with chosen word sequences; production uses the keyed stream words.
template <class Words>
PFA_FN int bernoulli_digits(const Words& w, double p) {
    if (pfa_bits(p) == 0) return 0;
    for (int j = 0; j < PFA_FRACTION_DIGITS; ++j) {
        const pfa_u64 r = w(j);
        const pfa_u64 c = fraction_digit(p, j);
        if (r < c) return 1;
        if (r > c) return 0;
    }
    return 0;
}
struct StreamWords {
    StreamKey s;
    pfa_u32 first;
    PFA_FN pfa_u64 operator()(int j) const { return stream_word(s, first + (pfa_u32)j); }
};
PFA_FN int bernoulli_exact(const StreamKey& s, double p, pfa_u32 first) {
    StreamWords w;
    w.s = s;
    w.first = first;
    return bernoulli_digits(w, p);
}

// Stochastic rounding of a finite binary64 x in [0, 2^52): floor(x) + Bernoulli(x - floor(x)).
// floor and the fraction are exact (x - n is exact by Sterbenz for n >= 1, and
// equals x for n = 0), so E[result] = x exactly for the given binary64 x.
PFA_FN pfa_u64 stochastic_round(const StreamKey& s, double x) {
    const pfa_u64 n = (pfa_u64)x;
    const double fraction = x - (double)n;
    return n + (pfa_u64)bernoulli_exact(s, fraction, 0u);
}

// Marks of one electron birth: open uniforms R(0..3) of its marks stream.
//   [0] normal_uniform  [1] transverse_uniform  [2] azimuth_uniform  [3] time_uniform
PFA_FN void electron_marks(const StreamKey& s, double* out4) {
    const Block a = stream_block(s, 0u), b = stream_block(s, 1u);
    out4[0] = u01_open(low_word(a));
    out4[1] = u01_open(high_word(a));
    out4[2] = u01_open(low_word(b));
    out4[3] = u01_open(high_word(b));
}
// Marks of one neutral leaf: open uniforms R(0..5) of its marks stream.
//   [0] normal [1] radius [2] angle [3] along [4] z [5] time
PFA_FN void neutral_marks(const StreamKey& s, double* out6) {
    for (pfa_u32 b = 0; b < 3u; ++b) {
        const Block x = stream_block(s, b);
        out6[2 * b] = u01_open(low_word(x));
        out6[2 * b + 1] = u01_open(high_word(x));
    }
}

}  // namespace pfa
