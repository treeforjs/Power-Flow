// forms_core.h -- native v7 producer forms, portable core.
//
// Compiled unchanged by the host C++ compiler (C++20, -ffp-contract=off) and by
// NVRTC (--fmad=false --prec-div=true --prec-sqrt=true --ftz=false). Integer
// arithmetic only, except the kinetic-energy formula (IEEE + - * / sqrt in a
// fixed order, no contraction), so the host and the device produce the same
// bits. No standard headers are needed under NVRTC.
//
// Contents
//  * exact superaccumulator for sums of exact products of binary64 values
//    (radix 2^32 digits held in int64, offset 2^-2148) and its correctly
//    rounded (round to nearest, ties to even) conversion to binary64;
//  * the kinetic energy of pic_kinetic_energy_J and the per-carrier change;
//  * SHA-256 (FIPS 180-4) streaming core over 32-bit big-endian words;
//  * record sizes and block-list geometry of the v7 Merkle forms.
// The definitions are in API-007.txt section 1.
#pragma once

#if defined(__CUDACC_RTC__) || defined(__CUDACC__)
#define PF7_HD __host__ __device__ inline
#define PF7_DEVICE 1
#define PF7_UNROLL _Pragma("unroll")
typedef unsigned long long pf7_u64;
typedef long long pf7_i64;
typedef unsigned int pf7_u32;
typedef unsigned char pf7_u8;
#else
#include <bit>
#include <cmath>
#include <cstdint>
#define PF7_HD inline
#define PF7_DEVICE 0
#define PF7_UNROLL
typedef std::uint64_t pf7_u64;
typedef std::int64_t pf7_i64;
typedef std::uint32_t pf7_u32;
typedef std::uint8_t pf7_u8;
#endif

namespace pf7 {

// ---------------------------------------------------------------- constants
// Digits: value = sum_k digit[k] * 2^(32k - 2148). 2^-2148 is the weight of
// the least significant bit of the product of two binary64 subnormals; the
// largest product is below 2^2048, i.e. below bit 4196; 134 digits (4288
// bits) leave 92 bits for carries and the sign.
constexpr int acc_digits = 134;
constexpr int acc_offset = 2148;
// A device or host accumulator absorbs at most 2^31 terms between
// normalisations (each term adds < 2^32 to a digit).
constexpr pf7_u64 acc_terms_per_normalisation = 1ull << 30;
// RN(299792458.0 * 299792458.0) (exact product 89875517873681764).
constexpr double c_squared = 89875517873681760.0;

// Block lists (API-007.txt 1.2): 32 records per chunk, 32 chunks (1024
// records) per leaf.
constexpr pf7_u64 chunk_records = 32;
constexpr pf7_u64 leaf_records = 1024;
constexpr pf7_u64 leaf_chunks = 32;

// Canonical record sizes (bytes, little-endian, packed).
constexpr pf7_u64 endpoint_bytes = 80; // u64 id, 3 f64 x, 3 f64 u, f64 q, f64 m, f64 w
constexpr pf7_u64 path_bytes = 64;     // u64 id, 3 f64 start, 3 f64 displacement, f64 flight
constexpr pf7_u64 sink_bytes = 136;    // u64 id, kind, surface, axis, f64 plane, fraction,
                                       // 3 f64 segment, 3 f64 hit, 3 f64 incident u,
                                       // f64 flight, f64 discarded flight
constexpr pf7_u64 birth_bytes = 24;    // u64 id, f64 birth time, f64 flight
constexpr pf7_u64 census_bytes = 96;   // u64 id, 3 f64 x, 3 f64 u, f64 q, f64 m, f64 w,
                                       // u64 packet key, u64 state index

// Status bits of an exact sum.
constexpr pf7_u32 sum_nonfinite_input = 1u;   // a factor was inf or NaN
constexpr pf7_u32 sum_overflow = 2u;          // |rounded sum| >= 2^1024
constexpr pf7_u32 sum_too_many_terms = 4u;    // launch above 2^31 terms
constexpr pf7_u32 sum_census_missing = 8u;    // split kinetic: ID not in the census

// ---------------------------------------------------------------- bit casts
PF7_HD pf7_u64 bits_of(double x) {
#if PF7_DEVICE
  return (pf7_u64)__double_as_longlong(x);
#else
  return std::bit_cast<pf7_u64>(x);
#endif
}
PF7_HD double double_of(pf7_u64 b) {
#if PF7_DEVICE
  return __longlong_as_double((long long)b);
#else
  return std::bit_cast<double>(b);
#endif
}
PF7_HD pf7_u64 mul_hi(pf7_u64 a, pf7_u64 b) {
#if PF7_DEVICE
  return __umul64hi(a, b);
#else
  __extension__ typedef unsigned __int128 wide;
  return (pf7_u64)(((wide)a * b) >> 64);
#endif
}
PF7_HD double fsqrt(double x) {
#if PF7_DEVICE
  return sqrt(x);
#else
  return std::sqrt(x);
#endif
}

// ---------------------------------------------------------------- exact terms
// One exact product a*b placed into the digit grid: digits d..d+4 receive the
// 32-bit chunks c[0..4] (with sign neg). zero: the product is exactly zero.
struct Term {
  pf7_u32 c[5];
  int d;
  int neg;
  int zero;
  int bad; // nonfinite factor
};
// Decomposes a finite binary64 x = (-1)^s m 2^e (m < 2^53, e = LSB exponent).
PF7_HD int decompose(double x, pf7_u64 &m, int &e, int &s) {
  const pf7_u64 b = bits_of(x);
  const int E = (int)((b >> 52) & 0x7ff);
  s = (int)(b >> 63);
  m = b & ((1ull << 52) - 1);
  if (E == 0x7ff)
    return 0;
  if (E) {
    m |= 1ull << 52;
    e = E - 1075;
  } else
    e = -1074;
  return 1;
}
PF7_HD Term product_term(double a, double b) {
  Term t;
  for (int k = 0; k < 5; ++k)
    t.c[k] = 0;
  t.d = 0;
  t.neg = 0;
  t.zero = 1;
  t.bad = 0;
  pf7_u64 ma, mb;
  int ea, eb, sa, sb;
  const int fa = decompose(a, ma, ea, sa), fb = decompose(b, mb, eb, sb);
  if (!fa || !fb) {
    t.bad = 1;
    return t;
  }
  if (ma == 0 || mb == 0)
    return t;
  const pf7_u64 lo = ma * mb, hi = mul_hi(ma, mb); // hi < 2^42
  const int pos = ea + eb + acc_offset;            // 0 .. 4090
  const int r = pos & 31;
  t.d = pos >> 5;
  const pf7_u64 v0 = lo << r;
  const pf7_u64 v1 = r ? ((hi << r) | (lo >> (64 - r))) : hi;
  const pf7_u64 v2 = r ? (hi >> (64 - r)) : 0;
  t.c[0] = (pf7_u32)(v0 & 0xffffffffull);
  t.c[1] = (pf7_u32)(v0 >> 32);
  t.c[2] = (pf7_u32)(v1 & 0xffffffffull);
  t.c[3] = (pf7_u32)(v1 >> 32);
  t.c[4] = (pf7_u32)v2;
  t.neg = sa ^ sb;
  t.zero = 0;
  return t;
}
PF7_HD void add_term(pf7_i64 *digit, const Term &t) {
  if (t.zero || t.bad)
    return;
  for (int k = 0; k < 5; ++k)
    digit[t.d + k] += t.neg ? -(pf7_i64)t.c[k] : (pf7_i64)t.c[k];
}
// Carry propagation: digits 0..132 in [0, 2^32), digit 133 signed.
PF7_HD void normalise(pf7_i64 *digit) {
  for (int k = 0; k + 1 < acc_digits; ++k) {
    const pf7_i64 carry = digit[k] >> 32; // arithmetic shift = floor division
    digit[k] &= 0xffffffffll;
    digit[k + 1] += carry;
  }
}
PF7_HD int clz32(pf7_u32 x) {
#if PF7_DEVICE
  return __clz((int)x);
#else
  return x ? std::countl_zero(x) : 32;
#endif
}
PF7_HD pf7_u64 bit_at(const pf7_i64 *digit, int i) {
  return ((pf7_u64)digit[i >> 5] >> (i & 31)) & 1ull;
}
// Correctly rounded (nearest, ties to even) binary64 value of the exact sum.
// digit is normalised in place (and negated when the sum is negative). An
// exactly zero sum gives +0.0. Overflow gives +-inf and sets *status bit
// sum_overflow.
PF7_HD double round_sum(pf7_i64 *digit, pf7_u32 *status) {
  normalise(digit);
  int neg = 0;
  if (digit[acc_digits - 1] < 0) {
    neg = 1;
    for (int k = 0; k < acc_digits; ++k)
      digit[k] = -digit[k];
    normalise(digit);
  }
  int h = acc_digits - 1;
  while (h >= 0 && digit[h] == 0)
    --h;
  if (h < 0)
    return 0.0;
  const pf7_u64 sign = (pf7_u64)neg << 63;
  const pf7_u64 inf = sign | (0x7ffull << 52);
  if (h == acc_digits - 1 || h >= 100) { // >= 2^(3200 - 2148)
    *status |= sum_overflow;
    return double_of(inf);
  }
  const int P = 32 * h + 31 - clz32((pf7_u32)digit[h]);
  if (P >= 1024 + acc_offset) {
    *status |= sum_overflow;
    return double_of(inf);
  }
  int L = P - 52;
  if (L < 1074) // binary64 LSB 2^-1074 sits at grid position 1074
    L = 1074;
  pf7_u64 m = 0;
  for (int i = P; i >= L; --i)
    m = (m << 1) | bit_at(digit, i);
  const int g = L - 1;
  const pf7_u64 guard = bit_at(digit, g);
  int sticky = 0;
  {
    const int dg = g >> 5, rb = g & 31;
    if (rb && ((pf7_u64)digit[dg] & ((1ull << rb) - 1)))
      sticky = 1;
    for (int k = 0; k < dg && !sticky; ++k)
      if (digit[k])
        sticky = 1;
  }
  if (guard && (sticky || (m & 1ull))) {
    ++m;
    if (m == (1ull << 53)) {
      m = 1ull << 52;
      ++L;
    }
  }
  if (m >= (1ull << 52)) {
    const pf7_u64 E = (pf7_u64)(L - 1073);
    if (E >= 2047) {
      *status |= sum_overflow;
      return double_of(inf);
    }
    return double_of(sign | (E << 52) | (m - (1ull << 52)));
  }
  return double_of(sign | m); // subnormal (L == 1074)
}

// ---------------------------------------------------------------- kinetic energy
// pic_kinetic_energy_J (native mfem_gap_pic.cpp): m |u|^2 / (gamma + 1),
// gamma = sqrt(1 + |u|^2 / c^2), |u|^2 = (u0 u0 + u1 u1) + u2 u2.
PF7_HD double kinetic_energy(double mass, double u0, double u1, double u2) {
  const double n = u0 * u0 + u1 * u1 + u2 * u2;
  const double g = fsqrt(1 + n / c_squared);
  return mass * n / (g + 1);
}
// RN(KE(end) - KE(start)) of one carrier.
PF7_HD double kinetic_change(double mass, const double *u_start, const double *u_end) {
  return kinetic_energy(mass, u_end[0], u_end[1], u_end[2]) -
         kinetic_energy(mass, u_start[0], u_start[1], u_start[2]);
}

// ---------------------------------------------------------------- SHA-256
struct Sha {
  pf7_u32 h[8];
  pf7_u32 w[16];
  pf7_u32 pos;   // bytes in the current block (0..63)
  pf7_u64 total; // bytes fed
};
PF7_HD pf7_u32 rotr(pf7_u32 x, int n) { return (x >> n) | (x << (32 - n)); }
// FIPS 180-4 round constants (first 32 bits of the fractional parts of the
// cube roots of the first 64 primes).
#if PF7_DEVICE
__constant__ const pf7_u32 sha_k_table[64] = {
      0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
      0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
      0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
      0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
      0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
      0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
      0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
      0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
#else
inline constexpr pf7_u32 sha_k_table[64] = {
      0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
      0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
      0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
      0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
      0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
      0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
      0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
      0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
#endif
PF7_HD void sha_compress(Sha &s) {
  pf7_u32 W[64];
  PF7_UNROLL
  for (int i = 0; i < 16; ++i)
    W[i] = s.w[i];
  PF7_UNROLL
  for (int i = 16; i < 64; ++i) {
    const pf7_u32 s0 = rotr(W[i - 15], 7) ^ rotr(W[i - 15], 18) ^ (W[i - 15] >> 3);
    const pf7_u32 s1 = rotr(W[i - 2], 17) ^ rotr(W[i - 2], 19) ^ (W[i - 2] >> 10);
    W[i] = W[i - 16] + s0 + W[i - 7] + s1;
  }
  pf7_u32 a = s.h[0], b = s.h[1], c = s.h[2], d = s.h[3], e = s.h[4], f = s.h[5], g = s.h[6], h = s.h[7];
  PF7_UNROLL
  for (int i = 0; i < 64; ++i) {
    const pf7_u32 S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
    const pf7_u32 ch = (e & f) ^ (~e & g);
    const pf7_u32 t1 = h + S1 + ch + sha_k_table[i] + W[i];
    const pf7_u32 S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
    const pf7_u32 mj = (a & b) ^ (a & c) ^ (b & c);
    const pf7_u32 t2 = S0 + mj;
    h = g;
    g = f;
    f = e;
    e = d + t1;
    d = c;
    c = b;
    b = a;
    a = t1 + t2;
  }
  s.h[0] += a;
  s.h[1] += b;
  s.h[2] += c;
  s.h[3] += d;
  s.h[4] += e;
  s.h[5] += f;
  s.h[6] += g;
  s.h[7] += h;
  PF7_UNROLL
  for (int i = 0; i < 16; ++i)
    s.w[i] = 0;
}
PF7_HD void sha_init(Sha &s) {
  s.h[0] = 0x6a09e667u;
  s.h[1] = 0xbb67ae85u;
  s.h[2] = 0x3c6ef372u;
  s.h[3] = 0xa54ff53au;
  s.h[4] = 0x510e527fu;
  s.h[5] = 0x9b05688cu;
  s.h[6] = 0x1f83d9abu;
  s.h[7] = 0x5be0cd19u;
  for (int i = 0; i < 16; ++i)
    s.w[i] = 0;
  s.pos = 0;
  s.total = 0;
}
// Invariant: message words at and after the current byte position are zero.
PF7_HD void sha_byte(Sha &s, pf7_u32 b) {
  s.w[s.pos >> 2] |= (b & 0xffu) << (24 - 8 * (s.pos & 3));
  ++s.total;
  if (++s.pos == 64) {
    sha_compress(s);
    s.pos = 0;
  }
}
// Four message bytes given as one big-endian word, at any byte alignment.
PF7_HD void sha_word(Sha &s, pf7_u32 x) {
  const pf7_u32 o = s.pos & 3, i = s.pos >> 2;
  s.total += 4;
  if (o == 0) {
    s.w[i] = x;
    s.pos += 4;
    if (s.pos == 64) {
      sha_compress(s);
      s.pos = 0;
    }
    return;
  }
  const int sh = 8 * (int)o;
  s.w[i] |= x >> sh;
  const pf7_u32 carry = x << (32 - sh);
  if (i == 15) {
    sha_compress(s);
    s.w[0] = carry;
    s.pos = s.pos + 4 - 64;
  } else {
    s.w[i + 1] = carry;
    s.pos += 4;
  }
}
PF7_HD pf7_u32 bswap32(pf7_u32 x) {
  return (x >> 24) | ((x >> 8) & 0xff00u) | ((x << 8) & 0xff0000u) | (x << 24);
}
// u64 in canonical little-endian byte order.
PF7_HD void sha_u64(Sha &s, pf7_u64 v) {
  sha_word(s, bswap32((pf7_u32)(v & 0xffffffffull)));
  sha_word(s, bswap32((pf7_u32)(v >> 32)));
}
PF7_HD void sha_f64(Sha &s, double x) { sha_u64(s, bits_of(x)); }
// Canonical text: u64 byte length, then the bytes.
PF7_HD void sha_text(Sha &s, const char *t, pf7_u32 n) {
  sha_u64(s, n);
  for (pf7_u32 i = 0; i < n; ++i)
    sha_byte(s, (pf7_u8)t[i]);
}
// 32 raw digest bytes given as the 8 state words of another hash.
PF7_HD void sha_digest_words(Sha &s, const pf7_u32 *d) {
  for (int i = 0; i < 8; ++i)
    sha_word(s, d[i]);
}
PF7_HD void sha_finish(Sha &s, pf7_u32 *out) {
  const pf7_u64 bits = s.total * 8;
  sha_byte(s, 0x80);
  while (s.pos != 56)
    sha_byte(s, 0);
  sha_word(s, (pf7_u32)(bits >> 32));
  sha_word(s, (pf7_u32)(bits & 0xffffffffull));
  for (int i = 0; i < 8; ++i)
    out[i] = s.h[i];
}

// ---------------------------------------------------------------- tags
// Section tags of the v7 block lists; level tags are "<tag>/chunk",
// "<tag>/block" and "<tag>/root" (API-007.txt 1.2).
enum Section : pf7_u32 {
  section_endpoints = 0,
  section_paths = 1,
  section_sinks = 2,
  section_births = 3,
  section_census = 4
};
constexpr int level_chunk = 0, level_block = 1, level_root = 2;
// Full tag text and its length (no terminator is hashed).
PF7_HD const char *section_tag(pf7_u32 section, int level, pf7_u32 &n) {
  // 5 sections x 3 levels.
  const char *const tags[15] = {
      "native-v7/motion-endpoints/chunk", "native-v7/motion-endpoints/block", "native-v7/motion-endpoints/root",
      "native-v7/motion-paths/chunk",     "native-v7/motion-paths/block",     "native-v7/motion-paths/root",
      "native-v7/motion-sinks/chunk",     "native-v7/motion-sinks/block",     "native-v7/motion-sinks/root",
      "native-v7/motion-births/chunk",    "native-v7/motion-births/block",    "native-v7/motion-births/root",
      "native-v7/census-carriers/chunk",  "native-v7/census-carriers/block",  "native-v7/census-carriers/root"};
  const char *t = tags[section * 3 + (pf7_u32)level];
  n = 0;
  while (t[n])
    ++n;
  return t;
}
PF7_HD pf7_u64 section_record_bytes(pf7_u32 section) {
  return section == section_endpoints ? endpoint_bytes
         : section == section_paths   ? path_bytes
         : section == section_sinks   ? sink_bytes
         : section == section_births  ? birth_bytes
                                      : census_bytes;
}
// chunk_j = SHA256(text <tag>/chunk, u64 j, u64 m, m records), records given
// as little-endian u64 words (record size is a multiple of 8).
PF7_HD void chunk_digest(pf7_u32 section, pf7_u64 j, pf7_u64 m, const pf7_u64 *words,
                         pf7_u32 *out) {
  Sha s;
  sha_init(s);
  pf7_u32 n;
  const char *t = section_tag(section, level_chunk, n);
  sha_text(s, t, n);
  sha_u64(s, j);
  sha_u64(s, m);
  const pf7_u64 count = m * (section_record_bytes(section) / 8);
  for (pf7_u64 i = 0; i < count; ++i)
    sha_u64(s, words[i]);
  sha_finish(s, out);
}
// block_k = SHA256(text <tag>/block, u64 k, u64 n_k, u64 q_k, q_k chunk digests).
PF7_HD void block_digest(pf7_u32 section, pf7_u64 k, pf7_u64 records, pf7_u64 chunks,
                         const pf7_u32 *chunk_words, pf7_u32 *out) {
  Sha s;
  sha_init(s);
  pf7_u32 n;
  const char *t = section_tag(section, level_block, n);
  sha_text(s, t, n);
  sha_u64(s, k);
  sha_u64(s, records);
  sha_u64(s, chunks);
  for (pf7_u64 q = 0; q < chunks; ++q)
    sha_digest_words(s, chunk_words + 8 * q);
  sha_finish(s, out);
}
// native v5 tree node and empty tree (mfem_pic_canonical.cpp).
PF7_HD void merkle_node(const pf7_u32 *a, const pf7_u32 *b, pf7_u32 *out) {
  Sha s;
  sha_init(s);
  sha_text(s, "native-v5/merkle-node", 21);
  sha_digest_words(s, a);
  sha_digest_words(s, b);
  sha_finish(s, out);
}
PF7_HD void merkle_empty(pf7_u32 *out) {
  Sha s;
  sha_init(s);
  sha_text(s, "native-v5/merkle-empty", 22);
  sha_finish(s, out);
}
// root = SHA256(text <tag>/root, u64 1024, u64 records, u64 leaves, tree).
PF7_HD void section_root(pf7_u32 section, pf7_u64 records, pf7_u64 leaves, const pf7_u32 *tree,
                         pf7_u32 *out) {
  Sha s;
  sha_init(s);
  pf7_u32 n;
  const char *t = section_tag(section, level_root, n);
  sha_text(s, t, n);
  sha_u64(s, leaf_records);
  sha_u64(s, records);
  sha_u64(s, leaves);
  sha_digest_words(s, tree);
  sha_finish(s, out);
}

// ---------------------------------------------------------------- census
// Failure codes of the canonical census check (first failing record wins).
enum CensusCode : pf7_u32 {
  census_ok = 0,
  census_id_zero = 1,
  census_not_ascending = 2,   // id <= previous id (unsorted or duplicate)
  census_key_invalid = 3,     // packet key 0 or above the carrier's own id
  census_key_absent = 4,      // no census carrier with id == packet key
  census_key_not_first = 5,   // the key carrier belongs to another packet
  census_divergent = 6,       // differs from its packet's key carrier
  census_state_range = 7      // state index >= number of states
};
// Record words: 0 id, 1-3 x, 4-6 u, 7 q, 8 m, 9 w, 10 packet key, 11 state.
PF7_HD pf7_u64 census_find(const pf7_u64 *rec, pf7_u64 n, pf7_u64 id) {
  pf7_u64 lo = 0, hi = n;
  while (lo < hi) {
    const pf7_u64 mid = lo + (hi - lo) / 2;
    if (rec[12 * mid] < id)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo;
}
// Check of record i of an ascending census (word view, 12 words per record).
PF7_HD pf7_u32 census_check(const pf7_u64 *rec, pf7_u64 n, pf7_u64 i, pf7_u64 states) {
  const pf7_u64 *r = rec + 12 * i;
  const pf7_u64 id = r[0], key = r[10];
  if (id == 0)
    return census_id_zero;
  if (i > 0 && rec[12 * (i - 1)] >= id)
    return census_not_ascending;
  if (key == 0 || key > id)
    return census_key_invalid;
  if (r[11] >= states)
    return census_state_range;
  if (key == id)
    return census_ok;
  const pf7_u64 k = census_find(rec, n, key);
  if (k >= n || rec[12 * k] != key)
    return census_key_absent;
  const pf7_u64 *f = rec + 12 * k;
  if (f[10] != key)
    return census_key_not_first;
  for (int w = 1; w <= 8; ++w) // position, proper velocity, charge, mass
    if (f[w] != r[w])
      return census_divergent;
  if (f[11] != r[11])
    return census_divergent;
  return census_ok;
}
// ---------------------------------------------------------------- neutral stream
// Native v7 port of the sealed resident-neutral stream (Hall 8dddec0
// kinetic_boundary.cpp prepare_neutral_stream and first_kinetic_hit; native
// v6 native_first_hit): same operands, operation order, tie rule and check
// order. The only change: the interval bound dt <= 1e-14 is the parameter
// dt_max (native v7 configuration, default 1e-14; API-007.txt 1.4).
// Policy foils: left, right, bottom, top.
struct Policy {
  double width, gap, short_z, feed_z;
  double lower[4], upper[4];
};
struct Hit {
  int sink;        // hall_driver::KineticSink value, -1 none
  double fraction; // 1 when none
  int axis;
  double plane;
};
PF7_HD int finite(double x) { return ((bits_of(x) >> 52) & 0x7ff) != 0x7ff; }
PF7_HD double plus_inf() { return double_of(0x7ff0000000000000ull); }
// require_kinetic_interior without the policy guard: 1 when x is strictly
// inside the kinetic box and not in or on a foil.
PF7_HD int kinetic_interior(const Policy &p, const double *x) {
  if (!(x[0] > 0 && x[0] < p.width && x[1] > 0 && x[1] < p.gap && x[2] > p.short_z && x[2] < p.feed_z))
    return 0;
  for (int b = 0; b < 2; ++b) {
    const double *f = b ? p.upper : p.lower;
    if (x[0] >= f[0] && x[0] <= f[1] && x[1] >= f[2] && x[1] <= f[3])
      return 0;
  }
  return 1;
}
PF7_HD Hit first_hit(const Policy &p, const double *x, const double *D) {
  Hit hit;
  hit.sink = -1;
  hit.fraction = 1;
  hit.axis = -1;
  hit.plane = 0;
  for (int axis = 0; axis < 2; ++axis)
    if (D[axis] != 0) {
      const double plane = D[axis] > 0 ? (axis ? p.gap : p.width) : 0;
      const double f = (plane - x[axis]) / D[axis];
      if (!(f < 0 || f > hit.fraction)) {
        hit.sink = 2;
        hit.fraction = f;
        hit.axis = axis;
        hit.plane = plane;
      }
    }
  for (int body = 0; body < 2; ++body) {
    const double *foil = body ? p.upper : p.lower;
    double enter = -plus_inf(), leave = plus_inf(), plane = 0;
    int enter_axis = -1, possible = 1;
    for (int axis = 0; axis < 2; ++axis) {
      const double low = axis ? foil[2] : foil[0], high = axis ? foil[3] : foil[1];
      if (D[axis] == 0) {
        if (x[axis] < low || x[axis] > high)
          possible = 0;
        continue;
      }
      double a = (low - x[axis]) / D[axis], b = (high - x[axis]) / D[axis], entry_plane = low;
      if (a > b) {
        const double t = a;
        a = b;
        b = t;
        entry_plane = high;
      }
      if (a >= enter) {
        enter = a;
        enter_axis = axis;
        plane = entry_plane;
      }
      leave = b < leave ? b : leave; // std::min(leave, b)
    }
    if (possible && leave > 0 && enter >= 0 && enter <= leave && enter_axis >= 0 &&
        !(enter < 0 || enter > hit.fraction)) {
      hit.sink = body ? 1 : 0;
      hit.fraction = enter;
      hit.axis = enter_axis;
      hit.plane = plane;
    }
  }
  if (D[2] != 0) {
    const double plane = D[2] > 0 ? p.feed_z : p.short_z, f = (plane - x[2]) / D[2];
    const double xx = x[0] + f * D[0], yy = x[1] + f * D[1];
    const int bridge = D[2] < 0 && xx >= p.lower[0] && xx <= p.lower[1] && yy >= p.lower[3] && yy <= p.upper[2];
    if (!(f < 0 || f > hit.fraction)) {
      hit.sink = bridge ? 3 : 4;
      hit.fraction = f;
      hit.axis = 2;
      hit.plane = plane;
    }
  }
  return hit;
}
enum NeutralCode : pf7_u32 {
  neutral_ok = 0,
  neutral_interval = 1,        // "neutral authoritative interval"
  neutral_identity = 2,        // "neutral stream identity/species"
  neutral_momentum = 3,        // "neutral momentum finite"
  neutral_resolved = 4,        // "neutral momentum outside resolved range"
  neutral_boundary_owner = 5,  // "particle requires native boundary owner"
  neutral_in_metal = 6,        // "particle in/on native metal"
  neutral_displacement = 7,    // "nonfinite kinetic displacement"
  neutral_incident = 8         // "neutral incident kinetic operand"
};
struct NeutralStream {
  double endpoint[3]; // native position at the end of the flight (hit point when hit)
  double D[3];        // dt u / gamma
  Hit hit;
  double flight_s, remaining_s, incident_kinetic_J;
};
// x, u: native position and proper velocity. The kinetic policy itself must
// have passed the sealed policy guard (checked on the host).
PF7_HD pf7_u32 neutral_stream(const Policy &p, pf7_u64 id, const double *x, const double *u, double charge,
                              double mass, double weight, double begin, double end, double dt, double dt_max,
                              NeutralStream &out) {
  if (!(finite(begin) && begin >= 0 && finite(end) && finite(dt) && dt > 0 && dt <= dt_max && end > begin &&
        end == begin + dt))
    return neutral_interval;
  if (!(id > 0 && charge == 0 && finite(weight) && weight > 0 && finite(mass) && mass > 0))
    return neutral_identity;
  double u2 = 0;
  for (int d = 0; d < 3; ++d) {
    if (!finite(u[d]))
      return neutral_momentum;
    u2 += u[d] * u[d];
  }
  const double gamma = fsqrt(1 + u2 / c_squared);
  if (!finite(gamma))
    return neutral_resolved;
  for (int d = 0; d < 3; ++d)
    out.D[d] = dt * u[d] / gamma;
  if (!(x[0] > 0 && x[0] < p.width && x[1] > 0 && x[1] < p.gap && x[2] > p.short_z && x[2] < p.feed_z))
    return neutral_boundary_owner;
  if (!kinetic_interior(p, x))
    return neutral_in_metal;
  for (int d = 0; d < 3; ++d)
    if (!finite(out.D[d]))
      return neutral_displacement;
  out.hit = first_hit(p, x, out.D);
  out.flight_s = dt * out.hit.fraction;
  out.remaining_s = dt - out.flight_s;
  for (int d = 0; d < 3; ++d)
    out.endpoint[d] = x[d] + out.hit.fraction * out.D[d];
  if (out.hit.sink >= 0)
    out.endpoint[out.hit.axis] = out.hit.plane;
  out.incident_kinetic_J = weight * kinetic_energy(mass, u[0], u[1], u[2]);
  if (!finite(out.incident_kinetic_J))
    return neutral_incident;
  return neutral_ok;
}
} // namespace pf7
