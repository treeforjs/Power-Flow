#pragma once
// Exact sum of nonnegative finite binary64 values (coordinator v7.1): an integer in units of
// 2^-1074 held in 35 little-endian 64-bit limbs (2240 bits: the binary64 range, 2098 bits, plus
// 142 bits of carry room). add() is exact (one shifted 53-bit significand, carries propagated),
// equality is exact, value() rounds once to nearest, ties to even (gradual underflow), the same
// value as native's ExactParticleCount::rounded_for_rate() of the same exact sum. It replaces
// the per-carrier ArenaCount credits (BIGNUM allocations) of the coordinator ledger totals.
#include "raikou/mfem_pic_full_time.h"
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace raikou::native_water::v7_1 {
class ExactWeightSum {
  static constexpr int L = 35;
  std::array<std::uint64_t, L> limb_{};
  void add_at(int k, std::uint64_t v) {
    for (; v && k < L; ++k) {
      const std::uint64_t s = limb_[std::size_t(k)] + v;
      v = s < v ? 1 : 0;
      limb_[std::size_t(k)] = s;
    }
    if (v) throw std::runtime_error("exact weight sum overflow");
  }

 public:
  void add(double x) {
    if (!(std::isfinite(x) && x >= 0)) throw std::runtime_error("exact weight sum: negative or nonfinite term");
    if (x == 0) return;
    const std::uint64_t bits = std::bit_cast<std::uint64_t>(x);
    const int be = int(bits >> 52) & 0x7ff;
    std::uint64_t m = bits & ((std::uint64_t(1) << 52) - 1);
    int shift = 0;
    if (be) { m |= std::uint64_t(1) << 52; shift = be - 1; }
    const int k = shift / 64, r = shift % 64;
    add_at(k, m << r);
    if (r) add_at(k + 1, m >> (64 - r));
  }
  void add(const ExactWeightSum& o) {
    std::uint64_t carry = 0;
    for (int k = 0; k < L; ++k) {
      const std::uint64_t a = limb_[std::size_t(k)], b = o.limb_[std::size_t(k)];
      const std::uint64_t s = a + b;
      const std::uint64_t c1 = s < a;
      const std::uint64_t t = s + carry;
      const std::uint64_t c2 = t < s;
      limb_[std::size_t(k)] = t;
      carry = c1 + c2;
    }
    if (carry) throw std::runtime_error("exact weight sum overflow");
  }
  // Exact value of a native exact count (its positive carriers sum to it exactly).
  void add(const m::ArenaCount& c) { for (double w : c.positive_carriers()) add(w); }
  bool operator==(const ExactWeightSum& o) const noexcept { return limb_ == o.limb_; }
  bool zero() const noexcept {
    for (auto x : limb_) if (x) return false;
    return true;
  }
  // Correctly rounded (to nearest, ties to even) binary64 value.
  double value() const {
    int top = L - 1;
    while (top >= 0 && !limb_[std::size_t(top)]) --top;
    if (top < 0) return 0.;
    const int H = 64 * top + 63 - std::countl_zero(limb_[std::size_t(top)]);   // highest set bit (units 2^-1074)
    const auto bit = [&](int i) { return i < 0 ? 0u : unsigned((limb_[std::size_t(i / 64)] >> (i % 64)) & 1u); };
    const auto field = [&](int lo, int n) {   // n <= 64 bits starting at lo (lo may be negative: zeros)
      std::uint64_t v = 0;
      for (int i = n - 1; i >= 0; --i) v = (v << 1) | bit(lo + i);
      return v;
    };
    if (H <= 52) return std::ldexp(double(field(0, H + 1)), -1074);
    std::uint64_t M = field(H - 52, 53);
    const unsigned guard = bit(H - 53);
    const auto any_below = [&](int n) {   // any set bit with index < n
      if (n <= 0) return false;
      const int full = n / 64, r = n % 64;
      for (int k = 0; k < full; ++k) if (limb_[std::size_t(k)]) return true;
      return r ? (limb_[std::size_t(full)] & ((std::uint64_t(1) << r) - 1)) != 0 : false;
    };
    const bool sticky = any_below(H - 53);
    int e = H - 52;
    if (guard && (sticky || (M & 1u))) {
      ++M;
      if (M == (std::uint64_t(1) << 53)) { M >>= 1; ++e; }
    }
    return std::ldexp(double(M), e - 1074);
  }
};
}  // namespace raikou::native_water::v7_1
