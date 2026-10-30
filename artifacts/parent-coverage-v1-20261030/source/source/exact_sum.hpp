#pragma once
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Exact, order-independent sum of finite nonnegative binary64 values.
// The value is held as one fixed-point integer in units of 2^-1074 (the
// smallest subnormal), so every binary64 addend is represented without
// rounding and the stored sum never depends on addition order. Rounding to
// binary64 happens only when a value is reported (round half to even).
namespace raikou::native_water::coverage_v1 {
class ExactSum {
public:
  // 2098 bits hold every finite binary64; 36 limbs leave 206 bits of carry
  // headroom, i.e. more than 2^200 maximal addends before overflow.
  static constexpr int bias = 1074;
  static constexpr std::size_t limbs = 36;
  static constexpr std::size_t maximum_hex_digits = limbs * 16;

  ExactSum() = default;
  explicit ExactSum(double x) { add(x); }

  void add(double x) {
    if (!(std::isfinite(x) && x >= 0))
      throw std::invalid_argument("exact sum accepts only finite nonnegative binary64 addends");
    std::uint64_t bits;
    std::memcpy(&bits, &x, sizeof bits);
    const std::uint64_t exponent = (bits >> 52) & 0x7ff, fraction = bits & ((std::uint64_t(1) << 52) - 1);
    if (exponent == 0 && fraction == 0) return; // +0 and -0 add nothing
    const std::uint64_t mantissa = exponent ? (fraction | (std::uint64_t(1) << 52)) : fraction;
    const std::size_t offset = exponent ? static_cast<std::size_t>(exponent - 1) : 0;
    add_shifted(mantissa, offset);
  }
  void add(const ExactSum& other) {
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < limbs; ++i) carry = add_with_carry(limb_[i], other.limb_[i], carry);
    if (carry) throw std::overflow_error("exact sum capacity exhausted");
  }
  bool zero() const noexcept {
    for (auto v : limb_) if (v) return false;
    return true;
  }
  bool operator==(const ExactSum& o) const noexcept { return limb_ == o.limb_; }
  bool operator!=(const ExactSum& o) const noexcept { return !(*this == o); }
  bool operator<(const ExactSum& o) const noexcept {
    for (std::size_t i = limbs; i-- > 0;)
      if (limb_[i] != o.limb_[i]) return limb_[i] < o.limb_[i];
    return false;
  }

  // Correctly rounded (nearest, ties to even) binary64 value; +inf when the
  // exact sum exceeds the binary64 range.
  double value() const noexcept {
    const long top = highest_bit();
    if (top < 0) return 0;
    if (top <= 52) return std::ldexp(static_cast<double>(limb_[0]), -bias);
    std::uint64_t mantissa = bits(static_cast<std::size_t>(top - 52), 53);
    const bool round = bit(static_cast<std::size_t>(top - 53));
    const bool sticky = top >= 54 && any_below(static_cast<std::size_t>(top - 53));
    long exponent = top - 52 - bias;
    if (round && (sticky || (mantissa & 1))) {
      ++mantissa;
      if (mantissa == (std::uint64_t(1) << 53)) { mantissa >>= 1; ++exponent; }
    }
    if (exponent > 1023 - 52) return std::numeric_limits<double>::infinity();
    return std::ldexp(static_cast<double>(mantissa), static_cast<int>(exponent));
  }

  // Canonical exact form: value = M * 2^E with M odd (or M = 0, E = 0).
  // M is lowercase hexadecimal without leading zeros.
  std::string mantissa_hex() const {
    const long low = lowest_bit();
    if (low < 0) return "0";
    const long top = highest_bit();
    const std::size_t width = static_cast<std::size_t>(top - low + 1);
    std::string out;
    const std::size_t digits = (width + 3) / 4;
    out.reserve(digits);
    for (std::size_t d = digits; d-- > 0;) {
      const std::size_t start = static_cast<std::size_t>(low) + 4 * d;
      const std::size_t count = std::min<std::size_t>(4, static_cast<std::size_t>(top) + 1 - start);
      const unsigned nibble = static_cast<unsigned>(bits(start, count));
      out.push_back("0123456789abcdef"[nibble]);
    }
    return out;
  }
  long exponent2() const noexcept {
    const long low = lowest_bit();
    return low < 0 ? 0 : low - bias;
  }
  static ExactSum from_canonical(const std::string& hex, long exponent) {
    ExactSum out;
    if (hex == "0") {
      if (exponent != 0) throw std::invalid_argument("noncanonical exact zero");
      return out;
    }
    if (hex.empty() || hex.size() > maximum_hex_digits || hex[0] == '0')
      throw std::invalid_argument("noncanonical exact mantissa");
    const long shift = exponent + bias;
    if (shift < 0 || shift > static_cast<long>(limbs * 64))
      throw std::invalid_argument("exact exponent outside capacity");
    for (std::size_t i = 0; i < hex.size(); ++i) {
      const char c = hex[hex.size() - 1 - i];
      unsigned v;
      if (c >= '0' && c <= '9') v = static_cast<unsigned>(c - '0');
      else if (c >= 'a' && c <= 'f') v = static_cast<unsigned>(c - 'a' + 10);
      else throw std::invalid_argument("exact mantissa is not lowercase hexadecimal");
      if (i == 0 && !(v & 1)) throw std::invalid_argument("noncanonical even exact mantissa");
      for (unsigned b = 0; b < 4; ++b)
        if (v >> b & 1) {
          const std::size_t position = static_cast<std::size_t>(shift) + 4 * i + b;
          if (position >= limbs * 64 - 64) throw std::invalid_argument("exact value exceeds capacity");
          out.limb_[position / 64] |= std::uint64_t(1) << (position % 64);
        }
    }
    return out;
  }

private:
  std::array<std::uint64_t, limbs> limb_{};

  // limb += a + carry_in (carry_in <= 1); returns the carry out (0 or 1).
  static std::uint64_t add_with_carry(std::uint64_t& limb, std::uint64_t a, std::uint64_t carry_in) noexcept {
    const std::uint64_t s1 = limb + a;
    const std::uint64_t c1 = s1 < limb;
    const std::uint64_t s2 = s1 + carry_in;
    const std::uint64_t c2 = s2 < s1;
    limb = s2;
    return c1 + c2;
  }
  void add_shifted(std::uint64_t mantissa, std::size_t offset) {
    const std::size_t q = offset / 64, r = offset % 64;
    const std::uint64_t lo = mantissa << r;
    const std::uint64_t hi = r ? mantissa >> (64 - r) : 0;
    std::uint64_t carry = add_with_carry(limb_[q], lo, 0);
    std::size_t i = q + 1;
    if (i < limbs) { carry = add_with_carry(limb_[i], hi, carry); ++i; }
    else if (hi || carry) throw std::overflow_error("exact sum capacity exhausted");
    for (; carry && i < limbs; ++i) carry = add_with_carry(limb_[i], 0, carry);
    if (carry) throw std::overflow_error("exact sum capacity exhausted");
  }
  long highest_bit() const noexcept {
    for (std::size_t i = limbs; i-- > 0;)
      if (limb_[i]) return static_cast<long>(64 * i) + 63 - std::countl_zero(limb_[i]);
    return -1;
  }
  long lowest_bit() const noexcept {
    for (std::size_t i = 0; i < limbs; ++i)
      if (limb_[i]) return static_cast<long>(64 * i) + std::countr_zero(limb_[i]);
    return -1;
  }
  bool bit(std::size_t position) const noexcept {
    return position < limbs * 64 && (limb_[position / 64] >> (position % 64) & 1);
  }
  // count <= 64 bits starting at position, least significant first.
  std::uint64_t bits(std::size_t position, std::size_t count) const noexcept {
    std::uint64_t out = 0;
    for (std::size_t k = 0; k < count; ++k)
      if (bit(position + k)) out |= std::uint64_t(1) << k;
    return out;
  }
  bool any_below(std::size_t position) const noexcept {
    const std::size_t full = position / 64;
    for (std::size_t i = 0; i < full; ++i) if (limb_[i]) return true;
    const std::size_t r = position % 64;
    return r && (limb_[full] & ((std::uint64_t(1) << r) - 1));
  }
};
} // namespace raikou::native_water::coverage_v1
