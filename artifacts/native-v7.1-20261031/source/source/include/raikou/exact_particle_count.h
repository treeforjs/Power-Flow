#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <openssl/bn.h>
#include <stdexcept>
#include <utility>
#include <vector>

namespace raikou::mfem_backend {
// Exact represented number, in binary64's existing2^-1074 units. This is a
// numerical encoding of admitted weights, not a physical source quantum/floor.
// Values are immutable and bounded to2098 bits (the finite binary64 range).
class ExactParticleCount {
public:
  using Integer = std::unique_ptr<BIGNUM, decltype(&BN_free)>;
  static constexpr int maximum_bits = 2098;
  static constexpr std::size_t maximum_carriers = 42;
  ExactParticleCount() noexcept = default;
  // Count2 interoperability: retain its normalized hi/lo input exactly.
  ExactParticleCount(double hi, double lo = 0) {
    need(std::isfinite(hi) && std::isfinite(lo) && std::isfinite(hi + lo) &&
             hi >= 0,
         "nonfinite/negative Count2 input");
    if (hi == 0)
      need(lo == 0, "Count2 zero normalization");
    else
      need(std::abs(lo) <= .5 * (std::nextafter(hi, INFINITY) - hi),
           "Count2 normalization");
    auto n = integer(hi), low = integer(lo);
    need(BN_add(n.get(), n.get(), low.get()), "Count2 import");
    *this = adopt(std::move(n));
  }
  template <class Count2>
  static ExactParticleCount from_two_word(const Count2 &n) {
    return {n.hi, n.lo};
  }
  template <class Count2> Count2 require_two_word() const {
    const auto w = require_two_limb_values();
    return {w[0], w[1]};
  }
  bool zero() const noexcept { return !value_; }
  bool operator==(const ExactParticleCount &other) const noexcept {
    if (!value_ || !other.value_)
      return value_ == other.value_;
    return BN_cmp(value_->units.get(), other.value_->units.get()) == 0;
  }
  bool operator!=(const ExactParticleCount &other) const noexcept {
    return !(*this == other);
  }
  Integer copy_units() const {
    auto n = allocate();
    if (value_)
      need(BN_copy(n.get(), value_->units.get()), "count copy");
    return n;
  }
  ExactParticleCount debit(double amount) const {
    need(std::isfinite(amount) && amount > 0, "positive finite count debit");
    auto n = copy_units(), w = integer(amount);
    need(BN_cmp(n.get(), w.get()) >= 0, "insufficient exact inventory");
    need(BN_sub(n.get(), n.get(), w.get()), "exact debit");
    return adopt(std::move(n));
  }
  ExactParticleCount credited(double amount) const {
    need(std::isfinite(amount) && amount >= 0,
         "nonnegative finite count credit");
    auto n = copy_units(), w = integer(amount);
    need(BN_add(n.get(), n.get(), w.get()), "exact credit");
    return adopt(std::move(n));
  }
  double rounded_for_rate() const { return rounded(copy_units().get()); }
  std::array<double, 2> require_two_limb_values() const {
    const double hi = rounded_for_rate();
    auto remainder = copy_units(), high = integer(hi);
    need(BN_sub(remainder.get(), remainder.get(), high.get()),
         "Count2 export subtraction");
    bool negative = BN_is_negative(remainder.get());
    BN_set_negative(remainder.get(), 0);
    double lo = rounded(remainder.get());
    auto encoded = integer(lo);
    need(BN_cmp(encoded.get(), remainder.get()) == 0,
         "exact gas count cannot be represented in two limbs; conversion "
         "refused");
    if (negative)
      lo = -lo;
    return {hi, lo};
  }
  std::vector<double> positive_carriers() const {
    auto remaining = copy_units();
    std::vector<double> out;
    while (!BN_is_zero(remaining.get())) {
      need(out.size() < maximum_carriers,
           "bounded positive-carrier representation");
      const int shift = std::max(0, BN_num_bits(remaining.get()) - 53);
      auto top = allocate();
      need(BN_rshift(top.get(), remaining.get(), shift), "carrier shift");
      const double weight =
          std::ldexp(double(BN_get_word(top.get())), shift - 1074);
      need(std::isfinite(weight) && weight > 0,
           "native positive carrier unrepresentable");
      auto word = integer(weight);
      need(BN_cmp(word.get(), remaining.get()) <= 0 &&
               BN_sub(remaining.get(), remaining.get(), word.get()),
           "exact positive carrier decomposition");
      out.push_back(weight);
    }
    return out;
  }
  static ExactParticleCount
  from_positive_carriers(const std::vector<double> &weights) {
    need(weights.size() <= maximum_carriers, "bounded count input limbs");
    auto n = allocate();
    for (double weight : weights) {
      need(std::isfinite(weight) && weight > 0,
           "positive finite count input limb");
      auto w = integer(weight);
      need(BN_add(n.get(), n.get(), w.get()), "count input limb sum");
    }
    return adopt(std::move(n));
  }

private:
  struct Payload {
    Integer units;
    explicit Payload(Integer n) : units(std::move(n)) {};
  };
  std::shared_ptr<const Payload> value_;
  static void need(bool ok, const char *message) {
    if (!ok)
      throw std::invalid_argument(message);
  }
  static Integer allocate() {
    Integer n(BN_new(), BN_free);
    if (!n)
      throw std::bad_alloc();
    return n;
  }
  static Integer integer(double x) {
    static_assert(sizeof(double) == 8 &&
                  std::numeric_limits<double>::is_iec559 &&
                  std::numeric_limits<double>::digits == 53);
    static_assert(sizeof(BN_ULONG) >= 8);
    need(std::isfinite(x), "finite binary64 count operand");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, 8);
    const unsigned exponent = (bits >> 52) & 2047;
    std::uint64_t mantissa = bits & ((1ULL << 52) - 1);
    if (exponent)
      mantissa |= 1ULL << 52;
    auto n = allocate();
    need(BN_set_word(n.get(), mantissa), "binary64 integer import");
    if (exponent)
      need(BN_lshift(n.get(), n.get(), int(exponent) - 1),
           "binary64 integer shift");
    BN_set_negative(n.get(), (bits >> 63) && mantissa);
    return n;
  }
  static ExactParticleCount adopt(Integer n) {
    need(!BN_is_negative(n.get()) && BN_num_bits(n.get()) <= maximum_bits,
         "exact count outside bounded nonnegative binary64 range");
    ExactParticleCount out;
    if (!BN_is_zero(n.get()))
      out.value_ = std::make_shared<const Payload>(std::move(n));
    return out;
  }
  static double rounded(const BIGNUM *n) {
    need(!BN_is_negative(n) && BN_num_bits(n) <= maximum_bits,
         "bounded rate-count operand");
    if (BN_is_zero(n))
      return 0;
    const int shift = std::max(0, BN_num_bits(n) - 53);
    auto top = allocate();
    need(BN_rshift(top.get(), n, shift), "rate count shift");
    BN_ULONG word = BN_get_word(top.get());
    if (shift) {
      auto remainder = allocate(), half = allocate();
      need(BN_copy(remainder.get(), n) &&
               BN_mask_bits(remainder.get(), shift) &&
               BN_set_bit(half.get(), shift - 1),
           "rate count rounding operands");
      const int compare = BN_cmp(remainder.get(), half.get());
      if (compare > 0 || (compare == 0 && (word & 1)))
        ++word;
    }
    const double result = std::ldexp(double(word), shift - 1074);
    need(std::isfinite(result),
         "rounded diagnostic count overflow; use exact positive carriers");
    return result;
  }
};
} // namespace raikou::mfem_backend
