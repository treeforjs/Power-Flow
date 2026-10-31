#pragma once
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

// Binary-encounter-Bethe (BEB) and relativistic BEB (RBEB) electron-impact
// ionization of one molecular orbital, and the singly differential cross
// sections (SDCS) whose integrals over the ejected-electron energy reproduce
// them exactly.
//
// Notation (per orbital): B binding energy, U mean orbital kinetic energy,
// N occupation, T incident kinetic energy in the target rest frame,
// t = T/B, u = U/B, w = W/B with W the kinetic energy of the slower (ejected)
// outgoing electron, 0 <= w <= (t-1)/2.
//
// BEB, Kim & Rudd, Phys. Rev. A 50, 3954 (1994), as used by the NIST
// electron-impact ionization database and implemented by the dataset
// (collision_data.h, model "nist_BEB", optional orbital_weights):
//   sigma = S/(t+u+1) [ Q ln t/2 (1 - 1/t^2) + (2-Q)(1 - 1/t - ln t/(t+1)) ],
//   S = 4 pi a0^2 N (R/B)^2.
// beb_weighted_sum_m2 repeats the dataset arithmetic operation by operation,
// so its binary64 values are identical to the dataset implementation.
//
// RBEB, Kim, Santos & Parente, Phys. Rev. A 62, 052710 (2000) (Q = 1):
//   sigma = 4 pi a0^2 alpha^4 N / ((bt^2 + bu^2 + bb^2) 2b') *
//           { 1/2 [ ln(bt^2/(1-bt^2)) - bt^2 - ln(2b') ] (1 - 1/t^2)
//             + 1 - 1/t - ln t/(t+1) (1+2t')/(1+t'/2)^2
//             + b'^2/(1+t'/2)^2 (t-1)/2 },
//   t' = T/mc^2, b' = B/mc^2, u' = U/mc^2, bx^2 = 1 - 1/(1+x')^2.
// alpha^4 is written (2R/mc^2)^2 (alpha^2 mc^2 = 2R), so the nonrelativistic
// limit (t', b', u' -> 0) is the BEB formula above with Q = 1 and the same a0
// and R. bx^2 = x'(x'+2)/(1+x')^2 and bt^2/(1-bt^2) = t'(t'+2) are used to
// avoid cancellation.
//
// SDCS. Both totals are P * G(t), G(t) = integral_0^{(t-1)/2} g(w) dw, with
//   g(w) = -phi/(t+1) [1/(w+1) + 1/(t-w)] + [1/(w+1)^2 + 1/(t-w)^2] + c_b
//          + L [1/(w+1)^3 + 1/(t-w)^3],
//   dsigma/dW = (P/B) g(W/B).
// RBEB: P as above, phi = (1+2t')/(1+t'/2)^2, c_b = b'^2/(1+t'/2)^2,
//       L = ln(bt^2/(1-bt^2)) - bt^2 - ln(2b')  (KSP 2000 RBED SDCS with the
//       BEB dipole substitution and Mott interference factor as in RBEB).
// BEB (Q = 1): P = S/(t+u+1), phi = 1, c_b = 0, L = ln t (Kim & Rudd SDCS with
//       the BEB dipole substitution). The dipole term is written symmetric under
//       the exchange w <-> t-1-w, which makes its integral over the ejected
//       half-range equal to the BEB term Q ln t/2 (1 - 1/t^2).
// The analytic integral G(w) (cumulative distribution) is
//   G(w) = -phi/(t+1) [log1p(w) - log1p(-w/t)] + w/(w+1) + w/(t(t-w)) + c_b w
//          + L/2 [w(w+2)/(w+1)^2 + w(2t-w)/(t^2 (t-w)^2)].
namespace raikou::native_water::beb_v2 {

struct Orbital {
  std::string label;
  double B_eV = 0, U_eV = 0, N = 0, Q = 0;
};
struct Constants {
  double bohr_radius_m = 0, rydberg_eV = 0, rest_energy_eV = 0;
};
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(what);
}

// Dataset arithmetic (collision_data.h, orbital-weighted nist_BEB), one
// operation at a time in the same order. weights may be empty (all one).
inline double beb_weighted_sum_m2(const std::vector<Orbital>& orbitals, const std::vector<double>& weights,
                                  const Constants& k, double energy) {
  double sum = 0;
  const double a = k.bohr_radius_m, R = k.rydberg_eV;
  std::size_t index = 0;
  for (const auto& orbital : orbitals) {
    const double weight = weights.empty() ? 1. : weights[index];
    ++index;
    const double B = orbital.B_eV;
    if (energy <= B) continue;
    const double t = energy / B, u = orbital.U_eV / B;
    const double Q = orbital.Q, N = orbital.N;
    const double S = 4 * std::numbers::pi * a * a * N * (R / B) * (R / B), l = std::log(t);
    sum += weight * S / (t + u + 1) * (Q * l / 2 * (1 - 1 / (t * t)) + (2 - Q) * (1 - 1 / t - l / (t + 1)));
  }
  return sum;
}
// One unweighted dataset-arithmetic orbital term.
inline double beb_orbital_m2(const Orbital& o, const Constants& k, double energy) {
  return beb_weighted_sum_m2({o}, {}, k, energy);
}

// Reduced SDCS shape of one orbital at one incident energy.
struct Shape {
  long double t = 0, phi = 1, cb = 0, L = 0, P = 0; // P in m^2
  long double B = 0;                                 // eV
  bool relativistic = false;
  long double w_max() const { return (t - 1) / 2; }
  long double g(long double w) const {
    const long double a = 1 / (w + 1), b = 1 / (t - w);
    return -phi / (t + 1) * (a + b) + (a * a + b * b) + cb + L * (a * a * a + b * b * b);
  }
  // Integral of g from 0 to w, analytic, without cancellation for small w.
  long double G(long double w) const {
    const long double tw = t - w;
    return -phi / (t + 1) * (std::log1p(w) - std::log1p(-w / t)) + w / (w + 1) + w / (t * tw) + cb * w +
           L / 2 * (w * (w + 2) / ((w + 1) * (w + 1)) + w * (2 * t - w) / (t * t * tw * tw));
  }
  long double total_m2() const { return P * G(w_max()); }
  // dsigma/dW in m^2/eV.
  long double sdcs_m2_per_eV(long double W) const { return P / B * g(W / B); }
};

inline Shape beb_shape(const Orbital& o, const Constants& k, long double T) {
  require(o.Q == 1, "the BEB SDCS is implemented for Q = 1 only");
  require(std::isfinite(static_cast<double>(T)) && T > o.B_eV, "BEB shape needs an incident energy above the binding energy");
  Shape s;
  s.B = o.B_eV;s.t = T / s.B;s.phi = 1;s.cb = 0;s.L = std::log(s.t);
  const long double a = k.bohr_radius_m, R = k.rydberg_eV, u = o.U_eV / s.B;
  const long double S = 4 * std::numbers::pi_v<long double> * a * a * o.N * (R / s.B) * (R / s.B);
  s.P = S / (s.t + u + 1);
  return s;
}
inline long double beta2(long double x) { return x * (x + 2) / ((1 + x) * (1 + x)); }
inline Shape rbeb_shape(const Orbital& o, const Constants& k, long double T) {
  require(o.Q == 1, "RBEB (Kim, Santos & Parente 2000) is defined with Q = 1");
  require(std::isfinite(static_cast<double>(T)) && T > o.B_eV, "RBEB shape needs an incident energy above the binding energy");
  require(k.rest_energy_eV > 0 && k.bohr_radius_m > 0 && k.rydberg_eV > 0, "RBEB constants");
  Shape s;s.relativistic = true;
  const long double mc2 = k.rest_energy_eV, B = o.B_eV, U = o.U_eV;
  s.B = B;s.t = T / B;
  const long double tp = T / mc2, bp = B / mc2, up = U / mc2;
  const long double bt2 = beta2(tp), bb2 = beta2(bp), bu2 = beta2(up);
  const long double alpha2 = 2 * static_cast<long double>(k.rydberg_eV) / mc2;
  const long double a = k.bohr_radius_m;
  s.P = 4 * std::numbers::pi_v<long double> * a * a * alpha2 * alpha2 * o.N / ((bt2 + bu2 + bb2) * 2 * bp);
  s.L = std::log(tp * (tp + 2)) - bt2 - std::log(2 * bp);
  const long double h = 1 + tp / 2;
  s.phi = (1 + 2 * tp) / (h * h);
  s.cb = bp * bp / (h * h);
  return s;
}
// RBEB total of one orbital written as the closed form of KSP 2000 Eq. (RBEB),
// independently of the SDCS integral.
inline long double rbeb_orbital_m2(const Orbital& o, const Constants& k, long double T) {
  if (!(T > o.B_eV)) return 0;
  const auto s = rbeb_shape(o, k, T);
  const long double t = s.t;
  return s.P * (s.L / 2 * (1 - 1 / (t * t)) + 1 - 1 / t - s.phi * std::log(t) / (t + 1) + s.cb * (t - 1) / 2);
}
// Nonrelativistic BEB (Q = 1) in extended precision, same closed form.
inline long double beb_orbital_ext_m2(const Orbital& o, const Constants& k, long double T) {
  if (!(T > o.B_eV)) return 0;
  const auto s = beb_shape(o, k, T);
  const long double t = s.t;
  return s.P * (s.L / 2 * (1 - 1 / (t * t)) + 1 - 1 / t - std::log(t) / (t + 1));
}

// Inverse of the normalized cumulative distribution by bisection on w. The
// result is the smallest representable w (to bisection resolution) with
// G(w) >= u G(w_max); u must lie in (0,1).
inline long double sample_reduced_w(const Shape& s, long double u) {
  require(u > 0 && u < 1, "SDCS variate must lie in (0,1)");
  const long double wm = s.w_max(), target = u * s.G(wm);
  long double lo = 0, hi = wm;
  for (int i = 0; i < 512; ++i) {
    const long double mid = lo + (hi - lo) / 2;
    if (!(mid > lo && mid < hi)) break;
    if (s.G(mid) < target) lo = mid;
    else hi = mid;
  }
  return hi;
}

// Orbital selection with probability proportional to weight_i sigma_i.
// sigma values are passed in, so the band's own rate terms are reused.
inline std::size_t select_orbital(const std::vector<long double>& weighted_terms, long double u) {
  require(u > 0 && u < 1, "orbital variate must lie in (0,1)");
  long double total = 0;
  for (auto x : weighted_terms) { require(std::isfinite(static_cast<double>(x)) && x >= 0, "orbital term"); total += x; }
  require(total > 0, "no positive orbital term");
  const long double choice = u * total;
  long double cumulative = 0;
  std::size_t last = weighted_terms.size();
  for (std::size_t i = 0; i < weighted_terms.size(); ++i) {
    if (!(weighted_terms[i] > 0)) continue;
    last = i;
    cumulative += weighted_terms[i];
    if (choice < cumulative) return i;
  }
  return last; // rounding at the upper end of the cumulative sum
}
} // namespace raikou::native_water::beb_v2
