#pragma once
// Nominal, conditional 304 thermal scenario. Not a compressible EOS or mechanics.
// Source and approximation contract: ../data/scenario.json and ../README.md.
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ss304 {
enum class Qualification { evaluated_solid_1500 = 0, nominal_synthesis_1800 = 1 };
enum class DensityConvention { fixed_reference_volume_cp_as_cv_approximation = 0 };
enum class Phase { solid = 0, mushy = 1, liquid = 2 };
struct Values {
  double temperature_k, h_j_kg, u_model_j_kg;
  double cp_sensible_j_kg_k, apparent_heat_capacity_j_kg_k;
  double k_w_m_k, resistivity_ohm_m, sigma_s_m, liquid_fraction;
  double energy_j_m3_reference, capacity_j_m3_reference_k;
  int phase, synthesized_high_temperature;
};
struct DepositResult { Values state; double signed_numerical_energy_residual_j; };

class Model {
 public:
  static constexpr double t_reference = 300.0;
  static constexpr double t_solidus = 1672.0;
  static constexpr double t_liquidus = 1727.0;
  static constexpr double cp_liquid = 800.0;
  static constexpr double latent_fusion = 160.0 * t_liquidus;

  // No implicit scenario or density selection. rho_ref is mass/reference volume;
  // density is NOT evolved as a function of T or liquid fraction.
  Model(Qualification q, DensityConvention convention, double rho_ref)
      : qualification_(q), rho_ref_(rho_ref) {
    if (q != Qualification::evaluated_solid_1500 && q != Qualification::nominal_synthesis_1800)
      throw std::invalid_argument("unknown qualification");
    if (convention != DensityConvention::fixed_reference_volume_cp_as_cv_approximation)
      throw std::invalid_argument("unsupported density/work convention");
    if (!std::isfinite(rho_ref) || rho_ref <= 0)
      throw std::invalid_argument("rho_ref must be finite and positive");
    // Guarantee finite forward endpoints for the declared scalar domain.
    if (!std::isfinite(rho_ref * h_unchecked(max_temperature_k())) ||
        !std::isfinite(rho_ref * (cp_liquid + latent_fusion / (t_liquidus - t_solidus))) ||
        rho_ref * cp_solid(t_reference) == 0)
      throw std::invalid_argument("reference density exceeds representable energy/capacity range");
  }
  double max_temperature_k() const {
    return qualification_ == Qualification::evaluated_solid_1500 ? 1500.0 : 1800.0;
  }
  double reference_density_kg_m3() const { return rho_ref_; }
  double max_enthalpy_j_kg() const { return h_unchecked(max_temperature_k()); }

  Values evaluate(double t) const {
    temperature_guard(t);
    double f = 0, cp = cp_solid(t), c = cp;
    Phase phase = Phase::solid;
    // Right-sided derivative at the two phase boundaries. h and f are continuous.
    if (t >= t_liquidus) {
      f = 1; cp = cp_liquid; c = cp; phase = Phase::liquid;
    } else if (t >= t_solidus) {
      f = (t - t_solidus) / (t_liquidus - t_solidus);
      cp = cp_solid(t_solidus) + (cp_liquid - cp_solid(t_solidus)) * f;
      c = cp + latent_fusion / (t_liquidus - t_solidus);
      phase = Phase::mushy;
    }
    const auto transport = transport_unchecked(t);
    const double h = h_unchecked(t);
    return {t, h, h, cp, c, transport[1], transport[0], 1.0 / transport[0], f,
            rho_ref_ * h, rho_ref_ * c, static_cast<int>(phase), t > 1500.0};
  }

  // Monotone bounded inversion, including exactly the rounded forward endpoints.
  // Bisection terminates at adjacent floating-point temperatures, not a physical
  // tolerance. Choose the adjacent state with the smaller forward-energy error.
  double temperature_from_enthalpy_k(double h) const {
    if (!std::isfinite(h) || h < 0 || h > max_enthalpy_j_kg())
      throw std::out_of_range("specific energy outside qualified range; no clamp/extrapolation");
    if (h == 0) return t_reference;
    if (h == max_enthalpy_j_kg()) return max_temperature_k();
    double lo = t_reference, hi = max_temperature_k();
    for (int i = 0; i < 128; ++i) {
      const double mid = lo + (hi - lo) * 0.5;
      if (mid == lo || mid == hi) break;
      const double hm = h_unchecked(mid);
      if (hm == h) return mid;
      if (hm < h) lo = mid; else hi = mid;
    }
    return std::abs(h_unchecked(lo) - h) <= std::abs(h_unchecked(hi) - h) ? lo : hi;
  }
  double temperature_from_u_model_k(double u) const { return temperature_from_enthalpy_k(u); }

  double temperature_from_reference_energy_k(double e) const {
    const double max_e = rho_ref_ * max_enthalpy_j_kg();
    if (!std::isfinite(e) || e < 0 || e > max_e)
      throw std::out_of_range("reference-volume energy outside qualified range");
    if (e == 0) return t_reference;
    if (e == max_e) return max_temperature_k();
    return temperature_from_enthalpy_k(e / rho_ref_);
  }

  // Fixed-reference-volume transaction: accepted state is never modified.
  // specific_energy + signed deposited J / accepted mass is the target state.
  DepositResult after_deposited_energy(double specific_energy, double mass_kg, double energy_j) const {
    if (!std::isfinite(mass_kg) || mass_kg <= 0 || !std::isfinite(energy_j))
      throw std::invalid_argument("mass must be positive; deposited energy must be finite");
    (void)temperature_from_u_model_k(specific_energy);
    const long double remaining = static_cast<long double>(mass_kg) *
      (static_cast<long double>(max_enthalpy_j_kg()) - specific_energy);
    const long double available = static_cast<long double>(mass_kg) * specific_energy;
    if (static_cast<long double>(energy_j) > remaining || -static_cast<long double>(energy_j) > available)
      throw std::out_of_range("deposited energy leaves qualified caloric range");
    const double target = specific_energy + energy_j / mass_kg;
    const Values state = evaluate(temperature_from_u_model_k(target));
    return {state, std::fma(mass_kg, state.u_model_j_kg - specific_energy, -energy_j)};
  }

 private:
  Qualification qualification_;
  double rho_ref_;
  void temperature_guard(double t) const {
    if (!std::isfinite(t) || t < t_reference || t > max_temperature_k())
      throw std::out_of_range("temperature outside qualified range; no clamp/extrapolation");
  }
  static double cp_solid(double t) { return 472.0 + 0.136 * t - 2820000.0 / (t * t); }
  static double h_solid(double t) {
    // Exact analytic integral of Mills Eq.10, re-zeroed at 300 K; factorization
    // retains the small energy increment near the reference temperature.
    return (t - t_reference) *
           (472.0 + 0.068 * (t + t_reference) - 2820000.0 / (t * t_reference));
  }
  static double h_mushy_sensible(double dt) {
    const double cs = cp_solid(t_solidus);
    return h_solid(t_solidus) + cs * dt +
           0.5 * (cp_liquid - cs) * dt * dt / (t_liquidus - t_solidus);
  }
  static double h_unchecked(double t) {
    const double width = t_liquidus - t_solidus;
    if (t < t_solidus) return h_solid(t);
    if (t < t_liquidus) {
      const double dt = t - t_solidus;
      return h_mushy_sensible(dt) + latent_fusion * (dt / width);
    }
    return h_mushy_sensible(width) + latent_fusion + cp_liquid * (t - t_liquidus);
  }
  static std::array<double, 2> transport_unchecked(double t) {
    // Ho & Chu (1977), Table2 printed p15/PDF22. Units: K, 1e-8 ohm m, W/m/K.
    // Published synthesized values; above1500K explicitly extrapolated by source.
    // PWL in resistivity and k (sigma is reciprocal AFTER interpolation).
    constexpr std::array<double, 18> ts = {
      300,350,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1672,1727,1800};
    constexpr std::array<double, 18> rs = {
      71.9,76.0,79.8,86.8,93.3,99.2,104.3,108.6,112.5,115.8,118.7,121.4,124.0,126.4,128.7,130.4,150.0,151.7};
    constexpr std::array<double, 18> ks = {
      14.89,15.79,16.61,18.28,19.77,21.21,22.59,23.99,25.33,26.58,27.81,29.18,30.34,31.55,32.70,33.53,28.15,28.99};
    const auto upper = std::upper_bound(ts.begin(), ts.end(), t);
    const std::size_t i = std::min<std::size_t>(static_cast<std::size_t>(upper - ts.begin() - 1), ts.size() - 2);
    const double f = (t - ts[i]) / (ts[i+1] - ts[i]);
    return {(rs[i] + f * (rs[i+1] - rs[i])) * 1e-8, ks[i] + f * (ks[i+1] - ks[i])};
  }
};
} // namespace ss304
