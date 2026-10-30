#include "include/ss304_thermal.hpp"
#include <exception>
#include <string>

namespace {
thread_local std::string last_error;
ss304::Model model(int mode, double density) {
  return ss304::Model(static_cast<ss304::Qualification>(mode),
    ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation, density);
}
template<class F> int guarded(F f) {
  try { f(); last_error.clear(); return 0; }
  catch (const std::exception& e) { last_error = e.what(); return 1; }
  catch (...) { last_error = "unknown property error"; return 2; }
}
}
extern "C" {
const char* ss304_error() { return last_error.c_str(); }
int ss304_evaluate(int mode, double density, double t, ss304::Values* out) {
  return guarded([&] { if (!out) throw std::invalid_argument("null output");
    const auto v = model(mode,density).evaluate(t); *out = v; });
}
int ss304_inverse(int mode, double density, double energy, int reference_volume, double* out) {
  return guarded([&] { if (!out) throw std::invalid_argument("null output");
    if (reference_volume != 0 && reference_volume != 1) throw std::invalid_argument("unknown energy measure");
    const auto m = model(mode,density);
    const double t = reference_volume ? m.temperature_from_reference_energy_k(energy) : m.temperature_from_u_model_k(energy);
    *out = t; });
}
int ss304_deposit(int mode, double density, double u, double mass, double joules, ss304::Values* out, double* residual) {
  return guarded([&] { if (!out || !residual) throw std::invalid_argument("null output");
    const auto v = model(mode,density).after_deposited_energy(u,mass,joules);
    *out = v.state; *residual = v.signed_numerical_energy_residual_j; });
}
}
