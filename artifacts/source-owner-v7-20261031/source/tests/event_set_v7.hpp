#pragma once
// Deterministic wall-event sets of the v7 qualification (shared by test_surface_v7 and test_wall_gpu_v7).
// Integer generator (splitmix64) with exact integer -> binary64 conversions only, so the set is the same bits
// on every build. Mix: electrons 0.01 eV - 1 MeV, ions 0.1 eV - 100 keV, neutrals 1 meV - 10 eV (log-uniform),
// cos(incidence) uniform in (0, 1], normals +-x / +-y, four patch states, three weights, unique 64-bit ids.
#include "surface_data_v7.hpp"
#include <cmath>
#include <cstdint>
#include <vector>

namespace v7test {
struct SplitMix {
  std::uint64_t s;
  std::uint64_t next() { std::uint64_t z = (s += 0x9E3779B97F4A7C15ull); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
  double u() { return double((next() >> 11) | 1ull) * 0x1p-53; }   // open (0, 1)
};
// Patch states (work function, temperature, coverage, adsorbed reference) of the four test patches.
inline std::vector<SW7Params> patch_params(const SW7Params& base) {
  std::vector<SW7Params> v(4, base);
  const double phi[4] = {4.35, 4.3, 3.6, 5.1}, T[4] = {300, 600, 1200, 300}, th[4] = {0.33, 0.1, 0.0, 0.9}, ads[4] = {0.87, 0.8, 0.0, 0.6};
  for (int i = 0; i < 4; ++i) { v[i].phi_eV = phi[i]; v[i].T_K = T[i]; v[i].theta = th[i]; v[i].adsorption_eV = ads[i]; }
  v[2].enable_sticking = 0;   // a patch without water stock
  return v;
}
inline std::vector<SW7Hit> event_set(std::size_t n, std::uint64_t seed, const std::array<SW7Species, SW7_NSPECIES>& table) {
  SplitMix g{seed};
  std::vector<SW7Hit> out(n);
  const double c = 299792458.0;
  for (std::size_t i = 0; i < n; ++i) {
    SW7Hit h{};
    h.id = (std::uint64_t(i) * 0x9E3779B97F4A7C15ull) ^ 0x0123456789ABCDEFull;
    const double r = g.u();
    unsigned sp;
    if (r < .5) sp = SW7_E;
    else if (r < .7) sp = SW7_H2O;
    else if (r < .85) sp = SW7_H2OP;
    else if (r < .9) sp = SW7_H3OP;
    else sp = 4 + unsigned(g.next() % 8);   // OH+, O+, H+, H2+, OH, O, H, H2
    h.species = sp;
    h.patch = unsigned(g.next() % 4);
    const double wsel = g.u();
    h.weight = wsel < .6 ? 1e6 : (wsel < .9 ? 2.5e5 : 1.0);
    const auto& s = table[sp];
    double lo, hi;
    if (s.kind == 0) { lo = 1e-2; hi = 1e6; } else if (s.kind == 1) { lo = 0.1; hi = 1e5; } else { lo = 1e-3; hi = 10; }
    const double E = std::exp(std::log(lo) + (std::log(hi) - std::log(lo)) * g.u());
    const double K = E * 1.602176634e-19, m = s.mass_kg;
    const double us = std::sqrt(K * (K + 2 * m * c * c)) / (m * c);
    const double ct = g.u(), st = std::sqrt(1 - ct * ct), az = 6.283185307179586 * g.u();
    const int face = int(g.next() % 4);
    double n3[3] = {0, 0, 0};
    n3[face / 2] = (face & 1) ? 1 : -1;
    // tangent directions orthogonal to the normal axis
    const int a1 = (face / 2 == 0) ? 1 : 0, a2 = 2;
    double d[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k) d[k] = -ct * n3[k];
    d[a1] += st * std::cos(az);
    d[a2] += st * std::sin(az);
    for (int k = 0; k < 3; ++k) { h.u[k] = us * d[k]; h.n[k] = n3[k]; }
    if (!(-(h.u[0] * h.n[0] + h.u[1] * h.n[1] + h.u[2] * h.n[2]) > 0)) { h.u[face / 2] = -n3[face / 2] * us; }   // grazing guard
    out[i] = h;
  }
  return out;
}
}  // namespace v7test
