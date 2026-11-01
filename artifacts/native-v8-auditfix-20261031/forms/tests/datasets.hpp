// datasets.hpp -- deterministic test data for the v7 forms (CPU and GPU tests).
#pragma once
#include "pf_v7/forms.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

namespace v7test {
using namespace pf_v7;
inline int failures = 0, checks = 0;
inline void check(bool ok, const std::string &what) {
  ++checks;
  if (!ok) {
    ++failures;
    std::printf("FAIL %s\n", what.c_str());
  }
}
inline std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }
// Distance in units in the last place between two finite doubles (ordered
// integer distance).
inline std::uint64_t ulps(double a, double b) {
  auto ord = [](double x) -> std::int64_t {
    const auto u = static_cast<std::int64_t>(std::bit_cast<std::uint64_t>(x));
    return u < 0 ? std::int64_t(0x8000000000000000ull) - u : u;
  };
  const auto d = ord(a) - ord(b);
  return d < 0 ? std::uint64_t(-d) : std::uint64_t(d);
}
// v6 form: sequential x87 long-double sum of (long double)a * b.
inline double x87_dot(std::span<const double> a, std::span<const double> b) {
  long double s = 0;
  for (std::size_t i = 0; i < a.size(); ++i)
    s += static_cast<long double>(a[i]) * b[i];
  return double(s);
}
inline double x87_kinetic(std::span<const Carrier> census, std::span<const double> u_end) {
  long double s = 0;
  for (std::size_t i = 0; i < census.size(); ++i)
    s += static_cast<long double>(census[i].weight) *
         (kinetic_energy_J(census[i].mass_kg, &u_end[3 * i]) -
          kinetic_energy_J(census[i].mass_kg, census[i].proper_velocity_m_s));
  return double(s);
}

struct Motion {
  std::vector<Carrier> census;    // ascending ID
  std::vector<double> u_end;      // 3 per census carrier
  std::vector<Carrier> endpoints; // survivors
  std::vector<Sink> sinks;        // removed
  std::vector<Path> paths;
  std::vector<Birth> births;
  std::vector<CensusRecord> census_records;
  std::vector<std::string> states;
};
// Electrons and neutrals with physical weights and velocities; ~1 % sinks,
// multi-carrier packets (2 co-located carriers) every 7th packet, sparse IDs.
inline Motion motion(std::size_t n, std::uint64_t seed) {
  std::mt19937_64 g(seed);
  std::uniform_real_distribution<double> U(0, 1);
  std::normal_distribution<double> N(0, 1);
  Motion m;
  m.states = {"H2O", "e", "n"};
  std::uint64_t id = 1 + (g() % 1000);
  std::uint64_t key = id;
  for (std::size_t i = 0; i < n; ++i) {
    const bool second = i > 0 && (i % 7 == 1); // second carrier of a 2-carrier packet
    Carrier c{};
    if (second) {
      c = m.census.back();
      c.id = id;
      c.weight = 1e9 * (1 + 999 * U(g));
    } else {
      key = id;
      c.id = id;
      const bool electron = g() % 3 != 0;
      for (int d = 0; d < 3; ++d) {
        c.position_m[d] = 0.001 + 0.004 * U(g);
        c.proper_velocity_m_s[d] = (electron ? 3e6 : 600.0) * N(g);
      }
      c.charge_C = electron ? -1.602176634e-19 : 0;
      c.mass_kg = electron ? 9.1093837015e-31 : 2.99e-26;
      c.weight = 1e9 * (1 + 999 * U(g));
    }
    m.census.push_back(c);
    const std::uint64_t state = c.charge_C ? 1 : (c.mass_kg > 1e-26 ? 0 : 2);
    m.census_records.push_back({c, key, state});
    id += 1 + (g() % 3 == 0 ? g() % 50 : 0);
  }
  // End states: survivors move; every 97th packet sinks (all its carriers).
  bool sink_packet = false;
  for (std::size_t i = 0; i < n; ++i) {
    const auto &c = m.census[i];
    if (m.census_records[i].packet_key == c.id)
      sink_packet = (i % 97) == 5;
    double u1[3];
    if (i > 0 && m.census_records[i].packet_key != c.id) { // co-located carrier: same end state
      for (int d = 0; d < 3; ++d)
        u1[d] = m.u_end[3 * (i - 1) + d];
    } else
      for (int d = 0; d < 3; ++d)
        u1[d] = c.proper_velocity_m_s[d] * (1 + 1e-4 * N(g)) + (c.charge_C ? 1e3 * N(g) : 0);
    for (int d = 0; d < 3; ++d)
      m.u_end.push_back(u1[d]);
    Path p{c.id, {}, {}, 1e-14};
    for (int d = 0; d < 3; ++d) {
      p.field_start_m[d] = c.position_m[d] - 0.0025;
      p.field_displacement_m[d] = 1e-14 * u1[d];
    }
    if (sink_packet) {
      Sink k{};
      k.id = c.id;
      k.kind = 1;
      k.surface = 0;
      k.axis = 1;
      k.plane_m = 0.0023;
      k.fraction = U(g);
      for (int d = 0; d < 3; ++d) {
        k.segment_m[d] = 1e-14 * u1[d];
        k.hit_point_m[d] = c.position_m[d] + k.fraction * k.segment_m[d];
        k.incident_proper_velocity_m_s[d] = u1[d];
      }
      k.flight_s = 1e-14 * k.fraction;
      k.discarded_flight_s = 1e-14 - k.flight_s;
      m.sinks.push_back(k);
      p.flight_s = k.flight_s;
    } else {
      Carrier e = c;
      for (int d = 0; d < 3; ++d) {
        e.position_m[d] = c.position_m[d] + p.field_displacement_m[d];
        e.proper_velocity_m_s[d] = u1[d];
      }
      m.endpoints.push_back(e);
    }
    m.paths.push_back(p);
    if (i % 101 == 3)
      m.births.push_back({c.id, 5e-15 * U(g), 1e-14 - 5e-15 * U(g)});
  }
  return m;
}
// Adversarial exact-sum sets (a_i, b_i): massive cancellation, wide exponent
// range, overflowing products, subnormal products, ties.
struct ProductSet {
  std::string name;
  std::vector<double> a, b;
};
inline std::vector<ProductSet> adversarial(std::uint64_t seed) {
  std::mt19937_64 g(seed);
  std::uniform_real_distribution<double> U(-1, 1);
  std::vector<ProductSet> out;
  auto add_pair = [](ProductSet &s, double a, double b) {
    s.a.push_back(a);
    s.b.push_back(b);
  };
  { // x, -x pairs over 600 decades plus one tiny residue
    ProductSet s{"cancel-pairs-wide", {}, {}};
    for (int i = 0; i < 2000; ++i) {
      const double a = U(g) * std::ldexp(1.0, int(g() % 1200) - 600), b = U(g) * std::ldexp(1.0, int(g() % 400) - 200);
      add_pair(s, a, b);
      add_pair(s, -a, b);
    }
    add_pair(s, 3e-300, 7e-20);
    out.push_back(s);
  }
  { // overflowing products that cancel to a moderate result
    ProductSet s{"overflow-cancel", {}, {}};
    add_pair(s, 1e308, 1e10);
    add_pair(s, 1e308, 1e10);
    add_pair(s, -1e308, 2e10);
    add_pair(s, 1.5, 2.5);
    out.push_back(s);
  }
  { // subnormal products: 3 x 2^-1075 = 1.5 x 2^-1074 (tie -> even 2 x 2^-1074)
    ProductSet s{"subnormal-tie", {}, {}};
    for (int i = 0; i < 3; ++i)
      add_pair(s, std::ldexp(1.0, -1074), 0.5);
    out.push_back(s);
  }
  { // 1 + 2^-53 (tie -> 1) and 1 + 2^-53 + 2^-200 (-> 1 + 2^-52)
    ProductSet s{"tie-even", {}, {}};
    add_pair(s, 1, 1);
    add_pair(s, std::ldexp(1.0, -53), 1);
    out.push_back(s);
    ProductSet t{"tie-sticky", {}, {}};
    add_pair(t, 1, 1);
    add_pair(t, std::ldexp(1.0, -53), 1);
    add_pair(t, std::ldexp(1.0, -100), std::ldexp(1.0, -100));
    out.push_back(t);
  }
  { // negative zero only -> +0
    ProductSet s{"negative-zero", {}, {}};
    add_pair(s, -0.0, 1);
    add_pair(s, 0.0, -5);
    out.push_back(s);
  }
  { // large sum with a tiny exact residue hidden behind cancellation chains
    ProductSet s{"chain-1e5", {}, {}};
    for (int i = 0; i < 50000; ++i) {
      const double a = U(g) * 1e20, b = U(g);
      add_pair(s, a, b);
      add_pair(s, -b, a);
    }
    add_pair(s, 1e-30, 1e-30);
    out.push_back(s);
  }
  { // result near the largest finite value
    ProductSet s{"near-max", {}, {}};
    add_pair(s, 1.7976931348623157e308, 1);
    add_pair(s, std::ldexp(1.0, 970), 0.49999999999);
    out.push_back(s);
  }
  { // just overflowing after rounding
    ProductSet s{"round-overflow", {}, {}};
    add_pair(s, 1.7976931348623157e308, 1);
    add_pair(s, std::ldexp(1.0, 970), 1);
    out.push_back(s);
  }
  return out;
}
// Random product sets with heavy cancellation of n terms.
inline ProductSet random_cancel(std::size_t n, std::uint64_t seed) {
  std::mt19937_64 g(seed);
  std::uniform_real_distribution<double> U(-1, 1);
  ProductSet s{"random-cancel-" + std::to_string(n), {}, {}};
  for (std::size_t i = 0; i < n; ++i) {
    const double scale = std::ldexp(1.0, int(g() % 80) - 40);
    s.a.push_back(U(g) * scale * 1e9);
    s.b.push_back(U(g) * 1e-18);
    if (i % 2 == 1) { // nearly cancel the previous term
      s.a.back() = -s.a[i - 1] * (1 + 1e-12 * U(g));
      s.b.back() = s.b[i - 1];
    }
  }
  return s;
}
} // namespace v7test
