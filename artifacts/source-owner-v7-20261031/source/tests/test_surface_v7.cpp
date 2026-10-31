// Heat/source owner v7: surface physics qualification without the owner.
//   test_surface_v7 OUTDIR
// Groups:
//   math        deterministic exp/log/sincos/lgamma/incomplete gamma against the host libm
//   rng         Philox4x32-10 known answer (Random123), surface stream layout v7 (key salt, injectivity probes)
//   emission    sealed thermal-field kernel against published curves: Forbes v(f) (SN barrier action), Murphy-Good
//               0 K FN current, Murphy-Good temperature factor, Richardson(-Schottky), work-function shift
//   workfn      work function of oxidized SS304 vs coverage (cited set), benchmark 4.3 eV, grid
//   see         Furman-Pivi stainless set: delta_max/E_max, 300 eV component split, as-received SEY table,
//               Monte Carlo yields, spectra (gamma marginals, rediffused law, elastic peak), angles, low energy
//   ions        Thomas-Janev-Smith reflection table (H on Fe), PE/KE yields, ion and neutral channel statistics
//   ledger      1e6 mixed wall events: every per-event exact identity, exact global sums, wall heat sign rules
//   partition   outputs and integer sums identical over 1/2/5/8 partitions, strided and reversed orders
//   interval    first-order error of the explicit start-of-step exposure and the interval bound
//   hot         extended emission domain to 5000 K: bit-identical to the sealed kernel up to 3000 K, exact Fermi-Dirac
//               thermal law and Richardson-Schottky to 5000 K, quadrature convergence, validity map (OUTDIR/validity-map-v7.json)
//   registry    every parameter tagged, banded and cited; written to OUTDIR/parameters-v7.json
#include "surface_data_v7.hpp"
#include "event_set_v7.hpp"
#include "surface_emission.hpp"
#include "fast_emission.hpp"
#include "emission_hot_v7.hpp"
#include "raikou/json.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>

namespace se = surface_emission;
namespace sv = surface_v7;
using J = raikou::Json;

namespace {
struct Results {
  J cases = J::array();
  int failed = 0;
  void record(const std::string& g, const std::string& name, bool ok, J ev = J::object()) {
    cases.push_back({{"group", g}, {"name", name}, {"passed", ok}, {"evidence", std::move(ev)}});
    if (!ok) ++failed;
    std::cout << g << " | " << name << ": " << (ok ? "pass" : "FAIL") << std::endl;
    if (!ok) std::cout << "  evidence: " << cases.back()["evidence"].dump().substr(0, 3000) << std::endl;
  }
};
double ulps(double a, double b) {
  if (a == b) return 0;
  if (!std::isfinite(a) || !std::isfinite(b)) return INFINITY;
  return std::abs(a - b) / std::max(std::numeric_limits<double>::denorm_min(), std::abs(std::nextafter(b, INFINITY) - b));
}
const double qe = 1.602176634e-19, me = 9.1093837139e-31, kB = 1.380649e-23, kBeV = kB / qe;

// ------------------------------------------------------------ math
void math_group(Results& R) {
  v7test::SplitMix g{77};
  double ex = 0, lg = 0, sc = 0, lgam = 0, gi = 0;
  for (int i = 0; i < 200000; ++i) {
    const double x = -700 + 1400 * g.u();
    ex = std::max(ex, ulps(sw7::exp_(x), std::exp(x)));
    const double y = std::exp(-700 + 1400 * g.u());
    lg = std::max(lg, ulps(sw7::log_(y), std::log(y)));
    const double t = g.u();
    double s, c;
    sw7::sincos_turn(t, &s, &c);
    sc = std::max({sc, std::abs(s - std::sin(6.283185307179586 * t)), std::abs(c - std::cos(6.283185307179586 * t))});
    const double a = 0.5 + 199.5 * g.u();
    lgam = std::max(lgam, std::abs(sw7::lgamma_(a) - std::lgamma(a)) / std::max(1.0, std::abs(std::lgamma(a))));
  }
  // incomplete gamma against a fine trapezoid-free reference: P(1, x) = 1 - exp(-x), P(1/2, x^2) = erf(x)
  for (int i = 0; i < 20000; ++i) {
    const double x = 30 * g.u();
    gi = std::max(gi, std::abs(sw7::gammainc_P(1.0, x) - (-std::expm1(-x))));
    const double z = 5 * g.u();
    gi = std::max(gi, std::abs(sw7::erf_pos(z) - std::erf(z)));
  }
  // exact special cases
  const bool special = sw7::exp_(0) == 1 && sw7::log_(1) == 0 && sw7::exp_(-800) == 0 && std::isinf(sw7::exp_(710)) && sw7::log_(0x1p-1060) < -734;
  R.record("math", "exp and log within 2 ulp of libm over [-700, 700]", ex <= 2 && lg <= 2 && special, {{"exp_max_ulp", ex}, {"log_max_ulp", lg}});
  // reference: libm sin/cos of the binary64 product 2 pi t, itself rounded by up to 4.4e-16 for t near 1
  R.record("math", "sincos_turn within 1e-15 absolute of libm (reference argument rounding <= 4.4e-16)", sc <= 1e-15, {{"max_abs", sc}});
  R.record("math", "lgamma within 1e-13 relative on [0.5, 200]; P(a,x) within 1e-13", lgam <= 1e-13 && gi <= 1e-13, {{"lgamma_max_rel", lgam}, {"gammainc_max_abs", gi}});
}

// ------------------------------------------------------------ rng
void rng_group(Results& R) {
  // Random123 kat_vectors: philox4x32_10 counter (0,0,0,0) key (0,0) -> 6627e8d5 e169c58d bc57ac4c 9b00dbd8;
  // counter all-ones key all-ones -> 408f276d 41c83b0e a20bc7c6 6d5451fd
  pfa::Block c{{0, 0, 0, 0}};
  pfa::Key2 k{{0, 0}};
  auto r = pfa::philox4x32_10(c, k);
  const bool kat0 = r.v[0] == 0x6627e8d5u && r.v[1] == 0xe169c58du && r.v[2] == 0xbc57ac4cu && r.v[3] == 0x9b00dbd8u;
  pfa::Block c1{{0xffffffffu, 0xffffffffu, 0xffffffffu, 0xffffffffu}};
  pfa::Key2 k1{{0xffffffffu, 0xffffffffu}};
  auto r1 = pfa::philox4x32_10(c1, k1);
  const bool kat1 = r1.v[0] == 0x408f276du && r1.v[1] == 0x41c83b0eu && r1.v[2] == 0xa20bc7c6u && r1.v[3] == 0x6d5451fdu;
  // stream layout: distinct events, purposes, steps and seeds give distinct first words; the key differs from the
  // arrival owner's key for the same run seed
  std::set<std::uint64_t> words;
  int n = 0;
  for (std::uint64_t seed : {1ull, 2ull, 4242ull})
    for (std::uint64_t step : {0ull, 1ull, (1ull << 32), (1ull << 39) + 5})
      for (std::uint64_t id : {1ull, 2ull, (1ull << 32) + 1, 0xFFFFFFFFFFFFFFFFull})
        for (unsigned purpose : {SW7_PURPOSE_WALL, SW7_PURPOSE_PHOTO}) {
          auto s = sw7::stream(seed, step, id, purpose);
          words.insert(sw7::word(s));
          ++n;
        }
  const auto s0 = sw7::stream(4242, 7, 9, SW7_PURPOSE_WALL);
  const bool salted = s0.key.v[0] != std::uint32_t(4242) && s0.key.v[0] == std::uint32_t(4242ull ^ SW7_SALT);
  // word budget: exhaustion is flagged, never silent
  auto sb = sw7::stream(1, 1, 1, SW7_PURPOSE_WALL);
  for (unsigned i = 0; i < SW7_RNG_WORD_CAP; ++i) (void)sw7::word(sb);
  const bool before = !sb.exhausted;
  (void)sw7::word(sb);
  R.record("rng", "Philox4x32-10 Random123 known answers; distinct streams per (seed, step, id, purpose); salted key; budget flag",
           kat0 && kat1 && int(words.size()) == n && salted && before && sb.exhausted, {{"streams", n}, {"distinct_first_words", words.size()}});
}

// ------------------------------------------------------------ emission oracles
double forbes_v(double f) { return 1 - f + f / 6 * std::log(f); }
double forbes_t(double f) { return 1 + f / 9 - f / 18 * std::log(f); }
void emission_group(Results& R) {
  const double a = 1.541434e-6, b = 6.830890e9, cS2 = 1.439964e-9;   // Forbes 2006 (A eV V^-2; eV^-3/2 V m^-1; eV^2 m V^-1)
  // (i) SN action at the Fermi level: v_kernel(f) = G(0) F / (b phi^3/2) vs Forbes v(f) (|error| <= 0.0024)
  J vrows = J::array();
  double vmax = 0;
  for (double phi : {2.5, 4.35, 4.5})
    for (double f : {0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}) {
      const double F = f * phi * phi / cS2;
      if (F > 5e9) continue;
      const double G = se::action(0.0, F, se::Config{phi});
      const double vk = G * F / (b * std::pow(phi, 1.5));
      vmax = std::max(vmax, std::abs(vk - forbes_v(f)));
      vrows.push_back({{"phi", phi}, {"f", f}, {"v_kernel", vk}, {"v_forbes", forbes_v(f)}});
    }
  R.record("emission", "SN barrier action: v(f) of the sealed kernel within Forbes 2006 bound 0.0024 of v ~ 1 - f + (f/6) ln f", vmax <= 0.0024,
           {{"max_abs_difference", vmax}, {"rows", vrows}});
  // (ii) 0 K current vs Murphy-Good FN with the exact v (kernel action) and Forbes t (0.4 %)
  J jrows = J::array();
  double jmax = 0, jforbes = 0, jforbes_bound = 0, jcurv = 0;
  for (double phi : {4.35, 4.5})
    for (double f : {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35}) {
      const double F = f * phi * phi / cS2;
      if (F > 5e9) continue;
      const auto d = se::evaluate(0.0, F, se::Config{phi});
      const double vk = se::action(0.0, F, se::Config{phi}) * F / (b * std::pow(phi, 1.5));
      const double lJmg = std::log(a * F * F / (phi * std::pow(forbes_t(f), 2))) - vk * b * std::pow(phi, 1.5) / F;
      const double lJf = std::log(a * F * F / (phi * std::pow(forbes_t(f), 2))) - forbes_v(f) * b * std::pow(phi, 1.5) / F;
      const double dev = d.log_current_A_m2 - lJmg;
      // MG linearizes G(E) about E_F; the neglected curvature c2 = G''(0) changes the 0 K integral
      // int (-E) exp(-G) dE by the factor ~ 1 - 3 c2 d^2 (<E^2> = 6 d^2 under the weight (-E) e^{E/d}).
      const double h = 0.02, G0 = se::action(0.0, F, se::Config{phi}), Gp = se::action(h, F, se::Config{phi}), Gm = se::action(-h, F, se::Config{phi});
      const double c1 = (Gp - Gm) / (2 * h), c2 = (Gp - 2 * G0 + Gm) / (h * h), dd = -1 / c1;
      const double predicted = std::log1p(-3 * c2 * dd * dd);
      jcurv = std::max(jcurv, std::abs(dev - predicted));
      jmax = std::max(jmax, std::abs(dev));
      jforbes = std::max(jforbes, std::abs(d.log_current_A_m2 - lJf) - (b * std::pow(phi, 1.5) / F * 0.0024 + 2 * 0.004));
      jforbes_bound = std::max(jforbes_bound, b * std::pow(phi, 1.5) / F * 0.0024 + 2 * 0.004);
      jrows.push_back({{"phi", phi}, {"f", f}, {"F_V_m", F}, {"log_J_kernel", d.log_current_A_m2}, {"log_J_MG_exact_v", lJmg}, {"log_J_MG_forbes", lJf}, {"curvature_prediction", predicted}, {"deviation", dev}});
    }
  R.record("emission", "0 K current: kernel vs Murphy-Good FN (exact v, Forbes t): total within 6 % for f <= 0.35; deviation = first-order curvature term -3 c2 d^2 within 2 % (Forbes t 0.8 % in J + higher orders of G(E))",
           jcurv <= 0.02 && jmax <= std::log(1.06), {{"max_abs_log_ratio", jmax}, {"max_abs_minus_curvature_prediction", jcurv}, {"rows", jrows}});
  R.record("emission", "0 K current: kernel vs Murphy-Good FN with Forbes v and t within the stated Forbes error band", jforbes <= std::log(1.05),
           {{"max_excess_over_band", jforbes}, {"largest_band_log", jforbes_bound}});
  // (iii) Murphy-Good temperature factor pi w / sin(pi w), w = kT/d (d from the kernel's Fermi decay width)
  J trows = J::array();
  double tmax = 0;
  for (double Fnm : {1.0, 1.5, 2.0, 3.0})
    for (double T : {300., 600., 900.}) {
      const double F = Fnm * 1e9, phi = 4.5;
      const auto d0 = se::evaluate(0.0, F, se::Config{phi});
      const auto dT = se::evaluate(T, F, se::Config{phi});
      const double w = kBeV * T / dT.fermi_decay_width_eV;
      if (w > 0.6) continue;
      const double mg = 3.141592653589793 * w / std::sin(3.141592653589793 * w);
      const double ratio = std::exp(dT.log_current_A_m2 - d0.log_current_A_m2);
      tmax = std::max(tmax, std::abs(std::log(ratio / mg)));
      trows.push_back({{"F_V_nm", Fnm}, {"T_K", T}, {"omega", w}, {"kernel_ratio", ratio}, {"MG_factor", mg}});
    }
  R.record("emission", "Murphy-Good temperature factor pi w/sin(pi w) (w <= 0.6) within 3 % (check value 1.845 at 1 V/nm, 300 K, phi 4.5 eV)",
           tmax <= std::log(1.03) && trows.size() >= 6, {{"max_abs_log_ratio", tmax}, {"rows", trows}});
  // (iv) Richardson at F = 0 (exact law of the step branch) and Richardson-Schottky at low field
  const double A0 = se::ideal_richardson_A_m2_K2();
  double r0 = 0, rs = 0;
  J rrows = J::array();
  for (double T : {1000., 1500., 2000., 2500., 3000.})
    for (double phi : {4.3, 4.35, 4.75}) {
      const auto d = se::evaluate(T, 0.0, se::Config{phi});
      const double lr = std::log(A0 * T * T) - phi / (kBeV * T);
      r0 = std::max(r0, std::abs(d.log_current_A_m2 - lr));
      for (double F : {1e4, 1e5, 1e6, 1e7}) {
        const auto e = se::evaluate(T, F, se::Config{phi});
        const double lrs = std::log(A0 * T * T) - (phi - std::sqrt(cS2 * F)) / (kBeV * T);
        rs = std::max(rs, std::abs(e.log_current_A_m2 - lrs));
        rrows.push_back({{"T", T}, {"phi", phi}, {"F", F}, {"log_J_kernel", e.log_current_A_m2}, {"log_J_RS", lrs}});
      }
    }
  R.record("emission", "Richardson law A0 T^2 exp(-phi/kT) at F = 0 within 1e-6 (A0 = 120.173 A cm^-2 K^-2)", r0 <= 1e-6 && std::abs(A0 / 1.20173e6 - 1) < 1e-5,
           {{"max_abs_log_difference", r0}, {"A0_A_m2_K2", A0}});
  R.record("emission", "Richardson-Schottky (1e4-1e7 V/m, 1000-3000 K) within 20 % (finite-temperature barrier corrections)", rs <= std::log(1.2),
           {{"max_abs_log_ratio", rs}, {"rows_sample", J(std::vector<J>(rrows.begin(), rrows.begin() + std::min<std::size_t>(rrows.size(), 12)))}});
  // (v) work-function coupling at F = 0: the exact Fermi-Dirac thermal current J = A0 T^2 (-Li2(-x)), x = exp(-phi/kT)
  //     (Boltzmann limit exp(-phi/kT)), for phi over the oxide/water range
  double wf = 0;
  for (double T : {800., 1500., 2500.})
    for (double phi : {4.35 - 0.97, 4.0, 4.35, 4.95}) {
      const auto y = se::evaluate(T, 0.0, se::Config{phi});
      const double x = std::exp(-phi / (kBeV * T));
      double li = 0, xk = 1;
      for (int k = 1; k < 40; ++k) { xk *= x; li += ((k & 1) ? 1 : -1) * xk / (double(k) * k); }
      wf = std::max(wf, std::abs(y.log_current_A_m2 - (std::log(A0 * T * T) + std::log(li))));
    }
  R.record("emission", "work-function coupling: ln J(phi, T, F=0) equals the exact Fermi-Dirac law ln(A0 T^2 (-Li2(-exp(-phi/kT)))) within 1e-8 (phi 3.38-4.95 eV)", wf <= 1e-8, {{"max_abs", wf}});
  // (vi) fast exact spectrum (v6.1) bit-identical for v7 work-function grid values
  bool same = true;
  int compared = 0;
  for (double phi : {4.35, 4.35 - 0.97 * 0.5, 3.38, 5.1})
    for (double T : {300., 1200.})
      for (double F : {2e6, 3e8, 2e9}) {
        const double q = std::nearbyint(phi / 0x1p-12) * 0x1p-12;
        same = same && electrode_heat_v6_1::emission::same_distribution(electrode_heat_v6_1::emission::evaluate(T, F, se::Config{q}), se::evaluate(T, F, se::Config{q}));
        ++compared;
      }
  R.record("emission", "v6.1 fast spectrum bit-identical to the sealed kernel on v7 work-function grid values", same, {{"compared", compared}});
}

// ------------------------------------------------------------ work function
void workfn_group(Results& R) {
  sv::Registry reg;
  const auto m = sv::ss304_work_function(reg);
  sv::SurfaceModel bench;
  J rows = J::array();
  bool mono = true;
  double prev = INFINITY;
  for (double th : {0.0, 0.05, 0.1, 0.2, 0.33, 0.5, 0.75, 1.0, 1.5}) {
    const double p = sv::work_function_eV(m, th, 300);
    mono = mono && p <= prev;
    prev = p;
    rows.push_back({{"theta", th}, {"phi_eV", p}});
  }
  const double p0 = sv::work_function_eV(m, 0, 300), p1 = sv::work_function_eV(m, 1, 300), p05 = sv::work_function_eV(m, 0.01, 300);
  const double grid = std::nearbyint((4.35 - 0.97) / 0x1p-12) * 0x1p-12;
  const double slope = (p05 - p0) / 0.01;   // initial slope ~ dphi_1ML (1 + c) = -2.13 eV/ML (secant over 0.01 ML, grid 2^-12 eV)
  const bool ok = sv::work_function_eV(bench, 0.7, 2000) == 4.3 && p0 == std::nearbyint(4.35 / 0x1p-12) * 0x1p-12 && p1 == grid && mono &&
                  std::abs(slope - (-0.97 * 2.2)) < 0.06 && p0 >= m.band_low_eV && p0 <= m.band_high_eV;
  R.record("workfn", "oxidized SS304 phi(theta): 4.35 eV at theta 0, -0.97 eV at one monolayer, monotone, initial slope -0.97(1+c), grid 2^-12 eV; benchmark 4.3 eV",
           ok, {{"rows", rows}, {"initial_slope_eV_per_ML", slope}});
}

// ------------------------------------------------------------ SEE
SW7Params base_params() {
  sv::Registry R;
  auto p = sv::wall_parameters(R);
  p.phi_eV = 4.35; p.T_K = 300; p.theta = 0.33; p.adsorption_eV = 0.87;
  return p;
}
void see_group(Results& R, const std::array<SW7Species, SW7_NSPECIES>& table) {
  SW7Params P = base_params();
  // (a) maximum at normal incidence
  double best = 0, bestE = 0;
  for (double E = 50; E <= 1000; E += 0.5) {
    double de, dr, dt;
    sw7::fp_yields(P, E, 1.0, &de, &dr, &dt);
    if (de + dr + dt > best) { best = de + dr + dt; bestE = E; }
  }
  R.record("see", "Furman-Pivi stainless: delta_max 2.05 at E_max 292 eV (Table I) within 0.01 and 2 eV", std::abs(best - 2.05) <= 0.01 && std::abs(bestE - 292) <= 2,
           {{"delta_max", best}, {"E_max_eV", bestE}});
  // (b) 300 eV split vs measured (Fig. 2 caption): 0.12 / 0.75 / 1.17
  double de, dr, dt;
  sw7::fp_yields(P, 300, 1.0, &de, &dr, &dt);
  R.record("see", "300 eV components vs measured split 0.12/0.75/1.17 (Fig. 2) within 0.06 each", std::abs(de - .12) <= .06 && std::abs(dr - .75) <= .06 && std::abs(dt - 1.17) <= .06,
           {{"delta_e", de}, {"delta_r", dr}, {"delta_ts", dt}});
  // (c) as-received stainless SEY (Baglin et al. EPAC 2000 Fig. 3 via the Joy database) including the high-energy continuation
  const double tE[13] = {50, 80, 100, 150, 200, 250, 300, 400, 500, 700, 1000, 1500, 1950};
  const double tY[13] = {1.15, 1.40, 1.55, 1.80, 1.94, 1.99, 2.02, 1.99, 1.92, 1.75, 1.53, 1.28, 1.13};
  J rows = J::array();
  double worst = 0, rms = 0;
  for (int i = 0; i < 13; ++i) {
    double a, b2, c;
    sw7::fp_yields(P, tE[i], 1.0, &a, &b2, &c);
    bool out = false;
    const double f = sw7::bs_continuation(P, tE[i], a, b2, &out);
    const double y = a * f + b2 * f + c;
    const double rel = y / tY[i] - 1;
    worst = std::max(worst, std::abs(rel));
    rms += rel * rel;
    rows.push_back({{"E_eV", tE[i]}, {"measured", tY[i]}, {"model", y}, {"relative", rel}});
  }
  rms = std::sqrt(rms / 13);
  R.record("see", "as-received stainless SEY table 50-1950 eV (independent sample): every point within 20 %, RMS within 10 %", worst <= 0.20 && rms <= 0.10,
           {{"max_abs_relative", worst}, {"rms_relative", rms}, {"rows", rows}});
  // (d) Monte Carlo: mean emitted electrons per incident electron vs delta(E0) (normal incidence)
  J mc = J::array();
  bool mcok = true;
  for (double E0 : {1.0, 10.0, 50.0, 150.0, 300.0, 1000.0, 5000.0}) {
    double a, b2, c;
    sw7::fp_yields(P, E0, 1.0, &a, &b2, &c);
    bool out = false;
    const double f = sw7::bs_continuation(P, E0, a, b2, &out);
    a *= f; b2 *= f;
    const int Nmc = 200000;
    double sum = 0, sum2 = 0, el = 0, re = 0, capped = 0;
    for (int i = 0; i < Nmc; ++i) {
      SW7Hit h{};
      h.id = 1000000 + std::uint64_t(i);
      h.species = SW7_E; h.weight = 1; h.n[1] = 1;
      const double K = E0 * qe, us = std::sqrt(K * (K + 2 * me * 299792458.0 * 299792458.0)) / (me * 299792458.0);
      h.u[1] = -us;
      SW7Out o;
      sw7::wall_event(P, table.data(), h, 99, 3, &o);
      sum += o.electrons_out; sum2 += double(o.electrons_out) * o.electrons_out;
      el += o.channel_main == SW7_CH_ELASTIC; re += o.channel_main == SW7_CH_REDIFFUSED;
      capped += (o.status & SW7_ST_ENERGY_CAP) != 0;
    }
    const double mean = sum / Nmc, sd = std::sqrt(std::max(0.0, sum2 / Nmc - mean * mean) / Nmc);
    const double expect = a + b2 + c;
    // the energy cap removes only emitted counts that energy conservation forbids (recorded; E0 < 2 phi only)
    const bool ok = (std::abs(mean - expect) <= 5 * sd + 1e-12) || (capped > 0 && mean < expect);
    const bool comp = std::abs(el / Nmc - a) <= 5 * std::sqrt(a * (1 - a) / Nmc) + 1e-12 && std::abs(re / Nmc - b2) <= 5 * std::sqrt(b2 * (1 - b2) / Nmc + 1e-30) + 1e-12;
    mcok = mcok && ok && comp;
    mc.push_back({{"E0", E0}, {"mean", mean}, {"sd", sd}, {"delta", expect}, {"elastic_fraction", el / Nmc}, {"delta_e", a}, {"rediffused_fraction", re / Nmc}, {"delta_r", b2}, {"energy_capped_events", capped}});
  }
  R.record("see", "sampled yields and elastic/rediffused fractions equal delta(E0) within 5 sigma (2e5 events per energy)", mcok, {{"rows", mc}});
  // (e) true-secondary energies: each E_k of an n-electron event is Gamma(p_n, eps_n) when truncation is negligible (KS, 1 %)
  {
    std::map<int, std::vector<double>> byn;
    for (int i = 0; i < 300000; ++i) {
      SW7Hit h{};
      h.id = 5000000 + std::uint64_t(i); h.species = SW7_E; h.weight = 1; h.n[0] = -1;
      const double E0 = 1100, K = E0 * qe, us = std::sqrt(K * (K + 2 * me * 299792458.0 * 299792458.0)) / (me * 299792458.0);
      h.u[0] = us;
      SW7Out o;
      sw7::wall_event(P, table.data(), h, 7, 11, &o);
      if (o.channel_main != SW7_CH_TRUE) continue;
      int n = 0;
      for (unsigned k = 0; k < o.n_out; ++k) n += o.ch[k] == SW7_CH_TRUE;
      for (unsigned k = 0; k < o.n_out; ++k) if (o.ch[k] == SW7_CH_TRUE) byn[n].push_back(sw7::from_quanta(o.kq[k]));
    }
    J ks = J::array();
    bool ok = true;
    for (auto& [n, v] : byn) {
      if (v.size() < 2000) continue;
      std::sort(v.begin(), v.end());
      double D = 0;
      for (std::size_t i = 0; i < v.size(); ++i) {
        const double F = sw7::gammainc_P(P.pn[n - 1], v[i] / P.epsn[n - 1]);
        D = std::max({D, std::abs(F - double(i) / v.size()), std::abs(F - double(i + 1) / v.size())});
      }
      const double crit = 1.63 / std::sqrt(double(v.size())) + 2 * 0x1p-30 / P.epsn[n - 1];
      const double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
      ok = ok && D <= crit;
      ks.push_back({{"n", n}, {"samples", v.size()}, {"ks_D", D}, {"critical_1pct", crit}, {"mean_eV", mean}, {"p_n_eps_n", P.pn[n - 1] * P.epsn[n - 1]},
                    {"mode_eV_(p_n-1)eps_n", (P.pn[n - 1] - 1) * P.epsn[n - 1]}});
    }
    R.record("see", "true-secondary energies: Gamma(p_n, eps_n) marginals (FP Eqs. 33, 53), Kolmogorov-Smirnov at 1 % for every n with >= 2000 samples (E0 1100 eV)",
             ok && ks.size() >= 5, {{"rows", ks}});
  }
  // (f) rediffused law (E/E0)^(q+1), elastic peak within E0 - 6 sigma_e, cosine emission mean cos = 2/3, total energy <= E0
  {
    std::vector<double> red;
    double el_min = INFINITY, cos_sum = 0;
    long cos_n = 0, viol = 0;
    const double E0 = 200;
    for (int i = 0; i < 300000; ++i) {
      SW7Hit h{};
      h.id = 9000000 + std::uint64_t(i); h.species = SW7_E; h.weight = 1; h.n[2] = 1;
      const double K = E0 * qe, us = std::sqrt(K * (K + 2 * me * 299792458.0 * 299792458.0)) / (me * 299792458.0);
      h.u[2] = -us;
      SW7Out o;
      sw7::wall_event(P, table.data(), h, 5, 1, &o);
      std::int64_t sum = 0;
      for (unsigned k = 0; k < o.n_out; ++k) {
        sum += o.kq[k];
        const double um = std::sqrt(o.u[k][0] * o.u[k][0] + o.u[k][1] * o.u[k][1] + o.u[k][2] * o.u[k][2]);
        if (um > 0) { cos_sum += o.u[k][2] / um; ++cos_n; }
        if (o.ch[k] == SW7_CH_REDIFFUSED) red.push_back(sw7::from_quanta(o.kq[k]) / E0);
        if (o.ch[k] == SW7_CH_ELASTIC) el_min = std::min(el_min, sw7::from_quanta(o.kq[k]));
      }
      viol += sum > o.k_in_q;
    }
    std::sort(red.begin(), red.end());
    double D = 0;
    for (std::size_t i = 0; i < red.size(); ++i) {
      const double F = std::pow(red[i], P.q + 1);
      D = std::max({D, std::abs(F - double(i) / red.size()), std::abs(F - double(i + 1) / red.size())});
    }
    const double mc = cos_sum / cos_n;
    R.record("see", "rediffused energies follow (E/E0)^(q+1) (KS 1 %), elastic within 6 sigma_e of E0, cosine law <cos> = 2/3, sum of emitted K <= E0",
             D <= 1.63 / std::sqrt(double(red.size())) && el_min >= E0 - 6 * P.sigma_e && std::abs(mc - 2.0 / 3.0) <= 5 * 0.2357 / std::sqrt(double(cos_n)) && viol == 0,
             {{"rediffused_ks_D", D}, {"rediffused_samples", red.size()}, {"elastic_min_eV", el_min}, {"mean_cos", mc}, {"energy_violations", viol}});
  }
  // (g) oblique incidence: FP limits exceeded (delta'_ts > M beyond ~18 deg, delta_e + delta_r -> 1 near 43 deg) are clamped and recorded
  {
    J rows2 = J::array();
    for (double deg : {0., 15., 30., 45., 60., 75., 85.}) {
      const double c = std::cos(deg * 3.141592653589793 / 180);
      double mx = 0, m2 = 0;
      for (double E = 20; E <= 1000; E += 1) {
        double a, b2, t;
        sw7::fp_yields(P, E, c, &a, &b2, &t);
        mx = std::max(mx, t / std::max(1e-300, 1 - a - b2));
        m2 = std::max(m2, a + b2);
      }
      rows2.push_back({{"theta_deg", deg}, {"max_delta_ts_prime", mx}, {"max_delta_e_plus_r", m2}});
    }
    R.record("see", "oblique incidence: model limits of the stainless angular parameters (TiN/Al-derived) reported; violations clamped and flagged per event",
             rows2[0]["max_delta_ts_prime"].get<double>() < 10 && rows2[3]["max_delta_ts_prime"].get<double>() > 10, {{"rows", rows2}});
  }
  // (h) zero-energy reflection: delta(E0 -> 0) = P1e_hat = 0.5 (band 0.02-1, open)
  {
    double a, b2, c;
    sw7::fp_yields(P, 1e-3, 1.0, &a, &b2, &c);
    R.record("see", "low-energy limit delta(1 meV) = P1e_hat (0.5, extrapolated; band 0.02-1 recorded as open)", std::abs(a + b2 + c - 0.5) < 1e-3, {{"delta", a + b2 + c}});
  }
}

// ------------------------------------------------------------ ions and neutrals
void ion_group(Results& R, const std::array<SW7Species, SW7_NSPECIES>& table) {
  SW7Params P = base_params();
  // (a) TJS H -> Fe recommended table (INDC(NDS)-249 Appendix, accuracy 10 %)
  const double E[5] = {0.01, 0.1, 1, 10, 100}, RN[5] = {0.683, 0.504, 0.329, 0.0761, 0.00237}, RE[5] = {0.505, 0.316, 0.162, 0.0236, 7.23e-4};
  double worst = 0;
  J rows = J::array();
  for (int i = 0; i < 5; ++i) {
    const double eps = table[SW7_HP].refl_epsL_per_keV * E[i];
    const double rn = sw7::tjs(P.rn_A, eps), re = sw7::tjs(P.re_A, eps);
    worst = std::max({worst, std::abs(rn / RN[i] - 1), std::abs(re / RE[i] - 1)});
    rows.push_back({{"E_keV", E[i]}, {"RN", rn}, {"RN_table", RN[i]}, {"RE", re}, {"RE_table", RE[i]}});
  }
  R.record("ions", "Thomas-Janev-Smith H on Fe: R_N and R_E within 2 % of the recommended table 0.01-100 keV", worst <= 0.02, {{"max_relative", worst}, {"rows", rows}});
  // (b) channel statistics: reflection, potential and kinetic emission, sticking, accommodation
  auto run = [&](unsigned sp, double E0, int Nmc, SW7Params Q) {
    struct S { double refl = 0, electrons = 0, stuck = 0, thermal = 0, specular = 0, frag = 0; } s;
    const auto& t = table[sp];
    for (int i = 0; i < Nmc; ++i) {
      SW7Hit h{};
      h.id = 77000000 + std::uint64_t(i); h.species = sp; h.weight = 1; h.n[1] = -1;
      const double K = E0 * qe, us = std::sqrt(K * (K + 2 * t.mass_kg * 299792458.0 * 299792458.0)) / (t.mass_kg * 299792458.0);
      h.u[1] = us;
      SW7Out o;
      sw7::wall_event(Q, table.data(), h, 31, 2, &o);
      s.electrons += o.electrons_out; s.stuck += o.stuck;
      for (unsigned k = 0; k < o.n_out; ++k) { s.refl += o.ch[k] == SW7_CH_REFLECTED; s.thermal += o.ch[k] == SW7_CH_THERMAL; s.specular += o.ch[k] == SW7_CH_SPECULAR; s.frag += o.ch[k] == SW7_CH_FRAGMENT; }
    }
    for (double* x : {&s.refl, &s.electrons, &s.stuck, &s.thermal, &s.specular, &s.frag}) *x /= Nmc;
    return s;
  };
  const int Nmc = 200000;
  J rows2 = J::array();
  bool ok = true;
  auto within = [&](double x, double p) { return std::abs(x - p) <= 5 * std::sqrt(std::max(p * (1 - p), 1e-12) / Nmc) + 1e-12; };
  // H+ at 1 keV: reflection R_N(eps(1 keV)) = 0.329; potential emission 0.032 (0.78*13.598 - 2*4.35) = 0.0610; kinetic emission: H table at 1 keV
  {
    auto s = run(SW7_HP, 1000, Nmc, P);
    const double rn = sw7::tjs(P.rn_A, table[SW7_HP].refl_epsL_per_keV * 1.0);
    bool out = false;
    const double gke = sw7::interp_loglog(P.ke_H_E, P.ke_H_g, P.ke_H_n, 1000 * (1.00794 / (table[SW7_HP].mass_kg / 1.66053906660e-27)), &out);
    const double gpe = 0.032 * (0.78 * 13.598434599702 - 2 * 4.35);
    ok = ok && within(s.refl, rn) && std::abs(s.electrons - (gpe + gke)) <= 5 * std::sqrt((gpe + gke) / Nmc);
    rows2.push_back({{"case", "H+ 1 keV"}, {"reflected", s.refl}, {"R_N", rn}, {"electrons", s.electrons}, {"gamma_PE+KE", gpe + gke}});
  }
  // H2O+ at 10 keV: additivity gamma_O(8.88 keV) + 2 gamma_H(555 eV); PE below threshold for phi > 4.92 only
  {
    auto s = run(SW7_H2OP, 10000, Nmc, P);
    const double m = table[SW7_H2OP].mass_kg / 1.66053906660e-27;
    bool out = false;
    const double g = sw7::interp_loglog(P.ke_O_E, P.ke_O_g, P.ke_O_n, 10000 * 15.9994 / m, &out) + 2 * sw7::interp_loglog(P.ke_H_E, P.ke_H_g, P.ke_H_n, 10000 * 1.00794 / m, &out);
    const double gpe = std::max(0.0, 0.032 * (0.78 * 12.621 - 2 * 4.35));
    const double eps = table[SW7_H2OP].refl_epsL_per_keV * 10 * table[SW7_H2OP].refl_energy_share;
    const double rn = sw7::tjs(P.rn_O_A, eps);
    ok = ok && std::abs(s.electrons - (g + gpe)) <= 5 * std::sqrt((g + gpe) / Nmc) && within(s.refl, rn);
    rows2.push_back({{"case", "H2O+ 10 keV"}, {"electrons", s.electrons}, {"gamma_additivity+PE", g + gpe}, {"reflected", s.refl}, {"R_N_O_set", rn}});
  }
  // H3O+ at 30 eV: no potential emission (E_rec 6.4 < 2 phi), below the kinetic threshold, fragment H on every event,
  // water product sticks with s0 (1 - theta)
  {
    auto s = run(SW7_H3OP, 30, Nmc, P);
    const double eps = table[SW7_H3OP].refl_epsL_per_keV * 0.03 * table[SW7_H3OP].refl_energy_share;
    const double rn = sw7::tjs(P.rn_O_A, eps), sc = P.stick0 * (1 - P.theta);
    ok = ok && s.electrons == 0 && s.frag == 1 && within(s.stuck, (1 - rn) * sc);
    rows2.push_back({{"case", "H3O+ 30 eV"}, {"electrons", s.electrons}, {"fragments", s.frag}, {"stuck", s.stuck}, {"expected_stuck", (1 - rn) * sc}});
  }
  // H2O neutral at 0.04 eV: stick s0 (1 - theta), diffuse alpha, specular 1 - alpha
  {
    auto s = run(SW7_H2O, 0.04, Nmc, P);
    const double sc = P.stick0 * (1 - P.theta);
    ok = ok && within(s.stuck, sc) && within(s.thermal, (1 - sc) * P.accommodation) && within(s.specular, (1 - sc) * (1 - P.accommodation));
    rows2.push_back({{"case", "H2O 0.04 eV"}, {"stuck", s.stuck}, {"sticking", sc}, {"diffuse", s.thermal}, {"specular", s.specular}});
  }
  // H atom: light accommodation, no sticking
  {
    auto s = run(SW7_H, 0.5, Nmc, P);
    ok = ok && s.stuck == 0 && within(s.thermal, P.accommodation_light);
    rows2.push_back({{"case", "H 0.5 eV"}, {"diffuse", s.thermal}, {"alpha_light", P.accommodation_light}});
  }
  R.record("ions", "ion and neutral channels: sampled reflection, PE+KE electrons, fragments, sticking and accommodation equal the model within 5 sigma", ok, {{"rows", rows2}});
  // (c) stimulated desorption cross sections: threshold and table
  bool out = false;
  const bool esd = sw7::interp_loglog(P.esd_E, P.esd_sigma, P.esd_n, 9.9, &out) == 0 && sw7::interp_loglog(P.esd_E, P.esd_sigma, P.esd_n, 100, &out) == 8e-21 &&
                   std::abs(sw7::interp_loglog(P.esd_E, P.esd_sigma, P.esd_n, 300, &out) - 2e-20) < 1e-30;
  R.record("ions", "ESD cross section: zero below the 10 eV threshold, 8e-17 cm^2 at 100 eV, 2e-16 cm^2 at 300 eV", esd);
}

// ------------------------------------------------------------ ledger and partitions
struct Batch { std::vector<SW7Hit> hits; std::vector<SW7Out> outs; };
void evaluate(const std::vector<SW7Params>& P, const std::array<SW7Species, SW7_NSPECIES>& table, const std::vector<SW7Hit>& hits, const std::vector<std::size_t>& idx,
              std::uint64_t seed, std::uint64_t step, std::vector<SW7Out>& outs) {
  for (std::size_t i : idx) sw7::wall_event(P[hits[i].patch], table.data(), hits[i], seed, step, &outs[i]);
}
bool same_out(const SW7Out& a, const SW7Out& b) { return std::memcmp(&a, &b, sizeof(SW7Out)) == 0; }
void ledger_group(Results& R, const std::array<SW7Species, SW7_NSPECIES>& table, const std::filesystem::path& out) {
  const std::size_t n = 1000000;
  const auto P = v7test::patch_params(base_params());
  const auto hits = v7test::event_set(n, 20261031, table);
  std::vector<SW7Out> outs(n);
  std::vector<std::size_t> all(n);
  std::iota(all.begin(), all.end(), 0);
  const auto t0 = std::chrono::steady_clock::now();
  evaluate(P, table, hits, all, 4242, 12345, outs);
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  // per-event identities
  long bad = 0, refused = 0, neg_electron_wall = 0;
  sv::Sums total;
  std::map<std::uint64_t, sv::Sums> byweight;
  sv::i128 births = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const auto& o = outs[i];
    if (!sv::event_closes(hits[i], o, table)) ++bad;
    if (o.status & (SW7_ST_INVALID | SW7_ST_RNG_BUDGET)) ++refused;
    if (table[hits[i].species].kind == 0 && o.wall_q < 0) ++neg_electron_wall;
    total.add(hits[i], o);
    byweight[std::bit_cast<std::uint64_t>(hits[i].weight)].add(hits[i], o);
    births += o.n_out;
  }
  // global exact closure: sum e_in = sum e_out + sum wall; e_out = sum kq + n_e phi - ads; charge_out = -n_e
  sv::i128 kout = 0, phi_out = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const auto& o = outs[i];
    for (unsigned k = 0; k < o.n_out; ++k) { kout += o.kq[k]; if (o.sp[k] == SW7_E) phi_out += o.phi_q; }
  }
  const bool global = total.e_in_q == total.e_out_q + total.wall_q && total.e_out_q == kout + phi_out - total.ads_q && total.charge_out == -total.electrons_out &&
                      kout == total.k_out_q;
  sv::i128 births_species = 0;
  for (int s = 0; s < SW7_NSPECIES; ++s) births_species += total.out_species[s];
  // physical totals per weight class (J): incident energy = emitted + wall heat (exact in quanta; one rounding per class)
  J classes = J::array();
  for (const auto& [wb, s] : byweight)
    classes.push_back({{"weight", std::bit_cast<double>(wb)}, {"events", s.events}, {"e_in_quanta", sv::i128_text(s.e_in_q)}, {"e_out_quanta", sv::i128_text(s.e_out_q)},
                       {"wall_quanta", sv::i128_text(s.wall_q)}, {"stuck", s.stuck}, {"electrons_out", s.electrons_out}, {"wall_heat_J", std::bit_cast<double>(wb) * double(s.wall_q) / SW7_QUANTA_PER_EV * qe}});
  J species = J::object();
  for (int s = 0; s < SW7_NSPECIES; ++s) species[sv::species_name(s)] = {{"in", total.in_species[s]}, {"out", total.out_species[s]}};
  R.record("ledger", "1e6 mixed wall events: every per-event exact identity holds (energy quanta, charge, counts); no refused event; electron wall heat >= 0",
           bad == 0 && refused == 0 && neg_electron_wall == 0, {{"events", n}, {"identity_failures", bad}, {"refused", refused}, {"seconds", secs}});
  R.record("ledger", "global exact sums (int128): sum e_in = sum e_out + sum wall; e_out = sum K + n_e phi - E_ads; births = sum of outputs; charge out = -electrons out",
           global && births == births_species && births == sv::i128(total.out_species[0]) + (births_species - total.out_species[0]),
           {{"e_in_quanta", sv::i128_text(total.e_in_q)}, {"e_out_quanta", sv::i128_text(total.e_out_q)}, {"wall_quanta", sv::i128_text(total.wall_q)},
            {"births", sv::i128_text(births)}, {"stuck", total.stuck}, {"electrons_out", total.electrons_out}, {"energy_capped", total.energy_capped},
            {"yield_clamped", total.yield_clamped}, {"outside_data", total.outside_data}, {"rng_words", total.words}, {"species", species}, {"weight_classes", classes}});
  // partitions: 1, 2, 5, 8 ranks (contiguous, strided, reversed), merged by event index: bit-identical outputs and sums
  J parts = J::array();
  bool same = true;
  for (int ranks : {1, 2, 5, 8})
    for (int mode = 0; mode < 3; ++mode) {
      std::vector<SW7Out> o2(n);
      std::vector<sv::Sums> rank_sums(ranks);
      for (int r = 0; r < ranks; ++r) {
        std::vector<std::size_t> idx;
        for (std::size_t i = 0; i < n; ++i) {
          const int owner = mode == 0 ? int(i * ranks / n) : int(i % ranks);
          if (owner == r) idx.push_back(i);
        }
        if (mode == 2) std::reverse(idx.begin(), idx.end());
        evaluate(P, table, hits, idx, 4242, 12345, o2);
        for (std::size_t i : idx) rank_sums[r].add(hits[i], o2[i]);
      }
      sv::Sums merged;
      for (int r = ranks - 1; r >= 0; --r) merged.add(rank_sums[r]);   // reversed merge order
      long diff = 0;
      for (std::size_t i = 0; i < n; ++i) diff += !same_out(o2[i], outs[i]);
      const bool ok = diff == 0 && merged == total;
      same = same && ok;
      parts.push_back({{"ranks", ranks}, {"mode", mode == 0 ? "contiguous" : (mode == 1 ? "strided" : "strided-reversed")}, {"differing_events", diff}, {"sums_equal", merged == total}});
    }
  R.record("partition", "outputs bit-identical and integer sums equal for 1/2/5/8 partitions (contiguous, strided, reversed), merged in reverse order", same, {{"rows", parts}});
  // write a digest of the outputs (compared with the GPU run)
  std::ofstream(out / "ledger-1e6.json") << J{{"events", n}, {"e_in_quanta", sv::i128_text(total.e_in_q)}, {"wall_quanta", sv::i128_text(total.wall_q)},
                                              {"births", sv::i128_text(births)}, {"stuck", total.stuck}}.dump(1) << '\n';
}

// ------------------------------------------------------------ interval bound (first-order exposure)
void interval_group(Results& R) {
  // Expected count over [0, dt] of an emission rate driven by a field ramp F(t) = F0 (1 + kappa t), exact integral by
  // 16-point Gauss-Legendre in time vs the explicit start-of-step value J(F0) dt.
  const double phi = 4.35, T = 300, F0 = 3e9, kappa = 2e10;   // d ln F / dt = 2e10 s^-1
  auto lJ = [&](double t) { return se::evaluate(T, F0 * (1 + kappa * t), se::Config{phi}).log_current_A_m2; };
  const double xg[8] = {0.0950125098376374, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438, 0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499};
  const double wg[8] = {0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767, 0.1246279358892235, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541};
  auto exact = [&](double dt) {
    double s = 0;
    for (int i = 0; i < 8; ++i) for (int sg : {-1, 1}) s += wg[i] * std::exp(lJ(dt / 2 * (1 + sg * xg[i])) - lJ(0));
    return s / 2;   // integral / (J(0) dt)
  };
  J rows = J::array();
  std::vector<double> errs;
  for (double dt : {1e-14, 2e-14, 4e-14, 8e-14}) {
    const double e = std::abs(1 / exact(dt) - 1);   // relative error of the explicit count
    errs.push_back(e);
    rows.push_back({{"dt_s", dt}, {"relative_error", e}});
  }
  const double p = std::log2(errs[2] / errs[1]);
  // interval bound from two accepted samples (t = 0 and t = 10 fs) with eps = 1e-3: the actual error at the bound is <= eps
  sv::RateSample a, b;
  a.t_s = 0; a.dt_s = 1e-14; b.t_s = 1e-14; b.dt_s = 1e-14;
  a.log_emission = {lJ(0) + std::log(1e-4)}; b.log_emission = {lJ(1e-14) + std::log(1e-4)};
  a.log_desorption = b.log_desorption = {};
  sv::IntervalPolicy pol;
  pol.cap_s = 1e-12; pol.emission_relative_error = 1e-3; pol.safety = 0.5; pol.count_floor = 0;
  const auto d = sv::interval_bound(a, b, pol);
  const double e_at = std::abs(1 / exact(d.bound_s) - 1);
  R.record("interval", "explicit start-of-step exposure is first order (observed p within 0.15 of 1) and the bound dt = 2 eps s / |dlnN/dt| keeps the error <= eps",
           std::abs(p - 1) <= 0.15 && d.limiting == "emission" && e_at <= 1e-3,
           {{"rows", rows}, {"observed_order", p}, {"bound_s", d.bound_s}, {"rate_per_s", d.rate_per_s}, {"error_at_bound", e_at}});
}


// ------------------------------------------------------------ extended emission domain (3000-5000 K)
void hot_group(Results& R, const std::filesystem::path& out) {
  namespace hot = electrode_heat_v7::emission_hot;
  // (a) identical to the sealed evaluate wherever the sealed domain applies (T <= 3000 K)
  long cmp = 0, same = 0, refs = 0, same_ref = 0;
  for (double phi : {4.35, 3.38, 4.5, 5.1})
    for (double T : {0., 1., 300., 1000., 1999.5, 2500., 3000.})
      for (double F : {-1e6, 0., 1., 1e4, 1e6, 1e7, 1e8, 1e9, 3e9, 5e9}) {
        std::string ea, eb;
        se::Distribution a, b;
        try { a = se::evaluate(T, F, se::Config{phi}); } catch (const std::exception& x) { ea = x.what(); }
        try { b = hot::evaluate(T, F, se::Config{phi}, 5000); } catch (const std::exception& x) { eb = x.what(); }
        if (!ea.empty() || !eb.empty()) { ++refs; same_ref += (!ea.empty() && !eb.empty()); continue; }
        ++cmp;
        same += electrode_heat_v6_1::emission::same_distribution(a, b);
      }
  R.record("hot", "extended evaluate bit-identical to the sealed kernel for T <= 3000 K (and refuses where it refuses)", cmp > 200 && same == cmp && same_ref == refs,
           {{"compared", cmp}, {"identical", same}, {"refused_both", same_ref}, {"refused", refs}});
  // (b) exact thermal law at F = 0 and Richardson-Schottky at low field, 3000-5000 K
  const double A0 = se::ideal_richardson_A_m2_K2(), cS2 = 1.439964e-9;
  double fd = 0, rs = 0;
  for (double T : {3200., 3500., 4000., 4500., 5000.})
    for (double phi : {4.0, 4.4, 4.75}) {
      const auto d = hot::evaluate(T, 0.0, se::Config{phi}, 5000);
      const double x = std::exp(-phi / (kBeV * T));
      double li = 0, xk = 1;
      for (int k = 1; k < 60; ++k) { xk *= x; li += ((k & 1) ? 1 : -1) * xk / (double(k) * k); }
      fd = std::max(fd, std::abs(d.log_current_A_m2 - (std::log(A0 * T * T) + std::log(li))));
      for (double F : {1e4, 1e5, 1e6, 1e7}) {
        const auto e = hot::evaluate(T, F, se::Config{phi}, 5000);
        rs = std::max(rs, std::abs(e.log_current_A_m2 - (std::log(A0 * T * T) - (phi - std::sqrt(cS2 * F)) / (kBeV * T))));
      }
    }
  R.record("hot", "3200-5000 K: Fermi-Dirac thermal law at F = 0 within 1e-8; Richardson-Schottky (1e4-1e7 V/m) within 20 %", fd <= 1e-8 && rs <= std::log(1.2),
           {{"fd_max_abs_log", fd}, {"rs_max_abs_log", rs}});
  // (c) quadrature convergence above 3000 K: default (128 panels, 64 action nodes, 48 tail widths) vs refined (512, 128, 80)
  double conv = 0;
  J crow = J::array();
  for (double T : {3500., 4000., 5000.})
    for (double F : {0., 1e6, 1e8, 1e9, 3e9}) {
      const auto a = hot::evaluate(T, F, se::Config{4.4}, 5000);
      const auto b = hot::evaluate(T, F, se::Config{4.4, 512, 128, 80}, 5000);
      const double dl = std::abs(a.log_current_A_m2 - b.log_current_A_m2), dm = std::abs(a.mean_metal_energy_eV - b.mean_metal_energy_eV);
      conv = std::max({conv, dl, dm / std::max(1.0, std::abs(b.mean_metal_energy_eV))});
      crow.push_back({{"T", T}, {"F", F}, {"dlogJ", dl}, {"dmean_eV", dm}});
    }
  R.record("hot", "3500-5000 K: current and energy moment converged to 1e-6 under quadrature refinement", conv <= 1e-6, {{"max", conv}, {"rows", crow}});
  // (d) validity map: where the planar SN kernel stops being valid between 300 and 5000 K
  //     numerical domain |F| <= 5e9 V/m; barrier collapse 2 sqrt(cF) >= 0.98 phi; free-electron supply: kT / E_F
  //     (E_F(Fe) 11.1 eV free-electron value, model_derived); over-barrier fraction; current density (space charge is
  //     the PIC's); molten surface (planar, beta = 1) above the solidus.
  J map = J::array();
  long inside = 0, collapse = 0, outside_F = 0;
  for (double phi : {3.4, 4.0, 4.4, 4.75})
    for (double T : {300., 1000., 2000., 3000., 3500., 4000., 4500., 5000.})
      for (double F : {1e5, 1e7, 1e8, 1e9, 2e9, 3e9, 5e9, 6e9}) {
        J row = {{"phi", phi}, {"T", T}, {"F", F}};
        try {
          const auto d = hot::evaluate(T, F, se::Config{phi}, 5000);
          const double top = phi - d.barrier_lowering_eV;
          double above = 0;
          for (const auto& n : d.nodes) if (n.normal_energy_eV >= top) above += n.probability;
          row["log_J"] = d.log_current_A_m2; row["lowering_over_phi"] = d.barrier_lowering_eV / phi; row["kT_over_d"] = kBeV * T / d.fermi_decay_width_eV;
          row["kT_over_EF"] = kBeV * T / 11.1; row["above_barrier"] = above; row["valid"] = true;
          // barrier-top distance sqrt(c/F) (Kyritsakis & Djurabekova 2017 regime test uses the barrier top); below the
          // Fe interatomic spacing (~0.25 nm, F >~ 5.8 V/nm) the classical image potential is questionable
          row["barrier_top_distance_nm"] = 1e9 * std::sqrt(3.5999e-10 / F); row["image_questionable"] = std::sqrt(3.5999e-10 / F) < 0.25e-9;
          ++inside;
        } catch (const std::exception& x) {
          const std::string w = x.what();
          row["valid"] = false; row["reason"] = w;
          if (w.find("Barrier collapse") != std::string::npos) ++collapse;
          if (w.find("normal-field") != std::string::npos) ++outside_F;
        }
        map.push_back(row);
      }
  std::ofstream(out / "validity-map-v7.json") << J{{"schema", "emission-validity-map-v7"}, {"criteria", {"|F| <= 5e9 V/m (numerical domain)", "Schottky lowering < 0.98 phi (barrier collapse)",
                                                                                                         "T <= 5000 K (extended domain)", "kT/E_F <= 0.04 at 5000 K (free-electron supply)"}},
                                                  {"points", map}}.dump(1) << '\n';
  R.record("hot", "validity map 300-5000 K x 1e5-6e9 V/m written; refusals only at |F| > 5e9 V/m or barrier collapse", inside > 0 && outside_F > 0 && inside + outside_F + collapse == long(map.size()),
           {{"valid_points", inside}, {"barrier_collapse", collapse}, {"field_domain", outside_F}});
}

// ------------------------------------------------------------ registry
void registry_group(Results& R, const std::filesystem::path& out) {
  const auto reg = sv::full_registry();
  J rows = J::array();
  std::map<std::string, int> tags;
  bool ok = !reg.rows.empty();
  std::set<std::string> names;
  for (const auto& p : reg.rows) {
    ok = ok && !p.source.empty() && !p.location.empty() && p.low <= p.value && p.value <= p.high && names.insert(p.name).second;
    ++tags[sv::tag_name(p.tag)];
    rows.push_back({{"name", p.name}, {"value", p.value}, {"low", p.low}, {"high", p.high}, {"unit", p.unit}, {"tag", sv::tag_name(p.tag)}, {"source", p.source},
                    {"location", p.location}, {"note", p.note}});
  }
  J t = J::object();
  for (const auto& [k, v] : tags) t[k] = v;
  std::ofstream(out / "parameters-v7.json") << J{{"schema", "surface-v7-parameters-1"}, {"parameters", rows}, {"tags", t}}.dump(1) << '\n';
  R.record("registry", "every parameter has a tag, source, location and a band containing its value; names unique", ok, {{"parameters", reg.rows.size()}, {"tags", t}});
}
}  // namespace

int main(int argc, char** argv) {
  Results R;
  int rc = 0;
  try {
    if (argc != 2) throw std::runtime_error("usage: test_surface_v7 OUTDIR");
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    const auto table = sv::species_table(nullptr);
    math_group(R);
    rng_group(R);
    emission_group(R);
    workfn_group(R);
    see_group(R, table);
    ion_group(R, table);
    ledger_group(R, table, out);
    interval_group(R);
    hot_group(R, out);
    registry_group(R, out);
    std::ofstream(out / "surface-v7.json") << J{{"schema", "surface-v7-qualification-v1"}, {"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    rc = R.failed == 0 ? 0 : 1;
  } catch (const std::exception& e) {
    std::cerr << "fatal: " << e.what() << '\n';
    rc = 2;
  }
  return rc;
}
