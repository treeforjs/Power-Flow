#include "parent_coverage_v2_subflow.hpp"
#include "parent_coverage_subflow.hpp" // coverage v1, for parity comparisons
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

// Coverage-scope v2 checks for the water parent-ionization channel: formulas,
// kernel conservation, classification, ledger categories, refusal limit and
// parity with v1. Conditional parent-only operator at a frozen particle stage;
// no physical trajectory time advances and no Hall coupling is constructed.
namespace nw = raikou::native_water;
namespace cv = raikou::native_water::coverage_v2;
namespace c1 = raikou::native_water::coverage_v1;
namespace bb = raikou::native_water::beb_v2;
namespace m = raikou::mfem_backend;
namespace pc = raikou::pic_collisions;
namespace w = raikou::pic_water;
using J = nw::Json;
J checks = J::array();
void check(bool ok, const std::string& name) {
  checks.push_back({{"name", name}, {"passed", ok}});
  if (!ok) throw std::runtime_error("check failed: " + name);
}
template<class F> void reject(F f, const std::string& s) {
  bool threw = false;try { f(); } catch (const std::exception&) { threw = true; }check(threw, s);
}
template<class F> cv::CoverageRefusal refusal(F f, const std::string& s) {
  try { f(); } catch (const cv::CoverageRefusal& r) { check(true, s);return r; } catch (const std::exception& e) {
    check(false, s + " (other error: " + e.what() + ")");
  }
  check(false, s + " (no refusal)");throw std::logic_error("unreachable");
}
constexpr double mn = 2.9915e-26;
const double inf = std::numeric_limits<double>::infinity();
const std::string outer3_sha = "c436020d3f2d07902fb059cfc0a269b2e914a3cbb2789300093d3631df896551";
m::GapPicOptions options() {
  m::GapPicOptions o;o.nx = o.ny = 45;o.width_m = o.gap_m = .009;o.depth_m = .025;
  o.finite_foils = true;o.finite_axial_extent = true;o.short_z_m = 0;o.feed_z_m = .025;
  o.lower_foil = {.002, .007, .0018, .002};o.upper_foil = {.002, .007, .007, .0072};return o;
}
m::PicVector electron_u(double E) { double x = E * pc::qe / (w::me * pc::c * pc::c);return {pc::c * std::sqrt(x * (x + 2)), 0, 0}; }
std::vector<m::ArenaStateSpec> states() {
  return {{"electron", 0, 0, -1, w::me, 0, 0}, {"H2O", 2, 1, 0, mn, 0, 0},
          {"H2O+", 2, 1, 1, mn - w::me + w::threshold_eV * pc::qe / (pc::c * pc::c), 0, w::threshold_eV * pc::qe}};
}
std::vector<m::ArenaChannel> channels(const cv::ArenaBinding& a) {
  return {{w::parent_channel, "electron", "H2O", {"electron", "electron", "H2O+"}, a.dataset_sha256, a.kernel_manifest_sha256,
           a.kernel_binding_file_sha256, a.model}};
}
const m::PicVector here{.004, .004, .0125}, water_at{.00403, .004, .0125}, elsewhere{.005, .004, .0125};
struct Electron { double weight, energy_eV;m::PicVector position; };
struct Fixture {
  m::GapPic pic{options()};std::unique_ptr<m::PicArena> arena;
  Fixture(const cv::ArenaBinding& a, std::vector<Electron> lead, double water_weight, std::vector<Electron> extra = {}) {
    std::vector<m::GapParticle> batch;
    for (const auto& x : lead) { m::GapParticle e;e.weight = x.weight;e.mass_kg = w::me;e.position_m = x.position;e.proper_velocity_m_s = electron_u(x.energy_eV);batch.push_back(e); }
    m::GapParticle n;n.weight = water_weight;n.mass_kg = mn;n.charge_C = 0;n.position_m = water_at;batch.push_back(n);
    for (const auto& x : extra) { m::GapParticle e;e.weight = x.weight;e.mass_kg = w::me;e.position_m = x.position;e.proper_velocity_m_s = electron_u(x.energy_eV);batch.push_back(e); }
    auto ids = pic.add_authored_batch(batch);pic.solve_fields();
    std::vector<m::ArenaPacket> packets;std::size_t k = 0;
    for (std::size_t i = 0; i < lead.size(); ++i, ++k)
      packets.push_back({i ? "e" + std::to_string(i) : "e", "electron", "authored_fixture", "", {lead[i].weight, 0}, {ids[k]}});
    packets.push_back({"n", "H2O", "authored_fixture", "", {water_weight, 0}, {ids[k++]}});
    for (std::size_t i = 0; i < extra.size(); ++i, ++k)
      packets.push_back({"x" + std::to_string(i), "electron", "authored_fixture", "", {extra[i].weight, 0}, {ids[k]}});
    arena = std::make_unique<m::PicArena>(pic, states(), channels(a), packets, "native-water-subflow-fixture/v1");
  }
};
template<class Flow>
nw::Interval interval(m::PicArena& arena, const Flow& flow, double duration, std::uint64_t sequence = 0) {
  nw::Interval i;i.sequence = sequence;i.end_s = duration;i.common_clock_sha256 = nw::sha("frozen-native-parent-fixture/" + std::to_string(sequence));
  i.phase.state_prefix = arena.state_prefix();i.phase.mesh_prefix = flow.grid_identity();i.phase.clock_prefix = i.common_clock_sha256;return i;
}
long double weight(const m::PicArena& a, const std::string& state) {
  long double total = 0;
  for (const auto& [_, p] : a.snapshot()->packets) if (p.state == state) for (double part : m::PicArena::positive_carriers(p.count)) total += part;
  return total;
}
J events_without_native(const J& audit) { J out = J::array();for (auto e : audit.at("events")) { e.erase("native");out.push_back(e); }return out; }
J rng_fields(const J& cp) { return {cp.at("scheduler"), cp.at("product_rng"), cp.at("next_event"), cp.at("next_interval")}; }
J census(const m::PicArena& a, const std::set<std::string>& moved = {}) {
  J out = J::array();
  for (const auto& r : a.rate_population(100000)) {
    J row{{"packet", r.packet_id}, {"state", r.state}, {"id", r.native.id}, {"weight", r.native.weight}, {"mass", r.native.mass_kg},
          {"u", r.native.proper_velocity_m_s}};
    if (!moved.count(r.packet_id)) row["x"] = r.native.position_m;
    out.push_back(row);
  }
  return out;
}
double value(const J& exact) { return exact.at("value").get<double>(); }
template<class P> P record_only() { P p;p.max_exposure_unresolved_events = inf;p.max_cumulative_unresolved_events = inf;p.max_unresolved_reaction_fraction = 1;p.admit_unbounded = true;return p; }
std::string hex(double x) { char b[64];std::snprintf(b, sizeof b, "%a", x);return b; }
const J& category(const J& coverage, cv::Category c) { return coverage.at("categories").at(static_cast<std::size_t>(c)); }
// Electron kinetic energy from a proper velocity.
long double kinetic_eV(const m::PicVector& u) { return w::kinetic(w::me, u) / pc::qe; }

// Composite Gauss-Legendre (8 points) of the reduced SDCS on log(1+w) panels.
long double quadrature(const bb::Shape& s, int panels) {
  static const long double x8[4] = {0.18343464249564980493947614236018398L, 0.52553240991632898581773904918924635L,
                                    0.79666647741362673959155393647583044L, 0.96028985649753623168356086856947299L};
  static const long double w8[4] = {0.36268378337836198296515044927719561L, 0.31370664587788728733796220198660131L,
                                    0.22238103445337447054435599442624088L, 0.10122853629037625915253135430996219L};
  const long double top = std::log1p(s.w_max());long double sum = 0;
  for (int k = 0; k < panels; ++k) {
    const long double a = top * k / panels, b = top * (k + 1) / panels, mid = (a + b) / 2, half = (b - a) / 2;
    for (int i = 0; i < 4; ++i)
      for (int sign : {-1, 1}) {
        const long double x = mid + sign * half * x8[i], wv = std::expm1(x);
        sum += half * w8[i] * s.g(wv) * (1 + wv);
      }
  }
  return sum;
}

int main(int argc, char** argv) {
  try {
    if (argc != 6) throw std::runtime_error("dataset binding outer3-dataset kernel-binding output-directory");
    std::ifstream f(argv[2]);J binding;f >> binding;
    nw::Binding b{binding.at("dataset_sha256"), binding.at("kernel_sha256"), binding.at("kernel_binding_file_sha256")};
    const std::string outer3_path = argv[3], kernel_binding_path = argv[4], dir = argv[5];
    auto data = std::make_shared<const w::Data>(w::Data::read(argv[1]));
    w::ParentBridge bridge(data, w::parent_model, b.dataset_sha256, b.kernel_sha256, mn);
    nw::Limits limits;nw::Grid grid(options(), {45, 45, 5}, limits.cells);
    const double V = grid.cells[0].volume_m3;
    auto make_bands = [&](cv::ProvisionalOptions o = {}) {
      return std::make_shared<const cv::HighEnergyBands>(cv::HighEnergyBands::read(outer3_path, outer3_sha, *data, b.dataset_sha256,
          kernel_binding_path, b.kernel_sha256, b.kernel_binding_file_sha256, o));
    };
    const auto bands = make_bands();
    const auto& he = *bands;const auto& K = he.constants();const auto& orb = he.orbitals();
    J receipts = J::array(), report;

    // ======================= A. formulas =======================
    {
      // A1. Dataset arithmetic: the unweighted v2 evaluation is the bound v1
      // H2O:nist_BEB_valence value bit for bit.
      bool same = true;
      for (double E : {0., 12.61, 12.62, 13., 15.57, 20., 50., 100., 999.999, 1000., 1000.5, 5000., 8084.470604111418, 9999.999, 10000.})
        same = same && he.valence_beb_m2(E) == data->observable_value_SI("H2O:nist_BEB_valence", E);
      check(same, "dataset BEB arithmetic reproduced bit for bit (bound H2O:nist_BEB_valence)");
      bool weights_ok = true;
      for (double E : {1000., 1500., 3000., 8000., 10000.}) {
        const double s = he.outer3_parent_m2(E);
        const long double parts = bb::beb_orbital_m2(orb[3], K, E) + bb::beb_orbital_m2(orb[2], K, E) + .08L * bb::beb_orbital_m2(orb[1], K, E);
        weights_ok = weights_ok && std::abs(s - parts) <= 8 * std::numeric_limits<double>::epsilon() * s && s < he.valence_beb_m2(E);
      }
      check(weights_ok, "outer3 parent rate = 1.00 s(1b1) + 1.00 s(3a1) + 0.08 s(1b2), 2a1 excluded (Murakami 2012 Eq. 27)");
      check(he.weights() == std::vector<double>{0, .08, 1, 1} && orb[0].label == "2a1" && orb[3].label == "1b1" &&
                orb[3].B_eV == 12.61 && orb[3].U_eV == 61.91 && orb[2].B_eV == 15.57 && orb[1].B_eV == 19.83 && orb[0].B_eV == 36.88 &&
                K.bohr_radius_m == 5.2918e-11 && K.rydberg_eV == 13.6057 && K.rest_energy_eV == 510998.95069,
            "orbital constants, parent weights and RBEB constants bound from the frozen dataset");
      // A2. Threshold limits.
      bool thr = true;
      for (const auto& o : orb) {
        thr = thr && bb::rbeb_orbital_m2(o, K, o.B_eV) == 0 && bb::beb_orbital_m2(o, K, o.B_eV) == 0 &&
              bb::rbeb_orbital_m2(o, K, o.B_eV * (1 - 1e-12)) == 0 && bb::rbeb_orbital_m2(o, K, o.B_eV * (1 + 1e-9)) > 0;
        for (bool rel : {false, true}) {
          auto s = [&](long double d) { const long double T = o.B_eV * (1 + d);return (rel ? bb::rbeb_orbital_m2(o, K, T) : bb::beb_orbital_ext_m2(o, K, T)) / d; };
          thr = thr && std::abs(s(1e-6L) / s(1e-7L) - 1) < 1e-5L && s(1e-7L) > 0;
        }
      }
      check(thr, "threshold: zero at and below B, linear positive onset (BEB and RBEB)");
      // A3. Nonrelativistic recovery.
      bb::Constants big = K;big.rest_energy_eV = K.rest_energy_eV * 1e8;
      long double worst_nr = 0;
      for (const auto& o : orb)
        for (double T : {13., 20., 50., 100., 1000., 5000., 10000.})
          if (T > o.B_eV) worst_nr = std::max(worst_nr, std::abs(bb::rbeb_orbital_m2(o, big, T) / bb::beb_orbital_ext_m2(o, K, T) - 1));
      check(worst_nr < 1e-9L, "RBEB with mc^2 -> 10^8 mc^2 reproduces BEB (Q=1) to 1e-9 (13 eV-10 keV)");
      bool small = true;J nr_rows = J::array();
      for (double T : {100., 1000., 2000., 5000., 10000.}) {
        const long double r = he.rbeb_parent_any(T) / he.beb_parent_any(T) - 1, x = T / K.rest_energy_eV;
        small = small && r > 0 && r / x > .5L && r / x < 2;
        nr_rows.push_back({{"energy_eV", T}, {"relative_correction", static_cast<double>(r)}, {"per_T_over_mc2", static_cast<double>(r / x)}});
      }
      check(small, "relativistic correction vanishes linearly in T/mc^2 at 0.1-10 keV");
      bool shape_ok = true;
      for (const auto& o : orb) for (double T : {50., 1000., 8000.}) {
        const auto s = bb::beb_shape(o, K, T);
        shape_ok = shape_ok && std::abs(s.total_m2() / bb::beb_orbital_m2(o, K, T) - 1) < 1e-15L;
      }
      check(shape_ok, "BEB SDCS integral equals the dataset BEB orbital term to 1e-15");
      // A4. High-energy Bethe-like behavior.
      long double prev = 1;bool asym = true;J asym_rows = J::array();
      for (double T : {1e5, 1e6, 1e7, 1e8, 1e9}) {
        long double worst = 0;
        for (std::size_t i = 1; i < 4; ++i) {
          const auto s = bb::rbeb_shape(orb[i], K, T);
          worst = std::max(worst, std::abs(bb::rbeb_orbital_m2(orb[i], K, T) / (s.P * (s.L / 2 + 1)) - 1));
        }
        asym = asym && worst < prev && (T < 1e6 || worst < 1e-4L);prev = worst;
        asym_rows.push_back({{"energy_eV", T}, {"relative_residual", static_cast<double>(worst)}});
      }
      check(asym, "RBEB approaches the relativistic Bethe form P [ (ln(b^2/(1-b^2)) - b^2 - ln 2b')/2 + 1 ] monotonically");
      bool fano = true;
      for (std::size_t i = 1; i < 4; ++i) {
        auto XY = [&](long double T) {
          const auto s = bb::rbeb_shape(orb[i], K, T);const long double tp = T / K.rest_energy_eV;
          return std::pair<long double, long double>{std::log(tp * (tp + 2)) - bb::beta2(tp), bb::rbeb_orbital_m2(orb[i], K, T) / s.P};
        };
        const auto a = XY(1e8), c = XY(1e9);
        fano = fano && std::abs((c.second - a.second) / (c.first - a.first) - .5L) < 1e-3L;
      }
      check(fano, "Fano-plot slope d[sigma/P]/d[ln(b^2/(1-b^2)) - b^2] -> 1/2 at 0.1-1 GeV");
      bool nrb = true;long double prevb = 1;
      for (double T : {1e5, 1e6, 1e7}) {
        long double worst = 0;
        for (std::size_t i = 1; i < 4; ++i) {
          const auto s = bb::beb_shape(orb[i], K, T);
          worst = std::max(worst, std::abs(bb::beb_orbital_ext_m2(orb[i], K, T) / (s.P * (s.L / 2 + 1)) - 1));
        }
        nrb = nrb && worst < prevb;prevb = worst;
      }
      check(nrb && prevb < 1e-4L, "nonrelativistic BEB approaches the Bethe form S/(t+u+1) (ln t/2 + 1)");
      const long double s1 = he.rbeb_parent_any(1e6), s2 = he.rbeb_parent_any(2e6), s5 = he.rbeb_parent_any(5e6), s8 = he.rbeb_parent_any(1e8);
      check(s2 < s1 && s2 < s5 && s5 < s8, "RBEB minimum between 1 and 5 MeV followed by the relativistic rise");
      report["formula_limits"] = {{"nonrelativistic_recovery_worst_relative", static_cast<double>(worst_nr)}, {"relativistic_correction_low_energy", nr_rows},
                                  {"bethe_asymptote_residual", asym_rows}};
      // A5. SDCS.
      bool sd = true, sym = true, pos = true, inv = true;long double worst_q = 0;
      for (double T : {1000.5, 5000., 10000., 10000.5, 1e5, 1e6})
        for (std::size_t i = 1; i < 4; ++i) {
          const bool rel = T > 10000;
          const auto s = rel ? bb::rbeb_shape(orb[i], K, T) : bb::beb_shape(orb[i], K, T);
          const long double closed = rel ? bb::rbeb_orbital_m2(orb[i], K, T) : bb::beb_orbital_ext_m2(orb[i], K, T);
          sd = sd && std::abs(s.total_m2() / closed - 1) < 1e-15L;
          worst_q = std::max(worst_q, std::abs(quadrature(s, 400) / s.G(s.w_max()) - 1));
          for (int k = 0; k <= 1000; ++k) {
            const long double wv = s.w_max() * k / 1000;
            pos = pos && s.g(wv) > 0;
            sym = sym && std::abs(s.g(wv) / s.g(s.t - 1 - wv) - 1) < 1e-15L;
          }
          for (int k = 1; k < 100; ++k) {
            const long double u = k / 100.L, x = bb::sample_reduced_w(s, u);
            inv = inv && std::abs(s.G(x) / s.G(s.w_max()) - u) < 1e-15L && x >= 0 && x <= s.w_max();
          }
        }
      check(sd, "SDCS analytic integral over 0..(T-B)/2 equals the closed-form BEB/RBEB total to 1e-15");
      check(worst_q < 1e-12L, "SDCS numerical quadrature equals its analytic integral to 1e-12");
      check(pos && sym, "SDCS positive and symmetric under exchange of the two outgoing electrons");
      check(inv, "SDCS inverse-CDF sampling inverts the analytic integral to 1e-15");
      std::mt19937_64 rng(2026);
      auto uni = [&] { return (static_cast<long double>(rng() >> 11) + .5L) * 0x1p-53L; };
      bool ks = true;J ks_rows = J::array();
      for (const auto& [T, rel] : std::vector<std::pair<double, bool>>{{8000., false}, {1e6, true}}) {
        const auto s = rel ? bb::rbeb_shape(orb[3], K, T) : bb::beb_shape(orb[3], K, T);
        const int n = 20000;std::vector<long double> x(n);
        for (auto& v : x) v = bb::sample_reduced_w(s, uni());
        std::sort(x.begin(), x.end());
        long double D = 0;const long double Gm = s.G(s.w_max());
        for (int k = 0; k < n; ++k) { const long double F = s.G(x[k]) / Gm;D = std::max({D, std::abs(F - (long double)k / n), std::abs(F - (long double)(k + 1) / n)}); }
        ks = ks && D < 1.628L / std::sqrt((long double)n);
        ks_rows.push_back({{"energy_eV", T}, {"orbital", "1b1"}, {"ks_statistic", static_cast<double>(D)}, {"critical_99", 1.628 / std::sqrt(n)}});
      }
      check(ks, "sampled ejected energies pass a Kolmogorov-Smirnov test against the analytic SDCS (99%)");
      bool chi = true;J chi_rows = J::array();
      for (double T : {5000., 1e5}) {
        const auto terms = T <= 10000 ? he.outer3_terms(T) : he.provisional_terms(T);
        long double tot = 0;for (auto x : terms) tot += x;
        std::array<int, 4> counts{};const int n = 30000;
        for (int k = 0; k < n; ++k) ++counts[bb::select_orbital(terms, uni())];
        long double x2 = 0;
        for (int i = 0; i < 4; ++i) {
          if (terms[i] == 0) { chi = chi && counts[i] == 0;continue; }
          const long double e = n * terms[i] / tot;x2 += (counts[i] - e) * (counts[i] - e) / e;
        }
        chi = chi && x2 < 13.8L;
        chi_rows.push_back({{"energy_eV", T}, {"counts_2a1_1b2_3a1_1b1", counts}, {"chi2_2dof", static_cast<double>(x2)}});
      }
      check(chi, "orbital selection frequencies follow w_i sigma_i (chi-square, 2 dof, 99.9%); 2a1 never selected");
      report["sdcs_sampling_tests"] = {{"ks", ks_rows}, {"orbital_chi2", chi_rows}, {"quadrature_worst_relative", static_cast<double>(worst_q)}};
    }

    // ======================= B. kernel conservation =======================
    {
      std::mt19937_64 rng(77);
      auto uni = [&] { return (static_cast<long double>(rng() >> 11) + .5L) * 0x1p-53L; };
      J rows = J::array();
      struct Case { const char* label;double T;cv::Band band; };
      for (const auto& c : {Case{"1 keV (outer3 band closed edge)", 1000., cv::Band::outer3_beb}, Case{"10 keV (outer3)", 10000., cv::Band::outer3_beb},
                            Case{"10 keV + 1 ulp (provisional)", std::nextafter(10000., 2e4), cv::Band::provisional_rbeb},
                            Case{"100 keV (provisional)", 1e5, cv::Band::provisional_rbeb}, Case{"1 MeV (provisional)", 1e6, cv::Band::provisional_rbeb}}) {
        const std::string L = c.label;
        bool energy = true, momentum = true, range = true, order = true, shell = true;
        long double worst_e = 0, worst_p = 0, max_recoil = 0, sumW = 0, max_defl = 0;
        const int n = 4000;
        std::array<double, 3> axis{0.48, -0.6, 0.64};
        for (int k = 0; k < n; ++k) {
          const auto p = cv::beb_sdcs_parent_products(he, c.band, c.T, mn, axis, uni(), uni(), uni() * .999999L);
          const long double T = c.T, p0 = std::sqrt(T * (T + 2 * 510998.95L));
          worst_e = std::max(worst_e, std::abs(p.energy_defect_eV) / T);worst_p = std::max(worst_p, p.momentum_defect_eV_c / p0);
          energy = energy && std::abs(p.energy_defect_eV) <= 1e-16L * T;momentum = momentum && p.momentum_defect_eV_c <= 1e-16L * p0;
          const double B = orb[p.orbital].B_eV;
          range = range && p.ejected_eV >= 0 && p.ejected_eV <= (c.T - B) / 2 * (1 + 1e-15) && p.orbital >= 1;
          order = order && p.primary_kinetic_J > p.secondary_kinetic_J;
          const long double recoil = p.ion_kinetic_J / pc::qe;max_recoil = std::max(max_recoil, recoil);
          for (const auto& [mom, mass, kin] : {std::tuple{p.primary_momentum_kg_m_s, w::me, p.primary_kinetic_J},
                                               std::tuple{p.secondary_momentum_kg_m_s, w::me, p.secondary_kinetic_J},
                                               std::tuple{p.ion_momentum_kg_m_s, p.ion_mass_kg, p.ion_kinetic_J}}) {
            w::Vec u{};for (int d = 0; d < 3; ++d) u[d] = mom[d] / mass;
            const long double k2 = w::kinetic(mass, u);
            shell = shell && std::abs(k2 - kin) <= 1e-12L * kin + 1e-40L;
          }
          sumW += p.ejected_eV;
          long double dot = 0, pn = 0;for (int d = 0; d < 3; ++d) { dot += p.primary_momentum_kg_m_s[d] * axis[d];pn += p.primary_momentum_kg_m_s[d] * p.primary_momentum_kg_m_s[d]; }
          max_defl = std::max(max_defl, std::acos(std::min(1.L, dot / std::sqrt(pn))));
        }
        check(energy && momentum, L + ": exact relativistic energy (with I) and momentum closure in the target rest frame");
        check(range && order, L + ": ejected energy within 0..(T-B)/2 from orbitals 1b2/3a1/1b1; scattered electron faster");
        check(shell && max_recoil < 1e-4L, L + ": products on shell; ground-ion recoil below 1e-4 eV");
        rows.push_back({{"case", L}, {"energy_eV", c.T}, {"band", cv::name(c.band)}, {"samples", n},
                        {"worst_relative_energy_defect", static_cast<double>(worst_e)}, {"worst_relative_momentum_defect", static_cast<double>(worst_p)},
                        {"max_ion_recoil_eV", static_cast<double>(max_recoil)}, {"mean_ejected_eV", static_cast<double>(sumW / n)},
                        {"equal_pair_share_eV_for_comparison", (c.T - w::threshold_eV) / 2}, {"max_scattered_deflection_rad", static_cast<double>(max_defl)}});
      }
      // Through the bridge with moving targets: laboratory closure is enforced
      // by the bridge (same 256-epsilon receipts as the sealed v1 bridge).
      std::uint64_t ignored = 0;
      const cv::CoverageScope scope(data, b.dataset_sha256, bridge, bands);
      const auto composite = scope.channel(bridge, &ignored, b.kernel_binding_file_sha256);
      bool lab = true;int tried = 0;
      for (double T : {1100., 9000., 11000., 1e5, 1e6})
        for (double speed : {0., 600., 2e5})
          for (int k = 0; k < 40; ++k) {
            const auto ax = std::array<double, 3>{std::cos(.3 * k), std::sin(.3 * k) * .6, std::sin(.3 * k) * .8};
            pc::Packet e{"e", "electron", 7, 1, w::me, electron_u(T)}, n{"n", "H2O", 7, 3, mn, {speed * ax[0], speed * ax[1], speed * ax[2]}};
            const auto cand = pc::candidates({e}, {n}, {{7, V}}, {composite}, 1, 1, 1);
            if (cand.size() != 1) { lab = false;continue; }
            const auto cat = scope.category(cand[0].energy_eV);
            if (cat == cv::Category::song_table) { lab = false;continue; }
            const auto band = cat == cv::Category::outer3_beb ? cv::Band::outer3_beb : cv::Band::provisional_rbeb;
            const auto s = cv::high_energy_parent_sample(he, band, composite, mn, cand[0], e, n, uni(), uni(), uni() * .999999L);
            long double before = w::kinetic(w::me, e.proper_velocity_m_s) + w::kinetic(mn, n.proper_velocity_m_s), after = s.products.particles[2].reference_J;
            for (const auto& q : s.products.particles) after += w::kinetic(q.mass_kg, q.proper_velocity_m_s);
            lab = lab && std::abs(after - before) <= 256 * std::numeric_limits<double>::epsilon() * (after + before) &&
                  s.products.particles[0].weight == 1 && s.products.particles[2].state == "H2O+";
            ++tried;
          }
      check(lab && tried == 600 && ignored == 0, "bridge: laboratory energy/momentum receipts close for moving targets (0, 600, 2e5 m/s)");
      // The Song band below 1 keV keeps the sealed equal-pair kernel.
      bool song = true;
      for (double T : {14., 100., 999.9}) {
        pc::Packet e{"e", "electron", 7, 1, w::me, electron_u(T)}, n{"n", "H2O", 7, 3, mn, {}};
        const auto cand = pc::candidates({e}, {n}, {{7, V}}, {bridge.channel()}, 1, 1, 1);
        if (cand.size() != 1) { song = false;continue; }
        const auto pr = bridge.sample(cand[0], e, n, {0, 0, 1});
        const long double k1 = kinetic_eV(pr.particles[0].proper_velocity_m_s), k2 = kinetic_eV(pr.particles[1].proper_velocity_m_s);
        const long double ki = w::kinetic(pr.particles[2].mass_kg, pr.particles[2].proper_velocity_m_s) / pc::qe;
        song = song && std::abs(k1 + k2 + ki - (cand[0].energy_eV - w::threshold_eV)) < 1e-9L * T && std::abs(k1 - k2) < 1e-3L * T;
      }
      check(song, "Song band keeps parent_ground_equal_pair_CM_v1 (sealed bridge) up to and including 1000 eV");
      receipts.push_back({{"case", "kernel conservation"}, {"rows", rows}});
      report["kernel_conservation"] = rows;
    }

    // ======================= C. classification and seams =======================
    const cv::CoverageScope scope(data, b.dataset_sha256, bridge, bands);
    const cv::CoverageScope open(data, b.dataset_sha256, bridge, bands, record_only<cv::CoveragePolicy>());
    const cv::CoverageScope v1eq(data, b.dataset_sha256, bridge);
    const auto no_prov = make_bands({false, 1e7});
    const cv::CoverageScope noprov(data, b.dataset_sha256, bridge, no_prov);
    const auto arena_v2 = scope.arena_binding(b.kernel_binding_file_sha256);
    const auto arena_v1 = v1eq.arena_binding(b.kernel_binding_file_sha256);
    {
      check(scope.identity() != open.identity() && scope.identity() != v1eq.identity() && scope.identity() != noprov.identity() &&
                make_bands({true, 2e7})->arena_binding().dataset_sha256 != arena_v2.dataset_sha256,
            "scope identity binds policy, high-energy bands, provisional admission and support");
      check(arena_v1.dataset_sha256 == b.dataset_sha256 && arena_v1.kernel_manifest_sha256 == b.kernel_sha256 &&
                arena_v1.kernel_binding_file_sha256 == b.kernel_binding_file_sha256 && arena_v1.model == w::parent_model &&
                arena_v2.dataset_sha256 != b.dataset_sha256 && arena_v2.model == cv::banded_model && pc::digest(arena_v2.dataset_sha256) &&
                pc::digest(arena_v2.kernel_manifest_sha256) && pc::digest(arena_v2.kernel_binding_file_sha256),
            "registered channel binding: v1 without bands, composite v2 identity with bands");
      const double s_inel = 2.800285205390781e-21;
      struct Row { double E;cv::Kind kind;cv::Sigma sigma;double sigma_m2;double upper; };
      const double up1 = std::nextafter(1000., 2000.), up10 = std::nextafter(10000., 2e4), up7 = std::nextafter(1e7, 2e7);
      const std::vector<Row> rows{
        {0, cv::Kind::known_zero, cv::Sigma::none, 0, 0}, {12.621, cv::Kind::known_zero, cv::Sigma::none, 0, 0},
        {std::nextafter(12.621, 13.), cv::Kind::onset_gap, cv::Sigma::bounded, 0, 1.6660000000000001 * s_inel},
        {13., cv::Kind::onset_gap, cv::Sigma::bounded, 0, 1.6660000000000001 * s_inel},
        {13.5, cv::Kind::supported, cv::Sigma::tabulated, data->observable_value_SI(w::parent_channel, 13.5), 0},
        {1000., cv::Kind::supported, cv::Sigma::tabulated, data->observable_value_SI(w::parent_channel, 1000.), 0},
        {up1, cv::Kind::beb_supported, cv::Sigma::dataset_formula, he.outer3_parent_m2(up1), 0},
        {8000., cv::Kind::beb_supported, cv::Sigma::dataset_formula, he.outer3_parent_m2(8000.), 0},
        {10000., cv::Kind::beb_supported, cv::Sigma::dataset_formula, he.outer3_parent_m2(10000.), 0},
        {up10, cv::Kind::provisional_supported, cv::Sigma::provisional_model, static_cast<double>(he.rbeb_parent_any(up10)), 0},
        {1e6, cv::Kind::provisional_supported, cv::Sigma::provisional_model, static_cast<double>(he.rbeb_parent_any(1e6)), 0},
        {1e7, cv::Kind::provisional_supported, cv::Sigma::provisional_model, static_cast<double>(he.rbeb_parent_any(1e7)), 0},
        {up7, cv::Kind::above_table, cv::Sigma::unbounded, 0, 0}, {1e9, cv::Kind::above_table, cv::Sigma::unbounded, 0, 0}};
      std::uint64_t counted = 0;const auto ch = scope.channel(bridge, &counted, b.kernel_binding_file_sha256);
      std::uint64_t expected_count = 0;
      for (const auto& r : rows) {
        const auto c = scope.classify(r.E);
        const bool excl = cv::unresolved(r.kind);expected_count += excl;
        const double s = ch.sigma_m2(r.E);
        check(c.kind == r.kind && c.sigma == r.sigma && c.sigma_m2 == r.sigma_m2 && c.upper_m2 == r.upper && scope.excluded(r.E) == excl &&
                  s == (excl ? 0. : r.sigma_m2), std::string("v2 coverage class and composite rate at ") + hex(r.E) + " eV");
      }
      check(counted == expected_count, "composite channel counts exactly the excluded energies");
      check(v1eq.classify(up1).kind == cv::Kind::above_table && v1eq.classify(8000.).sigma == cv::Sigma::bounded &&
                v1eq.classify(20000.).sigma == cv::Sigma::unbounded && v1eq.excluded(8000.) && !v1eq.excluded(1000.),
            "scope without bands keeps the v1 map above 1 keV");
      check(noprov.classify(8000.).kind == cv::Kind::beb_supported && noprov.classify(up10).kind == cv::Kind::above_table &&
                noprov.classify(up10).sigma == cv::Sigma::unbounded, "provisional band not admitted: above 10 keV unresolved without bound");
      const auto& seams = scope.description().at("seams");
      const double song1000 = data->observable_value_SI(w::parent_channel, 1000.);
      const double jump1 = (he.outer3_parent_m2(up1) - song1000) / song1000;
      const double jump10 = (static_cast<double>(he.rbeb_parent_any(up10)) - he.outer3_parent_m2(10000.)) / he.outer3_parent_m2(10000.);
      check(seams.at("song_to_outer3_BEB").at("relative_jump") == jump1 && jump1 > .2675 && jump1 < .2677 &&
                seams.at("outer3_BEB_to_provisional_RBEB").at("relative_jump") == jump10 && jump10 > .0257 && jump10 < .0259,
            "seams recorded in the scope identity: 1 keV Song->outer3 +26.76 %, 10 keV outer3->RBEB +2.58 %");
      reject([&] { (void)cv::HighEnergyBands::read(outer3_path, nw::sha("x"), *data, b.dataset_sha256, kernel_binding_path, b.kernel_sha256, b.kernel_binding_file_sha256); },
             "outer3 dataset bound by digest");
      reject([&] { (void)cv::HighEnergyBands(cv::read_bytes(outer3_path), outer3_sha, *data, b.dataset_sha256,
                                             "{\"schema\":\"parent-coverage-v2-kernel-binding\"}", b.kernel_sha256, b.kernel_binding_file_sha256); },
             "malformed kernel binding refused");
      reject([&] { (void)cv::HighEnergyBands::read(outer3_path, outer3_sha, *data, b.dataset_sha256, kernel_binding_path, nw::sha("other kernel"), b.kernel_binding_file_sha256); },
             "kernel binding bound to the Song-band kernel identity");
      reject([&] { (void)make_bands({true, 5000.}); }, "provisional support must lie above 10 keV");
      // Bit-level outer3 grid for the dataset-reference comparison.
      J grid_rows = J::array();
      for (double E : [] {
             std::vector<double> e{1000., std::nextafter(1000., 2000.), 1000.5, 1500., 2000., 3000., 5000., 8000., 8084.470604111418,
                                   9999.999, std::nextafter(10000., 0.), 10000.};
             for (int i = 0; i <= 200; ++i) e.push_back(1000. * std::pow(10., i / 200.));
             return e;
           }()) {
        if (E > 10000.) E = 10000.;
        grid_rows.push_back({{"energy", hex(E)}, {"parent", hex(he.outer3_parent_m2(E))}, {"valence", hex(he.valence_beb_m2(E))}});
      }
      std::ofstream(dir + "/outer3-grid.json") << grid_rows.dump(1) << '\n';
      report["seams"] = seams;
    }

    // ======================= D. subflow =======================
    // D1. Single high-energy electron, record-only policy, long exposure: the
    // whole cascade (BEB, RBEB and Song events) runs in one SSA.
    for (const double E : {1000.5, 8000., 20000., 1e5, 1e6}) {
      char label[64];std::snprintf(label, sizeof label, "single %g eV electron", E);const std::string L = label;
      Fixture c(arena_v2, {{1, E, here}}, 16);
      cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, open, 123, 456);
      auto trial = flow.prepare(interval(*c.arena, flow, 1e8));
      const auto& ev = trial.audit().at("events");const auto& cov = trial.audit().at("coverage");
      std::array<std::uint64_t, 3> by_cat{};bool tags = true, extras = true;
      for (const auto& x : ev) {
        const double energy = x.at("energy_eV");const auto cat = scope.category(energy);
        ++by_cat[static_cast<std::size_t>(cat)];
        tags = tags && x.at("event").get<std::string>().rfind(cv::event_prefix(cat), 0) == 0;
        if (cat == cv::Category::song_table) extras = extras && !x.contains("category") && x.size() == 10;
        else extras = extras && x.at("category") == cv::name(cat) && x.at("model_derived") == (cat == cv::Category::provisional_rbeb) &&
                      x.at("ejected_energy_target_rest_eV").get<double>() <= (energy - 12.61) / 2 && x.at("kernel") == cv::sdcs_kernel;
      }
      bool ledger = true;
      for (std::size_t k = 0; k < 3; ++k)
        ledger = ledger && category(cov, static_cast<cv::Category>(k)).at("executed_events") == by_cat[k] &&
                 value(category(cov, static_cast<cv::Category>(k)).at("executed_event_weight_sum")) == static_cast<double>(by_cat[k]);
      const auto first = scope.category(ev.at(0).at("energy_eV"));
      check(ev.size() >= 2 && first == scope.category(E) && tags && extras, L + ": events tagged by category; first event in the incident band");
      check(ledger && cov.at("exposure_limit").at("refused") == false, L + ": exact per-category executed counts and weights");
      check(flow.try_publish(trial) && weight(*c.arena, "H2O+") == ev.size() && weight(*c.arena, "electron") == 1 + ev.size() &&
                weight(*c.arena, "H2O") == 16 - ev.size(), L + ": native inventory: one H2O+ and one new electron per event");
      receipts.push_back({{"case", L}, {"events", events_without_native(trial.audit())}, {"coverage", cov}});
    }
    // D2. Mixed population across both seams, zero-event exposure: exact
    // per-category applied and sigma-free exposures against an independent
    // pair computation.
    {
      const std::vector<double> energies{13., 100., 999.9, 1000.1, 9999.9, 10000.1, 2e5};
      std::vector<Electron> lead;for (double E : energies) lead.push_back({1, E, here});
      Fixture c(arena_v2, lead, 16);
      cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, scope, 1, 2);
      const double dt = 1e-6;
      auto trial = flow.prepare(interval(*c.arena, flow, dt));
      const auto& cov = trial.audit().at("coverage");
      std::array<cv::ExactSum, 3> phys, ev_sig, free_p, free_e;std::array<std::uint64_t, 3> n{};std::array<std::set<std::size_t>, 3> bins;
      std::array<double, 3> lo{inf, inf, inf}, hi{0, 0, 0};
      for (double E : energies) {
        pc::Packet e{"e", "electron", 0, 1, w::me, electron_u(E)}, nn{"n", "H2O", 0, 16, mn, {}};
        const auto kin = pc::relative(e, nn);const double x = kin.target_frame_energy_eV;
        if (scope.excluded(x)) continue;
        const auto cat = scope.category(x);const auto k = static_cast<std::size_t>(cat);
        const double sigma = cat == cv::Category::song_table ? data->observable_value_SI(w::parent_channel, x)
                           : cat == cv::Category::outer3_beb ? bb::beb_weighted_sum_m2(orb, he.weights(), K, x)
                                                             : static_cast<double>(he.rbeb_parent_any(x));
        const double rate = static_cast<double>(16.L * sigma * kin.moller_flux_m_s / V);
        ++n[k];ev_sig[k].add(static_cast<double>(static_cast<long double>(rate) * dt));
        phys[k].add(static_cast<double>(static_cast<long double>(rate) * 1 * dt));
        free_p[k].add(static_cast<double>(1.L * 16 * kin.moller_flux_m_s / V * dt));free_e[k].add(static_cast<double>(16.L * kin.moller_flux_m_s / V * dt));
        bins[k].insert(scope.bin_of(x));lo[k] = std::min(lo[k], x);hi[k] = std::max(hi[k], x);
      }
      bool exact = true;
      for (std::size_t k = 0; k < 3; ++k) {
        const auto& r = category(cov, static_cast<cv::Category>(k));
        exact = exact && r.at("applied").at("candidate_evaluations") == n[k] && r.at("applied").at("expected_physical_reactions") == cv::exact_json(phys[k]) &&
                r.at("applied").at("expected_events") == cv::exact_json(ev_sig[k]) && r.at("sigma_free_physical_exposure_inv_m2") == cv::exact_json(free_p[k]) &&
                r.at("sigma_free_event_exposure_inv_m2") == cv::exact_json(free_e[k]) && r.at("bins").size() == bins[k].size() &&
                r.at("energy_min_eV") == lo[k] && r.at("energy_max_eV") == hi[k] && r.at("executed_events") == 0;
      }
      check(trial.audit().at("events").empty() && n[0] == 2 && n[1] == 2 && n[2] == 2, "mixed seams: zero-event exposure with two pairs per supported band");
      check(exact, "mixed seams: per-category applied and sigma-free exposures equal the independent pair integrals exactly");
      cv::Applied total;for (std::size_t k = 0; k < 3; ++k) total.merge(cv::Applied::from(category(cov, static_cast<cv::Category>(k)).at("applied")));
      check(total.json() == cov.at("applied") && cov.at("unresolved").at("pair_evaluations") == 1 && cov.at("exposure_limit").at("refused") == false,
            "mixed seams: applied total is the sum of the categories; the 13 eV pair stays unresolved");
      check(category(cov, cv::Category::provisional_rbeb).at("model_derived") == true && category(cov, cv::Category::outer3_beb).at("basis") == "dataset" &&
                category(cov, cv::Category::song_table).at("basis") == "dataset", "ledger categories: dataset Song, dataset outer3 BEB, model-derived RBEB");
      check(flow.try_publish(trial), "mixed seams: exposure publishes");
      receipts.push_back({{"case", "mixed populations across both seams"}, {"coverage", cov}});
      // Same population over a long exposure with record-only policy: events
      // in all three categories, each dispatched to its kernel.
      Fixture d(arena_v2, lead, 1e20);
      auto big = limits;big.reactions = 1000;
      cv::ParentSubflow f2(*d.arena, grid.counts, b, big, bridge, open, 5, 6);
      auto t2 = f2.prepare(interval(*d.arena, f2, 3e-16));
      std::array<std::uint64_t, 3> seen{};
      for (const auto& x : t2.audit().at("events")) ++seen[static_cast<std::size_t>(scope.category(x.at("energy_eV")))];
      bool consistent = true;
      for (std::size_t k = 0; k < 3; ++k) consistent = consistent && category(t2.audit().at("coverage"), static_cast<cv::Category>(k)).at("executed_events") == seen[k];
      check(seen[0] > 0 && seen[1] > 0 && seen[2] > 0 && consistent && f2.try_publish(t2), "mixed seams: events of all three categories in one SSA, ledger consistent");
      receipts.push_back({{"case", "mixed populations, all categories executed"}, {"events", events_without_native(t2.audit())}, {"coverage", t2.audit().at("coverage")}});
    }
    // D3. Parity with v1 below 1 keV.
    {
      c1::CoverageScope s1(data, b.dataset_sha256, bridge);
      for (const auto& [we, nwgt, duration] : std::vector<std::tuple<double, double, double>>{{1, 16, 1e8}, {.25, 1e20, 1e-5}}) {
        const std::string L = "v1 parity w_e=" + std::to_string(we);
        Fixture a(arena_v1, {{we, 100, here}}, nwgt), c(arena_v2, {{we, 100, here}}, nwgt), e(arena_v1, {{we, 100, here}}, nwgt);
        c1::ParentSubflow old(*a.arena, grid.counts, b, limits, bridge, s1, 123, 456);
        cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, scope, 123, 456);
        cv::ParentSubflow eq(*e.arena, grid.counts, b, limits, bridge, v1eq, 123, 456);
        auto ot = old.prepare(interval(*a.arena, old, duration));auto nt = flow.prepare(interval(*c.arena, flow, duration));
        auto qt = eq.prepare(interval(*e.arena, eq, duration));
        check(events_without_native(ot.audit()) == events_without_native(nt.audit()) && events_without_native(ot.audit()) == events_without_native(qt.audit()) &&
                  nt.audit().at("events").size() >= 2, L + ": event physics identical to v1 (bands present and absent)");
        check(old.try_publish(ot) && flow.try_publish(nt) && eq.try_publish(qt) && rng_fields(old.checkpoint()) == rng_fields(flow.checkpoint()) &&
                  rng_fields(old.checkpoint()) == rng_fields(eq.checkpoint()) && census(*a.arena) == census(*c.arena) && census(*a.arena) == census(*e.arena),
              L + ": scheduler/product RNG streams and native carriers identical to v1");
        const auto& cov = nt.audit().at("coverage");
        check(category(cov, cv::Category::outer3_beb).at("applied").at("candidate_evaluations") == 0 &&
                  category(cov, cv::Category::provisional_rbeb).at("applied").at("candidate_evaluations") == 0 &&
                  cov.at("applied") == ot.audit().at("coverage").at("applied") &&
                  category(cov, cv::Category::song_table).at("executed_events") == nt.audit().at("events").size(),
              L + ": applied measure equals v1; all executed events in the Song category");
      }
      // A high-energy electron in a disjoint cell changes nothing below 1 keV.
      Fixture ref(arena_v1, {{1, 100, here}}, 16, {{.5, 8000., elsewhere}});
      c1::ParentSubflow old(*ref.arena, grid.counts, b, limits, bridge, s1, 123, 456);
      auto ot = old.prepare(interval(*ref.arena, old, 1e8));check(old.try_publish(ot), "v1 reference with a disjoint 8 keV electron");
      Fixture mix(arena_v2, {{1, 100, here}}, 16, {{.5, 8000., elsewhere}});
      cv::ParentSubflow flow(*mix.arena, grid.counts, b, limits, bridge, scope, 123, 456);
      auto nt = flow.prepare(interval(*mix.arena, flow, 1e8));
      check(events_without_native(ot.audit()) == events_without_native(nt.audit()) && flow.try_publish(nt) &&
                rng_fields(old.checkpoint()) == rng_fields(flow.checkpoint()) && census(*mix.arena) == census(*ref.arena),
            "v1 parity with a disjoint-cell 8 keV electron: events, RNG and carriers identical");
    }
    // D4. Refusal limit now that the bands above 1 keV are covered.
    {
      struct Case { const char* label;double E;bool refused;const char* criterion;cv::Category cat; };
      for (const auto& r : {Case{"1000.5 eV (outer3)", 1000.5, false, "none", cv::Category::outer3_beb},
                            Case{"8000 eV (outer3)", 8000., false, "none", cv::Category::outer3_beb},
                            Case{"20000 eV (provisional)", 20000., false, "none", cv::Category::provisional_rbeb},
                            Case{"1 MeV (provisional)", 1e6, false, "none", cv::Category::provisional_rbeb},
                            Case{"20 MeV (above provisional support)", 2e7, true, "unbounded_unresolved_sigma", cv::Category::song_table},
                            Case{"13 eV (onset gap)", 13., true, "events_and_fraction_upper_bounds", cv::Category::song_table}}) {
        const std::string L = std::string("default limit, ") + r.label;
        Fixture c(arena_v2, {{1, r.E, here}}, 16);cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, scope, 1, 2);
        const auto cp = flow.checkpoint();
        if (r.refused) {
          const auto x = refusal([&] { flow.prepare(interval(*c.arena, flow, 1e3)); }, L + ": refused");
          check(x.criterion == r.criterion && flow.checkpoint() == cp, L + ": criterion recorded, owners unchanged");
        } else {
          auto t = flow.prepare(interval(*c.arena, flow, 1.));
          const auto& cov = t.audit().at("coverage");
          check(cov.at("exposure_limit").at("refused") == false && cov.at("unresolved").at("pair_evaluations") == 0 &&
                    category(cov, r.cat).at("applied").at("candidate_evaluations").get<std::uint64_t>() >= 1 && flow.try_publish(t),
                L + ": covered (v1 refused), applied in its category");
        }
      }
      // v1 behavior above 10 keV when the provisional band is not admitted.
      {
        const auto a = noprov.arena_binding(b.kernel_binding_file_sha256);
        Fixture c(a, {{1, 20000., here}}, 16);cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, noprov, 1, 2);
        const auto x = refusal([&] { flow.prepare(interval(*c.arena, flow, 1.)); }, "provisional band not admitted: 20 keV refused as unbounded");
        check(x.criterion == "unbounded_unresolved_sigma", "provisional band not admitted: criterion is the unbounded sigma");
      }
      // Covered high-energy hazard enters the limit denominator: a small gap
      // weight beside a large outer3/RBEB applied measure is admitted, a
      // comparable one refused.
      cv::CoveragePolicy tight;tight.max_exposure_unresolved_events = 1e-4;
      const cv::CoverageScope tscope(data, b.dataset_sha256, bridge, bands, tight);
      for (const auto& [gap_w, refused] : std::vector<std::pair<double, bool>>{{1e-3, false}, {1., true}}) {
        Fixture c(arena_v2, {{1, 8000., here}, {1, 1e5, here}}, 1e20, {{gap_w, 13., here}});
        cv::ParentSubflow flow(*c.arena, grid.counts, b, limits, bridge, tscope, 3, 4);
        if (refused) {
          const auto x = refusal([&] { flow.prepare(interval(*c.arena, flow, 1e-18)); }, "limit: comparable gap weight beside covered keV electrons refused");
          check(x.exposure_ledger.at("exposure_limit").at("fraction_exceeded") == true, "limit: refusal through the fraction criterion");
        } else {
          auto t = flow.prepare(interval(*c.arena, flow, 1e-18));const auto& d = t.audit().at("coverage").at("exposure_limit");
          check(d.at("refused") == false && d.at("absolute_exceeded") == true && d.at("fraction_exceeded") == false &&
                    d.at("applied_physical").get<double>() > 0, "limit: high-energy applied measure admits a small gap fraction");
        }
      }
    }
    // D5. Restart, rollback and tampering with categories.
    {
      Fixture c(arena_v2, {{1, 8000., here}, {1, 3e5, here}, {1, 3e5, here}}, 1e20);
      auto big = limits;big.reactions = 1000;
      cv::ParentSubflow flow(*c.arena, grid.counts, b, big, bridge, open, 77, 78);
      auto t0 = flow.prepare(interval(*c.arena, flow, 4e-16, 0));
      check(!t0.audit().at("events").empty() && flow.try_publish(t0), "restart: first exposure with high-energy events");
      c.arena->save(dir + "/native-k1");const auto k1 = flow.checkpoint();std::ofstream(dir + "/subflow-k1.json") << k1.dump(2);
      auto t1 = flow.prepare(interval(*c.arena, flow, 4e-16, 1));check(flow.try_publish(t1), "restart: second exposure");
      const auto k2 = flow.checkpoint();
      m::GapPic back{options()};
      auto a1 = m::PicArena::restore(back, dir + "/native-k1", states(), channels(arena_v2), "native-water-subflow-fixture/v1");
      cv::ParentSubflow r1(*a1, grid.counts, b, big, bridge, open, 9, 9);
      reject([&] { r1.restore(k2); }, "rollback: newer ledger refused against the older native history");
      r1.restore(k1);check(r1.checkpoint() == k1, "rollback: checkpoint with category ledger restored exactly");
      auto tampered = k1;
      tampered["coverage"]["cumulative"]["categories"][2]["executed_events"] = 0;
      reject([&] { r1.restore(tampered); }, "restart: tampered model-derived event count refused");
      tampered = k1;tampered["coverage"]["cumulative"]["categories"][2]["model_derived"] = false;
      reject([&] { r1.restore(tampered); }, "restart: relabeled model-derived category refused");
      tampered = k1;tampered["coverage"]["cumulative"]["categories"][1]["applied"]["candidate_evaluations"] = 0;
      reject([&] { r1.restore(tampered); }, "restart: tampered outer3 applied count refused");
      cv::ParentSubflow other(*a1, grid.counts, b, big, bridge, scope, 9, 9);
      reject([&] { other.restore(k1); }, "restart: ledger of another policy refused");
      auto replay = r1.prepare(interval(*a1, r1, 4e-16, 1));
      check(replay.audit().at("coverage") == t1.audit().at("coverage") && events_without_native(replay.audit()) == events_without_native(t1.audit()),
            "restart: next exposure, its high-energy events and category ledger byte-equivalent after restart");
      check(r1.try_publish(replay) && r1.checkpoint().at("coverage") == k2.at("coverage") && rng_fields(r1.checkpoint()) == rng_fields(k2),
            "restart: accepted ledger and RNG identical after replay");
      const auto cum = k2.at("coverage").at("cumulative");
      check(cum.at("categories").at(2).at("executed_events").get<std::uint64_t>() > 0 && cum.at("categories").at(1).at("executed_events").get<std::uint64_t>() > 0,
            "restart: cumulative ledger carries executed dataset-BEB and model-derived events");
      receipts.push_back({{"case", "restart with categories"}, {"cumulative", cum}});
    }
    // D6. Registration and configuration identity.
    {
      Fixture c(arena_v1, {{1, 100, here}}, 16);
      reject([&] { cv::ParentSubflow bad(*c.arena, grid.counts, b, limits, bridge, scope, 1, 2); }, "v2 bands refused against a v1 channel registration");
      Fixture d(arena_v2, {{1, 100, here}}, 16);
      reject([&] { cv::ParentSubflow bad(*d.arena, grid.counts, b, limits, bridge, v1eq, 1, 2); }, "v1-equivalent scope refused against the composite registration");
      auto wrong = b;wrong.kernel_binding_file_sha256 = nw::sha("wrong binding");
      reject([&] { cv::ParentSubflow bad(*d.arena, grid.counts, wrong, limits, bridge, scope, 1, 2); }, "Song-band binding file enforced with bands");
      const auto r_v2 = cv::ParentSubflow::initial_full_time_rng_sha256(options(), grid.counts, b, limits, scope, 5, 6);
      const auto r_eq = cv::ParentSubflow::initial_full_time_rng_sha256(options(), grid.counts, b, limits, v1eq, 5, 6);
      const auto r_np = cv::ParentSubflow::initial_full_time_rng_sha256(options(), grid.counts, b, limits, noprov, 5, 6);
      const auto r_v1 = c1::ParentSubflow::initial_full_time_rng_sha256(options(), grid.counts, b, limits, c1::CoverageScope(data, b.dataset_sha256, bridge), 5, 6);
      check(r_v2 != r_eq && r_v2 != r_np && r_eq != r_np && r_eq != r_v1 &&
                cv::ParentSubflow::full_time_configuration_sha256(options(), grid.counts, b, limits, scope) !=
                    cv::ParentSubflow::full_time_configuration_sha256(options(), grid.counts, b, limits, noprov),
            "configuration and genesis RNG identities bind the dataset/model bands");
    }
    // D7. Full-time first exposure with a high-energy electron.
    {
      auto make = [&](const cv::CoverageScope& s, double E) {
        m::FullTimeGenesisBinding g;g.time_s = 7.4978805549199997e-10;g.origin_provenance = "authored/coverage-fixture-no-simulation";
        g.field_sha256 = nw::sha("coverage-fixture/field");g.mesh_sha256 = nw::sha("coverage-fixture/mesh");
        g.material_sha256 = nw::sha("coverage-fixture/material");g.frame_sha256 = nw::sha("coverage-fixture/frame");
        g.common_clock_sha256 = nw::sha("coverage-fixture/clock");g.heat_sha256 = nw::sha("coverage-fixture/heat");
        g.source_sha256 = nw::sha("coverage-fixture/source");
        auto big = limits;big.reactions = 1000;
        g.rng_sha256 = cv::ParentSubflow::initial_full_time_rng_sha256(options(), grid.counts, b, big, s, 5, 6);
        return m::PicArena::create_full_time(options(), states(), channels(s.arena_binding(b.kernel_binding_file_sha256)),
          {{"e", "electron", "authored/fixture", {1, 0}, here, electron_u(E)}, {"n", "H2O", "authored/fixture", {3e19, 0}, water_at, {0, 0, 0}}}, g);
      };
      auto a = make(open, 5e4);
      auto big = limits;big.reactions = 1000;
      cv::ParentSubflow flow(*a, grid.counts, b, big, bridge, open, 5, 6, true);
      const auto cp = flow.checkpoint();const auto rng = flow.accepted_rng_sha256();
      {
        const auto bd = a->full_time_boundary();
        m::FullTimeStepBinding common{bd->accepted_step, bd->accepted_arena_revision, bd->time_s, bd->time_s + 4e-15, 4e-15,
                                      bd->accepted_arena_sha256, bd->field_sha256, bd->mesh_sha256, bd->material_sha256, nw::sha("coverage-common")};
        auto outer = a->begin_full_time_step(common, 100000, 4096);
        m::ExposureBinding e{common, m::FullTimeSlot::before_transport, common.begin_s, .5 * common.dt_s, 0, common.accepted_field_sha256};
        auto first_scope = outer.first_exposure(e);
        auto stage = flow.prepare_on(first_scope);
        const auto& cov = stage.audit().at("coverage");
        std::uint64_t model = 0;bool tagged = true;
        for (const auto& x : stage.audit().at("events")) {
          const bool md = x.at("event").get<std::string>().rfind("model-derived/", 0) == 0;model += md;
          tagged = tagged && md == (scope.category(x.at("energy_eV")) == cv::Category::provisional_rbeb);
        }
        check(model > 0 && tagged && cov.at("slot") == "before_transport" &&
                  category(cov, cv::Category::provisional_rbeb).at("executed_events") == model,
              "full-time first exposure: model-derived events staged with their category ledger");
        check(flow.staged_rng_sha256(stage) != rng && pc::digest(flow.staged_rng_sha256(stage)), "full-time first exposure: staged ledger bound into the staged RNG identity");
        receipts.push_back({{"case", "full-time first exposure, 50 keV"}, {"coverage", cov}});
      }
      check(flow.checkpoint() == cp && flow.accepted_rng_sha256() == rng, "full-time: abandoned private step leaves the accepted ledger unchanged");
    }

    // ======================= E. report tables =======================
    {
      J corr = J::array(), cmp = J::array(), meas = J::array(), gap = J::array();
      for (double E : {1000., 1500., 2000., 3000., 5000., 8000., 10000., 20000., 50000., 1e5, 2e5, 5e5, 1e6}) {
        const long double nr = he.beb_parent_any(E), r = he.rbeb_parent_any(E), nv = he.beb_valence_any(E), rv = he.rbeb_valence_any(E);
        corr.push_back({{"energy_eV", E}, {"nonrelativistic_BEB_parent_m2", static_cast<double>(nr)}, {"RBEB_parent_m2", static_cast<double>(r)},
                        {"relative_correction_parent", static_cast<double>(r / nr - 1)}, {"nonrelativistic_BEB_valence_m2", static_cast<double>(nv)},
                        {"RBEB_valence_m2", static_cast<double>(rv)}, {"relative_correction_valence", static_cast<double>(rv / nv - 1)},
                        {"parent_fraction_RBEB", static_cast<double>(r / rv)}, {"parent_fraction_BEB", static_cast<double>(nr / nv)}});
        if (E <= 10000) cmp.push_back({{"energy_eV", E}, {"outer3_BEB_Tan_parent_v1_m2", he.outer3_parent_m2(E)},
                                       {"RBEB_x_Tan_parent_m2", static_cast<double>(r)}, {"relative_difference", static_cast<double>(r / he.outer3_parent_m2(E) - 1)}});
      }
      for (double E : {1000., 1500., 2000., 3000., 4000., 5000.}) {
        const double m = data->observable_value_SI("H2O:munoz2007:measured_ionization", E);
        meas.push_back({{"energy_eV", E}, {"munoz2007_measured_gross_ionization_m2", m}, {"BEB_valence_total_m2", he.valence_beb_m2(E)},
                        {"RBEB_valence_total_m2", static_cast<double>(he.rbeb_valence_any(E))},
                        {"RBEB_over_measured", static_cast<double>(he.rbeb_valence_any(E) / m)}, {"BEB_over_measured", he.valence_beb_m2(E) / m}});
      }
      const double song_total = data->observable_value_SI("H2O:ion_yield:total", 1000.);
      // Probability that an ejected electron starts in the onset gap
      // (12.621, 13.5) eV or at/below 12.621 eV, per incident energy.
      for (double T : {1000.5, 5000., 10000., 1e5, 1e6}) {
        const auto terms = T <= 10000 ? he.outer3_terms(T) : he.provisional_terms(T);
        long double tot = 0, pg = 0, pz = 0, mean = 0;for (auto x : terms) tot += x;
        for (std::size_t i = 1; i < 4; ++i) {
          const auto s = T <= 10000 ? bb::beb_shape(orb[i], K, T) : bb::rbeb_shape(orb[i], K, T);
          const long double Gm = s.G(s.w_max()), p = terms[i] / tot;
          pz += p * s.G(w::threshold_eV / s.B) / Gm;pg += p * (s.G(13.5L / s.B) - s.G(w::threshold_eV / s.B)) / Gm;
          // Mean ejected energy by quadrature of W g(w).
          const long double top = std::log1p(s.w_max());long double num = 0;const int P = 2000;
          for (int k = 0; k < P; ++k) { const long double x = top * (k + .5L) / P, wv = std::expm1(x);num += s.g(wv) * (1 + wv) * wv * s.B * top / P; }
          mean += p * num / Gm;
        }
        gap.push_back({{"energy_eV", T}, {"P_ejected_at_or_below_12.621_eV", static_cast<double>(pz)}, {"P_ejected_in_onset_gap", static_cast<double>(pg)},
                       {"mean_ejected_eV", static_cast<double>(mean)}});
      }
      report["relativistic_correction"] = corr;report["outer3_vs_RBEB_1_10keV"] = cmp;
      report["total_ionization_vs_measured"] = {{"munoz2007_measured", meas}, {"song2021_total_at_1000eV_m2", song_total},
                                                {"RBEB_valence_at_1000eV_m2", static_cast<double>(he.rbeb_valence_any(1000.))},
                                                {"note", "Munoz 2007 Table I is a gross (ion-current, charge-weighted) observable including O 1s and "
                                                         "multiple ionization; the four-orbital BEB/RBEB total is a counting valence cross section."}};
      report["ejected_electron_thresholds"] = gap;
      // Values compared with the independent 50-digit oracle by the runner.
      J oracle = J::array();
      for (double T : {13., 20., 50., 100., 1000., 1500., 2000., 3000., 5000., 8000., 10000., 20000., 50000., 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 1e8, 1e9}) {
        J per = J::array();
        for (const auto& o : orb)
          per.push_back({{"label", o.label}, {"beb_dataset_arithmetic_m2", bb::beb_orbital_m2(o, K, T)},
                         {"beb_extended_m2", static_cast<double>(bb::beb_orbital_ext_m2(o, K, T))},
                         {"rbeb_m2", static_cast<double>(bb::rbeb_orbital_m2(o, K, T))}});
        oracle.push_back({{"energy_eV", T}, {"orbitals", per}, {"beb_parent_m2", static_cast<double>(he.beb_parent_any(T))},
                          {"rbeb_parent_m2", static_cast<double>(he.rbeb_parent_any(T))}, {"beb_valence_m2", static_cast<double>(he.beb_valence_any(T))},
                          {"rbeb_valence_m2", static_cast<double>(he.rbeb_valence_any(T))}});
      }
      report["oracle_values"] = oracle;
      report["song_parent_m2_at_1000eV"] = data->observable_value_SI(w::parent_channel, 1000.);
    }

    J result{{"status", "passed"}, {"checks", checks}, {"check_count", checks.size()}, {"coverage_scope", scope.description()},
             {"coverage_scope_sha256", scope.identity()}, {"record_only_scope_sha256", open.identity()},
             {"v1_equivalent_scope_sha256", v1eq.identity()}, {"arena_channel_binding", arena_v2.json()},
             {"report", report}, {"receipts", receipts}, {"conditional_parent_only", true},
             {"physical_trajectory_advanced", false}, {"full_coupled_run", false}};
    std::ofstream(dir + "/result.json") << result.dump(2) << '\n';
    std::cout << "Passed " << checks.size() << " checks\n";
  } catch (const std::exception& e) { std::cerr << e.what() << '\n';return 1; }
}
