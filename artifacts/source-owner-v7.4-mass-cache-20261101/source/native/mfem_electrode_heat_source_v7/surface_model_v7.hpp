#pragma once
// Heat/source owner v7 -- electrode surface physics (host side).
//
//   Registry        every parameter value with its tag (measured | recommended | model_derived |
//                   scenario | benchmark), band [low, high], unit and citation (DATA-v7.txt).
//   SurfaceModel    work function of oxidized SS304 with adsorbed water:
//                     phi(theta, T) = phi_ox + dphi_1ML g(theta) + dphi_dT (T - T_ref),
//                     g(theta) = theta (1 + c) / (1 + c theta) for theta <= 1, 1 above,
//                   theta = exact adsorbed stock / (patch area x monolayer density); phi is put on
//                   the grid 2^-12 eV (the fast spectrum caches one configuration per value).
//   Enhancement     field-enhancement classes (beta_k, area fraction a_k, local extent r_k): each
//                   class is one exposure of the sealed planar-barrier kernel at beta_k F with area
//                   a_k A; phase map per patch and class.
//   WallModel       parameters of the wall kernel device/sw7_wall.h (SW7Params, species table).
//   Ledger          exact integer ledgers of wall events (energy quanta 2^-30 eV, charges, counts).
//   IntervalBound   the source-owner step bound (API-v7.txt section 3).
//   Photo           optional, bounded photo-emission hook.
// Host only (the kernel header is shared with NVRTC).
#include "device/sw7_wall.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace surface_v7 {

// ------------------------------------------------------------------ registry
enum class Tag { measured, recommended, model_derived, scenario, benchmark };
inline const char* tag_name(Tag t) {
  switch (t) {
    case Tag::measured: return "measured";
    case Tag::recommended: return "recommended";
    case Tag::model_derived: return "model_derived";
    case Tag::scenario: return "scenario";
    case Tag::benchmark: return "benchmark";
  }
  return "?";
}
struct Param {
  std::string name;
  double value = 0, low = 0, high = 0;
  std::string unit;
  Tag tag = Tag::scenario;
  std::string source, location, note;
};
struct Registry {
  std::vector<Param> rows;
  double use(const std::string& name, double value, double low, double high, const std::string& unit, Tag tag, const std::string& source,
             const std::string& location, const std::string& note = "") {
    if (!(std::isfinite(value) && low <= value && value <= high)) throw std::invalid_argument("v7 parameter outside its own band: " + name);
    rows.push_back({name, value, low, high, unit, tag, source, location, note});
    return value;
  }
  double fixed(const std::string& name, double value, const std::string& unit, Tag tag, const std::string& source, const std::string& location,
               const std::string& note = "") {
    return use(name, value, value, value, unit, tag, source, location, note);
  }
};

// Sources (full citations in DATA-v7.txt).
inline constexpr const char* FP2002 = "furman_pivi_2002";          // PRST-AB 5, 124404, doi:10.1103/PhysRevSTAB.5.124404
inline constexpr const char* JOY_DB = "joy_database_2004";         // D. C. Joy, electron-solid interaction database (Fe eta: Hunger & Kuchler 1979; Palluel 1947)
inline constexpr const char* HUNGER1979 = "hunger_kuchler_1979";   // phys. stat. sol. (a) 56, K45 (1979)
inline constexpr const char* GRARD1973 = "grard_1973";             // J. Geophys. Res. 78, 2885, doi:10.1029/JA078i016p02885 (Feuerbacher & Fitton 1972 data)
inline constexpr const char* BAGLIN2000 = "baglin_epac2000";       // Baglin et al., EPAC 2000, p. 217 (as-received stainless SEY)
inline constexpr const char* MG1956 = "murphy_good_1956";          // Phys. Rev. 102, 1464, doi:10.1103/PhysRev.102.1464
inline constexpr const char* FORBES2006 = "forbes_2006";           // Appl. Phys. Lett. 89, 113122, doi:10.1063/1.2354582
inline constexpr const char* V61 = "source_owner_v6.1";            // runtimes/pic-electrode-heat-source-v6.1-20261031
inline constexpr const char* CAL_V7 = "coordinator_v7_calibration"; // runtimes/pic-water-common-commit-v7-20261031 README 6.2
inline constexpr const char* THIS = "this_component";              // derived here; derivation in API-v7.txt

// ------------------------------------------------------------------ work function
struct SurfaceModel {
  bool benchmark = true;          // phi = fixed_phi (v6.1: 4.3 eV clean-surface benchmark)
  double fixed_phi_eV = 4.3;
  double phi_ox_eV = 4.3;         // oxidized SS304 without adsorbed water
  double dphi_1ML_eV = 0;         // work-function change at one adsorbed water monolayer (negative lowers)
  double depolarization_c = 0;    // c of g(theta)
  double dphi_dT_eV_K = 0;
  double T_ref_K = 300;
  double grid_eV = 0x1p-12;       // phi grid
  double band_low_eV = 4.3, band_high_eV = 4.3;   // uncertainty band of phi at theta = 0 (oxide)
  double phi_liquid_eV = 0;       // molten (oxide-free) steel surface; 0: not configured (liquid class refused)
  double band_liquid_low_eV = 0, band_liquid_high_eV = 0;
  bool operator==(const SurfaceModel&) const = default;
};
inline double coverage_g(double theta, double c) {
  if (!(theta > 0)) return 0;
  if (theta >= 1) return 1;
  return theta * (1 + c) / (1 + c * theta);
}
inline double work_function_eV(const SurfaceModel& m, double theta, double T) {
  if (m.benchmark) return m.fixed_phi_eV;
  if (!(std::isfinite(theta) && theta >= 0 && std::isfinite(T) && T > 0)) throw std::invalid_argument("v7 surface state: finite coverage and temperature required");
  const double phi = m.phi_ox_eV + m.dphi_1ML_eV * coverage_g(theta, m.depolarization_c) + m.dphi_dT_eV_K * (T - m.T_ref_K);
  const double q = std::nearbyint(phi / m.grid_eV) * m.grid_eV;   // exact grid value
  if (!(q >= 2 && q <= 6)) throw std::invalid_argument("v7 surface state: work function outside the kernel domain 2..6 eV");
  return q;
}

// ------------------------------------------------------------------ field enhancement
struct EnhancementClass {
  double beta = 1, area_fraction = 1;
  double local_extent_m = INFINITY;   // emitter scale bounding the local field region (tip radius); inf for the flat class
  std::string label = "flat";
  bool operator==(const EnhancementClass&) const = default;
};
struct Enhancement {
  std::vector<EnhancementClass> classes{EnhancementClass{}};
  bool operator==(const Enhancement&) const = default;
};
inline void validate(const Enhancement& e) {
  if (e.classes.empty() || e.classes.size() > 32) throw std::invalid_argument("v7 enhancement: 1..32 classes required");
  double sum = 0;
  for (const auto& c : e.classes) {
    if (!(std::isfinite(c.beta) && c.beta >= 1 && c.beta <= 1e4 && std::isfinite(c.area_fraction) && c.area_fraction > 0 && c.area_fraction <= 1 &&
          c.local_extent_m > 0))
      throw std::invalid_argument("v7 enhancement class: beta in [1, 1e4], area fraction in (0, 1], positive extent");
    sum += c.area_fraction;
  }
  if (!(std::abs(sum - 1) <= 1e-12)) throw std::invalid_argument("v7 enhancement: area fractions must sum to 1");
}

// ------------------------------------------------------------------ wall model
struct WallModel {
  bool enabled = false;               // false: v6.1 behaviour (every first hit absorbed by the Hall sink)
  SW7Params p{};
  std::array<SW7Species, SW7_NSPECIES> species{};
  bool hall_applies_incident_kinetic = true;   // Hall's sink wall heat (incident K) stays; v7 loads carry the rest
  bool operator==(const WallModel& o) const;
};
// Bytewise identity of the parameter blocks (used for configuration digests and equality).
template <class T> std::string bytes_of(const T& x) { return std::string(reinterpret_cast<const char*>(&x), sizeof(T)); }
inline bool WallModel::operator==(const WallModel& o) const {
  return enabled == o.enabled && hall_applies_incident_kinetic == o.hall_applies_incident_kinetic && bytes_of(p) == bytes_of(o.p) && bytes_of(species) == bytes_of(o.species);
}

// ------------------------------------------------------------------ photo hook
struct PhotoModel {
  bool enabled = false;
  std::vector<double> E_eV, yield;    // quantum yield table (electrons per photon), log-log interpolation
  double scale = 1;
  bool operator==(const PhotoModel&) const = default;
};
struct PhotoExposure { double expected_electrons = 0, yield = 0, max_kinetic_eV = 0; bool outside_table = false; };
inline PhotoExposure photo_exposure(const PhotoModel& m, double photon_flux_m2_s, double photon_eV, double phi_eV, double area_m2, double dt_s) {
  PhotoExposure x;
  if (!m.enabled || !(photon_flux_m2_s > 0) || !(photon_eV > phi_eV)) return x;
  bool out = false;
  x.yield = sw7::interp_loglog(m.E_eV.data(), m.yield.data(), int(m.E_eV.size()), photon_eV, &out) * m.scale;
  x.outside_table = out;
  x.expected_electrons = x.yield * photon_flux_m2_s * area_m2 * dt_s;
  x.max_kinetic_eV = photon_eV - phi_eV;   // Einstein bound
  return x;
}

// ------------------------------------------------------------------ surface options per patch
struct SurfaceOptions {
  SurfaceModel work_function;
  Enhancement enhancement;
  WallModel wall;
  PhotoModel photo;
  bool benchmark() const { return work_function.benchmark && enhancement == Enhancement{} && !wall.enabled && !photo.enabled; }
};

// ------------------------------------------------------------------ species table
inline std::array<SW7Species, SW7_NSPECIES> species_table(Registry* R) {
  // Recombination energies: ionization energies of the neutral products (NIST Chemistry WebBook;
  // H3O+: recombination to H2O + H, see DATA-v7.txt).
  auto ie = [&](const char* n, double v, double lo, double hi, const char* src, const char* loc) { return R ? R->use(std::string("species.") + n, v, lo, hi, "eV", Tag::measured, src, loc) : v; };
  const double amu = 1.66053906660e-27;
  std::array<SW7Species, SW7_NSPECIES> t{};
  // Reflection: Thomas, Janev & Smith (IAEA INDC(NDS)-249, 1991; NIM B 69 (1992) 427) reduced energies on Fe,
  // eps = eps_L E_atom[keV]: H->Fe 0.3934 keV^-1 (Table 1), O->Fe 0.0340 keV^-1 (Eq. (2), model_derived);
  // a molecular projectile reflects as its heaviest atom with that atom's energy share (model_derived).
  const double epsH = R ? R->use("reflection.epsL.H_on_Fe", 0.3934, 0.3575, 0.3934, "1/keV", Tag::recommended, "thomas_janev_smith_1991", "INDC(NDS)-249 Table 1 p. 28 (Ni 0.3575)") : 0.3934;
  const double epsO = R ? R->use("reflection.epsL.O_on_Fe", 0.0340, 0.0306, 0.0374, "1/keV", Tag::model_derived, "thomas_janev_smith_1991", "Eq. (2) with M2=55.85, Z2=26 (Thomas-Fermi; Firsov within 10%)") : 0.0340;
  auto set = [&](int i, int kind, int product, int fragment, int water, double mass_amu, double Erec, double nH, double nO, int rset, double epsL, double share) {
    t[i] = SW7Species{kind, product, fragment, water, mass_amu * amu, Erec, nH, nO, rset, epsL, share};
  };
  set(SW7_E, 0, -1, -1, 0, 9.1093837139e-31 / amu, 0, 0, 0, 0, 0, 0);
  set(SW7_H2O, 2, SW7_H2O, -1, 1, 18.01528, 0, 2, 1, 1, epsO, 15.9994 / 18.01528);
  set(SW7_H2OP, 1, SW7_H2O, -1, 1, 18.01473, ie("H2O+.recombination", 12.621, 12.619, 12.623, "nist_webbook", "IE(H2O) 12.621 +- 0.002 eV (Lias evaluation)"), 2, 1, 1, epsO, 15.9994 / 18.01528);
  set(SW7_H3OP, 1, SW7_H2O, SW7_H, 1, 19.02267, ie("H3O+.recombination", 6.4, 4.3, 6.44, "neau_2000", "H3O+ + e -> H2O + H exothermicity 6.4 eV (J. Chem. Phys. 113, 1762); IE(H)-PA(H2O)=6.44; H3O radical IE 4.3 (Gellene & Porter 1984) as lower band"), 3, 1, 1, epsO, 15.9994 / 19.02322);
  set(SW7_OHP, 1, SW7_OH, -1, 0, 17.00680, ie("OH+.recombination", 13.017, 13.015, 13.019, "nist_webbook", "IE(OH)"), 1, 1, 1, epsO, 15.9994 / 17.00734);
  set(SW7_OP, 1, SW7_O, -1, 0, 15.99885, ie("O+.recombination", 13.618055, 13.618048, 13.618062, "nist_asd", "IE(O)"), 0, 1, 1, epsO, 1.0);
  set(SW7_HP, 1, SW7_H, -1, 0, 1.00728, ie("H+.recombination", 13.598434599702, 13.598434599702, 13.598434599702, "nist_asd", "IE(H)"), 1, 0, 0, epsH, 1.0);
  set(SW7_H2P, 1, SW7_H2, -1, 0, 2.01533, ie("H2+.recombination", 15.42593, 15.42588, 15.42598, "nist_webbook", "IE(H2)"), 2, 0, 0, epsH, 0.5);
  set(SW7_OH, 2, SW7_OH, -1, 0, 17.00734, 0, 1, 1, 1, epsO, 15.9994 / 17.00734);
  set(SW7_O, 2, SW7_O, -1, 0, 15.9994, 0, 0, 1, 1, epsO, 1.0);
  set(SW7_H, 2, SW7_H, -1, 0, 1.00794, 0, 1, 0, 0, epsH, 1.0);
  set(SW7_H2, 2, SW7_H2, -1, 0, 2.01588, 0, 2, 0, 0, epsH, 0.5);
  return t;
}
inline const char* species_name(unsigned i) {
  static const char* n[SW7_NSPECIES] = {"e", "H2O", "H2O+", "H3O+", "OH+", "O+", "H+", "H2+", "OH", "O", "H", "H2"};
  return i < SW7_NSPECIES ? n[i] : "?";
}
inline int species_index(const std::string& s) {
  for (unsigned i = 0; i < SW7_NSPECIES; ++i)
    if (s == species_name(i)) return int(i);
  if (s == "electron" || s == "e-") return SW7_E;
  return -1;
}

// ------------------------------------------------------------------ exact ledger
using i128 = __int128;
struct Sums {   // exact integer sums over wall events of one (patch, weight) class
  std::int64_t events = 0, invalid = 0, energy_capped = 0, yield_clamped = 0, outside_data = 0;
  std::int64_t in_species[SW7_NSPECIES]{}, out_species[SW7_NSPECIES]{};
  std::int64_t out_channel[16]{};
  std::int64_t stuck = 0, electrons_out = 0, charge_in = 0, charge_out = 0;
  i128 k_in_q = 0, e_in_q = 0, e_out_q = 0, wall_q = 0, ads_q = 0, k_out_q = 0, sigma_q = 0;
  std::int64_t words = 0;
  void add(const SW7Hit& h, const SW7Out& o) {
    ++events;
    if (o.status & SW7_ST_INVALID) ++invalid;
    if (o.status & SW7_ST_ENERGY_CAP) ++energy_capped;
    if (o.status & SW7_ST_YIELD_CLAMP) ++yield_clamped;
    if (o.status & SW7_ST_OUTSIDE_DATA) ++outside_data;
    if (h.species < SW7_NSPECIES) ++in_species[h.species];
    for (unsigned i = 0; i < o.n_out; ++i) { ++out_species[o.sp[i]]; ++out_channel[o.ch[i] & 15]; k_out_q += o.kq[i]; }
    stuck += o.stuck; electrons_out += o.electrons_out; charge_in += o.charge_in; charge_out += o.charge_out;
    k_in_q += o.k_in_q; e_in_q += o.e_in_q; e_out_q += o.e_out_q; wall_q += o.wall_q; ads_q += o.ads_q; sigma_q += o.sigma_q;
    words += o.words;
  }
  void add(const Sums& s) {
    events += s.events; invalid += s.invalid; energy_capped += s.energy_capped; yield_clamped += s.yield_clamped; outside_data += s.outside_data;
    for (int i = 0; i < SW7_NSPECIES; ++i) { in_species[i] += s.in_species[i]; out_species[i] += s.out_species[i]; }
    for (int i = 0; i < 16; ++i) out_channel[i] += s.out_channel[i];
    stuck += s.stuck; electrons_out += s.electrons_out; charge_in += s.charge_in; charge_out += s.charge_out;
    k_in_q += s.k_in_q; e_in_q += s.e_in_q; e_out_q += s.e_out_q; wall_q += s.wall_q; ads_q += s.ads_q; k_out_q += s.k_out_q; sigma_q += s.sigma_q;
    words += s.words;
  }
  bool operator==(const Sums& o) const {
    if (events != o.events || invalid != o.invalid || energy_capped != o.energy_capped || yield_clamped != o.yield_clamped || outside_data != o.outside_data) return false;
    for (int i = 0; i < SW7_NSPECIES; ++i) if (in_species[i] != o.in_species[i] || out_species[i] != o.out_species[i]) return false;
    for (int i = 0; i < 16; ++i) if (out_channel[i] != o.out_channel[i]) return false;
    return stuck == o.stuck && electrons_out == o.electrons_out && charge_in == o.charge_in && charge_out == o.charge_out && k_in_q == o.k_in_q &&
           e_in_q == o.e_in_q && e_out_q == o.e_out_q && wall_q == o.wall_q && ads_q == o.ads_q && k_out_q == o.k_out_q && sigma_q == o.sigma_q && words == o.words;
  }
};
inline std::string i128_text(i128 v) {
  if (v == 0) return "0";
  bool neg = v < 0;
  unsigned __int128 u = neg ? (unsigned __int128)(-(v + 1)) + 1 : (unsigned __int128)v;
  std::string s;
  while (u) { s.push_back(char('0' + int(u % 10))); u /= 10; }
  if (neg) s.push_back('-');
  return {s.rbegin(), s.rend()};
}
// Per-event exact identities (integers): energy  e_in = e_out + wall,  e_out = sum kq + n_e phi_q - ads_q,
// electron e_in = k_in + phi_q, neutral e_in = k_in; charge  charge_out = -n_e; counts n_out <= MAXOUT.
inline bool event_closes(const SW7Hit& h, const SW7Out& o, const std::array<SW7Species, SW7_NSPECIES>& table) {
  if (o.status & (SW7_ST_INVALID | SW7_ST_RNG_BUDGET)) return o.n_out <= SW7_MAXOUT;
  if (o.n_out > SW7_MAXOUT || o.e_in_q != o.e_out_q + o.wall_q) return false;
  std::int64_t out = -o.ads_q, ne = 0;
  for (unsigned i = 0; i < o.n_out; ++i) {
    if (o.kq[i] < 0 || o.sp[i] >= SW7_NSPECIES) return false;
    out += o.kq[i];
    if (o.sp[i] == SW7_E) { out += o.phi_q; ++ne; }
  }
  if (o.e_out_q != out || ne != std::int64_t(o.electrons_out) || o.charge_out != -int(ne)) return false;
  const auto& sp = table[h.species];
  if (sp.kind == 0) return o.e_in_q == o.k_in_q + o.phi_q && o.charge_in == -1;
  if (sp.kind == 2) return o.e_in_q == o.k_in_q && o.charge_in == 0;
  return o.e_in_q == o.k_in_q + sw7::quanta(sp.recombination_eV) - o.phi_q && o.charge_in == 1;
}

// ------------------------------------------------------------------ interval bound (API-v7.txt 3)
struct IntervalPolicy {
  double cap_s = 1e-14;                 // sealed trial interval bound (1e-14: v6.1)
  double emission_relative_error = 0;   // per-step relative error target of expected emission counts; 0: no bound (v6.1)
  double desorption_relative_error = 0; // same for desorption (thermal + stimulated) releases
  double safety = 0.5;                  // factor on the first-order bound (curvature margin)
  double count_floor = 1e-3;            // expected count per step below which a stream does not constrain dt
  bool operator==(const IntervalPolicy&) const = default;
};
struct RateSample {   // one accepted state, supplied by the owner (rate_sample) to the caller's history
  double t_s = 0, dt_s = 0;
  std::vector<double> log_emission;     // per (patch, class): log expected births per second (-inf: none)
  std::vector<double> log_desorption;   // per patch: log expected thermal release per second
  std::vector<double> T_K;              // per patch
};
struct IntervalDecision {
  double bound_s = INFINITY;
  std::string limiting = "none";
  double rate_per_s = 0;                // |d ln R / dt| of the limiting quantity
  std::size_t index = 0;
};
// First-order explicit exposure: relative error of a step's expected count ~ dt |d ln R/dt| / 2
// (observed order 1.00, coordinator v7 calibration). dt <= 2 eps safety / |d ln R/dt|.
inline IntervalDecision interval_bound(const RateSample& a, const RateSample& b, const IntervalPolicy& p) {
  IntervalDecision d;
  d.bound_s = p.cap_s;
  d.limiting = "cap";
  const double span = b.t_s - a.t_s;
  if (!(span > 0) || a.log_emission.size() != b.log_emission.size() || a.log_desorption.size() != b.log_desorption.size()) return d;
  auto consider = [&](const std::vector<double>& x, const std::vector<double>& y, double eps, const char* what) {
    if (!(eps > 0)) return;
    for (std::size_t i = 0; i < x.size(); ++i) {
      if (!std::isfinite(x[i]) || !std::isfinite(y[i])) continue;
      const double n_step = std::exp(std::max(x[i], y[i])) * b.dt_s;
      if (!(n_step >= p.count_floor)) continue;
      const double r = std::abs(y[i] - x[i]) / span;
      if (!(r > 0)) continue;
      const double bound = 2 * eps * p.safety / r;
      if (bound < d.bound_s) { d.bound_s = bound; d.limiting = what; d.rate_per_s = r; d.index = i; }
    }
  };
  consider(a.log_emission, b.log_emission, p.emission_relative_error, "emission");
  consider(a.log_desorption, b.log_desorption, p.desorption_relative_error, "desorption");
  return d;
}

}  // namespace surface_v7
