#pragma once
#include "high_energy_bands.hpp"

// Product kernel parent_ground_BEB_SDCS_binary_encounter_v2 for e + H2O ->
// e + e + H2O+(ground) above 1 keV. Pure function: no RNG, no inventories.
//
// Frame conventions as in the v1 kernel parent_ground_equal_pair_CM_v1: the
// input is the incident kinetic energy T and unit direction in the target rest
// frame (the stationary-target laboratory); the output momenta are in that
// frame and the host boosts them exactly to the laboratory; the ion is the
// ground H2O+ of mass M - m_e + I/c^2 with the chemical reference I = 12.621 eV
// counted once; relativistic energy and momentum are conserved exactly.
//
// Energy sharing (replaces the equal-pair split above 1 keV):
//  1. orbital i in {1b1, 3a1, 1b2} with probability w_i sigma_i(T) / sum_j
//     w_j sigma_j(T), using the band's own orbital rate terms (BEB in
//     1-10 keV, RBEB above 10 keV);
//  2. ejected (slower) electron kinetic energy W = B_i w, w sampled from the
//     orbital's singly differential cross section on 0 <= w <= (t-1)/2 by
//     exact inversion of its analytic integral (beb_rbeb.hpp); the BEB SDCS
//     in the 1-10 keV band, the RBEB SDCS above 10 keV;
//  3. ejected-electron direction from the relativistic binary-encounter
//     relation cos^2(theta) = W (T + 2mc^2) / (T (W + 2mc^2)) about the
//     incident direction, azimuth 2 pi u_azimuth;
//  4. the scattered electron and the ion share the remaining four-momentum
//     (T - W, p0 - p_ej) as an exact two-body system of invariant mass M_r:
//     in its rest frame they are back-to-back along the remaining momentum,
//     the electron forward. With the binary-encounter angle the free-electron
//     limit leaves the ion at rest; the ion recoil carries only the momentum
//     mismatch from the binding energy (order I/v).
// The scattered electron then has T - W - I - K_ion. The orbital binding
// B_i - I above the ion ground state is not retained as ion excitation (the
// same ground-ion convention as outer3_BEB_Tan_parent_v1 and v1).
// SDCS energies and the binary-encounter angle are target-rest quantities,
// the frame in which BEB/RBEB are defined. The angular law is a declared
// approximation; the SDCS supplies only the energy distribution.
namespace raikou::native_water::coverage_v2 {

enum class Band { outer3_beb, provisional_rbeb };
inline const char* name(Band b) { return b == Band::outer3_beb ? "outer3_BEB" : "provisional_RBEB"; }

struct SdcsParentProducts {
  std::array<double, 3> primary_momentum_kg_m_s{}, secondary_momentum_kg_m_s{}, ion_momentum_kg_m_s{};
  double primary_kinetic_J = 0, secondary_kinetic_J = 0, ion_kinetic_J = 0;
  double ion_mass_kg = 0, chemical_reference_gain_J = 0, cross_section_m2 = 0;
  std::size_t orbital = 0;
  std::string orbital_label;
  double ejected_eV = 0, cos_theta_ejected = 0;
  // Internal target-rest closure in extended precision (diagnostic).
  long double energy_defect_eV = 0, momentum_defect_eV_c = 0;
};

namespace kernel_detail {
using V = std::array<long double, 3>;
inline V basis_perpendicular(const V& n, V& e2) {
  // Deterministic orthonormal pair perpendicular to n.
  const int k = std::abs(n[0]) <= std::abs(n[1]) && std::abs(n[0]) <= std::abs(n[2]) ? 0 : (std::abs(n[1]) <= std::abs(n[2]) ? 1 : 2);
  V a{0, 0, 0};a[k] = 1;
  V e1{a[1] * n[2] - a[2] * n[1], a[2] * n[0] - a[0] * n[2], a[0] * n[1] - a[1] * n[0]};
  const long double s = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  for (auto& x : e1) x /= s;
  e2 = {n[1] * e1[2] - n[2] * e1[1], n[2] * e1[0] - n[0] * e1[2], n[0] * e1[1] - n[1] * e1[0]};
  return e1;
}
} // namespace kernel_detail

// Kinematics for a given ejected energy W (eV) and azimuth variate. Exposed
// for deterministic tests; the sampling wrapper is below.
inline SdcsParentProducts beb_sdcs_parent_kinematics(double T_eV, double W_eV, double neutral_mass_kg,
                                                     const std::array<double, 3>& incident_unit, long double u_azimuth) {
  using namespace kernel_detail;
  const auto require = [](bool ok, const char* what) { if (!ok) throw std::invalid_argument(what); };
  require(std::isfinite(T_eV) && T_eV > 0 && std::isfinite(W_eV) && W_eV >= 0 && W_eV < T_eV, "kernel energies");
  require(std::isfinite(neutral_mass_kg) && neutral_mass_kg > 1000 * w::me, "kernel neutral mass");
  require(u_azimuth >= 0 && u_azimuth < 1, "azimuth variate");
  long double nn = 0;for (double x : incident_unit) { require(std::isfinite(x), "incident axis");nn += static_cast<long double>(x) * x; }
  require(std::abs(nn - 1) < 1e-12L, "incident direction must be a unit vector");
  const long double qe = pc::qe, c = pc::c, I = w::threshold_eV;
  const long double me = static_cast<long double>(w::me) * c * c / qe;              // eV
  const long double M = static_cast<long double>(neutral_mass_kg) * c * c / qe;      // eV
  const long double Mion = M - me + I;                                              // eV, ground H2O+
  const long double K = T_eV, W = W_eV;
  const long double p0 = std::sqrt(K * (K + 2 * me)), p2 = std::sqrt(W * (W + 2 * me)); // eV/c
  const V n{incident_unit[0], incident_unit[1], incident_unit[2]};
  V e2{};const V e1 = basis_perpendicular(n, e2);
  // Binary-encounter polar angle; sin^2 written without cancellation.
  const long double den = K * (W + 2 * me);
  const long double cos_t = std::sqrt(W * (K + 2 * me) / den), sin_t = std::sqrt(2 * me * (K - W) / den);
  const long double phi = 2 * std::numbers::pi_v<long double> * u_azimuth;
  const long double cp = std::cos(phi), sp = std::sin(phi);
  V d2{};for (int i = 0; i < 3; ++i) d2[i] = cos_t * n[i] + sin_t * (cp * e1[i] + sp * e2[i]);
  // Remaining two-body system: energy above M, momentum P_r.
  const long double D = K - W;
  const long double par = p0 - p2 * cos_t, perp = p2 * sin_t;
  const long double Pr = std::sqrt(par * par + perp * perp);
  V dr{};for (int i = 0; i < 3; ++i) dr[i] = (p0 * n[i] - p2 * d2[i]) / Pr;
  const long double Sm = M + I; // m_e + M_ion
  const long double X = (D - I) * (2 * M + D + I) - Pr * Pr; // (M_r^2 - Sm^2) c^4
  require(X > 0, "ejected energy leaves no invariant mass for the scattered electron and ion");
  const long double Mr = Sm + X / (std::sqrt(Sm * Sm + X) + Sm);
  const long double Y = X + 4 * Mion * me;
  const long double ps = std::sqrt(X * Y) / (2 * Mr);
  const long double E1s = std::sqrt(ps * ps + me * me), Eis = std::sqrt(ps * ps + Mion * Mion);
  const long double g = (M + D) / Mr, gb = Pr / Mr;
  const long double p1 = g * ps + gb * E1s, pion = -g * ps + gb * Eis;
  const long double kin_e = me, kin_i = Mion;
  auto kinetic = [](long double p, long double m) { return p * p / (std::sqrt(p * p + m * m) + m); };
  SdcsParentProducts out;
  const long double to_SI = qe / c;
  V P1{}, P2{}, Pi{};
  for (int i = 0; i < 3; ++i) { P1[i] = p1 * dr[i];P2[i] = p2 * d2[i];Pi[i] = pion * dr[i]; }
  for (int i = 0; i < 3; ++i) {
    out.primary_momentum_kg_m_s[i] = static_cast<double>(P1[i] * to_SI);
    out.secondary_momentum_kg_m_s[i] = static_cast<double>(P2[i] * to_SI);
    out.ion_momentum_kg_m_s[i] = static_cast<double>(Pi[i] * to_SI);
  }
  const long double K1 = kinetic(p1, kin_e), K2 = kinetic(p2, kin_e), Ki = kinetic(std::abs(pion), kin_i);
  out.primary_kinetic_J = static_cast<double>(K1 * qe);out.secondary_kinetic_J = static_cast<double>(K2 * qe);
  out.ion_kinetic_J = static_cast<double>(Ki * qe);
  out.ion_mass_kg = neutral_mass_kg - w::me + w::threshold_eV * pc::qe / (pc::c * pc::c);
  out.chemical_reference_gain_J = w::threshold_eV * pc::qe;
  out.ejected_eV = W_eV;out.cos_theta_ejected = static_cast<double>(cos_t);
  out.energy_defect_eV = (K1 + K2 + Ki + I) - K;
  long double dp = 0;
  for (int i = 0; i < 3; ++i) { const long double r = P1[i] + P2[i] + Pi[i] - p0 * n[i];dp += r * r; }
  out.momentum_defect_eV_c = std::sqrt(dp);
  return out;
}

// Sampling wrapper. u_orbital, u_energy in (0,1), u_azimuth in [0,1).
inline SdcsParentProducts beb_sdcs_parent_products(const HighEnergyBands& bands, Band band, double T_eV,
                                                   double neutral_mass_kg, const std::array<double, 3>& incident_unit,
                                                   long double u_orbital, long double u_energy, long double u_azimuth) {
  const auto terms = band == Band::outer3_beb ? bands.outer3_terms(T_eV) : bands.provisional_terms(T_eV);
  const std::size_t i = beb_v2::select_orbital(terms, u_orbital);
  const auto& o = bands.orbitals()[i];
  const auto shape = band == Band::outer3_beb ? beb_v2::beb_shape(o, bands.constants(), T_eV)
                                              : beb_v2::rbeb_shape(o, bands.constants(), T_eV);
  const long double wr = beb_v2::sample_reduced_w(shape, u_energy);
  const double W = static_cast<double>(wr * shape.B);
  auto out = beb_sdcs_parent_kinematics(T_eV, W, neutral_mass_kg, incident_unit, u_azimuth);
  out.orbital = i;out.orbital_label = o.label;
  out.cross_section_m2 = band == Band::outer3_beb ? bands.outer3_parent_m2(T_eV) : bands.provisional_parent_m2(T_eV);
  return out;
}
// Event bridge for the two high-energy bands: candidate re-validation against
// the composite channel, target-rest kernel, exact boost to the laboratory and
// the same weighted energy/momentum receipts as the sealed v1 bridge.
struct HighEnergySample {
  w::Products products;
  SdcsParentProducts rest;
  Band band = Band::outer3_beb;
};
inline HighEnergySample high_energy_parent_sample(const HighEnergyBands& bands, Band band, const pc::Channel& composite,
                                                  double neutral_mass_kg, const pc::Candidate& event,
                                                  const pc::Packet& electron, const pc::Packet& target,
                                                  long double u_orbital, long double u_energy, long double u_azimuth) {
  pc::require(electron.state == "electron" && target.state == "H2O" && electron.mass_kg == w::me &&
                  target.mass_kg == neutral_mass_kg, "high-energy parent bridge species/catalog mass mismatch");
  pc::require(event.projectile_id == electron.id && event.target_id == target.id && electron.cell == target.cell &&
                  event.cell == electron.cell, "stale high-energy parent identity/cell");
  pc::require(event.projectile_before == electron.weight && event.target_before == target.weight, "stale high-energy parent weights");
  pc::require(event.channel_id == composite.id && event.data_sha256 == composite.data_sha256 &&
                  event.kernel_sha256 == composite.kernel_sha256, "foreign high-energy event data/kernel");
  const auto expected = pc::candidates({electron}, {target}, {{event.cell, event.volume_m3}}, {composite}, 1, 1, 1);
  pc::require(expected.size() == 1, "no positive supported high-energy parent event");
  const auto& a = expected[0];
  pc::require(event.event_weight == a.event_weight && event.energy_eV == a.energy_eV && event.rate_s == a.rate_s &&
                  event.flux_m_s == a.flux_m_s && event.sigma_m2 == a.sigma_m2, "stale or inconsistent high-energy rate candidate");
  pc::require(event.projectile_remainder.high == a.projectile_remainder.high && event.projectile_remainder.low == a.projectile_remainder.low &&
                  event.target_remainder.high == a.target_remainder.high && event.target_remainder.low == a.target_remainder.low,
              "changed compensated remainder");
  const bool in_beb = a.energy_eV > bands.beb_low_eV() && a.energy_eV <= bands.beb_high_eV();
  pc::require(band == (in_beb ? Band::outer3_beb : Band::provisional_rbeb), "high-energy band differs from the candidate energy");
  const auto direction = w::incident_direction(electron, target, a.energy_eV);
  HighEnergySample s;s.band = band;
  s.rest = beb_sdcs_parent_products(bands, band, a.energy_eV, target.mass_kg, direction, u_orbital, u_energy, u_azimuth);
  pc::require(s.rest.cross_section_m2 == a.sigma_m2, "kernel and candidate cross sections differ");
  auto& out = s.products;
  out.projectile_remainder = a.projectile_remainder;out.target_remainder = a.target_remainder;
  out.target_frame_energy_eV = a.energy_eV;out.cross_section_m2 = s.rest.cross_section_m2;
  out.particles = {w::Product{"electron", -1, w::me, a.event_weight, 0, {}}, w::Product{"electron", -1, w::me, a.event_weight, 0, {}},
                   w::Product{"H2O+", 1, s.rest.ion_mass_kg, a.event_weight, s.rest.chemical_reference_gain_J, {}}};
  const std::array<std::array<double, 3>, 3> momenta{s.rest.primary_momentum_kg_m_s, s.rest.secondary_momentum_kg_m_s,
                                                     s.rest.ion_momentum_kg_m_s};
  long double before = w::kinetic(electron.mass_kg, electron.proper_velocity_m_s) + w::kinetic(target.mass_kg, target.proper_velocity_m_s),
              after = s.rest.chemical_reference_gain_J;
  std::array<long double, 3> dp{}, abs_momentum{};
  for (int d = 0; d < 3; ++d) {
    dp[d] = -static_cast<long double>(electron.mass_kg) * electron.proper_velocity_m_s[d] -
            static_cast<long double>(target.mass_kg) * target.proper_velocity_m_s[d];
    abs_momentum[d] = std::abs(static_cast<long double>(electron.mass_kg) * electron.proper_velocity_m_s[d]) +
                      std::abs(static_cast<long double>(target.mass_kg) * target.proper_velocity_m_s[d]);
  }
  for (int i = 0; i < 3; ++i) {
    auto& product = out.particles[i];
    w::Vec u{};for (int d = 0; d < 3; ++d) u[d] = momenta[i][d] / product.mass_kg;
    product.proper_velocity_m_s = w::to_lab(u, target.proper_velocity_m_s);
    after += w::kinetic(product.mass_kg, product.proper_velocity_m_s);
    for (int d = 0; d < 3; ++d) {
      const long double p = static_cast<long double>(product.mass_kg) * product.proper_velocity_m_s[d];
      dp[d] += p;abs_momentum[d] += std::abs(p);
    }
  }
  const long double factor = 256 * std::numeric_limits<double>::epsilon();
  const long double energy_bound = factor * (std::abs(after) + std::abs(before));
  pc::require(std::abs(after - before) <= energy_bound, "unweighted high-energy parent lab energy closure");
  pc::require(energy_bound <= w::maximum_lab_roundoff_fraction_of_reaction_energy * a.energy_eV * pc::qe,
              "lab representation cannot resolve reaction energy");
  out.energy_defect_J = static_cast<double>((after - before) * a.event_weight);
  out.energy_roundoff_limit_J = static_cast<double>(energy_bound * a.event_weight);
  long double pscale = 0;
  for (int d = 0; d < 3; ++d) { out.momentum_defect_kg_m_s[d] = static_cast<double>(dp[d] * a.event_weight);pscale += abs_momentum[d]; }
  for (long double v : dp) pc::require(std::abs(v) <= factor * pscale, "unweighted high-energy parent lab momentum closure");
  out.momentum_roundoff_limit_kg_m_s = static_cast<double>(factor * pscale * a.event_weight);
  pc::require(std::isfinite(out.energy_defect_J) && std::isfinite(out.energy_roundoff_limit_J) && out.energy_roundoff_limit_J > 0 &&
                  std::abs(out.energy_defect_J) <= out.energy_roundoff_limit_J, "high-energy parent weighted energy receipt");
  pc::require(std::isfinite(out.momentum_roundoff_limit_kg_m_s) && out.momentum_roundoff_limit_kg_m_s > 0,
              "high-energy parent weighted momentum receipt");
  for (double v : out.momentum_defect_kg_m_s)
    pc::require(std::isfinite(v) && std::abs(v) <= out.momentum_roundoff_limit_kg_m_s, "high-energy parent lab momentum closure");
  return s;
}
} // namespace raikou::native_water::coverage_v2
