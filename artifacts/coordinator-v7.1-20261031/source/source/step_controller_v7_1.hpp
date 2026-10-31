#pragma once
// Coordinator v7 step controller: chooses each step's dt from the accepted state.
//
// Criteria (each a predicate on a candidate dt evaluated with + - * / sqrt only,
// so the same inputs give the same decision bits on any IEEE-754 machine):
//   plasma        omega_p,max dt <= c_plasma. omega_p^2 per Hall node from the actual
//                 charge density: the Hall motion guard's own operand (point deposit
//                 of w q^2/(eps0 m) of every charged carrier, divided by the node
//                 volume); per node the terms are summed in ascending order of value,
//                 so the sum does not depend on how the census is partitioned.
//   gyro          omega_c,max dt <= c_gyro, omega_c = |q/m| |B| with |B| bounded per
//                 cell by the face maxima of the accepted B over the 3x3x3 cell
//                 neighbourhood of every charged carrier's cell.
//   displacement  dt min(c, v_max + a_max dt) <= c_disp h_min (v = |u|/gamma of every
//                 carrier, a = |q/m| |E| with the same neighbourhood bound of E).
//   chemistry     predicted unresolved upper-bound parent events per exposure (rate
//                 of the last accepted step's exposures times dt/2) <= c_chem times
//                 the coverage policy's per-exposure limit.
//   drive         |I(t+dt) - I(t)| <= c_drive max(|I(t)|, drive_floor) (C19 drive).
//   lte           particle-free steps only: trapezoidal local truncation estimate
//                 (dt^3/12) max|x'''|/max|x| over the accepted E and B, x''' = 6 times
//                 the third divided difference of the last four accepted states,
//                 <= lte_tol.
//   births        particle-free steps only: expected charged births in the step (rate
//                 of the last step) <= birth_threshold, otherwise dt <= particle cap.
// Bounds: dt <= cap (configured; <= the field owner's admitted maximum), dt <=
// particle_cap while charged carriers exist, dt <= drive end - t. dt is the largest
// rung of the ladder particle_cap * 2^(-j/8) (j may be negative, up to the largest rung <= cap)
// admitted by every criterion; growth by at most
// growth_rungs rungs per step (8 rungs = x2), a rung above the previous one needs
// every criterion at grow_margin * dt; after a refused attempt no growth for
// hold_after_retry accepted steps. Optional landing times (sync_every_s): a step
// that would pass a landing time ends exactly on it.
// Guard refusal: an attempt refused by a step-size guard of an owner is retried
// retry_rungs rungs lower inside the same step transaction (owners are unchanged by
// a refused preparation); every attempt is recorded.
#include "joint_step_v7_1.hpp"
#include <algorithm>
#include <cstring>

namespace raikou::native_water::v7_1 {

struct ControllerConfig {
  bool adaptive = false;          // false: fixed dt (fixed_dt_s), no retries
  bool fixed_grid_retry = false;  // fixed dt only: a step refused by a step-size guard is retried at half its
                                  // length and the next step completes the grid interval (calibration runs)
  double fixed_dt_s = 1e-14;
  double cap_s = 1e-14;           // configured step cap (<= field owner maximum)
  double particle_cap_s = 1e-14;  // cap while carriers exist (see neutrals_need_particle_cap) or charged births are expected
  // A resident neutral carrier moves through the sealed Hall 8dddec0 neutral stream (kinetic_boundary.cpp,
  // prepare_neutral_stream), which requires dt <= 1e-14: while that holds, a step with any carrier is a
  // particle step. false: only charged carriers make a particle step.
  bool neutrals_need_particle_cap = true;
  // Calibrated defaults (README.txt, RESULTS): c_plasma 0.03 is below the onset of Hall midpoint fixed-point
  // refusals (observed from omega_p dt = 0.04) with a frequency error (omega_p dt)^2/12 <= 7.5e-5 per radian.
  double c_plasma = .03, c_gyro = .025, c_disp = .1, c_chem = .5, c_drive = 1e-3, drive_floor_A = 1e3;
  double lte_tol = 1e-8, birth_threshold = 1e-3;
  int growth_rungs = 8, hold_after_retry = 4, retry_rungs = 8, max_retries = 6;
  double grow_margin = 1.05, dt_min_s = 1e-18;
  double sync_every_s = 0, sync_origin_s = 0;
  std::string text() const {
    std::ostringstream s;
    s << std::hexfloat << "coordinator-v7-step-controller-v1;adaptive=" << adaptive << ";fixed_grid_retry=" << fixed_grid_retry << ";fixed_dt=" << fixed_dt_s << ";cap=" << cap_s
      << ";particle_cap=" << particle_cap_s << ";neutrals_need_particle_cap=" << neutrals_need_particle_cap << ";c_plasma=" << c_plasma << ";c_gyro=" << c_gyro << ";c_disp=" << c_disp << ";c_chem=" << c_chem
      << ";c_drive=" << c_drive << ";drive_floor=" << drive_floor_A << ";lte_tol=" << lte_tol << ";birth_threshold=" << birth_threshold
      << ";ladder=cap*2^(-j/8);growth_rungs=" << growth_rungs << ";grow_margin=" << grow_margin << ";hold_after_retry=" << hold_after_retry
      << ";retry_rungs=" << retry_rungs << ";max_retries=" << max_retries << ";dt_min=" << dt_min_s << ";sync_every=" << sync_every_s
      << ";sync_origin=" << sync_origin_s;
    return s.str();
  }
  std::string sha256() const { return sha(text()); }
  void validate(double owner_max_s) const {
    need(std::isfinite(cap_s) && cap_s > 0 && cap_s <= owner_max_s, "controller cap must lie in (0, largest step of the coupled owners]");
    need(std::isfinite(particle_cap_s) && particle_cap_s > 0 && particle_cap_s <= cap_s, "particle cap must lie in (0, cap]");
    need(std::isfinite(fixed_dt_s) && fixed_dt_s > 0 && fixed_dt_s <= cap_s, "fixed dt must lie in (0, cap]");
    for (double c : {c_plasma, c_gyro, c_disp, c_chem, c_drive, lte_tol, birth_threshold, drive_floor_A, dt_min_s})
      need(std::isfinite(c) && c > 0, "controller factors must be positive and finite");
    need(c_plasma < .1 && c_gyro < .05 && c_disp < .2 && c_chem <= 1, "controller factors must lie below the owners' guards");
    need(growth_rungs >= 0 && hold_after_retry >= 0 && retry_rungs > 0 && max_retries >= 0 && grow_margin >= 1, "controller growth/retry settings");
    need(std::isfinite(sync_every_s) && sync_every_s >= 0 && std::isfinite(sync_origin_s), "controller landing times");
  }
};

// 2^(-r/8), r = 0..7 (fixed binary64 constants).
inline constexpr double ladder_mantissa[8] = {0x1p+0, 0x1.d5818dcfba487p-1, 0x1.ae89f995ad3adp-1, 0x1.8ace5422aa0dbp-1,
                                              0x1.6a09e667f3bcdp-1, 0x1.4bfdad5362a27p-1, 0x1.306fe0a31b715p-1, 0x1.172b83c7d517bp-1};
// Rung j of the ladder anchored at a (j may be negative: rungs above the anchor).
inline double rung_dt(double anchor, int j) {
  const int q = j >= 0 ? j / 8 : -((-j + 7) / 8);
  return std::ldexp(anchor * ladder_mantissa[j - 8 * q], -q);
}
// The ladder is anchored at the particle cap (so 1e-14 stays a rung when the cap is larger);
// the top rung is the largest rung <= cap (ladder_dt, top_rung below).

inline double ladder_dt(const ControllerConfig& C, int j) { return rung_dt(C.particle_cap_s, j); }
inline int top_rung(const ControllerConfig& C) {
  int j = 0;
  while (rung_dt(C.particle_cap_s, j - 1) <= C.cap_s) --j;
  return j;
}
enum class Criterion : int { fixed = 0, cap, growth, hold, particle_cap, plasma, gyro, displacement, chemistry, drive, lte, births, sync, drive_end, retry };
inline const char* criterion_name(Criterion c) {
  static const char* n[] = {"fixed", "cap", "growth", "hold", "particle_cap", "plasma", "gyro", "displacement", "chemistry", "drive", "lte",
                            "births", "sync", "drive_end", "retry"};
  return n[int(c)];
}

// Inputs of one decision (accepted start-of-step state and the last accepted step).
struct ControlInputs {
  double t_s = 0, dt_prev_s = 0;
  int rung_prev = 0, hold = 0;
  std::uint64_t carriers = 0, charged = 0, electrons = 0;
  double omega_p_max = 0, omega_c_max = 0, v_max = 0, a_max = 0, h_min = 0;
  double chem_rate = 0, chem_limit = .01;        // unresolved upper-bound events per second of exposure; per-exposure limit
  double drive_I = 0, drive_end_s = 0;
  double lte_coeff = 0;                           // max|x'''|/max|x|
  bool lte_valid = false;
  double births_rate = 0;                         // expected charged births per second (last step)
  std::uint64_t plasma_terms = 0;
};

struct Decision {
  double dt_s = 0;
  int rung = 0;
  bool landing = false, refused = false;
  Criterion criterion = Criterion::fixed;
  double value = 0, limit = 0;                    // guard quantity of the limiting criterion at dt, and its limit
  // continuous bounds per criterion (diagnostic; +inf when not applicable)
  double bound_plasma = INFINITY, bound_gyro = INFINITY, bound_disp = INFINITY, bound_chem = INFINITY, bound_drive = INFINITY,
         bound_lte = INFINITY;
  std::string sha256;
};

// ------------------------------------------------------------ deterministic inputs
struct NodeTerm {
  std::int64_t node;
  double value;
};
inline bool operator<(const NodeTerm& a, const NodeTerm& b) { return a.node != b.node ? a.node < b.node : a.value < b.value; }
// Hall motion plasma operand of one carrier, as the eight CIC node terms of
// axial::point_deposit (same arithmetic, same order of factors).
inline void plasma_terms(const axial::Mesh& g, const m::PicVector& shift, const m::GapParticle& p, std::vector<NodeTerm>& out) {
  if (p.charge_C == 0) return;
  const double load = p.weight * p.charge_C * p.charge_C / (axial::eps * p.mass_kg);
  const double spacing[3] = {g.dx, g.dy, g.dz}, origin[3] = {-g.side / 2, -g.side / 2, g.zmin};
  const int counts[3] = {g.nx, g.ny, g.nz};
  int cell[3];
  double a[3];
  for (int d = 0; d < 3; ++d) {
    double v = ((p.position_m[std::size_t(d)] + shift[std::size_t(d)]) - origin[d]) / spacing[d];
    v = std::clamp(v, 0., double(counts[d]));
    cell[d] = std::min(int(std::floor(v)), counts[d] - 1);
    a[d] = v - cell[d];
  }
  for (int k : {0, 1})
    for (int j : {0, 1})
      for (int i : {0, 1})
        out.push_back({g.node(cell[0] + i, cell[1] + j, cell[2] + k), load * (i ? a[0] : 1 - a[0]) * (j ? a[1] : 1 - a[1]) * (k ? a[2] : 1 - a[2])});
}
// max over nodes of sqrt(sum/VN), each node's terms summed in ascending order of value.
// v7.1: terms are bucketed by node (counting sort, O(terms + nodes)) and only each node's
// terms are sorted; the order within a node, hence every sum and the maximum, are those of
// the global (node, value) sort of coordinator v7.
inline double omega_p_max(const axial::Mesh& g, const std::vector<NodeTerm>& terms) {
  if (terms.empty()) return 0.;
  std::vector<std::uint32_t> start(std::size_t(g.nn) + 1, 0);
  for (const auto& t : terms) ++start[std::size_t(t.node) + 1];
  for (std::size_t n = 0; n < std::size_t(g.nn); ++n) start[n + 1] += start[n];
  std::vector<double> v(terms.size());
  {
    auto at = start;
    for (const auto& t : terms) v[at[std::size_t(t.node)]++] = t.value;
  }
  double best = 0;
  for (std::size_t n = 0; n < std::size_t(g.nn); ++n) {
    const auto a = v.begin() + start[n], b = v.begin() + start[n + 1];
    if (a == b) continue;
    std::sort(a, b);
    double s = 0;
    for (auto it = a; it != b; ++it) s += *it;
    best = std::max(best, std::sqrt(s / g.VN[n]));
  }
  return best;
}
// Per-cell component maxima of B (faces) and E (edges) over the 3x3x3 neighbourhood.
struct CellBounds {
  int nx = 0, ny = 0, nz = 0;
  std::vector<double> B, E;   // |B| and |E| bounds per cell
  bool legacy_27_loop = false;   // --host-work v7: the coordinator v7 27-cell loop (same values)
  void build(const axial::Mesh& g, const hc::Vec& Ef, const hc::Vec& Bf) {
    nx = g.nx; ny = g.ny; nz = g.nz;
    const std::size_t nc = std::size_t(nx) * ny * nz;
    std::vector<double> bx(nc), by(nc), bz(nc), ex(nc), ey(nc), ez(nc);
    auto c = [&](int i, int j, int k) { return std::size_t(i) + std::size_t(nx) * (std::size_t(j) + std::size_t(ny) * std::size_t(k)); };
    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
          const auto x = c(i, j, k);
          bx[x] = std::max(std::abs(Bf[std::size_t(g.bx(i, j, k))]), std::abs(Bf[std::size_t(g.bx(i + 1, j, k))]));
          by[x] = std::max(std::abs(Bf[std::size_t(g.by(i, j, k))]), std::abs(Bf[std::size_t(g.by(i, j + 1, k))]));
          bz[x] = std::max(std::abs(Bf[std::size_t(g.bz(i, j, k))]), std::abs(Bf[std::size_t(g.bz(i, j, k + 1))]));
          double mx = 0, my = 0, mz = 0;
          for (int b : {0, 1})
            for (int a : {0, 1}) {
              mx = std::max(mx, std::abs(Ef[std::size_t(g.ex(i, j + a, k + b))]));
              my = std::max(my, std::abs(Ef[std::size_t(g.ey(i + a, j, k + b))]));
              mz = std::max(mz, std::abs(Ef[std::size_t(g.ez(i + a, j + b, k))]));
            }
          ex[x] = mx; ey[x] = my; ez[x] = mz;
        }
    if (legacy_27_loop) {
      B.assign(nc, 0); E.assign(nc, 0);
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
          for (int i = 0; i < nx; ++i) {
            double b[3] = {0, 0, 0}, e[3] = {0, 0, 0};
            for (int dk = -1; dk <= 1; ++dk)
              for (int dj = -1; dj <= 1; ++dj)
                for (int di = -1; di <= 1; ++di) {
                  const int I = i + di, Jj = j + dj, K = k + dk;
                  if (I < 0 || Jj < 0 || K < 0 || I >= nx || Jj >= ny || K >= nz) continue;
                  const auto y = c(I, Jj, K);
                  b[0] = std::max(b[0], bx[y]); b[1] = std::max(b[1], by[y]); b[2] = std::max(b[2], bz[y]);
                  e[0] = std::max(e[0], ex[y]); e[1] = std::max(e[1], ey[y]); e[2] = std::max(e[2], ez[y]);
                }
            B[c(i, j, k)] = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
            E[c(i, j, k)] = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
          }
      return;
    }
    // 3x3x3 neighbourhood maxima per component as three 1-D passes (exact: maxima commute);
    // the same values as the 27-cell loop of coordinator v7.
    const auto box_max = [&](std::vector<double>& a) {
      std::vector<double> t(nc);
      for (int pass = 0; pass < 3; ++pass) {
        for (int k = 0; k < nz; ++k)
          for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
              const int n = pass == 0 ? nx : pass == 1 ? ny : nz, x = pass == 0 ? i : pass == 1 ? j : k;
              const std::size_t stride = pass == 0 ? 1 : pass == 1 ? std::size_t(nx) : std::size_t(nx) * std::size_t(ny);
              const std::size_t at = c(i, j, k);
              double m = a[at];
              if (x > 0) m = std::max(m, a[at - stride]);
              if (x < n - 1) m = std::max(m, a[at + stride]);
              t[at] = m;
            }
        a.swap(t);
      }
    };
    for (auto* v : {&bx, &by, &bz, &ex, &ey, &ez}) box_max(*v);
    B.assign(nc, 0); E.assign(nc, 0);
    for (std::size_t x = 0; x < nc; ++x) {
      B[x] = std::sqrt(bx[x] * bx[x] + by[x] * by[x] + bz[x] * bz[x]);
      E[x] = std::sqrt(ex[x] * ex[x] + ey[x] * ey[x] + ez[x] * ez[x]);
    }
  }
  std::size_t cell_of(const axial::Mesh& g, const m::PicVector& shift, const m::PicVector& x) const {
    const double spacing[3] = {g.dx, g.dy, g.dz}, origin[3] = {-g.side / 2, -g.side / 2, g.zmin};
    const int counts[3] = {g.nx, g.ny, g.nz};
    int cell[3];
    for (int d = 0; d < 3; ++d) {
      const double v = std::clamp(((x[std::size_t(d)] + shift[std::size_t(d)]) - origin[d]) / spacing[d], 0., double(counts[d]));
      cell[d] = std::min(int(std::floor(v)), counts[d] - 1);
    }
    return std::size_t(cell[0]) + std::size_t(nx) * (std::size_t(cell[1]) + std::size_t(ny) * std::size_t(cell[2]));
  }
};
// Particle inputs from one part of the census (any partition): plasma node terms
// and the order-independent maxima.
struct ParticlePart {
  std::vector<NodeTerm> terms;
  std::uint64_t carriers = 0, charged = 0, electrons = 0;
  double omega_c_max = 0, v_max = 0, a_max = 0;
};
inline void particle_part(const axial::Mesh& g, const m::PicVector& shift, const CellBounds& cb, std::span<const m::ArenaRateEntry> rows,
                          ParticlePart& out) {
  constexpr double c = 299792458.;
  for (const auto& r : rows) {
    const auto& p = r.native;
    ++out.carriers;
    double u2 = 0;
    for (double x : p.proper_velocity_m_s) u2 += x * x;
    const double v = std::sqrt(u2) / std::sqrt(1 + u2 / (c * c));
    out.v_max = std::max(out.v_max, v);
    if (p.charge_C == 0) continue;
    ++out.charged;
    if (p.charge_C < 0) ++out.electrons;
    plasma_terms(g, shift, p, out.terms);
    const auto cell = cb.cell_of(g, shift, p.position_m);
    const double qm = std::abs(p.charge_C / p.mass_kg);
    out.omega_c_max = std::max(out.omega_c_max, qm * cb.B[cell]);
    out.a_max = std::max(out.a_max, qm * cb.E[cell]);
  }
}
inline void merge_parts(ControlInputs& in, const axial::Mesh& g, std::vector<ParticlePart>& parts) {
  std::vector<NodeTerm> merged;
  const bool single = parts.size() == 1;
  for (auto& p : parts) {
    if (!single) merged.insert(merged.end(), p.terms.begin(), p.terms.end());
    in.carriers += p.carriers; in.charged += p.charged; in.electrons += p.electrons;
    in.omega_c_max = std::max(in.omega_c_max, p.omega_c_max);
    in.v_max = std::max(in.v_max, p.v_max);
    in.a_max = std::max(in.a_max, p.a_max);
  }
  const auto& terms = single ? parts[0].terms : merged;
  in.plasma_terms = terms.size();
  in.omega_p_max = omega_p_max(g, terms);
}

// Accepted-state history of the controller (per run).
struct ControllerState {
  double dt_prev_s = 0;
  int rung_prev = 0, hold = 0;
  double chem_rate = 0, births_rate = 0;
  struct Snap { double t; hc::Vec E, B; };
  std::vector<Snap> ring;   // last <= 4 accepted states, oldest first
  void record_state(double t, const hc::Vec& E, const hc::Vec& B) {
    if (ring.size() == 4) ring.erase(ring.begin());
    ring.push_back({t, E, B});
  }
  // max over E and B of max_i |6 f[t0..t3]_i| / max_i |x_i| (third divided difference).
  bool lte_coefficient(double& out) const {
    if (ring.size() < 4) return false;
    const double t0 = ring[0].t, t1 = ring[1].t, t2 = ring[2].t, t3 = ring[3].t;
    double best = 0;
    for (int which = 0; which < 2; ++which) {
      const auto& a = which ? ring[0].B : ring[0].E;
      const auto& b = which ? ring[1].B : ring[1].E;
      const auto& c = which ? ring[2].B : ring[2].E;
      const auto& d = which ? ring[3].B : ring[3].E;
      double xmax = 0, dmax = 0;
      for (std::size_t i = 0; i < d.size(); ++i) {
        const double f01 = (b[i] - a[i]) / (t1 - t0), f12 = (c[i] - b[i]) / (t2 - t1), f23 = (d[i] - c[i]) / (t3 - t2);
        const double f012 = (f12 - f01) / (t2 - t0), f123 = (f23 - f12) / (t3 - t1);
        dmax = std::max(dmax, std::abs(6 * (f123 - f012) / (t3 - t0)));
        xmax = std::max(xmax, std::abs(d[i]));
      }
      if (xmax > 0) best = std::max(best, dmax / xmax);
    }
    out = best;
    return true;
  }
};

// ------------------------------------------------------------ decision
struct Verdict {
  bool ok = true;
  Criterion failing = Criterion::cap;
  double value = 0, limit = 0;
};
// All step-size criteria at dt. The caps apply to dt; the physical criteria to
// margin * dt (margin > 1 is the growth hysteresis).
inline Verdict admissible(const ControllerConfig& C, const ControlInputs& in, double dt, double margin) {
  const double x = dt * margin;
  constexpr double c = 299792458.;
  const bool particles = in.charged > 0 || (C.neutrals_need_particle_cap && in.carriers > 0);
  auto fail = [](Criterion k, double v, double l) { return Verdict{false, k, v, l}; };
  if (dt > C.cap_s) return fail(Criterion::cap, dt, C.cap_s);
  if (particles && dt > C.particle_cap_s) return fail(Criterion::particle_cap, dt, C.particle_cap_s);
  if (in.omega_p_max * x > C.c_plasma) return fail(Criterion::plasma, in.omega_p_max * x, C.c_plasma);
  if (in.omega_c_max * x > C.c_gyro) return fail(Criterion::gyro, in.omega_c_max * x, C.c_gyro);
  const double disp = x * std::min(c, in.v_max + in.a_max * x) / in.h_min;
  if (disp > C.c_disp) return fail(Criterion::displacement, disp, C.c_disp);
  if (in.chem_rate * (.5 * x) > C.c_chem * in.chem_limit) return fail(Criterion::chemistry, in.chem_rate * (.5 * x), C.c_chem * in.chem_limit);
  if (!particles) {
    if (in.births_rate * x > C.birth_threshold && dt > C.particle_cap_s) return fail(Criterion::births, in.births_rate * x, C.birth_threshold);
    if (!in.lte_valid) {
      if (dt > C.particle_cap_s) return fail(Criterion::lte, dt, C.particle_cap_s);
    } else if (in.lte_coeff * x * x * x / 12 > C.lte_tol)
      return fail(Criterion::lte, in.lte_coeff * x * x * x / 12, C.lte_tol);
  }
  return {};
}
// Guard value and limit of criterion k at dt (for the record of the chosen step).
inline Verdict measure(const ControllerConfig& C, const ControlInputs& in, Criterion k, double dt) {
  constexpr double c = 299792458.;
  switch (k) {
    case Criterion::plasma: return {true, k, in.omega_p_max * dt, C.c_plasma};
    case Criterion::gyro: return {true, k, in.omega_c_max * dt, C.c_gyro};
    case Criterion::displacement: return {true, k, dt * std::min(c, in.v_max + in.a_max * dt) / in.h_min, C.c_disp};
    case Criterion::chemistry: return {true, k, in.chem_rate * (.5 * dt), C.c_chem * in.chem_limit};
    case Criterion::lte: return in.lte_valid ? Verdict{true, k, in.lte_coeff * dt * dt * dt / 12, C.lte_tol} : Verdict{true, k, dt, C.particle_cap_s};
    case Criterion::births: return {true, k, in.births_rate * dt, C.birth_threshold};
    case Criterion::particle_cap: return {true, k, dt, C.particle_cap_s};
    default: return {true, k, dt, C.cap_s};
  }
}
inline void continuous_bounds(const ControllerConfig& C, const ControlInputs& in, Decision& d) {
  if (in.omega_p_max > 0) d.bound_plasma = C.c_plasma / in.omega_p_max;
  if (in.omega_c_max > 0) d.bound_gyro = C.c_gyro / in.omega_c_max;
  const double s = C.c_disp * in.h_min;
  if (in.v_max > 0 || in.a_max > 0) d.bound_disp = 2 * s / (in.v_max + std::sqrt(in.v_max * in.v_max + 4 * in.a_max * s));
  if (in.chem_rate > 0) d.bound_chem = 2 * C.c_chem * in.chem_limit / in.chem_rate;
  if (in.lte_valid && in.lte_coeff > 0 && in.charged == 0 && !(C.neutrals_need_particle_cap && in.carriers > 0))
    d.bound_lte = std::cbrt(12 * C.lte_tol / in.lte_coeff);
}
// Drive criterion needs the model (current at t and t+dt).
template <class Drive>
inline bool drive_ok(const ControllerConfig& C, const ControlInputs& in, const Drive& drive, double dt, double margin, double* rel = nullptr) {
  const double x = dt * margin;
  if (!(in.t_s + x <= drive.end())) return margin > 1 ? false : true;   // the end of the drive is handled as a landing time
  const double r = std::abs(drive.current(in.t_s + x) - in.drive_I) / std::max(std::abs(in.drive_I), C.drive_floor_A);
  if (rel) *rel = r;
  return r <= C.c_drive;
}

template <class Drive>
inline Decision decide(const ControllerConfig& C, const ControlInputs& in, const Drive& drive, bool first_step) {
  Decision d;
  continuous_bounds(C, in, d);
  {   // drive bound (diagnostic): slope over [t, t + cap]
    const double t1 = std::min(in.t_s + C.cap_s, drive.end());
    if (t1 > in.t_s) {
      const double slope = std::abs(drive.current(t1) - in.drive_I) / (t1 - in.t_s);
      if (slope > 0) d.bound_drive = C.c_drive * std::max(std::abs(in.drive_I), C.drive_floor_A) / slope;
    }
  }
  if (!C.adaptive) {
    d.dt_s = C.fixed_dt_s;
    d.criterion = Criterion::fixed;
    d.value = d.dt_s; d.limit = C.fixed_dt_s;
    if (C.fixed_grid_retry) {   // after a substep: complete the grid interval t0 + k dt
      const double k = std::floor((in.t_s + 1e-6 * C.fixed_dt_s - C.sync_origin_s) / C.fixed_dt_s) + 1;
      const double r = (C.sync_origin_s + k * C.fixed_dt_s) - in.t_s;
      if (r < C.fixed_dt_s * (1 - 1e-6)) {
        d.dt_s = r;
        d.landing = true;
        d.criterion = Criterion::sync;
        d.value = r;
      }
    }
    return d;
  }
  const int jtop = top_rung(C);
  const int jmin = first_step ? jtop : std::max(jtop, in.rung_prev - (in.hold > 0 ? 0 : C.growth_rungs));
  int j = jmin;
  Verdict above{};
  bool found = false;
  for (;; ++j) {
    const double dt = ladder_dt(C, j);
    if (dt < C.dt_min_s) break;
    const double margin = (!first_step && j < in.rung_prev) ? C.grow_margin : 1.;
    auto v = admissible(C, in, dt, margin);
    if (v.ok && !drive_ok(C, in, drive, dt, margin)) v = {false, Criterion::drive, 0, C.c_drive};
    if (v.ok) { found = true; break; }
    above = v;
  }
  if (!found) {
    d.refused = true;
    d.criterion = above.failing;
    d.value = above.value; d.limit = above.limit;
    return d;
  }
  d.rung = j;
  d.dt_s = ladder_dt(C, j);
  if (j == jmin) {
    d.criterion = jmin == jtop ? Criterion::cap : (in.hold > 0 ? Criterion::hold : Criterion::growth);
    d.value = d.dt_s; d.limit = d.criterion == Criterion::cap ? C.cap_s : ladder_dt(C, in.rung_prev) * (in.hold > 0 ? 1. : 2.);
  } else {
    d.criterion = above.failing;
    if (above.failing == Criterion::drive) {
      double rel = 0;
      (void)drive_ok(C, in, drive, d.dt_s, 1., &rel);
      d.value = rel; d.limit = C.c_drive;
    } else {
      const auto m = measure(C, in, above.failing, d.dt_s);
      d.value = m.value; d.limit = m.limit;
    }
  }
  // Landing times: the end of the drive and the optional sync grid. A time within
  // 1e-6 of a step of a landing time counts as on it (sums of steps round).
  double land = drive.end();
  if (C.sync_every_s > 0) {
    const double tau = 1e-6 * d.dt_s;
    const double k = std::floor((in.t_s + tau - C.sync_origin_s) / C.sync_every_s) + 1;
    land = std::min(land, C.sync_origin_s + k * C.sync_every_s);
  }
  const double r = land - in.t_s;
  if (r <= d.dt_s) {
    d.dt_s = r;
    d.landing = true;
  } else if (r <= d.dt_s * (1 + 1e-6)) {
    d.landing = true;                  // ends within 1e-6 dt of the landing time
  } else if (C.sync_every_s > 0 && r < 1.5 * d.dt_s) {
    d.dt_s = .5 * r;                   // two equal steps to the landing time (no sliver step)
    d.landing = true;
  }
  if (d.landing && d.dt_s != ladder_dt(C, d.rung)) {
    d.criterion = land == drive.end() ? Criterion::drive_end : Criterion::sync;
    d.value = d.dt_s; d.limit = land;
  }
  return d;
}

// Canonical record of a decision (inputs, configuration digest, decision) and its digest.
inline Json decision_json(const ControllerConfig& C, const ControlInputs& in, const Decision& d) {
  auto hex = [](double x) { std::ostringstream s; s << std::hexfloat << x; return s.str(); };
  auto num = [](double x) { return std::isfinite(x) ? Json(x) : Json(x > 0 ? "inf" : "-inf"); };
  return Json{{"schema", "coordinator-v7-step-decision-v1"}, {"controller_sha256", C.sha256()},
              {"inputs", {{"t_s", hex(in.t_s)}, {"dt_prev_s", hex(in.dt_prev_s)}, {"rung_prev", in.rung_prev}, {"hold", in.hold},
                          {"carriers", in.carriers}, {"charged", in.charged}, {"electrons", in.electrons}, {"omega_p_max", hex(in.omega_p_max)},
                          {"omega_c_max", hex(in.omega_c_max)}, {"v_max", hex(in.v_max)}, {"a_max", hex(in.a_max)}, {"h_min", hex(in.h_min)},
                          {"chem_rate", hex(in.chem_rate)}, {"chem_limit", hex(in.chem_limit)}, {"drive_I", hex(in.drive_I)},
                          {"lte_coeff", hex(in.lte_coeff)}, {"lte_valid", in.lte_valid}, {"births_rate", hex(in.births_rate)},
                          {"plasma_terms", in.plasma_terms}}},
              {"decision", {{"dt_s", hex(d.dt_s)}, {"rung", d.rung}, {"landing", d.landing}, {"refused", d.refused},
                            {"criterion", criterion_name(d.criterion)}, {"value", num(d.value)}, {"limit", num(d.limit)}}},
              {"bounds_s", {{"plasma", num(d.bound_plasma)}, {"gyro", num(d.bound_gyro)}, {"displacement", num(d.bound_disp)},
                            {"chemistry", num(d.bound_chem)}, {"drive", num(d.bound_drive)}, {"lte", num(d.bound_lte)}}}};
}

// Step-size guards of the owners: a refused preparation whose text names one of
// these is retried at a smaller dt (the owners are unchanged by a refused preparation).
inline bool step_size_refusal(const std::string& what) {
  static const char* guards[] = {"plasma guard", "gyro/displacement guard", "displacement guard", "motion budget outside the sealed guards",
                                 "did not converge", "excessive axial path", "coverage limit refused the exposure",
                                 "reaction/event resource limit", "exceeds-capacity", "did not reach the identical authoritative midpoint"};
  for (const char* g : guards)
    if (what.find(g) != std::string::npos) return true;
  return false;
}

// Inputs of the next decision from the accepted owners and the controller history.
// Inputs of the next decision from the accepted owners, the accepted population (ascending
// native ID, as rate_population returns it) and the controller history.
inline ControlInputs control_inputs_from(const OwnersV7& o, std::span<const m::ArenaRateEntry> census, const ControllerState& st,
                                         const ControllerConfig& C, CellBounds& cb, int parts = 1, bool always_bounds = false) {
  ControlInputs in;
  const auto& g = o.hall->mesh();
  in.t_s = o.hall->time_s();
  in.dt_prev_s = st.dt_prev_s;
  in.rung_prev = st.rung_prev;
  in.hold = st.hold;
  in.h_min = std::min({g.dx, g.dy, g.dz});
  // The field bounds enter only through charged carriers (omega_c, a_max).
  if (always_bounds || std::any_of(census.begin(), census.end(), [](const m::ArenaRateEntry& r) { return r.native.charge_C != 0; }))
    cb.build(g, o.hall->field().E, o.hall->field().B);
  const auto shift = native_to_field_shift(o.arena->native_options());
  need(parts >= 1, "census partition count");
  std::vector<ParticlePart> P(static_cast<std::size_t>(parts));
  if (parts == 1) particle_part(g, shift, cb, census, P[0]);
  else {
    // Partition by native ID modulo parts (as ranks would hold the census), parts in reverse order.
    std::vector<std::vector<m::ArenaRateEntry>> rows(static_cast<std::size_t>(parts));
    for (const auto& r : census) rows[std::size_t(r.native.id % std::uint64_t(parts))].push_back(r);
    for (int k = parts - 1; k >= 0; --k) {
      std::reverse(rows[std::size_t(k)].begin(), rows[std::size_t(k)].end());
      particle_part(g, shift, cb, rows[std::size_t(k)], P[std::size_t(parts - 1 - k)]);
    }
  }
  merge_parts(in, g, P);
  in.chem_rate = st.chem_rate;
  in.chem_limit = cv::CoveragePolicy{}.max_exposure_unresolved_events;
  in.drive_I = o.hall->model().drive().current(in.t_s);
  in.drive_end_s = o.hall->model().drive().end();
  in.lte_valid = st.lte_coefficient(in.lte_coeff);
  in.births_rate = st.births_rate;
  (void)C;
  return in;
}
// Coordinator v7 form (its own population traversal, bounds always built).
inline ControlInputs control_inputs(const OwnersV7& o, const ControllerState& st, const ControllerConfig& C, CellBounds& cb, int parts = 1) {
  const auto census = o.arena->rate_population(o.inputs.maximum_carriers);
  return control_inputs_from(o, census, st, C, cb, parts, true);
}
}  // namespace raikou::native_water::v7_1
