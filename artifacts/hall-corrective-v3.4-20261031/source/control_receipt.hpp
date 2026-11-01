#pragma once
// Hall v3.1: per-step quantities for a step-size controller (diagnostic receipt; no identity,
// no effect on any accepted value). The device-resident step fills it on the device every step
// (gpu::ResidentField); host_control_receipt computes the same quantities on the host for a v3
// Owner proposal on request.
//   omega_p_dt_max    max over nodes of sqrt(e * gross_n / (eps0 * m_e * VN_n)) * dt: plasma frequency of
//                     the actual free charge (gross charge of all carriers counted as electrons, so an upper
//                     bound), node dual volumes VN
//   B_max_T           max over cells of |B| at the cell centre (face values of the authoritative midpoint
//                     averaged per axis); omega_c_dt_max = e B_max dt / m_e
//   drive_rel_change  |I(t+dt) - I(t)| / max(|I(t+dt/2)|, 1 A) of the C19 drive; drive_knot_interval_s the
//                     C19 knot spacing at t+dt/2 (the drive is linear between knots)
//   increment_rel     ||E(t+dt) - E(t)||_M / ||E(t+dt)||_M (M = edge permittivity mass)
//   third_difference_rel  ||dE_n - 2 dE_n-1 + dE_n-2||_M / (12 ||E||_M): local truncation estimate of the
//                     trapezoidal (implicit midpoint) update, dt^3 E''' / 12, from the last three accepted
//                     increments (device path; 0 until two accepted increments exist)
// The particle quantities (maximum gyro angle omega_c dt and maximum displacement in cells over the
// carriers) are the motion record budget values maximum_gyro and maximum_displacement_cells.
#include "field_contract.hpp"

namespace hall_corrective_v3 {
struct ControlReceipt {
    double omega_p_dt_max = 0;
    double B_max_T = 0, omega_c_dt_max = 0;
    double drive_current_A = 0, drive_rel_change = 0;
    double drive_knot_interval_s = 0;
    double increment_rel = 0;
    double third_difference_rel = 0;
    int history = 0;
};
inline constexpr double control_electron_charge_C = 1.602176634e-19, control_electron_mass_kg = 9.1093837015e-31;
// Host form (same per-node and per-cell expressions as the device kernels k_omega_p and k_cells).
inline ControlReceipt host_control_receipt(const Mesh& g, const CurrentDrive& drive, const axial::State& old, const axial::State& next,
                                           const FreeCharge& free_next, const Vec& midB, double dt) {
    ControlReceipt c;
    const double factor = control_electron_charge_C / (hall::eps * control_electron_mass_kg);
    for (int n = 0; n < g.nn; ++n) {
        const double v = std::sqrt(factor * free_next.gross_C[n] / g.VN[n]) * dt;
        c.omega_p_dt_max = std::max(c.omega_p_dt_max, std::isfinite(v) ? v : INFINITY);
    }
    for (int k = 0; k < g.nz; ++k) for (int j = 0; j < g.ny; ++j) for (int i = 0; i < g.nx; ++i) {
        const double bx = .5 * (midB[g.bx(i, j, k)] + midB[g.bx(i + 1, j, k)]);
        const double by = .5 * (midB[g.by(i, j, k)] + midB[g.by(i, j + 1, k)]);
        const double bz = .5 * (midB[g.bz(i, j, k)] + midB[g.bz(i, j, k + 1)]);
        const double bm = std::sqrt(bx * bx + by * by + bz * bz);
        c.B_max_T = std::max(c.B_max_T, std::isfinite(bm) ? bm : INFINITY);
    }
    c.omega_c_dt_max = control_electron_charge_C * c.B_max_T / control_electron_mass_kg * dt;
    c.drive_current_A = drive.current(old.time_s + .5 * dt);
    c.drive_rel_change = std::abs(drive.current(old.time_s + dt) - drive.current(old.time_s)) / std::max(std::abs(c.drive_current_A), 1.);
    c.drive_knot_interval_s = drive.knot_interval(old.time_s + .5 * dt);
    long double sd = 0, se = 0;
    for (int e = 0; e < g.ne; ++e) {
        const long double a = (long double)old.E[e] + old.E_low[e], b = (long double)next.E[e] + next.E_low[e];
        sd += g.M[e] * (b - a) * (b - a); se += g.M[e] * b * b;
    }
    c.increment_rel = se > 0 ? double(std::sqrt(sd / se)) : 0.;
    return c;
}
}  // namespace hall_corrective_v3
