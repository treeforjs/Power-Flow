#pragma once
// c19-circuit-drive-v1 (2026-10-31): circuit drive for the Hall port of the Power-Flow coupled model.
//
// "C19" names the drive trace of Mykonos shot 17133 recovered from dissertation Figure C.19; the driver is the
// Mykonos linear transformer driver (LTD) at Sandia. This header supplies:
//   PrescribedCurrentDrive  the regression mode: the C19 trace with the exact parse and interpolation arithmetic of
//                           Hall corrective v3/v3.1 CurrentDrive (bit for bit), plus a source-energy ledger;
//   LumpedLtdCircuit        the lumped Mykonos equivalent (supplied machine data where available):
//                             C (144 nF erected, charged to 800 kV) - R_s(t) - L_ltd (36 nH) - node A (dI/dt
//                             monitor / water-line output) - [core branch L_m || R_c to return] - L_feed (40 nH)
//                             - R_feed - Hall port;
//                           implicit-midpoint step (exact discrete energy identity), per-step port relation
//                           I_mid = a - Y V_mid, and an energy ledger (stored, dissipated, delivered);
//   PortCoupler             the per-step circuit/port iteration (Newton with a secant port-impedance estimate) for a
//                           port evaluated by a field step at a trial current (Hall prepare), joined ledger terms.
// CPU reference, a few scalars per step, deterministic (fixed operation order, no threads, no reductions over
// carriers). Nothing here is called per carrier.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace c19_circuit_v1 {

inline void require(bool ok, const char* what) { if (!ok) throw std::runtime_error(what); }
inline constexpr double eps = std::numeric_limits<double>::epsilon();

// Error-free accumulation (two-sum), fixed order. Used for every ledger total.
struct Sum2 {
    double hi = 0, lo = 0;
    void add(double x) { const double s = hi + x, v = s - hi; lo += (hi - (s - v)) + (x - v); hi = s; }
    double value() const { return hi + lo; }
};

// ---------------------------------------------------------------------------------------------------------------
// Regression mode. Same bytes, same parse (std::stod of CSV columns 4 = time_s and 8 = current_A, header skipped),
// same strictly increasing support check and the same interpolation expression as
// hall_corrective_v3::CurrentDrive (runtimes/hall-corrective-v3.1-20261031/source/field_contract.hpp). The sha256
// check of the bytes stays in the Hall field model; the regression test binds this class to the Hall class directly.
class PrescribedCurrentDrive {
    std::vector<std::pair<double, double>> knots_;
public:
    explicit PrescribedCurrentDrive(const std::filesystem::path& path) {
        std::ifstream in(path, std::ios::binary);
        require(bool(in), "missing C19 current");
        std::string bytes((std::istreambuf_iterator<char>(in)), {});
        std::istringstream all(bytes); std::string line; std::getline(all, line);
        while (std::getline(all, line)) {
            std::istringstream row(line); std::string cell; std::vector<std::string> v;
            while (std::getline(row, cell, ',')) v.push_back(cell);
            require(v.size() == 12, "C19 row shape");
            double t = std::stod(v[4]), I = std::stod(v[8]);
            require(std::isfinite(t) && std::isfinite(I) && (knots_.empty() || t > knots_.back().first), "C19 finite monotone support");
            knots_.emplace_back(t, I);
        }
        require(knots_.size() > 1, "C19 support empty");
    }
    double start() const { return knots_.front().first; }
    double end() const { return knots_.back().first; }
    const std::vector<std::pair<double, double>>& knots() const { return knots_; }
    double current(double t) const {
        require(std::isfinite(t) && t >= start() && t <= end(), "C19 query outside original support");
        if (t == end()) return knots_.back().second;
        auto high = std::upper_bound(knots_.begin(), knots_.end(), t, [](double x, const auto& p) { return x < p.first; });
        auto low = high - 1;
        return low->second + (t - low->first) / (high->first - low->first) * (high->second - low->second);
    }
    // The Hall step queries the drive at old.time_s + .5*dt; the expression is reproduced literally.
    double step_current(double old_time_s, double dt) const { return current(old_time_s + .5 * dt); }
};

// Ledger of the prescribed (ideal current source) mode: the source supplies exactly what the field step receives.
struct PrescribedLedger {
    std::uint64_t steps = 0;
    Sum2 source_supplied_J;   // = sum of Hall circuit_work_J
    Sum2 charge_C;            // sum dt * I
    void accept(double dt, double current, double hall_circuit_work_J) {
        source_supplied_J.add(hall_circuit_work_J); charge_C.add(dt * current); ++steps;
    }
};

// ---------------------------------------------------------------------------------------------------------------
// Lumped LTD equivalent.
struct Parameters {
    // Supplied machine data (confirmed 2026-10-31): 5 stages x 36 bricks x 2 x 40 nF, +/-80 kV.
    double C_F = 144e-9;            // erected capacitance: 20 nF/brick, 720 nF/stage, 5 stages in series
    double V0_V = 800e3;            // open-circuit erected voltage (160 kV per stage)
    double L_ltd_H = 36e-9;         // total LTD inductance, all 5 stages (supplied)
    double L_feed_H = 40e-9;        // dI/dt monitor to target entrance (supplied; 36 nH alternative recalled)
    // Model-derived (fitted to shots 17133 and 15210; see calibration/).
    double R_inf_ohm = 0.0;         // late-time series resistance (bricks, switches, capacitor ESR, line)
    double R_sw_ohm = 0.0;          // switch resistive-phase amplitude: R_s(t) = R_inf + R_sw exp(-(t-t_trigger)/tau_sw)
    double tau_sw_s = 1e-9;
    double inv_L_core_per_H = 0.0;  // 1/L_m of the core magnetizing branch at node A (0 = absent)
    double G_core_S = 0.0;          // 1/R_c of the core loss branch at node A (0 = absent)
    double R_feed_ohm = 0.0;
    double t_trigger_s = 0.0;       // switch closure on the circuit clock; before it no current flows
    void check() const {
        for (double v : {C_F, V0_V, L_ltd_H, L_feed_H, R_inf_ohm, R_sw_ohm, tau_sw_s, inv_L_core_per_H, G_core_S, R_feed_ohm, t_trigger_s})
            require(std::isfinite(v), "nonfinite circuit parameter");
        require(C_F > 0 && L_ltd_H > 0 && L_feed_H > 0 && tau_sw_s > 0 && R_inf_ohm >= 0 && R_sw_ohm >= 0 &&
                inv_L_core_per_H >= 0 && G_core_S >= 0 && R_feed_ohm >= 0, "circuit parameter outside its domain");
    }
    double series_R(double t) const {
        if (t < t_trigger_s) return 0;
        return R_inf_ohm + R_sw_ohm * std::exp(-(t - t_trigger_s) / tau_sw_s);
    }
};

struct State {
    double time_s = 0;
    double v_C = 0;      // capacitor (erected) voltage
    double i_s = 0;      // series (LTD) inductor current
    double i_m = 0;      // core magnetizing current
    double i_f = 0;      // feed current = port current = dI/dt-monitor current
    std::uint64_t revision = 0;
};

// Midpoint solution of one step.
struct Midpoint {
    double v_C, i_s, i_m, i_f, v_A, V_port, R_s;
};

struct StepLedger {
    double dE_stored_J = 0;        // increment form: C dv v_mid + L_s di_s i_s,mid + L_m di_m i_m,mid + L_f di_f i_f,mid
    double dissipated_J = 0;       // dt (R_s i_s^2 + G_c v_A^2 + R_f i_f^2) at the midpoint
    double delivered_J = 0;        // dt i_f,mid V_port (circuit side of the port)
    double residual_J = 0;         // dE + dissipated + delivered (exact arithmetic: 0)
    double rounding_bound_J = 0;   // 16 eps * sum of |terms|
};

struct Ledger {
    std::uint64_t steps = 0;
    double initial_stored_J = 0;
    Sum2 dissipated_J, delivered_J, residual_J, dE_J;
    double max_step_residual_over_bound = 0;
    // Joined to the Hall step:
    Sum2 hall_circuit_work_J;           // sum of Hall circuit_work_J (authoritative energy received by the field)
    Sum2 port_iteration_defect_J;       // sum dt I (V_circuit - V_T): circuit/port voltage mismatch after the iteration
    Sum2 hall_dot_rounding_J;           // sum (dt I V_T - circuit_work_J): Hall's per-edge sum vs terminal-voltage product
};

class LumpedLtdCircuit {
    Parameters p_;
    State s_;
    Ledger ledger_;
public:
    explicit LumpedLtdCircuit(Parameters p, double t0 = 0) : p_(p) {
        p_.check();
        s_.time_s = t0;
        // Before switch closure the erected capacitor holds V0 and no current flows.
        s_.v_C = p_.V0_V;
        ledger_.initial_stored_J = stored_energy(s_);
    }
    // Start from an explicit state (e.g. a checkpoint or a calibration-derived state at the coupled genesis).
    LumpedLtdCircuit(Parameters p, const State& s) : p_(p), s_(s) { p_.check(); ledger_.initial_stored_J = stored_energy(s_); }

    const Parameters& parameters() const { return p_; }
    const State& state() const { return s_; }
    const Ledger& ledger() const { return ledger_; }

    double stored_energy(const State& s) const {
        double E = .5 * p_.C_F * s.v_C * s.v_C + .5 * p_.L_ltd_H * s.i_s * s.i_s + .5 * p_.L_feed_H * s.i_f * s.i_f;
        if (p_.inv_L_core_per_H > 0) E += .5 * s.i_m * s.i_m / p_.inv_L_core_per_H;
        return E;
    }

    // Solves the midpoint system of one step for a given port voltage (port_current == nullptr) or a given port
    // current (port_voltage ignored). Unknowns u = (v, i_s, i_m, i_f, v_A) at the midpoint; with the current given,
    // i_f is replaced by V_port. Gaussian elimination with partial pivoting in 80-bit long double (deterministic on
    // x86-64), rounded once to double, so the midpoint values carry ~1 ulp error and the energy identity closes at the
    // rounding level of its terms.
    Midpoint solve(double dt, double port_voltage, const double* port_current) const {
        require(std::isfinite(dt) && dt > 0, "circuit dt");
        const double tm = s_.time_s + .5 * dt;
        const bool closed = tm >= p_.t_trigger_s;
        const double Rs = p_.series_R(tm);
        const double C = p_.C_F, Ls = p_.L_ltd_H, Lf = p_.L_feed_H, Gm = p_.inv_L_core_per_H, Gc = p_.G_core_S, Rf = p_.R_feed_ohm;
        long double A[5][6] = {};
        // row 0: (2C/dt) v + i_s = (2C/dt) v0
        A[0][0] = 2.L * C / dt; A[0][1] = 1; A[0][5] = 2.L * C / dt * s_.v_C;
        // row 1: -v + (2Ls/dt + Rs) i_s + v_A = (2Ls/dt) i_s0      (open switch: i_s = 0)
        if (closed) { A[1][0] = -1; A[1][1] = 2.L * Ls / dt + Rs; A[1][4] = 1; A[1][5] = 2.L * Ls / dt * s_.i_s; }
        else { A[1][1] = 1; A[1][5] = 0; }
        // row 2: (2/dt) i_m - Gm v_A = (2/dt) i_m0
        A[2][2] = 2.L / dt; A[2][4] = -Gm; A[2][5] = 2.L / dt * s_.i_m;
        // row 3: (2Lf/dt + Rf) i_f - v_A + V_port = (2Lf/dt) i_f0
        // row 4: i_s - i_m - i_f - Gc v_A = 0
        A[4][1] = 1; A[4][2] = -1; A[4][4] = -Gc;
        if (!port_current) {
            A[3][3] = 2.L * Lf / dt + Rf; A[3][4] = -1; A[3][5] = 2.L * Lf / dt * s_.i_f - port_voltage;
            A[4][3] = -1;
        } else {
            // column 3 now holds V_port; i_f = I is known.
            const double I = *port_current;
            A[3][3] = 1; A[3][4] = -1; A[3][5] = 2.L * Lf / dt * s_.i_f - (2.L * Lf / dt + Rf) * I;
            A[4][5] = I;
        }
        for (int c = 0; c < 5; ++c) {
            int piv = c;
            for (int r = c + 1; r < 5; ++r) if (std::abs(A[r][c]) > std::abs(A[piv][c])) piv = r;
            require(A[piv][c] != 0, "singular circuit midpoint system");
            if (piv != c) for (int k = 0; k < 6; ++k) std::swap(A[c][k], A[piv][k]);
            for (int r = c + 1; r < 5; ++r) {
                const long double f = A[r][c] / A[c][c];
                if (f != 0) for (int k = c; k < 6; ++k) A[r][k] -= f * A[c][k];
            }
        }
        long double u[5];
        for (int r = 4; r >= 0; --r) {
            long double s = A[r][5];
            for (int k = r + 1; k < 5; ++k) s -= A[r][k] * u[k];
            u[r] = s / A[r][r];
        }
        Midpoint m;
        m.v_C = double(u[0]); m.i_s = double(u[1]); m.i_m = double(u[2]); m.v_A = double(u[4]); m.R_s = Rs;
        if (!port_current) { m.i_f = double(u[3]); m.V_port = port_voltage; }
        else { m.i_f = *port_current; m.V_port = double(u[3]); }
        return m;
    }

    // Affine port relation of the step: I_mid = a - Y V_mid (exact up to rounding: the system is linear).
    std::pair<double, double> port_relation(double dt) const {
        const auto m0 = solve(dt, 0., nullptr), m1 = solve(dt, 1., nullptr);
        return {m0.i_f, m0.i_f - m1.i_f};
    }

    // Accept a step with the field step's terminal voltage V_T as the port voltage (the circuit solves its own midpoint
    // current i_c = a - Y V_T), given the impressed current I_port the field step used and its circuit_work_J. The
    // circuit ledger closes with its own delivered energy dt i_c V_T; the join terms are recorded:
    //   port_iteration_defect = dt (i_c - I_port) V_T   (vanishes with the iteration residual F = I_port - i_c)
    //   hall_dot_rounding     = dt I_port V_T - circuit_work_J (Hall's per-edge sum against the terminal-voltage product)
    // The voltage is pinned (not the current) because the circuit's step impedance 2 L / dt is large at dt = 10 fs:
    // a current pinned to rounding still leaves a large voltage error, while a pinned voltage leaves a current error of
    // order Y |delta V| with Y ~ dt / (2 L), i.e. at rounding level.
    StepLedger accept(double dt, double I_port, double V_T, double hall_circuit_work_J) {
        const auto m = solve(dt, V_T, nullptr);
        State n = s_;
        n.v_C = 2 * m.v_C - s_.v_C; n.i_s = 2 * m.i_s - s_.i_s; n.i_m = 2 * m.i_m - s_.i_m; n.i_f = 2 * m.i_f - s_.i_f;
        n.time_s = s_.time_s + dt; ++n.revision;
        for (double v : {n.v_C, n.i_s, n.i_m, n.i_f}) require(std::isfinite(v), "nonfinite circuit state");
        StepLedger L;
        // Increment form of the stored-energy change (uses the actual new state, so it includes the midpoint
        // reconstruction rounding).
        const double dv = n.v_C - s_.v_C, dis = n.i_s - s_.i_s, dim = n.i_m - s_.i_m, dif = n.i_f - s_.i_f;
        const double vm = .5 * (n.v_C + s_.v_C), ism = .5 * (n.i_s + s_.i_s), imm = .5 * (n.i_m + s_.i_m), ifm = .5 * (n.i_f + s_.i_f);
        const double tC = p_.C_F * dv * vm, tS = p_.L_ltd_H * dis * ism, tF = p_.L_feed_H * dif * ifm;
        const double tM = p_.inv_L_core_per_H > 0 ? dim * imm / p_.inv_L_core_per_H : 0.;
        const double dRs = dt * m.R_s * m.i_s * m.i_s, dGc = dt * p_.G_core_S * m.v_A * m.v_A, dRf = dt * p_.R_feed_ohm * m.i_f * m.i_f;
        const double pw = dt * m.i_f * m.V_port;
        Sum2 de; de.add(tC); de.add(tS); de.add(tM); de.add(tF);
        Sum2 di; di.add(dRs); di.add(dGc); di.add(dRf);
        Sum2 r; r.add(tC); r.add(tS); r.add(tM); r.add(tF); r.add(dRs); r.add(dGc); r.add(dRf); r.add(pw);
        L.dE_stored_J = de.value(); L.dissipated_J = di.value(); L.delivered_J = pw; L.residual_J = r.value();
        // Rounding bound: arithmetic of the terms plus the representation of the stored state (each state value is a
        // double, so its energy is known to ~eps |x| |x_mid| Q even when the per-step change is much smaller).
        const double Lm = p_.inv_L_core_per_H > 0 ? 1 / p_.inv_L_core_per_H : 0.;
        const double state_rep = p_.C_F * (std::abs(n.v_C) + std::abs(s_.v_C)) * std::abs(vm) + p_.L_ltd_H * (std::abs(n.i_s) + std::abs(s_.i_s)) * std::abs(ism) +
                                 Lm * (std::abs(n.i_m) + std::abs(s_.i_m)) * std::abs(imm) + p_.L_feed_H * (std::abs(n.i_f) + std::abs(s_.i_f)) * std::abs(ifm);
        L.rounding_bound_J = 16 * eps * (std::abs(tC) + std::abs(tS) + std::abs(tM) + std::abs(tF) + std::abs(dRs) + std::abs(dGc) +
                                         std::abs(dRf) + std::abs(pw)) + 4 * eps * state_rep + std::numeric_limits<double>::denorm_min();
        ledger_.dE_J.add(L.dE_stored_J); ledger_.dissipated_J.add(L.dissipated_J); ledger_.delivered_J.add(pw);
        ledger_.residual_J.add(L.residual_J);
        ledger_.max_step_residual_over_bound = std::max(ledger_.max_step_residual_over_bound, std::abs(L.residual_J) / L.rounding_bound_J);
        ledger_.hall_circuit_work_J.add(hall_circuit_work_J);
        ledger_.port_iteration_defect_J.add(dt * (m.i_f - I_port) * V_T);
        ledger_.hall_dot_rounding_J.add(dt * I_port * V_T - hall_circuit_work_J);
        ++ledger_.steps;
        s_ = n;
        return L;
    }
    // Circuit-only advance against a given port voltage (fixtures).
    Midpoint advance_with_voltage(double dt, double V_port) {
        const auto m = solve(dt, V_port, nullptr);
        accept(dt, m.i_f, V_port, dt * m.i_f * V_port);
        return m;
    }
};

// ---------------------------------------------------------------------------------------------------------------
// Per-step circuit/port coupling. `Port` provides evaluate(I) -> terminal voltage of the field step prepared at the
// trial impressed current I (it keeps the prepared candidate of the last evaluation). The iteration solves
// F(I) = I - (a - Y V_T(I)) = 0 by Newton steps with a port-impedance estimate Z (secant within the step; carried
// across steps). For a field step that is affine in I at fixed plasma current and conductivity (Hall), the secant
// is exact up to rounding, so the iteration ends after two or three evaluations.
struct CouplerOptions {
    int max_evaluations = 8;
    double tol_rel_I = 4e-14;      // |F| <= tol_rel_I |I| + tol_abs_I
    double tol_abs_I = 1e-9;       // A
    double initial_Z_ohm = 0.0;    // port impedance estimate for the first step (0 = plain fixed point first)
};

struct CouplerStep {
    int evaluations = 0;
    double I = 0, V_T = 0, F = 0, a = 0, Y = 0, Z = 0;
    bool converged = false;
};

class PortCoupler {
    CouplerOptions o_;
    double Z_ = 0, last_V_ = 0;
    bool have_V_ = false;
public:
    explicit PortCoupler(CouplerOptions o = {}) : o_(o), Z_(o.initial_Z_ohm) {}
    double impedance_estimate() const { return Z_; }
    // evaluate(I) must run the field step at trial current I and return its terminal voltage.
    // others_converged (optional) reports whether the other unknowns of a joint fixed point (plasma current, motion)
    // have converged at the last evaluation; the step then ends only when both are converged. While they are not,
    // the secant is not used to update Z (the voltage change would mix in their change).
    CouplerStep step(const LumpedLtdCircuit& c, double dt, const std::function<double(double)>& evaluate,
                     const std::function<bool()>& others_converged = {}) {
        CouplerStep r;
        const auto [a, Y] = c.port_relation(dt);
        r.a = a; r.Y = Y;
        auto others = [&] { return !others_converged || others_converged(); };
        auto ok = [&](double I, double F) { return std::abs(F) <= o_.tol_rel_I * std::abs(I) + o_.tol_abs_I; };
        double I = have_V_ ? a - Y * last_V_ : a;
        double V = evaluate(I); ++r.evaluations;
        double F = I - (a - Y * V);
        bool prev_clean = others();
        while (!(ok(I, F) && others()) && r.evaluations < o_.max_evaluations) {
            const double Ip = I, Vp = V;
            I = I - F / (1 + Y * Z_);
            V = evaluate(I); ++r.evaluations;
            F = I - (a - Y * V);
            const bool clean = others();
            // Secant update only from a current change well above rounding (a secant over a few ulps is noise).
            if (prev_clean && clean && std::abs(I - Ip) > 1e-11 * std::max(1., std::abs(I))) {
                const double Zs = (V - Vp) / (I - Ip); if (std::isfinite(Zs) && 1 + Y * Zs > 0) Z_ = Zs; }
            prev_clean = clean;
        }
        r.I = I; r.V_T = V; r.F = F; r.Z = Z_;
        r.converged = ok(I, F) && others();
        last_V_ = V; have_V_ = true;
        return r;
    }
};

}  // namespace c19_circuit_v1
