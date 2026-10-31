// test_v12.cpp -- v1.2 benchmarks and ledgers (CPU reference): depth-averaged melt layer (dam break against
// Stoker's solution, moving bore against the Rankine-Hugoniot speed, flow along z = flow along x, a flat surface over
// a varying melt depth stays at rest, Poiseuille and Couette limits at low Reynolds number against the lubrication
// model), Lorentz and surface-force work and dissipation, ejected-melt reservoir (open edge, departure limit),
// thermal sub-cycling (exact window Joule, accuracy bound), full v1.2 ledger, checkpoint and layer-flow digest.
//   test_v12 OUT_DIR [suite ...]   suites: stoker bore zflow balance limits work lubwork eject departure subcycle
//                                          ledger checkpoint
#include "check.hpp"
#include "phase_owner.hpp"
#include <chrono>
#include <cstring>
#include <functional>
#include <set>

using namespace ss304_phase_v1_2;
using Clock = std::chrono::steady_clock;

static std::vector<double> uniform_rows(double H, int n) { std::vector<double> e(n + 1); for (int j = 0; j <= n; ++j) e[j] = H * j / n; e.back() = H; return e; }
static Config base(int S, double W, double H, double dz, int nx_em, int ny_em, int rx, std::vector<double> rows) {
    Config c;
    for (int s = 0; s < S; ++s) c.slices.push_back({0, s, -0.5 * W, 0.0});
    c.width = W; c.thickness = H; c.height = dz; c.nx_em = nx_em; c.ny_em = ny_em; c.rx = rx; c.y_edges = std::move(rows);
    return c;
}
static LawOptions liquid_law(double mu, double sigma, double dsdT) {
    LawOptions o; o.kind = 1; o.c_T0 = 300; o.c_cps = o.c_cpl = 500; o.c_ks = o.c_kl = 1e-12; o.c_Tm = 1000; o.c_dTm = 1; o.c_L = 0; o.T_min = 200; o.T_max = 5000;
    o.c_rho = 7000; o.c_mu = mu; o.c_sigma = sigma; o.c_dsdT = dsdT; o.vapor = false;
    return o;
}
static const double kRho = 7000;

// ------------------------------------------------------------------------------------------------ shallow water
// Inviscid layer (mu 1e-9 Pa s, no surface tension) of the whole HighY half column (two uniform rows, all liquid),
// pressed onto the wall by the inward normal Lorentz load (J_z B_x < 0 on HighY): rho g = 7e6 N/m3, g = 1000 m/s2.
struct Sw { int nx; double W, H, g, dz; };
static Config sw_config(const Sw& p, int S, bool along_z) {
    Config c = along_z ? base(S, p.dz, p.H, p.dz, 1, 1, 1, uniform_rows(p.H, 2)) : base(1, p.W, p.H, p.dz, p.nx, 1, 1, uniform_rows(p.H, 2));
    c.law = liquid_law(1e-9, 0.0, 0.0); c.rho_ref = kRho; c.T0 = 1200; c.evap_mode = 0; c.flow = true; c.flow_model = 2;
    c.max_substeps = 100000000;
    return c;
}
static std::vector<double> sw_JB(const Owner& o, double g) {
    std::vector<double> JB(6 * o.dims().n_em, 0.0);
    for (int e = 0; e < o.dims().n_em; ++e) { JB[6 * e + 2] = 1e6; JB[6 * e + 3] = -kRho * g / 1e6; }   // (J x B)_y = J_z B_x = -rho g
    return JB;
}
// layer thickness and velocity of the HighY columns from the state (layer mass of the last sub-step)
static void sw_profile(const Owner& o, std::vector<double>& h, std::vector<double>& u, bool along_z) {
    const auto st = o.state();
    const int nx = o.dims().nx, S = o.dims().S;
    const int n = along_z ? S : nx;
    h.assign(n, 0); u.assign(n, 0);
    for (int k = 0; k < n; ++k) {
        const int fc = along_z ? (k * 2 + 1) * nx : nx + k;
        const double m = st.fcV[3 * fc + 2];
        h[k] = m / (kRho * o.dims().dx * o.dims().dz);
        u[k] = m > 0 ? st.fcV[3 * fc + (along_z ? 1 : 0)] / m : 0.0;
    }
}
static void sw_run(Owner& o, const std::vector<double>& JB, double t_end, double dt, std::int64_t* nsub, std::int64_t* retries) {
    const int steps = int(std::llround(t_end / dt));
    for (int k = 0; k < steps; ++k) {
        StepInput in; in.dt = dt; in.em_power_W.assign(o.dims().n_em, 0.0); in.em_JB = JB;
        const auto r = o.advance(in);
        *nsub += r.n_sub; *retries += r.cfl_retries;
    }
}
// Stoker (1957) wet-bed dam break: middle state from the rarefaction and the shock (Rankine-Hugoniot) relations
struct Stoker { double hL, hR, g, hm, um, s; };
static Stoker stoker(double hL, double hR, double g) {
    auto F = [&](double hm) { const double s = std::sqrt(g * hm * (hm + hR) / (2 * hR)); return 2 * (std::sqrt(g * hL) - std::sqrt(g * hm)) - s * (1 - hR / hm); };
    double a = hR, b = hL;
    for (int k = 0; k < 200; ++k) { const double m = 0.5 * (a + b); if (F(m) > 0) a = m; else b = m; }
    const double hm = 0.5 * (a + b), s = std::sqrt(g * hm * (hm + hR) / (2 * hR));
    return {hL, hR, g, hm, s * (1 - hR / hm), s};
}
static double stoker_h(const Stoker& q, double xi) {   // xi = (x - x0) / t
    const double cL = std::sqrt(q.g * q.hL), cm = std::sqrt(q.g * q.hm);
    if (xi < -cL) return q.hL;
    if (xi < q.um - cm) { const double c = (2 * cL - xi) / 3; return c * c / q.g; }
    if (xi < q.s) return q.hm;
    return q.hR;
}
static double dam_break(int nx, const char* tag, std::vector<double>* hz_out) {
    const Sw p{nx, 0.2, 0.002, 1000.0, 2.5e-4};
    Config c = sw_config(p, 1, false);
    Owner o(c);
    const double h0 = 0.5 * p.H, hL = 1.4e-3, hR = 0.6e-3;
    o.initialize_displacement([&](int, int f, double x) { return f == 1 ? (x < 0 ? hL - h0 : hR - h0) : 0.0; });
    const auto JB = sw_JB(o, p.g);
    std::int64_t nsub = 0, retries = 0;
    const double t = 0.04;
    sw_run(o, JB, t, 1e-3, &nsub, &retries);
    std::vector<double> h, u; sw_profile(o, h, u, false);
    const Stoker q = stoker(hL, hR, p.g);
    const double dx = p.W / nx;
    double l1 = 0, ref = 0, plateau = 0; int np = 0; double xs = 0;
    for (int i = 0; i < nx; ++i) {
        const double x = -0.5 * p.W + (i + 0.5) * dx, he = stoker_h(q, x / t);
        l1 += std::abs(h[i] - he) * dx; ref += std::abs(he - (x < 0 ? hL : hR)) * dx;
        const double cm = std::sqrt(p.g * q.hm);
        if (x / t > q.um - cm + 0.2 * (q.s - q.um + cm) && x / t < q.s - 0.2 * (q.s - q.um + cm)) { plateau += h[i]; ++np; }
        if (i > 0 && h[i - 1] >= 0.5 * (q.hm + hR) && h[i] < 0.5 * (q.hm + hR) && x > 0) xs = x - 0.5 * dx;
    }
    plateau /= np;
    const std::string T = tag;
    chk::put(T + "_nx", nx); chk::put(T + "_substeps", double(nsub)); chk::put(T + "_cfl_retries", double(retries));
    chk::put(T + "_L1_error_over_disturbance", l1 / ref);
    chk::put(T + "_plateau_h_m", plateau); chk::put(T + "_stoker_h_m", q.hm);
    chk::put(T + "_shock_x_m", xs); chk::put(T + "_stoker_shock_x_m", q.s * t);
    chk::check(std::abs(plateau / q.hm - 1) < 0.005, ("dam break " + T + ": middle depth within 0.5 % of Stoker's").c_str());
    chk::check(std::abs(xs - q.s * t) < 3 * dx, ("dam break " + T + ": shock position within 3 columns of the Rankine-Hugoniot shock").c_str());
    chk::check(o.ledger_closes(), ("exact ledgers in the dam break " + T).c_str());
    if (hz_out) *hz_out = h;
    return l1 / ref;
}
static void stoker_suite() {
    const double e1 = dam_break(200, "stoker200", nullptr), e2 = dam_break(400, "stoker400", nullptr);
    chk::put("stoker_L1_ratio_200_to_400", e1 / e2);
    chk::check(e2 < e1 && e2 < 0.05, "dam break: L1 error decreases with resolution and is below 5 % of the disturbance");
}
// moving bore: left state (hL, uL) from the Rankine-Hugoniot relations for a bore of speed s into still water
static void bore() {
    const int nx = 400;
    const Sw p{nx, 0.2, 0.002, 1000.0, 2.5e-4};
    Config c = sw_config(p, 1, false);
    Owner o(c);
    const double h0 = 0.5 * p.H, hL = 1.4e-3, hR = 0.6e-3, xb = -0.05;
    const double s = std::sqrt(p.g * hL * (hL + hR) / (2 * hR)), uL = s * (1 - hR / hL);
    o.initialize_displacement([&](int, int f, double x) { return f == 1 ? (x < xb ? hL - h0 : hR - h0) : 0.0; });
    o.initialize_flow([&](int, int f, double x) { return std::array<double, 2>{f == 1 && x < xb ? uL : 0.0, 0.0}; });
    const auto JB = sw_JB(o, p.g);
    std::int64_t nsub = 0, retries = 0;
    const double t = 0.025;
    sw_run(o, JB, t, 1e-3, &nsub, &retries);
    std::vector<double> h, u; sw_profile(o, h, u, false);
    const double dx = p.W / nx;
    double xs = 0, hl = 0, ul = 0; int nl = 0;
    for (int i = 1; i < nx; ++i) {
        const double x = -0.5 * p.W + (i + 0.5) * dx;
        if (h[i - 1] >= 0.5 * (hL + hR) && h[i] < 0.5 * (hL + hR) && x > xb) xs = x - 0.5 * dx;
        if (x > -0.04 && x < xb + s * t - 0.01) { hl += h[i]; ul += u[i]; ++nl; }
    }
    hl /= nl; ul /= nl;
    chk::put("bore_speed_rankine_hugoniot_m_s", s); chk::put("bore_speed_measured_m_s", (xs - xb) / t);
    chk::put("bore_left_h_m", hl); chk::put("bore_left_u_m_s", ul); chk::put("bore_substeps", double(nsub));
    chk::check(std::abs(xs - (xb + s * t)) < 3 * dx, "moving bore travels at the Rankine-Hugoniot speed (within 3 columns over 0.025 s)");
    chk::check(std::abs(hl / hL - 1) < 0.005 && std::abs(ul / uL - 1) < 0.01, "moving bore: the upstream state is kept (0.5 % depth, 1 % velocity)");
    chk::check(o.ledger_closes(), "exact ledgers with the moving bore (initial kinetic energy booked)");
}
// the same dam break along z (400 slices of one body, one column each) and along x
static void zflow() {
    const int n = 200;
    const Sw p{n, n * 2.5e-4, 0.002, 1000.0, 2.5e-4};
    const double h0 = 0.5 * p.H, hL = 1.4e-3, hR = 0.6e-3, t = 0.01;
    std::vector<double> hx, ux, hz, uz;
    {
        Config c = sw_config(p, 1, false);
        Owner o(c);
        o.initialize_displacement([&](int, int f, double x) { return f == 1 ? (x < 0 ? hL - h0 : hR - h0) : 0.0; });
        std::int64_t a = 0, b = 0; sw_run(o, sw_JB(o, p.g), t, 1e-3, &a, &b);
        sw_profile(o, hx, ux, false);
    }
    {
        Config c = sw_config(p, n, true);
        Owner o(c);
        o.initialize_displacement([&](int s, int f, double) { return f == 1 ? (s < n / 2 ? hL - h0 : hR - h0) : 0.0; });
        std::int64_t a = 0, b = 0; sw_run(o, sw_JB(o, p.g), t, 1e-3, &a, &b);
        sw_profile(o, hz, uz, true);
        chk::check(o.ledger_closes(), "exact ledgers with flow along z");
    }
    double dh = 0, du = 0;
    for (int k = 0; k < n; ++k) { dh = std::max(dh, std::abs(hx[k] - hz[k])); du = std::max(du, std::abs(ux[k] - uz[k])); }
    chk::put("zflow_max_depth_difference_m", dh); chk::put("zflow_max_velocity_difference_m_s", du);
    chk::check(dh <= 1e-12 * hL && du <= 1e-9, "flow along z (through the slices of one body) = flow along x");
}
// flat surface over a melt depth that varies along x (well-balanced hydrostatic reconstruction): no flow
static void balance() {
    const int nx = 40; const double W = 0.02, H = 0.002, g = 1000;
    Config c = base(1, W, H, 2.5e-4, nx, 1, 1, uniform_rows(H, 40));
    c.law = liquid_law(1e-9, 0.0, 0.0); c.rho_ref = kRho; c.T0 = 800; c.evap_mode = 0; c.flow = true; c.flow_model = 2;
    Owner o(c);
    // HighY half: rows 20..39; liquid from the surface down to a depth of 2..18 rows varying with x
    o.initialize_temperature([&](int, double x, double y) { const double depth = 1e-4 * (2 + 16 * (0.5 + 0.5 * std::sin(2 * M_PI * x / W))); return y > H - depth ? 1200.0 : 800.0; });
    const auto JB = sw_JB(o, g);
    double umax = 0;
    for (int k = 0; k < 200; ++k) { StepInput in; in.dt = 1e-4; in.em_power_W.assign(o.dims().n_em, 0.0); in.em_JB = JB; o.advance(in); }
    std::vector<double> h, u; sw_profile(o, h, u, false);
    double hmin = 1, hmax = 0;
    for (int i = 0; i < nx; ++i) { umax = std::max(umax, std::abs(u[i])); hmin = std::min(hmin, h[i]); hmax = std::max(hmax, h[i]); }
    chk::put("balance_max_speed_m_s", umax); chk::put("balance_depth_min_m", hmin); chk::put("balance_depth_max_m", hmax);
    chk::check(hmax > 4 * hmin && umax < 1e-12, "flat surface over a varying melt depth stays at rest (hydrostatic reconstruction; |u| < 1e-12 m/s)");
}

// ------------------------------------------------------------------------------------------------ low-Re limits
// Thin liquid layer (4 rows, 20 um) on a solid substrate (k ~ 0, fixed thickness); f_x from J_z B_y, Marangoni from
// a temperature gradient. Steady depth-averaged velocity at the centre column versus f h^2/(3 mu) and tau h/(2 mu),
// and the depth-averaged flux versus the lubrication flux of the v1.1 model on the same configuration.
static Config thin_layer(int nx, double W, double delta, double H, double sigma, double dsdT, int model) {
    std::vector<double> rows;
    for (int j = 0; j <= 4; ++j) rows.push_back((H - delta) * j / 4);
    for (int j = 1; j <= 4; ++j) rows.push_back(H - delta + delta * j / 4);
    rows.back() = H;
    Config c = base(1, W, H, 1e-4, nx, 1, 1, rows);
    c.law = liquid_law(6e-3, sigma, dsdT); c.rho_ref = kRho; c.evap_mode = 0; c.flow = true; c.T0 = 800; c.flow_model = model;
    return c;
}
static double fc_speed(const Owner& o, int fc) { const auto st = o.state(); return st.fcV[3 * fc + 2] > 0 ? st.fcV[3 * fc] / st.fcV[3 * fc + 2] : 0.0; }
static void limits() {
    const int nx = 20; const double W = 2e-3, H = 1e-4, delta = 2e-5, mu = 6e-3;
    const double f = 1e6;                       // J_z = 1e6 A/m2, B_y = -1 T
    // Poiseuille
    {
        Config c = thin_layer(nx, W, delta, H, 0.0, 0.0, 2);
        Owner o(c);
        o.initialize_temperature([&](int, double, double y) { return y > H - delta ? 1200.0 : 800.0; });
        std::vector<double> JB(6 * o.dims().n_em, 0.0);
        for (int e = 0; e < o.dims().n_em; ++e) { JB[6 * e + 2] = 1e6; JB[6 * e + 4] = -1.0; }
        for (int k = 0; k < 200; ++k) { StepInput in; in.dt = 2e-5; in.em_power_W.assign(o.dims().n_em, 0.0); in.em_JB = JB; o.advance(in); }   // 25 viscous times
        const double u = fc_speed(o, nx + nx / 2), ue = f * delta * delta / (3 * mu), re = kRho * u * delta / mu;
        chk::put("limits_poiseuille_u_m_s", u); chk::put("limits_poiseuille_expected_m_s", ue); chk::put("limits_poiseuille_reynolds", re);
        chk::check(std::abs(u / ue - 1) < 1e-6, "depth-averaged layer, low Re: steady Poiseuille velocity f h^2/(3 mu)");
        chk::check(o.ledger_closes(), "exact ledgers (depth-averaged Poiseuille)");
        // lubrication model on the same configuration: one-step flux at the centre interface
        Config c1 = thin_layer(nx, W, delta, H, 0.0, 0.0, 1);
        Owner o1(c1);
        o1.initialize_temperature([&](int, double, double y) { return y > H - delta ? 1200.0 : 800.0; });
        StepInput in; in.dt = 1e-9; in.em_power_W.assign(o1.dims().n_em, 0.0); in.em_JB = JB; o1.advance(in);
        const double q1 = o1.patch_reports()[1].melt_speed_m_s * delta;   // mean |q| / delta x delta
        chk::put("limits_lubrication_flux_m2_s", q1); chk::put("limits_depth_averaged_flux_m2_s", u * delta);
        chk::check(std::abs(u * delta / q1 - 1) < 1e-6, "depth-averaged flux = lubrication flux (Poiseuille) at low Re");
    }
    // Couette (Marangoni)
    {
        const double dsdT = -4e-4, G = 1e4;   // 10 K/mm
        Config c = thin_layer(nx, W, delta, H, 0.0, dsdT, 2);
        Owner o(c);
        o.initialize_temperature([&](int, double x, double y) { return y > H - delta ? 2000.0 + G * x : 800.0; });
        for (int k = 0; k < 200; ++k) { StepInput in; in.dt = 2e-5; in.em_power_W.assign(o.dims().n_em, 0.0); o.advance(in); }
        const double u = fc_speed(o, nx + nx / 2), ue = dsdT * G * delta / (2 * mu);
        chk::put("limits_couette_u_m_s", u); chk::put("limits_couette_expected_m_s", ue); chk::put("limits_couette_reynolds", kRho * std::abs(u) * delta / mu);
        chk::check(std::abs(u / ue - 1) < 1e-3, "depth-averaged layer, low Re: steady Marangoni (Couette) velocity tau h/(2 mu) (0.1 %, slow advection of the temperature field)");
        Config c1 = thin_layer(nx, W, delta, H, 0.0, dsdT, 1);
        Owner o1(c1);
        o1.initialize_temperature([&](int, double x, double y) { return y > H - delta ? 2000.0 + G * x : 800.0; });
        StepInput in; in.dt = 1e-9; in.em_power_W.assign(o1.dims().n_em, 0.0); o1.advance(in);
        const double q1 = o1.patch_reports()[1].melt_speed_m_s * delta;
        chk::put("limits_couette_lubrication_flux_m2_s", q1);
        chk::check(std::abs(std::abs(u) * delta / q1 - 1) < 1e-3, "depth-averaged flux = lubrication flux (Marangoni) at low Re");
        chk::check(o.ledger_closes(), "exact ledgers (depth-averaged Marangoni)");
    }
}

// ------------------------------------------------------------------------------------------------ work and dissipation
// Steady Poiseuille and Couette layers: the booked Lorentz work equals sum f u h A dt over the columns; the booked
// dissipation equals the viscous dissipation of the Poiseuille-Couette profile, per column 3 mu u^2 / h + tau^2 h/(4 mu)
// (from the column's own depth and velocity; the end columns are not in the uniform state).
static double profile_dissipation(const Owner& o, double tau, double mu, double A, double dt, int i0, int i1) {
    const auto st = o.state(); const int nx = o.dims().nx;
    double D = 0;
    for (int i = i0; i < i1; ++i) {
        const int fc = nx + i; const double m = st.fcV[3 * fc + 2];
        if (m <= 0) continue;
        const double u = st.fcV[3 * fc] / m, h = m / (kRho * A);
        D += (3 * mu * u * u / h + tau * tau * h / (4 * mu)) * A * dt;
    }
    return D;
}
// HighY face of the single slice split into patches [0, 0.1), [0.1, 0.9), [0.9, 1]: the middle patch holds the
// columns 2..17 of 20, away from the end columns (drain and pile-up against the walls).
static int middle_patch(const Owner& o) {
    for (int p = 0; p < o.n_patches(); ++p) { const auto& ps = o.patch_specs()[p]; if (ps.face == 1 && ps.begin == 0.1) return p; }
    return -1;
}
static void work() {
    const int nx = 20; const double W = 2e-3, H = 1e-4, delta = 2e-5, mu = 6e-3, f = 1e6, A = (W / nx) * 1e-4;
    {
        Config c = thin_layer(nx, W, delta, H, 0.0, 0.0, 2);
        c.patches = {{0, 1, 0.0, 0.1}, {0, 1, 0.1, 0.9}, {0, 1, 0.9, 1.0}};
        Owner o(c);
        o.initialize_temperature([&](int, double, double y) { return y > H - delta ? 1200.0 : 800.0; });
        std::vector<double> JB(6 * o.dims().n_em, 0.0);
        for (int e = 0; e < o.dims().n_em; ++e) { JB[6 * e + 2] = 1e6; JB[6 * e + 4] = -1.0; }
        StepRecord r;
        for (int k = 0; k < 200; ++k) { StepInput in; in.dt = 2e-5; in.em_power_W.assign(o.dims().n_em, 0.0); in.em_JB = JB; r = o.advance(in); }
        const auto st = o.state();
        double wl = 0;   // independent: sum f u h A dt over the HighY columns
        for (int i = 0; i < nx; ++i) { const int fc = nx + i; const double m = st.fcV[3 * fc + 2]; if (m > 0) wl += f * (st.fcV[3 * fc] / m) * (m / (kRho * A)) * A * r.dt; }
        const int pm = middle_patch(o);
        const double Qm = r.J(r.patch_dissipation[pm]), Dm = profile_dissipation(o, 0.0, mu, A, r.dt, 2, 18);
        const double Du = f * f * delta * delta * delta / (3 * mu) * A * r.dt * 16;
        chk::put("work_poiseuille_lorentz_J", r.J(r.lorentz)); chk::put("work_poiseuille_independent_J", wl);
        chk::put("work_poiseuille_dissipation_J", r.J(r.dissipation)); chk::put("work_poiseuille_middle_dissipation_J", Qm);
        chk::put("work_poiseuille_middle_profile_dissipation_J", Dm); chk::put("work_poiseuille_middle_uniform_layer_J", Du);
        chk::check(std::abs(r.J(r.lorentz) / wl - 1) < 1e-3, "booked Lorentz work = sum f u h A dt (steady Poiseuille, 0.1 %)");
        chk::check(std::abs(Qm / Dm - 1) < 1e-3, "booked dissipation (interior columns) = viscous dissipation 3 mu u^2/h of the Poiseuille profile of each column (0.1 %)");
        chk::check(r.J(r.dissipation) >= profile_dissipation(o, 0.0, mu, A, r.dt, 0, nx), "total booked dissipation >= friction dissipation (the wall stagnation adds)");
        chk::check(px_is_zero(r.eject_E) && o.ledger_closes(), "exact ledgers with the Lorentz work");
    }
    {
        const double dsdT = -4e-4, G = 1e4, tau = dsdT * G;
        Config c = thin_layer(nx, W, delta, H, 0.0, dsdT, 2);
        c.patches = {{0, 1, 0.0, 0.1}, {0, 1, 0.1, 0.9}, {0, 1, 0.9, 1.0}};
        Owner o(c);
        o.initialize_temperature([&](int, double x, double y) { return y > H - delta ? 2000.0 + G * x : 800.0; });
        StepRecord r;
        for (int k = 0; k < 200; ++k) { StepInput in; in.dt = 2e-5; in.em_power_W.assign(o.dims().n_em, 0.0); r = o.advance(in); }
        const int pm = middle_patch(o);
        const double Qm = r.J(r.patch_dissipation[pm]), Dm = profile_dissipation(o, tau, mu, A, r.dt, 2, 18), Du = tau * tau * delta / mu * A * r.dt * 16;
        chk::put("work_couette_surface_J", r.J(r.surface)); chk::put("work_couette_dissipation_J", r.J(r.dissipation));
        chk::put("work_couette_middle_dissipation_J", Qm); chk::put("work_couette_middle_profile_dissipation_J", Dm); chk::put("work_couette_middle_uniform_layer_J", Du);
        chk::check(std::abs(Qm / Dm - 1) < 5e-3, "booked dissipation (interior columns) of the Marangoni layer = 3 mu u^2/h + tau^2 h/(4 mu) of each column (Couette profile, 0.5 %)");
        chk::check(o.ledger_closes(), "exact ledgers with the surface-force work");
    }
}
// lubrication model with book_work: Lorentz work f q dx dz dt booked as heat, energy rises by exactly W_L + W_S
static void lubwork() {
    const double W = 2e-4, H = 1e-4, delta = 2e-5, mu = 6e-3, f = 1e6, dz = 1e-4, dt = 1e-9;
    Config c = thin_layer(2, W, delta, H, 0.0, 0.0, 1);
    c.book_work = true;
    Owner o(c);
    o.initialize_temperature([&](int, double, double y) { return y > H - delta ? 1200.0 : 800.0; });
    std::vector<double> JB(6 * o.dims().n_em, 0.0);
    for (int e = 0; e < o.dims().n_em; ++e) { JB[6 * e + 2] = 1e6; JB[6 * e + 4] = -1.0; }
    StepInput in; in.dt = dt; in.em_power_W.assign(o.dims().n_em, 0.0); in.em_JB = JB;
    const auto r = o.advance(in);
    const double q = delta * delta * delta * f / (3 * mu), wl = f * q * (W / 2) * dz * dt;
    chk::put("lubwork_lorentz_J", r.J(r.lorentz)); chk::put("lubwork_expected_J", wl); chk::put("lubwork_dE_J", r.J(r.dE));
    chk::check(std::abs(r.J(r.lorentz) / wl - 1) < 1e-9, "lubrication model: booked Lorentz work = f q dx dz dt");
    int ov = 0;
    chk::check(px_eq(r.dE, px_add(r.lorentz, r.surface, &ov)) && !ov && o.ledger_closes(), "lubrication model: the work is dissipated as heat (dE = W_L + W_S exactly)");
}

// ------------------------------------------------------------------------------------------------ ejection
// Layer driven to the open +x edge: the outflow leaves through the ejected-melt reservoir (atoms and energy, thermal
// and kinetic), the ledgers stay exact; ejected mass = outflow rho u h dz dt in the steady state.
static void eject() {
    const int nx = 20; const double W = 2e-3, H = 1e-4, delta = 2e-5, mu = 6e-3, f = 1e6, dz = 1e-4;
    Config c = thin_layer(nx, W, delta, H, 0.0, 0.0, 2);
    c.eject = true; c.open_x = true; c.eject_fraction = 1e9;   // no departure ejection here
    Owner o(c);
    o.initialize_temperature([&](int, double, double y) { return y > H - delta ? 1200.0 : 800.0; });
    std::vector<double> JB(6 * o.dims().n_em, 0.0);
    for (int e = 0; e < o.dims().n_em; ++e) { JB[6 * e + 2] = 1e6; JB[6 * e + 4] = -1.0; }
    StepRecord r;
    double ej_kg = 0;
    for (int k = 0; k < 200; ++k) { StepInput in; in.dt = 2e-5; in.em_power_W.assign(o.dims().n_em, 0.0); in.em_JB = JB; r = o.advance(in); }
    for (int el = 0; el < PX_NEL; ++el) ej_kg += px_to_double(r.eject_N[el], PX_N_INV) * o.law().atom_mass[el];
    const double u = f * delta * delta / (3 * mu), expect = kRho * u * delta * dz * r.dt;
    const auto& L = o.ledger();
    double tot = 0; for (int el = 0; el < PX_NEL; ++el) tot += px_to_double(L.eject_N[el], PX_N_INV) * o.law().atom_mass[el];
    const auto rep = o.patch_reports();
    double rep_tot = 0; for (const auto& p : rep) for (int el = 0; el < PX_NEL; ++el) rep_tot += p.cumulative_ejected_atoms[el] * o.law().atom_mass[el];
    chk::put("eject_last_step_kg", ej_kg); chk::put("eject_expected_step_kg", expect); chk::put("eject_total_kg", tot);
    chk::put("eject_total_J", px_to_double(L.eject_E, PX_E_INV)); chk::put("eject_kinetic_J", px_to_double(L.eject_K, PX_E_INV));
    chk::check(std::abs(ej_kg / expect - 1) < 0.02, "edge outflow into the ejected-melt reservoir = rho u h dz dt (2 %, end-column depth)");
    chk::check(std::abs(rep_tot / tot - 1) < 1e-12, "per-patch ejecta add up to the ledger");
    chk::check(o.ledger_closes(), "exact ledgers with ejection (atoms and energy leave through the reservoir)");
}
// Liquid piled against the +x wall by the Lorentz force rises above the departure limit and is ejected; the
// displacement stays near the limit; ledgers exact.
static void departure() {
    const int nx = 10; const double W = 1e-3, H = 1e-4, delta = 2e-5, dz = 1e-4;
    Config c = thin_layer(nx, W, delta, H, 0.0, 0.0, 2);
    c.eject = true; c.open_x = false; c.eject_fraction = 0.02;   // 2 um of the 100 um EM cell
    Owner o(c);
    o.initialize_temperature([&](int, double, double y) { return y > H - delta ? 1200.0 : 800.0; });
    std::vector<double> JB(6 * o.dims().n_em, 0.0);
    for (int e = 0; e < o.dims().n_em; ++e) { JB[6 * e + 2] = 1e7; JB[6 * e + 4] = -1.0; }
    double emax = 0; StepRecord r;
    for (int k = 0; k < 300; ++k) {
        StepInput in; in.dt = 2e-5; in.em_power_W.assign(o.dims().n_em, 0.0); in.em_JB = JB; r = o.advance(in);
        for (double e : o.face_displacements()) emax = std::max(emax, e / o.em_cell_thickness());   // outward displacement
    }
    const auto& L = o.ledger();
    double tot = 0; for (int el = 0; el < PX_NEL; ++el) tot += px_to_double(L.eject_N[el], PX_N_INV) * o.law().atom_mass[el];
    chk::put("departure_ejected_kg", tot); chk::put("departure_max_outward_displacement_fraction", emax); chk::put("departure_max_crater_fraction", r.max_crater_fraction); chk::put("departure_limit_fraction", 0.02);
    chk::check(tot > 0 && emax < 0.04, "liquid above the departure limit is ejected; the displacement stays near the limit");
    chk::check(o.ledger_closes(), "exact ledgers with departure ejection");
    (void)delta; (void)dz;
}

// ------------------------------------------------------------------------------------------------ sub-cycling
// SS304 v1.2 foil, Joule power varying by field step: owner every field step versus every K field steps with
// accumulated Joule words. Window Joule words equal the sum of the per-step words exactly; temperatures agree within
// the reported window heating estimate; conductivity change over the window reported (the field-solve lag bound).
static Config foil(int S, bool v12) {
    const double W = 0.005, H = 0.0002;
    Config c;
    for (int s = 0; s < S; ++s) c.slices.push_back({s % 2, s / 2, -0.5 * W, s % 2 ? 0.0025 : -0.0027});
    c.width = W; c.thickness = H; c.height = 2e-4; c.nx_em = 5; c.ny_em = 1; c.rx = 2; c.y_edges = graded_rows(H, 0.5e-6, 1.35, 10e-6, 1);
    if (v12) Config::enable_v1_2(c); else Config::enable_v1_1(c);
    return c;
}
static void subcycle() {
    for (int K : {4, 16}) {
        Config c = foil(4, true);
        Owner a(c), b(c);
        auto prof = [&](int s, double x, double) { return 2400.0 + 2.0e5 * (x + 0.0025) + 30.0 * s; };
        a.initialize_temperature(prof); b.initialize_temperature(prof);
        const int nem = a.dims().n_em;
        std::vector<double> JB(6 * nem);
        for (int e = 0; e < nem; ++e) { JB[6 * e + 2] = 5e8; JB[6 * e + 3] = 20.0; JB[6 * e + 4] = -0.05; }
        const double dtf = 2.5e-12;
        px_i128 ja = px_zero(); int ov = 0;
        double lag = 0, heat = 0;
        long windows = 0;
        for (int k = 0; k < 64 * 16 / K * K; ++k) {
            std::vector<double> P(nem); for (int e = 0; e < nem; ++e) P[e] = 3e6 * (1 + 0.3 * std::sin(0.37 * k + e));
            StepInput in; in.dt = dtf; in.em_power_W = P; in.em_JB = JB;
            ja = px_add(ja, a.advance(in).joule, &ov);
            b.accumulate(P, JB, dtf);
            if ((k + 1) % K == 0) {
                StepInput w; w.use_accumulated = true;
                const auto r = b.advance(w);
                lag = std::max(lag, r.sigma_change_max); heat = std::max(heat, r.window_heating_K); ++windows;
            }
        }
        const auto Ta = a.temperatures(), Tb = b.temperatures();
        double dT = 0; for (std::size_t q = 0; q < Ta.size(); ++q) dT = std::max(dT, std::abs(Ta[q] - Tb[q]));
        const std::string T = "subcycle_K" + std::to_string(K);
        chk::put(T + "_windows", double(windows)); chk::put(T + "_max_T_difference_K", dT); chk::put(T + "_window_heating_bound_K", heat);
        chk::put(T + "_sigma_change_max", lag);
        chk::put(T + "_time_difference_s", a.time_s() - b.time_s());
        chk::check(px_eq(ja, a.ledger().joule) && px_eq(a.ledger().joule, b.ledger().joule), (T + ": window Joule words = sum of the per-field-step words, exactly").c_str());
        chk::check(dT <= heat, (T + ": temperatures agree within the window heating estimate").c_str());
        chk::check(a.ledger_closes() && b.ledger_closes(), (T + ": exact ledgers with and without sub-cycling").c_str());
    }
}

// ------------------------------------------------------------------------------------------------ full v1.2 ledger
// SS304 v1.2 law and options on four slices (two bodies, z neighbours), Hall-like J x B, back pressures, loads,
// deposits, sub-cycling K = 4: exact ledgers every step; Reynolds number of the inertial layer.
static void ledger() {
    Config c = foil(4, true);
    Owner o(c);
    const double W = 0.005, H = 0.0002;
    o.initialize_temperature([&](int s, double x, double) { return 2400.0 + 2.0e5 * (x + 0.5 * W) + 30.0 * s; });
    const int nem = o.dims().n_em;
    std::vector<double> JB(6 * nem);
    for (int e = 0; e < nem; ++e) { const double j = 5e8 * (1 + 0.02 * (e % 5)); JB[6 * e + 2] = j; JB[6 * e + 3] = 20.0; JB[6 * e + 4] = -0.05; }
    std::vector<double> pb(o.n_patches(), 0.0);
    for (int p = 0; p < o.n_patches(); ++p) pb[p] = (p % 3) * 2e3;
    const auto t0 = Clock::now();
    bool exact = true; long steps = 0; std::uint64_t seq = 0; double maxre = 0, maxspeed = 0, maxeta = 0, lag = 0;
    const int K = 4;
    for (int k = 0; k < 3000; ++k) {
        for (int q = 0; q < K; ++q) o.accumulate(std::vector<double>(nem, 3e6), JB, 1e-11 / K);
        StepInput in; in.use_accumulated = true; in.back_pressure_Pa = pb;
        for (int s = 0; s < 4; ++s) {
            ++seq; const double y1 = c.slices[s].y0 + H;
            in.loads.push_back({s, -0.5 * W + 0.0011, y1, 2e-9, "t/wall", "wall", seq});
            in.deposits.push_back({s, -0.5 * W + 0.0021, y1, 3, 1e6, 1e-13, "t/mn-return", seq});
        }
        StepRecord r;
        try { r = o.advance(in); } catch (const PhaseRefusal& e) { chk::put_text("ledger12_stop", e.what()); break; }
        ++steps;
        maxre = std::max(maxre, r.max_reynolds); maxeta = std::max(maxeta, r.max_displacement_fraction); lag = std::max(lag, r.sigma_change_max);
        if (k % 100 == 0) { exact = exact && o.ledger_closes(); for (const auto& p : o.patch_reports()) maxspeed = std::max(maxspeed, p.melt_speed_m_s); }
    }
    exact = exact && o.ledger_closes();
    const auto& L = o.ledger();
    chk::check(exact && steps == 3000, "exact energy and atom ledgers with all v1.2 options (inertial layer, ejection, sub-cycling K = 4), 3000 windows");
    chk::put("ledger12_steps", double(steps)); chk::put("ledger12_seconds", std::chrono::duration<double>(Clock::now() - t0).count());
    chk::put("ledger12_max_T", o.max_temperature_K()); chk::put("ledger12_max_reynolds", maxre); chk::put("ledger12_max_speed_m_s", maxspeed);
    chk::put("ledger12_max_displacement_fraction", maxeta); chk::put("ledger12_sigma_change_max", lag);
    chk::put("ledger12_lorentz_work_J", px_to_double(L.lorentz, PX_E_INV)); chk::put("ledger12_surface_work_J", px_to_double(L.surface, PX_E_INV));
    chk::put("ledger12_dissipation_J", px_to_double(L.dissipation, PX_E_INV)); chk::put("ledger12_kinetic_J", px_to_double(o.kinetic_word(), PX_E_INV));
    chk::put("ledger12_ejected_J", px_to_double(L.eject_E, PX_E_INV)); chk::put("ledger12_joule_J", px_to_double(L.joule, PX_E_INV));
    double ej = 0; for (int el = 0; el < PX_NEL; ++el) ej += px_to_double(L.eject_N[el], PX_N_INV) * o.law().atom_mass[el];
    chk::put("ledger12_ejected_kg", ej);
    chk::put_text("ledger12_state_sha256", o.state_sha256());
    chk::check(o.verify_state_leaves(), "Merkle leaves of the state with the layer-flow section: engine = host recomputation");
}
// checkpoint with a moving layer: save, restore into a fresh owner, continue both: identical digests
static void checkpoint(const std::string& out) {
    Config c = foil(2, true);
    Owner a(c);
    a.initialize_temperature([&](int s, double x, double) { return 2600.0 + 2.0e5 * (x + 0.0025) + 30.0 * s; });
    const int nem = a.dims().n_em;
    std::vector<double> JB(6 * nem);
    for (int e = 0; e < nem; ++e) { JB[6 * e + 2] = 5e8; JB[6 * e + 3] = 20.0; JB[6 * e + 4] = -0.05; }
    auto step = [&](Owner& o) { StepInput in; in.dt = 1e-11; in.em_power_W.assign(nem, 3e6); in.em_JB = JB; o.advance(in); };
    for (int k = 0; k < 200; ++k) step(a);
    const std::string path = out + "/v12-checkpoint.bin";
    std::filesystem::remove(path);
    a.save(path);
    Owner b(c);
    b.load(path);
    bool same = a.state_sha256() == b.state_sha256();
    for (int k = 0; k < 100; ++k) { step(a); step(b); same = same && a.state_sha256() == b.state_sha256(); }
    const auto st = a.state(); double kin = 0; for (auto w : st.fcK) kin += px_to_double(w, PX_E_INV);
    chk::put("checkpoint_layer_kinetic_J", kin);
    chk::check(kin > 0 && same, "checkpoint with a moving layer: restored owner continues bit for bit");
    std::filesystem::remove(path);
}

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : ".";
    std::set<std::string> want; for (int i = 2; i < argc; ++i) want.insert(argv[i]);
    auto timed = [&](const char* s, const std::function<void()>& f) {
        if (!want.empty() && !want.count(s)) return;
        const auto t = Clock::now(); f();
        chk::put(std::string("seconds_") + s, std::chrono::duration<double>(Clock::now() - t).count());
        std::printf("suite %s done\n", s); std::fflush(stdout);
    };
    timed("stoker", stoker_suite);
    timed("bore", bore);
    timed("zflow", zflow);
    timed("balance", balance);
    timed("limits", limits);
    timed("work", work);
    timed("lubwork", lubwork);
    timed("eject", eject);
    timed("departure", departure);
    timed("subcycle", subcycle);
    timed("checkpoint", [&] { checkpoint(out); });
    timed("ledger", ledger);
    std::string name = "test_v12"; for (const auto& w : want) name += "-" + w;
    return chk::finish(name, out);
}
