// test_v11.cpp -- v1.1 owner benchmarks and ledgers (CPU reference).
//   test_v11 OUT_DIR [suite ...]   suites: pulse knudsen poiseuille marangoni leveling ledger geometry explosion
//                                          radiation expansion merkle
#include "check.hpp"
#include "phase_owner.hpp"
#include <chrono>
#include <functional>
#include <set>

using namespace ss304_phase_v1_1;
using Clock = std::chrono::steady_clock;

static std::vector<double> uniform_rows(double H, int n) { std::vector<double> e(n + 1); for (int j = 0; j <= n; ++j) e[j] = H * j / n; e.back() = H; return e; }
static Config base(int S, double W, double H, double dz, int nx_em, int ny_em, int rx, std::vector<double> rows) {
    Config c;
    for (int s = 0; s < S; ++s) c.slices.push_back({0, s, -0.5 * W, 0.0});
    c.width = W; c.thickness = H; c.height = dz; c.nx_em = nx_em; c.ny_em = ny_em; c.rx = rx; c.y_edges = std::move(rows);
    return c;
}
static std::string jarr(const std::vector<double>& v) { std::string s = "["; for (std::size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + chk::num(v[i]); return s + "]"; }

// ------------------------------------------------------------------------------------------------ pulse heating
// Ohmic pulse heating emulated on the owner (uniform Joule power, no evaporation, adiabatic): owner enthalpy
// (admitted Joule / mass, relative to 298.15 K) and resistivity (1/sigma per EM cell) at the owner temperature,
// versus Pichler et al. 2020 Table 6, for the v1.1 and the v1 laws.
static void pulse() {
    struct Row { double T, H, rho; } rows[] = {{1800, 1190, 1.495}, {2000, 1360, 1.532}, {2200, 1529, 1.570}, {2400, 1699, 1.608}, {2600, 1868, 1.646}, {2800, 2037, 1.684}};
    for (int mode = 0; mode < 2; ++mode) {
        Config c = base(1, 0.002, 0.0002, 2e-4, 2, 1, 1, uniform_rows(0.0002, 4));
        if (mode) Config::enable_v1_1(c);
        c.evap_mode = 0; c.law.radiation = 0; c.flow = false;
        Owner o(c);
        double mass = 0; { const auto N = o.atom_words(); for (int el = 0; el < PQ_NEL; ++el) mass += pq_to_double(N[el], PQ_N_INV) * o.law().atom_mass[el]; }
        const double p = 1e17, V = 0.001 * 0.0002 * 2e-4, dt = 2e-11;   // 1e17 W/m3: about 1e10 K/s, the pulse-heating regime
        double maxH = 0, maxR = 0; std::size_t k = 0;
        std::vector<double> hs, rs;
        while (k < 6) {
            StepInput in; in.dt = dt; in.em_power_W = {p * V, p * V};
            o.advance(in);
            const double T = o.max_temperature_K();
            if (T >= rows[k].T) {
                const double H = (pq_to_double(o.ledger().joule, PQ_E_INV) / mass + pq_h(&o.law(), 300.0) - pq_ss_h_solid(298.15)) / 1000.0;
                // data at the owner temperature (Pichler fits, the source of Table 6)
                const double Hd = -335 + 0.847 * T, Rd = 1.154 + 1.893e-4 * T;
                const double R = 1e6 / o.sigma_em()[0];
                maxH = std::max(maxH, std::abs(H / Hd - 1)); maxR = std::max(maxR, std::abs(R / Rd - 1));
                hs.push_back(H); rs.push_back(R);
                ++k;
            }
        }
        chk::put_raw(std::string("pulse_owner_H_kJ_kg_") + (mode ? "v11" : "v1"), jarr(hs));
        chk::put_raw(std::string("pulse_owner_rho_uOhm_m_") + (mode ? "v11" : "v1"), jarr(rs));
        chk::put(std::string("pulse_owner_H_maxdev_") + (mode ? "v11" : "v1"), maxH);
        chk::put(std::string("pulse_owner_rho_maxdev_") + (mode ? "v11" : "v1"), maxR);
        if (mode) chk::check(maxH < 0.029 * 0.1 && maxR < 0.045 * 0.1, "owner reproduces the pulse-heating enthalpy and resistivity (within a tenth of their uncertainty)");
        chk::check(o.ledger_closes(), "exact ledgers (pulse heating)");
    }
}

// ------------------------------------------------------------------------------------------------ Knudsen layer
static void knudsen() {
    std::vector<double> evap;
    double ps = 0;
    for (double frac : {0.0, 0.6, 1.5}) {
        Config c = base(1, 0.001, 0.0002, 2e-4, 1, 1, 1, uniform_rows(0.0002, 4));
        Config::enable_v1_1(c); c.flow = false; c.law.radiation = 0;
        c.evap_face = {0, 1, 0, 0};
        Owner o(c);
        o.initialize_temperature([](int, double, double) { return 3500.0; });
        if (ps == 0) { const auto st = o.state(); double a[PQ_NEL], x[PQ_NEL]; for (int el = 0; el < PQ_NEL; ++el) a[el] = pq_to_double(st.N[el], PQ_N_INV); pq_fractions(a, x); int s2 = 0; ps = pq_psat_total(&o.law(), 3500.0, x, &s2); }
        std::vector<double> pb(o.n_patches(), frac * ps);
        StepInput in; in.dt = 1e-13; in.em_power_W = {0.0}; in.back_pressure_Pa = pb;
        const auto r = o.advance(in);
        double n = 0; for (int el = 0; el < PQ_NEL; ++el) n += pq_to_double(r.evap_N[el], PQ_N_INV);
        evap.push_back(n);
        chk::check(o.ledger_closes(), "exact ledgers (Knudsen layer)");
        if (frac == 0.0) chk::put("knudsen_recoil_vacuum_Pa", o.patch_reports()[1].recoil_pressure_Pa);
    }
    double f0, f6, pr;
    const pq_law X = make_law(LawOptions::v1_1());
    pq_knudsen(&X, 0.0, &f0, &pr); pq_knudsen(&X, 0.6, &f6, &pr);
    chk::put_raw("knudsen_owner_evaporated_atoms_pb_0_0.6_1.5", jarr(evap));
    chk::put("knudsen_psat_3500K_Pa", ps);
    chk::check(evap[2] == 0.0, "no evaporation for a back pressure above saturation");
    chk::check(std::abs(evap[1] / evap[0] - f6 / f0) < 1e-6, "owner evaporation scales with the Knight net fraction");
}

// ------------------------------------------------------------------------------------------------ melt-layer benchmarks
// Constant-property liquid layer (thickness delta, T > Tm) on a solid substrate (T < Tm), heat conduction switched
// off (k ~ 0) so that the layer thickness is fixed; the face HighY carries the layer.
static Config layer_config(int nx, double W, double delta, double H, int rows_liq, int rows_sol, double sigma, double dsdT) {
    std::vector<double> rows;
    for (int j = 0; j <= rows_sol; ++j) rows.push_back((H - delta) * j / rows_sol);
    for (int j = 1; j <= rows_liq; ++j) rows.push_back(H - delta + delta * j / rows_liq);
    rows.back() = H;
    Config c = base(1, W, H, 1e-4, nx, 1, 1, rows);
    LawOptions o; o.kind = 1; o.c_T0 = 300; o.c_cps = o.c_cpl = 500; o.c_ks = o.c_kl = 1e-12; o.c_Tm = 1000; o.c_dTm = 1; o.c_L = 0; o.T_min = 200; o.T_max = 5000;
    o.c_rho = 7000; o.c_mu = 6e-3; o.c_sigma = sigma; o.c_dsdT = dsdT;
    c.law = o; c.rho_ref = 7000; c.evap_mode = 0; c.flow = true; c.T0 = 800;
    return c;
}
static void layer_profile(Owner& o, double H, double delta, std::function<double(double)> Tliq) {
    o.initialize_temperature([&](int, double x, double y) { return y > H - delta ? Tliq(x) : 800.0; });
}
// Mass gained by column i between two states, from exact atom-word differences.
static double moved_mass(const Owner& o, const HostState& s0, const HostState& s1, int i) {
    double m = 0;
    for (int el = 0; el < PQ_NEL; ++el) {
        pq_i128 d = pq_zero(); int ov = 0;
        for (int jj = 0; jj < o.dims().ny; ++jj) { const int cc = jj * o.dims().nx + i; d = pq_add(d, pq_sub(s1.N[PQ_NEL * cc + el], s0.N[PQ_NEL * cc + el], &ov), &ov); }
        m += pq_to_double(d, PQ_N_INV) * o.law().atom_mass[el];
    }
    return m;
}
static void poiseuille() {
    const double W = 2e-4, H = 1e-4, delta = 2e-5, j = 1e9, b = 10.0, mu = 6e-3, rho = 7000, dz = 1e-4;
    Config c = layer_config(2, W, delta, H, 4, 4, 0.0, 0.0);
    Owner o(c);
    layer_profile(o, H, delta, [](double) { return 1200.0; });
    const auto s0 = o.state();
    const double dt = 1e-9; const int steps = 1;   // one step: the layer has its initial thickness
    std::vector<double> JB(6 * o.dims().n_em, 0.0);
    for (int e = 0; e < o.dims().n_em; ++e) { JB[6 * e + 2] = j; JB[6 * e + 4] = -b; }   // J_z = j, B_y = -b: (J x B)_x = j b
    for (int k = 0; k < steps; ++k) { StepInput in; in.dt = dt; in.em_power_W.assign(o.dims().n_em, 0.0); in.em_JB = JB; o.advance(in); }
    const auto s1 = o.state();
    for (int k = 0; k < 200; ++k) { StepInput in; in.dt = dt; in.em_power_W.assign(o.dims().n_em, 0.0); in.em_JB = JB; o.advance(in); }
    chk::put("poiseuille_speed_after_201_steps_m_s", o.patch_reports()[1].melt_speed_m_s);
    const double moved = moved_mass(o, s0, s1, 1);
    const double expected = rho * (delta * delta * delta * j * b / (3 * mu)) * dz * dt * steps;
    chk::put("poiseuille_moved_kg", moved); chk::put("poiseuille_expected_kg", expected);
    chk::check(std::abs(moved / expected - 1) < 1e-9, "body-force melt flux q = delta^3 f/(3 mu) (Poiseuille layer)");
    chk::check(o.ledger_closes(), "exact ledgers with melt flow (Poiseuille)");
}
static void marangoni() {
    const double W = 2e-4, H = 1e-4, delta = 2e-5, mu = 6e-3, rho = 7000, dz = 1e-4, dsdT = -4e-4, G = 1e5;   // 100 K/mm
    Config c = layer_config(2, W, delta, H, 4, 4, 0.0, dsdT);
    Owner o(c);
    layer_profile(o, H, delta, [&](double x) { return 2000.0 + G * x; });
    const auto s0 = o.state();
    const double dt = 1e-10; const int steps = 1;   // one step: the surface temperatures are the initial ones
    for (int k = 0; k < steps; ++k) { StepInput in; in.dt = dt; in.em_power_W.assign(o.dims().n_em, 0.0); o.advance(in); }
    const auto s1 = o.state();
    for (int k = 0; k < 49; ++k) { StepInput in; in.dt = dt; in.em_power_W.assign(o.dims().n_em, 0.0); o.advance(in); }
    chk::put("marangoni_speed_after_50_steps_m_s", o.patch_reports()[1].melt_speed_m_s);
    const double moved = moved_mass(o, s0, s1, 1);
    const double tau = dsdT * G * (W / 2) / (W / 2);   // dsigma/dT x dT/dx
    const double expected = rho * (delta * delta * tau / (2 * mu)) * dz * dt * steps;
    chk::put("marangoni_moved_kg", moved); chk::put("marangoni_expected_kg", expected);
    chk::check(std::abs(moved / expected - 1) < 1e-9, "Marangoni (Couette) melt flux q = delta^2 tau/(2 mu)");
    chk::check(o.ledger_closes(), "exact ledgers with melt flow (Marangoni)");
}
// Capillary leveling of a cosine surface mode on the melt layer: discrete eigen-rate sigma delta^3 kappa^2/(3 mu)
// (kappa from the cell-centred Neumann second difference) and the continuum rate with kappa = (pi/W)^2.
static void leveling() {
    const int nx = 20; const double W = 2e-3, H = 1e-4, delta = 2e-5, sigma = 1.8, mu = 6e-3, A = 1e-8;
    Config c = layer_config(nx, W, delta, H, 4, 4, sigma, 0.0);
    Owner o(c);
    layer_profile(o, H, delta, [](double) { return 1200.0; });
    const double x0 = -0.5 * W;
    o.initialize_displacement([&](int, int f, double x) { return f == 1 ? A * std::cos(M_PI * (x - x0) / W) : 0.0; });
    const double dx = W / nx, kd = (2 - 2 * std::cos(M_PI / nx)) / (dx * dx), kc = M_PI * M_PI / (W * W);
    const double lam_d = sigma * delta * delta * delta * kd * kd / (3 * mu), lam_c = sigma * delta * delta * delta * kc * kc / (3 * mu);
    auto amplitude = [&]() {
        const auto eta = o.face_displacements();
        double num = 0, den = 0;
        for (int i = 0; i < nx; ++i) { const double sh = std::cos(M_PI * (i + 0.5) / nx); num += eta[nx + i] * sh; den += sh * sh; }   // face 1 columns
        return num / den;
    };
    const double a0 = amplitude();
    const double dt = 1e-3; const int steps = 100;
    std::int64_t nsub = 0;
    for (int k = 0; k < steps; ++k) { StepInput in; in.dt = dt; in.em_power_W.assign(o.dims().n_em, 0.0); nsub += o.advance(in).n_sub; }
    const double tau = dt / double(o.substeps_for(dt)), t = dt * steps;
    const double a1 = amplitude(), disc = a0 * std::pow(1 - lam_d * tau, double(nsub)), cont = a0 * std::exp(-lam_c * t);
    chk::put("leveling_initial_amplitude_m", a0); chk::put("leveling_final_amplitude_m", a1);
    chk::put("leveling_discrete_prediction_m", disc); chk::put("leveling_continuum_prediction_m", cont);
    chk::put("leveling_rate_discrete_per_s", lam_d); chk::put("leveling_rate_continuum_per_s", lam_c);
    chk::check(std::abs(a0 / A - 1) < 1e-6, "initial displacement mode recovered");
    chk::check(std::abs(a1 / disc - 1) < 1e-4, "capillary leveling follows the discrete eigen-rate (linear regime, A/delta = 5e-4)");
    chk::check(std::abs(a1 / cont - 1) < 0.01, "capillary leveling within 1 % of the lubrication-theory rate");
    chk::check(o.ledger_closes(), "exact ledgers with capillary flow");
}

// ------------------------------------------------------------------------------------------------ full v1.1 ledger
// SS304 v1.1 law, Mn, activities, Knudsen layer with back pressures, radiation, internal-energy channel, melt flow
// driven by Hall J x B, loads and deposits: exact ledgers every step; Mn depletion of the surface layer.
static void ledger() {
    const double W = 0.005, H = 0.0002;
    Config c = base(2, W, H, 2e-4, 5, 1, 2, graded_rows(H, 0.5e-6, 1.35, 10e-6, 1));
    Config::enable_v1_1(c);
    Owner o(c);
    o.initialize_temperature([&](int s, double x, double) { return 2400.0 + 2.0e5 * (x + 0.5 * W) + 30.0 * s; });
    const int nem = o.dims().n_em;
    std::vector<double> JB(6 * nem);
    // J along the foil (z), B mostly tangential (x) with a small normal component: tangential force J_z |B_y|,
    // normal magnetic load J_z B_x varying by 2 % across the width
    for (int e = 0; e < nem; ++e) { const double j = 5e8 * (1 + 0.02 * (e % 5)); JB[6 * e + 2] = j; JB[6 * e + 3] = 20.0; JB[6 * e + 4] = -0.05; }
    std::vector<double> pb(o.n_patches(), 0.0);
    for (int p = 0; p < o.n_patches(); ++p) pb[p] = (p % 3) * 2e3;
    double mn_surface0 = -1, mn_surface = 0, mn_bulk = 0;
    const auto t0 = Clock::now();
    bool exact = true; long steps = 0; std::uint64_t seq = 0; double maxspeed = 0, maxeta = 0, maxre = 0;
    for (int k = 0; k < 3000; ++k) {
        StepInput in; in.dt = 1e-11; in.em_power_W.assign(nem, 3e6); in.em_JB = JB; in.back_pressure_Pa = pb;
        for (int s = 0; s < 2; ++s) {
            ++seq;
            in.loads.push_back({s, -0.5 * W + 0.0011, H, 2e-9, "t/wall", "wall", seq});
            in.deposits.push_back({s, -0.5 * W + 0.0021, H, 3, 1e6, 1e-13, "t/mn-return", seq});
        }
        StepRecord r;
        try { r = o.advance(in); } catch (const PhaseRefusal& e) { chk::put_text("ledger_stop", e.what()); break; }
        ++steps;
        if (k % 100 == 0) exact = exact && o.ledger_closes();
        const auto rep = o.patch_reports();
        for (const auto& p : rep) { maxspeed = std::max(maxspeed, p.melt_speed_m_s); maxeta = std::max(maxeta, p.displacement_max_abs_m); maxre = std::max(maxre, p.reynolds_max); }
        if (mn_surface0 < 0) mn_surface0 = rep[1].surface_atom_fraction[3];
        mn_surface = rep[1].surface_atom_fraction[3];
    }
    { const auto N = o.atom_words(); double tot = 0; for (int el = 0; el < PQ_NEL; ++el) tot += pq_to_double(N[el], PQ_N_INV); mn_bulk = pq_to_double(N[3], PQ_N_INV) / tot; }
    exact = exact && o.ledger_closes();
    const auto& L = o.ledger();
    chk::check(exact, "exact energy and atom ledgers with all v1.1 options (melt flow, Mn, Knudsen, radiation, internal energy)");
    chk::check(steps == 3000, "3000 steps without refusal");
    chk::check(mn_surface < mn_surface0, "preferential Mn depletion of the evaporating surface cell");
    chk::put("ledger11_steps", double(steps)); chk::put("ledger11_seconds", std::chrono::duration<double>(Clock::now() - t0).count());
    chk::put("ledger11_max_T", o.max_temperature_K());
    chk::put("ledger11_Mn_surface_fraction_start", mn_surface0); chk::put("ledger11_Mn_surface_fraction_end", mn_surface); chk::put("ledger11_Mn_bulk_fraction", mn_bulk);
    double evap = 0, evap_mn = 0; for (int el = 0; el < PQ_NEL; ++el) { evap += pq_to_double(L.evap_N[el], PQ_N_INV); } evap_mn = pq_to_double(L.evap_N[3], PQ_N_INV);
    chk::put("ledger11_Mn_share_of_evaporated_atoms", evap_mn / evap);
    chk::put("ledger11_evap_internal_J", pq_to_double(L.evap_int, PQ_E_INV)); chk::put("ledger11_radiation_J", pq_to_double(L.radiation, PQ_E_INV));
    chk::put("ledger11_evap_reference_J", pq_to_double(L.evap_ref, PQ_E_INV)); chk::put("ledger11_joule_J", pq_to_double(L.joule, PQ_E_INV));
    chk::put("ledger11_max_melt_speed_m_s", maxspeed); chk::put("ledger11_max_displacement_m", maxeta); chk::put("ledger11_max_reynolds", maxre);
    chk::put_text("ledger11_state_sha256", o.state_sha256());
}

// ------------------------------------------------------------------------------------------------ refusals
static void geometry() {
    const double W = 2e-4, H = 1e-4, delta = 2e-5;
    Config c = layer_config(2, W, delta, H, 4, 4, 0.0, 0.0);
    c.geometry_refusal_fraction = 0.002;   // 0.2 um of the 100 um EM cell
    Owner o(c);
    layer_profile(o, H, delta, [](double) { return 1200.0; });
    std::vector<double> JB(6 * o.dims().n_em, 0.0);
    for (int e = 0; e < o.dims().n_em; ++e) { JB[6 * e + 2] = 1e9; JB[6 * e + 4] = -10.0; }
    std::string kind; int k = 0; double frac = 0;
    for (; k < 100000; ++k) {
        StepInput in; in.dt = 1e-9; in.em_power_W.assign(o.dims().n_em, 0.0); in.em_JB = JB;
        try { frac = o.advance(in).max_displacement_fraction; } catch (const PhaseRefusal& e) { kind = e.kind; break; }
    }
    chk::check(kind == "geometry-departure", "typed refusal when the displacement exceeds the stated fraction of the EM cell");
    chk::put("geometry_refusal_after_steps", k); chk::put("geometry_last_accepted_fraction", frac); chk::put_text("geometry_refusal_kind", kind);
}
static void explosion() {
    Config c = base(1, 0.001, 0.0002, 2e-4, 1, 1, 1, uniform_rows(0.0002, 2));
    Config::enable_v1_1(c); c.evap_mode = 0; c.law.radiation = 0; c.flow = false;
    c.law.T_critical = 4000.0;   // test value: phase-explosion temperature 3600 K
    Owner o(c);
    std::string kind; double T = 0;
    for (int k = 0; k < 100000; ++k) {
        StepInput in; in.dt = 1e-10; in.em_power_W = {2e5};
        try { const auto r = o.advance(in); T = r.max_T; } catch (const PhaseRefusal& e) { kind = e.kind; break; }
    }
    chk::check(kind == "phase-explosion" && T < 3600 && T > 3500, "typed refusal at 0.9 T_c");
    chk::put("explosion_last_accepted_T", T);
}
static void radiation() {
    Config c = base(1, 0.001, 0.0002, 2e-4, 1, 1, 1, uniform_rows(0.0002, 2));
    Config::enable_v1_1(c); c.evap_mode = 0; c.flow = false;
    Owner o(c);
    o.initialize_temperature([](int, double, double) { return 2000.0; });
    const double dt = 1e-12;
    StepInput in; in.dt = dt; in.em_power_W = {0.0};
    const auto r = o.advance(in);
    const double A = 2 * (0.001 * 2e-4) + 2 * (0.0002 * 2e-4);
    const double expect = 0.37 * 5.670374419e-8 * (std::pow(2000.0, 4) - std::pow(300.0, 4)) * A * dt;
    chk::put("radiation_J", r.J(r.radiation)); chk::put("radiation_expected_J", expect);
    chk::check(std::abs(r.J(r.radiation) / expect - 1) < 1e-9, "radiation loss eps sigma (T^4 - T_env^4) A dt");
    chk::check(o.ledger_closes(), "exact ledgers with radiation");
}
static void expansion() {
    Config c = base(1, 0.001, 0.0002, 2e-4, 1, 1, 1, uniform_rows(0.0002, 4));
    Config::enable_v1_1(c); c.evap_mode = 0; c.law.radiation = 0; c.flow = false;
    Owner o(c);
    o.initialize_temperature([](int, double, double) { return 2500.0; });
    const auto eta = o.face_displacements();
    const double rho0 = pq_density(&o.law(), 300.0);
    double m = 0; { const auto N = o.atom_words(); for (int el = 0; el < PQ_NEL; ++el) m += pq_to_double(N[el], PQ_N_INV) * o.law().atom_mass[el]; }
    const double expect = 0.5 * m * (1 / pq_density(&o.law(), 2500.0) - 1 / rho0) / (0.001 * 2e-4);
    chk::put("expansion_eta_m", eta[0]); chk::put("expansion_expected_m", expect);
    chk::check(std::abs(eta[0] / expect - 1) < 1e-9 && std::abs(eta[1] / expect - 1) < 1e-9, "thermal-expansion displacement from the density law");
}
static void merkle() {
    Config c = base(2, 0.002, 0.0002, 2e-4, 4, 1, 1, graded_rows(0.0002, 1e-6, 1.5, 20e-6, 1));
    Config::enable_v1_1(c);
    Owner o(c);
    const std::string d0 = o.state_sha256();
    chk::check(o.verify_state_leaves(), "Merkle leaves recomputed from the state words");
    StepInput in; in.dt = 1e-12; in.em_power_W.assign(o.dims().n_em, 1e5);
    o.advance(in);
    chk::check(o.state_sha256() != d0, "state digest changes after a step");
    chk::put("merkle_cell_leaves", double(o.state_roots().cell_leaves));
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
    timed("pulse", pulse);
    timed("knudsen", knudsen);
    timed("poiseuille", poiseuille);
    timed("marangoni", marangoni);
    timed("leveling", leveling);
    timed("radiation", radiation);
    timed("expansion", expansion);
    timed("merkle", merkle);
    timed("explosion", explosion);
    timed("geometry", geometry);
    timed("ledger", ledger);
    std::string name = "test_v11"; for (const auto& w : want) name += "-" + w;
    return chk::finish(name, out);
}
