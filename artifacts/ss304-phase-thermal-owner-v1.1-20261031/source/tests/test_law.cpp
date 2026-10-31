// test_law.cpp -- SS304 phase law v1: bit identity with the sealed law, source spot values, vapor cross-checks,
// deterministic math and the water-desorption cap analysis.
//   test_law OUT_DIR DATA_DIR
#include "check.hpp"
#include "phase_law_host.hpp"
#include "ss304_thermal_v2.hpp"   // sealed heat-v2 material (read in place)
#include <cfloat>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>

using namespace ss304_phase_v1_1;

static std::uint64_t ulp_distance(double a, double b) {
    if (a == b) return 0;
    std::int64_t ia, ib; std::memcpy(&ia, &a, 8); std::memcpy(&ib, &b, 8);
    if (ia < 0) ia = INT64_MIN - ia;
    if (ib < 0) ib = INT64_MIN - ib;
    return ia > ib ? std::uint64_t(ia - ib) : std::uint64_t(ib - ia);
}
static bool same_bits(double a, double b) { return std::memcmp(&a, &b, 8) == 0; }

struct JRow { double T, Cp, S, H, dfH, dfG, logKf; bool has_formation; };
static std::vector<JRow> read_janaf(const std::string& path) {
    std::ifstream f(path); std::string line; std::vector<JRow> out; int n = 0;
    while (std::getline(f, line)) {
        if (++n <= 2) continue;
        std::istringstream s(line); std::vector<std::string> p; std::string x; while (s >> x) p.push_back(x);
        if (p.size() < 5) continue;
        JRow r{}; try { r.T = std::stod(p[0]); r.Cp = std::stod(p[1]); r.S = std::stod(p[2]); r.H = std::stod(p[4]); } catch (...) { continue; }
        r.has_formation = false;
        if (p.size() >= 8) { try { r.dfH = std::stod(p[5]); r.dfG = std::stod(p[6]); r.logKf = std::stod(p[7]); r.has_formation = true; } catch (...) {} }
        out.push_back(r);
    }
    return out;
}

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : "", data = argc > 2 ? argv[2] : "data";
    LawOptions opt;
    const pq_law L = make_law(opt);
    chk::put_text("law_sha256", law_sha256(L));

    // ---------------------------------------------------------------- 1. bit identity with the sealed law
    {
        const ss304v2::Model q0(ss304v2::Qualification::evaluated_solid_1500, ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation, 8000.);
        const ss304v2::Model q1(ss304v2::Qualification::nominal_1800, ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation, 8000.);
        std::vector<double> Ts;
        for (int i = 0; i <= 400000; ++i) Ts.push_back(298.2 + (1800.0 - 298.2) * i / 400000.0);
        for (double t : {298.2, 299.0, 299.999999, 300.0, 350.0, 1500.0, 1600.0, 1672.0, 1700.0, 1727.0, 1799.5, 1800.0}) Ts.push_back(t);
        std::mt19937_64 rng(20261031);
        std::uniform_real_distribution<double> U(298.2, 1800.0);
        for (int i = 0; i < 200000; ++i) Ts.push_back(U(rng));
        long n1500 = 0, n1800 = 0, bad1500 = 0, bad1800 = 0;
        for (double T : Ts) {
            if (T > 1800.0) continue;
            const auto& m = T <= 1500.0 ? q0 : q1;
            const auto v = m.evaluate(T);
            double rho, k; pq_ss_transport(T, &rho, &k);
            const bool same = same_bits(pq_ss_h(T), v.h_j_kg) && same_bits(pq_ss_cap(T), v.apparent_heat_capacity_j_kg_k) && same_bits(rho, v.resistivity_ohm_m) &&
                              same_bits(k, v.k_w_m_k) && same_bits(1.0 / rho, v.sigma_s_m) && same_bits(pq_ss_liquid_fraction(T), v.liquid_fraction) &&
                              same_bits(8000.0 * pq_ss_h(T), v.energy_j_m3_reference);
            if (T <= 1500.0) { ++n1500; if (!same) ++bad1500; }
            // the same expressions also hold on 1500-1800 K against the sealed nominal_1800 qualification
            const auto w = q1.evaluate(T);
            const bool same1 = same_bits(pq_ss_h(T), w.h_j_kg) && same_bits(pq_ss_cap(T), w.apparent_heat_capacity_j_kg_k) && same_bits(rho, w.resistivity_ohm_m) && same_bits(k, w.k_w_m_k);
            ++n1800; if (!same1) ++bad1800;
        }
        chk::check(bad1500 == 0, "bit identity with sealed evaluated_solid_1500 on 298.2-1500 K");
        chk::check(bad1800 == 0, "bit identity with sealed nominal_1800 on 298.2-1800 K");
        chk::put("bit_identity_points_298_1500", double(n1500)); chk::put("bit_identity_mismatches_298_1500", double(bad1500));
        chk::put("bit_identity_points_298_1800", double(n1800)); chk::put("bit_identity_mismatches_298_1800", double(bad1800));
        // inversion: new Newton/quadratic/linear inverse versus the sealed bisection inverse
        std::uint64_t max_ulp = 0; long inv_points = 0, worse = 0;
        for (std::size_t i = 0; i < Ts.size(); i += 3) {
            const double T = Ts[i];
            const double h = pq_ss_h(T);
            int st = 0; const double Tn = pq_T_of_h(&L, h, &st);
            const double Ts_ = (T <= 1500.0 ? q0 : q1).temperature_from_u(h);
            max_ulp = std::max(max_ulp, ulp_distance(Tn, Ts_));
            if (std::abs(pq_ss_h(Tn) - h) > std::abs(pq_ss_h(Ts_) - h)) ++worse;
            ++inv_points;
            chk::check(st == 0, "inverse status inside the domain");
        }
        chk::check(worse == 0, "inverse is never further from h than the sealed bisection");
        chk::put("inverse_points", double(inv_points)); chk::put("inverse_max_ulp_vs_sealed_bisection", double(max_ulp)); chk::put("inverse_worse_than_sealed", double(worse));
        chk::check(max_ulp <= 2, "inverse within 2 ulp of the sealed bisection");
    }
    // ---------------------------------------------------------------- 2. domain, continuity, inverse over 298.2-5000 K
    {
        std::mt19937_64 rng(7);
        std::uniform_real_distribution<double> U(298.2, 5000.0);
        long bad = 0;
        for (int i = 0; i < 300000; ++i) {
            const double T = i < 4702 ? 298.2 + i : U(rng);
            if (T > 5000.0) continue;
            int st = 0; const double h = pq_ss_h(T), Tn = pq_T_of_h(&L, h, &st);
            if (st != 0 || std::abs(pq_ss_h(Tn) - h) > std::abs(pq_ss_h(T) - h)) ++bad;
        }
        chk::check(bad == 0, "inverse returns the closest representable temperature over 298.2-5000 K");
        int st = 0;
        (void)pq_T_of_h(&L, pq_ss_h(5000.0) * (1 + 1e-12), &st); chk::check(st == PQ_ERR_T_HIGH, "typed status above 5000 K");
        st = 0; (void)pq_T_of_h(&L, pq_ss_h(298.2) - 1.0, &st); chk::check(st == PQ_ERR_T_LOW, "typed status below 298.2 K");
        st = 0; (void)pq_T_of_h(&L, NAN, &st); chk::check(st == PQ_ERR_NONFINITE, "typed status for nonfinite energy");
        for (double Tb : {1672.0, 1727.0, 1800.0}) {
            const double a = pq_ss_h(std::nextafter(Tb, 0.0)), b = pq_ss_h(Tb);
            chk::check(std::abs(b - a) < 1e-6, "enthalpy continuous at " + std::to_string(Tb));
        }
        double r1, k1, r2, k2; pq_ss_transport(1800.0, &r1, &k1); pq_ss_transport(std::nextafter(1800.0, 1e9), &r2, &k2);
        chk::check(std::abs(r2 / r1 - 1) < 1e-12 && std::abs(k2 / k1 - 1) < 1e-12, "transport continuous at 1800 K");
        chk::put("lorenz_ratio_anchor", (28.99 * 151.7e-8) / 1800.0);
        // monotone h
        long nonmono = 0; double prev = -1e300;
        for (int i = 0; i <= 470180; ++i) { const double T = 298.2 + i * 0.01; if (T > 5000) break; const double h = pq_ss_h(T); if (!(h > prev)) ++nonmono; prev = h; }
        chk::check(nonmono == 0, "enthalpy strictly increasing on a 0.01 K grid");
        std::ostringstream tab; tab << "[";
        bool first = true;
        for (double T : {298.2, 300.0, 500.0, 1000.0, 1500.0, 1672.0, 1700.0, 1727.0, 1800.0, 1900.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0}) {
            double r, k; pq_ss_transport(T, &r, &k);
            tab << (first ? "" : ",") << "{\"T\":" << chk::num(T) << ",\"h_J_kg\":" << chk::num(pq_ss_h(T)) << ",\"cap_J_kg_K\":" << chk::num(pq_ss_cap(T))
                << ",\"k_W_m_K\":" << chk::num(k) << ",\"rho_e_Ohm_m\":" << chk::num(r) << ",\"liquid_fraction\":" << chk::num(pq_ss_liquid_fraction(T))
                << ",\"density_kg_m3\":" << chk::num(pq_ss_density(T)) << "}";
            first = false;
        }
        tab << "]";
        chk::put_raw("law_table", tab.str());
    }
    // ---------------------------------------------------------------- 3. source spot values (Mills 2004 Table 3)
    {
        const double h298 = pq_ss_h(298.2);
        struct Row { double T, H_kJ; } solid[] = {{400, 50.7}, {500, 103}, {600, 156}, {700, 212}, {800, 269}, {900, 327}, {1000, 387}, {1100, 448}, {1200, 511}, {1300, 575}, {1400, 640}, {1500, 707}, {1600, 775}};
        double maxrel = 0;
        for (auto r : solid) maxrel = std::max(maxrel, std::abs((pq_ss_h(r.T) - h298) / 1000.0 - r.H_kJ) / r.H_kJ);
        chk::check(maxrel < 0.01, "solid enthalpy within 1% of Mills Table 3(a) (table rounding)");
        chk::put("mills_table3a_H_max_rel_dev", maxrel);
        // liquid: Mills uses T_liq = 1723 K; this law uses the producer liquidus 1727 K. Compare at 1773-1900 K.
        struct Row2 { double T, H_kJ; } liquid[] = {{1773, 1180}, {1800, 1202}, {1900, 1283}};
        double maxrel_l = 0;
        for (auto r : liquid) maxrel_l = std::max(maxrel_l, std::abs((pq_ss_h(r.T) - h298) / 1000.0 - r.H_kJ) / r.H_kJ);
        chk::check(maxrel_l < 0.05, "liquid enthalpy within the Mills 5% uncertainty (Table 3(b))");
        chk::put("mills_table3b_H_max_rel_dev", maxrel_l);
        chk::check(std::abs(pq_ss_density(298.2) - 7921) < 1e-9 && std::abs(pq_ss_density(1500) - 7320) < 0.2, "solid density Mills Table 3(a)");
        chk::check(std::abs(pq_ss_density(1773) - 6865) < 0.5 && std::abs(pq_ss_density(1900) - 6756) < 0.5, "liquid density Mills Table 3(b)");
        double r, k; pq_ss_transport(1900, &r, &k);
        chk::check(std::abs(k - 30) < 0.02, "liquid conductivity Mills Table 3(b) 1900 K");
        chk::check(std::abs(PQ_SS_LATENT - 276320.0) == 0, "latent heat 160 J/K/kg x 1727 K");
    }
    // ---------------------------------------------------------------- 4. vapor law
    {
        // transcription: embedded rows equal the JANAF text files
        const std::pair<const janaf::Table*, std::string> tabs[] = {{&janaf::fe_gas, "Fe-008"}, {&janaf::fe_liquid, "Fe-006"}, {&janaf::cr_gas, "Cr-005"},
                                                                    {&janaf::cr_liquid, "Cr-003"}, {&janaf::ni_gas, "Ni-005"}, {&janaf::ni_liquid, "Ni-003"},
                                                                    {&janaf::mn_gas, "Mn-005"}, {&janaf::mn_liquid, "Mn-003"}, {&janaf::si_gas, "Si-005"},
                                                                    {&janaf::si_liquid, "Si-003"}, {&janaf::c_gas, "C-003"}, {&janaf::c_graphite, "C-002"}};
        for (const auto& [t, id] : tabs) {
            const auto rows = read_janaf(data + "/janaf/" + id + ".txt");
            long found = 0;
            for (int i = 0; i < t->n; ++i) {
                bool ok = false;
                for (const auto& r : rows) if (r.T == t->rows[i].T && r.S == t->rows[i].S && r.H == t->rows[i].HmH298) { ok = true; break; }
                if (ok) ++found;
            }
            chk::check(found == t->n, "JANAF transcription " + id);
        }
        // above the normal boiling point the liquid table's log Kf is -log10(p/bar): compare
        const char* liq[3] = {"Fe-006", "Cr-003", "Ni-003"};
        double max_dlog_bp = 0;
        for (int el = 0; el < 3; ++el) {
            for (const auto& r : read_janaf(data + "/janaf/" + liq[el] + ".txt")) {
                if (!r.has_formation || r.logKf >= 0 || r.T < 3200 || r.T > 5000) continue;   // above the boiling point only
                int st = 0;
                const double lp = pq_lnp(&L, el, r.T, &st) / std::log(10.0) - 5.0;   // log10(p/bar)
                max_dlog_bp = std::max(max_dlog_bp, std::abs(lp + r.logKf));
            }
        }
        chk::check(max_dlog_bp < 0.002, "vapor pressure above the boiling point equals -logKf(liquid) (3-decimal table)");
        chk::put("max_abs_dlog10p_vs_janaf_logKf_above_bp", max_dlog_bp);
        // Alcock, Itkin, Horrigan 1984 (transcription: mmrc.caltech.edu table): log10(p/atm) = A + B/T + C log10 T + D T/1000
        auto alcock = [](double A, double B, double C, double D, double T) { return A + B / T + C * std::log10(T) + D * T / 1000.0 + std::log10(101325.0); };
        double maxfe = 0, maxni = 0;
        for (double T = 1810; T <= 2090; T += 10) { int st = 0; maxfe = std::max(maxfe, std::abs(pq_lnp(&L, 0, T, &st) / std::log(10.0) - alcock(6.347, -19574, 0, 0, T))); }
        for (double T = 1730; T <= 2140; T += 10) { int st = 0; maxni = std::max(maxni, std::abs(pq_lnp(&L, 2, T, &st) / std::log(10.0) - alcock(6.666, -20765, 0, 0, T))); }
        chk::check(maxfe < std::log10(1.10), "Fe liquid vapor pressure within 10% of Alcock et al. 1984 (fit range)");
        chk::check(maxni < std::log10(1.10), "Ni liquid vapor pressure within 10% of Alcock et al. 1984 (fit range)");
        chk::put("alcock_fe_liquid_max_abs_dlog10p", maxfe); chk::put("alcock_ni_liquid_max_abs_dlog10p", maxni);
        // Cr: Alcock gives the solid only (to 2130 K); the pure-liquid standard state lies above it by the fusion
        // Gibbs energy. Check the ratio p_liquid/p_crystal against JANAF Cr(cr,l) and Cr(l).
        double cr_ratio_min = 1e9, cr_ratio_max = 0;
        for (double T = 1700; T <= 2100; T += 100) {
            int st = 0;
            const double ratio = std::pow(10.0, pq_lnp(&L, 1, T, &st) / std::log(10.0) - alcock(6.8, -20733, 0.4391, -0.4094, T));
            cr_ratio_min = std::min(cr_ratio_min, ratio); cr_ratio_max = std::max(cr_ratio_max, ratio);
        }
        chk::put("cr_liquid_over_alcock_solid_ratio_min", cr_ratio_min); chk::put("cr_liquid_over_alcock_solid_ratio_max", cr_ratio_max);
        double cr_logkf_ratio_max = 0;
        {
            const auto lq = read_janaf(data + "/janaf/Cr-003.txt");
            for (const auto& r : lq) if (r.has_formation && r.T >= 1700 && r.T <= 2100) cr_logkf_ratio_max = std::max(cr_logkf_ratio_max, std::pow(10.0, -r.logKf));
        }
        chk::put("cr_janaf_liquid_over_crystal_ratio_max_1700_2100", cr_logkf_ratio_max);
        chk::check(cr_ratio_min > 0.9 && cr_ratio_max < 1.4, "Cr liquid-standard pressure close to Alcock solid (fusion Gibbs energy offset)");
        // latent heat implied by the constant eps0 versus JANAF dHvap
        std::ostringstream lt; lt << "[";
        double max_lat_dev = 0, max_lat_dev_4000 = 0, max_lat_dev_2000_3500 = 0;
        for (int el = 0; el < 3; ++el) {
            const auto k = janaf_knots(gas_table(el), condensed_table(el), 298.15, 5000.0);
            for (const auto& x : k) {
                if (x.T < 1700 || x.T > 5000) continue;
                const double lam = (L.eps0[el] + 2.5 * boltzmann * x.T - L.atom_mass[el] * pq_ss_h(x.T)) * avogadro;
                const double dev = lam / x.dHvap_J_mol - 1;
                max_lat_dev = std::max(max_lat_dev, std::abs(dev));
                if (x.T <= 4000) max_lat_dev_4000 = std::max(max_lat_dev_4000, std::abs(dev));
                if (x.T >= 2000 && x.T <= 3500) max_lat_dev_2000_3500 = std::max(max_lat_dev_2000_3500, std::abs(dev));
                if (int(x.T) % 500 == 0) lt << (lt.tellp() > 1 ? "," : "") << "{\"el\":\"" << element_name[el] << "\",\"T\":" << x.T << ",\"model_J_mol\":" << chk::num(lam) << ",\"janaf_J_mol\":" << chk::num(x.dHvap_J_mol) << ",\"rel_dev\":" << chk::num(dev) << "}";
            }
        }
        lt << "]";
        chk::put_raw("latent_heat_model_vs_janaf", lt.str());
        chk::put("latent_heat_max_rel_dev_1700_5000", max_lat_dev);
        chk::check(max_lat_dev_4000 < 0.06, "constant-eps0 latent heat within 6% of JANAF over 1700-4000 K");
        chk::check(max_lat_dev_2000_3500 < 0.03, "constant-eps0 latent heat within 3% of JANAF over 2000-3500 K");
        chk::put("latent_heat_max_rel_dev_2000_3500", max_lat_dev_2000_3500);
        chk::put("latent_heat_max_rel_dev_1700_4000", max_lat_dev_4000);
        for (int el = 0; el < 3; ++el) chk::put(std::string("eps0_eV_") + element_name[el], L.eps0[el] / 1.602176634e-19);
        // fluxes for the nominal composition
        const double xs[3] = {0.7102, 0.2021, 0.0877};
        std::ostringstream ft; ft << "[";
        for (double T : {1700.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0}) {
            double mf = 0, q = 0; int st = 0;
            ft << (T > 1700 ? "," : "") << "{\"T\":" << T;
            for (int el = 0; el < 3; ++el) {
                const double j = pq_flux(&L, el, T, xs[el], &st);
                mf += j * L.atom_mass[el]; q += j * (L.eps0[el] + 2 * boltzmann * T - L.atom_mass[el] * pq_ss_h(T));
                ft << ",\"p_" << element_name[el] << "_Pa\":" << chk::num(std::exp(pq_lnp(&L, el, T, &st))) << ",\"j_" << element_name[el] << "_m2s\":" << chk::num(j);
            }
            ft << ",\"mass_flux_kg_m2_s\":" << chk::num(mf) << ",\"evaporative_cooling_W_m2\":" << chk::num(q) << ",\"recession_m_s\":" << chk::num(mf / pq_ss_density(T)) << "}";
        }
        ft << "]";
        chk::put_raw("hkl_nominal_304_beta_0.18", ft.str());
    }
    // ---------------------------------------------------------------- 5. deterministic math and fixed point
    {
        std::mt19937_64 rng(11);
        std::uniform_real_distribution<double> U(-740, 709), V(-60, 60);
        std::uint64_t me = 0, ml = 0;
        for (int i = 0; i < 2000000; ++i) {
            const double x = i % 2 ? U(rng) : V(rng);
            me = std::max(me, ulp_distance(pq_exp(x), std::exp(x)));
            const double y = std::exp(V(rng) / 2);
            ml = std::max(ml, ulp_distance(pq_log(y), std::log(y)));
        }
        chk::put("pq_exp_max_ulp_vs_libm", double(me)); chk::put("pq_log_max_ulp_vs_libm", double(ml));
        chk::check(me <= 2 && ml <= 2, "pq_exp/pq_log within 2 ulp of glibc");
        int bad = 0, ov = 0;
        std::uniform_real_distribution<double> W(-1e6, 1e6);
        long rt = 0;
        for (int i = 0; i < 1000000; ++i) {
            const double x = W(rng) * std::pow(2.0, double(int(rng() % 45)) - 40);   // |x| < 2^24 J, inside the 2^31 J word range
            const pq_i128 q = pq_from_double(x, PQ_E_POW2, &bad);
            const double xs = std::ldexp(x, 96);
            const double back_d = pq_to_double(q, PQ_E_INV);
            if (xs == std::floor(xs) ? back_d != x : std::abs(std::ldexp(back_d, 96) - xs) > 0.5) ++rt;
            // divmod identity
            pq_u32 r = 0; const pq_u32 d = 1 + pq_u32(rng() % 1000000);
            const pq_i128 qq = pq_divmod_u32(q, d, &r, &ov);
            pq_i128 back = pq_zero();
            for (pq_u32 k = 0; k < 1; ++k) {}
            // back = qq*d + r computed by repeated doubling
            pq_i128 acc = pq_zero(), base = qq; pq_u32 dd = d;
            while (dd) { if (dd & 1) acc = pq_add(acc, base, &ov); base = pq_add(base, base, &ov); dd >>= 1; }
            pq_i128 rr; rr.lo = r; rr.hi = 0;
            back = pq_add(acc, rr, &ov);
            if (!pq_eq(back, q) || r >= d) ++rt;
        }
        chk::check(rt == 0 && !bad && !ov, "fixed-point round trip and exact floor division");
        int bad2 = 0; (void)pq_from_double(0x1p31, PQ_E_POW2, &bad2);
        chk::check(bad2 == 1, "fixed-point range refusal at 2^31 J");
        // rounding: ties to even
        chk::check(pq_from_double(0x1p-97, PQ_E_POW2, &bad).lo == 0 && pq_from_double(3 * 0x1p-97, PQ_E_POW2, &bad).lo == 2, "fixed-point ties to even");
    }
    // ---------------------------------------------------------------- 6. water desorption: cap at 3000 K
    {
        // Source law (heat-v2 water_scenario.hpp, ConditionalAnalytic298To3000): k = (kT/h) exp(-0.87 (1 - 0.5 theta) eV / kT).
        const double kB = 1.380649e-23, hP = 6.62607015e-34, qe = 1.602176634e-19;
        auto rate = [&](double T, double th) { return kB * T / hP * std::exp(-0.87 * (1 - 0.5 * th) * qe / (kB * T)); };
        std::ostringstream o; o << "[";
        bool first = true;
        for (double beta : {1e9, 1e10, 1e11, 1e12, 1e13, 1e14}) {
            double lnth = 0.0, T = 300.0;     // theta0 = 1 (largest coverage, lowest barrier)
            const double dT = 0.01;
            while (T < 3000.0 && lnth > -800) { const double Tm = T + 0.5 * dT; lnth -= rate(Tm, std::exp(lnth)) * dT / beta; T += dT; }
            o << (first ? "" : ",") << "{\"ramp_K_per_s\":" << chk::num(beta) << ",\"log10_theta_ratio_at_3000K\":" << chk::num(lnth > -800 ? lnth / std::log(10.0) : -347.0) << ",\"reached_floor\":" << (lnth <= -800 ? "true" : "false") << "}";
            first = false;
        }
        o << "]";
        chk::put_raw("water_survival_linear_ramps_from_300K_theta0_1", o.str());
        chk::put("water_rate_3000K_theta0_per_s", rate(3000, 0));
    }
    // ---------------------------------------------------------------- 7. v1.1 data comparisons
    {
        const pq_law X = make_law(LawOptions::v1_1());
        chk::put_text("law_v1_1_sha256", law_sha256(X));
        // 7a. pulse-heating liquid (Pichler et al. 2020, Table 6 as published; H relative to 298 K in kJ/kg,
        //     resistivity with volume correction in uOhm m, density in kg/m3; uncertainties 2.9 %, 4.5 %, 3 % at k = 2)
        struct Row { double T, H, rho_corr, D; } rows[] = {
            {1800, 1190, 1.495, 6857}, {1900, 1275, 1.514, 6791}, {2000, 1360, 1.532, 6725}, {2100, 1444, 1.551, 6659}, {2200, 1529, 1.570, 6593},
            {2300, 1614, 1.589, 6527}, {2400, 1699, 1.608, 6461}, {2500, 1783, 1.627, 6395}, {2600, 1868, 1.646, 6329}, {2700, 1953, 1.665, 6262},
            {2800, 2037, 1.684, 6196}};
        const pq_law V = make_law(LawOptions{});
        double dH11 = 0, dR11 = 0, dD11 = 0, dH1 = 0, dR1 = 0, dD1 = 0;
        std::ostringstream tab; tab << "[";
        for (const auto& r : rows) {
            const double H11 = (pq_h(&X, r.T) - pq_ss_h_solid(298.15)) / 1000.0, H1 = (pq_h(&V, r.T) - pq_ss_h_solid(298.15)) / 1000.0;
            const double R11 = pq_rhoe(&X, r.T) * 1e6, R1 = pq_rhoe(&V, r.T) * 1e6;
            const double D11 = pq_density(&X, r.T), D1 = pq_density(&V, r.T);
            dH11 = std::max(dH11, std::abs(H11 / r.H - 1)); dR11 = std::max(dR11, std::abs(R11 / r.rho_corr - 1)); dD11 = std::max(dD11, std::abs(D11 / r.D - 1));
            dH1 = std::max(dH1, std::abs(H1 / r.H - 1)); dR1 = std::max(dR1, std::abs(R1 / r.rho_corr - 1)); dD1 = std::max(dD1, std::abs(D1 / r.D - 1));
            tab << (tab.tellp() > 1 ? "," : "") << "{\"T\":" << r.T << ",\"H_data\":" << r.H << ",\"H_v11\":" << chk::num(H11) << ",\"H_v1\":" << chk::num(H1)
                << ",\"rho_data\":" << r.rho_corr << ",\"rho_v11\":" << chk::num(R11) << ",\"rho_v1\":" << chk::num(R1)
                << ",\"D_data\":" << r.D << ",\"D_v11\":" << chk::num(D11) << ",\"D_v1\":" << chk::num(D1) << ",\"k_v11\":" << chk::num(pq_k(&X, r.T)) << "}";
        }
        tab << "]";
        chk::put_raw("pulse_heating_comparison", tab.str());
        chk::put("pulse_heating_H_maxdev_v11", dH11); chk::put("pulse_heating_rho_maxdev_v11", dR11); chk::put("pulse_heating_D_maxdev_v11", dD11);
        chk::put("pulse_heating_H_maxdev_v1", dH1); chk::put("pulse_heating_rho_maxdev_v1", dR1); chk::put("pulse_heating_D_maxdev_v1", dD1);
        chk::check(dH11 < 0.029 && dR11 < 0.045 && dD11 < 0.03, "v1.1 liquid law within the pulse-heating uncertainties (1800-2800 K)");
        chk::check(std::abs(pq_k(&X, 1800.0) / 29.0 - 1) < 0.10, "Wiedemann-Franz liquid k within the Mills 10 % at 1800 K");
        chk::put("lorenz_sommerfeld", X.ph_lorenz);
        // 7b. vapor pressure equals 1 bar at the JANAF normal boiling points (FUGACITY = 1 bar rows)
        struct Bp { int el; double T; } bps[] = {{0, 3133.345}, {1, 2952.078}, {2, 3156.584}, {3, 2334.526}, {4, 3504.616}};
        double maxbp = 0;
        for (const auto& b : bps) { int st = 0; maxbp = std::max(maxbp, std::abs(pq_lnp(&X, b.el, b.T, &st) / std::log(10.0) - 5.0)); }
        chk::put("boiling_point_max_abs_dlog10p", maxbp);
        chk::check(maxbp < 0.004, "vapor pressure 1 bar at the JANAF boiling points of Fe, Cr, Ni, Mn, Si (within 1 %)");
        // Mn against Alcock et al. 1984 (solid Mn, 1300-1500 K) with the JANAF liquid/crystal fusion ratio
        // 7c. activity: gamma_Ni at infinite dilution in liquid Fe at 1873 K = 0.66 (Zellars et al. 1960)
        double x0[PQ_NEL] = {1, 0, 0, 0, 0, 0};
        chk::check(std::abs(pq_gamma(&X, 2, 1873.0, x0) - 0.66) < 1e-12, "regular solution reproduces gamma_Ni infinite dilution 0.66 at 1873 K");
        double x7[PQ_NEL] = {0.7, 0, 0.3, 0, 0, 0};
        chk::put("gamma_Fe_at_xFe_0.7_1873K", pq_gamma(&X, 0, 1873.0, x7));
        // 7d. internal-energy channel: eps0 + E_int + 5/2 kT - m h(T) equals the JANAF vaporization enthalpy at the rows
        double max_int_dev = 0, min_eint = 1e300;
        for (int el = 0; el < PQ_NEL; ++el) {
            const auto k = janaf_knots(gas_table(el), condensed_table(el), 298.15, 6000.0);
            for (std::size_t i = 0; i < k.size(); ++i) {
                if (k[i].T < X.T_min || k[i].T > X.T_max) continue;
                const double lam = (X.eps0[el] + pq_eint(&X, el, k[i].T) + 2.5 * boltzmann * k[i].T - X.atom_mass[el] * pq_h(&X, k[i].T)) * avogadro;
                max_int_dev = std::max(max_int_dev, std::abs(lam / k[i].dHvap_J_mol - 1));
                min_eint = std::min(min_eint, pq_eint(&X, el, k[i].T));
            }
        }
        chk::put("internal_channel_latent_max_rel_dev", max_int_dev); chk::put("internal_channel_min_eint_J", min_eint);
        chk::check(max_int_dev < 1e-9 && min_eint >= 0, "with the internal-energy channel the latent heat equals JANAF at every row, E_int >= 0");
        // 7e. Knudsen layer (Knight 1979): sonic values and limits
        const auto s1 = knight(1.0);
        chk::put("knight_M1_Tv_Ts", s1.Tr); chk::put("knight_M1_rho_ratio", s1.rhor); chk::put("knight_M1_p_ratio", s1.pr);
        chk::put("knight_M1_net_fraction", s1.fnet); chk::put("knight_M1_recoil_ratio", s1.prec);
        chk::check(std::abs(s1.Tr - 0.669) < 0.001 && std::abs(s1.rhor - 0.308) < 0.001 && std::abs(s1.pr - 0.206) < 0.001, "Knight sonic jump 0.669 / 0.308 / 0.206");
        chk::check(std::abs(s1.fnet - 0.82) < 0.01, "sonic net flux fraction 0.82 (18 % recondensation)");
        double f, pr;
        pq_knudsen(&X, 0.0, &f, &pr); const double fv = f;
        pq_knudsen(&X, 1.0, &f, &pr); const double fh = f;
        pq_knudsen(&X, 0.6, &f, &pr); const double fm = f;
        chk::put("knudsen_fnet_vacuum", fv); chk::put("knudsen_fnet_back_0.6", fm); chk::put("knudsen_fnet_back_1.0", fh);
        chk::check(fh == 0.0 && fm > 0 && fm < fv, "net flux vanishes at saturation back pressure and decreases monotonically");
        long nonmono = 0; for (int i = 1; i < PQ_KL_N; ++i) if (!(X.kl_pv[i] < X.kl_pv[i - 1])) ++nonmono;
        chk::check(nonmono == 0, "Knudsen-layer pressure ratio strictly decreasing in Mach number");
        // 7f. Mills 2004 viscosity (Table 3(b)) and surface tension (Eq. 23)
        chk::check(std::abs(pq_viscosity(1773) * 1e3 - 8.07) < 0.01 && std::abs(pq_viscosity(1800) * 1e3 - 7.54) < 0.01 && std::abs(pq_viscosity(1900) * 1e3 - 5.98) < 0.01, "viscosity reproduces Mills Table 3(b)");
        double g, dg; pq_surface_tension(&X, 1823.0, &g, &dg);
        chk::check(std::abs(g - 1.840) < 1e-12 && dg == -0.4e-3, "surface tension Mills Eq. 23 (S = 0)");
        // 7g. HKL with Mn and activities (nominal 304 with Mn, vacuum, Knight sonic)
        std::array<double, PQ_NEL> xs{};
        { double tot = 0; for (int el = 0; el < PQ_NEL; ++el) { xs[el] = nominal_mass_fraction[el] / molar_mass_g[el]; tot += xs[el]; } for (auto& v : xs) v /= tot; }
        std::ostringstream ft; ft << "[";
        for (double T : {1700.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 5000.0, 6000.0}) {
            int st = 0; const double fn = pq_flux_fraction(&X, T, xs.data(), 0.0, &st);
            double tot = 0, jmn = 0, mf = 0;
            for (int el = 0; el < PQ_NEL; ++el) { const double j = pq_flux_f(&X, el, T, xs.data(), fn, &st); tot += j; mf += j * X.atom_mass[el]; if (el == 3) jmn = j; }
            ft << (T > 1700 ? "," : "") << "{\"T\":" << T << ",\"mass_flux_kg_m2_s\":" << chk::num(mf) << ",\"Mn_atom_share\":" << chk::num(jmn / tot)
               << ",\"recoil_Pa\":" << chk::num(pq_recoil(&X, T, xs.data(), 0.0, &st)) << ",\"psat_Pa\":" << chk::num(pq_psat_total(&X, T, xs.data(), &st)) << "}";
        }
        ft << "]";
        chk::put_raw("hkl_v11_nominal_with_Mn", ft.str());
        chk::put("phase_explosion_T_K", X.T_explosion);
    }
    return chk::finish("test_law", out);
}
