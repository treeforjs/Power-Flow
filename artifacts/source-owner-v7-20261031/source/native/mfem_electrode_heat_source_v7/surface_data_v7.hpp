#pragma once
// Heat/source owner v7 -- the cited SS304 surface parameter set ("ss304-oxide-water-v7") and the v6.1
// benchmark set. Every number is entered once through Registry::use with its tag, band, unit, source id
// and location; DATA-v7.txt lists the full citations of the source ids and the access level. Values read
// from published figures carry "fig" in the location (reading error stated there).
#include "surface_model_v7.hpp"

namespace surface_v7 {

// ------------------------------------------------------------------ electron wall model (Furman-Pivi)
inline void furman_pivi_ss(Registry& R, SW7Params& p) {
  const char* T1 = "PRST-AB 5 124404 Table I p.124404-10 (stainless column); fit to R. Kirby data, SLAC 304 rolled sheet, etched + passivated, unconditioned, 0-1100 eV";
  p.P1e_inf = R.fixed("fp.P1e_inf", 0.07, "1", Tag::recommended, FP2002, T1);
  p.P1e_hat = R.use("fp.P1e_hat", 0.5, 0.02, 1.0, "1", Tag::recommended, FP2002, T1,
                    "zero-energy reflectivity extrapolated; band: 0.02-0.07 (image-potential theory, Andronov/Kaganovich) to 1 (Cimino et al. PRL 93 014801, scrubbed Cu, 10 K); no stainless data below 20 eV");
  p.Ee_hat = R.fixed("fp.Ee_hat", 0, "eV", Tag::recommended, FP2002, T1);
  p.W = R.fixed("fp.W", 100, "eV", Tag::recommended, FP2002, T1);
  p.p = R.fixed("fp.p", 0.9, "1", Tag::recommended, FP2002, T1);
  p.sigma_e = R.fixed("fp.sigma_e", 1.9, "eV", Tag::recommended, FP2002, T1);
  p.e1 = R.fixed("fp.e1", 0.26, "1", Tag::recommended, FP2002, "Table I; angular parameters from TiN/Al data (p. 9), not stainless");
  p.e2 = R.fixed("fp.e2", 2, "1", Tag::recommended, FP2002, "Table I; from TiN/Al data");
  p.P1r_inf = R.fixed("fp.P1r_inf", 0.74, "1", Tag::recommended, FP2002, T1);
  p.Er = R.fixed("fp.Er", 40, "eV", Tag::recommended, FP2002, T1);
  p.r = R.fixed("fp.r", 1, "1", Tag::recommended, FP2002, T1);
  p.q = R.fixed("fp.q", 0.4, "1", Tag::recommended, FP2002, T1);
  p.r1 = R.fixed("fp.r1", 0.26, "1", Tag::recommended, FP2002, "Table I; from TiN/Al data");
  p.r2 = R.fixed("fp.r2", 2, "1", Tag::recommended, FP2002, "Table I; from TiN/Al data");
  p.dts_hat = R.fixed("fp.dts_hat", 1.22, "1", Tag::recommended, FP2002, T1);
  p.Ets_hat = R.fixed("fp.Ets_hat", 310, "eV", Tag::recommended, FP2002, T1);
  p.s = R.fixed("fp.s", 1.813, "1", Tag::recommended, FP2002, T1);
  p.t1 = R.fixed("fp.t1", 0.66, "1", Tag::recommended, FP2002, "Table I; from TiN/Al data");
  p.t2 = R.fixed("fp.t2", 0.8, "1", Tag::recommended, FP2002, "Table I; from TiN/Al data");
  p.t3 = R.fixed("fp.t3", 0.7, "1", Tag::recommended, FP2002, "Table I; from TiN/Al data");
  p.t4 = R.fixed("fp.t4", 1, "1", Tag::recommended, FP2002, "Table I; from TiN/Al data");
  const double pn[10] = {1.6, 2, 1.8, 4.7, 1.8, 2.4, 1.8, 1.8, 2.3, 1.8};
  const double en[10] = {3.9, 6.2, 13, 8.8, 6.25, 2.25, 9.2, 5.3, 17.8, 10};
  for (int n = 0; n < 10; ++n) {
    p.pn[n] = R.fixed("fp.p_n" + std::to_string(n + 1), pn[n], "1", Tag::recommended, FP2002, "Table II p.124404-11 (stainless); pinned only for 4 <= n <= 8 (p. 11)");
    p.epsn[n] = R.fixed("fp.eps_n" + std::to_string(n + 1), en[n], "eV", Tag::recommended, FP2002, "Table II p.124404-11 (stainless)");
  }
  p.M = int(R.fixed("fp.M", 10, "1", Tag::recommended, FP2002, "Eq. (46), Table II"));
  p.see_scale = R.use("see.scale_true_secondary", 1.0, 0.5, 1.2, "1", Tag::scenario, BAGLIN2000,
                      "1 = as-received fit; conditioned/scrubbed SS reaches delta_max 1.35-1.5 (Ji et al. NAPAC2016 Table 1; Hilleret ICFA 2000 slide 12)");
  p.fp_fit_E_max = R.fixed("fp.fit_E_max", 1100, "eV", Tag::recommended, FP2002, "Sec. VI, Figs. 4-5: data 0-1100 eV");
  p.bs_blend_E = R.fixed("backscatter.blend_E", 3930, "eV", Tag::model_derived, THIS,
                         "FP plateau P1e_inf+P1r_inf=0.81 vs measured eta_Fe~0.28: log-linear blend 1.1 keV -> first Hunger-Kuchler point 3.93 keV");
  const double bE[11] = {1000, 1500, 2000, 3000, 3930, 5000, 7500, 10000, 15200, 20200, 31000};
  const double bY[11] = {0.237, 0.263, 0.270, 0.277, 0.294, 0.284, 0.287, 0.285, 0.278, 0.276, 0.277};
  const char* bsrc[11] = {"Palluel 1947", "Palluel 1947", "Palluel 1947", "Palluel 1947", "Hunger & Kuchler 1979", "Hunger & Kuchler 1979", "Hunger & Kuchler 1979",
                          "Hunger & Kuchler 1979", "Hunger & Kuchler 1979", "Hunger & Kuchler 1979", "Hunger & Kuchler 1979"};
  p.bs_n = 11;
  for (int i = 0; i < 11; ++i) {
    p.bs_E[i] = bE[i];
    p.bs_eta[i] = R.use("backscatter.eta_Fe@" + std::to_string(int(bE[i])) + "eV", bY[i], 0.219, 0.307, "1", Tag::measured, JOY_DB,
                        std::string("Section 2 Iron (Z=26): ") + bsrc[i] + "; spread of eight data sets 0.251-0.307 at 4-30 keV; Fe proxy for 304");
  }
}

// ------------------------------------------------------------------ ion wall model
inline void ion_model(Registry& R, SW7Params& p) {
  p.pe_a = R.fixed("iiee.pe.a", 0.032, "1/eV", Tag::recommended, "baragiola_1979", "Surf. Sci. 90, 240, p. 241 and Fig. 1: gamma_p = 0.032 (0.78 E_i - 2 phi); noble-gas fit, upper bound for molecular ions (p. 248: molecular PEE 'very small')");
  p.pe_b = R.fixed("iiee.pe.b", 0.78, "1", Tag::recommended, "baragiola_1979", "p. 241");
  p.pe_c = R.fixed("iiee.pe.c", 2, "1", Tag::recommended, "baragiola_1979", "p. 241 (requires E_i > 2 phi)");
  // Kinetic emission per atom vs energy per atom: Alonso et al. Radiat. Eff. 45 (1979) 119, Fig. 3 (fig, +-5-10%),
  // Inconel 625 with native oxide (proxy for oxidized 304). Threshold 40 eV (20-60 eV for gas-covered surfaces, p. 119).
  const double hE[7] = {40, 1200, 2500, 5000, 10000, 20000, 30000};
  const double hG[7] = {0, 0.72, 0.97, 1.5, 1.8, 2.3, 2.5};
  const char* hS[7] = {"threshold (model_derived): 20-60 eV gas-covered, p. 119", "H2+ per atom at E/2, oxidized (fig)", "H2+ per atom, oxidized (fig)",
                       "H+ oxidized (fig)", "H+ oxidized (fig)", "H+ oxidized (fig)", "H+ oxidized (fig)"};
  p.ke_H_n = 7;
  for (int i = 0; i < 7; ++i) {
    p.ke_H_E[i] = R.use("iiee.ke.H.E" + std::to_string(i), hE[i], i ? hE[i] : 20, i ? hE[i] : 60, "eV", i ? Tag::measured : Tag::model_derived, "alonso_1979", hS[i]);
    p.ke_H_g[i] = R.use("iiee.ke.H.gamma" + std::to_string(i), hG[i], hG[i] * 0.9, hG[i] * 1.1, "1", i ? Tag::measured : Tag::model_derived, "alonso_1979", hS[i]);
  }
  const double oE[6] = {40, 2500, 5000, 10000, 20000, 30000};
  const double oG[6] = {0, 1.35, 2.05, 3.0, 4.2, 5.3};
  p.ke_O_n = 6;
  for (int i = 0; i < 6; ++i) {
    p.ke_O_E[i] = R.use("iiee.ke.O.E" + std::to_string(i), oE[i], i ? oE[i] : 20, i ? oE[i] : 60, "eV", i ? Tag::measured : Tag::model_derived, "alonso_1979",
                        i ? "O+ oxidized Inconel 625 (fig)" : "threshold (model_derived)");
    p.ke_O_g[i] = R.use("iiee.ke.O.gamma" + std::to_string(i), oG[i], oG[i] * 0.9, oG[i] * 1.1, "1", i ? Tag::measured : Tag::model_derived, "alonso_1979",
                        i ? "O+ oxidized Inconel 625 (fig)" : "threshold");
  }
  p.ke_scale = R.use("iiee.ke.scale", 1.0, 0.5, 1.0, "1", Tag::scenario, "baragiola_1979",
                     "molecular additivity gamma(H2O+,E) = gamma_O(16E/18) + 2 gamma_H(E/18) is an upper bound (~20%, pp. 247-248); clean surfaces about half (Alonso Fig. 4)");
  const double rn[6] = {0.8250, 21.41, 8.606, 0.6425, 1.907, 1.927}, re[6] = {0.6831, 27.16, 15.66, 0.6598, 7.967, 1.822};
  const double ron[6] = {0.3680, 2.985, 7.122, 0.5802, 4.211, 1.597}, roe[6] = {0.2058, 3.848, 19.07, 0.4872, 15.13, 1.638};
  for (int i = 0; i < 6; ++i) {
    p.rn_A[i] = R.fixed("reflection.H.RN.A" + std::to_string(i + 1), rn[i], "1", Tag::recommended, "thomas_janev_smith_1991", "INDC(NDS)-249 Table 2 p. 29, mu >= 20; RMS 12 %, data within +-25 %");
    p.re_A[i] = R.fixed("reflection.H.RE.A" + std::to_string(i + 1), re[i], "1", Tag::recommended, "thomas_janev_smith_1991", "Table 2 p. 29, mu >= 20");
    p.rn_O_A[i] = R.fixed("reflection.O.RN.A" + std::to_string(i + 1), ron[i], "1", Tag::model_derived, "thomas_janev_smith_1991", "Table 2 p. 29, mu ~ 3 set (database 'very poor'); O on Fe mu = 3.49");
    p.re_O_A[i] = R.fixed("reflection.O.RE.A" + std::to_string(i + 1), roe[i], "1", Tag::model_derived, "thomas_janev_smith_1991", "Table 2 p. 29, mu ~ 3 set");
  }
  p.reflect_eps_min = R.fixed("reflection.eps_min", 1e-3, "1", Tag::recommended, "thomas_janev_smith_1991", "validity 1e-3 <= eps <= 40");
  p.reflect_eps_max = R.fixed("reflection.eps_max", 40, "1", Tag::recommended, "thomas_janev_smith_1991", "validity range");
  p.reflect_scale = R.use("reflection.scale", 1.0, 0.5, 1.0, "1", Tag::scenario, "thomas_janev_smith_1991", "oxide layers and retained hydrogen lower R (p. 5, 5.3 p. 23)");
  // Ion-stimulated desorption of adsorbed water: no H2O data for any projectile (open). Scenario: eta ~ 1 molecule
  // per keV ion at one monolayer (unbaked SS ~10x baked, attributed to ~1 ML water: Edwards IEEE TNS 26 (1979) 4089
  // Table I), sigma = eta / N_ML with N_ML = 1e19 m^-2; energy shape of the CO yield (Mathewson CERN-ISR-VA/76-5 via
  // Malyshev CAS 2017 slide 33, fig): 0.5/1/2/3 keV -> 0.71/1/1.375/1.5; threshold 10 eV (linear below 500 eV).
  const double iE[5] = {10, 500, 1000, 2000, 3000};
  const double iS[5] = {0, 0.71e-19, 1.0e-19, 1.375e-19, 1.5e-19};
  p.isd_n = 5;
  for (int i = 0; i < 5; ++i) {
    p.isd_E[i] = iE[i];
    p.isd_sigma[i] = R.use("isd.sigma@" + std::to_string(int(iE[i])) + "eV", iS[i], iS[i] * 0.01, i ? iS[i] * 10 : 0, "m^2", Tag::scenario, "edwards_1979+mathewson_1976",
                           "scenario: no measured ISD yield of H2O; band two decades down, one up");
  }
  p.isd_scale = R.use("isd.scale", 1.0, 0.0, 10.0, "1", Tag::scenario, THIS, "sensitivity factor on the ISD scenario");
}

// ------------------------------------------------------------------ neutral wall and ESD
inline void neutral_esd_model(Registry& R, SW7Params& p) {
  p.stick0 = R.use("neutral.H2O.sticking", 0.01, 9e-4, 0.2, "1", Tag::recommended, "tuzi_1997",
                   "Shinku 40, 377, Fig. 3 p. 380: electropolished SS 0.005-0.01 during RT pumping; as-received glass-bead SS 9e-4; degassed 0.1-0.2 (Fig. 4)");
  p.stick_theta_max = R.fixed("neutral.H2O.sticking_theta_max", 1.0, "ML", Tag::model_derived, THIS, "Langmuir site blocking s = s0 (1 - theta/theta_max); no multilayer adsorption at T >= 298 K");
  p.accommodation = R.use("neutral.accommodation.molecular", 0.9, 0.8, 1.0, "1", Tag::recommended, "trott_2011",
                          "SAND2011-5640P p. 17 / RSI 82 035120: 304 SS untreated N2 0.87, Ar 0.95; water 'slightly larger' than the noble-gas trend (IJHMT 183 122195, abstract)");
  p.accommodation_light = R.use("neutral.accommodation.H_H2", 0.3, 0.15, 0.6, "1", Tag::recommended, "chiggiato_2015",
                                "arXiv:1501.07162 Table 7 p. 11: technical Pt H2 0.29, clean Pt 0.15; adsorbed layers raise alpha");
  // ESD cross section of adsorbed H2O per molecule, sigma(E) (eta = sigma N_ads; first-order loss).
  // Threshold ~10 eV (Lane et al. JCP 127 224706; Drinkwine & Lichtman p. 7); 100 eV: >= 8e-17 cm^2 (Petrik & Kimmel
  // PCCP 2018, alpha-Al2O3 D2O monolayer depleted at 3e16 cm^-2; first-order inference, model_derived); 300 eV:
  // 2e-16 cm^2 (ice 1.6-5.2e-16 cm^2 at 200-300 eV, Abdulgalil 2017; technical SS initial total 1e-16 cm^2, McCracken
  // 1967 via Drinkwine & Lichtman pp. 18-20); held above 300 eV (SS all-gas yields keep rising to 6.5 keV, ice peaks
  // at 300-400 eV: outside the H2O data, flagged).
  const double eE[4] = {10, 50, 100, 300};
  const double eS[4] = {0, 3.6e-21, 8e-21, 2e-20};
  const char* eL[4] = {"threshold ~10 eV", "linear between threshold and 100 eV (TiO2 monotonic rise 5-50 eV)", "Al2O3 monolayer inference (model_derived)", "ice / technical SS"};
  p.esd_n = 4;
  for (int i = 0; i < 4; ++i) {
    p.esd_E[i] = eE[i];
    p.esd_sigma[i] = R.use("esd.H2O.sigma@" + std::to_string(int(eE[i])) + "eV", eS[i], eS[i] * 1e-4, eS[i] * 2.6, "m^2", Tag::model_derived, "esd_compilation_v7", eL[i]);
  }
  p.esd_scale = R.use("esd.scale", 1.0, 1e-4, 2.6, "1", Tag::scenario, "drinkwine_lichtman",
                      "long bombardment lowers sigma to ~1e-20 cm^2 (1e-4 of initial); ice up to 5.2e-16 cm^2");
}

// ------------------------------------------------------------------ photo hook (optional, off by default)
inline PhotoModel photo_model(Registry& R) {
  PhotoModel m;
  const double E[9] = {5.1, 6.9, 9.7, 10.3, 12.1, 14, 19, 20, 25};
  const double Y[9] = {5e-7, 1e-4, 1e-2, 0.025, 0.11, 0.19, 0.19, 0.155, 0.125};
  for (int i = 0; i < 9; ++i) {
    m.E_eV.push_back(E[i]);
    m.yield.push_back(R.use("photo.yield@" + std::to_string(E[i]).substr(0, 4) + "eV", Y[i], Y[i] / 1.12, Y[i] * 1.12, "e/photon", Tag::measured, GRARD1973,
                            "Fig. 5a stainless (Feuerbacher & Fitton, J. Appl. Phys. 43, 1563 (1972)); fig, +-0.05 decade; surface preparation not verified"));
  }
  m.enabled = false;
  return m;
}

// ------------------------------------------------------------------ assembled sets
inline SW7Params wall_parameters(Registry& R) {
  SW7Params p{};
  furman_pivi_ss(R, p);
  ion_model(R, p);
  neutral_esd_model(R, p);
  p.enable_see = p.enable_iiee = p.enable_reflection = p.enable_sticking = p.enable_esd = p.enable_isd = 1;
  return p;
}
inline WallModel ss304_wall_model(Registry& R) {
  WallModel w;
  w.enabled = true;
  w.p = wall_parameters(R);
  w.species = species_table(&R);
  return w;
}


// ------------------------------------------------------------------ work function of oxidized SS304 with adsorbed water
inline SurfaceModel ss304_work_function(Registry& R) {
  SurfaceModel m;
  m.benchmark = false;
  m.fixed_phi_eV = R.fixed("wf.benchmark_clean", 4.3, "eV", Tag::benchmark, V61, "v6.1 clean-surface benchmark (selectable: SurfaceModel::benchmark)");
  m.phi_ox_eV = R.use("wf.oxide", 4.35, 4.19, 5.3, "eV", Tag::recommended, "walton_1999",
                      "J. Appl. Phys. 85, 1832, p. 1833 Fig. 1: clean 302 SS 4.75 eV (photoemission at 23 eV, sputter + anneal), -0.4 eV after 60 L O2 (initial oxide); "
                      "band: 4.19 +- 0.02 eV sputtered 304 (Marlow et al. JAP 133 085104 (2023) Table II, UPS) to 5.3 eV (air-exposed literature range 4.20-5.3, secondary via Marlow refs 14-17); "
                      "native passive film (2-5 nm) not measured");
  m.band_low_eV = 4.19;
  m.band_high_eV = 5.3;
  m.dphi_1ML_eV = R.use("wf.water_dphi_1ML", -0.97, -1.20, -0.60, "eV", Tag::recommended, "thiel_madey_1987",
                        "Surf. Sci. Rep. 7, 211, Table 5 p. 250: molecular H2O on polycrystalline Fe, Delta phi_max -0.97 +- 0.08 eV (77 K); band: Co -1.20 to Au -0.60 (same table); "
                        "no data for water on Fe/Cr oxides (open); sign negative in every reported case (p. 249)");
  m.depolarization_c = R.use("wf.water_depolarization_c", 1.2, 0.0, 3.0, "1", Tag::model_derived, "thiel_madey_1987",
                             "initial slope d(Delta phi)/d(theta) = -N_ML mu0/eps0 (Helmholtz) with mu0 = 0.56 D (Ni polycr., Table 5) and N_ML = 1e19 m^-2 gives -2.1 eV/ML = "
                             "Delta phi_1ML (1 + c); band: c = 0 (linear) to 3 (Cu(110) 0.87 D, N_ML 1.3e19)");
  m.dphi_dT_eV_K = R.use("wf.dphi_dT", 0.0, -4.6e-4, 4.6e-4, "eV/K", Tag::scenario, "tolias_2022",
                         "no SS/Fe/Cr/oxide value; 'a few tens to few hundred ueV/K, either sign, generally neglected' (Tolias, Nucl. Mater. Energy 31 101195 (2022) sec. 2.9); Ag -4.58e-4 (Akbi 2014, abstract); "
                         "extrapolation to 5000 K not supported: temperature acts through the coverage (desorption) and the phase only");
  m.phi_liquid_eV = R.use("wf.liquid_steel", 4.3, 4.04, 4.6, "eV", Tag::model_derived, "wilson_1966+gordon_chapman_1972",
                          "no measured phi of liquid Fe/Ni/Cr/steel (open); outgassed 304 SS thermionic phi 4.3 eV (Wilson, JAP 37 2261 (1966), abstract) carried into the melt: "
                          "molten refractory metals at T_m within ~0.1 V of solid values (Gordon & Chapman, NASA TN D-6888 (1972) Table I, FT), simple liquids 0-0.17 eV; "
                          "band = Fe thermionic spread 4.04-4.6 eV (Fomenko 1966 p. 20)");
  m.band_liquid_low_eV = 4.04;
  m.band_liquid_high_eV = 4.6;
  m.T_ref_K = 300;
  m.grid_eV = R.fixed("wf.grid", 0x1p-12, "eV", Tag::model_derived, THIS, "phi on the 2^-12 eV grid (one fast-spectrum configuration per value); 0.12 meV half-step << band");
  return m;
}

// ------------------------------------------------------------------ field-enhancement classes (scenario)
// No published beta or emitter-area distribution exists for broad-area stainless steel or rolled foil (open item).
// Scenario "ss304-foil-sites-v7": a flat class, a geometric-roughness class (hemisphere on a plane, beta = 3,
// Forbes, Edgcombe & Valdre, Ultramicroscopy 95 (2003) 57) and two field-emission-site classes in the beta range
// of FN plots of polished 304 cathodes (BastaniNejad et al., PRSTAB 15 083502 (2012) Tables II and IV: 214-1156)
// with area fractions of one emitter of 1e-16..1e-12 m^2 (Latham 1995, typical emitter areas, secondary) on a
// 3.1e-5 m^2 patch. Every value is a scenario parameter with its sensitivity band.
inline Enhancement ss304_enhancement(Registry& R) {
  Enhancement e;
  e.classes.clear();
  const double a3 = R.use("beta.roughness.area_fraction", 0.1, 0.0, 0.5, "1", Tag::scenario, "forbes_edgcombe_valdre_2003", "hemispherical asperities on the rolled finish; fraction not measured");
  const double b3 = R.use("beta.roughness.beta", 3.0, 1.0, 3.0, "1", Tag::model_derived, "forbes_edgcombe_valdre_2003", "hemisphere on a plane: beta = 3");
  const double r3 = R.use("beta.roughness.extent", 1e-6, 1e-7, 1e-5, "m", Tag::scenario, THIS, "asperity radius scale of a rolled foil");
  const double b100 = R.use("beta.site100.beta", 100.0, 50.0, 300.0, "1", Tag::scenario, "descoeudres_2009_calatroni_2010", "Cu after conditioning 77 +- 28 (LINAC2010 MOP070 Table 1); FN-plot betas are often spuriously large (Forbes et al. 2015)");
  const double a100 = R.use("beta.site100.area_fraction", 1e-9, 3e-12, 3e-8, "1", Tag::scenario, "latham_1995", "one emitter of 3e-14 m^2 on 3.1e-5 m^2 (band: 1e-16..1e-12 m^2 per emitter)");
  const double r100 = R.use("beta.site100.extent", 1e-7, 1e-8, 1e-6, "m", Tag::scenario, THIS, "tip radius ~ height/beta");
  const double b500 = R.use("beta.site500.beta", 500.0, 214.0, 1156.0, "1", Tag::measured, "bastaninejad_2012", "PRSTAB 15 083502 Tables II, IV: diamond-polished 304, 214-1156 (before/after Kr processing)");
  const double a500 = R.use("beta.site500.area_fraction", 1e-12, 3e-21, 3e-8, "1", Tag::scenario, "bastaninejad_2012", "FN-intercept single-emitter areas 1e-25..1e-22 m^2 (Table III) to 1e-12 m^2 (Latham) on 3.1e-5 m^2");
  const double r500 = R.use("beta.site500.extent", 2e-8, 1e-9, 1e-7, "m", Tag::scenario, THIS, "tip radius ~ height/beta");
  e.classes.push_back({1.0, 1.0 - a3 - a100 - a500, INFINITY, "flat"});
  e.classes.push_back({b3, a3, r3, "roughness-hemisphere"});
  e.classes.push_back({b100, a100, r100, "fn-site-100"});
  e.classes.push_back({b500, a500, r500, "fn-site-500"});
  validate(e);
  return e;
}

// ------------------------------------------------------------------ assembled surface options
inline SurfaceOptions ss304_surface(Registry& R, bool with_enhancement = true, bool with_wall = true) {
  SurfaceOptions s;
  s.work_function = ss304_work_function(R);
  if (with_enhancement) s.enhancement = ss304_enhancement(R);
  if (with_wall) s.wall = ss304_wall_model(R);
  s.photo = photo_model(R);
  return s;
}
// Richardson factor (OwnerOptions::richardson_lambda): the unified integral uses the ideal free-electron A0;
// measured effective constants A_R of the alloy elements are below A0 and absorb the temperature coefficient of phi
// (A_R = lambda A0 exp(-alpha/k)), so a lambda < 1 must be paired with the matching phi (scenario only).
inline double richardson_lambda(Registry& R) {
  return R.use("emission.richardson_lambda", 1.0, 0.0125, 1.0, "1", Tag::scenario, "fomenko_1966",
               "1 = ideal A0 = 120.17 A cm^-2 K^-2 (model_derived, CODATA); band: A_R/A0 of polycrystalline Fe (beta 26, gamma 1.5), Ni 30, Co 41, Cr 48 A cm^-2 K^-2 "
               "(Wahlin, Phys. Rev. 61 509 (1942) / 73 1458 (1948) via Fomenko pp. 17-22)");
}
inline Registry full_registry() {
  Registry R;
  (void)ss304_surface(R);
  (void)richardson_lambda(R);
  R.fixed("emission.T_max", 5000, "K", Tag::model_derived, THIS, "extended kernel domain (API-v7.txt 10.2); the phase owner's domain 298.2-5000 K");
  R.fixed("water.law_T_max", 3000, "K", Tag::recommended, "ss304_phase_owner_v1", "INTEGRATION.txt 2.2: water desorption law capped at 3000 K; stock above: typed refusal");
  R.fixed("free_electron.E_F_Fe", 11.1, "eV", Tag::recommended, "ashcroft_mermin_1976", "Table 2.1 (via HyperPhysics): free-electron Fermi energy of Fe; kT/E_F = 0.039 at 5000 K");
  return R;
}
}  // namespace surface_v7
