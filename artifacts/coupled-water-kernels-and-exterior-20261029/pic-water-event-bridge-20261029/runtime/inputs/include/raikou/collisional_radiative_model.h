#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <algorithm>
#include <stdexcept>

namespace raikou {

/**
 * @brief Approximate five-level hydrogen collisional-radiative model.
 *
 * Solves non-LTE excited state population kinetics (n = 1..5) in the Quasi-Steady-State (QSS)
 * regime, including:
 *   1. Electron-impact excitation C_ij(Te) and ionization S_i(Te) (Johnson 1972 / Drawin).
 *   2. Detailed balance de-excitation C_ji(Te).
 *   3. Einstein A_ij transition probabilities (NIST ASD).
 *   4. Radiative and three-body recombination (alpha_rad, alpha_3b).
 *   5. Molecular hydrogen (H2) dissociative excitation feeding.
 *   6. Self-consistent line-center optical depth tau_0 and Holstein escape factors Lambda_1p(tau_0)
 *      for resonant Lyman transitions (Ly-alpha 121.57 nm, Ly-beta 102.57 nm).
 *   7. Synthetic VUV & visible spectral radiance calculation.
 *   8. Objective chi^2 inverse parameter fitting engine.
 */
class CollisionalRadiativeModel {
public:
    static constexpr int NUM_LEVELS = 5; // Principal quantum numbers n = 1, 2, 3, 4, 5

    // Fundamental Physical Constants (SI Units)
    static constexpr double QE = 1.602176634e-19;       // Elementary charge (C)
    static constexpr double ME = 9.1093837015e-31;      // Electron mass (kg)
    static constexpr double H_PLANCK = 6.62607015e-34;  // Planck constant (J s)
    static constexpr double C_LIGHT = 299792458.0;      // Speed of light (m/s)
    static constexpr double KB = 1.380649e-23;          // Boltzmann constant (J/K)
    static constexpr double EPS0 = 8.8541878128e-12;    // Vacuum permittivity (F/m)
    static constexpr double A0 = 5.29177210903e-11;     // Bohr radius (m)
    static constexpr double IH_EV = 13.605693122994;    // Rydberg energy (eV)
    static constexpr double AMU = 1.66053906660e-27;    // Atomic mass unit (kg)
    static constexpr double M_H = 1.00784 * AMU;        // Hydrogen atom mass (kg)
    static constexpr double PI = 3.14159265358979323846;

    struct PlasmaState {
        double Te_eV{2.5};         // Electron temperature (eV)
        double ne_m3{1.0e19};       // Electron density (m^-3)
        double Tgas_K{5000.0};      // Neutral gas temperature (K)
        double nH_ground_m3{1.0e20};// Atomic hydrogen ground state n(1) (m^-3)
        double nH2_m3{5.0e19};      // Molecular hydrogen density (m^-3)
        double path_length_m{5.0e-3}; // Line-of-sight path length (m)
    };

    struct CRRadianceResult {
        std::array<double, NUM_LEVELS> populations_m3{}; // n(1) to n(5)
        double lyman_alpha_radiance_W_m2_sr{0.0};        // Ly-alpha (121.57 nm, 2p -> 1s)
        double lyman_beta_radiance_W_m2_sr{0.0};         // Ly-beta (102.57 nm, 3p -> 1s)
        double balmer_alpha_radiance_W_m2_sr{0.0};       // H-alpha (656.28 nm, 3 -> 2)
        double balmer_beta_radiance_W_m2_sr{0.0};        // H-beta (486.13 nm, 4 -> 2)
        double molecular_lyman_proxy_W_m2_sr{0.0};       // H2 Lyman band proxy (~160 nm)
        double tau_0_lyman_alpha{0.0};                   // Line-center Ly-alpha optical depth
        double escape_factor_lyman_alpha{1.0};           // Holstein factor Lambda_21
        double tau_0_lyman_beta{0.0};                    // Line-center Ly-beta optical depth
        double escape_factor_lyman_beta{1.0};            // Holstein factor Lambda_31
    };

    struct FitPriors {
        double min_Te_eV{0.5};
        double max_Te_eV{20.0};
        double min_ne_m3{1.0e17};
        double max_ne_m3{1.0e22};
        double min_nH_m3{1.0e17};
        double max_nH_m3{1.0e23};
        double min_nH2_m3{1.0e17};
        double max_nH2_m3{1.0e23};
    };

    struct FitResult {
        double Te_eV{0.0};
        double Te_err_eV{0.0};
        double ne_m3{0.0};
        double ne_err_m3{0.0};
        double nH_m3{0.0};
        double nH_err_m3{0.0};
        double nH2_m3{0.0};
        double nH2_err_m3{0.0};
        double tau_0_lyman_alpha{0.0};
        double escape_factor_lyman_alpha{1.0};
        double chi2_reduced{0.0};
        bool converged{false};
        int iterations{0};
    };

    CollisionalRadiativeModel();

    /// Evaluates excited populations and spectral radiances for a given plasma state.
    /// Temperatures must be positive; densities and path length nonnegative; all finite.
    /// Throws std::invalid_argument for invalid input, std::overflow_error on overflow.
    CRRadianceResult evaluate(const PlasmaState& state) const;

    /// Line-center Doppler optical depth and Holstein escape factor
    static double compute_lyman_alpha_cross_section(double Tgas_K);
    static double compute_lyman_beta_cross_section(double Tgas_K);
    static double compute_optical_depth(double ground_density_m3, double path_length_m, double cross_section_m2);
    static double compute_holstein_escape_factor(double tau_0);

    /// State-to-state rate coefficients
    static double excitation_rate(int i, int j, double Te_eV);    // i < j (m^3/s)
    static double deexcitation_rate(int j, int i, double Te_eV);  // j > i (m^3/s)
    static double ionization_rate(int n, double Te_eV);           // (m^3/s)
    static double radiative_recombination(int n, double Te_eV);   // (m^3/s)
    static double three_body_recombination(int n, double Te_eV);  // (m^6/s)
    static double h2_dissociation_to_level(int p, double Te_eV);  // (m^3/s)

    /// Regularized log-radiance fit against two VUV channels with four bounded parameters.
    /// The fit is underdetermined. Error fields are simplex-span heuristics, not
    /// statistical confidence intervals; chi2_reduced is the legacy objective / 2.
    FitResult fit_vuv_shot(
        double measured_ly_alpha_radiance,
        double measured_ly_alpha_sigma,
        double measured_molecular_radiance,
        double measured_molecular_sigma,
        double path_length_m,
        double fixed_Tgas_K,
        const FitPriors& priors
    ) const;

    /// Same as above with the default gas temperature (5000 K) and default priors. (A default
    /// argument of a nested aggregate with member initializers does not compile on GCC/Clang
    /// before the class is complete: CWG 1397.)
    FitResult fit_vuv_shot(
        double measured_ly_alpha_radiance,
        double measured_ly_alpha_sigma,
        double measured_molecular_radiance,
        double measured_molecular_sigma,
        double path_length_m,
        double fixed_Tgas_K = 5000.0
    ) const;

private:
    // Atomic level energies (eV) and statistical weights
    std::array<double, NUM_LEVELS> level_energies_eV_{};
    std::array<double, NUM_LEVELS> stat_weights_{};

    // Einstein spontaneous emission coefficients A_ji (s^-1) [upper][lower]
    std::array<std::array<double, NUM_LEVELS>, NUM_LEVELS> A_einstein_{};

    // Absorption oscillator strengths f_ij [lower][upper]
    std::array<std::array<double, NUM_LEVELS>, NUM_LEVELS> f_oscillator_{};

    void init_atomic_data();

    /// Solves linear system M * x = b for 4x4 excited populations (n=2..5)
    static void solve_4x4(const double M[4][4], const double b[4], double x[4]);
};

inline CollisionalRadiativeModel::FitResult CollisionalRadiativeModel::fit_vuv_shot(
    double measured_ly_alpha_radiance, double measured_ly_alpha_sigma, double measured_molecular_radiance,
    double measured_molecular_sigma, double path_length_m, double fixed_Tgas_K) const {
    return fit_vuv_shot(measured_ly_alpha_radiance, measured_ly_alpha_sigma, measured_molecular_radiance,
                        measured_molecular_sigma, path_length_m, fixed_Tgas_K, FitPriors{});
}

} // namespace raikou
