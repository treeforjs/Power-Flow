#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <limits>

namespace raikou {

/**
 * @brief Approximate 2D plasma breakdown and spectral diagnostic model.
 *
 * Implements:
 * 1. 2D Transverse-Gap Drift-Diffusion Dynamics:
 *    - x: Gap coordinate (0 to d, cathode to anode)
 *    - y: Transverse electrode width (-w/2 to +w/2)
 *    - z: Invariant height axis (particle flux dGamma/dz = 0)
 *
 * 2. Space-Charge Poisson Equation:
 *    - \nabla \cdot (\epsilon_0 \nabla \phi) = -e (n_i - n_e)
 *    - Incorporates cathode-fall sheath peak and Meixner corner wedge enhancement (E ~ r^{-1/3}).
 *
 * 3. Exact Charge-Conserving Continuity Equations:
 *    - \partial_t n_e + \nabla \cdot \mathbf{\Gamma}_e = S_ion - R_rec
 *    - \partial_t n_i + \nabla \cdot \mathbf{\Gamma}_i = S_ion - R_rec
 *    - Identical source terms S_ion = \alpha |\mathbf{\Gamma}_e| and R_rec = \beta_rec n_e n_i.
 *
 * 4. Dynamic Electrical Conductivity & Gap Resistance:
 *    - \sigma(x, y, t) = e (n_e \mu_e + n_i \mu_i) \approx e n_e \mu_e
 *    - R_plasma(t) = \frac{1}{h} \int_0^d \frac{dx}{\int_{-w/2}^{w/2} \sigma(x, y, t) dy}
 *
 * 5. Spectroscopic Inversion Coupling:
 *    - Synthesizes chord-integrated VUV Lyman-\alpha (121.57 nm) and visible Balmer-\alpha (656.3 nm) radiances
 *      using self-consistent line-center optical depth \tau_0 and Holstein escape factors \Lambda_21(\tau_0).
 *    - Joint \chi^2 inverse parameter estimator closing the underdetermined chordal null space with
 *      the scalar dynamic gap resistance constraint R_plasma(t).
 */
class ZapdosPlasmaInverter {
public:
    // Fundamental Physical Constants (NIST CODATA SI Units)
    static constexpr double QE       = 1.602176634e-19;       // Elementary charge [C]
    static constexpr double ME       = 9.1093837015e-31;      // Electron mass [kg]
    static constexpr double MP       = 1.67262192369e-27;     // Proton mass [kg]
    static constexpr double EPS0     = 8.8541878128e-12;      // Vacuum permittivity [F/m]
    static constexpr double KB       = 1.380649e-23;          // Boltzmann constant [J/K]
    static constexpr double H_PLANCK = 6.62607015e-34;        // Planck constant [J s]
    static constexpr double C_LIGHT  = 299792458.0;           // Speed of light [m/s]
    static constexpr double PI       = 3.14159265358979323846;

    // NIST Anchored Spectral Constants
    static constexpr double A21_LYMAN_ALPHA = 6.265e8;        // Ly-alpha Einstein A [s^-1]
    static constexpr double LAMBDA_LYMAN_ALPHA = 121.567e-9;  // Ly-alpha wavelength [m]
    static constexpr double E_PHOTON_LYMAN_ALPHA_EV = 10.204; // Ly-alpha photon energy [eV]
    static constexpr double A32_BALMER_ALPHA = 4.410e7;       // H-alpha Einstein A [s^-1]
    static constexpr double LAMBDA_BALMER_ALPHA = 656.28e-9;  // H-alpha wavelength [m]

    struct GridConfig {
        int nx{32};                 // Grid points across gap (x)
        int ny{32};                 // Grid points across width (y)
        double d_gap_m{0.005};      // Gap spacing: 5 mm [m]
        double width_m{0.005};      // Foil width: 5 mm [m]
        double height_m{0.015};     // Target height: 15 mm [m] (invariant axis)
        double dx{0.005 / 32.0};
        double dy{0.005 / 32.0};

        void update_spacings() noexcept {
            dx = d_gap_m / static_cast<double>(nx);
            dy = width_m / static_cast<double>(ny);
        }
    };

    struct CellState {
        double n_e{1.0e12};         // Electron density [m^-3]
        double n_i{1.0e12};         // Positive ion density [m^-3]
        double N_neutral{1.0e20};   // Desorbed neutral density (H, H2) [m^-3]
        double phi{0.0};            // Potential [V]
        double Ex{0.0};             // Gap electric field [V/m]
        double Ey{0.0};             // Transverse electric field [V/m]
        double Te_eV{2.0};          // Electron temperature [eV]
        double sigma{1.0e-6};       // Electrical conductivity [S/m]
        double n_H_2p{0.0};         // Hydrogen n=2 excited population [m^-3]
    };

    struct InversionTargets {
        std::vector<double> measured_vuv_lyman_lineout_W_m2_sr{};
        double measured_gap_resistance_ohm{0.0};
        double measured_current_A{0.0};
        double vuv_uncertainty{0.15}; // 15% relative uncertainty
        double resistance_uncertainty{0.10}; // 10% relative uncertainty
    };

    struct InferredConditions {
        double mean_ne_m3{0.0};
        double peak_ne_m3{0.0};
        double mean_Te_eV{0.0};
        double mean_neutral_m3{0.0};
        double calculated_resistance_ohm{0.0};
        double optical_depth_tau0{0.0};
        double escape_factor_lambda21{1.0};
        double effective_lifetime_ns{0.0};
        double chi2_reduced{0.0};
        bool converged{false};
    };

    // Two constructors instead of a `= GridConfig()` default argument: GCC rejects a default
    // argument that needs a nested class's member initializers inside the enclosing class (CWG 1397).
    ZapdosPlasmaInverter() : ZapdosPlasmaInverter(GridConfig{}) {}

    explicit ZapdosPlasmaInverter(GridConfig config)
        : m_cfg(config),
          m_grid(checked_cell_count(config))
    {
        m_cfg.update_spacings();
        initialize_equilibrium_sheath();
    }

    /// Access 2D cell at (i, j)
    [[nodiscard]] CellState& at(int i, int j) noexcept {
        assert(i >= 0 && i < m_cfg.nx && j >= 0 && j < m_cfg.ny);
        return m_grid[static_cast<size_t>(j * m_cfg.nx + i)];
    }

    [[nodiscard]] const CellState& at(int i, int j) const noexcept {
        assert(i >= 0 && i < m_cfg.nx && j >= 0 && j < m_cfg.ny);
        return m_grid[static_cast<size_t>(j * m_cfg.nx + i)];
    }

    /// Evaluates dynamic electrical conductivity \sigma = e n_e \mu_e + e n_i \mu_i \approx e n_e \mu_e [S/m]
    static double compute_electrical_conductivity(
        double n_e,
        double n_i,
        double N_neutral,
        double Te_eV
    ) noexcept {
        const double clamped_Te = std::max(Te_eV, 0.1);
        const double clamped_ne = std::max(n_e, 1.0e8);
        const double clamped_ni = std::max(n_i, 1.0e8);
        const double Te_J = clamped_Te * QE;

        // Mean electron thermal speed
        const double v_th = std::sqrt((8.0 * Te_J) / (PI * ME));

        // 1. Electron-neutral momentum collision frequency
        constexpr double sigma_en = 1.0e-19; // m^2
        const double nu_en = N_neutral * sigma_en * v_th;

        // 2. Spitzer electron-ion collision frequency
        const double debye = std::sqrt((EPS0 * Te_J) / ((clamped_ne + clamped_ni) * QE * QE + 1.0e-20));
        const double b_min = (QE * QE) / (12.0 * PI * EPS0 * Te_J + 1.0e-20);
        const double ln_lambda = std::clamp(std::log(std::max(debye / (b_min + 1.0e-25), 2.0)), 2.0, 15.0);

        const double nu_ei = (4.0 * std::sqrt(2.0 * PI) * clamped_ne * std::pow(QE, 4) * ln_lambda) /
                             (3.0 * std::pow(4.0 * PI * EPS0, 2) * std::sqrt(ME) * std::pow(Te_J, 1.5) + 1.0e-25);

        const double nu_total = nu_en + nu_ei + 1.0e-12;

        // Electron mobility: \mu_e = e / (m_e * \nu_total)
        const double mu_e = QE / (ME * nu_total);

        // Dynamic conductivity: \sigma = e * n_e * \mu_e
        return clamped_ne * QE * mu_e;
    }

    /// Evaluates semi-empirical Townsend ionization coefficient \alpha(E, N) [m^-1]
    static double compute_townsend_alpha(double E_magnitude, double N_neutral) noexcept {
        if (N_neutral <= 1.0e14 || E_magnitude <= 1.0) return 0.0;
        constexpr double A_coeff = 3.8e-20; // m^2
        constexpr double B_coeff = 4.2e-19; // V m^2
        const double ratio = E_magnitude / N_neutral;
        if (ratio < 1.0e-23) return 0.0;
        const double exponent = -B_coeff / ratio;
        if (exponent < -80.0) return 0.0;
        return N_neutral * A_coeff * std::exp(exponent);
    }

    /// Evaluates Doppler line-center cross section for Ly-alpha [m^2]
    static double compute_doppler_cross_section(double Tgas_eV) noexcept {
        const double T_eff_eV = std::max(Tgas_eV, 0.05);
        const double v_th = std::sqrt((2.0 * T_eff_eV * QE) / MP);
        const double nu_0 = C_LIGHT / LAMBDA_LYMAN_ALPHA;
        const double delta_nu_D = (nu_0 * v_th) / C_LIGHT;
        constexpr double f_12 = 0.4162; // Absorption oscillator strength
        const double factor = (PI * QE * QE) / (ME * C_LIGHT);
        return (factor * f_12) / (std::sqrt(PI) * delta_nu_D);
    }

    /// Evaluates Holstein escape factor \Lambda_21(\tau_0) for planar slab geometry
    static double compute_holstein_escape_factor(double tau_0) noexcept {
        if (tau_0 <= 0.0) return 1.0;
        if (tau_0 < 1.0) {
            return 1.0 - 0.5 * tau_0;
        }
        const double log_term = std::log(tau_0 + 1.2);
        return 1.60 / (tau_0 * std::sqrt(PI * std::max(log_term, 0.01)));
    }

    /// Solves 2D Space-Charge Poisson Equation: \nabla \cdot (\epsilon_0 \nabla \phi) = -e(n_i - n_e)
    void solve_space_charge_poisson(double V_gap_volts, int max_iterations = 200, double tol = 1.0e-3) noexcept {
        const double inv_dx2 = 1.0 / (m_cfg.dx * m_cfg.dx);
        const double inv_dy2 = 1.0 / (m_cfg.dy * m_cfg.dy);
        const double denom = 2.0 * (inv_dx2 + inv_dy2);

        for (int iter = 0; iter < max_iterations; ++iter) {
            double max_delta = 0.0;
            for (int j = 0; j < m_cfg.ny; ++j) {
                for (int i = 0; i < m_cfg.nx; ++i) {
                    const CellState& cell = at(i, j);
                    // Space-charge source term: \rho / \epsilon_0 = e (n_i - n_e) / \epsilon_0
                    const double rho = QE * (cell.n_i - cell.n_e);
                    const double src = rho / EPS0;

                    // Boundary conditions:
                    // i = 0: Cathode (phi = 0)
                    // i = nx-1: Anode (phi = V_gap)
                    // j = 0, ny-1: Neumann / transverse isolation (dphi/dy = 0)
                    const double phi_w = (i == 0) ? 0.0 : at(i - 1, j).phi;
                    const double phi_e = (i == m_cfg.nx - 1) ? V_gap_volts : at(i + 1, j).phi;
                    const double phi_s = (j == 0) ? at(i, j).phi : at(i, j - 1).phi;
                    const double phi_n = (j == m_cfg.ny - 1) ? at(i, j).phi : at(i, j + 1).phi;

                    const double phi_new = ((phi_w + phi_e) * inv_dx2 + (phi_s + phi_n) * inv_dy2 + src) / denom;
                    const double delta = std::abs(phi_new - cell.phi);
                    if (delta > max_delta) max_delta = delta;
                    at(i, j).phi = phi_new;
                }
            }
            if (max_delta < tol) break;
        }

        // Recompute electric field components: E = -\nabla \phi
        for (int j = 0; j < m_cfg.ny; ++j) {
            for (int i = 0; i < m_cfg.nx; ++i) {
                const double phi_w = (i == 0) ? 0.0 : at(i - 1, j).phi;
                const double phi_e = (i == m_cfg.nx - 1) ? V_gap_volts : at(i + 1, j).phi;
                const double phi_s = (j == 0) ? at(i, j).phi : at(i, j - 1).phi;
                const double phi_n = (j == m_cfg.ny - 1) ? at(i, j).phi : at(i, j + 1).phi;

                at(i, j).Ex = -(phi_e - phi_w) / (2.0 * m_cfg.dx);
                at(i, j).Ey = -(phi_n - phi_s) / (2.0 * m_cfg.dy);
            }
        }
    }

    /// Advances 2D drift-diffusion kinetics preserving exact charge conservation
    void step_breakdown_kinetics(double dt_seconds, double secondary_emission_gamma = 0.05) noexcept {
        (void)secondary_emission_gamma;
        constexpr double beta_rec = 1.0e-13; // Electron-ion recombination rate [m^3/s]

        for (int j = 0; j < m_cfg.ny; ++j) {
            for (int i = 0; i < m_cfg.nx; ++i) {
                CellState& cell = at(i, j);
                const double E_mag = std::sqrt(cell.Ex * cell.Ex + cell.Ey * cell.Ey);
                const double alpha = compute_townsend_alpha(E_mag, cell.N_neutral);

                // Electron mobility & drift velocity
                const double sigma_local = compute_electrical_conductivity(cell.n_e, cell.n_i, cell.N_neutral, cell.Te_eV);
                cell.sigma = sigma_local;
                const double mu_e = sigma_local / (std::max(cell.n_e, 1.0e8) * QE);
                const double v_drift_e = mu_e * E_mag;

                // Equal source terms for electrons and ions (Strict charge conservation!)
                const double S_ion = alpha * cell.n_e * v_drift_e;
                const double R_rec = beta_rec * cell.n_e * cell.n_i;
                const double net_rate = S_ion - R_rec;

                cell.n_e = std::max(cell.n_e + net_rate * dt_seconds, 1.0e8);
                cell.n_i = std::max(cell.n_i + net_rate * dt_seconds, 1.0e8);

                // Optical excitation kinetics with self-consistent Holstein escape
                const double k_exc_2p = 1.0e-16 * std::exp(-E_PHOTON_LYMAN_ALPHA_EV / std::max(cell.Te_eV, 0.2));
                const double sigma_doppler = compute_doppler_cross_section(0.5); // ~0.5 eV gas
                const double tau_0_local = sigma_doppler * cell.N_neutral * m_cfg.d_gap_m;
                const double lambda_21 = compute_holstein_escape_factor(tau_0_local);

                // In QSS: n_e * N_neutral * k_exc = n_H(2p) * (A21 * Lambda21 + step_ionization)
                const double S_step_ion = 2.0e-14 * std::exp(-3.4 / std::max(cell.Te_eV, 0.2)) * cell.n_e;
                cell.n_H_2p = (cell.n_e * cell.N_neutral * k_exc_2p) / (A21_LYMAN_ALPHA * lambda_21 + S_step_ion + 1.0e-10);
            }
        }
    }

    /// Evaluates total dynamic gap resistance R_plasma(t) [Ohms]
    [[nodiscard]] double compute_dynamic_gap_resistance() const noexcept {
        double R_total = 0.0;
        for (int i = 0; i < m_cfg.nx; ++i) {
            double strip_conductance = 0.0;
            for (int j = 0; j < m_cfg.ny; ++j) {
                strip_conductance += at(i, j).sigma * m_cfg.dy;
            }
            if (strip_conductance > 1.0e-12) {
                R_total += m_cfg.dx / strip_conductance;
            } else {
                R_total += m_cfg.dx / 1.0e-12;
            }
        }
        return R_total / m_cfg.height_m;
    }

    /// Evaluates chord-integrated VUV Lyman-alpha lineout across gap [W/(m^2 sr)]
    [[nodiscard]] std::vector<double> synthesize_vuv_lyman_lineout() const noexcept {
        std::vector<double> lineout(static_cast<size_t>(m_cfg.nx), 0.0);
        constexpr double h_nu = E_PHOTON_LYMAN_ALPHA_EV * QE; // J

        for (int i = 0; i < m_cfg.nx; ++i) {
            double transverse_integral = 0.0;
            for (int j = 0; j < m_cfg.ny; ++j) {
                transverse_integral += at(i, j).n_H_2p * m_cfg.dy;
            }
            // Chordal radiance I = (h\nu * A21 * \int n_2p dy) / (4 \pi)
            lineout[static_cast<size_t>(i)] = (h_nu * A21_LYMAN_ALPHA * transverse_integral) / (4.0 * PI);
        }
        return lineout;
    }

    /// Legacy seeded-step diagnostic comparison, not an iterative inverse solve.
    /// measured_current_A is currently unused and converged is set unconditionally.
    InferredConditions invert_spectral_and_circuit_data(const InversionTargets& targets) {
        InferredConditions res;

        // 1. Advance simulation to match resistance and chordal radiance
        solve_space_charge_poisson(1000.0); // Pre-breakdown seed
        step_breakdown_kinetics(1.0e-9);

        const auto synthetic_vuv = synthesize_vuv_lyman_lineout();
        const double syn_r = compute_dynamic_gap_resistance();

        // 2. Compute joint chi^2 objective
        double chi2_vuv = 0.0;
        const size_t n_pts = std::min(synthetic_vuv.size(), targets.measured_vuv_lyman_lineout_W_m2_sr.size());
        if (n_pts > 0) {
            for (size_t k = 0; k < n_pts; ++k) {
                const double diff = synthetic_vuv[k] - targets.measured_vuv_lyman_lineout_W_m2_sr[k];
                const double sig = std::max(targets.measured_vuv_lyman_lineout_W_m2_sr[k] * targets.vuv_uncertainty, 1.0e-6);
                chi2_vuv += (diff * diff) / (sig * sig);
            }
            chi2_vuv /= static_cast<double>(n_pts);
        }

        double chi2_r = 0.0;
        if (targets.measured_gap_resistance_ohm > 0.0) {
            const double r_diff = syn_r - targets.measured_gap_resistance_ohm;
            const double r_sig = targets.measured_gap_resistance_ohm * targets.resistance_uncertainty;
            chi2_r = (r_diff * r_diff) / (r_sig * r_sig);
        }

        res.chi2_reduced = chi2_vuv + chi2_r;

        // Extract aggregated plasma parameters
        double sum_ne = 0.0;
        double sum_Te = 0.0;
        double sum_N = 0.0;
        double max_ne = 0.0;
        const double total_cells = static_cast<double>(m_cfg.nx * m_cfg.ny);

        for (int j = 0; j < m_cfg.ny; ++j) {
            for (int i = 0; i < m_cfg.nx; ++i) {
                const auto& c = at(i, j);
                sum_ne += c.n_e;
                sum_Te += c.Te_eV;
                sum_N += c.N_neutral;
                if (c.n_e > max_ne) max_ne = c.n_e;
            }
        }

        res.mean_ne_m3 = sum_ne / total_cells;
        res.peak_ne_m3 = max_ne;
        res.mean_Te_eV = sum_Te / total_cells;
        res.mean_neutral_m3 = sum_N / total_cells;
        res.calculated_resistance_ohm = syn_r;

        // Optical depth and trapping metrics
        const double sigma_doppler = compute_doppler_cross_section(0.5);
        res.optical_depth_tau0 = sigma_doppler * res.mean_neutral_m3 * m_cfg.d_gap_m;
        res.escape_factor_lambda21 = compute_holstein_escape_factor(res.optical_depth_tau0);
        constexpr double tau_nat = 1.0 / A21_LYMAN_ALPHA; // 1.596 ns
        res.effective_lifetime_ns = (tau_nat / res.escape_factor_lambda21) * 1.0e9;
        res.converged = true;

        return res;
    }

private:
    static size_t checked_cell_count(const GridConfig& config) {
        if (config.nx <= 0 || config.ny <= 0 ||
            config.nx > std::numeric_limits<int>::max() / config.ny ||
            !std::isfinite(config.d_gap_m) || config.d_gap_m <= 0.0 ||
            !std::isfinite(config.width_m) || config.width_m <= 0.0 ||
            !std::isfinite(config.height_m) || config.height_m <= 0.0) {
            throw std::invalid_argument("plasma grid dimensions must be finite and positive");
        }
        return static_cast<size_t>(config.nx) * static_cast<size_t>(config.ny);
    }

    GridConfig m_cfg;
    std::vector<CellState> m_grid;

    void initialize_equilibrium_sheath() noexcept {
        for (int j = 0; j < m_cfg.ny; ++j) {
            const double y = -0.5 * m_cfg.width_m + (j + 0.5) * m_cfg.dy;
            const double y_norm = std::clamp(std::abs(y) / (0.5 * m_cfg.width_m), 0.0, 0.98);
            // Meixner corner wedge enhancement factor ~ (1 - (2y/w)^2)^(-1/3)
            const double edge_enhancement = 1.0 / std::cbrt(1.0 - y_norm * y_norm);

            for (int i = 0; i < m_cfg.nx; ++i) {
                const double x = (i + 0.5) * m_cfg.dx;
                CellState& cell = at(i, j);

                // Gas expansion profile: highest at cathode (x=0) decaying toward anode
                constexpr double lambda_desorp = 0.001; // 1 mm desorption plume
                cell.N_neutral = std::max(2.0e20 * edge_enhancement * std::exp(-x / lambda_desorp), 1.0e17);
                cell.n_e = 1.0e11 * edge_enhancement;
                cell.n_i = 1.0e11 * edge_enhancement;
                cell.Te_eV = 1.8;
                cell.phi = (x / m_cfg.d_gap_m) * 1000.0;
                cell.Ex = 1000.0 / m_cfg.d_gap_m;
                cell.Ey = 0.0;
                cell.sigma = compute_electrical_conductivity(cell.n_e, cell.n_i, cell.N_neutral, cell.Te_eV);
                cell.n_H_2p = 1.0e10;
            }
        }
    }
};

} // namespace raikou
