#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <algorithm>
#include <stdexcept>

namespace raikou {

/**
 * @brief Approximate VUV spectrometer response and branching-ratio calibration model.
 *
 * Models the McPherson 1-meter normal-incidence VUV spectrometer and Photek iCMOS diagnostic
 * deployed on the Mykonos 850 kA pulsed-power accelerator.
 *
 * Features:
 *   1. Optical collection geometry & throughput: solid angle Omega, entrance slit area,
 *      magnification, and limiting etendue G = A_slit * Omega.
 *   2. Component-level efficiency modeling:
 *      - MgF2 vacuum window cutoff & transmission T_MgF2(lambda) (115 - 200 nm).
 *      - Al + MgF2 turning mirror reflectances R_mirrors(lambda).
 *      - Grating diffraction efficiency eta_grating(lambda).
 *      - Photek iCMOS MCP photocathode quantum efficiency QE(lambda), MCP electron gain,
 *        phosphor screen conversion, fiber optic relay, and ADC counts-per-photon conversion.
 *   3. Theoretical relative radiometric calibration via Molecular Nitrogen (N2)
 *      Lyman-Birge-Hopfield (LBH, a^1Pi_g -> X^1Sigma_g+) branching ratios across 120 - 180 nm.
 *   4. Direct two-way conversion between raw detector counts S(lambda) and physical
 *      plasma emission radiance I(lambda) [W/(m^2 sr)].
 */
class SpectrometerCalibration {
public:
    static constexpr double H_PLANCK = 6.62607015e-34;  // J s
    static constexpr double C_LIGHT = 299792458.0;      // m/s
    static constexpr double PI = 3.14159265358979323846;

    struct CollectionGeometry {
        double mirror_diameter_m{0.050};   // 50 mm concave turning mirror
        double source_distance_m{0.350};   // 350 mm distance from plasma AK gap to mirror
        double slit_width_m{25.0e-6};       // 25 um spectrometer entrance slit
        double slit_height_m{3.0e-3};      // 3.0 mm slit height
        double optical_magnification{1.0}; // 1:1 relay imaging

        [[nodiscard]] double solid_angle_sr() const noexcept {
            const double r = mirror_diameter_m / 2.0;
            return (PI * r * r) / (source_distance_m * source_distance_m);
        }

        [[nodiscard]] double slit_area_m2() const noexcept {
            return slit_width_m * slit_height_m;
        }

        [[nodiscard]] double etendue_m2_sr() const noexcept {
            return slit_area_m2() * solid_angle_sr();
        }
    };

    struct DetectorConfig {
        double mcp_gain{2500.0};                        // MCP electron multiplication gain
        double phosphor_photons_per_electron{30.0};    // P43 phosphor screen conversion yield
        double fiber_coupling_efficiency{0.25};        // Fiber optic taper transmission
        double cmos_quantum_efficiency{0.55};          // CMOS green (530 nm) sensor QE
        double adc_counts_per_photoelectron{0.75};     // Camera digital ADC sensitivity factor
    };

    struct N2LBHTransition {
        int v_upper{0};             // Upper vibrational state v' in a^1Pi_g
        int v_lower{0};             // Lower vibrational state v'' in X^1Sigma_g+
        double wavelength_nm{0.0};  // Transition wavelength (nm)
        double einstein_A_s{0.0};   // Spontaneous decay rate A(v', v'') (s^-1)
        double franck_condon{0.0};  // Franck-Condon factor q(v', v'')
    };

    struct MeasuredBand {
        int v_upper{0};
        int v_lower{0};
        double wavelength_nm{0.0};
        double raw_counts{0.0};
        double counts_sigma{0.0};
    };

    struct RelativeCalibrationPoint {
        double wavelength_nm{0.0};
        double relative_efficiency{1.0}; // Pairwise efficiency ratio to the reference band
        double uncertainty{0.0};
    };

    SpectrometerCalibration(const CollectionGeometry& geom, const DetectorConfig& det);
    /// Default detector configuration. (Default arguments of nested aggregates with member
    /// initializers do not compile on GCC/Clang before the class is complete: CWG 1397.)
    explicit SpectrometerCalibration(const CollectionGeometry& geom);
    SpectrometerCalibration();

    // Component transmission functions
    [[nodiscard]] static double transmission_mgf2(double wavelength_nm) noexcept;
    [[nodiscard]] static double reflectance_mirrors(double wavelength_nm, int num_mirrors = 2) noexcept;
    [[nodiscard]] static double efficiency_grating(double wavelength_nm) noexcept;
    [[nodiscard]] static double photocathode_qe(double wavelength_nm) noexcept;

    /// Master absolute system efficiency E(lambda) [counts per incident photon at slit]
    [[nodiscard]] double master_efficiency(double wavelength_nm) const noexcept;

    /// Converts raw detector counts into physical emission radiance I(lambda) [W/(m^2 sr)]
    [[nodiscard]] double counts_to_radiance(double raw_counts, double wavelength_nm, double gate_time_s) const;

    /// Converts physical emission radiance into expected raw camera counts
    [[nodiscard]] double radiance_to_counts(double radiance_W_m2_sr, double wavelength_nm, double gate_time_s) const;

    /// Reference N2 Lyman-Birge-Hopfield (LBH) branching ratio table
    [[nodiscard]] const std::vector<N2LBHTransition>& get_n2_lbh_database() const noexcept {
        return n2_lbh_table_;
    }

    /// Theoretical branching ratio between two transitions originating from the same v'
    [[nodiscard]] double theoretical_branching_ratio(int v_prime, int v_double_prime_1, int v_double_prime_2) const;

    /// Evaluates relative spectrometer calibration R(lambda) from in-situ N2 band measurements
    [[nodiscard]] std::vector<RelativeCalibrationPoint> evaluate_branching_calibration(
        const std::vector<MeasuredBand>& measured_bands
    ) const;

private:
    CollectionGeometry geom_;
    DetectorConfig det_;
    std::vector<N2LBHTransition> n2_lbh_table_;

    void init_n2_lbh_table();
};

inline SpectrometerCalibration::SpectrometerCalibration(const CollectionGeometry& geom)
    : SpectrometerCalibration(geom, DetectorConfig{}) {}
inline SpectrometerCalibration::SpectrometerCalibration()
    : SpectrometerCalibration(CollectionGeometry{}, DetectorConfig{}) {}

} // namespace raikou
