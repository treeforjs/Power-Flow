#pragma once
#include "raikou/electron_boltzmann.h"
#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace raikou::kinetics {
using Vector3 = std::array<double, 3>;
// Active, proper Cartesian rotation in row-major order.
struct AngularRotation {
    std::array<double, 9> matrix{1, 0, 0, 0, 1, 0, 0, 0, 1};
};
[[nodiscard]] AngularRotation axis_angle_rotation(const Vector3& axis, double angle_rad);
[[nodiscard]] AngularRotation rotation_to_axis(const Vector3& axis);
[[nodiscard]] AngularRotation inverse_rotation(const AngularRotation& rotation);
[[nodiscard]] AngularRotation compose_rotations(const AngularRotation& left, const AngularRotation& right);
[[nodiscard]] Vector3 rotate_vector(const AngularRotation& rotation, const Vector3& vector);
[[nodiscard]] std::size_t spherical_index(int degree, int order);
// Real basis orthonormal under the spherical average: R00=1, R11=sqrt(3)x,
// R1,-1=sqrt(3)y, R10=sqrt(3)z. Negative orders are sine, positive orders cosine.
// No Condon-Shortley phase. Coefficients therefore have the units of F0.
[[nodiscard]] std::vector<double> real_spherical_harmonics(int terms, const Vector3& direction);
struct AngularNode {
    Vector3 direction;
    double weight{0}; // weights sum to one (spherical average)
    std::vector<double> basis;
};
[[nodiscard]] std::vector<AngularNode> spherical_quadrature(int terms, int polar_points = 0, int azimuthal_points = 0);
class HarmonicRotation {
public:
    HarmonicRotation(int terms, const AngularRotation& rotation);
    void apply(std::span<const double> input_modes, std::span<double> output_modes) const;
    [[nodiscard]] int terms() const noexcept { return terms_; }
    [[nodiscard]] const std::vector<std::size_t>& offsets() const noexcept { return offsets_; }
    [[nodiscard]] const std::vector<double>& blocks() const noexcept { return blocks_; }
private:
    int terms_{0};
    std::vector<std::size_t> offsets_;
    std::vector<double> blocks_;
};

// Full angular state in a fixed laboratory basis. Per energy slot, all terms^2 modes
// use index l*l+l+m, -l<=m<=l. Even degrees live on energy cells, odd degrees at
// energy face i+1; the last odd-degree slots are zero boundary padding. Rotation mixes
// only modes of one degree and thus never interpolates between radial locations.
struct SphericalEedfState {
    EedfGrid grid;
    double time_s{0};
    std::vector<double> density_modes;
    void validate() const;
    [[nodiscard]] double electron_density_m3() const;
    [[nodiscard]] double electron_energy_ev_m3() const;
    [[nodiscard]] double mean_energy_ev() const;
    [[nodiscard]] Vector3 particle_flux_m2_s() const;
    [[nodiscard]] Vector3 energy_flux_ev_m2_s() const;
};
[[nodiscard]] SphericalEedfState spherical_from_legendre(const EedfState& state,
                                                       const Vector3& symmetry_axis = {0, 0, 1});
// Refuses to discard non-axisymmetric modes. Explicit angular projections must be
// performed by a caller choosing an approximation, not hidden in compatibility code.
[[nodiscard]] EedfState legendre_from_spherical(const SphericalEedfState& state,
                                               const Vector3& symmetry_axis = {0, 0, 1},
                                               double angular_tolerance = 1e-10);
[[nodiscard]] SphericalEedfState rotate_spherical_eedf(const SphericalEedfState& state,
                                                      const AngularRotation& rotation);
struct SphericalEedfStepResult {
    SphericalEedfState state;
    SwarmResult swarm;
    Vector3 flux_velocity_m_s{};
    int nonlinear_iterations{0};
    double linear_residual{0};
    double energy_balance_relative{0};
};
struct SphericalEedfEvaluation {
    std::vector<double> density_mode_rates;
    SwarmResult swarm;
    double electron_production_m3_s{0};
    double electron_energy_rate_ev_m3_s{0};
};
// Instantaneous semidiscrete collision/field RHS for implicit spatial/time residuals.
[[nodiscard]] SphericalEedfEvaluation evaluate_spherical_eedf(const std::vector<CollisionProcess>& processes,
    const EedfConditions& conditions, const SphericalEedfState& state, const Vector3& electric_field_V_m);
// Physical electric field: electron acceleration is along -E. Arbitrary time-varying
// 2D/3D field direction. All real azimuthal modes are
// retained while aligning the velocity operator to the instantaneous field axis.
[[nodiscard]] SphericalEedfStepResult advance_spherical_eedf(const std::vector<CollisionProcess>& processes,
    const EedfConditions& conditions, const SphericalEedfState& initial,
    const Vector3& electric_field_V_m, double dt_s);
[[nodiscard]] raikou::Json spherical_state_to_json(const SphericalEedfState&);
[[nodiscard]] SphericalEedfState spherical_state_from_json(const raikou::Json&);
[[nodiscard]] raikou::Json evaluate_spherical_eedf_json(const raikou::Json&);
[[nodiscard]] raikou::Json advance_spherical_eedf_json(const raikou::Json&);
} // namespace raikou::kinetics
