#pragma once
#include "raikou/spherical_harmonics.h"
#include "raikou/distributed.h"
#include <limits>
#include <optional>

namespace raikou::kinetics {
enum class KineticWall { Absorbing, Specular, PrescribedIncoming };
struct KineticFace {
    std::size_t left{0};
    std::size_t right{std::numeric_limits<std::size_t>::max()};
    Vector3 normal{1, 0, 0}; // unit normal from left to right/outside
    double area_m2{1};
    KineticWall wall{KineticWall::Absorbing};
    std::optional<SphericalEedfState> incoming{};
};
struct KineticMesh {
    std::vector<double> cell_volumes_m3;
    // Each physical face appears exactly once. Periodic faces connect their two cells.
    std::vector<KineticFace> faces;
};
struct KineticStreamingResult {
    // Native staggered energy/angular layout, one rate vector per spatial cell.
    std::vector<std::vector<double>> density_mode_rates;
    std::vector<double> face_particle_flux_m2_s;
    std::vector<double> face_energy_flux_ev_m2_s;
    double outward_particle_rate_s{0};
    double outward_energy_rate_ev_s{0};
    double maximum_explicit_step_s{0};
};
// Conservative P_N upwind flux in Cartesian space, retaining all real m modes. Even/odd
// energy grids use a shared overlap operator and its mass adjoint; the overlap preserves
// both discrete number-current and energy-current moments. No drift-diffusion closure is
// substituted for distribution transport. Suitable as an implicit spatial residual term.
[[nodiscard]] KineticStreamingResult kinetic_streaming_rhs(const KineticMesh& mesh,
    const std::vector<SphericalEedfState>& states, std::size_t threads = 1);
struct KineticStreamingOperator {
    backend::DistributedCsrMatrix matrix;
    std::vector<double> affine_rates; // prescribed incoming boundary contribution
};
// Exact row-owned assembly of the same conservative streaming operator. Static
// geometry may replicate; only owned matrix rows and boundary source are stored.
[[nodiscard]] KineticStreamingOperator assemble_kinetic_streaming(const KineticMesh& mesh,
    const EedfGrid& grid, backend::RowPartition mode_rows);
struct KineticStreamingStep {
    std::vector<SphericalEedfState> states;
    double outward_particles{0};
    double outward_energy_ev{0};
};
// Transactional explicit Euler reference/transport step with a conservative CFL bound.
// Implicit coupled solvers call kinetic_streaming_rhs directly at their candidate state.
[[nodiscard]] KineticStreamingStep advance_kinetic_streaming(const KineticMesh& mesh,
    const std::vector<SphericalEedfState>& states, double dt_s, std::size_t threads = 1);
} // namespace raikou::kinetics
