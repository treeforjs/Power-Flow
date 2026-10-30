#pragma once
#include "raikou/electron_coupling.h"
#include "raikou/kinetic_transport.h"
#include "raikou/distributed.h"
#include <span>

namespace raikou::spatial::mesh { struct Geometry; }
namespace raikou::kinetics {
struct SpatialKineticCell {
    // Null coupling means an absent electron phase. Its full state must be zero;
    // every active/absent contact must be supplied as an explicit kinetic wall.
    const chemistry::Mechanism* mechanism{nullptr};
    const ElectronCoupling* coupling{nullptr};
    chemistry::State candidate_state;
    chemistry::Environment environment;
    Vector3 electric_field_V_m{};
};
struct SpatialKineticBoundary {
    KineticWall kind{KineticWall::Absorbing};
    std::optional<SphericalEedfState> incoming{};
};
// Geometry is static and may be replicated. Phase volumes and contact apertures
// are authored physical data. Each active/absent interior contact requires a wall
// entry keyed by Geometry::faces index; walls follow Geometry::walls order.
[[nodiscard]] KineticMesh make_kinetic_mesh(const spatial::mesh::Geometry& geometry,
    std::span<const double> gas_cell_volumes_m3,
    std::span<const double> gas_interior_areas_m2,
    std::span<const SpatialKineticBoundary> walls,
    const std::map<std::size_t, SpatialKineticBoundary>& phase_contacts = {});
struct SpatialKineticOptions {
    EedfGrid grid;
    double dt_s{0};
    // a0*f - a1*f_old + a2*f_older = dt*C(f). BE=(1,1,0).
    double a0{1}, a1{1}, a2{0};
    double relative_tolerance{1e-8};
    int maximum_newton_iterations{15};
    int maximum_krylov_iterations{200};
    int krylov_restart{30};
    // Independent CPU Newton/JVP oracle. Production uses assembled row-owned
    // frozen-coefficient CSR with resident MPI/CPU/CUDA linear algebra.
    bool matrix_free_reference{false};
    backend::ExecutionOptions execution;
    // Optional cylindrical geometry coefficient 2*pi*A_meridional/V. Empty means
    // Cartesian. RZ basis is (radial, azimuthal, axial).
    std::vector<double> inverse_radius_m1;
};
struct SpatialKineticResult {
    // Rank-local cells in the distributed overload; complete cells in serial.
    std::vector<SphericalEedfState> states;
    std::vector<ElectronFeedback> feedback;
    std::vector<std::size_t> face_global_ids;
    std::vector<double> face_particle_flux_m2_s;
    std::vector<double> face_energy_flux_ev_m2_s;
    double outward_particle_rate_s{0};
    double outward_energy_rate_ev_s{0};
    double nonlinear_residual{0};
    int nonlinear_iterations{0};
    int krylov_iterations{0};
    double linear_residual{0};
    std::string linear_device;
    std::size_t local_unknowns{0};
    std::size_t halo_unknowns{0};
};
// Collision, vector-field acceleration, and conservative spatial streaming enter
// one implicit block. Heavy inventories/temperatures are the outer Newton trial.
// Electron number/energy equations are the F00 moments of this block, never a
// second fluid equation. Inputs are immutable on nonlinear/linear failure.
[[nodiscard]] SpatialKineticResult advance_spatial_kinetics(const KineticMesh& mesh,
    std::span<const SpatialKineticCell> cells,
    std::span<const SphericalEedfState> previous,
    const SpatialKineticOptions& options,
    std::span<const SphericalEedfState> older = {});
// Cell-contiguous ownership. No prognostic vector or Krylov basis is gathered.
// Only neighboring cell distributions cross rank boundaries. Empty ranks work.
[[nodiscard]] SpatialKineticResult advance_distributed_spatial_kinetics(const KineticMesh& mesh,
    backend::RowPartition cell_partition,
    std::span<const SpatialKineticCell> local_cells,
    std::span<const SphericalEedfState> local_previous,
    const SpatialKineticOptions& options,
    std::span<const SphericalEedfState> local_older = {});
} // namespace raikou::kinetics
