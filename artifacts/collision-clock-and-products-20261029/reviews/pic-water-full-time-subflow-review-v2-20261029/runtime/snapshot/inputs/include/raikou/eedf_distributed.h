#pragma once
#include "raikou/distributed.h"
#include "raikou/electron_boltzmann.h"

namespace raikou::kinetics {
// The partition is over packed mode rows and must start/end on an energy-cell
// boundary. Local storage uses the same even/odd mode ordering as EedfState.
struct DistributedEedfState {
    EedfGrid grid;
    backend::RowPartition partition;
    double time_s{};
    std::vector<double> density_modes;
};
struct FrozenEedfOperator {
    backend::DistributedCsrMatrix matrix;
    std::vector<double> mass;
    std::vector<unsigned char> padding;
    std::size_t ghost_modes{};
};
struct DistributedEedfResult {
    DistributedEedfState state;
    // Scalar swarm properties and process rates are global; energy_ev/eepf
    // contain only owned energy cells, never a replicated global distribution.
    SwarmResult swarm;
    // Hydrodynamic gradient diffusion uses the owned CPU matrix-free solver,
    // including when the main collision/field linear executor is CUDA.
    int nonlinear_iterations{};
    double linear_residual{}, energy_balance_relative{};
    std::size_t ghost_modes{}, local_matrix_entries{};
};
[[nodiscard]] backend::RowPartition distributed_eedf_partition(const EedfGrid&, const backend::ExecutionOptions& = {});
[[nodiscard]] DistributedEedfState maxwellian_eedf_distributed(double density_m3, double mean_energy_ev,
    const EedfGrid&, const backend::ExecutionOptions& = {});
// Linear operator with one physical frozen Landau flux (not its doubled Newton
// derivative). Frozen coefficients are computed from isotropic_state using
// local faces and scalar prefix/suffix moments. mass_normalized=true returns
// continuous physical rates in 1/s; false returns the legacy per-N mass form.
[[nodiscard]] FrozenEedfOperator freeze_eedf_operator(const std::vector<CollisionProcess>&,
    const EedfConditions&, const DistributedEedfState& isotropic_state,
    int azimuthal_order = 0, bool mass_normalized = true, bool normalized_real_harmonics = true);
[[nodiscard]] DistributedEedfResult advance_eedf_distributed(const std::vector<CollisionProcess>&,
    const EedfConditions&, const DistributedEedfState&, double dt_s);
[[nodiscard]] DistributedEedfResult solve_eedf_distributed(const std::vector<CollisionProcess>&,
    const EedfConditions&, const EedfGrid& = {});
[[nodiscard]] raikou::Json advance_distributed_eedf_json(const raikou::Json&);
[[nodiscard]] raikou::Json solve_distributed_eedf_json(const raikou::Json&);
}
