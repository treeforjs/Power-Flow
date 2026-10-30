#pragma once
#include "raikou/spatial.h"
#include "raikou/spatial_partition.h"
#include <functional>
#include <memory>

namespace raikou::spatial {
struct DistributedInitialValues {
    // Called only for owned cells. Empty callbacks use unexpanded constants in
    // the physical model (one species vector / scalar thermal/photon fields).
    std::function<double(std::size_t cell, std::size_t species)> density;
    std::function<double(std::size_t cell, std::size_t field)> temperature;
    std::function<double(std::size_t cell, std::size_t group)> photon;
    std::function<double(std::size_t cell)> potential;
};
struct DistributedState {
    double time_s = 0;
    backend::RowPartition rows;
    partition::CellInterval cells;
    std::vector<double> owned_values;
    std::vector<kinetics::SphericalEedfState> spherical_electrons{};
    std::vector<kinetics::EedfState> electrons{}; // exact owned legacy local/fixed-axis history
};
struct DistributedCheckpoint {
    DistributedState state, previous_state;
    double previous_step_s = 0, next_step_s = 0;
    std::size_t accepted_steps = 0;
    std::vector<double> local_boundary_counts;
    std::vector<double> previous_local_boundary_counts;
    std::map<std::string, double> local_energy_integrals_J;
    std::map<std::string, double> previous_local_energy_integrals_J;
};
struct DistributedProblem {
    std::shared_ptr<const Problem> model{};
    DistributedInitialValues initial{};
    std::optional<partition::CellInterval> owned_cells{};
    std::optional<DistributedCheckpoint> checkpoint{};
    std::optional<bool> adaptive_time{}; // default: model.options.adaptive_time
};
struct DistributedResult {
    DistributedState state;
    bool success = false;
    std::string message;
    std::size_t accepted_steps = 0, rejected_steps = 0;
    std::vector<Sample> samples;
    std::vector<double> particle_counts; // reduced species inventories
    std::vector<double> particles_collected_at_boundaries;
    std::map<std::string, double> chemical_energy_integrals_J;
    DistributedCheckpoint checkpoint;
    std::size_t local_unknowns = 0, halo_unknowns = 0;
    std::size_t local_kinetic_unknowns = 0, kinetic_halo_unknowns = 0;
    std::size_t peak_local_matrix_entries = 0, peak_local_state_entries = 0;
};
// Cell blocks and full angular states remain owned. Geometry, mechanism and
// static material coefficient metadata may be shared; this API never gathers
// a prognostic vector, global residual, matrix or Krylov basis.
[[nodiscard]] DistributedResult solve_distributed(const DistributedProblem& problem);
// Shared scalar physical metadata plus optional rank-local initial fields or
// a rank-local checkpoint. No global prognostic vector is constructed/gathered.
[[nodiscard]] Json solve_distributed_json(const Json& request);
// Existing replicated API boundary: partition before advancing, gather completed
// current/previous states only when constructing the historical Result object.
[[nodiscard]] Result solve_mpi_compat(const Problem& problem);
} // namespace raikou::spatial
