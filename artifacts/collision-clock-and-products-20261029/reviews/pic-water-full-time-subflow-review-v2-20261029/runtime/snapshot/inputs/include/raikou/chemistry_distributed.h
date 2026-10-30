#pragma once

#include "raikou/chemistry.h"
#include "raikou/distributed.h"
#include <map>
#include <string>
#include <vector>

namespace raikou::chemistry {
// Species IDs are mechanism indices. Each rank stores its contiguous species
// interval only; static mechanism/rate metadata may be shared on every rank.
struct DistributedState {
    backend::RowPartition partition;
    std::vector<double> densities_m3;
};
struct DistributedEvaluation {
    std::vector<double> production_rates_m3_s;
    backend::DistributedCsrMatrix jacobian;
    // Global reservoir totals, summed once from uniquely owned reactions.
    std::map<std::string, double> energy_sources_W_m3;
    std::size_t ghost_species{};
};
struct DistributedAdvanceResult {
    DistributedState state;
    bool success{false};
    double advanced_time_s{};
    std::size_t accepted_steps{}, rejected_steps{};
    std::string message;
    std::map<std::string, double> integrated_energy_J_m3;
    std::size_t ghost_species{};
};
// Collective on execution.communicator_f. Host initializes MPI; rank ownership
// must cover the mechanism in communicator order. Numerical vectors/Jacobians
// remain owned+reactant halo throughout evaluation, Newton, and adaptation.
[[nodiscard]] DistributedEvaluation evaluate_distributed(
    const Mechanism&, const DistributedState&, const Environment& = {},
    const backend::ExecutionOptions& = {});
// Transactional across all ranks: failure returns each original local state,
// zero advanced_time_s and empty integrated energy. GPU execution uses resident
// distributed linear algebra; reaction assembly currently executes on the CPU.
[[nodiscard]] DistributedAdvanceResult advance_distributed(
    const Mechanism&, const DistributedState&, const Environment&, double interval_s,
    const AdvanceOptions& = {});
}
