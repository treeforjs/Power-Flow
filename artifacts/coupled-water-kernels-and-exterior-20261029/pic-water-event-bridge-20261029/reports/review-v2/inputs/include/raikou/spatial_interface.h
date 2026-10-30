#pragma once
#include "raikou/spatial.h"

namespace raikou::spatial::phase_interface {
struct ReactionEvaluation {
    chemistry::Evaluation own_volume_sources;
    std::vector<double> particle_rates_s;
    std::vector<double> particle_jacobian_s_m3;
    std::map<std::string, double> energy_powers_W;
    double electric_power_W = 0;
    double phase_a_charge_rate_A = 0;
    double phase_b_charge_rate_A = 0;
};
// State entries are own-phase concentrations read from the cell owning each
// species phase. Environment supplies the two actual positive phase volumes and
// their respective temperatures. The reaction must be authored per interface
// area. No transfer coefficient, equilibrium law, or energy partition is guessed.
[[nodiscard]] ReactionEvaluation evaluate_reaction(const chemistry::Mechanism& mechanism,
    const FaceReaction& binding, const chemistry::State& phase_owned_concentrations,
    const chemistry::Environment& environment, double area_m2,
    double phase_a_potential_V, double phase_b_potential_V);
} // namespace raikou::spatial::phase_interface
