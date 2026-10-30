#pragma once
#include "raikou/chemistry.h"
#include "raikou/electron_boltzmann.h"
#include "raikou/spherical_harmonics.h"

namespace raikou::kinetics {

struct ElectronCoupling {
    std::vector<CollisionProcess> processes;
    EedfGrid grid;
    EedfConditions conditions;
    std::string electron_species{"e"};
    // A collision target can name one chemistry species or an explicit sum of resolved
    // states. No element/formula heuristics are used to invent populations.
    std::map<std::string, std::vector<std::string>> targets;
    // Zero uses the current gas heavy-particle count. A positive value is a fixed reference
    // N, with every target fraction still recalculated from the evolving species inventory.
    double gas_density_m3{0};
    bool strict_event_mapping{false};
    static ElectronCoupling from_json(const raikou::Json& value);
};

struct ElectronFeedback {
    EedfConditions conditions;
    SwarmResult swarm;
    chemistry::Environment environment;
    double mobility_m2_V_s{0};
    double diffusion_m2_s{0};
    double kinetic_electron_source_m3_s{0};
    double mapped_electron_source_m3_s{0};
    bool collision_event_mapping_complete{false};
    // The reservoirs are the operator's energy moments, not additional chemistry sources.
    // Add each physical transfer to an energy equation exactly once.
    std::map<std::string, double> collision_energy_sources_W_m3;
    std::vector<std::string> diagnostics;
};

[[nodiscard]] EedfConditions coupled_eedf_conditions(const chemistry::Mechanism& mechanism,
    const chemistry::State& state, const chemistry::Environment& environment,
    double field_magnitude_V_m, const ElectronCoupling& coupling);

// Recompute all populations, EEDF rates and transport inside an outer nonlinear residual.
// Reactions bind to collision labels; label+".reverse" binds superelastic excitation.
[[nodiscard]] ElectronFeedback evaluate_electron_feedback(const chemistry::Mechanism& mechanism,
    const chemistry::State& state, const chemistry::Environment& environment,
    double field_magnitude_V_m, const ElectronCoupling& coupling);
// The same, continuing from the last steady solve of this cell (solve_eedf's SteadyEedfReuse).
[[nodiscard]] ElectronFeedback evaluate_electron_feedback(const chemistry::Mechanism& mechanism,
    const chemistry::State& state, const chemistry::Environment& environment,
    double field_magnitude_V_m, const ElectronCoupling& coupling, SteadyEedfReuse& reuse);

// The local field approximation: steady swarms solved once at declared reduced fields, target fractions
// and gas temperature, then read at each cell's own E/N. Rates and transport between two tabulated fields
// are interpolated geometrically in E/N (linearly where a value changes sign); below the first field the
// first swarm holds and above the last the last one, which `at` reports. Electron-electron collisions,
// the cell's composition and its gas temperature do not enter the tabulated distributions.
struct SteadyEedfTable {
    std::vector<double> reduced_field_td;  // strictly increasing, positive
    EedfConditions conditions;             // as tabulated; reduced_field_td ignored
    std::vector<SwarmResult> swarms;
    static SteadyEedfTable build(const std::vector<CollisionProcess>& processes, const EedfConditions& conditions,
                                 std::vector<double> reduced_fields_td, const EedfGrid& grid);
    [[nodiscard]] SwarmResult at(double reduced_field_td, bool* above = nullptr) const;
};
// The steady feedback read from the table at the cell's E/N; ionization and attachment per N are summed
// again over the cell's own target fractions so mapped and kinetic electron sources agree.
[[nodiscard]] ElectronFeedback evaluate_electron_feedback(const chemistry::Mechanism& mechanism,
    const chemistry::State& state, const chemistry::Environment& environment,
    double field_magnitude_V_m, const ElectronCoupling& coupling, const SteadyEedfTable& table);

// Bind rates/energy diagnostics of an already evaluated full angular distribution.
// This performs no additional kinetic time step or steady-state approximation.
[[nodiscard]] ElectronFeedback spherical_electron_feedback(const chemistry::Mechanism& mechanism,
    const chemistry::State& candidate_state, const chemistry::Environment& environment,
    const Vector3& electric_field_V_m, const ElectronCoupling& coupling,
    const SphericalEedfState& electrons, const SphericalEedfEvaluation& evaluation);

struct TransientElectronFeedback {
    EedfStepResult kinetic_step;
    ElectronFeedback feedback;
};
// Trial state for a fully coupled backward-Euler residual. The caller supplies candidate
// heavy populations; the resulting F0 moments replace the electron density/energy equations.
// The caller must not also integrate the mapped electron collision source in chemistry.
[[nodiscard]] TransientElectronFeedback advance_electron_feedback(const chemistry::Mechanism& mechanism,
    const chemistry::State& candidate_state, const chemistry::Environment& environment,
    double field_magnitude_V_m, const ElectronCoupling& coupling,
    const EedfState& previous_electrons, double dt_s);

struct CoupledCellOptions {
    double relative_tolerance{1e-7};
    double absolute_density_tolerance_m3{1};
    double temperature_tolerance_K{1e-6};
    std::size_t maximum_newton_iterations{25};
    // A phase omitted here retains its prescribed temperature. Capacities are volumetric
    // J/(m^3 K) in that phase. gas and either liquid or aqueous are supported; liquid
    // and aqueous share Environment::liquid_temperature_K and cannot both be unknowns.
    std::map<std::string, double> heat_capacity_J_m3_K;
    // Explicit chemistry energy reservoirs feeding each phase's thermal equation. The gas
    // elastic EEDF transfer is included automatically, exactly once. Electron reservoirs
    // cannot be selected: electron energy belongs exclusively to kinetic moments.
    std::map<std::string, std::vector<std::string>> thermal_reservoirs;
};

struct CoupledCellResult {
    chemistry::State state;
    chemistry::Environment environment;
    EedfState electrons;
    ElectronFeedback feedback;
    chemistry::Evaluation chemistry;
    std::map<std::string, double> thermal_energy_change_J_m3;
    std::size_t nonlinear_iterations{0};
    double nonlinear_residual{0};
};

// One monolithic backward-Euler step. Heavy species, kinetic electron number/energy,
// and requested phase temperatures share one convergence/acceptance decision. The input
// objects remain unchanged on every failure. Other electron sources (photoemission etc.)
// require explicit distribution source operators and are rejected here.
[[nodiscard]] CoupledCellResult advance_coupled_cell(const chemistry::Mechanism& mechanism,
    const chemistry::State& initial, const chemistry::Environment& environment,
    const ElectronCoupling& coupling, const EedfState& previous_electrons,
    double field_magnitude_V_m, double dt_s, const CoupledCellOptions& options = {});
[[nodiscard]] raikou::Json advance_coupled_cell_json(const raikou::Json& request);

} // namespace raikou::kinetics
