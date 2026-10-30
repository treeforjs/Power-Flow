#pragma once

#include "raikou/json.h"
#include "raikou/backend.h"

#include <cstddef>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace raikou::chemistry {

// All concentrations are particles/m^3 of the species' own phase. Rate
// coefficients use the corresponding SI particle-number convention, not mol/L.
struct Species {
    std::string id;
    std::string phase = "gas";
    std::map<std::string, int> elements;
    int charge = 0;
    double mass_kg = 0.0;
    bool nonnegative = true;
};

struct Term {
    std::size_t species = 0;
    double coefficient = 1.0;
};

enum class RateKind {
    Constant,
    Arrhenius,
    ElectronTemperatureTable,
    ElectronEedf,
    Expression,
    Parameter
};

struct RateLaw {
    RateKind kind = RateKind::Constant;
    double A = 0.0;
    double exponent = 0.0;
    double activation_temperature_K = 0.0;
    double reference_temperature_K = 1.0;
    // "gas", "liquid", or "electron" (converted from eV to K).
    std::string temperature = "gas";
    std::vector<double> temperatures_eV;
    std::vector<double> coefficients;
    std::string key;
    std::string expression;
    std::map<std::string, double> constants;
    double minimum_temperature_K = 0.0;
    double maximum_temperature_K = 0.0; // zero means unspecified
};

struct Reaction {
    std::string id;
    std::vector<Term> reactants;
    std::vector<Term> products;
    // Empty means reactant stoichiometry. Explicit orders may include catalysts.
    std::vector<Term> orders;
    bool explicit_orders = false;
    RateLaw rate;
    // Binds an electron rate to exactly the collision used by the EEDF solver.
    std::string collision_label;
    // Volumetric rates are per m^3 of basis_phase; area rates per m^2 of
    // basis_interface. Area/volume factors are caller-supplied, never guessed.
    std::string basis_phase = "gas";
    std::string basis_interface;
    // Signed energy change per event, keyed "phase/reservoir". A positive value
    // adds energy. If checked, the sum must equal external_energy_eV_per_event.
    std::map<std::string, double> energy_eV_per_event;
    bool check_energy_balance = false;
    double external_energy_eV_per_event = 0.0;
    std::string source;
    std::string status = "research";
    Json metadata = Json::object();
};

struct Mechanism {
    std::string id;
    std::string status = "research";
    std::vector<Species> species;
    std::vector<Reaction> reactions;
    Json provenance = Json::array();
    Json coverage = Json::object();
    static Mechanism from_json(const Json& value);
    [[nodiscard]] Json to_json() const;
    [[nodiscard]] std::vector<std::string> validate() const;
    [[nodiscard]] std::size_t species_index(std::string_view id_value) const;
};

struct State {
    std::vector<double> densities_m3;
};

struct Environment {
    double gas_temperature_K = 300.0;
    double electron_temperature_eV = 1.0;
    double liquid_temperature_K = 300.0;
    std::map<std::string, double> electron_rates;
    std::map<std::string, double> parameters;
    std::map<std::string, double> phase_volumes_m3;
    std::map<std::string, double> interface_areas_m2;
};

struct Evaluation {
    std::vector<double> production_rates_m3_s;
    std::vector<double> jacobian_row_major;
    // Reaction progress per m^3/s or m^2/s according to its declared basis.
    std::vector<double> reaction_rates_m3_s;
    std::vector<double> rate_coefficients;
    std::map<std::string, double> energy_sources_W_m3;
    std::vector<std::string> diagnostics;
};

struct AdvanceOptions {
    double relative_tolerance = 1.0e-5;
    double absolute_tolerance_m3 = 1.0e-6;
    double minimum_step_s = 1.0e-20;
    double maximum_step_s = 0.0; // zero means the requested interval
    std::size_t maximum_steps = 10000;
    std::size_t maximum_newton_iterations = 30;
    // Empty uses Species::nonnegative. A supplied mask must match the species.
    std::vector<bool> nonnegative;
    // One local CR system can distribute its Newton linear operator. A host
    // using cell worker threads must leave mpi=false for those local calls.
    backend::ExecutionOptions execution;
};

struct AdvanceResult {
    State state;
    bool success = false;
    double advanced_time_s = 0.0;
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    std::string message;
    std::map<std::string, double> integrated_energy_J_m3;
};

[[nodiscard]] Evaluation evaluate(const Mechanism& mechanism, const State& state,
                                  const Environment& environment = {});
// Compiled/distributed assemblers validate the mechanism and environment once,
// then request coefficients without evaluating a global state or dense Jacobian.
void validate_environment(const Environment& environment);
[[nodiscard]] double rate_coefficient(const Reaction& reaction,const Environment& environment);
// Transactional: failure returns the original state and zero advanced_time_s.
[[nodiscard]] AdvanceResult advance(const Mechanism& mechanism, const State& state,
                                    const Environment& environment, double interval_s,
                                    const AdvanceOptions& options = {});
[[nodiscard]] Json coverage_report(const Mechanism& mechanism,
                                   const Environment& environment = {});
[[nodiscard]] Json validate_mechanism(const Json& value);
// Request shape: {mechanism, state:{densities_m3:{id:number}}, environment:{...}}.
[[nodiscard]] Json evaluate_json(const Json& request);
[[nodiscard]] Json advance_json(const Json& request);
// Shared mechanism with cells:[{state, environment?}], common environment/options.
// Cells advance independently; per-cell numerical failures retain original state.
[[nodiscard]] Json evaluate_batch_json(const Json& request);
[[nodiscard]] Json advance_batch_json(const Json& request);
[[nodiscard]] Json evaluate_distributed_json(const Json& request);
[[nodiscard]] Json advance_distributed_json(const Json& request);

} // namespace raikou::chemistry
