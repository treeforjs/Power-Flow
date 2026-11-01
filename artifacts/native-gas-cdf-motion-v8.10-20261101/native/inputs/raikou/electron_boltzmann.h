#pragma once
// Electron energy distribution (EEDF) from the Boltzmann equation, expanded in Legendre
// polynomials to any order (multi-term; two terms is the classical approximation), state to state, with the processes that steady swarm codes leave out: superelastic collisions
// from excited states (each rotational and vibrational level is its own target with its own
// population), electron-electron and electron-ion Coulomb collisions, recombination with
// positive ions as an energy-dependent electron sink, attachment, and ionization with either
// equal energy sharing or the secondary born at rest.
//
// With f(v, mu) = sum_l F_l(e) P_l(mu) (mu the cosine to the acceleration), the spatially
// uniform equation divided by N v is, for each l (gamma = sqrt(2q/m), E/N in V m^2),
//   (1/(N v)) dF_l/dt + (E/N) [ l/(2l-1) (F_{l-1}' - (l-1)/(2e) F_{l-1})
//                              + (l+1)/(2l+3) (F_{l+1}' + (l+2)/(2e) F_{l+1}) ] = C_l / (N v)
// with C_0 the elastic recoil (Davydov) operator, the inelastic, superelastic, ionization,
// attachment and recombination gains and losses, and electron-electron collisions; and C_l
// (l >= 1) the total loss -sigma_0 F_l plus in-scattering weighted by the l-th Legendre moment
// sigma_l of each differential cross section. Dipole (Born) rotational collisions and screened
// elastic scattering have closed-form moments: sigma_l / sigma_0 = Q_l(a) / Q_0(a), Q_l the
// Legendre functions of the second kind. Superelastic collisions follow from detailed balance.
// Recombination, attachment and ionization change the electron number, so the steady solution
// is a population growing (or decaying) at its net rate, a term in every F_l. Energies are in
// eV, cross sections in m^2, rate coefficients in m^3/s, E/N in Td.
//
// Numerics: even l on the cells and odd l on the faces of a grid that is fine at low energy
// (faces at max (k/n)^p); a Scharfetter-Gummel flux for recoil and electron-electron drift and
// diffusion in energy; exact, particle-conserving overlap integrals for the energy jumps of
// inelastic and superelastic collisions in every l; GMRES with a banded LU preconditioner; the
// steady swarm is the dominant eigenvector, found by implicit inverse iteration.
#include <map>
#include <string>
#include <vector>
#include "raikou/json.h"
#include "raikou/backend.h"

namespace raikou::kinetics {

enum class CollisionKind { Elastic, Effective, Excitation, Ionization, Attachment, Recombination };

struct CollisionProcess {
    CollisionKind kind{CollisionKind::Elastic};
    std::string target;              // the species or state struck, e.g. "H2O" or "H2O(000,1_01)"
    std::string product;             // excitation: the state produced (superelastic from it)
    std::string label;
    double mass_ratio{0};            // elastic and effective: m_e / M
    double threshold_ev{0};          // excitation and ionization energy loss
    double weight_ratio{0};          // excitation: g_target / g_product; > 0 enables superelastic
    enum class Group { Automatic, Rotational, Vibrational, Electronic } group{Group::Automatic};  // for the energy budget;
                                     // automatic: below 0.05 eV rotational, below 2 eV vibrational
    std::vector<double> energy_ev;   // increasing; the cross section is linear between points,
    std::vector<double> cross_section_m2;  // zero below the first point and constant past the last
                                     // (elastic and effective: the momentum transfer)
    // Angular distribution, for every term beyond the first: isotropic; Born dipole (the
    // point-dipole rotational form, from the threshold); or, for elastic scattering, a screened
    // forward-peaked form fitted at each energy to sigma_m / sigma_0 from integral_cross_section_m2.
    enum class Angular { Isotropic, BornDipole, Screened } angular{Angular::Isotropic};
    std::vector<double> integral_cross_section_m2;  // Screened elastic: sigma_0 on energy_ev
};
// The LXCat text format: ELASTIC, EFFECTIVE, EXCITATION, IONIZATION and ATTACHMENT blocks, plus
// RECOMBINATION (target is the ion; no parameter line value is needed). An EXCITATION parameter
// line "threshold weight_ratio" enables superelastic collisions; "A <-> B" names the product.
[[nodiscard]] std::vector<CollisionProcess> read_lxcat(const std::string& text);

enum class IonizationSharing { Equal, OneTakesAll };

struct EedfGrid {
    int cells{1000};
    double max_energy_ev{100};
    double power{2};                 // faces at max (k/cells)^power; 1 is uniform
    int legendre_terms{2};           // F_0 .. F_{terms-1}; 2 is the two-term approximation
};

struct EedfConditions {
    double reduced_field_td{100};
    double gas_temperature_k{300};
    std::map<std::string, double> fractions;  // target densities over N; missing targets are 0
    double gas_density_m3{2.414e25};          // N, for the Coulomb terms
    double electron_fraction{0};              // n_e / N; electron-electron collisions when > 0
    double ion_fraction{-1};                  // sum z_i^2 n_i / N for positive-ion Coulomb momentum transfer
                                              // (-1: equal to electron_fraction)
    IonizationSharing sharing{IonizationSharing::Equal};
    bool temporal_growth{true};               // the population grows (or decays) at its net rate
    // MPI steady/advance calls use owned energy/mode rows and scalar Landau
    // prefixes; their linear solves run on the selected CPU/CUDA executor.
    // Auxiliary gradient diffusion remains on the owned CPU matrix-free path.
    // Serial calls retain the banded kinetic solver. MPI lifecycle belongs to
    // the host/CLI. Default one thread preserves small-solve cost.
    backend::ExecutionOptions execution{.threads = 1};
};

struct ProcessRate {
    std::string label;
    CollisionKind kind{CollisionKind::Elastic};
    std::string target;
    double fraction{0};
    double rate_m3_s{0};              // forward rate coefficient per target density
    double superelastic_rate_m3_s{0}; // the reverse, per density of the product state
    double power_ev_m3_s{0};          // net energy lost per electron per N: x u k - x* u k_s
};

// Energy per electron per unit gas density, eV m^3/s. In a steady swarm the field input equals
// elastic + inelastic (net of superelastic) + the energy carried off by attached and recombined
// electrons + the energy given to the population's own growth (net growth rate x mean energy);
// `residual` is what is left, relative to the field input.
struct EnergyBudget {
    double field{0};
    double elastic{0};
    double rotational{0}, vibrational{0}, electronic{0}, ionization{0};  // by threshold class
    double superelastic_gain{0};      // included (negative) in the classes above
    double attachment{0}, recombination{0};
    double growth{0};
    double electron_electron{0};      // exchanged among electrons; zero up to discretisation
    double residual{0};
};

struct SwarmResult {
    bool converged{false};
    int iterations{0};
    double mean_energy_ev{0};
    double mobility_n{0};             // mu N, 1/(V m s)
    double diffusion_n{0};            // D N, 1/(m s): isotropic two-term estimate
    double drift_velocity_m_s{0};     // flux drift velocity
    // From the density-gradient terms (first order in the gradient expansion, to every Legendre
    // order): flux longitudinal and transverse diffusion and the bulk drift velocity, which differs
    // from the flux drift when ionization, attachment or recombination depend on energy.
    double longitudinal_diffusion_n{0}, transverse_diffusion_n{0};  // 1/(m s)
    double bulk_drift_velocity_m_s{0};
    double gradient_residual{0};      // relative residual of the gradient solves
    double ionization_rate_m3_s{0}, attachment_rate_m3_s{0}, recombination_rate_m3_s{0};  // per N: sum of fraction x rate
    double alpha_n_m2{0}, eta_n_m2{0};  // Townsend coefficients over N from the flux
    double net_rate_m3_s{0};          // (ionization - attachment - recombination) over N
    std::vector<ProcessRate> rates;
    EnergyBudget budget;
    std::vector<double> energy_ev;    // cell centres
    std::vector<double> eepf;         // F(e), normalised: integral sqrt(e) F de = 1
};

// Restartable axisymmetric velocity-space state. Every Legendre mode is density weighted;
// electron number and electron energy are moments, never independently evolved variables.
// Cell i contains the even modes followed by the odd modes at face i+1. The last odd
// entries are zero boundary padding. The field direction must remain fixed during a step;
// this representation does not retain non-axisymmetric spherical-harmonic modes.
struct EedfState {
    EedfGrid grid;
    double time_s{0};
    std::vector<double> density_modes;
    [[nodiscard]] double electron_density_m3() const;
    [[nodiscard]] double electron_energy_ev_m3() const;
    [[nodiscard]] double mean_energy_ev() const;
};

struct EedfStepResult {
    EedfState state;
    SwarmResult swarm;
    int nonlinear_iterations{0};
    double linear_residual{0};
    double energy_balance_relative{0};
};

[[nodiscard]] EedfState maxwellian_eedf(double density_m3, double mean_energy_ev,
                                      const EedfGrid& grid = {});
// One backward-Euler step, implicit in every mode and in electron-electron collisions.
// Fractions are the current actual target/ion populations divided by gas_density_m3.
// All channels are retained even when their current population is zero. The caller's state
// is unchanged on failure; reduce dt or update the outer coupled nonlinear iterate and retry.
[[nodiscard]] EedfStepResult advance_eedf(const std::vector<CollisionProcess>& processes,
                                        const EedfConditions& conditions, const EedfState& initial,
                                        double dt_s);

[[nodiscard]] SwarmResult solve_eedf(const std::vector<CollisionProcess>& processes, const EedfConditions& conditions,
                                     const EedfGrid& grid = {});

// For a caller that solves nearby states again and again, such as a fluid closure solving every
// cell at every Newton assembly. `modes` is the last converged distribution (every Legendre mode)
// on the same grid: the solve continues from it, starting at the pseudo time step a cold solve
// converges at, so a nearby state takes a few iterations instead of the whole ramp, and a
// continuation that fails falls back to the cold solve. Empty starts cold; either way it returns
// holding this solve's distribution.
// gradient_terms false skips the density-gradient (bulk) transport, which leaves
// longitudinal_diffusion_n, transverse_diffusion_n, bulk_drift_velocity_m_s and gradient_residual 0.
struct SteadyEedfReuse {
    std::vector<double> modes;
    bool gradient_terms{true};
};
[[nodiscard]] SwarmResult solve_eedf(const std::vector<CollisionProcess>& processes, const EedfConditions& conditions,
                                     const EedfGrid& grid, SteadyEedfReuse& reuse);

// Several reduced fields on one grid: the collision operator is assembled once and each field
// starts from the previous one's distribution, so a sweep in increasing E/N costs little more
// than its solves. conditions.reduced_field_td is ignored.
[[nodiscard]] std::vector<SwarmResult> solve_eedf_sweep(const std::vector<CollisionProcess>& processes, const EedfConditions& conditions,
                                                        const std::vector<double>& reduced_fields_td, const EedfGrid& grid = {});

// Steady-state Townsend (spatial growth) coefficients, what a Townsend gap and a Paschen curve
// see: a steady stream n ~ exp(alpha_eff z) obeys D_L alpha^2 - W alpha + nu = 0 with the bulk
// drift W, longitudinal diffusion D_L and net rate nu, so
//   alpha_eff / N = (W - sqrt(W^2 - 4 (N D_L) (nu / N))) / (2 N D_L),
// split into ionization and loss in proportion to their rates (nu / W when D_L is unknown).
struct SteadyTownsend { double alpha_n_m2{0}, eta_n_m2{0}; };
[[nodiscard]] SteadyTownsend steady_townsend(const SwarmResult& swarm);

// Time-dependent breakdown of a uniform gap from a few seed electrons, in one volume (0D): the
// distribution evolves in time from the seed shape (multi-term), the electron number grows by ionization and falls by attachment,
// recombination and drift to the anode, positive and negative ions drift to their electrodes,
// ions reaching the cathode release secondary electrons, and the gap is a capacitor fed from
// the source through a series resistance, so the voltage collapses as the current grows.
//   dn_e/dt = (nu_i - nu_a - nu_r) n_e - (W / d) n_e + gamma (W_i / d) n_+
//   dn_+/dt = nu_i n_e - (W_i / d) n_+ - k_r n_e n_+,   dn_-/dt = nu_a n_e - (W_i / d) n_-
//   C dV/dt = (V0 - V) / R - I,   I = e A (n_e W + (n_+ + n_-) W_i),   C = eps0 A / d
// Every Legendre term carries its time derivative (at low pressure momentum relaxes on the time
// scale of the breakdown). Time stepping is second order (variable-step BDF2), implicit in the
// distribution and the circuit. Every energy channel is integrated with the same scheme as the distribution, so the
// electrons' books close to the solver tolerance: field work = elastic + rotational +
// vibrational + electronic + ionization + attachment + recombination + electron-electron +
// energy carried to the anode - energy brought by secondaries + the change in electron energy.
// The gas composition is held fixed (excited and ionized fractions small).
struct BreakdownSetup {
    double gap_m{0.01};
    double area_m2{1e-4};
    double voltage_v{0};                  // source voltage; the gap starts charged to it
    double series_resistance_ohm{50};     // 0: an ideal source holds the gap voltage
    double seed_density_m3{1e9};
    double seed_mean_energy_ev{1};        // seed and secondary electrons: Maxwellian of this mean
    double ion_mobility_n{0};             // N mu_i, 1/(V m s), both ion signs; 0: ions stay put
    double secondary_yield{0};            // gamma, electrons per ion reaching the cathode
    bool drift_loss{true};                // electrons leave through the anode at W / d
    double end_time_s{1e-6};
    int outputs{200};                     // samples evenly spaced in time, plus one whenever the electron
                                          // density grows by a quarter decade or the voltage moves 1%
    int max_steps{50000};
    double collapse_fraction{0.5};        // stop once the gap voltage falls below this share of the source
    double step_tolerance{1};             // scales the per-step change limits (smaller: more accurate)
};
struct BreakdownSample {
    double time_s{0}, gap_voltage_v{0}, reduced_field_td{0};
    double electron_density_m3{0}, positive_ion_density_m3{0}, negative_ion_density_m3{0};
    double current_a{0}, mean_energy_ev{0}, drift_velocity_m_s{0};
    double ionization_rate_m3_s{0}, attachment_rate_m3_s{0};  // per N
};
// Cumulative energies, J. The electron channels are those of EnergyBudget (net of superelastic
// collisions); residuals are relative to the largest flow in each set of books.
struct BreakdownEnergy {
    double source{0}, resistor{0}, capacitor_change{0}, gap{0};      // circuit: source = resistor + capacitor change + gap
    double electron_field_work{0}, ion_joule{0};                    // gap = electron field work + ion Joule heating
    double elastic{0}, rotational{0}, vibrational{0}, electronic{0}, ionization{0};
    double attachment{0}, recombination{0}, electron_electron{0}, superelastic_gain{0};
    double anode_loss{0}, secondary_input{0}, electron_energy_change{0};
    double electron_residual{0}, gap_residual{0}, circuit_residual{0};
};
struct BreakdownResult {
    std::vector<BreakdownSample> samples;
    BreakdownEnergy energy;
    double breakdown_time_s{-1};          // the gap voltage first 10% below the source (-1: not reached)
    std::string stop;                     // "end time", "collapse" or "step limit"
    int steps{0}, rejected_steps{0};
    std::vector<double> energy_ev, eepf;  // the final distribution, normalised
};
[[nodiscard]] BreakdownResult solve_breakdown(const std::vector<CollisionProcess>& processes, const EedfConditions& conditions,
                                              const BreakdownSetup& setup, const EedfGrid& grid = {});

// JSON entrypoints used by the common CLI, C ABI, and Python runtime.
[[nodiscard]] std::vector<CollisionProcess> collision_processes_from_json(const raikou::Json& value);
[[nodiscard]] EedfConditions eedf_conditions_from_json(const raikou::Json& value);
[[nodiscard]] EedfGrid eedf_grid_from_json(const raikou::Json& value);
[[nodiscard]] raikou::Json solve_eedf_json(const raikou::Json& request);
[[nodiscard]] raikou::Json advance_eedf_json(const raikou::Json& request);
[[nodiscard]] raikou::Json water_levels_json(const raikou::Json& request);
[[nodiscard]] raikou::Json solve_breakdown_json(const raikou::Json& request);

} // namespace raikou::kinetics
