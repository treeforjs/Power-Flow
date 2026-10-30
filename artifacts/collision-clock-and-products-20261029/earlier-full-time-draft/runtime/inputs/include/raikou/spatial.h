#pragma once

#include "raikou/chemistry.h"
#include "raikou/electron_coupling.h"
#include "raikou/spatial_kinetics.h"
#include "raikou/json.h"
#include <array>
#include <cstddef>
#include <string>
#include <optional>
#include <vector>

namespace raikou::spatial {

// Cell-centred finite volumes. Planar 2D uses extent_m[2] as physical depth;
// axisymmetry uses x=r and y=z and integrates the full revolution.
struct Mesh {
    struct Cell { std::string type; std::vector<std::size_t> nodes; };
    struct BoundaryFace { std::vector<std::size_t> nodes; int side = 0; };
    int dimension = 2;
    bool axisymmetric = false;
    std::array<std::size_t, 3> shape{1, 1, 1};
    std::array<double, 3> extent_m{1, 1, 1};
    // Supplied convex SI tri3/quad4/tet4/hex8 cells with planar faces.
    // Conservative least-squares multipoint flux correction accepts skew cells.
    std::vector<std::array<double, 3>> nodes_m;
    std::vector<Cell> cells;
    std::vector<BoundaryFace> boundary_faces;
};
struct Transport {
    double mobility_m2_V_s = 0;
    double diffusion_m2_s = 0;
};
enum class BoundaryKind { Insulating, Electrode, Dielectric };
struct Boundary {
    struct Sheath {
        std::string species;
        std::string model = "maxwellian"; // maxwellian or bohm (positive ions)
        double temperature_eV = 0; // prescribed wall-closure temperature
        double electron_temperature_eV = 0; // required by Bohm closure
        double transmission = 1;
    };
    struct SecondaryEmission {
        std::string incident_species;
        std::string emitted_species;
        double yield = 0; // authored emitted particles / incident ion
    };
    BoundaryKind kind = BoundaryKind::Insulating;
    double potential_V = 0;
    bool absorbing = false;
    double capacitance_F_m2 = 0; // dielectric backing capacitance per area
    double initial_surface_charge_C_m2 = 0;
    std::string kinetic_electrons = "auto"; // auto: absorbing flag, otherwise specular
    std::optional<kinetics::SphericalEedfState> incoming_electrons;
    std::vector<Sheath> sheaths;
    std::vector<SecondaryEmission> secondary_emission;
    // Fluid electrons use this explicitly supplied wall energy per emission.
    // Full kinetics uses the exact energy moment of incoming_electrons instead.
    double secondary_energy_eV = 0;
};
struct Circuit {
    bool enabled = false;
    int driven_side = 1; // xmin,xmax,ymin,ymax,zmin,zmax
    double source_voltage_V = 0;
    double initial_voltage_V = 0;
    double series_resistance_ohm = 0; // zero: ideal voltage source
    double external_capacitance_F = 0;
};
struct Options {
    double time_step_s = 1e-9;
    double minimum_time_step_s = 1e-18;
    double end_time_s = 1e-8;
    double relative_tolerance = 1e-8;
    double absolute_density_tolerance_m3 = 1;
    double potential_tolerance_V = 1e-8;
    double surface_charge_tolerance_C_m2 = 1e-12;
    double circuit_charge_tolerance_C = 1e-20;
    std::size_t maximum_newton_iterations = 30;
    std::size_t maximum_steps = 100000;
    std::size_t threads = 1;
    bool mpi = false;
    int communicator_f = -1;
    bool bdf2 = true;
    bool adaptive_time = true;
    double temporal_relative_tolerance = 1e-3;
    double temporal_density_tolerance_m3 = 1;
    double temporal_temperature_tolerance_K = 1e-3;
    double temporal_potential_tolerance_V = 1e-5;
    double temporal_surface_charge_tolerance_C_m2 = 1e-12;
    std::string device = "cpu";
    int cuda_device = -1;
    std::string preconditioner = "auto"; // block ILU(0) on CPU and CUDA
    int linear_restart = 40;
    int linear_max_iterations = 4000;
    // "krylov": restarted GMRES with the preconditioner above; "direct": one sparse LU per Newton
    // iteration (serial; needs a build with Eigen), whose cost grows far more slowly with the mesh.
    std::string linear_solver = "krylov";
    double nonorthogonal_correction_limiter = 2; // zero explicitly selects unlimited correction
    // The steady Boltzmann closure starts each cell's solve from that cell's last one (the same
    // converged distribution to the solve tolerance, in a few iterations instead of the full ramp).
    bool eedf_warm_start = true;
};
struct Thermal {
    struct Phase {
        std::vector<double> heat_capacity_J_m3_K;
        std::vector<double> conductivity_W_m_K;
        std::vector<double> initial_temperature_K;
        std::vector<double> heating_W_m3;
        std::vector<std::string> heat_reservoirs;
    };
    bool enabled = false;
    std::vector<double> heat_capacity_J_m3_K;
    std::vector<double> conductivity_W_m_K;
    std::vector<double> initial_temperature_K;
    std::vector<double> heating_W_m3;
    std::vector<std::string> heat_reservoirs;
    std::vector<std::string> joule_species;
    bool use_eedf_elastic_heating = false;
    // Alternative to common temperature: per-own-phase energy equations.
    std::map<std::string, Phase> phases;
};
struct FaceReaction {
    std::string reaction_id;
    std::string phase_a = "gas";
    std::string phase_b = "aqueous";
    std::string electrical_work_to_phase;
};
struct InterfaceHeatTransfer {
    std::string phase_a = "gas";
    std::string phase_b = "aqueous";
    double conductance_W_m2_K = 0;
};
struct KineticPhaseContact {
    kinetics::SpatialKineticBoundary boundary;
    std::string absorbed_species;
    std::string energy_to_phase;
};
struct RadiationGroup {
    struct Absorber { std::string reaction_id; double cross_section_m2 = 0; };
    std::string id;
    double photon_energy_eV = 0;
    double speed_m_s = 299792458;
    double diffusion_m2_s = 0;
    double absorption_m1 = 0;
    std::map<std::string, double> absorption_m1_by_phase;
    std::vector<double> initial_density_m3;
    std::vector<double> emission_W_m3;
    std::string energy_reservoir;
    std::vector<Absorber> absorbers;
    std::array<bool, 6> escaping_sides{};
};
struct State {
    double time_s = 0;
    std::vector<double> densities_m3;
    std::vector<double> potential_V;
    // One entry per dielectric boundary face, in cell/side traversal order.
    std::vector<double> surface_charge_C_m2;
    std::vector<double> dielectric_potential_V; // one independent Robin face potential per surface
    double driven_voltage_V = 0;
    std::vector<double> temperature_K;
    std::map<std::string, std::vector<double>> phase_temperatures_K;
    std::vector<double> photon_densities_m3; // cell-major, group-minor
    // Full density-weighted Legendre modes, one restartable state per cell.
    std::vector<kinetics::EedfState> electrons;
    std::vector<kinetics::SphericalEedfState> spherical_electrons;
    std::vector<double> ohmic_charge_C_m3; // signed bulk charge, separate from transported species
};
struct Checkpoint {
    State state;
    State previous_state;
    double previous_step_s = 0;
    double next_step_s = 0;
    std::size_t accepted_steps = 0;
    std::vector<double> boundary_counts;
    std::vector<double> previous_boundary_counts;
    std::map<std::string, double> energy_integrals_J;
    std::map<std::string, double> previous_energy_integrals_J;
    std::vector<double> ohmic_heating_J, previous_ohmic_heating_J;
    double ohmic_boundary_charge_C = 0, previous_ohmic_boundary_charge_C = 0;
    // Per-cell physical J accumulated from the assembled thermal RHS, not a
    // difference of baseline temperatures. Empty phase key is common thermal.
    std::map<std::string, std::vector<double>> thermal_energy_integrals_J, previous_thermal_energy_integrals_J;
    std::map<std::string, std::vector<double>> thermal_energy_residual_bound_J, previous_thermal_energy_residual_bound_J;
    // False after an advance on a backend without per-cell thermal history.
    // Continuation retains this status rather than reconstructing lost heat.
    bool thermal_energy_ledger_available = true;
};
struct Electrical {
    std::string phase; // owner of background carriers, distinct from charged DD transport
    std::vector<double> ohmic_conductivity_S_m; // effective bulk coefficient; not phase-weighted again
    std::vector<double> initial_ohmic_charge_C_m3;
};
struct HeavyBarycentricPair {
    std::string neutral, ion;
};
struct Problem {
    Mesh mesh;
    chemistry::Mechanism mechanism;
    chemistry::Environment environment;
    std::optional<kinetics::ElectronCoupling> electron_closure;
    std::string electron_mode = "steady_boltzmann";
    // local_field_table: the reduced fields and the conditions its steady swarms are solved at.
    std::vector<double> electron_table_reduced_fields_td;
    kinetics::EedfConditions electron_table_conditions;
    std::string electron_table_key;  // the closure as requested: solves with the same key share one table in this process
    [[nodiscard]] bool steady_electrons() const { return electron_mode == "steady_boltzmann" || electron_mode == "local_field_table"; }
    std::array<double, 3> electron_field_axis{1, 0, 0};
    double initial_electron_mean_energy_eV = 0;
    std::vector<Transport> transport; // mechanism species order
    std::array<Boundary, 6> boundaries;
    std::vector<double> permittivity_F_m; // one value or one per cell
    std::vector<double> initial_densities_m3; // cell-major, species-minor
    std::map<std::string, std::vector<double>> phase_fractions;
    std::vector<FaceReaction> face_reactions;
    std::vector<InterfaceHeatTransfer> interface_heat_transfer;
    std::map<std::size_t, KineticPhaseContact> electron_phase_contacts;
    Thermal thermal;
    std::vector<RadiationGroup> radiation;
    Circuit circuit;
    Options options;
    std::optional<Checkpoint> initial_checkpoint;
    std::optional<Electrical> electrical;
    std::optional<HeavyBarycentricPair> heavy_barycentric_pair;
};
struct Sample {
    double time_s = 0;
    double driven_voltage_V = 0;
    double source_current_A = 0;
    double gas_charge_C = 0;
    double surface_charge_C = 0;
    double field_energy_J = 0;
    double electric_work_W = 0; // signed J.E from the conservative face flux
    double thermal_energy_J = 0;
    double radiation_energy_J = 0;
    double electron_energy_J = 0;
    double nonlinear_residual = 0; // <= 1 meets all declared tolerances
    double temporal_error = 0; // <=1 for an accepted adaptive timestep
    int order = 1;
    std::size_t newton_iterations = 0;
    double ohmic_charge_C = 0, ohmic_work_W = 0, ohmic_terminal_current_A = 0;
    double ohmic_reconstructed_work_W = 0;
    double barycentric_frame_work_W = 0, neutral_frame_particle_work_W = 0;
    double collisional_heavy_heating_W = 0;
};
struct Result {
    State state;
    bool success = false;
    std::string message;
    std::size_t accepted_steps = 0;
    std::size_t rejected_steps = 0;
    std::vector<Sample> samples;
    std::vector<double> particle_counts;
    std::vector<double> particles_collected_at_boundaries;
    std::map<std::string, double> chemical_energy_integrals_J;
    Checkpoint checkpoint;
    std::vector<double> ohmic_heating_J, ohmic_heating_W_m3;
    double ohmic_charge_collected_at_boundaries_C = 0;
    // Keys are thermal model phase names. Check checkpoint availability: a
    // serial continuation cannot recover an unavailable distributed history.
    std::map<std::string, std::vector<double>> thermal_energy_integrals_J, thermal_energy_residual_bound_J;
};

[[nodiscard]] double bernoulli(double x);
[[nodiscard]] double scharfetter_gummel_flux(double left_density_m3,
    double right_density_m3, double left_potential_V, double right_potential_V,
    double distance_m, int charge_number, const Transport& transport);
[[nodiscard]] std::vector<double> cell_volumes(const Mesh& mesh);
// RHS is charge density C/m3. Dirichlet electrode and dielectric Robin data
// are used; a boundary anchor is required. Throws on a failed linear solve.
[[nodiscard]] std::vector<double> solve_poisson(const Problem& problem,
    const std::vector<double>& charge_density_C_m3,
    const std::vector<double>& surface_charge_C_m2 = {},
    double driven_voltage_V = 0);
[[nodiscard]] Result solve(const Problem& problem);
// Metadata-only parsing keeps scalar species initial values unexpanded and leaves
// supplied arrays for the owned caller to distribute; no checkpoint is decoded.
[[nodiscard]] Problem problem_from_json(const Json& request, bool expand_initial_state = true);
[[nodiscard]] Json solve_json(const Json& request);
// Zero-time host material/charge update; returns a new local-time checkpoint and
// an electrical work audit, never a heat source. Serial steady-electron path.
[[nodiscard]] Json reconfigure_json(const Json& request);

} // namespace raikou::spatial
