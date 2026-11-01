#pragma once
#include <array>
#include <cstdint>
#include <vector>

namespace surface_emission {
inline constexpr double qe=1.602176634e-19,me=9.1093837139e-31,kb=1.380649e-23;
inline constexpr double planck=6.62607015e-34,eps0=8.8541878188e-12,c_light=299792458.;
using Vector=std::array<double,3>;
struct Config {
    double work_function_eV; // mandatory explicit input;4.3 is only a clean-surface benchmark
    int energy_panels=128,action_nodes=64;
    double tail_widths=48;
};
struct Node {double normal_energy_eV,probability,cdf,mean_transverse_eV;};
struct Distribution {
    Config config;
    double T_K=0,extraction_field_V_m=0,current_A_m2=0,log_current_A_m2=0;
    bool has_log_current=false,current_underflow=false;
    double mean_metal_energy_eV=0,wall_heat_W_m2=0,barrier_lowering_eV=0;
    double fermi_decay_width_eV=0;
    std::vector<Node> nodes; // quadrature-discrete normal spectrum; caller treats trial as immutable
};
double ideal_richardson_A_m2_K2();
double extraction_field(const Vector& solved_one_sided_total_E,const Vector& material_outward_unit_normal);
double schottky_lowering_eV(double positive_field_V_m);
double action(double normal_energy_eV,double positive_field_V_m,const Config&);
double transmission(double normal_energy_eV,double positive_field_V_m,const Config&);
Distribution evaluate(double temperature_K,double signed_extraction_field_V_m,const Config&);
double transverse_quantile_eV(double normal_energy_eV,double T_K,double open_uniform);

struct Stage {
    std::uint64_t accepted_state=0,mesh_revision=0,material_revision=0,patch_id=0;
    std::uint64_t field_accepted_state=0,thermal_accepted_state=0;
    double t=0,dt=0,stage_time=0,field_time=0,thermal_time=0;
    Vector total_E{},outward_normal{};
    double T_K=0,area_m2=0;
    bool material_temperature_phase_admitted=false;
    bool locally_planar_barrier_admitted=false;
    bool local_equilibrium_electron_supply_admitted=false;
    Vector total_B_T{}; // same field_time/accepted revision as solved total E
};
struct PreparedFlux {
    Stage stage;
    Distribution spectrum;
    double expected_physical_electrons=0;
    double mean_wall_heat_diagnostic_J=0; // never debit this AND sampled packet heat
};
PreparedFlux prepare(const Stage&,const Config&);

enum class PlaneKind {ThermalVacuumSurfaceStep,ExplicitSchottkyNordheimHandoff};
struct MatchingPlane {
    PlaneKind kind=PlaneKind::ThermalVacuumSurfaceStep;
    double outward_distance_m=0;
    bool microscopic_matching_admitted=false;
    bool field_path_work_owned_once=false;
    bool image_reference_port_admitted=false;
};
double required_SN_matching_distance_m(const Distribution&);
struct Ledger {
    double physical_electron_number=0,particle_charge_C=0,electrode_charge_C=0;
    double wall_heat_J=0,outgoing_K_J=0,incoming_K_J=0;
    double electron_reference_change_J=0,image_reference_change_J=0;
    double resolved_matching_field_work_J=0;
    double residual_J() const;
};
struct Packet {
    std::uint64_t event_id=0,accepted_state=0,patch_id=0;
    double weight=0,normal_metal_energy_eV=0,transverse_energy_eV=0;
    Vector proper_velocity_m_s{};
    Ledger ledger;
};
// Caller owns all draws, event IDs, geometric position sampling and atomic commit.
// FE particles never pass through the old half-Maxwellian sampler.
Packet sample(const PreparedFlux&,const MatchingPlane&,double physical_weight,
              double normal_open_uniform,double transverse_open_uniform,double azimuth_open_uniform,
              std::uint64_t event_id);
struct ConditionalBoundary {
    // Explicit conditional model: locally planar, stationary SN barrier,
    // asymptotic image-tail restoration before PIC-only collective transport.
    double local_linear_vacuum_extent_m=0; // computed from actual field/geometry cell
    double source_log_current_timescale_s=0; // from accepted source variation, or declared steady max-double
    double image_force_fraction=.03,image_energy_fraction=.03;
    double maximum_acceleration_time_fraction=.05;
    bool PIC_already_applies_microscopic_single_electron_image_force=false;
    double maximum_straight_projection_gyroangle=.05;
};
struct ConditionalPacket {
    Packet particle;
    Vector outward_displacement_m{},boundary_charge_displacement_C_m{};
    double matching_distance_m=0,classical_outer_barrier_distance_m=0;
    double image_force_fraction=0,image_tail_restored_eV=0;
    double acceleration_timescale_s=0,source_variation_fraction=0;
    double projected_gyroangle=0;
};
// Usable conditional subgrid source. All distances follow a calculated barrier/
// image-tail criterion. Boundary current and work must be deposited once by
// the parent, together with the electron/electrode charge and heat transaction.
ConditionalPacket sample_conditional(const PreparedFlux&,const ConditionalBoundary&,double physical_weight,
    double normal_open_uniform,double transverse_open_uniform,double azimuth_open_uniform,std::uint64_t event_id);
double source_log_current_timescale_s(const Distribution& previous,const Distribution& current,double dt_s);
Ledger collect_at_material_surface(double work_function_eV,double incoming_K_J,double physical_weight,
                                   bool electronic_accommodation_admitted);
enum class LayerStatus {Active,PICHandedOff,Recollected};
struct LayerPolicy {
    double local_linear_vacuum_extent_m=0; // actual geometry/field-valid region, including lateral motion
    double image_force_fraction=.03,image_energy_fraction=.03;
    double maximum_quantum_magnetic_energy_fraction=.2; // conditional unmagnetized barrier criterion
    double maximum_gyroangle=.02,maximum_relative_normal_step=.02;
    int maximum_substeps=4096;
    bool PIC_already_applies_microscopic_single_electron_image_force=false;
    bool own_return_unit_accommodation=true; // explicit conditional absorbing electronic reservoir
};
struct CurrentPath {
    Vector from_m{},to_m{};double start_time_s=0,end_time_s=0,charge_C=0,field_work_J=0;
};
struct LayerExchange {
    double physical_electron_number_change=0,particle_charge_change_C=0,electrode_charge_change_C=0;
    double wall_heat_J=0,particle_K_change_J=0,electron_reference_change_J=0,image_energy_change_J=0,field_work_J=0;
    double residual_J() const;
};
struct LayerState {
    Packet particle;
    Vector position_from_surface_m{},outward_normal{};
    double current_time_s=0,return_overlap_distance_m=0,handoff_distance_m=0;
    double work_function_eV=0,source_energy_scale_eV=0,quantum_magnetic_energy_fraction=0;
    std::uint64_t mesh_revision=0,material_revision=0;
    LayerStatus status=LayerStatus::Active;
};
struct LayerTrial {
    LayerState state;
    LayerExchange exchange;
    std::vector<CurrentPath> paths;
    int substeps=0;
};
// Birth is a stationary quantum-to-classical overlap projection. Its calculated
// displacement/current/work is explicit; subsequent Lorentz+image flight has
// actual time and path, can remain active across accepted intervals or return.
LayerTrial begin_layer(const PreparedFlux&,const LayerPolicy&,double physical_weight,
    double normal_open_uniform,double transverse_open_uniform,double azimuth_open_uniform,std::uint64_t event_id);
LayerTrial advance_layer(const LayerState&,const Stage& frozen_field_stage,const LayerPolicy&,double duration_s);
Ledger collect_absorbed(const Distribution&,const MatchingPlane&,double incoming_K_J,
                       double physical_weight,bool electronic_accommodation_admitted);
Ledger sum_packets(const std::vector<Packet>&);
double kinetic_energy_J(const Vector& proper_velocity);
}
