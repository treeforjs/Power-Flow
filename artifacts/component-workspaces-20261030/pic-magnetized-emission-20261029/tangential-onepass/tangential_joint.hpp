#pragma once
#include "surface_emission.hpp"
#include <vector>

namespace tangential_onepass {
struct Config {
    double work_function_eV;
    int transverse_nodes=16,azimuth_nodes=16,normal_panels=4,action_nodes=32;
    double normal_lower_eV=-1,normal_tail_kT=48,transverse_tail_kT=120;
    double unresolved_relative_budget=1e-8;
    bool continuum_free_electron_reservoir_selected=false;
    bool incoherent_onepass_matching_selected=false;
};
struct Node {double normal_eV,transverse_eV,azimuth_rad,probability,cdf,log_transmission;};
struct Distribution {
    surface_emission::Stage stage;
    Config config;
    double known_component_current_A_m2=0,known_component_log_current_A_m2=0;
    double known_component_mean_metal_energy_eV=0,known_component_complete_heat_W_m2=0;
    double unresolved_current_upper_A_m2=0,unresolved_absolute_heat_upper_W_m2=0;
    double magnetic_reservoir_energy_scale_over_kT=0;
    double quantum_region_extent_bound_m=0;
    bool truncation_bounds_within_declared_budget=false;
    bool physical_foil_model_admitted=false;
    std::size_t unresolved_quadrature_phase_points=0;
    std::vector<Node> known_component_nodes;
};
// Conditional canonical one-pass prototype. It rejects non-normal E,
// non-tangential B, retarding fields and undeclared reservoir assumptions.
// Unknown/tail bounds remain explicit; total physical current is not inferred.
Distribution evaluate(const surface_emission::Stage&,const Config&,
    const surface_emission::LayerPolicy&);
// Samples ONLY the separately named known one-pass component, at its own
// current, never the old nominal current renormalized onto accepted phases.
// No packet rejection/resampling. Native integration remains parent-owned.
surface_emission::LayerTrial sample_known_component(const Distribution&,
    const surface_emission::LayerPolicy&,double physical_weight,double joint_uniform,
    std::uint64_t event_id);
struct QuantumReturn {
    surface_emission::LayerTrial trial;
    bool absorbed=false;
    double transmission=0,log_transmission=0;
    double reconstructed_normal_eV=0,reconstructed_transverse_eV=0;
};
// Apply the SAME canonical one-pass transmission to an inward quantum event.
// The library's collected trial is still provisional. Reflection undoes its
// quantum material bridge, preserves the preceding real trajectory, reverses
// normal momentum at overlap and retains Active ownership, without heat/work.
// Caller owns each distinct event draw and rolls back the whole trial on error.
QuantumReturn resolve_quantum_return(const surface_emission::Stage&,const Config&,
    const surface_emission::LayerPolicy&,const surface_emission::LayerTrial& provisional_collected,
    double return_uniform);
}
