#pragma once
#include "surface_emission.hpp"

namespace magnetized_emission {
struct RejectedSupply {
    bool SchottkyNordheim_guard_applies=false;
    double rejected_fraction=0,log_rejected_fraction=0;
    bool has_log_rejected_fraction=false;
    // Fermi-relative (normal+transverse) first moment of rejected phases,
    // weighted by their fraction of the NOMINAL unmagnetized emitted flux.
    double rejected_energy_moment_eV=0;
    double rejected_current_A_m2=0,rejected_complete_wall_heat_W_m2=0;
    double log_rejected_current_A_m2=0;
    bool has_log_rejected_current=false;
};
// Exact analytic transverse integration of the unchanged v1 nominal0.2
// applicability predicate, followed by its normal-spectrum quadrature.
// This quantifies nominal MG supply rejected by that predicate. It does not
// infer actual magnetized emission or authorize dropping rejected phases.
RejectedSupply assess_rejected_supply(const surface_emission::Distribution&,
    const surface_emission::Vector& B_T,const surface_emission::Vector& material_outward_normal,
    double image_energy_fraction=.03);
}
