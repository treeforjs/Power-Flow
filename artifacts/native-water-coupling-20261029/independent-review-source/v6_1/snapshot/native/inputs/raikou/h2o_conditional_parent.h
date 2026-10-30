#pragma once
// Explicit conditional kinematic comparison model, not measured water branching.
#include "raikou/collision_data.h"
#include <array>

namespace raikou::kinetics {
using CollisionMomentum=std::array<double,3>;
struct ConditionalParentProducts {
    CollisionMomentum primary_momentum_kg_m_s{},secondary_momentum_kg_m_s{},ion_momentum_kg_m_s{};
    double primary_kinetic_J{0},secondary_kinetic_J{0},ion_kinetic_J{0};
    double ion_mass_kg{0},chemical_reference_gain_J{0},physical_weight_each{0},cross_section_m2{0};
};
[[nodiscard]] inline ConditionalParentProducts conditional_h2o_parent_products(
    const StrictCollisionDataset& data,const std::string& explicitly_selected_model,
    double incident_target_rest_eV,double physical_weight,double neutral_mass_kg,
    const CollisionMomentum& incident_unit_direction,const CollisionMomentum& isotropic_CM_pair_unit_axis) {
    constexpr double e=1.602176634e-19,me=9.1093837139e-31,c=299792458.,I=12.621;
    const auto require=[](bool value,const char* text) {if (!value) throw std::invalid_argument(text);};
    require(explicitly_selected_model=="parent_ground_equal_pair_CM_v1","Explicit conditional parent model opt-in required");
    data.require_MCC_event_kernel("H2O:ion_yield:H2O+",explicitly_selected_model,incident_target_rest_eV);
    require(std::isfinite(physical_weight)&&physical_weight>0,"Invalid represented event weight");
    require(std::isfinite(neutral_mass_kg)&&neutral_mass_kg>1000*me,"Invalid neutral molecular mass");
    for (const auto& axis:{incident_unit_direction,isotropic_CM_pair_unit_axis}) {
        double norm=0;for (const auto v:axis) {require(std::isfinite(v),"Nonfinite model axis");norm+=v*v;}
        require(std::abs(norm-1)<1e-12,"Model axes must already be unit vectors from the host");
    }
    // This checks the original parent-only cross-section domain, including its
    // onset gap. The BEB/total or partial-yield curves cannot be substituted.
    const double sigma=data.observable_value_SI("H2O:ion_yield:H2O+",incident_target_rest_eV);
    require(sigma>0,"No conditional parent event below its supported onset");
    const double K=incident_target_rest_eV*e,D=(neutral_mass_kg+me)*c*c;
    const double root_s=std::sqrt(D*D+2*neutral_mass_kg*c*c*K);
    const double available=2*neutral_mass_kg*c*c*K/(root_s+D),Q=available-I*e;
    require(Q>0&&std::isfinite(Q),"Insufficient invariant energy for parent model");
    const double p0=std::sqrt(K*(K+2*me*c*c))/c;
    const double beta=p0*c/(D+K),gamma=1/std::sqrt(1-beta*beta);
    const double gamma_minus_one=beta*beta*gamma*gamma/(gamma+1);
    const double each_K=Q/2,each_E=me*c*c+each_K;
    const double each_p=std::sqrt(each_K*(each_K+2*me*c*c))/c;
    double dot=0;for (int d=0;d<3;++d)dot+=isotropic_CM_pair_unit_axis[d]*incident_unit_direction[d];
    ConditionalParentProducts out;out.ion_mass_kg=neutral_mass_kg-me+I*e/(c*c);
    out.chemical_reference_gain_J=I*e;out.physical_weight_each=physical_weight;out.cross_section_m2=sigma;
    out.primary_kinetic_J=each_K+gamma_minus_one*each_E+gamma*beta*c*each_p*dot;
    out.secondary_kinetic_J=each_K+gamma_minus_one*each_E-gamma*beta*c*each_p*dot;
    out.ion_kinetic_J=gamma_minus_one*out.ion_mass_kg*c*c;
    for (int d=0;d<3;++d) {
        out.primary_momentum_kg_m_s[d]=each_p*isotropic_CM_pair_unit_axis[d]+
            (gamma_minus_one*each_p*dot+gamma*beta*each_E/c)*incident_unit_direction[d];
        out.secondary_momentum_kg_m_s[d]=-each_p*isotropic_CM_pair_unit_axis[d]+
            (-gamma_minus_one*each_p*dot+gamma*beta*each_E/c)*incident_unit_direction[d];
        out.ion_momentum_kg_m_s[d]=gamma*beta*out.ion_mass_kg*c*incident_unit_direction[d];
    }
    return out;
}
} // namespace raikou::kinetics
