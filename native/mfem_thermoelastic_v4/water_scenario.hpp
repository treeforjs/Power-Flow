#pragma once
#include "electrode_heat_feedback.hpp"
namespace electrode_heat_v2 {
// A conditional ideal Fe2O3(0001)-like oxide surrogate, not measured SS304 kinetics.
// No default physical inventory, mass, chemical state identity or regime selection.
enum class WaterRange {Reviewed300To1300,ConditionalAnalytic298To3000};
inline Patch oxide_water_patch(const Geometry& g,std::string id,Face face,double begin,double end,double assumed_initial_coverage,double parent_registry_mass_kg,std::string gas_state_id,std::string parent_id,std::string scenario_hash,WaterRange range) {
 need(range==WaterRange::Reviewed300To1300||range==WaterRange::ConditionalAnalytic298To3000,"explicit supported water-range selection required");
 need(finite(assumed_initial_coverage)&&assumed_initial_coverage>0&&assumed_initial_coverage<=1,"explicit submonolayer initial coverage required");
 Patch p;p.id=std::move(id);p.face=face;p.begin=begin;p.end=end;p.desorption_enabled=true;p.state_id=std::move(gas_state_id);p.parent_id=std::move(parent_id);p.source_hash=std::move(scenario_hash);p.mass_kg=parent_registry_mass_kg;p.coverage_temkin=true;p.monolayer_m2=1.3e19;p.rate_alpha=.5;p.binding_alpha=.5;p.activation_eV=.87;p.desorption_eV=.87;p.internal_eV=0;
 // Canonical Eyring kBT/h, kappa=1; chosen model, not a fitted material prefactor.
 p.attempt_s=kb*300./6.62607015e-34;p.prefactor=Prefactor::EyringScaled;p.prefactor_reference_T=300;p.source_T_min=range==WaterRange::Reviewed300To1300?300:298.2;p.source_T_max=range==WaterRange::Reviewed300To1300?1300:3000;
 const double A=(end-begin)*(face==Face::LowY||face==Face::HighY?g.width:g.thickness)*g.height;p.initial_particles=assumed_initial_coverage*p.monolayer_m2*A;return p;
}
}
