#pragma once
#include "surface_emission.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace surface_emission_integration {
// Pure inverse of usable-api-v1's outward image-tail reference restoration.
// Caller resolves the first inward handoff-plane crossing in ordinary PIC,
// commits this exchange once, then advances the returned Active LayerState
// over the remaining common interval. This transition has no displacement.
inline surface_emission::LayerTrial reenter_image_layer(
    const surface_emission::LayerState& initial,
    const surface_emission::Stage& field,
    const surface_emission::LayerPolicy& policy) {
    namespace se=surface_emission;
    const auto require=[](bool b,const char* text){if(!b)throw std::invalid_argument(text);};
    const auto dot=[](const se::Vector& a,const se::Vector& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];};
    require(initial.status==se::LayerStatus::PICHandedOff,"Re-entry requires a PIC-owned handoff state; no duplicate restore");
    require(field.accepted_state>0 && field.accepted_state==field.field_accepted_state &&
        field.accepted_state==field.thermal_accepted_state && field.field_time==field.stage_time &&
        field.thermal_time==field.stage_time && field.mesh_revision==initial.mesh_revision &&
        field.material_revision==initial.material_revision && field.patch_id==initial.particle.patch_id &&
        field.outward_normal==initial.outward_normal,"Image-layer re-entry field/material/stage mismatch");
    require(std::isfinite(initial.current_time_s) && initial.current_time_s>=field.t &&
        initial.current_time_s<=field.t+field.dt+1e-24,"Re-entry exceeds common accepted field interval");
    const double x=dot(initial.position_from_surface_m,initial.outward_normal),un=dot(initial.particle.proper_velocity_m_s,initial.outward_normal);
    require(x>0 && std::isfinite(x) && initial.handoff_distance_m>0 &&
        x<=initial.handoff_distance_m*(1+1e-12) && x>=initial.return_overlap_distance_m &&
        std::sqrt(dot(initial.position_from_surface_m,initial.position_from_surface_m))<=policy.local_linear_vacuum_extent_m &&
        un<0 && std::isfinite(un) && !policy.PIC_already_applies_microscopic_single_electron_image_force,
        "Re-entry must resolve inward crossing in the admitted local image region");
    constexpr double pi=3.14159265358979323846;
    const double image=-se::qe*se::qe/(16*pi*se::eps0*x),K=se::kinetic_energy_J(initial.particle.proper_velocity_m_s),target=K-image;
    const double perpendicular2=std::max(0.,dot(initial.particle.proper_velocity_m_s,initial.particle.proper_velocity_m_s)-un*un);
    const double total_u2=target*(target+2*se::me*se::c_light*se::c_light)/(se::me*se::me*se::c_light*se::c_light);
    require(std::isfinite(total_u2) && total_u2>=perpendicular2 && initial.particle.weight>0,
        "Re-entry kinetic target is not representable");
    se::LayerTrial trial;trial.state=initial;
    const double normal_change=-std::sqrt(total_u2-perpendicular2)-un;
    for(int i=0;i<3;++i)trial.state.particle.proper_velocity_m_s[i]+=normal_change*initial.outward_normal[i];
    trial.state.status=se::LayerStatus::Active;
    const double weight=initial.particle.weight;
    trial.exchange.particle_K_change_J=-image*weight;
    trial.exchange.image_energy_change_J=image*weight;
    require(std::abs(se::kinetic_energy_J(trial.state.particle.proper_velocity_m_s)-target)<1e-11*std::max(target,1e-30),
        "Image re-entry momentum failed energy closure");
    return trial;
}
}
