#pragma once
#include "raikou/mfem_pic_full_time.h"
#include "axial_joint_v2.hpp"

namespace hall_driver {
// Native identity envelope; owner_audit remains the actual qualified provider
// receipt. Native acceptance does not replace field/heat physics preflight.
inline native::ArenaJson native_motion_envelope(const native::FullTimeMotionBinding&admission,
    const MotionReceipt&receipt){
    const auto&common=admission.common;const auto&binding=receipt.binding();
    hall::require(admission.transport_trial_sha256==receipt.seal_sha256(),"native envelope requires actual provider trial digest");
    hall::require(common.accepted_step==binding.accepted_step&&common.begin_s==binding.begin_s&&
        common.end_s==binding.end_s&&common.dt_s==binding.dt_s&&common.common_clock_sha256==binding.common_clock_sha256&&
        common.mesh_sha256==binding.mesh_sha256&&admission.field_before_sha256==common.accepted_field_sha256,
        "native envelope differs from provider common stage");
    return {{"schema","arena_total_field_motion_receipt_v1"},
        {"transport_trial_sha256",admission.transport_trial_sha256},
        {"common_clock_sha256",common.common_clock_sha256},{"mesh_sha256",common.mesh_sha256},
        {"field_before_sha256",admission.field_before_sha256},{"field_after_sha256",admission.field_after_sha256},
        {"begin_s",common.begin_s},{"end_s",common.end_s},{"dt_s",common.dt_s},
        {"owner_audit",motion_record(receipt)}};
}
// A transient proposal owned by the root common transaction. No accepted
// particle population, parent RNG or source state is held by this helper.
struct PreparedNativeJoint {
    AxialJointTrial field_heat;
    native::PreparedPrivateMotion native_motion;
    native::FullTimeMotionBinding admission;
    std::vector<native::ArenaRateEntry> before_census;
};
// Root owns first/second SSA exposure, outer trial selection and common commit.
// This helper fills all Hall-owned motion/field/heat seals and calls the native
// owner's ACTUAL typed motion admission on its private transport trial.
inline PreparedNativeJoint prepare_native_joint(native::PrivateTransportTrial&transport,
    const native::FullTimeStepBinding&common,const Mesh&mesh,const axial::Fields&fields,
    const axial::State&accepted_field,const native::GapPicOptions&kinetic_geometry,
    Thermal&thermal,const CurrentTrace&current,const std::string&species_sha256,
    const std::string&inactive_source_sha256,std::size_t maximum_carriers){
    hall::require(common.accepted_field_sha256==field_state_sha256(accepted_field),"native joint stale accepted field binding");
    auto census=transport.native_census(maximum_carriers);
    MotionBinding binding;binding.phase=ParticlePhase::full_time_position_and_proper_velocity;
    binding.accepted_step=common.accepted_step;binding.staged_arena_revision=transport.staged_arena_revision();
    binding.begin_s=common.begin_s;binding.end_s=common.end_s;binding.dt_s=common.dt_s;
    binding.staged_arena_sha256=transport.state_prefix();binding.mesh_sha256=common.mesh_sha256;
    binding.material_sha256=common.material_sha256;binding.species_sha256=species_sha256;
    binding.common_clock_sha256=common.common_clock_sha256;
    binding.kinetic_policy_sha256=kinetic_policy_sha256(declared_kinetic_policy(kinetic_geometry));
    auto joint=prepare_axial_joint(mesh,fields,accepted_field,kinetic_geometry,census,binding,thermal,current);
    native::FullTimeMotionBinding admission;admission.common=common;
    admission.transport_trial_sha256=joint.motion->seal_sha256();
    admission.field_before_sha256=common.accepted_field_sha256;admission.field_after_sha256=field_state_sha256(joint.field);
    admission.heat_after_sha256=thermal_trial_state_sha256(joint.heat);admission.source_after_sha256=inactive_source_sha256;
    std::vector<native::GapParticle>endpoints(joint.motion->endpoints().begin(),joint.motion->endpoints().end());
    auto native_motion=transport.prepare_total_field_motion(admission,endpoints,native_motion_envelope(admission,*joint.motion));
    return {std::move(joint),std::move(native_motion),std::move(admission),std::move(census)};
}
inline bool preflight_native_joint(const axial::State&accepted_field,const Thermal&thermal,const PreparedNativeJoint&p){
    return preflight_axial_joint(accepted_field,thermal,p.before_census,p.field_heat);
}
} // namespace hall_driver
