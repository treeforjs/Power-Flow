#pragma once
#include "raikou/mfem_pic_full_time.h"
#include "native_joint_binding.hpp"

namespace hall_driver {
// Minimal integration example using the native owner's actual DRAFT signatures
// and root ParentSubflow::prepare_on/seal API. This file is a caller example,
// not a second arena implementation or a claim of linked driver qualification.
// The caller must provide an already-authored FULL-TIME genesis and the exact
// current source/species/common-clock identities and parent interval sequence.
// Final native/parent tokens MUST certify exact first/second private identities,
// transport trial identity and accepted genesis/RNG seal before using this flow.
template<class ParentSubflow>
void root_two_exposure_example(native::PicArena&arena,ParentSubflow&parent,
    const Mesh&mesh,const axial::Fields&fields,axial::State&accepted_field,
    Thermal&thermal,const CurrentTrace&current,double dt,
    const std::string&species_sha256,const std::string&inactive_source_sha256,
    const std::string&common_clock_sha256,std::uint64_t first_owner_sequence,
    std::size_t maximum_carriers,std::size_t maximum_records){
    hall::require(first_owner_sequence<std::numeric_limits<std::uint64_t>::max()-1,"two exposure sequence exhaustion");
    native::FullTimeStepBinding common;
    common.accepted_step=arena.snapshot()->accepted_steps;common.accepted_arena_revision=arena.snapshot()->revision;
    common.begin_s=accepted_field.time_s;common.end_s=common.begin_s+dt;common.dt_s=dt;
    common.accepted_arena_sha256=arena.state_prefix();common.accepted_field_sha256=field_state_sha256(accepted_field);
    common.mesh_sha256=motion_mesh_sha256(mesh,arena.native_options());common.material_sha256=material_binding_sha256();
    common.common_clock_sha256=common_clock_sha256;
    auto outer=arena.begin_full_time_step(common,maximum_carriers,maximum_records);
    native::ExposureBinding before{common,native::FullTimeSlot::before_transport,common.begin_s,dt/2,
        first_owner_sequence,common.accepted_field_sha256};
    auto first_scope=outer.first_exposure(before);
    auto first=parent.prepare_on(first_scope);
    auto transport=outer.fork_transport_trial();
    auto proposal=prepare_native_joint(transport,common,mesh,fields,accepted_field,arena.native_options(),
        thermal,current,species_sha256,inactive_source_sha256,maximum_carriers);
    outer.accept_transport_trial(std::move(proposal.native_motion));
    const auto&admission=proposal.admission;
    native::ExposureBinding after{common,native::FullTimeSlot::after_transport,common.end_s,dt/2,
        first_owner_sequence+1,admission.field_after_sha256};
    auto second_scope=outer.second_exposure(after);
    auto second=parent.prepare_on(second_scope,&first); // BOTH staged RNG streams
    native::FullTimeCommitSeals seals{admission.field_after_sha256,admission.heat_after_sha256,
        inactive_source_sha256,parent.staged_rng_sha256(second)};
    auto native_token=outer.seal(seals);auto parent_token=parent.seal(first,second,native_token);
    hall::require(preflight_native_joint(accepted_field,thermal,proposal)&&
        arena.can_publish(native_token)&&parent.can_publish(parent_token,native_token),"all-owner preflight");
    // No fallible computation, allocations, callbacks or I/O below this point.
    if(!arena.try_publish(native_token))std::terminate();
    parent.publish_prevalidated(parent_token);
    publish_field_heat_prevalidated(accepted_field,thermal,std::move(proposal.field_heat));
}
} // namespace hall_driver
