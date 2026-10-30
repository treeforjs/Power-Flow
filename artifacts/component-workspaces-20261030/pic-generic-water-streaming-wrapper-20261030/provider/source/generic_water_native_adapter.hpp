#pragma once
#include "generic_water_bridge.hpp"
#include "raikou/mfem_pic_arena.h"
namespace raikou::generic_water {
namespace m=mfem_backend;
inline std::vector<m::ArenaStateSpec> native_catalog(const Bridge& b){
    std::vector<m::ArenaStateSpec> out;for(const auto&s:b.catalog())out.push_back({s.id,s.H,s.O,s.q,s.mass_kg,s.internal_J,s.reference_J});return out;
}
// Existing fixed-product registry adapter for forced admission comparisons.
// This is NOT a replacement for a mechanism/channel-set native registry.
inline m::ArenaChannel fixed_native_comparison_channel(const Outcome&o,const std::string&binding_file){
    pc::require(pc::digest(binding_file),"unbound native comparison file");std::vector<std::string> states;
    for(const auto&p:o.products)states.push_back(p.state.id);
    return {o.channel,"electron","H2O",states,o.data_sha256,o.kernel_sha256,binding_file,o.model};
}
inline m::ArenaEvent native_event(const Outcome&o,const std::string&event_id,const std::string&primary,
    const std::string&target,m::ArenaCount primary_count,m::ArenaPhase phase,m::ArenaProvenance provenance){
    pc::require(!event_id.empty()&&!primary.empty()&&!target.empty()&&primary!=target,"native outcome identity");
    pc::require(provenance.dataset_sha256==o.data_sha256&&provenance.kernel_manifest_sha256==o.kernel_sha256&&provenance.model==o.model,"native outcome provenance differs");
    m::ArenaEvent e;e.id=event_id;e.channel=o.channel;e.primary=primary;e.target=target;e.event_weight=o.event_weight;
    if(!o.captures_primary)e.scattered_primary_id=primary;
    if(!m::PicArena::positive_carriers(m::PicArena::debit(primary_count,o.event_weight)).empty())e.untouched_primary_id=event_id+"/untouched-primary";
    e.phase=std::move(phase);e.provenance=std::move(provenance);
    for(std::size_t i=0;i<o.products.size();++i){const auto&p=o.products[i];e.products.push_back({p.at_primary?primary:event_id+"/product/"+std::to_string(i),p.state.id,p.proper_velocity_m_s,p.at_primary?m::ArenaProduct::At::Primary:m::ArenaProduct::At::Target});}
    return e;
}
inline Json mechanism_contract(const Bridge& b,const std::string&binding_file){
    pc::require(pc::digest(binding_file),"unbound mechanism binding file");
    Json channels=Json::array();for(const auto&s:b.specs())channels.push_back({{"channel",s.channel},{"model",s.model},{"dataset_sha256",b.dataset_sha256()},{"kernel_manifest_sha256",b.kernel_sha256()},{"kernel_binding_file_sha256",binding_file},{"data_support_eV",{s.low_eV,s.high_eV}},{"captures_primary",s.captures_primary},{"allowed_product_multisets",s.allowed_product_states},{"rate_multiplicity",1},{"product_variates",s.variates}});
    Json states=Json::array();for(const auto&s:b.catalog())states.push_back({{"state",s.id},{"H",s.H},{"O",s.O},{"charge_e",s.q},{"mass_kg",s.mass_kg},{"internal_J",s.internal_J},{"reference_J",s.reference_J}});
    return {{"schema","generic-water-comparison-mechanism-contract-v1"},{"complete_mechanism",false},{"competing_SSA_candidate_set",true},{"conditional_OHv_family_one_base_rate",true},{"channels",channels},{"common_state_catalog",states},{"OH_vibrational_family",b.OH_family_metadata()},
      {"native_completion_requirement","immutable mechanism/channel-set digest; exact per-child registered channel/data/kernel/model membership plus existing SSA stage/duration/order/state-prefix checks"},
      {"native_product_family_requirement","Hminus permits exactly one OH_X_v0..4 draw in a single captured-primary event; no base-rate duplication"},
      {"current_native_integration_scope","fixed-product forced admission comparisons only; existing single-channel completion remains strict"},
      {"collision_site_policy","scattered primary stays at primary; other products at target; capture requires exact full-3D reactant coincidence; no relocation or reaction-current construction"}};
}
} // namespace raikou::generic_water
