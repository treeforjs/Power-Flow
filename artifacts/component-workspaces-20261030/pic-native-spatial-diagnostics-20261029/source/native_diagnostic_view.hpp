#pragma once
#include "raikou/mfem_pic_arena.h"

namespace native_diagnostics {
using namespace raikou::mfem_backend;
// Caller supplies its bound physical geometry/state registry and freezes the
// accepted owner for this read. This function never selects or publishes events.
inline ArenaJson export_view(const PicArena&arena,const std::vector<ArenaStateSpec>&states,
 ArenaJson geometry,ArenaJson energy_edges,double physical_time_s,
 std::size_t maximum_entries,std::size_t maximum_sink_receipts){
 const auto prefix=arena.state_prefix();const auto snapshot=arena.snapshot();
 if(maximum_entries>2000000||maximum_sink_receipts>100000||snapshot->sinks.size()>maximum_sink_receipts)throw std::invalid_argument("Diagnostic census/sink cap");
 const auto entries=arena.rate_population(maximum_entries);ArenaJson rows=ArenaJson::array(),registry=ArenaJson::array(),packets=ArenaJson::object();
 for(const auto&s:states)registry.push_back({{"id",s.id},{"hydrogen",s.hydrogen},{"oxygen",s.oxygen},{"charge_e",s.charge_e},{"mass_kg",s.mass_kg},{"internal_J",s.internal_J},{"reference_J",s.reference_J}});
 for(const auto&[id,p]:snapshot->packets)packets[id]={{"state",p.state},{"count",{{"hi",p.count.hi},{"lo",p.count.lo}}},{"native_ids",p.native_ids},{"origin",p.origin},{"last_event",p.last_event}};
 for(const auto&e:entries){const auto&n=e.native;rows.push_back({{"packet_id",e.packet_id},{"state",e.state},{"kinetic_J_per_particle",pic_kinetic_energy_J(n)},
   {"native",{{"id",n.id},{"position_m",n.position_m},{"proper_velocity_m_s",n.proper_velocity_m_s},{"charge_C",n.charge_C},{"mass_kg",n.mass_kg},{"weight",n.weight}}}});}
 if(arena.state_prefix()!=prefix||arena.snapshot()!=snapshot)throw std::runtime_error("Arena changed during diagnostic census");
 return {{"schema","native_spatial_diagnostics_input/v1"},{"state_prefix",prefix},{"snapshot_chain",snapshot->chain},{"revision",snapshot->revision},{"time_s",physical_time_s},{"maximum_entries",maximum_entries},{"arena_v5_source_identity",PicArena::source_identity()},
  {"states",registry},{"canonical_packets",packets},{"entries",rows},{"sinks",snapshot->sinks},{"geometry",geometry},{"energy_edges_eV",energy_edges},{"full_mechanism_qualification",false}};
}
} // namespace native_diagnostics
