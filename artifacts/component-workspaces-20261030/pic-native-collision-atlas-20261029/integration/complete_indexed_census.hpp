#pragma once
#include "physical_collision_atlas.hpp"
#include <utility>

namespace raikou::collision_geometry {
template<class Entries>struct CompleteIndexedCensus {
 Entries entries;
 OccupiedCellView occupied;
};
// No chemistry, source-status filter, RNG, pusher, or publication. The caller
// must synchronize the common accepted owner during this bounded read.
template<class NativeArena>auto complete_indexed_census(const NativeArena&arena,
 const CollisionAtlas&atlas,CensusLimits limits={}){
 pc::require(limits.carriers<=1000000&&limits.occupied_cells<=1000000&&limits.string_bytes<=256*1024*1024,"complete census resource caps");
 const auto snapshot=arena.snapshot();pc::require(bool(snapshot),"missing accepted canonical snapshot");
 const auto prefix=arena.state_prefix();CanonicalCensusBinding binding{prefix,{}};
 for(const auto&[id,packet]:snapshot->packets){
  (void)id;pc::require(packet.native_ids.size()<=limits.carriers-binding.expected_native_ids.size(),"canonical native ID census cap");
  binding.expected_native_ids.insert(binding.expected_native_ids.end(),packet.native_ids.begin(),packet.native_ids.end());
 }
 auto entries=arena.rate_population(limits.carriers);
 auto occupied=atlas.occupied_cells(entries,binding,limits);
 pc::require(arena.snapshot()==snapshot&&arena.state_prefix()==prefix,"accepted native owner changed during complete indexed census");
 return CompleteIndexedCensus<decltype(entries)>{std::move(entries),std::move(occupied)};
}
} // namespace raikou::collision_geometry
