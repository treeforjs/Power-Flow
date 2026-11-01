#pragma once
#include "census_charge.hpp"
#include <algorithm>
#include <stdexcept>

// Root-side numerical predicate; it does not mint a publication capability.
// The coordinator must obtain both arguments from its live authenticated Hall
// proposal and actual post-stage2 native candidate, retain their identities,
// and recheck those capabilities immediately before no-fail publication.
namespace powerflow::final_charge_guard_v1 {
namespace h=::hall_corrective_v3;
namespace n=::raikou::mfem_backend;
struct Receipt {
  std::string final_native_census_sha256;
  std::size_t final_carriers=0;
  double maximum_admitted_node_difference_C=0;
};
inline Receipt require_consistent(const axial::Mesh &mesh,const n::GapPicOptions &options,
                                 std::span<const n::GapParticle> motion_endpoints,
                                 std::span<const n::ArenaRateEntry> final_native_census,
                                 std::size_t maximum_carriers) {
  if(maximum_carriers==0)throw std::invalid_argument("final charge: positive carrier cap required");
  const auto expected=h::charge_of_endpoints(mesh,options,motion_endpoints,maximum_carriers);
  const auto actual=h::charge_of_census(mesh,options,final_native_census,maximum_carriers);
  // Reuse Hall's existing per-node gross-charge rounding bound exactly. Global
  // charge conservation cannot substitute for this local Maxwell constraint.
  if(!h::same_net_charge(expected,actual))
    throw std::invalid_argument("final charge: post-chemistry native charge differs from Hall motion endpoints; reaction current/work is required");
  Receipt receipt;receipt.final_native_census_sha256=h::canonical_census_sha256(final_native_census,maximum_carriers);
  receipt.final_carriers=final_native_census.size();
  for(std::size_t k=0;k<expected.net_C.size();++k)
    receipt.maximum_admitted_node_difference_C=std::max(receipt.maximum_admitted_node_difference_C,
      double(std::abs((actual.net_C[k]-expected.net_C[k]).wide())));
  return receipt;
}
} // namespace powerflow::final_charge_guard_v1
