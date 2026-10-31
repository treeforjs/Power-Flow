#pragma once
// Native v6 common predicate for coordinators (COORDINATOR-V6-LINK.txt
// section C). Include with -I pic-water-common-commit-v3-20261030/source; the
// V3 contract header is used unchanged. V3 accepts_closed_step requires an
// unchanged source seal; a step with births changes it. This predicate
// accepts exactly that change when the sealed native step carries a birth
// receipt binding the same transition (source_before = accepted seal,
// source_after = candidate seal) and the arrival owner's staged token, plan
// digest, B_s and birth count; every other V3 condition is applied unchanged.
#include "common_commit_contract.hpp"
#include "raikou/mfem_pic_full_time.h"

namespace raikou::native_water::common_commit_v6 {
namespace cc = ::raikou::native_water::common_commit;
namespace nm = ::raikou::mfem_backend;
// Arrival owner evidence of the step (pf_arrival::Staged of the step).
struct ArrivalEvidence {
  bool staged = false;
  std::string plan_sha256, token_sha256;
  std::uint64_t base = 0, births = 0;
};
// flows = PreparedArena::candidate_full_time_flows() of the sealed native
// step (null when the step has no births and no sinks).
inline bool accepts_closed_or_birth_step(const cc::ExpectedStep &e,
                                         const cc::ObservedOwners &current,
                                         const cc::ProposedOwners &candidate,
                                         const nm::FullTimeStepFlows *flows,
                                         const ArrivalEvidence &arrival) noexcept {
  try {
    const auto &b = e.before;
    const auto &a = candidate.after.native;
    const nm::FullTimeBirthReceipt *births =
        flows && flows->births ? flows->births.get() : nullptr;
    if (!cc::digest(candidate.after.source) ||
        candidate.after.source != a.source_sha256)
      return false;
    if (births) {
      if (births->accepted_step != b.accepted_step ||
          births->source_before_sha256 != b.source_sha256 ||
          births->source_after_sha256 != a.source_sha256 ||
          births->next_id_after < births->base ||
          births->next_id_after - births->base != births->entries.size())
        return false;
      if (!arrival.staged || arrival.plan_sha256 != births->plan_sha256 ||
          arrival.token_sha256 != births->staged_token_sha256 ||
          arrival.base != births->base ||
          arrival.births != births->entries.size())
        return false;
    } else if (arrival.staged || a.source_sha256 != b.source_sha256)
      return false;
    // Every V3 condition, with the (bound) source transition substituted.
    auto same = candidate;
    same.after.source = b.source_sha256;
    same.after.native.source_sha256 = b.source_sha256;
    return cc::accepts_closed_step(e, current, same);
  } catch (...) {
    return false;
  }
}
} // namespace raikou::native_water::common_commit_v6
