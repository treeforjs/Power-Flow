#pragma once
// Common-commit predicate v7: the V3 closed-step contract
// (pic-water-common-commit-v3-20261030/source/common_commit_contract.hpp) and the
// native v6 birth-step admission (candidate006 adapters/common_commit_v6.hpp)
// with a variable step length. V3 hard-codes 0 < dt <= 1e-14; v7 admits
// 0 < dt <= cap, where cap is the configured step cap of the run, itself bounded
// by the largest step the coupled owners admit (Hall: 1e-14 for v3, Model::dt_max() for
// v3.1; heat/source owner: 1e-14), and dt
// must be bit-identical to the step controller's decision of the step. Every
// other condition is the V3/native-v6 condition, in the same order, unchanged.
// With cap = 1e-14 and dt <= 1e-14 the decision equals the V3/v6 predicate's.
#include "common_commit_contract.hpp"
#include "raikou/mfem_pic_full_time.h"
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>

namespace raikou::native_water::common_commit_v7 {
namespace cc = ::raikou::native_water::common_commit;
namespace nm = ::raikou::mfem_backend;

// Clock admission of one step: the configured cap, the largest step the coupled owners
// admit (field owner cap and the other owners' bounds) and the controller decision
// (dt bits and the decision record digest).
struct StepClock {
  double cap_s = 0, owner_max_s = 0, decided_dt_s = 0;
  std::string decision_sha256;
};

inline bool admits_clock(const cc::ExpectedStep& e, const StepClock& c) noexcept {
  return std::isfinite(c.owner_max_s) && c.owner_max_s > 0 && std::isfinite(c.cap_s) && c.cap_s > 0 &&
         c.cap_s <= c.owner_max_s && cc::digest(c.decision_sha256) &&
         std::bit_cast<std::uint64_t>(c.decided_dt_s) == std::bit_cast<std::uint64_t>(e.dt_s) && std::isfinite(e.dt_s) && e.dt_s > 0 &&
         e.dt_s <= c.cap_s && .5 * e.dt_s > 0;
}

// V3 accepts_closed_step with dt <= cap (StepClock) in place of dt <= 1e-14.
inline bool accepts_closed_step(const cc::ExpectedStep& e, const cc::ObservedOwners& current, const cc::ProposedOwners& candidate,
                                const StepClock& clock) noexcept {
  const auto& b = e.before;
  const auto& a = candidate.after.native;
  if (!cc::valid_boundary(b) || !cc::owner_values_match(current) || !cc::same_boundary(b, current.native) ||
      !cc::owner_values_match(candidate.after))
    return false;
  if (!admits_clock(e, clock) || !std::isfinite(e.end_s) || e.end_s != b.time_s + e.dt_s || e.end_s <= b.time_s || !cc::digest(e.common_clock) ||
      !cc::digest(e.configuration) || current.configuration != e.configuration || candidate.after.configuration != e.configuration)
    return false;
  constexpr auto max = std::numeric_limits<std::uint64_t>::max();
  if (b.accepted_step >= max - 1 || b.accepted_arena_revision >= max - 1 || a.accepted_step != b.accepted_step + 1 ||
      a.accepted_arena_revision <= b.accepted_arena_revision || a.accepted_arena_revision == max || a.time_s != e.end_s ||
      a.authored_full_time_genesis)
    return false;
  if (a.common_clock_sha256 != e.common_clock || a.mesh_sha256 != b.mesh_sha256 || a.material_sha256 != b.material_sha256 ||
      a.frame_sha256 != b.frame_sha256 || a.source_sha256 != b.source_sha256)
    return false;
  return cc::digest(candidate.motion_midpoint) && candidate.motion_midpoint == candidate.field_midpoint;
}

// Arrival owner evidence of the step (as native v6 common_commit_v6::ArrivalEvidence).
struct ArrivalEvidence {
  bool staged = false;
  std::string plan_sha256, token_sha256;
  std::uint64_t base = 0, births = 0;
};

// Native v6 accepts_closed_or_birth_step with the v7 closed-step predicate.
inline bool accepts_closed_or_birth_step(const cc::ExpectedStep& e, const cc::ObservedOwners& current, const cc::ProposedOwners& candidate,
                                         const nm::FullTimeStepFlows* flows, const ArrivalEvidence& arrival, const StepClock& clock) noexcept {
  try {
    const auto& b = e.before;
    const auto& a = candidate.after.native;
    const nm::FullTimeBirthReceipt* births = flows && flows->births ? flows->births.get() : nullptr;
    if (!cc::digest(candidate.after.source) || candidate.after.source != a.source_sha256) return false;
    if (births) {
      if (births->accepted_step != b.accepted_step || births->source_before_sha256 != b.source_sha256 ||
          births->source_after_sha256 != a.source_sha256 || births->next_id_after < births->base ||
          births->next_id_after - births->base != births->entries.size())
        return false;
      if (!arrival.staged || arrival.plan_sha256 != births->plan_sha256 || arrival.token_sha256 != births->staged_token_sha256 ||
          arrival.base != births->base || arrival.births != births->entries.size())
        return false;
    } else if (arrival.staged || a.source_sha256 != b.source_sha256)
      return false;
    auto same = candidate;
    same.after.source = b.source_sha256;
    same.after.native.source_sha256 = b.source_sha256;
    return accepts_closed_step(e, current, same, clock);
  } catch (...) {
    return false;
  }
}
}  // namespace raikou::native_water::common_commit_v7
