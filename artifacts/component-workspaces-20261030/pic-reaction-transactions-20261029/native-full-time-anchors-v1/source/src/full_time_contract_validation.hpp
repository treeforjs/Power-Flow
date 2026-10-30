#pragma once
#include "raikou/mfem_pic_full_time.h"
#include <limits>
#include <stdexcept>
namespace raikou::mfem_backend::full_time_detail {
// Pure contracts shared by live admission, serialized replay and component
// tests.
inline void
completion_chain_guard(const ExposureBinding &e,
                       const ArenaStochasticCompletion &c,
                       const FullTimeAcceptedBoundary &accepted,
                       const std::optional<ArenaStochasticCompletion> &first) {
  auto need = [](bool ok, const char *message) {
    if (!ok)
      throw std::invalid_argument(message);
  };
  need(c.common_clock_sha256 == e.common.common_clock_sha256 &&
           c.stage_sequence == e.owner_stage_sequence &&
           c.common_interval_begin_s == e.common.begin_s &&
           c.common_interval_end_s == e.common.end_s &&
           c.frozen_physical_stage_s == e.frozen_physical_stage_s,
       "native completion exact common exposure binding");
  if (e.slot == FullTimeSlot::before_transport)
    need(!first && c.before_rng_clock_sha256 == accepted.rng_sha256,
         "first native completion differs from accepted chemistry RNG seal");
  else
    need(
        first && c.integration_owner_id == first->integration_owner_id &&
            c.channel_id == first->channel_id &&
            c.dataset_sha256 == first->dataset_sha256 &&
            c.kernel_manifest_sha256 == first->kernel_manifest_sha256 &&
            c.kernel_binding_file_sha256 == first->kernel_binding_file_sha256 &&
            c.model == first->model &&
            first->stage_sequence < std::numeric_limits<std::uint64_t>::max() &&
            c.stage_sequence == first->stage_sequence + 1 &&
            c.before_rng_clock_sha256 == first->after_rng_clock_sha256,
        "native second completion must continue same first "
        "owner/channel/data/kernel/model/RNG");
}
inline void final_rng_guard(const FullTimeCommitSeals &seals,
                            const ArenaStochasticCompletion &second) {
  if (seals.rng_sha256 != second.after_rng_clock_sha256)
    throw std::invalid_argument(
        "final native RNG seal differs from actual second completion");
}
} // namespace raikou::mfem_backend::full_time_detail
