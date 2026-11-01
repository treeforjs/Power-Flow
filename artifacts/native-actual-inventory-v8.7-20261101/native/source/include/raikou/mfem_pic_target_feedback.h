#pragma once
#include "raikou/mfem_pic_full_time.h"
#include "chem_event.h" // frozen C3 ABI; this interface is not a portable wire format
namespace raikou::mfem_backend {
// Authentic native-owned immutable inputs. Every preparation creates a fresh
// deep snapshot; public data and digests are facts, never owner authority.
struct FullTimeHostTargetInputs {
  ExposureBinding exposure;
  std::uint64_t stage = 0, native_rows = 0;
  std::string contract_sha256, configuration_sha256, registry_sha256, blob_sha256,
      grid_sha256, native_rows_root, census_sha256, water_scope_binding_sha256,
      frame_snapshot_sha256, frame_history_sha256, carriers_root, targets_root, sha256;
  std::vector<std::uint32_t> atlas;
  std::vector<pfc::CellTargets> cells;
};
struct FullTimeHostTargetDraft {
  std::vector<pfc::CellTargets> cells;
  std::string input_sha256, owner_configuration_sha256, owner_accepted_sha256,
      owner_candidate_sha256, owner_feedback_sha256, before_targets_root, after_targets_root;
};
class FullTimeHostTargetProvider {
public:
  virtual ~FullTimeHostTargetProvider() = default;
  virtual FullTimeHostTargetDraft apply(std::shared_ptr<const FullTimeHostTargetInputs>) = 0;
  // Preparation-time operation; may allocate/throw. No publication authority.
  virtual void validate_current(std::shared_ptr<const FullTimeHostTargetInputs>,
                                const FullTimeHostTargetDraft &) const = 0;
};
struct FullTimeHostTargetReceipt {
  std::shared_ptr<const FullTimeHostTargetInputs> input;
  FullTimeHostTargetDraft output; // native deep copy, sealed after validation
  std::string sha256;
};
std::string full_time_host_target_cells_sha256(std::span<const pfc::CellTargets>);
} // namespace raikou::mfem_backend
