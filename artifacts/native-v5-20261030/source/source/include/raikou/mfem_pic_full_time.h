#pragma once
#include "raikou/mfem_pic_arena.h"
#include <functional>
#include <memory>
#include <stdexcept>
namespace hall_driver {
class MotionReceipt;
}

namespace raikou::mfem_backend {
enum class FullTimeSlot { before_transport, after_transport };
struct FullTimeStepBinding {
  std::uint64_t accepted_step = 0, accepted_arena_revision = 0;
  double begin_s = 0, end_s = 0, dt_s = 0;
  std::string accepted_arena_sha256, accepted_field_sha256;
  std::string mesh_sha256, material_sha256, common_clock_sha256;
};
struct ExposureBinding {
  FullTimeStepBinding common;
  FullTimeSlot slot = FullTimeSlot::before_transport;
  double frozen_physical_stage_s = 0, collision_duration_s = 0;
  std::uint64_t owner_stage_sequence = 0;
  std::string frozen_field_sha256;
};
struct FullTimeGenesisBinding {
  double time_s = 0;
  std::string origin_provenance, field_sha256, mesh_sha256, material_sha256;
  std::string frame_sha256, common_clock_sha256, heat_sha256, source_sha256,
      rng_sha256;
};
struct FullTimeAuthoredPacket {
  std::string id, state, origin;
  ArenaCount count;
  PicVector position_m{}, proper_velocity_m_s{};
};
struct FullTimeMotionBinding {
  FullTimeStepBinding common;
  std::string transport_trial_sha256, field_before_sha256, field_after_sha256;
  std::string heat_after_sha256, source_after_sha256;
};
// rng_sha256 binds chemistry scheduler/product RNG + configuration only.
// Future source arrival RNG/count/layer state belongs in source_sha256.
struct FullTimeCommitSeals {
  std::string field_sha256, heat_sha256, source_sha256, rng_sha256;
};
class FullTimeExposureKey;
struct FullTimeExposureReceipt {
  ExposureBinding exposure;
  ArenaStochasticCompletion completion;
  std::string completed_state_prefix, completion_sha256;
  // Live-only opaque identity; not serialized as an accepted checkpoint key.
  std::shared_ptr<const FullTimeExposureKey> key;
};
struct FullTimeCandidateExposures {
  FullTimeStepBinding common;
  std::string outer_step_sha256, transport_trial_sha256;
  FullTimeExposureReceipt first, second;
};
struct FullTimeAcceptedBoundary {
  std::uint64_t accepted_step = 0, accepted_arena_revision = 0;
  double time_s = 0;
  std::string accepted_arena_sha256, field_sha256, heat_sha256, source_sha256,
      rng_sha256;
  std::string mesh_sha256, material_sha256, frame_sha256, common_clock_sha256;
  bool authored_full_time_genesis = false;
};
// ---- native v5 configuration, bounded history and typed refusals ----
// Typed resource refusal: thrown before any callback/mutation when an explicit
// cap or budget would be exceeded. History is never truncated silently.
class FullTimeResourceRefusal : public std::invalid_argument {
public:
  enum class Kind {
    carrier_cap,
    memory_budget,
    configuration,
    anchor_count,
    anchor_bytes,
    window_rows,
    window_bytes,
    identity_count,
    step_records,
    restore_deferred
  };
  FullTimeResourceRefusal(Kind k, const char *message)
      : std::invalid_argument(message), kind_(k) {}
  Kind kind() const noexcept { return kind_; }

private:
  Kind kind_;
};
// Resident carrier cap and the stated host-memory budget it is checked
// against. Both are bound into the native configuration identity (initial
// prefix), so different caps give different accepted-state identities.
struct FullTimeNativeLimits {
  std::size_t maximum_carriers = 100000;
  std::uint64_t memory_budget_bytes = 2ull << 30;
  // Logical/identity entries (ever packet IDs + event IDs) retained exactly.
  std::size_t maximum_identity_entries = 1u << 22;
};
// Anchor every K accepted steps: the validated window is compacted into one
// chained commitment. K is recorded in every anchor and in the policy digest,
// but NOT in accepted-state identities (any K gives identical boundaries).
struct FullTimeHistoryPolicy {
  std::uint64_t anchor_interval_steps = 1024;
  std::size_t maximum_anchors = 1u << 16;
  std::size_t maximum_anchor_bytes = 64u << 20;
  std::size_t maximum_window_rows = 1u << 20;
  std::size_t maximum_window_bytes = 256u << 20;
};
struct FullTimeConfiguration {
  FullTimeNativeLimits limits;
  FullTimeHistoryPolicy history;
};
// Conservative per-carrier host model used by the budget check (bytes):
// resident row copies in accepted/outer/frame/trial/motion/token states, packet
// and identity records, census/index entries and digest leaves.
inline constexpr std::uint64_t full_time_bytes_per_carrier = 2048;
inline constexpr std::uint64_t full_time_fixed_bytes = 64ull << 20;
inline constexpr std::size_t full_time_absolute_carrier_cap = 1u << 27;
struct FullTimeHistoryStatus {
  std::uint64_t anchor_interval_steps = 0, anchors = 0, anchor_bytes = 0;
  std::uint64_t window_rows = 0, window_bytes = 0;
  std::uint64_t cumulative_rows = 0, cumulative_bytes = 0;
  std::uint64_t base_step = 0, accepted_step = 0;
  std::string history_policy_sha256, anchor_root_sha256, anchor_head_sha256;
};
// Decoded anchor fields (round trip of the canonical anchor record).
struct FullTimeAnchorFields {
  std::string source_sha256, history_policy_sha256, previous_anchor_sha256;
  std::uint64_t anchor_interval_steps = 0, sequence = 0;
  std::uint64_t first_step = 0, last_step = 0;
  std::uint64_t base_revision = 0, tail_revision = 0;
  std::string base_history_chain, tail_history_chain;
  std::uint64_t window_rows = 0, window_bytes = 0, cumulative_bytes = 0;
  std::uint64_t window_reactions = 0, window_motions = 0;
  std::string accepted_arena_sha256, snapshot_sha256, native_rows_sha256;
  std::string native_ledger_sha256, completions_sha256, event_ids_sha256,
      packet_ids_sha256;
  std::uint64_t next_native_id = 0;
  double previous_dt_s = 0;
  // Accepted boundary at the anchor: arena/field/heat/source/RNG/mesh/
  // material/frame/clock seals, time, step and revision.
  FullTimeAcceptedBoundary boundary;
  std::string sha256; // digest of the canonical bytes
};
FullTimeAnchorFields decode_full_time_anchor(const std::string &canonical);
// Verifies a serialized chain against an independently held root and head:
// tags/lengths, digests, previous links, contiguous step/revision/chain ranges,
// constant policy and interval, cumulative totals. Throws on any mismatch.
std::vector<FullTimeAnchorFields>
verify_full_time_anchor_chain(const std::string &serialized,
                              const std::string &expected_root,
                              const std::string &expected_head);
// Moved-from handles support destruction and move assignment only, except
// publication predicates below, which return false for an empty handle.
class CollisionScope {
public:
  CollisionScope(CollisionScope &&) noexcept;
  CollisionScope &operator=(CollisionScope &&) noexcept;
  ~CollisionScope();
  CollisionScope(const CollisionScope &) = delete;
  CollisionScope &operator=(const CollisionScope &) = delete;
  const ExposureBinding &exposure() const;
  std::shared_ptr<const FullTimeExposureKey> exposure_key() const;
  ArenaPhase phase() const;
  const GapPicOptions &native_options() const noexcept;
  std::shared_ptr<const ArenaSnapshot> snapshot() const;
  std::vector<ArenaRateEntry>
  rate_population(std::size_t maximum_entries) const;
  const ArenaChannel &channel_binding(const std::string &) const;
  std::optional<ArenaStochasticCompletion>
  stochastic_completion(const std::string &) const;
  std::string state_prefix() const;
  PreparedArena prepare(const ArenaEvent &, int manufactured_fault_point = 0);
  PreparedArena prepare_compound_with_completion(
      std::size_t maximum_reactions,
      const std::function<ArenaStochasticCompletion(CollisionScope &)> &);
  bool can_publish(const PreparedArena &) const noexcept;
  bool try_publish(PreparedArena &) noexcept;

private:
  friend class FullTimeAccess;
  friend class PicArena;
  friend class PrivateFullTimeStep;
  struct Impl;
  explicit CollisionScope(std::unique_ptr<Impl>);
  std::unique_ptr<Impl> impl_;
};
class PreparedPrivateMotion {
public:
  PreparedPrivateMotion(PreparedPrivateMotion &&) noexcept;
  PreparedPrivateMotion &operator=(PreparedPrivateMotion &&) noexcept;
  ~PreparedPrivateMotion();
  PreparedPrivateMotion(const PreparedPrivateMotion &) = delete;
  PreparedPrivateMotion &operator=(const PreparedPrivateMotion &) = delete;
  const ArenaJson &audit() const;

private:
  friend class PrivateTransportTrial;
  friend class FullTimeAccess;
  friend class PicArena;
  friend class PrivateFullTimeStep;
  struct Impl;
  explicit PreparedPrivateMotion(std::unique_ptr<Impl>);
  std::unique_ptr<Impl> impl_;
};
// Only a genuine immutable Hall MotionReceipt can create live admission. JSON
// is a restart receipt, not an authority to propose new motion.
class FullTimeMotionAdmission {
public:
  ~FullTimeMotionAdmission() noexcept;
  FullTimeMotionAdmission(FullTimeMotionAdmission &&) noexcept = default;
  FullTimeMotionAdmission &operator=(FullTimeMotionAdmission &&) noexcept;
  FullTimeMotionAdmission(const FullTimeMotionAdmission &) = delete;
  FullTimeMotionAdmission &operator=(const FullTimeMotionAdmission &) = delete;

private:
  friend FullTimeMotionAdmission
  bind_hall_total_field_motion(const hall_driver::MotionReceipt &);
  friend class PrivateTransportTrial;
  explicit FullTimeMotionAdmission(ArenaJson record)
      : record_(std::move(record)) {}
  ArenaJson record_;
};
FullTimeMotionAdmission
bind_hall_total_field_motion(const hall_driver::MotionReceipt &);

// Moved-from trial/outer handles support destruction and move assignment only.
class PrivateTransportTrial {
public:
  PrivateTransportTrial(PrivateTransportTrial &&) noexcept;
  PrivateTransportTrial &operator=(PrivateTransportTrial &&) noexcept;
  ~PrivateTransportTrial();
  PrivateTransportTrial(const PrivateTransportTrial &) = delete;
  PrivateTransportTrial &operator=(const PrivateTransportTrial &) = delete;
  std::vector<ArenaRateEntry> native_census(std::size_t maximum_entries) const;
  std::shared_ptr<const ArenaSnapshot> snapshot() const;
  std::uint64_t staged_arena_revision() const;
  std::string state_prefix() const;
  PreparedPrivateMotion
  prepare_total_field_motion(const FullTimeMotionBinding &,
                             FullTimeMotionAdmission);
  // Compatibility refusal only: JSON cannot authorize live motion.
  PreparedPrivateMotion prepare_total_field_motion(
      const FullTimeMotionBinding &, const std::vector<GapParticle> &endpoints,
      const ArenaJson &transport_receipts,
      const ArenaJson &source_receipts = ArenaJson::array(),
      const ArenaJson &boundary_receipts = ArenaJson::array());

private:
  friend class FullTimeAccess;
  friend class PicArena;
  friend class PrivateFullTimeStep;
  struct Impl;
  explicit PrivateTransportTrial(std::unique_ptr<Impl>);
  std::unique_ptr<Impl> impl_;
};
class PrivateFullTimeStep {
public:
  PrivateFullTimeStep(PrivateFullTimeStep &&) noexcept;
  PrivateFullTimeStep &operator=(PrivateFullTimeStep &&) noexcept;
  ~PrivateFullTimeStep();
  PrivateFullTimeStep(const PrivateFullTimeStep &) = delete;
  PrivateFullTimeStep &operator=(const PrivateFullTimeStep &) = delete;
  CollisionScope first_exposure(const ExposureBinding &);
  // Four total forks per outer reservation. Each starts the exact first frame.
  PrivateTransportTrial fork_transport_trial();
  void accept_transport_trial(PreparedPrivateMotion &&);
  CollisionScope second_exposure(const ExposureBinding &);
  PreparedArena seal(const FullTimeCommitSeals &);

private:
  friend class PicArena;
  friend class FullTimeAccess;
  struct Impl;
  explicit PrivateFullTimeStep(std::unique_ptr<Impl>);
  std::unique_ptr<Impl> impl_;
};
} // namespace raikou::mfem_backend
