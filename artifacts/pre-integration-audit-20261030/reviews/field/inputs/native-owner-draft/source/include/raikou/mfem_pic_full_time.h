#pragma once
#include "raikou/mfem_pic_arena.h"
#include <functional>
#include <memory>
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
  std::string outer_step_sha256;
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
  friend class PrivateFullTimeStep;
  struct Impl;
  explicit PreparedPrivateMotion(std::unique_ptr<Impl>);
  std::unique_ptr<Impl> impl_;
};
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
  PreparedPrivateMotion prepare_total_field_motion(
      const FullTimeMotionBinding &, const std::vector<GapParticle> &endpoints,
      const ArenaJson &transport_receipts,
      const ArenaJson &source_receipts = ArenaJson::array(),
      const ArenaJson &boundary_receipts = ArenaJson::array());

private:
  friend class FullTimeAccess;
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
