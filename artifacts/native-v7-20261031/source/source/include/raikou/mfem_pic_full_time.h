#pragma once
#include "raikou/mfem_pic_arena.h"
#include <array>
#include <functional>
#include <memory>
#include <span>
#include <stdexcept>

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
// ---- native v6: births (H2), first-hit sinks (H3), binary motion v3 ----
// Typed refusal of a birth credit, sink record or binary motion record. The
// owner and every private frame are unchanged when it is thrown.
class FullTimeFlowRefusal : public std::invalid_argument {
public:
  enum class Kind {
    birth_scope,       // scope reserved twice / wrong phase / stale outer
    birth_binding,     // arrival step, base B_s, plan/token/source digests
    birth_order,       // not canonical (stream id, index) order
    birth_identity,    // foreign/duplicate packet ID, native ID != B_s + rank
    birth_namespace,   // reserved "arrival-v1/" ID from another producer
    birth_exhaustion,  // native ID space exhausted
    birth_count,       // native credit != source debit (exact), weights
    birth_kinematics,  // position/velocity/birth time/remaining flight
    sink_record,       // inconsistent with the native kinetic policy
    sink_outside_step, // hit fraction/flight/time outside the carrier flight
    sink_not_crossing, // carrier segment has no (or another) first hit
    sink_duplicate,    // carrier removed twice / also an endpoint / unknown
    sink_partial_packet, // only some carriers of one logical packet removed
    motion_record      // binary record schema/digest/census/bijection/budget
  };
  FullTimeFlowRefusal(Kind k, const std::string &message)
      : std::invalid_argument(message), kind_(k) {}
  Kind kind() const noexcept { return kind_; }

private:
  Kind kind_;
};
// Reserved logical packet namespace of the arrival owner v1. Only the birth
// credit creates IDs in it; genesis and reaction products refuse it.
inline constexpr char full_time_birth_namespace[] = "arrival-v1/";
// Canonical ID "arrival-v1/s<step>/p<stream>/k<index>" (decimal, no leading
// zeros) of an arrival birth key.
std::string full_time_birth_packet_id(std::uint64_t step,
                                      std::uint32_t stream_id,
                                      std::uint32_t index);
bool full_time_reserved_packet_id(const std::string &id) noexcept;
// One staged arrival birth (arrival owner v1 Birth after Staged::stage and
// bind_source_times) with the source owner's kinematics: neutral = source v5
// packet x and proper velocity at its release time; electron = source layer
// state at its handoff (or birth) point. Position must be a strictly interior
// vacuum point of the native kinetic box.
struct FullTimeBirth {
  std::uint64_t step = 0;
  std::uint32_t stream_id = 0, index = 0;
  std::string packet_id, state, origin;
  double weight = 0; // exactly one positive binary64 carrier
  PicVector position_m{}, proper_velocity_m_s{};
  double birth_time_s = 0, remaining_flight_s = 0; // tau, end - tau
  std::uint64_t native_id = 0;                      // staged B_s + rank
};
// One arrival stream of the step: its species state and the exact source
// debit (neutral: source v5 ReleaseQuote released = Trial inventory decrease;
// electron: n * w0 = sum of the source candidates' weights).
struct FullTimeBirthStream {
  std::uint32_t stream_id = 0;
  std::string state;
  ArenaCount source_debit;
};
struct FullTimeBirthCredit {
  std::uint64_t arrival_step = 0, base_native_id = 0;
  std::string plan_sha256, staged_token_sha256;
  std::string source_before_sha256, source_after_sha256;
  std::vector<FullTimeBirthStream> streams; // ascending stream id
  std::vector<FullTimeBirth> births;        // canonical order
};
struct FullTimeBirthEntry {
  std::string packet_id, state;
  std::uint64_t native_id = 0;
  double weight = 0;
};
// Native birth receipt (arrival NativeBirthReceipt + the journal binding).
struct FullTimeBirthReceipt {
  std::uint64_t accepted_step = 0, base = 0, next_id_after = 0;
  std::string plan_sha256, staged_token_sha256;
  std::string source_before_sha256, source_after_sha256;
  std::string births_root_sha256, row_sha256;
  std::vector<FullTimeBirthEntry> entries;  // ascending native ID
  std::vector<FullTimeBirthStream> streams; // native credit == source debit
};
// Partial flight of a newborn carrier: it moves only over [tau, end].
struct FullTimeBirthFlight {
  std::uint64_t native_id = 0;
  double birth_time_s = 0, flight_s = 0;
};
// First-hit sink of one carrier (values of hall_driver::KineticSink for the
// surface). collected_at_conductor: lower/upper foil, short bridge;
// escaped_exterior: outer x/y faces, axial (feed or non-bridge short) plane.
enum class FullTimeSinkKind : std::uint8_t {
  collected_at_conductor = 1,
  escaped_exterior = 2
};
enum class FullTimeSinkSurface : std::uint8_t {
  lower_foil = 0,
  upper_foil = 1,
  outer_xy = 2,
  short_bridge = 3,
  axial = 4
};
struct FullTimeSinkRecord {
  std::uint64_t native_id = 0;
  FullTimeSinkKind kind = FullTimeSinkKind::escaped_exterior;
  FullTimeSinkSurface surface = FullTimeSinkSurface::outer_xy;
  std::uint32_t axis = 0;  // hit plane normal axis (0 x, 1 y, 2 z)
  double plane_m = 0;      // native coordinate of the hit plane
  double fraction = 0;     // first-hit fraction of the straight segment
  PicVector segment_m{};   // full native displacement of the carrier flight
  PicVector hit_point_m{}; // start + fraction*segment, axis coordinate = plane
  PicVector incident_proper_velocity_m_s{};
  double flight_s = 0;           // carrier flight to the hit = F * fraction
  double discarded_flight_s = 0; // F - flight_s (remaining flight discarded)
};
// Exact per (surface, species) totals of removed carriers.
struct FullTimeSinkTally {
  FullTimeSinkSurface surface = FullTimeSinkSurface::outer_xy;
  FullTimeSinkKind kind = FullTimeSinkKind::escaped_exterior;
  std::string state;
  int charge_e = 0; // species charge number; charge = charge_e * e * weight
  std::uint64_t carriers = 0;
  ArenaCount weight;           // exact sum of carrier weights
  ArenaCount kinetic_energy_J; // exact sum of fl(w * KE(incident u))
};
// Exact per-surface totals (all five surfaces, ascending surface value).
struct FullTimeSinkSurfaceTotals {
  FullTimeSinkSurface surface = FullTimeSinkSurface::outer_xy;
  FullTimeSinkKind kind = FullTimeSinkKind::escaped_exterior;
  std::uint64_t carriers = 0;
  ArenaCount weight, kinetic_energy_J;
  ArenaCount positive_charge_e, negative_charge_e; // elementary charges
};
struct FullTimeSinkLedger {
  std::uint64_t accepted_step = 0, carriers = 0;
  std::string transport_trial_sha256, records_root_sha256, ledger_sha256,
      row_sha256;
  std::vector<FullTimeSinkRecord> records; // ascending native ID
  std::vector<std::string> removed_packets;
  std::vector<FullTimeSinkTally> tallies;  // (surface, state) ascending
  std::vector<FullTimeSinkSurfaceTotals> surfaces;
};
struct FullTimeStepFlows {
  std::shared_ptr<const FullTimeBirthReceipt> births;
  std::shared_ptr<const FullTimeSinkLedger> sinks;
};
struct FullTimeFlowTotals {
  std::uint64_t accepted_step = 0, next_native_id = 0;
  std::uint64_t births = 0, birth_rows = 0, sunk_carriers = 0, sink_rows = 0;
  bool has_last_birth = false;
  std::uint64_t last_birth_step = 0, last_birth_base = 0,
                last_birth_count = 0;
  std::string last_birth_plan_sha256, last_birth_token_sha256;
  std::vector<FullTimeSinkSurfaceTotals> surfaces; // cumulative
  std::string flows_sha256; // empty when the state has no flows
};
// Native port of the sealed Hall kinetic first-hit classification
// (kinetic_boundary.cpp first_kinetic_hit; surface = hall_driver::KineticSink
// value, -1 none). Same operands and operation order; it also accepts a start
// point on a foil face (a carrier leaving that face reports no hit there).
struct FullTimeKineticHit {
  int surface = -1;
  double fraction = 1;
  int axis = -1;
  double plane_m = 0;
};
FullTimeKineticHit full_time_first_kinetic_hit(const GapPicOptions &,
                                               const PicVector &start_m,
                                               const PicVector &segment_m);
// ---- native v7: device-computable motion forms (API-007.txt) ----
// Motion record v7: the v6 binary record content with block-Merkle section
// roots (1024-record leaves of 32-record chunks), the exact kinetic/work sums
// and the census v2 digest. Spans reference caller memory that must outlive
// the admission built from it (consumed by prepare_total_field_motion).
inline constexpr char motion_v7_schema[] =
    "canonical-full-time-total-field-motion-v7-merkle";
inline constexpr char motion_v7_coverage[] =
    "all_native_carriers_births_and_first_hit_sinks";
struct FullTimeMotionPath {
  std::uint64_t native_id = 0;
  PicVector field_start_m{}, field_displacement_m{};
  double flight_s = 0;
};
struct FullTimeMotionRecordV7 {
  std::string schema = motion_v7_schema;
  std::uint64_t accepted_step = 0, staged_arena_revision = 0, trial_stage = 0;
  double begin_s = 0, end_s = 0, dt_s = 0;
  std::string staged_arena_sha256, trial_field_sha256, mesh_sha256,
      material_sha256, species_sha256, common_clock_sha256,
      kinetic_policy_sha256;
  // kinetic_change (exact form), field_work (exact form), work_defect,
  // continuity, plasma_dt, maximum_gyro, maximum_displacement_cells
  std::array<double, 7> budget{};
  std::string census_sha256;            // census v2 (API-007.txt 1.3)
  std::string coverage = motion_v7_coverage;
  std::string current_A_m_sha256, endpoint_charge_C_sha256,
      initial_charge_C_sha256;
  std::span<const GapParticle> endpoints; // survivors, ascending native ID
  // width, gap, short_z, feed_z, lower l/r/b/t, upper l/r/b/t
  std::array<double, 12> kinetic_policy{};
  PicVector native_to_field_m{};
  std::string neutral_motion_owner = "same_full_time_q0_drift";
  bool native_DKD_called = false, native_Poisson_called = false;
  double neutral_stream_dt_max_s = 1e-14; // = the native configuration value
  std::span<const FullTimeMotionPath> paths;   // every census carrier
  std::span<const FullTimeSinkRecord> sinks;   // ascending native ID
  std::span<const FullTimeBirthFlight> births; // newborns of the trial
  std::string record_sha256; // producer's digest of the record
  // Admission inputs bound by the record's digests (not hashed into the record
  // digest itself): the trial midpoint total field (edge E, face B values of
  // the Hall axial mesh; trial_field_sha256 = total-midpoint-E-B-binary64-v1)
  // and the plasma current (current_A_m_sha256 = motion-array-binary64-v1).
  // Native recomputes field_work_J and the motion-work-gross-v1 tolerance
  // from them (API-007.txt 1.1 amendment).
  std::span<const double> trial_E, trial_B, current_A_m;
};
// Canonical digest of a motion record v7 (API-007.txt 1.2).
std::string full_time_motion_record_v7_sha256(const FullTimeMotionRecordV7 &);
// Census v2 of a transport trial (API-007.txt 1.3), computed once per trial.
// Spans reference the trial's storage (valid while the trial lives).
struct FullTimeCensusV7 {
  std::span<const GapParticle> carriers;       // ascending native ID
  std::span<const std::uint64_t> packet_keys;  // smallest native ID of the packet
  std::span<const std::uint64_t> state_indices; // rank in states
  std::span<const std::string> states;         // registered states, ascending
  std::string states_sha256, carriers_root, sha256;
};
// Resident-neutral interval bound (API-007.txt 1.4). Bound into the native
// configuration identity only when it differs from the default.
struct FullTimeMotionPolicy {
  double neutral_stream_dt_max_s = 1e-14;
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
  FullTimeMotionPolicy motion; // v7
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
  // v6: birth/sink rows and carriers in the window; cumulative flows digest
  // at the anchor (empty when the state has no births or sinks).
  std::uint64_t window_birth_rows = 0, window_births = 0, window_sink_rows = 0,
                window_sunk_carriers = 0;
  std::string flows_sha256;
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
// Live motion admission is created only from a typed motion record v7. JSON
// cannot authorize live motion.
class FullTimeMotionAdmission {
public:
  ~FullTimeMotionAdmission() noexcept;
  FullTimeMotionAdmission(FullTimeMotionAdmission &&) noexcept = default;
  FullTimeMotionAdmission &operator=(FullTimeMotionAdmission &&) noexcept;
  FullTimeMotionAdmission(const FullTimeMotionAdmission &) = delete;
  FullTimeMotionAdmission &operator=(const FullTimeMotionAdmission &) = delete;

private:
  friend FullTimeMotionAdmission
  bind_total_field_motion_v7(const FullTimeMotionRecordV7 &,
                             std::shared_ptr<const void>);
  friend class PrivateTransportTrial;
  FullTimeMotionAdmission() = default;
  std::shared_ptr<const FullTimeMotionRecordV7> v7_;
  // Keeps producer-owned arrays alive until the admission is consumed.
  std::shared_ptr<const void> v7_storage_;
};
// Typed motion v7 admission. The record digest must equal
// record.record_sha256 and the motion binding's transport_trial_sha256.
FullTimeMotionAdmission
bind_total_field_motion_v7(const FullTimeMotionRecordV7 &,
                           std::shared_ptr<const void> storage = {});

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
  // v7: census v2 of this trial (API-007.txt 1.3), computed at most once per
  // trial; census_sha256() = census_v7().sha256.
  const FullTimeCensusV7 &census_v7() const;
  const std::string &census_sha256() const;
  // v7: packet ID and species state of a census carrier (producer bookkeeping
  // of sink and newborn records).
  std::pair<std::string, std::string> carrier_packet(std::uint64_t native_id) const;
  // v6: newborn carriers of this trial (ascending native ID) and the source
  // seal after the births (accepted source seal when the trial has none).
  std::vector<FullTimeBirthFlight> birth_flights() const;
  std::string source_after_sha256() const;
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
// Prepared births: the private post-first-exposure frame plus one birth
// credit. Pure with respect to the outer step; a re-evaluated source step
// prepares another one from the same base. Consumed only by forking a
// transport trial from it (any number of forks within the outer fork cap).
class PreparedBirths {
public:
  PreparedBirths(PreparedBirths &&) noexcept;
  PreparedBirths &operator=(PreparedBirths &&) noexcept;
  ~PreparedBirths();
  PreparedBirths(const PreparedBirths &) = delete;
  PreparedBirths &operator=(const PreparedBirths &) = delete;
  const FullTimeBirthReceipt &receipt() const;
  const std::vector<FullTimeBirthFlight> &flights() const;

private:
  friend class BirthScope;
  friend class PrivateFullTimeStep;
  friend class FullTimeAccess;
  struct Impl;
  explicit PreparedBirths(std::unique_ptr<Impl>);
  std::unique_ptr<Impl> impl_;
};
// Birth scope of one outer step (reserved once, after the first exposure and
// before transport). base_native_id() is B_s: native's next ID at the birth
// insertion point, fixed for every re-evaluation of the step.
class BirthScope {
public:
  BirthScope(BirthScope &&) noexcept;
  BirthScope &operator=(BirthScope &&) noexcept;
  ~BirthScope();
  BirthScope(const BirthScope &) = delete;
  BirthScope &operator=(const BirthScope &) = delete;
  std::uint64_t base_native_id() const;
  std::uint64_t accepted_step() const;
  const std::string &source_before_sha256() const;
  PreparedBirths prepare(const FullTimeBirthCredit &) const;

private:
  friend class PrivateFullTimeStep;
  friend class FullTimeAccess;
  struct Impl;
  explicit BirthScope(std::unique_ptr<Impl>);
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
  // v6: birth scope (once per outer) and a fork that starts from the
  // post-birth frame of a prepared birth credit (counts toward the same cap).
  BirthScope birth_scope();
  PrivateTransportTrial fork_transport_trial(const PreparedBirths &);
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
