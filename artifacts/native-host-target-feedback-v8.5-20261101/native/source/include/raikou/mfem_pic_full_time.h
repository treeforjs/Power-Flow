#pragma once
#include "raikou/mfem_pic_arena.h"
#include <array>
#include <functional>
#include <memory>
#include <span>
#include <stdexcept>

struct MsParams; // merge-split-transaction-v1 policy (resolution-series-prep merge-split ms_core.h)
namespace raikou::mfem_backend {
class FullTimeHostTargetProvider;
struct FullTimeHostTargetReceipt;
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
    motion_record,     // binary record schema/digest/census/bijection/budget
    // v8.1 merge/split (merge-split INTEGRATION 2.4)
    merge_scope,       // not a transaction step, reserved twice, wrong phase, foreign/stale transaction
    merge_binding,     // producer's transaction digest differs from native's recomputation
    merge_kinematics,  // an output not strictly interior (v5 rule)
    merge_count,       // exact count equality of a species fails
    merge_exhaustion   // native ID space exhausted by the outputs
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
  // v8.1 (identity schema v8 only; API-008.txt 8.1-8.3):
  // matching-plane emission birth (source owner v7 B3): born at bridge_start_m +
  // L n (n the outward normal of the emitting foil face, L > 0), with a
  // zero-duration bridge path from the face point carrying charge q w (Hall
  // deposits it; its field work is the motion record's extra work).
  bool bridge = false;
  PicVector bridge_start_m{};
  // born-and-collected birth (source owner v7 B4, ReturnedSubcell): credited and
  // collected at its emission point (on an emitting foil face) in this step; it
  // takes a native ID but never becomes a carrier.
  bool returned = false;
  // internal energy carried by the macro (SS304 phase owner v1.1 3.2: released
  // internal energy / n per (patch, element)); neutral states only.
  double internal_energy_J = 0;
};
// v8.1 (source owner v7 INTEGRATION B2): a birth authored from a wall
// response. Its parent is a first-hit sink published by the previous accepted
// step (FullTimeSinkExposure); weight and position are exactly the parent's
// weight and hit point; it flies the whole step (birth time = step begin; the
// one-step lag is recorded by parent_hit_time_s). Identity: (wall response
// content digest, parent native ID, ordinal).
inline constexpr char full_time_wall_birth_namespace[] = "surface-v7/";
// Canonical ID "surface-v7/s<step>/c<parent native ID>/k<ordinal>".
std::string full_time_wall_birth_packet_id(std::uint64_t step, std::uint64_t parent_native_id,
                                           std::uint32_t ordinal);
struct FullTimeWallBirth {
  std::string wall_response_sha256;
  std::uint64_t parent_native_id = 0;
  std::uint32_t ordinal = 0; // 0..k-1 per parent (k <= 64)
  std::string packet_id, state, origin;
  double weight = 0;              // = the parent's weight (bitwise)
  PicVector position_m{};         // = the parent's hit point (native frame, bitwise)
  PicVector proper_velocity_m_s{};
  double parent_hit_time_s = 0;   // parent's hit time (<= step begin)
  std::uint64_t native_id = 0;    // B_s + arrival births + rank
};
// Zero-duration bridge path of a matching-plane birth (field frame, as motion
// paths): start = face point + shift, displacement = birth point - face point.
struct FullTimeBirthBridge {
  std::uint64_t native_id = 0;
  PicVector field_start_m{}, field_displacement_m{};
  double charge_C = 0; // q w
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
  // v8.1: wall-response births, ascending (parent native ID, ordinal); native
  // IDs follow the arrival births.
  std::vector<FullTimeWallBirth> wall;
};
struct FullTimeBirthEntry {
  std::string packet_id, state;
  std::uint64_t native_id = 0;
  double weight = 0;
  bool returned = false; // v8.1: born and collected (no carrier)
};
// Native birth receipt (arrival NativeBirthReceipt + the journal binding).
struct FullTimeBirthReceipt {
  std::uint64_t accepted_step = 0, base = 0, next_id_after = 0;
  std::string plan_sha256, staged_token_sha256;
  std::string source_before_sha256, source_after_sha256;
  std::string births_root_sha256, row_sha256;
  std::vector<FullTimeBirthEntry> entries;  // ascending native ID
  std::vector<FullTimeBirthStream> streams; // native credit == source debit
  // v8.1: next ID after the arrival births (= next_id_after without wall
  // births; the arrival owner's receipt uses it), wall-response births
  // (ascending native ID) and their exact credit per state (stream_id 0),
  // bridge paths, and the exact internal energy carried by the births (J).
  std::uint64_t arrival_next_id_after = 0;
  std::vector<FullTimeBirthEntry> wall_entries;
  std::vector<FullTimeBirthStream> wall_credits;
  std::vector<FullTimeBirthBridge> bridges;
  std::string internal_energy_J_digits; // exact sum (forms digits); empty when zero
  double internal_energy_J = 0;
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
// v8 (source owner v7 INTEGRATION B1): one published first-hit sink with the
// fields the next step's wall response needs. Audit schema v8.2 binds the
// complete exposure into the accepted snapshot, including at history anchors.
struct FullTimeSinkExposure {
  FullTimeSinkRecord record;
  std::string packet_id, state;
  double weight = 0, charge_C = 0, mass_kg = 0;
  PicVector start_position_m{}; // carrier position at the step start
  // v8.1: internal energy carried by the sunk macro (metal vapour births; the
  // coordinator returns it as wall heat with the deposit, phase owner 3.2).
  double internal_energy_J = 0;
};
struct FullTimeSinkLedger {
  std::uint64_t accepted_step = 0, carriers = 0;
  std::string transport_trial_sha256, records_root_sha256, ledger_sha256,
      row_sha256;
  std::vector<FullTimeSinkRecord> records; // ascending native ID
  std::vector<std::string> removed_packets;
  std::vector<FullTimeSinkTally> tallies;  // (surface, state) ascending
  std::vector<FullTimeSinkSurfaceTotals> surfaces;
  std::vector<FullTimeSinkExposure> exposed; // v8: ascending native ID (as records)
};
struct FullTimeMergeSplitReceipt;
struct FullTimeStepFlows {
  std::shared_ptr<const FullTimeBirthReceipt> births;
  std::shared_ptr<const FullTimeSinkLedger> sinks;
  std::shared_ptr<const FullTimeMergeSplitReceipt> merge; // v8.1 (transaction steps)
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
// v8.1 motion record with extras (API-008.txt 8.4): the v7 record plus the
// zero-time current sources that share the plasma current but are not motion
// work: merge/split relocations (merge-split INTEGRATION 2.9) and bridge paths
// of matching-plane births (source owner v7 B3/C4). record_sha256 =
// SHA256(text "native-v8.1/motion-extras", text v7 digest, u64 relocations,
// text relocations root, u64 bridges, text bridges root, text
// extra_current_A_m_sha256, f64 extra_work_J).
inline constexpr char motion_v81_schema[] = "canonical-full-time-total-field-motion-v8.1-merkle-extras";
// One relocation record (72 bytes, field frame): input ID, group rank, start,
// displacement D = (c + shift) - (x + shift), charge w q.
struct FullTimeMotionRelocation {
  std::uint64_t native_id = 0, group = 0;
  PicVector field_start_m{}, field_displacement_m{};
  double charge_C = 0;
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
  // v8.1 extras (schema motion_v81_schema only; empty in v7 records):
  // relocations (transaction order), bridges (ascending native ID), the extra
  // current J_extra of both on the Hall edges (included in current_A_m with one
  // rounding per edge; digest motion-array-binary64-v1 or the block form) and
  // extra_work_J = RN(dt x RN(exact sum E J_extra)). Admission: work defect =
  // kinetic_change - (field_work - extra_work), motion-work-gross-v1 tolerance.
  std::span<const FullTimeMotionRelocation> relocations;
  std::span<const struct FullTimeBirthBridge> bridges;
  std::span<const double> extra_current_A_m;
  std::string extra_current_A_m_sha256;
  double extra_work_J = 0;
  // v7.1 device lane: the same arrays in device memory on the shared session
  // (addresses; 0 when absent). Canonical layouts of API-007.txt 1.2 (sinks as
  // 136-byte records); E/B/J binary64. When given, the device lane reads them
  // in place; the host spans are then needed only by the host path.
  struct DeviceArrays {
    std::uint64_t endpoints = 0, paths = 0, sinks = 0, births = 0, E = 0, B = 0, J = 0;
    std::uint64_t extra_J = 0; // v8.1: extra_current_A_m on the device (n_E values)
    std::uint64_t n_endpoints = 0, n_paths = 0, n_sinks = 0, n_births = 0, n_E = 0, n_B = 0;
    bool present() const noexcept { return endpoints || paths || E || J; }
  } device;
};
// Canonical digest of a motion record v7 (API-007.txt 1.2).
std::string full_time_motion_record_v7_sha256(const FullTimeMotionRecordV7 &);
// v8.1: digest of a record with the extras schema (motion_v81_schema).
std::string full_time_motion_record_v81_sha256(const FullTimeMotionRecordV7 &);
// Census v2 of a transport trial (API-007.txt 1.3), computed once per trial.
// Spans reference the trial's storage (valid while the trial lives).
struct FullTimeCensusV7 {
  std::span<const GapParticle> carriers;       // ascending native ID
  std::span<const std::uint64_t> packet_keys;  // smallest native ID of the packet
  std::span<const std::uint64_t> state_indices; // rank in states
  std::span<const std::string> states;         // registered states, ascending
  std::string states_sha256, carriers_root, sha256;
};
// v7.1: string-free view of the accepted population (API-071.txt 3): the
// native rows (ascending native ID), the packet key (smallest native ID of the
// logical packet) and the state index per row, the state table, and per-state
// totals (carriers; exact sums of weights and of weight x charge, rounded
// once). Recomputed only when the packet table changes.
struct FullTimePopulationView {
  std::span<const GapParticle> rows;
  std::span<const std::uint64_t> packet_keys, state_indices;
  std::span<const std::string> states;
  struct StateTotals {
    std::uint64_t carriers = 0;
    double weight = 0, charge_C = 0;
  };
  std::vector<StateTotals> totals; // per state index
  std::shared_ptr<const void> storage; // keeps the arrays alive
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
// ---- native v8 (API-008.txt) ----
// Accepted-state identity schema. v5: the v7.1 identities (streamed packet
// table and ID sets; compatible mode). v8: content-chunked Merkle packet table
// and ID-set digests (API-008.txt 2), required by chemistry v3.
enum class FullTimeIdentitySchema { v5 = 0, v8 = 1 };
// Where a chemistry v3 stage runs: host (CPU reference), device (attached
// device lane; host otherwise), device_verified (both; any difference in
// records, ledgers or the applied frame is an internal error).
enum class FullTimeChemistryExecution { host = 0, device = 1, device_verified = 2 };
// coverage-ledger-v3 refusal policy (schema pfc-chem-v3-coverage-policy-1):
// refuse a step iff a bounded region has B_r > max(event_limit,
// event_fraction * A) (A = admitted events of the step) or an unbounded region
// has reachable carrier support > 0 and is not admitted (including zero-speed
// encounters with physical X_r = B_r = 0). Region names are those of the blob.
// A = the applied physical events of the step (exact sum of the admitted
// events' weights W over both stages), as chemistry v3.1's policy file states.
// policy_file (schema pfc-chem-v3.1-coverage-policy-1, sha256 policy_file_sha256)
// supplies event_limit, event_fraction and the admitted regions instead of the
// three fields (which must then be left empty/default); its blob sha256, region
// count, names and bounded set must match the loaded blob.
struct FullTimeChemistryPolicy {
  double event_limit = 1.0, event_fraction = 0.0;
  std::vector<std::string> admitted_regions;
  std::string policy_file, policy_file_sha256;
};
// Chemistry v3 (pic-water-chemistry-v3-20261031, INTEGRATION.txt): two
// collision stages per step (dt/2 each, purposes 0x30 and 0x34) on the native
// census, products, target debits, the quantized neutral release and exact
// ledgers. Requires identity schema v8.
enum class FullTimeWaterScopeMode : std::uint32_t { legacy = 0, minimum_water = 1 };
struct FullTimeChemistryConfiguration {
  bool enabled = false;
  // Diagnostic field-target MCC has no owner for target energy/momentum.
  // Default refuses events. Opt-in is bound into configuration/checkpoint identity;
  // it never establishes whole-system energy closure.
  bool diagnostic_field_target_approximation = false;
  // Optional minimum encounter gate; passing never establishes network completeness.
  bool host_target_feedback = false;
  std::string host_target_contract_sha256;
  FullTimeWaterScopeMode water_scope_mode = FullTimeWaterScopeMode::legacy;
  std::string water_scope_contract, water_scope_contract_sha256;
  std::string blob_path, blob_sha256;  // table blob and its expected identity
  std::uint64_t run_seed = 0;          // the arrival owner's run seed
  double release_weight = 0;           // w0 of the quantized release
  double trunc_sigmas = 7;
  std::uint32_t max_candidates = 4096, max_events = 64;
  std::array<std::uint32_t, 3> grid_counts{}; // collision grid (x, y, z cells)
  FullTimeChemistryPolicy policy;
  std::string preset; // named preset ("production-v1"); bound into the identity when set
  // Not bound into identities (outputs are identical for every value):
  std::uint32_t host_partitions = 1;   // carrier partitions of the host stage
  FullTimeChemistryExecution execution = FullTimeChemistryExecution::host;
  bool retain_stage_inputs = false;    // keep CarrierIn/CellTargets bytes in the receipt
  std::uint32_t device_partitions = 0; // device stage launches (0: by the 512 MiB event budget)
};
// Typed chemistry refusal ("chemistry_v3/<kind>"); nothing changes when thrown.
class FullTimeChemistryRefusal : public std::invalid_argument {
public:
  enum class Kind {
    status,            // a carrier record status other than OK/CONSUMED
    cell_targets,      // more than 12 target species in one collision cell
    coverage_policy,   // coverage-ledger-v3 refusal policy
    projectile_packet, // projectile carrier of a multi-carrier packet
    carrier_domain,    // native ID >= 2^48, carrier outside the grid
    debit_conflict,    // admitted debits exceed the carriers left by the stage
    configuration,     // blob/registry/configuration mismatch
    target_energy_owner, // event requires a target-energy owner; diagnostic opt-in only
    target_fields,     // invalid/unqualified external target input or stale owner
    network_scope,     // absent/unsupported required encounter or unqualified census lane
    event              // event record fails conservation/destination checks
  };
  FullTimeChemistryRefusal(Kind k, const std::string &message)
      : std::invalid_argument(message), kind_(k) {}
  Kind kind() const noexcept { return kind_; }

private:
  Kind kind_;
};
// Registered chemistry species (state "pfc-v3/<name>"), from the blob.
struct FullTimeChemistrySpecies {
  std::string state;
  double mass_kg = 0, mass_eV = 0, reference_eV = 0;
  int charge_e = 0;
  std::array<std::int32_t, 16> nuclei{}; // nuclei ABI v4 (0 H, 1 O, 2 Fe, 3 Cr, 4 Ni, ...)
  bool projectile = false, target = false;
};
struct FullTimeChemistryRegistry {
  std::string blob_sha256, registry_sha256, configuration_sha256, policy_sha256, grid_sha256;
  std::string initial_rng_sha256; // chemistry state digest of a fresh run
  std::vector<FullTimeChemistrySpecies> species; // blob index order
  ArenaChannel channel;                          // channel to register
  std::vector<std::string> region_names;
  std::vector<bool> region_bounded;
};
// Loads and validates the blob; registry, identities and the channel binding.
FullTimeChemistryRegistry full_time_chemistry_registry(const FullTimeChemistryConfiguration &,
                                                       const GapPicOptions &);
// Build the canonical minimum contract against the base physics configuration.
// Input must have legacy mode and empty scope fields; callers explicitly opt in.
std::string full_time_minimum_water_scope_contract(const FullTimeChemistryConfiguration &, const GapPicOptions &);
// Arena state specs of the 22 species (hydrogen/oxygen counts, charge, mass).
std::vector<ArenaStateSpec> full_time_chemistry_states(const FullTimeChemistryRegistry &);
// Recommended admitted unbounded regions (SCOPE.md section 5 placeholders and
// every region of the 1e7 eV decade).
std::vector<std::string> full_time_chemistry_recommended_admitted(const FullTimeChemistryRegistry &);
// Named configuration presets (API-008.txt 3.1). "production-v1": chemistry
// v3.1 (blob and coverage-policy file bound by sha256 below; refused when they
// differ), max_events 64, max_candidates 4096, the policy file's admitted list,
// bounded regions refused when B_r > max(1, 0.01 A). An empty blob_path selects
// the default v3.1 location.
inline constexpr char full_time_chemistry_production_v1_blob_sha256[] =
    "70fcc2adc2fdaa62abed7ed02d68ffc09b1b18554eda8fa99c5451120ae454b8";
inline constexpr char full_time_chemistry_production_v1_policy_sha256[] =
    "15d44a91b8af93e0b14beb89434c55dc402e7196f10ca4bdaa2357564474b2d3";
inline constexpr char full_time_chemistry_v31_tables[] =
    "/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031/data/tables/";
FullTimeChemistryConfiguration full_time_chemistry_preset(const std::string &name, const GapPicOptions &,
                                                          const std::string &blob_path, std::uint64_t run_seed,
                                                          double release_weight,
                                                          const std::array<std::uint32_t, 3> &grid_counts);
// Exact value of a ledger quantity: correctly rounded binary64 and the
// normalised forms_core digits (hexadecimal text of the 134 digits).
struct FullTimeExactValue {
  double value = 0;
  std::string digits;
};
// Native v8.3 diagnostic export of admitted events and actual packet debits.
// These immutable facts are not publication authority: a consumer must also
// preflight the live PreparedArena candidate in a common transaction.
struct FullTimeChemistryAdmission {
  std::uint64_t carrier = 0, raw_begin = 0, raw_count = 0, admitted_count = 0;
  std::uint32_t compact_cell = 0, atlas_cell = 0;
};
struct FullTimeChemistryAdmittedEvent {
  std::uint64_t raw_index = 0;
  std::uint32_t atlas_cell = 0, target_species = 0;
  bool consumes_tracked_target = false;
  double weight = 0;
  std::string event_bytes; // exact pinned C3 EventRec ABI, not a native struct dump
};
struct FullTimeChemistryPacketSnapshot {
  std::string packet_id, state;
  ArenaCount count;
  std::vector<GapParticle> carriers; // ascending native ID; exact positive weights
  double carried_internal_energy_J = 0;
};
struct FullTimeChemistryTargetDebit {
  std::uint32_t atlas_cell = 0, species = 0;
  ArenaCount amount;
  FullTimeChemistryPacketSnapshot before, remainder;
  bool also_projectile = false;
  bool internal_energy_deferred = false;
};
struct FullTimeChemistryStageExport {
  std::string schema = "native-v8.5/admitted-target-export-v2";
  std::uint64_t accepted_step = 0, stage = 0, purpose = 0;
  std::uint64_t base_native_id = 0, next_native_id = 0;
  std::uint64_t owner_stage_sequence = 0;
  double begin_s = 0, end_s = 0, step_duration_s = 0, collision_duration_s = 0;
  double frozen_physical_stage_s = 0;
  std::string accepted_native_sha256, common_clock_sha256, frame_snapshot_sha256,
      frame_history_sha256, configuration_sha256, registry_sha256,
      inputs_sha256, chemistry_before_sha256, chemistry_after_sha256,
      frozen_field_sha256, ledger_sha256, births_root, raw_records_root, host_target_feedback_sha256;
  std::vector<FullTimeChemistryAdmission> admissions;
  std::vector<FullTimeChemistryAdmittedEvent> events;
  std::vector<FullTimeChemistryTargetDebit> target_debits;
  std::vector<FullTimeChemistryPacketSnapshot> releases;
  std::string admitted_root, target_debits_root, releases_root, sha256;
};
// Checks canonical commitments and exact inventory identities. The default also
// refuses target/projectile overlap or deferred internal-energy attribution.
// Passing false permits diagnostic inspection only, never gas publication.
void validate_full_time_chemistry_stage_export(const FullTimeChemistryStageExport &,
                                              bool require_unambiguous_targets = true);
// One chemistry stage (API-008.txt 3): counts, digests, exact ledgers.
struct FullTimeChemistryStageReceipt {
  bool diagnostic_field_target_approximation = false;
  bool full_system_energy_closure = false; // never asserted by this component
  std::uint64_t accepted_step = 0, stage = 0, purpose = 0;
  double dt_s = 0, nu_dt_max = 0;
  std::uint64_t projectiles = 0, cells = 0, candidates = 0, events = 0, admitted_events = 0,
                refused_debits = 0, births = 0, consumed = 0, debited_packets = 0,
                removed_packets = 0, released_packets = 0, remainder_carriers = 0;
  std::uint64_t base_native_id = 0, next_native_id = 0;
  std::string inputs_sha256, records_sha256, records_root, ledger_sha256,
      coverage_sha256, state_before_sha256, state_after_sha256, births_root;
  // Empty in legacy mode; native-created pre-filter census commitment otherwise.
  std::string water_scope_contract_sha256, water_scope_census_sha256, water_scope_binding_sha256;
  std::shared_ptr<const FullTimeHostTargetReceipt> host_target_feedback;
  // exact energy ledger (eV): ke_in_projectile, ke_in_target, ke_out_survivor,
  // ke_out_birth, ke_out_reservoir, ke_out_gas, q_reference, internal[6], defect
  std::vector<std::pair<std::string, FullTimeExactValue>> energy;
  std::vector<std::pair<std::string, FullTimeExactValue>> counts; // created/destroyed per species, charge, nuclei
  std::uint64_t coverage_regions_exposed = 0;
  bool policy_refused = false;
  // Canonical stage inputs (retain_stage_inputs): StepConfig, CarrierIn and
  // CellTargets bytes in the C3 layouts, for replay.
  std::string step_config_bytes, carrier_bytes, cell_bytes;
  // Device accounting (device modes): device wall time (ms) of the stage.
  double device_ms = 0;
  bool device = false;
  std::shared_ptr<const FullTimeChemistryStageExport> admitted_export;
};
// Ejected melt of one phase-owner step (StepRecord::eject_N[6] Fe, Cr, Ni, Mn,
// Si, C atoms; eject_E total and eject_K kinetic energy, J) and the producer's
// digest of the step record it comes from.
struct FullTimeEjectedMelt {
  std::array<double, 6> atoms{}; // Fe, Cr, Ni, Mn, Si, C (nuclei ABI v4 slots 2..7)
  double energy_J = 0, kinetic_J = 0;
  std::string phase_record_sha256;
};
struct FullTimeMeltReservoir {
  std::uint64_t credits = 0;
  std::array<FullTimeExactValue, 6> atoms; // cumulative, exact
  FullTimeExactValue energy_J, kinetic_J;
  std::string last_record_sha256, sha256;
};
// Accepted chemistry state summary (API-008.txt 3.6).
struct FullTimeChemistryStatus {
  std::string water_scope_census_sha256, water_scope_binding_sha256, water_scope_history_sha256;
  std::string host_target_feedback_sha256, host_target_history_sha256;
  bool enabled = false;
  bool diagnostic_field_target_approximation = false;
  bool full_system_energy_closure = false;
  std::string state_sha256, configuration_sha256, registry_sha256, blob_sha256;
  std::uint64_t stages = 0, reservoirs = 0, events = 0;
  std::vector<std::pair<std::string, FullTimeExactValue>> cumulative_energy, cumulative_counts;
  // exact reservoir stock per species state (sum over cells)
  std::vector<std::pair<std::string, FullTimeExactValue>> reservoir_counts;
  std::string cumulative_ledger_sha256, coverage_json;
  FullTimeMeltReservoir melt; // v8: ejected-melt reservoir input (credits = 0: none)
};
// ---- native v8.1: merge/split transaction (merge-split INTEGRATION 2) ----
// Policy of the merge-split-transaction-v1 component (ms_core.h MsParams:
// species, energy bins, direction and spatial sub-bins, merge/split ranges,
// Hall mesh, native -> field shift; seed = the arrival run seed; native sets
// step and base_id per transaction). species_states[k] is the native state of
// policy species k (mass and charge must equal the registered state). A
// transaction step is one with accepted_step % cadence_steps == 0. Bound into
// the configuration identity (policy digest, cadence, states); schema v8 only.
struct FullTimeMergeSplitConfiguration {
  bool enabled = false;
  std::uint64_t cadence_steps = 1;
  std::shared_ptr<const ::MsParams> params;
  std::vector<std::string> species_states;
};
struct FullTimeMergeSplitReceipt {
  std::uint64_t accepted_step = 0, base_native_id = 0, outputs = 0, next_native_id = 0;
  std::uint64_t segments = 0, groups = 0, merge_groups = 0, split_groups = 0, accepted_merge = 0,
                accepted_split = 0, inputs_consumed = 0, relocations = 0, formation_refused = 0,
                centroid_fallback = 0, isotropic = 0;
  std::array<std::uint64_t, 4> status{}; // MS_OK, no real solution, operand range, RNG exhausted
  std::string params_sha256, policy_sha256, census_sha256, transaction_sha256, records_root, relocation_root,
      ledger_sha256, row_sha256;
  std::vector<FullTimeMotionRelocation> relocation_records; // transaction order (Hall J_reloc)
  std::vector<std::uint64_t> removed;                       // input native IDs, ascending
  // per state: exact residuals (value; digits = signed hexadecimal integer of
  // units 2^-1074): momentum x, y, z (kg m/s), kinetic energy (J); in + residual = out
  std::vector<std::pair<std::string, std::array<FullTimeExactValue, 4>>> residuals;
};
struct FullTimeConfiguration {
  FullTimeNativeLimits limits;
  FullTimeHistoryPolicy history;
  FullTimeMotionPolicy motion; // v7
  FullTimeIdentitySchema identity = FullTimeIdentitySchema::v5; // v8
  FullTimeChemistryConfiguration chemistry;                    // v8
  FullTimeMergeSplitConfiguration merge_split;                 // v8.1
};
// v8.1 accepted flows summary (API-008.txt 8.5).
struct FullTimeFlowsV81Status {
  bool present = false;
  std::uint64_t wall_births = 0, returned_births = 0, bridges = 0, carrying_packets = 0;
  FullTimeExactValue internal_born_J, internal_sunk_J, internal_released_J, internal_in_flight_J, extra_work_J;
  std::vector<std::pair<std::string, FullTimeExactValue>> wall_weight, returned_weight;
  std::uint64_t merge_transactions = 0, merge_groups = 0, merge_inputs = 0, merge_outputs = 0, merge_relocations = 0,
                merge_formation_refused = 0, merge_centroid_fallback = 0;
  std::array<std::uint64_t, 4> merge_status{};
  std::vector<std::pair<std::string, std::array<FullTimeExactValue, 4>>> merge_residual;
  bool has_last_merge = false;
  std::uint64_t last_merge_step = 0, last_merge_base = 0, last_merge_outputs = 0;
  std::string last_merge_sha256, merge_policy_sha256, sha256;
};
// v8.1 metal vapour states (SS304 phase owner v1.1 3.1): neutral Mn, Si, C
// ("metal-v1/Mn", "metal-v1/Si", "metal-v1/C"; masses 54.938043, 28.085,
// 12.011 u) and, without chemistry v3 (whose registry has pfc-v3/Fe, Cr, Ni),
// also Fe, Cr, Ni (55.845, 51.9961, 58.6934 u). Mass = u x 1.66053906660e-27 kg.
std::vector<ArenaStateSpec> full_time_metal_vapour_states(bool with_fe_cr_ni);
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
  // v8: runs this exposure's chemistry v3 stage (stage 1 before transport,
  // stage 2 with the quantized release after it) on the private frame and
  // returns the completed exposure token (publish it with try_publish).
  // Throws FullTimeChemistryRefusal (typed) or FullTimeResourceRefusal.
  PreparedArena prepare_chemistry_v3();
  PreparedArena prepare_chemistry_v3_host_targets(std::shared_ptr<FullTimeHostTargetProvider>);
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
  // v7.1: census v2 on the device (requires an attached device; computed at
  // most once per trial; no host row access). Declared in
  // mfem_pic_full_time_device.h.
  const struct FullTimeCensusDevice &census_v7_device() const;
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
// v8.1 prepared merge/split transaction: the private post-transaction frame
// (survivors + outputs, IDs M_s + rank) and its receipt. Consumed by
// birth_scope(const PreparedMergeSplit &) or fork_transport_trial(const
// PreparedMergeSplit &); enters the outer only with the accepted motion.
class PreparedMergeSplit {
public:
  PreparedMergeSplit(PreparedMergeSplit &&) noexcept;
  PreparedMergeSplit &operator=(PreparedMergeSplit &&) noexcept;
  ~PreparedMergeSplit();
  PreparedMergeSplit(const PreparedMergeSplit &) = delete;
  PreparedMergeSplit &operator=(const PreparedMergeSplit &) = delete;
  const FullTimeMergeSplitReceipt &receipt() const;

private:
  friend class MergeScope;
  friend class PrivateFullTimeStep;
  friend class FullTimeAccess;
  struct Impl;
  explicit PreparedMergeSplit(std::unique_ptr<Impl>);
  std::unique_ptr<Impl> impl_;
};
// v8.1 merge scope (once per transaction step, after the first exposure and
// before the birth scope). base_native_id() = M_s (native next ID after the
// first-exposure products). prepare() runs the merge-split-transaction-v1 CPU
// reference on the frame census (deterministic, rank independent); a producer
// that computed the transaction elsewhere (device) passes its digest, which
// must equal native's (merge_binding refusal otherwise).
class MergeScope {
public:
  MergeScope(MergeScope &&) noexcept;
  MergeScope &operator=(MergeScope &&) noexcept;
  ~MergeScope();
  MergeScope(const MergeScope &) = delete;
  MergeScope &operator=(const MergeScope &) = delete;
  std::uint64_t base_native_id() const;
  std::uint64_t accepted_step() const;
  PreparedMergeSplit prepare(const std::string &expected_transaction_sha256 = {}) const;

private:
  friend class PrivateFullTimeStep;
  friend class FullTimeAccess;
  struct Impl;
  explicit MergeScope(std::unique_ptr<Impl>);
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
  // v8.1: merge scope (transaction steps), birth scope on a prepared
  // transaction (B_s = M_s + n_out) and a trial from the post-transaction frame.
  MergeScope merge_scope();
  BirthScope birth_scope(const PreparedMergeSplit &);
  PrivateTransportTrial fork_transport_trial(const PreparedMergeSplit &);
  void accept_transport_trial(PreparedPrivateMotion &&);
  CollisionScope second_exposure(const ExposureBinding &);
  // v8 (SS304 phase owner v1.2 INTEGRATION 2.3): the phase owner's ejected melt of
  // this step as a reservoir input (no droplet species yet). At most once per
  // outer step, before seal; schema v8 only. Exact cumulative atoms per element
  // and energies enter the accepted state (identity, checkpoint).
  void credit_ejected_melt(const struct FullTimeEjectedMelt &);
  PreparedArena seal(const FullTimeCommitSeals &);

private:
  friend class PicArena;
  friend class FullTimeAccess;
  struct Impl;
  explicit PrivateFullTimeStep(std::unique_ptr<Impl>);
  std::unique_ptr<Impl> impl_;
};
} // namespace raikou::mfem_backend
