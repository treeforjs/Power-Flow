// gas_owner.hpp -- gas target owner v1: host API.
//
// The owner holds, per (native collision-grid cell, gas species), the exact
// energy, momentum and internal-energy deposit of the neutral gas on top of
// the moments of native's carriers of that species (README 2), every exact
// ledger of the target side of chemistry v3 events (README 3), the population
// flow ledgers (README 5) and the deposit transport (README 5.2). It feeds the
// local gas state back to chemistry as the target Maxwellian of each
// CellTargets slot (README 4). All mutations go into a step candidate that is
// committed (bound to its digest) or abandoned with the coordinator's step.
//
// Execution: host (CPU reference) or device (GPU-resident state, NVRTC
// kernels compiled from gas_core.h and forms_core.h) or device_verified (both,
// every digest compared). CPU/device equality requires qualification evidence.
#pragma once
#include "pf_gas_v1/gas_core.h"

#include <array>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace pf_gas_v1 {

using pfg::CellTargetsAbi;
using pfg::FlowRecord;
using pfg::MoverRecord;
using pfg::RecoilRecord;

enum class Transport : std::uint32_t { none = 0, carrier_shares = 1 };
enum class Execution { host, device, device_verified };

struct Species {
  std::string state;        // native state name ("pfc-v3/H2O", "metal-v1/Mn", ...)
  double mass_kg = 0;       // native registry mass (bitwise)
  int charge = 0;           // elementary charges
  bool chemistry_target = false;
};
struct Process {
  std::uint32_t projectile = 0, n_out = 0, internal_carrier = pfg::kNone, kind = 0;
};
// native-v8/collision-grid-v1 geometry (API-008 3.3): volume as native
// computes it ((hx / nx) * (hy / ny)) * (hz / nz).
struct Grid {
  std::array<double, 3> low{}, high{};
  std::array<std::uint32_t, 3> n{};
  double volume = 0;
  std::uint32_t locate(const double *x) const noexcept;   // native formula; kNone outside
};

struct Configuration {
  std::string preset = "custom";
  bool feedback = false;                  // CellTargets from the local gas state
  Transport transport = Transport::none;  // deposit transport with the carriers
  std::vector<Species> species;           // owner species index = chemistry registry index first
  std::vector<Process> processes;         // chemistry process table (blob order)
  std::string chemistry_blob_sha256, chemistry_registry_sha256;
  Grid grid;
  std::uint64_t max_entries = 1ull << 20;
  // runtime only (not in the identity)
  Execution execution = Execution::host;
  std::string device_uuid, nvrtc_library;
};
// Presets (README 7, INTEGRATION.txt 4):
//   "production-gas-feedback": feedback on, carrier_shares transport (historical preset name; unqualified)
//   "production-gas-ledger":   feedback off, carrier_shares transport (chemistry results identical to the
//                              field-target approximation; closure requires independent qualification)
// Returns base with the preset fields set and the preset name bound.
Configuration preset(const std::string &name, Configuration base);
std::string configuration_sha256(const Configuration &);

class Refusal : public std::runtime_error {
public:
  enum class Kind { configuration, input, exactness, capacity, transaction, restore, device };
  Refusal(Kind k, const std::string &what) : std::runtime_error("gas_target_owner_v1/" + what), kind_(k) {}
  Kind kind() const noexcept { return kind_; }

private:
  Kind kind_;
};

// Correctly rounded value and the normalised forms digits (134, hex text
// "d133:...:d0" of nonzero span) of one exact quantity.
struct ExactValue {
  double value = 0;
  std::string digits;
};
struct FeedbackReceipt {
  std::uint64_t slots = 0, eligible = 0, changed = 0, clamped = 0;
  std::string cells_before_sha256, cells_after_sha256;
};
struct EntryValue {
  std::uint32_t cell = 0, species = 0;
  std::array<double, 5> value{};   // E (J), Px, Py, Pz (kg m/s), U (J)
};
struct Ledgers {
  std::string state_sha256;
  std::array<std::uint64_t, pfg::C_COUNTERS> counters{};
  std::vector<std::pair<std::string, ExactValue>> fixed;            // G_* by name
  std::vector<std::pair<std::string, ExactValue>> per_species;      // flow/<kind>/<state>/<q>, chg/..., wall/...
  std::vector<std::pair<std::string, ExactValue>> per_process;      // proc/<index>/<q>
  std::vector<EntryValue> entries;                                  // nonzero, ascending key
  std::array<ExactValue, 5> entry_totals;                           // exact sums over the entries
};

// Exact accumulator on the forms digit grid (host).
class Exact {
public:
  Exact() { d_.fill(0); }
  void add(const pfg::Chunks &c, int negate = 0) {
    if (c.bad) throw Refusal(Refusal::Kind::exactness, "exactness: unrepresentable exact term");
    pfg::add_chunks(d_.data(), c, negate);
  }
  void add_product(double a, double b, int negate = 0) { add(pfg::chunks_product(a, b), negate); }
  void add_triple(double a, double b, double c, int negate = 0) { add(pfg::chunks_triple(a, b, c), negate); }
  void add_digits(const std::int64_t *d, int negate = 0) {
    for (int k = 0; k < pfg::kDigits; ++k)
      d_[k] += negate ? -d[k] : d[k];
  }
  void merge(const Exact &o, int negate = 0) { add_digits(o.d_.data(), negate); }
  std::array<std::int64_t, pfg::kDigits> normalised() const {
    auto x = d_;
    pf7::normalise(x.data());
    return x;
  }
  double value() const {
    auto x = d_;
    pf7_u32 st = 0;
    const double result = pf7::round_sum(x.data(), &st);
    if (st) throw Refusal(Refusal::Kind::exactness, "exactness: rounded exact value overflow");
    return result;
  }
  bool zero() const {
    auto x = normalised();
    return pfg::digits_zero(x.data(), pfg::kDigits);
  }
  bool operator==(const Exact &o) const { return normalised() == o.normalised(); }
  const std::int64_t *data() const { return d_.data(); }

private:
  std::array<std::int64_t, pfg::kDigits> d_;
};
std::string digits_text(const std::int64_t *normalised);
ExactValue exact_value(const std::int64_t *digits);

// Partial accumulation of a subset of inputs over the candidate's key table
// (rank-count independence: any split of the inputs, merged in any order,
// gives the same state).
struct Partial {
  std::vector<std::int64_t> entry, global;   // digits (not normalised)
  std::array<std::uint64_t, pfg::C_COUNTERS> counter{};
  std::uint64_t keys = 0;                    // key-table size it was built on
  std::uint64_t accepted_step = 0, stage = 0;
  std::string key_table_sha256, configuration_sha256, payload_sha256, candidate_sha256;
  std::vector<std::string> input_sha256;
};

// Census of a population (host utility): exact count, kinetic energy and
// momentum per (cell, gas species) and per species. Rows of non-gas species
// are ignored in the sums, but every row must have a unique ID, a known
// species, positive finite weight, and finite phase space inside the grid.
struct CensusRow {
  std::uint64_t native_id = 0;
  std::uint32_t species = 0, pad = 0;
  double position_m[3]{}, proper_velocity_m_s[3]{}, weight = 0;
};
struct CensusCell {
  std::uint32_t cell = 0, species = 0;
  Exact count, kinetic, momentum[3];
};
struct Census {
  std::vector<CensusCell> cells;                       // ascending (cell, species)
  std::vector<std::array<Exact, 5>> species_totals;    // per species: count, K, P3
  const CensusCell *find(std::uint32_t cell, std::uint32_t species) const;
};
Census census(const Configuration &, std::span<const CensusRow> rows);
// Movers of a motion phase: rows before and after (any order; joined by
// native ID), the IDs sunk at walls. Gas species only; cell_count from the
// before census. After and sunk must partition the before inventory exactly;
// motion preserves species and binary weight. Birth/debit/release/merge/split
// operators must be represented separately. Output in ascending native ID.
std::vector<MoverRecord> movers(const Configuration &, std::span<const CensusRow> before,
                                std::span<const CensusRow> after, std::span<const std::uint64_t> sunk);

class DeviceLane;   // gas_device.hpp

class Owner {
public:
  // Rollback authority for one begin attempt, including its mutated, failed or
  // prepared state. Unlike CandidateBinding this intentionally survives changes
  // within that attempt. It cannot authorize reads or publication and is never
  // serialized. Calls on one Owner must remain serial.
  class AttemptBinding {
    friend class Owner;
    std::shared_ptr<const void> lifetime_;
    std::uint64_t generation_ = 0, step_ = 0;
  public:
    AttemptBinding() = default;
    AttemptBinding(const AttemptBinding &) = delete;
    AttemptBinding &operator=(const AttemptBinding &) = delete;
    AttemptBinding(AttemptBinding &&) noexcept = default;
    AttemptBinding &operator=(AttemptBinding &&) noexcept = default;
  };
  // If initialization fails after opening this attempt, roll it back before
  // propagating the error. A pre-existing attempt is never discarded.
  AttemptBinding begin_attempt(std::uint64_t accepted_step);
  bool matches_attempt(const AttemptBinding &) const noexcept;
  bool abandon_attempt(AttemptBinding &) noexcept;

  // Read-only binding to one open candidate lifetime and mutation revision.
  // This confers no mutation/publication authority; contents still need digest
  // binding. It is transient and is never encoded in checkpoints.
  class CandidateBinding {
    friend class Owner;
    std::shared_ptr<const void> lifetime_;
    std::uint64_t generation_ = 0, revision_ = 0, step_ = 0;
  public:
    CandidateBinding() = default;
    CandidateBinding(const CandidateBinding &) = delete;
    CandidateBinding &operator=(const CandidateBinding &) = delete;
    CandidateBinding(CandidateBinding &&) noexcept = default;
    CandidateBinding &operator=(CandidateBinding &&) noexcept = default;
    std::uint64_t accepted_step() const noexcept { return step_; }
  };
  CandidateBinding candidate_binding() const;
  bool matches_candidate(const CandidateBinding &) const noexcept;
  bool matches_candidate(const CandidateBinding &, std::uint64_t expected_step) const noexcept;

  // A capability for this exact live candidate, never a transferable digest.
  class PreparedCommit {
    friend class Owner;
    std::shared_ptr<const void> lifetime_;
    std::uint64_t generation_ = 0;
    std::string digest_;
  public:
    PreparedCommit() = default;
    PreparedCommit(const PreparedCommit &) = delete;
    PreparedCommit &operator=(const PreparedCommit &) = delete;
    PreparedCommit(PreparedCommit &&) noexcept = default;
    PreparedCommit &operator=(PreparedCommit &&) noexcept = default;
    const std::string &candidate_sha256() const noexcept { return digest_; }
  };
  explicit Owner(const Configuration &);
  ~Owner();
  Owner(const Owner &) = delete;
  Owner &operator=(const Owner &) = delete;

  const Configuration &configuration() const noexcept { return cfg_; }
  const std::string &configuration_sha256() const noexcept { return cfg_sha_; }
  pfg::Layout layout() const noexcept { return L_; }

  // ---- step transaction
  void begin_step(std::uint64_t accepted_step);
  // One chemistry stage: its admitted events' recoil records (canonical order)
  // and its target-population flows (debits; stage 2 releases). Returns the
  // candidate state digest (bound into native's stage receipt).
  std::string take(std::uint64_t accepted_step, std::uint64_t stage, std::span<const RecoilRecord>,
                   std::span<const FlowRecord>);
  // Flows outside the chemistry stages (births, sinks, transaction outputs).
  std::string apply_flows(std::uint64_t accepted_step, std::span<const FlowRecord>);
  // Deposit transport of one motion phase (carrier_shares; refused when the
  // configuration has transport none).
  std::string transport(std::uint64_t accepted_step, std::span<const MoverRecord>);
  // Feedback: rewrites drift and s of eligible slots in place (no-op when the
  // configuration has feedback off). Reads the candidate when a step is open.
  FeedbackReceipt feedback(std::uint64_t accepted_step, std::uint64_t stage, std::span<const std::uint32_t> atlas,
                           std::span<CellTargetsAbi> cells) const;
  void commit(const std::string &candidate_sha256);
  // All hashing, validation, compaction and allocation finish here. Success
  // freezes the candidate until publication or abandon. The common coordinator
  // must prepare and preflight every owner before publishing any owner.
  PreparedCommit prepare_commit(const std::string &candidate_sha256);
  bool can_publish(const PreparedCommit &) const noexcept;
  bool try_publish(PreparedCommit &) noexcept;
  void abandon() noexcept;
  bool step_open() const noexcept { return open_; }

  // ---- partial accumulation (rank independence)
  void prepare_keys(std::span<const RecoilRecord>, std::span<const FlowRecord>, std::span<const MoverRecord>);
  Partial accumulate(std::span<const RecoilRecord>, std::span<const FlowRecord>) const;
  std::string merge(std::uint64_t accepted_step, std::uint64_t stage, std::span<const Partial>);

  // ---- state
  std::string accepted_sha256() const;
  std::string candidate_sha256() const;
  std::uint64_t accepted_step() const noexcept;
  Ledgers ledgers(bool candidate = false) const;
  std::array<std::int64_t, pfg::kDigits> global_digits(std::uint32_t g, bool candidate = false) const;
  std::vector<EntryValue> entries(bool candidate = false) const;

  // ---- persistence: canonical codec of the accepted state (README 6)
  std::vector<std::uint8_t> encode() const;
  static std::unique_ptr<Owner> restore(std::span<const std::uint8_t>, const Configuration &,
                                        const std::string &expected_state_sha256);

  struct State;   // gas_owner.cpp

  // Exact sum over the entries of quantity q (E, Px, Py, Pz, U).
  Exact entry_total(std::uint32_t q, bool candidate = false) const;
  std::string device_json() const;
  double device_kernel_ms() const;

private:
  Configuration cfg_;
  std::string cfg_sha_;
  pfg::Layout L_{};
  std::vector<pfg::SpeciesInfo> species_;
  std::vector<pfg::ProcessInfo> procs_;
  std::unique_ptr<State> acc_, cand_;
  bool open_ = false;
  std::shared_ptr<const void> lifetime_ = std::make_shared<const int>(0);
  std::uint64_t generation_ = 0, candidate_revision_ = 0;
  bool publication_ready_ = false;
  bool failed_ = false;  // a failed mutating operation requires abandon
  std::uint64_t next_stage_ = 1;
  bool partial_batch_ready_ = false;
  std::vector<std::string> partial_batch_inputs_;
  std::uint64_t open_step_ = 0;
  std::unique_ptr<DeviceLane> dev_;
  bool host() const noexcept;
  void ensure_open(std::uint64_t step, const char *what) const;
  void invalidate_candidate_bindings();
  void insert_keys(State &, std::vector<std::uint64_t> keys);
  void prepare_keys_impl(std::span<const RecoilRecord>, std::span<const FlowRecord>, std::span<const MoverRecord>, bool bind);
  Partial accumulate_impl(std::span<const RecoilRecord>, std::span<const FlowRecord>, bool bind) const;
  std::string finish(State &, const char *what);
  std::string run(std::uint64_t stage, std::span<const RecoilRecord>, std::span<const FlowRecord>, const char *what);
  std::array<std::uint64_t, pfg::C_COUNTERS> device_apply(std::span<const RecoilRecord>, std::span<const FlowRecord>,
                                                         std::span<const MoverRecord>,
                                                         const std::array<std::uint64_t, pfg::C_COUNTERS> *host);
  std::string digest(const State &, bool candidate) const;
  const State &view(bool candidate, State &tmp) const;
  void verify_device(const char *what, const std::string &host_digest) const;
};

}  // namespace pf_gas_v1
