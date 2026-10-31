#pragma once
#include "raikou/exact_particle_count.h"
#include "raikou/json.hpp"
#include "raikou/mfem_gap_pic.h"
#include <map>
#include <memory>
#include <optional>
#include <set>

namespace raikou::mfem_backend {
namespace canonical {
struct BlockDigest;
}
class PrivateFullTimeStep;
class FullTimeAccess;
struct FullTimeStepBinding;
struct FullTimeGenesisBinding;
struct FullTimeAuthoredPacket;
struct FullTimeAcceptedBoundary;
struct FullTimeCandidateExposures;
struct FullTimeConfiguration;
struct FullTimeHistory;
struct FullTimeHistoryStatus;
using ArenaJson = nlohmann::json;
// Structural sharing: immutable shared payload, copy-on-write on mutation.
// Copying a snapshot copies one pointer per container; the first mutation of
// a shared payload copies it once. The payload carries a lazily cached
// canonical digest that every mutation clears.
template <class T> class SharedValue {
public:
  SharedValue() : p_(empty()) {}
  SharedValue(T value) : p_(std::make_shared<Payload>(std::move(value))) {}
  const T &get() const noexcept { return p_->value; }
  T &mutate() {
    if (p_.use_count() != 1)
      p_ = std::make_shared<Payload>(p_->value);
    p_->digest.clear();
    return p_->value;
  }
  bool shares(const SharedValue &other) const noexcept {
    return p_ == other.p_;
  }
  // Keeps the immutable payload alive (memo keys cannot be reused).
  std::shared_ptr<const void> identity() const noexcept { return p_; }
  const std::string &cached_digest() const noexcept { return p_->digest; }
  void cache_digest(std::string d) const { p_->digest = std::move(d); }

private:
  struct Payload {
    explicit Payload(T v) : value(std::move(v)) {}
    T value;
    mutable std::string digest;
  };
  static const std::shared_ptr<Payload> &empty() {
    static const auto e = std::make_shared<Payload>(T{});
    return e;
  }
  std::shared_ptr<Payload> p_;
};
template <class K, class V> class SharedMap : public SharedValue<std::map<K, V>> {
public:
  using map_type = std::map<K, V>;
  using const_iterator = typename map_type::const_iterator;
  using iterator = const_iterator;
  using size_type = typename map_type::size_type;
  using SharedValue<map_type>::SharedValue;
  using SharedValue<map_type>::get;
  const V &at(const K &k) const { return get().at(k); }
  size_type count(const K &k) const { return get().count(k); }
  const_iterator find(const K &k) const { return get().find(k); }
  const_iterator begin() const noexcept { return get().begin(); }
  const_iterator end() const noexcept { return get().end(); }
  size_type size() const noexcept { return get().size(); }
  bool empty() const noexcept { return get().empty(); }
  template <class... A> std::pair<const_iterator, bool> emplace(A &&...a) {
    auto r = this->mutate().emplace(std::forward<A>(a)...);
    return {const_iterator(r.first), r.second};
  }
  size_type erase(const K &k) {
    return get().count(k) ? this->mutate().erase(k) : 0;
  }
};
template <class K> class SharedSet : public SharedValue<std::set<K>> {
public:
  using set_type = std::set<K>;
  using const_iterator = typename set_type::const_iterator;
  using iterator = const_iterator;
  using size_type = typename set_type::size_type;
  using SharedValue<set_type>::SharedValue;
  using SharedValue<set_type>::get;
  size_type count(const K &k) const { return get().count(k); }
  const_iterator begin() const noexcept { return get().begin(); }
  const_iterator end() const noexcept { return get().end(); }
  size_type size() const noexcept { return get().size(); }
  bool empty() const noexcept { return get().empty(); }
  std::pair<const_iterator, bool> insert(const K &k) {
    if (auto i = get().find(k); i != get().end())
      return {i, false};
    auto r = this->mutate().insert(k);
    return {const_iterator(r.first), r.second};
  }
  bool operator==(const SharedSet &o) const {
    return this->shares(o) || get() == o.get();
  }
};
// Immutable exact dyadic gas inventory. Surface/solid inventory has one
// separate owner and is transferred exactly; no rounded conversion is implied
// here.
using ArenaCount = ExactParticleCount;
struct ArenaStateSpec {
  std::string id;
  int hydrogen = 0, oxygen = 0, charge_e = 0;
  double mass_kg = 0, internal_J = 0, reference_J = 0;
};
struct ArenaChannel {
  std::string id, primary_state, target_state;
  std::vector<std::string> product_states;
  std::string dataset_sha256, kernel_manifest_sha256,
      kernel_binding_file_sha256, model;
};
struct ArenaPacket {
  std::string id, state, origin, last_event;
  ArenaCount count;
  std::vector<std::uint64_t> native_ids;
};
struct FullTimePhaseKey;
struct ArenaPhase {
  enum class Kind {
    InitialFieldProof,
    ProvisionalEqualDtMidX,
    FullTimeEndpoint
  };
  Kind kind = Kind::InitialFieldProof;
  std::uint64_t accepted_step = 0;
  double physical_stage_s = 0, previous_dt_s = 0, next_dt_s = 0;
  std::string state_prefix, mesh_prefix, clock_prefix;
  std::shared_ptr<const FullTimePhaseKey> full_time_key;
};
struct ArenaProvenance {
  enum class Kind { ManufacturedForced, ConditionalPhysicalData };
  Kind kind = Kind::ManufacturedForced;
  std::string dataset_sha256, kernel_manifest_sha256,
      kernel_binding_file_sha256;
  std::string model, selection_method, selection_state_prefix;
  // Collision subflow waiting clock. Never advances native trajectory time.
  double collision_elapsed_s = 0;
};
struct ArenaProduct {
  enum class At { Primary, Target };
  std::string id, state;
  PicVector proper_velocity_m_s{};
  At location = At::Target;
};
struct ArenaEvent {
  std::string id, channel, primary, target, untouched_primary_id;
  double event_weight = 0;
  // A scattered primary is one of products and is identified explicitly.
  // All products have event_weight; untouched reactant remainders do not.
  std::string scattered_primary_id;
  std::vector<ArenaProduct> products;
  ArenaPhase phase;
  ArenaProvenance provenance;
};
struct ArenaRateEntry {
  std::string packet_id, state;
  GapParticle native;
};
struct ArenaStochasticCompletion {
  std::string integration_owner_id, before_rng_clock_sha256,
      after_rng_clock_sha256;
  std::string common_clock_sha256;
  std::string channel_id, dataset_sha256, kernel_manifest_sha256,
      kernel_binding_file_sha256, model;
  std::uint64_t stage_sequence = 0;
  double frozen_physical_stage_s = 0, common_interval_begin_s = 0,
         common_interval_end_s = 0;
};
struct ArenaSnapshot {
  ArenaSnapshot() = default;
  ArenaSnapshot(const ArenaSnapshot &);
  ArenaSnapshot &operator=(const ArenaSnapshot &) = delete;
  ~ArenaSnapshot() noexcept;
  SharedMap<std::string, ArenaPacket> packets;
  SharedSet<std::string> ever_packet_ids, event_ids;
  // Legacy JSON journal. Always empty in full-time mode, whose bounded typed
  // history lives in full_time_history (anchors + unanchored window).
  ArenaJson journal = ArenaJson::array(), sinks = ArenaJson::array();
  ArenaJson stochastic_completions = ArenaJson::object();
  ArenaJson full_time_context;
  // Full-time: cumulative canonical row bytes ever appended (never reset).
  std::size_t journal_bytes = 0;
  std::string chain;
  // Full-time: revision = cumulative rows ever appended (never reset).
  std::uint64_t revision = 0, accepted_steps = 0;
  double previous_dt_s = 0;
  // Storage only: excluded from the logical snapshot digest, so anchoring
  // cadence never changes an accepted state identity.
  std::shared_ptr<const FullTimeHistory> full_time_history;
  // Lazily cached canonical digest of this immutable snapshot (full-time).
  mutable std::string digest_cache;
};
class PreparedArena {
public:
  PreparedArena(PreparedArena &&) noexcept;
  PreparedArena &operator=(PreparedArena &&) noexcept;
  ~PreparedArena();
  PreparedArena(const PreparedArena &) = delete;
  PreparedArena &operator=(const PreparedArena &) = delete;
  const ArenaJson &audit() const;
  std::string candidate_state_prefix() const;
  std::optional<ArenaStochasticCompletion>
  candidate_stochastic_completion(const std::string &owner) const;
  std::shared_ptr<const FullTimeCandidateExposures>
  candidate_full_time_exposures() const;
  std::shared_ptr<const FullTimeAcceptedBoundary>
  candidate_full_time_boundary() const;

private:
  friend class PicArena;
  friend class FullTimeAccess;
  friend class CollisionScope;
  friend class PrivateFullTimeStep;
  struct Impl;
  explicit PreparedArena(std::unique_ptr<Impl>);
  std::unique_ptr<Impl> impl_;
};
// Exclusive-owner prototype; no own collision selection, RNG, cross sections,
// external energy reservoir or independent neutral trajectories.
class PicArena {
public:
  struct NativeState;
  PicArena(GapPic &, std::vector<ArenaStateSpec>, std::vector<ArenaChannel>,
           std::vector<ArenaPacket>, std::string initial_coupling_prefix);
  PicArena(const PicArena &) = delete;
  PicArena &operator=(const PicArena &) = delete;
  PicArena(PicArena &&) = delete;
  PicArena &operator=(PicArena &&) = delete;
  std::vector<ArenaRateEntry>
  rate_population(std::size_t maximum_entries) const;
  // Default configuration: 100000 carriers, default budget and history policy.
  static std::unique_ptr<PicArena> create_full_time(
      GapPicOptions, std::vector<ArenaStateSpec>, std::vector<ArenaChannel>,
      std::vector<FullTimeAuthoredPacket>, const FullTimeGenesisBinding &);
  static std::unique_ptr<PicArena> create_full_time(
      GapPicOptions, std::vector<ArenaStateSpec>, std::vector<ArenaChannel>,
      std::vector<FullTimeAuthoredPacket>, const FullTimeGenesisBinding &,
      const FullTimeConfiguration &);
  // Bounded history (v5). Status/export never mutate the owner.
  FullTimeHistoryStatus full_time_history_status() const;
  // Serialized anchor chain: text tag, u64 count, then each anchor's
  // canonical bytes as text. Disk restore from it is not implemented.
  std::string export_full_time_anchor_chain() const;
  // Revalidates the anchor chain and the unanchored window from the last
  // anchor (grammar, hash chain, counters, RNG chain, owner seals).
  void validate_full_time_history() const;
  std::size_t full_time_carrier_cap() const noexcept { return carrier_cap_; }
  static std::unique_ptr<PicArena>
  restore_full_time(const std::string &path, GapPicOptions,
                    std::vector<ArenaStateSpec>, std::vector<ArenaChannel>,
                    const FullTimeAcceptedBoundary &expected);
  std::shared_ptr<const FullTimeAcceptedBoundary> full_time_boundary() const;
  PrivateFullTimeStep begin_full_time_step(const FullTimeStepBinding &,
                                           std::size_t maximum_carriers,
                                           std::size_t maximum_records);
  std::shared_ptr<const ArenaSnapshot> snapshot() const noexcept {
    return accepted_;
  }
  // Kinematics view has rounded diagnostic weight; do not select collisions
  // from it after compensated debits. Use positive rate_carriers instead.
  const GapPicOptions &native_options() const noexcept {
    return pic_->options();
  }
  GapParticle view(const std::string &) const;
  std::vector<GapParticle> rate_carriers(const std::string &) const;
  std::string state_prefix() const;
  PreparedArena prepare(const ArenaEvent &, int manufactured_fault_point = 0);
  PreparedArena
  prepare_compound(std::size_t maximum_reactions,
                   const std::function<void(PicArena &)> &builder);
  PreparedArena prepare_compound_with_completion(
      std::size_t maximum_reactions,
      const std::function<ArenaStochasticCompletion(PicArena &)> &builder);
  std::optional<ArenaStochasticCompletion>
  stochastic_completion(const std::string &owner) const;
  const ArenaChannel &channel_binding(const std::string &id) const {
    return channels_.at(id);
  }
  PreparedArena prepare_transport(double dt_s,
                                  int manufactured_fault_point = 0);
  bool can_publish(const PreparedArena &) const noexcept;
  bool try_publish(PreparedArena &) noexcept;
  bool try_rollback(PreparedArena &) noexcept;
  void save(const std::string &) const;
  static std::unique_ptr<PicArena>
  restore(GapPic &, const std::string &, const std::vector<ArenaStateSpec> &,
          const std::vector<ArenaChannel> &,
          const std::string &expected_initial_coupling_prefix);
  static std::vector<double> positive_carriers(ArenaCount);
  static ArenaCount debit(ArenaCount, double);
  static bool exact_count_equal(ArenaCount, ArenaCount);
  static std::string fingerprint(const GapPic &);
  static const char *source_identity() noexcept;

private:
  friend class FullTimeAccess;
  friend class CollisionScope;
  friend class PrivateFullTimeStep;
  struct ReactionBuild;
  ReactionBuild build_reaction(const ArenaEvent &, const ArenaSnapshot &,
                               const NativeState &, const std::string &,
                               bool) const;
  PicArena(GapPic &, std::vector<ArenaStateSpec>, std::vector<ArenaChannel>,
           std::vector<ArenaPacket>, std::string, bool,
           std::size_t carrier_cap = 100000);
  void save_full_time(const std::string &) const;
  std::shared_ptr<void> full_time_lifetime_ = std::make_shared<int>(0);
  std::size_t carrier_cap_ = 100000, identity_cap_ = 1u << 22;
  std::uint64_t memory_budget_bytes_ = 0;
  std::shared_ptr<const FullTimeConfiguration> configuration_;
  // Exact-count validation memo: (shared packets payload, SHA-256 of the
  // ascending (native ID, weight) sequence) pairs already proven exact.
  mutable std::vector<std::pair<std::shared_ptr<const void>, std::string>>
      count_memo_;
  std::vector<ArenaRateEntry> full_time_census(const ArenaSnapshot &,
                                               const std::vector<GapParticle> &,
                                               std::size_t cap,
                                               bool entries) const;
  std::string registry_v5_sha_, options_v5_sha_;
  std::shared_ptr<const canonical::BlockDigest> accepted_native_digest_;
  std::weak_ptr<void> active_full_time_outer_;
  bool full_time_builder_active_ = false;
  std::unique_ptr<GapPic> owned_pic_;
  std::shared_ptr<const FullTimeAcceptedBoundary> full_time_boundary_;
  GapPic *pic_;
  std::map<std::string, ArenaStateSpec> states_;
  std::map<std::string, ArenaChannel> channels_;
  std::shared_ptr<const ArenaSnapshot> accepted_, genesis_;
  std::vector<GapParticle> genesis_native_;
  std::uint64_t genesis_next_id_ = 0;
  double genesis_time_s_ = 0;
  void validate_history() const;
  std::string initial_prefix_, registry_sha_;
  std::uint64_t scope_counter_ = 0, active_scope_ = 0,
                compound_base_revision_ = 0;
  std::size_t compound_limit_ = 0;
  bool completion_compound_ = false;
  std::string compound_base_chain_;
  PreparedArena compound_impl(std::size_t,
                              const std::function<void(PicArena &)> &, bool);
  void complete_stochastic(const ArenaStochasticCompletion &);
  static NativeState capture(const GapPic &);
  static void swap_native(GapPic &, NativeState &) noexcept;
  static bool matches(const GapPic &, const NativeState &) noexcept;
  static std::string sha256(const std::string &);
  ArenaJson registry() const;
  void validate_population(const ArenaSnapshot &) const;
  PreparedArena stage(std::shared_ptr<ArenaSnapshot>, NativeState &, ArenaJson,
                      int);
};
} // namespace raikou::mfem_backend
