#pragma once
#include "electrode_heat_feedback.hpp"
#include "raikou/json.hpp"
#include "raikou/mfem_gap_pic.h"
#include <map>
#include <memory>
#include <optional>
#include <set>

namespace raikou::mfem_backend {
class PrivateFullTimeStep;
struct FullTimeStepBinding;
struct FullTimeGenesisBinding;
struct FullTimeAuthoredPacket;
struct FullTimeAcceptedBoundary;
struct FullTimeCandidateExposures;
using ArenaJson = nlohmann::json;
// Reuses the heat owner's public inventory type. This is the gas inventory
// after an accepted surface transfer, not a second owner of solid inventory.
using ArenaCount = electrode_heat::Count;
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
  enum class Kind { InitialFieldProof, ProvisionalEqualDtMidX, FullTimeEndpoint };
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
  std::map<std::string, ArenaPacket> packets;
  std::set<std::string> ever_packet_ids, event_ids;
  ArenaJson journal = ArenaJson::array(), sinks = ArenaJson::array();
  ArenaJson stochastic_completions = ArenaJson::object();
  std::size_t journal_bytes = 0;
  std::string chain;
  std::uint64_t revision = 0, accepted_steps = 0;
  double previous_dt_s = 0;
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
  std::optional<ArenaStochasticCompletion> candidate_stochastic_completion(const std::string& owner) const;
  std::shared_ptr<const FullTimeCandidateExposures> candidate_full_time_exposures() const;

private:
  friend class PicArena;
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
  static std::unique_ptr<PicArena> create_full_time(GapPicOptions,
      std::vector<ArenaStateSpec>,std::vector<ArenaChannel>,
      std::vector<FullTimeAuthoredPacket>,const FullTimeGenesisBinding&);
  static std::unique_ptr<PicArena> restore_full_time(const std::string& path,
      GapPicOptions,std::vector<ArenaStateSpec>,std::vector<ArenaChannel>,
      const FullTimeAcceptedBoundary& expected);
  std::shared_ptr<const FullTimeAcceptedBoundary> full_time_boundary() const;
  PrivateFullTimeStep begin_full_time_step(const FullTimeStepBinding&,
      std::size_t maximum_carriers,std::size_t maximum_records);
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
