#pragma once
#include "parent_coverage_v2_scope.hpp"

// Parent-ionization subflow with coverage scope v2. Derived from the coverage
// v1 subflow (sealed v3 native_parent_subflow.hpp included unchanged for
// Limits, Binding, Interval, Grid and sha). Differences from v1:
//  - one composite channel: the sealed Song bridge sigma for E <= 1 keV (same
//    calls and values as v1), the dataset outer3 BEB parent rate on
//    (1, 10] keV and the provisional RBEB parent rate above 10 keV;
//  - events above 1 keV use the BEB/RBEB SDCS kernel and three product-RNG
//    variates (orbital, ejected energy, azimuth); Song-band events use the
//    sealed v1 kernel and one isotropic axis exactly as v1;
//  - per-category applied/executed measures in the exact ledger; event ids
//    carry the category (model-derived events: "model-derived/...");
//  - the registered channel binding is the composite v2 identity when the
//    scope has high-energy bands, and the v1 registration otherwise.
// With no pair above 1 keV the candidate list, scheduler draws, product RNG,
// event identities and products are those of v1 (and v3).
namespace raikou::native_water::coverage_v2 {
using ::raikou::native_water::Binding;
using ::raikou::native_water::Grid;
using ::raikou::native_water::Interval;
using ::raikou::native_water::Limits;

struct State {
  std::string scheduler, product_rng;
  std::uint64_t next_event = 1, next_interval = 0;
  Ledger coverage;
};
class ParentSubflow;
class Prepared {
  friend class ParentSubflow;
  std::shared_ptr<const char> lifetime;
  std::uint64_t generation = 0;
  std::unique_ptr<m::PreparedArena> arena;
  std::unique_ptr<State> after;
  Json audit_;
public:
  const Json& audit() const { return audit_; }
};
class PreparedParentStage {
  friend class ParentSubflow;
  std::shared_ptr<const char> lifetime_;
  std::uint64_t generation_ = 0;
  m::FullTimeSlot slot_ = m::FullTimeSlot::before_transport;
  std::string common_binding_, common_clock_, accepted_arena_, before_hash_, after_hash_, private_prefix_;
  std::string exposure_binding_;
  std::shared_ptr<const m::FullTimeExposureKey> exposure_key_;
  std::unique_ptr<State> after_;
  Json audit_;
public:
  PreparedParentStage() = default;
  PreparedParentStage(PreparedParentStage&&) noexcept = default;
  PreparedParentStage& operator=(PreparedParentStage&&) noexcept = default;
  PreparedParentStage(const PreparedParentStage&) = delete;
  PreparedParentStage& operator=(const PreparedParentStage&) = delete;
  const Json& audit() const { return audit_; }
};
class PreparedParentState {
  friend class ParentSubflow;
  std::shared_ptr<const char> lifetime_;
  std::uint64_t generation_ = 0;
  std::string native_prefix_, after_hash_, common_clock_;
  std::unique_ptr<State> after_;
  Json coverage_;
public:
  PreparedParentState() = default;
  PreparedParentState(PreparedParentState&&) noexcept = default;
  PreparedParentState& operator=(PreparedParentState&&) noexcept = default;
  PreparedParentState(const PreparedParentState&) = delete;
  PreparedParentState& operator=(const PreparedParentState&) = delete;
  // The sealed step's unresolved-coverage ledger (both exposures, totals,
  // cumulative), for the coordinator's per-step ledgers.
  const Json& coverage() const { return coverage_; }
};
class ParentSubflow {
  m::PicArena& arena_ref_;
  const Grid grid_;
  const Binding binding_;
  const Limits limits_;
  const w::ParentBridge bridge_;
  const CoverageScope scope_;
  const ArenaBinding arena_;
  const bool full_time_;
  std::string integration_config_;
  std::unordered_map<std::uint64_t, double> volumes_;
  std::unique_ptr<State> accepted_;
  const std::shared_ptr<const char> lifetime_ = std::make_shared<const char>(0);
  std::uint64_t generation_ = 0;
  static constexpr const char* owner_id = "native_water_parent_subflow/coverage-v2";
  static constexpr const char* schema_legacy = "native-water-parent-subflow-coverage-v2";
  static constexpr const char* schema_full_time = "native-water-parent-full-time-coverage-v2";
  Json payload(const State& s) const {
    return {{"scheduler", s.scheduler}, {"product_rng", s.product_rng}, {"next_event", s.next_event},
            {"next_interval", s.next_interval}, {"coverage", scope_.ledger_json(s.coverage)}};
  }
  static std::string configuration_identity(const Grid& grid, const Binding& binding, const Limits& limits,
                                            const CoverageScope& scope) {
    return sha(Json{{"mode", "full-time-parent-subflow-coverage-v2"}, {"grid", grid.identity},
      {"dataset", binding.dataset_sha256}, {"kernel", binding.kernel_sha256}, {"binding_file", binding.kernel_binding_file_sha256},
      {"arena_channel_binding", scope.arena_binding(binding.kernel_binding_file_sha256).json()},
      {"limits", {limits.reactions, limits.carriers, limits.cells, limits.pairs, limits.queries, limits.candidates}},
      {"coverage_scope", scope.identity()}}.dump());
  }
  static std::string hash_with(const std::string& configuration, const Json& payload, bool full_time) {
    if (!full_time) return sha(Json{{"configuration", configuration}, {"state", payload}}.dump());
    return sha(Json{{"configuration", configuration}, {"state", payload}, {"full_time", true}}.dump());
  }
  // Unlike v3, legacy mode also binds the configuration (and so the coverage
  // scope and policy) into every RNG/clock state hash.
  std::string state_hash(const State& s) const { return hash_with(integration_config_, payload(s), full_time_); }
  template<class Arena>
  void check_completion_on(const Arena& arena, const State& s) const {
    auto prior = arena.stochastic_completion(owner_id);
    if (prior) {
      pc::require(prior->stage_sequence != UINT64_MAX && prior->stage_sequence + 1 == s.next_interval &&
        prior->after_rng_clock_sha256 == state_hash(s), "RNG/clock state does not match accepted native exposure");
    } else pc::require(s.next_interval == 0, "missing accepted native exposure history");
  }
  void check_completion(const State& s) const {
    check_completion_on(arena_ref_, s);
    pc::require(!s.coverage.pending, "accepted parent state has an unpublished exposure ledger");
    if (full_time_) {
      const auto boundary = arena_ref_.full_time_boundary();
      pc::require(boundary && boundary->rng_sha256 == state_hash(s) &&
        boundary->accepted_arena_sha256 == arena_ref_.state_prefix() &&
        boundary->accepted_arena_revision == arena_ref_.snapshot()->revision &&
        boundary->accepted_step == arena_ref_.snapshot()->accepted_steps,
        "accepted full-time RNG/configuration seal differs from native boundary");
    }
  }
  static std::string common_binding(const m::FullTimeStepBinding& b) {
    return sha(Json{{"step", b.accepted_step}, {"revision", b.accepted_arena_revision},
      {"begin", b.begin_s}, {"end", b.end_s}, {"dt", b.dt_s}, {"arena", b.accepted_arena_sha256},
      {"field", b.accepted_field_sha256}, {"mesh", b.mesh_sha256}, {"material", b.material_sha256},
      {"clock", b.common_clock_sha256}}.dump());
  }
  static std::string exposure_binding(const m::ExposureBinding& e) {
    return sha(Json{{"common", common_binding(e.common)},
      {"slot", e.slot == m::FullTimeSlot::before_transport ? "before_transport" : "after_transport"},
      {"physical_stage", e.frozen_physical_stage_s}, {"duration", e.collision_duration_s},
      {"sequence", e.owner_stage_sequence}, {"field", e.frozen_field_sha256}}.dump());
  }
  static std::string encode(const std::mt19937_64& rng) { std::ostringstream s;s << rng;return s.str(); }
  static std::mt19937_64 decode(const std::string& text) {
    pc::require(text.size() < 10000, "product RNG checkpoint cap");
    std::istringstream s(text);std::mt19937_64 r;
    pc::require(bool(s >> r), "product RNG checkpoint");s >> std::ws;
    pc::require(s.eof(), "product RNG trailing bytes");return r;
  }
  static double uniform(std::mt19937_64& rng) { return (static_cast<double>(rng() >> 12) + .5) * 0x1p-52; }
  static pc::Vec direction(std::mt19937_64& rng) {
    const double z = 2 * uniform(rng) - 1, phi = 6.283185307179586476925286766559 * uniform(rng);
    const double r = std::sqrt((1 - z) * (1 + z));return {r * std::cos(phi), r * std::sin(phi), z};
  }
  struct Population {
    std::vector<pc::Packet> electrons, water;
    std::vector<std::uint64_t> electron_ids, water_ids;
    std::map<std::string, std::string> parent;
    std::map<std::string, pc::Packet> by_id;
  };
  template<class Arena>
  Population population(const Arena& a) const {
    Population p;
    for (const auto& entry : a.rate_population(limits_.carriers)) {
      if (entry.state != "electron" && entry.state != "H2O") continue;
      const auto& carrier = entry.native;
      const auto name = "native-carrier/" + std::to_string(carrier.id);
      pc::Packet q{name, entry.state, grid_.locate(carrier.position_m), carrier.weight,
                   carrier.mass_kg, carrier.proper_velocity_m_s};
      pc::require(p.parent.emplace(name, entry.packet_id).second, "duplicate native rate carrier");
      p.by_id.emplace(name, q);
      if (entry.state == "electron") { p.electrons.push_back(q);p.electron_ids.push_back(carrier.id); }
      else { p.water.push_back(q);p.water_ids.push_back(carrier.id); }
    }
    return p;
  }
  // The sealed bridge sigma is called unchanged unless the scope classifies
  // the energy as unresolved; each exclusion is counted and must equal the
  // ledger pass count for the same population.
  pc::Channel covered_channel(std::uint64_t* excluded) const {
    return scope_.channel(bridge_, excluded, binding_.kernel_binding_file_sha256);
  }
  struct Built {
    std::unique_ptr<m::PreparedArena> arena;
    std::unique_ptr<State> after;
    Json audit;
  };
  template<class Arena>
  Built build_events(Arena& arena, const Interval& interval, const State& before,
                     const m::ExposureBinding* exposure = nullptr) const {
    pc::require(before.next_interval < UINT64_MAX - 1,
      "collision exposure would leave an unrestorable exhausted interval counter");
    auto clock = std::unique_ptr<pc::Scheduler>(new pc::Scheduler(pc::Scheduler::restore(before.scheduler)));
    auto rng = decode(before.product_rng);
    auto next = std::make_unique<State>(before);
    Json events = Json::array();std::size_t pairs = 0, queries = 0;
    double elapsed = 0;bool terminal = false;
    const auto before_hash = state_hash(before);
    const Slot slot = !exposure ? Slot::legacy
      : (exposure->slot == m::FullTimeSlot::before_transport ? Slot::before_transport : Slot::after_transport);
    // Identical to the duration each scheduler path uses internally.
    const double duration = exposure ? exposure->collision_duration_s : interval.end_s - interval.begin_s;
    ExposureAccumulator ledger(scope_, slot, interval.sequence, duration);
    Json coverage_audit;
    auto token = arena.prepare_compound_with_completion(limits_.reactions, [&](auto& scoped) {
      while (!terminal) {
        const auto p = population(scoped);
        pc::CandidateInputLimits caps;caps.projectiles = caps.targets = limits_.carriers;caps.cells = limits_.cells;
        pc::CandidateWork work;std::uint64_t excluded = 0;
        auto candidates = pc::candidates(p.electrons, p.water, grid_.cells, {covered_channel(&excluded)},
          limits_.pairs - pairs, limits_.queries - queries, limits_.candidates, caps, &work);
        pairs += work.pair_visits;queries += work.channel_queries;
        const auto segment = ledger.segment(p.electrons, p.electron_ids, p.water, p.water_ids, volumes_, candidates);
        pc::require(segment.pair_visits == work.pair_visits && segment.rows.size() == excluded,
          "coverage ledger pass differs from the candidate query exclusions");
        const auto prefix = scoped.state_prefix();
        pc::Clock stage{interval.phase.accepted_step, scoped.snapshot()->revision, interval.begin_s,
          interval.end_s, prefix, grid_.identity, arena_.kernel_manifest_sha256, interval.phase.physical_stage_s, elapsed};
        const auto proposal = exposure
          ? clock->prepare_exposure({stage, exposure->common.dt_s, exposure->collision_duration_s}, candidates, limits_.candidates)
          : clock->prepare(stage, candidates, limits_.candidates);
        if (!proposal.event()) {
          pc::require(clock->try_publish(proposal), "terminal scheduler publication");
          ledger.close(segment, elapsed, duration);
          terminal = true;break;
        }
        pc::require(events.size() < limits_.reactions && next->next_event < UINT64_MAX - 1, "native reaction/event resource limit");
        const auto& c = *proposal.event();
        const auto category = scope_.category(c.energy_eV);
        Json extra;
        w::Products products;
        if (category == Category::song_table) {
          // The sealed v1 bridge validates against its own Song binding.
          auto song = c;
          song.data_sha256 = bridge_.channel().data_sha256;song.kernel_sha256 = bridge_.channel().kernel_sha256;
          products = bridge_.sample(song, p.by_id.at(c.projectile_id), p.by_id.at(c.target_id), direction(rng));
        } else {
          const auto band = category == Category::outer3_beb ? Band::outer3_beb : Band::provisional_rbeb;
          const long double u_orbital = uniform(rng), u_energy = uniform(rng), u_azimuth = uniform(rng);
          std::uint64_t ignored = 0;
          const auto s = high_energy_parent_sample(*scope_.bands(), band, covered_channel(&ignored), p.by_id.at(c.target_id).mass_kg, c,
                                                   p.by_id.at(c.projectile_id), p.by_id.at(c.target_id), u_orbital, u_energy, u_azimuth);
          pc::require(ignored == 0, "high-energy candidate re-validation excluded its own pair");
          products = s.products;
          extra = {{"category", name(category)}, {"basis", basis(category)}, {"model_derived", category == Category::provisional_rbeb},
                   {"kernel", sdcs_kernel}, {"orbital", s.rest.orbital_label},
                   {"ejected_energy_target_rest_eV", s.rest.ejected_eV},
                   {"scattered_energy_target_rest_eV", s.rest.primary_kinetic_J / pc::qe},
                   {"ion_recoil_energy_target_rest_eV", s.rest.ion_kinetic_J / pc::qe},
                   {"cos_theta_ejected", s.rest.cos_theta_ejected}};
        }
        const std::string id = event_prefix(category) + std::to_string(next->next_event++);
        m::ArenaEvent e;e.id = id;e.channel = c.channel_id;e.primary = p.parent.at(c.projectile_id);
        e.target = p.parent.at(c.target_id);e.event_weight = c.event_weight;e.scattered_primary_id = e.primary;
        const auto remaining = m::PicArena::debit(scoped.snapshot()->packets.at(e.primary).count, e.event_weight);
        if (!m::PicArena::positive_carriers(remaining).empty()) e.untouched_primary_id = id + "/unscattered-primary";
        if constexpr (std::is_same_v<Arena, m::CollisionScope>) e.phase = scoped.phase();
        else e.phase = interval.phase;
        e.phase.state_prefix = prefix;
        e.provenance = {m::ArenaProvenance::Kind::ConditionalPhysicalData, arena_.dataset_sha256,
          arena_.kernel_manifest_sha256, arena_.kernel_binding_file_sha256, arena_.model,
          "ssa_weighted_pair_clock_v1", prefix, proposal.event_subflow_elapsed_s()};
        for (int i = 0; i < 3; ++i) e.products.push_back({i == 0 ? e.primary : id + (i == 1 ? "/secondary" : "/ion"),
          products.particles[i].state, products.particles[i].proper_velocity_m_s,
          i == 0 ? m::ArenaProduct::At::Primary : m::ArenaProduct::At::Target});
        auto reaction = scoped.prepare(e);
        pc::require(clock->can_publish(proposal) && scoped.can_publish(reaction), "reaction common preflight");
        pc::require(scoped.try_publish(reaction), "native reaction publication");
        pc::require(clock->try_publish(proposal), "reaction RNG publication");
        ledger.close(segment, elapsed, proposal.event_subflow_elapsed_s());
        ledger.executed(c.energy_eV, c.event_weight);
        elapsed = proposal.event_subflow_elapsed_s();
        Json record{{"event", id}, {"cell", c.cell}, {"physical_volume_m3", c.volume_m3},
          {"weight", c.event_weight}, {"rate_s", c.rate_s}, {"total_rate_s", proposal.total_rate_s()},
          {"energy_eV", c.energy_eV}, {"physical_stage_s", stage.frozen_physical_stage_s},
          {"collision_elapsed_s", elapsed}, {"native", reaction.audit()}};
        if (!extra.is_null()) for (auto& [k, v] : extra.items()) record[k] = v;
        events.push_back(std::move(record));
      }
      next->scheduler = clock->checkpoint();next->product_rng = encode(rng);++next->next_interval;
      // Throws CoverageRefusal before any completion is formed; the native
      // compound/private frame and both RNG streams are then discarded.
      auto exposure_ledger = ledger.finish();
      next->coverage = scope_.advance(before.coverage, exposure_ledger);
      const auto& recorded = slot == Slot::before_transport ? *next->coverage.pending : next->coverage.last_step.back();
      coverage_audit = scope_.exposure_json(recorded);
      m::ArenaStochasticCompletion completion;
      completion.integration_owner_id = owner_id;
      completion.before_rng_clock_sha256 = before_hash;
      completion.after_rng_clock_sha256 = state_hash(*next);
      completion.channel_id = w::parent_channel;
      completion.dataset_sha256 = arena_.dataset_sha256;
      completion.kernel_manifest_sha256 = arena_.kernel_manifest_sha256;
      completion.kernel_binding_file_sha256 = arena_.kernel_binding_file_sha256;
      completion.model = arena_.model;
      completion.common_clock_sha256 = interval.common_clock_sha256;
      completion.stage_sequence = interval.sequence;
      completion.frozen_physical_stage_s = interval.phase.physical_stage_s;
      completion.common_interval_begin_s = interval.begin_s;
      completion.common_interval_end_s = interval.end_s;
      return completion;
    });
    Built out;
    out.arena = std::make_unique<m::PreparedArena>(std::move(token));out.after = std::move(next);
    out.audit = {{"schema", "native-parent-subflow-coverage-receipt-v2"}, {"conditional_parent_only", true},
      {"full_mechanism", false}, {"physical_time_advanced", false}, {"interval_sequence", interval.sequence},
      {"interval_begin_s", interval.begin_s}, {"interval_end_s", interval.end_s},
      {"events", std::move(events)}, {"pairs", pairs}, {"queries", queries}, {"terminal_clock_sampled", terminal},
      {"coverage_scope", scope_.identity()}, {"coverage", std::move(coverage_audit)},
      {"compound", out.arena->audit()}};
    return out;
  }
  State restored_state(const Json& j) const {
    pc::require(j.at("schema") == (full_time_ ? schema_full_time : schema_legacy) && j.at("grid") == grid_.identity &&
      j.at("dataset") == binding_.dataset_sha256 && j.at("kernel") == binding_.kernel_sha256 &&
      j.at("binding_file") == binding_.kernel_binding_file_sha256 && j.at("arena") == arena_ref_.state_prefix() &&
      j.at("integration_configuration") == integration_config_ && j.at("coverage_scope") == scope_.identity(),
      "subflow restart binding mismatch");
    State s;s.scheduler = j.at("scheduler");s.product_rng = j.at("product_rng");
    pc::require(s.scheduler.size() < 10000, "scheduler restart size cap");
    (void)pc::Scheduler::restore(s.scheduler);(void)decode(s.product_rng);
    pc::require(j.at("next_event").is_number_unsigned() && j.at("next_interval").is_number_unsigned(), "restart sequence type");
    s.next_event = j.at("next_event");s.next_interval = j.at("next_interval");
    pc::require(s.next_event > 0 && s.next_event != UINT64_MAX && s.next_interval != UINT64_MAX, "restart event/interval sequence exhausted");
    s.coverage = scope_.ledger_from(j.at("coverage"));
    pc::require(!s.coverage.pending, "restart cannot carry an unpublished exposure ledger");
    return s;
  }

public:
  ParentSubflow(m::PicArena& arena, std::array<std::uint32_t, 3> grid_counts, Binding binding, Limits limits,
                const w::ParentBridge& bridge, const CoverageScope& scope, std::uint64_t clock_seed,
                std::uint64_t product_seed, bool full_time = false)
      : arena_ref_(arena), grid_(arena.native_options(), grid_counts, limits.cells), binding_(std::move(binding)), limits_(limits),
        bridge_(bridge), scope_(scope), arena_(scope_.arena_binding(binding_.kernel_binding_file_sha256)), full_time_(full_time),
        accepted_(std::make_unique<State>()) {
    for (const auto& s : {binding_.dataset_sha256, binding_.kernel_sha256, binding_.kernel_binding_file_sha256})
      pc::require(pc::digest(s), "unbound native subflow inputs");
    const auto c = bridge_.channel();
    pc::require(c.data_sha256 == binding_.dataset_sha256 && c.kernel_sha256 == binding_.kernel_sha256, "bridge/binding mismatch");
    pc::require(scope_.dataset_sha256() == binding_.dataset_sha256 && scope_.kernel_sha256() == binding_.kernel_sha256,
      "coverage scope/binding mismatch");
    scope_.require_bridge(bridge_);
    const auto& catalog = arena_ref_.channel_binding(w::parent_channel);
    pc::require(catalog.dataset_sha256 == arena_.dataset_sha256 && catalog.kernel_manifest_sha256 == arena_.kernel_manifest_sha256 &&
      catalog.kernel_binding_file_sha256 == arena_.kernel_binding_file_sha256 && catalog.model == arena_.model,
      "native registered channel differs from subflow binding");
    pc::require(limits_.reactions > 0 && limits_.carriers > 0 && limits_.cells >= grid_.cells.size(), "native subflow limits");
    integration_config_ = configuration_identity(grid_, binding_, limits_, scope_);
    for (const auto& cell : grid_.cells) volumes_.emplace(cell.id, cell.volume_m3);
    accepted_->scheduler = pc::Scheduler(clock_seed).checkpoint();
    accepted_->product_rng = encode(std::mt19937_64(product_seed));
  }
  static std::string initial_full_time_rng_sha256(const m::GapPicOptions& options,
      std::array<std::uint32_t, 3> grid_counts, const Binding& binding, const Limits& limits,
      const CoverageScope& scope, std::uint64_t clock_seed, std::uint64_t product_seed) {
    const Grid grid(options, grid_counts, limits.cells);State state;
    state.scheduler = pc::Scheduler(clock_seed).checkpoint();state.product_rng = encode(std::mt19937_64(product_seed));
    const Json payload{{"scheduler", state.scheduler}, {"product_rng", state.product_rng}, {"next_event", state.next_event},
                       {"next_interval", state.next_interval}, {"coverage", scope.ledger_json(state.coverage)}};
    return hash_with(configuration_identity(grid, binding, limits, scope), payload, true);
  }
  static std::string full_time_configuration_sha256(const m::GapPicOptions& options,
      std::array<std::uint32_t, 3> grid_counts, const Binding& binding, const Limits& limits, const CoverageScope& scope) {
    return configuration_identity(Grid(options, grid_counts, limits.cells), binding, limits, scope);
  }
  ParentSubflow(const ParentSubflow&) = delete;
  ParentSubflow& operator=(const ParentSubflow&) = delete;
  Json checkpoint() const {
    check_completion(*accepted_);
    return {{"schema", full_time_ ? schema_full_time : schema_legacy}, {"grid", grid_.identity},
      {"dataset", binding_.dataset_sha256}, {"kernel", binding_.kernel_sha256},
      {"binding_file", binding_.kernel_binding_file_sha256}, {"arena", arena_ref_.state_prefix()},
      {"scheduler", accepted_->scheduler}, {"product_rng", accepted_->product_rng},
      {"next_event", accepted_->next_event}, {"next_interval", accepted_->next_interval},
      {"integration_configuration", integration_config_}, {"coverage_scope", scope_.identity()},
      {"coverage", scope_.ledger_json(accepted_->coverage)}};
  }
  // Called only beside a matching restored native arena. Validation precedes
  // swap; the ledger must hash to the native accepted completion.
  void restore(const Json& j) {
    pc::require(generation_ != UINT64_MAX, "subflow revision exhausted");
    auto s = std::make_unique<State>(restored_state(j));
    check_completion(*s);accepted_.swap(s);++generation_;
  }
  Prepared prepare(const Interval& interval) {
    check_completion(*accepted_);
    pc::require(generation_ != UINT64_MAX && interval.sequence == accepted_->next_interval && interval.sequence != UINT64_MAX, "duplicate/stale subflow interval");
    pc::require(pc::digest(interval.common_clock_sha256) && interval.phase.mesh_prefix == grid_.identity &&
      interval.phase.clock_prefix == interval.common_clock_sha256 && interval.phase.state_prefix == arena_ref_.state_prefix(), "native subflow common stage binding");
    pc::require(!full_time_, "legacy preparation is not a full-time endpoint stage");
    auto built = build_events(arena_ref_, interval, *accepted_);
    Prepared out;out.lifetime = lifetime_;out.generation = generation_;
    out.arena = std::move(built.arena);out.after = std::move(built.after);out.audit_ = std::move(built.audit);
    return out;
  }
  bool can_publish(const Prepared& p) const noexcept {
    return p.lifetime == lifetime_ && p.generation == generation_ && generation_ != UINT64_MAX && p.after && p.arena && arena_ref_.can_publish(*p.arena);
  }
  bool try_publish(Prepared& p) noexcept {
    if (!can_publish(p)) return false;
    if (!arena_ref_.try_publish(*p.arena)) return false;
    accepted_.swap(p.after);++generation_;return true;
  }
  PreparedParentStage prepare_on(m::CollisionScope& scope, const PreparedParentStage* predecessor = nullptr) {
    pc::require(full_time_ && generation_ != UINT64_MAX, "not a full-time parent owner or exhausted revision");
    check_completion(*accepted_);
    const auto& exposure = scope.exposure();const auto& common = exposure.common;
    pc::require(common.accepted_arena_sha256 == arena_ref_.state_prefix() &&
      common.accepted_arena_revision == arena_ref_.snapshot()->revision &&
      common.accepted_step == arena_ref_.snapshot()->accepted_steps, "full-time exposure has stale accepted arena");
    pc::require(pc::digest(common.common_clock_sha256) && common.dt_s > 0 &&
      std::isfinite(common.dt_s) && common.end_s == common.begin_s + common.dt_s && common.end_s > common.begin_s &&
      exposure.collision_duration_s == .5 * common.dt_s && exposure.collision_duration_s > 0,
      "full-time common step or exact half-exposure differs");
    const bool first = exposure.slot == m::FullTimeSlot::before_transport;
    pc::require(exposure.frozen_physical_stage_s == (first ? common.begin_s : common.end_s), "full-time endpoint stage differs");
    const State* before = accepted_.get();
    if (first) {
      pc::require(predecessor == nullptr, "first exposure cannot use a predecessor stage");
      pc::require(before->next_interval <= UINT64_MAX - 3,
        "common step lacks restorable counter capacity for both collision exposures");
    } else {
      pc::require(predecessor && predecessor->lifetime_ == lifetime_ && predecessor->generation_ == generation_ &&
        predecessor->slot_ == m::FullTimeSlot::before_transport && predecessor->common_binding_ == common_binding(common) &&
        predecessor->accepted_arena_ == common.accepted_arena_sha256 && predecessor->after_,
        "second exposure requires matching first staged random streams");
      before = predecessor->after_.get();
    }
    check_completion_on(scope, *before);
    pc::require(exposure.owner_stage_sequence == before->next_interval && before->next_interval != UINT64_MAX,
      "full-time exposure sequence differs");
    Interval request;request.sequence = before->next_interval;request.begin_s = common.begin_s;request.end_s = common.end_s;
    request.phase = scope.phase();request.common_clock_sha256 = common.common_clock_sha256;
    pc::require(request.phase.state_prefix == scope.state_prefix() && request.phase.clock_prefix == common.common_clock_sha256 &&
      request.phase.physical_stage_s == exposure.frozen_physical_stage_s, "full-time scope phase mismatch");
    auto built = build_events(scope, request, *before, &exposure);
    PreparedParentStage stage;stage.lifetime_ = lifetime_;stage.generation_ = generation_;stage.slot_ = exposure.slot;
    stage.exposure_key_ = scope.exposure_key();
    pc::require(bool(stage.exposure_key_), "full-time exposure has no native identity key");
    stage.exposure_binding_ = exposure_binding(exposure);
    stage.common_binding_ = common_binding(common);stage.common_clock_ = common.common_clock_sha256;stage.accepted_arena_ = common.accepted_arena_sha256;
    stage.before_hash_ = state_hash(*before);stage.after_hash_ = state_hash(*built.after);
    stage.after_ = std::move(built.after);stage.audit_ = std::move(built.audit);
    stage.audit_["schema"] = "native-parent-full-time-exposure-coverage-v2";
    stage.audit_["slot"] = first ? "before_transport" : "after_transport";
    stage.audit_["authoritative_common_dt_s"] = common.dt_s;
    stage.audit_["collision_exposure_s"] = exposure.collision_duration_s;
    stage.audit_["integration_configuration"] = integration_config_;
    stage.private_prefix_ = built.arena->candidate_state_prefix();
    pc::require(scope.can_publish(*built.arena), "private exposure publication preflight");
    pc::require(scope.try_publish(*built.arena), "private exposure publication");
    return stage;
  }
  std::string staged_rng_sha256(const PreparedParentStage& stage) const {
    pc::require(full_time_ && stage.lifetime_ == lifetime_ && stage.generation_ == generation_ && stage.after_, "foreign or stale staged RNG");
    const auto actual = state_hash(*stage.after_);
    pc::require(actual == stage.after_hash_, "staged RNG payload differs from its sealed identity");
    return actual;
  }
  PreparedParentState seal(const PreparedParentStage& first, const PreparedParentStage& last,
                           const m::PreparedArena& outer) const {
    pc::require(full_time_ && generation_ != UINT64_MAX, "not a full-time RNG owner");
    pc::require(first.lifetime_ == lifetime_ && last.lifetime_ == lifetime_ && first.generation_ == generation_ &&
      last.generation_ == generation_ && first.after_ && last.after_ && first.slot_ == m::FullTimeSlot::before_transport &&
      last.slot_ == m::FullTimeSlot::after_transport && first.common_binding_ == last.common_binding_ && first.common_clock_ == last.common_clock_ &&
      first.accepted_arena_ == arena_ref_.state_prefix() && last.accepted_arena_ == first.accepted_arena_ &&
      first.before_hash_ == state_hash(*accepted_) && first.after_hash_ == state_hash(*first.after_) &&
      last.after_hash_ == state_hash(*last.after_) && last.before_hash_ == first.after_hash_ &&
      first.after_->next_interval != UINT64_MAX && last.after_->next_interval == first.after_->next_interval + 1,
      "first/second staged RNG chain is incomplete or foreign");
    pc::require(first.after_->coverage.pending && !last.after_->coverage.pending && last.after_->coverage.last_step.size() == 2,
      "staged coverage ledgers do not close one two-exposure step");
    const auto native_exposures = outer.candidate_full_time_exposures();
    pc::require(native_exposures && pc::digest(native_exposures->outer_step_sha256) &&
      common_binding(native_exposures->common) == first.common_binding_,
      "outer native token has no matching full-time exposure history");
    const auto matches = [&](const PreparedParentStage& stage, const m::FullTimeExposureReceipt& receipt) {
      const auto& c = receipt.completion;
      return stage.exposure_key_ && stage.exposure_key_ == receipt.key &&
        stage.private_prefix_ == receipt.completed_state_prefix &&
        stage.exposure_binding_ == exposure_binding(receipt.exposure) &&
        pc::digest(receipt.completion_sha256) && c.integration_owner_id == owner_id &&
        c.before_rng_clock_sha256 == stage.before_hash_ && c.after_rng_clock_sha256 == stage.after_hash_ &&
        c.stage_sequence != UINT64_MAX && c.stage_sequence + 1 == stage.after_->next_interval;
    };
    pc::require(matches(first, native_exposures->first) && matches(last, native_exposures->second),
      "staged random streams belong to different private native exposures");
    const auto completion = outer.candidate_stochastic_completion(owner_id);
    pc::require(completion && completion->after_rng_clock_sha256 == last.after_hash_ &&
      completion->before_rng_clock_sha256 == last.before_hash_ && completion->stage_sequence != UINT64_MAX && completion->stage_sequence + 1 == last.after_->next_interval &&
      completion->common_clock_sha256 == last.common_clock_ && completion->channel_id == w::parent_channel && completion->model == arena_.model &&
      completion->dataset_sha256 == arena_.dataset_sha256 && completion->kernel_manifest_sha256 == arena_.kernel_manifest_sha256 &&
      completion->kernel_binding_file_sha256 == arena_.kernel_binding_file_sha256,
      "outer native token and final staged RNG completion differ");
    pc::require(arena_ref_.can_publish(outer), "outer native token is not publishable");
    PreparedParentState state;state.lifetime_ = lifetime_;state.generation_ = generation_;
    state.after_ = std::make_unique<State>(*last.after_);state.after_hash_ = last.after_hash_;
    state.native_prefix_ = outer.candidate_state_prefix();state.common_clock_ = completion->common_clock_sha256;
    state.coverage_ = scope_.ledger_json(state.after_->coverage);
    return state;
  }
  bool can_publish(const PreparedParentState& state, const m::PreparedArena& outer) const noexcept {
    try {
      if (!full_time_ || state.lifetime_ != lifetime_ || state.generation_ != generation_ || generation_ == UINT64_MAX ||
         !state.after_ || !arena_ref_.can_publish(outer) || state.native_prefix_ != outer.candidate_state_prefix()) return false;
      const auto completion = outer.candidate_stochastic_completion(owner_id);
      return completion && state.after_hash_ == state_hash(*state.after_) &&
        completion->after_rng_clock_sha256 == state.after_hash_ &&
        completion->common_clock_sha256 == state.common_clock_ && completion->stage_sequence != UINT64_MAX &&
        completion->stage_sequence + 1 == state.after_->next_interval;
    } catch (...) { return false; }
  }
  std::string accepted_rng_sha256() const {
    pc::require(full_time_, "accepted full-time RNG observation requires full-time mode");
    check_completion(*accepted_);return state_hash(*accepted_);
  }
  std::string candidate_rng_sha256(const PreparedParentState& state, const m::PreparedArena& outer) const {
    pc::require(can_publish(state, outer), "private RNG candidate is stale or differs from native completion");
    return state_hash(*state.after_);
  }
  const std::string& configuration_sha256() const noexcept { return integration_config_; }
  void publish_prevalidated(PreparedParentState& state) noexcept {
    if (state.lifetime_ != lifetime_ || state.generation_ != generation_ || generation_ == UINT64_MAX || !state.after_) std::terminate();
    accepted_.swap(state.after_);++generation_;
  }
  const std::string& grid_identity() const { return grid_.identity; }
  const CoverageScope& coverage_scope() const noexcept { return scope_; }
  // Accepted unresolved-coverage ledger (last accepted step and cumulative).
  Json accepted_coverage() const { check_completion(*accepted_);return scope_.ledger_json(accepted_->coverage); }
  static constexpr const char* integration_owner() noexcept { return owner_id; }
};
} // namespace raikou::native_water::coverage_v2
