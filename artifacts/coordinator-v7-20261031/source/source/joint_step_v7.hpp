#pragma once
#include "identity_v7.hpp"
#include "common_commit_v7.hpp"

// Coordinator v7: the coordinator v6.1 step (pic-water-common-commit-v6.1-20261031)
// with a variable step length. Changes relative to v6.1 (README.txt):
//   * prepare_v7 admits 0 < dt <= the configured cap of the run (OwnersV7::dt_cap_s,
//     itself <= the largest step the coupled owners admit, owner_dt_max)
//     in place of the literal dt <= 1e-14, and binds the step controller's decision
//     (dt bits, decision digest) into the trial;
//   * can_publish_v7 applies common-commit predicate v7 (common_commit_v7.hpp: the
//     V3/native-v6 conditions with dt <= cap and dt equal to the decision).
// Every owner operation, identity text and refusal point is the v6.1 one, so a run
// at dt = 1e-14 reproduces v6.1 bit for bit.

// Coordinator v6.1: the coordinator v6 step (pic-water-common-commit-v6-20261031,
// sealed) with heat/source owner v6.1 and less host work per step; accepted states,
// ledgers and digests are those of v6. Changes (README.txt):
//   * source field: Field from the immutable mesh identity of the Hall mesh (hashed once,
//     reused while the live mesh is exactly equal to it); the one-sided extraction field
//     is sampled from that Field (v6 built a second trial field from a mesh copy);
//   * emission exposures: source v6.1 prepare_patch (exact fast spectrum); the prepared
//     flux of each exposure is reused for the electron births of that stream (v6
//     re-prepared the patch on a third trial field);
//   * observe(): digests of immutable objects (Hall model motion mesh and frame) are
//     computed once per model object and options;
//   * sub-phase timings of the source field, emission exposures and can_publish.
// v6 description:
// one accepted step with surface births, first-hit sinks and the
// one thermal owner.
//   first exposure (parent coverage v2)
//   -> source v6 epochs, desorption quotes and emission exposures on the
//      accepted start-of-step Hall field -> arrival v1.1 plan and staging at the
//      native birth base B_s -> source v6 trials (quantized neutral leaves,
//      emission surface loads) -> native v6 birth credit and fork
//   -> Hall v3 closed proposal (newborns, source loads, first-hit sinks, banks,
//      ThermalV3) with the selected motion provider and field backend
//   -> native admission of the binary flows motion record -> second exposure
//   -> seal -> preflight of every owner -> allocation-free publication
//      native -> chemistry -> Hall -> source owners -> coordinator source seal,
//      then the arrival ledger commit (its checks are replayed in preflight).
// Ions are ordinary moving carriers of the Hall motion; screening can only come
// from the Hall EM solve. Test instrumentation is reachable only through
// TrialTestAccessV7, declared here and defined only in test code.
namespace raikou::native_water::v7 {
namespace c7 = ::raikou::native_water::common_commit_v7;
// Admitted maximum step of the field owner: Model::dt_max() of a Hall owner with a
// configurable cap (v3.1), else the v3 literal 1e-14 s (bindings_v7.hpp).
inline double field_owner_dt_max(const hc::Model& m) { return bind::field_owner_dt_max(m); }
// Largest step every coupled owner admits (bindings_v7.hpp).
inline double owner_dt_max(const hc::Model& m) { return bind::coupled_dt_max(m); }
// Decision record digest of a step whose dt is given by the caller (no controller).
inline std::string direct_decision_sha256(double dt_s, double cap_s) {
  Encoding e;
  e.text("coordinator-v7-direct-dt-decision-v1");
  e.f64(dt_s); e.f64(cap_s);
  return e.finish();
}
struct TrialTestAccessV7;
struct GenesisTestAccessV7;

struct InputsV7 {
  std::vector<m::ArenaStateSpec> states;
  std::vector<m::ArenaChannel> channels;
  std::string species_sha256, frame_sha256, source_configuration_sha256, settings_text;
  std::size_t maximum_native_records = 4096, maximum_carriers = 100000;
};

// Configuration identity of the coupled step.
inline std::string configuration_sha256(const std::string& chemistry, const std::string& field, const InputsV7& in,
                                        const std::string& arrival_policy) {
  need(cc::digest(chemistry) && cc::digest(field) && cc::digest(in.species_sha256) && cc::digest(in.frame_sha256) &&
           cc::digest(in.source_configuration_sha256) && cc::digest(arrival_policy) && in.maximum_native_records > 0 && in.maximum_carriers > 0,
       "explicit identities and finite resource caps required");
  need(species_sha256(in.states, in.channels) == in.species_sha256, "declared species identity differs from its catalog");
  Encoding e;
  e.text("coordinator-v6-configuration-v1");
  e.text(field); e.text(chemistry); e.text(in.species_sha256); e.text(in.frame_sha256);
  e.u64(in.maximum_native_records); e.u64(in.maximum_carriers);
  e.text(in.source_configuration_sha256); e.text(arrival_policy); e.text(in.settings_text);
  e.text("surface births: source v6 quantized desorption (arrival v1.1) and face-born thermal-field emission; first-hit sinks; "
         "one thermal owner (Hall ThermalV3); source exposures on the accepted start-of-step Hall field");
  return e.finish();
}
inline std::string source_configuration_sha256(const std::vector<hb::SliceSource>& slices, const std::vector<std::string>& owner_signatures) {
  need(slices.size() == owner_signatures.size(), "one source owner per slice");
  Encoding e;
  e.text("coordinator-v6-source-configuration-v1");
  e.u64(slices.size());
  for (std::size_t i = 0; i < slices.size(); ++i) { e.i64(slices[i].body); e.i64(slices[i].slice); e.text(owner_signatures[i]); }
  return e.finish();
}
inline std::string genesis_source_seal(const std::string& configuration, const std::vector<std::string>& owner_states,
                                       const std::string& arrival_state, const std::string& arrival_policy) {
  Encoding e;
  e.text("coordinator-v6-genesis-source-v1");
  e.text(configuration);
  e.u64(owner_states.size());
  for (const auto& s : owner_states) e.text(s);
  e.text(arrival_state); e.text(arrival_policy);
  return e.finish();
}

// Owners of one run. Member order fixes destruction order: arrival, source owners
// (const views into the Hall thermal owner), chemistry (references the arena),
// Hall owner, arena, motion provider, model.
struct OwnersV7 {
  std::shared_ptr<const hc::Model> model;
  std::shared_ptr<hc::MotionProvider> motion = std::make_shared<hc::CpuMotionProvider>();
  std::unique_ptr<m::PicArena> arena;
  std::unique_ptr<hc::Owner> hall;
  std::unique_ptr<cv::ParentSubflow> chemistry;
  std::vector<hb::SliceSource> slices;
  std::vector<PatchRef> patches;
  std::vector<std::unique_ptr<s6::Owner>> sources;
  std::shared_ptr<const pa::Policy> arrival_policy;
  std::unique_ptr<pa::Owner> arrival;
  InputsV7 inputs;
  PilotSettings settings;
  // Accepted source seal (equal to the native boundary source seal) and the live
  // owner digests recorded at the last publication.
  std::string accepted_source_seal;
  std::vector<std::string> accepted_source_states;
  // Immutable identity of the Hall mesh for source-side trial fields (hashed once).
  std::shared_ptr<const s6::MeshIdentity> mesh_identity;
  // Digests of immutable objects: the Hall model's motion mesh digest and the
  // native/Hall frame digest, for one (model object, native options) pair.
  struct ImmutableDigests { const hc::Model* model = nullptr; const m::GapPicOptions* options = nullptr; std::string options_bytes, mesh, frame; };
  mutable ImmutableDigests immutable;
  std::string accepted_arrival_state;
  Json genesis;
  // Configured step cap of the run (v7): 0 < dt <= dt_cap_s <= owner_dt_max(model).
  double dt_cap_s = bind::literal_dt_max_s;
};

inline std::string arrival_configuration(const OwnersV7& o) { return o.arrival_policy->configuration_sha256(); }

// Live source seal: the coordinator-held accepted seal, valid only while every
// live source owner and the arrival owner still have the recorded digests.
inline std::string live_source_seal(const OwnersV7& o) {
  need(o.sources.size() == o.accepted_source_states.size(), "source owner set changed");
  for (std::size_t i = 0; i < o.sources.size(); ++i)
    need(o.sources[i]->state_sha256() == o.accepted_source_states[i], "live source owner differs from its accepted state");
  need(o.arrival->state_sha256() == o.accepted_arrival_state, "live arrival owner differs from its accepted state");
  return o.accepted_source_seal;
}

// Byte image of every native options field: a different or changed options object
// recomputes the cached digests.
inline std::string options_bytes(const m::GapPicOptions& g) {
  Encoding e;
  e.i64(g.nx); e.i64(g.ny);
  for (double x : {g.width_m, g.gap_m, g.depth_m, g.lower_potential_V, g.upper_potential_V}) e.f64(x);
  for (double x : g.magnetic_field_T) e.f64(x);
  e.u64(g.distributed_fields); e.u64(g.finite_foils);
  for (const auto* f : {&g.lower_foil, &g.upper_foil}) for (double x : {f->left_m, f->right_m, f->bottom_m, f->top_m}) e.f64(x);
  e.u64(g.finite_axial_extent); e.f64(g.short_z_m); e.f64(g.feed_z_m);
  return e.bytes;
}
// The Hall model is immutable (shared const object); its motion mesh digest and the
// frame digest are recomputed only for another model object or changed native options.
inline const OwnersV7::ImmutableDigests& immutable_digests(const OwnersV7& o) {
  const auto& model = o.hall->model();
  const auto& options = o.arena->native_options();
  auto& c = o.immutable;
  if (c.model != &model || c.options != &options || c.options_bytes != options_bytes(options)) {
    c.mesh = model.motion_mesh_sha256(options);
    c.frame = frame_sha256(options, model);
    c.model = &model; c.options = &options; c.options_bytes = options_bytes(options);
  }
  return c;
}
inline cc::ObservedOwners observe(const OwnersV7& o) {
  const auto boundary = o.arena->full_time_boundary();
  need(bool(boundary), "actual full-time native boundary required");
  cc::ObservedOwners out;
  out.native = *boundary;
  out.actual_native_prefix = o.arena->state_prefix();
  const auto snapshot = o.arena->snapshot();
  out.actual_native_revision = snapshot->revision;
  out.actual_native_step = snapshot->accepted_steps;
  out.field = o.hall->field_sha256();
  out.heat = o.hall->heat_sha256();
  out.rng = o.chemistry->accepted_rng_sha256();
  const auto& id = immutable_digests(o);
  out.mesh = id.mesh;
  need(out.mesh == o.hall->motion_mesh_sha256(), "Hall owner geometry differs from the live native geometry");
  out.material = o.hall->model().material_sha256();
  require_arena_species(*o.arena, o.inputs.states, o.inputs.channels, o.inputs.maximum_carriers);
  out.frame = id.frame;
  need(out.frame == o.inputs.frame_sha256, "declared frame differs from the live native/Hall translation frame");
  std::vector<std::string> signatures;
  for (const auto& s : o.sources) signatures.push_back(s->configuration_sha256());
  need(source_configuration_sha256(o.slices, signatures) == o.inputs.source_configuration_sha256, "live source configuration differs");
  out.configuration = configuration_sha256(o.chemistry->configuration_sha256(), o.hall->configuration_sha256(), o.inputs, arrival_configuration(o));
  out.source = live_source_seal(o);
  out.field_time_s = out.heat_time_s = o.hall->time_s();
  out.field_revision = out.heat_revision = o.hall->revision();
  return out;
}

// v6.1 sub-phase timings (written to subphase-timings.bin by the driver).
struct SubTimesV61 {
  double start_field_digest_s = 0, source_field_build_s = 0;
  double emission_extracting_s = 0, emission_retarding_s = 0, emission_other_s = 0;
  unsigned extracting_patches = 0, retarding_patches = 0, unresolved_patches = 0;
  double cp_native_s = 0, cp_chemistry_s = 0, cp_observe_s = 0, cp_predicate_s = 0, cp_sources_s = 0, cp_arrival_s = 0, cp_hall_s = 0;
  unsigned cp_source_preflights = 0;
};
struct PhaseTimesV7 {
  double observe_s = 0, begin_s = 0, first_exposure_s = 0, source_field_s = 0, source_quotes_s = 0, emission_exposure_s = 0,
         arrival_plan_s = 0, emission_births_s = 0, source_trials_s = 0, birth_credit_s = 0, census_s = 0, hall_proposal_s = 0,
         binding_checks_s = 0, motion_admission_s = 0, second_exposure_s = 0, seal_s = 0, ledger_s = 0;
  unsigned hall_iterations = 0;
  hc::ProposalTimings hall;
};
// Numbers of one prepared step for the per-step ledger (filled before publication).
struct StepNumbersV7 {
  std::uint64_t electron_births = 0, neutral_births = 0, sinks = 0, census = 0, newborns = 0;
  double expected_electrons = 0, expected_neutral_molecules = 0;
  double emission_nottingham_J = 0, emission_reference_J = 0, emission_matching_work_J = 0, emission_kinetic_J = 0;
  double neutral_kinetic_J = 0, neutral_internal_J = 0, neutral_adsorbed_reference_J = 0, source_load_J = 0;
  double max_extraction_field_V_m = -INFINITY, min_extraction_field_V_m = INFINITY, max_log10_expected_electrons = -INFINITY;
  double max_patch_T_K = 0, min_patch_T_K = INFINITY;
  std::uint64_t unresolved_emission_patches = 0;
  double unresolved_emission_lower_electrons = 0, unresolved_emission_upper_electrons = 0;
  m::ArenaCount source_debit;          // exact molecules released this step
  m::ArenaCount sunk_h2o, sunk_h2o_ion;  // exact weights removed by sinks this step
  std::uint64_t parent_events = 0;
};
inline constexpr int prepare_fault_points_v7 = 12;

class TrialV7 {
  friend TrialV7 prepare_v7(OwnersV7&, double, int, const std::string*);
  friend bool can_publish_v7(const OwnersV7&, const TrialV7&) noexcept;
  friend bool try_publish_v7(OwnersV7&, TrialV7&) noexcept;
  friend struct TrialTestAccessV7;
  std::unique_ptr<m::PrivateFullTimeStep> outer_;
  std::optional<m::BirthScope> scope_;
  std::optional<m::PreparedBirths> births_;
  std::unique_ptr<m::PreparedArena> native_;
  cv::PreparedParentState random_;
  std::unique_ptr<hc::JointProposal> field_heat_;
  std::vector<m::ArenaRateEntry> census_;
  std::vector<hc::Newborn> newborns_;
  std::vector<hc::SurfaceLoad> loads_;
  std::vector<std::optional<s6::Trial>> source_trials_;
  std::shared_ptr<const s6::Field> source_field_;
  v3k::TrialFrame source_frame_;
  std::optional<pa::Staged> staged_;
  pa::NativeBirthReceipt arrival_receipt_;
  std::string source_after_;
  InputsV7 inputs_;
  cc::ExpectedStep expected_;
  c7::StepClock clock_;
  StepNumbersV7 numbers_;
  PhaseTimesV7 times_{};
  mutable SubTimesV61 sub_{};
  std::vector<std::optional<s6::PreparedPatch>> prepared_;   // emission exposure per patch (index of OwnersV7::patches)
  bool consumed_ = false;
  // Test instrumentation (defaults leave behavior unchanged).
  int can_publish_fault_ = 0;
  mutable int refusal_point_ = 0;
  // v7 test instrumentation: evaluate the native v6 (V3) predicate next to predicate v7.
  bool compare_v6_predicate_ = false;
  mutable int v6_predicate_ = -1, v7_predicate_ = -1;
  void (*publish_probe_)(int) = nullptr;
  TrialV7() = default;

 public:
  TrialV7(TrialV7&&) noexcept = default;
  TrialV7& operator=(TrialV7&&) noexcept = default;
  TrialV7(const TrialV7&) = delete;
  TrialV7& operator=(const TrialV7&) = delete;
  const PhaseTimesV7& times() const noexcept { return times_; }
  const SubTimesV61& sub_times() const noexcept { return sub_; }
  const StepNumbersV7& numbers() const noexcept { return numbers_; }
  const hc::JointProposal& proposal() const { need(bool(field_heat_), "trial has no proposal"); return *field_heat_; }
  const cv::PreparedParentState& chemistry_state() const noexcept { return random_; }
  std::shared_ptr<const m::FullTimeStepFlows> flows() const { return native_ ? native_->candidate_full_time_flows() : nullptr; }
  const std::string& source_after() const noexcept { return source_after_; }
  const std::vector<m::ArenaRateEntry>& census() const noexcept { return census_; }
  int refusal_point() const noexcept { return refusal_point_; }
  bool consumed() const noexcept { return consumed_; }
  const c7::StepClock& clock() const noexcept { return clock_; }
  const std::string& common_clock_sha256() const noexcept { return expected_.common_clock; }
  double dt_s() const noexcept { return expected_.dt_s; }
};

// Arrival commit checks of pf_arrival_v11::Owner::commit, replayed before any
// publication (the commit itself allocates its ledger record and therefore runs
// after the allocation-free swap window).
inline bool arrival_commit_admissible(const pa::Owner& owner, const pa::Staged& s, const pa::NativeBirthReceipt& r) noexcept {
  try {
    const auto& plan = s.plan();
    if (s.owner_base_sha256() != owner.state_sha256() || plan.step() != owner.next_step()) return false;

    const std::uint64_t n = plan.births().size();
    if (r.base != s.native_base() || r.next_id_after != s.native_base() + n || r.entries.size() != n) return false;
    for (std::uint64_t k = 0; k < n; ++k) {
      const auto& b = plan.births()[k];
      const auto& e = r.entries[k];
      if (!(std::isfinite(b.time_s) && std::isfinite(b.remaining_s))) return false;
      if (e.packet_id != b.packet_id() || e.state != owner.policy().stream(b.stream_id).species_state || e.native_id != b.native_id ||
          std::bit_cast<std::uint64_t>(e.weight) != std::bit_cast<std::uint64_t>(b.weight))
        return false;
    }
    return true;
  } catch (...) {
    return false;
  }
}

inline TrialV7 prepare_v7(OwnersV7& o, double dt_s, int failure_point = 0, const std::string* decision_sha256 = nullptr) {
  using clock = std::chrono::steady_clock;
  need((failure_point >= 0 && failure_point <= prepare_fault_points_v7) || (failure_point > 100 && failure_point <= 100 + prepare_fault_points_v7),
       "unknown manufactured preparation fault point");
  // Points 1..12: the v6 manufactured fault; 101..112 (v7 tests): a manufactured
  // refusal with the text of an owner step-size guard at the same point.
  const auto fail = [&](int point) {
    if (failure_point == point) throw std::runtime_error("manufactured v6 preparation fault");
    if (failure_point == 100 + point) throw std::runtime_error("manufactured step-size refusal: motion full-census plasma guard");
  };
  TrialV7 out;
  auto t = clock::now();
  const auto before = observe(o);
  need(cc::owner_values_match(before), "accepted native/field/heat/RNG/source boundaries differ");
  const double owner_max = owner_dt_max(o.hall->model());
  need(std::isfinite(o.dt_cap_s) && o.dt_cap_s > 0 && o.dt_cap_s <= owner_max, "configured step cap within the field owner maximum");
  need(std::isfinite(dt_s) && dt_s > 0 && dt_s <= o.dt_cap_s && .5 * dt_s > 0, "dt within the configured step cap");
  const auto& boundary = before.native;
  const auto first_sequence = o.chemistry->checkpoint().at("next_interval").get<std::uint64_t>();
  need(first_sequence <= UINT64_MAX - 3, "two restorable collision exposure counters required");
  m::FullTimeStepBinding common;
  common.accepted_step = boundary.accepted_step;
  common.accepted_arena_revision = boundary.accepted_arena_revision;
  common.begin_s = boundary.time_s;
  common.dt_s = dt_s;
  common.end_s = common.begin_s + dt_s;
  need(std::isfinite(common.end_s) && common.end_s > common.begin_s, "common endpoint is not representable");
  common.accepted_arena_sha256 = before.actual_native_prefix;
  common.accepted_field_sha256 = before.field;
  common.mesh_sha256 = before.mesh;
  common.material_sha256 = before.material;
  {
    Encoding e;
    e.text("coordinator-v6-common-clock-v1");
    e.u64(common.accepted_step); e.u64(common.accepted_arena_revision);
    e.f64(common.begin_s); e.f64(common.end_s); e.f64(common.dt_s);
    for (const std::string* x : std::initializer_list<const std::string*>{&before.actual_native_prefix, &before.field, &before.heat, &before.source, &before.rng, &before.mesh,
                                 &before.material, &before.frame, &before.configuration, &o.inputs.species_sha256})
      e.text(*x);
    e.text(o.accepted_arrival_state); e.text(o.arrival->ledger_head());
    common.common_clock_sha256 = e.finish();
  }
  out.inputs_ = o.inputs;
  out.expected_ = {boundary, common.dt_s, common.end_s, common.common_clock_sha256, before.configuration};
  out.clock_ = {o.dt_cap_s, owner_max, dt_s, decision_sha256 ? *decision_sha256 : direct_decision_sha256(dt_s, o.dt_cap_s)};
  out.times_.observe_s = seconds_since(t);
  fail(1);
  t = clock::now();
  out.outer_ = std::make_unique<m::PrivateFullTimeStep>(
      o.arena->begin_full_time_step(common, o.inputs.maximum_carriers, o.inputs.maximum_native_records));
  out.times_.begin_s = seconds_since(t);
  fail(2);
  t = clock::now();
  m::ExposureBinding exposure{common, m::FullTimeSlot::before_transport, common.begin_s, .5 * common.dt_s, first_sequence, before.field};
  auto first_scope = out.outer_->first_exposure(exposure);
  auto first = o.chemistry->prepare_on(first_scope);
  out.times_.first_exposure_s = seconds_since(t);
  fail(3);

  // ---- births: source v6 + arrival v1.1 through the native v6 birth scope
  t = clock::now();
  out.scope_.emplace(out.outer_->birth_scope());
  auto& scope = *out.scope_;
  const std::uint64_t B = scope.base_native_id();
  need(scope.accepted_step() == common.accepted_step && scope.source_before_sha256() == before.source,
       "native birth scope differs from the accepted step/source seal");
  const auto& mesh = o.hall->mesh();
  const auto& options = o.arena->native_options();
  // Source field: the accepted start-of-step Hall arrays (explicit source coupling).
  const auto& accepted = o.hall->field();
  out.source_frame_.role = v3k::FieldRole::TotalMidpoint;
  out.source_frame_.accepted_revision = o.hall->revision();
  out.source_frame_.trial_stage = 1;
  out.source_frame_.mesh_revision = 1;
  out.source_frame_.material_revision = 2;
  out.source_frame_.begin_s = common.begin_s;
  out.source_frame_.end_s = common.end_s;
  out.source_frame_.dt_s = common.dt_s;
  out.source_frame_.field_stage_s = common.begin_s + .5 * common.dt_s;
  auto ts = clock::now();
  const std::string start_field_digest = h::total_field_sha256(accepted.E, accepted.B);
  out.sub_.start_field_digest_s = seconds_since(ts);
  ts = clock::now();
  out.source_frame_.arena_hash = before.actual_native_prefix;
  out.source_frame_.field_hash = "hall-v3-accepted-start-field/" + start_field_digest;
  out.source_frame_.mesh_hash = before.mesh;
  out.source_frame_.material_hash = before.material;
  out.source_frame_.species_hash = o.inputs.species_sha256;
  // Immutable identity of the live Hall mesh: reused while the mesh is exactly equal to it.
  o.mesh_identity = s6::MeshIdentity::of(mesh);
  out.source_field_ = std::make_shared<const s6::Field>(o.mesh_identity, out.source_frame_, accepted.E, accepted.B);
  const s6::FinalFieldView view{out.source_field_->frame(), mesh, accepted.E, accepted.B};
  out.sub_.source_field_build_s = seconds_since(ts);
  out.times_.source_field_s = seconds_since(t);
  t = clock::now();
  const std::size_t nslices = o.sources.size();
  std::vector<eh::Epoch> epochs;
  std::vector<std::vector<s6::ReleaseQuote>> quotes;
  for (std::size_t i = 0; i < nslices; ++i) {
    epochs.push_back(o.sources[i]->make_epoch(dt_s, "coordinator-v6/accepted-start-field/" + start_field_digest));
    quotes.push_back(o.sources[i]->release_quotes(epochs.back()));
  }
  out.times_.source_quotes_s = seconds_since(t);
  t = clock::now();
  pa::StepInput in;
  in.step = common.accepted_step;
  const v3k::Exposure ex{common.begin_s, dt_s, 1, "physical-quadrature"};
  // One-sided extraction field at every patch origin: the same surface sample and
  // normal-field expression as the kernel exposure (v3 prepare_patch), from the source
  // field's own trial field (same mesh values, frame and arrays as v6's second field).
  out.prepared_.assign(o.patches.size(), std::nullopt);
  for (std::size_t gi = 0; gi < o.patches.size(); ++gi) {
    const auto& r = o.patches[gi];
    const auto tp = clock::now();
    auto& src = *o.sources[r.slice];
    const auto& heat = o.hall->thermal().slice(r.body, r.k);
    const double F = se::extraction_field(out.source_field_->surface(patch_origin(mesh, r), heat.normal(r.local)).E, heat.normal(r.local));
    out.numbers_.max_extraction_field_V_m = std::max(out.numbers_.max_extraction_field_V_m, F);
    out.numbers_.min_extraction_field_V_m = std::min(out.numbers_.min_extraction_field_V_m, F);
    const double T = heat.surface_temperatures().at(r.local);
    out.numbers_.max_patch_T_K = std::max(out.numbers_.max_patch_T_K, T);
    out.numbers_.min_patch_T_K = std::min(out.numbers_.min_patch_T_K, T);
    if (F > 0 && F < 1) {
      // The emission kernel's numerical domain is F <= 0 or 1 <= F <= 5e9 V/m. A
      // field in (0, 1) V/m is an unresolved, bounded exposure: no electron
      // decision this step; its expectation lies between the kernel values at
      // F = 0 and F = 1 V/m (monotone in F), both recorded.
      const auto& cfg = o.slices[r.slice].settings[r.local].spectrum;
      const double scale = heat.area(r.local) * ex.duration_s / se::qe;
      out.numbers_.unresolved_emission_patches += 1;
      out.numbers_.unresolved_emission_lower_electrons += se::evaluate(T, 0., cfg).current_A_m2 * scale;
      out.numbers_.unresolved_emission_upper_electrons += se::evaluate(T, 1., cfg).current_A_m2 * scale;
      out.sub_.unresolved_patches += 1;
    } else {
      const auto pp = src.prepare_patch(epochs[r.slice], r.local, patch_binding(heat, mesh, r), out.source_field_, patch_origin(mesh, r), ex);
      out.prepared_[gi] = pp;
      need(pp.temperature_K() == T, "emission exposure temperature differs from the accepted patch mean");
      in.electrons.push_back(pas::electron_input(r.electron_stream, pp, ex, common.end_s));
      out.numbers_.expected_electrons += pp.expected_number();
      if (pp.expected_number() > 0)
        out.numbers_.max_log10_expected_electrons = std::max(out.numbers_.max_log10_expected_electrons, pp.log_expected_number() / std::log(10.));
      (F > 0 ? out.sub_.emission_extracting_s : out.sub_.emission_retarding_s) += seconds_since(tp);
      (F > 0 ? out.sub_.extracting_patches : out.sub_.retarding_patches) += 1;
    }
    const auto& q = quotes[r.slice][r.local];
    pa::NeutralRelease nr;
    nr.stream_id = r.neutral_stream;
    nr.debit = q.released;
    nr.log_debit = q.log_released;
    nr.inventory = src.inventory(r.local);
    nr.begin_s = epochs[r.slice].t;
    nr.end_s = epochs[r.slice].t + epochs[r.slice].dt;
    in.neutrals.push_back(std::move(nr));
    if (!q.released.zero()) out.numbers_.expected_neutral_molecules += q.released.rounded_for_rate();
  }
  out.times_.emission_exposure_s = seconds_since(t);
  out.sub_.emission_other_s = out.times_.emission_exposure_s - out.sub_.emission_extracting_s - out.sub_.emission_retarding_s;
  t = clock::now();
  auto plan = o.arrival->plan(in);
  out.staged_.emplace(o.arrival->stage(std::move(plan), B));
  auto& staged = *out.staged_;
  out.times_.arrival_plan_s = seconds_since(t);
  fail(4);

  // Electron births: face points with the source kernel's emitted phase space.
  t = clock::now();
  std::map<std::uint32_t, const PatchRef*> by_electron, by_neutral;
  for (const auto& r : o.patches) { by_electron[r.electron_stream] = &r; by_neutral[r.neutral_stream] = &r; }
  std::vector<std::vector<eh::Transfer>> other(nslices);
  std::map<std::uint64_t, av::BirthKinematics> electron_kinematics;   // by ordinal
  {
    for (const auto& b : staged.births()) {
      if (b.kind != pa::Kind::electron) continue;
      const auto& r = *by_electron.at(b.stream_id);
      const auto& heat = o.hall->thermal().slice(r.body, r.k);
      // The flux of this stream's emission exposure (the same thermal epoch, binding,
      // origin, exposure, spectrum configuration and field values as v6's re-preparation).
      const auto gi = std::size_t(&r - o.patches.data());
      need(gi < out.prepared_.size() && out.prepared_[gi].has_value(), "electron birth without its emission exposure");
      const auto& flux = out.prepared_[gi]->flux();
      // Kernel phase space of the emitted packet: at the surface step for a
      // retarding/zero field, at the Schottky-Nordheim matching plane (outer
      // turning point of the lowest occupied energy, widened by 1e-12) for an
      // extracting field. The packet is born on the face; the sub-cell matching
      // displacement is not represented (its field work is recorded).
      se::MatchingPlane plane{se::PlaneKind::ThermalVacuumSurfaceStep, 0, true, true, true};
      if (flux.spectrum.extraction_field_V_m > 0)
        plane = {se::PlaneKind::ExplicitSchottkyNordheimHandoff, se::required_SN_matching_distance_m(flux.spectrum) * (1 + 1e-12), true, true, true};
      const se::ConditionalPacket c{se::sample(flux, plane, b.weight, b.marks[0], b.marks[1], b.marks[2], b.ordinal + 1)};
      const auto uv = electron_position_uniforms(*o.arrival_policy, b);
      const auto& gg = heat.geometry();
      const auto local = face_point_local(gg, heat.patch(r.local), uv[0], uv[1]);
      av::BirthKinematics k;
      k.position_m = native_face_point(options, r, local, hb::slice_center_z(mesh, r.k));
      k.proper_velocity_m_s = c.particle.proper_velocity_m_s;
      k.origin = "coordinator-v6/face-emission/" + r.id + "/ordinal-" + std::to_string(b.ordinal);
      electron_kinematics.emplace(b.ordinal, std::move(k));
      // Nottingham energy of the emitted packet: a surface load of the one thermal owner.
      eh::Transfer x;
      x.producer = "coordinator-v6-emission/" + r.id + "/" + std::to_string(b.index);
      x.energy_owner = x.producer;
      x.state_id = "electron";
      x.parent_id = "electrode-electron-reservoir";
      x.channel_id = "electron-emission-nottingham";
      x.sequence = x.owner_sequence = heat.state().revision + 1;
      x.event_time_s = b.time_s;
      x.surface = true;
      x.x = local[0];
      x.y = local[1];
      x.authored_heat_J = c.particle.ledger.wall_heat_J;
      other[r.slice].push_back(std::move(x));
      out.numbers_.emission_nottingham_J += c.particle.ledger.wall_heat_J;
      out.numbers_.emission_reference_J += c.particle.ledger.electron_reference_change_J;
      out.numbers_.emission_matching_work_J += c.particle.ledger.resolved_matching_field_work_J;
      out.numbers_.emission_kinetic_J += c.particle.ledger.outgoing_K_J;
      ++out.numbers_.electron_births;
    }
  }
  out.times_.emission_births_s = seconds_since(t);
  t = clock::now();
  // Source trials: quantized neutral leaves of exactly w0 per patch (canonical
  // order). A slice owner with no leaf and no emission load this step has no
  // trial: its exact stock, layers and cumulative energies are unchanged (the
  // quantized remainder stays in the stock) and it is neither preflighted nor
  // published; its live digest is still re-recorded after the step.
  out.source_trials_.resize(nslices);
  for (std::size_t i = 0; i < nslices; ++i) {
    std::vector<std::vector<s6::NeutralLeaf>> leaves(o.slices[i].patches.size());
    bool any = !other[i].empty();
    for (const auto& b : staged.births()) {
      if (b.kind != pa::Kind::neutral) continue;
      const auto& r = *by_neutral.at(b.stream_id);
      if (r.slice != i) continue;
      leaves[r.local].push_back({b.weight, eh::Draw{b.marks[0], b.marks[1], b.marks[2], b.marks[3], b.marks[4], b.marks[5]}});
      any = true;
    }
    if (any) out.source_trials_[i].emplace(o.sources[i]->prepare(epochs[i], out.source_field_, view, {}, {}, leaves, other[i]));
  }
  for (const auto& r : o.patches) {
    std::uint64_t n = 0;
    for (const auto& b : staged.births()) n += b.stream_id == r.neutral_stream;
    if (n) staged.bind_source_times(r.neutral_stream, pas::neutral_times(*out.source_trials_[r.slice], staged, r.neutral_stream, r.id));
  }
  // Exact source debits = arrival decisions (n w0); neutral: the source inventory decrease.
  for (const auto& d : staged.plan().decisions()) {
    if (d.kind != pa::Kind::neutral) continue;
    const auto& r = *by_neutral.at(d.stream_id);
    if (d.debit.zero()) {
      need(d.births == 0, "nonzero release without a debit");
      continue;
    }
    need(out.source_trials_[r.slice].has_value(), "released molecules without a source trial");
    auto after = out.source_trials_[r.slice]->inventory(r.local);
    for (double w : d.debit.positive_carriers()) after = after.credited(w);
    need(after == o.sources[r.slice]->inventory(r.local), "source inventory decrease differs from the arrival debit (exact)");
    out.numbers_.source_debit = exact_add(out.numbers_.source_debit, d.debit);
  }
  out.times_.source_trials_s = seconds_since(t);
  fail(5);
  t = clock::now();
  // Kinematics in staged canonical order.
  std::vector<av::BirthKinematics> kinematics;
  {
    std::map<std::uint32_t, std::size_t> next_packet;
    for (const auto& b : staged.births()) {
      if (b.kind == pa::Kind::electron) { kinematics.push_back(electron_kinematics.at(b.ordinal)); continue; }
      const auto& r = *by_neutral.at(b.stream_id);
      const auto& packets = out.source_trials_[r.slice]->packets();
      std::size_t& at = next_packet[b.stream_id];
      std::size_t seen = 0;
      const eh::Packet* p = nullptr;
      for (const auto& x : packets)
        if (x.patch == r.id) { if (seen == at) { p = &x; break; } ++seen; }
      need(p && std::bit_cast<std::uint64_t>(p->weight) == std::bit_cast<std::uint64_t>(b.weight) &&
               std::bit_cast<std::uint64_t>(p->time) == std::bit_cast<std::uint64_t>(b.time_s),
           "source neutral packet differs from its staged leaf");
      ++at;
      av::BirthKinematics k;
      k.position_m = native_face_point(options, r, p->x, hb::slice_center_z(mesh, r.k));
      k.proper_velocity_m_s = p->proper_velocity_m_s;
      k.origin = "source-v6/" + p->event_id;
      kinematics.push_back(std::move(k));
      out.numbers_.neutral_kinetic_J += p->kinetic_J;
      out.numbers_.neutral_internal_J += p->internal_J;
      out.numbers_.neutral_adsorbed_reference_J += p->adsorbed_reference_change_J;
      ++out.numbers_.neutral_births;
    }
  }
  // Source transition seal: accepted seal, step, arrival plan/token, every slice trial.
  {
    Encoding e;
    e.text("coordinator-v6-source-transition-v1");
    e.text(before.source); e.u64(common.accepted_step); e.text(staged.plan().sha256()); e.text(staged.token());
    e.u64(nslices);
    for (std::size_t i = 0; i < nslices; ++i) {
      if (out.source_trials_[i]) e.text(out.source_trials_[i]->content_sha256());
      else { e.text("no-trial"); e.text(o.accepted_source_states[i]); }
    }
    out.source_after_ = e.finish();
  }
  const auto credit = av::credit(*o.arrival_policy, staged, kinematics, av::plan_debits(staged.plan()), before.source, out.source_after_);
  out.births_.emplace(scope.prepare(credit));
  auto motion_trial = out.outer_->fork_transport_trial(*out.births_);
  out.times_.birth_credit_s = seconds_since(t);
  fail(6);
  t = clock::now();
  out.census_ = motion_trial.native_census(o.inputs.maximum_carriers);
  require_census_species(out.census_, o.inputs.states);
  for (const auto& f : motion_trial.birth_flights()) out.newborns_.push_back({f.native_id, f.birth_time_s});
  for (std::size_t i = 0; i < nslices; ++i) {
    if (!out.source_trials_[i]) continue;
    auto L = hb::loads(*out.source_trials_[i], o.slices[i]);
    for (auto& l : L) {
      out.numbers_.source_load_J += hc::transfer_energy_J(l.transfer);
      out.loads_.push_back(std::move(l));
    }
  }
  h::MotionBinding motion;
  motion.phase = h::ParticlePhase::full_time_position_and_proper_velocity;
  motion.accepted_step = common.accepted_step;
  motion.staged_arena_revision = motion_trial.staged_arena_revision();
  motion.begin_s = common.begin_s;
  motion.end_s = common.end_s;
  motion.dt_s = common.dt_s;
  motion.staged_arena_sha256 = motion_trial.state_prefix();
  motion.mesh_sha256 = before.mesh;
  motion.material_sha256 = before.material;
  motion.species_sha256 = o.inputs.species_sha256;
  motion.common_clock_sha256 = common.common_clock_sha256;
  motion.kinetic_policy_sha256 = o.hall->motion_context().policy_sha256;
  out.numbers_.census = out.census_.size();
  out.numbers_.newborns = out.newborns_.size();
  out.times_.census_s = seconds_since(t);
  fail(7);
  t = clock::now();
  out.field_heat_ = std::make_unique<hc::JointProposal>(hc::prepare_closed_proposal(
      *o.hall, options, out.census_, motion, out.source_after_, *o.motion, hc::StepInputs{out.newborns_, out.loads_, true}));
  out.times_.hall_proposal_s = seconds_since(t);
  out.times_.hall_iterations = out.field_heat_->iterations();
  out.times_.hall = out.field_heat_->timings();
  fail(8);
  t = clock::now();
  const auto& proposal = *out.field_heat_;
  need(proposal.motion().frame().native_to_field_m == native_to_field_shift(options), "Hall motion translation differs from the derived frame");
  const auto midpoint_digest = h::total_field_sha256(proposal.midpoint().E(), proposal.midpoint().B());
  need(proposal.motion().binding().trial_field_sha256 == midpoint_digest && midpoint_digest != proposal.midpoint().seal(),
       "motion is not bound to the actual authoritative midpoint digest");
  need(proposal.births().size() == out.newborns_.size(), "Hall motion did not move every newborn");
  const auto& field_after = proposal.candidate_field_sha256();
  const auto& heat_after = proposal.candidate_heat_sha256();
  out.times_.binding_checks_s = seconds_since(t);
  t = clock::now();
  m::FullTimeMotionBinding native_motion{common, proposal.motion().record_v3_sha256(), before.field, field_after, heat_after, out.source_after_};
  auto moved = motion_trial.prepare_total_field_motion(native_motion, hn::bind_hall_v3_motion(proposal.motion()));
  fail(9);
  out.outer_->accept_transport_trial(std::move(moved));
  out.times_.motion_admission_s = seconds_since(t);
  t = clock::now();
  exposure.slot = m::FullTimeSlot::after_transport;
  exposure.frozen_physical_stage_s = common.end_s;
  exposure.owner_stage_sequence = first_sequence + 1;
  exposure.frozen_field_sha256 = field_after;
  auto second_scope = out.outer_->second_exposure(exposure);
  auto second = o.chemistry->prepare_on(second_scope, &first);
  out.times_.second_exposure_s = seconds_since(t);
  fail(10);
  t = clock::now();
  out.native_ = std::make_unique<m::PreparedArena>(out.outer_->seal({field_after, heat_after, out.source_after_, o.chemistry->staged_rng_sha256(second)}));
  out.random_ = o.chemistry->seal(first, second, *out.native_);
  const auto flows = out.native_->candidate_full_time_flows();
  need(flows && flows->births, "sealed native step without its birth receipt");
  out.arrival_receipt_ = av::receipt(*flows->births);
  out.times_.seal_s = seconds_since(t);
  fail(11);
  t = clock::now();
  out.numbers_.sinks = proposal.sinks().size();
  if (flows->sinks)
    for (const auto& tally : flows->sinks->tallies) {
      if (tally.state == "H2O") out.numbers_.sunk_h2o = exact_add(out.numbers_.sunk_h2o, tally.weight);
      if (tally.state == "H2O+") out.numbers_.sunk_h2o_ion = exact_add(out.numbers_.sunk_h2o_ion, tally.weight);
    }
  for (const auto* x : {&first.audit(), &second.audit()})
    if (x->contains("events") && x->at("events").is_array()) out.numbers_.parent_events += x->at("events").size();
  out.times_.ledger_s = seconds_since(t);
  fail(12);
  return out;
}

inline bool can_publish_v7(const OwnersV7& o, const TrialV7& trial) noexcept {
  int point = 0;
  const auto fault = [&] { ++point; if (trial.can_publish_fault_ == point) throw std::runtime_error("manufactured v6 can_publish fault"); };
  const auto refuse = [&](int code) { trial.refusal_point_ = code; return false; };
  trial.refusal_point_ = 0;
  using clock = std::chrono::steady_clock;
  auto& st = trial.sub_;
  st.cp_native_s = st.cp_chemistry_s = st.cp_observe_s = st.cp_predicate_s = st.cp_sources_s = st.cp_arrival_s = st.cp_hall_s = 0;
  st.cp_source_preflights = 0;
  auto t = clock::now();
  try {
    if (trial.consumed_ || !trial.native_ || !trial.field_heat_ || !trial.staged_ || trial.source_trials_.size() != o.sources.size()) return refuse(1);
    fault();
    if (!o.arena->can_publish(*trial.native_)) return refuse(2);
    st.cp_native_s = seconds_since(t); t = clock::now();
    fault();
    if (!o.chemistry->can_publish(trial.random_, *trial.native_)) return refuse(3);
    st.cp_chemistry_s = seconds_since(t); t = clock::now();
    fault();
    const auto boundary = trial.native_->candidate_full_time_boundary();
    const auto exposures = trial.native_->candidate_full_time_exposures();
    const auto flows = trial.native_->candidate_full_time_flows();
    const auto& proposal = *trial.field_heat_;
    if (!boundary || !exposures || !flows || !flows->births || exposures->transport_trial_sha256 != proposal.motion().record_v3_sha256())
      return refuse(4);
    fault();
    // Fresh observations of every owner (digests recomputed, live source/arrival states).
    const auto current = observe(o);
    st.cp_observe_s = seconds_since(t); t = clock::now();
    fault();
    cc::ProposedOwners candidate;
    auto& after = candidate.after;
    after.native = *boundary;
    after.actual_native_prefix = trial.native_->candidate_state_prefix();
    after.actual_native_revision = boundary->accepted_arena_revision;
    after.actual_native_step = boundary->accepted_step;
    after.field = proposal.candidate_field_sha256();
    after.heat = proposal.candidate_heat_sha256();
    after.source = trial.source_after_;
    after.rng = o.chemistry->candidate_rng_sha256(trial.random_, *trial.native_);
    after.mesh = current.mesh;
    after.material = current.material;
    after.frame = current.frame;
    after.configuration = current.configuration;
    after.field_time_s = after.heat_time_s = proposal.end_s();
    after.field_revision = after.heat_revision = proposal.endpoint_revision();
    fault();
    candidate.motion_midpoint = proposal.motion().binding().trial_field_sha256;
    candidate.field_midpoint = h::total_field_sha256(proposal.midpoint().E(), proposal.midpoint().B());
    if (candidate.field_midpoint == proposal.midpoint().seal()) return refuse(5);
    const auto& motion = proposal.motion().binding();
    if (proposal.motion().frame().native_to_field_m != native_to_field_shift(o.arena->native_options()) ||
        motion.species_sha256 != trial.inputs_.species_sha256 || motion.mesh_sha256 != current.mesh ||
        motion.kinetic_policy_sha256 != o.hall->motion_context().policy_sha256)
      return refuse(6);
    fault();
    c7::ArrivalEvidence evidence;
    evidence.staged = true;
    evidence.plan_sha256 = trial.staged_->plan().sha256();
    evidence.token_sha256 = trial.staged_->token();
    evidence.base = trial.staged_->native_base();
    evidence.births = trial.staged_->births().size();
    const bool predicate = c7::accepts_closed_or_birth_step(trial.expected_, current, candidate, flows.get(), evidence, trial.clock_);
    if (trial.compare_v6_predicate_) {
      c6::ArrivalEvidence e6;
      e6.staged = evidence.staged; e6.plan_sha256 = evidence.plan_sha256; e6.token_sha256 = evidence.token_sha256;
      e6.base = evidence.base; e6.births = evidence.births;
      trial.v6_predicate_ = c6::accepts_closed_or_birth_step(trial.expected_, current, candidate, flows.get(), e6) ? 1 : 0;
      trial.v7_predicate_ = predicate ? 1 : 0;
    }
    if (proposal.configuration_sha256() != o.hall->configuration_sha256() || !predicate) return refuse(7);
    st.cp_predicate_s = seconds_since(t); t = clock::now();
    fault();
    // Source owners: each trial at the accepted thermal clock and bound to the
    // live accepted Hall start-of-step arrays.
    const s6::FinalFieldView live{trial.source_frame_, o.hall->mesh(), o.hall->field().E, o.hall->field().B};
    for (std::size_t i = 0; i < o.sources.size(); ++i)
      if (trial.source_trials_[i]) {
        ++st.cp_source_preflights;
        if (!o.sources[i]->preflight(*trial.source_trials_[i], live)) return refuse(8);
      }
    st.cp_sources_s = seconds_since(t); t = clock::now();
    fault();
    if (!arrival_commit_admissible(*o.arrival, *trial.staged_, trial.arrival_receipt_)) return refuse(9);
    st.cp_arrival_s = seconds_since(t); t = clock::now();
    fault();
    // Final Hall approval is last.
    if (!o.hall->preflight_final(proposal, trial.census_, *boundary)) return refuse(10);
    st.cp_hall_s = seconds_since(t);
    return true;
  } catch (...) {
    trial.refusal_point_ = 100 + point;
    return false;
  }
}

// Publication: native -> chemistry -> Hall -> source owners -> coordinator source
// seal, with no I/O, callbacks or allocation between the first and last swap; then
// the arrival ledger commit and the record of the live owner digests.
inline bool try_publish_v7(OwnersV7& o, TrialV7& trial) noexcept {
  if (!can_publish_v7(o, trial)) return false;
  if (trial.publish_probe_) trial.publish_probe_(0);
  if (!o.arena->try_publish(*trial.native_)) { trial.refusal_point_ = 90; return false; }
  o.chemistry->publish_prevalidated(trial.random_);
  o.hall->publish_prevalidated(std::move(*trial.field_heat_));
  for (std::size_t i = 0; i < o.sources.size(); ++i)
    if (trial.source_trials_[i]) o.sources[i]->publish_prevalidated(std::move(*trial.source_trials_[i]));
  o.accepted_source_seal.swap(trial.source_after_);
  if (trial.publish_probe_) trial.publish_probe_(1);
  try {
    (void)o.arrival->commit(std::move(*trial.staged_), trial.arrival_receipt_);
    for (std::size_t i = 0; i < o.sources.size(); ++i) o.accepted_source_states[i] = o.sources[i]->state_sha256();
    o.accepted_arrival_state = o.arrival->state_sha256();
  } catch (...) {
    std::terminate();
  }
  trial.consumed_ = true;
  trial.native_.reset();
  trial.field_heat_.reset();
  trial.staged_.reset();
  trial.source_trials_.clear();
  trial.births_.reset();
  trial.scope_.reset();
  trial.outer_.reset();
  return true;
}
}  // namespace raikou::native_water::v7
