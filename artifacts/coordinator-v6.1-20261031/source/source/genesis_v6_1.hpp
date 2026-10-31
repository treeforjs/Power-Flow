#pragma once
#include "joint_step_v6_1.hpp"

// Joint fresh genesis with A09 admission (coordinator v6). Order:
//  1 model/time/geometry admission       6 genesis binding (field, heat, mesh,
//  2 predicted native census               material, frame, genesis source
//  3 declared DC ledger + census charge    seal, RNG, genesis clock)
//  4 provisional Hall owner (absolute      7 native create_full_time (caps)
//    Gauss/divB against the prediction),   8 actual census == prediction
//    provisional source owners bound to    9 real Hall owner from the actual census
//    its thermal slices, arrival owner    10 real source owners (same configuration
//  5 chemistry configuration and RNG         and state digests), arrival owner,
//    (parent coverage v2 scope)              chemistry; all owners agree
// Nothing is repaired: every inconsistency refuses before an owner is kept.
namespace raikou::native_water::v6_1 {
struct GenesisRequestV6 {
  std::shared_ptr<const hc::Model> model;
  m::GapPicOptions options;
  axial::State seed;
  std::vector<m::ArenaStateSpec> states;
  std::vector<m::ArenaChannel> channels;
  std::vector<m::FullTimeAuthoredPacket> packets;
  std::array<std::uint32_t, 3> grid_counts{};
  Binding binding;
  Limits limits;
  std::shared_ptr<const w::ParentBridge> bridge;
  std::shared_ptr<const cv::CoverageScope> scope;
  std::uint64_t clock_seed = 0, product_seed = 0;
  std::string origin_provenance, conductor_declaration, exterior_declaration;
  PilotSettings settings;
  SourceSetup sources;
};

class GenesisV6 {
  friend struct GenesisTestAccessV6;
 public:
  using PacketHook = void (*)(std::vector<m::FullTimeAuthoredPacket>&);
  using LedgerHook = void (*)(hc::ChargeLedger&);
 private:
  static std::vector<m::ArenaRateEntry> predicted_census(const std::vector<m::ArenaStateSpec>& states,
                                                        const std::vector<m::FullTimeAuthoredPacket>& packets, std::size_t cap) {
    const auto index = state_index(states);
    std::vector<m::ArenaRateEntry> out;
    std::uint64_t id = 1;
    for (const auto& p : packets) {
      const auto it = index.find(p.state);
      need(it != index.end(), "genesis packet has an undeclared species");
      for (double wgt : m::PicArena::positive_carriers(p.count)) {
        need(out.size() < cap, "genesis census cap");
        m::GapParticle n;
        n.id = id++; n.position_m = p.position_m; n.proper_velocity_m_s = p.proper_velocity_m_s;
        n.charge_C = it->second->charge_e * m::pic_e; n.mass_kg = it->second->mass_kg; n.weight = wgt;
        out.push_back({p.id, p.state, n});
      }
    }
    return out;
  }
  static hc::ChargeLedger declared_ledger(const hc::Model& model, const axial::State& seed, const std::string& conductor,
                                          const std::string& exterior) {
    need(!conductor.empty() && !exterior.empty(), "explicit conductor/exterior bank declarations required");
    const auto& g = model.mesh();
    hc::ChargeLedger q;
    q.particles.net_C.resize(g.nn); q.particles.gross_C.resize(g.nn); q.conductor_C.resize(g.nn); q.exterior_C.resize(g.nn);
    q.conductor_declaration = hc::sha256(conductor); q.exterior_declaration = hc::sha256(exterior);
    const auto field = hc::field_charge(g, seed);
    for (int node = 0; node < g.nn; ++node) {
      if (g.node_body[node] >= 0 || g.short_node[node]) q.conductor_C[node] = field[node];
      else if (g.outer_node[node]) q.exterior_C[node] = field[node];
    }
    return q;
  }
  struct SourceProbe {
    std::vector<std::string> signatures, states;
    std::string arrival_state;
  };
  static SourceProbe probe_sources(const hc::Owner& hall, const GenesisRequestV6& rq) {
    SourceProbe p;
    std::vector<std::unique_ptr<s6::Owner>> owners;
    for (const auto& s : rq.sources.slices) owners.push_back(hb::bind(hall, s));
    for (const auto& x : owners) { p.signatures.push_back(x->configuration_sha256()); p.states.push_back(x->state_sha256()); }
    pa::Owner arrival(rq.sources.arrival, 0);
    p.arrival_state = arrival.state_sha256();
    return p;
  }
  static OwnersV6 admit_with_hooks(const GenesisRequestV6& rq, PacketHook packet_hook, LedgerHook ledger_hook, int* reached) {
    const auto stage = [&](int s) { if (reached) *reached = s; };
    stage(1);
    need(bool(rq.model) && bool(rq.bridge) && bool(rq.scope) && bool(rq.sources.arrival), "genesis requires model, bridge, scope and arrival policy");
    const auto& model = *rq.model;
    need(rq.seed.time_s == model.drive().start(), "genesis time differs from the first supported C19 knot");
    need(rq.seed.revision == 0, "genesis field revision must be zero");
    const auto& S = rq.settings;
    need(S.maximum_carriers > 0 && S.maximum_carriers <= model.limits().maximum_carriers && S.maximum_native_records >= 2 * rq.limits.reactions + 6 &&
             S.maximum_native_records <= 4096,
         "genesis native resource caps (maximum_records >= 2M+6 for births/sinks)");
    need(rq.sources.slices.size() == std::size_t(2 * model.mesh().nz_foil), "one source owner per Hall foil slice");
    const auto mesh = model.motion_mesh_sha256(rq.options);
    stage(2);
    const auto predicted = predicted_census(rq.states, rq.packets, S.maximum_carriers);
    require_census_species(predicted, rq.states);
    const auto predicted_sha = hc::canonical_census_sha256(predicted, S.maximum_carriers);
    stage(3);
    auto ledger = declared_ledger(model, rq.seed, rq.conductor_declaration, rq.exterior_declaration);
    ledger.particles = hc::charge_of_census(model.mesh(), rq.options, predicted, model.limits().maximum_carriers);
    if (ledger_hook) ledger_hook(ledger);
    stage(4);
    std::string field, heat;
    SourceProbe probe;
    {
      hc::Owner provisional(rq.model, rq.seed, ledger, rq.options, predicted);
      field = provisional.field_sha256();
      heat = provisional.heat_sha256();
      probe = probe_sources(provisional, rq);
    }
    stage(5);
    const auto chemistry_configuration = cv::ParentSubflow::full_time_configuration_sha256(rq.options, rq.grid_counts, rq.binding, rq.limits, *rq.scope);
    const auto rng = cv::ParentSubflow::initial_full_time_rng_sha256(rq.options, rq.grid_counts, rq.binding, rq.limits, *rq.scope, rq.clock_seed,
                                                                     rq.product_seed);
    InputsV6 inputs;
    inputs.states = rq.states;
    inputs.channels = rq.channels;
    inputs.species_sha256 = species_sha256(rq.states, rq.channels);
    inputs.frame_sha256 = frame_sha256(rq.options, model);
    inputs.source_configuration_sha256 = source_configuration_sha256(rq.sources.slices, probe.signatures);
    inputs.settings_text = S.json_text();
    inputs.maximum_native_records = S.maximum_native_records;
    inputs.maximum_carriers = S.maximum_carriers;
    const auto config = configuration_sha256(chemistry_configuration, model.configuration_sha256(), inputs, rq.sources.arrival->configuration_sha256());
    const auto source_seal = genesis_source_seal(config, probe.states, probe.arrival_state, rq.sources.arrival->configuration_sha256());
    stage(6);
    m::FullTimeGenesisBinding g;
    g.time_s = rq.seed.time_s;
    g.origin_provenance = rq.origin_provenance;
    g.field_sha256 = field; g.heat_sha256 = heat; g.mesh_sha256 = mesh;
    g.material_sha256 = model.material_sha256(); g.frame_sha256 = inputs.frame_sha256;
    g.source_sha256 = source_seal; g.rng_sha256 = rng;
    {
      Encoding e;
      e.text("coordinator-v6-genesis-clock-v1");
      e.f64(g.time_s); e.f64(model.drive().start());
      for (const std::string* x : std::initializer_list<const std::string*>{&config, &field, &heat, &mesh, &g.material_sha256, &g.frame_sha256, &g.source_sha256, &rng, &predicted_sha})
        e.text(*x);
      g.common_clock_sha256 = e.finish();
    }
    stage(7);
    auto packets = rq.packets;
    if (packet_hook) packet_hook(packets);
    OwnersV6 out;
    out.model = rq.model;
    {
      m::FullTimeConfiguration configuration;
      configuration.limits.maximum_carriers = S.maximum_carriers;
      configuration.limits.memory_budget_bytes = S.native_memory_budget;
      out.arena = m::PicArena::create_full_time(rq.options, rq.states, rq.channels, packets, g, configuration);
    }
    stage(8);
    const auto actual = out.arena->rate_population(S.maximum_carriers);
    need(hc::canonical_census_sha256(actual, S.maximum_carriers) == predicted_sha, "actual native genesis census differs from the provisional Hall prediction");
    stage(9);
    auto actual_ledger = declared_ledger(model, rq.seed, rq.conductor_declaration, rq.exterior_declaration);
    actual_ledger.particles = hc::charge_of_census(model.mesh(), rq.options, actual, model.limits().maximum_carriers);
    const auto constraints = model.admit(rq.seed, actual_ledger);
    out.hall = std::make_unique<hc::Owner>(rq.model, rq.seed, actual_ledger, rq.options, actual);
    need(out.hall->field_sha256() == field && out.hall->heat_sha256() == heat, "actual-census Hall bootstrap differs from the provisional genesis seals");
    stage(10);
    out.slices = rq.sources.slices;
    out.patches = rq.sources.patches;
    for (const auto& s : out.slices) out.sources.push_back(hb::bind(*out.hall, s));
    for (std::size_t i = 0; i < out.sources.size(); ++i) {
      need(out.sources[i]->configuration_sha256() == probe.signatures[i] && out.sources[i]->state_sha256() == probe.states[i],
           "real source owner differs from its provisional genesis digest");
      out.accepted_source_states.push_back(out.sources[i]->state_sha256());
    }
    out.arrival_policy = rq.sources.arrival;
    out.arrival = std::make_unique<pa::Owner>(out.arrival_policy, out.arena->full_time_boundary()->accepted_step);
    need(out.arrival->state_sha256() == probe.arrival_state, "arrival owner differs from its provisional genesis digest");
    out.accepted_arrival_state = out.arrival->state_sha256();
    out.accepted_source_seal = source_seal;
    out.chemistry = std::make_unique<cv::ParentSubflow>(*out.arena, rq.grid_counts, rq.binding, rq.limits, *rq.bridge, *rq.scope, rq.clock_seed,
                                                        rq.product_seed, true);
    out.inputs = inputs;
    out.settings = S;
    const auto observed = observe(out);
    need(cc::owner_values_match(observed) && observed.configuration == config && observed.native.authored_full_time_genesis,
         "genesis owners disagree with the accepted native boundary");
    stage(11);
    Json j = Json::object();
    j["schema"] = "coordinator-v6-genesis-v1";
    j["time_s"] = g.time_s; j["origin_provenance"] = g.origin_provenance; j["field_sha256"] = field; j["heat_sha256"] = heat;
    j["mesh_sha256"] = mesh; j["material_sha256"] = g.material_sha256; j["frame_sha256"] = g.frame_sha256;
    j["species_sha256"] = inputs.species_sha256; j["source_sha256"] = source_seal; j["rng_sha256"] = rng;
    j["common_clock_sha256"] = g.common_clock_sha256; j["configuration_sha256"] = config;
    j["chemistry_configuration_sha256"] = chemistry_configuration; j["field_configuration_sha256"] = model.configuration_sha256();
    j["source_configuration_sha256"] = inputs.source_configuration_sha256;
    j["arrival_policy_sha256"] = rq.sources.arrival->configuration_sha256();
    j["arrival_sampling_policy_identity"] = rq.sources.arrival->sampling_policy_identity();
    j["settings"] = S.json_text();
    j["field_solve_backend"] = model.solver().name();
    j["field_solve_backend_descriptor_sha256"] = model.solver().descriptor_sha256();
    j["predicted_census_sha256"] = predicted_sha;
    j["census_rows"] = actual.size();
    j["native_state_prefix"] = out.arena->state_prefix();
    j["native_source_identity"] = m::PicArena::source_identity();
    j["source_owners"] = out.sources.size();
    j["source_patches"] = out.patches.size();
    j["admitted_constraints"] = {{"max_gauss_defect_C", constraints.max_gauss_C}, {"max_gauss_bound_C", constraints.max_gauss_bound_C},
                                 {"max_divB_T_m", constraints.max_divB_T_m}, {"max_divB_bound_T_m", constraints.max_divB_bound_T_m},
                                 {"initial_gauss_solver_bound_C", model.initial_gauss_solve_bound_C()}};
    out.genesis = std::move(j);
    return out;
  }

 public:
  static OwnersV6 admit(const GenesisRequestV6& rq) { return admit_with_hooks(rq, nullptr, nullptr, nullptr); }
};
}  // namespace raikou::native_water::v6_1
