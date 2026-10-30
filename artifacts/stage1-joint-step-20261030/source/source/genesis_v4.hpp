#pragma once
#include "closed_parent_joint_step_v4.hpp"

// Joint fresh genesis with A09 admission. Order:
//  1 model/time/geometry admission     6 genesis binding (field, heat, mesh,
//  2 predicted native census             material, frame, inactive source,
//  3 declared DC ledger + census charge  RNG, genesis clock)
//  4 provisional Hall owner (absolute    7 native create_full_time
//    Gauss/divB against the prediction)  8 actual census == prediction
//  5 chemistry configuration and RNG     9 real Hall owner from actual census
//                                        10 chemistry owner; all owners agree
// Nothing is repaired: every inconsistency refuses before an owner is kept.
namespace raikou::native_water::v4 {
struct GenesisTestAccess;

struct GenesisRequestV4 {
  std::shared_ptr<const hc::Model> model;
  m::GapPicOptions options;
  axial::State seed;
  std::vector<m::ArenaStateSpec> states;
  std::vector<m::ArenaChannel> channels;
  std::vector<m::FullTimeAuthoredPacket> packets;
  std::array<std::uint32_t,3> grid_counts{};
  Binding binding;
  Limits limits;
  std::uint64_t clock_seed=0,product_seed=0;
  std::string origin_provenance,conductor_declaration,exterior_declaration;
  std::size_t maximum_native_records=4096,maximum_carriers=100000;
};

// Member order fixes destruction order: chemistry (which references the
// arena) is destroyed before the arena.
struct JointOwnersV4 {
  std::shared_ptr<const hc::Model> model;
  std::unique_ptr<m::PicArena> arena;
  std::unique_ptr<hc::Owner> owner;
  std::unique_ptr<ParentSubflow> chemistry;
  ClosedParentInputsV4 inputs;
  Json genesis;
};

class JointGenesisV4 {
  friend struct GenesisTestAccess;
  using PacketHook=void(*)(std::vector<m::FullTimeAuthoredPacket>&);
  using LedgerHook=void(*)(hc::ChargeLedger&);

  // Exactly the native create_full_time carrier assignment: IDs from 1 over
  // positive_carriers in packet order, charge = charge_e x pic_e.
  static std::vector<m::ArenaRateEntry> predicted_census(const std::vector<m::ArenaStateSpec>& states,
      const std::vector<m::FullTimeAuthoredPacket>& packets,std::size_t cap){
    const auto index=state_index(states);
    std::vector<m::ArenaRateEntry> out;std::uint64_t id=1;
    for(const auto& p:packets){
      const auto it=index.find(p.state);
      pc::require(it!=index.end(),"genesis packet has an undeclared species");
      for(double w:m::PicArena::positive_carriers(p.count)){
        pc::require(out.size()<cap,"genesis census cap");
        m::GapParticle n;n.id=id++;n.position_m=p.position_m;n.proper_velocity_m_s=p.proper_velocity_m_s;
        n.charge_C=it->second->charge_e*m::pic_e;n.mass_kg=it->second->mass_kg;n.weight=w;
        out.push_back({p.id,p.state,n});
      }
    }
    return out;
  }
  // Declared DC banks: the seed's field charge on conductor/short nodes is
  // conductor charge and on the EM exterior is exterior termination charge.
  static hc::ChargeLedger declared_ledger(const hc::Model& model,const axial::State& seed,
      const std::string& conductor,const std::string& exterior){
    pc::require(!conductor.empty()&&!exterior.empty(),"explicit conductor/exterior bank declarations required");
    const auto& g=model.mesh();hc::ChargeLedger q;
    q.particles.net_C.resize(g.nn);q.particles.gross_C.resize(g.nn);q.conductor_C.resize(g.nn);q.exterior_C.resize(g.nn);
    q.conductor_declaration=hc::sha256(conductor);q.exterior_declaration=hc::sha256(exterior);
    const auto field=hc::field_charge(g,seed);
    for(int node=0;node<g.nn;++node){
      if(g.node_body[node]>=0||g.short_node[node])q.conductor_C[node]=field[node];
      else if(g.outer_node[node])q.exterior_C[node]=field[node];
    }
    return q;
  }
  static JointOwnersV4 admit_with_hooks(const GenesisRequestV4& rq,const w::ParentBridge& bridge,
      PacketHook packet_hook,LedgerHook ledger_hook,int* reached){
    const auto stage=[&](int s){if(reached)*reached=s;};
    stage(1);
    pc::require(bool(rq.model),"genesis requires an immutable Hall model");
    const auto& model=*rq.model;
    pc::require(rq.seed.time_s==model.drive().start(),"genesis time differs from the first supported C19 knot");
    pc::require(rq.seed.revision==0,"genesis field revision must be zero");
    pc::require(rq.maximum_carriers>0&&rq.maximum_carriers<=100000&&
      rq.maximum_native_records>=7&&rq.maximum_native_records<=4096,"genesis native resource caps");
    // Geometry guard: native options must describe the Hall field geometry.
    const auto mesh=h::motion_mesh_sha256(model.mesh(),rq.options);
    stage(2);
    const auto predicted=predicted_census(rq.states,rq.packets,rq.maximum_carriers);
    require_census_species(predicted,rq.states);
    const auto predicted_sha=h::canonical_census_sha256(predicted);
    stage(3);
    auto ledger=declared_ledger(model,rq.seed,rq.conductor_declaration,rq.exterior_declaration);
    ledger.particles=hc::charge_of_census(model.mesh(),rq.options,predicted);
    if(ledger_hook)ledger_hook(ledger);
    stage(4);
    std::string field,heat;
    {
      hc::Owner provisional(rq.model,rq.seed,ledger,rq.options,predicted);
      field=provisional.field_sha256();heat=provisional.heat_sha256();
    }
    stage(5);
    const auto chemistry_configuration=ParentSubflow::full_time_configuration_sha256(
      rq.options,rq.grid_counts,rq.binding,rq.limits);
    const auto rng=ParentSubflow::initial_full_time_rng_sha256(
      rq.options,rq.grid_counts,rq.binding,rq.limits,rq.clock_seed,rq.product_seed);
    ClosedParentInputsV4 inputs{rq.states,rq.channels,species_sha256(rq.states,rq.channels),
      frame_sha256(rq.options,model.mesh()),rq.maximum_native_records,rq.maximum_carriers};
    const auto config=closed_parent_configuration_sha256(chemistry_configuration,model.configuration_sha256(),inputs);
    stage(6);
    m::FullTimeGenesisBinding g;
    g.time_s=rq.seed.time_s;g.origin_provenance=rq.origin_provenance;
    g.field_sha256=field;g.heat_sha256=heat;g.mesh_sha256=mesh;
    g.material_sha256=h::material_binding_sha256();g.frame_sha256=inputs.frame_sha256;
    g.source_sha256=closed_parent_inactive_source_sha256(config);g.rng_sha256=rng;
    {
      Encoding e;e.text("closed-parent-genesis-clock-v4");e.f64(g.time_s);e.f64(model.drive().start());
      for(const std::string* x:std::initializer_list<const std::string*>{&config,&field,&heat,&mesh,&g.material_sha256,
          &g.frame_sha256,&g.source_sha256,&rng,&predicted_sha})e.text(*x);
      g.common_clock_sha256=e.finish();
    }
    stage(7);
    auto packets=rq.packets;
    if(packet_hook)packet_hook(packets);
    JointOwnersV4 out;out.model=rq.model;
    out.arena=m::PicArena::create_full_time(rq.options,rq.states,rq.channels,packets,g);
    stage(8);
    const auto actual=out.arena->rate_population(rq.maximum_carriers);
    pc::require(h::canonical_census_sha256(actual)==predicted_sha,
      "actual native genesis census differs from the provisional Hall prediction");
    stage(9);
    auto actual_ledger=declared_ledger(model,rq.seed,rq.conductor_declaration,rq.exterior_declaration);
    actual_ledger.particles=hc::charge_of_census(model.mesh(),rq.options,actual);
    const auto constraints=model.admit(rq.seed,actual_ledger);
    out.owner=std::make_unique<hc::Owner>(rq.model,rq.seed,actual_ledger,rq.options,actual);
    pc::require(out.owner->field_sha256()==field&&out.owner->heat_sha256()==heat,
      "actual-census Hall bootstrap differs from the provisional genesis seals");
    stage(10);
    out.chemistry=std::make_unique<ParentSubflow>(*out.arena,rq.grid_counts,rq.binding,rq.limits,
      bridge,rq.clock_seed,rq.product_seed,true);
    out.inputs=inputs;
    const auto observed=observe_closed_accepted(*out.arena,*out.chemistry,*out.owner,out.inputs);
    pc::require(common_commit::owner_values_match(observed)&&observed.configuration==config&&
      observed.native.authored_full_time_genesis,"genesis owners disagree with the accepted native boundary");
    stage(11);
    out.genesis={{"schema","closed-parent-joint-genesis-v4"},{"time_s",g.time_s},
      {"origin_provenance",g.origin_provenance},{"field_sha256",field},{"heat_sha256",heat},
      {"mesh_sha256",mesh},{"material_sha256",g.material_sha256},{"frame_sha256",g.frame_sha256},
      {"species_sha256",inputs.species_sha256},{"source_sha256",g.source_sha256},{"rng_sha256",rng},
      {"common_clock_sha256",g.common_clock_sha256},{"configuration_sha256",config},
      {"chemistry_configuration_sha256",chemistry_configuration},
      {"field_configuration_sha256",model.configuration_sha256()},
      {"predicted_census_sha256",predicted_sha},{"actual_census_sha256",h::canonical_census_sha256(actual)},
      {"census_rows",actual.size()},{"native_state_prefix",out.arena->state_prefix()},
      {"admitted_constraints",{{"max_gauss_defect_C",constraints.max_gauss_C},
        {"max_gauss_bound_C",constraints.max_gauss_bound_C},{"max_divB_T_m",constraints.max_divB_T_m},
        {"max_divB_bound_T_m",constraints.max_divB_bound_T_m},
        {"initial_gauss_solver_bound_C",model.initial_gauss_solve_bound_C()}}}};
    return out;
  }
public:
  static JointOwnersV4 admit(const GenesisRequestV4& rq,const w::ParentBridge& bridge){
    return admit_with_hooks(rq,bridge,nullptr,nullptr,nullptr);
  }
};
} // namespace raikou::native_water::v4
