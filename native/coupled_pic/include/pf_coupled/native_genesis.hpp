#pragma once
// Development integration: requires the new native-v8 SourceStage interface.
// This builder never aliases the legacy parent-chemistry species catalogue.
#include "source_stage_v8.hpp"
#include "pf_gas_v1/gas_owner.hpp"
#include "native_gas_stage.hpp"
#include "raikou/mfem_pic_canonical.h"
#include <bit>
#include <map>
#include <set>

namespace powerflow::coupled_genesis_v9 {
namespace ss = ::powerflow::source_stage_v8;
namespace n = ::raikou::mfem_backend;
namespace h = ::hall_corrective_v3;
namespace g = ::pf_gas_v1;
inline void require(bool x,const char *why) {
  if (!x) throw std::invalid_argument(std::string("coupled_genesis_v9: ")+why);
}
struct Encoding {
  std::string bytes;
  void u64(std::uint64_t x) { for(int k=0;k<8;++k)bytes.push_back(char((x>>(8*k))&255)); }
  void number(double x) { u64(std::bit_cast<std::uint64_t>(x)); }
  void text(const std::string &s) { u64(s.size());bytes+=s; }
  std::string hash()const {return n::canonical::sha256_hex(bytes);}
};
struct PortSpec {
  ss::hb::SliceSource slice;
  ss::s7::OwnerOptions options;
  std::vector<ss::sr::PatchRegistration> registration;
};
struct Request {
  std::shared_ptr<const h::Model> model;
  n::GapPicOptions options;
  n::FullTimeConfiguration native_configuration;
  g::Configuration gas_configuration;
  axial::State seed;
  std::vector<n::FullTimeAuthoredPacket> packets;
  std::vector<PortSpec> ports;
  std::shared_ptr<const ss::pa::Policy> arrival_policy;
  std::string origin_provenance,conductor_declaration,exterior_declaration;
};
struct Owners {
  std::shared_ptr<const h::Model> model;
  std::shared_ptr<h::Owner> hall;
  ss::OwnerSet sources;
  std::shared_ptr<g::Owner> gas;
  std::shared_ptr<const native_gas_stage::Binding> gas_binding;
  std::unique_ptr<n::PicArena> native;
  n::FullTimeConfiguration native_configuration;
  n::FullTimeChemistryRegistry registry;
  std::string configuration_sha256,source_sha256,arrival_state_sha256;
  std::vector<std::string> source_states;
};
inline std::vector<n::ArenaRateEntry> predict(const n::FullTimeChemistryRegistry &registry,
                                             const std::vector<n::FullTimeAuthoredPacket> &packets,
                                             std::size_t cap) {
  const auto states=n::full_time_chemistry_states(registry);
  std::map<std::string,const n::ArenaStateSpec*> by_state;
  for(const auto&s:states)require(by_state.emplace(s.id,&s).second,"duplicate actual registry state");
  std::set<std::string> ids;std::vector<n::ArenaRateEntry> result;std::uint64_t next=1;
  for(const auto&p:packets) {
    const auto it=by_state.find(p.state);
    require(it!=by_state.end() && ids.insert(p.id).second,"unknown/duplicate authored packet");
    for(double w:n::PicArena::positive_carriers(p.count)) {
      require(result.size()<cap && next<UINT64_MAX,"genesis carrier capacity");
      n::GapParticle r;r.id=next++;r.position_m=p.position_m;r.proper_velocity_m_s=p.proper_velocity_m_s;
      r.weight=w;r.mass_kg=it->second->mass_kg;r.charge_C=it->second->charge_e*n::pic_e;
      result.push_back({p.id,p.state,r});
    }
  }
  return result;
}
inline h::ChargeLedger charge(const Request &r,std::span<const n::ArenaRateEntry> predicted) {
  require(!r.conductor_declaration.empty() && !r.exterior_declaration.empty(),"explicit boundary charge declarations required");
  const auto&m=r.model->mesh();h::ChargeLedger q;
  q.particles=h::charge_of_census(m,r.options,predicted,r.model->limits().maximum_carriers);
  q.conductor_C.resize(m.nn);q.exterior_C.resize(m.nn);
  q.conductor_declaration=h::sha256(r.conductor_declaration);q.exterior_declaration=h::sha256(r.exterior_declaration);
  const auto fq=h::field_charge(m,r.seed);
  for(int i=0;i<m.nn;++i) {
    if(m.node_body[i]>=0 || m.short_node[i])q.conductor_C[i]=fq[i];
    else if(m.outer_node[i])q.exterior_C[i]=fq[i];
  }
  return q;
}
inline Owners admit(const Request &r) {
  require(r.model && r.arrival_policy,"model and arrival policy required");
  require(r.seed.time_s==r.model->drive().start() && r.seed.revision==0,
          "fresh genesis must begin at the supported drive start with revision zero");
  require(r.native_configuration.identity==n::FullTimeIdentitySchema::v8 && r.native_configuration.chemistry.enabled,
          "actual full-native chemistry and identity v8 required");
  require(!r.native_configuration.chemistry.diagnostic_field_target_approximation,
          "unowned diagnostic target energy cannot enter common genesis");
  require(!r.native_configuration.merge_split.enabled,
          "merge/split gas transport is not qualified for this common genesis");
  require(r.arrival_policy->configuration().run_seed==r.native_configuration.chemistry.run_seed,
          "source and collision RNG run seed mismatch");
  require(r.gas_configuration.feedback && r.gas_configuration.transport==g::Transport::census_conservative &&
          r.gas_configuration.execution==g::Execution::host,
          "qualified host census transport and gas feedback required");
  require(!r.origin_provenance.empty(),"declared genesis provenance required");
  const auto cap=r.native_configuration.limits.maximum_carriers;
  require(cap>0 && cap<=r.model->limits().maximum_carriers,"joint carrier capacity");
  Owners out;out.model=r.model;out.native_configuration=r.native_configuration;
  out.registry=n::full_time_chemistry_registry(r.native_configuration.chemistry,r.options);
  require(r.gas_configuration.chemistry_registry_sha256==out.registry.registry_sha256 &&
          r.gas_configuration.chemistry_blob_sha256==out.registry.blob_sha256,
          "actual native/gas registry and blob mismatch");
  require(r.gas_configuration.grid.n==r.native_configuration.chemistry.grid_counts &&
          r.gas_configuration.grid.low==std::array<double,3>{0,0,r.options.short_z_m} &&
          r.gas_configuration.grid.high==std::array<double,3>{r.options.width_m,r.options.gap_m,r.options.feed_z_m},
          "actual native/gas collision grid mismatch");
  require(r.gas_configuration.species.size()>=out.registry.species.size(),"gas species catalogue is incomplete");
  for(std::size_t i=0;i<out.registry.species.size();++i) {
    const auto&a=out.registry.species[i];const auto&b=r.gas_configuration.species[i];
    require(a.state==b.state && std::bit_cast<std::uint64_t>(a.mass_kg)==std::bit_cast<std::uint64_t>(b.mass_kg) &&
            a.charge_e==b.charge && a.target==b.chemistry_target,"gas species differ from actual native operands");
  }
  out.gas=std::make_shared<g::Owner>(r.gas_configuration);
  const auto blob=pfc::load_blob(r.native_configuration.chemistry.blob_path,out.registry.blob_sha256);
  out.gas_binding=std::make_shared<const native_gas_stage::Binding>(*out.gas,out.registry,blob,r.options);
  const auto predicted=predict(out.registry,r.packets,cap);
  const auto census_hash=h::canonical_census_sha256(predicted,cap);
  // Hall constructor performs the actual absolute Gauss/divB admission. No
  // seed field or particle charge is repaired to make an inconsistent seed fit.
  out.hall=std::make_shared<h::Owner>(r.model,r.seed,charge(r,predicted),r.options,predicted);
  require(out.hall->motion_mesh_sha256()==r.model->motion_mesh_sha256(r.options),"Hall/native geometry mismatch");
  std::set<std::pair<int,int>> slices;
  for(const auto&p:r.ports) {
    require(slices.emplace(p.slice.body,p.slice.slice).second,"duplicate source slice");
    out.sources.ports.push_back(ss::bind_port(out.hall,p.slice,p.options,r.options,
                                             r.native_configuration.chemistry,p.registration));
  }
  require(r.ports.size()==std::size_t(2*r.model->mesh().nz_foil),"one source port per foil slice required");
  out.sources.arrival_policy=r.arrival_policy;
  out.sources.arrival=std::make_shared<ss::pa::Owner>(r.arrival_policy,0);
  out.arrival_state_sha256=out.sources.arrival->state_sha256();
  Encoding config;config.text("joint-native-gas-source-Hall-genesis-v9/configuration");
  for(const auto*s:{&out.registry.configuration_sha256,&out.registry.registry_sha256,&out.registry.blob_sha256})config.text(*s);
  config.text(r.model->configuration_sha256());config.text(out.gas->configuration_sha256());
  config.text(r.arrival_policy->configuration_sha256());config.u64(cap);config.u64(r.ports.size());
  // Bind common-run resource policies as well as physical inputs. Native's
  // accepted-state hash intentionally omits history cadence, but common restart
  // admission must still retain the requested bounded-history policy.
  const auto &nc=r.native_configuration;
  config.u64(nc.limits.memory_budget_bytes);config.u64(nc.limits.maximum_identity_entries);
  config.u64(nc.history.anchor_interval_steps);config.u64(nc.history.maximum_anchors);
  config.u64(nc.history.maximum_anchor_bytes);config.u64(nc.history.maximum_window_rows);
  config.u64(nc.history.maximum_window_bytes);config.number(nc.motion.neutral_stream_dt_max_s);
  for(const auto&p:out.sources.ports) {
    config.text(p.source()->configuration_sha256());config.text(p.reference()->configuration_sha256());
    out.source_states.push_back(p.source()->state_sha256());
  }
  out.configuration_sha256=config.hash();
  Encoding source;source.text("joint-native-gas-source-Hall-genesis-v9/source");source.text(out.configuration_sha256);
  for(const auto&s:out.source_states)source.text(s);
  source.text(out.arrival_state_sha256);out.source_sha256=source.hash();
  n::FullTimeGenesisBinding b;b.time_s=r.seed.time_s;b.origin_provenance=r.origin_provenance;
  b.field_sha256=out.hall->field_sha256();b.heat_sha256=out.hall->heat_sha256();
  b.mesh_sha256=out.hall->motion_mesh_sha256();b.material_sha256=r.model->material_sha256();
  // Frame translation is independently retained in model geometry and native
  // options; bind all actual geometry operands, not a legacy species digest.
  Encoding frame;frame.text("joint-native-Hall-frame-v9");frame.text(b.mesh_sha256);
  frame.number(r.options.width_m);frame.number(r.options.gap_m);frame.number(r.options.depth_m);
  frame.number(r.options.short_z_m);frame.number(r.options.feed_z_m);b.frame_sha256=frame.hash();
  b.source_sha256=out.source_sha256;b.rng_sha256=out.registry.initial_rng_sha256;
  Encoding clock;clock.text("joint-native-gas-source-Hall-genesis-v9/clock");clock.number(b.time_s);
  const std::string*clock_items[]={&out.configuration_sha256,&b.field_sha256,&b.heat_sha256,&b.mesh_sha256,&b.material_sha256,
                                    &b.frame_sha256,&b.source_sha256,&b.rng_sha256,&census_hash};
  for(const auto*s:clock_items)clock.text(*s);
  clock.text(out.gas->accepted_sha256());b.common_clock_sha256=clock.hash();
  out.native=n::PicArena::create_full_time(r.options,n::full_time_chemistry_states(out.registry),
                                          {out.registry.channel},r.packets,b,r.native_configuration);
  const auto actual=out.native->rate_population(cap);
  require(h::canonical_census_sha256(actual,cap)==census_hash,"actual native census differs from Hall genesis");
  // No provisional source is replaced: every port retains this exact Hall
  // object's thermal view and lifetime throughout the accepted run.
  for(std::size_t i=0;i<out.sources.ports.size();++i)
    require(out.sources.ports[i].hall()==out.hall && out.sources.ports[i].source()->state_sha256()==out.source_states[i],
            "source changed during native genesis");
  require(out.sources.arrival->state_sha256()==out.arrival_state_sha256 && !out.gas->step_open(),
          "arrival/gas changed during native genesis");
  return out;
}
} // namespace powerflow::coupled_genesis_v9
