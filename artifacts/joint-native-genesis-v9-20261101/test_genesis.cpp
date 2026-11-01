#include "native_genesis.hpp"
#include "chem_v8_fixture.hpp"
#include "raikou/json.h"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace cg=powerflow::coupled_genesis_v9;
namespace ss=powerflow::source_stage_v8;
namespace h=hall_corrective_v3;
namespace n=raikou::mfem_backend;
namespace g=pf_gas_v1;
using J=raikou::Json;

cg::Request request(const char *trace) {
  cg::Request r;r.options=chemfix::options();
  auto &c=r.native_configuration;c.identity=n::FullTimeIdentitySchema::v8;
  c.chemistry=chemfix::chemistry_configuration(r.options);
  c.chemistry.diagnostic_field_target_approximation=false;
  c.limits.maximum_carriers=1000;c.limits.memory_budget_bytes=1000*n::full_time_bytes_per_carrier+n::full_time_fixed_bytes;
  c.limits.maximum_identity_entries=1u<<22;c.motion.neutral_stream_dt_max_s=1e-11;
  const auto reg=n::full_time_chemistry_registry(c.chemistry,r.options);
  const auto states=n::full_time_chemistry_states(reg);
  const auto water=std::find_if(states.begin(),states.end(),[](const auto&s){return s.id=="pfc-v3/H2O";});
  const auto spectrum=native_source_v7::electron_mass_configuration(surface_emission::Config{4.3},chemfix::blob_path(),chemfix::blob_sha(),reg,states);
  axial::Mesh geometry(45,4,.009,1);std::vector<ss::hb::SliceSource> slices;
  ss::pa::Configuration arrivals;arrivals.run_seed=c.chemistry.run_seed;
  for(int k=geometry.kshort;k<geometry.kfeed;++k)for(int body=0;body<2;++body) {
    const auto index=std::uint32_t(slices.size());cg::PortSpec port;auto &s=port.slice;
    s.body=body;s.slice=k;s.source_hash="manufactured/joint-genesis/finite-water-stock";
    ss::eh::Patch p;p.id="water/"+std::to_string(index);p.face=body?ss::eh::Face::LowY:ss::eh::Face::HighY;
    p.begin=0;p.end=1;p.state_id=p.parent_id="pfc-v3/H2O";p.mass_kg=water->mass_kg;p.initial_particles=64;
    p.desorption_enabled=true;p.attempt_s=1e13;p.activation_eV=p.desorption_eV=p.internal_eV=0;p.source_hash=s.source_hash;
    s.patches.push_back(p);ss::s7::SourceSettings settings;settings.spectrum=spectrum;
    settings.layer.local_linear_vacuum_extent_m=25e-6;settings.quantized_neutral_weight=2;
    settings.sampling_policy_identity="manufactured/joint-genesis/neutral-weight2";
    s.settings.push_back(settings);s.identities.push_back({p.id,"benchmark-water-layer",.005*geometry.dz});
    arrivals.streams.push_back({2*index,ss::pa::Kind::electron,index,"pfc-v3/e",1e6,64});
    arrivals.streams.push_back({2*index+1,ss::pa::Kind::neutral,index,"pfc-v3/H2O",2,64});
    port.registration.push_back({0,2*index+1,2*index,body==1,ss::hb::slice_center_z(geometry,k)});
    slices.push_back(s);r.ports.push_back(port);
  }
  h::MotionForms forms;forms.native_v7=forms.native_v81_bridges=true;forms.neutral_stream_dt_max_s=1e-11;
  r.model=std::make_shared<h::Model>(h::MeshSpec{},trace,h::make_cpu_pcg(),h::CarrierLimits{},ss::hb::declaration(slices),1e-14,forms);
  r.seed=r.model->dc_seed();r.arrival_policy=std::make_shared<const ss::pa::Policy>(arrivals);
  auto &gas=r.gas_configuration;gas.preset="manufactured/joint-genesis/host-census";
  gas.feedback=true;gas.transport=g::Transport::census_conservative;
  const auto grid=n::chem8::make_grid(r.options,c.chemistry.grid_counts);gas.grid={grid.low,grid.high,grid.n,grid.volume};
  gas.chemistry_blob_sha256=reg.blob_sha256;gas.chemistry_registry_sha256=reg.registry_sha256;
  for(const auto&s:reg.species)gas.species.push_back({s.state,s.mass_kg,s.charge_e,s.target});
  const auto blob=pfc::load_blob(chemfix::blob_path(),chemfix::blob_sha());
  for(std::uint32_t k=0;k<blob.n_proc;++k){const auto&p=blob.view.procs[k];gas.processes.push_back({blob.view.groups[p.group].projectile,p.n_out,native_gas_stage::internal_carrier(p.kind,p.n_out),p.kind});}
  r.origin_provenance="manufactured full-surface startup fixture at supported C19 drive start; no physical evolution";
  r.conductor_declaration="manufactured C19 DC-seed conductor charge";r.exterior_declaration="manufactured C19 DC-seed exterior charge";
  return r;
}
int main(int argc,char**argv) {
 J cases=J::array();int failed=0;
 auto check=[&](const std::string&name,bool pass,J evidence=J::object()){cases.push_back({{"name",name},{"passed",pass},{"evidence",evidence}});failed+=!pass;std::cout<<(pass?"PASS ":"FAIL ")<<name<<'\n';};
 try {
  if(argc!=3)throw std::runtime_error("TRACE OUT required");std::filesystem::create_directories(argv[2]);
  auto r=request(argv[1]);auto own=cg::admit(r);const auto b=own.native->full_time_boundary();
  check("native and Hall start at actual supported drive time",b&&b->time_s==own.hall->time_s()&&b->time_s==r.model->drive().start()&&b->accepted_step==0,{{"time_s",b->time_s}});
  check("all foil slices retain exact Hall and thermal owners",own.sources.ports.size()==std::size_t(2*r.model->mesh().nz_foil)&&std::all_of(own.sources.ports.begin(),own.sources.ports.end(),[&](const auto&p){return p.hall()==own.hall&&&p.source()->thermal_state()==&own.hall->thermal().slice(p.slice().body,p.slice().slice).state();}),{{"ports",own.sources.ports.size()}});
  check("actual catalogue and source roots bind native genesis",b->source_sha256==own.source_sha256&&b->field_sha256==own.hall->field_sha256()&&b->heat_sha256==own.hall->heat_sha256()&&own.gas_binding->registry().registry_sha256==own.registry.registry_sha256);
  check("gas starts accepted and empty without open attempt",!own.gas->step_open()&&own.gas->ledgers().entries.empty());
  auto refused=[&](const std::string&name,auto mutate){auto bad=r;mutate(bad);bool caught=false;std::string why;try{auto unused=cg::admit(bad);}catch(const std::exception&e){caught=true;why=e.what();}check(name,caught,{{"reason",why}});};
  refused("later cold start refused",[](auto&x){x.seed.time_s+=1e-14;});
  refused("diagnostic target-energy approximation refused",[](auto&x){x.native_configuration.chemistry.diagnostic_field_target_approximation=true;});
  refused("mismatched source collision seed refused",[](auto&x){auto p=x.arrival_policy->configuration();++p.run_seed;x.arrival_policy=std::make_shared<const ss::pa::Policy>(p);});
  refused("old leaking gas transport mode refused",[](auto&x){x.gas_configuration.transport=g::Transport::carrier_shares;});
  refused("unsupported gas device mode refused",[](auto&x){x.gas_configuration.execution=g::Execution::device;});
  refused("wrong process projectile despite matching hashes refused",[](auto&x){auto&p=x.gas_configuration.processes[0];p.projectile=(p.projectile+1)%x.gas_configuration.species.size();});
  refused("wrong process internal destination refused",[](auto&x){x.gas_configuration.processes[0].internal_carrier=1;});
  refused("wrong native gas cell volume refused",[](auto&x){x.gas_configuration.grid.volume=std::nextafter(x.gas_configuration.grid.volume,INFINITY);});
  refused("wrong exact native mass refused",[](auto&x){x.gas_configuration.species[0].mass_kg=std::nextafter(x.gas_configuration.species[0].mass_kg,INFINITY);});
  refused("missing surface slice refused",[](auto&x){x.ports.pop_back();});
  refused("duplicate surface slice refused",[](auto&x){x.ports[1]=x.ports[0];});
  refused("unsupported merge transport refused",[](auto&x){x.native_configuration.merge_split.enabled=true;});
  auto populated=r;n::FullTimeAuthoredPacket p;p.id="manufactured-neutral";p.origin="manufactured/joint-genesis/neutral";p.state="pfc-v3/H2O";p.count=n::ArenaCount{}.credited(2.5);p.position_m={.0045,.0045,0};p.proper_velocity_m_s={10,20,30};populated.packets.push_back(p);
  auto second=cg::admit(populated);const auto actual=second.native->rate_population(1000);const auto expected=cg::predict(second.registry,populated.packets,1000);
  check("populated actual native census equals Hall startup prediction",!actual.empty()&&h::canonical_census_sha256(actual,1000)==h::canonical_census_sha256(expected,1000),{{"carriers",actual.size()}});
  auto resource=r;++resource.native_configuration.limits.memory_budget_bytes;auto limited=cg::admit(resource);
  check("common configuration binds exact native resource budget",limited.configuration_sha256!=own.configuration_sha256);
  auto history=r;++history.native_configuration.history.anchor_interval_steps;auto archived=cg::admit(history);
  check("common configuration binds restart history policy",archived.configuration_sha256!=own.configuration_sha256);
  auto repeat=cg::admit(r);check("failed admissions leave deterministic retry roots unchanged",repeat.native->state_prefix()==own.native->state_prefix()&&repeat.source_sha256==own.source_sha256&&repeat.gas->accepted_sha256()==own.gas->accepted_sha256());
 }catch(const std::exception&e){check(std::string("unexpected: ")+e.what(),false);}
 std::ofstream(std::filesystem::path(argv[2])/"genesis.json")<<J{{"failed",failed},{"scope","manufactured startup admission only; no coupled evolution or common restart"},{"cases",cases}}.dump(1)<<'\n';return failed?1:0;
}
