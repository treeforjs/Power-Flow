#include "source_hall_iteration.hpp"
#include "positive_source_fixture.hpp"
#include "final_charge_guard.hpp"
#include "chem_v8_fixture.hpp"
#include "raikou/json.h"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace ss=powerflow::source_stage_v8;namespace h=hall_corrective_v3;namespace n=raikou::mfem_backend;
namespace it=powerflow::source_hall_iteration_v1;namespace fx=powerflow::positive_source_fixture_v1;namespace cg=powerflow::coupled_genesis_v9;
using J=raikou::Json;
struct Fixture {
 n::GapPicOptions options=chemfix::options();n::FullTimeChemistryConfiguration chemistry=chemfix::chemistry_configuration(options);
 n::FullTimeChemistryRegistry registry=n::full_time_chemistry_registry(chemistry,options);
 std::vector<n::ArenaStateSpec> states=n::full_time_chemistry_states(registry);
 std::shared_ptr<h::Model> model;std::shared_ptr<h::Owner>hall;ss::OwnerSet owners;ss::AcceptedSourceState base;
 std::unique_ptr<n::PicArena>arena;n::FullTimeStepBinding common;fx::Receipt declared;
 std::string native_before,field_before,heat_before;double dt=1e-14;
 Fixture(const char*trace,bool water_enabled,double weight){
  axial::Mesh geometry(45,4,.009,1);std::vector<cg::PortSpec>ports;ss::pa::Configuration arrival;arrival.run_seed=chemistry.run_seed;
  const auto water=std::find_if(states.begin(),states.end(),[](const auto&s){return s.id=="pfc-v3/H2O";});
  const auto config=native_source_v7::electron_mass_configuration(surface_emission::Config{4.3},chemfix::blob_path(),chemfix::blob_sha(),registry,states);
  for(int k=geometry.kshort;k<geometry.kfeed;++k)for(int body=0;body<2;++body){
   const auto index=std::uint32_t(ports.size());cg::PortSpec port;auto&s=port.slice;s.body=body;s.slice=k;s.source_hash="manufactured/positive-source-Hall/C19-start";
   ss::eh::Patch p;p.id="gap/"+std::to_string(index);p.face=body?ss::eh::Face::LowY:ss::eh::Face::HighY;p.begin=0;p.end=1;p.state_id=p.parent_id="pfc-v3/H2O";p.mass_kg=water->mass_kg;p.initial_particles=64;p.desorption_enabled=true;p.attempt_s=1e13;p.activation_eV=p.desorption_eV=p.internal_eV=0;p.source_hash=s.source_hash;s.patches.push_back(p);
   ss::s7::SourceSettings set;set.spectrum=config;set.layer.local_linear_vacuum_extent_m=25e-6;set.quantized_neutral_weight=2;set.sampling_policy_identity="manufactured/positive-source-Hall";s.settings.push_back(set);s.identities.push_back({p.id,"manufactured-layer",.005*geometry.dz});
   arrival.streams.push_back({2*index,ss::pa::Kind::electron,index,"pfc-v3/e",1e6,64});arrival.streams.push_back({2*index+1,ss::pa::Kind::neutral,index,"pfc-v3/H2O",2,64});port.registration.push_back({0,2*index+1,2*index,body==1,ss::hb::slice_center_z(geometry,k)});ports.push_back(port);
  }
  fx::Spec spec;spec.water_enabled=water_enabled;spec.electron_weight=weight;declared=fx::apply(ports,arrival,options,chemistry,spec);
  std::vector<ss::hb::SliceSource>slices;for(const auto&p:ports)slices.push_back(p.slice);h::MotionForms forms;forms.native_v7=forms.native_v81_bridges=true;forms.neutral_stream_dt_max_s=1e-11;
  model=std::make_shared<h::Model>(h::MeshSpec{},trace,h::make_cpu_pcg(),h::CarrierLimits{},ss::hb::declaration(slices),dt,forms);auto seed=model->dc_seed();const auto&g=model->mesh();h::ChargeLedger charge;charge.particles=h::charge_of_census(g,options,{});charge.conductor_C.resize(g.nn);charge.exterior_C.resize(g.nn);charge.conductor_declaration=chemfix::sha("positive/actual-dc-conductor");charge.exterior_declaration=chemfix::sha("positive/actual-dc-exterior");const auto Q=h::field_charge(g,seed);for(int k=0;k<g.nn;++k){if(g.node_body[k]>=0||g.short_node[k])charge.conductor_C[k]=Q[k];else if(g.outer_node[k])charge.exterior_C[k]=Q[k];}
  hall=std::make_shared<h::Owner>(model,seed,charge,options,std::span<const n::ArenaRateEntry>{});owners.arrival_policy=std::make_shared<const ss::pa::Policy>(arrival);owners.arrival=std::make_shared<ss::pa::Owner>(owners.arrival_policy,0);for(auto&p:ports)owners.ports.push_back(ss::bind_port(hall,p.slice,p.options,options,chemistry,p.registration));
  base.arrival_state=owners.arrival->state_sha256();h::Digest digest("manufactured/positive-source-Hall/genesis");for(const auto&p:owners.ports){base.source_states.push_back(p.source()->state_sha256());digest.text(base.source_states.back());}digest.text(base.arrival_state);base.source_seal=digest.hex();
  n::FullTimeConfiguration nc;nc.identity=n::FullTimeIdentitySchema::v8;nc.chemistry=chemistry;nc.limits.maximum_carriers=1000;nc.limits.memory_budget_bytes=1000*n::full_time_bytes_per_carrier+n::full_time_fixed_bytes;nc.limits.maximum_identity_entries=1u<<22;nc.motion.neutral_stream_dt_max_s=1e-11;
  n::FullTimeGenesisBinding genesis;genesis.time_s=hall->time_s();genesis.field_sha256=hall->field_sha256();genesis.heat_sha256=hall->heat_sha256();genesis.source_sha256=base.source_seal;genesis.mesh_sha256=hall->motion_mesh_sha256();genesis.material_sha256=model->material_sha256();genesis.frame_sha256=chemfix::sha("positive/frame");genesis.common_clock_sha256=chemfix::sha("positive/genesis-clock");genesis.rng_sha256=registry.initial_rng_sha256;genesis.origin_provenance="manufactured phi2eV cold C19-start source-Hall component; diagnostic chemistry, not plasma calibration";
  arena=n::PicArena::create_full_time(options,states,{registry.channel},{},genesis,nc);const auto b=arena->full_time_boundary();common={b->accepted_step,b->accepted_arena_revision,b->time_s,b->time_s+dt,dt,b->accepted_arena_sha256,b->field_sha256,b->mesh_sha256,b->material_sha256,chemfix::sha("positive/fixed-common-clock")};native_before=arena->state_prefix();field_before=hall->field_sha256();heat_before=hall->heat_sha256();
 }
 bool unchanged()const {if(arena->state_prefix()!=native_before||hall->field_sha256()!=field_before||hall->heat_sha256()!=heat_before||owners.arrival->state_sha256()!=base.arrival_state)return false;for(std::size_t i=0;i<owners.ports.size();++i)if(owners.ports[i].source()->state_sha256()!=base.source_states[i])return false;return true;}
};
