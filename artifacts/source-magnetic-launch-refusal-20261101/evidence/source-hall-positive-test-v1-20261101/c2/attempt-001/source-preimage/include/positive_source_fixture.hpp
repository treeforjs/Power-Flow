#pragma once
// Manufactured qualification inputs only. Never a default SS304 material law.
#include "native_genesis.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <map>
#include <set>
namespace powerflow::positive_source_fixture_v1 {
namespace cg=::powerflow::coupled_genesis_v9;
namespace ss=::powerflow::source_stage_v8;
namespace n=::raikou::mfem_backend;
struct Spec {
 double work_function_eV=.75;
 double electron_weight=.01;
 bool water_enabled=false;
 // Explicit finite stock/packet inputs for the simultaneous-source variant.
 double water_stock_per_patch=64.,water_packet_weight=2.;
 double water_attempt_s=1e13,water_activation_eV=0.,water_desorption_eV=0.;
 std::uint32_t maximum_electrons_per_stream=32,maximum_electrons_total=128;
};
struct Receipt {
 Spec declared;
 std::uint64_t electron_weight_bits=0,work_function_bits=0,run_seed=0;
 std::string registry_sha256,blob_sha256,electron_mass_identity;
 double electron_mass_kg=0,water_mass_kg=0;
 std::size_t lower_patches=0,upper_patches=0;
 const char* scope="manufactured300K_C19-start_source-Hall_component_not_SS304_emission_calibration";
};
// Call before constructing the Hall Model/thermal declaration or any owner.
// Rebuild the model from hb::declaration(the modified port slices), use its real
// dc_seed at drive.start(), and let actual genesis enforce Gauss and identities.
inline Receipt apply(std::vector<cg::PortSpec>&ports,ss::pa::Configuration&arrival,
                     const n::GapPicOptions&options,const n::FullTimeChemistryConfiguration&chemistry,
                     const Spec&spec={}) {
 auto need=[](bool ok,const char*why){if(!ok)throw std::invalid_argument(why);};
 need(!ports.empty()&&std::isfinite(spec.work_function_eV)&&spec.work_function_eV>0&&
      std::isnormal(spec.electron_weight)&&spec.electron_weight>0,"positive fixture needs explicit finite work function and macro weight");
 need(spec.maximum_electrons_per_stream>0&&spec.maximum_electrons_per_stream<=32&&
      spec.maximum_electrons_total>0&&spec.maximum_electrons_total<=128,"positive fixture must retain bounded birth caps");
 need(spec.water_stock_per_patch==64&&spec.water_packet_weight==2&&spec.water_attempt_s==1e13&&
      spec.water_activation_eV==0&&spec.water_desorption_eV==0,"water component fixture has fixed declared finite inputs");
 const auto registry=n::full_time_chemistry_registry(chemistry,options);const auto states=n::full_time_chemistry_states(registry);
 const auto water=std::find_if(states.begin(),states.end(),[](const auto&s){return s.id=="pfc-v3/H2O";});
 need(water!=states.end()&&water->charge_e==0&&water->internal_J==0,"actual ground-state registry water required");
 ss::se::Config spectrum{spec.work_function_eV};spectrum=native_source_v7::electron_mass_configuration(spectrum,chemistry.blob_path,chemistry.blob_sha256,registry,states);
 Receipt receipt;receipt.declared=spec;receipt.electron_weight_bits=std::bit_cast<std::uint64_t>(spec.electron_weight);receipt.work_function_bits=std::bit_cast<std::uint64_t>(spec.work_function_eV);receipt.run_seed=chemistry.run_seed;receipt.registry_sha256=registry.registry_sha256;receipt.blob_sha256=registry.blob_sha256;receipt.electron_mass_identity=spectrum.electron_mass_identity;receipt.electron_mass_kg=spectrum.electron_mass_kg;receipt.water_mass_kg=water->mass_kg;
 std::set<std::uint32_t>keep;
 for(auto&port:ports){
  need(!port.options.phase&&port.options.water_writeoff_T_K==0,"fixture cannot remove enabled phase/writeoff physics");
  auto&slice=port.slice;need(slice.body==0||slice.body==1,"fixture electrode body");
  need(slice.patches.size()==slice.settings.size()&&slice.patches.size()==port.registration.size(),"fixture requires complete existing patch registration");
  for(auto&reg:port.registration){
   need(reg.patch<slice.patches.size()&&reg.electron_stream.has_value(),"fixture requires each gap-face electron stream");
   auto&p=slice.patches[reg.patch];auto&setting=slice.settings[reg.patch];
   need(p.face==(slice.body?ss::eh::Face::LowY:ss::eh::Face::HighY),"positive fixture accepts actual gap-facing foil patches only");
   need(setting.surface.benchmark()&&setting.surface.enhancement.classes.size()==1,"fixture does not override enabled enhancement/wall/photo models");
   setting.spectrum=spectrum;setting.sampling_policy_identity="manufactured/positive-source-B3/explicit-phi-and-weight";
   p.state_id=p.parent_id="pfc-v3/H2O";p.mass_kg=water->mass_kg;p.internal_eV=0;p.coverage_temkin=false;
   p.desorption_enabled=spec.water_enabled;p.initial_particles=spec.water_enabled?spec.water_stock_per_patch:0.;p.initial_low_particles=0.;p.attempt_s=spec.water_attempt_s;p.activation_eV=spec.water_activation_eV;p.desorption_eV=spec.water_desorption_eV;
   setting.quantized_neutral_weight=spec.water_packet_weight;
   if(slice.body)++receipt.upper_patches;else ++receipt.lower_patches;
   const auto electron=*reg.electron_stream;need(keep.insert(electron).second,"duplicate fixture electron stream");
   auto es=std::find_if(arrival.streams.begin(),arrival.streams.end(),[&](const auto&s){return s.stream_id==electron;});need(es!=arrival.streams.end()&&es->kind==ss::pa::Kind::electron&&es->species_state=="pfc-v3/e","fixture electron stream catalogue");es->target_weight=spec.electron_weight;es->maximum_births=spec.maximum_electrons_per_stream;
   if(spec.water_enabled){need(reg.neutral_stream.has_value(),"water fixture needs predeclared neutral stream");const auto neutral=*reg.neutral_stream;need(keep.insert(neutral).second,"duplicate water fixture stream");auto ns=std::find_if(arrival.streams.begin(),arrival.streams.end(),[&](const auto&s){return s.stream_id==neutral;});need(ns!=arrival.streams.end()&&ns->kind==ss::pa::Kind::neutral&&ns->species_state=="pfc-v3/H2O","fixture neutral stream catalogue");ns->target_weight=spec.water_packet_weight;}
   else reg.neutral_stream.reset();
  }
 }
 need(receipt.lower_patches&&receipt.upper_patches,"fixture requires both electrodes");
 arrival.streams.erase(std::remove_if(arrival.streams.begin(),arrival.streams.end(),[&](const auto&s){return !keep.contains(s.stream_id);}),arrival.streams.end());
 arrival.run_seed=chemistry.run_seed;arrival.maximum_electron_births_per_step=spec.maximum_electrons_total;
 // Constructor validates exact positive binary64 weights and all original
 // policy/native caps. No accepted state, field, charge or temperature changed.
 (void)ss::pa::Policy(arrival);return receipt;
}
} // namespace powerflow::positive_source_fixture_v1
