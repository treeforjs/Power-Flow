#pragma once
// Boundary adapter only. Source core does not depend on native chemistry types.
#include "surface_emission.hpp"
#include "chem_host.hpp"
#include <algorithm>
#include <bit>
#include <set>
#include <string>
namespace native_source_v7 {
inline constexpr char supported_blob_sha256[]="70fcc2adc2fdaa62abed7ed02d68ffc09b1b18554eda8fa99c5451120ae454b8";
inline bool same(double a,double b){return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);}
inline bool digest(const std::string&s){return s.size()==64&&std::all_of(s.begin(),s.end(),[](char c){return (c>='0'&&c<='9')||(c>='a'&&c<='f');});}
// Registry and States are actual native full_time_chemistry_registry/states
// results. Templates avoid importing a second native header version into the
// source owner. Verify the blob independently and every registry mass/nucleus
// operand; callers cannot substitute a remembered rounded electron mass.
template<class Registry,class States>
surface_emission::Config electron_mass_configuration(surface_emission::Config cfg,const std::string&blob_path,
 const std::string&expected_blob_sha256,const Registry&reg,const States&states){
 using pfc::require;
 require(expected_blob_sha256==supported_blob_sha256&&reg.blob_sha256==expected_blob_sha256,"source mass: unsupported or mismatched pinned blob");
 require(digest(reg.registry_sha256)&&digest(reg.configuration_sha256)&&digest(reg.policy_sha256)&&digest(reg.grid_sha256),"source mass: missing native configuration provenance");
 const auto b=pfc::load_blob(blob_path,expected_blob_sha256);
 require(b.n_species==22&&reg.species.size()==b.n_species&&states.size()==b.n_species,"source mass: incomplete native catalogue");
 pfc::Sha256 registry;registry.text("native-v8/chemistry-v3-registry-v1");registry.text(b.sha256);registry.text("nuclei-abi-v4");registry.pod(std::uint64_t(b.n_species));
 std::set<std::string> names,state_names;double mass=0;bool electron=false;
 for(const auto&s:states)require(state_names.insert(s.id).second,"source mass: duplicate arena state");
 for(std::size_t i=0;i<b.n_species;++i){
  const auto&sp=b.view.species[i];const std::string name="pfc-v3/"+b.species_names[i];
  const auto it=std::find_if(reg.species.begin(),reg.species.end(),[&](const auto&s){return s.state==name;});
  require(it!=reg.species.end()&&names.insert(name).second,"source mass: missing or duplicate native species");const auto&r=*it;
  const double m=(sp.mass_eV*surface_emission::qe)/(surface_emission::c_light*surface_emission::c_light);
  require(same(r.mass_kg,m)&&same(r.mass_eV,sp.mass_eV)&&same(r.reference_eV,sp.ref_eV)&&r.charge_e==sp.charge,"source mass: native species operand differs from pinned blob");
  registry.text(name);registry.pod(std::bit_cast<std::uint64_t>(m));registry.pod(std::uint64_t(std::int64_t(sp.charge)));
  require(sp.nuclei[5]==0,"source mass: nonzero spare blob nucleus slot");
  for(std::size_t k=0;k<r.nuclei.size();++k){const auto n=k<5?std::int32_t(sp.nuclei[k]):0;require(r.nuclei[k]==n,"source mass: native nucleus operand differs from blob");registry.pod(std::uint32_t(n));}
  const auto a=std::find_if(states.begin(),states.end(),[&](const auto&s){return s.id==name;});
  require(a!=states.end()&&same(a->mass_kg,m)&&a->charge_e==sp.charge&&a->hydrogen==int(sp.nuclei[0])&&a->oxygen==int(sp.nuclei[1])&&same(a->internal_J,0.0)&&same(a->reference_J,sp.ref_eV*surface_emission::qe),"source mass: arena mass/charge/nucleus/energy convention differs from blob");
  if(name=="pfc-v3/e"){require(!electron&&sp.charge==-1&&same(sp.ref_eV,0.0),"source mass: electron charge/reference convention");for(auto n:r.nuclei)require(n==0,"source mass: electron has nuclei");mass=m;electron=true;}
 }
 require(electron&&registry.hex()==reg.registry_sha256,"source mass: native registry digest mismatch");
 pfc::Sha256 policy;policy.text("source-v7.1/blob-electron-mass-v1");policy.text(b.sha256);policy.text(reg.registry_sha256);policy.text(reg.configuration_sha256);policy.text(reg.policy_sha256);policy.text(reg.grid_sha256);policy.text("pfc-v3/e");policy.pod(std::bit_cast<std::uint64_t>(mass));
 cfg.electron_mass_kg=mass;cfg.electron_mass_identity="source-v7.1/blob-electron-mass-v1/"+policy.hex();surface_emission::validate_electron_mass(cfg);return cfg;
}
}
