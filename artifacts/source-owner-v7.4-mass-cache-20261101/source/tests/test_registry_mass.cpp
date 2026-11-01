#include "chem_v8_fixture.hpp"
#include "native_source_v7_mass.hpp"
#include "fast_emission.hpp"
#include "emission_hot_v7.hpp"
#include "raikou/json.h"
#include <filesystem>
#include <fstream>
#include <iostream>
namespace se=surface_emission;namespace m=raikou::mfem_backend;using J=raikou::Json;
int main(int argc,char**argv){int fail=0;J cases=J::array();auto check=[&](std::string n,bool b){cases.push_back({{"name",n},{"group","registry-mass"},{"passed",b}});if(!b)++fail;std::cout<<n<<": "<<b<<'\n';};auto refuses=[&](auto f){try{f();return false;}catch(const std::exception&){return true;}};
 try{if(argc!=2)throw std::runtime_error("output directory required");std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
 auto opt=chemfix::options();auto cc=chemfix::chemistry_configuration(opt);auto reg=m::full_time_chemistry_registry(cc,opt);auto states=m::full_time_chemistry_states(reg);
 const auto config=native_source_v7::electron_mass_configuration(se::Config{4.3},chemfix::blob_path(),chemfix::blob_sha(),reg,states);
 const auto electron=std::find_if(reg.species.begin(),reg.species.end(),[](auto&s){return s.state=="pfc-v3/e";});
 check("actual pinned registry mass chosen by exact state",electron!=reg.species.end()&&native_source_v7::same(config.electron_mass_kg,electron->mass_kg)&&config.electron_mass_kg!=se::me&&!se::legacy_electron_mass(config));
 auto permuted=reg;std::reverse(permuted.species.begin(),permuted.species.end());auto ps=states;std::reverse(ps.begin(),ps.end());auto pc=native_source_v7::electron_mass_configuration(se::Config{4.3},chemfix::blob_path(),chemfix::blob_sha(),permuted,ps);check("lookup independent of input ordering",pc.electron_mass_kg==config.electron_mass_kg&&pc.electron_mass_identity==config.electron_mass_identity);
 for(int k=0;k<9;++k){auto r=reg;auto s=states;std::string sha=chemfix::blob_sha();
  auto e=std::find_if(r.species.begin(),r.species.end(),[](auto&x){return x.state=="pfc-v3/e";});auto a=std::find_if(s.begin(),s.end(),[](auto&x){return x.id=="pfc-v3/e";});
  if(k==0)e->mass_kg=se::me;if(k==1)e->charge_e=0;if(k==2)e->nuclei[15]=1;if(k==3)a->internal_J=se::qe;if(k==4)a->reference_J=se::qe;if(k==5)r.registry_sha256[0]=r.registry_sha256[0]=='0'?'1':'0';if(k==6)sha[0]='0';if(k==7)r.configuration_sha256.clear();if(k==8)r.species.pop_back();
  check("altered registry refusal "+std::to_string(k),refuses([&]{native_source_v7::electron_mass_configuration(se::Config{4.3},chemfix::blob_path(),sha,r,s);}));}
 se::Stage stage;stage.accepted_state=stage.field_accepted_state=stage.thermal_accepted_state=1;stage.mesh_revision=stage.material_revision=stage.patch_id=1;stage.dt=1e-14;stage.outward_normal={0,1,0};stage.area_m2=1e-6;stage.material_temperature_phase_admitted=stage.locally_planar_barrier_admitted=stage.local_equilibrium_electron_supply_admitted=true;
 se::LayerPolicy lp;lp.local_linear_vacuum_extent_m=25e-6;
 for(double T:{300.,1500.,3000.})for(double F:{-1e8,0.,1e8}){stage.T_K=T;stage.total_E={0,-F,0};auto f=se::prepare(stage,config);auto old=se::prepare(stage,se::Config{4.3});auto p=se::sample(f,{F>0?se::PlaneKind::ExplicitSchottkyNordheimHandoff:se::PlaneKind::ThermalVacuumSurfaceStep,F>0?se::required_SN_matching_distance_m(f.spectrum)*1.1:0,true,true,true},1,.6,.35,.1,1);double K=se::kinetic_energy_J(p.proper_velocity_m_s,config.electron_mass_kg);check("birth uses registry mass and exact energy "+std::to_string(T)+"/"+std::to_string(F),native_source_v7::same(p.electron_mass_kg,config.electron_mass_kg)&&std::abs(K-p.ledger.outgoing_K_J)<=2e-12*std::abs(K));
  auto l=se::begin_layer(f,lp,1,.6,.35,.1,2);check("layer retains declared mass "+std::to_string(T)+"/"+std::to_string(F),native_source_v7::same(l.state.particle.electron_mass_kg,config.electron_mass_kg));
  if(l.state.status==se::LayerStatus::Active){auto a=se::advance_layer(l.state,stage,lp,1e-16);check("active layer mass and work closure "+std::to_string(T),native_source_v7::same(a.state.particle.electron_mass_kg,config.electron_mass_kg)&&std::abs(a.exchange.residual_J())<=2e-10*(std::abs(a.exchange.particle_K_change_J)+std::abs(a.exchange.field_work_J)+std::abs(a.exchange.image_energy_change_J)));}
  auto fast=electrode_heat_v6_1::emission::evaluate(T,F,config);auto normal=electrode_heat_v6_1::emission::evaluate(T,F,se::Config{4.3});check("nonlegacy cache isolation "+std::to_string(T)+"/"+std::to_string(F),fast.config.electron_mass_identity==config.electron_mass_identity&&normal.config.electron_mass_identity.empty()&&native_source_v7::same(fast.current_A_m2,f.spectrum.current_A_m2));
 }
 se::Config twice{4.3};twice.electron_mass_kg=2*se::me;twice.electron_mass_identity="manufactured-mass-scaling-test";auto doubled=se::evaluate(1500,0,twice),single=se::evaluate(1500,0,se::Config{4.3});check("manufactured thermal supply scales with declared mass",std::abs(doubled.current_A_m2/single.current_A_m2-2)<2e-12);check("manufactured tunneling action scales with square root mass",std::abs(se::action(-1,1e8,twice)/se::action(-1,1e8,se::Config{4.3})-std::sqrt(2.))<2e-15);
 auto invalid=config;invalid.electron_mass_kg=0;check("zero mass refused",refuses([&]{se::evaluate(1500,0,invalid);}));invalid=config;invalid.electron_mass_identity.clear();check("unidentified nonlegacy mass refused",refuses([&]{se::evaluate(1500,0,invalid);}));
 auto before_stats=electrode_heat_v6_1::emission::statistics();
 auto cached=electrode_heat_v6_1::emission::evaluate(1500,3e9,config);
 auto after_stats=electrode_heat_v6_1::emission::statistics();
 check("actual registry mass takes qualified extracting fast path",after_stats.extracting_fast>before_stats.extracting_fast&&electrode_heat_v6_1::emission::same_distribution(cached,se::evaluate(1500,3e9,config)));
 auto hot=electrode_heat_v7::emission_hot::evaluate(3500,1e8,config,4000);check("hot law retains declared mass",hot.config.electron_mass_kg==config.electron_mass_kg&&hot.current_A_m2>0);
 std::ofstream(out/"mass-configuration.json")<<J{{"mass_kg",config.electron_mass_kg},{"identity",config.electron_mass_identity},{"blob_sha256",reg.blob_sha256},{"registry_sha256",reg.registry_sha256},{"configuration_sha256",reg.configuration_sha256}}.dump(2)<<'\n';
 }catch(const std::exception&e){check(std::string("unexpected: ")+e.what(),false);}
 std::ofstream(std::filesystem::path(argv[1])/"registry-mass.json")<<J{{"failed",fail},{"cases",cases}}.dump(1)<<'\n';return fail?1:0;}
