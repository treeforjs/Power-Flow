#pragma once
// Read-only operand reconstruction, never an emission acceptance API.
// The actual sealed begin_layer remains the sole source of each outcome.
#include "surface_emission.hpp"
#include "raikou/json.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace powerflow::magnetic_source_census_v1 {
using J=raikou::Json;
namespace se=surface_emission;
inline double dot(const se::Vector&a,const se::Vector&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline J spectrum(const se::PreparedFlux&f,const se::LayerPolicy&p){
 if(f.spectrum.nodes.size()>8192)throw std::runtime_error("bounded spectrum census exceeded");
 J nodes=J::array();for(const auto&n:f.spectrum.nodes)nodes.push_back({n.normal_energy_eV,n.probability,n.cdf,n.mean_transverse_eV});
 return {{"T_K",f.spectrum.T_K},{"stage_T_K",f.stage.T_K},{"F_V_m",f.spectrum.extraction_field_V_m},
 {"E_V_m",f.stage.total_E},{"B_T",f.stage.total_B_T},{"outward_normal",f.stage.outward_normal},{"area_m2",f.stage.area_m2},
 {"field_time_s",f.stage.field_time},{"thermal_time_s",f.stage.thermal_time},{"stage_time_s",f.stage.stage_time},
 {"field_accepted_state",f.stage.field_accepted_state},{"thermal_accepted_state",f.stage.thermal_accepted_state},
 {"accepted_state",f.stage.accepted_state},{"mesh_revision",f.stage.mesh_revision},{"material_revision",f.stage.material_revision},
 {"work_function_eV",f.spectrum.config.work_function_eV},{"electron_mass_kg",f.spectrum.config.electron_mass_kg},
 {"electron_mass_identity",f.spectrum.config.electron_mass_identity},{"current_A_m2",f.spectrum.current_A_m2},
 {"log_current_A_m2",f.spectrum.log_current_A_m2},{"barrier_lowering_eV",f.spectrum.barrier_lowering_eV},
 {"fermi_decay_width_eV",f.spectrum.fermi_decay_width_eV},{"energy_panels",f.spectrum.config.energy_panels},
 {"action_nodes",f.spectrum.config.action_nodes},{"tail_widths",f.spectrum.config.tail_widths},
 {"image_energy_fraction",p.image_energy_fraction},{"image_force_fraction",p.image_force_fraction},
 {"quantum_fraction_limit",p.maximum_quantum_magnetic_energy_fraction},{"local_extent_m",p.local_linear_vacuum_extent_m},
 {"node_columns",{"normal_energy_eV","probability","cdf","mean_transverse_eV"}},{"nodes",std::move(nodes)}};
}
inline J operands(const se::PreparedFlux&f,const se::LayerPolicy&p,double un,double ut,double ua){
 if(!(un>0&&un<1&&ut>0&&ut<1&&ua>0&&ua<1))throw std::invalid_argument("invalid preserved sample mark");
 if(f.spectrum.nodes.empty()||f.spectrum.nodes.size()>8192)throw std::invalid_argument("invalid bounded source spectrum");
 const auto it=std::lower_bound(f.spectrum.nodes.begin(),f.spectrum.nodes.end(),un,[](const auto&n,double u){return n.cdf<u;});
 if(it==f.spectrum.nodes.end())throw std::invalid_argument("sample outside source CDF");
 const double E=it->normal_energy_eV,t=se::transverse_quantile_eV(E,f.spectrum.T_K,ut),F=f.spectrum.extraction_field_V_m;
 J out={{"marks",{un,ut,ua}},{"node_index",std::size_t(it-f.spectrum.nodes.begin())},{"normal_metal_energy_eV",E},{"transverse_energy_eV",t},
 {"magnetic_barrier_checked",F>0},{"mode",F>0?"extracting-Schottky-overlap":"thermal-surface-step"}};
 if(F<=0)return out;
 constexpr double pi=3.1415926535897932384626433832795;
 constexpr double image_c=se::qe/(16*pi*se::eps0),hbar=se::planck/(2*pi);
 const double mass=f.spectrum.config.electron_mass_kg,phi=f.spectrum.config.work_function_eV;
 const double scale=std::max(se::kb*f.spectrum.T_K/se::qe,f.spectrum.fermi_decay_width_eV),overlap=p.image_energy_fraction*scale;
 const double A=phi-E+overlap,lowering=f.spectrum.barrier_lowering_eV;
 const double x=A>lowering?(A+std::sqrt((A-lowering)*(A+lowering)))/(2*F):std::sqrt(image_c/F);
 const auto&n=f.stage.outward_normal;const double Bn=dot(f.stage.total_B_T,n);
 se::Vector Bt=f.stage.total_B_T;for(int k=0;k<3;++k)Bt[k]-=n[k]*Bn;
 const double bt=std::sqrt(dot(Bt,Bt)),vt=std::sqrt(2*se::qe*t/mass);
 const double linear=bt*x*vt,quadratic=se::qe*dot(Bt,Bt)*x*x/(2*mass),normal=.5*hbar*std::abs(Bn)/mass;
 const double fraction=(linear+quadratic+normal)/scale;
 out.update({{"source_scale_eV",scale},{"overlap_allowance_eV",overlap},{"barrier_root_A_eV",A},{"overlap_distance_m",x},
 {"barrier_maximum_distance_m",std::sqrt(image_c/F)},{"Bnormal_T",Bn},{"Btangent_T",bt},{"transverse_speed_m_s",vt},
 {"linear_bound_eV",linear},{"quadratic_eV",quadratic},{"normal_zero_point_eV",normal},{"magnetic_fraction",fraction},
 {"quantum_fraction_limit",p.maximum_quantum_magnetic_energy_fraction},{"magnetic_guard_would_refuse",fraction>p.maximum_quantum_magnetic_energy_fraction},
 {"geometry_guard_would_refuse",!(x>0&&x<p.local_linear_vacuum_extent_m)},
 {"handoff_distance_m",std::max(std::sqrt(image_c/(p.image_force_fraction*F)),image_c/(p.image_energy_fraction*scale))}});
 return out;
}
inline void compare_begin(const J&d,const se::LayerState&s){
 if(d.at("magnetic_barrier_checked").get<bool>()&&d.at("magnetic_fraction").get<double>()!=s.quantum_magnetic_energy_fraction)
  throw std::runtime_error("diagnostic magnetic reconstruction differs from sealed source outcome");
}
} // namespace powerflow::magnetic_source_census_v1
