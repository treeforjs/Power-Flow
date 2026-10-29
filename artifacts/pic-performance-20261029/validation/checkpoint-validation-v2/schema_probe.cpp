#include "raikou/mfem_gap_pic.h"
#include <cstddef>
#include <bit>
#include <iostream>
#include <type_traits>
using namespace raikou::mfem_backend;
int main(){
static_assert(sizeof(double)==8 && sizeof(int)==4 && sizeof(bool)==1);
static_assert(std::is_standard_layout_v<GapPicOptions> && std::is_standard_layout_v<GapPicDiagnostics> && std::is_standard_layout_v<GapParticle>);
std::cout<<"{\"endian\":\""<<(std::endian::native==std::endian::little?"little":"big")<<"\",\"options_size\":"<<sizeof(GapPicOptions)<<",\"diagnostics_size\":"<<sizeof(GapPicDiagnostics)<<",\"particle_size\":"<<sizeof(GapParticle)<<",";
std::cout<<"\"options\":[";
std::cout<<"{\"name\":\"nx\",\"kind\":\"i32\",\"offset\":"<<(offsetof(GapPicOptions,nx))<<"}";
std::cout<<",{\"name\":\"ny\",\"kind\":\"i32\",\"offset\":"<<(offsetof(GapPicOptions,ny))<<"}";
std::cout<<",{\"name\":\"width_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,width_m))<<"}";
std::cout<<",{\"name\":\"gap_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,gap_m))<<"}";
std::cout<<",{\"name\":\"depth_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,depth_m))<<"}";
std::cout<<",{\"name\":\"lower_potential_V\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,lower_potential_V))<<"}";
std::cout<<",{\"name\":\"upper_potential_V\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,upper_potential_V))<<"}";
std::cout<<",{\"name\":\"magnetic_field_T[0]\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,magnetic_field_T)+0*sizeof(double))<<"}";
std::cout<<",{\"name\":\"magnetic_field_T[1]\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,magnetic_field_T)+1*sizeof(double))<<"}";
std::cout<<",{\"name\":\"magnetic_field_T[2]\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,magnetic_field_T)+2*sizeof(double))<<"}";
std::cout<<",{\"name\":\"distributed_fields\",\"kind\":\"bool\",\"offset\":"<<(offsetof(GapPicOptions,distributed_fields))<<"}";
std::cout<<",{\"name\":\"finite_foils\",\"kind\":\"bool\",\"offset\":"<<(offsetof(GapPicOptions,finite_foils))<<"}";
std::cout<<",{\"name\":\"lower_foil.left_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,lower_foil)+offsetof(GapFoil,left_m))<<"}";
std::cout<<",{\"name\":\"lower_foil.right_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,lower_foil)+offsetof(GapFoil,right_m))<<"}";
std::cout<<",{\"name\":\"lower_foil.bottom_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,lower_foil)+offsetof(GapFoil,bottom_m))<<"}";
std::cout<<",{\"name\":\"lower_foil.top_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,lower_foil)+offsetof(GapFoil,top_m))<<"}";
std::cout<<",{\"name\":\"upper_foil.left_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,upper_foil)+offsetof(GapFoil,left_m))<<"}";
std::cout<<",{\"name\":\"upper_foil.right_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,upper_foil)+offsetof(GapFoil,right_m))<<"}";
std::cout<<",{\"name\":\"upper_foil.bottom_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,upper_foil)+offsetof(GapFoil,bottom_m))<<"}";
std::cout<<",{\"name\":\"upper_foil.top_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,upper_foil)+offsetof(GapFoil,top_m))<<"}";
std::cout<<",{\"name\":\"finite_axial_extent\",\"kind\":\"bool\",\"offset\":"<<(offsetof(GapPicOptions,finite_axial_extent))<<"}";
std::cout<<",{\"name\":\"short_z_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,short_z_m))<<"}";
std::cout<<",{\"name\":\"feed_z_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicOptions,feed_z_m))<<"}";
std::cout<<"],";
std::cout<<"\"diagnostics\":[";
std::cout<<"{\"name\":\"time_s\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,time_s))<<"}";
std::cout<<",{\"name\":\"particle_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,particle_charge_C))<<"}";
std::cout<<",{\"name\":\"deposited_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,deposited_charge_C))<<"}";
std::cout<<",{\"name\":\"lower_field_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,lower_field_charge_C))<<"}";
std::cout<<",{\"name\":\"upper_field_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,upper_field_charge_C))<<"}";
std::cout<<",{\"name\":\"gauss_defect_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,gauss_defect_C))<<"}";
std::cout<<",{\"name\":\"particle_charge_defect_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,particle_charge_defect_C))<<"}";
std::cout<<",{\"name\":\"field_energy_J\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,field_energy_J))<<"}";
std::cout<<",{\"name\":\"particle_kinetic_energy_J\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,particle_kinetic_energy_J))<<"}";
std::cout<<",{\"name\":\"lower_absorbed_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,lower_absorbed_charge_C))<<"}";
std::cout<<",{\"name\":\"upper_absorbed_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,upper_absorbed_charge_C))<<"}";
std::cout<<",{\"name\":\"lower_absorbed_energy_J\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,lower_absorbed_energy_J))<<"}";
std::cout<<",{\"name\":\"upper_absorbed_energy_J\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,upper_absorbed_energy_J))<<"}";
std::cout<<",{\"name\":\"lower_injected_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,lower_injected_charge_C))<<"}";
std::cout<<",{\"name\":\"upper_injected_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,upper_injected_charge_C))<<"}";
std::cout<<",{\"name\":\"injected_energy_J\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,injected_energy_J))<<"}";
std::cout<<",{\"name\":\"authored_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,authored_charge_C))<<"}";
std::cout<<",{\"name\":\"authored_kinetic_energy_J\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,authored_kinetic_energy_J))<<"}";
std::cout<<",{\"name\":\"maximum_nodal_continuity_defect_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,maximum_nodal_continuity_defect_C))<<"}";
std::cout<<",{\"name\":\"escaped_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,escaped_charge_C))<<"}";
std::cout<<",{\"name\":\"escaped_kinetic_energy_J\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,escaped_kinetic_energy_J))<<"}";
std::cout<<",{\"name\":\"short_absorbed_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,short_absorbed_charge_C))<<"}";
std::cout<<",{\"name\":\"short_absorbed_energy_J\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,short_absorbed_energy_J))<<"}";
std::cout<<",{\"name\":\"axial_escaped_charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,axial_escaped_charge_C))<<"}";
std::cout<<",{\"name\":\"poisson_residual_C_per_m\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapPicDiagnostics,poisson_residual_C_per_m))<<"}";
std::cout<<"],";
std::cout<<"\"particle\":[";
std::cout<<"{\"name\":\"id\",\"kind\":\"u64\",\"offset\":"<<(offsetof(GapParticle,id))<<"}";
std::cout<<",{\"name\":\"position_m[0]\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapParticle,position_m)+0*sizeof(double))<<"}";
std::cout<<",{\"name\":\"position_m[1]\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapParticle,position_m)+1*sizeof(double))<<"}";
std::cout<<",{\"name\":\"position_m[2]\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapParticle,position_m)+2*sizeof(double))<<"}";
std::cout<<",{\"name\":\"proper_velocity_m_s[0]\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapParticle,proper_velocity_m_s)+0*sizeof(double))<<"}";
std::cout<<",{\"name\":\"proper_velocity_m_s[1]\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapParticle,proper_velocity_m_s)+1*sizeof(double))<<"}";
std::cout<<",{\"name\":\"proper_velocity_m_s[2]\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapParticle,proper_velocity_m_s)+2*sizeof(double))<<"}";
std::cout<<",{\"name\":\"charge_C\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapParticle,charge_C))<<"}";
std::cout<<",{\"name\":\"mass_kg\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapParticle,mass_kg))<<"}";
std::cout<<",{\"name\":\"weight\",\"kind\":\"f64\",\"offset\":"<<(offsetof(GapParticle,weight))<<"}";
std::cout<<"]}";
}
