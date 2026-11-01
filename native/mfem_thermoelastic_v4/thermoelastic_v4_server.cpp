#include "thermoelastic_v4.h"
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
namespace {
void require(bool b,const std::string& s){if(!b)throw std::invalid_argument(s);}
void call(int code){if(code)throw std::runtime_error(tfv4_last_error());}
template<class T>void read(std::istream& in,T& x){
 if constexpr(std::is_unsigned_v<T>){
  std::string token;require(bool(in>>token)&&!token.empty()&&token.front()!='-',"unsigned field must not be negative");
  std::istringstream number(token);require(bool(number>>x)&&number.peek()==std::char_traits<char>::eof(),"malformed unsigned field");
 }else require(bool(in>>x),"malformed numeric field");
}
std::string quoted(std::istream& in){std::string s;require(bool(in>>std::quoted(s)),"malformed quoted field");return s;}
void complete(std::istream& in){in>>std::ws;require(in.eof(),"unexpected trailing fields");}
std::size_t count(std::istream& in){std::size_t n;read(in,n);require(n<=10000000,"count too large");return n;}
template<class T>void read_array(std::istream& in,T& a){for(auto& v:a)read(in,v);}
std::vector<double> values(std::istream& in,std::size_t n){std::vector<double> v(n);read_array(in,v);return v;}
template<class T>void output(std::ostream& o,const T& a){for(auto x:a)o<<' '<<x;}
tfv4_state_info info(tfv4_handle h,int trial=0){tfv4_state_info i{};call(tfv4_get_state_info(h,trial,&i));return i;}
void budget(std::ostream& o,const tfv4_budget& b){
 for(double x:{b.thermal_delta_J,b.mechanical_exchange_J,b.mechanical_absolute_activity_J,b.joule_J,b.authored_heat_J,b.wall_net_J,b.residual_J,b.roundoff_bound_J,b.accepted_reference_energy_J,b.trial_reference_energy_J,b.nodal_residual_max_abs_J,b.mechanical_load_conservation_residual_J})o<<' '<<x;o<<' '<<b.newton_iterations;
 for(double x:{b.source_wall_J,b.source_particle_K_change_J,b.source_electron_reference_change_J,b.source_image_energy_change_J,b.source_field_work_J,b.neutral_K_J,b.neutral_internal_J,b.adsorbed_reference_change_J,b.gas_reference_J,b.particle_charge_change_C,b.electrode_charge_change_C,b.physical_electron_number_change})o<<' '<<x;o<<' '<<b.new_gas_packet_count<<' '<<b.new_path_count;
}
void action(std::istream& in,tfv4_source_action& a){
 read(in,a.kind);read(in,a.layer_id);read(in,a.patch_index);read(in,a.patch_id);read(in,a.weight);read(in,a.normal_uniform);read(in,a.transverse_uniform);read(in,a.azimuth_uniform);read_array(in,a.origin_m);read(in,a.duration_s);read(in,a.work_function_eV);read(in,a.energy_panels);read(in,a.action_nodes);read(in,a.tail_widths);read(in,a.local_linear_vacuum_extent_m);read(in,a.image_force_fraction);read(in,a.image_energy_fraction);read(in,a.maximum_quantum_magnetic_energy_fraction);read(in,a.maximum_gyroangle);read(in,a.maximum_relative_normal_step);read(in,a.maximum_substeps);read(in,a.PIC_applies_microscopic_image_force);read(in,a.own_return_unit_accommodation);read(in,a.incoming_weighted_K_J);read(in,a.collection_accommodation_admitted);
}
}
int main(){tfv4_handle h=nullptr;std::string line;while(std::getline(std::cin,line)){bool exit=false;try{
 std::istringstream in(line);std::string cmd;read(in,cmd);std::ostringstream out;out<<std::setprecision(17)<<"OK";
 if(cmd=="CREATE"){
  require(!h,"CREATE requires fresh server");tfv4_config c{};read(in,c.abi_version);read(in,c.nx);read(in,c.ny);read(in,c.x0_m);read(in,c.y0_m);read(in,c.width_m);read(in,c.thickness_m);read(in,c.height_m);read(in,c.reference_density_kg_m3);read(in,c.initial_temperature_K);read(in,c.initial_time_s);auto mesh=quoted(in),material=quoted(in),source=quoted(in),electrode=quoted(in);auto n=count(in);require(c.nx>0&&c.ny>0&&(std::uint64_t(c.nx)+1)*(std::uint64_t(c.ny)+1)==n,"nodal grid count");auto T=values(in,n);auto np=count(in);std::vector<tfv4_patch> patches(np);std::vector<std::array<std::string,4>> strings(np);
  for(size_t i=0;i<np;++i){auto& p=patches[i];for(auto& s:strings[i])s=quoted(in);p.id=strings[i][0].c_str();p.gas_state_id=strings[i][1].c_str();p.parent_id=strings[i][2].c_str();p.scenario_sha256=strings[i][3].c_str();read(in,p.face);read(in,p.begin);read(in,p.end);read(in,p.water_enabled);read(in,p.water_range);read(in,p.initial_coverage);read(in,p.registry_mass_kg);}
  auto electron=quoted(in),parent=quoted(in);complete(in);c.mesh_revision=mesh.c_str();c.material_sha256=material.c_str();c.source_sha256=source.c_str();c.electrode_id=electrode.c_str();c.initial_nodal_temperature_K=T.data();c.initial_node_count=T.size();c.patches=patches.data();c.patch_count=patches.size();c.electron_state_id=electron.c_str();c.electron_parent_id=parent.c_str();call(tfv4_create(&c,&h));
 }else if(cmd=="INFO"){
  int trial;read(in,trial);complete(in);auto i=info(h,trial);out<<' '<<i.revision<<' '<<i.time_s<<' '<<i.node_count<<' '<<i.cell_count<<' '<<i.reference_density_kg_m3<<' '<<i.stored_reference_energy_J<<' '<<i.cumulative_mechanical_exchange_J<<' '<<i.minimum_supported_temperature_K<<' '<<i.maximum_supported_temperature_K<<' '<<i.state_sha256<<' '<<i.patch_count<<' '<<i.layer_count<<' '<<i.gas_packet_count<<' '<<i.path_count;
 }else if(cmd=="TEMP"||cmd=="CELL"){
  int trial;read(in,trial);complete(in);auto i=info(h,trial);std::vector<double> v(cmd=="TEMP"?i.node_count:i.cell_count);if(cmd=="TEMP")call(tfv4_copy_nodal_temperature(h,trial,v.data(),v.size()));else call(tfv4_copy_cell_energy(h,trial,v.data(),v.size()));output(out,v);
 }else if(cmd=="VOL"||cmd=="COORD"){
  complete(in);auto i=info(h);std::vector<double> x(i.node_count),y(i.node_count);if(cmd=="VOL"){call(tfv4_copy_nodal_reference_volume(h,x.data(),x.size()));output(out,x);}else{call(tfv4_copy_node_coordinates(h,x.data(),y.data(),x.size()));for(size_t k=0;k<x.size();++k)out<<' '<<x[k]<<' '<<y[k];}
 }else if(cmd=="PREP"){
  tfv4_epoch e{};read(in,e.base_revision);read(in,e.mechanical_sequence);read(in,e.start_time_s);read(in,e.dt_s);read(in,e.stage_time_s);auto mesh=quoted(in),material=quoted(in),source=quoted(in),em=quoted(in),base=quoted(in),tx=quoted(in);read(in,e.common_state);read(in,e.field_state);read(in,e.source_mesh_revision);read(in,e.source_material_revision);read(in,e.field_time_s);read_array(in,e.E_V_m);read_array(in,e.B_T);read(in,e.material_phase_admitted);read(in,e.planar_barrier_admitted);read(in,e.electron_supply_admitted);
  auto n=count(in);require(n==info(h).cell_count,"mechanical cell count");auto mech=values(in,n);auto na=count(in);std::vector<tfv4_source_action> actions(na);for(auto& a:actions)action(in,a);auto nd=count(in);std::vector<tfv4_water_draw> draws(nd);for(auto& d:draws){read(in,d.patch_index);read(in,d.normal);read(in,d.radius);read(in,d.angle);read(in,d.along);read(in,d.z);read(in,d.time);}complete(in);e.mesh_revision=mesh.c_str();e.material_sha256=material.c_str();e.source_sha256=source.c_str();e.em_state_sha256=em.c_str();e.base_state_sha256=base.c_str();e.transaction_id=tx.c_str();tfv4_budget b{};call(tfv4_prepare(h,&e,mech.data(),mech.size(),actions.data(),actions.size(),draws.data(),draws.size(),&b));budget(out,b);
 }else if(cmd=="LEDGER"){
  complete(in);auto n=info(h,1).node_count;std::array<std::vector<double>,7> v;for(auto& a:v)a.resize(n);call(tfv4_copy_trial_nodal_ledger(h,v[0].data(),v[1].data(),v[2].data(),v[3].data(),v[4].data(),v[5].data(),v[6].data(),n));for(const auto& a:v)output(out,a);
 }else if(cmd=="MCELL"){
  complete(in);std::vector<double> v(info(h,1).cell_count);call(tfv4_copy_trial_cell_mechanical(h,v.data(),v.size()));output(out,v);
 }else if(cmd=="LAYERS"||cmd=="PATHS"||cmd=="GAS"||cmd=="RESERVOIRS"||cmd=="CUMULATIVE"){
  int trial;read(in,trial);complete(in);auto i=info(h,trial);
  if(cmd=="LAYERS"){std::vector<tfv4_layer_info> v(i.layer_count);call(tfv4_copy_layers(h,trial,v.data(),v.size()));for(const auto& l:v){out<<' '<<l.layer_id<<' '<<l.patch_index<<' '<<l.status<<' '<<l.weight<<' '<<l.time_s;output(out,l.origin_m);output(out,l.position_m);output(out,l.proper_velocity_m_s);out<<' '<<l.current_weighted_K_J<<' '<<l.work_function_eV;output(out,l.outward_normal);}}
  else if(cmd=="PATHS"){std::vector<tfv4_path_info> v(i.path_count);call(tfv4_copy_paths(h,trial,v.data(),v.size()));for(const auto& p:v){out<<' '<<p.layer_id;output(out,p.from_m);output(out,p.to_m);out<<' '<<p.start_time_s<<' '<<p.end_time_s<<' '<<p.charge_C<<' '<<p.field_work_J;}}
  else if(cmd=="GAS"){std::vector<tfv4_gas_info> v(i.gas_packet_count);call(tfv4_copy_gas_packets(h,trial,v.data(),v.size()));for(const auto& p:v){out<<' '<<p.patch_index<<' '<<p.source_revision<<' '<<p.weight<<' '<<p.T_K<<' '<<p.time_s;output(out,p.position_m);output(out,p.proper_velocity_m_s);output(out,p.wall_support_impulse_N_s);for(double d:{p.mass_kg,p.charge_C,p.kinetic_J,p.internal_J,p.reference_J,p.adsorbed_reference_change_J})out<<' '<<d;}}
  else{std::vector<double> v(cmd=="RESERVOIRS"?i.patch_count*5:23);if(cmd=="RESERVOIRS")call(tfv4_copy_reservoirs(h,trial,v.data(),v.size()));else call(tfv4_copy_cumulative(h,trial,v.data(),v.size()));output(out,v);}
 }else if(cmd=="PREFLIGHT"||cmd=="COMMIT"){
  auto tx=quoted(in),base=quoted(in);complete(in);if(cmd=="PREFLIGHT")call(tfv4_preflight(h,tx.c_str(),base.c_str()));else call(tfv4_commit(h,tx.c_str(),base.c_str()));
 }else if(cmd=="DISCARD"){auto tx=quoted(in);complete(in);call(tfv4_discard(h,tx.c_str()));}
 else if(cmd=="SAVE"||cmd=="LOAD"){auto path=quoted(in);complete(in);if(cmd=="SAVE")call(tfv4_save(h,path.c_str()));else call(tfv4_load(h,path.c_str()));}
 else if(cmd=="EXIT"){complete(in);exit=true;}
 else throw std::invalid_argument("unknown native v4 command");
 std::cout<<out.str()<<'\n'<<std::flush;
 }catch(const std::exception& e){std::string s=e.what();for(char& c:s)if(c=='\n'||c=='\r')c=' ';std::cout<<"ERR "<<s<<'\n'<<std::flush;}catch(...){std::cout<<"ERR unknown native v4 failure\n"<<std::flush;}
 if(exit)break;
 }tfv4_destroy(h);return 0;
}
