// Persistent one-electrode native process. The frozen MFEM archive is non-PIC;
// this executable preserves it unchanged and carries the internal C ABI over a
// strict, deterministic line protocol. It performs no scheduling or launch.
#include "thermoelastic_v3.h"
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool value,const std::string& message){if(!value)throw std::invalid_argument(message);}
void call(int code){if(code)throw std::runtime_error(tfv3_last_error());}
template<class T>void read(std::istream& input,T& value){require(bool(input>>value),"malformed numeric request field");}
std::string quoted(std::istream& input){std::string value;require(bool(input>>std::quoted(value)),"malformed quoted request field");return value;}
void complete(std::istream& input){input>>std::ws;require(input.eof(),"unexpected trailing request fields");}
std::vector<double> values(std::istream& input,std::size_t count){require(count<=10000000,"request vector too large");std::vector<double> result(count);for(double& x:result)read(input,x);return result;}
void output_values(std::ostream& out,const std::vector<double>& values){for(double x:values)out<<' '<<x;}
tfv3_state_info info(tfv3_handle handle){tfv3_state_info result{};call(tfv3_get_state_info(handle,&result));return result;}
void output_budget(std::ostream& out,const tfv3_budget& b){
 out<<' '<<b.thermal_delta_J<<' '<<b.mechanical_exchange_J<<' '<<b.mechanical_absolute_activity_J
    <<' '<<b.joule_J<<' '<<b.authored_heat_J<<' '<<b.wall_net_J<<' '<<b.residual_J<<' '<<b.roundoff_bound_J
    <<' '<<b.accepted_reference_energy_J<<' '<<b.trial_reference_energy_J<<' '<<b.nodal_residual_max_abs_J
    <<' '<<b.mechanical_load_conservation_residual_J<<' '<<b.newton_iterations;
}
}

int main(){
 tfv3_handle handle=nullptr;std::string line;
 while(std::getline(std::cin,line)){
  bool exit=false;
  try{
   std::istringstream request(line);std::string command;require(bool(request>>command),"empty request");
   std::ostringstream response;response<<std::setprecision(17)<<"OK";
   if(command=="CREATE"){
    require(!handle,"CREATE requires a fresh server");tfv3_config c{};
    read(request,c.abi_version);read(request,c.nx);read(request,c.ny);read(request,c.x0_m);read(request,c.y0_m);
    read(request,c.width_m);read(request,c.thickness_m);read(request,c.height_m);read(request,c.reference_density_kg_m3);
    read(request,c.initial_temperature_K);read(request,c.initial_time_s);
    auto mesh=quoted(request),material=quoted(request),source=quoted(request),electrode=quoted(request);
    std::size_t n=0;read(request,n);require(c.nx>0&&c.ny>0&&(std::uint64_t(c.nx)+1)*(std::uint64_t(c.ny)+1)==n,"CREATE nodal count/grid mismatch");
    auto nodes=values(request,n);complete(request);c.mesh_revision=mesh.c_str();c.material_sha256=material.c_str();c.source_sha256=source.c_str();c.electrode_id=electrode.c_str();c.initial_nodal_temperature_K=nodes.data();c.initial_node_count=n;call(tfv3_create(&c,&handle));
   }else if(command=="INFO"){
    complete(request);auto i=info(handle);response<<' '<<i.revision<<' '<<i.time_s<<' '<<i.node_count<<' '<<i.cell_count<<' '<<i.reference_density_kg_m3<<' '<<i.stored_reference_energy_J<<' '<<i.cumulative_mechanical_exchange_J<<' '<<i.minimum_supported_temperature_K<<' '<<i.maximum_supported_temperature_K<<' '<<i.state_sha256;
   }else if(command=="TEMP"||command=="CELL"){
    int trial=0;read(request,trial);complete(request);auto i=info(handle);std::vector<double> out(command=="TEMP"?i.node_count:i.cell_count);
    if(command=="TEMP")call(tfv3_copy_nodal_temperature(handle,trial,out.data(),out.size()));else call(tfv3_copy_cell_energy(handle,trial,out.data(),out.size()));output_values(response,out);
   }else if(command=="COORD"){
    complete(request);auto i=info(handle);std::vector<double>x(i.node_count),y(i.node_count);call(tfv3_copy_node_coordinates(handle,x.data(),y.data(),x.size()));for(std::size_t k=0;k<x.size();++k)response<<' '<<x[k]<<' '<<y[k];
   }else if(command=="VOL"){
    complete(request);std::vector<double>out(info(handle).node_count);call(tfv3_copy_nodal_reference_volume(handle,out.data(),out.size()));output_values(response,out);
   }else if(command=="PREP"){
    tfv3_epoch e{};read(request,e.base_revision);read(request,e.mechanical_sequence);read(request,e.start_time_s);read(request,e.dt_s);read(request,e.stage_time_s);
    auto mesh=quoted(request),material=quoted(request),source=quoted(request),em=quoted(request),base=quoted(request),transaction=quoted(request);
    std::size_t n=0;read(request,n);require(n==info(handle).cell_count,"PREP cell count/grid mismatch");auto energy=values(request,n);complete(request);
    e.mesh_revision=mesh.c_str();e.material_sha256=material.c_str();e.source_sha256=source.c_str();e.em_state_sha256=em.c_str();e.base_state_sha256=base.c_str();e.transaction_id=transaction.c_str();tfv3_budget b{};call(tfv3_prepare(handle,&e,energy.data(),energy.size(),&b));output_budget(response,b);
   }else if(command=="LEDGER"){
    complete(request);auto n=info(handle).node_count;std::vector<double>d(n),m(n),q(n),r(n);call(tfv3_copy_trial_nodal_ledger(handle,d.data(),m.data(),q.data(),r.data(),n));for(const auto* v:{&d,&m,&q,&r})output_values(response,*v);
   }else if(command=="MCELL"){
    complete(request);std::vector<double>out(info(handle).cell_count);call(tfv3_copy_trial_cell_mechanical(handle,out.data(),out.size()));output_values(response,out);
   }else if(command=="PREFLIGHT"||command=="COMMIT"){
    auto transaction=quoted(request),base=quoted(request);complete(request);if(command=="PREFLIGHT")call(tfv3_preflight(handle,transaction.c_str(),base.c_str()));else call(tfv3_commit(handle,transaction.c_str(),base.c_str()));
   }else if(command=="DISCARD"){
    auto transaction=quoted(request);complete(request);call(tfv3_discard(handle,transaction.c_str()));
   }else if(command=="SAVE"||command=="LOAD"){
    auto path=quoted(request);complete(request);if(command=="SAVE")call(tfv3_save(handle,path.c_str()));else call(tfv3_load(handle,path.c_str()));
   }else if(command=="EXIT"){
    complete(request);exit=true;
   }else throw std::invalid_argument("unknown native thermoelastic command");
   std::cout<<response.str()<<'\n'<<std::flush;
  }catch(const std::exception& e){std::string error=e.what();for(char& c:error)if(c=='\n'||c=='\r')c=' ';std::cout<<"ERR "<<error<<'\n'<<std::flush;}
  catch(...){std::cout<<"ERR unknown native thermoelastic failure\n"<<std::flush;}
  if(exit)break;
 }
 tfv3_destroy(handle);return 0;
}
