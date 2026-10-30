#include "thermoelastic_v4.h"
#include "ss304_heat_material_v2.hpp"
#include "surface_source_binding.hpp"
#include "water_scenario.hpp"
#include <cstring>
#include <numeric>
#include <set>
#include <type_traits>
namespace {
namespace eh=electrode_heat_v2;namespace se=surface_emission;
using eh::need;
thread_local std::string last_error;
struct Path {std::uint64_t layer_id=0;se::CurrentPath native;};
struct Aux {
 std::map<std::uint64_t,eh::BoundLayerState> layers;
 std::vector<Path> paths;std::vector<eh::Packet> gas;std::vector<eh::Count> gas_inventory;
 std::array<double,3> cumulative_charge_number{};std::uint64_t source_sequence=0;
};
struct Pending {eh::Trial thermal;Aux aux;};
struct Holder {
 std::unique_ptr<eh::Coupler> owner;std::unique_ptr<Pending> pending;Aux aux;
 std::string transaction,base_digest,mesh,material,source,electrode,electron_state,electron_parent;
 double rho=0,min_T=0,max_T=0;std::size_t patch_count=0;
};
std::set<Holder*> live;
template<class F>int guarded(F&& fn){try{fn();last_error.clear();return 0;}catch(const std::exception& e){last_error=e.what();return -1;}catch(...){last_error="unknown native v4 failure";return -1;}}
Holder& get(tfv4_handle x){auto* h=static_cast<Holder*>(x);need(h&&live.count(h),"invalid native v4 handle");return *h;}
std::string text(const char* p,const char* label){need(p&&*p,std::string("missing ")+label);return p;}
std::string hash(const char* p,const char* label){auto s=text(p,label);need(s.size()==64&&std::all_of(s.begin(),s.end(),[](char c){return(c>='0'&&c<='9')||(c>='a'&&c<='f');}),std::string("invalid SHA256 ")+label);return s;}
void boolean(int value){need(value==0||value==1,"explicit Boolean must be0 or1");}
const eh::State& state(Holder& h,int which){need(which==0||which==1,"state selector accepted0 or trial1");need(!which||bool(h.pending),"no pending native v4 trial");return which?h.pending->thermal.state():h.owner->state();}
const Aux& auxiliary(Holder& h,int which){state(h,which);return which?h.pending->aux:h.aux;}
void copy(const std::vector<double>& v,double* out,std::size_t n){need(out&&n==v.size(),"output dimension mismatch");std::copy(v.begin(),v.end(),out);}
std::size_t patch_index(const Holder& h,const std::string& id){for(std::size_t i=0;i<h.patch_count;++i)if(h.owner->patch(i).id==id)return i;throw std::invalid_argument("gas packet patch missing");}
#include "joint_state_v4.inc"
std::string digest(Holder& h,int trial){return thermal_v4_detail::sha256(payload(h,trial));}
void preflight(Holder& h,const char* id,const char* base){need(h.pending&&text(id,"transaction_id")==h.transaction&&hash(base,"base_state_sha256")==h.base_digest&&digest(h,0)==h.base_digest&&h.owner->preflight(h.pending->thermal),"stale/mismatched native joint transaction");}
eh::SourceStageBinding binding(const Holder& h,const tfv4_epoch& e,const tfv4_source_action& a){
 eh::SourceStageBinding b;b.common_state=e.common_state;b.thermal_revision=e.base_revision;b.field_state=e.field_state;b.mesh_revision=e.source_mesh_revision;b.material_revision=e.source_material_revision;b.patch_id=a.patch_id;b.field_time=e.field_time_s;b.mesh_identity=h.mesh;b.material_hash=h.material;
 for(int i=0;i<3;++i){b.E[i]=e.E_V_m[i];b.B[i]=e.B_T[i];}
 for(int v:{e.material_phase_admitted,e.planar_barrier_admitted,e.electron_supply_admitted})boolean(v);
 b.material_phase_admitted=e.material_phase_admitted;b.planar_barrier_admitted=e.planar_barrier_admitted;b.electron_supply_admitted=e.electron_supply_admitted;return b;
}
se::LayerPolicy policy(const tfv4_source_action& a){
 boolean(a.PIC_applies_microscopic_image_force);boolean(a.own_return_unit_accommodation);
 se::LayerPolicy p;p.local_linear_vacuum_extent_m=a.local_linear_vacuum_extent_m;p.image_force_fraction=a.image_force_fraction;p.image_energy_fraction=a.image_energy_fraction;p.maximum_quantum_magnetic_energy_fraction=a.maximum_quantum_magnetic_energy_fraction;p.maximum_gyroangle=a.maximum_gyroangle;p.maximum_relative_normal_step=a.maximum_relative_normal_step;p.maximum_substeps=a.maximum_substeps;p.PIC_already_applies_microscopic_single_electron_image_force=a.PIC_applies_microscopic_image_force;p.own_return_unit_accommodation=a.own_return_unit_accommodation;return p;
}
void source_action(Holder& h,Aux& candidate,const eh::Epoch& epoch,const tfv4_epoch& e,const tfv4_source_action& a,std::vector<eh::Transfer>& transfers,std::array<double,3>& charge){
 need(a.layer_id>0&&a.patch_index<h.patch_count&&a.patch_id>0,"source action identity/patch");
 need(candidate.source_sequence<std::numeric_limits<std::uint64_t>::max(),"source sequence overflow");
 const auto seq=++candidate.source_sequence;eh::SourceEventIdentity id{"native-v4-source/"+h.electrode,h.electron_state,h.electron_parent,"sealed-source-layer","electron/"+h.electrode+"/"+std::to_string(a.layer_id),seq,seq};
 auto bind=binding(h,e,a);const auto stage=eh::accepted_source_stage(*h.owner,epoch,a.patch_index,bind);
 auto found=candidate.layers.find(a.layer_id);eh::BoundLayerTrial result;
 if(a.kind==1){
  need(found==candidate.layers.end(),"duplicate electron birth layer ID");se::Config cfg{a.work_function_eV,a.energy_panels,a.action_nodes,a.tail_widths};
  se::Vector origin{a.origin_m[0],a.origin_m[1],a.origin_m[2]};
  result=eh::begin_bound_source(*h.owner,epoch,a.patch_index,bind,cfg,policy(a),a.weight,a.normal_uniform,a.transverse_uniform,a.azimuth_uniform,a.layer_id,origin,id);
 }else if(a.kind==2){
  need(found!=candidate.layers.end()&&found->second.native.status==se::LayerStatus::Active,"advance requires existing active same electron");
  need(found->second.native.particle.patch_id==a.patch_id,"source numeric patch identity changed");
  result=eh::advance_bound_source(*h.owner,epoch,a.patch_index,bind,found->second,policy(a),a.duration_s,id);
 }else if(a.kind==3){
  need(found!=candidate.layers.end()&&found->second.native.status==se::LayerStatus::PICHandedOff,"collection requires existing handed-off electron");
  auto& layer=found->second;need(layer.patch_index==a.patch_index&&layer.native.particle.patch_id==a.patch_id,"collection patch changed");
  need(a.weight==layer.native.particle.weight&&a.work_function_eV==layer.native.work_function_eV,"collection same particle weight/work function mismatch");
  const double K=se::kinetic_energy_J(layer.native.particle.proper_velocity_m_s)*a.weight;
  need(eh::finite(a.incoming_weighted_K_J)&&std::abs(K-a.incoming_weighted_K_J)<=256*eh::eps*(std::abs(K)+std::abs(a.incoming_weighted_K_J)),"collection K differs from same particle; intervening PIC work is not modeled");
  boolean(a.collection_accommodation_admitted);se::Vector wall{a.origin_m[0],a.origin_m[1],a.origin_m[2]};eh::admit_origin(*h.owner,a.patch_index,wall);
  transfers.push_back(eh::collected_source_heat(a.work_function_eV,a.incoming_weighted_K_J,a.weight,a.collection_accommodation_admitted,id,wall,epoch.t+epoch.dt));
  layer.native.status=se::LayerStatus::Recollected;layer.native.current_time_s=epoch.t+epoch.dt;for(int i=0;i<3;++i)layer.native.position_from_surface_m[i]=wall[i]-layer.surface_origin_m[i];
  charge[0]+=a.weight*se::qe;charge[1]-=a.weight*se::qe;charge[2]-=a.weight;return;
 }else throw std::invalid_argument("unsupported source action kind");
 transfers.push_back(result.heat);charge[0]+=result.exchange.particle_charge_change_C;charge[1]+=result.exchange.electrode_charge_change_C;charge[2]+=result.exchange.physical_electron_number_change;
 for(const auto& path:result.global_paths)candidate.paths.push_back({a.layer_id,path});candidate.layers[a.layer_id]=std::move(result.state);
}
}
extern "C" {
const char* tfv4_last_error(void){return last_error.c_str();}
int tfv4_create(const tfv4_config* c,tfv4_handle* out){if(out)*out=nullptr;return guarded([&]{
 need(c&&out&&c->abi_version==4,"native thermoelastic ABI4 required");need(mfem::GetVersionMajor()==4&&mfem::GetVersionMinor()==10&&mfem::GetVersionPatch()==0,"frozen MFEM4.10.0 required");
 need(c->nx>0&&c->ny>0&&(std::int64_t(c->nx)+1)*(std::int64_t(c->ny)+1)<std::numeric_limits<int>::max(),"native grid dimensions");
 const auto n=std::size_t(c->nx+1)*std::size_t(c->ny+1);need(c->initial_nodal_temperature_K&&c->initial_node_count==n,"full canonical initial nodal temperature required");need(c->patch_count<100000&&(!c->patch_count||c->patches),"patch array missing");
 auto h=std::make_unique<Holder>();h->mesh=text(c->mesh_revision,"mesh_revision");h->material=hash(c->material_sha256,"material_sha256");h->source=hash(c->source_sha256,"source_sha256");h->electrode=text(c->electrode_id,"electrode_id");h->electron_state=text(c->electron_state_id,"electron_state_id");h->electron_parent=text(c->electron_parent_id,"electron_parent_id");h->rho=c->reference_density_kg_m3;h->patch_count=c->patch_count;
 need(h->material=="1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d","material owner provenance seal mismatch");
 auto law=eh::conditional_ss304_v2(ss304v2::Qualification::nominal_1800,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,h->rho,h->material);h->min_T=law.knots.front().T;h->max_T=law.knots.back().T;
 eh::Geometry g;g.nx=c->nx;g.ny=c->ny;g.x0=c->x0_m;g.y0=c->y0_m;g.width=c->width_m;g.thickness=c->thickness_m;g.height=c->height_m;g.revision=h->mesh;g.electrode=h->electrode;
 std::vector<eh::Patch> patches;for(std::size_t i=0;i<c->patch_count;++i){const auto& p=c->patches[i];need(p.face>=0&&p.face<=3,"invalid physical face");boolean(p.water_enabled);eh::Patch patch;
  if(p.water_enabled){need(p.water_range==0||p.water_range==1,"explicit supported water range required");patch=eh::oxide_water_patch(g,text(p.id,"patch id"),eh::Face(p.face),p.begin,p.end,p.initial_coverage,p.registry_mass_kg,text(p.gas_state_id,"gas state"),text(p.parent_id,"water parent"),hash(p.scenario_sha256,"water scenario"),eh::WaterRange(p.water_range));}
  else{patch.id=text(p.id,"patch id");patch.face=eh::Face(p.face);patch.begin=p.begin;patch.end=p.end;patch.source_hash=hash(p.scenario_sha256,"patch scenario");}
  patches.push_back(std::move(patch));
 }
 for(size_t i=0;i<patches.size();++i)for(size_t j=0;j<i;++j)if(patches[i].desorption_enabled&&patches[j].desorption_enabled&&patches[i].face==patches[j].face)need(std::max(patches[i].begin,patches[j].begin)>=std::min(patches[i].end,patches[j].end),"overlapping water patches duplicate a physical reservoir");
 h->owner=std::make_unique<eh::Coupler>(g,std::move(law),std::move(patches),c->initial_temperature_K,h->source,c->initial_time_s);h->owner->initialize_temperature_cartesian(std::vector<double>(c->initial_nodal_temperature_K,c->initial_nodal_temperature_K+n));for(size_t i=0;i<h->patch_count;++i){const auto& p=h->owner->patch(i);if(p.desorption_enabled){const double T=h->owner->surface_temperatures()[i];need(T>=p.source_T_min&&T<=p.source_T_max,"initial remaining-water surface outside selected range");}}h->aux.gas_inventory.resize(c->patch_count);live.insert(h.get());*out=h.release();
 });}
void tfv4_destroy(tfv4_handle x){auto* h=static_cast<Holder*>(x);if(h&&live.erase(h))delete h;}
int tfv4_get_state_info(tfv4_handle x,int which,tfv4_state_info* out){return guarded([&]{auto& h=get(x);need(out,"missing info output");const auto& s=state(h,which);const auto& a=auxiliary(h,which);const auto& g=h.owner->geometry();tfv4_state_info i{};i.revision=s.revision;i.time_s=s.time;i.node_count=s.T.size();i.cell_count=std::size_t(g.nx)*g.ny;i.reference_density_kg_m3=h.rho;i.stored_reference_energy_J=0;for(double e:h.owner->cell_energies(s))i.stored_reference_energy_J+=e;i.cumulative_mechanical_exchange_J=s.cumulative_mechanical_J;i.minimum_supported_temperature_K=h.min_T;i.maximum_supported_temperature_K=h.max_T;i.patch_count=h.patch_count;i.layer_count=a.layers.size();i.gas_packet_count=a.gas.size();i.path_count=a.paths.size();auto d=digest(h,which);std::memcpy(i.state_sha256,d.c_str(),65);*out=i;});}
int tfv4_copy_nodal_temperature(tfv4_handle x,int w,double* out,size_t n){return guarded([&]{auto& h=get(x);copy(h.owner->canonical_nodal_values(state(h,w).T),out,n);});}
int tfv4_copy_cell_energy(tfv4_handle x,int w,double* out,size_t n){return guarded([&]{auto& h=get(x);copy(h.owner->cell_energies(state(h,w)),out,n);});}
int tfv4_copy_nodal_reference_volume(tfv4_handle x,double* out,size_t n){return guarded([&]{copy(get(x).owner->nodal_reference_volumes(),out,n);});}
int tfv4_copy_node_coordinates(tfv4_handle x,double* a,double* b,size_t n){return guarded([&]{auto v=get(x).owner->reference_node_coordinates();need(a&&b&&n==v.size(),"coordinate dimensions");for(size_t i=0;i<n;++i){a[i]=v[i][0];b[i]=v[i][1];}});}
int tfv4_prepare(tfv4_handle x,const tfv4_epoch* e,const double* cell_J,size_t n,const tfv4_source_action* actions,size_t na,const tfv4_water_draw* draws,size_t nd,tfv4_budget* out){return guarded([&]{
 auto& h=get(x);need(e&&out&&cell_J&&!h.pending&&(!na||actions)&&(!nd||draws),"prepare inputs/pending state");
 need(text(e->mesh_revision,"mesh")==h.mesh&&hash(e->material_sha256,"material")==h.material&&hash(e->source_sha256,"source")==h.source,"immutable mesh/material/source mismatch");auto base=hash(e->base_state_sha256,"base");need(base==digest(h,0),"joint base state mismatch");
 need(e->mechanical_sequence==e->base_revision+1&&text(e->transaction_id,"transaction")==base+"/mechanical/"+std::to_string(e->mechanical_sequence),"joint transaction identity mismatch");
 const auto& g=h.owner->geometry();need(n==std::size_t(g.nx)*g.ny,"complete mechanical cell grid required");eh::Epoch epoch{e->start_time_s,e->dt_s,e->base_revision,h.mesh,h.material,h.source,hash(e->em_state_sha256,"EM state")};
 for(const auto& [id,l]:h.aux.layers)if(l.native.status==se::LayerStatus::Active)need(l.native.current_time_s==h.owner->state().time,"accepted active layer clock mismatch");
 Aux candidate=h.aux;std::vector<eh::Transfer> transfers;std::array<double,3> charge{};for(size_t i=0;i<na;++i)source_action(h,candidate,epoch,*e,actions[i],transfers,charge);
 for(const auto& [id,layer]:candidate.layers)if(layer.native.status==se::LayerStatus::Active)need(layer.native.current_time_s==epoch.t+epoch.dt,"active layer must advance to common accepted interval end");
 std::vector<std::vector<eh::Draw>> water(h.patch_count);for(size_t i=0;i<nd;++i){const auto& d=draws[i];need(d.patch_index<h.patch_count,"water draw patch index");water[d.patch_index].push_back({d.normal,d.radius,d.angle,d.along,d.z,d.time});}
 eh::MechanicalCellExchange m;m.energy_owner="reversible-mechanics/"+h.electrode;m.base_state_sha256=h.owner->state_digest();m.sequence=e->mechanical_sequence;m.stage_time_s=e->stage_time_s;m.transaction_id=m.base_state_sha256+"/mechanical/"+std::to_string(m.sequence);m.cell_J.assign(cell_J,cell_J+n);
 mfem::ConstantCoefficient zero(0);auto thermal=h.owner->prepare(epoch,zero,transfers,water,&m);const auto& b=thermal.budget();for(size_t i=0;i<h.patch_count;++i){const auto& p=h.owner->patch(i);if(p.desorption_enabled&&thermal.state().inventory[i].rounded_for_rate()>0){const double T=thermal.next_surface_T()[i];need(T>=p.source_T_min&&T<=p.source_T_max,"candidate remaining-water surface outside selected range");}}need(b.joule_J==0&&b.authored_heat_J==0,"source/mechanics must not become Joule/authored heat");
 for(const auto& p:thermal.packets()){need(p.reference_J==0,"water gas reference must remain zero");candidate.gas.push_back(p);auto i=patch_index(h,p.patch);candidate.gas_inventory[i]=eh::add(candidate.gas_inventory[i],p.weight);need(eh::canonical(candidate.gas_inventory[i]),"gas canonical inventory overflow");}
 for(int i=0;i<3;++i){candidate.cumulative_charge_number[i]+=charge[i];need(eh::finite(candidate.cumulative_charge_number[i]),"charge/number ledger overflow");}
 tfv4_budget z{};z.thermal_delta_J=b.thermal_delta_J;z.mechanical_exchange_J=b.mechanical_exchange_J;z.mechanical_absolute_activity_J=b.mechanical_absolute_activity_J;z.joule_J=b.joule_J;z.authored_heat_J=b.authored_heat_J;z.wall_net_J=b.wall_net_J;z.residual_J=b.residual_J;z.roundoff_bound_J=b.roundoff_bound_J;z.accepted_reference_energy_J=h.owner->stored_energy_J();for(double energy:h.owner->cell_energies(thermal.state()))z.trial_reference_energy_J+=energy;z.nodal_residual_max_abs_J=b.nodal_residual_max_abs_J;z.mechanical_load_conservation_residual_J=b.mechanical_load_conservation_residual_J;z.newton_iterations=b.newton_iterations;
 z.source_wall_J=b.source_wall_J;z.source_particle_K_change_J=b.source_particle_K_change_J;z.source_electron_reference_change_J=b.source_electron_reference_change_J;z.source_image_energy_change_J=b.source_image_energy_change_J;z.source_field_work_J=b.source_field_work_J;z.neutral_K_J=b.neutral_K_J;z.neutral_internal_J=b.neutral_internal_J;z.adsorbed_reference_change_J=b.adsorbed_reference_change_J;z.particle_charge_change_C=charge[0];z.electrode_charge_change_C=charge[1];z.physical_electron_number_change=charge[2];z.new_gas_packet_count=thermal.packets().size();z.new_path_count=candidate.paths.size()-h.aux.paths.size();
 auto pending=std::make_unique<Pending>(Pending{std::move(thermal),std::move(candidate)});h.transaction=e->transaction_id;h.base_digest=base;h.pending=std::move(pending);*out=z;
 });}
int tfv4_copy_trial_nodal_ledger(tfv4_handle x,double* d,double* m,double* c,double* r,double* ext,double* electron,double* water,size_t n){return guarded([&]{auto& h=get(x);state(h,1);const auto& t=h.pending->thermal;copy(h.owner->canonical_nodal_values(t.node_ledger(0)),d,n);copy(h.owner->canonical_nodal_values(t.node_ledger(1)),m,n);copy(h.owner->canonical_nodal_values(t.node_ledger(2)),c,n);copy(h.owner->canonical_nodal_values(t.node_ledger(4)),r,n);copy(h.owner->canonical_nodal_values(t.node_ledger(5)),electron,n);copy(h.owner->canonical_nodal_values(t.node_ledger(6)),water,n);need(ext,"external nodal output missing");for(size_t i=0;i<n;++i)ext[i]=electron[i]+water[i];});}
int tfv4_copy_trial_cell_mechanical(tfv4_handle x,double* out,size_t n){return guarded([&]{auto& h=get(x);state(h,1);copy(h.pending->thermal.mechanical_cell_J(),out,n);});}
int tfv4_copy_layers(tfv4_handle x,int w,tfv4_layer_info* out,size_t n){return guarded([&]{auto& h=get(x);const auto& a=auxiliary(h,w);need(n==a.layers.size()&&(!n||out),"layer output dimensions");size_t k=0;for(const auto& [id,l]:a.layers){tfv4_layer_info v{};v.layer_id=id;v.patch_index=l.patch_index;v.status=int(l.native.status);v.weight=l.native.particle.weight;v.time_s=l.native.current_time_s;v.work_function_eV=l.native.work_function_eV;for(int j=0;j<3;++j){v.origin_m[j]=l.surface_origin_m[j];v.position_m[j]=l.surface_origin_m[j]+l.native.position_from_surface_m[j];v.proper_velocity_m_s[j]=l.native.particle.proper_velocity_m_s[j];v.outward_normal[j]=l.native.outward_normal[j];}v.current_weighted_K_J=l.native.status==se::LayerStatus::Recollected?0:se::kinetic_energy_J(l.native.particle.proper_velocity_m_s)*v.weight;out[k++]=v;}});}
int tfv4_copy_paths(tfv4_handle x,int w,tfv4_path_info* out,size_t n){return guarded([&]{auto& h=get(x);const auto& a=auxiliary(h,w);need(n==a.paths.size()&&(!n||out),"path output dimensions");for(size_t i=0;i<n;++i){tfv4_path_info v{};v.layer_id=a.paths[i].layer_id;const auto& p=a.paths[i].native;for(int j=0;j<3;++j){v.from_m[j]=p.from_m[j];v.to_m[j]=p.to_m[j];}v.start_time_s=p.start_time_s;v.end_time_s=p.end_time_s;v.charge_C=p.charge_C;v.field_work_J=p.field_work_J;out[i]=v;}});}
int tfv4_copy_gas_packets(tfv4_handle x,int w,tfv4_gas_info* out,size_t n){return guarded([&]{auto& h=get(x);const auto& a=auxiliary(h,w);need(n==a.gas.size()&&(!n||out),"gas output dimensions");for(size_t i=0;i<n;++i){const auto& p=a.gas[i];tfv4_gas_info v{};v.patch_index=patch_index(h,p.patch);v.source_revision=p.source_revision;v.weight=p.weight;v.T_K=p.T;v.time_s=p.time;for(int j=0;j<3;++j){v.position_m[j]=p.x[j];v.proper_velocity_m_s[j]=p.proper_velocity_m_s[j];v.wall_support_impulse_N_s[j]=p.wall_support_impulse_N_s[j];}v.mass_kg=p.mass_kg;v.charge_C=p.charge_C;v.kinetic_J=p.kinetic_J;v.internal_J=p.internal_J;v.reference_J=p.reference_J;v.adsorbed_reference_change_J=p.adsorbed_reference_change_J;out[i]=v;}});}
int tfv4_copy_reservoirs(tfv4_handle x,int w,double* out,size_t n){return guarded([&]{auto& h=get(x);const auto& s=state(h,w);const auto& a=auxiliary(h,w);need(n==5*h.patch_count&&(!n||out),"reservoir dimensions5 perpatch");for(size_t i=0;i<h.patch_count;++i){const auto& p=h.owner->patch(i);auto N=s.inventory[i];long double total=(long double)N.hi+N.lo;long double correction=p.coverage_temkin?p.binding_alpha*total/(2*(long double)h.owner->area(i)*p.monolayer_m2):0;out[5*i]=N.hi;out[5*i+1]=N.lo;out[5*i+2]=a.gas_inventory[i].hi;out[5*i+3]=a.gas_inventory[i].lo;out[5*i+4]=p.desorption_enabled?-p.desorption_eV*eh::qe*total*(1-correction):0;}});}
int tfv4_copy_cumulative(tfv4_handle x,int w,double* out,size_t n){return guarded([&]{auto& h=get(x);const auto& s=state(h,w);const auto& a=auxiliary(h,w);need(out&&n==23,"cumulative dimensions23");std::copy(s.cumulative.begin(),s.cumulative.end(),out);out[19]=s.cumulative_mechanical_J;for(int i=0;i<3;++i)out[20+i]=a.cumulative_charge_number[i];});}
int tfv4_preflight(tfv4_handle x,const char* id,const char* base){return guarded([&]{preflight(get(x),id,base);});}
int tfv4_commit(tfv4_handle x,const char* id,const char* base){return guarded([&]{auto& h=get(x);preflight(h,id,base);static_assert(std::is_nothrow_swappable<Aux>::value&&std::is_nothrow_destructible<Pending>::value,"native acceptance cleanup must not throw");if(!h.owner->commit(std::move(h.pending->thermal)))throw std::invalid_argument("native thermal commit rejected");using std::swap;swap(h.aux,h.pending->aux);h.pending.reset();h.transaction.clear();h.base_digest.clear();});}
int tfv4_discard(tfv4_handle x,const char* id){return guarded([&]{auto& h=get(x);need(h.pending&&text(id,"transaction")==h.transaction,"discard matching pending transaction required");h.pending.reset();h.transaction.clear();h.base_digest.clear();});}
int tfv4_save(tfv4_handle x,const char* path){return guarded([&]{auto& h=get(x);need(!h.pending,"checkpoint cannot include pending trial");auto name=text(path,"checkpoint path");auto bytes=payload(h,0);auto checksum=thermal_v4_detail::sha256(bytes);std::ofstream f(name+".tmp",std::ios::binary);f.write(bytes.data(),bytes.size());f.write(checksum.data(),checksum.size());f.close();need(bool(f),"joint checkpoint write failed");std::filesystem::rename(name+".tmp",name);});}
int tfv4_load(tfv4_handle x,const char* path){return guarded([&]{auto& h=get(x);need(!h.pending&&h.owner->state().revision==0,"restart requires fresh native owner");std::ifstream f(text(path,"checkpoint path"),std::ios::binary);need(bool(f),"checkpoint open failed");std::string bytes((std::istreambuf_iterator<char>(f)),{});need(bytes.size()>=64,"checkpoint truncated");auto data=bytes.substr(0,bytes.size()-64);need(bytes.substr(bytes.size()-64)==thermal_v4_detail::sha256(data),"joint checkpoint SHA256 mismatch");restore(h,data);});}
}
