#include "thermoelastic_v3.h"
#include "ss304_heat_material_v2.hpp"
#include <cstring>
#include <numeric>
#include <set>
#include <type_traits>

namespace {
thread_local std::string last_error;
struct Holder {
 std::unique_ptr<electrode_heat::Coupler> owner;
 std::unique_ptr<electrode_heat::Trial> trial;
 std::string transaction,base_digest,mesh,material,source,electrode;
 double rho=0,min_T=0,max_T=0;
};
std::set<Holder*> live;
using electrode_heat::need;
template<class F> int guarded(F&& fn){try{fn();last_error.clear();return 0;}catch(const std::exception& e){last_error=e.what();return -1;}catch(...){last_error="unknown native thermoelastic failure";return -1;}}
Holder& get(tfv3_handle handle){auto* h=static_cast<Holder*>(handle);need(h&&live.count(h),"invalid native v3 handle");return *h;}
std::string text(const char* x,const char* label){need(x&&*x,std::string("missing ")+label);return x;}
std::string hash(const char* x,const char* label){auto s=text(x,label);need(s.size()==64&&std::all_of(s.begin(),s.end(),[](char c){return(c>='0'&&c<='9')||(c>='a'&&c<='f');}),std::string("invalid SHA256 ")+label);return s;}
void copy(const std::vector<double>& values,double* out,std::size_t count){need(out&&count==values.size(),"output array dimension mismatch");std::copy(values.begin(),values.end(),out);}
const electrode_heat::Trial& trial(Holder& h){need(bool(h.trial),"no prepared native thermal trial");return *h.trial;}
const electrode_heat::State& state(Holder& h,int which){need(which==0||which==1,"state selector must be accepted0 or trial1");return which?trial(h).state():h.owner->state();}
void preflight(Holder& h,const char* id,const char* base){need(h.trial&&text(id,"transaction_id")==h.transaction&&hash(base,"base_state_sha256")==h.base_digest&&h.owner->state_digest()==h.base_digest&&h.owner->preflight(*h.trial),"stale/mismatched prepared mechanical transaction");}
}

extern "C" {
const char* tfv3_last_error(void){return last_error.c_str();}
int tfv3_create(const tfv3_config* c,tfv3_handle* out){if(out)*out=nullptr;return guarded([&]{
 need(c&&out&&c->abi_version==3,"native thermoelastic ABI version must be3");
 need(mfem::GetVersionMajor()==4&&mfem::GetVersionMinor()==10&&mfem::GetVersionPatch()==0,"frozen MFEM4.10.0 required");
 need(c->nx>0&&c->ny>0&&(std::int64_t(c->nx)+1)*(std::int64_t(c->ny)+1)<std::numeric_limits<int>::max(),"native reference grid dimensions");
 const auto count=std::size_t(c->nx+1)*std::size_t(c->ny+1);
 need(c->initial_nodal_temperature_K&&c->initial_node_count==count,"full initial canonical nodal temperatures required");
 auto h=std::make_unique<Holder>();h->mesh=text(c->mesh_revision,"mesh_revision");h->material=hash(c->material_sha256,"material_sha256");h->source=hash(c->source_sha256,"source_sha256");h->electrode=text(c->electrode_id,"electrode_id");h->rho=c->reference_density_kg_m3;
 auto law=electrode_heat::conditional_ss304_v2(ss304v2::Qualification::nominal_1800,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,h->rho,h->material);
 h->min_T=law.knots.front().T;h->max_T=law.knots.back().T;
 electrode_heat::Geometry g;g.nx=c->nx;g.ny=c->ny;g.x0=c->x0_m;g.y0=c->y0_m;g.width=c->width_m;g.thickness=c->thickness_m;g.height=c->height_m;g.revision=h->mesh;g.electrode=h->electrode;
 h->owner=std::make_unique<electrode_heat::Coupler>(g,std::move(law),std::vector<electrode_heat::Patch>{},c->initial_temperature_K,h->source,c->initial_time_s);
 h->owner->initialize_temperature_cartesian(std::vector<double>(c->initial_nodal_temperature_K,c->initial_nodal_temperature_K+count));
 live.insert(h.get());*out=h.release();
 });}
void tfv3_destroy(tfv3_handle handle){auto* h=static_cast<Holder*>(handle);if(h&&live.erase(h))delete h;}
int tfv3_get_state_info(tfv3_handle handle,tfv3_state_info* out){return guarded([&]{auto& h=get(handle);need(out,"missing state-info output");const auto& s=h.owner->state();const auto& g=h.owner->geometry();tfv3_state_info info{};info.revision=s.revision;info.time_s=s.time;info.node_count=s.T.size();info.cell_count=std::size_t(g.nx)*g.ny;info.reference_density_kg_m3=h.rho;info.stored_reference_energy_J=h.owner->stored_energy_J();info.cumulative_mechanical_exchange_J=s.cumulative_mechanical_J;info.minimum_supported_temperature_K=h.min_T;info.maximum_supported_temperature_K=h.max_T;auto digest=h.owner->state_digest();std::memcpy(info.state_sha256,digest.c_str(),65);*out=info;});}
int tfv3_copy_nodal_temperature(tfv3_handle handle,int which,double* out,size_t n){return guarded([&]{auto& h=get(handle);copy(h.owner->canonical_nodal_values(state(h,which).T),out,n);});}
int tfv3_copy_cell_energy(tfv3_handle handle,int which,double* out,size_t n){return guarded([&]{auto& h=get(handle);copy(h.owner->cell_energies(state(h,which)),out,n);});}
int tfv3_copy_nodal_reference_volume(tfv3_handle handle,double* out,size_t n){return guarded([&]{copy(get(handle).owner->nodal_reference_volumes(),out,n);});}
int tfv3_copy_node_coordinates(tfv3_handle handle,double* x,double* y,size_t n){return guarded([&]{auto coords=get(handle).owner->reference_node_coordinates();need(x&&y&&n==coords.size(),"coordinate output dimension mismatch");for(std::size_t i=0;i<n;++i){x[i]=coords[i][0];y[i]=coords[i][1];}});}
int tfv3_prepare(tfv3_handle handle,const tfv3_epoch* e,const double* cell_J,size_t n,tfv3_budget* out){return guarded([&]{
 auto& h=get(handle);need(e&&out&&cell_J&&!h.trial,"prepare requires inputs and no pending trial");
 need(text(e->mesh_revision,"mesh_revision")==h.mesh&&hash(e->material_sha256,"material_sha256")==h.material&&hash(e->source_sha256,"source_sha256")==h.source,"C ABI immutable mesh/material/source mismatch");
 const auto base=hash(e->base_state_sha256,"base_state_sha256");need(base==h.owner->state_digest(),"C ABI stale base state digest");
 const auto& g=h.owner->geometry();need(n==std::size_t(g.nx)*g.ny,"mechanical cell array must cover complete native reference grid");
 electrode_heat::Epoch epoch{e->start_time_s,e->dt_s,e->base_revision,h.mesh,h.material,h.source,hash(e->em_state_sha256,"em_state_sha256")};
 electrode_heat::MechanicalCellExchange mechanical;mechanical.energy_owner="reversible-mechanics/"+h.electrode;mechanical.transaction_id=text(e->transaction_id,"transaction_id");mechanical.base_state_sha256=base;mechanical.sequence=e->mechanical_sequence;mechanical.stage_time_s=e->stage_time_s;mechanical.cell_J.assign(cell_J,cell_J+n);
 mfem::ConstantCoefficient zero_joule(0.);auto candidate=std::make_unique<electrode_heat::Trial>(h.owner->prepare(epoch,zero_joule,{}, {},&mechanical));
 const auto& b=candidate->budget();need(b.joule_J==0&&b.authored_heat_J==0&&b.wall_net_J==0,"internal mechanics cannot masquerade as external/Joule/wall heat");
 tfv3_budget budget{};budget.thermal_delta_J=b.thermal_delta_J;budget.mechanical_exchange_J=b.mechanical_exchange_J;budget.mechanical_absolute_activity_J=b.mechanical_absolute_activity_J;budget.joule_J=b.joule_J;budget.authored_heat_J=b.authored_heat_J;budget.wall_net_J=b.wall_net_J;budget.residual_J=b.residual_J;budget.roundoff_bound_J=b.roundoff_bound_J;budget.accepted_reference_energy_J=h.owner->stored_energy_J();const auto energy=h.owner->cell_energies(candidate->state());budget.trial_reference_energy_J=std::accumulate(energy.begin(),energy.end(),0.);budget.nodal_residual_max_abs_J=b.nodal_residual_max_abs_J;budget.mechanical_load_conservation_residual_J=b.mechanical_load_conservation_residual_J;budget.newton_iterations=b.newton_iterations;
 h.transaction=mechanical.transaction_id;h.base_digest=base;h.trial=std::move(candidate);*out=budget;
 });}
int tfv3_copy_trial_nodal_ledger(tfv3_handle handle,double* delta,double* mechanical,double* conduction,double* residual,size_t n){return guarded([&]{auto& h=get(handle);const auto& t=trial(h);copy(h.owner->canonical_nodal_values(t.node_ledger(0)),delta,n);copy(h.owner->canonical_nodal_values(t.node_ledger(1)),mechanical,n);copy(h.owner->canonical_nodal_values(t.node_ledger(2)),conduction,n);copy(h.owner->canonical_nodal_values(t.node_ledger(4)),residual,n);});}
int tfv3_copy_trial_cell_mechanical(tfv3_handle handle,double* out,size_t n){return guarded([&]{copy(trial(get(handle)).mechanical_cell_J(),out,n);});}
int tfv3_preflight(tfv3_handle handle,const char* id,const char* base){return guarded([&]{preflight(get(handle),id,base);});}
int tfv3_commit(tfv3_handle handle,const char* id,const char* base){return guarded([&]{
 auto& h=get(handle);preflight(h,id,base);
 static_assert(std::is_nothrow_destructible<electrode_heat::Trial>::value,
               "accepted commit cleanup must not throw");
 // Construct a rejection exception only after commit has returned false. Passing
 // commit and an allocating error string to need() would leave argument order
 // unspecified, permitting an allocation failure after the accepted state swap.
 if(!h.owner->commit(std::move(*h.trial)))
  throw std::invalid_argument("native commit rejected stale proposal");
 h.trial.reset();h.transaction.clear();h.base_digest.clear();
});}
int tfv3_discard(tfv3_handle handle,const char* id){return guarded([&]{auto& h=get(handle);need(h.trial&&text(id,"transaction_id")==h.transaction,"discard requires matching pending transaction");h.trial.reset();h.transaction.clear();h.base_digest.clear();});}
int tfv3_save(tfv3_handle handle,const char* path){return guarded([&]{auto& h=get(handle);need(!h.trial,"checkpoint cannot include an unaccepted proposal");h.owner->save(text(path,"checkpoint path"));});}
int tfv3_load(tfv3_handle handle,const char* path){return guarded([&]{auto& h=get(handle);need(!h.trial,"restart cannot replace a prepared proposal");h.owner->load(text(path,"checkpoint path"));});}
}
