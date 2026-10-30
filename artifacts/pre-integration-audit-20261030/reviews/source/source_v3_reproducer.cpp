#include "hall-owner/development/hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "axial_trial_field.hpp"
#include "ss304_heat_material_v2.hpp"
#include <iostream>
#include <sched.h>
#include <sys/resource.h>
using namespace electrode_heat_v3;
namespace {
int checks=0;void check(bool good,const std::string&why){++checks;if(!good)throw std::runtime_error(why);}
template<class F>void rejects(F&&f,const std::string&why){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,why);}
std::string bytes(const eh::BoundLayerState&s){std::ostringstream out(std::ios::binary);eh::write_bound_layer_state(out,s);return out.str();}
eh::Geometry geom(){eh::Geometry g;g.nx=4;g.ny=2;g.revision="thermal-mesh-v3";g.electrode="upper";return g;}
eh::Material material(){return eh::conditional_ss304_v2(ss304v2::Qualification::source_conditioned_liquid_3000,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000,eh::material_headers_sha256);}
eh::Patch patch(){eh::Patch p;p.id="lower-face";p.face=eh::Face::LowY;p.source_hash="v3-source-scenario";return p;}
eh::SourceEventIdentity event(std::uint64_t seq){return {"source-v3","e-","electrode-electron-reservoir","native-layer","canonical-electron-42",seq,seq};}
V3 origin{.0001,.0025,.0002};
TrialFrame frame(const eh::Coupler&h,std::uint64_t stage=1,double dt=1e-14){TrialFrame f;f.role=FieldRole::TotalMidpoint;f.accepted_revision=h.state().revision;f.trial_stage=stage;f.mesh_revision=1;f.material_revision=2;f.begin_s=h.state().time;f.end_s=f.begin_s+dt;f.dt_s=dt;f.field_stage_s=f.begin_s+.5*dt;f.arena_hash="accepted-arena-revision-"+std::to_string(h.state().revision);f.field_hash="immutable-midpoint-trial-"+std::to_string(stage);f.mesh_hash="Hall-axial-sealed-mesh";f.material_hash="global-material-v2-map";f.species_hash="actual-source-v1-electron";return f;}
PatchBinding binding(const eh::Coupler&h){PatchBinding b;b.thermal_revision=h.state().revision;b.source_accepted_state=h.state().revision+1;b.patch_id=7;b.thermal_mesh_identity=h.geometry().revision;b.thermal_material_hash=eh::material_headers_sha256;b.material_phase_admitted=b.planar_barrier_admitted=b.electron_supply_admitted=true;b.temperature=HeldTemperature::AcceptedPatchMean;b.extrusion_center_z_m=0;b.extrusion_height_m=h.geometry().height;return b;}
struct Arrays{foil::Vec E,B;};
Arrays arrays(const axial::Mesh&g,double F,double B,bool affine=false){Arrays a{foil::Vec(g.ne),foil::Vec(g.nb)};std::vector<V3> nodes(g.nn);for(int k=0;k<=g.nz;++k)for(int j=0;j<=g.ny;++j)for(int i=0;i<=g.nx;++i)nodes[g.node(i,j,k)]={i*g.dx-g.side/2,j*g.dy-g.side/2,g.zmin+k*g.dz};for(int e=0;e<g.ne;++e){const auto&edge=g.edges[e];V3 x{};for(int d=0;d<3;++d)x[d]=.5*(nodes[edge.a][d]+nodes[edge.b][d]);if(edge.axis==0)a.E[e]=affine?1e6+1e13*(x[1]-origin[1]):0;if(edge.axis==1)a.E[e]=F+(affine?1e13*(x[0]-origin[0]):0);if(edge.axis==2)a.E[e]=affine?2e6+1e13*(x[1]-origin[1]):0;}
 for(int k=0;k<g.nz;++k)for(int j=0;j<g.ny;++j)for(int i=0;i<=g.nx;++i)a.B[g.bx(i,j,k)]=B+(affine?1e6*(i*g.dx-g.side/2-origin[0]):0);
 for(int k=0;k<g.nz;++k)for(int j=0;j<=g.ny;++j)for(int i=0;i<g.nx;++i)a.B[g.by(i,j,k)]=affine?-1e6*(j*g.dy-g.side/2-origin[1]):0;
 return a;}
std::shared_ptr<const AxialTrialField> field(const std::shared_ptr<const axial::Mesh>&g,const eh::Coupler&h,std::uint64_t stage,double F,double B=0,bool affine=false){auto a=arrays(*g,F,B,affine);return std::make_shared<const AxialTrialField>(g,frame(h,stage),a.E,a.B);}
PreparedPatch prepared(eh::Coupler&h,std::shared_ptr<const TrialField>f){return prepare_patch(h,h.make_epoch(f->frame().dt_s,"trial-midpoint-field"),0,binding(h),f,origin,Exposure{f->frame().begin_s,f->frame().dt_s,1,"native-surface-quadrature"},se::Config{4.3});}
BirthSample sample(const PreparedPatch&p,double weight=1e9,bool qualification=true){BirthSample b;b.physical_weight=weight;b.normal_uniform=.6;b.transverse_uniform=.35;b.azimuth_uniform=.1;b.time_uniform=.23;b.event_id=123;b.canonical_particle_id=42;b.count_owner_revision=p.field->frame().accepted_revision;b.count_owner_base_hash="unchanged-accepted-source-count-owner";b.count_owner_token=qualification?"qualification-weight-not-an-event-rate":"staged-canonical-fractional-count-receipt";b.exposure_key=p.exposure_key;b.time_policy=BirthTime::ExposureBeginLumped;b.weight_provenance=qualification?WeightProvenance::QualificationOnly:WeightProvenance::StagedCanonicalOwner;return b;}
se::LayerPolicy policy(){se::LayerPolicy p;p.local_linear_vacuum_extent_m=25e-6;return p;}
}

int main(){
 try {
  std::cout<<std::setprecision(17)<<std::boolalpha;
  cpu_set_t affinity; CPU_ZERO(&affinity);sched_getaffinity(0,sizeof(affinity),&affinity);
  std::cout<<"cpu7="<<bool(CPU_ISSET(7,&affinity))<<" affinity_count="<<CPU_COUNT(&affinity)<<" nice="<<getpriority(PRIO_PROCESS,0)<<"\n";
  auto g=std::make_shared<axial::Mesh>(45,4,.009,1);
  eh::Coupler h(geom(),material(),{patch()},1500,"v3-source-scenario",5e-9);
  auto f=field(g,h,1,1e8);auto p=prepared(h,f);auto b=begin_source(p,policy(),sample(p,1e9,false),event(1));
  auto mutated=b;mutated.state.native.particle.weight*=2;
  bool import_accepted=true;try{(void)f->import(std::span<const Candidate>(&mutated,1),{42});}catch(...){import_accepted=false;}
  std::cout<<"mutated_weight_preflight="<<preflight(mutated,h,f)<<" import_accepted="<<import_accepted<<" stored_weight="<<mutated.state.native.particle.weight<<" path_charge="<<mutated.paths[0].charge_C<<"\n";
  auto altered=p;altered.native.spectrum=se::evaluate(1600,1e8,se::Config{4.3});
  auto other=begin_source(altered,policy(),sample(p,1e9,false),event(1));
  std::cout<<"mutated_spectrum_old_exposure_key_accepted="<<(other.exposure_key==p.exposure_key)<<" particle_changed="<<(other.state.native.particle.proper_velocity_m_s!=b.state.native.particle.proper_velocity_m_s)<<"\n";
  eh::Coupler cold(geom(),material(),{patch()},300,"v3-source-scenario",5e-9);
  auto fc=field(g,cold,2,0);auto pc=prepared(cold,fc);
  auto tiny=begin_source(pc,policy(),sample(pc,1e-310,false),event(1));
  std::cout<<"positive_weight="<<tiny.exchange.physical_electron_number_change<<" weighted_charge="<<tiny.exchange.particle_charge_change_C<<" weighted_wall_heat="<<tiny.exchange.wall_heat_J<<" preflight="<<preflight(tiny,cold,fc)<<"\n";
  auto count=eh::add(eh::Count{1e15,-.01},-1e-30);
  std::cout<<"source_count_debit_lost="<<(count.hi==1e15&&count.lo==-.01)<<" hi="<<std::hexfloat<<count.hi<<" lo="<<count.lo<<std::defaultfloat<<"\n";
  auto sub=prepare_patch(cold,cold.make_epoch(fc->frame().dt_s,"trial-midpoint-field"),0,binding(cold),fc,origin,Exposure{fc->frame().begin_s,fc->frame().dt_s,1e-250,"tiny-positive-area"},se::Config{4.3});
  std::cout<<"tiny_area_scalar_count="<<sub.native.expected_physical_electrons<<" log_count="<<sub.log_expected_number<<" reconstructed_count="<<std::exp(sub.log_expected_number)<<" current_underflow_flag="<<sub.current_underflow<<"\n";
  GapParticle canonical;canonical.id=42;canonical.mass_kg=se::me;canonical.charge_C=-se::qe;canonical.weight=1;
  auto thermal=begin_source(pc,policy(),sample(pc,1,false),event(1));
  try{auto moved=fc->handoff_remainder(thermal,canonical);std::cout<<"thermal_surface_handoff_succeeds=true boundary="<<moved.boundary.surface<<" duration="<<moved.remaining_duration_s<<" endpoint_time="<<moved.endpoint_time_s<<"\n";}catch(const std::exception&e){std::cout<<"thermal_surface_handoff_succeeds=false reason="<<e.what()<<"\n";}
  alignas(eh::Coupler) std::byte storage[sizeof(eh::Coupler)];
  auto* old=new(storage)eh::Coupler(geom(),material(),{patch()},1500,"v3-source-scenario",5e-9);
  auto fa=field(g,*old,3,1e8);auto pa=prepared(*old,fa);auto ca=begin_source(pa,policy(),sample(pa,1,false),event(1));
  std::destroy_at(old);auto* replacement=new(storage)eh::Coupler(geom(),material(),{patch()},1500,"v3-source-scenario",5e-9);
  std::cout<<"replacement_owner_accepts_stale_candidate="<<preflight(ca,*replacement,fa)<<"\n";std::destroy_at(replacement);

  auto water_patch=[](double high,double low,double rate){auto p=patch();p.desorption_enabled=true;p.state_id="H2O";p.parent_id="H2O";p.initial_particles=high;p.initial_low_particles=low;p.mass_kg=18*1.66053906660e-27;p.attempt_s=rate;p.activation_eV=p.desorption_eV=p.internal_eV=0;return p;};
  mfem::ConstantCoefficient zero(0);eh::Draw draw{.5,.5,.5,.5,.5,.5};
  eh::Coupler reservoir(geom(),material(),{water_patch(1e15,-.01,1e-31)},1500,"v3-source-scenario",5e-9);
  auto release=reservoir.prepare(reservoir.make_epoch(1e-14,"manufactured-source-audit"),zero,{},{{draw}});
  std::cout<<"actual_positive_packet_weight="<<release.packets()[0].weight<<" actual_source_inventory_unchanged="<<(release.state().inventory[0].hi==reservoir.state().inventory[0].hi&&release.state().inventory[0].lo==reservoir.state().inventory[0].lo)<<"\n";
  eh::Coupler split(geom(),material(),{water_patch(1,0,1e300)},1500,"v3-source-scenario",5e-9);
  auto allocated=split.prepare(split.make_epoch(1e-14,"manufactured-source-audit"),zero,{}, {std::vector<eh::Draw>(10,draw)});
  std::cout<<"allocation_source_remaining_hi="<<allocated.state().inventory[0].hi<<" remaining_lo="<<allocated.state().inventory[0].lo<<" packet_weights_hex=";
  for(const auto&p:allocated.packets())std::cout<<std::hexfloat<<p.weight<<' ';std::cout<<std::defaultfloat<<"\n";
  auto gm=std::make_shared<axial::Mesh>(45,4,.009,1);auto fm=field(gm,h,4,1e8,0,true);auto before=fm->surface(origin,h.normal(0));gm->dx*=1.001;
  bool changed_mesh_sample=false;try{changed_mesh_sample=fm->surface(origin,h.normal(0)).E!=before.E;}catch(...){changed_mesh_sample=true;}
  std::cout<<"mutable_mesh_changes_field_gather="<<changed_mesh_sample<<" matches_trial_still_true="<<fm->matches_trial(fm->frame(),fm->total_mid_E(),fm->total_mid_B())<<"\n";
  return 0;
 } catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}
