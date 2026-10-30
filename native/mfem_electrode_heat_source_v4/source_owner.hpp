#pragma once
// Include the sealed Hall axial_surface_path_port.hpp first, as required by v3.
#include "axial_trial_field.hpp"
#include "admission.hpp"
#include HEAT_SOURCE_EXACT_COUNT_HEADER
#include <map>
#include <span>
namespace electrode_heat_v4 {
namespace v3=electrode_heat_v3;namespace eh=electrode_heat_v2;namespace se=surface_emission;
using Count=raikou::mfem_backend::ExactParticleCount;
using V3=se::Vector;using Particle=raikou::mfem_backend::GapParticle;
class Owner;
struct PatchIdentity {std::string physical_surface_id,material_layer_id;double physical_area_m2=0;bool operator==(const PatchIdentity&)const=default;};
namespace detail {
inline std::string mesh_hash(const axial::Mesh&g){Encoding e;e.s("source-v4-actual-axial-mesh");
 for(int x:{g.nx,g.ny,g.nz,g.nz_foil,g.kshort,g.kfeed,g.nn,g.ne,g.nb,g.nex,g.ney,g.nbx,g.nby})e.u(x);
 for(double x:{g.side,g.depth,g.dx,g.dy,g.dz,g.zmin})e.d(x);
 for(const auto*v:{&g.M,&g.H,&g.Y,&g.VE,&g.VB,&g.VN})e.doubles(*v);
 e.u(g.edges.size());for(const auto&x:g.edges){e.u(x.a);e.u(x.b);e.u(x.axis);e.d(x.length);}
 e.u(g.C.size());for(const auto&x:g.C){for(auto i:x.index)e.u(i);e.doubles(x.coefficient);e.u(x.count);}
 e.u(g.cells.size());for(const auto&x:g.cells){for(auto i:x.edges)e.u(i);for(int i:{x.i,x.j,x.k,x.body})e.u(i);e.d(x.volume);}
 e.u(g.node_body.size());for(auto i:g.node_body)e.u(i);
 for(const auto*v:{&g.short_node,&g.fixed_E,&g.fixed_A,&g.outer_node}){e.u(v->size());for(bool x:*v)e.u(x);}
 return e.sha();}
inline void policy(Encoding&e,const se::LayerPolicy&p){for(double x:{p.local_linear_vacuum_extent_m,p.image_force_fraction,p.image_energy_fraction,p.maximum_quantum_magnetic_energy_fraction,p.maximum_gyroangle,p.maximum_relative_normal_step})e.d(x);e.u(p.maximum_substeps);e.u(p.PIC_already_applies_microscopic_single_electron_image_force);e.u(p.own_return_unit_accommodation);}
inline std::string policy_hash(const se::LayerPolicy&p){Encoding e;policy(e,p);return e.sha();}
inline void exchange(Encoding&e,const se::LayerExchange&x){for(double a:{x.physical_electron_number_change,x.particle_charge_change_C,x.electrode_charge_change_C,x.wall_heat_J,x.particle_K_change_J,x.electron_reference_change_J,x.image_energy_change_J,x.field_work_J})e.d(a);}
inline void source_exchange_represented(const se::LayerExchange&unit,const se::LayerExchange&actual,double w){
 const double a[]={unit.physical_electron_number_change,unit.particle_charge_change_C,unit.electrode_charge_change_C,unit.wall_heat_J,unit.particle_K_change_J,unit.electron_reference_change_J,unit.image_energy_change_J,unit.field_work_J};
 const double b[]={actual.physical_electron_number_change,actual.particle_charge_change_C,actual.electrode_charge_change_C,actual.wall_heat_J,actual.particle_K_change_J,actual.electron_reference_change_J,actual.image_energy_change_J,actual.field_work_J};
 const char* names[]={"particle count","particle charge C","electrode charge C","wall heat J","particle kinetic increment J","electron reference J","image energy J","field work J"};
 for(int i=0;i<8;++i)require_product(a[i],w,b[i],names[i]);}
inline void admit_weighted_exchange(const se::LayerExchange&u,double w){se::LayerExchange predicted{u.physical_electron_number_change*w,u.particle_charge_change_C*w,u.electrode_charge_change_C*w,u.wall_heat_J*w,u.particle_K_change_J*w,u.electron_reference_change_J*w,u.image_energy_change_J*w,u.field_work_J*w};source_exchange_represented(u,predicted,w);}
inline void configuration(Encoding&e,const eh::Geometry&g,const eh::Material&m,const std::vector<eh::Patch>&patches,double initial_T,double start,const std::string&source){e.u(g.nx);e.u(g.ny);for(double x:{g.x0,g.y0,g.width,g.thickness,g.height,initial_T,start})e.d(x);e.s(g.revision);e.s(g.electrode);for(const auto*p:{&m.id,&m.phase,&m.data_hash,&m.volume_convention,&m.native_callback_contract})e.s(*p);e.s(source);e.u(m.knots.size());for(const auto&k:m.knots)for(double x:{k.T,k.h,k.k,k.sigma})e.d(x);e.u(patches.size());for(const auto&p:patches){for(const auto*q:{&p.id,&p.state_id,&p.parent_id,&p.source_hash})e.s(*q);e.u(p.desorption_enabled);e.u(int(p.face));for(double x:{p.begin,p.end,p.initial_particles,p.initial_low_particles,p.mass_kg,p.attempt_s,p.activation_eV,p.desorption_eV,p.internal_eV,p.monolayer_m2,p.rate_alpha,p.binding_alpha,p.source_T_min,p.source_T_max,p.prefactor_reference_T})e.d(x);e.u(p.coverage_temkin);e.u(int(p.prefactor));}}
inline void budget(Encoding&e,const eh::Budget&b){for(double x:{b.source_wall_J,b.source_particle_K_change_J,b.source_electron_reference_change_J,b.source_image_energy_change_J,b.source_field_work_J,b.neutral_relativistic_conversion_J,b.work_function_reference_part_J,b.joule_J,b.incoming_K_J,b.outgoing_K_J,b.reference_out_J,b.nottingham_wall_J,b.solid_collision_J,b.authored_heat_J,b.neutral_K_J,b.neutral_internal_J,b.adsorbed_reference_change_J,b.expected_neutral_K_J,b.thermal_delta_J,b.wall_net_J,b.residual_J,b.roundoff_bound_J})e.d(x);e.u(b.newton_iterations);}
inline void count(Encoding&e,const Count&n){e.doubles(n.positive_carriers());}
inline void layer(Encoding&e,const eh::BoundLayerState&s){std::ostringstream f(std::ios::binary);eh::write_bound_layer_state(f,s);e.s(f.str());}
inline void identity(Encoding&e,const PatchIdentity&p){e.s(p.physical_surface_id);e.s(p.material_layer_id);e.d(p.physical_area_m2);}
inline void transfer(Encoding&e,const eh::Transfer&t){for(const auto*s:{&t.producer,&t.state_id,&t.parent_id,&t.channel_id,&t.energy_owner})e.s(*s);e.u(t.sequence);e.u(t.owner_sequence);for(double x:{t.event_time_s,t.x,t.y,t.incoming_K_J,t.outgoing_K_J,t.reference_out_J,t.solid_collision_J,t.authored_heat_J,t.nottingham_wall_J,t.work_function_reference_part_J,t.source_wall_J,t.source_particle_K_change_J,t.source_electron_reference_change_J,t.source_image_energy_change_J,t.source_field_work_J})e.d(x);e.u(t.surface);e.u(t.complete_nottingham);e.u(t.certified_source_exchange);}
}
struct FinalFieldView {const v3::TrialFrame&frame;const axial::Mesh&mesh;const foil::Vec&E;const foil::Vec&B;};
class Field final {
 friend class Owner;
 const std::shared_ptr<const axial::Mesh>mesh_;
 const std::shared_ptr<const v3::AxialTrialField>native_;
 const std::string mesh_sha_,content_sha_;
 static std::string content(const v3::TrialFrame&f,const std::string&m,const foil::Vec&E,const foil::Vec&B){detail::Encoding e;e.s("source-v4-total-midpoint");e.s(m);e.u(f.accepted_revision);e.u(f.trial_stage);e.u(f.mesh_revision);e.u(f.material_revision);for(double x:{f.begin_s,f.end_s,f.dt_s,f.field_stage_s})e.d(x);for(const auto*s:{&f.arena_hash,&f.field_hash,&f.mesh_hash,&f.material_hash,&f.species_hash})e.s(*s);e.doubles(E);e.doubles(B);return e.sha();}
public:
 Field(const axial::Mesh&g,v3::TrialFrame f,foil::Vec E,foil::Vec B):mesh_(std::make_shared<const axial::Mesh>(g)),native_(std::make_shared<const v3::AxialTrialField>(mesh_,f,E,B)),mesh_sha_(detail::mesh_hash(*mesh_)),content_sha_(content(f,mesh_sha_,E,B)){}
 const v3::TrialFrame&frame()const{return native_->frame();}
 const std::string&content_sha256()const{return content_sha_;}const std::string&mesh_sha256()const{return mesh_sha_;}
 bool matches(const FinalFieldView&v)const{return frame()==v.frame&&detail::mesh_hash(v.mesh)==mesh_sha_&&content(v.frame,mesh_sha_,v.E,v.B)==content_sha_;}
 v3::FieldSample surface(V3 x,V3 n)const{return native_->surface(x,n);}
};
struct SourceSettings {se::Config spectrum;se::LayerPolicy layer;std::string sampling_policy_identity;};
struct NeutralLeaf {double physical_weight=0;eh::Draw draw{};}; // explicit authored leaves; never adjusted
struct ReleaseQuote {std::size_t patch_index=0;PatchIdentity identity;Count released;double log_released=-INFINITY;double T_K=0,fraction=0,hazard=0,log_hazard=0,reference_eV=0;};
class PreparedPatch {
 friend class Owner;
 struct Data {v3::PreparedPatch native;std::weak_ptr<const char>life;std::uint64_t generation=0;std::shared_ptr<const Field>field;PatchIdentity identity;se::LayerPolicy policy;std::string sampling_policy,owner_signature,content;};
 std::shared_ptr<const Data>d_;explicit PreparedPatch(Data d):d_(std::make_shared<const Data>(std::move(d))){}
public:
 const std::string&exposure_key()const{return d_->native.exposure_key;}double expected_number()const{return d_->native.native.expected_physical_electrons;}double log_expected_number()const{return d_->native.log_expected_number;}const PatchIdentity&identity()const{return d_->identity;}double temperature_K()const{return d_->native.native.stage.T_K;}
};
class Candidate {
 friend class Owner;
 struct Data {v3::Candidate native;std::vector<eh::Transfer>heat_events;std::weak_ptr<const char>life;std::uint64_t generation=0;std::shared_ptr<const Field>field;PatchIdentity identity;se::LayerPolicy policy;std::size_t patch_index=0;double slice_center_z_m=0,extrusion_height_m=0;bool birth=false;std::string owner_signature,prepared_content,content;};
 std::shared_ptr<const Data>d_;explicit Candidate(Data d):d_(std::make_shared<const Data>(std::move(d))){}
public:
 Candidate()=default;
 const eh::BoundLayerState&state()const{return d_->native.state;}const se::LayerExchange&exchange()const{return d_->native.exchange;}const eh::Transfer&aggregate_exchange_receipt()const{return d_->native.heat;}const std::vector<eh::Transfer>&heat_events()const{return d_->heat_events;}const std::vector<v3::PrecisionPath>&paths()const{return d_->native.paths;}std::uint64_t canonical_id()const{return d_->native.canonical_particle_id;}const std::string&content_sha256()const{return d_->content;}const std::optional<v3::BirthSample>&birth_sample()const{return d_->native.birth_sample;}
};
class Handoff {
 friend class Owner;
 struct Data {v3::HandoffTrial native;std::string source_content,field_content,content;};std::shared_ptr<const Data>d_;explicit Handoff(Data d):d_(std::make_shared<const Data>(std::move(d))){}
public:
 Handoff()=default;const Particle&particle()const{return d_->native.particle;}const foil::SurfacePathReceipt&current()const{return d_->native.current;}double handoff_time_s()const{return d_->native.handoff_time_s;}double remaining_duration_s()const{return d_->native.remaining_duration_s;}int boundary_surface()const{return d_->native.boundary.surface;}const std::string&content_sha256()const{return d_->content;}
};
struct LayerRecord {eh::BoundLayerState state;PatchIdentity identity;se::LayerPolicy policy;double slice_center_z_m=0,extrusion_height_m=0;bool qualification=false;std::string energy_owner;};
class Trial {
 friend class Owner;
 struct Data {eh::Trial heat;std::shared_ptr<const Field>field;std::vector<Candidate>sources;std::vector<Handoff>handoffs;bool qualification=false;std::vector<Count>inventory;std::map<std::uint64_t,LayerRecord>layers;std::vector<eh::Packet>packets;std::weak_ptr<const char>life;std::uint64_t generation=0;std::string source_base,content;eh::Budget reported_budget;std::array<double,3>neutral_cumulative{};double neutral_K_J=0,internal_J=0,adsorbed_reference_J=0,expected_neutral_K_J=0,neutral_relativistic_conversion_J=0;};
 std::unique_ptr<Data>d_;explicit Trial(std::unique_ptr<Data>d):d_(std::move(d)){}
public:
 Trial(Trial&&)=default;Trial&operator=(Trial&&)=default;Trial(const Trial&)=delete;Trial&operator=(const Trial&)=delete;
 const std::vector<eh::Packet>&packets()const{return d_->packets;}const Count&inventory(std::size_t j)const{return d_->inventory.at(j);}const eh::Budget&heat_budget()const{return d_->reported_budget;}const std::string&content_sha256()const{return d_->content;}double neutral_kinetic_J()const{return d_->neutral_K_J;}
};
class Owner {
 eh::Geometry geometry_;eh::Material material_;std::vector<eh::Patch>patches_;std::vector<PatchIdentity>identities_;std::vector<SourceSettings>settings_;double initial_T_,start_;std::string source_hash_,signature_;
 std::unique_ptr<eh::Coupler>heat_;std::vector<Count>inventory_;std::map<std::uint64_t,LayerRecord>layers_;
 std::array<double,3>neutral_cumulative_{};
 std::shared_ptr<const char>life_=std::make_shared<const char>(0);std::uint64_t generation_=0;
 static std::vector<eh::Patch>heat_only(std::vector<eh::Patch>p){for(auto&x:p){x.desorption_enabled=false;x.coverage_temkin=false;x.initial_particles=x.initial_low_particles=0;}return p;}
 std::unique_ptr<eh::Coupler>new_heat()const{return std::make_unique<eh::Coupler>(geometry_,material_,heat_only(patches_),initial_T_,source_hash_,start_);}
 bool live(const std::weak_ptr<const char>&w,std::uint64_t g)const{return w.lock()==life_&&g==generation_;}
 void epoch(const eh::Epoch&e)const{eh::need(e.t==heat_->state().time&&e.revision==heat_->state().revision&&e.mesh_revision==geometry_.revision&&e.material_hash==material_.data_hash&&e.source_hash==source_hash_&&std::isfinite(e.dt)&&e.dt>0&&e.t+e.dt>e.t,"v4 accepted epoch mismatch");}
 static std::string candidate_content(const Candidate::Data&d){detail::Encoding e;e.s("source-v4-candidate");e.s(d.owner_signature);e.s(d.prepared_content);e.s(d.field->content_sha256());detail::identity(e,d.identity);detail::policy(e,d.policy);e.u(d.patch_index);e.d(d.slice_center_z_m);e.d(d.extrusion_height_m);e.u(d.birth);e.u(d.native.canonical_particle_id);e.s(d.native.exposure_key);e.s(d.native.count_owner_token);e.u(d.native.qualification_only);detail::layer(e,d.native.state);detail::exchange(e,d.native.exchange);detail::transfer(e,d.native.heat);e.u(d.heat_events.size());for(const auto&t:d.heat_events)detail::transfer(e,t);e.u(d.native.paths.size());for(const auto&p:d.native.paths){e.doubles(p.origin_m);e.doubles(p.from_relative_m);e.doubles(p.to_relative_m);e.doubles(p.displacement_m);for(double x:{p.start_s,p.end_s,p.charge_C,p.field_work_J})e.d(x);}if(d.native.birth_sample){e.u(1);const auto&s=*d.native.birth_sample;for(double x:{s.physical_weight,s.normal_uniform,s.transverse_uniform,s.azimuth_uniform,s.time_uniform,s.specified_birth_s})e.d(x);e.u(s.event_id);e.u(s.canonical_particle_id);e.u(s.count_owner_revision);e.s(s.count_owner_base_hash);e.s(s.count_owner_token);e.s(s.exposure_key);e.u(int(s.time_policy));e.u(int(s.weight_provenance));}else e.u(0);return e.sha();}
 void prepared(const PreparedPatch&p)const{eh::need(p.d_&&live(p.d_->life,p.d_->generation)&&p.d_->owner_signature==signature_,"v4 stale/foreign prepared source");}
 void candidate(const Candidate&c,bool qualification=false)const{eh::need(c.d_&&live(c.d_->life,c.d_->generation)&&c.d_->owner_signature==signature_&&c.d_->content==candidate_content(*c.d_)&&v3::preflight(c.d_->native,*heat_,c.d_->field->native_,qualification),"v4 stale/foreign/changed source candidate");}
 Candidate seal(v3::Candidate c,const PreparedPatch&p,bool birth)const{double w=c.state.native.particle.weight;detail::require_product(-se::qe,w,-se::qe*w,"source charge C");double K=se::kinetic_energy_J(c.state.native.particle.proper_velocity_m_s);detail::require_product(K,w,K*w,"source endpoint kinetic J");for(const auto&path:c.paths){auto E=p.d_->field->native_->sample(path).E;double unitW=-se::qe*v3::dot(E,path.displacement_m);detail::require_product(unitW,w,path.field_work_J,"source path work J");detail::require_product(-se::qe,w,path.charge_C,"source path charge C");}Candidate::Data d;d.heat_events.push_back(c.heat);d.native=std::move(c);d.life=life_;d.generation=generation_;d.field=p.d_->field;d.identity=p.d_->identity;d.policy=p.d_->policy;d.patch_index=p.d_->native.patch_index;d.slice_center_z_m=p.d_->native.binding.extrusion_center_z_m;d.extrusion_height_m=p.d_->native.binding.extrusion_height_m;d.birth=birth;d.owner_signature=signature_;d.prepared_content=p.d_->content;d.content=candidate_content(d);return Candidate(std::move(d));}
 static std::string handoff_content(const Handoff::Data&d){detail::Encoding e;e.s(d.source_content);e.s(d.field_content);const auto&p=d.native.particle;e.u(p.id);e.doubles(p.position_m);e.doubles(p.proper_velocity_m_s);for(double x:{p.mass_kg,p.charge_C,p.weight,d.native.handoff_time_s,d.native.remaining_duration_s,d.native.endpoint_time_s,d.native.kinetic_change_J})e.d(x);e.u(d.native.boundary.surface);e.d(d.native.boundary.fraction);e.doubles(d.native.displacement_m);e.doubles(d.native.current.current_load_C_m_s);e.doubles(d.native.current.endpoint_charge_impulse_C);for(double x:{d.native.current.field_work_J,d.native.current.claimed_field_work_J,d.native.current.work_defect_J,d.native.current.continuity_C})e.d(x);e.u(d.native.current.zero_duration_bridges);return e.sha();}
 std::vector<ReleaseQuote>quotes(const eh::Epoch&e)const;
 static std::string trial_content(const Trial::Data&d){detail::Encoding e;e.s(d.source_base);e.s(d.field->content_sha256());e.u(d.qualification);for(const auto&c:d.sources)e.s(c.content_sha256());for(const auto&h:d.handoffs)e.s(h.content_sha256());for(const auto&n:d.inventory)detail::count(e,n);for(const auto&[id,r]:d.layers){e.u(id);detail::layer(e,r.state);detail::identity(e,r.identity);detail::policy(e,r.policy);e.d(r.slice_center_z_m);e.d(r.extrusion_height_m);e.u(r.qualification);e.s(r.energy_owner);}for(const auto&p:d.packets){e.s(p.event_id);e.s(p.state_id);e.s(p.parent_id);e.s(p.channel_id);e.s(p.patch);e.u(p.source_revision);for(double x:{p.weight,p.T,p.time,p.mass_kg,p.charge_C,p.kinetic_J,p.internal_J,p.reference_J,p.adsorbed_reference_change_J,p.classical_sample_K_J,p.relativistic_conversion_J})e.d(x);e.doubles(p.x);e.doubles(p.v);e.doubles(p.proper_velocity_m_s);e.doubles(p.wall_support_impulse_N_s);}detail::budget(e,d.heat.budget());detail::budget(e,d.reported_budget);const auto&s=d.heat.state();e.doubles(s.T);e.d(s.time);e.u(s.revision);e.doubles(s.cumulative);for(const auto&[k,v]:s.watermarks){e.s(k);e.u(v);}e.d(d.neutral_K_J);e.d(d.internal_J);e.d(d.adsorbed_reference_J);e.d(d.expected_neutral_K_J);e.d(d.neutral_relativistic_conversion_J);e.doubles(d.neutral_cumulative);return e.sha();}
public:
 static constexpr std::size_t maximum_leaves_per_patch=4096;
 Owner(eh::Geometry g,eh::Material m,std::vector<eh::Patch>p,std::vector<PatchIdentity>ids,std::vector<SourceSettings>settings,double T,std::string source,double start=0):geometry_(std::move(g)),material_(std::move(m)),patches_(std::move(p)),identities_(std::move(ids)),settings_(std::move(settings)),initial_T_(T),start_(start),source_hash_(std::move(source)){
  eh::need(geometry_.nx>0&&geometry_.ny>0&&geometry_.nx<=1024&&geometry_.ny<=1024&&patches_.size()<=256,"v4 component resource geometry/patch cap");eh::Coupler validation(geometry_,material_,patches_,T,source_hash_,start);eh::need(identities_.size()==patches_.size()&&settings_.size()==patches_.size(),"v4 physical surface/layer identity required per patch");detail::Encoding enc;enc.s("heat-source-owner-v4");enc.s(validation.config_signature());detail::configuration(enc,geometry_,material_,patches_,T,start,source_hash_);enc.s("46e01fe22e67cc9deb0da9320edef4041f4b3ef433de8ce1eca484463b405501");for(std::size_t j=0;j<patches_.size();++j){const auto&i=identities_[j];eh::need(!i.physical_surface_id.empty()&&!i.material_layer_id.empty()&&std::isfinite(i.physical_area_m2)&&i.physical_area_m2>0&&i.physical_area_m2==validation.area(j),"v4 physical patch area/layer identity mismatch");detail::identity(enc,i);const auto&setting=settings_[j];eh::need(!setting.sampling_policy_identity.empty(),"v4 explicit external sampling policy identity required");v3::admit_policy(setting.layer);enc.d(setting.spectrum.work_function_eV);enc.u(setting.spectrum.energy_panels);enc.u(setting.spectrum.action_nodes);enc.d(setting.spectrum.tail_widths);detail::policy(enc,setting.layer);enc.s(setting.sampling_policy_identity);inventory_.emplace_back(patches_[j].initial_particles,patches_[j].initial_low_particles);}signature_=enc.sha();heat_=new_heat();}
 Owner(const Owner&)=delete;Owner&operator=(const Owner&)=delete;Owner(Owner&&)=delete;Owner&operator=(Owner&&)=delete;
 const eh::State&thermal_state()const{return heat_->state();}const std::array<double,3>&neutral_cumulative_J()const{return neutral_cumulative_;}const Count&inventory(std::size_t j)const{return inventory_.at(j);}const std::map<std::uint64_t,LayerRecord>&active_layers()const{return layers_;}const PatchIdentity&patch_identity(std::size_t j)const{return identities_.at(j);}eh::Epoch make_epoch(double dt,const std::string&field)const{return heat_->make_epoch(dt,field);}const std::string&configuration_sha256()const{return signature_;}
 std::string state_sha256()const{detail::Encoding e;e.s(signature_);e.doubles(neutral_cumulative_);const auto&s=heat_->state();e.doubles(s.T);e.d(s.time);e.u(s.revision);e.doubles(s.cumulative);for(const auto&[k,v]:s.watermarks){e.s(k);e.u(v);}for(const auto&n:inventory_)detail::count(e,n);for(const auto&[id,r]:layers_){e.u(id);detail::layer(e,r.state);detail::identity(e,r.identity);detail::policy(e,r.policy);e.d(r.slice_center_z_m);e.d(r.extrusion_height_m);e.u(r.qualification);e.s(r.energy_owner);}return e.sha();}
 std::vector<ReleaseQuote>release_quotes(const eh::Epoch&e)const{return quotes(e);}
 PreparedPatch prepare_patch(const eh::Epoch&e,std::size_t j,v3::PatchBinding binding,std::shared_ptr<const Field>field,V3 origin,const v3::Exposure&exposure)const;
 Candidate begin_source(const PreparedPatch&p,const v3::BirthSample&sample,const eh::SourceEventIdentity&id)const;
 Candidate advance_active(const PreparedPatch&p,std::uint64_t canonical_id,const eh::SourceEventIdentity&id,const v3::AdvanceControls&controls={})const;
 Candidate advance_birth(const PreparedPatch&p,const Candidate&birth,const eh::SourceEventIdentity&continuation_event,const v3::AdvanceControls&controls={})const;
 bool preflight(const Candidate&c,const std::shared_ptr<const Field>&field,const FinalFieldView&actual,bool qualification=false)const noexcept{try{candidate(c,qualification);return c.d_->field==field&&field&&field->matches(actual);}catch(...){return false;}}
 foil::SurfacePathReceipt import(const std::shared_ptr<const Field>&field,const FinalFieldView&actual,std::span<const Candidate>sources,const std::set<std::uint64_t>&excluded,bool qualification=false)const{eh::need(field&&field->matches(actual),"v4 actual final field/mesh changed");std::vector<v3::Candidate>raw;raw.reserve(sources.size());for(const auto&c:sources){candidate(c,qualification);eh::need(c.d_->field==field,"v4 source final field identity refused");raw.push_back(c.d_->native);}return field->native_->import(raw,excluded,qualification);}
 Handoff handoff(const Candidate&c,const FinalFieldView&actual,const Particle&canonical,bool qualification=false)const{candidate(c,qualification);eh::need(c.d_->field->matches(actual),"v4 actual handoff field/mesh changed");Handoff::Data d;d.native=c.d_->field->native_->handoff_remainder(c.d_->native,canonical,qualification);auto unit_source=c.d_->native;unit_source.state.native.particle.weight=1;auto unit_particle=canonical;unit_particle.weight=1;auto unit=c.d_->field->native_->handoff_remainder(unit_source,unit_particle,qualification);detail::require_product(unit.kinetic_change_J,canonical.weight,d.native.kinetic_change_J,"ordinary handoff kinetic increment J");detail::require_product(unit.current.field_work_J,canonical.weight,d.native.current.field_work_J,"ordinary handoff work J");double K=raikou::mfem_backend::pic_kinetic_energy_J(d.native.particle);detail::require_product(K,canonical.weight,K*canonical.weight,"ordinary handoff endpoint kinetic J");d.source_content=c.d_->content;d.field_content=c.d_->field->content_sha256();d.content=handoff_content(d);return Handoff(std::move(d));}
 Trial prepare(const eh::Epoch&e,mfem::Coefficient&joule,const std::shared_ptr<const Field>&field,const FinalFieldView&actual,std::span<const Candidate>sources,std::span<const Handoff>handoffs,const std::vector<std::vector<NeutralLeaf>>&leaves,const std::vector<eh::Transfer>&other={},bool qualification=false);
 bool preflight(const Trial&t,const FinalFieldView&actual)const noexcept;
 bool commit(Trial&&t,const FinalFieldView&actual)noexcept;
 void publish_prevalidated(Trial&&t)noexcept; // sole coordinator, after all-owner preflight in one critical section
 void save(const std::string&path)const;void load(const std::string&path);
};

inline std::vector<ReleaseQuote>Owner::quotes(const eh::Epoch&e)const {
 epoch(e);std::vector<ReleaseQuote>out;
 for(std::size_t j=0;j<patches_.size();++j){const auto&p=patches_[j];ReleaseQuote q;q.patch_index=j;q.identity=identities_[j];
  if(!p.desorption_enabled||inventory_[j].zero()){out.push_back(q);continue;}
  double N=inventory_[j].rounded_for_rate(),T=heat_->surface_temperatures().at(j);q.T_K=T;double theta=0;
  if(p.coverage_temkin){eh::need(T>=p.source_T_min&&T<=p.source_T_max,"remaining water outside explicitly selected source-law domain");theta=N/(q.identity.physical_area_m2*p.monolayer_m2);eh::need(theta>=0&&theta<=1,"coverage outside submonolayer model");}
  double prefactor=p.attempt_s*(p.prefactor==eh::Prefactor::EyringScaled?T/p.prefactor_reference_T:1.);double barrier=p.activation_eV*(1-p.rate_alpha*theta);
  q.log_hazard=std::log(prefactor)+std::log(e.dt)-barrier*eh::qe/(eh::kb*T);double released;
  if(q.log_hazard<-36){q.log_released=std::log(N)+q.log_hazard;released=std::exp(q.log_released);q.fraction=std::exp(q.log_hazard);}
  else{q.hazard=q.log_hazard>700?INFINITY:std::exp(q.log_hazard);q.fraction=-std::expm1(-q.hazard);released=N*q.fraction;q.log_released=std::log(N)+std::log(q.fraction);}
  if(!(std::isfinite(released)&&released>0))throw Refusal("neutral release count",q.log_released);
  q.released=Count(released);auto remaining=inventory_[j].debit(released); // exact, no nextafter stock/weight adjustment
  double theta_after=p.coverage_temkin?remaining.rounded_for_rate()/(q.identity.physical_area_m2*p.monolayer_m2):0;
  q.reference_eV=p.desorption_eV*(1-p.binding_alpha*.5*(theta+theta_after));out.push_back(std::move(q));
 }
 return out;
}
inline PreparedPatch Owner::prepare_patch(const eh::Epoch&e,std::size_t j,v3::PatchBinding binding,std::shared_ptr<const Field>field,V3 origin,const v3::Exposure&exposure)const {
 const auto&setting=settings_.at(j);const auto&config=setting.spectrum;const auto&policy=setting.layer;std::string sampling_policy=setting.sampling_policy_identity;epoch(e);eh::need(field&&!sampling_policy.empty(),"v4 immutable field and declared sampling-policy identity required");v3::admit_policy(policy);eh::need(field->frame().accepted_revision==e.revision&&field->frame().begin_s==e.t&&field->frame().dt_s==e.dt,"v4 common source/field clock mismatch");
 PreparedPatch::Data d;d.native=v3::prepare_patch(*heat_,e,j,std::move(binding),field->native_,origin,exposure,config);d.life=life_;d.generation=generation_;d.field=std::move(field);d.identity=identities_.at(j);d.policy=policy;d.sampling_policy=std::move(sampling_policy);d.owner_signature=signature_;
 if(d.native.has_log_expected_number&&(!(d.native.native.expected_physical_electrons>0)||!std::isfinite(d.native.log_expected_number)))throw Refusal("electron expected exposure count",d.native.log_expected_number);
 detail::Encoding enc;enc.s(d.native.exposure_key);enc.s(signature_);enc.s(d.field->content_sha256());detail::identity(enc,d.identity);detail::policy(enc,policy);enc.s(d.sampling_policy);enc.d(d.native.log_expected_number);enc.d(d.native.native.expected_physical_electrons);enc.d(d.native.native.mean_wall_heat_diagnostic_J);d.content=enc.sha();d.native.exposure_key+="/v4-immutable/"+d.content;return PreparedPatch(std::move(d));
}
inline Candidate Owner::begin_source(const PreparedPatch&p,const v3::BirthSample&sample,const eh::SourceEventIdentity&id)const {
 prepared(p);eh::need(std::isfinite(sample.physical_weight)&&sample.physical_weight>0,"v4 explicit positive representable electron weight required");eh::need(!layers_.contains(sample.canonical_particle_id),"v4 birth ID already owns accepted source layer");
 // Identical marks, unit weight only for numerical admission; never replacement sampling.
 auto unit_sample=sample;unit_sample.physical_weight=1;auto unit=v3::begin_source(p.d_->native,p.d_->policy,unit_sample,id);
 detail::source_exchange_represented(unit.exchange,unit.exchange,1);
 detail::admit_weighted_exchange(unit.exchange,sample.physical_weight);
 auto raw=v3::begin_source(p.d_->native,p.d_->policy,sample,id);detail::source_exchange_represented(unit.exchange,raw.exchange,sample.physical_weight);return seal(std::move(raw),p,true);
}
inline Candidate Owner::advance_active(const PreparedPatch&p,std::uint64_t id,const eh::SourceEventIdentity&event,const v3::AdvanceControls&controls)const {
 prepared(p);const auto&record=layers_.at(id);eh::need(event.energy_owner==record.energy_owner,"v4 accepted canonical energy owner changed");eh::need(record.identity==p.d_->identity&&record.slice_center_z_m==p.d_->native.binding.extrusion_center_z_m&&record.extrusion_height_m==p.d_->native.binding.extrusion_height_m&&detail::policy_hash(record.policy)==detail::policy_hash(p.d_->policy),"v4 accepted layer area/material/policy changed");
 auto unit=record.state;double w=unit.native.particle.weight;unit.native.particle.weight=1;
 auto check=v3::advance_source(p.d_->native,unit,id,p.d_->policy,p.d_->field->frame().dt_s,event,controls);
 detail::admit_weighted_exchange(check.exchange,w);
 auto raw=v3::advance_source(p.d_->native,record.state,id,p.d_->policy,p.d_->field->frame().dt_s,event,controls);detail::source_exchange_represented(check.exchange,raw.exchange,w);raw.qualification_only=record.qualification;return seal(std::move(raw),p,false);
}
inline Candidate Owner::advance_birth(const PreparedPatch&p,const Candidate&birth,const eh::SourceEventIdentity&event,const v3::AdvanceControls&controls)const {
 prepared(p);candidate(birth,true);eh::need(event.energy_owner==birth.aggregate_exchange_receipt().energy_owner,"v4 newborn canonical energy owner changed");eh::need(birth.d_->birth&&birth.d_->field==p.d_->field&&birth.d_->prepared_content==p.d_->content&&birth.state().native.status==se::LayerStatus::Active,"v4 same-trial Active birth required");
 double duration=p.d_->field->frame().end_s-birth.state().native.current_time_s;eh::need(duration>0,"v4 birth has no active remaining interval");
 auto unit=birth.state();double w=unit.native.particle.weight;unit.native.particle.weight=1;auto check=v3::advance_source(p.d_->native,unit,birth.canonical_id(),p.d_->policy,duration,event,controls);
 detail::admit_weighted_exchange(check.exchange,w);
 auto move=v3::advance_source(p.d_->native,birth.state(),birth.canonical_id(),p.d_->policy,duration,event,controls);detail::source_exchange_represented(check.exchange,move.exchange,w);
 auto raw=birth.d_->native;raw.state=move.state;v3::add(raw.exchange,move.exchange);raw.paths.insert(raw.paths.end(),move.paths.begin(),move.paths.end());v3::finalize(raw,p.d_->native,event,raw.state.native.current_time_s);auto combined=seal(std::move(raw),p,true);auto data=*combined.d_;data.heat_events=birth.d_->heat_events;data.heat_events.push_back(move.heat);data.content=candidate_content(data);return Candidate(std::move(data));
}
inline Trial Owner::prepare(const eh::Epoch&e,mfem::Coefficient&joule,const std::shared_ptr<const Field>&field,const FinalFieldView&actual,std::span<const Candidate>sources,std::span<const Handoff>handoffs,const std::vector<std::vector<NeutralLeaf>>&leaves,const std::vector<eh::Transfer>&other,bool qualification){
 epoch(e);eh::need(field&&field->frame().accepted_revision==e.revision&&field->frame().begin_s==e.t&&field->frame().dt_s==e.dt,"v4 common source/field clock mismatch");eh::need(sources.size()<=4096&&handoffs.size()<=4096,"v4 source proposal cap; whole proposal refused");eh::need(generation_<UINT64_MAX&&field&&field->matches(actual)&&leaves.size()==patches_.size(),"v4 common field/mesh or leaf dimensions invalid");
 auto d=std::make_unique<Trial::Data>();d->life=life_;d->generation=generation_;d->source_base=state_sha256();d->field=field;d->sources.assign(sources.begin(),sources.end());d->handoffs.assign(handoffs.begin(),handoffs.end());d->qualification=qualification;d->inventory=inventory_;d->layers=layers_;auto transfers=other;
 std::set<std::uint64_t>seen;std::map<std::uint64_t,const Handoff*>hand;for(const auto&h:handoffs){eh::need(h.d_&&h.d_->content==handoff_content(*h.d_)&&hand.emplace(h.particle().id,&h).second,"v4 changed/duplicate handoff");}
 for(const auto&c:sources){candidate(c,qualification);eh::need(c.d_->field==field&&seen.insert(c.canonical_id()).second,"v4 source final identity/duplicate");auto id=c.canonical_id();eh::need(c.d_->birth?!layers_.contains(id):layers_.contains(id),"v4 source birth/accepted-layer ownership mismatch");
  if(c.state().native.status==se::LayerStatus::Active){eh::need(c.state().native.current_time_s==field->frame().end_s&&!hand.contains(id),"v4 Active layer not at common endpoint");d->layers[id]={c.state(),c.d_->identity,c.d_->policy,c.d_->slice_center_z_m,c.d_->extrusion_height_m,c.d_->native.qualification_only,c.aggregate_exchange_receipt().energy_owner};}
  else{d->layers.erase(id);if(c.state().native.status==se::LayerStatus::PICHandedOff){auto it=hand.find(id);eh::need(it!=hand.end()&&it->second->d_->source_content==c.content_sha256()&&it->second->d_->field_content==field->content_sha256(),"v4 missing or foreign ordinary handoff remainder");}else eh::need(!hand.contains(id),"v4 recollected source also handed off");}
  transfers.insert(transfers.end(),c.d_->heat_events.begin(),c.d_->heat_events.end());
 }
 std::set<std::string>energy_owners;for(const auto&[id,r]:d->layers)eh::need(!r.energy_owner.empty()&&energy_owners.insert(r.energy_owner).second,"v4 multiple Active IDs share an energy owner");for(const auto&[id,_]:layers_)eh::need(seen.contains(id),"v4 accepted Active layer omitted from interval");for(const auto&[id,_]:hand)eh::need(seen.contains(id),"v4 orphan handoff");
 auto q=quotes(e);for(std::size_t j=0;j<patches_.size();++j){const auto&p=patches_[j];const auto&quote=q[j];eh::need(leaves[j].size()<=maximum_leaves_per_patch,"v4 neutral leaf cap; whole proposal refused");Count total;for(const auto&leaf:leaves[j]){eh::need(std::isfinite(leaf.physical_weight)&&leaf.physical_weight>0,"v4 explicit positive neutral weight required");total=total.credited(leaf.physical_weight);}eh::need(total==quote.released,"v4 exact newborn leaf sum differs from exact source debit");if(total.zero())continue;
  for(double w:total.positive_carriers())d->inventory[j]=d->inventory[j].debit(w);
  for(std::size_t k=0;k<leaves[j].size();++k){const auto&leaf=leaves[j][k];const auto&a=leaf.draw;for(double u:{a.normal,a.radius,a.angle,a.along,a.z,a.time})eh::need(std::isfinite(u)&&u>0&&u<1,"v4 caller-owned open neutral draws required");eh::Packet packet;packet.patch=p.id;packet.state_id=p.state_id;packet.parent_id=p.parent_id;packet.channel_id=p.coverage_temkin?"conditional-oxide-water-Temkin-frozen-coverage":"authored-first-order-desorption";packet.source_revision=e.revision+1;packet.event_id="v4/desorption/"+identities_[j].physical_surface_id+"/"+identities_[j].material_layer_id+"/"+p.id+"/"+std::to_string(packet.source_revision)+"/"+std::to_string(k);packet.T=quote.T_K;packet.weight=leaf.physical_weight;
   double along=p.begin+(p.end-p.begin)*a.along;if(p.face==eh::Face::LowY||p.face==eh::Face::HighY)packet.x={geometry_.x0+along*geometry_.width,p.face==eh::Face::LowY?geometry_.y0:geometry_.y0+geometry_.thickness,(a.z-.5)*geometry_.height};else packet.x={p.face==eh::Face::LowX?geometry_.x0:geometry_.x0+geometry_.width,geometry_.y0+along*geometry_.thickness,(a.z-.5)*geometry_.height};
   double vn=std::sqrt(-2*eh::kb*packet.T/p.mass_kg*std::log(a.normal)),vt=std::sqrt(-2*eh::kb*packet.T/p.mass_kg*std::log(a.radius)),angle=2*3.14159265358979323846*a.angle;if(p.face==eh::Face::LowY||p.face==eh::Face::HighY)packet.v={vt*std::cos(angle),(p.face==eh::Face::LowY?-1:1)*vn,vt*std::sin(angle)};else packet.v={(p.face==eh::Face::LowX?-1:1)*vn,vt*std::cos(angle),vt*std::sin(angle)};
   auto unit=packet;unit.weight=1;eh::assign_native_kinematics(unit,p.mass_kg);eh::assign_native_kinematics(packet,p.mass_kg);detail::require_product(unit.kinetic_J,packet.weight,packet.kinetic_J,"neutral kinetic J");
   packet.internal_J=p.internal_eV*eh::qe*packet.weight;packet.reference_J=0;packet.adsorbed_reference_change_J=quote.reference_eV*eh::qe*packet.weight;if(p.internal_eV>0)detail::require_log_product(std::log(p.internal_eV)+std::log(eh::qe)+std::log(packet.weight),packet.internal_J,"neutral internal J");if(quote.reference_eV>0)detail::require_log_product(std::log(quote.reference_eV)+std::log(eh::qe)+std::log(packet.weight),packet.adsorbed_reference_change_J,"adsorbed reference J");for(int axis=0;axis<3;++axis)packet.wall_support_impulse_N_s[axis]=-packet.weight*packet.mass_kg*packet.proper_velocity_m_s[axis];
   double offset=quote.log_hazard<-36?e.dt*a.time:-std::log1p(-a.time*quote.fraction)/quote.hazard*e.dt;packet.time=e.t+offset;double energy=packet.kinetic_J+packet.internal_J+packet.adsorbed_reference_change_J;eh::need(std::isfinite(energy)&&std::isfinite(packet.time),"v4 neutral energy/time overflow");
   eh::Transfer t;t.producer="v4-neutral/"+identities_[j].physical_surface_id+"/"+identities_[j].material_layer_id+"/"+p.id+"/"+std::to_string(k);t.state_id=p.state_id;t.parent_id=p.parent_id;t.channel_id=packet.channel_id;t.energy_owner=t.producer;t.sequence=t.owner_sequence=e.revision+1;t.event_time_s=packet.time;t.x=packet.x[0];t.y=packet.x[1];t.authored_heat_J=-energy;transfers.push_back(t);d->neutral_K_J+=packet.kinetic_J;d->expected_neutral_K_J+=2*eh::kb*packet.T*packet.weight;d->neutral_relativistic_conversion_J+=packet.relativistic_conversion_J;d->internal_J+=packet.internal_J;d->adsorbed_reference_J+=packet.adsorbed_reference_change_J;eh::need(d->packets.size()<4096,"v4 total neutral proposal cap; whole proposal refused");d->packets.push_back(std::move(packet));
  }
 }
 eh::need(d->layers.size()<=4096,"v4 Active layer cap; whole proposal refused");d->heat=heat_->prepare(e,joule,transfers,std::vector<std::vector<eh::Draw>>(patches_.size()));d->reported_budget=d->heat.budget();d->reported_budget.authored_heat_J=0;for(const auto&t:other)d->reported_budget.authored_heat_J+=t.authored_heat_J;d->reported_budget.neutral_K_J=d->neutral_K_J;d->reported_budget.expected_neutral_K_J=d->expected_neutral_K_J;d->reported_budget.neutral_relativistic_conversion_J=d->neutral_relativistic_conversion_J;d->reported_budget.neutral_internal_J=d->internal_J;d->reported_budget.adsorbed_reference_change_J=d->adsorbed_reference_J;d->neutral_cumulative=neutral_cumulative_;for(int i=0;i<3;++i){double amount=i==0?d->neutral_K_J:i==1?d->internal_J:d->adsorbed_reference_J;d->neutral_cumulative[i]+=amount;eh::need(std::isfinite(d->neutral_cumulative[i]),"v4 neutral cumulative energy overflow");}d->content=trial_content(*d);return Trial(std::move(d));
}
inline bool Owner::preflight(const Trial&t,const FinalFieldView&actual)const noexcept{try{if(!t.d_||!live(t.d_->life,t.d_->generation)||t.d_->source_base!=state_sha256()||!heat_->preflight(t.d_->heat)||!t.d_->field->matches(actual)||t.d_->content!=trial_content(*t.d_))return false;for(const auto&c:t.d_->sources){candidate(c,t.d_->qualification);if(c.d_->field!=t.d_->field)return false;}for(const auto&h:t.d_->handoffs)if(!h.d_||h.d_->content!=handoff_content(*h.d_))return false;return true;}catch(...){return false;}}
inline bool Owner::commit(Trial&&t,const FinalFieldView&actual)noexcept{if(!preflight(t,actual))return false;publish_prevalidated(std::move(t));return true;}
inline void Owner::publish_prevalidated(Trial&&t)noexcept{if(!t.d_||!live(t.d_->life,t.d_->generation)||!heat_->preflight(t.d_->heat))std::terminate();if(!heat_->commit(std::move(t.d_->heat)))std::terminate();static_assert(std::is_nothrow_swappable_v<decltype(inventory_)>&&std::is_nothrow_swappable_v<decltype(layers_)>);inventory_.swap(t.d_->inventory);layers_.swap(t.d_->layers);neutral_cumulative_.swap(t.d_->neutral_cumulative);++generation_;t.d_.reset();}
}

#include "component_checkpoint.hpp"
