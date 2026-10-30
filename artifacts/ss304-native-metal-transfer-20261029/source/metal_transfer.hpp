#pragma once
#include "dyadic.hpp"
#include "material/ss304_thermal_v2.hpp"
#include "raikou/mfem_pic_arena.h"
#include <openssl/sha.h>
#include <iomanip>
#include <sstream>

namespace metal_transfer {
namespace m=raikou::mfem_backend;
namespace h=electrode_heat;
using J=nlohmann::json;
inline std::string hash(const std::string& s){unsigned char out[SHA256_DIGEST_LENGTH];SHA256(reinterpret_cast<const unsigned char*>(s.data()),s.size(),out);std::ostringstream o;for(auto x:out)o<<std::hex<<std::setfill('0')<<std::setw(2)<<unsigned(x);return o.str();}
inline D count(h::Count n){return D(n.hi)+D(n.lo);}
inline h::Count credit(h::Count n,double w){auto r=h::add(n,w);need(h::canonical(r)&&count(r).compare(count(n)+D(w))==0,"exact credit exceeds Count support");return r;}
struct Config {
 std::string material_seal,initial_prefix;
 h::Geometry geometry;
 double rho_reference=8000,internal_per_atom_J=0;
 double axial_center_m=0;
 // Only a fixed upper face, uniform thermal block and immobile support here.
 bool conditional_congruent_neutral=true,fixed_support=true;
 bool manufactured_internal_state=false;
};
struct Snapshot {
 h::Count solid_atoms;
 D energy;
 std::array<D,3> support_impulse;
 std::set<std::string> events;
 std::set<std::string> consumed_sinks;
 J journal=J::array();
 std::uint64_t revision=0;
};
struct Event {
 std::string id,state_prefix;
 std::uint64_t base_revision=0;
 std::vector<m::MetalSurfacePacket> packets;
 double other_heat_J=0,dt_s=1e-12,physical_stage_s=0;
 std::string heat_owner;
 double field_work_J=0,extra_latent_J=0;
 std::vector<std::string> native_sink_receipts;
};
class Adapter;
class Prepared {
 friend class Adapter;
 Adapter* owner_=nullptr;
 std::uint64_t generation_=0;
 std::shared_ptr<const Snapshot> before_,after_;
 m::PreparedArena gas_;
 std::unique_ptr<h::Coupler> heat_cache_;
 h::Trial heat_trial_;
 J audit_;
 enum class Status {Staged,Published,RolledBack};Status status_=Status::Staged;
 Prepared(m::PreparedArena gas):gas_(std::move(gas)){}
 public:
 Prepared(Prepared&&) noexcept=default;Prepared&operator=(Prepared&&) noexcept=default;
 Prepared(const Prepared&)=delete;Prepared&operator=(const Prepared&)=delete;
 const J& audit()const{return audit_;}
};
class Adapter {
 Config cfg_;
 m::PicArena* arena_;
 ss304v2::Model material_;
 ss304v2::NeutralVapor vapor_;
 std::shared_ptr<const Snapshot> accepted_;
 std::uint64_t generation_=0;
 double atom_mass_,reference_per_atom_;
 std::string signature_;
 double volume()const{return cfg_.geometry.width*cfg_.geometry.thickness*cfg_.geometry.height;}
 D mass(h::Count n)const{return count(n)*D(atom_mass_);}
 double temperature(const Snapshot& s)const{
   need(count(s.solid_atoms).compare(D(0))>0,"empty solid requires geometry handoff");
   const D kg=mass(s.solid_atoms);
   need(s.energy.compare(kg*D(material_.min_u()))>=0&&s.energy.compare(kg*D(material_.max_u()))<=0,"exact energy outside298.2..3000K domain");
   return material_.temperature_from_u(static_cast<double>(s.energy.ld()/kg.ld()));
 }
 h::Material law(double kg)const {
   const auto material=material_;const double density=kg/volume();
   h::Material a;a.id="conditional-congruent-SS304";a.phase="bounded-solid-mushy-liquid";a.data_hash=cfg_.material_seal;a.volume_convention="current-mass-per-fixed-reference-volume";a.native_callback_contract="uniform-metal-block/v1/"+signature_;
   for(double t:{298.2,300.,1500.,1672.,1727.,1800.,3000.}){auto v=material.evaluate(t);a.knots.push_back({t,density*v.h_j_kg,v.k_w_m_k,v.sigma_s_m});}
   a.h_fn=[=](double t){return density*material.evaluate(t).h_j_kg;};
   a.cv_fn=[=](double t){return density*material.evaluate(t).apparent_heat_capacity_j_kg_k;};
   a.k_fn=[=](double t){return material.evaluate(t).k_w_m_k;};a.sigma_fn=[=](double t){return material.evaluate(t).sigma_s_m;};
   a.phase_fn=[=](double t){return material.evaluate(t).phase;};a.liquid_fraction_fn=[=](double t){return material.evaluate(t).liquid_fraction;};
   a.delta_h_fn=[=](double t,double old){return density*(material.evaluate(t).h_j_kg-material.evaluate(old).h_j_kg);};return a;
 }
 J state_json(const Snapshot& s)const{J impulse=J::array();for(auto& x:s.support_impulse)impulse.push_back(x.json());return {{"solid_atoms",{s.solid_atoms.hi,s.solid_atoms.lo}},{"energy",s.energy.json()},{"support_impulse",impulse},{"events",s.events},{"consumed_sinks",s.consumed_sinks},{"journal",s.journal},{"revision",s.revision}};}
 public:
 Adapter(Config cfg,m::PicArena& arena,h::Count initial_atoms,double initial_T)
  :cfg_(std::move(cfg)),arena_(&arena),material_(ss304v2::Qualification::source_conditioned_liquid_3000,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,cfg_.rho_reference),vapor_(material_) {
  need(cfg_.conditional_congruent_neutral&&cfg_.fixed_support&&cfg_.material_seal=="1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d"&&!cfg_.initial_prefix.empty(),"explicit sealed-v2/conditional/fixed-support binding required");
  cfg_.geometry=h::admit_geometry(cfg_.geometry);need(cfg_.geometry.nx==1&&cfg_.geometry.ny==1,"only bounded uniform1x1 thermal block qualified");
  need(std::isfinite(cfg_.axial_center_m)&&std::isfinite(cfg_.axial_center_m-cfg_.geometry.height/2)&&std::isfinite(cfg_.axial_center_m+cfg_.geometry.height/2)&&cfg_.axial_center_m-cfg_.geometry.height/2<cfg_.axial_center_m+cfg_.geometry.height/2,"unresolved axial source frame");
  need(h::canonical(initial_atoms)&&count(initial_atoms).compare(D(0))>0,"initial metal inventory");
  atom_mass_=ss304v2::NeutralVapor::molar_mass/ss304v2::NeutralVapor::avogadro;
  reference_per_atom_=atom_mass_*vapor_.atomic_u0();
  auto s=std::make_shared<Snapshot>();s->solid_atoms=initial_atoms;s->energy=mass(initial_atoms)*D(material_.evaluate(initial_T).u_model_j_kg);accepted_=s;
  J config={{"native_source_sha256",m::PicArena::source_identity()},{"material",cfg_.material_seal},{"initial_prefix",cfg_.initial_prefix},{"rho",cfg_.rho_reference},{"internal_per_atom_J",cfg_.internal_per_atom_J},{"geometry",{cfg_.geometry.nx,cfg_.geometry.ny,cfg_.geometry.x0,cfg_.geometry.y0,cfg_.geometry.width,cfg_.geometry.thickness,cfg_.geometry.height}},{"mesh_revision",cfg_.geometry.revision},{"electrode",cfg_.geometry.electrode},{"initial_count",{initial_atoms.hi,initial_atoms.lo}},{"initial_T",initial_T}};
  need(std::isfinite(cfg_.internal_per_atom_J)&&cfg_.internal_per_atom_J>=0&&(cfg_.internal_per_atom_J==0||cfg_.manufactured_internal_state),"nonzero internal states require explicit manufactured qualification; no physical excitation data supplied");config["manufactured_internal_state"]=cfg_.manufactured_internal_state;config["axial_center_m"]=cfg_.axial_center_m;signature_=hash(config.dump());
 }
 m::ArenaStateSpec registry_state()const{return {"SS304_CONGRUENT_NEUTRAL",0,0,0,atom_mass_,cfg_.internal_per_atom_J,reference_per_atom_};}
 static m::ArenaStateSpec required_registry_state(double internal=0){ss304v2::Model m0(ss304v2::Qualification::source_conditioned_liquid_3000,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000);ss304v2::NeutralVapor v(m0);const double mass0=ss304v2::NeutralVapor::molar_mass/ss304v2::NeutralVapor::avogadro;return {"SS304_CONGRUENT_NEUTRAL",0,0,0,mass0,internal,mass0*v.atomic_u0()};}
 const std::shared_ptr<const Snapshot>& snapshot()const{return accepted_;}
 double temperature()const{return temperature(*accepted_);}
 std::string prefix()const{return hash(signature_+state_json(*accepted_).dump()+arena_->state_prefix());}
 // Instantaneous upper-face proposal from five caller-owned open-unit draws.
 // The heat API performs physical-v to native proper-u/kinetic conversion.
 // No inventory debit, RNG or mean-energy replacement occurs in this helper.
 m::MetalSurfacePacket sample(std::string id,h::Count atoms,const std::array<double,5>& draws,h::Packet* per_atom=nullptr)const {
  const double T=temperature();need(T>=1800,"metal sample below vapor scenario support");
  for(double u:draws)need(std::isfinite(u)&&u>0&&u<1,"open caller-owned metal sampling draws required");
  need(!m::PicArena::positive_carriers(atoms).empty(),"positive sampled count required");
  h::Packet p;p.weight=1;p.T=T;p.x={cfg_.geometry.x0+draws[3]*cfg_.geometry.width,std::nextafter(cfg_.geometry.y0+cfg_.geometry.thickness,INFINITY),cfg_.axial_center_m+(draws[4]-.5)*cfg_.geometry.height};
  const double vn=std::sqrt(-2*h::kb*T/atom_mass_*std::log(draws[0])),vt=std::sqrt(-2*h::kb*T/atom_mass_*std::log(draws[1])),angle=2*std::acos(-1.)*draws[2];
  p.v={vt*std::cos(angle),vn,vt*std::sin(angle)};h::assign_native_kinematics(p,atom_mass_);p.internal_J=cfg_.internal_per_atom_J;p.reference_J=reference_per_atom_;
  m::GapParticle pose;pose.position_m=p.x;pose.proper_velocity_m_s=p.proper_velocity_m_s;pose.mass_kg=atom_mass_;pose.charge_C=0;pose.weight=1;
  if(per_atom)*per_atom=p;
  return {true,std::move(id),atoms,pose};
 }
 Prepared prepare(const Event& e,int fault=0) {need(e.native_sink_receipts.empty(),"native receipts require typed staged transport");return prepare_impl(e,fault,nullptr);}
 Prepared prepare_native_collection(Event e,m::PreparedArena&& transport,int fault=0) {
  need(e.packets.empty()&&e.native_sink_receipts.empty()&&arena_->can_publish(transport),"native collection requires fresh actual transport token");
  const auto& audit=transport.audit();need(audit.at("kind")=="native_transport"&&audit.at("state_prefix")==arena_->state_prefix()&&audit.at("dt_s")==e.dt_s&&audit.at("start_time_s")==e.physical_stage_s,"native collection clock/prefix");
  need(!audit.at("species_sink_receipts").empty(),"no actual native collection receipts");
  // Bounded one-block owner: all new sinks must belong here. Root must provide
  // a multi-owner splitter before using this primitive for other sinks/species.
  for(const auto& r:audit.at("species_sink_receipts")) {
   const std::string receipt=r.at("receipt_id"),id=r.at("packet_id");
   need(!accepted_->consumed_sinks.count(receipt)&&r.at("state")=="SS304_CONGRUENT_NEUTRAL"&&r.at("sink")=="upper_foil"&&r.at("reference_per_particle_J")==reference_per_atom_&&r.at("internal_per_particle_J")==cfg_.internal_per_atom_J,"foreign/replayed/mismatched native sink");
   need(arena_->view(id).mass_kg==atom_mass_&&arena_->view(id).charge_C==0,"actual native sink species binding");
   m::MetalSurfacePacket op;op.emit=false;op.id=id;op.count={r.at("count").at(0),r.at("count").at(1)};op.pose.position_m=r.at("position_m").get<m::PicVector>();op.pose.proper_velocity_m_s=r.at("proper_velocity_m_s").get<m::PicVector>();op.pose.mass_kg=atom_mass_;op.pose.charge_C=0;
   e.packets.push_back(op);e.native_sink_receipts.push_back(receipt);
  }
  return prepare_impl(e,fault,&transport);
 }
 private:
 Prepared prepare_impl(const Event& e,int fault,m::PreparedArena* transport) {
  need(e.base_revision==accepted_->revision&&e.state_prefix==prefix()&&!e.id.empty()&&!accepted_->events.count(e.id)&&std::isfinite(e.other_heat_J)&&std::isfinite(e.dt_s)&&e.dt_s>0,"stale/invalid metal event");
  need(e.field_work_J==0&&e.extra_latent_J==0&&(e.other_heat_J==0||!e.heat_owner.empty()),"duplicate field/latent work or missing canonical solid heat owner");
  need(accepted_->revision<std::numeric_limits<std::uint64_t>::max()&&generation_<std::numeric_limits<std::uint64_t>::max(),"metal revision/generation exhaustion");
  const double oldT=temperature();auto next=std::make_shared<Snapshot>(*accepted_);D outward,K,I,R;std::array<D,3> momentum;J packet_audit=J::array();
  m::MetalSurfaceEvent surface;surface.id=e.id;surface.state="SS304_CONGRUENT_NEUTRAL";surface.state_prefix=arena_->state_prefix();surface.material_seal_sha256=cfg_.material_seal;surface.atom_mass_kg=atom_mass_;surface.reference_per_atom_J=reference_per_atom_;surface.internal_per_atom_J=cfg_.internal_per_atom_J;surface.physical_stage_s=e.physical_stage_s;surface.thermal_foil={cfg_.geometry.x0,cfg_.geometry.x0+cfg_.geometry.width,cfg_.geometry.y0,cfg_.geometry.y0+cfg_.geometry.thickness};surface.thermal_depth_m=cfg_.geometry.height;surface.thermal_short_z_m=cfg_.axial_center_m-cfg_.geometry.height/2;surface.thermal_feed_z_m=cfg_.axial_center_m+cfg_.geometry.height/2;surface.manufactured_internal_state=cfg_.manufactured_internal_state;surface.packets=e.packets;
  std::set<std::string> operation_ids;
  for(const auto& op:e.packets) {
   need(operation_ids.insert(op.id).second,"duplicate metal packet");auto weights=m::PicArena::positive_carriers(op.count);
   auto pose=(op.emit||transport)?op.pose:arena_->view(op.id);
   const double tol=64*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(cfg_.geometry.y0));
   need(pose.position_m[0]>=cfg_.geometry.x0&&pose.position_m[0]<=cfg_.geometry.x0+cfg_.geometry.width&&std::abs(pose.position_m[1]-(cfg_.geometry.y0+cfg_.geometry.thickness))<=tol&&std::abs(pose.position_m[2]-cfg_.axial_center_m)<=cfg_.geometry.height/2,"metal source/collection not on qualified upper face");
   need(op.emit?pose.proper_velocity_m_s[1]>=0:pose.proper_velocity_m_s[1]<=0,"wrong surface momentum orientation");
   if(op.emit)need(oldT>=1800,"vapor emission below source range");
   if(!op.emit)need(m::PicArena::exact_count_equal(op.count,arena_->snapshot()->packets.at(op.id).count),"partial condensation is not implemented");
   for(double weight:weights) {
    need(std::isfinite(atom_mass_*weight)&&atom_mass_*weight>0&&std::isfinite(reference_per_atom_*weight)&&reference_per_atom_*weight!=0&&(cfg_.internal_per_atom_J==0||(std::isfinite(cfg_.internal_per_atom_J*weight)&&cfg_.internal_per_atom_J*weight>0)),"weighted particle mass/state energy not representable by native receiver");
    if(op.emit)next->solid_atoms=m::PicArena::debit(next->solid_atoms,weight);else next->solid_atoms=credit(next->solid_atoms,weight);
    const D sign(op.emit?1.:-1.);outward=outward+sign*D(weight);
    pose.mass_kg=atom_mass_;pose.weight=1;pose.charge_C=0;
    K=K+sign*D(weight)*D(m::pic_kinetic_energy_J(pose));
    I=I+sign*D(weight)*D(cfg_.internal_per_atom_J);R=R+sign*D(weight)*D(reference_per_atom_);
    for(int d=0;d<3;++d)momentum[d]=momentum[d]+sign*D(weight)*D(atom_mass_)*D(pose.proper_velocity_m_s[d]);
   }
   packet_audit.push_back({{"id",op.id},{"emit",op.emit},{"weights",weights},{"position",pose.position_m},{"proper_velocity",pose.proper_velocity_m_s},{"native_per_atom_kinetic_J",m::pic_kinetic_energy_J(pose)}});
  }
  const D packets=K+I+R;next->energy=accepted_->energy+D(e.other_heat_J)-packets;
  for(auto& id:e.native_sink_receipts)need(next->consumed_sinks.insert(id).second,"duplicate native sink receipt");
  for(int d=0;d<3;++d)next->support_impulse[d]=next->support_impulse[d]-momentum[d];
  const double target=temperature(*next);next->events.insert(e.id);++next->revision;
  // Native heat owns its normal temporary Trial. Accepted metal state remains
  // the one canonical count+exact-energy snapshot above, never a second pool.
  const D new_mass=mass(next->solid_atoms);
  const double kg=new_mass.value();
  constexpr double capacity_bound=ss304::Model::cp_liquid+ss304::Model::latent_fusion/(ss304::Model::t_liquidus-ss304::Model::t_solidus);
  const double native_activity_scale=kg*capacity_bound*3000/(4*cfg_.geometry.height);
  need(kg>0&&std::isfinite(native_activity_scale)&&native_activity_scale<=std::sqrt(std::numeric_limits<double>::max())/4,"native heat activity norm is not representable");
  auto cache=std::make_unique<h::Coupler>(cfg_.geometry,law(kg),std::vector<h::Patch>{},oldT,"metal-transfer/"+signature_,e.physical_stage_s);
  const D correction=next->energy-new_mass*D(material_.evaluate(oldT).u_model_j_kg);
  const double load=correction.value();std::vector<h::Transfer> transfers;
  int sequence=0;
  for(double x:{cfg_.geometry.x0,cfg_.geometry.x0+cfg_.geometry.width})for(double y:{cfg_.geometry.y0,cfg_.geometry.y0+cfg_.geometry.thickness}) {
   h::Transfer t;t.producer="metal-thermal-native-proof";t.energy_owner="metal-canonical-energy";t.sequence=t.owner_sequence=++sequence;t.state_id="SS304_CONGRUENT_NEUTRAL";t.parent_id=e.id;t.channel_id="conditional-metal-surface";t.event_time_s=e.physical_stage_s;t.x=x;t.y=y;t.authored_heat_J=load/4;transfers.push_back(t);
  }
  mfem::ConstantCoefficient zero(0);auto ht=cache->prepare(cache->make_epoch(e.dt_s,"no-metal-field-work"),zero,transfers,{});
  const double average=(ht.state().T[0]+ht.state().T[1]+ht.state().T[2]+ht.state().T[3])/4;
  const D native_projection=new_mass*D(material_.evaluate(average).u_model_j_kg)-next->energy;
  const double native_bound=ht.budget().roundoff_bound_J+128*std::numeric_limits<double>::epsilon()*std::max(std::abs(next->energy.value()),kg*material_.max_u());
  need(std::isfinite(native_bound)&&native_bound>=0,"nonfinite native heat proof bound");
  need(std::abs(native_projection.ld())<=native_bound,"native heat/canonical mass energy mismatch");
  for(double t:ht.state().T)need(std::abs(t-average)<=1e-8,"uniform thermal block lost symmetry");
  auto gas=[&]() -> m::PreparedArena {
   if(!transport)return arena_->prepare_metal_surface(surface,fault<=3?fault:0);
   // Validate the same sealed registry/geometry against the actual base. This
   // empty binding trial is discarded, never published or counted as a source.
   auto binding=surface;binding.packets.clear();auto admission=arena_->prepare_metal_surface(binding);
   need(arena_->can_publish(admission)&&arena_->can_publish(*transport),"native receipt staging changed base");
   return std::move(*transport);
  }();
  J impulse=J::array();for(auto& x:momentum)impulse.push_back(x.json());
  next->journal.push_back({{"event",e.id},{"atom_delta",outward.json()},{"K",K.json()},{"I",I.json()},{"R",R.json()},{"other_heat",e.other_heat_J},{"heat_owner",e.heat_owner},{"gas_momentum",impulse},{"native_sink_receipts",e.native_sink_receipts},{"arena_stage_sha256",hash(gas.audit().dump())}});
  Prepared p(std::move(gas));p.owner_=this;p.generation_=generation_;p.before_=accepted_;p.after_=next;p.heat_cache_=std::move(cache);p.heat_trial_=std::move(ht);
  p.audit_={{"event",e.id},{"material_seal",cfg_.material_seal},{"signature",signature_},{"atom_mass_kg",atom_mass_},{"reference_per_atom_J",reference_per_atom_},{"internal_per_atom_J",cfg_.internal_per_atom_J},{"actual_packets",packet_audit},{"other_heat_J",e.other_heat_J},{"solid_before",state_json(*accepted_)},{"solid_after",state_json(*next)},{"outward_atoms",outward.json()},{"actual_native_kinetic_J",K.json()},{"actual_state_internal_J",I.json()},{"actual_state_reference_J",R.json()},{"outward_packet_energy_J",packets.json()},{"gas_momentum",impulse},{"fixed_support_work_J",0},{"field_work_J",0},{"extra_latent_debit_J",0},{"native_heat_load_J",load},{"native_heat_projection_residual",native_projection.json()},{"native_heat_roundoff_bound_J",native_bound},{"native_heat_temperature_K",average},{"canonical_temperature_K",target},{"gas_stage",p.gas_.audit()}};
  p.audit_["canonical_before_temperature_K"]=oldT;p.audit_["u_before_J_kg"]=material_.evaluate(oldT).u_model_j_kg;p.audit_["native_trial_current_mass_kg"]=kg;
  p.audit_["consumed_native_sink_receipts"]=e.native_sink_receipts;p.audit_["other_heat_owner"]=e.heat_owner;
  if(fault==4)throw std::runtime_error("fault after both native stages");
  return p;
 }
 public:
 bool can_publish(const Prepared& p)const noexcept{return generation_<std::numeric_limits<std::uint64_t>::max()&&p.owner_==this&&p.status_==Prepared::Status::Staged&&p.generation_==generation_&&accepted_==p.before_&&p.heat_cache_->preflight(p.heat_trial_)&&arena_->can_publish(p.gas_);}
 bool try_publish(Prepared& p,bool fail_after_gas=false)noexcept{
  if(!can_publish(p)||!arena_->try_publish(p.gas_))return false;
  if(fail_after_gas){if(!arena_->try_rollback(p.gas_))std::terminate();p.status_=Prepared::Status::RolledBack;++generation_;return false;}
  accepted_=p.after_;p.status_=Prepared::Status::Published;++generation_;return true;
 }
 bool try_rollback(Prepared& p)noexcept{
  if(p.owner_!=this||p.status_!=Prepared::Status::Published||accepted_!=p.after_||!arena_->try_rollback(p.gas_))return false;
  accepted_=p.before_;p.status_=Prepared::Status::RolledBack;if(generation_<std::numeric_limits<std::uint64_t>::max())++generation_;return true;
 }
 J checkpoint()const{J body={{"schema","uniform-native-metal-transfer-v1"},{"signature",signature_},{"arena_prefix",arena_->state_prefix()},{"state",state_json(*accepted_)}};return {{"body",body},{"sha256",hash(body.dump())}};}
 void restore(const J& record){need(accepted_->revision==0&&generation_==0&&record.size()==2&&record.at("sha256")==hash(record.at("body").dump()),"metal restart digest/freshness");const auto& b=record.at("body");need(b.size()==4&&b.at("schema")=="uniform-native-metal-transfer-v1"&&b.at("signature")==signature_&&b.at("arena_prefix")==arena_->state_prefix(),"metal restart source/arena binding");const auto& j=b.at("state");need(j.size()==7,"metal state schema");auto n=std::make_shared<Snapshot>();n->solid_atoms={j.at("solid_atoms").at(0),j.at("solid_atoms").at(1)};need(h::canonical(n->solid_atoms),"metal restart count");n->energy=D::read(j.at("energy"));for(int d=0;d<3;++d)n->support_impulse[d]=D::read(j.at("support_impulse").at(d));n->events=j.at("events").get<std::set<std::string>>();n->consumed_sinks=j.at("consumed_sinks").get<std::set<std::string>>();n->journal=j.at("journal");n->revision=j.at("revision");need(n->events.size()==n->revision&&n->journal.is_array()&&n->journal.size()==n->revision,"metal restart event count");
  std::map<std::string,J> known_stages;std::set<std::string>used_stages,events,sinks;for(auto row:arena_->snapshot()->journal){J a=row.at("payload");a.erase("previous_chain");known_stages.emplace(hash(a.dump()),a);}
  D atoms=count(accepted_->solid_atoms),energy=accepted_->energy;std::array<D,3>support;
  for(const auto& row:n->journal){std::string id=row.at("event"),stage=row.at("arena_stage_sha256");need(events.insert(id).second&&known_stages.count(stage)&&used_stages.insert(stage).second,"metal restart arena-stage journal binding");D delta=D::read(row.at("atom_delta")),k=D::read(row.at("K")),i=D::read(row.at("I")),r=D::read(row.at("R"));need(i.compare(delta*D(cfg_.internal_per_atom_J))==0&&r.compare(delta*D(reference_per_atom_))==0,"metal restart phase/internal reference ledger");
   const auto& actual=known_stages.at(stage);J packets=J::array();std::set<std::string>stage_sinks;
   if(actual.at("kind")=="metal_surface"){need(actual.at("event_id")==id&&actual.at("material_seal_sha256")==cfg_.material_seal,"metal source restart identity");packets=actual.at("packets");}
   else {need(actual.at("kind")=="native_transport","foreign arena stage");for(const auto& s:actual.at("species_sink_receipts")){need(s.at("state")=="SS304_CONGRUENT_NEUTRAL"&&s.at("sink")=="upper_foil"&&s.at("internal_per_particle_J")==cfg_.internal_per_atom_J&&s.at("reference_per_particle_J")==reference_per_atom_,"foreign native sink during restart");J w=J::array();for(const auto& carrier:s.at("native_carrier_receipts"))w.push_back(carrier.at("weight"));packets.push_back({{"emit",false},{"weights",w},{"proper_velocity",s.at("proper_velocity_m_s")}});stage_sinks.insert(s.at("receipt_id").get<std::string>());}}
   D expected_delta,expected_k;std::array<D,3>expected_momentum;
   for(const auto& p:packets){m::GapParticle pose;pose.mass_kg=atom_mass_;pose.weight=1;pose.proper_velocity_m_s=p.at("proper_velocity").get<m::PicVector>();const D sign(p.at("emit").get<bool>()?1.:-1.);for(double w:p.at("weights").get<std::vector<double>>()){expected_delta=expected_delta+sign*D(w);expected_k=expected_k+sign*D(w)*D(m::pic_kinetic_energy_J(pose));for(int d=0;d<3;++d)expected_momentum[d]=expected_momentum[d]+sign*D(w)*D(atom_mass_)*D(pose.proper_velocity_m_s[d]);}}
   need(expected_delta.compare(delta)==0&&expected_k.compare(k)==0,"metal restart native count/kinetic receipt mismatch");for(int d=0;d<3;++d)need(expected_momentum[d].compare(D::read(row.at("gas_momentum").at(d)))==0,"metal restart native momentum mismatch");
   std::set<std::string>declared=row.at("native_sink_receipts").get<std::set<std::string>>();need(declared==stage_sinks,"metal restart consumed native sink set mismatch");
   atoms=atoms-delta;energy=energy+D(row.at("other_heat").get<double>())-k-i-r;need(atoms.compare(D(0))>0&&energy.compare(atoms*D(atom_mass_)*D(material_.min_u()))>=0&&energy.compare(atoms*D(atom_mass_)*D(material_.max_u()))<=0,"metal restart journal domain");for(int d=0;d<3;++d)support[d]=support[d]-D::read(row.at("gas_momentum").at(d));for(const auto& x:row.at("native_sink_receipts"))need(sinks.insert(x.get<std::string>()).second,"metal restart sink replay");}
  need(atoms.compare(count(n->solid_atoms))==0&&energy.compare(n->energy)==0&&events==n->events&&sinks==n->consumed_sinks,"metal restart canonical count/energy/event ledger mismatch");for(int d=0;d<3;++d)need(support[d].compare(n->support_impulse[d])==0,"metal restart support momentum ledger mismatch");
  (void)temperature(*n);accepted_=n;++generation_;}
};
} // namespace metal_transfer
