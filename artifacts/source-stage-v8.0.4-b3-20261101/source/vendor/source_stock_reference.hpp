#pragma once
// Immutable chemical-reference transfer. Source owns stock; native owns births.
// This component never modifies a source/thermal/native accepted state.
#include "raikou/mfem_pic_full_time.h"
#include "hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "pic-electrode-heat-source-v7.4-mass-cache-20261101/native/mfem_electrode_heat_source_v7/source_owner.hpp"
#include "pic-electrode-heat-source-v7.4-mass-cache-20261101/adapters/native_source_v7_mass.hpp"
#include <map>
#include <set>
#include <optional>
#include <tuple>
namespace source_stock_reference_v2 {
namespace m=raikou::mfem_backend;namespace s=electrode_heat_v7;namespace eh=electrode_heat_v2;namespace se=surface_emission;
using V=m::PicVector;using Count=m::ArenaCount;using pfc::require;
inline bool bits(double a,double b){return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);}
inline bool bits(const V&a,const V&b){for(int k=0;k<3;++k)if(!bits(a[k],b[k]))return false;return true;}
inline void number(pfc::Sha256&h,double x){h.pod(std::bit_cast<std::uint64_t>(x));}
inline void vector(pfc::Sha256&h,const V&v){for(double x:v)number(h,x);}
inline void count(pfc::Sha256&h,const Count&n){auto w=n.positive_carriers();h.pod(std::uint64_t(w.size()));for(double x:w)number(h,x);}
inline pfc::ExactSum reference(const Count&n,double R){pfc::ExactSum x;for(double w:n.positive_carriers())x.add_product(w,R);return x;}
struct PatchRegistration {
 std::size_t patch=0;
 std::optional<std::uint32_t> neutral_stream,electron_stream;
 bool upper_foil=true;
 double z_center_m=0;
};
// Explicit coordinate-domain wrappers: neutral packet coordinates are local to
// the thermal slice; exposure origins already use the global field frame.
struct LocalThermalPoint {V value;};
struct GlobalFieldPoint {V value;};
struct NativePoint {V value;};
struct NeutralLink {std::size_t source_packet=0;std::string native_packet_id;};
struct ElectronLink {std::uint64_t candidate_native_id=0;std::string native_packet_id;};
struct StockEntry {std::size_t patch=0;std::string patch_id,state;Count inventory;double reference_J=0;std::string exact_reference_J;};
class StockSnapshot {
 friend class Adapter;
 std::string configuration_,source_state_,identity_;std::vector<StockEntry>stock_;
public:
 const std::string&sha256()const noexcept{return identity_;}
 const std::string&source_state_sha256()const noexcept{return source_state_;}
 const std::vector<StockEntry>&entries()const noexcept{return stock_;}
};
struct SurfaceTransfer {
 std::size_t source_port=0,patch=0;std::uint64_t native_id=0;std::string source_token,packet_id;
 V field_origin{},declared_displacement{},native_origin{},native_position{};
 double charge_C=0,wall_heat_J=0,electron_reference_J=0,matching_work_J=0;
 bool bridge=false;
};
struct PatchTransfer {
 std::size_t source_port=0,patch=0;std::string patch_id,state;Count before,after,birth_count;
 double reference_J=0;
 std::string before_reference_J,after_reference_J,gas_reference_J,closure_J;
};
class PreparedReferenceTransfer {
 friend class Adapter;
 struct Data {
  m::FullTimeBirthCredit credit;
  m::PreparedBirths native;
  std::vector<PatchTransfer> transfers;std::vector<SurfaceTransfer>surface;
  std::string configuration,source_before,source_trial,credit_sha,identity,source_load_sha;
  std::size_t neutral_births=0,electron_births=0;
  Data(m::FullTimeBirthCredit c,m::PreparedBirths n):credit(std::move(c)),native(std::move(n)){}
 };
 std::unique_ptr<const Data>d_;
 explicit PreparedReferenceTransfer(std::unique_ptr<const Data>d):d_(std::move(d)){}
public:
 PreparedReferenceTransfer(PreparedReferenceTransfer&&)noexcept=default;
 PreparedReferenceTransfer&operator=(PreparedReferenceTransfer&&)noexcept=default;
 PreparedReferenceTransfer(const PreparedReferenceTransfer&)=delete;
 PreparedReferenceTransfer&operator=(const PreparedReferenceTransfer&)=delete;
 // Borrowed immutable views expire on move/destruction. Native retains its own
 // original outer-step lifetime/phase validation; this wrapper grants no publish.
 const m::PreparedBirths&native_births()const{return d_->native;}
 const m::FullTimeBirthCredit&credit()const{return d_->credit;}
 const std::vector<PatchTransfer>&transfers()const{return d_->transfers;}
 const std::vector<SurfaceTransfer>&surface_transfers()const{return d_->surface;}
 const std::string&sha256()const{return d_->identity;}
 const std::string&credit_sha256()const{return d_->credit_sha;}
 const std::string&source_load_sha256()const{return d_->source_load_sha;}
 std::size_t neutral_births()const{return d_->neutral_births;}
 std::size_t electron_births()const{return d_->electron_births;}
};
class Adapter {
 struct PortValidation {std::string configuration,source_before,source_trial,credit_sha,load_sha;std::vector<PatchTransfer>transfers;std::vector<SurfaceTransfer>surface;};
 m::GapPicOptions options_;m::FullTimeChemistryRegistry registry_;std::vector<m::ArenaStateSpec>states_;
 se::Config mass_policy_;std::vector<PatchRegistration>patches_;std::string source_configuration_,configuration_;
 const m::ArenaStateSpec&species(const std::string&id)const{auto p=std::find_if(states_.begin(),states_.end(),[&](const auto&x){return x.id==id;});require(p!=states_.end(),"source reference: unregistered source species");return *p;}
 const PatchRegistration&registration(std::size_t j)const{auto p=std::find_if(patches_.begin(),patches_.end(),[&](const auto&x){return x.patch==j;});require(p!=patches_.end(),"source reference: unmapped patch");return *p;}
 std::size_t patch_index(const s::Owner&o,const std::string&id)const{for(std::size_t j=0;j<o.source_patches();++j)if(o.source_patch(j).id==id)return j;throw std::invalid_argument("source reference: packet has unknown patch");}
 void admission(const s::Owner&o)const{
  require(o.configuration_sha256()==source_configuration_,"source reference: foreign source configuration");
  require(!o.options().phase&&o.options().water_writeoff_T_K==0,"source reference: phase/metal/writeoff mapping unsupported");
  require(patches_.size()==o.source_patches(),"source reference: complete patch registration required");
  for(std::size_t j=0;j<o.source_patches();++j){const auto p=o.source_patch(j);const auto st=o.source_settings(j);const auto&r=registration(j);
   require(!st.surface.wall.enabled,"source reference: every enabled wall path requires a separate species/reference mapping");require(!st.surface.photo.enabled,"source reference: photoemission mapping is not qualified");
   require(bits(st.spectrum.electron_mass_kg,mass_policy_.electron_mass_kg)&&st.spectrum.electron_mass_identity==mass_policy_.electron_mass_identity,"source reference: electron mass/provenance differs from pinned native registry");
   if(p.desorption_enabled){require(r.neutral_stream.has_value(),"source reference: enabled desorption patch has no neutral stream");const auto&sp=species(p.state_id);
    require(sp.charge_e==0&&bits(sp.internal_J,0.)&&bits(p.internal_eV,0.)&&bits(p.mass_kg,sp.mass_kg),"source reference: unsupported charged/internal source or mismatched kinetic mass");
    auto q=std::find_if(registry_.species.begin(),registry_.species.end(),[&](const auto&x){return x.state==p.state_id;});require(q!=registry_.species.end(),"source reference: missing blob species");
    for(std::size_t k=2;k<q->nuclei.size();++k)require(q->nuclei[k]==0,"source reference: metal/other-nucleus source mapping unsupported");
   }else require(o.inventory(j).zero()&&!r.neutral_stream,"source reference: disabled desorption patch must have no stock or neutral stream");
  }
 }
 static std::string loads_sha(const s::Trial&t){s::detail::Encoding e;e.s("source-reference/unmodified-thermal-loads-v1");e.u(t.loads().size());for(const auto&x:t.loads())s::detail::transfer(e,x);return e.sha();}
 static std::string credit_sha(const m::FullTimeBirthCredit&c){
  pfc::Sha256 h;h.text("source-reference/full-native-credit-v1");h.pod(c.arrival_step);h.pod(c.base_native_id);
  for(const auto*p:{&c.plan_sha256,&c.staged_token_sha256,&c.source_before_sha256,&c.source_after_sha256})h.text(*p);
  h.pod(std::uint64_t(c.streams.size()));for(const auto&x:c.streams){h.pod(x.stream_id);h.text(x.state);count(h,x.source_debit);}
  h.pod(std::uint64_t(c.births.size()));for(const auto&b:c.births){h.pod(b.step);h.pod(b.stream_id);h.pod(b.index);h.text(b.packet_id);h.text(b.state);h.text(b.origin);number(h,b.weight);vector(h,b.position_m);vector(h,b.proper_velocity_m_s);number(h,b.birth_time_s);number(h,b.remaining_flight_s);h.pod(b.native_id);h.pod(std::uint8_t(b.bridge));vector(h,b.bridge_start_m);h.pod(std::uint8_t(b.returned));number(h,b.internal_energy_J);}
  // Unsupported wall records are refused, never omitted from a supported digest.
  require(c.wall.empty(),"source reference: native wall births unsupported");h.pod(std::uint64_t(0));return h.hex();
 }
public:
 Adapter(const s::Owner&o,const m::GapPicOptions&options,const m::FullTimeChemistryConfiguration&chem,std::vector<PatchRegistration>patches):options_(options),registry_(m::full_time_chemistry_registry(chem,options)),states_(m::full_time_chemistry_states(registry_)),patches_(std::move(patches)),source_configuration_(o.configuration_sha256()){
  require(chem.blob_sha256==native_source_v7::supported_blob_sha256,"source reference: unsupported chemistry blob");
  mass_policy_=native_source_v7::electron_mass_configuration(se::Config{4.3},chem.blob_path,chem.blob_sha256,registry_,states_);
  require(options.finite_foils&&options.finite_axial_extent,"source reference: finite foil geometry required");
  std::set<std::size_t>index;std::set<std::uint32_t>streams;for(const auto&r:patches_){require(r.patch<o.source_patches()&&index.insert(r.patch).second&&std::isfinite(r.z_center_m),"source reference: duplicate/invalid patch registration");for(const auto*x:{&r.neutral_stream,&r.electron_stream})if(*x)require(streams.insert(**x).second,"source reference: duplicate stream registration");}
  admission(o);pfc::Sha256 h;h.text("source-reference/registry-stock-surface-adapter-v2");h.text(registry_.blob_sha256);h.text(registry_.registry_sha256);h.text(registry_.configuration_sha256);h.text(source_configuration_);
  for(double x:{options.width_m,options.gap_m,options.depth_m,options.short_z_m,options.feed_z_m,options.lower_foil.left_m,options.lower_foil.right_m,options.lower_foil.bottom_m,options.lower_foil.top_m,options.upper_foil.left_m,options.upper_foil.right_m,options.upper_foil.bottom_m,options.upper_foil.top_m})number(h,x);
  h.text("native-face-translation-v2;typed-local-neutral/global-exposure;normal-snap-tolerance=1e-12m");h.pod(std::uint64_t(patches_.size()));for(const auto&r:patches_){h.pod(std::uint64_t(r.patch));h.pod(std::uint8_t(r.upper_foil));number(h,r.z_center_m);for(const auto*x:{&r.neutral_stream,&r.electron_stream}){h.pod(std::uint8_t(x->has_value()));if(*x)h.pod(**x);}}
  configuration_=h.hex();
 }
 const std::string&configuration_sha256()const noexcept{return configuration_;}
 const m::FullTimeChemistryRegistry&registry()const noexcept{return registry_;}
private:
 V native_face_point(const s::Owner&o,std::size_t patch,const V&x)const{
  const auto&r=registration(patch);const auto p=o.source_patch(patch);const auto&f=r.upper_foil?options_.upper_foil:options_.lower_foil;
  V v{x[0]+.5*options_.width_m,x[1]+.5*options_.gap_m,r.z_center_m+x[2]};int axis=1;double plane=0;
  switch(p.face){case eh::Face::LowY:plane=f.bottom_m;break;case eh::Face::HighY:plane=f.top_m;break;case eh::Face::LowX:axis=0;plane=f.left_m;break;case eh::Face::HighX:axis=0;plane=f.right_m;break;}
  require(std::isfinite(v[axis])&&std::abs(v[axis]-plane)<=1e-12,"source reference: field point does not map onto registered native foil face");v[axis]=plane;return v;
 }
public:
 V native_neutral_face_point(const s::Owner&o,std::size_t patch,LocalThermalPoint x)const{return native_face_point(o,patch,x.value);}
 V native_surface_point(const s::Owner&o,std::size_t patch,GlobalFieldPoint x)const{
  const double global_z=x.value[2];x.value[2]-=registration(patch).z_center_m;auto point=native_face_point(o,patch,x.value);point[2]=global_z;return point;
 }
 GlobalFieldPoint field_point(NativePoint x)const{return {{x.value[0]-.5*options_.width_m,x.value[1]-.5*options_.gap_m,x.value[2]}};}
 StockSnapshot stock_checkpoint(const s::Owner&o)const{
  admission(o);StockSnapshot a;a.configuration_=configuration_;a.source_state_=o.state_sha256();pfc::Sha256 h;h.text("source-reference/derived-stock-checkpoint-v1");h.text(a.configuration_);h.text(a.source_state_);
  for(std::size_t j=0;j<o.source_patches();++j){const auto p=o.source_patch(j);if(!p.desorption_enabled)continue;const auto&sp=species(p.state_id);StockEntry x{j,p.id,p.state_id,o.inventory(j),sp.reference_J,reference(o.inventory(j),sp.reference_J).hex()};h.pod(std::uint64_t(j));h.text(x.patch_id);h.text(x.state);count(h,x.inventory);number(h,x.reference_J);h.text(x.exact_reference_J);a.stock_.push_back(std::move(x));}a.identity_=h.hex();return a;
 }
 bool restored_stock_matches(const s::Owner&o,const StockSnapshot&old)const{return stock_checkpoint(o).sha256()==old.sha256();}
 private:
 PortValidation inspect(const s::Owner&o,const s::Trial&t,const std::shared_ptr<const s::Field>&field,const s::FinalFieldView&actual,const m::FullTimeBirthCredit&credit,std::span<const NeutralLink>neutral,std::span<const ElectronLink>electrons)const{
  admission(o);require(field&&field->matches(actual)&&o.preflight(t,actual),"source reference: stale/foreign source trial or field");
  require(!t.wall()&&t.water_writeoff().empty()&&t.phase_reports().empty()&&credit.wall.empty(),"source reference: unsupported wall/phase/writeoff flow");
  const auto before=o.state_sha256(),trial_sha=t.content_sha256(),load_sha=loads_sha(t);
  require(credit.births.size()==neutral.size()+electrons.size()&&neutral.size()==t.packets().size(),"source reference: full credit/source packet partition is incomplete");
  std::map<std::string,const m::FullTimeBirth*>born;std::map<std::uint32_t,Count>stream_counts;for(const auto&b:credit.births){require(born.emplace(b.packet_id,&b).second,"source reference: duplicate native birth");require(!b.returned&&bits(b.internal_energy_J,0.),"source reference: returned/internal paths require qualified work mapping");}
  std::map<std::uint32_t,std::string>streams;for(const auto&r:patches_){const auto p=o.source_patch(r.patch);if(r.neutral_stream)streams.emplace(*r.neutral_stream,p.state_id);if(r.electron_stream)streams.emplace(*r.electron_stream,"pfc-v3/e");}
  require(credit.streams.size()==streams.size(),"source reference: native stream declaration differs from complete source mapping");
  std::set<std::uint32_t>seen_streams;for(const auto&x:credit.streams)require(streams.contains(x.stream_id)&&streams.at(x.stream_id)==x.state&&seen_streams.insert(x.stream_id).second,"source reference: duplicate/unmapped native stream");
  std::set<std::string>used_births;std::set<std::size_t>used_packets;std::vector<Count>neutral_counts(o.source_patches());std::vector<pfc::ExactSum>gas_ref(o.source_patches());
  for(const auto&link:neutral){require(link.source_packet<t.packets().size()&&used_packets.insert(link.source_packet).second&&born.contains(link.native_packet_id)&&used_births.insert(link.native_packet_id).second,"source reference: missing/duplicate neutral packet mapping");const auto&p=t.packets()[link.source_packet];const auto&b=*born.at(link.native_packet_id);auto j=patch_index(o,p.patch);const auto&r=registration(j);const auto&sp=species(p.state_id);
   require(!b.bridge&&r.neutral_stream&&b.stream_id==*r.neutral_stream&&b.state==p.state_id&&b.origin==p.event_id&&bits(b.weight,p.weight)&&bits(b.birth_time_s,p.time)&&bits(b.position_m,native_neutral_face_point(o,j,LocalThermalPoint{p.x}))&&bits(b.proper_velocity_m_s,p.proper_velocity_m_s),"source reference: native neutral does not match actual sealed source packet");
   require(bits(p.mass_kg,sp.mass_kg)&&bits(p.charge_C,0.)&&bits(p.internal_J,0.)&&bits(p.reference_J,0.)&&sp.charge_e==0&&bits(sp.internal_J,0.),"source reference: neutral mass/charge/internal/reference convention mismatch");
   require(bits(b.remaining_flight_s,field->frame().end_s-b.birth_time_s),"source reference: neutral birth flight mismatch");
   neutral_counts[j]=neutral_counts[j].credited(p.weight);stream_counts[b.stream_id]=stream_counts[b.stream_id].credited(p.weight);gas_ref[j].add_product(p.weight,sp.reference_J);
  }
  // Every electron comes from a genuine candidate AND the actual handoff stored
  // in this trial. Restrict this first adapter to the zero-displacement thermal
  // surface handoff; it does not silently drop active/returned/bridge currents.
  require(electrons.size()==t.source_candidates().size()+t.surface_births().size()&&t.source_candidates().size()==t.source_handoffs().size(),"source reference: source electron candidate/handoff partition incomplete");
  std::set<std::uint64_t>used_candidates;std::vector<SurfaceTransfer>surface;
  for(const auto&link:electrons){require(used_candidates.insert(link.candidate_native_id).second&&born.contains(link.native_packet_id)&&used_births.insert(link.native_packet_id).second,"source reference: missing/duplicate electron mapping");
   const auto surface_it=std::find_if(t.surface_births().begin(),t.surface_births().end(),[&](const auto&x){return x.canonical_id()==link.candidate_native_id;});
   if(surface_it!=t.surface_births().end()){
    const auto&x=*surface_it;require(o.preflight(x,field,actual),"source reference: stale/foreign instantaneous surface token");
    const auto&p=x.phase().packet.particle;const auto&b=*born.at(link.native_packet_id);const auto&r=registration(x.patch_index());const auto&sp=species("pfc-v3/e");
    require(x.heat().state_id==sp.id&&x.phase().outcome==s::BirthOutcome::Handoff&&bits(x.policy().return_extent_m,0.)&&r.electron_stream&&b.stream_id==*r.electron_stream&&b.state==sp.id&&b.native_id==x.canonical_id()&&bits(b.weight,p.weight)&&bits(p.electron_mass_kg,sp.mass_kg)&&bits(sp.internal_J,0.)&&bits(sp.reference_J,0.)&&sp.charge_e==-1,"source reference: surface birth ID/count/mass/reference/return mismatch");
    V face=native_surface_point(o,x.patch_index(),GlobalFieldPoint{x.origin_m()}),point=face;for(int k=0;k<3;++k)point[k]+=x.phase().displacement_m[k];
    const bool bridge=x.phase().packet.matching_distance_m>0;
    require(b.origin=="source-v7.3/surface/"+x.content_sha256()&&bits(b.position_m,point)&&bits(b.proper_velocity_m_s,p.proper_velocity_m_s)&&bits(b.birth_time_s,x.birth_time_s())&&bits(b.remaining_flight_s,field->frame().end_s-x.birth_time_s())&&b.bridge==bridge&&(!bridge||bits(b.bridge_start_m,face)),"source reference: native surface birth or bridge differs from actual sealed matching plane");
    if(!bridge)require(b.bridge_start_m==V{},"source reference: unused bridge start must be zero");
    surface.push_back({0,x.patch_index(),x.canonical_id(),x.content_sha256(),b.packet_id,x.origin_m(),x.phase().displacement_m,face,point,-se::qe*p.weight,p.ledger.wall_heat_J,p.ledger.electron_reference_change_J,p.ledger.resolved_matching_field_work_J,bridge});
    stream_counts[b.stream_id]=stream_counts[b.stream_id].credited(p.weight);continue;
   }
   auto c=std::find_if(t.source_candidates().begin(),t.source_candidates().end(),[&](const auto&x){return x.canonical_id()==link.candidate_native_id;});auto h=std::find_if(t.source_handoffs().begin(),t.source_handoffs().end(),[&](const auto&x){return x.particle().id==link.candidate_native_id;});
   require(c!=t.source_candidates().end()&&h!=t.source_handoffs().end()&&o.preflight(*c,field,actual),"source reference: electron is not an actual live trial source/handoff");const auto&bs=c->birth_sample();const auto&state=c->state();const auto&p=state.native.particle;const auto&b=*born.at(link.native_packet_id);const auto&r=registration(state.patch_index);const auto&sp=species("pfc-v3/e");
   require(!b.bridge&&bs&&state.native.status==se::LayerStatus::PICHandedOff&&c->paths().empty()&&state.native.position_from_surface_m==V{},"source reference: source microscopic path requires qualified field-work mapping");
   require(r.electron_stream&&b.stream_id==*r.electron_stream&&b.state==sp.id&&b.native_id==c->canonical_id()&&bits(b.weight,p.weight)&&bits(p.electron_mass_kg,sp.mass_kg)&&bits(sp.internal_J,0.)&&bits(sp.reference_J,0.)&&sp.charge_e==-1,"source reference: electron ID/count/mass/reference mismatch");
   require(b.origin=="source-v7/candidate/"+c->content_sha256()&&bits(b.position_m,native_surface_point(o,state.patch_index,GlobalFieldPoint{state.surface_origin_m}))&&bits(b.proper_velocity_m_s,p.proper_velocity_m_s)&&bits(b.birth_time_s,state.native.current_time_s)&&bits(b.remaining_flight_s,field->frame().end_s-state.native.current_time_s)&&bits(b.birth_time_s,h->handoff_time_s())&&bits(b.remaining_flight_s,h->remaining_duration_s()),"source reference: native electron birth kinematics differ from actual source handoff start");
   auto verified=o.handoff(*c,actual,h->particle());require(verified.content_sha256()==h->content_sha256(),"source reference: changed or foreign handoff");
   stream_counts[b.stream_id]=stream_counts[b.stream_id].credited(p.weight);
  }
  require(used_births.size()==credit.births.size(),"source reference: unpartitioned native birth");
  for(const auto&x:credit.streams)require(stream_counts[x.stream_id]==x.source_debit,"source reference: native debit differs from exact mapped source counts");
  std::vector<PatchTransfer>transfers;for(std::size_t j=0;j<o.source_patches();++j){const auto p=o.source_patch(j);if(!p.desorption_enabled){require(neutral_counts[j].zero()&&t.inventory(j)==o.inventory(j),"source reference: disabled source inventory changed");continue;}const auto&sp=species(p.state_id);Count remaining=o.inventory(j);for(double w:neutral_counts[j].positive_carriers())remaining=remaining.debit(w);require(remaining==t.inventory(j),"source reference: exact source stock decrease differs from gas birth count");auto old=reference(o.inventory(j),sp.reference_J),next=reference(t.inventory(j),sp.reference_J);auto closure=next;closure.sub(old);closure.add(gas_ref[j]);require(closure.is_zero(),"source reference: stock plus gas chemical-reference transfer does not close exactly");transfers.push_back({0,j,p.id,p.state_id,o.inventory(j),t.inventory(j),neutral_counts[j],sp.reference_J,old.hex(),next.hex(),gas_ref[j].hex(),closure.hex()});}
  return {configuration_,before,trial_sha,credit_sha(credit),load_sha,std::move(transfers),std::move(surface)};
 }
public:
 struct Port {
  const Adapter&adapter;const s::Owner&owner;const s::Trial&trial;
  std::shared_ptr<const s::Field>field;s::FinalFieldView actual;
  std::span<const NeutralLink>neutral;std::span<const ElectronLink>electrons;
 };
 static PreparedReferenceTransfer prepare_group(std::span<const Port>ports,m::BirthScope&scope,m::FullTimeBirthCredit credit){
  require(!ports.empty()&&credit.wall.empty(),"source reference: source ports required and wall paths unsupported");
  std::set<const s::Owner*>owners;std::set<std::tuple<std::string,std::string,std::string>>physical_patches;std::map<std::uint32_t,std::size_t>stream_owner;
  for(std::size_t j=0;j<ports.size();++j){const auto&p=ports[j];require(owners.insert(&p.owner).second,"source reference: duplicate source owner port");require(p.adapter.registry_.registry_sha256==ports[0].adapter.registry_.registry_sha256&&p.adapter.registry_.configuration_sha256==ports[0].adapter.registry_.configuration_sha256,"source reference: source ports have different native registry/configuration bindings");for(const auto&r:p.adapter.patches_){const auto&id=p.owner.patch_identity(r.patch);require(physical_patches.emplace(id.physical_surface_id,id.material_layer_id,p.owner.source_patch(r.patch).id).second,"source reference: physical patch represented by multiple source owners");for(const auto*x:{&r.neutral_stream,&r.electron_stream})if(*x)require(stream_owner.emplace(**x,j).second,"source reference: native stream belongs to multiple source ports");}}
  require(stream_owner.size()==credit.streams.size(),"source reference: full native stream partition incomplete");
  // Copy only common metadata. Copying all births into every port before
  // filtering would retain O(ports*births) capacity on a many-slice surface.
  std::vector<m::FullTimeBirthCredit>local(ports.size());
  for(auto&x:local){x.arrival_step=credit.arrival_step;x.base_native_id=credit.base_native_id;x.plan_sha256=credit.plan_sha256;x.staged_token_sha256=credit.staged_token_sha256;x.source_before_sha256=credit.source_before_sha256;x.source_after_sha256=credit.source_after_sha256;}
  std::set<std::uint32_t>declared;
  for(const auto&x:credit.streams){require(stream_owner.contains(x.stream_id)&&declared.insert(x.stream_id).second,"source reference: duplicate/unmapped full native stream");local[stream_owner.at(x.stream_id)].streams.push_back(x);}
  for(const auto&x:credit.births){require(stream_owner.contains(x.stream_id),"source reference: birth outside all source ports");local[stream_owner.at(x.stream_id)].births.push_back(x);}
  std::vector<PortValidation>checked;checked.reserve(ports.size());std::size_t neutrals=0,electrons=0;
  for(std::size_t j=0;j<ports.size();++j){const auto&p=ports[j];checked.push_back(p.adapter.inspect(p.owner,p.trial,p.field,p.actual,local[j],p.neutral,p.electrons));neutrals+=p.neutral.size();electrons+=p.electrons.size();}
  require(neutrals+electrons==credit.births.size(),"source reference: whole native birth partition incomplete");
  const auto whole_credit_sha=credit_sha(credit);
  // Exactly one native prepare, with the full owned credit. The authentic token
  // and its credit travel together; native original-scope lifetime checks remain.
  auto native=scope.prepare(credit);auto d=std::make_unique<PreparedReferenceTransfer::Data>(std::move(credit),std::move(native));d->credit_sha=whole_credit_sha;d->neutral_births=neutrals;d->electron_births=electrons;
  pfc::Sha256 hash,loads,configuration,bases,trials;hash.text("source-reference/prepared-port-set-transfer-v2");hash.text(whole_credit_sha);hash.text(d->native.receipt().row_sha256);hash.text(d->native.receipt().births_root_sha256);hash.pod(std::uint64_t(ports.size()));loads.text("source-reference/port-set-thermal-loads-v1");configuration.text("source-reference/port-set-configuration-v1");bases.text("source-reference/port-set-source-bases-v1");trials.text("source-reference/port-set-source-trials-v1");
  for(std::size_t j=0;j<ports.size();++j){const auto&p=ports[j];auto&v=checked[j];require(p.owner.state_sha256()==v.source_before&&p.owner.preflight(p.trial,p.actual)&&loads_sha(p.trial)==v.load_sha,"source reference: source changed during native preparation");
   for(const auto*x:{&v.configuration,&v.source_before,&v.source_trial,&v.credit_sha,&v.load_sha})hash.text(*x);loads.text(v.load_sha);configuration.text(v.configuration);bases.text(v.source_before);trials.text(v.source_trial);
   hash.pod(std::uint64_t(p.neutral.size()));for(const auto&l:p.neutral){hash.pod(std::uint64_t(l.source_packet));hash.text(l.native_packet_id);}hash.pod(std::uint64_t(p.electrons.size()));for(const auto&l:p.electrons){hash.pod(l.candidate_native_id);hash.text(l.native_packet_id);}
   for(auto&x:v.surface){x.source_port=j;d->surface.push_back(std::move(x));}
   hash.pod(std::uint64_t(v.transfers.size()));for(auto&x:v.transfers){x.source_port=j;hash.pod(std::uint64_t(j));hash.pod(std::uint64_t(x.patch));hash.text(x.patch_id);hash.text(x.state);count(hash,x.before);count(hash,x.after);count(hash,x.birth_count);number(hash,x.reference_J);hash.text(x.before_reference_J);hash.text(x.after_reference_J);hash.text(x.gas_reference_J);hash.text(x.closure_J);d->transfers.push_back(std::move(x));}
  }
  // Bind the actual native bridge receipt. Parent must still deposit these exact
  // field-frame paths and compare its measured E.J work to matching_work_J.
  std::size_t bridge_count=0;
  for(const auto&x:d->surface){if(!x.bridge)continue;++bridge_count;const auto b=std::find_if(d->native.receipt().bridges.begin(),d->native.receipt().bridges.end(),[&](const auto&y){return y.native_id==x.native_id;});require(b!=d->native.receipt().bridges.end(),"source reference: native omitted actual surface bridge");
   const auto&adapter=ports[x.source_port].adapter;const auto start=adapter.field_point(NativePoint{x.native_origin}).value,end=adapter.field_point(NativePoint{x.native_position}).value;V displacement{};for(int k=0;k<3;++k)displacement[k]=end[k]-start[k];
   require(bits(b->field_start_m,start)&&bits(b->field_displacement_m,displacement)&&bits(b->charge_C,x.charge_C),"source reference: native bridge translation/charge mismatch");
  }
  require(bridge_count==d->native.receipt().bridges.size(),"source reference: unowned native bridge receipt");
  d->configuration=configuration.hex();d->source_before=bases.hex();d->source_trial=trials.hex();d->source_load_sha=loads.hex();d->identity=hash.hex();return PreparedReferenceTransfer(std::move(d));
 }
 PreparedReferenceTransfer prepare(const s::Owner&o,const s::Trial&t,const std::shared_ptr<const s::Field>&field,const s::FinalFieldView&actual,m::BirthScope&scope,m::FullTimeBirthCredit credit,std::span<const NeutralLink>neutral,std::span<const ElectronLink>electrons)const{
  const Port port{*this,o,t,field,actual,neutral,electrons};return prepare_group(std::span<const Port>(&port,1),scope,std::move(credit));
 }
};
}
