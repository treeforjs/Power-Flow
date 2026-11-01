#pragma once
// Native-v8 source stage. It owns no second thermal store and no chemistry.
#include "root_api.hpp"
#include "source_stock_reference.hpp"
#include "pf_arrival_v11/arrival_owner.hpp"
#include "pf_arrival_v11/native_v6_births.hpp"
#include "pic-electrode-heat-source-v7.4-mass-cache-20261101/binding/hall_source_binding_v7.hpp"
#include <map>
#include <memory>
#include <optional>
#include <set>

namespace powerflow::source_stage_v8 {
namespace hc=hall_corrective_v3;
namespace hb=hall_source_binding_v7;
namespace s7=electrode_heat_v7;
namespace sr=source_stock_reference_v2;
namespace pa=pf_arrival_v11;
namespace m=raikou::mfem_backend;
namespace av=m::arrival_v11;
namespace eh=electrode_heat_v2;
namespace se=surface_emission;
using hc::need;
class PreparedSourceStage;
struct SourceStageTestAccess;

// Only bind_port can create a port. Its source reads this exact retained Hall
// owner's thermal slice; equal-content foreign thermal owners cannot be used.
class Port {
  friend class PreparedSourceStage;
  friend Port bind_port(std::shared_ptr<hc::Owner>,hb::SliceSource,const s7::OwnerOptions&,
                        const m::GapPicOptions&,const m::FullTimeChemistryConfiguration&,
                        std::vector<sr::PatchRegistration>);
  std::shared_ptr<hc::Owner> hall_;
  hb::SliceSource slice_;
  std::string registry_sha256_,chemistry_configuration_sha256_,blob_sha256_;
  std::uint64_t run_seed_=0;
  std::shared_ptr<s7::Owner> source_;
  std::vector<sr::PatchRegistration> registration_;
  std::shared_ptr<const sr::Adapter> reference_;
  Port()=default;
 public:
  const auto& hall()const noexcept{return hall_;}
  const auto& registry_sha256()const noexcept{return registry_sha256_;}
  const auto& chemistry_configuration_sha256()const noexcept{return chemistry_configuration_sha256_;}
  const auto& blob_sha256()const noexcept{return blob_sha256_;}
  std::uint64_t run_seed()const noexcept{return run_seed_;}
  const auto& slice()const noexcept{return slice_;}
  const auto& source()const noexcept{return source_;}
  const auto& registration()const noexcept{return registration_;}
  const auto& reference()const noexcept{return reference_;}
};
inline Port bind_port(std::shared_ptr<hc::Owner> hall,hb::SliceSource slice,const s7::OwnerOptions& source_options,
                      const m::GapPicOptions& options,const m::FullTimeChemistryConfiguration& chemistry,
                      std::vector<sr::PatchRegistration> registrations) {
  need(bool(hall),"source port requires actual shared Hall owner");
  need(slice.body>=0&&slice.body<2&&slice.slice>=hall->mesh().kshort&&slice.slice<hall->mesh().kfeed,"source port outside actual Hall foil");
  std::set<std::size_t> patches;
  for(const auto&r:registrations){
    need(r.patch<slice.patches.size()&&patches.insert(r.patch).second,"source port patch registration duplicates or exceeds slice");
    need(r.upper_foil==(slice.body==1)&&sr::bits(r.z_center_m,hb::slice_center_z(hall->mesh(),slice.slice)),"source registration does not match actual Hall slice");
  }
  need(patches.size()==slice.patches.size(),"every configured source patch must be registered");
  Port p;p.hall_=std::move(hall);p.slice_=std::move(slice);p.registration_=std::move(registrations);
  p.source_=hb::bind(*p.hall_,p.slice_,source_options);
  need(&p.source_->thermal_state()==&p.hall_->thermal().slice(p.slice_.body,p.slice_.slice).state(),"source bound to foreign thermal state");
  const auto registry=m::full_time_chemistry_registry(chemistry,options);p.registry_sha256_=registry.registry_sha256;p.chemistry_configuration_sha256_=registry.configuration_sha256;p.blob_sha256_=registry.blob_sha256;p.run_seed_=chemistry.run_seed;
  p.reference_=std::make_shared<const sr::Adapter>(*p.source_,options,chemistry,p.registration_);
  return p;
}
struct OwnerSet {
  std::vector<Port> ports;
  std::shared_ptr<const pa::Policy> arrival_policy;
  std::shared_ptr<pa::Owner> arrival;
};
struct AcceptedSourceState {
  std::vector<std::string> source_states;
  std::string arrival_state,source_seal;
};
struct PublicationState {
  std::vector<std::string> source_states;
  std::string arrival_state,source_transition,receipt_sha256;
  double time_s=0;
  std::uint64_t thermal_revision=0;
};

class PreparedSourceStage {
  friend struct SourceStageTestAccess;
  struct Data {
    OwnerSet owners;
    m::FullTimeStepBinding common;
    std::string field_before,heat_before,arrival_before,source_transition;
    std::vector<std::string> source_before;
    std::shared_ptr<const s7::Field> field;
    s7::v3::TrialFrame frame;
    std::shared_ptr<const s7::MeshIdentity> mesh;
    foil::Vec E,B;
    std::vector<eh::Epoch> epochs;
    std::vector<s7::Trial> trials;
    std::vector<hc::SurfaceLoad> loads;
    std::optional<sr::PreparedReferenceTransfer> reference;
    std::optional<pa::Owner::PreparedCommit> arrival;
    double expected_electrons=0,expected_neutrals=0;
    std::optional<PublicationState> publication;
    std::string hall_after_field,hall_after_heat;
    bool consumed=false;
    s7::FinalFieldView view()const{return {frame,*mesh->shared(),E,B};}
  };
  std::unique_ptr<Data>d_;
  explicit PreparedSourceStage(std::unique_ptr<Data>d):d_(std::move(d)){}
  struct Patch {
    std::size_t port=0,local=0;
    const sr::PatchRegistration* registration=nullptr;
    std::optional<s7::PreparedPatch> exposure;
  };
  static se::Vector patch_origin(const eh::Coupler&heat,std::size_t local,double center) {
    const auto&g=heat.geometry();const auto&p=heat.patch(local);const double u=.5*(p.begin+p.end);
    se::Vector x{g.x0+.5*g.width,g.y0+.5*g.thickness,center};
    if(p.face==eh::Face::LowY||p.face==eh::Face::HighY){x[0]=g.x0+u*g.width;x[1]=p.face==eh::Face::LowY?g.y0:g.y0+g.thickness;}
    else {x[0]=p.face==eh::Face::LowX?g.x0:g.x0+g.width;x[1]=g.y0+u*g.thickness;}
    return x;
  }
  static s7::v3::PatchBinding patch_binding(const eh::Coupler&heat,std::size_t local,double center) {
    s7::v3::PatchBinding b;b.thermal_revision=heat.state().revision;b.source_accepted_state=b.thermal_revision+1;b.patch_id=local+1;
    b.thermal_mesh_identity=heat.geometry().revision;b.thermal_material_hash=eh::material_headers_sha256;
    b.material_phase_admitted=b.planar_barrier_admitted=b.electron_supply_admitted=true;
    b.temperature=s7::v3::HeldTemperature::AcceptedPatchMean;b.extrusion_center_z_m=center;b.extrusion_height_m=heat.geometry().height;return b;
  }
  static std::unique_ptr<Data> prepare_impl(const OwnerSet&owners,const AcceptedSourceState&accepted,const m::FullTimeStepBinding&common,
      m::BirthScope&scope,const s7::BirthPolicy&birth_policy,
      std::span<const s7::IntervalEvidence>intervals,const hc::JointProposal* field_proposal,
      std::span<const m::ArenaRateEntry> proposal_census);
 public:
  PreparedSourceStage(const PreparedSourceStage&)=delete;
  PreparedSourceStage&operator=(const PreparedSourceStage&)=delete;
  PreparedSourceStage(PreparedSourceStage&&)noexcept=default;
  PreparedSourceStage&operator=(PreparedSourceStage&&)noexcept=default;
  // Fixed patch-centre spatial quadrature, accepted patch-mean temperature and
  // accepted-start Hall fields. No random relocation after exposure sampling.
  static PreparedSourceStage prepare(const OwnerSet&owners,const AcceptedSourceState&accepted,const m::FullTimeStepBinding&common,
      m::BirthScope&scope,const s7::BirthPolicy&birth_policy,
      std::span<const s7::IntervalEvidence>intervals={}) {
    return PreparedSourceStage(prepare_impl(owners,accepted,common,scope,birth_policy,intervals,nullptr,{}));
  }
  // Retry using only an authentic, live proposal made by this exact Hall owner.
  // Source temperature/stock and arrival bases remain accepted; the caller must
  // abandon/recreate the native outer if its bounded transport-fork cap is hit.
  static PreparedSourceStage prepare_against_proposal(const OwnerSet&owners,const AcceptedSourceState&accepted,
      const m::FullTimeStepBinding&common,m::BirthScope&scope,const hc::JointProposal&proposal,
      std::span<const m::ArenaRateEntry> proposal_census,const s7::BirthPolicy&birth_policy,
      std::span<const s7::IntervalEvidence>intervals={}) {
    return PreparedSourceStage(prepare_impl(owners,accepted,common,scope,birth_policy,intervals,&proposal,proposal_census));
  }
  const m::PreparedBirths& native_births()const{return d_->reference->native_births();}
  const sr::PreparedReferenceTransfer& reference()const{return *d_->reference;}
  const std::vector<hc::SurfaceLoad>& loads()const noexcept{return d_->loads;}
  const std::vector<s7::Trial>& source_trials()const noexcept{return d_->trials;}
  const pa::Owner::PreparedCommit& prepared_arrival()const{return *d_->arrival;}
  const std::string& source_after_sha256()const noexcept{return d_->source_transition;}
  const std::shared_ptr<const s7::Field>& source_field()const noexcept{return d_->field;}
  const m::FullTimeStepBinding& common()const noexcept{return d_->common;}
  double expected_electrons()const noexcept{return d_->expected_electrons;}
  double expected_neutrals()const noexcept{return d_->expected_neutrals;}
  bool preflight_sources(const OwnerSet&)const noexcept;
  const PublicationState& prepare_publication(const OwnerSet&,const hc::JointProposal&);
  // Root must first preflight its native, gas and actual Hall capabilities.
  // Hall publication occurs before this call. No allocation occurs here.
  void publish_sources_prevalidated(const OwnerSet&)noexcept;
};

inline std::unique_ptr<PreparedSourceStage::Data> PreparedSourceStage::prepare_impl(
    const OwnerSet&owners,const AcceptedSourceState&base,const m::FullTimeStepBinding&common,m::BirthScope&scope,
    const s7::BirthPolicy&birth_policy,std::span<const s7::IntervalEvidence>intervals,
    const hc::JointProposal* field_proposal,std::span<const m::ArenaRateEntry> proposal_census) {
  need(!owners.ports.empty()&&owners.arrival&&owners.arrival_policy,"source stage needs ports and arrival owners");
  need(&owners.arrival->policy()==owners.arrival_policy.get(),"arrival owner does not retain the declared exact policy");
  need(owners.arrival_policy->configuration().run_seed==owners.ports.front().run_seed_,"source arrival run seed differs from actual native chemistry configuration");
  need(base.source_states.size()==owners.ports.size()&&base.arrival_state==owners.arrival->state_sha256()&&base.source_seal==scope.source_before_sha256(),"source/arrival accepted root differs from native source boundary");
  for(std::size_t i=0;i<owners.ports.size();++i)need(base.source_states[i]==owners.ports[i].source_->state_sha256(),"source accepted component state changed between steps");
  need(birth_policy.return_extent_m==0,"source stage requires explicit ordinary-native return collection");
  need(scope.accepted_step()==common.accepted_step&&owners.arrival->next_step()==common.accepted_step,"native/source/arrival common step mismatch");
  need(intervals.empty()||intervals.size()==owners.ports.size(),"interval evidence must cover every source port");
  auto out=std::make_unique<Data>();out->owners=owners;out->common=common;
  const auto&hall=*owners.ports.front().hall_;const auto&g=hall.mesh();
  need(common.begin_s==hall.time_s()&&common.accepted_step==hall.revision()&&common.end_s==common.begin_s+common.dt_s&&common.dt_s>0&&common.dt_s<=hall.model().dt_max(),"actual Hall/source common clock mismatch");
  need(common.accepted_field_sha256==hall.field_sha256()&&common.mesh_sha256==hall.motion_mesh_sha256()&&common.material_sha256==hall.model().material_sha256(),"source stage common binding differs from actual Hall");
  out->field_before=hall.field_sha256();out->heat_before=hall.heat_sha256();out->arrival_before=owners.arrival->state_sha256();
  std::string held_field_hash;
  if(field_proposal){
    need(hall.preflight(*field_proposal,proposal_census),"source retry requires an authentic live proposal of this exact Hall owner and census");
    const auto&binding=field_proposal->motion().binding();
    need(binding.accepted_step==common.accepted_step&&binding.begin_s==common.begin_s&&binding.end_s==common.end_s&&binding.dt_s==common.dt_s&&binding.common_clock_sha256==common.common_clock_sha256&&binding.mesh_sha256==common.mesh_sha256&&binding.material_sha256==common.material_sha256,"source retry Hall proposal differs from the actual common step");
    out->E=field_proposal->trial_E();out->B=field_proposal->trial_B();
    held_field_hash="source-stage-v8/proposal-midpoint/"+field_proposal->midpoint().total_field_sha256();
  } else {
    need(proposal_census.empty(),"source accepted-start preparation cannot carry a proposal census");
    const auto accepted=hall.model().at(hall.field());out->E=accepted.E();out->B=accepted.B();held_field_hash="source-stage-v8/accepted-start/"+accepted.total_field_sha256();
  }
  out->mesh=s7::MeshIdentity::of(g);auto&frame=out->frame;
  frame.role=s7::v3::FieldRole::TotalMidpoint;frame.accepted_revision=hall.revision();frame.trial_stage=1;frame.mesh_revision=1;frame.material_revision=2;
  frame.begin_s=common.begin_s;frame.end_s=common.end_s;frame.dt_s=common.dt_s;frame.field_stage_s=common.begin_s+.5*common.dt_s;
  frame.arena_hash=common.accepted_arena_sha256;frame.field_hash=std::move(held_field_hash);frame.mesh_hash=common.mesh_sha256;frame.material_hash=common.material_sha256;
  frame.species_hash=owners.ports.front().registry_sha256_;
  out->field=std::make_shared<const s7::Field>(out->mesh,frame,out->E,out->B,hall.model().dt_max());
  std::vector<Patch> patches;std::map<std::uint32_t,std::size_t> streams;std::set<std::pair<int,int>> slices;
  pa::StepInput input;input.step=common.accepted_step;
  const s7::v3::Exposure exposure{common.begin_s,common.dt_s,1,"source-stage-v8/fixed-patch-centre-quadrature"};
  for(std::size_t i=0;i<owners.ports.size();++i){
    const auto&p=owners.ports[i];need(p.hall_.get()==&hall&&p.source_&&p.reference_,"source ports belong to different Hall owners");
    const auto&first=owners.ports.front();need(p.registry_sha256_==first.registry_sha256_&&p.chemistry_configuration_sha256_==first.chemistry_configuration_sha256_&&p.blob_sha256_==first.blob_sha256_&&p.run_seed_==first.run_seed_,"source ports use different actual native catalogues/configurations");
    need(slices.insert({p.slice_.body,p.slice_.slice}).second,"duplicate Hall source slice");
    const auto&heat=hall.thermal().slice(p.slice_.body,p.slice_.slice);
    need(&p.source_->thermal_state()==&heat.state(),"source port thermal lifetime differs from retained Hall");
    need(p.source_->active_layers().empty(),"source stage has no retained microscopic-layer handoff implementation");
    out->source_before.push_back(p.source_->state_sha256());out->epochs.push_back(p.source_->make_epoch(common.dt_s,frame.field_hash));
    const auto quotes=p.source_->release_quotes(out->epochs.back());
    for(const auto&r:p.registration_){
      Patch work;work.port=i;work.local=r.patch;work.registration=&r;
      const auto patch=p.source_->source_patch(r.patch);const auto settings=p.source_->source_settings(r.patch);
      need(!patch.desorption_enabled||r.neutral_stream.has_value(),"enabled neutral desorption has no arrival stream");
      if(r.electron_stream){
        need(p.source_->enhancement_classes(r.patch)==1,"source stage requires a single emission enhancement class per patch; class partitioning is not implemented");
        const auto stream=*r.electron_stream;const auto&sc=owners.arrival_policy->stream(stream);
        need(sc.kind==pa::Kind::electron&&sc.species_state=="pfc-v3/e","electron arrival stream is not actual registry electron");
        need(streams.emplace(stream,patches.size()).second,"source stream mapped twice");
        const auto origin=patch_origin(heat,r.patch,r.z_center_m);const double F=se::extraction_field(out->field->surface(origin,heat.normal(r.patch)).E,heat.normal(r.patch));
        if(F>0&&F<1)throw std::runtime_error("source/emission_field_domain: positive extraction below1V/m has no qualified launch model");
        work.exposure.emplace(p.source_->prepare_patch(out->epochs.back(),r.patch,patch_binding(heat,r.patch,r.z_center_m),out->field,origin,exposure));
        const auto&pp=*work.exposure;input.electrons.push_back({stream,pp.expected_number(),pp.log_expected_number(),std::isfinite(pp.log_expected_number()),common.begin_s,common.dt_s,common.end_s});out->expected_electrons+=pp.expected_number();
      }
      if(r.neutral_stream){
        const auto stream=*r.neutral_stream;const auto&sc=owners.arrival_policy->stream(stream);
        need(sc.kind==pa::Kind::neutral&&sc.species_state==patch.state_id&&sr::bits(sc.target_weight,settings.quantized_neutral_weight),"neutral arrival/source state or exact weight differs");
        need(streams.emplace(stream,patches.size()).second,"source stream mapped twice");const auto&q=quotes.at(r.patch);
        input.neutrals.push_back({stream,q.released,q.log_released,p.source_->inventory(r.patch),common.begin_s,common.end_s});out->expected_neutrals+=q.released.rounded_for_rate();
      }
      patches.push_back(std::move(work));
    }
  }
  need(streams.size()==owners.arrival_policy->configuration().streams.size(),"arrival policy contains unmapped source streams");
  auto staged=owners.arrival->stage(owners.arrival->plan(input),scope.base_native_id());
  const auto n=owners.ports.size();std::vector<std::vector<s7::SurfaceBirth>>surface(n);
  std::vector<std::vector<std::vector<s7::NeutralLeaf>>>leaves(n);
  for(std::size_t i=0;i<n;++i)leaves[i].resize(owners.ports[i].source_->source_patches());
  for(const auto&b:staged.births()){
    const auto&patch=patches.at(streams.at(b.stream_id));const auto&p=owners.ports[patch.port];
    if(b.kind==pa::Kind::neutral){leaves[patch.port][patch.local].push_back({b.weight,eh::Draw{b.marks[0],b.marks[1],b.marks[2],b.marks[3],b.marks[4],b.marks[5]}});continue;}
    need(patch.exposure.has_value(),"electron birth has no genuine source exposure");
    eh::SourceEventIdentity event{"source-stage-v8/emission/"+b.packet_id(),"pfc-v3/e","electrode-electron-reservoir","source-stage-v8/surface-birth",b.packet_id(),b.ordinal+1,b.ordinal+1};
    surface[patch.port].push_back(p.source_->prepare_surface_birth(*patch.exposure,b.weight,b.marks[0],b.marks[1],b.marks[2],b.native_id,b.time_s,event,birth_policy));
  }
  out->trials.reserve(n);
  for(std::size_t i=0;i<n;++i){const auto&p=owners.ports[i];out->trials.push_back(p.source_->prepare_with_surface_births(out->epochs[i],out->field,out->view(),surface[i],leaves[i],{},intervals.empty()?nullptr:&intervals[i]));auto loads=hb::loads(out->trials.back(),p.slice_);out->loads.insert(out->loads.end(),loads.begin(),loads.end());}
  // Source owns neutral release times; bind them before native credit or the
  // prepared arrival ledger is made. Exact packet weights/order must agree.
  for(const auto&patch:patches)if(patch.registration->neutral_stream){
    const auto stream=*patch.registration->neutral_stream;const auto id=owners.ports[patch.port].source_->source_patch(patch.local).id;std::vector<double>times,weights,expected;
    for(const auto&p:out->trials[patch.port].packets())if(p.patch==id){times.push_back(p.time);weights.push_back(p.weight);}
    for(const auto&b:staged.births())if(b.stream_id==stream)expected.push_back(b.weight);
    need(weights==expected,"neutral source packet weights differ from arrival leaves");if(!times.empty())staged.bind_source_times(stream,times);
  }
  hc::Digest transition("source-stage-v8-transition-v1");transition.text(scope.source_before_sha256());transition.u64(common.accepted_step);transition.text(staged.plan().sha256());transition.text(staged.token());transition.u64(n);
  for(std::size_t i=0;i<n;++i){transition.text(out->source_before[i]);transition.text(out->trials[i].content_sha256());}
  out->source_transition=transition.hex();
  std::vector<av::BirthKinematics>kinematics;std::vector<std::vector<sr::NeutralLink>>neutral(n);std::vector<std::vector<sr::ElectronLink>>electrons(n);std::vector<const s7::SurfaceBirth*>electron_tokens;
  std::map<std::uint32_t,std::size_t>next_packet;
  for(const auto&b:staged.births()){
    const auto&patch=patches.at(streams.at(b.stream_id));const auto&p=owners.ports[patch.port];av::BirthKinematics k;
    if(b.kind==pa::Kind::electron){
      const auto&v=out->trials[patch.port].surface_births();const auto it=std::find_if(v.begin(),v.end(),[&](const auto&x){return x.canonical_id()==b.native_id;});need(it!=v.end(),"staged electron lost actual SurfaceBirth");
      k.position_m=p.reference_->native_surface_point(*p.source_,patch.local,sr::GlobalFieldPoint{it->origin_m()});for(int d=0;d<3;++d)k.position_m[d]+=it->phase().displacement_m[d];k.proper_velocity_m_s=it->phase().packet.particle.proper_velocity_m_s;k.origin="source-v7.3/surface/"+it->content_sha256();electrons[patch.port].push_back({b.native_id,b.packet_id()});electron_tokens.push_back(&*it);
    } else {
      const auto&packets=out->trials[patch.port].packets();const auto id=p.source_->source_patch(patch.local).id;const auto wanted=next_packet[b.stream_id]++;std::size_t seen=0,index=0;const eh::Packet*packet=nullptr;
      for(;index<packets.size();++index)if(packets[index].patch==id&&seen++==wanted){packet=&packets[index];break;}
      need(packet&&sr::bits(packet->weight,b.weight)&&sr::bits(packet->time,b.time_s),"neutral source packet differs from staged exact leaf/time");
      k.position_m=p.reference_->native_neutral_face_point(*p.source_,patch.local,sr::LocalThermalPoint{packet->x});k.proper_velocity_m_s=packet->proper_velocity_m_s;k.origin=packet->event_id;neutral[patch.port].push_back({index,b.packet_id()});electron_tokens.push_back(nullptr);
    }
    kinematics.push_back(std::move(k));
  }
  auto credit=av::credit(*owners.arrival_policy,staged,kinematics,av::plan_debits(staged.plan()),scope.source_before_sha256(),out->source_transition);
  for(std::size_t i=0;i<credit.births.size();++i)if(electron_tokens[i]){const auto&t=*electron_tokens[i];auto&b=credit.births[i];b.bridge=t.phase().packet.matching_distance_m>0;if(b.bridge){const auto&patch=patches.at(streams.at(b.stream_id));const auto&p=owners.ports[patch.port];b.bridge_start_m=p.reference_->native_surface_point(*p.source_,patch.local,sr::GlobalFieldPoint{t.origin_m()});}}
  std::vector<sr::Adapter::Port>ports;ports.reserve(n);
  for(std::size_t i=0;i<n;++i){const auto&p=owners.ports[i];ports.push_back({*p.reference_,*p.source_,out->trials[i],out->field,out->view(),neutral[i],electrons[i]});}
  out->reference.emplace(sr::Adapter::prepare_group(ports,scope,std::move(credit)));
  out->arrival.emplace(owners.arrival->prepare_commit(std::move(staged),av::receipt(out->reference->native_births().receipt())));
  return out;
}

inline bool PreparedSourceStage::preflight_sources(const OwnerSet&owners)const noexcept {
  try{
    if(!d_||d_->consumed||owners.ports.size()!=d_->owners.ports.size()||owners.arrival!=d_->owners.arrival||owners.arrival_policy!=d_->owners.arrival_policy||!owners.arrival->can_publish(*d_->arrival))return false;
    const auto&hall=*d_->owners.ports.front().hall_;
    if(hall.field_sha256()!=d_->field_before||hall.heat_sha256()!=d_->heat_before||owners.arrival->state_sha256()!=d_->arrival_before)return false;
    for(std::size_t i=0;i<owners.ports.size();++i){const auto&p=owners.ports[i];const auto&original=d_->owners.ports[i];if(p.hall_!=original.hall_||p.source_!=original.source_||p.reference_!=original.reference_||p.source_->state_sha256()!=d_->source_before[i]||!p.source_->preflight(d_->trials[i],d_->view()))return false;}
    return true;
  }catch(...){return false;}
}
inline const PublicationState& PreparedSourceStage::prepare_publication(const OwnerSet&owners,const hc::JointProposal&p) {
  need(d_&&!d_->publication,"source publication already prepared or moved");
  need(preflight_sources(owners),"source stage accepted owners changed before publication preparation");
  const auto&hall=*d_->owners.ports.front().hall_;
  need(p.end_s()==d_->common.end_s&&p.endpoint_revision()==hall.revision()+1&&p.configuration_sha256()==hall.configuration_sha256(),"source stage Hall proposal clock/configuration mismatch");
  const auto&expected=d_->reference->native_births().receipt().bridges;const auto&actual=p.motion().bridges();const auto&work=p.motion().bridge_work_J();
  need(expected.size()==actual.size()&&actual.size()==work.size(),"source stage Hall bridge coverage mismatch");
  std::map<std::uint64_t,const sr::SurfaceTransfer*> transfers;
  for(const auto&s:d_->reference->surface_transfers())if(s.bridge)need(transfers.emplace(s.native_id,&s).second,"duplicate source bridge native ID");
  need(transfers.size()==actual.size(),"source stage unmatched Hall bridge");
  for(std::size_t j=0;j<actual.size();++j){
    need(transfers.contains(actual[j].native_id),"Hall bridge has no source transfer");const auto&s=*transfers.at(actual[j].native_id);
    need(expected[j].native_id==actual[j].native_id&&sr::bits(expected[j].field_start_m,actual[j].field_start_m)&&sr::bits(expected[j].field_displacement_m,actual[j].field_displacement_m)&&sr::bits(expected[j].charge_C,actual[j].charge_C),"source stage Hall bridge operand mismatch");
    const double bound=1e-10*(std::abs(work[j])+std::abs(s.matching_work_J));
    need(std::isfinite(work[j])&&std::abs(work[j]-s.matching_work_J)<=bound,"source/Hall matching-field work mismatch; common field/source iteration required");
  }
  PublicationState state;state.source_transition=d_->source_transition;state.arrival_state=d_->arrival->candidate_state_sha256();state.time_s=p.end_s();state.thermal_revision=p.endpoint_revision();
  hc::Digest seal("source-stage-v8-future-state-v1");seal.text(state.source_transition);seal.text(state.arrival_state);seal.text(p.candidate_field_sha256());seal.text(p.candidate_heat_sha256());
  for(std::size_t i=0;i<owners.ports.size();++i){const auto&source=*owners.ports[i].source_;state.source_states.push_back(source.candidate_state_sha256(d_->trials[i],d_->view(),{source.bound_thermal_configuration_signature(),state.time_s,state.thermal_revision}));seal.text(state.source_states.back());}
  state.receipt_sha256=seal.hex();d_->hall_after_field=p.candidate_field_sha256();d_->hall_after_heat=p.candidate_heat_sha256();d_->publication.emplace(std::move(state));return *d_->publication;
}
inline void PreparedSourceStage::publish_sources_prevalidated(const OwnerSet&owners)noexcept {
  if(!d_||d_->consumed||owners.ports.size()!=d_->owners.ports.size()||owners.arrival!=d_->owners.arrival||!owners.arrival->can_publish(*d_->arrival))std::terminate();
  if(!d_->publication)std::terminate();const auto&hall=*d_->owners.ports.front().hall_;if(hall.field_sha256()!=d_->hall_after_field||hall.heat_sha256()!=d_->hall_after_heat||hall.time_s()!=d_->publication->time_s||hall.revision()!=d_->publication->thermal_revision)std::terminate();
  for(std::size_t i=0;i<owners.ports.size();++i){if(owners.ports[i].source_!=d_->owners.ports[i].source_||owners.ports[i].hall_!=d_->owners.ports[i].hall_)std::terminate();owners.ports[i].source_->publish_prevalidated(std::move(d_->trials[i]));}
  if(!owners.arrival->try_publish(*d_->arrival))std::terminate();d_->consumed=true;
}
} // namespace powerflow::source_stage_v8
