#pragma once
#include "surface_source_binding.hpp"
#include <bit>
#include <span>
#include <optional>
namespace electrode_heat_v3 {
namespace eh=electrode_heat_v2;namespace se=surface_emission;
using V3=se::Vector;
inline double dot(V3 a,V3 b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline double norm(V3 a){return std::sqrt(dot(a,a));}
inline bool finite(V3 a){return std::all_of(a.begin(),a.end(),[](double x){return std::isfinite(x);});}
enum class FieldRole {Unbound,TotalMidpoint};
struct TrialFrame {
 FieldRole role=FieldRole::Unbound;
 std::uint64_t accepted_revision=0,trial_stage=0,mesh_revision=0,material_revision=0;
 double begin_s=0,end_s=0,dt_s=0,field_stage_s=0;
 std::string arena_hash,field_hash,mesh_hash,material_hash,species_hash;
 void validate()const{eh::need(role==FieldRole::TotalMidpoint&&trial_stage>0&&mesh_revision>0&&material_revision>0&&std::isfinite(begin_s)&&std::isfinite(end_s)&&std::isfinite(dt_s)&&dt_s>0&&dt_s<=1e-14&&end_s==begin_s+dt_s&&end_s>begin_s&&field_stage_s==begin_s+.5*dt_s,"invalid immutable midpoint/common interval");eh::need(!arena_hash.empty()&&!field_hash.empty()&&!mesh_hash.empty()&&!material_hash.empty()&&!species_hash.empty(),"unsealed trial identities");}
 bool operator==(const TrialFrame&)const=default;
};
struct PrecisionPath {
 V3 origin_m{},from_relative_m{},to_relative_m{},displacement_m{};
 double start_s=0,end_s=0,charge_C=0,field_work_J=0;
 V3 absolute_from()const{V3 x=origin_m;for(int i=0;i<3;++i)x[i]+=from_relative_m[i];return x;}
 V3 absolute_to()const{V3 x=origin_m;for(int i=0;i<3;++i)x[i]+=to_relative_m[i];return x;}
};
inline PrecisionPath precision_path(const se::CurrentPath&p,V3 origin){PrecisionPath q;q.origin_m=origin;q.from_relative_m=p.from_m;q.to_relative_m=p.to_m;for(int i=0;i<3;++i)q.displacement_m[i]=p.to_m[i]-p.from_m[i];q.start_s=p.start_time_s;q.end_s=p.end_time_s;q.charge_C=p.charge_C;q.field_work_J=p.field_work_J;return q;}
struct FieldSample {V3 E{},B{};};
class TrialField {
 const TrialFrame frame_;
public:
 explicit TrialField(TrialFrame f):frame_(std::move(f)){frame_.validate();}virtual ~TrialField()=default;
 const TrialFrame& frame()const{return frame_;}
 virtual FieldSample sample(const PrecisionPath&)const=0;
 virtual FieldSample surface(V3 origin,V3 outward_normal)const=0;
 virtual bool admits_normal_quantum_bridge(const PrecisionPath&,V3 outward_normal)const=0;
};
enum class HeldTemperature {Unbound,AcceptedPatchMean,AcceptedPointTrace};
struct PatchBinding {
 std::uint64_t thermal_revision=0,source_accepted_state=0,patch_id=0;
 std::string thermal_mesh_identity,thermal_material_hash;double extrusion_center_z_m=std::numeric_limits<double>::quiet_NaN(),extrusion_height_m=std::numeric_limits<double>::quiet_NaN();
 bool material_phase_admitted=false,planar_barrier_admitted=false,electron_supply_admitted=false;
 HeldTemperature temperature=HeldTemperature::Unbound;
};
struct Exposure {double begin_s=0,duration_s=0,area_fraction=0;std::string quadrature_id;};
struct PreparedPatch {
 const eh::Coupler* heat=nullptr;eh::Epoch epoch;std::size_t patch_index=0;PatchBinding binding;
 std::shared_ptr<const TrialField> field;V3 origin_m{};se::PreparedFlux native;
 Exposure exposure;std::string exposure_key,thermal_signature;
 double log_expected_number=-INFINITY;bool has_log_expected_number=false,current_underflow=false;
};
inline std::string exposure_key(const TrialFrame&f,const PatchBinding&b,const Exposure&e,const V3&origin,double T,const se::Config&c){std::ostringstream s;auto str=[&](const std::string&x){s<<x.size()<<':'<<x;};s<<"source-exposure-v3/"<<f.accepted_revision<<'/'<<f.trial_stage<<'/'<<f.mesh_revision<<'/'<<f.material_revision<<'/'<<b.patch_id<<'/'<<b.thermal_revision<<'/'<<b.source_accepted_state<<' ';for(const auto&x:{f.field_hash,f.arena_hash,f.mesh_hash,f.material_hash,f.species_hash,b.thermal_mesh_identity,b.thermal_material_hash,e.quadrature_id})str(x);s<<std::hexfloat<<'/'<<f.begin_s<<'/'<<f.end_s<<'/'<<f.dt_s<<'/'<<f.field_stage_s<<'/'<<e.begin_s<<'/'<<e.duration_s<<'/'<<e.area_fraction<<'/'<<T<<'/'<<c.work_function_eV<<'/'<<c.energy_panels<<'/'<<c.action_nodes<<'/'<<c.tail_widths<<'/'<<int(b.temperature)<<'/'<<b.extrusion_center_z_m<<'/'<<b.extrusion_height_m;for(double x:origin)s<<'/'<<x;return s.str();}
inline PreparedPatch prepare_patch(const eh::Coupler&heat,const eh::Epoch&epoch,std::size_t patch,PatchBinding bind,std::shared_ptr<const TrialField>field,V3 origin,const Exposure&exposure,const se::Config&config){
 eh::need(bool(field),"missing immutable trial field");const auto&f=field->frame();f.validate();
 eh::need(epoch.t==heat.state().time&&epoch.revision==heat.state().revision&&epoch.dt==f.dt_s&&epoch.t==f.begin_s&&bind.thermal_revision==heat.state().revision&&bind.source_accepted_state>0&&bind.patch_id>0&&bind.thermal_mesh_identity==epoch.mesh_revision&&bind.thermal_material_hash==epoch.material_hash,"accepted thermal versus trial-field mapping mismatch");
 eh::need(bind.material_phase_admitted&&bind.planar_barrier_admitted&&bind.electron_supply_admitted&&(bind.temperature==HeldTemperature::AcceptedPatchMean||bind.temperature==HeldTemperature::AcceptedPointTrace),"explicit source/held-temperature regime admission required");eh::need(std::isfinite(bind.extrusion_center_z_m)&&std::isfinite(bind.extrusion_height_m)&&bind.extrusion_height_m==heat.geometry().height,"explicit physical axial slice mapping must use thermal height once");V3 local_origin=origin;local_origin[2]-=bind.extrusion_center_z_m;eh::admit_origin(heat,patch,local_origin);
 eh::need(std::isfinite(exposure.begin_s)&&std::isfinite(exposure.duration_s)&&exposure.duration_s>0&&exposure.duration_s<=f.dt_s&&exposure.begin_s>=f.begin_s&&exposure.begin_s+exposure.duration_s<=f.end_s&&exposure.area_fraction>0&&exposure.area_fraction<=1&&std::isfinite(exposure.area_fraction)&&!exposure.quadrature_id.empty(),"explicit local exposure/quadrature required");
 PreparedPatch p;p.heat=&heat;p.epoch=epoch;p.patch_index=patch;p.binding=std::move(bind);p.field=std::move(field);p.origin_m=origin;p.exposure=exposure;p.thermal_signature=heat.config_signature();auto n=heat.normal(patch);auto sample=p.field->surface(origin,n);eh::need(f.dt_s*se::qe/se::me*norm(sample.B)<=.05,"source full-common-interval gyro guard failed");eh::need(finite(sample.E)&&finite(sample.B),"nonfinite one-sided trial fields");
 se::Stage s;s.accepted_state=p.binding.source_accepted_state;s.field_accepted_state=s.accepted_state;s.thermal_accepted_state=s.accepted_state;s.mesh_revision=f.mesh_revision;s.material_revision=f.material_revision;s.patch_id=p.binding.patch_id;s.t=f.begin_s;s.dt=f.dt_s;s.stage_time=exposure.begin_s;s.field_time=f.field_stage_s;s.thermal_time=heat.state().time;s.total_E=sample.E;s.total_B_T=sample.B;s.outward_normal=n;s.area_m2=heat.area(patch)*exposure.area_fraction;s.T_K=p.binding.temperature==HeldTemperature::AcceptedPatchMean?heat.surface_temperatures().at(patch):heat.temperature_at(origin[0],origin[1]);s.material_temperature_phase_admitted=s.locally_planar_barrier_admitted=s.local_equilibrium_electron_supply_admitted=true;
 // Reuse the unchanged source-v1 supply/transmission law. Its legacy prepare()
 // conflates timestamps; v3 validates the distinct roles above, without relabeling.
 auto d=se::evaluate(s.T_K,se::extraction_field(s.total_E,s.outward_normal),config);
 double number=d.current_A_m2*s.area_m2*exposure.duration_s/se::qe,diagnostic=d.wall_heat_W_m2*s.area_m2*exposure.duration_s;eh::need(std::isfinite(number)&&std::isfinite(diagnostic),"exposure product overflow");p.current_underflow=d.current_underflow;p.has_log_expected_number=d.has_log_current;if(d.has_log_current)p.log_expected_number=d.log_current_A_m2+std::log(s.area_m2)+std::log(exposure.duration_s)-std::log(se::qe);p.native={s,std::move(d),number,diagnostic};p.exposure_key=exposure_key(f,p.binding,exposure,origin,s.T_K,config);for(const std::string&x:{p.thermal_signature,epoch.source_hash,std::string(eh::source_kernel_sha256)})p.exposure_key+=std::to_string(x.size())+":"+x;return p;
}
enum class BirthTime {Unbound,UniformHeldExposure,ExposureBeginLumped,SpecifiedPhysicalTime};
enum class WeightProvenance {Unbound,StagedCanonicalOwner,QualificationOnly};
struct BirthSample {
 double physical_weight=0,normal_uniform=0,transverse_uniform=0,azimuth_uniform=0,time_uniform=0,specified_birth_s=0;
 std::uint64_t event_id=0,canonical_particle_id=0,count_owner_revision=0;
 std::string count_owner_base_hash,count_owner_token,exposure_key;
 BirthTime time_policy=BirthTime::Unbound;WeightProvenance weight_provenance=WeightProvenance::Unbound;
};
struct Candidate {
 std::shared_ptr<const TrialField> field;const eh::Coupler* heat_owner=nullptr;
 std::uint64_t thermal_revision=0,canonical_particle_id=0;
 double thermal_time_s=0;std::string thermal_signature,exposure_key,count_owner_token;
 eh::BoundLayerState state;se::LayerExchange exchange;eh::Transfer heat;std::optional<BirthSample> birth_sample;
 std::vector<PrecisionPath> paths;bool qualification_only=false;int kernel_calls=0,field_iterations=0,refinements=0;
 double maximum_quantum_magnetic_fraction=0;
};
inline bool preflight(const Candidate&c,const eh::Coupler&h,const std::shared_ptr<const TrialField>&field,bool qualification=false)noexcept{return c.field==field&&c.heat_owner==&h&&c.thermal_revision==h.state().revision&&c.thermal_time_s==h.state().time&&c.thermal_signature==h.config_signature()&&(!c.qualification_only||qualification);}
inline void admit_policy(const se::LayerPolicy&p){eh::need(p.maximum_quantum_magnetic_energy_fraction>0&&p.maximum_quantum_magnetic_energy_fraction<=.2,"v3 forbids relaxed quantum magnetic applicability");}
inline bool path_work_matches(const TrialField&field,const PrecisionPath&p){auto a=field.sample(p);eh::need(finite(a.E)&&finite(a.B)&&finite(p.displacement_m)&&finite(p.origin_m),"nonfinite path gather");double W=p.charge_C*dot(a.E,p.displacement_m),scale=std::abs(p.field_work_J);for(int i=0;i<3;++i)scale+=std::abs(p.charge_C*a.E[i]*p.displacement_m[i]);return std::abs(W-p.field_work_J)<=1e-10*scale;}
inline void add(se::LayerExchange&a,const se::LayerExchange&b){a.physical_electron_number_change+=b.physical_electron_number_change;a.particle_charge_change_C+=b.particle_charge_change_C;a.electrode_charge_change_C+=b.electrode_charge_change_C;a.wall_heat_J+=b.wall_heat_J;a.particle_K_change_J+=b.particle_K_change_J;a.electron_reference_change_J+=b.electron_reference_change_J;a.image_energy_change_J+=b.image_energy_change_J;a.field_work_J+=b.field_work_J;}
inline Candidate candidate_base(const PreparedPatch&p,std::uint64_t particle_id){eh::need(p.heat&&p.field&&p.heat->state().revision==p.binding.thermal_revision&&p.heat->state().time==p.epoch.t&&p.heat->config_signature()==p.thermal_signature,"prepared patch no longer at accepted thermal base");eh::need(particle_id>0,"positive canonical source particle ID required");Candidate c;c.field=p.field;c.heat_owner=p.heat;c.thermal_revision=p.binding.thermal_revision;c.thermal_time_s=p.epoch.t;c.thermal_signature=p.thermal_signature;c.canonical_particle_id=particle_id;c.exposure_key=p.exposure_key;return c;}
inline void finalize(Candidate&c,const PreparedPatch&p,const eh::SourceEventIdentity&id,double time){V3 wall=c.state.surface_origin_m;if(c.state.native.status==se::LayerStatus::Recollected)for(int i=0;i<3;++i)wall[i]+=c.state.native.position_from_surface_m[i];c.heat=eh::source_heat(c.exchange,id,wall,time);double scale=0;for(double x:{c.exchange.wall_heat_J,c.exchange.particle_K_change_J,c.exchange.electron_reference_change_J,c.exchange.image_energy_change_J,c.exchange.field_work_J})scale+=std::abs(x);eh::need(std::abs(c.exchange.residual_J())<=1e-10*scale,"source incremental aggregate energy failed");for(const auto&path:c.paths)eh::need(path.start_s>=p.field->frame().begin_s&&path.end_s<=p.field->frame().end_s&&path.end_s>=path.start_s&&path_work_matches(*p.field,path),"source path/common-stage work unresolved");}
inline Candidate begin_source(const PreparedPatch&p,const se::LayerPolicy&policy,const BirthSample&sample,const eh::SourceEventIdentity&id){
 admit_policy(policy);const auto&frame=p.field->frame();eh::need((sample.weight_provenance==WeightProvenance::StagedCanonicalOwner||sample.weight_provenance==WeightProvenance::QualificationOnly)&&sample.event_id>0&&sample.canonical_particle_id>0&&std::isfinite(sample.physical_weight)&&sample.physical_weight>0&&sample.exposure_key==p.exposure_key&&!sample.count_owner_token.empty()&&!sample.count_owner_base_hash.empty()&&sample.count_owner_revision==frame.accepted_revision,"explicit staged count-owner receipt/positive representable packet weight required");double birth=0;if(sample.time_policy==BirthTime::UniformHeldExposure){eh::need(sample.time_uniform>0&&sample.time_uniform<1,"caller-owned open birth-time draw required");birth=p.exposure.begin_s+sample.time_uniform*p.exposure.duration_s;}else if(sample.time_policy==BirthTime::ExposureBeginLumped)birth=p.exposure.begin_s;else if(sample.time_policy==BirthTime::SpecifiedPhysicalTime)birth=sample.specified_birth_s;else eh::need(false,"explicit microscopic birth-time policy required");eh::need(std::isfinite(birth)&&birth>=p.exposure.begin_s&&birth<=p.exposure.begin_s+p.exposure.duration_s,"birth outside sealed exposure");auto flux=p.native;flux.stage.stage_time=birth;auto trial=se::begin_layer(flux,policy,sample.physical_weight,sample.normal_uniform,sample.transverse_uniform,sample.azimuth_uniform,sample.event_id);Candidate c=candidate_base(p,sample.canonical_particle_id);c.kernel_calls=1;c.qualification_only=sample.weight_provenance==WeightProvenance::QualificationOnly;c.count_owner_token=sample.count_owner_token;c.birth_sample=sample;c.state={trial.state,p.origin_m,p.patch_index,p.binding.thermal_mesh_identity,p.binding.thermal_material_hash};c.exchange=trial.exchange;c.maximum_quantum_magnetic_fraction=trial.state.quantum_magnetic_energy_fraction;for(const auto&path:trial.paths)c.paths.push_back(precision_path(path,p.origin_m));
 // Native birth uses B only for its quantum applicability guard. Re-run that
 // SAME guard at the other end of the normal bridge, without another birth debit.
 if(!c.paths.empty()){eh::need(p.field->admits_normal_quantum_bridge(c.paths.front(),p.native.stage.outward_normal),"source-v1 constant-normal-field barrier not admitted across trial-field cells");auto endpoint=c.paths.back();endpoint.from_relative_m=endpoint.to_relative_m;endpoint.displacement_m={};auto far=p.field->sample(endpoint);auto check_flux=flux;check_flux.stage.total_B_T=far.B;auto guard=se::begin_layer(check_flux,policy,sample.physical_weight,sample.normal_uniform,sample.transverse_uniform,sample.azimuth_uniform,sample.event_id);++c.kernel_calls;c.maximum_quantum_magnetic_fraction=std::max(c.maximum_quantum_magnetic_fraction,guard.state.quantum_magnetic_energy_fraction);}
 finalize(c,p,id,birth);return c;
}
struct AdvanceControls {int maximum_field_iterations=24,maximum_refinements=32,maximum_kernel_calls=4096;};
inline bool field_close(V3 a,V3 b){V3 d{};for(int i=0;i<3;++i)d[i]=a[i]-b[i];return norm(d)<=2e-13*(norm(a)+norm(b));}
inline Candidate advance_source(const PreparedPatch&p,const eh::BoundLayerState&accepted_layer,std::uint64_t canonical_id,const se::LayerPolicy&policy,double duration_s,const eh::SourceEventIdentity&id,const AdvanceControls&controls={}){
 admit_policy(policy);eh::need(controls.maximum_field_iterations>0&&controls.maximum_field_iterations<=24&&controls.maximum_refinements>=0&&controls.maximum_refinements<=32&&controls.maximum_kernel_calls>0&&controls.maximum_kernel_calls<=4096,"bounded path-field iteration controls required");eh::need(accepted_layer.kernel_hash==eh::source_kernel_sha256&&accepted_layer.patch_index==p.patch_index&&accepted_layer.mesh_identity==p.binding.thermal_mesh_identity&&accepted_layer.material_hash==p.binding.thermal_material_hash&&accepted_layer.surface_origin_m==p.origin_m&&accepted_layer.native.status==se::LayerStatus::Active,"accepted active source geometry/material/origin mismatch");Candidate out=candidate_base(p,canonical_id);out.state=accepted_layer;double target=accepted_layer.native.current_time_s+duration_s;const auto&frame=p.field->frame();eh::need(std::isfinite(duration_s)&&duration_s>0&&duration_s<=frame.dt_s&&accepted_layer.native.current_time_s>=frame.begin_s&&target<=frame.end_s&&target>accepted_layer.native.current_time_s,"active source duration outside common interval");
 while(out.state.native.status==se::LayerStatus::Active&&out.state.native.current_time_s<target){const auto base=out.state.native;double h=target-base.current_time_s;bool admitted=false;for(int refine=0;refine<=controls.maximum_refinements&&!admitted;++refine){PrecisionPath at;at.origin_m=p.origin_m;at.from_relative_m=at.to_relative_m=base.position_from_surface_m;auto guess=p.field->sample(at);bool shrink=false;for(int iteration=0;iteration<controls.maximum_field_iterations;++iteration){eh::need(out.kernel_calls<controls.maximum_kernel_calls,"source path-field kernel-call budget exceeded; reject/refine common trial");auto stage=p.native.stage;stage.stage_time=stage.field_time=frame.field_stage_s;stage.total_E=guess.E;stage.total_B_T=guess.B;auto trial=se::advance_layer(base,stage,policy,h);++out.kernel_calls;++out.field_iterations;const se::CurrentPath* moving=nullptr;int moves=0;for(const auto&path:trial.paths)if(path.end_time_s>path.start_time_s){moving=&path;++moves;}if(moves!=1){shrink=true;break;}auto actual=p.field->sample(precision_path(*moving,p.origin_m));eh::need(frame.dt_s*se::qe/se::me*norm(actual.B)<=.05,"active source full-common-interval gyro guard failed");if(!field_close(guess.E,actual.E)||!field_close(guess.B,actual.B)){guess=actual;continue;}std::vector<PrecisionPath> paths;bool work_ok=true;for(const auto&path:trial.paths){paths.push_back(precision_path(path,p.origin_m));work_ok=work_ok&&path_work_matches(*p.field,paths.back());}if(!work_ok){shrink=true;break;}out.state.native=trial.state;add(out.exchange,trial.exchange);out.paths.insert(out.paths.end(),paths.begin(),paths.end());admitted=true;break;}
 if(!admitted){++out.refinements;h*=.5;eh::need(h>0&&base.current_time_s+h>base.current_time_s,"unresolved spatial source path at native clock resolution; no work correction");if(refine==controls.maximum_refinements)eh::need(false,"unresolved spatial bridge/path work or field iteration; reject common trial");}(void)shrink;}
 eh::need(admitted,"source path-field trial not admitted");}
 finalize(out,p,id,out.state.native.current_time_s);return out;
}
inline eh::Trial prepare_heat(eh::Coupler&heat,const eh::Epoch&epoch,const std::shared_ptr<const TrialField>&field,mfem::Coefficient&joule,std::span<const Candidate>sources,const std::vector<eh::Transfer>&other,const std::vector<std::vector<eh::Draw>>&draws,bool qualification=false){auto transfers=other;for(const auto&c:sources){eh::need(preflight(c,heat,field,qualification),"source receipt does not bind this accepted thermal state and nonlinear trial");transfers.push_back(c.heat);}return heat.prepare(epoch,joule,transfers,draws);}
}
