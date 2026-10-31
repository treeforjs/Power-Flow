#pragma once
// Heat/source owner v7 (namespace electrode_heat_v7): electrode surface physics on top of owner v6.1.
// A versioned copy of owner v6.1 (pic-electrode-heat-source-v6.1-20261031, sealed, text below) with
// the v7 additions of API-v7.txt:
//  * configurable trial interval bound (OwnerOptions::interval, default 1e-14 s = v6.1) through the
//    versioned trial layer electrode_heat_v7k (trial_binding_v7.hpp), and the source-owner step bound;
//  * surface state per patch (SourceSettings::surface): work function of oxidized SS304 with the
//    adsorbed-water coverage of the exact stock; field-enhancement classes (prepare_class), phase map;
//  * emission births at the macroscopic matching plane (sample_birth: matching displacement, image tail
//    restored, sub-cell return in retarding fields);
//  * wall responses to first hits (respond): secondary/backscattered electrons, ion-induced emission,
//    ion and neutral reflection, sticking into the stock, electron/ion-stimulated desorption through
//    the quantized release, exact ledgers.
// With every v7 option at its benchmark value the owner is the v6.1 owner: same digests, values,
// decisions and checkpoint bytes (tests/test_source_v7.cpp, the sealed v6 qualification via compat-v6).
// Original v6.1 header text:
// Heat/source owner v6.1: a performance revision of owner v6 (pic-electrode-heat-source-
// v6-20261030, sealed; API-v6.txt). Every configuration, state, prepared-patch, candidate,
// trial, field-content and mesh digest is computed from the same bytes as in v6 and has the
// same value; every accepted/refused outcome is the same. Changes (API-v6.1.txt):
//  * Field: the mesh identity is an immutable snapshot hashed once (MeshIdentity); the field
//    content digest is computed once at construction; matches() compares the view's frame,
//    mesh (identity of the immutable snapshot, else exact member-wise equality) and E/B
//    arrays (bit for bit) instead of re-hashing the whole mesh and both arrays per call.
//  * prepare_patch: the sealed v3 prepare_patch steps with the emission spectrum from
//    emission::evaluate (bit-identical to surface_emission::evaluate, fast_emission.hpp);
//    PreparedPatch::flux() exposes the prepared flux (for sampling without re-preparing).
//  * state_sha256(): cached per (owner generation, thermal revision, thermal time).
//  * detail::Encoding appends 8-byte words in bulk (same bytes).
// Namespace electrode_heat_v6_1, so v6 and v6.1 can be compared in one program.
// Include the sealed Hall axial_surface_path_port.hpp first, as required by v3.
#include "trial_binding_v7.hpp"
#include "admission_v7.hpp"
#include "mesh_identity_v7.hpp"
#include "fast_emission.hpp"
#include "surface_data_v7.hpp"
#include "emission_hot_v7.hpp"
#include <numeric>
#include HEAT_SOURCE_EXACT_COUNT_HEADER
#include <map>
#include <optional>
#include <span>
namespace electrode_heat_v7 {
namespace v3=electrode_heat_v7k;namespace eh=electrode_heat_v2;namespace se=surface_emission;namespace emission=electrode_heat_v6_1::emission;
using Count=raikou::mfem_backend::ExactParticleCount;
using V3=se::Vector;using Particle=raikou::mfem_backend::GapParticle;
class Owner;
// v7 emission dispatch: the v6.1 exact fast path up to 3000 K (bit-identical to the sealed kernel), the extended
// evaluate above (only when the owner admits T_max > 3000 K).
inline se::Distribution evaluate_emission(double T,double F,const se::Config&c,double T_max){
 if(T>3000&&T_max>3000)return emission_hot::evaluate(T,F,c,T_max);
 return emission::evaluate(T,F,c);}
struct PatchIdentity {std::string physical_surface_id,material_layer_id;double physical_area_m2=0;bool operator==(const PatchIdentity&)const=default;};
namespace detail {
inline void policy(Encoding&e,const se::LayerPolicy&p){for(double x:{p.local_linear_vacuum_extent_m,p.image_force_fraction,p.image_energy_fraction,p.maximum_quantum_magnetic_energy_fraction,p.maximum_gyroangle,p.maximum_relative_normal_step})e.d(x);e.u(p.maximum_substeps);e.u(p.PIC_already_applies_microscopic_single_electron_image_force);e.u(p.own_return_unit_accommodation);}
inline std::string policy_hash(const se::LayerPolicy&p){Encoding e;policy(e,p);return e.sha();}
inline void exchange(Encoding&e,const se::LayerExchange&x){for(double a:{x.physical_electron_number_change,x.particle_charge_change_C,x.electrode_charge_change_C,x.wall_heat_J,x.particle_K_change_J,x.electron_reference_change_J,x.image_energy_change_J,x.field_work_J})e.d(a);}

// Weighted values of one kernel evaluation, in a fixed order with fixed names.
inline std::array<double,8> operands(const se::LayerExchange&x){return {x.physical_electron_number_change,x.particle_charge_change_C,x.electrode_charge_change_C,x.wall_heat_J,x.particle_K_change_J,x.electron_reference_change_J,x.image_energy_change_J,x.field_work_J};}
inline constexpr std::array<const char*,8> operand_names{"particle count","particle charge C","electrode charge C","wall heat J","particle kinetic increment J","electron reference J","image energy J","field work J"};
inline std::array<double,9> ledger_values(const se::Ledger&l){return {l.physical_electron_number,l.particle_charge_C,l.electrode_charge_C,l.wall_heat_J,l.outgoing_K_J,l.incoming_K_J,l.electron_reference_change_J,l.image_reference_change_J,l.resolved_matching_field_work_J};}
inline constexpr std::array<const char*,9> ledger_names{"ledger particle count","ledger particle charge C","ledger electrode charge C","ledger wall heat J","ledger outgoing kinetic J","ledger incoming kinetic J","ledger electron reference J","ledger image reference J","ledger matching field work J"};
inline std::array<double,5> heat_values(const eh::Transfer&t){return {t.source_wall_J,t.source_particle_K_change_J,t.source_electron_reference_change_J,t.source_image_energy_change_J,t.source_field_work_J};}
inline constexpr std::array<const char*,5> heat_names{"heat receipt source wall J","heat receipt particle kinetic increment J","heat receipt electron reference J","heat receipt image energy J","heat receipt field work J"};
inline double energy_scale(const se::LayerExchange&x){double s=std::abs(x.wall_heat_J)+std::abs(x.particle_K_change_J)+std::abs(x.electron_reference_change_J)+std::abs(x.image_energy_change_J)+std::abs(x.field_work_J);return std::isfinite(s)&&s>0?s:se::qe;}

inline void source_exchange_represented(const se::LayerExchange&unit,const se::LayerExchange&actual,double w,const std::string&id){
 auto a=operands(unit),b=operands(actual);for(int i=0;i<8;++i)require_product(a[i],w,b[i],operand_names[i],id);}
inline void admit_weighted_exchange(const se::LayerExchange&u,double w,const std::string&id){se::LayerExchange predicted{u.physical_electron_number_change*w,u.particle_charge_change_C*w,u.electrode_charge_change_C*w,u.wall_heat_J*w,u.particle_K_change_J*w,u.electron_reference_change_J*w,u.image_energy_change_J*w,u.field_work_J*w};source_exchange_represented(u,predicted,w,id);}
template<std::size_t N>void predicted_values(const Weighted&W,const std::array<const char*,N>&names,const std::array<double,N>&unit,const std::array<double,N>&tilde){for(std::size_t i=0;i<N;++i)W.predicted(names[i],unit[i],tilde[i]);}
template<std::size_t N>void consistent_values(const Weighted&W,const std::array<const char*,N>&names,const std::array<double,N>&unit,const std::array<double,N>&tilde,const std::array<double,N>&actual){for(std::size_t i=0;i<N;++i)W.consistent(names[i],unit[i],tilde[i],actual[i]);}

inline void configuration(Encoding&e,const eh::Geometry&g,const eh::Material&m,const std::vector<eh::Patch>&patches,double initial_T,double start,const std::string&source){e.u(g.nx);e.u(g.ny);for(double x:{g.x0,g.y0,g.width,g.thickness,g.height,initial_T,start})e.d(x);e.s(g.revision);e.s(g.electrode);for(const auto*p:{&m.id,&m.phase,&m.data_hash,&m.volume_convention,&m.native_callback_contract})e.s(*p);e.s(source);e.u(m.knots.size());for(const auto&k:m.knots)for(double x:{k.T,k.h,k.k,k.sigma})e.d(x);e.u(patches.size());for(const auto&p:patches){for(const auto*q:{&p.id,&p.state_id,&p.parent_id,&p.source_hash})e.s(*q);e.u(p.desorption_enabled);e.u(int(p.face));for(double x:{p.begin,p.end,p.initial_particles,p.initial_low_particles,p.mass_kg,p.attempt_s,p.activation_eV,p.desorption_eV,p.internal_eV,p.monolayer_m2,p.rate_alpha,p.binding_alpha,p.source_T_min,p.source_T_max,p.prefactor_reference_T})e.d(x);e.u(p.coverage_temkin);e.u(int(p.prefactor));}}
inline void budget(Encoding&e,const eh::Budget&b){for(double x:{b.source_wall_J,b.source_particle_K_change_J,b.source_electron_reference_change_J,b.source_image_energy_change_J,b.source_field_work_J,b.neutral_relativistic_conversion_J,b.work_function_reference_part_J,b.joule_J,b.incoming_K_J,b.outgoing_K_J,b.reference_out_J,b.nottingham_wall_J,b.solid_collision_J,b.authored_heat_J,b.neutral_K_J,b.neutral_internal_J,b.adsorbed_reference_change_J,b.expected_neutral_K_J,b.thermal_delta_J,b.wall_net_J,b.residual_J,b.roundoff_bound_J})e.d(x);e.u(b.newton_iterations);}
inline void count(Encoding&e,const Count&n){e.doubles(n.positive_carriers());}
inline void layer(Encoding&e,const eh::BoundLayerState&s){std::ostringstream f(std::ios::binary);eh::write_bound_layer_state(f,s);e.s(f.str());}
inline void identity(Encoding&e,const PatchIdentity&p){e.s(p.physical_surface_id);e.s(p.material_layer_id);e.d(p.physical_area_m2);}
inline void transfer(Encoding&e,const eh::Transfer&t){for(const auto*s:{&t.producer,&t.state_id,&t.parent_id,&t.channel_id,&t.energy_owner})e.s(*s);e.u(t.sequence);e.u(t.owner_sequence);for(double x:{t.event_time_s,t.x,t.y,t.incoming_K_J,t.outgoing_K_J,t.reference_out_J,t.solid_collision_J,t.authored_heat_J,t.nottingham_wall_J,t.work_function_reference_part_J,t.source_wall_J,t.source_particle_K_change_J,t.source_electron_reference_change_J,t.source_image_energy_change_J,t.source_field_work_J})e.d(x);e.u(t.surface);e.u(t.complete_nottingham);e.u(t.certified_source_exchange);}
// Exact comparison stock >= amount: the first positive carrier is the stock
// truncated to 53 significant bits; a binary64 amount above it exceeds the stock.
inline bool covers(const Count&stock,double amount){if(stock.zero())return amount==0;return stock.positive_carriers().front()>=amount;}
// The sealed electrode_heat_v3::prepare_patch (trial_source_binding.hpp), verbatim except that
// the spectrum is emission::evaluate (bit-identical to surface_emission::evaluate).
// v7: beta (enhancement class) multiplies the one-sided extraction field seen by the kernel (beta == 1: unchanged).
// v7: T_override (phase-owner patch temperature; NaN: the bound thermal view), T_max (extended emission domain),
// lambda (Richardson factor on the current; 1: unchanged).
inline v3::PreparedPatch prepare_patch(const eh::Coupler&heat,const eh::Epoch&epoch,std::size_t patch,v3::PatchBinding bind,std::shared_ptr<const v3::TrialField>field,V3 origin,const v3::Exposure&exposure,const se::Config&config,double beta=1,double T_override=std::numeric_limits<double>::quiet_NaN(),double T_max=3000,double lambda=1){
 eh::need(bool(field),"missing immutable trial field");const auto&f=field->frame();v3::validate(f,field->interval_cap_s());
 eh::need(epoch.t==heat.state().time&&epoch.revision==heat.state().revision&&epoch.dt==f.dt_s&&epoch.t==f.begin_s&&bind.thermal_revision==heat.state().revision&&bind.source_accepted_state>0&&bind.patch_id>0&&bind.thermal_mesh_identity==epoch.mesh_revision&&bind.thermal_material_hash==epoch.material_hash,"accepted thermal versus trial-field mapping mismatch");
 eh::need(bind.material_phase_admitted&&bind.planar_barrier_admitted&&bind.electron_supply_admitted&&(bind.temperature==v3::HeldTemperature::AcceptedPatchMean||bind.temperature==v3::HeldTemperature::AcceptedPointTrace),"explicit source/held-temperature regime admission required");eh::need(std::isfinite(bind.extrusion_center_z_m)&&std::isfinite(bind.extrusion_height_m)&&bind.extrusion_height_m==heat.geometry().height,"explicit physical axial slice mapping must use thermal height once");V3 local_origin=origin;local_origin[2]-=bind.extrusion_center_z_m;eh::admit_origin(heat,patch,local_origin);
 eh::need(std::isfinite(exposure.begin_s)&&std::isfinite(exposure.duration_s)&&exposure.duration_s>0&&exposure.duration_s<=f.dt_s&&exposure.begin_s>=f.begin_s&&exposure.begin_s+exposure.duration_s<=f.end_s&&exposure.area_fraction>0&&exposure.area_fraction<=1&&std::isfinite(exposure.area_fraction)&&!exposure.quadrature_id.empty(),"explicit local exposure/quadrature required");
 v3::PreparedPatch p;p.heat=&heat;p.epoch=epoch;p.patch_index=patch;p.binding=std::move(bind);p.field=std::move(field);p.origin_m=origin;p.exposure=exposure;p.thermal_signature=heat.config_signature();auto n=heat.normal(patch);auto sample=p.field->surface(origin,n);eh::need(f.dt_s*se::qe/se::me*v3::norm(sample.B)<=.05,"source full-common-interval gyro guard failed");eh::need(v3::finite(sample.E)&&v3::finite(sample.B),"nonfinite one-sided trial fields");
 se::Stage s;s.accepted_state=p.binding.source_accepted_state;s.field_accepted_state=s.accepted_state;s.thermal_accepted_state=s.accepted_state;s.mesh_revision=f.mesh_revision;s.material_revision=f.material_revision;s.patch_id=p.binding.patch_id;s.t=f.begin_s;s.dt=f.dt_s;s.stage_time=exposure.begin_s;s.field_time=f.field_stage_s;s.thermal_time=heat.state().time;s.total_E=sample.E;s.total_B_T=sample.B;s.outward_normal=n;s.area_m2=heat.area(patch)*exposure.area_fraction;s.T_K=p.binding.temperature==v3::HeldTemperature::AcceptedPatchMean?heat.surface_temperatures().at(patch):heat.temperature_at(origin[0],origin[1]);if(!std::isnan(T_override))s.T_K=T_override;s.material_temperature_phase_admitted=s.locally_planar_barrier_admitted=s.local_equilibrium_electron_supply_admitted=true;
 double F=se::extraction_field(s.total_E,s.outward_normal);if(beta!=1)F*=beta;
 auto d=evaluate_emission(s.T_K,F,config,T_max);
 if(lambda!=1&&d.has_log_current){d.current_A_m2*=lambda;d.log_current_A_m2+=std::log(lambda);d.wall_heat_W_m2*=lambda;d.current_underflow=d.current_A_m2==0;}
 double number=d.current_A_m2*s.area_m2*exposure.duration_s/se::qe,diagnostic=d.wall_heat_W_m2*s.area_m2*exposure.duration_s;eh::need(std::isfinite(number)&&std::isfinite(diagnostic),"exposure product overflow");p.current_underflow=d.current_underflow;p.has_log_expected_number=d.has_log_current;if(d.has_log_current)p.log_expected_number=d.log_current_A_m2+std::log(s.area_m2)+std::log(exposure.duration_s)-std::log(se::qe);p.native={s,std::move(d),number,diagnostic};p.exposure_key=v3::exposure_key(f,p.binding,exposure,origin,s.T_K,config);for(const std::string&x:{p.thermal_signature,epoch.source_hash,std::string(eh::source_kernel_sha256)})p.exposure_key+=std::to_string(x.size())+":"+x;return p;
}
}
struct FinalFieldView {const v3::TrialFrame&frame;const axial::Mesh&mesh;const foil::Vec&E;const foil::Vec&B;};
class Field final {
 friend class Owner;
 const std::shared_ptr<const MeshIdentity>identity_;
 const std::shared_ptr<const axial::Mesh>mesh_;
 const std::shared_ptr<const v3::AxialTrialField>native_;
 const std::string mesh_sha_,content_sha_;
 static std::string content(const v3::TrialFrame&f,const std::string&m,const foil::Vec&E,const foil::Vec&B,double cap=v3::v61_trial_interval_s){detail::Encoding e;e.bytes.reserve(16*(E.size()+B.size())+1024);e.s("source-v6-total-midpoint");e.s(m);e.u(f.accepted_revision);e.u(f.trial_stage);e.u(f.mesh_revision);e.u(f.material_revision);for(double x:{f.begin_s,f.end_s,f.dt_s,f.field_stage_s})e.d(x);for(const auto*s:{&f.arena_hash,&f.field_hash,&f.mesh_hash,&f.material_hash,&f.species_hash})e.s(*s);e.doubles(E);e.doubles(B);if(cap!=v3::v61_trial_interval_s){e.s("v7-trial-interval-cap");e.d(cap);}return e.sha();}
 static std::shared_ptr<const MeshIdentity>checked(std::shared_ptr<const MeshIdentity>id){eh::need(bool(id),"v6.1 immutable mesh identity required");return id;}
public:
 // v6 constructor: the mesh identity of g (an immutable snapshot, reused while g is exactly equal to it).
 // v7: interval_cap_s is the trial interval bound (API-v7.txt 3); 1e-14 is the v6.1 value and digest.
 Field(const axial::Mesh&g,v3::TrialFrame f,foil::Vec E,foil::Vec B,double interval_cap_s=v3::v61_trial_interval_s):Field(MeshIdentity::of(g),std::move(f),std::move(E),std::move(B),interval_cap_s){}
 // v6.1: an explicit immutable mesh identity (no copy, no hash).
 Field(std::shared_ptr<const MeshIdentity>id,v3::TrialFrame f,foil::Vec E,foil::Vec B,double interval_cap_s=v3::v61_trial_interval_s):identity_(checked(std::move(id))),mesh_(identity_->shared()),native_(std::make_shared<const v3::AxialTrialField>(mesh_,std::move(f),std::move(E),std::move(B),interval_cap_s)),mesh_sha_(identity_->sha256()),content_sha_(content(native_->frame(),mesh_sha_,native_->total_mid_E(),native_->total_mid_B(),interval_cap_s)){}
 double interval_cap_s()const{return native_->interval_cap_s();}
 const v3::TrialFrame&frame()const{return native_->frame();}
 const std::string&content_sha256()const{return content_sha_;}const std::string&mesh_sha256()const{return mesh_sha_;}
 const std::shared_ptr<const MeshIdentity>&mesh_identity()const noexcept{return identity_;}
 // Same decision as v6 (frame equal, mesh digest equal, content digest equal), by exact comparison.
 bool matches(const FinalFieldView&v)const{return frame()==v.frame&&identity_->same_as(v.mesh)&&detail::same_doubles(v.E,native_->total_mid_E())&&detail::same_doubles(v.B,native_->total_mid_B());}
 v3::FieldSample surface(V3 x,V3 n)const{return native_->surface(x,n);}
};
// v6: quantized_neutral_weight > 0 selects the quantized neutral release mode
// (arrival policy v1.1): every neutral leaf of the patch has weight exactly w0, the
// exact stock is debited by the leaves' exact sum n*w0 (<= stock), and the quote D
// is only the expectation (E[n w0] = D); 0 keeps the v5 exact-quote debit.
// v7: surface options per patch (benchmark by default: the v6.1 patch).
struct SourceSettings {se::Config spectrum;se::LayerPolicy layer;std::string sampling_policy_identity;double quantized_neutral_weight=0;surface_v7::SurfaceOptions surface{};};
// v7: per-patch report of the SS304 phase thermal owner (ss304-phase-thermal-owner-v1-20261031, INTEGRATION.txt 2.1):
// area-weighted surface temperature, maximum, liquid fraction, melt depth, recession and cumulative evaporated atoms.
struct PatchPhase {double mean_T=0,max_T=0,liquid_fraction=0,melt_depth_m=0,recession_m=0;std::array<double,3>evaporated_atoms{};};   // Fe, Cr, Ni
// Const view of the accepted phase-owner state (the coordinator adapts thermal.patch_reports()); bound by revision,
// time and digest; the owner reads every surface temperature from it when set.
class PhaseView {
public:
 virtual ~PhaseView()=default;
 virtual std::uint64_t revision()const=0;virtual double time_s()const=0;virtual std::string sha256()const=0;
 virtual std::size_t patches()const=0;virtual PatchPhase patch(std::size_t j)const=0;
};
// v7 owner-wide options (default: v6.1).
struct OwnerOptions {
 surface_v7::IntervalPolicy interval{};
 double emission_T_max_K=3000;     // > 3000: the extended kernel (emission_hot_v7) above 3000 K, <= 5000 K (API-v7.txt 10)
 double water_law_T_max_K=3000;    // water desorption/adsorption law domain (phase owner INTEGRATION 2.2); stock above: typed refusal
 double water_writeoff_T_K=0;      // > 0: a sub-w0 quantized residue is written off explicitly at T >= this (before 3000 K)
 double richardson_lambda=1;       // multiplies every model-mode emission current (scenario; 1: ideal free-electron A0)
 const PhaseView*phase=nullptr;    // temperature source (not part of equality; its presence is)
 bool operator==(const OwnerOptions&o)const{return interval==o.interval&&emission_T_max_K==o.emission_T_max_K&&water_law_T_max_K==o.water_law_T_max_K&&water_writeoff_T_K==o.water_writeoff_T_K&&richardson_lambda==o.richardson_lambda&&(phase!=nullptr)==(o.phase!=nullptr);}
};

// v7 first-hit wall event (from the Hall/native sink record of the previous step).
struct WallHit {std::uint64_t carrier_id=0;std::uint32_t species=0;std::size_t patch=0;double weight=0;V3 proper_velocity_m_s{};V3 hit_point_m{};double hit_time_s=0;};
// v7 particle born by a wall event: same weight as the incident macro-particle.
struct WallBirth {std::uint64_t parent_carrier=0;std::uint32_t ordinal=0,species=0,channel=0;std::size_t patch=0;double weight=0;V3 position_m{},proper_velocity_m_s{};double kinetic_eV=0,parent_hit_time_s=0;};
// v7 emission birth policy (sample_birth) and result.
struct BirthPolicy {double subcell_extent_m=0,return_extent_m=0,source_log_current_timescale_s=std::numeric_limits<double>::max();double image_force_fraction=.03,image_energy_fraction=.03,maximum_acceleration_time_fraction=.05,maximum_straight_projection_gyroangle=.05;};
enum class BirthOutcome {Handoff,ReturnedSubcell};
struct BirthPhase {BirthOutcome outcome=BirthOutcome::Handoff;se::ConditionalPacket packet;V3 displacement_m{};double turning_distance_m=INFINITY,return_time_s=0;se::Ledger cycle;};
// v7 phase map entry (one enhancement class of one patch).
struct PhaseEntry {std::size_t enhancement_class=0;double beta=1,area_fraction=1,macroscopic_field_V_m=0,local_field_V_m=0,work_function_eV=0,coverage=0,T_K=0;std::string regime,refusal;double log_current_A_m2=-INFINITY,current_A_m2=0,barrier_lowering_eV=0,fermi_decay_width_eV=0,kT_over_d=0,above_barrier_fraction=0,mean_metal_energy_eV=0;bool in_domain=false;};
// v7 interval evidence: two accepted rate samples (the caller's history); the owner recomputes the bound.
struct IntervalEvidence {surface_v7::RateSample previous,current;};
struct NeutralLeaf {double physical_weight=0;eh::Draw draw{};}; // explicit authored leaves; never adjusted
struct ReleaseQuote {std::size_t patch_index=0;PatchIdentity identity;Count released;double log_released=-INFINITY;double T_K=0,fraction=0,hazard=0,log_hazard=0,reference_eV=0;double stimulated_hazard=0;};
class PreparedPatch {
 friend class Owner;
 struct Data {v3::PreparedPatch native;std::weak_ptr<const char>life;std::uint64_t generation=0;std::shared_ptr<const Field>field;PatchIdentity identity;se::LayerPolicy policy;std::string sampling_policy,owner_signature,content;std::size_t enhancement_class=0;double beta=1,macroscopic_field_V_m=0,work_function_eV=0,coverage=0;bool v7=false;};
 std::shared_ptr<const Data>d_;explicit PreparedPatch(Data d):d_(std::make_shared<const Data>(std::move(d))){}
public:
 const std::string&exposure_key()const{return d_->native.exposure_key;}double expected_number()const{return d_->native.native.expected_physical_electrons;}double log_expected_number()const{return d_->native.log_expected_number;}const PatchIdentity&identity()const{return d_->identity;}double temperature_K()const{return d_->native.native.stage.T_K;}
 double mean_wall_heat_diagnostic_J()const{return d_->native.native.mean_wall_heat_diagnostic_J;}
 // v6.1: the prepared flux (stage, spectrum, expected count), e.g. for surface_emission::sample.
 const se::PreparedFlux&flux()const{return d_->native.native;}
 double log_current_A_m2()const{return d_->native.native.spectrum.log_current_A_m2;}double mean_metal_energy_eV()const{return d_->native.native.spectrum.mean_metal_energy_eV;}
 // v7: enhancement class, its beta, the macroscopic one-sided extraction field and the work function used.
 std::size_t enhancement_class()const{return d_->enhancement_class;}double beta()const{return d_->beta;}double macroscopic_field_V_m()const{return d_->macroscopic_field_V_m;}double work_function_eV()const{return d_->work_function_eV;}double coverage()const{return d_->coverage;}
};
// v7 wall response of one step (pure; bound to the owner state it was computed on).
class WallResponse {
 friend class Owner;
public:
 struct Data {std::weak_ptr<const char>life;std::uint64_t generation=0;std::string owner_signature,source_base;double epoch_t=0,epoch_dt=0;std::uint64_t epoch_revision=0,seed=0,step=0;
  std::vector<SW7Hit>hits;std::vector<SW7Out>outs;std::vector<double>hit_time;std::vector<V3>hit_point;std::vector<WallBirth>births;
  std::vector<std::map<std::uint64_t,surface_v7::Sums>>sums;std::vector<Count>stuck;std::vector<double>hazard;std::vector<eh::Transfer>loads;std::string content;};
private:
 std::shared_ptr<const Data>d_;explicit WallResponse(Data d):d_(std::make_shared<const Data>(std::move(d))){}
public:
 WallResponse()=default;
 const std::vector<WallBirth>&births()const{return d_->births;}const std::vector<SW7Hit>&hits()const{return d_->hits;}const std::vector<SW7Out>&outputs()const{return d_->outs;}
 const std::map<std::uint64_t,surface_v7::Sums>&sums(std::size_t j)const{return d_->sums.at(j);}const Count&stuck(std::size_t j)const{return d_->stuck.at(j);}double hazard(std::size_t j)const{return d_->hazard.at(j);}
 const std::vector<eh::Transfer>&loads()const{return d_->loads;}const std::string&content_sha256()const{return d_->content;}std::uint64_t seed()const{return d_->seed;}std::uint64_t step()const{return d_->step;}
 surface_v7::Sums total()const{surface_v7::Sums s;for(const auto&m:d_->sums)for(const auto&[w,x]:m)s.add(x);return s;}
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
 const V3&displacement_m()const{return d_->native.displacement_m;}double kinetic_change_J()const{return d_->native.kinetic_change_J;}
};
struct LayerRecord {eh::BoundLayerState state;PatchIdentity identity;se::LayerPolicy policy;double slice_center_z_m=0,extrusion_height_m=0;bool qualification=false;std::string energy_owner;};
class Trial {
 friend class Owner;
 struct Data {std::vector<eh::Transfer> loads;std::uint64_t thermal_revision=0;double thermal_time=0;std::string thermal_signature;std::shared_ptr<const Field>field;std::vector<Candidate>sources;std::vector<Handoff>handoffs;bool qualification=false;std::vector<Count>inventory;std::map<std::uint64_t,LayerRecord>layers;std::vector<eh::Packet>packets;std::weak_ptr<const char>life;std::uint64_t generation=0;std::string source_base,content;eh::Budget reported_budget;std::array<double,3>neutral_cumulative{};double neutral_K_J=0,internal_J=0,adsorbed_reference_J=0,expected_neutral_K_J=0,neutral_relativistic_conversion_J=0;std::shared_ptr<const WallResponse::Data>wall;std::string interval_record;std::vector<Count>writeoff;std::vector<PatchPhase>phase;std::string phase_sha;};
 std::unique_ptr<Data>d_;explicit Trial(std::unique_ptr<Data>d):d_(std::move(d)){}
public:
 Trial(Trial&&)=default;Trial&operator=(Trial&&)=default;Trial(const Trial&)=delete;Trial&operator=(const Trial&)=delete;
 // v7: the wall response carried by this trial (null without one) and the interval receipt text.
 const WallResponse::Data*wall()const{return d_->wall.get();}const std::string&interval_record()const{return d_->interval_record;}
 // v7: water residue written off per patch (empty without write-offs); phase-owner reports the trial was built on
 // (metal vapor: cumulative evaporated atoms per element -- owned and paid by the phase owner, never a source load).
 const std::vector<Count>&water_writeoff()const{return d_->writeoff;}const std::vector<PatchPhase>&phase_reports()const{return d_->phase;}
 const std::vector<eh::Packet>&packets()const{return d_->packets;}const Count&inventory(std::size_t j)const{return d_->inventory.at(j);}const eh::Budget&heat_budget()const{return d_->reported_budget;}const std::vector<eh::Transfer>&loads()const{return d_->loads;}std::uint64_t thermal_revision()const{return d_->thermal_revision;}const std::string&content_sha256()const{return d_->content;}double neutral_kinetic_J()const{return d_->neutral_K_J;}
};
class Owner {
 eh::Geometry geometry_;eh::Material material_;std::vector<eh::Patch>patches_;std::vector<PatchIdentity>identities_;std::vector<SourceSettings>settings_;double initial_T_,start_;std::string source_hash_,signature_;
 const eh::Coupler*heat_=nullptr;std::vector<Count>inventory_;std::map<std::uint64_t,LayerRecord>layers_;
 std::array<double,3>neutral_cumulative_{};
 std::shared_ptr<const char>life_=std::make_shared<const char>(0);std::uint64_t generation_=0;
 OwnerOptions options_;bool v7_identity_=false; // v7: owner options; true when any v7 option differs from v6.1
 // v6.1: state digest cache; the digest depends only on members that change with generation_
 // (inventory, layers, neutral cumulative) and on the thermal clock (revision, time).
 struct StateCache {bool valid=false;std::uint64_t generation=0,revision=0,time_bits=0;std::string digest;};mutable StateCache state_cache_;
 static std::vector<eh::Patch>heat_only(std::vector<eh::Patch>p){for(auto&x:p){x.desorption_enabled=false;x.coverage_temkin=false;x.initial_particles=x.initial_low_particles=0;}return p;}
 // v6: the thermal view must be a Coupler of exactly this geometry, material, initial
 // temperature, source hash and heat-only patch list (same configuration signature).
 const eh::Coupler*checked_view(const eh::Coupler&view)const{eh::Coupler shadow(geometry_,material_,heat_only(patches_),initial_T_,source_hash_,start_);
  eh::need(view.config_signature()==shadow.config_signature(),"v6 thermal view is not the declared heat-only slice of this source configuration");return &view;}
 bool live(const std::weak_ptr<const char>&w,std::uint64_t g)const{return w.lock()==life_&&g==generation_;}
 void epoch(const eh::Epoch&e)const{eh::need(e.t==heat_->state().time&&e.revision==heat_->state().revision&&e.mesh_revision==geometry_.revision&&e.material_hash==material_.data_hash&&e.source_hash==source_hash_&&std::isfinite(e.dt)&&e.dt>0&&e.t+e.dt>e.t,"v6 accepted epoch mismatch");}
 static std::string candidate_content(const Candidate::Data&d){detail::Encoding e;e.s("source-v6-candidate");e.s(d.owner_signature);e.s(d.prepared_content);e.s(d.field->content_sha256());detail::identity(e,d.identity);detail::policy(e,d.policy);e.u(d.patch_index);e.d(d.slice_center_z_m);e.d(d.extrusion_height_m);e.u(d.birth);e.u(d.native.canonical_particle_id);e.s(d.native.exposure_key);e.s(d.native.count_owner_token);e.u(d.native.qualification_only);detail::layer(e,d.native.state);detail::exchange(e,d.native.exchange);detail::transfer(e,d.native.heat);e.u(d.heat_events.size());for(const auto&t:d.heat_events)detail::transfer(e,t);e.u(d.native.paths.size());for(const auto&p:d.native.paths){e.doubles(p.origin_m);e.doubles(p.from_relative_m);e.doubles(p.to_relative_m);e.doubles(p.displacement_m);for(double x:{p.start_s,p.end_s,p.charge_C,p.field_work_J})e.d(x);}if(d.native.birth_sample){e.u(1);const auto&s=*d.native.birth_sample;for(double x:{s.physical_weight,s.normal_uniform,s.transverse_uniform,s.azimuth_uniform,s.time_uniform,s.specified_birth_s})e.d(x);e.u(s.event_id);e.u(s.canonical_particle_id);e.u(s.count_owner_revision);e.s(s.count_owner_base_hash);e.s(s.count_owner_token);e.s(s.exposure_key);e.u(int(s.time_policy));e.u(int(s.weight_provenance));}else e.u(0);return e.sha();}
 void prepared(const PreparedPatch&p)const{eh::need(p.d_&&live(p.d_->life,p.d_->generation)&&p.d_->owner_signature==signature_,"v6 stale/foreign prepared source");}
 void candidate(const Candidate&c,bool qualification=false)const{eh::need(c.d_&&live(c.d_->life,c.d_->generation)&&c.d_->owner_signature==signature_&&c.d_->content==candidate_content(*c.d_)&&v3::preflight(c.d_->native,*heat_,c.d_->field->native_,qualification),"v6 stale/foreign/changed source candidate");}

 // ---- identities carried by typed refusals
 std::string patch_text(std::size_t j)const{const auto&i=identities_.at(j);std::ostringstream s;s<<"patch-index="<<j<<";patch="<<patches_.at(j).id<<";surface="<<i.physical_surface_id<<";layer="<<i.material_layer_id;return s.str();}
 std::string stock_text(std::size_t j,const Count&n)const{detail::Encoding e;e.s("source-v6-exact-stock");detail::count(e,n);return patch_text(j)+";accepted-revision="+std::to_string(heat_->state().revision)+";exact-stock-sha256="+e.sha();}
 std::string particle_text(std::size_t j,std::uint64_t id,const std::string&event,double w)const{return patch_text(j)+";canonical-particle="+std::to_string(id)+";event="+event+";weight="+detail::hex(w);}
 static std::string event_text(const eh::SourceEventIdentity&id){return id.producer+"/"+id.energy_owner+"/"+std::to_string(id.sequence);}

 // ---- v5 weighted admission of one sealed kernel evaluation (API-v5.txt)
 struct Deposit {foil::Vec J,Q,R;};
 // The same per-path deposit operations as the sealed AxialTrialField::import.
 static Deposit deposit(const Field&f,std::span<const v3::PrecisionPath>paths){
  const auto&g=*f.mesh_;const double dt=f.frame().dt_s;Deposit d{foil::Vec(g.ne),foil::Vec(g.nn),foil::Vec(g.nn)};
  for(const auto&p:paths){axial::Path path(g,p.absolute_from(),p.displacement_m);path.deposit(d.J,p.charge_C,p.displacement_m,dt);axial::point_deposit(g,d.Q,p.absolute_from(),-p.charge_C);axial::point_deposit(g,d.Q,p.absolute_to(),p.charge_C);}
  auto div=g.gradient_transpose(d.J);for(int n=0;n<g.nn;++n)d.R[n]=d.Q[n]-dt*div[n];return d;}
 static Particle canonical_for(const v3::Candidate&c,double w){Particle p;p.id=c.canonical_particle_id;p.mass_kg=se::me;p.charge_C=-se::qe;p.weight=w;return p;}
 static v3::HandoffTrial remainder(const Field&f,const v3::Candidate&c,const Particle&canonical,double w,bool qualification){auto s=c;s.state.native.particle.weight=w;auto p=canonical;p.weight=w;return f.native_->handoff_remainder(s,p,qualification);}
 static void predicted_handoff(const detail::Weighted&W,const v3::HandoffTrial&unit,const v3::HandoffTrial&tilde){
  W.predicted_array("ordinary handoff current load C m/s",tilde.current.current_load_C_m_s);
  W.predicted_array("ordinary handoff endpoint charge impulse C",tilde.current.endpoint_charge_impulse_C);
  W.predicted("ordinary handoff work J",unit.current.field_work_J,tilde.current.field_work_J);
  W.predicted("ordinary handoff claimed work J",unit.current.claimed_field_work_J,tilde.current.claimed_field_work_J);
  W.predicted("ordinary handoff work defect J",unit.current.work_defect_J,tilde.current.work_defect_J);
  W.predicted("ordinary handoff continuity C",unit.current.continuity_C,tilde.current.continuity_C);
  W.predicted("ordinary handoff kinetic increment J",unit.kinetic_change_J,tilde.kinetic_change_J);
  const double K=raikou::mfem_backend::pic_kinetic_energy_J(tilde.particle);W.predicted("ordinary handoff endpoint kinetic J",K,K*W.tilde());}
 static void consistent_handoff(const detail::Weighted&W,const v3::HandoffTrial&unit,const v3::HandoffTrial&t,const v3::HandoffTrial&a){
  const auto&x=a.particle,&y=t.particle;
  W.same("ordinary handoff phase space",x.id==y.id&&x.position_m==y.position_m&&x.proper_velocity_m_s==y.proper_velocity_m_s&&x.mass_kg==y.mass_kg&&x.charge_C==y.charge_C&&x.weight==W.w()&&y.weight==W.tilde());
  W.same("ordinary handoff path and times",a.boundary.surface==t.boundary.surface&&a.boundary.fraction==t.boundary.fraction&&a.displacement_m==t.displacement_m&&a.handoff_time_s==t.handoff_time_s&&a.remaining_duration_s==t.remaining_duration_s&&a.endpoint_time_s==t.endpoint_time_s&&a.current.zero_duration_bridges==t.current.zero_duration_bridges);
  W.consistent_array("ordinary handoff current load C m/s",t.current.current_load_C_m_s,a.current.current_load_C_m_s);
  W.consistent_array("ordinary handoff endpoint charge impulse C",t.current.endpoint_charge_impulse_C,a.current.endpoint_charge_impulse_C);
  W.consistent("ordinary handoff work J",unit.current.field_work_J,t.current.field_work_J,a.current.field_work_J);
  W.consistent("ordinary handoff claimed work J",unit.current.claimed_field_work_J,t.current.claimed_field_work_J,a.current.claimed_field_work_J);
  W.consistent("ordinary handoff work defect J",unit.current.work_defect_J,t.current.work_defect_J,a.current.work_defect_J);
  W.consistent("ordinary handoff continuity C",unit.current.continuity_C,t.current.continuity_C,a.current.continuity_C);
  W.consistent("ordinary handoff kinetic increment J",unit.kinetic_change_J,t.kinetic_change_J,a.kinetic_change_J);
  const double K=raikou::mfem_backend::pic_kinetic_energy_J(a.particle);W.consistent("ordinary handoff endpoint kinetic J",K,K*W.tilde(),K*W.w());}
 // Ordinary-handoff remainder at the actual weight, admitted only when every
 // weighted current/charge array entry, continuity, work and energy value is the
 // exact binary rescaling of the remainder at the binade-normalized weight.
 v3::HandoffTrial handoff_trial(const v3::Candidate&raw,const Field&f,const Particle&canonical,bool qualification,const detail::Weighted&W)const{
  const auto&n=raw.state.native;
  eh::need(canonical.id==raw.canonical_particle_id&&canonical.mass_kg==se::me&&canonical.charge_C==-se::qe&&canonical.weight==n.particle.weight,"canonical/source electron species, ID or weight mismatch");
  auto unit=remainder(f,raw,canonical,1,qualification);
  std::optional<v3::HandoffTrial>tilde;if(!W.trivial()){tilde=remainder(f,raw,canonical,W.tilde(),qualification);predicted_handoff(W,unit,*tilde);}
  const double scale=std::abs(unit.current.field_work_J)+std::abs(unit.kinetic_change_J)+raikou::mfem_backend::pic_kinetic_energy_J(unit.particle);
  auto actual=W.at_weight([&]{return f.native_->handoff_remainder(raw,canonical,qualification);},scale);
  if(!W.trivial())consistent_handoff(W,unit,*tilde,actual);
  const double w=canonical.weight;
  detail::require_product(unit.kinetic_change_J,w,actual.kinetic_change_J,"ordinary handoff kinetic increment J",W.identity());
  detail::require_product(unit.current.field_work_J,w,actual.current.field_work_J,"ordinary handoff work J",W.identity());
  const double K=raikou::mfem_backend::pic_kinetic_energy_J(actual.particle);detail::require_product(K,w,K*w,"ordinary handoff endpoint kinetic J",W.identity());
  return actual;}
 // Before the actual-weight evaluation: every weighted value of the binade-normalized
 // evaluation must rescale exactly to the actual weight.
 // birth: the particle ledger is produced by this evaluation (an advance carries it unchanged).
 void predict(const detail::Weighted&W,const v3::Candidate&unit,const v3::Candidate&tilde,const Field&f,std::size_t first,bool birth)const{
  if(W.trivial())return;
  detail::predicted_values(W,detail::operand_names,detail::operands(unit.exchange),detail::operands(tilde.exchange));
  W.predicted("source charge C",-se::qe,-se::qe*W.tilde());
  const double K=se::kinetic_energy_J(tilde.state.native.particle.proper_velocity_m_s);W.predicted("source endpoint kinetic J",K,K*W.tilde());
  if(birth)detail::predicted_values(W,detail::ledger_names,detail::ledger_values(unit.state.native.particle.ledger),detail::ledger_values(tilde.state.native.particle.ledger));
  detail::predicted_values(W,detail::heat_names,detail::heat_values(unit.heat),detail::heat_values(tilde.heat));
  const bool paired=unit.paths.size()==tilde.paths.size();
  for(std::size_t k=first;k<tilde.paths.size();++k){const auto&p=tilde.paths[k];W.predicted("source path charge C",paired?unit.paths[k].charge_C:NAN,p.charge_C);W.predicted("source path work J",paired?unit.paths[k].field_work_J:NAN,p.field_work_J);}
  if(first<tilde.paths.size()){auto d=deposit(f,std::span(tilde.paths).subspan(first));W.predicted_array("source path current load C m/s",d.J);W.predicted_array("source path endpoint charge impulse C",d.Q);W.predicted_array("source path continuity residual C",d.R);}
  if(tilde.state.native.status==se::LayerStatus::PICHandedOff){const auto c=canonical_for(tilde,W.tilde());auto u=remainder(f,tilde,c,1,true);auto t=remainder(f,tilde,c,W.tilde(),true);predicted_handoff(W,u,t);}
 }
 // After the actual-weight evaluation: every weighted value must equal that rescaling exactly.
 void confirm(const detail::Weighted&W,const v3::Candidate&unit,const v3::Candidate&tilde,const v3::Candidate&actual,const Field&f,std::size_t first,bool birth)const{
  if(W.trivial())return;
  const auto&a=actual.state.native,&t=tilde.state.native;
  W.same("source state phase space",a.status==t.status&&a.current_time_s==t.current_time_s&&a.position_from_surface_m==t.position_from_surface_m&&a.particle.proper_velocity_m_s==t.particle.proper_velocity_m_s&&a.particle.normal_metal_energy_eV==t.particle.normal_metal_energy_eV&&a.particle.transverse_energy_eV==t.particle.transverse_energy_eV&&a.return_overlap_distance_m==t.return_overlap_distance_m&&a.handoff_distance_m==t.handoff_distance_m&&a.particle.weight==W.w()&&t.particle.weight==W.tilde());
  detail::consistent_values(W,detail::operand_names,detail::operands(unit.exchange),detail::operands(tilde.exchange),detail::operands(actual.exchange));
  W.consistent("source charge C",-se::qe,-se::qe*W.tilde(),-se::qe*W.w());
  const double K=se::kinetic_energy_J(a.particle.proper_velocity_m_s);W.consistent("source endpoint kinetic J",K,K*W.tilde(),K*W.w());
  if(birth)detail::consistent_values(W,detail::ledger_names,detail::ledger_values(unit.state.native.particle.ledger),detail::ledger_values(t.particle.ledger),detail::ledger_values(a.particle.ledger));
  detail::consistent_values(W,detail::heat_names,detail::heat_values(unit.heat),detail::heat_values(tilde.heat),detail::heat_values(actual.heat));
  W.same("source path count",actual.paths.size()==tilde.paths.size());
  const bool paired=unit.paths.size()==tilde.paths.size();
  for(std::size_t k=first;k<actual.paths.size();++k){const auto&p=actual.paths[k],&q=tilde.paths[k];
   W.same("source path geometry",p.origin_m==q.origin_m&&p.from_relative_m==q.from_relative_m&&p.to_relative_m==q.to_relative_m&&p.displacement_m==q.displacement_m&&p.start_s==q.start_s&&p.end_s==q.end_s);
   W.consistent("source path charge C",paired?unit.paths[k].charge_C:NAN,q.charge_C,p.charge_C);W.consistent("source path work J",paired?unit.paths[k].field_work_J:NAN,q.field_work_J,p.field_work_J);}
  if(first<actual.paths.size()){auto dt_=deposit(f,std::span(tilde.paths).subspan(first)),da=deposit(f,std::span(actual.paths).subspan(first));
   W.consistent_array("source path current load C m/s",dt_.J,da.J);W.consistent_array("source path endpoint charge impulse C",dt_.Q,da.Q);W.consistent_array("source path continuity residual C",dt_.R,da.R);}
  if(a.status==se::LayerStatus::PICHandedOff)(void)handoff_trial(actual,f,canonical_for(actual,W.w()),true,W);
 }
 // The w~ counterpart of an admitted candidate (its weighted values are exact rescalings).
 static v3::Candidate scaled_down(const v3::Candidate&c,const detail::Weighted&W){
  auto x=c;const int e=W.scale().e;
  auto down=[&](double&v,const char*q){double r=std::scalbn(v,-e);W.same(q,std::scalbn(r,e)==v);v=r;};
  auto&s=x.exchange;for(double*v:{&s.physical_electron_number_change,&s.particle_charge_change_C,&s.electrode_charge_change_C,&s.wall_heat_J,&s.particle_K_change_J,&s.electron_reference_change_J,&s.image_energy_change_J,&s.field_work_J})down(*v,"admitted birth exchange rescaling");
  for(auto&p:x.paths){down(p.charge_C,"admitted birth path rescaling");down(p.field_work_J,"admitted birth path rescaling");}
  auto&l=x.state.native.particle.ledger;for(double*v:{&l.physical_electron_number,&l.particle_charge_C,&l.electrode_charge_C,&l.wall_heat_J,&l.outgoing_K_J,&l.incoming_K_J,&l.electron_reference_change_J,&l.image_reference_change_J,&l.resolved_matching_field_work_J})down(*v,"admitted birth ledger rescaling");
  x.state.native.particle.weight=W.tilde();return x;}

 Candidate seal(v3::Candidate c,const PreparedPatch&p,bool birth,const std::string&id)const{double w=c.state.native.particle.weight;detail::require_product(-se::qe,w,-se::qe*w,"source charge C",id);double K=se::kinetic_energy_J(c.state.native.particle.proper_velocity_m_s);detail::require_product(K,w,K*w,"source endpoint kinetic J",id);for(const auto&path:c.paths){auto E=p.d_->field->native_->sample(path).E;double unitW=-se::qe*v3::dot(E,path.displacement_m);detail::require_product(unitW,w,path.field_work_J,"source path work J",id);detail::require_product(-se::qe,w,path.charge_C,"source path charge C",id);}Candidate::Data d;d.heat_events.push_back(c.heat);d.native=std::move(c);d.life=life_;d.generation=generation_;d.field=p.d_->field;d.identity=p.d_->identity;d.policy=p.d_->policy;d.patch_index=p.d_->native.patch_index;d.slice_center_z_m=p.d_->native.binding.extrusion_center_z_m;d.extrusion_height_m=p.d_->native.binding.extrusion_height_m;d.birth=birth;d.owner_signature=signature_;d.prepared_content=p.d_->content;d.content=candidate_content(d);return Candidate(std::move(d));}
 static std::string handoff_content(const Handoff::Data&d){detail::Encoding e;e.s(d.source_content);e.s(d.field_content);const auto&p=d.native.particle;e.u(p.id);e.doubles(p.position_m);e.doubles(p.proper_velocity_m_s);for(double x:{p.mass_kg,p.charge_C,p.weight,d.native.handoff_time_s,d.native.remaining_duration_s,d.native.endpoint_time_s,d.native.kinetic_change_J})e.d(x);e.u(d.native.boundary.surface);e.d(d.native.boundary.fraction);e.doubles(d.native.displacement_m);e.doubles(d.native.current.current_load_C_m_s);e.doubles(d.native.current.endpoint_charge_impulse_C);for(double x:{d.native.current.field_work_J,d.native.current.claimed_field_work_J,d.native.current.work_defect_J,d.native.current.continuity_C})e.d(x);e.u(d.native.current.zero_duration_bridges);return e.sha();}
 std::vector<ReleaseQuote>quotes(const eh::Epoch&e,const WallResponse*wall=nullptr)const;
 // v7 helpers
 double coverage_of(std::size_t j)const{const auto&p=patches_[j];if(!p.desorption_enabled||!(p.monolayer_m2>0))return 0;return inventory_[j].rounded_for_rate()/(identities_[j].physical_area_m2*p.monolayer_m2);}
 // v7 temperature source: the phase-owner patch report when bound, else the thermal view's patch mean.
 double T_of(std::size_t j)const{return options_.phase?options_.phase->patch(j).mean_T:heat_->surface_temperatures().at(j);}
 double liquid_fraction_of(std::size_t j)const{if(!options_.phase)return 0;const double f=options_.phase->patch(j).liquid_fraction;eh::need(std::isfinite(f)&&f>=0&&f<=1,"v7 phase report liquid fraction outside [0, 1]");return f;}
 void phase_bound()const{if(options_.phase)eh::need(options_.phase->time_s()==heat_->state().time&&options_.phase->patches()==patches_.size(),"v7 phase-owner report does not match the accepted clock or patch layout");}
 // water law domain (phase owner INTEGRATION 2.2): a nonzero water stock above water_law_T_max_K is a typed refusal
 void water_domain(std::size_t j)const{const double T=T_of(j);if(!inventory_[j].zero()&&T>options_.water_law_T_max_K)throw Refusal("adsorbed water above the desorption-law domain",std::log(inventory_[j].rounded_for_rate()),patch_text(j)+";T="+detail::hex(T),RefusalKind::surface_domain,std::numeric_limits<double>::quiet_NaN(),0,"v7: nonzero water inventory reported above "+std::to_string(options_.water_law_T_max_K)+" K",T,options_.water_law_T_max_K);}
 double work_function_of(std::size_t j,const se::Config&c)const{const auto&m=settings_[j].surface.work_function;return m.benchmark?c.work_function_eV:surface_v7::work_function_eV(m,coverage_of(j),T_of(j));}
 double liquid_work_function_of(std::size_t j)const{const auto&m=settings_[j].surface.work_function;eh::need(!m.benchmark&&std::isfinite(m.phi_liquid_eV),"v7 liquid work function not configured");const double q=std::nearbyint(m.phi_liquid_eV/m.grid_eV)*m.grid_eV;eh::need(q>=2&&q<=6,"v7 liquid work function outside 2..6 eV");return q;}
 // class k of patch j: configured (solid) classes scaled by (1 - liquid fraction); with a phase view one more
 // class (index K) for the molten area: beta 1, area = liquid fraction, liquid work function.
 bool liquid_class(std::size_t j,std::size_t k)const{return options_.phase&&!settings_[j].surface.benchmark()&&k==settings_[j].surface.enhancement.classes.size();}
 double adsorbed_reference_eV(std::size_t j)const{const auto&p=patches_[j];return p.desorption_eV*(1-(p.coverage_temkin?p.binding_alpha*coverage_of(j):0.));}
 SW7Params wall_parameters(std::size_t j)const;
 V3 patch_midpoint(std::size_t j)const;
 void field_cap(const Field&f)const{eh::need(f.interval_cap_s()==options_.interval.cap_s,"v7 field trial interval bound differs from the owner configuration");}
 void admit_exposure(const v3::PreparedPatch&n,std::size_t j,const v3::Exposure&x)const;
 static std::string trial_content(const Trial::Data&d){detail::Encoding e;e.s("source-v6-trial-loads");e.s(d.source_base);e.s(d.field->content_sha256());e.u(d.qualification);for(const auto&c:d.sources)e.s(c.content_sha256());for(const auto&h:d.handoffs)e.s(h.content_sha256());for(const auto&n:d.inventory)detail::count(e,n);for(const auto&[id,r]:d.layers){e.u(id);detail::layer(e,r.state);detail::identity(e,r.identity);detail::policy(e,r.policy);e.d(r.slice_center_z_m);e.d(r.extrusion_height_m);e.u(r.qualification);e.s(r.energy_owner);}for(const auto&p:d.packets){e.s(p.event_id);e.s(p.state_id);e.s(p.parent_id);e.s(p.channel_id);e.s(p.patch);e.u(p.source_revision);for(double x:{p.weight,p.T,p.time,p.mass_kg,p.charge_C,p.kinetic_J,p.internal_J,p.reference_J,p.adsorbed_reference_change_J,p.classical_sample_K_J,p.relativistic_conversion_J})e.d(x);e.doubles(p.x);e.doubles(p.v);e.doubles(p.proper_velocity_m_s);e.doubles(p.wall_support_impulse_N_s);}detail::budget(e,d.reported_budget);e.u(d.loads.size());for(const auto&t:d.loads)detail::transfer(e,t);e.u(d.thermal_revision);e.d(d.thermal_time);e.s(d.thermal_signature);e.d(d.neutral_K_J);e.d(d.internal_J);e.d(d.adsorbed_reference_J);e.d(d.expected_neutral_K_J);e.d(d.neutral_relativistic_conversion_J);e.doubles(d.neutral_cumulative);if(d.wall){e.s("v7-wall-response");e.s(d.wall->content);}if(!d.interval_record.empty()){e.s("v7-interval");e.s(d.interval_record);}
  if(!d.writeoff.empty()){e.s("v7-water-writeoff");for(const auto&c:d.writeoff)detail::count(e,c);}
  if(!d.phase_sha.empty()){e.s("v7-phase-view");e.s(d.phase_sha);for(const auto&r:d.phase){for(double x:{r.mean_T,r.max_T,r.liquid_fraction,r.melt_depth_m,r.recession_m})e.d(x);e.doubles(r.evaporated_atoms);}}return e.sha();}
public:
 static constexpr std::size_t maximum_leaves_per_patch=4096;
 Owner(eh::Geometry g,eh::Material m,std::vector<eh::Patch>p,std::vector<PatchIdentity>ids,std::vector<SourceSettings>settings,double T,std::string source,double start,const eh::Coupler&thermal_view,OwnerOptions options={}):geometry_(std::move(g)),material_(std::move(m)),patches_(std::move(p)),identities_(std::move(ids)),settings_(std::move(settings)),initial_T_(T),start_(start),source_hash_(std::move(source)),options_(std::move(options)){
  eh::need(geometry_.nx>0&&geometry_.ny>0&&geometry_.nx<=1024&&geometry_.ny<=1024&&patches_.size()<=256,"v6 component resource geometry/patch cap");eh::Coupler validation(geometry_,material_,patches_,T,source_hash_,start);eh::need(identities_.size()==patches_.size()&&settings_.size()==patches_.size(),"v6 physical surface/layer identity required per patch");detail::Encoding enc;enc.s("heat-source-owner-v6");enc.s("thermal-view;no-private-engine;no-Joule-input");enc.s(validation.config_signature());detail::configuration(enc,geometry_,material_,patches_,T,start,source_hash_);enc.s("46e01fe22e67cc9deb0da9320edef4041f4b3ef433de8ce1eca484463b405501");for(std::size_t j=0;j<patches_.size();++j){const auto&i=identities_[j];eh::need(!i.physical_surface_id.empty()&&!i.material_layer_id.empty()&&std::isfinite(i.physical_area_m2)&&i.physical_area_m2>0&&i.physical_area_m2==validation.area(j),"v6 physical patch area/layer identity mismatch");detail::identity(enc,i);const auto&setting=settings_[j];eh::need(!setting.sampling_policy_identity.empty(),"v6 explicit external sampling policy identity required");v3::admit_policy(setting.layer);enc.d(setting.spectrum.work_function_eV);enc.u(setting.spectrum.energy_panels);enc.u(setting.spectrum.action_nodes);enc.d(setting.spectrum.tail_widths);detail::policy(enc,setting.layer);enc.s(setting.sampling_policy_identity);eh::need(setting.quantized_neutral_weight==0||(std::isfinite(setting.quantized_neutral_weight)&&setting.quantized_neutral_weight>=DBL_MIN),"v6 quantized neutral weight must be 0 or a positive normal binary64");enc.s("quantized-neutral-release-v1");enc.d(setting.quantized_neutral_weight);inventory_.emplace_back(patches_[j].initial_particles,patches_[j].initial_low_particles);}
  // v7 configuration records, appended only when an option differs from v6.1 (benchmark identity unchanged).
  {const auto&iv=options_.interval;eh::need(std::isfinite(iv.cap_s)&&iv.cap_s>0&&iv.cap_s<=v3::hard_trial_interval_cap_s&&iv.emission_relative_error>=0&&iv.desorption_relative_error>=0&&iv.safety>0&&iv.safety<=1&&iv.count_floor>=0,"v7 interval policy outside its domain");
   v7_identity_=!(options_==OwnerOptions{});for(const auto&s:settings_)v7_identity_=v7_identity_||!s.surface.benchmark();
   eh::need(options_.emission_T_max_K>=3000&&options_.emission_T_max_K<=5000&&options_.water_law_T_max_K>0&&options_.water_law_T_max_K<=3000&&options_.water_writeoff_T_K>=0&&(options_.water_writeoff_T_K==0||options_.water_writeoff_T_K<=options_.water_law_T_max_K)&&std::isfinite(options_.richardson_lambda)&&options_.richardson_lambda>0&&options_.richardson_lambda<=1,"v7 owner options outside their domains (T_max 3000..5000 K, water law <= 3000 K, write-off <= law cap, 0 < lambda <= 1)");
   if(v7_identity_){enc.s("heat-source-owner-v7");enc.s("surface-physics-v7;interval-bound-v7");for(double x:{iv.cap_s,iv.emission_relative_error,iv.desorption_relative_error,iv.safety,iv.count_floor})enc.d(x);
    enc.s("high-temperature-v7");for(double x:{options_.emission_T_max_K,options_.water_law_T_max_K,options_.water_writeoff_T_K,options_.richardson_lambda})enc.d(x);enc.u(options_.phase!=nullptr);
    for(std::size_t j=0;j<settings_.size();++j){const auto&s=settings_[j].surface;surface_v7::validate(s.enhancement);const auto&w=s.work_function;
     if(!w.benchmark)eh::need(std::isfinite(w.phi_ox_eV)&&std::isfinite(w.dphi_1ML_eV)&&std::isfinite(w.dphi_dT_eV_K)&&w.depolarization_c>=0&&w.grid_eV>0&&w.band_low_eV<=w.phi_ox_eV&&w.phi_ox_eV<=w.band_high_eV,"v7 work-function model invalid");
     enc.u(w.benchmark);for(double x:{w.fixed_phi_eV,w.phi_ox_eV,w.dphi_1ML_eV,w.depolarization_c,w.dphi_dT_eV_K,w.T_ref_K,w.grid_eV,w.band_low_eV,w.band_high_eV,w.phi_liquid_eV,w.band_liquid_low_eV,w.band_liquid_high_eV})enc.d(x);
     enc.u(s.enhancement.classes.size());for(const auto&c:s.enhancement.classes){enc.d(c.beta);enc.d(c.area_fraction);enc.d(c.local_extent_m);enc.s(c.label);}
     enc.u(s.wall.enabled);enc.u(s.wall.hall_applies_incident_kinetic);enc.s(surface_v7::bytes_of(s.wall.p));enc.s(surface_v7::bytes_of(s.wall.species));
     enc.u(s.photo.enabled);enc.doubles(s.photo.E_eV);enc.doubles(s.photo.yield);enc.d(s.photo.scale);}}}
  signature_=enc.sha();heat_=checked_view(thermal_view);}
 Owner(const Owner&)=delete;Owner&operator=(const Owner&)=delete;Owner(Owner&&)=delete;Owner&operator=(Owner&&)=delete;
 const eh::State&thermal_state()const{return heat_->state();}const std::array<double,3>&neutral_cumulative_J()const{return neutral_cumulative_;}const Count&inventory(std::size_t j)const{return inventory_.at(j);}const std::map<std::uint64_t,LayerRecord>&active_layers()const{return layers_;}const PatchIdentity&patch_identity(std::size_t j)const{return identities_.at(j);}eh::Epoch make_epoch(double dt,const std::string&field)const{return heat_->make_epoch(dt,field);}const std::string&configuration_sha256()const{return signature_;}
 std::string state_sha256()const{const auto&st=heat_->state();const auto tb=std::bit_cast<std::uint64_t>(st.time);
  if(state_cache_.valid&&state_cache_.generation==generation_&&state_cache_.revision==st.revision&&state_cache_.time_bits==tb)return state_cache_.digest;
  auto digest=state_sha256_uncached();state_cache_={true,generation_,st.revision,tb,digest};return digest;}
 // v6 state digest (no cache).
 std::string state_sha256_uncached()const{detail::Encoding e;e.s(signature_);e.doubles(neutral_cumulative_);const auto&s=heat_->state();e.d(s.time);e.u(s.revision);for(const auto&n:inventory_)detail::count(e,n);for(const auto&[id,r]:layers_){e.u(id);detail::layer(e,r.state);detail::identity(e,r.identity);detail::policy(e,r.policy);e.d(r.slice_center_z_m);e.d(r.extrusion_height_m);e.u(r.qualification);e.s(r.energy_owner);}return e.sha();}
 std::vector<ReleaseQuote>release_quotes(const eh::Epoch&e)const{return quotes(e);}
 PreparedPatch prepare_patch(const eh::Epoch&e,std::size_t j,v3::PatchBinding binding,std::shared_ptr<const Field>field,V3 origin,const v3::Exposure&exposure)const;
 // v7: exposure of enhancement class k of patch j (class 0 = prepare_patch). The class's area fraction
 // multiplies exposure.area_fraction; the kernel sees beta_k times the one-sided macroscopic field.
 PreparedPatch prepare_class(const eh::Epoch&e,std::size_t j,std::size_t k,v3::PatchBinding binding,std::shared_ptr<const Field>field,V3 origin,const v3::Exposure&exposure)const;
 std::size_t enhancement_classes(std::size_t j)const{const auto&s=settings_.at(j).surface;return s.enhancement.classes.size()+(options_.phase&&!s.benchmark()?1:0);}
 // area fraction of class k at the accepted state (0: inactive; the caller skips it)
 double class_area_fraction(std::size_t j,std::size_t k)const{const auto&s=settings_.at(j).surface;const double fl=s.benchmark()?0:liquid_fraction_of(j);if(liquid_class(j,k))return fl;eh::need(k<s.enhancement.classes.size(),"v7 enhancement class index");return s.enhancement.classes[k].area_fraction*(1-fl);}
 // v7: emission phase map of patch j at the given macroscopic one-sided field and temperature (pure).
 std::vector<PhaseEntry>phase_map(std::size_t j,double macroscopic_field_V_m,double T_K)const;
 // v7: emitted phase space of one birth at the macroscopic matching plane (pure).
 BirthPhase sample_birth(const PreparedPatch&p,double weight,double normal_uniform,double transverse_uniform,double azimuth_uniform,std::uint64_t event_id,const BirthPolicy&policy)const;
 // v7: wall responses of first hits (pure; bound to the current accepted state).
 WallResponse respond(const eh::Epoch&e,std::span<const WallHit>hits,std::uint64_t seed,std::uint64_t step)const;
 std::vector<ReleaseQuote>release_quotes(const eh::Epoch&e,const WallResponse&wall)const{return quotes(e,&wall);}
 // v7: rate sample of the accepted state (log expected births per second per (patch, class); log thermal
 // release per second per patch), for the caller's interval history; and the bound from two samples.
 surface_v7::RateSample rate_sample(const eh::Epoch&e,const std::shared_ptr<const Field>&field,std::span<const V3>origins)const;
 surface_v7::IntervalDecision interval_bound(const IntervalEvidence&x)const{return surface_v7::interval_bound(x.previous,x.current,options_.interval);}
 const OwnerOptions&options()const{return options_;}bool v7_identity()const{return v7_identity_;}double surface_work_function_eV(std::size_t j)const{return work_function_of(j,settings_.at(j).spectrum);}double surface_coverage(std::size_t j)const{return coverage_of(j);}
 Candidate begin_source(const PreparedPatch&p,const v3::BirthSample&sample,const eh::SourceEventIdentity&id)const;
 Candidate advance_active(const PreparedPatch&p,std::uint64_t canonical_id,const eh::SourceEventIdentity&id,const v3::AdvanceControls&controls={})const;
 Candidate advance_birth(const PreparedPatch&p,const Candidate&birth,const eh::SourceEventIdentity&continuation_event,const v3::AdvanceControls&controls={})const;
 bool preflight(const Candidate&c,const std::shared_ptr<const Field>&field,const FinalFieldView&actual,bool qualification=false)const noexcept{try{candidate(c,qualification);return c.d_->field==field&&field&&field->matches(actual);}catch(...){return false;}}
 foil::SurfacePathReceipt import(const std::shared_ptr<const Field>&field,const FinalFieldView&actual,std::span<const Candidate>sources,const std::set<std::uint64_t>&excluded,bool qualification=false)const{eh::need(field&&field->matches(actual),"v6 actual final field/mesh changed");std::vector<v3::Candidate>raw;raw.reserve(sources.size());for(const auto&c:sources){candidate(c,qualification);eh::need(c.d_->field==field,"v6 source final field identity refused");raw.push_back(c.d_->native);}return field->native_->import(raw,excluded,qualification);}
 Handoff handoff(const Candidate&c,const FinalFieldView&actual,const Particle&canonical,bool qualification=false)const{
  candidate(c,qualification);eh::need(c.d_->field->matches(actual),"v6 actual handoff field/mesh changed");
  eh::need(std::isfinite(canonical.weight)&&canonical.weight>0,"canonical/source electron species, ID or weight mismatch");
  const detail::Weighted W(canonical.weight,particle_text(c.d_->patch_index,c.canonical_id(),"ordinary-handoff",canonical.weight));
  Handoff::Data d;d.native=handoff_trial(c.d_->native,*c.d_->field,canonical,qualification,W);d.source_content=c.d_->content;d.field_content=c.d_->field->content_sha256();d.content=handoff_content(d);return Handoff(std::move(d));}
 Trial prepare(const eh::Epoch&e,const std::shared_ptr<const Field>&field,const FinalFieldView&actual,std::span<const Candidate>sources,std::span<const Handoff>handoffs,const std::vector<std::vector<NeutralLeaf>>&leaves,const std::vector<eh::Transfer>&other={},bool qualification=false){return prepare(e,field,actual,sources,handoffs,leaves,other,qualification,nullptr,nullptr);}
 // v7: with the step's wall response (stuck molecules credited, stimulated release, wall loads) and the
 // interval evidence (required when the interval policy has error targets).
 Trial prepare(const eh::Epoch&e,const std::shared_ptr<const Field>&field,const FinalFieldView&actual,std::span<const Candidate>sources,std::span<const Handoff>handoffs,const std::vector<std::vector<NeutralLeaf>>&leaves,const std::vector<eh::Transfer>&other,bool qualification,const WallResponse*wall,const IntervalEvidence*interval);
 bool preflight(const Trial&t,const FinalFieldView&actual)const noexcept;
 bool commit(Trial&&t,const FinalFieldView&actual)noexcept;
 void publish_prevalidated(Trial&&t)noexcept; // sole coordinator, after all-owner preflight in one critical section
 void save(const std::string&path)const;void load(const std::string&path);
};

inline std::vector<ReleaseQuote>Owner::quotes(const eh::Epoch&e,const WallResponse*wall)const {
 phase_bound();
 if(wall)eh::need(wall->d_&&wall->d_->life.lock()==life_&&wall->d_->generation==generation_&&wall->d_->owner_signature==signature_&&wall->d_->epoch_t==e.t&&wall->d_->epoch_dt==e.dt&&wall->d_->epoch_revision==e.revision,"v7 wall response does not bind this accepted state and epoch");
 epoch(e);std::vector<ReleaseQuote>out;
 for(std::size_t j=0;j<patches_.size();++j){const auto&p=patches_[j];ReleaseQuote q;q.patch_index=j;q.identity=identities_[j];
  if(!p.desorption_enabled||inventory_[j].zero()){out.push_back(q);continue;}
  double N=inventory_[j].rounded_for_rate(),T=T_of(j);q.T_K=T;double theta=0;water_domain(j);
  if(p.coverage_temkin){eh::need(T>=p.source_T_min&&T<=p.source_T_max,"remaining water outside explicitly selected source-law domain");theta=N/(q.identity.physical_area_m2*p.monolayer_m2);eh::need(theta>=0&&theta<=1,"coverage outside submonolayer model");}
  double prefactor=p.attempt_s*(p.prefactor==eh::Prefactor::EyringScaled?T/p.prefactor_reference_T:1.);double barrier=p.activation_eV*(1-p.rate_alpha*theta);
  q.log_hazard=std::log(prefactor)+std::log(e.dt)-barrier*eh::qe/(eh::kb*T);
  // v7: electron/ion-stimulated desorption hazard of the step's wall response (first-order, same stock).
  if(wall&&wall->d_->hazard[j]>0){q.stimulated_hazard=wall->d_->hazard[j];q.log_hazard=std::log(std::exp(q.log_hazard)+q.stimulated_hazard);}
  double released;
  if(q.log_hazard<-36){q.log_released=std::log(N)+q.log_hazard;released=std::exp(q.log_released);q.fraction=std::exp(q.log_hazard);}
  else{q.hazard=q.log_hazard>700?INFINITY:std::exp(q.log_hazard);q.fraction=-std::expm1(-q.hazard);released=N*q.fraction;q.log_released=std::log(N)+std::log(q.fraction);}
  if(!(std::isfinite(released)&&released>0))throw Refusal("neutral release count",q.log_released,stock_text(j,inventory_[j]),RefusalKind::release_scalar,released==0?INFINITY:std::numeric_limits<double>::quiet_NaN());
  // v5: a rounded release above the exact stock is a typed, scoped refusal (v4: untyped
  // "insufficient exact inventory"). The exact stock is never rounded or retreated.
  if(!detail::covers(inventory_[j],released))throw Refusal("neutral release count",q.log_released,stock_text(j,inventory_[j]),RefusalKind::stock_drain,std::numeric_limits<double>::quiet_NaN(),0,"rounded release "+detail::hex(released)+" exceeds the exact accepted stock (rate view "+detail::hex(N)+"); no stock or weight adjustment",released,N);
  q.released=Count(released);auto remaining=inventory_[j].debit(released); // exact, no nextafter stock/weight adjustment
  double theta_after=p.coverage_temkin?remaining.rounded_for_rate()/(q.identity.physical_area_m2*p.monolayer_m2):0;
  q.reference_eV=p.desorption_eV*(1-p.binding_alpha*.5*(theta+theta_after));out.push_back(std::move(q));
 }
 return out;
}
// v5 finite-log exposure admission (API-v5.txt): the expected count scalar and its
// wall-heat diagnostic must lie within relative_tolerance of the exact value defined
// by the kernel's finite log current, else a typed refusal carries the log scope.
inline void Owner::admit_exposure(const v3::PreparedPatch&n,std::size_t j,const v3::Exposure&x)const{
 if(!n.has_log_expected_number)return; // no positive flux defined (T=0 and F<=0): zero is the exact value
 const auto&spectrum=n.native.spectrum;const double area=heat_->area(j),fraction=x.area_fraction,dt=x.duration_s,count=n.native.expected_physical_electrons,L=n.log_expected_number;
 using LD=long double;
 const LD log_current=spectrum.log_current_A_m2;
 const LD log_exact=log_current+std::log(LD(area))+std::log(LD(fraction))+std::log(LD(dt))-std::log(LD(se::qe));
 const double scope=(n.native.stage.area_m2>=DBL_MIN&&std::isfinite(L))?L:double(log_exact);
 std::ostringstream id;id<<patch_text(j)<<";exposure-begin="<<detail::hex(x.begin_s)<<";duration="<<detail::hex(dt)<<";area-fraction="<<detail::hex(fraction)<<";quadrature="<<x.quadrature_id<<";log-current="<<detail::hex(spectrum.log_current_A_m2);
 const LD exact=std::exp(log_current)*LD(area)*LD(fraction)*LD(dt)/LD(se::qe);
 if(!(count>0)||!std::isfinite(L))throw Refusal("electron expected exposure count",scope,id.str(),RefusalKind::exposure_scalar,count==0?1.:std::numeric_limits<double>::quiet_NaN(),relative_tolerance,"finite log count with an unrepresented scalar",count,double(exact));
 const LD rel=exact>0&&std::isfinite(double(exact))?std::fabs((LD(count)-exact)/exact):LD(INFINITY);
 if(!(rel<=LD(admission_bound)))throw Refusal("electron expected exposure count",scope,id.str(),RefusalKind::exposure_scalar,double(rel),relative_tolerance,"scalar "+detail::hex(count)+" differs from exp(log count)",count,double(exact));
 const double diagnostic=n.native.mean_wall_heat_diagnostic_J;
 const LD heat=-std::exp(log_current)*LD(spectrum.mean_metal_energy_eV)*LD(area)*LD(fraction)*LD(dt);
 if(heat==0){if(diagnostic!=0)throw Refusal("electron mean wall heat exposure diagnostic J",double(log_exact),id.str(),RefusalKind::exposure_scalar,INFINITY,relative_tolerance,"nonzero diagnostic for a zero exact value",diagnostic,0.);}
 else{const LD rh=std::fabs((LD(diagnostic)-heat)/heat);if(!(rh<=LD(admission_bound)))throw Refusal("electron mean wall heat exposure diagnostic J",double(std::log(std::fabs(heat))),id.str(),RefusalKind::exposure_scalar,double(rh),relative_tolerance,"diagnostic "+detail::hex(diagnostic),diagnostic,double(heat));}
}
inline PreparedPatch Owner::prepare_patch(const eh::Epoch&e,std::size_t j,v3::PatchBinding binding,std::shared_ptr<const Field>field,V3 origin,const v3::Exposure&exposure)const {
 return prepare_class(e,j,0,std::move(binding),std::move(field),origin,exposure);
}
inline PreparedPatch Owner::prepare_class(const eh::Epoch&e,std::size_t j,std::size_t k,v3::PatchBinding binding,std::shared_ptr<const Field>field,V3 origin,const v3::Exposure&exposure_in)const {
 const auto&setting=settings_.at(j);auto config=setting.spectrum;const auto&policy=setting.layer;std::string sampling_policy=setting.sampling_policy_identity;epoch(e);eh::need(field&&!sampling_policy.empty(),"v6 immutable field and declared sampling-policy identity required");v3::admit_policy(policy);eh::need(field->frame().accepted_revision==e.revision&&field->frame().begin_s==e.t&&field->frame().dt_s==e.dt,"v6 common source/field clock mismatch");
 field_cap(*field);
 // v7 surface state and enhancement class
 const auto&surface=setting.surface;const bool bench=surface.benchmark();phase_bound();if(!bench)water_domain(j);
 const bool liquid=liquid_class(j,k);eh::need(liquid||k<surface.enhancement.classes.size(),"v7 enhancement class index");
 const surface_v7::EnhancementClass cls=liquid?surface_v7::EnhancementClass{1.,1.,INFINITY,"liquid"}:surface.enhancement.classes[k];
 const double fraction=bench?1.:class_area_fraction(j,k);eh::need(fraction>0,"v7 class has zero area at the accepted state (skip it)");
 config.work_function_eV=liquid?liquid_work_function_of(j):work_function_of(j,setting.spectrum);auto exposure=exposure_in;if(fraction!=1)exposure.area_fraction=exposure_in.area_fraction*fraction;
 const double T_over=options_.phase?T_of(j):std::numeric_limits<double>::quiet_NaN();
 PreparedPatch::Data d;d.native=detail::prepare_patch(*heat_,e,j,std::move(binding),field->native_,origin,exposure,config,cls.beta,T_over,options_.emission_T_max_K,bench?1.:options_.richardson_lambda);d.life=life_;d.generation=generation_;d.field=std::move(field);d.identity=identities_.at(j);d.policy=policy;d.sampling_policy=std::move(sampling_policy);d.owner_signature=signature_;
 d.enhancement_class=k;d.beta=cls.beta;d.macroscopic_field_V_m=se::extraction_field(d.native.native.stage.total_E,d.native.native.stage.outward_normal);d.work_function_eV=config.work_function_eV;d.coverage=coverage_of(j);d.v7=!bench;
 admit_exposure(d.native,j,exposure);
 if(!bench){std::ostringstream s;s<<std::hexfloat<<"/v7-class/"<<k<<'/'<<cls.beta<<'/'<<cls.area_fraction<<'/'<<config.work_function_eV<<'/'<<d.coverage;d.native.exposure_key+=s.str();}
 detail::Encoding enc;enc.s(d.native.exposure_key);enc.s(signature_);enc.s(d.field->content_sha256());detail::identity(enc,d.identity);detail::policy(enc,policy);enc.s(d.sampling_policy);enc.d(d.native.log_expected_number);enc.d(d.native.native.expected_physical_electrons);enc.d(d.native.native.mean_wall_heat_diagnostic_J);
 if(!bench){enc.s("v7-class");enc.u(k);enc.d(d.beta);enc.d(d.macroscopic_field_V_m);enc.d(d.work_function_eV);enc.d(d.coverage);enc.u(liquid);enc.d(fraction);}
 if(options_.phase){enc.s("v7-phase-view");enc.s(options_.phase->sha256());}
 d.content=enc.sha();d.native.exposure_key+="/v6-immutable/"+d.content;return PreparedPatch(std::move(d));
}
inline Candidate Owner::begin_source(const PreparedPatch&p,const v3::BirthSample&sample,const eh::SourceEventIdentity&id)const {
 prepared(p);eh::need(p.d_->beta==1,"v7: the time-resolved image layer is admitted for the beta = 1 class only (enhanced classes: sample_birth)");eh::need(std::isfinite(sample.physical_weight)&&sample.physical_weight>0,"v6 explicit positive representable electron weight required");eh::need(!layers_.contains(sample.canonical_particle_id),"v6 birth ID already owns accepted source layer");
 const detail::Weighted W(sample.physical_weight,particle_text(p.d_->native.patch_index,sample.canonical_particle_id,"birth/"+std::to_string(sample.event_id)+"/"+event_text(id),sample.physical_weight));
 auto at=[&](double w){auto s=sample;s.physical_weight=w;return v3::begin_source(p.d_->native,p.d_->policy,s,id);};
 // Identical marks; weight 1 and the binade-normalized weight are numerical admission
 // evaluations only, never replacement sampling.
 auto unit=at(1);detail::source_exchange_represented(unit.exchange,unit.exchange,1,W.identity());
 std::optional<v3::Candidate>tilde;if(!W.trivial()){tilde=at(W.tilde());predict(W,unit,*tilde,*p.d_->field,0,true);}
 detail::admit_weighted_exchange(unit.exchange,sample.physical_weight,W.identity());
 auto raw=W.at_weight([&]{return v3::begin_source(p.d_->native,p.d_->policy,sample,id);},detail::energy_scale(unit.exchange));
 if(!W.trivial())confirm(W,unit,*tilde,raw,*p.d_->field,0,true);
 detail::source_exchange_represented(unit.exchange,raw.exchange,sample.physical_weight,W.identity());return seal(std::move(raw),p,true,W.identity());
}
inline Candidate Owner::advance_active(const PreparedPatch&p,std::uint64_t id,const eh::SourceEventIdentity&event,const v3::AdvanceControls&controls)const {
 prepared(p);const auto&record=layers_.at(id);eh::need(event.energy_owner==record.energy_owner,"v6 accepted canonical energy owner changed");eh::need(record.identity==p.d_->identity&&record.slice_center_z_m==p.d_->native.binding.extrusion_center_z_m&&record.extrusion_height_m==p.d_->native.binding.extrusion_height_m&&detail::policy_hash(record.policy)==detail::policy_hash(p.d_->policy),"v6 accepted layer area/material/policy changed");
 const double w=record.state.native.particle.weight;const detail::Weighted W(w,particle_text(p.d_->native.patch_index,id,"advance/"+event_text(event),w));
 auto at=[&](double weight){auto s=record.state;s.native.particle.weight=weight;return v3::advance_source(p.d_->native,s,id,p.d_->policy,p.d_->field->frame().dt_s,event,controls);};
 auto check=at(1);
 std::optional<v3::Candidate>tilde;if(!W.trivial()){tilde=at(W.tilde());predict(W,check,*tilde,*p.d_->field,0,false);}
 detail::admit_weighted_exchange(check.exchange,w,W.identity());
 auto raw=W.at_weight([&]{return v3::advance_source(p.d_->native,record.state,id,p.d_->policy,p.d_->field->frame().dt_s,event,controls);},detail::energy_scale(check.exchange));
 if(!W.trivial())confirm(W,check,*tilde,raw,*p.d_->field,0,false);
 detail::source_exchange_represented(check.exchange,raw.exchange,w,W.identity());raw.qualification_only=record.qualification;return seal(std::move(raw),p,false,W.identity());
}
inline Candidate Owner::advance_birth(const PreparedPatch&p,const Candidate&birth,const eh::SourceEventIdentity&event,const v3::AdvanceControls&controls)const {
 prepared(p);candidate(birth,true);eh::need(event.energy_owner==birth.aggregate_exchange_receipt().energy_owner,"v6 newborn canonical energy owner changed");eh::need(birth.d_->birth&&birth.d_->field==p.d_->field&&birth.d_->prepared_content==p.d_->content&&birth.state().native.status==se::LayerStatus::Active,"v6 same-trial Active birth required");
 double duration=p.d_->field->frame().end_s-birth.state().native.current_time_s;eh::need(duration>0,"v6 birth has no active remaining interval");
 const double w=birth.state().native.particle.weight;const detail::Weighted W(w,particle_text(birth.d_->patch_index,birth.canonical_id(),"continuation/"+event_text(event),w));
 auto at=[&](double weight){auto s=birth.state();s.native.particle.weight=weight;return v3::advance_source(p.d_->native,s,birth.canonical_id(),p.d_->policy,duration,event,controls);};
 auto check=at(1);
 std::optional<v3::Candidate>tilde;if(!W.trivial()){tilde=at(W.tilde());predict(W,check,*tilde,*p.d_->field,0,false);}
 detail::admit_weighted_exchange(check.exchange,w,W.identity());
 auto move=W.at_weight([&]{return v3::advance_source(p.d_->native,birth.state(),birth.canonical_id(),p.d_->policy,duration,event,controls);},detail::energy_scale(check.exchange));
 if(!W.trivial())confirm(W,check,*tilde,move,*p.d_->field,0,false);
 detail::source_exchange_represented(check.exchange,move.exchange,w,W.identity());
 auto raw=birth.d_->native;raw.state=move.state;v3::add(raw.exchange,move.exchange);raw.paths.insert(raw.paths.end(),move.paths.begin(),move.paths.end());
 if(!W.trivial()){
  // The combined birth+continuation at w~: the admitted birth rescaled down exactly, plus the w~ continuation.
  auto combined=scaled_down(birth.d_->native,W);combined.state=tilde->state;v3::add(combined.exchange,tilde->exchange);combined.paths.insert(combined.paths.end(),tilde->paths.begin(),tilde->paths.end());
  v3::finalize(combined,p.d_->native,event,combined.state.native.current_time_s);
  W.at_weight([&]{v3::finalize(raw,p.d_->native,event,raw.state.native.current_time_s);return 0;},detail::energy_scale(check.exchange));
  detail::consistent_values(W,detail::operand_names,detail::operands(check.exchange),detail::operands(combined.exchange),detail::operands(raw.exchange));
  detail::consistent_values(W,detail::heat_names,detail::heat_values(check.heat),detail::heat_values(combined.heat),detail::heat_values(raw.heat));
 }else v3::finalize(raw,p.d_->native,event,raw.state.native.current_time_s);
 auto combined=seal(std::move(raw),p,true,W.identity());auto data=*combined.d_;data.heat_events=birth.d_->heat_events;data.heat_events.push_back(move.heat);data.content=candidate_content(data);return Candidate(std::move(data));
}
inline Trial Owner::prepare(const eh::Epoch&e,const std::shared_ptr<const Field>&field,const FinalFieldView&actual,std::span<const Candidate>sources,std::span<const Handoff>handoffs,const std::vector<std::vector<NeutralLeaf>>&leaves,const std::vector<eh::Transfer>&other,bool qualification,const WallResponse*wall,const IntervalEvidence*interval){
 if(field)field_cap(*field);
 phase_bound();
 std::string interval_record;
 if(options_.interval.emission_relative_error>0||options_.interval.desorption_relative_error>0){eh::need(interval&&interval->current.t_s==e.t&&interval->previous.t_s<interval->current.t_s,"v7 interval evidence (two accepted rate samples ending at this epoch) required by the interval policy");
  const auto b=interval_bound(*interval);std::ostringstream s;s<<std::hexfloat<<"bound="<<b.bound_s<<";limiting="<<b.limiting<<";rate="<<b.rate_per_s<<";index="<<b.index<<";dt="<<e.dt;interval_record=s.str();
  if(!(e.dt<=b.bound_s))throw Refusal("source interval bound",std::log(e.dt),interval_record,RefusalKind::interval_bound,e.dt/b.bound_s,1,"v7 step exceeds the source-owner interval bound",e.dt,b.bound_s);}
 epoch(e);eh::need(field&&field->frame().accepted_revision==e.revision&&field->frame().begin_s==e.t&&field->frame().dt_s==e.dt,"v6 common source/field clock mismatch");eh::need(sources.size()<=4096&&handoffs.size()<=4096,"v6 source proposal cap; whole proposal refused");eh::need(generation_<UINT64_MAX&&field&&field->matches(actual)&&leaves.size()==patches_.size(),"v6 common field/mesh or leaf dimensions invalid");
 auto d=std::make_unique<Trial::Data>();d->life=life_;d->generation=generation_;d->source_base=state_sha256();d->field=field;d->sources.assign(sources.begin(),sources.end());d->handoffs.assign(handoffs.begin(),handoffs.end());d->qualification=qualification;d->inventory=inventory_;d->layers=layers_;auto transfers=other;
 std::set<std::uint64_t>seen;std::map<std::uint64_t,const Handoff*>hand;for(const auto&h:handoffs){eh::need(h.d_&&h.d_->content==handoff_content(*h.d_)&&hand.emplace(h.particle().id,&h).second,"v6 changed/duplicate handoff");}
 for(const auto&c:sources){candidate(c,qualification);eh::need(c.d_->field==field&&seen.insert(c.canonical_id()).second,"v6 source final identity/duplicate");auto id=c.canonical_id();eh::need(c.d_->birth?!layers_.contains(id):layers_.contains(id),"v6 source birth/accepted-layer ownership mismatch");
  if(c.state().native.status==se::LayerStatus::Active){eh::need(c.state().native.current_time_s==field->frame().end_s&&!hand.contains(id),"v6 Active layer not at common endpoint");d->layers[id]={c.state(),c.d_->identity,c.d_->policy,c.d_->slice_center_z_m,c.d_->extrusion_height_m,c.d_->native.qualification_only,c.aggregate_exchange_receipt().energy_owner};}
  else{d->layers.erase(id);if(c.state().native.status==se::LayerStatus::PICHandedOff){auto it=hand.find(id);eh::need(it!=hand.end()&&it->second->d_->source_content==c.content_sha256()&&it->second->d_->field_content==field->content_sha256(),"v6 missing or foreign ordinary handoff remainder");}else eh::need(!hand.contains(id),"v6 recollected source also handed off");}
  transfers.insert(transfers.end(),c.d_->heat_events.begin(),c.d_->heat_events.end());
 }
 std::set<std::string>energy_owners;for(const auto&[id,r]:d->layers)eh::need(!r.energy_owner.empty()&&energy_owners.insert(r.energy_owner).second,"v6 multiple Active IDs share an energy owner");for(const auto&[id,_]:layers_)eh::need(seen.contains(id),"v6 accepted Active layer omitted from interval");for(const auto&[id,_]:hand)eh::need(seen.contains(id),"v6 orphan handoff");
 auto q=quotes(e,wall);for(std::size_t j=0;j<patches_.size();++j){const auto&p=patches_[j];const auto&quote=q[j];eh::need(leaves[j].size()<=maximum_leaves_per_patch,"v6 neutral leaf cap; whole proposal refused");Count total;for(const auto&leaf:leaves[j]){eh::need(std::isfinite(leaf.physical_weight)&&leaf.physical_weight>0,"v6 explicit positive neutral weight required");total=total.credited(leaf.physical_weight);}
  const double w0=settings_[j].quantized_neutral_weight;double reference_eV=quote.reference_eV;
  if(w0>0){
   // quantized release: leaves of exactly w0, exact debit n*w0 within the exact stock
   for(const auto&leaf:leaves[j])eh::need(leaf.physical_weight==w0,"v6 quantized release: every neutral leaf has weight exactly w0");
   if(total.zero())continue;
   {Count rest=inventory_[j];try{for(double w:total.positive_carriers())rest=rest.debit(w);}catch(const std::exception&){
     throw Refusal("neutral release count",std::log(total.rounded_for_rate()),stock_text(j,inventory_[j]),RefusalKind::stock_drain,std::numeric_limits<double>::quiet_NaN(),0,
       "quantized release n*w0 exceeds the exact accepted stock; no stock or weight adjustment",total.rounded_for_rate(),inventory_[j].rounded_for_rate());}}
   for(double w:total.positive_carriers())d->inventory[j]=d->inventory[j].debit(w);
   if(p.coverage_temkin){const double N=inventory_[j].rounded_for_rate(),A=identities_[j].physical_area_m2;
     const double theta=N/(A*p.monolayer_m2),theta_after=d->inventory[j].rounded_for_rate()/(A*p.monolayer_m2);
     reference_eV=p.desorption_eV*(1-p.binding_alpha*.5*(theta+theta_after));}
  }else{
   eh::need(total==quote.released,"v6 exact newborn leaf sum differs from exact source debit");if(total.zero())continue;
   for(double w:total.positive_carriers())d->inventory[j]=d->inventory[j].debit(w);
  }
  for(std::size_t k=0;k<leaves[j].size();++k){const auto&leaf=leaves[j][k];const auto&a=leaf.draw;for(double u:{a.normal,a.radius,a.angle,a.along,a.z,a.time})eh::need(std::isfinite(u)&&u>0&&u<1,"v6 caller-owned open neutral draws required");eh::Packet packet;packet.patch=p.id;packet.state_id=p.state_id;packet.parent_id=p.parent_id;packet.channel_id=p.coverage_temkin?"conditional-oxide-water-Temkin-frozen-coverage":"authored-first-order-desorption";packet.source_revision=e.revision+1;packet.event_id="v6/desorption/"+identities_[j].physical_surface_id+"/"+identities_[j].material_layer_id+"/"+p.id+"/"+std::to_string(packet.source_revision)+"/"+std::to_string(k);packet.T=quote.T_K;packet.weight=leaf.physical_weight;
   double along=p.begin+(p.end-p.begin)*a.along;if(p.face==eh::Face::LowY||p.face==eh::Face::HighY)packet.x={geometry_.x0+along*geometry_.width,p.face==eh::Face::LowY?geometry_.y0:geometry_.y0+geometry_.thickness,(a.z-.5)*geometry_.height};else packet.x={p.face==eh::Face::LowX?geometry_.x0:geometry_.x0+geometry_.width,geometry_.y0+along*geometry_.thickness,(a.z-.5)*geometry_.height};
   double vn=std::sqrt(-2*eh::kb*packet.T/p.mass_kg*std::log(a.normal)),vt=std::sqrt(-2*eh::kb*packet.T/p.mass_kg*std::log(a.radius)),angle=2*3.14159265358979323846*a.angle;if(p.face==eh::Face::LowY||p.face==eh::Face::HighY)packet.v={vt*std::cos(angle),(p.face==eh::Face::LowY?-1:1)*vn,vt*std::sin(angle)};else packet.v={(p.face==eh::Face::LowX?-1:1)*vn,vt*std::cos(angle),vt*std::sin(angle)};
   // v5: every weighted neutral value must be the exact binary rescaling of its
   // binade-normalized evaluation (same draw, same formulas).
   const detail::Weighted W(packet.weight,stock_text(j,inventory_[j])+";neutral-leaf="+std::to_string(k)+";weight="+detail::hex(packet.weight));
   auto weighted=[&](eh::Packet&x){x.internal_J=p.internal_eV*eh::qe*x.weight;x.reference_J=0;x.adsorbed_reference_change_J=reference_eV*eh::qe*x.weight;for(int axis=0;axis<3;++axis)x.wall_support_impulse_N_s[axis]=-x.weight*x.mass_kg*x.proper_velocity_m_s[axis];};
   auto values=[&](const eh::Packet&x)->std::array<double,9>{return {x.kinetic_J,x.classical_sample_K_J,x.relativistic_conversion_J,x.internal_J,x.adsorbed_reference_change_J,x.wall_support_impulse_N_s[0],x.wall_support_impulse_N_s[1],x.wall_support_impulse_N_s[2],2*eh::kb*x.T*x.weight};};
   static constexpr std::array<const char*,9>names{"neutral kinetic J","neutral classical sample kinetic J","neutral relativistic conversion J","neutral internal J","adsorbed reference J","neutral wall support impulse x N s","neutral wall support impulse y N s","neutral wall support impulse z N s","expected neutral kinetic J"};
   auto unit=packet;unit.weight=1;eh::assign_native_kinematics(unit,p.mass_kg);weighted(unit);
   std::optional<eh::Packet>tilde;if(!W.trivial()){tilde=packet;tilde->weight=W.tilde();eh::assign_native_kinematics(*tilde,p.mass_kg);weighted(*tilde);detail::predicted_values(W,names,values(unit),values(*tilde));}
   W.at_weight([&]{eh::assign_native_kinematics(packet,p.mass_kg);return 0;},unit.kinetic_J);detail::require_product(unit.kinetic_J,packet.weight,packet.kinetic_J,"neutral kinetic J",W.identity());
   weighted(packet);if(p.internal_eV>0)detail::require_log_product(std::log(p.internal_eV)+std::log(eh::qe)+std::log(packet.weight),packet.internal_J,"neutral internal J",W.identity());if(reference_eV>0)detail::require_log_product(std::log(reference_eV)+std::log(eh::qe)+std::log(packet.weight),packet.adsorbed_reference_change_J,"adsorbed reference J",W.identity());
   if(!W.trivial()){W.same("neutral phase space",packet.proper_velocity_m_s==tilde->proper_velocity_m_s&&packet.v==tilde->v&&packet.x==tilde->x);detail::consistent_values(W,names,values(unit),values(*tilde),values(packet));}
   double offset=quote.log_hazard<-36?e.dt*a.time:-std::log1p(-a.time*quote.fraction)/quote.hazard*e.dt;packet.time=e.t+offset;double energy=packet.kinetic_J+packet.internal_J+packet.adsorbed_reference_change_J;eh::need(std::isfinite(energy)&&std::isfinite(packet.time),"v6 neutral energy/time overflow");
   eh::Transfer t;t.producer="v6-neutral/"+identities_[j].physical_surface_id+"/"+identities_[j].material_layer_id+"/"+p.id+"/"+std::to_string(k);t.state_id=p.state_id;t.parent_id=p.parent_id;t.channel_id=packet.channel_id;t.energy_owner=t.producer;t.sequence=t.owner_sequence=e.revision+1;t.event_time_s=packet.time;t.x=packet.x[0];t.y=packet.x[1];t.authored_heat_J=-energy;transfers.push_back(t);d->neutral_K_J+=packet.kinetic_J;d->expected_neutral_K_J+=2*eh::kb*packet.T*packet.weight;d->neutral_relativistic_conversion_J+=packet.relativistic_conversion_J;d->internal_J+=packet.internal_J;d->adsorbed_reference_J+=packet.adsorbed_reference_change_J;eh::need(d->packets.size()<4096,"v6 total neutral proposal cap; whole proposal refused");d->packets.push_back(std::move(packet));
  }
 }
 // v7: molecules stuck by the step's wall events are credited to the exact stock (after the release debit);
 // the wall loads (energy left in the foil beyond the Hall sink's incident kinetic energy) join the loads.
 if(wall){eh::need(wall->d_&&wall->d_->life.lock()==life_&&wall->d_->generation==generation_&&wall->d_->owner_signature==signature_&&wall->d_->source_base==d->source_base&&wall->d_->epoch_t==e.t&&wall->d_->epoch_dt==e.dt&&wall->d_->epoch_revision==e.revision,"v7 wall response does not bind this accepted state and epoch");
  for(std::size_t j=0;j<patches_.size();++j)for(double w:wall->d_->stuck[j].positive_carriers())d->inventory[j]=d->inventory[j].credited(w);
  transfers.insert(transfers.end(),wall->d_->loads.begin(),wall->d_->loads.end());d->wall=wall->d_;}
 d->interval_record=std::move(interval_record);
 // v7: explicit write-off of a sub-w0 quantized water residue at T >= water_writeoff_T_K (phase owner INTEGRATION
 // 2.2: the residue must leave before the patch passes the 3000 K law cap). The molecules leave without PIC
 // particles (ledger column 'water residue written off'); their adsorbed reference energy is debited from the foil.
 if(options_.water_writeoff_T_K>0)for(std::size_t j=0;j<patches_.size();++j){const double w0=settings_[j].quantized_neutral_weight;
  if(!(w0>0)||d->inventory[j].zero()||T_of(j)<options_.water_writeoff_T_K)continue;const double rest=d->inventory[j].rounded_for_rate();if(!(rest<w0))continue;
  if(d->writeoff.empty())d->writeoff.resize(patches_.size());d->writeoff[j]=d->inventory[j];d->inventory[j]=Count();
  eh::Transfer x;x.producer="v7-water-writeoff/"+identities_[j].physical_surface_id+"/"+identities_[j].material_layer_id+"/"+patches_[j].id;x.energy_owner=x.producer;x.state_id=patches_[j].state_id;x.parent_id=patches_[j].parent_id;
  x.channel_id="surface-v7-water-writeoff";x.sequence=x.owner_sequence=e.revision+1;x.event_time_s=e.t;x.surface=true;const auto m=patch_midpoint(j);x.x=m[0];x.y=m[1];x.authored_heat_J=-adsorbed_reference_eV(j)*eh::qe*rest;transfers.push_back(x);}
 if(options_.phase){d->phase_sha=options_.phase->sha256();for(std::size_t j=0;j<patches_.size();++j)d->phase.push_back(options_.phase->patch(j));}
 eh::need(d->layers.size()<=4096,"v6 Active layer cap; whole proposal refused");
 // v6: no thermal solve here. The loads go to the thermal owner, which applies them
 // exactly once; the reported budget is the sum of the loads (no Joule term).
 d->loads=std::move(transfers);d->thermal_revision=heat_->state().revision;d->thermal_time=heat_->state().time;d->thermal_signature=heat_->config_signature();
 {auto&b=d->reported_budget;b=eh::Budget{};for(const auto&x:d->loads){
   for(double v:{x.incoming_K_J,x.outgoing_K_J,x.reference_out_J,x.solid_collision_J,x.authored_heat_J,x.nottingham_wall_J,x.work_function_reference_part_J,x.source_wall_J,x.source_particle_K_change_J,x.source_electron_reference_change_J,x.source_image_energy_change_J,x.source_field_work_J})eh::need(std::isfinite(v),"v6 nonfinite surface load operand");
   if(x.certified_source_exchange){b.source_wall_J+=x.source_wall_J;b.source_particle_K_change_J+=x.source_particle_K_change_J;b.source_electron_reference_change_J+=x.source_electron_reference_change_J;b.source_image_energy_change_J+=x.source_image_energy_change_J;b.source_field_work_J+=x.source_field_work_J;b.wall_net_J+=x.source_wall_J;continue;}
   b.incoming_K_J+=x.incoming_K_J;b.outgoing_K_J+=x.outgoing_K_J;b.reference_out_J+=x.reference_out_J;b.nottingham_wall_J+=x.nottingham_wall_J;b.solid_collision_J+=x.solid_collision_J;b.work_function_reference_part_J+=x.work_function_reference_part_J;
   const double boundary=x.complete_nottingham?x.nottingham_wall_J:x.incoming_K_J-x.outgoing_K_J-x.reference_out_J;b.wall_net_J+=boundary+x.solid_collision_J+x.authored_heat_J;}
  b.joule_J=0;b.thermal_delta_J=0;b.residual_J=0;b.roundoff_bound_J=0;b.authored_heat_J=0;for(const auto&t:other)b.authored_heat_J+=t.authored_heat_J;
  b.neutral_K_J=d->neutral_K_J;b.expected_neutral_K_J=d->expected_neutral_K_J;b.neutral_relativistic_conversion_J=d->neutral_relativistic_conversion_J;b.neutral_internal_J=d->internal_J;b.adsorbed_reference_change_J=d->adsorbed_reference_J;}
 d->neutral_cumulative=neutral_cumulative_;for(int i=0;i<3;++i){double amount=i==0?d->neutral_K_J:i==1?d->internal_J:d->adsorbed_reference_J;d->neutral_cumulative[i]+=amount;eh::need(std::isfinite(d->neutral_cumulative[i]),"v6 neutral cumulative energy overflow");}d->content=trial_content(*d);return Trial(std::move(d));
}
// ---------------------------------------------------------------- v7 surface physics (owner methods)
inline V3 Owner::patch_midpoint(std::size_t j)const{const auto&p=patches_.at(j);const double along=p.begin+(p.end-p.begin)*.5;
 if(p.face==eh::Face::LowY||p.face==eh::Face::HighY)return {geometry_.x0+along*geometry_.width,p.face==eh::Face::LowY?geometry_.y0:geometry_.y0+geometry_.thickness,0.};
 return {p.face==eh::Face::LowX?geometry_.x0:geometry_.x0+geometry_.width,geometry_.y0+along*geometry_.thickness,0.};}
// Kernel parameters of patch j at the accepted state: the patch's wall model with the surface state
// (work function, temperature, coverage, adsorbed reference) of the exact stock.
inline SW7Params Owner::wall_parameters(std::size_t j)const{
 const auto&s=settings_.at(j);eh::need(s.surface.wall.enabled,"v7 wall model not selected for this patch");SW7Params p=s.surface.wall.p;
 water_domain(j);const double fl=liquid_fraction_of(j);
 p.phi_eV=fl>=.5?liquid_work_function_of(j):work_function_of(j,s.spectrum);p.T_K=T_of(j);p.theta=coverage_of(j);   // dominant phase of the patch
 if(patches_[j].desorption_enabled){p.adsorption_eV=adsorbed_reference_eV(j);eh::need(std::isfinite(p.adsorption_eV)&&p.adsorption_eV>=0,"v7 adsorbed reference energy");}
 else{p.adsorption_eV=0;p.enable_sticking=0;p.enable_esd=0;p.enable_isd=0;}   // no water stock on this patch
 if(p.T_K>options_.water_law_T_max_K)p.enable_sticking=0;   // no adsorption above the water-law domain
 return p;}
inline WallResponse Owner::respond(const eh::Epoch&e,std::span<const WallHit>hits,std::uint64_t seed,std::uint64_t step)const{
 epoch(e);phase_bound();eh::need(hits.size()<=(std::size_t(1)<<20),"v7 wall hit cap (2^20 per step)");eh::need(step<(std::uint64_t(1)<<40),"v7 wall step index outside the stream layout");
 WallResponse::Data d;d.life=life_;d.generation=generation_;d.owner_signature=signature_;d.source_base=state_sha256();d.epoch_t=e.t;d.epoch_dt=e.dt;d.epoch_revision=e.revision;d.seed=seed;d.step=step;
 const std::size_t P=patches_.size();d.sums.resize(P);d.stuck.resize(P);d.hazard.assign(P,0.);
 // canonical order (patch, carrier id): the response does not depend on the order or partition of the input
 std::vector<std::size_t>order(hits.size());std::iota(order.begin(),order.end(),std::size_t(0));
 for(const auto&h:hits)eh::need(h.patch<P&&h.species<SW7_NSPECIES&&std::isfinite(h.weight)&&h.weight>0&&std::isfinite(h.hit_time_s)&&v3::finite(h.proper_velocity_m_s)&&v3::finite(h.hit_point_m),"v7 wall hit: patch, species, positive weight, finite kinematics required");
 std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){return std::pair(hits[a].patch,hits[a].carrier_id)<std::pair(hits[b].patch,hits[b].carrier_id);});
 {std::set<std::uint64_t>ids;for(const auto&h:hits)eh::need(ids.insert(h.carrier_id).second,"v7 wall hit: one first hit per carrier and step");}
 std::vector<std::optional<SW7Params>>params(P);
 for(std::size_t i:order){const auto&h=hits[i];const std::size_t j=h.patch;if(!params[j])params[j]=wall_parameters(j);const auto&wm=settings_[j].surface.wall;
  SW7Hit x{};x.id=h.carrier_id;x.species=h.species;x.patch=std::uint32_t(j);x.weight=h.weight;const auto n=heat_->normal(j);for(int a=0;a<3;++a){x.u[a]=h.proper_velocity_m_s[a];x.n[a]=n[a];}
  SW7Out o{};sw7::wall_event(*params[j],wm.species.data(),x,seed,step,&o);
  eh::need(!(o.status&(SW7_ST_INVALID|SW7_ST_RNG_BUDGET)),"v7 wall event refused (invalid incidence or RNG word budget); never resampled");
  eh::need(surface_v7::event_closes(x,o,wm.species),"v7 wall event exact ledger does not close");
  d.sums[j][std::bit_cast<std::uint64_t>(h.weight)].add(x,o);
  if(o.stuck)d.stuck[j]=d.stuck[j].credited(h.weight);
  for(std::uint32_t k=0;k<o.n_out;++k)d.births.push_back({h.carrier_id,k,o.sp[k],o.ch[k],j,h.weight,h.hit_point_m,{o.u[k][0],o.u[k][1],o.u[k][2]},sw7::from_quanta(o.kq[k]),h.hit_time_s});
  // energy left in the foil: wall quanta, less the incident kinetic energy that the Hall sink already applied
  const std::int64_t q=o.wall_q-(wm.hall_applies_incident_kinetic?o.k_in_q:0);
  if(q!=0){eh::Transfer t;t.producer="v7-wall/"+identities_[j].physical_surface_id+"/"+identities_[j].material_layer_id+"/"+patches_[j].id+"/"+std::to_string(h.carrier_id);t.energy_owner=t.producer;
   t.state_id=surface_v7::species_name(h.species);t.parent_id="surface-v7-wall";t.channel_id="surface-v7-wall-response";t.sequence=t.owner_sequence=e.revision+1;t.event_time_s=e.t;t.surface=true;t.x=h.hit_point_m[0];t.y=h.hit_point_m[1];
   t.authored_heat_J=h.weight*(double(q)/SW7_QUANTA_PER_EV)*se::qe;eh::need(std::isfinite(t.authored_heat_J),"v7 wall load overflow");d.loads.push_back(std::move(t));}
  d.hits.push_back(x);d.outs.push_back(o);d.hit_time.push_back(h.hit_time_s);d.hit_point.push_back(h.hit_point_m);}
 // stimulated desorption hazard per patch: sum over weight classes of w * sigma_q 2^-100 / A (weight-bit order)
 for(std::size_t j=0;j<P;++j){double H=0;for(const auto&[wb,s]:d.sums[j])if(s.sigma_q!=0)H+=std::bit_cast<double>(wb)*(double(s.sigma_q)/SW7_SIGMA_QUANTA_PER_M2)/identities_[j].physical_area_m2;eh::need(std::isfinite(H)&&H>=0,"v7 stimulated hazard");d.hazard[j]=H;}
 detail::Encoding enc;enc.s("surface-v7-wall-response");enc.s(d.owner_signature);enc.s(d.source_base);enc.d(e.t);enc.d(e.dt);enc.u(e.revision);enc.u(seed);enc.u(step);enc.u(d.hits.size());
 for(std::size_t i=0;i<d.hits.size();++i){enc.s(surface_v7::bytes_of(d.hits[i]));enc.s(surface_v7::bytes_of(d.outs[i]));enc.d(d.hit_time[i]);enc.doubles(d.hit_point[i]);}
 enc.doubles(d.hazard);for(const auto&c:d.stuck)detail::count(enc,c);for(const auto&t:d.loads)detail::transfer(enc,t);if(options_.phase){enc.s("v7-phase-view");enc.s(options_.phase->sha256());}d.content=enc.sha();
 return WallResponse(std::move(d));
}
inline std::vector<PhaseEntry>Owner::phase_map(std::size_t j,double F,double T)const{
 const auto&s=settings_.at(j);std::vector<PhaseEntry>out;auto config=s.spectrum;double phi=NAN;
 try{config.work_function_eV=work_function_of(j,s.spectrum);phi=config.work_function_eV;}catch(const std::exception&x){PhaseEntry e;e.refusal=x.what();out.push_back(e);return out;}
 for(std::size_t k=0;k<enhancement_classes(j);++k){const bool liq=liquid_class(j,k);const auto c=liq?surface_v7::EnhancementClass{1.,1.,INFINITY,"liquid"}:s.surface.enhancement.classes[k];PhaseEntry e;e.enhancement_class=k;e.beta=c.beta;e.area_fraction=s.surface.benchmark()?1.:class_area_fraction(j,k);e.macroscopic_field_V_m=F;e.local_field_V_m=F>0?F*c.beta:F;e.work_function_eV=phi;e.coverage=coverage_of(j);e.T_K=T;
  auto cfg=config;if(liq){try{cfg.work_function_eV=e.work_function_eV=liquid_work_function_of(j);}catch(const std::exception&x){e.refusal=x.what();e.regime="liquid-work-function-missing";out.push_back(std::move(e));continue;}}
  try{const auto d=evaluate_emission(T,e.local_field_V_m,cfg,options_.emission_T_max_K);e.in_domain=true;e.log_current_A_m2=d.has_log_current?d.log_current_A_m2:-INFINITY;e.current_A_m2=d.current_A_m2;e.barrier_lowering_eV=d.barrier_lowering_eV;e.fermi_decay_width_eV=d.fermi_decay_width_eV;e.mean_metal_energy_eV=d.mean_metal_energy_eV;
   const double kT=se::kb*T/se::qe;e.kT_over_d=d.fermi_decay_width_eV>0?kT/d.fermi_decay_width_eV:INFINITY;
   const double top=cfg.work_function_eV-d.barrier_lowering_eV;double above=0;for(const auto&n:d.nodes)if(n.normal_energy_eV>=top)above+=n.probability;e.above_barrier_fraction=above;
   // regime (diagnostic classification of the one source integral; never a second population)
   if(!(e.local_field_V_m>0))e.regime=T>0?"retarding-thermal":"none";
   else if(!d.has_log_current||d.current_A_m2==0)e.regime=d.has_log_current?"underflow":"none";
   else if(above>=.5)e.regime="thermionic-schottky";
   else if(e.kT_over_d<.2)e.regime="field-emission";
   else e.regime="thermal-field";
  }catch(const std::exception&x){e.refusal=x.what();e.regime="outside-kernel-domain";}
  out.push_back(std::move(e));}
 return out;
}
inline BirthPhase Owner::sample_birth(const PreparedPatch&p,double weight,double un,double ut,double ua,std::uint64_t id,const BirthPolicy&b)const{
 prepared(p);eh::need(std::isfinite(weight)&&weight>0&&id>0,"v7 birth: positive weight and event id");
 const double F=p.d_->macroscopic_field_V_m;const auto&flux=p.d_->native.native;
 // handoff view: the class spectrum with the macroscopic one-sided field at the matching plane (potential of the PIC field)
 se::PreparedFlux view=flux;view.spectrum.extraction_field_V_m=F;view.spectrum.barrier_lowering_eV=F>0?se::schottky_lowering_eV(F):0;
 if(F>0)eh::need(view.spectrum.barrier_lowering_eV<view.spectrum.config.work_function_eV*.98,"v7 birth: barrier collapse at the macroscopic field");
 se::ConditionalBoundary cb;cb.local_linear_vacuum_extent_m=b.subcell_extent_m;cb.source_log_current_timescale_s=b.source_log_current_timescale_s;cb.image_force_fraction=b.image_force_fraction;cb.image_energy_fraction=b.image_energy_fraction;cb.maximum_acceleration_time_fraction=b.maximum_acceleration_time_fraction;cb.PIC_already_applies_microscopic_single_electron_image_force=false;cb.maximum_straight_projection_gyroangle=b.maximum_straight_projection_gyroangle;
 BirthPhase r;r.packet=se::sample_conditional(view,cb,weight,un,ut,ua,id);r.displacement_m=r.packet.outward_displacement_m;r.cycle=r.packet.particle.ledger;
 if(F<0){const double K=se::kinetic_energy_J(r.packet.particle.proper_velocity_m_s);r.turning_distance_m=K/(se::qe*(-F));   // bound with any B (B does no work)
  if(r.turning_distance_m<=b.return_extent_m){r.outcome=BirthOutcome::ReturnedSubcell;
   const auto&u=r.packet.particle.proper_velocity_m_s;const auto&n=flux.stage.outward_normal;r.return_time_s=2*se::me*v3::dot(u,n)/(se::qe*(-F));   // B = 0 value (diagnostic)
   const auto c=se::collect_at_material_surface(view.spectrum.config.work_function_eV,K,weight,true);auto&l=r.cycle;
   l.physical_electron_number+=c.physical_electron_number;l.particle_charge_C+=c.particle_charge_C;l.electrode_charge_C+=c.electrode_charge_C;l.wall_heat_J+=c.wall_heat_J;l.outgoing_K_J+=c.outgoing_K_J;l.incoming_K_J+=c.incoming_K_J;l.electron_reference_change_J+=c.electron_reference_change_J;l.image_reference_change_J+=c.image_reference_change_J;l.resolved_matching_field_work_J+=c.resolved_matching_field_work_J;}}
 return r;
}
inline surface_v7::RateSample Owner::rate_sample(const eh::Epoch&e,const std::shared_ptr<const Field>&field,std::span<const V3>origins)const{
 epoch(e);eh::need(field&&origins.size()==patches_.size(),"v7 rate sample: field and one origin per patch");surface_v7::RateSample r;r.t_s=e.t;r.dt_s=e.dt;
 phase_bound();
 for(std::size_t j=0;j<patches_.size();++j){const auto&s=settings_[j];auto config=s.spectrum;config.work_function_eV=work_function_of(j,s.spectrum);const double T=T_of(j);r.T_K.push_back(T);
  const auto n=heat_->normal(j);const double F=se::extraction_field(field->surface(origins[j],n).E,n);const double A=identities_[j].physical_area_m2;
  for(std::size_t k=0;k<enhancement_classes(j);++k){const bool liq=liquid_class(j,k);const double beta=liq?1.:s.surface.enhancement.classes[k].beta;const double a=s.surface.benchmark()?1.:class_area_fraction(j,k);
   double L=-INFINITY;if(a>0){try{auto cfg=config;if(liq)cfg.work_function_eV=liquid_work_function_of(j);const auto d=evaluate_emission(T,F>0?F*beta:F,cfg,options_.emission_T_max_K);if(d.has_log_current)L=d.log_current_A_m2+std::log(A*a)-std::log(se::qe)+(s.surface.benchmark()?0.:std::log(options_.richardson_lambda));}catch(const std::exception&){L=NAN;}}r.log_emission.push_back(L);}
  const auto&p=patches_[j];double D=-INFINITY;if(p.desorption_enabled&&!inventory_[j].zero()){const double N=inventory_[j].rounded_for_rate(),th=coverage_of(j);const double pre=p.attempt_s*(p.prefactor==eh::Prefactor::EyringScaled?T/p.prefactor_reference_T:1.);D=std::log(N)+std::log(pre)-p.activation_eV*(1-p.rate_alpha*(p.coverage_temkin?th:0.))*eh::qe/(eh::kb*T);}
  r.log_desorption.push_back(D);}
 return r;
}
inline bool Owner::preflight(const Trial&t,const FinalFieldView&actual)const noexcept{try{if(t.d_&&!t.d_->phase_sha.empty()&&(!options_.phase||options_.phase->sha256()!=t.d_->phase_sha))return false;if(!t.d_||!live(t.d_->life,t.d_->generation)||t.d_->source_base!=state_sha256()||t.d_->thermal_revision!=heat_->state().revision||t.d_->thermal_time!=heat_->state().time||t.d_->thermal_signature!=heat_->config_signature()||!t.d_->field->matches(actual)||t.d_->content!=trial_content(*t.d_))return false;for(const auto&c:t.d_->sources){candidate(c,t.d_->qualification);if(c.d_->field!=t.d_->field)return false;}for(const auto&h:t.d_->handoffs)if(!h.d_||h.d_->content!=handoff_content(*h.d_))return false;return true;}catch(...){return false;}}
inline bool Owner::commit(Trial&&t,const FinalFieldView&actual)noexcept{if(!preflight(t,actual))return false;publish_prevalidated(std::move(t));return true;}
inline void Owner::publish_prevalidated(Trial&&t)noexcept{if(!t.d_||!live(t.d_->life,t.d_->generation))std::terminate();static_assert(std::is_nothrow_swappable_v<decltype(inventory_)>&&std::is_nothrow_swappable_v<decltype(layers_)>);inventory_.swap(t.d_->inventory);layers_.swap(t.d_->layers);neutral_cumulative_.swap(t.d_->neutral_cumulative);++generation_;t.d_.reset();}
}

#include "component_checkpoint_v7.hpp"
