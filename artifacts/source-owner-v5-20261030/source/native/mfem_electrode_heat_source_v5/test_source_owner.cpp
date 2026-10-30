// Focused tests of heat/source owner v5. The first part is the v4 focused suite
// (97 checks) ported unchanged except for the namespace; the v5 checks follow.
#include <cstdlib>
#include <new>
thread_local bool forbid_allocation=false;
void*operator new(std::size_t n){if(forbid_allocation)throw std::bad_alloc();if(void*p=std::malloc(n?n:1))return p;throw std::bad_alloc();}
void operator delete(void*p)noexcept{std::free(p);}
void operator delete(void*p,std::size_t)noexcept{std::free(p);}
void*operator new[](std::size_t n){return ::operator new(n);}
void operator delete[](void*p)noexcept{::operator delete(p);}
void operator delete[](void*p,std::size_t)noexcept{::operator delete(p);}
#include "hall-owner/development/hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "source_owner.hpp"
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

namespace s5=electrode_heat_v5;
namespace {
using Owner=s5::Owner;
eh::Patch water(double n=1e15,double low=-.01,double rate=1e-31){auto p=patch();p.desorption_enabled=true;p.state_id=p.parent_id="H2O";p.initial_particles=n;p.initial_low_particles=low;p.mass_kg=18*1.66053906660e-27;p.attempt_s=rate;p.activation_eV=p.desorption_eV=p.internal_eV=0;return p;}
s5::SourceSettings settings(){return {se::Config{4.3},policy(),"component-authored-marks-no-production-arrival-policy"};}
std::unique_ptr<Owner> owner(eh::Patch p=patch(),double T=1500,std::string layer="conditional-single-layer"){return std::make_unique<Owner>(geom(),material(),std::vector<eh::Patch>{p},std::vector<s5::PatchIdentity>{{"upper-foil-lower-surface",layer,.005*.025}},std::vector<s5::SourceSettings>{settings()},T,"v3-source-scenario",5e-9);}
TrialFrame frame4(const Owner&o,std::uint64_t stage=1){TrialFrame f;f.role=FieldRole::TotalMidpoint;f.accepted_revision=o.thermal_state().revision;f.trial_stage=stage;f.mesh_revision=1;f.material_revision=2;f.begin_s=o.thermal_state().time;f.dt_s=1e-14;f.end_s=f.begin_s+f.dt_s;f.field_stage_s=f.begin_s+.5*f.dt_s;f.arena_hash=o.state_sha256();f.field_hash="declared-parent-field";f.mesh_hash="declared-parent-mesh";f.material_hash="global-material-v2-map";f.species_hash="actual-source-v1-electron";return f;}
PatchBinding bind4(const Owner&o){PatchBinding b;b.thermal_revision=o.thermal_state().revision;b.source_accepted_state=b.thermal_revision+1;b.patch_id=7;b.thermal_mesh_identity=geom().revision;b.thermal_material_hash=eh::material_headers_sha256;b.material_phase_admitted=b.planar_barrier_admitted=b.electron_supply_admitted=true;b.temperature=HeldTemperature::AcceptedPatchMean;b.extrusion_center_z_m=0;b.extrusion_height_m=geom().height;return b;}
std::shared_ptr<const s5::Field>field4(const axial::Mesh&g,const Owner&o,std::uint64_t stage,double F){auto a=arrays(g,F,0);return std::make_shared<const s5::Field>(g,frame4(o,stage),a.E,a.B);}
s5::FinalFieldView view(const s5::Field&f,const axial::Mesh&g,const Arrays&a){return {f.frame(),g,a.E,a.B};}
s5::PreparedPatch prep(Owner&o,const std::shared_ptr<const s5::Field>&f,double area=1){return o.prepare_patch(o.make_epoch(1e-14,"actual-parent-midpoint"),0,bind4(o),f,origin,{f->frame().begin_s,1e-14,area,"physical-quadrature"});}
BirthSample sample4(const s5::PreparedPatch&p,const Owner&o,double w=1e9){BirthSample b;b.physical_weight=w;b.normal_uniform=.6;b.transverse_uniform=.35;b.azimuth_uniform=.1;b.time_uniform=.23;b.event_id=123;b.canonical_particle_id=42;b.count_owner_revision=o.thermal_state().revision;b.count_owner_base_hash="external-accepted-count-base";b.count_owner_token="external-staged-ticket-needs-parent-preflight";b.exposure_key=p.exposure_key();b.time_policy=BirthTime::UniformHeldExposure;b.weight_provenance=WeightProvenance::StagedCanonicalOwner;return b;}
std::vector<std::vector<s5::NeutralLeaf>>leaves(const std::vector<s5::ReleaseQuote>&q){std::vector<std::vector<s5::NeutralLeaf>>out;for(const auto&r:q){out.emplace_back();for(double w:r.released.positive_carriers())out.back().push_back({w,eh::Draw{.5,.5,.5,.5,.5,.5}});}return out;}
void dump(std::ofstream&f,std::string label,const s5::Count&n){f<<label;for(double w:n.positive_carriers())f<<','<<std::hexfloat<<w;f<<std::defaultfloat<<'\n';}

void legacy_parity(const axial::Mesh&g){
 auto wp=water(1e15,0,1e10);wp.desorption_eV=.87;wp.internal_eV=.2;auto o=owner(wp);eh::Coupler old(geom(),material(),{wp},1500,"v3-source-scenario",5e-9);auto f=field4(g,*o,1,0);auto a=arrays(g,0,0);auto v=view(*f,g,a);auto e=o->make_epoch(1e-14,"actual-parent-midpoint");auto l=leaves(o->release_quotes(e));mfem::ConstantCoefficient zero(0);auto corrected=o->prepare(e,zero,f,v,{}, {},l);auto prior=old.prepare(old.make_epoch(1e-14,"actual-parent-midpoint"),zero,{},{{l[0][0].draw}});const auto&x=corrected.packets()[0];const auto&y=prior.packets()[0];check(x.weight==y.weight&&x.T==y.T&&x.time==y.time&&x.x==y.x&&x.v==y.v&&x.proper_velocity_m_s==y.proper_velocity_m_s&&x.kinetic_J==y.kinetic_J&&x.internal_J==y.internal_J&&x.adsorbed_reference_change_J==y.adsorbed_reference_change_J,"admissible neutral source law/phase-space changed");check(corrected.heat_budget().thermal_delta_J==prior.budget().thermal_delta_J,"admissible neutral heat solve changed");
 auto metal=owner();auto extracting=field4(g,*metal,1,1e8);auto p=prep(*metal,extracting);auto birth=metal->begin_source(p,sample4(p,*metal),event(1));eh::Coupler reference(geom(),material(),{patch()},1500,"v3-source-scenario",5e-9);auto data=arrays(g,1e8,0);auto rawfield=std::make_shared<const AxialTrialField>(std::make_shared<const axial::Mesh>(g),extracting->frame(),data.E,data.B);auto rp=prepare_patch(reference,reference.make_epoch(1e-14,"actual-parent-midpoint"),0,binding(reference),rawfield,origin,{rawfield->frame().begin_s,1e-14,1,"physical-quadrature"},se::Config{4.3});auto marks=sample4(p,*metal);marks.exposure_key=rp.exposure_key;auto raw=begin_source(rp,policy(),marks,event(1));check(bytes(raw.state)==bytes(birth.state())&&raw.exchange.wall_heat_J==birth.exchange().wall_heat_J&&raw.exchange.field_work_J==birth.exchange().field_work_J,"admissible electron kernel/phase-space changed");
 auto other_arrays=arrays(g,1.1e8,0);auto changed=std::make_shared<const s5::Field>(g,extracting->frame(),other_arrays.E,other_arrays.B);auto pp=prep(*metal,changed);check(p.exposure_key()!=pp.exposure_key(),"actual field bits absent from immutable exposure key");auto stale=sample4(p,*metal);rejects([&]{metal->begin_source(pp,stale,event(1));},"old count receipt admitted with different actual field bits");
 auto ids=std::vector<s5::PatchIdentity>{{"upper-foil-lower-surface","conditional-single-layer",.000124}};rejects([&]{Owner invalid(geom(),material(),{patch()},ids,{settings()},1500,"v3-source-scenario",5e-9);},"wrong physical patch area admitted");auto relaxed=settings();relaxed.layer.maximum_quantum_magnetic_energy_fraction=.3;rejects([&]{Owner invalid(geom(),material(),{patch()},{{"surface","layer",.000125}},{relaxed},1500,"v3-source-scenario",5e-9);},"quantum applicability guard relaxed");
}
void exact_stock(const axial::Mesh&g,const std::filesystem::path&out){auto o=owner(water());mfem::ConstantCoefficient zero(0);std::ofstream oracle(out/"count-oracle.csv");dump(oracle,"initial",o->inventory(0));s5::Count gas;
 for(int step=0;step<12;++step){auto f=field4(g,*o,step+1,0);auto a=arrays(g,0,0);auto v=view(*f,g,a);auto epoch=o->make_epoch(1e-14,"actual-parent-midpoint");auto l=leaves(o->release_quotes(epoch));auto before=o->state_sha256();auto trial=o->prepare(epoch,zero,f,v,{}, {},l);check(o->state_sha256()==before,"neutral preparation mutated accepted state");for(const auto&p:trial.packets())gas=gas.credited(p.weight);auto expected=o->inventory(0).debit(l[0][0].physical_weight);check(trial.inventory(0)==expected&&trial.inventory(0)!=o->inventory(0),"C1 positive exact debit lost");if(step==0){check(o->preflight(trial,v),"common barrier preflight failed");forbid_allocation=true;o->publish_prevalidated(std::move(trial));forbid_allocation=false;}else check(o->commit(std::move(trial),v),"exact source common component commit failed");dump(oracle,"source-"+std::to_string(step+1),o->inventory(0));dump(oracle,"gas-"+std::to_string(step+1),gas);}
 bool two_word_refused=false;try{(void)o->inventory(0).require_two_limb_values();}catch(...){two_word_refused=true;}check(two_word_refused,"multi-scale inventory silently exported through two words");o->save((out/"exact-stock.bin").string());auto restored=owner(water());restored->load((out/"exact-stock.bin").string());check(restored->state_sha256()==o->state_sha256(),"exact stock component restart differs");
}

void configuration_binding(){auto a=material(),b=material();a.id="a b";a.phase="c";b.id="a";b.phase="b c";eh::Coupler olda(geom(),a,{patch()},1500,"v3-source-scenario",5e-9),oldb(geom(),b,{patch()},1500,"v3-source-scenario",5e-9);check(olda.config_signature()==oldb.config_signature(),"configuration counterexample does not collide in legacy text");Owner x(geom(),a,{patch()},{{"surface","layer",.000125}},{settings()},1500,"v3-source-scenario",5e-9),y(geom(),b,{patch()},{{"surface","layer",.000125}},{settings()},1500,"v3-source-scenario",5e-9);check(x.configuration_sha256()!=y.configuration_sha256(),"authored configuration strings lack typed length binding");}
void fresh_rate_view(const axial::Mesh&g){auto o=owner(water(1e15,-.01,1e5));mfem::ConstantCoefficient zero(0);auto f=field4(g,*o,1,0);auto a=arrays(g,0,0);auto v=view(*f,g,a);auto e=o->make_epoch(1e-14,"actual-parent-midpoint");auto first=o->release_quotes(e);auto t=o->prepare(e,zero,f,v,{}, {},leaves(first));check(o->commit(std::move(t),v),"fresh-stock fixture commit failed");auto next=o->release_quotes(o->make_epoch(1e-14,"actual-parent-midpoint"));double fraction=-std::expm1(-std::exp(std::log(1e5)+std::log(1e-14)));double expected=o->inventory(0).rounded_for_rate()*fraction;check(next[0].released==s5::Count(expected)&&next[0].released!=first[0].released,"next rate used a stale or independent stock");check(o->thermal_state().inventory[0].hi==0&&o->thermal_state().inventory[0].lo==0,"heat engine carries a second source stock");}
void exact_allocation(const axial::Mesh&g){auto o=owner(water(1,0,1e300));auto f=field4(g,*o,1,0);auto a=arrays(g,0,0);auto v=view(*f,g,a);auto e=o->make_epoch(1e-14,"actual-parent-midpoint");auto base=o->state_sha256();auto draw=eh::Draw{.5,.5,.5,.5,.5,.5};std::vector<std::vector<s5::NeutralLeaf>>bad(1);double assigned=0;for(int i=0;i<10;++i){double w=i==9?1-assigned:.1;assigned+=w;bad[0].push_back({w,draw});}check(assigned==1,"old counterexample no longer exercises rounded check");mfem::ConstantCoefficient zero(0);rejects([&]{o->prepare(e,zero,f,v,{}, {},bad);},"C2 excess exact leaf sum accepted");check(o->state_sha256()==base,"C2 refusal mutated source or heat");std::vector<std::vector<s5::NeutralLeaf>>good(1,std::vector<s5::NeutralLeaf>(8,{.125,draw}));auto t=o->prepare(e,zero,f,v,{}, {},good);s5::Count total;for(const auto&p:t.packets()){check(p.weight==.125,"explicit leaf weight changed");total=total.credited(p.weight);}check(total==s5::Count(1)&&t.inventory(0).zero(),"exact newborn sum/debit do not match");bad[0].assign(4097,{1./8192,draw});rejects([&]{o->prepare(e,zero,f,v,{}, {},bad);},"proposal cap truncated instead of rejecting");check(o->state_sha256()==base,"cap rejection changed accepted state");}
void immutable_and_underflow(axial::Mesh&g){auto o=owner(patch(),300);auto f=field4(g,*o,1,0);auto a=arrays(g,0,0);auto v=view(*f,g,a);auto p=prep(*o,f);auto base=o->state_sha256();
 static_assert(std::is_const_v<std::remove_reference_t<decltype(std::declval<const s5::Candidate&>().state())>>);static_assert(!std::is_constructible_v<s5::Candidate,Candidate>);static_assert(!std::is_constructible_v<s5::PreparedPatch,PreparedPatch>);
 for(double w:{1e-310,1e-305}){bool refused=false;try{o->begin_source(p,sample4(p,*o,w),event(1));}catch(const s5::Refusal&e){refused=std::isfinite(e.log_abs_scope)&&(w!=1e-310||e.quantity=="particle charge C");}check(refused,"C4 weighted charge/energy was zeroed or refused without scope");}
 bool log_refused=false;try{prep(*o,f,1e-250);}catch(const s5::Refusal&e){log_refused=e.quantity=="electron expected exposure count"&&std::abs(e.log_abs_scope+714.51649451502885)<1e-10;}check(log_refused,"finite positive log exposure silently admitted as scalar zero");
 auto birth=o->begin_source(p,sample4(p,*o,p.expected_number()),event(1));check(birth.exchange().particle_charge_change_C!=0&&birth.exchange().wall_heat_J!=0,"representable cold fractional source rejected or zeroed");check(o->preflight(birth,f,v),"valid immutable source preflight failed");auto again=o->begin_source(p,sample4(p,*o,p.expected_number()),event(1));check(again.content_sha256()==birth.content_sha256(),"identical nonlinear retry changed source");
 auto affine=arrays(g,1e8,0,true);auto frozen=std::make_shared<const s5::Field>(g,frame4(*o,9),affine.E,affine.B);auto original=frozen->surface(origin,{0,-1,0});double dx=g.dx;g.dx*=1.001;check(frozen->surface(origin,{0,-1,0}).E==original.E,"affine field borrowed a mutable mesh");check(!f->matches(v)&&!o->preflight(birth,f,v),"mutated final mesh accepted");g.dx=dx;auto changed=a;changed.E[0]=1;auto changed_view=view(*f,g,changed);check(!o->preflight(birth,f,changed_view),"mutated final field accepted");
 s5::Particle canonical;canonical.id=42;canonical.mass_kg=se::me;canonical.charge_C=-se::qe;canonical.weight=birth.state().native.particle.weight;auto hand=o->handoff(birth,v,canonical);check(hand.remaining_duration_s()>0&&hand.boundary_surface()<0,"actual thermal remainder failed");mfem::ConstantCoefficient zero(0);std::vector<s5::Candidate>sources{birth};std::vector<s5::Handoff>hands{hand};auto trial=o->prepare(o->make_epoch(1e-14,"actual-parent-midpoint"),zero,f,v,sources,hands,{{}});check(!o->commit(std::move(trial),changed_view)&&o->state_sha256()==base,"changed field published source/heat");check(o->commit(std::move(trial),v),"valid source/heat publication failed");check(!o->preflight(birth,f,v),"stale source survived accepted commit");}
void owner_aba(const axial::Mesh&g){alignas(Owner)std::byte storage[sizeof(Owner)];auto make=[&](){return new(storage)Owner(geom(),material(),{patch()},{{"upper-foil-lower-surface","conditional-single-layer",.000125}},{settings()},1500,"v3-source-scenario",5e-9);};auto*o=make();auto f=field4(g,*o,1,1e8);auto a=arrays(g,1e8,0);auto v=view(*f,g,a);auto p=prep(*o,f);auto b=o->begin_source(p,sample4(p,*o),event(1));std::destroy_at(o);auto*n=make();check(!n->preflight(b,f,v),"C3 owner address ABA accepted old source");rejects([&]{n->begin_source(p,sample4(p,*n),event(1));},"destroyed-owner patch reused");std::destroy_at(n);}
void active_restart(const axial::Mesh&g,const std::filesystem::path&out){auto o=owner(water());auto f=field4(g,*o,1,1e8);auto a=arrays(g,1e8,0);auto v=view(*f,g,a);auto p=prep(*o,f);auto b=o->begin_source(p,sample4(p,*o),event(1));auto c=o->advance_birth(p,b,event(2));auto again=o->advance_birth(p,b,event(2));check(c.heat_events().size()==2&&c.heat_events()[0].event_time_s==b.state().native.current_time_s&&c.heat_events()[1].event_time_s==c.state().native.current_time_s,"source birth/continuation heat times lost");check(c.content_sha256()==again.content_sha256(),"newborn same-base retry differs");check(c.state().native.status==se::LayerStatus::Active&&c.state().native.current_time_s==f->frame().end_s,"Active birth was not advanced to common endpoint");std::vector<s5::Candidate>sources{c};auto receipt=o->import(f,v,sources,{42});check(std::abs(receipt.work_defect_J)<=1e-10*std::abs(receipt.field_work_J),"actual source current/work does not close");mfem::ConstantCoefficient zero(0);auto epoch=o->make_epoch(1e-14,"actual-parent-midpoint");auto l=leaves(o->release_quotes(epoch));auto t=o->prepare(epoch,zero,f,v,sources,{},l);check(o->commit(std::move(t),v)&&o->active_layers().size()==1,"Active layer acceptance failed");o->save((out/"active-source.bin").string());auto restored=owner(water());restored->load((out/"active-source.bin").string());check(o->state_sha256()==restored->state_sha256(),"accepted layer/count/heat restart differs");
 auto fn=field4(g,*o,2,1.1e8),fr=field4(g,*restored,2,1.1e8);auto an=arrays(g,1.1e8,0);auto vn=view(*fn,g,an),vr=view(*fr,g,an);auto pn=prep(*o,fn),pr=prep(*restored,fr);auto omitted_epoch=o->make_epoch(1e-14,"actual-parent-midpoint");auto before_omission=o->state_sha256();rejects([&]{o->prepare(omitted_epoch,zero,fn,vn,{}, {},leaves(o->release_quotes(omitted_epoch)));},"accepted Active layer silently omitted");check(o->state_sha256()==before_omission,"omitted layer refusal changed accepted owners");auto wrong_event=event(3);wrong_event.energy_owner="other-canonical-energy-owner";rejects([&]{o->advance_active(pn,42,wrong_event);},"accepted canonical energy-owner identity changed");auto cn=o->advance_active(pn,42,event(3)),cr=restored->advance_active(pr,42,event(3));check(cn.content_sha256()==cr.content_sha256(),"exact next Active proposal differs after restart");std::vector<s5::Candidate>sn{cn},sr{cr};auto in=o->import(fn,vn,sn,{42}),ir=restored->import(fr,vr,sr,{42});check(in.current_load_C_m_s==ir.current_load_C_m_s&&in.endpoint_charge_impulse_C==ir.endpoint_charge_impulse_C,"restarted source currents differ");auto en=o->make_epoch(1e-14,"actual-parent-midpoint"),er=restored->make_epoch(1e-14,"actual-parent-midpoint");auto tn=o->prepare(en,zero,fn,vn,sn,{},leaves(o->release_quotes(en)));auto tr=restored->prepare(er,zero,fr,vr,sr,{},leaves(restored->release_quotes(er)));check(tn.content_sha256()==tr.content_sha256(),"next full component proposal differs");check(o->commit(std::move(tn),vn)&&restored->commit(std::move(tr),vr)&&o->state_sha256()==restored->state_sha256(),"next component accepted state differs");
 auto data=s5::detail::read_file((out/"active-source.bin").string(),64000000);data[data.size()/2]^=1;std::ofstream corrupt(out/"corrupt.bin",std::ios::binary);corrupt.write(data.data(),data.size());corrupt.close();auto fresh=owner(water());auto initial=fresh->state_sha256();rejects([&]{fresh->load((out/"corrupt.bin").string());},"corrupt component restart accepted");check(fresh->state_sha256()==initial,"failed component restart changed accepted owner");auto wrong=owner(water(),1500,"different-material-layer");rejects([&]{wrong->load((out/"active-source.bin").string());},"changed layer identity checkpoint accepted");}
}

// ============================================================================
// v5 checks: A08 G1 (weighted current/charge arrays at handoff and import),
// G2 (finite-log exposure scalar), subnormal weighted values on births and
// neutral leaves, typed scoped refusals for stock drain and sealed guards.
// Every admitted/refused case is also written to precision-dump.json for the
// independent exact-arithmetic oracle (precision_oracle.py).
// ============================================================================
namespace {
std::string hx(double x){std::ostringstream s;s<<std::hexfloat<<x;return s.str();}
std::string jq(const std::string&x){std::string r="\"";for(char c:x){if(c=='"'||c=='\\'){r+='\\';r+=c;}else if(c=='\n')r+="\\n";else if(static_cast<unsigned char>(c)<0x20){char b[8];std::snprintf(b,sizeof b,"\\u%04x",unsigned(c));r+=b;}else r+=c;}return r+"\"";}
struct J {
 std::vector<std::string>items;
 J&s(const std::string&k,const std::string&v){items.push_back(jq(k)+":"+jq(v));return *this;}
 J&b(const std::string&k,bool v){items.push_back(jq(k)+":"+(v?"true":"false"));return *this;}
 J&i(const std::string&k,long long v){items.push_back(jq(k)+":"+std::to_string(v));return *this;}
 J&h(const std::string&k,double v){items.push_back(jq(k)+":"+jq(hx(v)));return *this;}
 J&raw(const std::string&k,const std::string&v){items.push_back(jq(k)+":"+v);return *this;}
 std::string str()const{std::string r="{";for(std::size_t n=0;n<items.size();++n){if(n)r+=",";r+=items[n];}return r+"}";}
};
std::string arr(const std::vector<std::string>&v){std::string r="[";for(std::size_t n=0;n<v.size();++n){if(n)r+=",";r+=v[n];}return r+"]";}
std::string refusal_json(const s5::Refusal&r){J j;j.s("quantity",r.quantity).s("kind",s5::kind_name(r.kind)).h("log_abs_scope",r.log_abs_scope).h("relative_error",r.relative_error).h("tolerance",r.tolerance).s("identity",r.identity).s("detail",r.detail).h("observed",r.observed).h("reference",r.reference).i("binary_exponent",r.binary_exponent);return j.str();}
bool typed_scoped(const s5::Refusal&r){return !r.quantity.empty()&&!r.identity.empty()&&std::isfinite(r.log_abs_scope);}
constexpr long double T=s5::relative_tolerance;

// One weighted value against its unit-weight counterpart: |a-w*u| <= T*w*g (g=|u| for scalars).
struct Tally {long values=0,zeroed=0;long double worst=0;std::string worst_name;};
std::string value_json(const std::string&name,double u,double a,double g){return J().s("name",name).h("unit",u).h("actual",a).h("gross",g).str();}
void tally(Tally&t,const std::string&name,double u,double a,double w,double g){
 ++t.values;if(u!=0&&a==0)++t.zeroed;
 long double err;if(g==0)err=a==0?0:INFINITY;else err=std::fabs((long double)a-(long double)u*w)/((long double)w*g);
 if(!(err<=t.worst)){t.worst=err;t.worst_name=name;}}
struct Gross{foil::Vec J,Q;};
Gross path_gross(const axial::Mesh&g,const std::vector<PrecisionPath>&paths,double dt){Gross r{foil::Vec(g.ne),foil::Vec(g.nn)};
 for(const auto&p:paths){axial::Path path(g,p.absolute_from(),p.displacement_m);const double q=std::abs(p.charge_C);for(const auto&w:path.electric)r.J[w.index]+=q*std::abs(p.displacement_m[w.axis])*w.value/dt;axial::point_deposit(g,r.Q,p.absolute_from(),q);axial::point_deposit(g,r.Q,p.absolute_to(),q);}return r;}
Gross handoff_gross(const axial::Mesh&g,const s5::Candidate&c,const s5::Handoff&h,double dt){Gross r{foil::Vec(g.ne),foil::Vec(g.nn)};if(!(h.remaining_duration_s()>0))return r;
 V3 start=c.state().surface_origin_m;for(int d=0;d<3;++d)start[d]+=c.state().native.position_from_surface_m[d];const auto D=h.displacement_m();axial::Path path(g,start,D);
 for(const auto&w:path.electric)r.J[w.index]+=se::qe*std::abs(D[w.axis])*w.value/dt;axial::point_deposit(g,r.Q,start,se::qe);axial::point_deposit(g,r.Q,h.particle().position_m,se::qe);return r;}
std::string array_json(const std::string&name,const foil::Vec&u,const foil::Vec&a,const foil::Vec&gross,double w,Tally&t){std::vector<std::string>e;
 for(std::size_t i=0;i<u.size();++i)if(u[i]!=0||a[i]!=0||gross[i]!=0){e.push_back("["+std::to_string(i)+","+jq(hx(u[i]))+","+jq(hx(a[i]))+","+jq(hx(gross[i]))+"]");tally(t,name,u[i],a[i],w,gross[i]);}
 return J().s("name",name).raw("entries",arr(e)).str();}
s5::Particle electron(double w){s5::Particle c;c.id=42;c.mass_kg=se::me;c.charge_C=-se::qe;c.weight=w;return c;}
std::vector<std::pair<std::string,double>>exchange_named(const se::LayerExchange&x){return {{"particle count",x.physical_electron_number_change},{"particle charge C",x.particle_charge_change_C},{"electrode charge C",x.electrode_charge_change_C},{"wall heat J",x.wall_heat_J},{"particle kinetic increment J",x.particle_K_change_J},{"electron reference J",x.electron_reference_change_J},{"image energy J",x.image_energy_change_J},{"field work J",x.field_work_J}};}
std::vector<std::pair<std::string,double>>ledger_named(const se::Ledger&l){return {{"ledger particle count",l.physical_electron_number},{"ledger particle charge C",l.particle_charge_C},{"ledger electrode charge C",l.electrode_charge_C},{"ledger wall heat J",l.wall_heat_J},{"ledger outgoing kinetic J",l.outgoing_K_J},{"ledger incoming kinetic J",l.incoming_K_J},{"ledger electron reference J",l.electron_reference_change_J},{"ledger image reference J",l.image_reference_change_J},{"ledger matching field work J",l.resolved_matching_field_work_J}};}

struct SweepStats {long admitted=0,refused=0,untyped=0,admitted_subnormal_value=0;Tally t;std::vector<std::string>refused_quantities;double smallest_admitted=INFINITY,largest_refused=0;};
// Births (and continuation, import, handoff) at each weight against the unit-weight chain.
std::string birth_sweep(const axial::Mesh&g,const std::string&name,double Tk,double F,const std::vector<double>&weights,SweepStats&st){
 auto o=owner(patch(),Tk);auto f=field4(g,*o,1,F);auto a=arrays(g,F,0);auto v=view(*f,g,a);auto p=prep(*o,f);const double dt=f->frame().dt_s;const auto base=o->state_sha256();
 auto unit=o->begin_source(p,sample4(p,*o,1),event(1));const bool active=unit.state().native.status==se::LayerStatus::Active;
 auto unit_final=active?o->advance_birth(p,unit,event(2)):unit;
 const bool handed=unit_final.state().native.status==se::LayerStatus::PICHandedOff;
 std::optional<s5::Handoff>unit_hand;if(handed)unit_hand=o->handoff(unit_final,v,electron(1));
 std::optional<foil::SurfacePathReceipt>unit_import;if(!unit_final.paths().empty())unit_import=o->import(f,v,std::span<const s5::Candidate>(&unit_final,1),{42});
 const auto import_gross=path_gross(g,unit_final.paths(),dt);std::optional<Gross>hand_gross;if(unit_hand)hand_gross=handoff_gross(g,unit_final,*unit_hand,dt);
 std::vector<std::string>rows;
 for(double w:weights){J row;row.h("weight",w);
  try{
   auto c=o->begin_source(p,sample4(p,*o,w),event(1));
   auto fin=active?o->advance_birth(p,c,event(2)):c;
   std::vector<std::string>values,arrays_;Tally local;
   auto scalar=[&](const std::string&n,double u,double x){values.push_back(value_json(n,u,x,std::abs(u)));tally(local,n,u,x,w,std::abs(u));if(x!=0&&std::abs(x)<DBL_MIN)++st.admitted_subnormal_value;};
   for(auto [n,u]:exchange_named(unit.exchange())){double x=0;for(auto [m,y]:exchange_named(c.exchange()))if(m==n)x=y;scalar("birth "+n,u,x);}
   for(auto [n,u]:ledger_named(unit.state().native.particle.ledger)){double x=0;for(auto [m,y]:ledger_named(c.state().native.particle.ledger))if(m==n)x=y;scalar("birth "+n,u,x);}
   check(fin.paths().size()==unit_final.paths().size(),"v5 weighted path count differs from unit weight");
   for(std::size_t k=0;k<fin.paths().size();++k){scalar("path charge C",unit_final.paths()[k].charge_C,fin.paths()[k].charge_C);scalar("path work J",unit_final.paths()[k].field_work_J,fin.paths()[k].field_work_J);}
   if(active)for(auto [n,u]:exchange_named(unit_final.exchange())){double x=0;for(auto [m,y]:exchange_named(fin.exchange()))if(m==n)x=y;scalar("combined "+n,u,x);}
   if(unit_import){auto r=o->import(f,v,std::span<const s5::Candidate>(&fin,1),{42});arrays_.push_back(array_json("import current load C m/s",unit_import->current_load_C_m_s,r.current_load_C_m_s,import_gross.J,w,local));arrays_.push_back(array_json("import endpoint charge impulse C",unit_import->endpoint_charge_impulse_C,r.endpoint_charge_impulse_C,import_gross.Q,w,local));scalar("import field work J",unit_import->field_work_J,r.field_work_J);row.h("import_continuity_C",r.continuity_C).h("unit_import_continuity_C",unit_import->continuity_C);}
   if(handed){auto h=o->handoff(fin,v,electron(w));arrays_.push_back(array_json("handoff current load C m/s",unit_hand->current().current_load_C_m_s,h.current().current_load_C_m_s,hand_gross->J,w,local));arrays_.push_back(array_json("handoff endpoint charge impulse C",unit_hand->current().endpoint_charge_impulse_C,h.current().endpoint_charge_impulse_C,hand_gross->Q,w,local));
    scalar("handoff work J",unit_hand->current().field_work_J,h.current().field_work_J);scalar("handoff kinetic increment J",unit_hand->kinetic_change_J(),h.kinetic_change_J());const double K=raikou::mfem_backend::pic_kinetic_energy_J(h.particle());scalar("handoff endpoint kinetic J",K,K*w);
    row.h("handoff_continuity_C",h.current().continuity_C).h("unit_handoff_continuity_C",unit_hand->current().continuity_C);
    check(h.current().continuity_C<=1e-12*se::qe*w,"v5 admitted handoff continuity above its guard");}
   ++st.admitted;st.smallest_admitted=std::min(st.smallest_admitted,w);
   st.t.values+=local.values;st.t.zeroed+=local.zeroed;if(!(local.worst<=st.t.worst)){st.t.worst=local.worst;st.t.worst_name=local.worst_name;}
   check(local.zeroed==0,"v5 admitted weighted value silently zeroed at weight "+hx(w));
   check(local.worst<=T,"v5 admitted weighted value outside tolerance at weight "+hx(w)+": "+local.worst_name);
   row.b("admitted",true).raw("values",arr(values)).raw("arrays",arr(arrays_)).h("worst_relative_error",double(local.worst)).s("worst_quantity",local.worst_name);
  }catch(const s5::Refusal&r){++st.refused;st.largest_refused=std::max(st.largest_refused,w);st.refused_quantities.push_back(r.quantity);row.b("admitted",false).raw("refusal",refusal_json(r));check(typed_scoped(r),"v5 refusal without quantity/scope/identity");check(o->state_sha256()==base,"v5 refused birth changed the accepted owner");}
  catch(const std::exception&e){++st.untyped;row.b("admitted",false).s("untyped_error",e.what());}
  rows.push_back(row.str());}
 check(st.untyped==0,"v5 tiny-weight sweep produced an untyped error ("+name+")");
 return J().s("fixture",name).h("T_K",Tk).h("F_V_m",F).b("active_birth",active).b("handed_off",handed).raw("rows",arr(rows)).i("admitted",st.admitted).i("refused",st.refused).i("untyped",st.untyped).i("admitted_values",st.t.values).i("admitted_values_zeroed",st.t.zeroed).i("admitted_subnormal_values",st.admitted_subnormal_value).h("worst_admitted_relative_error",double(st.t.worst)).s("worst_admitted_quantity",st.t.worst_name).h("smallest_admitted_weight",st.smallest_admitted).h("largest_refused_weight",st.largest_refused).str();
}

std::string exposure_sweep(const axial::Mesh&g,long&admitted,long&refused){
 auto o=owner(patch(),300);auto f=field4(g,*o,1,0);const auto base=o->state_sha256();auto one=prep(*o,f,1);
 const V3 n{0,-1,0};auto spectrum=se::evaluate(300,se::extraction_field(f->surface(origin,n).E,n),se::Config{4.3});
 check(spectrum.log_current_A_m2==one.log_current_A_m2()&&spectrum.mean_metal_energy_eV==one.mean_metal_energy_eV()&&one.temperature_K()==300,"v5 exposure oracle inputs differ from the prepared spectrum");
 std::vector<double>areas{1,1e-50,1e-100,1e-200,1e-220,1e-228,1e-230};for(int k=0;k<=60;++k)areas.push_back(std::pow(10.,-225.-.5*k));
 for(double x:{3.1622776601683792e-245,1e-250,1e-260,1e-280,1e-300,1e-305,1e-310,1e-320})areas.push_back(x);
 std::vector<std::string>rows;
 for(double x:areas){J row;row.h("area_fraction",x);
  try{auto p=prep(*o,f,x);++admitted;row.b("admitted",true).h("count",p.expected_number()).h("log_count",p.log_expected_number()).h("diagnostic_J",p.mean_wall_heat_diagnostic_J());
   const long double exact=std::exp((long double)spectrum.log_current_A_m2)*(long double)(.005*.025)*x*(long double)1e-14/(long double)se::qe;check(std::fabs((p.expected_number()-exact)/exact)<=T,"v5 admitted exposure scalar outside tolerance");}
  catch(const s5::Refusal&r){++refused;row.b("admitted",false).raw("refusal",refusal_json(r));check(typed_scoped(r)&&r.kind==s5::RefusalKind::exposure_scalar&&o->state_sha256()==base,"v5 exposure refusal untyped/unscoped or changed owner");}
  rows.push_back(row.str());}
 return J().h("T_K",300).h("area_m2",.005*.025).h("duration_s",1e-14).h("qe",se::qe).h("log_current_A_m2",spectrum.log_current_A_m2).h("mean_metal_energy_eV",spectrum.mean_metal_energy_eV).raw("rows",arr(rows)).str();
}

// Neutral leaves: an exact partition of a unit release reaching subnormal leaf weights.
std::string neutral_leaves(const axial::Mesh&g){
 auto o=owner(water(1,0,1e300));auto f=field4(g,*o,1,0);auto a=arrays(g,0,0);auto v=view(*f,g,a);auto e=o->make_epoch(1e-14,"actual-parent-midpoint");mfem::ConstantCoefficient zero(0);const eh::Draw draw{.5,.5,.5,.5,.5,.5};const auto base=o->state_sha256();
 auto unit_trial=o->prepare(e,zero,f,v,{},{},{{{1.,draw}}});const auto unit=unit_trial.packets().at(0);
 auto row=[&](const std::vector<double>&w,const std::string&label){J j;j.s("partition",label);std::vector<std::string>ws;s5::Count total;for(double x:w){ws.push_back(jq(hx(x)));total=total.credited(x);}j.raw("leaves",arr(ws)).b("exact_partition_of_one",total==s5::Count(1.));
  std::vector<s5::NeutralLeaf>l;for(double x:w)l.push_back({x,draw});
  try{auto t=o->prepare(e,zero,f,v,{},{},{l});std::vector<std::string>vals;Tally tl;
   for(const auto&p:t.packets()){const double x=p.weight;auto add=[&](const std::string&n,double u,double y,double gross){vals.push_back(J().h("weight",x).s("name",n).h("unit",u).h("actual",y).h("gross",gross).str());tally(tl,n,u,y,x,gross);};
    add("neutral kinetic J",unit.kinetic_J,p.kinetic_J,std::abs(unit.kinetic_J));add("neutral classical sample kinetic J",unit.classical_sample_K_J,p.classical_sample_K_J,std::abs(unit.classical_sample_K_J));
    add("neutral relativistic conversion J",unit.relativistic_conversion_J,p.relativistic_conversion_J,std::abs(unit.kinetic_J)+std::abs(unit.classical_sample_K_J));
    for(int k=0;k<3;++k)add("neutral wall support impulse N s",unit.wall_support_impulse_N_s[k],p.wall_support_impulse_N_s[k],std::abs(unit.wall_support_impulse_N_s[k]));}
   check(tl.zeroed==0&&tl.worst<=T,"v5 admitted neutral leaf value zeroed or outside tolerance");
   j.b("admitted",true).raw("values",arr(vals)).h("worst_relative_error",double(tl.worst)).s("worst_quantity",tl.worst_name);}
  catch(const s5::Refusal&r){j.b("admitted",false).raw("refusal",refusal_json(r));check(typed_scoped(r)&&o->state_sha256()==base,"v5 neutral refusal untyped/unscoped or changed owner");}
  return j.str();};
 std::vector<double>deep{1-0x1p-53};for(int k=1;k<20;++k)deep.push_back(std::ldexp(1.,-53*k)-std::ldexp(1.,-53*(k+1)));deep.push_back(std::ldexp(1.,-1060));
 std::vector<double>normal{.5,.25,.25-0x1p-54,0x1p-54};
 std::vector<double>edge{1-0x1p-53};for(int k=1;k<18;++k)edge.push_back(std::ldexp(1.,-53*k)-std::ldexp(1.,-53*(k+1)));edge.push_back(std::ldexp(1.,-53*18));
 auto r1=row(deep,"chain to 2^-1060"),r2=row(normal,"normal-range leaves"),r3=row(edge,"chain to 2^-954");
 check(r1.find("\"admitted\":false")!=std::string::npos,"v5 subnormal neutral leaf admitted");check(r2.find("\"admitted\":true")!=std::string::npos,"v5 normal-range exact partition refused");
 return J().h("unit_kinetic_J",unit.kinetic_J).raw("cases",arr({r1,r2,r3})).str();
}

std::string stock_drain(const axial::Mesh&g){
 auto o=owner(water(1e15,-.01,1e300));auto f=field4(g,*o,1,0);auto a=arrays(g,0,0);auto v=view(*f,g,a);auto e=o->make_epoch(1e-14,"actual-parent-midpoint");mfem::ConstantCoefficient zero(0);const auto base=o->state_sha256();
 std::optional<s5::Refusal>quote,prep_;try{(void)o->release_quotes(e);}catch(const s5::Refusal&r){quote.emplace(r);}
 try{(void)o->prepare(e,zero,f,v,{},{},{{{1e15,eh::Draw{.5,.5,.5,.5,.5,.5}}}});}catch(const s5::Refusal&r){prep_.emplace(r);}
 check(quote&&quote->kind==s5::RefusalKind::stock_drain&&quote->quantity=="neutral release count"&&typed_scoped(*quote)&&quote->identity.find("exact-stock-sha256=")!=std::string::npos,"v5 stock drain not a typed scoped refusal");
 check(prep_&&prep_->kind==s5::RefusalKind::stock_drain&&o->state_sha256()==base,"v5 stock drain refusal in prepare missing or changed owner");
 std::vector<std::string>carriers;for(double c:o->inventory(0).positive_carriers())carriers.push_back(jq(hx(c)));
 // Positive control: a full release of an exactly representable stock is covered.
 auto exact=owner(water(1,0,1e300));auto q=exact->release_quotes(exact->make_epoch(1e-14,"actual-parent-midpoint"));check(q[0].released==s5::Count(1.),"v5 exact full release refused");
 return J().raw("stock_carriers",arr(carriers)).raw("quote_refusal",quote?refusal_json(*quote):"null").raw("prepare_refusal",prep_?refusal_json(*prep_):"null").b("owner_unchanged",o->state_sha256()==base).str();
}

std::string guard_conversion(){
 const s5::detail::Weighted W(1e-300,"unit-test;canonical-particle=42");std::vector<std::string>rows;
 auto convert=[&](const std::string&m){try{W.at_weight([&]()->int{throw std::invalid_argument(m);},1e-20);}catch(const s5::Refusal&r){return std::optional<s5::Refusal>(r);}catch(...){}return std::optional<s5::Refusal>();};
 for(std::string m:{"electrode heat: source incremental aggregate energy failed","Micro-layer aggregate energy closure failed","electrode heat: actual handoff remainder work/continuity failed W=1"}){auto r=convert(m);check(r&&r->kind==s5::RefusalKind::sealed_guard_at_weight&&r->quantity!="weighted sealed numerical guard"&&std::abs(r->log_abs_scope-(std::log(1e-20)+std::log(1e-300)))<1e-12&&r->detail.find(m)!=std::string::npos,"v5 sealed energy-closure guard not converted to a typed refusal");rows.push_back(refusal_json(*r));}
 bool passthrough=false;try{W.at_weight([&]()->int{throw s5::Refusal("x",1,"id",s5::RefusalKind::weighted_representation);},1);}catch(const s5::Refusal&r){passthrough=r.quantity=="x";}check(passthrough,"v5 typed refusal altered by guard conversion");
 bool alloc=false;try{W.at_weight([&]()->int{throw std::bad_alloc();},1);}catch(const std::bad_alloc&){alloc=true;}catch(...){}check(alloc,"v5 allocation failure reclassified as a numerical refusal");
 return arr(rows);
}

// Exact binary rescaling is satisfied for normal-range weights: no new refusal.
std::string normal_range(const axial::Mesh&g){
 std::vector<std::string>rows;
 for(double Tk:{300.,1500.}){const double F=Tk==300?0:1e8;auto o=owner(patch(),Tk);auto f=field4(g,*o,1,F);auto a=arrays(g,F,0);auto v=view(*f,g,a);auto p=prep(*o,f);
  for(double w:{1e-250,1e-200,1e-100,1e-61,0x1p-60,.75,1.,1.5,2.,3.,1e9,1e100,1e250}){bool ok=true;std::string why;
   try{auto c=o->begin_source(p,sample4(p,*o,w),event(1));auto fin=c.state().native.status==se::LayerStatus::Active?o->advance_birth(p,c,event(2)):c;if(fin.state().native.status==se::LayerStatus::PICHandedOff)(void)o->handoff(fin,v,electron(w));if(!fin.paths().empty())(void)o->import(f,v,std::span<const s5::Candidate>(&fin,1),{42});}
   catch(const std::exception&e){ok=false;why=e.what();}
   check(ok,"v5 refused a normal-range weight "+hx(w)+": "+why);rows.push_back(J().h("T_K",Tk).h("weight",w).b("admitted",ok).str());}}
 return arr(rows);
}
}

int main(int argc,char**argv){try{
 check(argc==2,"output directory required");std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
 cpu_set_t cpus;CPU_ZERO(&cpus);sched_getaffinity(0,sizeof(cpus),&cpus);const int cpu=CPU_ISSET(6,&cpus)?6:7;
 check(CPU_COUNT(&cpus)==1&&(CPU_ISSET(6,&cpus)||CPU_ISSET(7,&cpus))&&getpriority(PRIO_PROCESS,0)==15,"CPU6-or-7/nice15 resource contract");
 axial::Mesh g(45,4,.009,1);
 // ---- the 97 v4 focused checks, ported unchanged except for the namespace
 configuration_binding();legacy_parity(g);exact_stock(g,out);fresh_rate_view(g);exact_allocation(g);immutable_and_underflow(g);owner_aba(g);active_restart(g,out);
 const int ported=checks;check(ported==97,"ported v4 check count changed");
 // ---- v5 checks
 auto g1_owner=owner(patch(),300);auto g1f=field4(g,*g1_owner,1,0);auto g1p=prep(*g1_owner,g1f);const auto g1base=g1_owner->state_sha256();
 std::optional<s5::Refusal>g1;try{(void)g1_owner->begin_source(g1p,sample4(g1p,*g1_owner,1e-301),event(1));}catch(const s5::Refusal&r){g1.emplace(r);}
 check(g1&&typed_scoped(*g1)&&g1->identity.find("canonical-particle=42")!=std::string::npos&&g1_owner->state_sha256()==g1base,"v5 G1 weight 1e-301 birth not refused with typed scope");
 std::optional<s5::Refusal>g2;try{(void)prep(*g1_owner,g1f,3.1622776601683792e-245);}catch(const s5::Refusal&r){g2.emplace(r);}
 check(g2&&g2->kind==s5::RefusalKind::exposure_scalar&&g2->quantity=="electron expected exposure count"&&g2->relative_error>.9&&std::isfinite(g2->log_abs_scope),"v5 G2 imprecise exposure scalar not refused with typed scope");
 std::vector<double>cold{1e300,1e250,1e100,1e9,2,1.5,1,.75,1e-9,1e-61,1e-100,1e-200,1e-250,1e-301,1e-310};for(int k=0;k<=60;++k)cold.push_back(std::pow(10.,-260.-k));
 for(int k=0;k<=40;++k)cold.push_back(std::pow(10.,-270.-.25*k));
 std::vector<double>hot{1e250,1e100,1e9,1,1e-100,1e-200};for(int k=0;k<=35;++k)hot.push_back(std::pow(10.,-250.-2*k));for(int k=0;k<=40;++k)hot.push_back(std::pow(10.,-270.-.5*k));
 SweepStats cs,hs;auto cold_json=birth_sweep(g,"300K-F0-cold-handoff",300,0,cold,cs);auto hot_json=birth_sweep(g,"1500K-F1e8-active-continuation",1500,1e8,hot,hs);
 check(cs.admitted>0&&cs.refused>0&&hs.admitted>0&&hs.refused>0,"v5 sweeps do not cover both admitted and refused weights");
 long ea=0,er=0;auto exposure_json=exposure_sweep(g,ea,er);check(ea>0&&er>0,"v5 exposure sweep does not cover both outcomes");
 auto neutral_json=neutral_leaves(g);auto drain_json=stock_drain(g);auto guard_json=guard_conversion();auto normal_json=normal_range(g);
 {std::ofstream d(out/"precision-dump.json");d<<J().s("schema","source-v5-precision-dump-v1").h("relative_tolerance",s5::relative_tolerance).h("evaluation_bound",s5::evaluation_bound)
   .raw("g1_refusal",refusal_json(*g1)).raw("g2_refusal",refusal_json(*g2)).raw("birth_sweeps",arr({cold_json,hot_json})).raw("exposure_sweep",exposure_json)
   .raw("neutral_leaves",neutral_json).raw("stock_drain",drain_json).raw("sealed_guard_conversion",guard_json).raw("normal_range_weights",normal_json).str()<<"\n";}
 const int added=checks-ported;
 std::ofstream result(out/"results.json");result<<"{\"status\":\"passed\",\"checks\":"<<checks<<",\"ported_v4_checks\":"<<ported<<",\"new_v5_checks\":"<<added<<",\"nice\":15,\"full_joint_test\":false,\"arrival_policy_implemented\":false}\n";
 std::cout<<checks<<" focused source v5 checks passed ("<<ported<<" ported v4, "<<added<<" new v5) on CPU"<<cpu<<"\n";
 std::cout<<"cold sweep admitted="<<cs.admitted<<" refused="<<cs.refused<<" smallest_admitted="<<cs.smallest_admitted<<" largest_refused="<<cs.largest_refused<<" worst="<<double(cs.t.worst)<<"\n";
 std::cout<<"hot sweep admitted="<<hs.admitted<<" refused="<<hs.refused<<" smallest_admitted="<<hs.smallest_admitted<<" largest_refused="<<hs.largest_refused<<" worst="<<double(hs.t.worst)<<"\n";
 std::cout<<"exposure sweep admitted="<<ea<<" refused="<<er<<"\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
