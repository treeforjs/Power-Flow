// Port of the pre-integration audit reproducer
//   Power-Flow/artifacts/pre-integration-audit-20261030/reviews/source/source_v3_reproducer.cpp
// to heat/source owner v5 (runtimes/pic-electrode-heat-source-v5-20261030), derived
// from the source-v4 port (artifacts/audit-reproducer-ports-20261030/source-v4/port-source).
// Same rows, inputs and recorded observations as the source-v4 port; the A08 extension
// sweeps are widened so admitted and refused weights/areas both appear, and the exact
// operands of every admitted/refused A08 case are written for the Fraction oracle.
//
// Every audit input (inventories, rates, temperatures, fields, weights, areas,
// draws, event identities, frames and bindings) is reused unchanged. Only the
// API calls differ. Outcomes are RECORDED rather than asserted, so a defective
// result is reported instead of aborting the run. Probes whose id starts with
// "EXT-" are additional probes, not audit counterexamples.
//
// usage: source_v5_port <output-directory>
#include "hall-owner/development/hall-finite-foil-em-20261029/axial_surface_path_port.hpp"
#include "source_owner.hpp"
#include "ss304_heat_material_v2.hpp"
#include <cfloat>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sched.h>
#include <sys/resource.h>
using namespace electrode_heat_v3;
namespace s5=electrode_heat_v5;

namespace {
// ---------------------------------------------------------------- recording
std::string hexs(double x){std::ostringstream s;s<<std::hexfloat<<x;return s.str();}
std::string decs(double x){std::ostringstream s;s<<std::setprecision(17)<<x;return s.str();}
std::string quote(const std::string&x){std::string r="\"";for(char c:x){if(c=='"'||c=='\\'){r+='\\';r+=c;}else if(c=='\n')r+="\\n";else if(static_cast<unsigned char>(c)<0x20){char b[8];std::snprintf(b,sizeof b,"\\u%04x",unsigned(c));r+=b;}else r+=c;}return r+"\"";}
struct Obj {
 std::vector<std::string>items;
 Obj&s(const std::string&k,const std::string&v){items.push_back(quote(k)+":"+quote(v));return *this;}
 Obj&b(const std::string&k,bool v){items.push_back(quote(k)+":"+(v?"true":"false"));return *this;}
 Obj&i(const std::string&k,long long v){items.push_back(quote(k)+":"+std::to_string(v));return *this;}
 Obj&d(const std::string&k,double v){items.push_back(quote(k)+":"+quote(decs(v)));items.push_back(quote(k+"_hex")+":"+quote(hexs(v)));return *this;}
 Obj&raw(const std::string&k,const std::string&v){items.push_back(quote(k)+":"+v);return *this;}
 std::string str()const{std::string r="{";for(std::size_t n=0;n<items.size();++n){if(n)r+=",";r+=items[n];}return r+"}";}
};
std::string hexlist(const std::vector<double>&v){std::string r="[";for(std::size_t n=0;n<v.size();++n){if(n)r+=",";r+=quote(hexs(v[n]));}return r+"]";}
std::string jarray(const std::vector<std::string>&v){std::string r="[";for(std::size_t n=0;n<v.size();++n){if(n)r+=",";r+=v[n];}return r+"]";}
struct Outcome {bool threw=false,typed=false;std::string what,quantity,kind,identity;double log_abs_scope=NAN,relative_error=NAN,observed=NAN,reference=NAN;int exponent=0;};
template<class F>Outcome attempt(F&&f){Outcome o;try{f();}catch(const s5::Refusal&e){o.threw=o.typed=true;o.what=e.what();o.quantity=e.quantity;o.log_abs_scope=e.log_abs_scope;o.kind=s5::kind_name(e.kind);o.identity=e.identity;o.relative_error=e.relative_error;o.observed=e.observed;o.reference=e.reference;o.exponent=e.binary_exponent;}catch(const std::exception&e){o.threw=true;o.what=e.what();}return o;}
std::string record(const Outcome&o){Obj x;x.b("refused",o.threw).b("typed_refusal_with_scope",o.typed).s("message",o.what);if(o.typed)x.s("quantity",o.quantity).d("log_abs_scope",o.log_abs_scope).s("kind",o.kind).s("identity",o.identity).d("relative_error",o.relative_error).d("observed",o.observed).d("reference",o.reference).i("binary_exponent",o.exponent);return x.str();}
bool scoped(const Outcome&o){return o.typed&&std::isfinite(o.log_abs_scope)&&!o.identity.empty()&&!o.quantity.empty();}
constexpr long double T5=s5::relative_tolerance;
Obj exact_dump; // hexadecimal exact-count operands for the independent Fraction oracle
std::vector<std::string>summary;

// ------------------------------------------- audit fixture (unchanged inputs)
eh::Geometry geom(){eh::Geometry g;g.nx=4;g.ny=2;g.revision="thermal-mesh-v3";g.electrode="upper";return g;}
eh::Material material(){return eh::conditional_ss304_v2(ss304v2::Qualification::source_conditioned_liquid_3000,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000,eh::material_headers_sha256);}
eh::Patch patch(){eh::Patch p;p.id="lower-face";p.face=eh::Face::LowY;p.source_hash="v3-source-scenario";return p;}
eh::SourceEventIdentity event(std::uint64_t seq){return {"source-v3","e-","electrode-electron-reservoir","native-layer","canonical-electron-42",seq,seq};}
V3 origin{.0001,.0025,.0002};
se::LayerPolicy policy(){se::LayerPolicy p;p.local_linear_vacuum_extent_m=25e-6;return p;}
eh::Patch water_patch(double high,double low,double rate){auto p=patch();p.desorption_enabled=true;p.state_id="H2O";p.parent_id="H2O";p.initial_particles=high;p.initial_low_particles=low;p.mass_kg=18*1.66053906660e-27;p.attempt_s=rate;p.activation_eV=p.desorption_eV=p.internal_eV=0;return p;}
struct Arrays{foil::Vec E,B;};
Arrays arrays(const axial::Mesh&g,double F,double B,bool affine=false){Arrays a{foil::Vec(g.ne),foil::Vec(g.nb)};std::vector<V3> nodes(g.nn);for(int k=0;k<=g.nz;++k)for(int j=0;j<=g.ny;++j)for(int i=0;i<=g.nx;++i)nodes[g.node(i,j,k)]={i*g.dx-g.side/2,j*g.dy-g.side/2,g.zmin+k*g.dz};for(int e=0;e<g.ne;++e){const auto&edge=g.edges[e];V3 x{};for(int d=0;d<3;++d)x[d]=.5*(nodes[edge.a][d]+nodes[edge.b][d]);if(edge.axis==0)a.E[e]=affine?1e6+1e13*(x[1]-origin[1]):0;if(edge.axis==1)a.E[e]=F+(affine?1e13*(x[0]-origin[0]):0);if(edge.axis==2)a.E[e]=affine?2e6+1e13*(x[1]-origin[1]):0;}
 for(int k=0;k<g.nz;++k)for(int j=0;j<g.ny;++j)for(int i=0;i<=g.nx;++i)a.B[g.bx(i,j,k)]=B+(affine?1e6*(i*g.dx-g.side/2-origin[0]):0);
 for(int k=0;k<g.nz;++k)for(int j=0;j<=g.ny;++j)for(int i=0;i<g.nx;++i)a.B[g.by(i,j,k)]=affine?-1e6*(j*g.dy-g.side/2-origin[1]):0;
 return a;}

// ------------------------------------------------- v5 API adaptation only (as for v4)
// v5 (as v4) requires a physical surface/layer identity, frozen source settings and an
// explicit external sampling-policy identity per patch. These are identity
// fields; the physical area is the actual geometry patch area.
s5::SourceSettings settings(se::Config c=se::Config{4.3}){return {c,policy(),"audit-port-caller-authored-marks-no-arrival-policy"};}
std::unique_ptr<s5::Owner> owner(eh::Patch p,double T,se::Config c=se::Config{4.3}){return std::make_unique<s5::Owner>(geom(),material(),std::vector<eh::Patch>{p},std::vector<s5::PatchIdentity>{{"upper-foil-lower-surface","conditional-single-layer",.005*.025}},std::vector<s5::SourceSettings>{settings(c)},T,"v3-source-scenario",5e-9);}
// Same frame/binding values as the v3 reproducer; accepted revision/time come from the v5 owner.
TrialFrame frame(const s5::Owner&h,std::uint64_t stage=1,double dt=1e-14){TrialFrame f;f.role=FieldRole::TotalMidpoint;f.accepted_revision=h.thermal_state().revision;f.trial_stage=stage;f.mesh_revision=1;f.material_revision=2;f.begin_s=h.thermal_state().time;f.end_s=f.begin_s+dt;f.dt_s=dt;f.field_stage_s=f.begin_s+.5*dt;f.arena_hash="accepted-arena-revision-"+std::to_string(h.thermal_state().revision);f.field_hash="immutable-midpoint-trial-"+std::to_string(stage);f.mesh_hash="Hall-axial-sealed-mesh";f.material_hash="global-material-v2-map";f.species_hash="actual-source-v1-electron";return f;}
PatchBinding binding(const s5::Owner&h){PatchBinding b;b.thermal_revision=h.thermal_state().revision;b.source_accepted_state=h.thermal_state().revision+1;b.patch_id=7;b.thermal_mesh_identity=geom().revision;b.thermal_material_hash=eh::material_headers_sha256;b.material_phase_admitted=b.planar_barrier_admitted=b.electron_supply_admitted=true;b.temperature=HeldTemperature::AcceptedPatchMean;b.extrusion_center_z_m=0;b.extrusion_height_m=geom().height;return b;}
// The parent's actual final arrays are retained so FinalFieldView can reference them.
struct ParentField{std::shared_ptr<const s5::Field>f;Arrays a;};
ParentField field(const axial::Mesh&g,const s5::Owner&h,std::uint64_t stage,double F,double B=0,bool affine=false){auto a=arrays(g,F,B,affine);auto f=std::make_shared<const s5::Field>(g,frame(h,stage),a.E,a.B);return {f,std::move(a)};}
s5::FinalFieldView view(const ParentField&p,const axial::Mesh&g){return {p.f->frame(),g,p.a.E,p.a.B};}
s5::PreparedPatch prepared(s5::Owner&h,const std::shared_ptr<const s5::Field>&f,double area=1,const char*quadrature="native-surface-quadrature"){return h.prepare_patch(h.make_epoch(f->frame().dt_s,"trial-midpoint-field"),0,binding(h),f,origin,Exposure{f->frame().begin_s,f->frame().dt_s,area,quadrature});}
BirthSample sample(const s5::PreparedPatch&p,const s5::Owner&h,double weight=1e9,bool qualification=true){BirthSample b;b.physical_weight=weight;b.normal_uniform=.6;b.transverse_uniform=.35;b.azimuth_uniform=.1;b.time_uniform=.23;b.event_id=123;b.canonical_particle_id=42;b.count_owner_revision=h.thermal_state().revision;b.count_owner_base_hash="unchanged-accepted-source-count-owner";b.count_owner_token=qualification?"qualification-weight-not-an-event-rate":"staged-canonical-fractional-count-receipt";b.exposure_key=p.exposure_key();b.time_policy=BirthTime::ExposureBeginLumped;b.weight_provenance=qualification?WeightProvenance::QualificationOnly:WeightProvenance::StagedCanonicalOwner;return b;}
s5::Particle canonical_electron(double weight){s5::Particle c;c.id=42;c.mass_kg=se::me;c.charge_C=-se::qe;c.weight=weight;return c;}
std::vector<double>exchange_values(const se::LayerExchange&x){return {x.physical_electron_number_change,x.particle_charge_change_C,x.electrode_charge_change_C,x.wall_heat_J,x.particle_K_change_J,x.electron_reference_change_J,x.image_energy_change_J,x.field_work_J};}
const char* exchange_names[]={"particle count","particle charge C","electrode charge C","wall heat J","particle kinetic increment J","electron reference J","image energy J","field work J"};
s5::Count exact_sum(const std::vector<double>&w){s5::Count c;for(double x:w)c=c.credited(x);return c;}
s5::Count exact_plus(s5::Count a,const s5::Count&b){for(double w:b.positive_carriers())a=a.credited(w);return a;}
void finish(Obj&j,const std::string&id,const std::string&verdict){j.s("v5_verdict",verdict);summary.push_back(id+" -> "+verdict);std::cout<<id<<" v5_verdict="<<verdict<<"\n";}

// ======================================================= A07 / source C1
std::string c1_release(const axial::Mesh&g){
 Obj j;const std::string id="A07-C1-actual-release";j.s("id",id).s("audit_finding","A07 (source review C1)");
 j.s("audit_input","heat-v2 Coupler water patch inventory {hi=1e15, lo=-0.01}, attempt 1e-31 1/s, activation/desorption/internal 0 eV, T=1500 K, dt=1e-14 s, epoch label manufactured-source-audit, one draw {0.5 x6}");
 j.s("v3_behavior","Coupler::prepare emitted packet weight 1.0000000000000024e-30 while both candidate inventory limbs stayed {1e15,-0.01}");
 const double v3_weight=std::strtod("1.0000000000000024e-30",nullptr);j.d("v3_packet_weight",v3_weight);
 auto o=owner(water_patch(1e15,-.01,1e-31),1500);auto pf=field(g,*o,1,0);auto v=view(pf,g);mfem::ConstantCoefficient zero(0);const eh::Draw draw{.5,.5,.5,.5,.5,.5};
 const auto initial=o->inventory(0);const auto base=o->state_sha256();
 j.b("initial_stock_equals_count2_input_exactly",initial==s5::Count(1e15,-.01));
 j.b("thermal_engine_legacy_inventory_zero",o->thermal_state().inventory.at(0).hi==0&&o->thermal_state().inventory.at(0).lo==0);
 auto e=o->make_epoch(1e-14,"manufactured-source-audit");
 std::vector<s5::ReleaseQuote>q;auto qo=attempt([&]{q=o->release_quotes(e);});j.raw("quote_outcome",record(qo));
 if(q.empty()){finish(j,id,"harness-error-no-quote");return j.str();}
 j.raw("quote_released_carriers",hexlist(q[0].released.positive_carriers())).d("quote_log_released",q[0].log_released).d("quote_log_hazard",q[0].log_hazard);
 const bool quote_is_v3=q[0].released==s5::Count(v3_weight);j.b("quote_equals_v3_emitted_weight_exactly",quote_is_v3);
 // A caller leaf that misses the exact debit by one ulp must refuse.
 const double below=std::nextafter(v3_weight,0.);
 auto mo=attempt([&]{(void)o->prepare(e,zero,pf.f,v,{},{},{{{below,draw}}});});j.raw("one_ulp_low_leaf",record(mo)).b("one_ulp_low_leaf_owner_unchanged",o->state_sha256()==base);
 // The same single leaf the v3 owner returned.
 std::optional<s5::Trial>t;auto po=attempt([&]{t.emplace(o->prepare(e,zero,pf.f,v,{},{},{{{v3_weight,draw}}}));});j.raw("v3_leaf_prepare",record(po)).b("prepare_left_accepted_owner_unchanged",o->state_sha256()==base);
 exact_dump.raw("c1_initial",hexlist(initial.positive_carriers())).raw("c1_count2_input",hexlist({1e15,-.01})).raw("c1_leaves",hexlist({v3_weight}));
 if(!t){finish(j,id,mo.threw?"refused":"still-defective");return j.str();}
 const auto cand=t->inventory(0);std::vector<double>pw;for(const auto&p:t->packets())pw.push_back(p.weight);const auto gas=exact_sum(pw);
 j.i("packet_count",long(pw.size())).raw("returned_packet_weights",hexlist(pw)).b("returned_weight_verbatim",pw.size()==1&&pw[0]==v3_weight);
 j.b("candidate_stock_changed",cand!=initial).b("candidate_stock_equals_exact_debit",cand==initial.debit(v3_weight)).b("candidate_plus_returned_leaves_equals_initial",exact_plus(cand,gas)==initial);
 j.raw("candidate_stock_carriers",hexlist(cand.positive_carriers()));
 const bool committed=o->commit(std::move(*t),v);j.b("committed",committed);
 const auto accepted=o->inventory(0);
 j.b("accepted_stock_equals_candidate",accepted==cand).b("accepted_stock_plus_gas_credit_equals_initial",exact_plus(accepted,gas)==initial).b("thermal_engine_legacy_inventory_zero_after_commit",o->thermal_state().inventory.at(0).hi==0&&o->thermal_state().inventory.at(0).lo==0);
 auto two=attempt([&]{(void)accepted.require_two_limb_values();});j.raw("two_limb_export_of_accepted_stock",record(two)).i("accepted_stock_carrier_count",long(accepted.positive_carriers().size()));
 exact_dump.raw("c1_accepted",hexlist(accepted.positive_carriers())).raw("c1_returned",hexlist(pw));
 const bool exact=quote_is_v3&&mo.threw&&pw.size()==1&&pw[0]==v3_weight&&cand!=initial&&exact_plus(cand,gas)==initial&&committed&&accepted==cand&&exact_plus(accepted,gas)==initial;
 finish(j,id,exact?"exact":"still-defective");return j.str();
}
std::string c1_arithmetic(){
 Obj j;const std::string id="A07-C1-count-arithmetic";j.s("id",id).s("audit_finding","A07 (source review C1, reproducer line source_count_debit_lost)");
 j.s("audit_input","eh::add(eh::Count{1e15,-0.01}, -1e-30)").s("v3_behavior","two-word debit returned hi=0x1.c6bf52634p+49 lo=-0x1.47ae147ae147bp-7, both limbs unchanged");
 const auto legacy=eh::add(eh::Count{1e15,-.01},-1e-30);
 j.b("legacy_v2_two_word_add_still_loses_debit",legacy.hi==1e15&&legacy.lo==-.01).s("legacy_note","sealed heat-v2 function is unchanged; v5 (as v4) keeps the source stock in ExactParticleCount and holds the heat-v2 legacy inventory at zero (see A07-C1-actual-release)");
 const s5::Count c(1e15,-.01);const auto d=c.debit(1e-30);
 j.b("exact_debit_changed_stock",d!=c).b("exact_debit_plus_amount_equals_initial",d.credited(1e-30)==c).raw("exact_after_carriers",hexlist(d.positive_carriers()));
 exact_dump.raw("c1_arith_after",hexlist(d.positive_carriers())).raw("c1_arith_amount",hexlist({1e-30}));
 finish(j,id,d!=c&&d.credited(1e-30)==c?"exact (legacy two-word function unchanged but not the v5 stock)":"still-defective");return j.str();
}

// ======================================================= A07 / source C2
std::string c2_allocation(const axial::Mesh&g){
 Obj j;const std::string id="A07-C2-ten-packet-allocation";j.s("id",id).s("audit_finding","A07 (source review C2)");
 j.s("audit_input","heat-v2 Coupler water patch inventory {1,0}, attempt 1e300 1/s, T=1500 K, dt=1e-14 s, ten draws {0.5 x6}");
 j.s("v3_behavior","full release of 1 into ten packets passed assigned==released; exact packet sum 1+5/2^55");
 // The ten weights exactly as recorded in reviews/source/exact-allocation-evidence.json.
 const char* v3_hex[10]={"0x1.999999999999ap-4","0x1.999999999999ap-4","0x1.999999999999ap-4","0x1.999999999999ap-4","0x1.999999999999ap-4","0x1.999999999999ap-4","0x1.999999999999ap-4","0x1.999999999999ap-4","0x1.999999999999ap-4","0x1.99999999999a0p-4"};
 std::vector<double>v3w;for(auto*h:v3_hex)v3w.push_back(std::strtod(h,nullptr));
 auto o=owner(water_patch(1,0,1e300),1500);auto pf=field(g,*o,1,0);auto v=view(pf,g);mfem::ConstantCoefficient zero(0);const eh::Draw draw{.5,.5,.5,.5,.5,.5};
 auto e=o->make_epoch(1e-14,"manufactured-source-audit");const auto base=o->state_sha256();
 std::vector<s5::ReleaseQuote>q;auto qo=attempt([&]{q=o->release_quotes(e);});j.raw("quote_outcome",record(qo));
 if(q.empty()){finish(j,id,"harness-error-no-quote");return j.str();}
 j.raw("quote_released_carriers",hexlist(q[0].released.positive_carriers())).b("quote_equals_one",q[0].released==s5::Count(1.));
 // The v3 allocation rule (equal split, last = released - assigned) reproduces the audit weights.
 const double released=q[0].released.rounded_for_rate();double assigned=0;std::vector<double>rule;for(int k=0;k<10;++k){double w=k+1==10?released-assigned:released/10;assigned+=w;rule.push_back(w);}
 j.b("v3_rule_reproduces_audit_weights",rule==v3w).b("v3_rule_floating_assigned_equals_released",assigned==released);
 const auto exact_v3=exact_sum(v3w);j.b("audit_weights_exact_sum_equals_debit",exact_v3==q[0].released);
 auto excess=attempt([&]{j.raw("audit_weights_exact_excess_carriers",hexlist(exact_v3.debit(1.).positive_carriers()));});j.raw("excess_computation",record(excess));
 auto leaves=[&](const std::vector<double>&w){std::vector<std::vector<s5::NeutralLeaf>>l(1);for(double x:w)l[0].push_back({x,draw});return l;};
 auto ro=attempt([&]{(void)o->prepare(e,zero,pf.f,v,{},{},leaves(v3w));});j.raw("audit_ten_leaves",record(ro)).b("audit_ten_leaves_owner_unchanged",o->state_sha256()==base);
 const std::vector<double>tenth(10,.1);auto eo=attempt([&]{(void)o->prepare(e,zero,pf.f,v,{},{},leaves(tenth));});j.raw("ten_equal_0p1_leaves",record(eo)).b("ten_equal_leaves_exact_sum_equals_debit",exact_sum(tenth)==q[0].released).b("ten_equal_leaves_owner_unchanged",o->state_sha256()==base);
 // An exact ten-leaf partition near 0.1 (nine binary64 0.1 = 7205759403792794*2^-56 and one 7205759403792790*2^-56) must be admitted verbatim.
 std::vector<double>exact10(9,.1);exact10.push_back(std::strtod("0x1.9999999999996p-4",nullptr));
 j.b("exact_partition_sums_to_one",exact_sum(exact10)==s5::Count(1.));
 std::optional<s5::Trial>t;auto ao=attempt([&]{t.emplace(o->prepare(e,zero,pf.f,v,{},{},leaves(exact10)));});j.raw("exact_ten_leaves",record(ao));
 exact_dump.raw("c2_audit_weights",hexlist(v3w)).raw("c2_exact_partition",hexlist(exact10));
 bool admitted_exact=false;
 if(t){std::vector<double>pw;for(const auto&p:t->packets())pw.push_back(p.weight);j.raw("exact_ten_returned_weights",hexlist(pw)).b("exact_ten_returned_verbatim",pw==exact10).b("exact_ten_returned_sum_equals_debit",exact_sum(pw)==q[0].released).b("exact_ten_candidate_stock_zero",t->inventory(0).zero());
  const bool committed=o->commit(std::move(*t),v);j.b("exact_ten_committed",committed).b("exact_ten_accepted_stock_zero",o->inventory(0).zero());exact_dump.raw("c2_exact_returned",hexlist(pw));
  admitted_exact=pw==exact10&&exact_sum(pw)==q[0].released&&committed&&o->inventory(0).zero();}
 const bool refused=ro.threw&&eo.threw&&!(exact_v3==q[0].released);
 finish(j,id,refused&&admitted_exact?"refused (audit leaves); exact (explicit exact partition)":"still-defective");return j.str();
}
std::string ext_overdraw(const axial::Mesh&g){
 Obj j;const std::string id="EXT-A07-full-release-of-multiscale-stock";j.s("id",id).s("audit_finding","A07 extension: audit inventory {1e15,-0.01} with the audit full-release rate 1e300");
 j.s("v3_behavior","not part of the audit reproducer; heat-v2 would retreat the released amount with nextafter").s("v4_behavior","refused with the untyped error 'insufficient exact inventory'");
 auto o=owner(water_patch(1e15,-.01,1e300),1500);auto pf=field(g,*o,1,0);auto v=view(pf,g);mfem::ConstantCoefficient zero(0);auto e=o->make_epoch(1e-14,"manufactured-source-audit");const auto base=o->state_sha256();
 auto qo=attempt([&]{(void)o->release_quotes(e);});j.raw("quote_outcome",record(qo)).b("owner_unchanged",o->state_sha256()==base);
 auto po=attempt([&]{(void)o->prepare(e,zero,pf.f,v,{},{},{{{1e15,eh::Draw{.5,.5,.5,.5,.5,.5}}}});});j.raw("prepare_with_rounded_full_leaf",record(po)).b("owner_unchanged_after_prepare",o->state_sha256()==base);
 exact_dump.raw("drain_stock",hexlist(o->inventory(0).positive_carriers())).s("drain_released",hexs(qo.observed));
 const bool typed=scoped(qo)&&scoped(po)&&qo.kind=="stock-drain"&&po.kind=="stock-drain"&&qo.identity.find("exact-stock-sha256=")!=std::string::npos;
 finish(j,id,typed&&o->state_sha256()==base?"refused (typed stock-drain refusal with log scope and stock identity; stock unchanged)":qo.threw&&po.threw?"refused (untyped)":"admitted");return j.str();
}
// ======================================================= A08 / source C4, L1
std::string c4_tiny_birth(const axial::Mesh&g){
 Obj j;const std::string id="A08-C4-tiny-weight-birth";j.s("id",id).s("audit_finding","A08 (source review C4)");
 j.s("audit_input","300 K owner, field stage 2 F=0, exposure area fraction 1, staged weight 1e-310, ExposureBeginLumped, event 1, canonical 42");
 j.s("v3_behavior","begin_source succeeded with count 9.9999999999999694e-311, particle charge -0, wall heat -0, preflight true");
 auto cold=owner(patch(),300);auto pc_=field(g,*cold,2,0);auto v=view(pc_,g);auto pc=prepared(*cold,pc_.f);const auto base=cold->state_sha256();
 auto unit=cold->begin_source(pc,sample(pc,*cold,1,false),event(1));j.raw("unit_weight_exchange",hexlist(exchange_values(unit.exchange())));
 std::optional<s5::Candidate>c;auto r=attempt([&]{c=cold->begin_source(pc,sample(pc,*cold,1e-310,false),event(1));});j.raw("weight_1e-310",record(r)).b("owner_unchanged",cold->state_sha256()==base);
 const double expected=std::log(se::qe)+std::log(1e-310);j.d("expected_log_abs_charge_scope",expected);
 if(c)j.raw("admitted_exchange",hexlist(exchange_values(c->exchange())));
 finish(j,id,r.typed&&std::isfinite(r.log_abs_scope)&&cold->state_sha256()==base?"refused (typed, with log scope)":r.threw?"refused (untyped)":"still-defective");return j.str();
}
std::string l1_exposure(const axial::Mesh&g){
 Obj j;const std::string id="A08-L1-exposure-underflow";j.s("id",id).s("audit_finding","A08 (source review L1)");
 j.s("audit_input","300 K owner, field stage 2 F=0, exposure area fraction 1e-250, dt 1e-14 s, quadrature id tiny-positive-area");
 j.s("v3_behavior","scalar expected count 0, log count -714.51649451502885 (exp gives 4.8913552196320394e-311), current_underflow=false; admitted");
 auto cold=owner(patch(),300);auto pc_=field(g,*cold,2,0);const auto base=cold->state_sha256();
 auto r=attempt([&]{(void)prepared(*cold,pc_.f,1e-250,"tiny-positive-area");});j.raw("area_1e-250",record(r)).b("owner_unchanged",cold->state_sha256()==base);
 const double v3_log=std::strtod("-714.51649451502885",nullptr);j.d("v3_log_count",v3_log);
 if(r.typed)j.b("refusal_scope_bit_identical_to_v3_log",r.log_abs_scope==v3_log).d("refusal_scope_minus_v3_log",r.log_abs_scope-v3_log);
 auto full=prepared(*cold,pc_.f,1,"tiny-positive-area");j.d("area_1_expected_count",full.expected_number()).d("area_1_log_count",full.log_expected_number());
 finish(j,id,r.typed&&r.quantity=="electron expected exposure count"&&std::isfinite(r.log_abs_scope)?"refused (typed, with log scope; finite-log admission refuses the unrepresented scalar)":r.threw?"refused (untyped)":"still-defective");return j.str();
}
std::string handoff_control(const axial::Mesh&g){
 Obj j;const std::string id="CONTROL-thermal-surface-handoff";j.s("id",id).s("audit_finding","positive control from the audit reproducer (not a defect)");
 j.s("audit_input","300 K owner, field stage 2 F=0, area 1, staged weight 1, ExposureBeginLumped, canonical electron 42 weight 1");
 j.s("v3_behavior","handoff succeeded, boundary=-1, remaining duration 9.9999999998789903e-15, endpoint 5.00001e-09");
 auto cold=owner(patch(),300);auto pc_=field(g,*cold,2,0);auto v=view(pc_,g);auto pc=prepared(*cold,pc_.f);
 auto thermal=cold->begin_source(pc,sample(pc,*cold,1,false),event(1));j.i("birth_status",int(thermal.state().native.status));
 std::optional<s5::Handoff>h;auto r=attempt([&]{h=cold->handoff(thermal,v,canonical_electron(1));});j.raw("handoff",record(r));
 bool same=false;if(h){j.i("boundary",h->boundary_surface()).d("remaining_duration_s",h->remaining_duration_s()).d("handoff_time_s",h->handoff_time_s());same=h->boundary_surface()==-1&&h->remaining_duration_s()==std::strtod("9.9999999998789903e-15",nullptr);}
 j.b("matches_v3_boundary_and_duration",same);
 finish(j,id,same?"positive-control-passed (identical to v3)":"positive-control-differs");return j.str();
}

// ======================================================= A02 / source C3
std::string a02_weight(const axial::Mesh&g){
 Obj j;const std::string id="A02-C3-changed-particle-weight";j.s("id",id).s("audit_finding","A02 (source review C3: changed particle weight)");
 j.s("audit_input","1500 K owner, field stage 1 F=1e8, area 1, staged weight 1e9; candidate particle weight then doubled");
 j.s("v3_behavior","preflight stayed true and import accepted the doubled weight with the original path charge -1.6021766339999998e-10");
 j.s("v5_direct_mutation","not expressible: Candidate::state() is const and the payload is private const (see compile-reject cases 1, 2, 3, 4)");
 auto o=owner(patch(),1500);auto pf=field(g,*o,1,1e8);auto v=view(pf,g);auto p=prepared(*o,pf.f);
 auto b=o->begin_source(p,sample(p,*o,1e9,false),event(1));
 j.d("candidate_weight",b.state().native.particle.weight).d("path0_charge_C",b.paths().empty()?NAN:b.paths()[0].charge_C).b("genuine_candidate_preflight",o->preflight(b,pf.f,v));
 auto io=attempt([&]{(void)o->import(pf.f,v,std::span<const s5::Candidate>(&b,1),{42});});j.raw("genuine_candidate_import",record(io));
 // Runtime analogue that IS expressible: the canonical particle presented at handoff carries a changed weight.
 auto cold=owner(patch(),300);auto cf=field(g,*cold,2,0);auto cv=view(cf,g);auto cp=prepared(*cold,cf.f);auto thermal=cold->begin_source(cp,sample(cp,*cold,1,false),event(1));const auto base=cold->state_sha256();
 auto ho=attempt([&]{(void)cold->handoff(thermal,cv,canonical_electron(2));});j.raw("handoff_with_doubled_canonical_weight",record(ho)).b("owner_unchanged",cold->state_sha256()==base);
 finish(j,id,ho.threw&&cold->state_sha256()==base?"not-expressible (compile-reject) + runtime analogue refused":"still-defective");return j.str();
}
std::string a02_spectrum(const axial::Mesh&g){
 Obj j;const std::string id="A02-C3-changed-spectrum";j.s("id",id).s("audit_finding","A02 (source review C3: substituted spectrum with old exposure key)");
 j.s("audit_input","PreparedPatch at 1500 K, F=1e8; spectrum replaced by se::evaluate(1600,1e8,Config{4.3}); the old exposure key reused");
 j.s("v3_behavior","candidate built from the altered patch kept the 1500 K exposure key and a changed particle");
 j.s("v5_direct_mutation","not expressible: PreparedPatch has no public payload and no conversion from v3 (see compile-reject cases 5, 6, 7)");
 auto o15=owner(patch(),1500),o16=owner(patch(),1600);auto f15=field(g,*o15,1,1e8),f16=field(g,*o16,1,1e8);auto p15=prepared(*o15,f15.f),p16=prepared(*o16,f16.f);
 const auto base15=o15->state_sha256(),base16=o16->state_sha256();
 j.b("exposure_keys_differ_1500_vs_1600",p15.exposure_key()!=p16.exposure_key()).d("T_1500_patch",p15.temperature_K()).d("T_1600_patch",p16.temperature_K());
 auto old_key=sample(p15,*o15,1e9,false);
 auto r1=attempt([&]{(void)o16->begin_source(p16,old_key,event(1));});j.raw("1600K_patch_with_1500K_count_receipt",record(r1));
 auto r2=attempt([&]{(void)o15->begin_source(p16,sample(p16,*o16,1e9,false),event(1));});j.raw("1600K_patch_presented_to_1500K_owner",record(r2));
 auto ctrl=attempt([&]{(void)o16->begin_source(p16,sample(p16,*o16,1e9,false),event(1));});j.raw("control_matching_1600K_receipt",record(ctrl));
 auto o43=owner(patch(),1500,se::Config{4.5});auto f43=field(g,*o43,1,1e8);auto p43=prepared(*o43,f43.f);
 j.b("exposure_keys_differ_work_function_4p3_vs_4p5",p43.exposure_key()!=p15.exposure_key());
 auto r3=attempt([&]{(void)o43->begin_source(p43,old_key,event(1));});j.raw("4p5eV_patch_with_4p3eV_count_receipt",record(r3));
 const bool unchanged=o15->state_sha256()==base15&&o16->state_sha256()==base16;j.b("owners_unchanged",unchanged);
 finish(j,id,r1.threw&&r2.threw&&r3.threw&&!ctrl.threw&&unchanged?"not-expressible (compile-reject) + runtime analogues refused":"still-defective");return j.str();
}
std::string a02_aba(const axial::Mesh&g){
 Obj j;const std::string id="A02-C3-reused-owner-address";j.s("id",id).s("audit_finding","A02 (source review C3: reused owner address)");
 j.s("audit_input","owner at 1500 K in fixed storage, field stage 3 F=1e8, staged weight 1 candidate; owner destroyed and an identical owner constructed at the same address");
 j.s("v3_behavior","replacement owner preflight accepted the predecessor candidate");
 alignas(s5::Owner) std::byte storage[sizeof(s5::Owner)];
 auto make=[&]{return new(storage)s5::Owner(geom(),material(),{patch()},{{"upper-foil-lower-surface","conditional-single-layer",.005*.025}},{settings()},1500,"v3-source-scenario",5e-9);};
 auto*old=make();auto fa=field(g,*old,3,1e8);auto va=view(fa,g);auto pa=prepared(*old,fa.f);auto ca=old->begin_source(pa,sample(pa,*old,1,false),event(1));
 const bool before=old->preflight(ca,fa.f,va);const auto old_address=reinterpret_cast<std::uintptr_t>(old);std::destroy_at(old);
 auto*replacement=make();const auto base=replacement->state_sha256();
 j.b("same_address",reinterpret_cast<std::uintptr_t>(replacement)==old_address).b("predecessor_preflight_before_destruction",before);
 const bool accepted=replacement->preflight(ca,fa.f,va);j.b("replacement_accepts_stale_candidate",accepted);
 auto r1=attempt([&]{(void)replacement->begin_source(pa,sample(pa,*replacement,1,false),event(1));});j.raw("replacement_begin_source_with_stale_patch",record(r1));
 mfem::ConstantCoefficient zero(0);
 auto r2=attempt([&]{(void)replacement->prepare(replacement->make_epoch(1e-14,"trial-midpoint-field"),zero,fa.f,va,std::span<const s5::Candidate>(&ca,1),{},{{}});});j.raw("replacement_prepare_with_stale_candidate",record(r2));
 auto r3=attempt([&]{(void)replacement->import(fa.f,va,std::span<const s5::Candidate>(&ca,1),{42});});j.raw("replacement_import_with_stale_candidate",record(r3));
 const bool unchanged=replacement->state_sha256()==base;j.b("replacement_owner_unchanged",unchanged);std::destroy_at(replacement);
 finish(j,id,!accepted&&r1.threw&&r2.threw&&r3.threw&&unchanged?"refused":"still-defective");return j.str();
}
std::string a02_mesh(){
 Obj j;const std::string id="A02-C3-mutated-mesh-old-metadata";j.s("id",id).s("audit_finding","A02 (source review C3: mutated mesh with old metadata)");
 j.s("audit_input","nonconst axial::Mesh(45,4,0.009,1), affine field stage 4 F=1e8; mesh dx *= 1.001 after the trial field is built; frame labels unchanged");
 j.s("v3_behavior","field gathering changed while matches_trial stayed true");
 auto o=owner(patch(),1500);axial::Mesh gm(45,4,.009,1);
 auto fm=field(gm,*o,4,1e8,0,true);auto fu=field(gm,*o,5,1e8);auto vm=view(fm,gm);auto vu=view(fu,gm);
 const auto before=fm.f->surface(origin,{0,-1,0});
 auto pu=prepared(*o,fu.f);auto cu=o->begin_source(pu,sample(pu,*o,1e9,false),event(1));
 auto e=o->make_epoch(1e-14,"trial-midpoint-field");mfem::ConstantCoefficient zero(0);
 std::optional<s5::Trial>t;t.emplace(o->prepare(e,zero,fm.f,vm,{},{},{{}}));
 j.b("pre_mutation_field_matches",fm.f->matches(vm)).b("pre_mutation_candidate_preflight",o->preflight(cu,fu.f,vu)).b("pre_mutation_trial_preflight",o->preflight(*t,vm));
 const auto base=o->state_sha256();const double dx=gm.dx;gm.dx*=1.001;
 const auto after=fm.f->surface(origin,{0,-1,0});
 j.b("sealed_field_gather_unchanged_by_alias_mutation",after.E==before.E&&after.B==before.B);
 const bool m1=fm.f->matches(vm),c1=o->preflight(cu,fu.f,vu),t1=o->preflight(*t,vm);
 j.b("mutated_mesh_field_matches",m1).b("mutated_mesh_candidate_preflight",c1).b("mutated_mesh_trial_preflight",t1);
 auto io=attempt([&]{(void)o->import(fu.f,vu,std::span<const s5::Candidate>(&cu,1),{42});});j.raw("mutated_mesh_import",record(io));
 auto ho=attempt([&]{(void)o->handoff(cu,vu,canonical_electron(1e9));});j.raw("mutated_mesh_handoff",record(ho));
 auto po=attempt([&]{(void)o->prepare(e,zero,fm.f,vm,{},{},{{}});});j.raw("mutated_mesh_new_prepare",record(po));
 const bool committed_mutated=o->commit(std::move(*t),vm);j.b("mutated_mesh_commit",committed_mutated).b("owner_unchanged_after_refusals",o->state_sha256()==base);
 gm.dx=dx;
 // Other mesh members are also bound (zmin, one dual volume, one edge length).
 const double zmin=gm.zmin;gm.zmin+=1e-9;const bool mz=fm.f->matches(vm);gm.zmin=zmin;
 const double ve=gm.VE[0];gm.VE[0]*=2;const bool mv=fm.f->matches(vm);gm.VE[0]=ve;
 const double len=gm.edges[0].length;gm.edges[0].length*=1+1e-15;const bool ml=fm.f->matches(vm);gm.edges[0].length=len;
 j.b("mutated_zmin_matches",mz).b("mutated_VE0_matches",mv).b("mutated_edge0_length_matches",ml);
 const bool restored=fm.f->matches(vm);const bool committed=o->commit(std::move(*t),vm);
 j.b("restored_mesh_matches",restored).b("restored_mesh_commit_positive_control",committed).b("predecessor_candidate_stale_after_commit",!o->preflight(cu,fu.f,vu));
 const bool ok=!m1&&!c1&&!t1&&io.threw&&ho.threw&&po.threw&&!committed_mutated&&!mz&&!mv&&!ml&&restored&&committed;
 finish(j,id,ok?"refused (mutation detected, owner unchanged; positive control commits)":"still-defective");return j.str();
}
std::string ext_publish_gap(){
 Obj j;const std::string id="EXT-A02-mutation-between-preflight-and-publish";j.s("id",id).s("audit_finding","A02 extension: scope boundary of the final check");
 j.s("v3_behavior","not part of the audit reproducer");
 auto o=owner(patch(),1500);axial::Mesh g(45,4,.009,1);auto pf=field(g,*o,1,1e8);auto v=view(pf,g);auto e=o->make_epoch(1e-14,"trial-midpoint-field");mfem::ConstantCoefficient zero(0);
 std::optional<s5::Trial>t;t.emplace(o->prepare(e,zero,pf.f,v,{},{},{{}}));const bool pre=o->preflight(*t,v);
 pf.a.E[0]+=1;const bool still=pf.f->matches(v);const auto base=o->state_sha256();
 o->publish_prevalidated(std::move(*t));const bool published=o->state_sha256()!=base;
 j.b("preflight_before_mutation",pre).b("field_matches_after_parent_array_mutation",still).b("publish_prevalidated_published",published);
 j.s("contract","API-v5 (unchanged from v4): preflight every owner, then publish_prevalidated in the same root critical section; publish_prevalidated performs no field recheck by design");
 finish(j,id,published&&!still?"observation: undetected after preflight (root critical-section contract, not a component check)":"observation: detected");return j.str();
}

// Entry-wise comparison of a weighted receipt array against unit-weight entries scaled exactly by w.
std::string compare_array(const foil::Vec&unitv,const foil::Vec&act,double w){
 long unit_nz=0,act_nz=0,zero=0;long double min_zeroed=INFINITY,max_zeroed=0,exact_sum=0,act_sum=0,abs_exact=0;double worst_nonzero=0;
 for(std::size_t i=0;i<unitv.size()&&i<act.size();++i){const long double ex=(long double)unitv[i]*(long double)w;exact_sum+=ex;act_sum+=act[i];abs_exact+=std::fabs(ex);if(unitv[i]!=0)++unit_nz;if(act[i]!=0)++act_nz;
  if(unitv[i]!=0&&act[i]==0){++zero;min_zeroed=std::min(min_zeroed,std::fabs(ex));max_zeroed=std::max(max_zeroed,std::fabs(ex));}
  else if(unitv[i]!=0&&ex!=0)worst_nonzero=std::max(worst_nonzero,double(std::fabs(((long double)act[i]-ex)/ex)));}
 Obj o;o.i("entries",long(std::min(unitv.size(),act.size()))).i("unit_nonzero",unit_nz).i("weighted_nonzero",act_nz).i("zeroed",zero).d("max_relative_error_of_nonzero_entries",worst_nonzero);
 if(zero)o.d("smallest_exact_magnitude_zeroed",double(min_zeroed)).d("largest_exact_magnitude_zeroed",double(max_zeroed)).d("half_smallest_subnormal",std::numeric_limits<double>::denorm_min()/2);
 o.d("exact_sum_unit_times_w",double(exact_sum)).d("weighted_sum",double(act_sum)).d("sum_of_exact_magnitudes",double(abs_exact));return o.str();}
std::string ext_handoff_current(const axial::Mesh&g){
 Obj j;const std::string id="EXT-A08-handoff-current-underflow";j.s("id",id).s("audit_finding","A08 extension (G1): weighted ordinary-handoff current of a tiny-weight birth (same 300 K, F=0 audit fixture)");
 j.s("v4_behavior","birth admitted; handoff zeroed 12 current-load entries, 8 endpoint charge impulses and the continuity value without refusal");
 const double w=std::strtod("1e-301",nullptr);j.d("weight",w);
 auto cold=owner(patch(),300);auto pf=field(g,*cold,2,0);auto v=view(pf,g);auto pc=prepared(*cold,pf.f);const auto base=cold->state_sha256();
 auto unit=cold->begin_source(pc,sample(pc,*cold,1,false),event(1));auto uh=cold->handoff(unit,v,canonical_electron(1));
 std::optional<s5::Candidate>c;auto bo=attempt([&]{c=cold->begin_source(pc,sample(pc,*cold,w,false),event(1));});j.raw("begin_source",record(bo)).b("owner_unchanged",cold->state_sha256()==base);
 exact_dump.raw("g1_refusal",Obj().s("quantity",bo.quantity).s("kind",bo.kind).s("observed",hexs(bo.observed)).s("reference",hexs(bo.reference)).i("binary_exponent",bo.exponent).s("weight",hexs(w)).str());
 if(!c){finish(j,id,scoped(bo)&&cold->state_sha256()==base?"refused at birth (typed, with log scope and identity; no handoff admitted)":"birth refused (untyped)");return j.str();}
 j.b("birth_preflight",cold->preflight(*c,pf.f,v));
 std::optional<s5::Handoff>h;auto ho=attempt([&]{h=cold->handoff(*c,v,canonical_electron(w));});j.raw("handoff",record(ho));
 if(!h){finish(j,id,scoped(ho)?"handoff refused (typed)":"handoff refused (untyped)");return j.str();}
 j.raw("current_load_C_m_s",compare_array(uh.current().current_load_C_m_s,h->current().current_load_C_m_s,w)).raw("endpoint_charge_impulse_C",compare_array(uh.current().endpoint_charge_impulse_C,h->current().endpoint_charge_impulse_C,w));
 long z=0;
 for(std::size_t i=0;i<uh.current().current_load_C_m_s.size();++i)if(uh.current().current_load_C_m_s[i]!=0&&h->current().current_load_C_m_s[i]==0)++z;
 for(std::size_t i=0;i<uh.current().endpoint_charge_impulse_C.size();++i)if(uh.current().endpoint_charge_impulse_C[i]!=0&&h->current().endpoint_charge_impulse_C[i]==0)++z;
 j.i("entries_zeroed",z);
 finish(j,id,z?"defect: admitted handoff with weighted current entries silently zeroed":"admitted without zeroing");return j.str();
}
std::string ext_exposure_scalar(const axial::Mesh&g){
 Obj j;const std::string id="EXT-A08-exposure-scalar-admitted-imprecise";j.s("id",id).s("audit_finding","A08 extension (G2): single area fraction from the sweep (same 300 K, F=0 audit fixture)");
 j.s("v4_behavior","admitted scalar 3.08e-305 against exact 1.547e-305 (ratio 1.994)");
 const double a=std::strtod("3.1622776601683792e-245",nullptr);j.d("area_fraction",a);
 auto cold=owner(patch(),300);auto pf=field(g,*cold,2,0);auto one=prepared(*cold,pf.f,1,"tiny-positive-area");const auto base=cold->state_sha256();
 std::optional<s5::PreparedPatch>p;auto o=attempt([&]{p=prepared(*cold,pf.f,a,"tiny-positive-area");});j.raw("prepare_patch",record(o)).b("owner_unchanged",cold->state_sha256()==base);
 exact_dump.raw("g2_refusal",Obj().s("quantity",o.quantity).s("kind",o.kind).s("observed",hexs(o.observed)).s("area_fraction",hexs(a)).str());
 if(!p){finish(j,id,scoped(o)&&o.kind=="exposure-scalar"?"refused (typed exposure-scalar refusal with finite log scope)":"refused (untyped)");return j.str();}
 const double n=p->expected_number();const long double ref=(long double)one.expected_number()*(long double)a;
 j.d("admitted_scalar_count",n).d("scalar_over_reference",double((long double)n/ref));
 finish(j,id,std::fabs((long double)n/ref-1)>T5?"gap: admitted positive scalar count differs from the exact finite-log exposure":"admitted within tolerance");return j.str();
}
// ======================================================= A08 extensions
std::string ext_weight_sweep(const axial::Mesh&g){
 Obj j;const std::string id="EXT-A08-weight-precision-sweep";j.s("id",id).s("audit_finding","A08 extension around the audit 1e-310 input (same 300 K, F=0 fixture); v4 rows 1e-300..1e-311 plus 1e-250..1e-300");
 auto cold=owner(patch(),300);auto pf=field(g,*cold,2,0);auto v=view(pf,g);auto pc=prepared(*cold,pf.f);const auto base=cold->state_sha256();
 auto unit=cold->begin_source(pc,sample(pc,*cold,1,false),event(1));const auto u=exchange_values(unit.exchange());
 std::optional<foil::SurfacePathReceipt>ui;std::optional<s5::Handoff>uh;
 j.raw("unit_import",record(attempt([&]{ui=cold->import(pf.f,v,std::span<const s5::Candidate>(&unit,1),{42});}))).raw("unit_handoff",record(attempt([&]{uh=cold->handoff(unit,v,canonical_electron(1));})));
 // Gross (sum of absolute contributions) of the unit handoff deposit, same operations as the sealed deposit.
 foil::Vec gJ(g.ne),gQ(g.nn);const double dt=pf.f->frame().dt_s;
 if(uh&&uh->remaining_duration_s()>0){V3 start=unit.state().surface_origin_m;for(int d=0;d<3;++d)start[d]+=unit.state().native.position_from_surface_m[d];const auto D=uh->displacement_m();axial::Path path(g,start,D);for(const auto&w:path.electric)gJ[w.index]+=se::qe*std::abs(D[w.axis])*w.value/dt;axial::point_deposit(g,gQ,start,se::qe);axial::point_deposit(g,gQ,uh->particle().position_m,se::qe);}
 std::vector<std::string>rows,dump;double worst_exchange=0,worst_current=0;long zeroed_total=0;long admitted=0,admitted_subnormal=0,refused=0,typed=0,owner_changed=0;std::string worst_w,worst_q;double smallest_admitted=INFINITY;
 auto entries=[&](const foil::Vec&uv,const foil::Vec&av,const foil::Vec&gv){std::vector<std::string>e;for(std::size_t i=0;i<uv.size();++i)if(uv[i]!=0||av[i]!=0||gv[i]!=0)e.push_back("["+std::to_string(i)+","+quote(hexs(uv[i]))+","+quote(hexs(av[i]))+","+quote(hexs(gv[i]))+"]");return jarray(e);};
 auto zeroed=[](const foil::Vec&unitv,const foil::Vec&act,const foil::Vec&gross,double w,double&worst){long z=0;for(std::size_t i=0;i<unitv.size()&&i<act.size();++i){if(unitv[i]!=0&&act[i]==0)++z;if(gross[i]!=0){long double rel=std::fabs(((long double)act[i]-(long double)unitv[i]*w)/((long double)gross[i]*w));worst=std::max(worst,double(rel));}else if(act[i]!=0)worst=INFINITY;}return z;};
 for(int k=-200;k<=44;++k){const double w=std::pow(10.,-300.-.25*k);Obj r;r.d("weight",w);
  std::optional<s5::Candidate>c;auto o=attempt([&]{c=cold->begin_source(pc,sample(pc,*cold,w,false),event(1));});r.raw("begin_source",record(o));
  Obj dumprow;dumprow.s("weight",hexs(w)).b("admitted",bool(c));
  if(!c){++refused;if(scoped(o))++typed;if(cold->state_sha256()!=base)++owner_changed;dumprow.s("quantity",o.quantity).s("kind",o.kind).s("observed",hexs(o.observed)).s("reference",hexs(o.reference)).i("binary_exponent",o.exponent);}
  if(c){++admitted;smallest_admitted=std::min(smallest_admitted,w);const auto a=exchange_values(c->exchange());double mx=0;std::string q;bool sub=false;
   for(int n=0;n<8;++n){if(a[n]!=0&&std::abs(a[n])<DBL_MIN)sub=true;long double ex=(long double)u[n]*(long double)w;if(ex!=0){double rel=double(std::fabs(((long double)a[n]-ex)/ex));if(rel>mx){mx=rel;q=exchange_names[n];}}else if(a[n]!=0)mx=INFINITY;}
   if(sub)++admitted_subnormal;r.d("max_relative_exchange_error",mx).s("worst_quantity",q).b("subnormal_weighted_operand",sub).raw("exchange",hexlist(a));
   dumprow.raw("unit_exchange",hexlist(u)).raw("exchange",hexlist(a));
   if(mx>worst_exchange){worst_exchange=mx;worst_w=hexs(w);worst_q=q;}
   std::optional<s5::Handoff>h;auto ho=attempt([&]{h=cold->handoff(*c,v,canonical_electron(w));});r.raw("handoff",record(ho));
   if(ho.threw)dumprow.s("handoff_refusal",ho.quantity);
   if(h&&uh){double wc=0;long z=zeroed(uh->current().current_load_C_m_s,h->current().current_load_C_m_s,gJ,w,wc)+zeroed(uh->current().endpoint_charge_impulse_C,h->current().endpoint_charge_impulse_C,gQ,w,wc);r.i("handoff_entries_zeroed",z).d("handoff_max_entry_error_relative_to_gross",wc).d("handoff_continuity_C",h->current().continuity_C).d("unit_handoff_continuity_C",uh->current().continuity_C);zeroed_total+=z;worst_current=std::max(worst_current,wc);
    dumprow.raw("handoff_current_load",entries(uh->current().current_load_C_m_s,h->current().current_load_C_m_s,gJ)).raw("handoff_endpoint_charge_impulse",entries(uh->current().endpoint_charge_impulse_C,h->current().endpoint_charge_impulse_C,gQ));}
  }
  rows.push_back(r.str());dump.push_back(dumprow.str());}
 exact_dump.raw("weight_sweep",jarray(dump));
 j.raw("unit_exchange",hexlist(u)).raw("rows",jarray(rows)).i("admitted_weights",admitted).i("refused_weights",refused).i("typed_scoped_refusals",typed).i("refusals_changing_owner",owner_changed).d("smallest_admitted_weight",smallest_admitted).i("admitted_with_subnormal_weighted_operand",admitted_subnormal).d("worst_admitted_relative_exchange_error",worst_exchange).s("worst_weight",worst_w).s("worst_quantity",worst_q).i("admitted_current_entries_zeroed",zeroed_total).d("worst_admitted_current_entry_error_relative_to_gross",worst_current);
 std::string verdict=zeroed_total?"defect: admitted proposal with silently zeroed weighted current entries":worst_exchange>T5||worst_current>T5?"gap: admitted value outside tolerance":typed!=refused||owner_changed?"refusal without typed scope or with owner change":"no admitted zeroing or value outside tolerance; every refusal typed with scope; owner unchanged";
 finish(j,id,verdict);return j.str();
}
std::string ext_exposure_sweep(const axial::Mesh&g){
 Obj j;const std::string id="EXT-A08-exposure-scalar-precision-sweep";j.s("id",id).s("audit_finding","A08 extension around the audit area 1e-250 input (same 300 K, F=0 fixture); v4 rows 1e-230..1e-252 plus 1e-200..1e-230");
 auto cold=owner(patch(),300);auto pf=field(g,*cold,2,0);const auto base=cold->state_sha256();
 auto one=prepared(*cold,pf.f,1,"tiny-positive-area");const double n1=one.expected_number();j.d("area_1_count",n1).d("area_1_log_count",one.log_expected_number()).d("log_current_A_m2",one.log_current_A_m2()).d("mean_metal_energy_eV",one.mean_metal_energy_eV());
 exact_dump.s("exposure_log_current",hexs(one.log_current_A_m2())).s("exposure_mean_eV",hexs(one.mean_metal_energy_eV())).s("exposure_area",hexs(.005*.025)).s("exposure_dt",hexs(1e-14)).s("qe",hexs(se::qe));
 std::vector<std::string>rows,dump;double worst=0;std::string worst_a;long admitted=0,admitted_bad=0,refused=0,typed=0;
 for(int k=-60;k<=44;++k){const double a=std::pow(10.,-230.-.5*k);Obj r;r.d("area_fraction",a);
  std::optional<s5::PreparedPatch>p;auto o=attempt([&]{p=prepared(*cold,pf.f,a,"tiny-positive-area");});r.raw("prepare_patch",record(o));
  Obj d;d.s("area_fraction",hexs(a)).b("admitted",bool(p));
  if(p){++admitted;const double n=p->expected_number(),l=p->log_expected_number();
   const long double exact=std::exp((long double)one.log_current_A_m2())*(long double)(.005*.025)*(long double)a*(long double)1e-14/(long double)se::qe;
   const double rel=double(std::fabs(((long double)n-exact)/exact)),rel_lin=double(std::fabs(((long double)n-(long double)n1*a)/((long double)n1*a)));
   r.d("scalar_count",n).d("log_count",l).d("relative_error_vs_exact_finite_log",rel).d("relative_error_vs_linear_scaling",rel_lin).d("diagnostic_J",p->mean_wall_heat_diagnostic_J());
   d.s("count",hexs(n)).s("diagnostic",hexs(p->mean_wall_heat_diagnostic_J()));
   if(rel>worst){worst=rel;worst_a=hexs(a);}if(rel>T5)++admitted_bad;}
  else{++refused;if(scoped(o)&&o.kind=="exposure-scalar")++typed;d.s("quantity",o.quantity).s("observed",hexs(o.observed));}
  rows.push_back(r.str());dump.push_back(d.str());}
 exact_dump.raw("exposure_sweep",jarray(dump));
 j.raw("rows",jarray(rows)).i("admitted_area_fractions",admitted).i("refused_area_fractions",refused).i("typed_scoped_refusals",typed).i("admitted_with_relative_error_above_tolerance",admitted_bad).d("worst_admitted_relative_error",worst).s("worst_area_fraction",worst_a).d("tolerance",double(T5)).b("owner_unchanged",cold->state_sha256()==base);
 finish(j,id,admitted_bad?"gap: nonzero scalar count admitted outside tolerance":typed!=refused?"refusal without typed scope":"every admitted scalar within tolerance of the exact finite-log count; every refusal typed with scope");return j.str();
}
}

int main(int argc,char**argv){
 try{
  if(argc!=2){std::cerr<<"usage: source_v5_port <output-directory>\n";return 2;}
  std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
  std::cout<<std::setprecision(17)<<std::boolalpha;
  cpu_set_t cpus;CPU_ZERO(&cpus);sched_getaffinity(0,sizeof(cpus),&cpus);const int nice=getpriority(PRIO_PROCESS,0);
  std::cout<<"cpu="<<(CPU_ISSET(6,&cpus)?6:CPU_ISSET(7,&cpus)?7:-1)<<" affinity_count="<<CPU_COUNT(&cpus)<<" nice="<<nice<<"\n";
  if(!(CPU_COUNT(&cpus)==1&&(CPU_ISSET(6,&cpus)||CPU_ISSET(7,&cpus))&&nice==15)){std::cerr<<"CPU6-or-7/nice15 resource contract not met\n";return 3;}
  axial::Mesh g(45,4,.009,1);
  std::vector<std::string>probes{c1_release(g),c1_arithmetic(),c2_allocation(g),ext_overdraw(g),c4_tiny_birth(g),l1_exposure(g),handoff_control(g),a02_weight(g),a02_spectrum(g),a02_aba(g),a02_mesh(),ext_publish_gap(),ext_weight_sweep(g),ext_exposure_sweep(g),ext_handoff_current(g),ext_exposure_scalar(g)};
  std::vector<std::string>quoted;for(auto&s:summary)quoted.push_back(quote(s));
  {std::ofstream f(out/"port-results.json");f<<Obj().s("schema","source-v5-audit-reproducer-port-v1").s("reproducer","reviews/source/source_v3_reproducer.cpp").s("derived_from","artifacts/audit-reproducer-ports-20261030/source-v4/port-source/source_v4_port.cpp").s("target","electrode-heat-source-v5 (runtimes/pic-electrode-heat-source-v5-20261030)").d("relative_tolerance",s5::relative_tolerance).b("full_joint_test",false).raw("summary",jarray(quoted)).raw("probes",jarray(probes)).str()<<"\n";}
  {std::ofstream f(out/"exact-counts.json");f<<exact_dump.str()<<"\n";}
  std::cout<<"probes="<<probes.size()<<"\n";return 0;
 }catch(const std::exception&e){std::cerr<<"harness error: "<<e.what()<<"\n";return 1;}
}
