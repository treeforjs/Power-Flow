#pragma once
#include "dyadic.hpp"
#include "inputs/exact_particle_count.h"
#include "material/ss304_thermal_v2.hpp"
#include <openssl/sha.h>
#include <array>
#include <functional>
#include <iomanip>
#include <map>
#include <optional>
#include <set>
#include <sstream>

namespace ablation_proposal {
using metal_transfer::D;
using metal_transfer::need;
using Count=raikou::mfem_backend::ExactParticleCount;
using J=nlohmann::json;
inline D zero(){return D(0.);}
inline bool same(const D&a,const D&b){return a.compare(b)==0;}
inline bool finite(double x){return std::isfinite(x);}
inline std::string hash(const std::string&s){unsigned char a[SHA256_DIGEST_LENGTH];SHA256(reinterpret_cast<const unsigned char*>(s.data()),s.size(),a);std::ostringstream o;for(auto c:a)o<<std::hex<<std::setfill('0')<<std::setw(2)<<unsigned(c);return o.str();}
inline bool digest(const std::string&s){return s.size()==64&&s.find_first_not_of("0123456789abcdef")==std::string::npos;}
inline D atoms(const Count&c){D n;for(double w:c.positive_carriers())n=n+D(w);return n;}
inline D abs(const D&d){return d.compare(zero())<0?zero()-d:d;}
struct Caloric {
 std::string id,source_sha256;
 double Tmin=0,Tmax=0;
 std::function<double(double)> u_J_kg,cp_J_kg_K;
};
inline Caloric sealed_ss304(ss304v2::Qualification q,double rho,const std::string&seal){
 need(seal=="1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d","sealed material identity");
 auto m=std::make_shared<const ss304v2::Model>(q,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,rho);
 return {"sealed-SS304-q"+std::to_string(int(q)),seal,298.2,m->max_temperature_k(),[m](double T){return m->evaluate(T).u_model_j_kg;},[m](double T){return m->evaluate(T).apparent_heat_capacity_j_kg_k;}};
}
enum class Kind {Condensed,Adsorbate,Water};
enum class Trigger {MeltBudget,IndependentDesorption};
struct Model {
 std::string id,material_sha256,reference_id,reference_law_sha256;
 Kind kind=Kind::Condensed;
 Trigger trigger=Trigger::MeltBudget;
 std::optional<double> unit_mass_kg,reference_density_kg_m3,condensed_reference_J_kg;
 std::optional<double> trigger_K,removal_J_kg_excluding_actual_K;
 std::optional<double> source_Tmin,source_Tmax;
 std::optional<Caloric> caloric;
 // Derived bound-reference energy from the external canonical count, not a
 // second reservoir. The common owner supplies the immutable admitted law.
 std::function<D(const Count&)> adsorbed_reference_J;
 std::string support_id;
 std::map<std::string,double> parent_elements,product_elements;
 std::string product_state,product_reference_id;
 std::optional<double> product_mass_kg,product_R_J,product_I_J;
 std::optional<int> product_charge_e;
};
struct Scenario {
 std::string id,common_reference_id,mesh_sha256,geometry_sha256;
 std::string qualification; // e.g. manufactured_verification_only, explicit_reduction
 bool stationary_fixed_support=false;
 std::vector<Model> models;
 std::vector<std::string> condensed_order; // substrate -> outer covering layers
};
// Borrowed values from the ONE common owner/pre-removal trial. Condensed mass
// and adsorbate count are mutually exclusive encodings, never independent pools.
struct Stock {std::string id;D mass_kg;Count count;};
struct View {
 std::string accepted_owner_sha256,common_trial_sha256,field_sha256,source_sha256;
 std::uint64_t accepted_revision=0;
 double begin_s=0,end_s=0,dt_s=0;
 D thermal_energy_J; // already includes other admitted heat/work, exactly once
 std::vector<Stock> stocks;
 std::set<std::string> consumed_event_ids;
};
struct Carrier {
 std::string planned_id;
 double weight=0,kinetic_per_unit_J=0;
 std::array<double,3> proper_velocity_m_s{};
};
struct Action {
 std::string event_id,stock_id;
 std::vector<Carrier> carriers;
 D other_energy_out_J;
 std::string other_energy_owner;
};
struct Request {std::string expected_view_sha256,expected_exposure_sha256;std::vector<Action> actions;};
struct BudgetRefusal:std::invalid_argument {
 J evidence;
 explicit BudgetRefusal(J e):std::invalid_argument("insufficient local removal budget; proposed weights unchanged"),evidence(std::move(e)){}
};
class Proposal {
 friend class Evaluator;
 std::string base_view_sha256,scenario_sha256;
 std::vector<Stock> candidate_stocks;
 D candidate_thermal_energy_J;
 J receipt;
 Proposal()=default;
 public:
 Proposal(const Proposal&)=default;Proposal(Proposal&&) noexcept=default;
 const std::vector<Stock>& stocks()const{return candidate_stocks;}
 const D& thermal_energy()const{return candidate_thermal_energy_J;}
 const J& evidence()const{return receipt;}
};
inline double kinetic(double mass,const std::array<double,3>&u){double u2=0;for(double x:u){need(finite(x),"nonfinite proposed proper velocity");u2+=x*x;}need(finite(u2),"proposed proper velocity overflow");double K=mass*u2/(std::sqrt(1+u2/(299792458.*299792458.))+1);need(finite(K)&&K>=0,"proposed kinetic energy overflow");return K;}

// Stateless, bounded evaluator. No accepted state, publisher, RNG, gas ID
// allocator, temperature cache, or persistent count owner is implemented here.
class Evaluator {
 const Scenario c_;
 std::map<std::string,const Model*> models_;
 double Tmin_=0,Tmax_=INFINITY;
 const Model& model(const std::string&id)const{auto i=models_.find(id);need(i!=models_.end(),"unknown layer model");return *i->second;}
 static J stock_json(const Stock&s){return {{"id",s.id},{"mass_kg",s.mass_kg.json()},{"count_positive_carriers",s.count.positive_carriers()}};}
 static std::map<std::string,std::size_t> index(const std::vector<Stock>&stocks){std::map<std::string,std::size_t> r;for(std::size_t i=0;i<stocks.size();++i)need(r.emplace(stocks[i].id,i).second,"duplicate stock id");return r;}
 bool positive(const Stock&s)const{return model(s.id).kind==Kind::Condensed?s.mass_kg.compare(zero())>0:!s.count.zero();}
 D reference(const Stock&s)const{const auto&m=model(s.id);return m.kind==Kind::Condensed?s.mass_kg*D(*m.condensed_reference_J_kg):m.adsorbed_reference_J(s.count);}
 D thermal_at(const std::vector<Stock>&stocks,double T)const{
  need(finite(T)&&T>=Tmin_&&T<=Tmax_,"common caloric temperature outside admitted domain");D E;
  for(auto&s:stocks){const auto&m=model(s.id);if(m.kind!=Kind::Condensed)continue;double u=m.caloric->u_J_kg(T),cp=m.caloric->cp_J_kg_K(T);need(finite(u)&&finite(cp)&&cp>0,"invalid admitted caloric callback");E=E+s.mass_kg*D(u);}return E;
 }
 double temperature(const std::vector<Stock>&stocks,const D&E)const{
  D low=thermal_at(stocks,Tmin_),high=thermal_at(stocks,Tmax_);need(high.compare(low)>0&&E.compare(low)>=0&&E.compare(high)<=0,"local exact thermal energy outside admitted range");if(same(E,low))return Tmin_;if(same(E,high))return Tmax_;double a=Tmin_,b=Tmax_;
  for(int i=0;i<80;++i){double mid=a+(b-a)/2;if(mid==a||mid==b)break;if(thermal_at(stocks,mid).compare(E)<=0)a=mid;else b=mid;}
  return abs(thermal_at(stocks,a)-E).compare(abs(thermal_at(stocks,b)-E))<=0?a:b;
 }
 double admit_temperature(const std::vector<Stock>&stocks,const D&E)const{double T=temperature(stocks,E);for(auto&s:stocks){const auto&m=model(s.id);if(m.kind!=Kind::Condensed&&!s.count.zero())need(T>=*m.source_Tmin&&T<=*m.source_Tmax,"remaining adsorbate outside explicit source range");}return T;}
 void admit_stocks(const std::vector<Stock>&stocks)const{
  need(stocks.size()==models_.size(),"complete external layer inventory required");auto ix=index(stocks);
  for(auto&s:stocks){const auto&m=model(s.id);if(m.kind==Kind::Condensed){need(s.count.zero()&&s.mass_kg.compare(zero())>=0,"condensed inventory must be canonical mass only");}else{need(same(s.mass_kg,zero()),"adsorbate inventory must be canonical count only");(void)s.count.require_two_limb_values();need(ix.count(m.support_id)&&model(m.support_id).kind==Kind::Condensed,"adsorbate support identity");if(!s.count.zero())need(positive(stocks.at(ix.at(m.support_id))),"retained adsorbate has no material support");}}
  need(positive(stocks.at(ix.at(c_.condensed_order.front()))),"empty substrate requires common geometry handoff");
 }
 J model_json(const Model&m)const{return {{"id",m.id},{"kind",int(m.kind)},{"trigger",int(m.trigger)},{"material_sha256",m.material_sha256},{"reference_id",m.reference_id},{"reference_law_sha256",m.reference_law_sha256},{"unit_mass_kg",*m.unit_mass_kg},{"rho",m.reference_density_kg_m3?J(*m.reference_density_kg_m3):J()},{"condensed_reference_J_kg",m.condensed_reference_J_kg?J(*m.condensed_reference_J_kg):J()},{"trigger_K",m.trigger_K?J(*m.trigger_K):J()},{"removal_J_kg_excluding_actual_K",m.removal_J_kg_excluding_actual_K?J(*m.removal_J_kg_excluding_actual_K):J()},{"source_Tmin",m.source_Tmin?J(*m.source_Tmin):J()},{"source_Tmax",m.source_Tmax?J(*m.source_Tmax):J()},{"support",m.support_id},{"parent_elements",m.parent_elements},{"product_elements",m.product_elements},{"caloric_id",m.caloric?J(m.caloric->id):J()},{"caloric_source",m.caloric?J(m.caloric->source_sha256):J()},{"caloric_range",m.caloric?J::array({m.caloric->Tmin,m.caloric->Tmax}):J()},{"product",{m.product_state,m.product_reference_id,*m.product_mass_kg,*m.product_R_J,*m.product_I_J,*m.product_charge_e}}};}
 public:
 Evaluator(const Evaluator&)=delete;Evaluator&operator=(const Evaluator&)=delete;Evaluator(Evaluator&&)=delete;Evaluator&operator=(Evaluator&&)=delete;
 explicit Evaluator(const Scenario&c):c_(c){
  need(!c.id.empty()&&!c.common_reference_id.empty()&&!c.qualification.empty()&&digest(c.mesh_sha256)&&digest(c.geometry_sha256)&&c.stationary_fixed_support,"explicit scenario/reference/mesh and fixed stationary support required");need(!c.models.empty()&&c.models.size()<=64&&!c.condensed_order.empty()&&c.condensed_order.size()<=c.models.size(),"bounded explicit material stack required");
  for(const auto&m:c_.models){need(!m.id.empty()&&m.id.size()<=128&&models_.emplace(m.id,&m).second&&digest(m.material_sha256)&&digest(m.reference_law_sha256)&&m.reference_id==c.common_reference_id,"explicit per-material/reference identities required");
   need((m.kind==Kind::Condensed||m.kind==Kind::Adsorbate||m.kind==Kind::Water)&&!m.parent_elements.empty()&&m.parent_elements==m.product_elements,"explicit identity-preserving elemental composition required; reactions unsupported");for(auto&[element,n]:m.parent_elements)need(!element.empty()&&finite(n)&&n>0,"invalid elemental composition");
   need(m.unit_mass_kg&&finite(*m.unit_mass_kg)&&*m.unit_mass_kg>0&&!m.product_state.empty()&&m.product_reference_id==c.common_reference_id&&m.product_mass_kg&&*m.product_mass_kg==*m.unit_mass_kg&&m.product_R_J&&finite(*m.product_R_J)&&m.product_I_J&&finite(*m.product_I_J)&&*m.product_I_J>=0&&m.product_charge_e&&*m.product_charge_e==0,"explicit compatible neutral product mass/reference/internal state required");
   if(m.kind==Kind::Condensed){need(m.trigger==Trigger::MeltBudget&&m.caloric&&digest(m.caloric->source_sha256)&&!m.caloric->id.empty()&&m.caloric->u_J_kg&&m.caloric->cp_J_kg_K&&finite(m.caloric->Tmin)&&finite(m.caloric->Tmax)&&m.caloric->Tmin<m.caloric->Tmax,"explicit condensed caloric/trigger contract required");need(m.reference_density_kg_m3&&finite(*m.reference_density_kg_m3)&&*m.reference_density_kg_m3>0&&m.condensed_reference_J_kg&&finite(*m.condensed_reference_J_kg)&&m.trigger_K&&finite(*m.trigger_K)&&*m.trigger_K>=m.caloric->Tmin&&*m.trigger_K<=m.caloric->Tmax&&m.removal_J_kg_excluding_actual_K&&finite(*m.removal_J_kg_excluding_actual_K)&&*m.removal_J_kg_excluding_actual_K>0,"explicit density/reference/trigger/positive removal cost required");Tmin_=std::max(Tmin_,m.caloric->Tmin);Tmax_=std::min(Tmax_,m.caloric->Tmax);
    const double u=m.caloric->u_J_kg(*m.trigger_K);D declared=D(*m.unit_mass_kg)*D(*m.removal_J_kg_excluding_actual_K);D actual=D(*m.product_R_J)+D(*m.product_I_J)-D(*m.unit_mass_kg)*(D(*m.condensed_reference_J_kg)+D(u));D activity=abs(D(*m.product_R_J))+abs(D(*m.product_I_J))+abs(D(*m.unit_mass_kg)*D(*m.condensed_reference_J_kg))+abs(D(*m.unit_mass_kg)*D(u))+abs(declared);need(actual.compare(zero())>0&&abs(actual-declared).compare(activity*D(8*std::numeric_limits<double>::epsilon()))<=0,"removal cost and product/reference states disagree; no latent retuning");
   }else{need(m.trigger==Trigger::IndependentDesorption&&!m.trigger_K&&!m.removal_J_kg_excluding_actual_K&&!m.caloric&&!m.condensed_reference_J_kg&&!m.reference_density_kg_m3&&m.adsorbed_reference_J&&!m.support_id.empty()&&m.source_Tmin&&m.source_Tmax&&finite(*m.source_Tmin)&&finite(*m.source_Tmax)&&*m.source_Tmin<*m.source_Tmax,"independent desorption/reference/range required; no melt gate or duplicate latent cost");need(same(m.adsorbed_reference_J(Count()),zero()),"empty adsorbate reference must be zero");}
  }
  need(finite(Tmin_)&&finite(Tmax_)&&Tmin_<Tmax_,"empty common caloric domain");std::set<std::string>order;for(auto&id:c.condensed_order)need(order.insert(id).second&&model(id).kind==Kind::Condensed,"unique explicit condensed order");for(auto&[id,m]:models_)if(m->kind==Kind::Condensed)need(order.count(id),"condensed layer absent from exposure stack");
 }
 J exposure(const std::vector<Stock>&stocks)const{
  admit_stocks(stocks);auto ix=index(stocks);std::string outer;for(auto&id:c_.condensed_order)if(positive(stocks.at(ix.at(id))))outer=id;J bound=J::array();for(auto&s:stocks){const auto&m=model(s.id);if(m.kind!=Kind::Condensed&&!s.count.zero()&&m.support_id==outer)bound.push_back({{"id",s.id},{"count",s.count.positive_carriers()}});}return {{"outer_structural_material",outer},{"bound_adsorbates",bound}};
 }
 std::string scenario_sha256()const{J models=J::array();for(auto&m:c_.models)models.push_back(model_json(m));return hash(J({{"source",ABLATION_PROPOSAL_SOURCE_SHA256},{"id",c_.id},{"qualification",c_.qualification},{"reference",c_.common_reference_id},{"mesh",c_.mesh_sha256},{"geometry",c_.geometry_sha256},{"fixed_stationary_support",true},{"models",models},{"stack",c_.condensed_order}}).dump());}
 std::string view_sha256(const View&v)const{J stocks=J::array();for(auto&s:v.stocks)stocks.push_back(stock_json(s));return hash(J({{"scenario",scenario_sha256()},{"owner",v.accepted_owner_sha256},{"trial",v.common_trial_sha256},{"field",v.field_sha256},{"source",v.source_sha256},{"revision",v.accepted_revision},{"clock",{v.begin_s,v.end_s,v.dt_s}},{"thermal_energy_J",v.thermal_energy_J.json()},{"stocks",stocks},{"consumed",v.consumed_event_ids}}).dump());}
 D energy_at(const std::vector<Stock>&s,double T)const{admit_stocks(s);return thermal_at(s,T);}
 double projected_temperature(const View&v)const{admit_stocks(v.stocks);return admit_temperature(v.stocks,v.thermal_energy_J);}
 Proposal evaluate(const View&v,const Request&r)const{
  admit_stocks(v.stocks);for(auto*p:{&v.accepted_owner_sha256,&v.common_trial_sha256,&v.field_sha256,&v.source_sha256})need(digest(*p),"missing common-owner trial binding");need(v.accepted_revision<UINT64_MAX&&finite(v.begin_s)&&finite(v.end_s)&&finite(v.dt_s)&&v.dt_s>0&&v.end_s==v.begin_s+v.dt_s&&v.end_s>v.begin_s,"invalid authoritative common clock");const auto base=view_sha256(v);need(r.expected_view_sha256==base&&r.expected_exposure_sha256==hash(exposure(v.stocks).dump())&&r.actions.size()<=256,"stale trial/exposure or proposal resource cap");
  Proposal p;p.base_view_sha256=base;p.scenario_sha256=scenario_sha256();p.candidate_stocks=v.stocks;p.candidate_thermal_energy_J=v.thermal_energy_J;auto ix=index(v.stocks);const double initialT=admit_temperature(v.stocks,v.thermal_energy_J);D gasE,totalR,totalX,totalMass;std::array<D,3>support;J rows=J::array();auto events=v.consumed_event_ids;std::set<std::string>carrierids;std::size_t carrier_count=0;
  for(auto&a:r.actions){need(!a.event_id.empty()&&a.event_id.size()<=256&&events.insert(a.event_id).second&&!a.carriers.empty()&&a.carriers.size()<=Count::maximum_carriers,"duplicate/missing release event or carrier cap");need(carrier_count<=4096-a.carriers.size(),"total proposal carrier cap");carrier_count+=a.carriers.size();auto mi=ix.find(a.stock_id);need(mi!=ix.end(),"unknown external source inventory");auto&stock=p.candidate_stocks[mi->second];const auto&m=model(stock.id);J exp=exposure(p.candidate_stocks);std::string outer=exp.at("outer_structural_material");need(m.kind==Kind::Condensed?stock.id==outer:m.support_id==outer,"buried layer/source release requires an explicit transport or exposure model");need(a.other_energy_out_J.compare(zero())>=0&&(same(a.other_energy_out_J,zero())||!a.other_energy_owner.empty()),"unowned/negative additional removal energy");
   D beforeR=reference(stock),n,K,E;std::vector<double>weights;J carriers=J::array();for(auto&c:a.carriers){need(!c.planned_id.empty()&&c.planned_id.size()<=256&&carrierids.insert(c.planned_id).second&&finite(c.weight)&&c.weight>0&&finite(c.kinetic_per_unit_J)&&c.kinetic_per_unit_J>=0,"invalid/duplicate proposed carrier");double expectedK=kinetic(*m.unit_mass_kg,c.proper_velocity_m_s);D scale=abs(D(expectedK))+abs(D(c.kinetic_per_unit_J));need(abs(D(expectedK)-D(c.kinetic_per_unit_J)).compare(D(8*std::numeric_limits<double>::epsilon())*scale)<=0,"proposed K/proper momentum mismatch");weights.push_back(c.weight);n=n+D(c.weight);K=K+D(c.weight)*D(c.kinetic_per_unit_J);E=E+D(c.weight)*(D(*m.product_R_J)+D(*m.product_I_J)+D(c.kinetic_per_unit_J));for(int d=0;d<3;++d)support[d]=support[d]-D(c.weight)*D(*m.unit_mass_kg)*D(c.proper_velocity_m_s[d]);carriers.push_back({{"planned_id",c.planned_id},{"weight",c.weight},{"kinetic_per_unit_J",c.kinetic_per_unit_J},{"proper_velocity_m_s",c.proper_velocity_m_s}});}
   const auto requestedCount=Count::from_positive_carriers(weights);need(same(atoms(requestedCount),n),"exact proposed positive leaf sum");D dm=n*D(*m.unit_mass_kg);D threshold,available,incremental;
   if(m.kind==Kind::Condensed){threshold=thermal_at(p.candidate_stocks,*m.trigger_K);available=p.candidate_thermal_energy_J-threshold;need(available.compare(zero())>=0,"melt trigger not reached in local net trial energy");D remaining=stock.mass_kg-dm;need(remaining.compare(zero())>=0,"finite condensed inventory overdraw");if(stock.id==c_.condensed_order.front())need(remaining.compare(zero())>0,"full substrate depletion requires geometry handoff");if(same(remaining,zero()))for(auto&other:p.candidate_stocks)if(model(other.id).kind!=Kind::Condensed&&model(other.id).support_id==stock.id)need(other.count.zero(),"cannot remove support with retained adsorbate; no silent entrainment/deletion");stock.mass_kg=remaining;
   }else{need(initialT>=*m.source_Tmin&&initialT<=*m.source_Tmax,"independent desorption outside declared source range");for(double w:weights)stock.count=stock.count.debit(w);(void)stock.count.require_two_limb_values();}
   D dR=reference(stock)-beforeR;D debit=E+dR+a.other_energy_out_J;
   if(m.kind==Kind::Condensed){incremental=debit-dm*D(m.caloric->u_J_kg(*m.trigger_K));need(incremental.compare(zero())>0,"nonpositive actual incremental removal cost");if(available.compare(incremental)<0)throw BudgetRefusal({{"stock",stock.id},{"available_J",available.json()},{"requested_incremental_J",incremental.json()},{"plan_scale_limit_ratio",{{"numerator",available.json()},{"denominator",incremental.json()}}},{"proposed_weights_unchanged",true}});}
   p.candidate_thermal_energy_J=p.candidate_thermal_energy_J-debit;admit_stocks(p.candidate_stocks);(void)admit_temperature(p.candidate_stocks,p.candidate_thermal_energy_J);gasE=gasE+E;totalR=totalR+dR;totalX=totalX+a.other_energy_out_J;totalMass=totalMass+dm;
   rows.push_back({{"event",a.event_id},{"stock",stock.id},{"trigger",m.kind==Kind::Condensed?"explicit_melt_energy_budget":"independent_desorption_no_melt_gate"},{"count",n.json()},{"mass_removed_kg",dm.json()},{"gas_energy_J",E.json()},{"reference_change_J",dR.json()},{"other_energy_out_J",a.other_energy_out_J.json()},{"other_energy_owner",a.other_energy_owner},{"incremental_cost_including_actual_K_J",incremental.json()},{"available_above_trigger_J",available.json()},{"plan_scale_limit_ratio",m.kind==Kind::Condensed?J({{"numerator",available.json()},{"denominator",incremental.json()}}):J()},{"trigger_thermal_energy_J",threshold.json()},{"product",m.product_state},{"product_R_J",*m.product_R_J},{"product_I_J",*m.product_I_J},{"unit_mass_kg",*m.unit_mass_kg},{"carriers",carriers},{"stock_after",stock_json(stock)},{"removed_reference_volume_ratio",m.kind==Kind::Condensed?J({{"mass_numerator_kg",dm.json()},{"density_denominator_kg_m3",*m.reference_density_kg_m3}}):J()}});
  }
  D beforeMass,afterMass;D beforeRef,afterRef;for(auto&s:v.stocks){const auto&m=model(s.id);beforeMass=beforeMass+(m.kind==Kind::Condensed?s.mass_kg:atoms(s.count)*D(*m.unit_mass_kg));beforeRef=beforeRef+reference(s);}for(auto&s:p.candidate_stocks){const auto&m=model(s.id);afterMass=afterMass+(m.kind==Kind::Condensed?s.mass_kg:atoms(s.count)*D(*m.unit_mass_kg));afterRef=afterRef+reference(s);}
  D massResidual=afterMass-beforeMass+totalMass;D energyResidual=p.candidate_thermal_energy_J-v.thermal_energy_J+(afterRef-beforeRef)+gasE+totalX;need(same(massResidual,zero())&&same(energyResidual,zero())&&same(afterRef-beforeRef,totalR),"exact proposal mass/reference/energy closure");J before=J::array(),after=J::array(),impulse=J::array();for(auto&s:v.stocks)before.push_back(stock_json(s));for(auto&s:p.candidate_stocks)after.push_back(stock_json(s));for(auto&x:support)impulse.push_back(x.json());
  p.receipt={{"schema","melt-ablation-layer-proposal-v1"},{"source",ABLATION_PROPOSAL_SOURCE_SHA256},{"scope","stateless proposal only; no accepted thermal/count owner or native birth"},{"qualification",c_.qualification},{"base_view_sha256",base},{"scenario_sha256",p.scenario_sha256},{"reference_id",c_.common_reference_id},{"stocks_before",before},{"stocks_after",after},{"thermal_before_J",v.thermal_energy_J.json()},{"thermal_after_J",p.candidate_thermal_energy_J.json()},{"reference_before_J",beforeRef.json()},{"reference_after_J",afterRef.json()},{"outgoing_energy_J",gasE.json()},{"other_energy_out_J",totalX.json()},{"outgoing_mass_kg",totalMass.json()},{"exact_mass_residual",massResidual.json()},{"exact_energy_residual",energyResidual.json()},{"initial_projected_T_K",initialT},{"final_projected_T_K",admit_temperature(p.candidate_stocks,p.candidate_thermal_energy_J)},{"exposure_before",exposure(v.stocks)},{"exposure_after",exposure(p.candidate_stocks)},{"support_impulse_N_s",impulse},{"fixed_support_work_J",0},{"charge_change_C",0},{"actions",rows}};return p;
 }
 // Freshness only: the common owner must still bind actual registry/particle,
 // source, geometry, field, thermal and RNG candidates before atomic acceptance.
 bool matches_base(const View&v,const Proposal&p)const{return p.scenario_sha256==scenario_sha256()&&p.base_view_sha256==view_sha256(v);}
 // Configuration evidence only. Callback implementations remain admitted by the
 // common owner under their supplied hashes; this JSON grants no publication.
 J configuration_evidence()const{J out=J::array();for(auto&m:c_.models){J j=model_json(m);if(m.kind==Kind::Condensed)j["u_at_trigger_J_kg"]=m.caloric->u_J_kg(*m.trigger_K);out.push_back(j);}return out;}
};
} // namespace ablation_proposal
