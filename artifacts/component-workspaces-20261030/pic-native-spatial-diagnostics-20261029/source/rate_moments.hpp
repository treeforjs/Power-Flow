#pragma once
#include "raikou/mfem_pic_arena.h"
#include "collision_scheduler.hpp"
#include <iomanip>
#include <sstream>
#include <unordered_map>

namespace native_diagnostics {
using namespace raikou::mfem_backend;
namespace pc=raikou::pic_collisions;
using Json=ArenaJson;
struct AssignedRateEntry { ArenaRateEntry entry;std::uint64_t physical_cell=0; };
struct ReactionCell {std::uint64_t id;double physical_volume_m3;};
struct SigmaCoverage {
 bool supported=false;double sigma_m2=0;std::string reason;
 // Inclusive yield or aggregate observables are not exclusive MCC events.
 bool exclusive_event_kernel_supported=false;
};
struct RateQuery {
 std::string channel,projectile_state,target_state,dataset_sha256,kernel_sha256,model,observable;
 std::function<SigmaCoverage(double target_frame_energy_eV)> evaluate;
};
struct Sum {
 long double high=0,low=0;
 void add(long double v){const auto t=high+v;low+=std::abs(high)>=std::abs(v)?(high-t)+v:(v-t)+high;high=t;}
 long double value()const{return high+low;}
};
inline Json numeric(long double v){
 pc::require(std::isfinite(v),"diagnostic moment overflow");const double d=double(v);std::ostringstream s;s<<std::setprecision(std::numeric_limits<long double>::max_digits10)<<v;
 return {{"value",std::isfinite(d)&&(v==0||d!=0)?Json(d):Json(nullptr)},{"long_double_decimal",s.str()},{"positive_below_double",v>0&&d==0}};
}
inline Json rate_moments(const std::vector<AssignedRateEntry>& entries,
 const std::vector<ReactionCell>& cells,const std::vector<RateQuery>& queries,
 const std::string& state_prefix,const std::string& physical_geometry_sha256,
 std::size_t maximum_pairs){
 pc::require(pc::digest(state_prefix)&&pc::digest(physical_geometry_sha256),"unbound rate diagnostic state/physical volume");
 pc::require(entries.size()<=2000000&&cells.size()<=1010000&&queries.size()<=64&&maximum_pairs<=1000000,"rate diagnostic resource cap");
 std::map<std::uint64_t,double> volumes;for(const auto&cell:cells){pc::positive(cell.physical_volume_m3,"physical3Dreaction volume");pc::require(volumes.emplace(cell.id,cell.physical_volume_m3).second,"duplicate reaction cell");}
 std::map<std::uint64_t,std::vector<const AssignedRateEntry*>> by_cell;std::set<std::uint64_t> ids;
 for(const auto&e:entries){const auto&p=e.entry.native;pc::positive(p.weight,"actual positive rate-carrier weight");pc::positive(p.mass_kg,"native mass");pc::require(p.id>0&&!e.entry.packet_id.empty()&&!e.entry.state.empty()&&ids.insert(p.id).second&&volumes.count(e.physical_cell),"invalid/duplicate native identity or unknown physical reaction cell");for(double v:p.position_m)pc::require(std::isfinite(v),"nonfinite native position");for(double v:p.proper_velocity_m_s)pc::require(std::isfinite(v),"nonfinite native proper velocity");by_cell[e.physical_cell].push_back(&e);}
 std::size_t pair_count=0;
 for(const auto&q:queries){pc::require(q.evaluate&&pc::digest(q.dataset_sha256)&&pc::digest(q.kernel_sha256)&&!q.channel.empty()&&!q.model.empty()&&!q.observable.empty()&&!q.projectile_state.empty()&&!q.target_state.empty()&&q.projectile_state!=q.target_state,"unbound query or unsupported identical-species pair convention");
  for(const auto&[cell,population]:by_cell){std::size_t a=0,b=0;for(auto*e:population){a+=e->entry.state==q.projectile_state;b+=e->entry.state==q.target_state;}pc::require(!a||b<=(maximum_pairs-pair_count)/a,"pair cap; no partial/truncated rate report");pair_count+=a*b;}}
 Json rows=Json::array();
 for(const auto&q:queries)for(const auto&[cell,population]:by_cell){
  std::vector<const AssignedRateEntry*> projectiles,targets;Sum ne,nt,total_pair,supported_pair,unsupported_pair,event_kernel_missing,raw;
  for(auto*e:population){if(e->entry.state==q.projectile_state){projectiles.push_back(e);ne.add(e->entry.native.weight);}if(e->entry.state==q.target_state){targets.push_back(e);nt.add(e->entry.native.weight);}}
  if(projectiles.empty())continue;
  std::map<std::string,Sum> missing;std::set<std::uint64_t> unsupported_projectiles,unsupported_targets;double minE=INFINITY,maxE=0;
  for(auto*pe:projectiles)for(auto*pn:targets){const auto&a=pe->entry.native;const auto&b=pn->entry.native;const long double pair_weight=static_cast<long double>(a.weight)*b.weight;total_pair.add(pair_weight);
   pc::Packet p{pe->entry.packet_id,pe->entry.state,cell,a.weight,a.mass_kg,a.proper_velocity_m_s},n{pn->entry.packet_id,pn->entry.state,cell,b.weight,b.mass_kg,b.proper_velocity_m_s};pc::Kinematics k;
   try{k=pc::relative(p,n);}catch(const pc::Error&error){unsupported_pair.add(pair_weight);missing[std::string("relative_kinematics: ")+error.what()].add(pair_weight);unsupported_projectiles.insert(a.id);unsupported_targets.insert(b.id);continue;}
   minE=std::min(minE,k.target_frame_energy_eV);maxE=std::max(maxE,k.target_frame_energy_eV);const auto coverage=q.evaluate(k.target_frame_energy_eV);
   if(!coverage.supported){pc::require(!coverage.reason.empty(),"unsupported positive weight needs a reason");unsupported_pair.add(pair_weight);missing[coverage.reason].add(pair_weight);unsupported_projectiles.insert(a.id);unsupported_targets.insert(b.id);continue;}
   pc::require(std::isfinite(coverage.sigma_m2)&&coverage.sigma_m2>=0,"invalid supported cross section");supported_pair.add(pair_weight);raw.add(pair_weight*coverage.sigma_m2*k.moller_flux_m_s);if(coverage.sigma_m2>0&&!coverage.exclusive_event_kernel_supported)event_kernel_missing.add(pair_weight);
  }
  Sum affected,target_affected;for(auto*e:projectiles)if(unsupported_projectiles.count(e->entry.native.id))affected.add(e->entry.native.weight);for(auto*e:targets)if(unsupported_targets.count(e->entry.native.id))target_affected.add(e->entry.native.weight);
  Json reasons=Json::object();for(const auto&[reason,w]:missing)reasons[reason]=numeric(w.value());
  rows.push_back({{"cell_id",cell},{"channel",q.channel},{"observable",q.observable},{"model",q.model},{"dataset_sha256",q.dataset_sha256},{"kernel_sha256",q.kernel_sha256},{"physical_reaction_volume_m3",volumes.at(cell)},
   {"projectile_weight",numeric(ne.value())},{"target_weight",numeric(nt.value())},{"pair_weight",numeric(total_pair.value())},{"supported_pair_weight",numeric(supported_pair.value())},{"unsupported_pair_weight",numeric(unsupported_pair.value())},{"unsupported_pair_fraction",total_pair.value()>0?numeric(unsupported_pair.value()/total_pair.value()):Json(nullptr)},
   {"unsupported_projectile_weight",numeric(affected.value())},{"unsupported_projectile_fraction",numeric(affected.value()/ne.value())},{"unsupported_target_weight",numeric(target_affected.value())},{"unsupported_target_fraction",nt.value()>0?numeric(target_affected.value()/nt.value()):Json(nullptr)},{"missing_exclusive_event_kernel_pair_weight",numeric(event_kernel_missing.value())},{"missing_exclusive_event_kernel_pair_fraction",total_pair.value()>0?numeric(event_kernel_missing.value()/total_pair.value()):Json(nullptr)},{"unsupported_reasons",reasons},
   {"supported_partial_sum_WeWn_sigma_relative_flux_m3_s",numeric(raw.value())},{"supported_partial_physical_number_rate_s",numeric(raw.value()/volumes.at(cell))},{"supported_partial_sum_We_target_averaged_sigma_flux_m3_s",nt.value()>0?numeric(raw.value()/nt.value()):Json(nullptr)},
   {"target_frame_energy_min_eV",std::isfinite(minE)?Json(minE):Json(nullptr)},{"target_frame_energy_max_eV",std::isfinite(minE)?Json(maxE):Json(nullptr)},{"complete_sigma_coverage",unsupported_pair.value()==0&&nt.value()>0},{"status",nt.value()>0?"actual positive-weight pair moments":"no actual target population; no target velocity inferred"}});
 }
 return {{"schema","native_rate_diagnostics/v1"},{"state_prefix",state_prefix},{"physical_geometry_sha256",physical_geometry_sha256},{"actual_pair_count",pair_count},{"maximum_pairs",maximum_pairs},{"rows",rows},{"sum_convention","long-double compensated positive pair quadrature; canonical inventory counts are separately exact"},{"no_support_renormalization",true},{"no_observation_slab_in_rate_volume",true},{"full_mechanism_qualification",false},{"new_evolution_or_GPU",false}};
}
} // namespace native_diagnostics
