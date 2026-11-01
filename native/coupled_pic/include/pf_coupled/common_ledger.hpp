#pragma once
#include "common_balance.hpp"
#include "native_genesis.hpp"

namespace powerflow::common_step_v9 {
namespace cm=::powerflow::common_material_stocks_v1;
namespace g=::pf_gas_v1;
namespace n=::raikou::mfem_backend;
namespace cg=::powerflow::coupled_genesis_v9;
inline void ledger_need(bool ok,const char*why){if(!ok)throw std::invalid_argument(std::string("common ledger: ")+why);}
inline g::Exact canonical_exact(const g::Exact&value){g::Exact out;const auto digits=value.normalised();out.add_digits(digits.data());return out;}
inline void accumulate(g::Exact&total,const g::Exact&delta){total.merge(delta);total=canonical_exact(total);}
inline void encode_exact(cg::Encoding&encoding,const g::Exact&value){for(const auto digit:value.normalised())encoding.u64(std::bit_cast<std::uint64_t>(digit));}
inline void encode_inventory(cg::Encoding&e,const cm::Inventory&v){for(const auto*p:{&v.kinetic_J,&v.reference_J,&v.internal_J,&v.binding_J,&v.charge_e,&v.hydrogen,&v.oxygen})encode_exact(e,*p);}
inline void accumulate_inventory(cm::Inventory&a,const cm::Inventory&b){
 for(auto member:{&cm::Inventory::kinetic_J,&cm::Inventory::reference_J,&cm::Inventory::internal_J,&cm::Inventory::binding_J,&cm::Inventory::charge_e,&cm::Inventory::hydrogen,&cm::Inventory::oxygen})accumulate(a.*member,b.*member);
}
inline bool same_inventory(const cm::Inventory&a,const cm::Inventory&b){for(auto member:{&cm::Inventory::kinetic_J,&cm::Inventory::reference_J,&cm::Inventory::internal_J,&cm::Inventory::binding_J,&cm::Inventory::charge_e,&cm::Inventory::hydrogen,&cm::Inventory::oxygen})if(!(a.*member==b.*member))return false;return true;}
struct CommonLedger {
 std::uint64_t accepted_step=0,events=0;
 double begin_s=0,time_s=0;
 cm::Inventory initial_material,accepted_material;
 g::Exact initial_energy_J,accepted_energy_J;
 g::Exact emitted_electrode_reference_J,emitted_electrode_charge_e;
 cm::BoundaryStock boundary;
 std::map<std::string,g::Exact> conductor_species,exterior_species;
 g::Exact external_net_J,event_defect_J,event_reference_conversion_J,thermal_reference_rounding_J;
 g::Exact source_kinematic_rounding_J,source_binding_conversion_J;
 g::Exact raw_residual_J,accounted_residual_J,admitted_bound_sum_J;
 std::string sha256;
};
inline std::string ledger_hash(const CommonLedger&v){
 cg::Encoding e;e.text("common-accepted-material-energy-ledger-v1");e.u64(v.accepted_step);e.u64(v.events);e.number(v.begin_s);e.number(v.time_s);
 encode_inventory(e,v.initial_material);encode_inventory(e,v.accepted_material);encode_inventory(e,v.boundary.conductor);encode_inventory(e,v.boundary.exterior);
 for(const auto*p:{&v.initial_energy_J,&v.accepted_energy_J,&v.emitted_electrode_reference_J,&v.emitted_electrode_charge_e,
  &v.boundary.short_kinetic_J,&v.boundary.exterior_kinetic_J,&v.boundary.exterior_correction_J,&v.boundary.foil_kinetic_J,
  &v.external_net_J,&v.event_defect_J,&v.event_reference_conversion_J,&v.thermal_reference_rounding_J,
  &v.source_kinematic_rounding_J,&v.source_binding_conversion_J,&v.raw_residual_J,&v.accounted_residual_J,&v.admitted_bound_sum_J})encode_exact(e,*p);
 for(const auto*m:{&v.conductor_species,&v.exterior_species}){e.u64(m->size());for(const auto&[name,count]:*m){e.text(name);encode_exact(e,count);}}
 return e.hash();
}
// This record is prepared with the actual sealed candidate before publication.
// Exact cumulative digits are normalized each step so history does not grow
// unbounded integer limbs; no accepted account changes during preparation.
inline std::shared_ptr<const CommonLedger> prepare_ledger(const std::shared_ptr<const CommonLedger>&before,
 const n::FullTimeStepBinding&common,const cm::CombinedAssessment&step,const cm::SourceAssessment&source,
 const cm::BoundaryStock&boundary,const cm::ChemistryOperands&chem,const n::FullTimeSinkLedger*sinks,
 const n::FullTimeCandidatePopulation&catalogue,double admitted_energy_bound_J){
 ledger_need(common.accepted_step<UINT64_MAX&&std::isfinite(admitted_energy_bound_J)&&admitted_energy_bound_J>=0,"valid bounded transition required");
 auto next=std::make_shared<CommonLedger>();
 if(before){
  ledger_need(before->accepted_step==common.accepted_step&&before->time_s==common.begin_s&&before->sha256==ledger_hash(*before),"stale or altered accepted common ledger");
  ledger_need(before->accepted_energy_J==step.physical_stored_before_J&&same_inventory(before->accepted_material,step.material_before),"actual before stocks do not equal previous accepted stocks");
  *next=*before;
 }else{
  ledger_need(common.accepted_step==0,"missing non-genesis common ledger");next->begin_s=common.begin_s;
  next->initial_energy_J=step.physical_stored_before_J;next->initial_material=step.material_before;
  ledger_need(!catalogue.states.empty()&&catalogue.states.size()<=64,"bounded actual species catalogue");
  for(const auto&s:catalogue.states){ledger_need(next->conductor_species.emplace(s.id,g::Exact{}).second,"duplicate species");next->exterior_species.emplace(s.id,g::Exact{});}
 }
 next->accepted_step=common.accepted_step+1;next->time_s=common.end_s;
 next->accepted_energy_J=step.physical_stored_after_J;next->accepted_material=step.material_after;
 ledger_need(chem.events<=UINT64_MAX-next->events,"common event count overflow");next->events+=chem.events;
 accumulate(next->emitted_electrode_reference_J,source.emitted_electrode_reference_J);accumulate(next->emitted_electrode_charge_e,source.emitted_electrode_charge_e);
 accumulate_inventory(next->boundary.conductor,boundary.conductor);accumulate_inventory(next->boundary.exterior,boundary.exterior);
 for(auto member:{&cm::BoundaryStock::short_kinetic_J,&cm::BoundaryStock::exterior_kinetic_J,&cm::BoundaryStock::exterior_correction_J,&cm::BoundaryStock::foil_kinetic_J})accumulate(next->boundary.*member,boundary.*member);
 if(sinks)for(const auto&sink:sinks->exposed){auto&map=sink.record.kind==n::FullTimeSinkKind::escaped_exterior?next->exterior_species:next->conductor_species;
  auto found=map.find(sink.state);ledger_need(found!=map.end(),"sink species absent from accepted catalogue");g::Exact count;count.add_product(sink.weight,1);accumulate(found->second,count);}
 accumulate(next->external_net_J,step.external_net_J);accumulate(next->event_defect_J,step.event_numerical_defect_J);
 accumulate(next->event_reference_conversion_J,step.event_reference_conversion_J);accumulate(next->thermal_reference_rounding_J,step.thermal_reference_rounding_J);
 accumulate(next->source_kinematic_rounding_J,step.source_kinematic_rounding_J);accumulate(next->source_binding_conversion_J,step.source_binding_operand_conversion_J);
 g::Exact bound;bound.add_product(admitted_energy_bound_J,1);accumulate(next->admitted_bound_sum_J,bound);
 const auto previous_raw=next->raw_residual_J,previous_accounted=next->accounted_residual_J;
 next->raw_residual_J=cm::diff(next->accepted_energy_J,next->initial_energy_J);next->raw_residual_J.merge(next->emitted_electrode_reference_J);
 next->raw_residual_J.merge(next->boundary.conductor.energy());next->raw_residual_J.merge(next->boundary.exterior.energy());
 next->raw_residual_J.merge(next->boundary.short_kinetic_J);next->raw_residual_J.merge(next->boundary.exterior_kinetic_J);next->raw_residual_J.merge(next->boundary.exterior_correction_J);next->raw_residual_J.merge(next->external_net_J,1);
 next->accounted_residual_J=next->raw_residual_J;next->accounted_residual_J.merge(next->event_defect_J);next->accounted_residual_J.merge(next->event_reference_conversion_J,1);
 next->accounted_residual_J.merge(next->thermal_reference_rounding_J,1);next->accounted_residual_J.merge(next->source_kinematic_rounding_J);next->accounted_residual_J.merge(next->source_binding_conversion_J,1);
 auto raw_expected=previous_raw;raw_expected.merge(step.raw_residual_J);auto accounted_expected=previous_accounted;accounted_expected.merge(step.accounted_residual_J);
 ledger_need(next->raw_residual_J==raw_expected&&next->accounted_residual_J==accounted_expected,"whole-history energy accounts do not telescope");
 auto H=cm::diff(next->accepted_material.hydrogen,next->initial_material.hydrogen);H.merge(next->boundary.conductor.hydrogen);H.merge(next->boundary.exterior.hydrogen);
 auto O=cm::diff(next->accepted_material.oxygen,next->initial_material.oxygen);O.merge(next->boundary.conductor.oxygen);O.merge(next->boundary.exterior.oxygen);
 auto Q=cm::diff(next->accepted_material.charge_e,next->initial_material.charge_e);Q.merge(next->boundary.conductor.charge_e);Q.merge(next->boundary.exterior.charge_e);Q.merge(next->emitted_electrode_charge_e);
 ledger_need(H.zero()&&O.zero()&&Q.zero(),"whole-history material nuclei/charge imbalance");next->sha256=ledger_hash(*next);return next;
}
} // namespace powerflow::common_step_v9
