#pragma once
#include "common_ledger.hpp"
#include <string_view>

namespace powerflow::common_ledger_codec_v1 {
namespace cs=::powerflow::common_step_v9;
namespace cm=cs::cm;
namespace g=cs::g;
namespace n=cs::n;
inline constexpr std::size_t byte_cap=1024*1024;
inline void need(bool ok,const char*why){if(!ok)throw std::invalid_argument(std::string("common ledger codec: ")+why);}
inline bool nonnegative(const g::Exact&v){return v.normalised().back()>=0;}
// Arithmetic and canonical-byte checks only. Restore authority still requires
// the root coordinator's genuine native/gas continuation and restored owners.
inline void validate(const cs::CommonLedger&v){
 need(v.accepted_step>0&&std::isfinite(v.begin_s)&&std::isfinite(v.time_s)&&v.time_s>v.begin_s,"accepted clock");
 need(v.sha256.size()==64&&v.sha256==cs::ledger_hash(v),"content seal");
 need(v.conductor_species.size()>0&&v.conductor_species.size()<=64&&v.conductor_species.size()==v.exterior_species.size(),"bounded species accounts");
 for(const auto&[id,count]:v.conductor_species){auto found=v.exterior_species.find(id);need(!id.empty()&&id.size()<=4096&&found!=v.exterior_species.end()&&nonnegative(count)&&nonnegative(found->second),"species inventory shape/sign");}
 auto raw=cm::diff(v.accepted_energy_J,v.initial_energy_J);raw.merge(v.emitted_electrode_reference_J);
 raw.merge(v.boundary.conductor.energy());raw.merge(v.boundary.exterior.energy());raw.merge(v.boundary.short_kinetic_J);raw.merge(v.boundary.exterior_kinetic_J);raw.merge(v.boundary.exterior_correction_J);raw.merge(v.external_net_J,1);
 auto accounted=raw;accounted.merge(v.event_defect_J);accounted.merge(v.event_reference_conversion_J,1);accounted.merge(v.thermal_reference_rounding_J,1);accounted.merge(v.source_kinematic_rounding_J);accounted.merge(v.source_binding_conversion_J,1);
 need(raw==v.raw_residual_J&&accounted==v.accounted_residual_J,"cumulative energy accounts");
 auto H=cm::diff(v.accepted_material.hydrogen,v.initial_material.hydrogen);H.merge(v.boundary.conductor.hydrogen);H.merge(v.boundary.exterior.hydrogen);
 auto O=cm::diff(v.accepted_material.oxygen,v.initial_material.oxygen);O.merge(v.boundary.conductor.oxygen);O.merge(v.boundary.exterior.oxygen);
 auto Q=cm::diff(v.accepted_material.charge_e,v.initial_material.charge_e);Q.merge(v.boundary.conductor.charge_e);Q.merge(v.boundary.exterior.charge_e);Q.merge(v.emitted_electrode_charge_e);
 need(H.zero()&&O.zero()&&Q.zero(),"cumulative nuclei/charge accounts");
 need(nonnegative(v.admitted_bound_sum_J),"cumulative bound sign");
 auto lo=accounted,hi=accounted;lo.merge(v.admitted_bound_sum_J);hi.merge(v.admitted_bound_sum_J,1);
 need(nonnegative(lo)&&(hi.normalised().back()<0||hi.zero()),"cumulative derived energy bound");
}
inline std::string encode(const cs::CommonLedger&v){
 validate(v);cs::cg::Encoding e;e.text("common-ledger-checkpoint-v1");e.u64(v.accepted_step);e.u64(v.events);e.number(v.begin_s);e.number(v.time_s);
 for(const auto*p:{&v.initial_material,&v.accepted_material,&v.boundary.conductor,&v.boundary.exterior})cs::encode_inventory(e,*p);
 for(const auto*p:{&v.initial_energy_J,&v.accepted_energy_J,&v.emitted_electrode_reference_J,&v.emitted_electrode_charge_e,
  &v.boundary.short_kinetic_J,&v.boundary.exterior_kinetic_J,&v.boundary.exterior_correction_J,&v.boundary.foil_kinetic_J,
  &v.external_net_J,&v.event_defect_J,&v.event_reference_conversion_J,&v.thermal_reference_rounding_J,
  &v.source_kinematic_rounding_J,&v.source_binding_conversion_J,&v.raw_residual_J,&v.accounted_residual_J,&v.admitted_bound_sum_J})cs::encode_exact(e,*p);
 for(const auto*m:{&v.conductor_species,&v.exterior_species}){e.u64(m->size());for(const auto&[id,count]:*m){e.text(id);cs::encode_exact(e,count);}}
 e.text(v.sha256);need(e.bytes.size()<=byte_cap,"encoded byte cap");return e.bytes;
}
class Reader {
 std::string_view bytes_;std::size_t pos_=0;
public:
 explicit Reader(std::string_view bytes):bytes_(bytes){need(bytes.size()<=byte_cap,"input byte cap");}
 std::uint64_t u64(){need(pos_<=bytes_.size()&&bytes_.size()-pos_>=8,"truncated integer");std::uint64_t value=0;for(unsigned k=0;k<8;++k)value|=std::uint64_t(static_cast<unsigned char>(bytes_[pos_++]))<<(8*k);return value;}
 double number(){const auto value=std::bit_cast<double>(u64());need(std::isfinite(value),"nonfinite scalar");return value;}
 std::string text(std::size_t cap){const auto count=u64();need(count<=cap&&pos_<=bytes_.size()&&count<=bytes_.size()-pos_,"bounded text");auto result=std::string(bytes_.substr(pos_,count));pos_+=count;return result;}
 g::Exact exact(){
  std::array<std::int64_t,pfg::kDigits> digits{};
  for(std::size_t k=0;k<digits.size();++k){digits[k]=std::bit_cast<std::int64_t>(u64());need(k+1==digits.size()?(digits[k]==0||digits[k]==-1):(digits[k]>=0&&digits[k]<=INT64_C(0xffffffff)),"noncanonical or overflowing exact digit");}
  g::Exact value;value.add_digits(digits.data());need(value.normalised()==digits&&std::isfinite(value.value()),"finite canonical exact value");return value;
 }
 cm::Inventory inventory(){cm::Inventory v;for(auto member:{&cm::Inventory::kinetic_J,&cm::Inventory::reference_J,&cm::Inventory::internal_J,&cm::Inventory::binding_J,&cm::Inventory::charge_e,&cm::Inventory::hydrogen,&cm::Inventory::oxygen})v.*member=exact();return v;}
 bool ended()const noexcept{return pos_==bytes_.size();}
};
inline std::shared_ptr<const cs::CommonLedger> decode(std::string_view bytes,const std::string&expected_sha256,
 std::uint64_t expected_step,double expected_time_s,std::span<const n::ArenaStateSpec>states){
 Reader r(bytes);need(r.text(128)=="common-ledger-checkpoint-v1","schema");auto v=std::make_shared<cs::CommonLedger>();
 v->accepted_step=r.u64();v->events=r.u64();v->begin_s=r.number();v->time_s=r.number();
 for(auto*p:{&v->initial_material,&v->accepted_material,&v->boundary.conductor,&v->boundary.exterior})*p=r.inventory();
 for(auto*p:{&v->initial_energy_J,&v->accepted_energy_J,&v->emitted_electrode_reference_J,&v->emitted_electrode_charge_e,
  &v->boundary.short_kinetic_J,&v->boundary.exterior_kinetic_J,&v->boundary.exterior_correction_J,&v->boundary.foil_kinetic_J,
  &v->external_net_J,&v->event_defect_J,&v->event_reference_conversion_J,&v->thermal_reference_rounding_J,
  &v->source_kinematic_rounding_J,&v->source_binding_conversion_J,&v->raw_residual_J,&v->accounted_residual_J,&v->admitted_bound_sum_J})*p=r.exact();
 std::set<std::string>actual;need(!states.empty()&&states.size()<=64,"bounded actual catalogue");for(const auto&s:states)need(actual.insert(s.id).second,"duplicate actual species");
 for(auto*m:{&v->conductor_species,&v->exterior_species}){
  need(r.u64()==actual.size(),"species account count");std::string previous;
  for(std::size_t k=0;k<actual.size();++k){auto id=r.text(4096);need(actual.contains(id)&&(k==0||id>previous),"canonical actual species order");auto count=r.exact();need(nonnegative(count),"negative species count");previous=id;m->emplace(std::move(id),std::move(count));}
 }
 v->sha256=r.text(64);need(r.ended()&&v->sha256==expected_sha256&&v->accepted_step==expected_step&&v->time_s==expected_time_s,"actual common root/clock binding or trailing data");validate(*v);return v;
}
} // namespace powerflow::common_ledger_codec_v1
