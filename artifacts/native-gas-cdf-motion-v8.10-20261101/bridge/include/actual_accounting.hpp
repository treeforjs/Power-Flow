#pragma once
#include "native_gas_targets.hpp"
#include "json_guard.hpp"
#include <map>
namespace native_gas_common {
namespace n=raikou::mfem_backend;namespace g=pf_gas_v1;namespace a=native_gas_stage;
inline void need(bool b,const char*s){if(!b)throw std::invalid_argument(std::string("native_gas_common: ")+s);}
struct Inventory {
  g::Exact kinetic,reference,packet_internal,charge,hydrogen,oxygen;
  std::array<g::Exact,3> momentum;
  std::map<std::string,g::Exact> species_count;
};
inline void carrier(Inventory&v,const n::GapParticle&p,const n::ArenaStateSpec&s) {
  need(p.weight>0 && std::isfinite(p.weight) && p.mass_kg==s.mass_kg && p.charge_C==s.charge_e*n::pic_e,"actual inventory operands");
  v.kinetic.add_product(p.weight,pfg::kinetic(s.mass_kg,p.proper_velocity_m_s.data()));
  v.reference.add_product(p.weight,s.reference_J);v.charge.add_product(p.weight,double(s.charge_e));
  v.hydrogen.add_product(p.weight,double(s.hydrogen));v.oxygen.add_product(p.weight,double(s.oxygen));
  v.species_count[s.id].add_product(p.weight,1);
  for(int k=0;k<3;++k)v.momentum[k].add_triple(p.weight,s.mass_kg,p.proper_velocity_m_s[k]);
}
inline std::map<std::string,n::ArenaStateSpec> catalogue(const n::FullTimeCandidatePopulation&v){std::map<std::string,n::ArenaStateSpec>x;for(const auto&s:v.states)need(x.emplace(s.id,s).second,"duplicate state");return x;}
inline Inventory before_inventory(std::span<const n::ArenaRateEntry>rows,const n::FullTimeCandidatePopulation&post,const n::FullTimeChemistryStockFrame&stock) {
  const auto states=catalogue(post);Inventory v;for(const auto&r:rows){const auto i=states.find(r.state);need(i!=states.end(),"before state absent");carrier(v,r.native,i->second);}
  for(const auto&p:stock.packet_internal){v.packet_internal.add_product(p.energy_J,1);}
  return v;
}
inline Inventory after_inventory(const n::FullTimeCandidatePopulation&post,const n::FullTimeChemistryStockFrame&stock) {
  const auto states=catalogue(post);Inventory v;for(const auto&p:post.packets){const auto i=states.find(p.population.state);need(i!=states.end(),"after state absent");for(const auto&r:p.population.carriers)carrier(v,r,i->second);}
  for(const auto&p:stock.packet_internal){v.packet_internal.add_product(p.energy_J,1);}
  return v;
}
struct Pending {Inventory si;g::Exact per_parcel_converted_ke_J;std::array<g::Exact,3>per_parcel_converted_p;};
inline Pending pending(const n::FullTimeChemistryStockFrame&stock,const n::FullTimeCandidatePopulation&post,const pfc::Blob&blob) {
  const auto states=catalogue(post);Pending out;
  for(const auto&r:stock.reservoirs){need(r.species<blob.n_species && r.state=="pfc-v3/"+blob.species_names[r.species] && r.rest_energy_eV==blob.view.species[r.species].mass_eV && r.reference_eV==blob.view.species[r.species].ref_eV,"stock actual blob operands");const auto i=states.find(r.state);need(i!=states.end(),"stock state catalogue");
    for(const auto&p:r.parcels){n::GapParticle row{};row.weight=p.weight;row.mass_kg=i->second.mass_kg;row.charge_C=i->second.charge_e*n::pic_e;row.position_m=p.position_m;row.proper_velocity_m_s=p.proper_velocity_m_s;carrier(out.si,row,i->second);
      auto u=pfc::v3(p.proper_velocity_m_s[0],p.proper_velocity_m_s[1],p.proper_velocity_m_s[2]);const auto mom=pfc::p_from_u(r.rest_energy_eV,u);
      out.per_parcel_converted_ke_J.add_product(p.weight,pfc::ke_from_u(r.rest_energy_eV,u)*n::pic_e);
      const double convert=n::pic_e/n::pic_c;const double vv[3]={mom.x,mom.y,mom.z};for(int k=0;k<3;++k)out.per_parcel_converted_p[k].add_product(p.weight,vv[k]*convert);
    }
  }return out;
}
struct OwnerStock {g::Exact corrected_energy,unheld_internal,reference_transfer,event_defect,reservoir_energy;std::array<g::Exact,3>corrected_momentum,unheld_internal_momentum,event_defect_p,reservoir_momentum;};
inline g::Exact global(const g::Owner&o,std::uint32_t k){g::Exact v;const auto d=o.global_digits(k,true);v.add_digits(d.data());return v;}
inline OwnerStock owner_stock(const g::Owner&o){OwnerStock v;
 v.corrected_energy=o.entry_total(pfg::kQE,true);v.corrected_energy.merge(o.entry_total(pfg::kQU,true));for(int k=0;k<3;++k)v.corrected_momentum[k]=o.entry_total(pfg::kQPX+k,true);
 for(std::uint32_t s=0;s<o.configuration().species.size();++s)if(o.configuration().species[s].charge!=0){v.corrected_energy.merge(global(o,o.layout().chg(s,0)));for(int k=0;k<3;++k)v.corrected_momentum[k].merge(global(o,o.layout().chg(s,k+1)));}
 for(int k=0;k<6;++k)v.unheld_internal.merge(global(o,pfg::G_I_UN+k));
 v.reference_transfer=global(o,pfg::G_Q);v.event_defect=global(o,pfg::G_DEF_E);v.reservoir_energy=global(o,pfg::G_RES_E);
 for(int k=0;k<3;++k){v.unheld_internal_momentum[k]=global(o,pfg::G_PINT_UN+k);v.event_defect_p[k]=global(o,pfg::G_DEF_P+k);v.reservoir_momentum[k]=global(o,pfg::G_RES_P+k);}return v;
}
inline g::Exact difference(const g::Exact&a,const g::Exact&b){auto x=a;x.merge(b,1);return x;}
inline n::ArenaJson exact_json(const g::Exact&x){const auto d=x.normalised();GuardedProofJson out;out.get()["rounded"]=x.value();out.get()["digits"]=g::digits_text(d.data());return out.release();}
struct ReferenceOperands {g::Exact actual_species_delta_J,process_Q_J,operand_conversion_J;};
inline ReferenceOperands reference_operands(const a::Batch&batch,const n::FullTimeCandidatePopulation&population,const pfc::Blob&blob) {
 const auto states=catalogue(population);ReferenceOperands out;
 auto state=[&](std::uint32_t sp)->const n::ArenaStateSpec&{need(sp<blob.n_species,"reference species index");const auto it=states.find("pfc-v3/"+blob.species_names[sp]);need(it!=states.end() && it->second.reference_J==blob.view.species[sp].ref_eV*n::pic_e && it->second.internal_J==0,"reference actual pinned catalogue");return it->second;};
 for(const auto&r:batch.recoils()) {
  const auto&e=r.ev;need(e.process<blob.n_proc,"reference process index");const auto&pr=blob.view.procs[e.process];need(pr.group<blob.n_groups,"reference process group");const auto&group=blob.view.groups[pr.group];
  need(e.target_species==group.target && e.n_out==pr.n_out && e.q_eV==pr.q_eV,"actual reference process operands");
  out.actual_species_delta_J.add_product(e.weight,state(group.projectile).reference_J,1);
  out.actual_species_delta_J.add_product(e.weight,state(group.target).reference_J,1);
  for(std::uint32_t k=0;k<e.n_out;++k){need(e.out_species[k]==pr.out_species[k] && e.out_dest[k]==pr.out_dest[k],"actual reference product operands");
   if(e.out_dest[k]==pfg::kDestGas)need(e.out_species[k]==group.target,"changed-species GAS target reference owner unsupported");
   out.actual_species_delta_J.add_product(e.weight,state(e.out_species[k]).reference_J);
  }
  out.process_Q_J.add_product(e.weight,e.q_eV*n::pic_e);
 }
 out.operand_conversion_J=difference(out.actual_species_delta_J,out.process_Q_J);return out;
}
inline n::ArenaJson closure(const Inventory&b,const Inventory&a,const Pending&pb,const Pending&pa,const OwnerStock&gb,const OwnerStock&ga,const ReferenceOperands&ref){
 // This independently compares actual stored parcels with the owner's event/release
 // mirror. No event record is substituted for the native pending population.
 need(pb.si.kinetic==gb.reservoir_energy && pa.si.kinetic==ga.reservoir_energy,"actual reservoir energy differs from gas mirror");
 for(int k=0;k<3;++k)need(pb.si.momentum[k]==gb.reservoir_momentum[k] && pa.si.momentum[k]==ga.reservoir_momentum[k],"actual reservoir momentum differs from gas mirror");
 auto delta=[&](const g::Exact&x,const g::Exact&y,const g::Exact&xr,const g::Exact&yr){auto v=difference(x,y);v.merge(xr);v.merge(yr,1);return v;};
 auto E=delta(a.kinetic,b.kinetic,pa.si.kinetic,pb.si.kinetic);E.merge(difference(a.packet_internal,b.packet_internal));
 E.merge(difference(ga.corrected_energy,gb.corrected_energy));E.merge(difference(ga.unheld_internal,gb.unheld_internal));E.merge(difference(ga.event_defect,gb.event_defect));
 auto actual_reference=delta(a.reference,b.reference,pa.si.reference,pb.si.reference);
 need(actual_reference==ref.actual_species_delta_J,"actual stored chemical reference differs from admitted stoichiometry");
 need(difference(ga.reference_transfer,gb.reference_transfer)==ref.process_Q_J,"gas reference ledger differs from pinned admitted process Q");
 E.merge(actual_reference);E.merge(ref.operand_conversion_J,1);
 auto charge=delta(a.charge,b.charge,pa.si.charge,pb.si.charge),H=delta(a.hydrogen,b.hydrogen,pa.si.hydrogen,pb.si.hydrogen),O=delta(a.oxygen,b.oxygen,pa.si.oxygen,pb.si.oxygen);
 need(charge.zero() && H.zero() && O.zero(),"actual carrier/reservoir charge or nucleus imbalance");
 GuardedProofJson Pg,rawPg,defectPg;auto&P=Pg.get();auto&rawP=rawPg.get();auto&defectP=defectPg.get();P=n::ArenaJson::array();rawP=n::ArenaJson::array();defectP=n::ArenaJson::array();bool momentum_zero=true;
 for(int k=0;k<3;++k){auto v=delta(a.momentum[k],b.momentum[k],pa.si.momentum[k],pb.si.momentum[k]);v.merge(difference(ga.corrected_momentum[k],gb.corrected_momentum[k]));v.merge(difference(ga.unheld_internal_momentum[k],gb.unheld_internal_momentum[k]));v.merge(difference(ga.event_defect_p[k],gb.event_defect_p[k]));P.push_back(nullptr);P.back()=exact_json(v);momentum_zero&=v.zero();const auto defect=difference(ga.event_defect_p[k],gb.event_defect_p[k]);defectP.push_back(nullptr);defectP.back()=exact_json(defect);v.merge(defect,1);rawP.push_back(nullptr);rawP.back()=exact_json(v);}
 auto reference=delta(a.reference,b.reference,pa.si.reference,pb.si.reference);reference.merge(difference(ga.reference_transfer,gb.reference_transfer),1);
 auto rawE=E;rawE.merge(difference(ga.event_defect,gb.event_defect),1);rawE.merge(ref.operand_conversion_J);
 GuardedProofJson result;auto&r=result.get();
 auto exact=[&](const char*key,const g::Exact&x){auto&v=r[key];v=exact_json(x);};
 exact("energy_residual_J",E);
 r["residual_definition"]="exact accounting after explicit per-event numerical-defect and reference-operand conversion ledgers; not raw physical sums";
 exact("raw_stored_energy_balance_J",rawE);
 exact("recorded_numerical_energy_defect_J",difference(ga.event_defect,gb.event_defect));
 r["raw_stored_momentum_balance"].swap(rawP);r["recorded_numerical_momentum_defect"].swap(defectP);r["momentum_residual_kg_m_s"].swap(P);
 exact("charge_e",charge);exact("H",H);exact("O",O);exact("actual_reference_minus_event_Q_J",reference);
 exact("independent_stoichiometric_reference_delta_J",ref.actual_species_delta_J);
 exact("independent_reference_operand_conversion_J",ref.operand_conversion_J);
 r["reference_operands_independently_verified"]=true;
 exact("pending_after_C3_per_parcel_conversion_minus_native_J",difference(pa.per_parcel_converted_ke_J,pa.si.kinetic));
 r["energy_exact_zero"]=E.zero();r["momentum_exact_zero"]=momentum_zero;r["production_qualified"]=false;
 r["charged_target_correction"]="ledger only; no charged target distribution feedback qualified";
 return result.release();
}
} // namespace native_gas_common
