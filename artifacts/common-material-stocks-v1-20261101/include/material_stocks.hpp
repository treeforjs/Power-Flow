#pragma once
#include "pf_gas_v1/gas_owner.hpp"
#include "chem_host.hpp"
#include "source_stage_v8.hpp"
#include "conductor_correction_guard.hpp"
#include <cstring>
// Stock readers only. Immutable snapshots are evidence, never live authority.
// Root must acquire them from the exact coupled candidate and preflight all
// owners immediately before publication. This helper does not publish.
namespace powerflow::common_material_stocks_v1 {
namespace n=::raikou::mfem_backend;namespace g=::pf_gas_v1;namespace ss=::powerflow::source_stage_v8;
inline void need(bool b,const char*s){if(!b)throw std::invalid_argument(std::string("common material stocks: ")+s);}
inline g::Exact diff(const g::Exact&a,const g::Exact&b){auto x=a;x.merge(b,1);return x;}
struct Inventory {
 g::Exact kinetic_J,reference_J,internal_J,binding_J,charge_e,hydrogen,oxygen;
 g::Exact energy()const {auto v=kinetic_J;v.merge(reference_J);v.merge(internal_J);v.merge(binding_J);return v;}
 void add(const Inventory&v){kinetic_J.merge(v.kinetic_J);reference_J.merge(v.reference_J);internal_J.merge(v.internal_J);binding_J.merge(v.binding_J);charge_e.merge(v.charge_e);hydrogen.merge(v.hydrogen);oxygen.merge(v.oxygen);}
};
inline std::map<std::string,n::ArenaStateSpec> species_catalogue(const n::FullTimeCandidatePopulation&population){std::map<std::string,n::ArenaStateSpec> out;for(const auto&sp:population.states)need(out.emplace(sp.id,sp).second,"duplicate species");return out;}
inline void add_particle(Inventory&out,const n::GapParticle&p,const n::ArenaStateSpec&s){need(p.weight>0&&std::isfinite(p.weight)&&p.mass_kg==s.mass_kg&&p.charge_C==s.charge_e*n::pic_e&&s.internal_J==0,"actual particle operands and explicit packet-owned internal energy");out.kinetic_J.add_product(p.weight,pfg::kinetic(s.mass_kg,p.proper_velocity_m_s.data()));out.reference_J.add_product(p.weight,s.reference_J);out.charge_e.add_product(p.weight,double(s.charge_e));out.hydrogen.add_product(p.weight,double(s.hydrogen));out.oxygen.add_product(p.weight,double(s.oxygen));}
inline Inventory pending_inventory(const n::FullTimeChemistryStockFrame&stock,const n::FullTimeCandidatePopulation&population,const pfc::Blob&blob){const auto states=species_catalogue(population);Inventory out;for(const auto&r:stock.reservoirs){need(r.species<blob.n_species&&r.state=="pfc-v3/"+blob.species_names[r.species]&&r.rest_energy_eV==blob.view.species[r.species].mass_eV&&r.reference_eV==blob.view.species[r.species].ref_eV,"pending actual blob operands");auto it=states.find(r.state);need(it!=states.end(),"pending species catalogue");for(const auto&p:r.parcels){n::GapParticle row{};row.weight=p.weight;row.mass_kg=it->second.mass_kg;row.charge_C=it->second.charge_e*n::pic_e;row.proper_velocity_m_s=p.proper_velocity_m_s;add_particle(out,row,it->second);}}return out;}
inline Inventory native_before(std::span<const n::ArenaRateEntry>rows,const n::FullTimeCandidatePopulation&population,const n::FullTimeChemistryStockFrame&stock,const pfc::Blob&blob){const auto states=species_catalogue(population);Inventory out=pending_inventory(stock,population,blob);std::set<std::uint64_t>ids;for(const auto&r:rows){auto it=states.find(r.state);need(it!=states.end()&&ids.insert(r.native.id).second,"before particle species/identity");add_particle(out,r.native,it->second);}for(const auto&p:stock.packet_internal)out.internal_J.add_product(p.energy_J,1);return out;}
inline Inventory native_after(const n::FullTimeCandidatePopulation&population,const n::FullTimeChemistryStockFrame&stock,const pfc::Blob&blob){std::vector<n::ArenaRateEntry>rows;rows.reserve(population.carriers);for(const auto&p:population.packets)for(const auto&c:p.population.carriers)rows.push_back({p.population.packet_id,p.population.state,c});need(rows.size()==population.carriers,"actual population cardinality");return native_before(rows,population,stock,blob);}
inline g::Exact global(const g::Owner&o,std::uint32_t k,bool candidate){g::Exact v;const auto d=o.global_digits(k,candidate);v.add_digits(d.data());return v;}
struct GasStock {g::Exact correction_J,unheld_internal_J,pending_mirror_J,numerical_defect_J,process_Q_J;g::Exact energy()const{auto v=correction_J;v.merge(unheld_internal_J);return v;}};
inline GasStock gas_stock(const g::Owner&o,bool candidate){
 GasStock v;v.correction_J=o.entry_total(pfg::kQE,candidate);v.correction_J.merge(o.entry_total(pfg::kQU,candidate));
 for(std::uint32_t sp=0;sp<o.configuration().species.size();++sp)if(o.configuration().species[sp].charge!=0)v.correction_J.merge(global(o,o.layout().chg(sp,0),candidate));
 for(std::uint32_t k=0;k<6;++k)v.unheld_internal_J.merge(global(o,pfg::G_I_UN+k,candidate));
 v.pending_mirror_J=global(o,pfg::G_RES_E,candidate);v.numerical_defect_J=global(o,pfg::G_DEF_E,candidate);v.process_Q_J=global(o,pfg::G_Q,candidate);return v;
}
// This is the finite-stock potential of the declared constant/Temkin law.
// Linear and quadratic coefficients are explicit binary64 SI operands; the
// count polynomial is accumulated exactly, including low count digits.
// The difference from the actually released rounded binding operand is exposed
// separately by SourceAssessment, never silently removed from raw energy.
template<class Count>inline g::Exact binding_stock(const Count&count,const ss::eh::Patch&p,double area){
 need(p.desorption_enabled&&std::isfinite(p.desorption_eV)&&p.desorption_eV>=0,"declared adsorption binding");
 const auto terms=count.positive_carriers();need(terms.size()<=256,"finite exact count decomposition cap");const double D=p.desorption_eV*n::pic_e;need(std::isfinite(D),"binding coefficient");
 g::Exact result;for(double w:terms)result.add_product(w,D,1);
 if(p.coverage_temkin){need(std::isfinite(area)&&area>0&&std::isfinite(p.monolayer_m2)&&p.monolayer_m2>0&&std::isfinite(p.binding_alpha)&&p.binding_alpha>=0&&p.binding_alpha<=1,"Temkin binding operands");const double capacity=area*p.monolayer_m2;need(std::isfinite(capacity)&&capacity>0,"finite monolayer capacity");const double coefficient=D*p.binding_alpha/(2*capacity);need(std::isfinite(coefficient),"Temkin coefficient");for(double a:terms)for(double b:terms)result.add_triple(coefficient,a,b);}
 return result;
}
inline Inventory patch_stock(const ss::s7::Count&count,const ss::eh::Patch&p,double area,const n::ArenaStateSpec&species){
 need(species.id==p.state_id&&species.charge_e==0&&species.internal_J==0&&species.mass_kg==p.mass_kg,"source stock species operands");Inventory v;
 for(double w:count.positive_carriers()){v.reference_J.add_product(w,species.reference_J);v.hydrogen.add_product(w,double(species.hydrogen));v.oxygen.add_product(w,double(species.oxygen));}
 v.binding_J=binding_stock(count,p,area);return v;
}
struct SourceAssessment {
 Inventory before,after;
 g::Exact emitted_electrode_reference_J,emitted_electrode_charge_e,source_wall_J,declared_neutral_binding_J;
 g::Exact neutral_kinematic_rounding_J,electron_kinematic_rounding_J,binding_operand_conversion_J;
 std::size_t patches=0,neutral_births=0,electron_births=0;
};
inline SourceAssessment source_stocks(const ss::OwnerSet&owners,const ss::PreparedSourceStage&stage,const n::FullTimeCandidatePopulation&catalogue){
 need(stage.preflight_sources(owners)&&stage.source_trials().size()==owners.ports.size(),"actual source candidate preflight");const auto states=species_catalogue(catalogue);SourceAssessment out;
 for(std::size_t i=0;i<owners.ports.size();++i){const auto&port=owners.ports[i];const auto&owner=*port.source();const auto&trial=stage.source_trials()[i];
  // Active-layer source storage requires the separate held-layer owner reader.
  // Current surface-birth stage explicitly excludes these paths.
  need(owner.active_layers().empty()&&trial.source_candidates().empty()&&trial.source_handoffs().empty(),"held/recollected source layers require an actual before/after layer stock reader");
  for(std::size_t j=0;j<owner.source_patches();++j){const auto&p=owner.source_patch(j);if(!p.desorption_enabled){need(owner.inventory(j)==trial.inventory(j),"disabled source changed");continue;}auto it=states.find(p.state_id);need(it!=states.end(),"source species absent from actual catalogue");const double area=owner.patch_identity(j).physical_area_m2;out.before.add(patch_stock(owner.inventory(j),p,area,it->second));out.after.add(patch_stock(trial.inventory(j),p,area,it->second));++out.patches;}
  for(const auto&p:trial.packets()){auto it=states.find(p.state_id);need(it!=states.end()&&p.mass_kg==it->second.mass_kg&&p.weight>0,"actual neutral packet operands");out.declared_neutral_binding_J.add_product(p.adsorbed_reference_change_J,1);out.neutral_kinematic_rounding_J.add_product(p.kinetic_J,1);out.neutral_kinematic_rounding_J.add_product(p.weight,pfg::kinetic(p.mass_kg,p.proper_velocity_m_s.data()),1);++out.neutral_births;}
  for(const auto&b:trial.surface_births()){const auto&packet=b.phase().packet.particle;const auto&ledger=packet.ledger;need(ledger.image_reference_change_J==0&&ledger.incoming_K_J==0&&packet.weight>0,"native handoff must restore image tail and have no incoming particle");out.emitted_electrode_reference_J.add_product(ledger.electron_reference_change_J,1);out.emitted_electrode_charge_e.add_product(ledger.physical_electron_number,1);out.electron_kinematic_rounding_J.add_product(ledger.outgoing_K_J,1);out.electron_kinematic_rounding_J.add_product(packet.weight,pfg::kinetic(packet.electron_mass_kg,packet.proper_velocity_m_s.data()),1);++out.electron_births;}
  for(const auto&load:trial.loads())out.source_wall_J.add_product(hall_corrective_v3::transfer_energy_J(load),1);
 }
 out.binding_operand_conversion_J=diff(diff(out.after.binding_J,out.before.binding_J),out.declared_neutral_binding_J);return out;
}
struct ChemistryOperands {g::Exact species_reference_delta_J,process_Q_J,numerical_defect_J,reference_conversion_J;std::uint64_t events=0;};
inline ChemistryOperands admitted_event_operands(std::span<const std::shared_ptr<const n::FullTimeChemistryStageExport>>exports,const n::FullTimeCandidatePopulation&catalogue,const pfc::Blob&blob,std::uint64_t cap){
 const auto states=species_catalogue(catalogue);ChemistryOperands out;std::set<std::pair<std::uint64_t,std::uint64_t>>seen;
 auto state=[&](std::uint32_t sp)->const n::ArenaStateSpec&{need(sp<blob.n_species,"event species");auto it=states.find("pfc-v3/"+blob.species_names[sp]);need(it!=states.end()&&it->second.reference_J==blob.view.species[sp].ref_eV*n::pic_e,"event actual species reference");return it->second;};
 for(const auto&x:exports){need(bool(x)&&seen.emplace(x->accepted_step,x->owner_stage_sequence).second,"missing or duplicate actual stage export");need(x->events.size()<=cap-out.events,"bounded event operands");
  for(const auto&record:x->events){need(record.event_bytes.size()==sizeof(pfc::EventRec),"event ABI");pfc::EventRec e;std::memcpy(&e,record.event_bytes.data(),sizeof e);need(e.process<blob.n_proc&&e.weight==record.weight&&e.weight>0,"admitted event operands");const auto&pr=blob.view.procs[e.process];need(pr.group<blob.n_groups&&e.n_out==pr.n_out,"event process");const auto&group=blob.view.groups[pr.group];need(group.target==e.target_species&&e.q_eV==pr.q_eV,"event target/reference");const auto&sp=state(group.projectile);const auto&st=state(group.target);
   out.species_reference_delta_J.add_product(e.weight,sp.reference_J,1);out.species_reference_delta_J.add_product(e.weight,st.reference_J,1);out.numerical_defect_J.add_product(e.weight,pfg::kinetic(sp.mass_kg,e.u_in));out.numerical_defect_J.add_product(e.weight,pfg::kinetic(st.mass_kg,e.ut));
   for(std::uint32_t j=0;j<e.n_out;++j){need(e.out_species[j]==pr.out_species[j]&&e.out_dest[j]==pr.out_dest[j],"event product identity");const auto&so=state(e.out_species[j]);out.species_reference_delta_J.add_product(e.weight,so.reference_J);out.numerical_defect_J.add_product(e.weight,pfg::kinetic(so.mass_kg,e.out_u[j]),1);}
   const double Q=e.q_eV*n::pic_e,I=e.e_internal*n::pic_e;out.process_Q_J.add_product(e.weight,Q);out.numerical_defect_J.add_product(e.weight,Q,1);out.numerical_defect_J.add_product(e.weight,I,1);++out.events;
  }
 }
 out.reference_conversion_J=diff(out.species_reference_delta_J,out.process_Q_J);return out;
}
struct BoundaryStock {Inventory conductor,exterior;g::Exact short_kinetic_J,exterior_kinetic_J,exterior_correction_J,foil_kinetic_J;};
inline BoundaryStock boundary_stocks(const n::FullTimeSinkLedger*ledger,const g::TransportReceipt&transport,const n::FullTimeCandidatePopulation&catalogue,const g::Configuration&configuration,std::uint64_t expected_step){
 ::powerflow::conductor_correction_guard_v1::require_no_unrouted_conductor_correction(transport,expected_step,configuration);
 const auto states=species_catalogue(catalogue);BoundaryStock out;std::map<std::uint64_t,const n::FullTimeSinkExposure*>sink_ids;
 if(ledger){need(ledger->accepted_step==expected_step&&ledger->exposed.size()==ledger->records.size(),"complete actual sink exposure");std::size_t record_index=0;for(const auto&x:ledger->exposed){const auto&canonical=ledger->records.at(record_index++);need(x.record.native_id==canonical.native_id&&x.record.kind==canonical.kind&&x.record.surface==canonical.surface&&x.record.incident_proper_velocity_m_s==canonical.incident_proper_velocity_m_s,"sink exposure differs from canonical actual record");need(sink_ids.emplace(x.record.native_id,&x).second,"duplicate sink");auto it=states.find(x.state);need(it!=states.end()&&x.mass_kg==it->second.mass_kg&&x.charge_C==it->second.charge_e*n::pic_e&&x.weight>0,"sink species operands");Inventory one;one.reference_J.add_product(x.weight,it->second.reference_J);one.internal_J.add_product(x.internal_energy_J,1);one.charge_e.add_product(x.weight,double(it->second.charge_e));one.hydrogen.add_product(x.weight,double(it->second.hydrogen));one.oxygen.add_product(x.weight,double(it->second.oxygen));g::Exact K;K.add_product(x.weight,pfg::kinetic(x.mass_kg,x.record.incident_proper_velocity_m_s.data()));
  const auto surface=x.record.surface;const bool foil=surface==n::FullTimeSinkSurface::lower_foil||surface==n::FullTimeSinkSurface::upper_foil,shorted=surface==n::FullTimeSinkSurface::short_bridge,exterior=surface==n::FullTimeSinkSurface::outer_xy||surface==n::FullTimeSinkSurface::axial;
  need((foil||shorted||exterior)&&((foil||shorted)==(x.record.kind==n::FullTimeSinkKind::collected_at_conductor)),"sink destination partition");
  if(foil){need(x.internal_energy_J==0,"foil packet internal energy requires qualified thermal routing");out.foil_kinetic_J.merge(K);}else if(shorted)out.short_kinetic_J.merge(K);else out.exterior_kinetic_J.merge(K);
  (exterior?out.exterior:out.conductor).add(one);
 }}
 for(const auto&exit:transport.exits){need(sink_ids.count(exit.native_id),"gas correction exit lacks actual native sink");const auto&actual=*sink_ids.at(exit.native_id);need(configuration.species[exit.species].state==actual.state&&((exit.destination==g::ExitDestination::exterior)==(actual.record.kind==n::FullTimeSinkKind::escaped_exterior)),"gas correction sink species/destination");if(exit.destination==g::ExitDestination::exterior){out.exterior_correction_J.add_digits(exit.digits[0].data());out.exterior_correction_J.add_digits(exit.digits[4].data());}}
 return out;
}
} // namespace powerflow::common_material_stocks_v1
