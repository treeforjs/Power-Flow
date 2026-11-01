#pragma once
#include "raikou/mfem_pic_full_time.h"
#include "pf_gas_v1/gas_owner.hpp"
#include "chem_host.hpp"
#include "chem_v3_native.hpp"
#include <bit>
#include <cstring>
#include <limits>
#include <set>

// Development-only records adapter. It never publishes either owner, injects
// targets, enables chemistry, or establishes whole-system energy closure.
namespace native_gas_stage {
namespace n = raikou::mfem_backend;
namespace g = pf_gas_v1;
inline void need(bool value,const char *message) {
  if(!value) throw std::invalid_argument(std::string("native_gas_stage: ")+message);
}
inline bool same(double a,double b) {return std::bit_cast<std::uint64_t>(a)==std::bit_cast<std::uint64_t>(b);}
static_assert(std::endian::native==std::endian::little);
static_assert(sizeof(pfg::ChemEvent)==sizeof(pfc::EventRec));
#define NGS_FIELD(field) static_assert(offsetof(pfg::ChemEvent,field)==offsetof(pfc::EventRec,field))
NGS_FIELD(carrier);NGS_FIELD(candidate);NGS_FIELD(process);NGS_FIELD(cell);NGS_FIELD(target_species);
NGS_FIELD(n_out);NGS_FIELD(cover_class);NGS_FIELD(internal_sink);NGS_FIELD(flags);NGS_FIELD(weight);
NGS_FIELD(time_s);NGS_FIELD(u_in);NGS_FIELD(ut);NGS_FIELD(E_tf);NGS_FIELD(ke_in_p);NGS_FIELD(ke_in_t);
NGS_FIELD(q_eV);NGS_FIELD(e_internal);NGS_FIELD(out_species);NGS_FIELD(out_dest);NGS_FIELD(out_u);NGS_FIELD(ke_out);
#undef NGS_FIELD

inline std::uint32_t internal_carrier(std::uint32_t kind,std::uint32_t count) {
  switch(kind) {
    case 1:return pfg::kNone;
    case 2:need(count>=2,"excitation product count");return 1;
    case 3:need(count>=2,"dissociation product count");return count-1;
    case 4:need(count>=3,"ionization product count");return 2;
    case 5:need(count>=1,"capture product count");return 0;
    case 6:return count>2?count-1:pfg::kNone;
    case 7:need(count>=3,"ion-impact ionization product count");return 2;
    case 8:need(count>=2,"breakup product count");return count-2;
    default:throw std::invalid_argument("native_gas_stage: unsupported process kind");
  }
}

// Verify actual blob-derived process operands rather than accepting a registry
// hash alongside a different process table. Extra non-chemistry gas species may
// follow the registry, but cannot replace/reorder its indices.
inline void configuration(const g::Configuration &cfg,const n::FullTimeChemistryRegistry &registry,
                          const pfc::Blob &blob) {
  need(blob.bytes>0 && blob.bytes<=blob.storage.size()*sizeof(std::uint64_t) &&
       pfc::sha256_hex(blob.base(),blob.bytes)==blob.sha256,"actual blob bytes");
  std::size_t ns=0,ng=0,np=0;
  need(blob.view.species==blob.at<pfc::SpeciesRec>("species",&ns) &&
       blob.view.groups==blob.at<pfc::GroupRec>("groups",&ng) &&
       blob.view.procs==blob.at<pfc::ProcessRec>("procs",&np) && ns==blob.n_species &&
       ng==blob.n_groups && np==blob.n_proc && blob.species_names.size()==ns,"blob table views");
  need(cfg.chemistry_blob_sha256==blob.sha256 && registry.blob_sha256==blob.sha256 &&
       cfg.chemistry_registry_sha256==registry.registry_sha256,"blob/registry identity");
  need(registry.species.size()==blob.n_species && cfg.species.size()>=registry.species.size() &&
       cfg.processes.size()==blob.n_proc,"registry/process cardinality");
  std::vector<bool> targets(blob.n_species,false);
  for(std::uint32_t k=0;k<blob.n_groups;++k) {
    need(blob.view.groups[k].target<blob.n_species && blob.view.groups[k].projectile<blob.n_species,"group species");
    targets[blob.view.groups[k].target]=true;
  }
  for(std::size_t k=0;k<registry.species.size();++k) {
    const auto &s=registry.species[k];const auto &c=cfg.species[k];
    const double mass=(blob.view.species[k].mass_eV*n::pic_e)/(n::pic_c*n::pic_c);
    need(s.state=="pfc-v3/"+blob.species_names[k] && c.state==s.state && same(c.mass_kg,s.mass_kg) &&
         same(s.mass_kg,mass) && c.charge==s.charge_e && c.charge==blob.view.species[k].charge &&
         c.chemistry_target==targets[k] && s.target==targets[k],"species operands");
  }
  for(std::uint32_t k=0;k<blob.n_proc;++k) {
    const auto &p=blob.view.procs[k];const auto &c=cfg.processes[k];
    need(p.group<blob.n_groups && p.n_out>=1 && p.n_out<=pfg::kMaxOut,"process layout");
    need(c.projectile==blob.view.groups[p.group].projectile && c.n_out==p.n_out && c.kind==p.kind &&
         c.internal_carrier==internal_carrier(p.kind,p.n_out),"process operands");
  }
}

// Construct once per immutable native/gas configuration. The native grid and
// solid-cell mask are reused, not rebuilt on every collision stage.
class Binding {
  g::Configuration gas_;
  n::FullTimeChemistryRegistry registry_;
  n::chem8::Grid grid_;
  std::string gas_identity_;
  std::map<std::string,std::uint32_t> species_;
public:
  Binding(const g::Owner &gas,const n::FullTimeChemistryRegistry &registry,const pfc::Blob &blob,
          const n::GapPicOptions &options) : gas_(gas.configuration()),registry_(registry) {
    configuration(gas_,registry_,blob);
    grid_=n::chem8::make_grid(options,gas_.grid.n);
    need(grid_.sha256==registry_.grid_sha256 && grid_.low==gas_.grid.low && grid_.high==gas_.grid.high &&
         grid_.n==gas_.grid.n && same(grid_.volume,gas_.grid.volume),"native/gas collision geometry");
    gas_identity_=gas.configuration_sha256();
    for(std::uint32_t k=0;k<gas_.species.size();++k)
      need(species_.emplace(gas_.species[k].state,k).second,"duplicate gas species");
  }
  const g::Configuration &gas() const noexcept {return gas_;}
  const n::FullTimeChemistryRegistry &registry() const noexcept {return registry_;}
  const n::chem8::Grid &grid() const noexcept {return grid_;}
  const std::string &gas_identity() const noexcept {return gas_identity_;}
  const std::map<std::string,std::uint32_t> &species() const noexcept {return species_;}
};

class Batch {
  friend Batch prepare(const n::CollisionScope &,const n::PreparedArena &,const g::Owner &,
                       const Binding &,std::size_t);
  std::shared_ptr<const n::FullTimeChemistryStageExport> export_;
  std::string gas_configuration_,gas_accepted_,gas_candidate_;
  g::Owner::CandidateBinding gas_attempt_;
  std::vector<pfg::RecoilRecord> recoil_;
  std::vector<pfg::FlowRecord> flows_;
  Batch()=default;
public:
  Batch(const Batch &)=delete;
  Batch &operator=(const Batch &)=delete;
  Batch(Batch &&) noexcept=default;
  Batch &operator=(Batch &&) noexcept=default;
  std::span<const pfg::RecoilRecord> recoils() const noexcept {return recoil_;}
  std::span<const pfg::FlowRecord> flows() const noexcept {return flows_;}
  const n::FullTimeChemistryStageExport &source() const {need(bool(export_),"moved batch");return *export_;}
  // This is a preparation-time check, not a no-fail publication operation.
  // A digest alone never authorizes an abandoned deterministic retry.
  void validate_current(const n::CollisionScope &scope,const n::PreparedArena &token,const g::Owner &gas) const {
    need(export_ && scope.can_inspect(token) && gas.matches_candidate(gas_attempt_,export_->accepted_step),"dead native/gas candidate");
    // The pinned native implementation creates and validates one immutable
    // deep snapshot per attempt and retains no mutable alias. Rechecking its
    // complete event/count hashes on every preflight adds work but no new
    // information. Authority still comes from both live serial-owner tokens.
    const auto stage=token.candidate_chemistry_stage();
    need(stage && stage->admitted_export==export_,"different native attempt");
    need(gas.configuration_sha256()==gas_configuration_ && gas.accepted_sha256()==gas_accepted_ &&
         gas.candidate_sha256()==gas_candidate_,"changed gas base/candidate");
    need(scope.can_inspect(token),"native candidate changed during validation");
  }
};

inline Batch prepare(const n::CollisionScope &scope,const n::PreparedArena &token,const g::Owner &gas,
                     const Binding &binding,std::size_t cap) {
  need(cap>0 && scope.can_inspect(token) && gas.step_open(),"live candidate and positive capacity required");
  need(gas.configuration_sha256()==binding.gas_identity(),"gas configuration changed");
  const auto &cfg=binding.gas();const auto &registry=binding.registry();
  Batch b;b.export_=token.candidate_admitted_export(true);need(bool(b.export_),"missing admitted export");
  // candidate_admitted_export(true) above already performs full strict validation.
  const auto &x=*b.export_;
  const auto &exposure=scope.exposure();const auto &common=exposure.common;
  need(x.configuration_sha256==registry.configuration_sha256 && x.registry_sha256==registry.registry_sha256 &&
       x.accepted_native_sha256==common.accepted_arena_sha256 && x.common_clock_sha256==common.common_clock_sha256 &&
       x.accepted_step==common.accepted_step && same(x.begin_s,common.begin_s) && same(x.end_s,common.end_s) &&
       same(x.step_duration_s,common.dt_s) && same(x.collision_duration_s,exposure.collision_duration_s) &&
       same(x.frozen_physical_stage_s,exposure.frozen_physical_stage_s) &&
       x.frozen_field_sha256==exposure.frozen_field_sha256 && x.owner_stage_sequence==exposure.owner_stage_sequence &&
       x.stage==(exposure.slot==n::FullTimeSlot::before_transport?1u:2u),"native exposure/registry binding");
  b.gas_attempt_=gas.candidate_binding();
  need(gas.matches_candidate(b.gas_attempt_,x.accepted_step),"gas/native step differs");
  b.gas_configuration_=gas.configuration_sha256();b.gas_accepted_=gas.accepted_sha256();b.gas_candidate_=gas.candidate_sha256();
  need(x.events.size()<=cap,"event capacity");
  std::size_t remaining=cap-x.events.size();
  auto account=[&](std::size_t count){need(count<=remaining,"total record capacity");remaining-=count;};
  for(const auto &d:x.target_debits){account(d.before.carriers.size());account(d.remainder.carriers.size());}
  for(const auto &r:x.releases)account(r.carriers.size());
  b.recoil_.reserve(x.events.size());b.flows_.reserve(cap-x.events.size()-remaining);
  const std::uint64_t cells=std::uint64_t(cfg.grid.n[0])*cfg.grid.n[1]*cfg.grid.n[2];
  for(const auto &event:x.events) {
    pfg::RecoilRecord r{};r.accepted_step=x.accepted_step;r.stage=std::uint32_t(x.stage);r.grid_cell=event.atlas_cell;
    need(event.event_bytes.size()==sizeof(r.ev) && r.grid_cell<cells && binding.grid().vacuum[r.grid_cell],"event ABI/cell");
    std::memcpy(&r.ev,event.event_bytes.data(),sizeof(r.ev));
    need(r.ev.process<cfg.processes.size() && r.ev.n_out>=1 && r.ev.n_out<=pfg::kMaxOut &&
         r.ev.target_species<registry.species.size() && same(r.ev.weight,event.weight),"event process/species");
    bool gas_return=false;
    for(std::uint32_t k=0;k<r.ev.n_out;++k)gas_return|=r.ev.out_dest[k]==pfg::kDestGas;
    need(gas_return?!event.consumes_tracked_target:event.consumes_tracked_target,
         "untracked consumed target or GAS-return debit");
    b.recoil_.push_back(r);
  }
  const auto &species=binding.species();
  auto packet=[&](const n::FullTimeChemistryPacketSnapshot &p,std::uint32_t sp,std::uint32_t kind,
                  bool negative,std::uint32_t expected_cell) {
    need(p.carried_internal_energy_J==0,"packet internal-energy attribution not supported");
    if(p.packet_id.empty()) {need(p.count.zero() && p.carriers.empty(),"absent packet count");return;}
    need(sp<cfg.species.size() && p.state==cfg.species[sp].state,"packet species mapping");
    n::ArenaCount total;
    for(const auto &particle:p.carriers) {
      need(b.flows_.size()<cap-b.recoil_.size() && std::isfinite(particle.weight) && particle.weight>0 &&
           same(particle.mass_kg,cfg.species[sp].mass_kg) && particle.charge_C==cfg.species[sp].charge*n::pic_e &&
           pfg::finite3(particle.position_m.data()) && pfg::finite3(particle.proper_velocity_m_s.data()),"carrier operands/capacity");
      const auto cell=binding.grid().locate(particle.position_m.data());
      need(cell!=pfg::kNone && (expected_cell==pfg::kNone || cell==expected_cell),"carrier atlas cell");
      pfg::FlowRecord f{};f.kind=kind;f.negative=negative?1u:0u;f.grid_cell=cell;f.species=sp;
      f.native_id=particle.id;f.weight=particle.weight;
      std::copy(particle.proper_velocity_m_s.begin(),particle.proper_velocity_m_s.end(),f.u);
      b.flows_.push_back(f);total=total.credited(particle.weight);
    }
    need(total==p.count,"carrier count projection differs from exact packet count");
  };
  for(const auto &debit:x.target_debits) {
    need(!debit.also_projectile && !debit.internal_energy_deferred,"ambiguous target debit");
    packet(debit.before,debit.species,pfg::kFlowDebit,true,debit.atlas_cell);
    packet(debit.remainder,debit.species,pfg::kFlowDebit,false,debit.atlas_cell);
  }
  for(const auto &release:x.releases) {
    need(x.stage==2 && species.count(release.state),"release stage/species");
    packet(release,species.at(release.state),pfg::kFlowRelease,false,pfg::kNone);
  }
  b.validate_current(scope,token,gas);
  return b;
}
} // namespace native_gas_stage
