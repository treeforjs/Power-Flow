#pragma once
#include "actual_accounting.hpp"
#include "raikou/mfem_pic_coupled_targets.h"
#include "pf_coupled/gas_step_scope.hpp"

namespace native_gas_common {
// Strict target-owner bookkeeping, bounded CPU reference. Unknown-rate/network
// policy, native motion and materials retain their independent admission gates.
inline std::string coupled_contract(const n::FullTimeChemistryConfiguration&configuration,
                                   const n::GapPicOptions&options,const g::Owner&gas) {
 auto base=configuration;base.host_target_feedback=false;base.host_target_contract_sha256.clear();
 base.coupled_target_contract_sha256.clear();
 const auto reg=n::full_time_chemistry_registry(base,options);
 return native_gas_targets::sha(n::ArenaJson{{"schema","native-gas-provisional-owner-cpu-v1"},
  {"native_source",n::PicArena::source_identity()},{"base_configuration",reg.configuration_sha256},
  {"gas_configuration",gas.configuration_sha256()},
  {"accounting","actual-stocks;water-nuclei;explicit-numerical-defect-and-reference-conversion"},
  {"motion_scope","same-carriers-species-weights-and-cells;no-source-or-sink"}}.dump());
}
class CoupledTrial final {
public:
 class Continuation final {
   friend class CoupledTrial;
   std::shared_ptr<const n::FullTimeCoupledAccess::Continuation> native_;
   std::shared_ptr<g::Owner> gas_;
   std::string gas_after_,common_receipt_;
   std::uint64_t next_step_=0;
   Continuation()=default;
 public:
   Continuation(const Continuation&)=delete;Continuation&operator=(const Continuation&)=delete;
   const std::string&common_receipt_sha256()const noexcept{return common_receipt_;}
   std::uint64_t next_step()const noexcept{return next_step_;}
 };
private:
 struct StageWitness {
  std::shared_ptr<g::Owner> gas;g::Owner::CandidateBinding binding;std::uint64_t step=0;
  static bool valid(const void*p) noexcept {const auto &w=*static_cast<const StageWitness*>(p);return w.gas && w.gas->matches_candidate(w.binding,w.step);}
 };
 struct FinalWitness {
  std::shared_ptr<g::Owner> gas;g::Owner::PreparedCommit prepared;
  static bool valid(const void*p) noexcept {const auto &w=*static_cast<const FinalWitness*>(p);return w.gas && w.gas->can_publish(w.prepared);}
 };
 struct State {
  // Exact native outer and gas attempt ownership survive wrapper moves.
  n::PrivateFullTimeStep outer;std::shared_ptr<g::Owner> gas;pf_coupled::GasStepScope rollback;
  n::FullTimeChemistryConfiguration cfg;n::GapPicOptions options;n::FullTimeChemistryRegistry registry;
  pfc::Blob blob;a::Binding binding;
  std::string initial_gas,initial_native,contract,parent_receipt;
  std::uint64_t step=0;std::size_t cap=0;unsigned completed=0;bool failed=false,published=false;
  std::shared_ptr<const n::FullTimeCandidatePopulation> first_population;
  n::ArenaJson stages=n::ArenaJson::array(),receipt;
  std::unique_ptr<n::PreparedArena> final_native;std::shared_ptr<FinalWitness> final_gas;
  std::shared_ptr<const Continuation> continuation;
  g::Owner::CandidateBinding initial_binding;std::shared_ptr<StageWitness> last_witness;
  ~State()noexcept{drain_proof_json(stages);drain_proof_json(receipt);}
  State(n::PrivateFullTimeStep o,std::shared_ptr<g::Owner> gowner,const n::FullTimeChemistryConfiguration&c,
    const n::GapPicOptions&opts,std::uint64_t accepted_step,std::size_t maximum_carriers)
   :outer(std::move(o)),gas(std::move(gowner)),cfg(c),options(opts),registry(n::full_time_chemistry_registry(c,opts)),
    blob(pfc::load_blob(c.blob_path,c.blob_sha256)),binding(*(need(bool(gas),"missing gas owner"),gas),registry,blob,opts),step(accepted_step),cap(maximum_carriers) {
    need(cap>0 && cap<=1000000,"coupled carrier capacity");
    need(!cfg.diagnostic_field_target_approximation && cfg.host_target_feedback &&
      !cfg.coupled_target_contract_sha256.empty() && cfg.execution==n::FullTimeChemistryExecution::host &&
      gas->configuration().execution==g::Execution::host && gas->configuration().transport==g::Transport::none,
      "bounded coupled CPU strict mode/transport configuration");
    contract=coupled_contract(cfg,options,*gas);need(contract==cfg.coupled_target_contract_sha256,"coupled contract differs");
    initial_gas=gas->accepted_sha256();rollback=pf_coupled::GasStepScope::begin(gas,step);initial_binding=gas->candidate_binding();
  }
 };
 std::unique_ptr<State>s_;
 void active() const {need(s_ && !s_->failed && !s_->published && !s_->final_native && s_->rollback.active(),"closed/failed coupled trial");}
 void fail() noexcept {if(s_){s_->failed=true;s_->rollback.reset();}}
 void water_scope(std::span<const n::ArenaRateEntry>rows) const {
  for(const auto&r:rows){auto i=s_->binding.species().find(r.state);need(i!=s_->binding.species().end() && i->second<s_->registry.species.size(),"unregistered actual state");
   const auto&sp=s_->registry.species[i->second];for(std::size_t k=2;k<sp.nuclei.size();++k)need(sp.nuclei[k]==0,"independent nucleus accounting currently supports H/O only");}
 }
 void motion_scope(std::span<const n::ArenaRateEntry>rows) const {
  need(s_->first_population && rows.size()==s_->first_population->carriers,"bounded motion carrier count changed");
  std::map<std::uint64_t,std::pair<std::string,n::GapParticle>>before;
  for(const auto&p:s_->first_population->packets)for(const auto&r:p.population.carriers)
    need(before.emplace(r.id,std::make_pair(p.population.state,r)).second,"duplicate stage1 carrier");
  const auto &cfg=s_->gas->configuration();
  auto cell=[&](const n::GapParticle&r){std::uint32_t ix[3];for(int d=0;d<3;++d){const double x=(r.position_m[d]-cfg.grid.low[d])/(cfg.grid.high[d]-cfg.grid.low[d])*cfg.grid.n[d];need(std::isfinite(x)&&x>=0 && x<cfg.grid.n[d],"motion outside gas atlas");ix[d]=std::uint32_t(x);}return std::uint64_t(ix[0])+std::uint64_t(cfg.grid.n[0])*(ix[1]+std::uint64_t(cfg.grid.n[1])*ix[2]);};
  for(const auto&r:rows){auto i=before.find(r.native.id);need(i!=before.end() && i->second.first==r.state &&
    a::same(i->second.second.weight,r.native.weight) && cell(i->second.second)==cell(r.native),
    "gas motion adapter required for new/sunk/changed/cell-crossing carriers");before.erase(i);}
  need(before.empty(),"missing moved carriers");
 }
public:
 CoupledTrial(n::PrivateFullTimeStep outer,std::shared_ptr<g::Owner>gas,const n::FullTimeChemistryConfiguration&cfg,
              const n::GapPicOptions&options,std::uint64_t step,std::size_t maximum_carriers,
              std::shared_ptr<const Continuation> previous={})
  :s_(std::make_unique<State>(std::move(outer),std::move(gas),cfg,options,step,maximum_carriers)) {
  if(previous)need(previous->gas_==s_->gas && previous->gas_after_==s_->initial_gas &&
                    previous->next_step_==step,"foreign/stale gas accepted-pair continuation");
  else {
    g::Owner pristine(s_->gas->configuration());
    need(step==0 && pristine.accepted_sha256()==s_->initial_gas,"unattested non-genesis gas pair");
  }
  s_->parent_receipt=previous?previous->common_receipt_:std::string{};
  n::FullTimeCoupledAccess::attach(s_->outer,s_->contract,s_->initial_gas,step,
                                 previous?previous->native_:nullptr);
 }
 CoupledTrial(CoupledTrial&&) noexcept=default;
 CoupledTrial&operator=(CoupledTrial&&) noexcept=default;
 CoupledTrial(const CoupledTrial&)=delete;CoupledTrial&operator=(const CoupledTrial&)=delete;
 n::PrivateFullTimeStep&outer(){active();return s_->outer;}
 const n::ArenaJson&stage_receipts()const {need(bool(s_),"moved coupled trial");return s_->stages;}
 std::shared_ptr<const Continuation> continuation()const {need(s_ && s_->published && s_->continuation,"common pair not published");return s_->continuation;}
 const n::ArenaJson&receipt()const {need(s_ && s_->final_native,"common publication not prepared");return s_->receipt;}
 std::shared_ptr<const n::FullTimeChemistryStageReceipt> run_stage(n::CollisionScope&scope) {
  try {
   active();const unsigned expected=s_->completed+1;need(expected<=2,"both stages already complete");
   const auto&ex=scope.exposure();need(ex.common.accepted_step==s_->step &&
     (ex.slot==n::FullTimeSlot::before_transport?1u:2u)==expected,"wrong common collision stage");
   if(expected==1)s_->initial_native=ex.common.accepted_arena_sha256;
   else need(ex.common.accepted_arena_sha256==s_->initial_native,"native common accepted base changed");
   need(s_->gas->accepted_sha256()==s_->initial_gas,"gas accepted base changed");
   need(expected==1?s_->gas->matches_candidate(s_->initial_binding,s_->step):(s_->last_witness && StageWitness::valid(s_->last_witness.get())),"gas candidate changed outside prior stage transcript");
   const auto rows=scope.rate_population(s_->cap);water_scope(rows);if(expected==2)motion_scope(rows);
   const auto gb=owner_stock(*s_->gas);
   auto provider=std::make_shared<native_gas_targets::Provider>(s_->gas,s_->registry,s_->blob,s_->cfg,s_->options);
   auto token=n::FullTimeCoupledAccess::prepare(s_->outer,scope,provider);
   need(scope.can_inspect(token) && !scope.can_publish(token),"provisional token escaped quarantine");
   token.validate_candidate_host_targets();
   auto population=scope.candidate_population(token,s_->cap,1u<<28);
   auto stock=scope.candidate_stocks(token,1000000,1u<<28);
   auto batch=a::prepare(scope,token,*s_->gas,s_->binding,1000000);batch.validate_current(scope,token,*s_->gas);
   const auto pre=before_inventory(rows,*population,stock->before),post=after_inventory(*population,stock->after);
   const auto pb=pending(stock->before,*population,s_->blob),pa=pending(stock->after,*population,s_->blob);
   const auto refs=reference_operands(batch,*population,s_->blob);
   auto stage=token.candidate_chemistry_stage();need(stage && stage->stage==expected,"native stage identity");
   const auto gas_before=s_->gas->candidate_sha256();
   const auto gas_after=s_->gas->take(s_->step,expected,batch.recoils(),batch.flows());
   GuardedProofJson proof_guard(closure(pre,post,pb,pa,gb,owner_stock(*s_->gas),refs));auto&proof=proof_guard.get();
   need(proof.at("energy_exact_zero").get<bool>() && proof.at("momentum_exact_zero").get<bool>(),"actual native/gas accounted ledger imbalance");
   proof["stage"]=expected;proof["accepted_step"]=s_->step;proof["initial_native"]=s_->initial_native;proof["initial_gas"]=s_->initial_gas;
   proof["gas_before"]=gas_before;proof["gas_after"]=gas_after;proof["native_population"]=population->sha256;
   proof["native_stock"]=stock->sha256;proof["admitted_export"]=stage->admitted_export->sha256;
   proof["native_candidate"]=token.candidate_state_prefix();proof["contract"]=s_->contract;
   proof["events"]=stage->admitted_events;proof["released_packets"]=stage->released_packets;
   proof["quarantine_verified"]=true;proof["feedback_changed_slots"]=provider->feedback().changed;
   const auto hash=native_gas_targets::sha(proof.dump());proof["sha256"]=hash;
   auto witness=std::make_shared<StageWitness>();witness->gas=s_->gas;witness->binding=s_->gas->candidate_binding();witness->step=s_->step;
   n::FullTimeCoupledAccess::authorize_stage(scope,token,witness,&StageWitness::valid,hash);
   // Prepare all record allocations before any private publication.
   GuardedProofJson next_guard;auto&next=next_guard.get();clone_proof_json(next,s_->stages);next.push_back(nullptr);clone_proof_json(next.back(),proof);
   need(scope.can_publish(token) && StageWitness::valid(witness.get()),"dead common stage witness");
   need(scope.try_publish(token),"private coupled stage publication failed");
   s_->stages.swap(next);s_->last_witness=std::move(witness);++s_->completed;if(expected==1)s_->first_population=std::move(population);
   return stage;
  }catch(...){fail();throw;}
 }
 void seal(const n::FullTimeCommitSeals&seals) {
  try {
   active();need(s_->completed==2 && s_->stages.size()==2,"missing stage transcript");
   need(s_->last_witness && StageWitness::valid(s_->last_witness.get()),"stage2 gas candidate changed before final preparation");
   const auto candidate=s_->gas->candidate_sha256();
   need(s_->stages.at(1).at("gas_after")==candidate && s_->gas->accepted_sha256()==s_->initial_gas,"gas state changed after closure");
   auto gas=std::make_shared<FinalWitness>();gas->gas=s_->gas;gas->prepared=s_->gas->prepare_commit(candidate);
   auto native=std::make_unique<n::PreparedArena>(s_->outer.seal(seals));
   GuardedProofJson record_guard;auto&record=record_guard.get();
   record["schema"]="native-gas-coupled-step-v1";record["contract"]=s_->contract;record["accepted_step"]=s_->step;
   record["parent_common_receipt"]=s_->parent_receipt;record["native_before"]=s_->initial_native;record["gas_before"]=s_->initial_gas;
   record["native_after"]=native->candidate_state_prefix();record["gas_after"]=candidate;
   clone_proof_json(record["stages"],s_->stages);record["production_qualified"]=false;
   record["scope"]="CPU bounded same-cell inventory motion; complete two-stage accounting, not full network/transport qualification";
   const auto hash=native_gas_targets::sha(record.dump());record["sha256"]=hash;
   n::FullTimeCoupledAccess::authorize_final(s_->outer,*native,gas,&FinalWitness::valid,hash);
   auto continuation=std::shared_ptr<Continuation>(new Continuation);
   continuation->native_=n::FullTimeCoupledAccess::prepare_continuation(s_->outer,*native);
   continuation->gas_=s_->gas;continuation->gas_after_=candidate;continuation->common_receipt_=hash;
   continuation->next_step_=s_->step+1;
   s_->receipt.swap(record);s_->continuation=std::move(continuation);
   s_->final_gas=std::move(gas);s_->final_native=std::move(native);
  }catch(...){fail();throw;}
 }
 bool can_publish()const noexcept {
  return s_ && !s_->failed && !s_->published && s_->final_native && s_->final_gas &&
   s_->rollback.active() && FinalWitness::valid(s_->final_gas.get()) &&
   n::FullTimeCoupledAccess::can_publish_final(s_->outer,*s_->final_native);
 }
 bool try_publish()noexcept {
  if(!can_publish())return false;
  if(!n::FullTimeCoupledAccess::publish_final(s_->outer,*s_->final_native))return false;
  if(!s_->gas->try_publish(s_->final_gas->prepared))std::terminate();
  s_->published=true;return true;
 }
};
} // namespace native_gas_common
