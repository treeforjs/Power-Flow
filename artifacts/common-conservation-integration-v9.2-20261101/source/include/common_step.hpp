#pragma once
#include "coupled_trial.hpp"
#include "native_genesis.hpp"
#include "source_hall_iteration.hpp"
#include "final_charge_guard.hpp"
#include "common_ledger.hpp"
#include "energy_bound.hpp"

namespace powerflow::common_step_v9 {
namespace cg=::powerflow::coupled_genesis_v9;
namespace sh=::powerflow::source_hall_iteration_v1;
namespace ss=::powerflow::source_stage_v8;
namespace nc=::native_gas_common;
namespace n=::raikou::mfem_backend;
namespace g=::pf_gas_v1;
namespace h=::hall_corrective_v3;
inline void require(bool ok,const char *why){if(!ok)throw std::invalid_argument(std::string("common step v9: ")+why);}
class Step;
struct StepTestAccess;
class Run {
  friend class Step;
  cg::Owners owners_;
  std::shared_ptr<const nc::CoupledTrial::Continuation> continuation_;
  std::string receipt_root_;
  std::string configuration_sha256_;
  std::shared_ptr<const CommonLedger> ledger_;
  explicit Run(cg::Owners o):owners_(std::move(o)){
    // Root policy identity is downstream of the concrete owner graph. It is
    // never fed back into chemistry configuration (which would form a cycle).
    cg::Encoding e;e.text("common-native-gas-source-Hall-v9.2/independent-conservation-v1");
    e.text(owners_.configuration_sha256);e.text("exact-cumulative-material-energy-v1;no-conductor-sinks;held-layers-refused");
    configuration_sha256_=e.hash();receipt_root_=configuration_sha256_;
  }
public:
  static std::shared_ptr<Run> create(cg::Request request){
    auto &chem=request.native_configuration.chemistry;
    require(!chem.diagnostic_field_target_approximation,"unowned target-energy approximation refused");
    require(!chem.host_target_feedback && chem.host_target_contract_sha256.empty() && chem.coupled_target_contract_sha256.empty(),
            "fresh common run requires an unbound base chemistry configuration");
    // Validate declared base provenance before rebinding only mass metadata to
    // the final coupled registry. Do not repair an authored mass/provenance error.
    const auto base_registry=n::full_time_chemistry_registry(chem,request.options);
    const auto base_states=n::full_time_chemistry_states(base_registry);
    for(const auto&port:request.ports)for(const auto&setting:port.slice.settings){
      const auto expected=::native_source_v7::electron_mass_configuration(setting.spectrum,chem.blob_path,chem.blob_sha256,base_registry,base_states);
      require(std::bit_cast<std::uint64_t>(expected.electron_mass_kg)==std::bit_cast<std::uint64_t>(setting.spectrum.electron_mass_kg)&&
        expected.electron_mass_identity==setting.spectrum.electron_mass_identity,"base source electron mass/provenance differs from actual base registry");
    }
    g::Owner prototype(request.gas_configuration);
    chem.coupled_target_contract_sha256=nc::coupled_contract(chem,request.options,prototype);
    chem.host_target_feedback=true;
    chem.host_target_contract_sha256=native_gas_targets::contract(chem,request.options,prototype);
    const auto coupled_registry=n::full_time_chemistry_registry(chem,request.options);
    const auto coupled_states=n::full_time_chemistry_states(coupled_registry);
    for(auto&port:request.ports)for(auto&setting:port.slice.settings)
      setting.spectrum=::native_source_v7::electron_mass_configuration(setting.spectrum,chem.blob_path,chem.blob_sha256,coupled_registry,coupled_states);
    return std::shared_ptr<Run>(new Run(cg::admit(request)));
  }
  const cg::Owners& owners()const noexcept{return owners_;}
  const std::string& receipt_root()const noexcept{return receipt_root_;}
  const std::string& configuration_sha256()const noexcept{return configuration_sha256_;}
  const std::shared_ptr<const CommonLedger>& ledger()const noexcept{return ledger_;}
};
struct PhaseTime {std::string phase;sh::Time time;};
struct Receipt {
  std::uint64_t step=0;
  double begin_s=0,end_s=0;
  unsigned attempts=0;
  std::size_t final_carriers=0;
  std::uint64_t stage1_events=0,stage2_events=0;
  std::string parent,common_clock,native_after,gas_after,field_after,heat_after,source_after,
    native_gas_receipt,source_publication_receipt,final_census,conservation_sha256,sha256;
  double raw_energy_residual_J=0,accounted_energy_residual_J=0,energy_bound_J=0;
  std::vector<sh::IterationReceipt> iterations;
  std::vector<PhaseTime> phases;
};
inline std::vector<n::ArenaRateEntry> rows(const n::FullTimeCandidatePopulation &population,std::size_t cap){
  require(population.carriers<=cap,"actual final native census exceeds cap");
  std::vector<n::ArenaRateEntry> result;result.reserve(population.carriers);
  for(const auto &packet:population.packets)for(const auto &particle:packet.population.carriers){
    require(result.size()<cap,"actual final native census cardinality");
    result.push_back({packet.population.packet_id,packet.population.state,particle});
  }
  require(result.size()==population.carriers,"actual final native census mismatch");
  std::sort(result.begin(),result.end(),[](const auto&a,const auto&b){return a.native.id<b.native.id;});return result;
}
// First integrated reference path refuses conductor correction until its heat
// destination is qualified. It does not silently discard or thermalize it.
inline void require_no_unrouted_conductor_correction(const g::Owner &gas){
  const auto layout=gas.layout();
  for(std::uint32_t sp=0;sp<gas.configuration().species.size();++sp)for(std::uint32_t q=0;q<5;++q)
    require(gas.global_digits(layout.conductor(sp,q),true)==gas.global_digits(layout.conductor(sp,q),false),
            "nonzero conductor gas correction needs a qualified thermal/momentum destination");
}
class Step {
  friend struct StepTestAccess;
  struct State {
    std::shared_ptr<Run> run;
    std::unique_ptr<nc::CoupledTrial> coupled;
    std::unique_ptr<sh::Ready> source_field;
    std::shared_ptr<const n::FullTimeCandidatePopulation> final_population;
    std::shared_ptr<const n::FullTimeAcceptedBoundary> final_boundary;
    std::vector<n::ArenaRateEntry> final_rows;
    ss::AcceptedSourceState before;
    n::FullTimeStepBinding common;
    std::string before_run_root,gas_before,future_run_root;
    Receipt receipt;
    std::shared_ptr<const CommonLedger> before_ledger,future_ledger;
    std::shared_ptr<const cm::CombinedAssessment> conservation;
    bool prevalidated=false,published=false;
  };
  std::unique_ptr<State> state_;
  explicit Step(std::unique_ptr<State>s):state_(std::move(s)){}
  void publish_prevalidated()noexcept{
    auto&s=*state_;auto&o=s.run->owners_;
    if(!s.prevalidated||s.published||!s.coupled->can_publish())std::terminate();
    if(!s.coupled->try_publish())std::terminate();
    o.hall->publish_prevalidated(std::move(s.source_field->proposal));
    s.source_field->source.publish_sources_prevalidated(o.sources);
    // All strings, vectors and the native/gas continuation were prepared before
    // the first swap. This is one sole-writer critical section with no callbacks.
    o.source_states.swap(s.source_field->future.source_states);
    o.source_sha256.swap(s.source_field->future.source_transition);
    o.arrival_state_sha256.swap(s.source_field->future.arrival_state);
    s.run->receipt_root_.swap(s.future_run_root);
    s.run->ledger_.swap(s.future_ledger);
    s.run->continuation_=s.coupled->continuation();
    s.published=true;s.prevalidated=false;
  }
public:
  Step(const Step&)=delete;Step&operator=(const Step&)=delete;
  Step(Step&&)noexcept=default;Step&operator=(Step&&)noexcept=default;
  static Step prepare(std::shared_ptr<Run> run,double dt,const ss::s7::BirthPolicy &policy,
                      unsigned maximum_attempts=4,std::span<const ss::s7::IntervalEvidence>intervals={}){
    require(run&&std::isfinite(dt)&&dt>0&&maximum_attempts>=1&&maximum_attempts<=4,"run, timestep and bounded source retry policy required");
    auto&o=run->owners_;const auto cap=o.native_configuration.limits.maximum_carriers;
    require(dt<=o.model->dt_max(),"Hall timestep cap");
    const auto accepted=o.native->full_time_boundary();require(bool(accepted),"missing native accepted boundary");
    require(accepted->time_s==o.hall->time_s()&&accepted->accepted_step==o.hall->revision()&&
      accepted->field_sha256==o.hall->field_sha256()&&accepted->heat_sha256==o.hall->heat_sha256()&&
      accepted->source_sha256==o.source_sha256&&!o.gas->step_open()&&o.gas->accepted_step()==(accepted->accepted_step?accepted->accepted_step-1:0),
      "accepted native/gas/Hall/source clocks or roots differ");
    ss::AcceptedSourceState source_base{o.source_states,o.arrival_state_sha256,o.source_sha256};
    cg::Encoding clock;clock.text("common-native-gas-source-Hall-step-v9.2");clock.text(run->receipt_root_);
    clock.text(run->configuration_sha256_);clock.text(o.configuration_sha256);clock.text(accepted->accepted_arena_sha256);clock.text(o.gas->accepted_sha256());
    clock.text(accepted->field_sha256);clock.text(accepted->heat_sha256);clock.text(accepted->source_sha256);
    clock.u64(accepted->accepted_step);clock.number(accepted->time_s);clock.number(dt);
    n::FullTimeStepBinding common{accepted->accepted_step,accepted->accepted_arena_revision,accepted->time_s,accepted->time_s+dt,dt,
      accepted->accepted_arena_sha256,accepted->field_sha256,accepted->mesh_sha256,accepted->material_sha256,clock.hash()};
    const auto options=o.native->native_options();const auto stages=o.native->full_time_chemistry().stages;
    const auto accepted_rows=o.native->rate_population(cap);
    const auto gas_stock_before=cm::gas_stock(*o.gas,false);
    const auto blob=pfc::load_blob(o.native_configuration.chemistry.blob_path,o.registry.blob_sha256);
    std::unique_ptr<sh::RetrySeed> seed;std::vector<sh::IterationReceipt> iterations;std::vector<PhaseTime> phases;
    for(unsigned attempt=0;attempt<maximum_attempts;++attempt){
      auto state=std::make_unique<State>();state->run=run;state->before=source_base;state->common=common;
      state->before_run_root=run->receipt_root_;state->gas_before=o.gas->accepted_sha256();state->before_ledger=run->ledger_;
      sh::Clock begin_clock;
      auto outer=o.native->begin_full_time_step(common,cap,64);
      state->coupled=std::make_unique<nc::CoupledTrial>(std::move(outer),o.gas,o.native_configuration.chemistry,options,
        common.accepted_step,cap,run->continuation_);
      phases.push_back({"begin",begin_clock.elapsed()});
      auto&private_outer=state->coupled->outer();
      n::ExposureBinding first{common,n::FullTimeSlot::before_transport,common.begin_s,.5*dt,stages,common.accepted_field_sha256};
      sh::Clock first_clock;auto first_scope=private_outer.first_exposure(first);auto stage1=state->coupled->run_stage(first_scope);
      phases.push_back({"chemistry1",first_clock.elapsed()});
      auto source=sh::prepare(private_outer,o.sources,source_base,common,options,o.registry,policy,cap,4,seed.get(),intervals);
      iterations.insert(iterations.end(),source.iterations.begin(),source.iterations.end());
      if(!source.ready){require(bool(source.retry),"source iteration produced neither ready step nor retry seed");seed=std::move(source.retry);continue;}
      state->source_field=std::move(source.ready);auto&ready=*state->source_field;
      sh::Clock motion_clock;auto moved=ready.prepare_motion(common);
      // v1.3 retains only aggregate correction exits. Until the successor's
      // authentic per-exit receipt is adopted, this initial integration subset
      // refuses every conductor sink, even when its correction might be zero.
      const auto actual_motion=state->coupled->motion_snapshot(moved);
      require(actual_motion && state->coupled->can_inspect_motion(moved,*actual_motion),"missing authentic motion sink evidence");
      if(actual_motion->sinks)for(const auto &exit:actual_motion->sinks->records)
        require(exit.kind!=n::FullTimeSinkKind::collected_at_conductor,"initial v1.3 integration subset refuses conductor sinks pending individual correction routing guard");
      state->coupled->accept_transport(std::move(moved));
      phases.push_back({"native_gas_motion",motion_clock.elapsed()});
      auto second=first;second.slot=n::FullTimeSlot::after_transport;second.frozen_physical_stage_s=common.end_s;
      second.owner_stage_sequence=stages+1;second.frozen_field_sha256=ready.proposal.candidate_field_sha256();
      sh::Clock second_clock;auto second_scope=private_outer.second_exposure(second);auto stage2=state->coupled->run_stage(second_scope);
      phases.push_back({"chemistry2",second_clock.elapsed()});
      require_no_unrouted_conductor_correction(*o.gas);
      sh::Clock seal_clock;
      state->coupled->seal({ready.proposal.candidate_field_sha256(),ready.proposal.candidate_heat_sha256(),ready.source.source_after_sha256(),stage2->state_after_sha256});
      state->final_population=state->coupled->final_candidate_population();state->final_boundary=state->coupled->candidate_boundary();
      require(state->final_population&&state->final_boundary,"missing actual final native population/boundary");
      state->final_rows=rows(*state->final_population,cap);
      const auto charge=::powerflow::final_charge_guard_v1::require_consistent(o.hall->mesh(),options,
        ready.proposal.motion().endpoints(),state->final_rows,cap);
      phases.push_back({"seal_final_charge",seal_clock.elapsed()});
      sh::Clock conservation_clock;
      // These pointers are obtained from this exact retained native/gas trial,
      // never reconstructed from caller-provided receipt hashes.
      const auto first_stocks=state->coupled->stage_stocks(1),second_stocks=state->coupled->stage_stocks(2);
      const auto first_population=state->coupled->stage_population(1);
      const auto movement=state->coupled->accepted_motion_snapshot();
      const auto transport_receipt=state->coupled->motion_transport_receipt();
      const std::array exports{state->coupled->stage_export(1),state->coupled->stage_export(2)};
      require(first_stocks&&second_stocks&&first_population&&movement&&transport_receipt&&state->coupled->can_publish(),"actual conservation operands unavailable or stale");
      const auto native_before=cm::native_before(accepted_rows,*first_population,first_stocks->before,blob);
      const auto native_after=cm::native_after(*state->final_population,second_stocks->after,blob);
      const auto source_stock=cm::source_stocks(o.sources,ready.source,*state->final_population);
      const auto gas_stock_after=cm::gas_stock(*o.gas,true);
      require(cm::pending_inventory(first_stocks->before,*first_population,blob).kinetic_J==gas_stock_before.pending_mirror_J&&
        cm::pending_inventory(second_stocks->after,*state->final_population,blob).kinetic_J==gas_stock_after.pending_mirror_J,
        "actual pending native kinetic stock differs from gas mirror");
      const auto chemistry=cm::admitted_event_operands(exports,*state->final_population,blob,1000000);
      const auto boundary_stock=cm::boundary_stocks(movement->sinks.get(),*transport_receipt,*state->final_population,o.gas->configuration(),common.accepted_step);
      const auto field=::powerflow::common_conservation_v1::field_thermal_assessment(*o.hall,ready.proposal);
      auto assessment=std::make_shared<const cm::CombinedAssessment>(cm::combine(native_before,native_after,source_stock,gas_stock_before,gas_stock_after,boundary_stock,chemistry,field));
      // Actual independent operands define the allowed unchanged component
      // tolerance. No self-reported component residual sets this bound.
      const auto energy_bound=cm::derive_energy_bound(*o.hall,ready.proposal,ready.source,native_before,native_after,boundary_stock,field);
      cm::require_energy_closure(*assessment,energy_bound);
      state->future_ledger=prepare_ledger(state->before_ledger,common,*assessment,source_stock,boundary_stock,chemistry,movement->sinks.get(),*state->final_population,energy_bound.total_J);
      state->conservation=std::move(assessment);
      phases.push_back({"independent_conservation",conservation_clock.elapsed()});
      auto&r=state->receipt;r.step=common.accepted_step;r.begin_s=common.begin_s;r.end_s=common.end_s;r.attempts=attempt+1;
      r.final_carriers=state->final_rows.size();r.stage1_events=stage1->admitted_events;r.stage2_events=stage2->admitted_events;
      r.parent=state->before_run_root;r.common_clock=common.common_clock_sha256;r.native_after=state->final_boundary->accepted_arena_sha256;
      r.gas_after=o.gas->candidate_sha256();r.field_after=ready.proposal.candidate_field_sha256();r.heat_after=ready.proposal.candidate_heat_sha256();
      r.source_after=ready.source.source_after_sha256();r.source_publication_receipt=ready.future.receipt_sha256;
      r.native_gas_receipt=state->coupled->receipt().at("sha256").get<std::string>();r.final_census=charge.final_native_census_sha256;
      r.conservation_sha256=state->future_ledger->sha256;r.raw_energy_residual_J=state->conservation->raw_residual_J.value();
      r.accounted_energy_residual_J=state->conservation->accounted_residual_J.value();r.energy_bound_J=energy_bound.total_J;
      r.iterations=std::move(iterations);r.phases=std::move(phases);
      cg::Encoding seal;seal.text("common-native-gas-source-Hall-receipt-v9.2");seal.text(run->configuration_sha256_);
      for(const auto *v:{&r.parent,&r.common_clock,&r.native_after,&r.gas_after,&r.field_after,&r.heat_after,&r.source_after,&r.source_publication_receipt,&r.native_gas_receipt,&r.final_census,&r.conservation_sha256})seal.text(*v);
      seal.u64(r.step);seal.number(r.begin_s);seal.number(r.end_s);r.sha256=seal.hash();state->future_run_root=r.sha256;
      return Step(std::move(state));
    }
    throw std::runtime_error("source/Hall fixed point did not converge within bounded fresh outer attempts");
  }
  const Receipt& receipt()const{require(bool(state_),"moved step");return state_->receipt;}
  bool preflight()noexcept{
    if(!state_||state_->published)return false;auto&s=*state_;s.prevalidated=false;auto&o=s.run->owners_;
    try{
      if(s.run->ledger_!=s.before_ledger||s.run->receipt_root_!=s.before_run_root||o.source_sha256!=s.before.source_seal||o.source_states!=s.before.source_states||
        o.arrival_state_sha256!=s.before.arrival_state||!s.coupled->can_publish()||
        !o.hall->preflight_final(s.source_field->proposal,s.source_field->census,*s.final_boundary)||
        !s.source_field->source.preflight_sources(o.sources))return false;
      s.prevalidated=true;return true;
    }catch(...){return false;}
  }
  bool try_publish()noexcept{if(!preflight())return false;publish_prevalidated();return true;}
};
} // namespace powerflow::common_step_v9
