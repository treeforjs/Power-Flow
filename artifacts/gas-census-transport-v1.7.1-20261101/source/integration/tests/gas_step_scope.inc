// Included by the gas-owner qualification driver, using its actual owners and
// stage/flow fixtures. These controls cover the coordinator's rollback scope;
// they do not qualify a jointly coupled physical step.
#include <pf_coupled/gas_step_scope.hpp>
#include <type_traits>

using pf_coupled::GasStepScope;
static_assert(!std::is_copy_constructible_v<GasStepScope>);
static_assert(!std::is_copy_assignable_v<GasStepScope>);
static_assert(std::is_nothrow_move_constructible_v<GasStepScope>);
static_assert(std::is_nothrow_move_assignable_v<GasStepScope>);
static_assert(std::is_nothrow_destructible_v<GasStepScope>);
static_assert(!std::is_copy_constructible_v<Owner::AttemptBinding>);
static_assert(noexcept(std::declval<Owner&>().abandon_attempt(std::declval<Owner::AttemptBinding&>())));

void scope_nonempty(Owner &o, std::uint64_t step) {
  auto f=debit(0,2,true);
  o.take(step,1,{},{&f,1});
  o.take(step,2,{},{});
}

void scope_tests(Execution execution) {
  const auto cfg=config(execution);
  auto owner=std::make_shared<Owner>(cfg);
  const auto genesis=owner->encode();
  Owner::AttemptBinding empty;
  check(!owner->matches_attempt(empty) && !owner->abandon_attempt(empty),"empty rollback authority refuses");
  bool null_refused=false;
  try{auto invalid=GasStepScope::begin({},0);}catch(const std::invalid_argument&){null_refused=true;}
  check(null_refused,"null scope owner refused");

  // Failure after actual gas mutation unwinds the current candidate, including
  // the prepared-publication state; no accepted byte may change.
  for(bool seal:{false,true}) {
    bool raised=false;
    try {
      auto scope=GasStepScope::begin(owner,0);
      scope_nonempty(*owner,0);
      auto candidate=owner->candidate_sha256();
      Owner::PreparedCommit prepared;
      if(seal) prepared=owner->prepare_commit(candidate);
      check(scope.active(),"rollback scope survives mutations and preparation");
      throw std::runtime_error("manufactured later-owner refusal");
    } catch(const std::runtime_error& e) {
      raised=std::string(e.what())=="manufactured later-owner refusal";
    }
    check(raised && !owner->step_open() && owner->encode()==genesis,"later-owner refusal rolls back exact accepted bytes");
  }

  // Old trial destruction must not cancel a same-digest retry, even if both
  // belong to the same accepted step and all numerical bytes match.
  auto old=GasStepScope::begin(owner,0);
  scope_nonempty(*owner,0);
  const auto original=owner->candidate_sha256();
  owner->abandon();
  auto retry=GasStepScope::begin(owner,0);
  scope_nonempty(*owner,0);
  check(owner->candidate_sha256()==original && !old.active() && retry.active(),"equal-content retry has a distinct rollback lifetime");
  old.reset();
  check(retry.active() && owner->candidate_sha256()==original,"old scope reset preserves newer retry");
  retry.reset();
  check(!owner->step_open() && owner->encode()==genesis,"retry itself can roll back");

  // Beginning a second scope must fail without disturbing the first.
  auto active=GasStepScope::begin(owner,0);
  auto read_binding=owner->candidate_binding();
  refuses([&]{auto wrong=GasStepScope::begin(owner,0);},"nested scope refused");
  check(active.active() && owner->matches_candidate(read_binding),"nested begin refusal preserves live attempt");
  auto moved=std::move(active);
  check(!active.active() && moved.active(),"scope move construction transfers cleanup");
  active.reset();
  check(moved.active(),"moved-from destruction preserves live attempt");
  auto assign=[](GasStepScope &to, GasStepScope &from){to=std::move(from);};
  assign(moved,moved);
  check(moved.active(),"self move preserves scope");

  auto other=std::make_shared<Owner>(cfg);
  const auto other_genesis=other->encode();
  auto destination=GasStepScope::begin(other,0);
  scope_nonempty(*other,0);
  destination=std::move(moved);
  check(!other->step_open() && other->encode()==other_genesis && destination.active() && !moved.active(),"move assignment first rolls back its previous owner");
  destination.reset();

  // A foreign owner and an orphaned lifetime cannot consume the right owner's
  // token. Failed attempt matching must not mutate either token or state.
  auto attempt=owner->begin_attempt(0);
  auto foreign=other->begin_attempt(0);
  check(!other->matches_attempt(attempt) && !other->abandon_attempt(attempt),"foreign rollback authority refused");
  check(owner->matches_attempt(attempt) && other->matches_attempt(foreign),"foreign refusal leaves both candidates live");
  check(owner->abandon_attempt(attempt) && !owner->abandon_attempt(attempt),"exact attempt can roll back only once");
  other->abandon_attempt(foreign);
  Owner::AttemptBinding orphan;
  {Owner destroyed(cfg);orphan=destroyed.begin_attempt(0);}
  auto fresh=owner->begin_attempt(0);
  check(!owner->matches_attempt(orphan) && !owner->abandon_attempt(orphan) && owner->matches_attempt(fresh),"destroyed-owner token cannot affect a fresh owner");
  owner->abandon_attempt(fresh);

  // A poisoned key-setup candidate still belongs to the trial and must be
  // abandoned during unwinding. Content/read bindings intentionally refuse it.
  auto limited_cfg=cfg;limited_cfg.max_entries=1;
  auto limited=std::make_shared<Owner>(limited_cfg);
  const auto limited_base=limited->encode();
  refuses([&]{
    auto scope=GasStepScope::begin(limited,0);
    std::array<FlowRecord,2> batch={debit(0,2,true),debit(1,2,true)};
    limited->prepare_keys({},batch,{});
  },"key capacity refuses inside scope");
  check(!limited->step_open() && limited->encode()==limited_base,"poisoned candidate cleaned up on failure");

  // Cleanup itself must allocate nothing, even when it destroys the last owner
  // reference. This tests exception-unwinding safety under memory pressure.
  auto cleanup=GasStepScope::begin(owner,0);
  scope_nonempty(*owner,0);
  const auto allocations=allocation_probe::allocations;
  allocation_probe::fail_after=0;
  cleanup.reset();
  allocation_probe::fail_after=-1;
  check(allocations==allocation_probe::allocations && !owner->step_open() && owner->encode()==genesis,"scope cleanup allocates nothing and preserves accepted bytes");
  std::weak_ptr<Owner> weak;
  GasStepScope sole;
  {
    auto temporary=std::make_shared<Owner>(cfg);weak=temporary;
    sole=GasStepScope::begin(temporary,0);
    scope_nonempty(*temporary,0);
  }
  check(!weak.expired() && sole.active(),"scope retains owner lifetime");
  sole.reset();
  check(weak.expired(),"scope releases last owner safely");

  // Compare real populated accepted history and restored continuation against
  // the existing explicit begin/publish path. Scope destruction after publish
  // must be harmless, including when the next candidate is already open.
  Owner reference(cfg);
  for(std::uint64_t step=0;step<2;++step) {
    auto scope=GasStepScope::begin(owner,step);
    scope_nonempty(*owner,step);
    auto prepared=owner->prepare_commit(owner->candidate_sha256());
    stage_nonempty(reference,step);reference.commit(reference.candidate_sha256());
    check(owner->try_publish(prepared) && !scope.active(),"published attempt no longer has rollback authority");
    check(owner->encode()==reference.encode(),"scoped populated state matches explicit publication");
    auto next=GasStepScope::begin(owner,step+1);
    scope.reset();
    check(next.active(),"published old scope cannot abandon following step");
  }
  auto restored=std::shared_ptr<Owner>(Owner::restore(owner->encode(),cfg,owner->accepted_sha256()));
  for(auto o:{owner,restored}) {
    auto scope=GasStepScope::begin(o,2);scope_nonempty(*o,2);
    auto prepared=o->prepare_commit(o->candidate_sha256());check(o->try_publish(prepared),"scoped continuation published");
  }
  check(owner->encode()==restored->encode(),"nonempty scoped restored continuation exact");

  if(execution==Execution::host) {
    bool complete=false;unsigned failures=0;
    for(long fault=0;fault<128;++fault) {
      const auto accepted=owner->encode();bool failed=false;
      allocation_probe::fail_after=fault;
      try{auto scope=GasStepScope::begin(owner,3);}catch(const std::bad_alloc&){failed=true;}
      allocation_probe::fail_after=-1;
      check(!owner->step_open() && owner->encode()==accepted,"begin allocation failure leaves no open attempt or changed accepted state");
      if(!failed){complete=true;break;}++failures;
    }
    check(complete && failures>0,"scoped begin allocation failure sweep exercised");
    std::cout<<"scope_begin_allocation_failures "<<failures<<"\n";
  }
#ifdef PF_GAS_BEGIN_FAULT_WRAPPER
  if(execution!=Execution::host) {
    for(int mode:{1,2}) {
      const auto accepted=owner->encode();bool failed=false;
      device_begin_failure_mode=mode;
      try{auto scope=GasStepScope::begin(owner,3);}catch(const std::runtime_error&){failed=true;}
      device_begin_failure_mode=0;
      check(failed && !owner->step_open() && owner->encode()==accepted,"device begin failure automatically rolls back its own attempt");
      auto recovered=GasStepScope::begin(owner,3);check(recovered.active(),"device begin failure permits fresh retry");
    }
  }
#endif
  std::cout<<"scope_state "<<owner->accepted_sha256()<<"\n";
}
