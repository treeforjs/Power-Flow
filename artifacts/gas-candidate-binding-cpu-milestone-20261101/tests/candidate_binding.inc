// Lifetime capability controls; no physics or numerical-identity changes.
#include <type_traits>
static_assert(!std::is_copy_constructible_v<Owner::CandidateBinding>);
static_assert(std::is_nothrow_move_constructible_v<Owner::CandidateBinding>);
static_assert(noexcept(std::declval<const Owner&>().matches_candidate(std::declval<const Owner::CandidateBinding&>())));
static_assert(noexcept(std::declval<const Owner&>().matches_candidate(std::declval<const Owner::CandidateBinding&>(),0)));
#ifdef PF_GAS_BEGIN_FAULT_WRAPPER
extern int device_begin_failure_mode;
#endif
void binding_tests(Execution execution) {
  const auto cfg=config(execution);
  Owner a(cfg),equal(cfg);
  Owner::CandidateBinding empty;
  check(!a.matches_candidate(empty),"default binding closed refusal");
  refuses([&]{(void)a.candidate_binding();},"capture closed candidate");
  a.begin_step(0);equal.begin_step(0);
  auto token=a.candidate_binding();const auto base=a.accepted_sha256(),candidate=a.candidate_sha256();
  check(token.accepted_step()==0 && a.matches_candidate(token,0) && !a.matches_candidate(token,1),"expected native step binding");
  check(!a.matches_candidate(empty) && equal.candidate_sha256()==candidate && !equal.matches_candidate(token),"equal-content foreign owner refused");
  auto moved=std::move(token);check(!a.matches_candidate(token) && a.matches_candidate(moved),"move invalidates original binding");
  (void)a.entries(true);(void)a.ledgers(true);
  std::array<CellTargetsAbi,1> one={cell()};std::array<unsigned,1> atlas={0};a.feedback(0,1,atlas,one);
  check(a.matches_candidate(moved) && a.accepted_sha256()==base && a.candidate_sha256()==candidate,"read-only queries preserve binding and contents");
  const auto allocations=allocation_probe::allocations;
  allocation_probe::fail_after=0;
  const bool match=a.matches_candidate(moved),correct=a.matches_candidate(moved,0),wrong=a.matches_candidate(moved,1),foreign=equal.matches_candidate(moved);
  allocation_probe::fail_after=-1;
  check(match && correct && !wrong && !foreign && allocations==allocation_probe::allocations,"preflight allocates nothing");
  a.prepare_keys({},{},{});
  check(a.candidate_sha256()==candidate && !a.matches_candidate(moved),"same-digest key preparation invalidates old binding");
  auto keyed=a.candidate_binding();
  if(execution==Execution::host) {
    auto part=a.accumulate({},{});check(a.matches_candidate(keyed),"const accumulation preserves binding");
    a.merge(0,1,{&part,1});check(!a.matches_candidate(keyed),"partial merge invalidates binding");
  } else {
    a.take(0,1,{},{});check(!a.matches_candidate(keyed),"take invalidates binding");
  }
  auto stage1=a.candidate_binding();
  refuses([&]{a.take(0,1,{},{});},"wrong stage refused");
  check(!a.matches_candidate(stage1),"failed attempted stage invalidates old binding");
  auto beforeflows=a.candidate_binding();a.apply_flows(0,{});
  check(!a.matches_candidate(beforeflows),"external flow call invalidates binding");
  auto beforetransport=a.candidate_binding();refuses([&]{a.transport(0,{});},"disabled transport refuses");
  check(!a.matches_candidate(beforetransport),"rejected transport invalidates old binding");
  a.take(0,2,{},{});auto beforeprepare=a.candidate_binding();const auto done=a.candidate_sha256();
  refuses([&]{a.prepare_commit(std::string(64,'0'));},"wrong digest prepare refuses");
  check(!a.matches_candidate(beforeprepare) && a.candidate_sha256()==done,"failed prepare invalidates binding without changing contents");
  auto beforeseal=a.candidate_binding();auto publication=a.prepare_commit(done);
  check(!a.matches_candidate(beforeseal) && a.can_publish(publication),"prepared publication closes read-only binding");
  refuses([&]{(void)a.candidate_binding();},"capture sealed candidate");
  check(a.try_publish(publication) && !a.matches_candidate(beforeseal),"publication closes candidate binding");
  refuses([&]{(void)a.candidate_binding();},"capture published candidate");
  auto retry_old=equal.candidate_binding();equal.abandon();equal.begin_step(0);
  check(equal.candidate_sha256()==candidate && !equal.matches_candidate(retry_old),"same-digest abandon/rebegin invalidates binding");
  check(equal.matches_candidate(equal.candidate_binding()),"fresh retry binding accepted");
  Owner::CandidateBinding orphan;
  {Owner dead(cfg);dead.begin_step(0);orphan=dead.candidate_binding();}
  check(!equal.matches_candidate(orphan),"destroyed equal-owner binding refused");
  auto restored=Owner::restore(a.encode(),cfg,a.accepted_sha256());restored->begin_step(1);
  check(!restored->matches_candidate(beforeseal) && restored->matches_candidate(restored->candidate_binding(),1),"restore creates new capability lifetime");
  {
    auto limited=cfg;limited.max_entries=1;Owner o(limited);o.begin_step(0);
    auto old=o.candidate_binding();const auto accepted=o.accepted_sha256();
    std::array<FlowRecord,2> batch={debit(0,2,true),debit(1,2,true)};
    refuses([&]{o.prepare_keys({},batch,{});},"key capacity refusal");
    check(!o.matches_candidate(old) && o.accepted_sha256()==accepted,"key capacity preserves accepted and invalidates binding");
    refuses([&]{(void)o.candidate_binding();},"key capacity failure requires abandon");
    o.abandon();o.begin_step(0);check(o.matches_candidate(o.candidate_binding()),"key capacity fresh retry binding");
  }
  if(execution==Execution::host) {
    unsigned poisoned=0,failures=0;bool complete=false;
    for(long fault=0;fault<128;++fault) {
      Owner o(cfg);o.begin_step(0);auto old=o.candidate_binding();const auto accepted=o.accepted_sha256();
      auto f=debit(0,2,true);bool bad=false;
      allocation_probe::fail_after=fault;
      try{o.prepare_keys({},{&f,1},{});}catch(const std::bad_alloc&){bad=true;}
      allocation_probe::fail_after=-1;
      check(o.accepted_sha256()==accepted,"key allocation failure preserves accepted state");
      if(!bad){complete=true;check(!o.matches_candidate(old),"successful keys invalidate binding");break;}
      ++failures;bool capture_refused=false;
      try{(void)o.candidate_binding();}catch(const Refusal&){capture_refused=true;}
      if(capture_refused){++poisoned;check(!o.matches_candidate(old),"poisoned setup refuses old binding");}
      o.abandon();o.begin_step(0);check(o.matches_candidate(o.candidate_binding()),"failed key setup can retry only on fresh lifetime");
    }
    check(complete && failures>0 && poisoned>0,"allocation sweep exercises poisoned key preparation");
    std::cout<<"binding_key_allocation_failures "<<failures<<" poisoned "<<poisoned<<"\n";
  }
#ifdef PF_GAS_BEGIN_FAULT_WRAPPER
  if(execution!=Execution::host) {
    Owner o(cfg);stage_nonempty(o,1);o.commit(o.candidate_sha256());const auto accepted=o.accepted_sha256();
    for(int mode:{1,2}) {
      bool rejected=false;device_begin_failure_mode=mode;
      try{o.begin_step(2);}catch(const std::runtime_error&){rejected=true;}
      device_begin_failure_mode=0;
      check(rejected && o.accepted_sha256()==accepted,"device begin fault preserves accepted state");
      check(!o.matches_candidate(empty),"failed device begin cannot preflight");
      refuses([&]{(void)o.candidate_binding();},"failed device begin cannot capture");
      refuses([&]{o.take(2,1,{},{});},"failed device begin must abandon");
      o.abandon();o.begin_step(2);check(o.matches_candidate(o.candidate_binding(),2),"device begin fresh retry is live");o.abandon();
    }
  }
#endif
  std::cout<<"binding_state "<<a.accepted_sha256()<<"\n";
}
