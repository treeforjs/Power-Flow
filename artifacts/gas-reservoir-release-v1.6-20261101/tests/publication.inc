void stage_nonempty(Owner &o,std::uint64_t step) {
  o.begin_step(step);auto f=debit(0,2,true);o.take(step,1,{},{&f,1});o.take(step,2,{},{});
}
void publication_tests(Execution execution) {
  const auto cfg=config(execution);Owner a(cfg),b(cfg);
  stage_nonempty(a,1);stage_nonempty(b,1);
  const auto accepted=a.accepted_sha256(),candidate=a.candidate_sha256();
  refuses([&]{a.prepare_commit(std::string(64,'0'));},"wrong publication digest");
  auto token=a.prepare_commit(candidate);
  check(token.candidate_sha256()==candidate && a.accepted_sha256()==accepted,"prepare leaves accepted state");
  check(a.can_publish(token) && !b.can_publish(token),"token owner identity");
  Owner::PreparedCommit empty;
  check(!a.try_publish(empty) && !b.try_publish(token),"empty/foreign tokens refused");
  auto moved=std::move(token);
  check(!a.can_publish(token) && a.can_publish(moved),"moved token invalidated");
  refuses([&]{a.prepare_commit(candidate);},"duplicate preparation");
  refuses([&]{a.apply_flows(1,{});},"sealed candidate mutation");
  refuses([&]{a.take(1,2,{},{});},"sealed candidate stage");
  refuses([&]{a.prepare_keys({},{},{});},"sealed candidate keys");
  const auto alloc_before=allocation_probe::allocations;
  allocation_probe::fail_after=0;
  const bool ready=a.can_publish(moved),published=a.try_publish(moved);
  allocation_probe::fail_after=-1;
  check(ready && published && alloc_before==allocation_probe::allocations,"publication allocates nothing");
  check(a.accepted_sha256()==candidate && !a.step_open(),"prepared state published exactly");
  check(!a.try_publish(moved),"token replay refused");
  auto restored=Owner::restore(a.encode(),cfg,a.accepted_sha256());
  check(restored->accepted_sha256()==candidate,"published restore identity");
  auto abandoned=b.prepare_commit(b.candidate_sha256());b.abandon();
  stage_nonempty(b,1);
  check(b.candidate_sha256()==candidate && !b.can_publish(abandoned),"same-digest retry rejects stale capability");
  auto retry=b.prepare_commit(candidate);check(b.try_publish(retry),"fresh retry token accepted");
  check(b.accepted_sha256()==a.accepted_sha256(),"retry accepted identity");
  Owner::PreparedCommit orphan;
  {Owner destroyed(cfg);stage_nonempty(destroyed,1);orphan=destroyed.prepare_commit(destroyed.candidate_sha256());}
  check(!a.can_publish(orphan),"destroyed-owner capability refused");
  for(Owner *owner : {&a,restored.get()}) {
    stage_nonempty(*owner,2);auto next=owner->prepare_commit(owner->candidate_sha256());
    check(owner->try_publish(next),"continued prepared publication");
  }
  check(a.accepted_sha256()==restored->accepted_sha256(),"continued restore identity");
  if(execution==Execution::host) {
    bool reached_success=false;unsigned failures=0;
    for(long fault=0;fault<128;++fault) {
      Owner owner(cfg);stage_nonempty(owner,1);const auto base=owner.accepted_sha256(),wanted=owner.candidate_sha256();
      Owner::PreparedCommit prepared;bool failed=false;
      allocation_probe::fail_after=fault;
      try {prepared=owner.prepare_commit(wanted);}catch(const std::bad_alloc &){failed=true;}
      allocation_probe::fail_after=-1;
      check(owner.accepted_sha256()==base,"allocation failure preserves accepted state");
      check(owner.candidate_sha256()==wanted,"preparation preserves candidate identity");
      if(failed) {
        ++failures;check(!owner.can_publish(prepared),"failed preparation gives no capability");
        auto recover=owner.prepare_commit(wanted);check(owner.try_publish(recover),"retry after allocation failure");
      }else {check(owner.try_publish(prepared),"fault sweep reaches successful preparation");reached_success=true;break;}
    }
    check(reached_success && failures>0,"preparation allocation-fault coverage");
    std::cout<<"publication_allocation_failures "<<failures<<"\n";
  }
  std::cout<<"publication_state "<<a.accepted_sha256()<<"\n";
}
