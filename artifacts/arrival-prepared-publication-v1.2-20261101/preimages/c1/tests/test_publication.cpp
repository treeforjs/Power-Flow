#include <pf_arrival_v11/arrival_owner.hpp>
#include "allocation_probe.hpp"
#include <iostream>
#include <type_traits>

namespace pa=pf_arrival_v11;
int checks=0;
void check(bool ok,const char* label){if(!ok)throw std::runtime_error(label);++checks;}
template<class F> void refuses(F f,const char* label){bool caught=false;try{f();}catch(const pa::Refusal&){caught=true;}check(caught,label);}
std::shared_ptr<const pa::Policy> policy(){
  pa::Configuration c;c.run_seed=98325;
  c.streams={{17,pa::Kind::electron,0,"electron",2.5,128},{18,pa::Kind::neutral,1,"H2O",.75,128}};
  return std::make_shared<const pa::Policy>(c);
}
pa::Staged stage(pa::Owner& o,bool empty=false){
  const auto k=o.next_step();const double t=double(k)*1e-12,dt=1e-12;
  pa::StepInput input;input.step=k;
  const double e=empty?0:9.25+.125*double(k%3),n=empty?0:4.125+.125*double(k%4);
  input.electrons.push_back({17,e,e>0?std::log(e):-INFINITY,true,t,dt,t+dt});
  input.neutrals.push_back({18,pa::Count(n),n>0?std::log(n):-INFINITY,pa::Count(1e6),t,t+dt});
  auto s=o.stage(o.plan(input),o.next_ordinal()+1);
  std::vector<double> times;
  for(const auto& b:s.births())if(b.stream_id==18)times.push_back(t+.5*dt);
  s.bind_source_times(18,times);
  return s;
}
pa::NativeBirthReceipt receipt(const pa::Owner& o,const pa::Staged& s){
  pa::NativeBirthReceipt r;r.base=s.native_base();r.next_id_after=r.base+s.births().size();
  for(const auto& b:s.births())r.entries.push_back({b.packet_id(),o.policy().stream(b.stream_id).species_state,b.native_id,b.weight});
  return r;
}
void parity_history(){
  auto p=policy();pa::Owner a(p,0);
  std::string records;
  for(int k=0;k<24;++k){auto s=stage(a,k%7==0);auto r=receipt(a,s);records+=a.commit(std::move(s),r).bytes;}
  auto restored=pa::Owner::restore(p,a.save());
  for(int k=0;k<8;++k){
    for(pa::Owner* o:{&a,restored.get()}){auto s=stage(*o);auto r=receipt(*o,s);o->commit(std::move(s),r);}
    check(a.save()==restored->save(),"nonempty continuation after restore");
  }
  std::cout<<"history "<<a.state_sha256()<<" "<<pa::sha256_hex(records)<<" "<<pa::sha256_hex(a.save())<<"\n";
}
#ifdef PF_TEST_PREPARED
static_assert(!std::is_copy_constructible_v<pa::Owner::PreparedCommit>);
static_assert(std::is_nothrow_move_constructible_v<pa::Owner::PreparedCommit>);
static_assert(std::is_nothrow_move_assignable_v<pa::Owner::PreparedCommit>);
static_assert(noexcept(std::declval<pa::Owner&>().try_publish(std::declval<pa::Owner::PreparedCommit&>())));
void prepared_tests(){
  auto p=policy();pa::Owner a(p,0),b(p,0),reference(p,0);
  const auto initial=a.save();
  pa::Owner::PreparedCommit empty;
  check(!a.can_publish(empty)&&!a.try_publish(empty)&&!empty.staged(),"default prepared token refuses");
  auto s=stage(a);auto r=receipt(a,s);const auto plan=s.plan().sha256();
  auto wrong=r;wrong.entries.front().weight+=1;
  refuses([&]{a.prepare_commit(std::move(s),wrong);},"wrong native weight refuses preparation");
  check(a.save()==initial && s.plan().sha256()==plan && !s.births().empty(),"failed preparation preserves owner and staged input");
  auto token=a.prepare_commit(std::move(s),r);
  check(a.save()==initial && a.can_publish(token),"preparation does not publish");
  check(token.staged() && token.staged()->births().size()==r.entries.size(),"prepared token retains immutable staged mapping");
  check(!b.can_publish(token)&&!b.try_publish(token)&&a.can_publish(token),"equal-state foreign owner cannot publish token");
  auto moved=std::move(token);
  check(!a.can_publish(token)&&a.can_publish(moved),"move invalidates former publication token");
  auto sr=stage(reference);auto rr=receipt(reference,sr);const auto ledger=reference.commit(std::move(sr),rr);
  check(moved.record().bytes==ledger.bytes && moved.record().record_sha256==ledger.record_sha256 && moved.record().chain_head==ledger.chain_head,"prepared ledger equals compatibility commit");
  const auto expected=moved.candidate_state_sha256();
  const auto allocations=allocation_probe::allocations;
  allocation_probe::fail_after=0;
  const bool ready=a.can_publish(moved),published=a.try_publish(moved),again=a.try_publish(moved);
  allocation_probe::fail_after=-1;
  check(ready&&published&&!again&&allocations==allocation_probe::allocations,"whole arrival publication allocates nothing");
  check(a.save()==reference.save()&&a.state_sha256()==expected,"prepared state exactly equals actual accepted state");

  // Multiple immutable proposals may coexist for one accepted state. The first
  // publication advances generation and invalidates every other proposal.
  auto s1=stage(a),s2=stage(a);auto r1=receipt(a,s1),r2=receipt(a,s2);
  auto first=a.prepare_commit(std::move(s1),r1),second=a.prepare_commit(std::move(s2),r2);
  check(a.try_publish(first)&&!a.can_publish(second)&&!a.try_publish(second),"accepted advance invalidates competing prepared state");
  auto restored=pa::Owner::restore(p,a.save());auto s3=stage(a);auto r3=receipt(a,s3);auto live=a.prepare_commit(std::move(s3),r3);
  check(!restored->can_publish(live),"restored equal state has a distinct owner lifetime");
  auto before=a.save();
  {auto discarded=std::move(live);check(a.can_publish(discarded),"discarded private proposal was live");}
  check(a.save()==before,"discarding prepared proposal preserves accepted state");
  pa::Owner::PreparedCommit orphan;
  {pa::Owner dead(p,0);auto sd=stage(dead);auto rd=receipt(dead,sd);orphan=dead.prepare_commit(std::move(sd),rd);}
  check(!b.can_publish(orphan)&&!b.try_publish(orphan),"destroyed owner's token cannot publish into another owner");

  // Exact prepared hashes/bytes across empty and populated steps and restart.
  for(int k=0;k<10;++k){
    for(pa::Owner* o:{&a,restored.get()}){
      auto st=stage(*o,k%3==0);auto re=receipt(*o,st);auto t=o->prepare_commit(std::move(st),re);
      const auto wanted=t.candidate_state_sha256();check(o->try_publish(t),"prepared continuation publishes");
      check(o->state_sha256()==wanted,"prepared hash equals accepted hash through continuation");
    }
    check(a.save()==restored->save(),"prepared nonempty restored continuation matches");
  }

  // Enumerate every C++ allocation during preparation. Accepted bytes and the
  // staged input stay intact on each failure; successful retry yields one state.
  bool reached=false;unsigned failures=0;
  for(long fault=0;fault<512;++fault){
    auto st=stage(a);auto re=receipt(a,st);const auto accepted=a.save(),staged=st.plan().encode();
    const auto n=st.births().size();bool failed=false;pa::Owner::PreparedCommit t;
    allocation_probe::fail_after=fault;
    try{t=a.prepare_commit(std::move(st),re);}catch(const std::bad_alloc&){failed=true;}
    allocation_probe::fail_after=-1;
    check(a.save()==accepted,"allocation refusal preserves accepted owner bytes");
    if(!failed){check(a.can_publish(t),"allocation sweep reaches prepared state");reached=true;break;}
    ++failures;check(!a.can_publish(t)&&st.births().size()==n&&st.plan().encode()==staged,"allocation refusal preserves staged input and mints no token");
    auto retry=a.prepare_commit(std::move(st),re);check(a.can_publish(retry),"preparation can retry after allocation failure");
  }
  check(reached&&failures>0,"preparation allocation sweep completed");
  std::cout<<"prepare_allocation_failures "<<failures<<"\n";
}
#endif
int main(){try{
  parity_history();
#ifdef PF_TEST_PREPARED
  prepared_tests();
#endif
  std::cout<<"PASS "<<checks<<" checks\n";return 0;
}catch(const std::exception& e){allocation_probe::fail_after=-1;std::cerr<<"FAIL "<<e.what()<<" after "<<checks<<"\n";return 1;}}
