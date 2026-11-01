// End-to-end publication allocation control. This runs in its own bounded
// process. The baseline's terminate handler emits an explicit mixed-state
// preimage and exits without a core dump; it is never a production callback.
v8::OwnersV8<hc::Owner>* publication_observed_owner=nullptr;
std::uint64_t publication_hall_before=0,publication_arrival_before=0;
void publication_terminated() noexcept {
  stage1_set_allocation_countdown(-1);
  const auto* o=publication_observed_owner;
  const bool hall_advanced=o && o->hall->revision()==publication_hall_before+1;
  const bool arrival_advanced=o && o->arrival->next_step()==publication_arrival_before+1;
  if(hall_advanced && !arrival_advanced){
    constexpr char message[]="PUBLICATION_PREIMAGE: Hall accepted step advanced; arrival accepted step unchanged after allocation failure\n";
    (void)::write(2,message,sizeof(message)-1);std::_Exit(77);
  }
  constexpr char message[]="PUBLICATION_FAILURE: unexpected accepted-state combination\n";
  (void)::write(2,message,sizeof(message)-1);std::_Exit(78);
}
void publication_fail_next_allocation(int phase){
  if(phase==0)stage1_set_allocation_countdown(0);
}
int publication_guard_mode(const fs::path& out){
  const auto in=load_inputs("qual");auto o=genesis<hc::Owner>(in,"qual");
  auto trial=v8::prepare_v8(o,opt.dt_s);
  if(!v8::can_publish_v8(o,trial))throw std::runtime_error("publication guard requires a valid populated trial");
  publication_observed_owner=&o;publication_hall_before=o.hall->revision();publication_arrival_before=o.arrival->next_step();
  const auto previous=std::set_terminate(&publication_terminated);
  TA::set_publish_probe(trial,&publication_fail_next_allocation);
  const bool published=v8::try_publish_v8(o,trial);
  const auto remaining=stage1_allocation_countdown();
  stage1_set_allocation_countdown(-1);
  std::set_terminate(previous);publication_observed_owner=nullptr;
  if(!published || remaining!=0 || o.hall->revision()!=publication_hall_before+1 || o.arrival->next_step()!=publication_arrival_before+1)
    throw std::runtime_error("whole publication path allocated or did not advance both owners");
  // Compare each prepared bookkeeping digest to its actual accepted owner.
  for(std::size_t i=0;i<o.sources.size();++i)
    if(o.accepted_source_states[i]!=o.sources[i]->state_sha256())throw std::runtime_error("prepared source digest differs after publication");
  if(o.accepted_arrival_state!=o.arrival->state_sha256())throw std::runtime_error("prepared arrival digest differs after publication");
  const J result={{"passed",true},{"whole_publication_allocations",0},{"hall_revision",o.hall->revision()},
                  {"arrival_next_step",o.arrival->next_step()},{"source_digests_match",true},{"arrival_digest_matches",true},
                  {"fixture","manufactured qual population; not a physical breakdown run"}};
  write_text(out/"publication-guard.json",result.dump(2)+"\n");
  return 0;
}
