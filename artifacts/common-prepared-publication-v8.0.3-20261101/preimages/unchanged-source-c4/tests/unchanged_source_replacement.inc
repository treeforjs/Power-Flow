// Include after the driver's genesis/load_inputs/TestAccess definitions.
// Dispatch as a separate test mode after finish_configuration(); no production run.
// This uses actual accepted owners/prepare, never edits a prepared source trial.
int unchanged_source_replacement_mode(const fs::path& out) {
  const auto in=load_inputs("qual");
  J attempts=J::array();
  for (double divisor : {1.0,16.0,256.0,4096.0}) {
    auto o=genesis<hc::Owner>(in,"qual");
    auto twin=genesis<hc::Owner>(in,"qual");
    const double dt=opt.dt_s/divisor;
    auto trial=v8::prepare_v8(o,dt);
    if (!v8::can_publish_v8(o,trial))
      throw std::runtime_error("replacement fixture requires initially publishable trial");
    std::size_t index=o.sources.size();
    const auto& source_trials=TA::source_trials(trial);
    for (std::size_t i=0;i<source_trials.size();++i)
      if (!source_trials[i]) { index=i;break; }
    attempts.push_back(J{{"dt_s",dt},{"source_count",o.sources.size()},
                         {"no_trial_found",index<o.sources.size()}});
    if (index==o.sources.size()) continue;
    if (o.sources[index].get()==twin.sources[index].get() ||
        o.sources[index]->state_sha256()!=twin.sources[index]->state_sha256() ||
        o.sources[index]->configuration_sha256()!=twin.sources[index]->configuration_sha256())
      throw std::runtime_error("replacement fixture lacks equal-content distinct source owners");
    const auto initial_state=o.sources[index]->state_sha256();
    const auto hall_before=o.hall->revision();
    const auto twin_before=twin.hall->revision();
    // Both halls and both source owners remain alive throughout observation.
    // Explicit swap-back also keeps each source bound to its own Hall at cleanup.
    auto original=std::move(o.sources[index]);
    o.sources[index]=std::move(twin.sources[index]);
    bool admissible=false,published=false,mismatch=false;
    J result;
    try {
      admissible=v8::can_publish_v8(o,trial);
      published=admissible && v8::try_publish_v8(o,trial);
      const auto actual=o.sources[index]->state_sha256();
      mismatch=published && actual!=o.accepted_source_states[index];
      result=J{{"source_index",index},{"dt_s",dt},{"source_trial_present",false},
               {"distinct_equal_content_source",true},{"can_publish_after_replacement",admissible},
               {"published",published},{"actual_source_state",actual},
               {"recorded_source_state",o.accepted_source_states[index]},
               {"original_accepted_source_state",initial_state},
               {"actual_vs_recorded_mismatch",mismatch},
               {"hall_advanced",o.hall->revision()==hall_before+1},
               {"foreign_thermal_unchanged",twin.hall->revision()==twin_before},
               {"attempts",attempts},{"fixture","equal-content foreign thermal owner; no physical run"}};
    } catch (...) {
      twin.sources[index]=std::move(o.sources[index]);
      o.sources[index]=std::move(original);
      throw;
    }
    twin.sources[index]=std::move(o.sources[index]);
    o.sources[index]=std::move(original);
    write_text(out/"unchanged-source-replacement.json",result.dump(2)+"\n");
    // Both baseline preimage and corrected fail-closed result are recognizable.
    // The parent harness should assert the expected branch for its pinned archive.
    if (admissible && (!published || !mismatch))
      throw std::runtime_error("replacement was admitted without the expected baseline preimage");
    return 0;
  }
  write_text(out/"unchanged-source-replacement-attempts.json",attempts.dump(2)+"\n");
  throw std::runtime_error("no unchanged source slice found; use benchmark interval policy without changing private trial state");
}
