// Qualification of read-only source state projection, not a Hall capability.
void prepared_state_group(Results&R,const std::filesystem::path&out){
 using Clock=s7::CandidateThermalClock;axial::Mesh g(45,4,.009,1);
 auto clock=[](const s7::Owner&o,const eh::Epoch&e){return Clock{o.bound_thermal_configuration_signature(),e.t+e.dt,e.revision+1};};
 auto restore_check=[&](const eh::Patch&p,double T,double w,const se::Config&cfg,eh::Coupler&thermal,const s7::Owner&o,const std::string&hash,const std::string&tag){
  auto thermal_path=out/(tag+"-thermal.bin"),source_path=out/(tag+"-source.bin");thermal.save(thermal_path.string());o.save(source_path.string());
  auto hv=thermal_view(p,T);hv->load(thermal_path.string());auto restored=owner<VB>(p,T,*hv,cfg,w);restored->load(source_path.string());
  return restored->state_sha256()==hash&&restored->state_sha256_uncached()==hash;
 };
 // An empty source slice: thermal publication only, with actual prepared thermal state.
 {
  auto hv=thermal_view(patch(),1500);auto o=owner<VB>(patch(),1500,*hv);auto e=o->make_epoch(1e-14,"projection-field");auto ck=clock(*o,e);const auto before=o->state_sha256();
  auto prediction=o->unchanged_state_sha256_after_thermal_step(e,ck);const bool pure=o->state_sha256()==before&&o->thermal_state().revision==e.revision;
  mfem::ConstantCoefficient zero(0);auto heat=hv->prepare(e,zero,{},{{}});const bool thermal_binding=heat.state().time==ck.time_s&&heat.state().revision==ck.revision&&ck.configuration_signature==hv->config_signature();
  const bool commit=hv->commit(std::move(heat));
  R.record("prepared-state","unchanged slice prediction equals actual thermal publication and restored state",pure&&thermal_binding&&commit&&prediction==o->state_sha256()&&restore_check(patch(),1500,0,se::Config{4.3},*hv,*o,prediction,"projection-unchanged"));
  R.record("prepared-state","stale accepted epoch cannot project a second step",attempt([&]{o->unchanged_state_sha256_after_thermal_step(e,ck);}).threw);
 }
 // A live empty Trial exercises lifetime, exact field/config/clock checks and move/ABA behavior.
 {
  auto hv=thermal_view(patch(),1500);auto o=owner<VB>(patch(),1500,*hv);auto foreign=owner<VB>(patch(),1500,*hv);auto pf=field<VB>(g,*o,1,0);auto v=view<VB>(pf,g);auto e=o->make_epoch(1e-14,"projection-field");auto ck=clock(*o,e);
  auto trial=o->prepare(e,pf.f,v,{}, {},{{}});const auto before=o->state_sha256(),content=trial.content_sha256();
  auto prediction=o->candidate_state_sha256(trial,v,ck);const auto unchanged=o->unchanged_state_sha256_after_thermal_step(e,ck);
  R.record("prepared-state","empty trial and unchanged slice share the exact state encoding without mutation",prediction==unchanged&&o->state_sha256()==before&&trial.content_sha256()==content);
  J bad=J::object();auto wrong=ck;wrong.configuration_signature+="wrong";bad["configuration"]=attempt([&]{o->candidate_state_sha256(trial,v,wrong);}).threw;
  wrong=ck;++wrong.revision;bad["revision"]=attempt([&]{o->candidate_state_sha256(trial,v,wrong);}).threw;
  wrong=ck;wrong.time_s=std::nextafter(ck.time_s,INFINITY);bad["one_ulp_endpoint"]=attempt([&]{o->candidate_state_sha256(trial,v,wrong);}).threw;
  wrong=ck;wrong.time_s=NAN;bad["nonfinite_endpoint"]=attempt([&]{o->candidate_state_sha256(trial,v,wrong);}).threw;
  wrong=ck;wrong.revision=UINT64_MAX;bad["max_revision"]=attempt([&]{o->candidate_state_sha256(trial,v,wrong);}).threw;
  bad["foreign_owner"]=attempt([&]{foreign->candidate_state_sha256(trial,v,ck);}).threw;
  pf.a.E[0]=std::nextafter(pf.a.E[0],INFINITY);bad["changed_actual_field"]=attempt([&]{o->candidate_state_sha256(trial,v,ck);}).threw;pf.a.E[0]=0;
  auto changed=e;changed.dt=2e-14;auto changed_clock=clock(*o,changed);bad["interval_cap"]=attempt([&]{o->unchanged_state_sha256_after_thermal_step(changed,changed_clock);}).threw;
  changed=e;changed.source_hash+="wrong";bad["source_epoch"]=attempt([&]{o->unchanged_state_sha256_after_thermal_step(changed,ck);}).threw;
  changed=e;changed.material_hash+="wrong";bad["material_epoch"]=attempt([&]{o->unchanged_state_sha256_after_thermal_step(changed,ck);}).threw;
  bool all=true;for(const auto&x:bad)all=all&&x.get<bool>();
  R.record("prepared-state","foreign trial and invalid candidate clock/config/field/epoch refuse with accepted state intact",all&&o->state_sha256()==before&&foreign->state_sha256()==before&&trial.content_sha256()==content,bad);
  // Sweep every ordinary C++ allocation in this projection, until the first
  // complete success; failures of private hash caches do not change owner state.
  std::size_t refused=0;bool success=false,intact=true;std::string projected;
  for(std::ptrdiff_t k=0;k<4096;++k){bool failed=false;source_projection_fail_after=k;try{projected=o->candidate_state_sha256(trial,v,ck);}catch(...){failed=true;}source_projection_fail_after=-1;
   intact=intact&&o->state_sha256()==before&&trial.content_sha256()==content&&o->thermal_state().revision==e.revision;
   if(!failed){success=projected==prediction;break;}++refused;
  }
  R.record("prepared-state","ordinary C++ allocation failures in prepublication projection preserve accepted state and trial",intact&&success&&refused>0,{{"failure_sites",refused},{"eventual_success",success},{"scope","operator new/new[]; external C-library allocation is not injected"}});
  mfem::ConstantCoefficient zero(0);auto heat=hv->prepare(e,zero,trial.loads(),{{}});const bool admitted=o->preflight(trial,v)&&hv->preflight(heat);auto moved=std::move(trial);
  const bool moved_refused=attempt([&]{o->candidate_state_sha256(trial,v,ck);}).threw&&trial.source_candidates().empty()&&trial.source_handoffs().empty();
  // Thermal first then source: exactly the coordinator no-fail publication order.
  const bool thermal=hv->commit(std::move(heat));o->publish_prevalidated(std::move(moved));
  R.record("prepared-state","precomputed digest survives thermal-first source publication and moved-from trial refusal",admitted&&moved_refused&&thermal&&o->state_sha256()==prediction&&restore_check(patch(),1500,0,se::Config{4.3},*hv,*o,prediction,"projection-empty"));
 }
 // Nonempty exact source inventory + cumulative ledgers change in a real trial.
 {
  auto wp=water_patch(64,0x1p-50,1e13);auto hv=thermal_view(wp,1500);auto o=owner<VB>(wp,1500,*hv,se::Config{4.3},2.);
  auto pf=field<VB>(g,*o,1,0);auto v=view<VB>(pf,g);auto e=o->make_epoch(1e-14,"projection-field");auto ck=clock(*o,e);
  std::vector<std::vector<VB::Leaf>> leaves(1);for(int k=0;k<3;++k)leaves[0].push_back({2.,eh::Draw{.3,.4,.5,.6,.7,.1+.1*k}});
  auto trial=o->prepare(e,pf.f,v,{}, {},leaves);const auto before=o->state_sha256();const auto prediction=o->candidate_state_sha256(trial,v,ck);const auto untouched=o->unchanged_state_sha256_after_thermal_step(e,ck);
  const auto inventory_expected=o->inventory(0).debit(6.).positive_carriers();const bool has_packet=trial.packets().size()==3&&trial.inventory(0).positive_carriers()==inventory_expected;
  mfem::ConstantCoefficient zero(0);auto heat=hv->prepare(e,zero,trial.loads(),{{}});bool admitted=o->preflight(trial,v)&&hv->preflight(heat),pure=o->state_sha256()==before;
  bool thermal=hv->commit(std::move(heat));o->publish_prevalidated(std::move(trial));
  R.record("prepared-state","neutral debit and cumulative energy are predicted exactly and survive combined restore",admitted&&pure&&has_packet&&prediction!=untouched&&thermal&&o->state_sha256()==prediction&&o->inventory(0).positive_carriers()==inventory_expected&&restore_check(wp,1500,2.,se::Config{4.3},*hv,*o,prediction,"projection-neutral"));
 }
 // Real nonlegacy-mass microscopic layer; complete advance, persist, restore,
 // advance a second interval and refuse the unchanged-layer shortcut.
 {
  const char*file=std::getenv("PF_SOURCE_MASS_CONFIG");if(!file)throw std::runtime_error("actual registry source mass is required");std::ifstream in(file);J bound;in>>bound;se::Config cfg{4.3};cfg.electron_mass_kg=bound.at("mass_kg").get<double>();cfg.electron_mass_identity=bound.at("identity").get<std::string>();
  auto hv=thermal_view(patch(),1500);auto o=owner<VB>(patch(),1500,*hv,cfg);auto pf=field<VB>(g,*o,1,1e8);auto v=view<VB>(pf,g);auto e=o->make_epoch(1e-14,"projection-field");auto ck=clock(*o,e);
  auto p=prepared<VB>(*o,pf.f);auto born=o->begin_source(p,sample<VB>(p,*o,1,false),event(1));auto next=o->advance_birth(p,born,event(2));
  bool active=next.state().native.status==se::LayerStatus::Active;auto trial=o->prepare(e,pf.f,v,std::span<const VB::Candidate>(&next,1),{},{{}});const auto prediction=o->candidate_state_sha256(trial,v,ck);
  R.record("prepared-state","trial observational spans expose its exact immutable candidate and handoff set",trial.source_candidates().size()==1&&trial.source_candidates()[0].content_sha256()==next.content_sha256()&&trial.source_handoffs().empty()&&o->preflight(trial,v));
  const auto layer_base=o->state_sha256();bool layer_success=false,layer_intact=true;std::size_t layer_failures=0;std::string layer_projection;
  for(std::ptrdiff_t k=0;k<4096;++k){bool failed=false;source_projection_fail_after=k;try{layer_projection=o->candidate_state_sha256(trial,v,ck);}catch(...){failed=true;}source_projection_fail_after=-1;layer_intact=layer_intact&&o->state_sha256()==layer_base&&o->thermal_state().revision==e.revision;if(!failed){layer_success=layer_projection==prediction;break;}++layer_failures;}
  R.record("prepared-state","allocation faults in actual microscopic layer serialization refuse or return the exact full digest",layer_success&&layer_intact&&layer_failures>0,{{"failure_sites",layer_failures},{"eventual_success",layer_success}});

  mfem::ConstantCoefficient zero(0);auto heat=hv->prepare(e,zero,trial.loads(),{{}});bool admitted=o->preflight(trial,v)&&hv->preflight(heat),thermal=hv->commit(std::move(heat));o->publish_prevalidated(std::move(trial));
  R.record("prepared-state","active microscopic layer future digest equals accepted/restored declared-mass state",active&&admitted&&thermal&&o->state_sha256()==prediction&&restore_check(patch(),1500,0,cfg,*hv,*o,prediction,"projection-layer"));
  auto e2=o->make_epoch(1e-14,"projection-field-2");auto ck2=clock(*o,e2);const auto before=o->state_sha256();
  R.record("prepared-state","unchanged projection refuses a retained active microscopic layer",attempt([&]{o->unchanged_state_sha256_after_thermal_step(e2,ck2);}).threw&&o->state_sha256()==before);
  auto pf2=field<VB>(g,*o,2,1e8);auto v2=view<VB>(pf2,g);auto p2=prepared<VB>(*o,pf2.f);auto adv=o->advance_active(p2,42,event(3));
  std::vector<s7::Handoff>handoffs;if(adv.state().native.status==se::LayerStatus::PICHandedOff){auto canonical=canonical_electron<VB>(1);canonical.mass_kg=cfg.electron_mass_kg;handoffs.push_back(o->handoff(adv,v2,canonical));}
  auto tr2=o->prepare(e2,pf2.f,v2,std::span<const VB::Candidate>(&adv,1),handoffs,{{}});auto hash2=o->candidate_state_sha256(tr2,v2,ck2);auto heat2=hv->prepare(e2,zero,tr2.loads(),{{}});
  bool ready=o->preflight(tr2,v2)&&hv->preflight(heat2),heat_ok=hv->commit(std::move(heat2));o->publish_prevalidated(std::move(tr2));
  R.record("prepared-state","persisted layer real next-interval trial predicts its exact state",ready&&heat_ok&&o->state_sha256()==hash2&&restore_check(patch(),1500,0,cfg,*hv,*o,hash2,"projection-layer2"));
 }
 // A same-state empty publication invalidates previous trial generation even
 // while thermal time and all serialized source bytes remain unchanged.
 {
  auto hv=thermal_view(patch(),1500);auto o=owner<VB>(patch(),1500,*hv);auto pf=field<VB>(g,*o,1,0);auto v=view<VB>(pf,g);auto e=o->make_epoch(1e-14,"projection-field");auto ck=clock(*o,e);
  auto first=o->prepare(e,pf.f,v,{}, {},{{}}),second=o->prepare(e,pf.f,v,{}, {},{{}});const auto before=o->state_sha256();o->publish_prevalidated(std::move(first));
  R.record("prepared-state","same-state publication invalidates previous trial generation",o->state_sha256()==before&&attempt([&]{o->candidate_state_sha256(second,v,ck);}).threw);
 }
}
