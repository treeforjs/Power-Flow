// Private integration entry points; no public arbitrary authority factory.
struct FullTimeCoupledAccess::Continuation {
  PicArena *owner=nullptr;std::weak_ptr<void> lifetime;
  std::shared_ptr<const ArenaSnapshot> accepted;
  std::shared_ptr<const FullTimeAcceptedBoundary> boundary;
  std::string contract,proof;
};
void FullTimeCoupledAccess::attach(PrivateFullTimeStep &step,const std::string &contract,
                                  const std::string &gas,std::uint64_t gas_step,
                                  const std::shared_ptr<const Continuation>&previous) {
  require(step.impl_ && step.impl_->outer,"empty coupled outer");
  auto &o=*step.impl_->outer;FullTimeAccess::live(o);
  const auto setup=std::static_pointer_cast<const chem8::Setup>(o.owner->chemistry_);
  require(setup && !contract.empty() && contract==setup->config.coupled_target_contract_sha256 &&
          digest(gas) && gas_step==o.common.accepted_step && o.phase==0 &&
          !o.first_reserved && o.coupled_contract.empty(),"coupled owner attachment/base mismatch");
  require(o.owner->full_time_execution()==FullTimeExecution::host &&
          setup->config.execution==FullTimeChemistryExecution::host,"coupled owner CPU path only");
  if(previous)
    require(!previous->lifetime.expired() && previous->owner==o.owner &&
            previous->accepted==o.before.snapshot && previous->boundary==o.boundary &&
            previous->contract==contract && digest(previous->proof),"foreign/stale coupled continuation");
  else
    require(o.boundary->authored_full_time_genesis && o.common.accepted_step==0,
            "coupled continuation requires an authenticated common accepted pair");
  // Allocate all content before committing attachment.
  std::string c=contract,g=gas;o.coupled_contract.swap(c);o.coupled_initial_gas.swap(g);
}
PreparedArena FullTimeCoupledAccess::prepare(PrivateFullTimeStep &step,CollisionScope &scope,
                                           std::shared_ptr<FullTimeHostTargetProvider> provider) {
  require(scope.impl_ && scope.impl_->scope && provider,"missing provisional scope/provider");
  require(step.impl_ && scope.impl_->scope->outer==step.impl_->outer,"foreign coupled outer/scope");
  auto token=FullTimeAccess::chemistry(*scope.impl_->scope,std::move(provider),true);
  token.impl_->target_scope=scope.impl_->scope;
  return token;
}
void FullTimeCoupledAccess::authorize_stage(CollisionScope &scope,PreparedArena &token,
    std::shared_ptr<const void> retained,Check check,const std::string &proof) {
  require(scope.impl_ && scope.impl_->scope,"empty coupled scope");
  auto &s=*scope.impl_->scope;
  require(FullTimeAccess::inspectable(s,token) && token.impl_->coupled_target_required &&
          !token.impl_->coupled_witness && token.impl_->chemistry && token.impl_->population_view && digest(proof),"wrong provisional token/proof");
  auto w=std::make_shared<CoupledTargetWitness>();w->retained=std::move(retained);w->check=check;w->proof=proof;
  w->population=token.impl_->population_view;
  w->native_candidate=token.impl_->candidate_prefix;
  require(w->valid(),"stale gas stage witness");
  token.impl_->coupled_witness=std::move(w);
}
void FullTimeCoupledAccess::authorize_final(PrivateFullTimeStep &step,PreparedArena &token,
    std::shared_ptr<const void> retained,Check check,const std::string &proof) {
  require(step.impl_ && step.impl_->outer && token.impl_,"empty coupled final token");
  auto &o=*step.impl_->outer;
  require(o.sealed && token.impl_->full_outer_identity==step.impl_->outer &&
          token.impl_->coupled_target_required && !token.impl_->private_identity &&
          !token.impl_->coupled_witness && o.coupled_stages[0] && o.coupled_stages[1] && digest(proof),
          "coupled final requires exact sealed outer and complete transcript");
  auto w=std::make_shared<CoupledTargetWitness>();w->retained=std::move(retained);w->check=check;w->proof=proof;
  require(w->valid(),"stale gas final witness");token.impl_->coupled_witness=std::move(w);
}

bool FullTimeCoupledAccess::can_publish_final(const PrivateFullTimeStep &step,const PreparedArena &token) noexcept {
  if(!step.impl_ || !step.impl_->outer || !token.impl_) return false;
  const auto &o=*step.impl_->outer;
  return !o.lifetime.expired() && o.owner && o.sealed &&
    token.impl_->full_outer_identity==step.impl_->outer && token.impl_->coupled_target_required &&
    o.owner->can_publish(token);
}
bool FullTimeCoupledAccess::publish_final(PrivateFullTimeStep &step,PreparedArena &token) noexcept {
  return can_publish_final(step,token) && step.impl_->outer->owner->try_publish(token);
}

std::shared_ptr<const FullTimeCoupledAccess::Continuation>
FullTimeCoupledAccess::prepare_continuation(const PrivateFullTimeStep&step,const PreparedArena&token) {
  require(can_publish_final(step,token),"continuation requires valid prepared common final token");
  const auto&o=*step.impl_->outer;auto p=std::make_shared<Continuation>();
  p->owner=o.owner;p->lifetime=o.lifetime;p->accepted=token.impl_->after;p->boundary=token.impl_->boundary_after;
  p->contract=o.coupled_contract;p->proof=token.impl_->coupled_witness->proof;return p;
}

bool FullTimeCoupledAccess::motion_live(const PrivateFullTimeStep&step,const PreparedPrivateMotion&token) noexcept {
  try {
    if(!step.impl_ || !step.impl_->outer || !token.impl_ || !token.impl_->trial) return false;
    const auto&o=*step.impl_->outer;const auto&t=*token.impl_;
    FullTimeAccess::live(o);
    return !o.coupled_contract.empty() && o.phase==1 && !o.builder_active &&
      t.trial->outer==step.impl_->outer && !t.consumed && !t.trial->consumed &&
      o.frame.snapshot==t.trial->fork_frame && t.after.snapshot && t.trial->base.snapshot;
  }catch(...){return false;}
}
bool FullTimeCoupledAccess::can_inspect_motion(const PrivateFullTimeStep&step,const PreparedPrivateMotion&token,
    const FullTimePreparedMotionSnapshot&view) noexcept {
  return motion_live(step,token) && token.impl_->motion_view.get()==&view &&
    step.impl_->outer->owner->full_time_execution()==FullTimeExecution::host &&
    step.impl_->outer->frame.native.particles.host_present() &&
    token.impl_->trial->base.native.particles.host_present() && token.impl_->after.native.particles.host_present();
}
namespace full_time_detail {
std::vector<FullTimeStoredPacketEnergy> motion_internal_stock(const ArenaSnapshot&snap) {
  std::vector<FullTimeStoredPacketEnergy> out;
  if(snap.full_time_v81) {
    out.reserve(snap.full_time_v81->internal.size());
    for(const auto&[id,p]:snap.full_time_v81->internal) {
      require(snap.packets.count(id) && std::isfinite(p.weight) && p.weight>0 &&
              std::isfinite(p.energy_J) && p.energy_J>=0,"motion packet energy owner");
      out.push_back({id,p.weight,p.energy_J});
    }
  }
  return out;
}
std::string motion_inventory_digest(const FullTimePreparedMotionSnapshot&v) {
  Writer w;w.tag(v.schema);auto common=v.binding.common;codec_common(w,common);
  for(const auto*s:{&v.binding.transport_trial_sha256,&v.binding.field_before_sha256,&v.binding.field_after_sha256,
      &v.binding.heat_after_sha256,&v.binding.source_after_sha256,&v.source_sha256,&v.configuration_sha256,
      &v.registry_sha256,&v.outer_step_sha256,&v.motion_record_sha256,&v.before_prefix,&v.start_prefix,&v.after_prefix})w.text(*s);
  w.boolean(v.has_merge);w.u64(v.payload_bound_bytes);
  for(const auto*rows:{&v.before,&v.start,&v.after}) {
    w.u64(rows->size());for(const auto&r:*rows) {
      w.text(r.packet_id);w.text(r.state);w.u64(r.native.id);w.vec(r.native.position_m);w.vec(r.native.proper_velocity_m_s);
      w.f64(r.native.charge_C);w.f64(r.native.mass_kg);w.f64(r.native.weight);
    }
  }
  for(const auto*stock:{&v.before_internal,&v.start_internal,&v.after_internal}) {
    w.u64(stock->size());for(const auto&p:*stock){w.text(p.packet_id);w.f64(p.original_birth_weight);w.f64(p.energy_J);}
  }
  w.boolean(bool(v.births));
  if(v.births) {
    const auto&b=*v.births;w.u64(b.accepted_step);w.u64(b.base);w.u64(b.next_id_after);w.u64(b.arrival_next_id_after);
    for(const auto*s:{&b.plan_sha256,&b.staged_token_sha256,&b.source_before_sha256,&b.source_after_sha256,
                     &b.births_root_sha256,&b.row_sha256,&b.internal_energy_J_digits})w.text(*s);
    w.f64(b.internal_energy_J);
    for(const auto*entries:{&b.entries,&b.wall_entries}){w.u64(entries->size());for(const auto&e:*entries){w.text(e.packet_id);w.text(e.state);w.u64(e.native_id);w.f64(e.weight);w.boolean(e.returned);}}
    for(const auto*streams:{&b.streams,&b.wall_credits}){w.u64(streams->size());for(const auto&s:*streams){w.u64(s.stream_id);w.text(s.state);codec_count(w,const_cast<ArenaCount&>(s.source_debit));}}
    w.u64(b.bridges.size());for(const auto&p:b.bridges){w.u64(p.native_id);w.vec(p.field_start_m);w.vec(p.field_displacement_m);w.f64(p.charge_C);}
  }
  w.boolean(bool(v.sinks));if(v.sinks)w.text(published_sinks_digest(*v.sinks));
  return cn::sha256_hex(w.bytes);
}
}
std::shared_ptr<const FullTimePreparedMotionSnapshot> FullTimeCoupledAccess::prepare_motion_snapshot(
    PrivateFullTimeStep&step,PreparedPrivateMotion&token,std::size_t cap,std::size_t bytes) {
  require(motion_live(step,token),"motion view requires exact live coupled outer/token");
  auto&o=*step.impl_->outer;auto&t=*token.impl_;auto&owner=*o.owner;
  const FullFrame*frames[]={&o.frame,&t.trial->base,&t.after};
  if(owner.full_time_execution()!=FullTimeExecution::host ||
     std::any_of(std::begin(frames),std::end(frames),[](const FullFrame*f){return !f->native.particles.host_present();}))
    throw FullTimeChemistryRefusal(FullTimeChemistryRefusal::Kind::target_fields,"motion inventory has no qualified device path");
  if(!cap || cap>owner.carrier_cap_ || !bytes || bytes>(std::size_t(1)<<30))
    throw refusal(RK::memory_budget,"motion view explicit carrier/byte capacity");
  for(const auto*f:frames)if(f->native.particles.size()>cap || f->snapshot->packets.size()>cap)
    throw refusal(RK::carrier_cap,"motion view per-frame carrier capacity");
  const auto birth=t.trial->births?t.trial->births->receipt:nullptr;
  // The same resource checks apply to fresh and cached views, even when a
  // smaller population accompanies many zero-count streams or sink records.
  if(birth && (birth->entries.size()>cap || birth->streams.size()>cap || birth->bridges.size()>cap))
    throw refusal(RK::carrier_cap,"motion birth receipt capacity");
  if(t.sinks && (t.sinks->records.size()>cap || t.sinks->exposed.size()>cap ||
     t.sinks->removed_packets.size()>cap || t.sinks->tallies.size()>cap || t.sinks->surfaces.size()>cap))
    throw refusal(RK::carrier_cap,"motion sink receipt capacity");
  if(t.motion_view) {
    if(t.motion_view->payload_bound_bytes>bytes)throw refusal(RK::memory_budget,"cached motion view byte capacity");
    return t.motion_view;
  }
  require(!t.trial->merge,"motion inventory merge needs a qualified common gas mapping");
  if(birth) {
    require(birth->wall_entries.empty() && birth->wall_credits.empty(),"motion inventory wall-response mapping unavailable");
    for(const auto&e:birth->entries)require(!e.returned,"motion inventory returned-birth mapping unavailable");
  }
  // Logical retained payload, not an RSS claim. Bound every row string and
  // retained receipt before allocating row vectors or copying any population.
  std::uint64_t payload=sizeof(FullTimePreparedMotionSnapshot)+4096;
  auto add=[&](std::uint64_t n){payload=checked_add(payload,n,"motion inventory payload overflow");
    if(payload>bytes)throw refusal(RK::memory_budget,"motion inventory payload capacity");};
  auto mult=[&](std::uint64_t n,std::uint64_t width){require(!width || n<=std::numeric_limits<std::uint64_t>::max()/width,"motion inventory payload product overflow");add(n*width);};
  add(0);
  for(const auto*f:frames) {
    mult(f->native.particles.size(),sizeof(ArenaRateEntry));
    for(const auto&[id,p]:f->snapshot->packets){(void)id;mult(p.native_ids.size(),p.id.size()+p.state.size());}
    if(f->snapshot->full_time_v81)for(const auto&[id,p]:f->snapshot->full_time_v81->internal){(void)p;add(sizeof(FullTimeStoredPacketEnergy)+id.size());}
  }
  constexpr std::uint64_t count_bytes=(ArenaCount::maximum_bits+7)/8;
  if(birth) {
    if(birth->entries.size()>cap || birth->streams.size()>cap || birth->bridges.size()>cap)throw refusal(RK::carrier_cap,"motion birth receipt capacity");
    add(sizeof(FullTimeBirthReceipt));for(const auto*s:{&birth->plan_sha256,&birth->staged_token_sha256,&birth->source_before_sha256,
      &birth->source_after_sha256,&birth->births_root_sha256,&birth->row_sha256,&birth->internal_energy_J_digits})add(s->size());
    for(const auto&e:birth->entries)add(sizeof e+e.packet_id.size()+e.state.size());
    for(const auto&s:birth->streams)add(sizeof s+s.state.size()+count_bytes);
    mult(birth->bridges.size(),sizeof(FullTimeBirthBridge));
  }
  if(t.sinks) {
    const auto&s=*t.sinks;
    if(s.records.size()>cap || s.exposed.size()>cap || s.removed_packets.size()>cap || s.tallies.size()>cap || s.surfaces.size()>cap)
      throw refusal(RK::carrier_cap,"motion sink receipt capacity");
    add(sizeof s);for(const auto*x:{&s.transport_trial_sha256,&s.records_root_sha256,&s.ledger_sha256,&s.row_sha256})add(x->size());
    mult(s.records.size(),sizeof(FullTimeSinkRecord));
    for(const auto&x:s.exposed)add(sizeof x+x.packet_id.size()+x.state.size());
    for(const auto&x:s.removed_packets)add(sizeof x+x.size());
    for(const auto&x:s.tallies)add(sizeof x+x.state.size()+2*count_bytes);
    mult(s.surfaces.size(),sizeof(FullTimeSinkSurfaceTotals)+4*count_bytes);
  }
  const auto setup=std::static_pointer_cast<const chem8::Setup>(owner.chemistry_);require(bool(setup),"motion inventory chemistry registry");
  auto v=std::make_shared<FullTimePreparedMotionSnapshot>();v->binding=t.binding;v->source_sha256=PicArena::source_identity();
  v->configuration_sha256=setup->configuration_sha256;v->registry_sha256=setup->reg->sha256;
  v->outer_step_sha256=o.outer_sha;v->motion_record_sha256=t.v7_sha256;
  require(digest(v->motion_record_sha256),"motion inventory requires typed native motion record");
  v->before_prefix=FullTimeAccess::prefix(o,o.frame);v->start_prefix=t.trial->prefix;v->after_prefix=FullTimeAccess::prefix(o,t.after);
  v->before=FullTimeAccess::population(o,o.frame,cap);v->start=FullTimeAccess::population(o,t.trial->base,cap);v->after=FullTimeAccess::population(o,t.after,cap);
  v->before_internal=motion_internal_stock(*o.frame.snapshot);v->start_internal=motion_internal_stock(*t.trial->base.snapshot);v->after_internal=motion_internal_stock(*t.after.snapshot);
  v->births=birth;v->sinks=t.sinks;v->has_merge=false;v->payload_bound_bytes=payload;v->sha256=motion_inventory_digest(*v);
  require(motion_live(step,token),"motion view token changed during construction");t.motion_view=v;return v;
}
std::shared_ptr<const FullTimeCandidatePopulation> FullTimeCoupledAccess::final_population(
    const PrivateFullTimeStep&step,const PreparedArena&token) {
  require(can_publish_final(step,token),"final population requires live common publication token");
  const auto&w=step.impl_->outer->coupled_stages[1];
  require(w && w->population && token.impl_->native_after_digest &&
          w->population->native_rows_root==token.impl_->native_after_digest->root &&
          w->population->carriers==token.impl_->native_after.particles.size(),
          "final population differs from actual native publication rows");
  return w->population;
}
