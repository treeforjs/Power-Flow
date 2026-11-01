// CPU-only actual candidate census. It is computed once per prepared token;
// all external uses must retain/preflight the live collision scope and token.
namespace full_time_detail {
std::uint64_t population_payload_bytes(const FullTimeCandidatePopulation &x) {
  std::uint64_t n=sizeof x;
  auto add=[&](std::uint64_t v){n=checked_add(n,v,"candidate population payload");};
  add(x.states.size()*sizeof(ArenaStateSpec));
  add(x.packets.size()*sizeof(FullTimeCandidatePacket));
  add(x.carriers*sizeof(GapParticle));
  for(const auto *s:{&x.schema,&x.source_sha256,&x.configuration_sha256,&x.registry_sha256,
                    &x.before_state_prefix,&x.candidate_state_prefix,&x.candidate_snapshot_sha256,
                    &x.native_rows_root,&x.admitted_export_sha256,&x.sha256})add(s->size());
  const auto &e=x.exposure;const auto &b=e.common;
  for(const auto *s:{&b.accepted_arena_sha256,&b.accepted_field_sha256,&b.mesh_sha256,&b.material_sha256,
                    &b.common_clock_sha256,&e.frozen_field_sha256})add(s->size());
  for(const auto &s:x.states)add(s.id.size());
  for(const auto &p:x.packets)
    for(const auto *s:{&p.population.packet_id,&p.population.state,&p.origin,&p.last_event})add(s->size());
  return n;
}
std::string candidate_population_digest(const FullTimeCandidatePopulation &x) {
  Writer w;w.tag(x.schema);auto exposure=x.exposure;codec_exposure(w,exposure);
  w.u64(x.carriers);w.u64(x.payload_bytes);
  for(const auto *s:{&x.source_sha256,&x.configuration_sha256,&x.registry_sha256,&x.before_state_prefix,
                    &x.candidate_state_prefix,&x.candidate_snapshot_sha256,&x.native_rows_root,&x.admitted_export_sha256})w.text(*s);
  w.u64(x.states.size());
  for(const auto &s:x.states) {
    w.text(s.id);w.u64(std::uint64_t(s.hydrogen));w.u64(std::uint64_t(s.oxygen));
    w.u64(std::uint64_t(std::int64_t(s.charge_e)));w.f64(s.mass_kg);w.f64(s.internal_J);w.f64(s.reference_J);
  }
  w.u64(x.packets.size());
  for(const auto &p:x.packets) {
    export_packet(w,p.population);w.text(p.origin);w.text(p.last_event);
    w.boolean(p.has_carried_internal_energy);w.f64(p.internal_birth_weight);
  }
  return cn::sha256_hex(w.bytes);
}
} // namespace full_time_detail
void validate_full_time_candidate_population(const FullTimeCandidatePopulation &x) {
  using namespace full_time_detail;
  require(x.schema=="native-v8.6/chemistry-candidate-population-v1","candidate population schema");
  for(const auto *s:{&x.source_sha256,&x.configuration_sha256,&x.registry_sha256,&x.before_state_prefix,
                    &x.candidate_state_prefix,&x.candidate_snapshot_sha256,&x.native_rows_root,&x.admitted_export_sha256})
    require(digest(*s),"candidate population binding digest");
  require(x.carriers<=full_time_absolute_carrier_cap && x.packets.size()<=x.carriers,
          "candidate bounded shape");
  std::uint64_t actual=0;
  for(const auto &p:x.packets) {
    actual=checked_add(actual,p.population.carriers.size(),"candidate vector cardinality");
    require(actual<=x.carriers,"candidate declared cardinality underflow");
  }
  require(actual==x.carriers,"candidate declared cardinality mismatch");
  require(x.payload_bytes==population_payload_bytes(x),"candidate population payload count");
  require(x.sha256==candidate_population_digest(x),"candidate population content digest");
  std::map<std::string,const ArenaStateSpec *> states;
  for(const auto &s:x.states) {
    require(!s.id.empty() && (states.empty() || states.rbegin()->first<s.id) &&
            s.hydrogen>=0 && s.oxygen>=0 && s.charge_e>=-8 && s.charge_e<=8 &&
            std::isfinite(s.mass_kg) && s.mass_kg>0 && std::isfinite(s.internal_J) && s.internal_J>=0 &&
            std::isfinite(s.reference_J),"candidate state catalogue");
    states.emplace(s.id,&s);
  }
  std::set<std::uint64_t> seen;std::string previous;
  std::vector<GapParticle> ordered;ordered.reserve(x.carriers);
  for(const auto &entry:x.packets) {
    const auto &p=entry.population;
    require(!p.packet_id.empty() && (previous.empty() || previous<p.packet_id) && !p.count.zero() &&
            !p.carriers.empty() && states.count(p.state),"candidate packet identity/count");
    previous=p.packet_id;const auto &state=*states.at(p.state);ArenaCount count;std::uint64_t id=0;
    for(const auto &r:p.carriers) {
      native_particle_guard(r);
      require(r.id>id && seen.insert(r.id).second && r.mass_kg==state.mass_kg &&
              r.charge_C==state.charge_e*pic_e,"candidate carrier identity/catalogue");
      id=r.id;count=count.credited(r.weight);ordered.push_back(r);
    }
    require(count==p.count,"candidate exact packet count differs from positive carriers");
    require(std::isfinite(p.carried_internal_energy_J) && p.carried_internal_energy_J>=0 &&
            std::isfinite(entry.internal_birth_weight),"candidate packet-owned energy");
    if(entry.has_carried_internal_energy)require(entry.internal_birth_weight>0,"candidate internal birth weight");
    else require(entry.internal_birth_weight==0 && p.carried_internal_energy_J==0,"unexpected candidate internal owner");
  }
  require(seen.size()==x.carriers,"candidate carrier cardinality");
  std::sort(ordered.begin(),ordered.end(),[](const auto&a,const auto&b){return a.id<b.id;});
  require(cn::carrier_rows_digest(ordered)->root==x.native_rows_root,"candidate actual row commitment");
}
std::shared_ptr<const FullTimeCandidatePopulation> FullTimeAccess::candidate_population(
    const FullCollision &s,const PreparedArena &t,std::size_t cap,std::size_t bytes) {
  require(inspectable(s,t) && t.impl_->chemistry && t.impl_->full_time_exposure,
          "candidate population requires this live chemistry token");
  require(!s.building && !s.outer->builder_active,"candidate population callback reentry");
  auto &owner=*s.outer->owner;
  if(owner.full_time_execution()!=FullTimeExecution::host || t.impl_->chemistry->device ||
     !t.impl_->native_after.particles.host_present())
    throw FullTimeChemistryRefusal(FullTimeChemistryRefusal::Kind::target_fields,
                                  "candidate population has no qualified device view");
  const auto &snap=*t.impl_->after;const auto &rows=t.impl_->native_after.particles.host();
  if(cap==0 || cap>owner.carrier_cap_ || rows.size()>cap || snap.packets.size()>cap)
    throw refusal(RK::carrier_cap,"candidate population carrier/packet capacity");
  if(t.impl_->population_view) {
    if(t.impl_->population_view->payload_bytes>bytes)
      throw refusal(RK::memory_budget,"candidate population payload capacity");
    return t.impl_->population_view;
  }
  // Count every resident carrier and owned string before allocating the view.
  // payload_bytes is logical owned data, not an assertion of allocator/RSS cost.
  FullTimeCandidatePopulation header;
  header.exposure=s.exposure;header.carriers=rows.size();header.source_sha256=PicArena::source_identity();
  const auto setup=std::static_pointer_cast<const chem8::Setup>(owner.chemistry_);
  header.configuration_sha256=setup->configuration_sha256;header.registry_sha256=setup->reg->sha256;
  header.before_state_prefix=prefix(*s.outer,s.outer->frame);header.candidate_state_prefix=t.impl_->candidate_prefix;
  header.candidate_snapshot_sha256=snapshot_digest(snap);header.native_rows_root=t.impl_->native_after_digest->root;
  header.admitted_export_sha256=t.impl_->chemistry->admitted_export->sha256;header.sha256=std::string(64,'0');
  std::uint64_t required=population_payload_bytes(header),carriers=0;
  auto add=[&](std::uint64_t n){required=checked_add(required,n,"candidate population preallocation");};
  add(owner.states_.size()*sizeof(ArenaStateSpec));add(snap.packets.size()*sizeof(FullTimeCandidatePacket));
  for(const auto &[id,state]:owner.states_){(void)state;add(id.size());}
  for(const auto &[id,p]:snap.packets) {
    carriers=checked_add(carriers,p.native_ids.size(),"candidate population carrier inventory");
    require(carriers<=rows.size(),"candidate population duplicated inventory");
    for(const auto *v:{&id,&p.state,&p.origin,&p.last_event})add(v->size());
  }
  require(carriers==rows.size(),"candidate population missing inventory");
  if(required>bytes)throw refusal(RK::memory_budget,"candidate population payload capacity");
  auto view=std::make_shared<FullTimeCandidatePopulation>(std::move(header));
  view->payload_bytes=required;view->states.reserve(owner.states_.size());view->packets.reserve(snap.packets.size());
  for(const auto &[id,state]:owner.states_){(void)id;view->states.push_back(state);}
  std::vector<unsigned char> seen(rows.size(),0);
  for(const auto &[id,p]:snap.packets) {
    FullTimeCandidatePacket out;out.population.packet_id=id;out.population.state=p.state;
    out.population.count=p.count;out.origin=p.origin;out.last_event=p.last_event;
    out.population.carriers.reserve(p.native_ids.size());
    for(auto nid:p.native_ids) {
      const auto it=std::lower_bound(rows.begin(),rows.end(),nid,[](const auto&r,auto v){return r.id<v;});
      require(it!=rows.end() && it->id==nid,"candidate native ID missing");
      const auto index=std::size_t(it-rows.begin());require(!seen[index],"candidate native ID reused");seen[index]=1;
      out.population.carriers.push_back(*it);
    }
    std::sort(out.population.carriers.begin(),out.population.carriers.end(),[](const auto&a,const auto&b){return a.id<b.id;});
    if(snap.full_time_v81) {
      const auto it=snap.full_time_v81->internal.find(id);
      if(it!=snap.full_time_v81->internal.end()) {
        out.has_carried_internal_energy=true;out.internal_birth_weight=it->second.weight;
        out.population.carried_internal_energy_J=it->second.energy_J;
      }
    }
    view->packets.push_back(std::move(out));
  }
  view->sha256=candidate_population_digest(*view);
  validate_full_time_candidate_population(*view);
  require(inspectable(s,t),"candidate population changed during construction");
  t.impl_->population_view=std::shared_ptr<const FullTimeCandidatePopulation>(std::move(view));
  return t.impl_->population_view;
}
