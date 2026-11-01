namespace full_time_detail {
void stock_exact(Writer&w,const FullTimeExactValue&v){w.f64(v.value);w.text(v.digits);}
void stock_frame(Writer&w,const FullTimeChemistryStockFrame&f) {
  w.text(f.chemistry_state_sha256);w.u64(f.stages);w.boolean(f.pending);w.u64(f.reservoirs.size());
  for(const auto&r:f.reservoirs){w.u64(r.atlas_cell);w.u64(r.species);w.text(r.state);w.f64(r.rest_energy_eV);w.f64(r.reference_eV);
    stock_exact(w,r.count);stock_exact(w,r.kinetic_eV);for(const auto&p:r.momentum_eV_c)stock_exact(w,p);
    w.u64(r.parcels.size());for(const auto&p:r.parcels){w.f64(p.weight);for(double v:p.proper_velocity_m_s)w.f64(v);for(double v:p.position_m)w.f64(v);}}
  for(const auto&v:f.cumulative_internal_sink_eV)stock_exact(w,v);
  w.u64(f.cumulative_gas_transfer.size());for(const auto&v:f.cumulative_gas_transfer){w.u64(v.atlas_cell);w.u64(v.species);stock_exact(w,v.energy_eV);}
  w.u64(f.packet_internal.size());for(const auto&p:f.packet_internal){w.text(p.packet_id);w.f64(p.original_birth_weight);w.f64(p.energy_J);}
  for(const auto*v:{&f.cumulative_internal_born_J,&f.cumulative_internal_sunk_J,&f.cumulative_internal_released_J})stock_exact(w,*v);
}
std::string candidate_stocks_digest(const FullTimeCandidateStocks&v) {
  Writer w;w.tag(v.schema);auto e=v.exposure;codec_exposure(w,e);w.u64(v.parcels);w.u64(v.entries_and_parcels);w.u64(v.payload_bound_bytes);
  for(const auto*s:{&v.source_sha256,&v.configuration_sha256,&v.registry_sha256,&v.blob_sha256,
                    &v.before_state_prefix,&v.candidate_state_prefix,&v.admitted_export_sha256})w.text(*s);
  stock_frame(w,v.before);stock_frame(w,v.after);return cn::sha256_hex(w.bytes);
}
FullTimeExactValue stored_exact(const chem8::ESum&s) {
  auto v=chem8::exact_value(s);
  require(v.digits.size()<=4096,"stock exact digit bound");return v;
}
FullTimeChemistryStockFrame actual_stock_frame(const ArenaSnapshot&snap,const chem8::Setup&setup) {
  require(bool(snap.full_time_chemistry),"chemistry stock missing");const auto&st=*snap.full_time_chemistry;
  const auto&blob=*setup.reg->blob;FullTimeChemistryStockFrame f;
  f.chemistry_state_sha256=chem8::state_digest(st);f.stages=st.stages;f.pending=st.pending;
  f.reservoirs.reserve(st.reservoirs.size());
  for(const auto&[key,r]:st.reservoirs) {
    require(key.cell<setup.grid.cells() && setup.grid.vacuum[key.cell] && key.species<blob.n_species &&
            blob.view.species[key.species].charge==0,"reservoir stock key");
    // Validate once from actual persistent positive parcels, not event records.
    pfc::validate_reservoir(blob,{key.cell,key.species},r);
    FullTimeReservoirStock x;x.atlas_cell=key.cell;x.species=key.species;x.state=setup.reg->states[key.species];
    x.rest_energy_eV=blob.view.species[key.species].mass_eV;x.reference_eV=blob.view.species[key.species].ref_eV;
    x.count=stored_exact(r.count);x.kinetic_eV=stored_exact(r.kinetic);
    for(int a=0;a<3;++a)x.momentum_eV_c[a]=stored_exact(r.momentum[a]);
    x.parcels.reserve(r.parcels.size());for(const auto&p:r.parcels) {
      require(setup.grid.locate(p.x)==key.cell,"reservoir parcel cell differs from stored key");
      FullTimeReservoirParcel q;q.weight=p.weight;std::copy(p.u,p.u+3,q.proper_velocity_m_s.begin());
      std::copy(p.x,p.x+3,q.position_m.begin());x.parcels.push_back(q);
    }
    f.reservoirs.push_back(std::move(x));
  }
  static_assert(PFC_SINK_COUNT==6);
  for(std::size_t k=0;k<6;++k)f.cumulative_internal_sink_eV[k]=stored_exact(st.cumulative.c3.internal[k]);
  f.cumulative_gas_transfer.reserve(st.cumulative.gas_heating.size());
  for(const auto&[key,e]:st.cumulative.gas_heating)f.cumulative_gas_transfer.push_back({key.cell,key.species,stored_exact(e)});
  if(snap.full_time_v81) {
    const auto&flow=*snap.full_time_v81;f.packet_internal.reserve(flow.internal.size());
    for(const auto&[id,p]:flow.internal) {
      require(snap.packets.count(id) && std::isfinite(p.weight) && p.weight>0 && std::isfinite(p.energy_J) && p.energy_J>=0,"stored packet energy owner");
      f.packet_internal.push_back({id,p.weight,p.energy_J});
    }
    f.cumulative_internal_born_J=stored_exact(flow.internal_born_J);f.cumulative_internal_sunk_J=stored_exact(flow.internal_sunk_J);
    f.cumulative_internal_released_J=stored_exact(flow.internal_released_J);
  }else{
    const auto zero=stored_exact(chem8::ESum{});f.cumulative_internal_born_J=zero;f.cumulative_internal_sunk_J=zero;f.cumulative_internal_released_J=zero;
  }
  return f;
}
} // namespace full_time_detail
std::shared_ptr<const FullTimeCandidateStocks> FullTimeAccess::candidate_stocks(
    const FullCollision&s,const PreparedArena&t,std::size_t cap,std::size_t bytes) {
  require(inspectable(s,t) && t.impl_->chemistry && t.impl_->full_time_exposure,"candidate stocks require live chemistry token");
  require(!s.building && !s.outer->builder_active,"candidate stocks callback reentry");
  auto&owner=*s.outer->owner;
  if(owner.full_time_execution()!=FullTimeExecution::host || t.impl_->chemistry->device || !t.impl_->native_after.particles.host_present())
    throw FullTimeChemistryRefusal(FullTimeChemistryRefusal::Kind::target_fields,"actual stock export has no qualified device path");
  if(cap==0 || cap>full_time_absolute_carrier_cap)throw refusal(RK::carrier_cap,"candidate stock capacity range");
  if(t.impl_->stock_view) {
    if(t.impl_->stock_view->entries_and_parcels>cap)throw refusal(RK::carrier_cap,"candidate stock cached capacity");
    if(t.impl_->stock_view->payload_bound_bytes>bytes)throw refusal(RK::memory_budget,"candidate stock cached payload capacity");
    return t.impl_->stock_view;
  }
  const auto setup=std::static_pointer_cast<const chem8::Setup>(owner.chemistry_);
  std::uint64_t entries=0,parcels=0,bound=262144; // fixed fields and bounded 134-limb exact operands
  auto add=[&](std::uint64_t n){bound=checked_add(bound,n,"stock payload bound");};
  for(const auto*snap:{t.impl_->before.get(),t.impl_->after.get()}) {
    require(snap && snap->full_time_chemistry,"candidate stock state missing");const auto&st=*snap->full_time_chemistry;
    entries=checked_add(entries,st.reservoirs.size(),"stock entries");entries=checked_add(entries,st.cumulative.gas_heating.size(),"stock transfers");
    add(st.cumulative.gas_heating.size()*(sizeof(FullTimeCumulativeGasTransfer)+4096));
    for(const auto&[key,r]:st.reservoirs) {
      require(key.species<setup->reg->states.size(),"stock species bounds");parcels=checked_add(parcels,r.parcels.size(),"stock parcels");
      add(sizeof(FullTimeReservoirStock)+5*4096+setup->reg->states[key.species].size());add(r.parcels.size()*sizeof(FullTimeReservoirParcel));
    }
    if(snap->full_time_v81) {
      entries=checked_add(entries,snap->full_time_v81->internal.size(),"stock packet energy entries");
      for(const auto&[id,p]:snap->full_time_v81->internal){(void)p;add(sizeof(FullTimeStoredPacketEnergy)+id.size());}
    }
  }
  if(checked_add(entries,parcels,"stock combined capacity")>cap)throw refusal(RK::carrier_cap,"candidate stock entry/parcel capacity");
  if(bound>bytes)throw refusal(RK::memory_budget,"candidate stock payload capacity");
  auto v=std::make_shared<FullTimeCandidateStocks>();v->exposure=s.exposure;v->parcels=parcels;v->entries_and_parcels=entries+parcels;v->payload_bound_bytes=bound;
  v->source_sha256=PicArena::source_identity();v->configuration_sha256=setup->configuration_sha256;v->registry_sha256=setup->reg->sha256;
  v->blob_sha256=setup->reg->blob->sha256;v->before_state_prefix=prefix(*s.outer,s.outer->frame);v->candidate_state_prefix=t.impl_->candidate_prefix;
  v->admitted_export_sha256=t.impl_->chemistry->admitted_export->sha256;
  v->before=actual_stock_frame(*t.impl_->before,*setup);v->after=actual_stock_frame(*t.impl_->after,*setup);
  require(v->before.chemistry_state_sha256==t.impl_->chemistry->state_before_sha256 &&
          v->after.chemistry_state_sha256==t.impl_->chemistry->state_after_sha256,"actual stock/receipt chemistry roots");
  v->sha256=candidate_stocks_digest(*v);require(inspectable(s,t),"candidate stock changed during construction");
  t.impl_->stock_view=std::shared_ptr<const FullTimeCandidateStocks>(std::move(v));return t.impl_->stock_view;
}
