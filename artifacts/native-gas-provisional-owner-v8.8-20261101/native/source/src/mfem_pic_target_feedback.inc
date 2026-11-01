// Host target overrides are facts from an external owner, never energy-closure
// or publication authority. Native owns and seals every input/output byte here.
std::string full_time_host_target_cells_sha256(std::span<const pfc::CellTargets> cells) {
  return chem8::block_root("native-v8/chemistry-cells",cells.data(),sizeof(pfc::CellTargets),cells.size(),64);
}
namespace full_time_detail {
std::string target_input_digest(const FullTimeHostTargetInputs &x) {
  Writer w;w.tag("native-v8.5/host-target-input-v1");
  const auto &e=x.exposure;const auto &b=e.common;
  for(auto v:{b.accepted_step,b.accepted_arena_revision,x.stage,x.native_rows,e.owner_stage_sequence})w.u64(v);
  for(double v:{b.begin_s,b.end_s,b.dt_s,e.frozen_physical_stage_s,e.collision_duration_s})w.f64(v);
  w.u64(std::uint64_t(e.slot));
  for(const auto *v:{&b.accepted_arena_sha256,&b.accepted_field_sha256,&b.mesh_sha256,&b.material_sha256,
      &b.common_clock_sha256,&e.frozen_field_sha256,&x.contract_sha256,&x.configuration_sha256,&x.registry_sha256,
      &x.blob_sha256,&x.grid_sha256,&x.native_rows_root,&x.census_sha256,&x.water_scope_binding_sha256,
      &x.frame_snapshot_sha256,&x.frame_history_sha256,&x.carriers_root,&x.targets_root})w.text(*v);
  w.u64(x.atlas.size());for(auto c:x.atlas)w.u64(c);
  return cn::sha256_hex(w.bytes);
}
std::string target_receipt_digest(const FullTimeHostTargetInputs &in,const FullTimeHostTargetDraft &out) {
  Writer w;w.tag("native-v8.5/host-target-receipt-v1");w.text(in.sha256);
  for(const auto *p:{&out.input_sha256,&out.owner_configuration_sha256,&out.owner_accepted_sha256,
      &out.owner_candidate_sha256,&out.owner_feedback_sha256,&out.before_targets_root,&out.after_targets_root})w.text(*p);
  return cn::sha256_hex(w.bytes);
}
std::shared_ptr<const FullTimeHostTargetInputs> make_target_inputs(
    const FullCollision &s,const chem8::Setup &setup,const chem8::StageInputs &in,
    const std::string &census,const std::string &scope_binding) {
  auto x=std::make_shared<FullTimeHostTargetInputs>();const auto &f=s.outer->frame;
  x->exposure=s.exposure;x->stage=s.exposure.slot==FullTimeSlot::before_transport?1:2;
  x->native_rows=f.native.particles.size();x->contract_sha256=setup.config.host_target_contract_sha256;
  x->configuration_sha256=setup.configuration_sha256;x->registry_sha256=setup.reg->sha256;
  x->blob_sha256=setup.reg->blob->sha256;x->grid_sha256=setup.grid.sha256;
  x->native_rows_root=cn::carrier_rows_digest(f.native.particles)->root;x->census_sha256=census;
  x->water_scope_binding_sha256=scope_binding;x->frame_snapshot_sha256=snapshot_digest(*f.snapshot);
  x->frame_history_sha256=f.snapshot->chain;
  x->carriers_root=chem8::block_root("native-v8/chemistry-carriers",in.carriers.data(),sizeof(pfc::CarrierIn),in.carriers.size(),1024);
  x->atlas=in.cell_atlas;x->cells=in.cells;x->targets_root=full_time_host_target_cells_sha256(x->cells);
  x->sha256=target_input_digest(*x);
  return std::shared_ptr<const FullTimeHostTargetInputs>(std::move(x));
}
void validate_target_draft(const chem8::Setup &setup,const FullTimeHostTargetInputs &in,
                           const FullTimeHostTargetDraft &out) {
  using K=FullTimeChemistryRefusal::Kind;
  auto check=[](bool yes,const char *why){if(!yes)throw FullTimeChemistryRefusal(K::target_fields,why);};
  check(in.sha256==target_input_digest(in) && in.targets_root==full_time_host_target_cells_sha256(in.cells),
        "host target input changed");
  check(in.contract_sha256==setup.config.host_target_contract_sha256 &&
        in.configuration_sha256==setup.configuration_sha256 && in.registry_sha256==setup.reg->sha256 &&
        in.blob_sha256==setup.reg->blob->sha256 && in.grid_sha256==setup.grid.sha256,"host target native configuration");
  for(const auto *p:{&out.input_sha256,&out.owner_configuration_sha256,&out.owner_accepted_sha256,
      &out.owner_candidate_sha256,&out.owner_feedback_sha256,&out.before_targets_root,&out.after_targets_root})
    check(digest(*p),"host target owner/receipt digest");
  check(out.input_sha256==in.sha256 && out.before_targets_root==in.targets_root && out.cells.size()==in.cells.size() &&
        out.after_targets_root==full_time_host_target_cells_sha256(out.cells),"host target byte binding");
  for(std::size_t i=0;i<in.cells.size();++i) {
    const auto &base=in.cells[i];auto restored=out.cells[i];
    check(base.n_targets<=PFC_MAX_TARGETS,"host target slot count");
    for(pfc_u32 k=0;k<base.n_targets;++k) {
      const auto sp=base.species[k];check(sp<setup.reg->states.size(),"host target species");
      if(setup.reg->charge[sp]==0 && setup.reg->target[sp] && base.n[k]>0) {
        check(std::isfinite(restored.s[k]) && restored.s[k]>=0,"host target thermal width");
        for(int a=0;a<3;++a){check(std::isfinite(restored.drift[k][a]),"host target drift");restored.drift[k][a]=base.drift[k][a];}
        restored.s[k]=base.s[k];
      }
    }
    check(std::memcmp(&restored,&base,sizeof base)==0,"host target changed density/species/layout/noneligible slot");
  }
}
struct TargetCallbackGuard {
  FullCollision &s;
  explicit TargetCallbackGuard(FullCollision &x):s(x){
    FullTimeAccess::enter_target_callback(s);
  }
  ~TargetCallbackGuard(){FullTimeAccess::leave_target_callback(s);}
};
std::shared_ptr<const FullTimeHostTargetReceipt> apply_target_provider(
    FullCollision &s,const chem8::Setup &setup,chem8::StageInputs &in,
    const std::string &census,const std::string &scope_binding,
    const std::shared_ptr<FullTimeHostTargetProvider> &provider) {
  const auto input=make_target_inputs(s,setup,in,census,scope_binding);
  FullTimeHostTargetDraft draft;
  {
    TargetCallbackGuard guard(s);
    draft=provider->apply(input);FullTimeAccess::live(s);
    provider->validate_current(input,draft);FullTimeAccess::live(s);
  }
  validate_target_draft(setup,*input,draft);
  auto receipt=std::make_shared<FullTimeHostTargetReceipt>();receipt->input=input;
  receipt->output=draft; // deep copy; no returned mutable storage survives here
  receipt->sha256=target_receipt_digest(*input,receipt->output);
  in.cells=receipt->output.cells;chem8::rehash_inputs(in,receipt->sha256);
  return std::shared_ptr<const FullTimeHostTargetReceipt>(std::move(receipt));
}
void validate_target_provider(FullCollision &s,const chem8::Setup &setup,
    const std::shared_ptr<FullTimeHostTargetProvider> &provider,const FullTimeHostTargetReceipt &receipt) {
  validate_target_draft(setup,*receipt.input,receipt.output);
  require(receipt.sha256==target_receipt_digest(*receipt.input,receipt.output),"host target receipt changed");
  TargetCallbackGuard guard(s);provider->validate_current(receipt.input,receipt.output);FullTimeAccess::live(s);
}
} // namespace full_time_detail
