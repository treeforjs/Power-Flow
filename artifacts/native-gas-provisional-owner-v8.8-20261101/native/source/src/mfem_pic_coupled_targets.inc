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
          !token.impl_->coupled_witness && token.impl_->chemistry && digest(proof),"wrong provisional token/proof");
  auto w=std::make_shared<CoupledTargetWitness>();w->retained=std::move(retained);w->check=check;w->proof=proof;
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
