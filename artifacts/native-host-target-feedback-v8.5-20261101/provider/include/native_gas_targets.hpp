#pragma once
#include "raikou/mfem_pic_target_feedback.h"
#include "native_gas_stage.hpp"

// CPU reference preparation only. This never authorizes native/gas publication
// or opens the strict eventful target-energy-owner gate.
namespace native_gas_targets {
namespace n=raikou::mfem_backend;
namespace g=pf_gas_v1;
namespace a=native_gas_stage;
inline void need(bool x,const char *why){if(!x)throw std::invalid_argument(std::string("native_gas_targets: ")+why);}
static_assert(sizeof(pfc::CellTargets)==sizeof(g::CellTargetsAbi));
#define NGT_FIELD(x) static_assert(offsetof(pfc::CellTargets,x)==offsetof(g::CellTargetsAbi,x))
NGT_FIELD(n);NGT_FIELD(drift);NGT_FIELD(s);NGT_FIELD(species);NGT_FIELD(n_targets);NGT_FIELD(pad);
#undef NGT_FIELD
inline std::string sha(const std::string &s){return pfc::sha256_hex(s.data(),s.size());}
inline std::string contract(const n::FullTimeChemistryConfiguration &native,const n::GapPicOptions &options,
                            const g::Owner &gas) {
  auto base=native;base.host_target_feedback=false;base.host_target_contract_sha256.clear();
  const auto registry=n::full_time_chemistry_registry(base,options);
  return sha(n::ArenaJson{{"schema","native-gas-target-feedback-cpu-v1"},
    {"native_base_configuration",registry.configuration_sha256},{"registry",registry.registry_sha256},
    {"gas_configuration",gas.configuration_sha256()},{"model","neutral-Maxwellian-moment-correction"}}.dump());
}
inline std::string feedback_digest(const g::FeedbackReceipt &r) {
  return sha(n::ArenaJson{{"schema","native-gas-feedback-receipt-v1"},{"slots",r.slots},
    {"eligible",r.eligible},{"changed",r.changed},{"clamped",r.clamped},
    {"before",r.cells_before_sha256},{"after",r.cells_after_sha256}}.dump());
}
inline bool equal(const n::FullTimeHostTargetDraft &x,const n::FullTimeHostTargetDraft &y) {
  return x.input_sha256==y.input_sha256 && x.owner_configuration_sha256==y.owner_configuration_sha256 &&
    x.owner_accepted_sha256==y.owner_accepted_sha256 && x.owner_candidate_sha256==y.owner_candidate_sha256 &&
    x.owner_feedback_sha256==y.owner_feedback_sha256 && x.before_targets_root==y.before_targets_root &&
    x.after_targets_root==y.after_targets_root && x.cells.size()==y.cells.size() &&
    (x.cells.empty() || std::memcmp(x.cells.data(),y.cells.data(),x.cells.size()*sizeof(pfc::CellTargets))==0);
}
class Provider final:public n::FullTimeHostTargetProvider {
  // Shared ownership prevents a callback retaining a dangling gas owner. Its
  // private capability additionally rejects mutation/retry/equal-state owners.
  std::shared_ptr<g::Owner> gas_;
  a::Binding binding_;
  std::string contract_;
  g::Owner::CandidateBinding attempt_;
  std::shared_ptr<const n::FullTimeHostTargetInputs> input_;
  n::FullTimeHostTargetDraft output_;
  g::FeedbackReceipt feedback_;
  bool attempted_=false,ready_=false;
public:
  Provider(std::shared_ptr<g::Owner> gas,const n::FullTimeChemistryRegistry &registry,const pfc::Blob &blob,
           const n::FullTimeChemistryConfiguration &configuration,const n::GapPicOptions &options)
    :gas_(std::move(gas)),binding_(*(need(bool(gas_),"missing gas owner"),gas_),registry,blob,options),
     contract_(contract(configuration,options,*gas_)) {
    need(configuration.host_target_feedback && configuration.host_target_contract_sha256==contract_,"native contract mode");
    need(binding_.gas().feedback && binding_.gas().execution==g::Execution::host &&
      configuration.execution==n::FullTimeChemistryExecution::host,"explicit CPU gas feedback required");
    const auto actual=n::full_time_chemistry_registry(configuration,options);
    need(actual.configuration_sha256==registry.configuration_sha256 && actual.registry_sha256==registry.registry_sha256,
      "actual native configuration differs");
  }
  n::FullTimeHostTargetDraft apply(std::shared_ptr<const n::FullTimeHostTargetInputs> input) override {
    need(!attempted_,"provider is one-shot per native stage attempt");attempted_=true;
    need(bool(input),"missing authentic native input");input_=std::move(input);
    const auto &x=*input_;const auto &r=binding_.registry();
    need(x.configuration_sha256==r.configuration_sha256 && x.registry_sha256==r.registry_sha256 &&
      x.blob_sha256==r.blob_sha256 && x.grid_sha256==binding_.grid().sha256 && x.contract_sha256==contract_,
      "native/gas configuration binding");
    need(x.stage==(x.exposure.slot==n::FullTimeSlot::before_transport?1u:2u) &&
      x.atlas.size()==x.cells.size() && x.targets_root==n::full_time_host_target_cells_sha256(x.cells),"native stage/targets");
    std::uint32_t prev=0;bool first=true;
    for(auto cell:x.atlas){need(cell<binding_.grid().vacuum.size() && binding_.grid().vacuum[cell] &&
      (first || cell>prev),"native vacuum atlas ordering");first=false;prev=cell;}
    attempt_=gas_->candidate_binding();
    need(gas_->matches_candidate(attempt_,x.exposure.common.accepted_step),"gas/native step or candidate mismatch");
    output_.input_sha256=x.sha256;output_.owner_configuration_sha256=gas_->configuration_sha256();
    output_.owner_accepted_sha256=gas_->accepted_sha256();output_.owner_candidate_sha256=gas_->candidate_sha256();
    output_.before_targets_root=x.targets_root;
    std::vector<g::CellTargetsAbi> cells(x.cells.size());
    if(!cells.empty())std::memcpy(cells.data(),x.cells.data(),cells.size()*sizeof(g::CellTargetsAbi));
    feedback_=gas_->feedback(x.exposure.common.accepted_step,x.stage,x.atlas,cells);
    need(feedback_.clamped==0,"gas feedback reported clamping");
    output_.cells.resize(cells.size());
    if(!cells.empty())std::memcpy(output_.cells.data(),cells.data(),cells.size()*sizeof(g::CellTargetsAbi));
    output_.after_targets_root=n::full_time_host_target_cells_sha256(output_.cells);
    output_.owner_feedback_sha256=feedback_digest(feedback_);ready_=true;
    validate_current(input_,output_);return output_;
  }
  void validate_current(std::shared_ptr<const n::FullTimeHostTargetInputs> input,
                        const n::FullTimeHostTargetDraft &draft) const override {
    need(ready_ && input && input==input_ && equal(draft,output_),"changed feedback attempt/input/result");
    need(gas_->matches_candidate(attempt_,input_->exposure.common.accepted_step),"gas lifetime/step/revision changed");
    need(gas_->configuration_sha256()==output_.owner_configuration_sha256 && gas_->accepted_sha256()==output_.owner_accepted_sha256 &&
      gas_->candidate_sha256()==output_.owner_candidate_sha256,"gas feedback state changed");
  }
  const g::FeedbackReceipt &feedback() const {need(ready_,"feedback not prepared");return feedback_;}
};
} // namespace native_gas_targets
