#pragma once
// REVIEW CONTRACT ONLY. No definitions/bindings/runtime backend are supplied.
// CMake rejects enabling this bridge. Do not add a v2 reference-free overload.
#include "raikou/mfem_gap_pic.h"
#include <memory>
#include <span>

namespace raikou::mfem_backend::experimental_force {
enum class Mode { NativeCpu, QualificationShadow, QualifiedRuntime };
struct ForceSample {
    PicVector E_V_m{}, B_T{}, proper_velocity_m_s{};
    int status=0;
};
struct MidpointView {
    const GapPicOptions& options;
    const mfem::GridFunction& recovered_phi;
    std::span<const GapParticle> particles;
    double relative_midpoint_s,dt_s;
};
// An application-owned typed provider must lease exact unit-psi, native frame
// IDs/times/topology/axes/Az/rate arrays, fraction/voltage/origin and explicit
// layout/unit/native epochs. This type and provider remain UNIMPLEMENTED.
// Callback identity or rounded/reconstructed time is not a valid epoch.
class NativeMidpointState;
class QualifiedRuntimeApproval; // opaque; no issuer exists in this proposal

class PreparedForce {
public:
    virtual ~PreparedForce()=default; // discard uncommitted work safely
    virtual std::span<const ForceSample> output() const=0;
    // Whole-batch finite/status/shape/physical-input/domain/capacity/epoch and
    // unchanged-input/identity validation precedes velocity-only mutation.
    // Context/thread/resource requirements and safe cleanup are mandatory
    // without reference comparisons. No partial or implicit fallback commit.
    virtual void validate_and_commit(std::span<GapParticle> live)=0;
};
class RuntimeExecutor {
public:
    virtual ~RuntimeExecutor()=default;
    // Normal execution has NO independent reference argument or callback.
    // No implementation is supplied; the reviewed v2 API stays unchanged.
    virtual std::unique_ptr<PreparedForce> prepare(
        const MidpointView&,const NativeMidpointState&)=0;
};
using ReferenceProducer=std::function<std::vector<ForceSample>()>;

// Future bindings must be per GapPic owner/thread, explicit and scoped. The
// default is NativeCpu; a requested unsupported mode fails closed. Shared
// cuDSS-primary context, immutable leases, resource budgets and safe draining
// remain mandatory. No definitions/registration facility exist here.
Mode selected_mode(const GapPic&);
void execute_qualified_runtime(GapPic&,const MidpointView&,
                              std::span<GapParticle> live);
// The separate qualification wrapper generates fresh independent original
// E/B/u, then executes the runtime, applies ALL original componentwise/weighted
// reference gates, and only then permits transactional commit. The runtime
// path above cannot access or silently suppress this reference producer.
void execute_qualification_shadow(GapPic&,const MidpointView&,
    const ReferenceProducer&,std::span<GapParticle> live);
} // namespace raikou::mfem_backend::experimental_force
