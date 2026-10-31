// Declaration-only qualification of the v5 coordinator (production view: the
// test access structures are declared by the sources and defined nowhere).
// Without REJECT_CASE this translation unit must compile; with
// -DREJECT_CASE=k (k = 1..9) it must fail to compile for the stated reason.
#include "genesis_v5.hpp"
#include <type_traits>

namespace nw = raikou::native_water;
namespace v5 = raikou::native_water::v5;
namespace m = raikou::mfem_backend;
namespace hc = ::hall_corrective_v2;
using Trial = v5::ClosedParentTrialV5;

static_assert(!std::is_default_constructible_v<Trial>);
static_assert(!std::is_copy_constructible_v<Trial> && !std::is_copy_assignable_v<Trial>);
static_assert(std::is_nothrow_move_constructible_v<Trial> && std::is_nothrow_move_assignable_v<Trial>);
template <class T> concept ExposesNative = requires(T& t) { t.native_; };
template <class T> concept ExposesField = requires(T& t) { t.field_heat_; };
template <class T> concept ExposesRng = requires(T& t) { t.random_; };
template <class T> concept ExposesStep = requires(T& t) { t.private_step_; };
template <class T> concept ExposesInputs = requires(T& t) { t.inputs_; };
template <class T> concept ExposesHooks = requires(T& t) { t.publish_probe_; t.can_publish_fault_; t.retain_native_; };
static_assert(!ExposesNative<Trial> && !ExposesField<Trial> && !ExposesRng<Trial> &&
              !ExposesStep<Trial> && !ExposesInputs<Trial> && !ExposesHooks<Trial>);
static_assert(std::is_same_v<decltype(std::declval<const Trial&>().audit()), const nw::Json&>);
// The test access structures are incomplete in a production translation unit.
template <class T> concept Complete = requires { sizeof(T); };
static_assert(!Complete<v5::TrialTestAccess> && !Complete<v5::GenesisTestAccess>);

// T5: v5 never exposes the native rollback. No member, no ADL-visible free
// function taking a trial, and the trial holds no public native token.
template <class T> concept MemberRollback = requires(T& t) { t.try_rollback(); };
template <class T> concept AdlRollback = requires(T& t, m::PicArena& a) { try_rollback(a, t); };
template <class T> concept AdlRollbackV5 = requires(T& t, m::PicArena& a, nw::ParentSubflow& c, hc::Owner& o) {
  try_rollback_closed_parent_joint_v5(a, c, o, t);
};
static_assert(!MemberRollback<Trial> && !AdlRollback<Trial> && !AdlRollbackV5<Trial>);
template <class T> concept MemberRollbackOwners = requires(T& t) { t.try_rollback(); };
static_assert(!MemberRollbackOwners<v5::JointOwnersV5>);

// Geometry has no non-const accessor on any owner reachable by the coordinator.
static_assert(std::is_same_v<decltype(std::declval<m::PicArena&>().native_options()), const m::GapPicOptions&>);
static_assert(std::is_same_v<decltype(std::declval<hc::Owner&>().mesh()), const axial::Mesh&>);
static_assert(std::is_same_v<decltype(std::declval<hc::Model&>().mesh()), const axial::Mesh&>);
static_assert(std::is_same_v<decltype(std::declval<m::CollisionScope&>().native_options()), const m::GapPicOptions&>);
template <class T> concept MutableGeometry = requires(T& t) { t.native_options().width_m = 1.; };
static_assert(!MutableGeometry<m::PicArena> && !MutableGeometry<m::CollisionScope>);
template <class T> concept MutableMesh = requires(T& t) { t.mesh().side = 1.; };
static_assert(!MutableMesh<hc::Owner> && !MutableMesh<hc::Model>);

// Publication predicates are nonthrowing; preparation may throw.
static_assert(noexcept(v5::can_publish_closed_parent_joint_v5(std::declval<const m::PicArena&>(),
    std::declval<const nw::ParentSubflow&>(), std::declval<const hc::Owner&>(), std::declval<const Trial&>())));
static_assert(noexcept(v5::try_publish_closed_parent_joint_v5(std::declval<m::PicArena&>(),
    std::declval<nw::ParentSubflow&>(), std::declval<hc::Owner&>(), std::declval<Trial&>())));

#if defined(REJECT_CASE)
#if REJECT_CASE == 1
void reject(m::PicArena& a) { a.native_options().width_m = 1; }            // const native geometry
#elif REJECT_CASE == 2
void reject(hc::Owner& o) { o.mesh().side = 1; }                           // const Hall mesh
#elif REJECT_CASE == 3
void reject(const hc::Model& model) { model.mesh().dx = 1; }                // const Hall model mesh
#elif REJECT_CASE == 4
void reject(Trial& t) { (void)t.native_; }                                 // private native token
#elif REJECT_CASE == 5
Trial reject() { return Trial(); }                                          // private construction
#elif REJECT_CASE == 6
void reject(const Trial& t) { Trial copy(t); }                              // copy deleted
#elif REJECT_CASE == 7
void reject(m::PicArena& a, Trial& t) { v5::try_rollback_closed_parent_joint_v5(a, t); } // no rollback
#elif REJECT_CASE == 8
void reject(Trial& t) { v5::TrialTestAccess::set_publish_probe(t, nullptr); } // test access undefined
#elif REJECT_CASE == 9
void reject(const nw::v5::GenesisRequestV5& rq, const raikou::pic_water::ParentBridge& b) {
  (void)v5::JointGenesisV5::admit_with_hooks(rq, b, nullptr, nullptr, nullptr);  // genesis hooks private
}
#else
#error unknown rejection case
#endif
#endif
