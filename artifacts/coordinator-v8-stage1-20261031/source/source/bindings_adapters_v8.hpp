#pragma once
// Owner bindings of coordinator v8 (continued): one section per admitted Hall/native pair, chosen by
// run_v8.py with the pair (PF_BINDING_PAIR_*). Every section provides the same names, the only
// pair-specific surface the step code uses (templates over the Hall proposal type, so the host Owner's
// JointProposal and the ResidentOwner's ResidentJointProposal bind through the same section):
//   native_admission            namespace of the native birth-step admission predicate (V3 types)
//   pair_name, native_v7_forms
//   motion_record_sha256(p)     record digest of an accepted Hall proposal that native admits
//   iterate_record_sha256(r)    record digest of one motion receipt (GPU push parity checks)
//   motion_admission(p)         native admission of the proposal's motion record
//   census_states(trial)        native's census v2 state table for Hall
//   step_inputs(newborns, loads, trial)   Hall StepInputs of a step with a birth scope
//   make_hall_model(dt_max, neutral_dt_max, spec, trace, solver, limits, declaration)
#include <memory>
#include <span>
#include <string>
#include <utility>
#if defined(PF_BINDING_PAIR_HALL_V32_NATIVE_V7) || defined(PF_BINDING_PAIR_HALL_V33_NATIVE_V71)
// Hall corrective v3.2 / v3.3 (MotionForms::native_v7) with native candidate007 (v7) / candidate007.1
// (v7.1, host mode: the v7 admission path; API-071.txt 3): the native adapters/common_commit_v6.hpp
// (unchanged from v6) and the coordinator-side record adapter pair_hall_native_v7_forms.hpp (motion record
// v7 with the A1 admission inputs). The device lane of native v7.1 (FullTimeDevice on Hall's session and
// forms) needs the Hall v3.4 device step; it is not attached by this section.
#include "common_commit_v6.hpp"
#include "pair_hall_native_v7_forms.hpp"
namespace raikou::native_water::v8::bind {
namespace native_admission = ::raikou::native_water::common_commit_v6;
#if defined(PF_BINDING_PAIR_HALL_V32_NATIVE_V7)
inline constexpr char pair_name[] = "hall-v3.2+native-v7";
#else
inline constexpr char pair_name[] = "hall-v3.3+native-v7.1";
#endif
inline constexpr bool native_v7_forms = true;
template <class P> const std::string& motion_record_sha256(const P& p) { return p.motion().record_v7_sha256(); }
inline std::string iterate_record_sha256(const ::hall_driver::MotionReceipt& r) { return hall::MotionAssembly::record_v7(r); }
template <class P> native::FullTimeMotionAdmission motion_admission(const P& p) { return pair_v7_forms::admission(p); }
template <class Trial> std::span<const std::string> census_states(const Trial& t) { return t.census_v7().states; }
template <class Trial>
hall::StepInputs step_inputs(std::span<const hall::Newborn> n, std::span<const hall::SurfaceLoad> l, const Trial& t) {
  return {n, l, true, census_states(t)};
}
// Hall model with the configured step cap and the resident-neutral interval bound shared by the Hall
// motion forms and the native configuration (FullTimeMotionPolicy).
template <class... A> std::shared_ptr<const hall::Model> make_hall_model(double dt_max, double neutral_dt_max, A&&... a) {
  return std::make_shared<const hall::Model>(std::forward<A>(a)..., dt_max, hall::MotionForms{true, neutral_dt_max});
}
}  // namespace raikou::native_water::v8::bind
#else
#error "no adapter section for this Hall/native binding pair (bindings_adapters_v8.hpp)"
#endif
