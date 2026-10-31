#pragma once
// Owner bindings of coordinator v7.1 (continued): one section per admitted Hall/native pair, chosen
// by run_v7_1.py with the pair (PF_BINDING_PAIR_*). Every section provides the same names, the only
// pair-specific surface the step code uses:
//   native_admission            namespace of the native birth-step admission predicate (V3 types)
//   pair_name, native_v7_forms
//   motion_record_sha256(p)     record digest of an accepted Hall proposal that native admits
//   iterate_record_sha256(r)    record digest of one motion receipt (GPU push parity checks)
//   motion_admission(p)         native admission of the proposal's motion record
//   census_states(trial)        native's census v2 state table for Hall (empty: Hall sorts the states)
//   step_inputs(newborns, loads, trial)   Hall StepInputs of a step with a birth scope
//   make_hall_model(dt_max, spec, trace, solver, limits, declaration)
#include <memory>
#include <span>
#include <string>
#include <utility>
#if defined(PF_BINDING_PAIR_HALL_V3X_NATIVE_V6)
// Hall corrective v3 / v3.1 with native candidate006: candidate006 adapters/common_commit_v6.hpp
// and hall-corrective-v3 native-v6-link/hall_v3_flows_adapter.hpp (binary v3 / flows record).
#include "common_commit_v6.hpp"
#include "hall_v3_flows_adapter.hpp"
namespace raikou::native_water::v7_1::bind {
namespace native_admission = ::raikou::native_water::common_commit_v6;
inline constexpr char pair_name[] = "hall-v3x+native-v6";
inline constexpr bool native_v7_forms = false;
inline const std::string& motion_record_sha256(const hall::JointProposal& p) { return p.motion().record_v3_sha256(); }
inline std::string iterate_record_sha256(const ::hall_driver::MotionReceipt& r) { return r.record_v3_sha256(); }
inline native::FullTimeMotionAdmission motion_admission(const hall::JointProposal& p) {
  return ::hall_v3_native_v6::bind_hall_v3_motion(p.motion());
}
template <class Trial> std::span<const std::string> census_states(const Trial&) { return {}; }
template <class Trial>
hall::StepInputs step_inputs(std::span<const hall::Newborn> n, std::span<const hall::SurfaceLoad> l, const Trial&) { return {n, l, true}; }
template <class... A> std::shared_ptr<const hall::Model> make_hall_model(double dt_max, A&&... a) {
  if constexpr (std::is_constructible_v<hall::Model, A..., double>) {
    return std::make_shared<const hall::Model>(std::forward<A>(a)..., dt_max);
  } else {
    if (dt_max != literal_dt_max_s) throw std::runtime_error(std::string("coordinator v7.1: the Hall binding ") + hall_binding +
                                                             " has no configurable step cap (only 1e-14 s)");
    return std::make_shared<const hall::Model>(std::forward<A>(a)...);
  }
}
}  // namespace raikou::native_water::v7_1::bind
#elif defined(PF_BINDING_PAIR_HALL_V32_NATIVE_V7)
// Hall corrective v3.2 (MotionForms::native_v7) with native candidate007: candidate007
// adapters/common_commit_v6.hpp (unchanged from v6) and the coordinator-side record adapter
// pair_hall_v32_native_v7.hpp (motion record v7 with the A1 admission inputs).
#include "common_commit_v6.hpp"
#include "pair_hall_v32_native_v7.hpp"
namespace raikou::native_water::v7_1::bind {
namespace native_admission = ::raikou::native_water::common_commit_v6;
inline constexpr char pair_name[] = "hall-v3.2+native-v7";
inline constexpr bool native_v7_forms = true;
// Resident-neutral interval bound shared by the Hall motion forms and the native configuration (default).
inline constexpr double neutral_stream_dt_max_s = 1e-14;
inline const std::string& motion_record_sha256(const hall::JointProposal& p) { return p.motion().record_v7_sha256(); }
inline std::string iterate_record_sha256(const ::hall_driver::MotionReceipt& r) { return hall::MotionAssembly::record_v7(r); }
inline native::FullTimeMotionAdmission motion_admission(const hall::JointProposal& p) { return pair_v32_v7::admission(p); }
template <class Trial> std::span<const std::string> census_states(const Trial& t) { return t.census_v7().states; }
template <class Trial>
hall::StepInputs step_inputs(std::span<const hall::Newborn> n, std::span<const hall::SurfaceLoad> l, const Trial& t) {
  return {n, l, true, census_states(t)};
}
template <class... A> std::shared_ptr<const hall::Model> make_hall_model(double dt_max, A&&... a) {
  return std::make_shared<const hall::Model>(std::forward<A>(a)..., dt_max, hall::MotionForms{true, neutral_stream_dt_max_s});
}
}  // namespace raikou::native_water::v7_1::bind
#else
#error "no adapter section for this Hall/native binding pair (bindings_adapters_v7_1.hpp)"
#endif
