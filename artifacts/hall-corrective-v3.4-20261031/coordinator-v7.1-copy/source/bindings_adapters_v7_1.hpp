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
#elif defined(PF_BINDING_PAIR_HALL_V34_NATIVE_V7)
// Copy for Hall corrective v3.4 (resident owner, device coupled step) with native candidate007: the record adapter of
// the v3.2 pair (templated on the proposal type); StepInputs carry native's census v2 tables (packet keys, state
// indices) so Hall's device census needs no packet/state texts; the Model uses the resident backend.
#include "common_commit_v6.hpp"
#include "pair_hall_v32_native_v7.hpp"
namespace raikou::native_water::v7_1::bind {
namespace native_admission = ::raikou::native_water::common_commit_v6;
inline constexpr char pair_name[] = "hall-v3.4-resident+native-v7";
inline constexpr bool native_v7_forms = true;
inline constexpr double neutral_stream_dt_max_s = 1e-14;
inline const std::string& motion_record_sha256(const FieldProposal& p) { return p.motion().record_v7_sha256(); }
inline std::string iterate_record_sha256(const ::hall_driver::MotionReceipt& r) { return hall::MotionAssembly::record_v7(r); }
inline native::FullTimeMotionAdmission motion_admission(const FieldProposal& p) { return pair_v32_v7::admission(p); }
template <class Trial> std::span<const std::string> census_states(const Trial& t) { return t.census_v7().states; }
template <class Trial>
hall::StepInputs step_inputs(std::span<const hall::Newborn> n, std::span<const hall::SurfaceLoad> l, const Trial& t) {
  hall::StepInputs in{n, l, true, census_states(t)};
  const auto& c = t.census_v7();
  in.census_packet_keys = c.packet_keys; in.census_state_indices = c.state_indices;
  return in;
}
// The resident backend record of the device: descriptor of the resident kernels compiled for this device (a probe
// Model of the same mesh constructs one resident field), then the Model with that backend; the given solver is used by
// the probe only.
template <class Spec, class Trace, class Solver, class... A>
std::shared_ptr<const hall::Model> make_hall_model(double dt_max, Spec spec, Trace trace, Solver solver, A&&... a) {
  auto probe = std::make_shared<const hall::Model>(spec, trace, solver);
  const auto descriptor = hall::resident_descriptor_sha256(field_device(), *probe);
  return std::make_shared<const hall::Model>(spec, trace, std::make_shared<hall::ResidentBackend>(descriptor), std::forward<A>(a)..., dt_max,
                                             hall::MotionForms{true, neutral_stream_dt_max_s});
}
}  // namespace raikou::native_water::v7_1::bind
#elif defined(PF_BINDING_PAIR_HALL_V34_NATIVE_V71)
// Copy for Hall corrective v3.4 (resident owner, device coupled step) with native candidate007.1 (native v7.1 device lane,
// HALL-LINK-7.1.txt): native runs on Hall's session and forms (FullTimeDevice, attached by the driver); Hall reads the
// trial census in place (trial.census_v7_device()); the motion record carries the device addresses of the record arrays
// and of E, B, J (FullTimeMotionRecordV7::device) and the block-form field/current digests (pair_hall_v34_native_v71.hpp).
#include "common_commit_v6.hpp"
#include "pair_hall_v34_native_v71.hpp"
#include "raikou/mfem_pic_full_time_device.h"
#include <algorithm>
namespace raikou::native_water::v7_1::bind {
namespace native_admission = ::raikou::native_water::common_commit_v6;
inline constexpr char pair_name[] = "hall-v3.4-resident+native-v7.1-device-lane";
inline constexpr bool native_v7_forms = true;
inline constexpr double neutral_stream_dt_max_s = 1e-14;
// Native v7.1 lane of this process (set by the driver at genesis; counters in the run summary).
inline std::shared_ptr<native::FullTimeDevice>& native_lane() { static std::shared_ptr<native::FullTimeDevice> d; return d; }
inline native::FullTimeExecution& native_execution() { static native::FullTimeExecution e = native::FullTimeExecution::device; return e; }
// The device step of the field owner (record view of the accepted iterate; set by the driver at genesis).
inline hall::gpu::DeviceStep*& field_step() { static hall::gpu::DeviceStep* s = nullptr; return s; }
// Census v2 state table: native's registered states, ascending (API-007.txt 1.3), from the declared catalog; checked
// against the device census's states digest every step.
inline std::vector<std::string>& census_state_table() { static std::vector<std::string> s; return s; }
inline void set_registered_states(const std::vector<native::ArenaStateSpec>& states) {
  auto& t = census_state_table(); t.clear();
  for (const auto& s : states) t.push_back(s.id);
  std::sort(t.begin(), t.end());
  t.erase(std::unique(t.begin(), t.end()), t.end());
}
inline const std::string& motion_record_sha256(const FieldProposal& p) { return p.motion().record_v7_sha256(); }
inline std::string iterate_record_sha256(const ::hall_driver::MotionReceipt& r) { return r.record_v7_sha256(); }
inline native::FullTimeMotionAdmission motion_admission(const FieldProposal& p) {
  if (!field_step()) throw std::runtime_error("coordinator v7.1 copy: native v7.1 admission without the field owner's device step");
  return pair_v34_v71::admission(p, *field_step());
}
template <class Trial> std::span<const std::string> census_states(const Trial&) { return census_state_table(); }
template <class Trial>
hall::StepInputs step_inputs(std::span<const hall::Newborn> n, std::span<const hall::SurfaceLoad> l, const Trial& t) {
  const auto& d = t.census_v7_device();
  if (pf_v7::census_states_sha256(census_state_table()) != d.states_sha256)
    throw std::runtime_error("coordinator v7.1 copy: native census state table differs from the declared catalog");
  hall::StepInputs in{n, l, true, census_state_table()};
  in.device_rows = d.rows; in.device_keys = d.packet_keys; in.device_states = d.state_indices; in.device_count = d.carriers;
  in.packet_of = [&t](std::uint64_t id) { return t.carrier_packet(id); };
  return in;
}
// The midpoint digest the coordinator binds: native v7.1 recomputes the block-form digest from the device E/B arrays at
// admission (no host copy of the midpoint exists in the device lane).
inline std::string midpoint_digest(const FieldProposal& p) { return p.midpoint().total_field_sha256(); }
template <class Spec, class Trace, class Solver, class... A>
std::shared_ptr<const hall::Model> make_hall_model(double dt_max, Spec spec, Trace trace, Solver solver, A&&... a) {
  auto probe = std::make_shared<const hall::Model>(spec, trace, solver);
  const auto descriptor = hall::resident_descriptor_sha256(field_device(), *probe);
  return std::make_shared<const hall::Model>(spec, trace, std::make_shared<hall::ResidentBackend>(descriptor), std::forward<A>(a)..., dt_max,
                                             hall::MotionForms{true, neutral_stream_dt_max_s});
}
}  // namespace raikou::native_water::v7_1::bind
#else
#error "no adapter section for this Hall/native binding pair (bindings_adapters_v7_1.hpp)"
#endif
