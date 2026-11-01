#pragma once
// Owner bindings of coordinator v8: the one place in the coordinator sources that names the field
// (Hall) owner build and the native owner build; they are selected as a pair by run_v8.py
// (build LABEL --hall NAME --native NAME), which sets the include paths, owner objects and libraries
// and records the choice (PF_HALL_BINDING, PF_NATIVE_BINDING, PF_BINDING_PAIR_*) in the binary and the
// build receipt.
//   Hall:   namespace hall_corrective_v3 of Hall corrective v3.2 or v3.3 (both: host Owner and
//           device-resident ResidentOwner, configurable step cap, native v7 motion forms). Which of the
//           two owner types drives a run is chosen at launch (HallTraits, hall_owner_v8.hpp); a later Hall
//           version with the same protocol (v3.4: device step of the ResidentOwner) is one more
//           HALL_BINDINGS entry and one pair section in bindings_adapters_v8.hpp.
//   Native: raikou::mfem_backend of native candidate007 (v7) or candidate007.1 (v7.1, host mode).
// The heat/source owner is source owner v7 (electrode_heat_v7) for every pair; its trial interval
// bound is a run option (source_binding_v8.hpp).
// The coordinator reaches the field-owner step cap only through field_owner_dt_max(), which detects
// the configurable-cap interface (Model::dt_max()); without it the cap is the v3 literal.
#include "raikou/mfem_pic_full_time.h"
#include "root_api.hpp"
#include "resident_owner.hpp"
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

#ifndef PF_HALL_BINDING
#error "PF_HALL_BINDING (Hall owner build) must be set by the build"
#endif
#ifndef PF_NATIVE_BINDING
#error "PF_NATIVE_BINDING (native owner build) must be set by the build"
#endif

namespace raikou::native_water::v8::bind {
namespace hall = ::hall_corrective_v3;
namespace native = ::raikou::mfem_backend;
inline constexpr char hall_binding[] = PF_HALL_BINDING;
inline constexpr char native_binding[] = PF_NATIVE_BINDING;
// Step bound of a Hall owner without a configurable cap (Hall corrective v3 literal); also the
// default of every configurable bound (Hall dt_max, source trial interval, native neutral stream).
inline constexpr double literal_dt_max_s = 1e-14;
template <class M> constexpr bool configurable_cap() { return requires(const M& m) { { m.dt_max() } -> std::convertible_to<double>; }; }
inline constexpr bool hall_configurable_cap = configurable_cap<hall::Model>();
// Admitted maximum step of the field owner's model.
inline double field_owner_dt_max(const hall::Model& m) {
  if constexpr (hall_configurable_cap) return m.dt_max();
  else return literal_dt_max_s;
}
// Largest step every coupled owner admits: the Hall model's cap and the source owner's configured
// trial interval bound (source v7 OwnerOptions::interval.cap_s, <= 1e-12). The resident-neutral stream
// bound (Hall MotionForms / native FullTimeMotionPolicy) applies only while neutral carriers are
// resident; the step controller applies it through its particle cap (neutrals_need_particle_cap).
inline double coupled_dt_max(const hall::Model& m, double source_cap_s) { return std::min(field_owner_dt_max(m), source_cap_s); }
}  // namespace raikou::native_water::v8::bind
