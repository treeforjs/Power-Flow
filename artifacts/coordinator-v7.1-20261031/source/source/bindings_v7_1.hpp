#pragma once
// Owner bindings of coordinator v7.1: the one place in the coordinator sources that
// names the field (Hall) owner and the native owner; they are selected as a pair. Which owner build is compiled
// in is chosen by run_v7.py (build --hall NAME --native NAME), which sets the
// include paths, the owner objects and the libraries of that build and records the
// choice (PF_HALL_BINDING, PF_NATIVE_BINDING) in the binary and the build receipt.
//   Hall:   namespace hall_corrective_v3 from v3 (hall-corrective-v3-20261030) or v3.1
//           (hall-corrective-v3.1-20261031, configurable step cap); a later Hall owner with
//           the same namespace and API is one more entry in run_v7.py.
//   Native: raikou::mfem_backend from native candidate006 (native v6).
// The coordinator reaches the field-owner step cap only through field_owner_dt_max(), which
// detects the configurable-cap interface (Model::dt_max()); without it the cap is the v3 literal.
#include "raikou/mfem_pic_full_time.h"
#include "root_api.hpp"
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

namespace raikou::native_water::v7_1::bind {
namespace hall = ::hall_corrective_v3;
namespace native = ::raikou::mfem_backend;
inline constexpr char hall_binding[] = PF_HALL_BINDING;
inline constexpr char native_binding[] = PF_NATIVE_BINDING;
// Step bound of a Hall owner without a configurable cap (Hall corrective v3 literal).
inline constexpr double literal_dt_max_s = 1e-14;
template <class M> constexpr bool configurable_cap() { return requires(const M& m) { { m.dt_max() } -> std::convertible_to<double>; }; }
inline constexpr bool hall_configurable_cap = configurable_cap<hall::Model>();
// Admitted maximum step of the field owner's model.
inline double field_owner_dt_max(const hall::Model& m) {
  if constexpr (hall_configurable_cap) return m.dt_max();
  else return literal_dt_max_s;
}
// Step bound of the heat/source owner (not a swappable binding): source v6.1 uses the
// electrode-heat v3 kernel TrialFrame (pic-electrode-heat-feedback-20261029,
// mfem_electrode_heat_source_v3/trial_source_binding.hpp:18), which requires 0 < dt <= 1e-14.
inline constexpr double source_owner_dt_max_s = 1e-14;
// Largest step every coupled owner admits: the Hall model's cap and the source owner's bound.
// (Resident neutral carriers add the sealed neutral-stream bound 1e-14; the step controller
// applies it through its particle cap.)
inline double coupled_dt_max(const hall::Model& m) { return std::min(field_owner_dt_max(m), source_owner_dt_max_s); }
}  // namespace raikou::native_water::v7_1::bind
