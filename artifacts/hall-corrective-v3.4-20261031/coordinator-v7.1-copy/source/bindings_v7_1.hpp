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
// copy: the Hall v3.4 resident owner with native v7 (host admission) or native v7.1 (device lane)
#if defined(PF_BINDING_PAIR_HALL_V34_NATIVE_V7) || defined(PF_BINDING_PAIR_HALL_V34_NATIVE_V71)
#define PF_HALL_V34_RESIDENT 1
#else
#define PF_HALL_V34_RESIDENT 0
#endif
#if PF_HALL_V34_RESIDENT
#include "resident_owner.hpp"   // Hall v3.4: the resident owner (field state and the coupled step on the device)
#include "gpu/gpu_backends.hpp"
#else
#include "root_api.hpp"
#endif
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
// Field owner type of the pair (copy for Hall v3.4: the resident owner; otherwise the host Owner of v7.1).
#if PF_HALL_V34_RESIDENT
using FieldOwner = hall::ResidentOwner;
using FieldProposal = hall::ResidentJointProposal;
// The device shared by the resident field owner, its device step and the coordinator's GPU backends (set by the driver).
inline std::shared_ptr<hall::gpu::Device>& field_device() { static std::shared_ptr<hall::gpu::Device> d; return d; }
// Options of the resident owner: native v7 admits from host arrays (streaming digests); the native v7.1 device lane reads
// Hall's device arrays in place (no host record arrays, block-form field/current digests: HALL-LINK-7.1.txt 3-5).
inline hall::gpu::ResidentOptions field_owner_options() {
  hall::gpu::ResidentOptions r;
#if defined(PF_BINDING_PAIR_HALL_V34_NATIVE_V71)
  r.host_record_arrays = false;
  r.field_digest_form = hall::gpu::FieldDigestForm::blocks_v71;
#endif
  return r;
}
inline std::unique_ptr<FieldOwner> make_field_owner(std::shared_ptr<const hall::Model> m, axial::State seed, hall::ChargeLedger q,
                                                    const native::GapPicOptions& o, std::span<const native::ArenaRateEntry> census) {
  if (!field_device()) throw std::runtime_error("coordinator v7.1 copy: the Hall v3.4 resident owner needs the GPU device (field_device)");
  return std::make_unique<FieldOwner>(std::move(m), field_device(), std::move(seed), std::move(q), o, census, field_owner_options());
}
#else
using FieldOwner = hall::Owner;
using FieldProposal = hall::JointProposal;
inline std::unique_ptr<FieldOwner> make_field_owner(std::shared_ptr<const hall::Model> m, axial::State seed, hall::ChargeLedger q,
                                                    const native::GapPicOptions& o, std::span<const native::ArenaRateEntry> census) {
  return std::make_unique<FieldOwner>(std::move(m), std::move(seed), std::move(q), o, census);
}
#endif
}  // namespace raikou::native_water::v7_1::bind
