#pragma once
#include "raikou/mfem_pic_full_time.h"
namespace native_gas_common { class CoupledTrial; }
namespace raikou::mfem_backend {
// Private native entry points for the pinned concrete CPU gas bridge. A hash,
// diagnostic token, or caller-supplied public callback cannot authorize this path.
class FullTimeCoupledAccess final {
  friend class ::native_gas_common::CoupledTrial;
  struct Continuation;
  using Check = bool (*)(const void *) noexcept;
  static void attach(PrivateFullTimeStep &, const std::string &contract,
                     const std::string &initial_gas, std::uint64_t gas_step,
                     const std::shared_ptr<const Continuation> &previous);
  static PreparedArena prepare(PrivateFullTimeStep &,CollisionScope &,std::shared_ptr<FullTimeHostTargetProvider>);
  static void authorize_stage(CollisionScope &,PreparedArena &,
      std::shared_ptr<const void>,Check,const std::string &proof);
  static std::shared_ptr<const Continuation> prepare_continuation(
      const PrivateFullTimeStep &,const PreparedArena &);
  static bool can_publish_final(const PrivateFullTimeStep &,const PreparedArena &) noexcept;
  static bool publish_final(PrivateFullTimeStep &,PreparedArena &) noexcept;
  static void authorize_final(PrivateFullTimeStep &,PreparedArena &,
      std::shared_ptr<const void>,Check,const std::string &proof);
};
}
