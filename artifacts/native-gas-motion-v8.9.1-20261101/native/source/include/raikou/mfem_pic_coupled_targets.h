#pragma once
#include "raikou/mfem_pic_full_time.h"
namespace native_gas_common { class CoupledTrial; }
namespace raikou::mfem_backend {
// Actual host populations of one genuine, still-private native motion.
// These records are evidence, not publication authority. The concrete common
// owner must retain and preflight the exact prepared motion and outer scope.
struct FullTimePreparedMotionSnapshot {
  std::string schema="native-v8.9/prepared-motion-inventory-v1";
  FullTimeMotionBinding binding;
  std::string source_sha256,configuration_sha256,registry_sha256,outer_step_sha256,
      motion_record_sha256,before_prefix,start_prefix,after_prefix,sha256;
  std::vector<ArenaRateEntry> before,start,after;
  std::vector<FullTimeStoredPacketEnergy> before_internal,start_internal,after_internal;
  std::shared_ptr<const FullTimeBirthReceipt> births;
  std::shared_ptr<const FullTimeSinkLedger> sinks;
  bool has_merge=false;
  std::uint64_t payload_bound_bytes=0;
};
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
  static std::shared_ptr<const FullTimePreparedMotionSnapshot> prepare_motion_snapshot(
      PrivateFullTimeStep &,PreparedPrivateMotion &,std::size_t carrier_cap,std::size_t byte_cap);
  static bool motion_live(const PrivateFullTimeStep &,const PreparedPrivateMotion &) noexcept;
  static bool can_inspect_motion(const PrivateFullTimeStep &,const PreparedPrivateMotion &,
      const FullTimePreparedMotionSnapshot &) noexcept;
  static std::shared_ptr<const FullTimeCandidatePopulation> final_population(
      const PrivateFullTimeStep &,const PreparedArena &);
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
