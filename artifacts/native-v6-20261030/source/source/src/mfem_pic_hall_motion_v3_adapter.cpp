// Native v6 adapter: Hall corrective v2 MotionReceipt -> typed binary motion
// v3 admission. Compiled against hall-corrective-v2-20261030/source/
// motion_receipt.hpp (which provides record_v3_sha256()); linked only with
// that Hall motion provider. No JSON value is built.
#include "raikou/mfem_pic_hall_motion_v3.h"
#include "raikou/mfem_pic_canonical.h"
#include "motion_receipt.hpp"
namespace raikou::mfem_backend {
namespace {
// motion-array-binary64-v1 digest (Hall array_hash encoding).
std::string array_digest(const std::vector<double> &v) {
  canonical::Hasher h;
  h.w.tag("motion-array-binary64-v1");
  h.w.u64(v.size());
  for (double x : v) {
    h.w.f64(x);
    h.flush_if_large();
  }
  return canonical::hex(h.finish());
}
struct Storage {
  std::vector<FullTimeMotionPath> paths;
};
} // namespace
FullTimeMotionAdmission
bind_hall_total_field_motion_v3(const hall_driver::MotionReceipt &r) {
  auto storage = std::make_shared<Storage>();
  const auto paths = r.paths();
  storage->paths.resize(paths.size());
  for (std::size_t i = 0; i < paths.size(); ++i)
    storage->paths[i] = {paths[i].native_id, paths[i].field_start_m,
                         paths[i].field_displacement_m, paths[i].flight_s};
  const auto &b = r.binding();
  FullTimeMotionRecordV3 x;
  x.schema = motion_v3_schema;
  x.accepted_step = b.accepted_step;
  x.staged_arena_revision = b.staged_arena_revision;
  x.trial_stage = b.trial_stage;
  x.begin_s = b.begin_s;
  x.end_s = b.end_s;
  x.dt_s = b.dt_s;
  x.staged_arena_sha256 = b.staged_arena_sha256;
  x.trial_field_sha256 = b.trial_field_sha256;
  x.mesh_sha256 = b.mesh_sha256;
  x.material_sha256 = b.material_sha256;
  x.species_sha256 = b.species_sha256;
  x.common_clock_sha256 = b.common_clock_sha256;
  x.kinetic_policy_sha256 = b.kinetic_policy_sha256;
  const auto &u = r.budget();
  x.budget = {u.kinetic_change_J, u.field_work_J,  u.work_defect_J,
              u.continuity_C,     u.plasma_dt,     u.maximum_gyro,
              u.maximum_displacement_cells};
  x.census_sha256 = r.census_sha256();
  x.coverage = "all_native_carriers_no_sources_no_boundaries";
  x.current_A_m_sha256 = array_digest(r.current_A_m());
  x.endpoint_charge_C_sha256 = array_digest(r.endpoint_charge_C());
  x.initial_charge_C_sha256 = array_digest(r.initial_charge_C());
  x.endpoints = r.endpoints();
  const auto &k = r.kinetic_policy();
  x.kinetic_policy = {k.width_m,       k.gap_m,          k.short_z_m,
                      k.feed_z_m,      k.lower.left_m,   k.lower.right_m,
                      k.lower.bottom_m, k.lower.top_m,   k.upper.left_m,
                      k.upper.right_m, k.upper.bottom_m, k.upper.top_m};
  x.native_to_field_m = r.frame().native_to_field_m;
  x.neutral_motion_owner = "same_full_time_q0_drift";
  x.paths = storage->paths;
  x.record_sha256 = r.record_v3_sha256();
  return bind_total_field_motion_v3(x, storage);
}
} // namespace raikou::mfem_backend
