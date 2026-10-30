#include "raikou/mfem_pic_full_time.h"
// This separate adapter is linked only with the exact qualified Hall owner.
// The core/legacy arena has no Hall runtime or field-solver dependency.
namespace hall_driver {
raikou::mfem_backend::ArenaJson motion_record(const MotionReceipt &);
}
namespace raikou::mfem_backend {
FullTimeMotionAdmission
bind_hall_total_field_motion(const hall_driver::MotionReceipt &receipt) {
  return FullTimeMotionAdmission(hall_driver::motion_record(receipt));
}
} // namespace raikou::mfem_backend
