#pragma once
#include "raikou/mfem_pic_full_time.h"
namespace hall_driver {
class MotionReceipt;
}
namespace raikou::mfem_backend {
// Binary motion v3 admission of a Hall corrective v2 receipt
// (hall-corrective-v2-20261030/source/motion_receipt.hpp). Separate archive;
// the core has no Hall dependency. The receipt must outlive the admission
// (its endpoint array is referenced, not copied).
FullTimeMotionAdmission
bind_hall_total_field_motion_v3(const hall_driver::MotionReceipt &);
} // namespace raikou::mfem_backend
