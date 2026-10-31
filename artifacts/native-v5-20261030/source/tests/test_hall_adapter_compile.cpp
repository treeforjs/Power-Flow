#include "raikou/mfem_pic_full_time.h"
#include "hall-motion-v2/motion_receipt.hpp"
#include <type_traits>
namespace m=raikou::mfem_backend;
static_assert(!std::is_default_constructible_v<hall_driver::MotionReceipt>);
static_assert(std::is_same_v<decltype(m::bind_hall_total_field_motion(std::declval<const hall_driver::MotionReceipt&>())),m::FullTimeMotionAdmission>);
static_assert(!std::is_copy_constructible_v<m::FullTimeMotionAdmission>);
static_assert(!std::is_constructible_v<m::FullTimeMotionAdmission,m::ArenaJson>);
