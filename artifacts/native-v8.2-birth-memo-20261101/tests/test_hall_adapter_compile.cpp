// Native v7: compile-time traits of the motion admission and of the
// transitional receipt adapter against the sealed Hall 8dddec0 receipt and
// the Hall corrective v2 receipt types.
#include "raikou/mfem_pic_full_time.h"
#include "hall_receipt_v7.hpp"
#include "motion_receipt.hpp"
#include <type_traits>
namespace m = raikou::mfem_backend;
static_assert(!std::is_default_constructible_v<hall_driver::MotionReceipt>);
static_assert(std::is_same_v<decltype(m::v7_adapter::from_hall_receipt(
                                 std::declval<const hall_driver::MotionReceipt &>(),
                                 std::declval<const m::PrivateTransportTrial &>(),
                                 std::declval<std::span<const double>>(),
                                 std::declval<std::span<const double>>(), 1e-14)),
                             std::shared_ptr<m::v7_adapter::MotionV7>>);
static_assert(std::is_same_v<decltype(m::bind_total_field_motion_v7(
                                 std::declval<const m::FullTimeMotionRecordV7 &>(), {})),
                             m::FullTimeMotionAdmission>);
static_assert(!std::is_default_constructible_v<m::FullTimeMotionAdmission>);
static_assert(!std::is_copy_constructible_v<m::FullTimeMotionAdmission>);
static_assert(!std::is_constructible_v<m::FullTimeMotionAdmission, m::ArenaJson>);
static_assert(sizeof(m::FullTimeMotionPath) == sizeof(hall_driver::CarrierPath));
