// work_scales_v32.hpp -- Hall v3.2 work scales of the device motion fold: motion-work-gross-v1
// (../motion_work_gross_v1.hpp, the same code as the CPU motion assembly).
#pragma once
#include "../motion_work_gross_v1.hpp"

namespace hall_gpu::push {
using WorkScalesV32 = motion_work_v1::Scales;
inline double work_tolerance_v32(const WorkScalesV32& s) { return motion_work_v1::tolerance(s); }
}  // namespace hall_gpu::push
