#pragma once
#include <stdexcept>

namespace raikou::native_water::v8 {
// The source kernel's positive-field quadrature and microscopic launch have
// not been qualified below 1 V/m. Recording a bound while omitting its stream
// is not an admissible source approximation, even for a short timestep.
inline void require_emission_field_domain(double extraction_field_V_m) {
  if (extraction_field_V_m > 0 && extraction_field_V_m < 1)
    throw std::runtime_error("source/emission_field_domain: 0 < extraction field < 1 V/m requires a qualified low-field emission and launch model");
}
}  // namespace raikou::native_water::v8
