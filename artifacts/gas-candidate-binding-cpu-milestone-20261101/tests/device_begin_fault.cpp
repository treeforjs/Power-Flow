// Test-only linker wrapper: explicitly injected before/after the real device
// initialization call. Normal operations always delegate to the actual method.
#include <cstdint>
#include <stdexcept>
int device_begin_failure_mode=0;
extern "C" void __real__ZN9pf_gas_v110DeviceLane10begin_stepEm(void*,std::uint64_t);
extern "C" void __wrap__ZN9pf_gas_v110DeviceLane10begin_stepEm(void* self,std::uint64_t step) {
  if(device_begin_failure_mode==1)throw std::runtime_error("manufactured device initialization failure before copy");
  __real__ZN9pf_gas_v110DeviceLane10begin_stepEm(self,step);
  if(device_begin_failure_mode==2)throw std::runtime_error("manufactured device initialization failure after copy");
}
