#pragma once
// Owner bindings of coordinator v7 (continued): the adapters of the selected Hall/native
// pair, chosen by run_v7.py with the pair (PF_BINDING_PAIR_*): the native-side admission
// predicate of the native owner and the Hall-to-native motion-record adapter. The
// coordinator uses them only through bind::native_admission and bind::motion_adapter.
#if defined(PF_BINDING_PAIR_HALL_V3X_NATIVE_V6)
// Hall corrective v3 / v3.1 with native candidate006: candidate006 adapters/common_commit_v6.hpp
// (V3 contract types, birth-step admission) and hall-corrective-v3 native-v6-link/hall_v3_flows_adapter.hpp.
#include "common_commit_v6.hpp"
#include "hall_v3_flows_adapter.hpp"
namespace raikou::native_water::v7::bind {
namespace native_admission = ::raikou::native_water::common_commit_v6;
namespace motion_adapter = ::hall_v3_native_v6;
}  // namespace raikou::native_water::v7::bind
#else
#error "no adapter section for this Hall/native binding pair (bindings_adapters_v7.hpp)"
#endif
