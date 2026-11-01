#pragma once
// Builds an unchanged Hall corrective v2 client (coordinator v5, read in place) against
// Hall corrective v3 for the Stage 1/v5 acceptance parity run: the client's namespace
// aliases (namespace hc = ::hall_corrective_v2) resolve to the v3 owner. Force-included
// with -include; the v3 public API is a superset of the v2 API.
// Angle includes: the v3 headers are found through the client's -isystem path, so
// the client's warning flags apply to its own code only (as for the v2 headers).
#pragma GCC system_header
#include <root_api.hpp>
#include <gpu/gpu_backends.hpp>
namespace hall_corrective_v2 = ::hall_corrective_v3;
