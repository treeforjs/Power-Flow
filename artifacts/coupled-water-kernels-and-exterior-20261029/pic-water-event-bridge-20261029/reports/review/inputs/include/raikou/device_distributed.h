#pragma once
#include "raikou/distributed.h"
#include "raikou/device.h"

namespace raikou::backend {
// Every bulk Krylov vector remains on the GPU. Only neighbor send/receive
// entries and reduction scalars cross the device/host boundary during solves.
DistributedLinearResult solve_distributed_cuda(const DistributedLinearOperator&,
                                               std::span<const double> local_rhs);
}
