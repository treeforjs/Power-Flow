#pragma once
#include "raikou/chemistry.h"
#include "raikou/device.h"
#include <span>

namespace raikou::chemistry {
struct DeviceEvaluationBatch {
    std::vector<Evaluation> cells;
    backend::DeviceTransferMetrics transfers;
};
struct DeviceAdvanceBatch {
    std::vector<AdvanceResult> cells;
    backend::DeviceTransferMetrics transfers;
};
// One shared mechanism; either one environment for all cells, or one per cell.
// Temperatures/external rate inputs are prescribed during this local interval.
// Reaction progress, Jacobians, implicit Newton steps and adaptive control run
// on the GPU. Each cell is transactional if its timestep/Newton budget fails.
DeviceEvaluationBatch evaluate_cuda_batch(const Mechanism&,std::span<const State>,
                                          std::span<const Environment>,int device_id=-1);
DeviceAdvanceBatch advance_cuda_batch(const Mechanism&,std::span<const State>,
                                     std::span<const Environment>,double interval_s,
                                     const AdvanceOptions& = {});
}
