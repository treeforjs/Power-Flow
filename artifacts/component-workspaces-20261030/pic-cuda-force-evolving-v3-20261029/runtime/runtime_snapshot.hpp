#pragma once
#include "batch_snapshot.hpp"
namespace pic_batch {
// Typed live static state. It contains no independent reference or saved
// particle batch. Dynamic particles/phi/parameters are checked by evaluate.
struct RuntimeSnapshot {
    Parameters parameters;
    std::vector<gather_push::Cell> self_cells,native_cells;
    std::vector<double> phi,psi,az0,az1,rate0,rate1,xs,ys;
    void validate() const;
    static RuntimeSnapshot from_host(const HostSnapshot&);
};
}
