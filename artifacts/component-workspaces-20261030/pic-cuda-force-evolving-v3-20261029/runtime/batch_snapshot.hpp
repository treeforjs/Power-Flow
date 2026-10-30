#pragma once
#include "snapshot_io.hpp"
#include <span>

namespace pic_batch {
inline constexpr std::uint64_t snapshot_magic=0x5049434241544348ULL;
struct Parameters {
    std::uint64_t magic=snapshot_magic,flags=0;
    int self_nx=0,self_ny=0,native_nx=0,native_ny=0;
    double width=0,height=0,depth=0,short_z=0,fraction=0,voltage=0,dt=0;
};
static_assert(sizeof(Parameters)==88);
// CPU-owned immutable-at-use snapshot. Static topology/unit field and native
// frame pair are explicitly uploaded once by CudaBatchedForce::set_snapshot;
// callers must provide a new snapshot when layout/unit/native frames change.
struct HostSnapshot {
    Parameters parameters;
    std::vector<gather_push::Cell> self_cells,native_cells;
    std::vector<double> phi,psi,az0,az1,rate0,rate1,xs,ys;
    std::vector<raikou::mfem_backend::GapParticle> particles;
    std::vector<gather_push::Inputs> inputs;
    std::vector<gather_push::Outputs> reference;
    gather_push::Tables tables() const;
    void validate(bool require_real_midpoint=true) const;
    static HostSnapshot load(const std::filesystem::path&,bool require_real_midpoint=true);
};
}
