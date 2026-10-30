#include "runtime_snapshot.hpp"
namespace pic_batch {
namespace io=snapshot_io;
namespace gp=gather_push;
void RuntimeSnapshot::validate() const {
    const auto& p=parameters;
    io::require(p.magic==snapshot_magic && (p.flags&~1ULL)==0 && (p.flags&1),"Unqualified or incompatible midpoint snapshot");
    for(int n:{p.self_nx,p.self_ny,p.native_nx,p.native_ny})io::require(n>=1 && n<=4096,"Invalid snapshot dimensions");
    for(double x:{p.width,p.height,p.depth,p.dt})io::require(std::isfinite(x) && x>0,"Invalid geometry/dt");
    io::require(std::isfinite(p.short_z) && std::isfinite(p.voltage) && std::isfinite(p.fraction) && p.fraction>=0 && p.fraction<=1,"Invalid provider parameters");
    io::require(self_cells.size()==std::size_t(p.self_nx*p.self_ny) && native_cells.size()==std::size_t(p.native_nx*p.native_ny),"Cell count mismatch");
    io::require(!phi.empty() && phi.size()==psi.size() && !az0.empty() && az0.size()==az1.size() && az0.size()==rate0.size() && az0.size()==rate1.size(),"Nodal field count mismatch");
    for(const auto* data:{&phi,&psi,&az0,&az1,&rate0,&rate1,&xs,&ys})
        for(double value:*data)io::require(std::isfinite(value),"Nonfinite field or axis");
    io::require(xs.size()==std::size_t(p.native_nx+1) && ys.size()==std::size_t(p.native_ny+1),"Native axes count mismatch");
    for(const auto* axis:{&xs,&ys})for(std::size_t i=1;i<axis->size();++i)io::require((*axis)[i]>(*axis)[i-1],"Native axes not strictly increasing");
    for(const auto pair:{std::pair{&self_cells,phi.size()},std::pair{&native_cells,az0.size()}}) {
        for(const auto& cell:*pair.first) {
            for(int k=0;k<4;++k)io::require(cell.dof[k]>=0 && std::size_t(cell.dof[k])<pair.second &&
                std::isfinite(cell.x[k]) && std::isfinite(cell.y[k]),"Bad connectivity/geometry");
        }
    }
}
RuntimeSnapshot RuntimeSnapshot::from_host(const HostSnapshot& s) {
 return {s.parameters,s.self_cells,s.native_cells,s.phi,s.psi,s.az0,s.az1,s.rate0,s.rate1,s.xs,s.ys};
}
}
