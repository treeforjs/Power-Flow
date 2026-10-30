#include "batch_snapshot.hpp"
#include <set>

namespace pic_batch {
namespace io=snapshot_io;
namespace gp=gather_push;
gp::Tables HostSnapshot::tables() const {
    const auto& p=parameters;
    return {{self_cells.data(),p.self_nx,p.self_ny},{native_cells.data(),p.native_nx,p.native_ny},
        phi.data(),psi.data(),az0.data(),az1.data(),rate0.data(),rate1.data(),xs.data(),ys.data(),
        p.width,p.height,p.depth,p.short_z,p.fraction,p.voltage,p.dt};
}
void HostSnapshot::validate(bool require_real) const {
    const auto& p=parameters;
    io::require(p.magic==snapshot_magic && (p.flags&~1ULL)==0 && (!require_real || (p.flags&1)),"Unqualified or incompatible midpoint snapshot");
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
    io::require(!particles.empty() && particles.size()<=2000000 && particles.size()==inputs.size() && particles.size()==reference.size(),"Invalid batch sizes");
    std::set<std::uint64_t> ids;
    for(std::size_t i=0;i<particles.size();++i) {
        const auto& q=particles[i];const auto& v=inputs[i];const auto& ref=reference[i];
        io::require(ids.insert(q.id).second && gp::finite(io::vec(q.position_m)) && gp::finite(io::vec(q.proper_velocity_m_s)) &&
            std::isfinite(q.charge_C) && std::isfinite(q.mass_kg) && q.mass_kg>0 && std::isfinite(q.weight) && q.weight>0 &&
            std::isfinite(q.charge_C*q.weight) && std::isfinite(q.charge_C/q.mass_kg),"Invalid particle metadata");
        const gp::Inputs expected{q.position_m[0],q.position_m[1],q.position_m[2],io::vec(q.proper_velocity_m_s),q.charge_C/q.mass_kg};
        io::require(std::memcmp(&expected,&v,sizeof(v))==0,"Packed input differs from saved particle");
        io::require(v.x>=0 && v.x<=p.width && v.y>=0 && v.y<=p.height && gp::finite(gp::norm2(v.u)),"Particle outside domain or overflow");
        io::require(ref.status==0 && gp::finite(ref.E) && gp::finite(ref.B) && gp::finite(ref.u) && gp::finite(gp::norm2(ref.u)),"Invalid CPU reference");
    }
}
HostSnapshot HostSnapshot::load(const std::filesystem::path& path,bool require_real) {
    HostSnapshot s;const auto params=io::read<Parameters>(path/"parameters.bin");io::require(params.size()==1,"Expected one parameter record");s.parameters=params[0];
    s.self_cells=io::read<gp::Cell>(path/"self_cells.bin");s.native_cells=io::read<gp::Cell>(path/"native_cells.bin");
    s.phi=io::read<double>(path/"phi.bin");s.psi=io::read<double>(path/"psi.bin");
    s.az0=io::read<double>(path/"az0.bin");s.az1=io::read<double>(path/"az1.bin");s.rate0=io::read<double>(path/"rate0.bin");s.rate1=io::read<double>(path/"rate1.bin");
    s.xs=io::read<double>(path/"xs.bin");s.ys=io::read<double>(path/"ys.bin");
    s.particles=io::read<raikou::mfem_backend::GapParticle>(path/"particles.bin");s.inputs=io::read<gp::Inputs>(path/"inputs.bin");
    s.reference=io::read<gp::Outputs>(path/"reference.bin");s.validate(require_real);return s;
}
}
