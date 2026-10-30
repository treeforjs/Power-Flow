#pragma once
#include "batch_snapshot.hpp"
namespace pic_batch {
inline void validate_runtime_outputs(std::span<const gather_push::Outputs> output,
                                    std::span<const raikou::mfem_backend::GapParticle> particles) {
    namespace gp=gather_push;namespace io=snapshot_io;
    io::require(output.size()==particles.size(),"GPU runtime output shape mismatch");
    long double weighted_energy=0,weighted_u[3]{};
    for(std::size_t i=0;i<output.size();++i) {
        const auto& r=output[i];
        io::require(r.status==0 && gp::finite(r.E) && gp::finite(r.B) && gp::finite(r.u) && gp::finite(gp::norm2(r.u)),
                    "GPU runtime output finite/status guard failed");
        auto p=particles[i];p.proper_velocity_m_s=io::vec(r.u);
        const double energy=raikou::mfem_backend::pic_kinetic_energy_J(p);
        io::require(std::isfinite(energy),"GPU runtime output kinetic energy guard failed");
        weighted_energy+=static_cast<long double>(p.weight)*energy;
        for(int k=0;k<3;++k)weighted_u[k]+=static_cast<long double>(p.weight)*r.u.v[k];
    }
    io::require(std::isfinite(weighted_energy),"GPU runtime output weighted energy guard failed");
    for(auto value:weighted_u)io::require(std::isfinite(value),"GPU runtime output weighted momentum guard failed");
}
}
