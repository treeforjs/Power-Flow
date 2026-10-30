#pragma once
#include "snapshot_io.hpp"
#include <algorithm>
#include <span>
#include <iostream>
#include <iomanip>
namespace batch_compare {
namespace gp=gather_push;
using namespace raikou::mfem_backend;
using snapshot_io::vec;using snapshot_io::require;
inline constexpr double E_abs=1e-6,B_abs=1e-13,u_abs=1e-5,field_rel=5e-12,u_rel=1e-11,weighted_limit=1e-12;
struct Comparison {
    double max_E=0,max_B=0,max_u=0,weighted_u=0,weighted_energy=0;
    std::size_t failures=0;
};
inline Comparison compare(std::span<const gp::Outputs> ref,std::span<const gp::Outputs> actual,std::span<const GapParticle> particles) {
    require(ref.size()==actual.size() && ref.size()==particles.size(),"Comparison size mismatch");
    Comparison c;long double momentum_error[3]{},momentum_scale[3]{},energy_error=0,energy_scale=0;
    for(std::size_t i=0;i<ref.size();++i) {
        const auto& a=actual[i];const auto& r=ref[i];
        if(a.status!=0 || r.status!=0 || !gp::finite(r.E) || !gp::finite(r.B) || !gp::finite(r.u) ||
           !gp::finite(a.E) || !gp::finite(a.B) || !gp::finite(a.u) || !gp::finite(gp::norm2(a.u))) {
            ++c.failures;continue; // Explicitly reject NaN before comparisons.
        }
        for(int k=0;k<3;++k) {
            const double e=std::abs(a.E.v[k]-r.E.v[k]),b=std::abs(a.B.v[k]-r.B.v[k]),u=std::abs(a.u.v[k]-r.u.v[k]);
            c.max_E=std::max(c.max_E,e);c.max_B=std::max(c.max_B,b);c.max_u=std::max(c.max_u,u);
            if(e>E_abs+field_rel*std::abs(r.E.v[k]) || b>B_abs+field_rel*std::abs(r.B.v[k]) || u>u_abs+u_rel*std::abs(r.u.v[k]))++c.failures;
            momentum_error[k]+=static_cast<long double>(particles[i].weight)*(a.u.v[k]-r.u.v[k]);
            momentum_scale[k]+=static_cast<long double>(particles[i].weight)*std::abs(r.u.v[k]);
        }
        auto p=particles[i],q=p;p.proper_velocity_m_s=vec(r.u);q.proper_velocity_m_s=vec(a.u);
        const auto er=pic_kinetic_energy_J(p),ea=pic_kinetic_energy_J(q);
        if(!std::isfinite(ea) || !std::isfinite(er)) {++c.failures;continue;}
        energy_error+=static_cast<long double>(p.weight)*(ea-er);energy_scale+=static_cast<long double>(p.weight)*std::abs(er);
    }
    for(int k=0;k<3;++k)c.weighted_u=std::max(c.weighted_u,double(std::abs(momentum_error[k])/std::max(momentum_scale[k],1e-300L)));
    c.weighted_energy=double(std::abs(energy_error)/std::max(energy_scale,1e-300L));
    if(c.weighted_u>weighted_limit || c.weighted_energy>weighted_limit)++c.failures;
    return c;
}
inline void print_comparison(const Comparison& c,std::ostream& output=std::cout) {
    output<<std::setprecision(17)<<"{\"passed\":"<<(c.failures==0?"true":"false")<<",\"failures\":"<<c.failures
        <<",\"max_E_error_V_m\":"<<c.max_E<<",\"max_B_error_T\":"<<c.max_B<<",\"max_u_error_m_s\":"<<c.max_u
        <<",\"weighted_u_error_scaled_by_sum_abs\":"<<c.weighted_u<<",\"relative_weighted_energy_error\":"<<c.weighted_energy<<"}\n";
}
}
