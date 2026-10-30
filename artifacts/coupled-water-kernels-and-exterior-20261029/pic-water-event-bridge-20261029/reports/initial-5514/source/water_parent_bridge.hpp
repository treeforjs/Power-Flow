#pragma once
#include "collision_scheduler.hpp"
#include "raikou/h2o_conditional_parent.h"
#include <memory>

// Conditional parent-only event bridge. It owns no particles, neutral inventory,
// random generator, field state, heat, or accepted-step mutation.
namespace raikou::pic_water {
namespace pc = raikou::pic_collisions;
using Vec = pc::Vec;
constexpr double me = 9.1093837139e-31;
constexpr double threshold_eV = 12.621;
inline const std::string parent_model="parent_ground_equal_pair_CM_v1";
inline const std::string parent_channel="H2O:ion_yield:H2O+";
using Data = raikou::kinetics::StrictCollisionDataset;

inline long double gamma(const Vec& u) {
    long double u2=0;for(double v:u){pc::require(std::isfinite(v),"nonfinite proper velocity");u2+=static_cast<long double>(v)*v;}
    return std::sqrt(1+u2/(static_cast<long double>(pc::c)*pc::c));
}
inline long double kinetic(double mass,const Vec& u) {
    pc::positive(mass,"mass");long double u2=0;for(double v:u)u2+=static_cast<long double>(v)*v;
    return mass*u2/(gamma(u)+1);
}
// Boost from a target-rest frame to the lab. target_u is the target's LAB
// proper velocity; no velocity or time is borrowed from the SSA waiting clock.
inline Vec to_lab(const Vec& u,const Vec& target_u) {
    long double dot=0;for(int d=0;d<3;++d)dot+=static_cast<long double>(u[d])*target_u[d];
    const long double a=gamma(u)+dot/(static_cast<long double>(pc::c)*pc::c*(gamma(target_u)+1));
    Vec out{};for(int d=0;d<3;++d){out[d]=static_cast<double>(u[d]+a*target_u[d]);pc::require(std::isfinite(out[d]),"lab boost overflow");}
    return out;
}
inline Vec incident_direction(const pc::Packet& e,const pc::Packet& target,double energy_eV) {
    long double dot=0;for(int d=0;d<3;++d)dot+=static_cast<long double>(e.proper_velocity_m_s[d])*target.proper_velocity_m_s[d];
    const long double a=dot/(static_cast<long double>(pc::c)*pc::c*(gamma(target.proper_velocity_m_s)+1))-gamma(e.proper_velocity_m_s);
    std::array<long double,3> rest{};long double norm2=0;
    for(int d=0;d<3;++d){rest[d]=e.proper_velocity_m_s[d]+a*target.proper_velocity_m_s[d];norm2+=rest[d]*rest[d];}
    const long double x=energy_eV*static_cast<long double>(pc::qe)/(me*static_cast<long double>(pc::c)*pc::c);
    const long double expected=pc::c*std::sqrt(x*(x+2)),norm=std::sqrt(norm2);
    pc::require(norm>0&&std::isfinite(norm),"undefined incident direction");
    // If subtraction has lost direction information under a huge common boost,
    // refuse the event. Never normalize an inaccurate direction silently.
    pc::require(std::abs(norm-expected)<=256*std::numeric_limits<double>::epsilon()*expected,"target-rest direction precision unavailable");
    Vec direction{};for(int d=0;d<3;++d)direction[d]=static_cast<double>(rest[d]/norm);
    return direction;
}
struct Product {
    std::string state;
    int charge_e=0;
    double mass_kg=0,weight=0,reference_J=0;
    Vec proper_velocity_m_s{};
};
struct Products {
    std::array<Product,3> particles;
    pc::Remainder projectile_remainder,target_remainder;
    double energy_defect_J=0,energy_roundoff_limit_J=0;
    Vec momentum_defect_kg_m_s{};
    double momentum_roundoff_limit_kg_m_s=0;
    double target_frame_energy_eV=0,cross_section_m2=0;
};

class ParentBridge {
    std::shared_ptr<const Data> data_;
    std::string data_sha_,kernel_sha_;
public:
    // Host must bind these digests to exact bytes before construction. They are
    // propagated to every candidate, not interpreted as mutable file paths.
    ParentBridge(std::shared_ptr<const Data> data,std::string selected,
                 std::string data_sha,std::string kernel_sha)
        :data_(std::move(data)),data_sha_(std::move(data_sha)),kernel_sha_(std::move(kernel_sha)) {
        pc::require(data_&&selected==parent_model,"explicit conditional water parent model required");
        pc::require(pc::digest(data_sha_)&&pc::digest(kernel_sha_),"unbound parent model inputs");
        pc::require(data_->kernel_status(parent_channel,parent_model).at("status")=="conditional_available","parent model unavailable");
    }
    pc::Channel channel()const {
        auto data=data_;
        return {parent_channel,"electron","H2O",data_sha_,kernel_sha_,pc::Observable::ExclusiveEvent,
            [data](double E){
                const double sigma=data->observable_value_SI(parent_channel,E);
                if(sigma>0)data->require_MCC_event_kernel(parent_channel,parent_model,E);
                return sigma;
            }};
    }
    Products sample(const pc::Candidate& event,const pc::Packet& electron,
                    const pc::Packet& target,const Vec& isotropic_axis)const {
        pc::require(electron.state=="electron"&&target.state=="H2O"&&electron.mass_kg==me,"parent bridge species mismatch");
        pc::require(event.projectile_id==electron.id&&event.target_id==target.id&&electron.cell==target.cell&&event.cell==electron.cell,"stale parent identity/cell");
        pc::require(event.projectile_before==electron.weight&&event.target_before==target.weight,"stale parent weights");
        pc::require(event.channel_id==parent_channel&&event.data_sha256==data_sha_&&event.kernel_sha256==kernel_sha_,"foreign event data/kernel");
        const auto expected=pc::candidates({electron},{target},{{event.cell,event.volume_m3}},{channel()},1,1,1);
        pc::require(expected.size()==1,"no positive supported parent event");const auto& a=expected[0];
        pc::require(event.event_weight==a.event_weight&&event.energy_eV==a.energy_eV&&event.rate_s==a.rate_s&&event.flux_m_s==a.flux_m_s&&event.sigma_m2==a.sigma_m2,"stale or inconsistent rate candidate");
        pc::require(event.projectile_remainder.high==a.projectile_remainder.high&&event.projectile_remainder.low==a.projectile_remainder.low&&event.target_remainder.high==a.target_remainder.high&&event.target_remainder.low==a.target_remainder.low,"changed compensated remainder");
        const auto direction=incident_direction(electron,target,a.energy_eV);
        const auto rest=kinetics::conditional_h2o_parent_products(*data_,parent_model,a.energy_eV,a.event_weight,target.mass_kg,direction,isotropic_axis);
        Products out;out.projectile_remainder=a.projectile_remainder;out.target_remainder=a.target_remainder;
        out.target_frame_energy_eV=a.energy_eV;out.cross_section_m2=rest.cross_section_m2;
        out.particles={Product{"electron",-1,me,a.event_weight,0,{}},Product{"electron",-1,me,a.event_weight,0,{}},Product{"H2O+",1,rest.ion_mass_kg,a.event_weight,rest.chemical_reference_gain_J,{}}};
        const std::array<Vec,3> momenta{rest.primary_momentum_kg_m_s,rest.secondary_momentum_kg_m_s,rest.ion_momentum_kg_m_s};
        long double before=kinetic(electron.mass_kg,electron.proper_velocity_m_s)+kinetic(target.mass_kg,target.proper_velocity_m_s),after=rest.chemical_reference_gain_J;
        std::array<long double,3> dp{},abs_momentum{};
        for(int d=0;d<3;++d){dp[d]=-static_cast<long double>(electron.mass_kg)*electron.proper_velocity_m_s[d]-static_cast<long double>(target.mass_kg)*target.proper_velocity_m_s[d];abs_momentum[d]=std::abs(static_cast<long double>(electron.mass_kg)*electron.proper_velocity_m_s[d])+std::abs(static_cast<long double>(target.mass_kg)*target.proper_velocity_m_s[d]);}
        for(int i=0;i<3;++i){auto& product=out.particles[i];Vec u{};for(int d=0;d<3;++d)u[d]=momenta[i][d]/product.mass_kg;product.proper_velocity_m_s=to_lab(u,target.proper_velocity_m_s);after+=kinetic(product.mass_kg,product.proper_velocity_m_s);for(int d=0;d<3;++d){const long double p=static_cast<long double>(product.mass_kg)*product.proper_velocity_m_s[d];dp[d]+=p;abs_momentum[d]+=std::abs(p);}}
        const long double factor=256*std::numeric_limits<double>::epsilon();
        out.energy_defect_J=static_cast<double>((after-before)*a.event_weight);
        out.energy_roundoff_limit_J=static_cast<double>(factor*(std::abs(after)+std::abs(before))*a.event_weight);
        long double pscale=0;for(int d=0;d<3;++d){out.momentum_defect_kg_m_s[d]=static_cast<double>(dp[d]*a.event_weight);pscale+=abs_momentum[d];}
        out.momentum_roundoff_limit_kg_m_s=static_cast<double>(factor*pscale*a.event_weight);
        pc::require(std::isfinite(out.energy_defect_J)&&std::isfinite(out.energy_roundoff_limit_J)&&std::abs(out.energy_defect_J)<=out.energy_roundoff_limit_J,"parent lab energy closure");
        pc::require(std::isfinite(out.momentum_roundoff_limit_kg_m_s),"parent lab momentum bound overflow");
        for(double v:out.momentum_defect_kg_m_s)pc::require(std::isfinite(v)&&std::abs(v)<=out.momentum_roundoff_limit_kg_m_s,"parent lab momentum closure");
        return out;
    }
};
} // namespace raikou::pic_water
