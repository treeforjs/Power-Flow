#include "kinetic_boundary.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <iomanip>
#include <limits>
#include <openssl/sha.h>
#include <sstream>
#include <stdexcept>
namespace hall_driver {
namespace {
void need(bool ok,const char*message){if(!ok)throw std::invalid_argument(message);}
void policy_guard(const KineticPolicy&p){
    for(double v:{p.width_m,p.gap_m,p.short_z_m,p.feed_z_m})need(std::isfinite(v),"nonfinite kinetic boundary");
    need(p.width_m>0&&p.gap_m>0&&p.feed_z_m>p.short_z_m,"invalid kinetic extents");
    for(const auto&f:{p.lower,p.upper})need(std::isfinite(f.left_m)&&std::isfinite(f.right_m)&&std::isfinite(f.bottom_m)&&std::isfinite(f.top_m)&&
        f.left_m>=0&&f.right_m<=p.width_m&&f.left_m<f.right_m&&f.bottom_m>=0&&f.top_m<=p.gap_m&&f.bottom_m<f.top_m,"invalid kinetic foil");
    need(p.lower.top_m<p.upper.bottom_m&&p.lower.left_m==p.upper.left_m&&p.lower.right_m==p.upper.right_m,"kinetic short bridge geometry");
}
}
KineticPolicy declared_kinetic_policy(const kinetic_native::GapPicOptions&o){
    need(o.finite_foils&&o.finite_axial_extent,"explicit finite-target kinetic policy required");
    KineticPolicy p{o.width_m,o.gap_m,o.short_z_m,o.feed_z_m,o.lower_foil,o.upper_foil};policy_guard(p);
    need(std::isfinite(o.depth_m)&&o.depth_m==o.feed_z_m-o.short_z_m,"kinetic depth must match target-end planes");return p;
}
std::string kinetic_policy_sha256(const KineticPolicy&p){policy_guard(p);std::string bytes=kinetic_boundary_api;
    bytes.push_back('\0');for(double v:{p.width_m,p.gap_m,p.short_z_m,p.feed_z_m,p.lower.left_m,p.lower.right_m,
        p.lower.bottom_m,p.lower.top_m,p.upper.left_m,p.upper.right_m,p.upper.bottom_m,p.upper.top_m}){
        auto bits=std::bit_cast<std::uint64_t>(v);for(int i=0;i<8;++i)bytes.push_back(char((bits>>(8*i))&255));}
    unsigned char digest[SHA256_DIGEST_LENGTH];SHA256(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size(),digest);
    std::ostringstream out;for(auto b:digest)out<<std::hex<<std::setfill('0')<<std::setw(2)<<unsigned(b);return out.str();
}
void require_kinetic_interior(const KineticPolicy&p,const kinetic_native::PicVector&x){policy_guard(p);
    need(x[0]>0&&x[0]<p.width_m&&x[1]>0&&x[1]<p.gap_m&&x[2]>p.short_z_m&&x[2]<p.feed_z_m,"particle requires native boundary owner");
    for(const auto&f:{p.lower,p.upper})need(!(x[0]>=f.left_m&&x[0]<=f.right_m&&x[1]>=f.bottom_m&&x[1]<=f.top_m),"particle in/on native metal");
}
KineticHit first_kinetic_hit(const KineticPolicy&p,const kinetic_native::PicVector&x,const kinetic_native::PicVector&D){
    require_kinetic_interior(p,x);for(double v:D)need(std::isfinite(v),"nonfinite kinetic displacement");KineticHit hit;
    auto take=[&](double f,KineticSink sink,int axis,double plane){if(f<0||f>hit.fraction)return;
        hit={sink,f,axis,plane,{}};hit.outward_from_vacuum_normal[axis]=D[axis]>0?1:-1;};
    for(int axis=0;axis<2;++axis)if(D[axis]!=0){double plane=D[axis]>0?(axis?p.gap_m:p.width_m):0;
        take((plane-x[axis])/D[axis],KineticSink::outer_xy_escape,axis,plane);}
    int body=0;for(const auto&foil:{p.lower,p.upper}){
        double enter=-std::numeric_limits<double>::infinity(),leave=std::numeric_limits<double>::infinity(),plane=0;
        int enter_axis=-1;bool possible=true;for(int axis=0;axis<2;++axis){double low=axis?foil.bottom_m:foil.left_m,high=axis?foil.top_m:foil.right_m;
            if(D[axis]==0){if(x[axis]<low||x[axis]>high)possible=false;continue;}
            double a=(low-x[axis])/D[axis],b=(high-x[axis])/D[axis],entry_plane=low;if(a>b){std::swap(a,b);entry_plane=high;}
            if(a>=enter){enter=a;enter_axis=axis;plane=entry_plane;}leave=std::min(leave,b);}
        if(possible&&leave>0&&enter>=0&&enter<=leave&&enter_axis>=0)take(enter,body?KineticSink::upper_foil:KineticSink::lower_foil,enter_axis,plane);++body;
    }
    if(D[2]!=0){double plane=D[2]>0?p.feed_z_m:p.short_z_m,f=(plane-x[2])/D[2];
        const double xx=x[0]+f*D[0],yy=x[1]+f*D[1];bool bridge=D[2]<0&&xx>=p.lower.left_m&&xx<=p.lower.right_m&&yy>=p.lower.top_m&&yy<=p.upper.bottom_m;
        take(f,bridge?KineticSink::short_bridge:KineticSink::axial_escape,2,plane);}
    return hit;
}
NeutralStreamProposal prepare_neutral_stream(const KineticPolicy&policy,const kinetic_native::GapParticle&p,double begin,double end,double dt){
    need(std::isfinite(begin)&&begin>=0&&std::isfinite(end)&&std::isfinite(dt)&&dt>0&&dt<=1e-14&&end>begin&&end==begin+dt,"neutral authoritative interval");
    need(p.id>0&&p.charge_C==0&&std::isfinite(p.weight)&&p.weight>0&&std::isfinite(p.mass_kg)&&p.mass_kg>0,"neutral stream identity/species");
    double u2=0;for(double v:p.proper_velocity_m_s){need(std::isfinite(v),"neutral momentum finite");u2+=v*v;}
    const double gamma=std::sqrt(1+u2/(kinetic_native::pic_c*kinetic_native::pic_c));need(std::isfinite(gamma),"neutral momentum outside resolved range");
    kinetic_native::PicVector D{};for(int d=0;d<3;++d)D[d]=dt*p.proper_velocity_m_s[d]/gamma;
    NeutralStreamProposal out;out.endpoint=p;out.hit=first_kinetic_hit(policy,p.position_m,D);out.begin_s=begin;out.end_s=end;out.dt_s=dt;
    out.flight_s=dt*out.hit.fraction;out.remaining_s=dt-out.flight_s;out.kinetic_policy_sha256=kinetic_policy_sha256(policy);
    for(int d=0;d<3;++d)out.endpoint.position_m[d]=p.position_m[d]+out.hit.fraction*D[d];
    if(out.hit.sink!=KineticSink::none)out.endpoint.position_m[out.hit.axis]=out.hit.plane_m;
    out.incident_kinetic_J=p.weight*kinetic_native::pic_kinetic_energy_J(p);need(std::isfinite(out.incident_kinetic_J),"neutral incident kinetic operand");return out;
}
} // namespace hall_driver
