#pragma once
#include "axial_fields.hpp"
#include "native_total_field_transport.hpp"

namespace axial {
struct Path {
    std::vector<foil::Weight> electric,magnetic;
    Path(const Mesh&g,V3 x,V3 D){
        V3 spacing{g.dx,g.dy,g.dz},origin{-g.side/2,-g.side/2,g.zmin};
        std::array<int,3>counts{g.nx,g.ny,g.nz};std::vector<double>cuts{0,1};
        for(int d=0;d<3;++d)if(D[d]!=0){double a=(x[d]-origin[d])/spacing[d],b=a+D[d]/spacing[d];
            hall::require(std::abs(b-a)<10000,"excessive axial path");
            for(int n=std::floor(std::min(a,b))+1;n<std::max(a,b);++n){double t=(n-a)/(b-a);if(t>0&&t<1)cuts.push_back(t);}}
        std::sort(cuts.begin(),cuts.end());
        for(std::size_t segment=1;segment<cuts.size();++segment){double midpoint=.5*(cuts[segment]+cuts[segment-1]),span=cuts[segment]-cuts[segment-1];
            for(double q:{-.5773502691896257645,.5773502691896257645}){double t=midpoint+.5*span*q,w=.5*span;
                std::array<int,3>cell{};V3 a{};
                for(int d=0;d<3;++d){double v=(x[d]+t*D[d]-origin[d])/spacing[d];
                    hall::require(v>=-1e-10&&v<=counts[d]+1e-10,"axial path outside mesh");v=std::clamp(v,0.,double(counts[d]));
                    cell[d]=std::min(int(std::floor(v)),counts[d]-1);a[d]=v-cell[d];}
                int i=cell[0],j=cell[1],k=cell[2];
                for(int r:{0,1})for(int s:{0,1}){
                    electric.push_back({g.ex(i,j+r,k+s),0,w*(r?a[1]:1-a[1])*(s?a[2]:1-a[2])});
                    electric.push_back({g.ey(i+r,j,k+s),1,w*(r?a[0]:1-a[0])*(s?a[2]:1-a[2])});
                    electric.push_back({g.ez(i+r,j+s,k),2,w*(r?a[0]:1-a[0])*(s?a[1]:1-a[1])});
                }
                for(int r:{0,1}){
                    magnetic.push_back({g.bx(i+r,j,k),0,w*(r?a[0]:1-a[0])});
                    magnetic.push_back({g.by(i,j+r,k),1,w*(r?a[1]:1-a[1])});
                    magnetic.push_back({g.bz(i,j,k+r),2,w*(r?a[2]:1-a[2])});
                }
            }
        }
    }
    V3 gather_E(const Vec&E)const{V3 v{};for(const auto&w:electric)v[w.axis]+=w.value*E[w.index];return v;}
    V3 gather_B(const Vec&B)const{V3 v{};for(const auto&w:magnetic)v[w.axis]+=w.value*B[w.index];return v;}
    void deposit(Vec&J,double charge,V3 D,double dt)const{for(const auto&w:electric)J[w.index]+=charge*D[w.axis]*w.value/dt;}
};
inline void point_deposit(const Mesh&g,Vec&Q,V3 x,double charge){
    V3 spacing{g.dx,g.dy,g.dz},origin{-g.side/2,-g.side/2,g.zmin},a{};
    std::array<int,3>counts{g.nx,g.ny,g.nz},cell{};
    for(int d=0;d<3;++d){double v=(x[d]-origin[d])/spacing[d];hall::require(v>=-1e-10&&v<=counts[d]+1e-10,"axial charge outside mesh");
        v=std::clamp(v,0.,double(counts[d]));cell[d]=std::min(int(std::floor(v)),counts[d]-1);a[d]=v-cell[d];}
    for(int k:{0,1})for(int j:{0,1})for(int i:{0,1})Q[g.node(cell[0]+i,cell[1]+j,cell[2]+k)]+=charge*(i?a[0]:1-a[0])*(j?a[1]:1-a[1])*(k?a[2]:1-a[2]);
}
inline foil::Hit intersection(const Mesh&g,V3 x,V3 D){
    foil::Hit hit;V3 lo{-g.side/2,-g.side/2,g.zmin},hi{g.side/2,g.side/2,g.zmin+g.nz*g.dz};
    for(int d=0;d<3;++d)if(D[d]!=0){double t=((D[d]>0?hi[d]:lo[d])-x[d])/D[d];if(t>=0&&t<=hit.fraction)hit={t,2};}
    for(int body=0;body<2;++body){V3 a{-.0025,body?.0025:-.0027,-g.depth/2},b{.0025,body?.0027:-.0025,g.depth/2};
        double enter=-1e100,leave=1e100;bool possible=true;
        for(int d=0;d<3;++d){if(D[d]==0){if(x[d]<a[d]||x[d]>b[d])possible=false;}else{
            double u=(a[d]-x[d])/D[d],v=(b[d]-x[d])/D[d];if(u>v)std::swap(u,v);enter=std::max(enter,u);leave=std::min(leave,v);}}
        if(possible&&leave>0&&enter>=0&&enter<=leave&&enter<=hit.fraction)hit={enter,body};
    }
    if(D[2]!=0){double t=(-g.depth/2-x[2])/D[2];if(t>=0&&t<=hit.fraction&&std::abs(x[0]+t*D[0])<=.0025&&std::abs(x[1]+t*D[1])<=.0027)hit={t,3};}
    return hit;
}
struct Push {GapParticle particle;V3 D;foil::Hit hit;double flight_s,gyro,displacement_cells;};
inline Push push(const Mesh&g,const GapParticle&p,const Vec&E,const Vec&B,double dt){
    auto trajectory=[&](double tau){V3 D=p.proper_velocity_m_s;for(double&v:D)v*=tau/foil::gamma(p.proper_velocity_m_s);
        V3 u{},Bg{};bool converged=false;
        for(int iteration=0;iteration<12;++iteration){auto hit=intersection(g,p.position_m,D);V3 gatherD=D;
            if(hit.surface>=0)for(double&v:gatherD)v*=hit.fraction;
            Path path(g,p.position_m,gatherD);auto Eg=path.gather_E(E);Bg=path.gather_B(B);
            auto step=foil::work_kick(p.proper_velocity_m_s,p.charge_C/p.mass_kg,Eg,Bg,tau);double diff=0;
            for(int d=0;d<3;++d)diff=std::max(diff,std::abs(step.second[d]-D[d]));u=step.first;D=step.second;
            if(diff<=1e-18){converged=true;break;}}
        hall::require(converged,"axial native work path did not converge");auto q=p;q.proper_velocity_m_s=u;
        for(int d=0;d<3;++d)q.position_m[d]+=D[d];
        return Push{q,D,intersection(g,p.position_m,D),tau,dt*std::abs(p.charge_C/p.mass_kg)*std::sqrt(hall::dot(Bg,Bg)),
            std::sqrt(hall::dot(D,D))/std::min({g.dx,g.dy,g.dz})};
    };
    auto result=trajectory(dt);
    hall::require(result.gyro<=.05&&result.displacement_cells<=.2,"axial full-interval gyro/displacement guard");
    if(result.hit.surface>=0){double lo=0,hi=dt;for(int i=0;i<48;++i){double t=.5*(lo+hi);auto trial=trajectory(t);
            if(trial.hit.surface>=0)hi=t;else lo=t;}
        result=trajectory(hi);hall::require(result.hit.surface>=0,"lost axial surface intersection");
        for(int d=0;d<3;++d){result.D[d]*=result.hit.fraction;result.particle.position_m[d]=p.position_m[d]+result.D[d];}}
    return result;
}
inline foil::NativeTransportTrial prepare_native_total_field_transport(
    const Mesh&g,const foil::NativeInterval&interval,std::span<const GapParticle>arena,
    const std::set<std::uint64_t>&layer_owned_ids,const Vec&E,const Vec&B){
    interval.validate();hall::require(int(E.size())==g.ne&&int(B.size())==g.nb&&finite(E)&&finite(B),"axial total-field binding");
    double dt=interval.dt_s;foil::NativeTransportTrial out;out.interval=interval;
    out.current_load_A_m.resize(g.ne);out.initial_charge_C.resize(g.nn);out.endpoint_charge_C.resize(g.nn);
    Vec gross(g.nn),omega(g.nn);std::set<std::uint64_t>ids;
    for(const auto&p:arena){hall::require(p.id!=0&&ids.insert(p.id).second&&std::isfinite(p.weight)&&p.weight>0&&
        std::isfinite(p.mass_kg)&&p.mass_kg>0&&std::isfinite(p.charge_C)&&p.charge_C!=0,"axial canonical particle invalid");
        for(const auto&v:{p.position_m,p.proper_velocity_m_s})for(double x:v)hall::require(std::isfinite(x),"nonfinite axial particle");
        point_deposit(g,omega,p.position_m,p.weight*p.charge_C*p.charge_C/(eps*p.mass_kg));
    }
    for(auto id:layer_owned_ids)hall::require(ids.count(id),"axial layer ID absent from canonical arena");
    for(int n=0;n<g.nn;++n)out.plasma_dt=std::max(out.plasma_dt,dt*std::sqrt(omega[n]/g.VN[n]));
    hall::require(out.plasma_dt<=.1,"axial full-arena plasma guard");
    for(const auto&p:arena){if(layer_owned_ids.count(p.id))continue;
        auto moved=push(g,p,E,B,dt);Path path(g,p.position_m,moved.D);
        path.deposit(out.current_load_A_m,p.weight*p.charge_C,moved.D,dt);
        point_deposit(g,out.initial_charge_C,p.position_m,p.weight*p.charge_C);
        point_deposit(g,out.endpoint_charge_C,moved.particle.position_m,p.weight*p.charge_C);
        point_deposit(g,gross,p.position_m,p.weight*std::abs(p.charge_C));
        out.kinetic_change_J+=p.weight*(pic_kinetic_energy_J(moved.particle)-pic_kinetic_energy_J(p));
        if(moved.hit.surface<0)out.updates.push_back(moved.particle);
        else out.boundaries.push_back({moved.particle,moved.hit.surface,p.weight*pic_kinetic_energy_J(moved.particle),moved.flight_s});
    }
    out.field_work_J=dt*dot(E,out.current_load_A_m);out.work_defect_J=out.kinetic_change_J-out.field_work_J;
    auto div=g.gradient_transpose(out.current_load_A_m);
    for(int n=0;n<g.nn;++n)out.continuity_C=std::max(out.continuity_C,std::abs(out.endpoint_charge_C[n]-out.initial_charge_C[n]-dt*div[n]));
    hall::require(out.continuity_C<=1e-12*std::max(maximum(gross),1e-30),"axial path continuity mismatch");
    hall::require(std::abs(out.work_defect_J)<=1e-20+1e-10*(std::abs(out.kinetic_change_J)+std::abs(out.field_work_J)),"axial native work mismatch");
    return out;
}
} // namespace axial
