#pragma once
#include "foil_em.hpp"

namespace foil {
struct Weight {int index,axis;double value;};
struct Path {
    std::vector<Weight>electric,magnetic;
    Path(const Mesh&g,V3 x,V3 D){
        std::vector<double>cuts{0,1};for(int d=0;d<2;++d){double spacing=d?g.dy:g.dx,origin=(d?g.height:g.width)/2,a=(x[d]+origin)/spacing,b=(x[d]+D[d]+origin)/spacing;if(D[d]!=0){hall::require(std::abs(b-a)<10000,"excessive geometric path");for(int k=int(std::floor(std::min(a,b)))+1;k<std::max(a,b);++k){double t=(k-a)/(b-a);if(t>0&&t<1)cuts.push_back(t);}}}
        std::sort(cuts.begin(),cuts.end());
        for(std::size_t k=1;k<cuts.size();++k){double mid=.5*(cuts[k]+cuts[k-1]),span=cuts[k]-cuts[k-1];
            // Two-point Gauss exactly integrates each bilinear path basis.
            for(double q:{-.5773502691896257645,.5773502691896257645}){double t=mid+.5*span*q;double X=(x[0]+D[0]*t+g.width/2)/g.dx,Y=(x[1]+D[1]*t+g.height/2)/g.dy;hall::require(X>=-1e-10&&X<=g.nx+1e-10&&Y>=-1e-10&&Y<=g.ny+1e-10,"path outside field mesh");X=std::clamp(X,0.,double(g.nx));Y=std::clamp(Y,0.,double(g.ny));int i=std::min(int(std::floor(X)),g.nx-1),j=std::min(int(std::floor(Y)),g.ny-1);double a=X-i,b=Y-j,w=.5*span;
                electric.push_back({g.ex(i,j),0,w*(1-b)});electric.push_back({g.ex(i,j+1),0,w*b});electric.push_back({g.ey(i,j),1,w*(1-a)});electric.push_back({g.ey(i+1,j),1,w*a});
                for(int dj:{0,1})for(int di:{0,1})electric.push_back({g.ez(i+di,j+dj),2,w*(di?a:1-a)*(dj?b:1-b)});
                magnetic.push_back({g.bx(i,j),0,w*(1-a)});magnetic.push_back({g.bx(i+1,j),0,w*a});magnetic.push_back({g.by(i,j),1,w*(1-b)});magnetic.push_back({g.by(i,j+1),1,w*b});magnetic.push_back({g.bz(i,j),2,w});
            }
        }
    }
    V3 gather_E(const Vec&E)const{V3 v{};for(const auto&w:electric)v[w.axis]+=w.value*E[w.index];return v;}
    V3 gather_B(const Vec&B)const{V3 v{};for(const auto&w:magnetic)v[w.axis]+=w.value*B[w.index];return v;}
    void deposit(Vec&J,double charge,V3 displacement,double dt)const{for(const auto&w:electric)J[w.index]+=charge*displacement[w.axis]*w.value/dt;}
};
struct Species {std::string id;double charge_C,mass_kg;};
struct Particle {std::uint64_t id;int species;V3 x,u;double count;};
inline double gamma(V3 u){return std::sqrt(1+hall::dot(u,u)/(pic_c*pic_c));}
inline double kinetic(const Particle&p,const std::vector<Species>&registry){double v=hall::dot(p.u,p.u);return p.count*registry.at(p.species).mass_kg*v/(std::sqrt(1+v/(pic_c*pic_c))+1);}
inline double kinetic(const std::vector<Particle>&p,const std::vector<Species>&r){double K=0;for(const auto&v:p)K+=kinetic(v,r);return K;}
inline void point_deposit(const Mesh&g,Vec&q,V3 x,double charge){double X=(x[0]+g.width/2)/g.dx,Y=(x[1]+g.height/2)/g.dy;hall::require(X>=-1e-10&&X<=g.nx+1e-10&&Y>=-1e-10&&Y<=g.ny+1e-10,"point outside mesh");X=std::clamp(X,0.,double(g.nx));Y=std::clamp(Y,0.,double(g.ny));int i=std::min(int(std::floor(X)),g.nx-1),j=std::min(int(std::floor(Y)),g.ny-1);for(int dj:{0,1})for(int di:{0,1})q[g.node(i+di,j+dj)]+=charge*(di?X-i:1-X+i)*(dj?Y-j:1-Y+j);}
inline Vec charge(const Mesh&g,const std::vector<Particle>&p,const std::vector<Species>&r){Vec q(g.nn);for(const auto&v:p)point_deposit(g,q,v.x,v.count*r.at(v.species).charge_C);return q;}
struct Hit {double fraction=1;int surface=-1;}; // 0 lower,1 upper,2 exterior
inline Hit intersection(const Mesh&g,V3 x,V3 D){Hit hit;
    for(int d=0;d<2;++d)if(D[d]!=0){double width=d?g.height:g.width;double f=((D[d]>0?1:-1)*width/2-x[d])/D[d];if(f>=0&&f<=hit.fraction)hit={f,2};}
    if(g.foils)for(int who=0;who<2;++who){double lo[2]={-.0025,who?.0025:-.0027},hi[2]={.0025,who?.0027:-.0025};double enter=-1e100,leave=1e100;bool possible=true;for(int d=0;d<2;++d){if(D[d]==0){if(x[d]<lo[d]||x[d]>hi[d])possible=false;}else{double a=(lo[d]-x[d])/D[d],b=(hi[d]-x[d])/D[d];if(a>b)std::swap(a,b);enter=std::max(enter,a);leave=std::min(leave,b);}}if(possible&&leave>0&&enter<=leave&&enter>=0&&enter<=hit.fraction)hit={enter,who};}
    return hit;
}
inline std::pair<V3,V3>work_kick(V3 u,double qm,V3 E,V3 B,double dt){
    auto next=pic_boris(u,qm,E,B,dt);V3 minus=u,plus=next,D{};for(int d=0;d<3;++d){double a=.5*dt*qm*E[d];minus[d]+=a;plus[d]-=a;}
    double g0=gamma(u),gm=gamma(minus),gp=gamma(plus),g1=gamma(next);for(int d=0;d<3;++d)D[d]=.5*dt*((minus[d]+u[d])/(gm+g0)+(next[d]+plus[d])/(g1+gp));return {next,D};
}
struct BoundaryParticle {std::uint64_t id;int species,surface;V3 point_m;double count,charge_C,kinetic_J,flight_s;};
struct Push {Particle particle;V3 displacement;Hit hit;double flight_s,gyro,displacement_cells;};
inline Push push(const Mesh&g,const Particle&p,const Species&s,const Vec&E,const Vec&B,double dt){
    auto trajectory=[&](double tau){V3 D=p.u;for(double&v:D)v*=tau/gamma(p.u);V3 u{},Bg{};bool converged=false;
        for(int k=0;k<12;++k){auto estimate=intersection(g,p.x,D);V3 gatherD=D;if(estimate.surface>=0)for(double&v:gatherD)v*=estimate.fraction;Path path(g,p.x,gatherD);auto Eg=path.gather_E(E);Bg=path.gather_B(B);auto step=work_kick(p.u,s.charge_C/s.mass_kg,Eg,Bg,tau);double diff=0;for(int d=0;d<3;++d)diff=std::max(diff,std::abs(step.second[d]-D[d]));u=step.first;D=step.second;if(diff<=1e-18){converged=true;break;}}
        hall::require(converged,"finite foil particle path did not converge");auto q=p;q.u=u;for(int d=0;d<3;++d)q.x[d]+=D[d];return Push{q,D,intersection(g,p.x,D),tau,dt*std::abs(s.charge_C)/s.mass_kg*std::sqrt(hall::dot(Bg,Bg)),std::sqrt(hall::dot(D,D))/std::min(g.dx,g.dy)};
    };
    auto result=trajectory(dt);
    if(result.hit.surface>=0){double lo=0,hi=dt;for(int k=0;k<48;++k){double t=.5*(lo+hi);auto trial=trajectory(t);if(trial.hit.surface>=0)hi=t;else lo=t;}result=trajectory(hi);hall::require(result.hit.surface>=0,"lost boundary intersection");
        for(int d=0;d<3;++d){result.displacement[d]*=result.hit.fraction;result.particle.x[d]=p.x[d]+result.displacement[d];}
    }
    hall::require(result.gyro<=.05&&result.displacement_cells<=.2,"finite foil gyro/displacement guard");return result;
}
struct ParticleState {
    State field;std::vector<Particle>particles;std::uint64_t next_id=1;Vec last_current;
    double collected_energy_J=0,escaped_energy_J=0,collected_charge_C=0,escaped_charge_C=0;
    double cumulative_particle_work_J=0,cumulative_circuit_work_J=0,cumulative_background_work_J=0,cumulative_outward_energy_J=0;
};
struct ParticleReceipt {
    FieldReceipt field;std::vector<BoundaryParticle>boundary;
    double particle_work_defect_J=0,gauss_increment_C=0,path_continuity_C=0,plasma_dt=0,nonlinear_relative=0;
    int iterations=0;
    Json json()const{return {{"field",field.json()},{"particle_work_defect_J",particle_work_defect_J},{"gauss_increment_C",gauss_increment_C},{"path_continuity_C",path_continuity_C},{"plasma_dt",plasma_dt},{"nonlinear_relative",nonlinear_relative},{"iterations",iterations},{"boundary_particles",boundary.size()}};}
};
class System {
    const Mesh&g;Fields field;double dt;
public:
    std::vector<Species>registry;
    System(const Mesh&mesh,double step,double voltage,double dc_current,std::vector<Species>species):g(mesh),field(mesh,step,voltage,dc_current),dt(step),registry(std::move(species)){
        std::set<std::string>ids;for(const auto&s:registry)hall::require(!s.id.empty()&&ids.insert(s.id).second&&s.charge_C!=0&&std::isfinite(s.charge_C)&&std::isfinite(s.mass_kg)&&s.mass_kg>0,"invalid charged species registry");
    }
    ParticleState initial(double voltage)const{return {field.initial(voltage),{},1,Vec(g.ne)};}
    const Fields&fields()const{return field;}
    void add(ParticleState&state,Particle p)const{hall::require(p.species>=0&&p.species<int(registry.size())&&p.count>0&&std::isfinite(p.count),"invalid particle species/count");for(const auto&v:{p.x,p.u})for(double x:v)hall::require(std::isfinite(x),"nonfinite particle");hall::require(std::abs(p.x[0])<g.width/2&&std::abs(p.x[1])<g.height/2,"particle outside domain");if(g.foils)hall::require(!(std::abs(p.x[0])<=.0025&&std::abs(p.x[1])>=.0025&&std::abs(p.x[1])<=.0027),"particle inside electrode");p.id=state.next_id++;state.particles.push_back(p);}
    ParticleReceipt advance(ParticleState&accepted,double voltage)const{
        const auto old=accepted;auto guess=old.field;ParticleState next=old;ParticleReceipt receipt;
        Vec plasma(g.nn);for(const auto&p:old.particles){const auto&s=registry.at(p.species);point_deposit(g,plasma,p.x,p.count*s.charge_C*s.charge_C/(eps*s.mass_kg));}
        for(int j=0;j<=g.ny;++j)for(int i=0;i<=g.nx;++i)receipt.plasma_dt=std::max(receipt.plasma_dt,dt*std::sqrt(plasma[g.node(i,j)]/(g.dual_x(i)*g.dual_y(j)*g.depth)));
        hall::require(receipt.plasma_dt<=.1,"sum-species full-interval plasma guard");const auto q0=charge(g,old.particles,registry);Vec gross_charge(g.nn);for(const auto&p:old.particles)point_deposit(g,gross_charge,p.x,p.count*std::abs(registry.at(p.species).charge_C));bool converged=false;double wall_energy=0,escape_energy=0;
        for(int iteration=0;iteration<16;++iteration){Vec E(g.ne),B(g.nb);for(int e=0;e<g.ne;++e)E[e]=.5*(old.field.E[e]+guess.E[e]);for(int b=0;b<g.nb;++b)B[b]=.5*(old.field.B[b]+guess.B[b]);
            Vec J(g.ne),end_charge(g.nn);next.particles.clear();receipt.boundary.clear();wall_energy=escape_energy=0;
            for(const auto&p:old.particles){const auto&s=registry.at(p.species);auto pushed=push(g,p,s,E,B,dt);Path path(g,p.x,pushed.displacement);path.deposit(J,p.count*s.charge_C,pushed.displacement,dt);point_deposit(g,end_charge,pushed.particle.x,p.count*s.charge_C);
                if(pushed.hit.surface<0)next.particles.push_back(pushed.particle);else{double K=kinetic(pushed.particle,registry);receipt.boundary.push_back({p.id,p.species,pushed.hit.surface,pushed.particle.x,p.count,p.count*s.charge_C,K,pushed.flight_s});if(pushed.hit.surface<2)wall_energy+=K;else escape_energy+=K;}
            }
            auto div=g.gradient_transpose(J);receipt.path_continuity_C=0;for(int n=0;n<g.nn;++n)receipt.path_continuity_C=std::max(receipt.path_continuity_C,std::abs(end_charge[n]-q0[n]-dt*div[n]));
            auto solved=field.advance(old.field,J,voltage);next.field=std::move(solved.first);receipt.field=solved.second;next.last_current=std::move(J);
            double difference=0,scale=std::max(1.,maximum(next.field.E));for(int e=0;e<g.ne;++e)difference=std::max(difference,std::abs(next.field.E[e]-guess.E[e]));for(int b=0;b<g.nb;++b){difference=std::max(difference,pic_c*std::abs(next.field.B[b]-guess.B[b]));scale=std::max(scale,pic_c*std::abs(next.field.B[b]));}
            receipt.nonlinear_relative=difference/scale;receipt.iterations=iteration+1;guess=next.field;if(receipt.nonlinear_relative<=1e-13){converged=true;break;}
        }
        hall::require(converged,"finite foil field/particle iteration failed");
        double dK=kinetic(next.particles,registry)-kinetic(old.particles,registry)+wall_energy+escape_energy;
        receipt.particle_work_defect_J=dK-receipt.field.particle_work_J;
        const double particle_scale=std::abs(dK)+std::abs(receipt.field.particle_work_J)+wall_energy+escape_energy;
        hall::require(std::abs(receipt.particle_work_defect_J)<=1e-20+1e-10*particle_scale,"particle/field work mismatch");
        const auto Q0=g.field_charge(old.field.E),Q1=g.field_charge(next.field.E),q1=charge(g,next.particles,registry);
        double charge_scale=0;for(int n=0;n<g.nn;++n){charge_scale=std::max(charge_scale,std::abs(q1[n]-q0[n]));if(g.owner[n]<0&&!g.outer_node[n])receipt.gauss_increment_C=std::max(receipt.gauss_increment_C,std::abs((Q1[n]-Q0[n])-(q1[n]-q0[n])));}
        hall::require(receipt.path_continuity_C<=1e-12*std::max(maximum(gross_charge),1e-30),"finite foil path continuity mismatch");
        hall::require(receipt.gauss_increment_C<=1e-22+2e-10*charge_scale,"finite foil Gauss propagation mismatch");
        const auto&r=receipt.field;double field_scale=std::abs(r.field_change_J)+std::abs(r.particle_work_J)+std::abs(r.outward_boundary_energy_J)+std::abs(r.background_supply_work_J)+std::abs(r.circuit_supply_work_J);
        hall::require(std::abs(r.balance_J)<=1e-20+1e-10*field_scale,"finite foil field/circuit/boundary energy mismatch");
        next.collected_energy_J+=wall_energy;next.escaped_energy_J+=escape_energy;for(const auto&b:receipt.boundary){if(b.surface<2)next.collected_charge_C+=b.charge_C;else next.escaped_charge_C+=b.charge_C;}
        next.cumulative_particle_work_J+=r.particle_work_J;next.cumulative_circuit_work_J+=r.circuit_supply_work_J;next.cumulative_background_work_J+=r.background_supply_work_J;next.cumulative_outward_energy_J+=r.outward_boundary_energy_J;accepted=std::move(next);return receipt;
    }
};
}
