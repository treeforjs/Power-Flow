#include "tangential_joint.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace tangential_onepass {
namespace se=surface_emission;
namespace {
constexpr double pi=3.14159265358979323846,hbar=se::planck/(2*pi),c=se::qe/(16*pi*se::eps0),kb_eV=se::kb/se::qe;
const double kappa=2*std::sqrt(2*se::me*se::qe)/hbar;
void need(bool b,const char* s){if(!b)throw std::invalid_argument(s);}
double dot(const se::Vector& a,const se::Vector& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
se::Vector cross(const se::Vector& a,const se::Vector& b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
se::Vector mul(se::Vector a,double s){for(double& x:a)x*=s;return a;}
double norm(const se::Vector& a){return std::sqrt(dot(a,a));}
double softplus(double x){if(x>40)return x;if(x< -40)return std::exp(x);return std::log1p(std::exp(x));}
double logD(double G){return -softplus(G);}
void validate_policy(const se::LayerPolicy& p) {
    need(std::isfinite(p.local_linear_vacuum_extent_m) && p.local_linear_vacuum_extent_m>0 &&
        p.image_energy_fraction>0 && p.image_energy_fraction<=.2 && p.image_force_fraction>0 && p.image_force_fraction<=.2 &&
        p.maximum_quantum_magnetic_energy_fraction==.2 && !p.PIC_already_applies_microscopic_single_electron_image_force,
        "Tangential prototype requires nominal0.2 policy and explicit separate image ownership");
}
struct Rule {std::vector<double> x,w;};
Rule legendre(int n) {
    Rule q;q.x.resize(n);q.w.resize(n);
    for(int i=0;i<(n+1)/2;++i) {
        double z=std::cos(pi*(i+.75)/(n+.5)),derivative=0;
        for(int it=0;it<30;++it){double a=1,b=0;for(int j=1;j<=n;++j){const double old=b;b=a;a=((2*j-1)*z*b-(j-1)*old)/j;}derivative=n*(z*a-b)/(z*z-1);const double next=z-a/derivative;if(std::abs(next-z)<2e-15){z=next;break;}z=next;}
        q.x[i]=-z;q.x[n-1-i]=z;const double w=2/((1-z*z)*derivative*derivative);q.w[i]=q.w[n-1-i]=w;
    }return q;
}
Rule laguerre(int n) {
    Rule q;double z=0;
    for(int i=1;i<=n;++i) {
        if(i==1)z=3./(1+2.4*n);else if(i==2)z+=15./(1+2.5*n);else{const double a=i-2;z+=(1+2.55*a)/(1.9*a)*(z-q.x[i-3]);}
        double derivative=0,previous=0;
        for(int it=0;it<40;++it){double a=1,b=0;for(int j=1;j<=n;++j){const double old=b;b=a;a=((2*j-1-z)*b-(j-1)*old)/j;}derivative=n*(a-b)/z;previous=b;const double next=z-a/derivative;if(std::abs(next-z)<2e-14*std::max(1.,std::abs(z))){z=next;break;}z=next;}
        const double w=-1/(derivative*n*previous);need(z>0 && w>0 && std::isfinite(w),"Invalid Laguerre rule");q.x.push_back(z);q.w.push_back(w);
    }
    long double sum=0,mean=0;for(int i=0;i<n;++i){sum+=q.w[i];mean+=q.w[i]*q.x[i];}
    need(std::abs(sum-1)<1e-11 && std::abs(mean-1)<1e-11,"Laguerre supply normalization failed");return q;
}
template<class Function> double root(Function f,double a,double b) {
    const bool sign=f(a)>0;need(sign!=(f(b)>0),"Unbracketed canonical root");
    for(int i=0;i<72;++i){const double m=(a+b)/2;if((f(m)>0)==sign)a=m;else b=m;}return (a+b)/2;
}
struct Barrier {
    double phi,F,B,v,alpha,Feff,xs=0,xmin=0,top=0,bottom=0,width=0;bool has_pocket=false;
    double U(double x)const{return phi-Feff*x-c/x+alpha*x*x;}
    Barrier(double work,double electric,double magnetic,double canonical_v):phi(work),F(electric),B(magnetic),v(canonical_v),alpha(se::qe*B*B/(2*se::me)),Feff(F+B*v) {
        if(alpha==0){if(Feff<=0)return;xs=std::sqrt(c/Feff);xmin=std::numeric_limits<double>::infinity();top=U(xs);bottom=-std::numeric_limits<double>::infinity();has_pocket=true;}
        else {
            const double critical=std::cbrt(c/alpha);if(Feff<=3*std::cbrt(alpha*alpha*c))return;
            auto derivative=[&](double x){return -Feff+c/(x*x)+2*alpha*x;};
            xs=root(derivative,critical*1e-12,critical);xmin=root(derivative,critical,std::max(2*critical,Feff/(2*alpha)));top=U(xs);bottom=U(xmin);has_pocket=true;
        }
        const double curvature=-2*c/(xs*xs*xs)+2*alpha;need(curvature<0,"Canonical local saddle curvature");
        width=hbar*std::sqrt(-se::qe*curvature/se::me)/(2*pi*se::qe);
    }
    double outer(double E)const {
        need(has_pocket && E>bottom,"No resolved canonical overlap channel");if(E>=top)return xs;
        double limit=xmin;if(!std::isfinite(limit)){limit=2*xs;while(U(limit)>E)limit*=2;}
        return root([&](double x){return U(x)-E;},xs,limit);
    }
    double gamow(double E,const Rule& rule)const {
        if(E>=top)return (top-E)/width;
        const double a=root([&](double x){return U(x)-E;},xs*1e-12,xs),b=outer(E),mid=(a+b)/2,half=(b-a)/2;
        long double sum=0;
        for(std::size_t i=0;i<rule.x.size();++i){const double angle=pi*(rule.x[i]+1)/2,x=mid-half*std::cos(angle);sum+=rule.w[i]*std::sqrt(std::max(0.,U(x)-E))*half*std::sin(angle);}
        return double(sum)*pi/2*kappa;
    }
};
struct Bounds {double current=0,heat=0;};
Bounds unknown_bounds(double phi,double F,double B,double t,double emin,double emax,double r,double overlap) {
    const double C=se::ideal_richardson_A_m2_K2()/(kb_eV*kb_eV),alpha=se::qe*B*B/(2*se::me),k=B*std::sqrt(2*se::qe/se::me);
    Bounds out;
    if(B>0) {
        const double threshold=std::max(3*std::cbrt(alpha*alpha*c),2*std::sqrt(alpha*(phi-emin+overlap)));
        const double etcrit=F>threshold?std::pow((F-threshold)/k,2):0,a=emin+etcrit;
        if(a>0){const double exponential=std::exp(-a/t);out.current+=C*t*t*exponential;out.heat+=C*t*t*(a+2*t)*exponential;}
        else {out.current+=C*(a*a/2+t*t-a*t);out.heat+=C*(-a*a*a/6+2*t*t*t-a*t*t);}
    }
    // Deep-normal bound on a fixed SN barrier slab for all Et<=r and angles.
    const double Fmax=F+k*std::sqrt(r),top=phi-2*std::sqrt(c*Fmax),m=(top+emin)/2;
    need(top>emin,"No deep-normal barrier slab for the declared prototype domain");
    const double w=phi-m,discriminant=w*w-4*c*Fmax,L=std::sqrt(discriminant)/Fmax,a0=kappa*L,z=std::sqrt(m-emin),factor=2*std::exp(-a0*z);
    const double I0=factor*(z/a0+1/(a0*a0));
    const double I2=factor*(z*z*z/a0+3*z*z/(a0*a0)+6*z/(a0*a0*a0)+6/std::pow(a0,4));
    out.current+=C*r*I0;out.heat+=C*r*((std::abs(m)+r)*I0+I2);
    // Et tail: split at E=-Et/2. FD bounds the upper half; a triangular
    // SN subinterval bounds the lower half of the one-pass JWKB barrier.
    need(r>=2*phi/3 && std::pow(phi+r/2,2)>=16*c*Fmax,"Transverse-tail action bound not established");
    const double a=kappa/(8*std::sqrt(2.)*(F/std::sqrt(r)+k));
    const double expA=std::exp(-r/(2*t)),expB=std::exp(-a*(phi+r/2));
    out.current+=2*C*t*t*expA+2*C/(a*a)*expB;
    out.heat+=C*t*t*(r+4*t)*expA+C*((4*phi+3*r)/(a*a)+8/(a*a*a))*expB;
    const double high=std::exp(-emax/t);out.current+=C*t*t*high;out.heat+=C*t*t*(emax+2*t)*high;
    return out;
}
struct Weighted {Node node;double log_weight;};
}
Distribution evaluate(const se::Stage& stage,const Config& cfg,const se::LayerPolicy& policy) {
    validate_policy(policy);
    need(cfg.continuum_free_electron_reservoir_selected && cfg.incoherent_onepass_matching_selected,"Conditional continuum-reservoir/one-pass assumptions must be selected");
    need(stage.T_K>=298.2 && stage.T_K<=1500 && cfg.work_function_eV>=2 && cfg.work_function_eV<=6 &&
        cfg.transverse_nodes>=8 && cfg.transverse_nodes<=64 && cfg.azimuth_nodes>=8 && cfg.azimuth_nodes<=64 &&
        cfg.normal_panels>=1 && cfg.normal_panels<=32 && cfg.action_nodes>=16 && cfg.action_nodes<=96 &&
        cfg.normal_lower_eV<0 && cfg.normal_lower_eV>=-4 && cfg.normal_tail_kT>=40 && cfg.transverse_tail_kT>=80,
        "Tangential prototype numerical/thermal domain");
    const auto nominal=se::prepare(stage,se::Config{cfg.work_function_eV});const double F=nominal.spectrum.extraction_field_V_m,B=norm(stage.total_B_T),Bn=dot(stage.total_B_T,stage.outward_normal);
    need(F>=5e7 && F<=1e8 && std::abs(Bn)<=1e-12*std::max(1.,B) && B<=40,"Tangent-only prototype does not cover oblique B or this field range");
    se::Vector Etangent=stage.total_E;for(int j=0;j<3;++j)Etangent[j]+=F*stage.outward_normal[j];
    need(norm(Etangent)<=1e-12*F,"Tangent-only prototype requires normal E");
    const double t=kb_eV*stage.T_K,scale=std::max(t,nominal.spectrum.fermi_decay_width_eV),delta=policy.image_energy_fraction*scale;
    const double emin=cfg.normal_lower_eV,emax=cfg.work_function_eV+cfg.normal_tail_kT*t,r=2*cfg.work_function_eV+cfg.transverse_tail_kT*t;
    const double vmax=std::sqrt(2*se::qe*r/se::me),alpha=se::qe*B*B/(2*se::me);
    const double quantum_extent=B*vmax<=F/2?2*(cfg.work_function_eV-emin+delta)/(F-B*vmax):(F+B*vmax)/(2*alpha);
    need(policy.local_linear_vacuum_extent_m>=quantum_extent,"Entire resolved quantum region is not certified inside actual local field/geometry extent");
    const auto transverse=laguerre(cfg.transverse_nodes),normal=legendre(8),action_rule=legendre(cfg.action_nodes);
    const double C=se::ideal_richardson_A_m2_K2()/(kb_eV*kb_eV);Distribution out;out.stage=stage;out.config=cfg;out.magnetic_reservoir_energy_scale_over_kT=hbar*se::qe*B/(se::me*se::qe*t);out.quantum_region_extent_bound_m=quantum_extent;
    std::vector<Weighted> values;double peak=-std::numeric_limits<double>::infinity();
    for(int j=0;j<cfg.transverse_nodes;++j) {
        const double transverse_energy=t*transverse.x[j];if(transverse_energy>r)continue;
        const double speed=std::sqrt(2*se::qe*transverse_energy/se::me);
        for(int k=0;k<cfg.azimuth_nodes;++k) {
            const double theta=2*pi*(k+.5)/cfg.azimuth_nodes;const Barrier barrier(cfg.work_function_eV,F,B,speed*std::cos(theta));
            if(!barrier.has_pocket){++out.unresolved_quadrature_phase_points;continue;}
            std::vector<double> edges{emin,emax};
            for(double point:{-transverse_energy,barrier.bottom+delta,barrier.top,barrier.top-8*barrier.width,barrier.top+8*barrier.width,barrier.top-8*t,barrier.top+8*t})if(point>emin && point<emax && std::isfinite(point))edges.push_back(point);
            std::sort(edges.begin(),edges.end());edges.erase(std::unique(edges.begin(),edges.end()),edges.end());
            for(std::size_t segment=1;segment<edges.size();++segment)for(int panel=0;panel<cfg.normal_panels;++panel) {
                const double a=edges[segment-1]+(edges[segment]-edges[segment-1])*panel/cfg.normal_panels,b=edges[segment-1]+(edges[segment]-edges[segment-1])*(panel+1)/cfg.normal_panels;
                for(std::size_t n=0;n<normal.x.size();++n) {
                    const double E=(a+b)/2+(b-a)*normal.x[n]/2;
                    if(E-delta<=barrier.bottom){++out.unresolved_quadrature_phase_points;continue;}
                    const double transmission=logD(barrier.gamow(E,action_rule));
                    const double log_weight=std::log(C*t*transverse.w[j]*normal.w[n]*(b-a)/(2*cfg.azimuth_nodes))+transverse.x[j]-softplus((E+transverse_energy)/t)+transmission;
                    peak=std::max(peak,log_weight);values.push_back({{E,transverse_energy,theta,0,0,transmission},log_weight});
                }
            }
        }
    }
    need(std::isfinite(peak) && !values.empty(),"No resolved known-component source spectrum");
    long double sum=0,moment=0;for(const auto& v:values){const double w=std::exp(v.log_weight-peak);sum+=w;moment+=w*(v.node.normal_eV+v.node.transverse_eV);}
    out.known_component_log_current_A_m2=peak+std::log(double(sum));out.known_component_current_A_m2=std::exp(out.known_component_log_current_A_m2);out.known_component_mean_metal_energy_eV=double(moment/sum);
    out.known_component_complete_heat_W_m2=-out.known_component_current_A_m2*out.known_component_mean_metal_energy_eV;
    long double cdf=0;out.known_component_nodes.reserve(values.size());
    for(auto& v:values){v.node.probability=std::exp(v.log_weight-peak)/double(sum);cdf+=v.node.probability;v.node.cdf=double(cdf);out.known_component_nodes.push_back(v.node);}out.known_component_nodes.back().cdf=1;
    const auto bound=unknown_bounds(cfg.work_function_eV,F,B,t,emin,emax,r,delta);out.unresolved_current_upper_A_m2=bound.current;out.unresolved_absolute_heat_upper_W_m2=bound.heat;
    out.truncation_bounds_within_declared_budget=bound.current<=cfg.unresolved_relative_budget*out.known_component_current_A_m2 && bound.heat<=cfg.unresolved_relative_budget*out.known_component_current_A_m2*std::max(1.,std::abs(out.known_component_mean_metal_energy_eV));
    return out;
}
se::LayerTrial sample_known_component(const Distribution& d,const se::LayerPolicy& policy,double weight,double u,std::uint64_t id) {
    validate_policy(policy);
    need(d.truncation_bounds_within_declared_budget && u>0 && u<1 && weight>0 && std::isfinite(weight) && id>0,"Known component sampling requires explicit small unresolved bound and caller draws/weight");
    const auto it=std::lower_bound(d.known_component_nodes.begin(),d.known_component_nodes.end(),u,[](const Node& n,double q){return n.cdf<q;});need(it!=d.known_component_nodes.end(),"Invalid joint CDF");
    const auto& n=d.stage.outward_normal;const double F=se::extraction_field(d.stage.total_E,n),B=norm(d.stage.total_B_T),scale=kb_eV*d.stage.T_K;
    se::Vector b;if(B>0)b=mul(d.stage.total_B_T,1/B);else {const se::Vector axis=std::abs(n[0])<.8?se::Vector{1,0,0}:se::Vector{0,1,0};b=cross(n,axis);b=mul(b,1/norm(b));}
    const auto direction=cross(n,b);const double ut=std::sqrt(2*se::qe*it->transverse_eV/se::me),vz=ut*std::cos(it->azimuth_rad);
    const Barrier barrier(d.config.work_function_eV,F,B,vz);const double x=barrier.outer(it->normal_eV-policy.image_energy_fraction*scale);
    need(x>0 && x<policy.local_linear_vacuum_extent_m && !policy.PIC_already_applies_microscopic_single_electron_image_force,"Joint overlap leaves local geometry or duplicates microscopic image force");
    se::Vector proper{};for(int j=0;j<3;++j)proper[j]=direction[j]*(vz-se::qe*B*x/se::me)+b[j]*ut*std::sin(it->azimuth_rad);
    const double image=-c/x,K_eV=it->normal_eV+it->transverse_eV-d.config.work_function_eV+F*x-image,K=K_eV*se::qe;
    const double u2=K*(K+2*se::me*se::c_light*se::c_light)/(se::me*se::me*se::c_light*se::c_light),normal2=u2-dot(proper,proper);need(K>0 && normal2>=0,"Joint canonical overlap normal momentum unavailable");
    for(int j=0;j<3;++j)proper[j]+=n[j]*std::sqrt(normal2);
    se::LayerTrial trial;auto& s=trial.state;const double q=se::qe*weight;
    se::Ledger ledger{weight,-q,q,-(it->normal_eV+it->transverse_eV)*q,K*weight,0,d.config.work_function_eV*q,image*q,F*x*q};
    s.particle={id,d.stage.accepted_state,d.stage.patch_id,weight,it->normal_eV,it->transverse_eV,proper,ledger};s.position_from_surface_m=mul(n,x);s.outward_normal=n;s.current_time_s=d.stage.stage_time;
    s.return_overlap_distance_m=x;s.handoff_distance_m=std::max(std::sqrt(c/(policy.image_force_fraction*F)),c/(policy.image_energy_fraction*scale));s.work_function_eV=d.config.work_function_eV;s.source_energy_scale_eV=scale;
    s.quantum_magnetic_energy_fraction=(B*x*ut+se::qe*B*B*x*x/(2*se::me))/scale;s.mesh_revision=d.stage.mesh_revision;s.material_revision=d.stage.material_revision;
    trial.exchange={weight,-q,q,ledger.wall_heat_J,K*weight,ledger.electron_reference_change_J,ledger.image_reference_change_J,ledger.resolved_matching_field_work_J};
    trial.paths.push_back({{},s.position_from_surface_m,s.current_time_s,s.current_time_s,-q,F*x*q});return trial;
}
QuantumReturn resolve_quantum_return(const se::Stage& stage,const Config& cfg,const se::LayerPolicy& policy,const se::LayerTrial& incoming,double uniform) {
    validate_policy(policy);
    need(cfg.continuum_free_electron_reservoir_selected && cfg.incoherent_onepass_matching_selected && uniform>0 && uniform<1 &&
        incoming.state.status==se::LayerStatus::Recollected && !incoming.paths.empty(),"Quantum return requires a provisional collected trial and distinct caller draw");
    const auto& state=incoming.state;const auto& path=incoming.paths.back();const auto& n=state.outward_normal;
    need(stage.accepted_state>0 && stage.accepted_state==stage.field_accepted_state && stage.accepted_state==stage.thermal_accepted_state &&
        stage.field_time==stage.stage_time && stage.thermal_time==stage.stage_time && stage.mesh_revision==state.mesh_revision &&
        stage.material_revision==state.material_revision && stage.patch_id==state.particle.patch_id && stage.outward_normal==n &&
        cfg.work_function_eV==state.work_function_eV,"Quantum return stage/material/source identity mismatch");
    const double F=se::extraction_field(stage.total_E,n),B=norm(stage.total_B_T),Bn=dot(stage.total_B_T,n),x=dot(path.from_m,n);
    se::Vector tangentE=stage.total_E;for(int j=0;j<3;++j)tangentE[j]+=F*n[j];
    need(F>=5e7 && F<=1e8 && B<=40 && std::abs(Bn)<=1e-12*std::max(1.,B) && norm(tangentE)<=1e-12*F &&
        x>0 && norm(path.from_m)<=policy.local_linear_vacuum_extent_m && path.start_time_s==path.end_time_s &&
        path.end_time_s==state.current_time_s && state.current_time_s>=stage.t && state.current_time_s<=stage.t+stage.dt+1e-24 &&
        norm(se::Vector{path.to_m[0]-state.position_from_surface_m[0],path.to_m[1]-state.position_from_surface_m[1],path.to_m[2]-state.position_from_surface_m[2]})<1e-20,
        "Quantum return has no resolved inward bridge in the supported geometry/epoch");
    const auto& u=state.particle.proper_velocity_m_s;const double un=dot(u,n);need(un<0,"Quantum return momentum is not inward");
    se::Vector b;if(B>0)b=mul(stage.total_B_T,1/B);else {const se::Vector axis=std::abs(n[0])<.8?se::Vector{1,0,0}:se::Vector{0,1,0};b=cross(n,axis);b=mul(b,1/norm(b));}
    const auto direction=cross(n,b);auto metal=u;
    for(int j=0;j<3;++j)metal[j]+=-un*n[j]+se::qe*B*x/se::me*direction[j];
    const double Et=se::me*dot(metal,metal)/(2*se::qe),K=se::kinetic_energy_J(u),image=-se::qe*c/x;
    const double total_metal=K/se::qe+cfg.work_function_eV+image/se::qe-F*x,E=total_metal-Et;
    const Barrier barrier(cfg.work_function_eV,F,B,dot(metal,direction));
    need(barrier.has_pocket && E>barrier.bottom,"Returning channel requires an unresolved quantum treatment; rollback instead of absorption0");
    const double logarithm=logD(barrier.gamow(E,legendre(cfg.action_nodes)));
    QuantumReturn result;result.trial=incoming;result.log_transmission=logarithm;result.transmission=std::exp(logarithm);result.reconstructed_normal_eV=E;result.reconstructed_transverse_eV=Et;result.absorbed=std::log(uniform)<logarithm;
    if(result.absorbed)return result;
    // Subtract the provisional reverse-quantum bridge exactly, retaining all
    // real image/Lorentz steps and their q dx and work. No accepted state exists
    // between this undo and the returned reflected trial.
    auto& e=result.trial.exchange;const double weight=state.particle.weight,charge=-se::qe*weight,phi=se::qe*cfg.work_function_eV;
    e.physical_electron_number_change+=weight;e.particle_charge_change_C+=charge;e.electrode_charge_change_C-=charge;
    e.wall_heat_J-=(K+phi+image)*weight+path.field_work_J;e.particle_K_change_J+=K*weight;
    e.electron_reference_change_J+=phi*weight;e.image_energy_change_J+=image*weight;e.field_work_J-=path.field_work_J;
    auto& s=result.trial.state;s.position_from_surface_m=path.from_m;s.status=se::LayerStatus::Active;
    for(int j=0;j<3;++j)s.particle.proper_velocity_m_s[j]-=2*un*n[j];
    result.trial.paths.pop_back();return result;
}
}
