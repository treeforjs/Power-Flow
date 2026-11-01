#include "surface_emission.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace surface_emission {
namespace {
constexpr double pi=3.1415926535897932384626433832795;
constexpr double hbar=planck/(2*pi),kb_eV=kb/qe,image_c=qe/(16*pi*eps0);
void require(bool okay,const char* message){if(!okay)throw std::invalid_argument(message);}
double dot(const Vector& a,const Vector& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
Vector cross(const Vector& a,const Vector& b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
bool finite(const Vector& v){return std::all_of(v.begin(),v.end(),[](double x){return std::isfinite(x);});}
void unit(const Vector& n){require(finite(n) && std::abs(dot(n,n)-1)<1e-12,"Material-outward normal must be finite and unit length");}
void config(const Config& c) {
    require(std::isfinite(c.work_function_eV) && c.work_function_eV>=2 && c.work_function_eV<=6,"Explicit work function outside numerical domain2..6eV");
    require(c.energy_panels>=16 && c.energy_panels<=1024 && c.action_nodes>=16 && c.action_nodes<=256 && c.action_nodes%2==0,
            "Bounded quadrature order required");
    require(std::isfinite(c.tail_widths) && c.tail_widths>=32 && c.tail_widths<=80,"Tail width outside32..80");
}
std::pair<std::vector<double>,std::vector<double>> legendre(int n) {
    std::vector<double> x(n),w(n);
    for(int i=0;i<(n+1)/2;++i) {
        double z=std::cos(pi*(i+.75)/(n+.5)),derivative=0;
        for(int j=0;j<100;++j) {
            double p=1,q=0;
            for(int k=1;k<=n;++k){const double next=((2*k-1)*z*p-(k-1)*q)/k;q=p;p=next;}
            derivative=n*(z*p-q)/(z*z-1);const double next=z-p/derivative;
            if(std::abs(next-z)<2e-15){z=next;break;}z=next;
        }
        // Re-evaluate derivative at converged root.
        double p=1,q=0;for(int k=1;k<=n;++k){const double next=((2*k-1)*z*p-(k-1)*q)/k;q=p;p=next;}
        derivative=n*(z*p-q)/(z*z-1);
        x[i]=-z;x[n-1-i]=z;w[i]=w[n-1-i]=2/((1-z*z)*derivative*derivative);
    }
    return {x,w};
}
struct Barrier {
    double F,phi,lowering,top,top_slope;
    std::vector<double> cosine,sine2,weights;
    Barrier(double field,const Config& c):F(field),phi(c.work_function_eV),lowering(2*std::sqrt(image_c*field)),top(phi-lowering) {
        const double z=std::sqrt(image_c/F),omega=std::sqrt(2*qe*image_c/(me*z*z*z));
        top_slope=2*pi*qe/(hbar*omega);
        const auto [x,w]=legendre(c.action_nodes);
        for(std::size_t i=0;i<x.size();++i){const double angle=pi*(x[i]+1)/2;cosine.push_back(std::cos(angle));sine2.push_back(std::pow(std::sin(angle),2));weights.push_back(w[i]*pi/2);}
    }
    double gamow(double E) const {
        if(E>=top)return -top_slope*(E-top); // local inverted-parabola continuation above saddle
        const double A=phi-E,disc=std::sqrt((A-lowering)*(A+lowering));
        const double z1=2*image_c/(A+disc),z2=(A+disc)/(2*F),mid=(z1+z2)/2,half=(z2-z1)/2;
        long double integral=0;
        for(std::size_t i=0;i<weights.size();++i)integral+=weights[i]*sine2[i]/std::sqrt(mid+half*cosine[i]);
        const double G=2*std::sqrt(2*me*qe)/hbar*std::sqrt(F)*half*half*double(integral);
        require(std::isfinite(G) && G>=0,"Nonfinite Schottky-Nordheim action");return G;
    }
};
double softplus(double x){if(x>40)return x+std::log1p(std::exp(-x));return std::log1p(std::exp(x));}
double log_transmission(double G){return -softplus(G);}
double log_supply(double E,double kT) {
    if(kT==0)return E<0?std::log(-E):-std::numeric_limits<double>::infinity();
    const double eta=-E/kT;
    return std::log(kT)+(eta< -36?eta:std::log(softplus(eta)));
}
double fermi1_negative(double eta) {
    const double z=1/(1+std::exp(-eta)),L=softplus(eta);double series=0,power=z;
    for(int k=1;k<=100;++k){const double term=power/(double(k)*k);series+=term;if(term<1e-17*std::max(series,1e-300))break;power*=z;}
    return series+.5*L*L;
}
double mean_transverse(double E,double kT) {
    if(kT==0)return -.5*E;
    const double eta=-E/kT;
    if(eta< -36)return kT;
    const double L=softplus(eta),f1=eta<=0?fermi1_negative(eta):.5*eta*eta+pi*pi/6-fermi1_negative(-eta);
    return kT*f1/L;
}
void uniform(double u){require(std::isfinite(u) && u>0 && u<1,"Caller-owned draws must lie strictly inside(0,1)");}
struct PlaneState {double U,field_work,image_energy;};
PlaneState plane(const Distribution& d,const MatchingPlane& p) {
    require(p.microscopic_matching_admitted && p.field_path_work_owned_once && p.image_reference_port_admitted,
            "Unresolved microscopic matching/field/image energy ownership: launch refused");
    if(d.extraction_field_V_m<=0) {
        require(p.kind==PlaneKind::ThermalVacuumSurfaceStep && p.outward_distance_m==0,"Retarding/zero-field branch requires explicit local vacuum surface step");
        return {d.config.work_function_eV,0,0};
    }
    require(p.kind==PlaneKind::ExplicitSchottkyNordheimHandoff && std::isfinite(p.outward_distance_m) && p.outward_distance_m>0,
            "Extracting field requires a declared SN matching surface; no automatic launch offset");
    require(p.outward_distance_m>=required_SN_matching_distance_m(d),"Matching plane is inside the emission barrier or cuts the occupied emitted spectrum");
    const double field=d.extraction_field_V_m*p.outward_distance_m,image=-image_c/p.outward_distance_m;
    const double U=d.config.work_function_eV-field+image;
    require(std::isfinite(U) && std::isfinite(field) && std::isfinite(image),"Nonfinite matching-plane energy");
    require(d.nodes.empty() || d.nodes.front().normal_energy_eV>=U,"Matching surface cuts occupied emitted spectrum; no rejection/resampling/current renormalization allowed");
    return {U,field,image};
}
}
double ideal_richardson_A_m2_K2(){return 4*pi*me*qe*kb*kb/(planck*planck*planck);}
double extraction_field(const Vector& E,const Vector& n){unit(n);require(finite(E),"Solved total one-sided vacuum E must be finite");const double F=-dot(E,n);require(std::isfinite(F),"Normal field overflow");return F;}
double schottky_lowering_eV(double F){require(std::isfinite(F) && F>=0,"Positive extraction field required");return 2*std::sqrt(image_c*F);}
double action(double E,double F,const Config& c){config(c);require(std::isfinite(E) && F>=1 && F<=5e9,"Action domain requires1..5e9V/m and finite energy");return Barrier(F,c).gamow(E);}
double transmission(double E,double F,const Config& c){return std::exp(log_transmission(action(E,F,c)));}
Distribution evaluate(double T,double F,const Config& c) {
    config(c);require(std::isfinite(T) && (T==0 || (T>=1 && T<=3000)),"Numerical free-electron law domain is T=0 or1..3000K; no extrapolation/clamp");
    require(std::isfinite(F) && std::abs(F)<=5e9 && (F<=0 || F>=1),"Signed normal-field numerical domain exceeded");
    const double lowering=F>0?schottky_lowering_eV(F):0;
    require(lowering<c.work_function_eV*.98,"Barrier collapse outside admitted numerical domain");
    Distribution d{};d.config=c;d.T_K=T;d.extraction_field_V_m=F;d.barrier_lowering_eV=lowering;
    if(T==0 && F<=0)return d;
    const double kT=kb_eV*T;
    const double decay=F>0?F/(1.5*6.830890e9*std::sqrt(c.work_function_eV)):0;
    const double low=F>0?-std::max({2.,c.tail_widths*decay,c.tail_widths*kT}):c.work_function_eV;
    const double high=T>0?c.work_function_eV+c.tail_widths*kT:0;
    const auto barrier=F>0?std::make_unique<Barrier>(F,c):nullptr;
    if(barrier)d.fermi_decay_width_eV=-2e-5/(barrier->gamow(1e-5)-barrier->gamow(-1e-5));
    std::vector<double> breaks{low,high};
    if(F>0) {
        const double top=c.work_function_eV-lowering,top_width=1/barrier->top_slope;
        for(double scale:{-48.,-8.,-1.,0.,1.,8.,48.})for(double value:{scale*decay,top+scale*top_width})
            if(value>low && value<high)breaks.push_back(value);
    }
    std::sort(breaks.begin(),breaks.end());breaks.erase(std::unique(breaks.begin(),breaks.end()),breaks.end());
    const auto [x,w]=legendre(8);std::vector<double> logweights,energies;
    for(std::size_t segment=1;segment<breaks.size();++segment)for(int cell=0;cell<c.energy_panels;++cell) {
        const double a=breaks[segment-1]+(breaks[segment]-breaks[segment-1])*cell/c.energy_panels;
        const double b=breaks[segment-1]+(breaks[segment]-breaks[segment-1])*(cell+1)/c.energy_panels;
        for(std::size_t q=0;q<x.size();++q) {
            const double E=(a+b)/2+(b-a)/2*x[q],supply=log_supply(E,kT);
            if(!std::isfinite(supply))continue;
            const double logD=barrier?log_transmission(barrier->gamow(E)):0;
            energies.push_back(E);logweights.push_back(std::log((b-a)/2*w[q])+supply+logD);
        }
    }
    require(!energies.empty(),"Empty numerical emission quadrature");
    const double largest=*std::max_element(logweights.begin(),logweights.end());long double sum=0;
    for(double l:logweights)sum+=std::exp(l-largest);
    const double logIntegral=largest+std::log(double(sum));
    d.log_current_A_m2=std::log(ideal_richardson_A_m2_K2()/(kb_eV*kb_eV))+logIntegral;d.has_log_current=true;
    d.current_A_m2=std::exp(d.log_current_A_m2);d.current_underflow=d.current_A_m2==0;
    require(std::isfinite(d.current_A_m2),"Emission current overflow");
    long double cdf=0,moment=0;
    for(std::size_t i=0;i<energies.size();++i) {
        const double p=std::exp(logweights[i]-logIntegral);if(p==0)continue;
        const double transverse=mean_transverse(energies[i],kT);cdf+=p;moment+=static_cast<long double>(p)*(energies[i]+transverse);
        d.nodes.push_back({energies[i],p,double(cdf),transverse});
    }
    require(!d.nodes.empty() && std::abs(double(cdf)-1)<1e-10,"Spectrum normalization failed");d.nodes.back().cdf=1;
    d.mean_metal_energy_eV=double(moment);d.wall_heat_W_m2=-d.current_A_m2*d.mean_metal_energy_eV;
    require(std::isfinite(d.mean_metal_energy_eV) && std::isfinite(d.wall_heat_W_m2),"Emission energy moment overflow");return d;
}
double transverse_quantile_eV(double E,double T,double u) {
    uniform(u);require(std::isfinite(E) && std::isfinite(T) && T>=0,"Finite transverse distribution inputs required");
    if(T==0){require(E<0,"Cold spectrum must lie below Fermi level");return -E*u;}
    const double kT=kb_eV*T,eta=-E/kT;
    if(eta< -36)return -kT*std::log1p(-u);
    const double target=(1-u)*softplus(eta);
    const double log_expm1=target>36?target+std::log1p(-std::exp(-target)):std::log(std::expm1(target));
    const double transverse=-kT*log_expm1-E;
    require(std::isfinite(transverse) && transverse>=-1e-12,"Invalid emitted transverse energy");return std::max(0.,transverse);
}
PreparedFlux prepare(const Stage& s,const Config& c) {
    require(s.accepted_state>0 && s.mesh_revision>0 && s.material_revision>0 && s.patch_id>0 &&
            s.field_accepted_state==s.accepted_state && s.thermal_accepted_state==s.accepted_state,"Field/T/geometry must bind one immutable accepted revision");
    require(std::isfinite(s.t) && std::isfinite(s.dt) && s.dt>0 && std::isfinite(s.t+s.dt) &&
            s.stage_time>=s.t && s.stage_time<=s.t+s.dt && s.field_time==s.stage_time && s.thermal_time==s.stage_time,"Emission feedback clock/stage mismatch");
    require(s.material_temperature_phase_admitted && std::isfinite(s.area_m2) && s.area_m2>0,"Actual evolved T/phase and physical patch area must be admitted");
    require(s.locally_planar_barrier_admitted && s.local_equilibrium_electron_supply_admitted,"Planar barrier/local equilibrium supply approximation must be admitted");
    require(finite(s.total_B_T),"Solved one-sided local B must be finite");
    auto d=evaluate(s.T_K,extraction_field(s.total_E,s.outward_normal),c);
    const double number=d.current_A_m2*s.area_m2*s.dt/qe,heat=d.wall_heat_W_m2*s.area_m2*s.dt;
    require(std::isfinite(number) && std::isfinite(heat),"Integrated emission flux overflow");return {s,std::move(d),number,heat};
}
double required_SN_matching_distance_m(const Distribution& d) {
    require(d.extraction_field_V_m>0,"SN matching distance only exists for extracting field");double distance=std::sqrt(image_c/d.extraction_field_V_m);
    if(!d.nodes.empty()){const double A=d.config.work_function_eV-d.nodes.front().normal_energy_eV;
        if(A>d.barrier_lowering_eV)distance=std::max(distance,(A+std::sqrt((A-d.barrier_lowering_eV)*(A+d.barrier_lowering_eV)))/(2*d.extraction_field_V_m));}
    return distance;
}
double Ledger::residual_J() const{return wall_heat_J+outgoing_K_J-incoming_K_J+electron_reference_change_J+image_reference_change_J-resolved_matching_field_work_J;}
double kinetic_energy_J(const Vector& u){require(finite(u),"Nonfinite particle momentum");const double u2=dot(u,u);return me*u2/(std::sqrt(1+u2/(c_light*c_light))+1);}
Packet sample(const PreparedFlux& f,const MatchingPlane& p,double weight,double un,double ut,double ua,std::uint64_t id) {
    uniform(un);uniform(ut);uniform(ua);require(id>0 && std::isfinite(weight) && weight>0 && !f.spectrum.nodes.empty(),"Finite positive packet weight, ID and nonempty spectrum required");
    const auto state=plane(f.spectrum,p);
    const auto it=std::lower_bound(f.spectrum.nodes.begin(),f.spectrum.nodes.end(),un,[](const Node& n,double u){return n.cdf<u;});
    require(it!=f.spectrum.nodes.end(),"Invalid spectrum CDF");const double E=it->normal_energy_eV,transverse=transverse_quantile_eV(E,f.spectrum.T_K,ut),normal=E-state.U;
    require(normal>=0 && std::isfinite(normal+transverse) && normal+transverse>0,"Nonpositive/nonfinite kinetic energy at matching surface");
    const double K=(normal+transverse)*qe,proper_speed=std::sqrt(K*(K+2*me*c_light*c_light))/(me*c_light);
    const auto& n=f.stage.outward_normal;unit(n);Vector axis=std::abs(n[0])<.8?Vector{1,0,0}:Vector{0,1,0};auto tangent=cross(n,axis);
    const double length=std::sqrt(dot(tangent,tangent));for(double& v:tangent)v/=length;const auto other=cross(n,tangent);
    Vector u{};for(int k=0;k<3;++k)u[k]=proper_speed*(std::sqrt(normal/(normal+transverse))*n[k]+
        std::sqrt(transverse/(normal+transverse))*(std::cos(2*pi*ua)*tangent[k]+std::sin(2*pi*ua)*other[k]));
    const double scale=qe*weight;
    Ledger ledger{weight,-scale,scale,-(E+transverse)*scale,K*weight,0,f.spectrum.config.work_function_eV*scale,state.image_energy*scale,state.field_work*scale};
    require(finite(u) && dot(u,n)>=0 && std::isfinite(ledger.residual_J()),"Invalid emitted phase-space or energy ledger");
    return {id,f.stage.accepted_state,f.stage.patch_id,weight,E,transverse,u,ledger};
}
double source_log_current_timescale_s(const Distribution& previous,const Distribution& current,double dt) {
    require(std::isfinite(dt) && dt>0,"Finite positive source-history interval required");
    if(!previous.has_log_current && !current.has_log_current)return std::numeric_limits<double>::max();
    require(previous.has_log_current && current.has_log_current,"Unresolved source onset: quasi-steady history criterion not established");
    const double change=std::abs(current.log_current_A_m2-previous.log_current_A_m2);
    return change==0?std::numeric_limits<double>::max():dt/change;
}
ConditionalPacket sample_conditional(const PreparedFlux& f,const ConditionalBoundary& policy,double weight,double un,double ut,double ua,std::uint64_t id) {
    uniform(un);uniform(ut);uniform(ua);
    require(!policy.PIC_already_applies_microscopic_single_electron_image_force,"Conditional restored-image model cannot duplicate an existing microscopic self-image force");
    if(f.spectrum.extraction_field_V_m<=0) {
        const MatchingPlane surface{PlaneKind::ThermalVacuumSurfaceStep,0,true,true,true};
        return {sample(f,surface,weight,un,ut,ua,id)};
    }
    require(std::isfinite(policy.local_linear_vacuum_extent_m) && policy.local_linear_vacuum_extent_m>0 &&
            std::isfinite(policy.source_log_current_timescale_s) && policy.source_log_current_timescale_s>0 &&
            std::isfinite(policy.image_force_fraction) && policy.image_force_fraction>0 && policy.image_force_fraction<=.2 &&
            std::isfinite(policy.image_energy_fraction) && policy.image_energy_fraction>0 && policy.image_energy_fraction<=.2 &&
            std::isfinite(policy.maximum_acceleration_time_fraction) && policy.maximum_acceleration_time_fraction>0 && policy.maximum_acceleration_time_fraction<=.2,
            "Conditional boundary requires actual linear-vacuum extent, source timescale and bounded image tolerances");
    require(!f.spectrum.nodes.empty() && std::isfinite(weight) && weight>0 && id>0,"Positive packet weight/ID and nonempty spectrum required");
    const auto it=std::lower_bound(f.spectrum.nodes.begin(),f.spectrum.nodes.end(),un,[](const Node& n,double u){return n.cdf<u;});
    require(it!=f.spectrum.nodes.end(),"Invalid normal spectrum CDF");
    const double E=it->normal_energy_eV,transverse=transverse_quantile_eV(E,f.spectrum.T_K,ut);
    const double F=f.spectrum.extraction_field_V_m,phi=f.spectrum.config.work_function_eV,A=phi-E;
    const double zmax=std::sqrt(image_c/F),lowering=f.spectrum.barrier_lowering_eV;
    const double exit=A>lowering?(A+std::sqrt((A-lowering)*(A+lowering)))/(2*F):zmax;
    const double scale=std::max(kb_eV*f.spectrum.T_K,f.spectrum.fermi_decay_width_eV);
    require(std::isfinite(scale) && scale>0,"No resolved source energy scale for image matching");
    const double tolerance=policy.image_energy_fraction*scale;
    const double L=std::max({exit,std::max(0.,A)/F,std::sqrt(image_c/(policy.image_force_fraction*F)),image_c/tolerance});
    require(std::isfinite(L) && L<=policy.local_linear_vacuum_extent_m,"Calculated subgrid layer exceeds actual locally linear vacuum; refine field/geometry or conditional tolerance");
    const double timescale=std::sqrt(2*me*L/(qe*F)); // acceleration scale, NOT a quantum tunnelling delay/upper bound
    const double fraction=timescale/policy.source_log_current_timescale_s;
    require(std::isfinite(timescale) && fraction<=policy.maximum_acceleration_time_fraction,"Quasi-steady matching acceleration/source-variation criterion failed");
    const double gyroangle=qe/me*std::sqrt(dot(f.stage.total_B_T,f.stage.total_B_T))*timescale;
    require(std::isfinite(policy.maximum_straight_projection_gyroangle) && policy.maximum_straight_projection_gyroangle>0 &&
            gyroangle<=policy.maximum_straight_projection_gyroangle,"Straight matching projection fails local B gyro/deflection criterion; use time-resolved Lorentz-image layer");
    // D(E) already contains the image barrier. Restore its conservative tail
    // before entering an image-free PIC micro-layer representation: no image
    // heat term or second explicit self-image force is added to collective E.
    const double normal=E-phi+F*L,K=(normal+transverse)*qe;
    require(normal>=-1e-12 && K>0 && std::isfinite(K),"Asymptotic matching produced invalid kinetic energy");
    const double utangent=std::sqrt(2*qe*transverse/me);
    const double unormal=std::sqrt(std::max(0.,2*qe*normal/me)+std::pow(K/(me*c_light),2));
    const auto& n=f.stage.outward_normal;unit(n);Vector axis=std::abs(n[0])<.8?Vector{1,0,0}:Vector{0,1,0};auto tangent=cross(n,axis);
    const double length=std::sqrt(dot(tangent,tangent));for(double& v:tangent)v/=length;const auto other=cross(n,tangent);
    Vector u{},displacement{},charge_displacement{};
    for(int j=0;j<3;++j){u[j]=unormal*n[j]+utangent*(std::cos(2*pi*ua)*tangent[j]+std::sin(2*pi*ua)*other[j]);displacement[j]=L*n[j];charge_displacement[j]=-qe*weight*displacement[j];}
    const double charge=qe*weight;
    Ledger ledger{weight,-charge,charge,-(E+transverse)*charge,K*weight,0,phi*charge,0,F*L*charge};
    require(finite(u) && std::isfinite(ledger.residual_J()),"Conditional packet/ledger overflow");
    Packet packet{id,f.stage.accepted_state,f.stage.patch_id,weight,E,transverse,u,ledger};
    return {packet,displacement,charge_displacement,L,exit,image_c/(L*L*F),image_c/L,timescale,fraction,gyroangle};
}
Ledger collect_at_material_surface(double phi,double K,double weight,bool admitted) {
    require(std::isfinite(phi) && phi>0 && std::isfinite(K) && K>=0 && std::isfinite(weight) && weight>0 && admitted,
            "Actual fully absorbed material collection and finite reference/kinetic energies required");
    const double charge=qe*weight;
    return {-weight,charge,-charge,K*weight+phi*charge,0,K*weight,-phi*charge,0,0};
}
Ledger collect_absorbed(const Distribution& d,const MatchingPlane& p,double incoming_K_J,double weight,bool admitted) {
    require(admitted && std::isfinite(incoming_K_J) && incoming_K_J>=0 && std::isfinite(weight) && weight>0,
            "Collection accommodation/reflection/secondary energy ownership must be explicit");
    const auto state=plane(d,p);const double scale=qe*weight;
    Ledger ledger{-weight,scale,-scale,(incoming_K_J+state.U*qe)*weight,0,incoming_K_J*weight,-d.config.work_function_eV*scale,-state.image_energy*scale,-state.field_work*scale};
    require(std::isfinite(ledger.residual_J()),"Collection energy overflow");return ledger;
}
Ledger sum_packets(const std::vector<Packet>& packets) {
    Ledger sum;std::set<std::uint64_t> ids;std::uint64_t revision=0;
    for(const auto& p:packets){require(ids.insert(p.event_id).second && p.event_id>0 && p.accepted_state>0 && (!revision || revision==p.accepted_state),"Duplicate source event or mixed accepted revision");revision=p.accepted_state;
        const auto& a=p.ledger;sum.physical_electron_number+=a.physical_electron_number;sum.particle_charge_C+=a.particle_charge_C;sum.electrode_charge_C+=a.electrode_charge_C;
        sum.wall_heat_J+=a.wall_heat_J;sum.outgoing_K_J+=a.outgoing_K_J;sum.incoming_K_J+=a.incoming_K_J;sum.electron_reference_change_J+=a.electron_reference_change_J;
        sum.image_reference_change_J+=a.image_reference_change_J;sum.resolved_matching_field_work_J+=a.resolved_matching_field_work_J;}
    require(std::isfinite(sum.residual_J()) && std::isfinite(sum.physical_electron_number),"Batch emission ledger overflow");return sum;
}
}
