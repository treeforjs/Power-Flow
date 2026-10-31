// Heat/source owner v7: the sealed thermal-field emission evaluate with the temperature domain extended to
// T_max <= 5000 K (API-v7.txt section 10). Generated from the sealed source
//   hall-coupled-transport-20261029/.../bindings/heat-v2/source-kernel/src/surface_emission.cpp
// (helpers of its anonymous namespace and evaluate(), verbatim) with one change: the domain check
// "T=0 or 1..3000 K" reads "T=0 or 1..T_max K" (3000 <= T_max <= 5000). For T <= 3000 K the result is the
// sealed evaluate's, bit for bit (same expressions and order; checked by test_surface_v7). Compile with the
// sealed kernel's flags (-O2 -ffp-contract=off -fno-fast-math).
#include "emission_hot_v7.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace electrode_heat_v7::emission_hot {
using namespace surface_emission;
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
}
Distribution evaluate(double T,double F,const Config& c,double T_max) {
    config(c);require(std::isfinite(T_max) && T_max>=3000 && T_max<=5000,"v7 extended emission domain bound must lie in 3000..5000 K");
    require(std::isfinite(T) && (T==0 || (T>=1 && T<=T_max)),"Numerical free-electron law domain is T=0 or1..T_max K (v7 extended, T_max<=5000K); no extrapolation/clamp");
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
}  // namespace electrode_heat_v7::emission_hot
