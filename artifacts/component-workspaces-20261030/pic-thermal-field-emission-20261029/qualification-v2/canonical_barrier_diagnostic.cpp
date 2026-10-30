#include "surface_emission.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <sched.h>
#include <sys/resource.h>

// Exploratory locally tangential-B first-barrier calculation, not an admitted
// source law. From minimal coupling (pz-eBx)^2/(2m), retaining the cross term.
// Source: Lebedynskyi et al., arXiv:1702.02713 Eqs.4,20,28. No current, supply,
// general oblique-B scattering, magnetic work, or terminal escape is inferred.
namespace se=surface_emission;
namespace {
constexpr double pi=3.14159265358979323846,hbar=se::planck/(2*pi),c=se::qe/(16*pi*se::eps0),phi=4.3;
void need(bool b,const char* s){if(!b)throw std::runtime_error(s);}
template<class Function> double root(Function f,double a,double b) {
    const bool sign=f(a)>0;
    need(sign!=(f(b)>0),"Unbracketed canonical barrier root");
    for(int i=0;i<90;++i){const double m=(a+b)/2;if((f(m)>0)==sign)a=m;else b=m;}
    return (a+b)/2;
}
struct Barrier {
    double F,B,v,alpha,Feff,xs,xmin,top,bottom;
    bool has_pocket=false;
    double U(double x)const{return phi-Feff*x-c/x+alpha*x*x;}
    double curvature(double x)const{return -2*c/(x*x*x)+2*alpha;}
    Barrier(double electric,double magnetic,double canonical_v):F(electric),B(magnetic),v(canonical_v),alpha(se::qe*B*B/(2*se::me)),Feff(F+B*v) {
        if(alpha==0){if(Feff<=0)return;xs=std::sqrt(c/Feff);xmin=std::numeric_limits<double>::infinity();top=U(xs);bottom=-std::numeric_limits<double>::infinity();has_pocket=true;return;}
        const double critical=std::cbrt(c/alpha);
        if(Feff<=3*std::cbrt(alpha*alpha*c))return;
        auto derivative=[&](double x){return -Feff+c/(x*x)+2*alpha*x;};
        xs=root(derivative,critical*1e-12,critical);
        xmin=root(derivative,critical,std::max(2*critical,Feff/(2*alpha)));
        top=U(xs);bottom=U(xmin);has_pocket=true;
    }
    double signed_action(double E,int panels,double& outer)const {
        if(!has_pocket || E<=bottom){outer=0;return std::numeric_limits<double>::infinity();}
        if(E>=top){outer=xs;return 2*pi*se::qe*(top-E)/(hbar*std::sqrt(se::qe*(-curvature(xs))/se::me));}
        const double left=root([&](double x){return U(x)-E;},xs*1e-12,xs);
        double right_limit=xmin;
        if(!std::isfinite(right_limit)){right_limit=2*xs;while(U(right_limit)>E)right_limit*=2;}
        outer=root([&](double x){return U(x)-E;},xs,right_limit);
        const double mid=(outer+left)/2,half=(outer-left)/2;
        long double sum=0;
        for(int i=1;i<panels;++i){const double t=pi*i/panels,x=mid-half*std::cos(t);const double integrand=std::sqrt(std::max(0.,U(x)-E))*half*std::sin(t);sum+=(i%2?4:2)*integrand;}
        return double(sum)*(pi/panels)/3*2*std::sqrt(2*se::me*se::qe)/hbar;
    }
};
double logD(double G){if(G>40)return -G;if(G< -40)return -std::exp(G);return -std::log1p(std::exp(G));}
}
int main(int argc,char** argv) {
    try {
        need(argc==2,"Expected diagnostic output directory");cpu_set_t mask;CPU_ZERO(&mask);
        need(sched_getaffinity(0,sizeof(mask),&mask)==0 && CPU_COUNT(&mask)==1 && CPU_ISSET(7,&mask) && getpriority(PRIO_PROCESS,0)>=15,"CPU7/nice15 required");
        const std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
        std::ofstream table(out/"canonical-tangential-B-diagnostic.csv");table<<std::setprecision(17)<<"F_V_m,B_tangent_T,metal_canonical_v_z_m_s,E_normal_eV,local_allowed_pocket,local_saddle_eV,pocket_minimum_eV,first_outer_root_m,G_SN,G_canonical,logD_SN,logD_canonical\n";
        double max_B0_relative_error=0,max_action_refinement=0;int pockets=0,missing=0;
        for(double F:{213175.33729099843,1e7,97724679.40458213})for(double B:{0.,30.})for(double v:{-3e5,0.,3e5}) {
            const double E=phi-2*std::sqrt(c*F)-.05,G0=se::action(E,F,se::Config{phi});
            Barrier b(F,B,v);double x=0;const double G=b.signed_action(E,2048,x);double fine_x=0;const double fine=b.signed_action(E,4096,fine_x);
            const bool pocket=std::isfinite(G);if(pocket){++pockets;max_action_refinement=std::max(max_action_refinement,std::abs(G-fine)/std::max(1.,std::abs(fine)));}else ++missing;
            if(B==0){max_B0_relative_error=std::max(max_B0_relative_error,std::abs(G-G0)/std::max(1.,std::abs(G0)));need(std::abs(G-G0)<1e-7*std::max(1.,std::abs(G0)),"Canonical B0 failed sealed SN action comparison");}
            table<<F<<','<<B<<','<<v<<','<<E<<','<<pocket<<',';
            if(b.has_pocket)table<<b.top<<','<<b.bottom;else table<<"nan,nan";
            table<<','<<x<<','<<G0<<','<<G<<','<<logD(G0)<<','<<logD(G)<<'\n';
        }
        need(max_action_refinement<1e-7,"Canonical action refinement failed");
        need(missing>0,"Diagnostic did not exercise the absence of a local allowed pocket");
        std::ofstream report(out/"canonical-barrier-diagnostic.json");
        report<<std::setprecision(17)<<"{\"diagnostic_checks_passed\":true,\"source_law_admitted\":false,\"stable_API_modified\":false,\"nominal_guard_widened\":false,\"B0_action_max_relative_error\":"<<max_B0_relative_error<<",\"action_refinement_max_relative_error\":"<<max_action_refinement<<",\"local_pocket_cases\":"<<pockets<<",\"no_local_pocket_cases\":"<<missing<<",\"scope\":\"Minimal-coupling first-barrier sensitivity for purely tangential B only. Not finite-T supply or joint sampling; not a total escaping current or a model for oblique B. Local allowed pockets require resolved subsequent Lorentz/image returns.\",\"primary_source\":\"https://arxiv.org/abs/1702.02713\"}\n";
        std::cout<<"Canonical barrier diagnostic checked; not an admitted source law\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
