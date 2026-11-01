#include "magnetic_guard_assessment.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace magnetized_emission {
namespace se=surface_emission;
namespace {
constexpr double pi=3.14159265358979323846,hbar=se::planck/(2*pi),image_c=se::qe/(16*pi*se::eps0),kb_eV=se::kb/se::qe;
double dot(const se::Vector& a,const se::Vector& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double loglog1pexp(double y){if(y < -40)return y;if(y>40)return std::log(y);return std::log(std::log1p(std::exp(y)));}
double logS(double E,double t){return std::log(t)+loglog1pexp(-E/t);}
double F1negative(double y) {
    // -Li2(-exp(y)) = Li2(exp(y)/(1+exp(y))) + log(1+exp(y))^2/2.
    // The series argument is at most1/2, including the degenerate endpoint.
    const double q=std::exp(y),r=q/(1+q),l=std::log1p(q);long double sum=0,power=r;
    for(int k=1;k<=80;++k){sum+=power/(k*k);power*=r;}
    return double(sum)+l*l/2;
}
double logF1(double y){
    if(y < -40)return y;
    if(y<=0)return std::log(F1negative(y));
    return std::log(y*y/2+pi*pi/6-F1negative(-y));
}
double logadd(double a,double b){if(a==-std::numeric_limits<double>::infinity())return b;if(b==-std::numeric_limits<double>::infinity())return a;const double hi=std::max(a,b);return hi+std::log1p(std::exp(std::min(a,b)-hi));}
}
RejectedSupply assess_rejected_supply(const se::Distribution& d,const se::Vector& B,const se::Vector& n,double eta) {
    if(!(eta>0 && eta<=.2) || !std::all_of(B.begin(),B.end(),[](double x){return std::isfinite(x);}) || std::abs(dot(n,n)-1)>1e-12)
        throw std::invalid_argument("Invalid magnetic guard assessment inputs");
    RejectedSupply out;
    if(d.extraction_field_V_m<=0)return out; // thermal-step branch is distinct, not proof of zero magnetic effect
    out.SchottkyNordheim_guard_applies=true;
    const double scale=std::max(kb_eV*d.T_K,d.fermi_decay_width_eV),Bn=dot(B,n);
    double Bt2=0;for(int j=0;j<3;++j)Bt2+=std::pow(B[j]-Bn*n[j],2);
    const double Bt=std::sqrt(Bt2),t=kb_eV*d.T_K,F=d.extraction_field_V_m,lowering=d.barrier_lowering_eV;
    if(!(scale>0))throw std::invalid_argument("Missing v1 source energy scale");
    long double rejected=0,energy=0,probability=0;double log_probability=-std::numeric_limits<double>::infinity();
    for(const auto& node:d.nodes) {
        probability+=node.probability;if(node.probability==0)continue;
        const double E=node.normal_energy_eV,A=d.config.work_function_eV-E+eta*scale;
        const double x=A>lowering?(A+std::sqrt((A-lowering)*(A+lowering)))/(2*F):std::sqrt(image_c/F);
        const double room=.2*scale-se::qe*Bt2*x*x/(2*se::me)-.5*hbar*std::abs(Bn)/se::me;
        double cut=0;
        if(room>0){if(Bt==0)continue;const double vmax=room/(Bt*x);cut=se::me*vmax*vmax/(2*se::qe);}
        double p=0,moment=0,logp=-std::numeric_limits<double>::infinity();
        if(d.T_K==0) {
            if(E<0 && cut< -E){p=1-cut/(-E);logp=std::log(p);moment=-(E+cut)*(E+cut)/(2*(-E));}
        } else {
            const double a=E+cut,denominator=logS(E,t);logp=std::min(0.,logS(a,t)-denominator);p=std::exp(logp);
            moment=a*p+std::exp(2*std::log(t)+logF1(-a/t)-denominator);
        }
        rejected+=node.probability*p;energy+=node.probability*moment;
        log_probability=logadd(log_probability,std::log(node.probability)+logp);
    }
    if(std::abs(probability-1)>1e-10)throw std::runtime_error("Unnormalized immutable v1 normal spectrum");
    out.rejected_fraction=std::clamp(double(rejected),0.,1.);out.rejected_energy_moment_eV=double(energy);
    out.has_log_rejected_fraction=std::isfinite(log_probability);out.log_rejected_fraction=log_probability;
    out.rejected_current_A_m2=d.current_A_m2*out.rejected_fraction;
    out.rejected_complete_wall_heat_W_m2=-d.current_A_m2*out.rejected_energy_moment_eV;
    if(d.has_log_current && out.has_log_rejected_fraction){out.has_log_rejected_current=true;out.log_rejected_current_A_m2=d.log_current_A_m2+log_probability;}
    return out;
}
}
