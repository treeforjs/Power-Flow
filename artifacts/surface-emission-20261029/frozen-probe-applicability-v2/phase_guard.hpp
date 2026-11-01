#pragma once
#include "surface_emission.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
namespace original_v1_guard {
namespace se=surface_emission;
constexpr double pi=3.14159265358979323846;
struct PhaseBoundary { double x_m=0,scale_eV=0,Bn_T=0,Bt_T=0,room_eV=0,maximum_transverse_eV=0; bool thermal_step=false,all_transverse=false; };
inline double dot(const se::Vector&a,const se::Vector&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline PhaseBoundary boundary(const se::Distribution&d,const se::Vector&B,const se::Vector&n,double E){
 PhaseBoundary p;p.thermal_step=d.extraction_field_V_m<=0;if(p.thermal_step){p.all_transverse=true;p.maximum_transverse_eV=INFINITY;return p;}
 p.scale_eV=std::max(se::kb*d.T_K/se::qe,d.fermi_decay_width_eV);p.Bn_T=dot(B,n);double bt2=0;for(int i=0;i<3;++i)bt2+=std::pow(B[i]-p.Bn_T*n[i],2);p.Bt_T=std::sqrt(bt2);
 const double F=d.extraction_field_V_m,A=d.config.work_function_eV-E+.03*p.scale_eV,L=d.barrier_lowering_eV,c=se::qe/(16*pi*se::eps0);
 p.x_m=A>L?(A+std::sqrt((A-L)*(A+L)))/(2*F):std::sqrt(c/F);
 p.room_eV=.2*p.scale_eV-se::qe*bt2*p.x_m*p.x_m/(2*se::me)-.5*(se::planck/(2*pi))*std::abs(p.Bn_T)/se::me;
 if(p.room_eV>=0&&p.Bt_T==0){p.all_transverse=true;p.maximum_transverse_eV=INFINITY;}
 else if(p.room_eV>0)p.maximum_transverse_eV=se::me/(2*se::qe)*std::pow(p.room_eV/(p.Bt_T*p.x_m),2);
 return p;
}
inline double logsoftplus(double y){if(y< -40)return y;if(y>40)return std::log(y);return std::log(std::log1p(std::exp(y)));}
// Conditional FD fraction admitted at a fixed normal energy. Finite-temperature
// tails remain nonzero; no rejected phase is dropped or renormalized in a source.
inline double admitted_fraction(const PhaseBoundary&p,double E,double T){
 if(p.all_transverse)return 1;if(p.room_eV<=0)return 0;
 if(T==0)return std::min(1.,p.maximum_transverse_eV/(-E));
 const double t=se::kb*T/se::qe,l=std::min(0.,logsoftplus(-(E+p.maximum_transverse_eV)/t)-logsoftplus(-E/t));
 return -std::expm1(l);
}
inline double ratio(const PhaseBoundary&p,double Et){if(p.thermal_step)return 0;return (p.Bt_T*p.x_m*std::sqrt(2*se::qe*Et/se::me)+se::qe*p.Bt_T*p.Bt_T*p.x_m*p.x_m/(2*se::me)+.5*(se::planck/(2*pi))*std::abs(p.Bn_T)/se::me)/p.scale_eV;}
}
