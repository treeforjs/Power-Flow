#include "include/ss304_thermal_v2.hpp"
#include <string>
namespace {
thread_local std::string error;
ss304v2::Model material(int q,double rho) {
  return ss304v2::Model(static_cast<ss304v2::Qualification>(q),
    ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,rho);
}
template<class F> int guarded(F f) { try { f();error.clear();return 0; }
  catch(const std::exception& e) {error=e.what();return 1;} }
}
extern "C" {
const char* ss304v2_error() {return error.c_str();}
int ss304v2_evaluate(int q,double rho,double t,ss304v2::Values* out) {
 return guarded([&]{if(!out)throw std::invalid_argument("null output");const auto v=material(q,rho).evaluate(t);*out=v;});
}
int ss304v2_inverse(int q,double rho,double e,int vol,double* out) {
 return guarded([&]{if(!out||vol<0||vol>1)throw std::invalid_argument("invalid inverse output/measure");const auto m=material(q,rho);const double t=vol?m.temperature_from_reference_energy(e):m.temperature_from_u(e);*out=t;});
}
int ss304v2_vapor(double t,ss304v2::VaporValues* out) {
 return guarded([&]{if(!out)throw std::invalid_argument("null output");const auto m=material(2,8000);const auto v=ss304v2::NeutralVapor(m).evaluate(t);*out=v;});
}
int ss304v2_flux(double t,double back,double alpha,double* out) {
 return guarded([&]{if(!out)throw std::invalid_argument("null output");const auto m=material(2,8000);const double f=ss304v2::NeutralVapor(m).net_mass_flux(t,back,alpha);*out=f;});
}
int ss304v2_gas_k(double t,double diameter,double ionized_fraction,double* out) {
 return guarded([&]{if(!out)throw std::invalid_argument("null output");const auto m=material(2,8000);const double k=ss304v2::NeutralVapor(m).neutral_gas_k(t,diameter,ionized_fraction);*out=k;});
}
int ss304v2_transfer(double rho,const ss304v2::MetalState* base,std::uint64_t revision,std::uint64_t event,
 double dm,double kinetic,double heat,ss304v2::TransferTrial* out) {
 return guarded([&]{if(!base||!out)throw std::invalid_argument("null transfer state");const auto m=material(2,rho);const ss304v2::NeutralVapor v(m);const auto trial=ss304v2::prepare_transfer(m,v,*base,revision,event,dm,kinetic,heat);*out=trial;});
}
int ss304v2_energy_support(int q,double rho,double mass,double energy,ss304v2::EnergySupport* out) {
 return guarded([&]{if(!out)throw std::invalid_argument("null output");const auto v=material(q,rho).assess_energy(mass,energy);*out=v;});
}
int ss304v2_evaluate_inventory(int q,double rho,double t,double mass,double volume,ss304v2::Values* out) {
 return guarded([&]{if(!out)throw std::invalid_argument("null output");const auto v=material(q,rho).evaluate_inventory(t,mass,volume);*out=v;});
}
}
