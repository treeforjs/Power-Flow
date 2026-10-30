#pragma once
#include "electrode_heat_feedback.hpp"
#include "binding_seals.hpp"
#include "material-v2/ss304_thermal_v2.hpp"
namespace electrode_heat_v2 {
inline Material conditional_ss304_v2(ss304v2::Qualification mode,ss304::DensityConvention convention,double rho_ref,const std::string& seal) {
 auto model=std::make_shared<const ss304v2::Model>(mode,convention,rho_ref);Material m;std::ostringstream id;id<<"nominal304_v2_q="<<int(mode)<<" rho_ref="<<std::hexfloat<<rho_ref;m.id=id.str();m.phase="explicit_conditional_equilibrium_fixed_geometry";m.data_hash=seal;m.volume_convention="fixed_reference_volume_cp_as_cv_approximation";m.native_callback_contract=std::string("ss304_owner_v2_direct_callbacks_stable_increment_v2/")+material_headers_sha256;
 for(double T:{298.2,300.,1500.,1672.,1727.,1800.,3000.})if(T<=model->max_temperature_k()){auto v=model->evaluate(T);m.knots.push_back({T,v.energy_j_m3_reference,v.k_w_m_k,v.sigma_s_m});}
 m.cv_fn=[model](double T){return model->evaluate(T).capacity_j_m3_reference_k;};m.k_fn=[model](double T){return model->evaluate(T).k_w_m_k;};m.sigma_fn=[model](double T){return model->evaluate(T).sigma_s_m;};m.h_fn=[model](double T){return model->evaluate(T).energy_j_m3_reference;};m.phase_fn=[model](double T){return model->evaluate(T).phase;};m.liquid_fraction_fn=[model](double T){return model->evaluate(T).liquid_fraction;};
 m.delta_h_fn=[model](double T,double old){model->evaluate(T);model->evaluate(old);if(T==old)return 0.;double sign=T>old?1.:-1.,a=std::min(T,old),end=std::max(T,old),result=0.;while(a<end){double b=end;for(double cut:{300.,1672.,1727.,1800.})if(a<cut){b=std::min(b,cut);break;}double width=b-a;if(width>1e-4*std::max(a,b))result+=model->evaluate(b).energy_j_m3_reference-model->evaluate(a).energy_j_m3_reference;else{double mid=a+.5*width,off=.5*width/std::sqrt(3.);double x=mid-off,y=mid+off;if(x>=b)x=std::nextafter(b,a);if(y>=b)y=std::nextafter(b,a);result+=.5*width*(model->evaluate(x).capacity_j_m3_reference_k+model->evaluate(y).capacity_j_m3_reference_k);}a=b;}return sign*result;};m.validate();return m;
}
}
