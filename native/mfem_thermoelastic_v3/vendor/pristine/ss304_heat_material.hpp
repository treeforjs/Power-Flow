#pragma once
#include "electrode_heat_feedback.hpp"
#include "material/ss304_thermal.hpp"
namespace electrode_heat {
// No default qualification or density: all assumptions are caller selections.
inline Material conditional_ss304(ss304::Qualification mode,ss304::DensityConvention convention,double rho_ref,const std::string& source_seal_hash){
 auto model=std::make_shared<const ss304::Model>(mode,convention,rho_ref);Material m;std::ostringstream identity;identity<<"nominal304_mode="<<int(mode)<<" rho_ref="<<std::hexfloat<<rho_ref;m.id=identity.str();m.phase=mode==ss304::Qualification::nominal_synthesis_1800?"conditional_equilibrium_solid_mushy_liquid":"evaluated_solid";m.data_hash=source_seal_hash;m.volume_convention="fixed_reference_volume_cp_as_cv_approximation";m.native_callback_contract="ss304_owner_v1_direct_callbacks_stable_small_increment_v1";
 for(double T:{300.,ss304::Model::t_solidus,ss304::Model::t_liquidus,model->max_temperature_k()})if(T<=model->max_temperature_k()&&(m.knots.empty()||T>m.knots.back().T)){auto v=model->evaluate(T);m.knots.push_back({T,v.energy_j_m3_reference,v.k_w_m_k,v.sigma_s_m});}
 // The low-temperature mode ends before the phase markers.
 if(m.knots.back().T!=model->max_temperature_k()){auto v=model->evaluate(model->max_temperature_k());m.knots.push_back({v.temperature_k,v.energy_j_m3_reference,v.k_w_m_k,v.sigma_s_m});}
 m.cv_fn=[model](double T){return model->evaluate(T).capacity_j_m3_reference_k;};m.k_fn=[model](double T){return model->evaluate(T).k_w_m_k;};m.sigma_fn=[model](double T){return model->evaluate(T).sigma_s_m;};m.h_fn=[model](double T){return model->evaluate(T).energy_j_m3_reference;};m.phase_fn=[model](double T){return model->evaluate(T).phase;};m.liquid_fraction_fn=[model](double T){return model->evaluate(T).liquid_fraction;};
 m.delta_h_fn=[model](double T,double old){model->evaluate(T);model->evaluate(old);if(T==old)return 0.;double sign=T>old?1.:-1.,a=std::min(T,old),stop=std::max(T,old),sum=0.;while(a<stop){double b=stop;if(a<ss304::Model::t_solidus)b=std::min(b,ss304::Model::t_solidus);else if(a<ss304::Model::t_liquidus)b=std::min(b,ss304::Model::t_liquidus);double width=b-a;
   if(width>1e-4*std::max(a,b))sum+=model->evaluate(b).energy_j_m3_reference-model->evaluate(a).energy_j_m3_reference;
   else{double mid=a+.5*width,offset=.5*width/std::sqrt(3.);double x=mid-offset,y=mid+offset;if(b==ss304::Model::t_solidus||b==ss304::Model::t_liquidus){if(x>=b)x=std::nextafter(b,a);if(y>=b)y=std::nextafter(b,a);}sum+=.5*width*(model->evaluate(x).capacity_j_m3_reference_k+model->evaluate(y).capacity_j_m3_reference_k);}
   a=b;
  }return sign*sum;};
 m.validate();return m;
}
// Direct reusable conductivity callback for the field owner's NEXT accepted stage.
class AcceptedSolidConductivity:public mfem::Coefficient {
 const Coupler& thermal;
public:explicit AcceptedSolidConductivity(const Coupler& state):thermal(state){}
 mfem::real_t Eval(mfem::ElementTransformation& tr,const mfem::IntegrationPoint& ip)override{mfem::Vector x(2);tr.Transform(ip,x);return thermal.electrical_conductivity_at(x[0],x[1]);}
};
// The field owner may assemble its mass/current constraint with this coefficient.
// It returns exactly zero in vacuum; there is no spurious vacuum conductivity.
class AcceptedMetalConductivity:public mfem::Coefficient {
 std::vector<const Coupler*> bodies;
public:
 explicit AcceptedMetalConductivity(std::vector<const Coupler*> thermal):bodies(std::move(thermal)){need(!bodies.empty(),"no solid conductivity owners");for(auto* p:bodies)need(p!=nullptr,"null thermal conductivity owner");}
 mfem::real_t Eval(mfem::ElementTransformation&tr,const mfem::IntegrationPoint&ip)override{mfem::Vector x(2);tr.Transform(ip,x);double result=0;bool found=false;const auto& clock=bodies[0]->state();for(auto* body:bodies){need(body->state().revision==clock.revision&&body->state().time==clock.time,"solid conductivity owners not at common accepted stage");const auto& g=body->geometry();if(x[0]>=g.x0&&x[0]<=g.x0+g.width&&x[1]>=g.y0&&x[1]<=g.y0+g.thickness){need(!found,"overlapping solid conductivity owners");found=true;result=body->electrical_conductivity_at(x[0],x[1]);}}return result;}
};

}
