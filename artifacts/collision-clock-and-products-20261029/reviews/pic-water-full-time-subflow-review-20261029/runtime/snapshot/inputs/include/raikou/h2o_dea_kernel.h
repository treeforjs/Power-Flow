#pragma once
// Source-grounded, explicitly conditional6.5eV DEA comparison. No trajectories,
// inventories, random engine, rate scheduling or accepted-state commit here.
#include "raikou/h2o_conditional_parent.h"
#include <numeric>

namespace raikou::kinetics {
struct HminusDeaProducts {
    CollisionMomentum Hminus_momentum_kg_m_s{},OH_momentum_kg_m_s{};
    double Hminus_kinetic_J{0},OH_kinetic_J{0},OH_state_mass_kg{0};
    double chemical_reference_gain_J{0},OH_internal_gain_J{0},physical_weight_each{0},cross_section_m2{0};
    int OH_vibrational_level{0};double fragment_cosine{0};
};

namespace h2o_dea_detail {
inline const Json& model(const StrictCollisionDataset& data) {
    const auto s=data.kernel_status("H2O:attachment_yield:Hminus","Hminus_Ram65_factorized_window_v1");
    if (s.at("status")!="conditional_available") throw CollisionCoverageError("Ram6.5eV DEA model not present in this frozen release");
    // Return a reference into the original immutable dataset, not the status copy.
    for (const auto& p:data.metadata().at("conditional_event_models")) if (p.at("id")=="Hminus_Ram65_factorized_window_v1")return p;
    throw CollisionCoverageError("Missing DEA model");
}
inline void require(bool value,const char* text) {if (!value)throw std::invalid_argument(text);}
inline std::array<double,3> coefficients(const Json& p) {
    const auto& a=p.at("angular");const double x=a.at("a").get<double>(),y=a.at("b").get<double>(),d=a.at("phase_rad").get<double>();
    return {2*x*x+.96*y*y,1.2*y*y,4.2*x*y*std::cos(d)};
}
inline double primitive(double mu,const std::array<double,3>& q) {
    const double m2=mu*mu,m3=m2*mu,m4=m2*m2,m5=m4*mu;
    return q[0]*(mu-m3/3)+q[1]*(m3/3-m5/5)+q[2]*(m2/2-m4/4);
}
} // namespace h2o_dea_detail

[[nodiscard]] inline double hminus_Ram65_cosine_density(const StrictCollisionDataset& data,double mu) {
    h2o_dea_detail::require(std::isfinite(mu)&&mu>=-1&&mu<=1,"Invalid fragment cosine");
    const auto q=h2o_dea_detail::coefficients(h2o_dea_detail::model(data));
    const double norm=4*q[0]/3+4*q[1]/15;
    return (1-mu*mu)*(q[0]+q[1]*mu*mu+q[2]*mu)/norm;
}
[[nodiscard]] inline double hminus_Ram65_cosine_cdf(const StrictCollisionDataset& data,double mu) {
    h2o_dea_detail::require(std::isfinite(mu)&&mu>=-1&&mu<=1,"Invalid fragment cosine");
    const auto q=h2o_dea_detail::coefficients(h2o_dea_detail::model(data));
    const double norm=4*q[0]/3+4*q[1]/15;
    return (h2o_dea_detail::primitive(mu,q)-h2o_dea_detail::primitive(-1,q))/norm;
}
[[nodiscard]] inline double hminus_Ram65_vibrational_eV(const StrictCollisionDataset& data,int level) {
    h2o_dea_detail::require(level>=0&&level<=4,"Unsupported OH vibrational state");
    const auto& p=h2o_dea_detail::model(data);const auto a=p.at("OH_X_G_cm1_coefficients").get<std::vector<double>>();
    const auto G=[&](double q) {double answer=0,power=q;for (const double v:a) {answer+=v*power;power*=q;}return answer;};
    return (G(level+.5)-G(.5))*p.at("h_c_eV_cm").get<double>();
}
[[nodiscard]] inline HminusDeaProducts conditional_hminus_Ram65_products(
    const StrictCollisionDataset& data,const std::string& explicitly_selected_model,
    double incident_target_rest_eV,double represented_weight,double neutral_mass_kg,double Hminus_mass_kg,
    const CollisionMomentum& incident_unit_direction,double level_uniform,double cosine_uniform,double azimuth_uniform) {
    using h2o_dea_detail::require;
    constexpr double e=1.602176634e-19,me=9.1093837139e-31,c=299792458.;
    require(explicitly_selected_model=="Hminus_Ram65_factorized_window_v1","Explicit Ram6.5eV model opt-in required");
    data.require_MCC_event_kernel("H2O:attachment_yield:Hminus",explicitly_selected_model,incident_target_rest_eV);
    require(std::isfinite(represented_weight)&&represented_weight>0,"Invalid represented weight");
    require(std::isfinite(neutral_mass_kg)&&std::isfinite(Hminus_mass_kg)&&neutral_mass_kg>Hminus_mass_kg&&Hminus_mass_kg>100*me,"Invalid fragment mass basis");
    for (const double u:{level_uniform,cosine_uniform,azimuth_uniform})require(std::isfinite(u)&&u>=0&&u<1,"Host RNG variates must be in[0,1)");
    double norm=0;for (const auto v:incident_unit_direction) {require(std::isfinite(v),"Invalid incident axis");norm+=v*v;}require(std::abs(norm-1)<1e-12,"Incident direction must be unit length");
    const auto& p=h2o_dea_detail::model(data);const auto weights=p.at("OH_vibrational_relative_intensities").get<std::vector<double>>();
    const double total=std::accumulate(weights.begin(),weights.end(),0.);double cumulative=0;int level=4;
    for (int v=0;v<5;++v) {cumulative+=weights[v]/total;if (level_uniform<cumulative) {level=v;break;}}
    const auto angular=h2o_dea_detail::coefficients(p);const double angular_norm=4*angular[0]/3+4*angular[1]/15;
    const auto F=[&](double mu) {return (h2o_dea_detail::primitive(mu,angular)-h2o_dea_detail::primitive(-1,angular))/angular_norm;};
    double lo=-1,hi=1;for (int i=0;i<56;++i) {const double mid=(lo+hi)/2;if (F(mid)<cosine_uniform)lo=mid;else hi=mid;}
    const double mu=(lo+hi)/2,phi=2*std::numbers::pi*azimuth_uniform;
    CollisionMomentum first{};const int component=std::abs(incident_unit_direction[0])<.8?0:1;
    first[component]=1;double dot=incident_unit_direction[component];norm=0;
    for (int d=0;d<3;++d) {first[d]-=dot*incident_unit_direction[d];norm+=first[d]*first[d];}
    for (auto& v:first)v/=std::sqrt(norm);
    const CollisionMomentum second{incident_unit_direction[1]*first[2]-incident_unit_direction[2]*first[1],incident_unit_direction[2]*first[0]-incident_unit_direction[0]*first[2],incident_unit_direction[0]*first[1]-incident_unit_direction[1]*first[0]};
    const double reference=p.at("chemical_reference_gain_eV").get<double>()*e,internal=hminus_Ram65_vibrational_eV(data,level)*e;
    const double K=incident_target_rest_eV*e,D=(neutral_mass_kg+me)*c*c;
    const double W=std::sqrt(D*D+2*neutral_mass_kg*c*c*K),available=2*neutral_mass_kg*c*c*K/(W+D);
    const double Q=available-reference-internal;require(Q>0,"Selected source state not energetically available; no resampling/renormalization");
    const double mass2=neutral_mass_kg+me+(reference+internal)/(c*c)-Hminus_mass_kg;
    require(mass2>0,"Invalid OH state mass");
    const double k1=Q*(Q+2*mass2*c*c)/(2*W),k2=Q*(Q+2*Hminus_mass_kg*c*c)/(2*W);
    const double mom=std::sqrt(k1*(k1+2*Hminus_mass_kg*c*c))/c;
    const double p0=std::sqrt(K*(K+2*me*c*c))/c,beta=p0*c/(D+K),gamma=1/std::sqrt(1-beta*beta),gm1=beta*beta*gamma*gamma/(gamma+1);
    const double E1=Hminus_mass_kg*c*c+k1,E2=mass2*c*c+k2;
    // The measured fit is a LAB fragment angle. Invert the two-body boost
    // instead of silently treating that fit as a center-of-mass angle.
    require(mom>beta*E1/c,"Lab angular inversion is not one-to-one in this mass/energy basis");
    const auto lab_cosine=[&](double x) {const double parallel=gamma*(mom*x+beta*E1/c);
        return parallel/std::sqrt(parallel*parallel+mom*mom*(1-x*x));};
    lo=-1;hi=1;for (int i=0;i<56;++i) {const double mid=(lo+hi)/2;if (lab_cosine(mid)<mu)lo=mid;else hi=mid;}
    const double cm_mu=(lo+hi)/2,cm_sine=std::sqrt(std::max(0.,1-cm_mu*cm_mu));
    CollisionMomentum axis{};for (int d=0;d<3;++d)axis[d]=cm_mu*incident_unit_direction[d]+cm_sine*(std::cos(phi)*first[d]+std::sin(phi)*second[d]);
    HminusDeaProducts out;out.OH_vibrational_level=level;out.fragment_cosine=mu;out.OH_state_mass_kg=mass2;
    out.chemical_reference_gain_J=reference;out.OH_internal_gain_J=internal;out.physical_weight_each=represented_weight;
    out.cross_section_m2=data.observable_value_SI("H2O:attachment_yield:Hminus",incident_target_rest_eV);
    out.Hminus_kinetic_J=k1+gm1*E1+gamma*beta*c*mom*cm_mu;
    out.OH_kinetic_J=k2+gm1*E2-gamma*beta*c*mom*cm_mu;
    for (int d=0;d<3;++d) {
        out.Hminus_momentum_kg_m_s[d]=mom*axis[d]+(gm1*mom*cm_mu+gamma*beta*E1/c)*incident_unit_direction[d];
        out.OH_momentum_kg_m_s[d]=-mom*axis[d]+(-gm1*mom*cm_mu+gamma*beta*E2/c)*incident_unit_direction[d];
    }
    return out;
}
} // namespace raikou::kinetics
