#pragma once
#include "electrode_heat_feedback_v3.hpp"
#include "material_v2/ss304_thermal_v2.hpp"

namespace electrode_heat {

inline constexpr const char* ss304_v2_source_seal_sha256() {
  return "1a1845e04f1f6fda33338d9bd284f5e1ead9cbc4b7747d0c6355e32cb5661c2d";
}

// This adapter is deliberately q1 only. No q2 liquid/vapor scenario is selected.
// The caller supplies accepted mass/reference volume as rho_ref. Material
// callbacks are volumetric already; the heat solver applies reference volume
// once. Mechanical773.15K and smooth-segment guards remain independent.
inline Material conditional_ss304_v2(ss304v2::Qualification mode,
    ss304::DensityConvention convention, double rho_ref,
    const std::string& source_seal_hash) {
  need(mode == ss304v2::Qualification::nominal_1800,
       "v2 heat adapter requires explicit q1 nominal_1800");
  need(convention == ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,
       "v2 heat adapter requires explicit reference-volume caloric convention");
  need(source_seal_hash == ss304_v2_source_seal_sha256(),
       "v2 heat adapter source seal differs from reviewed compiled model");
  const auto model = std::make_shared<const ss304v2::Model>(mode, convention, rho_ref);
  Material material;
  std::ostringstream identity;
  identity << "nominal304_v2_q1_cold298p2_reference300 rho_ref=" << std::hexfloat << rho_ref;
  material.id = identity.str();
  material.phase = "conditional_equilibrium_solid_mushy_liquid";
  material.data_hash = source_seal_hash;
  material.volume_convention = "fixed_reference_volume_cp_as_cv_approximation";
  material.native_callback_contract = "ss304_owner_v2_q1_direct_callbacks_exact_increment_v1";
  for (const double t : {ss304v2::Model::minimum_temperature_k, 300.,
                        ss304::Model::t_solidus, ss304::Model::t_liquidus, 1800.}) {
    const auto v = model->evaluate(t);
    material.knots.push_back({t, v.energy_j_m3_reference, v.k_w_m_k, v.sigma_s_m});
  }
  material.cv_fn = [model](double t) { return model->evaluate(t).capacity_j_m3_reference_k; };
  material.k_fn = [model](double t) { return model->evaluate(t).k_w_m_k; };
  material.sigma_fn = [model](double t) { return model->evaluate(t).sigma_s_m; };
  material.h_fn = [model](double t) { return model->evaluate(t).energy_j_m3_reference; };
  material.phase_fn = [model](double t) { return model->evaluate(t).phase; };
  material.liquid_fraction_fn = [model](double t) { return model->evaluate(t).liquid_fraction; };

  // Stable exact integral of the sealed q1 capacity. In the solid branch:
  // Delta h=(b-a)[472+.068(a+b)-2820000/(a*b)]. This includes298.2-300K
  // and avoids subtracting large reference energies for tiny increments.
  // The mushy branch integrates its affine sensible capacity plus fusion once;
  // the liquid branch has constant800J/(kg K). No duplicated latent debit.
  material.delta_h_fn = [model](double t, double old) {
    (void)model->evaluate(t); (void)model->evaluate(old);
    if (t == old) return 0.;
    const double sign = t > old ? 1. : -1.;
    const double solidus = ss304::Model::t_solidus;
    const double liquidus = ss304::Model::t_liquidus;
    const double width = liquidus-solidus;
    const double cp_solidus = model->evaluate(solidus).cp_sensible_j_kg_k;
    double a = std::min(t, old), stop = std::max(t, old), sum = 0.;
    while (a < stop) {
      double b = stop, increment = 0.;
      if (a < solidus) {
        b = std::min(b, solidus);
        increment = (b-a)*(472.+.068*(a+b)-2820000./(a*b));
      } else if (a < liquidus) {
        b = std::min(b, liquidus);
        const double slope = (ss304::Model::cp_liquid-cp_solidus)/width;
        increment = (b-a)*(cp_solidus+.5*slope*((a-solidus)+(b-solidus))+
                           ss304::Model::latent_fusion/width);
      } else {
        increment = (b-a)*ss304::Model::cp_liquid;
      }
      sum += increment;
      a = b;
    }
    const double answer = sign*model->density()*sum;
    need(finite(answer), "nonfinite v2 reference-volume enthalpy increment");
    return answer;
  };
  material.validate();
  return material;
}

} // namespace electrode_heat
