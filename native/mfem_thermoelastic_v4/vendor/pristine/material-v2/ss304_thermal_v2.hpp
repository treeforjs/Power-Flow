#pragma once
#include "sealed1800.hpp"
#include <cstdint>

namespace ss304v2 {
using Values = ss304::Values;
enum class Qualification { evaluated_solid_1500=0, nominal_1800=1, source_conditioned_liquid_3000=2 };
struct EnergySupport { double minimum_j,maximum_j,cold_margin_j,hot_margin_j; int range_status; };
class Model {
  Qualification q_;
  ss304::Model base_;
  static double cold_h(double t) {
    return (t-300.0)*(472.0+.068*(t+300.0)-2820000.0/(300.0*t));
  }
 public:
  static constexpr double minimum_temperature_k=298.2;
  static constexpr double kim_liquid_cp=0.190*4184.0;
  Model(Qualification q, ss304::DensityConvention convention, double rho_ref)
      : q_(q), base_(q==Qualification::evaluated_solid_1500 ? ss304::Qualification::evaluated_solid_1500 : ss304::Qualification::nominal_synthesis_1800, convention,rho_ref) {
    if (q!=Qualification::evaluated_solid_1500 && q!=Qualification::nominal_1800 && q!=Qualification::source_conditioned_liquid_3000)
      throw std::invalid_argument("unknown qualification");
    if (!std::isfinite(rho_ref*max_u())) throw std::invalid_argument("density exceeds representable caloric domain");
  }
  double max_temperature_k() const { return q_==Qualification::evaluated_solid_1500 ? 1500 : q_==Qualification::nominal_1800 ? 1800 : 3000; }
  double density() const { return base_.reference_density_kg_m3(); }
  double min_u() const { return cold_h(minimum_temperature_k); }
  double max_u() const {
    if (q_!=Qualification::source_conditioned_liquid_3000) return base_.max_enthalpy_j_kg();
    return base_.max_enthalpy_j_kg()+kim_liquid_cp*1200.0;
  }
  Values evaluate(double t) const {
    if (!std::isfinite(t) || t<minimum_temperature_k || t>max_temperature_k())
      throw std::out_of_range("temperature outside selected reviewed domain");
    if (t>=300 && t<=std::min(1800.0,max_temperature_k())) return base_.evaluate(t);
    Values v{};v.temperature_k=t;
    if (t<300) {
      v.h_j_kg=cold_h(t);v.cp_sensible_j_kg_k=472+.136*t-2820000/(t*t);
      // Source Table2 brackets298.2K with293K/300K; no extrapolation.
      const double f=(t-293)/7;
      v.resistivity_ohm_m=(71.3+f*(71.9-71.3))*1e-8;
      v.k_w_m_k=14.76+f*(14.89-14.76);
    } else {
      const auto a=base_.evaluate(1800);
      v.h_j_kg=a.h_j_kg+kim_liquid_cp*(t-1800);
      v.cp_sensible_j_kg_k=kim_liquid_cp;
      // Source-conditioned shape: Kim Eq29, W/(cm K)->W/(m K),
      // normalized to the admitted Ho-Chu1800K value. Explicit model choice.
      v.k_w_m_k=a.k_w_m_k*(.1229+.00003248*t)/(.1229+.00003248*1800);
      const double lorenz=a.k_w_m_k*a.resistivity_ohm_m/1800;
      v.resistivity_ohm_m=lorenz*t/v.k_w_m_k;
      v.liquid_fraction=1;v.phase=2;v.synthesized_high_temperature=1;
    }
    v.u_model_j_kg=v.h_j_kg;v.apparent_heat_capacity_j_kg_k=v.cp_sensible_j_kg_k;
    v.sigma_s_m=1/v.resistivity_ohm_m;
    v.energy_j_m3_reference=density()*v.u_model_j_kg;
    v.capacity_j_m3_reference_k=density()*v.apparent_heat_capacity_j_kg_k;
    return v;
  }
  double temperature_from_u(double u) const {
    if (!std::isfinite(u) || u<min_u() || u>max_u()) throw std::out_of_range("specific energy outside reviewed domain");
    if (u==min_u()) return minimum_temperature_k;
    if (u==max_u()) return max_temperature_k();
    if (u>=0 && u<=base_.max_enthalpy_j_kg()) return base_.temperature_from_u_model_k(u);
    double lo=minimum_temperature_k,hi=max_temperature_k();
    for (int i=0;i<128;++i) {
      const double t=lo+(hi-lo)/2;
      if (t==lo || t==hi) break;
      const double h=evaluate(t).u_model_j_kg;
      if (h==u) return t;
      if (h<u) lo=t;else hi=t;
    }
    return std::abs(evaluate(lo).u_model_j_kg-u)<=std::abs(evaluate(hi).u_model_j_kg-u) ? lo : hi;
  }
  Values evaluate_inventory(double t,double mass,double reference_volume) const {
    if (!std::isfinite(mass)||mass<=0||!std::isfinite(reference_volume)||reference_volume<=0)
      throw std::invalid_argument("invalid current condensed mass/reference volume");
    Values v=evaluate(t);
    const double mass_per_reference_volume=mass/reference_volume;
    v.energy_j_m3_reference=mass_per_reference_volume*v.u_model_j_kg;
    v.capacity_j_m3_reference_k=mass_per_reference_volume*v.apparent_heat_capacity_j_kg_k;
    if (!std::isfinite(v.energy_j_m3_reference)||!std::isfinite(v.capacity_j_m3_reference_k)||v.capacity_j_m3_reference_k<=0)
      throw std::out_of_range("inventory caloric values not representable");
    // k and sigma remain INTRINSIC condensed-phase properties. Geometry/mixture
    // owner must supply any effective-medium/remap closure; no implicit scaling.
    return v;
  }
  double temperature_from_reference_energy(double e) const {
    if (!std::isfinite(e) || e<density()*min_u() || e>density()*max_u()) throw std::out_of_range("reference-volume energy outside reviewed domain");
    if (e==density()*min_u()) return minimum_temperature_k;
    if (e==density()*max_u()) return max_temperature_k();
    return temperature_from_u(e/density());
  }
  EnergySupport assess_energy(double mass,double total_energy) const {
    if (!std::isfinite(mass)||mass<=0||!std::isfinite(total_energy)) throw std::invalid_argument("invalid energy-budget input");
    const double lo=mass*min_u(),hi=mass*max_u();
    if (!std::isfinite(lo)||!std::isfinite(hi)) throw std::out_of_range("energy budget overflow");
    return {lo,hi,total_energy-lo,hi-total_energy,total_energy<lo ? -1 : total_energy>hi ? 1 : 0};
  }
};

// Lumping/congruent evaporation and a neutral monatomic vapor are explicit
// scenario approximations. No bulk liquid conductivity is assigned to vapor.
struct VaporValues {
  double saturation_pressure_pa, kim_pressure_pa;
  double vapor_h_j_kg, vapor_u_j_kg, latent_h_j_kg;
  double gas_cp_j_kg_k, gas_cv_j_kg_k, atomic_reference_u0_j_kg;
};
class NeutralVapor {
  Model material_;
 public:
  static constexpr double ru=8.31446261815324;
  static constexpr double avogadro=6.02214076e23;
  static constexpr double kb=1.380649e-23;
  static constexpr double molar_mass=.713*.055845+.202*.0519961+.085*.0586934;
  static constexpr double r=ru/molar_mass;
  static constexpr double cp=2.5*r;
  static constexpr double latent_star=1770.0*4184.0;
  static constexpr double tstar=3000.0;
  explicit NeutralVapor(const Model& material):material_(material) {
    if (material.max_temperature_k()!=3000) throw std::invalid_argument("vapor requires explicit3000K source-conditioned scenario");
  }
  static double kim_pressure(double t) {
    if (!std::isfinite(t) || t<1800 || t>3000) throw std::out_of_range("vapor source/model domain1800-3000K");
    return 101325.0*std::pow(10.0,6.1210-18836.0/t);
  }
  double atomic_u0() const { return material_.evaluate(tstar).h_j_kg+latent_star-cp*tstar; }
  double material_reference_density() const { return material_.density(); }
  VaporValues evaluate(double t) const {
    const double raw=kim_pressure(t);
    const double dc=cp-Model::kim_liquid_cp;
    const double latent=latent_star+dc*(t-tstar);
    const double log_ratio=dc/r*std::log(t/tstar)-(latent/t-latent_star/tstar)/r;
    const double p=kim_pressure(tstar)*std::exp(log_ratio);
    const double h=material_.evaluate(tstar).h_j_kg+latent_star+cp*(t-tstar);
    return {p,raw,h,h-r*t,latent,cp,1.5*r,atomic_u0()};
  }
  // Net Hertz-Knudsen flux, same surface/background temperature assumption.
  // Signed return; negative flux is condensation, never clipped to evaporation.
  double net_mass_flux(double t,double back_pressure_pa,double accommodation) const {
    const auto v=evaluate(t);
    if (!std::isfinite(back_pressure_pa)||back_pressure_pa<0||!std::isfinite(accommodation)||accommodation<0||accommodation>1)
      throw std::invalid_argument("invalid vapor pressure/accommodation");
    return accommodation*(v.saturation_pressure_pa-back_pressure_pa)/std::sqrt(2*std::acos(-1.0)*r*t);
  }
  // A declared neutral hard-sphere alternative, with caller-supplied diameter.
  // No fitted collision diameter or ionized/plasma conductivity is invented.
  double neutral_gas_k(double t,double diameter_m,double ionized_fraction) const {
    (void)evaluate(t);
    if (!std::isfinite(diameter_m)||diameter_m<=0||ionized_fraction!=0)
      throw std::invalid_argument("requires positive supplied collision diameter and strictly neutral gas; plasma owner supplies charged transport");
    const double m=molar_mass/avogadro;
    const double k=75*kb/(64*diameter_m*diameter_m)*std::sqrt(kb*t/(std::acos(-1.0)*m));
    if (!std::isfinite(k)||k<=0) throw std::out_of_range("gas conductivity outside representable range");
    return k;
  }
};

struct MetalState {
  double mass_kg, energy_j, reference_volume_m3;
  std::uint64_t revision,last_event;
};
struct TransferTrial {
  MetalState next;
  double temperature_k, emitted_reference_plus_kinetic_j;
  double signed_energy_residual_j,caloric_projection_residual_j;
  double remaining_mass_fraction,removed_reference_volume_m3;
};
// Canonical accepted metal inventory belongs to heat owner. This pure prepare
// function does not mutate state, emit packets, advance RNG or commit revisions.
inline TransferTrial prepare_transfer(const Model& material,const NeutralVapor& vapor,
    const MetalState& base,std::uint64_t expected_revision,std::uint64_t event,
    double outward_mass_kg,double outward_actual_kinetic_j,double other_heat_j) {
  if (material.max_temperature_k()!=3000 || material.density()!=vapor.material_reference_density())
    throw std::invalid_argument("material/vapor scenario mismatch");
  if (base.revision!=expected_revision||event<=base.last_event||base.revision==std::numeric_limits<std::uint64_t>::max())
    throw std::invalid_argument("stale/replayed/overflowed metal transaction");
  if (!std::isfinite(base.mass_kg)||base.mass_kg<=0||!std::isfinite(base.energy_j)||
      !std::isfinite(base.reference_volume_m3)||base.reference_volume_m3<=0||
      !std::isfinite(outward_mass_kg)||!std::isfinite(outward_actual_kinetic_j)||!std::isfinite(other_heat_j)||
      (outward_mass_kg>0&&outward_actual_kinetic_j<0)||(outward_mass_kg<0&&outward_actual_kinetic_j>0)||
      (outward_mass_kg==0&&outward_actual_kinetic_j!=0))
    throw std::invalid_argument("invalid signed mass/energy/reference volume");
  const double old_t=material.temperature_from_u(base.energy_j/base.mass_kg);
  if (outward_mass_kg>0&&old_t<1800) throw std::out_of_range("evaporation source model starts1800K");
  const double mass=base.mass_kg-outward_mass_kg;
  if (!std::isfinite(mass)||mass<=0) throw std::out_of_range("depletion needs empty-cell/geometry handoff; no mass clamp");
  const double packet=std::fma(outward_mass_kg,vapor.atomic_u0(),outward_actual_kinetic_j);
  const double energy=base.energy_j+other_heat_j-packet;
  const double t=material.temperature_from_u(energy/mass);
  const double projected=mass*material.evaluate(t).u_model_j_kg;
  const long double residual=static_cast<long double>(base.energy_j)+other_heat_j-energy-packet;
  const double fraction=mass/(material.density()*base.reference_volume_m3);
  if (!std::isfinite(fraction)) throw std::out_of_range("reference mass fraction overflow");
  return {{mass,energy,base.reference_volume_m3,base.revision+1,event},t,packet,
          static_cast<double>(residual),projected-energy,fraction,outward_mass_kg/material.density()};
}
} // namespace ss304v2
