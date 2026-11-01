#pragma once
#include "common_balance.hpp"
namespace powerflow::common_material_stocks_v1 {
namespace fs=::powerflow::common_conservation_v1;namespace h=::hall_corrective_v3;
struct EnergyBound {
 double hall_field_heat_J=0,motion_work_J=0,source_partition_J=0,reduction_roundoff_J=0,compensated_field_roundoff_J=0,total_J=0;
 std::uint64_t reduction_operations=0;
 double reduction_operand_norm_J=0;
};
inline double upper_bound(long double x){need(std::isfinite(x)&&x>=0&&x<=std::numeric_limits<double>::max(),"finite derived energy bound");if(x==0)return 0;const double d=static_cast<double>(x);return static_cast<long double>(d)<x?std::nextafter(d,INFINITY):d;}
// Reuses every original component guard, then bounds the change of arithmetic
// from their ordered rounded reductions to independent exact stock sums. The
// norm is formed from actual transition operands, not reported residuals and
// not the large absolute thermal reference stock. Signed physical conversion
// accounts remain separate; they are not added as discretionary tolerances.
inline EnergyBound derive_energy_bound(const h::Owner&owner,const h::JointProposal&proposal,
 std::span<const ss::s7::Trial>source_trials,const Inventory&native_before,const Inventory&native_after,
 const BoundaryStock&boundary,const fs::FieldThermalAssessment&field){
 need(proposal.configuration_sha256()==owner.configuration_sha256()&&proposal.endpoint_revision()==owner.revision()+1,"energy bound candidate identity");
 const auto&r=proposal.field_heat().field_receipt();const auto&thermal=proposal.field_heat().heat();const auto&motion=proposal.motion();const auto&mesh=owner.mesh();const double dt=motion.binding().dt_s;
 EnergyBound out;
 // Exactly the FieldHeatProposal guard, including its load activity terms.
 double transfers=std::abs(r.field_delta_J)+std::abs(r.plasma_work_J)+std::abs(thermal.thermal_delta_J)+std::abs(r.outer_loss_J)+std::abs(r.outer_supply_J)+std::abs(r.circuit_work_J);
 if(thermal.loads)transfers+=std::abs(thermal.wall_heat_J)+thermal.external_heating_J+thermal.external_cooling_J;
 out.hall_field_heat_J=2e-20+1e-10*transfers+thermal.roundoff_bound_J;
 out.motion_work_J=h::motion_work_tolerance(motion.work_scales());
 const long double eps=std::numeric_limits<double>::epsilon();long double norm=0,source_bound=0;std::uint64_t operands=0;
 auto add=[&](long double v){need(std::isfinite(v),"finite arithmetic norm operand");norm+=std::abs(v);need(++operands<=100000000,"bounded arithmetic norm");};
 const auto&before=owner.field();const auto&after=proposal.field_heat().candidate_field();const auto&E=proposal.trial_E();const auto&initial=owner.model().initial_fields();const double I=owner.model().drive().current(motion.binding().begin_s+.5*dt);
 for(int i=0;i<mesh.ne;++i){const h::Words a(before.E[i],before.E_low[i]),b(after.E[i],after.E_low[i]);add(mesh.M[i]*((a+b)*.5).wide()*(b-a).wide());add(static_cast<long double>(dt)*I*initial.feed[i]*E[i]);add(static_cast<long double>(dt)*initial.background_force[i]*E[i]);add(static_cast<long double>(dt)*mesh.Y[i]*E[i]*E[i]);}
 for(int i=0;i<mesh.nb;++i){const h::Words a(before.B[i],before.B_low[i]),b(after.B[i],after.B_low[i]);add(mesh.H[i]*((a+b)*.5).wide()*(b-a).wide());}
 const auto law=fs::eh::conditional_ss304_v2(ss304v2::Qualification::evaluated_solid_1500,ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,8000.,fs::eh::material_headers_sha256);
 for(int k=mesh.kshort;k<mesh.kfeed;++k)for(int body=0;body<2;++body){const auto&accepted=owner.thermal().slice(body,k);const auto&candidate=thermal.bodies.at(owner.thermal().index(body,k)).state();const auto areas=fs::lumped_control_areas(accepted.geometry());need(areas.size()==candidate.T.size(),"actual thermal norm shape");for(std::size_t j=0;j<areas.size();++j)add(static_cast<long double>(areas[j])*law.delta_h(candidate.T[j],accepted.state().T[j])*accepted.geometry().height);}
 // Actual endpoint kinetic stocks and sink operands bound weighted-KE sum
 // conversion. Their independent exact inventories are already reconstructed.
 for(const auto*x:{&native_before.kinetic_J,&native_after.kinetic_J,&boundary.foil_kinetic_J,&boundary.short_kinetic_J,&boundary.exterior_kinetic_J})add(x->value());
 for(const auto&t:source_trials){
  need(t.source_candidates().empty()&&t.source_handoffs().empty(),"source layer bound needs qualified held-layer operands");
  for(const auto&p:t.packets()){add(p.kinetic_J);add(p.internal_J);add(p.adsorbed_reference_change_J);}
  for(const auto&b:t.surface_births()){const auto&l=b.phase().packet.particle.ledger;long double activity=0;for(double v:{l.wall_heat_J,l.outgoing_K_J,l.incoming_K_J,l.electron_reference_change_J,l.image_reference_change_J,l.resolved_matching_field_work_J}){add(v);activity+=std::abs(static_cast<long double>(v));}source_bound+=static_cast<long double>(1e-10+128*std::numeric_limits<double>::epsilon())*activity;}
  for(const auto&load:t.loads())add(h::transfer_energy_J(load));
 }
 // At most four arithmetic operations per scalar terminal/thermal/source term,
 // plus weighted kinetic reduction. 32 covers the fixed final composition.
 // Exact accumulators introduce no reduction rounding themselves.
 need(motion.paths().size()<=1000000,"bounded actual motion norm count");operands+=8*motion.paths().size();out.reduction_operations=4*operands+32;const long double neps=out.reduction_operations*eps;need(neps<.01L,"bounded floating reduction length");out.reduction_operand_norm_J=upper_bound(norm);out.reduction_roundoff_J=upper_bound(neps/(1-neps)*norm);
 // Two-word endpoint add/subtract/midpoint error is second order in epsilon.
 // Unlike thermal reference energy, this small compensation error depends on
 // EM endpoint magnitudes even when their net transition cancels.
 const long double em_norm=std::abs(static_cast<long double>(field.field_before.total().value()))+std::abs(static_cast<long double>(field.field_after.total().value()));out.compensated_field_roundoff_J=upper_bound(64*eps*eps*em_norm);
 out.source_partition_J=upper_bound(source_bound);out.total_J=upper_bound(static_cast<long double>(out.hall_field_heat_J)+out.motion_work_J+out.source_partition_J+out.reduction_roundoff_J+out.compensated_field_roundoff_J);return out;
}
inline EnergyBound derive_energy_bound(const h::Owner&owner,const h::JointProposal&proposal,const ss::PreparedSourceStage&source,const Inventory&before,const Inventory&after,const BoundaryStock&boundary,const fs::FieldThermalAssessment&field){return derive_energy_bound(owner,proposal,std::span<const ss::s7::Trial>(source.source_trials()),before,after,boundary,field);}
inline void require_energy_closure(const CombinedAssessment&assessment,const EnergyBound&bound){
 for(double x:{bound.hall_field_heat_J,bound.motion_work_J,bound.source_partition_J,bound.reduction_roundoff_J,bound.compensated_field_roundoff_J,bound.total_J})need(std::isfinite(x)&&x>=0,"invalid derived energy tolerance");
 const double expected=upper_bound(static_cast<long double>(bound.hall_field_heat_J)+bound.motion_work_J+bound.source_partition_J+bound.reduction_roundoff_J+bound.compensated_field_roundoff_J);need(bound.total_J==expected,"derived energy bound altered");
 g::Exact allowed;allowed.add_product(bound.total_J,1);const auto r=assessment.accounted_residual_J;auto positive=r,negative=r;negative.merge(allowed);positive.merge(allowed,1);
 // value() can round tiny exact digits to zero. Compare signed canonical top
 // digits of (r-bound) and (r+bound), never only the rounded residual.
 const auto hi=positive.normalised(),lo=negative.normalised();need(hi.back()<0||positive.zero(),"whole-system energy exceeds positive derived bound");need(lo.back()>=0,"whole-system energy exceeds negative derived bound");
}
} // namespace powerflow::common_material_stocks_v1
