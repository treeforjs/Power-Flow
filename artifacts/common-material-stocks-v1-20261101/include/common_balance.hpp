#pragma once
#include "material_stocks.hpp"
#include "field_thermal_stocks.hpp"
namespace powerflow::common_material_stocks_v1 {
struct CombinedAssessment {
 Inventory material_before,material_after;
 g::Exact physical_stored_before_J,physical_stored_after_J,material_boundary_J,external_net_J;
 g::Exact raw_residual_J,accounted_residual_J,hydrogen_residual,oxygen_residual,charge_residual_e;
 g::Exact event_numerical_defect_J,event_reference_conversion_J,thermal_reference_rounding_J;
 g::Exact source_kinematic_rounding_J,source_binding_operand_conversion_J;
};
// Composition of independently recovered actual operands. This helper does not
// claim the supplied values are authoritative: root acquires/pins their owners,
// verifies the exact stage/common bindings, and performs all live preflights.
// No supplied component residual enters either result. Boundary increments must
// be retained in the common checkpoint so successive transitions telescope.
inline CombinedAssessment combine(const Inventory&native_before,const Inventory&native_after,
 const SourceAssessment&source,const GasStock&gas_before,const GasStock&gas_after,
 const BoundaryStock&boundary,const ChemistryOperands&chem,
 const ::powerflow::common_conservation_v1::FieldThermalAssessment&field){
 need(diff(gas_after.numerical_defect_J,gas_before.numerical_defect_J)==chem.numerical_defect_J,"actual event defect differs from gas ledger");
 need(diff(gas_after.process_Q_J,gas_before.process_Q_J)==chem.process_Q_J,"actual process Q differs from gas ledger");
 CombinedAssessment out;out.material_before=native_before;out.material_before.add(source.before);out.material_after=native_after;out.material_after.add(source.after);
 out.physical_stored_before_J=out.material_before.energy();out.physical_stored_before_J.merge(gas_before.energy());out.physical_stored_before_J.merge(field.field_before.total());out.physical_stored_before_J.merge(field.thermal_before_J);
 out.physical_stored_after_J=out.material_after.energy();out.physical_stored_after_J.merge(gas_after.energy());out.physical_stored_after_J.merge(field.field_after.total());out.physical_stored_after_J.merge(field.thermal_after_J);
 out.material_boundary_J=boundary.conductor.energy();out.material_boundary_J.merge(boundary.exterior.energy());out.material_boundary_J.merge(boundary.short_kinetic_J);out.material_boundary_J.merge(boundary.exterior_kinetic_J);out.material_boundary_J.merge(boundary.exterior_correction_J);
 out.external_net_J=field.circuit_input_J;out.external_net_J.merge(field.outer_input_J);out.external_net_J.merge(field.outer_loss_J,1);
 out.raw_residual_J=diff(out.physical_stored_after_J,out.physical_stored_before_J);out.raw_residual_J.merge(source.emitted_electrode_reference_J);out.raw_residual_J.merge(out.material_boundary_J);out.raw_residual_J.merge(out.external_net_J,1);
 out.hydrogen_residual=diff(out.material_after.hydrogen,out.material_before.hydrogen);out.hydrogen_residual.merge(boundary.conductor.hydrogen);out.hydrogen_residual.merge(boundary.exterior.hydrogen);
 out.oxygen_residual=diff(out.material_after.oxygen,out.material_before.oxygen);out.oxygen_residual.merge(boundary.conductor.oxygen);out.oxygen_residual.merge(boundary.exterior.oxygen);
 out.charge_residual_e=diff(out.material_after.charge_e,out.material_before.charge_e);out.charge_residual_e.merge(source.emitted_electrode_charge_e);out.charge_residual_e.merge(boundary.conductor.charge_e);out.charge_residual_e.merge(boundary.exterior.charge_e);
 need(out.hydrogen_residual.zero()&&out.oxygen_residual.zero()&&out.charge_residual_e.zero(),"whole material nuclei/charge imbalance");
 auto actual_reference=diff(out.material_after.reference_J,out.material_before.reference_J);actual_reference.merge(boundary.conductor.reference_J);actual_reference.merge(boundary.exterior.reference_J);need(actual_reference==chem.species_reference_delta_J,"whole chemical-reference stock differs from actual admitted events");
 out.event_numerical_defect_J=chem.numerical_defect_J;out.event_reference_conversion_J=chem.reference_conversion_J;out.thermal_reference_rounding_J=field.thermal_reference_rounding_J;out.source_binding_operand_conversion_J=source.binding_operand_conversion_J;
 out.source_kinematic_rounding_J=source.neutral_kinematic_rounding_J;out.source_kinematic_rounding_J.merge(source.electron_kinematic_rounding_J);
 out.accounted_residual_J=out.raw_residual_J;out.accounted_residual_J.merge(out.event_numerical_defect_J);out.accounted_residual_J.merge(out.event_reference_conversion_J,1);out.accounted_residual_J.merge(out.thermal_reference_rounding_J,1);out.accounted_residual_J.merge(out.source_kinematic_rounding_J);out.accounted_residual_J.merge(out.source_binding_operand_conversion_J,1);
 return out;
}
} // namespace powerflow::common_material_stocks_v1
