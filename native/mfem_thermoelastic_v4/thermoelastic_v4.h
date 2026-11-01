#ifndef SS304_THERMOELASTIC_V4_H
#define SS304_THERMOELASTIC_V4_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Single persistent native transaction participant. All energies are full-height
 * weighted J. Native nodes/cells are canonical (y,x), x fastest. All source laws
 * are exact copied sealed-v2 calls; source action weights are caller events, not
 * predicted event counts. Return0 success,-1 rejection; see last_error. */
typedef void* tfv4_handle;
typedef struct tfv4_patch {
 const char *id,*gas_state_id,*parent_id,*scenario_sha256;
 int32_t face; /* LowY0,HighY1,LowX2,HighX3 */
 double begin,end;
 int32_t water_enabled,water_range; /* reviewed0:300..1300, conditional1:298.2..3000 */
 double initial_coverage,registry_mass_kg;
} tfv4_patch;
typedef struct tfv4_config {
 uint32_t abi_version; /*4*/
 int32_t nx,ny;
 double x0_m,y0_m,width_m,thickness_m,height_m,reference_density_kg_m3;
 double initial_temperature_K,initial_time_s;
 const char *mesh_revision,*material_sha256,*source_sha256,*electrode_id;
 const double* initial_nodal_temperature_K;size_t initial_node_count;
 const tfv4_patch* patches;size_t patch_count;
 const char *electron_state_id,*electron_parent_id;
} tfv4_config;
typedef struct tfv4_state_info {
 uint64_t revision;double time_s;size_t node_count,cell_count;
 double reference_density_kg_m3,stored_reference_energy_J,cumulative_mechanical_exchange_J;
 double minimum_supported_temperature_K,maximum_supported_temperature_K;
 char state_sha256[65]; /* SHA256 complete accepted heat+source+gas/path payload */
 size_t patch_count,layer_count,gas_packet_count,path_count;
} tfv4_state_info;
typedef struct tfv4_epoch {
 uint64_t base_revision,mechanical_sequence;double start_time_s,dt_s,stage_time_s;
 const char *mesh_revision,*material_sha256,*source_sha256,*em_state_sha256;
 const char *base_state_sha256,*transaction_id;
 /* Explicit source accepted-start field/thermal binding; positive IDs required
  * whenever source actions exist. field_state must equal common_state. */
 uint64_t common_state,field_state,source_mesh_revision,source_material_revision;
 double field_time_s,E_V_m[3],B_T[3];
 int32_t material_phase_admitted,planar_barrier_admitted,electron_supply_admitted;
} tfv4_epoch;
typedef struct tfv4_source_action {
 int32_t kind; /*1 begin,2 advance,3 collect existing PICHandedOff electron*/
 uint64_t layer_id;size_t patch_index;uint64_t patch_id;
 double weight,normal_uniform,transverse_uniform,azimuth_uniform,origin_m[3];
 double duration_s,work_function_eV;int32_t energy_panels,action_nodes;double tail_widths;
 double local_linear_vacuum_extent_m,image_force_fraction,image_energy_fraction;
 double maximum_quantum_magnetic_energy_fraction,maximum_gyroangle,maximum_relative_normal_step;
 int32_t maximum_substeps,PIC_applies_microscopic_image_force,own_return_unit_accommodation;
 /* Collection uses this ACTUAL weighted K, matching the same handed-off native
  * particle K in this no-PIC-transport adapter; arbitrary intervening PIC work
  * requires a separate admitted particle owner, not an invented energy change. */
 double incoming_weighted_K_J;int32_t collection_accommodation_admitted;
} tfv4_source_action;
typedef struct tfv4_water_draw {size_t patch_index;double normal,radius,angle,along,z,time;} tfv4_water_draw;
typedef struct tfv4_budget {
 /* v3 prefix retained; wall_net is now effective external thermal input. */
 double thermal_delta_J,mechanical_exchange_J,mechanical_absolute_activity_J;
 double joule_J,authored_heat_J,wall_net_J,residual_J,roundoff_bound_J;
 double accepted_reference_energy_J,trial_reference_energy_J,nodal_residual_max_abs_J;
 double mechanical_load_conservation_residual_J;int32_t newton_iterations;
 double source_wall_J,source_particle_K_change_J,source_electron_reference_change_J;
 double source_image_energy_change_J,source_field_work_J;
 double neutral_K_J,neutral_internal_J,adsorbed_reference_change_J,gas_reference_J;
 double particle_charge_change_C,electrode_charge_change_C,physical_electron_number_change;
 size_t new_gas_packet_count,new_path_count;
} tfv4_budget;
/* Rows preserve native accepted source/gas state, not a second population. */
typedef struct tfv4_layer_info {
 uint64_t layer_id;size_t patch_index;int32_t status; /*Active0,PICHandedOff1,Recollected2*/
 double weight,time_s,origin_m[3],position_m[3],proper_velocity_m_s[3];
 double current_weighted_K_J,work_function_eV,outward_normal[3];
} tfv4_layer_info;
typedef struct tfv4_path_info {uint64_t layer_id;double from_m[3],to_m[3],start_time_s,end_time_s,charge_C,field_work_J;} tfv4_path_info;
typedef struct tfv4_gas_info {
 size_t patch_index;uint64_t source_revision;double weight,T_K,time_s;
 double position_m[3],proper_velocity_m_s[3],wall_support_impulse_N_s[3];
 double mass_kg,charge_C,kinetic_J,internal_J,reference_J,adsorbed_reference_change_J;
} tfv4_gas_info;
const char* tfv4_last_error(void);
int tfv4_create(const tfv4_config*,tfv4_handle*);void tfv4_destroy(tfv4_handle);
int tfv4_get_state_info(tfv4_handle,int trial,tfv4_state_info*);
int tfv4_copy_nodal_temperature(tfv4_handle,int trial,double*,size_t);
int tfv4_copy_cell_energy(tfv4_handle,int trial,double*,size_t);
int tfv4_copy_nodal_reference_volume(tfv4_handle,double*,size_t);
int tfv4_copy_node_coordinates(tfv4_handle,double*,double*,size_t);
int tfv4_prepare(tfv4_handle,const tfv4_epoch*,const double* cell_J,size_t,
 const tfv4_source_action*,size_t,const tfv4_water_draw*,size_t,tfv4_budget*);
/* Seven arrays: deltaH,mechanical,conduction,residual,external,electron,water.
 * deltaH - mechanical - external + conduction = residual;
 * external=electron+water, water=-(gasK+internal+dUads), all gas reference0. */
int tfv4_copy_trial_nodal_ledger(tfv4_handle,double*,double*,double*,double*,double*,double*,double*,size_t);
int tfv4_copy_trial_cell_mechanical(tfv4_handle,double*,size_t);
int tfv4_copy_layers(tfv4_handle,int trial,tfv4_layer_info*,size_t);
int tfv4_copy_paths(tfv4_handle,int trial,tfv4_path_info*,size_t);
int tfv4_copy_gas_packets(tfv4_handle,int trial,tfv4_gas_info*,size_t);
/* Per patch rows: accepted/candidate adsorbed Nhi,Nlo,gas Nhi,Nlo,
 * negative adsorbed reference diagnostic (stable per-step debit is budget).
 * Gas totals and19 sealed thermal cumulative terms are checkpointed. */
int tfv4_copy_reservoirs(tfv4_handle,int trial,double*,size_t);
/*19 sealed cumulative terms + mechanical + particle/electrode charge +electronN*/
int tfv4_copy_cumulative(tfv4_handle,int trial,double*,size_t);
int tfv4_preflight(tfv4_handle,const char*,const char*);
int tfv4_commit(tfv4_handle,const char*,const char*);
int tfv4_discard(tfv4_handle,const char*);
int tfv4_save(tfv4_handle,const char*);int tfv4_load(tfv4_handle,const char*);
#ifdef __cplusplus
}
#endif
#endif
