#ifndef SS304_THERMOELASTIC_V3_H
#define SS304_THERMOELASTIC_V3_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* All functions return 0 on success, -1 on failure, unless stated otherwise.
 * Failure leaves accepted state unchanged; tfv3_last_error() describes it.
 * Native arrays use canonical row-major (y,x), x fastest. Cells are nx*ny,
 * nodes are (nx+1)*(ny+1). All energies already include the full height once.
 * One handle is one disconnected electrode. CPU-only native MFEM conduction.
 */
typedef void* tfv3_handle;

typedef struct tfv3_config {
    uint32_t abi_version; /* exactly 3 */
    int32_t nx, ny;
    double x0_m, y0_m, width_m, thickness_m, height_m;
    double reference_density_kg_m3; /* explicit fixed-reference-volume density */
    double initial_temperature_K; /* scalar baseline/reference used by native solver */
    double initial_time_s;
    const char* mesh_revision;
    const char* material_sha256;
    const char* source_sha256;
    const char* electrode_id;
    /* Mandatory full nodal state; no implicit temperature/mirroring defaults. */
    const double* initial_nodal_temperature_K;
    size_t initial_node_count;
} tfv3_config;

typedef struct tfv3_state_info {
    uint64_t revision;
    double time_s;
    size_t node_count, cell_count;
    double reference_density_kg_m3;
    double stored_reference_energy_J;
    double cumulative_mechanical_exchange_J;
    double minimum_supported_temperature_K, maximum_supported_temperature_K;
    char state_sha256[65];
} tfv3_state_info;

typedef struct tfv3_epoch {
    uint64_t base_revision;
    uint64_t mechanical_sequence; /* exactly base_revision+1 */
    double start_time_s, dt_s, stage_time_s; /* stage_time_s == start_time_s+dt_s */
    const char* mesh_revision;
    const char* material_sha256;
    const char* source_sha256;
    const char* em_state_sha256; /* matching declared state; zero-Joule still explicit */
    const char* base_state_sha256;
    /* Exactly base_state_sha256 + "/mechanical/" + decimal(mechanical_sequence).
     * Discarded corrector trials may retry this identity; committed replay fails. */
    const char* transaction_id;
} tfv3_epoch;

typedef struct tfv3_budget {
    double thermal_delta_J;
    double mechanical_exchange_J; /* positive heats native caloric state; signed */
    double mechanical_absolute_activity_J;
    double joule_J, authored_heat_J, wall_net_J; /* all exactly zero in this ABI */
    double residual_J, roundoff_bound_J;
    double accepted_reference_energy_J, trial_reference_energy_J;
    double nodal_residual_max_abs_J;
    double mechanical_load_conservation_residual_J;
    int32_t newton_iterations;
} tfv3_budget;

const char* tfv3_last_error(void);
int tfv3_create(const tfv3_config* config, tfv3_handle* out_handle);
void tfv3_destroy(tfv3_handle handle);
int tfv3_get_state_info(tfv3_handle handle, tfv3_state_info* info);
/* trial=0 accepted, trial=1 prepared; getters never interpolate/extrapolate. */
int tfv3_copy_nodal_temperature(tfv3_handle handle, int trial, double* output, size_t count);
/* Discrete native cell energy is Vcell/4 * sum(h_vol(T_corner)). Summed
 * cell energies equal native lumped nodal enthalpy, not Gauss h(T_Q1). */
int tfv3_copy_cell_energy(tfv3_handle handle, int trial, double* output, size_t count);
int tfv3_copy_nodal_reference_volume(tfv3_handle handle, double* output, size_t count);
/* x and y arrays are separate canonical coordinate vectors, full node count. */
int tfv3_copy_node_coordinates(tfv3_handle handle, double* x_m, double* y_m, size_t count);
int tfv3_prepare(tfv3_handle handle, const tfv3_epoch* epoch,
                 const double* signed_mechanical_cell_J, size_t cell_count,
                 tfv3_budget* budget);
/* Dedicated per-node ledger, full-height J, canonical order:
 * delta_H - mechanical_load + conduction_out = residual (external sources zero).
 * Uniform-cell mechanical load distributes each cell's exact J to four corners.
 * Do not falsely identify reaggregated cell delta_H with its original cell debit. */
int tfv3_copy_trial_nodal_ledger(tfv3_handle handle, double* delta_H_J,
                                double* mechanical_load_J, double* conduction_out_J,
                                double* residual_J, size_t node_count);
int tfv3_copy_trial_cell_mechanical(tfv3_handle handle, double* output, size_t cell_count);
/* preflight validates the same prepared transaction and accepted base without mutation. */
int tfv3_preflight(tfv3_handle handle, const char* transaction_id, const char* base_state_sha256);
int tfv3_commit(tfv3_handle handle, const char* transaction_id, const char* base_state_sha256);
int tfv3_discard(tfv3_handle handle, const char* transaction_id);
/* Versioned checkpoint includes numerical/material/source identity, nodal T,
 * cumulative mechanical ledger and replay watermarks. Save refuses pending
 * trials; load requires fresh revision-zero handle with identical initial config.
 * Parent owns atomic multi-component checkpoint/commit orchestration. */
int tfv3_save(tfv3_handle handle, const char* path);
int tfv3_load(tfv3_handle handle, const char* path);

#ifdef __cplusplus
}
#endif
#endif
