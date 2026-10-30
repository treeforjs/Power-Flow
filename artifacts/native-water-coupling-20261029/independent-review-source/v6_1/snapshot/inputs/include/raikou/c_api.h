#ifndef RAIKOU_C_API_H
#define RAIKOU_C_API_H

/** Raikou's versioned C ABI. All quantities use SI units except Te_eV. */
#define RAIKOU_VERSION "1.0.0"
#define RAIKOU_VERSION_MAJOR 1
#define RAIKOU_VERSION_MINOR 0
#define RAIKOU_VERSION_PATCH 0
#define RAIKOU_NUM_LEVELS 5

#if defined(_WIN32) && !defined(RAIKOU_STATIC_DEFINE)
#  if defined(RAIKOU_BUILDING_LIBRARY)
#    define RAIKOU_API __declspec(dllexport)
#  else
#    define RAIKOU_API __declspec(dllimport)
#  endif
#else
#  define RAIKOU_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef enum raikou_status {
    RAIKOU_OK = 0,
    RAIKOU_INVALID_ARGUMENT = 1,
    RAIKOU_NUMERICAL_ERROR = 2,
    RAIKOU_INTERNAL_ERROR = 3
} raikou_status;

typedef struct raikou_plasma_state {
    double Te_eV;
    double ne_m3;
    double Tgas_K;
    double nH_ground_m3;
    double nH2_m3;
    double path_length_m;
} raikou_plasma_state;

typedef struct raikou_radiance_result {
    double populations_m3[RAIKOU_NUM_LEVELS];
    double lyman_alpha_radiance_W_m2_sr;
    double lyman_beta_radiance_W_m2_sr;
    double balmer_alpha_radiance_W_m2_sr;
    double balmer_beta_radiance_W_m2_sr;
    double molecular_lyman_proxy_W_m2_sr;
    double tau_0_lyman_alpha;
    double escape_factor_lyman_alpha;
    double tau_0_lyman_beta;
    double escape_factor_lyman_beta;
} raikou_radiance_result;

/** Matrices index [source level - 1][destination level - 1].
 * Diagonal and inapplicable transitions are zero; vectors index level - 1. */
typedef struct raikou_rate_result {
    double excitation_m3_s[RAIKOU_NUM_LEVELS][RAIKOU_NUM_LEVELS];
    double deexcitation_m3_s[RAIKOU_NUM_LEVELS][RAIKOU_NUM_LEVELS];
    double ionization_m3_s[RAIKOU_NUM_LEVELS];
    double radiative_recombination_m3_s[RAIKOU_NUM_LEVELS];
    double three_body_recombination_m6_s[RAIKOU_NUM_LEVELS];
    double h2_dissociation_m3_s[RAIKOU_NUM_LEVELS];
} raikou_rate_result;

typedef struct raikou_fit_priors {
    double min_Te_eV;
    double max_Te_eV;
    double min_ne_m3;
    double max_ne_m3;
    double min_nH_m3;
    double max_nH_m3;
    double min_nH2_m3;
    double max_nH2_m3;
} raikou_fit_priors;

typedef struct raikou_fit_result {
    double Te_eV;
    double Te_err_eV;
    double ne_m3;
    double ne_err_m3;
    double nH_m3;
    double nH_err_m3;
    double nH2_m3;
    double nH2_err_m3;
    double tau_0_lyman_alpha;
    double escape_factor_lyman_alpha;
    double chi2_reduced;
    int converged;
    int iterations;
} raikou_fit_result;

typedef struct raikou_spark_config {
    double peak_current;
    double pulse_duration_s;
    double initial_density;
    double plasma_conductivity;
    double channel_length;
} raikou_spark_config;

typedef struct raikou_spark_state {
    double time_s;
    double current_A;
    double channel_radius_m;
    double pressure_Pa;
    double temperature_K;
    double deposited_energy_J;
    double accumulated_action;
} raikou_spark_state;

RAIKOU_API const char* raikou_version(void);
/** JSON runtime entrypoint. Allocated output is released with raikou_free_string.
 * Response is {ok:true,result:...} or {ok:false,error:...}; status is zero on success.
 * The embedding host owns MPI initialization; see backend.h for communicator selection. */
RAIKOU_API int raikou_execute_json(const char* operation, const char* request, char** output);
RAIKOU_API void raikou_free_string(char* output);
/** Error description local to the calling thread, valid until its next API call.
 * Successful operations clear the error. No C++ exception crosses this ABI.
 * Output pointers are required and remain unchanged on failure. */
RAIKOU_API const char* raikou_last_error(void);
RAIKOU_API raikou_status raikou_default_plasma_state(raikou_plasma_state* out);
RAIKOU_API raikou_status raikou_evaluate(const raikou_plasma_state* state, raikou_radiance_result* out);
RAIKOU_API raikou_status raikou_rate_coefficients(double Te_eV, raikou_rate_result* out);
RAIKOU_API raikou_status raikou_default_fit_priors(raikou_fit_priors* out);
/** A null priors pointer selects the default fit priors. A successful call can
 * return converged == 0; callers should check this diagnostic before using a fit. */
RAIKOU_API raikou_status raikou_fit_vuv_shot(
    double measured_ly_alpha_radiance, double measured_ly_alpha_sigma,
    double measured_molecular_radiance, double measured_molecular_sigma,
    double path_length_m, double fixed_Tgas_K,
    const raikou_fit_priors* priors, raikou_fit_result* out);
RAIKOU_API raikou_status raikou_default_spark_config(raikou_spark_config* out);
RAIKOU_API raikou_status raikou_default_spark_state(raikou_spark_state* out);
/** current and out may point to the same state, permitting in-place stepping. */
RAIKOU_API raikou_status raikou_step_spark(
    const raikou_spark_state* current, double dt_s,
    const raikou_spark_config* config, raikou_spark_state* out);

#ifdef __cplusplus
}
#endif
#endif
