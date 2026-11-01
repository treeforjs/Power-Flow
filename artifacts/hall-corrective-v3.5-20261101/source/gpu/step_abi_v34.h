// step_abi_v34.h -- plain-data layout shared by the host engine (device_step_v34.cpp) and the NVRTC
// kernels (step_kernels_v34.cu) of the Hall v3.4 device coupled step. Device addresses are 64-bit
// integers; no CUDA type appears.
#ifndef HALL_V34_STEP_ABI_H
#define HALL_V34_STEP_ABI_H

typedef unsigned long long s34_u64;
typedef unsigned int s34_u32;

#define S34_TPB 256
#define S34_NONE 0xffffffffu

// Exact-sum slots of one motion fold (index into the XsOut array).
#define S34_SUM_K 0      // kinetic_change_J = RN(sum w_i d_i), all census carriers
#define S34_SUM_W 1      // field_work_J = RN(dt RN(sum E_e J_e))
#define S34_SUM_GK 2     // G_K = RN(sum_{q != 0} w |d|)
#define S34_SUM_CK 3     // C_K = RN(sum_{q != 0} w K0 + w RN(K0 + d))
#define S34_SUM_Q 4      // Q = RN(sum_{q != 0} |q| w)
#define S34_SUM_GW 5     // G_W = RN(dt RN(sum |E| |J|))
#define S34_SUMS 6

// Operands of the exact sums (kind selects the term).
struct S34SumArgs {
    s34_u64 n;
    s34_u64 w, q, dke, k0;       // carriers (census order)
    s34_u64 a, b;                // field: E, J (or r, r / b, b for the residual norms)
};

// Results of one device motion iterate (downloaded once per iterate).
struct S34Results {
    s34_u32 status[8];           // push: census rank/code, motion rank/code, charge rank/code, edge records, 0
    double plasma_dt, max_gyro, max_dcells;
    s34_u32 field_nonfinite;     // first nonfinite E/B index (S34_NONE: none)
    s34_u32 current_nonfinite;   // first nonfinite J index
    double continuity_C, gross_max, field_bound;
    s34_u32 free_guard_fail;     // first node failing free_charge_guard (endpoint free charge)
    s34_u32 assemble_fail;       // first carrier failing the record operand checks
    s34_u32 support_fail;        // first carrier failing the deposit support (rank); code below
    s34_u32 support_code;
    s34_u32 support_records;
    s34_u32 survivors, sinks, births;
    s34_u32 field_nonfinite_b;   // first nonfinite B index
    s34_u32 pad;
};
// Preflight recomputation results.
struct S34Preflight {
    s34_u32 endpoint_charge_fail;   // first endpoint failing charge_of_endpoints
    s34_u32 differ[3];              // first node where the recomputed endpoint free charge (hi, lo, gross) differs from the candidate slot
    s34_u32 pad[4];
};

// Census admission results (begin of a proposal).
struct S34Census {
    s34_u32 charge_fail;         // first resident failing charge_of_particles (rank); code below
    s34_u32 charge_code;
    s34_u32 guard_fail;          // first node failing free_charge_guard of the census charge
    s34_u32 same_fail;           // first node where the census charge differs from the accepted free charge
    s34_u32 census_index_lo, census_index_hi, census_code, pad0;   // census v2 check (first failing record, code)
    s34_u32 neutral_fail;        // first resident neutral carrier (native v7 neutral interval bound)
    s34_u32 newborn_fail_rank;   // smallest census index of a newborn failing the flight/start rule; code below
    s34_u32 newborn_code;
    s34_u32 newborn_absent;      // number of newborn IDs absent from the census
};

#endif
