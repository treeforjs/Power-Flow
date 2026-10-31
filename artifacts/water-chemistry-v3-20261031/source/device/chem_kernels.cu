// chem_kernels.cu -- NVRTC program of chemistry v3. Compiled by the shared GPU
// layer (hall-gpu-resident-v1-20261030/common) with --fmad=false
// --prec-div=true --prec-sqrt=true --ftz=false and the architecture of the
// actual device, from the same header bytes the host compiles. One thread per
// carrier; every output slot is owned by exactly one thread (no atomics, no
// reductions): carrier i writes outs[i] and events[i*max_events ...].
#include "chem_event.h"

extern "C" __global__ void pfc_mcc_step(pfc::TableView t, pfc::CoverView cv, pfc::StepConfig cfg, const pfc::CarrierIn* carriers,
                                        const pfc::CellTargets* cells, unsigned n_cells, unsigned n, unsigned first,
                                        pfc::EventRec* events, pfc::CarrierOut* outs) {
    const unsigned i = first + blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const pfc::CarrierIn c = carriers[i];
    if (c.cell >= n_cells) {
        pfc::CarrierOut o;
        o.id = c.id;
        o.status = PFC_ST_BAD_INPUT;
        o.n_events = 0;
        o.n_candidates = 0;
        o.species_final = c.species;
        o.u_final[0] = c.u[0];
        o.u_final[1] = c.u[1];
        o.u_final[2] = c.u[2];
        o.nu_dt_max = 0.0;
        o.n_cover = 0;
        o.pad = 0;
        for (int k = 0; k < PFC_MAX_COVER_SLOTS; ++k) {
            o.cover_region[k] = 0xFFFFFFFFu;
            o.cover_exposure[k] = 0.0;
            o.cover_bound[k] = 0.0;
        }
        outs[i] = o;
        return;
    }
    pfc::mcc_carrier(t, cv, cfg, c, cells[c.cell], events + (unsigned long long)i * cfg.max_events, outs + i);
}

// Elementary-function parity: out[4i..4i+3] = log(x), log1p(y), cos(2 pi u), sin(2 pi u).
extern "C" __global__ void pfc_math_parity(const double* x, const double* y, const double* u, unsigned n, double* out) {
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[4ull * i] = pfc::log(x[i]);
    out[4ull * i + 1] = pfc::log1p(y[i]);
    double c, s;
    pfc::sincos_2pi(u[i], &c, &s);
    out[4ull * i + 2] = c;
    out[4ull * i + 3] = s;
}
