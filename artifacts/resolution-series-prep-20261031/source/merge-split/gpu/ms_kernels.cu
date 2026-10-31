// ms_kernels.cu -- NVRTC kernels of the merge/split transaction (schema merge-split-transaction-v1).
// Compiled by the hall-gpu-resident-v1 common session with --fmad=false --prec-div=true
// --prec-sqrt=true --ftz=false for the actual device. All floating-point work is in
// ms_core.h (the same source as the CPU reference). No atomics; every sum has a fixed
// order (one thread per segment or group, sequential), and the ledger uses exact integer
// superaccumulators, so reduction order cannot change any bit.
#include "ms_core.h"
#include "deterministic_primitives.cuh"

#define MS_LEDGER_Q 10u   // quantities per species: W_in W_out P_in[3] P_out[3] E_in E_out

extern "C" __global__ void ms_classify(const MsParams* P, ms_u32 n, const ms_u64* id, const ms_u32* sp, const double* x0,
                                       const double* x1, const double* x2, const double* u0, const double* u1,
                                       const double* u2, const double* w, ms_u32* cell, ms_u32* bin, double* ke,
                                       ms_u32* key_lo, ms_u32* status) {
    const ms_u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    ms_u32 c = 0, b = 0;
    double k = 0;
    const int ok = ms_classify(P, sp[i], x0[i], x1[i], x2[i], u0[i], u1[i], u2[i], w[i], &c, &b, &k);
    int bad = !ok || id[i] == 0ull || (i > 0 && !(id[i - 1] < id[i]));
    cell[i] = c;
    bin[i] = b;
    ke[i] = k;
    key_lo[i] = ((sp[i] & 0xffu) << 24) | (b & 0xffffffu);
    status[i] = bad ? 1u : 0u;
}

extern "C" __global__ void ms_gather_u32(ms_u32 n, const ms_u32* perm, const ms_u32* src, ms_u32* dst) {
    const ms_u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) dst[i] = src[perm[i]];
}

extern "C" __global__ void ms_gather_sorted(ms_u32 n, ms_u32 n_species, const ms_u32* perm, const ms_u64* id,
                                            const ms_u32* sp, const ms_u32* fl, const double* x0, const double* x1,
                                            const double* x2, const double* u0, const double* u1, const double* u2,
                                            const double* w, const ms_u32* cell, const ms_u32* bin, const double* ke,
                                            ms_u64* sid, ms_u32* ssp, ms_u32* sfl, ms_u32* sbin, double* sx0, double* sx1,
                                            double* sx2, double* su0, double* su1, double* su2, double* sw, double* ske,
                                            ms_u32* sseg) {
    const ms_u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const ms_u32 o = perm[i];
    sid[i] = id[o];
    ssp[i] = sp[o];
    sfl[i] = fl[o];
    sbin[i] = bin[o];
    sx0[i] = x0[o];
    sx1[i] = x1[o];
    sx2[i] = x2[o];
    su0[i] = u0[o];
    su1[i] = u1[o];
    su2[i] = u2[o];
    sw[i] = w[o];
    ske[i] = ke[o];
    sseg[i] = cell[o] * n_species + sp[o];
}

extern "C" __global__ void ms_segment_flags(ms_u32 nseg, const ms_u32* begin, const ms_u32* end, ms_u32* nonempty) {
    const ms_u32 g = blockIdx.x * blockDim.x + threadIdx.x;
    if (g < nseg) nonempty[g] = end[g] > begin[g] ? 1u : 0u;
}

extern "C" __global__ void ms_form_count(const MsParams* P, ms_u32 nseg, const ms_u32* begin, const ms_u32* end, MsView v,
                                         ms_u32* counts, ms_u32* refused) {
    const ms_u32 g = blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= nseg) return;
    ms_u32 c = 0, r = 0;
    if (end[g] > begin[g]) c = ms_form_segment(P, v, begin[g], end[g], g / P->n_species, g % P->n_species, (MsGroup*)0, &r);
    counts[g] = c;
    refused[g] = r;
}

extern "C" __global__ void ms_form_write(const MsParams* P, ms_u32 nseg, const ms_u32* begin, const ms_u32* end, MsView v,
                                         const ms_u32* offset, MsGroup* groups) {
    const ms_u32 g = blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= nseg || !(end[g] > begin[g])) return;
    ms_u32 r = 0;
    (void)ms_form_segment(P, v, begin[g], end[g], g / P->n_species, g % P->n_species, groups + offset[g], &r);
}

extern "C" __global__ void ms_kinematics(const MsParams* P, ms_u32 G, MsGroup* groups, MsView v, MsMergeResult* merge,
                                         ms_u32* n_out, ms_u32* n_reloc) {
    const ms_u32 g = blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= G) return;
    MsGroup gr = groups[g];
    MsMergeResult r;
    for (int d = 0; d < 3; ++d) { r.c[d] = 0; r.ua[d] = 0; r.ub[d] = 0; }
    ms_u32 gf = 0;
    if (gr.kind == MS_KIND_MERGE) {
        gr.status = ms_merge_group(P, v, &gr, &r, &gf);
        gr.n_out = gr.status == MS_OK ? 2u : 0u;
    } else {
        gr.status = ms_split_group(P, v, &gr);
        gr.n_out = gr.status == MS_OK ? (1u << P->species[gr.species].split_levels) : 0u;
    }
    gr.gflags = gf;
    if (gr.status != MS_OK)
        for (int d = 0; d < 3; ++d) { r.c[d] = 0; r.ua[d] = 0; r.ub[d] = 0; }
    groups[g] = gr;
    merge[g] = r;
    n_out[g] = gr.n_out;
    n_reloc[g] = (gr.kind == MS_KIND_MERGE && gr.status == MS_OK) ? gr.k : 0u;
}

struct MsOutputs {
    ms_u64* id;
    ms_u32* sp;
    ms_u32* fl;
    double* x[3];
    double* u[3];
    double* w;
};
struct MsReloc {
    ms_u64* id;
    ms_u32* group;
    double* start[3];
    double* D[3];
    double* charge;
};

extern "C" __global__ void ms_emit(const MsParams* P, ms_u32 G, MsGroup* groups, const MsMergeResult* merge,
                                   const ms_u32* out_off, const ms_u32* reloc_off, MsView v, ms_u64* first_id,
                                   MsOutputs o, MsReloc rl, ms_u32* consumed) {
    const ms_u32 g = blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= G) return;
    MsGroup gr = groups[g];
    gr.out_offset = out_off[g];
    gr.reloc_offset = reloc_off[g];
    groups[g] = gr;
    first_id[g] = v.id[gr.start];
    if (gr.status != MS_OK) return;
    const MsSpecies* sp = &P->species[gr.species];
    const ms_u64 base = P->base_id + (ms_u64)gr.out_offset;
    if (gr.kind == MS_KIND_MERGE) {
        const MsMergeResult m = merge[g];
        const double Wh = gr.W * 0.5;
        for (ms_u32 j = 0; j < 2; ++j) {
            const ms_u32 k = gr.out_offset + j;
            o.id[k] = base + j;
            o.sp[k] = gr.species;
            o.fl[k] = MS_FLAG_MERGE | MS_FLAG_SPLIT;
            for (int d = 0; d < 3; ++d) {
                o.x[d][k] = m.c[d];
                o.u[d][k] = j ? m.ub[d] : m.ua[d];
            }
            o.w[k] = Wh;
        }
        double cf[3];
        for (int d = 0; d < 3; ++d) cf[d] = m.c[d] + P->shift[d];
        ms_u32 r = gr.reloc_offset;
        for (ms_u32 s = gr.start; s < gr.end; ++s) {
            if (!(v.flags[s] & MS_FLAG_MERGE)) continue;
            consumed[s] = 1u;
            rl.id[r] = v.id[s];
            rl.group[r] = g;
            for (int d = 0; d < 3; ++d) {
                const double st = v.x[d][s] + P->shift[d];
                rl.start[d][r] = st;
                rl.D[d][r] = cf[d] - st;
            }
            rl.charge[r] = v.w[s] * sp->charge_C;
            ++r;
        }
    } else {
        const ms_u32 s = gr.start;
        consumed[s] = 1u;
        const ms_u32 clones = 1u << sp->split_levels;
        const double wm = v.w[s] * (1.0 / (double)clones);
        for (ms_u32 j = 0; j < clones; ++j) {
            const ms_u32 k = gr.out_offset + j;
            o.id[k] = base + j;
            o.sp[k] = gr.species;
            o.fl[k] = MS_FLAG_MERGE | MS_FLAG_SPLIT;
            for (int d = 0; d < 3; ++d) {
                o.x[d][k] = v.x[d][s];
                o.u[d][k] = v.u[d][s];
            }
            o.w[k] = wm;
        }
    }
}

extern "C" __global__ void ms_scatter_u32(ms_u32 n, const ms_u32* perm, const ms_u32* src, ms_u32* dst) {
    const ms_u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) dst[perm[i]] = src[i];
}

extern "C" __global__ void ms_compact_removed(ms_u32 n, const ms_u32* flag, const ms_u32* offset, const ms_u64* id, ms_u64* removed) {
    const ms_u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n && flag[i]) removed[offset[i]] = id[i];
}

// Exact ledger, stage 1: thread (chunk c, quantity q) accumulates its chunk of the
// consumed inputs (sorted order) or of the outputs into one superaccumulator.
extern "C" __global__ void ms_ledger_chunks(const MsParams* P, ms_u32 n, MsView v, const ms_u32* consumed, ms_u32 nout,
                                            const ms_u32* osp, const double* ow, const double* ou0, const double* ou1,
                                            const double* ou2, ms_u32 chunks, ms_i64* partial) {
    const ms_u32 c = blockIdx.x * blockDim.x + threadIdx.x, q = blockIdx.y;
    if (c >= chunks) return;
    const ms_u32 s = q / MS_LEDGER_Q, k = q % MS_LEDGER_Q;
    const double m = P->species[s].mass_kg;
    const int out = k == 1u || (k >= 5u && k <= 7u) || k == 9u;
    const ms_u32 total = out ? nout : n;
    const ms_u32 per = (total + chunks - 1u) / chunks;
    const ms_u32 lo = c * per;
    const ms_u32 hi = lo + per < total ? lo + per : total;
    MsAcc acc;
    ms_acc_zero(&acc);
    for (ms_u32 i = lo; i < hi; ++i) {
        if (out) {
            if (osp[i] != s) continue;
            const double w = ow[i];
            if (k == 1u) ms_acc_add(&acc, w);
            else if (k <= 7u) ms_acc_add_product3(&acc, w, k == 5u ? ou0[i] : (k == 6u ? ou1[i] : ou2[i]), m);
            else ms_acc_add_product2(&acc, w, ms_kinetic_energy(m, ou0[i], ou1[i], ou2[i]));
        } else {
            if (!consumed[i] || v.species[i] != s) continue;
            const double w = v.w[i];
            if (k == 0u) ms_acc_add(&acc, w);
            else if (k <= 4u) ms_acc_add_product3(&acc, w, v.u[k - 2u][i], m);
            else ms_acc_add_product2(&acc, w, v.ke[i]);
        }
    }
    ms_acc_normalize(&acc);
    ms_i64* dst = partial + ((ms_u64)q * chunks + c) * MS_ACC_DIGITS;
    for (int d = 0; d < MS_ACC_DIGITS; ++d) dst[d] = acc.d[d];
}

// Exact ledger, stage 2: digit (q, d) = sum over chunks (integer, exact).
extern "C" __global__ void ms_ledger_final(ms_u32 nq, ms_u32 chunks, const ms_i64* partial, ms_i64* total) {
    const ms_u32 t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= nq * MS_ACC_DIGITS) return;
    const ms_u32 q = t / MS_ACC_DIGITS, d = t % MS_ACC_DIGITS;
    ms_i64 sum = 0;
    for (ms_u32 c = 0; c < chunks; ++c) sum += partial[((ms_u64)q * chunks + c) * MS_ACC_DIGITS + d];
    total[t] = sum;
}
