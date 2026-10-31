// phase_kernels.h -- the per-cell, per-EM-cell, per-element and per-patch arithmetic of one owner step.
// The CPU reference calls these functions in index loops; the NVRTC program wraps the same functions in
// extern "C" kernels (one thread per item). No atomics: every output word is written by exactly one item, and
// every sum runs over a fixed index order (energy and atom sums are exact 128-bit integer sums).
//
// Layout. Cells c = (s*ny + j)*nx + i (slice s, thermal row j through the thickness, column i across the width).
// Faces f: 0 LowY (-y), 1 HighY (+y), 2 LowX (-x), 3 HighX (+x) (the electrode_heat_v2::Face order).
// Surface elements e: exposed (cell, face) pairs that carry loads, deposits, releases, evaporation and
// prescribed-temperature exchange; sorted by patch, then cell, then face; every element belongs to one patch.
#pragma once
#include "phase_law.h"

struct pt_dims {
    int S, nx, ny, nx_em, ny_em, rx;
    int n_cells, n_em, n_el, n_patches;
    double dx, dz;
    int n_sub, sub;               // sub-step count of this step and the current sub-step index
    double tau;                   // sub-step length (s)
    double dt;                    // whole step (s), for quotes
    int joule_mode;               // 0: uniform power density per EM cell (Hall v3), 1: sigma-weighted (parallel axial conductors)
    int evap_mode;                // 0: none, 1: expected HKL flux evaluated per sub-step, 2: quantized releases supplied by the caller
    int bc_kind[4];               // per face kind: 0 adiabatic/loads only, 1 prescribed temperature (benchmarks)
    double bc_T[4];
    int evap_face[4];             // per face kind: evaporating
    double phi_floor;             // admitted minimum mass fraction of the reference mass during a step
    double phi_merge;             // end-of-step merge threshold (owner policy)
    int n_scatter;                // number of surface elements with loads/deposits this step
};

struct pt_buffers {
    const pt_law* law;
    // geometry
    const double* row_dy;         // [ny]
    const int* row_em;            // [ny] EM row of each thermal row
    const int* em_row_begin;      // [ny_em+1] first thermal row of each EM row
    const double* gx;             // [ny] x-face conductance factor dy_j*dz/dx
    const double* gy;             // [ny-1] y-face factor dx*dz/((dy_j+dy_j+1)/2)
    const double* vref;           // [n_cells] reference volume
    const double* mref;           // [n_cells] reference mass
    // state (current and next slot)
    const pt_i128* E;  const pt_i128* N;  const int* flags;      // E [n_cells], N [3 n_cells], flags bit0 = void
    pt_i128* E_out;    pt_i128* N_out;
    // per-cell scratch
    double* T; double* k; double* sig; double* m; int* status;
    pt_i128* joule_total;         // [n_cells] Joule energy of the step
    const double* em_energy;      // [n_em] Hall Joule energy of the step (J)
    // surface elements
    const int* el_of;             // [4 n_cells] element of (cell, face) or -1
    const int* el_cell; const int* el_face; const double* el_area; const int* el_patch;
    const int* patch_begin;       // [n_patches+1] element range of each patch
    const pt_i128* in_load_E;     // [n_el] step totals
    const pt_i128* in_dep_N; const pt_i128* in_dep_E;
    pt_i128* in_rel_N; pt_i128* in_rel_E;     // written by the release split (quantized mode)
    pt_i128* acc_ev_N; pt_i128* acc_ev_Eref; pt_i128* acc_ev_Ekin; pt_i128* acc_bnd_E;   // [n_el*3], [n_el] x3
    double* el_quote;             // [n_el*3] expected atoms over dt at the accepted state
    // patch outputs
    const pt_i128* rel_total_N;   // [n_patches*3] released atoms (quantized mode), units 2^-32 atom
    const double* rel_kin_J;      // [n_patches*3] actual kinetic energy of those atoms (J)
    double* patch_out;            // [n_patches*PT_PATCH_OUT]
    pt_i128* patch_sum;           // [n_patches*PT_PATCH_SUM]
    pt_i128* slice_sum;           // [S*PT_SLICE_SUM]
    int* slice_status;            // [S*4]: first failing cell, code, merge-candidate count, first candidate
    double* slice_maxT;           // [S]
    double* sigma_em;             // [n_em] effective axial conductivity of each EM cell (S/m)
    const int* sc_el;             // [n_scatter] element of each scatter entry (distinct)
    const pt_i128* sc_val;        // [n_scatter*5]: load_E, dep_N x3, dep_E
    pt_i128* in_load_E_w; pt_i128* in_dep_N_w; pt_i128* in_dep_E_w;   // writable aliases for the scatter
};

#define PT_PATCH_OUT 12   // area, mean T, max T, mean liquid fraction, mean melt depth, D_Fe, D_Cr, D_Ni (atoms over dt), j_Fe, j_Cr, j_Ni (atoms/m2/s), mean surface density
#define PT_PATCH_SUM 9    // acc_ev_N x3, acc_ev_Eref, acc_ev_Ekin, acc_bnd_E, load_E, dep_E, rel_E   (sums over the patch elements)
#define PT_SLICE_SUM 4    // E, N_Fe, N_Cr, N_Ni

PT_HD int pt_neighbor(const pt_dims* d, int c, int f) {
    const int i = c % d->nx, j = (c / d->nx) % d->ny;
    if (f == 0) return j > 0 ? c - d->nx : -1;
    if (f == 1) return j + 1 < d->ny ? c + d->nx : -1;
    if (f == 2) return i > 0 ? c - 1 : -1;
    return i + 1 < d->nx ? c + 1 : -1;
}
// Geometric factor (area/distance) of the interior face between c and its neighbour across f.
PT_HD double pt_face_g(const pt_dims* d, const pt_buffers* b, int c, int f) {
    const int j = (c / d->nx) % d->ny;
    if (f == 0) return b->gy[j - 1];
    if (f == 1) return b->gy[j];
    return b->gx[j];
}
// Boundary face area and centre distance of cell c on face f.
PT_HD double pt_face_area(const pt_dims* d, const pt_buffers* b, int c, int f) {
    const int j = (c / d->nx) % d->ny;
    return f < 2 ? d->dx * d->dz : b->row_dy[j] * d->dz;
}
PT_HD double pt_face_half(const pt_dims* d, const pt_buffers* b, int c, int f) {
    const int j = (c / d->nx) % d->ny;
    return f < 2 ? 0.5 * b->row_dy[j] : 0.5 * d->dx;
}
// Energy (J) conducted from cell a to cell b over tau; always called with a < b so that both cells use the
// same bits.
PT_HD double pt_face_energy(double Ta, double ka, double Tb, double kb, double g, double tau) {
    const double kh = (ka == 0.0 || kb == 0.0) ? 0.0 : 2.0 * ka * kb / (ka + kb);
    return kh * g * (Ta - Tb) * tau;
}
PT_HD pt_i128 pt_part(pt_i128 total, int n, int s, int* ov) {
    pt_u32 r = 0;
    pt_i128 q = pt_divmod_u32(total, (pt_u32)n, &r, ov);
    if ((pt_u32)s < r) { pt_i128 one; one.lo = 1; one.hi = 0; q = pt_add(q, one, ov); }
    return q;
}
PT_HD void pt_set_status(int* st, int code) { if (*st == PT_OK) *st = code; }

// Atoms (double) and mass (kg) of a cell.
PT_HD double pt_cell_mass(const pt_law* L, const pt_i128* N, int c, double* atoms) {
    double m = 0.0;
    for (int el = 0; el < PT_NEL; ++el) { atoms[el] = pt_to_double(N[3 * c + el], PT_N_INV); m += atoms[el] * L->atom_mass[el]; }
    return m;
}

// 1. Thermodynamic state of one cell from (E, N).
PT_HD void pt_thermo(const pt_dims* d, const pt_buffers* b, int c) {
    (void)d;
    int st = PT_OK;
    if (b->flags[c] & 1) { b->T[c] = 0; b->k[c] = 0; b->sig[c] = 0; b->m[c] = 0; return; }
    double atoms[PT_NEL];
    const double m = pt_cell_mass(b->law, b->N, c, atoms);
    const double h = pt_to_double(b->E[c], PT_E_INV) / m;
    const double T = pt_T_of_h(b->law, h, &st);
    b->T[c] = T; b->m[c] = m;
    b->k[c] = pt_k(b->law, T);
    b->sig[c] = 1.0 / pt_rhoe(b->law, T);
    if (b->status[c] == PT_OK) b->status[c] = st;   // the first error of the step is kept
}

// 1b. Scatter of this step's surface loads and deposits (one entry per distinct element).
PT_HD void pt_scatter(const pt_dims* d, const pt_buffers* b, int q) {
    (void)d;
    const int e = b->sc_el[q];
    b->in_load_E_w[e] = b->sc_val[5 * q + 0];
    for (int el = 0; el < PT_NEL; ++el) b->in_dep_N_w[3 * e + el] = b->sc_val[5 * q + 1 + el];
    b->in_dep_E_w[e] = b->sc_val[5 * q + 4];
}

// 1c. Effective axial conductivity of one EM cell: sub-cells are parallel conductors for an axial field,
//     sigma_eff = sum(sigma V) / sum(V) over its non-void sub-cells (fixed order).
PT_HD void pt_sigma_em(const pt_dims* d, const pt_buffers* b, int e) {
    const int per_slice = d->nx_em * d->ny_em;
    const int s = e / per_slice, je = (e % per_slice) / d->nx_em, ie = e % d->nx_em;
    double sv = 0.0, v = 0.0;
    for (int j = b->em_row_begin[je]; j < b->em_row_begin[je + 1]; ++j) for (int i = ie * d->rx; i < (ie + 1) * d->rx; ++i) {
        const int c = (s * d->ny + j) * d->nx + i;
        if (b->flags[c] & 1) continue;
        sv += b->sig[c] * b->vref[c]; v += b->vref[c];
    }
    b->sigma_em[e] = v > 0.0 ? sv / v : 0.0;
}

// 2. Exact partition of one EM cell's Joule energy over its thermal cells.
PT_HD void pt_joule_split(const pt_dims* d, const pt_buffers* b, int e) {
    const int per_slice = d->nx_em * d->ny_em;
    const int s = e / per_slice, je = (e % per_slice) / d->nx_em, ie = e % d->nx_em;
    const int j0 = b->em_row_begin[je], j1 = b->em_row_begin[je + 1], i0 = ie * d->rx, i1 = i0 + d->rx;
    const double energy = b->em_energy[e];
    int bad = 0, ov = 0;
    double W = 0.0, wmax = -1.0; int cmax = -1;
    for (int j = j0; j < j1; ++j) for (int i = i0; i < i1; ++i) {
        const int c = (s * d->ny + j) * d->nx + i;
        b->joule_total[c] = pt_zero();
        if (b->flags[c] & 1) continue;
        const double w = d->joule_mode == 1 ? b->sig[c] * b->vref[c] : b->vref[c];
        W += w;
        if (w > wmax) { wmax = w; cmax = c; }
    }
    if (energy == 0.0) return;
    if (cmax < 0 || !(W > 0.0)) { b->status[(s * d->ny + j0) * d->nx + i0] = PT_ERR_NONFINITE; return; }
    const pt_i128 total = pt_from_double(energy, PT_E_POW2, &bad);
    pt_i128 used = pt_zero();
    for (int j = j0; j < j1; ++j) for (int i = i0; i < i1; ++i) {
        const int c = (s * d->ny + j) * d->nx + i;
        if ((b->flags[c] & 1) || c == cmax) continue;
        const double w = d->joule_mode == 1 ? b->sig[c] * b->vref[c] : b->vref[c];
        const pt_i128 q = pt_from_double(energy * (w / W), PT_E_POW2, &bad);
        b->joule_total[c] = q;
        used = pt_add(used, q, &ov);
    }
    const pt_i128 rest = pt_sub(total, used, &ov);
    b->joule_total[cmax] = rest;
    if (bad || ov || pt_is_neg(rest)) pt_set_status(&b->status[cmax], PT_ERR_FIXED_RANGE);
}

// 3. One explicit sub-step of one cell: conduction, prescribed-temperature exchange, Joule, surface loads,
//    deposits, releases and evaporation.
PT_HD void pt_step_cell(const pt_dims* d, const pt_buffers* b, int c) {
    const pt_law* L = b->law;
    for (int el = 0; el < PT_NEL; ++el) b->N_out[3 * c + el] = b->N[3 * c + el];
    b->E_out[c] = b->E[c];
    if (b->flags[c] & 1) return;
    int st = b->status[c], ov = 0, bad = 0;
    pt_i128 dE = pt_zero();
    pt_i128 dN[PT_NEL] = {pt_zero(), pt_zero(), pt_zero()};
    const double Tc = b->T[c], kc = b->k[c];
    double atoms[PT_NEL];
    double atoms_tot = 0.0;
    for (int el = 0; el < PT_NEL; ++el) { atoms[el] = pt_to_double(b->N[3 * c + el], PT_N_INV); atoms_tot += atoms[el]; }
    for (int f = 0; f < 4; ++f) {
        const int nb = pt_neighbor(d, c, f);
        if (nb >= 0 && !(b->flags[nb] & 1)) {
            const double g = pt_face_g(d, b, c, f);
            // lower index first: f 0/2 -> neighbour is lower; f 1/3 -> c is lower
            if (f == 0 || f == 2) {
                const double F = pt_face_energy(b->T[nb], b->k[nb], Tc, kc, g, d->tau);
                dE = pt_add(dE, pt_from_double(F, PT_E_POW2, &bad), &ov);
            } else {
                const double F = pt_face_energy(Tc, kc, b->T[nb], b->k[nb], g, d->tau);
                dE = pt_sub(dE, pt_from_double(F, PT_E_POW2, &bad), &ov);
            }
            continue;
        }
        const int e = b->el_of[4 * c + f];
        if (e < 0) continue;   // side wall toward a merged (void) cell: adiabatic
        if (nb < 0 && d->bc_kind[f] == 1) {
            const double F = kc * pt_face_area(d, b, c, f) / pt_face_half(d, b, c, f) * (d->bc_T[f] - Tc) * d->tau;
            const pt_i128 q = pt_from_double(F, PT_E_POW2, &bad);
            dE = pt_add(dE, q, &ov);
            b->acc_bnd_E[e] = pt_add(b->acc_bnd_E[e], q, &ov);
        }
        dE = pt_add(dE, pt_part(b->in_load_E[e], d->n_sub, d->sub, &ov), &ov);
        dE = pt_add(dE, pt_part(b->in_dep_E[e], d->n_sub, d->sub, &ov), &ov);
        for (int el = 0; el < PT_NEL; ++el) dN[el] = pt_add(dN[el], pt_part(b->in_dep_N[3 * e + el], d->n_sub, d->sub, &ov), &ov);
        if (d->evap_mode == 2) {
            dE = pt_add(dE, pt_part(b->in_rel_E[e], d->n_sub, d->sub, &ov), &ov);
            for (int el = 0; el < PT_NEL; ++el) dN[el] = pt_add(dN[el], pt_part(b->in_rel_N[3 * e + el], d->n_sub, d->sub, &ov), &ov);
        } else if (d->evap_mode == 1 && d->evap_face[f] && L->vapor_enabled) {
            const double A = b->el_area[e];
            for (int el = 0; el < PT_NEL; ++el) {
                const double D = pt_flux(L, el, Tc, atoms[el] / atoms_tot, &st) * A * d->tau;
                const pt_i128 qn = pt_from_double(D, PT_N_POW2, &bad);
                if (pt_is_zero(qn)) continue;
                const double a = pt_to_double(qn, PT_N_INV);
                const pt_i128 qr = pt_from_double(a * L->eps0[el], PT_E_POW2, &bad);
                const pt_i128 qk = pt_from_double(a * (2.0 * PT_KB * Tc), PT_E_POW2, &bad);
                dN[el] = pt_sub(dN[el], qn, &ov);
                dE = pt_sub(dE, pt_add(qr, qk, &ov), &ov);
                b->acc_ev_N[3 * e + el] = pt_add(b->acc_ev_N[3 * e + el], qn, &ov);
                b->acc_ev_Eref[e] = pt_add(b->acc_ev_Eref[e], qr, &ov);
                b->acc_ev_Ekin[e] = pt_add(b->acc_ev_Ekin[e], qk, &ov);
            }
        }
    }
    dE = pt_add(dE, pt_part(b->joule_total[c], d->n_sub, d->sub, &ov), &ov);
    const pt_i128 En = pt_add(b->E[c], dE, &ov);
    double m_new = 0.0;
    for (int el = 0; el < PT_NEL; ++el) {
        const pt_i128 Nn = pt_add(b->N[3 * c + el], dN[el], &ov);
        if (pt_is_neg(Nn)) pt_set_status(&st, PT_ERR_NEGATIVE_INV);
        b->N_out[3 * c + el] = Nn;
        m_new += pt_to_double(Nn, PT_N_INV) * L->atom_mass[el];
    }
    if (!(m_new >= d->phi_floor * b->mref[c])) pt_set_status(&st, PT_ERR_DEPLETED);
    if (bad || ov) pt_set_status(&st, PT_ERR_FIXED_RANGE);
    b->E_out[c] = En;
    b->status[c] = st;
}

// 4. Per-slice exact sums of the state and the first failing cell; merge candidates (end of step).
PT_HD void pt_slice_reduce(const pt_dims* d, const pt_buffers* b, int s) {
    int ov = 0;
    pt_i128 sE = pt_zero(), sN[PT_NEL] = {pt_zero(), pt_zero(), pt_zero()};
    int first = -1, code = 0, cand = 0, first_cand = -1;
    double maxT = 0.0;
    const int per = d->nx * d->ny;
    for (int c = s * per; c < (s + 1) * per; ++c) {
        if (first < 0 && b->status[c] != PT_OK) { first = c; code = b->status[c]; }
        if (b->flags[c] & 1) continue;
        if (b->T[c] > maxT) maxT = b->T[c];
        sE = pt_add(sE, b->E[c], &ov);
        for (int el = 0; el < PT_NEL; ++el) sN[el] = pt_add(sN[el], b->N[3 * c + el], &ov);
        bool exposed = false;
        for (int f = 0; f < 4; ++f) if (b->el_of[4 * c + f] >= 0) exposed = true;
        if (exposed) {
            double atoms[PT_NEL];
            const double m = pt_cell_mass(b->law, b->N, c, atoms);
            if (m < d->phi_merge * b->mref[c]) { ++cand; if (first_cand < 0) first_cand = c; }
        }
    }
    pt_i128* o = b->slice_sum + PT_SLICE_SUM * s;
    o[0] = sE; for (int el = 0; el < PT_NEL; ++el) o[1 + el] = sN[el];
    b->slice_status[4 * s + 0] = ov ? (first >= 0 ? first : s * per) : first;
    b->slice_status[4 * s + 1] = ov ? PT_ERR_FIXED_RANGE : code;
    b->slice_status[4 * s + 2] = cand;
    b->slice_status[4 * s + 3] = first_cand;
    b->slice_maxT[s] = maxT;
}

// 5. Expected atoms per element over dt at the current state (quantized-release quotes).
PT_HD void pt_quote_element(const pt_dims* d, const pt_buffers* b, int e) {
    const int c = b->el_cell[e], f = b->el_face[e];
    int st = 0;
    double atoms[PT_NEL], tot = 0.0;
    for (int el = 0; el < PT_NEL; ++el) { atoms[el] = pt_to_double(b->N[3 * c + el], PT_N_INV); tot += atoms[el]; }
    for (int el = 0; el < PT_NEL; ++el) {
        double D = 0.0;
        if (d->evap_face[f] && b->law->vapor_enabled) D = pt_flux(b->law, el, b->T[c], atoms[el] / tot, &st) * b->el_area[e] * d->dt;
        b->el_quote[3 * e + el] = D;
    }
}

// 6. Split the quantized release of one patch over its elements in proportion to the element quotes
//    (remainder to the first element with the largest quote): atoms exactly, reference energy per atom
//    eps0, actual kinetic energy split by the same weights.
PT_HD void pt_release_split(const pt_dims* d, const pt_buffers* b, int p) {
    (void)d;
    const int e0 = b->patch_begin[p], e1 = b->patch_begin[p + 1];
    int bad = 0, ov = 0;
    for (int e = e0; e < e1; ++e) { b->in_rel_E[e] = pt_zero(); for (int el = 0; el < PT_NEL; ++el) b->in_rel_N[3 * e + el] = pt_zero(); }
    for (int el = 0; el < PT_NEL; ++el) {
        const pt_i128 R = b->rel_total_N[3 * p + el];
        const double K = b->rel_kin_J[3 * p + el];
        if (pt_is_zero(R) && K == 0.0) continue;
        double W = 0.0, wmax = -1.0; int emax = -1;
        for (int e = e0; e < e1; ++e) { const double w = b->el_quote[3 * e + el]; W += w; if (w > wmax) { wmax = w; emax = e; } }
        if (!(W > 0.0)) { b->status[b->el_cell[e0]] = PT_ERR_NEGATIVE_INV; continue; }
        const double Rd = pt_to_double(R, PT_N_INV);
        const pt_i128 Kq = pt_from_double(K, PT_E_POW2, &bad);
        pt_i128 usedN = pt_zero(), usedK = pt_zero();
        for (int e = e0; e < e1; ++e) {
            if (e == emax) continue;
            const double share = b->el_quote[3 * e + el] / W;
            const pt_i128 qn = pt_from_double(Rd * share, PT_N_POW2, &bad);
            const pt_i128 qk = pt_from_double(K * share, PT_E_POW2, &bad);
            usedN = pt_add(usedN, qn, &ov); usedK = pt_add(usedK, qk, &ov);
            const pt_i128 qr = pt_from_double(pt_to_double(qn, PT_N_INV) * b->law->eps0[el], PT_E_POW2, &bad);
            b->in_rel_N[3 * e + el] = pt_neg(qn, &ov);
            b->in_rel_E[e] = pt_sub(b->in_rel_E[e], pt_add(qr, qk, &ov), &ov);
        }
        const pt_i128 qn = pt_sub(R, usedN, &ov), qk = pt_sub(Kq, usedK, &ov);
        const pt_i128 qr = pt_from_double(pt_to_double(qn, PT_N_INV) * b->law->eps0[el], PT_E_POW2, &bad);
        b->in_rel_N[3 * emax + el] = pt_neg(qn, &ov);
        b->in_rel_E[emax] = pt_sub(b->in_rel_E[emax], pt_add(qr, qk, &ov), &ov);
        if (pt_is_neg(qn) || pt_is_neg(qk)) bad = 1;
    }
    if (bad || ov) b->status[b->el_cell[e0]] = PT_ERR_FIXED_RANGE;
}

// 7. Patch report at the current state and exact sums of the step accumulators of its elements.
PT_HD void pt_patch_report(const pt_dims* d, const pt_buffers* b, int p) {
    const int e0 = b->patch_begin[p], e1 = b->patch_begin[p + 1];
    double A = 0, sT = 0, mT = 0, sf = 0, sdepth = 0, sD[PT_NEL] = {0, 0, 0}, sj[PT_NEL] = {0, 0, 0}, srho = 0;
    int ov = 0, st = 0;
    pt_i128 acc[PT_PATCH_SUM];
    for (int q = 0; q < PT_PATCH_SUM; ++q) acc[q] = pt_zero();
    for (int e = e0; e < e1; ++e) {
        const int c = b->el_cell[e], f = b->el_face[e];
        const double a = b->el_area[e], T = b->T[c];
        A += a; sT += a * T; if (T > mT) mT = T;
        sf += a * pt_fl(b->law, T);
        if (b->law->kind == 0) srho += a * pt_ss_density(T);
        // contiguous melt layer from this face inward along the normal
        double depth = 0.0;
        int cc = c;
        const int inward = f == 0 ? 1 : f == 1 ? 0 : f == 2 ? 3 : 2;
        while (cc >= 0 && !(b->flags[cc] & 1)) {
            const double fl = pt_fl(b->law, b->T[cc]);
            if (fl <= 0.0) break;
            depth += fl * (f < 2 ? b->row_dy[(cc / d->nx) % d->ny] : d->dx);
            cc = pt_neighbor(d, cc, inward);
        }
        sdepth += a * depth;
        double atoms[PT_NEL], tot = 0;
        for (int el = 0; el < PT_NEL; ++el) { atoms[el] = pt_to_double(b->N[3 * c + el], PT_N_INV); tot += atoms[el]; }
        for (int el = 0; el < PT_NEL; ++el) {
            const double j = (d->evap_face[f] && b->law->vapor_enabled) ? pt_flux(b->law, el, T, atoms[el] / tot, &st) : 0.0;
            sj[el] += a * j;
            sD[el] += j * a * d->dt;
        }
        for (int el = 0; el < PT_NEL; ++el) acc[el] = pt_add(acc[el], b->acc_ev_N[3 * e + el], &ov);
        acc[3] = pt_add(acc[3], b->acc_ev_Eref[e], &ov);
        acc[4] = pt_add(acc[4], b->acc_ev_Ekin[e], &ov);
        acc[5] = pt_add(acc[5], b->acc_bnd_E[e], &ov);
        acc[6] = pt_add(acc[6], b->in_load_E[e], &ov);
        acc[7] = pt_add(acc[7], b->in_dep_E[e], &ov);
        acc[8] = pt_add(acc[8], b->in_rel_E[e], &ov);
    }
    double* o = b->patch_out + PT_PATCH_OUT * p;
    o[0] = A; o[1] = A > 0 ? sT / A : 0; o[2] = mT; o[3] = A > 0 ? sf / A : 0; o[4] = A > 0 ? sdepth / A : 0;
    for (int el = 0; el < PT_NEL; ++el) { o[5 + el] = sD[el]; o[8 + el] = A > 0 ? sj[el] / A : 0; }
    o[11] = A > 0 ? srho / A : 0;
    for (int q = 0; q < PT_PATCH_SUM; ++q) b->patch_sum[PT_PATCH_SUM * p + q] = acc[q];
    (void)ov;
}

// Note on dN for deposits: in_dep_N is per element and per chemical element (3 per surface element).

#if defined(__CUDACC_RTC__)
extern "C" __global__ void k_thermo(pt_dims d, pt_buffers b) { const int c = blockIdx.x * blockDim.x + threadIdx.x; if (c < d.n_cells) pt_thermo(&d, &b, c); }
extern "C" __global__ void k_scatter(pt_dims d, pt_buffers b) { const int q = blockIdx.x * blockDim.x + threadIdx.x; if (q < d.n_scatter) pt_scatter(&d, &b, q); }
extern "C" __global__ void k_sigma(pt_dims d, pt_buffers b) { const int e = blockIdx.x * blockDim.x + threadIdx.x; if (e < d.n_em) pt_sigma_em(&d, &b, e); }
extern "C" __global__ void k_joule(pt_dims d, pt_buffers b) { const int e = blockIdx.x * blockDim.x + threadIdx.x; if (e < d.n_em) pt_joule_split(&d, &b, e); }
extern "C" __global__ void k_step(pt_dims d, pt_buffers b) { const int c = blockIdx.x * blockDim.x + threadIdx.x; if (c < d.n_cells) pt_step_cell(&d, &b, c); }
extern "C" __global__ void k_slice(pt_dims d, pt_buffers b) { const int s = blockIdx.x * blockDim.x + threadIdx.x; if (s < d.S) pt_slice_reduce(&d, &b, s); }
extern "C" __global__ void k_quote(pt_dims d, pt_buffers b) { const int e = blockIdx.x * blockDim.x + threadIdx.x; if (e < d.n_el) pt_quote_element(&d, &b, e); }
extern "C" __global__ void k_release(pt_dims d, pt_buffers b) { const int p = blockIdx.x * blockDim.x + threadIdx.x; if (p < d.n_patches) pt_release_split(&d, &b, p); }
extern "C" __global__ void k_patch(pt_dims d, pt_buffers b) { const int p = blockIdx.x * blockDim.x + threadIdx.x; if (p < d.n_patches) pt_patch_report(&d, &b, p); }
#endif
