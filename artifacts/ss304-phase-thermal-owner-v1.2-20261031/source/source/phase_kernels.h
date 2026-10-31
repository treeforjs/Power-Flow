// phase_kernels.h -- the per-cell, per-EM-cell, per-element, per-face-column, per-interface and per-patch
// arithmetic of one owner step (v1.2). The CPU reference calls these functions in index loops; the NVRTC program
// wraps the same functions in extern "C" kernels (one thread per item). No atomics: every output word is written
// by exactly one item, and every sum runs over a fixed index order (energy and atom sums are exact 128-bit sums).
//
// Layout. Cells c = (s*ny + j)*nx + i. Faces f: 0 LowY, 1 HighY, 2 LowX, 3 HighX. Surface elements e: exposed
// (cell, face) pairs, sorted by patch, cell, face. Face columns fc = (s*2 + f)*nx + i for the two large faces
// f = 0, 1 (the faces that carry the melt-layer flow along x); interfaces ic = (s*2 + f)*(nx-1) + i between
// face columns i and i+1.
#pragma once
#include "phase_law.h"

struct px_dims {
    int S, nx, ny, nx_em, ny_em, rx;
    int n_cells, n_em, n_el, n_patches;
    double dx, dz;
    int n_sub, sub;               // sub-step count of this step and the current sub-step index
    double tau;                   // sub-step length (s)
    double dt;                    // whole step (s), for quotes
    int joule_mode;               // 0: uniform power density per EM cell (Hall v3), 1: sigma-weighted
    int evap_mode;                // 0: none, 1: expected HKL flux per sub-step, 2: quantized releases from the caller
    int bc_kind[4];               // per face kind: 0 adiabatic/loads only, 1 prescribed temperature (benchmarks)
    double bc_T[4];
    int evap_face[4];             // per face kind: evaporating
    double phi_floor;             // admitted minimum mass fraction of the reference mass during a step
    double phi_merge;             // end-of-step merge threshold (owner policy)
    int n_scatter;                // number of surface elements with loads/deposits this step
    // v1.1
    int flow;                     // melt-layer flow on the two large faces
    double flow_max_fraction;     // admitted fraction of a donor cell moved by flow in one sub-step
    double rho0;                  // reference density for displacements (density law at T0)
    int n_fc, n_ic;               // face columns, interfaces
    double thickness;             // foil thickness (reference), caps a face's melt layer at half of it
    // v1.2
    int flow_model;               // 1 lubrication (v1.1), 2 depth-averaged layer with inertia
    int book_work;                // lubrication: book the Lorentz work and the surface-force work as dissipated heat
    int eject;                    // depth-averaged: open-edge outflow and displacement excess go to the ejected-melt reservoir
    int open_x, open_z;           // depth-averaged: the x edges / z ends of the face are open (outflow) rather than walls
    double eta_limit;             // displacement above which layer liquid is ejected (m, eject)
    double cfl;                   // admitted a tau / distance per interface and sub-step (depth-averaged)
    double k_cap;                 // largest resolved surface wavenumber pi / min(dx, dz) (capillary wave speed)
    int n_if;                     // depth-averaged interfaces (x, z, boundaries)
    int use_acc;                  // Joule energy and Lorentz force of the step from the accumulators (sub-cycling)
    int diag;                     // compute the step diagnostics (conductivity change, window heating estimate)
    double cp_min;                // smallest heat capacity of the law (window heating estimate)
};

#define PX_SC 8   // scatter values per entry: load_E, dep_N x6, dep_E
// interface directions of the depth-averaged layer: interior x, interior z, boundaries -x, +x, -z, +z
#define PX_IF_X 0
#define PX_IF_Z 1
#define PX_IF_MX 2
#define PX_IF_PX 3
#define PX_IF_MZ 4
#define PX_IF_PZ 5
// per-face-column step accumulators of the depth-averaged layer
#define FA_WL 0           // Lorentz work on the layer (taken from the field)
#define FA_WS 1           // work of the surface forces (recoil, capillary pressure, Marangoni shear)
#define FA_Q 2            // dissipation booked to the condensed phase (internal: kinetic -> thermal)
#define FA_EJN 3          // +el, ejected atoms
#define FA_EJE 9          // ejected energy (thermal + kinetic)
#define FA_EJK 10         // ejected kinetic energy
#define PX_FC_ACC 11
#define PX_ERR_CFL 9      // a melt-layer wave would cross more than the admitted fraction of a face column in one sub-step

struct px_buffers {
    const px_law* law;
    // geometry
    const double* row_dy; const int* row_em; const int* em_row_begin; const double* gx; const double* gy;
    const double* vref; const double* mref;
    // state (current and next slot)
    const px_i128* E; const px_i128* N; const int* flags; const px_i128* fcm;   // fcm [n_fc*PX_NEL]: face mass words
    px_i128* E_out; px_i128* N_out; px_i128* fcm_out;
    int* flags_out;               // v1.2: the step copies the flags of the next slot (null: the engine copies them)
    // per-cell scratch
    double* T; double* k; double* sig; double* m; int* status;
    px_i128* joule_total; const double* em_energy;
    const double* em_JB;          // [n_em*6]: J (A/m2) x,y,z then B (T) x,y,z at the EM cell (coupling input)
    // surface elements
    const int* el_of; const int* el_cell; const int* el_face; const double* el_area; const int* el_patch; const int* patch_begin;
    const px_i128* in_load_E; const px_i128* in_dep_N; const px_i128* in_dep_E;
    px_i128* in_rel_N; px_i128* in_rel_E;
    px_i128* acc_ev_N; px_i128* acc_ev_Eref; px_i128* acc_ev_Ekin; px_i128* acc_ev_Eint; px_i128* acc_bnd_E; px_i128* acc_rad;
    double* el_quote;             // [n_el*PX_NEL]
    const double* patch_pb;       // [n_patches] back pressure (Pa)
    // face columns and interfaces (melt flow, displacement)
    int* fc_cell; double* fc_delta; double* fc_T; double* fc_eta; double* fc_p; double* fc_x; // fc_x [n_fc*PX_NEL]
    px_i128* ic_N; px_i128* ic_E; double* ic_q;      // [n_ic*PX_NEL], [n_ic], [n_ic] volumetric flux (m2/s)
    // patch / slice outputs
    const px_i128* rel_total_N; const double* rel_kin_J;
    double* patch_out; px_i128* patch_sum; px_i128* slice_sum; int* slice_status; double* slice_maxT; double* sigma_em;
    const int* sc_el; const px_i128* sc_val;
    px_i128* in_load_E_w; px_i128* in_dep_N_w; px_i128* in_dep_E_w;
    // v1.2
    const px_i128* acc_E;         // [n_em] accumulated Joule energy words (use_acc)
    const double* acc_F;          // [3 n_em] accumulated (J x B) dt (use_acc)
    const int* acc_bad;           // [n_em] nonzero: an accumulated input was not finite or negative
    double* em_F;                 // [3 n_em] Lorentz force density (J x B) used by the step (N/m3)
    double* em_diag;              // [2 n_em] window heating estimate (K) at the start, conductivity at the start
    double* slice_diag;           // [4 S] largest relative conductivity change, largest window heating estimate, accumulator flag,
                                  //       largest layer signal speed |u| + a
    const double* fcV; double* fcV_out;     // [3 n_fc] layer momentum x, z (kg m/s) and layer mass (kg), state
    const px_i128* fcK; px_i128* fcK_out;   // [n_fc] layer kinetic energy word, state
    const int* if_a; const int* if_b; const int* if_dir; const int* fc_if;   // fc_if [4 n_fc]: -x, +x, -z, +z interfaces
    const int* slice_zn;          // [2 S] slice of the same body at k - 1, k + 1 (-1: none)
    int* fc_lay;                  // cells of the layer, -1: blocked column (surface cell outside the face's half)
    int* fc_nw;                   // cells of the column written by the layer pass this sub-step
    double* fc_m; double* fc_h; double* fc_rho; double* fc_mu; double* fc_dsdT; double* fc_sig;
    double* fc_f;                 // [3 n_fc] f_x, f_z, inward normal force (<= 0)
    double* fc_u;                 // [2 n_fc] velocity at the start of the sub-step
    double* fc_un;                // [2 n_fc] velocity at the end of the sub-step
    double* fc_a; double* fc_ps; double* fc_g; double* fc_re; double* fc_s;   // wave speed, surface pressure, load, Re, |u| + a
    double* if_dm; double* if_F; int* if_st; px_i128* if_N; px_i128* if_E; px_i128* if_K;
    px_i128* cflow_N; px_i128* cflow_E; px_i128* fc_netN;
    px_i128* acc_fc;              // [PX_FC_ACC n_fc]
    px_i128* ic_W;                // [n_ic] lubrication heat word of the sub-step
    px_i128* acc_ic;              // [2 n_ic] lubrication Lorentz work, surface-force work (accumulated)
    px_i128* el_dE; px_i128* el_dN; int* el_st; int* el_bad;   // [n_el] surface physics of the sub-step per element
    px_i128* row_sum; int* row_info; double* row_maxT;          // [S ny] row partial sums of the slice reduction
};

// patch_out layout
#define XO_AREA 0
#define XO_MEAN_T 1
#define XO_MAX_T 2
#define XO_FL 3
#define XO_MELT 4
#define XO_D 5            // +el, expected atoms over dt
#define XO_J 11           // +el, flux atoms/m2/s
#define XO_RHO 17         // mean surface density
#define XO_ETA_MEAN 18    // mean displacement (m, outward)
#define XO_ETA_MAX 19     // max |displacement|
#define XO_RECOIL 20      // mean recoil pressure (Pa)
#define XO_SPEED 21       // mean |melt speed| at the patch interfaces (m/s)
#define XO_X 22           // +el, mean surface atom fraction
#define XO_EINT 28        // released internal energy (quantized mode, J)
#define XO_RE 29          // largest melt-layer Reynolds number rho |q| / mu at the patch interfaces
#define XO_UX 30          // mean layer velocity along x (depth-averaged model, m/s)
#define XO_UZ 31          // mean layer velocity along z
#define XO_H 32           // mean layer thickness (depth-averaged model, m)
#define XO_ETA_NEG 33     // largest inward displacement (crater depth, m)
#define PX_PATCH_OUT 34
// patch_sum layout
#define XS_EVN 0          // +el
#define XS_EREF 6
#define XS_EKIN 7
#define XS_EINT 8
#define XS_BND 9
#define XS_LOAD 10
#define XS_DEP 11
#define XS_REL 12
#define XS_RAD 13
#define XS_EJN 14         // +el, ejected atoms
#define XS_EJE 20         // ejected energy (thermal + kinetic)
#define XS_EJK 21         // ejected kinetic energy
#define XS_WL 22          // Lorentz work on the melt (the field-energy exchange)
#define XS_WS 23          // surface-force work on the melt
#define XS_Q 24           // melt dissipation booked to the condensed phase
#define PX_PATCH_SUM 25
// slice sums: E, N x6, layer kinetic energy, accumulated Joule energy
#define SL_K (1 + PX_NEL)
#define SL_J (2 + PX_NEL)
#define PX_SLICE_SUM (3 + PX_NEL)

PX_HD int px_neighbor(const px_dims* d, int c, int f) {
    const int i = c % d->nx, j = (c / d->nx) % d->ny;
    if (f == 0) return j > 0 ? c - d->nx : -1;
    if (f == 1) return j + 1 < d->ny ? c + d->nx : -1;
    if (f == 2) return i > 0 ? c - 1 : -1;
    return i + 1 < d->nx ? c + 1 : -1;
}
PX_HD double px_face_g(const px_dims* d, const px_buffers* b, int c, int f) {
    const int j = (c / d->nx) % d->ny;
    if (f == 0) return b->gy[j - 1];
    if (f == 1) return b->gy[j];
    return b->gx[j];
}
PX_HD double px_face_area(const px_dims* d, const px_buffers* b, int c, int f) {
    const int j = (c / d->nx) % d->ny;
    return f < 2 ? d->dx * d->dz : b->row_dy[j] * d->dz;
}
PX_HD double px_face_half(const px_dims* d, const px_buffers* b, int c, int f) {
    const int j = (c / d->nx) % d->ny;
    return f < 2 ? 0.5 * b->row_dy[j] : 0.5 * d->dx;
}
PX_HD double px_face_energy(double Ta, double ka, double Tb, double kb, double g, double tau) {
    const double kh = (ka == 0.0 || kb == 0.0) ? 0.0 : 2.0 * ka * kb / (ka + kb);
    return kh * g * (Ta - Tb) * tau;
}
PX_HD px_i128 px_part(px_i128 total, int n, int s, int* ov) {
    px_u32 r = 0;
    px_i128 q = px_divmod_u32(total, (px_u32)n, &r, ov);
    if ((px_u32)s < r) { px_i128 one; one.lo = 1; one.hi = 0; q = px_add(q, one, ov); }
    return q;
}
PX_HD void px_set_status(int* st, int code) { if (*st == PX_OK) *st = code; }
PX_HD double px_cell_mass(const px_law* L, const px_i128* N, int c, double* atoms) {
    double m = 0.0;
    for (int el = 0; el < PX_NEL; ++el) { atoms[el] = px_to_double(N[PX_NEL * c + el], PX_N_INV); m += atoms[el] * L->atom_mass[el]; }
    return m;
}
PX_HD void px_fractions(const double* atoms, double* x) {
    double tot = 0.0;
    for (int el = 0; el < PX_NEL; ++el) tot += atoms[el];
    for (int el = 0; el < PX_NEL; ++el) x[el] = atoms[el] / tot;
}
PX_HD int px_fc_index(const px_dims* d, int s, int f, int i) { return (s * 2 + f) * d->nx + i; }

// 1. Thermodynamic state of one cell from (E, N).
PX_HD void px_thermo(const px_dims* d, const px_buffers* b, int c) {
    (void)d;
    int st = PX_OK;
    if (b->flags[c] & 1) { b->T[c] = 0; b->k[c] = 0; b->sig[c] = 0; b->m[c] = 0; return; }
    double atoms[PX_NEL];
    const double m = px_cell_mass(b->law, b->N, c, atoms);
    const double h = px_to_double(b->E[c], PX_E_INV) / m;
    const double T = px_T_of_h(b->law, h, &st);
    b->T[c] = T; b->m[c] = m;
    b->k[c] = px_k(b->law, T);
    b->sig[c] = 1.0 / px_rhoe(b->law, T);
    if (b->law->T_explosion > 0.0 && T >= b->law->T_explosion) px_set_status(&st, PX_ERR_PHASE_EXPLOSION);
    if (b->status[c] == PX_OK) b->status[c] = st;
}

PX_HD void px_scatter(const px_dims* d, const px_buffers* b, int q) {
    (void)d;
    const int e = b->sc_el[q];
    b->in_load_E_w[e] = b->sc_val[PX_SC * q + 0];
    for (int el = 0; el < PX_NEL; ++el) b->in_dep_N_w[PX_NEL * e + el] = b->sc_val[PX_SC * q + 1 + el];
    b->in_dep_E_w[e] = b->sc_val[PX_SC * q + 1 + PX_NEL];
}

PX_HD double px_sigma_value(const px_dims* d, const px_buffers* b, int e) {
    const int per_slice = d->nx_em * d->ny_em;
    const int s = e / per_slice, je = (e % per_slice) / d->nx_em, ie = e % d->nx_em;
    double sv = 0.0, v = 0.0;
    for (int j = b->em_row_begin[je]; j < b->em_row_begin[je + 1]; ++j) for (int i = ie * d->rx; i < (ie + 1) * d->rx; ++i) {
        const int c = (s * d->ny + j) * d->nx + i;
        if (b->flags[c] & 1) continue;
        sv += b->sig[c] * b->vref[c]; v += b->vref[c];
    }
    return v > 0.0 ? sv / v : 0.0;
}
PX_HD void px_sigma_em(const px_dims* d, const px_buffers* b, int e) { b->sigma_em[e] = px_sigma_value(d, b, e); }
// v1.2: conductivity of the EM cell at the start of the step (diagnostics: its change over the step bounds the error of
// the field solve that used the start value for the whole window).
PX_HD void px_sigma_start(const px_dims* d, const px_buffers* b, int e) { b->em_diag[2 * e + 1] = px_sigma_value(d, b, e); }
// v1.2: Lorentz force density of the step: instantaneous J x B (single step) or the window mean of the accumulated
// (J x B) dt (sub-cycling).
PX_HD void px_em_prepare(const px_dims* d, const px_buffers* b, int e) {
    double* o = b->em_F + 3 * e;
    if (d->use_acc) {
        for (int k = 0; k < 3; ++k) o[k] = b->acc_F[3 * e + k] / d->dt;
    } else {
        const double* v = b->em_JB + 6 * e;
        o[0] = v[1] * v[5] - v[2] * v[4];
        o[1] = v[2] * v[3] - v[0] * v[5];
        o[2] = v[0] * v[4] - v[1] * v[3];
    }
}
// v1.2: one field step into the accumulators of EM cell e: Joule energy word += P dt (exact sum of the converted
// field-step energies), (J x B) dt += (J x B) dt; bad[e] marks a nonfinite or negative input.
PX_HD void px_accumulate_one(int e, const double* P, const double* JB, double dt, px_i128* accE, double* accF, int* bad) {
    int b = 0, ov = 0;
    const double p = P[e];
    if (!(px_finite(p) && p >= 0.0)) bad[e] = 1;
    else accE[e] = px_add(accE[e], px_from_double(p * dt, PX_E_POW2, &b), &ov);
    if (JB) {
        const double* v = JB + 6 * e;
        const double f0 = v[1] * v[5] - v[2] * v[4], f1 = v[2] * v[3] - v[0] * v[5], f2 = v[0] * v[4] - v[1] * v[3];
        if (!(px_finite(f0) && px_finite(f1) && px_finite(f2))) bad[e] = 1;
        else { accF[3 * e] += f0 * dt; accF[3 * e + 1] += f1 * dt; accF[3 * e + 2] += f2 * dt; }
    }
    if (b || ov) bad[e] = 1;
}

// 2. Exact partition of one EM cell's Joule energy over its thermal cells.
PX_HD void px_joule_split(const px_dims* d, const px_buffers* b, int e) {
    const int per_slice = d->nx_em * d->ny_em;
    const int s = e / per_slice, je = (e % per_slice) / d->nx_em, ie = e % d->nx_em;
    const int j0 = b->em_row_begin[je], j1 = b->em_row_begin[je + 1], i0 = ie * d->rx, i1 = i0 + d->rx;
    const double energy = d->use_acc ? px_to_double(b->acc_E[e], PX_E_INV) : b->em_energy[e];
    int bad = 0, ov = 0;
    double W = 0.0, wmax = -1.0, msum = 0.0; int cmax = -1;
    for (int j = j0; j < j1; ++j) for (int i = i0; i < i1; ++i) {
        const int c = (s * d->ny + j) * d->nx + i;
        b->joule_total[c] = px_zero();
        if (b->flags[c] & 1) continue;
        const double w = d->joule_mode == 1 ? b->sig[c] * b->vref[c] : b->vref[c];
        W += w;
        msum += b->m[c];
        if (w > wmax) { wmax = w; cmax = c; }
    }
    if (d->diag) b->em_diag[2 * e] = msum > 0.0 ? energy / (msum * d->cp_min) : 0.0;
    if (energy == 0.0 && !(d->use_acc && !px_is_zero(b->acc_E[e]))) return;
    if (cmax < 0 || !(W > 0.0)) { b->status[(s * d->ny + j0) * d->nx + i0] = PX_ERR_NONFINITE; return; }
    const px_i128 total = d->use_acc ? b->acc_E[e] : px_from_double(energy, PX_E_POW2, &bad);
    px_i128 used = px_zero();
    for (int j = j0; j < j1; ++j) for (int i = i0; i < i1; ++i) {
        const int c = (s * d->ny + j) * d->nx + i;
        if ((b->flags[c] & 1) || c == cmax) continue;
        const double w = d->joule_mode == 1 ? b->sig[c] * b->vref[c] : b->vref[c];
        const px_i128 q = px_from_double(energy * (w / W), PX_E_POW2, &bad);
        b->joule_total[c] = q;
        used = px_add(used, q, &ov);
    }
    const px_i128 rest = px_sub(total, used, &ov);
    b->joule_total[cmax] = rest;
    if (bad || ov || px_is_neg(rest)) px_set_status(&b->status[cmax], PX_ERR_FIXED_RANGE);
}

// 2b. v1.2 start-of-step pass of one EM cell (fused): Lorentz force of the step, conductivity at the start
//     (diagnostics), Joule partition.
PX_HD void px_em_start(const px_dims* d, const px_buffers* b, int e) {
    px_em_prepare(d, b, e);
    if (d->diag) px_sigma_start(d, b, e);
    px_joule_split(d, b, e);
}
// v1.2: zero the per-step arrays (status, element inputs and accumulators, quotes, patch outputs, layer and
//    lubrication accumulators) in one pass over t < max(n_cells, n_el, n_patches, n_fc, n_ic).
PX_HD void px_zero_step(const px_dims* d, const px_buffers* b, int t) {
    if (t < d->n_cells) b->status[t] = 0;
    if (t < d->n_el) {
        b->in_load_E_w[t] = px_zero(); b->in_dep_E_w[t] = px_zero(); b->in_rel_E[t] = px_zero();
        b->acc_ev_Eref[t] = px_zero(); b->acc_ev_Ekin[t] = px_zero(); b->acc_ev_Eint[t] = px_zero(); b->acc_bnd_E[t] = px_zero(); b->acc_rad[t] = px_zero();
        for (int el = 0; el < PX_NEL; ++el) {
            b->in_dep_N_w[PX_NEL * t + el] = px_zero(); b->in_rel_N[PX_NEL * t + el] = px_zero(); b->acc_ev_N[PX_NEL * t + el] = px_zero();
            b->el_quote[PX_NEL * t + el] = 0.0;
        }
    }
    if (t < d->n_patches) {
        for (int q = 0; q < PX_PATCH_OUT; ++q) b->patch_out[PX_PATCH_OUT * t + q] = 0.0;
        for (int q = 0; q < PX_PATCH_SUM; ++q) b->patch_sum[PX_PATCH_SUM * t + q] = px_zero();
    }
    if (t < d->n_fc && d->flow && d->flow_model == 2) { for (int q = 0; q < PX_FC_ACC; ++q) b->acc_fc[PX_FC_ACC * t + q] = px_zero(); b->fc_s[t] = 0.0; }
    if (t < d->n_ic && d->flow && d->book_work) { b->acc_ic[2 * t] = px_zero(); b->acc_ic[2 * t + 1] = px_zero(); }
}

// 3a. Face column: surface cell, melt-layer thickness, surface temperature and composition, displacement.
//     eta = (face mass gained) / (rho0 A) + expansion / 2, expansion = sum_c m_c (1/rho(T_c) - 1/rho0) / A.
PX_HD void px_fc_prepare(const px_dims* d, const px_buffers* b, int fc) {
    const int i = fc % d->nx, sf = fc / d->nx, s = sf / 2, f = sf % 2;
    const px_law* L = b->law;
    int sc = -1;
    for (int jj = 0; jj < d->ny; ++jj) {
        const int j = f == 0 ? jj : d->ny - 1 - jj;
        const int c = (s * d->ny + j) * d->nx + i;
        if (!(b->flags[c] & 1)) { sc = c; break; }
    }
    b->fc_cell[fc] = sc;
    double delta = 0.0, expansion = 0.0, x[PX_NEL];
    for (int el = 0; el < PX_NEL; ++el) x[el] = 0.0;
    double Ts = 0.0;
    if (sc >= 0) {
        Ts = b->T[sc];
        double atoms[PX_NEL];
        (void)px_cell_mass(L, b->N, sc, atoms);
        px_fractions(atoms, x);
        int c = sc;
        const int inward = f == 0 ? 1 : 0;
        while (c >= 0 && !(b->flags[c] & 1)) {
            const double fl = px_fl(L, b->T[c]);
            if (fl <= 0.0) break;
            // physical thickness of the liquid in this cell from its mass and the liquid density
            const double rho_c = L->kind == 0 ? px_density(L, b->T[c]) : (L->c_rho > 0.0 ? L->c_rho : d->rho0);
            delta += fl * b->m[c] / (rho_c * d->dx * d->dz);
            c = px_neighbor(d, c, inward);
        }
        // a face's layer is at most half the foil: a fully molten column is split between its two faces
        if (delta > 0.5 * d->thickness) delta = 0.5 * d->thickness;
        if (L->kind == 0 && d->rho0 > 0.0) {
            for (int j = 0; j < d->ny; ++j) {
                const int cc = (s * d->ny + j) * d->nx + i;
                if (b->flags[cc] & 1) continue;
                expansion += b->m[cc] * (1.0 / px_density(L, b->T[cc]) - 1.0 / d->rho0);
            }
            expansion /= d->dx * d->dz;
        }
    }
    double gained = 0.0;
    for (int el = 0; el < PX_NEL; ++el) gained += px_to_double(b->fcm[PX_NEL * fc + el], PX_N_INV) * L->atom_mass[el];
    const double rho = d->rho0 > 0.0 ? d->rho0 : 1.0;
    b->fc_delta[fc] = delta; b->fc_T[fc] = Ts;
    b->fc_eta[fc] = gained / (rho * d->dx * d->dz) + 0.5 * expansion;
    for (int el = 0; el < PX_NEL; ++el) b->fc_x[PX_NEL * fc + el] = x[el];
}
// 3b. Face column pressure: recoil (Knudsen layer, back pressure of the patch) - surface tension x curvature
//     - mean normal Lorentz load (n . (J x B)) delta / 2.
PX_HD double px_patch_pb(const px_buffers* b, int c, int f) {
    const int e = b->el_of[4 * c + f];
    return e >= 0 ? b->patch_pb[b->el_patch[e]] : 0.0;
}
PX_HD void px_surface_props(const px_law* L, double T, double* sigma, double* dsdT, double* mu) {
    if (L->kind == 0) { px_surface_tension(L, T, sigma, dsdT); *mu = px_viscosity(T); }
    else { *sigma = L->c_sigma; *dsdT = L->c_dsdT; *mu = L->c_mu; }
}
// Lorentz force density of the EM cell at the surface row of face f (em_F: (J x B)_x, _y, _z, filled per step).
PX_HD void px_em_force(const px_dims* d, const px_buffers* b, int s, int f, int i, double* fx, double* fn) {
    const int je = f == 0 ? 0 : d->ny_em - 1, ie = i / d->rx;
    const double* v = b->em_F + 3 * ((s * d->ny_em + je) * d->nx_em + ie);
    *fx = v[0];
    *fn = f == 0 ? -v[1] : v[1];                     // outward-normal component
}
PX_HD double px_em_force_z(const px_dims* d, const px_buffers* b, int s, int f, int i) {
    const int je = f == 0 ? 0 : d->ny_em - 1, ie = i / d->rx;
    return b->em_F[3 * ((s * d->ny_em + je) * d->nx_em + ie) + 2];
}
PX_HD void px_fc_pressure(const px_dims* d, const px_buffers* b, int fc) {
    const int i = fc % d->nx, sf = fc / d->nx, s = sf / 2, f = sf % 2;
    const int sc = b->fc_cell[fc];
    if (sc < 0) { b->fc_p[fc] = 0.0; return; }
    int st = 0;
    const double T = b->fc_T[fc];
    const double pr = b->law->vapor_enabled ? px_recoil(b->law, T, b->fc_x + PX_NEL * fc, px_patch_pb(b, sc, f), &st) : 0.0;
    double sigma, dsdT, mu;
    px_surface_props(b->law, T, &sigma, &dsdT, &mu);
    const double em = b->fc_eta[fc];
    const double el = i > 0 ? b->fc_eta[fc - 1] : em, er = i + 1 < d->nx ? b->fc_eta[fc + 1] : em;
    const double curvature = (er - 2.0 * em + el) / (d->dx * d->dx);
    double fx, fn;
    px_em_force(d, b, s, f, i, &fx, &fn);
    b->fc_p[fc] = pr - sigma * curvature - fn * b->fc_delta[fc] * 0.5;
    if (d->book_work) b->fc_ps[fc] = -(fn * b->fc_delta[fc] * 0.5);   // normal-Lorentz part of the layer pressure
}
// 3c. Interface between face columns i and i+1: lubrication flux
//     q = delta^3/(3 mu) (f_x - dp/dx) + delta^2/(2 mu) dsigma/dT dT_s/dx   (m2/s, + toward +x), with delta the
//     mean of the two columns' liquid thicknesses (from mass and liquid density); no flow unless both are molten
//     and the exact transfer (donor composition and specific energy) over tau.
PX_HD void px_ic_flux(const px_dims* d, const px_buffers* b, int ic) {
    const int i = ic % (d->nx - 1), sf = ic / (d->nx - 1), s = sf / 2, f = sf % 2;
    const int fa = px_fc_index(d, s, f, i), fb = fa + 1;
    const int ca = b->fc_cell[fa], cb = b->fc_cell[fb];
    const px_law* L = b->law;
    for (int el = 0; el < PX_NEL; ++el) b->ic_N[PX_NEL * ic + el] = px_zero();
    b->ic_E[ic] = px_zero(); b->ic_q[ic] = 0.0;
    if (d->book_work) b->ic_W[ic] = px_zero();
    if (ca < 0 || cb < 0) return;
    if (!(b->fc_delta[fa] > 0.0) || !(b->fc_delta[fb] > 0.0) || px_fl(L, b->T[ca]) < 1.0 || px_fl(L, b->T[cb]) < 1.0) return;
    const double delta = 0.5 * (b->fc_delta[fa] + b->fc_delta[fb]);   // layer thickness at the interface
    double sa, dsa, ma, sb, dsb, mb;
    px_surface_props(L, b->fc_T[fa], &sa, &dsa, &ma);
    px_surface_props(L, b->fc_T[fb], &sb, &dsb, &mb);
    const double mu = 0.5 * (ma + mb), dsdT = 0.5 * (dsa + dsb);
    double fxa, fna, fxb, fnb;
    px_em_force(d, b, s, f, i, &fxa, &fna);
    px_em_force(d, b, s, f, i + 1, &fxb, &fnb);
    const double fx = 0.5 * (fxa + fxb);
    const double dpdx = (b->fc_p[fb] - b->fc_p[fa]) / d->dx;
    const double tau_m = dsdT * (b->fc_T[fb] - b->fc_T[fa]) / d->dx;
    const double q = delta * delta * delta / (3.0 * mu) * (fx - dpdx) + delta * delta / (2.0 * mu) * tau_m;
    b->ic_q[ic] = q;
    if (q == 0.0) return;
    const int donor = q > 0.0 ? ca : cb;
    const double rho = L->kind == 0 ? px_density(L, b->T[donor]) : L->c_rho;
    const double dm = rho * (q > 0.0 ? q : -q) * d->dz * d->tau;
    const double phi = dm / b->m[donor];
    int bad = 0;
    if (!(phi <= d->flow_max_fraction)) { b->status[donor] = PX_ERR_FLOW; return; }
    int ov = 0;
    for (int el = 0; el < PX_NEL; ++el) {
        px_i128 qn = px_from_double(px_to_double(b->N[PX_NEL * donor + el], PX_N_INV) * phi, PX_N_POW2, &bad);
        b->ic_N[PX_NEL * ic + el] = q > 0.0 ? qn : px_neg(qn, &ov);
    }
    px_i128 qe = px_from_double(px_to_double(b->E[donor], PX_E_INV) * phi, PX_E_POW2, &bad);
    b->ic_E[ic] = q > 0.0 ? qe : px_neg(qe, &ov);
    if (d->book_work) {
        // Stokes layer: all work is dissipated at once. Per unit area the dissipation is G q + tau_m u_s (G = f_x - dp/dx,
        // u_s the surface speed of the profile); its Lorentz part (tangential force and the normal-load pressure) is
        // taken from the field, the rest is the work of the recoil, capillary and Marangoni forces.
        const double dpn = (b->fc_ps[fb] - b->fc_ps[fa]) / d->dx;
        const double us = delta * delta / (2.0 * mu) * (fx - dpdx) + tau_m * delta / mu;
        const double area = d->dx * d->dz * d->tau;
        const px_i128 wl = px_from_double((fx - dpn) * q * area, PX_E_POW2, &bad);
        const px_i128 ws = px_from_double((-(dpdx - dpn) * q + tau_m * us) * area, PX_E_POW2, &bad);
        b->ic_W[ic] = px_add(wl, ws, &ov);
        b->acc_ic[2 * ic] = px_add(b->acc_ic[2 * ic], wl, &ov);
        b->acc_ic[2 * ic + 1] = px_add(b->acc_ic[2 * ic + 1], ws, &ov);
    }
    if (bad || ov) b->status[donor] = PX_ERR_FIXED_RANGE;
}

// ---------------------------------------------------------------------------------------------------------
// v1.2 depth-averaged melt layer (flow_model 2). The layer of face column fc = (s, f, i) is the contiguous run of
// fully liquid cells from the surface cell inward inside the face's half of the column (rows [0, ny/2) for LowY,
// [ny/2, ny) for HighY). Column state: momentum M = (M_x, M_z) (kg m/s), the layer mass at the end of the last
// sub-step and the kinetic-energy word K. Shallow-layer equations on the face (x along the width, z along the
// height through the slices of one body):
//   d(m)/dt + div(m u) = 0,
//   d(m u)/dt + div(m u u + P) = A [h f_t + 3/2 tau_s - h grad p_s] - 3 mu u A / h,
// with P = rho g h^2 / 2 the depth-integrated load of the inward normal Lorentz force (rho g = -min(f_n, 0)),
// p_s = recoil - sigma kappa (kappa from the displacement field), tau_s = dsigma/dT grad T_s and the bed friction
// 3 mu u / h - tau_s / 2 of the Poiseuille-Couette profile (steady limits: u = f h^2 / (3 mu), u = tau_s h / (2 mu)).
// Rusanov fluxes with hydrostatic reconstruction (bed b = eta - h), semi-implicit friction, explicit sub-steps
// limited by (|u| + a) tau <= cfl dx. Every transfer moves exact words: a donor gives the fraction dm / m_layer of
// each of its layer cells' atoms and energy.
PX_HD int px_in_half(const px_dims* d, int c, int f) { const int j = (c / d->nx) % d->ny; return f == 0 ? (j < d->ny / 2) : (j >= d->ny / 2); }
PX_HD int px_layer_step(const px_dims* d, int f) { return f == 0 ? d->nx : -d->nx; }
PX_HD double px_cell_rho(const px_dims* d, const px_law* L, double T) { return L->kind == 0 ? px_density(L, T) : (L->c_rho > 0.0 ? L->c_rho : d->rho0); }
PX_HD int px_fc_neighbor(const px_dims* d, const px_buffers* b, int fc, int side) {
    const int i = fc % d->nx, sf = fc / d->nx, s = sf / 2, f = sf % 2;
    if (side == 0) return i > 0 ? fc - 1 : -1;
    if (side == 1) return i + 1 < d->nx ? fc + 1 : -1;
    const int t = b->slice_zn[2 * s + (side == 2 ? 0 : 1)];
    return t >= 0 ? px_fc_index(d, t, f, i) : -1;
}
// central difference along x (side pair 0/1, spacing dx) or z (2/3, dz); one-sided at a missing neighbour
PX_HD double px_fc_grad(const px_dims* d, const px_buffers* b, int fc, const double* v, int stride, int comp, int dir) {
    const int m = px_fc_neighbor(d, b, fc, 2 * dir), p = px_fc_neighbor(d, b, fc, 2 * dir + 1);
    const double h = dir == 0 ? d->dx : d->dz, v0 = v[stride * fc + comp];
    if (m >= 0 && p >= 0) return (v[stride * p + comp] - v[stride * m + comp]) / (2.0 * h);
    if (p >= 0) return (v[stride * p + comp] - v0) / h;
    if (m >= 0) return (v0 - v[stride * m + comp]) / h;
    return 0.0;
}

// 3d. Layer of one face column at the start of a sub-step.
PX_HD void px_fc_prepare2(const px_dims* d, const px_buffers* b, int fc) {
    px_fc_prepare(d, b, fc);
    const int i = fc % d->nx, sf = fc / d->nx, s = sf / 2, f = sf % 2;
    const px_law* L = b->law;
    const int sc = b->fc_cell[fc];
    int n = 0, lay = -1;
    double m = 0.0, h = 0.0, mT = 0.0;
    if (sc >= 0 && px_in_half(d, sc, f)) {
        int c = sc;
        for (;;) {
            if ((b->flags[c] & 1) || px_fl(L, b->T[c]) < 1.0) break;
            m += b->m[c]; mT += b->m[c] * b->T[c];
            h += b->m[c] / (px_cell_rho(d, L, b->T[c]) * d->dx * d->dz);
            ++n;
            const int nxt = c + px_layer_step(d, f);
            if (!px_in_half(d, nxt, f)) break;
            c = nxt;
        }
        lay = n;
    }
    double sigma = 0.0, dsdT = 0.0, mu = 0.0;
    const double Ts = sc >= 0 ? b->T[sc] : L->T_min;
    px_surface_props(L, Ts, &sigma, &dsdT, &mu);
    if (n > 0 && L->kind == 0) mu = px_viscosity(mT / m);      // viscosity at the layer's mass-mean temperature
    const double rho = n > 0 ? m / (h * d->dx * d->dz) : px_cell_rho(d, L, Ts);
    double fx, fn;
    px_em_force(d, b, s, f, i, &fx, &fn);
    const double fz = px_em_force_z(d, b, s, f, i);
    const double fin = fn < 0.0 ? fn : 0.0;                    // inward normal load only
    double ux = 0.0, uz = 0.0;
    if (n > 0) {
        double Mx = b->fcV[3 * fc], Mz = b->fcV[3 * fc + 1];
        const double mref = b->fcV[3 * fc + 2];
        // liquid that left the layer by solidification or evaporation takes its share of the momentum
        if (mref > 0.0 && m < mref) { const double r = m / mref; Mx *= r; Mz *= r; }
        ux = Mx / m; uz = Mz / m;
    }
    b->fc_lay[fc] = lay; b->fc_m[fc] = m; b->fc_h[fc] = h; b->fc_rho[fc] = rho; b->fc_mu[fc] = mu;
    b->fc_sig[fc] = sigma; b->fc_dsdT[fc] = dsdT;
    b->fc_f[3 * fc] = fx; b->fc_f[3 * fc + 1] = fz; b->fc_f[3 * fc + 2] = fin;
    b->fc_g[fc] = -fin / rho;
    b->fc_u[2 * fc] = ux; b->fc_u[2 * fc + 1] = uz;
}
// 3e. Surface pressure (recoil - sigma kappa, curvature along x and z) and wave speed of one face column.
PX_HD void px_fc_pressure2(const px_dims* d, const px_buffers* b, int fc) {
    const int sf = fc / d->nx, f = sf % 2;
    const int sc = b->fc_cell[fc];
    double ps = 0.0;
    if (sc >= 0) {
        int st = 0;
        const double pr = b->law->vapor_enabled ? px_recoil(b->law, b->fc_T[fc], b->fc_x + PX_NEL * fc, px_patch_pb(b, sc, f), &st) : 0.0;
        const double em = b->fc_eta[fc];
        double e[4];
        for (int side = 0; side < 4; ++side) { const int nb = px_fc_neighbor(d, b, fc, side); e[side] = nb >= 0 ? b->fc_eta[nb] : em; }
        const double kappa = (e[1] - 2.0 * em + e[0]) / (d->dx * d->dx) + (e[3] - 2.0 * em + e[2]) / (d->dz * d->dz);
        ps = pr - b->fc_sig[fc] * kappa;
    }
    b->fc_ps[fc] = ps;
    const double h = b->fc_h[fc];
    b->fc_a[fc] = h > 0.0 ? px_sqrt(b->fc_g[fc] * h + b->fc_sig[fc] * h * d->k_cap * d->k_cap / b->fc_rho[fc]) : 0.0;
}
// 3f. Interface flux (interior: Rusanov with hydrostatic reconstruction; boundary: mirror wall, or transmissive
//     outflow at an open edge with ejection) and the exact words leaving the donor layer.
PX_HD void px_if_flux2(const px_dims* d, const px_buffers* b, int q) {
    const int a = b->if_a[q], c2 = b->if_b[q], dir = b->if_dir[q];
    b->if_dm[q] = 0.0; b->if_F[3 * q] = 0.0; b->if_F[3 * q + 1] = 0.0; b->if_F[3 * q + 2] = 0.0; b->if_st[q] = 0;
    for (int el = 0; el < PX_NEL; ++el) b->if_N[PX_NEL * q + el] = px_zero();
    b->if_E[q] = px_zero(); b->if_K[q] = px_zero();
    const int nd = (dir == PX_IF_X || dir == PX_IF_MX || dir == PX_IF_PX) ? 0 : 1, td = 1 - nd;
    const double len = nd == 0 ? d->dz : d->dx, dist = nd == 0 ? d->dx : d->dz, k = len * d->tau;
    const int la = b->fc_lay[a];
    double dm = 0.0;
    if (c2 >= 0 && (la < 0) != (b->fc_lay[c2] < 0)) {
        // one column blocked (its half of the column is gone): mirror wall for the other
        const int w = la < 0 ? c2 : a;
        const double hw = b->fc_h[w];
        if (!(hw > 0.0)) return;
        const double sgn = w == a ? 1.0 : -1.0;
        const double rw = b->fc_rho[w], un = sgn * b->fc_u[2 * w + nd], qw = rw * hw, sp = px_fabs(un) + b->fc_a[w];
        if (sp * d->tau > d->cfl * dist) b->if_st[q] = PX_ERR_CFL;
        const double Fn = (qw * un * un + 0.5 * rw * b->fc_g[w] * hw * hw + sp * qw * un) * k;
        if (w == a) b->if_F[3 * q] = Fn; else b->if_F[3 * q + 1] = Fn;   // pushes the wet column away from the wall
        return;
    }
    if (c2 >= 0) {
        if (la < 0) return;
        const double ha = b->fc_h[a], hb = b->fc_h[c2];
        if (!(ha > 0.0) && !(hb > 0.0)) return;
        const double ba = b->fc_eta[a] - ha, bb = b->fc_eta[c2] - hb, bm = ba > bb ? ba : bb;
        double hsa = ha + ba - bm, hsb = hb + bb - bm;
        if (hsa < 0.0) hsa = 0.0;
        if (hsa > ha) hsa = ha;
        if (hsb < 0.0) hsb = 0.0;
        if (hsb > hb) hsb = hb;
        const double g = 0.5 * (b->fc_g[a] + b->fc_g[c2]);
        const double ra = b->fc_rho[a], rb = b->fc_rho[c2];
        const double una = b->fc_u[2 * a + nd], uta = b->fc_u[2 * a + td], unb = b->fc_u[2 * c2 + nd], utb = b->fc_u[2 * c2 + td];
        const double qa = ra * hsa, qb = rb * hsb;
        const double sa = px_fabs(una) + b->fc_a[a], sb = px_fabs(unb) + b->fc_a[c2], sp = sa > sb ? sa : sb;
        if (sp * d->tau > d->cfl * dist) b->if_st[q] = PX_ERR_CFL;
        const double Fm = 0.5 * (qa * una + qb * unb) - 0.5 * sp * (qb - qa);
        const double pa = 0.5 * ra * g * hsa * hsa, pb = 0.5 * rb * g * hsb * hsb;
        const double Fn = 0.5 * (qa * una * una + pa + qb * unb * unb + pb) - 0.5 * sp * (qb * unb - qa * una);
        const double Ft = 0.5 * (qa * una * uta + qb * unb * utb) - 0.5 * sp * (qb * utb - qa * uta);
        b->if_F[3 * q] = (Fn + 0.5 * ra * g * (ha * ha - hsa * hsa)) * k;
        b->if_F[3 * q + 1] = (Fn + 0.5 * rb * g * (hb * hb - hsb * hsb)) * k;
        b->if_F[3 * q + 2] = Ft * k;
        dm = Fm * k;
    } else {
        if (la < 0) return;
        const double ha = b->fc_h[a];
        if (!(ha > 0.0)) return;
        const double sgn = (dir == PX_IF_MX || dir == PX_IF_MZ) ? -1.0 : 1.0;
        const double ra = b->fc_rho[a], g = b->fc_g[a];
        const double un = sgn * b->fc_u[2 * a + nd], ut = b->fc_u[2 * a + td];
        const double qa = ra * ha, sp = px_fabs(un) + b->fc_a[a];
        if (sp * d->tau > d->cfl * dist) b->if_st[q] = PX_ERR_CFL;
        const int open = d->eject && (nd == 0 ? d->open_x : d->open_z);
        const double pa = 0.5 * ra * g * ha * ha;
        double Fm = 0.0, Fn, Ft = 0.0;
        if (open && un > 0.0) { Fm = qa * un; Fn = qa * un * un + pa; Ft = qa * un * ut; }
        else Fn = qa * un * un + pa + sp * qa * un;
        b->if_F[3 * q] = Fn * k;        // outward normal impulse leaving column a
        b->if_F[3 * q + 2] = Ft * k;
        dm = Fm * k;
    }
    b->if_dm[q] = dm;
    if (dm == 0.0) return;
    const int donor = dm > 0.0 ? a : c2;
    const double phi = px_fabs(dm) / b->fc_m[donor];
    int bad = 0, ov = 0;
    px_i128 sN[PX_NEL], sE = px_zero();
    for (int el = 0; el < PX_NEL; ++el) sN[el] = px_zero();
    const int n = b->fc_lay[donor], stp = px_layer_step(d, (donor / d->nx) % 2);
    int c = b->fc_cell[donor];
    for (int kk = 0; kk < n; ++kk, c += stp) {
        for (int el = 0; el < PX_NEL; ++el) sN[el] = px_add(sN[el], px_from_double(px_to_double(b->N[PX_NEL * c + el], PX_N_INV) * phi, PX_N_POW2, &bad), &ov);
        sE = px_add(sE, px_from_double(px_to_double(b->E[c], PX_E_INV) * phi, PX_E_POW2, &bad), &ov);
    }
    for (int el = 0; el < PX_NEL; ++el) b->if_N[PX_NEL * q + el] = sN[el];
    b->if_E[q] = sE;
    const double u0 = b->fc_u[2 * donor], u1 = b->fc_u[2 * donor + 1];
    b->if_K[q] = px_from_double(0.5 * px_fabs(dm) * (u0 * u0 + u1 * u1), PX_E_POW2, &bad);
    if (bad || ov) b->if_st[q] = PX_ERR_FIXED_RANGE;
}
// 3g. Layer pass of one face column: momentum update, kinetic energy, Lorentz and surface-force work, dissipation,
//     departure ejection, and the exact per-cell changes (outflows, inflow and dissipation shares by mass).
//     Energy of the column: K_new - K_old = W_L + W_S + K_in - K_out - K_ejected - Q (Q booked as heat).
PX_HD void px_fc_pass(const px_dims* d, const px_buffers* b, int fc) {
    const int sf = fc / d->nx, f = sf % 2;
    const int lay = b->fc_lay[fc], sc = b->fc_cell[fc];
    for (int el = 0; el < PX_NEL; ++el) b->fc_netN[PX_NEL * fc + el] = px_zero();
    b->fc_un[2 * fc] = 0.0; b->fc_un[2 * fc + 1] = 0.0; b->fc_re[fc] = 0.0; b->fc_s[fc] = 0.0;
    if (lay < 0 || sc < 0) {
        // blocked column (its half of the column is gone): no flow, the column state is kept
        b->fc_nw[fc] = 0;
        for (int k = 0; k < 3; ++k) b->fcV_out[3 * fc + k] = b->fcV[3 * fc + k];
        b->fcK_out[fc] = b->fcK[fc];
        return;
    }
    int st = 0, bad = 0, ov = 0;
    const int n = lay, stp = px_layer_step(d, f);
    const double m = b->fc_m[fc], A = d->dx * d->dz, h = b->fc_h[fc];
    const double ux0 = b->fc_u[2 * fc], uz0 = b->fc_u[2 * fc + 1];
    double M[2] = {m * ux0, m * uz0};
    double phis[5]; int nphi = 0;
    double phi_out = 0.0, m_out = 0.0, m_in = 0.0;
    px_i128 inN[PX_NEL], inE = px_zero(), Kin = px_zero(), Kout = px_zero(), ejK = px_zero(), ejN[PX_NEL], ejE = px_zero();
    for (int el = 0; el < PX_NEL; ++el) { inN[el] = px_zero(); ejN[el] = px_zero(); }
    for (int side = 0; side < 4; ++side) {
        const int q = b->fc_if[4 * fc + side];
        if (q < 0) continue;
        if (b->if_st[q]) px_set_status(&st, b->if_st[q]);
        const int dir = b->if_dir[q], nd = (dir == PX_IF_X || dir == PX_IF_MX || dir == PX_IF_PX) ? 0 : 1, td = 1 - nd;
        const double dm = b->if_dm[q];
        if (b->if_b[q] < 0) {
            const double sgn = (dir == PX_IF_MX || dir == PX_IF_MZ) ? -1.0 : 1.0;
            M[nd] -= sgn * b->if_F[3 * q]; M[td] -= b->if_F[3 * q + 2];
            if (dm > 0.0) {
                phis[nphi++] = dm / m; phi_out += dm / m; m_out += dm;
                ejK = px_add(ejK, b->if_K[q], &ov);
                for (int el = 0; el < PX_NEL; ++el) ejN[el] = px_add(ejN[el], b->if_N[PX_NEL * q + el], &ov);
                ejE = px_add(ejE, b->if_E[q], &ov);
            }
            continue;
        }
        const int isa = b->if_a[q] == fc;
        if (isa) { M[nd] -= b->if_F[3 * q]; M[td] -= b->if_F[3 * q + 2]; }
        else { M[nd] += b->if_F[3 * q + 1]; M[td] += b->if_F[3 * q + 2]; }
        const double out = isa ? dm : -dm;
        if (out > 0.0) { phis[nphi++] = out / m; phi_out += out / m; m_out += out; Kout = px_add(Kout, b->if_K[q], &ov); }
        else if (out < 0.0) {
            m_in -= out; Kin = px_add(Kin, b->if_K[q], &ov);
            for (int el = 0; el < PX_NEL; ++el) inN[el] = px_add(inN[el], b->if_N[PX_NEL * q + el], &ov);
            inE = px_add(inE, b->if_E[q], &ov);
        }
    }
    // departure: liquid lifted beyond the admitted displacement leaves as ejected melt
    double phi_ej = 0.0;
    if (d->eject && n > 0 && b->fc_eta[fc] > d->eta_limit) {
        phi_ej = (b->fc_eta[fc] - d->eta_limit) * d->rho0 * A / m;
        if (phi_ej > 0.5 * d->flow_max_fraction) phi_ej = 0.5 * d->flow_max_fraction;
    }
    if (!(phi_out + phi_ej <= d->flow_max_fraction)) px_set_status(&st, PX_ERR_FLOW);
    // momentum: sources, semi-implicit friction
    double un[2] = {0.0, 0.0};
    const double m_new = m - m_out - phi_ej * m + m_in;
    px_i128 wl = px_zero(), ws = px_zero(), Kej = px_zero(), Knew = px_zero();
    if (n > 0) {
        const double fx = b->fc_f[3 * fc], fz = b->fc_f[3 * fc + 1], fin = b->fc_f[3 * fc + 2], mu = b->fc_mu[fc];
        const double dpx = px_fc_grad(d, b, fc, b->fc_ps, 1, 0, 0), dpz = px_fc_grad(d, b, fc, b->fc_ps, 1, 0, 1);
        const double tsx = b->fc_dsdT[fc] * px_fc_grad(d, b, fc, b->fc_T, 1, 0, 0), tsz = b->fc_dsdT[fc] * px_fc_grad(d, b, fc, b->fc_T, 1, 0, 1);
        const double div = px_fc_grad(d, b, fc, b->fc_u, 2, 0, 0) + px_fc_grad(d, b, fc, b->fc_u, 2, 1, 1);
        const double sx = h * fx + 1.5 * tsx - h * dpx, sz = h * fz + 1.5 * tsz - h * dpz;
        M[0] += d->tau * A * sx; M[1] += d->tau * A * sz;
        const double beta = 3.0 * mu * d->tau / (b->fc_rho[fc] * h * h);
        M[0] /= 1.0 + beta; M[1] /= 1.0 + beta;
        if (m_new > 0.0) { un[0] = M[0] / m_new; un[1] = M[1] / m_new; }
        else { M[0] = 0.0; M[1] = 0.0; }
        const double umx = 0.5 * (ux0 + un[0]), umz = 0.5 * (uz0 + un[1]);
        wl = px_from_double(d->tau * A * (h * (umx * fx + umz * fz) - 0.5 * fin * h * h * div), PX_E_POW2, &bad);
        ws = px_from_double(d->tau * A * (umx * (1.5 * tsx - h * dpx) + umz * (1.5 * tsz - h * dpz) + (tsx * tsx + tsz * tsz) * h / (4.0 * mu)), PX_E_POW2, &bad);
        Kej = px_from_double(phi_ej * 0.5 * m * (ux0 * ux0 + uz0 * uz0), PX_E_POW2, &bad);
        if (m_new > 0.0) Knew = px_from_double(0.5 * m_new * (un[0] * un[0] + un[1] * un[1]), PX_E_POW2, &bad);
        const double sp = px_sqrt(un[0] * un[0] + un[1] * un[1]);
        b->fc_re[fc] = b->fc_rho[fc] * sp * h / mu;
        b->fc_s[fc] = sp + b->fc_a[fc];
    } else { M[0] = 0.0; M[1] = 0.0; }   // inflow onto a column without liquid stops there (its kinetic energy is dissipated)
    // Q = W_L + W_S + K_in - K_out - K_ejected - (K_new - K_old)
    px_i128 Q = px_add(wl, ws, &ov);
    Q = px_add(Q, Kin, &ov); Q = px_sub(Q, Kout, &ov); Q = px_sub(Q, ejK, &ov); Q = px_sub(Q, Kej, &ov);
    Q = px_sub(Q, Knew, &ov); Q = px_add(Q, b->fcK[fc], &ov);
    // per-cell changes: removals (identical words to the interface sums), departure ejection, inflow and heat shares
    const px_i128 H = px_add(inE, Q, &ov);
    const int nw = n > 0 ? n : 1;
    px_i128 usedN[PX_NEL], usedE = px_zero(), netN[PX_NEL], depN[PX_NEL], depE = px_zero();
    for (int el = 0; el < PX_NEL; ++el) { usedN[el] = px_zero(); netN[el] = px_zero(); depN[el] = px_zero(); }
    for (int kk = 1; kk <= nw; ++kk) {
        const int k2 = kk == nw ? 0 : kk;           // the surface cell last (it takes the remainders)
        const int c = sc + k2 * stp;
        px_i128 dNc[PX_NEL], dEc = px_zero();
        for (int el = 0; el < PX_NEL; ++el) dNc[el] = px_zero();
        if (n > 0) {
            for (int el = 0; el < PX_NEL; ++el) {
                const double a = px_to_double(b->N[PX_NEL * c + el], PX_N_INV);
                for (int r = 0; r < nphi; ++r) dNc[el] = px_sub(dNc[el], px_from_double(a * phis[r], PX_N_POW2, &bad), &ov);
                if (phi_ej > 0.0) {
                    const px_i128 x = px_from_double(a * phi_ej, PX_N_POW2, &bad);
                    dNc[el] = px_sub(dNc[el], x, &ov); depN[el] = px_add(depN[el], x, &ov);
                }
            }
            const double e = px_to_double(b->E[c], PX_E_INV);
            for (int r = 0; r < nphi; ++r) dEc = px_sub(dEc, px_from_double(e * phis[r], PX_E_POW2, &bad), &ov);
            if (phi_ej > 0.0) { const px_i128 x = px_from_double(e * phi_ej, PX_E_POW2, &bad); dEc = px_sub(dEc, x, &ov); depE = px_add(depE, x, &ov); }
        }
        if (k2 != 0) {
            const double share = b->m[c] / m;
            for (int el = 0; el < PX_NEL; ++el) {
                const px_i128 x = px_from_double(px_to_double(inN[el], PX_N_INV) * share, PX_N_POW2, &bad);
                dNc[el] = px_add(dNc[el], x, &ov); usedN[el] = px_add(usedN[el], x, &ov);
            }
            const px_i128 x = px_from_double(px_to_double(H, PX_E_INV) * share, PX_E_POW2, &bad);
            dEc = px_add(dEc, x, &ov); usedE = px_add(usedE, x, &ov);
        } else {
            for (int el = 0; el < PX_NEL; ++el) dNc[el] = px_add(dNc[el], px_sub(inN[el], usedN[el], &ov), &ov);
            dEc = px_add(dEc, px_sub(H, usedE, &ov), &ov);
        }
        for (int el = 0; el < PX_NEL; ++el) { b->cflow_N[PX_NEL * c + el] = dNc[el]; netN[el] = px_add(netN[el], dNc[el], &ov); }
        b->cflow_E[c] = dEc;
    }
    b->fc_nw[fc] = nw;
    for (int el = 0; el < PX_NEL; ++el) b->fc_netN[PX_NEL * fc + el] = netN[el];
    b->fcV_out[3 * fc] = M[0]; b->fcV_out[3 * fc + 1] = M[1]; b->fcV_out[3 * fc + 2] = (n > 0 && m_new > 0.0) ? m_new : 0.0;
    b->fcK_out[fc] = Knew;
    b->fc_un[2 * fc] = un[0]; b->fc_un[2 * fc + 1] = un[1];
    // step accumulators
    px_i128* acc = b->acc_fc + PX_FC_ACC * fc;
    acc[FA_WL] = px_add(acc[FA_WL], wl, &ov);
    acc[FA_WS] = px_add(acc[FA_WS], ws, &ov);
    acc[FA_Q] = px_add(acc[FA_Q], Q, &ov);
    for (int el = 0; el < PX_NEL; ++el) acc[FA_EJN + el] = px_add(acc[FA_EJN + el], px_add(ejN[el], depN[el], &ov), &ov);
    const px_i128 kej = px_add(ejK, Kej, &ov);
    acc[FA_EJK] = px_add(acc[FA_EJK], kej, &ov);
    acc[FA_EJE] = px_add(acc[FA_EJE], px_add(px_add(ejE, depE, &ov), kej, &ov), &ov);
    if (bad || ov) px_set_status(&st, PX_ERR_FIXED_RANGE);
    if (st) px_set_status(&b->status[sc], st);
}

// 4a. v1.2: surface physics of one element (cell c, exposed face f) for the current sub-step: prescribed-temperature
//     exchange, loads, deposits, releases, evaporation (HKL with Knudsen layer), radiation. One thread per element
//     (full warps for the transcendental evaluations); the cell step adds el_dE / el_dN. Status codes keep the
//     overwrite semantics of the evaporation functions; range flags are carried to the cell.
PX_HD void px_surface_flux(const px_dims* d, const px_buffers* b, int e) {
    const px_law* L = b->law;
    const int c = b->el_cell[e], f = b->el_face[e];
    int st = 0, ov = 0, bad = 0;
    px_i128 dE = px_zero(), dN[PX_NEL], dNf[PX_NEL];
    for (int el = 0; el < PX_NEL; ++el) { dN[el] = px_zero(); dNf[el] = px_zero(); }
    if (!(b->flags[c] & 1)) {
        const double Tc = b->T[c], kc = b->k[c];
        double atoms[PX_NEL], x[PX_NEL];
        for (int el = 0; el < PX_NEL; ++el) atoms[el] = px_to_double(b->N[PX_NEL * c + el], PX_N_INV);
        px_fractions(atoms, x);
        const int nb = px_neighbor(d, c, f);
            if (nb < 0 && d->bc_kind[f] == 1) {
                const double F = kc * px_face_area(d, b, c, f) / px_face_half(d, b, c, f) * (d->bc_T[f] - Tc) * d->tau;
                const px_i128 q = px_from_double(F, PX_E_POW2, &bad);
                dE = px_add(dE, q, &ov);
                b->acc_bnd_E[e] = px_add(b->acc_bnd_E[e], q, &ov);
            }
            dE = px_add(dE, px_part(b->in_load_E[e], d->n_sub, d->sub, &ov), &ov);
            dE = px_add(dE, px_part(b->in_dep_E[e], d->n_sub, d->sub, &ov), &ov);
            for (int el = 0; el < PX_NEL; ++el) {
                const px_i128 p = px_part(b->in_dep_N[PX_NEL * e + el], d->n_sub, d->sub, &ov);
                dN[el] = px_add(dN[el], p, &ov); dNf[el] = px_add(dNf[el], p, &ov);
            }
            if (d->evap_mode == 2) {
                dE = px_add(dE, px_part(b->in_rel_E[e], d->n_sub, d->sub, &ov), &ov);
                for (int el = 0; el < PX_NEL; ++el) {
                    const px_i128 p = px_part(b->in_rel_N[PX_NEL * e + el], d->n_sub, d->sub, &ov);
                    dN[el] = px_add(dN[el], p, &ov); dNf[el] = px_add(dNf[el], p, &ov);
                }
            } else if (d->evap_mode == 1 && d->evap_face[f] && L->vapor_enabled) {
                const double A = b->el_area[e];
                const double fnet = px_flux_fraction(L, Tc, x, b->patch_pb[b->el_patch[e]], &st);
                for (int el = 0; el < PX_NEL; ++el) {
                    const double D = px_flux_f(L, el, Tc, x, fnet, &st) * A * d->tau;
                    const px_i128 qn = px_from_double(D, PX_N_POW2, &bad);
                    if (px_is_zero(qn)) continue;
                    const double a = px_to_double(qn, PX_N_INV);
                    const px_i128 qr = px_from_double(a * L->eps0[el], PX_E_POW2, &bad);
                    const px_i128 qk = px_from_double(a * (2.0 * PX_KB * Tc), PX_E_POW2, &bad);
                    dN[el] = px_sub(dN[el], qn, &ov); dNf[el] = px_sub(dNf[el], qn, &ov);
                    dE = px_sub(dE, px_add(qr, qk, &ov), &ov);
                    b->acc_ev_N[PX_NEL * e + el] = px_add(b->acc_ev_N[PX_NEL * e + el], qn, &ov);
                    b->acc_ev_Eref[e] = px_add(b->acc_ev_Eref[e], qr, &ov);
                    b->acc_ev_Ekin[e] = px_add(b->acc_ev_Ekin[e], qk, &ov);
                    if (L->gas) {
                        const px_i128 qi = px_from_double(a * px_eint(L, el, Tc), PX_E_POW2, &bad);
                        dE = px_sub(dE, qi, &ov);
                        b->acc_ev_Eint[e] = px_add(b->acc_ev_Eint[e], qi, &ov);
                    }
                }
            }
            if (L->radiation) {
                const double T4 = Tc * Tc * Tc * Tc, E4 = L->T_env * L->T_env * L->T_env * L->T_env;
                const px_i128 q = px_from_double(L->emissivity * PX_SIGMA_SB * (T4 - E4) * b->el_area[e] * d->tau, PX_E_POW2, &bad);
                dE = px_sub(dE, q, &ov);
                b->acc_rad[e] = px_add(b->acc_rad[e], q, &ov);
            }
    }
    b->el_dE[e] = dE;
    for (int el = 0; el < PX_NEL; ++el) b->el_dN[PX_NEL * e + el] = dN[el];
    b->el_st[e] = st; b->el_bad[e] = (bad || ov) ? 1 : 0;
    (void)dNf;
}

// 4. One explicit sub-step of one cell: conduction, prescribed-temperature exchange, Joule, surface loads,
//    deposits, releases, evaporation, radiation and melt-layer transfers.
PX_HD void px_step_cell(const px_dims* d, const px_buffers* b, int c) {
    const px_law* L = b->law;
    for (int el = 0; el < PX_NEL; ++el) b->N_out[PX_NEL * c + el] = b->N[PX_NEL * c + el];
    b->E_out[c] = b->E[c];
    if (b->flags_out) b->flags_out[c] = b->flags[c];
    if (b->flags[c] & 1) return;
    int st = b->status[c], ov = 0, bad = 0;
    px_i128 dE = px_zero();
    px_i128 dN[PX_NEL];
    for (int el = 0; el < PX_NEL; ++el) dN[el] = px_zero();
    const double Tc = b->T[c], kc = b->k[c];
    const int s = c / (d->nx * d->ny), i = c % d->nx;
    for (int f = 0; f < 4; ++f) {
        const int nb = px_neighbor(d, c, f);
        if (nb >= 0 && !(b->flags[nb] & 1)) {
            const double g = px_face_g(d, b, c, f);
            if (f == 0 || f == 2) {
                const double F = px_face_energy(b->T[nb], b->k[nb], Tc, kc, g, d->tau);
                dE = px_add(dE, px_from_double(F, PX_E_POW2, &bad), &ov);
            } else {
                const double F = px_face_energy(Tc, kc, b->T[nb], b->k[nb], g, d->tau);
                dE = px_sub(dE, px_from_double(F, PX_E_POW2, &bad), &ov);
            }
            continue;
        }
        const int e = b->el_of[4 * c + f];
        if (e < 0) continue;
        // surface physics of this face (boundary exchange, loads, deposits, releases, evaporation, radiation),
        // computed per element by px_surface_flux; the sums are exact, so the order does not change the words
        px_i128 dNf[PX_NEL];                 // atoms gained through this face (displacement bookkeeping)
        dE = px_add(dE, b->el_dE[e], &ov);
        for (int el = 0; el < PX_NEL; ++el) { dNf[el] = b->el_dN[PX_NEL * e + el]; dN[el] = px_add(dN[el], dNf[el], &ov); }
        if (b->el_st[e]) st = b->el_st[e];
        if (b->el_bad[e]) bad = 1;
        if (f < 2) {
            const int fc = px_fc_index(d, s, f, i);
            if (d->flow && d->flow_model == 1) {
                // melt-layer transfers through this face column's two interfaces (+ toward +x)
                const int base = (s * 2 + f) * (d->nx - 1);
                if (i > 0) {
                    const int ic = base + i - 1;
                    for (int el = 0; el < PX_NEL; ++el) { dN[el] = px_add(dN[el], b->ic_N[PX_NEL * ic + el], &ov); dNf[el] = px_add(dNf[el], b->ic_N[PX_NEL * ic + el], &ov); }
                    dE = px_add(dE, b->ic_E[ic], &ov);
                    if (d->book_work) { px_u32 r = 0; const px_i128 h = px_divmod_u32(b->ic_W[ic], 2u, &r, &ov); dE = px_add(dE, h, &ov); }
                }
                if (i + 1 < d->nx) {
                    const int ic = base + i;
                    for (int el = 0; el < PX_NEL; ++el) { dN[el] = px_sub(dN[el], b->ic_N[PX_NEL * ic + el], &ov); dNf[el] = px_sub(dNf[el], b->ic_N[PX_NEL * ic + el], &ov); }
                    dE = px_sub(dE, b->ic_E[ic], &ov);
                    if (d->book_work) {
                        // the lower-x cell takes the larger half (exact split of the work word)
                        px_u32 r = 0; px_i128 h = px_divmod_u32(b->ic_W[ic], 2u, &r, &ov);
                        if (r) { px_i128 one; one.lo = 1; one.hi = 0; h = px_add(h, one, &ov); }
                        dE = px_add(dE, h, &ov);
                    }
                }
            }
            if (d->flow && d->flow_model == 2)
                for (int el = 0; el < PX_NEL; ++el) dNf[el] = px_add(dNf[el], b->fc_netN[PX_NEL * fc + el], &ov);
            for (int el = 0; el < PX_NEL; ++el) b->fcm_out[PX_NEL * fc + el] = px_add(b->fcm[PX_NEL * fc + el], dNf[el], &ov);
        }
    }
    dE = px_add(dE, px_part(b->joule_total[c], d->n_sub, d->sub, &ov), &ov);
    if (d->flow && d->flow_model == 2) {
        // depth-averaged layer: transfers, inflows, ejection and booked dissipation of this cell (written by the layer
        // pass of the face column whose half of the column holds the cell)
        const int j = (c / d->nx) % d->ny, f = j < d->ny / 2 ? 0 : 1;
        const int fc = px_fc_index(d, s, f, i), sc = b->fc_cell[fc], nw = b->fc_nw[fc];
        if (sc >= 0 && nw > 0) {
            const int jsc = (sc / d->nx) % d->ny, depth = f == 0 ? j - jsc : jsc - j;
            if (depth >= 0 && depth < nw) {
                dE = px_add(dE, b->cflow_E[c], &ov);
                for (int el = 0; el < PX_NEL; ++el) dN[el] = px_add(dN[el], b->cflow_N[PX_NEL * c + el], &ov);
            }
        }
    }
    const px_i128 En = px_add(b->E[c], dE, &ov);
    double m_new = 0.0;
    for (int el = 0; el < PX_NEL; ++el) {
        const px_i128 Nn = px_add(b->N[PX_NEL * c + el], dN[el], &ov);
        if (px_is_neg(Nn)) px_set_status(&st, PX_ERR_NEGATIVE_INV);
        b->N_out[PX_NEL * c + el] = Nn;
        m_new += px_to_double(Nn, PX_N_INV) * L->atom_mass[el];
    }
    if (!(m_new >= d->phi_floor * b->mref[c])) px_set_status(&st, PX_ERR_DEPLETED);
    if (bad || ov) px_set_status(&st, PX_ERR_FIXED_RANGE);
    b->E_out[c] = En;
    b->status[c] = st;
}

// 5. Per-slice exact sums of the state and the first failing cell; merge candidates (end of step). Two levels:
//    px_row_reduce per (slice, row) over the nx cells of the row, px_slice_reduce over the rows of the slice in row
//    order (exact integer sums; the first failing cell and the first merge candidate in cell order; the largest T).
#define PX_ROW_SUM (1 + PX_NEL)
PX_HD void px_row_reduce(const px_dims* d, const px_buffers* b, int r) {
    int ov = 0;
    px_i128 sE = px_zero(), sN[PX_NEL];
    for (int el = 0; el < PX_NEL; ++el) sN[el] = px_zero();
    int first = -1, code = 0, cand = 0, first_cand = -1;
    double maxT = 0.0;
    for (int c = r * d->nx; c < (r + 1) * d->nx; ++c) {
        if (first < 0 && b->status[c] != PX_OK) { first = c; code = b->status[c]; }
        if (b->flags[c] & 1) continue;
        if (b->T[c] > maxT) maxT = b->T[c];
        sE = px_add(sE, b->E[c], &ov);
        for (int el = 0; el < PX_NEL; ++el) sN[el] = px_add(sN[el], b->N[PX_NEL * c + el], &ov);
        bool exposed = false;
        for (int f = 0; f < 4; ++f) if (b->el_of[4 * c + f] >= 0) exposed = true;
        if (exposed) {
            double atoms[PX_NEL];
            const double m = px_cell_mass(b->law, b->N, c, atoms);
            if (m < d->phi_merge * b->mref[c]) { ++cand; if (first_cand < 0) first_cand = c; }
        }
    }
    px_i128* o = b->row_sum + PX_ROW_SUM * r;
    o[0] = sE; for (int el = 0; el < PX_NEL; ++el) o[1 + el] = sN[el];
    int* q = b->row_info + 5 * r;
    q[0] = first; q[1] = code; q[2] = cand; q[3] = first_cand; q[4] = ov;
    b->row_maxT[r] = maxT;
}
PX_HD void px_slice_reduce(const px_dims* d, const px_buffers* b, int s) {
    int ov = 0;
    px_i128 sE = px_zero(), sN[PX_NEL];
    for (int el = 0; el < PX_NEL; ++el) sN[el] = px_zero();
    int first = -1, code = 0, cand = 0, first_cand = -1;
    double maxT = 0.0;
    const int per = d->nx * d->ny;
    for (int r = s * d->ny; r < (s + 1) * d->ny; ++r) {
        const px_i128* o = b->row_sum + PX_ROW_SUM * r;
        const int* q = b->row_info + 5 * r;
        if (first < 0 && q[0] >= 0) { first = q[0]; code = q[1]; }
        cand += q[2];
        if (first_cand < 0 && q[3] >= 0) first_cand = q[3];
        if (q[4]) ov = 1;
        if (b->row_maxT[r] > maxT) maxT = b->row_maxT[r];
        sE = px_add(sE, o[0], &ov);
        for (int el = 0; el < PX_NEL; ++el) sN[el] = px_add(sN[el], o[1 + el], &ov);
    }
    // v1.2: layer kinetic energy of the slice's face columns; accumulated Joule energy of its EM cells; diagnostics
    px_i128 sK = px_zero(), sJ = px_zero();
    for (int fc = 2 * s * d->nx; fc < (2 * s + 2) * d->nx; ++fc) sK = px_add(sK, b->fcK[fc], &ov);
    const int per_em = d->nx_em * d->ny_em;
    double lag = 0.0, heat = 0.0, accbad = 0.0, smax = 0.0;
    for (int e = s * per_em; e < (s + 1) * per_em; ++e) {
        if (d->use_acc) { sJ = px_add(sJ, b->acc_E[e], &ov); if (b->acc_bad[e]) accbad = 1.0; }
        if (d->diag) {
            const double s0 = b->em_diag[2 * e + 1];
            if (s0 > 0.0) { const double r = px_fabs(b->sigma_em[e] / s0 - 1.0); if (r > lag) lag = r; }
            if (b->em_diag[2 * e] > heat) heat = b->em_diag[2 * e];
        }
    }
    if (d->flow && d->flow_model == 2) for (int fc = 2 * s * d->nx; fc < (2 * s + 2) * d->nx; ++fc) if (b->fc_s[fc] > smax) smax = b->fc_s[fc];
    px_i128* o = b->slice_sum + PX_SLICE_SUM * s;
    o[0] = sE; for (int el = 0; el < PX_NEL; ++el) o[1 + el] = sN[el];
    o[SL_K] = sK; o[SL_J] = sJ;
    b->slice_diag[4 * s] = lag; b->slice_diag[4 * s + 1] = heat; b->slice_diag[4 * s + 2] = accbad; b->slice_diag[4 * s + 3] = smax;
    b->slice_status[4 * s + 0] = ov ? (first >= 0 ? first : s * per) : first;
    b->slice_status[4 * s + 1] = ov ? PX_ERR_FIXED_RANGE : code;
    b->slice_status[4 * s + 2] = cand;
    b->slice_status[4 * s + 3] = first_cand;
    b->slice_maxT[s] = maxT;
}

// 6. Expected atoms per element over dt at the current state (quantized-release quotes).
PX_HD void px_quote_element(const px_dims* d, const px_buffers* b, int e) {
    const int c = b->el_cell[e], f = b->el_face[e];
    int st = 0;
    double atoms[PX_NEL], x[PX_NEL];
    for (int el = 0; el < PX_NEL; ++el) atoms[el] = px_to_double(b->N[PX_NEL * c + el], PX_N_INV);
    px_fractions(atoms, x);
    const double fnet = px_flux_fraction(b->law, b->T[c], x, b->patch_pb[b->el_patch[e]], &st);
    for (int el = 0; el < PX_NEL; ++el) {
        double D = 0.0;
        if (d->evap_face[f] && b->law->vapor_enabled) D = px_flux_f(b->law, el, b->T[c], x, fnet, &st) * b->el_area[e] * d->dt;
        b->el_quote[PX_NEL * e + el] = D;
    }
}

// 7. Split the quantized release of one patch over its elements in proportion to the element quotes
//    (remainder to the first element with the largest quote): atoms exactly, reference energy eps0 per atom,
//    internal energy E_int(T_cell) per atom (gas = 1), actual kinetic energy split by the same weights.
PX_HD void px_release_split(const px_dims* d, const px_buffers* b, int p) {
    (void)d;
    const int e0 = b->patch_begin[p], e1 = b->patch_begin[p + 1];
    int bad = 0, ov = 0;
    for (int e = e0; e < e1; ++e) { b->in_rel_E[e] = px_zero(); for (int el = 0; el < PX_NEL; ++el) b->in_rel_N[PX_NEL * e + el] = px_zero(); }
    double eint_total = 0.0;
    for (int el = 0; el < PX_NEL; ++el) {
        const px_i128 R = b->rel_total_N[PX_NEL * p + el];
        const double K = b->rel_kin_J[PX_NEL * p + el];
        if (px_is_zero(R) && K == 0.0) continue;
        double W = 0.0, wmax = -1.0; int emax = -1;
        for (int e = e0; e < e1; ++e) { const double w = b->el_quote[PX_NEL * e + el]; W += w; if (w > wmax) { wmax = w; emax = e; } }
        if (!(W > 0.0)) { b->status[b->el_cell[e0]] = PX_ERR_NEGATIVE_INV; continue; }
        const double Rd = px_to_double(R, PX_N_INV);
        const px_i128 Kq = px_from_double(K, PX_E_POW2, &bad);
        px_i128 usedN = px_zero(), usedK = px_zero();
        for (int e = e0; e <= e1; ++e) {
            if (e == emax) continue;
            const bool last = e == e1;
            const int ee = last ? emax : e;
            px_i128 qn, qk;
            if (!last) {
                const double share = b->el_quote[PX_NEL * e + el] / W;
                qn = px_from_double(Rd * share, PX_N_POW2, &bad);
                qk = px_from_double(K * share, PX_E_POW2, &bad);
                usedN = px_add(usedN, qn, &ov); usedK = px_add(usedK, qk, &ov);
            } else {
                qn = px_sub(R, usedN, &ov); qk = px_sub(Kq, usedK, &ov);
                if (px_is_neg(qn) || px_is_neg(qk)) bad = 1;
            }
            const double an = px_to_double(qn, PX_N_INV);
            px_i128 qr = px_from_double(an * b->law->eps0[el], PX_E_POW2, &bad);
            if (b->law->gas) {
                const double ei = an * px_eint(b->law, el, b->T[b->el_cell[ee]]);
                qr = px_add(qr, px_from_double(ei, PX_E_POW2, &bad), &ov);
                eint_total += ei;
            }
            b->in_rel_N[PX_NEL * ee + el] = px_neg(qn, &ov);
            b->in_rel_E[ee] = px_sub(b->in_rel_E[ee], px_add(qr, qk, &ov), &ov);
        }
    }
    b->patch_out[PX_PATCH_OUT * p + XO_EINT] = eint_total;
    if (bad || ov) b->status[b->el_cell[e0]] = PX_ERR_FIXED_RANGE;
}

// 8. Patch report at the current state and exact sums of the step accumulators of its elements. px_patch_sums: the
//    step part (accumulator sums, displacement, melt speed and Reynolds number, layer means); px_patch_report: the
//    sums plus the surface state (temperatures, melt depth, evaporation flux, composition, recoil).
PX_HD void px_patch_sums(const px_dims* d, const px_buffers* b, int p) {
    const int e0 = b->patch_begin[p], e1 = b->patch_begin[p + 1];
    double A = 0, seta = 0, meta = 0, mneg = 0, sspeed = 0, mre = 0, sux = 0, suz = 0, sh = 0;
    int ov = 0, nspeed = 0, nfc = 0;
    px_i128 acc[PX_PATCH_SUM];
    for (int q = 0; q < PX_PATCH_SUM; ++q) acc[q] = px_zero();
    const px_law* L = b->law;
    for (int e = e0; e < e1; ++e) {
        const int c = b->el_cell[e], f = b->el_face[e];
        const double a = b->el_area[e];
        A += a;
        if (f < 2) {
            const int s = c / (d->nx * d->ny), i = c % d->nx;
            const int fc = px_fc_index(d, s, f, i);
            const double eta = b->fc_eta[fc];
            seta += a * eta; if (px_fabs(eta) > meta) meta = px_fabs(eta);
            if (-eta > mneg) mneg = -eta;
            if (d->flow && d->flow_model != 2 && i + 1 < d->nx) {
                const int ic = (s * 2 + f) * (d->nx - 1) + i;
                const double dl = 0.5 * (b->fc_delta[fc] + b->fc_delta[fc + 1]);
                if (dl > 0.0) {
                    sspeed += px_fabs(b->ic_q[ic]) / dl; ++nspeed;
                    double sg, dsg, mu; px_surface_props(L, b->fc_T[fc], &sg, &dsg, &mu);
                    const double rho = L->kind == 0 ? px_density(L, b->fc_T[fc]) : L->c_rho;
                    const double re = rho * px_fabs(b->ic_q[ic]) / mu;
                    if (re > mre) mre = re;
                }
                if (d->book_work) {
                    acc[XS_WL] = px_add(acc[XS_WL], b->acc_ic[2 * ic], &ov);
                    acc[XS_WS] = px_add(acc[XS_WS], b->acc_ic[2 * ic + 1], &ov);
                    acc[XS_Q] = px_add(acc[XS_Q], px_add(b->acc_ic[2 * ic], b->acc_ic[2 * ic + 1], &ov), &ov);
                }
            }
            if (d->flow && d->flow_model == 2) {
                const double ux = b->fc_u[2 * fc], uz = b->fc_u[2 * fc + 1], h = b->fc_h[fc];
                ++nfc; sux += ux; suz += uz; sh += h;
                if (h > 0.0) {
                    const double sp = px_sqrt(ux * ux + uz * uz);
                    sspeed += sp; ++nspeed;
                    const double re = b->fc_rho[fc] * sp * h / b->fc_mu[fc];
                    if (re > mre) mre = re;
                }
                const px_i128* fa = b->acc_fc + PX_FC_ACC * fc;
                acc[XS_WL] = px_add(acc[XS_WL], fa[FA_WL], &ov);
                acc[XS_WS] = px_add(acc[XS_WS], fa[FA_WS], &ov);
                acc[XS_Q] = px_add(acc[XS_Q], fa[FA_Q], &ov);
                for (int el = 0; el < PX_NEL; ++el) acc[XS_EJN + el] = px_add(acc[XS_EJN + el], fa[FA_EJN + el], &ov);
                acc[XS_EJE] = px_add(acc[XS_EJE], fa[FA_EJE], &ov);
                acc[XS_EJK] = px_add(acc[XS_EJK], fa[FA_EJK], &ov);
            }
        }
        for (int el = 0; el < PX_NEL; ++el) acc[XS_EVN + el] = px_add(acc[XS_EVN + el], b->acc_ev_N[PX_NEL * e + el], &ov);
        acc[XS_EREF] = px_add(acc[XS_EREF], b->acc_ev_Eref[e], &ov);
        acc[XS_EKIN] = px_add(acc[XS_EKIN], b->acc_ev_Ekin[e], &ov);
        acc[XS_EINT] = px_add(acc[XS_EINT], b->acc_ev_Eint[e], &ov);
        acc[XS_BND] = px_add(acc[XS_BND], b->acc_bnd_E[e], &ov);
        acc[XS_LOAD] = px_add(acc[XS_LOAD], b->in_load_E[e], &ov);
        acc[XS_DEP] = px_add(acc[XS_DEP], b->in_dep_E[e], &ov);
        acc[XS_REL] = px_add(acc[XS_REL], b->in_rel_E[e], &ov);
        acc[XS_RAD] = px_add(acc[XS_RAD], b->acc_rad[e], &ov);
    }
    double* o = b->patch_out + PX_PATCH_OUT * p;
    o[XO_ETA_MEAN] = A > 0 ? seta / A : 0; o[XO_ETA_MAX] = meta; o[XO_ETA_NEG] = mneg;
    o[XO_SPEED] = nspeed > 0 ? sspeed / nspeed : 0;
    o[XO_RE] = mre;
    o[XO_UX] = nfc > 0 ? sux / nfc : 0; o[XO_UZ] = nfc > 0 ? suz / nfc : 0; o[XO_H] = nfc > 0 ? sh / nfc : 0;
    for (int q = 0; q < PX_PATCH_SUM; ++q) b->patch_sum[PX_PATCH_SUM * p + q] = acc[q];
    (void)ov;
}
PX_HD void px_patch_report(const px_dims* d, const px_buffers* b, int p) {
    px_patch_sums(d, b, p);
    const int e0 = b->patch_begin[p], e1 = b->patch_begin[p + 1];
    double A = 0, sT = 0, mT = 0, sf = 0, sdepth = 0, srho = 0, srec = 0;
    double sD[PX_NEL], sj[PX_NEL], sx[PX_NEL];
    for (int el = 0; el < PX_NEL; ++el) { sD[el] = 0; sj[el] = 0; sx[el] = 0; }
    int st = 0;
    const px_law* L = b->law;
    for (int e = e0; e < e1; ++e) {
        const int c = b->el_cell[e], f = b->el_face[e];
        const double a = b->el_area[e], T = b->T[c];
        A += a; sT += a * T; if (T > mT) mT = T;
        sf += a * px_fl(L, T);
        if (L->kind == 0) srho += a * px_density(L, T);
        double depth = 0.0;
        int cc = c;
        const int inward = f == 0 ? 1 : f == 1 ? 0 : f == 2 ? 3 : 2;
        while (cc >= 0 && !(b->flags[cc] & 1)) {
            const double fl = px_fl(L, b->T[cc]);
            if (fl <= 0.0) break;
            depth += fl * (f < 2 ? b->row_dy[(cc / d->nx) % d->ny] : d->dx);
            cc = px_neighbor(d, cc, inward);
        }
        sdepth += a * depth;
        double atoms[PX_NEL], x[PX_NEL];
        for (int el = 0; el < PX_NEL; ++el) atoms[el] = px_to_double(b->N[PX_NEL * c + el], PX_N_INV);
        px_fractions(atoms, x);
        const double pb = b->patch_pb[p];
        const double fnet = px_flux_fraction(L, T, x, pb, &st);
        for (int el = 0; el < PX_NEL; ++el) {
            const double j = (d->evap_face[f] && L->vapor_enabled) ? px_flux_f(L, el, T, x, fnet, &st) : 0.0;
            sj[el] += a * j;
            sD[el] += j * a * d->dt;
            sx[el] += a * x[el];
        }
        if (L->vapor_enabled && d->evap_face[f]) srec += a * px_recoil(L, T, x, pb, &st);
    }
    double* o = b->patch_out + PX_PATCH_OUT * p;
    o[XO_AREA] = A; o[XO_MEAN_T] = A > 0 ? sT / A : 0; o[XO_MAX_T] = mT; o[XO_FL] = A > 0 ? sf / A : 0; o[XO_MELT] = A > 0 ? sdepth / A : 0;
    for (int el = 0; el < PX_NEL; ++el) { o[XO_D + el] = sD[el]; o[XO_J + el] = A > 0 ? sj[el] / A : 0; o[XO_X + el] = A > 0 ? sx[el] / A : 0; }
    o[XO_RHO] = A > 0 ? srho / A : 0;
    o[XO_RECOIL] = A > 0 ? srec / A : 0;
}

#if defined(__CUDACC_RTC__)
#define PX_KERNEL(name, count, call) extern "C" __global__ void name(px_dims d, px_buffers b) { const int t = blockIdx.x * blockDim.x + threadIdx.x; if (t < (count)) call(&d, &b, t); }
PX_KERNEL(k_thermo, d.n_cells, px_thermo)
PX_KERNEL(k_scatter, d.n_scatter, px_scatter)
PX_KERNEL(k_sigma, d.n_em, px_sigma_em)
PX_KERNEL(k_joule, d.n_em, px_joule_split)
PX_KERNEL(k_fc_prepare, d.n_fc, px_fc_prepare)
PX_KERNEL(k_fc_pressure, d.n_fc, px_fc_pressure)
PX_KERNEL(k_ic_flux, d.n_ic, px_ic_flux)
PX_KERNEL(k_step, d.n_cells, px_step_cell)
PX_KERNEL(k_slice, d.S, px_slice_reduce)
PX_KERNEL(k_rows, d.S * d.ny, px_row_reduce)
PX_KERNEL(k_surface, d.n_el, px_surface_flux)
PX_KERNEL(k_quote, d.n_el, px_quote_element)
PX_KERNEL(k_release, d.n_patches, px_release_split)
PX_KERNEL(k_patch, d.n_patches, px_patch_report)
PX_KERNEL(k_patch_sums, d.n_patches, px_patch_sums)
PX_KERNEL(k_em_start, d.n_em, px_em_start)
PX_KERNEL(k_fc_prepare2, d.n_fc, px_fc_prepare2)
PX_KERNEL(k_fc_pressure2, d.n_fc, px_fc_pressure2)
PX_KERNEL(k_if_flux2, d.n_if, px_if_flux2)
PX_KERNEL(k_fc_pass, d.n_fc, px_fc_pass)
extern "C" __global__ void k_zero_step(px_dims d, px_buffers b, int n) {
    const int t = blockIdx.x * blockDim.x + threadIdx.x; if (t < n) px_zero_step(&d, &b, t);
}
extern "C" __global__ void k_accumulate(int n, const double* P, const double* JB, double dt, px_i128* accE, double* accF, int* bad) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x; if (e < n) px_accumulate_one(e, P, JB, dt, accE, accF, bad);
}
#endif
