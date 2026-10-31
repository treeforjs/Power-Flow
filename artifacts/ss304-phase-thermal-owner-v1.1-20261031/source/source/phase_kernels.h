// phase_kernels.h -- the per-cell, per-EM-cell, per-element, per-face-column, per-interface and per-patch
// arithmetic of one owner step (v1.1). The CPU reference calls these functions in index loops; the NVRTC program
// wraps the same functions in extern "C" kernels (one thread per item). No atomics: every output word is written
// by exactly one item, and every sum runs over a fixed index order (energy and atom sums are exact 128-bit sums).
//
// Layout. Cells c = (s*ny + j)*nx + i. Faces f: 0 LowY, 1 HighY, 2 LowX, 3 HighX. Surface elements e: exposed
// (cell, face) pairs, sorted by patch, cell, face. Face columns fc = (s*2 + f)*nx + i for the two large faces
// f = 0, 1 (the faces that carry the melt-layer flow along x); interfaces ic = (s*2 + f)*(nx-1) + i between
// face columns i and i+1.
#pragma once
#include "phase_law.h"

struct pq_dims {
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
};

#define PQ_SC 8   // scatter values per entry: load_E, dep_N x6, dep_E

struct pq_buffers {
    const pq_law* law;
    // geometry
    const double* row_dy; const int* row_em; const int* em_row_begin; const double* gx; const double* gy;
    const double* vref; const double* mref;
    // state (current and next slot)
    const pq_i128* E; const pq_i128* N; const int* flags; const pq_i128* fcm;   // fcm [n_fc*PQ_NEL]: face mass words
    pq_i128* E_out; pq_i128* N_out; pq_i128* fcm_out;
    // per-cell scratch
    double* T; double* k; double* sig; double* m; int* status;
    pq_i128* joule_total; const double* em_energy;
    const double* em_JB;          // [n_em*6]: J (A/m2) x,y,z then B (T) x,y,z at the EM cell (coupling input)
    // surface elements
    const int* el_of; const int* el_cell; const int* el_face; const double* el_area; const int* el_patch; const int* patch_begin;
    const pq_i128* in_load_E; const pq_i128* in_dep_N; const pq_i128* in_dep_E;
    pq_i128* in_rel_N; pq_i128* in_rel_E;
    pq_i128* acc_ev_N; pq_i128* acc_ev_Eref; pq_i128* acc_ev_Ekin; pq_i128* acc_ev_Eint; pq_i128* acc_bnd_E; pq_i128* acc_rad;
    double* el_quote;             // [n_el*PQ_NEL]
    const double* patch_pb;       // [n_patches] back pressure (Pa)
    // face columns and interfaces (melt flow, displacement)
    int* fc_cell; double* fc_delta; double* fc_T; double* fc_eta; double* fc_p; double* fc_x; // fc_x [n_fc*PQ_NEL]
    pq_i128* ic_N; pq_i128* ic_E; double* ic_q;      // [n_ic*PQ_NEL], [n_ic], [n_ic] volumetric flux (m2/s)
    // patch / slice outputs
    const pq_i128* rel_total_N; const double* rel_kin_J;
    double* patch_out; pq_i128* patch_sum; pq_i128* slice_sum; int* slice_status; double* slice_maxT; double* sigma_em;
    const int* sc_el; const pq_i128* sc_val;
    pq_i128* in_load_E_w; pq_i128* in_dep_N_w; pq_i128* in_dep_E_w;
};

// patch_out layout
#define PO_AREA 0
#define PO_MEAN_T 1
#define PO_MAX_T 2
#define PO_FL 3
#define PO_MELT 4
#define PO_D 5            // +el, expected atoms over dt
#define PO_J 11           // +el, flux atoms/m2/s
#define PO_RHO 17         // mean surface density
#define PO_ETA_MEAN 18    // mean displacement (m, outward)
#define PO_ETA_MAX 19     // max |displacement|
#define PO_RECOIL 20      // mean recoil pressure (Pa)
#define PO_SPEED 21       // mean |melt speed| at the patch interfaces (m/s)
#define PO_X 22           // +el, mean surface atom fraction
#define PO_EINT 28        // released internal energy (quantized mode, J)
#define PO_RE 29          // largest melt-layer Reynolds number rho |q| / mu at the patch interfaces
#define PQ_PATCH_OUT 30
// patch_sum layout
#define PS_EVN 0          // +el
#define PS_EREF 6
#define PS_EKIN 7
#define PS_EINT 8
#define PS_BND 9
#define PS_LOAD 10
#define PS_DEP 11
#define PS_REL 12
#define PS_RAD 13
#define PQ_PATCH_SUM 14
#define PQ_SLICE_SUM (1 + PQ_NEL)

PQ_HD int pq_neighbor(const pq_dims* d, int c, int f) {
    const int i = c % d->nx, j = (c / d->nx) % d->ny;
    if (f == 0) return j > 0 ? c - d->nx : -1;
    if (f == 1) return j + 1 < d->ny ? c + d->nx : -1;
    if (f == 2) return i > 0 ? c - 1 : -1;
    return i + 1 < d->nx ? c + 1 : -1;
}
PQ_HD double pq_face_g(const pq_dims* d, const pq_buffers* b, int c, int f) {
    const int j = (c / d->nx) % d->ny;
    if (f == 0) return b->gy[j - 1];
    if (f == 1) return b->gy[j];
    return b->gx[j];
}
PQ_HD double pq_face_area(const pq_dims* d, const pq_buffers* b, int c, int f) {
    const int j = (c / d->nx) % d->ny;
    return f < 2 ? d->dx * d->dz : b->row_dy[j] * d->dz;
}
PQ_HD double pq_face_half(const pq_dims* d, const pq_buffers* b, int c, int f) {
    const int j = (c / d->nx) % d->ny;
    return f < 2 ? 0.5 * b->row_dy[j] : 0.5 * d->dx;
}
PQ_HD double pq_face_energy(double Ta, double ka, double Tb, double kb, double g, double tau) {
    const double kh = (ka == 0.0 || kb == 0.0) ? 0.0 : 2.0 * ka * kb / (ka + kb);
    return kh * g * (Ta - Tb) * tau;
}
PQ_HD pq_i128 pq_part(pq_i128 total, int n, int s, int* ov) {
    pq_u32 r = 0;
    pq_i128 q = pq_divmod_u32(total, (pq_u32)n, &r, ov);
    if ((pq_u32)s < r) { pq_i128 one; one.lo = 1; one.hi = 0; q = pq_add(q, one, ov); }
    return q;
}
PQ_HD void pq_set_status(int* st, int code) { if (*st == PQ_OK) *st = code; }
PQ_HD double pq_cell_mass(const pq_law* L, const pq_i128* N, int c, double* atoms) {
    double m = 0.0;
    for (int el = 0; el < PQ_NEL; ++el) { atoms[el] = pq_to_double(N[PQ_NEL * c + el], PQ_N_INV); m += atoms[el] * L->atom_mass[el]; }
    return m;
}
PQ_HD void pq_fractions(const double* atoms, double* x) {
    double tot = 0.0;
    for (int el = 0; el < PQ_NEL; ++el) tot += atoms[el];
    for (int el = 0; el < PQ_NEL; ++el) x[el] = atoms[el] / tot;
}
PQ_HD int pq_fc_index(const pq_dims* d, int s, int f, int i) { return (s * 2 + f) * d->nx + i; }

// 1. Thermodynamic state of one cell from (E, N).
PQ_HD void pq_thermo(const pq_dims* d, const pq_buffers* b, int c) {
    (void)d;
    int st = PQ_OK;
    if (b->flags[c] & 1) { b->T[c] = 0; b->k[c] = 0; b->sig[c] = 0; b->m[c] = 0; return; }
    double atoms[PQ_NEL];
    const double m = pq_cell_mass(b->law, b->N, c, atoms);
    const double h = pq_to_double(b->E[c], PQ_E_INV) / m;
    const double T = pq_T_of_h(b->law, h, &st);
    b->T[c] = T; b->m[c] = m;
    b->k[c] = pq_k(b->law, T);
    b->sig[c] = 1.0 / pq_rhoe(b->law, T);
    if (b->law->T_explosion > 0.0 && T >= b->law->T_explosion) pq_set_status(&st, PQ_ERR_PHASE_EXPLOSION);
    if (b->status[c] == PQ_OK) b->status[c] = st;
}

PQ_HD void pq_scatter(const pq_dims* d, const pq_buffers* b, int q) {
    (void)d;
    const int e = b->sc_el[q];
    b->in_load_E_w[e] = b->sc_val[PQ_SC * q + 0];
    for (int el = 0; el < PQ_NEL; ++el) b->in_dep_N_w[PQ_NEL * e + el] = b->sc_val[PQ_SC * q + 1 + el];
    b->in_dep_E_w[e] = b->sc_val[PQ_SC * q + 1 + PQ_NEL];
}

PQ_HD void pq_sigma_em(const pq_dims* d, const pq_buffers* b, int e) {
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
PQ_HD void pq_joule_split(const pq_dims* d, const pq_buffers* b, int e) {
    const int per_slice = d->nx_em * d->ny_em;
    const int s = e / per_slice, je = (e % per_slice) / d->nx_em, ie = e % d->nx_em;
    const int j0 = b->em_row_begin[je], j1 = b->em_row_begin[je + 1], i0 = ie * d->rx, i1 = i0 + d->rx;
    const double energy = b->em_energy[e];
    int bad = 0, ov = 0;
    double W = 0.0, wmax = -1.0; int cmax = -1;
    for (int j = j0; j < j1; ++j) for (int i = i0; i < i1; ++i) {
        const int c = (s * d->ny + j) * d->nx + i;
        b->joule_total[c] = pq_zero();
        if (b->flags[c] & 1) continue;
        const double w = d->joule_mode == 1 ? b->sig[c] * b->vref[c] : b->vref[c];
        W += w;
        if (w > wmax) { wmax = w; cmax = c; }
    }
    if (energy == 0.0) return;
    if (cmax < 0 || !(W > 0.0)) { b->status[(s * d->ny + j0) * d->nx + i0] = PQ_ERR_NONFINITE; return; }
    const pq_i128 total = pq_from_double(energy, PQ_E_POW2, &bad);
    pq_i128 used = pq_zero();
    for (int j = j0; j < j1; ++j) for (int i = i0; i < i1; ++i) {
        const int c = (s * d->ny + j) * d->nx + i;
        if ((b->flags[c] & 1) || c == cmax) continue;
        const double w = d->joule_mode == 1 ? b->sig[c] * b->vref[c] : b->vref[c];
        const pq_i128 q = pq_from_double(energy * (w / W), PQ_E_POW2, &bad);
        b->joule_total[c] = q;
        used = pq_add(used, q, &ov);
    }
    const pq_i128 rest = pq_sub(total, used, &ov);
    b->joule_total[cmax] = rest;
    if (bad || ov || pq_is_neg(rest)) pq_set_status(&b->status[cmax], PQ_ERR_FIXED_RANGE);
}

// 3a. Face column: surface cell, melt-layer thickness, surface temperature and composition, displacement.
//     eta = (face mass gained) / (rho0 A) + expansion / 2, expansion = sum_c m_c (1/rho(T_c) - 1/rho0) / A.
PQ_HD void pq_fc_prepare(const pq_dims* d, const pq_buffers* b, int fc) {
    const int i = fc % d->nx, sf = fc / d->nx, s = sf / 2, f = sf % 2;
    const pq_law* L = b->law;
    int sc = -1;
    for (int jj = 0; jj < d->ny; ++jj) {
        const int j = f == 0 ? jj : d->ny - 1 - jj;
        const int c = (s * d->ny + j) * d->nx + i;
        if (!(b->flags[c] & 1)) { sc = c; break; }
    }
    b->fc_cell[fc] = sc;
    double delta = 0.0, expansion = 0.0, x[PQ_NEL];
    for (int el = 0; el < PQ_NEL; ++el) x[el] = 0.0;
    double Ts = 0.0;
    if (sc >= 0) {
        Ts = b->T[sc];
        double atoms[PQ_NEL];
        (void)pq_cell_mass(L, b->N, sc, atoms);
        pq_fractions(atoms, x);
        int c = sc;
        const int inward = f == 0 ? 1 : 0;
        while (c >= 0 && !(b->flags[c] & 1)) {
            const double fl = pq_fl(L, b->T[c]);
            if (fl <= 0.0) break;
            // physical thickness of the liquid in this cell from its mass and the liquid density
            const double rho_c = L->kind == 0 ? pq_density(L, b->T[c]) : (L->c_rho > 0.0 ? L->c_rho : d->rho0);
            delta += fl * b->m[c] / (rho_c * d->dx * d->dz);
            c = pq_neighbor(d, c, inward);
        }
        // a face's layer is at most half the foil: a fully molten column is split between its two faces
        if (delta > 0.5 * d->thickness) delta = 0.5 * d->thickness;
        if (L->kind == 0 && d->rho0 > 0.0) {
            for (int j = 0; j < d->ny; ++j) {
                const int cc = (s * d->ny + j) * d->nx + i;
                if (b->flags[cc] & 1) continue;
                expansion += b->m[cc] * (1.0 / pq_density(L, b->T[cc]) - 1.0 / d->rho0);
            }
            expansion /= d->dx * d->dz;
        }
    }
    double gained = 0.0;
    for (int el = 0; el < PQ_NEL; ++el) gained += pq_to_double(b->fcm[PQ_NEL * fc + el], PQ_N_INV) * L->atom_mass[el];
    const double rho = d->rho0 > 0.0 ? d->rho0 : 1.0;
    b->fc_delta[fc] = delta; b->fc_T[fc] = Ts;
    b->fc_eta[fc] = gained / (rho * d->dx * d->dz) + 0.5 * expansion;
    for (int el = 0; el < PQ_NEL; ++el) b->fc_x[PQ_NEL * fc + el] = x[el];
}
// 3b. Face column pressure: recoil (Knudsen layer, back pressure of the patch) - surface tension x curvature
//     - mean normal Lorentz load (n . (J x B)) delta / 2.
PQ_HD double pq_patch_pb(const pq_buffers* b, int c, int f) {
    const int e = b->el_of[4 * c + f];
    return e >= 0 ? b->patch_pb[b->el_patch[e]] : 0.0;
}
PQ_HD void pq_surface_props(const pq_law* L, double T, double* sigma, double* dsdT, double* mu) {
    if (L->kind == 0) { pq_surface_tension(L, T, sigma, dsdT); *mu = pq_viscosity(T); }
    else { *sigma = L->c_sigma; *dsdT = L->c_dsdT; *mu = L->c_mu; }
}
PQ_HD void pq_em_force(const pq_dims* d, const pq_buffers* b, int s, int f, int i, double* fx, double* fn) {
    const int je = f == 0 ? 0 : d->ny_em - 1, ie = i / d->rx;
    const double* v = b->em_JB + 6 * ((s * d->ny_em + je) * d->nx_em + ie);
    const double fxv = v[1] * v[5] - v[2] * v[4];   // (J x B)_x
    const double fyv = v[2] * v[3] - v[0] * v[5];   // (J x B)_y
    *fx = fxv;
    *fn = f == 0 ? -fyv : fyv;                       // outward-normal component
}
PQ_HD void pq_fc_pressure(const pq_dims* d, const pq_buffers* b, int fc) {
    const int i = fc % d->nx, sf = fc / d->nx, s = sf / 2, f = sf % 2;
    const int sc = b->fc_cell[fc];
    if (sc < 0) { b->fc_p[fc] = 0.0; return; }
    int st = 0;
    const double T = b->fc_T[fc];
    const double pr = b->law->vapor_enabled ? pq_recoil(b->law, T, b->fc_x + PQ_NEL * fc, pq_patch_pb(b, sc, f), &st) : 0.0;
    double sigma, dsdT, mu;
    pq_surface_props(b->law, T, &sigma, &dsdT, &mu);
    const double em = b->fc_eta[fc];
    const double el = i > 0 ? b->fc_eta[fc - 1] : em, er = i + 1 < d->nx ? b->fc_eta[fc + 1] : em;
    const double curvature = (er - 2.0 * em + el) / (d->dx * d->dx);
    double fx, fn;
    pq_em_force(d, b, s, f, i, &fx, &fn);
    b->fc_p[fc] = pr - sigma * curvature - fn * b->fc_delta[fc] * 0.5;
}
// 3c. Interface between face columns i and i+1: lubrication flux
//     q = delta^3/(3 mu) (f_x - dp/dx) + delta^2/(2 mu) dsigma/dT dT_s/dx   (m2/s, + toward +x), with delta the
//     mean of the two columns' liquid thicknesses (from mass and liquid density); no flow unless both are molten
//     and the exact transfer (donor composition and specific energy) over tau.
PQ_HD void pq_ic_flux(const pq_dims* d, const pq_buffers* b, int ic) {
    const int i = ic % (d->nx - 1), sf = ic / (d->nx - 1), s = sf / 2, f = sf % 2;
    const int fa = pq_fc_index(d, s, f, i), fb = fa + 1;
    const int ca = b->fc_cell[fa], cb = b->fc_cell[fb];
    const pq_law* L = b->law;
    for (int el = 0; el < PQ_NEL; ++el) b->ic_N[PQ_NEL * ic + el] = pq_zero();
    b->ic_E[ic] = pq_zero(); b->ic_q[ic] = 0.0;
    if (ca < 0 || cb < 0) return;
    if (!(b->fc_delta[fa] > 0.0) || !(b->fc_delta[fb] > 0.0) || pq_fl(L, b->T[ca]) < 1.0 || pq_fl(L, b->T[cb]) < 1.0) return;
    const double delta = 0.5 * (b->fc_delta[fa] + b->fc_delta[fb]);   // layer thickness at the interface
    double sa, dsa, ma, sb, dsb, mb;
    pq_surface_props(L, b->fc_T[fa], &sa, &dsa, &ma);
    pq_surface_props(L, b->fc_T[fb], &sb, &dsb, &mb);
    const double mu = 0.5 * (ma + mb), dsdT = 0.5 * (dsa + dsb);
    double fxa, fna, fxb, fnb;
    pq_em_force(d, b, s, f, i, &fxa, &fna);
    pq_em_force(d, b, s, f, i + 1, &fxb, &fnb);
    const double fx = 0.5 * (fxa + fxb);
    const double dpdx = (b->fc_p[fb] - b->fc_p[fa]) / d->dx;
    const double tau_m = dsdT * (b->fc_T[fb] - b->fc_T[fa]) / d->dx;
    const double q = delta * delta * delta / (3.0 * mu) * (fx - dpdx) + delta * delta / (2.0 * mu) * tau_m;
    b->ic_q[ic] = q;
    if (q == 0.0) return;
    const int donor = q > 0.0 ? ca : cb;
    const double rho = L->kind == 0 ? pq_density(L, b->T[donor]) : L->c_rho;
    const double dm = rho * (q > 0.0 ? q : -q) * d->dz * d->tau;
    const double phi = dm / b->m[donor];
    int bad = 0;
    if (!(phi <= d->flow_max_fraction)) { b->status[donor] = PQ_ERR_FLOW; return; }
    int ov = 0;
    for (int el = 0; el < PQ_NEL; ++el) {
        pq_i128 qn = pq_from_double(pq_to_double(b->N[PQ_NEL * donor + el], PQ_N_INV) * phi, PQ_N_POW2, &bad);
        b->ic_N[PQ_NEL * ic + el] = q > 0.0 ? qn : pq_neg(qn, &ov);
    }
    pq_i128 qe = pq_from_double(pq_to_double(b->E[donor], PQ_E_INV) * phi, PQ_E_POW2, &bad);
    b->ic_E[ic] = q > 0.0 ? qe : pq_neg(qe, &ov);
    if (bad || ov) b->status[donor] = PQ_ERR_FIXED_RANGE;
}

// 4. One explicit sub-step of one cell: conduction, prescribed-temperature exchange, Joule, surface loads,
//    deposits, releases, evaporation, radiation and melt-layer transfers.
PQ_HD void pq_step_cell(const pq_dims* d, const pq_buffers* b, int c) {
    const pq_law* L = b->law;
    for (int el = 0; el < PQ_NEL; ++el) b->N_out[PQ_NEL * c + el] = b->N[PQ_NEL * c + el];
    b->E_out[c] = b->E[c];
    if (b->flags[c] & 1) return;
    int st = b->status[c], ov = 0, bad = 0;
    pq_i128 dE = pq_zero();
    pq_i128 dN[PQ_NEL];
    for (int el = 0; el < PQ_NEL; ++el) dN[el] = pq_zero();
    const double Tc = b->T[c], kc = b->k[c];
    double atoms[PQ_NEL], x[PQ_NEL];
    for (int el = 0; el < PQ_NEL; ++el) atoms[el] = pq_to_double(b->N[PQ_NEL * c + el], PQ_N_INV);
    pq_fractions(atoms, x);
    const int s = c / (d->nx * d->ny), i = c % d->nx;
    for (int f = 0; f < 4; ++f) {
        const int nb = pq_neighbor(d, c, f);
        if (nb >= 0 && !(b->flags[nb] & 1)) {
            const double g = pq_face_g(d, b, c, f);
            if (f == 0 || f == 2) {
                const double F = pq_face_energy(b->T[nb], b->k[nb], Tc, kc, g, d->tau);
                dE = pq_add(dE, pq_from_double(F, PQ_E_POW2, &bad), &ov);
            } else {
                const double F = pq_face_energy(Tc, kc, b->T[nb], b->k[nb], g, d->tau);
                dE = pq_sub(dE, pq_from_double(F, PQ_E_POW2, &bad), &ov);
            }
            continue;
        }
        const int e = b->el_of[4 * c + f];
        if (e < 0) continue;
        pq_i128 dNf[PQ_NEL];                 // atoms gained through this face (displacement bookkeeping)
        for (int el = 0; el < PQ_NEL; ++el) dNf[el] = pq_zero();
        if (nb < 0 && d->bc_kind[f] == 1) {
            const double F = kc * pq_face_area(d, b, c, f) / pq_face_half(d, b, c, f) * (d->bc_T[f] - Tc) * d->tau;
            const pq_i128 q = pq_from_double(F, PQ_E_POW2, &bad);
            dE = pq_add(dE, q, &ov);
            b->acc_bnd_E[e] = pq_add(b->acc_bnd_E[e], q, &ov);
        }
        dE = pq_add(dE, pq_part(b->in_load_E[e], d->n_sub, d->sub, &ov), &ov);
        dE = pq_add(dE, pq_part(b->in_dep_E[e], d->n_sub, d->sub, &ov), &ov);
        for (int el = 0; el < PQ_NEL; ++el) {
            const pq_i128 p = pq_part(b->in_dep_N[PQ_NEL * e + el], d->n_sub, d->sub, &ov);
            dN[el] = pq_add(dN[el], p, &ov); dNf[el] = pq_add(dNf[el], p, &ov);
        }
        if (d->evap_mode == 2) {
            dE = pq_add(dE, pq_part(b->in_rel_E[e], d->n_sub, d->sub, &ov), &ov);
            for (int el = 0; el < PQ_NEL; ++el) {
                const pq_i128 p = pq_part(b->in_rel_N[PQ_NEL * e + el], d->n_sub, d->sub, &ov);
                dN[el] = pq_add(dN[el], p, &ov); dNf[el] = pq_add(dNf[el], p, &ov);
            }
        } else if (d->evap_mode == 1 && d->evap_face[f] && L->vapor_enabled) {
            const double A = b->el_area[e];
            const double fnet = pq_flux_fraction(L, Tc, x, b->patch_pb[b->el_patch[e]], &st);
            for (int el = 0; el < PQ_NEL; ++el) {
                const double D = pq_flux_f(L, el, Tc, x, fnet, &st) * A * d->tau;
                const pq_i128 qn = pq_from_double(D, PQ_N_POW2, &bad);
                if (pq_is_zero(qn)) continue;
                const double a = pq_to_double(qn, PQ_N_INV);
                const pq_i128 qr = pq_from_double(a * L->eps0[el], PQ_E_POW2, &bad);
                const pq_i128 qk = pq_from_double(a * (2.0 * PQ_KB * Tc), PQ_E_POW2, &bad);
                dN[el] = pq_sub(dN[el], qn, &ov); dNf[el] = pq_sub(dNf[el], qn, &ov);
                dE = pq_sub(dE, pq_add(qr, qk, &ov), &ov);
                b->acc_ev_N[PQ_NEL * e + el] = pq_add(b->acc_ev_N[PQ_NEL * e + el], qn, &ov);
                b->acc_ev_Eref[e] = pq_add(b->acc_ev_Eref[e], qr, &ov);
                b->acc_ev_Ekin[e] = pq_add(b->acc_ev_Ekin[e], qk, &ov);
                if (L->gas) {
                    const pq_i128 qi = pq_from_double(a * pq_eint(L, el, Tc), PQ_E_POW2, &bad);
                    dE = pq_sub(dE, qi, &ov);
                    b->acc_ev_Eint[e] = pq_add(b->acc_ev_Eint[e], qi, &ov);
                }
            }
        }
        if (L->radiation) {
            const double T4 = Tc * Tc * Tc * Tc, E4 = L->T_env * L->T_env * L->T_env * L->T_env;
            const pq_i128 q = pq_from_double(L->emissivity * PQ_SIGMA_SB * (T4 - E4) * b->el_area[e] * d->tau, PQ_E_POW2, &bad);
            dE = pq_sub(dE, q, &ov);
            b->acc_rad[e] = pq_add(b->acc_rad[e], q, &ov);
        }
        if (f < 2) {
            const int fc = pq_fc_index(d, s, f, i);
            if (d->flow) {
                // melt-layer transfers through this face column's two interfaces (+ toward +x)
                const int base = (s * 2 + f) * (d->nx - 1);
                if (i > 0) {
                    const int ic = base + i - 1;
                    for (int el = 0; el < PQ_NEL; ++el) { dN[el] = pq_add(dN[el], b->ic_N[PQ_NEL * ic + el], &ov); dNf[el] = pq_add(dNf[el], b->ic_N[PQ_NEL * ic + el], &ov); }
                    dE = pq_add(dE, b->ic_E[ic], &ov);
                }
                if (i + 1 < d->nx) {
                    const int ic = base + i;
                    for (int el = 0; el < PQ_NEL; ++el) { dN[el] = pq_sub(dN[el], b->ic_N[PQ_NEL * ic + el], &ov); dNf[el] = pq_sub(dNf[el], b->ic_N[PQ_NEL * ic + el], &ov); }
                    dE = pq_sub(dE, b->ic_E[ic], &ov);
                }
            }
            for (int el = 0; el < PQ_NEL; ++el) b->fcm_out[PQ_NEL * fc + el] = pq_add(b->fcm[PQ_NEL * fc + el], dNf[el], &ov);
        }
    }
    dE = pq_add(dE, pq_part(b->joule_total[c], d->n_sub, d->sub, &ov), &ov);
    const pq_i128 En = pq_add(b->E[c], dE, &ov);
    double m_new = 0.0;
    for (int el = 0; el < PQ_NEL; ++el) {
        const pq_i128 Nn = pq_add(b->N[PQ_NEL * c + el], dN[el], &ov);
        if (pq_is_neg(Nn)) pq_set_status(&st, PQ_ERR_NEGATIVE_INV);
        b->N_out[PQ_NEL * c + el] = Nn;
        m_new += pq_to_double(Nn, PQ_N_INV) * L->atom_mass[el];
    }
    if (!(m_new >= d->phi_floor * b->mref[c])) pq_set_status(&st, PQ_ERR_DEPLETED);
    if (bad || ov) pq_set_status(&st, PQ_ERR_FIXED_RANGE);
    b->E_out[c] = En;
    b->status[c] = st;
}

// 5. Per-slice exact sums of the state and the first failing cell; merge candidates (end of step).
PQ_HD void pq_slice_reduce(const pq_dims* d, const pq_buffers* b, int s) {
    int ov = 0;
    pq_i128 sE = pq_zero(), sN[PQ_NEL];
    for (int el = 0; el < PQ_NEL; ++el) sN[el] = pq_zero();
    int first = -1, code = 0, cand = 0, first_cand = -1;
    double maxT = 0.0;
    const int per = d->nx * d->ny;
    for (int c = s * per; c < (s + 1) * per; ++c) {
        if (first < 0 && b->status[c] != PQ_OK) { first = c; code = b->status[c]; }
        if (b->flags[c] & 1) continue;
        if (b->T[c] > maxT) maxT = b->T[c];
        sE = pq_add(sE, b->E[c], &ov);
        for (int el = 0; el < PQ_NEL; ++el) sN[el] = pq_add(sN[el], b->N[PQ_NEL * c + el], &ov);
        bool exposed = false;
        for (int f = 0; f < 4; ++f) if (b->el_of[4 * c + f] >= 0) exposed = true;
        if (exposed) {
            double atoms[PQ_NEL];
            const double m = pq_cell_mass(b->law, b->N, c, atoms);
            if (m < d->phi_merge * b->mref[c]) { ++cand; if (first_cand < 0) first_cand = c; }
        }
    }
    pq_i128* o = b->slice_sum + PQ_SLICE_SUM * s;
    o[0] = sE; for (int el = 0; el < PQ_NEL; ++el) o[1 + el] = sN[el];
    b->slice_status[4 * s + 0] = ov ? (first >= 0 ? first : s * per) : first;
    b->slice_status[4 * s + 1] = ov ? PQ_ERR_FIXED_RANGE : code;
    b->slice_status[4 * s + 2] = cand;
    b->slice_status[4 * s + 3] = first_cand;
    b->slice_maxT[s] = maxT;
}

// 6. Expected atoms per element over dt at the current state (quantized-release quotes).
PQ_HD void pq_quote_element(const pq_dims* d, const pq_buffers* b, int e) {
    const int c = b->el_cell[e], f = b->el_face[e];
    int st = 0;
    double atoms[PQ_NEL], x[PQ_NEL];
    for (int el = 0; el < PQ_NEL; ++el) atoms[el] = pq_to_double(b->N[PQ_NEL * c + el], PQ_N_INV);
    pq_fractions(atoms, x);
    const double fnet = pq_flux_fraction(b->law, b->T[c], x, b->patch_pb[b->el_patch[e]], &st);
    for (int el = 0; el < PQ_NEL; ++el) {
        double D = 0.0;
        if (d->evap_face[f] && b->law->vapor_enabled) D = pq_flux_f(b->law, el, b->T[c], x, fnet, &st) * b->el_area[e] * d->dt;
        b->el_quote[PQ_NEL * e + el] = D;
    }
}

// 7. Split the quantized release of one patch over its elements in proportion to the element quotes
//    (remainder to the first element with the largest quote): atoms exactly, reference energy eps0 per atom,
//    internal energy E_int(T_cell) per atom (gas = 1), actual kinetic energy split by the same weights.
PQ_HD void pq_release_split(const pq_dims* d, const pq_buffers* b, int p) {
    (void)d;
    const int e0 = b->patch_begin[p], e1 = b->patch_begin[p + 1];
    int bad = 0, ov = 0;
    for (int e = e0; e < e1; ++e) { b->in_rel_E[e] = pq_zero(); for (int el = 0; el < PQ_NEL; ++el) b->in_rel_N[PQ_NEL * e + el] = pq_zero(); }
    double eint_total = 0.0;
    for (int el = 0; el < PQ_NEL; ++el) {
        const pq_i128 R = b->rel_total_N[PQ_NEL * p + el];
        const double K = b->rel_kin_J[PQ_NEL * p + el];
        if (pq_is_zero(R) && K == 0.0) continue;
        double W = 0.0, wmax = -1.0; int emax = -1;
        for (int e = e0; e < e1; ++e) { const double w = b->el_quote[PQ_NEL * e + el]; W += w; if (w > wmax) { wmax = w; emax = e; } }
        if (!(W > 0.0)) { b->status[b->el_cell[e0]] = PQ_ERR_NEGATIVE_INV; continue; }
        const double Rd = pq_to_double(R, PQ_N_INV);
        const pq_i128 Kq = pq_from_double(K, PQ_E_POW2, &bad);
        pq_i128 usedN = pq_zero(), usedK = pq_zero();
        for (int e = e0; e <= e1; ++e) {
            if (e == emax) continue;
            const bool last = e == e1;
            const int ee = last ? emax : e;
            pq_i128 qn, qk;
            if (!last) {
                const double share = b->el_quote[PQ_NEL * e + el] / W;
                qn = pq_from_double(Rd * share, PQ_N_POW2, &bad);
                qk = pq_from_double(K * share, PQ_E_POW2, &bad);
                usedN = pq_add(usedN, qn, &ov); usedK = pq_add(usedK, qk, &ov);
            } else {
                qn = pq_sub(R, usedN, &ov); qk = pq_sub(Kq, usedK, &ov);
                if (pq_is_neg(qn) || pq_is_neg(qk)) bad = 1;
            }
            const double an = pq_to_double(qn, PQ_N_INV);
            pq_i128 qr = pq_from_double(an * b->law->eps0[el], PQ_E_POW2, &bad);
            if (b->law->gas) {
                const double ei = an * pq_eint(b->law, el, b->T[b->el_cell[ee]]);
                qr = pq_add(qr, pq_from_double(ei, PQ_E_POW2, &bad), &ov);
                eint_total += ei;
            }
            b->in_rel_N[PQ_NEL * ee + el] = pq_neg(qn, &ov);
            b->in_rel_E[ee] = pq_sub(b->in_rel_E[ee], pq_add(qr, qk, &ov), &ov);
        }
    }
    b->patch_out[PQ_PATCH_OUT * p + PO_EINT] = eint_total;
    if (bad || ov) b->status[b->el_cell[e0]] = PQ_ERR_FIXED_RANGE;
}

// 8. Patch report at the current state and exact sums of the step accumulators of its elements.
PQ_HD void pq_patch_report(const pq_dims* d, const pq_buffers* b, int p) {
    const int e0 = b->patch_begin[p], e1 = b->patch_begin[p + 1];
    double A = 0, sT = 0, mT = 0, sf = 0, sdepth = 0, srho = 0, seta = 0, meta = 0, srec = 0, sspeed = 0, mre = 0;
    double sD[PQ_NEL], sj[PQ_NEL], sx[PQ_NEL];
    for (int el = 0; el < PQ_NEL; ++el) { sD[el] = 0; sj[el] = 0; sx[el] = 0; }
    int ov = 0, st = 0, nspeed = 0;
    pq_i128 acc[PQ_PATCH_SUM];
    for (int q = 0; q < PQ_PATCH_SUM; ++q) acc[q] = pq_zero();
    const pq_law* L = b->law;
    for (int e = e0; e < e1; ++e) {
        const int c = b->el_cell[e], f = b->el_face[e];
        const double a = b->el_area[e], T = b->T[c];
        A += a; sT += a * T; if (T > mT) mT = T;
        sf += a * pq_fl(L, T);
        if (L->kind == 0) srho += a * pq_density(L, T);
        double depth = 0.0;
        int cc = c;
        const int inward = f == 0 ? 1 : f == 1 ? 0 : f == 2 ? 3 : 2;
        while (cc >= 0 && !(b->flags[cc] & 1)) {
            const double fl = pq_fl(L, b->T[cc]);
            if (fl <= 0.0) break;
            depth += fl * (f < 2 ? b->row_dy[(cc / d->nx) % d->ny] : d->dx);
            cc = pq_neighbor(d, cc, inward);
        }
        sdepth += a * depth;
        double atoms[PQ_NEL], x[PQ_NEL];
        for (int el = 0; el < PQ_NEL; ++el) atoms[el] = pq_to_double(b->N[PQ_NEL * c + el], PQ_N_INV);
        pq_fractions(atoms, x);
        const double pb = b->patch_pb[p];
        const double fnet = pq_flux_fraction(L, T, x, pb, &st);
        for (int el = 0; el < PQ_NEL; ++el) {
            const double j = (d->evap_face[f] && L->vapor_enabled) ? pq_flux_f(L, el, T, x, fnet, &st) : 0.0;
            sj[el] += a * j;
            sD[el] += j * a * d->dt;
            sx[el] += a * x[el];
        }
        if (L->vapor_enabled && d->evap_face[f]) srec += a * pq_recoil(L, T, x, pb, &st);
        if (f < 2) {
            const int s = c / (d->nx * d->ny), i = c % d->nx;
            const int fc = pq_fc_index(d, s, f, i);
            const double eta = b->fc_eta[fc];
            seta += a * eta; if (pq_fabs(eta) > meta) meta = pq_fabs(eta);
            if (d->flow && i + 1 < d->nx) {
                const int ic = (s * 2 + f) * (d->nx - 1) + i;
                const double dl = 0.5 * (b->fc_delta[fc] + b->fc_delta[fc + 1]);
                if (dl > 0.0) {
                    sspeed += pq_fabs(b->ic_q[ic]) / dl; ++nspeed;
                    double sg, dsg, mu; pq_surface_props(L, b->fc_T[fc], &sg, &dsg, &mu);
                    const double rho = L->kind == 0 ? pq_density(L, b->fc_T[fc]) : L->c_rho;
                    const double re = rho * pq_fabs(b->ic_q[ic]) / mu;
                    if (re > mre) mre = re;
                }
            }
        }
        for (int el = 0; el < PQ_NEL; ++el) acc[PS_EVN + el] = pq_add(acc[PS_EVN + el], b->acc_ev_N[PQ_NEL * e + el], &ov);
        acc[PS_EREF] = pq_add(acc[PS_EREF], b->acc_ev_Eref[e], &ov);
        acc[PS_EKIN] = pq_add(acc[PS_EKIN], b->acc_ev_Ekin[e], &ov);
        acc[PS_EINT] = pq_add(acc[PS_EINT], b->acc_ev_Eint[e], &ov);
        acc[PS_BND] = pq_add(acc[PS_BND], b->acc_bnd_E[e], &ov);
        acc[PS_LOAD] = pq_add(acc[PS_LOAD], b->in_load_E[e], &ov);
        acc[PS_DEP] = pq_add(acc[PS_DEP], b->in_dep_E[e], &ov);
        acc[PS_REL] = pq_add(acc[PS_REL], b->in_rel_E[e], &ov);
        acc[PS_RAD] = pq_add(acc[PS_RAD], b->acc_rad[e], &ov);
    }
    double* o = b->patch_out + PQ_PATCH_OUT * p;
    o[PO_AREA] = A; o[PO_MEAN_T] = A > 0 ? sT / A : 0; o[PO_MAX_T] = mT; o[PO_FL] = A > 0 ? sf / A : 0; o[PO_MELT] = A > 0 ? sdepth / A : 0;
    for (int el = 0; el < PQ_NEL; ++el) { o[PO_D + el] = sD[el]; o[PO_J + el] = A > 0 ? sj[el] / A : 0; o[PO_X + el] = A > 0 ? sx[el] / A : 0; }
    o[PO_RHO] = A > 0 ? srho / A : 0;
    o[PO_ETA_MEAN] = A > 0 ? seta / A : 0; o[PO_ETA_MAX] = meta; o[PO_RECOIL] = A > 0 ? srec / A : 0;
    o[PO_SPEED] = nspeed > 0 ? sspeed / nspeed : 0;
    o[PO_RE] = mre;
    for (int q = 0; q < PQ_PATCH_SUM; ++q) b->patch_sum[PQ_PATCH_SUM * p + q] = acc[q];
    (void)ov;
}

#if defined(__CUDACC_RTC__)
#define PQ_KERNEL(name, count, call) extern "C" __global__ void name(pq_dims d, pq_buffers b) { const int t = blockIdx.x * blockDim.x + threadIdx.x; if (t < (count)) call(&d, &b, t); }
PQ_KERNEL(k_thermo, d.n_cells, pq_thermo)
PQ_KERNEL(k_scatter, d.n_scatter, pq_scatter)
PQ_KERNEL(k_sigma, d.n_em, pq_sigma_em)
PQ_KERNEL(k_joule, d.n_em, pq_joule_split)
PQ_KERNEL(k_fc_prepare, d.n_fc, pq_fc_prepare)
PQ_KERNEL(k_fc_pressure, d.n_fc, pq_fc_pressure)
PQ_KERNEL(k_ic_flux, d.n_ic, pq_ic_flux)
PQ_KERNEL(k_step, d.n_cells, pq_step_cell)
PQ_KERNEL(k_slice, d.S, pq_slice_reduce)
PQ_KERNEL(k_quote, d.n_el, pq_quote_element)
PQ_KERNEL(k_release, d.n_patches, pq_release_split)
PQ_KERNEL(k_patch, d.n_patches, pq_patch_report)
#endif
