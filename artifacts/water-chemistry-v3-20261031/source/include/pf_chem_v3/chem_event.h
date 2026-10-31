// chem_event.h -- chemistry v3 per-carrier collision step (null-collision
// Monte Carlo), shared by the CPU reference and the NVRTC kernel.
//
// Model of one step (operator split; the carrier sits at its post-push
// position in one cell during the collision phase):
//  * The carrier's collision process is a Poisson process whose rate depends
//    only on its current state. Candidate times are drawn with the majorant
//    rate nu_max(state) = sum_k n_k M_k(state), M_k an upper bound of
//    sigma_tot(E') v' over every target-frame energy E' reachable with a target
//    velocity inside the truncated thermal sphere |u_t - drift| <= trunc*s_k.
//    Waiting time tau = -ln(U0)/nu_max (fresh exponential after every
//    candidate; exact by memorylessness, and the majorant is re-evaluated
//    whenever a real event changes the state).
//  * At a candidate: target slot k with probability n_k M_k / nu_max (U1); a
//    target velocity from the truncated Maxwellian (Box-Muller, U4..U19);
//    real/null and channel with one uniform U2: channel j is chosen iff
//    U2 M_k falls in [sum_{i<j} sigma_i v', sum_{i<=j} sigma_i v'); null if
//    U2 M_k >= sigma_tot v'.
//  * Every number is drawn from the chemistry stream (seed, step, carrier id,
//    candidate index, word), so results are independent of the order and
//    partition in which carriers are processed.
// The kernel writes one CarrierOut and up to max_events EventRec per carrier.
// Ledgers, target-inventory debits and neutral releases are formed on the
// host from these records (chem_host.hpp).
#pragma once
#include "chem_kinematics.h"
#include "chem_rng.h"
#include "chem_tables.h"

namespace pfc {

// Carrier status codes.
#define PFC_ST_OK 0u
#define PFC_ST_CONSUMED 1u          // projectile removed by a real event (still OK)
#define PFC_ST_CAND_CAP 2u          // candidate cap reached before dt was exhausted (refuse step)
#define PFC_ST_EVENT_CAP 3u         // real-event record cap reached (refuse step)
#define PFC_ST_KEY_DOMAIN 4u        // id/step/event outside the stream domain (refuse step)
#define PFC_ST_MAJORANT 5u          // sigma_tot v' exceeded the majorant (refuse step; never expected)
#define PFC_ST_COVER_CAP 6u         // more distinct unresolved regions than slots (refuse step)
#define PFC_ST_BAD_INPUT 7u         // nonfinite input, unknown species or cell (refuse step)

#define PFC_MAX_COVER_SLOTS 16

struct StepConfig {
    pfc_u64 seed, step;
    double dt_s;
    double trunc_sigmas;      // truncation radius of the target Maxwellian in standard deviations
    pfc_u32 purpose;
    pfc_u32 max_candidates;   // per carrier per step
    pfc_u32 max_events;       // real events per carrier per step (record slots)
    pfc_u32 pad;
};
struct CarrierIn {           // 48 bytes
    pfc_u64 id;
    double u[3];             // proper velocity, m/s
    double weight;           // physical particles represented
    pfc_u32 species, cell;
};
struct CellTargets {         // targets present in one cell
    double n[PFC_MAX_TARGETS];           // m^-3
    double drift[PFC_MAX_TARGETS][3];    // mean proper velocity, m/s
    double s[PFC_MAX_TARGETS];           // 1D thermal standard deviation sqrt(kT/m), m/s
    pfc_u32 species[PFC_MAX_TARGETS];
    pfc_u32 n_targets, pad;
};
struct EventRec {            // one real event (sizes asserted on the host)
    pfc_u64 carrier;
    pfc_u32 candidate, process, cell, target_species;
    pfc_u32 n_out, cover_class, internal_sink, flags;
    double weight;
    double time_s;           // event time within the step
    double u_in[3];          // projectile proper velocity before (m/s)
    double ut[3];            // sampled target proper velocity (m/s)
    double E_tf;             // projectile kinetic energy in the target frame (eV)
    double ke_in_p, ke_in_t; // canonical kinetic energies (eV, per physical particle)
    double q_eV, e_internal; // reference-energy change and internal sink (eV)
    pfc_u32 out_species[PFC_MAX_OUT];
    pfc_u32 out_dest[PFC_MAX_OUT];
    double out_u[PFC_MAX_OUT][3];
    double ke_out[PFC_MAX_OUT];
};
struct CarrierOut {
    pfc_u64 id;
    pfc_u32 status, n_events, n_candidates, species_final;
    double u_final[3];
    double nu_dt_max;            // max over the step of nu_max * dt (time-step guard quantity)
    pfc_u32 n_cover, pad;
    pfc_u32 cover_region[PFC_MAX_COVER_SLOTS];
    double cover_exposure[PFC_MAX_COVER_SLOTS];   // sum W n v dt (sigma-free), m^-2
    double cover_bound[PFC_MAX_COVER_SLOTS];      // sum W n v sigma_bound dt (expected events bound)
};

// ---------------------------------------------------------------- BEB / RBEB
// One orbital, Q = 1 (Kim & Rudd 1994; Kim, Santos & Parente 2000).
struct SdcsShape {
    double t, phi, cb, L, P;    // P in m^2
    double B;
};
PFC_FN double beta2_of(double x) { return x * (x + 2.0) / ((1.0 + x) * (1.0 + x)); }
PFC_FN SdcsShape sdcs_shape(double B, double U, double N, double T, double a0, double R, double mc2, int rel) {
    SdcsShape s;
    s.B = B;
    s.t = T / B;
    const double pi = 3.141592653589793;
    if (!rel) {
        s.phi = 1.0;
        s.cb = 0.0;
        s.L = pfc::log(s.t);
        const double u = U / B;
        const double S = 4.0 * pi * a0 * a0 * N * (R / B) * (R / B);
        s.P = S / (s.t + u + 1.0);
    } else {
        const double tp = T / mc2, bp = B / mc2, up = U / mc2;
        const double bt2 = beta2_of(tp), bb2 = beta2_of(bp), bu2 = beta2_of(up);
        const double alpha2 = 2.0 * R / mc2;
        s.P = 4.0 * pi * a0 * a0 * alpha2 * alpha2 * N / ((bt2 + bu2 + bb2) * 2.0 * bp);
        s.L = pfc::log(tp * (tp + 2.0)) - bt2 - pfc::log(2.0 * bp);
        const double h = 1.0 + tp / 2.0;
        s.phi = (1.0 + 2.0 * tp) / (h * h);
        s.cb = bp * bp / (h * h);
    }
    return s;
}
PFC_FN double sdcs_G(const SdcsShape& s, double w) {
    const double t = s.t, tw = t - w;
    return -s.phi / (t + 1.0) * (pfc::log1p(w) - pfc::log1p(-w / t)) + w / (w + 1.0) + w / (t * tw) + s.cb * w +
           s.L / 2.0 * (w * (w + 2.0) / ((w + 1.0) * (w + 1.0)) + w * (2.0 * t - w) / (t * t * tw * tw));
}
PFC_FN double sdcs_total(const SdcsShape& s) { return s.P * sdcs_G(s, (s.t - 1.0) / 2.0); }
// Reduced ejected energy w in [0, (t-1)/2] with G(w) = u G(wmax); 64 bisections.
PFC_FN double sdcs_sample_w(const SdcsShape& s, double u) {
    const double wm = (s.t - 1.0) / 2.0;
    const double target = u * sdcs_G(s, wm);
    double lo = 0.0, hi = wm;
    for (int i = 0; i < 64; ++i) {
        const double mid = lo + (hi - lo) / 2.0;
        if (!(mid > lo && mid < hi)) break;
        if (sdcs_G(s, mid) < target) lo = mid;
        else hi = mid;
    }
    return hi;
}

// ------------------------------------------------------------- chain decay
// Decay a virtual particle of rest mass Mv (eV) with momentum Pv in frame X
// into the final species list f[0..n-1] (masses m[]), with total kinetic
// energy Kkin in its rest frame (flat sequential partition for n > 2,
// isotropic directions). Momenta in X are written to pout[]. Words are taken
// from word index *wp onward (3 per split). Returns 0 on kinematic failure.
PFC_FN int chain_decay(const ChemKey& key, pfc_u32 ev, pfc_u32* wp, double Mv, V3 Pv, const double* m, pfc_u32 n,
                       double Kkin, V3* pout) {
    if (n == 1) {
        pout[0] = Pv;
        return 1;
    }
    double rest = 0.0;
    for (pfc_u32 i = 0; i < n; ++i) rest += m[i];
    double Krem = Kkin;
    V3 P = Pv;
    double M = Mv;
    for (pfc_u32 i = 0; i + 1 < n; ++i) {
        rest -= m[i];   // rest mass of the remainder after removing fragment i
        double Ks;
        if (i + 2 == n) Ks = Krem;
        else Ks = Krem * chem_u01(key, ev, (*wp)++);
        if (!(Ks >= 0.0)) return 0;
        const double Mr = rest + (Krem - Ks);   // virtual mass of the remainder
        const double ps = pstar(Ks, m[i], Mr);
        const double ua = chem_u01(key, ev, (*wp)++), ub = chem_u01(key, ev, (*wp)++);
        const V3 d = iso_dir(ua, ub);
        const V3 w = scale(P, 1.0 / M);
        const double gw = pfc_sqrt(1.0 + norm2(w));
        const TwoBody tb = two_body_tf(w, gw, m[i], Mr, ps, d);
        pout[i] = tb.p3;
        P = tb.p4;
        M = Mr;
        Krem = Krem - Ks;
    }
    pout[n - 1] = P;
    return 1;
}

// ---------------------------------------------------------------- one event
// Inputs in the target frame (TF): projectile momentum pp (eV), its kinetic
// energy T, projectile mass mp, target mass mt, target proper velocity wt
// (dimensionless) with gamma gt (for the final boost to the lab). On success
// writes lab momenta plab[] and kinematic masses mk[] (rest energy plus the
// internal energy I_i the product carries) of the n_out final particles in
// record order, and the lab-frame internal-sink energy sum_i gamma_i I_i
// (the internal energy is released in each product's rest frame, so its lab
// energy is gamma_i I_i and the product keeps its velocity); returns 0 when
// the event is kinematically closed.
// Internal-energy carriers: EXCITE -> target (I = loss); DISSOC -> last
// fragment (I = (1 - f) (loss - Q_chain)); IONIZE -> ion, first chain product
// (I = (1 - f)(B_eff - q)); CAPTURE -> first product (I = (1 - f)(K_cm - q)).
PFC_FN int event_kinematics(const TableView& t, const ProcessRec& pr, const GroupRec& gr, int gi, const ChemKey& key,
                            pfc_u32 ev, double mp, double mt, V3 pp, double T, V3 wt, double gt, V3* plab, double* mk,
                            double* e_int) {
    pfc_u32 wp = PFC_W_KIN;
    V3 ptf[PFC_MAX_OUT];
    double I[PFC_MAX_OUT];
    for (pfc_u32 i = 0; i < pr.n_out; ++i) {
        mk[i] = t.species[pr.out_species[i]].mass_eV;
        I[i] = 0.0;
    }
    *e_int = 0.0;
    const double Kcm = kcm_from_T(mp, mt, T);
    const double sqrt_s = mp + mt + Kcm;
    const V3 wcm = scale(pp, 1.0 / sqrt_s);
    const double gcm = pfc_sqrt(1.0 + norm2(wcm));
    const V3 din = unit_or_z(pp);
    if (pr.kind == PFC_P_ELASTIC || pr.kind == PFC_P_EXCITE || pr.kind == PFC_P_DISSOC || pr.kind == PFC_P_CHARGE_X) {
        // Scatter family: (out0, out1-or-virtual) in CM.
        double m3, m4, Kout;
        if (pr.kind == PFC_P_CHARGE_X) {
            m3 = mk[0];
            m4 = mk[1];
            Kout = Kcm - pr.q_eV;
        } else {
            m3 = mp;
            m4 = mt + pr.loss_eV;
            Kout = Kcm - pr.loss_eV;
        }
        if (!(Kout >= 0.0)) return 0;
        const double ua = chem_u01(key, ev, wp++), ub = chem_u01(key, ev, wp++);
        V3 d;
        if (pr.angular == PFC_ANG_ISO) d = iso_dir(ua, ub);
        else if (pr.angular == PFC_ANG_FORWARD) d = din;
        else {
            const double xi = xi_eval(t, gr, gi, T);
            d = deflect(din, screened_mu(xi, ua), ub);
        }
        const double ps = pstar(Kout, m3, m4);
        const TwoBody tb = two_body_tf(wcm, gcm, m3, m4, ps, d);
        ptf[0] = tb.p3;
        if (pr.kind == PFC_P_DISSOC) {
            // virtual excited target -> fragments out[1..]; internal on the last fragment
            const double avail = pr.loss_eV - pr.chain_q_eV;
            if (!(avail >= 0.0)) return 0;
            const double kin = avail * pr.kinetic_fraction;
            const pfc_u32 nf = pr.n_out - 1;
            I[pr.n_out - 1] = avail - kin;
            mk[pr.n_out - 1] += I[pr.n_out - 1];
            double Mv = kin;
            for (pfc_u32 i = 1; i < pr.n_out; ++i) Mv += mk[i];
            if (!chain_decay(key, ev, &wp, Mv, tb.p4, mk + 1, nf, kin, ptf + 1)) return 0;
        } else {
            ptf[1] = tb.p4;
            if (pr.kind == PFC_P_EXCITE) {
                I[1] = pr.loss_eV;
                mk[1] = m4;
            }
        }
    } else if (pr.kind == PFC_P_CAPTURE) {
        const double avail = Kcm - pr.q_eV;
        if (!(avail >= 0.0)) return 0;
        const double kin = avail * pr.kinetic_fraction;
        I[0] = avail - kin;
        mk[0] += I[0];
        double Mv = kin;
        for (pfc_u32 i = 0; i < pr.n_out; ++i) Mv += mk[i];
        // compound at rest in CM; decay there, then boost CM -> TF
        V3 pc[PFC_MAX_OUT];
        if (!chain_decay(key, ev, &wp, Mv, v3(0.0, 0.0, 0.0), mk, pr.n_out, kin, pc)) return 0;
        for (pfc_u32 i = 0; i < pr.n_out; ++i) ptf[i] = boost_out(wcm, gcm, pc[i], energy_of(mk[i], pc[i]));
    } else if (pr.kind == PFC_P_IONIZE) {
        // out[0] scattered electron (survivor), out[1] ejected electron (birth),
        // out[2..] ion-system products (chain; the ion is out[2]).
        const double me = mp;
        // orbital selection by w_i sigma_i(T) with B_eff = max(B_i, q)
        double terms[PFC_MAX_ORB];
        double total = 0.0;
        const int rel = T > pr.relativistic_from_eV ? 1 : 0;
        for (pfc_u32 i = 0; i < pr.n_orb; ++i) {
            terms[i] = 0.0;
            const double Be = pr.orb_B[i] > pr.q_eV ? pr.orb_B[i] : pr.q_eV;
            if (pr.orb_w[i] > 0.0 && T > Be) {
                const SdcsShape s = sdcs_shape(Be, pr.orb_U[i], pr.orb_N[i], T, pr.bohr_m, pr.rydberg_eV, me, rel);
                const double v = pr.orb_w[i] * sdcs_total(s);
                terms[i] = v > 0.0 ? v : 0.0;
            }
            total += terms[i];
        }
        if (!(total > 0.0)) return 0;
        const double uo = chem_u01(key, ev, wp++) * total;
        pfc_u32 io = 0;
        double cum = 0.0;
        for (pfc_u32 i = 0; i < pr.n_orb; ++i) {
            if (!(terms[i] > 0.0)) continue;
            io = i;
            cum += terms[i];
            if (uo < cum) break;
        }
        const double Be = pr.orb_B[io] > pr.q_eV ? pr.orb_B[io] : pr.q_eV;
        const SdcsShape s = sdcs_shape(Be, pr.orb_U[io], pr.orb_N[io], T, pr.bohr_m, pr.rydberg_eV, me, rel);
        const double W = Be * sdcs_sample_w(s, chem_u01(key, ev, wp++));
        // binary-encounter ejection angle about the incident direction
        double c2 = W * (T + 2.0 * me) / (T * (W + 2.0 * me));
        if (c2 > 1.0) c2 = 1.0;
        const double mu = pfc_sqrt(c2);
        const V3 dw = deflect(din, mu, chem_u01(key, ev, wp++));
        const V3 pw = scale(dw, pfc_sqrt(W * (W + 2.0 * me)));
        const double Ex = Be - pr.q_eV;               // ion-system excitation above ground products
        const double Sigma = mt + Be;                // m_scattered + M_x (M_x = mt - me + q + Ex)
        const double D = T - W - Be;
        const V3 pr_ = sub(pp, pw);
        const double x = D * (D + 2.0 * Sigma) - norm2(pr_);
        if (!(x >= 0.0)) return 0;
        const double Kr = x / (pfc_sqrt(Sigma * Sigma + x) + Sigma);
        const double sr = Sigma + Kr;
        const V3 wr = scale(pr_, 1.0 / sr);
        const double gr_ = pfc_sqrt(1.0 + norm2(wr));
        const double Mx = Sigma - me;                // ion-system mass (with Ex)
        const double ps = pstar(Kr, me, Mx);
        const TwoBody tb = two_body_tf(wr, gr_, me, Mx, ps, unit_or_z(pr_));
        ptf[0] = tb.p3;
        ptf[1] = pw;
        // ion system: products out[2..]; internal energy on the ion
        const double kin = Ex * pr.kinetic_fraction;
        I[2] = Ex - kin;
        mk[2] += I[2];
        double Mv = kin;
        for (pfc_u32 i = 2; i < pr.n_out; ++i) Mv += mk[i];
        if (!chain_decay(key, ev, &wp, Mv, tb.p4, mk + 2, pr.n_out - 2, kin, ptf + 2)) return 0;
    } else {
        return 0;
    }
    double eint = 0.0;
    for (pfc_u32 i = 0; i < pr.n_out; ++i) {
        plab[i] = boost_out(wt, gt, ptf[i], energy_of(mk[i], ptf[i]));
        if (I[i] > 0.0) {
            const double g = pfc_sqrt(1.0 + norm2(plab[i]) / (mk[i] * mk[i]));
            eint += g * I[i];
        }
    }
    *e_int = eint;
    return 1;
}

// ------------------------------------------------------------- coverage add
PFC_FN int cover_add(CarrierOut* o, pfc_u32 region, double x, double b) {
    for (pfc_u32 i = 0; i < o->n_cover; ++i)
        if (o->cover_region[i] == region) {
            o->cover_exposure[i] += x;
            o->cover_bound[i] += b;
            return 1;
        }
    if (o->n_cover >= PFC_MAX_COVER_SLOTS) return 0;
    o->cover_region[o->n_cover] = region;
    o->cover_exposure[o->n_cover] = x;
    o->cover_bound[o->n_cover] = b;
    o->n_cover++;
    return 1;
}

// Coverage bins of a group are stored as rows of the bound table: the group's
// bin edges are its grid points; bin i (grid interval i) has region id
// bin_region[group.bin_off + i] (0xFFFFFFFF when every process is resolved
// there) and bound row bin_row[...] (0xFFFFFFFF when no bound exists). The
// two regions "below the grid" and "above the grid" have ids bin_low/bin_high.
struct CoverView {
    const pfc_u32* bin_region;   // per group: grid_n + 1 entries (below, intervals 0..n-2, above)
    const pfc_u32* bin_row;
    const pfc_u32* bin_off;      // per group offset
};

PFC_FN pfc_u32 bin_index(int i, pfc_u32 n) {
    if (i < 0) return 0u;
    if (i >= (int)n - 1) return n;   // at/above the last point
    return (pfc_u32)i + 1u;
}

// --------------------------------------------------------------- the step
PFC_FN void mcc_carrier(const TableView& t, const CoverView& cv, const StepConfig& cfg, const CarrierIn& cin,
                        const CellTargets& cell, EventRec* evs, CarrierOut* out) {
    const double c = PFC_C_M_S;
    out->id = cin.id;
    out->status = PFC_ST_OK;
    out->n_events = 0;
    out->n_candidates = 0;
    out->species_final = cin.species;
    out->u_final[0] = cin.u[0];
    out->u_final[1] = cin.u[1];
    out->u_final[2] = cin.u[2];
    out->nu_dt_max = 0.0;
    out->n_cover = 0;
    out->pad = 0;
    for (int i = 0; i < PFC_MAX_COVER_SLOTS; ++i) {
        out->cover_region[i] = 0xFFFFFFFFu;
        out->cover_exposure[i] = 0.0;
        out->cover_bound[i] = 0.0;
    }
    ChemKey key;
    key.seed = cfg.seed;
    key.step = cfg.step;
    key.carrier = cin.id;
    key.purpose = cfg.purpose;
    if (!chem_key_valid(key, 0u) || cin.species >= t.n_species) {
        out->status = cin.species >= t.n_species ? PFC_ST_BAD_INPUT : PFC_ST_KEY_DOMAIN;
        return;
    }
    const double m = t.species[cin.species].mass_eV;
    V3 u = v3(cin.u[0], cin.u[1], cin.u[2]);
    if (!(norm2(u) < 1e300) || !(cin.weight > 0.0) || !(cin.weight < 1e300)) {
        out->status = PFC_ST_BAD_INPUT;
        return;
    }
    const double W = cin.weight;
    double tacc = 0.0;
    pfc_u32 e = 0;
    for (;;) {
        const V3 p = p_from_u(m, u);
        const double pn = pfc_sqrt(norm2(p));
        const double K = ke_from_u(m, u);
        // majorant per target slot
        double nu[PFC_MAX_TARGETS], M[PFC_MAX_TARGETS];
        pfc_u32 gidx[PFC_MAX_TARGETS];
        double numax = 0.0;
        for (pfc_u32 k = 0; k < PFC_MAX_TARGETS; ++k) {
            nu[k] = 0.0;
            M[k] = 0.0;
            gidx[k] = 0xFFFFFFFFu;
        }
        for (pfc_u32 k = 0; k < cell.n_targets && k < PFC_MAX_TARGETS; ++k) {
            const pfc_u32 ts = cell.species[k];
            const pfc_u32 g = group_index(t, cin.species, ts);
            if (g == 0xFFFFFFFFu || !(cell.n[k] > 0.0)) continue;
            gidx[k] = g;
            const GroupRec& gr = t.groups[g];
            const double vcap = pfc_sqrt(norm2(v3(cell.drift[k][0], cell.drift[k][1], cell.drift[k][2]))) + cfg.trunc_sigmas * cell.s[k];
            const double wc = vcap / c;
            const double gc = pfc_sqrt(1.0 + wc * wc);
            const double base = gc * K + gamma_minus_one(wc * wc) * m;
            const double Ehi = (base + wc * pn) * (1.0 + 1e-12);
            // projectile at rest in some admissible target frame?  (beta_cap >= beta_p)
            const double wp = pn / m;   // projectile proper speed (dimensionless)
            double Elo = (wc >= wp) ? 0.0 : (base - wc * pn) * (1.0 - 1e-12);
            if (Elo < 0.0) Elo = 0.0;
            M[k] = majorant(t, gr, Elo, Ehi);
            nu[k] = cell.n[k] * M[k];
            numax += nu[k];
        }
        const double remaining = cfg.dt_s - tacc;
        // waiting time to the next candidate
        double tau;
        if (numax > 0.0) {
            const double g = numax * cfg.dt_s;
            if (g > out->nu_dt_max) out->nu_dt_max = g;
            if (e >= cfg.max_candidates) {
                out->status = PFC_ST_CAND_CAP;
                break;
            }
            if (!chem_key_valid(key, e)) {
                out->status = PFC_ST_KEY_DOMAIN;
                break;
            }
            tau = -pfc::log(chem_u01(key, e, PFC_W_TIME)) / numax;
        } else {
            tau = remaining;
        }
        const double seg = tau < remaining ? tau : remaining;
        // coverage exposure of the current state over seg (energy relative to the target drift)
        for (pfc_u32 k = 0; k < cell.n_targets && k < PFC_MAX_TARGETS; ++k) {
            const pfc_u32 g = group_index(t, cin.species, cell.species[k]);
            if (g == 0xFFFFFFFFu || !(cell.n[k] > 0.0) || !(seg > 0.0)) continue;
            const GroupRec& gr = t.groups[g];
            const V3 wd = w_from_u(v3(cell.drift[k][0], cell.drift[k][1], cell.drift[k][2]));
            const double gd = pfc_sqrt(1.0 + norm2(wd));
            const V3 pd = boost_in(wd, gd, p, m + K);
            const double p2d = norm2(pd);
            const double Ed = kinetic_from_p2(m, p2d);
            const double vd = pfc_sqrt(p2d) / (m + Ed) * c;
            const int gi = grid_locate(t.grid + gr.grid_off, gr.grid_n, Ed);
            const pfc_u32 b = bin_index(gi, gr.grid_n);
            const pfc_u32 region = cv.bin_region[cv.bin_off[g] + b];
            if (region == 0xFFFFFFFFu) continue;
            const pfc_u32 row = cv.bin_row[cv.bin_off[g] + b];
            const double x = W * cell.n[k] * vd * seg;
            double bnd = 0.0;
            if (row != 0xFFFFFFFFu) bnd = x * bound_eval(t, gr, row, gi, Ed);
            if (!cover_add(out, region, x, bnd)) {
                out->status = PFC_ST_COVER_CAP;
            }
        }
        if (out->status != PFC_ST_OK) break;
        if (!(numax > 0.0) || !(tau < remaining)) break;   // no further candidate in this step
        tacc += tau;
        // ---- candidate e
        double xs = chem_u01(key, e, PFC_W_TARGET) * numax;
        pfc_u32 k = 0xFFFFFFFFu, last = 0xFFFFFFFFu;
        double cum = 0.0;
        for (pfc_u32 j = 0; j < cell.n_targets && j < PFC_MAX_TARGETS; ++j) {
            if (!(nu[j] > 0.0)) continue;
            last = j;
            cum += nu[j];
            if (xs < cum) {
                k = j;
                break;
            }
        }
        if (k == 0xFFFFFFFFu) k = last;   // rounding at the top of the cumulative sum
        const GroupRec& gr = t.groups[gidx[k]];
        const pfc_u32 ts = cell.species[k];
        const double mt = t.species[ts].mass_eV;
        // truncated Maxwellian target velocity (up to 4 attempts, then the drift itself)
        V3 ut = v3(cell.drift[k][0], cell.drift[k][1], cell.drift[k][2]);
        {
            const double r2max = cfg.trunc_sigmas * cfg.trunc_sigmas;
            for (pfc_u32 a = 0; a < 4u; ++a) {
                const pfc_u32 w0 = PFC_W_TVEL + 4u * a;
                const double r1 = pfc_sqrt(-2.0 * pfc::log(chem_u01(key, e, w0)));
                const double r2 = pfc_sqrt(-2.0 * pfc::log(chem_u01(key, e, w0 + 2u)));
                double c1, s1, c2, s2;
                sincos_2pi(chem_u01(key, e, w0 + 1u), &c1, &s1);
                sincos_2pi(chem_u01(key, e, w0 + 3u), &c2, &s2);
                const double z0 = r1 * c1, z1 = r1 * s1, z2 = r2 * c2;
                if (z0 * z0 + z1 * z1 + z2 * z2 <= r2max) {
                    ut = add(ut, scale(v3(z0, z1, z2), cell.s[k]));
                    break;
                }
            }
        }
        const V3 wt = w_from_u(ut);
        const double gt = pfc_sqrt(1.0 + norm2(wt));
        const V3 pt = boost_in(wt, gt, p, m + K);
        const double T = kinetic_from_p2(m, norm2(pt));
        const double vt = pfc_sqrt(norm2(pt)) / (m + T) * c;
        const int gi = grid_locate(t.grid + gr.grid_off, gr.grid_n, T);
        const double y = chem_u01(key, e, PFC_W_CHANNEL) * M[k];
        pfc_u32 sel = 0xFFFFFFFFu;
        double F = 0.0;
        for (pfc_u32 j = 0; j < gr.proc_n; ++j) {
            const double sv = sigma_eval(t, gr, j, gi, T) * vt;
            F += sv;
            if (sel == 0xFFFFFFFFu && y < F) sel = j;
        }
        if (F > M[k]) {
            out->status = PFC_ST_MAJORANT;
            break;
        }
        if (sel == 0xFFFFFFFFu) {   // null candidate
            ++e;
            continue;
        }
        if (out->n_events >= cfg.max_events) {
            out->status = PFC_ST_EVENT_CAP;
            break;
        }
        const pfc_u32 pidx = gr.proc_off + sel;
        const ProcessRec& pr = t.procs[pidx];
        V3 plab[PFC_MAX_OUT];
        double mk[PFC_MAX_OUT];
        double eint = 0.0;
        const int ok = event_kinematics(t, pr, gr, gi, key, e, m, mt, pt, T, wt, gt, plab, mk, &eint);
        if (!ok) {   // kinematically closed at this energy: counted as a null candidate
            ++e;
            continue;
        }
        EventRec& r = evs[out->n_events];
        r.carrier = cin.id;
        r.candidate = e;
        r.process = pidx;
        r.cell = cin.cell;
        r.target_species = ts;
        r.n_out = pr.n_out;
        {
            pfc_u32 reg;
            r.cover_class = coverage_class(pr, T, &reg);
        }
        r.internal_sink = pr.internal_sink;
        r.flags = 0;
        r.weight = W;
        r.time_s = tacc;
        r.u_in[0] = u.x;
        r.u_in[1] = u.y;
        r.u_in[2] = u.z;
        r.ut[0] = ut.x;
        r.ut[1] = ut.y;
        r.ut[2] = ut.z;
        r.E_tf = T;
        r.ke_in_p = K;
        r.ke_in_t = ke_from_u(mt, ut);
        r.q_eV = pr.q_eV;
        r.e_internal = eint;
        int survivor = -1;
        for (pfc_u32 i = 0; i < PFC_MAX_OUT; ++i) {
            if (i < pr.n_out) {
                const double mo = t.species[pr.out_species[i]].mass_eV;
                const V3 uo = u_from_p(mk[i], plab[i]);   // velocity of the (possibly excited) product
                r.out_species[i] = pr.out_species[i];
                r.out_dest[i] = pr.out_dest[i];
                r.out_u[i][0] = uo.x;
                r.out_u[i][1] = uo.y;
                r.out_u[i][2] = uo.z;
                r.ke_out[i] = ke_from_u(mo, uo);
                if (pr.out_dest[i] == PFC_DEST_SURVIVOR) survivor = (int)i;
            } else {
                r.out_species[i] = 0xFFFFFFFFu;
                r.out_dest[i] = 0u;
                r.out_u[i][0] = r.out_u[i][1] = r.out_u[i][2] = 0.0;
                r.ke_out[i] = 0.0;
            }
        }
        out->n_events++;
        ++e;
        if (survivor < 0) {
            out->status = PFC_ST_CONSUMED;
            out->species_final = 0xFFFFFFFFu;
            break;
        }
        u = v3(r.out_u[survivor][0], r.out_u[survivor][1], r.out_u[survivor][2]);
    }
    out->n_candidates = e;
    out->u_final[0] = u.x;
    out->u_final[1] = u.y;
    out->u_final[2] = u.z;
}

}  // namespace pfc
