// resident_kernels.cu -- NVRTC device code of the Hall v3.2 device-resident field step
// (v3.1 kernels; v3.2: k_continuity keeps the per-node continuity defect words, k_charge
// computes the candidate field charge, k_gauss_increment accumulates the admitted Gauss
// ledger L += (delta - c), k_admit_nodes measures the absolute Gauss defect relative to L).
// Compiled by hall_gpu::Session::compile_and_load with --fmad=false --prec-div=true
// --prec-sqrt=true --ftz=false for the actual device. No atomics; every sum has a fixed
// order (per-element gathers in the ascending order of the sealed host scatter loops, or
// the fixed two-stage reduction tree below).
//
// Physical arithmetic (bitwise equal to Hall v3 Model::prepare, field_contract.hpp Words):
//   k_operator     conduction S (cells ascending), a = (.5 dt)(S+Y), Jacobi diagonal
//   k_face_rhs     curl words of E, estimated B, H B words
//   k_rhs          ampere words (faces ascending) and the compensated right-hand side
//   k_update_E     E words + increment, compensated midpoint E
//   k_update_B     curl(mid E), B words, compensated midpoint B, H mid B
//   k_currents     ampere(mid B), displacement/solid/exterior/metal/total/short currents
//   k_ledger       metal/exterior incidence words (edges ascending), conductor/exterior
//                  ledger words, bank transfer, KCL node sums, free charge copy
//   k_cells        cell Joule power (12 edges in cell order)
// Guard arithmetic (double-word; Hall v3 uses x87 long double): absolute Gauss, divB,
// continuity, Gauss increment, energy sums. Values are guard quantities, not state.
#include "resident_abi.h"

typedef rf_u64 u64;
typedef rf_u32 u32;
#define TPB RF_TPB
#define RBLK RF_RBLK
#define EPS 2.220446049250313e-16

struct W { double hi, lo; };
__device__ __forceinline__ W two_sum(double a, double b) {
    const double s = a + b, v = s - a;
    W r; r.hi = s; r.lo = (a - (s - v)) + (b - v); return r;
}
// Words(double): hi = value, lo = 0 (no normalization)
__device__ __forceinline__ W wd(double a) { W r; r.hi = a; r.lo = 0.0; return r; }
__device__ __forceinline__ W wpair(double hi, double lo) { W r; r.hi = hi; r.lo = lo; return r; }
// operator+(Words, Words)
__device__ __forceinline__ W wadd(W a, W b) { const W f = two_sum(a.hi, b.hi); return two_sum(f.hi, f.lo + a.lo + b.lo); }
// operator-(Words): Words(-hi, -lo) (two_sum constructor)
__device__ __forceinline__ W wneg(W a) { return two_sum(-a.hi, -a.lo); }
__device__ __forceinline__ W wsub(W a, W b) { return wadd(a, wneg(b)); }
// operator*(Words, double)
__device__ __forceinline__ W wmul(W a, double b) { const double p = a.hi * b; return two_sum(p, fma(a.hi, b, -p) + a.lo * b); }
// operator/(Words, double)
__device__ __forceinline__ W wdiv(W a, double b) {
    const double q = a.hi / b; const W r = wsub(a, wmul(wd(q), b)); return two_sum(q, (r.hi + r.lo) / b);
}
__device__ __forceinline__ double wround(W a) { return a.hi + a.lo; }
__device__ __forceinline__ bool wcanonical(W a) {
    const W n = two_sum(a.hi, a.lo);
    return isfinite(a.hi) && isfinite(a.lo) && n.hi == a.hi && n.lo == a.lo;
}
__device__ __forceinline__ bool wzero(W a) { return a.hi == 0.0 && a.lo == 0.0; }
// Guard arithmetic: double-word values (|x| as a pair, accumulation, comparison).
__device__ __forceinline__ W gabs(W a) { return (a.hi < 0.0 || (a.hi == 0.0 && a.lo < 0.0)) ? wpair(-a.hi, -a.lo) : a; }
__device__ __forceinline__ W gmag(W a) { return wadd(wd(fabs(a.hi)), wd(fabs(a.lo))); }   // |hi| + |lo|
__device__ __forceinline__ bool gle(W a, W b) {   // a <= b for double-word values
    const W d = wsub(b, a); return d.hi > 0.0 || (d.hi == 0.0 && d.lo >= 0.0);
}
__device__ __forceinline__ double gdbl(W a) { return a.hi + a.lo; }
__device__ __forceinline__ W gmul(W a, double b) { return wmul(a, b); }

#define ROWS1(i, n) for (int i = blockIdx.x * TPB + threadIdx.x; i < (n); i += gridDim.x * TPB)
#define PD(p) ((double*)(p))
#define PCD(p) ((const double*)(p))
#define PI(p) ((const int*)(p))
#define PU(p) ((u32*)(p))
#define ET(c, k) (PD((c).eterm) + (size_t)(k) * (c).ne)
#define FT(c, k) (PD((c).fterm) + (size_t)(k) * (c).nb)
#define NV(c, k) (PD((c).nval) + (size_t)(k) * (c).nn)
#define CV(c, k) (PD((c).cval) + (size_t)(k) * (c).ncell)

// ------------------------------------------------------------------ operator
extern "C" __global__ void k_operator(RfCtx c) {
    ROWS1(e, c.ne) {
        const int* ecp = PI(c.ecp); const int* ecc = PI(c.ecc); const double* qv = PCD(c.qv); const double* sigma = PCD(c.sigma);
        double s = 0.0;
        for (int j = ecp[e]; j < ecp[e + 1]; ++j) { const int k = ecc[j]; s += qv[k] * sigma[k]; }
        const double Y = PCD(c.Y)[e], M = PCD(c.M)[e];
        PD(c.S)[e] = s;
        PD(c.a)[e] = (.5 * c.dt) * (s + Y);
        double d = 1.0;
        if (!((const unsigned char*)c.fixed)[e]) {
            d += ((.5 * c.dt) * (s + Y)) / M;
            const int* ctp = PI(c.ctp); const int* ctf = PI(c.ctf); const double* pc2 = PCD(c.pc2); const double* H = PCD(c.H);
            for (int j = ctp[e]; j < ctp[e + 1]; ++j) d += ((((.25 * c.dt) * c.dt) * H[ctf[j]]) * pc2[j]) / M;
        }
        PD(c.diag)[e] = d;
    }
}

// ------------------------------------------------------------------ right-hand side
extern "C" __global__ void k_face_rhs(RfCtx c) {
    ROWS1(b, c.nb) {
        const int* fe = PI(c.fe); const double* fc = PCD(c.fc);
        W sum = wpair(0.0, 0.0);
        for (int k = 0; k < 4; ++k) {
            const int e = fe[4 * b + k];
            sum = wadd(sum, wmul(wpair(PCD(c.old_.E)[e], PCD(c.old_.El)[e]), fc[4 * b + k]));
        }
        const double curl = wround(sum);
        const W est = wsub(wpair(PCD(c.old_.B)[b], PCD(c.old_.Bl)[b]), wmul(wd(curl), .5 * c.dt));
        const W HB = wmul(est, PCD(c.H)[b]);
        PD(c.HBh)[b] = HB.hi; PD(c.HBl)[b] = HB.lo;
    }
}
extern "C" __global__ void k_rhs(RfCtx c) {
    ROWS1(e, c.ne) {
        const int* ctp = PI(c.ctp); const int* ctf = PI(c.ctf); const double* ctc = PCD(c.ctc);
        W r = wpair(0.0, 0.0);
        for (int j = ctp[e]; j < ctp[e + 1]; ++j) { const int b = ctf[j]; r = wadd(r, wmul(wpair(PCD(c.HBh)[b], PCD(c.HBl)[b]), ctc[j])); }
        const double amp = wround(r);
        double out = 0.0;
        if (!((const unsigned char*)c.fixed)[e]) {
            const W oldE = wpair(PCD(c.old_.E)[e], PCD(c.old_.El)[e]);
            W t = wsub(wd(amp), wd(PCD(c.plasma)[e]));
            t = wsub(t, wmul(wd(c.current), PCD(c.feed)[e]));
            t = wsub(t, wmul(oldE, PCD(c.S)[e] + PCD(c.Y)[e]));
            t = wadd(t, wd(PCD(c.bg)[e]));
            t = wmul(t, c.dt);
            t = wdiv(t, PCD(c.sqrtM)[e]);
            out = wround(t);
        }
        PD(c.rhs)[e] = out;
    }
}

// ------------------------------------------------------------------ update
// Edge/face/node/cell guard terms: see resident_abi.h (ET_*, FT_*, NV_*, CV_*).
extern "C" __global__ void k_update_E(RfCtx c) {
    ROWS1(e, c.ne) {
        const W oldE = wpair(PCD(c.old_.E)[e], PCD(c.old_.El)[e]);
        const double inc = PCD(c.x)[e] / PCD(c.sqrtM)[e];
        PD(c.inc0)[e] = inc;
        const W nextE = wadd(oldE, wd(inc));
        PD(c.new_.E)[e] = nextE.hi; PD(c.new_.El)[e] = nextE.lo;
        PD(c.midE)[e] = wround(wmul(wadd(oldE, nextE), .5));
        // field-accuracy indicator (diagnostic): increment history in field units
        const double M = PCD(c.M)[e];
        const double d0 = inc;
        double third = 0.0;
        if (c.history >= 2) third = d0 - 2.0 * PCD(c.inc1)[e] + PCD(c.inc2)[e];
        ET(c, ET_MDE2)[e] = M * d0 * d0;
        ET(c, ET_MTHIRD2)[e] = M * third * third;
        const double En = wround(nextE);
        ET(c, ET_ME2)[e] = M * En * En;
    }
}
extern "C" __global__ void k_update_B(RfCtx c) {
    ROWS1(b, c.nb) {
        const int* fe = PI(c.fe); const double* fc = PCD(c.fc);
        double s = 0.0;
        for (int k = 0; k < 4; ++k) s += fc[4 * b + k] * PCD(c.midE)[fe[4 * b + k]];
        PD(c.cmid)[b] = s;
        const W oldB = wpair(PCD(c.old_.B)[b], PCD(c.old_.Bl)[b]);
        const W nextB = wsub(oldB, wmul(wd(s), c.dt));
        PD(c.new_.B)[b] = nextB.hi; PD(c.new_.Bl)[b] = nextB.lo;
        const double mB = wround(wmul(wadd(oldB, nextB), .5));
        PD(c.midB)[b] = mB;
        const double H = PCD(c.H)[b];
        PD(c.HmidB)[b] = mB * H;
        // energy change term of B (guard): H * ((nextB + oldB) * .5) * (nextB - oldB)
        const W tm = wmul(wadd(nextB, oldB), .5), dl = wsub(nextB, oldB);
        const W p = wmul(tm, H);
        FT(c, FT_FDB)[b] = gdbl(wadd(wmul(p, dl.hi), wd(p.hi * dl.lo)));
    }
}
extern "C" __global__ void k_currents(RfCtx c) {
    ROWS1(e, c.ne) {
        const int* ctp = PI(c.ctp); const int* ctf = PI(c.ctf); const double* ctc = PCD(c.ctc);
        double amp = 0.0;
        for (int j = ctp[e]; j < ctp[e + 1]; ++j) amp += ctc[j] * PCD(c.HmidB)[ctf[j]];
        const W oldE = wpair(PCD(c.old_.E)[e], PCD(c.old_.El)[e]);
        const W nextE = wpair(PCD(c.new_.E)[e], PCD(c.new_.El)[e]);
        const W delta = wsub(nextE, oldE);
        const W tmid = wmul(wadd(nextE, oldE), .5);
        const double M = PCD(c.M)[e], Y = PCD(c.Y)[e], bg = PCD(c.bg)[e], feed = PCD(c.feed)[e], plasma = PCD(c.plasma)[e];
        const double mid = PCD(c.midE)[e];
        const double disp = wround(wmul(delta, M / c.dt));
        const double solid = PCD(c.S)[e] * mid;
        const double ext = Y * mid - bg;
        double metal = solid + c.current * feed;
        double total = metal + plasma + ext + disp;
        const double residual = amp - total;
        double shortc = 0.0, ampres = 0.0;
        if (((const unsigned char*)c.fixed)[e]) { shortc = residual; metal += residual; total += residual; }
        else ampres = fabs(residual);
        PD(c.disp)[e] = disp; PD(c.solid)[e] = solid; PD(c.exterior)[e] = ext; PD(c.metal)[e] = metal; PD(c.total)[e] = total;
        PD(c.short_cur)[e] = shortc;
        // guard terms
        const W p = wmul(tmid, M);
        ET(c, ET_FDE)[e] = gdbl(wadd(wmul(p, delta.hi), wd(p.hi * delta.lo)));
        ET(c, ET_PW)[e] = c.dt * mid * plasma;
        ET(c, ET_JOULE)[e] = c.dt * mid * solid;
        ET(c, ET_CIRCUIT)[e] = c.dt * c.current * feed * mid;
        ET(c, ET_OLOSS)[e] = c.dt * Y * mid * mid;
        ET(c, ET_OSUPPLY)[e] = c.dt * bg * mid;
        ET(c, ET_AMPRES)[e] = isfinite(residual) ? ampres : 1.0 / 0.0;
        ET(c, ET_TERMINAL)[e] = feed * mid;
        const double rg = PCD(c.region_grad)[e];
        ET(c, ET_REG_SOLID)[e] = rg * solid;
        ET(c, ET_REG_PLASMA)[e] = rg * plasma;
        ET(c, ET_REG_DISP)[e] = rg * disp;
        ET(c, ET_REG_EXT)[e] = rg * ext;
    }
}

// ------------------------------------------------------------------ continuity (before the solve)
// code_cont: endpoint continuity mismatch; NV_CONT = |defect|/tolerance (ratio, guard).
extern "C" __global__ void k_continuity(RfCtx c) {
    ROWS1(n, c.nn) {
        const int* ntp = PI(c.ntp); const int* nte = PI(c.nte); const int* nts = PI(c.nts);
        W dv = wpair(0.0, 0.0);
        for (int j = ntp[n]; j < ntp[n + 1]; ++j) {
            const int e = nte[j];
            const W q = wdiv(wd(PCD(c.plasma)[e]), PCD(c.elen)[e]);
            dv = nts[j] < 0 ? wsub(dv, q) : wadd(dv, q);
        }
        PD(c.div_hi)[n] = dv.hi; PD(c.div_lo)[n] = dv.lo;
        const W after = wpair(PCD(c.after_hi)[n], PCD(c.after_lo)[n]);
        const W before = wpair(PCD(c.old_.fnet_hi)[n], PCD(c.old_.fnet_lo)[n]);
        const W ddt = wmul(dv, c.dt);
        W defect;
        if (c.transfer)
            defect = wsub(wsub(wadd(wadd(after, wpair(PCD(c.bank_cond_hi)[n], PCD(c.bank_cond_lo)[n])), wpair(PCD(c.bank_ext_hi)[n], PCD(c.bank_ext_lo)[n])), before), ddt);
        else defect = wsub(wsub(after, before), ddt);
        W scale = wadd(wadd(wd(PCD(c.old_.fgross)[n]), wd(PCD(c.after_gross)[n])), gmag(ddt));
        if (c.transfer) scale = wadd(scale, wd(PCD(c.bank_gross)[n]));
        const W tol = wadd(wmul(scale, 1e-12), wd(PCD(c.support)[n]));
        const W size = gabs(defect);
        const bool ok = gle(size, tol);
        double ratio = 0.0;
        if (gdbl(size) > 0.0) ratio = gdbl(size) / gdbl(tol);
        PU(c.code_cont)[n] = ok ? 0u : 1u;
        NV(c, NV_CONT)[n] = isfinite(ratio) ? ratio : (ok ? 0.0 : 1.0 / 0.0);
        PD(c.cdef_hi)[n] = defect.hi; PD(c.cdef_lo)[n] = defect.lo;   // v3.2: c_n for the admitted Gauss ledger
    }
}

// ------------------------------------------------------------------ admission (absolute Gauss)
// Per node of slot s (0 old, 1 new): field charge words Q (stored), absolute Gauss defect and
// bound, support masks, canonical owned charge and the free-charge operand guard.
// Codes (first failing node, checks in the sealed order): 2 noncanonical owned charge,
// 3 conductor support, 4 exterior support, 5 absolute Gauss; free_out: free charge operand guard.
__device__ __forceinline__ const RfSlot& slot_of(const RfCtx& c, int s) { return s ? c.new_ : c.old_; }
// Field charge words of node n (field_charge: incidence of E words * (-M), edges ascending).
__device__ __forceinline__ W node_charge(const RfCtx& c, const RfSlot& z, int n) {
    const int* ntp = PI(c.ntp); const int* nte = PI(c.nte); const int* nts = PI(c.nts);
    W Q = wpair(0.0, 0.0);
    for (int j = ntp[n]; j < ntp[n + 1]; ++j) {
        const int e = nte[j];
        const W q = wdiv(wmul(wpair(PCD(z.E)[e], PCD(z.El)[e]), -PCD(c.M)[e]), PCD(c.elen)[e]);
        Q = nts[j] < 0 ? wsub(Q, q) : wadd(Q, q);
    }
    return Q;
}
// v3.2: field charge words of slot s (the candidate charge Q1 before the admitted Gauss ledger update).
extern "C" __global__ void k_charge(RfCtx c, int s, u64 Qh, u64 Ql) {
    const RfSlot& z = slot_of(c, s);
    ROWS1(n, c.nn) { const W Q = node_charge(c, z, n); PD(Qh)[n] = Q.hi; PD(Ql)[n] = Q.lo; }
}
// v3.2 outputs: defect_out |D - L| (the checked value), raw_out |D|, led_out |L|, rat_out |D - L| / bound,
// apr_out |D| / (bound + |L|);
// code 6: noncanonical admitted Gauss ledger (checked after the support codes, before the Gauss code 5).
extern "C" __global__ void k_admit_nodes(RfCtx c, int s, u64 Qh, u64 Ql, u64 code_out, u64 free_out, u64 defect_out, u64 bound_out, u64 act_out,
                                         u64 raw_out, u64 led_out, u64 rat_out, u64 apr_out) {
    const RfSlot& z = slot_of(c, s);
    ROWS1(n, c.nn) {
        const int* ntp = PI(c.ntp); const int* nte = PI(c.nte);
        W activity = wpair(0.0, 0.0);
        for (int j = ntp[n]; j < ntp[n + 1]; ++j) {
            const int e = nte[j];
            const W E = wpair(PCD(z.E)[e], PCD(z.El)[e]);
            const double M = PCD(c.M)[e], len = PCD(c.elen)[e];
            activity = wadd(activity, wmul(wmul(gmag(E), M), 1.0 / len));
        }
        const W Q = node_charge(c, z, n);
        PD(Qh)[n] = Q.hi; PD(Ql)[n] = Q.lo;
        const W net = wpair(PCD(z.fnet_hi)[n], PCD(z.fnet_lo)[n]);
        const double gross = PCD(z.fgross)[n];
        const W cond = wpair(PCD(z.cond_hi)[n], PCD(z.cond_lo)[n]), ext = wpair(PCD(z.ext_hi)[n], PCD(z.ext_lo)[n]);
        u32 code = 0;
        // free charge operand guard (runs over all nodes before the per-node loop in the sealed order)
        const bool free_ok = wcanonical(net) && isfinite(gross) && gross >= 0.0 && gle(gabs(net), wmul(wd(gross), 1.0 + 32.0 * EPS));
        PU(free_out)[n] = free_ok ? 0u : 1u;
        const double residue = PCD(z.residue)[n];
        if (!wcanonical(cond) || !wcanonical(ext)) code = 2;
        const bool cond_supported = ((const unsigned char*)c.node_body_or_short)[n] != 0;
        const bool ext_supported = ((const unsigned char*)c.exterior_support)[n] != 0;
        if (!code && !(cond_supported || wzero(cond) || gle(gmag(cond), wd(residue)))) code = 3;
        if (!code && !(ext_supported || wzero(ext) || gle(gmag(ext), wd(residue)))) code = 4;
        const W ledger = wpair(PCD(z.gl_hi)[n], PCD(z.gl_lo)[n]);
        if (!code && !wcanonical(ledger)) code = 6;
        const W defect = wsub(wsub(wsub(Q, net), cond), ext);
        const W centered = wsub(defect, ledger);
        const W bsum = wadd(wadd(wadd(activity, wd(gross)), gmag(cond)), gmag(ext));
        const W bound = wadd(wd(c.initial_gauss_bound), wmul(bsum, 256.0 * EPS));
        if (!code && !gle(gabs(centered), bound)) code = 5;
        PU(code_out)[n] = code;
        const double size = gdbl(gabs(centered)), bnd = gdbl(bound);
        PD(defect_out)[n] = size;
        PD(bound_out)[n] = bnd;
        PD(act_out)[n] = gdbl(activity);
        PD(raw_out)[n] = gdbl(gabs(defect));
        PD(led_out)[n] = gdbl(gabs(ledger));
        const double ratio = bnd > 0.0 ? size / bnd : (size > 0.0 ? 1.0 / 0.0 : 0.0);
        PD(rat_out)[n] = isfinite(ratio) ? ratio : 1.0 / 0.0;
        const double raw = gdbl(gabs(defect)), wide = gdbl(wadd(bound, gabs(ledger)));
        const double implied = wide > 0.0 ? raw / wide : (raw > 0.0 ? 1.0 / 0.0 : 0.0);
        PD(apr_out)[n] = isfinite(implied) ? implied : 1.0 / 0.0;
    }
}
// Per edge: canonical E words, zero PEC tangential field; per face: canonical B words (u32 flags).
extern "C" __global__ void k_admit_edges(RfCtx c, int s, u64 edge_canon, u64 edge_pec, u64 face_code) {
    const RfSlot& z = slot_of(c, s);
    ROWS1(e, c.ne) {
        const W E = wpair(PCD(z.E)[e], PCD(z.El)[e]);
        PU(edge_canon)[e] = wcanonical(E) ? 0u : 1u;
        PU(edge_pec)[e] = (((const unsigned char*)c.fixed)[e] && !wzero(E)) ? 1u : 0u;
    }
    ROWS1(b, c.nb) {
        const W B = wpair(PCD(z.B)[b], PCD(z.Bl)[b]);
        PU(face_code)[b] = wcanonical(B) ? 0u : 1u;
    }
}
// Per cell: absolute divB (guard) of slot s; cell_v0 |divB|, cell_v1 bound; code 1 on failure.
__device__ __forceinline__ int bx_(const RfCtx& c, int i, int j, int k) { return (k * c.ny + j) * (c.nx + 1) + i; }
__device__ __forceinline__ int by_(const RfCtx& c, int i, int j, int k) { return (c.nx + 1) * c.ny * c.nz + (k * (c.ny + 1) + j) * c.nx + i; }
__device__ __forceinline__ int bz_(const RfCtx& c, int i, int j, int k) { return (c.nx + 1) * c.ny * c.nz + c.nx * (c.ny + 1) * c.nz + (k * c.ny + j) * c.nx + i; }
extern "C" __global__ void k_admit_cells(RfCtx c, int s, u64 code_out, u64 value_out, u64 bound_out) {
    const RfSlot& z = slot_of(c, s);
    ROWS1(q, c.ncell) {
        const int i = q % c.nx, j = (q / c.nx) % c.ny, k = q / (c.nx * c.ny);
        const int fa[3] = {bx_(c, i + 1, j, k), by_(c, i, j + 1, k), bz_(c, i, j, k + 1)};
        const int fb[3] = {bx_(c, i, j, k), by_(c, i, j, k), bz_(c, i, j, k)};
        const double h[3] = {c.dx, c.dy, c.dz};
        W d = wpair(0.0, 0.0), scale = wpair(0.0, 0.0);
        for (int a = 0; a < 3; ++a) {
            const W Ba = wpair(PCD(z.B)[fa[a]], PCD(z.Bl)[fa[a]]), Bb = wpair(PCD(z.B)[fb[a]], PCD(z.Bl)[fb[a]]);
            d = wadd(d, wdiv(wsub(Ba, Bb), h[a]));
            scale = wadd(scale, wmul(wadd(gmag(Ba), gmag(Bb)), 1.0 / h[a]));
        }
        const W bound = wmul(scale, 128.0 * EPS);
        const bool ok = gle(gabs(d), bound);
        PU(code_out)[q] = ok ? 0u : 1u;
        PD(value_out)[q] = gdbl(gabs(d));
        PD(bound_out)[q] = gdbl(bound);
    }
}

// ------------------------------------------------------------------ ledger and node sums
// x87 extended emulation (non-negative operands) for the bank residue bound, as in
// hall-gpu-resident-v1 push_core.hpp (x80_add, x80_to_double with the upward step).
struct X80 { u64 m; int e; };
__device__ X80 x80_from_double(double d) {
    const u64 bits = (u64)__double_as_longlong(d) & 0x7fffffffffffffffull;
    const int exponent = int(bits >> 52);
    const u64 fraction = bits & 0x000fffffffffffffull;
    X80 r;
    if (exponent == 0) {
        if (fraction == 0) { r.m = 0; r.e = 0; return r; }
        const int shift = __clzll((long long)fraction);
        r.m = fraction << shift; r.e = -1074 - shift; return r;
    }
    r.m = (fraction | 0x0010000000000000ull) << 11; r.e = exponent - 1075 - 11; return r;
}
__device__ X80 x80_add(X80 a, X80 b) {
    if (a.m == 0) return b;
    if (b.m == 0) return a;
    if (a.e < b.e) { const X80 t = a; a = b; b = t; }
    const int d = a.e - b.e;
    u64 bhi, blo; bool sticky = false;
    if (d == 0) { bhi = b.m; blo = 0; }
    else if (d < 64) { bhi = b.m >> d; blo = b.m << (64 - d); }
    else if (d == 64) { bhi = 0; blo = b.m; }
    else if (d < 128) { bhi = 0; blo = b.m >> (d - 64); sticky = (b.m << (128 - d)) != 0; }
    else { bhi = 0; blo = 0; sticky = true; }
    const u64 lo = blo, hi = a.m + bhi; const bool carry = hi < a.m;
    X80 r; bool guard, rest; u64 m;
    if (carry) { m = 0x8000000000000000ull | (hi >> 1); guard = (hi & 1) != 0; rest = lo != 0 || sticky; r.e = a.e + 1; }
    else { m = hi; guard = (lo >> 63) != 0; rest = (lo << 1) != 0 || sticky; r.e = a.e; }
    if (guard && (rest || (m & 1))) { m += 1; if (m == 0) { m = 0x8000000000000000ull; r.e += 1; } }
    r.m = m; return r;
}
__device__ double x80_to_double_up(X80 x) {
    if (x.m == 0) return 0.0;
    const int top = x.e + 63;
    u64 q, rem, half, high_part;
    if (top >= -1022) { q = x.m >> 11; rem = x.m & 0x7ffull; half = 0x400ull; high_part = u64(top + 1022) << 52; }
    else {
        const int shift = -(x.e + 1074); high_part = 0;
        if (shift > 64) { q = 0; rem = x.m; half = 0; }
        else if (shift == 64) { q = 0; rem = x.m; half = 0x8000000000000000ull; }
        else { q = x.m >> shift; rem = x.m & ((u64(1) << shift) - 1); half = u64(1) << (shift - 1); }
    }
    bool up;
    if (half == 0) up = false; else up = rem > half || (rem == half && (q & 1));
    u64 bits = high_part + q + (up ? 1 : 0);
    if (bits >= 0x7ff0000000000000ull) bits = 0x7ff0000000000000ull;
    if (rem != 0 && !up && bits < 0x7ff0000000000000ull) bits += 1;
    return __longlong_as_double((long long)bits);
}
extern "C" __global__ void k_ledger(RfCtx c, u64 md_hi, u64 md_lo, u64 od_hi, u64 od_lo) {
    ROWS1(n, c.nn) {
        const int* ntp = PI(c.ntp); const int* nte = PI(c.nte); const int* nts = PI(c.nts);
        W md = wpair(0.0, 0.0), od = wpair(0.0, 0.0);
        double kcl = 0.0;
        for (int j = ntp[n]; j < ntp[n + 1]; ++j) {
            const int e = nte[j]; const double len = PCD(c.elen)[e];
            const W qm = wdiv(wd(PCD(c.metal)[e]), len), qo = wdiv(wd(PCD(c.exterior)[e]), len);
            const double qt = PCD(c.total)[e] / len;
            if (nts[j] < 0) { md = wsub(md, qm); od = wsub(od, qo); kcl -= qt; }
            else { md = wadd(md, qm); od = wadd(od, qo); kcl += qt; }
        }
        PD(md_hi)[n] = md.hi; PD(md_lo)[n] = md.lo; PD(od_hi)[n] = od.hi; PD(od_lo)[n] = od.lo;
        NV(c, NV_KCL)[n] = isfinite(kcl) ? fabs(kcl) : 1.0 / 0.0;
        W cond = wpair(PCD(c.old_.cond_hi)[n], PCD(c.old_.cond_lo)[n]);
        W ext = wpair(PCD(c.old_.ext_hi)[n], PCD(c.old_.ext_lo)[n]);
        double residue = PCD(c.old_.residue)[n];
        if (c.transfer) {
            const W bc = wpair(PCD(c.bank_cond_hi)[n], PCD(c.bank_cond_lo)[n]), be = wpair(PCD(c.bank_ext_hi)[n], PCD(c.bank_ext_lo)[n]);
            if (bc.hi != 0.0 || bc.lo != 0.0) cond = wadd(cond, bc);
            if (be.hi != 0.0 || be.lo != 0.0) ext = wadd(ext, be);
            if (c.residue_any) {
                const double r = PCD(c.bank_residue)[n];
                if (r != 0.0) residue = x80_to_double_up(x80_add(x80_from_double(residue), x80_from_double(r)));
            }
        }
        cond = wadd(cond, wmul(md, c.dt));
        ext = wadd(ext, wmul(od, c.dt));
        PD(c.new_.cond_hi)[n] = cond.hi; PD(c.new_.cond_lo)[n] = cond.lo;
        PD(c.new_.ext_hi)[n] = ext.hi; PD(c.new_.ext_lo)[n] = ext.lo;
        PD(c.new_.residue)[n] = residue;
        PD(c.new_.fnet_hi)[n] = PCD(c.after_hi)[n]; PD(c.new_.fnet_lo)[n] = PCD(c.after_lo)[n]; PD(c.new_.fgross)[n] = PCD(c.after_gross)[n];
    }
}
// Gauss increment (guard): |Q1 - Q0 - (div + metal + outer) dt| and |Q1 - Q0| per node.
// v3.2: admitted Gauss ledger of the candidate, L1 = L0 + (delta - c) (state words; the operations of
// Model::prepare: ledger[n] + (delta - continuity[n])).
extern "C" __global__ void k_gauss_increment(RfCtx c, u64 md_hi, u64 md_lo, u64 od_hi, u64 od_lo) {
    ROWS1(n, c.nn) {
        const W change = wsub(wpair(PCD(c.Q1_hi)[n], PCD(c.Q1_lo)[n]), wpair(PCD(c.Q0_hi)[n], PCD(c.Q0_lo)[n]));
        const W expected = wmul(wadd(wadd(wpair(PCD(c.div_hi)[n], PCD(c.div_lo)[n]), wpair(PCD(md_hi)[n], PCD(md_lo)[n])),
                                     wpair(PCD(od_hi)[n], PCD(od_lo)[n])), c.dt);
        const W delta = wsub(change, expected);
        const double v = gdbl(gabs(delta));
        const double a = gdbl(gmag(change));
        NV(c, NV_GINC)[n] = isfinite(v) ? v : 1.0 / 0.0;
        NV(c, NV_GACT)[n] = isfinite(a) ? a : 0.0;
        const W L = wadd(wpair(PCD(c.old_.gl_hi)[n], PCD(c.old_.gl_lo)[n]), wsub(delta, wpair(PCD(c.cdef_hi)[n], PCD(c.cdef_lo)[n])));
        PD(c.new_.gl_hi)[n] = L.hi; PD(c.new_.gl_lo)[n] = L.lo;
    }
}
// Cell Joule power (bitwise) and the controller quantity |B| at cell centres (mid B).
extern "C" __global__ void k_cells(RfCtx c) {
    ROWS1(q, c.ncell) {
        const int* ce = PI(c.cedge);
        const double w = PCD(c.qv)[q] * PCD(c.sigma)[q];
        double p = 0.0;
        for (int k = 0; k < 12; ++k) { const double E = PCD(c.midE)[ce[12 * q + k]]; p += w * E * E; }
        PD(c.joule_cell)[q] = p;
        const int i = q % c.nx, j = (q / c.nx) % c.ny, k = q / (c.nx * c.ny);
        const double* mB = PCD(c.midB);
        const double bx = .5 * (mB[bx_(c, i, j, k)] + mB[bx_(c, i + 1, j, k)]);
        const double by = .5 * (mB[by_(c, i, j, k)] + mB[by_(c, i, j + 1, k)]);
        const double bz = .5 * (mB[bz_(c, i, j, k)] + mB[bz_(c, i, j, k + 1)]);
        const double bm = sqrt(bx * bx + by * by + bz * bz);
        CV(c, CV_BMAG)[q] = isfinite(bm) ? bm : 1.0 / 0.0;
        // vacuum cells must carry exactly +0 Joule power (thermal owner rule)
        const bool body = ((const int*)c.cell_body)[q] >= 0;
        PU(c.code_joule)[q] = (isfinite(p) && p >= 0.0 && (body || p == 0.0)) ? 0u : 1u;
    }
}

// ------------------------------------------------------------------ reductions (fixed tree)
__device__ __forceinline__ W tree_w(W* sh, W v) {
    const int t = threadIdx.x;
    sh[t] = v; __syncthreads();
    for (int s = TPB / 2; s > 0; s >>= 1) { if (t < s) sh[t] = wadd(sh[t], sh[t + s]); __syncthreads(); }
    return sh[0];
}
__device__ __forceinline__ double tree_max(double* sh, double v) {
    const int t = threadIdx.x;
    sh[t] = v; __syncthreads();
    for (int s = TPB / 2; s > 0; s >>= 1) { if (t < s) sh[t] = fmax(sh[t], sh[t + s]); __syncthreads(); }
    return sh[0];
}
__device__ __forceinline__ u32 tree_min(u32* sh, u32 v) {
    const int t = threadIdx.x;
    sh[t] = v; __syncthreads();
    for (int s = TPB / 2; s > 0; s >>= 1) { if (t < s) sh[t] = min(sh[t], sh[t + s]); __syncthreads(); }
    return sh[0];
}
// Stage 1: grid (RBLK, rows). Partials: W per (row, block) at part + (row*RBLK + block)*16.
// kind 2 (min u32) reduces "first index with nonzero code": value = code ? index : 0xffffffff.
extern "C" __global__ void k_reduce1(RfRows r, u64 part) {
    const int row = blockIdx.y;
    if (row >= r.count) return;
    __shared__ W shw[TPB];
    const int n = r.n[row], kind = r.kind[row];
    W* out = (W*)part + (row * RBLK + blockIdx.x);
    if (kind == 0) {
        W acc = wpair(0.0, 0.0);
        for (int i = blockIdx.x * TPB + threadIdx.x; i < n; i += RBLK * TPB) acc = wadd(acc, wd(PCD(r.ptr[row])[i]));
        const W v = tree_w(shw, acc);
        if (threadIdx.x == 0) *out = v;
    } else if (kind == 1) {
        double acc = 0.0;
        for (int i = blockIdx.x * TPB + threadIdx.x; i < n; i += RBLK * TPB) { const double x = PCD(r.ptr[row])[i]; acc = fmax(acc, x != x ? 1.0 / 0.0 : x); }
        const double v = tree_max((double*)shw, acc);
        if (threadIdx.x == 0) *out = wpair(v, 0.0);
    } else {
        u32 acc = 0xffffffffu;
        for (int i = blockIdx.x * TPB + threadIdx.x; i < n; i += RBLK * TPB) if (((const u32*)r.ptr[row])[i] != 0u) acc = min(acc, (u32)i);
        const u32 v = tree_min((u32*)shw, acc);
        if (threadIdx.x == 0) *out = wpair((double)v, 0.0);
    }
}
// Stage 2: one block per row; result W per row at res + row*16.
extern "C" __global__ void k_reduce2(RfRows r, u64 part, u64 res) {
    const int row = blockIdx.x;
    if (row >= r.count) return;
    __shared__ W shw[TPB];
    const int kind = r.kind[row];
    const W* p = (const W*)part + row * RBLK;
    W v = p[threadIdx.x];   // RBLK == TPB
    W out;
    if (kind == 0) out = tree_w(shw, v);
    else if (kind == 1) out = wpair(tree_max((double*)shw, v.hi), 0.0);
    else out = wpair((double)tree_min((u32*)shw, (u32)v.hi), 0.0);
    if (threadIdx.x == 0) ((W*)res)[row] = out;
}

// ------------------------------------------------------------------ comparison
// flags[i] = 1 where a[i] != b[i] as doubles (std::vector<double> equality semantics).
extern "C" __global__ void k_differ(int n, u64 a, u64 b, u64 flags) {
    ROWS1(i, n) PU(flags)[i] = (PCD(a)[i] == PCD(b)[i]) ? 0u : 1u;
}
// Midpoint recomputation (preflight): mid E/B from the accepted and candidate words.
extern "C" __global__ void k_midpoint(RfCtx c, u64 outE, u64 outB) {
    ROWS1(e, c.ne) PD(outE)[e] = wround(wmul(wadd(wpair(PCD(c.old_.E)[e], PCD(c.old_.El)[e]), wpair(PCD(c.new_.E)[e], PCD(c.new_.El)[e])), .5));
    ROWS1(b, c.nb) PD(outB)[b] = wround(wmul(wadd(wpair(PCD(c.old_.B)[b], PCD(c.old_.Bl)[b]), wpair(PCD(c.new_.B)[b], PCD(c.new_.Bl)[b])), .5));
}
// Rounded accepted words (Model::at): E = hi + lo, B = hi + lo.
extern "C" __global__ void k_rounded(RfCtx c, u64 outE, u64 outB) {
    ROWS1(e, c.ne) PD(outE)[e] = PCD(c.old_.E)[e] + PCD(c.old_.El)[e];
    ROWS1(b, c.nb) PD(outB)[b] = PCD(c.old_.B)[b] + PCD(c.old_.Bl)[b];
}
extern "C" __global__ void k_copy_f64(int n, u64 src, u64 dst) { ROWS1(i, n) PD(dst)[i] = PCD(src)[i]; }
// Diagnostic accumulation: out = acc + add (per element).
extern "C" __global__ void k_add_f64(int n, u64 acc, u64 add, u64 out) { ROWS1(i, n) PD(out)[i] = PCD(acc)[i] + PCD(add)[i]; }
extern "C" __global__ void k_gather_f64(int n, u64 idx, u64 src, u64 dst) { ROWS1(i, n) PD(dst)[i] = PCD(src)[((const int*)idx)[i]]; }
// Controller: omega_p dt from the free gross charge of slot new (node volumes VN, electron mass).
extern "C" __global__ void k_omega_p(RfCtx c, u64 VN, double factor) {
    ROWS1(n, c.nn) {
        const double g = PCD(c.new_.fgross)[n];
        const double v = sqrt(factor * g / PCD(VN)[n]) * c.dt;   // factor = e / (eps0 m_e)
        NV(c, NV_OMEGAP)[n] = isfinite(v) ? v : 1.0 / 0.0;
    }
}

// ------------------------------------------------------------------ SHA-256 leaves
__constant__ u32 K256[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
__device__ __forceinline__ u32 rotr(u32 x, int n) { return (x >> n) | (x << (32 - n)); }
__device__ __forceinline__ u32 bswap(u32 x) { return (x >> 24) | ((x >> 8) & 0xff00u) | ((x << 8) & 0xff0000u) | (x << 24); }
__device__ void sha_compress(u32* h, const u32* w16) {
    u32 w[64];
    for (int i = 0; i < 16; ++i) w[i] = w16[i];
    for (int i = 16; i < 64; ++i) {
        const u32 s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
        const u32 s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    u32 a = h[0], b = h[1], cc = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
    for (int i = 0; i < 64; ++i) {
        const u32 S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        const u32 ch = (e & f) ^ (~e & g);
        const u32 t1 = hh + S1 + ch + K256[i] + w[i];
        const u32 S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        const u32 mj = (a & b) ^ (a & cc) ^ (b & cc);
        const u32 t2 = S0 + mj;
        hh = g; g = f; f = e; e = d + t1; d = cc; cc = b; b = a; a = t1 + t2;
    }
    h[0] += a; h[1] += b; h[2] += cc; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
}
// Leaf L of the table: SHA-256 of u64 tag || u64 block index || u64 byte count || payload
// (all u64 little-endian; payload = the binary64/u64 words of the block as stored).
extern "C" __global__ void k_leaves(RfLeafTable t) {
    const int L = blockIdx.x * blockDim.x + threadIdx.x;
    if (L >= t.first_leaf[t.arrays]) return;
    int a = 0;
    while (L >= t.first_leaf[a + 1]) ++a;
    const u64 k = (u64)(L - t.first_leaf[a]);
    const u64 total = t.bytes[a];
    const u64 begin = k * RF_BLOCK_BYTES;
    const u64 nbytes = (total - begin < (u64)RF_BLOCK_BYTES) ? total - begin : (u64)RF_BLOCK_BYTES;
    const u64* data = (const u64*)(t.ptr[a] + begin);
    const u64 words = 3 + nbytes / 8;            // message length in u64 words (multiple of 8 bytes)
    const u64 bitlen = words * 64;
    const u64 total_words = ((words + 1 + 1) + 7) / 8 * 8;   // + padding word + length word, rounded to 64 bytes
    u32 h[8] = {0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au, 0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
    u32 w[16];
    for (u64 base = 0; base < total_words; base += 8) {
        for (int j = 0; j < 8; ++j) {
            const u64 idx = base + j;
            u64 v;
            bool raw = false;   // raw: already SHA word order (length word)
            if (idx == 0) v = t.tag[a];
            else if (idx == 1) v = k;
            else if (idx == 2) v = nbytes;
            else if (idx < words) v = data[idx - 3];
            else if (idx == words) v = 0x80ull;
            else if (idx == total_words - 1) { v = bitlen; raw = true; }
            else v = 0ull;
            if (raw) { w[2 * j] = (u32)(v >> 32); w[2 * j + 1] = (u32)(v & 0xffffffffull); }
            else { w[2 * j] = bswap((u32)(v & 0xffffffffull)); w[2 * j + 1] = bswap((u32)(v >> 32)); }
        }
        sha_compress(h, w);
    }
    unsigned char* out = (unsigned char*)(t.out + 32ull * L);
    for (int i = 0; i < 8; ++i) { out[4 * i] = (unsigned char)(h[i] >> 24); out[4 * i + 1] = (unsigned char)(h[i] >> 16); out[4 * i + 2] = (unsigned char)(h[i] >> 8); out[4 * i + 3] = (unsigned char)h[i]; }
}
// Copies the payload blocks of sampled leaves into staging (RF_BLOCK_BYTES per sample).
extern "C" __global__ void k_gather_blocks(RfLeafTable t, int samples, u64 leaf_ids, u64 staging) {
    const int s = blockIdx.x;
    if (s >= samples) return;
    const int L = ((const int*)leaf_ids)[s];
    int a = 0;
    while (L >= t.first_leaf[a + 1]) ++a;
    const u64 k = (u64)(L - t.first_leaf[a]);
    const u64 begin = k * RF_BLOCK_BYTES, total = t.bytes[a];
    const u64 nbytes = (total - begin < (u64)RF_BLOCK_BYTES) ? total - begin : (u64)RF_BLOCK_BYTES;
    const u64* src = (const u64*)(t.ptr[a] + begin);
    u64* dst = (u64*)(staging + (u64)s * RF_BLOCK_BYTES);
    for (u64 i = threadIdx.x; i < nbytes / 8; i += blockDim.x) dst[i] = src[i];
}
