// step_kernels_v34.cu -- NVRTC program of the Hall v3.4 device coupled step (device_step_v34.cpp).
// Compiled with --fmad=false --prec-div=true --prec-sqrt=true --ftz=false for the actual device. No atomics:
// exact sums are integer digit sums (exact_sum_v34.cuh), node accumulations walk each node's records after a
// stable sort in the sealed carrier order, failure searches are fixed-tree minima (primitives).
//
//   census     k34_rows_to_soa, k34_find_ids, k34_flights_init, k34_newborn_place, k34_k0, k34_neutral_fail,
//              k34_records_charge (charge_of_particles per carrier record: the sealed checks and Words CIC weights),
//              k34_charge_nodes (net words in carrier order, x87 gross), k34_census_nodes (free_charge_guard,
//              same_net_charge against the accepted free charge, x87 emulated)
//   fold       k34_sum / k34_sum_finish (K, W, G_K, C_K, Q, G_W; norms of the independent residual),
//              k34_continuity, k34_node_marks / k34_edge_bound (E^), k34_free_guard, k34_fail_nonfinite
//   record     k34_assemble_check, k34_record_flags, k34_record_build (endpoints 80 B, paths 64 B, sinks 136 B,
//              births 24 B, canonical layouts of API-007.txt 1.2)
//   support    k34_support_carriers, k34_support_emit, k34_support_nodes (support_v34.hpp, x87 emulated)
//   blocks     k34_f64_chunks / k34_f64_blocks / k34_f64_tree: block lists of binary64 arrays with the native
//              v7.1 field/current tags (API-071.txt 1.1)
#include "forms_core.h"
#include "exact_sum_v34.cuh"
#include "sink_core.hpp"
#include "support_v34.hpp"
#include "step_abi_v34.h"

typedef s34_u64 u64;
typedef s34_u32 u32;
#define TPB S34_TPB
#define IDX() (blockIdx.x * (u64)TPB + threadIdx.x)
#define PD(p) ((double*)(p))
#define PCD(p) ((const double*)(p))
#define PU64(p) ((u64*)(p))
#define PCU64(p) ((const u64*)(p))
#define PU32(p) ((u32*)(p))
#define PCU32(p) ((const u32*)(p))
#define PI(p) ((int*)(p))
#define PCI(p) ((const int*)(p))

// ---------------------------------------------------------------- census
// Census rows (GapParticle, 10 words: id, x[3], u[3], q, m, w) -> push carriers (SoA, stride n) and IDs.
extern "C" __global__ void k34_rows_to_soa(u64 n, u64 rows, u64 ids, u64 soa) {
    const u64 i = IDX();
    if (i >= n) return;
    const u64* r = PCU64(rows) + 10 * i;
    PU64(ids)[i] = r[0];
    double* s = PD(soa);
    for (int k = 0; k < 9; ++k) s[k * n + i] = pf7::double_of(r[1 + k]);
}
// Ascending-ID position of each wanted ID (n when absent).
extern "C" __global__ void k34_find_ids(u64 count, u64 want, u64 ids, u64 n, u64 out) {
    const u64 j = IDX();
    if (j >= count) return;
    const u64 id = PCU64(want)[j];
    u64 lo = 0, hi = n;
    while (lo < hi) { const u64 mid = lo + (hi - lo) / 2; if (PCU64(ids)[mid] < id) lo = mid + 1; else hi = mid; }
    PU64(out)[j] = (lo < n && PCU64(ids)[lo] == id) ? lo : n;
}
// Push flights (flight F, newborn face, newborn flag), the newborn face axis and birth time per carrier.
extern "C" __global__ void k34_flights_init(u64 n, u64 flight, u64 face, u64 newborn, u64 axis, u64 tau, double dt, double begin_s) {
    const u64 i = IDX();
    if (i >= n) return;
    PD(flight)[i] = dt; PI(face)[i] = -1; PI(newborn)[i] = 0; PI(axis)[i] = -1; PD(tau)[i] = begin_s;
}
// Newborn j at census index idx[j] (n: absent): flight = end - tau, the newborn start rule (sink_core.hpp).
// code[j]: 0 ok, 1 no remaining flight, 100 + start-rule code; rank[j] = census index of a failing newborn.
extern "C" __global__ void k34_newborn_place(u64 count, u64 idx, u64 nb_tau, double end_s, hgp::Policy p, u64 soa, u64 n,
                                             u64 flight, u64 face, u64 newborn, u64 axis, u64 tau, u64 code, u64 rank, u64 nb_face) {
    const u64 j = IDX();
    if (j >= count) return;
    const u64 i = PCU64(idx)[j];
    PU32(code)[j] = 0; PU32(rank)[j] = S34_NONE; PI(nb_face)[2 * j] = -1; PI(nb_face)[2 * j + 1] = -1;
    if (i >= n) return;
    const double t = PCD(nb_tau)[j], F = end_s - t;
    PI(newborn)[i] = 1; PD(flight)[i] = F; PD(tau)[i] = t;
    if (!(F > 0)) { PU32(code)[j] = 1; PU32(rank)[j] = u32(i); return; }
    const double* s = PCD(soa);
    const hgp::V3 x = {{s[i], s[n + i], s[2 * n + i]}}, u = {{s[3 * n + i], s[4 * n + i], s[5 * n + i]}};
    int body = -1, ax = -1;
    const int c = hgp::newborn_start_rule(p, x, u, body, ax);
    if (c != hgp::ok) { PU32(code)[j] = 100u + u32(c); PU32(rank)[j] = u32(i); return; }
    PI(face)[i] = body; PI(axis)[i] = ax;
    PI(nb_face)[2 * j] = body; PI(nb_face)[2 * j + 1] = ax;
}
// K0 = KE(m, u_start) per carrier (pf7::kinetic_energy, the work scales of motion-work-gross-v1).
extern "C" __global__ void k34_k0(u64 n, u64 soa, u64 k0) {
    const u64 i = IDX();
    if (i >= n) return;
    const double* s = PCD(soa);
    PD(k0)[i] = pf7::kinetic_energy(s[7 * n + i], s[3 * n + i], s[4 * n + i], s[5 * n + i]);
}
// First resident neutral carrier (native v7 resident-neutral interval bound).
extern "C" __global__ void k34_neutral_fail(u64 n, u64 soa, u64 newborn, u64 fail) {
    const u64 i = IDX();
    if (i >= n) return;
    PU32(fail)[i] = (PCD(soa)[6 * n + i] == 0 && PCI(newborn)[i] == 0) ? u32(i) : S34_NONE;
}
// charge_of_particles (census_charge.hpp) of 10-word carrier records (GapParticle: id, x[3], u[3], q, m, w):
// the census residents (mask: newborn flags, masked carriers give sentinel records) or the endpoints (no mask).
// The sealed check order per carrier: operand, kinetic interior, unresolved words, deposited words, finite
// position, inside the mesh (push_core.hpp codes 40-44 and the interior codes).
extern "C" __global__ void k34_records_charge(u64 n, u64 rec, hgp::Mesh g, hgp::Policy p, u64 mask, u32 sentinel, u64 key, u64 pay,
                                             u64 hi, u64 lo, u64 code, u64 fail) {
    const u64 i = IDX();
    if (i >= n) return;
    const u64* r = PCU64(rec) + 10 * i;
    const double w = pf7::double_of(r[9]), q = pf7::double_of(r[7]);
    const hgp::V3 x = {{pf7::double_of(r[1]), pf7::double_of(r[2]), pf7::double_of(r[3])}};
    int c = hgp::ok;
    hgp::Words qw = {0, 0};
    hgp::NodeWeights cw;
    const bool counted = !mask || PCI(mask)[i] == 0;
    if (counted) {
        if (!(hgp::finite(w) && w > 0 && hgp::finite(q))) c = hgp::canonical_charge_operand;
        else c = hgp::require_kinetic_interior(p, x);
        if (c == hgp::ok) {
            hgp::V3 xc = x;
            xc.x[0] -= p.width / 2;
            xc.x[1] -= p.gap / 2;
            const hgp::Words ww = {w, 0};
            qw = hgp::words_mul(ww, q);
            if (!(hgp::words_canonical(qw) && (q == 0 || qw.hi != 0 || qw.lo != 0))) c = hgp::unresolved_canonical_charge;
            else if (!hgp::words_canonical(qw)) c = hgp::deposited_words;
            else c = hgp::point_cell(g, xc, cw, hgp::canonical_charge_outside_mesh);
        }
    }
    for (int k = 0; k < 8; ++k) {
        const u64 rr = 8 * i + k;
        PU32(pay)[rr] = u32(rr);
        if (counted && c == hgp::ok) {
            PU32(key)[rr] = u32(hgp::slot_node(g, cw, k));
            const hgp::Words v = hgp::slot_words(cw, k, qw);
            PD(hi)[rr] = v.hi; PD(lo)[rr] = v.lo;
        } else { PU32(key)[rr] = sentinel; PD(hi)[rr] = 0; PD(lo)[rr] = 0; }
    }
    PU32(code)[i] = u32(c);
    PU32(fail)[i] = c != hgp::ok ? u32(i) : S34_NONE;
}
// Constant fill of u32 words.
extern "C" __global__ void k34_fill_u32(u64 n, u64 p, u32 v) {
    const u64 i = IDX();
    if (i < n) PU32(p)[i] = v;
}
// Rows of a device census (10-word records) selected by census index (O(sinks + newborns) downloads).
extern "C" __global__ void k34_gather_rows(u64 count, u64 rows, u64 idx, u64 out) {
    const u64 j = IDX();
    if (j >= count) return;
    const u64 i = PCU64(idx)[j];
    for (int w = 0; w < 10; ++w) PU64(out)[10 * j + w] = PCU64(rows)[10 * i + w];
}
// First differing index of two u64 arrays (bitwise comparison).
extern "C" __global__ void k34_differ_u64(u64 n, u64 a, u64 b, u64 fail) {
    const u64 i = IDX();
    if (i >= n) return;
    PU32(fail)[i] = PCU64(a)[i] == PCU64(b)[i] ? S34_NONE : u32(i);
}
// Net words (sequential Words sums in carrier order) and x87 gross of each node's records.
extern "C" __global__ void k34_charge_nodes(u32 nn, u64 pay, u64 begin, u64 end, u64 hi, u64 lo, u64 out_hi, u64 out_lo, u64 out_gross) {
    const u64 k = IDX();
    if (k >= nn) return;
    hgp::Words net = {0, 0};
    double gross = 0;
    for (u32 r = PCU32(begin)[k]; r < PCU32(end)[k]; ++r) {
        const u32 p = PCU32(pay)[r];
        const hgp::Words v = {PCD(hi)[p], PCD(lo)[p]};
        net = hgp::words_add(net, v);
        gross = hgp::gross_step(gross, v);
    }
    PD(out_hi)[k] = net.hi; PD(out_lo)[k] = net.lo; PD(out_gross)[k] = gross;
}
// free_charge_guard of the census charge and same_net_charge(accepted free charge, census charge).
extern "C" __global__ void k34_census_nodes(u32 nn, u64 c_hi, u64 c_lo, u64 c_gross, u64 a_hi, u64 a_lo, u64 a_gross, u64 guard_fail, u64 same_fail) {
    const u64 k = IDX();
    if (k >= nn) return;
    const double ch = PCD(c_hi)[k], cl = PCD(c_lo)[k], cg = PCD(c_gross)[k];
    PU32(guard_fail)[k] = sup34::free_guard_node(ch, cl, cg) ? S34_NONE : u32(k);
    PU32(same_fail)[k] = sup34::same_net_node(PCD(a_hi)[k], PCD(a_lo)[k], PCD(a_gross)[k], ch, cl, cg) ? S34_NONE : u32(k);
}

// ---------------------------------------------------------------- fold
extern "C" __global__ void k34_fail_nonfinite(u64 n, u64 values, u64 fail) {
    const u64 i = IDX();
    if (i >= n) return;
    PU32(fail)[i] = pf7::finite(PCD(values)[i]) ? S34_NONE : u32(i);
}
// Exact sums. kind: 0 K (w d, every carrier), 1 G_K (w |d|, charged), 2 C_K (w K0 then w RN(K0 + d), charged;
// 2n terms), 3 Q (|q| w, charged), 4 W (E J), 5 G_W (|E| |J|), 6 sum of squares of a.
extern "C" __global__ void k34_sum(int kind, S34SumArgs s, u64 partial, u64 status) {
    const u64 terms = kind == 2 ? 2 * s.n : s.n;
    xs_body(terms, [&](u64 i) {
        pf7::Term t; t.zero = 1; t.bad = 0; t.d = 0; t.neg = 0; for (int k = 0; k < 5; ++k) t.c[k] = 0;
        if (kind == 4) return pf7::product_term(PCD(s.a)[i], PCD(s.b)[i]);
        if (kind == 5) return pf7::product_term(fabs(PCD(s.a)[i]), fabs(PCD(s.b)[i]));
        if (kind == 6) return pf7::product_term(PCD(s.a)[i], PCD(s.a)[i]);
        const u64 c = kind == 2 && i >= s.n ? i - s.n : i;
        const double w = PCD(s.w)[c], d = PCD(s.dke)[c];
        if (kind == 0) return pf7::product_term(w, d);
        const double q = PCD(s.q)[c];
        if (q == 0) return t;
        if (kind == 1) return pf7::product_term(w, fabs(d));
        if (kind == 3) return pf7::product_term(fabs(q), w);
        const double K0 = PCD(s.k0)[c];
        return pf7::product_term(w, i < s.n ? K0 : K0 + d);
    }, (pf7_i64*)partial, PU32(status));
}
extern "C" __global__ void k34_sum_finish(u32 parts, u64 partial, u64 status, u64 terms, double scale, u64 out) {
    xs_finish(parts, (const pf7_i64*)partial, PCU32(status), terms, scale, (XsOut*)out);
}
// Continuity of the motion fold: div in ascending edge order (plain binary64, the sealed gradient_transpose
// loop), |endpoint - initial - dt div| per node (NaN skipped by the max, as std::max).
extern "C" __global__ void k34_continuity(u32 nn, u64 ntp, u64 nte, u64 nts, u64 current, u64 elen, u64 endpoint, u64 initial, double dt, u64 out) {
    const u64 k = IDX();
    if (k >= nn) return;
    double div = 0.0;
    for (int j = PCI(ntp)[k]; j < PCI(ntp)[k + 1]; ++j) {
        const int e = PCI(nte)[j];
        const double q = PCD(current)[e] / PCD(elen)[e];
        div = PCI(nts)[j] < 0 ? div - q : div + q;
    }
    PD(out)[k] = fabs(PCD(endpoint)[k] - PCD(initial)[k] - dt * div);
}
// Nodes on an edge carrying J != 0 (motion-work-gross-v1 E^).
extern "C" __global__ void k34_node_marks(u32 nn, u64 ntp, u64 nte, u64 current, u64 mark) {
    const u64 k = IDX();
    if (k >= nn) return;
    int m = 0;
    for (int j = PCI(ntp)[k]; j < PCI(ntp)[k + 1]; ++j) if (PCD(current)[PCI(nte)[j]] != 0) m = 1;
    PI(mark)[k] = m;
}
extern "C" __global__ void k34_edge_bound(u32 ne, u64 ea, u64 eb, u64 mark, u64 E, u64 out) {
    const u64 e = IDX();
    if (e >= ne) return;
    PD(out)[e] = (PCI(mark)[PCI(ea)[e]] || PCI(mark)[PCI(eb)[e]]) ? fabs(PCD(E)[e]) : 0.0;
}
extern "C" __global__ void k34_free_guard(u32 nn, u64 hi, u64 lo, u64 gross, u64 fail) {
    const u64 k = IDX();
    if (k >= nn) return;
    PU32(fail)[k] = sup34::free_guard_node(PCD(hi)[k], PCD(lo)[k], PCD(gross)[k]) ? S34_NONE : u32(k);
}

// ---------------------------------------------------------------- record
// MotionAssembly::assemble operand checks: path start = census position + frame shift, displacement =
// (end + shift) - start, sink flight = F fraction and discarded = F - flight.
// outs: push carrier outputs (endpoint 6n, start 3n, disp 3n); sink_real: 7n (flight, discarded, fraction, plane, segment).
extern "C" __global__ void k34_assemble_check(u64 n, u64 soa, u64 outs, u64 sink_kind, u64 sink_real, u64 flight, double sx, double sy, u64 fail) {
    const u64 i = IDX();
    if (i >= n) return;
    const double* s = PCD(soa); const double* o = PCD(outs);
    const double shift[3] = {sx, sy, 0};
    bool ok = true;
    for (int d = 0; d < 3; ++d) {
        const double start = o[(6 + d) * n + i], disp = o[(9 + d) * n + i], end = o[d * n + i];
        if (!(start == s[d * n + i] + shift[d] && disp == (end + shift[d]) - start)) ok = false;
    }
    if (PCI(sink_kind)[i] >= 0) {
        const double F = PCD(flight)[i], f = PCD(sink_real)[i], disc = PCD(sink_real)[n + i], frac = PCD(sink_real)[2 * n + i];
        if (!(f == F * frac && disc == F - f)) ok = false;
    }
    PU32(fail)[i] = ok ? S34_NONE : u32(i);
}
extern "C" __global__ void k34_record_flags(u64 n, u64 sink_kind, u64 newborn, u64 surv, u64 sink, u64 born) {
    const u64 i = IDX();
    if (i >= n) return;
    const int k = PCI(sink_kind)[i];
    PU32(surv)[i] = k < 0 ? 1u : 0u; PU32(sink)[i] = k >= 0 ? 1u : 0u; PU32(born)[i] = PCI(newborn)[i] ? 1u : 0u;
}
// Canonical record arrays (API-007.txt 1.2), ascending native ID; sink_index/birth_index: census index.
extern "C" __global__ void k34_record_build(u64 n, u64 ids, u64 soa, u64 outs, u64 sink_kind, u64 sink_axis, u64 sink_real, u64 flight, u64 tau,
                                            u64 off_surv, u64 off_sink, u64 endpoints, u64 paths, u64 sinks, u64 sink_index, int main_records, u64 sink_extra) {
    const u64 i = IDX();
    if (i >= n) return;
    const double* s = PCD(soa); const double* o = PCD(outs); const double* sr = PCD(sink_real);
    const u64 id = PCU64(ids)[i];
    const int kind = PCI(sink_kind)[i];
    if (main_records) {
        u64* path = PU64(paths) + 8 * i;
        path[0] = id;
        for (int d = 0; d < 3; ++d) { path[1 + d] = pf7::bits_of(o[(6 + d) * n + i]); path[4 + d] = pf7::bits_of(o[(9 + d) * n + i]); }
        path[7] = pf7::bits_of(kind >= 0 ? sr[i] : PCD(flight)[i]);
    }
    if (kind < 0) {
        if (!main_records) return;
        u64* e = PU64(endpoints) + 10 * PCU32(off_surv)[i];
        e[0] = id;
        for (int d = 0; d < 6; ++d) e[1 + d] = pf7::bits_of(o[d * n + i]);
        for (int k = 0; k < 3; ++k) e[7 + k] = pf7::bits_of(s[(6 + k) * n + i]);
    } else {
        const u32 j = PCU32(off_sink)[i];
        u64* r = PU64(sinks) + 17 * j;
        const bool collected = kind == 0 || kind == 1 || kind == 3;
        r[0] = id; r[1] = collected ? 1 : 2; r[2] = u64(kind); r[3] = u64((long long)PCI(sink_axis)[i]);
        r[4] = pf7::bits_of(sr[3 * n + i]); r[5] = pf7::bits_of(sr[2 * n + i]);
        for (int d = 0; d < 3; ++d) { r[6 + d] = pf7::bits_of(sr[(4 + d) * n + i]); r[9 + d] = pf7::bits_of(o[d * n + i]); r[12 + d] = pf7::bits_of(o[(3 + d) * n + i]); }
        r[15] = pf7::bits_of(sr[i]); r[16] = pf7::bits_of(sr[n + i]);
        PU64(sink_index)[j] = i;
        PD(sink_extra)[2 * j] = o[13 * n + i]; PD(sink_extra)[2 * j + 1] = o[14 * n + i];   // gyro, displacement cells
    }
}
extern "C" __global__ void k34_birth_build(u64 n, u64 ids, u64 newborn, u64 flight, u64 tau, u64 off_born, u64 births, u64 birth_index) {
    const u64 i = IDX();
    if (i >= n || PCI(newborn)[i] == 0) return;
    const u32 j = PCU32(off_born)[i];
    u64* b = PU64(births) + 3 * j;
    b[0] = PCU64(ids)[i]; b[1] = pf7::bits_of(PCD(tau)[i]); b[2] = pf7::bits_of(PCD(flight)[i]);
    PU64(birth_index)[j] = i;
}

// ---------------------------------------------------------------- deposit support
// Per carrier: x87 terms t1, t2, the node box and its node count (0 without a contribution).
extern "C" __global__ void k34_support_carriers(u64 n, sup34::Geo g, u64 soa, u64 outs, u64 t1, u64 t2, u64 box, u64 count, u64 fail, u64 code) {
    const u64 i = IDX();
    if (i >= n) return;
    const double* s = PCD(soa); const double* o = PCD(outs);
    const double start[3] = {s[i], s[n + i], s[2 * n + i]}, end[3] = {o[i], o[n + i], o[2 * n + i]};
    sup34::Carrier c;
    sup34::carrier(g, s[8 * n + i], s[6 * n + i], start, end, c);
    ((x87::XBits*)t1)[i] = x87::pack(c.t1);
    ((x87::XBits*)t2)[i] = x87::pack(c.t2);
    for (int d = 0; d < 3; ++d) { PI(box)[6 * i + d] = c.lo[d]; PI(box)[6 * i + 3 + d] = c.hi[d]; }
    PU32(count)[i] = sup34::box_nodes(c);
    PU32(code)[i] = u32(c.code);
    PU32(fail)[i] = c.code != sup34::ok ? u32(i) : S34_NONE;
}
// (node, carrier) records in carrier order, nodes of a carrier in the sealed box loop order.
extern "C" __global__ void k34_support_emit(u64 n, int nx, int ny, u64 box, u64 count, u64 offset, u64 key, u64 pay) {
    const u64 i = IDX();
    if (i >= n) return;
    const u32 m = PCU32(count)[i];
    if (!m) return;
    sup34::Carrier c;
    c.contributes = 1;
    for (int d = 0; d < 3; ++d) { c.lo[d] = PCI(box)[6 * i + d]; c.hi[d] = PCI(box)[6 * i + 3 + d]; }
    const u32 base = PCU32(offset)[i];
    for (u32 k = 0; k < m; ++k) { PU32(key)[base + k] = u32(sup34::box_node(c, k, nx, ny)); PU32(pay)[base + k] = u32(i); }
}
// Per node: s1, s2, m in ascending carrier order (stable sort), then the sealed finish.
extern "C" __global__ void k34_support_nodes(u32 nn, sup34::Geo g, u64 pay, u64 begin, u64 end, u64 t1, u64 t2, u64 out) {
    const u64 k = IDX();
    if (k >= nn) return;
    x87::X s1 = x87::zero(), s2 = x87::zero();
    x87_u64 m = 0;
    for (u32 r = PCU32(begin)[k]; r < PCU32(end)[k]; ++r) {
        const u32 c = PCU32(pay)[r];
        s1 = x87::add(s1, x87::unpack(((const x87::XBits*)t1)[c]));
        s2 = x87::add(s2, x87::unpack(((const x87::XBits*)t2)[c]));
        ++m;
    }
    PD(out)[k] = sup34::finish(g, s1, s2, m);
}

// ---------------------------------------------------------------- block lists of binary64 arrays (native v7.1 tags)
// Same construction as the forms library block lists (API-007.txt 1.2) with 8-byte records and the tags
// native-v7.1/field-E, native-v7.1/field-B, native-v7.1/current-J (API-071.txt 1.1).
__device__ const char* f64_tag(u32 which, int level, u32& len) {
    const char* const t[9] = {"native-v7.1/field-E/chunk", "native-v7.1/field-E/block", "native-v7.1/field-E/root",
                              "native-v7.1/field-B/chunk", "native-v7.1/field-B/block", "native-v7.1/field-B/root",
                              "native-v7.1/current-J/chunk", "native-v7.1/current-J/block", "native-v7.1/current-J/root"};
    const char* s = t[which * 3 + level];
    len = 0; while (s[len]) ++len;
    return s;
}
extern "C" __global__ void k34_f64_chunks(u32 which, u64 n, u64 words, u64 out) {
    const u64 j = IDX();
    const u64 chunks = (n + pf7::chunk_records - 1) / pf7::chunk_records;
    if (j >= chunks) return;
    const u64 m = n - j * pf7::chunk_records < pf7::chunk_records ? n - j * pf7::chunk_records : pf7::chunk_records;
    u32 len; const char* tag = f64_tag(which, 0, len);
    pf7::Sha sh; pf7::sha_init(sh); pf7::sha_text(sh, tag, len); pf7::sha_u64(sh, j); pf7::sha_u64(sh, m);
    const u64* w = PCU64(words) + j * pf7::chunk_records;
    for (u64 r = 0; r < m; ++r) pf7::sha_u64(sh, w[r]);
    pf7::sha_finish(sh, PU32(out) + 8 * j);
}
extern "C" __global__ void k34_f64_blocks(u32 which, u64 n, u64 chunks_in, u64 out) {
    const u64 k = IDX();
    const u64 leaves = (n + pf7::leaf_records - 1) / pf7::leaf_records;
    if (k >= leaves) return;
    const u64 records = n - k * pf7::leaf_records < pf7::leaf_records ? n - k * pf7::leaf_records : pf7::leaf_records;
    const u64 q = (records + pf7::chunk_records - 1) / pf7::chunk_records;
    u32 len; const char* tag = f64_tag(which, 1, len);
    pf7::Sha sh; pf7::sha_init(sh); pf7::sha_text(sh, tag, len); pf7::sha_u64(sh, k); pf7::sha_u64(sh, records); pf7::sha_u64(sh, q);
    for (u64 c = 0; c < q; ++c) pf7::sha_digest_words(sh, PCU32(chunks_in) + 8 * (k * pf7::leaf_chunks + c));
    pf7::sha_finish(sh, PU32(out) + 8 * k);
}
// One block: native v5 tree over the leaves (pairs left to right, an odd last node carried up), then the root.
extern "C" __global__ void k34_f64_tree(u32 which, u64 n, u64 leaves, u64 leaf, u64 a, u64 b, u64 root_out) {
    if (threadIdx.x != 0) return;
    u32 tree[8];
    if (leaves == 0) pf7::merkle_empty(tree);
    else {
        const u32* src = PCU32(leaf);
        u32* bufs[2] = {PU32(a), PU32(b)};
        u64 count = leaves; int which_buf = 0;
        while (count > 1) {
            u32* dst = bufs[which_buf];
            u64 m = 0;
            for (u64 p = 0; p + 1 < count; p += 2) pf7::merkle_node(src + 8 * p, src + 8 * (p + 1), dst + 8 * m++);
            if (count & 1) { for (int w = 0; w < 8; ++w) dst[8 * m + w] = src[8 * (count - 1) + w]; ++m; }
            src = dst; count = m; which_buf ^= 1;
        }
        for (int w = 0; w < 8; ++w) tree[w] = src[w];
    }
    u32 len; const char* tag = f64_tag(which, 2, len);
    pf7::Sha sh; pf7::sha_init(sh); pf7::sha_text(sh, tag, len); pf7::sha_u64(sh, pf7::leaf_records); pf7::sha_u64(sh, n); pf7::sha_u64(sh, leaves);
    pf7::sha_digest_words(sh, tree);
    pf7::sha_finish(sh, PU32(root_out));
}
