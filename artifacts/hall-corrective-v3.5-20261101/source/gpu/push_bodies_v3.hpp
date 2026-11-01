// push_bodies_v3.hpp -- per-thread bodies of the Hall v3 motion kernels. The carrier
// body is the v1 body (hall-gpu-resident-v1-20261030/push/push_bodies.hpp) with the
// v3 first-hit truncation: a charged carrier whose represented path hits a field
// surface or the kinetic boundary is truncated by sink_truncate (sink_core.hpp, the
// same source as the CPU provider). Sink carriers keep their endpoint charge record
// (the bank deposit point) but leave the free-charge census. Compiled by NVRTC and
// by the host (port checks), like the v1 bodies.
#pragma once
#include "push_abi_v3.hpp"

namespace hgp {

HGP_BODY void carrier_body_v3(const CarrierArgsV3& args, hgp_u32 i) {
    const CarrierArgs& a = args.v1;
    if (i >= a.n) return;
    const Mesh& g = a.g;
    const hgp_u32 n = a.n;
    const V3 x = {{a.x[i], a.y[i], a.z[i]}};
    const V3 u = {{a.ux[i], a.uy[i], a.uz[i]}};
    const double q = a.q[i], m = a.m[i], w = a.w[i];
    const double F = args.flight ? args.flight[i] : a.dt;           // v3: newborn flight end - tau
    const int face = args.face ? args.face[i] : -1;                 // v3: newborn foil-face start
    const bool newborn = args.newborn && args.newborn[i] != 0;

    // ---- first census loop (v1; v3: a newborn face start skips the interior rule)
    int c1 = ok;
    if (!(finite(w) && w > 0 && finite(m) && m > 0 && finite(q))) c1 = invalid_carrier;
    else if (!finite3(x) || !finite3(u)) c1 = nonfinite_kinematics;
    else if (face < 0) c1 = require_kinetic_interior(a.p, x);
    V3 fx = {{0, 0, 0}};
    NodeWeights sw;
    double load = 0;
    if (c1 == ok) {
        for (int d = 0; d < 3; ++d) fx.x[d] = x.x[d] + a.off[d];
        load = w * q * q / (eps0 * m);
        if (!(finite(load) && (q == 0 || load > 0))) c1 = unresolved_plasma_operand;
        else c1 = point_cell(g, fx, sw, charge_outside_mesh);
    }
    const double charge = w * q;
    for (int s = 0; s < 8; ++s) {
        const hgp_u32 r = 8 * i + s;
        a.start_pay[r] = r;
        if (c1 == ok) {
            a.start_key[r] = hgp_u32(slot_node(g, sw, s));
            a.start_omega[r] = slot_value(sw, s, load);
            a.start_init[r] = slot_value(sw, s, charge);
            a.start_gross[r] = slot_value(sw, s, HGP_FABS(charge));
        } else {
            a.start_key[r] = a.node_sentinel;
            a.start_omega[r] = 0; a.start_init[r] = 0; a.start_gross[r] = 0;
        }
    }

    // ---- second loop: motion, first hit (v3: truncation), endpoint, canonical path
    int c2 = c1 == ok ? ok : skipped;
    bool moved_valid = false, sink = false;
    V3 movedD = {{0, 0, 0}}, movedU = {{0, 0, 0}};
    Hit mhit = {1, -1};
    double flight = 0, gyro = 0, dcells = 0;
    KineticHit kh = {sink_none, 1, -1, 0};
    SinkOut so;
    so.sink = sink_none; so.axis = -1; so.plane = 0; so.fraction = 0; so.flight = 0; so.discarded = 0;
    so.segment.x[0] = so.segment.x[1] = so.segment.x[2] = 0;
    if (c2 == ok && newborn) {
        // v3 newborn over [tau, end]: same code as MotionAssembly::cpu (ports, face-aware first hit)
        if (q == 0) {
            V3 end; double dc = 0;
            c2 = neutral_drift(g, args.sg, x, u, F, face, so, end, dc);
            if (c2 == ok) {
                if (so.sink != sink_none) { sink = true; movedU = u; flight = so.flight; gyro = 0; dcells = so.dcells; }
                else {
                    movedU = u; flight = F; for (int d = 0; d < 3; ++d) movedD.x[d] = end.x[d] - x.x[d]; dcells = dc; moved_valid = true;
                    if (!(dcells <= .2)) c2 = neutral_displacement;
                }
            }
        } else {
            PushResult pr;
            c2 = push(g, fx, u, q, m, a.E, a.B, F, pr);
            if (c2 == ok) {
                KineticHit k2;
                c2 = first_hit_from(a.p, x, pr.D, face, k2);
                if (c2 == ok && (!(pr.hit.surface < 0 && pr.flight_s == F) || k2.sink != sink_none)) {
                    c2 = sink_truncate(g, args.sg, x, fx, u, q, m, a.E, a.B, F, so, face);
                    sink = c2 == ok;
                    if (sink) { movedU = so.u; flight = so.flight; gyro = so.gyro; dcells = so.dcells; }
                } else if (c2 == ok) {
                    movedD = pr.D; movedU = pr.u; mhit = pr.hit; flight = pr.flight_s; gyro = pr.gyro; dcells = pr.displacement_cells; moved_valid = true;
                }
            }
        }
    } else if (c2 == ok) {
        if (q == 0) {
            V3 se;
            KineticHit sh;
            c2 = neutral_stream(a.p, x, u, q, m, w, a.dt, se, sh);
            if (c2 == ok) {
                if (sh.sink != sink_none) {   // v3: neutral first hit -> sink (v2 refused)
                    V3 end; double dc = 0;
                    c2 = neutral_drift(g, args.sg, x, u, a.dt, -1, so, end, dc);
                    if (c2 == ok && so.sink == sink_none) c2 = lost_kinetic_surface;
                    sink = c2 == ok;
                    if (sink) { movedU = u; flight = so.flight; gyro = 0; dcells = so.dcells; }
                }
                else {
                    movedU = u;
                    flight = a.dt;
                    for (int d = 0; d < 3; ++d) movedD.x[d] = se.x[d] - x.x[d];
                    dcells = HGP_SQRT(hall_dot(movedD, movedD)) / min_spacing(g);
                    moved_valid = true;
                    if (!(dcells <= .2)) c2 = neutral_displacement;
                }
            }
        } else {
            PushResult pr;
            c2 = push(g, fx, u, q, m, a.E, a.B, a.dt, pr);
            if (c2 == ok) {
                movedD = pr.D; movedU = pr.u; mhit = pr.hit; flight = pr.flight_s; gyro = pr.gyro; dcells = pr.displacement_cells;
                moved_valid = true;
            }
        }
    }
    int khc = -1;
    if (moved_valid && !newborn) {
        KineticHit k2;
        khc = first_kinetic_hit(a.p, x, movedD, k2);
        if (khc == ok) kh = k2;
    }
    if (c2 == ok && !sink && !newborn && q != 0 && (!(mhit.surface < 0 && flight == a.dt) || (khc == ok && kh.sink != sink_none))) {
        c2 = sink_truncate(g, args.sg, x, fx, u, q, m, a.E, a.B, a.dt, so);
        sink = c2 == ok;
        if (sink) { movedU = so.u; flight = so.flight; gyro = so.gyro; dcells = so.dcells; }
    }
    if (c2 == ok && !sink && !newborn) {
        if (!(mhit.surface < 0 && flight == a.dt)) c2 = boundary_owner;
        else if (khc != ok) c2 = khc;
        else if (kh.sink != sink_none) c2 = charged_first_hit;
    }
    if (c2 == ok && !sink && newborn && !(mhit.surface < 0 && flight == F)) c2 = boundary_owner;
    V3 e = {{0, 0, 0}}, Dp = {{0, 0, 0}};
    if (c2 == ok) {
        if (sink) e = so.point;
        else {
            for (int d = 0; d < 3; ++d) e.x[d] = x.x[d] + movedD.x[d];
            c2 = require_kinetic_interior(a.p, e);
        }
    }
    if (c2 == ok) {
        for (int d = 0; d < 3; ++d) Dp.x[d] = (e.x[d] + a.off[d]) - fx.x[d];
        if (!sink && !(intersection(g, fx, Dp).surface < 0)) c2 = rounded_endpoint_crosses;
    }
    int segments = 0;
    if (c2 == ok) c2 = path_check(g, fx, Dp, segments);
    if (c2 == ok && !(finite(charge) && (q == 0 || charge != 0))) c2 = unresolved_macrocharge;
    NodeWeights ew;
    if (c2 == ok) {
        V3 fe;
        for (int d = 0; d < 3; ++d) fe.x[d] = e.x[d] + a.off[d];
        c2 = point_cell(g, fe, ew, charge_outside_mesh);
    }

    // ---- corrected free-charge census of the surviving endpoints (charge_of_endpoints)
    int c3 = (c2 == ok && !sink) ? ok : skipped;
    Words qw = {0, 0};
    NodeWeights cw;
    if (c3 == ok) {
        if (!(finite(w) && w > 0 && finite(q))) c3 = canonical_charge_operand;
        else c3 = require_kinetic_interior(a.p, e);
    }
    if (c3 == ok) {
        V3 xc = e;
        xc.x[0] -= a.p.width / 2;
        xc.x[1] -= a.p.gap / 2;
        const Words ww = {w, 0};
        qw = words_mul(ww, q);
        if (!(words_canonical(qw) && (q == 0 || qw.hi != 0 || qw.lo != 0))) c3 = unresolved_canonical_charge;
        else if (!words_canonical(qw)) c3 = deposited_words;
        else c3 = point_cell(g, xc, cw, canonical_charge_outside_mesh);
        if (c3 == ok)
            for (int d = 0; d < 3; ++d)
                if (cw.cell[d] != ew.cell[d] || hgp_bits(cw.a[d]) != hgp_bits(ew.a[d])) c3 = gpu_frame_mismatch;
    }

    // ---- outputs
    const bool good = c2 == ok;
    for (int s = 0; s < 8; ++s) {
        const hgp_u32 r = 8 * i + s;
        a.end_pay[r] = r;
        if (good) {
            a.end_key[r] = hgp_u32(slot_node(g, ew, s));
            a.end_charge[r] = slot_value(ew, s, charge);
            if (c3 == ok) { const Words v = slot_words(cw, s, qw); a.end_hi[r] = v.hi; a.end_lo[r] = v.lo; }
            else { a.end_hi[r] = 0; a.end_lo[r] = 0; }
        } else {
            a.end_key[r] = a.node_sentinel;
            a.end_charge[r] = 0; a.end_hi[r] = 0; a.end_lo[r] = 0;
        }
    }
    a.code1[i] = c1; a.code2[i] = c2; a.code3[i] = c3;
    a.fail1[i] = c1 != ok ? i : 0xffffffffu;
    a.fail2[i] = (c1 == ok && c2 != ok) ? i : 0xffffffffu;
    a.fail3[i] = (c2 == ok && !sink && c3 != ok) ? i : 0xffffffffu;
    a.endpoint[i] = e.x[0]; a.endpoint[n + i] = e.x[1]; a.endpoint[2 * n + i] = e.x[2];
    a.endpoint[3 * n + i] = good ? movedU.x[0] : 0; a.endpoint[4 * n + i] = good ? movedU.x[1] : 0; a.endpoint[5 * n + i] = good ? movedU.x[2] : 0;
    for (int d = 0; d < 3; ++d) { a.start[d * n + i] = fx.x[d]; a.disp[d * n + i] = Dp.x[d]; }
    a.dke[i] = good ? kinetic_energy(m, movedU) - kinetic_energy(m, u) : 0;
    a.gyro[i] = good ? gyro : 0;
    a.dcells[i] = good ? dcells : 0;
    a.edge_count[i] = good ? hgp_u32(24 * segments) : 0u;
    args.sink_kind[i] = (good && sink) ? so.sink : -1;
    args.sink_axis[i] = (good && sink) ? so.axis : -1;
    args.sink_real[i] = (good && sink) ? so.flight : 0;
    args.sink_real[n + i] = (good && sink) ? so.discarded : 0;
    args.sink_real[2 * n + i] = (good && sink) ? so.fraction : 0;
    args.sink_real[3 * n + i] = (good && sink) ? so.plane : 0;
    for (int d = 0; d < 3; ++d) args.sink_real[(4 + d) * n + i] = (good && sink) ? so.segment.x[d] : 0;
    if (a.diagnostics) {
        double* r = a.diag_real;
        int* t = a.diag_int;
        for (int d = 0; d < 3; ++d) { r[(diag_moved_dx + d) * n + i] = movedD.x[d]; r[(diag_moved_ux + d) * n + i] = movedU.x[d]; }
        r[diag_hit_fraction * n + i] = mhit.fraction;
        r[diag_flight_s * n + i] = flight;
        r[diag_gyro * n + i] = gyro;
        r[diag_displacement_cells * n + i] = dcells;
        r[diag_kinetic_fraction * n + i] = kh.fraction;
        r[diag_kinetic_plane * n + i] = kh.plane;
        t[diag_moved_valid * n + i] = moved_valid ? 1 : 0;
        t[diag_hit_surface * n + i] = moved_valid ? mhit.surface : -1;
        t[diag_kinetic_valid * n + i] = khc == ok ? 1 : 0;
        t[diag_kinetic_sink * n + i] = kh.sink;
        t[diag_kinetic_axis * n + i] = kh.axis;
    }
}

// Endpoint node reduction: endpoint charge of every record (incl. sink carriers);
// free-charge words and x87 gross of the surviving carriers only, in record order.
HGP_BODY void reduce_end_body_v3(const EndNodeArgsV3& a, hgp_u32 k) {
    if (k >= a.nn) return;
    double acc = 0, gross = 0;
    Words net = {0, 0};
    for (hgp_u32 r = a.begin[k]; r < a.end[k]; ++r) {
        const hgp_u32 p = a.pay[r];
        acc += a.charge[p];
        if (a.sink_kind[p >> 3] >= 0) continue;
        const Words v = {a.hi[p], a.lo[p]};
        net = words_add(net, v);
        gross = gross_step(gross, v);
    }
    a.out_charge[k] = acc; a.out_hi[k] = net.hi; a.out_lo[k] = net.lo; a.out_gross[k] = gross;
}

}  // namespace hgp
