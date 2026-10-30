// push_bodies.hpp -- per-thread bodies of the push kernels. Included by the NVRTC
// program (push_kernels.cu, one thread per index) and by host code that runs the
// identical logic sequentially (port checks against the sealed CPU reference).
#pragma once
#include "push_abi.hpp"

namespace hgp {

HGP_BODY void carrier_body(const CarrierArgs& a, hgp_u32 i) {
    if (i >= a.n) return;
    const Mesh& g = a.g;
    const hgp_u32 n = a.n;
    const V3 x = {{a.x[i], a.y[i], a.z[i]}};
    const V3 u = {{a.ux[i], a.uy[i], a.uz[i]}};
    const double q = a.q[i], m = a.m[i], w = a.w[i];

    // ---- first census loop: carrier validity, interior, plasma operand, start deposit cell
    int c1 = ok;
    if (!(finite(w) && w > 0 && finite(m) && m > 0 && finite(q))) c1 = invalid_carrier;
    else if (!finite3(x) || !finite3(u)) c1 = nonfinite_kinematics;
    else c1 = require_kinetic_interior(a.p, x);
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

    // ---- second loop: motion, first hit, endpoint, canonical path, endpoint deposit cell
    int c2 = c1 == ok ? ok : skipped;
    bool moved_valid = false, kinetic_valid = false;
    V3 movedD = {{0, 0, 0}}, movedU = {{0, 0, 0}};
    Hit mhit = {1, -1};
    double flight = 0, gyro = 0, dcells = 0;
    KineticHit kh = {sink_none, 1, -1, 0};
    if (c2 == ok) {
        if (q == 0) {
            V3 se;
            KineticHit sh;
            c2 = neutral_stream(a.p, x, u, q, m, w, a.dt, se, sh);
            if (c2 == ok) {
                if (sh.sink != sink_none) { c2 = neutral_first_hit; kh = sh; kinetic_valid = true; }
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
    if (moved_valid) {
        KineticHit k2;
        khc = first_kinetic_hit(a.p, x, movedD, k2);
        if (khc == ok) { kh = k2; kinetic_valid = true; }
    }
    if (c2 == ok) {
        if (!(mhit.surface < 0 && flight == a.dt)) c2 = boundary_owner;
        else if (khc != ok) c2 = khc;
        else if (kh.sink != sink_none) c2 = charged_first_hit;
    }
    V3 e = {{0, 0, 0}}, Dp = {{0, 0, 0}};
    if (c2 == ok) {
        for (int d = 0; d < 3; ++d) e.x[d] = x.x[d] + movedD.x[d];
        c2 = require_kinetic_interior(a.p, e);
    }
    if (c2 == ok) {
        for (int d = 0; d < 3; ++d) Dp.x[d] = (e.x[d] + a.off[d]) - fx.x[d];
        if (!(intersection(g, fx, Dp).surface < 0)) c2 = rounded_endpoint_crosses;
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

    // ---- corrected free-charge census of the endpoint (charge_of_endpoints)
    int c3 = c2 == ok ? ok : skipped;
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
    a.fail3[i] = (c2 == ok && c3 != ok) ? i : 0xffffffffu;
    a.endpoint[i] = e.x[0]; a.endpoint[n + i] = e.x[1]; a.endpoint[2 * n + i] = e.x[2];
    a.endpoint[3 * n + i] = good ? movedU.x[0] : 0; a.endpoint[4 * n + i] = good ? movedU.x[1] : 0; a.endpoint[5 * n + i] = good ? movedU.x[2] : 0;
    for (int d = 0; d < 3; ++d) { a.start[d * n + i] = fx.x[d]; a.disp[d * n + i] = Dp.x[d]; }
    a.dke[i] = good ? kinetic_energy(m, movedU) - kinetic_energy(m, u) : 0;
    a.gyro[i] = good ? gyro : 0;
    a.dcells[i] = good ? dcells : 0;
    a.edge_count[i] = good ? hgp_u32(24 * segments) : 0u;
    if (a.diagnostics) {
        double* r = a.diag_real;
        int* t = a.diag_int;
        for (int d = 0; d < 3; ++d) { r[(diag_moved_dx + d) * n + i] = movedD.x[d]; r[(diag_moved_ux + d) * n + i] = movedU.x[d]; }
        r[diag_hit_fraction * n + i] = mhit.fraction;
        r[diag_flight_s * n + i] = flight;
        r[diag_gyro * n + i] = gyro;
        r[diag_displacement_cells * n + i] = dcells;
        r[diag_kinetic_fraction * n + i] = kinetic_valid ? kh.fraction : 0;
        r[diag_kinetic_plane * n + i] = kinetic_valid ? kh.plane : 0;
        t[diag_moved_valid * n + i] = moved_valid ? 1 : 0;
        t[diag_hit_surface * n + i] = moved_valid ? mhit.surface : -1;
        t[diag_kinetic_valid * n + i] = kinetic_valid ? 1 : 0;
        t[diag_kinetic_sink * n + i] = kinetic_valid ? kh.sink : -1;
        t[diag_kinetic_axis * n + i] = kinetic_valid ? kh.axis : -1;
    }
}

// Edge records of the canonical path deposit (axial::Path::deposit), in list order.
HGP_BODY void emit_body(const EmitArgs& a, hgp_u32 i) {
    if (i >= a.n) return;
    const hgp_u32 count = a.edge_count[i];
    if (count == 0) return;
    const Mesh& g = a.g;
    const hgp_u32 n = a.n;
    const V3 x = {{a.start[i], a.start[n + i], a.start[2 * n + i]}};
    const V3 D = {{a.disp[i], a.disp[n + i], a.disp[2 * n + i]}};
    const double charge = a.w[i] * a.q[i];
    Cuts c;
    if (path_cuts(g, x, D, c) != ok) return;   // validated by hgp_carrier
    hgp_u32 r = a.edge_offset[i];
    for (int segment = 1; segment < c.count; ++segment) {
        const double midpoint = .5 * (c.t[segment] + c.t[segment - 1]), span = c.t[segment] - c.t[segment - 1];
        for (int qi = 0; qi < 2; ++qi) {
            const double gq = qi ? gauss_node : -gauss_node;
            const double t = midpoint + .5 * span * gq, w = .5 * span;
            int cell[3];
            double f[3];
            if (path_point(g, x, D, t, cell, f) != ok) return;
            const int ci = cell[0], cj = cell[1], ck = cell[2];
            for (int rr = 0; rr < 2; ++rr)
                for (int s = 0; s < 2; ++s) {
                    const double v0 = w * (rr ? f[1] : 1 - f[1]) * (s ? f[2] : 1 - f[2]);
                    const double v1 = w * (rr ? f[0] : 1 - f[0]) * (s ? f[2] : 1 - f[2]);
                    const double v2 = w * (rr ? f[0] : 1 - f[0]) * (s ? f[1] : 1 - f[1]);
                    a.key[r] = hgp_u32(ex(g, ci, cj + rr, ck + s)); a.pay[r] = r; a.value[r] = charge * D.x[0] * v0 / a.dt; ++r;
                    a.key[r] = hgp_u32(ey(g, ci + rr, cj, ck + s)); a.pay[r] = r; a.value[r] = charge * D.x[1] * v1 / a.dt; ++r;
                    a.key[r] = hgp_u32(ez(g, ci + rr, cj + s, ck)); a.pay[r] = r; a.value[r] = charge * D.x[2] * v2 / a.dt; ++r;
                }
        }
    }
}

// current[e] = sum of the edge's records in (carrier ID, list position) order.
HGP_BODY void reduce_edges_body(const NodeReduceArgs& a, hgp_u32 e) {
    if (e >= a.nn) return;
    double acc = 0;
    for (hgp_u32 r = a.begin[e]; r < a.end[e]; ++r) acc += a.v0[a.pay[r]];
    a.out0[e] = acc;
}

// Start-position deposits: omega (plasma load), initial charge, gross |charge|.
HGP_BODY void reduce_start_body(const NodeReduceArgs& a, hgp_u32 k) {
    if (k >= a.nn) return;
    double omega = 0, init = 0, gross = 0;
    for (hgp_u32 r = a.begin[k]; r < a.end[k]; ++r) {
        const hgp_u32 p = a.pay[r];
        omega += a.v0[p];
        init += a.v1[p];
        gross += a.v2[p];
    }
    a.out0[k] = omega; a.out1[k] = init; a.out2[k] = gross;
}

// Endpoint deposits: endpoint charge (binary64) and the corrected free charge
// (compensated net words, x87-rounded gross) of charge_of_endpoints.
HGP_BODY void reduce_end_body(const EndNodeArgs& a, hgp_u32 k) {
    if (k >= a.nn) return;
    double acc = 0, gross = 0;
    Words net = {0, 0};
    for (hgp_u32 r = a.begin[k]; r < a.end[k]; ++r) {
        const hgp_u32 p = a.pay[r];
        acc += a.charge[p];
        const Words v = {a.hi[p], a.lo[p]};
        net = words_add(net, v);
        gross = gross_step(gross, v);
    }
    a.out_charge[k] = acc; a.out_hi[k] = net.hi; a.out_lo[k] = net.lo; a.out_gross[k] = gross;
}

}  // namespace hgp
