// sink_core.hpp -- first-hit truncation of one charged carrier (Hall corrective v3).
// Compiled twice from the same bytes: by the host compiler (-ffp-contract=off) for
// the CPU motion provider and by NVRTC (--fmad=false, IEEE division and square
// root) for the GPU motion provider, so both produce identical truncated paths
// and sink records. Uses only the transcribed sealed arithmetic of push_core.hpp
// (hall-gpu-resident-v1-20261030/push, unchanged): trajectory (the sealed
// axial::push work-path fixed point), first_kinetic_hit (sealed kinetic_boundary.cpp).
//
// Kinetic model (native-finite-target-kinetic-policy-v1) and the native v6 sink
// record (pic-reaction-transactions-20261029/native-full-time-candidate006,
// API-006.txt section 2): the first hit is evaluated on the integrator's
// represented straight segment D of the full flight F = dt (the converged sealed
// work path, trajectory(dt).D); f = first_kinetic_hit(x, D).fraction; hit point =
// x + f D componentwise with the hit coordinate set to the plane
// (kinetic_boundary.cpp:63 order); flight = F f; discarded flight = F - flight.
// The carrier ends at the hit point: its path (and current deposit) is the
// truncated segment from x to the hit point. Incident proper velocity: along D,
// with the kinetic energy KE(u0) + q (E_path . D_trunc), E_path gathered along the
// truncated segment with the deposit's quadrature (the work its current does), so
// the motion work guard closes with the truncated deposit.
#pragma once
#if defined(__CUDACC_RTC__)
#include "push_core.hpp"
#else
#include "push/push_core.hpp"
#endif

namespace hgp {

enum SinkStatus : int {
    sink_field_without_kinetic = 70,   // field surface crossed, no kinetic boundary hit on the represented segment
    lost_kinetic_surface = 71,         // the sealed neutral stream hit a surface that the straight drift does not reach
    sink_off_support = 72,             // (not produced by this version; code number kept)
    sink_unaligned = 73,               // kinetic planes not on field node planes
    sink_unsupported_bank = 74,        // the hit surface has no admitted bank node support
    sink_incident_energy = 75,         // incident kinetic energy not positive and finite
    newborn_start = 76                 // newborn start neither strictly interior nor on one foil face leaving it
};

// Node-index support of every sink kind (0 lower foil, 1 upper foil, 2 outer xy,
// 3 short bridge, 4 axial escape) and the native -> field translation.
struct SinkGeometry {
    Policy p;
    double off[3];
    int lo[5][3], hi[5][3];
    int supported[5];
    int aligned, pad_;
};

struct SinkOut {
    int sink, axis;
    double plane, fraction, flight, discarded, gyro, dcells;
    V3 segment, u, point, field_point;
};

HGP_FN double axis_origin(const Mesh& g, int d) { return d < 2 ? -g.side / 2 : g.zmin; }
HGP_FN double axis_spacing(const Mesh& g, int d) { return d == 0 ? g.dx : (d == 1 ? g.dy : g.dz); }
HGP_FN int axis_count(const Mesh& g, int d) { return d == 0 ? g.nx : (d == 1 ? g.ny : g.nz); }
// Cell coordinate of the field image of a native coordinate: the operations of
// axial::point_deposit / deposit_charge on x_field = x_native + offset.
HGP_FN double cell_of(const Mesh& g, int d, double native, double off) {
    const double field = native + off;
    return (field - axis_origin(g, d)) / axis_spacing(g, d);
}
HGP_FN int nearest_plane(const Mesh& g, int d, double native_plane, double off, int& p) {
    const double v = cell_of(g, d, native_plane, off);
    const double r = HGP_FLOOR(v + .5);
    if (!(HGP_FABS(v - r) <= 1e-8 && r >= 0 && r <= axis_count(g, d))) return sink_unaligned;
    p = int(r);
    return ok;
}
// first_kinetic_hit of push_core.hpp without the interior requirement of its start
// (same operands, operation order and tie rule). Used only for a newborn that starts
// on a foil face and leaves it: the face it leaves reports no hit (leave fraction 0).
HGP_FN int first_kinetic_hit_core(const Policy& p, const V3& x, const V3& D, KineticHit& hit) {
    for (int d = 0; d < 3; ++d) if (!finite(D.x[d])) return kinetic_displacement_nonfinite;
    hit.sink = sink_none; hit.fraction = 1; hit.axis = -1; hit.plane = 0;
    for (int axis = 0; axis < 2; ++axis)
        if (D.x[axis] != 0) {
            const double plane = D.x[axis] > 0 ? (axis ? p.gap : p.width) : 0;
            take_hit(hit, (plane - x.x[axis]) / D.x[axis], sink_outer_xy, axis, plane);
        }
    for (int body = 0; body < 2; ++body) {
        const Foil& foil = body ? p.upper : p.lower;
        double enter = hgp_from_bits(0xfff0000000000000ull), leave = hgp_from_bits(0x7ff0000000000000ull), plane = 0;
        int enter_axis = -1;
        bool possible = true;
        for (int axis = 0; axis < 2; ++axis) {
            const double low = axis ? foil.bottom : foil.left, high = axis ? foil.top : foil.right;
            if (D.x[axis] == 0) { if (x.x[axis] < low || x.x[axis] > high) possible = false; continue; }
            double a = (low - x.x[axis]) / D.x[axis], b = (high - x.x[axis]) / D.x[axis], entry_plane = low;
            if (a > b) { const double t = a; a = b; b = t; entry_plane = high; }
            if (a >= enter) { enter = a; enter_axis = axis; plane = entry_plane; }
            leave = dmin(leave, b);
        }
        if (possible && leave > 0 && enter >= 0 && enter <= leave && enter_axis >= 0)
            take_hit(hit, enter, body ? sink_upper_foil : sink_lower_foil, enter_axis, plane);
    }
    if (D.x[2] != 0) {
        const double plane = D.x[2] > 0 ? p.feed_z : p.short_z, f = (plane - x.x[2]) / D.x[2];
        const double xx = x.x[0] + f * D.x[0], yy = x.x[1] + f * D.x[1];
        const bool bridge = D.x[2] < 0 && xx >= p.lower.left && xx <= p.lower.right && yy >= p.lower.top && yy <= p.upper.bottom;
        take_hit(hit, f, bridge ? sink_short_bridge : sink_axial, 2, plane);
    }
    return ok;
}
HGP_FN int first_hit_from(const Policy& p, const V3& x, const V3& D, int face, KineticHit& hit) {
    return face >= 0 ? first_kinetic_hit_core(p, x, D, hit) : first_kinetic_hit(p, x, D, hit);
}
// Start rule of a newborn (native v6 birth_kinematics): a strictly interior vacuum
// point (body = -1), or a point on exactly one face of one foil, not an edge or a
// corner, inside the kinetic box, with the proper velocity pointing out of the metal
// (body = foil, axis = face normal axis).
HGP_FN int newborn_start_rule(const Policy& p, const V3& x, const V3& u, int& body, int& axis) {
    body = -1; axis = -1;
    if (!(x.x[0] > 0 && x.x[0] < p.width && x.x[1] > 0 && x.x[1] < p.gap && x.x[2] > p.short_z && x.x[2] < p.feed_z)) return newborn_start;
    for (int b = 0; b < 2; ++b) {
        const Foil& f = b ? p.upper : p.lower;
        if (!(x.x[0] >= f.left && x.x[0] <= f.right && x.x[1] >= f.bottom && x.x[1] <= f.top)) continue;
        const bool inx = x.x[0] > f.left && x.x[0] < f.right, iny = x.x[1] > f.bottom && x.x[1] < f.top;
        if (iny && (x.x[0] == f.left || x.x[0] == f.right)) {
            const double out = x.x[0] == f.left ? -1. : 1.;
            if (!(u.x[0] * out > 0)) return newborn_start;
            body = b; axis = 0; return ok;
        }
        if (inx && (x.x[1] == f.bottom || x.x[1] == f.top)) {
            const double out = x.x[1] == f.bottom ? -1. : 1.;
            if (!(u.x[1] * out > 0)) return newborn_start;
            body = b; axis = 1; return ok;
        }
        return newborn_start;   // inside the metal, on an edge or on a corner
    }
    return ok;   // interior vacuum point
}
// x: native start, fx: field start (x + off), u0: proper velocity. Returns ok and
// the truncated motion, or a status code (reference texts in sink_code_message).
HGP_HEAVY int sink_truncate(const Mesh& g, const SinkGeometry& s, const V3& x, const V3& fx, const V3& u0, double q, double m,
                            const double* E, const double* B, double dt, SinkOut& out, int face = -1) {
    PushResult full;
    int code = trajectory(g, fx, u0, q, m, E, B, dt, dt, full);
    if (code) return code;
    if (!(full.gyro <= .05 && full.displacement_cells <= .2)) return gyro_displacement_guard;
    KineticHit kh;
    code = first_hit_from(s.p, x, full.D, face, kh);
    if (code) return code;
    if (kh.sink == sink_none) return sink_field_without_kinetic;
    if (!s.aligned) return sink_unaligned;
    if (!s.supported[kh.sink]) return sink_unsupported_bank;
    V3 point, Dt;
    for (int d = 0; d < 3; ++d) point.x[d] = x.x[d] + kh.fraction * full.D.x[d];
    point.x[kh.axis] = kh.plane;
    for (int d = 0; d < 3; ++d) Dt.x[d] = (point.x[d] + s.off[d]) - fx.x[d];
    V3 Eg, Bg;
    code = path_gather(g, fx, Dt, E, B, Eg, Bg);
    if (code) return code;
    const double K0 = kinetic_energy(m, u0);
    const double K1 = K0 + q * hall_dot(Eg, Dt);
    const double Dn = HGP_SQRT(hall_dot(full.D, full.D));
    if (!(finite(K1) && K1 > 0 && Dn > 0)) return sink_incident_energy;
    const double gm1 = K1 / (m * pic_c * pic_c);
    const double un = pic_c * HGP_SQRT(gm1 * (gm1 + 2));
    if (!finite(un)) return sink_incident_energy;
    for (int d = 0; d < 3; ++d) out.u.x[d] = un * (full.D.x[d] / Dn);
    out.sink = kh.sink;
    out.axis = kh.axis;
    out.plane = kh.plane;
    out.fraction = kh.fraction;
    out.flight = dt * kh.fraction;
    out.discarded = dt - out.flight;
    out.gyro = full.gyro;
    out.dcells = full.displacement_cells;
    out.segment = full.D;
    out.point = point;
    for (int d = 0; d < 3; ++d) out.field_point.x[d] = point.x[d] + s.off[d];
    return ok;
}

// Neutral straight drift over the flight F (prepare_neutral_stream arithmetic,
// face-aware first hit): survivor endpoint x + D or a sink at x + f D (plane
// coordinate), incident u = u0. sink = sink_none for a survivor.
HGP_FN int neutral_drift(const Mesh& g, const SinkGeometry& s, const V3& x, const V3& u, double F, int face, SinkOut& out, V3& end, double& dcells) {
    double u2 = 0;
    for (int d = 0; d < 3; ++d) { if (!finite(u.x[d])) return neutral_momentum_finite; u2 += u.x[d] * u.x[d]; }
    const double gamma = HGP_SQRT(1 + u2 / (pic_c * pic_c));
    if (!finite(gamma)) return neutral_gamma;
    V3 D;
    for (int d = 0; d < 3; ++d) D.x[d] = F * u.x[d] / gamma;
    KineticHit hit;
    const int code = first_hit_from(s.p, x, D, face, hit);
    if (code) return code;
    for (int d = 0; d < 3; ++d) end.x[d] = x.x[d] + hit.fraction * D.x[d];
    out.sink = hit.sink; out.u = u; out.segment = D; out.gyro = 0;
    if (hit.sink == sink_none) { V3 m; for (int d = 0; d < 3; ++d) m.x[d] = end.x[d] - x.x[d]; dcells = HGP_SQRT(hall_dot(m, m)) / min_spacing(g); return ok; }
    end.x[hit.axis] = hit.plane;
    if (!s.aligned) return sink_unaligned;
    if (!s.supported[hit.sink]) return sink_unsupported_bank;
    out.axis = hit.axis; out.plane = hit.plane; out.fraction = hit.fraction; out.flight = F * hit.fraction; out.discarded = F - out.flight;
    out.dcells = HGP_SQRT(hall_dot(D, D)) / min_spacing(g); dcells = out.dcells;
    out.point = end;
    for (int d = 0; d < 3; ++d) out.field_point.x[d] = end.x[d] + s.off[d];
    return ok;
}

}  // namespace hgp

#if !defined(__CUDACC_RTC__)
namespace hgp {
// Reference texts of the codes that the truncation can return (host side).
inline const char* sink_code_message(int code) {
    switch (code) {
        case excessive_path: return "excessive axial path";
        case path_outside_mesh: return "axial path outside mesh";
        case boris_input: return "Boris update requires finite data and positive dt";
        case boris_overflow: return "Boris update overflow";
        case not_converged: return "axial native work path did not converge";
        case gyro_displacement_guard: return "axial full-interval gyro/displacement guard";
        case interior_boundary_owner: return "particle requires native boundary owner";
        case interior_metal: return "particle in/on native metal";
        case kinetic_displacement_nonfinite: return "nonfinite kinetic displacement";
        case gpu_path_capacity: return "device path cut capacity exceeded";
        case sink_field_without_kinetic: return "kinetic sink: field surface crossed without a kinetic boundary hit";
        case lost_kinetic_surface: return "kinetic sink: lost kinetic surface intersection";
        case sink_off_support: return "kinetic sink: hit point outside the conductor/escape node support";
        case sink_unaligned: return "kinetic sink: kinetic boundary planes are not aligned with field node planes";
        case sink_unsupported_bank: return "kinetic sink: the hit surface has no admitted bank node support";
        case sink_incident_energy: return "kinetic sink: incident kinetic energy not positive and finite";
        case newborn_start: return "newborn start neither a strictly interior vacuum point nor a foil-face point leaving the metal";
        case neutral_momentum_finite: return "neutral momentum finite";
        case neutral_gamma: return "neutral momentum outside resolved range";
        default: return "kinetic sink: unknown status";
    }
}
}  // namespace hgp
#endif
