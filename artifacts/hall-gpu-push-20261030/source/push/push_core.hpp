// push_core.hpp -- per-carrier arithmetic of the corrected Hall motion path,
// compiled twice: by NVRTC for the device and by the host compiler (with
// -ffp-contract=off) for CPU-side port checks. Each function is a line-by-line
// transcription of the sealed reference at Hall commit 8dddec0:
//   hall-finite-foil-em-20261029/axial_particles.hpp   Path, point_deposit, intersection, push
//   hall-finite-foil-em-20261029/foil_particles.hpp    gamma, work_kick
//   hall-coupled-transport-20261029/coupled_pic.hpp    hall::dot
//   native-v6.2 mfem_gap_pic.cpp                       pic_boris, pic_kinetic_energy_J
//   hall-canonical-driver-20261029/kinetic_boundary.cpp first_kinetic_hit, require_kinetic_interior,
//                                                       prepare_neutral_stream
//   hall-field-correction-20261030/field_contract.hpp  Words (two_sum, +, *)
//   hall-field-correction-20261030/census_charge.hpp   deposit_charge (net words, x87 gross)
// Operand order and association follow the C++ source exactly; no contraction
// (NVRTC --fmad=false, host -ffp-contract=off). Only + - * / sqrt fma floor fabs
// and comparisons are used, all correctly rounded in binary64.
// The x87 80-bit long double steps of the gross charge are emulated exactly
// (64-bit significand, round to nearest even) by X80 below.
#pragma once

#if defined(__CUDACC_RTC__)
#define HGP_FN __device__ inline
#define HGP_HEAVY __device__ __noinline__
#define HGP_BODY __device__ __forceinline__
typedef unsigned long long hgp_u64;
typedef unsigned int hgp_u32;
#define HGP_SQRT(x) sqrt(x)
#define HGP_FLOOR(x) floor(x)
#define HGP_FABS(x) fabs(x)
#define HGP_FMA(a, b, c) fma(a, b, c)
HGP_FN hgp_u64 hgp_bits(double x) { return (hgp_u64)__double_as_longlong(x); }
HGP_FN double hgp_from_bits(hgp_u64 b) { return __longlong_as_double((long long)b); }
HGP_FN int hgp_clz64(hgp_u64 x) { return __clzll((long long)x); }
#else
#include <bit>
#include <cmath>
#include <cstdint>
#define HGP_FN inline
#define HGP_HEAVY inline
#define HGP_BODY inline
typedef std::uint64_t hgp_u64;
typedef std::uint32_t hgp_u32;
#define HGP_SQRT(x) std::sqrt(x)
#define HGP_FLOOR(x) std::floor(x)
#define HGP_FABS(x) std::fabs(x)
#define HGP_FMA(a, b, c) std::fma(a, b, c)
HGP_FN hgp_u64 hgp_bits(double x) { return std::bit_cast<hgp_u64>(x); }
HGP_FN double hgp_from_bits(hgp_u64 b) { return std::bit_cast<double>(b); }
HGP_FN int hgp_clz64(hgp_u64 x) { return std::countl_zero(x); }
#endif

namespace hgp {

constexpr double pic_c = 299792458.;
constexpr double eps0 = 8.8541878188e-12;           // hall::eps (axial::eps)
constexpr double gauss_node = .5773502691896257645; // axial::Path two-point Gauss abscissa

// ---------------------------------------------------------------- status codes
// Values are stable; host code maps them to the reference exception texts.
enum Code : int {
    ok = 0,
    // first census loop of prepare_canonical_motion
    invalid_carrier = 1, nonfinite_kinematics = 2, interior_boundary_owner = 3, interior_metal = 4,
    unresolved_plasma_operand = 5, charge_outside_mesh = 6,
    // second loop
    neutral_identity = 11, neutral_momentum_finite = 12, neutral_gamma = 13, kinetic_displacement_nonfinite = 14,
    neutral_incident = 15, neutral_first_hit = 16, neutral_displacement = 17,
    excessive_path = 20, path_outside_mesh = 21, boris_input = 22, boris_overflow = 23, not_converged = 24,
    gyro_displacement_guard = 25, lost_surface = 26, boundary_owner = 27, charged_first_hit = 28,
    rounded_endpoint_crosses = 29, unresolved_macrocharge = 30,
    // charge_of_endpoints (corrected free-charge census of the endpoints)
    canonical_charge_operand = 40, unresolved_canonical_charge = 41, deposited_words = 42,
    nonfinite_charge_position = 43, canonical_charge_outside_mesh = 44,
    // device-only capacity/consistency limits (never produced by the reference)
    gpu_path_capacity = 60, gpu_frame_mismatch = 61,
    skipped = 99
};

// ---------------------------------------------------------------- data
struct V3 { double x[3]; };
struct Mesh {
    int nx, ny, nz, nex, ney, nbx, nby, nn, ne, nb;
    double side, depth, dx, dy, dz, zmin;
};
struct Foil { double left, right, bottom, top; };
struct Policy { double width, gap, short_z, feed_z; Foil lower, upper; };

HGP_FN bool finite(double v) { return (hgp_bits(v) & 0x7ff0000000000000ull) != 0x7ff0000000000000ull; }
HGP_FN bool finite3(const V3& v) { return finite(v.x[0]) && finite(v.x[1]) && finite(v.x[2]); }
HGP_FN double dmax(double a, double b) { return (a < b) ? b : a; }   // std::max
HGP_FN double dmin(double a, double b) { return (b < a) ? b : a; }   // std::min
HGP_FN int imin(int a, int b) { return (b < a) ? b : a; }

HGP_FN int node(const Mesh& g, int i, int j, int k) { return (k * (g.ny + 1) + j) * (g.nx + 1) + i; }
HGP_FN int ex(const Mesh& g, int i, int j, int k) { return (k * (g.ny + 1) + j) * g.nx + i; }
HGP_FN int ey(const Mesh& g, int i, int j, int k) { return g.nex + (k * g.ny + j) * (g.nx + 1) + i; }
HGP_FN int ez(const Mesh& g, int i, int j, int k) { return g.nex + g.ney + (k * (g.ny + 1) + j) * (g.nx + 1) + i; }
HGP_FN int bx(const Mesh& g, int i, int j, int k) { return (k * g.ny + j) * (g.nx + 1) + i; }
HGP_FN int by(const Mesh& g, int i, int j, int k) { return g.nbx + (k * (g.ny + 1) + j) * g.nx + i; }
HGP_FN int bz(const Mesh& g, int i, int j, int k) { return g.nbx + g.nby + (k * g.ny + j) * g.nx + i; }
HGP_FN double min_spacing(const Mesh& g) { double m = g.dx; if (g.dy < m) m = g.dy; if (g.dz < m) m = g.dz; return m; }

// hall::dot
HGP_FN double hall_dot(const V3& a, const V3& b) { double s = 0; for (int k = 0; k < 3; ++k) s += a.x[k] * b.x[k]; return s; }
// native norm2 / gamma (mfem_gap_pic.cpp)
HGP_FN double norm2(const V3& a) { return a.x[0] * a.x[0] + a.x[1] * a.x[1] + a.x[2] * a.x[2]; }
HGP_FN double native_gamma(const V3& u) { return HGP_SQRT(1 + norm2(u) / (pic_c * pic_c)); }
// foil::gamma
HGP_FN double foil_gamma(const V3& u) { return HGP_SQRT(1 + hall_dot(u, u) / (pic_c * pic_c)); }
HGP_FN V3 cross(const V3& a, const V3& b) {
    V3 r;
    r.x[0] = a.x[1] * b.x[2] - a.x[2] * b.x[1];
    r.x[1] = a.x[2] * b.x[0] - a.x[0] * b.x[2];
    r.x[2] = a.x[0] * b.x[1] - a.x[1] * b.x[0];
    return r;
}
// pic_kinetic_energy_J
HGP_FN double kinetic_energy(double mass, const V3& u) { return mass * norm2(u) / (native_gamma(u) + 1); }

// pic_boris
HGP_FN int pic_boris(const V3& u, double qm, const V3& E, const V3& B, double dt, V3& result) {
    if (!finite3(u) || !finite3(E) || !finite3(B) || !finite(qm) || !finite(dt) || dt <= 0) return boris_input;
    V3 minus, t, prime, plus, s;
    for (int i = 0; i < 3; ++i) minus.x[i] = u.x[i] + .5 * dt * qm * E.x[i];
    for (int i = 0; i < 3; ++i) t.x[i] = .5 * dt * qm * B.x[i] / native_gamma(minus);
    const V3 mt = cross(minus, t);
    for (int i = 0; i < 3; ++i) { prime.x[i] = minus.x[i] + mt.x[i]; s.x[i] = 2 * t.x[i] / (1 + norm2(t)); }
    const V3 ps = cross(prime, s);
    for (int i = 0; i < 3; ++i) { plus.x[i] = minus.x[i] + ps.x[i]; result.x[i] = plus.x[i] + .5 * dt * qm * E.x[i]; }
    if (!finite3(result) || !finite(norm2(result))) return boris_overflow;
    return ok;
}

// foil::work_kick
HGP_FN int work_kick(const V3& u, double qm, const V3& E, const V3& B, double dt, V3& next, V3& D) {
    const int c = pic_boris(u, qm, E, B, dt, next);
    if (c) return c;
    V3 minus = u, plus = next;
    for (int d = 0; d < 3; ++d) { const double a = .5 * dt * qm * E.x[d]; minus.x[d] += a; plus.x[d] -= a; }
    const double g0 = foil_gamma(u), gm = foil_gamma(minus), gp = foil_gamma(plus), g1 = foil_gamma(next);
    for (int d = 0; d < 3; ++d) D.x[d] = .5 * dt * ((minus.x[d] + u.x[d]) / (gm + g0) + (next.x[d] + plus.x[d]) / (g1 + gp));
    return ok;
}

// ---------------------------------------------------------------- axial::intersection
struct Hit { double fraction; int surface; };
HGP_FN Hit intersection(const Mesh& g, const V3& x, const V3& D) {
    Hit hit = {1, -1};
    const double lo[3] = {-g.side / 2, -g.side / 2, g.zmin};
    const double hi[3] = {g.side / 2, g.side / 2, g.zmin + g.nz * g.dz};
    for (int d = 0; d < 3; ++d)
        if (D.x[d] != 0) {
            const double t = ((D.x[d] > 0 ? hi[d] : lo[d]) - x.x[d]) / D.x[d];
            if (t >= 0 && t <= hit.fraction) { hit.fraction = t; hit.surface = 2; }
        }
    for (int body = 0; body < 2; ++body) {
        const double a[3] = {-.0025, body ? .0025 : -.0027, -g.depth / 2};
        const double b[3] = {.0025, body ? .0027 : -.0025, g.depth / 2};
        double enter = -1e100, leave = 1e100;
        bool possible = true;
        for (int d = 0; d < 3; ++d) {
            if (D.x[d] == 0) {
                if (x.x[d] < a[d] || x.x[d] > b[d]) possible = false;
            } else {
                double u = (a[d] - x.x[d]) / D.x[d], v = (b[d] - x.x[d]) / D.x[d];
                if (u > v) { const double w = u; u = v; v = w; }
                enter = dmax(enter, u);
                leave = dmin(leave, v);
            }
        }
        if (possible && leave > 0 && enter >= 0 && enter <= leave && enter <= hit.fraction) { hit.fraction = enter; hit.surface = body; }
    }
    if (D.x[2] != 0) {
        const double t = (-g.depth / 2 - x.x[2]) / D.x[2];
        if (t >= 0 && t <= hit.fraction && HGP_FABS(x.x[0] + t * D.x[0]) <= .0025 && HGP_FABS(x.x[1] + t * D.x[1]) <= .0027) {
            hit.fraction = t; hit.surface = 3;
        }
    }
    return hit;
}

// ---------------------------------------------------------------- axial::Path
#ifndef HGP_CUT_CAPACITY
#define HGP_CUT_CAPACITY 16
#endif
struct Cuts { int count; double t[HGP_CUT_CAPACITY]; };

HGP_FN int path_cuts(const Mesh& g, const V3& x, const V3& D, Cuts& c) {
    const double spacing[3] = {g.dx, g.dy, g.dz};
    const double origin[3] = {-g.side / 2, -g.side / 2, g.zmin};
    c.count = 0;
    c.t[c.count++] = 0;
    c.t[c.count++] = 1;
    for (int d = 0; d < 3; ++d)
        if (D.x[d] != 0) {
            const double a = (x.x[d] - origin[d]) / spacing[d], b = a + D.x[d] / spacing[d];
            if (!(HGP_FABS(b - a) < 10000)) return excessive_path;
            const double lo = dmin(a, b), hi = dmax(a, b);
            for (int n = int(HGP_FLOOR(lo) + 1); n < hi; ++n) {
                const double t = (n - a) / (b - a);
                if (t > 0 && t < 1) {
                    if (c.count == HGP_CUT_CAPACITY) return gpu_path_capacity;
                    c.t[c.count++] = t;
                }
            }
        }
    for (int i = 1; i < c.count; ++i) {   // ascending; equal values are bitwise equal (all > 0)
        const double v = c.t[i];
        int j = i - 1;
        while (j >= 0 && v < c.t[j]) { c.t[j + 1] = c.t[j]; --j; }
        c.t[j + 1] = v;
    }
    return ok;
}

// Cell and local coordinates of x + t*D (Path constructor inner loop).
HGP_FN int path_point(const Mesh& g, const V3& x, const V3& D, double t, int cell[3], double a[3]) {
    const double spacing[3] = {g.dx, g.dy, g.dz};
    const double origin[3] = {-g.side / 2, -g.side / 2, g.zmin};
    const int counts[3] = {g.nx, g.ny, g.nz};
    for (int d = 0; d < 3; ++d) {
        double v = (x.x[d] + t * D.x[d] - origin[d]) / spacing[d];
        if (!(v >= -1e-10 && v <= counts[d] + 1e-10)) return path_outside_mesh;
        const double top = double(counts[d]);
        v = (v < 0.) ? 0. : ((top < v) ? top : v);
        cell[d] = imin(int(HGP_FLOOR(v)), counts[d] - 1);
        a[d] = v - cell[d];
    }
    return ok;
}

// Gather of E and B along the path (Path ctor + gather_E/gather_B), fused, same order.
HGP_HEAVY int path_gather(const Mesh& g, const V3& x, const V3& D, const double* E, const double* B, V3& Eg, V3& Bg) {
    Cuts c;
    int code = path_cuts(g, x, D, c);
    if (code) return code;
    for (int k = 0; k < 3; ++k) { Eg.x[k] = 0; Bg.x[k] = 0; }
    for (int segment = 1; segment < c.count; ++segment) {
        const double midpoint = .5 * (c.t[segment] + c.t[segment - 1]), span = c.t[segment] - c.t[segment - 1];
        for (int qi = 0; qi < 2; ++qi) {
            const double q = qi ? gauss_node : -gauss_node;
            const double t = midpoint + .5 * span * q, w = .5 * span;
            int cell[3];
            double a[3];
            code = path_point(g, x, D, t, cell, a);
            if (code) return code;
            const int i = cell[0], j = cell[1], k = cell[2];
            for (int r = 0; r < 2; ++r)
                for (int s = 0; s < 2; ++s) {
                    Eg.x[0] += w * (r ? a[1] : 1 - a[1]) * (s ? a[2] : 1 - a[2]) * E[ex(g, i, j + r, k + s)];
                    Eg.x[1] += w * (r ? a[0] : 1 - a[0]) * (s ? a[2] : 1 - a[2]) * E[ey(g, i + r, j, k + s)];
                    Eg.x[2] += w * (r ? a[0] : 1 - a[0]) * (s ? a[1] : 1 - a[1]) * E[ez(g, i + r, j + s, k)];
                }
            for (int r = 0; r < 2; ++r) {
                Bg.x[0] += w * (r ? a[0] : 1 - a[0]) * B[bx(g, i + r, j, k)];
                Bg.x[1] += w * (r ? a[1] : 1 - a[1]) * B[by(g, i, j + r, k)];
                Bg.x[2] += w * (r ? a[2] : 1 - a[2]) * B[bz(g, i, j, k + r)];
            }
        }
    }
    return ok;
}

// Validates every path point (the checks of the Path constructor) and returns the segment count.
HGP_FN int path_check(const Mesh& g, const V3& x, const V3& D, int& segments) {
    Cuts c;
    int code = path_cuts(g, x, D, c);
    if (code) return code;
    for (int segment = 1; segment < c.count; ++segment) {
        const double midpoint = .5 * (c.t[segment] + c.t[segment - 1]), span = c.t[segment] - c.t[segment - 1];
        for (int qi = 0; qi < 2; ++qi) {
            const double q = qi ? gauss_node : -gauss_node;
            const double t = midpoint + .5 * span * q;
            int cell[3];
            double a[3];
            code = path_point(g, x, D, t, cell, a);
            if (code) return code;
        }
    }
    segments = c.count - 1;
    return ok;
}

// ---------------------------------------------------------------- axial::push
struct PushResult {
    V3 u, position, D;
    Hit hit;
    double flight_s, gyro, displacement_cells;
};

HGP_HEAVY int trajectory(const Mesh& g, const V3& x, const V3& u0, double charge, double mass, const double* E, const double* B,
                      double dt, double tau, PushResult& out) {
    V3 D = u0;
    for (int d = 0; d < 3; ++d) D.x[d] *= tau / foil_gamma(u0);
    V3 u = {{0, 0, 0}}, Bg = {{0, 0, 0}};
    bool converged = false;
    for (int iteration = 0; iteration < 12; ++iteration) {
        const Hit hit = intersection(g, x, D);
        V3 gatherD = D;
        if (hit.surface >= 0) for (int d = 0; d < 3; ++d) gatherD.x[d] *= hit.fraction;
        V3 Eg;
        int code = path_gather(g, x, gatherD, E, B, Eg, Bg);
        if (code) return code;
        V3 next, nextD;
        code = work_kick(u0, charge / mass, Eg, Bg, tau, next, nextD);
        if (code) return code;
        double diff = 0;
        for (int d = 0; d < 3; ++d) diff = dmax(diff, HGP_FABS(nextD.x[d] - D.x[d]));
        u = next;
        D = nextD;
        if (diff <= 1e-18) { converged = true; break; }
    }
    if (!converged) return not_converged;
    out.u = u;
    for (int d = 0; d < 3; ++d) out.position.x[d] = x.x[d] + D.x[d];
    out.D = D;
    out.hit = intersection(g, x, D);
    out.flight_s = tau;
    out.gyro = dt * HGP_FABS(charge / mass) * HGP_SQRT(hall_dot(Bg, Bg));
    out.displacement_cells = HGP_SQRT(hall_dot(D, D)) / min_spacing(g);
    return ok;
}

HGP_FN int push(const Mesh& g, const V3& x, const V3& u0, double charge, double mass, const double* E, const double* B,
                double dt, PushResult& result) {
    int code = trajectory(g, x, u0, charge, mass, E, B, dt, dt, result);
    if (code) return code;
    if (!(result.gyro <= .05 && result.displacement_cells <= .2)) return gyro_displacement_guard;
    if (result.hit.surface >= 0) {
        double lo = 0, hi = dt;
        for (int i = 0; i < 48; ++i) {
            const double t = .5 * (lo + hi);
            PushResult trial;
            code = trajectory(g, x, u0, charge, mass, E, B, dt, t, trial);
            if (code) return code;
            if (trial.hit.surface >= 0) hi = t; else lo = t;
        }
        code = trajectory(g, x, u0, charge, mass, E, B, dt, hi, result);
        if (code) return code;
        if (!(result.hit.surface >= 0)) return lost_surface;
        for (int d = 0; d < 3; ++d) { result.D.x[d] *= result.hit.fraction; result.position.x[d] = x.x[d] + result.D.x[d]; }
    }
    return ok;
}

// ---------------------------------------------------------------- kinetic_boundary.cpp
enum Sink : int { sink_none = -1, sink_lower_foil = 0, sink_upper_foil = 1, sink_outer_xy = 2, sink_short_bridge = 3, sink_axial = 4 };
struct KineticHit { int sink; double fraction; int axis; double plane; };

HGP_FN int require_kinetic_interior(const Policy& p, const V3& x) {
    if (!(x.x[0] > 0 && x.x[0] < p.width && x.x[1] > 0 && x.x[1] < p.gap && x.x[2] > p.short_z && x.x[2] < p.feed_z))
        return interior_boundary_owner;
    for (int body = 0; body < 2; ++body) {
        const Foil& f = body ? p.upper : p.lower;
        if (x.x[0] >= f.left && x.x[0] <= f.right && x.x[1] >= f.bottom && x.x[1] <= f.top) return interior_metal;
    }
    return ok;
}

HGP_FN void take_hit(KineticHit& hit, double f, int sink, int axis, double plane) {
    if (f < 0 || f > hit.fraction) return;
    hit.sink = sink; hit.fraction = f; hit.axis = axis; hit.plane = plane;
}

HGP_FN int first_kinetic_hit(const Policy& p, const V3& x, const V3& D, KineticHit& hit) {
    int code = require_kinetic_interior(p, x);
    if (code) return code;
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
            if (a > b) { const double w = a; a = b; b = w; entry_plane = high; }
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

// prepare_neutral_stream (interval and policy checks are host-side and global).
HGP_FN int neutral_stream(const Policy& p, const V3& x, const V3& u, double charge, double mass, double weight, double dt,
                          V3& endpoint, KineticHit& hit) {
    if (!(charge == 0 && finite(weight) && weight > 0 && finite(mass) && mass > 0)) return neutral_identity;
    double u2 = 0;
    for (int d = 0; d < 3; ++d) { if (!finite(u.x[d])) return neutral_momentum_finite; u2 += u.x[d] * u.x[d]; }
    const double gamma = HGP_SQRT(1 + u2 / (pic_c * pic_c));
    if (!finite(gamma)) return neutral_gamma;
    V3 D;
    for (int d = 0; d < 3; ++d) D.x[d] = dt * u.x[d] / gamma;
    const int code = first_kinetic_hit(p, x, D, hit);
    if (code) return code;
    for (int d = 0; d < 3; ++d) endpoint.x[d] = x.x[d] + hit.fraction * D.x[d];
    if (hit.sink != sink_none) endpoint.x[hit.axis] = hit.plane;
    if (!finite(weight * kinetic_energy(mass, u))) return neutral_incident;
    return ok;
}

// ---------------------------------------------------------------- axial::point_deposit geometry
struct NodeWeights { int cell[3]; double a[3]; };
// outside_code: axial point_deposit (charge_outside_mesh) or deposit_charge
// (canonical_charge_outside_mesh, which also checks finiteness per axis first).
HGP_FN int point_cell(const Mesh& g, const V3& x, NodeWeights& w, int outside_code) {
    const double spacing[3] = {g.dx, g.dy, g.dz};
    const double origin[3] = {-g.side / 2, -g.side / 2, g.zmin};
    const int counts[3] = {g.nx, g.ny, g.nz};
    for (int d = 0; d < 3; ++d) {
        if (outside_code == canonical_charge_outside_mesh && !finite(x.x[d])) return nonfinite_charge_position;
        double v = (x.x[d] - origin[d]) / spacing[d];
        if (!(v >= -1e-10 && v <= counts[d] + 1e-10)) return outside_code;
        const double top = double(counts[d]);
        v = (v < 0.) ? 0. : ((top < v) ? top : v);
        w.cell[d] = imin(int(HGP_FLOOR(v)), counts[d] - 1);
        w.a[d] = v - w.cell[d];
    }
    return ok;
}
// Slot s in 0..7 follows the reference loop order: for k, for j, for i.
HGP_FN int slot_node(const Mesh& g, const NodeWeights& w, int s) {
    const int k = (s >> 2) & 1, j = (s >> 1) & 1, i = s & 1;
    return node(g, w.cell[0] + i, w.cell[1] + j, w.cell[2] + k);
}
HGP_FN double slot_value(const NodeWeights& w, int s, double charge) {
    const int k = (s >> 2) & 1, j = (s >> 1) & 1, i = s & 1;
    return charge * (i ? w.a[0] : 1 - w.a[0]) * (j ? w.a[1] : 1 - w.a[1]) * (k ? w.a[2] : 1 - w.a[2]);
}

// ---------------------------------------------------------------- Words (field_contract.hpp)
struct Words { double hi, lo; };
HGP_FN Words two_sum_words(double a, double b) {   // Words(double,double) via axial::two_sum
    const double s = a + b, v = s - a;
    Words w; w.hi = s; w.lo = (a - (s - v)) + (b - v);
    return w;
}
HGP_FN Words words_add(Words a, Words b) {
    const Words first = two_sum_words(a.hi, b.hi);
    return two_sum_words(first.hi, first.lo + a.lo + b.lo);
}
HGP_FN Words words_mul(Words a, double b) {
    const double p = a.hi * b;
    return two_sum_words(p, HGP_FMA(a.hi, b, -p) + a.lo * b);
}
HGP_FN bool words_canonical(Words a) {
    const Words n = two_sum_words(a.hi, a.lo);
    return finite(a.hi) && finite(a.lo) && n.hi == a.hi && n.lo == a.lo;
}
HGP_FN Words slot_words(const NodeWeights& w, int s, Words charge) {
    const int k = (s >> 2) & 1, j = (s >> 1) & 1, i = s & 1;
    return words_mul(words_mul(words_mul(charge, i ? w.a[0] : 1 - w.a[0]), j ? w.a[1] : 1 - w.a[1]), k ? w.a[2] : 1 - w.a[2]);
}

// ---------------------------------------------------------------- x87 extended emulation
// Non-negative values m * 2^e with a 64-bit significand (top bit set) or zero.
struct X80 { hgp_u64 m; int e; };
HGP_FN X80 x80_from_double(double d) {   // d >= 0 and finite
    const hgp_u64 bits = hgp_bits(d) & 0x7fffffffffffffffull;
    const int exponent = int(bits >> 52);
    const hgp_u64 fraction = bits & 0x000fffffffffffffull;
    X80 r;
    if (exponent == 0) {
        if (fraction == 0) { r.m = 0; r.e = 0; return r; }
        const int shift = hgp_clz64(fraction);
        r.m = fraction << shift; r.e = -1074 - shift;
        return r;
    }
    r.m = (fraction | 0x0010000000000000ull) << 11;
    r.e = exponent - 1075 - 11;
    return r;
}
// Round-to-nearest-even sum of two non-negative values to a 64-bit significand.
HGP_FN X80 x80_add(X80 a, X80 b) {
    if (a.m == 0) return b;
    if (b.m == 0) return a;
    if (a.e < b.e) { const X80 t = a; a = b; b = t; }
    const int d = a.e - b.e;
    hgp_u64 bhi, blo;
    bool sticky = false;
    if (d == 0) { bhi = b.m; blo = 0; }
    else if (d < 64) { bhi = b.m >> d; blo = b.m << (64 - d); }
    else if (d == 64) { bhi = 0; blo = b.m; }
    else if (d < 128) { bhi = 0; blo = b.m >> (d - 64); sticky = (b.m << (128 - d)) != 0; }
    else { bhi = 0; blo = 0; sticky = true; }
    // (a.m, 0) + (bhi, blo) as a 129-bit sum
    const hgp_u64 lo = blo;                 // a's low word is zero
    const hgp_u64 hi = a.m + bhi;
    const bool carry = hi < a.m;
    X80 r;
    bool guard, rest;
    hgp_u64 m;
    if (carry) {
        m = 0x8000000000000000ull | (hi >> 1);
        guard = (hi & 1) != 0;
        rest = lo != 0 || sticky;
        r.e = a.e + 1;
    } else {
        m = hi;
        guard = (lo >> 63) != 0;
        rest = (lo << 1) != 0 || sticky;
        r.e = a.e;
    }
    if (guard && (rest || (m & 1))) {
        m += 1;
        if (m == 0) { m = 0x8000000000000000ull; r.e += 1; }
    }
    r.m = m;
    return r;
}
// Conversion to binary64: rounding to nearest even, then optionally the
// reference's "if the double is below the long double, nextafter(+inf)" step.
HGP_FN double x80_to_double(X80 x, bool upward_if_below) {
    if (x.m == 0) return 0.0;
    const int top = x.e + 63;   // value in [2^top, 2^(top+1))
    hgp_u64 q, rem, half;
    hgp_u64 high_part;
    if (top >= -1022) {
        q = x.m >> 11; rem = x.m & 0x7ffull; half = 0x400ull;
        high_part = hgp_u64(top + 1022) << 52;   // adding q (with its hidden bit) yields the biased exponent
    } else {
        const int shift = -(x.e + 1074);       // >= 12
        high_part = 0;
        if (shift > 64) { q = 0; rem = x.m; half = 0; }   // below half of the smallest subnormal
        else if (shift == 64) { q = 0; rem = x.m; half = 0x8000000000000000ull; }
        else { q = x.m >> shift; rem = x.m & ((hgp_u64(1) << shift) - 1); half = hgp_u64(1) << (shift - 1); }
    }
    bool up;
    if (half == 0) up = false;
    else up = rem > half || (rem == half && (q & 1));
    hgp_u64 bits = high_part + q + (up ? 1 : 0);
    if (bits >= 0x7ff0000000000000ull) bits = 0x7ff0000000000000ull;   // overflow to +inf
    if (upward_if_below && rem != 0 && !up && bits < 0x7ff0000000000000ull) bits += 1;
    return hgp_from_bits(bits);
}
// census_charge.hpp deposit_charge gross step:
//   long double gross = (long double)g + (|hi| + |lo| in long double);
//   g = double(gross); if ((long double)g < gross) g = nextafter(g, +inf);
HGP_FN double gross_step(double g, Words v) {
    const X80 magnitude = x80_add(x80_from_double(HGP_FABS(v.hi)), x80_from_double(HGP_FABS(v.lo)));
    return x80_to_double(x80_add(x80_from_double(g), magnitude), true);
}

}  // namespace hgp
