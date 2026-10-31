// relocation.cpp -- relocation current and Hall continuity check (see relocation.hpp).
#include "relocation.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace merge_split {
namespace {
void need(bool ok, const char* message) {
    if (!ok) throw Error(message);
}
// axial_fields.hpp:110 two_sum
std::pair<double, double> two_sum(double a, double b) {
    double s = a + b, v = s - a;
    return {s, (a - (s - v)) + (b - v)};
}
long double magnitude(HWords a) { return std::abs((long double)a.hi) + std::abs((long double)a.lo); }
}  // namespace

HWords::HWords(double a, double b) {
    auto p = two_sum(a, b);
    hi = p.first;
    lo = p.second;
}
HWords operator+(HWords a, HWords b) {
    auto first = two_sum(a.hi, b.hi);
    return HWords(first.first, first.second + a.lo + b.lo);
}
HWords operator-(HWords a) { return HWords(-a.hi, -a.lo); }
HWords operator-(HWords a, HWords b) { return a + (-b); }
HWords operator*(HWords a, double b) {
    const double p = a.hi * b;
    return HWords(p, std::fma(a.hi, b, -p) + a.lo * b);
}
HWords operator/(HWords a, double b) {
    const double q = a.hi / b;
    const auto r = a - HWords(q) * b;
    return HWords(q, (r.hi + r.lo) / b);
}

HallMesh::HallMesh(int n, int axial_cells, double width, int air_cells, double depth_m)
    : nx(n), ny(n), nz(axial_cells + 2 * air_cells), axial(axial_cells), air(air_cells),
      nn((n + 1) * (n + 1) * (nz + 1)),
      ne(n * (n + 1) * (nz + 1) + (n + 1) * n * (nz + 1) + (n + 1) * (n + 1) * nz),
      nex(n * (n + 1) * (nz + 1)), ney((n + 1) * n * (nz + 1)), side(width), depth(depth_m), dx(width / n),
      dy(width / n), dz(depth_m / axial_cells), zmin(-depth_m / 2 - air_cells * dz) {
    need(n >= 9 && n <= 500 && axial_cells >= 2 && air_cells >= 1 && width > .0054, "invalid axial mesh");
}

namespace {
struct Weight {
    int index, axis;
    double value;
};
// axial::Path electric weights (axial_particles.hpp:8-31).
void path_weights(const HallMesh& g, const double x[3], const double D[3], std::vector<double>& cuts, std::vector<Weight>& electric) {
    const double spacing[3] = {g.dx, g.dy, g.dz}, origin[3] = {-g.side / 2, -g.side / 2, g.zmin};
    const int counts[3] = {g.nx, g.ny, g.nz};
    cuts.assign({0, 1});
    electric.clear();
    for (int d = 0; d < 3; ++d)
        if (D[d] != 0) {
            double a = (x[d] - origin[d]) / spacing[d], b = a + D[d] / spacing[d];
            need(std::abs(b - a) < 10000, "excessive axial path");
            for (int n = std::floor(std::min(a, b)) + 1; n < std::max(a, b); ++n) {
                double t = (n - a) / (b - a);
                if (t > 0 && t < 1) cuts.push_back(t);
            }
        }
    std::sort(cuts.begin(), cuts.end());
    for (std::size_t segment = 1; segment < cuts.size(); ++segment) {
        double midpoint = .5 * (cuts[segment] + cuts[segment - 1]), span = cuts[segment] - cuts[segment - 1];
        for (double q : {-.5773502691896257645, .5773502691896257645}) {
            double t = midpoint + .5 * span * q, w = .5 * span;
            int cell[3];
            double a[3];
            for (int d = 0; d < 3; ++d) {
                double v = (x[d] + t * D[d] - origin[d]) / spacing[d];
                need(v >= -1e-10 && v <= counts[d] + 1e-10, "axial path outside mesh");
                v = std::clamp(v, 0., double(counts[d]));
                cell[d] = std::min(int(std::floor(v)), counts[d] - 1);
                a[d] = v - cell[d];
            }
            int i = cell[0], j = cell[1], k = cell[2];
            for (int r : {0, 1})
                for (int s : {0, 1}) {
                    electric.push_back({g.ex(i, j + r, k + s), 0, w * (r ? a[1] : 1 - a[1]) * (s ? a[2] : 1 - a[2])});
                    electric.push_back({g.ey(i + r, j, k + s), 1, w * (r ? a[0] : 1 - a[0]) * (s ? a[2] : 1 - a[2])});
                    electric.push_back({g.ez(i + r, j + s, k), 2, w * (r ? a[0] : 1 - a[0]) * (s ? a[1] : 1 - a[1])});
                }
        }
    }
}
// deposit_charge (census_charge.hpp:15-26)
void deposit_charge(const HallMesh& g, FreeChargeW& q, const double xin[3], HWords charge) {
    const double spacing[3] = {g.dx, g.dy, g.dz}, origin[3] = {-g.side / 2, -g.side / 2, g.zmin};
    const int counts[3] = {g.nx, g.ny, g.nz};
    int cell[3];
    double a[3];
    for (int d = 0; d < 3; ++d) {
        need(std::isfinite(xin[d]), "nonfinite charge position");
        double v = (xin[d] - origin[d]) / spacing[d];
        need(v >= -1e-10 && v <= counts[d] + 1e-10, "canonical charge outside mesh");
        v = std::clamp(v, 0., double(counts[d]));
        cell[d] = std::min(int(std::floor(v)), counts[d] - 1);
        a[d] = v - cell[d];
    }
    for (int k : {0, 1})
        for (int j : {0, 1})
            for (int i : {0, 1}) {
                const int n = g.node(cell[0] + i, cell[1] + j, cell[2] + k);
                const auto v = charge * (i ? a[0] : 1 - a[0]) * (j ? a[1] : 1 - a[1]) * (k ? a[2] : 1 - a[2]);
                q.net[n] = q.net[n] + v;
                const long double gross = (long double)q.gross[n] + magnitude(v);
                q.gross[n] = double(gross);
                if ((long double)q.gross[n] < gross) q.gross[n] = std::nextafter(q.gross[n], INFINITY);
            }
}
}  // namespace

std::vector<double> relocation_current(const HallMesh& g, const Relocation& r, double dt) {
    need(std::isfinite(dt) && dt > 0, "relocation current interval");
    std::vector<double> J(std::size_t(g.ne), 0.);
    std::vector<double> cuts;
    std::vector<Weight> electric;
    for (std::size_t i = 0; i < r.size(); ++i) {
        const double x[3] = {r.start[0][i], r.start[1][i], r.start[2][i]}, D[3] = {r.D[0][i], r.D[1][i], r.D[2][i]};
        path_weights(g, x, D, cuts, electric);
        const double charge = r.charge[i];
        for (const auto& w : electric) J[std::size_t(w.index)] += charge * D[w.axis] * w.value / dt;
    }
    return J;
}

FreeChargeW free_charge(const HallMesh& g, const MsParams& P, const Carriers& C) {
    FreeChargeW out{std::vector<HWords>(std::size_t(g.nn)), std::vector<double>(std::size_t(g.nn), 0.)};
    std::vector<std::size_t> order(C.size());
    std::iota(order.begin(), order.end(), std::size_t(0));
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return C.id[a] < C.id[b]; });
    for (auto i : order) {
        const double charge_C = P.species[C.species[i]].charge_C;
        double x[3] = {C.x[0][i], C.x[1][i], C.x[2][i]};
        x[0] -= g.side / 2;   // charge_of_particles: x[0] -= width/2, x[1] -= gap/2 (width = gap = side)
        x[1] -= g.side / 2;
        const auto q = HWords(C.w[i]) * charge_C;
        deposit_charge(g, out, x, q);
    }
    return out;
}

ContinuityReport relocation_continuity(const HallMesh& g, const MsParams& P, const Carriers& before, const Carriers& after,
                                       const std::vector<double>& J, double dt) {
    need(J.size() == std::size_t(g.ne), "relocation current shape");
    const auto b = free_charge(g, P, before), a = free_charge(g, P, after);
    // incidence (field_contract.hpp): q = Words(v[e])/length; out[a] -= q; out[b] += q, edges ascending.
    std::vector<HWords> div(std::size_t(g.nn));
    auto edge = [&](int e, int na, int nb, double length) {
        const auto q = HWords(J[std::size_t(e)]) / length;
        div[std::size_t(na)] = div[std::size_t(na)] - q;
        div[std::size_t(nb)] = div[std::size_t(nb)] + q;
    };
    for (int k = 0; k <= g.nz; ++k)
        for (int j = 0; j <= g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) edge(g.ex(i, j, k), g.node(i, j, k), g.node(i + 1, j, k), g.dx);
    for (int k = 0; k <= g.nz; ++k)
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i <= g.nx; ++i) edge(g.ey(i, j, k), g.node(i, j, k), g.node(i, j + 1, k), g.dy);
    for (int k = 0; k < g.nz; ++k)
        for (int j = 0; j <= g.ny; ++j)
            for (int i = 0; i <= g.nx; ++i) edge(g.ez(i, j, k), g.node(i, j, k), g.node(i, j, k + 1), g.dz);
    ContinuityReport rep;
    for (int n = 0; n < g.nn; ++n) {
        const auto change = a.net[std::size_t(n)] - b.net[std::size_t(n)];
        const auto defect = change - div[std::size_t(n)] * dt;
        const long double scale = (long double)b.gross[std::size_t(n)] + a.gross[std::size_t(n)] + magnitude(div[std::size_t(n)] * dt);
        const long double size = std::abs(defect.wide());
        rep.max_defect_C = std::max(rep.max_defect_C, double(size));
        if (size > 0) rep.max_ratio = std::max(rep.max_ratio, double(size / (1e-12L * scale)));
        const double moved = double(std::abs(change.wide()));
        rep.max_net_change_C = std::max(rep.max_net_change_C, moved);
        if (change.hi != 0 || change.lo != 0) ++rep.nodes_changed;
    }
    return rep;
}

double relocation_work_J(const std::vector<double>& E, const std::vector<double>& J, double dt) {
    need(E.size() == J.size(), "relocation work shapes");
    double w = 0;
    for (std::size_t e = 0; e < E.size(); ++e) w += dt * E[e] * J[e];
    return w;
}

std::vector<double> test_edge_field(const HallMesh& g, const double E0[3], const double G[3][3]) {
    std::vector<double> E(std::size_t(g.ne));
    auto value = [&](int c, double x, double y, double z) { return E0[c] + G[c][0] * x + G[c][1] * y + G[c][2] * z; };
    const double o0 = -g.side / 2, o1 = -g.side / 2, o2 = g.zmin;
    for (int k = 0; k <= g.nz; ++k)
        for (int j = 0; j <= g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) E[std::size_t(g.ex(i, j, k))] = value(0, o0 + (i + .5) * g.dx, o1 + j * g.dy, o2 + k * g.dz);
    for (int k = 0; k <= g.nz; ++k)
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i <= g.nx; ++i) E[std::size_t(g.ey(i, j, k))] = value(1, o0 + i * g.dx, o1 + (j + .5) * g.dy, o2 + k * g.dz);
    for (int k = 0; k < g.nz; ++k)
        for (int j = 0; j <= g.ny; ++j)
            for (int i = 0; i <= g.nx; ++i) E[std::size_t(g.ez(i, j, k))] = value(2, o0 + i * g.dx, o1 + j * g.dy, o2 + (k + .5) * g.dz);
    return E;
}

}  // namespace merge_split
