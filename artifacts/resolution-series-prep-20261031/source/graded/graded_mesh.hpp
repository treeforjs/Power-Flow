#pragma once
// Graded transverse Yee mesh with the Hall axial-mesh conventions (resolution-series prep).
//
// The sealed Hall mesh (hall-finite-foil-em-20261029/axial_mesh.hpp, 8dddec0) is uniform in x and y. This structure is the
// same complex with per-index transverse node coordinates: a uniform core of cell size h over the 9 mm kinetic box (foils,
// gap, every region carriers can reach), then cells growing geometrically by q per cell up to the box half-width; axial
// cells uniform as in Hall (4 foil cells + air). Every operator array has the meaning and the construction order of
// axial::Mesh (edge/face/node numbering, M = eps0 x dual edge volume, H = dual face volume / mu0, Y = boundary area / Z0,
// curl rows with coefficients +-1/length, PEC short, foil bodies, outer nodes), with the local spacing in place of dx, dy.
// With q = 1 and the 9 mm box it reproduces axial::Mesh(45, 4, .009, 1) (checked bit for bit in graded-dc mode).
// Hall v3.1 cannot run this mesh (its Model takes a uniform MeshSpec); it is used here for sizes, memory and the DC port
// inductance (the dc_initial uniform-sigma path, same operation order, with these operators).
#include "hall-finite-foil-em-20261029/axial_mesh.hpp"
#include <vector>

namespace graded {
using Vec = std::vector<double>;
// Transverse grid: a uniform core of ncells = round(2 core / h) cells of size dx = (2 core) / ncells with node coordinates
// i dx - core (the arithmetic of axial::Mesh for the core), then on each side cells h q, h q^2, ... (capped at hmax) until the
// box half-width; the last cell ends on the box (it absorbs the remainder when that is at least half a grown cell, otherwise
// the remainder is merged into it). Returns node coordinates and the cell sizes used for lengths and volumes.
struct Axis { Vec x, hcell; };
inline Axis axis(double h, double core, double box, double q, double hmax) {
    const int ncells = int(std::lround(2 * core / h));
    hall::require(std::abs(ncells * h - 2 * core) < 1e-12 && box >= core, "graded core must be a whole number of cells");
    const double dx = (2 * core) / ncells;
    std::vector<double> grow;   // outer cell sizes from the core outward
    double reach = core, step = dx;
    while (reach < box - 1e-12) {
        step = std::min(step * q, hmax);
        double next = reach + step;
        if (box - next < .5 * step) { grow.push_back(box - reach); reach = box; break; }
        grow.push_back(step); reach = next;
    }
    Axis a;
    for (std::size_t k = grow.size(); k-- > 0;) a.hcell.push_back(grow[k]);
    for (int i = 0; i < ncells; ++i) a.hcell.push_back(dx);
    for (double g : grow) a.hcell.push_back(g);
    // coordinates: core nodes exactly i dx - core; outer nodes by cumulative sums from the core edges
    const int nl = int(grow.size());
    a.x.assign(a.hcell.size() + 1, 0.);
    for (int i = 0; i <= ncells; ++i) a.x[nl + i] = i * dx - core;
    for (int k = 0; k < nl; ++k) { a.x[nl - 1 - k] = a.x[nl - k] - grow[k]; a.x[nl + ncells + 1 + k] = a.x[nl + ncells + k] + grow[k]; }
    if (nl) { a.x.front() = -box; a.x.back() = box; }
    return a;
}

struct Mesh {
    int nx, ny, nz, nz_foil, kshort, kfeed, nn, ne, nb, nex, ney, nbx, nby;
    double side, depth = .025, dz, zmin, h;
    Vec xs, ys, hx, hy;   // transverse node coordinates and cell sizes
    Vec M, H, Y, VE, VB, VN;
    std::vector<axial::Edge> edges; std::vector<foil::CurlRow> C; std::vector<axial::Cell> cells; std::vector<int> node_body;
    std::vector<bool> short_node, fixed_E, fixed_A, outer_node;
    double dxi(int i) const { return hx[i]; }
    double dyj(int j) const { return hy[j]; }
    // dual lengths (the uniform case gives .5 h at the ends and h inside, as axial::Mesh::dual)
    double dualx(int i) const { return i == 0 ? .5 * hx[0] : i == nx ? .5 * hx[nx - 1] : (hx[i - 1] == hx[i] ? hx[i] : .5 * hx[i - 1] + .5 * hx[i]); }
    double dualy(int j) const { return j == 0 ? .5 * hy[0] : j == ny ? .5 * hy[ny - 1] : (hy[j - 1] == hy[j] ? hy[j] : .5 * hy[j - 1] + .5 * hy[j]); }
    double dualz(int k) const { return (k == 0 || k == nz) ? .5 * dz : dz; }
    Mesh(const Axis& a, int axial_cells, int air_cells)
        : nx(int(a.x.size()) - 1), ny(int(a.x.size()) - 1), nz(axial_cells + 2 * air_cells), nz_foil(axial_cells), kshort(air_cells),
          kfeed(air_cells + axial_cells), side(a.x.back() - a.x.front()), dz(depth / axial_cells), zmin(-depth / 2 - air_cells * dz),
          h(*std::min_element(a.hcell.begin(), a.hcell.end())), xs(a.x), ys(a.x), hx(a.hcell), hy(a.hcell) {
        const int n = nx;
        nn = (n + 1) * (n + 1) * (nz + 1); ne = n * (n + 1) * (nz + 1) + (n + 1) * n * (nz + 1) + (n + 1) * (n + 1) * nz;
        nb = (n + 1) * n * nz + n * (n + 1) * nz + n * n * (nz + 1); nex = n * (n + 1) * (nz + 1); ney = (n + 1) * n * (nz + 1);
        nbx = (n + 1) * n * nz; nby = n * (n + 1) * nz;
        M.assign(ne, 0); H.assign(nb, 0); Y.assign(ne, 0); VE.assign(ne, 0); VB.assign(nb, 0); VN.assign(nn, 0);
        edges.resize(ne); C.resize(nb); node_body.assign(nn, -1); short_node.assign(nn, false); fixed_E.assign(ne, false); fixed_A.assign(ne, false);
        outer_node.assign(nn, false);
        for (double f : {-.0027, -.0025, .0025, .0027}) {
            bool found = false; for (double c : xs) found = found || std::abs(c - f) < 1e-12;
            hall::require(found, "unaligned physical foil (graded)");
        }
        using foil::eps; using foil::mu; using foil::Z0;
        for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i) {
            int n0 = node(i, j, k); double xx = xs[i], yy = ys[j];
            VN[n0] = dualx(i) * dualy(j) * dualz(k);
            outer_node[n0] = i == 0 || i == nx || j == 0 || j == ny || k == 0 || k == nz;
            bool inx = std::abs(xx) <= .0025 + 1e-12;
            if (inx && k >= kshort && k <= kfeed) {
                if (yy >= -.0027 - 1e-12 && yy <= -.0025 + 1e-12) node_body[n0] = 0;
                if (yy >= .0025 - 1e-12 && yy <= .0027 + 1e-12) node_body[n0] = 1;
            }
            short_node[n0] = k == kshort && inx && std::abs(yy) <= .0027 + 1e-12;
        }
        auto edge = [&](int e, int a, int b, int axis, double length, double volume, bool outer, double boundary_area) {
            edges[e] = {a, b, axis, length}; VE[e] = volume; M[e] = eps * volume; Y[e] = boundary_area / Z0; fixed_A[e] = outer;
            fixed_E[e] = axis != 2 && short_node[a] && short_node[b];
        };
        for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i < nx; ++i) {
            double ay = dualy(j), az = dualz(k), dx = dxi(i);
            edge(ex(i, j, k), node(i, j, k), node(i + 1, j, k), 0, dx, dx * ay * az, j == 0 || j == ny || k == 0 || k == nz,
                 (j == 0 || j == ny ? dx * az : 0) + (k == 0 || k == nz ? dx * ay : 0));
        }
        for (int k = 0; k <= nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i <= nx; ++i) {
            double ax = dualx(i), az = dualz(k), dy = dyj(j);
            edge(ey(i, j, k), node(i, j, k), node(i, j + 1, k), 1, dy, ax * dy * az, i == 0 || i == nx || k == 0 || k == nz,
                 (i == 0 || i == nx ? dy * az : 0) + (k == 0 || k == nz ? ax * dy : 0));
        }
        for (int k = 0; k < nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i) {
            double ax = dualx(i), ay = dualy(j);
            edge(ez(i, j, k), node(i, j, k), node(i, j, k + 1), 2, dz, ax * ay * dz, i == 0 || i == nx || j == 0 || j == ny,
                 (i == 0 || i == nx ? ay * dz : 0) + (j == 0 || j == ny ? ax * dz : 0));
        }
        for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i <= nx; ++i) {
            int b = bx(i, j, k); double dy = dyj(j); VB[b] = dualx(i) * dy * dz;
            C[b] = {{ez(i, j + 1, k), ez(i, j, k), ey(i, j, k + 1), ey(i, j, k)}, {1 / dy, -1 / dy, -1 / dz, 1 / dz}, 4};
        }
        for (int k = 0; k < nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i < nx; ++i) {
            int b = by(i, j, k); double dx = dxi(i); VB[b] = dx * dualy(j) * dz;
            C[b] = {{ex(i, j, k + 1), ex(i, j, k), ez(i + 1, j, k), ez(i, j, k)}, {1 / dz, -1 / dz, -1 / dx, 1 / dx}, 4};
        }
        for (int k = 0; k <= nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            int b = bz(i, j, k); double dx = dxi(i), dy = dyj(j); VB[b] = dx * dy * dualz(k);
            C[b] = {{ey(i + 1, j, k), ey(i, j, k), ex(i, j + 1, k), ex(i, j, k)}, {1 / dx, -1 / dx, -1 / dy, 1 / dy}, 4};
        }
        for (int b = 0; b < nb; ++b) H[b] = VB[b] / mu;
        for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            double x = .5 * (xs[i] + xs[i + 1]), y = .5 * (ys[j] + ys[j + 1]); int body = -1;
            if (k >= kshort && k < kfeed && std::abs(x) < .0025 && std::abs(y) > .0025 && std::abs(y) < .0027) body = y > 0 ? 1 : 0;
            cells.push_back({{ex(i, j, k), ex(i, j + 1, k), ex(i, j, k + 1), ex(i, j + 1, k + 1), ey(i, j, k), ey(i + 1, j, k), ey(i, j, k + 1),
                              ey(i + 1, j, k + 1), ez(i, j, k), ez(i + 1, j, k), ez(i, j + 1, k), ez(i + 1, j + 1, k)},
                             i, j, k, body, dxi(i) * dyj(j) * dz});
        }
    }
    int node(int i, int j, int k) const { return (k * (ny + 1) + j) * (nx + 1) + i; }
    int ex(int i, int j, int k) const { return (k * (ny + 1) + j) * nx + i; }
    int ey(int i, int j, int k) const { return nex + (k * ny + j) * (nx + 1) + i; }
    int ez(int i, int j, int k) const { return nex + ney + (k * (ny + 1) + j) * (nx + 1) + i; }
    int bx(int i, int j, int k) const { return (k * ny + j) * (nx + 1) + i; }
    int by(int i, int j, int k) const { return nbx + (k * (ny + 1) + j) * nx + i; }
    int bz(int i, int j, int k) const { return nbx + nby + (k * ny + j) * nx + i; }
    int index_of(const Vec& c, double v) const { for (std::size_t i = 0; i < c.size(); ++i) if (std::abs(c[i] - v) < 1e-12) return int(i); return -1; }
    Vec gradient(const Vec& p) const { Vec v(ne); for (int e = 0; e < ne; ++e) v[e] = (p[edges[e].b] - p[edges[e].a]) / edges[e].length; return v; }
    Vec gradient_transpose(const Vec& v) const { Vec p(nn); for (int e = 0; e < ne; ++e) { double q = v[e] / edges[e].length; p[edges[e].a] -= q; p[edges[e].b] += q; } return p; }
    Vec curl(const Vec& E) const { Vec B(nb); for (int b = 0; b < nb; ++b) for (int k = 0; k < 4; ++k) B[b] += C[b].coefficient[k] * E[C[b].index[k]]; return B; }
    Vec curl_transpose(const Vec& B) const { Vec E(ne); for (int b = 0; b < nb; ++b) for (int k = 0; k < 4; ++k) E[C[b].index[k]] += C[b].coefficient[k] * B[b]; return E; }
    Vec curl_adjoint(Vec B) const { for (int b = 0; b < nb; ++b) B[b] *= H[b]; return curl_transpose(B); }
    Vec conduction(const Vec& sigma_cells) const {
        Vec s(ne); for (std::size_t c = 0; c < cells.size(); ++c) for (int e : cells[c].edges) s[e] += .25 * cells[c].volume * sigma_cells[c]; return s; }
    // Same quadrature as axial::Mesh::feed_shape over the (uniform-core) foil contact at the feed plane.
    Vec feed_shape() const {
        Vec f(ne); double weight = 0;
        const int i0 = index_of(xs, -.0025), i1 = index_of(xs, .0025), j0 = index_of(ys, .0025), j1 = index_of(ys, .0027);
        hall::require(i0 >= 0 && i1 >= 0 && j0 >= 0 && j1 >= 0, "feed contact not on nodes");
        for (int j = j0; j <= j1; ++j) for (int i = i0; i <= i1; ++i) {
            const double w = ((i == i0 || i == i1) ? .5 : 1) * ((j == j0 || j == j1) ? .5 : 1) * hx[i0] * hy[j0] / (.005 * .0002);
            weight += w; for (int y = ny - j; y < j; ++y) f[ey(i, y, kfeed)] -= w * dyj(y);
        }
        hall::require(std::abs(weight - 1) <= 2e-14, "feed contact quadrature normalization (graded)");
        return f;
    }
};
}  // namespace graded
