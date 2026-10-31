#pragma once
// Deterministic geometric multigrid (GMG) preconditioner for the Hall implicit edge operator (resolution-series prep).
//
// Operator (Hall corrective v3.1 field_model.hpp, solver_backend.hpp): A = M^-1/2 K M^-1/2 on free edges, identity on PEC
// edges, with K = M + (dt/2)(S + Y) + (dt^2/4) C^T H C (M edge permittivity mass, S conduction, Y outer admittance, C the
// Yee curl, H face reluctance mass). The outer solver is CG with the same stopping rule as the Hall gpu1 CG (recursive
// residual <= 0.2 target, true residual recomputed, at most three correction cycles, target = max(1e-30, 2e-13 ||b||));
// acceptance is the independent check of Hall v3.1 (HostApply, bitwise the sealed apply; foil::norm, the sealed long-double
// norm).
//
// Preconditioner: one symmetric V-cycle on K.
//   hierarchy   semi-coarsening in x and y only (factor 2, or 3 when the transverse cell count is odd), the axial cells are
//               kept (the Hall mesh has dz/dx = 31-250, the strong couplings are transverse); lowest-order edge-element
//               prolongation (x-edge values constant along x within a coarse cell and linear in y, z-edges bilinear in x
//               and y), which commutes with the discrete gradient (P_e G_c = G_f P_n); PEC edges have zero rows of P;
//               Galerkin coarse operators P^T K P, kept as three parts (mass, conduction+admittance, curl-curl) so that a dt
//               change is a linear combination with the same sparsity.
//   smoother    Hiptmair hybrid smoother: l1-Jacobi on K, then l1-Jacobi on the nodal operator G^T (M + dt/2 (S+Y)) G mapped
//               back with G (the gradient fields are the near-kernel of the curl-curl term that edge Jacobi cannot smooth);
//               pre-smoothing (edge, node), post-smoothing (node, edge), so the V-cycle is symmetric and positive definite;
//               coarsest level: nu_c (edge, node) sweeps followed by nu_c (node, edge) sweeps.
//   determinism every product is a CSR row loop in ascending column order (gmg_core.h), compiled identically for the
//               host reference and NVRTC (--fmad=false --prec-div=true --prec-sqrt=true --ftz=false); dot products use a
//               fixed two-stage tree (TPB 256, NBLK 512) that the host reference reproduces; no atomics. The CPU reference and
//               the GPU produce the same bits.
#include "gmg_core.h"
#include "hall-finite-foil-em-20261029/axial_mesh.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace gmg {
using Vec = std::vector<double>;
using Bytes = std::vector<unsigned char>;
struct Csr { int n = 0, m = 0; std::vector<int> rp{0}, ci; Vec v; std::size_t nnz() const { return ci.size(); } };
struct Csr3 { int n = 0; std::vector<int> rp{0}, ci; Vec a, b, c; std::size_t nnz() const { return ci.size(); } };
inline void need(bool ok, const char* m) { if (!ok) throw std::runtime_error(m); }

struct Grid {
    int n = 0, nz = 0; double dx = 0, dz = 0;
    int nn = 0, nex = 0, ney = 0, ne = 0;
    Grid() = default;
    Grid(int n_, int nz_, double dx_, double dz_) : n(n_), nz(nz_), dx(dx_), dz(dz_) {
        nn = (n + 1) * (n + 1) * (nz + 1); nex = n * (n + 1) * (nz + 1); ney = (n + 1) * n * (nz + 1); ne = nex + ney + (n + 1) * (n + 1) * nz;
    }
    int node(int i, int j, int k) const { return (k * (n + 1) + j) * (n + 1) + i; }
    int ex(int i, int j, int k) const { return (k * (n + 1) + j) * n + i; }
    int ey(int i, int j, int k) const { return nex + (k * n + j) * (n + 1) + i; }
    int ez(int i, int j, int k) const { return nex + ney + (k * (n + 1) + j) * (n + 1) + i; }
};

struct Level {
    Grid g;
    Bytes fixed;               // per edge: excluded (identity row) at this level
    Bytes fixed_nodes;         // per node: phi = 0 in the auxiliary nodal space (PEC nodes)
    Csr3 K3;                   // a: mass, b: conduction + outer admittance, c: curl-curl (H-weighted)
    Csr K; Vec invd;           // combined for the current dt, l1-Jacobi inverse diagonal
    Csr G, GT;                 // gradient (edges x nodes; zero rows at fixed edges) and transpose
    Csr3 N3; Csr N; Vec invdn; // nodal auxiliary operator G^T K G (all three parts, so its l1 diagonal bounds G^T K G)
    Csr P, PT;                 // prolongation from level+1 (rows: this level's edges) and its transpose
};

inline Csr transpose(const Csr& A) {
    Csr T; T.n = A.m; T.m = A.n; T.rp.assign(A.m + 1, 0);
    for (int c : A.ci) ++T.rp[c + 1];
    for (int i = 0; i < A.m; ++i) T.rp[i + 1] += T.rp[i];
    T.ci.resize(A.nnz()); T.v.resize(A.nnz());
    std::vector<int> fill(T.rp.begin(), T.rp.end() - 1);
    for (int r = 0; r < A.n; ++r) for (int k = A.rp[r]; k < A.rp[r + 1]; ++k) { const int c = A.ci[k]; T.ci[fill[c]] = r; T.v[fill[c]] = A.v[k]; ++fill[c]; }
    return T;   // rows of T in ascending column order because rows of A are visited in ascending order
}

// Galerkin product R X P (R: coarse x fine, X: fine x fine, P: fine x coarse) of the parts a, b (and c when with_c),
// accumulated in a fixed order (R row entries, then X row entries, then P row entries, all ascending).
inline Csr3 rap(const Csr& R, const Csr3& X, const Csr& P, bool with_c) {
    Csr3 out; out.n = R.n; out.rp.assign(1, 0);
    const int nc = P.m;
    Vec acc_a(nc, 0.), acc_b(nc, 0.), acc_c(nc, 0.); std::vector<int> mark(nc, -1), touched;
    for (int I = 0; I < R.n; ++I) {
        touched.clear();
        for (int k = R.rp[I]; k < R.rp[I + 1]; ++k) {
            const int i = R.ci[k]; const double w1 = R.v[k];
            for (int kk = X.rp[i]; kk < X.rp[i + 1]; ++kk) {
                const int j = X.ci[kk];
                const double xa = w1 * X.a[kk], xb = w1 * X.b[kk], xc = with_c ? w1 * X.c[kk] : 0.;
                for (int q = P.rp[j]; q < P.rp[j + 1]; ++q) {
                    const int J = P.ci[q]; const double w2 = P.v[q];
                    if (mark[J] != I) { mark[J] = I; touched.push_back(J); acc_a[J] = 0; acc_b[J] = 0; acc_c[J] = 0; }
                    acc_a[J] = acc_a[J] + xa * w2; acc_b[J] = acc_b[J] + xb * w2;
                    if (with_c) acc_c[J] = acc_c[J] + xc * w2;
                }
            }
        }
        std::sort(touched.begin(), touched.end());
        for (int J : touched) { out.ci.push_back(J); out.a.push_back(acc_a[J]); out.b.push_back(acc_b[J]); out.c.push_back(acc_c[J]); }
        out.rp.push_back(int(out.ci.size()));
    }
    return out;
}

// Fine-level operator parts from the Hall mesh (free rows only couple free columns; fixed rows are identity in part a).
inline Csr3 fine_operator(const axial::Mesh& g, const Vec& S) {
    need(S.size() == std::size_t(g.ne), "conduction size");
    // edge -> faces (ascending face index) and position in the face row
    std::vector<int> cp(g.ne + 1, 0), cf, ck;
    for (int b = 0; b < g.nb; ++b) for (int q = 0; q < 4; ++q) ++cp[g.C[b].index[q] + 1];
    for (int e = 0; e < g.ne; ++e) cp[e + 1] += cp[e];
    cf.resize(cp[g.ne]); ck.resize(cp[g.ne]);
    { std::vector<int> fill(cp.begin(), cp.end() - 1);
      for (int b = 0; b < g.nb; ++b) for (int q = 0; q < 4; ++q) { const int e = g.C[b].index[q]; cf[fill[e]] = b; ck[fill[e]] = q; ++fill[e]; } }
    Csr3 K; K.n = g.ne; K.rp.assign(1, 0);
    std::vector<std::pair<int, double>> row;
    for (int e = 0; e < g.ne; ++e) {
        if (g.fixed_E[e]) { K.ci.push_back(e); K.a.push_back(1.); K.b.push_back(0.); K.c.push_back(0.); K.rp.push_back(int(K.ci.size())); continue; }
        row.clear(); row.emplace_back(e, 0.);
        for (int t = cp[e]; t < cp[e + 1]; ++t) {
            const int b = cf[t]; const double ce = g.C[b].coefficient[ck[t]];
            for (int q = 0; q < 4; ++q) {
                const int e2 = g.C[b].index[q]; if (g.fixed_E[e2]) continue;
                const double v = g.H[b] * ce * g.C[b].coefficient[q];
                auto it = std::find_if(row.begin(), row.end(), [&](const auto& p) { return p.first == e2; });
                if (it == row.end()) row.emplace_back(e2, v); else it->second = it->second + v;
            }
        }
        std::sort(row.begin(), row.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
        for (const auto& [c, v] : row) { K.ci.push_back(c); K.a.push_back(c == e ? g.M[e] : 0.); K.b.push_back(c == e ? S[e] + g.Y[e] : 0.); K.c.push_back(v); }
        K.rp.push_back(int(K.ci.size()));
    }
    return K;
}
// Discrete gradient (edges x nodes): zero rows at fixed edges, zero columns at fixed nodes (phi = 0 on the PEC nodes,
// so G phi is an exact gradient field with zero tangential component on the PEC sheet).
inline Csr gradient(const Grid& G, const Bytes& fixed, const Bytes& fixed_nodes) {
    Csr D; D.n = G.ne; D.m = G.nn; D.rp.assign(1, 0);
    auto add = [&](int e, int a, int b, double h) {
        if (!fixed[e]) { if (!fixed_nodes[a]) { D.ci.push_back(a); D.v.push_back(-1. / h); } if (!fixed_nodes[b]) { D.ci.push_back(b); D.v.push_back(1. / h); } }
        D.rp.push_back(int(D.ci.size())); };
    const int n = G.n;
    for (int k = 0; k <= G.nz; ++k) for (int j = 0; j <= n; ++j) for (int i = 0; i < n; ++i) add(G.ex(i, j, k), G.node(i, j, k), G.node(i + 1, j, k), G.dx);
    for (int k = 0; k <= G.nz; ++k) for (int j = 0; j < n; ++j) for (int i = 0; i <= n; ++i) add(G.ey(i, j, k), G.node(i, j, k), G.node(i, j + 1, k), G.dx);
    for (int k = 0; k < G.nz; ++k) for (int j = 0; j <= n; ++j) for (int i = 0; i <= n; ++i) add(G.ez(i, j, k), G.node(i, j, k), G.node(i, j, k + 1), G.dz);
    return D;   // rows are in edge index order because the loops follow the edge numbering
}
// Edge prolongation from the coarse grid (n/r) to the fine grid (n); zero rows at fixed fine edges.
inline Csr prolongation(const Grid& F, const Grid& Cg, int r, const Bytes& fixed) {
    Csr P; P.n = F.ne; P.m = Cg.ne; P.rp.assign(1, 0);
    std::vector<std::pair<int, double>> row;
    auto emit = [&](int e) {
        if (!fixed[e]) { std::sort(row.begin(), row.end()); for (auto& [c, w] : row) if (w != 0) { P.ci.push_back(c); P.v.push_back(w); } }
        P.rp.push_back(int(P.ci.size())); row.clear(); };
    const int n = F.n; const double rr = r;
    for (int k = 0; k <= F.nz; ++k) for (int j = 0; j <= n; ++j) for (int i = 0; i < n; ++i) {
        const int ic = i / r, jc = j / r, m = j % r;
        row.emplace_back(Cg.ex(ic, jc, k), (rr - m) / rr); if (m) row.emplace_back(Cg.ex(ic, jc + 1, k), m / rr);
        emit(F.ex(i, j, k)); }
    for (int k = 0; k <= F.nz; ++k) for (int j = 0; j < n; ++j) for (int i = 0; i <= n; ++i) {
        const int ic = i / r, jc = j / r, m = i % r;
        row.emplace_back(Cg.ey(ic, jc, k), (rr - m) / rr); if (m) row.emplace_back(Cg.ey(ic + 1, jc, k), m / rr);
        emit(F.ey(i, j, k)); }
    for (int k = 0; k < F.nz; ++k) for (int j = 0; j <= n; ++j) for (int i = 0; i <= n; ++i) {
        const int ic = i / r, jc = j / r, mi = i % r, mj = j % r;
        for (int b = 0; b <= (mj ? 1 : 0); ++b) for (int a = 0; a <= (mi ? 1 : 0); ++a) {
            const double wi = a ? mi / rr : (rr - mi) / rr, wj = b ? mj / rr : (rr - mj) / rr;
            row.emplace_back(Cg.ez(ic + a, jc + b, k), wi * wj); }
        emit(F.ez(i, j, k)); }
    return P;
}
inline Csr combine(const Csr3& X, double alpha, double beta, bool with_c) {
    Csr K; K.n = X.n; K.m = X.n; K.rp = X.rp; K.ci = X.ci; K.v.resize(X.nnz());
    for (std::size_t k = 0; k < X.nnz(); ++k) K.v[k] = X.a[k] + alpha * X.b[k] + (with_c ? beta * X.c[k] : 0.);
    return K;
}
inline Vec l1_inverse(const Csr& K) {
    Vec d(K.n, 0.);
    for (int i = 0; i < K.n; ++i) { double s = 0; for (int k = K.rp[i]; k < K.rp[i + 1]; ++k) s += std::abs(K.v[k]); d[i] = s > 0 ? 1. / s : 0.; }
    return d;
}

struct Options { int nu = 1; int nu_coarse = 4; int min_cells = 5; int max_levels = 8; };
struct Hierarchy {
    std::vector<Level> L; std::vector<int> ratio;   // ratio[l]: coarsening factor from level l to l+1
    Vec sqrtM, inv_sqrtM; Bytes fixed0; double dt = 0; Options opt;
    std::size_t nnz_total() const { std::size_t s = 0; for (const auto& l : L) s += l.K.nnz() + l.N.nnz() + l.P.nnz() + l.PT.nnz() + l.G.nnz() + l.GT.nnz(); return s; }
};
inline Hierarchy build(const axial::Mesh& g, const Vec& S, Options opt = {}) {
    need(g.nx == g.ny && g.dx == g.dy, "square transverse grid required");
    Hierarchy H; H.opt = opt;
    H.sqrtM.resize(g.ne); H.inv_sqrtM.resize(g.ne); H.fixed0.resize(g.ne);
    for (int e = 0; e < g.ne; ++e) { H.sqrtM[e] = std::sqrt(g.M[e]); H.fixed0[e] = g.fixed_E[e] ? 1 : 0; H.inv_sqrtM[e] = g.fixed_E[e] ? 0. : 1. / H.sqrtM[e]; }
    Level l0; l0.g = Grid(g.nx, g.nz, g.dx, g.dz); l0.fixed = H.fixed0; l0.K3 = fine_operator(g, S);
    l0.fixed_nodes.resize(g.nn); for (int q = 0; q < g.nn; ++q) l0.fixed_nodes[q] = g.short_node[q] ? 1 : 0;
    H.L.push_back(std::move(l0));
    while (int(H.L.size()) < opt.max_levels) {
        Level& F = H.L.back(); const int n = F.g.n;
        int r = 0; if (n % 2 == 0 && n / 2 >= opt.min_cells) r = 2; else if (n % 3 == 0 && n / 3 >= opt.min_cells) r = 3;
        if (!r) break;
        Level C; C.g = Grid(n / r, F.g.nz, F.g.dx * r, F.g.dz);
        F.P = prolongation(F.g, C.g, r, F.fixed); F.PT = transpose(F.P);
        C.K3 = rap(F.PT, F.K3, F.P, true);
        C.fixed.assign(C.g.ne, 0);
        for (int i = 0; i < C.K3.n; ++i) {   // a coarse edge that no free fine edge interpolates from is fixed (identity row)
            bool any = false; for (int k = C.K3.rp[i]; k < C.K3.rp[i + 1]; ++k) if (C.K3.ci[k] == i && C.K3.a[k] != 0) any = true;
            C.fixed[i] = !any;
        }
        {   // identity rows for fixed coarse edges (their P columns are zero, so the rows/columns are empty)
            Csr3 X; X.n = C.K3.n; X.rp.assign(1, 0);
            for (int i = 0; i < C.K3.n; ++i) {
                if (C.fixed[i]) { X.ci.push_back(i); X.a.push_back(1.); X.b.push_back(0.); X.c.push_back(0.); }
                else for (int k = C.K3.rp[i]; k < C.K3.rp[i + 1]; ++k) { X.ci.push_back(C.K3.ci[k]); X.a.push_back(C.K3.a[k]); X.b.push_back(C.K3.b[k]); X.c.push_back(C.K3.c[k]); }
                X.rp.push_back(int(X.ci.size()));
            }
            C.K3 = std::move(X);
        }
        H.ratio.push_back(r);
        H.L.push_back(std::move(C));
    }
    for (std::size_t li = 0; li < H.L.size(); ++li) {
        auto& l = H.L[li];
        if (li > 0) {   // coarse levels: nodes whose incident edges are all fixed carry phi = 0
            l.fixed_nodes.assign(l.g.nn, 1);
            Csr full = gradient(l.g, Bytes(l.g.ne, 0), Bytes(l.g.nn, 0));
            for (int e = 0; e < l.g.ne; ++e) if (!l.fixed[e]) for (int k = full.rp[e]; k < full.rp[e + 1]; ++k) l.fixed_nodes[full.ci[k]] = 0;
        }
        l.G = gradient(l.g, l.fixed, l.fixed_nodes); l.GT = transpose(l.G);
        Csr3 X = l.K3;   // fixed rows carry the identity in part a; G has zero rows there, exclude them explicitly as well
        for (int i = 0; i < X.n; ++i) if (l.fixed[i]) for (int k = X.rp[i]; k < X.rp[i + 1]; ++k) { X.a[k] = 0; X.b[k] = 0; X.c[k] = 0; }
        l.N3 = rap(l.GT, X, l.G, true);
    }
    return H;
}
inline void set_dt(Hierarchy& H, double dt) {
    H.dt = dt;
    for (auto& l : H.L) {
        l.K = combine(l.K3, .5 * dt, .25 * dt * dt, true); l.invd = l1_inverse(l.K);
        l.N = combine(l.N3, .5 * dt, .25 * dt * dt, true); l.invdn = l1_inverse(l.N);
    }
}
// Hall's Jacobi diagonal of A (Model::prepare), for the reference Jacobi PCG.
inline Vec jacobi_diagonal(const axial::Mesh& g, const Vec& S, double dt) {
    Vec d(g.ne, 1.);
    for (int e = 0; e < g.ne; ++e) if (!g.fixed_E[e]) d[e] += .5 * dt * (S[e] + g.Y[e]) / g.M[e];
    for (int b = 0; b < g.nb; ++b) for (int k = 0; k < 4; ++k) { int e = g.C[b].index[k]; if (!g.fixed_E[e]) d[e] += .25 * dt * dt * g.H[b] * std::pow(g.C[b].coefficient[k], 2) / g.M[e]; }
    return d;
}

// ------------------------------------------------------------------ CPU reference (same operation order as the GPU)
constexpr int TPB = 256, NBLK = 512;
inline double tree_dot(const Vec& a, const Vec& b) {
    std::vector<double> part(NBLK), sh(TPB);
    const int n = int(a.size());
    for (int blk = 0; blk < NBLK; ++blk) {
        for (int t = 0; t < TPB; ++t) { double s = 0.0; for (int i = blk * TPB + t; i < n; i += NBLK * TPB) s = s + a[i] * b[i]; sh[t] = s; }
        for (int s = TPB / 2; s > 0; s >>= 1) for (int t = 0; t < s; ++t) sh[t] = sh[t] + sh[t + s];
        part[blk] = sh[0];
    }
    for (int s = NBLK / 2; s > 0; s >>= 1) for (int t = 0; t < s; ++t) part[t] = part[t] + part[t + s];
    return part[0];
}
inline void residual(const Csr& K, const Vec& x, const Vec& b, Vec& r) { for (int i = 0; i < K.n; ++i) r[i] = gmg_row_residual(i, K.rp.data(), K.ci.data(), K.v.data(), x.data(), b.data()); }
inline void mul(const Csr& K, const Vec& x, Vec& y) { for (int i = 0; i < K.n; ++i) y[i] = gmg_row_mul(i, K.rp.data(), K.ci.data(), K.v.data(), x.data()); }
inline void mul_add(const Csr& K, const Vec& x, Vec& y) { for (int i = 0; i < K.n; ++i) y[i] = y[i] + gmg_row_mul(i, K.rp.data(), K.ci.data(), K.v.data(), x.data()); }
inline void jacobi(const Vec& d, const Vec& r, Vec& x) { for (std::size_t i = 0; i < x.size(); ++i) x[i] = x[i] + d[i] * r[i]; }
inline void scale(const Vec& d, const Vec& x, Vec& y) { for (std::size_t i = 0; i < y.size(); ++i) y[i] = d[i] * x[i]; }

struct CpuWork { std::vector<Vec> x, b, r, rn, phi; Vec t0, t1; };
inline CpuWork cpu_work(const Hierarchy& H) {
    CpuWork w; for (const auto& l : H.L) { w.x.emplace_back(l.g.ne); w.b.emplace_back(l.g.ne); w.r.emplace_back(l.g.ne); w.rn.emplace_back(l.g.nn); w.phi.emplace_back(l.g.nn); }
    w.t0.resize(H.L[0].g.ne); w.t1.resize(H.L[0].g.ne); return w;
}
inline void cpu_edge(const Level& l, CpuWork& w, int i) { residual(l.K, w.x[i], w.b[i], w.r[i]); jacobi(l.invd, w.r[i], w.x[i]); }
inline void cpu_node(const Level& l, CpuWork& w, int i) { residual(l.K, w.x[i], w.b[i], w.r[i]); mul(l.GT, w.r[i], w.rn[i]); scale(l.invdn, w.rn[i], w.phi[i]); mul_add(l.G, w.phi[i], w.x[i]); }
inline void cpu_vcycle(const Hierarchy& H, CpuWork& w, int i) {
    const Level& l = H.L[i];
    std::fill(w.x[i].begin(), w.x[i].end(), 0.);
    if (i + 1 == int(H.L.size())) {
        for (int s = 0; s < H.opt.nu_coarse; ++s) { cpu_edge(l, w, i); cpu_node(l, w, i); }
        for (int s = 0; s < H.opt.nu_coarse; ++s) { cpu_node(l, w, i); cpu_edge(l, w, i); }
        return;
    }
    for (int s = 0; s < H.opt.nu; ++s) { cpu_edge(l, w, i); cpu_node(l, w, i); }
    residual(l.K, w.x[i], w.b[i], w.r[i]); mul(l.PT, w.r[i], w.b[i + 1]);
    cpu_vcycle(H, w, i + 1);
    mul_add(l.P, w.x[i + 1], w.x[i]);
    for (int s = 0; s < H.opt.nu; ++s) { cpu_node(l, w, i); cpu_edge(l, w, i); }
}
// z = B r (preconditioner of A): Jacobi (jac non-empty) or the V-cycle
inline void cpu_prec(const Hierarchy& H, CpuWork& w, const Vec& jac_inv, const Vec& r, Vec& z) {
    if (!jac_inv.empty()) { scale(jac_inv, r, z); return; }
    const int n = H.L[0].g.ne;
    for (int e = 0; e < n; ++e) w.b[0][e] = H.fixed0[e] ? 0.0 : H.sqrtM[e] * r[e];
    cpu_vcycle(H, w, 0);
    for (int e = 0; e < n; ++e) z[e] = H.fixed0[e] ? r[e] : H.sqrtM[e] * w.x[0][e];
}
// q = A p = fixed ? p : (K (p / sqrtM)) / sqrtM
inline void cpu_apply(const Hierarchy& H, CpuWork& w, const Vec& p, Vec& q) {
    const int n = H.L[0].g.ne;
    scale(H.inv_sqrtM, p, w.t0); mul(H.L[0].K, w.t0, w.t1);
    for (int e = 0; e < n; ++e) q[e] = H.fixed0[e] ? p[e] : H.inv_sqrtM[e] * w.t1[e];
}
struct SolveStats { int iterations = 0, corrections = 0; double rhs_norm = 0, target = 0, backend_residual = 0; std::vector<int> cycle_iterations; };
// CG with the Hall gpu1 stopping rule; the same scalar sequence as GpuSolver::solve.
inline Vec cpu_solve(const Hierarchy& H, const Vec& jac_inv, const Vec& b, int maxit, SolveStats& st) {
    const int n = int(b.size()); CpuWork w = cpu_work(H);
    Vec x(n, 0.), r = b, z(n), p(n), q(n);
    st.rhs_norm = std::sqrt(tree_dot(b, b)); st.target = std::max(1e-30, 2e-13 * st.rhs_norm);
    for (int corr = 0; corr < 3; ++corr) {
        st.corrections = corr + 1;
        cpu_prec(H, w, jac_inv, r, z); p = z; double rz = tree_dot(r, z); int it = 0;
        double rr = tree_dot(r, r);
        while (it < maxit && std::sqrt(rr) > .2 * st.target) {
            cpu_apply(H, w, p, q); const double pq = tree_dot(p, q); const double alpha = rz / pq;
            for (int i = 0; i < n; ++i) { x[i] = x[i] + alpha * p[i]; r[i] = r[i] - alpha * q[i]; }
            rr = tree_dot(r, r); ++it;
            if (!(std::sqrt(rr) > .2 * st.target)) break;
            cpu_prec(H, w, jac_inv, r, z); const double rz2 = tree_dot(r, z); const double beta = rz2 / rz; rz = rz2;
            for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        }
        st.iterations += it; st.cycle_iterations.push_back(it);
        cpu_apply(H, w, x, q); for (int i = 0; i < n; ++i) r[i] = b[i] - q[i];
        st.backend_residual = std::sqrt(tree_dot(r, r));
        if (st.backend_residual <= st.target) return x;
    }
    return x;
}
}  // namespace gmg
