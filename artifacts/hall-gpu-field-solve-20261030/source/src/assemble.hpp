// Explicit assembly of the Hall implicit edge operator from an axial::Mesh,
// a step dt and a cell conductivity field. Every scalar that also appears in
// Model::prepare (field_model.hpp:79-85) is evaluated with the identical
// expression text, so the preconditioner diagonal, a = (.5*dt)*(S+Y),
// q = (.25*dt)*dt and sqrt(M) are bitwise equal to the CPU path. This is
// checked at run time against the diagonal handed to the solver.
#pragma once
#include "csr_operator.hpp"

namespace hgf {
template <class Mesh>
Operator assemble(const Mesh& g, double dt, const Vec& sigma) {
    Operator A;
    A.n = g.ne; A.nb = g.nb; A.dt = dt;
    A.q = .25 * dt * dt;
    const auto S = g.conduction(sigma);
    A.sqrtM.resize(g.ne); A.a.resize(g.ne); A.fixed.resize(g.ne);
    Vec diagonal(g.ne, 1);
    for (int e = 0; e < g.ne; ++e) {
        A.sqrtM[e] = std::sqrt(g.M[e]);
        A.fixed[e] = g.fixed_E[e] ? 1 : 0;
        A.a[e] = .5 * dt * (S[e] + g.Y[e]);
        if (!g.fixed_E[e]) diagonal[e] += .5 * dt * (S[e] + g.Y[e]) / g.M[e];
    }
    for (int b = 0; b < g.nb; ++b)
        for (int k = 0; k < 4; ++k) {
            int e = g.C[b].index[k];
            if (!g.fixed_E[e]) diagonal[e] += .25 * dt * dt * g.H[b] * std::pow(g.C[b].coefficient[k], 2) / g.M[e];
        }
    A.diag = diagonal;
    // factored rows
    A.face_edge.resize(4 * std::size_t(g.nb)); A.face_coef.resize(4 * std::size_t(g.nb)); A.H.resize(g.nb);
    std::vector<std::int32_t> count(g.ne + 1, 0);
    for (int b = 0; b < g.nb; ++b) {
        need(g.C[b].count == 4, "axial curl rows have four entries");
        A.H[b] = g.H[b];
        for (int k = 0; k < 4; ++k) {
            A.face_edge[4 * b + k] = g.C[b].index[k];
            A.face_coef[4 * b + k] = g.C[b].coefficient[k];
            ++count[g.C[b].index[k] + 1];
        }
    }
    A.ct_ptr.assign(g.ne + 1, 0);
    for (int e = 0; e < g.ne; ++e) A.ct_ptr[e + 1] = A.ct_ptr[e] + count[e + 1];
    A.ct_face.resize(A.ct_ptr[g.ne]); A.ct_coef.resize(A.ct_ptr[g.ne]);
    std::vector<std::int32_t> fill(A.ct_ptr.begin(), A.ct_ptr.end() - 1);
    for (int b = 0; b < g.nb; ++b)  // ascending face order per edge, as in curl_transpose
        for (int k = 0; k < 4; ++k) {
            const int e = g.C[b].index[k];
            A.ct_face[fill[e]] = b; A.ct_coef[fill[e]] = g.C[b].coefficient[k]; ++fill[e];
        }
    // assembled rows
    A.row_ptr.assign(g.ne + 1, 0);
    A.col.reserve(std::size_t(g.ne) * 13); A.val.reserve(std::size_t(g.ne) * 13);
    std::vector<std::pair<std::int32_t, double>> row;
    for (int e = 0; e < g.ne; ++e) {
        row.clear();
        row.emplace_back(e, A.diag[e]);
        if (!g.fixed_E[e]) {
            for (std::int32_t j = A.ct_ptr[e]; j < A.ct_ptr[e + 1]; ++j) {
                const std::int32_t b = A.ct_face[j];
                const double ce = A.ct_coef[j];
                for (int k = 0; k < 4; ++k) {
                    const std::int32_t f = A.face_edge[4 * b + k];
                    if (f == e || g.fixed_E[f]) continue;
                    const double cf = A.face_coef[4 * b + k];
                    const double v = (A.q * A.H[b]) * (ce * cf) / (A.sqrtM[e] * A.sqrtM[f]);
                    auto it = std::find_if(row.begin(), row.end(), [&](const auto& p) { return p.first == f; });
                    need(it == row.end(), "two faces share an edge pair");
                    row.emplace_back(f, v);
                }
            }
        }
        std::sort(row.begin(), row.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
        for (const auto& [c, v] : row) { A.col.push_back(c); A.val.push_back(v); }
        A.row_ptr[e + 1] = std::int32_t(A.col.size());
        if (S[e] != 0) {
            A.sigma_rows.push_back(e);
            for (std::int32_t j = A.row_ptr[e]; j < A.row_ptr[e + 1]; ++j) if (A.col[j] == e) A.sigma_diag_pos.push_back(j);
        }
    }
    need(A.sigma_rows.size() == A.sigma_diag_pos.size(), "sigma row diagonal positions");
    return A;
}

// Assembled-matrix symmetry check (the entry formula is symmetric in e and f).
inline bool csr_symmetric(const Operator& A) {
    for (std::int32_t i = 0; i < A.n; ++i)
        for (std::int32_t j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
            const std::int32_t c = A.col[j];
            const auto begin = A.col.begin() + A.row_ptr[c], end = A.col.begin() + A.row_ptr[c + 1];
            const auto it = std::lower_bound(begin, end, i);
            const bool present = it != end && *it == i;
            if (A.fixed[i] || A.fixed[c]) { if (c != i && present) return false; continue; }
            if (!present || A.val[std::size_t(it - A.col.begin())] != A.val[j]) return false;
        }
    return true;
}
}  // namespace hgf
