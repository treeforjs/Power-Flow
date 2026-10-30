// Hall implicit edge operator in explicit form, host reference arithmetic and
// fixture I/O. No Hall headers are included here, so the GPU executable can use
// this file without MFEM. Assembly from axial::Mesh lives in assemble.hpp.
//
// Unknown: v = sqrt(M) * dE on every edge (the scaled variable of Model::prepare).
// Free rows:  (A v)_e = v_e + ( a_e * E_e + q * K_e ) / sqrtM_e,
//             E_f = v_f / sqrtM_f (0 on PEC edges), K = C^T H C E,
//             a_e = (.5*dt)*(S_e+Y_e), q = (.25*dt)*dt.
// Fixed (PEC) rows: (A v)_e = v_e.
// Two device representations are kept:
//  * assembled: one scaled CSR matrix A (13 entries per interior row);
//  * factored:  curl rows (4 per face) and transposed curl rows (faces in
//    ascending order per edge) that reproduce the matrix-free apply exactly.
#pragma once
#include <openssl/evp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace hgf {
using Vec = std::vector<double>;

inline void need(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}

// ---------------------------------------------------------------- sha256
inline std::string hex(const unsigned char* d, std::size_t n) {
    std::ostringstream o;
    for (std::size_t i = 0; i < n; ++i) o << std::hex << std::setw(2) << std::setfill('0') << unsigned(d[i]);
    return o.str();
}
inline std::string sha256_bytes(const void* data, std::size_t size) {
    unsigned char out[EVP_MAX_MD_SIZE];
    unsigned int len = 0;
    need(EVP_Digest(data, size, out, &len, EVP_sha256(), nullptr) == 1, "sha256 digest failed");
    return hex(out, len);
}
template <class T>
inline std::string sha256_vec(const std::vector<T>& v) { return sha256_bytes(v.data(), v.size() * sizeof(T)); }
inline std::string sha256_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    need(bool(in), "cannot read " + path);
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    need(ctx && EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr) == 1, "sha256 init failed");
    std::vector<char> buf(1 << 20);
    while (in) {
        in.read(buf.data(), std::streamsize(buf.size()));
        if (in.gcount() > 0) EVP_DigestUpdate(ctx, buf.data(), std::size_t(in.gcount()));
    }
    unsigned char out[EVP_MAX_MD_SIZE];
    unsigned int len = 0;
    EVP_DigestFinal_ex(ctx, out, &len);
    EVP_MD_CTX_free(ctx);
    return hex(out, len);
}

// ---------------------------------------------------------------- operator
struct Operator {
    std::int32_t n = 0, nb = 0;
    double dt = 0, q = 0;
    // assembled scaled CSR
    std::vector<std::int32_t> row_ptr, col;
    Vec val;
    // factored representation
    std::vector<std::int32_t> face_edge;  // 4 per face
    Vec face_coef, H;                     // 4 per face, 1 per face
    std::vector<std::int32_t> ct_ptr, ct_face;
    Vec ct_coef;
    Vec sqrtM, a;
    std::vector<std::uint8_t> fixed;
    Vec diag;  // Jacobi preconditioner, bitwise equal to Model::prepare's diagonal
    // Rows whose entries depend on sigma (S_e > 0) and the position of their
    // diagonal entry in the assembled CSR: the only values a sigma(T) update changes.
    std::vector<std::int32_t> sigma_rows, sigma_diag_pos;
    std::size_t nnz() const { return col.size(); }
};

// Assembled CSR product. Each row sums from +0.0 in ascending column order;
// the device kernel uses the same order.
inline void spmv_csr(const Operator& A, const double* v, double* out) {
    for (std::int32_t i = 0; i < A.n; ++i) {
        double y = 0.0;
        for (std::int32_t j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) y += A.val[j] * v[A.col[j]];
        out[i] = y;
    }
}

// Factored product in the exact operation order of the matrix-free apply in
// Model::prepare (field_model.hpp:84-85) and axial::Mesh::curl /
// curl_adjoint / curl_transpose (axial_mesh.hpp:84-87).
inline void apply_factored(const Operator& A, const double* v, double* out, Vec& Bh) {
    Bh.assign(std::size_t(A.nb), 0.0);
    for (std::int32_t b = 0; b < A.nb; ++b) {
        double s = 0.0;
        for (int k = 0; k < 4; ++k) {
            const std::int32_t e = A.face_edge[4 * b + k];
            const double E = A.fixed[e] ? 0.0 : v[e] / A.sqrtM[e];
            s += A.face_coef[4 * b + k] * E;
        }
        Bh[b] = s * A.H[b];
    }
    for (std::int32_t e = 0; e < A.n; ++e) {
        if (A.fixed[e]) { out[e] = v[e]; continue; }
        const double E = v[e] / A.sqrtM[e];
        double K = 0.0;
        for (std::int32_t j = A.ct_ptr[e]; j < A.ct_ptr[e + 1]; ++j) K += A.ct_coef[j] * Bh[A.ct_face[j]];
        out[e] = v[e] + (A.a[e] * E + A.q * K) / A.sqrtM[e];
    }
}

// ---------------------------------------------------------------- fixed-tree reductions
// Device layout: NBLK blocks of TPB threads, grid-stride rows (thread slot
// gid = block*TPB + t takes rows gid, gid+NBLK*TPB, ...), shared-memory halving
// tree inside each block, then one block of NBLK threads halving the partials.
// The constants do not depend on the device, so the summation tree is the
// same on every GPU architecture.
constexpr int TPB = 256;
constexpr int NBLK = 512;
constexpr int SLOTS = TPB * NBLK;

inline double tree_halve(double* sh, int count) {
    for (int s = count / 2; s > 0; s >>= 1)
        for (int t = 0; t < s; ++t) sh[t] = sh[t] + sh[t + s];
    return sh[0];
}
// Emulates block trees and the final tree given per-slot accumulators.
inline double tree_finish(std::vector<double>& slot) {
    std::array<double, NBLK> partial{};
    std::array<double, TPB> sh{};
    for (int b = 0; b < NBLK; ++b) {
        for (int t = 0; t < TPB; ++t) sh[t] = slot[std::size_t(b) * TPB + t];
        partial[b] = tree_halve(sh.data(), TPB);
    }
    return tree_halve(partial.data(), NBLK);
}

// ---------------------------------------------------------------- CG trace
struct CgTrace {
    double bnorm = 0, target = 0, true_norm = 0;
    int corrections = 0;
    std::vector<int> iterations;               // per correction
    std::vector<std::vector<double>> alpha, beta, rnorm;  // per correction
    bool converged = false;
    int total_iterations() const { int s = 0; for (int i : iterations) s += i; return s; }
};

// Host emulation of the device CG: same kernels, same fixed-tree sums, same
// gating. Its result is expected to equal the device result bit for bit.
template <class Apply>
Vec emulate_device_cg(Apply&& apply, const Vec& b, const Vec& diag, int maxit, CgTrace& tr) {
    const std::size_t n = b.size();
    std::vector<double> s1(SLOTS), s2(SLOTS);
    auto reset = [&] { std::fill(s1.begin(), s1.end(), 0.0); std::fill(s2.begin(), s2.end(), 0.0); };
    reset();
    for (std::size_t i = 0; i < n; ++i) s1[i % SLOTS] += b[i] * b[i];
    tr.bnorm = std::sqrt(tree_finish(s1));
    tr.target = std::fmax(1e-30, 2e-13 * tr.bnorm);
    Vec x(n, 0.0), r = b, z(n), p(n), Ap(n);
    for (int corr = 0; corr < 3; ++corr) {
        tr.corrections = corr + 1;
        tr.alpha.emplace_back(); tr.beta.emplace_back(); tr.rnorm.emplace_back();
        reset();
        for (std::size_t i = 0; i < n; ++i) {
            z[i] = r[i] / diag[i]; p[i] = z[i];
            s1[i % SLOTS] += r[i] * z[i]; s2[i % SLOTS] += r[i] * r[i];
        }
        double rz = tree_finish(s1), rr = tree_finish(s2), rnorm = std::sqrt(rr);
        int it = 0;
        bool done = !(it < maxit && rnorm > tr.target * .2);
        tr.rnorm.back().push_back(rnorm);
        while (!done) {
            apply(p.data(), Ap.data());
            reset();
            for (std::size_t i = 0; i < n; ++i) s1[i % SLOTS] += p[i] * Ap[i];
            const double pAp = tree_finish(s1);
            need(pAp > 0 && std::isfinite(pAp), "non-SPD field operator");
            const double alpha = rz / pAp;
            reset();
            for (std::size_t i = 0; i < n; ++i) {
                x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; z[i] = r[i] / diag[i];
                s1[i % SLOTS] += r[i] * z[i]; s2[i % SLOTS] += r[i] * r[i];
            }
            const double nrz = tree_finish(s1);
            rr = tree_finish(s2);
            const double beta = nrz / rz;
            rz = nrz; rnorm = std::sqrt(rr); ++it;
            tr.alpha.back().push_back(alpha); tr.beta.back().push_back(beta); tr.rnorm.back().push_back(rnorm);
            done = !(it < maxit && rnorm > tr.target * .2);
            if (!done) for (std::size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        }
        tr.iterations.push_back(it);
        apply(x.data(), Ap.data());
        reset();
        for (std::size_t i = 0; i < n; ++i) { r[i] = b[i] - Ap[i]; s1[i % SLOTS] += r[i] * r[i]; }
        tr.true_norm = std::sqrt(tree_finish(s1));
        if (tr.true_norm <= tr.target) { tr.converged = true; return x; }
    }
    throw std::runtime_error("field solve independent residual failed");
}

// Plain single-thread CPU Jacobi PCG on a preassembled operator (sequential
// double dot products), same stopping rule. Used only as a timing baseline for
// an optimized CPU path; it is not a reference for bitwise comparisons.
template <class Apply>
Vec plain_cpu_pcg(Apply&& apply, const Vec& b, const Vec& diag, int maxit, CgTrace& tr) {
    const std::size_t n = b.size();
    auto dot = [&](const Vec& u, const Vec& w) { double s = 0; for (std::size_t i = 0; i < n; ++i) s += u[i] * w[i]; return s; };
    tr.bnorm = std::sqrt(dot(b, b));
    tr.target = std::max(1e-30, 2e-13 * tr.bnorm);
    Vec x(n, 0.0), r = b, z(n), p(n), Ap(n);
    for (int corr = 0; corr < 3; ++corr) {
        tr.corrections = corr + 1;
        double rz = 0, rr = 0;
        for (std::size_t i = 0; i < n; ++i) { z[i] = r[i] / diag[i]; p[i] = z[i]; rz += r[i] * z[i]; rr += r[i] * r[i]; }
        int it = 0;
        while (it < maxit && std::sqrt(rr) > tr.target * .2) {
            apply(p.data(), Ap.data());
            const double pAp = dot(p, Ap);
            need(pAp > 0 && std::isfinite(pAp), "non-SPD field operator");
            const double alpha = rz / pAp;
            double nrz = 0; rr = 0;
            for (std::size_t i = 0; i < n; ++i) {
                x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; z[i] = r[i] / diag[i]; nrz += r[i] * z[i]; rr += r[i] * r[i];
            }
            const double beta = nrz / rz; rz = nrz;
            for (std::size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
            ++it;
        }
        tr.iterations.push_back(it);
        apply(x.data(), Ap.data());
        double t = 0;
        for (std::size_t i = 0; i < n; ++i) { r[i] = b[i] - Ap[i]; t += r[i] * r[i]; }
        tr.true_norm = std::sqrt(t);
        if (tr.true_norm <= tr.target) { tr.converged = true; return x; }
    }
    throw std::runtime_error("field solve independent residual failed");
}

// ---------------------------------------------------------------- Lanczos estimate
// Eigenvalue bounds of the Lanczos tridiagonal built from CG coefficients of
// one correction cycle (Jacobi-preconditioned operator). Extreme Ritz values
// are inner estimates, so kappa = lmax/lmin is a lower estimate of the true
// condition number of D^{-1/2} A D^{-1/2}.
struct Ritz { double lmin = 0, lmax = 0, kappa = 0; int k = 0; };
inline Ritz lanczos_extremes(const std::vector<double>& alpha, const std::vector<double>& beta) {
    Ritz out;
    const int k = int(alpha.size());
    out.k = k;
    if (k == 0) return out;
    std::vector<double> d(k), e(std::max(0, k - 1));
    for (int j = 0; j < k; ++j) {
        d[j] = 1 / alpha[j] + (j > 0 ? beta[j - 1] / alpha[j - 1] : 0.0);
        if (j + 1 < k) e[j] = std::sqrt(beta[j]) / alpha[j];
    }
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (int j = 0; j < k; ++j) {
        double rad = (j > 0 ? std::abs(e[j - 1]) : 0) + (j + 1 < k ? std::abs(e[j]) : 0);
        lo = std::min(lo, d[j] - rad); hi = std::max(hi, d[j] + rad);
    }
    auto count_below = [&](double x) {
        int c = 0; double q = d[0] - x;
        if (q < 0) ++c;
        for (int j = 1; j < k; ++j) {
            if (q == 0) q = 1e-300;
            q = (d[j] - x) - e[j - 1] * e[j - 1] / q;
            if (q < 0) ++c;
        }
        return c;
    };
    auto bisect = [&](int index) {  // index-th smallest eigenvalue (0-based)
        double a = lo, b = hi;
        for (int it = 0; it < 200 && b - a > 1e-15 * std::max(std::abs(a), std::abs(b)); ++it) {
            double m = .5 * (a + b);
            if (count_below(m) > index) b = m; else a = m;
        }
        return .5 * (a + b);
    };
    out.lmin = bisect(0);
    out.lmax = bisect(k - 1);
    out.kappa = out.lmax / out.lmin;
    return out;
}

// ---------------------------------------------------------------- deterministic test vectors
inline std::uint64_t splitmix(std::uint64_t& s) {
    std::uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
inline Vec pseudo_random(std::size_t n, std::uint64_t seed) {
    Vec v(n);
    for (auto& x : v) x = double(splitmix(seed) >> 11) * 0x1.0p-53 * 2 - 1;
    return v;
}

struct Diff { double max_abs = 0, max_rel = 0, norm_rel = 0, max_scaled = 0; std::size_t differing = 0; bool bitwise = true; };
// |A| |v| with the assembled entries: the natural rounding scale of (A v)_e.
inline Vec abs_product(const Operator& A, const Vec& v) {
    Vec out(v.size());
    for (std::int32_t i = 0; i < A.n; ++i) {
        double y = 0.0;
        for (std::int32_t j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) y += std::abs(A.val[j]) * std::abs(v[A.col[j]]);
        out[i] = y;
    }
    return out;
}
inline Diff compare(const Vec& got, const Vec& ref, const Vec* scale = nullptr) {
    Diff d;
    long double num = 0, den = 0;
    for (std::size_t i = 0; i < ref.size(); ++i) {
        std::uint64_t a, b;
        std::memcpy(&a, &got[i], 8); std::memcpy(&b, &ref[i], 8);
        if (a != b) { d.bitwise = false; ++d.differing; }
        const double e = std::abs(got[i] - ref[i]);
        d.max_abs = std::max(d.max_abs, e);
        if (ref[i] != 0) d.max_rel = std::max(d.max_rel, e / std::abs(ref[i]));
        else if (e != 0) d.max_rel = std::numeric_limits<double>::infinity();
        if (scale && e != 0) d.max_scaled = std::max(d.max_scaled, (*scale)[i] > 0 ? e / (*scale)[i] : std::numeric_limits<double>::infinity());
        num += (long double)e * e; den += (long double)ref[i] * ref[i];
    }
    d.norm_rel = den > 0 ? double(std::sqrt(num / den)) : (num > 0 ? std::numeric_limits<double>::infinity() : 0.0);
    return d;
}

// ---------------------------------------------------------------- binary fixture I/O
// Record format: u32 name length, name, u8 type (1=f64, 2=i32, 3=u8),
// u64 element count, raw little-endian payload. A zero-length name ends the file.
class ArrayFile {
public:
    struct Entry { std::uint8_t type = 0; std::vector<char> bytes; };
    std::map<std::string, Entry> entries;
    static constexpr char magic[8] = {'H', 'G', 'F', 'A', 'R', 'R', '0', '1'};

    void put(const std::string& name, const Vec& v) { put_raw(name, 1, v.data(), v.size() * 8); }
    void put(const std::string& name, const std::vector<std::int32_t>& v) { put_raw(name, 2, v.data(), v.size() * 4); }
    void put(const std::string& name, const std::vector<std::uint8_t>& v) { put_raw(name, 3, v.data(), v.size()); }
    void put_scalar(const std::string& name, double x) { put(name, Vec{x}); }
    void put_raw(const std::string& name, std::uint8_t type, const void* p, std::size_t bytes) {
        Entry e; e.type = type; e.bytes.resize(bytes);
        if (bytes) std::memcpy(e.bytes.data(), p, bytes);
        entries[name] = std::move(e);
    }
    template <class T> std::vector<T> get(const std::string& name, std::uint8_t type) const {
        auto it = entries.find(name);
        need(it != entries.end(), "fixture array missing: " + name);
        need(it->second.type == type && it->second.bytes.size() % sizeof(T) == 0, "fixture array type mismatch: " + name);
        std::vector<T> out(it->second.bytes.size() / sizeof(T));
        if (!out.empty()) std::memcpy(out.data(), it->second.bytes.data(), it->second.bytes.size());
        return out;
    }
    Vec f64(const std::string& n) const { return get<double>(n, 1); }
    std::vector<std::int32_t> i32(const std::string& n) const { return get<std::int32_t>(n, 2); }
    std::vector<std::uint8_t> u8(const std::string& n) const { return get<std::uint8_t>(n, 3); }
    double scalar(const std::string& n) const { auto v = f64(n); need(v.size() == 1, "fixture scalar shape: " + n); return v[0]; }

    void write(const std::string& path) const {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        need(bool(out), "cannot write " + path);
        out.write(magic, 8);
        for (const auto& [name, e] : entries) {
            const std::uint32_t len = std::uint32_t(name.size());
            const std::uint64_t count = e.bytes.size() / (e.type == 1 ? 8 : e.type == 2 ? 4 : 1);
            out.write(reinterpret_cast<const char*>(&len), 4);
            out.write(name.data(), len);
            out.write(reinterpret_cast<const char*>(&e.type), 1);
            out.write(reinterpret_cast<const char*>(&count), 8);
            out.write(e.bytes.data(), std::streamsize(e.bytes.size()));
        }
        const std::uint32_t end = 0;
        out.write(reinterpret_cast<const char*>(&end), 4);
        need(bool(out), "write failed: " + path);
    }
    static ArrayFile read(const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        need(bool(in), "cannot read " + path);
        char m[8];
        in.read(m, 8);
        need(bool(in) && std::memcmp(m, magic, 8) == 0, "bad fixture magic: " + path);
        ArrayFile f;
        for (;;) {
            std::uint32_t len = 0;
            in.read(reinterpret_cast<char*>(&len), 4);
            need(bool(in), "truncated fixture: " + path);
            if (len == 0) break;
            need(len < 256, "fixture name length");
            std::string name(len, '\0');
            in.read(name.data(), len);
            Entry e;
            std::uint64_t count = 0;
            in.read(reinterpret_cast<char*>(&e.type), 1);
            in.read(reinterpret_cast<char*>(&count), 8);
            need(bool(in) && (e.type >= 1 && e.type <= 3), "fixture record header");
            const std::size_t size = count * (e.type == 1 ? 8 : e.type == 2 ? 4 : 1);
            e.bytes.resize(size);
            in.read(e.bytes.data(), std::streamsize(size));
            need(bool(in), "truncated fixture payload: " + name);
            f.entries[name] = std::move(e);
        }
        return f;
    }
};

inline void store_operator(ArrayFile& f, const Operator& A) {
    f.put("n", std::vector<std::int32_t>{A.n});
    f.put("nb", std::vector<std::int32_t>{A.nb});
    f.put_scalar("dt", A.dt);
    f.put_scalar("q", A.q);
    f.put("row_ptr", A.row_ptr); f.put("col", A.col); f.put("val", A.val);
    f.put("face_edge", A.face_edge); f.put("face_coef", A.face_coef); f.put("H", A.H);
    f.put("ct_ptr", A.ct_ptr); f.put("ct_face", A.ct_face); f.put("ct_coef", A.ct_coef);
    f.put("sqrtM", A.sqrtM); f.put("a", A.a); f.put("fixed", A.fixed); f.put("diag", A.diag);
    f.put("sigma_rows", A.sigma_rows); f.put("sigma_diag_pos", A.sigma_diag_pos);
}
inline Operator load_operator(const ArrayFile& f) {
    Operator A;
    A.n = f.i32("n").at(0); A.nb = f.i32("nb").at(0);
    A.dt = f.scalar("dt"); A.q = f.scalar("q");
    A.row_ptr = f.i32("row_ptr"); A.col = f.i32("col"); A.val = f.f64("val");
    A.face_edge = f.i32("face_edge"); A.face_coef = f.f64("face_coef"); A.H = f.f64("H");
    A.ct_ptr = f.i32("ct_ptr"); A.ct_face = f.i32("ct_face"); A.ct_coef = f.f64("ct_coef");
    A.sqrtM = f.f64("sqrtM"); A.a = f.f64("a"); A.fixed = f.u8("fixed"); A.diag = f.f64("diag");
    A.sigma_rows = f.i32("sigma_rows"); A.sigma_diag_pos = f.i32("sigma_diag_pos");
    const std::size_t n = std::size_t(A.n), nb = std::size_t(A.nb);
    need(A.row_ptr.size() == n + 1 && A.ct_ptr.size() == n + 1 && A.face_edge.size() == 4 * nb &&
         A.face_coef.size() == 4 * nb && A.H.size() == nb && A.sqrtM.size() == n && A.a.size() == n &&
         A.fixed.size() == n && A.diag.size() == n && A.val.size() == A.col.size() &&
         std::size_t(A.row_ptr[n]) == A.col.size() && A.ct_face.size() == A.ct_coef.size() &&
         std::size_t(A.ct_ptr[n]) == A.ct_face.size(), "operator fixture shape");
    for (auto c : A.col) need(c >= 0 && c < A.n, "operator column index");
    for (auto e : A.face_edge) need(e >= 0 && e < A.n, "operator face edge index");
    for (auto b : A.ct_face) need(b >= 0 && b < A.nb, "operator transposed face index");
    for (std::size_t i = 0; i < n; ++i) need(A.row_ptr[i] <= A.row_ptr[i + 1] && A.ct_ptr[i] <= A.ct_ptr[i + 1], "operator row pointers");
    need(A.sigma_rows.size() == A.sigma_diag_pos.size(), "sigma row list shape");
    for (std::size_t k = 0; k < A.sigma_rows.size(); ++k) {
        const auto e = A.sigma_rows[k], p = A.sigma_diag_pos[k];
        need(e >= 0 && e < A.n && p >= A.row_ptr[e] && p < A.row_ptr[e + 1] && A.col[p] == e, "sigma row diagonal position");
    }
    return A;
}
}  // namespace hgf
