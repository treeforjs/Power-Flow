#pragma once
// GPU execution of the deterministic GMG-preconditioned CG (gmg.hpp) on one hall_gpu::Session (hall-gpu-resident-v1
// common: one device, one ordered stream, NVRTC with --fmad=false --prec-div=true --prec-sqrt=true --ftz=false for the
// actual device architecture). The kernels perform exactly the operations of the CPU reference in gmg.hpp in the same
// order; the host reads only the scalar r.r per iteration for the stopping rule.
#include "gmg.hpp"
#include "common/hall_gpu_common.hpp"
#include <chrono>
#include <cstring>
#include <map>

namespace gmg {
namespace detail {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
inline const char kGmgKernels[] = {
#embed "gmg_kernels.cu"
    , 0};
inline const char kGmgCore[] = {
#embed "gmg_core.h"
    , 0};
#pragma GCC diagnostic pop
}  // namespace detail

class GpuSolver {
    hall_gpu::Session& s;
    hall_gpu::Module mod;
    std::map<std::string, hall_gpu::Kernel> k;
    hall_gpu::HostBuffer stage, scal;
    struct DCsr { hall_gpu::DeviceBuffer rp, ci, v; int n = 0; };
    struct DLevel { DCsr K, G, GT, N, P, PT; hall_gpu::DeviceBuffer invd, invdn, x, b, r, rn, phi; int ne = 0, nn = 0; };
    std::vector<DLevel> L;
    hall_gpu::DeviceBuffer sqrtM, inv_sqrtM, fixed0, jac, X, R, Z, Pd, Q, T0, T1, B, part, sc;
    int n0 = 0, nu = 1, nu_coarse = 4; bool have_jac = false;
    std::size_t stage_bytes = std::size_t(64) << 20;
    std::vector<hall_gpu::DeviceBuffer> owned;
    static unsigned blocks(int n) { return unsigned((n + 255) / 256); }
    void put(const hall_gpu::DeviceBuffer& d, const void* src, std::size_t bytes) {
        for (std::size_t off = 0; off < bytes; off += stage_bytes) {
            const std::size_t m = std::min(stage_bytes, bytes - off);
            std::memcpy(stage.pointer, static_cast<const char*>(src) + off, m);
            s.upload(d, off, stage, 0, m);
            s.synchronize("gmg upload");
        }
    }
    template <class T> hall_gpu::DeviceBuffer up(const std::vector<T>& v, const char* tag) {
        auto d = s.alloc(std::max<std::size_t>(v.size() * sizeof(T), 16), tag); owned.push_back(d);
        if (!v.empty()) put(d, v.data(), v.size() * sizeof(T));
        return d;
    }
    DCsr up(const Csr& A, const char* tag) { DCsr d; d.n = A.n; d.rp = up(A.rp, tag); d.ci = up(A.ci, tag); d.v = up(A.v, tag); return d; }
    hall_gpu::DeviceBuffer vec(std::size_t n, const char* tag) { auto d = s.alloc(std::max<std::size_t>(n * 8, 16), tag); owned.push_back(d); s.zero(d, 0, d.bytes); return d; }
    void launch(const char* name, int n, std::vector<void*> args) { s.launch(k.at(name), {blocks(n)}, {256}, 0, args.data()); }
    // kernels
    void residual(const DCsr& A, const hall_gpu::DeviceBuffer& x, const hall_gpu::DeviceBuffer& b, const hall_gpu::DeviceBuffer& r) {
        int n = A.n; auto rp = A.rp.address, ci = A.ci.address, v = A.v.address, xa = x.address, ba = b.address, ra = r.address;
        launch("g_residual", n, {&n, &rp, &ci, &v, &xa, &ba, &ra}); }
    void mul(const DCsr& A, const hall_gpu::DeviceBuffer& x, const hall_gpu::DeviceBuffer& y, bool add) {
        int n = A.n; auto rp = A.rp.address, ci = A.ci.address, v = A.v.address, xa = x.address, ya = y.address;
        launch(add ? "g_mul_add" : "g_mul", n, {&n, &rp, &ci, &v, &xa, &ya}); }
    void jacobi(int n, const hall_gpu::DeviceBuffer& d, const hall_gpu::DeviceBuffer& r, const hall_gpu::DeviceBuffer& x) {
        auto da = d.address, ra = r.address, xa = x.address; launch("g_jacobi", n, {&n, &da, &ra, &xa}); }
    void scale(int n, const hall_gpu::DeviceBuffer& d, const hall_gpu::DeviceBuffer& x, const hall_gpu::DeviceBuffer& y) {
        auto da = d.address, xa = x.address, ya = y.address; launch("g_scale", n, {&n, &da, &xa, &ya}); }
    void zero(int n, const hall_gpu::DeviceBuffer& x) { auto xa = x.address; launch("g_zero", n, {&n, &xa}); }
    void copy(int n, const hall_gpu::DeviceBuffer& a, const hall_gpu::DeviceBuffer& b) { auto aa = a.address, ba = b.address; launch("g_copy", n, {&n, &aa, &ba}); }
    void dot(int n, const hall_gpu::DeviceBuffer& a, const hall_gpu::DeviceBuffer& b, int slot) {
        auto aa = a.address, ba = b.address, pa = part.address, sa = sc.address;
        s.launch(k.at("g_dot"), {512}, {256}, 0, std::vector<void*>{&n, &aa, &ba, &pa}.data());
        s.launch(k.at("g_fin"), {1}, {512}, 0, std::vector<void*>{&pa, &sa, &slot}.data()); }
    void sdiv(int dst, int a, int b) { auto sa = sc.address; s.launch(k.at("g_div"), {1}, {1}, 0, std::vector<void*>{&sa, &dst, &a, &b}.data()); }
    void smov(int dst, int src) { auto sa = sc.address; s.launch(k.at("g_mov"), {1}, {1}, 0, std::vector<void*>{&sa, &dst, &src}.data()); }
    double read(int slot) { s.download(scal, 0, sc, std::size_t(slot) * 8, 8); s.synchronize("gmg scalar"); return *scal.as<double>(); }
    void edge(int i) { auto& l = L[i]; residual(l.K, l.x, l.b, l.r); jacobi(l.ne, l.invd, l.r, l.x); }
    void node(int i) { auto& l = L[i]; residual(l.K, l.x, l.b, l.r); mul(l.GT, l.r, l.rn, false); scale(l.nn, l.invdn, l.rn, l.phi); mul(l.G, l.phi, l.x, true); }
    void vcycle(int i) {
        auto& l = L[i];
        zero(l.ne, l.x);
        if (i + 1 == int(L.size())) {
            for (int q = 0; q < nu_coarse; ++q) { edge(i); node(i); }
            for (int q = 0; q < nu_coarse; ++q) { node(i); edge(i); }
            return;
        }
        for (int q = 0; q < nu; ++q) { edge(i); node(i); }
        residual(l.K, l.x, l.b, l.r); mul(l.PT, l.r, L[i + 1].b, false);
        vcycle(i + 1);
        mul(l.P, L[i + 1].x, l.x, true);
        for (int q = 0; q < nu; ++q) { node(i); edge(i); }
    }
    void prec(const hall_gpu::DeviceBuffer& r, const hall_gpu::DeviceBuffer& z, bool use_jacobi) {
        int n = n0;
        if (use_jacobi) { scale(n, jac, r, z); return; }
        auto f = fixed0.address, sq = sqrtM.address, ra = r.address, b0 = L[0].b.address, x0 = L[0].x.address, za = z.address;
        launch("g_start", n, {&n, &f, &sq, &ra, &b0});
        vcycle(0);
        launch("g_finish", n, {&n, &f, &sq, &ra, &x0, &za});
    }
    void apply(const hall_gpu::DeviceBuffer& p, const hall_gpu::DeviceBuffer& q) {
        int n = n0; scale(n, inv_sqrtM, p, T0); mul(L[0].K, T0, T1, false);
        auto f = fixed0.address, is = inv_sqrtM.address, pa = p.address, t1 = T1.address, qa = q.address;
        launch("g_finish", n, {&n, &f, &is, &pa, &t1, &qa});
    }
public:
    ~GpuSolver() { try { s.synchronize("gmg release"); for (auto& d : owned) s.release(d); s.release(stage); s.release(scal); } catch (...) {} }
    GpuSolver(const GpuSolver&) = delete;
    GpuSolver& operator=(const GpuSolver&) = delete;
    hall_gpu::CompileRecord compile_record;
    double upload_ms = 0;
    GpuSolver(hall_gpu::Session& session, const Hierarchy& H, const Vec* jacobi_inverse) : s(session), nu(H.opt.nu), nu_coarse(H.opt.nu_coarse) {
        const auto t0 = std::chrono::steady_clock::now();
        hall_gpu::CompileRequest rq; rq.main = {"gmg_kernels.cu", detail::kGmgKernels}; rq.headers = {{"gmg_core.h", detail::kGmgCore}};
        mod = s.compile_and_load(rq); compile_record = mod.record;
        for (const char* nm : {"g_residual", "g_mul", "g_mul_add", "g_jacobi", "g_scale", "g_finish", "g_start", "g_zero", "g_copy", "g_cg_update",
                               "g_cg_direction", "g_dot", "g_fin", "g_div", "g_sub", "g_mov"})
            k[nm] = s.kernel(mod, nm);
        stage = s.alloc_pinned(stage_bytes, "gmg staging"); scal = s.alloc_pinned(64, "gmg scalar");
        n0 = H.L[0].g.ne;
        for (const auto& l : H.L) {
            DLevel d; d.ne = l.g.ne; d.nn = l.g.nn;
            d.K = up(l.K, "K"); d.G = up(l.G, "G"); d.GT = up(l.GT, "GT"); d.N = up(l.N, "N");
            if (l.P.n) { d.P = up(l.P, "P"); d.PT = up(l.PT, "PT"); }
            d.invd = up(l.invd, "invd"); d.invdn = up(l.invdn, "invdn");
            d.x = vec(d.ne, "x"); d.b = vec(d.ne, "b"); d.r = vec(d.ne, "r"); d.rn = vec(d.nn, "rn"); d.phi = vec(d.nn, "phi");
            L.push_back(d);
        }
        sqrtM = up(H.sqrtM, "sqrtM"); inv_sqrtM = up(H.inv_sqrtM, "inv sqrtM"); fixed0 = up(H.fixed0, "fixed");
        if (jacobi_inverse) { jac = up(*jacobi_inverse, "jacobi"); have_jac = true; }
        for (auto* b : {&X, &R, &Z, &Pd, &Q, &T0, &T1, &B}) *b = vec(n0, "cg");
        part = vec(512, "partials"); sc = vec(16, "scalars");
        s.synchronize("gmg setup");
        upload_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    }
    // Same scalar sequence as gmg::cpu_solve.
    Vec solve(const Vec& b, int maxit, bool use_jacobi, SolveStats& st) {
        need(int(b.size()) == n0, "rhs size"); need(!use_jacobi || have_jac, "no Jacobi diagonal uploaded");
        int n = n0;
        put(B, b.data(), b.size() * 8);
        dot(n, B, B, 6); st.rhs_norm = std::sqrt(read(6)); st.target = std::max(1e-30, 2e-13 * st.rhs_norm);
        zero(n, X); copy(n, B, R);
        auto sa = sc.address, xa = X.address, ra = R.address, za = Z.address, pa = Pd.address, qa = Q.address, ba = B.address;
        for (int corr = 0; corr < 3; ++corr) {
            st.corrections = corr + 1;
            prec(R, Z, use_jacobi); copy(n, Z, Pd); dot(n, R, Z, 0); dot(n, R, R, 5);
            double rr = read(5); int it = 0;
            while (it < maxit && std::sqrt(rr) > .2 * st.target) {
                apply(Pd, Q); dot(n, Pd, Q, 1); sdiv(2, 0, 1);
                int slot = 2; launch("g_cg_update", n, {&n, &sa, &slot, &pa, &qa, &xa, &ra});
                dot(n, R, R, 5); rr = read(5); ++it;
                if (!(std::sqrt(rr) > .2 * st.target)) break;
                prec(R, Z, use_jacobi); dot(n, R, Z, 3); sdiv(4, 3, 0); smov(0, 3);
                int bslot = 4; launch("g_cg_direction", n, {&n, &sa, &bslot, &za, &pa});
            }
            st.iterations += it; st.cycle_iterations.push_back(it);
            apply(X, Q); launch("g_sub", n, {&n, &ba, &qa, &ra});
            dot(n, R, R, 5); st.backend_residual = std::sqrt(read(5));
            if (st.backend_residual <= st.target) break;
        }
        Vec x(n);
        for (std::size_t off = 0; off < std::size_t(n) * 8; off += stage_bytes) {
            const std::size_t m = std::min(stage_bytes, std::size_t(n) * 8 - off);
            s.download(stage, 0, X, off, m); s.synchronize("gmg download");
            std::memcpy(reinterpret_cast<char*>(x.data()) + off, stage.pointer, m);
        }
        return x;
    }
    // One V-cycle application (timing of the preconditioner alone).
    double time_vcycle(int reps) {
        s.synchronize("gmg"); const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < reps; ++i) prec(R, Z, false);
        s.synchronize("gmg vcycle"); return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() / reps;
    }
    double time_apply(int reps) {
        s.synchronize("gmg"); const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < reps; ++i) apply(Pd, Q);
        s.synchronize("gmg apply"); return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() / reps;
    }
};
}  // namespace gmg
