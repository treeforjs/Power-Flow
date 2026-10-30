// Device-resident deterministic FP64 Jacobi CG for the Hall edge operator.
// Stopping rule identical to foil::solve (foil_em.hpp:56-63): inner loop while
// iteration < maxit and ||r_recursive|| > 0.2 * target, then the true residual
// ||b - A x|| is recomputed on the device and compared with
// target = max(1e-30, 2e-13 ||b||); at most three correction cycles.
//
// Launch structure: one "head" (target, correction-0 init, BATCH gated
// iterations, true residual into a separate buffer) followed by one control
// readback. Further gated batches run only if the stopping rule has not fired;
// gated kernels are no-ops after it fires, so the result does not depend on
// BATCH or on graph use. The true residual never overwrites the recursive
// residual; a correction cycle starts from it explicitly (as foil::solve).
#pragma once
#include "cu_driver.hpp"
#include "kernels.hpp"
#include <map>
#include <tuple>

namespace hgf {
enum class OpMode { assembled, factored };
inline const char* name(OpMode m) { return m == OpMode::assembled ? "assembled" : "factored"; }

struct SolveStats {
    std::vector<int> iterations;
    int corrections = 0, readbacks = 0;
    double bnorm = 0, target = 0, true_norm = 0;
    bool converged = false;
};

class GpuSolver {
public:
    GpuSolver(cu::Session& s, const Operator& A, int maxit) : s_(s), n_(A.n), nb_(A.nb), maxit_(maxit), q_(A.q) {
        auto dev = [&](auto& ptr, const auto& v) {
            using T = typename std::decay_t<decltype(v)>::value_type;
            ptr = s_.device(v.size() * sizeof(T));
            if (!v.empty()) s_.upload_static(ptr, v.data(), v.size() * sizeof(T));
        };
        dev(rp_, A.row_ptr); dev(ci_, A.col); dev(val_, A.val);
        dev(fe_, A.face_edge); dev(fc_, A.face_coef); dev(H_, A.H);
        dev(ctp_, A.ct_ptr); dev(ctf_, A.ct_face); dev(ctc_, A.ct_coef);
        dev(sqrtM_, A.sqrtM); dev(a_, A.a); dev(fixed_, A.fixed); dev(diag_, A.diag);
        dev(srows_, A.sigma_rows); dev(spos_, A.sigma_diag_pos);
        m_ = int(A.sigma_rows.size());
        const std::size_t vb = std::size_t(n_) * 8;
        b_ = s_.device(vb); x_ = s_.device(vb); r_ = s_.device(vb); rt_ = s_.device(vb); z_ = s_.device(vb); p_ = s_.device(vb); Ap_ = s_.device(vb);
        Bh_ = s_.device(std::size_t(nb_) * 8);
        part1_ = s_.device(NBLK * 8); part2_ = s_.device(NBLK * 8);
        ctl_ = s_.device(sizeof(Ctl));
        ha_ = s_.device(std::size_t(3) * maxit_ * 8); hb_ = s_.device(std::size_t(3) * maxit_ * 8);
        hr_ = s_.device(std::size_t(3) * (maxit_ + 1) * 8);
        snew_diag_ = s_.device(std::size_t(m_) * 8); snew_a_ = s_.device(std::size_t(m_) * 8);
        s_.zero(ha_, std::size_t(3) * maxit_ * 8); s_.zero(hb_, std::size_t(3) * maxit_ * 8); s_.zero(hr_, std::size_t(3) * (maxit_ + 1) * 8);
        ctl_host_ = static_cast<Ctl*>(s_.pinned(sizeof(Ctl)));
        b_host_ = static_cast<double*>(s_.pinned(vb)); x_host_ = static_cast<double*>(s_.pinned(vb));
        hist_host_ = static_cast<double*>(s_.pinned(std::size_t(3) * (maxit_ + 1) * 8));
        val_host_ = static_cast<double*>(s_.pinned(A.val.size() * 8)); diag_host_ = static_cast<double*>(s_.pinned(vb));
        a_host_ = static_cast<double*>(s_.pinned(vb));
        sdiag_host_ = static_cast<double*>(s_.pinned(std::size_t(m_) * 8)); sa_host_ = static_cast<double*>(s_.pinned(std::size_t(m_) * 8));
        std::memcpy(val_host_, A.val.data(), A.val.size() * 8);
        std::memcpy(diag_host_, A.diag.data(), vb);
        std::memcpy(a_host_, A.a.data(), vb);
        for (int k = 0; k < m_; ++k) { sdiag_host_[k] = A.diag[A.sigma_rows[k]]; sa_host_[k] = A.a[A.sigma_rows[k]]; }
        nnz_ = A.val.size();
        for (const char* k : {"k_dot_self", "k_fin_target", "k_init", "k_fin_init", "k_spmv_csr", "k_face", "k_edge", "k_fin_alpha",
                              "k_update", "k_fin_beta", "k_direction", "k_residual", "k_fin_true", "k_refresh"})
            fn_[k] = s_.function(k);
        s_.sync();  // static operator resident and staging drained before use
    }

    // b staged in pinned memory (outside any timer), then copied on the stream.
    void stage_rhs(const Vec& b) { need(b.size() == std::size_t(n_), "rhs size"); std::memcpy(b_host_, b.data(), b.size() * 8); }
    void upload_rhs() { s_.h2d(b_, b_host_, std::size_t(n_) * 8); }
    void download_x() { s_.d2h(x_host_, x_, std::size_t(n_) * 8); }
    Vec x_host() const { return Vec(x_host_, x_host_ + n_); }
    // Upper bound of a sigma(T)/dt refresh: all CSR values, Jacobi diagonal and a.
    void upload_operator_values() {
        s_.h2d(val_, val_host_, nnz_ * 8);
        s_.h2d(diag_, diag_host_, std::size_t(n_) * 8);
        s_.h2d(a_, a_host_, std::size_t(n_) * 8);
    }
    std::size_t operator_value_bytes() const { return nnz_ * 8 + std::size_t(n_) * 16; }
    // sigma(T) refresh as in a time step: new diagonal and a for the conduction
    // rows only (host-computed values), scattered on the device.
    void refresh_sigma_rows() {
        if (m_ == 0) return;
        s_.h2d(snew_diag_, sdiag_host_, std::size_t(m_) * 8);
        s_.h2d(snew_a_, sa_host_, std::size_t(m_) * 8);
        int m = m_;
        launch("k_refresh", unsigned((m_ + TPB - 1) / TPB), TPB, {&m, &srows_, &spos_, &snew_diag_, &snew_a_, &val_, &diag_, &a_});
    }
    std::size_t sigma_refresh_bytes() const { return std::size_t(m_) * 16; }
    int sigma_rows() const { return m_; }

    // y = A v on the device for a host probe vector (verification only).
    Vec product(OpMode m, const Vec& v) {
        s_.upload_static(p_, v.data(), v.size() * 8);
        apply(m, p_, Ap_, 0, 0);
        s_.d2h(x_host_, Ap_, std::size_t(n_) * 8);
        s_.sync();
        return Vec(x_host_, x_host_ + n_);
    }

    SolveStats solve(OpMode m, int batch, bool graph) {
        SolveStats st;
        if (graph) s_.launch_graph(graph_for(m, batch, true));
        else enqueue_head(m, batch);
        read(st);
        for (int corr = 0; corr < 3; ++corr) {
            st.corrections = corr + 1;
            bool tail = false;
            if (corr > 0) {
                int n = n_;
                s_.d2d(r_, rt_, std::size_t(n_) * 8);
                launch("k_init", NBLK, TPB, {&n, &r_, &diag_, &z_, &p_, &part1_, &part2_});
                launch("k_fin_init", 1, NBLK, {&part1_, &part2_, &ctl_, &corr, &hr_});
                ctl_host_->done = 0;
                tail = true;
            }
            for (int guard = 0; !ctl_host_->done; ++guard) {
                need(guard <= maxit_ / batch + 2, "device CG did not stop within maxit");
                if (graph) s_.launch_graph(graph_for(m, batch, false));
                else for (int k = 0; k < batch; ++k) enqueue_iteration(m);
                read(st);
                tail = true;
            }
            if (tail) { enqueue_tail(m); read(st); }
            st.iterations.push_back(ctl_host_->it);
            st.bnorm = ctl_host_->bnorm; st.target = ctl_host_->target; st.true_norm = ctl_host_->tnorm;
            if (ctl_host_->tnorm <= ctl_host_->target) { st.converged = true; return st; }
        }
        throw std::runtime_error("field solve independent residual failed (device CG)");
    }

    // Coefficient and residual histories of the last solve.
    void histories(const SolveStats& st, std::vector<std::vector<double>>& alpha, std::vector<std::vector<double>>& beta,
                   std::vector<std::vector<double>>& rnorm) {
        Vec ha(std::size_t(3) * maxit_), hb(ha.size()), hr(std::size_t(3) * (maxit_ + 1));
        double* tmp = hist_host_;
        s_.d2h(tmp, ha_, ha.size() * 8); s_.sync(); std::memcpy(ha.data(), tmp, ha.size() * 8);
        s_.d2h(tmp, hb_, hb.size() * 8); s_.sync(); std::memcpy(hb.data(), tmp, hb.size() * 8);
        s_.d2h(tmp, hr_, hr.size() * 8); s_.sync(); std::memcpy(hr.data(), tmp, hr.size() * 8);
        alpha.clear(); beta.clear(); rnorm.clear();
        for (int c = 0; c < st.corrections; ++c) {
            const int k = st.iterations[c];
            alpha.emplace_back(ha.begin() + c * maxit_, ha.begin() + c * maxit_ + k);
            beta.emplace_back(hb.begin() + c * maxit_, hb.begin() + c * maxit_ + k);
            rnorm.emplace_back(hr.begin() + c * (maxit_ + 1), hr.begin() + c * (maxit_ + 1) + k + 1);
        }
    }

private:
    void read(SolveStats& st) {
        s_.d2h(ctl_host_, ctl_, sizeof(Ctl));
        s_.sync();
        ++st.readbacks;
        if (ctl_host_->error) throw std::runtime_error("non-SPD field operator (device CG)");
    }
    void launch(const char* k, unsigned grid, unsigned block, std::initializer_list<void*> args) {
        std::vector<void*> a(args);
        s_.launch(fn_.at(k), grid, block, a.data());
    }
    void apply(OpMode m, CUdeviceptr& v, CUdeviceptr& out, int gated, int with_dot) {
        int n = n_, nb = nb_;
        if (m == OpMode::assembled) {
            launch("k_spmv_csr", NBLK, TPB, {&n, &rp_, &ci_, &val_, &v, &out, &part1_, &ctl_, &gated, &with_dot});
        } else {
            launch("k_face", unsigned((nb_ + TPB - 1) / TPB), TPB, {&nb, &fe_, &fc_, &H_, &sqrtM_, &fixed_, &v, &Bh_, &ctl_, &gated});
            launch("k_edge", NBLK, TPB, {&n, &ctp_, &ctf_, &ctc_, &Bh_, &sqrtM_, &a_, &q_, &fixed_, &v, &out, &part1_, &ctl_, &gated, &with_dot});
        }
    }
    void enqueue_iteration(OpMode m) {
        int n = n_;
        apply(m, p_, Ap_, 1, 1);
        launch("k_fin_alpha", 1, NBLK, {&part1_, &ctl_, &ha_});
        launch("k_update", NBLK, TPB, {&n, &x_, &r_, &p_, &Ap_, &diag_, &z_, &part1_, &part2_, &ctl_});
        launch("k_fin_beta", 1, NBLK, {&part1_, &part2_, &ctl_, &hb_, &hr_});
        launch("k_direction", unsigned((n_ + TPB - 1) / TPB), TPB, {&n, &z_, &p_, &ctl_});
    }
    void enqueue_tail(OpMode m) {
        int n = n_;
        apply(m, x_, Ap_, 0, 0);
        launch("k_residual", NBLK, TPB, {&n, &b_, &Ap_, &rt_, &part1_});
        launch("k_fin_true", 1, NBLK, {&part1_, &ctl_});
    }
    void enqueue_head(OpMode m, int batch) {
        int n = n_, maxit = maxit_, corr = 0;
        launch("k_dot_self", NBLK, TPB, {&n, &b_, &part1_});
        launch("k_fin_target", 1, NBLK, {&part1_, &ctl_, &maxit});
        s_.zero(x_, std::size_t(n_) * 8);
        s_.d2d(r_, b_, std::size_t(n_) * 8);
        launch("k_init", NBLK, TPB, {&n, &r_, &diag_, &z_, &p_, &part1_, &part2_});
        launch("k_fin_init", 1, NBLK, {&part1_, &part2_, &ctl_, &corr, &hr_});
        for (int k = 0; k < batch; ++k) enqueue_iteration(m);
        enqueue_tail(m);
    }
    CUgraphExec graph_for(OpMode m, int batch, bool head) {
        const auto key = std::make_tuple(int(m), batch, head);
        auto it = graphs_.find(key);
        if (it != graphs_.end()) return it->second;
        s_.begin_capture();
        try {
            if (head) enqueue_head(m, batch);
            else for (int k = 0; k < batch; ++k) enqueue_iteration(m);
        } catch (...) {
            CUgraph g = nullptr;
            cuStreamEndCapture(s_.stream(), &g);
            if (g) cuGraphDestroy(g);
            throw;
        }
        auto e = s_.end_capture();
        graphs_[key] = e;
        return e;
    }

    cu::Session& s_;
    int n_, nb_, maxit_, m_ = 0;
    double q_;
    std::size_t nnz_ = 0;
    CUdeviceptr rp_ = 0, ci_ = 0, val_ = 0, fe_ = 0, fc_ = 0, H_ = 0, ctp_ = 0, ctf_ = 0, ctc_ = 0, sqrtM_ = 0, a_ = 0, fixed_ = 0, diag_ = 0;
    CUdeviceptr srows_ = 0, spos_ = 0, snew_diag_ = 0, snew_a_ = 0;
    CUdeviceptr b_ = 0, x_ = 0, r_ = 0, rt_ = 0, z_ = 0, p_ = 0, Ap_ = 0, Bh_ = 0, part1_ = 0, part2_ = 0, ctl_ = 0, ha_ = 0, hb_ = 0, hr_ = 0;
    Ctl* ctl_host_ = nullptr;
    double *b_host_ = nullptr, *x_host_ = nullptr, *hist_host_ = nullptr, *val_host_ = nullptr, *diag_host_ = nullptr, *a_host_ = nullptr;
    double *sdiag_host_ = nullptr, *sa_host_ = nullptr;
    std::map<std::string, CUfunction> fn_;
    std::map<std::tuple<int, int, bool>, CUgraphExec> graphs_;
};
}  // namespace hgf
