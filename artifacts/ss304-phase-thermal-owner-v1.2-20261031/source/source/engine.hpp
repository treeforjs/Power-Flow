// engine.hpp -- execution backends of the phase owner v1.2. The CPU reference runs the phase_kernels.h functions
// in index loops; the GPU engine (gpu_engine.hpp) runs the same functions as NVRTC kernels in the same sequence,
// so the trial state words are bit-identical.
#pragma once
#include "phase_law_host.hpp"
#include "state_digest.h"
#include <algorithm>
#include <memory>
#include <vector>

namespace ss304_phase_v1_2 {

struct HostState {
    std::vector<px_i128> E, N;   // [n_cells], [PX_NEL n_cells]
    std::vector<int> flags;      // bit0: void (merged away)
    std::vector<px_i128> fcm;    // [n_fc PX_NEL]: atoms gained through each large-face column (displacement)
    std::vector<double> fcV;     // [3 n_fc]: layer momentum x, z and layer mass (depth-averaged model)
    std::vector<px_i128> fcK;    // [n_fc]: layer kinetic energy word
};
struct HostGeometry {
    std::vector<double> row_dy, gx, gy, vref, mref;
    std::vector<int> row_em, em_row_begin;
    std::vector<int> if_a, if_b, if_dir, fc_if, slice_zn;   // depth-averaged layer topology
};
struct HostElements {
    std::vector<int> el_of, el_cell, el_face, el_patch, patch_begin;
    std::vector<double> el_area;
    bool operator==(const HostElements&) const = default;
};
struct StepData {
    bool use_acc = false;                     // Joule energy and Lorentz force from the engine's accumulators
    std::vector<double> em_energy;            // [n_em] J
    std::vector<double> em_JB;                // [6 n_em] J (A/m2) and B (T) per EM cell
    std::vector<double> patch_pb;             // [n_patches] back pressure (Pa)
    std::vector<int> sc_el;                   // distinct elements with loads/deposits
    std::vector<px_i128> sc_val;              // [PX_SC per entry]: load_E, dep_N x6, dep_E
    std::vector<px_i128> rel_total_N;         // [PX_NEL n_patches]
    std::vector<double> rel_kin_J;            // [PX_NEL n_patches]
};
struct EngineResult {
    int trial_slot = -1;
    std::vector<px_i128> slice_sum;           // [S*PX_SLICE_SUM]
    std::vector<int> slice_status;            // [S*4]
    std::vector<double> slice_maxT;           // [S]
    std::vector<px_i128> patch_sum;           // [P*PX_PATCH_SUM]
    std::vector<double> patch_out;            // [P*PX_PATCH_OUT]
    std::vector<double> slice_diag;           // [4 S]
    std::vector<double> sigma_em;             // [n_em] at the trial state (step) / evaluated state
};

class Engine {
public:
    virtual ~Engine() = default;
    virtual std::string name() const = 0;
    virtual void configure(const px_law& law, const px_dims& dims, const HostGeometry& g) = 0;
    virtual void set_elements(const HostElements& el) = 0;
    virtual void write_state(int slot, const HostState& s) = 0;
    virtual HostState read_state(int slot) = 0;
    virtual void evaluate(int slot, const px_dims& d, const std::vector<double>& patch_pb, EngineResult& r, std::vector<double>* sigma_em) = 0;
    virtual void step(int from, const px_dims& d, const StepData& in, EngineResult& r) = 0;
    virtual std::vector<double> temperatures(int slot, const px_dims& d) = 0;
    virtual std::vector<double> face_displacements(int slot, const px_dims& d) = 0;   // fc_eta of the slot
    // Block-Merkle state digest of a slot (state_digest.h); leaves optionally returned for host sampling.
    virtual PxStateRoots state_roots(int slot, const px_dims& d, std::vector<std::array<px_u32, 8>>* cell_leaves) = 0;
    // Sub-cycling accumulators: one field step (power W and J, B per EM cell; JB may be empty) into exact Joule words
    // and (J x B) dt; reset after an accepted owner step.
    virtual void accumulate(const std::vector<double>& power, const std::vector<double>& JB, double dt) = 0;
    virtual void reset_accumulator() = 0;
    // Field-step data already on the device of the engine's session (GPU engine only).
    virtual void accumulate_device(const double* power, const double* JB, double dt) { (void)power; (void)JB; (void)dt; throw std::logic_error("device accumulation needs the GPU engine"); }
};

inline std::pair<int, int> work_slots(int from) {
    int a = -1, b = -1;
    for (int s = 0; s < 3; ++s) if (s != from) { if (a < 0) a = s; else b = s; }
    return {a, b};
}

// ---------------------------------------------------------------------------------------------------------
class CpuEngine final : public Engine {
    px_law law_{};
    px_dims base_{};
    HostGeometry g_;
    HostElements el_;
    HostState slot_[3];
    std::vector<double> T_, k_, sig_, m_, quote_, patch_out_, slice_maxT_, sigma_em_, em_energy_, em_JB_, patch_pb_, rel_kin_;
    std::vector<double> fc_delta_, fc_T_, fc_eta_, fc_p_, fc_x_, ic_q_;
    std::vector<int> status_, slice_status_, fc_cell_;
    std::vector<px_i128> joule_, in_load_, in_dep_N_, in_dep_E_, in_rel_N_, in_rel_E_, acc_ev_N_, acc_ev_Eref_, acc_ev_Ekin_,
        acc_ev_Eint_, acc_bnd_, acc_rad_, patch_sum_, slice_sum_, rel_total_, ic_N_, ic_E_;
    // v1.2
    std::vector<px_i128> acc_E_, if_N_, if_E_, if_K_, cflow_N_, cflow_E_, fc_netN_, acc_fc_, ic_W_, acc_ic_;
    std::vector<double> acc_F_, em_F_, em_diag_, slice_diag_, fc_m_, fc_h_, fc_rho_, fc_mu_, fc_dsdT_, fc_sig_, fc_f_, fc_u_, fc_un_,
        fc_a_, fc_ps_, fc_g_, fc_re_, fc_s_, if_dm_, if_F_;
    std::vector<int> acc_bad_, fc_lay_, fc_nw_, if_st_;
    std::vector<px_i128> el_dE_, el_dN_, row_sum_;
    std::vector<int> el_st_, el_bad_, row_info_;
    std::vector<double> row_maxT_;

    px_buffers buffers(int cur, int next) {
        px_buffers b{};
        b.law = &law_;
        b.row_dy = g_.row_dy.data(); b.row_em = g_.row_em.data(); b.em_row_begin = g_.em_row_begin.data();
        b.gx = g_.gx.data(); b.gy = g_.gy.data(); b.vref = g_.vref.data(); b.mref = g_.mref.data();
        b.E = slot_[cur].E.data(); b.N = slot_[cur].N.data(); b.flags = slot_[cur].flags.data(); b.fcm = slot_[cur].fcm.data();
        if (next >= 0) { b.E_out = slot_[next].E.data(); b.N_out = slot_[next].N.data(); b.fcm_out = slot_[next].fcm.data(); }
        b.T = T_.data(); b.k = k_.data(); b.sig = sig_.data(); b.m = m_.data(); b.status = status_.data();
        b.joule_total = joule_.data(); b.em_energy = em_energy_.data(); b.em_JB = em_JB_.data();
        b.el_of = el_.el_of.data(); b.el_cell = el_.el_cell.data(); b.el_face = el_.el_face.data(); b.el_area = el_.el_area.data();
        b.el_patch = el_.el_patch.data(); b.patch_begin = el_.patch_begin.data();
        b.in_load_E = in_load_.data(); b.in_dep_N = in_dep_N_.data(); b.in_dep_E = in_dep_E_.data();
        b.in_rel_N = in_rel_N_.data(); b.in_rel_E = in_rel_E_.data();
        b.acc_ev_N = acc_ev_N_.data(); b.acc_ev_Eref = acc_ev_Eref_.data(); b.acc_ev_Ekin = acc_ev_Ekin_.data(); b.acc_ev_Eint = acc_ev_Eint_.data();
        b.acc_bnd_E = acc_bnd_.data(); b.acc_rad = acc_rad_.data();
        b.el_quote = quote_.data(); b.patch_pb = patch_pb_.data();
        b.fc_cell = fc_cell_.data(); b.fc_delta = fc_delta_.data(); b.fc_T = fc_T_.data(); b.fc_eta = fc_eta_.data(); b.fc_p = fc_p_.data(); b.fc_x = fc_x_.data();
        b.ic_N = ic_N_.data(); b.ic_E = ic_E_.data(); b.ic_q = ic_q_.data();
        b.rel_total_N = rel_total_.data(); b.rel_kin_J = rel_kin_.data();
        b.patch_out = patch_out_.data(); b.patch_sum = patch_sum_.data(); b.slice_sum = slice_sum_.data(); b.slice_status = slice_status_.data();
        b.slice_maxT = slice_maxT_.data(); b.sigma_em = sigma_em_.data();
        b.in_load_E_w = in_load_.data(); b.in_dep_N_w = in_dep_N_.data(); b.in_dep_E_w = in_dep_E_.data();
        b.fcV = slot_[cur].fcV.data(); b.fcK = slot_[cur].fcK.data();
        if (next >= 0) { b.fcV_out = slot_[next].fcV.data(); b.fcK_out = slot_[next].fcK.data(); }
        b.acc_E = acc_E_.data(); b.acc_F = acc_F_.data(); b.acc_bad = acc_bad_.data(); b.em_F = em_F_.data(); b.em_diag = em_diag_.data();
        b.slice_diag = slice_diag_.data();
        b.if_a = g_.if_a.data(); b.if_b = g_.if_b.data(); b.if_dir = g_.if_dir.data(); b.fc_if = g_.fc_if.data(); b.slice_zn = g_.slice_zn.data();
        b.fc_lay = fc_lay_.data(); b.fc_nw = fc_nw_.data(); b.fc_m = fc_m_.data(); b.fc_h = fc_h_.data(); b.fc_rho = fc_rho_.data(); b.fc_mu = fc_mu_.data();
        b.fc_dsdT = fc_dsdT_.data(); b.fc_sig = fc_sig_.data(); b.fc_f = fc_f_.data(); b.fc_u = fc_u_.data(); b.fc_un = fc_un_.data();
        b.fc_a = fc_a_.data(); b.fc_ps = fc_ps_.data(); b.fc_g = fc_g_.data(); b.fc_re = fc_re_.data(); b.fc_s = fc_s_.data();
        b.if_dm = if_dm_.data(); b.if_F = if_F_.data(); b.if_st = if_st_.data(); b.if_N = if_N_.data(); b.if_E = if_E_.data(); b.if_K = if_K_.data();
        b.cflow_N = cflow_N_.data(); b.cflow_E = cflow_E_.data(); b.fc_netN = fc_netN_.data(); b.acc_fc = acc_fc_.data();
        b.ic_W = ic_W_.data(); b.acc_ic = acc_ic_.data();
        b.el_dE = el_dE_.data(); b.el_dN = el_dN_.data(); b.el_st = el_st_.data(); b.el_bad = el_bad_.data();
        b.row_sum = row_sum_.data(); b.row_info = row_info_.data(); b.row_maxT = row_maxT_.data();
        return b;
    }
    void sized_elements() {
        const std::size_t n = el_.el_cell.size(), P = el_.patch_begin.size() - 1;
        for (auto* v : {&in_load_, &in_dep_E_, &in_rel_E_, &acc_ev_Eref_, &acc_ev_Ekin_, &acc_ev_Eint_, &acc_bnd_, &acc_rad_}) v->assign(std::max<std::size_t>(n, 1), px_zero());
        for (auto* v : {&in_dep_N_, &in_rel_N_, &acc_ev_N_}) v->assign(PX_NEL * std::max<std::size_t>(n, 1), px_zero());
        quote_.assign(PX_NEL * std::max<std::size_t>(n, 1), 0.0);
        el_dE_.assign(std::max<std::size_t>(n, 1), px_zero()); el_dN_.assign(PX_NEL * std::max<std::size_t>(n, 1), px_zero());
        el_st_.assign(std::max<std::size_t>(n, 1), 0); el_bad_.assign(std::max<std::size_t>(n, 1), 0);
        patch_out_.assign(PX_PATCH_OUT * P, 0.0);
        patch_sum_.assign(PX_PATCH_SUM * P, px_zero());
        rel_total_.assign(PX_NEL * P, px_zero());
        rel_kin_.assign(PX_NEL * P, 0.0);
    }
    void collect(EngineResult& r) {
        r.slice_sum = slice_sum_; r.slice_status = slice_status_; r.slice_maxT = slice_maxT_;
        r.patch_sum = patch_sum_; r.patch_out = patch_out_; r.slice_diag = slice_diag_; r.sigma_em = sigma_em_;
    }
    void flow_pass(const px_dims& d, const px_buffers& b) {
        if (d.flow_model == 2) {
            for (int fc = 0; fc < d.n_fc; ++fc) px_fc_prepare2(&d, &b, fc);
            for (int fc = 0; fc < d.n_fc; ++fc) px_fc_pressure2(&d, &b, fc);
            for (int q = 0; q < d.n_if; ++q) px_if_flux2(&d, &b, q);
            for (int fc = 0; fc < d.n_fc; ++fc) px_fc_pass(&d, &b, fc);
            return;
        }
        for (int fc = 0; fc < d.n_fc; ++fc) px_fc_prepare(&d, &b, fc);
        for (int fc = 0; fc < d.n_fc; ++fc) px_fc_pressure(&d, &b, fc);
        for (int ic = 0; ic < d.n_ic; ++ic) px_ic_flux(&d, &b, ic);
    }
    // full: patch report with the surface state (evaluate); otherwise the step sums only
    void finish(int slot, const px_dims& d, bool full) {
        px_buffers b = buffers(slot, -1);
        for (int c = 0; c < d.n_cells; ++c) px_thermo(&d, &b, c);
        if (d.flow && d.flow_model == 2) for (int fc = 0; fc < d.n_fc; ++fc) px_fc_prepare2(&d, &b, fc);
        else for (int fc = 0; fc < d.n_fc; ++fc) px_fc_prepare(&d, &b, fc);
        for (int e = 0; e < d.n_em; ++e) px_sigma_em(&d, &b, e);
        for (int r = 0; r < d.S * d.ny; ++r) px_row_reduce(&d, &b, r);
        for (int s = 0; s < d.S; ++s) px_slice_reduce(&d, &b, s);
        if (full) for (int p = 0; p < d.n_patches; ++p) px_patch_report(&d, &b, p);
        else for (int p = 0; p < d.n_patches; ++p) px_patch_sums(&d, &b, p);
    }

public:
    std::string name() const override { return "cpu-reference"; }
    void configure(const px_law& law, const px_dims& dims, const HostGeometry& g) override {
        law_ = law; base_ = dims; g_ = g;
        const std::size_t n = dims.n_cells, nfc = std::max(dims.n_fc, 1), nic = std::max(dims.n_ic, 1);
        for (auto& s : slot_) { s.E.assign(n, px_zero()); s.N.assign(PX_NEL * n, px_zero()); s.flags.assign(n, 0); s.fcm.assign(PX_NEL * nfc, px_zero()); s.fcV.assign(3 * nfc, 0.0); s.fcK.assign(nfc, px_zero()); }
        T_.assign(n, 0); k_.assign(n, 0); sig_.assign(n, 0); m_.assign(n, 0); status_.assign(n, 0); joule_.assign(n, px_zero());
        slice_sum_.assign(PX_SLICE_SUM * dims.S, px_zero()); slice_status_.assign(4 * dims.S, 0); slice_maxT_.assign(dims.S, 0);
        sigma_em_.assign(dims.n_em, 0); em_energy_.assign(dims.n_em, 0); em_JB_.assign(6 * dims.n_em, 0.0);
        patch_pb_.assign(std::max(dims.n_patches, 1), 0.0);
        fc_cell_.assign(nfc, -1); fc_delta_.assign(nfc, 0); fc_T_.assign(nfc, 0); fc_eta_.assign(nfc, 0); fc_p_.assign(nfc, 0); fc_x_.assign(PX_NEL * nfc, 0);
        ic_N_.assign(PX_NEL * nic, px_zero()); ic_E_.assign(nic, px_zero()); ic_q_.assign(nic, 0);
        const std::size_t nem = std::max(dims.n_em, 1), nif = std::max(dims.n_if, 1);
        acc_E_.assign(nem, px_zero()); acc_F_.assign(3 * nem, 0.0); acc_bad_.assign(nem, 0); em_F_.assign(3 * nem, 0.0); em_diag_.assign(2 * nem, 0.0);
        slice_diag_.assign(4 * std::size_t(dims.S), 0.0);
        for (auto* v : {&fc_m_, &fc_h_, &fc_rho_, &fc_mu_, &fc_dsdT_, &fc_sig_, &fc_a_, &fc_ps_, &fc_g_, &fc_re_, &fc_s_}) v->assign(nfc, 0.0);
        fc_f_.assign(3 * nfc, 0.0); fc_u_.assign(2 * nfc, 0.0); fc_un_.assign(2 * nfc, 0.0);
        fc_lay_.assign(nfc, 0); fc_nw_.assign(nfc, 0);
        if_dm_.assign(nif, 0.0); if_F_.assign(3 * nif, 0.0); if_st_.assign(nif, 0); if_N_.assign(PX_NEL * nif, px_zero()); if_E_.assign(nif, px_zero()); if_K_.assign(nif, px_zero());
        cflow_N_.assign(PX_NEL * n, px_zero()); cflow_E_.assign(n, px_zero()); fc_netN_.assign(PX_NEL * nfc, px_zero());
        acc_fc_.assign(PX_FC_ACC * nfc, px_zero()); ic_W_.assign(nic, px_zero()); acc_ic_.assign(2 * nic, px_zero());
        row_sum_.assign(PX_ROW_SUM * std::size_t(dims.S) * dims.ny, px_zero()); row_info_.assign(5 * std::size_t(dims.S) * dims.ny, 0);
        row_maxT_.assign(std::size_t(dims.S) * dims.ny, 0.0);
        if (g_.slice_zn.empty()) g_.slice_zn.assign(2 * std::size_t(dims.S), -1);
        if (g_.fc_if.empty()) g_.fc_if.assign(4 * nfc, -1);
        for (auto* v : {&g_.if_a, &g_.if_b, &g_.if_dir}) if (v->empty()) v->assign(1, -1);
    }
    void set_elements(const HostElements& el) override { el_ = el; sized_elements(); }
    void write_state(int slot, const HostState& s) override { slot_[slot] = s; }
    HostState read_state(int slot) override { return slot_[slot]; }
    void evaluate(int slot, const px_dims& d0, const std::vector<double>& pb, EngineResult& r, std::vector<double>* sigma_em) override {
        px_dims d = d0; d.use_acc = 0; d.diag = 0;
        std::fill(status_.begin(), status_.end(), 0);
        patch_pb_ = pb;
        finish(slot, d, true);
        if (sigma_em) *sigma_em = sigma_em_;
        collect(r);
    }
    void step(int from, const px_dims& d0, const StepData& in, EngineResult& r) override {
        px_dims d = d0;
        d.use_acc = in.use_acc ? 1 : 0;
        std::fill(status_.begin(), status_.end(), 0);
        sized_elements();
        std::fill(acc_fc_.begin(), acc_fc_.end(), px_zero()); std::fill(acc_ic_.begin(), acc_ic_.end(), px_zero());
        std::fill(fc_s_.begin(), fc_s_.end(), 0.0);
        if (!in.use_acc) { em_energy_ = in.em_energy; em_JB_ = in.em_JB; }
        patch_pb_ = in.patch_pb;
        rel_total_ = in.rel_total_N; rel_kin_ = in.rel_kin_J;
        {
            px_buffers b = buffers(from, -1);
            for (int c = 0; c < d.n_cells; ++c) px_thermo(&d, &b, c);
            px_dims ds = d; ds.n_scatter = int(in.sc_el.size());
            b.sc_el = in.sc_el.data(); b.sc_val = in.sc_val.data();
            for (int q = 0; q < ds.n_scatter; ++q) px_scatter(&ds, &b, q);
            if (d.evap_mode == 2) {
                for (int e = 0; e < d.n_el; ++e) px_quote_element(&d, &b, e);
                for (int p = 0; p < d.n_patches; ++p) px_release_split(&d, &b, p);
            }
            for (int e = 0; e < d.n_em; ++e) px_em_start(&d, &b, e);
        }
        const auto [w0, w1] = work_slots(from);
        int cur = from;
        for (int s = 0; s < d.n_sub; ++s) {
            const int next = (s % 2 == 0) ? w0 : w1;
            px_buffers b = buffers(cur, next);
            if (s > 0) for (int c = 0; c < d.n_cells; ++c) px_thermo(&d, &b, c);
            d.sub = s;
            if (d.flow) flow_pass(d, b);
            for (int e = 0; e < d.n_el; ++e) px_surface_flux(&d, &b, e);
            for (int c = 0; c < d.n_cells; ++c) px_step_cell(&d, &b, c);
            slot_[next].flags = slot_[cur].flags;
            cur = next;
        }
        finish(cur, d, false);
        collect(r);
        r.trial_slot = cur;
    }
    void accumulate(const std::vector<double>& power, const std::vector<double>& JB, double dt) override {
        for (int e = 0; e < base_.n_em; ++e) px_accumulate_one(e, power.data(), JB.empty() ? nullptr : JB.data(), dt, acc_E_.data(), acc_F_.data(), acc_bad_.data());
    }
    void reset_accumulator() override {
        std::fill(acc_E_.begin(), acc_E_.end(), px_zero()); std::fill(acc_F_.begin(), acc_F_.end(), 0.0); std::fill(acc_bad_.begin(), acc_bad_.end(), 0);
    }
    std::vector<double> temperatures(int slot, const px_dims& d) override {
        std::vector<int> keep = status_;
        px_buffers b = buffers(slot, -1);
        for (int c = 0; c < d.n_cells; ++c) px_thermo(&d, &b, c);
        status_ = keep;
        return T_;
    }
    std::vector<double> face_displacements(int slot, const px_dims& d) override {
        std::vector<int> keep = status_;
        px_buffers b = buffers(slot, -1);
        for (int c = 0; c < d.n_cells; ++c) px_thermo(&d, &b, c);
        for (int fc = 0; fc < d.n_fc; ++fc) px_fc_prepare(&d, &b, fc);
        status_ = keep;
        return std::vector<double>(fc_eta_.begin(), fc_eta_.begin() + d.n_fc);
    }
    PxStateRoots state_roots(int slot, const px_dims& d, std::vector<std::array<px_u32, 8>>* cell_leaves) override {
        const auto& st = slot_[slot];
        const px_digest_src g{st.E.data(), st.N.data(), st.flags.data(), st.fcm.data(), st.fcV.data(), st.fcK.data()};
        return px_host_state_roots(g, d.n_cells, d.n_fc, d.flow && d.flow_model == 2, cell_leaves);
    }
};

}  // namespace ss304_phase_v1_2
