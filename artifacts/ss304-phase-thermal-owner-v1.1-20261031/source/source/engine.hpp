// engine.hpp -- execution backends of the phase owner v1.1. The CPU reference runs the phase_kernels.h functions
// in index loops; the GPU engine (gpu_engine.hpp) runs the same functions as NVRTC kernels in the same sequence,
// so the trial state words are bit-identical.
#pragma once
#include "phase_law_host.hpp"
#include "state_digest.h"
#include <algorithm>
#include <memory>
#include <vector>

namespace ss304_phase_v1_1 {

struct HostState {
    std::vector<pq_i128> E, N;   // [n_cells], [PQ_NEL n_cells]
    std::vector<int> flags;      // bit0: void (merged away)
    std::vector<pq_i128> fcm;    // [n_fc PQ_NEL]: atoms gained through each large-face column (displacement)
};
struct HostGeometry {
    std::vector<double> row_dy, gx, gy, vref, mref;
    std::vector<int> row_em, em_row_begin;
};
struct HostElements {
    std::vector<int> el_of, el_cell, el_face, el_patch, patch_begin;
    std::vector<double> el_area;
    bool operator==(const HostElements&) const = default;
};
struct StepData {
    std::vector<double> em_energy;            // [n_em] J
    std::vector<double> em_JB;                // [6 n_em] J (A/m2) and B (T) per EM cell
    std::vector<double> patch_pb;             // [n_patches] back pressure (Pa)
    std::vector<int> sc_el;                   // distinct elements with loads/deposits
    std::vector<pq_i128> sc_val;              // [PQ_SC per entry]: load_E, dep_N x6, dep_E
    std::vector<pq_i128> rel_total_N;         // [PQ_NEL n_patches]
    std::vector<double> rel_kin_J;            // [PQ_NEL n_patches]
};
struct EngineResult {
    int trial_slot = -1;
    std::vector<pq_i128> slice_sum;           // [S*PQ_SLICE_SUM]
    std::vector<int> slice_status;            // [S*4]
    std::vector<double> slice_maxT;           // [S]
    std::vector<pq_i128> patch_sum;           // [P*PQ_PATCH_SUM]
    std::vector<double> patch_out;            // [P*PQ_PATCH_OUT]
};

class Engine {
public:
    virtual ~Engine() = default;
    virtual std::string name() const = 0;
    virtual void configure(const pq_law& law, const pq_dims& dims, const HostGeometry& g) = 0;
    virtual void set_elements(const HostElements& el) = 0;
    virtual void write_state(int slot, const HostState& s) = 0;
    virtual HostState read_state(int slot) = 0;
    virtual void evaluate(int slot, const pq_dims& d, const std::vector<double>& patch_pb, EngineResult& r, std::vector<double>* sigma_em) = 0;
    virtual void step(int from, const pq_dims& d, const StepData& in, EngineResult& r) = 0;
    virtual std::vector<double> temperatures(int slot, const pq_dims& d) = 0;
    virtual std::vector<double> face_displacements(int slot, const pq_dims& d) = 0;   // fc_eta of the slot
    // Block-Merkle state digest of a slot (state_digest.h); leaves optionally returned for host sampling.
    virtual StateRoots state_roots(int slot, const pq_dims& d, std::vector<std::array<pq_u32, 8>>* cell_leaves) = 0;
};

inline std::pair<int, int> work_slots(int from) {
    int a = -1, b = -1;
    for (int s = 0; s < 3; ++s) if (s != from) { if (a < 0) a = s; else b = s; }
    return {a, b};
}

// ---------------------------------------------------------------------------------------------------------
class CpuEngine final : public Engine {
    pq_law law_{};
    pq_dims base_{};
    HostGeometry g_;
    HostElements el_;
    HostState slot_[3];
    std::vector<double> T_, k_, sig_, m_, quote_, patch_out_, slice_maxT_, sigma_em_, em_energy_, em_JB_, patch_pb_, rel_kin_;
    std::vector<double> fc_delta_, fc_T_, fc_eta_, fc_p_, fc_x_, ic_q_;
    std::vector<int> status_, slice_status_, fc_cell_;
    std::vector<pq_i128> joule_, in_load_, in_dep_N_, in_dep_E_, in_rel_N_, in_rel_E_, acc_ev_N_, acc_ev_Eref_, acc_ev_Ekin_,
        acc_ev_Eint_, acc_bnd_, acc_rad_, patch_sum_, slice_sum_, rel_total_, ic_N_, ic_E_;

    pq_buffers buffers(int cur, int next) {
        pq_buffers b{};
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
        return b;
    }
    void sized_elements() {
        const std::size_t n = el_.el_cell.size(), P = el_.patch_begin.size() - 1;
        for (auto* v : {&in_load_, &in_dep_E_, &in_rel_E_, &acc_ev_Eref_, &acc_ev_Ekin_, &acc_ev_Eint_, &acc_bnd_, &acc_rad_}) v->assign(std::max<std::size_t>(n, 1), pq_zero());
        for (auto* v : {&in_dep_N_, &in_rel_N_, &acc_ev_N_}) v->assign(PQ_NEL * std::max<std::size_t>(n, 1), pq_zero());
        quote_.assign(PQ_NEL * std::max<std::size_t>(n, 1), 0.0);
        patch_out_.assign(PQ_PATCH_OUT * P, 0.0);
        patch_sum_.assign(PQ_PATCH_SUM * P, pq_zero());
        rel_total_.assign(PQ_NEL * P, pq_zero());
        rel_kin_.assign(PQ_NEL * P, 0.0);
    }
    void collect(EngineResult& r) {
        r.slice_sum = slice_sum_; r.slice_status = slice_status_; r.slice_maxT = slice_maxT_;
        r.patch_sum = patch_sum_; r.patch_out = patch_out_;
    }
    void flow_pass(const pq_dims& d, const pq_buffers& b) {
        for (int fc = 0; fc < d.n_fc; ++fc) pq_fc_prepare(&d, &b, fc);
        for (int fc = 0; fc < d.n_fc; ++fc) pq_fc_pressure(&d, &b, fc);
        for (int ic = 0; ic < d.n_ic; ++ic) pq_ic_flux(&d, &b, ic);
    }
    void finish(int slot, const pq_dims& d) {
        pq_buffers b = buffers(slot, -1);
        for (int c = 0; c < d.n_cells; ++c) pq_thermo(&d, &b, c);
        for (int fc = 0; fc < d.n_fc; ++fc) pq_fc_prepare(&d, &b, fc);
        for (int s = 0; s < d.S; ++s) pq_slice_reduce(&d, &b, s);
        for (int p = 0; p < d.n_patches; ++p) pq_patch_report(&d, &b, p);
    }

public:
    std::string name() const override { return "cpu-reference"; }
    void configure(const pq_law& law, const pq_dims& dims, const HostGeometry& g) override {
        law_ = law; base_ = dims; g_ = g;
        const std::size_t n = dims.n_cells, nfc = std::max(dims.n_fc, 1), nic = std::max(dims.n_ic, 1);
        for (auto& s : slot_) { s.E.assign(n, pq_zero()); s.N.assign(PQ_NEL * n, pq_zero()); s.flags.assign(n, 0); s.fcm.assign(PQ_NEL * nfc, pq_zero()); }
        T_.assign(n, 0); k_.assign(n, 0); sig_.assign(n, 0); m_.assign(n, 0); status_.assign(n, 0); joule_.assign(n, pq_zero());
        slice_sum_.assign(PQ_SLICE_SUM * dims.S, pq_zero()); slice_status_.assign(4 * dims.S, 0); slice_maxT_.assign(dims.S, 0);
        sigma_em_.assign(dims.n_em, 0); em_energy_.assign(dims.n_em, 0); em_JB_.assign(6 * dims.n_em, 0.0);
        patch_pb_.assign(std::max(dims.n_patches, 1), 0.0);
        fc_cell_.assign(nfc, -1); fc_delta_.assign(nfc, 0); fc_T_.assign(nfc, 0); fc_eta_.assign(nfc, 0); fc_p_.assign(nfc, 0); fc_x_.assign(PQ_NEL * nfc, 0);
        ic_N_.assign(PQ_NEL * nic, pq_zero()); ic_E_.assign(nic, pq_zero()); ic_q_.assign(nic, 0);
    }
    void set_elements(const HostElements& el) override { el_ = el; sized_elements(); }
    void write_state(int slot, const HostState& s) override { slot_[slot] = s; }
    HostState read_state(int slot) override { return slot_[slot]; }
    void evaluate(int slot, const pq_dims& d, const std::vector<double>& pb, EngineResult& r, std::vector<double>* sigma_em) override {
        std::fill(status_.begin(), status_.end(), 0);
        patch_pb_ = pb;
        finish(slot, d);
        if (sigma_em) {
            pq_buffers b = buffers(slot, -1);
            for (int e = 0; e < d.n_em; ++e) pq_sigma_em(&d, &b, e);
            *sigma_em = sigma_em_;
        }
        collect(r);
    }
    void step(int from, const pq_dims& d0, const StepData& in, EngineResult& r) override {
        pq_dims d = d0;
        std::fill(status_.begin(), status_.end(), 0);
        sized_elements();
        em_energy_ = in.em_energy; em_JB_ = in.em_JB; patch_pb_ = in.patch_pb;
        rel_total_ = in.rel_total_N; rel_kin_ = in.rel_kin_J;
        {
            pq_buffers b = buffers(from, -1);
            for (int c = 0; c < d.n_cells; ++c) pq_thermo(&d, &b, c);
            pq_dims ds = d; ds.n_scatter = int(in.sc_el.size());
            b.sc_el = in.sc_el.data(); b.sc_val = in.sc_val.data();
            for (int q = 0; q < ds.n_scatter; ++q) pq_scatter(&ds, &b, q);
            if (d.evap_mode == 2) {
                for (int e = 0; e < d.n_el; ++e) pq_quote_element(&d, &b, e);
                for (int p = 0; p < d.n_patches; ++p) pq_release_split(&d, &b, p);
            }
            for (int e = 0; e < d.n_em; ++e) pq_joule_split(&d, &b, e);
        }
        const auto [w0, w1] = work_slots(from);
        int cur = from;
        for (int s = 0; s < d.n_sub; ++s) {
            const int next = (s % 2 == 0) ? w0 : w1;
            pq_buffers b = buffers(cur, next);
            if (s > 0) for (int c = 0; c < d.n_cells; ++c) pq_thermo(&d, &b, c);
            d.sub = s;
            if (d.flow) flow_pass(d, b);
            for (int c = 0; c < d.n_cells; ++c) pq_step_cell(&d, &b, c);
            slot_[next].flags = slot_[cur].flags;
            cur = next;
        }
        finish(cur, d);
        collect(r);
        r.trial_slot = cur;
    }
    std::vector<double> temperatures(int slot, const pq_dims& d) override {
        std::vector<int> keep = status_;
        pq_buffers b = buffers(slot, -1);
        for (int c = 0; c < d.n_cells; ++c) pq_thermo(&d, &b, c);
        status_ = keep;
        return T_;
    }
    std::vector<double> face_displacements(int slot, const pq_dims& d) override {
        std::vector<int> keep = status_;
        pq_buffers b = buffers(slot, -1);
        for (int c = 0; c < d.n_cells; ++c) pq_thermo(&d, &b, c);
        for (int fc = 0; fc < d.n_fc; ++fc) pq_fc_prepare(&d, &b, fc);
        status_ = keep;
        return std::vector<double>(fc_eta_.begin(), fc_eta_.begin() + d.n_fc);
    }
    StateRoots state_roots(int slot, const pq_dims& d, std::vector<std::array<pq_u32, 8>>* cell_leaves) override {
        return host_state_roots(slot_[slot].E.data(), slot_[slot].N.data(), slot_[slot].flags.data(), d.n_cells, slot_[slot].fcm.data(), d.n_fc, cell_leaves);
    }
};

}  // namespace ss304_phase_v1_1
