// engine.hpp -- execution backends of the phase owner. The CPU reference runs the phase_kernels.h functions in
// index loops; the GPU engine (gpu_engine.hpp) runs the same functions as NVRTC kernels. Both execute the same
// sequence, so the trial state words are bit-identical.
#pragma once
#include "phase_law_host.hpp"
#include <algorithm>
#include <memory>
#include <vector>

namespace ss304_phase_v1 {

struct HostState {
    std::vector<pt_i128> E, N;   // [n_cells], [3 n_cells]
    std::vector<int> flags;      // bit0: void (merged away)
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
    std::vector<int> sc_el;                    // distinct elements with loads/deposits
    std::vector<pt_i128> sc_val;               // [5 per entry]: load_E, dep_N x3, dep_E
    std::vector<pt_i128> rel_total_N;          // [3 n_patches]
    std::vector<double> rel_kin_J;             // [3 n_patches]
};
struct EngineResult {
    int trial_slot = -1;
    std::vector<pt_i128> slice_sum;            // [S*PT_SLICE_SUM]
    std::vector<int> slice_status;             // [S*4]
    std::vector<double> slice_maxT;            // [S]
    std::vector<pt_i128> patch_sum;            // [P*PT_PATCH_SUM]
    std::vector<double> patch_out;             // [P*PT_PATCH_OUT]
};

class Engine {
public:
    virtual ~Engine() = default;
    virtual std::string name() const = 0;
    virtual void configure(const pt_law& law, const pt_dims& dims, const HostGeometry& g) = 0;
    virtual void set_elements(const HostElements& el) = 0;
    virtual void write_state(int slot, const HostState& s) = 0;
    virtual HostState read_state(int slot) = 0;
    // Thermodynamic evaluation of a slot: slice sums/statuses, patch report (quotes over d.dt), sigma per EM cell.
    virtual void evaluate(int slot, const pt_dims& d, EngineResult& r, std::vector<double>* sigma_em) = 0;
    // One owner step from slot `from` into one of the two other slots (r.trial_slot).
    virtual void step(int from, const pt_dims& d, const StepData& in, EngineResult& r) = 0;
    // Surface temperatures of the given cells at the last evaluated slot (diagnostics and tests).
    virtual std::vector<double> temperatures(int slot, const pt_dims& d) = 0;
};

inline std::pair<int, int> work_slots(int from) {
    int a = -1, b = -1;
    for (int s = 0; s < 3; ++s) if (s != from) { if (a < 0) a = s; else b = s; }
    return {a, b};
}

// ---------------------------------------------------------------------------------------------------------
class CpuEngine final : public Engine {
    pt_law law_{};
    pt_dims base_{};
    HostGeometry g_;
    HostElements el_;
    HostState slot_[3];
    std::vector<double> T_, k_, sig_, m_, quote_, patch_out_, slice_maxT_, sigma_em_, em_energy_, rel_kin_;
    std::vector<int> status_, slice_status_;
    std::vector<pt_i128> joule_, in_load_, in_dep_N_, in_dep_E_, in_rel_N_, in_rel_E_, acc_ev_N_, acc_ev_Eref_, acc_ev_Ekin_, acc_bnd_,
        patch_sum_, slice_sum_, rel_total_;

    pt_buffers buffers(int cur, int next) {
        pt_buffers b{};
        b.law = &law_;
        b.row_dy = g_.row_dy.data(); b.row_em = g_.row_em.data(); b.em_row_begin = g_.em_row_begin.data();
        b.gx = g_.gx.data(); b.gy = g_.gy.data(); b.vref = g_.vref.data(); b.mref = g_.mref.data();
        b.E = slot_[cur].E.data(); b.N = slot_[cur].N.data(); b.flags = slot_[cur].flags.data();
        if (next >= 0) { b.E_out = slot_[next].E.data(); b.N_out = slot_[next].N.data(); }
        b.T = T_.data(); b.k = k_.data(); b.sig = sig_.data(); b.m = m_.data(); b.status = status_.data();
        b.joule_total = joule_.data(); b.em_energy = em_energy_.data();
        b.el_of = el_.el_of.data(); b.el_cell = el_.el_cell.data(); b.el_face = el_.el_face.data(); b.el_area = el_.el_area.data();
        b.el_patch = el_.el_patch.data(); b.patch_begin = el_.patch_begin.data();
        b.in_load_E = in_load_.data(); b.in_dep_N = in_dep_N_.data(); b.in_dep_E = in_dep_E_.data();
        b.in_rel_N = in_rel_N_.data(); b.in_rel_E = in_rel_E_.data();
        b.acc_ev_N = acc_ev_N_.data(); b.acc_ev_Eref = acc_ev_Eref_.data(); b.acc_ev_Ekin = acc_ev_Ekin_.data(); b.acc_bnd_E = acc_bnd_.data();
        b.el_quote = quote_.data();
        b.rel_total_N = rel_total_.data(); b.rel_kin_J = rel_kin_.data();
        b.patch_out = patch_out_.data(); b.patch_sum = patch_sum_.data(); b.slice_sum = slice_sum_.data(); b.slice_status = slice_status_.data();
        b.slice_maxT = slice_maxT_.data(); b.sigma_em = sigma_em_.data();
        b.in_load_E_w = in_load_.data(); b.in_dep_N_w = in_dep_N_.data(); b.in_dep_E_w = in_dep_E_.data();
        return b;
    }
    void sized_elements() {
        const std::size_t n = el_.el_cell.size(), P = el_.patch_begin.size() - 1;
        for (auto* v : {&in_load_, &in_dep_E_, &in_rel_E_, &acc_ev_Eref_, &acc_ev_Ekin_, &acc_bnd_}) v->assign(n, pt_zero());
        for (auto* v : {&in_dep_N_, &in_rel_N_, &acc_ev_N_}) v->assign(3 * n, pt_zero());
        quote_.assign(3 * n, 0.0);
        patch_out_.assign(PT_PATCH_OUT * P, 0.0);
        patch_sum_.assign(PT_PATCH_SUM * P, pt_zero());
        rel_total_.assign(3 * P, pt_zero());
        rel_kin_.assign(3 * P, 0.0);
    }
    void collect(EngineResult& r) {
        r.slice_sum = slice_sum_; r.slice_status = slice_status_; r.slice_maxT = slice_maxT_;
        r.patch_sum = patch_sum_; r.patch_out = patch_out_;
    }
    void finish(int slot, const pt_dims& d) {
        pt_buffers b = buffers(slot, -1);
        for (int c = 0; c < d.n_cells; ++c) pt_thermo(&d, &b, c);
        for (int s = 0; s < d.S; ++s) pt_slice_reduce(&d, &b, s);
        for (int p = 0; p < d.n_patches; ++p) pt_patch_report(&d, &b, p);
    }

public:
    std::string name() const override { return "cpu-reference"; }
    void configure(const pt_law& law, const pt_dims& dims, const HostGeometry& g) override {
        law_ = law; base_ = dims; g_ = g;
        const std::size_t n = dims.n_cells;
        for (auto& s : slot_) { s.E.assign(n, pt_zero()); s.N.assign(3 * n, pt_zero()); s.flags.assign(n, 0); }
        T_.assign(n, 0); k_.assign(n, 0); sig_.assign(n, 0); m_.assign(n, 0); status_.assign(n, 0); joule_.assign(n, pt_zero());
        slice_sum_.assign(PT_SLICE_SUM * dims.S, pt_zero()); slice_status_.assign(4 * dims.S, 0); slice_maxT_.assign(dims.S, 0);
        sigma_em_.assign(dims.n_em, 0); em_energy_.assign(dims.n_em, 0);
    }
    void set_elements(const HostElements& el) override { el_ = el; sized_elements(); }
    void write_state(int slot, const HostState& s) override { slot_[slot] = s; }
    HostState read_state(int slot) override { return slot_[slot]; }
    void evaluate(int slot, const pt_dims& d, EngineResult& r, std::vector<double>* sigma_em) override {
        std::fill(status_.begin(), status_.end(), 0);
        finish(slot, d);
        if (sigma_em) {
            pt_buffers b = buffers(slot, -1);
            for (int e = 0; e < d.n_em; ++e) pt_sigma_em(&d, &b, e);
            *sigma_em = sigma_em_;
        }
        collect(r);
    }
    void step(int from, const pt_dims& d0, const StepData& in, EngineResult& r) override {
        pt_dims d = d0;
        std::fill(status_.begin(), status_.end(), 0);
        sized_elements();
        em_energy_ = in.em_energy;
        rel_total_ = in.rel_total_N; rel_kin_ = in.rel_kin_J;
        {
            pt_buffers b = buffers(from, -1);
            for (int c = 0; c < d.n_cells; ++c) pt_thermo(&d, &b, c);
            pt_dims ds = d; ds.n_scatter = int(in.sc_el.size());
            b.sc_el = in.sc_el.data(); b.sc_val = in.sc_val.data();
            for (int q = 0; q < ds.n_scatter; ++q) pt_scatter(&ds, &b, q);
            if (d.evap_mode == 2) {
                for (int e = 0; e < d.n_el; ++e) pt_quote_element(&d, &b, e);
                for (int p = 0; p < d.n_patches; ++p) pt_release_split(&d, &b, p);
            }
            for (int e = 0; e < d.n_em; ++e) pt_joule_split(&d, &b, e);
        }
        const auto [w0, w1] = work_slots(from);
        int cur = from;
        for (int s = 0; s < d.n_sub; ++s) {
            const int next = (s % 2 == 0) ? w0 : w1;
            pt_buffers b = buffers(cur, next);
            if (s > 0) for (int c = 0; c < d.n_cells; ++c) pt_thermo(&d, &b, c);
            d.sub = s;
            for (int c = 0; c < d.n_cells; ++c) pt_step_cell(&d, &b, c);
            slot_[next].flags = slot_[cur].flags;
            cur = next;
        }
        finish(cur, d);
        collect(r);
        r.trial_slot = cur;
    }
    std::vector<double> temperatures(int slot, const pt_dims& d) override {
        std::vector<int> keep = status_;
        pt_buffers b = buffers(slot, -1);
        for (int c = 0; c < d.n_cells; ++c) pt_thermo(&d, &b, c);
        status_ = keep;
        return T_;
    }
};

}  // namespace ss304_phase_v1
