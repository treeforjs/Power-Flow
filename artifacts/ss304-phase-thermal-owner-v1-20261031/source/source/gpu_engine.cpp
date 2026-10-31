// gpu_engine.cpp -- see gpu_engine.hpp.
#include "gpu_engine.hpp"
#include "common/hall_gpu_common.hpp"
#include <chrono>
#include <cstring>
#include <sstream>

namespace ss304_phase_v1 {
namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kKernels[] = {
#embed "phase_kernels.h"
    , 0};
const char kLaw[] = {
#embed "phase_law.h"
    , 0};
const char kMath[] = {
#embed "det_math.h"
    , 0};
#pragma GCC diagnostic pop
const char kMain[] =
    "#include \"phase_kernels.h\"\n"
    "extern \"C\" __global__ void k_copy_i32(const int* a, int* b, int n) {\n"
    "    const int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) b[i] = a[i];\n"
    "}\n";

hall_gpu::SessionOptions session_options(const GpuOptions& o) {
    hall_gpu::SessionOptions s;
    s.require_uuid = o.require_uuid;
    s.nvrtc_library = o.nvrtc_library;
    return s;
}
using hall_gpu::DeviceBuffer;
}  // namespace

struct GpuEngine::Impl {
    hall_gpu::Session s;
    hall_gpu::Module mod;
    hall_gpu::Kernel kthermo, kscatter, ksigma, kjoule, kstep, kslice, kquote, krelease, kpatch, kcopy;
    unsigned block;
    pt_dims base{};
    int n_el = 0, P = 0;
    std::size_t sc_cap = 0, el_cap = 0;
    DeviceBuffer law, row_dy, row_em, em_row_begin, gx, gy, vref, mref;
    DeviceBuffer E[3], N[3], flags[3];
    DeviceBuffer T, k, sig, m, status, joule, em_energy;
    DeviceBuffer el_of, el_cell, el_face, el_area, el_patch, patch_begin;
    DeviceBuffer in_load, in_dep_N, in_dep_E, in_rel_N, in_rel_E, acc_N, acc_Eref, acc_Ekin, acc_bnd, quote;
    DeviceBuffer rel_total, rel_kin, patch_out, patch_sum, slice_sum, slice_status, slice_maxT, sigma_em, sc_el, sc_val;
    hall_gpu::HostBuffer stage;

    explicit Impl(const GpuOptions& o) : s(session_options(o)), block(o.block) {
        hall_gpu::CompileRequest req;
        req.main = {"phase_owner_kernels.cu", kMain};
        req.headers = {{"phase_kernels.h", kKernels}, {"phase_law.h", kLaw}, {"det_math.h", kMath}};
        mod = s.compile_and_load(req);
        kthermo = s.kernel(mod, "k_thermo"); kscatter = s.kernel(mod, "k_scatter"); ksigma = s.kernel(mod, "k_sigma");
        kjoule = s.kernel(mod, "k_joule"); kstep = s.kernel(mod, "k_step"); kslice = s.kernel(mod, "k_slice");
        kquote = s.kernel(mod, "k_quote"); krelease = s.kernel(mod, "k_release"); kpatch = s.kernel(mod, "k_patch");
        kcopy = s.kernel(mod, "k_copy_i32");
    }
    void ensure(DeviceBuffer& b, std::size_t bytes, const char* tag) {
        if (bytes == 0) bytes = 16;
        if (b && b.bytes >= bytes) return;
        if (b) s.release(b);
        b = s.alloc(bytes, tag);
    }
    void ensure_stage(std::size_t bytes) {
        if (stage && stage.bytes >= bytes) return;
        s.synchronize("stage");
        if (stage) s.release(stage);
        stage = s.alloc_pinned(std::max<std::size_t>(bytes, 1 << 20), "stage");
    }
    void up(DeviceBuffer& d, const void* src, std::size_t bytes) {
        if (!bytes) return;
        ensure_stage(bytes);
        s.synchronize("upload staging");
        std::memcpy(stage.pointer, src, bytes);
        s.upload(d, 0, stage, 0, bytes);
        s.synchronize("upload");
    }
    void down(void* dst, const DeviceBuffer& d, std::size_t bytes) {
        if (!bytes) return;
        ensure_stage(bytes);
        s.download(stage, 0, d, 0, bytes);
        s.synchronize("download");
        std::memcpy(dst, stage.pointer, bytes);
    }
    template <class T> void upv(DeviceBuffer& d, const std::vector<T>& v, const char* tag) { ensure(d, v.size() * sizeof(T), tag); up(d, v.data(), v.size() * sizeof(T)); }
    template <class T> std::vector<T> downv(const DeviceBuffer& d, std::size_t n) { std::vector<T> v(n); down(v.data(), d, n * sizeof(T)); return v; }
    template <class T> static T* ptr(const DeviceBuffer& b) { return reinterpret_cast<T*>(b.address); }

    pt_buffers buffers(int cur, int next) const {
        pt_buffers b{};
        b.law = ptr<const pt_law>(law);
        b.row_dy = ptr<const double>(row_dy); b.row_em = ptr<const int>(row_em); b.em_row_begin = ptr<const int>(em_row_begin);
        b.gx = ptr<const double>(gx); b.gy = ptr<const double>(gy); b.vref = ptr<const double>(vref); b.mref = ptr<const double>(mref);
        b.E = ptr<const pt_i128>(E[cur]); b.N = ptr<const pt_i128>(N[cur]); b.flags = ptr<const int>(flags[cur]);
        if (next >= 0) { b.E_out = ptr<pt_i128>(E[next]); b.N_out = ptr<pt_i128>(N[next]); }
        b.T = ptr<double>(T); b.k = ptr<double>(k); b.sig = ptr<double>(sig); b.m = ptr<double>(m); b.status = ptr<int>(status);
        b.joule_total = ptr<pt_i128>(joule); b.em_energy = ptr<const double>(em_energy);
        b.el_of = ptr<const int>(el_of); b.el_cell = ptr<const int>(el_cell); b.el_face = ptr<const int>(el_face);
        b.el_area = ptr<const double>(el_area); b.el_patch = ptr<const int>(el_patch); b.patch_begin = ptr<const int>(patch_begin);
        b.in_load_E = ptr<const pt_i128>(in_load); b.in_dep_N = ptr<const pt_i128>(in_dep_N); b.in_dep_E = ptr<const pt_i128>(in_dep_E);
        b.in_rel_N = ptr<pt_i128>(in_rel_N); b.in_rel_E = ptr<pt_i128>(in_rel_E);
        b.acc_ev_N = ptr<pt_i128>(acc_N); b.acc_ev_Eref = ptr<pt_i128>(acc_Eref); b.acc_ev_Ekin = ptr<pt_i128>(acc_Ekin); b.acc_bnd_E = ptr<pt_i128>(acc_bnd);
        b.el_quote = ptr<double>(quote);
        b.rel_total_N = ptr<const pt_i128>(rel_total); b.rel_kin_J = ptr<const double>(rel_kin);
        b.patch_out = ptr<double>(patch_out); b.patch_sum = ptr<pt_i128>(patch_sum);
        b.slice_sum = ptr<pt_i128>(slice_sum); b.slice_status = ptr<int>(slice_status); b.slice_maxT = ptr<double>(slice_maxT);
        b.sigma_em = ptr<double>(sigma_em);
        b.sc_el = ptr<const int>(sc_el); b.sc_val = ptr<const pt_i128>(sc_val);
        b.in_load_E_w = ptr<pt_i128>(in_load); b.in_dep_N_w = ptr<pt_i128>(in_dep_N); b.in_dep_E_w = ptr<pt_i128>(in_dep_E);
        return b;
    }
    void launch(const hall_gpu::Kernel& kk, int n, pt_dims d, pt_buffers b) {
        if (n <= 0) return;
        void* args[] = {&d, &b};
        s.launch(kk, {unsigned((n + int(block) - 1) / int(block)), 1, 1}, {block, 1, 1}, 0, args);
    }
    void zero_elements() {
        const std::size_t n = std::size_t(std::max(n_el, 1));
        for (DeviceBuffer* b : {&in_load, &in_dep_E, &in_rel_E, &acc_Eref, &acc_Ekin, &acc_bnd}) s.zero(*b, 0, n * sizeof(pt_i128));
        for (DeviceBuffer* b : {&in_dep_N, &in_rel_N, &acc_N}) s.zero(*b, 0, 3 * n * sizeof(pt_i128));
        s.zero(quote, 0, 3 * n * sizeof(double));
        s.zero(patch_out, 0, std::size_t(P) * PT_PATCH_OUT * sizeof(double));
        s.zero(patch_sum, 0, std::size_t(P) * PT_PATCH_SUM * sizeof(pt_i128));
    }
    void finish(int slot, const pt_dims& d) {
        const pt_buffers b = buffers(slot, -1);
        launch(kthermo, d.n_cells, d, b);
        launch(kslice, d.S, d, b);
        launch(kpatch, d.n_patches, d, b);
    }
    void collect(const pt_dims& d, EngineResult& r) {
        r.slice_sum = downv<pt_i128>(slice_sum, std::size_t(d.S) * PT_SLICE_SUM);
        r.slice_status = downv<int>(slice_status, std::size_t(d.S) * 4);
        r.slice_maxT = downv<double>(slice_maxT, d.S);
        r.patch_sum = downv<pt_i128>(patch_sum, std::size_t(d.n_patches) * PT_PATCH_SUM);
        r.patch_out = downv<double>(patch_out, std::size_t(d.n_patches) * PT_PATCH_OUT);
    }
};

GpuEngine::GpuEngine(const GpuOptions& o) : p_(new Impl(o)) {}
GpuEngine::~GpuEngine() { delete p_; }

void GpuEngine::configure(const pt_law& law, const pt_dims& dims, const HostGeometry& g) {
    auto& I = *p_;
    I.base = dims; I.P = dims.n_patches;
    I.ensure(I.law, sizeof(pt_law), "law"); I.up(I.law, &law, sizeof(pt_law));
    I.upv(I.row_dy, g.row_dy, "row_dy"); I.upv(I.row_em, g.row_em, "row_em"); I.upv(I.em_row_begin, g.em_row_begin, "em_row_begin");
    I.upv(I.gx, g.gx, "gx"); I.upv(I.gy, g.gy, "gy"); I.upv(I.vref, g.vref, "vref"); I.upv(I.mref, g.mref, "mref");
    const std::size_t n = dims.n_cells;
    for (int s = 0; s < 3; ++s) {
        I.ensure(I.E[s], n * sizeof(pt_i128), "E"); I.ensure(I.N[s], 3 * n * sizeof(pt_i128), "N"); I.ensure(I.flags[s], n * sizeof(int), "flags");
        I.s.zero(I.E[s], 0, n * sizeof(pt_i128)); I.s.zero(I.N[s], 0, 3 * n * sizeof(pt_i128)); I.s.zero(I.flags[s], 0, n * sizeof(int));
    }
    for (DeviceBuffer* b : {&I.T, &I.k, &I.sig, &I.m}) { I.ensure(*b, n * sizeof(double), "scratch"); I.s.zero(*b, 0, n * sizeof(double)); }
    I.ensure(I.status, n * sizeof(int), "status"); I.s.zero(I.status, 0, n * sizeof(int));
    I.ensure(I.joule, n * sizeof(pt_i128), "joule"); I.s.zero(I.joule, 0, n * sizeof(pt_i128));
    I.ensure(I.em_energy, std::size_t(dims.n_em) * sizeof(double), "em_energy"); I.s.zero(I.em_energy, 0, std::size_t(dims.n_em) * sizeof(double));
    I.ensure(I.sigma_em, std::size_t(dims.n_em) * sizeof(double), "sigma_em");
    I.ensure(I.slice_sum, std::size_t(dims.S) * PT_SLICE_SUM * sizeof(pt_i128), "slice_sum");
    I.ensure(I.slice_status, std::size_t(dims.S) * 4 * sizeof(int), "slice_status");
    I.ensure(I.slice_maxT, std::size_t(dims.S) * sizeof(double), "slice_maxT");
    const std::size_t P = std::max(I.P, 1);
    I.ensure(I.rel_total, 3 * P * sizeof(pt_i128), "rel_total"); I.ensure(I.rel_kin, 3 * P * sizeof(double), "rel_kin");
    I.ensure(I.patch_out, P * PT_PATCH_OUT * sizeof(double), "patch_out"); I.ensure(I.patch_sum, P * PT_PATCH_SUM * sizeof(pt_i128), "patch_sum");
    I.ensure(I.sc_el, 64 * sizeof(int), "sc_el"); I.ensure(I.sc_val, 64 * 5 * sizeof(pt_i128), "sc_val");
    I.s.synchronize("configure");
}

void GpuEngine::set_elements(const HostElements& el) {
    auto& I = *p_;
    I.n_el = int(el.el_cell.size());
    I.upv(I.el_of, el.el_of, "el_of"); I.upv(I.el_cell, el.el_cell, "el_cell"); I.upv(I.el_face, el.el_face, "el_face");
    I.upv(I.el_area, el.el_area, "el_area"); I.upv(I.el_patch, el.el_patch, "el_patch"); I.upv(I.patch_begin, el.patch_begin, "patch_begin");
    const std::size_t n = std::size_t(std::max(I.n_el, 1));
    for (DeviceBuffer* b : {&I.in_load, &I.in_dep_E, &I.in_rel_E, &I.acc_Eref, &I.acc_Ekin, &I.acc_bnd}) I.ensure(*b, n * sizeof(pt_i128), "element");
    for (DeviceBuffer* b : {&I.in_dep_N, &I.in_rel_N, &I.acc_N}) I.ensure(*b, 3 * n * sizeof(pt_i128), "element3");
    I.ensure(I.quote, 3 * n * sizeof(double), "quote");
    I.zero_elements();
    const std::size_t P = std::max(I.P, 1);
    I.s.zero(I.rel_total, 0, 3 * P * sizeof(pt_i128)); I.s.zero(I.rel_kin, 0, 3 * P * sizeof(double));
    I.s.synchronize("set_elements");
}

void GpuEngine::write_state(int slot, const HostState& st) {
    auto& I = *p_;
    I.up(I.E[slot], st.E.data(), st.E.size() * sizeof(pt_i128));
    I.up(I.N[slot], st.N.data(), st.N.size() * sizeof(pt_i128));
    I.up(I.flags[slot], st.flags.data(), st.flags.size() * sizeof(int));
}
HostState GpuEngine::read_state(int slot) {
    auto& I = *p_;
    const std::size_t n = I.base.n_cells;
    HostState st;
    st.E = I.downv<pt_i128>(I.E[slot], n); st.N = I.downv<pt_i128>(I.N[slot], 3 * n); st.flags = I.downv<int>(I.flags[slot], n);
    return st;
}
void GpuEngine::evaluate(int slot, const pt_dims& d, EngineResult& r, std::vector<double>* sigma_em) {
    auto& I = *p_;
    I.s.zero(I.status, 0, std::size_t(d.n_cells) * sizeof(int));
    I.finish(slot, d);
    if (sigma_em) { I.launch(I.ksigma, d.n_em, d, I.buffers(slot, -1)); *sigma_em = I.downv<double>(I.sigma_em, d.n_em); }
    I.collect(d, r);
}
void GpuEngine::step(int from, const pt_dims& d0, const StepData& in, EngineResult& r) {
    auto& I = *p_;
    const auto t0 = std::chrono::steady_clock::now();
    pt_dims d = d0;
    I.s.zero(I.status, 0, std::size_t(d.n_cells) * sizeof(int));
    I.zero_elements();
    I.up(I.em_energy, in.em_energy.data(), in.em_energy.size() * sizeof(double));
    I.up(I.rel_total, in.rel_total_N.data(), in.rel_total_N.size() * sizeof(pt_i128));
    I.up(I.rel_kin, in.rel_kin_J.data(), in.rel_kin_J.size() * sizeof(double));
    if (!in.sc_el.empty()) { I.upv(I.sc_el, in.sc_el, "sc_el"); I.upv(I.sc_val, in.sc_val, "sc_val"); }
    {
        const pt_buffers b = I.buffers(from, -1);
        I.launch(I.kthermo, d.n_cells, d, b);
        pt_dims ds = d; ds.n_scatter = int(in.sc_el.size());
        I.launch(I.kscatter, ds.n_scatter, ds, b);
        if (d.evap_mode == 2) { I.launch(I.kquote, d.n_el, d, b); I.launch(I.krelease, d.n_patches, d, b); }
        I.launch(I.kjoule, d.n_em, d, b);
    }
    const auto [w0, w1] = work_slots(from);
    int cur = from;
    for (int s = 0; s < d.n_sub; ++s) {
        const int next = (s % 2 == 0) ? w0 : w1;
        const pt_buffers b = I.buffers(cur, next);
        if (s > 0) I.launch(I.kthermo, d.n_cells, d, b);
        d.sub = s;
        I.launch(I.kstep, d.n_cells, d, b);
        {
            const int* a = Impl::ptr<const int>(I.flags[cur]); int* o = Impl::ptr<int>(I.flags[next]); int n = d.n_cells;
            void* args[] = {&a, &o, &n};
            I.s.launch(I.kcopy, {unsigned((n + int(I.block) - 1) / int(I.block)), 1, 1}, {I.block, 1, 1}, 0, args);
        }
        cur = next;
    }
    I.finish(cur, d);
    I.collect(d, r);
    r.trial_slot = cur;
    last_step_ms_ = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
std::vector<double> GpuEngine::temperatures(int slot, const pt_dims& d) {
    auto& I = *p_;
    I.launch(I.kthermo, d.n_cells, d, I.buffers(slot, -1));
    return I.downv<double>(I.T, d.n_cells);
}
std::string GpuEngine::record_json() const {
    const auto& I = *p_;
    std::ostringstream o;
    o << "{\"device\":" << hall_gpu::device_json(I.s.device()) << ",\"compile\":" << hall_gpu::compile_json(I.mod.record) << "}";
    return o.str();
}

}  // namespace ss304_phase_v1
