// gpu_engine.cpp -- see gpu_engine.hpp.
#include "gpu_engine.hpp"
#include "common/hall_gpu_common.hpp"
#include <chrono>
#include <cstring>
#include <sstream>

namespace ss304_phase_v1_1 {
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
const char kDigest[] = {
#embed "state_digest.h"
    , 0};
#pragma GCC diagnostic pop
const char kMain[] =
    "#include \"phase_kernels.h\"\n"
    "#include \"state_digest.h\"\n"
    "extern \"C\" __global__ void k_copy_i32(const int* a, int* b, int n) {\n"
    "    const int i = blockIdx.x * blockDim.x + threadIdx.x; if (i < n) b[i] = a[i];\n"
    "}\n"
    "extern \"C\" __global__ void k_chunks(int section, pq_u64 n, const pq_i128* E, const pq_i128* N, const int* flags, const pq_i128* fcm, pq_u32* out) {\n"
    "    const pq_u64 j = (pq_u64)blockIdx.x * blockDim.x + threadIdx.x; if (j < (n + 31) / 32) pq_chunk_digest(section, j, n, E, N, flags, fcm, out + 8 * j);\n"
    "}\n"
    "extern \"C\" __global__ void k_leaves(int section, pq_u64 n, const pq_u32* chunks, pq_u32* out) {\n"
    "    const pq_u64 k = (pq_u64)blockIdx.x * blockDim.x + threadIdx.x; if (k < (n + 1023) / 1024) pq_leaf_digest(section, k, n, chunks, out + 8 * k);\n"
    "}\n";

using hall_gpu::DeviceBuffer;
}  // namespace

struct GpuEngine::Impl {
    std::unique_ptr<hall_gpu::Session> own;
    hall_gpu::Session* sp = nullptr;
    hall_gpu::Session& s() { return *sp; }
    hall_gpu::Module mod;
    hall_gpu::Kernel kthermo, kscatter, ksigma, kjoule, kfcp, kfcpr, kic, kstep, kslice, kquote, krelease, kpatch, kcopy, kchunks, kleaves;
    unsigned block;
    pq_dims base{};
    int n_el = 0, P = 0;
    DeviceBuffer law, row_dy, row_em, em_row_begin, gx, gy, vref, mref;
    DeviceBuffer E[3], N[3], flags[3], fcm[3];
    DeviceBuffer T, k, sig, m, status, joule, em_energy, em_JB, patch_pb;
    DeviceBuffer el_of, el_cell, el_face, el_area, el_patch, patch_begin;
    DeviceBuffer in_load, in_dep_N, in_dep_E, in_rel_N, in_rel_E, acc_N, acc_Eref, acc_Ekin, acc_Eint, acc_bnd, acc_rad, quote;
    DeviceBuffer fc_cell, fc_delta, fc_T, fc_eta, fc_p, fc_x, ic_N, ic_E, ic_q, chunks, leaves;
    DeviceBuffer rel_total, rel_kin, patch_out, patch_sum, slice_sum, slice_status, slice_maxT, sigma_em, sc_el, sc_val;
    hall_gpu::HostBuffer stage;

    explicit Impl(const GpuOptions& o) : block(o.block) {
        if (o.session) sp = o.session;
        else {
            hall_gpu::SessionOptions so; so.require_uuid = o.require_uuid; so.nvrtc_library = o.nvrtc_library;
            own = std::make_unique<hall_gpu::Session>(so); sp = own.get();
        }
        hall_gpu::CompileRequest req;
        req.main = {"phase_owner_kernels.cu", kMain};
        req.headers = {{"phase_kernels.h", kKernels}, {"phase_law.h", kLaw}, {"det_math.h", kMath}, {"state_digest.h", kDigest}};
        mod = s().compile_and_load(req);
        auto K = [&](const char* n) { return s().kernel(mod, n); };
        kthermo = K("k_thermo"); kscatter = K("k_scatter"); ksigma = K("k_sigma"); kjoule = K("k_joule");
        kfcp = K("k_fc_prepare"); kfcpr = K("k_fc_pressure"); kic = K("k_ic_flux"); kstep = K("k_step"); kslice = K("k_slice");
        kquote = K("k_quote"); krelease = K("k_release"); kpatch = K("k_patch"); kcopy = K("k_copy_i32"); kchunks = K("k_chunks"); kleaves = K("k_leaves");
    }
    ~Impl() {
        // buffers of a borrowed session are released explicitly; an own session releases everything at close
        if (!own && sp && !sp->failed()) {
            for (DeviceBuffer* b : {&law, &row_dy, &row_em, &em_row_begin, &gx, &gy, &vref, &mref, &E[0], &E[1], &E[2], &N[0], &N[1], &N[2], &flags[0], &flags[1], &flags[2],
                                    &fcm[0], &fcm[1], &fcm[2], &T, &k, &sig, &m, &status, &joule, &em_energy, &em_JB, &patch_pb, &el_of, &el_cell, &el_face, &el_area, &el_patch,
                                    &patch_begin, &in_load, &in_dep_N, &in_dep_E, &in_rel_N, &in_rel_E, &acc_N, &acc_Eref, &acc_Ekin, &acc_Eint, &acc_bnd, &acc_rad, &quote,
                                    &fc_cell, &fc_delta, &fc_T, &fc_eta, &fc_p, &fc_x, &ic_N, &ic_E, &ic_q, &chunks, &leaves, &rel_total, &rel_kin, &patch_out, &patch_sum,
                                    &slice_sum, &slice_status, &slice_maxT, &sigma_em, &sc_el, &sc_val})
                if (*b) sp->release(*b);
            if (stage) sp->release(stage);
        }
    }
    void ensure(DeviceBuffer& b, std::size_t bytes, const char* tag) {
        if (bytes == 0) bytes = 16;
        if (b && b.bytes >= bytes) return;
        if (b) s().release(b);
        b = s().alloc(bytes, tag);
    }
    void ensure_stage(std::size_t bytes) {
        if (stage && stage.bytes >= bytes) return;
        s().synchronize("stage");
        if (stage) s().release(stage);
        stage = s().alloc_pinned(std::max<std::size_t>(bytes, 1 << 20), "stage");
    }
    void up(DeviceBuffer& d, const void* src, std::size_t bytes) {
        if (!bytes) return;
        ensure_stage(bytes);
        s().synchronize("upload staging");
        std::memcpy(stage.pointer, src, bytes);
        s().upload(d, 0, stage, 0, bytes);
        s().synchronize("upload");
    }
    void down(void* dst, const DeviceBuffer& d, std::size_t bytes) {
        if (!bytes) return;
        ensure_stage(bytes);
        s().download(stage, 0, d, 0, bytes);
        s().synchronize("download");
        std::memcpy(dst, stage.pointer, bytes);
    }
    template <class T> void upv(DeviceBuffer& d, const std::vector<T>& v, const char* tag) { ensure(d, v.size() * sizeof(T), tag); up(d, v.data(), v.size() * sizeof(T)); }
    template <class T> std::vector<T> downv(const DeviceBuffer& d, std::size_t n) { std::vector<T> v(n); down(v.data(), d, n * sizeof(T)); return v; }
    template <class T> static T* ptr(const DeviceBuffer& b) { return reinterpret_cast<T*>(b.address); }

    pq_buffers buffers(int cur, int next) const {
        pq_buffers b{};
        b.law = ptr<const pq_law>(law);
        b.row_dy = ptr<const double>(row_dy); b.row_em = ptr<const int>(row_em); b.em_row_begin = ptr<const int>(em_row_begin);
        b.gx = ptr<const double>(gx); b.gy = ptr<const double>(gy); b.vref = ptr<const double>(vref); b.mref = ptr<const double>(mref);
        b.E = ptr<const pq_i128>(E[cur]); b.N = ptr<const pq_i128>(N[cur]); b.flags = ptr<const int>(flags[cur]); b.fcm = ptr<const pq_i128>(fcm[cur]);
        if (next >= 0) { b.E_out = ptr<pq_i128>(E[next]); b.N_out = ptr<pq_i128>(N[next]); b.fcm_out = ptr<pq_i128>(fcm[next]); }
        b.T = ptr<double>(T); b.k = ptr<double>(k); b.sig = ptr<double>(sig); b.m = ptr<double>(m); b.status = ptr<int>(status);
        b.joule_total = ptr<pq_i128>(joule); b.em_energy = ptr<const double>(em_energy); b.em_JB = ptr<const double>(em_JB);
        b.el_of = ptr<const int>(el_of); b.el_cell = ptr<const int>(el_cell); b.el_face = ptr<const int>(el_face);
        b.el_area = ptr<const double>(el_area); b.el_patch = ptr<const int>(el_patch); b.patch_begin = ptr<const int>(patch_begin);
        b.in_load_E = ptr<const pq_i128>(in_load); b.in_dep_N = ptr<const pq_i128>(in_dep_N); b.in_dep_E = ptr<const pq_i128>(in_dep_E);
        b.in_rel_N = ptr<pq_i128>(in_rel_N); b.in_rel_E = ptr<pq_i128>(in_rel_E);
        b.acc_ev_N = ptr<pq_i128>(acc_N); b.acc_ev_Eref = ptr<pq_i128>(acc_Eref); b.acc_ev_Ekin = ptr<pq_i128>(acc_Ekin); b.acc_ev_Eint = ptr<pq_i128>(acc_Eint);
        b.acc_bnd_E = ptr<pq_i128>(acc_bnd); b.acc_rad = ptr<pq_i128>(acc_rad);
        b.el_quote = ptr<double>(quote); b.patch_pb = ptr<const double>(patch_pb);
        b.fc_cell = ptr<int>(fc_cell); b.fc_delta = ptr<double>(fc_delta); b.fc_T = ptr<double>(fc_T); b.fc_eta = ptr<double>(fc_eta);
        b.fc_p = ptr<double>(fc_p); b.fc_x = ptr<double>(fc_x);
        b.ic_N = ptr<pq_i128>(ic_N); b.ic_E = ptr<pq_i128>(ic_E); b.ic_q = ptr<double>(ic_q);
        b.rel_total_N = ptr<const pq_i128>(rel_total); b.rel_kin_J = ptr<const double>(rel_kin);
        b.patch_out = ptr<double>(patch_out); b.patch_sum = ptr<pq_i128>(patch_sum);
        b.slice_sum = ptr<pq_i128>(slice_sum); b.slice_status = ptr<int>(slice_status); b.slice_maxT = ptr<double>(slice_maxT);
        b.sigma_em = ptr<double>(sigma_em);
        b.sc_el = ptr<const int>(sc_el); b.sc_val = ptr<const pq_i128>(sc_val);
        b.in_load_E_w = ptr<pq_i128>(in_load); b.in_dep_N_w = ptr<pq_i128>(in_dep_N); b.in_dep_E_w = ptr<pq_i128>(in_dep_E);
        return b;
    }
    void launch(const hall_gpu::Kernel& kk, int n, pq_dims d, pq_buffers b) {
        if (n <= 0) return;
        void* args[] = {&d, &b};
        s().launch(kk, {unsigned((n + int(block) - 1) / int(block)), 1, 1}, {block, 1, 1}, 0, args);
    }
    void zero_elements() {
        const std::size_t n = std::size_t(std::max(n_el, 1));
        for (DeviceBuffer* b : {&in_load, &in_dep_E, &in_rel_E, &acc_Eref, &acc_Ekin, &acc_Eint, &acc_bnd, &acc_rad}) s().zero(*b, 0, n * sizeof(pq_i128));
        for (DeviceBuffer* b : {&in_dep_N, &in_rel_N, &acc_N}) s().zero(*b, 0, PQ_NEL * n * sizeof(pq_i128));
        s().zero(quote, 0, PQ_NEL * n * sizeof(double));
        s().zero(patch_out, 0, std::size_t(std::max(P, 1)) * PQ_PATCH_OUT * sizeof(double));
        s().zero(patch_sum, 0, std::size_t(std::max(P, 1)) * PQ_PATCH_SUM * sizeof(pq_i128));
    }
    void flow_pass(const pq_dims& d, const pq_buffers& b) {
        launch(kfcp, d.n_fc, d, b);
        launch(kfcpr, d.n_fc, d, b);
        launch(kic, d.n_ic, d, b);
    }
    void finish(int slot, const pq_dims& d) {
        const pq_buffers b = buffers(slot, -1);
        launch(kthermo, d.n_cells, d, b);
        launch(kfcp, d.n_fc, d, b);
        launch(kslice, d.S, d, b);
        launch(kpatch, d.n_patches, d, b);
    }
    void collect(const pq_dims& d, EngineResult& r) {
        r.slice_sum = downv<pq_i128>(slice_sum, std::size_t(d.S) * PQ_SLICE_SUM);
        r.slice_status = downv<int>(slice_status, std::size_t(d.S) * 4);
        r.slice_maxT = downv<double>(slice_maxT, d.S);
        r.patch_sum = downv<pq_i128>(patch_sum, std::size_t(d.n_patches) * PQ_PATCH_SUM);
        r.patch_out = downv<double>(patch_out, std::size_t(d.n_patches) * PQ_PATCH_OUT);
    }
    std::vector<std::array<pq_u32, 8>> section_leaves(int section, pq_u64 n, int slot) {
        const pq_u64 nc = (n + 31) / 32, nl = (n + 1023) / 1024;
        if (nl == 0) return {};
        ensure(chunks, std::size_t(nc) * 32, "chunks");
        ensure(leaves, std::size_t(nl) * 32, "leaves");
        {
            const pq_i128* e = ptr<const pq_i128>(E[slot]); const pq_i128* nn = ptr<const pq_i128>(N[slot]); const int* fl = ptr<const int>(flags[slot]);
            const pq_i128* fm = ptr<const pq_i128>(fcm[slot]); pq_u32* out = ptr<pq_u32>(chunks); int sec = section; pq_u64 cnt = n;
            void* args[] = {&sec, &cnt, &e, &nn, &fl, &fm, &out};
            s().launch(kchunks, {unsigned((nc + 63) / 64), 1, 1}, {64, 1, 1}, 0, args);
        }
        {
            const pq_u32* ch = ptr<const pq_u32>(chunks); pq_u32* out = ptr<pq_u32>(leaves); int sec = section; pq_u64 cnt = n;
            void* args[] = {&sec, &cnt, &ch, &out};
            s().launch(kleaves, {unsigned((nl + 63) / 64), 1, 1}, {64, 1, 1}, 0, args);
        }
        const auto words = downv<pq_u32>(leaves, std::size_t(nl) * 8);
        std::vector<std::array<pq_u32, 8>> lv(nl);
        for (pq_u64 q = 0; q < nl; ++q) for (int w = 0; w < 8; ++w) lv[q][w] = words[8 * q + w];
        return lv;
    }
};

GpuEngine::GpuEngine(const GpuOptions& o) : p_(new Impl(o)) {}
GpuEngine::~GpuEngine() { delete p_; }
bool GpuEngine::borrowed_session() const { return !p_->own; }

void GpuEngine::configure(const pq_law& law, const pq_dims& dims, const HostGeometry& g) {
    auto& I = *p_;
    I.base = dims; I.P = dims.n_patches;
    I.ensure(I.law, sizeof(pq_law), "law"); I.up(I.law, &law, sizeof(pq_law));
    I.upv(I.row_dy, g.row_dy, "row_dy"); I.upv(I.row_em, g.row_em, "row_em"); I.upv(I.em_row_begin, g.em_row_begin, "em_row_begin");
    I.upv(I.gx, g.gx, "gx"); I.upv(I.gy, g.gy, "gy"); I.upv(I.vref, g.vref, "vref"); I.upv(I.mref, g.mref, "mref");
    const std::size_t n = dims.n_cells, nfc = std::max(dims.n_fc, 1), nic = std::max(dims.n_ic, 1), P = std::max(I.P, 1);
    for (int s = 0; s < 3; ++s) {
        I.ensure(I.E[s], n * sizeof(pq_i128), "E"); I.ensure(I.N[s], PQ_NEL * n * sizeof(pq_i128), "N"); I.ensure(I.flags[s], n * sizeof(int), "flags");
        I.ensure(I.fcm[s], PQ_NEL * nfc * sizeof(pq_i128), "fcm");
        I.s().zero(I.E[s], 0, n * sizeof(pq_i128)); I.s().zero(I.N[s], 0, PQ_NEL * n * sizeof(pq_i128)); I.s().zero(I.flags[s], 0, n * sizeof(int));
        I.s().zero(I.fcm[s], 0, PQ_NEL * nfc * sizeof(pq_i128));
    }
    for (DeviceBuffer* b : {&I.T, &I.k, &I.sig, &I.m}) { I.ensure(*b, n * sizeof(double), "scratch"); I.s().zero(*b, 0, n * sizeof(double)); }
    I.ensure(I.status, n * sizeof(int), "status"); I.s().zero(I.status, 0, n * sizeof(int));
    I.ensure(I.joule, n * sizeof(pq_i128), "joule"); I.s().zero(I.joule, 0, n * sizeof(pq_i128));
    I.ensure(I.em_energy, std::size_t(dims.n_em) * sizeof(double), "em_energy"); I.s().zero(I.em_energy, 0, std::size_t(dims.n_em) * sizeof(double));
    I.ensure(I.em_JB, 6 * std::size_t(dims.n_em) * sizeof(double), "em_JB"); I.s().zero(I.em_JB, 0, 6 * std::size_t(dims.n_em) * sizeof(double));
    I.ensure(I.patch_pb, P * sizeof(double), "patch_pb"); I.s().zero(I.patch_pb, 0, P * sizeof(double));
    I.ensure(I.sigma_em, std::size_t(dims.n_em) * sizeof(double), "sigma_em");
    I.ensure(I.slice_sum, std::size_t(dims.S) * PQ_SLICE_SUM * sizeof(pq_i128), "slice_sum");
    I.ensure(I.slice_status, std::size_t(dims.S) * 4 * sizeof(int), "slice_status");
    I.ensure(I.slice_maxT, std::size_t(dims.S) * sizeof(double), "slice_maxT");
    I.ensure(I.rel_total, PQ_NEL * P * sizeof(pq_i128), "rel_total"); I.ensure(I.rel_kin, PQ_NEL * P * sizeof(double), "rel_kin");
    I.ensure(I.patch_out, P * PQ_PATCH_OUT * sizeof(double), "patch_out"); I.ensure(I.patch_sum, P * PQ_PATCH_SUM * sizeof(pq_i128), "patch_sum");
    I.ensure(I.fc_cell, nfc * sizeof(int), "fc_cell");
    for (DeviceBuffer* b : {&I.fc_delta, &I.fc_T, &I.fc_eta, &I.fc_p}) { I.ensure(*b, nfc * sizeof(double), "fc"); I.s().zero(*b, 0, nfc * sizeof(double)); }
    I.ensure(I.fc_x, PQ_NEL * nfc * sizeof(double), "fc_x"); I.s().zero(I.fc_x, 0, PQ_NEL * nfc * sizeof(double));
    I.ensure(I.ic_N, PQ_NEL * nic * sizeof(pq_i128), "ic_N"); I.s().zero(I.ic_N, 0, PQ_NEL * nic * sizeof(pq_i128));
    I.ensure(I.ic_E, nic * sizeof(pq_i128), "ic_E"); I.s().zero(I.ic_E, 0, nic * sizeof(pq_i128));
    I.ensure(I.ic_q, nic * sizeof(double), "ic_q"); I.s().zero(I.ic_q, 0, nic * sizeof(double));
    I.ensure(I.sc_el, 64 * sizeof(int), "sc_el"); I.ensure(I.sc_val, 64 * PQ_SC * sizeof(pq_i128), "sc_val");
    I.s().synchronize("configure");
}

void GpuEngine::set_elements(const HostElements& el) {
    auto& I = *p_;
    I.n_el = int(el.el_cell.size());
    I.upv(I.el_of, el.el_of, "el_of"); I.upv(I.el_cell, el.el_cell, "el_cell"); I.upv(I.el_face, el.el_face, "el_face");
    I.upv(I.el_area, el.el_area, "el_area"); I.upv(I.el_patch, el.el_patch, "el_patch"); I.upv(I.patch_begin, el.patch_begin, "patch_begin");
    const std::size_t n = std::size_t(std::max(I.n_el, 1)), P = std::max(I.P, 1);
    for (DeviceBuffer* b : {&I.in_load, &I.in_dep_E, &I.in_rel_E, &I.acc_Eref, &I.acc_Ekin, &I.acc_Eint, &I.acc_bnd, &I.acc_rad}) I.ensure(*b, n * sizeof(pq_i128), "element");
    for (DeviceBuffer* b : {&I.in_dep_N, &I.in_rel_N, &I.acc_N}) I.ensure(*b, PQ_NEL * n * sizeof(pq_i128), "element6");
    I.ensure(I.quote, PQ_NEL * n * sizeof(double), "quote");
    I.zero_elements();
    I.s().zero(I.rel_total, 0, PQ_NEL * P * sizeof(pq_i128)); I.s().zero(I.rel_kin, 0, PQ_NEL * P * sizeof(double));
    I.s().synchronize("set_elements");
}

void GpuEngine::write_state(int slot, const HostState& st) {
    auto& I = *p_;
    I.up(I.E[slot], st.E.data(), st.E.size() * sizeof(pq_i128));
    I.up(I.N[slot], st.N.data(), st.N.size() * sizeof(pq_i128));
    I.up(I.flags[slot], st.flags.data(), st.flags.size() * sizeof(int));
    I.up(I.fcm[slot], st.fcm.data(), st.fcm.size() * sizeof(pq_i128));
}
HostState GpuEngine::read_state(int slot) {
    auto& I = *p_;
    const std::size_t n = I.base.n_cells, nfc = std::max(I.base.n_fc, 1);
    HostState st;
    st.E = I.downv<pq_i128>(I.E[slot], n); st.N = I.downv<pq_i128>(I.N[slot], PQ_NEL * n); st.flags = I.downv<int>(I.flags[slot], n);
    st.fcm = I.downv<pq_i128>(I.fcm[slot], PQ_NEL * nfc);
    return st;
}
void GpuEngine::evaluate(int slot, const pq_dims& d, const std::vector<double>& pb, EngineResult& r, std::vector<double>* sigma_em) {
    auto& I = *p_;
    I.s().zero(I.status, 0, std::size_t(d.n_cells) * sizeof(int));
    I.up(I.patch_pb, pb.data(), pb.size() * sizeof(double));
    I.finish(slot, d);
    if (sigma_em) { I.launch(I.ksigma, d.n_em, d, I.buffers(slot, -1)); *sigma_em = I.downv<double>(I.sigma_em, d.n_em); }
    I.collect(d, r);
}
void GpuEngine::step(int from, const pq_dims& d0, const StepData& in, EngineResult& r) {
    auto& I = *p_;
    const auto t0 = std::chrono::steady_clock::now();
    pq_dims d = d0;
    I.s().zero(I.status, 0, std::size_t(d.n_cells) * sizeof(int));
    I.zero_elements();
    I.up(I.em_energy, in.em_energy.data(), in.em_energy.size() * sizeof(double));
    I.up(I.em_JB, in.em_JB.data(), in.em_JB.size() * sizeof(double));
    I.up(I.patch_pb, in.patch_pb.data(), in.patch_pb.size() * sizeof(double));
    I.up(I.rel_total, in.rel_total_N.data(), in.rel_total_N.size() * sizeof(pq_i128));
    I.up(I.rel_kin, in.rel_kin_J.data(), in.rel_kin_J.size() * sizeof(double));
    if (!in.sc_el.empty()) { I.upv(I.sc_el, in.sc_el, "sc_el"); I.upv(I.sc_val, in.sc_val, "sc_val"); }
    {
        const pq_buffers b = I.buffers(from, -1);
        I.launch(I.kthermo, d.n_cells, d, b);
        pq_dims ds = d; ds.n_scatter = int(in.sc_el.size());
        I.launch(I.kscatter, ds.n_scatter, ds, b);
        if (d.evap_mode == 2) { I.launch(I.kquote, d.n_el, d, b); I.launch(I.krelease, d.n_patches, d, b); }
        I.launch(I.kjoule, d.n_em, d, b);
    }
    const auto [w0, w1] = work_slots(from);
    int cur = from;
    for (int s = 0; s < d.n_sub; ++s) {
        const int next = (s % 2 == 0) ? w0 : w1;
        const pq_buffers b = I.buffers(cur, next);
        if (s > 0) I.launch(I.kthermo, d.n_cells, d, b);
        d.sub = s;
        if (d.flow) I.flow_pass(d, b);
        I.launch(I.kstep, d.n_cells, d, b);
        {
            const int* a = Impl::ptr<const int>(I.flags[cur]); int* o = Impl::ptr<int>(I.flags[next]); int n = d.n_cells;
            void* args[] = {&a, &o, &n};
            I.s().launch(I.kcopy, {unsigned((n + int(I.block) - 1) / int(I.block)), 1, 1}, {I.block, 1, 1}, 0, args);
        }
        cur = next;
    }
    I.finish(cur, d);
    I.collect(d, r);
    r.trial_slot = cur;
    last_step_ms_ = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
std::vector<double> GpuEngine::temperatures(int slot, const pq_dims& d) {
    auto& I = *p_;
    I.launch(I.kthermo, d.n_cells, d, I.buffers(slot, -1));
    return I.downv<double>(I.T, d.n_cells);
}
std::vector<double> GpuEngine::face_displacements(int slot, const pq_dims& d) {
    auto& I = *p_;
    const pq_buffers b = I.buffers(slot, -1);
    I.launch(I.kthermo, d.n_cells, d, b);
    I.launch(I.kfcp, d.n_fc, d, b);
    return I.downv<double>(I.fc_eta, d.n_fc);
}
StateRoots GpuEngine::state_roots(int slot, const pq_dims& d, std::vector<std::array<pq_u32, 8>>* cell_leaves) {
    auto& I = *p_;
    StateRoots r;
    auto lc = I.section_leaves(PQ_SEC_CELLS, pq_u64(d.n_cells), slot);
    auto lf = I.section_leaves(PQ_SEC_FC, pq_u64(d.n_fc), slot);
    r.cell_leaves = lc.size(); r.fc_leaves = lf.size();
    if (cell_leaves) *cell_leaves = lc;
    r.cells = pq_section_root(PQ_SEC_CELLS, pq_u64(d.n_cells), lc);
    r.fcm = pq_section_root(PQ_SEC_FC, pq_u64(d.n_fc), lf);
    return r;
}
std::string GpuEngine::record_json() const {
    const auto& I = *p_;
    std::ostringstream o;
    o << "{\"device\":" << hall_gpu::device_json(I.sp->device()) << ",\"compile\":" << hall_gpu::compile_json(I.mod.record)
      << ",\"borrowed_session\":" << (I.own ? "false" : "true") << "}";
    return o.str();
}

}  // namespace ss304_phase_v1_1
