// gpu_engine.cpp -- see gpu_engine.hpp.
#include "gpu_engine.hpp"
#include "common/hall_gpu_common.hpp"
#include <chrono>
#include <cstring>
#include <map>
#include <sstream>

namespace ss304_phase_v1_2 {
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
}  // namespace
#include "gpu_program.h"
namespace {
const auto& kMain = kPxMain;

using hall_gpu::DeviceBuffer;
}  // namespace

struct GpuEngine::Impl {
    std::unique_ptr<hall_gpu::Session> own;
    hall_gpu::Session* sp = nullptr;
    hall_gpu::Session& s() { return *sp; }
    hall_gpu::Module mod;
    hall_gpu::Kernel kthermo, kscatter, kfcp, kfcpr, kic, kstep, kslice, kquote, krelease, kpatch, kchunks, kleaves;
    hall_gpu::Kernel kpsums, kemstart, kzero, kfcp2, kfcpr2, kif2, kfcpass, kaccum, ksurf, krows, ksigma;
    unsigned block;
    px_dims base{};
    int n_el = 0, P = 0;
    DeviceBuffer law, row_dy, row_em, em_row_begin, gx, gy, vref, mref;
    DeviceBuffer E[3], N[3], flags[3], fcm[3], fcV[3], fcK[3];
    DeviceBuffer T, k, sig, m, status, joule, em_energy, em_JB, patch_pb;
    DeviceBuffer el_of, el_cell, el_face, el_area, el_patch, patch_begin;
    DeviceBuffer in_load, in_dep_N, in_dep_E, in_rel_N, in_rel_E, acc_N, acc_Eref, acc_Ekin, acc_Eint, acc_bnd, acc_rad, quote;
    DeviceBuffer fc_cell, fc_delta, fc_T, fc_eta, fc_p, fc_x, ic_N, ic_E, ic_q, chunks, leaves;
    DeviceBuffer rel_total, rel_kin, patch_out, patch_sum, slice_sum, slice_status, slice_maxT, sigma_em, sc_el, sc_val;
    // v1.2
    DeviceBuffer acc_E, acc_F, acc_bad, acc_P, acc_JB, em_F, em_diag, slice_diag, if_a, if_b, if_dir, fc_if, slice_zn;
    DeviceBuffer fc_lay, fc_nw, fc_m, fc_h, fc_rho, fc_mu, fc_dsdT, fc_sig, fc_f, fc_u, fc_un, fc_a, fc_ps, fc_g, fc_re, fc_s;
    DeviceBuffer if_dm, if_F, if_st, if_N, if_E, if_K, cflow_N, cflow_E, fc_netN, acc_fc, ic_W, acc_ic;
    DeviceBuffer el_dE, el_dN, el_st, el_bad, row_sum, row_info, row_maxT;
    hall_gpu::HostBuffer stage;
    std::size_t up_end = 0;
    bool stage_busy = false;
    int thermo_slot = -1;   // slot whose T, k, sigma, m are in the scratch arrays (a deterministic function of its state)
    // profiling
    bool profile = false;
    std::vector<std::pair<std::string, int>> marks;
    std::map<std::string, double> prof_ms;
    long prof_steps = 0;
    void mark(const char* label) { if (profile) marks.push_back({label, s().mark(label)}); }
    void collect_profile() {
        if (!profile) return;
        for (std::size_t i = 1; i < marks.size(); ++i) prof_ms[marks[i].first] += s().elapsed_ms(marks[i - 1].second, marks[i].second);
        marks.clear(); s().clear_markers(); ++prof_steps;
    }

    explicit Impl(const GpuOptions& o) : block(o.block), profile(o.profile) {
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
        kthermo = K("k_thermo"); kscatter = K("k_scatter");
        kfcp = K("k_fc_prepare"); kfcpr = K("k_fc_pressure"); kic = K("k_ic_flux"); kstep = K("k_step"); kslice = K("k_slice");
        kquote = K("k_quote"); krelease = K("k_release"); kpatch = K("k_patch"); kchunks = K("k_chunks"); kleaves = K("k_leaves");
        kpsums = K("k_patch_sums"); kemstart = K("k_em_start"); kzero = K("k_zero_step"); kfcp2 = K("k_fc_prepare2"); kfcpr2 = K("k_fc_pressure2");
        kif2 = K("k_if_flux2"); kfcpass = K("k_fc_pass"); kaccum = K("k_accumulate");
        ksurf = K("k_surface"); krows = K("k_rows"); ksigma = K("k_sigma");
    }
    std::vector<DeviceBuffer*> all_buffers() {
        std::vector<DeviceBuffer*> v{&law, &row_dy, &row_em, &em_row_begin, &gx, &gy, &vref, &mref, &T, &k, &sig, &m, &status, &joule, &em_energy, &em_JB, &patch_pb,
            &el_of, &el_cell, &el_face, &el_area, &el_patch, &patch_begin, &in_load, &in_dep_N, &in_dep_E, &in_rel_N, &in_rel_E, &acc_N, &acc_Eref, &acc_Ekin,
            &acc_Eint, &acc_bnd, &acc_rad, &quote, &fc_cell, &fc_delta, &fc_T, &fc_eta, &fc_p, &fc_x, &ic_N, &ic_E, &ic_q, &chunks, &leaves, &rel_total, &rel_kin,
            &patch_out, &patch_sum, &slice_sum, &slice_status, &slice_maxT, &sigma_em, &sc_el, &sc_val,
            &acc_E, &acc_F, &acc_bad, &acc_P, &acc_JB, &em_F, &em_diag, &slice_diag, &if_a, &if_b, &if_dir, &fc_if, &slice_zn,
            &fc_lay, &fc_nw, &fc_m, &fc_h, &fc_rho, &fc_mu, &fc_dsdT, &fc_sig, &fc_f, &fc_u, &fc_un, &fc_a, &fc_ps, &fc_g, &fc_re, &fc_s,
            &if_dm, &if_F, &if_st, &if_N, &if_E, &if_K, &cflow_N, &cflow_E, &fc_netN, &acc_fc, &ic_W, &acc_ic,
            &el_dE, &el_dN, &el_st, &el_bad, &row_sum, &row_info, &row_maxT};
        for (int q = 0; q < 3; ++q) for (DeviceBuffer* b : {&E[q], &N[q], &flags[q], &fcm[q], &fcV[q], &fcK[q]}) v.push_back(b);
        return v;
    }
    ~Impl() {
        // buffers of a borrowed session are released explicitly; an own session releases everything at close
        if (!own && sp && !sp->failed()) {
            for (DeviceBuffer* b : all_buffers()) if (*b) sp->release(*b);
            if (stage) sp->release(stage);
        }
    }
    void ensure(DeviceBuffer& b, std::size_t bytes, const char* tag) {
        if (bytes == 0) bytes = 16;
        if (b && b.bytes >= bytes) return;
        if (b) { s().synchronize("reallocate"); s().release(b); }
        b = s().alloc(bytes, tag);
    }
    void ensure_stage(std::size_t bytes) {
        if (stage && stage.bytes >= bytes) return;
        s().synchronize("stage");
        if (stage) s().release(stage);
        stage = s().alloc_pinned(std::max<std::size_t>(bytes, 1 << 20), "stage");
    }
    void quiet() { if (stage_busy) { s().synchronize("stage reuse"); stage_busy = false; } }
    void up(DeviceBuffer& d, const void* src, std::size_t bytes) {
        if (!bytes) return;
        quiet();
        ensure_stage(bytes);
        std::memcpy(stage.pointer, src, bytes);
        s().upload(d, 0, stage, 0, bytes);
        s().synchronize("upload");
    }
    void down(void* dst, const DeviceBuffer& d, std::size_t bytes) {
        if (!bytes) return;
        quiet();
        ensure_stage(bytes);
        s().download(stage, 0, d, 0, bytes);
        s().synchronize("download");
        std::memcpy(dst, stage.pointer, bytes);
    }
    // Staged transfers of one step: every upload through one staging area without waiting, the downloads after the
    // kernels, one synchronization.
    static std::size_t al(std::size_t x) { return (x + 255) & ~std::size_t(255); }
    struct Up { DeviceBuffer* d; const void* src; std::size_t bytes; };
    struct Down { void* dst; const DeviceBuffer* d; std::size_t bytes; };
    void put(const std::vector<Up>& ups, std::size_t down_bytes) {
        quiet();
        std::size_t tot = 0;
        for (const auto& u : ups) tot += al(u.bytes);
        ensure_stage(tot + down_bytes);
        std::size_t off = 0;
        for (const auto& u : ups) {
            if (!u.bytes) continue;
            std::memcpy(static_cast<char*>(stage.pointer) + off, u.src, u.bytes);
            s().upload(*u.d, 0, stage, off, u.bytes);
            off += al(u.bytes);
        }
        up_end = off; stage_busy = true;
    }
    void get(const std::vector<Down>& downs) {
        std::size_t off = up_end;
        for (const auto& d : downs) { if (d.bytes) s().download(stage, off, *d.d, 0, d.bytes); off += al(d.bytes); }
        s().synchronize("step results");
        stage_busy = false;
        off = up_end;
        for (const auto& d : downs) { if (d.bytes) std::memcpy(d.dst, static_cast<char*>(stage.pointer) + off, d.bytes); off += al(d.bytes); }
        up_end = 0;
    }
    template <class T> void upv(DeviceBuffer& d, const std::vector<T>& v, const char* tag) { ensure(d, v.size() * sizeof(T), tag); up(d, v.data(), v.size() * sizeof(T)); }
    template <class T> std::vector<T> downv(const DeviceBuffer& d, std::size_t n) { std::vector<T> v(n); down(v.data(), d, n * sizeof(T)); return v; }
    template <class T> static T* ptr(const DeviceBuffer& b) { return reinterpret_cast<T*>(b.address); }

    px_buffers buffers(int cur, int next) const {
        px_buffers b{};
        b.law = ptr<const px_law>(law);
        b.row_dy = ptr<const double>(row_dy); b.row_em = ptr<const int>(row_em); b.em_row_begin = ptr<const int>(em_row_begin);
        b.gx = ptr<const double>(gx); b.gy = ptr<const double>(gy); b.vref = ptr<const double>(vref); b.mref = ptr<const double>(mref);
        b.E = ptr<const px_i128>(E[cur]); b.N = ptr<const px_i128>(N[cur]); b.flags = ptr<const int>(flags[cur]); b.fcm = ptr<const px_i128>(fcm[cur]);
        b.fcV = ptr<const double>(fcV[cur]); b.fcK = ptr<const px_i128>(fcK[cur]);
        if (next >= 0) {
            b.E_out = ptr<px_i128>(E[next]); b.N_out = ptr<px_i128>(N[next]); b.fcm_out = ptr<px_i128>(fcm[next]);
            b.fcV_out = ptr<double>(fcV[next]); b.fcK_out = ptr<px_i128>(fcK[next]);
            b.flags_out = ptr<int>(flags[next]);
        }
        b.T = ptr<double>(T); b.k = ptr<double>(k); b.sig = ptr<double>(sig); b.m = ptr<double>(m); b.status = ptr<int>(status);
        b.joule_total = ptr<px_i128>(joule); b.em_energy = ptr<const double>(em_energy); b.em_JB = ptr<const double>(em_JB);
        b.el_of = ptr<const int>(el_of); b.el_cell = ptr<const int>(el_cell); b.el_face = ptr<const int>(el_face);
        b.el_area = ptr<const double>(el_area); b.el_patch = ptr<const int>(el_patch); b.patch_begin = ptr<const int>(patch_begin);
        b.in_load_E = ptr<const px_i128>(in_load); b.in_dep_N = ptr<const px_i128>(in_dep_N); b.in_dep_E = ptr<const px_i128>(in_dep_E);
        b.in_rel_N = ptr<px_i128>(in_rel_N); b.in_rel_E = ptr<px_i128>(in_rel_E);
        b.acc_ev_N = ptr<px_i128>(acc_N); b.acc_ev_Eref = ptr<px_i128>(acc_Eref); b.acc_ev_Ekin = ptr<px_i128>(acc_Ekin); b.acc_ev_Eint = ptr<px_i128>(acc_Eint);
        b.acc_bnd_E = ptr<px_i128>(acc_bnd); b.acc_rad = ptr<px_i128>(acc_rad);
        b.el_quote = ptr<double>(quote); b.patch_pb = ptr<const double>(patch_pb);
        b.fc_cell = ptr<int>(fc_cell); b.fc_delta = ptr<double>(fc_delta); b.fc_T = ptr<double>(fc_T); b.fc_eta = ptr<double>(fc_eta);
        b.fc_p = ptr<double>(fc_p); b.fc_x = ptr<double>(fc_x);
        b.ic_N = ptr<px_i128>(ic_N); b.ic_E = ptr<px_i128>(ic_E); b.ic_q = ptr<double>(ic_q);
        b.rel_total_N = ptr<const px_i128>(rel_total); b.rel_kin_J = ptr<const double>(rel_kin);
        b.patch_out = ptr<double>(patch_out); b.patch_sum = ptr<px_i128>(patch_sum);
        b.slice_sum = ptr<px_i128>(slice_sum); b.slice_status = ptr<int>(slice_status); b.slice_maxT = ptr<double>(slice_maxT);
        b.sigma_em = ptr<double>(sigma_em);
        b.sc_el = ptr<const int>(sc_el); b.sc_val = ptr<const px_i128>(sc_val);
        b.in_load_E_w = ptr<px_i128>(in_load); b.in_dep_N_w = ptr<px_i128>(in_dep_N); b.in_dep_E_w = ptr<px_i128>(in_dep_E);
        b.acc_E = ptr<const px_i128>(acc_E); b.acc_F = ptr<const double>(acc_F); b.acc_bad = ptr<const int>(acc_bad);
        b.em_F = ptr<double>(em_F); b.em_diag = ptr<double>(em_diag); b.slice_diag = ptr<double>(slice_diag);
        b.if_a = ptr<const int>(if_a); b.if_b = ptr<const int>(if_b); b.if_dir = ptr<const int>(if_dir); b.fc_if = ptr<const int>(fc_if);
        b.slice_zn = ptr<const int>(slice_zn);
        b.fc_lay = ptr<int>(fc_lay); b.fc_nw = ptr<int>(fc_nw); b.fc_m = ptr<double>(fc_m); b.fc_h = ptr<double>(fc_h); b.fc_rho = ptr<double>(fc_rho);
        b.fc_mu = ptr<double>(fc_mu); b.fc_dsdT = ptr<double>(fc_dsdT); b.fc_sig = ptr<double>(fc_sig); b.fc_f = ptr<double>(fc_f);
        b.fc_u = ptr<double>(fc_u); b.fc_un = ptr<double>(fc_un); b.fc_a = ptr<double>(fc_a); b.fc_ps = ptr<double>(fc_ps); b.fc_g = ptr<double>(fc_g);
        b.fc_re = ptr<double>(fc_re); b.fc_s = ptr<double>(fc_s);
        b.if_dm = ptr<double>(if_dm); b.if_F = ptr<double>(if_F); b.if_st = ptr<int>(if_st); b.if_N = ptr<px_i128>(if_N); b.if_E = ptr<px_i128>(if_E);
        b.if_K = ptr<px_i128>(if_K); b.cflow_N = ptr<px_i128>(cflow_N); b.cflow_E = ptr<px_i128>(cflow_E); b.fc_netN = ptr<px_i128>(fc_netN);
        b.acc_fc = ptr<px_i128>(acc_fc); b.ic_W = ptr<px_i128>(ic_W); b.acc_ic = ptr<px_i128>(acc_ic);
        b.el_dE = ptr<px_i128>(el_dE); b.el_dN = ptr<px_i128>(el_dN); b.el_st = ptr<int>(el_st); b.el_bad = ptr<int>(el_bad);
        b.row_sum = ptr<px_i128>(row_sum); b.row_info = ptr<int>(row_info); b.row_maxT = ptr<double>(row_maxT);
        return b;
    }
    void launch(const hall_gpu::Kernel& kk, int n, px_dims d, px_buffers b) {
        if (n <= 0) return;
        void* args[] = {&d, &b};
        s().launch(kk, {unsigned((n + int(block) - 1) / int(block)), 1, 1}, {block, 1, 1}, 0, args);
    }
    void zero_elements(const px_dims& d) {
        const std::size_t n = std::size_t(std::max(n_el, 1));
        for (DeviceBuffer* b : {&in_load, &in_dep_E, &in_rel_E, &acc_Eref, &acc_Ekin, &acc_Eint, &acc_bnd, &acc_rad}) s().zero(*b, 0, n * sizeof(px_i128));
        for (DeviceBuffer* b : {&in_dep_N, &in_rel_N, &acc_N}) s().zero(*b, 0, PX_NEL * n * sizeof(px_i128));
        s().zero(quote, 0, PX_NEL * n * sizeof(double));
        s().zero(patch_out, 0, std::size_t(std::max(P, 1)) * PX_PATCH_OUT * sizeof(double));
        s().zero(patch_sum, 0, std::size_t(std::max(P, 1)) * PX_PATCH_SUM * sizeof(px_i128));
        const std::size_t nfc = std::max(d.n_fc, 1), nic = std::max(d.n_ic, 1);
        if (d.flow && d.flow_model == 2) { s().zero(acc_fc, 0, PX_FC_ACC * nfc * sizeof(px_i128)); s().zero(fc_s, 0, nfc * sizeof(double)); }
        if (d.flow && d.book_work) s().zero(acc_ic, 0, 2 * nic * sizeof(px_i128));
    }
    void flow_pass(const px_dims& d, const px_buffers& b) {
        if (d.flow_model == 2) {
            launch(kfcp2, d.n_fc, d, b);
            launch(kfcpr2, d.n_fc, d, b);
            launch(kif2, d.n_if, d, b);
            launch(kfcpass, d.n_fc, d, b);
            return;
        }
        launch(kfcp, d.n_fc, d, b);
        launch(kfcpr, d.n_fc, d, b);
        launch(kic, d.n_ic, d, b);
    }
    void finish(int slot, const px_dims& d, bool full) {
        const px_buffers b = buffers(slot, -1);
        launch(kthermo, d.n_cells, d, b);
        thermo_slot = slot;
        launch(d.flow && d.flow_model == 2 ? kfcp2 : kfcp, d.n_fc, d, b);
        launch(ksigma, d.n_em, d, b);
        launch(krows, d.S * d.ny, d, b);
        launch(kslice, d.S, d, b);
        launch(full ? kpatch : kpsums, d.n_patches, d, b);
    }
    std::vector<Down> results(const px_dims& d, EngineResult& r) {
        r.slice_sum.resize(std::size_t(d.S) * PX_SLICE_SUM); r.slice_status.resize(std::size_t(d.S) * 4); r.slice_maxT.resize(d.S);
        r.slice_diag.resize(std::size_t(d.S) * 4); r.patch_sum.resize(std::size_t(d.n_patches) * PX_PATCH_SUM);
        r.patch_out.resize(std::size_t(d.n_patches) * PX_PATCH_OUT); r.sigma_em.resize(d.n_em);
        return {{r.slice_sum.data(), &slice_sum, r.slice_sum.size() * sizeof(px_i128)}, {r.slice_status.data(), &slice_status, r.slice_status.size() * sizeof(int)},
                {r.slice_maxT.data(), &slice_maxT, r.slice_maxT.size() * sizeof(double)}, {r.slice_diag.data(), &slice_diag, r.slice_diag.size() * sizeof(double)},
                {r.patch_sum.data(), &patch_sum, r.patch_sum.size() * sizeof(px_i128)}, {r.patch_out.data(), &patch_out, r.patch_out.size() * sizeof(double)},
                {r.sigma_em.data(), &sigma_em, r.sigma_em.size() * sizeof(double)}};
    }
    static std::size_t bytes_of(const std::vector<Down>& v) { std::size_t t = 0; for (const auto& d : v) t += al(d.bytes); return t; }
    std::vector<std::array<px_u32, 8>> section_leaves(int section, px_u64 n, int slot) {
        const px_u64 nc = (n + 31) / 32, nl = (n + 1023) / 1024;
        if (nl == 0) return {};
        ensure(chunks, std::size_t(nc) * 32, "chunks");
        ensure(leaves, std::size_t(nl) * 32, "leaves");
        {
            px_digest_src g{ptr<const px_i128>(E[slot]), ptr<const px_i128>(N[slot]), ptr<const int>(flags[slot]), ptr<const px_i128>(fcm[slot]),
                            ptr<const double>(fcV[slot]), ptr<const px_i128>(fcK[slot])};
            px_u32* out = ptr<px_u32>(chunks); int sec = section; px_u64 cnt = n;
            void* args[] = {&sec, &cnt, &g, &out};
            s().launch(kchunks, {unsigned((nc + 63) / 64), 1, 1}, {64, 1, 1}, 0, args);
        }
        {
            const px_u32* ch = ptr<const px_u32>(chunks); px_u32* out = ptr<px_u32>(leaves); int sec = section; px_u64 cnt = n;
            void* args[] = {&sec, &cnt, &ch, &out};
            s().launch(kleaves, {unsigned((nl + 63) / 64), 1, 1}, {64, 1, 1}, 0, args);
        }
        const auto words = downv<px_u32>(leaves, std::size_t(nl) * 8);
        std::vector<std::array<px_u32, 8>> lv(nl);
        for (px_u64 q = 0; q < nl; ++q) for (int w = 0; w < 8; ++w) lv[q][w] = words[8 * q + w];
        return lv;
    }
};

GpuEngine::GpuEngine(const GpuOptions& o) : p_(new Impl(o)) {}
GpuEngine::~GpuEngine() { delete p_; }
bool GpuEngine::borrowed_session() const { return !p_->own; }

void GpuEngine::configure(const px_law& law, const px_dims& dims, const HostGeometry& g) {
    auto& I = *p_;
    I.base = dims; I.P = dims.n_patches;
    I.ensure(I.law, sizeof(px_law), "law"); I.up(I.law, &law, sizeof(px_law));
    I.upv(I.row_dy, g.row_dy, "row_dy"); I.upv(I.row_em, g.row_em, "row_em"); I.upv(I.em_row_begin, g.em_row_begin, "em_row_begin");
    I.upv(I.gx, g.gx, "gx"); I.upv(I.gy, g.gy, "gy"); I.upv(I.vref, g.vref, "vref"); I.upv(I.mref, g.mref, "mref");
    const std::size_t n = dims.n_cells, nfc = std::max(dims.n_fc, 1), nic = std::max(dims.n_ic, 1), P = std::max(I.P, 1);
    const std::size_t nem = std::max(dims.n_em, 1), nif = std::max(dims.n_if, 1), S = dims.S;
    auto Z = [&](DeviceBuffer& b, std::size_t bytes, const char* tag) { I.ensure(b, bytes, tag); I.s().zero(b, 0, bytes); };
    for (int s = 0; s < 3; ++s) {
        Z(I.E[s], n * sizeof(px_i128), "E"); Z(I.N[s], PX_NEL * n * sizeof(px_i128), "N"); Z(I.flags[s], n * sizeof(int), "flags");
        Z(I.fcm[s], PX_NEL * nfc * sizeof(px_i128), "fcm"); Z(I.fcV[s], 3 * nfc * sizeof(double), "fcV"); Z(I.fcK[s], nfc * sizeof(px_i128), "fcK");
    }
    for (DeviceBuffer* b : {&I.T, &I.k, &I.sig, &I.m}) Z(*b, n * sizeof(double), "scratch");
    Z(I.status, n * sizeof(int), "status");
    Z(I.joule, n * sizeof(px_i128), "joule");
    Z(I.em_energy, nem * sizeof(double), "em_energy");
    Z(I.em_JB, 6 * nem * sizeof(double), "em_JB");
    Z(I.patch_pb, P * sizeof(double), "patch_pb");
    Z(I.sigma_em, nem * sizeof(double), "sigma_em");
    Z(I.slice_sum, S * PX_SLICE_SUM * sizeof(px_i128), "slice_sum");
    Z(I.slice_status, S * 4 * sizeof(int), "slice_status");
    Z(I.slice_maxT, S * sizeof(double), "slice_maxT");
    Z(I.slice_diag, 4 * S * sizeof(double), "slice_diag");
    I.ensure(I.rel_total, PX_NEL * P * sizeof(px_i128), "rel_total"); I.ensure(I.rel_kin, PX_NEL * P * sizeof(double), "rel_kin");
    I.ensure(I.patch_out, P * PX_PATCH_OUT * sizeof(double), "patch_out"); I.ensure(I.patch_sum, P * PX_PATCH_SUM * sizeof(px_i128), "patch_sum");
    Z(I.fc_cell, nfc * sizeof(int), "fc_cell");
    for (DeviceBuffer* b : {&I.fc_delta, &I.fc_T, &I.fc_eta, &I.fc_p}) Z(*b, nfc * sizeof(double), "fc");
    Z(I.fc_x, PX_NEL * nfc * sizeof(double), "fc_x");
    Z(I.ic_N, PX_NEL * nic * sizeof(px_i128), "ic_N");
    Z(I.ic_E, nic * sizeof(px_i128), "ic_E");
    Z(I.ic_q, nic * sizeof(double), "ic_q");
    I.ensure(I.sc_el, 64 * sizeof(int), "sc_el"); I.ensure(I.sc_val, 64 * PX_SC * sizeof(px_i128), "sc_val");
    // v1.2
    Z(I.acc_E, nem * sizeof(px_i128), "acc_E"); Z(I.acc_F, 3 * nem * sizeof(double), "acc_F"); Z(I.acc_bad, nem * sizeof(int), "acc_bad");
    Z(I.acc_P, nem * sizeof(double), "acc_P"); Z(I.acc_JB, 6 * nem * sizeof(double), "acc_JB");
    Z(I.em_F, 3 * nem * sizeof(double), "em_F"); Z(I.em_diag, 2 * nem * sizeof(double), "em_diag");
    auto upi = [&](DeviceBuffer& b, const std::vector<int>& v, std::size_t count, const char* tag) {
        std::vector<int> w = v; if (w.size() < count) w.resize(count, -1); I.upv(b, w, tag);
    };
    upi(I.if_a, g.if_a, nif, "if_a"); upi(I.if_b, g.if_b, nif, "if_b"); upi(I.if_dir, g.if_dir, nif, "if_dir");
    upi(I.fc_if, g.fc_if, 4 * nfc, "fc_if"); upi(I.slice_zn, g.slice_zn, 2 * S, "slice_zn");
    Z(I.fc_lay, nfc * sizeof(int), "fc_lay"); Z(I.fc_nw, nfc * sizeof(int), "fc_nw");
    for (DeviceBuffer* b : {&I.fc_m, &I.fc_h, &I.fc_rho, &I.fc_mu, &I.fc_dsdT, &I.fc_sig, &I.fc_a, &I.fc_ps, &I.fc_g, &I.fc_re, &I.fc_s}) Z(*b, nfc * sizeof(double), "fc2");
    Z(I.fc_f, 3 * nfc * sizeof(double), "fc_f"); Z(I.fc_u, 2 * nfc * sizeof(double), "fc_u"); Z(I.fc_un, 2 * nfc * sizeof(double), "fc_un");
    Z(I.if_dm, nif * sizeof(double), "if_dm"); Z(I.if_F, 3 * nif * sizeof(double), "if_F"); Z(I.if_st, nif * sizeof(int), "if_st");
    Z(I.if_N, PX_NEL * nif * sizeof(px_i128), "if_N"); Z(I.if_E, nif * sizeof(px_i128), "if_E"); Z(I.if_K, nif * sizeof(px_i128), "if_K");
    Z(I.cflow_N, PX_NEL * n * sizeof(px_i128), "cflow_N"); Z(I.cflow_E, n * sizeof(px_i128), "cflow_E"); Z(I.fc_netN, PX_NEL * nfc * sizeof(px_i128), "fc_netN");
    Z(I.acc_fc, PX_FC_ACC * nfc * sizeof(px_i128), "acc_fc"); Z(I.ic_W, nic * sizeof(px_i128), "ic_W"); Z(I.acc_ic, 2 * nic * sizeof(px_i128), "acc_ic");
    const std::size_t rows = S * std::size_t(dims.ny);
    Z(I.row_sum, PX_ROW_SUM * rows * sizeof(px_i128), "row_sum"); Z(I.row_info, 5 * rows * sizeof(int), "row_info"); Z(I.row_maxT, rows * sizeof(double), "row_maxT");
    I.s().synchronize("configure");
}

void GpuEngine::set_elements(const HostElements& el) {
    auto& I = *p_;
    I.n_el = int(el.el_cell.size());
    I.upv(I.el_of, el.el_of, "el_of"); I.upv(I.el_cell, el.el_cell, "el_cell"); I.upv(I.el_face, el.el_face, "el_face");
    I.upv(I.el_area, el.el_area, "el_area"); I.upv(I.el_patch, el.el_patch, "el_patch"); I.upv(I.patch_begin, el.patch_begin, "patch_begin");
    const std::size_t n = std::size_t(std::max(I.n_el, 1)), P = std::max(I.P, 1);
    for (DeviceBuffer* b : {&I.in_load, &I.in_dep_E, &I.in_rel_E, &I.acc_Eref, &I.acc_Ekin, &I.acc_Eint, &I.acc_bnd, &I.acc_rad}) I.ensure(*b, n * sizeof(px_i128), "element");
    for (DeviceBuffer* b : {&I.in_dep_N, &I.in_rel_N, &I.acc_N}) I.ensure(*b, PX_NEL * n * sizeof(px_i128), "element6");
    I.ensure(I.quote, PX_NEL * n * sizeof(double), "quote");
    I.ensure(I.el_dE, n * sizeof(px_i128), "el_dE"); I.ensure(I.el_dN, PX_NEL * n * sizeof(px_i128), "el_dN");
    I.ensure(I.el_st, n * sizeof(int), "el_st"); I.ensure(I.el_bad, n * sizeof(int), "el_bad");
    I.zero_elements(I.base);
    I.s().zero(I.rel_total, 0, PX_NEL * P * sizeof(px_i128)); I.s().zero(I.rel_kin, 0, PX_NEL * P * sizeof(double));
    I.s().synchronize("set_elements");
}

void GpuEngine::write_state(int slot, const HostState& st) {
    auto& I = *p_;
    if (slot == I.thermo_slot) I.thermo_slot = -1;
    I.up(I.E[slot], st.E.data(), st.E.size() * sizeof(px_i128));
    I.up(I.N[slot], st.N.data(), st.N.size() * sizeof(px_i128));
    I.up(I.flags[slot], st.flags.data(), st.flags.size() * sizeof(int));
    I.up(I.fcm[slot], st.fcm.data(), st.fcm.size() * sizeof(px_i128));
    I.up(I.fcV[slot], st.fcV.data(), st.fcV.size() * sizeof(double));
    I.up(I.fcK[slot], st.fcK.data(), st.fcK.size() * sizeof(px_i128));
}
HostState GpuEngine::read_state(int slot) {
    auto& I = *p_;
    const std::size_t n = I.base.n_cells, nfc = std::max(I.base.n_fc, 1);
    HostState st;
    st.E = I.downv<px_i128>(I.E[slot], n); st.N = I.downv<px_i128>(I.N[slot], PX_NEL * n); st.flags = I.downv<int>(I.flags[slot], n);
    st.fcm = I.downv<px_i128>(I.fcm[slot], PX_NEL * nfc);
    st.fcV = I.downv<double>(I.fcV[slot], 3 * nfc); st.fcK = I.downv<px_i128>(I.fcK[slot], nfc);
    return st;
}
void GpuEngine::evaluate(int slot, const px_dims& d0, const std::vector<double>& pb, EngineResult& r, std::vector<double>* sigma_em) {
    auto& I = *p_;
    px_dims d = d0; d.use_acc = 0; d.diag = 0;
    I.s().zero(I.status, 0, std::size_t(d.n_cells) * sizeof(int));
    auto downs = I.results(d, r);
    I.put({{&I.patch_pb, pb.data(), pb.size() * sizeof(double)}}, Impl::bytes_of(downs));
    I.finish(slot, d, true);
    I.get(downs);
    if (sigma_em) *sigma_em = r.sigma_em;
}
void GpuEngine::step(int from, const px_dims& d0, const StepData& in, EngineResult& r) {
    auto& I = *p_;
    const auto t0 = std::chrono::steady_clock::now();
    px_dims d = d0;
    d.use_acc = in.use_acc ? 1 : 0;
    I.mark("begin");
    {
        // one pass zeroes the status, element, patch and accumulator arrays of the step
        px_dims dz = d; dz.n_el = I.n_el;
        int n = std::max({d.n_cells, I.n_el, d.n_patches, d.n_fc, d.n_ic});
        px_buffers bz = I.buffers(from, -1);
        void* args[] = {&dz, &bz, &n};
        I.s().launch(I.kzero, {unsigned((n + int(I.block) - 1) / int(I.block)), 1, 1}, {I.block, 1, 1}, 0, args);
    }
    if (!in.sc_el.empty()) {
        I.ensure(I.sc_el, in.sc_el.size() * sizeof(int), "sc_el"); I.ensure(I.sc_val, in.sc_val.size() * sizeof(px_i128), "sc_val");
    }
    auto downs = I.results(d, r);
    std::vector<Impl::Up> ups{{&I.patch_pb, in.patch_pb.data(), in.patch_pb.size() * sizeof(double)},
                              {&I.sc_el, in.sc_el.data(), in.sc_el.size() * sizeof(int)},
                              {&I.sc_val, in.sc_val.data(), in.sc_val.size() * sizeof(px_i128)}};
    if (d.evap_mode == 2) {   // quantized releases (read only by the release split)
        ups.push_back({&I.rel_total, in.rel_total_N.data(), in.rel_total_N.size() * sizeof(px_i128)});
        ups.push_back({&I.rel_kin, in.rel_kin_J.data(), in.rel_kin_J.size() * sizeof(double)});
    }
    if (!in.use_acc) {
        ups.push_back({&I.em_energy, in.em_energy.data(), in.em_energy.size() * sizeof(double)});
        ups.push_back({&I.em_JB, in.em_JB.data(), in.em_JB.size() * sizeof(double)});
    }
    I.put(ups, Impl::bytes_of(downs));
    I.mark("zero+upload");
    {
        const px_buffers b = I.buffers(from, -1);
        // the scratch thermodynamic state of the accepted slot is still valid when the last finish/evaluate computed it
        if (I.thermo_slot != from) I.launch(I.kthermo, d.n_cells, d, b);
        px_dims ds = d; ds.n_scatter = int(in.sc_el.size());
        I.launch(I.kscatter, ds.n_scatter, ds, b);
        if (d.evap_mode == 2) { I.launch(I.kquote, d.n_el, d, b); I.launch(I.krelease, d.n_patches, d, b); }
        I.launch(I.kemstart, d.n_em, d, b);   // force, start conductivity (diagnostics), Joule partition
    }
    I.mark("start: thermo, scatter, forces, joule");
    const auto [w0, w1] = work_slots(from);
    int cur = from;
    for (int s = 0; s < d.n_sub; ++s) {
        const int next = (s % 2 == 0) ? w0 : w1;
        const px_buffers b = I.buffers(cur, next);
        if (s > 0) { I.launch(I.kthermo, d.n_cells, d, b); I.thermo_slot = cur; }
        I.mark("sub-step thermo");
        d.sub = s;
        if (d.flow) I.flow_pass(d, b);
        I.mark("sub-step melt layer");
        I.launch(I.ksurf, d.n_el, d, b);      // surface physics per element
        I.launch(I.kstep, d.n_cells, d, b);   // includes the flags of the next slot
        I.mark("sub-step cells");
        cur = next;
    }
    {
        const px_buffers b = I.buffers(cur, -1);
        I.launch(I.kthermo, d.n_cells, d, b);
        I.thermo_slot = cur;
        I.launch(d.flow && d.flow_model == 2 ? I.kfcp2 : I.kfcp, d.n_fc, d, b);
        I.mark("finish thermo, face columns");
        I.launch(I.ksigma, d.n_em, d, b);
        I.launch(I.krows, d.S * d.ny, d, b);
        I.launch(I.kslice, d.S, d, b);
        I.mark("finish slice sums, conductivity");
        I.launch(I.kpsums, d.n_patches, d, b);
        I.mark("finish patch sums");
    }
    I.get(downs);
    I.mark("download");
    if (I.profile) { I.s().synchronize("profile"); I.collect_profile(); }
    r.trial_slot = cur;
    last_step_ms_ = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
void GpuEngine::accumulate(const std::vector<double>& power, const std::vector<double>& JB, double dt) {
    auto& I = *p_;
    const int n = I.base.n_em;
    std::vector<Impl::Up> ups{{&I.acc_P, power.data(), std::size_t(n) * sizeof(double)}};
    if (!JB.empty()) ups.push_back({&I.acc_JB, JB.data(), 6 * std::size_t(n) * sizeof(double)});
    I.put(ups, 0);
    accumulate_device(Impl::ptr<const double>(I.acc_P), JB.empty() ? nullptr : Impl::ptr<const double>(I.acc_JB), dt);
}
void GpuEngine::accumulate_device(const double* power, const double* JB, double dt) {
    auto& I = *p_;
    int n = I.base.n_em;
    px_i128* aE = Impl::ptr<px_i128>(I.acc_E); double* aF = Impl::ptr<double>(I.acc_F); int* bad = Impl::ptr<int>(I.acc_bad);
    void* args[] = {&n, &power, &JB, &dt, &aE, &aF, &bad};
    I.s().launch(I.kaccum, {unsigned((n + int(I.block) - 1) / int(I.block)), 1, 1}, {I.block, 1, 1}, 0, args);
}
void GpuEngine::reset_accumulator() {
    auto& I = *p_;
    const std::size_t nem = std::max(I.base.n_em, 1);
    I.s().zero(I.acc_E, 0, nem * sizeof(px_i128)); I.s().zero(I.acc_F, 0, 3 * nem * sizeof(double)); I.s().zero(I.acc_bad, 0, nem * sizeof(int));
}
std::vector<double> GpuEngine::temperatures(int slot, const px_dims& d) {
    auto& I = *p_;
    I.launch(I.kthermo, d.n_cells, d, I.buffers(slot, -1));
    I.thermo_slot = slot;
    return I.downv<double>(I.T, d.n_cells);
}
std::vector<double> GpuEngine::face_displacements(int slot, const px_dims& d) {
    auto& I = *p_;
    const px_buffers b = I.buffers(slot, -1);
    I.launch(I.kthermo, d.n_cells, d, b);
    I.thermo_slot = slot;
    I.launch(I.kfcp, d.n_fc, d, b);
    return I.downv<double>(I.fc_eta, d.n_fc);
}
PxStateRoots GpuEngine::state_roots(int slot, const px_dims& d, std::vector<std::array<px_u32, 8>>* cell_leaves) {
    auto& I = *p_;
    PxStateRoots r;
    auto lc = I.section_leaves(PX_SEC_CELLS, px_u64(d.n_cells), slot);
    auto lf = I.section_leaves(PX_SEC_FC, px_u64(d.n_fc), slot);
    r.cell_leaves = lc.size(); r.fc_leaves = lf.size();
    if (cell_leaves) *cell_leaves = lc;
    r.cells = px_section_root(PX_SEC_CELLS, px_u64(d.n_cells), lc);
    r.fcm = px_section_root(PX_SEC_FC, px_u64(d.n_fc), lf);
    if (d.flow && d.flow_model == 2) {
        auto lv = I.section_leaves(PX_SEC_FLOW, px_u64(d.n_fc), slot);
        r.flow_leaves = lv.size(); r.has_flow = true;
        r.flow = px_section_root(PX_SEC_FLOW, px_u64(d.n_fc), lv);
    }
    return r;
}
std::string GpuEngine::profile_json() const {
    const auto& I = *p_;
    std::ostringstream o;
    o << "{\"steps\":" << I.prof_steps << ",\"device_ms\":{";
    bool first = true;
    for (const auto& [k, v] : I.prof_ms) { o << (first ? "" : ",") << hall_gpu::json_string(k) << ":" << hall_gpu::json_number(v); first = false; }
    o << "}}";
    return o.str();
}
std::string GpuEngine::record_json() const {
    const auto& I = *p_;
    std::ostringstream o;
    o << "{\"device\":" << hall_gpu::device_json(I.sp->device()) << ",\"compile\":" << hall_gpu::compile_json(I.mod.record)
      << ",\"borrowed_session\":" << (I.own ? "false" : "true") << "}";
    return o.str();
}

}  // namespace ss304_phase_v1_2
