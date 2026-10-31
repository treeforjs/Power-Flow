// Hall v3.2 device-resident field step: host engine (see resident_field.hpp). v3.2: admitted Gauss
// ledger words in each slot (state digest, upload/download), k_charge before the ledger update, the
// absolute Gauss admission relative to the ledger, reduction rows of the v3.2 admission values.
#include "resident_field.hpp"
#include "common/hall_gpu_common.hpp"
#include "cg_kernels.hpp"
#include "resident_abi.h"
#include <openssl/sha.h>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <sstream>

namespace hall_corrective_v3::gpu {
namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kAbi[] = {
#embed "resident_abi.h"
    , 0};
const char kKernels[] = {
#embed "resident_kernels.cu"
    , 0};
#pragma GCC diagnostic pop

using Clock = std::chrono::steady_clock;
double ms_since(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }
constexpr int TPB = RF_TPB, CG_NBLK = 512, CG_TPB = 256, BATCH = 4;
unsigned blocks(std::size_t n, unsigned b = TPB) { return unsigned((n + b - 1) / b); }
struct Ctl {
    double rz, rr, rnorm, alpha, beta, pap, target, bnorm, tnorm, spare;
    std::int32_t it, done, error, maxit, corr, pad0, pad1, pad2;
};
static_assert(sizeof(Ctl) == 112, "control block layout");
struct Wd { double hi, lo; };
static_assert(sizeof(RfCtx) <= 4000, "kernel parameter size");
constexpr std::uint64_t state_tag_base = 0x484c333153000000ull, midpoint_tag_base = 0x484c33314d000000ull;
constexpr int state_arrays_max = 14, admit_rows = 14;
const char* const state_array_names[state_arrays_max] = {"E", "B", "E_low", "B_low", "free_net_hi", "free_net_lo", "free_gross",
                                           "conductor_hi", "conductor_lo", "exterior_hi", "exterior_lo", "gauss_admitted_hi",
                                           "gauss_admitted_lo", "residue"};
using Leaf = std::array<unsigned char, 32>;

const char* resident_rules =
    "hall-v3.2-resident-field-step-v1: admitted Gauss ledger words per slot (state), L1 = L0 + (delta - c) on the device "
    "with the Model::prepare word operations, absolute Gauss admission |D - L| <= bound; otherwise "
    "hall-v3.1-resident-field-step-v1: accepted and candidate field states device-resident (E/B words, charge ledger "
    "structure-of-arrays); Model::prepare arithmetic on the device in the order of the sealed host loops (per-element "
    "gathers in ascending scatter order, Words two-sum/fma operations, --fmad=false); gpu1 Jacobi CG on the device "
    "right-hand side; independent FP64 residual ||b-Ax||_2<=max(1e-30,2e-13||b||_2) on the host with the sealed apply "
    "arithmetic and the sealed long-double norm; guard quantities (absolute Gauss, divB, continuity, Gauss increment, "
    "energy sums) in double-word arithmetic with fixed reduction trees; identities: block-Merkle hash lists "
    "(2048-byte leaves hashed on the device, root and a random leaf sample re-hashed on the host)";
}  // namespace

// ------------------------------------------------------------------ host helpers
HostApply::HostApply(const Mesh& mesh) : g(mesh), sqrtM(mesh.ne), E(mesh.ne), B(mesh.nb), K(mesh.ne) {
    for (int e = 0; e < g.ne; ++e) sqrtM[e] = std::sqrt(g.M[e]);
}
// Same arithmetic as the apply lambda of Model::prepare: E = v/sqrtM (0 on PEC edges),
// K = curl_adjoint(curl(E)) (face sums in coefficient order, transpose scattered in face
// order), result = v + ((.5 dt (S+Y)) E + (.25 dt dt) K) / sqrtM (v on PEC edges).
void HostApply::apply(const Vec& S, double dt, const Vec& v, Vec& out) {
    for (int e = 0; e < g.ne; ++e) E[e] = g.fixed_E[e] ? 0. : v[e] / sqrtM[e];
    for (int b = 0; b < g.nb; ++b) {
        double s = 0;
        for (int k = 0; k < 4; ++k) s += g.C[b].coefficient[k] * E[g.C[b].index[k]];
        B[b] = s * g.H[b];
    }
    std::fill(K.begin(), K.end(), 0.);
    for (int b = 0; b < g.nb; ++b) for (int k = 0; k < 4; ++k) K[g.C[b].index[k]] += g.C[b].coefficient[k] * B[b];
    out.resize(g.ne);
    for (int e = 0; e < g.ne; ++e) out[e] = g.fixed_E[e] ? v[e] : v[e] + (.5 * dt * (S[e] + g.Y[e]) * E[e] + .25 * dt * dt * K[e]) / sqrtM[e];
}
void conduction_metal(const Mesh& g, const std::vector<int>& metal_cells, const Vec& sigma, Vec& S) {
    need(sigma.size() == g.cells.size(), "cell material binding mismatch");
    for (std::size_t c = 0; c < g.cells.size(); ++c)
        need(std::isfinite(sigma[c]) && sigma[c] >= 0 && (g.cells[c].body >= 0 || sigma[c] == 0), "nonfinite or vacuum conductivity");
    S.assign(g.ne, 0.);
    for (int c : metal_cells) for (int e : g.cells[c].edges) S[e] += .25 * g.cells[c].volume * sigma[c];
}

// ------------------------------------------------------------------ engine
struct DevSlot {
    hall_gpu::DeviceBuffer E, El, B, Bl, fnet_hi, fnet_lo, fgross, cond_hi, cond_lo, ext_hi, ext_lo, residue, gl_hi, gl_lo;
};
struct ResidentField::Impl {
    std::shared_ptr<Device> device;
    hall_gpu::Session& s;
    const Model& model;
    const Mesh& g;
    ResidentOptions opt;
    hall_gpu::Module mod_rf, mod_cg;
    std::map<std::string, hall_gpu::Kernel> k;
    std::string descriptor;
    // static
    hall_gpu::DeviceBuffer fe, fc, H, ctp, ctf, ctc, pc2, sqrtM, M, Y, fixed, ea, eb, elen, ecp, ecc, qv, ntp, nte, nts, cedge, cvol,
        feed, bg, cond_mask, ext_mask, cell_body, region_grad, VN, metal_idx;
    // slots and work
    DevSlot slot[2];
    int acc = 0;
    hall_gpu::DeviceBuffer sigma, plasma, zero_nn, zero_ne, after_hi, after_lo, after_gross, support, bank_ch, bank_cl, bank_eh, bank_el, bank_g, bank_r;
    hall_gpu::DeviceBuffer S, a, diag, rhs, x, HBh, HBl, midE, midB, HmidB, cmid, solid, short_cur, metal, exterior, total, disp;
    hall_gpu::DeviceBuffer div_hi, div_lo, Q0_hi, Q0_lo, Q1_hi, Q1_lo, md_hi, md_lo, od_hi, od_lo, joule_cell, joule_metal, cdef_hi, cdef_lo;
    hall_gpu::DeviceBuffer eterm, fterm, nval, cval, code_cont, code_joule, c_edge_canon, c_edge_pec, c_face, c_free, c_node, c_cell, c_diff, c_diff2;
    hall_gpu::DeviceBuffer inc[3];
    int inc0 = 0, inc1 = 1, inc2 = 2, history = 0;
    hall_gpu::DeviceBuffer prevE, prevB, chkE, chkB, part, res, leaves_dev, leaves_chk, sample_ids, sample_stage;
    // CG
    hall_gpu::DeviceBuffer r, rt, z, p, Ap, Bh, part1, part2, ctl, ha, hb, hr;
    // pinned
    hall_gpu::HostBuffer stage, ctl_host, res_host, rhs_host, x_host, joule_host, leaves_host, chk_host, samples_host, ids_host, inputs_host;
    // host caches
    Vec sigma_cached, S_host, cell_joule;
    std::vector<int> metal_cells;
    double dt_cached = -1;
    bool sigma_valid = false;
    RfCtx ctx{};
    HostApply host_apply;
    Vec Ax_host, r_host;
    // leaves: per physical slot (state) and the candidate midpoint
    std::vector<Leaf> slot_leaves[2];
    bool slot_leaves_valid[2] = {false, false};
    bool slot_residue[2] = {false, false};
    std::vector<Leaf> mid_leaves;
    bool mid_leaves_valid = false;
    std::mt19937_64 rng;
    DigestStats dstats;
    std::vector<int> h_ntp, h_nte;   // host node -> edges (ascending)
    Vec rq_acc, rq_cand;            // accumulated residual-charge bound (accepted, candidate)
    hall_gpu::DeviceBuffer ginc_acc[2]; int ginc_a = 0;   // accumulated per-node |Gauss increment| (diagnostic)
    std::map<std::string, double> device_ms;
    int maxit = 200;

    Impl(std::shared_ptr<Device> d, const Model& m, ResidentOptions o)
        : device(std::move(d)), s(device->session()), model(m), g(m.mesh()), opt(o), host_apply(m.mesh()) {
        rng.seed(o.sample_seed ? o.sample_seed : std::random_device{}());
        hall_gpu::CompileRequest rq;
        rq.main = {"resident_kernels.cu", kKernels};
        rq.headers = {{"resident_abi.h", kAbi}};
        mod_rf = s.compile_and_load(rq);
        hall_gpu::CompileRequest cq;
        cq.main = {"hall_v2_field_cg.cu", cg_kernel_source};
        mod_cg = s.compile_and_load(cq);
        for (const char* n : {"k_operator", "k_face_rhs", "k_rhs", "k_update_E", "k_update_B", "k_currents", "k_continuity", "k_admit_nodes", "k_charge",
                              "k_admit_edges", "k_admit_cells", "k_ledger", "k_gauss_increment", "k_cells", "k_reduce1", "k_reduce2", "k_differ",
                              "k_midpoint", "k_rounded", "k_copy_f64", "k_add_f64", "k_gather_f64", "k_omega_p", "k_leaves", "k_gather_blocks"})
            k[n] = s.kernel(mod_rf, n);
        for (const char* n : {"k_copy", "k_dot_self", "k_fin_target", "k_init", "k_fin_init", "k_face", "k_edge", "k_fin_alpha", "k_update",
                              "k_fin_beta", "k_direction", "k_residual", "k_fin_true"})
            k[std::string("cg.") + n] = s.kernel(mod_cg, n);
        Digest dg("hall-v3.2-resident-field-descriptor-v1");
        dg.text(resident_rules);
        dg.text(bytes_sha256(kKernels)); dg.text(bytes_sha256(kAbi)); dg.text(bytes_sha256(cg_kernel_source));
        for (const auto* rec : {&mod_rf.record, &mod_cg.record}) {
            dg.u64(rec->options.size()); for (const auto& x : rec->options) dg.text(x);
            dg.i64(rec->target_cc); dg.text(rec->image_kind); dg.i64(rec->nvrtc_major); dg.i64(rec->nvrtc_minor);
        }
        dg.u64(RF_BLOCK_BYTES); dg.u64(CG_NBLK); dg.u64(CG_TPB); dg.u64(BATCH); dg.u64(RF_RBLK);
        descriptor = dg.hex();
        build();
    }
    hall_gpu::DeviceBuffer alloc(std::size_t bytes, const char* tag) { return s.alloc(std::max<std::size_t>(bytes, 16), tag); }
    template <class T> hall_gpu::DeviceBuffer upload_static(const std::vector<T>& v, const char* tag) {
        auto d = alloc(v.size() * sizeof(T), tag);
        if (!v.empty()) {
            need(v.size() * sizeof(T) <= stage.bytes, "resident staging too small");
            std::memcpy(stage.pointer, v.data(), v.size() * sizeof(T));
            s.upload(d, 0, stage, 0, v.size() * sizeof(T));
            s.synchronize("resident static upload");
        }
        return d;
    }
    // Asynchronous upload through the step-input staging region (offset in doubles). The host
    // writes the region only at the start of a step, after the previous step's final synchronize.
    void upload_async(const hall_gpu::DeviceBuffer& d, std::size_t offset_doubles, const double* src, std::size_t count) {
        need((offset_doubles + count) * 8 <= inputs_host.bytes, "step input staging too small");
        std::memcpy(inputs_host.as<double>() + offset_doubles, src, count * 8);
        s.upload(d, 0, inputs_host, offset_doubles * 8, count * 8);
    }
    void upload_words_async(const hall_gpu::DeviceBuffer& hi, const hall_gpu::DeviceBuffer& lo, std::size_t offset_doubles, const WordVec& w) {
        need((offset_doubles + 2 * w.size()) * 8 <= inputs_host.bytes, "step input staging too small");
        double* h = inputs_host.as<double>() + offset_doubles; double* l = h + w.size();
        for (std::size_t i = 0; i < w.size(); ++i) { h[i] = w[i].hi; l[i] = w[i].lo; }
        s.upload(hi, 0, inputs_host, offset_doubles * 8, w.size() * 8);
        s.upload(lo, 0, inputs_host, (offset_doubles + w.size()) * 8, w.size() * 8);
    }
    void upload_now(const hall_gpu::DeviceBuffer& d, const void* src, std::size_t bytes) {
        need(bytes <= stage.bytes, "resident staging too small");
        std::memcpy(stage.pointer, src, bytes);
        s.upload(d, 0, stage, 0, bytes);
        s.synchronize("resident upload");
    }
    void build() {
        const int ne = g.ne, nb = g.nb, nn = g.nn, nc = int(g.cells.size());
        std::size_t big = std::size_t(std::max({ne, nb, nn, nc})) * 16 * 4;
        stage = s.alloc_pinned(std::max<std::size_t>(big, std::size_t(12) * nc * 4 + 64), "resident staging");
        // operator structure (same construction as the gpu1 backend)
        std::vector<std::int32_t> vfe(4 * std::size_t(nb)), vctp(ne + 1, 0), vctf, count(ne + 1, 0);
        Vec vfc(4 * std::size_t(nb)), vH(nb), vctc, vpc2, vsq(ne);
        std::vector<std::uint8_t> vfixed(ne);
        for (int b = 0; b < nb; ++b) {
            need(g.C[b].count == 4, "axial curl rows have four entries");
            vH[b] = g.H[b];
            for (int q = 0; q < 4; ++q) { vfe[4 * b + q] = g.C[b].index[q]; vfc[4 * b + q] = g.C[b].coefficient[q]; ++count[g.C[b].index[q] + 1]; }
        }
        for (int e = 0; e < ne; ++e) { vctp[e + 1] = vctp[e] + count[e + 1]; vsq[e] = std::sqrt(g.M[e]); vfixed[e] = g.fixed_E[e] ? 1 : 0; }
        vctf.resize(vctp[ne]); vctc.resize(vctp[ne]); vpc2.resize(vctp[ne]);
        std::vector<std::int32_t> fill(vctp.begin(), vctp.end() - 1);
        for (int b = 0; b < nb; ++b) for (int q = 0; q < 4; ++q) {
            const int e = g.C[b].index[q];
            vctf[fill[e]] = b; vctc[fill[e]] = g.C[b].coefficient[q]; vpc2[fill[e]] = std::pow(g.C[b].coefficient[q], 2); ++fill[e];
        }
        fe = upload_static(vfe, "fe"); fc = upload_static(vfc, "fc"); H = upload_static(vH, "H"); ctp = upload_static(vctp, "ctp");
        ctf = upload_static(vctf, "ctf"); ctc = upload_static(vctc, "ctc"); pc2 = upload_static(vpc2, "pc2"); sqrtM = upload_static(vsq, "sqrtM");
        fixed = upload_static(vfixed, "fixed"); M = upload_static(g.M, "M"); Y = upload_static(g.Y, "Y");
        std::vector<std::int32_t> va(ne), vb(ne); Vec vlen(ne);
        for (int e = 0; e < ne; ++e) { va[e] = g.edges[e].a; vb[e] = g.edges[e].b; vlen[e] = g.edges[e].length; }
        ea = upload_static(va, "ea"); eb = upload_static(vb, "eb"); elen = upload_static(vlen, "elen");
        // edge -> cells (ascending), quarter volumes
        std::vector<std::int32_t> vecp(ne + 1, 0), vecc, cnt(ne + 1, 0), vcedge(12 * std::size_t(nc)), vcbody(nc);
        Vec vqv(nc), vcvol(nc);
        for (int c = 0; c < nc; ++c) {
            vqv[c] = .25 * g.cells[c].volume; vcvol[c] = g.cells[c].volume; vcbody[c] = g.cells[c].body;
            for (int q = 0; q < 12; ++q) { vcedge[12 * c + q] = g.cells[c].edges[q]; ++cnt[g.cells[c].edges[q] + 1]; }
            if (g.cells[c].body >= 0) metal_cells.push_back(c);
        }
        for (int e = 0; e < ne; ++e) vecp[e + 1] = vecp[e] + cnt[e + 1];
        vecc.resize(vecp[ne]);
        std::vector<std::int32_t> f2(vecp.begin(), vecp.end() - 1);
        for (int c = 0; c < nc; ++c) for (int q = 0; q < 12; ++q) vecc[f2[g.cells[c].edges[q]]++] = c;
        ecp = upload_static(vecp, "ecp"); ecc = upload_static(vecc, "ecc"); qv = upload_static(vqv, "qv"); cvol = upload_static(vcvol, "cvol");
        cedge = upload_static(vcedge, "cedge"); cell_body = upload_static(vcbody, "cell body");
        metal_idx = upload_static(metal_cells, "metal cells");
        // node -> edges (ascending edge index, sign -1 for the a end, +1 for the b end)
        std::vector<std::int32_t> vntp(nn + 1, 0), vnte, vnts, nc2(nn + 1, 0);
        for (int e = 0; e < ne; ++e) { ++nc2[g.edges[e].a + 1]; ++nc2[g.edges[e].b + 1]; }
        for (int n = 0; n < nn; ++n) vntp[n + 1] = vntp[n] + nc2[n + 1];
        vnte.resize(vntp[nn]); vnts.resize(vntp[nn]);
        std::vector<std::int32_t> f3(vntp.begin(), vntp.end() - 1);
        for (int e = 0; e < ne; ++e) {
            need(g.edges[e].a != g.edges[e].b, "degenerate edge");
            vnte[f3[g.edges[e].a]] = e; vnts[f3[g.edges[e].a]++] = -1;
            vnte[f3[g.edges[e].b]] = e; vnts[f3[g.edges[e].b]++] = +1;
        }
        ntp = upload_static(vntp, "ntp"); nte = upload_static(vnte, "nte"); nts = upload_static(vnts, "nts");
        h_ntp.assign(vntp.begin(), vntp.end()); h_nte.assign(vnte.begin(), vnte.end()); rq_acc.assign(nn, 0.); rq_cand.assign(nn, 0.);
        for (auto& b : ginc_acc) { b = alloc(std::size_t(nn) * 8, "gauss increment accumulation"); s.zero(b, 0, b.bytes); }
        const auto& init = model.initial_fields();
        feed = upload_static(init.feed, "feed"); bg = upload_static(init.background_force, "background force");
        std::vector<std::uint8_t> vcm(nn), vem(nn);
        for (int n = 0; n < nn; ++n) { vcm[n] = (g.node_body[n] >= 0 || g.short_node[n]) ? 1 : 0; vem[n] = model.exterior_support()[n] ? 1 : 0; }
        cond_mask = upload_static(vcm, "conductor support"); ext_mask = upload_static(vem, "exterior support");
        {   // upper-region membership gradient (Mesh::upper_region_inflow)
            const int cut = g.kshort + g.nz_foil / 2 - 1;
            Vec membership(nn);
            for (int kk = cut + 1; kk <= g.kfeed; ++kk) for (int j = 0; j <= g.ny; ++j) if (j * g.dy - g.side / 2 > 0) for (int i = 0; i <= g.nx; ++i) membership[g.node(i, j, kk)] = 1;
            region_grad = upload_static(g.gradient(membership), "region gradient");
        }
        VN = upload_static(g.VN, "node volumes");
        // slots
        const std::size_t vne = std::size_t(ne) * 8, vnb = std::size_t(nb) * 8, vnn = std::size_t(nn) * 8;
        for (auto& z : slot) {
            z.E = alloc(vne, "E"); z.El = alloc(vne, "E_low"); z.B = alloc(vnb, "B"); z.Bl = alloc(vnb, "B_low");
            for (auto* b : {&z.fnet_hi, &z.fnet_lo, &z.fgross, &z.cond_hi, &z.cond_lo, &z.ext_hi, &z.ext_lo, &z.residue, &z.gl_hi, &z.gl_lo}) { *b = alloc(vnn, "ledger"); s.zero(*b, 0, b->bytes); }
        }
        sigma = alloc(std::size_t(nc) * 8, "sigma"); plasma = alloc(vne, "plasma"); zero_ne = alloc(vne, "zero ne"); zero_nn = alloc(vnn, "zero nn");
        s.zero(zero_ne, 0, zero_ne.bytes); s.zero(zero_nn, 0, zero_nn.bytes);
        for (auto* b : {&after_hi, &after_lo, &after_gross, &support, &bank_ch, &bank_cl, &bank_eh, &bank_el, &bank_g, &bank_r}) *b = alloc(vnn, "step input");
        for (auto* b : {&S, &a, &diag, &rhs, &x, &midE, &solid, &short_cur, &metal, &exterior, &total, &disp, &prevE, &chkE}) *b = alloc(vne, "edge work");
        for (auto* b : {&HBh, &HBl, &midB, &HmidB, &cmid, &prevB, &chkB}) *b = alloc(vnb, "face work");
        for (auto* b : {&div_hi, &div_lo, &Q0_hi, &Q0_lo, &Q1_hi, &Q1_lo, &md_hi, &md_lo, &od_hi, &od_lo, &cdef_hi, &cdef_lo}) *b = alloc(vnn, "node work");
        for (auto& b : inc) { b = alloc(vne, "increment history"); s.zero(b, 0, b.bytes); }
        joule_cell = alloc(std::size_t(nc) * 8, "cell joule"); joule_metal = alloc(metal_cells.size() * 8, "metal joule");
        eterm = alloc(vne * ET_COUNT, "edge terms"); fterm = alloc(vnb * FT_COUNT, "face terms");
        nval = alloc(vnn * NV_COUNT, "node values"); cval = alloc(std::size_t(nc) * 8 * CV_COUNT, "cell values");
        code_cont = alloc(std::size_t(nn) * 4, "codes"); code_joule = alloc(std::size_t(nc) * 4, "codes");
        c_edge_canon = alloc(std::size_t(ne) * 4, "codes"); c_edge_pec = alloc(std::size_t(ne) * 4, "codes"); c_face = alloc(std::size_t(nb) * 4, "codes");
        c_free = alloc(std::size_t(nn) * 4, "codes"); c_node = alloc(std::size_t(nn) * 4, "codes"); c_cell = alloc(std::size_t(nc) * 4, "codes");
        c_diff = alloc(std::size_t(std::max(ne, nb)) * 4, "codes"); c_diff2 = alloc(std::size_t(std::max(ne, nb)) * 4, "codes");
        part = alloc(std::size_t(RF_ROWS) * RF_RBLK * 16, "reduction partials"); res = alloc(std::size_t(4) * RF_ROWS * 16, "reduction results");
        const std::size_t max_leaves = std::size_t(6) * (vne / RF_BLOCK_BYTES + 2) + std::size_t(10) * (vnn / RF_BLOCK_BYTES + 2) + 16;   // state + midpoint
        leaves_dev = alloc(max_leaves * 32, "leaves"); leaves_chk = alloc(2 * max_leaves * 32, "leaves check");
        sample_ids = alloc(std::size_t(4096) * 4, "sample ids"); sample_stage = alloc(std::size_t(4096) * RF_BLOCK_BYTES / 16, "sample staging");
        // CG work (gpu1 layout)
        for (auto* b : {&r, &rt, &z, &p, &Ap}) *b = alloc(vne, "cg vector");
        Bh = alloc(vnb, "cg face work"); part1 = alloc(CG_NBLK * 8, "cg partials"); part2 = alloc(CG_NBLK * 8, "cg partials"); ctl = alloc(sizeof(Ctl), "cg control");
        ha = alloc(std::size_t(3) * maxit * 8, "cg alpha"); hb = alloc(std::size_t(3) * maxit * 8, "cg beta"); hr = alloc(std::size_t(3) * (maxit + 1) * 8, "cg rnorm");
        s.zero(ha, 0, ha.bytes); s.zero(hb, 0, hb.bytes); s.zero(hr, 0, hr.bytes);
        ctl_host = s.alloc_pinned(sizeof(Ctl), "cg control");
        res_host = s.alloc_pinned(std::size_t(4) * RF_ROWS * 16, "results");
        rhs_host = s.alloc_pinned(vne, "rhs"); x_host = s.alloc_pinned(vne, "x");
        joule_host = s.alloc_pinned(std::max<std::size_t>(metal_cells.size() * 8, 16), "joule");
        leaves_host = s.alloc_pinned(max_leaves * 32, "leaves"); chk_host = s.alloc_pinned(2 * max_leaves * 32, "leaves check");
        samples_host = s.alloc_pinned(std::size_t(4096) * RF_BLOCK_BYTES / 16, "samples"); ids_host = s.alloc_pinned(std::size_t(4096) * 4, "sample ids");
        inputs_host = s.alloc_pinned((std::size_t(nc) + std::size_t(ne) + std::size_t(10) * nn) * 8, "step inputs");
        s.synchronize("resident buffers");
        cell_joule.assign(nc, 0.);
        // context: static part
        RfCtx& c = ctx;
        c.ne = ne; c.nb = nb; c.nn = nn; c.ncell = nc; c.nx = g.nx; c.ny = g.ny; c.nz = g.nz; c.kshort = g.kshort; c.kfeed = g.kfeed;
        c.dx = g.dx; c.dy = g.dy; c.dz = g.dz; c.initial_gauss_bound = model.initial_gauss_solve_bound_C();
        c.fe = fe.address; c.fc = fc.address; c.H = H.address; c.ctp = ctp.address; c.ctf = ctf.address; c.ctc = ctc.address; c.pc2 = pc2.address;
        c.sqrtM = sqrtM.address; c.M = M.address; c.Y = Y.address; c.fixed = fixed.address; c.ea = ea.address; c.eb = eb.address; c.elen = elen.address;
        c.ecp = ecp.address; c.ecc = ecc.address; c.qv = qv.address; c.ntp = ntp.address; c.nte = nte.address; c.nts = nts.address;
        c.cedge = cedge.address; c.cvol = cvol.address; c.feed = feed.address; c.bg = bg.address; c.node_body_or_short = cond_mask.address;
        c.exterior_support = ext_mask.address; c.cell_body = cell_body.address; c.region_grad = region_grad.address;
        c.sigma = sigma.address; c.plasma = plasma.address; c.support = support.address;
        c.bank_cond_hi = bank_ch.address; c.bank_cond_lo = bank_cl.address; c.bank_ext_hi = bank_eh.address; c.bank_ext_lo = bank_el.address;
        c.bank_gross = bank_g.address; c.bank_residue = bank_r.address;
        c.S = S.address; c.a = a.address; c.diag = diag.address; c.rhs = rhs.address; c.x = x.address;
        c.HBh = HBh.address; c.HBl = HBl.address; c.midE = midE.address; c.midB = midB.address; c.HmidB = HmidB.address; c.cmid = cmid.address;
        c.solid = solid.address; c.short_cur = short_cur.address; c.metal = metal.address; c.exterior = exterior.address; c.total = total.address; c.disp = disp.address;
        c.cdef_hi = cdef_hi.address; c.cdef_lo = cdef_lo.address;
        c.div_hi = div_hi.address; c.div_lo = div_lo.address; c.Q0_hi = Q0_hi.address; c.Q0_lo = Q0_lo.address; c.Q1_hi = Q1_hi.address; c.Q1_lo = Q1_lo.address;
        c.joule_cell = joule_cell.address; c.eterm = eterm.address; c.fterm = fterm.address; c.nval = nval.address; c.cval = cval.address;
        c.code_cont = code_cont.address; c.code_joule = code_joule.address; c.metal_cells = metal_idx.address; c.joule_metal = joule_metal.address;
        c.nmetal = int(metal_cells.size());
        S_host.assign(ne, 0.);
    }
    RfSlot rf(const DevSlot& z) const {
        RfSlot r;
        r.E = z.E.address; r.El = z.El.address; r.B = z.B.address; r.Bl = z.Bl.address; r.fnet_hi = z.fnet_hi.address; r.fnet_lo = z.fnet_lo.address;
        r.fgross = z.fgross.address; r.cond_hi = z.cond_hi.address; r.cond_lo = z.cond_lo.address; r.ext_hi = z.ext_hi.address; r.ext_lo = z.ext_lo.address;
        r.residue = z.residue.address; r.gl_hi = z.gl_hi.address; r.gl_lo = z.gl_lo.address;
        return r;
    }
    DevSlot& phys(int role) { return slot[acc ^ role]; }
    void refresh_ctx() {
        ctx.old_ = rf(phys(0)); ctx.new_ = rf(phys(1));
        ctx.inc0 = inc[inc0].address; ctx.inc1 = inc[inc1].address; ctx.inc2 = inc[inc2].address; ctx.history = history;
    }
    void launch(const char* name, hall_gpu::Dim3 grid, hall_gpu::Dim3 block, std::vector<void*> args) {
        s.launch(k.at(name), grid, block, 0, args.data());
    }
    void launch_ctx(const char* name, unsigned n, std::vector<void*> extra = {}) {
        std::vector<void*> args{&ctx};
        args.insert(args.end(), extra.begin(), extra.end());
        launch(name, {blocks(n)}, {TPB}, args);
    }
    // Reduction of rows into res[row_base + row].
    void reduce(const RfRows& rows, int row_base) {
        RfRows rr = rows;
        std::uint64_t pa = part.address, ra = res.address + std::uint64_t(row_base) * 16;
        launch("k_reduce1", {RF_RBLK, unsigned(rr.count)}, {TPB}, {&rr, &pa});
        launch("k_reduce2", {unsigned(rr.count)}, {TPB}, {&rr, &pa, &ra});
    }
    static void row(RfRows& r, std::uint64_t ptr, int n, int kind) {
        need(r.count < RF_ROWS, "too many reduction rows");
        r.ptr[r.count] = ptr; r.n[r.count] = n; r.kind[r.count] = kind; ++r.count;
    }
    std::uint64_t etp(int t) const { return eterm.address + std::uint64_t(t) * g.ne * 8; }
    std::uint64_t ftp(int t) const { return fterm.address + std::uint64_t(t) * g.nb * 8; }
    std::uint64_t nvp(int t) const { return nval.address + std::uint64_t(t) * g.nn * 8; }
    std::uint64_t cvp(int t) const { return cval.address + std::uint64_t(t) * g.cells.size() * 8; }
    // Admission kernels of a role (0 old, 1 new) and their reduction rows (admit_rows = 14 rows from row_base):
    // 0-5 first failing index (edge canon, face canon, free, PEC, node, cell), 6 max |D-L|, 7 max bound, 8 max |divB|,
    // 9 max divB bound, 10 max |D|, 11 max |L|, 12 max |D-L|/bound, 13 max |D|/(bound+|L|).
    void enqueue_admit(int role, int row_base, std::uint64_t Qh, std::uint64_t Ql) {
        int sl = role;
        std::uint64_t cn = c_node.address, cf = c_free.address, ec = c_edge_canon.address, ep = c_edge_pec.address, fcd = c_face.address, cc = c_cell.address;
        std::uint64_t dv = nvp(role ? NV_DEF1 : NV_DEF0), bv = nvp(role ? NV_BND1 : NV_BND0), av = nvp(role ? NV_ACT1 : NV_ACT0);
        std::uint64_t rv = nvp(role ? NV_RAW1 : NV_RAW0), lv = nvp(role ? NV_LED1 : NV_LED0), qv2 = nvp(role ? NV_RAT1 : NV_RAT0), pv = nvp(role ? NV_APR1 : NV_APR0);
        std::uint64_t cv = cvp(role ? CV_DIV1 : CV_DIV0), cb = cvp(role ? CV_DBND1 : CV_DBND0);
        launch_ctx("k_admit_nodes", g.nn, {&sl, &Qh, &Ql, &cn, &cf, &dv, &bv, &av, &rv, &lv, &qv2, &pv});
        launch_ctx("k_admit_edges", std::max(g.ne, g.nb), {&sl, &ec, &ep, &fcd});
        launch_ctx("k_admit_cells", int(g.cells.size()), {&sl, &cc, &cv, &cb});
        RfRows rr{};
        row(rr, ec, g.ne, 2); row(rr, fcd, g.nb, 2); row(rr, cf, g.nn, 2); row(rr, ep, g.ne, 2); row(rr, cn, g.nn, 2); row(rr, cc, int(g.cells.size()), 2);
        row(rr, dv, g.nn, 1); row(rr, bv, g.nn, 1); row(rr, cv, int(g.cells.size()), 1); row(rr, cb, int(g.cells.size()), 1);
        row(rr, rv, g.nn, 1); row(rr, lv, g.nn, 1); row(rr, qv2, g.nn, 1); row(rr, pv, g.nn, 1);
        need(rr.count == admit_rows, "admission row layout");
        reduce(rr, row_base);
    }
    const Wd* results() const { return res_host.as<Wd>(); }
    static std::uint32_t idx(const Wd& w) { return std::uint32_t(w.hi); }
    // Host decoding of an admission (10 rows from base) in the sealed order.
    Constraints decode_admit(int role, int base, const SlotMeta& meta) {
        const Wd* R = results() + base;
        need(std::isfinite(meta.time_s) && meta.time_s >= 0 && meta.revision < UINT64_MAX, "field clock/revision");
        if (idx(R[0]) != 0xffffffffu || idx(R[1]) != 0xffffffffu) throw std::runtime_error("noncanonical field words");
        if (idx(R[2]) != 0xffffffffu) throw std::runtime_error("free charge/gross operand invalid");
        need(digest(meta.conductor_declaration) && digest(meta.exterior_declaration), "explicit conductor/exterior charge declarations required");
        if (idx(R[3]) != 0xffffffffu) throw std::runtime_error("nonzero PEC tangential field");
        const std::uint32_t n0 = idx(R[4]);
        if (n0 != 0xffffffffu) {
            std::vector<std::uint32_t> codes(g.nn);
            download_now(c_node, codes.data(), codes.size() * 4);
            const std::uint32_t code = codes.at(n0);
            if (code == 2) throw std::runtime_error("noncanonical owned charge");
            if (code == 3) throw std::runtime_error("conductor charge outside declared material/short");
            if (code == 4) throw std::runtime_error("exterior charge outside the EM boundary and the kinetic escape planes");
            if (code == 6) throw std::runtime_error("noncanonical admitted Gauss ledger");
            // |D - L| (the checked value); the signed values for the refusal text come from the slot words
            double dv = 0, bv = 0, av = 0, raw = 0, lh = 0, ll = 0, qh = 0, ql = 0;
            download_now_at(nval, (std::size_t(role ? NV_DEF1 : NV_DEF0) * g.nn + n0) * 8, &dv, 8);
            download_now_at(nval, (std::size_t(role ? NV_BND1 : NV_BND0) * g.nn + n0) * 8, &bv, 8);
            download_now_at(nval, (std::size_t(role ? NV_ACT1 : NV_ACT0) * g.nn + n0) * 8, &av, 8);
            download_now_at(nval, (std::size_t(role ? NV_RAW1 : NV_RAW0) * g.nn + n0) * 8, &raw, 8);
            const DevSlot& z = phys(role);
            download_now_at(z.gl_hi, std::size_t(n0) * 8, &lh, 8); download_now_at(z.gl_lo, std::size_t(n0) * 8, &ll, 8);
            (void)qh; (void)ql;
            std::ostringstream why;
            why << std::setprecision(17) << "absolute Gauss constraint disagrees with actual declared charge: node=" << n0 << " defect_C=" << dv
                << " bound_C=" << bv << " field_activity_C=" << av;
            if (lh != 0 || ll != 0) why << " admitted_gauss_ledger_C=" << double((long double)lh + ll) << " uncorrected_defect_C(abs)=" << raw;
            throw std::runtime_error(why.str());
        }
        if (idx(R[5]) != 0xffffffffu) throw std::runtime_error("absolute divB constraint failed");
        Constraints c;
        c.max_gauss_C = R[6].hi; c.max_gauss_bound_C = R[7].hi; c.max_divB_T_m = R[8].hi; c.max_divB_bound_T_m = R[9].hi;
        c.max_gauss_defect_C = R[10].hi; c.max_gauss_ledger_C = R[11].hi; c.max_gauss_node_ratio = R[12].hi; c.max_gauss_defect_over_bound_plus_ledger = R[13].hi;
        return c;
    }
    void download_now(const hall_gpu::DeviceBuffer& d, void* dst, std::size_t bytes) { download_now_at(d, 0, dst, bytes); }
    void download_now_at(const hall_gpu::DeviceBuffer& d, std::size_t off, void* dst, std::size_t bytes) {
        need(bytes <= stage.bytes, "resident staging too small");
        s.download(stage, 0, d, off, bytes);
        s.synchronize("resident download");
        std::memcpy(dst, stage.pointer, bytes);
    }
    // ---------------------------------------------------------------- CG (gpu1 kernels and control flow)
    void cg_apply(std::uint64_t v, std::uint64_t out, int gated, int with_dot, double q) {
        int n = g.ne, nbb = g.nb;
        auto fea = fe.address, fca = fc.address, Ha = H.address, sq = sqrtM.address, fx = fixed.address, Bha = Bh.address, ct = ctl.address,
             ctpa = ctp.address, ctfa = ctf.address, ctca = ctc.address, aa = a.address, pa = part1.address;
        launch("cg.k_face", {blocks(nbb, CG_TPB)}, {CG_TPB}, {&nbb, &fea, &fca, &Ha, &sq, &fx, &v, &Bha, &ct, &gated});
        launch("cg.k_edge", {CG_NBLK}, {CG_TPB}, {&n, &ctpa, &ctfa, &ctca, &Bha, &sq, &aa, &q, &fx, &v, &out, &pa, &ct, &gated, &with_dot});
    }
    void cg_iteration(double q) {
        int n = g.ne;
        auto xa = x.address, ra = r.address, pa = p.address, Apa = Ap.address, dg = diag.address, za = z.address, p1 = part1.address, p2 = part2.address,
             ct = ctl.address, haa = ha.address, hba = hb.address, hra = hr.address;
        cg_apply(pa, Apa, 1, 1, q);
        launch("cg.k_fin_alpha", {1}, {CG_NBLK}, {&p1, &ct, &haa});
        launch("cg.k_update", {CG_NBLK}, {CG_TPB}, {&n, &xa, &ra, &pa, &Apa, &dg, &za, &p1, &p2, &ct});
        launch("cg.k_fin_beta", {1}, {CG_NBLK}, {&p1, &p2, &ct, &hba, &hra});
        launch("cg.k_direction", {blocks(g.ne, CG_TPB)}, {CG_TPB}, {&n, &za, &pa, &ct});
    }
    void cg_tail(double q) {
        int n = g.ne;
        auto xa = x.address, Apa = Ap.address, ba = rhs.address, rta = rt.address, p1 = part1.address, ct = ctl.address;
        cg_apply(xa, Apa, 0, 0, q);
        launch("cg.k_residual", {CG_NBLK}, {CG_TPB}, {&n, &ba, &Apa, &rta, &p1});
        launch("cg.k_fin_true", {1}, {CG_NBLK}, {&p1, &ct});
    }
    void cg_init_from(std::uint64_t src, int corr) {
        int n = g.ne;
        auto ra = r.address, dg = diag.address, za = z.address, pa = p.address, p1 = part1.address, p2 = part2.address, ct = ctl.address, hra = hr.address;
        launch("cg.k_copy", {CG_NBLK}, {CG_TPB}, {&n, &src, &ra});
        launch("cg.k_init", {CG_NBLK}, {CG_TPB}, {&n, &ra, &dg, &za, &pa, &p1, &p2});
        launch("cg.k_fin_init", {1}, {CG_NBLK}, {&p1, &p2, &ct, &corr, &hra});
    }
    // pending_check: called after the first synchronize (refusals that precede the solve in the sealed order).
    template <class Check> void cg_read(Check&& pending_check) {
        s.download(ctl_host, 0, ctl, 0, sizeof(Ctl));
        s.synchronize("field control readback");
        pending_check();
        if (ctl_host.as<Ctl>()->error) throw std::runtime_error("non-SPD field operator (device CG)");
    }
    template <class Check> void cg_solve(double dt, SolveRecord& rec, Check&& pending_check) {
        const double q = .25 * dt * dt;
        int n = g.ne, mi = maxit, corr0 = 0;
        auto ba = rhs.address, p1 = part1.address, ct = ctl.address;
        launch("cg.k_dot_self", {CG_NBLK}, {CG_TPB}, {&n, &ba, &p1});
        launch("cg.k_fin_target", {1}, {CG_NBLK}, {&p1, &ct, &mi});
        s.zero(x, 0, x.bytes);
        cg_init_from(rhs.address, corr0);
        for (int i = 0; i < BATCH; ++i) cg_iteration(q);
        cg_tail(q);
        cg_read(pending_check);
        rec.iterations = 0; rec.corrections = 0;
        for (int corr = 0; corr < 3; ++corr) {
            rec.corrections = corr + 1;
            bool needs_tail = false;
            if (corr > 0) { cg_init_from(rt.address, corr); ctl_host.as<Ctl>()->done = 0; needs_tail = true; }
            for (int guard = 0; !ctl_host.as<Ctl>()->done; ++guard) {
                need(guard <= maxit / BATCH + 2, "device CG did not stop within maxit");
                for (int i = 0; i < BATCH; ++i) cg_iteration(q);
                cg_read([] {});
                needs_tail = true;
            }
            if (needs_tail) { cg_tail(q); cg_read([] {}); }
            const Ctl& C = *ctl_host.as<Ctl>();
            rec.iterations += C.it; rec.backend_residual = C.tnorm;
            if (C.tnorm <= C.target) return;
        }
        throw std::runtime_error("field solve independent residual failed (device CG)");
    }
    // ---------------------------------------------------------------- leaves
    RfLeafTable state_table(int role, std::uint64_t out) {
        DevSlot& z = phys(role);
        RfLeafTable t{};
        const hall_gpu::DeviceBuffer* arr[state_arrays_max] = {&z.E, &z.B, &z.El, &z.Bl, &z.fnet_hi, &z.fnet_lo, &z.fgross, &z.cond_hi, &z.cond_lo,
                                                               &z.ext_hi, &z.ext_lo, &z.gl_hi, &z.gl_lo, &z.residue};
        const std::uint64_t ne_ = std::uint64_t(g.ne), nb_ = std::uint64_t(g.nb), nn_ = std::uint64_t(g.nn);
        const std::uint64_t counts[state_arrays_max] = {ne_, nb_, ne_, nb_, nn_, nn_, nn_, nn_, nn_, nn_, nn_, nn_, nn_, nn_};
        const int arrays = slot_residue[acc ^ role] ? state_arrays_max : state_arrays_max - 1;
        int L = 0;
        for (int i = 0; i < arrays; ++i) {
            t.ptr[i] = arr[i]->address; t.bytes[i] = counts[i] * 8; t.tag[i] = state_tag_base + std::uint64_t(i); t.first_leaf[i] = L;
            L += int((t.bytes[i] + RF_BLOCK_BYTES - 1) / RF_BLOCK_BYTES);
        }
        t.first_leaf[arrays] = L; t.arrays = arrays; t.out = out;
        return t;
    }
    RfLeafTable midpoint_table(std::uint64_t out, std::uint64_t E, std::uint64_t B) {
        RfLeafTable t{};
        t.ptr[0] = E; t.bytes[0] = std::uint64_t(g.ne) * 8; t.tag[0] = midpoint_tag_base;
        t.ptr[1] = B; t.bytes[1] = std::uint64_t(g.nb) * 8; t.tag[1] = midpoint_tag_base + 1;
        t.first_leaf[0] = 0; t.first_leaf[1] = int((t.bytes[0] + RF_BLOCK_BYTES - 1) / RF_BLOCK_BYTES);
        t.first_leaf[2] = t.first_leaf[1] + int((t.bytes[1] + RF_BLOCK_BYTES - 1) / RF_BLOCK_BYTES);
        t.arrays = 2; t.out = out;
        return t;
    }
    static int leaf_count(const RfLeafTable& t) { return t.first_leaf[t.arrays]; }
    void enqueue_leaves(const RfLeafTable& t) {
        RfLeafTable tt = t;
        launch("k_leaves", {blocks(leaf_count(t), 64)}, {64}, {&tt});
    }
    // Host leaf of (table, leaf index) from a payload buffer.
    static Leaf host_leaf(std::uint64_t tag, std::uint64_t k, const unsigned char* payload, std::uint64_t nbytes) {
        SHA256_CTX c;
        SHA256_Init(&c);
        SHA256_Update(&c, &tag, 8); SHA256_Update(&c, &k, 8); SHA256_Update(&c, &nbytes, 8);
        SHA256_Update(&c, payload, nbytes);
        Leaf out;
        SHA256_Final(out.data(), &c);
        return out;
    }
    static void locate(const RfLeafTable& t, int L, int& array, std::uint64_t& k, std::uint64_t& nbytes) {
        array = 0;
        while (L >= t.first_leaf[array + 1]) ++array;
        k = std::uint64_t(L - t.first_leaf[array]);
        const std::uint64_t begin = k * RF_BLOCK_BYTES;
        nbytes = std::min<std::uint64_t>(t.bytes[array] - begin, RF_BLOCK_BYTES);
    }
    // Sample ids for a table (stored in ids_host from offset), enqueue the gather into samples staging.
    int enqueue_samples(const RfLeafTable& t, int count, int offset) {
        const int L = leaf_count(t);
        count = std::min(count, L);
        std::uniform_int_distribution<int> pick(0, L - 1);
        auto* ids = ids_host.as<std::int32_t>() + offset;
        for (int i = 0; i < count; ++i) ids[i] = pick(rng);
        s.upload(sample_ids, std::size_t(offset) * 4, ids_host, std::size_t(offset) * 4, std::size_t(count) * 4);
        RfLeafTable tt = t;
        std::uint64_t idp = sample_ids.address + std::uint64_t(offset) * 4, stp = sample_stage.address + std::uint64_t(offset) * RF_BLOCK_BYTES;
        int c = count;
        launch("k_gather_blocks", {unsigned(count)}, {128}, {&tt, &c, &idp, &stp});
        s.download(samples_host, std::size_t(offset) * RF_BLOCK_BYTES, sample_stage, std::size_t(offset) * RF_BLOCK_BYTES, std::size_t(count) * RF_BLOCK_BYTES);
        return count;
    }
    bool verify_samples(const RfLeafTable& t, const Leaf* leaves, int count, int offset) {
        const auto* ids = ids_host.as<std::int32_t>() + offset;
        for (int i = 0; i < count; ++i) {
            int arr; std::uint64_t kk, nbytes;
            locate(t, ids[i], arr, kk, nbytes);
            const auto* payload = samples_host.as<unsigned char>() + std::size_t(offset + i) * RF_BLOCK_BYTES;
            if (host_leaf(t.tag[arr], kk, payload, nbytes) != leaves[ids[i]]) return false;
        }
        dstats.samples += count; dstats.leaves_hashed_on_host += count;
        return true;
    }
    const hall_gpu::DeviceBuffer& buffer_of(std::uint64_t address) const {
        for (const auto& z : slot)
            for (const auto* d : {&z.E, &z.El, &z.B, &z.Bl, &z.fnet_hi, &z.fnet_lo, &z.fgross, &z.cond_hi, &z.cond_lo, &z.ext_hi, &z.ext_lo, &z.residue, &z.gl_hi, &z.gl_lo})
                if (d->address == address) return *d;
        for (const auto* d : {&midE, &midB}) if (d->address == address) return *d;
        throw std::runtime_error("leaf table array is not a resident buffer");
    }
    // Every leaf on the host (qualification mode): downloads the arrays of the table.
    bool verify_all_leaves(const RfLeafTable& t, const Leaf* leaves) {
        std::vector<unsigned char> buf;
        for (int a2 = 0; a2 < t.arrays; ++a2) {
            buf.resize(t.bytes[a2]);
            const hall_gpu::DeviceBuffer& view = buffer_of(t.ptr[a2]);
            for (std::size_t off = 0; off < t.bytes[a2]; off += stage.bytes) {
                const std::size_t nb2 = std::min<std::size_t>(stage.bytes, t.bytes[a2] - off);
                s.download(stage, 0, view, off, nb2); s.synchronize("verify download");
                std::memcpy(buf.data() + off, stage.pointer, nb2);
            }
            for (int L = t.first_leaf[a2]; L < t.first_leaf[a2 + 1]; ++L) {
                const std::uint64_t kk = std::uint64_t(L - t.first_leaf[a2]), begin = kk * RF_BLOCK_BYTES;
                const std::uint64_t nbytes = std::min<std::uint64_t>(t.bytes[a2] - begin, RF_BLOCK_BYTES);
                if (host_leaf(t.tag[a2], kk, buf.data() + begin, nbytes) != leaves[L]) return false;
                ++dstats.leaves_hashed_on_host;
            }
        }
        return true;
    }
    std::string state_root(int role, const std::string& api, const std::string& configuration, const SlotMeta& m, const std::vector<Leaf>& leaves) {
        const RfLeafTable t = state_table(role, 0);
        Digest d("hall-v3.2-field-state-merkle-v1");
        d.text(api); d.text(configuration); d.text("hall-v3.1-block-hash-list-v1"); d.u64(RF_BLOCK_BYTES);
        d.text("hall-v2-field-words-v1");
        for (int i = 0; i < t.arrays; ++i) {
            if (i == 4) d.text("hall-v3-charge-ledger-v1");
            if (i == 11) d.text("hall-v3.2-admitted-gauss-ledger-v1");
            d.text(state_array_names[i]); d.u64(t.bytes[i] / 8); d.u64(std::uint64_t(t.first_leaf[i + 1] - t.first_leaf[i]));
            d.bytes(leaves[t.first_leaf[i]].data(), std::size_t(t.first_leaf[i + 1] - t.first_leaf[i]) * 32);
        }
        d.f64(m.time_s); d.u64(m.revision);
        d.text(m.conductor_declaration); d.text(m.exterior_declaration); d.boolean(m.residue_present);
        return d.hex();
    }
    std::string midpoint_root(const std::vector<Leaf>& leaves) {
        const RfLeafTable t = midpoint_table(0, midE.address, midB.address);
        Digest d("hall-v3.1-authoritative-midpoint-merkle-v1");
        d.text("hall-v3.1-block-hash-list-v1"); d.u64(RF_BLOCK_BYTES);
        for (int i = 0; i < 2; ++i) {
            d.text(i ? "B" : "E"); d.u64(t.bytes[i] / 8); d.u64(std::uint64_t(t.first_leaf[i + 1] - t.first_leaf[i]));
            d.bytes(leaves[t.first_leaf[i]].data(), std::size_t(t.first_leaf[i + 1] - t.first_leaf[i]) * 32);
        }
        return d.hex();
    }
};

// ------------------------------------------------------------------ public API
ResidentField::ResidentField(std::shared_ptr<Device> d, const Model& m, ResidentOptions o)
    : p_(std::make_unique<Impl>(std::move(d), m, o)) {}
ResidentField::~ResidentField() = default;
std::string ResidentField::descriptor_text() { return resident_rules; }
const std::string& ResidentField::descriptor_sha256() const noexcept { return p_->descriptor; }
const ResidentOptions& ResidentField::options() const noexcept { return p_->opt; }
std::string ResidentField::compile_json() const {
    return "{\"resident\":" + hall_gpu::compile_json(p_->mod_rf.record) + ",\"cg\":" + hall_gpu::compile_json(p_->mod_cg.record) + "}";
}
std::string ResidentField::memory_json() { return p_->device->memory_json(); }
DigestStats ResidentField::digest_stats() const noexcept { return p_->dstats; }
std::map<std::string, double> ResidentField::device_profile() { return p_->device_ms; }
int ResidentField::history() const noexcept { return p_->history; }
const Vec& ResidentField::residual_charge_cumulative() const noexcept { return p_->rq_acc; }
void ResidentField::gauss_increment_cumulative(Vec& out) {
    auto& I = *p_;
    out.resize(I.g.nn);
    I.download_now(I.ginc_acc[I.ginc_a], out.data(), out.size() * 8);
}
Constraints ResidentField::gauss_profile(int role, const SlotMeta& m, Vec& defect, Vec& bound, Vec* uncorrected, Vec* ledger) {
    auto& I = *p_;
    I.refresh_ctx();
    I.enqueue_admit(role, 0, I.Q0_hi.address, I.Q0_lo.address);
    I.s.download(I.res_host, 0, I.res, 0, std::size_t(RF_ROWS) * 16);
    I.s.synchronize("gauss profile");
    const auto c = I.decode_admit(role, 0, m);
    defect.resize(I.g.nn); bound.resize(I.g.nn);
    I.download_now_at(I.nval, std::size_t(role ? NV_DEF1 : NV_DEF0) * I.g.nn * 8, defect.data(), defect.size() * 8);
    I.download_now_at(I.nval, std::size_t(role ? NV_BND1 : NV_BND0) * I.g.nn * 8, bound.data(), bound.size() * 8);
    if (uncorrected) { uncorrected->resize(I.g.nn); I.download_now_at(I.nval, std::size_t(role ? NV_RAW1 : NV_RAW0) * I.g.nn * 8, uncorrected->data(), I.g.nn * 8); }
    if (ledger) {   // signed rounded ledger hi + lo
        Vec hi(I.g.nn), lo(I.g.nn);
        const DevSlot& z = I.phys(role);
        I.download_now(z.gl_hi, hi.data(), I.g.nn * 8); I.download_now(z.gl_lo, lo.data(), I.g.nn * 8);
        ledger->resize(I.g.nn); for (int n = 0; n < I.g.nn; ++n) (*ledger)[n] = hi[n] + lo[n];
    }
    return c;
}

void ResidentField::upload(int role, const axial::State& st, const ChargeLedger& q) {
    auto& I = *p_;
    const auto& g = I.g;
    need(st.E.size() == std::size_t(g.ne) && st.B.size() == std::size_t(g.nb) && st.E_low.size() == st.E.size() && st.B_low.size() == st.B.size(), "field state shape");
    need(q.particles.net_C.size() == std::size_t(g.nn) && q.particles.gross_C.size() == std::size_t(g.nn) && q.conductor_C.size() == std::size_t(g.nn) &&
         q.exterior_C.size() == std::size_t(g.nn) && (q.residue_C.empty() || q.residue_C.size() == std::size_t(g.nn)) &&
         (q.gauss_admitted_C.empty() || q.gauss_admitted_C.size() == std::size_t(g.nn)), "charge ledger shape");
    DevSlot& z = I.phys(role);
    I.upload_now(z.E, st.E.data(), st.E.size() * 8); I.upload_now(z.El, st.E_low.data(), st.E_low.size() * 8);
    I.upload_now(z.B, st.B.data(), st.B.size() * 8); I.upload_now(z.Bl, st.B_low.data(), st.B_low.size() * 8);
    Vec hi(g.nn), lo(g.nn);
    auto put = [&](const WordVec& w, DevSlot& zz, hall_gpu::DeviceBuffer DevSlot::*h, hall_gpu::DeviceBuffer DevSlot::*l) {
        for (int n = 0; n < g.nn; ++n) { hi[n] = w[n].hi; lo[n] = w[n].lo; }
        I.upload_now(zz.*h, hi.data(), hi.size() * 8); I.upload_now(zz.*l, lo.data(), lo.size() * 8);
    };
    put(q.particles.net_C, z, &DevSlot::fnet_hi, &DevSlot::fnet_lo);
    I.upload_now(z.fgross, q.particles.gross_C.data(), g.nn * 8);
    put(q.conductor_C, z, &DevSlot::cond_hi, &DevSlot::cond_lo);
    put(q.exterior_C, z, &DevSlot::ext_hi, &DevSlot::ext_lo);
    if (q.residue_C.empty()) { I.s.zero(z.residue, 0, z.residue.bytes); I.s.synchronize("residue zero"); }
    else I.upload_now(z.residue, q.residue_C.data(), g.nn * 8);
    if (q.gauss_admitted_C.empty()) { I.s.zero(z.gl_hi, 0, z.gl_hi.bytes); I.s.zero(z.gl_lo, 0, z.gl_lo.bytes); I.s.synchronize("ledger zero"); }
    else put(q.gauss_admitted_C, z, &DevSlot::gl_hi, &DevSlot::gl_lo);
    I.slot_residue[I.acc ^ role] = !q.residue_C.empty();
    I.slot_leaves_valid[I.acc ^ role] = false;
}
void ResidentField::download(int role, axial::State& st, ChargeLedger& q, const SlotMeta& m) {
    auto& I = *p_;
    const auto& g = I.g;
    DevSlot& z = I.phys(role);
    st.E.resize(g.ne); st.E_low.resize(g.ne); st.B.resize(g.nb); st.B_low.resize(g.nb);
    I.download_now(z.E, st.E.data(), g.ne * 8); I.download_now(z.El, st.E_low.data(), g.ne * 8);
    I.download_now(z.B, st.B.data(), g.nb * 8); I.download_now(z.Bl, st.B_low.data(), g.nb * 8);
    st.time_s = m.time_s; st.revision = m.revision;
    Vec hi(g.nn), lo(g.nn);
    auto get = [&](WordVec& w, const hall_gpu::DeviceBuffer& h, const hall_gpu::DeviceBuffer& l) {
        I.download_now(h, hi.data(), g.nn * 8); I.download_now(l, lo.data(), g.nn * 8);
        w.resize(g.nn); for (int n = 0; n < g.nn; ++n) { w[n].hi = hi[n]; w[n].lo = lo[n]; }
    };
    get(q.particles.net_C, z.fnet_hi, z.fnet_lo);
    q.particles.gross_C.resize(g.nn); I.download_now(z.fgross, q.particles.gross_C.data(), g.nn * 8);
    get(q.conductor_C, z.cond_hi, z.cond_lo); get(q.exterior_C, z.ext_hi, z.ext_lo);
    get(q.gauss_admitted_C, z.gl_hi, z.gl_lo);   // v3.2: always nn words
    if (m.residue_present) { q.residue_C.resize(g.nn); I.download_now(z.residue, q.residue_C.data(), g.nn * 8); } else q.residue_C.clear();
    q.conductor_declaration = m.conductor_declaration; q.exterior_declaration = m.exterior_declaration;
}
Constraints ResidentField::admit(int role, const SlotMeta& m, const FreeCharge* census_charge) {
    auto& I = *p_;
    I.refresh_ctx();
    RfSlot keep = role ? I.ctx.new_ : I.ctx.old_;
    hall_gpu::DeviceBuffer th, tl, tg;
    if (census_charge) {   // admission with a census charge in place of the slot's free charge (admit_census)
        const auto& g = I.g;
        Vec hi(g.nn), lo(g.nn);
        for (int n = 0; n < g.nn; ++n) { hi[n] = census_charge->net_C[n].hi; lo[n] = census_charge->net_C[n].lo; }
        I.upload_now(I.after_hi, hi.data(), g.nn * 8); I.upload_now(I.after_lo, lo.data(), g.nn * 8);
        I.upload_now(I.after_gross, census_charge->gross_C.data(), g.nn * 8);
        RfSlot& sl = role ? I.ctx.new_ : I.ctx.old_;
        sl.fnet_hi = I.after_hi.address; sl.fnet_lo = I.after_lo.address; sl.fgross = I.after_gross.address;
    }
    const std::uint64_t Qh = I.Q0_hi.address, Ql = I.Q0_lo.address;
    I.enqueue_admit(role, 0, Qh, Ql);
    I.s.download(I.res_host, 0, I.res, 0, std::size_t(RF_ROWS) * 16);
    I.s.synchronize("resident admission");
    auto c = I.decode_admit(role, 0, m);
    (role ? I.ctx.new_ : I.ctx.old_) = keep;
    return c;
}

ResidentStep ResidentField::prepare(const ResidentInputs& in, const SlotMeta& accepted) {
    auto& I = *p_;
    const auto& g = I.g;
    ResidentStep out;
    auto t0 = Clock::now();
    need(in.sigma && in.sigma->size() == g.cells.size(), "cell material binding mismatch");
    need(std::isfinite(in.dt) && in.dt > 0, "authoritative field interval");
    // per-step inputs: staged once, uploaded asynchronously on the stream
    const std::size_t nc = g.cells.size(), o_sigma = 0, o_plasma = nc, o_after = o_plasma + g.ne, o_support = o_after + 3 * std::size_t(g.nn),
                      o_bank = o_support + g.nn;
    const bool sigma_changed = !I.sigma_valid || *in.sigma != I.sigma_cached;
    if (sigma_changed) {
        conduction_metal(g, I.metal_cells, *in.sigma, I.S_host);   // also checks the sealed conductivity rule
        I.upload_async(I.sigma, o_sigma, in.sigma->data(), nc);
        I.sigma_cached = *in.sigma; I.sigma_valid = true;
    }
    I.refresh_ctx();
    RfCtx& c = I.ctx;
    c.dt = in.dt; c.current = in.current;
    if (in.plasma) { need(in.plasma->size() == std::size_t(g.ne) && axial::finite(*in.plasma), "field current shape/finiteness"); I.upload_async(I.plasma, o_plasma, in.plasma->data(), g.ne); c.plasma = I.plasma.address; }
    else c.plasma = I.zero_ne.address;
    if (in.after) {
        need(in.after->net_C.size() == std::size_t(g.nn) && in.after->gross_C.size() == std::size_t(g.nn), "free charge shape");
        I.upload_words_async(I.after_hi, I.after_lo, o_after, in.after->net_C);
        I.upload_async(I.after_gross, o_after + 2 * std::size_t(g.nn), in.after->gross_C.data(), g.nn);
        c.after_hi = I.after_hi.address; c.after_lo = I.after_lo.address; c.after_gross = I.after_gross.address;
    } else { c.after_hi = c.old_.fnet_hi; c.after_lo = c.old_.fnet_lo; c.after_gross = c.old_.fgross; }
    if (in.support) { need(in.support->size() == std::size_t(g.nn), "deposit support shape"); I.upload_async(I.support, o_support, in.support->data(), g.nn); c.support = I.support.address; }
    else c.support = I.zero_nn.address;
    c.transfer = in.bank && !in.bank->empty ? 1 : 0; c.residue_any = 0;
    bool residue_present = accepted.residue_present;
    if (c.transfer) {
        const auto& b = *in.bank;
        need(b.conductor_C.size() == std::size_t(g.nn) && b.exterior_C.size() == std::size_t(g.nn) && b.gross_C.size() == std::size_t(g.nn), "bank transfer shape");
        I.upload_words_async(I.bank_ch, I.bank_cl, o_bank, b.conductor_C);
        I.upload_words_async(I.bank_eh, I.bank_el, o_bank + 2 * std::size_t(g.nn), b.exterior_C);
        I.upload_async(I.bank_g, o_bank + 4 * std::size_t(g.nn), b.gross_C.data(), g.nn);
        bool any = false; for (double v : b.residue_C) any = any || v != 0;
        if (any) {
            need(b.residue_C.size() == std::size_t(g.nn), "bank residue shape");
            I.upload_async(I.bank_r, o_bank + 5 * std::size_t(g.nn), b.residue_C.data(), g.nn); c.residue_any = 1; residue_present = true;
        }
    }
    out.host_ms["inputs"] = ms_since(t0); t0 = Clock::now();
    // operator (S, a, diagonal) when sigma or dt changed
    if (sigma_changed || in.dt != I.dt_cached) { I.launch_ctx("k_operator", g.ne); I.dt_cached = in.dt; }
    // group A: continuity, old admission
    I.launch_ctx("k_continuity", g.nn);
    I.enqueue_admit(0, 2, I.Q0_hi.address, I.Q0_lo.address);
    {
        RfRows rr{};
        Impl::row(rr, I.code_cont.address, g.nn, 2); Impl::row(rr, I.nvp(NV_CONT), g.nn, 1);
        I.reduce(rr, 0);
    }
    I.s.download(I.res_host, 0, I.res, 0, std::size_t(2 + admit_rows) * 16);
    // right-hand side
    I.launch_ctx("k_face_rhs", g.nb);
    I.launch_ctx("k_rhs", g.ne);
    auto check_A = [&] {
        const Wd* R = I.results();
        (void)I.decode_admit(0, 2, accepted);
        if (in.after) free_charge_guard(g, *in.after);
        if (Impl::idx(R[0]) != 0xffffffffu) throw std::runtime_error("actual particle endpoint charge/current mismatch");
        out.continuity_ratio = R[1].hi;
    };
    auto& sr = out.solve;
    sr.backend = "gpu1-resident";
    auto t = Clock::now();
    Vec rhs_v(g.ne), x_v(g.ne);
    if (I.opt.solve == ResidentOptions::Solve::device_cg) {
        I.cg_solve(in.dt, sr, check_A);
        I.s.download(I.rhs_host, 0, I.rhs, 0, g.ne * 8);
        I.s.download(I.x_host, 0, I.x, 0, g.ne * 8);
    } else {
        // parity mode: the sealed CPU PCG on the downloaded right-hand side
        I.s.download(I.rhs_host, 0, I.rhs, 0, g.ne * 8);
        I.s.synchronize("rhs download");
        check_A();
        std::memcpy(rhs_v.data(), I.rhs_host.pointer, g.ne * 8);
        const auto& S = I.S_host;
        Vec diagonal(g.ne, 1), sqrtM(g.ne);
        const double dt = in.dt;
        for (int e = 0; e < g.ne; ++e) { sqrtM[e] = std::sqrt(g.M[e]); if (!g.fixed_E[e]) diagonal[e] += .5 * dt * (S[e] + g.Y[e]) / g.M[e]; }
        for (int b = 0; b < g.nb; ++b) for (int q = 0; q < 4; ++q) { int e = g.C[b].index[q]; if (!g.fixed_E[e]) diagonal[e] += .25 * dt * dt * g.H[b] * std::pow(g.C[b].coefficient[q], 2) / g.M[e]; }
        auto apply = [&](const Vec& v) {
            Vec E(g.ne);
            for (int e = 0; e < g.ne; ++e) if (!g.fixed_E[e]) E[e] = v[e] / sqrtM[e];
            auto K = g.curl_adjoint(g.curl(E));
            Vec result(g.ne);
            for (int e = 0; e < g.ne; ++e) result[e] = g.fixed_E[e] ? v[e] : v[e] + (.5 * dt * (S[e] + g.Y[e]) * E[e] + .25 * dt * dt * K[e]) / sqrtM[e];
            return result;
        };
        const std::function<Vec(const Vec&)> f = apply;
        x_v = foil::solve(f, rhs_v, diagonal, I.maxit);
        sr.backend = "cpu-pcg-resident";
        std::memcpy(I.x_host.pointer, x_v.data(), g.ne * 8);
        I.s.upload(I.x, 0, I.x_host, 0, g.ne * 8);
    }
    const int marker = I.s.mark("rhs/x downloaded");
    sr.solve_s = std::chrono::duration<double>(Clock::now() - t).count();
    out.host_ms["solve"] = ms_since(t0); t0 = Clock::now();
    // post-solve arithmetic (queued; the host residual check runs while it executes)
    I.launch_ctx("k_update_E", g.ne);
    I.launch_ctx("k_update_B", g.nb);
    I.launch_ctx("k_currents", g.ne);
    {
        std::uint64_t mh = I.md_hi.address, ml = I.md_lo.address, oh = I.od_hi.address, ol = I.od_lo.address;
        I.launch_ctx("k_ledger", g.nn, {&mh, &ml, &oh, &ol});
        I.ctx.new_ = I.rf(I.phys(1));
        // v3.2: candidate field charge, then the Gauss increment and the admitted Gauss ledger of the
        // candidate, then the candidate admission (relative to that ledger), rows at 3 RF_ROWS
        int one = 1; std::uint64_t q1h = I.Q1_hi.address, q1l = I.Q1_lo.address;
        I.launch_ctx("k_charge", g.nn, {&one, &q1h, &q1l});
        I.launch_ctx("k_gauss_increment", g.nn, {&mh, &ml, &oh, &ol});
        I.enqueue_admit(1, 3 * RF_ROWS, I.Q1_hi.address, I.Q1_lo.address);
    }
    if (I.opt.trace_residual_charge) {
        int nn = g.nn; std::uint64_t a = I.ginc_acc[I.ginc_a].address, add = I.nvp(NV_GINC), o = I.ginc_acc[I.ginc_a ^ 1].address;
        I.launch("k_add_f64", {blocks(nn)}, {TPB}, {&nn, &a, &add, &o});
    }
    I.launch_ctx("k_cells", int(g.cells.size()));
    {
        std::uint64_t vn = I.VN.address;
        double factor = control_electron_charge_C / (hall::eps * control_electron_mass_kg);
        I.launch_ctx("k_omega_p", g.nn, {&vn, &factor});
        int nm = int(I.metal_cells.size());
        std::uint64_t mi = I.metal_idx.address, jc = I.joule_cell.address, jm = I.joule_metal.address;
        I.launch("k_gather_f64", {blocks(nm)}, {TPB}, {&nm, &mi, &jc, &jm});
    }
    {
        RfRows rr{};
        for (int t2 : {ET_FDE}) Impl::row(rr, I.etp(t2), g.ne, 0);
        Impl::row(rr, I.ftp(FT_FDB), g.nb, 0);
        for (int t2 : {ET_PW, ET_JOULE, ET_CIRCUIT, ET_OLOSS, ET_OSUPPLY, ET_MDE2, ET_MTHIRD2, ET_ME2, ET_TERMINAL, ET_REG_SOLID, ET_REG_PLASMA, ET_REG_DISP, ET_REG_EXT})
            Impl::row(rr, I.etp(t2), g.ne, 0);
        Impl::row(rr, I.etp(ET_AMPRES), g.ne, 1); Impl::row(rr, I.nvp(NV_KCL), g.nn, 1); Impl::row(rr, I.nvp(NV_GINC), g.nn, 1);
        Impl::row(rr, I.nvp(NV_GACT), g.nn, 1); Impl::row(rr, I.cvp(CV_BMAG), int(g.cells.size()), 1); Impl::row(rr, I.nvp(NV_OMEGAP), g.nn, 1);
        Impl::row(rr, I.code_joule.address, int(g.cells.size()), 2);
        need(rr.count == 22, "row layout");
        I.reduce(rr, RF_ROWS);
    }
    // identities of the candidate state and midpoint (device leaves + host sample)
    const int pslot = I.acc ^ 1;
    I.slot_residue[pslot] = residue_present;
    const RfLeafTable st = I.state_table(1, I.leaves_dev.address);
    const int nstate = Impl::leaf_count(st);
    const RfLeafTable mt = I.midpoint_table(I.leaves_dev.address + std::uint64_t(nstate) * 32, I.midE.address, I.midB.address);
    const int nmid = Impl::leaf_count(mt);
    I.enqueue_leaves(st); I.enqueue_leaves(mt);
    const int ns = I.opt.verify_all ? 0 : I.enqueue_samples(st, I.opt.verify_samples, 0);
    const int nms = I.opt.verify_all ? 0 : I.enqueue_samples(mt, std::max(1, I.opt.verify_samples / 4), ns);
    I.s.download(I.res_host, std::size_t(RF_ROWS) * 16, I.res, std::size_t(RF_ROWS) * 16, std::size_t(RF_ROWS) * 16);
    I.s.download(I.res_host, std::size_t(3 * RF_ROWS) * 16, I.res, std::size_t(3 * RF_ROWS) * 16, std::size_t(admit_rows) * 16);
    I.s.download(I.joule_host, 0, I.joule_metal, 0, I.metal_cells.size() * 8);
    I.s.download(I.leaves_host, 0, I.leaves_dev, 0, std::size_t(nstate + nmid) * 32);
    // independent residual on the host (overlaps the queued device work)
    I.s.wait_marker(marker, "rhs/x");
    if (I.opt.solve == ResidentOptions::Solve::device_cg) {
        std::memcpy(rhs_v.data(), I.rhs_host.pointer, g.ne * 8);
        std::memcpy(x_v.data(), I.x_host.pointer, g.ne * 8);
    }
    need(axial::finite(x_v), "field solve backend returned a misshaped or nonfinite increment");
    t = Clock::now();
    I.host_apply.apply(I.S_host, in.dt, x_v, I.Ax_host);
    I.r_host.resize(g.ne);
    for (int e = 0; e < g.ne; ++e) I.r_host[e] = rhs_v[e] - I.Ax_host[e];
    sr.rhs_norm = foil::norm(rhs_v); sr.target = std::max(1e-30, 2e-13 * sr.rhs_norm); sr.independent_residual = foil::norm(I.r_host);
    if (!(sr.independent_residual <= sr.target)) {
        std::ostringstream why;
        why << std::setprecision(17) << "independent FP64 field residual check failed: backend=" << sr.backend << " residual=" << sr.independent_residual
            << " target=" << sr.target;
        throw std::runtime_error(why.str());
    }
    sr.check_s = std::chrono::duration<double>(Clock::now() - t).count();
    if (I.opt.trace_residual_charge) {
        double smax = 0, cmax = 0;
        for (int n = 0; n < g.nn; ++n) {
            long double s = 0;
            for (int j = I.h_ntp[n]; j < I.h_ntp[n + 1]; ++j) { const int e = I.h_nte[j]; s += (long double)I.host_apply.sqrtM[e] * std::abs(I.r_host[e]) / g.edges[e].length; }
            I.rq_cand[n] = I.rq_acc[n] + double(s); smax = std::max(smax, double(s)); cmax = std::max(cmax, I.rq_cand[n]);
        }
        out.residual_charge_step_max_C = smax; out.residual_charge_cumulative_max_C = cmax;
    }
    out.host_ms["independent_residual"] = ms_since(t0); t0 = Clock::now();
    I.s.synchronize("resident step");
    I.s.clear_markers();
    out.host_ms["device_wait"] = ms_since(t0); t0 = Clock::now();
    // decode in the sealed order
    const Wd* R = I.results() + RF_ROWS;
    auto& r = out.receipt;
    r.current_A = in.current;
    r.field_delta_J = (R[0].hi + R[0].lo) + (R[1].hi + R[1].lo);
    r.plasma_work_J = R[2].hi + R[2].lo; r.joule_J = R[3].hi + R[3].lo; r.circuit_work_J = -(R[4].hi + R[4].lo);
    r.outer_loss_J = R[5].hi + R[5].lo; r.outer_supply_J = R[6].hi + R[6].lo;
    const double mde2 = R[7].hi + R[7].lo, mthird2 = R[8].hi + R[8].lo, me2 = R[9].hi + R[9].lo;
    r.terminal_V = -(R[10].hi + R[10].lo);
    r.region_conduction_A = R[11].hi + R[11].lo; r.region_plasma_A = R[12].hi + R[12].lo; r.region_displacement_A = R[13].hi + R[13].lo;
    r.region_boundary_A = R[14].hi + R[14].lo;
    r.ampere_residual_A_m = R[15].hi; r.kcl_residual_A = R[16].hi; r.gauss_increment_C = R[17].hi;
    const double activityQ = R[18].hi;
    const double Bmax = R[19].hi, omega_p_dt = R[20].hi;
    const bool joule_rule = Impl::idx(R[21]) == 0xffffffffu;
    r.balance_J = r.field_delta_J + r.plasma_work_J + r.joule_J + r.outer_loss_J - r.outer_supply_J - r.circuit_work_J;
    need(r.gauss_increment_C <= 1e-22 + 2e-10 * activityQ, "existing field Gauss increment bound");
    const double activity = std::abs(r.field_delta_J) + std::abs(r.plasma_work_J) + std::abs(r.joule_J) + std::abs(r.outer_loss_J) + std::abs(r.outer_supply_J) +
                            std::abs(r.circuit_work_J);
    need(std::isfinite(r.balance_J) && std::abs(r.balance_J) <= 1e-20 + 1e-10 * activity, "unchanged field/current/Joule energy bound");
    SlotMeta cand = accepted;
    cand.time_s = accepted.time_s + in.dt; cand.revision = accepted.revision + 1; cand.residue_present = residue_present;
    // candidate admission (rows at 3 RF_ROWS)
    out.constraints = I.decode_admit(1, 3 * RF_ROWS, cand);
    need(joule_rule, "v2 heat Joule power outside metal or invalid");
    for (std::size_t i = 0; i < I.metal_cells.size(); ++i) I.cell_joule[I.metal_cells[i]] = I.joule_host.as<double>()[i];
    r.cell_joule_W = I.cell_joule;
    // leaves
    auto& sl = I.slot_leaves[pslot];
    sl.resize(nstate);
    std::memcpy(sl.data(), I.leaves_host.pointer, std::size_t(nstate) * 32);
    I.mid_leaves.resize(nmid);
    std::memcpy(I.mid_leaves.data(), I.leaves_host.as<unsigned char>() + std::size_t(nstate) * 32, std::size_t(nmid) * 32);
    bool ok_leaves;
    if (I.opt.verify_all) ok_leaves = I.verify_all_leaves(st, sl.data()) && I.verify_all_leaves(mt, I.mid_leaves.data());
    else ok_leaves = I.verify_samples(st, sl.data(), ns, 0) && I.verify_samples(mt, I.mid_leaves.data(), nms, ns);
    need(ok_leaves, "device block digest differs from the host recomputation");
    I.slot_leaves_valid[pslot] = true; I.mid_leaves_valid = true;
    // controller receipt
    auto& ctl = out.control;
    ctl.omega_p_dt_max = omega_p_dt; ctl.B_max_T = Bmax; ctl.omega_c_dt_max = control_electron_charge_C * Bmax / control_electron_mass_kg * in.dt;
    ctl.drive_current_A = in.current; ctl.drive_rel_change = in.drive_rel_change; ctl.drive_knot_interval_s = in.drive_knot_interval_s;
    ctl.history = I.history;
    ctl.increment_rel = me2 > 0 ? std::sqrt(mde2 / me2) : 0.;
    ctl.third_difference_rel = (I.history >= 2 && me2 > 0) ? std::sqrt(mthird2 / me2) / 12. : 0.;
    out.residue_present = residue_present;
    out.host_ms["decode"] = ms_since(t0);
    return out;
}
std::string ResidentField::state_digest(int role, const std::string& api, const std::string& configuration, const SlotMeta& m) {
    auto& I = *p_;
    const auto t = Clock::now();
    const int ps = I.acc ^ role;
    need(I.slot_residue[ps] == m.residue_present, "slot residue presence differs from its metadata");
    if (!I.slot_leaves_valid[ps]) {
        const RfLeafTable st = I.state_table(role, I.leaves_dev.address);
        const int n = Impl::leaf_count(st);
        I.enqueue_leaves(st);
        const int ns = I.opt.verify_all ? 0 : I.enqueue_samples(st, I.opt.verify_samples, 0);
        I.s.download(I.leaves_host, 0, I.leaves_dev, 0, std::size_t(n) * 32);
        I.s.synchronize("state leaves");
        I.slot_leaves[ps].resize(n);
        std::memcpy(I.slot_leaves[ps].data(), I.leaves_host.pointer, std::size_t(n) * 32);
        const bool ok = I.opt.verify_all ? I.verify_all_leaves(st, I.slot_leaves[ps].data()) : I.verify_samples(st, I.slot_leaves[ps].data(), ns, 0);
        need(ok, "device block digest differs from the host recomputation");
        I.slot_leaves_valid[ps] = true;
    }
    auto d = I.state_root(role, api, configuration, m, I.slot_leaves[ps]);
    ++I.dstats.digests; I.dstats.host_ms += ms_since(t);
    return d;
}
std::string ResidentField::midpoint_seal() {
    auto& I = *p_;
    need(I.mid_leaves_valid, "midpoint leaves unavailable");
    return I.midpoint_root(I.mid_leaves);
}
void ResidentField::midpoint_download(Vec& E, Vec& B) {
    auto& I = *p_;
    E.resize(I.g.ne); B.resize(I.g.nb);
    I.download_now(I.midE, E.data(), I.g.ne * 8);
    I.download_now(I.midB, B.data(), I.g.nb * 8);
}
void ResidentField::accepted_rounded(Vec& E, Vec& B) {
    auto& I = *p_;
    I.refresh_ctx();
    std::uint64_t ce = I.chkE.address, cb = I.chkB.address;
    I.launch_ctx("k_rounded", std::max(I.g.ne, I.g.nb), {&ce, &cb});
    E.resize(I.g.ne); B.resize(I.g.nb);
    I.download_now(I.chkE, E.data(), I.g.ne * 8);
    I.download_now(I.chkB, B.data(), I.g.nb * 8);
}
void ResidentField::begin_fixed_point() {
    auto& I = *p_;
    I.refresh_ctx();
    std::uint64_t pe = I.prevE.address, pb = I.prevB.address;
    I.launch_ctx("k_rounded", std::max(I.g.ne, I.g.nb), {&pe, &pb});
}
bool ResidentField::midpoint_unchanged() {
    auto& I = *p_;
    const auto& g = I.g;
    int ne = g.ne, nb = g.nb;
    std::uint64_t a1 = I.midE.address, b1 = I.prevE.address, a2 = I.midB.address, b2 = I.prevB.address, f = I.c_diff.address, f2 = I.c_diff2.address;
    I.launch("k_differ", {blocks(ne)}, {TPB}, {&ne, &a1, &b1, &f});
    I.launch("k_differ", {blocks(nb)}, {TPB}, {&nb, &a2, &b2, &f2});
    RfRows rr{};
    Impl::row(rr, f, ne, 2); Impl::row(rr, f2, nb, 2);
    I.reduce(rr, 2 * RF_ROWS);
    I.s.download(I.res_host, std::size_t(2 * RF_ROWS) * 16, I.res, std::size_t(2 * RF_ROWS) * 16, 32);
    I.s.synchronize("midpoint comparison");
    const Wd* R = I.results() + 2 * RF_ROWS;
    const bool same = Impl::idx(R[0]) == 0xffffffffu && Impl::idx(R[1]) == 0xffffffffu;
    if (!same) {
        I.launch("k_copy_f64", {blocks(ne)}, {TPB}, {&ne, &a1, &b1});
        I.launch("k_copy_f64", {blocks(nb)}, {TPB}, {&nb, &a2, &b2});
    }
    return same;
}
bool ResidentField::preflight_recheck(const SlotMeta& accepted, const SlotMeta& candidate) {
    auto& I = *p_;
    const auto& g = I.g;
    try {
        if (!I.slot_leaves_valid[I.acc] || !I.slot_leaves_valid[I.acc ^ 1] || !I.mid_leaves_valid) return false;
        I.refresh_ctx();
        // midpoint recomputation from both slots, compared with the stored midpoint
        std::uint64_t ce = I.chkE.address, cb = I.chkB.address;
        I.launch_ctx("k_midpoint", std::max(g.ne, g.nb), {&ce, &cb});
        int ne = g.ne, nb = g.nb;
        std::uint64_t mE = I.midE.address, mB = I.midB.address, f = I.c_diff.address, f2 = I.c_diff2.address;
        I.launch("k_differ", {blocks(ne)}, {TPB}, {&ne, &ce, &mE, &f});
        I.launch("k_differ", {blocks(nb)}, {TPB}, {&nb, &cb, &mB, &f2});
        // leaves of both slots, recomputed
        const RfLeafTable ta = I.state_table(0, I.leaves_chk.address);
        const int na = Impl::leaf_count(ta);
        const RfLeafTable tc = I.state_table(1, I.leaves_chk.address + std::uint64_t(na) * 32);
        const int nc = Impl::leaf_count(tc);
        I.enqueue_leaves(ta); I.enqueue_leaves(tc);
        // admissions of both slots
        I.enqueue_admit(0, 0, I.Q0_hi.address, I.Q0_lo.address);
        {
            RfRows rr{};
            Impl::row(rr, f, ne, 2); Impl::row(rr, f2, nb, 2);
            I.reduce(rr, admit_rows);
        }
        I.enqueue_admit(1, RF_ROWS, I.Q1_hi.address, I.Q1_lo.address);
        const int s1 = I.enqueue_samples(ta, std::max(1, I.opt.verify_samples / 2), 0);
        const int s2 = I.enqueue_samples(tc, std::max(1, I.opt.verify_samples / 2), s1);
        I.s.download(I.res_host, 0, I.res, 0, std::size_t(2) * RF_ROWS * 16);
        I.s.download(I.chk_host, 0, I.leaves_chk, 0, std::size_t(na + nc) * 32);
        I.s.synchronize("preflight recheck");
        const Wd* R = I.results();
        (void)I.decode_admit(0, 0, accepted);
        (void)I.decode_admit(1, RF_ROWS, candidate);
        if (Impl::idx(R[admit_rows]) != 0xffffffffu || Impl::idx(R[admit_rows + 1]) != 0xffffffffu) return false;
        if (std::memcmp(I.chk_host.pointer, I.slot_leaves[I.acc].data(), std::size_t(na) * 32) != 0) return false;
        if (std::memcmp(I.chk_host.as<unsigned char>() + std::size_t(na) * 32, I.slot_leaves[I.acc ^ 1].data(), std::size_t(nc) * 32) != 0) return false;
        if (!I.verify_samples(ta, I.slot_leaves[I.acc].data(), s1, 0) || !I.verify_samples(tc, I.slot_leaves[I.acc ^ 1].data(), s2, s1)) return false;
        return true;
    } catch (const std::exception&) {
        return false;
    }
}
void ResidentField::publish() noexcept {
    auto& I = *p_;
    I.acc ^= 1;
    const int old2 = I.inc2;
    I.inc2 = I.inc1; I.inc1 = I.inc0; I.inc0 = old2;
    if (I.history < 2) ++I.history;
    I.mid_leaves_valid = false;
    I.slot_leaves_valid[I.acc ^ 1] = false;
    if (I.opt.trace_residual_charge) { I.rq_acc.swap(I.rq_cand); I.ginc_a ^= 1; }
}
}  // namespace hall_corrective_v3::gpu
