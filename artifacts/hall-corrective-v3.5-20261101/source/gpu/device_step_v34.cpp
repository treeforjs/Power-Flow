// Hall v3.4 device coupled step: host engine (see device_step_v34.hpp).
#include "device_step_v34.hpp"
#include "gpu_backends.hpp"
#include "common/hall_gpu_common.hpp"
#include "common/hall_gpu_primitives.hpp"
#include "push/hall_adapter.hpp"
#include "hall_gpu_push_v3.hpp"
#include "gpu/forms_gpu.hpp"
#include "support_host_v34.hpp"
#include "step_abi_v34.h"
#include "../sink_core.hpp"
#include <algorithm>
#include <chrono>
#include <ctime>
#include <cstring>
#include <set>

namespace hall_corrective_v3::gpu {
namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kKernels[] = {
#embed "step_kernels_v34.cu"
    , 0};
const char kAbi[] = {
#embed "step_abi_v34.h"
    , 0};
const char kExactSum[] = {
#embed "exact_sum_v34.cuh"
    , 0};
const char kX87[] = {
#embed "x87_v34.hpp"
    , 0};
const char kSupport[] = {
#embed "support_v34.hpp"
    , 0};
const char kSinkCore[] = {
#embed "../sink_core.hpp"
    , 0};
const char kPushCore[] = {
#embed "../../../hall-gpu-resident-v1-20261030/push/push_core.hpp"
    , 0};
const char kFormsCore[] = {
#embed "../../../pic-reaction-transactions-20261029/native-full-time-candidate007.1/forms/include/pf_v7/forms_core.h"
    , 0};
#pragma GCC diagnostic pop

using Clock = std::chrono::steady_clock;
double ms_since(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }
double cpu_ms() { timespec ts{}; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); return ts.tv_sec * 1e3 + ts.tv_nsec * 1e-6; }
struct CpuTimer { double& acc; double t0 = cpu_ms(); ~CpuTimer() { acc += cpu_ms() - t0; } };
constexpr unsigned TPB = S34_TPB, SUM_BLOCKS = 120;
constexpr std::size_t XS = 1104, SPEC = 1024;
unsigned grid(std::size_t n) { return unsigned(std::max<std::size_t>((n + TPB - 1) / TPB, 1)); }
unsigned key_bits_for(std::size_t max_key) { unsigned b = 1; while ((std::uint64_t(1) << b) <= max_key) ++b; return b; }
std::size_t a256(std::size_t b) { return (b + 255) & ~std::size_t(255); }
// result region offsets (bytes)
constexpr std::size_t R_ITER = 0, R_CENSUS = 256, R_CHECK = 512, R_CROOT = 544, R_SROOT = 576, R_FROOT = 7680, R_PUSH_STATUS = 768, R_PUSH_MAX = 800, R_PRE = 896,
                      R_SUMS = 1024, R_SINKS = 8192, R_SINK_INDEX = R_SINKS + SPEC * 136, R_SINK_EXTRA = R_SINK_INDEX + SPEC * 8,
                      R_BYTES = R_SINK_EXTRA + SPEC * 16;
struct XsHost { double value, scaled; std::uint32_t status, parts; std::uint64_t terms; std::int64_t digit[134]; };
static_assert(sizeof(XsHost) == XS);
const char* support_text(int code) {
    switch (code) {
        case sup34::charge_operand: return "deposit support charge operand";
        case sup34::nonfinite_coordinate: return "deposit support nonfinite coordinate";
        default: return "deposit support path exceeds the motion displacement guard";
    }
}
std::string hex_words(const std::uint32_t* w) { return pf_v7::hex(pf_v7::digest_of_words(w)); }
}  // namespace

struct DeviceStep::Impl {
    std::shared_ptr<Device> device;
    hall_gpu::Session& s;
    const Model& model;
    const Mesh& g;
    MotionContext ctx;
    std::size_t N;
    std::shared_ptr<pf_v7::gpu::Forms> forms;   // shared per device (Device::cached)
    hall_gpu::Module mod;
    std::map<std::string, hall_gpu::Kernel> k;
    std::unique_ptr<hall_gpu::push_v3::MotionEngineV3> push;
    sup34::Geo geo{};
    hgp::Mesh hm{};
    hgp::Policy pol{};
    unsigned node_bits = 0;
    ResidentField::Tables tab{};
    // census
    hall_gpu::HostBuffer h_rows, h_keys, h_states, h_nb, h_res, h_rec;
    hall_gpu::DeviceBuffer rows, keys, states, ids, crec, cscratch, chk_scratch, k0, tau_c, axis_c;
    hall_gpu::DeviceBuffer nb_ids, nb_tau, nb_idx, nb_code, nb_rank, nb_face;
    hall_gpu::DeviceBuffer ckey_a, cpay_a, ckey_b, cpay_b, chi, clo, ccode, cfail, bounds, sort_scratch, c_hi, c_lo, c_gross;
    hall_gpu::DeviceBuffer fail_nn, fail_nn2, fail_n, fail_ne, red_scratch, scan_scratch;
    // fold
    hall_gpu::DeviceBuffer sum_part, sum_stat, cont, marks, ebound, efail_b, e_hi, e_lo, e_gross, grows;
    // support
    hall_gpu::DeviceBuffer t1, t2, box, scount, soff, scode, sfail, skey_a, spay_a, skey_b, spay_b, support_out, ssort_scratch;
    std::size_t support_cap = 0;
    // record
    hall_gpu::DeviceBuffer surv, sinkf, bornf, off_surv, off_sink, off_born, rec_end, rec_path, rec_sink, rec_birth, sink_index, birth_index, sink_extra,
        bl_scratch[4], froot_scratch, flists_scratch;
    hall_gpu::DeviceBuffer d_res;
    // per proposal
    std::size_t n = 0;
    const native::ArenaRateEntry* rows_host = nullptr; std::size_t rows_host_size = 0;
    std::vector<const native::ArenaRateEntry*> order;   // sorted rows when the input was not ascending
    CensusSource src;
    double begin_s = 0, end_s = 0, dt = 0;
    std::string pending_identity, pending_census, pending_newborn;   // deferred refusals (motion admission order)
    bool any_resident_neutral = false;
    CensusV2 v2;
    std::vector<NewbornRecord> births;
    std::vector<std::uint64_t> nb_index_host;   // census index per newborn (ascending ID order of births)
    std::vector<int> newborn_flag_host;
    // per iterate
    hall_gpu::push_v3::MotionEngineV3::DeviceRun run;
    int support_code = 0; long support_fail = -1;
    std::size_t sinks = 0, survivors = 0, nbirths = 0;
    // seal
    std::string sealed_record, sealed_current_root;
    std::array<std::string, 3> sealed_array_digests;
    FieldDigestForm sealed_form = FieldDigestForm::streaming_v7;
    std::string sealed_trial_blocks;
    std::uint64_t seal_generation = 0;
    DeviceStep::DeviceRecordView view;
    DeviceStepTimes times;

    Impl(std::shared_ptr<Device> d, const Model& m, const MotionContext& c, std::size_t capacity)
        : device(std::move(d)), s(device->session()), model(m), g(m.mesh()), ctx(c), N(std::max<std::size_t>(capacity, 1)) {
        need(c.mesh == &g, "device step motion context belongs to another mesh");
        auto& dv = *device;
        forms = std::static_pointer_cast<pf_v7::gpu::Forms>(dv.cached("pf_v7-forms", [&] { return std::static_pointer_cast<void>(std::make_shared<pf_v7::gpu::Forms>(s)); }));
        mod = *std::static_pointer_cast<const hall_gpu::Module>(dv.cached("hall-v3.4-step-module", [&] {
            hall_gpu::CompileRequest rq;
            rq.main = {"step_kernels_v34.cu", kKernels};
            rq.headers = {{"step_abi_v34.h", kAbi}, {"exact_sum_v34.cuh", kExactSum}, {"x87_v34.hpp", kX87}, {"support_v34.hpp", kSupport},
                          {"sink_core.hpp", kSinkCore}, {"push_core.hpp", kPushCore}, {"forms_core.h", kFormsCore}};
            return std::static_pointer_cast<void>(std::make_shared<hall_gpu::Module>(s.compile_and_load(rq))); }));
        for (const char* nm : {"k34_rows_to_soa", "k34_find_ids", "k34_flights_init", "k34_newborn_place", "k34_k0", "k34_neutral_fail",
                               "k34_records_charge", "k34_differ_u64", "k34_charge_nodes", "k34_census_nodes", "k34_fail_nonfinite", "k34_sum",
                               "k34_sum_finish", "k34_continuity", "k34_node_marks", "k34_edge_bound", "k34_free_guard", "k34_assemble_check",
                               "k34_record_flags", "k34_record_build", "k34_birth_build", "k34_support_carriers", "k34_support_emit",
                               "k34_support_nodes", "k34_f64_chunks", "k34_f64_blocks", "k34_f64_tree", "k34_fill_u32", "k34_gather_rows"})
            k[nm] = s.kernel(mod, nm);
        // push engine (v3 kernels; device path)
        namespace hp = hall_gpu::push; namespace ha = hall_gpu::push::hall_adapter; namespace h3 = hall_gpu::push_v3;
        hp::EngineOptions eo; eo.capacity = N;
        h3::SinkSupport sup;
        for (int a = 0; a < 5; ++a) { for (int dd = 0; dd < 3; ++dd) { sup.lo[a][dd] = c.sink_lo[a][dd]; sup.hi[a][dd] = c.sink_hi[a][dd]; } sup.supported[a] = c.sink_supported[a]; }
        sup.aligned = c.sink_support_aligned;
        push = std::make_unique<h3::MotionEngineV3>(s, ha::mesh_tables(g), ha::geometry(c.options), sup, eo, h3::MotionEngineV3::shared_module(dv).get());
        geo = sup34::geo_of(g, c.options);
        hm = ha::mesh_tables(g).mesh;
        pol = ha::geometry(c.options).policy;
        node_bits = key_bits_for(std::size_t(g.nn));
        build();
    }
    // every allocation of the engine is released with it
    std::vector<hall_gpu::DeviceBuffer> owned_dev;
    std::vector<hall_gpu::HostBuffer> owned_host;
    hall_gpu::DeviceBuffer dev(std::size_t bytes, const char* tag) { auto b = s.alloc(std::max<std::size_t>(bytes, 16), tag); owned_dev.push_back(b); return b; }
    hall_gpu::HostBuffer pinned(std::size_t bytes, const char* tag) { auto h = s.alloc_pinned(bytes, tag); owned_host.push_back(h); return h; }
    void release_owned(hall_gpu::DeviceBuffer& b) {
        for (auto it = owned_dev.begin(); it != owned_dev.end(); ++it) if (it->address == b.address && it->slot == b.slot) { owned_dev.erase(it); break; }
        s.release(b);
    }
    ~Impl() {
        push.reset();   // its own buffers first
        if (s.failed()) return;
        try { for (auto& b : owned_dev) s.release(b); for (auto& h : owned_host) s.release(h); } catch (...) {}
    }
    void build() {
        const std::size_t nn = std::size_t(g.nn), ne = std::size_t(g.ne), nb = std::size_t(g.nb);
        h_rows = pinned(N * 80, "census rows"); h_keys = pinned(N * 8, "census keys"); h_states = pinned(N * 8, "census states");
        h_nb = pinned(N * 40 + 64, "newborns"); h_res = pinned(R_BYTES, "step results");
        h_rec = pinned((ne + 2 * nn) * 8 + N * 144 + 64, "record arrays");
        rows = dev(N * 80, "census rows"); keys = dev(N * 8, "census keys"); states = dev(N * 8, "census states"); ids = dev(N * 8, "ids");
        crec = dev(N * 96, "census records"); cscratch = dev(pf_v7::gpu::Forms::block_list_scratch_bytes(N), "census block list");
        chk_scratch = dev(pf_v7::gpu::Forms::census_check_scratch_bytes(N), "census check");
        k0 = dev(N * 8, "K0"); tau_c = dev(N * 8, "birth time"); axis_c = dev(N * 4, "face axis");
        for (auto* b : {&nb_ids, &nb_tau, &nb_idx}) *b = dev(N * 8, "newborns");
        nb_code = dev(N * 4, "newborn codes"); nb_rank = dev(N * 4, "newborn ranks"); nb_face = dev(N * 8, "newborn faces");
        const std::size_t R = 8 * N;
        for (auto* b : {&ckey_a, &cpay_a, &ckey_b, &cpay_b}) *b = dev(R * 4, "charge records");
        chi = dev(R * 8, "charge hi"); clo = dev(R * 8, "charge lo"); ccode = dev(N * 4, "charge codes"); cfail = dev(N * 4, "charge fails");
        bounds = dev(2 * (std::max(nn, ne) + 1) * 4, "segment bounds");
        sort_scratch = dev(hall_gpu::Primitives::sort_scratch(R) * 4, "sort scratch");
        c_hi = dev(nn * 8, "census charge hi"); c_lo = dev(nn * 8, "census charge lo"); c_gross = dev(nn * 8, "census charge gross");
        fail_nn = dev(nn * 4, "node fails"); fail_nn2 = dev(nn * 4, "node fails"); fail_n = dev(N * 4, "carrier fails");
        fail_ne = dev(std::max(ne, nb) * 4, "edge fails"); efail_b = dev(nb * 4, "face fails");
        red_scratch = dev(hall_gpu::Primitives::reduce_scratch(std::max({N, ne, nn, nb})) * 8 + 64, "reduce scratch");
        scan_scratch = dev(hall_gpu::Primitives::scan_scratch(N) * 4 + 64, "scan scratch");
        sum_part = dev(std::size_t(SUM_BLOCKS) * 8 * 134, "sum partials"); sum_stat = dev(std::size_t(SUM_BLOCKS) * 4, "sum status");
        cont = dev(nn * 8, "continuity"); marks = dev(nn * 4, "node marks"); ebound = dev(ne * 8, "edge bound");
        e_hi = dev(nn * 8, "endpoint charge hi"); e_lo = dev(nn * 8, "endpoint charge lo"); e_gross = dev(nn * 8, "endpoint charge gross");
        grows = dev(N * 80, "gathered rows");
        t1 = dev(N * 16, "support t1"); t2 = dev(N * 16, "support t2"); box = dev(N * 24, "support box");
        for (auto* b : {&scount, &soff, &scode, &sfail}) *b = dev(N * 4, "support carriers");
        grow_support(12 * N);
        support_out = dev(nn * 8, "support"); s.zero(support_out, 0, nn * 8);
        for (auto* b : {&surv, &sinkf, &bornf, &off_surv, &off_sink, &off_born}) *b = dev(N * 4, "record flags");
        rec_end = dev(N * 80, "endpoint records"); rec_path = dev(N * 64, "path records"); rec_sink = dev(N * 136, "sink records");
        rec_birth = dev(N * 24, "birth records"); sink_index = dev(N * 8, "sink index"); birth_index = dev(N * 8, "birth index");
        sink_extra = dev(N * 16, "sink extra");
        for (auto& b : bl_scratch) b = dev(pf_v7::gpu::Forms::block_list_scratch_bytes(N), "section block list");
        froot_scratch = dev(pf_v7::gpu::Forms::block_list_scratch_bytes(std::max(ne, nb)) + 4096, "field block list");
#if HALL_FORMS_V71
        {
            const pf_v7::gpu::Forms::BlockListJob jobs[3] = {{pf7::section_current_J, ne, 0, 0}, {pf7::section_field_E, ne, 0, 0},
                                                             {pf7::section_field_B, nb, 0, 0}};
            flists_scratch = dev(pf_v7::gpu::Forms::block_lists_scratch_bytes(jobs) + 4096, "field block lists");
        }
#endif
        d_res = dev(R_BYTES, "step results");
        s.synchronize("device step buffers");
    }
    void grow_support(std::size_t records) {
        if (records <= support_cap) return;
        for (auto* b : {&skey_a, &spay_a, &skey_b, &spay_b, &ssort_scratch}) if (*b) release_owned(*b);
        support_cap = records;
        for (auto* b : {&skey_a, &spay_a, &skey_b, &spay_b}) *b = dev(records * 4, "support records");
        ssort_scratch = dev(hall_gpu::Primitives::sort_scratch(records) * 4, "support sort scratch");
    }
    hall_gpu::Primitives& P() { return forms->primitives(); }
    void launch(const char* name, std::size_t n_threads, std::vector<void*> args) {
        s.launch(k.at(name), {grid(n_threads)}, {TPB}, 0, args.data());
    }
    std::uint64_t res(std::size_t off) const { return d_res.address + off; }
    void min_u32(std::uint64_t values, std::size_t n_, std::size_t off) {
        if (n_ == 0) { std::uint64_t one = 1, at = res(off); std::uint32_t v = S34_NONE; launch("k34_fill_u32", 1, {&one, &at, &v}); return; }
        P().reduce_min_u32(values, n_, red_scratch.address, res(off));
    }
    void max_f64(std::uint64_t values, std::size_t n_, std::size_t off) {
        if (n_ == 0) { s.zero(d_res, off, 8); return; }
        P().reduce_max_f64(values, n_, red_scratch.address, res(off));
    }
    // Sorts 8-node records (key, pay) by node and walks each node's records in record order.
    void charge_records(std::size_t count, std::uint64_t key_a, std::uint64_t pay_a, std::uint64_t key_b, std::uint64_t pay_b, std::uint64_t hi,
                        std::uint64_t lo, std::uint64_t out_hi, std::uint64_t out_lo, std::uint64_t out_gross) {
        std::uint32_t nn = std::uint32_t(g.nn);
        std::uint64_t keys_s = key_a, pays_s = pay_a;
        s.zero(bounds, 0, 2 * (std::size_t(g.nn) + 1) * 4);
        if (count) {
            const bool in_b = P().radix_sort(key_a, pay_a, key_b, pay_b, count, node_bits, sort_scratch.address);
            keys_s = in_b ? key_b : key_a; pays_s = in_b ? pay_b : pay_a;
            P().segment_bounds(keys_s, count, bounds.address, bounds.address + (std::size_t(g.nn) + 1) * 4);
        }
        std::uint64_t b0 = bounds.address, b1 = bounds.address + (std::size_t(g.nn) + 1) * 4;
        launch("k34_charge_nodes", g.nn, {&nn, &pays_s, &b0, &b1, &hi, &lo, &out_hi, &out_lo, &out_gross});
    }
    // One row of a device census (k34_gather_rows; O(1) per call, used for sinks and newborns only).
    native::GapParticle device_row(std::uint64_t i) {
        std::uint64_t one = 1, rows_d = src.device_rows, ia = grows.address + 80, out = grows.address;
        auto* h = h_rows.as<std::uint64_t>();
        h[0] = i;
        s.upload(grows, 80, h_rows, 0, 8);
        launch("k34_gather_rows", 1, {&one, &rows_d, &ia, &out});
        s.download(h_rows, 0, grows, 0, 80);
        s.synchronize("device census row");
        native::GapParticle p{};
        std::memcpy(&p, h_rows.pointer, 80);
        return p;
    }
    const native::ArenaRateEntry& row(std::size_t i) const { return order.empty() ? rows_host[i] : *order[i]; }
    // ---------------------------------------------------------------- sums
    void sum(int kind, const S34SumArgs& a, double scale, int slot) {
        std::uint64_t pa = sum_part.address, sa = sum_stat.address;
        S34SumArgs args = a;
        s.launch(k.at("k34_sum"), {SUM_BLOCKS}, {TPB}, 0, std::vector<void*>{&kind, &args, &pa, &sa}.data());
        unsigned parts = SUM_BLOCKS; std::uint64_t terms = kind == 2 ? 2 * a.n : a.n, out = res(R_SUMS + std::size_t(slot) * XS);
        s.launch(k.at("k34_sum_finish"), {1}, {TPB}, 0, std::vector<void*>{&parts, &pa, &sa, &terms, &scale, &out}.data());
    }
    // ---------------------------------------------------------------- block list of a binary64 array (v7.1 tags)
    void f64_root(std::uint32_t which, std::uint64_t values, std::uint64_t n_, std::size_t off) {
        std::uint64_t chunks = (n_ + pf7::chunk_records - 1) / pf7::chunk_records, leaves = (n_ + pf7::leaf_records - 1) / pf7::leaf_records;
        std::uint64_t chunk_buf = froot_scratch.address, leaf_buf = chunk_buf + a256(std::max<std::uint64_t>(chunks, 1) * 32);
        std::uint64_t a = leaf_buf + a256(std::max<std::uint64_t>(leaves, 1) * 32), b = a + a256((leaves + 1) / 2 * 32 + 32);
        std::uint64_t root = res(off);
        if (chunks) {
            launch("k34_f64_chunks", chunks, {&which, &n_, &values, &chunk_buf});
            launch("k34_f64_blocks", leaves, {&which, &n_, &chunk_buf, &leaf_buf});
        }
        s.launch(k.at("k34_f64_tree"), {1}, {32}, 0, std::vector<void*>{&which, &n_, &leaves, &leaf_buf, &a, &b, &root}.data());
    }
    // ---------------------------------------------------------------- field/current block roots (native v7.1 1.1)
    // Roots of the current J [ne] (R_FROOT), and of the trial E [ne] / B [nb] (R_FROOT + 32 / + 64) when given: the
    // fused Forms::block_lists of native v7.1 (one launch per level); without the v7.1 forms library the step's own
    // kernels (k34_f64_*, the same block lists and tags; qualification: test_v34 blocks).
    void field_roots(std::uint64_t J, std::uint64_t E, std::uint64_t B) {
        const std::uint64_t ne_ = std::uint64_t(g.ne), nb_ = std::uint64_t(g.nb);
#if HALL_FORMS_V71
        std::vector<pf_v7::gpu::Forms::BlockListJob> jobs{{pf7::section_current_J, ne_, J, res(R_FROOT)}};
        if (E && B) { jobs.push_back({pf7::section_field_E, ne_, E, res(R_FROOT + 32)}); jobs.push_back({pf7::section_field_B, nb_, B, res(R_FROOT + 64)}); }
        forms->block_lists(jobs, flists_scratch.address);
#else
        f64_root(2, J, ne_, R_FROOT);
        if (E && B) { f64_root(0, E, ne_, R_FROOT + 32); f64_root(1, B, nb_, R_FROOT + 64); }
#endif
    }
    template <class T> const T& at(std::size_t off) const { return *reinterpret_cast<const T*>(h_res.as<char>() + off); }
};
static std::string current_blocks_digest(const std::string& J_root) {
#if HALL_FORMS_V71
    return pf_v7::current_blocks_sha256(J_root);
#else
    Digest dj("motion-array-blocks-v1"); dj.text(J_root); return dj.hex();
#endif
}
static std::string trial_blocks_digest(const std::string& E_root, const std::string& B_root) {
#if HALL_FORMS_V71
    return pf_v7::trial_field_blocks_sha256(E_root, B_root);
#else
    Digest dt("total-midpoint-E-B-blocks-v1"); dt.text(E_root); dt.text(B_root); return dt.hex();
#endif
}

DeviceStep::DeviceStep(std::shared_ptr<Device> d, const Model& m, const MotionContext& c, std::size_t capacity)
    : p_(std::make_unique<Impl>(std::move(d), m, c, capacity)) {}
DeviceStep::~DeviceStep() = default;
pf_v7::gpu::Forms& DeviceStep::forms() { return *p_->forms; }
std::size_t DeviceStep::capacity() const noexcept { return p_->N; }
std::size_t DeviceStep::carriers() const noexcept { return p_->n; }
const CensusV2& DeviceStep::census_v2() const noexcept { return p_->v2; }
const DeviceStepTimes& DeviceStep::times() const noexcept { return p_->times; }
void DeviceStep::reset_times() noexcept { p_->times = DeviceStepTimes{}; }
std::string DeviceStep::compile_json() const {
    return "{\"step\":" + hall_gpu::compile_json(p_->mod.record) + ",\"forms\":" + hall_gpu::compile_json(p_->forms->compile_record()) +
           ",\"push\":" + hall_gpu::compile_json(p_->push->compile_record()) + "}";
}

// ------------------------------------------------------------------ census
void DeviceStep::begin(const CensusSource& src, std::span<const Newborn> newborns, double begin_s, double end_s, double dt, ResidentField& field,
                       const SlotMeta& accepted) {
    auto& I = *p_;
    CpuTimer cpu{I.times.census_cpu};
    auto t0 = Clock::now();
    I.src = src; I.begin_s = begin_s; I.end_s = end_s; I.dt = dt; I.tab = field.tables();
    I.pending_identity.clear(); I.pending_census.clear(); I.pending_newborn.clear(); I.births.clear(); I.nb_index_host.clear();
    I.order.clear(); I.rows_host = nullptr; I.rows_host_size = 0; I.v2 = CensusV2{};
    const bool device_rows = src.device_rows != 0;
    const std::size_t n = device_rows ? src.device_count : src.rows.size();
    const std::size_t cap = I.ctx.maximum_carriers;
    const char* cap_text = cap == 100000 ? "motion census exceeds sealed native carrier cap" : "motion census exceeds the configured Hall carrier cap";
    // newborn list checks of MotionAssembly::flights (list order; deferred to the motion admission)
    std::vector<std::uint64_t> nb_sorted;
    {
        std::set<std::uint64_t> seen;
        for (const auto& b : newborns) {
            if (I.pending_newborn.empty() && !seen.insert(b.native_id).second) I.pending_newborn = "newborn listed twice";
            if (I.pending_newborn.empty() && !(std::isfinite(b.birth_time_s) && b.birth_time_s >= begin_s && b.birth_time_s < end_s))
                I.pending_newborn = "newborn birth time outside the step or without remaining flight";
            seen.insert(b.native_id);
        }
        nb_sorted.assign(seen.begin(), seen.end());
    }
    // ---- host gather (host census): 80-byte rows in ascending ID, identity of all rows and of the residents
    std::size_t residents = 0;
    long bad_res = -1, bad_all = -1;
    if (!device_rows) {
        I.rows_host = src.rows.data(); I.rows_host_size = src.rows.size();
        bool sorted = true;
        for (std::size_t i = 1; i < n && sorted; ++i) sorted = src.rows[i - 1].native.id < src.rows[i].native.id;
        if (!sorted) {   // the v3.3 order(): stable ascending-ID order
            I.order.reserve(n);
            for (const auto& r : src.rows) I.order.push_back(&r);
            std::stable_sort(I.order.begin(), I.order.end(), [](auto a, auto b) { return a->native.id < b->native.id; });
        }
        need(n <= I.N, "census exceeds the device step capacity");
        auto* hr = I.h_rows.as<unsigned char>();
        std::uint64_t prior_all = 0, prior_res = 0;
        std::size_t nbp = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const auto& r = I.row(i);
            std::memcpy(hr + 80 * i, &r.native, 80);
            const std::uint64_t id = r.native.id;
            while (nbp < nb_sorted.size() && nb_sorted[nbp] < id) ++nbp;
            const bool resident = !(nbp < nb_sorted.size() && nb_sorted[nbp] == id);
            const bool texts = !r.packet_id.empty() && !r.state.empty();
            if (bad_all < 0 && !(id > prior_all && texts)) bad_all = long(i);
            prior_all = id;
            if (resident) {
                ++residents;
                if (bad_res < 0 && !(id > prior_res && texts)) bad_res = long(i);
                prior_res = id;
            }
        }
        if (!src.packet_keys.empty()) {   // native census_v7 tables: ascending native ID (the sorted row order)
            need(src.packet_keys.size() == n && src.state_indices.size() == n, "native census v2 tables misaligned with the census rows");
            std::memcpy(I.h_keys.pointer, src.packet_keys.data(), n * 8); std::memcpy(I.h_states.pointer, src.state_indices.data(), n * 8);
        }
    } else residents = n - nb_sorted.size();
    I.n = n;
    // census v2 state table (native's, else the sorted states of the census: the v3.3 rule)
    if (!src.states.empty()) I.v2.states.assign(src.states.begin(), src.states.end());
    else if (!device_rows) { std::set<std::string> u; for (std::size_t i = 0; i < n; ++i) u.insert(I.row(i).state); I.v2.states.assign(u.begin(), u.end()); }
    // packet keys and state indices from the texts when native's tables are absent (v3.3 rule, host)
    if (!device_rows && src.packet_keys.empty() && n) {
        try {
            std::vector<pf_v7::HostCensusEntry> entries; entries.reserve(n);
            for (std::size_t i = 0; i < n; ++i) { const auto& r = I.row(i); entries.push_back({r.packet_id, r.state, reinterpret_cast<const pf_v7::Carrier*>(&r.native)}); }
            const auto records = pf_v7::census_records(entries, I.v2.states);
            for (std::size_t i = 0; i < n; ++i) { I.h_keys.as<std::uint64_t>()[i] = records[i].packet_key; I.h_states.as<std::uint64_t>()[i] = records[i].state_index; }
        } catch (const std::exception& e) { I.pending_census = e.what(); }
    }
    if (n > cap) I.pending_identity = cap_text;
    else if (bad_all >= 0) I.pending_identity = "motion missing/duplicate native/packet identity";
    // step-1 identity refusals of the residents (charge_of_census -> require_census_identity)
    need(residents <= cap, cap_text);
    if (bad_res >= 0) throw std::runtime_error("motion missing/duplicate native/packet identity");
    I.times.census_gather += ms_since(t0); t0 = Clock::now();
    // ---- device: SoA carriers, newborn flights, census charge of the residents, census v2, K0
    const std::uint64_t rows_d = device_rows ? src.device_rows : I.rows.address;
    const std::uint64_t keys_d = device_rows ? src.device_keys : I.keys.address, states_d = device_rows ? src.device_states : I.states.address;
    if (!device_rows && n) {
        I.s.upload(I.rows, 0, I.h_rows, 0, n * 80);
        I.s.upload(I.keys, 0, I.h_keys, 0, n * 8); I.s.upload(I.states, 0, I.h_states, 0, n * 8);
    }
    const std::uint64_t soa = I.push->carriers_address(), fl = I.push->flights_address();
    std::uint64_t n64 = n, idsa = I.ids.address;
    if (n) I.launch("k34_rows_to_soa", n, {&n64, const_cast<std::uint64_t*>(&rows_d), &idsa, const_cast<std::uint64_t*>(&soa)});
    I.push->use_device_carriers(n, true);
    std::uint64_t flight = fl, face = fl + 8 * n, newborn = fl + 12 * n, axis = I.axis_c.address, tau = I.tau_c.address;
    double dtv = dt, bs = begin_s, es = end_s;
    if (n) I.launch("k34_flights_init", n, {&n64, &flight, &face, &newborn, &axis, &tau, &dtv, &bs});
    const std::size_t B = newborns.size();
    if (B) {
        need(B <= I.N, "newborns exceed the device step capacity");
        auto* hb = I.h_nb.as<std::uint64_t>();
        for (std::size_t j = 0; j < B; ++j) { hb[j] = newborns[j].native_id; std::memcpy(&hb[B + j], &newborns[j].birth_time_s, 8); }
        I.s.upload(I.nb_ids, 0, I.h_nb, 0, B * 8); I.s.upload(I.nb_tau, 0, I.h_nb, B * 8, B * 8);
        std::uint64_t B64 = B, nbi = I.nb_ids.address, nbx = I.nb_idx.address, nbt = I.nb_tau.address, nbc = I.nb_code.address, nbr = I.nb_rank.address,
                      nbf = I.nb_face.address;
        I.launch("k34_find_ids", B, {&B64, &nbi, &idsa, &n64, &nbx});
        I.launch("k34_newborn_place", B, {&B64, &nbx, &nbt, &es, &I.pol, const_cast<std::uint64_t*>(&soa), &n64, &flight, &face, &newborn, &axis, &tau,
                                          &nbc, &nbr, &nbf});
        I.s.download(I.h_nb, 2 * B * 8, I.nb_idx, 0, B * 8);
        I.s.download(I.h_nb, 3 * B * 8, I.nb_code, 0, B * 4);
        I.s.download(I.h_nb, 3 * B * 8 + B * 4, I.nb_face, 0, B * 8);
    }
    // census charge of the residents (newborn mask), against the accepted free charge (slot 0)
    {
        std::uint32_t sentinel = std::uint32_t(I.g.nn);
        std::uint64_t mask = newborn, ka = I.ckey_a.address, pa = I.cpay_a.address, hi = I.chi.address, lo = I.clo.address, cc = I.ccode.address,
                      cf = I.cfail.address;
        if (n) I.launch("k34_records_charge", n, {&n64, const_cast<std::uint64_t*>(&rows_d), &I.hm, &I.pol, &mask, &sentinel, &ka, &pa, &hi, &lo, &cc, &cf});
        I.charge_records(8 * n, I.ckey_a.address, I.cpay_a.address, I.ckey_b.address, I.cpay_b.address, I.chi.address, I.clo.address,
                         I.c_hi.address, I.c_lo.address, I.c_gross.address);
        const auto acc = field.free_addresses(0);
        std::uint32_t nn = std::uint32_t(I.g.nn);
        std::uint64_t ch = I.c_hi.address, cl = I.c_lo.address, cg = I.c_gross.address, ah = acc[0], al = acc[1], ag = acc[2], gf = I.fail_nn.address,
                      sf = I.fail_nn2.address;
        I.launch("k34_census_nodes", I.g.nn, {&nn, &ch, &cl, &cg, &ah, &al, &ag, &gf, &sf});
        I.min_u32(I.cfail.address, n, R_CENSUS + offsetof(S34Census, charge_fail));
        I.min_u32(I.fail_nn.address, std::size_t(I.g.nn), R_CENSUS + offsetof(S34Census, guard_fail));
        I.min_u32(I.fail_nn2.address, std::size_t(I.g.nn), R_CENSUS + offsetof(S34Census, same_fail));
    }
    // resident neutral carriers (native v7 neutral interval bound)
    {
        std::uint64_t f = I.fail_n.address;
        if (n) I.launch("k34_neutral_fail", n, {&n64, const_cast<std::uint64_t*>(&soa), &newborn, &f});
        I.min_u32(I.fail_n.address, n, R_CENSUS + offsetof(S34Census, neutral_fail));
    }
    // census v2 (forms): pack, check, carriers root
    if (n) {
        I.forms->census_pack(n, rows_d, 10, keys_d, states_d, I.crec.address);
        I.forms->census_check(n, I.crec.address, I.v2.states.size(), I.chk_scratch.address, I.res(R_CHECK));
    } else I.s.zero(I.d_res, R_CHECK, 16);
    I.forms->block_list(pf7::section_census, n, I.crec.address, I.cscratch.address, I.res(R_CROOT));
    { std::uint64_t k0a = I.k0.address; if (n) I.launch("k34_k0", n, {&n64, const_cast<std::uint64_t*>(&soa), &k0a}); }
    I.s.download(I.h_res, 0, I.d_res, 0, R_SUMS);
    I.s.synchronize("device census");
    I.times.census_device += ms_since(t0); t0 = Clock::now();
    if (!n) { std::memset(I.h_res.as<char>() + R_CHECK, 0xff, 8); }
    // ---- step-1 refusals in the v3.3 order
    const auto& C = I.at<S34Census>(R_CENSUS);
    if (C.charge_fail != S34_NONE) {
        std::uint32_t code = 0;
        I.s.download(I.h_res, R_BYTES - 8, I.ccode, std::size_t(C.charge_fail) * 4, 4);
        I.s.synchronize("charge code");
        std::memcpy(&code, I.h_res.as<char>() + R_BYTES - 8, 4);
        throw std::runtime_error(hall_gpu::push::code_message(int(code)));
    }
    if (C.guard_fail != S34_NONE) throw std::runtime_error("free charge/gross operand invalid");
    if (C.same_fail != S34_NONE) throw std::runtime_error("accepted charge does not match actual canonical census");
    (void)field.admit_device(0, accepted, I.c_hi.address, I.c_lo.address, I.c_gross.address);
    I.any_resident_neutral = C.neutral_fail != S34_NONE;
    // census v2 digest and deferred check refusal
    if (I.pending_census.empty()) {
        const std::uint64_t* chk = reinterpret_cast<const std::uint64_t*>(I.h_res.as<char>() + R_CHECK);
        if (n && chk[0] != ~std::uint64_t(0))
            I.pending_census = chk[1] == pf7::census_divergent ? "motion divergent Count2 packet carriers" : "motion missing/duplicate native/packet identity";
        I.v2.states_sha256 = pf_v7::census_states_sha256(I.v2.states);
        I.v2.carriers_root = hex_words(reinterpret_cast<const std::uint32_t*>(I.h_res.as<char>() + R_CROOT));
        I.v2.sha256 = pf_v7::census_v2_sha256(n, I.v2.states_sha256, I.v2.carriers_root);
    }
    // newborn records (ascending ID) and the deferred flight / start-rule refusals (census order)
    if (B) {
        const auto* idx = reinterpret_cast<const std::uint64_t*>(I.h_nb.as<char>() + 2 * B * 8);
        const auto* code = reinterpret_cast<const std::uint32_t*>(I.h_nb.as<char>() + 3 * B * 8);
        const auto* fa = reinterpret_cast<const std::int32_t*>(I.h_nb.as<char>() + 3 * B * 8 + B * 4);
        std::vector<std::size_t> order(B);
        for (std::size_t j = 0; j < B; ++j) order[j] = j;
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return idx[a] < idx[b]; });
        std::set<std::uint64_t> placed;
        std::size_t present = 0;
        std::string first_rule;
        for (std::size_t j : order) {
            if (idx[j] >= n) continue;
            if (!placed.insert(idx[j]).second) continue;   // a duplicate ID (refused above)
            ++present;
            if (first_rule.empty() && code[j] != 0)
                first_rule = code[j] == 1 ? "newborn birth time outside the step or without remaining flight" : hgp::sink_code_message(int(code[j]) - 100);
            NewbornRecord nb;
            nb.native_id = newborns[j].native_id; nb.birth_time_s = newborns[j].birth_time_s; nb.flight_s = end_s - newborns[j].birth_time_s;
            nb.origin = fa[2 * j]; nb.axis = fa[2 * j + 1];
            native::GapParticle p{};
            std::string packet, state;
            if (!device_rows) { const auto& r = I.row(idx[j]); p = r.native; packet = r.packet_id; state = r.state; }
            else { p = I.device_row(idx[j]); if (src.packet_of) std::tie(packet, state) = src.packet_of(p.id); }
            nb.start_m = p.position_m; nb.packet_id = packet; nb.state = state; nb.charge_C = p.charge_C; nb.weight = p.weight;
            nb.birth_weighted_kinetic_J = p.weight * native::pic_kinetic_energy_J(p);
            I.births.push_back(std::move(nb)); I.nb_index_host.push_back(idx[j]);
        }
        if (I.pending_newborn.empty()) {
            if (!first_rule.empty()) I.pending_newborn = first_rule;
            else if (present != nb_sorted.size()) I.pending_newborn = "newborn native ID absent from the census";
        }
    }
    I.times.census_gather += ms_since(t0);
}

// ------------------------------------------------------------------ motion
DeviceStep::Iterate DeviceStep::motion(const hall_driver::MotionBinding& binding, std::uint64_t E, std::uint64_t B, const VerifiedField& trial, bool birth_scope) {
    auto& I = *p_;
    CpuTimer cpu{I.times.motion_cpu};
    auto t0 = Clock::now();
    const auto& c = I.ctx; const auto& g = I.g;
    const std::size_t n = I.n, nn = std::size_t(g.nn), ne = std::size_t(g.ne), nb = std::size_t(g.nb);
    std::uint64_t n64 = n;
    // ---- device work of the iterate (refusals resolved after the synchronize, in the v3.3 order)
    std::uint64_t ea = E, ba = B, fe = I.fail_ne.address, fb = I.efail_b.address;
    std::uint64_t ne64 = ne, nb64 = nb;
    I.launch("k34_fail_nonfinite", ne, {&ne64, &ea, &fe});
    I.launch("k34_fail_nonfinite", nb, {&nb64, &ba, &fb});
    I.min_u32(fe, ne, R_ITER + offsetof(S34Results, field_nonfinite));
    I.min_u32(fb, nb, R_ITER + offsetof(S34Results, field_nonfinite_b));
    I.push->set_field_device(E, B);
    I.run = I.push->run_device(binding.dt_s);
    const auto& R = I.run;
    const std::uint64_t soa = I.push->carriers_address(), fl = I.push->flights_address();
    std::uint64_t flight = fl, newborn = fl + 12 * n;
    const std::uint64_t dke = R.outs + 12 * n * 8, current = R.current, initial = R.nodes + nn * 8, gross = R.nodes + 2 * nn * 8, endpoint = R.nodes + 3 * nn * 8;
    const std::uint64_t ahi = R.nodes + 4 * nn * 8, alo = R.nodes + 5 * nn * 8, agr = R.nodes + 6 * nn * 8;
    // fold: exact sums
    {
        S34SumArgs a{}; a.n = n; a.w = soa + 8 * n * 8; a.q = soa + 6 * n * 8; a.dke = dke; a.k0 = I.k0.address;
        I.sum(0, a, 1.0, S34_SUM_K);
        I.sum(1, a, 1.0, S34_SUM_GK);
        I.sum(2, a, 1.0, S34_SUM_CK);
        I.sum(3, a, 1.0, S34_SUM_Q);
        S34SumArgs f{}; f.n = ne; f.a = E; f.b = current;
        I.sum(4, f, binding.dt_s, S34_SUM_W);
        I.sum(5, f, binding.dt_s, S34_SUM_GW);
    }
    // continuity, gross maximum, current finiteness, E^, endpoint free-charge guard, record operands
    {
        std::uint32_t nn32 = std::uint32_t(nn), ne32 = std::uint32_t(ne);
        std::uint64_t ntp = I.tab.ntp, nte = I.tab.nte, nts = I.tab.nts, elen = I.tab.elen, ea2 = I.tab.ea, eb2 = I.tab.eb, cont = I.cont.address;
        double dtv = binding.dt_s;
        std::uint64_t cur = current, ini = initial, end = endpoint;
        I.launch("k34_continuity", nn, {&nn32, &ntp, &nte, &nts, &cur, &elen, &end, &ini, &dtv, &cont});
        I.max_f64(I.cont.address, nn, R_ITER + offsetof(S34Results, continuity_C));
        I.max_f64(gross, nn, R_ITER + offsetof(S34Results, gross_max));
        std::uint64_t f2 = I.fail_ne.address;
        I.launch("k34_fail_nonfinite", ne, {&ne64, &cur, &f2});
        I.min_u32(f2, ne, R_ITER + offsetof(S34Results, current_nonfinite));
        std::uint64_t mk = I.marks.address, eb3 = I.ebound.address, Ev = E;
        I.launch("k34_node_marks", nn, {&nn32, &ntp, &nte, &cur, &mk});
        I.launch("k34_edge_bound", ne, {&ne32, &ea2, &eb2, &mk, &Ev, &eb3});
        I.max_f64(I.ebound.address, ne, R_ITER + offsetof(S34Results, field_bound));
        std::uint64_t h = ahi, l = alo, gr = agr, ff = I.fail_nn.address;
        I.launch("k34_free_guard", nn, {&nn32, &h, &l, &gr, &ff});
        I.min_u32(I.fail_nn.address, nn, R_ITER + offsetof(S34Results, free_guard_fail));
        std::uint64_t outs = R.outs, sk = R.sinks_i, sr = R.sinks_r, af = I.fail_n.address;
        double sx = -c.options.width_m / 2, sy = -c.options.gap_m / 2;
        if (n) I.launch("k34_assemble_check", n, {&n64, const_cast<std::uint64_t*>(&soa), &outs, &sk, &sr, &flight, &sx, &sy, &af});
        I.min_u32(I.fail_n.address, n, R_ITER + offsetof(S34Results, assemble_fail));
    }
    // deposit support (x87 emulation, sealed carrier order per node)
    {
        std::uint64_t outs = R.outs, a1 = I.t1.address, a2 = I.t2.address, bx = I.box.address, sc = I.scount.address, sf = I.sfail.address, scd = I.scode.address;
        if (n) I.launch("k34_support_carriers", n, {&n64, &I.geo, const_cast<std::uint64_t*>(&soa), &outs, &a1, &a2, &bx, &sc, &sf, &scd});
        I.min_u32(I.sfail.address, n, R_ITER + offsetof(S34Results, support_fail));
        if (n) I.P().exclusive_scan(I.scount.address, n, I.soff.address, I.scan_scratch.address, I.res(R_ITER + offsetof(S34Results, support_records)));
        else I.s.zero(I.d_res, R_ITER + offsetof(S34Results, support_records), 4);
        I.s.download(I.h_res, R_ITER + offsetof(S34Results, support_records), I.d_res, R_ITER + offsetof(S34Results, support_records), 4);
    }
    const int total_mark = I.s.mark("support records");
    // record flags and the compacted sink records (O(sinks) download)
    {
        std::uint64_t sk = R.sinks_i, sv = I.surv.address, sn = I.sinkf.address, bn = I.bornf.address;
        if (n) {
            I.launch("k34_record_flags", n, {&n64, &sk, &newborn, &sv, &sn, &bn});
            I.P().exclusive_scan(I.surv.address, n, I.off_surv.address, I.scan_scratch.address, I.res(R_ITER + offsetof(S34Results, survivors)));
            I.P().exclusive_scan(I.sinkf.address, n, I.off_sink.address, I.scan_scratch.address, I.res(R_ITER + offsetof(S34Results, sinks)));
            I.P().exclusive_scan(I.bornf.address, n, I.off_born.address, I.scan_scratch.address, I.res(R_ITER + offsetof(S34Results, births)));
            std::uint64_t idsa = I.ids.address, outs = R.outs, sa = R.sinks_i + 4 * n, sr = R.sinks_r, tau = I.tau_c.address, os = I.off_surv.address,
                          ok = I.off_sink.address, re = I.rec_end.address, rp = I.rec_path.address, rs = I.rec_sink.address, si = I.sink_index.address,
                          sx = I.sink_extra.address;
            int main_records = 0;
            I.launch("k34_record_build", n, {&n64, &idsa, const_cast<std::uint64_t*>(&soa), &outs, &sk, &sa, &sr, &flight, &tau, &os, &ok, &re, &rp, &rs, &si,
                                             &main_records, &sx});
        } else I.s.zero(I.d_res, R_ITER + offsetof(S34Results, survivors), 12);
    }
    // downloads: results, sums, push status and maxima, the first SPEC sink records
    I.s.download(I.h_res, R_ITER, I.d_res, R_ITER, 256);
    I.s.download(I.h_res, R_SUMS, I.d_res, R_SUMS, S34_SUMS * XS);
    I.s.download(I.h_res, R_PUSH_STATUS, R.status_buf, 64, 32);
    I.s.download(I.h_res, R_PUSH_MAX, R.maxima_buf, 0, 24);
    const std::size_t spec = std::min(SPEC, n);
    if (spec) {
        I.s.download(I.h_res, R_SINKS, I.rec_sink, 0, spec * 136);
        I.s.download(I.h_res, R_SINK_INDEX, I.sink_index, 0, spec * 8);
        I.s.download(I.h_res, R_SINK_EXTRA, I.sink_extra, 0, spec * 16);
    }
    // support records: size known after the scan; emit, sort and the node walk
    I.s.synchronize("support records");   // blocking wait (an event wait busy-waits on the host)
    (void)total_mark;
    std::size_t total = I.at<S34Results>(R_ITER).support_records;
    if (total > I.support_cap) I.grow_support(total + total / 4 + 1024);
    {
        int nx = g.nx, ny = g.ny;
        std::uint64_t bx = I.box.address, sc = I.scount.address, so = I.soff.address, ka = I.skey_a.address, pa = I.spay_a.address;
        if (n) I.launch("k34_support_emit", n, {&n64, &nx, &ny, &bx, &sc, &so, &ka, &pa});
        std::uint64_t keys_s = ka, pays_s = pa;
        I.s.zero(I.bounds, 0, 2 * (nn + 1) * 4);
        if (total) {
            const bool in_b = I.P().radix_sort(ka, pa, I.skey_b.address, I.spay_b.address, total, I.node_bits, I.ssort_scratch.address);
            keys_s = in_b ? I.skey_b.address : ka; pays_s = in_b ? I.spay_b.address : pa;
            I.P().segment_bounds(keys_s, total, I.bounds.address, I.bounds.address + (nn + 1) * 4);
        }
        std::uint32_t nn32 = std::uint32_t(nn);
        std::uint64_t b0 = I.bounds.address, b1 = I.bounds.address + (nn + 1) * 4, a1 = I.t1.address, a2 = I.t2.address, out = I.support_out.address;
        I.launch("k34_support_nodes", nn, {&nn32, &I.geo, &pays_s, &b0, &b1, &a1, &a2, &out});
    }
    const double host_before_wait = ms_since(t0);
    auto tw = Clock::now();
    I.s.synchronize("device motion");
    const double wait = ms_since(tw);
    auto t1 = Clock::now();
    I.times.device_push_ms = 0;
    for (const auto& ph : I.push->device_phases()) if (ph.name == "gpu_total_run") I.times.device_push_ms += ph.ms;
    I.s.clear_markers();
    // ---- refusals in the order of MotionAssembly::admit, GpuMotionProvider::motion and fold_record
    const auto& Rr = I.at<S34Results>(R_ITER);
    const auto* st = reinterpret_cast<const std::uint32_t*>(I.h_res.as<char>() + R_PUSH_STATUS);
    const auto* mx = reinterpret_cast<const double*>(I.h_res.as<char>() + R_PUSH_MAX);
    const auto sums = reinterpret_cast<const XsHost*>(I.h_res.as<char>() + R_SUMS);
    hall::require(binding.phase == hall_driver::ParticlePhase::full_time_position_and_proper_velocity, "motion requires owner-certified full-time phase");
    hall::require(std::isfinite(binding.begin_s) && binding.begin_s >= 0 && std::isfinite(binding.end_s) && std::isfinite(binding.dt_s) && binding.dt_s > 0 &&
                      binding.dt_s <= c.dt_max && binding.end_s > binding.begin_s && binding.end_s == binding.begin_s + binding.dt_s,
                  "motion common interval/10 fs guard");
    for (const auto* x : {&binding.staged_arena_sha256, &binding.trial_field_sha256, &binding.mesh_sha256, &binding.material_sha256, &binding.species_sha256,
                          &binding.common_clock_sha256, &binding.kinetic_policy_sha256})
        hall::require(digest(*x), "motion requires sealed stage identities");
    hall::require(Rr.field_nonfinite == S34_NONE && Rr.field_nonfinite_b == S34_NONE, "motion total-field dimensions/finiteness");
    hall::require(trial.digest == binding.trial_field_sha256 && binding.mesh_sha256 == c.mesh_sha256, "motion changed trial field or mesh binding");
    hall::require(binding.kinetic_policy_sha256 == c.policy_sha256, "motion kinetic policy changed independently of field mesh");
    if (!I.pending_identity.empty()) throw std::runtime_error(I.pending_identity);
    if (!I.pending_census.empty()) throw std::runtime_error(I.pending_census);
    if (c.forms_v7 && !(binding.dt_s <= c.neutral_stream_dt_max) && I.any_resident_neutral) throw std::runtime_error("neutral authoritative interval");
    if (!I.pending_newborn.empty()) throw std::runtime_error(I.pending_newborn);
    // push and fold
    namespace hp = hall_gpu::push;
    auto carrier_text = [](std::uint32_t code) -> std::string { return hall_gpu::push_v3::code_message_v3(int(code)); };
    if (st[0] != S34_NONE) throw std::runtime_error(carrier_text(st[1]));
    const double plasma_dt = mx[0];
    if (!(std::isfinite(plasma_dt) && plasma_dt <= .1)) throw std::runtime_error("motion full-census plasma guard");
    if (st[2] != S34_NONE) throw std::runtime_error(carrier_text(st[3]));
    hall_driver::MotionBudget bu{};
    std::uint32_t form_status = sums[S34_SUM_K].status | sums[S34_SUM_W].status;
    bu.kinetic_change_J = sums[S34_SUM_K].value; bu.field_work_J = sums[S34_SUM_W].scaled;
    bu.work_defect_J = bu.kinetic_change_J - bu.field_work_J;
    if (form_status != 0) throw std::runtime_error("canonical motion work mismatch");
    WorkScales ws;
    ws.status |= sums[S34_SUM_GW].status;
    ws.field_work_gross_J = sums[S34_SUM_GW].scaled;
    ws.field_bound_V_m = Rr.field_bound;
    if (!std::isfinite(ws.field_bound_V_m)) ws.status |= pf7::sum_nonfinite_input;
    ws.kinetic_gross_J = sums[S34_SUM_GK].value; ws.status |= sums[S34_SUM_GK].status;
    ws.kinetic_content_J = sums[S34_SUM_CK].value; ws.status |= sums[S34_SUM_CK].status;
    ws.charge_abs_C = sums[S34_SUM_Q].value; ws.status |= sums[S34_SUM_Q].status;
    ws.tolerance_J = motion_work_v1::tolerance(ws);
    ws.v31_tolerance_J = 1e-20 + 1e-10 * (std::abs(bu.kinetic_change_J) + std::abs(bu.field_work_J));
    if (!motion_work_v1::admits(ws, bu.work_defect_J)) throw std::runtime_error("canonical motion work mismatch");
    bu.continuity_C = Rr.continuity_C;
    const double gross_max = Rr.gross_max;
    if (!(Rr.current_nonfinite == S34_NONE && std::isfinite(bu.continuity_C) && bu.continuity_C <= 1e-12 * std::max(gross_max, 1e-30)))
        throw std::runtime_error("canonical motion charge continuity mismatch");
    if (st[4] != S34_NONE) throw std::runtime_error(carrier_text(st[5]));
    if (Rr.free_guard_fail != S34_NONE) throw std::runtime_error("free charge/gross operand invalid");
    bu.plasma_dt = plasma_dt; bu.maximum_gyro = mx[1]; bu.maximum_displacement_cells = mx[2];
    if (Rr.assemble_fail != S34_NONE) throw std::runtime_error("device canonical path start/endpoint operands");
    hall::require(Rr.survivors + Rr.sinks == n, "device motion sink/endpoint partition");
    hall::require(std::isfinite(bu.work_defect_J) && bu.work_defect_J == bu.kinetic_change_J - bu.field_work_J && motion_work_v1::admits(ws, bu.work_defect_J) &&
                      ws.tolerance_J == motion_work_tolerance(ws) && bu.plasma_dt >= 0 && bu.plasma_dt <= .1 && bu.maximum_gyro >= 0 && bu.maximum_gyro <= .05 &&
                      bu.maximum_displacement_cells >= 0 && bu.maximum_displacement_cells <= .2 && std::isfinite(bu.continuity_C) && bu.continuity_C >= 0,
                  "device motion budget outside the sealed guards");
    I.support_fail = Rr.support_fail == S34_NONE ? -1 : long(Rr.support_fail);
    if (I.support_fail >= 0) {
        std::uint32_t code = 0;
        I.s.download(I.h_res, R_BYTES - 8, I.scode, std::size_t(I.support_fail) * 4, 4);
        I.s.synchronize("support code");
        std::memcpy(&code, I.h_res.as<char>() + R_BYTES - 8, 4);
        I.support_code = int(code);
    }
    I.survivors = Rr.survivors; I.sinks = Rr.sinks; I.nbirths = Rr.births;
    // ---- sinks (O(sinks)): the canonical records downloaded above, Hall-side values on the host
    const std::size_t S = I.sinks;
    const char* sink_base = I.h_res.as<char>() + R_SINKS;
    const char* index_base = I.h_res.as<char>() + R_SINK_INDEX;
    const char* extra_base = I.h_res.as<char>() + R_SINK_EXTRA;
    std::vector<std::uint64_t> more;
    if (S > spec) {
        more.resize(S * (17 + 1 + 2));
        hall_gpu::HostBuffer tmp = I.s.alloc_pinned(S * 160, "sink records");
        I.s.download(tmp, 0, I.rec_sink, 0, S * 136);
        I.s.download(tmp, S * 136, I.sink_index, 0, S * 8);
        I.s.download(tmp, S * 144, I.sink_extra, 0, S * 16);
        I.s.synchronize("sink records");
        std::memcpy(more.data(), tmp.pointer, S * 160);
        I.s.release(tmp);
        sink_base = reinterpret_cast<const char*>(more.data()); index_base = sink_base + S * 136; extra_base = sink_base + S * 144;
    }
    MotionAssembly::DeviceParts parts;
    parts.binding = binding; parts.census_sha256 = I.v2.sha256; parts.budget = bu; parts.work = ws; parts.birth_scope = birth_scope;
    const native::PicVector shift{-c.options.width_m / 2, -c.options.gap_m / 2, 0};
    parts.sinks.reserve(S);
    for (std::size_t j = 0; j < S; ++j) {
        std::uint64_t w[17]; std::memcpy(w, sink_base + 136 * j, 136);
        std::uint64_t ci; std::memcpy(&ci, index_base + 8 * j, 8);
        double extra[2]; std::memcpy(extra, extra_base + 16 * j, 16);
        auto dbl = [](std::uint64_t x) { return std::bit_cast<double>(x); };
        native::GapParticle census{};
        std::string packet, state;
        if (I.rows_host) { const auto& r = I.row(ci); census = r.native; packet = r.packet_id; state = r.state; }
        else { census = I.device_row(ci); if (I.src.packet_of) std::tie(packet, state) = I.src.packet_of(census.id); }
        SinkRecord kk;
        kk.native_id = w[0]; kk.packet_id = packet; kk.state = state; kk.boundary = int(w[2]);
        kk.collected = kk.boundary == 0 || kk.boundary == 1 || kk.boundary == 3; kk.kind = kk.collected ? 1 : 2;
        kk.conductor = kk.boundary == 3 ? 2 : (kk.collected ? kk.boundary : -1);
        kk.axis = int(std::int64_t(w[3])); kk.plane_m = dbl(w[4]); kk.fraction = dbl(w[5]);
        for (int d = 0; d < 3; ++d) {
            kk.segment_m[d] = dbl(w[6 + d]); kk.hit_point_m[d] = dbl(w[9 + d]); kk.incident_proper_velocity_m_s[d] = dbl(w[12 + d]);
            kk.field_point_m[d] = kk.hit_point_m[d] + shift[d];
        }
        kk.flight_s = dbl(w[15]); kk.discarded_flight_s = dbl(w[16]);
        native::GapParticle e = census; e.position_m = kk.hit_point_m; e.proper_velocity_m_s = kk.incident_proper_velocity_m_s;
        const double K1 = native::pic_kinetic_energy_J(e), K0 = native::pic_kinetic_energy_J(census);
        kk.charge_C = census.charge_C; kk.mass_kg = census.mass_kg; kk.weight = census.weight; kk.kinetic_J = K1; kk.weighted_kinetic_J = census.weight * K1;
        kk.initial_weighted_kinetic_J = census.weight * K0; kk.gyro = extra[0]; kk.displacement_cells = extra[1];
        parts.sinks.push_back(std::move(kk));
    }
    parts.births = I.births;
    Iterate it{MotionAssembly::device(c, std::move(parts))};
    it.current = R.current; it.after_hi = ahi; it.after_lo = alo; it.after_gross = agr; it.support = I.support_out.address;
    it.host_ms = host_before_wait + ms_since(t1); it.wait_ms = wait;
    I.times.motion_host += it.host_ms; I.times.motion_wait += wait; ++I.times.iterations;
    return it;
}
void DeviceStep::support_check() const {
    const auto& I = *p_;
    if (I.support_fail >= 0) throw std::runtime_error(support_text(I.support_code));
}

// ------------------------------------------------------------------ record of the accepted iterate
const DeviceStep::DeviceRecordView& DeviceStep::record_view() const noexcept { return p_->view; }
hall_gpu::Session& DeviceStep::session() { return p_->s; }
void DeviceStep::seal(hall_driver::MotionReceipt& r, bool host_arrays, std::uint64_t trial_E, std::uint64_t trial_B, FieldDigestForm form) {
    auto& I = *p_;
    const bool blocks = form == FieldDigestForm::blocks_v71;
#if !HALL_FORMS_V71
    need(!blocks, "the block-form field digests (native v7.1) need the forms library of native v7.1");
#endif
    need(!blocks || (trial_E && trial_B), "block-form trial-field digest without the trial arrays");
    CpuTimer cpu{I.times.record_cpu};
    auto t0 = Clock::now();
    const std::size_t n = I.n, nn = std::size_t(I.g.nn), ne = std::size_t(I.g.ne);
    const auto& R = I.run;
    std::uint64_t n64 = n;
    const std::uint64_t soa = I.push->carriers_address(), fl = I.push->flights_address();
    std::uint64_t flight = fl, newborn = fl + 12 * n;
    if (n) {
        std::uint64_t idsa = I.ids.address, outs = R.outs, sk = R.sinks_i, sa = R.sinks_i + 4 * n, sr = R.sinks_r, tau = I.tau_c.address, os = I.off_surv.address,
                      ok = I.off_sink.address, re = I.rec_end.address, rp = I.rec_path.address, rs = I.rec_sink.address, si = I.sink_index.address, sx = I.sink_extra.address;
        int main_records = 1;
        I.launch("k34_record_build", n, {&n64, &idsa, const_cast<std::uint64_t*>(&soa), &outs, &sk, &sa, &sr, &flight, &tau, &os, &ok, &re, &rp, &rs, &si,
                                         &main_records, &sx});
        std::uint64_t ob = I.off_born.address, rb = I.rec_birth.address, bi = I.birth_index.address;
        I.launch("k34_birth_build", n, {&n64, &idsa, &newborn, &flight, &tau, &ob, &rb, &bi});
    }
    const std::uint64_t counts[4] = {I.survivors, n, I.sinks, I.nbirths};
    const std::uint64_t arrays[4] = {I.rec_end.address, I.rec_path.address, I.rec_sink.address, I.rec_birth.address};
    for (std::uint32_t sct = 0; sct < 4; ++sct) I.forms->block_list(sct, counts[sct], arrays[sct], I.bl_scratch[sct].address, I.res(R_SROOT + 32 * sct));
    I.field_roots(R.current, trial_E, trial_B);
    I.s.download(I.h_res, R_SROOT, I.d_res, R_SROOT, 4 * 32);
    I.s.download(I.h_res, R_FROOT, I.d_res, R_FROOT, 96);
    // host arrays of the record header (array digests) and, on request, the endpoint and path records
    // the current is downloaded for its streaming digest or for host arrays (block form: its root suffices)
    const bool host_current = !blocks || host_arrays;
    std::vector<double> cur(host_current ? ne : 0), ini(nn), end(nn);
    const hall_gpu::HostBuffer& tmp = I.h_rec;
    if (host_current) I.s.download(tmp, 0, R.current_buf, 0, ne * 8);
    I.s.download(tmp, ne * 8, R.nodes_buf, nn * 8, nn * 8);
    I.s.download(tmp, (ne + nn) * 8, R.nodes_buf, 3 * nn * 8, nn * 8);
    const std::size_t o2 = (ne + 2 * nn) * 8;
    if (host_arrays && n) {
        if (I.survivors) I.s.download(tmp, o2, I.rec_end, 0, I.survivors * 80);
        I.s.download(tmp, o2 + I.survivors * 80, I.rec_path, 0, n * 64);
    }
    I.s.synchronize("record");
    if (host_current) std::memcpy(cur.data(), tmp.pointer, ne * 8);
    std::memcpy(ini.data(), tmp.as<char>() + ne * 8, nn * 8);
    std::memcpy(end.data(), tmp.as<char>() + (ne + nn) * 8, nn * 8);
    std::vector<native::GapParticle> eps; std::vector<hall_driver::CarrierPath> pts;
    if (host_arrays) {
        eps.resize(I.survivors); pts.resize(n);
        if (I.survivors) std::memcpy(eps.data(), tmp.as<char>() + o2, I.survivors * 80);
        if (n) std::memcpy(pts.data(), tmp.as<char>() + o2 + I.survivors * 80, n * 64);
    }
    MotionAssembly::set_device_arrays(r, std::move(cur), std::move(ini), std::move(end), host_arrays ? &eps : nullptr, host_arrays ? &pts : nullptr);
    const auto* froot = reinterpret_cast<const std::uint32_t*>(I.h_res.as<char>() + R_FROOT);
    const std::string J_root = hex_words(froot);
    std::string trial_blocks;
    if (trial_E && trial_B) trial_blocks = trial_blocks_digest(hex_words(froot + 8), hex_words(froot + 16));
    if (blocks) MotionAssembly::set_trial_field_sha256(r, trial_blocks);   // the binding's trial field in the block form
    auto h = MotionAssembly::header_v7(r);
    if (blocks) h.current_A_m_sha256 = current_blocks_digest(J_root);
    pf_v7::MotionSectionsV7 sec;
    sec.endpoints = counts[0]; sec.paths = counts[1]; sec.sinks = counts[2]; sec.births = counts[3];
    const auto* roots = reinterpret_cast<const std::uint32_t*>(I.h_res.as<char>() + R_SROOT);
    sec.endpoints_root = hex_words(roots); sec.paths_root = hex_words(roots + 8); sec.sinks_root = hex_words(roots + 16); sec.births_root = hex_words(roots + 24);
    I.sealed_record = pf_v7::motion_record_v7_sha256(h, sec);
    I.sealed_array_digests = {h.current_A_m_sha256, h.endpoint_charge_C_sha256, h.initial_charge_C_sha256};
    I.sealed_current_root = J_root;
    I.sealed_form = form; I.sealed_trial_blocks = trial_blocks;
    ++I.seal_generation;
    {   // native v7.1 device lane view
        auto& v = I.view;
        v = DeviceRecordView{};
        v.endpoints = I.rec_end.address; v.paths = I.rec_path.address; v.sinks = I.rec_sink.address; v.births = I.rec_birth.address;
        v.n_endpoints = counts[0]; v.n_paths = counts[1]; v.n_sinks = counts[2]; v.n_births = counts[3];
        v.E = trial_E; v.B = trial_B; v.J = R.current; v.n_E = ne; v.n_B = std::size_t(I.g.nb);
        v.J_root = I.sealed_current_root;
        v.current_blocks_sha256 = current_blocks_digest(v.J_root);
        if (trial_E && trial_B) {
            v.E_root = hex_words(froot + 8);
            v.B_root = hex_words(froot + 16);
            v.trial_field_blocks_sha256 = trial_blocks;
        }
        v.seal_generation = I.seal_generation;
    }
    MotionAssembly::set_record_v7(r, I.sealed_record, I.sealed_array_digests);
    I.times.record += ms_since(t0);
}

// ------------------------------------------------------------------ preflight
bool DeviceStep::preflight(const hall_driver::MotionReceipt& m, std::span<const native::ArenaRateEntry> census, ResidentField& field, const SlotMeta& accepted) {
    auto& I = *p_;
    CpuTimer cpu{I.times.preflight_cpu};
    try {
        const std::size_t n = I.n, nn = std::size_t(I.g.nn), ne = std::size_t(I.g.ne);
        // the census must be the one of the proposal (host span identity, or the device rows)
        if (I.src.device_rows == 0 && (census.data() != I.rows_host || census.size() != I.rows_host_size)) return false;
        std::uint64_t n64 = n;
        const std::uint64_t rows_d = I.src.device_rows ? I.src.device_rows : I.rows.address;
        const std::uint64_t keys_d = I.src.device_rows ? I.src.device_keys : I.keys.address, states_d = I.src.device_rows ? I.src.device_states : I.states.address;
        // census v2 root recomputed from the device census
        if (n) I.forms->census_pack(n, rows_d, 10, keys_d, states_d, I.crec.address);
        I.forms->block_list(pf7::section_census, n, I.crec.address, I.cscratch.address, I.res(R_CROOT));
        // section roots recomputed from the record arrays of the seal
        const std::uint64_t counts[4] = {I.survivors, n, I.sinks, I.nbirths};
        const std::uint64_t arrays[4] = {I.rec_end.address, I.rec_path.address, I.rec_sink.address, I.rec_birth.address};
        for (std::uint32_t sct = 0; sct < 4; ++sct) I.forms->block_list(sct, counts[sct], arrays[sct], I.bl_scratch[sct].address, I.res(R_SROOT + 32 * sct));
        const bool blocks = I.sealed_form == FieldDigestForm::blocks_v71;
        // current root (and, block form, the trial E/B roots) recomputed on the device
        I.field_roots(I.run.current, blocks ? I.view.E : 0, blocks ? I.view.B : 0);
        // census charge of the residents recomputed (device rows, newborn mask), against the accepted free charge
        const std::uint64_t fl = I.push->flights_address(), newborn = fl + 12 * n;
        std::uint32_t sentinel = std::uint32_t(nn);
        std::uint64_t mask = newborn, ka = I.ckey_a.address, pa = I.cpay_a.address, hi = I.chi.address, lo = I.clo.address, cc = I.ccode.address, cf = I.cfail.address;
        if (n) I.launch("k34_records_charge", n, {&n64, const_cast<std::uint64_t*>(&rows_d), &I.hm, &I.pol, &mask, &sentinel, &ka, &pa, &hi, &lo, &cc, &cf});
        I.charge_records(8 * n, I.ckey_a.address, I.cpay_a.address, I.ckey_b.address, I.cpay_b.address, I.chi.address, I.clo.address,
                         I.c_hi.address, I.c_lo.address, I.c_gross.address);
        const auto acc = field.free_addresses(0);
        std::uint32_t nn32 = std::uint32_t(nn);
        {
            std::uint64_t ch = I.c_hi.address, cl = I.c_lo.address, cg = I.c_gross.address, ah = acc[0], al = acc[1], ag = acc[2], gf = I.fail_nn.address,
                          sf = I.fail_nn2.address;
            I.launch("k34_census_nodes", nn, {&nn32, &ch, &cl, &cg, &ah, &al, &ag, &gf, &sf});
            I.min_u32(I.cfail.address, n, R_CENSUS + offsetof(S34Census, charge_fail));
            I.min_u32(I.fail_nn.address, nn, R_CENSUS + offsetof(S34Census, guard_fail));
            I.min_u32(I.fail_nn2.address, nn, R_CENSUS + offsetof(S34Census, same_fail));
        }
        // endpoint free charge recomputed from the endpoint records, compared bitwise with the candidate slot
        {
            const std::uint64_t ns = I.survivors;
            std::uint64_t rec = I.rec_end.address, none = 0;
            if (ns) I.launch("k34_records_charge", ns, {const_cast<std::uint64_t*>(&ns), &rec, &I.hm, &I.pol, &none, &sentinel, &ka, &pa, &hi, &lo, &cc, &cf});
            I.min_u32(I.cfail.address, ns, R_PRE + offsetof(S34Preflight, endpoint_charge_fail));
            I.charge_records(8 * ns, I.ckey_a.address, I.cpay_a.address, I.ckey_b.address, I.cpay_b.address, I.chi.address, I.clo.address,
                             I.e_hi.address, I.e_lo.address, I.e_gross.address);
            const auto cand = field.free_addresses(1);
            const std::uint64_t mine[3] = {I.e_hi.address, I.e_lo.address, I.e_gross.address};
            for (int q = 0; q < 3; ++q) {
                std::uint64_t a = mine[q], b = cand[q], f = q == 0 ? I.fail_nn.address : I.fail_nn2.address;
                std::uint64_t nn64 = nn;
                I.launch("k34_differ_u64", nn, {&nn64, &a, &b, &f});
                I.min_u32(f, nn, R_PRE + offsetof(S34Preflight, differ) + 4 * std::size_t(q));
            }
        }
        I.s.download(I.h_res, 0, I.d_res, 0, R_SUMS);
        I.s.download(I.h_res, R_FROOT, I.d_res, R_FROOT, 96);
        I.s.synchronize("device preflight");
        const auto& C = I.at<S34Census>(R_CENSUS);
        const auto& Pf = I.at<S34Preflight>(R_PRE);
        if (C.charge_fail != S34_NONE || C.guard_fail != S34_NONE || C.same_fail != S34_NONE) return false;
        if (Pf.endpoint_charge_fail != S34_NONE || Pf.differ[0] != S34_NONE || Pf.differ[1] != S34_NONE || Pf.differ[2] != S34_NONE) return false;
        (void)field.admit_device(0, accepted, I.c_hi.address, I.c_lo.address, I.c_gross.address);
        const auto census_root = hex_words(reinterpret_cast<const std::uint32_t*>(I.h_res.as<char>() + R_CROOT));
        if (pf_v7::census_v2_sha256(n, pf_v7::census_states_sha256(I.v2.states), census_root) != m.census_sha256()) return false;
        const auto* froot = reinterpret_cast<const std::uint32_t*>(I.h_res.as<char>() + R_FROOT);
        if (hex_words(froot) != I.sealed_current_root) return false;
        auto h = MotionAssembly::header_v7(m);
        if (blocks) {
            // block form: the trial field of the binding and the current from the device roots recomputed above; a host
            // copy of the current (host record arrays) must hash to the same block-form digest
            if (trial_blocks_digest(hex_words(froot + 8), hex_words(froot + 16)) != I.sealed_trial_blocks ||
                m.binding().trial_field_sha256 != I.sealed_trial_blocks)
                return false;
            h.current_A_m_sha256 = current_blocks_digest(hex_words(froot));
#if HALL_FORMS_V71
            if (!m.current_A_m().empty() && pf_v7::current_blocks_sha256(m.current_A_m()) != h.current_A_m_sha256) return false;
#endif
        }
        if (h.current_A_m_sha256 != I.sealed_array_digests[0] || h.endpoint_charge_C_sha256 != I.sealed_array_digests[1] ||
            h.initial_charge_C_sha256 != I.sealed_array_digests[2] || m.record_v7_array_digests() != I.sealed_array_digests)
            return false;
        pf_v7::MotionSectionsV7 sec;
        sec.endpoints = counts[0]; sec.paths = counts[1]; sec.sinks = counts[2]; sec.births = counts[3];
        const auto* roots = reinterpret_cast<const std::uint32_t*>(I.h_res.as<char>() + R_SROOT);
        sec.endpoints_root = hex_words(roots); sec.paths_root = hex_words(roots + 8); sec.sinks_root = hex_words(roots + 16); sec.births_root = hex_words(roots + 24);
        return pf_v7::motion_record_v7_sha256(h, sec) == m.record_v7_sha256();
    } catch (const std::exception&) {
        return false;
    }
}

// ------------------------------------------------------------------ block digests, qualification access
std::string DeviceStep::block_root(int which, std::uint64_t values, std::size_t n) {
    auto& I = *p_;
    I.f64_root(std::uint32_t(which), values, n, R_FROOT + 96);
    I.s.download(I.h_res, R_FROOT + 96, I.d_res, R_FROOT + 96, 32);
    I.s.synchronize("block root");
    return hex_words(reinterpret_cast<const std::uint32_t*>(I.h_res.as<char>() + R_FROOT + 96));
}
std::string DeviceStep::block_root_forms(int which, std::uint64_t values, std::size_t n) {
#if HALL_FORMS_V71
    auto& I = *p_;
    need(which >= 0 && which <= 2 && n <= std::size_t(std::max(I.g.ne, I.g.nb)), "block root: section 0..2 within the mesh");
    const std::uint32_t sct[3] = {pf7::section_field_E, pf7::section_field_B, pf7::section_current_J};
    const pf_v7::gpu::Forms::BlockListJob job{sct[which], n, values, I.res(R_FROOT + 128)};
    I.forms->block_lists({&job, 1}, I.flists_scratch.address);
    I.s.download(I.h_res, R_FROOT + 128, I.d_res, R_FROOT + 128, 32);
    I.s.synchronize("block root (forms)");
    return hex_words(reinterpret_cast<const std::uint32_t*>(I.h_res.as<char>() + R_FROOT + 128));
#else
    (void)which; (void)values; (void)n;
    throw std::runtime_error("Forms::block_lists needs the forms library of native v7.1");
#endif
}
Vec DeviceStep::support_download() {
    auto& I = *p_;
    Vec out(std::size_t(I.g.nn));
    hall_gpu::HostBuffer tmp = I.s.alloc_pinned(out.size() * 8, "support download");
    I.s.download(tmp, 0, I.support_out, 0, out.size() * 8);
    I.s.synchronize("support download");
    std::memcpy(out.data(), tmp.pointer, out.size() * 8);
    I.s.release(tmp);
    return out;
}
std::vector<std::uint64_t> DeviceStep::record_words(int section) {
    auto& I = *p_;
    need(section >= 0 && section <= 5, "record section 0..3, 4 (census rows) or 5 (current J)");
    const std::size_t words[6] = {10, 8, 17, 3, 10, 1};
    const std::size_t counts[6] = {I.survivors, I.n, I.sinks, I.nbirths, I.src.device_rows ? 0 : I.n, std::size_t(I.g.ne)};
    const hall_gpu::DeviceBuffer* arr[6] = {&I.rec_end, &I.rec_path, &I.rec_sink, &I.rec_birth, &I.rows, &I.run.current_buf};
    std::vector<std::uint64_t> out(counts[section] * words[section]);
    if (out.empty()) return out;
    hall_gpu::HostBuffer tmp = I.s.alloc_pinned(out.size() * 8, "record download");
    I.s.download(tmp, 0, *arr[section], 0, out.size() * 8);
    I.s.synchronize("record download");
    std::memcpy(out.data(), tmp.pointer, out.size() * 8);
    I.s.release(tmp);
    return out;
}
void DeviceStep::tamper_record_word(int section, std::size_t word, std::uint64_t xor_mask) {
    auto w = record_words(section);
    need(word < w.size(), "tamper word outside the record");
    auto& I = *p_;
    const hall_gpu::DeviceBuffer* arr[6] = {&I.rec_end, &I.rec_path, &I.rec_sink, &I.rec_birth, &I.rows, &I.run.current_buf};
    w[word] ^= xor_mask;
    hall_gpu::HostBuffer tmp = I.s.alloc_pinned(8, "tamper");
    std::memcpy(tmp.pointer, &w[word], 8);
    I.s.upload(*arr[section], word * 8, tmp, 0, 8);
    I.s.synchronize("tamper");
    I.s.release(tmp);
}

// ------------------------------------------------------------------ sparse bank transfer
SparseBank bank_transfer_sparse(const Mesh& g, const native::GapPicOptions& o, const hall_driver::MotionReceipt& r, const std::vector<bool>& exterior_support) {
    SparseBank t;
    if (r.sinks().empty() && r.births().empty()) return t;
    t.empty = false;
    struct Acc { Words cond, ext; double cg = 0, eg = 0, residue = 0; bool touched_c = false, touched_e = false; };
    std::map<int, Acc> nodes;
    auto conductor = [&](int n) { return g.node_body[n] >= 0 || g.short_node[n]; };
    // deposit_charge at the eight CIC nodes (field frame), the operations of census_charge.hpp in its order
    auto deposit = [&](bool to_conductor, const native::PicVector& xn, Words charge) {
        native::PicVector x = xn; x[0] -= o.width_m / 2; x[1] -= o.gap_m / 2;
        need(canonical(charge), "invalid deposited charge words");
        native::PicVector spacing{g.dx, g.dy, g.dz}, origin{-g.side / 2, -g.side / 2, g.zmin}, a{};
        std::array<int, 3> counts{g.nx, g.ny, g.nz}, cell{};
        for (int d = 0; d < 3; ++d) {
            need(std::isfinite(x[d]), "nonfinite charge position"); double v = (x[d] - origin[d]) / spacing[d];
            need(v >= -1e-10 && v <= counts[d] + 1e-10, "canonical charge outside mesh"); v = std::clamp(v, 0., double(counts[d]));
            cell[d] = std::min(int(std::floor(v)), counts[d] - 1); a[d] = v - cell[d];
        }
        const long double Q = magnitude(charge);
        std::array<std::pair<int, Words>, 8> slots;
        int s = 0;
        for (int kk : {0, 1}) for (int j : {0, 1}) for (int i : {0, 1}) {
            const int n = g.node(cell[0] + i, cell[1] + j, cell[2] + kk);
            slots[std::size_t(s++)] = {n, charge * (i ? a[0] : 1 - a[0]) * (j ? a[1] : 1 - a[1]) * (kk ? a[2] : 1 - a[2])};
        }
        std::sort(slots.begin(), slots.end(), [](const auto& x, const auto& y) { return x.first < y.first; });   // ascending node (the v3.3 scan)
        for (const auto& [n, v0] : slots) {
            const Words v = Words() + v0;   // one.net_C[n] = 0 + v
            if (v.hi == 0 && v.lo == 0) continue;
            const long double gross_x = (long double)0.0 + magnitude(v);
            double gross = double(gross_x); if ((long double)gross < gross_x) gross = std::nextafter(gross, INFINITY);
            auto& acc = nodes[n];
            if (to_conductor) { acc.cond = acc.cond + v; acc.cg += gross; acc.touched_c = true; }
            else { acc.ext = acc.ext + v; acc.eg += gross; acc.touched_e = true; }
            const bool supported = to_conductor ? conductor(n) : bool(exterior_support[std::size_t(n)]);
            if (!supported) {
                const long double m = magnitude(v);
                need(m <= 1e-12L * Q, to_conductor ? "kinetic sink/newborn conductor bank deposit outside the conductor nodes"
                                                   : "kinetic sink exterior bank deposit outside the exterior support");
                const long double x2 = (long double)acc.residue + m; double vv = double(x2); if ((long double)vv < x2) vv = std::nextafter(vv, INFINITY);
                acc.residue = vv;
            }
        }
    };
    for (const auto& kk : r.sinks()) {
        if (kk.charge_C == 0) continue;
        const auto q = Words(kk.weight) * kk.charge_C;
        need(canonical(q) && (q.hi != 0 || q.lo != 0), "kinetic sink unresolved macro-charge");
        deposit(kk.collected, kk.hit_point_m, q);
    }
    for (const auto& b : r.births()) {
        if (b.charge_C == 0) continue;
        need(b.origin >= 0, "charged interior newborn needs the source path current and charge impulses (not admitted by this Hall version)");
        const auto q = Words(b.weight) * b.charge_C;
        need(canonical(q) && (q.hi != 0 || q.lo != 0), "newborn unresolved macro-charge");
        deposit(true, b.start_m, -q);
    }
    for (const auto& [n, acc] : nodes) {
        t.nodes.push_back(n); t.conductor_C.push_back(acc.cond); t.exterior_C.push_back(acc.ext);
        const long double x = (long double)acc.cg + acc.eg; double v = double(x); if ((long double)v < x) v = std::nextafter(v, INFINITY);
        t.gross_C.push_back(v); t.residue_C.push_back(acc.residue); t.residue_any = t.residue_any || acc.residue != 0;
    }
    return t;
}
BankTransfer dense_of(const Mesh& g, const SparseBank& s) {
    BankTransfer t; t.empty = s.empty;
    if (s.empty) return t;
    t.conductor_C.assign(std::size_t(g.nn), Words()); t.exterior_C.assign(std::size_t(g.nn), Words()); t.gross_C.assign(std::size_t(g.nn), 0.);
    t.residue_C.assign(std::size_t(g.nn), 0.);
    for (std::size_t j = 0; j < s.nodes.size(); ++j) {
        const auto n = std::size_t(s.nodes[j]);
        t.conductor_C[n] = s.conductor_C[j]; t.exterior_C[n] = s.exterior_C[j]; t.gross_C[n] = s.gross_C[j]; t.residue_C[n] = s.residue_C[j];
    }
    return t;
}
}  // namespace hall_corrective_v3::gpu
