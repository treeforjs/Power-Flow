// ms_gpu.cpp -- host side of the device merge/split transaction (see ms_gpu.hpp).
#include "ms_gpu.hpp"

#include <algorithm>
#include <chrono>
#include <cstring>

namespace merge_split {
namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kKernels[] = {
#embed "ms_kernels.cu"
    , 0};
const char kCore[] = {
#embed "../include/ms_core.h"
    , 0};
const char kPhilox[] = {
#embed "philox_core.h"
    , 0};
#pragma GCC diagnostic pop

void need(bool ok, const char* message) {
    if (!ok) throw Error(message);
}
unsigned blocks(std::size_t n, unsigned b = 256) { return unsigned((n + b - 1) / b); }
constexpr unsigned NQ = 10;

// Device view types with the same layout as the kernel parameter structs.
struct DevOutputs {
    std::uint64_t id, sp, fl, x[3], u[3], w;
};
struct DevReloc {
    std::uint64_t id, group, start[3], D[3], charge;
};
static_assert(sizeof(MsView) == 12 * 8 && sizeof(DevOutputs) == 10 * 8 && sizeof(DevReloc) == 9 * 8);
static_assert(sizeof(MsGroup) == 64 && sizeof(MsMergeResult) == 72);
}  // namespace

struct GpuTransaction::Impl {
    hall_gpu::Session& s;
    hall_gpu::Module m;
    std::unique_ptr<hall_gpu::Primitives> prim;
    hall_gpu::Kernel classify, gather_u32, gather_sorted, seg_flags, form_count, form_write, kinematics, emit, scatter,
        compact, ledger_chunks, ledger_final;
    explicit Impl(hall_gpu::Session& session) : s(session) {
        hall_gpu::CompileRequest req;
        req.main = {"ms_kernels.cu", kKernels};
        req.headers = {{"ms_core.h", kCore}, {"philox_core.h", kPhilox}, hall_gpu::Primitives::device_source()};
        m = s.compile_and_load(req);
        prim = std::make_unique<hall_gpu::Primitives>(s, m);
        classify = s.kernel(m, "ms_classify");
        gather_u32 = s.kernel(m, "ms_gather_u32");
        gather_sorted = s.kernel(m, "ms_gather_sorted");
        seg_flags = s.kernel(m, "ms_segment_flags");
        form_count = s.kernel(m, "ms_form_count");
        form_write = s.kernel(m, "ms_form_write");
        kinematics = s.kernel(m, "ms_kinematics");
        emit = s.kernel(m, "ms_emit");
        scatter = s.kernel(m, "ms_scatter_u32");
        compact = s.kernel(m, "ms_compact_removed");
        ledger_chunks = s.kernel(m, "ms_ledger_chunks");
        ledger_final = s.kernel(m, "ms_ledger_final");
    }
};

GpuTransaction::GpuTransaction(hall_gpu::Session& s) : p_(std::make_unique<Impl>(s)) {}
GpuTransaction::~GpuTransaction() = default;
const hall_gpu::Module& GpuTransaction::module() const { return p_->m; }

Result GpuTransaction::run(const MsParams& P, const Carriers& C, GpuTimings* timings) {
    using Clock = std::chrono::steady_clock;
    const auto wall0 = Clock::now();
    validate(P);
    auto& s = p_->s;
    auto& prim = *p_->prim;
    const std::size_t n = C.size();
    need(n >= 1 && n < (std::size_t(1) << 31), "device transaction carrier count");
    need(C.species.size() == n && C.flags.size() == n && C.w.size() == n, "carrier array sizes");
    for (int d = 0; d < 3; ++d) need(C.x[d].size() == n && C.u[d].size() == n, "carrier array sizes");
    const unsigned S = P.n_species;
    const std::size_t ncells = std::size_t(P.nx) * P.ny * P.nz, nseg = ncells * S;
    const unsigned un = unsigned(n), unseg = unsigned(nseg);
    std::vector<hall_gpu::DeviceBuffer> owned;
    std::vector<hall_gpu::HostBuffer> pinned;
    std::size_t dev_bytes = 0, pin_bytes = 0;
    auto dalloc = [&](std::size_t bytes, const char* tag) {
        auto b = s.alloc(std::max<std::size_t>(bytes, 8), tag);
        owned.push_back(b);
        dev_bytes += b.bytes;
        return b.address;
    };
    // the allocation that contains a device address (sub-array views of the input block included)
    auto buf_of = [&](std::uint64_t addr) -> const hall_gpu::DeviceBuffer& {
        for (const auto& b : owned)
            if (addr >= b.address && addr < b.address + b.bytes) return b;
        throw Error("device address outside the transaction allocations");
    };
    auto palloc = [&](std::size_t bytes, const char* tag) {
        auto b = s.alloc_pinned(std::max<std::size_t>(bytes, 8), tag);
        pinned.push_back(b);
        pin_bytes += b.bytes;
        return b;
    };
    Result R;
    try {
        R.params_sha256 = params_sha256(P);
        R.counters.carriers = n;
        // ---- upload: params + 12 input arrays through one pinned staging buffer
        const std::size_t f8 = n * 8, f4 = n * 4;
        const std::size_t in_bytes = sizeof(MsParams) + f8 /*id*/ + 2 * f4 /*sp, fl*/ + 7 * f8 /*x,u,w*/;
        auto stage = palloc(in_bytes, "ms-upload");
        char* h = static_cast<char*>(stage.pointer);
        std::size_t off = 0;
        auto put = [&](const void* src, std::size_t bytes) { std::memcpy(h + off, src, bytes); off += bytes; };
        put(&P, sizeof(MsParams));
        put(C.id.data(), f8);
        put(C.species.data(), f4);
        put(C.flags.data(), f4);
        for (int d = 0; d < 3; ++d) put(C.x[d].data(), f8);
        for (int d = 0; d < 3; ++d) put(C.u[d].data(), f8);
        put(C.w.data(), f8);
        const std::uint64_t din = dalloc(in_bytes, "ms-input");
        std::uint64_t dP = din, did = din + sizeof(MsParams), dsp = did + f8, dfl = dsp + f4, dx[3], du[3];
        dx[0] = dfl + f4;
        dx[1] = dx[0] + f8;
        dx[2] = dx[1] + f8;
        du[0] = dx[2] + f8;
        du[1] = du[0] + f8;
        du[2] = du[1] + f8;
        std::uint64_t dw = du[2] + f8;
        // sizeof(MsParams) is a multiple of 8, so every double array stays 8-byte aligned.
        static_assert(sizeof(MsParams) % 8 == 0);
        s.synchronize("merge-split start");
        s.clear_markers();
        const int m0 = s.mark("upload");
        s.upload(owned.back(), 0, stage, 0, in_bytes);
        const int m1 = s.mark("uploaded");
        // ---- classify, stable sorts (species|bin, then cell), gather
        const std::uint64_t cell = dalloc(f4, "cell"), bin = dalloc(f4, "bin"), ke = dalloc(f8, "ke"), key_lo = dalloc(f4, "key_lo"),
                            status = dalloc(f4, "status"), keyA = dalloc(f4, "keyA"), keyB = dalloc(f4, "keyB"),
                            payA = dalloc(f4, "payA"), payB = dalloc(f4, "payB"),
                            scratch = dalloc(std::max(hall_gpu::Primitives::sort_scratch(n), hall_gpu::Primitives::scan_scratch(std::max(n, nseg))) * 4 + 64, "scratch"),
                            scalars = dalloc(64, "scalars");
        {
            void* a[] = {&dP, (void*)&un, &did, &dsp, &dx[0], &dx[1], &dx[2], &du[0], &du[1], &du[2], &dw, (void*)&cell, (void*)&bin,
                         (void*)&ke, (void*)&key_lo, (void*)&status};
            s.launch(p_->classify, {blocks(n)}, {256}, 0, a);
        }
        // status total (exclusive scan of 0/1 flags into keyB as scratch output; total in scalars[0])
        prim.exclusive_scan(status, n, keyB, scratch, scalars);
        prim.iota(payA, n);
        // copy key_lo into keyA (gather with the identity permutation)
        {
            void* a[] = {(void*)&un, (void*)&payA, (void*)&key_lo, (void*)&keyA};
            s.launch(p_->gather_u32, {blocks(n)}, {256}, 0, a);
        }
        const unsigned lo_bits = 32;
        bool in_b = prim.radix_sort(keyA, payA, keyB, payB, n, lo_bits, scratch);
        std::uint64_t perm = in_b ? payB : payA;
        std::uint64_t perm_other = in_b ? payA : payB;
        std::uint64_t kcur = in_b ? keyB : keyA, kother = in_b ? keyA : keyB;
        {
            void* a[] = {(void*)&un, (void*)&perm, (void*)&cell, (void*)&kcur};
            s.launch(p_->gather_u32, {blocks(n)}, {256}, 0, a);
        }
        const unsigned cell_bits = ((P.cell_bits + 7) / 8) * 8;
        const bool in_b2 = prim.radix_sort(kcur, perm, kother, perm_other, n, cell_bits, scratch);
        if (in_b2) perm = perm_other;
        // sorted arrays
        const std::uint64_t sid = dalloc(f8, "sid"), ssp = dalloc(f4, "ssp"), sfl = dalloc(f4, "sfl"), sbin = dalloc(f4, "sbin"),
                            sw = dalloc(f8, "sw"), ske = dalloc(f8, "ske"), sseg = dalloc(f4, "sseg");
        std::uint64_t sx[3], su[3];
        for (int d = 0; d < 3; ++d) { sx[d] = dalloc(f8, "sx"); su[d] = dalloc(f8, "su"); }
        {
            void* a[] = {(void*)&un, (void*)&S, &perm, &did, &dsp, &dfl, &dx[0], &dx[1], &dx[2], &du[0], &du[1], &du[2], &dw,
                         (void*)&cell, (void*)&bin, (void*)&ke, (void*)&sid, (void*)&ssp, (void*)&sfl, (void*)&sbin, &sx[0], &sx[1],
                         &sx[2], &su[0], &su[1], &su[2], (void*)&sw, (void*)&ske, (void*)&sseg};
            s.launch(p_->gather_sorted, {blocks(n)}, {256}, 0, a);
        }
        MsView view{reinterpret_cast<const ms_u64*>(sid), reinterpret_cast<const ms_u32*>(ssp), reinterpret_cast<const ms_u32*>(sfl),
                    reinterpret_cast<const ms_u32*>(sbin),
                    {reinterpret_cast<const double*>(sx[0]), reinterpret_cast<const double*>(sx[1]), reinterpret_cast<const double*>(sx[2])},
                    {reinterpret_cast<const double*>(su[0]), reinterpret_cast<const double*>(su[1]), reinterpret_cast<const double*>(su[2])},
                    reinterpret_cast<const double*>(sw), reinterpret_cast<const double*>(ske)};
        // segments (cell x species)
        const std::uint64_t seg_begin = dalloc(nseg * 4, "seg_begin"), seg_end = dalloc(nseg * 4, "seg_end"),
                            seg_aux = dalloc(nseg * 4, "seg_aux"), counts = dalloc(nseg * 4, "counts"),
                            refused = dalloc(nseg * 4, "refused"), goff = dalloc(nseg * 4, "goff");
        s.zero(buf_of(seg_begin), 0, nseg * 4);
        s.zero(buf_of(seg_end), 0, nseg * 4);
        prim.segment_bounds(sseg, n, seg_begin, seg_end);
        {
            void* a[] = {(void*)&unseg, (void*)&seg_begin, (void*)&seg_end, (void*)&seg_aux};
            s.launch(p_->seg_flags, {blocks(nseg)}, {256}, 0, a);
        }
        prim.exclusive_scan(seg_aux, nseg, seg_aux, scratch, scalars + 4);   // nonempty segments -> scalars[1]
        {
            void* a[] = {&dP, (void*)&unseg, (void*)&seg_begin, (void*)&seg_end, &view, (void*)&counts, (void*)&refused};
            s.launch(p_->form_count, {blocks(nseg)}, {256}, 0, a);
        }
        prim.exclusive_scan(refused, nseg, refused, scratch, scalars + 8);   // formation refusals -> scalars[2]
        prim.exclusive_scan(counts, nseg, goff, scratch, scalars + 12);      // groups -> scalars[3]
        const int m2 = s.mark("sorted-counted");
        auto small = palloc(64, "ms-scalars");
        const hall_gpu::DeviceBuffer scalar_buf = buf_of(scalars);
        s.download(small, 0, scalar_buf, 0, 16);
        s.synchronize("merge-split counts");
        std::uint32_t sc[4];
        std::memcpy(sc, small.pointer, 16);
        need(sc[0] == 0, "invalid census row or carriers not in strictly ascending native ID order");
        R.counters.segments = sc[1];
        R.counters.formation_refused = sc[2];
        const std::size_t G = sc[3];
        const unsigned uG = unsigned(G);
        R.counters.groups = G;
        // ---- groups, kinematics, scans
        const std::uint64_t groups = dalloc(G * sizeof(MsGroup), "groups"), merge = dalloc(G * sizeof(MsMergeResult), "merge"),
                            n_out = dalloc(G * 4, "n_out"), n_rel = dalloc(G * 4, "n_reloc"), out_off = dalloc(G * 4, "out_off"),
                            rel_off = dalloc(G * 4, "reloc_off"), first_id = dalloc(G * 8, "first_id");
        const int m3 = s.mark("kinematics");
        if (G) {
            {
                void* a[] = {&dP, (void*)&unseg, (void*)&seg_begin, (void*)&seg_end, &view, (void*)&goff, (void*)&groups};
                s.launch(p_->form_write, {blocks(nseg)}, {256}, 0, a);
            }
            {
                void* a[] = {&dP, (void*)&uG, (void*)&groups, &view, (void*)&merge, (void*)&n_out, (void*)&n_rel};
                s.launch(p_->kinematics, {blocks(G)}, {256}, 0, a);
            }
        }
        prim.exclusive_scan(n_out, G, out_off, scratch, scalars + 16);   // outputs -> scalars[4]
        prim.exclusive_scan(n_rel, G, rel_off, scratch, scalars + 20);   // relocations -> scalars[5]
        const int m4 = s.mark("kinematics-done");
        s.download(small, 0, scalar_buf, 0, 24);
        s.synchronize("merge-split outputs");
        std::uint32_t sc2[6];
        std::memcpy(sc2, small.pointer, 24);
        const std::size_t NO = sc2[4], NR = sc2[5];
        const unsigned uNO = unsigned(NO);
        // ---- emission, removed IDs, ledger
        DevOutputs o{};
        o.id = dalloc(NO * 8, "out_id");
        o.sp = dalloc(NO * 4, "out_sp");
        o.fl = dalloc(NO * 4, "out_fl");
        for (int d = 0; d < 3; ++d) { o.x[d] = dalloc(NO * 8, "out_x"); o.u[d] = dalloc(NO * 8, "out_u"); }
        o.w = dalloc(NO * 8, "out_w");
        DevReloc r{};
        r.id = dalloc(NR * 8, "rel_id");
        r.group = dalloc(NR * 4, "rel_group");
        for (int d = 0; d < 3; ++d) { r.start[d] = dalloc(NR * 8, "rel_start"); r.D[d] = dalloc(NR * 8, "rel_D"); }
        r.charge = dalloc(NR * 8, "rel_charge");
        const std::uint64_t consumed = dalloc(f4, "consumed"), consumed_orig = dalloc(f4, "consumed_orig"), rem_off = dalloc(f4, "rem_off"),
                            removed = dalloc(f8, "removed");
        const hall_gpu::DeviceBuffer consumed_buf = buf_of(consumed);
        const unsigned chunks = unsigned(std::clamp<std::size_t>((std::max(n, NO) + 2047) / 2048, 1, 1024));
        const unsigned nq = S * NQ;
        const std::uint64_t partial = dalloc(std::size_t(nq) * chunks * MS_ACC_DIGITS * 8, "ledger_partial"),
                            total = dalloc(std::size_t(nq) * MS_ACC_DIGITS * 8, "ledger_total");
        const int m5 = s.mark("emit");
        s.zero(consumed_buf, 0, f4);
        if (G) {
            void* a[] = {&dP, (void*)&uG, (void*)&groups, (void*)&merge, (void*)&out_off, (void*)&rel_off, &view, (void*)&first_id, &o, &r,
                         (void*)&consumed};
            s.launch(p_->emit, {blocks(G)}, {256}, 0, a);
        }
        {
            void* a[] = {(void*)&un, &perm, (void*)&consumed, (void*)&consumed_orig};
            s.launch(p_->scatter, {blocks(n)}, {256}, 0, a);
        }
        prim.exclusive_scan(consumed_orig, n, rem_off, scratch, scalars + 24);   // removed -> scalars[6]
        {
            void* a[] = {(void*)&un, (void*)&consumed_orig, (void*)&rem_off, &did, (void*)&removed};
            s.launch(p_->compact, {blocks(n)}, {256}, 0, a);
        }
        {
            void* a[] = {&dP, (void*)&un, &view, (void*)&consumed, (void*)&uNO, &o.sp, &o.w, &o.u[0], &o.u[1], &o.u[2], (void*)&chunks,
                         (void*)&partial};
            s.launch(p_->ledger_chunks, {blocks(chunks, 64), nq}, {64}, 0, a);
        }
        {
            void* a[] = {(void*)&nq, (void*)&chunks, (void*)&partial, (void*)&total};
            s.launch(p_->ledger_final, {blocks(std::size_t(nq) * MS_ACC_DIGITS)}, {256}, 0, a);
        }
        const int m6 = s.mark("emitted");
        // ---- download everything through one pinned buffer
        const std::size_t down = 28 + G * (sizeof(MsGroup) + sizeof(MsMergeResult) + 8) + NO * (8 + 4 + 4 + 7 * 8) +
                                 NR * (8 + 4 + 7 * 8) + f8 + std::size_t(nq) * MS_ACC_DIGITS * 8;
        auto back = palloc(down, "ms-download");
        std::size_t doff = 0;
        auto get = [&](std::uint64_t src, std::size_t bytes) {
            if (!bytes) return doff;
            const std::size_t at = doff;
            const hall_gpu::DeviceBuffer& owner = buf_of(src);
            s.download(back, at, owner, src - owner.address, bytes);
            doff += bytes;
            return at;
        };
        const int m7 = s.mark("download");
        const auto a_sc = get(scalars, 28);
        const auto a_groups = get(groups, G * sizeof(MsGroup));
        const auto a_merge = get(merge, G * sizeof(MsMergeResult));
        const auto a_first = get(first_id, G * 8);
        const auto a_oid = get(o.id, NO * 8), a_osp = get(o.sp, NO * 4), a_ofl = get(o.fl, NO * 4);
        std::size_t a_ox[3], a_ou[3], a_rs[3], a_rd[3];
        for (int d = 0; d < 3; ++d) a_ox[d] = get(o.x[d], NO * 8);
        for (int d = 0; d < 3; ++d) a_ou[d] = get(o.u[d], NO * 8);
        const auto a_ow = get(o.w, NO * 8);
        const auto a_rid = get(r.id, NR * 8), a_rg = get(r.group, NR * 4);
        for (int d = 0; d < 3; ++d) a_rs[d] = get(r.start[d], NR * 8);
        for (int d = 0; d < 3; ++d) a_rd[d] = get(r.D[d], NR * 8);
        const auto a_rc = get(r.charge, NR * 8);
        const auto a_tot = get(total, std::size_t(nq) * MS_ACC_DIGITS * 8);
        // removed IDs: count is known only after the scan; download the full-length array prefix later
        const int m8 = s.mark("downloaded");
        s.synchronize("merge-split results");
        const char* b = static_cast<const char*>(back.pointer);
        std::uint32_t sc3[7];
        std::memcpy(sc3, b + a_sc, 28);
        const std::size_t NREM = sc3[6];
        auto take = [&](auto& vec, std::size_t at, std::size_t count) {
            vec.resize(count);
            if (count) std::memcpy(vec.data(), b + at, count * sizeof(vec[0]));
        };
        take(R.groups, a_groups, G);
        take(R.merge, a_merge, G);
        take(R.group_first_id, a_first, G);
        take(R.outputs.id, a_oid, NO);
        take(R.outputs.species, a_osp, NO);
        take(R.outputs.flags, a_ofl, NO);
        for (int d = 0; d < 3; ++d) { take(R.outputs.x[d], a_ox[d], NO); take(R.outputs.u[d], a_ou[d], NO); }
        take(R.outputs.w, a_ow, NO);
        take(R.reloc.id, a_rid, NR);
        take(R.reloc.group, a_rg, NR);
        for (int d = 0; d < 3; ++d) { take(R.reloc.start[d], a_rs[d], NR); take(R.reloc.D[d], a_rd[d], NR); }
        take(R.reloc.charge, a_rc, NR);
        std::vector<std::int64_t> tot(std::size_t(nq) * MS_ACC_DIGITS);
        std::memcpy(tot.data(), b + a_tot, tot.size() * 8);
        // removed IDs (second small download of the compacted prefix)
        R.removed.resize(NREM);
        if (NREM) {
            const auto a_rem = get(removed, NREM * 8);
            s.synchronize("merge-split removed");
            std::memcpy(R.removed.data(), b + a_rem, NREM * 8);
        }
        R.counters.inputs_consumed = NREM;
        detail::fill_group_counters(R, P);
        need(R.counters.outputs == NO && R.counters.relocations == NR, "device output counts");
        std::vector<Exact> totals(nq);
        for (unsigned q = 0; q < nq; ++q) {
            MsAcc a;
            for (int d = 0; d < MS_ACC_DIGITS; ++d) a.d[d] = tot[std::size_t(q) * MS_ACC_DIGITS + d];
            totals[q] = Exact::of(a);
        }
        detail::fill_ledger(R, P, totals);
        R.digest_sha256 = sha256_hex(encode(R));
        if (timings) {
            timings->upload_ms = s.elapsed_ms(m0, m1);
            timings->sort_ms = s.elapsed_ms(m1, m2);
            timings->kinematics_ms = s.elapsed_ms(m3, m4);
            timings->emit_ledger_ms = s.elapsed_ms(m5, m6);
            timings->download_ms = s.elapsed_ms(m7, m8);
            timings->device_bytes = dev_bytes;
            timings->pinned_bytes = pin_bytes;
        }
    } catch (...) {
        for (auto& b : owned) if (b) s.release(b);
        for (auto& b : pinned) if (b) s.release(b);
        throw;
    }
    for (auto& b : owned) s.release(b);
    for (auto& b : pinned) s.release(b);
    if (timings) timings->host_wall_ms = std::chrono::duration<double, std::milli>(Clock::now() - wall0).count();
    return R;
}

}  // namespace merge_split
