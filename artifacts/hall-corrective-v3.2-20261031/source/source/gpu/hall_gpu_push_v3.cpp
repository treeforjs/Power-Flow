// Hall corrective v3 GPU motion engine. Copied from hall-gpu-resident-v1-20261030/
// push/hall_gpu_push.cpp (sealed, unchanged) with these changes only: the NVRTC
// program adds sink_core.hpp, push_abi_v3.hpp, push_bodies_v3.hpp and runs
// push_kernels_v3.cu (carrier body with first-hit truncation, endpoint reduction
// without sink records); the carrier kernel receives the sink geometry and writes
// the per-carrier sink outputs, which are downloaded with the carrier outputs.
#include "hall_gpu_push_v3.hpp"
#include "push/push_abi.hpp"
#include "push_fold_v32.hpp"   // v3.2: sealed fold with the v3.2 work tolerance
#include "push_abi_v3.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>

namespace hall_gpu::push_v3 {
namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kCore[] = {
#embed "../../../hall-gpu-resident-v1-20261030/push/push_core.hpp"
    , 0};
const char kAbi[] = {
#embed "../../../hall-gpu-resident-v1-20261030/push/push_abi.hpp"
    , 0};
const char kBodies[] = {
#embed "../../../hall-gpu-resident-v1-20261030/push/push_bodies.hpp"
    , 0};
const char kSinkCore[] = {
#embed "../sink_core.hpp"
    , 0};
const char kAbiV3[] = {
#embed "push_abi_v3.hpp"
    , 0};
const char kBodiesV3[] = {
#embed "push_bodies_v3.hpp"
    , 0};
const char kKernels[] = {
#embed "push_kernels_v3.cu"
    , 0};
#pragma GCC diagnostic pop

using Clock = std::chrono::steady_clock;
double ms_since(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }
unsigned blocks(std::size_t n, unsigned b) { return unsigned((n + b - 1) / b); }
unsigned key_bits_for(std::size_t max_key) {   // bits to represent keys 0..max_key
    unsigned bits = 1;
    while ((std::uint64_t(1) << bits) <= max_key) ++bits;
    return bits;
}
constexpr std::uint32_t none = 0xffffffffu;
}  // namespace

const char* code_message_v3(int code) { return code >= hgp::sink_field_without_kinetic && code <= hgp::sink_unsupported_bank ? hgp::sink_code_message(code) : code_message(code); }

struct MotionEngineV3::Impl {
    Session& s;
    MeshTables mt;
    Geometry geo;
    SinkSupport sup;
    EngineOptions opt;
    Module module;
    std::unique_ptr<Primitives> prim;
    Kernel k_carrier, k_emit, k_edges, k_start, k_end, k_plasma, k_final;
    std::size_t N = 0, nn = 0, ne = 0, nb = 0, n = 0, edge_cap = 0, sort_cap = 0;
    unsigned edge_bits = 0, node_bits = 0;
    // device
    DeviceBuffer carriers, field, VN, codes, fails, outs, counts, offsets;
    DeviceBuffer s_key_a, s_pay_a, s_key_b, s_pay_b, s_val, e_key_a, e_pay_a, e_key_b, e_pay_b, e_val;
    DeviceBuffer g_key_a, g_pay_a, g_key_b, g_pay_b, g_val;   // edge records
    DeviceBuffer nodes, current, bounds, scan_scratch, sort_scratch, red_scratch, small_u32, small_f64, diag_real, diag_int, sinks_i, sinks_r, flights;
    bool flights_set = false;
    // pinned
    HostBuffer h_carriers, h_field, h_out, h_nodes, h_current, h_status, h_maxes, h_codes, h_diag_real, h_diag_int, h_sinks_i, h_sinks_r, h_flights;
    std::vector<std::uint64_t> ids;
    std::vector<double> weights;
    std::vector<double> k0, charge;   // v3.2: initial kinetic energy and charge per carrier (work scales)
    bool forms_v7 = false;            // v3.2: native v7 exact sums in the fold
    bool carriers_set = false, field_set = false, field_ok = true;
    int mark_carriers0 = -1, mark_carriers1 = -1, mark_field0 = -1, mark_field1 = -1;
    std::size_t bytes = 0;

    Impl(Session& session, const MeshTables& m, const Geometry& g, const SinkSupport& u, const EngineOptions& o) : s(session), mt(m), geo(g), sup(u), opt(o) {}

    DeviceBuffer dev(std::size_t b, const char* tag) { bytes += b; return s.alloc(b, tag); }

    void allocate_edges(std::size_t records) {
        edge_cap = std::max<std::size_t>(records, 1024);
        g_key_a = dev(edge_cap * 4, "edge key a");
        g_pay_a = dev(edge_cap * 4, "edge payload a");
        g_key_b = dev(edge_cap * 4, "edge key b");
        g_pay_b = dev(edge_cap * 4, "edge payload b");
        g_val = dev(edge_cap * 8, "edge value");
        sort_cap = std::max(8 * N, edge_cap);
        sort_scratch = dev(Primitives::sort_scratch(sort_cap) * 4, "sort scratch");
    }
    void grow_edges(std::size_t records) {
        for (auto* b : {&g_key_a, &g_pay_a, &g_key_b, &g_pay_b, &g_val, &sort_scratch}) { bytes -= b->bytes; s.release(*b); }
        allocate_edges(records + records / 4 + 1024);
    }
};

MotionEngineV3::MotionEngineV3(Session& s, const MeshTables& m, const Geometry& g, const SinkSupport& u, const EngineOptions& o) : p_(new Impl(s, m, g, u, o)) {
    try {
        auto& I = *p_;
        I.N = o.capacity;
        I.nn = std::size_t(m.mesh.nn);
        I.ne = std::size_t(m.mesh.ne);
        I.nb = std::size_t(m.mesh.nb);
        if (I.N == 0 || I.N > 200000000 / 8) throw Error("MotionEngine capacity must be in 1..25e6 carriers");
        if (m.VN.size() != I.nn || m.edge_a.size() != I.ne || m.edge_b.size() != I.ne || m.edge_length.size() != I.ne)
            throw Error("MeshTables sizes disagree with the mesh");
        I.edge_bits = key_bits_for(I.ne ? I.ne - 1 : 0);
        I.node_bits = key_bits_for(I.nn);   // includes the sentinel key nn
        CompileRequest req;
        req.main = {"push_kernels_v3.cu", kKernels};
        req.headers = {{"push_core.hpp", kCore}, {"push_abi.hpp", kAbi}, {"push_bodies.hpp", kBodies}, {"sink_core.hpp", kSinkCore},
                       {"push_abi_v3.hpp", kAbiV3}, {"push_bodies_v3.hpp", kBodiesV3}, Primitives::device_source()};
        I.module = s.compile_and_load(req);
        I.prim = std::make_unique<Primitives>(s, I.module);
        I.k_carrier = s.kernel(I.module, "hgp_carrier_v3");
        I.k_emit = s.kernel(I.module, "hgp_emit_edges");
        I.k_edges = s.kernel(I.module, "hgp_reduce_edges");
        I.k_start = s.kernel(I.module, "hgp_reduce_start_nodes");
        I.k_end = s.kernel(I.module, "hgp_reduce_end_nodes_v3");
        I.k_plasma = s.kernel(I.module, "hgp_plasma");
        I.k_final = s.kernel(I.module, "hgp_finalize");
        const std::size_t N = I.N, R = 8 * N;
        I.carriers = I.dev(9 * N * 8, "carriers");
        I.field = I.dev((I.ne + I.nb) * 8, "midpoint E,B");
        I.VN = I.dev(I.nn * 8, "VN");
        I.codes = I.dev(3 * N * 4, "codes");
        I.fails = I.dev(3 * N * 4, "fail ranks");
        I.outs = I.dev(15 * N * 8, "carrier outputs");
        I.counts = I.dev(N * 4, "edge counts");
        I.offsets = I.dev(N * 4, "edge offsets");
        I.s_key_a = I.dev(R * 4, "start key a"); I.s_pay_a = I.dev(R * 4, "start payload a");
        I.s_key_b = I.dev(R * 4, "start key b"); I.s_pay_b = I.dev(R * 4, "start payload b");
        I.s_val = I.dev(3 * R * 8, "start values");
        I.e_key_a = I.dev(R * 4, "end key a"); I.e_pay_a = I.dev(R * 4, "end payload a");
        I.e_key_b = I.dev(R * 4, "end key b"); I.e_pay_b = I.dev(R * 4, "end payload b");
        I.e_val = I.dev(3 * R * 8, "end values");
        I.allocate_edges(std::size_t(double(N) * o.initial_edge_records_per_carrier));
        I.nodes = I.dev(8 * I.nn * 8, "node outputs");
        I.current = I.dev(std::max<std::size_t>(I.ne, 1) * 8, "current");
        I.bounds = I.dev(2 * std::max(I.ne, I.nn + 1) * 4, "segment bounds");
        I.scan_scratch = I.dev(Primitives::scan_scratch(N) * 4 + 64, "scan scratch");
        I.red_scratch = I.dev(Primitives::reduce_scratch(std::max(N, I.nn)) * 8 + 64, "reduce scratch");
        I.small_u32 = I.dev(64 * 4, "status words");
        I.small_f64 = I.dev(8 * 8, "status reals");
        I.sinks_i = I.dev(2 * N * 4, "sink kinds/axes");
        I.sinks_r = I.dev(7 * N * 8, "sink reals");
        I.h_sinks_i = s.alloc_pinned(2 * N * 4, "sink kinds/axes host");
        I.h_sinks_r = s.alloc_pinned(7 * N * 8, "sink reals host");
        I.flights = I.dev(N * 8 + 2 * N * 4, "newborn flights/faces/flags");
        I.h_flights = s.alloc_pinned(N * 8 + 2 * N * 4, "newborn flights host");
        if (o.diagnostics_capable) {
            I.diag_real = I.dev(std::size_t(hgp::diag_real_fields) * N * 8, "diagnostics real");
            I.diag_int = I.dev(std::size_t(hgp::diag_int_fields) * N * 4, "diagnostics int");
            I.h_codes = s.alloc_pinned(3 * N * 4, "codes host");
            I.h_diag_real = s.alloc_pinned(std::size_t(hgp::diag_real_fields) * N * 8, "diagnostics real host");
            I.h_diag_int = s.alloc_pinned(std::size_t(hgp::diag_int_fields) * N * 4, "diagnostics int host");
        }
        I.h_carriers = s.alloc_pinned(9 * N * 8, "carrier staging");
        I.h_field = s.alloc_pinned((I.ne + I.nb + I.nn) * 8, "field staging");
        I.h_out = s.alloc_pinned(15 * N * 8, "carrier outputs host");
        I.h_nodes = s.alloc_pinned(6 * I.nn * 8, "node outputs host");
        I.h_current = s.alloc_pinned(std::max<std::size_t>(I.ne, 1) * 8, "current host");
        I.h_status = s.alloc_pinned(64 * 4, "status host");
        I.h_maxes = s.alloc_pinned(8 * 8, "maxima host");
        // Static mesh data: staged in pinned memory, uploaded on the stream, fenced.
        std::memcpy(I.h_field.as<double>() + I.ne + I.nb, m.VN.data(), I.nn * 8);
        s.upload(I.VN, 0, I.h_field, (I.ne + I.nb) * 8, I.nn * 8);
        s.synchronize("static mesh upload");
    } catch (...) {
        delete p_;
        p_ = nullptr;
        throw;
    }
}

MotionEngineV3::~MotionEngineV3() {
    // Release this engine's buffers through the Session (each release drains the
    // stream first). After a session failure nothing is released here: the
    // Session's own teardown decides between orderly release and quarantine.
    if (p_ && !p_->s.failed()) {
        try {
            auto& I = *p_;
            for (auto* b : {&I.carriers, &I.field, &I.VN, &I.codes, &I.fails, &I.outs, &I.counts, &I.offsets, &I.s_key_a, &I.s_pay_a,
                            &I.s_key_b, &I.s_pay_b, &I.s_val, &I.e_key_a, &I.e_pay_a, &I.e_key_b, &I.e_pay_b, &I.e_val, &I.g_key_a,
                            &I.g_pay_a, &I.g_key_b, &I.g_pay_b, &I.g_val, &I.nodes, &I.current, &I.bounds, &I.scan_scratch,
                            &I.sort_scratch, &I.red_scratch, &I.small_u32, &I.small_f64, &I.diag_real, &I.diag_int, &I.sinks_i, &I.sinks_r, &I.flights})
                I.s.release(*b);
            for (auto* h : {&I.h_carriers, &I.h_field, &I.h_out, &I.h_nodes, &I.h_current, &I.h_status, &I.h_maxes, &I.h_codes,
                            &I.h_diag_real, &I.h_diag_int, &I.h_sinks_i, &I.h_sinks_r, &I.h_flights})
                I.s.release(*h);
        } catch (...) {
            // The session is now failed; its teardown handles the remaining resources.
        }
    }
    delete p_;
}

const CompileRecord& MotionEngineV3::compile_record() const noexcept { return p_->module.record; }
std::size_t MotionEngineV3::device_bytes() const noexcept { return p_->bytes; }

void MotionEngineV3::set_carriers(const CarrierView& c) {
    auto& I = *p_;
    const std::size_t n = c.size();
    for (auto span_size : {c.x.size(), c.y.size(), c.z.size(), c.ux.size(), c.uy.size(), c.uz.size(), c.charge.size(), c.mass.size(), c.weight.size()})
        if (span_size != n) throw Error("carrier arrays have different lengths");
    if (n > I.N) throw Error("carrier count exceeds engine capacity");
    std::uint64_t prior = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (!(c.id[i] > prior)) throw Error("carriers must be in strictly increasing native-ID order with IDs > 0");
        prior = c.id[i];
    }
    I.s.synchronize("before carrier staging");   // staging memory may still feed an earlier upload
    double* h = I.h_carriers.as<double>();
    const std::span<const double> parts[9] = {c.x, c.y, c.z, c.ux, c.uy, c.uz, c.charge, c.mass, c.weight};
    for (int k = 0; k < 9; ++k) if (n) std::memcpy(h + std::size_t(k) * n, parts[k].data(), n * 8);
    I.mark_carriers0 = I.s.mark("h2d_carriers_begin");
    I.s.upload(I.carriers, 0, I.h_carriers, 0, 9 * n * 8);
    I.mark_carriers1 = I.s.mark("h2d_carriers_end");
    I.ids.assign(c.id.begin(), c.id.end());
    I.weights.assign(c.weight.begin(), c.weight.end());
    I.charge.assign(c.charge.begin(), c.charge.end());
    I.k0.resize(n);
    for (std::size_t i = 0; i < n; ++i) I.k0[i] = pf7::kinetic_energy(c.mass[i], c.ux[i], c.uy[i], c.uz[i]);   // K0 (motion-work-gross-v1)
    I.n = n;
    I.carriers_set = true;
}

void MotionEngineV3::set_motion_forms(bool native_v7) { p_->forms_v7 = native_v7; }

void MotionEngineV3::set_flights(std::span<const double> flight, std::span<const int> face, std::span<const char> newborn) {
    auto& I = *p_;
    if (flight.empty()) { I.flights_set = false; return; }
    const std::size_t n = I.n;
    if (flight.size() != n || face.size() != n || newborn.size() != n) throw Error("newborn flight arrays must match the carriers");
    I.s.synchronize("before flight staging");
    char* h = static_cast<char*>(I.h_flights.pointer);
    std::memcpy(h, flight.data(), n * 8);
    int* f = reinterpret_cast<int*>(h + n * 8);
    for (std::size_t i = 0; i < n; ++i) { f[i] = face[i]; f[n + i] = newborn[i] ? 1 : 0; }
    I.s.upload(I.flights, 0, I.h_flights, 0, n * 8 + 2 * n * 4);
    I.flights_set = true;
}

void MotionEngineV3::set_field(std::span<const double> E, std::span<const double> B) {
    auto& I = *p_;
    I.field_ok = E.size() == I.ne && B.size() == I.nb;
    if (I.field_ok) for (double v : E) if (!std::isfinite(v)) { I.field_ok = false; break; }
    if (I.field_ok) for (double v : B) if (!std::isfinite(v)) { I.field_ok = false; break; }
    I.field_set = true;
    if (!I.field_ok) return;
    I.s.synchronize("before field staging");
    double* h = I.h_field.as<double>();
    std::memcpy(h, E.data(), I.ne * 8);
    std::memcpy(h + I.ne, B.data(), I.nb * 8);
    I.mark_field0 = I.s.mark("h2d_field_begin");
    I.s.upload(I.field, 0, I.h_field, 0, (I.ne + I.nb) * 8);
    I.mark_field1 = I.s.mark("h2d_field_end");
}

MotionViewV3 MotionEngineV3::run(double dt, bool diagnostics) {
    auto& I = *p_;
    auto& s = I.s;
    if (!I.carriers_set || !I.field_set) throw Error("run() requires carriers and a field");
    if (!(std::isfinite(dt) && dt > 0)) throw Error("run() requires a finite positive dt");
    if (diagnostics && !I.opt.diagnostics_capable) throw Error("engine was built without diagnostic outputs");
    const auto wall0 = Clock::now();
    MotionViewV3 v;
    const std::size_t n = I.n, nn = I.nn, ne = I.ne;
    v.n = n;
    v.diagnostics = diagnostics;
    if (!I.field_ok) {
        v.refusal = {true, "motion total-field dimensions/finiteness", "field", -1, 0, 0};
        return v;
    }
    const hgp::Mesh& g = I.mt.mesh;
    auto D = [](const DeviceBuffer& b, std::size_t byte_offset = 0) { return b.address + byte_offset; };
    const std::size_t R = 8 * I.N;

    // ---- carrier kernel
    hgp::CarrierArgsV3 av3{};
    hgp::CarrierArgs& a = av3.v1;
    a.g = g;
    a.p = I.geo.policy;
    for (int d = 0; d < 3; ++d) a.off[d] = I.geo.native_to_field[d];
    a.dt = dt;
    a.n = hgp_u32(n);
    a.node_sentinel = hgp_u32(nn);
    a.diagnostics = diagnostics ? 1 : 0;
    const double* cb = I.carriers.as<const double>();
    a.x = cb; a.y = cb + n; a.z = cb + 2 * n; a.ux = cb + 3 * n; a.uy = cb + 4 * n; a.uz = cb + 5 * n;
    a.q = cb + 6 * n; a.m = cb + 7 * n; a.w = cb + 8 * n;
    a.E = I.field.as<const double>();
    a.B = a.E + ne;
    int* codes = I.codes.as<int>();
    a.code1 = codes; a.code2 = codes + n; a.code3 = codes + 2 * n;
    hgp_u32* fails = I.fails.as<hgp_u32>();
    a.fail1 = fails; a.fail2 = fails + n; a.fail3 = fails + 2 * n;
    double* outs = I.outs.as<double>();
    a.endpoint = outs; a.start = outs + 6 * n; a.disp = outs + 9 * n; a.dke = outs + 12 * n; a.gyro = outs + 13 * n; a.dcells = outs + 14 * n;
    a.edge_count = I.counts.as<hgp_u32>();
    a.start_key = I.s_key_a.as<hgp_u32>(); a.start_pay = I.s_pay_a.as<hgp_u32>();
    double* sv = I.s_val.as<double>();
    a.start_omega = sv; a.start_init = sv + R; a.start_gross = sv + 2 * R;
    a.end_key = I.e_key_a.as<hgp_u32>(); a.end_pay = I.e_pay_a.as<hgp_u32>();
    double* ev = I.e_val.as<double>();
    a.end_charge = ev; a.end_hi = ev + R; a.end_lo = ev + 2 * R;
    a.diag_real = diagnostics ? I.diag_real.as<double>() : nullptr;
    a.diag_int = diagnostics ? I.diag_int.as<int>() : nullptr;
    {
        auto& sg = av3.sg;
        sg.p = I.geo.policy;
        for (int d = 0; d < 3; ++d) sg.off[d] = I.geo.native_to_field[d];
        for (int k = 0; k < 5; ++k) { for (int d = 0; d < 3; ++d) { sg.lo[k][d] = I.sup.lo[k][d]; sg.hi[k][d] = I.sup.hi[k][d]; } sg.supported[k] = I.sup.supported[k]; }
        sg.aligned = I.sup.aligned ? 1 : 0;
        av3.sink_kind = I.sinks_i.as<int>(); av3.sink_axis = I.sinks_i.as<int>() + n; av3.sink_real = I.sinks_r.as<double>();
        if (I.flights_set) {
            av3.flight = I.flights.as<const double>();
            av3.face = reinterpret_cast<const int*>(I.flights.as<const char>() + n * 8);
            av3.newborn = reinterpret_cast<const int*>(I.flights.as<const char>() + n * 8) + n;
        } else { av3.flight = nullptr; av3.face = nullptr; av3.newborn = nullptr; }
    }
    const int m0 = s.mark("run_begin");
    { void* args[] = {&av3}; s.launch(I.k_carrier, {blocks(n, 128)}, {128}, 0, args); }
    const int m1 = s.mark("carrier_kernel");

    // ---- edge record count (scan of 24 x segments); read back early
    const std::uint64_t edge_total = D(I.small_u32, 0), min_fail = D(I.small_u32, 16), status = D(I.small_u32, 64);
    I.prim->exclusive_scan(D(I.counts), n, D(I.offsets), D(I.scan_scratch), edge_total);
    s.download(I.h_status, 0, I.small_u32, 0, 4);
    const int m2 = s.mark("edge_scan");

    // ---- start-position node deposits (omega, initial charge, gross)
    const std::uint64_t begin = D(I.bounds), end = D(I.bounds, std::max(ne, nn + 1) * 4);
    auto node_sort = [&](const DeviceBuffer& ka, const DeviceBuffer& pa, const DeviceBuffer& kb, const DeviceBuffer& pb,
                         std::uint64_t& keys, std::uint64_t& pays) {
        const bool in_b = I.prim->radix_sort(D(ka), D(pa), D(kb), D(pb), 8 * n, I.node_bits, D(I.sort_scratch));
        keys = in_b ? D(kb) : D(ka);
        pays = in_b ? D(pb) : D(pa);
        s.zero(I.bounds, 0, 2 * std::max(ne, nn + 1) * 4);
        I.prim->segment_bounds(keys, 8 * n, begin, end);
    };
    double* node_out = I.nodes.as<double>();   // [8][nn]: omega init gross echarge hi lo fgross plasma
    {
        std::uint64_t keys = 0, pays = 0;
        node_sort(I.s_key_a, I.s_pay_a, I.s_key_b, I.s_pay_b, keys, pays);
        hgp::NodeReduceArgs r{hgp_u32(nn), 0, reinterpret_cast<const hgp_u32*>(pays), reinterpret_cast<const hgp_u32*>(begin),
                              reinterpret_cast<const hgp_u32*>(end), sv, sv + R, sv + 2 * R, node_out, node_out + nn, node_out + 2 * nn};
        void* args[] = {&r};
        s.launch(I.k_start, {blocks(nn, 128)}, {128}, 0, args);
    }
    const int m3 = s.mark("start_nodes");
    {
        std::uint64_t keys = 0, pays = 0;
        node_sort(I.e_key_a, I.e_pay_a, I.e_key_b, I.e_pay_b, keys, pays);
        hgp::EndNodeArgsV3 r{hgp_u32(nn), 0, reinterpret_cast<const hgp_u32*>(pays), reinterpret_cast<const hgp_u32*>(begin),
                             reinterpret_cast<const hgp_u32*>(end), ev, ev + R, ev + 2 * R, I.sinks_i.as<const int>(),
                             node_out + 3 * nn, node_out + 4 * nn, node_out + 5 * nn, node_out + 6 * nn};
        void* args[] = {&r};
        s.launch(I.k_end, {blocks(nn, 128)}, {128}, 0, args);
    }
    const int m4 = s.mark("end_nodes");

    // ---- exact maxima / first failures
    {
        const double* omega = node_out;
        const double* vn = I.VN.as<const double>();
        double* pd = node_out + 7 * nn;
        hgp_u32 nnu = hgp_u32(nn);
        void* args[] = {&omega, &vn, &nnu, &dt, &pd};
        s.launch(I.k_plasma, {blocks(nn, 256)}, {256}, 0, args);
        I.prim->reduce_max_f64(D(I.nodes, 7 * nn * 8), nn, D(I.red_scratch), D(I.small_f64, 0));
        I.prim->reduce_max_f64(D(I.outs, 13 * n * 8), n, D(I.red_scratch), D(I.small_f64, 8));
        I.prim->reduce_max_f64(D(I.outs, 14 * n * 8), n, D(I.red_scratch), D(I.small_f64, 16));
        for (int k = 0; k < 3; ++k) I.prim->reduce_min_u32(D(I.fails, std::size_t(k) * n * 4), n, D(I.red_scratch), min_fail + 4 * k);
    }
    const int m5 = s.mark("reductions");
    s.download(I.h_out, 0, I.outs, 0, 15 * n * 8);
    s.download(I.h_sinks_i, 0, I.sinks_i, 0, 2 * n * 4);
    s.download(I.h_sinks_r, 0, I.sinks_r, 0, 7 * n * 8);
    s.download(I.h_nodes, 0, I.nodes, nn * 8, 6 * nn * 8);
    const int m6 = s.mark("d2h_carriers_nodes");

    // ---- edge deposits: emit in carrier order, stable sort by edge, sequential segment sums
    s.wait_marker(m2, "edge record total");
    const std::uint32_t T = I.h_status.as<std::uint32_t>()[0];
    v.edge_records = T;
    if (T > I.edge_cap) I.grow_edges(T);
    {
        hgp::EmitArgs e{};
        e.g = g;
        e.dt = dt;
        e.n = hgp_u32(n);
        e.q = a.q; e.w = a.w; e.start = a.start; e.disp = a.disp;
        e.edge_count = I.counts.as<const hgp_u32>();
        e.edge_offset = I.offsets.as<const hgp_u32>();
        e.key = I.g_key_a.as<hgp_u32>(); e.pay = I.g_pay_a.as<hgp_u32>(); e.value = I.g_val.as<double>();
        void* args[] = {&e};
        s.launch(I.k_emit, {blocks(n, 128)}, {128}, 0, args);
    }
    const int m7 = s.mark("emit_edges");
    const bool in_b = I.prim->radix_sort(D(I.g_key_a), D(I.g_pay_a), D(I.g_key_b), D(I.g_pay_b), T, I.edge_bits, D(I.sort_scratch));
    const std::uint64_t ekeys = in_b ? D(I.g_key_b) : D(I.g_key_a), epays = in_b ? D(I.g_pay_b) : D(I.g_pay_a);
    const int m8 = s.mark("sort_edges");
    s.zero(I.bounds, 0, 2 * std::max(ne, nn + 1) * 4);
    I.prim->segment_bounds(ekeys, T, begin, end);
    {
        hgp::NodeReduceArgs r{hgp_u32(ne), 0, reinterpret_cast<const hgp_u32*>(epays), reinterpret_cast<const hgp_u32*>(begin),
                              reinterpret_cast<const hgp_u32*>(end), I.g_val.as<const double>(), nullptr, nullptr,
                              I.current.as<double>(), nullptr, nullptr};
        void* args[] = {&r};
        s.launch(I.k_edges, {blocks(ne, 128)}, {128}, 0, args);
    }
    const int m9 = s.mark("reduce_edges");
    {
        hgp::FinalizeArgs f{reinterpret_cast<const hgp_u32*>(min_fail), a.code1, a.code2, a.code3,
                            reinterpret_cast<const hgp_u32*>(edge_total), reinterpret_cast<hgp_u32*>(status)};
        void* args[] = {&f};
        s.launch(I.k_final, {1}, {32}, 0, args);
    }
    s.download(I.h_current, 0, I.current, 0, ne * 8);
    s.download(I.h_status, 0, I.small_u32, 64, 32);
    s.download(I.h_maxes, 0, I.small_f64, 0, 24);
    if (diagnostics) {
        s.download(I.h_codes, 0, I.codes, 0, 3 * n * 4);
        s.download(I.h_diag_real, 0, I.diag_real, 0, std::size_t(hgp::diag_real_fields) * n * 8);
        s.download(I.h_diag_int, 0, I.diag_int, 0, std::size_t(hgp::diag_int_fields) * n * 4);
    }
    const int m10 = s.mark("d2h_edges_status");
    s.synchronize("motion run");

    // ---- timings
    if (I.mark_carriers0 >= 0) v.gpu_phases.push_back({"h2d_carriers", s.elapsed_ms(I.mark_carriers0, I.mark_carriers1)});
    if (I.mark_field0 >= 0) v.gpu_phases.push_back({"h2d_field", s.elapsed_ms(I.mark_field0, I.mark_field1)});
    const int marks[] = {m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10};
    for (int k = 1; k < 11; ++k) v.gpu_phases.push_back({s.marker_label(marks[k]), s.elapsed_ms(marks[k - 1], marks[k])});
    v.gpu_phases.push_back({"gpu_total_run", s.elapsed_ms(m0, m10)});
    s.clear_markers();
    I.mark_carriers0 = I.mark_carriers1 = I.mark_field0 = I.mark_field1 = -1;

    // ---- views
    const double* ho = I.h_out.as<const double>();
    v.end_x = {ho, n}; v.end_y = {ho + n, n}; v.end_z = {ho + 2 * n, n};
    v.end_ux = {ho + 3 * n, n}; v.end_uy = {ho + 4 * n, n}; v.end_uz = {ho + 5 * n, n};
    v.start_x = {ho + 6 * n, n}; v.start_y = {ho + 7 * n, n}; v.start_z = {ho + 8 * n, n};
    v.disp_x = {ho + 9 * n, n}; v.disp_y = {ho + 10 * n, n}; v.disp_z = {ho + 11 * n, n};
    v.dke = {ho + 12 * n, n};
    {
        const int* si = I.h_sinks_i.as<const int>();
        const double* sr = I.h_sinks_r.as<const double>();
        v.sink_kind = {si, n}; v.sink_axis = {si + n, n};
        v.gyro = {ho + 13 * n, n}; v.dcells = {ho + 14 * n, n};
        v.sink_flight = {sr, n}; v.sink_remaining = {sr + n, n}; v.sink_fraction = {sr + 2 * n, n}; v.sink_plane = {sr + 3 * n, n};
        v.sink_segment_x = {sr + 4 * n, n}; v.sink_segment_y = {sr + 5 * n, n}; v.sink_segment_z = {sr + 6 * n, n};
    }
    const double* hn = I.h_nodes.as<const double>();
    v.initial_charge_C = {hn, nn}; v.gross_C = {hn + nn, nn}; v.endpoint_charge_C = {hn + 2 * nn, nn};
    v.free_net_hi = {hn + 3 * nn, nn}; v.free_net_lo = {hn + 4 * nn, nn}; v.free_gross_C = {hn + 5 * nn, nn};
    v.current_A_m = {I.h_current.as<const double>(), ne};
    if (diagnostics) {
        const int* hc = I.h_codes.as<const int>();
        v.code1 = {hc, n}; v.code2 = {hc + n, n}; v.code3 = {hc + 2 * n, n};
        v.diag_real = {I.h_diag_real.as<const double>(), std::size_t(hgp::diag_real_fields) * n};
        v.diag_int = {I.h_diag_int.as<const int>(), std::size_t(hgp::diag_int_fields) * n};
    }
    const std::uint32_t* st = I.h_status.as<const std::uint32_t>();
    const double* mx = I.h_maxes.as<const double>();
    v.budget.plasma_dt = mx[0];
    v.budget.maximum_gyro = mx[1];
    v.budget.maximum_displacement_cells = mx[2];

    // ---- refusal order of prepare_canonical_motion, then the host folds
    const auto fold0 = Clock::now();
    FoldInput f;
    f.n = n;
    f.ids = I.ids;
    f.weights = I.weights;
    f.dke = v.dke;
    f.k0 = I.k0;
    f.charge = I.charge;
    f.forms_v7 = I.forms_v7;
    f.E = {I.h_field.as<const double>(), ne};
    f.current = v.current_A_m;
    f.initial = v.initial_charge_C;
    f.endpoint = v.endpoint_charge_C;
    f.gross = v.gross_C;
    f.net_hi = v.free_net_hi;
    f.net_lo = v.free_net_lo;
    f.free_gross = v.free_gross_C;
    f.mt = &I.mt;
    f.dt = dt;
    for (int k = 0; k < 6; ++k) f.status[k] = st[k];
    fold_record(f, v.refusal, v.charge_refusal, v.budget, v.work);
    for (auto* rf : {&v.refusal, &v.charge_refusal})
        if (rf->refused && rf->code >= hgp::sink_field_without_kinetic && rf->code <= hgp::sink_unsupported_bank) rf->message = code_message_v3(rf->code);
    v.host_fold_ms = ms_since(fold0);
    v.run_wall_ms = ms_since(wall0);
    return v;
}

}  // namespace hall_gpu::push_v3
