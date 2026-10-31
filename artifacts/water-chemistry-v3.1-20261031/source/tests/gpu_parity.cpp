// gpu_parity.cpp -- CPU/GPU bit identity of chemistry v3.
// usage: gpu_parity BLOB OUTPUT_DIR NVRTC_LIBRARY DEVICE_UUID CARRIERS REPEATS
// The host runs the CPU reference (chem_event.h compiled by the host compiler);
// the device runs device/chem_kernels.cu compiled by NVRTC from the same header
// bytes. Compared bit for bit: every CarrierOut and EventRec, the canonical
// record digest and the exact ledger digest, for several step configurations,
// three block sizes and a split launch; plus elementary-function parity.
// exit: 0 identical, 4 mismatch, 2 repeat mismatch, 3 teardown quarantine, 1 error.
#include "common/hall_gpu_common.hpp"
#include "test_common.hpp"

#include <chrono>
#include <cstring>
#include <filesystem>
#include <iostream>

using namespace pfc;
using namespace pfc_test;

namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kKernel[] = {
#embed "../device/chem_kernels.cu"
    , 0};
const char kDet[] = {
#embed "../include/pf_chem_v3/det_math.h"
    , 0};
const char kRng[] = {
#embed "../include/pf_chem_v3/chem_rng.h"
    , 0};
const char kKin[] = {
#embed "../include/pf_chem_v3/chem_kinematics.h"
    , 0};
const char kTab[] = {
#embed "../include/pf_chem_v3/chem_tables.h"
    , 0};
const char kEvt[] = {
#embed "../include/pf_chem_v3/chem_event.h"
    , 0};
const char kPhilox[] = {
#embed "/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival/philox_core.h"
    , 0};
#pragma GCC diagnostic pop

template <class T>
std::size_t mismatches(const T* a, const T* b, std::size_t n) {
    std::size_t m = 0;
    for (std::size_t i = 0; i < n; ++i) m += std::memcmp(a + i, b + i, sizeof(T)) != 0;
    return m;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "usage: gpu_parity BLOB OUTPUT_DIR NVRTC_LIBRARY DEVICE_UUID CARRIERS REPEATS\n";
        return 1;
    }
    const std::string out = argv[2];
    std::filesystem::create_directories(out);
    const std::size_t n = std::stoull(argv[5]);
    const int repeats = std::stoi(argv[6]);
    int rc = 0;
    std::ostringstream js;
    hall_gpu::TeardownReport teardown;
    try {
        const Blob b = load_blob(argv[1]);
        hall_gpu::SessionOptions so;
        so.require_uuid = argv[4];
        so.nvrtc_library = argv[3];
        hall_gpu::Session session(so);
        hall_gpu::CompileRequest req;
        req.main = {"chem_kernels.cu", kKernel};
        req.headers = {{"det_math.h", kDet}, {"chem_rng.h", kRng}, {"chem_kinematics.h", kKin}, {"chem_tables.h", kTab},
                       {"chem_event.h", kEvt}, {"philox_core.h", kPhilox}};
        const auto module = session.compile_and_load(req);
        const auto k_step = session.kernel(module, "pfc_mcc_step");
        const auto k_math = session.kernel(module, "pfc_math_parity");
        // table blob on the device
        auto d_blob = session.alloc(b.bytes, "table blob");
        auto h_blob = session.alloc_pinned(b.bytes, "table staging");
        std::memcpy(h_blob.as<char>(), b.base(), b.bytes);
        session.upload(d_blob, 0, h_blob, 0, b.bytes);
        TableView tv = b.view_at(d_blob.address);
        CoverView cv = b.cover_at(d_blob.address);
        js << "{\"device\": " << hall_gpu::device_json(session.device()) << ", \"compile\": " << hall_gpu::compile_json(module.record)
           << ", \"blob_sha256\": \"" << b.sha256 << "\", \"carriers\": " << n << ", \"cases\": [";
        // ------------------------------------------------ math parity
        {
            const std::size_t m = 1 << 20;
            std::vector<double> x(m), y(m), u(m), host(4 * m), dev(4 * m);
            std::uint64_t s = 31337;
            for (std::size_t i = 0; i < m; ++i) {
                x[i] = (i % 3 == 0) ? std::ldexp(1.0 + u01(s), int(splitmix(s) % 2000) - 1000) : (i % 3 == 1 ? u01(s) : 1.0 + (u01(s) - 0.5) * 1e-7);
                y[i] = (u01(s) - 0.5) * ((i & 1) ? 1.9 : 1e-9);
                u[i] = u01(s);
                host[4 * i] = pfc::log(x[i]);
                host[4 * i + 1] = pfc::log1p(y[i]);
                pfc::sincos_2pi(u[i], &host[4 * i + 2], &host[4 * i + 3]);
            }
            auto d_in = session.alloc(3 * m * 8, "math in");
            auto d_out = session.alloc(4 * m * 8, "math out");
            auto h_in = session.alloc_pinned(3 * m * 8, "math in staging");
            auto h_out = session.alloc_pinned(4 * m * 8, "math out staging");
            std::memcpy(h_in.as<char>(), x.data(), m * 8);
            std::memcpy(h_in.as<char>() + m * 8, y.data(), m * 8);
            std::memcpy(h_in.as<char>() + 2 * m * 8, u.data(), m * 8);
            session.upload(d_in, 0, h_in, 0, 3 * m * 8);
            std::uint64_t ax = d_in.address, ay = ax + m * 8, au = ay + m * 8, ao = d_out.address;
            unsigned cnt = unsigned(m);
            void* args[] = {&ax, &ay, &au, &cnt, &ao};
            session.launch(k_math, {unsigned((m + 255) / 256)}, {256}, 0, args);
            session.download(h_out, 0, d_out, 0, 4 * m * 8);
            session.synchronize("math parity");
            std::memcpy(dev.data(), h_out.as<char>(), 4 * m * 8);
            const std::size_t mm = mismatches(host.data(), dev.data(), 4 * m);
            js << "{\"case\": \"math\", \"values\": " << 4 * m << ", \"mismatches\": " << mm << "}";
            std::cout << (mm ? "FAIL" : "PASS") << " math parity " << 4 * m << " values, " << mm << " mismatches\n";
            if (mm) rc = 4;
            session.release(d_in);
            session.release(d_out);
            session.release(h_in);
            session.release(h_out);
        }
        // ------------------------------------------------ step parity
        struct Case {
            std::uint64_t seed, step, pop_seed;
            double dt;
            std::size_t cells;
        };
        const Case cases[] = {{0x5EED5EED12345678ull, 4242, 1000, 4e-11, 97}, {0xABCDEF0123456789ull, 77, 4711, 5e-12, 31}, {3, (1ull << 39) + 5, 99, 2e-10, 7}};
        for (const auto& cs : cases) {
            StepConfig cfg{};
            cfg.seed = cs.seed;
            cfg.step = cs.step;
            cfg.dt_s = cs.dt;
            cfg.trunc_sigmas = 7.0;
            cfg.purpose = PFC_PURPOSE_MCC;
            cfg.max_candidates = 4096;
            cfg.max_events = 64;
            const Population P = make_population(b, n, cs.cells, cs.pop_seed);
            const auto t0 = std::chrono::steady_clock::now();
            StepResult ref;
            cpu_step(b, cfg, P.carriers, P.cells, ref);
            const double cpu_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            const std::string dref = step_digest(ref.outs, ref.events, cfg.max_events, P.carriers);
            const std::string lref = ledger_digest(ledger_of(b, ref, cfg.max_events));
            std::uint64_t n_events = 0;
            for (auto& o : ref.outs) n_events += o.n_events;
            const std::size_t in_bytes = n * sizeof(CarrierIn) + P.cells.size() * sizeof(CellTargets);
            const std::size_t ev_bytes = n * cfg.max_events * sizeof(EventRec), out_bytes = n * sizeof(CarrierOut);
            auto d_in = session.alloc(in_bytes, "carriers+cells");
            auto d_ev = session.alloc(ev_bytes, "events");
            auto d_out = session.alloc(out_bytes, "carrier outs");
            auto h_in = session.alloc_pinned(in_bytes, "in staging");
            auto h_ev = session.alloc_pinned(ev_bytes, "events staging");
            auto h_out = session.alloc_pinned(out_bytes, "outs staging");
            std::memcpy(h_in.as<char>(), P.carriers.data(), n * sizeof(CarrierIn));
            std::memcpy(h_in.as<char>() + n * sizeof(CarrierIn), P.cells.data(), P.cells.size() * sizeof(CellTargets));
            session.upload(d_in, 0, h_in, 0, in_bytes);
            std::uint64_t a_car = d_in.address, a_cell = a_car + n * sizeof(CarrierIn), a_ev = d_ev.address, a_out = d_out.address;
            unsigned n_cells = unsigned(P.cells.size()), nn = unsigned(n);
            struct Launch {
                unsigned block;
                bool split;
            };
            const Launch launches[] = {{128, false}, {64, false}, {256, true}};
            for (const auto& L : launches) {
                for (int r = 0; r < repeats; ++r) {
                    session.zero(d_ev, 0, ev_bytes);
                    session.zero(d_out, 0, out_bytes);
                    const int m0 = session.mark("step-begin");
                    if (!L.split) {
                        unsigned first = 0;
                        void* args[] = {&tv, &cv, &cfg, &a_car, &a_cell, &n_cells, &nn, &first, &a_ev, &a_out};
                        session.launch(k_step, {(nn + L.block - 1) / L.block}, {L.block}, 0, args);
                    } else {
                        // two launches over [n/2, n) then [0, n/2): order independence on the device
                        unsigned half = nn / 2, first = half, upto = nn;
                        void* a1[] = {&tv, &cv, &cfg, &a_car, &a_cell, &n_cells, &upto, &first, &a_ev, &a_out};
                        session.launch(k_step, {(nn - half + L.block - 1) / L.block}, {L.block}, 0, a1);
                        unsigned first0 = 0, upto0 = half;
                        void* a2[] = {&tv, &cv, &cfg, &a_car, &a_cell, &n_cells, &upto0, &first0, &a_ev, &a_out};
                        session.launch(k_step, {(half + L.block - 1) / L.block}, {L.block}, 0, a2);
                    }
                    const int m1 = session.mark("step-end");
                    session.download(h_ev, 0, d_ev, 0, ev_bytes);
                    session.download(h_out, 0, d_out, 0, out_bytes);
                    session.synchronize("step parity");
                    const double gpu_ms = session.elapsed_ms(m0, m1);
                    session.clear_markers();
                    StepResult dev;
                    dev.outs.resize(n);
                    dev.events.resize(n * cfg.max_events);
                    std::memcpy(dev.outs.data(), h_out.as<char>(), out_bytes);
                    std::memcpy(dev.events.data(), h_ev.as<char>(), ev_bytes);
                    // compare the populated parts (unused event slots are zero on the device and value-initialized on the host)
                    std::size_t mm_out = mismatches(ref.outs.data(), dev.outs.data(), n), mm_ev = 0;
                    for (std::size_t i = 0; i < n; ++i)
                        mm_ev += mismatches(ref.events.data() + i * cfg.max_events, dev.events.data() + i * cfg.max_events, ref.outs[i].n_events);
                    const std::string d = step_digest(dev.outs, dev.events, cfg.max_events, P.carriers);
                    const std::string l = ledger_digest(ledger_of(b, dev, cfg.max_events));
                    const bool same = mm_out == 0 && mm_ev == 0 && d == dref && l == lref;
                    if (!same) rc = 4;
                    js << ", {\"case\": \"step\", \"seed\": " << cs.seed << ", \"step\": " << cs.step << ", \"dt_s\": " << cs.dt << ", \"block\": " << L.block
                       << ", \"split_launch\": " << (L.split ? "true" : "false") << ", \"repeat\": " << r << ", \"events\": " << n_events
                       << ", \"carrier_mismatches\": " << mm_out << ", \"event_mismatches\": " << mm_ev << ", \"record_digest_cpu\": \"" << dref
                       << "\", \"record_digest_gpu\": \"" << d << "\", \"ledger_digest_cpu\": \"" << lref << "\", \"ledger_digest_gpu\": \"" << l
                       << "\", \"cpu_ms\": " << cpu_ms << ", \"gpu_kernel_ms\": " << gpu_ms << "}";
                    std::cout << (same ? "PASS" : "FAIL") << " step seed " << cs.seed << " step " << cs.step << " block " << L.block << (L.split ? " split" : "")
                              << " repeat " << r << ": " << n_events << " events, mismatches " << mm_out << "/" << mm_ev << ", cpu " << cpu_ms
                              << " ms, gpu kernel " << gpu_ms << " ms\n";
                }
            }
            session.release(d_in);
            session.release(d_ev);
            session.release(d_out);
            session.release(h_in);
            session.release(h_ev);
            session.release(h_out);
        }
        js << "], \"memory\": " << hall_gpu::memory_json(session.memory());
        teardown = session.close();
        js << ", \"teardown_quarantined\": " << (teardown.quarantined ? "true" : "false") << ", \"exit\": " << rc << "}";
        if (teardown.quarantined) rc = 3;
    } catch (const std::exception& e) {
        std::cerr << "ERROR " << e.what() << "\n";
        js << "], \"error\": \"" << e.what() << "\"}";
        rc = 1;
    }
    std::ofstream(out + "/gpu-parity.json") << js.str() << "\n";
    std::cout << (rc == 0 ? "ALL IDENTICAL" : "NOT IDENTICAL") << " exit " << rc << "\n";
    return rc;
}
