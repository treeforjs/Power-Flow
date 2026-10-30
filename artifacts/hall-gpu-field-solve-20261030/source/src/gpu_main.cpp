// GPU side of the Hall field backend qualification.
//
//   hall_field_gpu FIXROOT OUTDIR WARMUP REPEATS BATCH GRAPH BLOCKING SWEEP
//
// FIXROOT/prep.json lists the fixtures written by `hall_field_cpu prep`; every
// fixture file is checked against the sha256 recorded there. For each fixture
// and each operator representation the program checks device products against
// the host references bit for bit, solves repeatedly (device-resident; with
// right-hand-side/solution transfers; with operator-value refresh as well),
// checks repeat-solve bit identity and equality with the host emulation of the
// same algorithm, prechecks the residual on the CPU, and writes the solution
// for `hall_field_cpu guard`.
#include "gpu_cg.hpp"
#include "raikou/json.h"
#include <chrono>
#include <filesystem>
#include <iostream>

using Json = raikou::Json;
using Clock = std::chrono::steady_clock;
using hgf::Vec;

namespace {
double ms_since(Clock::time_point t0) { return std::chrono::duration<double, std::milli>(Clock::now() - t0).count(); }
double median(std::vector<double> v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    return v.size() % 2 ? v[v.size() / 2] : .5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
}
Json timing(const std::vector<double>& v) {
    return {{"samples_ms", v}, {"median_ms", median(v)},
            {"min_ms", v.empty() ? 0 : *std::min_element(v.begin(), v.end())},
            {"max_ms", v.empty() ? 0 : *std::max_element(v.begin(), v.end())}};
}
Json diff_json(const hgf::Diff& d) {
    return {{"bitwise_equal", d.bitwise}, {"differing_components", d.differing}, {"max_abs", d.max_abs},
            {"max_rel", d.max_rel}, {"normwise_rel", d.norm_rel}};
}
// Same arithmetic as foil::dot / foil::norm (foil_em.hpp:10-11).
double canonical_norm(const Vec& a) {
    long double s = 0;
    for (double x : a) s += (long double)x * x;
    return std::sqrt(double(s));
}
}  // namespace

int main(int argc, char** argv) {
    int rc = 0;
    Json report{{"argv", std::vector<std::string>(argv, argv + argc)}};
    hgf::cu::CleanupReport cleanup;
    std::filesystem::path out;
    try {
        hgf::need(argc == 9, "usage: hall_field_gpu FIXROOT OUTDIR WARMUP REPEATS BATCH GRAPH BLOCKING SWEEP");
        const std::filesystem::path root = argv[1];
        out = argv[2];
        const int warmup = std::stoi(argv[3]), repeats = std::stoi(argv[4]), batch = std::stoi(argv[5]);
        const bool graph = std::stoi(argv[6]) != 0, blocking = std::stoi(argv[7]) != 0, sweep = std::stoi(argv[8]) != 0;
        hgf::need(batch >= 1 && repeats >= 1 && warmup >= 0, "batch/repeats/warmup");
        std::filesystem::create_directories(out);
        Json prep;
        { std::ifstream in(root / "prep.json"); hgf::need(bool(in), "prep.json missing"); in >> prep; }
        report["prep_json_sha256"] = hgf::sha256_file((root / "prep.json").string());
        report["settings"] = {{"warmup", warmup}, {"repeats", repeats}, {"batch", batch}, {"graph", graph}, {"blocking_sync", blocking},
                              {"reduction_blocks", hgf::NBLK}, {"threads_per_block", hgf::TPB}, {"transfers", "pinned host buffers, session stream"}};
        {
            hgf::cu::Session session(0, blocking, cleanup);
            const auto& info = session.info();
            report["device"] = {{"ordinal", info.ordinal}, {"name", info.name}, {"uuid", info.uuid}, {"pci_bus_id", info.pci_bus_id},
                                {"compute_capability", std::to_string(info.major) + "." + std::to_string(info.minor)},
                                {"sm_count", info.sm_count}, {"total_memory_bytes", info.total_memory},
                                {"driver_api_version", info.driver_api_version}};
            auto t0 = Clock::now();
            const auto compiled = hgf::cu::compile(hgf::kernel_source, "hall_field_cg.cu", info.major, info.minor);
            report["nvrtc"] = {{"version", std::to_string(compiled.nvrtc_major) + "." + std::to_string(compiled.nvrtc_minor)},
                               {"options", compiled.options}, {"supported_archs", compiled.supported_archs}, {"arch_decision", compiled.arch_decision},
                               {"image_kind", compiled.ptx ? "PTX" : "CUBIN"}, {"image_sha256", compiled.image_sha256},
                               {"source_sha256", compiled.source_sha256}, {"log", compiled.log}, {"compile_ms", ms_since(t0)}};
            session.load(compiled);
            Json fixtures = Json::object();
            for (auto& [fname, rec] : prep["fixtures"].items()) {
                const std::string path = rec["fixture_file"];
                const std::string expect = rec["fixture_sha256"];
                const auto actual = hgf::sha256_file(path);
                hgf::need(actual == expect, "fixture sha256 differs from prep record: " + path);
                const auto F = hgf::ArrayFile::read(path);
                const auto A = hgf::load_operator(F);
                const Vec b = F.f64("b"), x_cpu = F.f64("x_cpu");
                const int maxit = F.i32("maxit").at(0);
                const double target_canonical = F.scalar("target_canonical");
                Json fr{{"fixture_file", path}, {"fixture_sha256", actual}, {"rows", A.n}, {"nnz_assembled", A.nnz()}, {"maxit", maxit}};
                auto ts = Clock::now();
                hgf::GpuSolver solver(session, A, maxit);
                fr["static_upload_ms"] = ms_since(ts);
                fr["operator_value_bytes"] = solver.operator_value_bytes();
                fr["sigma_rows"] = solver.sigma_rows();
                fr["sigma_refresh_bytes"] = solver.sigma_refresh_bytes();
                // Device products versus host references (bitwise expected).
                Json products = Json::array();
                Vec Bh, y(b.size());
                for (auto [label, v] : std::vector<std::pair<std::string, Vec>>{{"rhs", b}, {"x_cpu", x_cpu}, {"pseudo_random_seed_3", hgf::pseudo_random(b.size(), 3)}}) {
                    hgf::spmv_csr(A, v.data(), y.data());
                    const auto da = hgf::compare(solver.product(hgf::OpMode::assembled, v), y);
                    hgf::apply_factored(A, v.data(), y.data(), Bh);
                    const auto df = hgf::compare(solver.product(hgf::OpMode::factored, v), y);
                    products.push_back({{"vector", label}, {"device_assembled_vs_host_csr", diff_json(da)}, {"device_factored_vs_host_factored", diff_json(df)}});
                }
                fr["device_products"] = products;
                solver.stage_rhs(b);
                solver.upload_rhs();
                session.sync();
                Json modes = Json::object();
                for (auto m : {hgf::OpMode::assembled, hgf::OpMode::factored}) {
                    Json mr;
                    const Vec x_emu = F.f64(std::string("x_emu_") + hgf::name(m));
                    std::vector<double> resident, with_io, with_sigma, with_refresh, event_ms;
                    std::vector<std::string> digests;
                    hgf::SolveStats st;
                    Vec x;
                    auto ev0 = session.event(), ev1 = session.event();
                    auto keep = [&](const Vec& xi) {
                        digests.push_back(hgf::sha256_vec(xi));
                        if (x.empty()) x = xi;
                    };
                    for (int i = 0; i < warmup + repeats; ++i) {  // device-resident: b and x stay on the device
                        auto t = Clock::now();
                        session.record(ev0);
                        st = solver.solve(m, batch, graph);
                        session.record(ev1);
                        session.sync();
                        const double w = ms_since(t);
                        solver.download_x(); session.sync();
                        keep(solver.x_host());
                        if (i >= warmup) { resident.push_back(w); event_ms.push_back(session.elapsed_ms(ev0, ev1)); }
                    }
                    for (int i = 0; i < repeats; ++i) {  // + H2D of b and D2H of x
                        auto t = Clock::now();
                        solver.upload_rhs();
                        st = solver.solve(m, batch, graph);
                        solver.download_x();
                        session.sync();
                        with_io.push_back(ms_since(t));
                        keep(solver.x_host());
                    }
                    for (int i = 0; i < repeats; ++i) {  // + sigma(T) refresh of the conduction rows (per-step cost)
                        auto t = Clock::now();
                        solver.refresh_sigma_rows();
                        solver.upload_rhs();
                        st = solver.solve(m, batch, graph);
                        solver.download_x();
                        session.sync();
                        with_sigma.push_back(ms_since(t));
                        keep(solver.x_host());
                    }
                    for (int i = 0; i < repeats; ++i) {  // + upload of every operator value (upper bound, e.g. dt change)
                        auto t = Clock::now();
                        solver.upload_operator_values();
                        solver.upload_rhs();
                        st = solver.solve(m, batch, graph);
                        solver.download_x();
                        session.sync();
                        with_refresh.push_back(ms_since(t));
                        keep(solver.x_host());
                    }
                    bool identical = true;
                    for (const auto& d : digests) identical = identical && d == digests.front();
                    std::vector<std::vector<double>> ha, hb, hr;
                    solver.histories(st, ha, hb, hr);
                    const auto ritz = hgf::lanczos_extremes(ha.at(0), hb.at(0));
                    // CPU precheck with the host factored apply (bitwise equal to the
                    // sealed matrix-free apply per the prep records); the binding
                    // check is `hall_field_cpu guard`.
                    hgf::apply_factored(A, x.data(), y.data(), Bh);
                    Vec r(b.size());
                    for (std::size_t i = 0; i < b.size(); ++i) r[i] = b[i] - y[i];
                    const double res = canonical_norm(r);
                    int total = 0;
                    for (int k : st.iterations) total += k;
                    mr["solves"] = int(digests.size());
                    mr["repeat_solves_bitwise_identical"] = identical;
                    mr["x_sha256"] = digests.front();
                    mr["x_vs_host_emulation"] = diff_json(hgf::compare(x, x_emu));
                    mr["x_vs_cpu_sealed"] = diff_json(hgf::compare(x, x_cpu));
                    mr["iterations_per_correction"] = st.iterations;
                    mr["iterations_total"] = total;
                    mr["corrections"] = st.corrections;
                    mr["control_readbacks_last_solve"] = st.readbacks;
                    mr["device_bnorm"] = st.bnorm;
                    mr["device_target"] = st.target;
                    mr["device_true_residual"] = st.true_norm;
                    mr["residual_history_recursive"] = hr;
                    mr["lanczos"] = {{"steps", ritz.k}, {"lambda_min", ritz.lmin}, {"lambda_max", ritz.lmax}, {"kappa_estimate", ritz.kappa}};
                    mr["cpu_precheck"] = {{"residual", res}, {"target_canonical", target_canonical}, {"passed", res <= target_canonical},
                                          {"operator", "host factored apply"}, {"norm", "long double accumulation, canonical order"}};
                    mr["timing_resident"] = timing(resident);
                    mr["timing_resident_events"] = timing(event_ms);
                    mr["timing_with_rhs_solution_transfers"] = timing(with_io);
                    mr["timing_with_sigma_refresh_and_transfers"] = timing(with_sigma);
                    mr["timing_with_full_operator_upload_and_transfers"] = timing(with_refresh);
                    if (sweep) {
                        Json sw = Json::array();
                        for (int bsz : {1, 2, 4, 8, 16})
                            for (bool g : {false, true}) {
                                std::vector<double> ms;
                                bool same = true;
                                for (int i = 0; i < 1 + std::max(3, repeats / 2); ++i) {
                                    auto t = Clock::now();
                                    auto s2 = solver.solve(m, bsz, g);
                                    session.sync();
                                    const double w = ms_since(t);
                                    solver.download_x(); session.sync();
                                    same = same && hgf::sha256_vec(solver.x_host()) == digests.front();
                                    if (i) ms.push_back(w);
                                }
                                sw.push_back({{"batch", bsz}, {"graph", g}, {"median_ms", median(ms)}, {"min_ms", *std::min_element(ms.begin(), ms.end())},
                                              {"bitwise_identical_to_default", same}});
                            }
                        mr["batch_sweep_resident"] = sw;
                    }
                    const auto dir = out / fname;
                    std::filesystem::create_directories(dir);
                    hgf::ArrayFile X;
                    X.put("x", x);
                    const auto xp = (dir / (std::string("x_gpu_") + hgf::name(m) + ".bin")).string();
                    X.write(xp);
                    mr["x_file"] = xp;
                    mr["x_file_sha256"] = hgf::sha256_file(xp);
                    modes[hgf::name(m)] = mr;
                    std::cout << fname << " " << hgf::name(m) << " iterations=" << total << " resident_median_ms=" << median(resident)
                              << " identical=" << identical << " emulation_bitwise=" << hgf::compare(x, x_emu).bitwise << std::endl;
                }
                fr["modes"] = modes;
                fixtures[fname] = fr;
            }
            report["fixtures"] = fixtures;
            session.close();
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        report["error"] = e.what();
        rc = 1;
    }
    Json issues = Json::array();
    for (const auto& i : cleanup.issues) issues.push_back({{"operation", i.operation}, {"code", i.code}});
    report["cleanup"] = {{"closed", cleanup.closed}, {"drained", cleanup.drained}, {"quarantined", cleanup.quarantined},
                         {"primary_released", cleanup.primary_released}, {"context_restored", cleanup.context_restored},
                         {"freed_device_allocations", cleanup.freed_device}, {"freed_pinned_allocations", cleanup.freed_pinned}, {"issues", issues}};
    if (cleanup.quarantined || !cleanup.issues.empty()) rc = rc ? rc : 2;
    if (!out.empty()) std::ofstream(out / "gpu.json") << report.dump(2) << '\n';
    return rc;
}
