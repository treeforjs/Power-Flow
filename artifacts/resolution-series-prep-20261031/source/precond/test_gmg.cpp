// Deterministic GMG-preconditioned CG for the Hall implicit edge operator: measurement against the Hall gpu1 Jacobi PCG
// algorithm on the same right-hand sides, with the Hall v3.1 independent residual check (HostApply = the sealed apply
// arithmetic, foil::norm = the sealed long-double norm, target max(1e-30, 2e-13 ||b||)).
//
//   gmg_solve OUT TRACE NVRTC UUID N AXIAL AIR RHS CPUREF DT1 [DT2 ...]
//     RHS     random | smooth | physical (physical: the Hall vacuum-step right-hand side from the DC seed; needs the
//             mesh to be admitted by the Hall Model)
//     NVRTC   path of libnvrtc.so.12, or none: CPU reference only (sanitizer builds; no GPU, no lock)
//     CPUREF  1[,NU[,NU_COARSE]]: also run the CPU reference of the GMG solve and compare the increment bit for bit with the
//             GPU; NU smoothing sweeps per level (default 1), NU_COARSE coarsest-level sweeps (default 4)
// GPU runs execute under the GPU lock wrappers (GPU locks first, then the core lease); build_run.py checks both are held.
#include "resident_owner.hpp"
#include "hall-finite-foil-em-20261029/axial_particles.hpp"
namespace hc = hall_corrective_v3;
#include "motion_fixtures.hpp"
#include "gmg_gpu.hpp"
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <random>
#include <sys/file.h>
#include <unistd.h>

using J = hc::Json;
using steady = std::chrono::steady_clock;
namespace {
double ms_since(steady::time_point t) { return std::chrono::duration<double, std::milli>(steady::now() - t).count(); }
long rss_kb() { std::ifstream f("/proc/self/status"); std::string l; while (std::getline(f, l)) if (l.rfind("VmHWM:", 0) == 0) return std::stol(l.substr(6)); return -1; }
axial::Mesh make_mesh(int nT, int nA, int nAir) {
    gmg::need(nT >= 9 && nT <= 500 && nA >= 2 && nAir >= 1, "invalid immutable mesh specification");
    return axial::Mesh(nT, nA, .009, nAir);
}
J level_json(const gmg::Hierarchy& H) {
    J a = J::array();
    for (std::size_t i = 0; i < H.L.size(); ++i) {
        const auto& l = H.L[i];
        std::size_t fixed = 0; for (auto f : l.fixed) fixed += f;
        a.push_back({{"level", i}, {"n", l.g.n}, {"nz", l.g.nz}, {"edges", l.g.ne}, {"nodes", l.g.nn}, {"fixed_edges", fixed}, {"K_nnz", l.K.nnz()},
                     {"N_nnz", l.N.nnz()}, {"P_nnz", l.P.nnz()}, {"G_nnz", l.G.nnz()}, {"ratio_to_next", i < H.ratio.size() ? H.ratio[i] : 0}});
    }
    return a;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 11) { std::cerr << "usage: see the header of test_gmg.cpp\n"; return 64; }
        const std::filesystem::path out = argv[1]; const std::string trace = argv[2], nvrtc = argv[3], uuid = argv[4];
        const int nT = std::stoi(argv[5]), nA = std::stoi(argv[6]), nAir = std::stoi(argv[7]);
        const std::string rhs_kind = argv[8];
        // CPUREF[,NU[,NU_COARSE]]
        int cpuref_i = 0, nu = 1, nuc = 4; std::sscanf(argv[9], "%d,%d,%d", &cpuref_i, &nu, &nuc); const bool cpuref = cpuref_i == 1;
        std::vector<double> dts; for (int i = 10; i < argc; ++i) dts.push_back(std::stod(argv[i]));
        std::filesystem::create_directories(out);
        J summary{{"schema", "resolution-series-gmg-v1"}, {"spec", {nT, nA, nAir, .009}}, {"rhs", rhs_kind}, {"cpu_reference", cpuref}};
        auto g = make_mesh(nT, nA, nAir);
        auto law = electrode_heat_v2::conditional_ss304_v2(ss304v2::Qualification::evaluated_solid_1500,
            ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation, 8000., electrode_heat_v2::material_headers_sha256);
        gmg::Vec sigma(g.cells.size()); for (std::size_t c = 0; c < sigma.size(); ++c) if (g.cells[c].body >= 0) sigma[c] = law.sigma(300.);
        const auto S = g.conduction(sigma);
        std::unique_ptr<hc::Model> model;
        if (rhs_kind == "physical") model = std::make_unique<hc::Model>(hc::MeshSpec{nT, nA, nAir, .009}, trace);
        auto t = steady::now();
        gmg::Options opt; opt.nu = nu; opt.nu_coarse = nuc; auto H = gmg::build(g, S, opt);
        summary["options"] = {{"nu", nu}, {"nu_coarse", nuc}, {"min_cells", opt.min_cells}};
        summary["hierarchy_build_ms"] = ms_since(t);
        gmg::set_dt(H, dts.front());
        summary["levels"] = level_json(H);
        summary["host_rss_after_build_kB"] = rss_kb();
        hc::gpu::HostApply host_apply(g);
        // right-hand sides (per dt)
        auto make_rhs = [&](double dt) {
            gmg::Vec b(g.ne, 0.);
            if (rhs_kind == "random") {
                std::mt19937_64 rng(20261031); std::uniform_real_distribution<double> U(-1, 1);
                for (int e = 0; e < g.ne; ++e) { const double u = U(rng); if (!g.fixed_E[e]) b[e] = u; }
            } else if (rhs_kind == "smooth") {
                gmg::Vec xs(g.ne);
                for (int e = 0; e < g.ne; ++e) {
                    const auto pa = g.node_position(g.edges[e].a), pb = g.node_position(g.edges[e].b);
                    const double x = .5 * (pa[0] + pb[0]) / g.side, y = .5 * (pa[1] + pb[1]) / g.side, z = .5 * (pa[2] + pb[2]) / (g.nz * g.dz);
                    xs[e] = g.fixed_E[e] ? 0. : std::sqrt(g.M[e]) * (std::sin(3.14159265358979 * (x + .5)) * std::cos(6.28318530717959 * y) + .3 * std::sin(6.28318530717959 * z));
                }
                host_apply.apply(S, dt, xs, b);
                for (int e = 0; e < g.ne; ++e) if (g.fixed_E[e]) b[e] = 0;
            } else {
                const auto& in = model->initial_fields(); const auto seed = model->dc_seed();
                const double current = model->drive().current(seed.time_s + .5 * dt);
                auto cE = g.curl(seed.E); gmg::Vec Bm = seed.B; for (int f = 0; f < g.nb; ++f) Bm[f] -= .5 * dt * cE[f];
                const auto amp = g.curl_adjoint(Bm);
                for (int e = 0; e < g.ne; ++e) if (!g.fixed_E[e])
                    b[e] = dt * (amp[e] - current * in.feed[e] - (S[e] + g.Y[e]) * seed.E[e] + in.background_force[e]) / std::sqrt(g.M[e]);
            }
            return b;
        };
        if (nvrtc == "none") {   // CPU-only mode (sanitizer builds): CPU reference solves with the independent check
            J rows = J::array();
            for (double dt : dts) {
                gmg::set_dt(H, dt);
                auto jd = gmg::jacobi_diagonal(g, S, dt); gmg::Vec jinv(g.ne); for (int e = 0; e < g.ne; ++e) jinv[e] = 1. / jd[e];
                const auto b = make_rhs(dt);
                J row{{"dt_s", dt}};
                for (const char* mode : {"jacobi", "gmg"}) {
                    gmg::SolveStats st; auto x = gmg::cpu_solve(H, std::string(mode) == "jacobi" ? jinv : gmg::Vec{}, b, 2000, st);
                    gmg::Vec Ax; host_apply.apply(S, dt, x, Ax); gmg::Vec r(g.ne); for (int e = 0; e < g.ne; ++e) r[e] = b[e] - Ax[e];
                    const double target = std::max(1e-30, 2e-13 * foil::norm(b)), res = foil::norm(r);
                    row[mode] = {{"iterations", st.iterations}, {"independent_residual", res}, {"target", target}, {"passed", res <= target}};
                    std::cout << "cpu dt " << dt << " " << mode << " its " << st.iterations << " ratio " << res / target << std::endl;
                }
                rows.push_back(row);
            }
            summary["rows"] = rows; summary["cpu_only"] = true;
            std::ofstream(out / "gmg.json") << summary.dump(1) << '\n';
            bool ok = true; for (const auto& r : rows) ok = ok && r["jacobi"]["passed"].get<bool>() && r["gmg"]["passed"].get<bool>();
            return ok ? 0 : 2;
        }
            const auto t_gpu = steady::now();
        hall_gpu::SessionOptions so; so.require_uuid = uuid; so.nvrtc_library = nvrtc;
        hall_gpu::Session session(so);
        J rows = J::array();
        std::unique_ptr<gmg::GpuSolver> gpu;
        for (double dt : dts) {
            J row{{"dt_s", dt}, {"courant_transverse", 299792458. * dt / g.dx}};
            t = steady::now(); gmg::set_dt(H, dt); row["set_dt_ms"] = ms_since(t);
            auto jd = gmg::jacobi_diagonal(g, S, dt); gmg::Vec jinv(g.ne); for (int e = 0; e < g.ne; ++e) jinv[e] = 1. / jd[e];
            t = steady::now(); gpu.reset(); gpu = std::make_unique<gmg::GpuSolver>(session, H, &jinv); row["gpu_setup_ms"] = ms_since(t);
            if (!summary.contains("compile")) {
                J o = J::array(); for (const auto& x : gpu->compile_record.options) o.push_back(x);
                summary["compile"] = {{"options", o}, {"arch_option", gpu->compile_record.arch_option}, {"image_kind", gpu->compile_record.image_kind},
                                      {"nvrtc", std::to_string(gpu->compile_record.nvrtc_major) + "." + std::to_string(gpu->compile_record.nvrtc_minor)},
                                      {"device", session.device().name}};
            }
            const auto b = make_rhs(dt);
            auto check = [&](const gmg::Vec& x) {
                gmg::Vec Ax; host_apply.apply(S, dt, x, Ax); gmg::Vec r(g.ne); for (int e = 0; e < g.ne; ++e) r[e] = b[e] - Ax[e];
                const double nb = foil::norm(b), target = std::max(1e-30, 2e-13 * nb), res = foil::norm(r);
                return J{{"independent_residual", res}, {"target", target}, {"ratio", res / target}, {"passed", res <= target}};
            };
            for (const char* mode : {"jacobi", "gmg"}) {
                const bool jac = std::string(mode) == "jacobi";
                gmg::SolveStats st; t = steady::now(); auto x = gpu->solve(b, 2000, jac, st); const double ms = ms_since(t);
                J r{{"iterations", st.iterations}, {"cycle_iterations", st.cycle_iterations}, {"corrections", st.corrections}, {"solve_ms", ms},
                    {"backend_residual", st.backend_residual}, {"target", st.target}, {"check", check(x)}};
                if (!jac) {
                    gmg::SolveStats st2; auto x2 = gpu->solve(b, 2000, false, st2);
                    r["repeat_bitwise_identical"] = x2 == x && st2.iterations == st.iterations;
                    if (cpuref) {
                        gmg::SolveStats sc; t = steady::now(); auto xc = gmg::cpu_solve(H, {}, b, 2000, sc);
                        r["cpu_reference_ms"] = ms_since(t); r["cpu_gpu_bitwise_identical"] = std::memcmp(xc.data(), x.data(), x.size() * 8) == 0;
                        r["cpu_iterations"] = sc.iterations;
                    }
                    r["vcycle_ms"] = gpu->time_vcycle(5); r["apply_ms"] = gpu->time_apply(5);
                } else if (cpuref) {
                    gmg::SolveStats sc; auto xc = gmg::cpu_solve(H, jinv, b, 2000, sc);
                    r["cpu_gpu_bitwise_identical"] = std::memcmp(xc.data(), x.data(), x.size() * 8) == 0; r["cpu_iterations"] = sc.iterations;
                }
                row[mode] = r;
                std::cout << "dt " << dt << " " << mode << " its " << st.iterations << " ms " << ms << " ratio " << r["check"]["ratio"] << std::endl;
            }
            row["memory"] = {{"own_peak_device_bytes", session.memory().own_peak_bytes}, {"own_bytes_in_use", session.memory().own_bytes_in_use}};
            rows.push_back(row);
            std::ofstream(out / "gmg.json") << J{{"partial", rows}}.dump(1) << '\n';
        }
        gpu.reset();
        summary["rows"] = rows;
        summary["gpu_hold_ms"] = ms_since(t_gpu);
        summary["max_rss_kB"] = rss_kb();
        session.close();
        std::ofstream(out / "gmg.json") << summary.dump(1) << '\n';
        std::cout << summary.dump(1).substr(0, 4000) << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
}
