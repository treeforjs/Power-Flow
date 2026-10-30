// CPU side of the Hall GPU field backend qualification.
//
//   hall_field_cpu prep  TRACE OUTDIR TRANSVERSE AXIAL AIR WIDTH REPEATS
//   hall_field_cpu guard TRACE OUTDIR TRANSVERSE AXIAL AIR WIDTH GPUDIR
//
// prep:  builds the sealed Model, the DC-seed and heated fixtures and, inside
//        Model::prepare (through model_seam.hpp), assembles the operator,
//        compares it with the sealed matrix-free apply, times the sealed PCG
//        and writes one operator/RHS fixture per (state, dt).
// guard: rebuilds the same fixtures, re-enters Model::prepare once with the
//        sealed CPU increment and once per GPU result. The backend checks the
//        RHS/diagonal identity, recomputes the residual with the sealed apply
//        and returns the GPU increment; Model::prepare then evaluates its own
//        Gauss, divB, Gauss-increment and energy guards on it. The v2 thermal
//        owner evaluates the Joule transfer of every candidate.
#include "model_seam.hpp"
#include "assemble.hpp"
#include <mpi.h>
#include <chrono>
#include <iostream>

namespace h = hall_corrective;
using h::Json;
using axial::Vec;
using Clock = std::chrono::steady_clock;

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
            {"max_rel", d.max_rel}, {"normwise_rel", d.norm_rel}, {"max_abs_over_absA_absv", d.max_scaled}};
}

// Test fixture declaration (same construction as test_field.cpp): the DC seed's
// induced charge is declared as conductor/PEC charge on material/short nodes
// and as exterior charge on the EM boundary. It is a fixture, not a genesis.
h::ChargeLedger declared_dc(const h::Model& m, const axial::State& s) {
    const auto& g = m.mesh();
    h::ChargeLedger q;
    q.particles.net_C.resize(g.nn); q.particles.gross_C.resize(g.nn);
    q.conductor_C.resize(g.nn); q.exterior_C.resize(g.nn);
    q.conductor_declaration = h::sha256("GPU field fixture explicitly declares DC induced conductor/PEC surface charge");
    q.exterior_declaration = h::sha256("GPU field fixture explicitly declares DC exterior termination surface charge");
    auto field = h::field_charge(g, s);
    for (int n = 0; n < g.nn; ++n) {
        if (g.node_body[n] >= 0 || g.short_node[n]) q.conductor_C[n] = field[n];
        else if (g.outer_node[n]) q.exterior_C[n] = field[n];
    }
    return q;
}

// Synthetic nonuniform Joule load (W per cell) applied through the v2 thermal
// owner: stronger near the gap-facing foil surface, near the foil x-edges,
// toward the feed end, and on the upper foil. Only the temperature field it
// produces is used; sigma(T) then comes from the owner's own material law.
constexpr double heat_peak_W_m3 = 1.4e23, heat_dt_s = 1e-14;
constexpr int heat_steps = 4;
Vec heating_profile(const axial::Mesh& g) {
    Vec P(g.cells.size());
    for (std::size_t c = 0; c < g.cells.size(); ++c) {
        const auto& cell = g.cells[c];
        if (cell.body < 0) continue;
        const double x = (cell.i + .5) * g.dx - g.side / 2, y = (cell.j + .5) * g.dy - g.side / 2;
        const double depth = std::abs(y) - .0025, edge = .0025 - std::abs(x);
        const double kk = double(cell.k - g.kshort) / std::max(1, g.nz_foil - 1);
        const double shape = (.2 + .5 * std::exp(-depth / 50e-6) + .3 * std::exp(-edge / .3e-3)) * (.6 + .4 * kk) * (cell.body ? 1. : .85);
        P[c] = heat_peak_W_m3 * shape * cell.volume;
    }
    return P;
}
void heat(hall_driver::Thermal& thermal, const axial::Mesh& g) {
    const Vec P = heating_profile(g);
    for (int step = 0; step < heat_steps; ++step) {
        axial::Receipt r;
        r.cell_joule_W = P;
        double total = 0;
        for (double p : P) total += p;
        r.joule_J = heat_dt_s * total;
        auto trial = thermal.prepare(r, heat_dt_s, h::sha256("GPU field fixture synthetic Joule heating step " + std::to_string(step)));
        h::need(thermal.commit(std::move(trial)), "synthetic heating commit");
    }
}

struct Fixture {
    std::string name;
    const axial::State* old;
    const h::ChargeLedger* charge;
    const Vec* sigma;
    double dt;
    hall_driver::Thermal* thermal;
};

// Instrumented copy of foil::solve (foil_em.hpp:56-63), same arithmetic.
struct PcgRecord { std::vector<int> iterations; std::vector<std::vector<double>> alpha, beta, rnorm; double true_norm = 0, target = 0; };
Vec instrumented_solve(const hgf_seam::Apply& apply, const Vec& b, const Vec& diagonal, int maxit, PcgRecord& rec) {
    using foil::dot; using foil::norm;
    Vec x(b.size()), r = b;
    const double target = std::max(1e-30, 2e-13 * norm(b));
    rec.target = target;
    for (int correction = 0; correction < 3; ++correction) {
        rec.alpha.emplace_back(); rec.beta.emplace_back(); rec.rnorm.emplace_back();
        Vec z(r.size());
        for (std::size_t i = 0; i < r.size(); ++i) z[i] = r[i] / diagonal[i];
        Vec p = z;
        double rz = dot(r, z);
        int iteration = 0;
        for (; iteration < maxit && norm(r) > target * .2; ++iteration) {
            rec.rnorm.back().push_back(norm(r));
            auto Ap = apply(p);
            double denominator = dot(p, Ap);
            hall::require(denominator > 0 && std::isfinite(denominator), "non-SPD field operator");
            double alpha = rz / denominator;
            for (std::size_t i = 0; i < r.size(); ++i) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; z[i] = r[i] / diagonal[i]; }
            double nrz = dot(r, z), beta = nrz / rz;
            for (std::size_t i = 0; i < r.size(); ++i) p[i] = z[i] + beta * p[i];
            rz = nrz;
            rec.alpha.back().push_back(alpha); rec.beta.back().push_back(beta);
        }
        rec.rnorm.back().push_back(norm(r));
        rec.iterations.push_back(iteration);
        r = foil::subtract(b, apply(x));
        rec.true_norm = norm(r);
        if (norm(r) <= target) { hall::require(foil::finite(x), "nonfinite field solve"); return x; }
    }
    throw std::runtime_error("field solve independent residual failed");
}

Json pcg_json(const PcgRecord& r) {
    int total = 0;
    for (int i : r.iterations) total += i;
    auto ritz = hgf::lanczos_extremes(r.alpha.at(0), r.beta.at(0));
    return {{"iterations_per_correction", r.iterations}, {"iterations_total", total}, {"corrections", r.iterations.size()},
            {"residual_history_recursive", r.rnorm}, {"final_true_residual", r.true_norm}, {"target", r.target},
            {"lanczos", {{"steps", ritz.k}, {"lambda_min", ritz.lmin}, {"lambda_max", ritz.lmax}, {"kappa_estimate", ritz.kappa}}}};
}
Json trace_json(const hgf::CgTrace& t) {
    auto ritz = hgf::lanczos_extremes(t.alpha.at(0), t.beta.at(0));
    return {{"iterations_per_correction", t.iterations}, {"iterations_total", t.total_iterations()}, {"corrections", t.corrections},
            {"bnorm_fixed_tree", t.bnorm}, {"target_fixed_tree", t.target}, {"final_true_residual_fixed_tree", t.true_norm},
            {"residual_history_recursive", t.rnorm},
            {"lanczos", {{"steps", ritz.k}, {"lambda_min", ritz.lmin}, {"lambda_max", ritz.lmax}, {"kappa_estimate", ritz.kappa}}}};
}

// Guard record for one candidate: everything below is computed by the sealed
// model code (Model::prepare and admit_constraints ran already and would have
// thrown on failure); bounds are re-derived here only for reporting.
Json guard_record(const h::Model& model, const axial::State& old, const h::FieldCandidate& c, hall_driver::Thermal& thermal) {
    const auto& g = model.mesh();
    const auto& r = c.receipt();
    const auto Q0 = h::field_charge(g, old), Q1 = h::field_charge(g, c.field());
    double activityQ = 0;
    for (int n = 0; n < g.nn; ++n) activityQ = std::max(activityQ, double(h::magnitude(Q1[n] - Q0[n])));
    const double activity = std::abs(r.field_delta_J) + std::abs(r.plasma_work_J) + std::abs(r.joule_J) + std::abs(r.outer_loss_J) +
                            std::abs(r.outer_supply_J) + std::abs(r.circuit_work_J);
    const auto& k = c.constraints();
    // Joule transfer into the v2 thermal owner (axial_thermal_v2.hpp:103-104
    // rejects a mismatched Joule operand). The trial is discarded.
    auto trial = thermal.prepare(r, c.dt_s(), c.after_sha256());
    const double transfers = std::abs(r.field_delta_J) + std::abs(r.plasma_work_J) + std::abs(trial.thermal_delta_J) +
                             std::abs(r.outer_loss_J) + std::abs(r.outer_supply_J) + std::abs(r.circuit_work_J);
    const double field_heat = r.balance_J + trial.thermal_delta_J - r.joule_J;
    const double field_heat_bound = 2e-20 + 1e-10 * transfers + trial.roundoff_bound_J;
    return {{"model_prepare_returned", true},
            {"absolute_gauss", {{"max_defect_C", k.max_gauss_C}, {"max_bound_C", k.max_gauss_bound_C}, {"passed", k.max_gauss_C <= k.max_gauss_bound_C}}},
            {"absolute_divB", {{"max_T_m", k.max_divB_T_m}, {"max_bound_T_m", k.max_divB_bound_T_m}, {"passed", k.max_divB_T_m <= k.max_divB_bound_T_m}}},
            {"gauss_increment", {{"value_C", r.gauss_increment_C}, {"bound_C", 1e-22 + 2e-10 * activityQ}, {"passed", r.gauss_increment_C <= 1e-22 + 2e-10 * activityQ}}},
            {"field_energy", {{"balance_J", r.balance_J}, {"bound_J", 1e-20 + 1e-10 * activity}, {"passed", std::abs(r.balance_J) <= 1e-20 + 1e-10 * activity}}},
            {"joule_transfer", {{"field_joule_J", r.joule_J}, {"thermal_admitted_joule_J", trial.admitted_joule_J},
                                {"bound_J", 1e-22 + 2e-12 * std::abs(r.joule_J)}, {"thermal_owner_accepted", true}}},
            {"field_heat_energy", {{"balance_J", field_heat}, {"bound_J", field_heat_bound}, {"passed", std::abs(field_heat) <= field_heat_bound},
                                   {"formula", "owner.hpp:20-24 FieldHeatProposal bound, evaluated here"}}},
            {"receipt", r.json()},
            {"after_sha256", c.after_sha256()}};
}

struct Mesh0 { int transverse, axial, air; double width; };
}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rc = 0;
    try {
        h::need(argc == 9, "usage: hall_field_cpu prep|guard TRACE OUTDIR TRANSVERSE AXIAL AIR WIDTH REPEATS|GPUDIR");
        const std::string mode = argv[1];
        h::need(mode == "prep" || mode == "guard", "mode must be prep or guard");
        const std::filesystem::path trace = argv[2], out = argv[3];
        const Mesh0 spec{std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), std::stod(argv[7])};
        const int repeats = mode == "prep" ? std::stoi(argv[8]) : 0;
        const std::filesystem::path gpudir = mode == "guard" ? argv[8] : "";
        std::filesystem::create_directories(out);

        Json report{{"mode", mode}, {"argv", std::vector<std::string>(argv, argv + argc)},
                    {"mesh_spec", {{"transverse", spec.transverse}, {"axial", spec.axial}, {"air", spec.air}, {"width_m", spec.width}}}};
        auto t0 = Clock::now();
        h::Model model(h::MeshSpec{spec.transverse, spec.axial, spec.air, spec.width}, trace);
        report["model_construction_ms"] = ms_since(t0);
        const auto run_start = t0;
        const auto& g = model.mesh();
        int fixed = 0;
        for (int e = 0; e < g.ne; ++e) fixed += g.fixed_E[e];
        report["mesh"] = {{"edges", g.ne}, {"faces", g.nb}, {"nodes", g.nn}, {"cells", g.cells.size()}, {"pec_edges", fixed},
                          {"nx", g.nx}, {"nz", g.nz}, {"dx_m", g.dx}, {"dz_m", g.dz}};
        report["configuration_sha256"] = model.configuration_sha256();

        // States
        const axial::State dc = model.dc_seed();
        const auto charge = declared_dc(model, dc);
        const Vec sigma0 = model.initial_conductivity();
        hall_driver::Thermal thermal_dc(g, model.drive().start());
        hall_driver::Thermal thermal_hot(g, model.drive().start());
        t0 = Clock::now();
        heat(thermal_hot, g);
        report["heating_ms"] = ms_since(t0);
        const Vec sigma_hot = thermal_hot.conductivity();
        axial::State hot = dc;
        hot.time_s = thermal_hot.time_s();
        h::need(hot.time_s == model.drive().start() + heat_steps * heat_dt_s || hot.time_s > model.drive().start(), "thermal clock");
        {
            double smin = 1e300, smax = 0, s0min = 1e300, s0max = 0, Tmin = 1e300;
            for (std::size_t c = 0; c < g.cells.size(); ++c)
                if (g.cells[c].body >= 0) {
                    smin = std::min(smin, sigma_hot[c]); smax = std::max(smax, sigma_hot[c]);
                    s0min = std::min(s0min, sigma0[c]); s0max = std::max(s0max, sigma0[c]);
                }
            for (const auto& st : thermal_hot.states()) for (double T : st.T) Tmin = std::min(Tmin, T);
            report["heated_state"] = {{"heating_steps", heat_steps}, {"heating_dt_s", heat_dt_s}, {"peak_power_density_W_m3", heat_peak_W_m3},
                                      {"thermal_time_s", thermal_hot.time_s()}, {"thermal_revision", thermal_hot.revision()},
                                      {"T_min_K", Tmin}, {"T_max_K", thermal_hot.max_temperature_K()},
                                      {"sigma_min_S_m", smin}, {"sigma_max_S_m", smax},
                                      {"sigma_300K_min_S_m", s0min}, {"sigma_300K_max_S_m", s0max},
                                      {"sigma_sha256", h::sha256(Json(sigma_hot).dump())}};
        }
        const auto initial_constraints = model.admit(dc, charge);
        report["initial_constraints"] = {{"max_gauss_C", initial_constraints.max_gauss_C}, {"max_gauss_bound_C", initial_constraints.max_gauss_bound_C},
                                         {"max_divB_T_m", initial_constraints.max_divB_T_m}, {"max_divB_bound_T_m", initial_constraints.max_divB_bound_T_m}};

        std::vector<Fixture> fixtures{{"dc-10fs", &dc, &charge, &sigma0, 1e-14, &thermal_dc},
                                      {"dc-1fs", &dc, &charge, &sigma0, 1e-15, &thermal_dc},
                                      {"heated-10fs", &hot, &charge, &sigma_hot, 1e-14, &thermal_hot},
                                      {"heated-1fs", &hot, &charge, &sigma_hot, 1e-15, &thermal_hot}};
        Json cases = Json::object();
        for (const auto& fx : fixtures) {
            const std::filesystem::path dir = out / fx.name;
            std::filesystem::create_directories(dir);
            Json rec{{"fixture", fx.name}, {"dt_s", fx.dt}, {"state_time_s", fx.old->time_s}};
            if (mode == "prep") {
                hgf_seam::backend() = [&](const hgf_seam::Apply& apply, const Vec& b, const Vec& diagonal, int maxit) {
                    h::need(foil::finite(b) && b.size() == std::size_t(g.ne), "seam RHS");
                    rec["maxit"] = maxit;
                    auto ta = Clock::now();
                    auto A = hgf::assemble(g, fx.dt, *fx.sigma);
                    rec["assembly_ms"] = ms_since(ta);
                    const auto dd = hgf::compare(A.diag, diagonal);
                    h::need(dd.bitwise, "assembled preconditioner diagonal differs from Model::prepare's diagonal");
                    rec["operator"] = {{"rows", A.n}, {"faces", A.nb}, {"nnz_assembled", A.nnz()}, {"nnz_curl_transpose", A.ct_face.size()},
                                       {"symmetric_bitwise", hgf::csr_symmetric(A)}, {"diagonal_bitwise_equal_to_model", dd.bitwise},
                                       {"sigma_dependent_rows", A.sigma_rows.size()}};
                    // CSR / factored products versus the sealed matrix-free apply.
                    Json checks = Json::array();
                    Vec Bh, y(b.size());
                    std::vector<std::pair<std::string, Vec>> probes{{"rhs", b}, {"pseudo_random_seed_1", hgf::pseudo_random(b.size(), 1)},
                                                                     {"pseudo_random_seed_2", hgf::pseudo_random(b.size(), 2)}};
                    auto check_probe = [&](const std::string& label, const Vec& v) {
                        const Vec ref = apply(v), scale = hgf::abs_product(A, v);
                        hgf::spmv_csr(A, v.data(), y.data());
                        auto dcsr = hgf::compare(y, ref, &scale);
                        hgf::apply_factored(A, v.data(), y.data(), Bh);
                        auto dfac = hgf::compare(y, ref, &scale);
                        checks.push_back({{"vector", label}, {"assembled_csr_vs_matrix_free", diff_json(dcsr)},
                                          {"factored_vs_matrix_free", diff_json(dfac)}});
                        return dfac.bitwise;
                    };
                    bool factored_exact = true;
                    for (const auto& [label, v] : probes) factored_exact = check_probe(label, v) && factored_exact;
                    // Sealed CPU PCG, timed.
                    std::vector<double> sealed_ms;
                    Vec x;
                    for (int i = 0; i < repeats + 1; ++i) {
                        auto ts = Clock::now();
                        Vec xi = hgf_seam::sealed_solve(apply, b, diagonal, maxit);
                        const double t = ms_since(ts);
                        if (i == 0) x = xi; else { sealed_ms.push_back(t); h::need(hgf::compare(xi, x).bitwise, "sealed solve repeat differs"); }
                    }
                    factored_exact = check_probe("sealed_solution", x) && factored_exact;
                    rec["spmv_checks"] = checks;
                    rec["factored_bitwise_on_all_probes"] = factored_exact;
                    PcgRecord pr;
                    Vec xi = instrumented_solve(apply, b, diagonal, maxit, pr);
                    h::need(hgf::compare(xi, x).bitwise, "instrumented copy differs from sealed foil::solve");
                    const double bnorm = foil::norm(b), target = std::max(1e-30, 2e-13 * bnorm);
                    const double true_res = foil::norm(foil::subtract(b, apply(x)));
                    rec["cpu_sealed_pcg"] = {{"timing", timing(sealed_ms)}, {"record", pcg_json(pr)}, {"instrumented_copy_bitwise_equal", true},
                                             {"rhs_norm_canonical", bnorm}, {"target_canonical", target}, {"independent_residual", true_res},
                                             {"passed", true_res <= target}};
                    // Plain CPU PCG on the preassembled operators (timing baseline).
                    auto time_plain = [&](bool assembled) {
                        std::vector<double> ms;
                        hgf::CgTrace tr;
                        Vec xp;
                        for (int i = 0; i < repeats + 1; ++i) {
                            hgf::CgTrace t;
                            auto ts = Clock::now();
                            xp = assembled ? hgf::plain_cpu_pcg([&](const double* v, double* o) { hgf::spmv_csr(A, v, o); }, b, A.diag, maxit, t)
                                           : hgf::plain_cpu_pcg([&](const double* v, double* o) { hgf::apply_factored(A, v, o, Bh); }, b, A.diag, maxit, t);
                            if (i) ms.push_back(ms_since(ts));
                            tr = t;
                        }
                        const double res = foil::norm(foil::subtract(b, apply(xp)));
                        return Json{{"timing", timing(ms)}, {"iterations_per_correction", tr.iterations}, {"iterations_total", tr.total_iterations()},
                                    {"independent_residual_sealed_apply", res}, {"passed", res <= target}};
                    };
                    rec["cpu_plain_pcg_assembled"] = time_plain(true);
                    rec["cpu_plain_pcg_factored"] = time_plain(false);
                    // Host emulation of the device algorithm (expected device bits).
                    hgf::CgTrace ta_tr, tf_tr;
                    auto te = Clock::now();
                    Vec xa = hgf::emulate_device_cg([&](const double* v, double* o) { hgf::spmv_csr(A, v, o); }, b, A.diag, maxit, ta_tr);
                    Vec xf = hgf::emulate_device_cg([&](const double* v, double* o) { hgf::apply_factored(A, v, o, Bh); }, b, A.diag, maxit, tf_tr);
                    rec["emulation_ms"] = ms_since(te);
                    const double ra = foil::norm(foil::subtract(b, apply(xa))), rf = foil::norm(foil::subtract(b, apply(xf)));
                    rec["device_algorithm_emulation"] = {{"assembled", {{"trace", trace_json(ta_tr)}, {"independent_residual_sealed_apply", ra}, {"passed", ra <= target},
                                                                        {"x_vs_sealed", diff_json(hgf::compare(xa, x))}, {"x_sha256", hgf::sha256_vec(xa)}}},
                                                         {"factored", {{"trace", trace_json(tf_tr)}, {"independent_residual_sealed_apply", rf}, {"passed", rf <= target},
                                                                       {"x_vs_sealed", diff_json(hgf::compare(xf, x))}, {"x_sha256", hgf::sha256_vec(xf)}}}};
                    hgf::ArrayFile f;
                    hgf::store_operator(f, A);
                    f.put("b", b); f.put("x_cpu", x); f.put("x_emu_assembled", xa); f.put("x_emu_factored", xf);
                    f.put("maxit", std::vector<std::int32_t>{maxit});
                    f.put_scalar("target_canonical", target);
                    f.put_scalar("rhs_norm_canonical", bnorm);
                    const auto path = (dir / "fixture.bin").string();
                    f.write(path);
                    rec["fixture_file"] = path;
                    rec["fixture_sha256"] = hgf::sha256_file(path);
                    rec["rhs_sha256"] = hgf::sha256_vec(b);
                    rec["diagonal_sha256"] = hgf::sha256_vec(diagonal);
                    rec["x_cpu_sha256"] = hgf::sha256_vec(x);
                    return x;
                };
                auto c = model.prepare(*fx.old, *fx.charge, Vec(g.ne), fx.charge->particles, *fx.sigma, fx.dt);
                rec["cpu_path_guards"] = guard_record(model, *fx.old, c, *fx.thermal);
            } else {
                const auto fixture_path = (dir / "fixture.bin").string();
                const auto F = hgf::ArrayFile::read(fixture_path);
                rec["fixture_sha256"] = hgf::sha256_file(fixture_path);
                const Vec fb = F.f64("b"), fdiag = F.f64("diag"), xcpu = F.f64("x_cpu");
                Json variants = Json::object();
                std::vector<std::pair<std::string, std::string>> inputs{{"cpu_sealed", ""}, {"gpu_assembled", "x_gpu_assembled.bin"}, {"gpu_factored", "x_gpu_factored.bin"}};
                Vec x_ref;
                for (const auto& [label, file] : inputs) {
                    Vec xin;
                    Json v;
                    if (file.empty()) xin = xcpu;
                    else {
                        const auto p = (gpudir / fx.name / file).string();
                        const auto X = hgf::ArrayFile::read(p);
                        xin = X.f64("x");
                        v["input_file"] = p;
                        v["input_sha256"] = hgf::sha256_file(p);
                    }
                    v["x_sha256"] = hgf::sha256_vec(xin);
                    hgf_seam::backend() = [&](const hgf_seam::Apply& apply, const Vec& b, const Vec& diagonal, int) {
                        h::need(hgf::compare(b, fb).bitwise, "Model::prepare RHS differs from the stored fixture RHS");
                        h::need(hgf::compare(diagonal, fdiag).bitwise, "Model::prepare diagonal differs from the stored fixture");
                        h::need(xin.size() == b.size() && foil::finite(xin), "backend increment shape/finiteness");
                        const double bnorm = foil::norm(b), target = std::max(1e-30, 2e-13 * bnorm);
                        const double res = foil::norm(foil::subtract(b, apply(xin)));
                        v["independent_residual"] = {{"rhs_norm", bnorm}, {"target", target}, {"residual", res}, {"ratio_to_target", res / target},
                                                     {"passed", res <= target}, {"operator", "sealed matrix-free apply of Model::prepare"},
                                                     {"norm", "foil::norm (long double accumulation, canonical index order)"}};
                        h::need(res <= target, "GPU increment fails the independent FP64 residual bound");
                        return xin;
                    };
                    auto c = model.prepare(*fx.old, *fx.charge, Vec(g.ne), fx.charge->particles, *fx.sigma, fx.dt);
                    v["guards"] = guard_record(model, *fx.old, c, *fx.thermal);
                    if (label == "cpu_sealed") x_ref = xin;
                    v["increment_vs_cpu_sealed"] = diff_json(hgf::compare(xin, x_ref));
                    variants[label] = v;
                }
                rec["variants"] = variants;
                // Negative controls on the GPU (factored) increment.
                // (1) Mandatory: one component perturbed by 1e-6 relative must be
                //     refused by the independent residual gate.
                // (2) Characterization: the whole increment scaled by s with the
                //     residual gate bypassed; records whether Model::prepare's own
                //     guards refuse it (they are not designed as solver-accuracy checks).
                Json controls = Json::object();
                {
                    const auto X = hgf::ArrayFile::read((gpudir / fx.name / "x_gpu_factored.bin").string());
                    const Vec xg = X.f64("x");
                    std::size_t imax = 0;
                    for (std::size_t i = 0; i < xg.size(); ++i) if (std::abs(xg[i]) > std::abs(xg[imax])) imax = i;
                    Json scaled = Json::array();
                    for (double factor : {0.0, 1.001, 1.1, 2.0}) {
                        const bool gate = factor == 0.0;
                        Vec xin = xg;
                        if (gate) xin[imax] *= 1 + 1e-6; else for (double& v : xin) v *= factor;
                        Json c{{"perturbation", gate ? "largest component x (1+1e-6), residual gate active" : "whole increment scaled, residual gate bypassed"}};
                        if (!gate) c["scale"] = factor;
                        hgf_seam::backend() = [&](const hgf_seam::Apply& apply, const Vec& b, const Vec&, int) {
                            const double target = std::max(1e-30, 2e-13 * foil::norm(b));
                            const double res = foil::norm(foil::subtract(b, apply(xin)));
                            c["independent_residual"] = res; c["target"] = target; c["ratio_to_target"] = res / target;
                            if (gate) h::need(res <= target, "GPU increment fails the independent FP64 residual bound");
                            return xin;
                        };
                        try {
                            auto cand = model.prepare(*fx.old, *fx.charge, Vec(g.ne), fx.charge->particles, *fx.sigma, fx.dt);
                            const auto& r = cand.receipt();
                            const double activity = std::abs(r.field_delta_J) + std::abs(r.plasma_work_J) + std::abs(r.joule_J) +
                                                    std::abs(r.outer_loss_J) + std::abs(r.outer_supply_J) + std::abs(r.circuit_work_J);
                            c["rejected"] = false;
                            c["accepted_energy_balance_J"] = r.balance_J; c["energy_bound_J"] = 1e-20 + 1e-10 * activity;
                            c["accepted_gauss_increment_C"] = r.gauss_increment_C; c["accepted_ampere_residual_A_m"] = r.ampere_residual_A_m;
                        } catch (const std::exception& e) {
                            c["rejected"] = true; c["message"] = e.what();
                        }
                        if (gate) {
                            h::need(c["rejected"].get<bool>(), "residual-gate negative control was accepted");
                            controls["residual_gate"] = c;
                        } else scaled.push_back(c);
                    }
                    controls["model_guards_without_residual_gate"] = scaled;
                }
                rec["negative_controls"] = controls;
            }
            cases[fx.name] = rec;
            std::cout << fx.name << " done" << std::endl;
        }
        report["fixtures"] = cases;
        report["seam_calls"] = hgf_seam::calls();
        report["total_ms"] = ms_since(run_start);
        std::ofstream(out / (mode + ".json")) << report.dump(2) << '\n';
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        rc = 1;
    }
    MPI_Finalize();
    return rc;
}
