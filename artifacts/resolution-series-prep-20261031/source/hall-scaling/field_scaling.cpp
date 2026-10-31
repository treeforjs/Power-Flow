// Resolution-series preparation: mesh series and field-only scaling of the sealed Hall corrective v3.1
// (read in place, unchanged; this program only constructs its public types with other MeshSpec values).
//
//   field_scaling mesh-report OUT N AXIAL AIR [N AXIAL AIR ...]
//       geometry, identity and quality of Hall axial meshes (no field solve)
//   field_scaling admit       OUT TRACE N AXIAL AIR [WIDTH]
//       instrumented copy of axial::dc_initial (iteration counts, residual history, time), then the real
//       Model construction (Hall admission: mesh check, DC seed, identities); initial fields compared bitwise
//   field_scaling dc-floor    OUT TRACE N AXIAL AIR EXTRA
//       DC seed solves: sealed cycles, EXTRA further cycles, attainable accuracy of the double apply (long double reference)
//   field_scaling graded-dc   OUT TRACE H CORE_HALF BOX_HALF Q HMAX
//       graded transverse mesh (graded/graded_mesh.hpp): sizes, resident memory model, DC port inductance 2 W_B / I^2
//   field_scaling graded-sizes OUT H1,H2,.. BOXHALF1,.. Q1,.. HMAX
//   field_scaling cpu-vacuum  OUT TRACE N AXIAL AIR STEPS DT
//       sealed v3 Owner (cpu-pcg) field-only steps
//   field_scaling gpu-vacuum  OUT TRACE NVRTC UUID N AXIAL AIR STEPS DT DTMAX [EVERY [WIDTH]]
//       v3.1 ResidentOwner field-only steps (device CG), per-step timings, CG iterations, residuals, memory,
//       per-node absolute Gauss ratio drift
//   field_scaling gpu-dtscan  OUT TRACE NVRTC UUID N AXIAL AIR STEPS DT1 [DT2 ...]
//       iteration count and time per step versus dt (fresh owner from the DC seed for every dt)
//   field_scaling gpu-coupled OUT TRACE NVRTC UUID N AXIAL AIR STEPS CARRIERS DT
//       closed coupled Hall step (GPU push, no native/chemistry) with CARRIERS co-located electron/H2O+ pairs
//
// GPU modes run under the GPU lock wrappers (flock testing-gpu.lock flock testing-gpu-<UUID>.lock, then the core lease);
// build_run.py refuses a GPU mode unless both lock files are held. The resident descriptor is computed on a separate, closed probe session with the
// default 45-cell mesh (the descriptor depends only on the kernels, options and device), so the measured device
// memory is that of one resident field.
#include "resident_owner.hpp"
#include "gpu/gpu_backends.hpp"
#include "hall-finite-foil-em-20261029/axial_particles.hpp"
namespace hc = hall_corrective_v3;
#include "motion_fixtures.hpp"
#include "../graded/graded_mesh.hpp"
#include <algorithm>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <sys/file.h>
#include <unistd.h>

namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
namespace gp = hall_corrective_v3::gpu;
using J = hc::Json;
using steady = std::chrono::steady_clock;

namespace {
double ms_since(steady::time_point t) { return std::chrono::duration<double, std::milli>(steady::now() - t).count(); }
long rss_kb(const char* key) {
    std::ifstream f("/proc/self/status"); std::string line;
    while (std::getline(f, line)) if (line.rfind(key, 0) == 0) return std::stol(line.substr(std::strlen(key)));
    return -1;
}
struct Stat {
    std::vector<double> v;
    void add(double x) { v.push_back(x); }
    J json() const {
        if (v.empty()) return J::object();
        auto s = v; std::sort(s.begin(), s.end()); double sum = 0; for (double x : s) sum += x;
        return {{"n", s.size()}, {"median", s[s.size() / 2]}, {"mean", sum / s.size()}, {"p10", s[s.size() / 10]}, {"p90", s[s.size() * 9 / 10]},
                {"min", s.front()}, {"max", s.back()}};
    }
};
struct Stats {
    std::map<std::string, Stat> m;
    void add(const std::string& k, double x) { m[k].add(x); }
    J json() const { J j = J::object(); for (const auto& [k, s] : m) j[k] = s.json(); return j; }
};
struct Boundary { std::string field_sha256, heat_sha256, source_sha256, mesh_sha256, material_sha256, common_clock_sha256; double time_s = 0; std::uint64_t accepted_step = 0; };
// Declared DC ledger (same recipe as the Hall v3.1 tests: conductor/PEC and exterior induced charge).
hc::ChargeLedger declared(const hc::Model& m, const axial::State& s, hc::FreeCharge particles) {
    const auto& g = m.mesh(); hc::ChargeLedger q; q.particles = std::move(particles);
    q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
    q.conductor_declaration = hc::sha256("parity: declared DC conductor/PEC induced charge");
    q.exterior_declaration = hc::sha256("parity: declared DC exterior termination charge");
    const auto field = hc::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = field[k]; else if (g.outer_node[k]) q.exterior_C[k] = field[k]; }
    return q;
}
std::string hexbits(double x) { std::ostringstream o; o << std::hex << std::setw(16) << std::setfill('0') << std::bit_cast<std::uint64_t>(x); return o.str(); }
J dbl(double x) { return {{"value", x}, {"bits", hexbits(x)}}; }

// Same admission conditions as Model::checked_mesh (field_model.hpp), then the sealed Mesh constructor.
axial::Mesh make_mesh(int nT, int nA, int nAir, double width = .009) {
    hc::need(nT >= 9 && nT <= 500 && nA >= 2 && nAir >= 1 && std::isfinite(width) && width > .0054, "invalid immutable mesh specification");
    const std::uint64_t nn = nT, z = std::uint64_t(nA) + 2 * std::uint64_t(nAir);
    hc::need(2 * nn * (nn + 1) * (z + 1) + (nn + 1) * (nn + 1) * z <= std::uint64_t(INT_MAX), "mesh index product overflow");
    return axial::Mesh(nT, nA, width, nAir);
}

// ------------------------------------------------------------------ device memory model
// Bytes per edge, face, node, cell and fixed bytes of one gpu::ResidentField (allocation list of
// resident_field.cpp build(); every buffer rounded up to >= 16 bytes, leaf tables included).
template <class MeshT> std::uint64_t resident_bytes(const MeshT& g, std::size_t metal_cells) {
    const std::uint64_t ne = g.ne, nb = g.nb, nn = g.nn, nc = g.cells.size(), nm = metal_cells;
    auto r16 = [](std::uint64_t b) { return std::max<std::uint64_t>(b, 16); };
    std::uint64_t t = 0;
    // static operator structure
    t += r16(4 * nb * 4) + r16(4 * nb * 8) + r16(nb * 8) + r16((ne + 1) * 4) + r16(4 * nb * 4) + r16(4 * nb * 8) + r16(4 * nb * 8);  // fe fc H ctp ctf ctc pc2
    t += r16(ne * 8) + r16(ne) + r16(ne * 8) + r16(ne * 8);   // sqrtM fixed M Y
    t += r16(ne * 4) * 2 + r16(ne * 8);                        // ea eb elen
    t += r16((ne + 1) * 4) + r16(12 * nc * 4) + r16(nc * 8) * 2 + r16(12 * nc * 4) + r16(nc * 4) + r16(nm * 4);  // ecp ecc qv cvol cedge body metal
    t += r16((nn + 1) * 4) + r16(2 * ne * 4) * 2;              // ntp nte nts
    t += 2 * r16(nn * 8);                                       // gauss increment accumulation
    t += r16(ne * 8) * 2 + r16(nn) * 2 + r16(ne * 8) + r16(nn * 8);   // feed bg masks region VN
    t += 2 * (2 * r16(ne * 8) + 2 * r16(nb * 8) + 8 * r16(nn * 8));   // slots
    t += r16(nc * 8) + 2 * r16(ne * 8) + r16(nn * 8);           // sigma plasma zero_ne zero_nn
    t += 10 * r16(nn * 8) + 14 * r16(ne * 8) + 7 * r16(nb * 8) + 10 * r16(nn * 8) + 3 * r16(ne * 8);
    t += r16(nc * 8) + r16(nm * 8) + r16(ne * 8 * 15) + r16(nb * 8 * 1) + r16(nn * 8 * 11) + r16(nc * 8 * 5);
    t += r16(nn * 4) + r16(nc * 4) + 2 * r16(ne * 4) + r16(nb * 4) + 2 * r16(nn * 4) + r16(nc * 4) + 2 * r16(std::max(ne, nb) * 4);
    t += 32 * 256 * 16 + 3 * 32 * 16;
    const std::uint64_t leaves = 6 * (ne * 8 / 2048 + 2) + 8 * (nn * 8 / 2048 + 2) + 16;
    t += leaves * 32 * 3 + 4096 * 4 + 4096 * 2048 / 16;
    t += 5 * r16(ne * 8) + r16(nb * 8) + 2 * 512 * 8 + 112 + 3 * 200 * 8 * 2 + 3 * 201 * 8;
    return t;
}
template <class MeshT> std::uint64_t resident_pinned_bytes(const MeshT& g, std::size_t metal_cells) {
    const std::uint64_t ne = g.ne, nb = g.nb, nn = g.nn, nc = g.cells.size();
    const std::uint64_t big = std::max({ne, nb, nn, nc}) * 16 * 4;
    const std::uint64_t leaves = 6 * (ne * 8 / 2048 + 2) + 8 * (nn * 8 / 2048 + 2) + 16;
    return std::max<std::uint64_t>(big, 12 * nc * 4 + 64) + 112 + 3 * 32 * 16 + 2 * ne * 8 + std::max<std::uint64_t>(metal_cells * 8, 16) +
           leaves * 32 * 3 + 4096 * 2048 / 16 + 4096 * 4 + (nc + ne + 10 * nn) * 8;
}

// ------------------------------------------------------------------ mesh report
J mesh_record(const axial::Mesh& g) {
    J j = J::object();
    std::size_t metal = 0, metal0 = 0, metal1 = 0, metal_nodes = 0, short_nodes = 0, fixedE = 0, outer = 0;
    long double vol = 0, metal_vol = 0, vn = 0, ve = 0, vb = 0;
    for (const auto& c : g.cells) { vol += c.volume; if (c.body >= 0) { ++metal; metal_vol += c.volume; (c.body ? metal1 : metal0)++; } }
    for (int k = 0; k < g.nn; ++k) { metal_nodes += g.node_body[k] >= 0; short_nodes += g.short_node[k]; outer += g.outer_node[k]; vn += g.VN[k]; }
    for (int e = 0; e < g.ne; ++e) { fixedE += g.fixed_E[e]; ve += g.VE[e]; }
    for (int b = 0; b < g.nb; ++b) vb += g.VB[b];
    const double box = g.side * g.side * (g.nz * g.dz);
    j["spec"] = {{"transverse", g.nx}, {"axial", g.nz_foil}, {"air", g.kshort}, {"width_m", g.side}};
    j["counts"] = {{"nodes", g.nn}, {"edges", g.ne}, {"faces", g.nb}, {"cells", g.cells.size()}, {"nx", g.nx}, {"ny", g.ny}, {"nz", g.nz},
                   {"metal_cells", metal}, {"metal_cells_lower", metal0}, {"metal_cells_upper", metal1}, {"metal_nodes", metal_nodes},
                   {"short_nodes", short_nodes}, {"pec_edges", fixedE}, {"outer_nodes", outer}};
    j["spacing"] = {{"dx_m", dbl(g.dx)}, {"dy_m", dbl(g.dy)}, {"dz_m", dbl(g.dz)}, {"zmin_m", dbl(g.zmin)}, {"zmax_m", dbl(g.zmin + g.nz * g.dz)},
                    {"aspect_dz_over_dx", g.dz / g.dx}};
    // resolution of the physical features
    j["cells_across"] = {{"foil_thickness", .0002 / g.dy}, {"gap", .005 / g.dy}, {"foil_width", .005 / g.dx}, {"foil_height", g.nz_foil},
                         {"air_per_end", g.kshort}};
    // feature positions as Hall computes them (node_position of the bounding node planes)
    auto ix = [&](double x) { return int(std::lround((x + g.side / 2) / g.dx)); };
    auto iy = [&](double y) { return int(std::lround((y + g.side / 2) / g.dy)); };
    J faces = J::object();
    for (auto [name, v] : std::vector<std::pair<std::string, double>>{{"y_lower_outer", -.0027}, {"y_lower_inner", -.0025}, {"y_upper_inner", .0025}, {"y_upper_outer", .0027}}) {
        const int j0 = iy(v); const double pos = g.node_position(g.node(0, j0, 0))[1];
        faces[name] = {{"index", j0}, {"position", dbl(pos)}, {"literal", dbl(v)}, {"position_minus_literal_m", pos - v},
                       {"alignment_cells", std::abs((v + g.side / 2) / g.dy - std::round((v + g.side / 2) / g.dy))}};
    }
    for (auto [name, v] : std::vector<std::pair<std::string, double>>{{"x_foil_left", -.0025}, {"x_foil_right", .0025}}) {
        const int i0 = ix(v); const double pos = g.node_position(g.node(i0, 0, 0))[0];
        faces[name] = {{"index", i0}, {"position", dbl(pos)}, {"literal", dbl(v)}, {"position_minus_literal_m", pos - v},
                       {"alignment_cells", std::abs((v + g.side / 2) / g.dx - std::round((v + g.side / 2) / g.dx))}};
    }
    const double zs = g.node_position(g.node(0, 0, g.kshort))[2], zf = g.node_position(g.node(0, 0, g.kfeed))[2];
    faces["z_short_plane"] = {{"index", g.kshort}, {"position", dbl(zs)}, {"literal", dbl(-g.depth / 2)}, {"position_minus_literal_m", zs + g.depth / 2}};
    faces["z_feed_port"] = {{"index", g.kfeed}, {"position", dbl(zf)}, {"literal", dbl(g.depth / 2)}, {"position_minus_literal_m", zf - g.depth / 2}};
    faces["outer_box"] = {{"x_min", dbl(g.node_position(0)[0])}, {"x_max", dbl(g.node_position(g.node(g.nx, 0, 0))[0])},
                          {"z_min", dbl(g.node_position(0)[2])}, {"z_max", dbl(g.node_position(g.node(0, 0, g.nz))[2])}};
    j["feature_positions"] = faces;
    // the classification that Hall derives from the positions (metal/short/PEC masks) checked against the
    // nominal geometry: every metal cell lies inside a foil, the foil cross-section is fully covered.
    const long double nominal_metal = 2.L * .005L * .0002L * g.depth;
    j["volumes"] = {{"cells_m3", double(vol)}, {"box_m3", box}, {"cells_rel_box", double((vol - box) / box)},
                    {"node_dual_rel_box", double((vn - box) / box)}, {"edge_dual_rel_3box", double((ve - 3 * box) / (3 * box))},
                    {"face_dual_rel_3box", double((vb - 3 * box) / (3 * box))},
                    {"metal_m3", double(metal_vol)}, {"metal_rel_nominal", double((metal_vol - nominal_metal) / nominal_metal)}};
    // feed contact quadrature (must normalize to 1 within 2e-14) and edge count
    {
        const auto f = g.feed_shape(); std::size_t nfe = 0; long double s = 0;
        for (double x : f) if (x != 0) { ++nfe; s += x; }
        j["feed"] = {{"edges", nfe}, {"sum_shape", double(s)}};
    }
    // operator quality: positivity, discrete exactness of curl(grad) and div(curl) on random fields
    {
        double minM = INFINITY, minH = INFINITY, minVN = INFINITY, minlen = INFINITY, maxlen = 0;
        for (int e = 0; e < g.ne; ++e) { minM = std::min(minM, g.M[e]); minlen = std::min(minlen, g.edges[e].length); maxlen = std::max(maxlen, g.edges[e].length); }
        for (int b = 0; b < g.nb; ++b) minH = std::min(minH, g.H[b]);
        for (int k = 0; k < g.nn; ++k) minVN = std::min(minVN, g.VN[k]);
        std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(-1, 1);
        axial::Vec phi(g.nn); for (auto& x : phi) x = U(rng);
        const auto cg = g.curl(g.gradient(phi)); double cgmax = 0, gscale = 0;
        for (double x : cg) cgmax = std::max(cgmax, std::abs(x));
        for (double x : g.gradient(phi)) gscale = std::max(gscale, std::abs(x));
        axial::Vec E(g.ne); for (auto& x : E) x = U(rng);
        const auto B = g.curl(E); const auto dB = g.magnetic_divergence(B); double dbmax = 0, bscale = 0;
        for (double x : dB) dbmax = std::max(dbmax, std::abs(x));
        for (double x : B) bscale = std::max(bscale, std::abs(x));
        j["operator"] = {{"min_M", minM}, {"min_H", minH}, {"min_VN", minVN}, {"min_edge_m", minlen}, {"max_edge_m", maxlen},
                         {"curl_grad_max_rel", cgmax / (gscale / g.dx)}, {"div_curl_max_rel", dbmax / (bscale / g.dx)}};
    }
    // incidence sanity: every face 4 edges, every cell 12 distinct edges
    {
        bool ok = true;
        for (const auto& c : g.C) ok = ok && c.count == 4;
        for (const auto& c : g.cells) { auto e = c.edges; std::sort(e.begin(), e.end()); ok = ok && std::adjacent_find(e.begin(), e.end()) == e.end(); }
        j["incidence_ok"] = ok;
    }
    // identities
    {
        hc::Digest d("resolution-series-mesh-v1"); hc::encode_mesh(d, g); j["encode_mesh_sha256"] = d.hex();
        J mj = g.json(); j["mesh_json"] = mj; j["mesh_json_sha256"] = hc::sha256(mj.dump());
        const auto o = mf::geometry(g);
        J k = {{"lower_foil", {o.lower_foil.left_m, o.lower_foil.right_m, o.lower_foil.bottom_m, o.lower_foil.top_m}},
               {"upper_foil", {o.upper_foil.left_m, o.upper_foil.right_m, o.upper_foil.bottom_m, o.upper_foil.top_m}},
               {"short_z_m", o.short_z_m}, {"feed_z_m", o.feed_z_m}, {"width_m", o.width_m}, {"depth_m", o.depth_m}};
        hc::Digest kd("resolution-series-kinetic-geometry-v1");
        for (double v : {o.lower_foil.left_m, o.lower_foil.right_m, o.lower_foil.bottom_m, o.lower_foil.top_m, o.upper_foil.left_m, o.upper_foil.right_m, o.upper_foil.bottom_m,
                         o.upper_foil.top_m, o.short_z_m, o.feed_z_m, o.width_m, o.gap_m, o.depth_m}) kd.f64(v);
        j["kinetic_geometry"] = k; j["kinetic_geometry_sha256"] = kd.hex();
    }
    j["device_memory_model"] = {{"resident_field_bytes", resident_bytes(g, metal)}, {"pinned_host_bytes", resident_pinned_bytes(g, metal)}};
    return j;
}
int mesh_report(const std::filesystem::path& out, const std::vector<std::array<int, 3>>& specs) {
    J rows = J::array();
    for (const auto& s : specs) {
        auto t = steady::now();
        auto g = std::make_unique<axial::Mesh>(make_mesh(s[0], s[1], s[2]));
        J r = mesh_record(*g);
        r["construct_ms"] = ms_since(t);
        rows.push_back(r);
        std::cout << "mesh " << s[0] << "x" << s[1] << "+" << s[2] << " edges " << r["counts"]["edges"] << std::endl;
    }
    // cross-mesh identity of every feature position (bitwise) and of the kinetic geometry
    J id = J::object();
    for (const char* k : {"y_lower_outer", "y_lower_inner", "y_upper_inner", "y_upper_outer", "x_foil_left", "x_foil_right", "z_short_plane", "z_feed_port"}) {
        J bits = J::array(); bool same = true;
        for (const auto& r : rows) { bits.push_back(r["feature_positions"][k]["position"]["bits"]); same = same && bits.back() == bits.front(); }
        id[k] = {{"bitwise_identical", same}, {"bits", bits}};
    }
    {
        J d = J::array(); bool same = true;
        for (const auto& r : rows) { d.push_back(r["kinetic_geometry_sha256"]); same = same && d.back() == d.front(); }
        id["kinetic_geometry"] = {{"bitwise_identical", same}, {"sha256", d}};
    }
    J summary{{"schema", "resolution-series-mesh-report-v1"}, {"meshes", rows}, {"cross_mesh_identity", id}};
    std::filesystem::create_directories(out);
    std::ofstream(out / "mesh-report.json") << summary.dump(1) << '\n';
    return 0;
}

// ------------------------------------------------------------------ instrumented DC seed (copy of axial::dc_initial)
struct SolveTrace { std::string name; int unknowns = 0; std::vector<int> cycle_iterations; std::vector<double> cycle_true_residual; double target = 0, ms = 0; bool converged = false; };
// Copy of foil::solve (foil_em.hpp, 8dddec0) with iteration counters; the arithmetic is unchanged.
foil::Vec traced_solve(const std::function<foil::Vec(const foil::Vec&)>& apply, const foil::Vec& b, const foil::Vec& diagonal, int maxit, SolveTrace& tr) {
    using foil::Vec; auto t0 = steady::now();
    Vec x(b.size()), r = b; const double target = std::max(1e-30, 2e-13 * foil::norm(b)); tr.target = target; tr.unknowns = int(b.size());
    for (int correction = 0; correction < 3; ++correction) {
        Vec z(r.size()); for (std::size_t i = 0; i < r.size(); ++i) z[i] = r[i] / diagonal[i]; Vec p = z; double rz = foil::dot(r, z);
        int iteration = 0;
        for (; iteration < maxit && foil::norm(r) > target * .2; ++iteration) {
            auto Ap = apply(p); double denominator = foil::dot(p, Ap); hall::require(denominator > 0 && std::isfinite(denominator), "non-SPD field operator");
            double alpha = rz / denominator;
            for (std::size_t i = 0; i < r.size(); ++i) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; z[i] = r[i] / diagonal[i]; }
            double nrz = foil::dot(r, z), beta = nrz / rz; for (std::size_t i = 0; i < r.size(); ++i) p[i] = z[i] + beta * p[i]; rz = nrz;
        }
        r = foil::subtract(b, apply(x)); tr.cycle_iterations.push_back(iteration); tr.cycle_true_residual.push_back(foil::norm(r));
        if (foil::norm(r) <= target) { tr.converged = true; tr.ms = ms_since(t0); return x; }
    }
    tr.ms = ms_since(t0);
    return x;   // not converged: recorded (the sealed solve throws here)
}
foil::Vec traced_scalar_dirichlet(const axial::Mesh& g, const foil::Vec& weights, const foil::Vec& rhs, const std::vector<bool>& fixed, SolveTrace& tr) {
    using foil::Vec; Vec diagonal(g.nn), b = rhs;
    for (int e = 0; e < g.ne; ++e) { double q = weights[e] / std::pow(g.edges[e].length, 2); diagonal[g.edges[e].a] += q; diagonal[g.edges[e].b] += q; }
    for (int k = 0; k < g.nn; ++k) if (fixed[k]) { diagonal[k] = 1; b[k] = 0; }
    auto apply = [&](Vec v) { for (int k = 0; k < g.nn; ++k) if (fixed[k]) v[k] = 0;
        auto d = g.gradient(v); for (int e = 0; e < g.ne; ++e) d[e] *= weights[e]; auto a = g.gradient_transpose(d);
        for (int k = 0; k < g.nn; ++k) if (fixed[k]) a[k] = 0; return a; };
    return traced_solve(apply, b, diagonal, 4000, tr);
}
// Returns the traced initial E, B (uniform-sigma path of the sealed dc_initial, the only path the Model uses).
J dc_profile(const axial::Mesh& g, const foil::Vec& sigma_cells, double current, foil::Vec& E_out, foil::Vec& B_out) {
    using foil::Vec; J j = J::object(); std::vector<SolveTrace> traces;
    auto feed = g.feed_shape(); auto S = g.conduction(sigma_cells);
    Vec impressed = feed; for (double& v : impressed) v *= current;
    auto rhs = g.gradient_transpose(impressed); std::vector<bool> fixed(g.nn);
    for (int k = 0; k < g.nn; ++k) fixed[k] = g.node_body[k] < 0 || g.short_node[k];
    double uniform_sigma = 0; for (std::size_t c = 0; c < g.cells.size(); ++c) if (g.cells[c].body >= 0) uniform_sigma = sigma_cells[c];
    Vec phi(g.nn);
    double gradient = current / (.005 * .0002 * uniform_sigma);
    for (int k = g.kshort; k <= g.kfeed; ++k) for (int jj = 0; jj <= g.ny; ++jj) for (int i = 0; i <= g.nx; ++i) {
        int q = g.node(i, jj, k); if (g.node_body[q] >= 0) phi[q] = (g.node_body[q] ? -1 : 1) * gradient * (k - g.kshort) * g.dz; }
    for (int k = 0; k < g.nn; ++k) fixed[k] = g.node_body[k] >= 0 || g.short_node[k];
    auto known = g.gradient(phi); for (int e = 0; e < g.ne; ++e) known[e] *= g.VE[e];
    rhs = g.gradient_transpose(known); for (double& v : rhs) v = -v;
    SolveTrace vac; vac.name = "dc_vacuum_scalar";
    Vec vacuum = traced_scalar_dirichlet(g, g.VE, rhs, fixed, vac); traces.push_back(vac);
    for (int k = 0; k < g.nn; ++k) phi[k] += vacuum[k];
    Vec Ein = g.gradient(phi); for (double& v : Ein) v = -v;
    Vec solid(g.ne), closed(g.ne);
    for (int e = 0; e < g.ne; ++e) closed[e] = (solid[e] = S[e] * Ein[e]) + impressed[e];
    Vec short_current(g.ne);
    for (int jj = 0; jj < g.ny; ++jj) for (int i = 0; i <= g.nx; ++i) short_current[g.ey(i, jj, g.kshort)] = -impressed[g.ey(i, jj, g.kfeed)];
    for (int e = 0; e < g.ne; ++e) closed[e] += short_current[e];
    Vec diagonal(g.ne);
    for (int b = 0; b < g.nb; ++b) for (int q = 0; q < 4; ++q) { int e = g.C[b].index[q]; diagonal[e] += g.H[b] * std::pow(g.C[b].coefficient[q], 2); }
    for (int e = 0; e < g.ne; ++e) { const auto& v = g.edges[e];
        diagonal[e] += g.VE[e] * g.VE[e] / (foil::mu * v.length * v.length) * ((g.outer_node[v.a] ? 0 : 1 / g.VN[v.a]) + (g.outer_node[v.b] ? 0 : 1 / g.VN[v.b]));
        if (g.fixed_A[e]) { diagonal[e] = 1; closed[e] = 0; } }
    auto apply = [&](Vec A) { for (int e = 0; e < g.ne; ++e) if (g.fixed_A[e]) A[e] = 0;
        auto K = g.curl_adjoint(g.curl(A)); for (int e = 0; e < g.ne; ++e) A[e] *= g.VE[e];
        auto div = g.gradient_transpose(A); for (int k = 0; k < g.nn; ++k) div[k] = g.outer_node[k] ? 0 : div[k] / (foil::mu * g.VN[k]); auto grad = g.gradient(div);
        for (int e = 0; e < g.ne; ++e) K[e] = g.fixed_A[e] ? 0 : K[e] + g.VE[e] * grad[e]; return K; };
    SolveTrace vec; vec.name = "dc_coulomb_vector";
    Vec A = traced_solve(apply, closed, diagonal, 5000, vec); traces.push_back(vec);
    E_out = Ein; B_out = g.curl(A);
    J t = J::array();
    for (const auto& s : traces) t.push_back({{"name", s.name}, {"unknowns", s.unknowns}, {"cycle_iterations", s.cycle_iterations},
                                              {"cycle_true_residual", s.cycle_true_residual}, {"target", s.target}, {"ms", s.ms}, {"converged", s.converged},
                                              {"maxit_per_cycle", s.name == "dc_vacuum_scalar" ? 4000 : 5000}});
    j["solves"] = t;
    return j;
}

// ------------------------------------------------------------------ DC seed attainable-accuracy diagnostic
// For the two DC seed solves: the sealed three correction cycles, then EXTRA further cycles of the same algorithm, and at
// the end of each phase the residual of x evaluated (a) with the sealed double apply, (b) with the same operator in long
// double, and (c) the componentwise rounding scale || |A| |x| ||_2 (eps || |A||x| || / ||b|| is the size of the evaluation
// noise of the double apply relative to ||b||, to compare with the 2e-13 target).
using LD = long double;
struct LdOps {
    const axial::Mesh& g;
    std::vector<LD> grad(const std::vector<LD>& p, bool ab) const { std::vector<LD> v(g.ne);
        for (int e = 0; e < g.ne; ++e) { const auto& E = g.edges[e]; v[e] = ab ? (std::abs(p[E.b]) + std::abs(p[E.a])) / E.length : (p[E.b] - p[E.a]) / E.length; } return v; }
    std::vector<LD> gradT(const std::vector<LD>& v, bool ab) const { std::vector<LD> p(g.nn);
        for (int e = 0; e < g.ne; ++e) { const auto& E = g.edges[e]; LD q = v[e] / E.length; if (ab) { q = std::abs(q); p[E.a] += q; p[E.b] += q; } else { p[E.a] -= q; p[E.b] += q; } } return p; }
    std::vector<LD> curl(const std::vector<LD>& A, bool ab) const { std::vector<LD> B(g.nb);
        for (int b = 0; b < g.nb; ++b) for (int k = 0; k < 4; ++k) { const LD c = g.C[b].coefficient[k], a = A[g.C[b].index[k]]; B[b] += ab ? std::abs(c) * std::abs(a) : c * a; } return B; }
    std::vector<LD> curlT(const std::vector<LD>& B, bool ab) const { std::vector<LD> E(g.ne);
        for (int b = 0; b < g.nb; ++b) for (int k = 0; k < 4; ++k) { const LD c = g.C[b].coefficient[k] * (LD)g.H[b]; E[g.C[b].index[k]] += ab ? std::abs(c) * std::abs(B[b]) : c * B[b]; } return E; }
};
std::vector<LD> to_ld(const foil::Vec& v) { return std::vector<LD>(v.begin(), v.end()); }
LD norm_ld(const std::vector<LD>& v) { LD s = 0; for (auto x : v) s += x * x; return std::sqrt(s); }
J accuracy(const std::function<foil::Vec(const foil::Vec&)>& apply, const std::function<std::vector<LD>(const std::vector<LD>&, bool)>& apply_ld,
           const foil::Vec& b, const foil::Vec& x) {
    const auto Ax = apply(x); const auto xl = to_ld(x); const auto Axl = apply_ld(xl, false); const auto Aabs = apply_ld(xl, true);
    std::vector<LD> rd(b.size()), rl(b.size()), diff(b.size());
    for (std::size_t i = 0; i < b.size(); ++i) { rd[i] = (LD)b[i] - Ax[i]; rl[i] = (LD)b[i] - Axl[i]; diff[i] = (LD)Ax[i] - Axl[i]; }
    const LD nb = norm_ld(to_ld(b));
    return {{"residual_double_apply_rel", double(norm_ld(rd) / nb)}, {"residual_long_double_apply_rel", double(norm_ld(rl) / nb)},
            {"apply_evaluation_error_rel", double(norm_ld(diff) / nb)}, {"abs_scale_rel", double(norm_ld(Aabs) / nb)},
            {"eps_times_abs_scale_rel", double(std::numeric_limits<double>::epsilon() / 2 * norm_ld(Aabs) / nb)}, {"target_rel", 2e-13}};
}
foil::Vec more_cycles(const std::function<foil::Vec(const foil::Vec&)>& apply, const foil::Vec& b, const foil::Vec& diagonal, int maxit, foil::Vec x, int extra, J& hist) {
    using foil::Vec; const double target = std::max(1e-30, 2e-13 * foil::norm(b));
    Vec r = foil::subtract(b, apply(x));
    for (int c = 0; c < extra; ++c) {
        Vec z(r.size()); for (std::size_t i = 0; i < r.size(); ++i) z[i] = r[i] / diagonal[i]; Vec pp = z; double rz = foil::dot(r, z); int it = 0;
        for (; it < maxit && foil::norm(r) > target * .2; ++it) {
            auto Ap = apply(pp); double den = foil::dot(pp, Ap); double alpha = rz / den;
            for (std::size_t i = 0; i < r.size(); ++i) { x[i] += alpha * pp[i]; r[i] -= alpha * Ap[i]; z[i] = r[i] / diagonal[i]; }
            double nrz = foil::dot(r, z), beta = nrz / rz; for (std::size_t i = 0; i < r.size(); ++i) pp[i] = z[i] + beta * pp[i]; rz = nrz;
        }
        r = foil::subtract(b, apply(x)); hist.push_back({{"cycle", 4 + c}, {"iterations", it}, {"true_residual_rel", foil::norm(r) / foil::norm(b)}});
    }
    return x;
}
int dc_floor(const std::filesystem::path& out, const std::string& trace, int nT, int nA, int nAir, int extra) {
    using foil::Vec; std::filesystem::create_directories(out);
    auto g = make_mesh(nT, nA, nAir); LdOps L{g};
    auto law = electrode_heat_v2::conditional_ss304_v2(ss304v2::Qualification::evaluated_solid_1500,
        ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation, 8000., electrode_heat_v2::material_headers_sha256);
    Vec sigma(g.cells.size()); for (std::size_t c = 0; c < sigma.size(); ++c) if (g.cells[c].body >= 0) sigma[c] = law.sigma(300.);
    hc::CurrentDrive drive(trace); const double current = drive.current(drive.start());
    J summary{{"schema", "resolution-series-dc-floor-v1"}, {"spec", {nT, nA, nAir, .009}}, {"extra_cycles", extra}};
    // scalar (vacuum potential) solve, as in dc_initial
    auto feed = g.feed_shape(); auto S = g.conduction(sigma); Vec impressed = feed; for (double& v : impressed) v *= current;
    std::vector<bool> fixed(g.nn); for (int k = 0; k < g.nn; ++k) fixed[k] = g.node_body[k] >= 0 || g.short_node[k];
    double us = 0; for (std::size_t c = 0; c < g.cells.size(); ++c) if (g.cells[c].body >= 0) us = sigma[c];
    Vec phi(g.nn); const double grad = current / (.005 * .0002 * us);
    for (int k = g.kshort; k <= g.kfeed; ++k) for (int j = 0; j <= g.ny; ++j) for (int i = 0; i <= g.nx; ++i) { int q = g.node(i, j, k); if (g.node_body[q] >= 0) phi[q] = (g.node_body[q] ? -1 : 1) * grad * (k - g.kshort) * g.dz; }
    auto known = g.gradient(phi); for (int e = 0; e < g.ne; ++e) known[e] *= g.VE[e]; Vec rhs = g.gradient_transpose(known); for (double& v : rhs) v = -v;
    {
        Vec diagonal(g.nn), b = rhs;
        for (int e = 0; e < g.ne; ++e) { double q = g.VE[e] / std::pow(g.edges[e].length, 2); diagonal[g.edges[e].a] += q; diagonal[g.edges[e].b] += q; }
        for (int k = 0; k < g.nn; ++k) if (fixed[k]) { diagonal[k] = 1; b[k] = 0; }
        auto apply = [&](Vec v) { for (int k = 0; k < g.nn; ++k) if (fixed[k]) v[k] = 0; auto d = g.gradient(v); for (int e = 0; e < g.ne; ++e) d[e] *= g.VE[e];
            auto a = g.gradient_transpose(d); for (int k = 0; k < g.nn; ++k) if (fixed[k]) a[k] = 0; return a; };
        auto apply_ld = [&](std::vector<LD> v, bool ab) { for (int k = 0; k < g.nn; ++k) if (fixed[k]) v[k] = 0; auto d = L.grad(v, ab); for (int e = 0; e < g.ne; ++e) d[e] *= g.VE[e];
            auto a = L.gradT(d, ab); for (int k = 0; k < g.nn; ++k) if (fixed[k]) a[k] = 0; return a; };
        SolveTrace tr; auto t = steady::now(); Vec x = traced_solve(apply, b, diagonal, 4000, tr);
        J s{{"cycle_iterations", tr.cycle_iterations}, {"cycle_true_residual_rel", J::array()}, {"converged_sealed", tr.converged}, {"ms", ms_since(t)}};
        for (double r : tr.cycle_true_residual) s["cycle_true_residual_rel"].push_back(r / foil::norm(b));
        s["accuracy_after_sealed_cycles"] = accuracy(apply, apply_ld, b, x);
        summary["dc_vacuum_scalar"] = s;
        for (int k = 0; k < g.nn; ++k) phi[k] += x[k];
    }
    std::ofstream(out / "dc-floor.json") << summary.dump(1) << '\n';
    // vector (Coulomb gauge) solve, as in dc_initial (uniform-DC path)
    {
        Vec Ein = g.gradient(phi); for (double& v : Ein) v = -v; Vec closed(g.ne);
        for (int e = 0; e < g.ne; ++e) closed[e] = S[e] * Ein[e] + impressed[e];
        Vec shortc(g.ne); for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i) shortc[g.ey(i, j, g.kshort)] = -impressed[g.ey(i, j, g.kfeed)];
        for (int e = 0; e < g.ne; ++e) closed[e] += shortc[e];
        Vec diagonal(g.ne);
        for (int b = 0; b < g.nb; ++b) for (int q = 0; q < 4; ++q) { int e = g.C[b].index[q]; diagonal[e] += g.H[b] * std::pow(g.C[b].coefficient[q], 2); }
        for (int e = 0; e < g.ne; ++e) { const auto& v = g.edges[e];
            diagonal[e] += g.VE[e] * g.VE[e] / (foil::mu * v.length * v.length) * ((g.outer_node[v.a] ? 0 : 1 / g.VN[v.a]) + (g.outer_node[v.b] ? 0 : 1 / g.VN[v.b]));
            if (g.fixed_A[e]) { diagonal[e] = 1; closed[e] = 0; } }
        auto apply = [&](Vec A) { for (int e = 0; e < g.ne; ++e) if (g.fixed_A[e]) A[e] = 0;
            auto K = g.curl_adjoint(g.curl(A)); for (int e = 0; e < g.ne; ++e) A[e] *= g.VE[e];
            auto div = g.gradient_transpose(A); for (int k = 0; k < g.nn; ++k) div[k] = g.outer_node[k] ? 0 : div[k] / (foil::mu * g.VN[k]); auto gr = g.gradient(div);
            for (int e = 0; e < g.ne; ++e) K[e] = g.fixed_A[e] ? 0 : K[e] + g.VE[e] * gr[e]; return K; };
        auto apply_ld = [&](std::vector<LD> A, bool ab) { for (int e = 0; e < g.ne; ++e) if (g.fixed_A[e]) A[e] = 0;
            auto K = L.curlT(L.curl(A, ab), ab); for (int e = 0; e < g.ne; ++e) A[e] *= g.VE[e];
            auto div = L.gradT(A, ab); for (int k = 0; k < g.nn; ++k) div[k] = g.outer_node[k] ? 0 : div[k] / ((LD)foil::mu * g.VN[k]); auto gr = L.grad(div, ab);
            for (int e = 0; e < g.ne; ++e) K[e] = g.fixed_A[e] ? 0 : K[e] + (LD)g.VE[e] * gr[e]; return K; };
        SolveTrace tr; auto t = steady::now(); Vec x = traced_solve(apply, closed, diagonal, 5000, tr);
        J s{{"cycle_iterations", tr.cycle_iterations}, {"cycle_true_residual_rel", J::array()}, {"converged_sealed", tr.converged}, {"ms", ms_since(t)}};
        for (double r : tr.cycle_true_residual) s["cycle_true_residual_rel"].push_back(r / foil::norm(closed));
        s["accuracy_after_sealed_cycles"] = accuracy(apply, apply_ld, closed, x);
        J hist = J::array(); x = more_cycles(apply, closed, diagonal, 5000, x, extra, hist);
        s["extra_cycles"] = hist; s["accuracy_after_extra_cycles"] = accuracy(apply, apply_ld, closed, x);
        summary["dc_coulomb_vector"] = s;
    }
    summary["max_rss_kB"] = rss_kb("VmHWM:");
    std::ofstream(out / "dc-floor.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return 0;
}

int admit(const std::filesystem::path& out, const std::string& trace, int nT, int nA, int nAir, double width = .009) {
    std::filesystem::create_directories(out);
    J summary{{"schema", "resolution-series-admission-v1"}, {"spec", {nT, nA, nAir, width}}};
    auto g = make_mesh(nT, nA, nAir, width);
    // instrumented DC seed (same arithmetic as the sealed uniform path)
    auto law = electrode_heat_v2::conditional_ss304_v2(ss304v2::Qualification::evaluated_solid_1500,
        ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation, 8000., electrode_heat_v2::material_headers_sha256);
    foil::Vec sigma(g.cells.size()); for (std::size_t c = 0; c < sigma.size(); ++c) if (g.cells[c].body >= 0) sigma[c] = law.sigma(300.);
    hc::CurrentDrive drive(trace);
    foil::Vec E, B;
    auto t = steady::now();
    summary["dc_profile"] = dc_profile(g, sigma, drive.current(drive.start()), E, B);
    summary["dc_profile_ms"] = ms_since(t);
    std::ofstream(out / "admission.json") << summary.dump(1) << '\n';
    // the sealed Model (admission of record)
    t = steady::now();
    std::string error;
    try {
        hc::Model m(hc::MeshSpec{nT, nA, nAir, width}, trace);
        summary["model_construct_ms"] = ms_since(t);
        summary["admitted"] = true;
        summary["configuration_sha256"] = m.configuration_sha256();
        summary["motion_mesh_sha256"] = m.motion_mesh_sha256(mf::geometry(m.mesh()));
        summary["initial_gauss_solve_bound_C"] = m.initial_gauss_solve_bound_C();
        const auto& in = m.initial_fields();
        summary["dc"] = {{"terminal_V", in.terminal_V}, {"current_A", in.current_A}, {"dc_continuity_A", in.dc_continuity_A},
                         {"magnetic_residual_A_m", in.magnetic_residual_A_m}};
        summary["profile_equals_model_bitwise"] = {{"E", E == in.E}, {"B", B == in.B}};
        double bmax = 0; for (double x : in.B) bmax = std::max(bmax, std::abs(x));
        summary["dc_face_B_max_T"] = bmax;
        {   // DC magnetostatic port inductance L = 2 W_B / I^2 (W_B = sum_b H_b B_b^2 / 2 of the DC seed)
            long double W = 0; const auto& gg = m.mesh();
            for (int b = 0; b < gg.nb; ++b) W += .5L * gg.H[b] * (long double)in.B[b] * in.B[b];
            summary["dc_magnetic_energy_J"] = double(W);
            summary["dc_port_inductance_H"] = double(2 * W / ((long double)in.current_A * in.current_A));
            summary["dc_resistance_ohm"] = in.terminal_V / in.current_A;
        }
    } catch (const std::exception& e) {
        summary["model_construct_ms"] = ms_since(t); summary["admitted"] = false; summary["refusal"] = e.what();
    }
    summary["max_rss_kB"] = rss_kb("VmHWM:");
    std::ofstream(out / "admission.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return 0;
}

// ------------------------------------------------------------------ graded meshes (sizes, DC port inductance)
// The uniform-sigma path of axial::dc_initial (8dddec0) with the operators of graded::Mesh, the same operation order
// (analytic conductor potential, vacuum scalar Dirichlet solve, reflected short current, Coulomb-gauge vector potential),
// returning B; L = 2 W_B / I^2.
template <class G> foil::Vec dc_B(const G& g, const foil::Vec& sigma_cells, double current, J& info) {
    using foil::Vec;
    auto feed = g.feed_shape(); auto S = g.conduction(sigma_cells);
    Vec impressed = feed; for (double& v : impressed) v *= current;
    std::vector<bool> fixed(g.nn);
    double uniform_sigma = 0; for (std::size_t c = 0; c < g.cells.size(); ++c) if (g.cells[c].body >= 0) uniform_sigma = sigma_cells[c];
    Vec phi(g.nn); const double gradient = current / (.005 * .0002 * uniform_sigma);
    for (int k = g.kshort; k <= g.kfeed; ++k) for (int j = 0; j <= g.ny; ++j) for (int i = 0; i <= g.nx; ++i) {
        int q = g.node(i, j, k); if (g.node_body[q] >= 0) phi[q] = (g.node_body[q] ? -1 : 1) * gradient * (k - g.kshort) * g.dz; }
    for (int k = 0; k < g.nn; ++k) fixed[k] = g.node_body[k] >= 0 || g.short_node[k];
    auto known = g.gradient(phi); for (int e = 0; e < g.ne; ++e) known[e] *= g.VE[e];
    Vec rhs = g.gradient_transpose(known); for (double& v : rhs) v = -v;
    Vec vacuum;
    {   // scalar_dirichlet(g, VE, rhs, fixed), traced
        Vec diagonal(g.nn), b = rhs;
        for (int e = 0; e < g.ne; ++e) { double q = g.VE[e] / std::pow(g.edges[e].length, 2); diagonal[g.edges[e].a] += q; diagonal[g.edges[e].b] += q; }
        for (int k = 0; k < g.nn; ++k) if (fixed[k]) { diagonal[k] = 1; b[k] = 0; }
        auto apply = [&](Vec v) { for (int k = 0; k < g.nn; ++k) if (fixed[k]) v[k] = 0; auto d = g.gradient(v); for (int e = 0; e < g.ne; ++e) d[e] *= g.VE[e];
            auto a = g.gradient_transpose(d); for (int k = 0; k < g.nn; ++k) if (fixed[k]) a[k] = 0; return a; };
        SolveTrace tr; vacuum = traced_solve(apply, b, diagonal, 4000, tr);
        info["vacuum_scalar"] = {{"cycle_iterations", tr.cycle_iterations}, {"converged", tr.converged}, {"ms", tr.ms}};
    }
    for (int k = 0; k < g.nn; ++k) phi[k] += vacuum[k];
    Vec Ein = g.gradient(phi); for (double& v : Ein) v = -v;
    Vec closed(g.ne); for (int e = 0; e < g.ne; ++e) closed[e] = S[e] * Ein[e] + impressed[e];
    Vec shortc(g.ne); for (int j = 0; j < g.ny; ++j) for (int i = 0; i <= g.nx; ++i) shortc[g.ey(i, j, g.kshort)] = -impressed[g.ey(i, j, g.kfeed)];
    for (int e = 0; e < g.ne; ++e) closed[e] += shortc[e];
    Vec diagonal(g.ne);
    for (int b = 0; b < g.nb; ++b) for (int q = 0; q < 4; ++q) { int e = g.C[b].index[q]; diagonal[e] += g.H[b] * std::pow(g.C[b].coefficient[q], 2); }
    for (int e = 0; e < g.ne; ++e) { const auto& v = g.edges[e];
        diagonal[e] += g.VE[e] * g.VE[e] / (foil::mu * v.length * v.length) * ((g.outer_node[v.a] ? 0 : 1 / g.VN[v.a]) + (g.outer_node[v.b] ? 0 : 1 / g.VN[v.b]));
        if (g.fixed_A[e]) { diagonal[e] = 1; closed[e] = 0; } }
    auto apply = [&](Vec A) { for (int e = 0; e < g.ne; ++e) if (g.fixed_A[e]) A[e] = 0;
        auto K = g.curl_adjoint(g.curl(A)); for (int e = 0; e < g.ne; ++e) A[e] *= g.VE[e];
        auto div = g.gradient_transpose(A); for (int k = 0; k < g.nn; ++k) div[k] = g.outer_node[k] ? 0 : div[k] / (foil::mu * g.VN[k]); auto gr = g.gradient(div);
        for (int e = 0; e < g.ne; ++e) K[e] = g.fixed_A[e] ? 0 : K[e] + g.VE[e] * gr[e]; return K; };
    SolveTrace tr; Vec A = traced_solve(apply, closed, diagonal, 5000, tr);
    info["coulomb_vector"] = {{"cycle_iterations", tr.cycle_iterations}, {"converged", tr.converged}, {"ms", tr.ms},
                              {"final_true_residual_rel", tr.cycle_true_residual.empty() ? 0. : tr.cycle_true_residual.back() / foil::norm(closed)}};
    return g.curl(A);
}
template <class G> double inductance(const G& g, const foil::Vec& B, double I) {
    long double W = 0; for (int b = 0; b < g.nb; ++b) W += .5L * g.H[b] * (long double)B[b] * B[b];
    return double(2 * W / ((long double)I * I));
}
J graded_record(const graded::Mesh& g) {
    std::size_t metal = 0; for (const auto& c : g.cells) metal += c.body >= 0;
    double hmax = *std::max_element(g.hx.begin(), g.hx.end());
    double qmax = 1; for (std::size_t i = 1; i < g.hx.size(); ++i) qmax = std::max(qmax, std::max(g.hx[i] / g.hx[i - 1], g.hx[i - 1] / g.hx[i]));
    return {{"nx", g.nx}, {"nz", g.nz}, {"nodes", g.nn}, {"edges", g.ne}, {"faces", g.nb}, {"cells", g.cells.size()}, {"metal_cells", metal},
            {"box_m", g.side}, {"core_h_m", g.h}, {"largest_cell_m", hmax}, {"largest_adjacent_ratio", qmax},
            {"resident_field_bytes", resident_bytes(g, metal)}, {"pinned_host_bytes", resident_pinned_bytes(g, metal)}};
}
// graded-dc OUT TRACE H CORE BOX Q HMAX: sizes and the DC port inductance on one graded mesh; with BOX = CORE also the
// sealed dc_initial on axial::Mesh(n, 4, 2 CORE, 1) and the comparison of B bit for bit.
int graded_dc(const std::filesystem::path& out, const std::string& trace, double h, double core, double box, double q, double hmax) {
    std::filesystem::create_directories(out);
    auto t = steady::now();
    graded::Mesh g(graded::axis(h, core, box, q, hmax), 4, 1);
    J summary{{"schema", "resolution-series-graded-dc-v1"}, {"h_m", h}, {"core_half_m", core}, {"box_half_m", box}, {"q", q}, {"hmax_m", hmax},
              {"mesh", graded_record(g)}, {"construct_ms", ms_since(t)}};
    auto law = electrode_heat_v2::conditional_ss304_v2(ss304v2::Qualification::evaluated_solid_1500,
        ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation, 8000., electrode_heat_v2::material_headers_sha256);
    foil::Vec sigma(g.cells.size()); for (std::size_t c = 0; c < sigma.size(); ++c) if (g.cells[c].body >= 0) sigma[c] = law.sigma(300.);
    hc::CurrentDrive drive(trace); const double I = drive.current(drive.start());
    J info = J::object(); t = steady::now();
    const auto B = dc_B(g, sigma, I, info);
    summary["dc_ms"] = ms_since(t); summary["dc"] = info; summary["dc_port_inductance_H"] = inductance(g, B, I); summary["current_A"] = I;
    if (box == core) {   // reference: the sealed path on the uniform Hall mesh
        const int n = int(std::lround(2 * core / h));
        axial::Mesh u(n, 4, 2 * core, 1);
        foil::Vec su(u.cells.size()); for (std::size_t c = 0; c < su.size(); ++c) if (u.cells[c].body >= 0) su[c] = law.sigma(300.);
        const auto ini = axial::dc_initial(u, su, I);
        bool same_ops = u.ne == g.ne && u.nb == g.nb && u.M == g.M && u.H == g.H && u.VE == g.VE && u.VN == g.VN && u.Y == g.Y;
        for (int b = 0; b < u.nb && same_ops; ++b) for (int k = 0; k < 4; ++k) same_ops = same_ops && u.C[b].index[k] == g.C[b].index[k] && u.C[b].coefficient[k] == g.C[b].coefficient[k];
        summary["uniform_reference"] = {{"operators_bitwise_equal", same_ops}, {"feed_bitwise_equal", u.feed_shape() == g.feed_shape()},
                                        {"B_bitwise_equal", ini.B == B}, {"dc_port_inductance_H", inductance(u, ini.B, I)}};
    }
    summary["max_rss_kB"] = rss_kb("VmHWM:");
    std::ofstream(out / "graded-dc.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return 0;
}
// graded-sizes OUT H1,H2,.. BOX1,BOX2,.. Q1,Q2,.. HMAX: counts and resident memory of graded meshes (core = 4.5 mm)
int graded_sizes(const std::filesystem::path& out, std::vector<double> hs, std::vector<double> boxes, std::vector<double> qs, double hmax) {
    std::filesystem::create_directories(out);
    J rows = J::array();
    for (double h : hs) for (double box : boxes) for (double q : qs) {
        graded::Mesh g(graded::axis(h, .0045, box, q, hmax), 4, 1);
        J r = graded_record(g); r["h_m"] = h; r["box_half_m"] = box; r["q"] = q; r["hmax_m"] = hmax;
        rows.push_back(r);
        std::cout << "h " << h << " box " << 2 * box << " q " << q << " n " << g.nx << " edges " << g.ne << " bytes " << r["resident_field_bytes"] << std::endl;
    }
    std::ofstream(out / "graded-sizes.json") << J{{"schema", "resolution-series-graded-sizes-v1"}, {"rows", rows}}.dump(1) << '\n';
    return 0;
}

int cpu_vacuum(const std::filesystem::path& out, const std::string& trace, int nT, int nA, int nAir, long steps, double dt) {
    std::filesystem::create_directories(out);
    auto t = steady::now();
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{nT, nA, nAir, .009}, trace, hc::make_cpu_pcg(), hc::CarrierLimits{}, hc::ThermalDeclaration{},
                                                   std::max(dt, 1e-14));
    const double construct = ms_since(t);
    const auto& g = model->mesh(); const auto o = mf::geometry(g); const auto seed = model->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    hc::Owner owner(model, seed, declared(*model, seed, hc::charge_of_census(g, o, none)), o, none);
    Stats T; double max_resid = 0;
    for (long step = 1; step <= steps; ++step) {
        t = steady::now();
        auto p = owner.prepare_vacuum(dt);
        const double prep = ms_since(t); t = steady::now();
        Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
        hc::need(owner.preflight_vacuum_final(p, b), "vacuum preflight");
        const double pre = ms_since(t);
        const auto& sr = p.solve_record();
        T.add("prepare_ms", prep); T.add("preflight_ms", pre); T.add("iterations", sr.iterations); T.add("solve_ms", sr.solve_s * 1e3); T.add("check_ms", sr.check_s * 1e3);
        max_resid = std::max(max_resid, sr.independent_residual / sr.target);
        t = steady::now(); owner.publish_vacuum_prevalidated(std::move(p)); T.add("publish_ms", ms_since(t));
        T.add("step_ms", prep + pre + T.m["publish_ms"].v.back());
    }
    J summary{{"schema", "resolution-series-cpu-vacuum-v1"}, {"spec", {nT, nA, nAir}}, {"dt_s", dt}, {"steps", steps}, {"model_construct_ms", construct},
              {"phases", T.json()}, {"max_independent_residual_over_target", max_resid}, {"final_field_sha256", owner.field_sha256()},
              {"max_rss_kB", rss_kb("VmHWM:")}};
    std::ofstream(out / "cpu-vacuum.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return 0;
}

// ------------------------------------------------------------------ GPU
std::string probe_descriptor(const std::string& trace, const std::string& nvrtc, const std::string& uuid) {
    auto probe_model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    auto dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    const auto d = hc::resident_descriptor_sha256(dev, *probe_model);
    dev->close_json();
    return d;
}
// Per-node absolute Gauss ratio |defect|/bound of the accepted slot (device admission values).
J gauss_ratio(hc::ResidentOwner& res, std::vector<double>* ratio_out = nullptr) {
    hc::Vec d, b; res.device_field().gauss_profile(0, res.accepted_meta(), d, b);
    double mr = 0, md = 0, mb = 0; int node = -1;
    for (std::size_t k = 0; k < d.size(); ++k) { const double r = d[k] / b[k]; if (r > mr) { mr = r; node = int(k); } md = std::max(md, d[k]); mb = std::max(mb, b[k]); }
    if (ratio_out) { ratio_out->resize(d.size()); for (std::size_t k = 0; k < d.size(); ++k) (*ratio_out)[k] = d[k] / b[k]; }
    return {{"max_node_ratio", mr}, {"worst_node", node}, {"max_defect_C", md}, {"max_bound_C", mb},
            {"worst_defect_C", node >= 0 ? d[node] : 0.}, {"worst_bound_C", node >= 0 ? b[node] : 0.}};
}
// v7-style |B| bound: per cell, per component face maxima over the 3x3x3 cell neighbourhood.
J b_bounds(const axial::Mesh& g, const axial::State& s) {
    const int nx = g.nx, ny = g.ny, nz = g.nz; std::vector<double> mx(std::size_t(nx) * ny * nz), my(mx.size()), mz(mx.size());
    auto cid = [&](int i, int j, int k) { return (std::size_t(k) * ny + j) * nx + i; };
    double face_max = 0;
    for (int b = 0; b < g.nb; ++b) face_max = std::max(face_max, std::abs(s.B[b] + s.B_low[b]));
    double cell_centre = 0;
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        const auto c = cid(i, j, k);
        mx[c] = std::max(std::abs(s.B[g.bx(i, j, k)]), std::abs(s.B[g.bx(i + 1, j, k)]));
        my[c] = std::max(std::abs(s.B[g.by(i, j, k)]), std::abs(s.B[g.by(i, j + 1, k)]));
        mz[c] = std::max(std::abs(s.B[g.bz(i, j, k)]), std::abs(s.B[g.bz(i, j, k + 1)]));
        const double bx = .5 * (s.B[g.bx(i, j, k)] + s.B[g.bx(i + 1, j, k)]), by = .5 * (s.B[g.by(i, j, k)] + s.B[g.by(i, j + 1, k)]),
                     bz = .5 * (s.B[g.bz(i, j, k)] + s.B[g.bz(i, j, k + 1)]);
        cell_centre = std::max(cell_centre, std::sqrt(bx * bx + by * by + bz * bz));
    }
    double nb_max = 0; int wi = 0, wj = 0, wk = 0;
    for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
        double ax = 0, ay = 0, az = 0;
        for (int dk = -1; dk <= 1; ++dk) for (int dj = -1; dj <= 1; ++dj) for (int di = -1; di <= 1; ++di) {
            const int a = i + di, b = j + dj, c = k + dk; if (a < 0 || b < 0 || c < 0 || a >= nx || b >= ny || c >= nz) continue;
            const auto q = cid(a, b, c); ax = std::max(ax, mx[q]); ay = std::max(ay, my[q]); az = std::max(az, mz[q]); }
        const double v = std::sqrt(ax * ax + ay * ay + az * az);
        if (v > nb_max) { nb_max = v; wi = i; wj = j; wk = k; }
    }
    return {{"face_max_T", face_max}, {"cell_centre_max_T", cell_centre}, {"neighbourhood_bound_max_T", nb_max},
            {"neighbourhood_bound_cell", {wi, wj, wk}},
            {"neighbourhood_bound_cell_position_m", {(wi + .5) * g.dx - g.side / 2, (wj + .5) * g.dy - g.side / 2, g.zmin + (wk + .5) * g.dz}}};
}

int gpu_vacuum(const std::filesystem::path& out, const std::string& trace, const std::string& nvrtc, const std::string& uuid, int nT, int nA, int nAir,
               long steps, double dt, double dtmax, long every, double width = .009) {
    std::filesystem::create_directories(out);
    J summary{{"schema", "resolution-series-gpu-vacuum-v1"}, {"spec", {nT, nA, nAir, width}}, {"dt_s", dt}, {"dt_max_s", dtmax}, {"steps_requested", steps}};
    auto t = steady::now();
    auto cpu_model = std::make_shared<const hc::Model>(hc::MeshSpec{nT, nA, nAir, width}, trace);   // DC seed and ledger recipe (cpu-pcg identity)
    summary["cpu_model_construct_ms"] = ms_since(t);
    const auto& g = cpu_model->mesh(); const auto o = mf::geometry(g); const auto seed = cpu_model->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    const auto q0 = declared(*cpu_model, seed, hc::charge_of_census(g, o, none));
    summary["rss_after_cpu_model_kB"] = rss_kb("VmRSS:");
    const auto t_gpu = steady::now();
    t = steady::now();
    const auto desc = probe_descriptor(trace, nvrtc, uuid);
    summary["probe_descriptor_ms"] = ms_since(t);
    t = steady::now();
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{nT, nA, nAir, width}, trace, std::make_shared<hc::ResidentBackend>(desc), hc::CarrierLimits{},
                                                   hc::ThermalDeclaration{}, dtmax);
    summary["resident_model_construct_ms"] = ms_since(t);
    summary["configuration_sha256"] = model->configuration_sha256();
    auto dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    t = steady::now();
    gp::ResidentOptions opt; opt.sample_seed = 1;
    std::unique_ptr<hc::ResidentOwner> res;
    try {
        res = std::make_unique<hc::ResidentOwner>(model, dev, seed, q0, o, none, opt);
    } catch (const std::exception& e) {
        summary["owner_refusal"] = e.what(); summary["memory"] = J::parse(dev->memory_json());
        summary["teardown"] = J::parse(dev->close_json());
        std::ofstream(out / "gpu-vacuum.json") << summary.dump(1) << '\n';
        std::cout << summary.dump(1) << std::endl;
        return 2;
    }
    summary["owner_construct_ms"] = ms_since(t);
    summary["memory_after_owner"] = J::parse(dev->memory_json());
    summary["gauss_initial"] = gauss_ratio(*res);
    std::ofstream ledger(out / "gpu-vacuum-ledger.jsonl");
    Stats T; long done = 0; std::string refusal;
    double max_resid = 0, max_gauss_over = 0, max_divB_over = 0, max_energy = 0, max_heat = 0, max_bmax = 0, max_omega_c = 0, max_third = 0, max_incr = 0;
    J gauss_series = J::array();
    gauss_series.push_back({{"step", 0}, {"time_s", res->time_s()}, {"gauss", summary["gauss_initial"]}});
    const auto t_all = steady::now();
    for (long step = 1; step <= steps; ++step) {
        try {
            t = steady::now();
            auto p = res->prepare_vacuum(dt);
            const double prep = ms_since(t); t = steady::now();
            Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
            const bool ok = res->preflight_vacuum_final(p, b);
            const double pre = ms_since(t);
            hc::need(ok, "vacuum preflight");
            const auto& r = p.field_receipt(); const auto& c = p.candidate_constraints(); const auto& sr = p.solve_record(); const auto& cr = p.control();
            max_resid = std::max(max_resid, sr.independent_residual / sr.target);
            max_gauss_over = std::max(max_gauss_over, c.max_gauss_C / std::max(c.max_gauss_bound_C, 1e-300));
            max_divB_over = std::max(max_divB_over, c.max_divB_T_m / std::max(c.max_divB_bound_T_m, 1e-300));
            const double activity = std::abs(r.field_delta_J) + std::abs(r.plasma_work_J) + std::abs(r.joule_J) + std::abs(r.outer_loss_J) + std::abs(r.outer_supply_J) + std::abs(r.circuit_work_J);
            max_energy = std::max(max_energy, std::abs(r.balance_J) / (1e-20 + 1e-10 * activity));
            max_heat = std::max(max_heat, std::abs(p.field_heat_balance_J()));
            max_bmax = std::max(max_bmax, cr.B_max_T); max_omega_c = std::max(max_omega_c, cr.omega_c_dt_max);
            max_third = std::max(max_third, cr.third_difference_rel); max_incr = std::max(max_incr, cr.increment_rel);
            for (const auto& [k, v] : p.step().host_ms) T.add("prepare." + k, v);
            T.add("step.prepare_vacuum_ms", prep); T.add("step.preflight_vacuum_final_ms", pre);
            T.add("solve.iterations", sr.iterations); T.add("solve.corrections", sr.corrections); T.add("solve.residual_over_target", sr.independent_residual / sr.target);
            T.add("solve.rhs_norm", sr.rhs_norm);
            if (step == 1 || step % every == 0 || step == steps)
                ledger << J{{"step", step}, {"time_s", p.end_s()}, {"current_A", r.current_A}, {"terminal_V", r.terminal_V}, {"iterations", sr.iterations},
                            {"rhs_norm", sr.rhs_norm}, {"independent_residual", sr.independent_residual}, {"target", sr.target},
                            {"max_gauss_C", c.max_gauss_C}, {"max_gauss_bound_C", c.max_gauss_bound_C}, {"max_divB_T_m", c.max_divB_T_m},
                            {"B_max_T", cr.B_max_T}, {"omega_c_dt", cr.omega_c_dt_max}, {"increment_rel", cr.increment_rel},
                            {"third_difference_rel", cr.third_difference_rel}, {"prepare_ms", prep}, {"preflight_ms", pre}}.dump() << '\n';
            t = steady::now();
            res->publish_vacuum_prevalidated(std::move(p));
            T.add("step.publish_ms", ms_since(t));
            T.add("step.total_ms", prep + pre + T.m["step.publish_ms"].v.back());
            done = step;
            if (step % every == 0 || step == steps) {
                gauss_series.push_back({{"step", step}, {"time_s", res->time_s()}, {"gauss", gauss_ratio(*res)}});
                for (const auto& [k, v] : res->device_field().device_profile()) T.add("device." + k, v);
            }
        } catch (const std::exception& e) { refusal = e.what(); break; }
    }
    const double wall = ms_since(t_all);
    summary["steps_done"] = done; summary["refusal"] = refusal; summary["wall_ms"] = wall; summary["ms_per_step"] = done ? wall / done : 0.;
    summary["phases"] = T.json();
    summary["guards"] = {{"max_independent_residual_over_target", max_resid}, {"max_gauss_over_bound_maxima", max_gauss_over},
                         {"max_divB_over_bound", max_divB_over}, {"max_energy_over_bound", max_energy}, {"max_abs_field_heat_balance_J", max_heat}};
    summary["control"] = {{"B_max_T", max_bmax}, {"omega_c_dt_max", max_omega_c}, {"third_difference_rel_max", max_third}, {"increment_rel_max", max_incr}};
    summary["gauss_series"] = gauss_series;
    {
        // per-node absolute Gauss ratio drift per simulated second (linear fit of the worst ratio over the run)
        const auto& a = gauss_series.front(); const auto& z = gauss_series.back();
        const double r0 = a["gauss"]["max_node_ratio"], r1 = z["gauss"]["max_node_ratio"], t0 = a["time_s"], t1 = z["time_s"];
        summary["gauss_ratio_rate_per_s"] = t1 > t0 ? (r1 - r0) / (t1 - t0) : 0.;
        summary["gauss_ratio_projected_refusal_time_s"] = (t1 > t0 && r1 > r0) ? t1 + (1 - r1) / ((r1 - r0) / (t1 - t0)) : -1.;
    }
    {
        axial::State fs; hc::ChargeLedger fq; res->accepted_state(fs, fq);
        summary["B_bounds_final"] = b_bounds(g, fs);
        summary["final_time_s"] = fs.time_s;
        summary["final_current_A"] = model->drive().current(fs.time_s);
    }
    summary["final_field_sha256"] = res->field_sha256();
    summary["digests"] = {{"computed", res->device_field().digest_stats().digests}, {"leaves_hashed_on_host", res->device_field().digest_stats().leaves_hashed_on_host}};
    summary["device"] = J::parse(dev->device_json());
    summary["memory"] = J::parse(dev->memory_json());
    summary["memory_model_bytes"] = resident_bytes(g, [&] { std::size_t m = 0; for (const auto& c : g.cells) m += c.body >= 0; return m; }());
    summary["max_rss_kB"] = rss_kb("VmHWM:");
    summary["gpu_hold_ms"] = ms_since(t_gpu);
    res.reset();
    summary["teardown"] = J::parse(dev->close_json());
    std::ofstream(out / "gpu-vacuum.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return refusal.empty() ? 0 : 3;
}

int gpu_dtscan(const std::filesystem::path& out, const std::string& trace, const std::string& nvrtc, const std::string& uuid, int nT, int nA, int nAir,
               long steps, std::vector<double> dts) {
    std::filesystem::create_directories(out);
    J summary{{"schema", "resolution-series-gpu-dtscan-v1"}, {"spec", {nT, nA, nAir, .009}}, {"steps_per_dt", steps}};
    const double cap = std::max(1e-14, *std::max_element(dts.begin(), dts.end()));
    auto t = steady::now();
    auto cpu_model = std::make_shared<const hc::Model>(hc::MeshSpec{nT, nA, nAir, .009}, trace);
    summary["cpu_model_construct_ms"] = ms_since(t);
    const auto& g = cpu_model->mesh(); const auto o = mf::geometry(g); const auto seed = cpu_model->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    const auto q0 = declared(*cpu_model, seed, hc::charge_of_census(g, o, none));
    const auto t_gpu = steady::now();
    const auto desc = probe_descriptor(trace, nvrtc, uuid);
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{nT, nA, nAir, .009}, trace, std::make_shared<hc::ResidentBackend>(desc), hc::CarrierLimits{},
                                                   hc::ThermalDeclaration{}, cap);
    auto dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    J rows = J::array();
    for (double dt : dts) {
        Stats T; std::string refusal; long done = 0; double max_resid = 0;
        {
            hc::ResidentOwner res(model, dev, seed, q0, o, none, gp::ResidentOptions{});
            const auto t0 = steady::now();
            for (long k = 1; k <= steps; ++k) {
                try {
                    t = steady::now();
                    auto p = res.prepare_vacuum(dt);
                    Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
                    hc::need(res.preflight_vacuum_final(p, b), "dt scan preflight");
                    const auto& sr = p.solve_record();
                    T.add("iterations", sr.iterations); T.add("residual_over_target", sr.independent_residual / sr.target);
                    max_resid = std::max(max_resid, sr.independent_residual / sr.target);
                    for (const auto& [kk, v] : p.step().host_ms) T.add("prepare." + kk, v);
                    res.publish_vacuum_prevalidated(std::move(p));
                    T.add("step_ms", ms_since(t));
                    done = k;
                } catch (const std::exception& e) { refusal = e.what(); break; }
            }
            T.add("wall_ms", ms_since(t0));
        }
        const double courant = 299792458. * dt / g.dx;
        rows.push_back({{"dt_s", dt}, {"courant_transverse", courant}, {"steps_done", done}, {"refusal", refusal}, {"max_residual_over_target", max_resid},
                        {"stats", T.json()}});
        std::cout << "dt " << dt << " done " << done << " iters " << T.json()["iterations"].dump() << " " << refusal << std::endl;
        std::ofstream(out / "gpu-dtscan.json") << J{{"partial", rows}}.dump(1) << '\n';
    }
    summary["rows"] = rows;
    summary["memory"] = J::parse(dev->memory_json());
    summary["gpu_hold_ms"] = ms_since(t_gpu);
    summary["teardown"] = J::parse(dev->close_json());
    std::ofstream(out / "gpu-dtscan.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return 0;
}

// Closed coupled Hall step on the resident owner (Hall closed loop, no native/chemistry), the population of the Hall
// v3.1 coupled-run test (co-located electron/H2O+ pairs, 1-300 eV electrons, weight 1e3) on a chosen mesh.
hall_driver::MotionBinding binding_for(double t, std::uint64_t rev, const std::string& mesh, const std::string& material, const std::string& policy, long step,
                                       double dt) {
    hall_driver::MotionBinding b; b.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = rev; b.staged_arena_revision = rev + 1; b.begin_s = t; b.dt_s = dt; b.end_s = b.begin_s + b.dt_s;
    b.staged_arena_sha256 = hc::sha256("resident/staged/" + std::to_string(step)); b.mesh_sha256 = mesh;
    b.material_sha256 = material; b.species_sha256 = hc::sha256("resident/species");
    b.common_clock_sha256 = hc::sha256("resident/clock/" + std::to_string(step)); b.kinetic_policy_sha256 = policy;
    return b;
}
int gpu_coupled(const std::filesystem::path& out, const std::string& trace, const std::string& nvrtc, const std::string& uuid, int nT, int nA, int nAir,
                long steps, std::size_t carriers, double dt) {
    std::filesystem::create_directories(out);
    J summary{{"schema", "resolution-series-gpu-coupled-v1"}, {"spec", {nT, nA, nAir, .009}}, {"dt_s", dt}, {"steps_requested", steps}, {"carriers", carriers}};
    hc::CarrierLimits lim{std::max<std::size_t>(100000, carriers), 1ull << 34};
    auto t = steady::now();
    auto cpu_model = std::make_shared<const hc::Model>(hc::MeshSpec{nT, nA, nAir, .009}, trace, hc::make_cpu_pcg(), lim);
    summary["cpu_model_construct_ms"] = ms_since(t);
    const auto& g = cpu_model->mesh(); const auto o = mf::geometry(g); const auto seed = cpu_model->dc_seed();
    auto census = carriers == 2 ? mf::pair({g.side / 2, g.side / 2, 0.}, 1.) : mf::pairs(carriers, 7, 1e3);
    const auto q0 = declared(*cpu_model, seed, hc::charge_of_census(g, o, census, lim.maximum_carriers));
    const auto t_gpu = steady::now();
    const auto desc = probe_descriptor(trace, nvrtc, uuid);
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{nT, nA, nAir, .009}, trace, std::make_shared<hc::ResidentBackend>(desc), lim,
                                                   hc::ThermalDeclaration{}, std::max(dt, 1e-14));
    auto dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    std::string refusal; long done = 0; Stats T;
    {
        hc::ResidentOwner res(model, dev, seed, q0, o, census, gp::ResidentOptions{});
        auto provider = std::make_unique<gp::GpuMotionProvider>(dev, res.motion_context(), carriers);
        summary["memory_after_owner"] = J::parse(dev->memory_json());
        const std::string src = hc::sha256("resident/source");
        const auto t_all = steady::now();
        for (long step = 1; step <= steps; ++step) {
            try {
                t = steady::now();
                auto b = binding_for(res.time_s(), res.revision(), res.motion_mesh_sha256(), model->material_sha256(), res.motion_context().policy_sha256, step, dt);
                auto p = hc::prepare_closed_proposal(res, o, census, b, src, *provider, hc::StepInputs{});
                const double prep = ms_since(t); t = steady::now();
                const auto& bb = p.motion().binding();
                Boundary bd{p.candidate_field_sha256(), p.candidate_heat_sha256(), src, bb.mesh_sha256, bb.material_sha256, bb.common_clock_sha256, p.end_s(), p.endpoint_revision()};
                const bool ok = res.preflight_final(p, census, bd);
                const double pre = ms_since(t);
                hc::need(ok, "coupled preflight");
                const auto& TT = p.timings();
                T.add("step.prepare_closed_proposal_ms", prep); T.add("step.preflight_final_ms", pre);
                T.add("proposal.admit_census_ms", 1e3 * TT.admit_census_s); T.add("proposal.motion_ms", 1e3 * TT.motion_s); T.add("proposal.support_ms", 1e3 * TT.support_s);
                T.add("proposal.field_prepare_ms", 1e3 * TT.field_prepare_s); T.add("proposal.motion_seal_ms", 1e3 * TT.motion_v2_seal_s);
                T.add("proposal.heat_digests_ms", 1e3 * TT.heat_s); T.add("proposal.iterations", TT.iterations);
                T.add("solve.iterations", TT.last_solve.iterations);
                auto next = hc::next_census(census, p.motion());
                t = steady::now();
                res.publish_prevalidated(std::move(p));
                T.add("step.publish_ms", ms_since(t));
                T.add("step.total_ms", prep + pre + T.m["step.publish_ms"].v.back());
                census = std::move(next);
                done = step;
            } catch (const std::exception& e) { refusal = e.what(); break; }
        }
        summary["wall_ms"] = ms_since(t_all);
        summary["memory"] = J::parse(dev->memory_json());
        provider.reset();
    }
    summary["steps_done"] = done; summary["refusal"] = refusal; summary["phases"] = T.json();
    summary["ms_per_step"] = done ? double(summary["wall_ms"]) / done : 0.;
    summary["max_rss_kB"] = rss_kb("VmHWM:");
    summary["gpu_hold_ms"] = ms_since(t_gpu);
    summary["teardown"] = J::parse(dev->close_json());
    std::ofstream(out / "gpu-coupled.json") << summary.dump(1) << '\n';
    std::cout << summary.dump(1) << std::endl;
    return refusal.empty() ? 0 : 3;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout << std::setprecision(17);
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "mesh-report" && argc >= 6 && (argc - 3) % 3 == 0) {
            std::vector<std::array<int, 3>> specs;
            for (int i = 3; i + 2 < argc; i += 3) specs.push_back({std::stoi(argv[i]), std::stoi(argv[i + 1]), std::stoi(argv[i + 2])});
            return mesh_report(argv[2], specs);
        }
        if (mode == "graded-dc" && argc == 9)
            return graded_dc(argv[2], argv[3], std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]));
        if (mode == "graded-sizes" && argc == 7) {
            auto list = [](const std::string& s) { std::vector<double> v; std::stringstream ss(s); std::string x; while (std::getline(ss, x, ',')) v.push_back(std::stod(x)); return v; };
            return graded_sizes(argv[2], list(argv[3]), list(argv[4]), list(argv[5]), std::stod(argv[6]));
        }
        if (mode == "dc-floor" && argc == 8) return dc_floor(argv[2], argv[3], std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]));
        if (mode == "admit" && (argc == 7 || argc == 8))
            return admit(argv[2], argv[3], std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), argc == 8 ? std::stod(argv[7]) : .009);
        if (mode == "cpu-vacuum" && argc == 9)
            return cpu_vacuum(argv[2], argv[3], std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), std::stol(argv[7]), std::stod(argv[8]));
        if (mode == "gpu-vacuum" && argc >= 13)
            return gpu_vacuum(argv[2], argv[3], argv[4], argv[5], std::stoi(argv[6]), std::stoi(argv[7]), std::stoi(argv[8]), std::stol(argv[9]),
                              std::stod(argv[10]), std::stod(argv[11]), argc > 12 ? std::stol(argv[12]) : 100, argc > 13 ? std::stod(argv[13]) : .009);
        if (mode == "gpu-dtscan" && argc >= 11) {
            std::vector<double> dts; for (int i = 10; i < argc; ++i) dts.push_back(std::stod(argv[i]));
            return gpu_dtscan(argv[2], argv[3], argv[4], argv[5], std::stoi(argv[6]), std::stoi(argv[7]), std::stoi(argv[8]), std::stol(argv[9]), dts);
        }
        if (mode == "gpu-coupled" && argc == 12)
            return gpu_coupled(argv[2], argv[3], argv[4], argv[5], std::stoi(argv[6]), std::stoi(argv[7]), std::stoi(argv[8]), std::stol(argv[9]),
                               std::stoul(argv[10]), std::stod(argv[11]));
        std::cerr << "usage: see the header of field_scaling.cpp\n";
        return 64;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return 1;
    }
}
