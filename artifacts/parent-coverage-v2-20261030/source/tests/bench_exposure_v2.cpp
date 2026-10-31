#include "parent_coverage_v2_subflow.hpp"
#include "parent_coverage_subflow.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>

// Short exposure benchmark: cost of one legacy frozen-stage exposure through
// the coverage-v1 and coverage-v2 subflows on the same populations. One
// thread; fixtures are built outside the timed region; median of repeats.
namespace nw = raikou::native_water;
namespace cv = raikou::native_water::coverage_v2;
namespace c1 = raikou::native_water::coverage_v1;
namespace m = raikou::mfem_backend;
namespace pc = raikou::pic_collisions;
namespace w = raikou::pic_water;
using J = nw::Json;
constexpr double mn = 2.9915e-26;
m::GapPicOptions options() {
  m::GapPicOptions o;o.nx = o.ny = 45;o.width_m = o.gap_m = .009;o.depth_m = .025;
  o.finite_foils = true;o.finite_axial_extent = true;o.short_z_m = 0;o.feed_z_m = .025;
  o.lower_foil = {.002, .007, .0018, .002};o.upper_foil = {.002, .007, .007, .0072};return o;
}
m::PicVector electron_u(double E) { double x = E * pc::qe / (w::me * pc::c * pc::c);return {pc::c * std::sqrt(x * (x + 2)), 0, 0}; }
std::vector<m::ArenaStateSpec> states() {
  return {{"electron", 0, 0, -1, w::me, 0, 0}, {"H2O", 2, 1, 0, mn, 0, 0},
          {"H2O+", 2, 1, 1, mn - w::me + w::threshold_eV * pc::qe / (pc::c * pc::c), 0, w::threshold_eV * pc::qe}};
}
struct Fixture {
  m::GapPic pic{options()};std::unique_ptr<m::PicArena> arena;
  Fixture(const cv::ArenaBinding& a, const std::vector<double>& energies, double water) {
    std::vector<m::GapParticle> batch;
    for (double E : energies) { m::GapParticle e;e.weight = 1;e.mass_kg = w::me;e.position_m = {.004, .004, .0125};e.proper_velocity_m_s = electron_u(E);batch.push_back(e); }
    m::GapParticle n;n.weight = water;n.mass_kg = mn;n.charge_C = 0;n.position_m = {.00403, .004, .0125};batch.push_back(n);
    auto ids = pic.add_authored_batch(batch);pic.solve_fields();
    std::vector<m::ArenaPacket> packets;
    for (std::size_t i = 0; i < energies.size(); ++i) packets.push_back({"e" + std::to_string(i), "electron", "authored_fixture", "", {1, 0}, {ids[i]}});
    packets.push_back({"n", "H2O", "authored_fixture", "", {water, 0}, {ids.back()}});
    arena = std::make_unique<m::PicArena>(pic, states(), std::vector<m::ArenaChannel>{{w::parent_channel, "electron", "H2O", {"electron", "electron", "H2O+"},
      a.dataset_sha256, a.kernel_manifest_sha256, a.kernel_binding_file_sha256, a.model}}, packets, "native-water-subflow-fixture/v1");
  }
};
template<class Flow> nw::Interval interval(m::PicArena& arena, const Flow& flow, double duration) {
  nw::Interval i;i.end_s = duration;i.common_clock_sha256 = nw::sha("bench/0");
  i.phase.state_prefix = arena.state_prefix();i.phase.mesh_prefix = flow.grid_identity();i.phase.clock_prefix = i.common_clock_sha256;return i;
}
J strip(const J& audit) { J out = J::array();for (auto e : audit.at("events")) { e.erase("native");out.push_back(e); }return out; }

int main(int argc, char** argv) {
  try {
    if (argc != 6) throw std::runtime_error("dataset binding outer3 kernel-binding output-directory");
    std::ifstream f(argv[2]);J binding;f >> binding;
    nw::Binding b{binding.at("dataset_sha256"), binding.at("kernel_sha256"), binding.at("kernel_binding_file_sha256")};
    auto data = std::make_shared<const w::Data>(w::Data::read(argv[1]));
    w::ParentBridge bridge(data, w::parent_model, b.dataset_sha256, b.kernel_sha256, mn);
    auto bands = std::make_shared<const cv::HighEnergyBands>(cv::HighEnergyBands::read(argv[3],
        "c436020d3f2d07902fb059cfc0a269b2e914a3cbb2789300093d3631df896551", *data, b.dataset_sha256, argv[4], b.kernel_sha256,
        b.kernel_binding_file_sha256));
    cv::CoveragePolicy open2;open2.max_exposure_unresolved_events = open2.max_cumulative_unresolved_events = std::numeric_limits<double>::infinity();
    open2.max_unresolved_reaction_fraction = 1;open2.admit_unbounded = true;
    c1::CoveragePolicy open1;open1.max_exposure_unresolved_events = open1.max_cumulative_unresolved_events = std::numeric_limits<double>::infinity();
    open1.max_unresolved_reaction_fraction = 1;open1.admit_unbounded = true;
    const c1::CoverageScope s1(data, b.dataset_sha256, bridge, open1);
    const cv::CoverageScope s2(data, b.dataset_sha256, bridge, bands, open2), s2eq(data, b.dataset_sha256, bridge, open2);
    nw::Limits limits;limits.reactions = 1000;nw::Grid grid(options(), {45, 45, 5}, limits.cells);
    const int repeats = 7;
    // Sub-keV population (identical physics in v1 and v2) and a population
    // with 1 keV-1 MeV electrons (v1 excludes and records those pairs).
    std::vector<double> low, high;
    for (int i = 0; i < 64; ++i) low.push_back(20. + 14. * i);
    for (int i = 0; i < 64; ++i) high.push_back(i % 2 ? 20. + 28. * i : 1000. * std::pow(1000., i / 63.));
    J rows = J::array();
    for (const auto& [label, pop, dt] : std::vector<std::tuple<std::string, std::vector<double>, double>>{{"64 electrons 20-900 eV", low, 2e-17},
                                                                                                         {"32 electrons 1 keV-1 MeV + 32 below 1 keV", high, 2e-17}}) {
      auto run = [&](auto make_flow, const cv::ArenaBinding& reg) {
        std::vector<double> t;J first;std::size_t events = 0, pairs = 0;
        for (int r = 0; r < repeats; ++r) {
          Fixture fx(reg, pop, 1e20);auto flow = make_flow(*fx.arena);
          const auto i = interval(*fx.arena, *flow, dt);
          const auto t0 = std::chrono::steady_clock::now();
          auto trial = flow->prepare(i);
          const auto t1 = std::chrono::steady_clock::now();
          t.push_back(std::chrono::duration<double>(t1 - t0).count());
          if (!r) { first = strip(trial.audit());events = trial.audit().at("events").size();pairs = trial.audit().at("pairs"); }
        }
        std::sort(t.begin(), t.end());
        return std::tuple<double, J, std::size_t, std::size_t>{t[repeats / 2], first, events, pairs};
      };
      const auto reg1 = s2eq.arena_binding(b.kernel_binding_file_sha256), reg2 = s2.arena_binding(b.kernel_binding_file_sha256);
      const auto [tv1, ev1, n1, p1] = run([&](m::PicArena& a) { return std::make_unique<c1::ParentSubflow>(a, grid.counts, b, limits, bridge, s1, 11, 12); }, reg1);
      const auto [teq, eveq, neq, peq] = run([&](m::PicArena& a) { return std::make_unique<cv::ParentSubflow>(a, grid.counts, b, limits, bridge, s2eq, 11, 12); }, reg1);
      const auto [tv2, ev2, n2, p2] = run([&](m::PicArena& a) { return std::make_unique<cv::ParentSubflow>(a, grid.counts, b, limits, bridge, s2, 11, 12); }, reg2);
      rows.push_back({{"population", label}, {"exposure_s", dt}, {"repeats", repeats},
                      {"v1", {{"median_s", tv1}, {"events", n1}, {"pair_visits", p1}}},
                      {"v2_without_bands", {{"median_s", teq}, {"events", neq}, {"pair_visits", peq}, {"events_identical_to_v1", eveq == ev1}}},
                      {"v2", {{"median_s", tv2}, {"events", n2}, {"pair_visits", p2}, {"events_identical_to_v1", ev2 == ev1}}},
                      {"v2_over_v1_time", tv2 / tv1}, {"v2_without_bands_over_v1_time", teq / tv1},
                      {"v2_seconds_per_pair_visit", tv2 / static_cast<double>(p2)}, {"v1_seconds_per_pair_visit", tv1 / static_cast<double>(p1)}});
      std::cout << label << ": v1 " << tv1 << " s (" << n1 << " events), v2 " << tv2 << " s (" << n2 << " events), ratio " << tv2 / tv1 << "\n";
    }
    J out{{"benchmark", "legacy frozen-stage exposure, one thread, median of repeats"}, {"rows", rows}};
    std::ofstream(std::string(argv[5]) + "/benchmark.json") << out.dump(2) << '\n';
  } catch (const std::exception& e) { std::cerr << e.what() << '\n';return 1; }
}
