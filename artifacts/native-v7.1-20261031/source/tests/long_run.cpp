// Long synthetic resident run. Compiled unchanged against candidate003 (for the
// comparison baseline) and against native v5 (NATIVE_V5 defined). Each step is
// a complete closed full-time step with a genuine Hall v2 motion receipt at
// zero total field; no field solve, no GapPic advance/solve (link traps).
#include "resident_fixture.hpp"
#include <algorithm>
#include <cstdlib>
#include <numeric>
int forbidden_calls = 0;
extern "C" void
__wrap__ZN6raikou12mfem_backend6GapPic12solve_fieldsEv(m::GapPic *) {
  ++forbidden_calls;
  throw std::runtime_error("forbidden native solver trap");
}
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic7advanceEd(m::GapPic *,
                                                                 double) {
  ++forbidden_calls;
  throw std::runtime_error("forbidden native advance trap");
}
int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "usage: long_run OUTPUT.json STEPS SAMPLE_EVERY [K] [CARRIERS]\n";
    return 2;
  }
  const std::string output = argv[1];
  const std::uint64_t steps = std::strtoull(argv[2], nullptr, 10);
  const std::uint64_t sample = std::max<std::uint64_t>(1, std::strtoull(argv[3], nullptr, 10));
  J result = {{"label",
#ifdef NATIVE_V5
               "native-v5"
#else
               "candidate003"
#endif
              },
              {"source_identity", m::PicArena::source_identity()},
              {"requested_steps", steps},
              {"dt_s", 5e-15},
              {"fixture",
               "candidate003 storage-control species/channel/genesis on Hall "
               "geometry; electron u=(0,0,1e6) m/s, neutral at rest; genuine "
               "Hall v2 motion at identically zero total field; zero-event "
               "exposures; synthetic field/heat/RNG digests"}};
  const std::size_t carriers = argc > 5 ? std::strtoull(argv[5], nullptr, 10) : 2;
  const auto population =
      carriers == 2 ? resident::population() : resident::lattice(carriers);
  result["carriers"] = carriers;
  const auto start_rss = resident::status_kib("VmRSS:");
  std::unique_ptr<m::PicArena> arena;
#ifdef NATIVE_V5
  m::FullTimeConfiguration configuration;
  if (argc > 4)
    configuration.history.anchor_interval_steps = std::strtoull(argv[4], nullptr, 10);
  result["anchor_interval_steps"] = configuration.history.anchor_interval_steps;
  if (carriers > configuration.limits.maximum_carriers) {
    configuration.limits.maximum_carriers = carriers;
    configuration.limits.memory_budget_bytes =
        carriers * m::full_time_bytes_per_carrier + m::full_time_fixed_bytes;
  }
  arena = m::PicArena::create_full_time(resident::hall_options(), states(), channels(),
                                        population, resident::hall_genesis(),
                                        configuration);
#else
  arena = m::PicArena::create_full_time(resident::hall_options(), states(), channels(),
                                        population, resident::hall_genesis());
#endif
  resident::StepState state;
  resident::StepOptions options;
  options.maximum_carriers = std::max<std::size_t>(8, carriers);
  J samples = J::array();
  std::vector<double> block;
  std::uint64_t done = 0;
  std::string refusal;
  using clock = std::chrono::steady_clock;
  const auto begin = clock::now();
  double total_s = 0, first_block_mean = -1;
  try {
    for (; done < steps; ++done) {
      const auto t0 = clock::now();
      resident::full_step(*arena, options, state);
      const double dt = std::chrono::duration<double>(clock::now() - t0).count();
      total_s += dt;
      block.push_back(dt);
      if ((done + 1) % sample == 0) {
        std::sort(block.begin(), block.end());
        const double mean =
            std::accumulate(block.begin(), block.end(), 0.0) / block.size();
        if (first_block_mean < 0)
          first_block_mean = mean;
        J row = {{"step", done + 1},
                 {"mean_step_s", mean},
                 {"median_step_s", block[block.size() / 2]},
                 {"max_step_s", block.back()},
                 {"VmRSS_KiB", resident::status_kib("VmRSS:")},
                 {"VmHWM_KiB", resident::status_kib("VmHWM:")}};
#ifdef NATIVE_V5
        const auto h = arena->full_time_history_status();
        row["anchors"] = h.anchors;
        row["anchor_bytes"] = h.anchor_bytes;
        row["window_rows"] = h.window_rows;
        row["window_bytes"] = h.window_bytes;
#endif
        samples.push_back(row);
        block.clear();
        std::cerr << row.dump() << '\n';
      }
    }
  } catch (const std::exception &e) {
    refusal = e.what();
  }
  const double wall = std::chrono::duration<double>(clock::now() - begin).count();
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
  const auto boundary = arena->full_time_boundary();
  result["completed_steps"] = done;
  result["refusal"] = refusal;
  result["wall_s"] = wall;
  result["mean_step_s"] = done ? total_s / done : 0;
  result["samples"] = samples;
  result["start_VmRSS_KiB"] = start_rss;
  result["peak_VmHWM_KiB"] = resident::status_kib("VmHWM:");
  result["ru_maxrss_KiB"] = usage.ru_maxrss;
  result["final_accepted_step"] = boundary->accepted_step;
  result["final_accepted_arena_revision"] = boundary->accepted_arena_revision;
  result["final_time_s"] = boundary->time_s;
  result["final_accepted_arena_sha256"] = boundary->accepted_arena_sha256;
  result["forbidden_native_solver_or_advance_calls"] = forbidden_calls;
  result["phase_totals_s"] = {
      {"native_begin", state.begin_s},
      {"native_first_exposure", state.first_s},
      {"native_fork_census", state.census_s},
      {"hall_provider_motion_and_record", state.hall_s},
      {"native_motion_admission_accept", state.admit_s},
      {"native_second_exposure", state.second_s},
      {"native_seal", state.seal_s},
      {"native_preflight_publish", state.publish_s}};
  result["native_mean_step_s"] =
      done ? (total_s - state.hall_s) / done : 0;
  result["hall_provider_mean_step_s"] = done ? state.hall_s / done : 0;
#ifdef NATIVE_V5
  {
    const auto h = arena->full_time_history_status();
    result["history"] = {{"anchors", h.anchors},
                         {"anchor_bytes", h.anchor_bytes},
                         {"window_rows", h.window_rows},
                         {"window_bytes", h.window_bytes},
                         {"cumulative_rows", h.cumulative_rows},
                         {"cumulative_bytes", h.cumulative_bytes},
                         {"anchor_head_sha256", h.anchor_head_sha256},
                         {"history_policy_sha256", h.history_policy_sha256}};
    const auto t0 = clock::now();
    m::verify_full_time_anchor_chain(arena->export_full_time_anchor_chain(),
                                     h.anchor_root_sha256, h.anchor_head_sha256);
    result["anchor_chain_verification_s"] =
        std::chrono::duration<double>(clock::now() - t0).count();
    result["anchor_chain_verified"] = true;
  }
#endif
  std::ofstream(output) << result.dump(2) << '\n';
  std::cout << done << " steps; " << refusal << '\n';
  return forbidden_calls == 0 && (done == steps || !refusal.empty()) ? 0 : 1;
}
