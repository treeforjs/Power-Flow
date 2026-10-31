// Identity parity trace: the same resident run (genuine Hall v2 motion at
// zero total field, mixed forced-event exposures, K = 4) compiled unchanged
// against native v5 and native v6. Every accepted identity is written per
// step; without births or sinks v6 must reproduce v5 except the anchor
// records (source identity, anchor v6 layout).
// usage: parity_trace OUTPUT.json STEPS
#include "resident_fixture.hpp"
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
  try {
    if (argc != 3)
      throw std::runtime_error("usage: parity_trace OUTPUT.json STEPS");
    const std::uint64_t steps = std::strtoull(argv[2], nullptr, 10);
    m::FullTimeConfiguration c;
    c.history.anchor_interval_steps = 4;
    auto a = m::PicArena::create_full_time(resident::hall_options(), states(), channels(),
                                           resident::population(), resident::hall_genesis(), c);
    J rows = J::array();
    const auto g = a->full_time_boundary();
    J genesis_row = {{"accepted_arena_sha256", g->accepted_arena_sha256},
                     {"history_policy_sha256", a->full_time_history_status().history_policy_sha256}};
    resident::StepState st;
    for (std::uint64_t k = 0; k < steps; ++k) {
      resident::StepOptions o;
      o.first_events = k % 3 == 1 ? 1 : 0;
      o.second_events = k % 4 == 2 ? 2 : 0;
      resident::full_step(*a, o, st);
      const auto b = a->full_time_boundary();
      const auto h = a->full_time_history_status();
      rows.push_back({{"step", b->accepted_step},
                      {"accepted_arena_sha256", b->accepted_arena_sha256},
                      {"accepted_arena_revision", b->accepted_arena_revision},
                      {"rng_sha256", b->rng_sha256},
                      {"field_sha256", b->field_sha256},
                      {"time_s", b->time_s},
                      {"cumulative_rows", h.cumulative_rows},
                      {"cumulative_bytes", h.cumulative_bytes},
                      {"anchors", h.anchors},
                      {"anchor_head_sha256", h.anchor_head_sha256}});
    }
    a->validate_full_time_history();
    J out = {{"source_identity", m::PicArena::source_identity()},
             {"genesis", genesis_row},
             {"anchor_root_sha256", a->full_time_history_status().anchor_root_sha256},
             {"steps", rows},
             {"forbidden_calls", forbidden_calls}};
    std::ofstream(argv[1]) << out.dump(2) << '\n';
    std::cout << "parity trace " << steps << " steps\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
