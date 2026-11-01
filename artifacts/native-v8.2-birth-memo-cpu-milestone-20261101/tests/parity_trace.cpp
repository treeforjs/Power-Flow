// Identity parity trace: the same resident run (genuine Hall v2 motion at
// zero total field, mixed forced-event exposures, K = 4) compiled against the
// sealed native v6 (v6 JSON route, -DPARITY_NATIVE_V6) and native v7 (v7
// adapter: exact sums, Merkle record, census v2). Every accepted identity and
// the physical particle state (every native carrier with its packet ID and
// state, canonical bytes) is written per step. v7 must reproduce the v6
// physical state bit for bit; the accepted identities differ by the motion
// identities of API-007.txt section 5.
// usage: parity_trace OUTPUT.json STEPS
#ifndef PARITY_REFERENCE
#include "device_mode.hpp"
#endif
#include "resident_fixture.hpp"
#include "raikou/mfem_pic_canonical.h"
namespace {
std::string physical_state(const m::PicArena &a) {
  auto rows = a.rate_population(1u << 20);
  std::sort(rows.begin(), rows.end(), [](const auto &x, const auto &y) { return x.native.id < y.native.id; });
  m::canonical::Writer w;
  w.tag("parity-physical-state-v1");
  w.u64(rows.size());
  for (const auto &r : rows) {
    w.text(r.packet_id);
    w.text(r.state);
    w.raw(&r.native, sizeof r.native);
  }
  return m::canonical::sha256_hex(w.bytes);
}
} // namespace
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
#ifndef PARITY_REFERENCE
    v71::install_from_environment();
    struct ExecutionReport {
      std::string path;
      ~ExecutionReport() { v71::write_report(path); }
    } execution_report{std::string(argv[1]) + ".execution.json"};
#endif

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
                      {"physical_state_sha256", physical_state(*a)},
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
