// Motion admission timing: the native v5/v6 JSON route (Hall v2 record seal,
// JSON record build, native JSON admission) versus the native v6 binary
// motion v3 route (adapter bind, typed admission), for the same genuine Hall
// corrective v2 CPU motion of an N-carrier resident lattice at zero total
// field. Built against native v6 (both routes) and, with BENCH_V5, against
// the sealed native v5 archive (JSON route only).
// usage: bench_motion OUTPUT.json CARRIERS REPEATS
#include "raikou/mfem_pic_full_time.h"
#include "hall-finite-foil-em-20261029/axial_mesh.hpp"
#include "motion_receipt.hpp" // Hall corrective v2
#ifndef BENCH_V5
#include "raikou/mfem_pic_hall_motion_v3.h"
#endif
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <openssl/sha.h>

namespace m = raikou::mfem_backend;
using J = m::ArenaJson;
namespace {
std::string sha(const std::string &s) {
  unsigned char d[32];
  SHA256(reinterpret_cast<const unsigned char *>(s.data()), s.size(), d);
  static const char *h = "0123456789abcdef";
  std::string o;
  for (auto c : d) {
    o.push_back(h[c >> 4]);
    o.push_back(h[c & 15]);
  }
  return o;
}
long status_kib(const char *key) {
  FILE *f = std::fopen("/proc/self/status", "r");
  if (!f)
    return -1;
  char line[256];
  long value = -1;
  const auto n = std::strlen(key);
  while (std::fgets(line, sizeof line, f))
    if (!std::strncmp(line, key, n)) {
      value = std::strtol(line + n + 1, nullptr, 10);
      break;
    }
  std::fclose(f);
  return value;
}
m::GapPicOptions options() {
  m::GapPicOptions o;
  o.nx = o.ny = 45;
  o.width_m = o.gap_m = .009;
  o.depth_m = .025;
  o.finite_foils = o.finite_axial_extent = true;
  o.short_z_m = -.0125;
  o.feed_z_m = .0125;
  o.lower_foil = {.002, .007, .0018, .002};
  o.upper_foil = {.002, .007, .007, .0072};
  return o;
}
std::vector<m::FullTimeAuthoredPacket> lattice(std::size_t n) {
  std::vector<m::FullTimeAuthoredPacket> out;
  out.reserve(n);
  const m::ArenaCount four(4, 0);
  for (std::size_t i = 0; i < n; ++i) {
    const double x = .0021 + 4.6e-3 * double(i % 997) / 997.0;
    const double y = .0021 + 4.6e-3 * double((i / 997) % 1009) / 1009.0;
    const double z = -.01 + 1e-8 * double(i / (997 * 1009));
    char id[24];
    std::snprintf(id, sizeof id, "p%07zu", i);
    const bool electron = i % 2 == 0;
    out.push_back({id, electron ? "e" : "n", "manufactured/lattice", four, {x, y, z},
                   {0, 0, electron ? 1e6 : 0}});
  }
  return out;
}
double since(std::chrono::steady_clock::time_point &t) {
  const auto now = std::chrono::steady_clock::now();
  const double d = std::chrono::duration<double>(now - t).count();
  t = now;
  return d;
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 4)
      throw std::runtime_error("usage: bench_motion OUTPUT.json CARRIERS REPEATS");
    const std::size_t n = std::strtoull(argv[2], nullptr, 10);
    const int repeats = std::atoi(argv[3]);
    const axial::Mesh mesh;
    const auto o = options();
    m::FullTimeConfiguration c;
    c.limits.maximum_carriers = std::max<std::size_t>(n, 100000);
    c.limits.memory_budget_bytes = c.limits.maximum_carriers * m::full_time_bytes_per_carrier + m::full_time_fixed_bytes;
    c.limits.maximum_identity_entries = std::max<std::size_t>(2 * n, 1u << 22);
    m::FullTimeGenesisBinding g;
    g.time_s = 26.31715e-9;
    g.origin_provenance = "manufactured/admission-benchmark-no-simulation";
    g.field_sha256 = sha("field");
    g.mesh_sha256 = hall_driver::motion_mesh_sha256(mesh, o);
    g.material_sha256 = sha("material");
    g.frame_sha256 = sha("frame");
    g.common_clock_sha256 = sha("clock");
    g.heat_sha256 = sha("heat");
    g.source_sha256 = sha("source");
    g.rng_sha256 = sha("rng");
    const std::vector<m::ArenaStateSpec> states = {{"e", 0, 0, -1, m::pic_me, 0, 0}, {"n", 2, 1, 0, 2.9915e-26, 0, 0}};
    const std::vector<m::ArenaChannel> channels = {{"conditional-storage-control", "e", "n", {"e", "n"}, sha("data"),
                                                    sha("kernel"), sha("binding"), "manufactured-storage-control-only"}};
    auto t = std::chrono::steady_clock::now();
    auto arena = m::PicArena::create_full_time(o, states, channels, lattice(n), g, c);
    J out = {{"carriers", n},
#ifdef BENCH_V5
             {"native", "v5"},
#else
             {"native", "v6"},
#endif
             {"source_identity", m::PicArena::source_identity()},
             {"create_s", since(t)}};
    const std::string mesh_json = mesh.json().dump();
    const auto context = hall_corrective_v2::make_motion_context(mesh, o, mesh_json, std::max<std::size_t>(n, 100000));
    J runs = J::array();
    for (int r = 0; r < repeats; ++r) {
      const auto bnd = arena->full_time_boundary();
      m::FullTimeStepBinding b{bnd->accepted_step, bnd->accepted_arena_revision, bnd->time_s, bnd->time_s + 5e-15,
                               5e-15, bnd->accepted_arena_sha256, bnd->field_sha256, bnd->mesh_sha256,
                               bnd->material_sha256, sha("bench-clock/" + std::to_string(r))};
      auto outer = arena->begin_full_time_step(b, n, 32);
      m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s, .5 * b.dt_s, 0, b.accepted_field_sha256};
      auto s1 = outer.first_exposure(e1);
      auto tok = s1.prepare_compound_with_completion(1, [&](m::CollisionScope &) {
        m::ArenaStochasticCompletion x;
        const auto ch = channels[0];
        x.integration_owner_id = "bench";
        x.before_rng_clock_sha256 = bnd->rng_sha256;
        x.after_rng_clock_sha256 = sha("after");
        x.common_clock_sha256 = b.common_clock_sha256;
        x.channel_id = ch.id;
        x.dataset_sha256 = ch.dataset_sha256;
        x.kernel_manifest_sha256 = ch.kernel_manifest_sha256;
        x.kernel_binding_file_sha256 = ch.kernel_binding_file_sha256;
        x.model = ch.model;
        x.stage_sequence = 0;
        x.frozen_physical_stage_s = b.begin_s;
        x.common_interval_begin_s = b.begin_s;
        x.common_interval_end_s = b.end_s;
        return x;
      });
      if (!s1.try_publish(tok))
        throw std::runtime_error("first exposure");
      J row = J::object();
      t = std::chrono::steady_clock::now();
      auto trial1 = outer.fork_transport_trial();
      row["fork_s"] = since(t);
      const auto census = trial1.native_census(n);
      row["native_census_entries_s"] = since(t);
      static const std::vector<double> E(mesh.ne, 0.0), B(mesh.nb, 0.0);
      hall_driver::MotionBinding mb;
      mb.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
      mb.accepted_step = b.accepted_step;
      mb.staged_arena_revision = trial1.staged_arena_revision();
      mb.begin_s = b.begin_s;
      mb.end_s = b.end_s;
      mb.dt_s = b.dt_s;
      mb.staged_arena_sha256 = trial1.state_prefix();
      mb.trial_field_sha256 = hall_driver::total_field_sha256(E, B);
      mb.mesh_sha256 = b.mesh_sha256;
      mb.material_sha256 = b.material_sha256;
      mb.species_sha256 = sha("species");
      mb.common_clock_sha256 = b.common_clock_sha256;
      mb.kinetic_policy_sha256 = context.policy_sha256;
      t = std::chrono::steady_clock::now();
      auto receipt = hall_corrective_v2::MotionAssembly::cpu(context, census, mb, {}, E, B);
      row["hall_cpu_motion_incl_census_and_v3_digest_s"] = since(t);
      // JSON route (native v5/v6 v2 recognizer).
      hall_corrective_v2::MotionAssembly::seal_v2(receipt);
      row["hall_v2_json_seal_s"] = since(t);
      {
        auto admission = m::bind_hall_total_field_motion(receipt);
        row["v2_json_record_build_s"] = since(t);
        auto moved = trial1.prepare_total_field_motion(
            {b, receipt.seal_sha256(), b.accepted_field_sha256, sha("f"), sha("h"), bnd->source_sha256},
            std::move(admission));
        row["native_v2_json_admission_s"] = since(t);
      }
      row["v2_json_release_s"] = since(t);
#ifndef BENCH_V5
      {
        auto trial2 = outer.fork_transport_trial();
        t = std::chrono::steady_clock::now();
        auto admission = m::bind_hall_total_field_motion_v3(receipt);
        row["v3_adapter_bind_s"] = since(t);
        auto moved = trial2.prepare_total_field_motion(
            {b, receipt.record_v3_sha256(), b.accepted_field_sha256, sha("f"), sha("h"), bnd->source_sha256},
            std::move(admission));
        row["native_v3_admission_cold_census_s"] = since(t);
      }
      {
        auto trial3 = outer.fork_transport_trial();
        t = std::chrono::steady_clock::now();
        const auto census_sha = trial3.census_sha256();
        row["native_census_digest_once_per_trial_s"] = since(t);
        if (census_sha != receipt.census_sha256())
          throw std::runtime_error("census digest mismatch");
        auto admission = m::bind_hall_total_field_motion_v3(receipt);
        row["v3_adapter_bind_warm_s"] = since(t);
        auto moved = trial3.prepare_total_field_motion(
            {b, receipt.record_v3_sha256(), b.accepted_field_sha256, sha("f"), sha("h"), bnd->source_sha256},
            std::move(admission));
        row["native_v3_admission_cached_census_s"] = since(t);
        outer.accept_transport_trial(std::move(moved));
        row["accept_s"] = since(t);
      }
#endif
      row["VmRSS_KiB"] = status_kib("VmRSS:");
      runs.push_back(row);
    }
    out["runs"] = runs;
    out["VmHWM_KiB"] = status_kib("VmHWM:");
    std::ofstream(argv[1]) << out.dump(2) << '\n';
    std::cout << out.dump() << '\n';
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
