// Motion admission and host forms timing, native v6 vs native v7, for the same
// genuine Hall corrective v2 CPU motion of an N-carrier resident lattice at
// zero total field.
//  v6 build (BENCH_V6, sealed native v6 archives): Hall-side host forms of v6
//    (census v1 digest, v2 JSON seal, binary v3 record digest, x87 kinetic and
//    work sums) and native v6 admission (JSON route; binary v3 route with the
//    census computed in admission / cached).
//  v7 build: the same Hall-side quantities in the v7 forms computed on the host
//    with the CPU reference (the CPU fallback; the device path is timed by
//    forms/tests/test_forms_gpu), the transitional adapter, and native v7
//    admission (census v2 computed in admission / cached).
// usage: bench_motion OUTPUT.json CARRIERS REPEATS
#include "raikou/mfem_pic_full_time.h"
#include "hall-finite-foil-em-20261029/axial_mesh.hpp"
#include "motion_receipt.hpp" // Hall corrective v2
#ifdef BENCH_V6
#include "raikou/mfem_pic_hall_motion_v3.h"
#else
#include "hall_receipt_v7.hpp"
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
#ifdef BENCH_V6
             {"native", "v6"},
#else
             {"native", "v7"},
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
      row["hall_corrective_v2_cpu_motion_incl_census_v1_and_v3_digest_s"] = since(t);
      static const std::vector<double> E_mid(mesh.ne, 1.0e5);
      const auto &J_cur = receipt.current_A_m();
#ifdef BENCH_V6
      // Hall-side host forms (v6): census v1 digest, JSON seal, v3 digest, x87 sums.
      t = std::chrono::steady_clock::now();
      const auto census_v1 = hall_corrective_v2::canonical_census_sha256(census, std::max<std::size_t>(n, 100000));
      row["hall_census_v1_digest_s"] = since(t);
      hall_corrective_v2::MotionAssembly::seal_v2(receipt);
      row["hall_v2_json_seal_s"] = since(t);
      {
        long double k = 0;
        const auto ends = receipt.endpoints();
        for (std::size_t i = 0; i < ends.size(); ++i)
          k += static_cast<long double>(ends[i].weight) *
               (m::pic_kinetic_energy_J(ends[i]) - m::pic_kinetic_energy_J(census[i].native));
        long double w = 0;
        for (std::size_t e = 0; e < J_cur.size(); ++e)
          w += static_cast<long double>(E_mid[e]) * J_cur[e];
        row["x87_sums_s"] = since(t);
        row["x87_kinetic_J"] = double(k);
        row["x87_work_J"] = b.dt_s * double(w);
        // Binary v3 record digest (native v6 encoder = Hall record_v3).
        std::vector<m::FullTimeMotionPath> paths;
        for (const auto &p : receipt.paths())
          paths.push_back({p.native_id, p.field_start_m, p.field_displacement_m, p.flight_s});
        m::FullTimeMotionRecordV3 x;
        x.endpoints = receipt.endpoints();
        x.paths = paths;
        t = std::chrono::steady_clock::now();
        (void)m::full_time_motion_record_v3_sha256(x);
        row["v3_record_digest_s"] = since(t);
      }
      {
        auto admission = m::bind_hall_total_field_motion(receipt);
        row["v2_json_record_build_s"] = since(t);
        auto moved = trial1.prepare_total_field_motion(
            {b, receipt.seal_sha256(), b.accepted_field_sha256, sha("f"), sha("h"), bnd->source_sha256},
            std::move(admission));
        row["native_v2_json_admission_s"] = since(t);
      }
      row["v2_json_release_s"] = since(t);
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
        if (census_sha != census_v1)
          throw std::runtime_error("census digest mismatch");
        auto admission = m::bind_hall_total_field_motion_v3(receipt);
        auto moved = trial3.prepare_total_field_motion(
            {b, receipt.record_v3_sha256(), b.accepted_field_sha256, sha("f"), sha("h"), bnd->source_sha256},
            std::move(admission));
        row["native_v3_admission_cached_census_s"] = since(t);
        outer.accept_transport_trial(std::move(moved));
        row["accept_s"] = since(t);
      }
#else
      {
        // Hall-side host forms (v7) on the CPU (fallback path): census v2,
        // exact sums, the four section block lists and the header digest.
        auto trial0 = outer.fork_transport_trial();
        const auto &c7 = trial0.census_v7();
        t = std::chrono::steady_clock::now();
        const std::span<const pf_v7::Carrier> rows{reinterpret_cast<const pf_v7::Carrier *>(c7.carriers.data()),
                                                   c7.carriers.size()};
        const auto list = pf_v7::census_block_list(rows, c7.packet_keys, c7.state_indices);
        row["cpu_census_v2_carriers_root_s"] = since(t);
        if (pf_v7::census_v2_sha256(rows.size(), c7.states_sha256, list.root) != c7.sha256)
          throw std::runtime_error("census v2 mismatch");
        const auto ends = receipt.endpoints();
        const auto K = pf_v7::kinetic_change_J(rows, {reinterpret_cast<const pf_v7::Carrier *>(ends.data()), ends.size()}, {});
        row["cpu_exact_kinetic_s"] = since(t);
        const auto W = pf_v7::field_work_J(b.dt_s, E_mid, J_cur);
        row["cpu_exact_work_s"] = since(t);
        long double k = 0;
        for (std::size_t i = 0; i < ends.size(); ++i)
          k += static_cast<long double>(ends[i].weight) *
               (m::pic_kinetic_energy_J(ends[i]) - m::pic_kinetic_energy_J(c7.carriers[i]));
        long double w = 0;
        for (std::size_t e = 0; e < J_cur.size(); ++e)
          w += static_cast<long double>(E_mid[e]) * J_cur[e];
        row["x87_sums_s"] = since(t);
        row["exact_kinetic_J"] = K.value;
        row["x87_kinetic_J"] = double(k);
        row["exact_work_J"] = W.value;
        row["x87_work_J"] = b.dt_s * double(w);
        std::vector<pf_v7::Path> paths(receipt.paths().size());
        std::memcpy(paths.data(), receipt.paths().data(), paths.size() * sizeof(pf_v7::Path));
        t = std::chrono::steady_clock::now();
        const auto sec = pf_v7::motion_sections_v7({reinterpret_cast<const pf_v7::Carrier *>(ends.data()), ends.size()},
                                                   paths, {}, {});
        row["cpu_section_lists_s"] = since(t);
        pf_v7::MotionHeaderV7 h;
        h.census_sha256 = c7.sha256;
        (void)pf_v7::motion_record_v7_sha256(h, sec);
        row["record_header_digest_s"] = since(t);
      }
      {
        auto trial2 = outer.fork_transport_trial();
        t = std::chrono::steady_clock::now();
        auto v7 = m::v7_adapter::from_hall_receipt(receipt, trial2, E, B);
        row["v7_adapter_incl_census_v2_s"] = since(t);
        auto moved = trial2.prepare_total_field_motion(
            {b, v7->record.record_sha256, b.accepted_field_sha256, sha("f"), sha("h"), bnd->source_sha256},
            m::v7_adapter::admit(v7));
        row["native_v7_admission_cached_census_s"] = since(t);
      }
      {
        auto trial3 = outer.fork_transport_trial();
        t = std::chrono::steady_clock::now();
        const auto &census_sha = trial3.census_sha256();
        row["native_census_v2_once_per_trial_s"] = since(t);
        auto v7 = m::v7_adapter::from_hall_receipt(receipt, trial3, E, B);
        row["v7_adapter_s"] = since(t);
        auto moved = trial3.prepare_total_field_motion(
            {b, v7->record.record_sha256, b.accepted_field_sha256, sha("f"), sha("h"), bnd->source_sha256},
            m::v7_adapter::admit(v7));
        row["native_v7_admission_cached_census_s"] = since(t);
        (void)census_sha;
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
