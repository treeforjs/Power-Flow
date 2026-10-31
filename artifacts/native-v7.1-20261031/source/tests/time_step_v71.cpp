// Per-step native host cost, host path (v7 behaviour) vs device lane.
// usage: time_step_v71 OUT.json MODE CARRIERS STEPS MESH_N NVRTC UUID
//   MODE host | device | device_stream; MESH_N Hall axial mesh transverse cells (45: 200 um,
//   180: 50 um), nz = 6.
// Each step: begin, first exposure (no events), fork, census digest, motion
// admission, accept, second exposure, seal, publish (all native calls timed:
// thread CPU time and wall time). The producer (endpoints, record digests,
// array uploads) is not timed. Motion: straight drift at zero field (electron
// u = (0, 0, 1e6) m/s, neutrals at rest); E = 1e3 V/m, J = 0, B = 0 on the
// Hall mesh. Device mode: the record arrays and E, B, J are device-resident
// (uploaded by the producer before the timed calls) and the field/current
// digests are the v7.1 block forms; host mode uses the v7 streaming forms.
#ifndef TIMING_REFERENCE
#include "device_mode.hpp"
#endif
#include "raikou/mfem_pic_full_time.h"
#include "raikou/mfem_pic_canonical.h"
#include "pf_v7/forms.hpp"
#include "hall-canonical-driver-20261029/kinetic_boundary.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>

namespace m = raikou::mfem_backend;
using J = m::ArenaJson;
namespace {
std::string sha(const std::string &s) { return m::canonical::sha256_hex(s); }
double cpu_now() {
  timespec t;
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
  return double(t.tv_sec) + 1e-9 * double(t.tv_nsec);
}
double wall_now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
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
    out.push_back({id, electron ? "e" : "n", "manufactured/lattice", four, {x, y, z}, {0, 0, electron ? 1e6 : 0}});
  }
  return out;
}
struct Phase {
  double cpu = 0, wall = 0;
};
struct Timer {
  double c = 0, w = 0;
  void start() {
    c = cpu_now();
    w = wall_now();
  }
  void stop(Phase &p) {
    p.cpu += cpu_now() - c;
    p.wall += wall_now() - w;
  }
};
} // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 8)
      throw std::runtime_error("usage: time_step_v71 OUT.json MODE CARRIERS STEPS MESH_N NVRTC UUID");
    const std::string mode = argv[2];
    // device_stream: the device lane with the v7 streaming field/current digests (hashed by native on the host).
    const bool device = mode == "device" || mode == "device_stream", blocks = mode == "device";
    const std::size_t n = std::strtoull(argv[3], nullptr, 10);
    const int steps = std::atoi(argv[4]);
    const std::uint64_t mesh_n = std::strtoull(argv[5], nullptr, 10);
    const pf_v7::HallAxialMesh hm{mesh_n, 6};
#ifdef TIMING_REFERENCE
    if (device)
      throw std::runtime_error("the reference build (native v7) times the host path only");
    (void)argv[6];
    (void)blocks;
#else
    hall_gpu::SessionOptions so;
    so.require_uuid = argv[7];
    so.nvrtc_library = argv[6];
    std::unique_ptr<hall_gpu::Session> session;
    std::unique_ptr<pf_v7::gpu::Forms> forms;
    std::shared_ptr<m::FullTimeDevice> dev;
    if (device) {
      session = std::make_unique<hall_gpu::Session>(so);
      forms = std::make_unique<pf_v7::gpu::Forms>(*session);
      dev = std::make_shared<m::FullTimeDevice>(*session, *forms);
    }
#endif
    const auto o = options();
    m::FullTimeConfiguration c;
    c.limits.maximum_carriers = std::max<std::size_t>(n, 1000);
    c.limits.memory_budget_bytes = c.limits.maximum_carriers * m::full_time_bytes_per_carrier + m::full_time_fixed_bytes;
    c.limits.maximum_identity_entries = std::max<std::size_t>(2 * n, 1u << 22);
    m::FullTimeGenesisBinding g;
    g.time_s = 26.31715e-9;
    g.origin_provenance = "manufactured/step-timing-no-simulation";
    for (auto *x : {&g.field_sha256, &g.mesh_sha256, &g.material_sha256, &g.frame_sha256, &g.common_clock_sha256,
                    &g.heat_sha256, &g.source_sha256, &g.rng_sha256})
      *x = sha("timing/" + std::to_string(std::uintptr_t(x) - std::uintptr_t(&g)));
    const std::vector<m::ArenaStateSpec> states = {{"e", 0, 0, -1, m::pic_me, 0, 0}, {"n", 2, 1, 0, 2.9915e-26, 0, 0}};
    const std::vector<m::ArenaChannel> channels = {{"conditional-storage-control", "e", "n", {"e", "n"}, sha("data"),
                                                    sha("kernel"), sha("binding"), "manufactured-storage-control-only"}};
    auto arena = m::PicArena::create_full_time(o, states, channels, lattice(n), g, c);
    // Producer copy of the census rows (ascending native ID), kept by the
    // producer from step to step (not read back from native).
    std::vector<m::GapParticle> rows;
    for (const auto &r : arena->rate_population(n))
      rows.push_back(r.native);
    std::sort(rows.begin(), rows.end(), [](const auto &a, const auto &x) { return a.id < x.id; });
#ifndef TIMING_REFERENCE
    if (device)
      arena->attach_full_time_device(dev, m::FullTimeExecution::device);
#endif
    // Field arrays (host and device).
    std::vector<double> E(hm.ne(), 1e3), B(hm.nb(), 0.0), Jc(hm.ne(), 0.0);
#ifndef TIMING_REFERENCE
    auto up = [&](const void *p, std::size_t bytes) {
      auto d = session->alloc(bytes ? bytes : 8, "timing-array");
      if (bytes) {
        auto h = session->alloc_pinned(bytes, "timing-stage");
        std::memcpy(h.pointer, p, bytes);
        session->upload(d, 0, h, 0, bytes);
        session->synchronize("timing upload");
        session->release(h);
      }
      return d;
    };
    hall_gpu::DeviceBuffer dE, dB, dJ;
    if (device) {
      dE = up(E.data(), E.size() * 8);
      dB = up(B.data(), B.size() * 8);
      dJ = up(Jc.data(), Jc.size() * 8);
    }
    const std::string field_digest = blocks ? pf_v7::trial_field_blocks_sha256(E, B) : pf_v7::total_field_sha256(E, B);
    const std::string current_digest = blocks ? pf_v7::current_blocks_sha256(Jc) : pf_v7::array_sha256(Jc);
#else
    const std::string field_digest = pf_v7::total_field_sha256(E, B);
    const std::string current_digest = pf_v7::array_sha256(Jc);
#endif

    const auto policy_sha = hall_driver::kinetic_policy_sha256(hall_driver::declared_kinetic_policy(o));
    const m::PicVector shift{-.5 * o.width_m, -.5 * o.gap_m, 0};
    std::map<std::string, std::vector<Phase>> ph;
    Timer T;
    std::vector<double> step_cpu, step_wall;
#ifndef TIMING_REFERENCE
    const auto mat0 = m::FullTimeDevice::row_materialisations();
#endif
    for (int k = 0; k < steps; ++k) {
      double cpu0 = 0, wall0 = 0;
      std::map<std::string, Phase> one;
      const auto bnd = arena->full_time_boundary();
      const double dt = 5e-15;
      m::FullTimeStepBinding b{bnd->accepted_step, bnd->accepted_arena_revision, bnd->time_s, bnd->time_s + dt, dt,
                               bnd->accepted_arena_sha256, bnd->field_sha256, bnd->mesh_sha256, bnd->material_sha256,
                               sha("clock/" + std::to_string(k))};
      auto completion = [&](const m::ExposureBinding &e, const std::string &before, const std::string &after) {
        m::ArenaStochasticCompletion x;
        const auto ch = channels[0];
        x.integration_owner_id = "timing";
        x.before_rng_clock_sha256 = before;
        x.after_rng_clock_sha256 = after;
        x.common_clock_sha256 = b.common_clock_sha256;
        x.channel_id = ch.id;
        x.dataset_sha256 = ch.dataset_sha256;
        x.kernel_manifest_sha256 = ch.kernel_manifest_sha256;
        x.kernel_binding_file_sha256 = ch.kernel_binding_file_sha256;
        x.model = ch.model;
        x.stage_sequence = e.owner_stage_sequence;
        x.frozen_physical_stage_s = e.frozen_physical_stage_s;
        x.common_interval_begin_s = b.begin_s;
        x.common_interval_end_s = b.end_s;
        return x;
      };
      const std::string rng0 = bnd->rng_sha256, rng1 = sha("rng1/" + std::to_string(k)), rng2 = sha("rng2/" + std::to_string(k));
      T.start();
      auto outer = arena->begin_full_time_step(b, n, 32);
      T.stop(one["begin"]);
      m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s, .5 * dt, 2u * std::uint64_t(k), b.accepted_field_sha256};
      T.start();
      auto s1 = outer.first_exposure(e1);
      auto t1 = s1.prepare_compound_with_completion(1, [&](m::CollisionScope &) { return completion(e1, rng0, rng1); });
      if (!s1.try_publish(t1))
        throw std::runtime_error("first exposure");
      auto trial = outer.fork_transport_trial();
      T.stop(one["first_exposure_fork"]);
      T.start();
      const std::string census_sha = trial.census_sha256();
      T.stop(one["census"]);
      // ---- producer (untimed): straight drift record.
      std::vector<m::GapParticle> ends;
      std::vector<m::FullTimeMotionPath> paths;
      {
        ends = rows;
        paths.resize(rows.size());
        for (std::size_t i = 0; i < rows.size(); ++i) {
          const auto &p = rows[i];
          double u2 = 0;
          for (double v : p.proper_velocity_m_s)
            u2 += v * v;
          const double gamma = std::sqrt(1 + u2 / (m::pic_c * m::pic_c));
          auto &e = ends[i];
          for (int d = 0; d < 3; ++d)
            e.position_m[d] = p.position_m[d] + dt * p.proper_velocity_m_s[d] / gamma;
          auto &pa = paths[i];
          pa.native_id = p.id;
          pa.flight_s = dt;
          for (int d = 0; d < 3; ++d) {
            pa.field_start_m[d] = p.position_m[d] + shift[d];
            pa.field_displacement_m[d] = (e.position_m[d] + shift[d]) - pa.field_start_m[d];
          }
        }
      }
      m::FullTimeMotionRecordV7 x;
      x.accepted_step = b.accepted_step;
      x.staged_arena_revision = trial.staged_arena_revision();
      x.begin_s = b.begin_s;
      x.end_s = b.end_s;
      x.dt_s = dt;
      x.staged_arena_sha256 = trial.state_prefix();
      x.trial_field_sha256 = field_digest;
      x.mesh_sha256 = b.mesh_sha256;
      x.material_sha256 = b.material_sha256;
      x.species_sha256 = sha("species");
      x.common_clock_sha256 = b.common_clock_sha256;
      x.kinetic_policy_sha256 = policy_sha;
      std::vector<m::GapParticle> start_rows(ends.size());
      for (std::size_t i = 0; i < ends.size(); ++i) {
        start_rows[i] = ends[i];
        for (int d = 0; d < 3; ++d)
          start_rows[i].position_m[d] = paths[i].field_start_m[d] - shift[d];
      }
      const auto K = pf_v7::kinetic_change_J({reinterpret_cast<const pf_v7::Carrier *>(start_rows.data()), start_rows.size()},
                                             {reinterpret_cast<const pf_v7::Carrier *>(ends.data()), ends.size()}, {});
      const auto W = pf_v7::field_work_J(dt, E, Jc);
      x.budget = {K.value, W.value, K.value - W.value, 0, 0, 0, 0};
      x.census_sha256 = census_sha;
      x.current_A_m_sha256 = current_digest;
      x.endpoint_charge_C_sha256 = sha("endpoint-charge");
      x.initial_charge_C_sha256 = sha("initial-charge");
      x.kinetic_policy = {o.width_m, o.gap_m, o.short_z_m, o.feed_z_m, o.lower_foil.left_m, o.lower_foil.right_m,
                          o.lower_foil.bottom_m, o.lower_foil.top_m, o.upper_foil.left_m, o.upper_foil.right_m,
                          o.upper_foil.bottom_m, o.upper_foil.top_m};
      x.native_to_field_m = shift;
      x.record_sha256 = m::full_time_motion_record_v7_sha256([&] {
        auto y = x;
        y.endpoints = ends;
        y.paths = paths;
        return y;
      }());
#ifndef TIMING_REFERENCE
      hall_gpu::DeviceBuffer dEnds, dPaths;
      if (device) {
        dEnds = up(ends.data(), ends.size() * sizeof(m::GapParticle));
        dPaths = up(paths.data(), paths.size() * sizeof(m::FullTimeMotionPath));
        x.device.endpoints = dEnds.address;
        x.device.n_endpoints = ends.size();
        x.device.paths = dPaths.address;
        x.device.n_paths = paths.size();
        x.device.E = dE.address;
        x.device.B = dB.address;
        x.device.J = dJ.address;
        x.device.n_E = E.size();
        x.device.n_B = B.size();
      } else
#endif
      {
        x.endpoints = ends;
        x.paths = paths;
        x.trial_E = E;
        x.trial_B = B;
        x.current_A_m = Jc;
      }
      m::FullTimeMotionBinding mb{b, x.record_sha256, b.accepted_field_sha256, sha("f/" + std::to_string(k)),
                                  sha("h/" + std::to_string(k)), bnd->source_sha256};
      T.start();
      auto moved = trial.prepare_total_field_motion(mb, m::bind_total_field_motion_v7(x));
      outer.accept_transport_trial(std::move(moved));
      T.stop(one["admission_accept"]);
      m::ExposureBinding e2 = e1;
      e2.slot = m::FullTimeSlot::after_transport;
      e2.frozen_physical_stage_s = b.end_s;
      e2.owner_stage_sequence = 2u * std::uint64_t(k) + 1;
      e2.frozen_field_sha256 = mb.field_after_sha256;
      T.start();
      auto s2 = outer.second_exposure(e2);
      auto t2 = s2.prepare_compound_with_completion(1, [&](m::CollisionScope &) { return completion(e2, rng1, rng2); });
      if (!s2.try_publish(t2))
        throw std::runtime_error("second exposure");
      T.stop(one["second_exposure"]);
      T.start();
      auto token = outer.seal({mb.field_after_sha256, mb.heat_after_sha256, bnd->source_sha256, rng2});
      if (!arena->can_publish(token) || !arena->try_publish(token))
        throw std::runtime_error("publish");
      T.stop(one["seal_publish"]);
#ifndef TIMING_REFERENCE
      if (device) {
        session->release(dEnds);
        session->release(dPaths);
      }
#endif
      rows = ends;
      for (auto &[name, p] : one) {
        cpu0 += p.cpu;
        wall0 += p.wall;
        ph[name].push_back(p);
      }
      step_cpu.push_back(cpu0);
      step_wall.push_back(wall0);
    }
    // Steady state: medians over the steps after the first (the first step
    // also validates the created state and uploads the rows in device mode).
    auto after_first = [](std::vector<double> v) {
      if (v.size() > 1)
        v.erase(v.begin());
      return v;
    };
    auto median = [](std::vector<double> v) {
      std::sort(v.begin(), v.end());
      return v.empty() ? 0.0 : v[v.size() / 2];
    };
    J phases = J::object();
    for (const auto &[name, v] : ph) {
      std::vector<double> c, w;
      for (const auto &p : v) {
        c.push_back(p.cpu);
        w.push_back(p.wall);
      }
      phases[name] = {{"cpu_ms_per_step", 1e3 * median(after_first(c))}, {"wall_ms_per_step", 1e3 * median(after_first(w))},
                      {"first_step_cpu_ms", 1e3 * c.front()}};
    }
    std::vector<double> step_ms;
    for (double t : step_cpu)
      step_ms.push_back(1e3 * t);
    J out = {{"mode", mode},
             {"carriers", n},
             {"mesh_n", mesh_n},
             {"mesh_edges", hm.ne()},
             {"steps", steps},
             {"native_cpu_ms_per_step_median", 1e3 * median(after_first(step_cpu))},
             {"native_wall_ms_per_step_median", 1e3 * median(after_first(step_wall))},
             {"native_cpu_us_per_carrier", 1e6 * median(after_first(step_cpu)) / double(n)},
             {"phases", phases},
             {"step_cpu_ms", step_ms},
             {"source_identity", m::PicArena::source_identity()},
             {"final_accepted_arena_sha256", arena->full_time_boundary()->accepted_arena_sha256}};
#ifndef TIMING_REFERENCE
    out["host_row_materialisations"] = m::FullTimeDevice::row_materialisations() - mat0;
    if (dev) {
      const auto cnt = dev->counters();
      out["device_admissions"] = cnt.device_admissions;
      out["host_admissions"] = cnt.host_admissions;
      out["row_uploads"] = cnt.row_uploads;
      out["row_downloads"] = cnt.row_downloads;
      out["array_uploads"] = cnt.array_uploads;
      out["array_downloads"] = cnt.array_downloads;
      out["array_download_bytes"] = cnt.array_download_bytes;
      out["streaming_field_digests"] = cnt.streaming_field_digests;
    }
#endif
    std::ofstream(argv[1]) << out.dump(2) << '\n';
    std::cout << out.dump() << '\n';
    arena.reset();
#ifndef TIMING_REFERENCE
    if (device) {
      dev.reset();
      forms.reset();
      session->release(dE);
      session->release(dB);
      session->release(dJ);
      (void)session->close();
    }
#endif
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
