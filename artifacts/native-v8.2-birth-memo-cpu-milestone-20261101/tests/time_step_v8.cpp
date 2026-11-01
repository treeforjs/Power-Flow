// Native v8 per-step cost with chemistry v3 (API-008.txt 6).
// usage: time_step_v8 OUT.json MODE CARRIERS STEPS MESH_N WATER_WEIGHT NVRTC UUID
//   MODE host | device. MESH_N: Hall axial mesh transverse cells (45: 200 um,
//   180: 50 um), nz = 6; collision grid MESH_N x MESH_N x 5.
// Population (coordinator v7.1 load fixture): 49 % electrons (1-200 eV,
// isotropic), 49 % H2O+ co-located with them, 2 % H2O of weight WATER_WEIGHT,
// spread over the gap; charged weight 1. dt 10 fs.
// Each step: begin, stage 1 (prepare_chemistry_v3 + publish), fork + census,
// motion admission + accept, stage 2 + release, seal + publish. Every native
// call is timed (thread CPU = host cost, wall). The producer (straight drift
// endpoints from the census, record digests, device uploads) is not timed. In
// device mode the census is read from the device by the producer and the
// record arrays are device-resident; native performs no host row access.
#include "chem_v8_fixture.hpp"
#include "forms_gpu.hpp"
#include "hall_gpu_common.hpp"
#include "native_device.hpp"
#include "raikou/mfem_pic_full_time_device.h"
#include "raikou/json.hpp"
#include <fstream>
#include <iostream>
#include <sys/resource.h>

namespace m = raikou::mfem_backend;
using J = nlohmann::json;
int main(int argc, char **argv) {
  try {
    if (argc != 9)
      throw std::runtime_error("usage: time_step_v8 OUT.json MODE CARRIERS STEPS MESH_N WATER_WEIGHT NVRTC UUID");
    const std::string mode = argv[2];
    const bool device = mode == "device";
    const std::size_t n = std::strtoull(argv[3], nullptr, 10);
    const int steps = std::atoi(argv[4]);
    const int mesh_n = std::atoi(argv[5]);
    const double water_weight = std::strtod(argv[6], nullptr);
    std::unique_ptr<hall_gpu::Session> session;
    std::unique_ptr<pf_v7::gpu::Forms> forms;
    std::shared_ptr<m::FullTimeDevice> dev;
    if (device) {
      hall_gpu::SessionOptions so;
      so.require_uuid = argv[8];
      so.nvrtc_library = argv[7];
      session = std::make_unique<hall_gpu::Session>(so);
      forms = std::make_unique<pf_v7::gpu::Forms>(*session);
      dev = std::make_shared<m::FullTimeDevice>(*session, *forms);
    }
    const auto o = chemfix::options(mesh_n);
    auto chem = chemfix::chemistry_configuration(o);
    chem.execution = device ? m::FullTimeChemistryExecution::device : m::FullTimeChemistryExecution::host;
    chemfix::Population P;
    P.carriers = n;
    P.water_weight = water_weight;
    P.seed = 99;
    P.electron_eV_min = 1;
    P.electron_eV_max = 200;
    auto A = chemfix::make_arena(o, chem, chemfix::population(P), std::max<std::size_t>(2 * n, 4000));
    auto &arena = *A.arena;
    if (device)
      arena.attach_full_time_device(dev, m::FullTimeExecution::device);
    const pf_v7::HallAxialMesh hm{std::uint64_t(mesh_n), 6};
    std::vector<double> E(hm.ne(), 1e3), B(hm.nb(), 0.0), Jc(hm.ne(), 0.0);
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
    const std::string field_digest = pf_v7::trial_field_blocks_sha256(E, B);
    const std::string current_digest = pf_v7::current_blocks_sha256(Jc);
    const auto policy_sha = hall_driver::kinetic_policy_sha256(hall_driver::declared_kinetic_policy(o));
    const m::PicVector shift{-.5 * o.width_m, -.5 * o.gap_m, 0};
    const double dt = 1e-14;
    std::map<std::string, std::vector<chemfix::Phase>> ph;
    std::vector<double> step_cpu, step_wall;
    J per_step = J::array();
    const auto mat0 = m::FullTimeDevice::row_materialisations();
    std::uint64_t stages = 0;
    for (int k = 0; k < steps; ++k) {
      std::map<std::string, chemfix::Phase> one;
      chemfix::Timer T;
      const auto bnd = arena.full_time_boundary();
      m::FullTimeStepBinding b{bnd->accepted_step, bnd->accepted_arena_revision, bnd->time_s, bnd->time_s + dt, dt,
                               bnd->accepted_arena_sha256, bnd->field_sha256, bnd->mesh_sha256, bnd->material_sha256,
                               chemfix::sha("clock/" + std::to_string(k))};
      T.start();
      auto outer = arena.begin_full_time_step(b, A.cfg.limits.maximum_carriers, 64);
      T.stop(one["begin"]);
      m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s, .5 * dt, stages, b.accepted_field_sha256};
      T.start();
      auto s1 = outer.first_exposure(e1);
      auto t1 = s1.prepare_chemistry_v3();
      const auto st1 = t1.candidate_chemistry_stage();
      if (!s1.try_publish(t1))
        throw std::runtime_error("stage 1 publication");
      T.stop(one["stage1"]);
      T.start();
      auto trial = outer.fork_transport_trial();
      std::string census_sha;
      m::FullTimeCensusDevice cd;
      if (device) {
        cd = trial.census_v7_device();
        census_sha = cd.sha256;
      } else
        census_sha = trial.census_sha256();
      T.stop(one["fork_census"]);
      // ---- producer (untimed)
      std::vector<m::GapParticle> rows;
      if (device) {
        rows.resize(cd.carriers);
        m::device_read(*dev, cd.rows, rows.data(), rows.size() * sizeof(m::GapParticle));
      } else {
        const auto &c = trial.census_v7();
        rows.assign(c.carriers.begin(), c.carriers.end());
      }
      std::vector<m::GapParticle> ends = rows, start_rows;
      std::vector<m::FullTimeMotionPath> paths(rows.size());
      for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto &p = rows[i];
        double u2 = 0;
        for (double v : p.proper_velocity_m_s)
          u2 += v * v;
        const double gamma = std::sqrt(1 + u2 / (m::pic_c * m::pic_c));
        for (int d = 0; d < 3; ++d)
          ends[i].position_m[d] = p.position_m[d] + dt * p.proper_velocity_m_s[d] / gamma;
        auto &pa = paths[i];
        pa.native_id = p.id;
        pa.flight_s = dt;
        for (int d = 0; d < 3; ++d) {
          pa.field_start_m[d] = p.position_m[d] + shift[d];
          pa.field_displacement_m[d] = (ends[i].position_m[d] + shift[d]) - pa.field_start_m[d];
        }
      }
      start_rows = ends;
      for (std::size_t i = 0; i < ends.size(); ++i)
        for (int d = 0; d < 3; ++d)
          start_rows[i].position_m[d] = paths[i].field_start_m[d] - shift[d];
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
      x.species_sha256 = chemfix::sha("species");
      x.common_clock_sha256 = b.common_clock_sha256;
      x.kinetic_policy_sha256 = policy_sha;
      const auto K = pf_v7::kinetic_change_J({reinterpret_cast<const pf_v7::Carrier *>(start_rows.data()), start_rows.size()},
                                             {reinterpret_cast<const pf_v7::Carrier *>(ends.data()), ends.size()}, {});
      const auto W = pf_v7::field_work_J(dt, E, Jc);
      x.budget = {K.value, W.value, K.value - W.value, 0, 0, 0, 0};
      x.census_sha256 = census_sha;
      x.current_A_m_sha256 = current_digest;
      x.endpoint_charge_C_sha256 = chemfix::sha("endpoint-charge");
      x.initial_charge_C_sha256 = chemfix::sha("initial-charge");
      x.kinetic_policy = {o.width_m, o.gap_m, o.short_z_m, o.feed_z_m, o.lower_foil.left_m, o.lower_foil.right_m,
                          o.lower_foil.bottom_m, o.lower_foil.top_m, o.upper_foil.left_m, o.upper_foil.right_m,
                          o.upper_foil.bottom_m, o.upper_foil.top_m};
      x.native_to_field_m = shift;
      x.neutral_stream_dt_max_s = A.cfg.motion.neutral_stream_dt_max_s;
      x.record_sha256 = m::full_time_motion_record_v7_sha256([&] {
        auto y = x;
        y.endpoints = ends;
        y.paths = paths;
        return y;
      }());
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
      } else {
        x.endpoints = ends;
        x.paths = paths;
        x.trial_E = E;
        x.trial_B = B;
        x.current_A_m = Jc;
      }
      m::FullTimeMotionBinding mb{b, x.record_sha256, b.accepted_field_sha256, chemfix::sha("f/" + std::to_string(k)),
                                  chemfix::sha("h/" + std::to_string(k)), bnd->source_sha256};
      T.start();
      auto moved = trial.prepare_total_field_motion(mb, m::bind_total_field_motion_v7(x));
      outer.accept_transport_trial(std::move(moved));
      T.stop(one["admission_accept"]);
      m::ExposureBinding e2 = e1;
      e2.slot = m::FullTimeSlot::after_transport;
      e2.frozen_physical_stage_s = b.end_s;
      e2.owner_stage_sequence = stages + 1;
      e2.frozen_field_sha256 = mb.field_after_sha256;
      T.start();
      auto s2 = outer.second_exposure(e2);
      auto t2 = s2.prepare_chemistry_v3();
      const auto st2 = t2.candidate_chemistry_stage();
      if (!s2.try_publish(t2))
        throw std::runtime_error("stage 2 publication");
      T.stop(one["stage2"]);
      T.start();
      auto token = outer.seal({mb.field_after_sha256, mb.heat_after_sha256, bnd->source_sha256, st2->state_after_sha256});
      if (!arena.can_publish(token) || !arena.try_publish(token))
        throw std::runtime_error("publish");
      T.stop(one["seal_publish"]);
      if (device) {
        session->release(dEnds);
        session->release(dPaths);
      }
      stages += 2;
      double cpu = 0, wall = 0;
      for (auto &[name, p] : one) {
        cpu += p.cpu;
        wall += p.wall;
        ph[name].push_back(p);
      }
      step_cpu.push_back(cpu);
      step_wall.push_back(wall);
      per_step.push_back({{"cpu_ms", 1e3 * cpu}, {"wall_ms", 1e3 * wall}, {"census", rows.size()},
                          {"stage1", {{"projectiles", st1->projectiles}, {"events", st1->events},
                                      {"admitted", st1->admitted_events}, {"births", st1->births},
                                      {"device_ms", st1->device_ms}, {"device", st1->device}}},
                          {"stage2", {{"projectiles", st2->projectiles}, {"events", st2->events},
                                      {"admitted", st2->admitted_events}, {"births", st2->births},
                                      {"released", st2->released_packets}, {"device_ms", st2->device_ms},
                                      {"device", st2->device}}}});
    }
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
      phases[name] = {{"cpu_ms_per_step", 1e3 * median(after_first(c))},
                      {"wall_ms_per_step", 1e3 * median(after_first(w))}, {"first_step_cpu_ms", 1e3 * c.front()}};
    }
    rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    J out = {{"mode", mode},
             {"carriers", n},
             {"mesh_n", mesh_n},
             {"mesh_edges", hm.ne()},
             {"collision_cells", std::uint64_t(mesh_n) * mesh_n * 5},
             {"water_weight", water_weight},
             {"dt_s", dt},
             {"steps", steps},
             {"native_cpu_ms_per_step_median", 1e3 * median(after_first(step_cpu))},
             {"native_wall_ms_per_step_median", 1e3 * median(after_first(step_wall))},
             {"phases", phases},
             {"per_step", per_step},
             {"max_rss_kib", ru.ru_maxrss},
             {"host_row_materialisations", m::FullTimeDevice::row_materialisations() - mat0},
             {"source_identity", m::PicArena::source_identity()},
             {"final_accepted_arena_sha256", arena.full_time_boundary()->accepted_arena_sha256}};
    if (dev) {
      const auto cnt = dev->counters();
      out["device_admissions"] = cnt.device_admissions;
      out["host_admissions"] = cnt.host_admissions;
      out["row_uploads"] = cnt.row_uploads;
      out["row_downloads"] = cnt.row_downloads;
      out["row_download_bytes"] = cnt.row_download_bytes;
    }
    std::ofstream(argv[1]) << out.dump(2) << '\n';
    std::cout << out.dump() << '\n';
    A.arena.reset();
    if (device) {
      dev.reset();
      forms.reset();
      session->release(dE);
      session->release(dB);
      session->release(dJ);
      (void)session->close();
    }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
