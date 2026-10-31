// test_forms_gpu.cpp -- CPU vs GPU bit identity of the native v7 forms.
// Usage: test_forms_gpu NVRTC_LIBRARY DEVICE_UUID [max_n]
// Prints RESULT lines (parity, device timings) and the session receipts.
#include "datasets.hpp"
#include "forms_gpu.hpp"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <ctime>
#include <numeric>

using namespace v7test;
namespace g7 = pf_v7::gpu;
namespace {
bool same(const g7::SumOut &d, const SumValue &c, const ExactSum::Digits *digits, bool scaled = false) {
  bool ok = bits(scaled ? d.scaled : d.value) == bits(c.value) && d.status == c.status;
  if (digits)
    ok = ok && std::memcmp(d.digit, digits->data(), sizeof d.digit) == 0;
  return ok;
}
ExactSum::Digits digits_of(std::span<const double> a, std::span<const double> b) {
  ExactSum s;
  for (std::size_t i = 0; i < a.size(); ++i)
    s.add_product(a[i], b[i]);
  return s.digits();
}
double ms_since(std::chrono::steady_clock::time_point t) {
  return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count();
}
} // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: test_forms_gpu NVRTC_LIBRARY DEVICE_UUID [max_n]\n");
    return 2;
  }
  const std::size_t max_n = argc > 3 ? std::strtoull(argv[3], nullptr, 10) : 1000000;
  int rc = 0;
  {
    hall_gpu::SessionOptions so;
    so.require_uuid = argv[2];
    so.nvrtc_library = argv[1];
    so.check_each_launch = std::getenv("PF7_CHECK_EACH_LAUNCH") != nullptr;
    hall_gpu::Session s(so);
    g7::Forms forms(s);
    g7::HostRunner run(s, forms);
    std::printf("DEVICE %s\n", hall_gpu::device_json(s.device()).c_str());
    std::printf("COMPILE %s\n", hall_gpu::compile_json(forms.compile_record()).c_str());
    std::printf("RESULT kernel_source_sha256 %s\n", g7::Forms::kernel_source_sha256().c_str());
    const unsigned grids[] = {1, 2, 5, 8, 64, 240};

    // ---- (a) adversarial and random product sets, every grid size.
    {
      auto sets = adversarial(11);
      for (std::size_t n = 10000; n <= max_n; n *= 10)
        sets.push_back(random_cancel(n, 900 + n));
      for (const auto &set : sets) {
        const auto ref = exact_dot(set.a, set.b);
        const auto dig = digits_of(set.a, set.b);
        int bad = 0;
        double t = 0;
        for (unsigned gsz : grids) {
          const auto d = run.sum_products(set.a, set.b, 1.0, gsz);
          bad += !same(d, ref, &dig);
          if (gsz == 64)
            t = run.last_device_ms;
        }
        check(bad == 0, "GPU exact dot = CPU (" + set.name + ") on grids 1/2/5/8/64/240");
        std::printf("RESULT gpu_dot %s terms=%zu mismatches=%d device_ms_g64=%.4f\n", set.name.c_str(), set.a.size(),
                    bad, t);
      }
      // Field work over the default mesh edge count: dt x rounded sum.
      std::vector<double> E(41676), J(41676);
      std::mt19937_64 g(3);
      std::normal_distribution<double> N(0, 1);
      for (std::size_t e = 0; e < E.size(); ++e) {
        E[e] = 1e7 * N(g);
        J[e] = 1e3 * N(g);
      }
      const auto w = field_work_J(1e-14, E, J);
      const auto d = run.sum_products(E, J, 1e-14, 64);
      check(same(d, w, nullptr, true) && bits(d.value) == bits(w.sum), "GPU field work = CPU (41,676 edges)");
      std::printf("RESULT gpu_field_work edges=41676 value=%.17g identical=%d device_ms=%.4f\n", d.scaled,
                  int(same(d, w, nullptr, true)), run.last_device_ms);
    }

    // ---- (a)-(c) per motion size.
    for (std::size_t n = 10000; n <= max_n; n *= 10) {
      const auto m = motion(n, 31 + n);
      // (a) kinetic: aligned and split forms, every grid, = CPU.
      const auto ref = kinetic_change_J(m.census, m.endpoints, m.sinks);
      ExactSum dsum;
      for (std::size_t i = 0; i < m.census.size(); ++i)
        dsum.add_product(m.census[i].weight, pf7::kinetic_change(m.census[i].mass_kg, m.census[i].proper_velocity_m_s,
                                                                 &m.u_end[3 * i]));
      const auto dig = dsum.digits();
      int bad = 0;
      double t_aligned = 0, t_split = 0;
      for (unsigned gsz : grids) {
        const auto a = run.kinetic_aligned(m.census, m.u_end, gsz);
        if (gsz == 64)
          t_aligned = run.last_device_ms;
        const auto b = run.kinetic_split(m.census, m.endpoints, m.sinks, gsz);
        if (gsz == 64)
          t_split = run.last_device_ms;
        bad += !same(a, ref, &dig) + !same(b, ref, &dig);
      }
      check(bad == 0, "GPU kinetic (aligned and split) = CPU at n=" + std::to_string(n));
      std::printf("RESULT gpu_kinetic n=%zu value=%.17g mismatches=%d device_ms_aligned=%.4f device_ms_split=%.4f\n",
                  n, ref.value, bad, t_aligned, t_split);
      // (b) block lists of every section.
      struct Sec {
        const char *name;
        std::uint32_t s;
        const void *p;
        std::size_t n;
      };
      for (const auto &x : {Sec{"endpoints", pf7::section_endpoints, m.endpoints.data(), m.endpoints.size()},
                            Sec{"paths", pf7::section_paths, m.paths.data(), m.paths.size()},
                            Sec{"sinks", pf7::section_sinks, m.sinks.data(), m.sinks.size()},
                            Sec{"births", pf7::section_births, m.births.data(), m.births.size()},
                            Sec{"census", pf7::section_census, m.census_records.data(), m.census_records.size()}}) {
        const auto c = block_list(x.s, x.p, x.n);
        const auto t0 = std::chrono::steady_clock::now();
        const auto d = run.block_list(x.s, x.p, x.n);
        const double wall = ms_since(t0);
        const bool ok = d.root == c.root_digest && d.leaves == c.leaves;
        check(ok, std::string("GPU block list = CPU: ") + x.name + " n=" + std::to_string(x.n));
        std::printf("RESULT gpu_block_list section=%s records=%zu leaves=%zu identical=%d device_ms=%.4f wall_ms=%.3f root=%s\n",
                    x.name, x.n, c.leaves.size(), int(ok), run.last_device_ms, wall, c.root.c_str());
      }
      // Record digest from device section roots (host O(1) header).
      {
        MotionHeaderV7 h;
        h.census_sha256 = census_v2_sha256(m.census_records, m.states);
        const auto cpu_sections = motion_sections_v7(m.endpoints, m.paths, m.sinks, m.births);
        MotionSectionsV7 dev = cpu_sections;
        dev.endpoints_root = hex(run.block_list(pf7::section_endpoints, m.endpoints.data(), m.endpoints.size()).root);
        dev.paths_root = hex(run.block_list(pf7::section_paths, m.paths.data(), m.paths.size()).root);
        dev.sinks_root = hex(run.block_list(pf7::section_sinks, m.sinks.data(), m.sinks.size()).root);
        dev.births_root = hex(run.block_list(pf7::section_births, m.births.data(), m.births.size()).root);
        check(motion_record_v7_sha256(h, dev) == motion_record_v7_sha256(h, cpu_sections),
              "record digest from device section roots = CPU at n=" + std::to_string(n));
      }
      // (c) census check (and tamper), order, digest.
      {
        const auto k = run.census_check(m.census_records, m.states.size());
        const double t_check = run.last_device_ms;
        check(k.ok(), "GPU census check accepts the canonical census n=" + std::to_string(n));
        int tbad = 0;
        auto variant = [&](std::vector<CensusRecord> c) {
          const auto a = check_census(c, m.states.size());
          const auto b = run.census_check(c, m.states.size());
          tbad += a.code != b.code || a.index != b.index || a.ok();
        };
        auto c = m.census_records;
        std::swap(c[n / 2], c[n / 2 + 1]);
        variant(c);
        c = m.census_records;
        c.insert(c.begin() + n / 3, c[n / 3]);
        variant(c);
        c = m.census_records;
        c[1].carrier.proper_velocity_m_s[2] = std::nextafter(c[1].carrier.proper_velocity_m_s[2], 1e300);
        variant(c);
        c = m.census_records;
        c[n - 1].state_index = 99;
        variant(c);
        check(tbad == 0, "GPU census refuses reordered/duplicated/divergent/state tampering as the CPU (same index, code)");
        auto shuffled = m.census_records;
        std::mt19937_64 g(n);
        std::shuffle(shuffled.begin(), shuffled.end(), g);
        const auto ordered = run.census_order(shuffled);
        const double t_order = run.last_device_ms;
        const bool order_ok = ordered.size() == m.census_records.size() &&
                              std::memcmp(ordered.data(), m.census_records.data(), ordered.size() * sizeof(CensusRecord)) == 0;
        check(order_ok, "GPU canonical census order of a shuffled census n=" + std::to_string(n));
        const auto d = run.block_list(pf7::section_census, m.census_records.data(), m.census_records.size());
        const auto cpu = census_v2_sha256(m.census_records, m.states);
        const auto dev = census_v2_sha256(n, census_states_sha256(m.states), hex(d.root));
        check(dev == cpu, "census digest from device carriers root = CPU n=" + std::to_string(n));
        std::printf("RESULT gpu_census n=%zu check_ms=%.4f order_ms=%.4f hash_ms=%.4f tamper_mismatches=%d digest=%s\n", n,
                    t_check, t_order, run.last_device_ms, tbad, cpu.c_str());
      }
    }
    // Host cost of the device path with device-resident inputs (the Hall v3.2
    // situation): census check + census root, K (aligned), W (41,676 edges),
    // four section lists, one download of the results, host header digest.
    // Host wall and host CPU time (blocking sync yields the core).
    for (std::size_t n = 10000; n <= max_n; n *= 10) {
      const auto m = motion(n, 7 + n);
      std::vector<double> E(41676, 1e5), J(41676, 1e-3);
      auto up = [&](const void *p, std::size_t bytes, const char *tag) {
        auto d = s.alloc(bytes ? bytes : 8, tag);
        if (bytes) {
          auto h = s.alloc_pinned(bytes, tag);
          std::memcpy(h.pointer, p, bytes);
          s.upload(d, 0, h, 0, bytes);
          s.synchronize("resident upload");
          s.release(h);
        }
        return d;
      };
      auto census = up(m.census_records.data(), m.census_records.size() * sizeof(CensusRecord), "census");
      auto ue = up(m.u_end.data(), m.u_end.size() * 8, "u-end");
      auto ep = up(m.endpoints.data(), m.endpoints.size() * sizeof(Carrier), "endpoints");
      auto pa = up(m.paths.data(), m.paths.size() * sizeof(Path), "paths");
      auto sk = up(m.sinks.data(), m.sinks.size() * sizeof(Sink), "sinks");
      auto bi = up(m.births.data(), m.births.size() * sizeof(Birth), "births");
      auto dE = up(E.data(), E.size() * 8, "E");
      auto dJ = up(J.data(), J.size() * 8, "J");
      const unsigned blocks = 120;
      std::size_t bl = 0;
      for (std::size_t c : {m.census_records.size(), m.endpoints.size(), m.paths.size(), m.sinks.size(), m.births.size()})
        bl = std::max(bl, g7::Forms::block_list_scratch_bytes(c));
      auto scratch_sum = s.alloc(2 * g7::Forms::sum_scratch_bytes(blocks), "sum-scratch");
      auto scratch_census = s.alloc(g7::Forms::census_check_scratch_bytes(n), "census-scratch");
      std::vector<hall_gpu::DeviceBuffer> scratch_bl;
      for (int k = 0; k < 5; ++k)
        scratch_bl.push_back(s.alloc(bl, "block-scratch"));
      auto results = s.alloc(4096, "results");
      auto host = s.alloc_pinned(4096, "results-host");
      std::vector<double> wall, cpu, dev;
      std::string record;
      for (int rep = 0; rep < 21; ++rep) {
        timespec c0, c1;
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &c0);
        const auto t0 = std::chrono::steady_clock::now();
        s.clear_markers();
        const int k0 = s.mark("v7-forms");
        const auto R = results.address;
        forms.census_check(n, census.address, m.states.size(), scratch_census.address, R);
        forms.block_list(pf7::section_census, n, census.address, scratch_bl[4].address, R + 16);
        const auto cb = census.address;
        forms.kinetic_aligned(n, cb + 72, 12, cb + 64, 12, cb + 32, 12, ue.address, 3, blocks, scratch_sum.address, R + 64);
        forms.sum_products(E.size(), dE.address, 1, dJ.address, 1, 1e-14, blocks,
                           scratch_sum.address + g7::Forms::sum_scratch_bytes(blocks), R + 64 + sizeof(g7::SumOut));
        const std::uint64_t roots = R + 64 + 2 * sizeof(g7::SumOut);
        forms.block_list(pf7::section_endpoints, m.endpoints.size(), ep.address, scratch_bl[0].address, roots);
        forms.block_list(pf7::section_paths, m.paths.size(), pa.address, scratch_bl[1].address, roots + 32);
        forms.block_list(pf7::section_sinks, m.sinks.size(), sk.address, scratch_bl[2].address, roots + 64);
        forms.block_list(pf7::section_births, m.births.size(), bi.address, scratch_bl[3].address, roots + 96);
        s.download(host, 0, results, 0, 64 + 2 * sizeof(g7::SumOut) + 128);
        const int k1 = s.mark("v7-forms-end");
        s.synchronize("v7 forms resident");
        const auto *base = static_cast<const unsigned char *>(host.pointer);
        std::uint64_t fail[2];
        std::memcpy(fail, base, 16);
        g7::SumOut K, W;
        std::memcpy(&K, base + 64, sizeof K);
        std::memcpy(&W, base + 64 + sizeof K, sizeof W);
        std::uint32_t words[5][8];
        std::memcpy(words[4], base + 16, 32);
        std::memcpy(words, base + 64 + 2 * sizeof K, 4 * 32);
        MotionHeaderV7 h;
        h.budget = {K.value, W.scaled, K.value - W.scaled, 0, 0, 0, 0};
        h.census_sha256 = census_v2_sha256(n, census_states_sha256(m.states), hex(digest_of_words(words[4])));
        MotionSectionsV7 sec;
        sec.endpoints = m.endpoints.size();
        sec.paths = m.paths.size();
        sec.sinks = m.sinks.size();
        sec.births = m.births.size();
        sec.endpoints_root = hex(digest_of_words(words[0]));
        sec.paths_root = hex(digest_of_words(words[1]));
        sec.sinks_root = hex(digest_of_words(words[2]));
        sec.births_root = hex(digest_of_words(words[3]));
        record = motion_record_v7_sha256(h, sec);
        const auto t1 = std::chrono::steady_clock::now();
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &c1);
        if (fail[0] != ~0ull || K.status || W.status)
          check(false, "resident census/sums refused");
        if (rep) {
          wall.push_back(std::chrono::duration<double, std::milli>(t1 - t0).count());
          cpu.push_back(1e3 * (double(c1.tv_sec - c0.tv_sec) + 1e-9 * double(c1.tv_nsec - c0.tv_nsec)));
          dev.push_back(s.elapsed_ms(k0, k1));
        }
      }
      // Reference: the same record on the CPU.
      const auto Kc = kinetic_change_J_aligned(m.census, m.u_end);
      const auto Wc = field_work_J(1e-14, E, J);
      MotionHeaderV7 hc;
      hc.budget = {Kc.value, Wc.value, Kc.value - Wc.value, 0, 0, 0, 0};
      hc.census_sha256 = census_v2_sha256(m.census_records, m.states);
      check(record == motion_record_v7_sha256(hc, motion_sections_v7(m.endpoints, m.paths, m.sinks, m.births)),
            "resident device production = CPU record digest at n=" + std::to_string(n));
      auto med = [](std::vector<double> v) {
        std::sort(v.begin(), v.end());
        return v[v.size() / 2];
      };
      std::printf("RESULT resident_host n=%zu host_wall_ms=%.4f host_cpu_ms=%.4f device_ms=%.4f host_cpu_us_per_carrier=%.6f\n",
                  n, med(wall), med(cpu), med(dev), 1e3 * med(cpu) / double(n));
      for (auto *d : {&census, &ue, &ep, &pa, &sk, &bi, &dE, &dJ, &scratch_sum, &scratch_census, &results})
        s.release(*d);
      for (auto &d : scratch_bl)
        s.release(d);
      s.release(host);
    }
    // Determinism: repeat one large case three times in process.
    {
      const auto m = motion(std::min<std::size_t>(max_n, 100000), 99);
      std::string first;
      int bad = 0;
      for (int r = 0; r < 3; ++r) {
        const auto a = run.kinetic_split(m.census, m.endpoints, m.sinks, 64);
        const auto b = run.block_list(pf7::section_paths, m.paths.data(), m.paths.size());
        const std::string key = std::to_string(bits(a.value)) + hex(b.root);
        if (r == 0)
          first = key;
        bad += key != first;
      }
      check(bad == 0, "repeated device runs identical");
    }
    const auto mem = s.memory();
    std::printf("MEMORY %s\n", hall_gpu::memory_json(mem).c_str());
    const auto td = s.close();
    std::printf("TEARDOWN %s\n", hall_gpu::teardown_json(td).c_str());
    if (!td.released || td.quarantined)
      rc = 3;
  }
  std::printf("RESULT gpu_checks=%d failures=%d\n", checks, failures);
  return failures ? 1 : rc;
}
