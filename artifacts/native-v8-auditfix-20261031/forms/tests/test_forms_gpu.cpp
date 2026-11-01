// test_forms_gpu.cpp -- CPU vs GPU bit identity of the native v7 forms.
// Usage: test_forms_gpu NVRTC_LIBRARY DEVICE_UUID [max_n]
// Prints RESULT lines (parity, device timings) and the session receipts.
#include "datasets.hpp"
#include "forms_gpu.hpp"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <ctime>
#include <openssl/sha.h>
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
    // ---- native v7.1 additions (API-071.txt): fused block lists, v5 ID-block
    // roots, |E||J| exact sum, E^, streaming and block-form field digests.
    {
      auto up = [&](const void *p, std::size_t bytes) {
        auto d = s.alloc(bytes ? bytes : 8, "v71-test");
        if (bytes) {
          auto h = s.alloc_pinned(bytes, "v71-stage");
          std::memcpy(h.pointer, p, bytes);
          s.upload(d, 0, h, 0, bytes);
          s.synchronize("v71 upload");
          s.release(h);
        }
        return d;
      };
      auto down = [&](void *out, const hall_gpu::DeviceBuffer &d, std::size_t off, std::size_t bytes) {
        auto h = s.alloc_pinned(bytes, "v71-down");
        s.download(h, 0, d, off, bytes);
        s.synchronize("v71 download");
        std::memcpy(out, h.pointer, bytes);
        s.release(h);
      };
      for (std::size_t n = 10000; n <= max_n; n *= 10) {
        const auto m = motion(n, 71 + n);
        // Sparse IDs (v5 ID blocks with partial blocks) are in the generated census.
        auto dEnd = up(m.endpoints.data(), m.endpoints.size() * 80);
        auto dPath = up(m.paths.data(), m.paths.size() * 64);
        auto dSink = up(m.sinks.data(), m.sinks.size() * sizeof(Sink));
        auto dBirth = up(m.births.data(), m.births.size() * 24);
        auto dCen = up(m.census_records.data(), m.census_records.size() * 96);
        std::vector<g7::Forms::BlockListJob> jobs = {{pf7::section_endpoints, m.endpoints.size(), dEnd.address, 0},
                                                     {pf7::section_paths, m.paths.size(), dPath.address, 0},
                                                     {pf7::section_sinks, m.sinks.size(), dSink.address, 0},
                                                     {pf7::section_births, m.births.size(), dBirth.address, 0},
                                                     {pf7::section_census, m.census_records.size(), dCen.address, 0}};
        auto roots = s.alloc(jobs.size() * 32, "v71-roots");
        for (std::size_t k = 0; k < jobs.size(); ++k)
          jobs[k].root_out = roots.address + 32 * k;
        auto scr = s.alloc(g7::Forms::block_lists_scratch_bytes(jobs), "v71-scratch");
        auto sep = s.alloc(g7::Forms::block_list_scratch_bytes(std::max<std::size_t>(n, 1)), "v71-sep");
        auto sroots = s.alloc(jobs.size() * 32, "v71-sroots");
        s.synchronize("v71 alloc");
        // Three alternating repetitions (warm device), fastest of each.
        double fused_ms = 1e300, separate_ms = 1e300;
        for (int rep = 0; rep < 3; ++rep) {
          s.clear_markers();
          const int a0 = s.mark("fused");
          forms.block_lists(jobs, scr.address);
          const int a1 = s.mark("fused-end");
          const int b0 = s.mark("separate");
          for (std::size_t k = 0; k < jobs.size(); ++k)
            forms.block_list(jobs[k].section, jobs[k].n, jobs[k].records, sep.address, sroots.address + 32 * k);
          const int b1 = s.mark("separate-end");
          s.synchronize("v71 lists");
          fused_ms = std::min(fused_ms, s.elapsed_ms(a0, a1));
          separate_ms = std::min(separate_ms, s.elapsed_ms(b0, b1));
        }
        std::vector<unsigned char> f(jobs.size() * 32), q(jobs.size() * 32);
        down(f.data(), roots, 0, f.size());
        down(q.data(), sroots, 0, q.size());
        const Digest cpu_roots[5] = {block_list(pf7::section_endpoints, m.endpoints.data(), m.endpoints.size()).root_digest,
                                     block_list(pf7::section_paths, m.paths.data(), m.paths.size()).root_digest,
                                     block_list(pf7::section_sinks, m.sinks.data(), m.sinks.size()).root_digest,
                                     block_list(pf7::section_births, m.births.data(), m.births.size()).root_digest,
                                     block_list(pf7::section_census, m.census_records.data(), m.census_records.size()).root_digest};
        bool fused_ok = f == q;
        for (int k = 0; k < 5; ++k) {
          std::uint32_t w[8];
          std::memcpy(w, f.data() + 32 * k, 32);
          fused_ok = fused_ok && digest_of_words(w) == cpu_roots[k];
        }
        check(fused_ok, "fused block lists = separate = CPU at n=" + std::to_string(n));
        std::printf("RESULT gpu_fused_block_lists n=%zu sections=5 identical=%d fused_ms=%.4f separate_ms=%.4f\n", n,
                    int(fused_ok), fused_ms, separate_ms);
        // v5 ID-block roots: carrier rows (endpoints) and path rows.
        auto ids = s.alloc(g7::Forms::id_root_scratch_bytes(n), "v71-ids");
        auto idr = s.alloc(64, "v71-idr");
        s.synchronize("v71 alloc");
        s.clear_markers();
        const int c0 = s.mark("id");
        forms.id_root(0, m.endpoints.size(), dEnd.address, 10, ids.address, idr.address);
        forms.id_root(1, m.paths.size(), dPath.address, 8, ids.address, idr.address + 32);
        const int c1 = s.mark("id-end");
        s.synchronize("v71 id");
        unsigned char ir[64];
        down(ir, idr, 0, 64);
        // CPU v5 roots (canonical carrier_rows_digest form, written out here).
        auto v5_root = [&](const char *leaf_tag, const char *root_tag, const unsigned char *rec, std::size_t nrec,
                           std::size_t R) {
          std::vector<Digest> leaves;
          std::size_t i = 0;
          auto id_at = [&](std::size_t k) {
            std::uint64_t v;
            std::memcpy(&v, rec + k * R, 8);
            return v;
          };
          while (i < nrec) {
            const auto blk = (id_at(i) - 1) / 1024;
            std::size_t j = i;
            while (j < nrec && (id_at(j) - 1) / 1024 == blk)
              ++j;
            std::vector<unsigned char> msg;
            auto u64 = [&](std::uint64_t v) {
              for (int b = 0; b < 8; ++b)
                msg.push_back(static_cast<unsigned char>(v >> (8 * b)));
            };
            const std::string t = leaf_tag;
            u64(t.size());
            msg.insert(msg.end(), t.begin(), t.end());
            u64(blk);
            u64(j - i);
            msg.insert(msg.end(), rec + i * R, rec + j * R);
            Digest d{};
            SHA256(msg.data(), msg.size(), d.data());
            leaves.push_back(d);
            i = j;
          }
          // tree + root (v5 rule)
          std::vector<Digest> level = leaves;
          Digest tree{};
          auto node = [&](const Digest &a, const Digest &b) {
            std::vector<unsigned char> msg;
            const std::string t = "native-v5/merkle-node";
            for (int k = 0; k < 8; ++k)
              msg.push_back(static_cast<unsigned char>(std::uint64_t(t.size()) >> (8 * k)));
            msg.insert(msg.end(), t.begin(), t.end());
            msg.insert(msg.end(), a.begin(), a.end());
            msg.insert(msg.end(), b.begin(), b.end());
            Digest d{};
            SHA256(msg.data(), msg.size(), d.data());
            return d;
          };
          if (level.empty()) {
            std::vector<unsigned char> msg;
            const std::string t = "native-v5/merkle-empty";
            for (int k = 0; k < 8; ++k)
              msg.push_back(static_cast<unsigned char>(std::uint64_t(t.size()) >> (8 * k)));
            msg.insert(msg.end(), t.begin(), t.end());
            SHA256(msg.data(), msg.size(), tree.data());
          } else {
            while (level.size() > 1) {
              std::vector<Digest> upl;
              for (std::size_t k = 0; k + 1 < level.size(); k += 2)
                upl.push_back(node(level[k], level[k + 1]));
              if (level.size() % 2)
                upl.push_back(level.back());
              level.swap(upl);
            }
            tree = level[0];
          }
          std::vector<unsigned char> msg;
          auto u64 = [&](std::uint64_t v) {
            for (int b = 0; b < 8; ++b)
              msg.push_back(static_cast<unsigned char>(v >> (8 * b)));
          };
          const std::string t = root_tag;
          u64(t.size());
          msg.insert(msg.end(), t.begin(), t.end());
          u64(1024);
          u64(nrec);
          u64(leaves.size());
          msg.insert(msg.end(), tree.begin(), tree.end());
          Digest d{};
          SHA256(msg.data(), msg.size(), d.data());
          return d;
        };
        std::uint32_t w0[8], w1[8];
        std::memcpy(w0, ir, 32);
        std::memcpy(w1, ir + 32, 32);
        const bool id_ok = digest_of_words(w0) == v5_root("native-v5/carrier-block", "native-v5/carrier-root",
                                                           reinterpret_cast<const unsigned char *>(m.endpoints.data()),
                                                           m.endpoints.size(), 80) &&
                           digest_of_words(w1) == v5_root("native-v5/path-block", "native-v5/path-root",
                                                          reinterpret_cast<const unsigned char *>(m.paths.data()),
                                                          m.paths.size(), 64);
        check(id_ok, "v5 ID-block carrier/path roots on the device = CPU at n=" + std::to_string(n));
        std::printf("RESULT gpu_v5_id_roots n=%zu identical=%d device_ms=%.4f\n", n, int(id_ok), s.elapsed_ms(c0, c1));
        // Dense native IDs 1..n (full 1024-record v5 leaves, as in a native census): carrier rows (10 words) and
        // path rows (8 words) with the same IDs; the two roots separately and in one paired leaf launch.
        {
          std::vector<std::uint64_t> A(10 * n), Bw(8 * n);
          for (std::size_t i = 0; i < n; ++i) {
            for (int w = 0; w < 10; ++w)
              A[10 * i + w] = w ? 0x9e3779b97f4a7c15ull * (i * 16 + w) : i + 1;
            for (int w = 0; w < 8; ++w)
              Bw[8 * i + w] = w ? 0xc2b2ae3d27d4eb4full * (i * 16 + w) : i + 1;
          }
          auto dA = up(A.data(), A.size() * 8);
          auto dB = up(Bw.data(), Bw.size() * 8);
          auto pscr = s.alloc(g7::Forms::id_root_pair_scratch_bytes(n), "v71-pair");
          auto pr = s.alloc(128, "v71-pair-roots");
          s.synchronize("v71 pair alloc");
          double single_ms = 1e300, pair_ms = 1e300;
          for (int rep = 0; rep < 2; ++rep) {
            s.clear_markers();
            const int p0 = s.mark("single");
            forms.id_root(0, n, dA.address, 10, ids.address, pr.address);
            forms.id_root(1, n, dB.address, 8, ids.address, pr.address + 32);
            const int p1 = s.mark("pair");
            forms.id_root_pair(n, dA.address, 10, dB.address, 8, pscr.address, pr.address + 64, pr.address + 96);
            const int p2 = s.mark("end");
            s.synchronize("v71 pair");
            single_ms = std::min(single_ms, s.elapsed_ms(p0, p1));
            pair_ms = std::min(pair_ms, s.elapsed_ms(p1, p2));
          }
          unsigned char q[128];
          down(q, pr, 0, 128);
          std::uint32_t wa[8], wb[8];
          std::memcpy(wa, q, 32);
          std::memcpy(wb, q + 32, 32);
          const bool pair_ok =
              std::memcmp(q, q + 64, 64) == 0 &&
              digest_of_words(wa) == v5_root("native-v5/carrier-block", "native-v5/carrier-root",
                                             reinterpret_cast<const unsigned char *>(A.data()), n, 80) &&
              digest_of_words(wb) == v5_root("native-v5/path-block", "native-v5/path-root",
                                             reinterpret_cast<const unsigned char *>(Bw.data()), n, 64);
          check(pair_ok, "dense-ID v5 roots: separate = paired launch = CPU at n=" + std::to_string(n));
          std::printf("RESULT gpu_v5_id_roots_dense n=%zu identical=%d separate_ms=%.4f paired_ms=%.4f\n", n,
                      int(pair_ok), single_ms, pair_ms);
          for (auto *d : {&dA, &dB, &pscr, &pr})
            s.release(*d);
        }
        for (auto *d : {&dEnd, &dPath, &dSink, &dBirth, &dCen, &roots, &scr, &sep, &sroots, &ids, &idr})
          s.release(*d);
      }
      // Field side on both meshes: |E||J| sum (dt scaled), E^, streaming and block digests.
      for (std::uint64_t mn : {std::uint64_t(45), std::uint64_t(180)}) {
        const HallAxialMesh hm{mn, 6};
        std::vector<double> E(hm.ne()), B(hm.nb()), J(hm.ne(), 0.0);
        std::mt19937_64 g(mn);
        std::normal_distribution<double> N(0, 1);
        for (auto &x : E)
          x = 1e5 * N(g);
        for (auto &x : B)
          x = 1e-2 * N(g);
        for (std::size_t e = 0; e < J.size(); e += 97)
          J[e] = 1e2 * N(g);
        auto dE = up(E.data(), E.size() * 8), dB = up(B.data(), B.size() * 8), dJ = up(J.data(), J.size() * 8);
        auto out = s.alloc(2 * sizeof(g7::SumOut) + 16 + 64 + 3 * 32, "v71-field");
        auto ss = s.alloc(g7::Forms::sum_scratch_bytes(80), "v71-sum");
        auto fb = s.alloc(g7::Forms::field_bound_scratch_bytes(hm.ne(), hm.nn()), "v71-fb");
        std::vector<g7::Forms::BlockListJob> fj = {{pf7::section_field_E, hm.ne(), dE.address, out.address + 2 * sizeof(g7::SumOut) + 16 + 64},
                                                   {pf7::section_field_B, hm.nb(), dB.address, out.address + 2 * sizeof(g7::SumOut) + 16 + 96},
                                                   {pf7::section_current_J, hm.ne(), dJ.address, out.address + 2 * sizeof(g7::SumOut) + 16 + 128}};
        auto fs = s.alloc(g7::Forms::block_lists_scratch_bytes(fj), "v71-fs");
        s.synchronize("v71 alloc");
        s.clear_markers();
        const int d0 = s.mark("field");
        forms.sum_abs_products(hm.ne(), dE.address, 1, dJ.address, 1, 1e-14, 80, ss.address, out.address);
        forms.sum_products(hm.ne(), dE.address, 1, dJ.address, 1, 1e-14, 80, ss.address, out.address + sizeof(g7::SumOut));
        forms.field_bound(hm.ne(), dE.address, dJ.address, mn, 6, fb.address, out.address + 2 * sizeof(g7::SumOut));
        forms.block_lists(fj, fs.address);
        const int d1 = s.mark("field-end");
        forms.stream_digest(0, hm.ne(), dE.address, hm.nb(), dB.address, out.address + 2 * sizeof(g7::SumOut) + 16);
        forms.stream_digest(1, hm.ne(), dJ.address, 0, 0, out.address + 2 * sizeof(g7::SumOut) + 16 + 32);
        const int d2 = s.mark("stream-end");
        s.synchronize("v71 field");
        std::vector<unsigned char> o(2 * sizeof(g7::SumOut) + 16 + 64 + 96);
        down(o.data(), out, 0, o.size());
        g7::SumOut ga, wa;
        std::memcpy(&ga, o.data(), sizeof ga);
        std::memcpy(&wa, o.data() + sizeof ga, sizeof wa);
        double em;
        std::uint32_t bad;
        std::memcpy(&em, o.data() + 2 * sizeof ga, 8);
        std::memcpy(&bad, o.data() + 2 * sizeof ga + 8, 4);
        WorkScales ws;
        work_field_scales(ws, 1e-14, E, J, hm);
        const auto w = field_work_J(1e-14, E, J);
        auto word_hex = [&](std::size_t off) {
          std::uint32_t x[8];
          std::memcpy(x, o.data() + off, 32);
          return hex(digest_of_words(x));
        };
        const std::size_t base = 2 * sizeof(g7::SumOut) + 16;
        const bool ok = bits(ga.scaled) == bits(ws.field_work_gross_J) && bits(wa.scaled) == bits(w.value) &&
                        bits(em) == bits(ws.field_bound_V_m) && !bad &&
                        word_hex(base) == total_field_sha256(E, B) && word_hex(base + 32) == array_sha256(J) &&
                        trial_field_blocks_sha256(word_hex(base + 64), word_hex(base + 96)) == trial_field_blocks_sha256(E, B) &&
                        current_blocks_sha256(word_hex(base + 128)) == current_blocks_sha256(J);
        check(ok, "A1 field side and field/current digests on the device = CPU, mesh n=" + std::to_string(mn));
        std::printf("RESULT gpu_field_a1 mesh_n=%llu edges=%llu identical=%d scales_and_block_digests_ms=%.4f streaming_digests_ms=%.4f\n",
                    (unsigned long long)mn, (unsigned long long)hm.ne(), int(ok), s.elapsed_ms(d0, d1), s.elapsed_ms(d1, d2));
        for (auto *d : {&dE, &dB, &dJ, &out, &ss, &fb, &fs})
          s.release(*d);
      }
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
