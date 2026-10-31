// test_forms_cpu.cpp -- CPU tests of the native v7 producer forms.
// Usage: test_forms_cpu OUT_DIR [max_n]
// Writes OUT_DIR/adversarial.txt (terms and results, verified independently
// by verify_exact.py with exact rational arithmetic) and prints RESULT lines.
#include "datasets.hpp"
#include <chrono>
#include <cstring>
#include <fstream>
#include <numeric>
#include <openssl/sha.h>

using namespace v7test;
namespace {
double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
Digest openssl_sha(const std::vector<unsigned char> &m) {
  Digest d{};
  SHA256(m.data(), m.size(), d.data());
  return d;
}
// Portable core over an arbitrary byte string fed as a mix of bytes and words.
Digest core_sha(const std::vector<unsigned char> &m, std::size_t split) {
  pf7::Sha s;
  pf7::sha_init(s);
  std::size_t i = 0;
  for (; i < split && i < m.size(); ++i)
    pf7::sha_byte(s, m[i]);
  for (; i + 4 <= m.size(); i += 4)
    pf7::sha_word(s, (std::uint32_t(m[i]) << 24) | (std::uint32_t(m[i + 1]) << 16) | (std::uint32_t(m[i + 2]) << 8) |
                         m[i + 3]);
  for (; i < m.size(); ++i)
    pf7::sha_byte(s, m[i]);
  std::uint32_t w[8];
  pf7::sha_finish(s, w);
  return digest_of_words(w);
}
// Block list through the portable core (device arithmetic on the host).
Digest core_block_list(std::uint32_t section, const void *records, std::uint64_t n) {
  const std::uint64_t R = pf7::section_record_bytes(section), chunks = (n + 31) / 32, leaves = (n + 1023) / 1024;
  std::vector<std::uint32_t> cw(8 * chunks), lw(8 * (leaves ? leaves : 1));
  const auto *words = static_cast<const std::uint64_t *>(records);
  for (std::uint64_t j = 0; j < chunks; ++j)
    pf7::chunk_digest(section, j, std::min<std::uint64_t>(32, n - 32 * j), words + j * 32 * (R / 8), &cw[8 * j]);
  for (std::uint64_t k = 0; k < leaves; ++k) {
    const std::uint64_t nk = std::min<std::uint64_t>(1024, n - 1024 * k);
    pf7::block_digest(section, k, nk, (nk + 31) / 32, &cw[8 * 32 * k], &lw[8 * k]);
  }
  std::vector<std::uint32_t> level(lw.begin(), lw.begin() + 8 * leaves);
  std::uint32_t tree[8];
  if (!leaves)
    pf7::merkle_empty(tree);
  else {
    while (level.size() > 8) {
      std::vector<std::uint32_t> up;
      const std::size_t c = level.size() / 8;
      for (std::size_t i = 0; i + 1 < c; i += 2) {
        std::uint32_t o[8];
        pf7::merkle_node(&level[8 * i], &level[8 * i + 8], o);
        up.insert(up.end(), o, o + 8);
      }
      if (c % 2)
        up.insert(up.end(), level.end() - 8, level.end());
      level.swap(up);
    }
    std::copy(level.begin(), level.end(), tree);
  }
  std::uint32_t root[8];
  pf7::section_root(section, n, leaves, tree, root);
  return digest_of_words(root);
}
std::string hex64(double x) {
  char b[24];
  std::snprintf(b, sizeof b, "%016llx", (unsigned long long)bits(x));
  return b;
}
} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: test_forms_cpu OUT_DIR [max_n]\n");
    return 2;
  }
  const std::string out = argv[1];
  const std::size_t max_n = argc > 2 ? std::strtoull(argv[2], nullptr, 10) : 1000000;

  // ---- SHA-256 portable core = OpenSSL (lengths 0..300, every split).
  {
    std::mt19937_64 g(7);
    int bad = 0, n = 0;
    for (std::size_t len = 0; len <= 300; ++len) {
      std::vector<unsigned char> m(len);
      for (auto &c : m)
        c = static_cast<unsigned char>(g());
      for (std::size_t split : {std::size_t(0), std::size_t(1), std::size_t(2), std::size_t(3), len / 2}) {
        bad += core_sha(m, split) != openssl_sha(m);
        ++n;
      }
    }
    check(bad == 0, "portable SHA-256 equals OpenSSL on " + std::to_string(n) + " messages");
    std::printf("RESULT sha256_core_vs_openssl messages=%d mismatches=%d\n", n, bad);
  }

  // ---- (a) adversarial exact sums: written for independent verification.
  {
    std::ofstream f(out + "/adversarial.txt");
    auto sets = adversarial(11);
    for (std::size_t n : {std::size_t(1000), std::size_t(10000), std::size_t(100000)})
      sets.push_back(random_cancel(n, 100 + n));
    for (const auto &s : sets) {
      const auto v = exact_dot(s.a, s.b);
      const double x87 = x87_dot(s.a, s.b);
      f << "set " << s.name << " " << s.a.size() << " " << hex64(v.sum) << " " << v.status << "\n";
      for (std::size_t i = 0; i < s.a.size(); ++i)
        f << hex64(s.a[i]) << " " << hex64(s.b[i]) << "\n";
      std::printf("RESULT adversarial %s terms=%zu exact=%.17g status=%u x87=%.17g ulps_x87=%llu\n", s.name.c_str(),
                  s.a.size(), v.sum, v.status, x87,
                  (unsigned long long)(std::isfinite(v.sum) && std::isfinite(x87) ? ulps(v.sum, x87) : 0));
    }
    // Known values.
    auto value = [&](const char *name) {
      for (const auto &s : sets)
        if (s.name == name)
          return exact_dot(s.a, s.b);
      return SumValue{};
    };
    check(value("subnormal-tie").sum == std::ldexp(1.0, -1073), "3 x 2^-1075 rounds to even 2^-1073");
    check(value("tie-even").sum == 1.0, "1 + 2^-53 ties to even 1");
    check(value("tie-sticky").sum == 1.0 + std::ldexp(1.0, -52), "1 + 2^-53 + 2^-200 rounds up");
    check(bits(value("negative-zero").sum) == 0, "exact zero is +0.0");
    check(value("overflow-cancel").sum == 3.75 && value("overflow-cancel").ok(), "overflowing products cancel exactly");
    check(value("round-overflow").status & pf7::sum_overflow, "rounding above DBL_MAX reports overflow");
    check(value("near-max").ok() && value("near-max").sum == 1.7976931348623157e308, "near-max stays finite");
    const auto wide = value("cancel-pairs-wide");
    check(wide.sum == 3e-300 * 7e-20 || wide.sum == double((long double)3e-300 * 7e-20) ||
              wide.sum == std::fma(3e-300, 7e-20, 0.0),
          "wide cancellation leaves exactly the residue product");
    double nan_a[1] = {std::nan("")}, one[1] = {1};
    check(exact_dot(nan_a, one).status & pf7::sum_nonfinite_input, "NaN factor reported");
  }

  // ---- (a) permutation and partition invariance, 1e3 .. 1e6 terms.
  for (std::size_t n = 1000; n <= max_n; n *= 10) {
    const auto s = random_cancel(n, 900 + n);
    const auto ref = exact_dot(s.a, s.b);
    std::mt19937_64 g(n);
    int bad = 0, cases = 0;
    double t_exact = 0, t_x87 = 0;
    std::uint64_t max_ulps = 0, x87_differs = 0;
    for (int perm = 0; perm < 4; ++perm) {
      std::vector<std::size_t> p(n);
      std::iota(p.begin(), p.end(), 0);
      if (perm)
        std::shuffle(p.begin(), p.end(), g);
      std::vector<double> a(n), b(n);
      for (std::size_t i = 0; i < n; ++i) {
        a[i] = s.a[p[i]];
        b[i] = s.b[p[i]];
      }
      for (unsigned parts : {1u, 2u, 5u, 8u, 64u}) {
        std::vector<ExactSum> part(parts);
        const double t0 = now();
        for (std::size_t i = 0; i < n; ++i)
          part[(i * parts) / n].add_product(a[i], b[i]);
        // Merge through the interchange form in a shuffled order.
        std::vector<unsigned> order(parts);
        std::iota(order.begin(), order.end(), 0u);
        std::shuffle(order.begin(), order.end(), g);
        ExactSum total;
        for (auto k : order)
          total.merge(ExactSum::from_digits(part[k].digits(), part[k].status(), part[k].terms()));
        std::uint32_t st = 0;
        const double v = total.value(&st);
        if (parts == 1)
          t_exact += now() - t0;
        bad += bits(v) != bits(ref.sum) || st != ref.status || total.digits() != ExactSum::from_digits(total.digits()).digits();
        ++cases;
      }
      const double t1 = now();
      const double x = x87_dot(a, b);
      t_x87 += now() - t1;
      x87_differs += bits(x) != bits(ref.sum);
      max_ulps = std::max(max_ulps, ulps(x, ref.sum));
    }
    check(bad == 0, "exact sum identical over permutations and partitions 1/2/5/8/64 at n=" + std::to_string(n));
    std::printf("RESULT partition n=%zu cases=%d mismatches=%d exact=%.17g x87_permutations_differing=%llu/4 "
                "x87_max_ulps=%llu exact_ns_per_term=%.2f x87_ns_per_term=%.2f\n",
                n, cases, bad, ref.sum, (unsigned long long)x87_differs, (unsigned long long)max_ulps,
                1e9 * t_exact / (4.0 * n), 1e9 * t_x87 / (4.0 * n));
  }

  // ---- (a) kinetic change: split = aligned; vs x87 (native v6 form).
  for (std::size_t n = 1000; n <= max_n; n *= 10) {
    const auto m = motion(n, 31 + n);
    const auto split = kinetic_change_J(m.census, m.endpoints, m.sinks);
    const auto aligned = kinetic_change_J_aligned(m.census, m.u_end);
    const double x87 = x87_kinetic(m.census, m.u_end);
    check(split.ok() && bits(split.value) == bits(aligned.value), "kinetic split = aligned at n=" + std::to_string(n));
    std::printf("RESULT kinetic n=%zu sinks=%zu exact=%.17g x87=%.17g ulps=%llu\n", n, m.sinks.size(), split.value,
                x87, (unsigned long long)ulps(split.value, x87));
    auto bad_split = m.endpoints;
    if (!bad_split.empty()) {
      bad_split.pop_back();
      check(kinetic_change_J(m.census, bad_split, m.sinks).status & pf7::sum_census_missing,
            "kinetic: census carrier without end state reported");
    }
  }
  {
    const auto m = motion(4096, 5);
    const double dt = 1e-14;
    std::vector<double> E(41676), J(41676);
    std::mt19937_64 g(3);
    std::normal_distribution<double> N(0, 1);
    for (std::size_t e = 0; e < E.size(); ++e) {
      E[e] = 1e7 * N(g);
      J[e] = 1e3 * N(g);
    }
    const auto w = field_work_J(dt, E, J);
    const double x87 = dt * x87_dot(E, J);
    check(w.ok() && w.value == dt * w.sum, "field work = dt x rounded exact sum");
    std::printf("RESULT field_work edges=%zu exact=%.17g x87=%.17g ulps=%llu\n", E.size(), w.value, x87,
                (unsigned long long)ulps(w.value, x87));
    // One-ulp tamper of a sum is a different value (budget equality refuses).
    check(std::nextafter(w.value, 1e300) != w.value, "one-ulp change of the work sum differs");
  }

  // ---- (b) block lists: OpenSSL reference = portable core; tamper.
  {
    const auto m = motion(5000, 77);
    struct Sec {
      std::uint32_t s;
      const void *p;
      std::size_t n;
    };
    for (std::size_t len : {std::size_t(0), std::size_t(1), std::size_t(31), std::size_t(32), std::size_t(33),
                            std::size_t(1023), std::size_t(1024), std::size_t(1025), std::size_t(4097)}) {
      const std::size_t n = std::min(len, m.endpoints.size());
      for (const auto &x : {Sec{pf7::section_endpoints, m.endpoints.data(), n}, Sec{pf7::section_paths, m.paths.data(), std::min(len, m.paths.size())},
                            Sec{pf7::section_sinks, m.sinks.data(), std::min(len, m.sinks.size())},
                            Sec{pf7::section_births, m.births.data(), std::min(len, m.births.size())},
                            Sec{pf7::section_census, m.census_records.data(), std::min(len, m.census_records.size())}}) {
        const auto ref = block_list(x.s, x.p, x.n);
        check(ref.root_digest == core_block_list(x.s, x.p, x.n),
              "block list: OpenSSL and portable core agree (section " + std::to_string(x.s) + ", n=" + std::to_string(x.n) + ")");
        check(block_root(x.s, x.n, ref.leaves) == ref.root, "block list: root from leaves");
      }
    }
    auto ep = m.endpoints;
    const auto r0 = block_list(pf7::section_endpoints, ep.data(), ep.size());
    ep[2047].position_m[1] = std::nextafter(ep[2047].position_m[1], 1.0);
    const auto r1 = block_list(pf7::section_endpoints, ep.data(), ep.size());
    int changed = 0;
    for (std::size_t k = 0; k < r0.leaves.size(); ++k)
      changed += r0.leaves[k] != r1.leaves[k];
    check(r0.root != r1.root && changed == 1, "one-ulp change of one endpoint changes exactly one leaf and the root");
    auto leaves = r0.leaves;
    leaves[0][5] ^= 1;
    check(block_root(pf7::section_endpoints, ep.size(), leaves) != r0.root, "a tampered leaf changes the root");
    // Section tags separate equal bytes.
    check(block_list(pf7::section_paths, m.paths.data(), 10).root != block_list(pf7::section_births, m.paths.data(), 10).root,
          "section tags separate sections");
    // Record digest: one-ulp change of K, one section root, one header field.
    MotionHeaderV7 h;
    h.accepted_step = 3;
    h.dt_s = 1e-14;
    h.budget = {1e-15, 1e-15, 0, 0, .01, .01, .01};
    h.census_sha256 = std::string(64, 'a');
    const auto sec = motion_sections_v7(m.endpoints, m.paths, m.sinks, m.births);
    const auto d0 = motion_record_v7_sha256(h, sec);
    auto h1 = h;
    h1.budget[0] = std::nextafter(h.budget[0], 1.0);
    auto sec1 = sec;
    sec1.paths_root[0] = sec1.paths_root[0] == '0' ? '1' : '0';
    auto h2 = h;
    h2.neutral_stream_dt_max_s = 2e-14;
    check(motion_record_v7_sha256(h1, sec) != d0 && motion_record_v7_sha256(h, sec1) != d0 &&
              motion_record_v7_sha256(h2, sec) != d0,
          "record digest binds K to the ulp, every section root and the neutral bound");
    std::printf("RESULT motion_record_v7_example %s\n", d0.c_str());
  }

  // ---- (c) census: canonical check, host records, ordering, tamper.
  {
    const auto m = motion(20000, 9);
    const auto chk = check_census(m.census_records, m.states.size());
    check(chk.ok(), "generated census passes the canonical check");
    const auto digest = census_v2_sha256(m.census_records, m.states);
    // Host entries (packet ID text, state text) give the same records.
    std::vector<std::string> pid(m.census.size());
    std::vector<HostCensusEntry> entries;
    for (std::size_t i = 0; i < m.census.size(); ++i)
      pid[i] = "packet/" + std::to_string(m.census_records[i].packet_key);
    for (std::size_t i = m.census.size(); i-- > 0;) // reverse order input
      entries.push_back({pid[i], m.states[m.census_records[i].state_index], &m.census[i]});
    const auto rebuilt = census_records(entries, m.states);
    check(rebuilt.size() == m.census_records.size() &&
              std::memcmp(rebuilt.data(), m.census_records.data(), rebuilt.size() * sizeof(CensusRecord)) == 0,
          "census records from host packet/state texts in any order");
    std::vector<std::uint64_t> keys, states;
    for (const auto &r : m.census_records) {
      keys.push_back(r.packet_key);
      states.push_back(r.state_index);
    }
    check(census_v2_sha256(m.census.size(), census_states_sha256(m.states),
                           census_block_list(m.census, keys, states).root) == digest,
          "census digest from rows + keys equals the record form");
    auto expect = [&](std::vector<CensusRecord> c, std::uint32_t code, const char *what) {
      const auto k = check_census(c, m.states.size());
      check(k.code == code, std::string("census tamper refused: ") + what + " (code " + std::to_string(k.code) + ")");
    };
    auto c = m.census_records;
    std::swap(c[100], c[101]);
    expect(c, pf7::census_not_ascending, "reordered entry");
    c = m.census_records;
    c.insert(c.begin() + 51, c[50]);
    expect(c, pf7::census_not_ascending, "duplicated entry");
    c = m.census_records;
    c[1].carrier.position_m[0] = std::nextafter(c[1].carrier.position_m[0], 1.0); // 2nd carrier of packet 0
    expect(c, pf7::census_divergent, "divergent co-located carrier");
    c = m.census_records;
    c[1].packet_key = c[1].carrier.id + 1;
    expect(c, pf7::census_key_invalid, "packet key above the carrier ID");
    c = m.census_records;
    c[1].packet_key = c[0].carrier.id + 1 < c[1].carrier.id ? c[0].carrier.id + 1 : c[0].carrier.id;
    if (c[1].packet_key != c[0].carrier.id)
      expect(c, pf7::census_key_absent, "packet key of no census carrier");
    c = m.census_records;
    c[8].packet_key = c[7].carrier.id; // 7 is the second carrier of a packet
    if (c[7].packet_key != c[7].carrier.id)
      expect(c, pf7::census_key_not_first, "packet key of a non-first carrier");
    c = m.census_records;
    c[3].state_index = 3;
    expect(c, pf7::census_state_range, "state index out of range");
    c = m.census_records;
    c[0].carrier.id = 0;
    expect(c, pf7::census_id_zero, "native ID zero");
    // Ordering: shuffled census -> canonical order -> same digest.
    auto shuffled = m.census_records;
    std::mt19937_64 g(1);
    std::shuffle(shuffled.begin(), shuffled.end(), g);
    const auto p = census_order(shuffled);
    std::vector<CensusRecord> ordered;
    for (auto i : p)
      ordered.push_back(shuffled[i]);
    check(census_v2_sha256(ordered, m.states) == digest, "canonical census order restores the digest");
    std::printf("RESULT census carriers=%zu digest=%s\n", m.census.size(), digest.c_str());
  }

  // ---- neutral stream: bound parameter.
  {
    KineticPolicy k{0.005, 0.005, -0.0125, 0.0125, {0.0, 0.005, 0.0023, 0.0025}, {0.0, 0.005, 0.0027, 0.0029}};
    k.lower[0] = 0.0001;
    k.upper[0] = 0.0001;
    k.lower[1] = 0.0049;
    k.upper[1] = 0.0049;
    Carrier c{7, {0.002, 0.0026, 0.0}, {300, -200, 100}, 0, 2.99e-26, 1e9};
    bool refused = false;
    try {
      (void)neutral_stream(k, c, 0, 2e-14, 2e-14);
    } catch (const std::invalid_argument &e) {
      refused = std::string(e.what()) == "neutral authoritative interval";
    }
    check(refused, "default bound refuses dt = 2e-14 with the sealed text");
    const auto a = neutral_stream(k, c, 0, 1e-14, 1e-14);
    const auto b = neutral_stream(k, c, 0, 2e-14, 2e-14, 1e-13);
    check(a.hit.sink == -1 && b.hit.sink == -1 && b.displacement_m[0] == 2e-14 * 300 / std::sqrt(1 + (300.0 * 300 + 200.0 * 200 + 100.0 * 100) / pf7::c_squared),
          "configured bound admits 2e-14 with the sealed arithmetic");
  }

  std::printf("RESULT cpu_checks=%d failures=%d\n", checks, failures);
  return failures ? 1 : 0;
}
