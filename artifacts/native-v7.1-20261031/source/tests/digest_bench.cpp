// Per-step native row hashing cost versus carrier count (1e4, 1e5, 1e6):
// v5 block Merkle digest (full, one changed carrier, every carrier moved),
// host recomputation of the root from supplied leaves, sampled leaf
// verification, and the candidate003-style JSON row dump + SHA-256 baseline.
#include "raikou/mfem_pic_canonical.h"
#include "raikou/mfem_pic_arena.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <openssl/sha.h>
namespace m = raikou::mfem_backend;
namespace cn = m::canonical;
using J = m::ArenaJson;
template <class F> double median_s(int reps, F f) {
  std::vector<double> t;
  for (int r = 0; r < reps; ++r) {
    const auto t0 = std::chrono::steady_clock::now();
    f();
    t.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
  }
  std::sort(t.begin(), t.end());
  return t[t.size() / 2];
}
int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  J out = J::array();
  volatile std::size_t sink = 0;
  for (std::size_t n : {std::size_t(10000), std::size_t(100000), std::size_t(1000000)}) {
    std::vector<m::GapParticle> rows(n);
    for (std::size_t i = 0; i < n; ++i) {
      rows[i].id = i + 1;
      rows[i].position_m = {.002 + 1e-9 * double(i), .004, -.01};
      rows[i].proper_velocity_m_s = {0, 0, 1e6};
      rows[i].weight = 4;
    }
    const int reps = n >= 1000000 ? 5 : 15;
    std::shared_ptr<const cn::BlockDigest> base;
    J row = {{"carriers", n}, {"block", cn::carrier_block}};
    row["full_digest_s"] = median_s(reps, [&] { base = cn::carrier_rows_digest(rows); });
    auto one = rows;
    one[n / 2].position_m[2] += 1e-12;
    std::uint64_t rehashed = 0;
    row["one_changed_carrier_s"] = median_s(reps, [&] {
      sink += cn::carrier_rows_digest(one, &rows, base.get(), &rehashed)->root.size();
    });
    row["one_changed_rehashed_blocks"] = rehashed;
    auto all = rows;
    for (auto &p : all)
      p.position_m[2] += 5e-9;
    row["all_moved_with_parent_s"] = median_s(reps, [&] {
      sink += cn::carrier_rows_digest(all, &rows, base.get(), &rehashed)->root.size();
    });
    row["all_moved_rehashed_blocks"] = rehashed;
    row["root_from_supplied_leaves_s"] = median_s(reps, [&] {
      sink += cn::merkle_root("native-v5/carrier-root", base->leaves, n).size();
    });
    std::vector<std::size_t> sample;
    for (std::size_t k = 0; k < 16 && k < base->leaves.size(); ++k)
      sample.push_back((k * 2654435761u) % base->leaves.size());
    bool ok = true;
    row["host_verify_root_plus_16_leaves_s"] = median_s(reps, [&] {
      ok = ok && cn::verify_carrier_leaves(rows, base->leaves, base->root, sample);
    });
    row["host_verification_passed"] = ok;
    row["candidate003_json_rows_dump_sha_s"] = median_s(n >= 1000000 ? 3 : reps, [&] {
      J j = J::array();
      for (const auto &p : rows)
        j.push_back({{"id", p.id}, {"position_m", p.position_m},
                     {"u", p.proper_velocity_m_s}, {"q", p.charge_C},
                     {"m", p.mass_kg}, {"weight", p.weight}});
      const auto text = j.dump();
      unsigned char d[32];
      SHA256(reinterpret_cast<const unsigned char *>(text.data()), text.size(), d);
      sink += d[0];
    });
    row["hall_census_binary_s"] = median_s(reps, [&] {
      cn::Hasher h;
      h.w.tag("canonical-motion-census-v1");
      h.w.u64(n);
      for (const auto &p : rows) {
        h.w.text("packet");
        h.w.text("e");
        h.w.u64(p.id);
        h.w.vec(p.position_m);
        h.w.vec(p.proper_velocity_m_s);
        h.w.f64(p.charge_C);
        h.w.f64(p.mass_kg);
        h.w.f64(p.weight);
        h.flush_if_large();
      }
      sink += h.finish()[0];
    });
    std::cerr << row.dump() << '\n';
    out.push_back(row);
  }
  std::ofstream(argv[1]) << J{{"digest_scaling", out}}.dump(2) << '\n';
  return 0;
}
