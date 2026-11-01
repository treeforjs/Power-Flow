// Native v8 persistent containers (API-008.txt 2.1): randomized equivalence with
// std::map / std::set, structural sharing, content-defined chunk canonicity and
// Merkle digest equality for equal content built by different histories.
#include "raikou/persistent_containers.h"
#include <openssl/sha.h>
#include <cstdio>
#include <random>
#include <string>

using namespace raikou::mfem_backend;
static int failures = 0, checks = 0;
static void check(bool ok, const char *what) {
  ++checks;
  if (!ok) {
    ++failures;
    std::printf("FAIL %s\n", what);
  }
}
using Leaf = std::array<unsigned char, 32>;
template <class E, class F> static Leaf leaf_of(const std::vector<E> &v, F enc) {
  std::string bytes = "leaf";
  for (const auto &e : v)
    bytes += enc(e) + '\n';
  Leaf out;
  SHA256(reinterpret_cast<const unsigned char *>(bytes.data()), bytes.size(), out.data());
  return out;
}
static std::string root_of(std::uint64_t n, const std::vector<Leaf> &leaves) {
  std::string bytes = "root" + std::to_string(n);
  for (const auto &l : leaves)
    bytes.append(reinterpret_cast<const char *>(l.data()), 32);
  unsigned char d[32];
  SHA256(reinterpret_cast<const unsigned char *>(bytes.data()), bytes.size(), d);
  char hx[65];
  for (int i = 0; i < 32; ++i)
    std::snprintf(hx + 2 * i, 3, "%02x", d[i]);
  return std::string(hx, 64);
}
template <class M> static std::string map_digest(const M &m) {
  return m.merkle_digest(
      [](const std::pair<std::string, int> &e) {
        return leaf_of(std::vector<std::string>{e.first + "=" + std::to_string(e.second)},
                       [](const std::string &x) { return x; });
      },
      [](const std::vector<Leaf> &v) {
        std::string b = "leaf";
        for (const auto &l : v)
          b.append(reinterpret_cast<const char *>(l.data()), 32);
        Leaf out;
        SHA256(reinterpret_cast<const unsigned char *>(b.data()), b.size(), out.data());
        return out;
      },
      root_of);
}
int main() {
  std::mt19937_64 rng(20261031);
  // 1. randomized operations against std::map, with snapshots sharing chunks
  for (int round = 0; round < 6; ++round) {
    SharedMap<std::string, int> m;
    std::map<std::string, int> ref;
    std::vector<std::pair<SharedMap<std::string, int>, std::map<std::string, int>>> snaps;
    const int ops = round < 3 ? 4000 : 40000;
    const int space = round % 2 ? 500 : 20000;
    for (int k = 0; k < ops; ++k) {
      const std::string key = "packet/" + std::to_string(rng() % space);
      const int op = int(rng() % 10);
      if (op < 5) {
        const int v = int(rng() % 1000);
        auto r = m.emplace(key, v);
        auto q = ref.emplace(key, v);
        if (r.second != q.second || r.first->second != q.first->second)
          check(false, "emplace result");
      } else if (op < 8) {
        if (m.erase(key) != ref.erase(key))
          check(false, "erase count");
      } else if (op == 8) {
        if (ref.count(key)) {
          const int v = int(rng() % 1000);
          m.modify(key) = v;
          ref[key] = v;
        }
      } else if (snaps.size() < 16) {
        snaps.emplace_back(m, ref);
      }
      if (m.count(key) != ref.count(key))
        check(false, "count");
    }
    check(m.size() == ref.size(), "size");
    check(std::equal(m.begin(), m.end(), ref.begin(), ref.end(),
                     [](const auto &a, const auto &b) { return a.first == b.first && a.second == b.second; }),
          "iteration order and values");
    for (auto &[s, r] : snaps) {
      check(s.size() == r.size() && std::equal(s.begin(), s.end(), r.begin(), r.end(),
                                               [](const auto &a, const auto &b) {
                                                 return a.first == b.first && a.second == b.second;
                                               }),
            "snapshot unaffected by later mutation");
    }
    // canonical chunking: the same content built in bulk gives the same digest
    SharedMap<std::string, int> bulk(ref);
    check(bulk.chunk_count() == m.chunk_count(), "chunk count canonical");
    check(map_digest(bulk) == map_digest(m), "Merkle digest canonical");
    // one changed value changes the digest; restoring it restores the digest
    if (!ref.empty()) {
      const auto d0 = map_digest(m);
      auto copy = m;
      const auto &k0 = ref.begin()->first;
      copy.modify(k0) += 1;
      check(map_digest(copy) != d0, "digest sensitive");
      copy.modify(k0) -= 1;
      check(map_digest(copy) == d0, "digest restored");
      check(map_digest(m) == d0, "original digest unchanged");
    }
  }
  // 2. sets
  {
    SharedSet<std::string> s;
    std::set<std::string> ref;
    for (int k = 0; k < 30000; ++k) {
      const std::string key = "arrival-v1/s" + std::to_string(rng() % 100) + "/p0/k" + std::to_string(rng() % 400);
      auto a = s.insert(key);
      auto b = ref.insert(key);
      if (a.second != b.second)
        check(false, "set insert result");
    }
    check(s.get() == ref, "set content");
    SharedSet<std::string> bulk(ref);
    check(bulk == s, "set equality");
    check(bulk.chunk_count() == s.chunk_count(), "set chunk count canonical");
  }
  // 3. structural sharing: a copy shares; one mutation of a large map touches few chunks
  {
    std::map<std::string, int> ref;
    for (int k = 0; k < 200000; ++k)
      ref.emplace("p" + std::to_string(k), k);
    SharedMap<std::string, int> m(ref);
    const auto d0 = map_digest(m);
    auto c = m;
    check(c.shares(m), "copy shares");
    c.emplace("p-new", 1);
    check(!c.shares(m) && m.size() == 200000 && c.size() == 200001, "copy on write");
    check(map_digest(m) == d0, "parent digest cached and unchanged");
    std::map<std::string, int> ref2 = ref;
    ref2.emplace("p-new", 1);
    check(map_digest(c) == map_digest(SharedMap<std::string, int>(ref2)), "incremental digest equals bulk digest");
    std::printf("chunks %zu for %zu entries\n", m.chunk_count(), m.size());
  }
  std::printf("%d container checks, %d failures\n", checks, failures);
  return failures ? 1 : 0;
}
