// persistent_containers.h -- native v8 persistent sorted map/set with
// content-defined chunks (API-008.txt 2.1).
//
// A container holds its entries in ascending key order in immutable chunks.
// Chunk boundaries are content defined: an entry whose key is a boundary key
// (cdc_boundary) ends its chunk, every other entry does not; only the last
// chunk may end with a non-boundary key. The chunking is therefore a function
// of the content alone (two containers with equal entries have equal chunks,
// whatever their insertion history).
//
// Entries are held through shared pointers: copying a container copies one
// pointer; the first mutation of a shared container copies its chunk-pointer
// table (O(N / 64) pointers); a mutation copies the entry pointers of the
// chunks it touches (expected 64) and never an entry it does not change.
// Iteration order and lookups equal those of std::map / std::set with
// std::less<K>, so every v5 digest that streams the entries is unchanged.
//
// v8 Merkle digest: entry digest (caller-defined, cached per entry while the
// entry is unchanged), leaf = digest of one chunk's entry digests (cached per
// chunk while the chunk is unchanged), root = digest over the leaves. A
// mutation recomputes only the entries and leaves it touched.
#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace raikou::mfem_backend {
namespace persistent_detail {
// Deterministic 64-bit key hash (FNV-1a with the splitmix64 finaliser).
inline std::uint64_t key_hash(const std::string &k) noexcept {
  std::uint64_t h = 0xcbf29ce484222325ull;
  for (unsigned char c : k) {
    h ^= c;
    h *= 0x100000001b3ull;
  }
  h ^= h >> 30;
  h *= 0xbf58476d1ce4e5b9ull;
  h ^= h >> 27;
  h *= 0x94d049bb133111ebull;
  h ^= h >> 31;
  return h;
}
inline std::uint64_t key_hash(std::uint64_t k) noexcept {
  std::uint64_t h = k + 0x9e3779b97f4a7c15ull;
  h ^= h >> 30;
  h *= 0xbf58476d1ce4e5b9ull;
  h ^= h >> 27;
  h *= 0x94d049bb133111ebull;
  h ^= h >> 31;
  return h;
}
// Expected chunk size 64.
template <class K> inline bool cdc_boundary(const K &k) noexcept {
  return (key_hash(k) & 63u) == 0u;
}
using Digest32 = std::array<unsigned char, 32>;
template <class E> struct Slot {
  std::shared_ptr<E> e;
  mutable Digest32 digest{};
  mutable bool valid = false;
};
template <class E> struct Chunk {
  std::vector<Slot<E>> e;
  mutable Digest32 leaf{};
  mutable bool leaf_valid = false;
};
template <class E> using ChunkPtr = std::shared_ptr<Chunk<E>>;
template <class E> struct Table {
  std::vector<ChunkPtr<E>> chunks;
  std::size_t n = 0;
  mutable std::string digest; // v5 streaming digest cache
  mutable std::string merkle; // v8 Merkle digest cache
};
} // namespace persistent_detail

// E: entry type; KeyOf: functor returning the key of an entry.
template <class K, class E, class KeyOf> class PersistentSorted {
protected:
  using Slot = persistent_detail::Slot<E>;
  using Chunk = persistent_detail::Chunk<E>;
  using ChunkPtr = persistent_detail::ChunkPtr<E>;
  using Table = persistent_detail::Table<E>;
  std::shared_ptr<Table> t_;
  static const std::shared_ptr<Table> &empty_table() {
    static const auto e = std::make_shared<Table>();
    return e;
  }
  static const K &key(const E &e) { return KeyOf{}(e); }
  static const K &key(const Slot &s) { return KeyOf{}(*s.e); }

public:
  using size_type = std::size_t;
  using Digest32 = persistent_detail::Digest32;
  class const_iterator {
  public:
    using value_type = E;
    using reference = const E &;
    using pointer = const E *;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;
    const_iterator() = default;
    const_iterator(const Table *t, std::size_t c, std::size_t i) : t_(t), c_(c), i_(i) {}
    reference operator*() const { return *t_->chunks[c_]->e[i_].e; }
    pointer operator->() const { return t_->chunks[c_]->e[i_].e.get(); }
    const_iterator &operator++() {
      if (++i_ >= t_->chunks[c_]->e.size()) {
        ++c_;
        i_ = 0;
      }
      return *this;
    }
    const_iterator operator++(int) {
      auto r = *this;
      ++*this;
      return r;
    }
    bool operator==(const const_iterator &o) const noexcept { return c_ == o.c_ && i_ == o.i_; }
    bool operator!=(const const_iterator &o) const noexcept { return !(*this == o); }

  private:
    const Table *t_ = nullptr;
    std::size_t c_ = 0, i_ = 0;
  };
  using iterator = const_iterator;
  PersistentSorted() : t_(empty_table()) {}
  // Bulk construction from entries in strictly ascending key order.
  template <class It> void assign_sorted(It first, It last) {
    auto t = std::make_shared<Table>();
    ChunkPtr cur;
    bool has_prev = false;
    K prev{};
    for (; first != last; ++first) {
      const E &e = *first;
      if (has_prev && !(prev < key(e)))
        throw std::invalid_argument("persistent container: keys not strictly ascending");
      prev = key(e);
      has_prev = true;
      if (!cur)
        cur = std::make_shared<Chunk>();
      cur->e.push_back(Slot{std::make_shared<E>(e), {}, false});
      ++t->n;
      if (persistent_detail::cdc_boundary(key(e))) {
        t->chunks.push_back(std::move(cur));
        cur.reset();
      }
    }
    if (cur)
      t->chunks.push_back(std::move(cur));
    t_ = std::move(t);
  }
  size_type size() const noexcept { return t_->n; }
  bool empty() const noexcept { return t_->n == 0; }
  const_iterator begin() const noexcept { return const_iterator(t_.get(), 0, 0); }
  const_iterator end() const noexcept { return const_iterator(t_.get(), t_->chunks.size(), 0); }
  const_iterator cbegin() const noexcept { return begin(); }
  const_iterator cend() const noexcept { return end(); }
  const_iterator find(const K &k) const {
    std::size_t c, i;
    if (!locate(k, c, i))
      return end();
    return const_iterator(t_.get(), c, i);
  }
  size_type count(const K &k) const {
    std::size_t c, i;
    return locate(k, c, i) ? 1 : 0;
  }
  bool shares(const PersistentSorted &o) const noexcept { return t_ == o.t_; }
  std::shared_ptr<const void> identity() const noexcept { return t_; }
  const std::string &cached_digest() const noexcept { return t_->digest; }
  void cache_digest(std::string d) const { t_->digest = std::move(d); }
  std::size_t chunk_count() const noexcept { return t_->chunks.size(); }
  // v8 Merkle digest (API-008.txt 2.2): entry(E) -> 32-byte digest (cached per
  // entry); leaf(entry digests of one chunk) -> 32-byte digest (cached per
  // chunk); root(entries, leaf digests) -> hex digest.
  const std::string &merkle_digest(const std::function<Digest32(const E &)> &entry,
                                   const std::function<Digest32(const std::vector<Digest32> &)> &leaf,
                                   const std::function<std::string(std::uint64_t, const std::vector<Digest32> &)> &root)
      const {
    if (!t_->merkle.empty())
      return t_->merkle;
    std::vector<Digest32> leaves, entries;
    leaves.reserve(t_->chunks.size());
    for (const auto &c : t_->chunks) {
      if (!c->leaf_valid) {
        entries.clear();
        for (const auto &s : c->e) {
          if (!s.valid) {
            s.digest = entry(*s.e);
            s.valid = true;
          }
          entries.push_back(s.digest);
        }
        c->leaf = leaf(entries);
        c->leaf_valid = true;
      }
      leaves.push_back(c->leaf);
    }
    t_->merkle = root(t_->n, leaves);
    return t_->merkle;
  }

protected:
  // Finds key k: chunk c, index i. Returns false when absent (c, i then give
  // the insertion position: chunk c (possibly the last chunk when k is larger
  // than every key) and index i inside it).
  bool locate(const K &k, std::size_t &c, std::size_t &i) const {
    const auto &ch = t_->chunks;
    if (ch.empty()) {
      c = 0;
      i = 0;
      return false;
    }
    std::size_t lo = 0, hi = ch.size(); // first chunk whose last key >= k
    while (lo < hi) {
      const std::size_t mid = lo + (hi - lo) / 2;
      if (key(ch[mid]->e.back()) < k)
        lo = mid + 1;
      else
        hi = mid;
    }
    if (lo == ch.size()) {
      c = ch.size() - 1;
      i = ch[c]->e.size();
      return false;
    }
    c = lo;
    const auto &v = ch[c]->e;
    auto it = std::lower_bound(v.begin(), v.end(), k, [](const Slot &s, const K &x) { return key(s) < x; });
    i = std::size_t(it - v.begin());
    return it != v.end() && !(k < key(*it));
  }
  // Unique table (copy of the chunk-pointer table when shared); clears the
  // table digests.
  Table &own() {
    if (t_.use_count() != 1 || t_ == empty_table()) {
      auto t = std::make_shared<Table>();
      t->chunks = t_->chunks;
      t->n = t_->n;
      t_ = std::move(t);
    }
    t_->digest.clear();
    t_->merkle.clear();
    return *t_;
  }
  // Unique, mutable chunk c of an owned table (copies the entry pointers and
  // their cached digests when shared; invalidates the leaf).
  Chunk &own_chunk(Table &t, std::size_t c) {
    if (t.chunks[c].use_count() != 1)
      t.chunks[c] = std::make_shared<Chunk>(Chunk{t.chunks[c]->e, {}, false});
    t.chunks[c]->leaf_valid = false;
    return *t.chunks[c];
  }
  // Inserts e (key absent). Returns the position of the new entry.
  const_iterator insert_new(E e) {
    std::size_t c, i;
    if (locate(key(e), c, i))
      return const_iterator(t_.get(), c, i);
    auto &t = own();
    const bool boundary = persistent_detail::cdc_boundary(key(e));
    Slot slot{std::make_shared<E>(std::move(e)), {}, false};
    if (t.chunks.empty()) {
      auto ch = std::make_shared<Chunk>();
      ch->e.push_back(std::move(slot));
      t.chunks.push_back(std::move(ch));
      ++t.n;
      return const_iterator(t_.get(), 0, 0);
    }
    if (i == t.chunks[c]->e.size() && persistent_detail::cdc_boundary(key(t.chunks[c]->e.back()))) {
      // after a last chunk that ends with a boundary key: a new last chunk
      auto ch = std::make_shared<Chunk>();
      ch->e.push_back(std::move(slot));
      t.chunks.push_back(std::move(ch));
      ++t.n;
      return const_iterator(t_.get(), t.chunks.size() - 1, 0);
    }
    auto &ch = own_chunk(t, c);
    ch.e.insert(ch.e.begin() + std::ptrdiff_t(i), std::move(slot));
    ++t.n;
    if (boundary && i + 1 < ch.e.size()) {
      // split after the new boundary entry
      auto tail = std::make_shared<Chunk>();
      tail->e.assign(std::make_move_iterator(ch.e.begin() + std::ptrdiff_t(i + 1)),
                     std::make_move_iterator(ch.e.end()));
      ch.e.erase(ch.e.begin() + std::ptrdiff_t(i + 1), ch.e.end());
      t.chunks.insert(t.chunks.begin() + std::ptrdiff_t(c + 1), std::move(tail));
    }
    return const_iterator(t_.get(), c, i);
  }
  size_type erase_key(const K &k) {
    std::size_t c, i;
    if (!locate(k, c, i))
      return 0;
    auto &t = own();
    auto &ch = own_chunk(t, c);
    const bool was_last = i + 1 == ch.e.size();
    ch.e.erase(ch.e.begin() + std::ptrdiff_t(i));
    --t.n;
    if (was_last && c + 1 < t.chunks.size()) {
      // the erased entry ended this chunk (a boundary key): merge with the next
      const auto &next = t.chunks[c + 1]->e;
      ch.e.insert(ch.e.end(), next.begin(), next.end());
      t.chunks.erase(t.chunks.begin() + std::ptrdiff_t(c + 1));
    }
    if (ch.e.empty())
      t.chunks.erase(t.chunks.begin() + std::ptrdiff_t(c));
    return 1;
  }
  // Mutable entry: a private copy of the entry (the chunk is made unique, the
  // entry digest and the table digests are cleared).
  E &mutable_entry(const K &k) {
    std::size_t c, i;
    if (!locate(k, c, i))
      throw std::out_of_range("persistent container: key absent");
    auto &t = own();
    auto &s = own_chunk(t, c).e[i];
    s.e = std::make_shared<E>(*s.e);
    s.valid = false;
    return *s.e;
  }
};

namespace persistent_detail {
template <class K, class V> struct MapKey {
  const K &operator()(const std::pair<K, V> &e) const noexcept { return e.first; }
};
template <class K> struct SetKey {
  const K &operator()(const K &e) const noexcept { return e; }
};
} // namespace persistent_detail

// Persistent ordered map (interface of the v7.1 SharedMap).
template <class K, class V>
class SharedMap : public PersistentSorted<K, std::pair<K, V>, persistent_detail::MapKey<K, V>> {
  using Base = PersistentSorted<K, std::pair<K, V>, persistent_detail::MapKey<K, V>>;

public:
  using map_type = std::map<K, V>;
  using const_iterator = typename Base::const_iterator;
  using iterator = const_iterator;
  using size_type = std::size_t;
  SharedMap() = default;
  SharedMap(const map_type &m) { this->assign_sorted(m.begin(), m.end()); }
  const V &at(const K &k) const {
    auto it = this->find(k);
    if (it == this->end())
      throw std::out_of_range("SharedMap::at");
    return it->second;
  }
  // Materialised copy (O(N); serialisation only).
  map_type get() const { return map_type(this->begin(), this->end()); }
  template <class Key, class... A> std::pair<const_iterator, bool> emplace(Key &&k, A &&...a) {
    K key(std::forward<Key>(k));
    if (auto it = this->find(key); it != this->end())
      return {it, false};
    auto it = this->insert_new(std::pair<K, V>(std::move(key), V(std::forward<A>(a)...)));
    return {it, true};
  }
  size_type erase(const K &k) { return this->erase_key(k); }
  // v8: in-place update of an existing entry's value (private copy of the entry).
  V &modify(const K &k) { return this->mutable_entry(k).second; }
};

// Persistent ordered set (interface of the v7.1 SharedSet).
template <class K> class SharedSet : public PersistentSorted<K, K, persistent_detail::SetKey<K>> {
  using Base = PersistentSorted<K, K, persistent_detail::SetKey<K>>;

public:
  using set_type = std::set<K>;
  using const_iterator = typename Base::const_iterator;
  using iterator = const_iterator;
  using size_type = std::size_t;
  SharedSet() = default;
  SharedSet(const set_type &s) { this->assign_sorted(s.begin(), s.end()); }
  set_type get() const { return set_type(this->begin(), this->end()); }
  std::pair<const_iterator, bool> insert(const K &k) {
    if (auto it = this->find(k); it != this->end())
      return {it, false};
    return {this->insert_new(k), true};
  }
  bool operator==(const SharedSet &o) const {
    if (this->shares(o))
      return true;
    if (this->size() != o.size())
      return false;
    return std::equal(this->begin(), this->end(), o.begin());
  }
};
} // namespace raikou::mfem_backend
