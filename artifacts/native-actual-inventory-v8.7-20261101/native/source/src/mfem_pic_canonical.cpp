#include "raikou/mfem_pic_canonical.h"
#include <algorithm>
#include <openssl/evp.h>
#include <openssl/sha.h>

namespace raikou::mfem_backend::canonical {
namespace {
void need(bool ok, const char *message) {
  if (!ok)
    throw EncodingRefusal(message);
}
static_assert(sizeof(GapParticle) == 80 &&
              std::is_trivially_copyable_v<GapParticle>);
// Canonical carrier record is exactly the little-endian field sequence
// [u64 id, 9 x f64]; on this little-endian target GapParticle has that byte
// layout with no padding, so the record is the object representation.
static_assert(std::endian::native == std::endian::little);
static_assert(offsetof(GapParticle, id) == 0 &&
              offsetof(GapParticle, position_m) == 8 &&
              offsetof(GapParticle, proper_velocity_m_s) == 32 &&
              offsetof(GapParticle, charge_C) == 56 &&
              offsetof(GapParticle, mass_kg) == 64 &&
              offsetof(GapParticle, weight) == 72);
Digest node(const Digest &a, const Digest &b) {
  Writer w;
  w.tag("native-v5/merkle-node");
  w.raw(a.data(), a.size());
  w.raw(b.data(), b.size());
  return sha256_raw(w.bytes);
}
} // namespace
Digest sha256_raw(std::string_view bytes) {
  Digest d{};
  SHA256(reinterpret_cast<const unsigned char *>(bytes.data()), bytes.size(),
         d.data());
  return d;
}
std::string hex(const Digest &d) {
  static constexpr char digits[] = "0123456789abcdef";
  std::string out(64, '0');
  for (std::size_t i = 0; i < 32; ++i) {
    out[2 * i] = digits[d[i] >> 4];
    out[2 * i + 1] = digits[d[i] & 15];
  }
  return out;
}
std::string sha256_hex(std::string_view bytes) { return hex(sha256_raw(bytes)); }
bool is_digest(const std::string &s) noexcept {
  return s.size() == 64 &&
         std::all_of(s.begin(), s.end(), [](char c) {
           return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
         });
}
Digest from_hex(const std::string &s) {
  need(is_digest(s), "lower-case 64-digit hex digest required");
  Digest d{};
  auto v = [](char c) { return c <= '9' ? c - '0' : c - 'a' + 10; };
  for (std::size_t i = 0; i < 32; ++i)
    d[i] = static_cast<unsigned char>(16 * v(s[2 * i]) + v(s[2 * i + 1]));
  return d;
}
Hasher::Hasher() : ctx_(EVP_MD_CTX_new()) {
  if (!ctx_ || EVP_DigestInit_ex(static_cast<EVP_MD_CTX *>(ctx_), EVP_sha256(),
                                 nullptr) != 1) {
    EVP_MD_CTX_free(static_cast<EVP_MD_CTX *>(ctx_));
    throw std::bad_alloc();
  }
}
Hasher::~Hasher() { EVP_MD_CTX_free(static_cast<EVP_MD_CTX *>(ctx_)); }
void Hasher::flush_if_large(std::size_t threshold) {
  if (w.bytes.size() < threshold)
    return;
  need(EVP_DigestUpdate(static_cast<EVP_MD_CTX *>(ctx_), w.bytes.data(),
                        w.bytes.size()) == 1,
       "digest update");
  w.bytes.clear();
}
Digest Hasher::finish() {
  flush_if_large(0);
  Digest d{};
  unsigned n = 0;
  need(EVP_DigestFinal_ex(static_cast<EVP_MD_CTX *>(ctx_), d.data(), &n) == 1 &&
           n == 32,
       "digest final");
  return d;
}
Digest carrier_leaf(const GapParticle *rows, std::uint64_t block,
                    std::uint64_t count) {
  Hasher h;
  h.w.tag("native-v5/carrier-block");
  h.w.u64(block);
  h.w.u64(count);
  h.flush_if_large(0);
  h.w.raw(rows, count * sizeof(GapParticle));
  return h.finish();
}
std::string merkle_root(const std::string &root_tag,
                        const std::vector<BlockLeaf> &leaves,
                        std::uint64_t records) {
  Digest tree{};
  if (leaves.empty()) {
    Writer e;
    e.tag("native-v5/merkle-empty");
    tree = sha256_raw(e.bytes);
  } else {
    std::vector<Digest> level;
    level.reserve(leaves.size());
    for (std::size_t i = 0; i < leaves.size(); ++i) {
      need(i == 0 || leaves[i].block > leaves[i - 1].block,
           "Merkle leaves in strictly ascending block order");
      level.push_back(leaves[i].hash);
    }
    while (level.size() > 1) {
      std::vector<Digest> up;
      up.reserve((level.size() + 1) / 2);
      for (std::size_t i = 0; i + 1 < level.size(); i += 2)
        up.push_back(node(level[i], level[i + 1]));
      if (level.size() % 2)
        up.push_back(level.back());
      level.swap(up);
    }
    tree = level.front();
  }
  Writer w;
  w.tag(root_tag);
  w.u64(carrier_block);
  w.u64(records);
  w.u64(leaves.size());
  w.raw(tree.data(), tree.size());
  return sha256_hex(w.bytes);
}
std::shared_ptr<const BlockDigest>
carrier_rows_digest(const std::vector<GapParticle> &rows,
                    const std::vector<GapParticle> *parent,
                    const BlockDigest *parent_digest,
                    std::uint64_t *rehashed_blocks) {
  auto out = std::make_shared<BlockDigest>();
  out->records = rows.size();
  // Parent block lookup (both sequences are block-ascending).
  std::size_t pi = 0;
  const bool reuse = parent && parent_digest &&
                     parent_digest->records == parent->size();
  std::uint64_t rehashed = 0;
  std::size_t i = 0;
  while (i < rows.size()) {
    need(rows[i].id > 0, "native carrier ID zero");
    const std::uint64_t block = (rows[i].id - 1) / carrier_block;
    std::size_t j = i + 1;
    while (j < rows.size() && (rows[j].id - 1) / carrier_block == block) {
      need(rows[j].id > rows[j - 1].id,
           "native carrier rows must be in strictly ascending ID order");
      ++j;
    }
    need(j == rows.size() || rows[j].id > rows[j - 1].id,
         "native carrier rows must be in strictly ascending ID order");
    BlockLeaf leaf{block, i, j - i, {}};
    bool found = false;
    if (reuse) {
      const auto &pl = parent_digest->leaves;
      while (pi < pl.size() && pl[pi].block < block)
        ++pi;
      if (pi < pl.size() && pl[pi].block == block &&
          pl[pi].count == leaf.count &&
          std::memcmp(parent->data() + pl[pi].first, rows.data() + i,
                      leaf.count * sizeof(GapParticle)) == 0) {
        leaf.hash = pl[pi].hash;
        found = true;
      }
    }
    if (!found) {
      leaf.hash = carrier_leaf(rows.data() + i, block, leaf.count);
      ++rehashed;
    }
    out->leaves.push_back(leaf);
    i = j;
  }
  out->root = merkle_root("native-v5/carrier-root", out->leaves, rows.size());
  if (rehashed_blocks)
    *rehashed_blocks = rehashed;
  return out;
}
bool verify_carrier_leaves(const std::vector<GapParticle> &rows,
                           const std::vector<BlockLeaf> &leaves,
                           const std::string &root,
                           const std::vector<std::size_t> &sample) {
  std::uint64_t covered = 0;
  for (std::size_t k = 0; k < leaves.size(); ++k) {
    const auto &l = leaves[k];
    if (l.first != covered || l.count == 0 || l.first + l.count > rows.size() ||
        (rows[l.first].id - 1) / carrier_block != l.block ||
        (rows[l.first + l.count - 1].id - 1) / carrier_block != l.block)
      return false;
    covered += l.count;
  }
  if (covered != rows.size() || merkle_root("native-v5/carrier-root", leaves,
                                            rows.size()) != root)
    return false;
  for (auto k : sample) {
    if (k >= leaves.size())
      return false;
    const auto &l = leaves[k];
    if (carrier_leaf(rows.data() + l.first, l.block, l.count) != l.hash)
      return false;
  }
  return true;
}
} // namespace raikou::mfem_backend::canonical
