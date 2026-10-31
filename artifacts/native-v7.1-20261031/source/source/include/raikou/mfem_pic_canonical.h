#pragma once
// Canonical binary encodings for native v5 full-time identities.
//
// Primitive rules (every record):
//  * u64: 8 bytes little-endian. i64: two's-complement value as u64.
//  * f64: the IEEE-754 binary64 bit pattern as u64 (exact; -0.0 != +0.0).
//  * bool: one byte, 0x00 or 0x01 (any other byte refuses on decode).
//  * text: u64 byte length, then the bytes (no terminator).
//  * every record starts with its own version tag encoded as text, so two
//    records of different kinds or versions never share an encoding.
//  * sequences: u64 element count, then elements in a stated canonical order.
// Hashes are SHA-256 over these bytes; hex digests are lower-case.
#include "raikou/mfem_gap_pic.h"
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace raikou::mfem_backend::canonical {
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);

class EncodingRefusal : public std::invalid_argument {
public:
  using std::invalid_argument::invalid_argument;
};

class Writer {
public:
  std::string bytes;
  void u8(const std::uint8_t &v) { bytes.push_back(char(v)); }
  void u64(const std::uint64_t &v) {
    char b[8];
    for (int i = 0; i < 8; ++i)
      b[i] = char((v >> (8 * i)) & 255);
    bytes.append(b, 8);
  }
  void i64(const std::int64_t &v) { u64(static_cast<std::uint64_t>(v)); }
  void f64(const double &x) { u64(std::bit_cast<std::uint64_t>(x)); }
  void boolean(const bool &b) { u8(b ? 1 : 0); }
  void text(const std::string &s, std::size_t = 0) {
    u64(s.size());
    bytes.append(s);
  }
  void tag(std::string_view s) {
    u64(s.size());
    bytes.append(s);
  }
  void vec(const PicVector &v) {
    for (double x : v)
      f64(x);
  }
  void raw(const void *p, std::size_t n) {
    bytes.append(static_cast<const char *>(p), n);
  }
};

// Reads a caller-owned byte sequence; it never copies, so constructing it from
// a temporary string is rejected at compile time (no dangling view).
class Reader {
public:
  explicit Reader(std::string_view b) : b_(b) {}
  explicit Reader(const std::string &b) : b_(b) {}
  Reader(std::string &&) = delete;
  void u8(std::uint8_t &v) { v = static_cast<std::uint8_t>(take(1)[0]); }
  void u64(std::uint64_t &v) {
    const auto s = take(8);
    v = 0;
    for (int i = 0; i < 8; ++i)
      v |= std::uint64_t(static_cast<unsigned char>(s[i])) << (8 * i);
  }
  void i64(std::int64_t &v) {
    std::uint64_t u;
    u64(u);
    v = static_cast<std::int64_t>(u);
  }
  void f64(double &x) {
    std::uint64_t u;
    u64(u);
    x = std::bit_cast<double>(u);
  }
  void boolean(bool &b) {
    std::uint8_t v;
    u8(v);
    if (v > 1)
      throw EncodingRefusal("canonical boolean byte must be 0 or 1");
    b = v == 1;
  }
  void text(std::string &s, std::size_t maximum = 1u << 20) {
    std::uint64_t n;
    u64(n);
    if (n > maximum || n > remaining())
      throw EncodingRefusal("canonical text length exceeds bound/input");
    s.assign(take(std::size_t(n)));
  }
  void tag(std::string_view expected) {
    std::string actual;
    text(actual, 256);
    if (actual != expected)
      throw EncodingRefusal("canonical record version tag mismatch");
  }
  void vec(PicVector &v) {
    for (double &x : v)
      f64(x);
  }
  std::uint64_t count(std::uint64_t maximum) {
    std::uint64_t n;
    u64(n);
    if (n > maximum)
      throw EncodingRefusal("canonical sequence count exceeds bound");
    return n;
  }
  std::size_t remaining() const noexcept { return b_.size() - pos_; }
  void finish() const {
    if (pos_ != b_.size())
      throw EncodingRefusal("trailing bytes after canonical record");
  }

private:
  std::string_view take(std::size_t n) {
    if (n > remaining())
      throw EncodingRefusal("truncated canonical record");
    auto s = b_.substr(pos_, n);
    pos_ += n;
    return s;
  }
  std::string_view b_;
  std::size_t pos_ = 0;
};

using Digest = std::array<unsigned char, 32>;
Digest sha256_raw(std::string_view bytes);
std::string hex(const Digest &);
std::string sha256_hex(std::string_view bytes);
Digest from_hex(const std::string &);
bool is_digest(const std::string &) noexcept;

// Incremental SHA-256 for large records (bounded transient buffer).
class Hasher {
public:
  Hasher();
  ~Hasher();
  Hasher(const Hasher &) = delete;
  Hasher &operator=(const Hasher &) = delete;
  Writer w;
  void flush_if_large(std::size_t threshold = 1u << 16);
  Digest finish();

private:
  void *ctx_;
};

// Block Merkle digest over ID-sorted per-carrier records.
//  * Records are sorted by strictly increasing native ID (ID 0 refused).
//  * Block k holds IDs in [k*B + 1, (k+1)*B]; empty blocks are omitted.
//  * leaf(k) = SHA256(text leaf_tag, u64 k, u64 n, n canonical records).
//  * node(a,b) = SHA256(text "native-v5/merkle-node", a(32 raw), b(32 raw));
//    levels pair left to right; an odd last node is carried up unchanged.
//    Leaf and node tags differ, so a leaf can never be read as a node.
//  * root = SHA256(text root_tag, u64 B, u64 records, u64 leaves,
//    tree(32 raw) or, for no leaves, SHA256(text "native-v5/merkle-empty")).
// Leaves are independent: a device may compute them; the host recomputes the
// tree from leaves (O(records/B) hashes) and verifies a sample of leaves.
inline constexpr std::uint64_t carrier_block = 1024;
struct BlockLeaf {
  std::uint64_t block = 0, first = 0, count = 0; // first = record offset
  Digest hash{};
};
struct BlockDigest {
  std::string root;              // hex root commitment
  std::vector<BlockLeaf> leaves; // ascending block index
  std::uint64_t records = 0;
};
// Native carrier rows: [u64 id, 3 f64 position, 3 f64 proper velocity,
// f64 charge, f64 mass, f64 weight]. `parent` (optional) lets unchanged blocks
// reuse their leaf hash after a byte-for-byte comparison of the block.
std::shared_ptr<const BlockDigest>
carrier_rows_digest(const std::vector<GapParticle> &rows,
                    const std::vector<GapParticle> *parent = nullptr,
                    const BlockDigest *parent_digest = nullptr,
                    std::uint64_t *rehashed_blocks = nullptr);
Digest carrier_leaf(const GapParticle *rows, std::uint64_t block,
                    std::uint64_t count);
std::string merkle_root(const std::string &root_tag,
                        const std::vector<BlockLeaf> &leaves,
                        std::uint64_t records);
// Host verification of device-supplied leaves: recompute a sample of leaves
// from host rows and the root from all supplied leaves.
bool verify_carrier_leaves(const std::vector<GapParticle> &rows,
                           const std::vector<BlockLeaf> &leaves,
                           const std::string &root,
                           const std::vector<std::size_t> &sample);
} // namespace raikou::mfem_backend::canonical
