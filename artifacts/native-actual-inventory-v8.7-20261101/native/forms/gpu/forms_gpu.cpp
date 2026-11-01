// forms_gpu.cpp -- host launcher of the native v7 producer kernels.
#include "forms_gpu.hpp"
#include <cstring>
#include <openssl/sha.h>

namespace pf_v7::gpu {
namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kKernels[] = {
#embed "forms_kernels.cu"
    , 0};
const char kCore[] = {
#embed "../include/pf_v7/forms_core.h"
    , 0};
#pragma GCC diagnostic pop
using u64 = std::uint64_t;
unsigned grid_of(u64 n) {
  const u64 g = (n + Forms::block - 1) / Forms::block;
  if (g > 0x7fffffffu)
    throw hall_gpu::Error("v7 forms launch exceeds the grid range");
  return unsigned(g ? g : 1);
}
std::size_t align256(std::size_t b) { return (b + 255) & ~std::size_t(255); }
} // namespace

struct Forms::Impl {
  hall_gpu::Session &s;
  hall_gpu::Module module;
  std::unique_ptr<hall_gpu::Primitives> prim;
  hall_gpu::Kernel sum_products, kinetic_aligned, kinetic_split, sum_finish, chunks, blocks, tree, census_pack,
      census_check, census_finish, split_ids, gather_hi, census_gather;
  hall_gpu::Kernel chunks_multi, blocks_multi, tree_multi, id_flags, id_starts, id_leaves, id_leaves_pair, id_tree, sum_abs,
      node_marks, marked_max, marked_max_finish, stream, zero_bytes;
  hall_gpu::DeviceBuffer tags; // v5 ID-block tags: carrier-block, carrier-root, path-block, path-root
  std::uint64_t tag_addr[4] = {}, tag_len[4] = {};
  explicit Impl(hall_gpu::Session &session) : s(session) {
    hall_gpu::CompileRequest req;
    req.main = {"forms_kernels.cu", kKernels};
    req.headers = {{"forms_core.h", kCore}, hall_gpu::Primitives::device_source()};
    module = s.compile_and_load(req);
    prim = std::make_unique<hall_gpu::Primitives>(s, module);
    auto k = [&](const char *name) { return s.kernel(module, name); };
    sum_products = k("pf7_sum_products");
    kinetic_aligned = k("pf7_kinetic_aligned");
    kinetic_split = k("pf7_kinetic_split");
    sum_finish = k("pf7_sum_finish");
    chunks = k("pf7_chunks");
    blocks = k("pf7_blocks");
    tree = k("pf7_tree");
    census_pack = k("pf7_census_pack");
    census_check = k("pf7_census_check");
    census_finish = k("pf7_census_finish");
    split_ids = k("pf7_census_split_ids");
    gather_hi = k("pf7_gather_hi");
    census_gather = k("pf7_census_gather");
    chunks_multi = k("pf7_chunks_multi");
    blocks_multi = k("pf7_blocks_multi");
    tree_multi = k("pf7_tree_multi");
    id_flags = k("pf7_id_flags");
    id_starts = k("pf7_id_starts");
    id_leaves = k("pf7_id_leaves");
    id_leaves_pair = k("pf7_id_leaves_pair");
    id_tree = k("pf7_id_tree");
    sum_abs = k("pf7_sum_abs_products");
    node_marks = k("pf7_node_marks");
    marked_max = k("pf7_marked_max");
    marked_max_finish = k("pf7_marked_max_finish");
    stream = k("pf7_stream_digest");
    zero_bytes = k("pf7_zero_bytes");
    const char *names[4] = {"native-v5/carrier-block", "native-v5/carrier-root", "native-v5/path-block",
                            "native-v5/path-root"};
    std::string all;
    for (int i = 0; i < 4; ++i) {
      tag_addr[i] = all.size();
      tag_len[i] = std::strlen(names[i]);
      all += names[i];
    }
    tags = s.alloc(all.size(), "v7.1-id-tags");
    auto h = s.alloc_pinned(all.size(), "v7.1-id-tags");
    std::memcpy(h.pointer, all.data(), all.size());
    s.upload(tags, 0, h, 0, all.size());
    s.synchronize("v7.1 tag upload");
    s.release(h);
    for (auto &a : tag_addr)
      a += tags.address;
  }
  void finish(unsigned parts, u64 scratch, u64 terms, double scale, u64 out) {
    u64 partial = scratch, status = scratch + align256(std::size_t(parts) * 8 * pf7::acc_digits);
    void *a[] = {&parts, &partial, &status, &terms, &scale, &out};
    s.launch(sum_finish, {1}, {block}, 0, a);
  }
};

Forms::Forms(hall_gpu::Session &s) : p_(std::make_unique<Impl>(s)) {}
Forms::~Forms() = default;
const hall_gpu::CompileRecord &Forms::compile_record() const { return p_->module.record; }
hall_gpu::Primitives &Forms::primitives() { return *p_->prim; }
const std::string &Forms::kernel_source_sha256() {
  static const std::string h = [] {
    std::string both = std::string(kKernels) + std::string(1, '\0') + std::string(kCore);
    return hall_gpu::sha256_hex(both.data(), both.size());
  }();
  return h;
}

std::size_t Forms::sum_scratch_bytes(unsigned blocks) {
  return align256(std::size_t(blocks) * 8 * pf7::acc_digits) + align256(std::size_t(blocks) * 4);
}
void Forms::sum_products(u64 n, u64 a, u64 a_stride, u64 b, u64 b_stride, double scale, unsigned blocks,
                         u64 scratch, u64 out) {
  if (!blocks)
    throw hall_gpu::Error("v7 exact sum needs at least one block");
  u64 partial = scratch, status = scratch + align256(std::size_t(blocks) * 8 * pf7::acc_digits);
  void *args[] = {&n, &a, &a_stride, &b, &b_stride, &partial, &status};
  p_->s.launch(p_->sum_products, {blocks}, {block}, 0, args);
  p_->finish(blocks, scratch, n, scale, out);
}
void Forms::kinetic_aligned(u64 n, u64 w, u64 ws, u64 m, u64 ms, u64 u0, u64 u0s, u64 u1, u64 u1s, unsigned blocks,
                            u64 scratch, u64 out) {
  if (!blocks)
    throw hall_gpu::Error("v7 exact sum needs at least one block");
  u64 partial = scratch, status = scratch + align256(std::size_t(blocks) * 8 * pf7::acc_digits);
  void *args[] = {&n, &w, &ws, &m, &ms, &u0, &u0s, &u1, &u1s, &partial, &status};
  p_->s.launch(p_->kinetic_aligned, {blocks}, {block}, 0, args);
  p_->finish(blocks, scratch, n, 1.0, out);
}
void Forms::kinetic_split(u64 ne, u64 endpoints, u64 ns, u64 sinks, u64 nc, u64 census, u64 census_words,
                          unsigned blocks, u64 scratch, u64 out) {
  if (!blocks)
    throw hall_gpu::Error("v7 exact sum needs at least one block");
  u64 partial = scratch, status = scratch + align256(std::size_t(blocks) * 8 * pf7::acc_digits);
  void *args[] = {&ne, &endpoints, &ns, &sinks, &nc, &census, &census_words, &partial, &status};
  p_->s.launch(p_->kinetic_split, {blocks}, {block}, 0, args);
  p_->finish(blocks, scratch, ne + ns, 1.0, out);
}

// Block list scratch: chunk digests | leaf digests | two tree levels.
std::size_t Forms::leaf_offset(u64 n) {
  const u64 chunks = (n + pf7::chunk_records - 1) / pf7::chunk_records;
  return align256(std::size_t(chunks ? chunks : 1) * 32);
}
std::size_t Forms::block_list_scratch_bytes(u64 n) {
  const u64 leaves = (n + pf7::leaf_records - 1) / pf7::leaf_records;
  const std::size_t level = align256(std::size_t((leaves + 1) / 2 + 1) * 32);
  return leaf_offset(n) + align256(std::size_t(leaves ? leaves : 1) * 32) + 2 * level;
}
void Forms::block_list(std::uint32_t section, u64 n, u64 records, u64 scratch, u64 root_out) {
  if (section >= pf7::section_count)
    throw hall_gpu::Error("unknown v7 block-list section");
  const u64 chunks = (n + pf7::chunk_records - 1) / pf7::chunk_records;
  u64 leaves = (n + pf7::leaf_records - 1) / pf7::leaf_records;
  u64 chunk_buf = scratch, leaf_buf = scratch + leaf_offset(n);
  const std::size_t level = align256(std::size_t((leaves + 1) / 2 + 1) * 32);
  u64 a = leaf_buf + align256(std::size_t(leaves ? leaves : 1) * 32), b = a + level;
  if (chunks) {
    void *c[] = {&section, &n, &records, &chunk_buf};
    p_->s.launch(p_->chunks, {grid_of(chunks)}, {block}, 0, c);
    void *l[] = {&section, &n, &chunk_buf, &leaf_buf};
    p_->s.launch(p_->blocks, {grid_of(leaves)}, {block}, 0, l);
  }
  void *t[] = {&section, &n, &leaves, &leaf_buf, &a, &b, &root_out};
  p_->s.launch(p_->tree, {1}, {block}, 0, t);
}

void Forms::census_pack(u64 n, u64 rows, u64 row_words, u64 keys, u64 states, u64 out) {
  if (!n)
    return;
  void *a[] = {&n, &rows, &row_words, &keys, &states, &out};
  p_->s.launch(p_->census_pack, {grid_of(n)}, {block}, 0, a);
}
std::size_t Forms::census_check_scratch_bytes(u64 n) {
  const u64 g = grid_of(n);
  return align256(g * 8) + align256(g * 4);
}
void Forms::census_check(u64 n, u64 records, u64 states, u64 scratch, u64 out) {
  u64 g = grid_of(n);
  u64 idx = scratch, code = scratch + align256(g * 8);
  void *a[] = {&n, &records, &states, &idx, &code};
  p_->s.launch(p_->census_check, {unsigned(g)}, {block}, 0, a);
  void *f[] = {&g, &idx, &code, &out};
  p_->s.launch(p_->census_finish, {1}, {32}, 0, f);
}
std::size_t Forms::census_order_scratch_bytes(u64 n) {
  const std::size_t v = align256(std::size_t(n ? n : 1) * 4);
  return 5 * v + align256(hall_gpu::Primitives::sort_scratch(n ? n : 1) * 4);
}
void Forms::census_order(u64 n, u64 records, u64 rec_words, u64 scratch, u64 out) {
  if (!n)
    return;
  if (n >= (u64(1) << 31))
    throw hall_gpu::Error("v7 census order: more than 2^31 records");
  const std::size_t v = align256(std::size_t(n) * 4);
  u64 lo = scratch, hi = lo + v, pay_a = hi + v, key_b = pay_a + v, pay_b = key_b + v, sort = pay_b + v;
  auto &P = *p_->prim;
  void *s1[] = {&n, &records, &rec_words, &lo, &hi};
  p_->s.launch(p_->split_ids, {grid_of(n)}, {block}, 0, s1);
  P.iota(pay_a, n);
  // Pass 1: low halves (stable); 4 byte passes end in the a buffers.
  bool in_b = P.radix_sort(lo, pay_a, key_b, pay_b, n, 32, sort);
  u64 perm = in_b ? pay_b : pay_a;
  // Pass 2: high halves in pass-1 order (stable).
  u64 keys2 = in_b ? lo : key_b;
  void *g1[] = {&n, &hi, &perm, &keys2};
  p_->s.launch(p_->gather_hi, {grid_of(n)}, {block}, 0, g1);
  u64 other_key = keys2 == lo ? key_b : lo, other_pay = perm == pay_a ? pay_b : pay_a;
  in_b = P.radix_sort(keys2, perm, other_key, other_pay, n, 32, sort);
  u64 final_perm = in_b ? other_pay : perm;
  void *g2[] = {&n, &records, &rec_words, &final_perm, &out};
  p_->s.launch(p_->census_gather, {grid_of(n)}, {block}, 0, g2);
}

// ---------------------------------------------------------------- v7.1 additions
namespace {
struct MultiJobs { // layout of struct MultiJobs in forms_kernels.cu
  std::uint32_t count;
  std::uint32_t section[8];
  std::uint64_t n[8];
  std::uint64_t records[8], chunks[8], leaves[8], tree_a[8], tree_b[8], root[8];
  std::uint64_t chunk_start[9], leaf_start[9];
};
static_assert(sizeof(MultiJobs) == 40 + 64 + 6 * 64 + 2 * 72);
struct JobScratch {
  std::size_t chunks = 0, leaves = 0, level = 0, total = 0;
};
JobScratch job_scratch(u64 n) {
  const u64 chunks = (n + pf7::chunk_records - 1) / pf7::chunk_records;
  const u64 leaves = (n + pf7::leaf_records - 1) / pf7::leaf_records;
  JobScratch j;
  j.chunks = align256(std::size_t(chunks ? chunks : 1) * 32);
  j.leaves = align256(std::size_t(leaves ? leaves : 1) * 32);
  j.level = align256(std::size_t((leaves + 1) / 2 + 1) * 32);
  j.total = j.chunks + j.leaves + 2 * j.level;
  return j;
}
} // namespace
std::size_t Forms::block_lists_scratch_bytes(std::span<const BlockListJob> jobs) {
  std::size_t t = 0;
  for (const auto &j : jobs)
    t += job_scratch(j.n).total;
  return t;
}
void Forms::block_lists(std::span<const BlockListJob> jobs, u64 scratch) {
  if (jobs.empty())
    return;
  if (jobs.size() > 8)
    throw hall_gpu::Error("v7.1 fused block lists: at most 8 sections");
  MultiJobs J{};
  J.count = std::uint32_t(jobs.size());
  u64 at = scratch;
  for (std::size_t k = 0; k < jobs.size(); ++k) {
    const auto &j = jobs[k];
    if (j.section >= pf7::section_count)
      throw hall_gpu::Error("unknown v7 block-list section");
    const auto js = job_scratch(j.n);
    J.section[k] = j.section;
    J.n[k] = j.n;
    J.records[k] = j.records;
    J.chunks[k] = at;
    J.leaves[k] = at + js.chunks;
    J.tree_a[k] = J.leaves[k] + js.leaves;
    J.tree_b[k] = J.tree_a[k] + js.level;
    J.root[k] = j.root_out;
    at += js.total;
    J.chunk_start[k + 1] = J.chunk_start[k] + (j.n + pf7::chunk_records - 1) / pf7::chunk_records;
    J.leaf_start[k + 1] = J.leaf_start[k] + (j.n + pf7::leaf_records - 1) / pf7::leaf_records;
  }
  void *a[] = {&J};
  if (J.chunk_start[J.count]) {
    p_->s.launch(p_->chunks_multi, {grid_of(J.chunk_start[J.count])}, {block}, 0, a);
    p_->s.launch(p_->blocks_multi, {grid_of(J.leaf_start[J.count])}, {block}, 0, a);
  }
  p_->s.launch(p_->tree_multi, {J.count}, {block}, 0, a);
}
// ID-block root scratch: flags | scan | starts | scan scratch | total | leaves | 2 tree levels.
std::size_t Forms::id_root_count_offset(u64 n) {
  const std::size_t v = align256(std::size_t(n ? n : 1) * 4);
  return 3 * v + align256(hall_gpu::Primitives::scan_scratch(n ? n : 1) * 4);
}
std::size_t Forms::id_root_scratch_bytes(u64 n) {
  const std::size_t leaves = std::size_t(n ? n : 1);
  return id_root_count_offset(n) + 256 + align256(leaves * 32) + 2 * align256((leaves + 1) / 2 * 32 + 32);
}
void Forms::id_root(std::uint32_t kind, u64 n, u64 records, u64 words, u64 scratch, u64 root_out) {
  if (kind > 1)
    throw hall_gpu::Error("v7.1 ID-block root kind");
  if (n >= (u64(1) << 31))
    throw hall_gpu::Error("v7.1 ID-block root: more than 2^31 records");
  const std::size_t v = align256(std::size_t(n ? n : 1) * 4);
  u64 flags = scratch, scan = flags + v, starts = scan + v, scan_scratch = starts + v;
  u64 count = scratch + id_root_count_offset(n);
  u64 leaves = count + 256;
  const std::size_t L = std::size_t(n ? n : 1);
  u64 ta = leaves + align256(L * 32), tb = ta + align256((L + 1) / 2 * 32 + 32);
  u64 leaf_tag = p_->tag_addr[2 * kind], root_tag = p_->tag_addr[2 * kind + 1];
  std::uint32_t leaf_len = std::uint32_t(p_->tag_len[2 * kind]), root_len = std::uint32_t(p_->tag_len[2 * kind + 1]);
  if (n) {
    void *f[] = {&n, &records, &words, &flags};
    p_->s.launch(p_->id_flags, {grid_of(n)}, {block}, 0, f);
    p_->prim->exclusive_scan(flags, n, scan, scan_scratch, count);
    void *st[] = {&n, &flags, &scan, &starts};
    p_->s.launch(p_->id_starts, {grid_of(n)}, {block}, 0, st);
    void *l[] = {&n, &records, &words, &starts, &count, &leaf_tag, &leaf_len, &leaves};
    // Warp-sized blocks: the leaf threads (at most n / 1024 + 1 of them, each a long hash) spread over the SMs.
    p_->s.launch(p_->id_leaves, {unsigned((n + 31) / 32)}, {32}, 0, l);
  } else {
    u64 four = 4;
    void *z[] = {&four, &count};
    p_->s.launch(p_->zero_bytes, {1}, {block}, 0, z);
  }
  void *t[] = {&n, &count, &leaves, &ta, &tb, &root_tag, &root_len, &root_out};
  p_->s.launch(p_->id_tree, {1}, {block}, 0, t);
}
std::size_t Forms::id_root_pair_scratch_bytes(u64 n) {
  const std::size_t leaves = std::size_t(n ? n : 1);
  return id_root_scratch_bytes(n) + align256(leaves * 32) + 2 * align256((leaves + 1) / 2 * 32 + 32);
}
void Forms::id_root_pair(u64 n, u64 recA, u64 wordsA, u64 recB, u64 wordsB, u64 scratch, u64 rootA, u64 rootB) {
  if (n >= (u64(1) << 31))
    throw hall_gpu::Error("v7.1 ID-block root: more than 2^31 records");
  const std::size_t v = align256(std::size_t(n ? n : 1) * 4);
  u64 flags = scratch, scan = flags + v, starts = scan + v, scan_scratch = starts + v;
  u64 count = scratch + id_root_count_offset(n);
  const std::size_t L = std::size_t(n ? n : 1);
  u64 leavesA = count + 256;
  u64 taA = leavesA + align256(L * 32), tbA = taA + align256((L + 1) / 2 * 32 + 32);
  u64 leavesB = scratch + id_root_scratch_bytes(n);
  u64 taB = leavesB + align256(L * 32), tbB = taB + align256((L + 1) / 2 * 32 + 32);
  u64 tagA = p_->tag_addr[0], rtagA = p_->tag_addr[1], tagB = p_->tag_addr[2], rtagB = p_->tag_addr[3];
  std::uint32_t lenA = std::uint32_t(p_->tag_len[0]), rlenA = std::uint32_t(p_->tag_len[1]),
                lenB = std::uint32_t(p_->tag_len[2]), rlenB = std::uint32_t(p_->tag_len[3]);
  if (n) {
    void *f[] = {&n, &recA, &wordsA, &flags};
    p_->s.launch(p_->id_flags, {grid_of(n)}, {block}, 0, f);
    p_->prim->exclusive_scan(flags, n, scan, scan_scratch, count);
    void *st[] = {&n, &flags, &scan, &starts};
    p_->s.launch(p_->id_starts, {grid_of(n)}, {block}, 0, st);
    void *l[] = {&n, &recA, &wordsA, &recB, &wordsB, &starts, &count, &tagA, &lenA, &tagB, &lenB, &leavesA, &leavesB};
    p_->s.launch(p_->id_leaves_pair, {unsigned(2 * ((n + 31) / 32))}, {32}, 0, l);
  } else {
    u64 four = 4;
    void *z[] = {&four, &count};
    p_->s.launch(p_->zero_bytes, {1}, {block}, 0, z);
  }
  void *ta[] = {&n, &count, &leavesA, &taA, &tbA, &rtagA, &rlenA, &rootA};
  p_->s.launch(p_->id_tree, {1}, {block}, 0, ta);
  void *tb[] = {&n, &count, &leavesB, &taB, &tbB, &rtagB, &rlenB, &rootB};
  p_->s.launch(p_->id_tree, {1}, {block}, 0, tb);
}
void Forms::sum_abs_products(u64 n, u64 a, u64 a_stride, u64 b, u64 b_stride, double scale, unsigned blocks,
                             u64 scratch, u64 out) {
  if (!blocks)
    throw hall_gpu::Error("v7 exact sum needs at least one block");
  u64 partial = scratch, status = scratch + align256(std::size_t(blocks) * 8 * pf7::acc_digits);
  void *args[] = {&n, &a, &a_stride, &b, &b_stride, &partial, &status};
  p_->s.launch(p_->sum_abs, {blocks}, {block}, 0, args);
  p_->finish(blocks, scratch, n, scale, out);
}
std::size_t Forms::field_bound_scratch_bytes(u64 ne, u64 nn) {
  const u64 g = grid_of(ne);
  return align256(nn ? nn : 1) + align256(g * 8) + align256(g * 4);
}
void Forms::field_bound(u64 ne, u64 E, u64 J, u64 n, u64 nz, u64 scratch, u64 out) {
  const u64 nn = (n + 1) * (n + 1) * (nz + 1);
  u64 g = grid_of(ne);
  u64 marks = scratch, bmax = marks + align256(nn ? nn : 1), bbad = bmax + align256(g * 8);
  {
    u64 bytes = nn;
    void *z[] = {&bytes, &marks};
    p_->s.launch(p_->zero_bytes, {grid_of(nn)}, {block}, 0, z);
  }
  void *m[] = {&ne, &J, &n, &nz, &marks};
  p_->s.launch(p_->node_marks, {unsigned(g)}, {block}, 0, m);
  void *x[] = {&ne, &E, &n, &nz, &marks, &bmax, &bbad};
  p_->s.launch(p_->marked_max, {unsigned(g)}, {block}, 0, x);
  u64 omax = out, obad = out + 8;
  void *f[] = {&g, &bmax, &bbad, &omax, &obad};
  p_->s.launch(p_->marked_max_finish, {1}, {32}, 0, f);
}
void Forms::stream_digest(std::uint32_t kind, u64 ne, u64 E, u64 nb, u64 B, u64 out) {
  void *a[] = {&kind, &ne, &E, &nb, &B, &out};
  p_->s.launch(p_->stream, {1}, {32}, 0, a);
}

// ---------------------------------------------------------------- HostRunner
namespace {
struct Scope {
  hall_gpu::Session &s;
  std::vector<hall_gpu::DeviceBuffer> dev;
  std::vector<hall_gpu::HostBuffer> host;
  explicit Scope(hall_gpu::Session &x) : s(x) {}
  ~Scope() {
    for (auto &d : dev)
      if (d)
        s.release(d);
    for (auto &h : host)
      if (h)
        s.release(h);
  }
  hall_gpu::DeviceBuffer alloc(std::size_t bytes, const char *tag) {
    dev.push_back(s.alloc(bytes ? bytes : 8, tag));
    return dev.back();
  }
  hall_gpu::DeviceBuffer put(const void *data, std::size_t bytes, const char *tag) {
    auto d = alloc(bytes, tag);
    if (bytes) {
      host.push_back(s.alloc_pinned(bytes, tag));
      std::memcpy(host.back().pointer, data, bytes);
      s.upload(d, 0, host.back(), 0, bytes);
    }
    return d;
  }
  void get(void *data, const hall_gpu::DeviceBuffer &d, std::size_t offset, std::size_t bytes) {
    host.push_back(s.alloc_pinned(bytes, "v7-download"));
    s.download(host.back(), 0, d, offset, bytes);
    s.synchronize("v7 forms download");
    std::memcpy(data, host.back().pointer, bytes);
  }
};
} // namespace
SumOut HostRunner::sum_products(std::span<const double> a, std::span<const double> b, double scale,
                                unsigned blocks) {
  Scope sc(s_);
  auto da = sc.put(a.data(), a.size_bytes(), "v7-a");
  auto db = sc.put(b.data(), b.size_bytes(), "v7-b");
  auto scr = sc.alloc(Forms::sum_scratch_bytes(blocks), "v7-sum-scratch");
  auto out = sc.alloc(sizeof(SumOut), "v7-sum-out");
  s_.synchronize("v7 upload");
  s_.clear_markers();
  const int m0 = s_.mark("sum");
  f_.sum_products(a.size(), da.address, 1, db.address, 1, scale, blocks, scr.address, out.address);
  const int m1 = s_.mark("sum-end");
  SumOut r;
  sc.get(&r, out, 0, sizeof r);
  last_device_ms = s_.elapsed_ms(m0, m1);
  return r;
}
SumOut HostRunner::kinetic_aligned(std::span<const Carrier> census, std::span<const double> u_end,
                                   unsigned blocks) {
  Scope sc(s_);
  auto dc = sc.put(census.data(), census.size_bytes(), "v7-census");
  auto du = sc.put(u_end.data(), u_end.size_bytes(), "v7-u-end");
  auto scr = sc.alloc(Forms::sum_scratch_bytes(blocks), "v7-sum-scratch");
  auto out = sc.alloc(sizeof(SumOut), "v7-sum-out");
  s_.synchronize("v7 upload");
  s_.clear_markers();
  const int m0 = s_.mark("kinetic");
  const u64 base = dc.address;
  f_.kinetic_aligned(census.size(), base + 72, 10, base + 64, 10, base + 32, 10, du.address, 3, blocks, scr.address,
                     out.address);
  const int m1 = s_.mark("kinetic-end");
  SumOut r;
  sc.get(&r, out, 0, sizeof r);
  last_device_ms = s_.elapsed_ms(m0, m1);
  return r;
}
SumOut HostRunner::kinetic_split(std::span<const Carrier> census, std::span<const Carrier> endpoints,
                                 std::span<const Sink> sinks, unsigned blocks) {
  Scope sc(s_);
  auto dc = sc.put(census.data(), census.size_bytes(), "v7-census");
  auto de = sc.put(endpoints.data(), endpoints.size_bytes(), "v7-endpoints");
  auto ds = sc.put(sinks.data(), sinks.size_bytes(), "v7-sinks");
  auto scr = sc.alloc(Forms::sum_scratch_bytes(blocks), "v7-sum-scratch");
  auto out = sc.alloc(sizeof(SumOut), "v7-sum-out");
  s_.synchronize("v7 upload");
  s_.clear_markers();
  const int m0 = s_.mark("kinetic");
  f_.kinetic_split(endpoints.size(), de.address, sinks.size(), ds.address, census.size(), dc.address, 10, blocks,
                   scr.address, out.address);
  const int m1 = s_.mark("kinetic-end");
  SumOut r;
  sc.get(&r, out, 0, sizeof r);
  last_device_ms = s_.elapsed_ms(m0, m1);
  return r;
}
HostRunner::BlockListOut HostRunner::block_list(std::uint32_t section, const void *records, u64 n) {
  Scope sc(s_);
  auto dr = sc.put(records, n * pf7::section_record_bytes(section), "v7-records");
  auto scr = sc.alloc(Forms::block_list_scratch_bytes(n), "v7-block-scratch");
  auto root = sc.alloc(32, "v7-root");
  s_.synchronize("v7 upload");
  s_.clear_markers();
  const int m0 = s_.mark("block-list");
  f_.block_list(section, n, dr.address, scr.address, root.address);
  const int m1 = s_.mark("block-list-end");
  BlockListOut out;
  std::uint32_t w[8];
  sc.get(w, root, 0, 32);
  last_device_ms = s_.elapsed_ms(m0, m1);
  out.root = digest_of_words(w);
  const u64 leaves = (n + pf7::leaf_records - 1) / pf7::leaf_records;
  std::vector<std::uint32_t> lw(8 * leaves);
  if (leaves)
    sc.get(lw.data(), scr, Forms::leaf_offset(n), lw.size() * 4);
  for (u64 k = 0; k < leaves; ++k)
    out.leaves.push_back(digest_of_words(&lw[8 * k]));
  return out;
}
CensusCheck HostRunner::census_check(std::span<const CensusRecord> c, std::uint64_t states) {
  Scope sc(s_);
  auto dr = sc.put(c.data(), c.size_bytes(), "v7-census");
  auto scr = sc.alloc(Forms::census_check_scratch_bytes(c.size()), "v7-census-scratch");
  auto out = sc.alloc(16, "v7-census-out");
  s_.synchronize("v7 upload");
  s_.clear_markers();
  const int m0 = s_.mark("census");
  f_.census_check(c.size(), dr.address, states, scr.address, out.address);
  const int m1 = s_.mark("census-end");
  u64 r[2];
  sc.get(r, out, 0, 16);
  last_device_ms = s_.elapsed_ms(m0, m1);
  CensusCheck k;
  if (r[0] != ~u64(0)) {
    k.index = r[0];
    k.code = std::uint32_t(r[1]);
  }
  return k;
}
std::vector<CensusRecord> HostRunner::census_order(std::span<const CensusRecord> c) {
  Scope sc(s_);
  auto dr = sc.put(c.data(), c.size_bytes(), "v7-census");
  auto scr = sc.alloc(Forms::census_order_scratch_bytes(c.size()), "v7-order-scratch");
  auto out = sc.alloc(c.size_bytes(), "v7-census-ordered");
  s_.synchronize("v7 upload");
  s_.clear_markers();
  const int m0 = s_.mark("order");
  f_.census_order(c.size(), dr.address, 12, scr.address, out.address);
  const int m1 = s_.mark("order-end");
  std::vector<CensusRecord> r(c.size());
  if (!c.empty())
    sc.get(r.data(), out, 0, c.size_bytes());
  else
    s_.synchronize("v7 order");
  last_device_ms = s_.elapsed_ms(m0, m1);
  return r;
}
} // namespace pf_v7::gpu
