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
  if (section > pf7::section_census)
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
