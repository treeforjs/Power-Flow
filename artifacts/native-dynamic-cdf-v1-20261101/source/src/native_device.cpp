// native_device.cpp -- native v7.1 device lane (API-071.txt section 2).
// Runs on the caller's hall-gpu-resident-v1 Session (one stream) and pf_v7
// forms module; compiles native_kernels.cu with forms_core.h through NVRTC.
#include "native_device.hpp"
#include "forms_gpu.hpp"
#include "hall_gpu_common.hpp"
#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <map>
#include <mutex>
#include <stdexcept>

namespace raikou::mfem_backend {
namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kNativeKernels[] = {
#embed "native_kernels.cu"
    , 0};
const char kCore[] = {
#embed "../../forms/include/pf_v7/forms_core.h"
    , 0};
// Chemistry program and sealed RNG2/blob-v2 dynamic-CDF headers, read in place.
// Legacy blob-v1 loading remains explicit; model/configuration identities change.
const char kChem[] = {
#embed "chem_v3_kernels.cu"
    , 0};
const char kChemEvent[] = {
#embed "/home/treefor/projects/runtimes/chemistry-dynamic-cdf-v1.1-20261101/include/pf_chem_v3/chem_event.h"
    , 0};
const char kChemKin[] = {
#embed "/home/treefor/projects/runtimes/chemistry-dynamic-cdf-v1.1-20261101/include/pf_chem_v3/chem_kinematics.h"
    , 0};
const char kChemRng[] = {
#embed "/home/treefor/projects/runtimes/chemistry-dynamic-cdf-v1.1-20261101/include/pf_chem_v3/chem_rng.h"
    , 0};
const char kChemTables[] = {
#embed "/home/treefor/projects/runtimes/chemistry-dynamic-cdf-v1.1-20261101/include/pf_chem_v3/chem_tables.h"
    , 0};
const char kChemMath[] = {
#embed "/home/treefor/projects/runtimes/chemistry-dynamic-cdf-v1.1-20261101/include/pf_chem_v3/det_math.h"
    , 0};
const char kPhilox[] = {
#embed "/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival/philox_core.h"
    , 0};
#pragma GCC diagnostic pop
using u64 = std::uint64_t;
constexpr unsigned block = 256;
std::size_t align256(std::size_t b) { return (b + 255) & ~std::size_t(255); }
unsigned grid_of(u64 n) { return unsigned((n + block - 1) / block ? (n + block - 1) / block : 1); }
std::string hex_words(const unsigned char *bytes) {
  std::uint32_t w[8];
  std::memcpy(w, bytes, 32);
  return pf_v7::hex(pf_v7::digest_of_words(w));
}
// Device-side layouts (forms_kernels.cu SumOut, native_kernels.cu Geometry/AdmitArgs).
struct Geometry {
  double width, gap, short_z, feed_z;
  double lower[4], upper[4];
  double shift[3];
  double dt, neutral_bound;
  std::uint32_t finite_foils, finite_axial;
};
struct AdmitArgs {
  u64 n, base, ends, paths, keys;
  Geometry g;
  u64 k_w, k_d, g_w, g_d, c_w, c_k, q_a, block_fail, block_code;
  // v8: sinks and newborns
  u64 n_ends, sinks, n_sinks, births, n_births, states, after_keys, after_states, sunk_index, born_index, born_end;
};
} // namespace

// ---------------------------------------------------------------- pool
struct DevicePool {
  hall_gpu::Session &s;
  std::multimap<std::size_t, hall_gpu::DeviceBuffer> free;
  explicit DevicePool(hall_gpu::Session &x) : s(x) {}
  ~DevicePool() {
    for (auto &[bytes, b] : free) {
      (void)bytes;
      try {
        s.release(b);
      } catch (...) {
      }
    }
  }
  hall_gpu::DeviceBuffer take(std::size_t bytes) {
    const std::size_t want = std::max<std::size_t>(align256(bytes), 256);
    auto it = free.lower_bound(want);
    if (it != free.end() && it->first <= 2 * want + (std::size_t(1) << 20)) {
      auto b = it->second;
      free.erase(it);
      return b;
    }
    // Round to 64 KiB to improve reuse.
    const std::size_t cap = (want + 65535) & ~std::size_t(65535);
    return s.alloc(cap, "native-v7.1-pool");
  }
  void put(const hall_gpu::DeviceBuffer &b) { free.emplace(b.bytes, b); }
};
DeviceBlock::~DeviceBlock() {
  if (pool && address)
    pool->put(hall_gpu::DeviceBuffer{address, bytes, slot});
}
namespace {
std::shared_ptr<DeviceBlock> take(const std::shared_ptr<DevicePool> &pool, std::size_t bytes) {
  auto b = pool->take(bytes);
  auto d = std::make_shared<DeviceBlock>();
  d->pool = pool;
  d->address = b.address;
  d->bytes = b.bytes;
  d->slot = b.slot;
  return d;
}
hall_gpu::DeviceBuffer view(const DeviceBlock &b) { return {b.address, b.bytes, b.slot}; }
} // namespace

// ---------------------------------------------------------------- impl
struct FullTimeDevice::Impl {
  hall_gpu::Session &s;
  pf_v7::gpu::Forms &f;
  hall_gpu::Module module;
  hall_gpu::Kernel admit, fail_finish, copy, fill, gather;
  std::shared_ptr<DevicePool> pool;
  hall_gpu::HostBuffer staging;
  std::size_t staging_bytes = 0;
  FullTimeDeviceCounters counters;
  std::mutex m;
  // v8 chemistry module (a std::shared_ptr<ChemModule>, see chem_slot), declared
  // after the pool so that it is released first.
  std::shared_ptr<void> chem;
  Impl(hall_gpu::Session &session, pf_v7::gpu::Forms &forms) : s(session), f(forms) {
    hall_gpu::CompileRequest req;
    req.main = {"native_kernels.cu", kNativeKernels};
    req.headers = {{"forms_core.h", kCore}};
    module = s.compile_and_load(req);
    admit = s.kernel(module, "pf7n_admit");
    fail_finish = s.kernel(module, "pf7n_fail_finish");
    copy = s.kernel(module, "pf7n_copy");
    fill = s.kernel(module, "pf7n_fill");
    gather = s.kernel(module, "pf7n_gather");
    pool = std::make_shared<DevicePool>(s);
  }
  ~Impl() {
    try {
      if (staging)
        s.release(staging);
    } catch (...) {
    }
  }
  unsigned char *stage(std::size_t bytes) {
    if (bytes > staging_bytes) {
      if (staging)
        s.release(staging);
      staging_bytes = std::max<std::size_t>(align256(bytes), std::size_t(1) << 16);
      staging = s.alloc_pinned(staging_bytes, "native-v7.1-staging");
    }
    return static_cast<unsigned char *>(staging.pointer);
  }
  // Uploads bytes into a new pooled block (synchronous: the staging buffer is reused).
  std::shared_ptr<DeviceBlock> upload(const void *data, std::size_t bytes) {
    auto d = take(pool, bytes ? bytes : 8);
    if (bytes) {
      std::memcpy(stage(bytes), data, bytes);
      s.upload(view(*d), 0, staging, 0, bytes);
      s.synchronize("native v7.1 upload");
    }
    return d;
  }
  void download(void *out, const DeviceBlock &b, std::size_t offset, std::size_t bytes) {
    if (!bytes)
      return;
    stage(bytes);
    s.download(staging, 0, view(b), offset, bytes);
    s.synchronize("native v7.1 download");
    std::memcpy(out, staging.pointer, bytes);
  }
  void copy_words(u64 n, u64 src, u64 dst) {
    if (!n)
      return;
    void *a[] = {&n, &src, &dst};
    s.launch(copy, {grid_of(n)}, {block}, 0, a);
  }
  void fill_words(u64 n, u64 value, u64 dst) {
    if (!n)
      return;
    void *a[] = {&n, &value, &dst};
    s.launch(fill, {grid_of(n)}, {block}, 0, a);
  }
  // dst = m records of `words` words gathered from src by the device index array idx.
  void gather_words(u64 m, u64 idx, u64 n_src, u64 src, u64 words, u64 dst) {
    if (!m)
      return;
    void *a[] = {&m, &idx, &n_src, &src, &words, &dst};
    s.launch(gather, {grid_of(m * words)}, {block}, 0, a);
  }
};

FullTimeDevice::FullTimeDevice(hall_gpu::Session &s, pf_v7::gpu::Forms &f) : impl_(std::make_shared<Impl>(s, f)) {}
FullTimeDevice::~FullTimeDevice() = default;
FullTimeDeviceCounters FullTimeDevice::counters() const { return impl_->counters; }
std::uint64_t FullTimeDevice::row_materialisations() noexcept { return NativeRows::materialisations(); }
FullTimeDeviceCounters &device_counters(FullTimeDevice &d) { return d.impl_->counters; }

void device_read(FullTimeDevice &dev, std::uint64_t address, void *out, std::size_t bytes) {
  auto &I = *dev.impl_;
  if (!bytes)
    return;
  auto b = take(I.pool, (bytes + 7) & ~std::size_t(7));
  I.copy_words((bytes + 7) / 8, address, b->address);
  I.download(out, *b, 0, bytes);
  ++I.counters.row_downloads;
  I.counters.row_download_bytes += bytes;
}
std::vector<GapParticle> device_rows_download(const DeviceRows &r) {
  std::vector<GapParticle> out(r.n);
  auto &I = *r.dev;
  if (r.n) {
    const std::size_t bytes = r.n * sizeof(GapParticle);
    I.stage(bytes);
    I.s.download(I.staging, 0, view(*r.rows), 0, bytes);
    I.s.synchronize("native v7.1 row download");
    std::memcpy(out.data(), I.staging.pointer, bytes);
    ++I.counters.row_downloads;
    I.counters.row_download_bytes += bytes;
  }
  return out;
}

std::shared_ptr<const DeviceRows> device_upload_rows(FullTimeDevice &dev, const std::vector<GapParticle> &rows,
                                                     const std::vector<std::uint64_t> *keys,
                                                     const std::vector<std::uint64_t> *states) {
  auto &I = *dev.impl_;
  auto r = std::make_shared<DeviceRows>();
  r->dev = dev.impl_;
  r->n = rows.size();
  r->rows = I.upload(rows.data(), rows.size() * sizeof(GapParticle));
  ++I.counters.row_uploads;
  I.counters.row_upload_bytes += rows.size() * sizeof(GapParticle);
  if (keys && states) {
    auto c = std::make_shared<DeviceColumns>();
    c->n = rows.size();
    c->keys = I.upload(keys->data(), keys->size() * 8);
    c->states = I.upload(states->data(), states->size() * 8);
    I.counters.row_upload_bytes += 16 * rows.size();
    r->columns = c;
  }
  return r;
}

FullTimeCensusDevice device_census(FullTimeDevice &dev, const DeviceRows &rows, const std::string &states_sha256) {
  auto &I = *dev.impl_;
  if (!rows.columns)
    throw std::logic_error("native v7.1 device census without census columns");
  const u64 n = rows.n;
  auto rec = take(I.pool, std::max<std::size_t>(n * pf7::census_bytes, 8));
  I.f.census_pack(n, rows.rows->address, 10, rows.columns->keys->address, rows.columns->states->address,
                  rec->address);
  pf_v7::gpu::Forms::BlockListJob job{pf7::section_census, n, rec->address, 0};
  auto out = take(I.pool, 32);
  job.root_out = out->address;
  auto scratch = take(I.pool, pf_v7::gpu::Forms::block_lists_scratch_bytes({&job, 1}));
  I.f.block_lists({&job, 1}, scratch->address);
  unsigned char root[32];
  I.download(root, *out, 0, 32);
  FullTimeCensusDevice c;
  c.rows = rows.rows->address;
  c.packet_keys = rows.columns->keys->address;
  c.state_indices = rows.columns->states->address;
  c.carriers = n;
  c.states_sha256 = states_sha256;
  c.carriers_root = hex_words(root);
  c.sha256 = pf_v7::census_v2_sha256(n, states_sha256, c.carriers_root);
  return c;
}

DeviceMotionOutput device_admit_motion(FullTimeDevice &dev, const DeviceMotionInput &in) {
  auto &I = *dev.impl_;
  std::lock_guard<std::mutex> lock(I.m);
  const auto &r = *in.record;
  const auto &base = *in.base;
  if (!base.columns)
    throw std::logic_error("native v7.1 device admission without census columns");
  const u64 n = base.n;
  DeviceMotionOutput out;
  // ---- record arrays: in place when given on the device, uploaded otherwise.
  std::vector<std::shared_ptr<DeviceBlock>> held;
  auto array = [&](u64 device_address, const void *host, std::size_t bytes) -> u64 {
    if (device_address)
      return device_address;
    held.push_back(I.upload(host, bytes));
    ++I.counters.array_uploads;
    I.counters.array_upload_bytes += bytes;
    return held.back()->address;
  };
  const auto &D = r.device;
  const u64 n_e = D.endpoints ? D.n_endpoints : r.endpoints.size();
  const u64 n_p = D.paths ? D.n_paths : r.paths.size();
  const u64 n_E = D.E ? D.n_E : r.trial_E.size();
  const u64 n_B = D.B ? D.n_B : r.trial_B.size();
  const u64 n_J = D.J ? D.n_E : r.current_A_m.size();
  const u64 n_s = D.sinks ? D.n_sinks : in.sinks.size();
  const u64 n_b = D.births ? D.n_births : in.births.size();
  if (n_s > n || n_e != n - n_s || n_p != n || n_b > n)
    throw std::logic_error("native v8 device admission dimensions");
  const u64 ends = array(D.endpoints, r.endpoints.data(), n_e * sizeof(GapParticle));
  const u64 paths = array(D.paths, r.paths.data(), n * sizeof(FullTimeMotionPath));
  const u64 sinks = n_s ? array(D.sinks, in.sinks.data(), n_s * sizeof(pf_v7::Sink)) : 0;
  const u64 births = n_b ? array(D.births, in.births.data(), n_b * sizeof(pf_v7::Birth)) : 0;
  const bool extra = D.extra_J || !in.extra_J.empty();
  const u64 n_X = D.extra_J ? n_E : in.extra_J.size();
  const u64 Jx = extra ? array(D.extra_J, in.extra_J.data(), n_X * 8) : 0;
  out.extra_given = extra;
  const u64 E = array(D.E, r.trial_E.data(), n_E * 8);
  const u64 B = array(D.B, r.trial_B.data(), n_B * 8);
  const u64 J = array(D.J, r.current_A_m.data(), n_J * 8);
  pf_v7::HallAxialMesh mesh;
  out.mesh_ok = n_J == n_E && pf_v7::hall_axial_mesh_of(n_E, n_B, mesh);
  // ---- scratch
  const std::size_t nn = std::max<u64>(n, 1);
  auto terms = take(I.pool, 9 * nn * 8);
  const u64 k_w = terms->address, k_d = k_w + 8 * nn, g_w = k_d + 8 * nn, g_d = g_w + 8 * nn, c_w = g_d + 8 * nn,
            c_k = c_w + 16 * nn, q_a = c_k + 16 * nn;
  const unsigned g = grid_of(n);
  auto fails = take(I.pool, align256(g * 8) + align256(g * 4));
  const bool census_known = !in.census_root.empty();
  auto census_rec = take(I.pool, census_known ? 8 : nn * pf7::census_bytes);
  const unsigned sum_blocks = 80;
  auto sum_scratch = take(I.pool, pf_v7::gpu::Forms::sum_scratch_bytes(sum_blocks));
  auto after = take(I.pool, std::max<u64>(n_e, 1) * sizeof(GapParticle));
  // v8: sink/newborn indices (filled with ~0 so that a refused admission
  // gathers nothing out of range) and the compacted census columns.
  const bool flows = n_s || n_b;
  std::shared_ptr<DeviceBlock> idx, cols_keys, cols_states;
  if (flows) {
    idx = take(I.pool, 8 * (n_s + 2 * n_b + 1));
    I.fill_words(n_s + 2 * n_b, ~u64(0), idx->address);
  }
  if (n_s) {
    cols_keys = take(I.pool, 8 * std::max<u64>(n_e, 1));
    cols_states = take(I.pool, 8 * std::max<u64>(n_e, 1));
  }
  // results: fail[2] | 7 SumOut | E^ (16) | 9 roots | 2 id roots | v8 gathers
  constexpr std::size_t so = sizeof(pf_v7::gpu::SumOut);
  constexpr std::size_t off_sums = 16, off_bound = off_sums + 7 * so, off_roots = off_bound + 16,
                        off_ids = off_roots + 9 * 32, results_base = off_ids + 2 * 32;
  // v8 gathers: per sink row (10 words), path (8), key (1); per newborn row (10), endpoint (10), endpoint index (1)
  const std::size_t off_sunk = results_base, off_born = off_sunk + 8 * 19 * n_s,
                    results_bytes = off_born + 8 * 21 * n_b;
  auto results = take(I.pool, results_bytes);
  const u64 R = results->address;
  // ---- census records and per-carrier checks/terms
  if (!census_known)
    I.f.census_pack(n, base.rows->address, 10, base.columns->keys->address, base.columns->states->address,
                    census_rec->address);
  AdmitArgs a{};
  a.n = n;
  a.base = base.rows->address;
  a.ends = ends;
  a.paths = paths;
  a.keys = base.columns->keys->address;
  a.g.width = in.width;
  a.g.gap = in.gap;
  a.g.short_z = in.short_z;
  a.g.feed_z = in.feed_z;
  for (int k = 0; k < 4; ++k) {
    a.g.lower[k] = in.lower[k];
    a.g.upper[k] = in.upper[k];
  }
  for (int k = 0; k < 3; ++k)
    a.g.shift[k] = in.shift[k];
  a.g.dt = in.dt;
  a.g.neutral_bound = in.neutral_bound;
  a.g.finite_foils = in.finite_foils;
  a.g.finite_axial = in.finite_axial;
  a.k_w = k_w;
  a.k_d = k_d;
  a.g_w = g_w;
  a.g_d = g_d;
  a.c_w = c_w;
  a.c_k = c_k;
  a.q_a = q_a;
  a.block_fail = fails->address;
  a.block_code = fails->address + align256(g * 8);
  a.n_ends = n_e;
  a.sinks = sinks;
  a.n_sinks = n_s;
  a.births = births;
  a.n_births = n_b;
  a.states = base.columns->states->address;
  a.after_keys = n_s ? cols_keys->address : 0;
  a.after_states = n_s ? cols_states->address : 0;
  a.sunk_index = flows ? idx->address : 0;
  a.born_index = flows ? idx->address + 8 * n_s : 0;
  a.born_end = flows ? idx->address + 8 * (n_s + n_b) : 0;
  if (n) {
    void *args[] = {&a};
    I.s.launch(I.admit, {g}, {block}, 0, args);
  }
  {
    u64 blocks = n ? g : 0, bf = a.block_fail, bc = a.block_code, o = R;
    void *args[] = {&blocks, &bf, &bc, &o};
    I.s.launch(I.fail_finish, {1}, {32}, 0, args);
  }
  // ---- exact sums (one scratch, stream-ordered)
  const u64 S = sum_scratch->address;
  I.f.sum_products(n, k_w, 1, k_d, 1, 1.0, sum_blocks, S, R + off_sums + 0 * so);
  I.f.sum_products(n, g_w, 1, g_d, 1, 1.0, sum_blocks, S, R + off_sums + 1 * so);
  I.f.sum_products(2 * n, c_w, 1, c_k, 1, 1.0, sum_blocks, S, R + off_sums + 2 * so);
  I.f.sum_products(n, q_a, 1, k_w, 1, 1.0, sum_blocks, S, R + off_sums + 3 * so);
  const bool field = out.mesh_ok;
  if (field) {
    I.f.sum_products(n_E, E, 1, J, 1, in.dt, sum_blocks, S, R + off_sums + 4 * so);
    I.f.sum_abs_products(n_E, E, 1, J, 1, in.dt, sum_blocks, S, R + off_sums + 5 * so);
    auto fb = take(I.pool, pf_v7::gpu::Forms::field_bound_scratch_bytes(n_E, mesh.nn()));
    I.f.field_bound(n_E, E, J, mesh.n, mesh.nz, fb->address, R + off_bound);
    held.push_back(fb);
    if (extra && n_X == n_E)
      I.f.sum_products(n_E, E, 1, Jx, 1, in.dt, sum_blocks, S, R + off_sums + 6 * so);
  }
  // ---- fused block lists: endpoints, paths, sinks, births, census, E, B, J
  std::vector<pf_v7::gpu::Forms::BlockListJob> jobs = {
      {pf7::section_endpoints, n_e, ends, R + off_roots + 0 * 32},
      {pf7::section_paths, n, paths, R + off_roots + 1 * 32},
      {pf7::section_sinks, n_s, n_s ? sinks : ends, R + off_roots + 2 * 32},
      {pf7::section_births, n_b, n_b ? births : ends, R + off_roots + 3 * 32},
      {pf7::section_census, census_known ? 0 : n, census_rec->address, R + off_roots + 4 * 32},
      {pf7::section_field_E, n_E, E, R + off_roots + 5 * 32},
      {pf7::section_field_B, n_B, B, R + off_roots + 6 * 32},
      {pf7::section_current_J, n_J, J, R + off_roots + 7 * 32}};
  auto bl = take(I.pool, pf_v7::gpu::Forms::block_lists_scratch_bytes(jobs));
  I.f.block_lists(jobs, bl->address);
  std::shared_ptr<DeviceBlock> blx;
  if (extra) { // v8.1 extra current: its own block list (the fused launch takes at most 8 sections)
    const pf_v7::gpu::Forms::BlockListJob xj{pf7::section_current_J, n_X, Jx, R + off_roots + 8 * 32};
    blx = take(I.pool, pf_v7::gpu::Forms::block_lists_scratch_bytes({&xj, 1}));
    I.f.block_lists({&xj, 1}, blx->address);
  }
  // ---- after rows (copy of the endpoints) and the v5 ID-block roots
  I.copy_words(10 * n_e, ends, after->address);
  std::shared_ptr<DeviceBlock> ids, ids2;
  if (!n_s) {
    // Both v5 roots in one leaf launch (after rows and paths carry the same IDs whenever the per-carrier checks
    // pass; otherwise the admission is refused and the roots are not used).
    ids = take(I.pool, pf_v7::gpu::Forms::id_root_pair_scratch_bytes(n));
    I.f.id_root_pair(n, after->address, 10, paths, 8, ids->address, R + off_ids, R + off_ids + 32);
  } else {
    ids = take(I.pool, pf_v7::gpu::Forms::id_root_scratch_bytes(n_e));
    I.f.id_root(0, n_e, after->address, 10, ids->address, R + off_ids);
    ids2 = take(I.pool, pf_v7::gpu::Forms::id_root_scratch_bytes(n));
    I.f.id_root(1, n, paths, 8, ids2->address, R + off_ids + 32);
  }
  // ---- v8 gathers for the host's sink and newborn checks
  if (n_s) {
    const u64 si = idx->address;
    I.gather_words(n_s, si, n, base.rows->address, 10, R + off_sunk);
    I.gather_words(n_s, si, n, paths, 8, R + off_sunk + 80 * n_s);
    I.gather_words(n_s, si, n, base.columns->keys->address, 1, R + off_sunk + 144 * n_s);
  }
  if (n_b) {
    const u64 bi = idx->address + 8 * n_s, be = bi + 8 * n_b;
    I.gather_words(n_b, bi, n, base.rows->address, 10, R + off_born);
    I.gather_words(n_b, be, n_e, ends, 10, R + off_born + 80 * n_b);
    I.copy_words(n_b, be, R + off_born + 160 * n_b);
  }
  // ---- one download
  std::vector<unsigned char> res(results_bytes);
  I.download(res.data(), *results, 0, results_bytes);
  u64 fail[2];
  std::memcpy(fail, res.data(), 16);
  out.first_fail = fail[0];
  out.fail_code = std::uint32_t(fail[1]);
  out.carriers_pass = fail[0] == ~u64(0);
  pf_v7::gpu::SumOut sums[7];
  std::memcpy(sums, res.data() + off_sums, sizeof sums);
  out.kinetic.value = out.kinetic.sum = sums[0].value;
  out.kinetic.status = sums[0].status;
  out.scales.kinetic_gross_J = sums[1].value;
  out.scales.kinetic_content_J = sums[2].value;
  out.scales.charge_abs_C = sums[3].value;
  out.scales.status = sums[1].status | sums[2].status | sums[3].status;
  if (field) {
    out.work.sum = sums[4].value;
    out.work.value = sums[4].scaled;
    out.work.status = sums[4].status;
    out.scales.field_work_gross_J = sums[5].scaled;
    out.scales.status |= sums[5].status;
    double em;
    std::uint32_t bad;
    std::memcpy(&em, res.data() + off_bound, 8);
    std::memcpy(&bad, res.data() + off_bound + 8, 4);
    out.scales.field_bound_V_m = em;
    if (bad)
      out.scales.status |= pf7::sum_nonfinite_input;
    out.scales.tolerance_J = pf_v7::work_tolerance(out.scales);
    if (extra && n_X == n_E) {
      out.extra_work.sum = sums[6].value;
      out.extra_work.value = sums[6].scaled;
      out.extra_work.status = sums[6].status;
    } else if (extra)
      out.extra_work.status = pf7::sum_nonfinite_input; // not one Hall edge array
  }
  std::string roots[9];
  for (int k = 0; k < 9; ++k)
    roots[k] = hex_words(res.data() + off_roots + 32 * k);
  out.sections = {n_e, n, n_s, n_b, roots[0], roots[1], roots[2], roots[3]};
  if (n_s) {
    out.sunk_rows.resize(n_s);
    out.sunk_paths.resize(n_s);
    out.sunk_keys.resize(n_s);
    std::memcpy(out.sunk_rows.data(), res.data() + off_sunk, 80 * n_s);
    std::memcpy(out.sunk_paths.data(), res.data() + off_sunk + 80 * n_s, 64 * n_s);
    std::memcpy(out.sunk_keys.data(), res.data() + off_sunk + 144 * n_s, 8 * n_s);
  }
  if (n_b) {
    out.born_rows.resize(n_b);
    out.born_ends.resize(n_b);
    out.born_end_index.resize(n_b);
    std::memcpy(out.born_rows.data(), res.data() + off_born, 80 * n_b);
    std::memcpy(out.born_ends.data(), res.data() + off_born + 80 * n_b, 80 * n_b);
    std::memcpy(out.born_end_index.data(), res.data() + off_born + 160 * n_b, 8 * n_b);
  }
  out.census_root = census_known ? in.census_root : roots[4];
  out.census_sha256 = pf_v7::census_v2_sha256(n, in.states_sha256, out.census_root);
  out.E_root = roots[5];
  out.B_root = roots[6];
  out.J_root = roots[7];
  out.after_carrier_root = hex_words(res.data() + off_ids);
  out.paths_v5_root = hex_words(res.data() + off_ids + 32);
  // ---- field/current digests: block form, otherwise the v7 streaming form.
  if (field) {
    out.field_digest_ok = r.trial_field_sha256 == pf_v7::trial_field_blocks_sha256(out.E_root, out.B_root);
    out.current_digest_ok = r.current_A_m_sha256 == pf_v7::current_blocks_sha256(out.J_root);
    if (!out.field_digest_ok || !out.current_digest_ok) {
      // v7 streaming forms (one sequential SHA-256 stream each): hashed on the host from the host spans, or from
      // the device arrays read back (O(mesh) bytes). One SHA-256 stream is much faster on the host than on one
      // device thread and does not occupy the shared stream.
      std::vector<double> tmp_E, tmp_B, tmp_J;
      auto host_array = [&](u64 device_address, std::span<const double> host, u64 count,
                            std::vector<double> &tmp) -> std::span<const double> {
        if (!device_address)
          return host;
        tmp.resize(count);
        if (count) {
          auto b = take(I.pool, count * 8);
          I.copy_words(count, device_address, b->address);
          I.download(tmp.data(), *b, 0, count * 8);
          ++I.counters.array_downloads;
          I.counters.array_download_bytes += count * 8;
        }
        return tmp;
      };
      if (!out.field_digest_ok) {
        out.field_digest_ok = r.trial_field_sha256 == pf_v7::total_field_sha256(host_array(D.E, r.trial_E, n_E, tmp_E),
                                                                                 host_array(D.B, r.trial_B, n_B, tmp_B));
        ++I.counters.streaming_field_digests;
      }
      if (!out.current_digest_ok) {
        out.current_digest_ok = r.current_A_m_sha256 == pf_v7::array_sha256(host_array(D.J, r.current_A_m, n_J, tmp_J));
        ++I.counters.streaming_field_digests;
      }
    }
    if (extra && n_X == n_E) {
      out.extra_digest_ok = r.extra_current_A_m_sha256 == pf_v7::current_blocks_sha256(roots[8]);
      if (!out.extra_digest_ok) {
        std::vector<double> tmp_X;
        if (D.extra_J) {
          tmp_X.resize(n_X);
          auto b = take(I.pool, n_X * 8);
          I.copy_words(n_X, D.extra_J, b->address);
          I.download(tmp_X.data(), *b, 0, n_X * 8);
        }
        out.extra_digest_ok = r.extra_current_A_m_sha256 ==
                              pf_v7::array_sha256(D.extra_J ? std::span<const double>(tmp_X) : in.extra_J);
        ++I.counters.streaming_field_digests;
      }
    }
  }
  auto rows = std::make_shared<DeviceRows>();
  rows->dev = dev.impl_;
  rows->rows = after;
  if (n_s) {
    auto c = std::make_shared<DeviceColumns>();
    c->keys = cols_keys;
    c->states = cols_states;
    c->n = n_e;
    rows->columns = c;
  } else
    rows->columns = base.columns;
  rows->n = n_e;
  rows->carrier_root = out.after_carrier_root;
  out.after = rows;
  ++I.counters.device_admissions;
  return out;
}

// ---------------------------------------------------------------- v8 chemistry device stage
namespace {
double full_time_detail_now_ms() {
  timespec t;
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
  return 1e3 * double(t.tv_sec) + 1e-6 * double(t.tv_nsec);
}
struct ChemModule {
  hall_gpu::Module module;
  std::map<std::string, hall_gpu::Kernel> k;
  std::string blob_sha;
  std::shared_ptr<DeviceBlock> blob, vacuum, proj, targ;
  pfc::TableView tv{};
  pfc::CoverView cv{};
};
std::mutex chem_modules_m;
// The chemistry module lives in the device lane (FullTimeDevice::Impl::chem):
// released with it, before its session closes.
std::shared_ptr<ChemModule> chem_slot(FullTimeDevice::Impl &I) { return std::static_pointer_cast<ChemModule>(I.chem); }
// The compiled program (no table data): row edits, compaction and scans are
// also used without chemistry (v8 device births).
std::shared_ptr<ChemModule> chem_program_locked(FullTimeDevice::Impl &I) {
  auto m = chem_slot(I);
  if (!m) {
    m = std::make_shared<ChemModule>();
    I.chem = m;
    hall_gpu::CompileRequest req;
    req.main = {"chem_v3_kernels.cu", kChem};
    req.headers = {{"forms_core.h", kCore},      {"chem_event.h", kChemEvent}, {"chem_kinematics.h", kChemKin},
                   {"chem_rng.h", kChemRng},     {"chem_tables.h", kChemTables}, {"det_math.h", kChemMath},
                   {"philox_core.h", kPhilox}};
    m->module = I.s.compile_and_load(req);
    for (const char *name : {"ch_zero_words", "ch_scan_blocks", "ch_scan_totals", "ch_scan_add", "ch_fill_u32", "ch_classify", "ch_multi",
                             "ch_mark", "ch_cells", "ch_carriers", "ch_seg_count", "ch_seg_fill", "ch_seg_sum",
                             "ch_cell_count", "ch_mcc", "ch_summary", "ch_leaves", "ch_events_copy", "ch_add_total",
                             "ch_active", "ch_gather_inv", "ch_gather_seg", "ch_nonzero", "ch_record_leaves",
                             "ch_widen", "ch_edit_velocity", "ch_edit_weight", "ch_drop", "ch_compact"})
      m->k[name] = I.s.kernel(m->module, name);
  }
  return m;
}
ChemModule &chem_module(FullTimeDevice::Impl &I, const chem8::Setup &st) {
  std::lock_guard<std::mutex> lock(chem_modules_m);
  auto m = chem_program_locked(I);
  const auto &b = *st.reg->blob;
  if (m->blob_sha != b.sha256) {
    m->blob = I.upload(b.base(), b.bytes);
    m->tv = b.view_at(m->blob->address);
    m->cv = b.cover_at(m->blob->address);
    std::vector<std::uint32_t> vac(st.grid.vacuum.begin(), st.grid.vacuum.end()), pj(st.reg->projectile.begin(),
                                                                                       st.reg->projectile.end()),
        tg(st.reg->target.begin(), st.reg->target.end());
    m->vacuum = I.upload(vac.data(), vac.size() * 4);
    m->proj = I.upload(pj.data(), pj.size() * 4);
    m->targ = I.upload(tg.data(), tg.size() * 4);
    m->blob_sha = b.sha256;
  }
  return *m;
}
// Device-side state of one stage kept for the callbacks.
struct ChemKeep {
  std::shared_ptr<FullTimeDevice::Impl> dev;
  std::shared_ptr<ChemModule> mod;
  std::vector<std::shared_ptr<DeviceBlock>> blocks;
  u64 rows = 0, keys = 0, flags = 0, seg_rows = 0, inv = 0;
  u64 seg_count = 0, seg_off = 0, seg_rank = 0, n_seg = 0; // device arrays (compact cell x 22)
};
struct ClassArgs {
  u64 n, rows, states, keys, state_species, n_states, proj, targ;
  double low[3], high[3];
  std::uint32_t cnt[3], pad;
  u64 vacuum, species, cell, flags, pflag, occupied, err;
};
struct SegArgs {
  u64 n_seg, rows, seg_count, seg_rows, seg_off, seg_rank;
  double volume;
  u64 cells, inv, err;
};
struct SumArgs {
  std::uint32_t begin, count;
  u64 outs, nev, active, bad, candidates, status, numax, region, support;
};
struct ActiveRec {
  u64 id;
  std::uint32_t index, row, cell, species, nev, status;
  double x[3];
  double weight;
};
static_assert(sizeof(ActiveRec) == 64);
std::array<unsigned char, 32> digest_words(const std::uint32_t *w) {
  const auto d = pf_v7::digest_of_words(w);
  std::array<unsigned char, 32> out;
  std::memcpy(out.data(), d.data(), 32);
  return out;
}
chem8::ESum from_digits(const std::int64_t *d) {
  // value = sum_k d_k 2^(32k - 2148): digits 0..132 in [0, 2^32), 133 signed.
  chem8::ESum e;
  for (int k = 0; k < pf7::acc_digits; ++k) {
    if (!d[k])
      continue;
    const int ex = 32 * k - pf7::acc_offset;
    if (ex > 960 || ex < -1092)
      throw std::logic_error("native v8 chemistry exact digits outside the ledger range");
    e.add(std::ldexp(double(d[k]), ex));
  }
  return e;
}
} // namespace
namespace {
struct Scan {
  FullTimeDevice::Impl &I;
  ChemModule &M;
  ChemKeep &K;
  // exclusive scan of n u32 -> u64 (block) ; grand total at *grand
  void operator()(u64 in, u64 n, u64 out, u64 grand) {
    const u64 blocks = n ? (n + 255) / 256 : 1;
    auto totals = take(I.pool, blocks * 8 + 8);
    K.blocks.push_back(totals);
    u64 t = totals->address;
    if (n) {
      void *a[] = {&in, &n, &out, &t};
      I.s.launch(M.k["ch_scan_blocks"], {unsigned(blocks)}, {256}, 0, a);
    } else {
      u64 w = 1;
      void *z[] = {&t, &w};
      I.s.launch(M.k["ch_zero_words"], {1}, {256}, 0, z);
    }
    u64 nb = n ? blocks : 0;
    void *b[] = {&t, &nb, &grand};
    I.s.launch(M.k["ch_scan_totals"], {1}, {1}, 0, b);
    if (n) {
      void *c[] = {&out, &n, &t};
      I.s.launch(M.k["ch_scan_add"], {unsigned(blocks)}, {256}, 0, c);
    }
  }
};
} // namespace

ChemDeviceStage device_chem_stage(FullTimeDevice &dev, const ChemDeviceRequest &rq) {
  auto &I = *dev.impl_;
  std::lock_guard<std::mutex> lock(I.m);
  const auto t0 = std::chrono::steady_clock::now();
  const auto &st = *rq.setup;
  auto &M = chem_module(I, st);
  auto keep = std::make_shared<ChemKeep>();
  auto &K = *keep;
  K.dev = dev.impl_;
  {
    std::lock_guard<std::mutex> l2(chem_modules_m);
    K.mod = chem_slot(I);
  }
  auto blk = [&](std::size_t bytes) {
    auto b = take(I.pool, std::max<std::size_t>(bytes, 8));
    K.blocks.push_back(b);
    return b->address;
  };
  auto zero = [&](u64 addr, std::size_t bytes) {
    u64 a = addr, w = (bytes + 7) / 8;
    void *args[] = {&a, &w};
    if (w)
      I.s.launch(M.k["ch_zero_words"], {unsigned((w + 255) / 256)}, {256}, 0, args);
  };
  auto launch = [&](const char *name, u64 n, void **args, unsigned bs = 256) {
    if (!n)
      return;
    I.s.launch(M.k[name], {unsigned((n + bs - 1) / bs)}, {bs}, 0, args);
  };
  Scan scan{I, M, K};
  ChemDeviceStage out;
  const auto &rows = *rq.rows;
  if (!rows.columns)
    throw std::logic_error("native v8 chemistry device stage without census columns");
  const u64 n = rows.n;
  const auto &g = st.grid;
  const u64 n_cells = g.cells();
  const std::uint32_t S = 22, R = st.reg->blob->n_regions;
  // ---- classification, compaction, targets
  ClassArgs a{};
  a.n = n;
  a.rows = rows.rows->address;
  a.states = rows.columns->states->address;
  a.keys = rows.columns->keys->address;
  a.n_states = rq.state_species.size();
  a.state_species = blk(a.n_states * 4);
  {
    auto up = I.upload(rq.state_species.data(), rq.state_species.size() * 4);
    K.blocks.push_back(up);
    a.state_species = up->address;
  }
  a.proj = M.proj->address;
  a.targ = M.targ->address;
  for (int d = 0; d < 3; ++d) {
    a.low[d] = g.low[d];
    a.high[d] = g.high[d];
    a.cnt[d] = g.n[d];
  }
  a.vacuum = M.vacuum->address;
  a.species = blk(n * 4);
  a.cell = blk(n * 4);
  a.flags = blk(n * 4);
  a.pflag = blk(n * 4);
  a.occupied = blk(n_cells * 4);
  zero(a.occupied, n_cells * 4);
  a.err = blk(4 * 8);
  {
    const u64 init[4] = {~u64(0), ~u64(0), ~u64(0), ~u64(0)};
    auto up = I.upload(init, 32);
    K.blocks.push_back(up);
    I.copy_words(4, up->address, a.err);
  }
  {
    void *args[] = {&a};
    launch("ch_classify", n, args);
    launch("ch_multi", n, args);
    launch("ch_mark", n, args);
  }
  const u64 crank = blk(n_cells * 8), prank = blk(n * 8), grands = blk(8 * 8);
  scan(a.occupied, n_cells, crank, grands + 0);
  scan(a.pflag, n, prank, grands + 8);
  u64 head[6];
  {
    I.copy_words(4, a.err, grands + 16);
    std::vector<unsigned char> b(48);
    auto tmp = take(I.pool, 64);
    I.copy_words(6, grands, tmp->address);
    I.download(head, *tmp, 0, 48);
  }
  const u64 n_occ = head[0], n_proj = head[1];
  for (int c = 0; c < 3; ++c)
    if (head[2 + c] != ~u64(0)) {
      out.error = c + 1;
      out.error_index = head[2 + c];
      return out;
    }
  if (n_proj > 0xffffffffull)
    throw std::length_error("native chemistry projectile index exceeds packed refusal domain");
  out.projectiles = n_proj;
  const u64 cell_atlas = blk(n_occ * 4);
  {
    u64 occ = a.occupied, cr = crank, nc = n_cells, ca = cell_atlas;
    void *args[] = {&occ, &cr, &nc, &ca};
    launch("ch_cells", n_cells, args);
  }
  const u64 carriers = blk(n_proj * sizeof(pfc::CarrierIn)), carrier_row = blk(n_proj * 4);
  {
    u64 pr = prank, cr = crank, c = carriers, row = carrier_row;
    void *args[] = {&a, &pr, &cr, &c, &row};
    launch("ch_carriers", n, args);
  }
  const u64 n_seg = n_occ * S;
  const u64 seg_count = blk(n_seg * 4), seg_off = blk(n_seg * 8), seg_fill = blk(n_seg * 4), present = blk(n_seg * 4),
            seg_rank = blk(n_seg * 8);
  zero(seg_count, n_seg * 4);
  zero(seg_fill, n_seg * 4);
  {
    u64 cr = crank, sc = seg_count;
    void *args[] = {&a, &cr, &sc};
    launch("ch_seg_count", n, args);
  }
  scan(seg_count, n_seg, seg_off, grands + 40);
  {
    u64 sc = seg_count, ns = n_seg, pr = present;
    void *args[] = {&sc, &ns, &pr};
    launch("ch_nonzero", n_seg, args);
  }
  scan(present, n_seg, seg_rank, grands + 48);
  u64 tot[2];
  {
    auto tmp = take(I.pool, 16);
    I.copy_words(2, grands + 40, tmp->address);
    I.download(tot, *tmp, 0, 16);
  }
  const u64 T = tot[0], n_present = tot[1];
  const u64 seg_rows = blk(T * 4);
  {
    u64 cr = crank, so = seg_off, sf = seg_fill, sr = seg_rows;
    void *args[] = {&a, &cr, &so, &sf, &sr};
    launch("ch_seg_fill", n, args);
  }
  const u64 cells = blk(n_occ * sizeof(pfc::CellTargets)), inv = blk(n_present * pf7::acc_digits * 8);
  zero(cells, n_occ * sizeof(pfc::CellTargets));
  {
    SegArgs sa{n_seg, rows.rows->address, seg_count, seg_rows, seg_off, seg_rank, g.volume, cells, inv, a.err};
    void *args[] = {&sa};
    launch("ch_seg_sum", n_seg, args, 64);
    u64 sc = seg_count, no = n_occ, ce = cells;
    void *b2[] = {&sc, &no, &ce};
    launch("ch_cell_count", n_occ, b2);
  }
  {
    u64 e4;
    auto tmp = take(I.pool, 8);
    I.copy_words(1, a.err + 24, tmp->address);
    I.download(&e4, *tmp, 0, 8);
    if (e4 != ~u64(0)) {
      out.error = 4;
      out.error_index = e4;
      return out;
    }
  }
  K.rows = rows.rows->address;
  K.keys = rows.columns->keys->address;
  K.flags = a.flags;
  K.seg_rows = seg_rows;
  K.inv = inv;
  K.seg_count = seg_count;
  K.seg_off = seg_off;
  K.seg_rank = seg_rank;
  K.n_seg = n_seg;
  // ---- input leaves (carriers 1024/leaf, cells 64/leaf, atlas 4096/leaf)
  auto leaves_of = [&](const char *tag, u64 data, u64 words, u64 count, u64 per_leaf) {
    const u64 nl = (count + per_leaf - 1) / per_leaf;
    const std::string t = std::string(tag) + "/leaf";
    auto tg = I.upload(t.data(), t.size());
    K.blocks.push_back(tg);
    const u64 out_l = blk(nl * 32);
    u64 d = data, w = words, c = count, pl = per_leaf, ta = tg->address, o = out_l;
    unsigned tl = unsigned(t.size());
    void *args[] = {&c, &w, &pl, &d, &ta, &tl, &o};
    if (nl)
      I.s.launch(M.k["ch_record_leaves"], {unsigned((nl + 63) / 64)}, {64}, 0, args);
    return std::pair<u64, u64>{out_l, nl};
  };
  const u64 atlas64 = blk(n_occ * 8);
  {
    u64 no = n_occ, ca = cell_atlas, a64 = atlas64;
    void *args[] = {&no, &ca, &a64};
    launch("ch_widen", n_occ, args);
  }
  const auto Lc = leaves_of("native-v8/chemistry-carriers", carriers, 6, n_proj, 1024);
  const auto Lt = leaves_of("native-v8/chemistry-cells", cells, sizeof(pfc::CellTargets) / 8, n_occ, 64);
  const auto La = leaves_of("native-v8/chemistry-atlas", atlas64, 1, n_occ, 4096);
  // ---- the C3 stage in chunks
  const auto &cfg = rq.cfg;
  const u64 maxe = cfg.max_events;
  u64 chunk = (std::size_t(512) << 20) / (maxe * sizeof(pfc::EventRec) + sizeof(pfc::CarrierOut));
  chunk = std::max<u64>(256, std::min<u64>(chunk / 256 * 256, 1u << 20));
  if (rq.partitions > 0) // explicit partition count (qualification): ceil(n / P) rounded up to 256
    chunk = std::min<u64>(chunk, std::max<u64>(256, ((n_proj + rq.partitions - 1) / rq.partitions + 255) / 256 * 256));
  const u64 ev_chunk = blk(chunk * maxe * sizeof(pfc::EventRec)), out_chunk = blk(chunk * sizeof(pfc::CarrierOut));
  const u64 nev = blk(n_proj * 4), activef = blk(n_proj * 4), offs = blk(chunk * 8);
  const u64 counters = blk(16 * 8); // historical scalar offsets; status counters now separately allocated
  const u64 statuses = blk(9 * 8); // explicit OK..CDF_SUPPORT, never mask status8 to OK
  const u64 region = blk(u64(R) * 2 * pf7::acc_digits * 8);
  const u64 support = blk(u64(R) * 8);
  const u64 n_leaves = (n_proj + 255) / 256;
  const u64 rec_leaves = blk(n_leaves * 32);
  const std::string rtag = "native-v8/chemistry-records-leaf";
  auto rtag_dev = I.upload(rtag.data(), rtag.size());
  K.blocks.push_back(rtag_dev);
  u64 cap = std::max<u64>(1u << 16, n_proj / 4);
  u64 compact = 0;
  u64 ctr[16], status_counts[9];
  for (int attempt = 0; attempt < 2; ++attempt) {
    compact = blk(cap * sizeof(pfc::EventRec));
    zero(counters, 16 * 8);
    zero(statuses, 9 * 8);
    {
      const u64 bad = ~u64(0);
      auto up = I.upload(&bad, 8);
      K.blocks.push_back(up);
      I.copy_words(1, up->address, counters);
    }
    zero(region, u64(R) * 2 * pf7::acc_digits * 8);
    zero(support, u64(R) * 8);
    for (u64 b0 = 0; b0 < n_proj; b0 += chunk) {
      const unsigned begin = unsigned(b0), count = unsigned(std::min<u64>(chunk, n_proj - b0));
      {
        pfc::TableView tv = M.tv;
        pfc::CoverView cv = M.cv;
        pfc::StepConfig c = cfg;
        u64 ca = carriers, ce = cells, ev = ev_chunk, oc = out_chunk;
        unsigned nc = unsigned(n_occ), bg = begin, ct = count;
        void *args[] = {&tv, &cv, &c, &ca, &ce, &nc, &bg, &ct, &ev, &oc};
        I.s.launch(M.k["ch_mcc"], {(count + 127) / 128}, {128}, 0, args);
      }
      {
        SumArgs sa{begin, count, out_chunk, nev, activef, counters, counters + 8, statuses, counters + 80, region, support};
        void *args[] = {&sa};
        I.s.launch(M.k["ch_summary"], {(count + 255) / 256}, {256}, 0, args);
      }
      {
        unsigned bg = begin, ct = count, me = unsigned(maxe), tl = unsigned(rtag.size());
        u64 oc = out_chunk, ev = ev_chunk, ta = rtag_dev->address, lv = rec_leaves;
        void *args[] = {&bg, &ct, &me, &oc, &ev, &ta, &tl, &lv};
        const unsigned nl = (count + 255) / 256;
        I.s.launch(M.k["ch_leaves"], {(nl + 31) / 32}, {32}, 0, args);
      }
      scan(nev + 4 * u64(begin), count, offs, counters + 104);
      {
        unsigned bg = begin, ct = count, me = unsigned(maxe);
        u64 nv = nev, of = offs, base = counters + 88, cp = cap, ev = ev_chunk, co = compact, ov = counters + 96;
        void *args[] = {&bg, &ct, &me, &nv, &of, &base, &cp, &ev, &co, &ov};
        I.s.launch(M.k["ch_events_copy"], {(count + 255) / 256}, {256}, 0, args);
        u64 b2 = counters + 88, gr = counters + 104;
        void *args2[] = {&b2, &gr};
        I.s.launch(M.k["ch_add_total"], {1}, {1}, 0, args2);
      }
    }
    // active carriers: rank
    const u64 arank = blk(n_proj * 8);
    scan(activef, n_proj, arank, counters + 112);
    {
      auto tmp = take(I.pool, 128);
      I.copy_words(16, counters, tmp->address);
      I.download(ctr, *tmp, 0, 128);
    }
    if (ctr[12]) { // overflow: rerun with the exact capacity
      cap = ctr[11] + 1;
      ++out.retries;
      continue;
    }
    const u64 n_act = ctr[14];
    const u64 act = blk(n_act * sizeof(ActiveRec));
    {
      u64 np = n_proj, af = activef, ar = arank, nv = nev, ca = carriers, cr = carrier_row, rw = rows.rows->address,
          o = act;
      void *args[] = {&np, &af, &ar, &nv, &ca, &cr, &rw, &o};
      launch("ch_active", n_proj, args);
    }
    // ---- downloads
    const u64 n_ev = ctr[11];
    std::vector<ActiveRec> ar(n_act);
    out.active.events.resize(n_ev);
    std::vector<std::int64_t> reg(u64(R) * 2 * pf7::acc_digits);
    out.active.totals.support.resize(R);
    std::vector<std::uint32_t> lw(n_leaves * 8), lc(Lc.second * 8), lt(Lt.second * 8), la(La.second * 8);
    out.active.cell_atlas.resize(n_occ);
    auto dl = [&](void *dst, u64 addr, std::size_t bytes) {
      if (!bytes)
        return;
      auto b = take(I.pool, (bytes + 7) & ~std::size_t(7));
      I.copy_words((bytes + 7) / 8, addr, b->address);
      I.download(dst, *b, 0, bytes);
    };
    dl(status_counts, statuses, sizeof status_counts);
    dl(ar.data(), act, n_act * sizeof(ActiveRec));
    dl(out.active.events.data(), compact, n_ev * sizeof(pfc::EventRec));
    dl(reg.data(), region, reg.size() * 8);
    dl(out.active.totals.support.data(), support, u64(R) * 8);
    dl(lw.data(), rec_leaves, lw.size() * 4);
    dl(lc.data(), Lc.first, lc.size() * 4);
    dl(lt.data(), Lt.first, lt.size() * 4);
    dl(la.data(), La.first, la.size() * 4);
    dl(out.active.cell_atlas.data(), cell_atlas, n_occ * 4);
    I.counters.row_download_bytes += n_ev * sizeof(pfc::EventRec);
    const double th0 = full_time_detail_now_ms();
    // ---- host assembly
    auto &A = out.active;
    u64 at = 0;
    for (const auto &r : ar) {
      chem8::ActiveCarrier c;
      c.row = r.row;
      c.cell = r.cell;
      c.species = r.species;
      c.id = r.id;
      c.begin = at;
      at += r.nev;
      c.end = at;
      for (int d = 0; d < 3; ++d)
        c.x[d] = r.x[d];
      c.weight = r.weight;
      A.carriers.push_back(c);
    }
    if (at != n_ev)
      throw std::logic_error("native v8 chemistry device event count");
    auto &T2 = A.totals;
    T2.carriers = n_proj;
    T2.candidates = ctr[1];
    for (std::uint32_t k = 0; k < 9; ++k)
      if (status_counts[k])
        T2.status_counts[k] = status_counts[k];
    T2.nu_dt_max = std::bit_cast<double>(ctr[10]);
    T2.exposure.resize(R);
    T2.bound.resize(R);
    for (std::uint32_t r = 0; r < R; ++r) {
      std::array<std::int64_t, pf7::acc_digits> dx, db;
      std::memcpy(dx.data(), reg.data() + (u64(r) * 2) * pf7::acc_digits, sizeof dx);
      std::memcpy(db.data(), reg.data() + (u64(r) * 2 + 1) * pf7::acc_digits, sizeof db);
      // normalise (carries), then convert
      auto norm = [](std::array<std::int64_t, pf7::acc_digits> &d) {
        for (int k = 0; k + 1 < pf7::acc_digits; ++k) {
          const std::int64_t carry = d[k] >> 32;
          d[k] &= 0xffffffffll;
          d[k + 1] += carry;
        }
      };
      norm(dx);
      norm(db);
      T2.exposure[r] = from_digits(dx.data());
      T2.bound[r] = from_digits(db.data());
    }
    if (ctr[0] != ~u64(0)) {
      pfc::CarrierIn c;
      // ch_summary atomically minimizes (uint32 projectile index, uint32 status).
      // The status belongs to that exact first bad carrier, even with mixed errors.
      const u64 index = ctr[0] >> 32;
      const auto status = std::uint32_t(ctr[0]);
      if (index >= n_proj || status < PFC_ST_CAND_CAP || status > PFC_ST_CDF_SUPPORT)
        throw std::logic_error("native chemistry invalid packed refusal");
      dl(&c, carriers + index * sizeof(pfc::CarrierIn), sizeof c);
      T2.bad_carrier = c.id;
      T2.bad_status = status;
    }
    auto roots = [&](const std::string &tag, u64 count, const std::vector<std::uint32_t> &w) {
      std::vector<std::array<unsigned char, 32>> l(w.size() / 8);
      for (std::size_t k = 0; k < l.size(); ++k)
        l[k] = digest_words(w.data() + 8 * k);
      return chem8::block_root_of(tag, count, l);
    };
    {
      std::vector<std::array<unsigned char, 32>> l(n_leaves);
      for (u64 k = 0; k < n_leaves; ++k)
        l[k] = digest_words(lw.data() + 8 * k);
      pfc::Sha256 h;
      h.text("native-v8/chemistry-records-root");
      h.pod(n_proj);
      h.pod(n_leaves);
      for (const auto &x : l)
        h.update(x.data(), 32);
      out.records_root = h.hex();
    }
    if (std::getenv("PF_V8_PROFILE"))
      std::fprintf(stderr, "chem-profile device_host_assembly %.3f\n", full_time_detail_now_ms() - th0);
    out.inputs_sha256 = chem8::inputs_digest(cfg, roots("native-v8/chemistry-carriers", n_proj, lc), n_proj,
                                             roots("native-v8/chemistry-cells", n_occ, lt), n_occ,
                                             roots("native-v8/chemistry-atlas", n_occ, la));
    break;
  }
  out.keep = keep;
  out.device_wall_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  return out;
}
namespace {
// (count, offset, rank) of one segment (three small reads).
struct SegInfo {
  std::uint32_t count = 0;
  u64 offset = 0, rank = 0;
};
SegInfo seg_info(FullTimeDevice::Impl &I, const ChemKeep &K, u64 seg) {
  SegInfo r;
  if (seg >= K.n_seg)
    return r;
  u64 w[3] = {0, 0, 0};
  auto b = take(I.pool, 24);
  // seg_count is u32: copy the aligned u64 word holding it
  I.copy_words(1, K.seg_count + (seg / 2) * 8, b->address);
  I.copy_words(1, K.seg_off + seg * 8, b->address + 8);
  I.copy_words(1, K.seg_rank + seg * 8, b->address + 16);
  I.download(w, *b, 0, 24);
  r.count = (seg & 1) ? std::uint32_t(w[0] >> 32) : std::uint32_t(w[0] & 0xffffffffu);
  r.offset = w[1];
  r.rank = w[2];
  return r;
}
} // namespace
chem8::ESum device_chem_inventory(FullTimeDevice &dev, const ChemDeviceStage &stg, std::uint32_t cell,
                                  std::uint32_t species) {
  auto &I = *dev.impl_;
  std::lock_guard<std::mutex> lock(I.m);
  const auto &K = *std::static_pointer_cast<ChemKeep>(stg.keep);
  const u64 seg = u64(cell) * 22 + species;
  const auto si = seg_info(I, K, seg);
  if (!si.count)
    return {};
  std::array<std::int64_t, pf7::acc_digits> d;
  auto b = take(I.pool, sizeof d);
  I.copy_words(pf7::acc_digits, K.inv + si.rank * pf7::acc_digits * 8, b->address);
  I.download(d.data(), *b, 0, sizeof d);
  return from_digits(d.data());
}
std::vector<ChemSegmentRow> device_chem_segment(FullTimeDevice &dev, const ChemDeviceStage &stg, std::uint32_t cell,
                                                std::uint32_t species) {
  auto &I = *dev.impl_;
  std::lock_guard<std::mutex> lock(I.m);
  const auto &K = *std::static_pointer_cast<ChemKeep>(stg.keep);
  const u64 seg = u64(cell) * 22 + species;
  std::vector<ChemSegmentRow> out;
  const auto si = seg_info(I, K, seg);
  if (!si.count)
    return out;
  std::uint32_t cnt = si.count;
  auto b = take(I.pool, u64(cnt) * 88);
  u64 sr = K.seg_rows + 4 * si.offset, rw = K.rows, ky = K.keys, fl = K.flags, o = b->address;
  void *args[] = {&cnt, &sr, &rw, &ky, &fl, &o};
  I.s.launch(K.mod->k["ch_gather_seg"], {(cnt + 255) / 256}, {256}, 0, args);
  std::vector<u64> w(u64(cnt) * 11);
  I.download(w.data(), *b, 0, w.size() * 8);
  for (std::uint32_t j = 0; j < cnt; ++j)
    out.push_back({std::uint32_t(w[11 * j]), std::uint32_t(w[11 * j + 3]), w[11 * j + 1], w[11 * j + 2],
                   std::bit_cast<double>(w[11 * j + 4]),
                   {std::bit_cast<double>(w[11*j+5]), std::bit_cast<double>(w[11*j+6]), std::bit_cast<double>(w[11*j+7])},
                   {std::bit_cast<double>(w[11*j+8]), std::bit_cast<double>(w[11*j+9]), std::bit_cast<double>(w[11*j+10])}});
  std::sort(out.begin(), out.end(), [](const auto &a, const auto &c) { return a.id < c.id; });
  return out;
}
std::shared_ptr<const DeviceRows> device_chem_apply(FullTimeDevice &dev, const DeviceRows &base, const ChemDeviceEdits &e) {
  auto &I = *dev.impl_;
  std::lock_guard<std::mutex> lock(I.m);
  std::shared_ptr<ChemModule> mod;
  {
    std::lock_guard<std::mutex> l2(chem_modules_m);
    mod = chem_program_locked(I);
  }
  auto &M = *mod;
  if (!base.columns)
    throw std::logic_error("native v8 chemistry apply without census columns");
  ChemKeep K;
  const u64 n = base.n;
  auto blk = [&](std::size_t bytes) {
    auto b = take(I.pool, std::max<std::size_t>(bytes, 8));
    K.blocks.push_back(b);
    return b->address;
  };
  auto launch = [&](const char *name, u64 cnt, void **args) {
    if (cnt)
      I.s.launch(M.k[name], {unsigned((cnt + 255) / 256)}, {256}, 0, args);
  };
  Scan scan{I, M, K};
  const u64 work = blk(n * 80);
  I.copy_words(10 * n, base.rows->address, work);
  if (!e.velocity.empty()) {
    std::vector<u64> v;
    for (const auto &[r, u] : e.velocity) {
      v.push_back(r);
      for (double x : u)
        v.push_back(std::bit_cast<u64>(x));
    }
    auto up = I.upload(v.data(), v.size() * 8);
    K.blocks.push_back(up);
    u64 m = e.velocity.size(), ed = up->address, w = work;
    void *args[] = {&m, &ed, &w};
    launch("ch_edit_velocity", m, args);
  }
  if (!e.weight.empty()) {
    std::vector<u64> v;
    for (const auto &[r, w] : e.weight) {
      v.push_back(r);
      v.push_back(std::bit_cast<u64>(w));
    }
    auto up = I.upload(v.data(), v.size() * 8);
    K.blocks.push_back(up);
    u64 m = e.weight.size(), ed = up->address, w = work;
    void *args[] = {&m, &ed, &w};
    launch("ch_edit_weight", m, args);
  }
  const u64 keep = blk(n * 4);
  {
    u64 k = keep, nn = n;
    std::uint32_t one = 1;
    void *args[] = {&k, &nn, &one};
    launch("ch_fill_u32", n, args);
  }
  if (!e.drop.empty()) {
    std::vector<u64> v(e.drop.begin(), e.drop.end());
    auto up = I.upload(v.data(), v.size() * 8);
    K.blocks.push_back(up);
    u64 m = v.size(), l = up->address, k = keep;
    void *args[] = {&m, &l, &k};
    launch("ch_drop", m, args);
  }
  const u64 rank = blk(n * 8), grand = blk(8);
  scan(keep, n, rank, grand);
  u64 n_keep = 0;
  {
    auto tmp = take(I.pool, 8);
    I.copy_words(1, grand, tmp->address);
    I.download(&n_keep, *tmp, 0, 8);
  }
  const u64 n_out = n_keep + e.append.size();
  auto rows_out = take(I.pool, std::max<u64>(n_out * 80, 8));
  auto c = std::make_shared<DeviceColumns>();
  c->keys = take(I.pool, std::max<u64>(n_out * 8, 8));
  c->states = take(I.pool, std::max<u64>(n_out * 8, 8));
  c->n = n_out;
  {
    u64 nn = n, k = keep, rk = rank, rw = work, ky = base.columns->keys->address, stt = base.columns->states->address,
        orow = rows_out->address, okey = c->keys->address, ost = c->states->address;
    void *args[] = {&nn, &k, &rk, &rw, &ky, &stt, &orow, &okey, &ost};
    launch("ch_compact", n, args);
  }
  if (!e.append.empty()) {
    auto put = [&](const std::shared_ptr<DeviceBlock> &dst, u64 offset, const void *src, std::size_t bytes) {
      std::memcpy(I.stage(bytes), src, bytes);
      I.s.upload(view(*dst), offset, I.staging, 0, bytes);
      I.s.synchronize("native v8 chemistry append");
    };
    put(rows_out, n_keep * 80, e.append.data(), e.append.size() * 80);
    put(c->keys, n_keep * 8, e.append_keys.data(), e.append_keys.size() * 8);
    put(c->states, n_keep * 8, e.append_states.data(), e.append_states.size() * 8);
  }
  auto ids = take(I.pool, pf_v7::gpu::Forms::id_root_scratch_bytes(std::max<u64>(n_out, 1)));
  auto root = take(I.pool, 32);
  I.f.id_root(0, n_out, rows_out->address, 10, ids->address, root->address);
  unsigned char rb[32];
  I.download(rb, *root, 0, 32);
  auto r = std::make_shared<DeviceRows>();
  r->dev = dev.impl_;
  r->rows = rows_out;
  r->columns = c;
  r->n = n_out;
  r->carrier_root = hex_words(rb);
  return r;
}
} // namespace raikou::mfem_backend
