// native_device.cpp -- native v7.1 device lane (API-071.txt section 2).
// Runs on the caller's hall-gpu-resident-v1 Session (one stream) and pf_v7
// forms module; compiles native_kernels.cu with forms_core.h through NVRTC.
#include "native_device.hpp"
#include "forms_gpu.hpp"
#include "hall_gpu_common.hpp"
#include <algorithm>
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
  hall_gpu::Kernel admit, fail_finish, copy;
  std::shared_ptr<DevicePool> pool;
  hall_gpu::HostBuffer staging;
  std::size_t staging_bytes = 0;
  FullTimeDeviceCounters counters;
  std::mutex m;
  Impl(hall_gpu::Session &session, pf_v7::gpu::Forms &forms) : s(session), f(forms) {
    hall_gpu::CompileRequest req;
    req.main = {"native_kernels.cu", kNativeKernels};
    req.headers = {{"forms_core.h", kCore}};
    module = s.compile_and_load(req);
    admit = s.kernel(module, "pf7n_admit");
    fail_finish = s.kernel(module, "pf7n_fail_finish");
    copy = s.kernel(module, "pf7n_copy");
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
  if (n_e != n || n_p != n)
    throw std::logic_error("native v7.1 device admission dimensions");
  const u64 ends = array(D.endpoints, r.endpoints.data(), n * sizeof(GapParticle));
  const u64 paths = array(D.paths, r.paths.data(), n * sizeof(FullTimeMotionPath));
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
  auto after = take(I.pool, nn * sizeof(GapParticle));
  // results: fail[2] | 6 SumOut | E^ (16) | 8 roots | 2 id roots
  constexpr std::size_t so = sizeof(pf_v7::gpu::SumOut);
  constexpr std::size_t off_sums = 16, off_bound = off_sums + 6 * so, off_roots = off_bound + 16,
                        off_ids = off_roots + 8 * 32, results_bytes = off_ids + 2 * 32;
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
  }
  // ---- fused block lists: endpoints, paths, sinks, births, census, E, B, J
  std::vector<pf_v7::gpu::Forms::BlockListJob> jobs = {
      {pf7::section_endpoints, n, ends, R + off_roots + 0 * 32},
      {pf7::section_paths, n, paths, R + off_roots + 1 * 32},
      {pf7::section_sinks, 0, ends, R + off_roots + 2 * 32},
      {pf7::section_births, 0, ends, R + off_roots + 3 * 32},
      {pf7::section_census, census_known ? 0 : n, census_rec->address, R + off_roots + 4 * 32},
      {pf7::section_field_E, n_E, E, R + off_roots + 5 * 32},
      {pf7::section_field_B, n_B, B, R + off_roots + 6 * 32},
      {pf7::section_current_J, n_J, J, R + off_roots + 7 * 32}};
  auto bl = take(I.pool, pf_v7::gpu::Forms::block_lists_scratch_bytes(jobs));
  I.f.block_lists(jobs, bl->address);
  // ---- after rows (copy of the endpoints) and the v5 ID-block roots
  I.copy_words(10 * n, ends, after->address);
  // Both v5 roots in one leaf launch (after rows and paths carry the same IDs whenever the per-carrier checks
  // pass; otherwise the admission is refused and the roots are not used).
  auto ids = take(I.pool, pf_v7::gpu::Forms::id_root_pair_scratch_bytes(n));
  I.f.id_root_pair(n, after->address, 10, paths, 8, ids->address, R + off_ids, R + off_ids + 32);
  // ---- one download
  std::vector<unsigned char> res(results_bytes);
  I.download(res.data(), *results, 0, results_bytes);
  u64 fail[2];
  std::memcpy(fail, res.data(), 16);
  out.first_fail = fail[0];
  out.fail_code = std::uint32_t(fail[1]);
  out.carriers_pass = fail[0] == ~u64(0);
  pf_v7::gpu::SumOut sums[6];
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
  }
  std::string roots[8];
  for (int k = 0; k < 8; ++k)
    roots[k] = hex_words(res.data() + off_roots + 32 * k);
  out.sections = {n, n, 0, 0, roots[0], roots[1], roots[2], roots[3]};
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
  }
  auto rows = std::make_shared<DeviceRows>();
  rows->dev = dev.impl_;
  rows->rows = after;
  rows->columns = base.columns;
  rows->n = n;
  rows->carrier_root = out.after_carrier_root;
  out.after = rows;
  ++I.counters.device_admissions;
  return out;
}
} // namespace raikou::mfem_backend
