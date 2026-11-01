// gas_device.cpp -- gas target owner v1 device lane (hall_gpu::Session).
#include "pf_gas_v1/gas_device.hpp"

#include "common/hall_gpu_common.hpp"

#include <algorithm>
#include <cstring>

namespace pf_gas_v1 {
using pfg::i64;
using pfg::u32;
using pfg::u64;

namespace {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc++26-extensions"
const char kKernels[] = {
#embed "../device/gas_kernels.cu"
    , 0};
const char kGasCore[] = {
#embed "../include/pf_gas_v1/gas_core.h"
    , 0};
const char kFormsCore[] = {
#embed PF_GAS_FORMS_CORE
    , 0};
#pragma GCC diagnostic pop
constexpr u64 D = pfg::kDigits;
constexpr u64 PER = pfg::kEntryQ * D;
constexpr unsigned kBlock = 128;
hall_gpu::Dim3 grid_for(u64 n) {
  hall_gpu::Dim3 g;
  g.x = unsigned((n + kBlock - 1) / kBlock);
  if (g.x == 0)
    g.x = 1;
  return g;
}
hall_gpu::Dim3 block() {
  hall_gpu::Dim3 b;
  b.x = kBlock;
  return b;
}
}  // namespace

struct DeviceLane::Impl {
  std::unique_ptr<hall_gpu::Session> s;
  hall_gpu::Module module;
  hall_gpu::Kernel k_rec, k_flow, k_mov, k_add, k_copy, k_rekey, k_norm, k_fb, k_eleaf, k_gleaf;
  pfg::Layout L{};
  u32 S = 0, NP = 0;
  u64 G = 0;
  hall_gpu::DeviceBuffer species, procs, counters;   // counters: C_COUNTERS u64 delta
  // key table (superset) and state
  std::vector<u64> keys;
  u64 cap = 0;
  hall_gpu::DeviceBuffer d_keys, cand_e, acc_e, delta_e, cand_v, acc_v, cand_g, acc_g, delta_g;
  hall_gpu::DeviceBuffer in;   // inputs scratch
  hall_gpu::HostBuffer stage;  // pinned staging
  double kernel_ms = 0;

  void stage_bytes(std::size_t n) {
    if (stage.bytes >= n)
      return;
    if (stage)
      s->release(stage);
    stage = s->alloc_pinned(std::max<std::size_t>(n, 1 << 20), "gas staging");
  }
  void in_bytes(std::size_t n) {
    if (in.bytes >= n)
      return;
    if (in)
      s->release(in);
    in = s->alloc(std::max<std::size_t>(n, 1 << 20), "gas inputs");
  }
  void up(const hall_gpu::DeviceBuffer &dst, std::size_t off, const void *src, std::size_t n) {
    if (!n)
      return;
    stage_bytes(n);
    std::memcpy(stage.pointer, src, n);
    s->upload(dst, off, stage, 0, n);
    s->synchronize("gas upload");
  }
  void down(void *dst, const hall_gpu::DeviceBuffer &src, std::size_t off, std::size_t n) {
    if (!n)
      return;
    stage_bytes(n);
    s->download(stage, 0, src, off, n);
    s->synchronize("gas download");
    std::memcpy(dst, stage.pointer, n);
  }
  void launch(const hall_gpu::Kernel &k, u64 n, void **args) {
    const int a = s->mark("gas kernel");
    s->launch(k, grid_for(n), block(), 0, args);
    const int b = s->mark("gas kernel end");
    s->synchronize("gas kernel");
    kernel_ms += s->elapsed_ms(a, b);
    s->clear_markers();
  }
  pfg::Tables tables() const {
    pfg::Tables T;
    T.species = species.as<const pfg::SpeciesInfo>();
    T.n_species = S;
    T.procs = procs.as<const pfg::ProcessInfo>();
    T.n_proc = NP;
    T.keys = d_keys.as<const u64>();
    T.n_keys = keys.size();
    T.L = L;
    return T;
  }
  void copy(const hall_gpu::DeviceBuffer &dst, const hall_gpu::DeviceBuffer &src, u64 n) {
    if (!n)
      return;
    i64 *d = dst.as<i64>();
    const i64 *x = src.as<const i64>();
    void *args[] = {&d, &x, &n};
    launch(k_copy, n, args);
  }
  void normalise(const hall_gpu::DeviceBuffer &dig, const hall_gpu::DeviceBuffer *val, u64 n) {
    if (!n)
      return;
    i64 *d = dig.as<i64>();
    double *v = val ? val->as<double>() : nullptr;
    unsigned long long *c = counters.as<unsigned long long>();
    void *args[] = {&d, &v, &n, &c};
    launch(k_norm, n, args);
  }
  void grow(u64 need) {
    if (need <= cap)
      return;
    u64 nc = std::max<u64>(need, cap ? 2 * cap : 64);
    auto realloc_digits = [&](hall_gpu::DeviceBuffer &b, u64 per, const char *tag, bool keep) {
      hall_gpu::DeviceBuffer nb = s->alloc(nc * per * 8, tag);
      s->zero(nb, 0, nb.bytes);
      if (keep && b && cap)
        copy(nb, b, cap * per);
      if (b)
        s->release(b);
      b = nb;
    };
    realloc_digits(cand_e, PER, "gas candidate entries", true);
    realloc_digits(acc_e, PER, "gas accepted entries", true);
    realloc_digits(delta_e, PER, "gas delta entries", false);
    realloc_digits(cand_v, pfg::kEntryQ, "gas candidate values", true);
    realloc_digits(acc_v, pfg::kEntryQ, "gas accepted values", true);
    hall_gpu::DeviceBuffer nk = s->alloc(nc * 8, "gas keys");
    if (d_keys)
      s->release(d_keys);
    d_keys = nk;
    cap = nc;
  }
};

DeviceLane::DeviceLane(const Configuration &cfg, const std::vector<pfg::SpeciesInfo> &species,
                       const std::vector<pfg::ProcessInfo> &procs, pfg::Layout L)
    : p_(std::make_unique<Impl>()) {
  hall_gpu::SessionOptions so;
  so.require_uuid = cfg.device_uuid;
  so.nvrtc_library = cfg.nvrtc_library;
  p_->s = std::make_unique<hall_gpu::Session>(so);
  hall_gpu::CompileRequest req;
  req.main = {"gas_kernels.cu", kKernels};
  req.headers = {{"gas_core.h", kGasCore}, {"forms_core.h", kFormsCore}};
  p_->module = p_->s->compile_and_load(req);
  auto &k = *p_;
  k.k_rec = k.s->kernel(k.module, "gas_apply_records");
  k.k_flow = k.s->kernel(k.module, "gas_apply_flows");
  k.k_mov = k.s->kernel(k.module, "gas_apply_movers");
  k.k_add = k.s->kernel(k.module, "gas_add_digits");
  k.k_copy = k.s->kernel(k.module, "gas_copy_digits");
  k.k_rekey = k.s->kernel(k.module, "gas_rekey");
  k.k_norm = k.s->kernel(k.module, "gas_normalise");
  k.k_fb = k.s->kernel(k.module, "gas_feedback");
  k.k_eleaf = k.s->kernel(k.module, "gas_entry_leaves");
  k.k_gleaf = k.s->kernel(k.module, "gas_global_leaves");
  k.L = L;
  k.S = u32(species.size());
  k.NP = L.NP;
  k.G = L.total();
  k.species = k.s->alloc(species.size() * sizeof(pfg::SpeciesInfo), "gas species");
  k.up(k.species, 0, species.data(), species.size() * sizeof(pfg::SpeciesInfo));
  k.procs = k.s->alloc(procs.size() * sizeof(pfg::ProcessInfo), "gas processes");
  k.up(k.procs, 0, procs.data(), procs.size() * sizeof(pfg::ProcessInfo));
  k.counters = k.s->alloc(pfg::C_COUNTERS * 8, "gas counters");
  k.cand_g = k.s->alloc(k.G * D * 8, "gas candidate globals");
  k.acc_g = k.s->alloc(k.G * D * 8, "gas accepted globals");
  k.delta_g = k.s->alloc(k.G * D * 8, "gas delta globals");
  k.s->zero(k.cand_g, 0, k.cand_g.bytes);
  k.s->zero(k.acc_g, 0, k.acc_g.bytes);
  k.grow(64);
  k.s->synchronize("gas init");
}
DeviceLane::~DeviceLane() {
  if (p_ && p_->s) {
    const auto t = p_->s->close();
    (void)t;
  }
}

void DeviceLane::begin_step(std::uint64_t) {
  auto &k = *p_;
  k.copy(k.cand_e, k.acc_e, k.keys.size() * PER);
  k.copy(k.cand_v, k.acc_v, k.keys.size() * pfg::kEntryQ);
  k.copy(k.cand_g, k.acc_g, k.G * D);
}

void DeviceLane::ensure_keys(const std::vector<std::uint64_t> &want) {
  auto &k = *p_;
  std::vector<u64> fresh;
  for (u64 x : want)
    if (!std::binary_search(k.keys.begin(), k.keys.end(), x))
      fresh.push_back(x);
  if (fresh.empty())
    return;
  std::sort(fresh.begin(), fresh.end());
  fresh.erase(std::unique(fresh.begin(), fresh.end()), fresh.end());
  std::vector<u64> keys;
  std::merge(k.keys.begin(), k.keys.end(), fresh.begin(), fresh.end(), std::back_inserter(keys));
  std::vector<long long> from(keys.size(), -1);
  u64 j = 0;
  for (u64 i = 0; i < keys.size(); ++i)
    if (j < k.keys.size() && k.keys[j] == keys[i])
      from[i] = (long long)(j++);
  k.grow(keys.size());
  // rekey candidate and accepted through the delta buffers
  hall_gpu::DeviceBuffer d_from = k.s->alloc(from.size() * 8, "gas rekey map");
  k.up(d_from, 0, from.data(), from.size() * 8);
  u64 nk = keys.size();
  const long long *fp = d_from.as<const long long>();
  auto rekey = [&](hall_gpu::DeviceBuffer &digits, hall_gpu::DeviceBuffer &values) {
    i64 *dst = k.delta_e.as<i64>();
    const i64 *src = digits.as<const i64>();
    void *args[] = {&dst, &src, &fp, &nk};
    k.launch(k.k_rekey, nk * PER, args);
    k.copy(digits, k.delta_e, nk * PER);
    // values: zero for new keys, moved for old ones (recomputed by normalise)
    k.normalise(digits, &values, nk * pfg::kEntryQ);
  };
  rekey(k.cand_e, k.cand_v);
  rekey(k.acc_e, k.acc_v);
  k.s->release(d_from);
  k.keys = std::move(keys);
  k.up(k.d_keys, 0, k.keys.data(), k.keys.size() * 8);
}

std::array<std::uint64_t, pfg::C_COUNTERS> DeviceLane::apply(std::span<const RecoilRecord> recs,
                                                            std::span<const FlowRecord> flows,
                                                            std::span<const MoverRecord> movs) {
  auto &k = *p_;
  const u64 ne = k.keys.size() * PER;
  k.s->zero(k.delta_e, 0, ne * 8);
  k.s->zero(k.delta_g, 0, k.G * D * 8);
  k.s->zero(k.counters, 0, pfg::C_COUNTERS * 8);
  pfg::Tables T = k.tables();
  i64 *de = k.delta_e.as<i64>(), *dg = k.delta_g.as<i64>();
  unsigned long long *dc = k.counters.as<unsigned long long>();
  if (!recs.empty()) {
    k.in_bytes(recs.size_bytes());
    k.up(k.in, 0, recs.data(), recs.size_bytes());
    const pfg::RecoilRecord *r = k.in.as<const pfg::RecoilRecord>();
    u64 n = recs.size();
    void *args[] = {&T, &r, &n, &de, &dg, &dc};
    k.launch(k.k_rec, n, args);
  }
  if (!flows.empty()) {
    k.in_bytes(flows.size_bytes());
    k.up(k.in, 0, flows.data(), flows.size_bytes());
    const pfg::FlowRecord *f = k.in.as<const pfg::FlowRecord>();
    u64 n = flows.size();
    void *args[] = {&T, &f, &n, &de, &dg, &dc};
    k.launch(k.k_flow, n, args);
  }
  if (!movs.empty()) {
    k.in_bytes(movs.size_bytes());
    k.up(k.in, 0, movs.data(), movs.size_bytes());
    const pfg::MoverRecord *m = k.in.as<const pfg::MoverRecord>();
    u64 n = movs.size();
    const double *v = k.cand_v.as<const double>();
    void *args[] = {&T, &m, &n, &v, &de, &dg, &dc};
    k.launch(k.k_mov, n, args);
  }
  std::array<u64, pfg::C_COUNTERS> c{};
  k.down(c.data(), k.counters, 0, pfg::C_COUNTERS * 8);
  if (c[pfg::C_STATUS])
    throw Refusal(Refusal::Kind::input, "input: device status " + std::to_string(c[pfg::C_STATUS]));
  {
    i64 *dst = k.cand_e.as<i64>();
    const i64 *src = k.delta_e.as<const i64>();
    u64 n = ne;
    void *args[] = {&dst, &src, &n};
    if (n)
      k.launch(k.k_add, n, args);
  }
  {
    i64 *dst = k.cand_g.as<i64>();
    const i64 *src = k.delta_g.as<const i64>();
    u64 n = k.G * D;
    void *args[] = {&dst, &src, &n};
    k.launch(k.k_add, n, args);
  }
  k.s->zero(k.counters, 0, pfg::C_COUNTERS * 8);
  k.normalise(k.cand_e, &k.cand_v, k.keys.size() * pfg::kEntryQ);
  k.normalise(k.cand_g, nullptr, k.G);
  std::array<u64, pfg::C_COUNTERS> st{};
  k.down(st.data(), k.counters, 0, pfg::C_COUNTERS * 8);
  c[pfg::C_STATUS] |= st[pfg::C_STATUS];
  return c;
}

FeedbackReceipt DeviceLane::feedback(std::span<const std::uint32_t> atlas, std::vector<CellTargetsAbi> &cells,
                                     double volume, bool candidate) {
  auto &k = *p_;
  FeedbackReceipt r;
  if (cells.empty())
    return r;
  const std::size_t ab = atlas.size() * 4, cb = cells.size() * sizeof(CellTargetsAbi);
  const std::size_t off = (ab + 255) / 256 * 256;
  k.in_bytes(off + cb);
  k.up(k.in, 0, atlas.data(), ab);
  k.up(k.in, off, cells.data(), cb);
  k.s->zero(k.counters, 0, pfg::C_COUNTERS * 8);
  const u32 *a = k.in.as<const u32>();
  pfg::CellTargetsAbi *c = reinterpret_cast<pfg::CellTargetsAbi *>(k.in.address + off);
  u64 n = cells.size(), nk = k.keys.size();
  const u64 *keys = k.d_keys.as<const u64>();
  const double *v = (candidate ? k.cand_v : k.acc_v).as<const double>();
  const pfg::SpeciesInfo *sp = k.species.as<const pfg::SpeciesInfo>();
  u32 S = k.S;
  unsigned long long *out = k.counters.as<unsigned long long>();
  void *args[] = {&a, &c, &n, &keys, &nk, &v, &sp, &S, &volume, &out};
  k.launch(k.k_fb, n, args);
  std::array<u64, 4> o{};
  k.down(o.data(), k.counters, 0, 32);
  if (o[3])
    throw Refusal(Refusal::Kind::exactness, "exactness: nonrealizable gas feedback");
  k.down(cells.data(), k.in, off, cb);
  r.eligible = o[0];
  r.changed = o[1];
  r.clamped = o[2];
  return r;
}

void DeviceLane::commit() {
  auto &k = *p_;
  std::swap(k.cand_e, k.acc_e);
  std::swap(k.cand_v, k.acc_v);
  std::swap(k.cand_g, k.acc_g);
}
void DeviceLane::abandon() {}

void DeviceLane::leaves(bool candidate, std::vector<Leaf> &el, std::vector<Leaf> &gl) {
  auto &k = *p_;
  el.clear();
  gl.clear();
  const u64 nk = k.keys.size(), G = k.G;
  const std::size_t lb = (nk + G) * 32, fb = (nk + G) * 4;
  k.in_bytes(lb + fb);
  unsigned int *leaf = k.in.as<unsigned int>();
  unsigned int *nz = reinterpret_cast<unsigned int *>(k.in.address + lb);
  if (nk) {
    const i64 *d = (candidate ? k.cand_e : k.acc_e).as<const i64>();
    const u64 *keys = k.d_keys.as<const u64>();
    u64 n = nk;
    void *args[] = {&d, &keys, &n, &leaf, &nz};
    k.launch(k.k_eleaf, n, args);
  }
  {
    const i64 *d = (candidate ? k.cand_g : k.acc_g).as<const i64>();
    unsigned int *gleaf = leaf + 8 * nk;
    unsigned int *gnz = nz + nk;
    u64 n = G;
    void *args[] = {&d, &n, &gleaf, &gnz};
    k.launch(k.k_gleaf, n, args);
  }
  std::vector<unsigned int> L((nk + G) * 8), F(nk + G);
  k.down(L.data(), k.in, 0, lb);
  k.down(F.data(), k.in, lb, fb);
  for (u64 i = 0; i < nk + G; ++i) {
    if (!F[i])
      continue;
    Leaf w;
    std::memcpy(w.data(), L.data() + 8 * i, 32);
    (i < nk ? el : gl).push_back(w);
  }
}

void DeviceLane::download(bool candidate, std::vector<std::uint64_t> &keys, std::vector<std::int64_t> &entry,
                          std::vector<std::int64_t> &global) {
  auto &k = *p_;
  std::vector<i64> all(k.keys.size() * PER);
  k.down(all.data(), candidate ? k.cand_e : k.acc_e, 0, all.size() * 8);
  keys.clear();
  entry.clear();
  for (u64 i = 0; i < k.keys.size(); ++i) {
    const i64 *d = all.data() + i * PER;
    if (pfg::digits_zero(d, int(PER)))
      continue;
    keys.push_back(k.keys[i]);
    entry.insert(entry.end(), d, d + PER);
  }
  global.assign(k.G * D, 0);
  k.down(global.data(), candidate ? k.cand_g : k.acc_g, 0, global.size() * 8);
}

void DeviceLane::load(const std::vector<std::uint64_t> &keys, const std::vector<std::int64_t> &entry,
                      const std::vector<std::int64_t> &global) {
  auto &k = *p_;
  k.keys.clear();
  k.grow(std::max<u64>(keys.size(), 1));
  k.keys = keys;
  k.up(k.d_keys, 0, k.keys.data(), k.keys.size() * 8);
  k.up(k.acc_e, 0, entry.data(), entry.size() * 8);
  k.up(k.acc_g, 0, global.data(), global.size() * 8);
  k.normalise(k.acc_e, &k.acc_v, k.keys.size() * pfg::kEntryQ);
  k.normalise(k.acc_g, nullptr, k.G);
}

std::string DeviceLane::session_json() const {
  return "{\"device\":" + hall_gpu::device_json(p_->s->device()) + ",\"compile\":" +
         hall_gpu::compile_json(p_->module.record) + "}";
}
double DeviceLane::kernel_ms_total() const noexcept { return p_->kernel_ms; }

}  // namespace pf_gas_v1
