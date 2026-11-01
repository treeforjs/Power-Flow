// Heat/source owner v7: deterministic NVRTC wall kernel against the CPU reference, bit for bit.
//   test_wall_gpu_v7 NVRTC_LIBRARY DEVICE_UUID V7_NATIVE_DIR ARRIVAL_V1_INCLUDE OUTDIR
// The kernel source (device/sw7_kernels.cu, sw7_wall.h, sw7_math.h) and pf_arrival/philox_core.h are read from disk
// and compiled by hall_gpu::Session (hall-gpu-resident-v1 common layer: NVRTC dlopen'ed, --fmad=false --prec-div=true
// --prec-sqrt=true --ftz=false, architecture selected from the device at run time).
// Groups:
//   layout     sizeof/offsets of SW7Params, SW7Hit, SW7Out, SW7Species equal on host and device
//   parity     1e6 wall events (the test_surface_v7 event set and patch states): every SW7Out byte equal
//   reduction  fixed-tree block sums + sequential final pass (integer limbs, no atomics) equal the CPU int128 sums
//   repeat     a second run and a run with 128-thread blocks give identical outputs and sums
#include "common/hall_gpu_common.hpp"
#include "surface_data_v7.hpp"
#include "event_set_v7.hpp"
#include "raikou/json.h"
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

namespace sv = surface_v7;
using J = raikou::Json;

namespace {
struct Results {
  J cases = J::array();
  int failed = 0;
  void record(const std::string& g, const std::string& name, bool ok, J ev = J::object()) {
    cases.push_back({{"group", g}, {"name", name}, {"passed", ok}, {"evidence", std::move(ev)}});
    if (!ok) ++failed;
    std::cout << g << " | " << name << ": " << (ok ? "pass" : "FAIL") << std::endl;
    if (!ok) std::cout << "  evidence: " << cases.back()["evidence"].dump().substr(0, 3000) << std::endl;
  }
};
std::string slurp(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error("cannot read " + p.string());
  return {std::istreambuf_iterator<char>(f), {}};
}
constexpr int NF = 11;
void host_fields(const SW7Out& o, sv::i128* f) {
  f[0] = o.e_in_q; f[1] = o.e_out_q; f[2] = o.wall_q; f[3] = o.k_in_q; f[4] = o.ads_q; f[5] = o.sigma_q;
  f[6] = o.n_out; f[7] = o.electrons_out; f[8] = o.stuck; f[9] = o.charge_in; f[10] = o.charge_out;
}
}  // namespace

int main(int argc, char** argv) {
  Results R;
  int rc = 0;
  try {
    if (argc != 6) throw std::runtime_error("usage: test_wall_gpu_v7 NVRTC_LIBRARY DEVICE_UUID V7_NATIVE_DIR ARRIVAL_V1_INCLUDE OUTDIR");
    const std::filesystem::path v7 = argv[3], arr = argv[4], out = argv[5];
    std::filesystem::create_directories(out);
    hall_gpu::SessionOptions so;
    so.require_uuid = argv[2];
    so.nvrtc_library = argv[1];
    hall_gpu::Session s(so);
    hall_gpu::CompileRequest req;
    req.main = {"sw7_kernels.cu", slurp(v7 / "device/sw7_kernels.cu")};
    req.headers = {{"sw7_wall.h", slurp(v7 / "device/sw7_wall.h")}, {"sw7_math.h", slurp(v7 / "device/sw7_math.h")},
                   {"pf_arrival/philox_core.h", slurp(arr / "pf_arrival/philox_core.h")}};
    auto mod = s.compile_and_load(req);
    const auto& cr = mod.record;
    J dev = {{"name", s.device().name}, {"uuid", s.device().uuid}, {"cc", std::to_string(s.device().cc_major) + "." + std::to_string(s.device().cc_minor)},
             {"driver", s.device().kernel_module_version}};
    J comp = {{"nvrtc", std::to_string(cr.nvrtc_major) + "." + std::to_string(cr.nvrtc_minor)}, {"arch_option", cr.arch_option}, {"image_kind", cr.image_kind},
              {"options", cr.options}, {"image_sha256", cr.image_sha256}, {"nvrtc_library_sha256", cr.nvrtc_library_sha256}};
    auto k_sizes = s.kernel(mod, "sw7_sizes"), k_events = s.kernel(mod, "sw7_wall_events"), k_red = s.kernel(mod, "sw7_reduce_blocks"), k_fin = s.kernel(mod, "sw7_reduce_final");

    // ---- layout
    {
      auto d = s.alloc(16 * 8, "sizes");
      auto h = s.alloc_pinned(16 * 8, "sizes-host");
      s.zero(d, 0, 16 * 8);
      void* args[] = {&d.address};
      s.launch(k_sizes, {1, 1, 1}, {1, 1, 1}, 0, args);
      s.download(h, 0, d, 0, 16 * 8);
      s.synchronize("sizes");
      const auto* x = h.as<unsigned long long>();
      SW7Out o{};
      SW7Params p{};
      const unsigned long long host[9] = {sizeof(SW7Params), sizeof(SW7Hit), sizeof(SW7Out), sizeof(SW7Species),
                                          (unsigned long long)((char*)&o.k_in_q - (char*)&o), (unsigned long long)((char*)&o.sigma_q - (char*)&o),
                                          (unsigned long long)((char*)&p.pn[0] - (char*)&p), (unsigned long long)((char*)&p.enable_isd - (char*)&p),
                                          (unsigned long long)((char*)&o.u[0][0] - (char*)&o)};
      bool same = true;
      J rows = J::array();
      for (int i = 0; i < 9; ++i) { same = same && x[i] == host[i]; rows.push_back({x[i], host[i]}); }
      R.record("layout", "struct sizes and field offsets equal on host and device", same, {{"device_host", rows}});
      s.release(d);
      s.release(h);
    }

    // ---- inputs and CPU reference (the test_surface_v7 ledger set)
    const std::size_t n = 1000000;
    sv::Registry reg;
    SW7Params base = sv::wall_parameters(reg);
    base.phi_eV = 4.35; base.T_K = 300; base.theta = 0.33; base.adsorption_eV = 0.87;
    const auto params = v7test::patch_params(base);
    const auto table = sv::species_table(nullptr);
    const auto hits = v7test::event_set(n, 20261031, table);
    const std::uint64_t seed = 4242, step = 12345;
    std::vector<SW7Out> cpu(n);
    auto t0 = std::chrono::steady_clock::now();
    for (std::size_t i = 0; i < n; ++i) sw7::wall_event(params[hits[i].patch], table.data(), hits[i], seed, step, &cpu[i]);
    const double cpu_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    sv::i128 cpu_tot[NF] = {};
    for (const auto& o : cpu) { sv::i128 f[NF]; host_fields(o, f); for (int k = 0; k < NF; ++k) cpu_tot[k] += f[k]; }

    // ---- device buffers (batches of 2^17 events)
    const std::size_t B = std::size_t(1) << 17;
    auto d_par = s.alloc(sizeof(SW7Params) * params.size(), "params");
    auto d_tab = s.alloc(sizeof(SW7Species) * table.size(), "species");
    auto d_hit = s.alloc(sizeof(SW7Hit) * B, "hits");
    auto d_out = s.alloc(sizeof(SW7Out) * B, "outputs");
    const unsigned max_blocks = unsigned((B + 255) / 256);
    auto d_part = s.alloc(sizeof(long long) * 2 * NF * max_blocks, "partials");
    auto d_tot = s.alloc(sizeof(long long) * 2 * NF, "totals");
    auto h_par = s.alloc_pinned(sizeof(SW7Params) * params.size(), "params-host");
    auto h_tab = s.alloc_pinned(sizeof(SW7Species) * table.size(), "species-host");
    auto h_hit = s.alloc_pinned(sizeof(SW7Hit) * B, "hits-host");
    auto h_out = s.alloc_pinned(sizeof(SW7Out) * B, "outputs-host");
    auto h_tot = s.alloc_pinned(sizeof(long long) * 2 * NF, "totals-host");
    std::memcpy(h_par.pointer, params.data(), sizeof(SW7Params) * params.size());
    std::memcpy(h_tab.pointer, table.data(), sizeof(SW7Species) * table.size());
    s.upload(d_par, 0, h_par, 0, sizeof(SW7Params) * params.size());
    s.upload(d_tab, 0, h_tab, 0, sizeof(SW7Species) * table.size());
    s.synchronize("static uploads");

    auto run = [&](unsigned block, std::vector<SW7Out>& gpu, sv::i128* tot, double& ms) {
      gpu.assign(n, SW7Out{});
      for (int k = 0; k < NF; ++k) tot[k] = 0;
      ms = 0;
      for (std::size_t first = 0; first < n; first += B) {
        const unsigned m = unsigned(std::min(B, n - first));
        std::memcpy(h_hit.pointer, hits.data() + first, sizeof(SW7Hit) * m);
        s.upload(d_hit, 0, h_hit, 0, sizeof(SW7Hit) * m);
        s.zero(d_out, 0, sizeof(SW7Out) * B);
        const int a = s.mark("events");
        unsigned mm = m;
        unsigned long long sd = seed, st = step;
        void* args[] = {&d_par.address, &d_tab.address, &d_hit.address, &mm, &sd, &st, &d_out.address};
        s.launch(k_events, {(m + block - 1) / block, 1, 1}, {block, 1, 1}, 0, args);
        const int b = s.mark("events-done");
        const unsigned blocks = (m + 255) / 256;
        void* rargs[] = {&d_out.address, &mm, &d_part.address};
        s.launch(k_red, {blocks, 1, 1}, {256, 1, 1}, 0, rargs);
        unsigned bb = blocks;
        void* fargs[] = {&d_part.address, &bb, &d_tot.address};
        s.launch(k_fin, {1, 1, 1}, {1, 1, 1}, 0, fargs);
        s.download(h_out, 0, d_out, 0, sizeof(SW7Out) * m);
        s.download(h_tot, 0, d_tot, 0, sizeof(long long) * 2 * NF);
        s.synchronize("batch");
        ms += s.elapsed_ms(a, b);
        s.clear_markers();
        std::memcpy(gpu.data() + first, h_out.pointer, sizeof(SW7Out) * m);
        const auto* t = h_tot.as<long long>();
        for (int k = 0; k < NF; ++k) tot[k] += sv::i128(t[2 * k]) * (sv::i128(1) << 32) + sv::i128(t[2 * k + 1]);
      }
    };
    std::vector<SW7Out> g1, g2, g3;
    sv::i128 t1[NF], t2[NF], t3[NF];
    double ms1 = 0, ms2 = 0, ms3 = 0;
    run(256, g1, t1, ms1);
    run(256, g2, t2, ms2);
    run(128, g3, t3, ms3);
    long diff = 0, first_diff = -1;
    for (std::size_t i = 0; i < n; ++i)
      if (std::memcmp(&g1[i], &cpu[i], sizeof(SW7Out)) != 0) { if (first_diff < 0) first_diff = long(i); ++diff; }
    J pe = {{"events", n}, {"differing_events", diff}, {"first_difference", first_diff}, {"cpu_seconds", cpu_s}, {"gpu_kernel_ms", ms1}, {"device", dev}, {"compile", comp}};
    if (first_diff >= 0) {
      const auto& a = cpu[first_diff];
      const auto& b = g1[first_diff];
      pe["cpu_status_nout_wall"] = {a.status, a.n_out, sv::i128_text(a.wall_q)};
      pe["gpu_status_nout_wall"] = {b.status, b.n_out, sv::i128_text(b.wall_q)};
      pe["species"] = hits[first_diff].species;
    }
    R.record("parity", "1e6 wall events: every SW7Out byte of the NVRTC kernel equals the CPU reference", diff == 0, pe);
    bool sums = true;
    J sv_rows = J::array();
    for (int k = 0; k < NF; ++k) { sums = sums && t1[k] == cpu_tot[k]; sv_rows.push_back({sv::i128_text(t1[k]), sv::i128_text(cpu_tot[k])}); }
    R.record("reduction", "fixed-tree GPU sums (integer limbs, no atomics) of 11 ledger fields equal the CPU int128 sums", sums, {{"gpu_cpu", sv_rows}});
    long rep = 0;
    for (std::size_t i = 0; i < n; ++i) rep += std::memcmp(&g1[i], &g2[i], sizeof(SW7Out)) != 0 || std::memcmp(&g1[i], &g3[i], sizeof(SW7Out)) != 0;
    bool same_tot = true;
    for (int k = 0; k < NF; ++k) same_tot = same_tot && t1[k] == t2[k] && t1[k] == t3[k];
    R.record("repeat", "repeated run and 128-thread blocks: identical outputs and sums", rep == 0 && same_tot, {{"differing", rep}, {"ms_256", ms1}, {"ms_256_again", ms2}, {"ms_128", ms3}});
    for (auto* b : {&d_par, &d_tab, &d_hit, &d_out, &d_part, &d_tot}) s.release(*b);
    for (auto* b : {&h_par, &h_tab, &h_hit, &h_out, &h_tot}) s.release(*b);
    const auto td = s.close();
    R.record("teardown", "session teardown drained and released", td.released && !td.quarantined, {{"first_error", td.first_error}});
    std::ofstream(out / "wall-gpu-v7.json") << J{{"schema", "wall-gpu-v7-qualification-v1"}, {"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    rc = R.failed == 0 ? 0 : 1;
  } catch (const std::exception& e) {
    std::cerr << "fatal: " << e.what() << '\n';
    rc = 2;
  }
  return rc;
}
