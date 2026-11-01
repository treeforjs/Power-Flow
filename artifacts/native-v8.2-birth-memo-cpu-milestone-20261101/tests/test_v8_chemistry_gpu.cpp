// Native v8 chemistry v3 device stage controls (API-008.txt 3.7):
//   CPU/GPU bit identity  host, device_verified and device runs end in the same
//                         records roots, ledgers, chemistry states and accepted identities
//   partitions            device launches over 1/2/5/8 partitions give identical results
//   refusal parity        a coverage-policy refusal has the same typed text on the host and the device
//   regression            device_verified admission of a device-only record (v7.1 defect fixed in v8:
//                         tests/time_step_v71.cpp device_verified mode)
// usage: test_v8_chemistry_gpu OUT_DIR NVRTC UUID
#include "chem_v8_fixture.hpp"
#include "audit_cold_support.hpp"
#include "forms_gpu.hpp"
#include "hall_gpu_common.hpp"
#include "raikou/mfem_pic_full_time_device.h"
#include "raikou/json.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

namespace m = raikou::mfem_backend;
using J = nlohmann::json;
namespace {
int checks = 0, failures = 0;
void check(bool ok, const std::string &what, const J &d = J::object()) {
  ++checks;
  if (!ok) {
    ++failures;
    std::cout << "FAIL " << what << " " << d.dump() << "\n";
  }
}
chemfix::Population dense(std::size_t n, std::uint64_t seed) {
  chemfix::Population P;
  P.carriers = n;
  P.seed = seed;
  P.water_weight = 2e13;
  P.clustered = true;
  P.electron_eV_min = 5;
  P.electron_eV_max = 300;
  return P;
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 4)
      throw std::runtime_error("usage: test_v8_chemistry_gpu OUT_DIR NVRTC UUID");
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    hall_gpu::SessionOptions so;
    so.require_uuid = argv[3];
    so.nvrtc_library = argv[2];
    auto session = std::make_unique<hall_gpu::Session>(so);
    auto forms = std::make_unique<pf_v7::gpu::Forms>(*session);
    auto dev = std::make_shared<m::FullTimeDevice>(*session, *forms);
    const auto o = chemfix::options(45);
    const double dt = 1e-12;
    J res = J::object();
    auto run = [&](m::FullTimeChemistryExecution ex, m::FullTimeExecution lane, unsigned parts, std::size_t n,
                   std::uint64_t seed, int steps) {
      auto chem = chemfix::chemistry_configuration(o);
      chem.execution = ex;
      chem.device_partitions = parts;
      auto A = chemfix::make_arena(o, chem, chemfix::population(dense(n, seed)), 40 * n);
      if (ex != m::FullTimeChemistryExecution::host)
        A.arena->attach_full_time_device(dev, lane);
      std::vector<std::string> trace;
      J st = J::array();
      for (int k = 0; k < steps; ++k) {
        const auto r = chemfix::step(A, dt, k);
        trace.push_back(r.stage1->records_root + "/" + r.stage1->ledger_sha256 + "/" + r.stage2->records_root + "/" +
                        r.stage2->ledger_sha256 + "/" + r.stage2->state_after_sha256 + "/" + r.boundary_sha256);
        st.push_back({{"events", r.stage1->admitted_events + r.stage2->admitted_events},
                      {"device", r.stage1->device}, {"device_ms", r.stage1->device_ms + r.stage2->device_ms},
                      {"census", r.census}});
      }
      return std::pair<std::vector<std::string>, J>{trace, st};
    };
    // 1. host vs device_verified vs device
    const auto H = run(m::FullTimeChemistryExecution::host, m::FullTimeExecution::host, 0, 3000, 11, 3);
    const auto c0 = dev->counters().verified_comparisons;
    const auto V = run(m::FullTimeChemistryExecution::device_verified, m::FullTimeExecution::device_verified, 0, 3000, 11, 3);
    const auto c1 = dev->counters().verified_comparisons;
    const auto D = run(m::FullTimeChemistryExecution::device, m::FullTimeExecution::device, 0, 3000, 11, 3);
    check(H.first == V.first && H.first == D.first, "CPU/GPU: host, device_verified and device traces identical",
          {{"host", H.first}, {"device", D.first}});
    check(c1 - c0 >= 6, "device_verified: every stage compared with the host path", {{"comparisons", c1 - c0}});
    bool on_device = true;
    for (const auto &x : D.second)
      on_device = on_device && x.at("device").get<bool>();
    check(on_device, "device mode: every stage ran on the device");
    res["host"] = H.second;
    res["device"] = D.second;
    res["trace"] = H.first;
    // 2. device partitions 1/2/5/8 (larger population: several 256-carrier launches)
    std::vector<std::vector<std::string>> T;
    for (unsigned P : {1u, 2u, 5u, 8u})
      T.push_back(run(m::FullTimeChemistryExecution::device, m::FullTimeExecution::device, P, 6000, 23, 2).first);
    const auto Hp = run(m::FullTimeChemistryExecution::host, m::FullTimeExecution::host, 0, 6000, 23, 2).first;
    bool same = true;
    for (const auto &t : T)
      same = same && t == Hp;
    check(same, "device partitions 1/2/5/8 equal the host path");
    // 3. refusal parity: coverage policy (nothing admitted)
    std::string texts[2];
    for (int k = 0; k < 2; ++k) {
      auto chem = chemfix::chemistry_configuration(o);
      chem.policy.admitted_regions.clear();
      chem.execution = k ? m::FullTimeChemistryExecution::device : m::FullTimeChemistryExecution::host;
      auto A = chemfix::make_arena(o, chem, chemfix::population(dense(800, 51)), 8000);
      if (k)
        A.arena->attach_full_time_device(dev, m::FullTimeExecution::device);
      try {
        (void)chemfix::step(A, dt, 0);
        texts[k] = "none";
      } catch (const m::FullTimeChemistryRefusal &e) {
        texts[k] = e.what();
      }
    }
    check(texts[0] == texts[1] && texts[0].find("chemistry_v3/coverage_policy") == 0,
          "refusal parity: device and host give the same typed coverage refusal", {{"host", texts[0]}, {"device", texts[1]}});
    check(audit::cold_support(out/"cold-device",m::FullTimeChemistryExecution::device,
          [&](m::PicArena &a){a.attach_full_time_device(dev,m::FullTimeExecution::device);})==7,
          "cold reachable-support strict refusal, no mutation, diagnostic zero rate and restore");
    const auto cnt = dev->counters();
    res["counters"] = {{"verified_comparisons", cnt.verified_comparisons}, {"row_uploads", cnt.row_uploads},
                       {"row_downloads", cnt.row_downloads}, {"device_admissions", cnt.device_admissions}};
    res["checks"] = checks;
    res["failures"] = failures;
    std::ofstream(out / "v8-chemistry-gpu-tests.json") << res.dump(1) << "\n";
    std::cout << (failures ? "FAILED " : "") << checks - failures << " native v8 chemistry device controls passed\n";
    dev.reset();
    forms.reset();
    (void)session->close();
    return failures ? 1 : 0;
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 2;
  }
}
