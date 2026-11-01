// Native v8 nonempty full-time restore (API-008.txt 5).
// usage: test_v8_restore OUT_DIR unbroken|save|restore|tamper [STEPS]
//   unbroken  2K steps from genesis; trace of every accepted boundary
//   save      K steps, then the checkpoint (OUT_DIR/ckpt.json) and its boundary
//   restore   (fresh process) restore the checkpoint, K more steps; trace
//   tamper    truncated, flipped byte, mismatched window row, mismatched anchor,
//             wrong schema version (the last three with a recomputed envelope
//             digest so the deeper checks are reached): every one refused
// Chemistry v3 (host path), anchors every 2 steps.
#include "chem_v8_fixture.hpp"
#include "raikou/json.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <openssl/sha.h>
#include <sstream>

namespace m = raikou::mfem_backend;
using J = nlohmann::json;
namespace {
chemfix::Population pop() {
  chemfix::Population P;
  P.carriers = 1500;
  P.seed = 71;
  P.water_weight = 2e13;
  P.clustered = true;
  P.electron_eV_min = 5;
  P.electron_eV_max = 300;
  return P;
}
J boundary_json(const m::FullTimeAcceptedBoundary &b) {
  return {{"accepted_step", b.accepted_step}, {"revision", b.accepted_arena_revision}, {"time_s", b.time_s},
          {"arena", b.accepted_arena_sha256},  {"field", b.field_sha256},             {"heat", b.heat_sha256},
          {"source", b.source_sha256},         {"rng", b.rng_sha256},                 {"mesh", b.mesh_sha256},
          {"material", b.material_sha256},     {"frame", b.frame_sha256},             {"clock", b.common_clock_sha256},
          {"genesis", b.authored_full_time_genesis}};
}
m::FullTimeAcceptedBoundary boundary_of(const J &j) {
  m::FullTimeAcceptedBoundary b;
  b.accepted_step = j.at("accepted_step");
  b.accepted_arena_revision = j.at("revision");
  b.time_s = j.at("time_s");
  b.accepted_arena_sha256 = j.at("arena");
  b.field_sha256 = j.at("field");
  b.heat_sha256 = j.at("heat");
  b.source_sha256 = j.at("source");
  b.rng_sha256 = j.at("rng");
  b.mesh_sha256 = j.at("mesh");
  b.material_sha256 = j.at("material");
  b.frame_sha256 = j.at("frame");
  b.common_clock_sha256 = j.at("clock");
  b.authored_full_time_genesis = j.at("genesis");
  return b;
}
std::string sha_hex(const std::string &s) {
  unsigned char d[32];
  SHA256(reinterpret_cast<const unsigned char *>(s.data()), s.size(), d);
  std::ostringstream o;
  for (auto x : d)
    o << std::hex << std::setfill('0') << std::setw(2) << unsigned(x);
  return o.str();
}
const double dt = 1e-12;
// ejected-melt credit of step k (phase owner v1.2 eject_N/E/K; exact cumulative reservoir)
void melt(m::PrivateFullTimeStep &outer, int k) {
  m::FullTimeEjectedMelt e;
  e.atoms = {1e12 * (k + 1), 3e11, 1.5e11, 2e10, 1e10, 5e9};
  e.energy_J = 1e-3 * (k + 1);
  e.kinetic_J = 1e-5;
  e.phase_record_sha256 = chemfix::sha("phase-record/" + std::to_string(k));
  outer.credit_ejected_melt(e);
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc < 3)
      throw std::runtime_error("usage: test_v8_restore OUT_DIR unbroken|save|restore|tamper [STEPS]");
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    const std::string mode = argv[2];
    const int K = argc > 3 ? std::atoi(argv[3]) : 3;
    const auto o = chemfix::options(45);
    const auto chem = chemfix::chemistry_configuration(o);
    if (mode == "unbroken" || mode == "save") {
      auto A = chemfix::make_arena(o, chem, chemfix::population(pop()), 40000, 2);
      J trace = J::array();
      const int steps = mode == "unbroken" ? 2 * K : K;
      for (int k = 0; k < steps; ++k) {
        const auto r = chemfix::step(A, dt, k, {}, [&](m::PrivateFullTimeStep &o2) { melt(o2, k); });
        trace.push_back({{"step", k}, {"boundary", r.boundary_sha256}, {"rng", r.stage2->state_after_sha256},
                         {"events", r.stage1->admitted_events + r.stage2->admitted_events},
                         {"melt", A.arena->full_time_chemistry().melt.sha256},
                         {"stage1_export",r.stage1->admitted_export->sha256},
                         {"stage2_export",r.stage2->admitted_export->sha256}});
      }
      if (mode == "save") {
        A.arena->save((out / "ckpt.json").string());
        std::ofstream(out / "boundary.json") << boundary_json(*A.arena->full_time_boundary()).dump(1) << "\n";
      }
      std::ofstream(out / (mode + "-trace.json")) << trace.dump(1) << "\n";
      std::cout << mode << " " << steps << " steps\n";
      return 0;
    }
    J bj;
    std::ifstream(out / "boundary.json") >> bj;
    const auto expected = boundary_of(bj);
    auto A = chemfix::make_arena(o, chem, chemfix::population(pop()), 40000, 2); // catalogue/configuration only
    if (mode == "restore") {
      A.arena = m::PicArena::restore_full_time((out / "ckpt.json").string(), o, A.states, A.channels, expected);
      J trace = J::array();
      for (int k = K; k < 2 * K; ++k) {
        const auto r = chemfix::step(A, dt, k, {}, [&](m::PrivateFullTimeStep &o2) { melt(o2, k); });
        trace.push_back({{"step", k}, {"boundary", r.boundary_sha256}, {"rng", r.stage2->state_after_sha256},
                         {"events", r.stage1->admitted_events + r.stage2->admitted_events},
                         {"melt", A.arena->full_time_chemistry().melt.sha256},
                         {"stage1_export",r.stage1->admitted_export->sha256},
                         {"stage2_export",r.stage2->admitted_export->sha256}});
      }
      std::ofstream(out / "restore-trace.json") << trace.dump(1) << "\n";
      J unbroken;
      std::ifstream(out / "unbroken-trace.json") >> unbroken;
      bool same = true;
      for (int k = 0; k < K; ++k)
        same = same && trace[k] == unbroken[K + k];
      std::cout << (same ? "restore continuation identical to the unbroken run" : "FAIL restore continuation differs")
                << "\n";
      return same ? 0 : 1;
    }
    if (mode == "tamper") {
      std::string text;
      {
        std::ifstream f(out / "ckpt.json", std::ios::binary);
        std::stringstream ss;
        ss << f.rdbuf();
        text = ss.str();
      }
      const auto doc = J::parse(text);
      auto resign = [&](J payload) {
        return J{{"payload", payload}, {"sha256", sha_hex(payload.dump())}}.dump();
      };
      struct Case {
        std::string name, bytes;
      };
      std::vector<Case> cases;
      cases.push_back({"truncated", text.substr(0, text.size() / 2)});
      {
        auto t = text;
        t[t.size() / 3] = t[t.size() / 3] == '0' ? '1' : '0';
        cases.push_back({"flipped byte", t});
      }
      {
        auto p = doc.at("payload");
        auto &w = p["metadata"]["history"]["window"];
        if (!w.empty()) {
          auto h = w[0]["canonical_hex"].get<std::string>();
          h[h.size() - 2] = h[h.size() - 2] == '0' ? '1' : '0';
          w[0]["canonical_hex"] = h;
        }
        cases.push_back({"mismatched window row (re-signed)", resign(p)});
      }
      {
        auto p = doc.at("payload");
        auto &a = p["metadata"]["history"]["anchors"];
        if (!a.empty()) {
          auto h = a[0]["canonical_hex"].get<std::string>();
          h[h.size() - 2] = h[h.size() - 2] == '0' ? '1' : '0';
          a[0]["canonical_hex"] = h;
        }
        cases.push_back({"mismatched anchor (re-signed)", resign(p)});
      }
      {
        auto p = doc.at("payload");
        p["configuration"]["v8"]["identity"] = "v8.2-auditfix2";
        cases.push_back({"wrong schema version (re-signed)", resign(p)});
      }
      {
        auto p = doc.at("payload");
        auto &c = p["metadata"]["chemistry"]["state_hex"];
        auto h = c.get<std::string>();
        h[h.size() / 2] = h[h.size() / 2] == '0' ? '1' : '0';
        c = h;
        cases.push_back({"chemistry reservoir/ledger state changed (re-signed)", resign(p)});
      }
      {
        auto p = doc.at("payload");
        auto &rows = p["native"];
        cases.push_back({"native rows changed (re-signed)", [&] {
                           if (rows.is_object() && rows.contains("native_bytes_hex")) {
                             auto h = rows["native_bytes_hex"].get<std::string>();
                             h[h.size() / 2] = h[h.size() / 2] == '0' ? '1' : '0';
                             rows["native_bytes_hex"] = h;
                           } else if (rows.is_array() && !rows.empty())
                             rows.erase(rows.begin());
                           return resign(p);
                         }()});
      }
      {
        auto p = doc.at("payload");
        J saved_trace;
        { std::ifstream f(out / "save-trace.json"); f >> saved_trace; }
        const auto committed = saved_trace.back().at("stage1_export").get<std::string>();
        std::string hex;
        constexpr char digits[]="0123456789abcdef";
        for(unsigned char ch:committed) {hex+=digits[ch>>4];hex+=digits[ch&15];}
        bool changed=false;
        for(auto &row:p["metadata"]["history"]["window"]) {
          auto h=row.at("canonical_hex").get<std::string>();
          auto pos=h.find(hex);
          if(pos!=std::string::npos) {
            h.replace(pos,2,hex.substr(0,2)=="30"?"31":"30");
            row["canonical_hex"]=h;changed=true;break;
          }
        }
        if(!changed) throw std::runtime_error("retained chemistry export commitment fixture absent");
        cases.push_back({"changed admitted-export commitment (re-signed envelope)",resign(p)});
      }
      J res = J::array();
      int refused = 0;
      for (const auto &c : cases) {
        const auto path = out / "tampered.json";
        std::ofstream(path, std::ios::binary) << c.bytes;
        std::string what = "none";
        try {
          (void)m::PicArena::restore_full_time(path.string(), o, A.states, A.channels, expected);
        } catch (const std::exception &e) {
          what = e.what();
        }
        refused += what != "none";
        res.push_back({{"case", c.name}, {"refusal", what}});
        std::filesystem::remove(path);
      }
      std::ofstream(out / "tamper.json") << res.dump(1) << "\n";
      std::cout << refused << "/" << cases.size() << " tampered checkpoints refused\n";
      return refused == int(cases.size()) ? 0 : 1;
    }
    throw std::runtime_error("unknown mode");
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 2;
  }
}
