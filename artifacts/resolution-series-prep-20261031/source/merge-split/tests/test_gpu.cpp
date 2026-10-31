// test_gpu.cpp -- device transaction: CPU/GPU bit identity (every output word, every
// ledger term, the digest), repeat-run identity on the device, timing per 1e6 carriers.
// Usage: test_gpu --nvrtc LIB --uuid UUID [--n N] [--meshes 45,90,180,360] [--repeat R]
//                 [--corner GROUPS] [--check-each-launch 0|1] [--out FILE.json]
#include "fixtures.hpp"
#include "merge_split.hpp"
#include "ms_gpu.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace merge_split;
using Clock = std::chrono::steady_clock;

namespace {
int failures = 0;
void check(bool ok, const std::string& what) {
    if (!ok) {
        ++failures;
        std::cout << "FAIL " << what << std::endl;
    }
}
std::string num(double x) {
    char b[64];
    std::snprintf(b, sizeof b, "%.17g", x);
    return b;
}
std::string q(const std::string& s) { return "\"" + s + "\""; }
template <class T> bool same_bytes(const std::vector<T>& a, const std::vector<T>& b) {
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(T)) == 0);
}
// Field-by-field comparison; returns the names of differing parts.
std::string differences(const Result& a, const Result& b) {
    std::string d;
    auto add = [&](bool ok, const char* name) { if (!ok) d += std::string(d.empty() ? "" : ",") + name; };
    add(same_bytes(a.groups, b.groups), "groups");
    add(same_bytes(a.group_first_id, b.group_first_id), "group_first_id");
    add(same_bytes(a.merge, b.merge), "merge_results");
    add(same_bytes(a.outputs.id, b.outputs.id), "out_id");
    add(same_bytes(a.outputs.species, b.outputs.species), "out_species");
    add(same_bytes(a.outputs.flags, b.outputs.flags), "out_flags");
    for (int k = 0; k < 3; ++k) {
        add(same_bytes(a.outputs.x[k], b.outputs.x[k]), "out_x");
        add(same_bytes(a.outputs.u[k], b.outputs.u[k]), "out_u");
        add(same_bytes(a.reloc.start[k], b.reloc.start[k]), "reloc_start");
        add(same_bytes(a.reloc.D[k], b.reloc.D[k]), "reloc_D");
    }
    add(same_bytes(a.outputs.w, b.outputs.w), "out_w");
    add(same_bytes(a.removed, b.removed), "removed");
    add(same_bytes(a.reloc.id, b.reloc.id), "reloc_id");
    add(same_bytes(a.reloc.group, b.reloc.group), "reloc_group");
    add(same_bytes(a.reloc.charge, b.reloc.charge), "reloc_charge");
    bool led = a.ledger.size() == b.ledger.size();
    for (std::size_t s = 0; led && s < a.ledger.size(); ++s) {
        const auto &x = a.ledger[s], &y = b.ledger[s];
        led = x.w_in == y.w_in && x.w_out == y.w_out && x.e_in == y.e_in && x.e_out == y.e_out && x.e_res == y.e_res &&
              std::memcmp(&x.e_res_rn, &y.e_res_rn, 8) == 0 && x.e_res_expansion == y.e_res_expansion;
        for (int k = 0; k < 3 && led; ++k)
            led = x.p_in[k] == y.p_in[k] && x.p_out[k] == y.p_out[k] && x.p_res[k] == y.p_res[k] &&
                  std::memcmp(&x.p_res_rn[k], &y.p_res_rn[k], 8) == 0 && x.p_res_expansion[k] == y.p_res_expansion[k];
    }
    add(led, "ledger");
    add(std::memcmp(&a.counters, &b.counters, sizeof(Counters)) == 0, "counters");
    add(a.next_id_after == b.next_id_after, "next_id");
    add(encode(a) == encode(b), "encoding");
    add(a.digest_sha256 == b.digest_sha256, "digest");
    return d;
}
std::vector<unsigned> list(const char* s) {
    std::vector<unsigned> v;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, ',')) v.push_back(unsigned(std::stoul(t)));
    return v;
}
std::string timing_json(const GpuTimings& t, std::size_t n) {
    const double per = 1e6 / double(n);
    std::ostringstream o;
    o << "{\"upload_ms\": " << num(t.upload_ms) << ", \"sort_count_ms\": " << num(t.sort_ms) << ", \"form_kinematics_ms\": " << num(t.kinematics_ms)
      << ", \"emit_ledger_ms\": " << num(t.emit_ledger_ms) << ", \"download_ms\": " << num(t.download_ms)
      << ", \"device_compute_ms\": " << num(t.device_compute_ms()) << ", \"device_compute_ms_per_1e6\": " << num(t.device_compute_ms() * per)
      << ", \"with_transfers_ms_per_1e6\": " << num((t.device_compute_ms() + t.upload_ms + t.download_ms) * per)
      << ", \"host_wall_ms\": " << num(t.host_wall_ms) << ", \"host_wall_ms_per_1e6\": " << num(t.host_wall_ms * per)
      << ", \"device_bytes\": " << t.device_bytes << ", \"pinned_bytes\": " << t.pinned_bytes << "}";
    return o.str();
}
}  // namespace

int main(int argc, char** argv) {
    std::string nvrtc, uuid, out_path;
    std::size_t n = 1000000;
    unsigned repeat = 3, corner = 2000;
    bool check_each_launch = false;
    std::vector<unsigned> meshes{45, 90, 180, 360};
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string a = argv[i];
        if (a == "--nvrtc") nvrtc = argv[i + 1];
        else if (a == "--uuid") uuid = argv[i + 1];
        else if (a == "--n") n = std::stoull(argv[i + 1]);
        else if (a == "--meshes") meshes = list(argv[i + 1]);
        else if (a == "--repeat") repeat = unsigned(std::stoul(argv[i + 1]));
        else if (a == "--corner") corner = unsigned(std::stoul(argv[i + 1]));
        else if (a == "--check-each-launch") check_each_launch = std::stoul(argv[i + 1]) != 0;
        else if (a == "--out") out_path = argv[i + 1];
        else { std::cerr << "unknown option " << a << std::endl; return 2; }
    }
    std::ostringstream J;
    J << "{\n \"schema\": \"merge-split-gpu-tests-v1\", \"carriers_per_species\": " << n << ", \"check_each_launch\": " << (check_each_launch ? "true" : "false");
    const auto start = Clock::now();
    int rc = 0;
    {
        hall_gpu::SessionOptions so;
        so.require_uuid = uuid;
        so.nvrtc_library = nvrtc;
        so.check_each_launch = check_each_launch;   // synchronize after every launch (fault localisation)
        hall_gpu::Session session(so);
        try {
            GpuTransaction gpu(session);
            J << ",\n \"device\": " << hall_gpu::device_json(session.device()) << ",\n \"compile\": " << hall_gpu::compile_json(gpu.module().record);
            std::cout << "compiled: " << gpu.module().record.arch_option << " " << gpu.module().record.image_kind << " "
                      << gpu.module().record.compile_ms << " ms" << std::endl;
            auto run_case = [&](const std::string& name, const MsParams& P, const Carriers& C) {
                const auto t0 = Clock::now();
                const auto ref = transact(P, C);
                const double cpu_s = std::chrono::duration<double>(Clock::now() - t0).count();
                std::ostringstream o;
                o << "{\"case\": " << q(name) << ", \"carriers\": " << C.size() << ", \"groups\": " << ref.counters.groups
                  << ", \"accepted_merge\": " << ref.counters.accepted_merge << ", \"accepted_split\": " << ref.counters.accepted_split
                  << ", \"skipped\": " << (ref.counters.groups - ref.counters.accepted_merge - ref.counters.accepted_split)
                  << ", \"outputs\": " << ref.counters.outputs << ", \"cpu_digest\": " << q(ref.digest_sha256)
                  << ", \"cpu_s\": " << num(cpu_s) << ", \"cpu_s_per_1e6\": " << num(cpu_s * 1e6 / double(C.size())) << ", \"runs\": [";
                std::string first_digest;
                bool all_same = true;
                for (unsigned r = 0; r < repeat; ++r) {
                    GpuTimings t;
                    const auto R = gpu.run(P, C, &t);
                    const std::string diff = differences(ref, R);
                    all_same = all_same && diff.empty();
                    if (r == 0) first_digest = R.digest_sha256;
                    check(diff.empty(), name + " CPU/GPU differ in: " + diff);
                    check(R.digest_sha256 == first_digest, name + " device repeat-run digest differs");
                    o << (r ? ", " : "") << "{\"gpu_digest\": " << q(R.digest_sha256) << ", \"identical_to_cpu\": " << (diff.empty() ? "true" : "false")
                      << ", \"differences\": " << q(diff) << ", \"timing\": " << timing_json(t, C.size()) << "}";
                    std::cout << name << " run " << r << ": identical " << diff.empty() << " device " << t.device_compute_ms() << " ms, wall "
                              << t.host_wall_ms << " ms" << std::endl;
                }
                o << "], \"all_identical\": " << (all_same ? "true" : "false") << "}";
                return o.str();
            };
            std::string cases = "[";
            {
                const auto C = fixtures::make_corner_set(0xC0C0, corner);
                cases += run_case("corner-groups", fixtures::corner_policy(0xBADA55ull), C);
            }
            {
                const auto C = fixtures::make_set(0x1234 + 45, n);
                MsParams P = fixtures::policy(45, 0x9e3779b97f4a7c15ull ^ 45ull, 4242, 2 * n + 1);
                P.sub[0] = 2;
                P.sub[1] = 2;
                P.sub[2] = 8;
                cases += ", " + run_case("mesh-45-sub-2-2-8", P, C);
            }
            for (unsigned mesh : meshes) {
                const auto C = fixtures::make_set(0x1234 + std::uint64_t(mesh), n);
                const MsParams P = fixtures::policy(int(mesh), 0x9e3779b97f4a7c15ull ^ std::uint64_t(mesh), 4242, 2 * n + 1);
                cases += ", " + run_case("mesh-" + std::to_string(mesh), P, C);
            }
            J << ",\n \"cases\": " << cases << "]";
            J << ",\n \"memory\": " << hall_gpu::memory_json(session.memory());
        } catch (const std::exception& ex) {
            ++failures;
            std::cout << "FAIL exception: " << ex.what() << std::endl;
            J << ",\n \"exception\": " << q(ex.what());
        }
        const auto td = session.close();
        J << ",\n \"teardown\": " << hall_gpu::teardown_json(td);
        if (td.quarantined) rc = 3;
    }
    J << ",\n \"failures\": " << failures << ", \"wall_s\": " << num(std::chrono::duration<double>(Clock::now() - start).count()) << "\n}\n";
    if (!out_path.empty()) std::ofstream(out_path) << J.str();
    std::cout << (failures ? "FAILED " : "PASSED ") << failures << " failures" << std::endl;
    return rc ? rc : (failures ? 1 : 0);
}
