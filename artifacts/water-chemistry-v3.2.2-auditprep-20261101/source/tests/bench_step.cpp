// bench_step.cpp -- CPU sampling cost of the chemistry collision stage on fixed inputs.
//
//   bench_step gen BLOB N CELLS SEED OUT.bin      write a carrier population (tests/test_common.hpp) to OUT.bin
//   bench_step run BLOB IN.bin SEED STEP DT REPS  run cpu_step REPS times on IN.bin; prints the minimum wall time,
//                                                 real events and candidates (records are identical across reps)
// The run mode uses only chem_host.hpp, so the same source compiled against the v3.1 headers measures the v3.1 code on
// the same inputs (CarrierIn and CellTargets layouts are unchanged between v3.1 and v3.2).
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "chem_host.hpp"
#ifdef BENCH_GEN
#include "../tests/test_common.hpp"
using namespace pfc_test;
#endif

using namespace pfc;

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    const std::string mode = argv[1];
#ifdef BENCH_GEN
    if (mode == "gen" && argc == 7) {
        const Blob b = load_blob(argv[2]);
        const Population P = make_population(b, std::stoull(argv[3]), std::stoull(argv[4]), std::stoull(argv[5]));
        std::ofstream f(argv[6], std::ios::binary);
        const std::uint64_t nc = P.carriers.size(), ne = P.cells.size();
        f.write(reinterpret_cast<const char*>(&nc), 8);
        f.write(reinterpret_cast<const char*>(&ne), 8);
        f.write(reinterpret_cast<const char*>(P.carriers.data()), std::streamsize(nc * sizeof(CarrierIn)));
        f.write(reinterpret_cast<const char*>(P.cells.data()), std::streamsize(ne * sizeof(CellTargets)));
        std::cout << "wrote " << nc << " carriers, " << ne << " cells\n";
        return 0;
    }
#endif
    if (mode == "run" && argc == 8) {
        const Blob b = load_blob(argv[2]);
        std::ifstream f(argv[3], std::ios::binary);
        std::uint64_t nc = 0, ne = 0;
        f.read(reinterpret_cast<char*>(&nc), 8);
        f.read(reinterpret_cast<char*>(&ne), 8);
        std::vector<CarrierIn> cs(nc);
        std::vector<CellTargets> cells(ne);
        f.read(reinterpret_cast<char*>(cs.data()), std::streamsize(nc * sizeof(CarrierIn)));
        f.read(reinterpret_cast<char*>(cells.data()), std::streamsize(ne * sizeof(CellTargets)));
        StepConfig cfg{};
        cfg.seed = std::stoull(argv[4]);
        cfg.step = std::stoull(argv[5]);
        cfg.dt_s = std::stod(argv[6]);
        cfg.trunc_sigmas = 7.0;
        cfg.purpose = PFC_PURPOSE_MCC;
        cfg.max_candidates = 4096;
        cfg.max_events = 64;
        const int reps = std::stoi(argv[7]);
        double best = 1e300;
        std::uint64_t events = 0, cands = 0;
        for (int r = 0; r < reps; ++r) {
            StepResult res;
            const auto t0 = std::chrono::steady_clock::now();
            cpu_step(b, cfg, cs, cells, res);
            const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            if (ms < best) best = ms;
            events = 0;
            cands = 0;
            for (const auto& o : res.outs) {
                events += o.n_events;
                cands += o.n_candidates;
            }
        }
        std::printf("{\"blob\": \"%s\", \"carriers\": %llu, \"events\": %llu, \"candidates\": %llu, \"min_ms\": %.6g, \"us_per_event\": %.6g, \"ns_per_candidate\": %.6g}\n",
                    b.sha256.c_str(), (unsigned long long)nc, (unsigned long long)events, (unsigned long long)cands, best,
                    events ? 1e3 * best / double(events) : 0.0, cands ? 1e6 * best / double(cands) : 0.0);
        return 0;
    }
    std::cerr << "usage: bench_step gen BLOB N CELLS SEED OUT.bin | run BLOB IN.bin SEED STEP DT REPS\n";
    return 2;
}
