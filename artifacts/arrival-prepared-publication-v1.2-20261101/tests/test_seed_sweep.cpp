// Seed sweep of the v1.1 neutral count decision at a fixed expectation:
// n = SR(D / w0) per step over STEPS consecutive steps for SEEDS run seeds; prints the
// standardized deviation z = (sum n - STEPS x) / sqrt(STEPS f (1 - f)) per seed and the
// mean and variance of z over the seeds (expected 0 and 1).
//   test_seed_sweep D W0 STEPS SEEDS FIRST_SEED
#include "pf_arrival_v11/arrival_owner.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
namespace pa = pf_arrival_v11;
int main(int argc, char** argv) {
    if (argc != 6) { std::fprintf(stderr, "usage: test_seed_sweep D W0 STEPS SEEDS FIRST_SEED\n"); return 2; }
    const double D = std::atof(argv[1]), w0 = std::atof(argv[2]);
    const long steps = std::atol(argv[3]), seeds = std::atol(argv[4]);
    const unsigned long long first = std::strtoull(argv[5], nullptr, 10);
    const double x = D / w0, f = x - std::floor(x);
    double sz = 0, szz = 0;
    std::printf("{\"D\": %.17g, \"w0\": %.17g, \"x\": %.17g, \"steps\": %ld, \"seeds\": [", D, w0, x, steps);
    for (long s = 0; s < seeds; ++s) {
        pa::Configuration cfg; cfg.run_seed = first + s; cfg.streams = {{1, pa::Kind::neutral, 0, "H2O", w0, 4096}};
        const pa::Policy P(cfg);
        long double sum = 0;
        for (long st = 0; st < steps; ++st) {
            pa::NeutralRelease r; r.stream_id = 1; r.debit = pa::Count(D); r.log_debit = std::log(D); r.inventory = pa::Count(1e20);
            r.begin_s = 0; r.end_s = 1e-14;
            sum += pa::decide(P, std::uint64_t(st), r).births;
        }
        const double z = double((sum - (long double)steps * x) / std::sqrt(double(steps) * f * (1 - f)));
        sz += z; szz += z * z;
        std::printf("%s{\"seed\": %llu, \"n\": %.0Lf, \"z\": %.4f}", s ? ", " : "", first + s, sum, z);
    }
    const double m = sz / seeds, v = szz / seeds - m * m;
    std::printf("], \"mean_z\": %.4f, \"var_z\": %.4f}\n", m, v);
    return 0;
}
