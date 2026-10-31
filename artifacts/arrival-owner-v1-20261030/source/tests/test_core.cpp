// Core qualification of the arrival owner v1 (no source/native libraries):
// Philox known-answer vectors, stream layout, exact conversions, exact Bernoulli,
// stochastic-rounding unbiasedness, exact debit partition sweeps, order/partition
// independence, IDs, ledger, checkpoint and refusals.
// usage: test_core OUTPUT_DIRECTORY
#include "pf_arrival/arrival_owner.hpp"
#include <openssl/bn.h>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <set>
#include <iostream>
#include <random>

namespace fs = std::filesystem;
using namespace pf_arrival;

namespace {
long checks = 0;
std::vector<std::string> notes;   // JSON fragments for the results file
void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) throw std::runtime_error("check failed: " + what);
}
template <class F>
void refuses(F&& f, RefusalKind kind, const std::string& what) {
    bool hit = false;
    try {
        f();
    } catch (const Refusal& r) {
        hit = r.kind == kind;
        if (!hit) std::cerr << "unexpected refusal kind " << kind_name(r.kind) << ": " << r.what() << "\n";
    }
    check(hit, what + " (expected refusal kind " + kind_name(kind) + ")");
}
std::string hx(double x) { return detail::hexfloat(x); }
std::string u64hex(std::uint64_t v) {
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", (unsigned long long)v);
    return b;
}
std::string jstr(const std::string& s) {
    std::string r = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') r += '\\';
        r += c;
    }
    return r + "\"";
}
std::uint64_t splitmix(std::uint64_t& s) {
    std::uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

// ------------------------------------------------------------ 1. known answers
void known_answers() {
    // Random123 kat_vectors, philox4x32 R=10: counter[4] key[2] -> output[4].
    struct V { pfa_u32 c[4], k[2], out[4]; };
    const V vs[] = {
        {{0x00000000, 0x00000000, 0x00000000, 0x00000000}, {0x00000000, 0x00000000}, {0x6627e8d5, 0xe169c58d, 0xbc57ac4c, 0x9b00dbd8}},
        {{0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff}, {0xffffffff, 0xffffffff}, {0x408f276d, 0x41c83b0e, 0xa20bc7c6, 0x6d5451fd}},
        {{0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344}, {0xa4093822, 0x299f31d0}, {0xd16cfe09, 0x94fdcceb, 0x5001e420, 0x24126ea1}},
    };
    for (const auto& v : vs) {
        pfa::Block c{{v.c[0], v.c[1], v.c[2], v.c[3]}};
        pfa::Key2 k{{v.k[0], v.k[1]}};
        const auto o = pfa::philox4x32_10(c, k);
        for (int i = 0; i < 4; ++i) check(o.v[i] == v.out[i], "Philox4x32-10 known-answer vector word " + std::to_string(i));
    }
    notes.push_back("\"known_answer_vectors\":3");
}

// ------------------------------------------------------------ 2. stream layout
void stream_layout() {
    // The counter/key map is injective: each field moves only its own bits.
    pfa::StreamKey a{0x0123456789abcdefull, 0xA5A5A5A5A5ull, 0xBEEF, 0xCAFEBABE, 0x7F};
    auto ctr = [](const pfa::StreamKey& s, pfa_u32 b) {
        std::array<pfa_u32, 6> r{b, s.candidate, (s.stream & 0xFFFFu) | ((s.purpose & 0xFFu) << 16) | (pfa_u32((s.step >> 32) & 0xFF) << 24),
                                 pfa_u32(s.step), pfa_u32(s.seed), pfa_u32(s.seed >> 32)};
        return r;
    };
    const auto base = ctr(a, 3);
    const auto blk = pfa::stream_block(a, 3);
    pfa::Block c{{base[0], base[1], base[2], base[3]}};
    pfa::Key2 k{{base[4], base[5]}};
    const auto direct = pfa::philox4x32_10(c, k);
    for (int i = 0; i < 4; ++i) check(blk.v[i] == direct.v[i], "stream block equals Philox of the documented counter/key");
    check(base[2] == (0xBEEFu | (0x7Fu << 16) | (0xA5u << 24)) && base[3] == 0xA5A5A5A5u, "counter word layout");
    // Words: R(2b), R(2b+1).
    check(pfa::stream_word(a, 6) == ((pfa_u64(blk.v[1]) << 32) | blk.v[0]) && pfa::stream_word(a, 7) == ((pfa_u64(blk.v[3]) << 32) | blk.v[2]),
          "word order within a block");
    // Distinct fields give distinct counters at the domain edges.
    std::set<std::array<pfa_u32, 6>> seen;
    const pfa_u64 steps[] = {0, 1, 0xFFFFFFFFull, 0x100000000ull, PFA_STEP_LIMIT - 1};
    const pfa_u32 streams[] = {0, 1, 0xFFFF}, cands[] = {0, 1, 0xFFFFFFFFu}, purposes[] = {1, 2, 3, 255};
    std::size_t total = 0;
    for (auto st : steps)
        for (auto sp : streams)
            for (auto cd : cands)
                for (auto pu : purposes)
                    for (pfa_u32 b : {0u, 1u, 0xFFFFFFFFu}) {
                        seen.insert(ctr({7, st, sp, cd, pu}, b));
                        ++total;
                    }
    check(seen.size() == total, "stream layout injective over the edge grid (" + std::to_string(total) + " tuples)");
    notes.push_back("\"layout_edge_tuples\":" + std::to_string(total));
}

// ------------------------------------------------------------ 3. conversions
void conversions() {
    check(pfa::u01_open(0) == 0x1p-53, "u01_open minimum 2^-53");
    check(pfa::u01_open(~0ull) == 1 - 0x1p-53, "u01_open maximum 1-2^-53");
    check(pfa::u01_half_open(0) == 0 && pfa::u01_half_open(~0ull) == 1 - 0x1p-53, "u01_half_open range [0, 1-2^-53]");
    std::uint64_t s = 1;
    for (int i = 0; i < 200000; ++i) {
        const auto r = splitmix(s);
        const double u = pfa::u01_open(r);
        check(u > 0 && u < 1, "u01_open open interval");
        // Exactness: u * 2^53 is the odd integer 2*(r>>12)+1.
        check(u * 0x1p53 == double(((r >> 12) << 1) | 1), "u01_open exact odd lattice");
        check(pfa::u01_open(~r) == 1 - u, "u01_open symmetric under complement");
        const double h = pfa::u01_half_open(r);
        check(h >= 0 && h < 1 && h * 0x1p53 == double(r >> 11), "u01_half_open exact lattice");
    }
}

// ------------------------------------------------------------ 4. fraction digits and exact Bernoulli
BIGNUM* units_shifted(double p) {   // p * 2^1088 as an integer (p in [0,1))
    auto c = Count(p).copy_units();
    BIGNUM* n = BN_dup(c.get());
    BN_lshift(n, n, 14);
    return n;
}
BIGNUM* digits_value(const std::function<pfa_u64(int)>& digit) {
    BIGNUM* n = BN_new();
    BN_zero(n);
    for (int j = 0; j < PFA_FRACTION_DIGITS; ++j) {
        BN_lshift(n, n, 64);
        BIGNUM* d = BN_new();
        BN_set_word(d, digit(j));
        BN_add(n, n, d);
        BN_free(d);
    }
    return n;
}
void bernoulli_exactness() {
    std::vector<double> ps = {0.5, 0.25, 0.75, 1 - 0x1p-53, 0x1p-1074, 3 * 0x1p-1074, 0x1p-1022 - 0x1p-1074, 0x1p-1022, 1e-300, 1e-20, 0.1, 0.3};
    std::uint64_t s = 99;
    for (int i = 0; i < 3000; ++i) {
        double p;
        do {
            p = std::bit_cast<double>(splitmix(s) & 0x3FEFFFFFFFFFFFFFull);   // [0, 1)
        } while (!(p >= 0 && p < 1));
        ps.push_back(p);
        ps.push_back(std::bit_cast<double>(splitmix(s) & 0x000FFFFFFFFFFFFFull));   // subnormal
    }
    for (double p : ps) {
        BIGNUM* a = units_shifted(p);
        BIGNUM* b = digits_value([&](int j) { return pfa::fraction_digit(p, j); });
        check(BN_cmp(a, b) == 0, "fraction digits represent p exactly: " + hx(p));
        BN_free(a);
        BN_free(b);
    }
    // Digit comparison against the definition U < p, with words chosen to match
    // the first m digits of p (deep paths).
    std::size_t deep = 0;
    for (double p : ps) {
        for (int m : {0, 1, 2, 5, 16, 17}) {
            for (int variant = 0; variant < 3; ++variant) {
                std::array<pfa_u64, PFA_FRACTION_DIGITS> w{};
                for (int j = 0; j < PFA_FRACTION_DIGITS; ++j) w[j] = j < m ? pfa::fraction_digit(p, j) : splitmix(s);
                if (m < PFA_FRACTION_DIGITS) {
                    const auto c = pfa::fraction_digit(p, m);
                    if (variant == 1 && c > 0) w[m] = c - 1;
                    if (variant == 2 && c < ~0ull) w[m] = c + 1;
                }
                auto words = [&](int j) { return w[j]; };
                BIGNUM* U = digits_value([&](int j) { return w[j]; });
                BIGNUM* P = units_shifted(p);
                const int expected = (p != 0 && BN_cmp(U, P) < 0) ? 1 : 0;
                check(pfa::bernoulli_digits(words, p) == expected, "Bernoulli digit comparison equals U < p");
                BN_free(U);
                BN_free(P);
                ++deep;
            }
        }
    }
    notes.push_back("\"bernoulli_definition_cases\":" + std::to_string(deep));
}

// ------------------------------------------------------------ 5. statistics
struct SrStat { double x; std::uint64_t draws; double mean, z; };
std::vector<SrStat> sr_stats;
void stochastic_rounding_statistics() {
    // Unbiasedness: z = (mean - x) / sqrt(f(1-f)/M) with |z| <= 5 for each x.
    // Family-wise false-alarm probability <= 8 * 5.73e-7 = 4.6e-6 (Bonferroni),
    // i.e. confidence >= 99.9995 % (normal approximation; M f (1-f) >= 400).
    const double xs[] = {0.5, 0.3, 0.999, 1e-4, 1.25, 7.875, 4095.4, 0.0625};
    const std::uint64_t M = 1ull << 22;
    for (double x : xs) {
        std::uint64_t sum = 0;
        const pfa_u64 floor_x = pfa_u64(x);
        bool range = true;
        for (std::uint64_t i = 0; i < M; ++i) {
            const auto n = pfa::stochastic_round({0x5eed0001ull, i, 17, 0, PFA_PURPOSE_COUNT}, x);
            range = range && (n == floor_x || n == floor_x + 1);
            sum += n;
        }
        check(range, "SR result is floor(x) or floor(x)+1");
        const double f = x - double(floor_x);
        const double mean = double(sum) / double(M);
        const double z = (mean - x) / std::sqrt(f * (1 - f) / double(M));
        sr_stats.push_back({x, M, mean, z});
        check(std::abs(z) <= 5, "SR unbiased at x=" + hx(x) + " z=" + std::to_string(z));
    }
    // Integers are reproduced exactly.
    for (double x : {0., 1., 2., 4096., 12345.})
        for (std::uint64_t i = 0; i < 1000; ++i)
            check(pfa::stochastic_round({1, i, 2, 0, PFA_PURPOSE_COUNT}, x) == pfa_u64(x), "SR of an integer is exact");
    // u01_open uniformity: 256-bin chi-square, threshold = 1 - 1e-6 quantile
    // (Wilson-Hilferty), and the mean within 5 sigma.
    const std::uint64_t U = 1ull << 22;
    std::vector<std::uint64_t> bins(256, 0);
    long double total = 0;
    for (std::uint64_t i = 0; i < U; ++i) {
        const double u = pfa::u01_open(pfa::stream_word({0x5eed0002ull, i >> 1, 3, 0, PFA_PURPOSE_ELECTRON_MARKS}, pfa_u32(i & 1)));
        bins[std::size_t(u * 256)]++;
        total += u;
    }
    double chi2 = 0;
    const double e = double(U) / 256;
    for (auto b : bins) chi2 += (double(b) - e) * (double(b) - e) / e;
    const double df = 255, zq = 4.753424;
    const double threshold = df * std::pow(1 - 2 / (9 * df) + zq * std::sqrt(2 / (9 * df)), 3);
    const double mean = double(total / U), zmean = (mean - 0.5) / std::sqrt(1.0 / 12 / double(U));
    check(chi2 <= threshold, "u01_open 256-bin chi-square " + std::to_string(chi2) + " <= " + std::to_string(threshold));
    check(std::abs(zmean) <= 5, "u01_open mean within 5 sigma");
    std::ostringstream o;
    o << "\"stochastic_rounding\":[";
    for (std::size_t i = 0; i < sr_stats.size(); ++i)
        o << (i ? "," : "") << "{\"x\":" << jstr(hx(sr_stats[i].x)) << ",\"draws\":" << sr_stats[i].draws << ",\"mean\":" << std::setprecision(17)
          << sr_stats[i].mean << ",\"z\":" << sr_stats[i].z << "}";
    o << "],\"stochastic_rounding_test\":\"|z|<=5 per x, 8 values, family-wise false alarm <= 4.6e-6\"";
    o << ",\"u01_open_chi2\":" << chi2 << ",\"u01_open_chi2_threshold_1e-6\":" << threshold << ",\"u01_open_mean_z\":" << zmean
      << ",\"u01_open_draws\":" << U;
    notes.push_back(o.str());
}

// ------------------------------------------------------------ 6. policy, decisions, exact debits
std::shared_ptr<const Policy> policy(std::uint64_t seed = 20261030) {
    Configuration c;
    c.run_seed = seed;
    c.streams = {{0, Kind::electron, 0, "e-", 1e5, 4096},     {1, Kind::electron, 0, "e-", 3e4, 4096},
                 {2, Kind::electron, 1, "e-", 2.5e3, 4096},   {7, Kind::neutral, 0, "H2O", 1e6, 4096},
                 {9, Kind::neutral, 1, "H2O", 7.3e4, 4096},   {300, Kind::neutral, 2, "H2O", 1e-300, 4096},
                 {65535, Kind::electron, 3, "e-", 1e9, 64}};
    return std::make_shared<const Policy>(c);
}
ElectronExposure electron(std::uint32_t id, double N, double begin = 7.4978805549199997e-10, double dt = 1e-14) {
    return {id, N, N > 0 ? std::log(N) : -INFINITY, N > 0, begin, dt, begin + dt};
}
NeutralRelease neutral(std::uint32_t id, double D, Count inventory = Count(1e300), double begin = 7.4978805549199997e-10, double dt = 1e-14) {
    return {id, D > 0 ? Count(D) : Count(), D > 0 ? std::log(D) : -INFINITY, inventory, begin, begin + dt};
}

void configuration_refusals() {
    auto base = [] {
        Configuration c;
        c.run_seed = 1;
        c.streams = {{0, Kind::electron, 0, "e-", 1e5, 4096}};
        return c;
    };
    auto with = [&](auto&& f) {
        auto c = base();
        f(c);
        return c;
    };
    refuses([&] { Policy p(with([](auto& c) { c.streams.clear(); })); }, RefusalKind::configuration, "empty stream list");
    refuses([&] { Policy p(with([](auto& c) { c.streams[0].stream_id = 65536; })); }, RefusalKind::configuration, "stream id >= 2^16");
    refuses([&] { Policy p(with([](auto& c) { c.streams.push_back(c.streams[0]); })); }, RefusalKind::configuration, "duplicate stream id");
    refuses([&] { Policy p(with([](auto& c) { c.streams[0].target_weight = 0; })); }, RefusalKind::configuration, "zero target weight");
    refuses([&] { Policy p(with([](auto& c) { c.streams[0].target_weight = 1e-310; })); }, RefusalKind::configuration, "subnormal target weight");
    refuses([&] { Policy p(with([](auto& c) { c.streams[0].target_weight = INFINITY; })); }, RefusalKind::configuration, "infinite target weight");
    refuses([&] { Policy p(with([](auto& c) { c.streams[0].maximum_births = 4097; })); }, RefusalKind::configuration, "per-stream cap above 4096");
    refuses([&] { Policy p(with([](auto& c) { c.streams[0].maximum_births = 0; })); }, RefusalKind::configuration, "zero per-stream cap");
    refuses([&] { Policy p(with([](auto& c) { c.streams[0].patch_index = 256; })); }, RefusalKind::configuration, "patch index above source cap");
    refuses([&] { Policy p(with([](auto& c) { c.streams[0].species_state.clear(); })); }, RefusalKind::configuration, "empty species state");
    refuses([&] {
        Policy p(with([](auto& c) {
            c.streams = {{0, Kind::neutral, 4, "H2O", 1, 8}, {1, Kind::neutral, 4, "H2O", 2, 8}};
        }));
    }, RefusalKind::configuration, "two neutral streams on one source patch");
    refuses([&] { Policy p(with([](auto& c) { c.maximum_neutral_leaves_per_step = 5000; })); }, RefusalKind::configuration, "step cap above 4096");
    // Configuration identity binds every field.
    const auto a = Policy(base()).configuration_sha256();
    check(Policy(with([](auto& c) { c.run_seed = 2; })).configuration_sha256() != a, "seed bound into identity");
    check(Policy(with([](auto& c) { c.streams[0].target_weight = std::nextafter(1e5, 2e5); })).configuration_sha256() != a, "w0 bits bound");
    check(Policy(with([](auto& c) { c.streams[0].species_state = "e"; })).configuration_sha256() != a, "species bound");
    check(Policy(with([](auto& c) { c.streams[0].maximum_births = 4095; })).configuration_sha256() != a, "cap bound");
    check(Policy(base()).configuration_sha256() == a, "identity deterministic");
}

std::size_t debit_cases = 0, merged_cases = 0, floored_cases = 0;
std::vector<std::string> oracle_partitions;
void neutral_partition_case(double D, double w0, std::uint64_t step) {
    if (!detail::positive_normal(w0)) return;
    Configuration c;
    c.run_seed = 77;
    c.streams = {{5, Kind::neutral, 0, "H2O", w0, 4096}};
    Policy P(c);
    const double x = D / w0;
    if (!(x <= 4096)) {
        refuses([&] { (void)decide(P, step, neutral(5, D, Count(DBL_MAX))); }, RefusalKind::macroparticle_cap, "D/w0 above cap refused");
        return;
    }
    const auto d = decide(P, step, neutral(5, D, Count(DBL_MAX)));
    check(d.debit == Count(D) && d.births >= 1 && d.births <= 4096, "neutral decision bounds");
    check(d.target == std::max<std::uint64_t>(1, d.rounded) && d.births == std::min(d.target, d.quanta) && d.merged == (d.target > d.quanta),
          "neutral merge/floor rule");
    check(std::ldexp(double(d.quanta), int(d.exponent)) == D && d.quanta < (1ull << 53) && d.exponent >= -1074, "d = M * 2^E");
    Count total;
    double wmin = INFINITY, wmax = 0;
    for (std::uint32_t k = 0; k < d.births; ++k) {
        const auto b = evaluate_birth(P, d, k);
        check(b.weight > 0 && std::isfinite(b.weight), "positive finite leaf");
        check(Count(b.weight).positive_carriers().size() == 1, "one native carrier per leaf");
        total = total.credited(b.weight);
        wmin = std::min(wmin, b.weight);
        wmax = std::max(wmax, b.weight);
        for (double m : b.marks) check(m > 0 && m < 1, "open neutral marks");
    }
    check(total == Count(D), "leaves sum exactly to the source debit: " + hx(D));
    check(d.births == 1 || std::nextafter(wmin, INFINITY) >= wmax || (wmax - wmin) == std::ldexp(1.0, int(d.exponent)),
          "leaf weights differ by at most one quantum");
    ++debit_cases;
    merged_cases += d.merged;
    floored_cases += d.floored;
    if (oracle_partitions.size() < 400) {
        std::ostringstream o;
        o << "{\"D\":" << jstr(hx(D)) << ",\"w0\":" << jstr(hx(w0)) << ",\"x\":" << jstr(hx(d.ratio)) << ",\"rounded\":" << d.rounded
          << ",\"births\":" << d.births << ",\"step\":" << step << ",\"weights\":[";
        for (std::uint32_t k = 0; k < d.births && k < 64; ++k) o << (k ? "," : "") << jstr(hx(evaluate_birth(P, d, k).weight));
        o << "],\"weights_listed\":" << std::min<std::uint64_t>(d.births, 64) << "}";
        if (d.births <= 64) oracle_partitions.push_back(o.str());
    }
}
std::size_t electron_cases = 0;
std::vector<std::string> oracle_electrons;
void electron_case(double N, double w0, std::uint64_t step) {
    if (!detail::positive_normal(w0)) return;
    Configuration c;
    c.run_seed = 78;
    c.streams = {{4, Kind::electron, 0, "e-", w0, 4096}};
    Policy P(c);
    const double x = N / w0;
    if (!(x <= 4096)) {
        refuses([&] { (void)decide(P, step, electron(4, N)); }, RefusalKind::macroparticle_cap, "N/w0 above cap refused");
        return;
    }
    const auto d = decide(P, step, electron(4, N));
    check(d.ratio == x && (d.births == pfa_u64(x) || d.births == pfa_u64(x) + 1), "electron count is SR(N/w0)");
    Count total;
    for (std::uint32_t k = 0; k < d.births; ++k) {
        const auto b = evaluate_birth(P, d, k);
        check(b.weight == w0, "electron weight is exactly w0");
        total = total.credited(b.weight);
        check(b.time_s >= d.begin_s && b.time_s <= d.begin_s + d.duration_s && b.remaining_s == d.end_s - b.time_s, "electron birth time");
        for (int i = 0; i < 4; ++i) check(b.marks[i] > 0 && b.marks[i] < 1, "open electron marks");
    }
    check(total == d.debit && total == detail::credit_n(w0, d.births), "electron debit = n*w0 exactly");
    ++electron_cases;
    if (oracle_electrons.size() < 300) {
        std::ostringstream o;
        o << "{\"N\":" << jstr(hx(N)) << ",\"w0\":" << jstr(hx(w0)) << ",\"step\":" << step << ",\"births\":" << d.births << "}";
        oracle_electrons.push_back(o.str());
    }
}
void debit_sweeps() {
    const std::vector<double> Ds = {0x1p-1074, 3 * 0x1p-1074, 0x1p-1022 - 0x1p-1074, 0x1p-1022, 1e-300, 1e-200, 1e-100, 1e-20, 0.1, 1,
                                    1 + 0x1p-52, 0.125, 7, 1e15, 0x1p53, 1e100, 1e300, DBL_MAX};
    const std::vector<double> xs = {1e-9, 0.3, 0.999999, 1, 2.5, 3, 100.5, 1023.75, 4095.9, 4096, 4097, 1e6};
    std::uint64_t step = 0;
    for (double D : Ds)
        for (double x : xs) neutral_partition_case(D, D / x, step++);
    // Boundary of the merge-rule proof: the smallest admissible w0 (2^-1022) with
    // subnormal debits gives x < 1, so target = 1 <= M; normal debits have M >= 2^52.
    for (double D : {0x1p-1074, 3 * 0x1p-1074, 7 * 0x1p-1074, 1000 * 0x1p-1074, 0x1p-1022 - 0x1p-1074}) neutral_partition_case(D, DBL_MIN, step++);
    std::uint64_t s = 4242;
    for (int i = 0; i < 4000; ++i) {
        double D = std::bit_cast<double>(splitmix(s) & 0x7FEFFFFFFFFFFFFFull);
        if (!(D > 0 && std::isfinite(D))) continue;
        const double x = std::ldexp(double(splitmix(s) >> 11) * 0x1p-53, int(splitmix(s) % 13));   // [0, 4096)
        if (x > 0) neutral_partition_case(D, D / x, step++);
    }
    const std::vector<double> Ns = {1e-300, 1e-100, 1e-12, 0.5, 1, 1e3, 1e9, 1e15, 1e20, 1e100, 1e300};
    for (double N : Ns)
        for (double x : xs) electron_case(N, N / x, step++);
    for (int i = 0; i < 2000; ++i) {
        const double N = std::ldexp(1 + double(splitmix(s) >> 12) * 0x1p-52, int(splitmix(s) % 1600) - 800);
        const double x = std::ldexp(double(splitmix(s) >> 11) * 0x1p-53, int(splitmix(s) % 13));
        if (x > 0) electron_case(N, N / x, step++);
    }
    check(merged_cases == 0, "merge rule never triggered (unreachable under v1 constraints, API-v1.txt)");
    check(floored_cases > 0, "sweep exercised the floor rule");
    notes.push_back("\"neutral_debit_cases\":" + std::to_string(debit_cases) + ",\"neutral_merged_cases\":" + std::to_string(merged_cases) +
                    ",\"neutral_floored_cases\":" + std::to_string(floored_cases) + ",\"electron_debit_cases\":" + std::to_string(electron_cases));
}

void decision_refusals() {
    auto P = policy();
    refuses([&] { (void)decide(*P, 0, electron(0, NAN)); }, RefusalKind::expected_count, "NaN expected count");
    refuses([&] { (void)decide(*P, 0, electron(0, -1)); }, RefusalKind::expected_count, "negative expected count");
    refuses([&] {
        auto e = electron(0, 0);
        e.has_log = true;
        e.log_expected_number = -800;
        (void)decide(*P, 0, e);
    }, RefusalKind::expected_count, "finite log with zero scalar");
    refuses([&] {
        auto e = electron(0, 1);
        e.log_expected_number = NAN;
        (void)decide(*P, 0, e);
    }, RefusalKind::expected_count, "log flagged but NaN");
    refuses([&] { (void)decide(*P, 0, electron(0, 1e5 * 4097)); }, RefusalKind::macroparticle_cap, "electron N/w0 above cap");
    refuses([&] { (void)decide(*P, 0, electron(65535, 1e9 * 65)); }, RefusalKind::macroparticle_cap, "electron N/w0 above a stream cap of 64");
    refuses([&] { (void)decide(*P, 0, electron(3, 1)); }, RefusalKind::stream_input, "unknown stream");
    refuses([&] { (void)decide(*P, 0, electron(7, 1)); }, RefusalKind::stream_input, "electron exposure on a neutral stream");
    refuses([&] { (void)decide(*P, 0, neutral(0, 1)); }, RefusalKind::stream_input, "neutral release on an electron stream");
    refuses([&] { (void)decide(*P, PFA_STEP_LIMIT, electron(0, 1)); }, RefusalKind::step_index, "step >= 2^40");
    refuses([&] {
        auto e = electron(0, 1);
        e.duration_s = 0;
        (void)decide(*P, 0, e);
    }, RefusalKind::clock, "zero exposure duration");
    refuses([&] {
        auto e = electron(0, 1);
        e.end_s = e.begin_s;
        (void)decide(*P, 0, e);
    }, RefusalKind::clock, "exposure outside the step");
    refuses([&] {
        auto r = neutral(7, 1);
        r.debit = Count::from_positive_carriers({1.0, 0x1p-80});
        (void)decide(*P, 0, r);
    }, RefusalKind::debit_representation, "multi-carrier debit");
    refuses([&] { (void)decide(*P, 0, neutral(7, 2, Count(1))); }, RefusalKind::inventory, "debit above the exact inventory");
    refuses([&] { (void)decide(*P, 0, neutral(7, 1e6 * 5000, Count(1e300))); }, RefusalKind::macroparticle_cap, "neutral D/w0 above cap");
    // Zero exposure/debit: no births, no draws.
    check(decide(*P, 0, electron(0, 0)).births == 0 && decide(*P, 0, neutral(7, 0)).births == 0, "zero expected count and debit");
    // Exact inventory boundary: D == inventory admitted.
    check(decide(*P, 0, neutral(7, 3e6, Count(3e6))).births >= 1, "debit equal to the exact inventory admitted");
    // Per-step caps (drawn totals).
    Configuration c;
    c.run_seed = 5;
    c.streams = {{0, Kind::electron, 0, "e-", 1, 4096}, {1, Kind::electron, 1, "e-", 1, 4096}};
    c.maximum_electron_births_per_step = 100;
    Policy small(c);
    refuses([&] { (void)plan_step(small, {0, {electron(0, 60), electron(1, 60)}, {}}); }, RefusalKind::step_cap, "per-step electron total cap");
    refuses([&] { (void)plan_step(small, {0, {electron(0, 6), electron(0, 6)}, {}}); }, RefusalKind::stream_input, "duplicate stream in one step");
}

// ------------------------------------------------------------ 7. order and partition independence
StepInput step_input(std::uint64_t step) {
    const double f = 1 + double(step % 5);
    StepInput in;
    in.step = step;
    in.electrons = {electron(0, 1.234e7 * f), electron(1, 3.3e5 / f), electron(2, 0.37 * 2.5e3 * f), electron(65535, 0)};
    in.neutrals = {neutral(7, 7.77e8 * f + 0.123), neutral(9, 1.0e3 * f), neutral(300, (step % 3 + 1) * 0x1p-1074)};
    return in;
}
// A worker evaluates a subset of births with nothing shared but the policy and
// the step input (decisions are recomputed from the key).
std::vector<Birth> worker(const Policy& P, const StepInput& in, std::span<const WorkItem> items) {
    std::map<std::uint32_t, Decision> local;
    std::vector<Birth> out;
    for (const auto& it : items) {
        if (!local.contains(it.stream_id)) {
            for (const auto& e : in.electrons)
                if (e.stream_id == it.stream_id) local.emplace(it.stream_id, decide(P, in.step, e));
            for (const auto& r : in.neutrals)
                if (r.stream_id == it.stream_id) local.emplace(it.stream_id, decide(P, in.step, r));
        }
        out.push_back(evaluate_birth(P, local.at(it.stream_id), it.index));
    }
    return out;
}
void partition_independence() {
    auto P = policy();
    std::mt19937_64 shuffle(12345);
    std::size_t compared = 0;
    for (std::uint64_t step : {0ull, 1ull, 2ull, 17ull, 1000ull, (1ull << 39) + 5}) {
        const auto in = step_input(step);
        const auto serial = plan_step(*P, in);
        // Reverse the input order too.
        auto reversed = in;
        std::reverse(reversed.electrons.begin(), reversed.electrons.end());
        std::reverse(reversed.neutrals.begin(), reversed.neutrals.end());
        check(plan_step(*P, reversed).encode() == serial.encode(), "input order does not change the plan");
        auto items = work_items(decide_all(*P, in));
        for (int chunks : {1, 2, 3, 5, 8}) {
            for (int rep = 0; rep < 3; ++rep) {
                auto shuffled = items;
                std::shuffle(shuffled.begin(), shuffled.end(), shuffle);
                std::vector<Birth> merged;
                // Uneven split boundaries.
                std::vector<std::size_t> cuts{0};
                for (int c = 1; c < chunks; ++c) cuts.push_back(std::size_t(shuffle() % (shuffled.size() + 1)));
                cuts.push_back(shuffled.size());
                std::sort(cuts.begin(), cuts.end());
                std::vector<std::vector<Birth>> parts;
                for (int c = 0; c < chunks; ++c)
                    parts.push_back(worker(*P, in, std::span(shuffled).subspan(cuts[c], cuts[c + 1] - cuts[c])));
                std::shuffle(parts.begin(), parts.end(), shuffle);
                for (auto& p : parts) merged.insert(merged.end(), p.begin(), p.end());
                auto plan = assemble(*P, step, decide_all(*P, in), merged);
                check(plan.encode() == serial.encode() && plan.sha256() == serial.sha256(),
                      "partitioned plan identical: step " + std::to_string(step) + " chunks " + std::to_string(chunks));
                // IDs after staging are identical too.
                Owner o(P, step);
                auto a = o.stage(plan, 1000);
                auto b = o.stage(serial, 1000);
                for (std::size_t k = 0; k < a.births().size(); ++k)
                    check(a.births()[k].native_id == b.births()[k].native_id && a.births()[k].ordinal == b.births()[k].ordinal &&
                              a.births()[k].packet_id() == b.births()[k].packet_id(),
                          "partitioned IDs identical");
                check(a.token() == b.token(), "staged token identical");
                ++compared;
            }
        }
        // Incomplete or duplicated coverage refuses.
        auto births = worker(*P, in, items);
        if (!births.empty()) {
            auto missing = births;
            missing.pop_back();
            refuses([&] { (void)assemble(*P, step, decide_all(*P, in), missing); }, RefusalKind::incomplete, "missing birth refused");
            auto dup = births;
            dup.push_back(births.front());
            refuses([&] { (void)assemble(*P, step, decide_all(*P, in), dup); }, RefusalKind::incomplete, "duplicated birth refused");
            auto altered = births;
            altered.front().weight = std::nextafter(altered.front().weight, INFINITY);
            refuses([&] { (void)assemble(*P, step, decide_all(*P, in), altered); }, RefusalKind::debit_representation,
                    "altered leaf weight refused by the exact sum");
        }
    }
    notes.push_back("\"partition_comparisons\":" + std::to_string(compared) + ",\"partition_chunks\":[1,2,3,5,8]");
}

// ------------------------------------------------------------ 8. owner, IDs, ledger, checkpoint
// Native stand-in: one monotone ID counter shared by genesis carriers, reaction
// products and birth blocks (the native v5 genesis/reaction rule).
struct MockNative {
    std::uint64_t next_id = 1;
    std::set<std::uint64_t> issued;
    void allocate(std::uint64_t n) {
        for (std::uint64_t i = 0; i < n; ++i) check(issued.insert(next_id++).second, "mock native ID unique");
    }
    NativeBirthReceipt births(const Policy& P, const Staged& s) {
        NativeBirthReceipt r;
        r.base = next_id;
        for (const auto& b : s.births()) {
            r.entries.push_back({b.packet_id(), P.stream(b.stream_id).species_state, next_id, b.weight});
            check(issued.insert(next_id++).second, "birth native ID unique");
        }
        r.next_id_after = next_id;
        return r;
    }
};
void bind_neutral_times(Staged& s) {
    for (const auto& d : s.plan().decisions()) {
        if (d.kind != Kind::neutral) continue;
        std::vector<double> t;
        for (const auto& b : s.births())
            if (b.stream_id == d.stream_id) t.push_back(d.begin_s + b.marks[5] * (d.end_s - d.begin_s));
        s.bind_source_times(d.stream_id, t);
    }
}
struct RunResult {
    std::vector<std::string> heads, states;
    std::vector<LedgerRecord> records;
    std::string checkpoint_at_20;
};
RunResult run_steps(const std::shared_ptr<const Policy>& P, std::uint64_t first, std::uint64_t count, int chunks,
                    std::unique_ptr<Owner> start = nullptr, MockNative* native_in = nullptr) {
    MockNative local;
    MockNative& native = native_in ? *native_in : local;
    if (!native_in) native.allocate(42);   // genesis carriers
    auto o = start ? std::move(start) : std::make_unique<Owner>(P, first);
    RunResult rr;
    std::mt19937_64 g(chunks);
    for (std::uint64_t i = 0; i < count; ++i) {
        const auto step = o->next_step();
        const auto in = step_input(step);
        native.allocate(step % 4);   // first-exposure reaction products of this step
        auto items = work_items(decide_all(*P, in));
        std::shuffle(items.begin(), items.end(), g);
        std::vector<Birth> merged;
        const std::size_t per = (items.size() + chunks - 1) / std::max(chunks, 1);
        for (std::size_t a = 0; a < items.size(); a += std::max<std::size_t>(per, 1)) {
            auto part = worker(*P, in, std::span(items).subspan(a, std::min(per, items.size() - a)));
            merged.insert(merged.end(), part.begin(), part.end());
        }
        auto plan = assemble(*P, step, decide_all(*P, in), merged);
        auto staged = o->stage(std::move(plan), native.next_id);
        bind_neutral_times(staged);
        auto receipt = native.births(*P, staged);
        rr.records.push_back(o->commit(std::move(staged), receipt));
        native.allocate((step + 1) % 3);   // second-exposure products
        rr.heads.push_back(o->ledger_head());
        rr.states.push_back(o->state_sha256());
        if (i + 1 == 20) rr.checkpoint_at_20 = o->save();
    }
    return rr;
}
void owner_ids_ledger_checkpoint(const fs::path& out) {
    auto P = policy();
    const auto a = run_steps(P, 100, 40, 1);
    for (int chunks : {2, 3, 5, 8}) {
        const auto b = run_steps(P, 100, 40, chunks);
        check(a.heads == b.heads && a.states == b.states, "multi-step ledger identical for " + std::to_string(chunks) + " workers");
    }
    {   // The 40 ledger step records as a canonical sequence (u64 count, then text records).
        Writer w;
        w.tag("pf-arrival-ledger-export-v1");
        w.u64(a.records.size());
        for (const auto& r : a.records) w.text(r.bytes);
        std::ofstream f(out / "ledger-steps-100-139.bin", std::ios::binary);
        f.write(w.bytes.data(), std::streamsize(w.bytes.size()));
    }
    // Chain audit from the genesis head.
    Owner genesis(P, 100);
    check(chain(genesis.ledger_head(), a.records) == a.heads.back(), "ledger chain recomputes from genesis");
    auto broken = a.records;
    broken[5].bytes[40] ^= 1;
    bool refused = false;
    try {
        (void)chain(genesis.ledger_head(), broken);
    } catch (const EncodingError&) {
        refused = true;
    }
    check(refused, "altered ledger record detected");
    // Checkpoint after 20 steps, restore, continue: identical to the uninterrupted run.
    MockNative native;
    native.allocate(42);
    auto first = run_steps(P, 100, 20, 1, nullptr, &native);
    check(first.checkpoint_at_20.size() < 4096, "checkpoint is small (" + std::to_string(first.checkpoint_at_20.size()) + " bytes)");
    auto restored = Owner::restore(P, first.checkpoint_at_20);
    check(restored->state_sha256() == a.states[19], "restored state equals the accepted state");
    auto rest = run_steps(P, 0, 20, 3, std::move(restored), &native);
    check(rest.heads.back() == a.heads.back() && rest.states.back() == a.states.back(), "restart continues bit-identically");
    {
        std::ofstream f(out / "checkpoint-step-120.bin", std::ios::binary);
        f.write(first.checkpoint_at_20.data(), std::streamsize(first.checkpoint_at_20.size()));
    }
    // Checkpoint refusals.
    auto cp = first.checkpoint_at_20;
    for (std::size_t pos : {std::size_t(10), cp.size() / 2, cp.size() - 3}) {
        auto bad = cp;
        bad[pos] ^= 0x10;
        refuses([&] { (void)Owner::restore(P, bad); }, RefusalKind::checkpoint, "corrupted checkpoint byte " + std::to_string(pos));
    }
    refuses([&] { (void)Owner::restore(P, cp.substr(0, cp.size() - 1)); }, RefusalKind::checkpoint, "truncated checkpoint");
    refuses([&] { (void)Owner::restore(P, cp + "x"); }, RefusalKind::checkpoint, "trailing checkpoint bytes");
    refuses([&] { (void)Owner::restore(policy(7), cp); }, RefusalKind::checkpoint, "foreign configuration checkpoint");
    // Every ID unique across steps and disjoint from reaction products (MockNative checks).
    std::set<std::string> packets;
    std::set<std::uint64_t> ordinals;
    {
        MockNative n2;
        n2.allocate(42);
        auto o = std::make_unique<Owner>(P, 100);
        for (int i = 0; i < 25; ++i) {
            const auto step = o->next_step();
            n2.allocate(step % 4);
            auto s = o->stage(o->plan(step_input(step)), n2.next_id);
            bind_neutral_times(s);
            for (const auto& b : s.births()) {
                check(packets.insert(b.packet_id()).second && ordinals.insert(b.ordinal).second, "packet ID and ordinal unique");
                check(b.packet_id().rfind("arrival-v1/s", 0) == 0, "reserved packet namespace");
            }
            auto r = n2.births(*P, s);
            (void)o->commit(std::move(s), r);
        }
        check(ordinals.size() == o->next_ordinal() && *ordinals.rbegin() + 1 == o->next_ordinal(), "ordinals contiguous from 0");
    }
    notes.push_back("\"checkpoint_bytes\":" + std::to_string(cp.size()) + ",\"ledger_steps\":" + std::to_string(a.records.size()) +
                    ",\"ledger_head\":" + jstr(a.heads.back()) + ",\"unique_birth_ids\":" + std::to_string(packets.size()));
}

void staging_and_commit_refusals() {
    auto P = policy();
    MockNative native;
    native.allocate(10);
    Owner o(P, 5);
    const auto before = o.state_sha256();
    RefusalLog log(*P);
    auto attempt = [&](auto&& f, RefusalKind k, const std::string& what) {
        try {
            f();
            check(false, what + " was not refused");
        } catch (const Refusal& r) {
            check(r.kind == k, what + " refusal kind");
            (void)log.append(o.next_step(), r);
        }
        check(o.state_sha256() == before, what + ": accepted state unchanged");
        ++checks;
    };
    attempt([&] { (void)o.plan(step_input(6)); }, RefusalKind::step_index, "skipped step");
    attempt([&] { (void)o.plan(step_input(4)); }, RefusalKind::step_index, "replayed step");
    attempt([&] { (void)o.stage(o.plan(step_input(5)), 0); }, RefusalKind::id_space, "native base 0");
    attempt([&] { (void)o.stage(o.plan(step_input(5)), ~0ull - 3); }, RefusalKind::id_space, "native ID exhaustion");
    const auto plan = o.plan(step_input(5));
    auto make = [&] {
        auto s = o.stage(plan, native.next_id);
        bind_neutral_times(s);
        return s;
    };
    auto good_receipt = [&](const Staged& s) {
        MockNative copy = native;
        return copy.births(*P, s);
    };
    attempt([&] {
        auto s = o.stage(plan, native.next_id);
        (void)o.commit(std::move(s), good_receipt(s));
    }, RefusalKind::incomplete, "neutral times not bound");
    attempt([&] {
        auto s = make();
        auto r = good_receipt(s);
        r.base += 1;
        (void)o.commit(std::move(s), r);
    }, RefusalKind::native_receipt, "receipt base differs");
    attempt([&] {
        auto s = make();
        auto r = good_receipt(s);
        std::swap(r.entries[0], r.entries[1]);
        (void)o.commit(std::move(s), r);
    }, RefusalKind::native_receipt, "receipt order differs from canonical order");
    attempt([&] {
        auto s = make();
        auto r = good_receipt(s);
        r.entries.back().weight = std::nextafter(r.entries.back().weight, 0.);
        (void)o.commit(std::move(s), r);
    }, RefusalKind::native_receipt, "receipt carrier weight differs");
    attempt([&] {
        auto s = make();
        auto r = good_receipt(s);
        r.entries[2].state = "e";
        (void)o.commit(std::move(s), r);
    }, RefusalKind::native_receipt, "receipt species differs");
    attempt([&] {
        auto s = make();
        auto r = good_receipt(s);
        r.next_id_after += 1;
        (void)o.commit(std::move(s), r);
    }, RefusalKind::native_receipt, "receipt next ID differs");
    attempt([&] {
        auto s = make();
        auto r = good_receipt(s);
        r.entries.pop_back();
        (void)o.commit(std::move(s), r);
    }, RefusalKind::native_receipt, "receipt missing a birth");
    attempt([&] {
        auto s = make();
        std::size_t n7 = 0;
        for (const auto& b : s.births()) n7 += b.stream_id == 7;
        s.bind_source_times(7, std::vector<double>(n7, 0.0));
    }, RefusalKind::clock, "source time outside the step");
    // A staged step from another accepted state is stale.
    auto stale = make();
    {
        auto s = make();
        auto r = native.births(*P, s);
        (void)o.commit(std::move(s), r);
    }
    const auto after = o.state_sha256();
    bool hit = false;
    try {
        (void)o.commit(std::move(stale), good_receipt(stale));
    } catch (const Refusal& r) {
        hit = r.kind == RefusalKind::stale;
    }
    check(hit && o.state_sha256() == after, "stale staged step refused, state unchanged");
    check(log.count() == 12, "refusal log appended every refusal");
    notes.push_back("\"refusal_log_head\":" + jstr(log.head()) + ",\"refusal_log_records\":" + std::to_string(log.count()));
}

// ------------------------------------------------------------ 9. oracle dump
void oracle_dump(const fs::path& out) {
    std::ofstream f(out / "oracle-input.json");
    f << "{\"blocks\":[";
    std::uint64_t s = 31337;
    for (int i = 0; i < 2000; ++i) {
        pfa::StreamKey k{splitmix(s), splitmix(s) & (PFA_STEP_LIMIT - 1), pfa_u32(splitmix(s) & 0xFFFF), pfa_u32(splitmix(s)),
                         pfa_u32(1 + splitmix(s) % 255)};
        if (i == 0) k = {0, 0, 0, 0, 1};
        if (i == 1) k = {~0ull, PFA_STEP_LIMIT - 1, 0xFFFF, 0xFFFFFFFFu, 255};
        const pfa_u32 b = i < 2 ? (i ? 0xFFFFFFFFu : 0) : pfa_u32(splitmix(s));
        const auto o = pfa::stream_block(k, b);
        double em[4], nm[6];
        pfa::electron_marks(k, em);
        pfa::neutral_marks(k, nm);
        const double x = std::ldexp(double(splitmix(s) >> 11) * 0x1p-53, int(splitmix(s) % 13));
        f << (i ? "," : "") << "{\"seed\":\"" << u64hex(k.seed) << "\",\"step\":" << k.step << ",\"stream\":" << k.stream
          << ",\"candidate\":" << k.candidate << ",\"purpose\":" << k.purpose << ",\"block\":" << b << ",\"out\":[" << o.v[0] << "," << o.v[1]
          << "," << o.v[2] << "," << o.v[3] << "],\"electron_marks\":[";
        for (int j = 0; j < 4; ++j) f << (j ? "," : "") << jstr(hx(em[j]));
        f << "],\"neutral_marks\":[";
        for (int j = 0; j < 6; ++j) f << (j ? "," : "") << jstr(hx(nm[j]));
        f << "],\"x\":" << jstr(hx(x)) << ",\"sr\":" << pfa::stochastic_round(k, x) << "}";
    }
    f << "],\"partitions\":[";
    for (std::size_t i = 0; i < oracle_partitions.size(); ++i) f << (i ? "," : "") << oracle_partitions[i];
    f << "],\"electrons\":[";
    for (std::size_t i = 0; i < oracle_electrons.size(); ++i) f << (i ? "," : "") << oracle_electrons[i];
    f << "],\"seeds\":{\"neutral\":77,\"electron\":78,\"neutral_stream\":5,\"electron_stream\":4}}\n";
}
}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: test_core OUTPUT_DIRECTORY\n";
        return 2;
    }
    const fs::path out = argv[1];
    fs::create_directories(out);
    const auto t0 = std::chrono::steady_clock::now();
    std::string failure;
    try {
        known_answers();
        stream_layout();
        conversions();
        bernoulli_exactness();
        stochastic_rounding_statistics();
        configuration_refusals();
        debit_sweeps();
        decision_refusals();
        partition_independence();
        owner_ids_ledger_checkpoint(out);
        staging_and_commit_refusals();
        oracle_dump(out);
    } catch (const std::exception& e) {
        failure = e.what();
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::ofstream r(out / "core-results.json");
    r << "{\"program\":\"test_core\",\"checks\":" << checks << ",\"passed\":" << (failure.empty() ? "true" : "false")
      << ",\"failure\":" << jstr(failure) << ",\"seconds\":" << seconds;
    for (const auto& n : notes) r << "," << n;
    r << "}\n";
    std::cout << "test_core checks=" << checks << " " << (failure.empty() ? "PASS" : "FAIL: " + failure) << "\n";
    return failure.empty() ? 0 : 1;
}
