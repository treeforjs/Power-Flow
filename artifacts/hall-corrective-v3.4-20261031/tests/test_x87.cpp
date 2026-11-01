// Hall v3.4 CPU qualification of the x87 extended emulation (gpu/x87_v34.hpp) and of the deposit support /
// census comparisons built on it (gpu/support_v34.hpp), against the host's own long double arithmetic and
// the sealed DepositSupportBuilder (census_charge.hpp).
//   test_x87 ops OUT [N]          random add/sub/mul/floor/compare/conversion vs long double (N per class)
//   test_x87 support OUT [ROUNDS] random carrier sets on the 200 and 50 um meshes vs DepositSupportBuilder
//   test_x87 compare OUT [N]      same_net_charge / free_charge_guard node predicates vs the host forms
#include "root_api.hpp"
#include "gpu/support_host_v34.hpp"
namespace hc = hall_corrective_v3;
#include "motion_fixtures.hpp"
#include <fstream>
#include <iostream>
#include <random>

namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
using J = hc::Json;

namespace {
struct Results {
    J cases = J::array();
    int failed = 0;
    void record(const std::string& group, const std::string& name, bool ok, J evidence = J::object()) {
        cases.push_back({{"group", group}, {"name", name}, {"passed", ok}, {"evidence", std::move(evidence)}});
        if (!ok) ++failed;
        std::cout << group << " | " << name << ": " << (ok ? "pass" : "FAIL") << std::endl;
        if (!ok) std::cout << "  evidence: " << cases.back()["evidence"].dump().substr(0, 3000) << std::endl;
    }
};
bool same_x(x87::X a, x87::X b) { return a.m == b.m && (a.m == 0 || (a.e == b.e && a.s == b.s)); }
std::string hexld(long double v) { char buf[64]; std::snprintf(buf, sizeof buf, "%La", v); return buf; }
std::uint64_t bits(double d) { return std::bit_cast<std::uint64_t>(d); }

// Random binary64 with a random exponent in [-e_span, e_span] around 2^e0, random sign.
double rand_double(std::mt19937_64& rng, int e0, int e_span, bool sign) {
    std::uniform_int_distribution<int> ex(-e_span, e_span);
    std::uniform_int_distribution<std::uint64_t> fr(0, (1ull << 52) - 1);
    const int e = e0 + ex(rng);
    double v = std::ldexp(1.0 + double(fr(rng)) * 0x1p-52, e);
    if (rng() & 7) {} else v = std::ldexp(double(fr(rng) | 1), e - 52);   // odd significands, fewer bits
    if (sign && (rng() & 1)) v = -v;
    return v;
}
long double rand_ld(std::mt19937_64& rng, int e0, int e_span, bool sign) {
    const double a = rand_double(rng, e0 / 2, e_span / 2, sign), b = rand_double(rng, e0 - e0 / 2, e_span - e_span / 2, false);
    switch (rng() % 4) {
        case 0: return (long double)a;
        case 1: return (long double)a * b;
        case 2: return (long double)a + (long double)b * 0x1p-30L;
        default: return ((long double)a * b) * (1 + 1e-9L);
    }
}

int ops(const std::filesystem::path& out, long N) {
    Results R;
    std::mt19937_64 rng(20261031);
    J ev = J::object();
    long bad_add = 0, bad_sub = 0, bad_mul = 0, bad_floor = 0, bad_le = 0, bad_rne = 0, bad_ru = 0, bad_conv = 0;
    std::string first;
    for (long t = 0; t < N; ++t) {
        const int span = int(rng() % 8) == 0 ? 200 : int(rng() % 140);
        const long double a = rand_ld(rng, int(rng() % 200) - 100, 4, true), b = rand_ld(rng, int(rng() % 200) - 100, span, true);
        long double c = a, d = b;
        if (rng() % 3 == 0) d = a * (1 + 0x1p-60L * (long double)(int(rng() % 9) - 4));   // near-cancellation
        if (rng() % 11 == 0) d = -a;
        const x87::X xa = sup34::from_long_double(c), xb = sup34::from_long_double(d);
        if (!same_x(sup34::from_long_double(sup34::to_long_double(xa)), xa)) ++bad_conv;
        const long double s = c + d, df = c - d, p = c * d;
        if (!same_x(x87::add(xa, xb), sup34::from_long_double(s))) { if (!bad_add++) first = "add " + hexld(c) + " " + hexld(d); }
        if (!same_x(x87::sub(xa, xb), sup34::from_long_double(df))) { if (!bad_sub++ && first.empty()) first = "sub " + hexld(c) + " " + hexld(d); }
        if (!same_x(x87::mul(xa, xb), sup34::from_long_double(p))) { if (!bad_mul++ && first.empty()) first = "mul " + hexld(c) + " " + hexld(d); }
        if (x87::le(xa, xb) != (c <= d)) ++bad_le;
        // floor of moderately sized values
        const long double f = c * 0x1p-0L;
        if (std::fabs(f) < 1e15L) { if (x87::floor_int(xa) != (long long)std::floor(f)) { if (!bad_floor++ && first.empty()) first = "floor " + hexld(f); } }
        // conversions to binary64 (positive values for RU)
        const long double pos = std::fabs(p);
        if (pos > 0x1p-1000L && pos < 0x1p1000L) {
            const x87::X xp = sup34::from_long_double(pos);
            if (bits(x87::to_double(xp, false)) != bits(double(pos))) ++bad_rne;
            double r = double(pos); if ((long double)r < pos) r = std::nextafter(r, INFINITY);
            if (bits(x87::to_double(xp, true)) != bits(r)) ++bad_ru;
        }
    }
    // doubles: exact conversion and the sealed difference/sum chains of the support (cell coordinates)
    long bad_cell = 0;
    for (long t = 0; t < N; ++t) {
        std::uniform_real_distribution<double> u(0, 400);
        const double a = u(rng), b = a + (rng() % 4 == 0 ? 0.0 : std::ldexp(double(int(rng() % 2001) - 1000), -int(rng() % 60)));
        const long double delta = std::numeric_limits<double>::epsilon() * (3.0L * double(rng() % 400) + 8);
        long double R = 0; R += std::abs((long double)b - a) + 2 * delta;
        const long double lo = (long double)std::min(a, b) - delta, hi = (long double)std::max(a, b) + delta;
        x87::X xR = x87::add(x87::zero(), x87::add(x87::abs(x87::sub(x87::from_double(b), x87::from_double(a))), sup34::from_long_double(2 * delta)));
        const x87::X xlo = x87::sub(x87::from_double(std::min(a, b)), sup34::from_long_double(delta));
        const x87::X xhi = x87::add(x87::from_double(std::max(a, b)), sup34::from_long_double(delta));
        if (!same_x(xR, sup34::from_long_double(R)) || !same_x(xlo, sup34::from_long_double(lo)) || !same_x(xhi, sup34::from_long_double(hi)) ||
            x87::floor_int(xlo) != (long long)std::floor(lo) || x87::floor_int(xhi) != (long long)std::floor(hi)) ++bad_cell;
    }
    ev = {{"trials", N}, {"add", bad_add}, {"sub", bad_sub}, {"mul", bad_mul}, {"floor", bad_floor}, {"le", bad_le}, {"to_double_rne", bad_rne},
          {"to_double_ru", bad_ru}, {"conversion", bad_conv}, {"cell_chain", bad_cell}, {"first_difference", first}};
    R.record("x87", "add/sub/mul/floor/le/to_double bit-identical to the host long double (random operands, near cancellation, wide exponent spread)",
             bad_add + bad_sub + bad_mul + bad_floor + bad_le + bad_rne + bad_ru + bad_conv + bad_cell == 0, ev);
    std::filesystem::create_directories(out);
    std::ofstream(out / "x87-ops.json") << J{{"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// Random carriers near the gap, small displacements (some exactly on node planes), endpoints.
struct Set { std::vector<n::ArenaRateEntry> census; std::vector<n::GapParticle> endpoints; };
Set random_set(const axial::Mesh& g, std::size_t count, std::uint64_t seed, int flavour) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> ux(.0035, .0055), uz(-.004, .004), ue(0, 1);
    const double offsets[] = {0., 1e-15, -1e-15, 1e-12, -1e-12, 4e-16, -7e-16};
    Set s;
    const double h[3] = {g.dx, g.dy, g.dz};
    for (std::size_t k = 0; k < count; ++k) {
        n::PicVector x{ux(rng), ux(rng), uz(rng)};
        if (flavour >= 1 && (rng() % 3 == 0)) {   // on a node plane (native frame: plane = side/2 + (i - n/2) h ... via field frame)
            const int axis = int(rng() % 3);
            const double origin = axis < 2 ? -g.side / 2 : g.zmin;
            const double shift = axis == 0 ? g.side / 2 : axis == 1 ? g.side / 2 : 0;
            const double cell = std::floor(((x[axis] - shift) - origin) / h[axis]);
            x[axis] = origin + cell * h[axis] + shift + offsets[rng() % 7];
        }
        double w = flavour == 2 ? std::ldexp(1.0 + ue(rng), int(rng() % 40)) : (rng() % 5 == 0 ? 1e9 : 1.0);
        double q = rng() % 7 == 0 ? 0.0 : (rng() & 1 ? -n::pic_e : n::pic_e);
        if (flavour == 2 && rng() % 9 == 0) q *= 2.5;
        n::PicVector e = x;
        const double frac = flavour == 3 ? .3 : .19;
        for (int d = 0; d < 3; ++d) e[d] += (ue(rng) - .5) * 2 * frac * h[d];
        if (rng() % 13 == 0) e = x;   // no displacement
        s.census.push_back(mf::row(k + 1, "p" + std::to_string(k), q == 0 ? "H2O" : (q < 0 ? "electron" : "H2O+"), x, {0, 0, 0}, q, 1e-30, w));
        n::GapParticle p = s.census.back().native; p.position_m = e; s.endpoints.push_back(p);
    }
    return s;
}

int support(const std::filesystem::path& out, long rounds) {
    Results R;
    J ev = J::array();
    long total = 0, nodes_compared = 0, mismatched = 0, refusal_mismatch = 0;
    for (int transverse : {45, 180}) {
        const axial::Mesh g(transverse, 4, .009, 1);
        const auto o = mf::geometry(g);
        const auto geo = sup34::geo_of(g, o);
        for (long r = 0; r < rounds; ++r) {
            const int flavour = int(r % 4);
            const std::size_t count = r % 5 == 0 ? 20000 : 3000;
            auto set = random_set(g, count, 1000 * transverse + r, flavour);
            std::string ref_error; std::vector<double> ref;
            try { ref = hc::DepositSupportBuilder::build(g, o, set.census, set.endpoints).absolute_C(); } catch (const std::exception& e) { ref_error = e.what(); }
            std::vector<double> w, q, xs, xe;
            for (std::size_t i = 0; i < set.census.size(); ++i) {
                w.push_back(set.census[i].native.weight); q.push_back(set.census[i].native.charge_C);
                for (int d = 0; d < 3; ++d) { xs.push_back(set.census[i].native.position_m[d]); xe.push_back(set.endpoints[i].position_m[d]); }
            }
            const auto em = sup34::emulate(geo, g.nn, g.nx, g.ny, w, q, xs, xe);
            ++total;
            const char* texts[] = {"", "deposit support charge operand", "deposit support nonfinite coordinate",
                                   "deposit support path exceeds the motion displacement guard"};
            if (!ref_error.empty() || em.code != sup34::ok) {
                if (ref_error != texts[em.code]) { ++refusal_mismatch; ev.push_back({{"mesh", transverse}, {"round", r}, {"ref", ref_error}, {"emulated", em.code}}); }
                continue;
            }
            long bad = 0; long first = -1;
            for (int k = 0; k < g.nn; ++k) if (bits(ref[std::size_t(k)]) != bits(em.absolute_C[std::size_t(k)])) { if (first < 0) first = k; ++bad; }
            nodes_compared += g.nn; mismatched += bad;
            if (bad) ev.push_back({{"mesh", transverse}, {"round", r}, {"nodes_differing", bad}, {"first", first},
                                   {"ref", ref[std::size_t(first)]}, {"emulated", em.absolute_C[std::size_t(first)]}});
        }
        // refusal: a path beyond one cell sum
        auto set = random_set(g, 50, 77, 0);
        set.endpoints[17].position_m[0] += 3 * g.dx;
        set.census[17].native.charge_C = -n::pic_e;
        std::string ref_error;
        try { (void)hc::DepositSupportBuilder::build(g, o, set.census, set.endpoints); } catch (const std::exception& e) { ref_error = e.what(); }
        std::vector<double> w, q, xs, xe;
        for (std::size_t i = 0; i < set.census.size(); ++i) {
            w.push_back(set.census[i].native.weight); q.push_back(set.census[i].native.charge_C);
            for (int d = 0; d < 3; ++d) { xs.push_back(set.census[i].native.position_m[d]); xe.push_back(set.endpoints[i].position_m[d]); }
        }
        const auto em = sup34::emulate(geo, g.nn, g.nx, g.ny, w, q, xs, xe);
        R.record("support", "mesh " + std::to_string(transverse) + ": path beyond the displacement guard refused with the sealed text at the same carrier",
                 ref_error == "deposit support path exceeds the motion displacement guard" && em.code == sup34::path_exceeds && em.first_failure == 17,
                 {{"ref", ref_error}, {"code", em.code}, {"first", em.first_failure}});
    }
    R.record("support", "emulated deposit support bit-identical to DepositSupportBuilder (x87) on every node of every random set",
             mismatched == 0 && refusal_mismatch == 0,
             {{"sets", total}, {"nodes_compared", nodes_compared}, {"nodes_differing", mismatched}, {"refusal_mismatch", refusal_mismatch}, {"details", ev}});
    std::filesystem::create_directories(out);
    std::ofstream(out / "x87-support.json") << J{{"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

int compare(const std::filesystem::path& out, long N) {
    Results R;
    std::mt19937_64 rng(4242);
    long bad_same = 0, bad_guard = 0, same_true = 0, guard_true = 0;
    for (long t = 0; t < N; ++t) {
        const double base = rand_double(rng, -60, 10, true);
        hc::Words a(base, base * 0x1p-60 * double(int(rng() % 7) - 3));
        hc::Words b = a;
        const int kind = int(rng() % 4);
        if (kind == 1) b = a + hc::Words(base * std::ldexp(1.0, -int(rng() % 60)) * (rng() & 1 ? 1 : -1));
        if (kind == 2) b = hc::Words(rand_double(rng, -60, 10, true));
        const double ga = std::fabs(a.rounded()) * (1 + double(rng() % 3)), gb = std::fabs(b.rounded()) * (1 + double(rng() % 3)) + (kind == 3 ? std::fabs(base) * 1e-12 : 0);
        hc::FreeCharge A{{a}, {ga}}, B{{b}, {gb}};
        const bool host = hc::same_net_charge(A, B);
        const bool emu = sup34::same_net_node(a.hi, a.lo, ga, b.hi, b.lo, gb);
        if (host != emu) ++bad_same;
        same_true += host;
        // free charge guard: gross slightly above/below |net|
        const double g2 = std::fabs(a.rounded()) * (1 + std::ldexp(double(int(rng() % 65) - 32), -52));
        const long double wide = a.wide();
        const bool hg = hc::canonical(a) && std::isfinite(g2) && g2 >= 0 && std::abs(wide) <= g2 * (1 + 32 * std::numeric_limits<double>::epsilon());
        if (hg != sup34::free_guard_node(a.hi, a.lo, g2)) ++bad_guard;
        guard_true += hg;
    }
    R.record("compare", "same_net_charge node predicate (x87) identical to the host form",
             bad_same == 0, {{"trials", N}, {"differing", bad_same}, {"host_true", same_true}});
    R.record("compare", "free_charge_guard node predicate (x87) identical to the host form",
             bad_guard == 0, {{"trials", N}, {"differing", bad_guard}, {"host_true", guard_true}});
    std::filesystem::create_directories(out);
    std::ofstream(out / "x87-compare.json") << J{{"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    return R.failed ? 1 : 0;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "ops" && argc >= 3) return ops(argv[2], argc > 3 ? std::stol(argv[3]) : 2000000);
        if (mode == "support" && argc >= 3) return support(argv[2], argc > 3 ? std::stol(argv[3]) : 40);
        if (mode == "compare" && argc >= 3) return compare(argv[2], argc > 3 ? std::stol(argv[3]) : 2000000);
        std::cerr << "usage: test_x87 ops|support|compare OUT [N]\n";
        return 2;
    } catch (const std::exception& e) { std::cerr << "fatal: " << e.what() << '\n'; return 2; }
}
