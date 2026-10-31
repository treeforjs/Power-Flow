// test_cpu.cpp -- CPU reference tests of the merge/split transaction.
//  unit        Philox known answers, exact TwoProduct and superaccumulator against
//              OpenSSL BIGNUM, round-to-nearest of exact values, expansions
//  kinematics  random groups (non-relativistic and keV), identical velocities,
//              symmetric groups (P = 0, isotropic axis), relative residuals
//  split       split-only transactions: residuals exactly zero
//  ledger      1e6 carriers per species on the 45/90/180/360 meshes: exact count
//              (superaccumulator and native ExactParticleCount), exact charge, booked
//              momentum/energy residuals closing exactly, group weights exact, ID
//              allocation, outputs inside the group cell, relocation continuity and
//              relocation work, CPU time
//  partitions  1, 2, 5, 8 partitions, contiguous and shuffled: identical bytes
// Usage: test_cpu [--n N] [--meshes 45,90,180,360] [--partition-meshes 45,360]
//                 [--parts 1,2,5,8] [--out FILE.json]
#include "fixtures.hpp"
#include "merge_split.hpp"
#include "relocation.hpp"
#include "raikou/exact_particle_count.h"

#include <openssl/bn.h>

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace merge_split;
using Clock = std::chrono::steady_clock;

namespace {
int failures = 0;
std::ostringstream J;   // JSON body
bool first_field = true;
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
void field(const std::string& key, const std::string& raw) {
    J << (first_field ? "" : ",\n") << " \"" << key << "\": " << raw;
    first_field = false;
}
std::string q(const std::string& s) { return "\"" + s + "\""; }
double seconds(Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); }

// ---------------------------------------------------------------- BIGNUM helpers (independent exact arithmetic)
using BN = std::unique_ptr<BIGNUM, decltype(&BN_free)>;
BN bn() { return BN(BN_new(), BN_free); }
BN bn_of_double(double x) {   // signed integer of units 2^-1074
    std::uint64_t b;
    std::memcpy(&b, &x, 8);
    const unsigned ex = (b >> 52) & 2047;
    std::uint64_t m = b & ((1ull << 52) - 1);
    if (ex) m |= 1ull << 52;
    auto n = bn();
    BN_set_word(n.get(), m);
    if (ex) BN_lshift(n.get(), n.get(), int(ex) - 1);
    BN_set_negative(n.get(), (b >> 63) && m);
    return n;
}
BN bn_of_exact(const Exact& e) {
    std::string h = e.hex();
    bool neg = h[0] == '-';
    std::string digits = h.substr(neg ? 3 : (h == "0" ? 0 : 2));
    if (h == "0") digits = "0";
    BIGNUM* p = nullptr;
    BN_hex2bn(&p, digits.c_str());
    BN out(p, BN_free);
    BN_set_negative(out.get(), neg && !BN_is_zero(out.get()));
    return out;
}
bool bn_eq(const BN& a, const BN& b) { return BN_cmp(a.get(), b.get()) == 0; }

// ---------------------------------------------------------------- unit tests
void test_unit() {
    // Philox4x32-10 known answers (Random123 kat_vectors, as in the arrival owner tests)
    struct Kat { pfa_u32 c[4], k[2], r[4]; };
    const Kat kats[] = {
        {{0, 0, 0, 0}, {0, 0}, {0x6627e8d5, 0xe169c58d, 0xbc57ac4c, 0x9b00dbd8}},
        {{0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff}, {0xffffffff, 0xffffffff}, {0x408f276d, 0x41c83b0e, 0xa20bc7c6, 0x6d5451fd}},
        {{0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344}, {0xa4093822, 0x299f31d0}, {0xd16cfe09, 0x94fdcceb, 0x5001e420, 0x24126ea1}}};
    for (const auto& k : kats) {
        pfa::Block c{{k.c[0], k.c[1], k.c[2], k.c[3]}};
        pfa::Key2 key{{k.k[0], k.k[1]}};
        const auto r = pfa::philox4x32_10(c, key);
        check(r.v[0] == k.r[0] && r.v[1] == k.r[1] && r.v[2] == k.r[2] && r.v[3] == k.r[3], "philox known answer");
    }
    // exact products and sums against BIGNUM
    fixtures::Rng r(0xC0FFEE);
    auto rnd = [&](int lo, int hi) {
        const double m = 1 + double(r.g() >> 12) * 0x1p-52;
        const int ex = lo + int(r.g() % std::uint64_t(hi - lo + 1));
        return (r.g() & 1 ? -1 : 1) * std::ldexp(m, ex);
    };
    std::size_t products = 0, sums = 0, rounds = 0;
    for (int t = 0; t < 20000; ++t) {
        const double a = rnd(-120, 60), b = rnd(-120, 60), m = rnd(-110, -60);
        if (!ms_product3_safe(a, b, m)) continue;
        MsAcc acc;
        ms_acc_zero(&acc);
        ms_acc_add_product3(&acc, a, b, m);
        // BIGNUM: a*b*m exactly as integers: (Ma 2^ea)(Mb 2^eb)(Mm 2^em) -> units 2^-1074 needs scaling
        // compare instead (a*b*m) * 2^(2*1074) = A*B*M with A, B, M integer units.
        auto A = bn_of_double(a), B = bn_of_double(b), M = bn_of_double(m), X = bn_of_exact(Exact::of(acc));
        BN_CTX* ctx = BN_CTX_new();
        auto prod = bn();
        BN_mul(prod.get(), A.get(), B.get(), ctx);
        BN_mul(prod.get(), prod.get(), M.get(), ctx);
        BN_lshift(X.get(), X.get(), 2 * 1074);
        BN_CTX_free(ctx);
        check(bn_eq(prod, X), "exact three-factor product");
        ++products;
    }
    for (int t = 0; t < 200; ++t) {
        MsAcc acc;
        ms_acc_zero(&acc);
        auto ref = bn();
        const int n = 1 + int(r.g() % 5000);
        for (int i = 0; i < n; ++i) {
            const double x = rnd(-1074 + 60, 900);
            ms_acc_add(&acc, x);
            auto bx = bn_of_double(x);
            BN_add(ref.get(), ref.get(), bx.get());
        }
        const Exact e = Exact::of(acc);
        check(bn_eq(bn_of_exact(e), ref), "superaccumulator sum");
        ++sums;
        // round to nearest: |v - RN(v)| <= ulp/2 with ties to even, and RN within one ulp of hi
        const double x0 = e.rounded();
        auto diff = bn();
        auto bx0 = bn_of_double(x0);
        BN_sub(diff.get(), ref.get(), bx0.get());
        // half ulp of x0 in units 2^-1074
        std::uint64_t bits;
        std::memcpy(&bits, &x0, 8);
        const int ex = int((bits >> 52) & 2047);
        auto half = bn();
        BN_one(half.get());
        BN_lshift(half.get(), half.get(), std::max(ex, 1) - 1);   // ulp = 2^(max(ex,1)-1) units; half ulp shift-1
        BN_rshift1(half.get(), half.get());
        BN_set_negative(diff.get(), 0);
        check(BN_cmp(diff.get(), half.get()) <= 0, "round to nearest within half an ulp");
        // expansion recomposes exactly
        auto sum = bn();
        for (double c : e.expansion()) {
            auto bc = bn_of_double(c);
            BN_add(sum.get(), sum.get(), bc.get());
        }
        check(bn_eq(sum, ref), "canonical expansion is exact");
        ++rounds;
    }
    field("unit", "{\"philox_known_answers\": 3, \"exact_products\": " + std::to_string(products) + ", \"exact_sums\": " +
                      std::to_string(sums) + ", \"rounding_checks\": " + std::to_string(rounds) + "}");
    std::cout << "unit: products " << products << " sums " << sums << " rounding " << rounds << std::endl;
}

// ---------------------------------------------------------------- kinematics on small groups
MsParams one_bin_policy(std::uint64_t seed) {
    MsParams P = fixtures::policy(45, seed, 7, 1u << 30);
    P.nsub = 0;   // energy bins only, and one energy bin: every group is its whole cell
    for (unsigned s = 0; s < 2; ++s) {
        auto& sp = P.species[s];
        sp.n_thresh = 0;   // one energy bin
        sp.merge_above = 2;
        sp.merge_target = 2;
        sp.k_min = 3;
        sp.k_max = 64;
        sp.split_below = 0;
    }
    return P;
}
struct ResStats {
    double max_e = 0, max_p = 0;
    std::size_t groups = 0, accepted = 0, skipped = 0, isotropic = 0, zero_p = 0;
};
// One cell per group: carriers at random positions inside cell (cx, cy, kz) of the 45 mesh.
void group_cell_positions(fixtures::Rng& r, const MsParams& P, unsigned cell, double x[3]) {
    const unsigned i = cell % P.nx, j = (cell / P.nx) % P.ny, k = cell / (P.nx * P.ny);
    x[0] = (i + .05 + .9 * r.uniform()) * P.h[0];
    x[1] = (j + .05 + .9 * r.uniform()) * P.h[1];
    x[2] = P.origin[2] + (k + .05 + .9 * r.uniform()) * P.h[2];
}
void residual_stats(const MsParams& P, const Result& R, ResStats& st, const Carriers& C) {
    for (unsigned s = 0; s < P.n_species; ++s) {
        const auto& L = R.ledger[s];
        if (L.w_in.zero()) continue;
        const double ein = L.e_in.rounded();
        if (ein > 0) st.max_e = std::max(st.max_e, std::abs(L.e_res_rn) / ein);
        long double act = 0;   // sum |w m u_d| of the inputs
        for (std::size_t i = 0; i < C.size(); ++i)
            if (C.species[i] == s && std::binary_search(R.removed.begin(), R.removed.end(), C.id[i]))
                for (int d = 0; d < 3; ++d) act += std::abs((long double)C.w[i] * P.species[s].mass_kg * C.u[d][i]);
        for (int d = 0; d < 3; ++d)
            if (act > 0) st.max_p = std::max(st.max_p, double(std::abs((long double)L.p_res_rn[d]) / act));
        bool zp = true;
        for (int d = 0; d < 3; ++d) zp = zp && L.p_res[d].zero();
        if (zp) ++st.zero_p;
    }
    st.groups += R.counters.groups;
    st.accepted += R.counters.accepted_merge;
    st.skipped += R.counters.groups - R.counters.accepted_merge;
    st.isotropic += R.counters.isotropic;
}
std::string stats_json(const ResStats& s) {
    return "{\"groups\": " + std::to_string(s.groups) + ", \"accepted\": " + std::to_string(s.accepted) + ", \"skipped\": " +
           std::to_string(s.skipped) + ", \"isotropic\": " + std::to_string(s.isotropic) + ", \"transactions_with_zero_momentum_residual\": " +
           std::to_string(s.zero_p) + ", \"max_rel_energy_residual\": " + num(s.max_e) + ", \"max_rel_momentum_residual\": " + num(s.max_p) + "}";
}
void test_kinematics() {
    std::string out = "{";
    for (int kind = 0; kind < 4; ++kind) {
        ResStats st;
        fixtures::Rng r(1000 + kind);
        for (int t = 0; t < 400; ++t) {
            MsParams P = one_bin_policy(0xABCDEF00ull + std::uint64_t(t));
            Carriers C;
            const unsigned species = unsigned(t & 1);
            const unsigned k = kind == 2 ? 4u : 3u + unsigned(r.g() % 14);
            const unsigned cell = unsigned(r.g() % (45 * 45)) + 45 * 45 * 1;   // a kinetic layer cell
            double same[3];
            fixtures::velocity(r, species, same);
            double v2[3];
            fixtures::velocity(r, species, v2);
            for (unsigned i = 0; i < k; ++i) {
                double x[3], u[3];
                group_cell_positions(r, P, cell, x);
                if (kind == 0) fixtures::velocity(r, species, u);
                else if (kind == 1) for (int d = 0; d < 3; ++d) u[d] = same[d];
                else if (kind == 2) for (int d = 0; d < 3; ++d) u[d] = (i < 2 ? same[d] : v2[d]) * (i % 2 ? -1 : 1);
                else {   // relativistic electrons 50 keV - 2 MeV
                    double mag;
                    fixtures::speed_to_u(fixtures::me, (5e4 + 2e6 * r.uniform()) * fixtures::e, mag);
                    fixtures::isotropic(r, mag, u, .5);
                }
                const double w = kind == 2 ? 1e6 : (i % 3 ? 1e6 : 5e5);
                C.push(100 + i, kind == 3 ? 0 : species, MS_FLAG_MERGE, x, u, w);
            }
            const auto R = transact(P, C);
            residual_stats(P, R, st, C);
            if (kind == 2) check(R.counters.isotropic == R.counters.accepted_merge, "symmetric group uses the isotropic axis");
        }
        const char* names[] = {"maxwellian_or_tail", "identical_velocities", "symmetric_P_zero", "relativistic_50keV_2MeV"};
        out += std::string(kind ? ", " : "") + q(names[kind]) + ": " + stats_json(st);
        std::cout << "kinematics " << names[kind] << ": " << stats_json(st) << std::endl;
        if (kind == 2) check(st.zero_p == st.accepted, "P = 0 groups: momentum residual exactly zero");
        if (kind != 1) check(st.max_e < 1e-14 && st.max_p < 1e-14, "relative residuals at rounding level");
    }
    field("kinematics", out + "}");
}

// ---------------------------------------------------------------- split-only
void test_split_only(std::size_t n) {
    MsParams P = fixtures::policy(360, 0x51ull, 3, 2 * n + 1);
    for (unsigned s = 0; s < 2; ++s) {
        auto& sp = P.species[s];
        sp.merge_above = 0xFFFFFFF0u;
        sp.split_below = 1000;
        sp.split_target = 2000;
        sp.split_levels = 2;
    }
    const auto C = fixtures::make_set(0x5150, n);
    const auto R = transact(P, C);
    bool zero = true;
    for (const auto& L : R.ledger) {
        zero = zero && L.e_res.zero() && L.count_equal;
        for (int d = 0; d < 3; ++d) zero = zero && L.p_res[d].zero();
    }
    check(R.counters.accepted_split > 0 && zero, "split-only transaction: exact zero residuals");
    field("split_only", "{\"carriers\": " + std::to_string(C.size()) + ", \"splits\": " + std::to_string(R.counters.accepted_split) +
                            ", \"outputs\": " + std::to_string(R.counters.outputs) + ", \"all_residuals_exactly_zero\": " + (zero ? "true" : "false") + "}");
    std::cout << "split-only: splits " << R.counters.accepted_split << " zero residuals " << zero << std::endl;
}

// ---------------------------------------------------------------- ledger closure on large sets
std::string exact_json(const Exact& e) { return "{\"rn\": " + num(e.rounded()) + ", \"units_2^-1074\": " + q(e.hex()) + "}"; }

void test_ledger(int mesh, std::size_t n, std::string& out) {
    const std::uint64_t seed = 0x9e3779b97f4a7c15ull ^ std::uint64_t(mesh);
    const auto C = fixtures::make_set(0x1234 + std::uint64_t(mesh), n);
    MsParams P = fixtures::policy(mesh, seed, 4242, 2 * n + 1);
    const auto t0 = Clock::now();
    const auto R = transact(P, C);
    const auto t1 = Clock::now();
    const auto post = apply(C, R);
    std::ostringstream o;
    o << "{\"mesh\": " << mesh << ", \"carriers\": " << C.size() << ", \"cpu_transaction_s\": " << num(seconds(t0, t1))
      << ", \"cpu_s_per_1e6_carriers\": " << num(seconds(t0, t1) / (C.size() / 1e6));
    const auto& c = R.counters;
    o << ", \"segments\": " << c.segments << ", \"groups\": " << c.groups << ", \"merge_groups\": " << c.merge_groups
      << ", \"split_groups\": " << c.split_groups << ", \"accepted_merge\": " << c.accepted_merge << ", \"accepted_split\": " << c.accepted_split
      << ", \"skip_no_real_solution\": " << c.status[MS_SKIP_NO_REAL_SOLUTION] << ", \"skip_operand_range\": " << c.status[MS_SKIP_OPERAND_RANGE]
      << ", \"skip_rng\": " << c.status[MS_SKIP_RNG_EXHAUSTED] << ", \"formation_refused\": " << c.formation_refused
      << ", \"centroid_fallback\": " << c.centroid_fallback << ", \"isotropic\": " << c.isotropic << ", \"inputs_consumed\": " << c.inputs_consumed
      << ", \"outputs\": " << c.outputs << ", \"relocations\": " << c.relocations << ", \"census_after\": " << post.size()
      << ", \"digest\": " << q(R.digest_sha256);
    // count, charge: exact per species over the whole census, superaccumulator and native ExactParticleCount
    const auto wb = exact_weights(P, C), wa = exact_weights(P, post);
    bool counts = true, native_counts = true, charges = true;
    for (unsigned s = 0; s < 2; ++s) {
        counts = counts && wb[s] == wa[s] && R.ledger[s].count_equal;
        raikou::mfem_backend::ExactParticleCount nb, na, gin, gout;
        MsAcc qb, qa;
        ms_acc_zero(&qb);
        ms_acc_zero(&qa);
        for (std::size_t i = 0; i < C.size(); ++i)
            if (C.species[i] == s) { nb = nb.credited(C.w[i]); ms_acc_add_product2(&qb, C.w[i], P.species[s].charge_C); }
        for (std::size_t i = 0; i < post.size(); ++i)
            if (post.species[i] == s) { na = na.credited(post.w[i]); ms_acc_add_product2(&qa, post.w[i], P.species[s].charge_C); }
        native_counts = native_counts && nb == na;
        charges = charges && Exact::of(qb) == Exact::of(qa);
    }
    check(counts, "exact count (superaccumulator) conserved");
    check(native_counts, "exact count (native ExactParticleCount) conserved");
    check(charges, "exact charge conserved");
    // every merge group's weight is the exact sum of its inputs; halves are exact; split clones exact
    std::size_t group_checks = 0;
    bool groups_exact = true;
    {
        // map removed inputs per group through the relocation records (merge) and the canonical first ID (split)
        std::vector<raikou::mfem_backend::ExactParticleCount> per(R.groups.size());
        std::vector<std::size_t> idx(C.size());
        std::vector<std::pair<std::uint64_t, std::size_t>> by_id;
        by_id.reserve(C.size());
        for (std::size_t i = 0; i < C.size(); ++i) by_id.push_back({C.id[i], i});
        std::sort(by_id.begin(), by_id.end());
        auto find = [&](std::uint64_t id) { return std::lower_bound(by_id.begin(), by_id.end(), std::make_pair(id, std::size_t(0)))->second; };
        for (std::size_t k = 0; k < R.reloc.size(); ++k) per[R.reloc.group[k]] = per[R.reloc.group[k]].credited(C.w[find(R.reloc.id[k])]);
        for (std::size_t g = 0; g < R.groups.size(); ++g) {
            const auto& gr = R.groups[g];
            if (gr.status != MS_OK) continue;
            raikou::mfem_backend::ExactParticleCount outw;
            for (std::uint32_t j = 0; j < gr.n_out; ++j) outw = outw.credited(R.outputs.w[gr.out_offset + j]);
            if (gr.kind == MS_KIND_MERGE) {
                groups_exact = groups_exact && per[g] == outw && per[g] == raikou::mfem_backend::ExactParticleCount(gr.W) &&
                               R.outputs.w[gr.out_offset] * 2 == gr.W;
            } else {
                groups_exact = groups_exact && outw == raikou::mfem_backend::ExactParticleCount(C.w[find(R.group_first_id[g])]);
            }
            ++group_checks;
        }
    }
    check(groups_exact, "per-group ExactParticleCount equality");
    // booked residuals close exactly: in + booked == out; booked expansion recomposes exactly
    bool closure = true, expansion_exact = true;
    double rel_e[2] = {0, 0}, rel_p[2] = {0, 0};
    std::ostringstream lj;
    for (unsigned s = 0; s < 2; ++s) {
        const auto& L = R.ledger[s];
        closure = closure && (L.e_in + L.e_res) == L.e_out;
        for (int d = 0; d < 3; ++d) closure = closure && (L.p_in[d] + L.p_res[d]) == L.p_out[d];
        auto recompose = [&](const std::vector<double>& xs, const Exact& v) {
            MsAcc a;
            ms_acc_zero(&a);
            for (double x : xs) ms_acc_add(&a, x);
            return Exact::of(a) == v;
        };
        expansion_exact = expansion_exact && recompose(L.e_res_expansion, L.e_res);
        for (int d = 0; d < 3; ++d) expansion_exact = expansion_exact && recompose(L.p_res_expansion[d], L.p_res[d]);
        // independent BIGNUM recomputation of the input/output energy totals
        rel_e[s] = L.e_in.zero() ? 0 : std::abs(L.e_res_rn) / L.e_in.rounded();
        long double act = 0;
        for (std::size_t k = 0; k < R.removed.size(); ++k) (void)k;
        for (std::size_t i = 0; i < C.size(); ++i)
            if (C.species[i] == s && std::binary_search(R.removed.begin(), R.removed.end(), C.id[i]))
                for (int d = 0; d < 3; ++d) act += std::abs((long double)C.w[i] * P.species[s].mass_kg * C.u[d][i]);
        double pm = 0;
        for (int d = 0; d < 3; ++d) pm = std::max(pm, std::abs(L.p_res_rn[d]));
        rel_p[s] = act > 0 ? double(pm / act) : 0;
        lj << (s ? ", " : "") << "{\"species\": " << s << ", \"w_in\": " << exact_json(L.w_in) << ", \"w_out\": " << exact_json(L.w_out)
           << ", \"e_in_J\": " << exact_json(L.e_in) << ", \"e_out_J\": " << exact_json(L.e_out)
           << ", \"merge_energy_residual_J\": " << exact_json(L.e_res) << ", \"merge_energy_residual_expansion_terms\": " << L.e_res_expansion.size()
           << ", \"merge_momentum_residual_kg_m_s\": [" << exact_json(L.p_res[0]) << ", " << exact_json(L.p_res[1]) << ", " << exact_json(L.p_res[2]) << "]"
           << ", \"p_in_kg_m_s\": [" << num(L.p_in[0].rounded()) << ", " << num(L.p_in[1].rounded()) << ", " << num(L.p_in[2].rounded()) << "]"
           << ", \"rel_energy_residual\": " << num(rel_e[s]) << ", \"max_momentum_residual_over_sum_abs\": " << num(rel_p[s]) << "}";
    }
    check(closure, "booked residuals close exactly (in + residual == out)");
    check(expansion_exact, "booked residual expansions are exact");
    // independent BIGNUM check of one species' energy totals
    bool bn_ok = true;
    {
        BN_CTX* ctx = BN_CTX_new();
        for (unsigned s = 0; s < 2; ++s) {
            auto ein = bn(), eout = bn();
            for (std::size_t i = 0; i < C.size(); ++i) {
                if (C.species[i] != s || !std::binary_search(R.removed.begin(), R.removed.end(), C.id[i])) continue;
                const double ke = ms_kinetic_energy(P.species[s].mass_kg, C.u[0][i], C.u[1][i], C.u[2][i]);
                auto a = bn_of_double(C.w[i]), b = bn_of_double(ke), pr = bn();
                BN_mul(pr.get(), a.get(), b.get(), ctx);
                BN_add(ein.get(), ein.get(), pr.get());
            }
            for (std::size_t i = 0; i < R.outputs.size(); ++i) {
                if (R.outputs.species[i] != s) continue;
                const double ke = ms_kinetic_energy(P.species[s].mass_kg, R.outputs.u[0][i], R.outputs.u[1][i], R.outputs.u[2][i]);
                auto a = bn_of_double(R.outputs.w[i]), b = bn_of_double(ke), pr = bn();
                BN_mul(pr.get(), a.get(), b.get(), ctx);
                BN_add(eout.get(), eout.get(), pr.get());
            }
            auto ei = bn_of_exact(R.ledger[s].e_in), eo = bn_of_exact(R.ledger[s].e_out);
            BN_lshift(ei.get(), ei.get(), 1074);
            BN_lshift(eo.get(), eo.get(), 1074);
            bn_ok = bn_ok && bn_eq(ein, ei) && bn_eq(eout, eo);
        }
        BN_CTX_free(ctx);
    }
    check(bn_ok, "energy totals equal the independent BIGNUM sums");
    // ID allocation and output placement
    bool ids = R.outputs.size() == R.counters.outputs && R.next_id_after == P.base_id + R.outputs.size();
    for (std::size_t i = 0; i < R.outputs.size() && ids; ++i) ids = R.outputs.id[i] == P.base_id + i;
    bool placed = true;
    for (std::size_t g = 0; g < R.groups.size(); ++g) {
        const auto& gr = R.groups[g];
        if (gr.status != MS_OK) continue;
        for (std::uint32_t j = 0; j < gr.n_out; ++j) {
            const std::size_t i = gr.out_offset + j;
            std::uint32_t cell = 0;
            int open = 0;
            placed = placed && ms_cell_of(&P, R.outputs.x[0][i], R.outputs.x[1][i], R.outputs.x[2][i], &cell, &open) && cell == gr.cell;
        }
    }
    check(ids, "output IDs are base + rank in canonical order");
    check(placed, "outputs lie in their group cell");
    // relocation current: continuity and work
    const HallMesh hm(mesh);
    const double dt = 1e-14;
    const auto t2 = Clock::now();
    const auto Jr = relocation_current(hm, R.reloc, dt);
    const auto t3 = Clock::now();
    const auto cont = relocation_continuity(hm, P, C, post, Jr, dt);
    const double E0[3] = {1e5, -2e5, 1e6}, Z[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    // symmetric (curl-free) gradient with off-diagonal terms: the edge (Whitney) interpolant of a
    // diagonal gradient is constant along the edge direction inside a cell, the off-diagonal part
    // is what produces second-order relocation work (sum Q (x - c)_a (x - c)_b G_ab type terms)
    const double E00[3] = {0, 0, 0}, G[3][3] = {{1e8, 5e7, 2e7}, {5e7, -2e8, 3e7}, {2e7, 3e7, 1e8}};
    const double w_uniform = relocation_work_J(test_edge_field(hm, E0, Z), Jr, dt);
    const double w_grad = relocation_work_J(test_edge_field(hm, E00, G), Jr, dt);
    long double first_order = 0;   // sum |Q E0.D|: the uniform-field work without the centroid cancellation
    for (std::size_t k = 0; k < R.reloc.size(); ++k)
        first_order += std::abs((long double)R.reloc.charge[k] * (E0[0] * R.reloc.D[0][k] + E0[1] * R.reloc.D[1][k] + E0[2] * R.reloc.D[2][k]));
    const double ein_total = R.ledger[0].e_in.rounded() + R.ledger[1].e_in.rounded();
    check(cont.max_ratio < 1e-2, "relocation continuity within 1% of the Hall relative tolerance term");
    o << ", \"counts_exact\": " << (counts ? "true" : "false") << ", \"native_exact_particle_count_equal\": " << (native_counts ? "true" : "false")
      << ", \"charge_exact\": " << (charges ? "true" : "false") << ", \"group_weight_checks\": " << group_checks
      << ", \"groups_exact\": " << (groups_exact ? "true" : "false") << ", \"residual_closure_exact\": " << (closure ? "true" : "false")
      << ", \"residual_expansions_exact\": " << (expansion_exact ? "true" : "false") << ", \"bignum_energy_totals_equal\": " << (bn_ok ? "true" : "false")
      << ", \"ids_base_plus_rank\": " << (ids ? "true" : "false") << ", \"outputs_in_group_cell\": " << (placed ? "true" : "false")
      << ", \"ledger\": [" << lj.str() << "]"
      << ", \"relocation\": {\"dt_s\": " << num(dt) << ", \"deposit_s\": " << num(seconds(t2, t3)) << ", \"max_continuity_defect_C\": " << num(cont.max_defect_C)
      << ", \"max_defect_over_1e-12_scale\": " << num(cont.max_ratio) << ", \"max_node_charge_moved_C\": " << num(cont.max_net_change_C)
      << ", \"nodes_changed\": " << cont.nodes_changed << ", \"work_uniform_field_J\": " << num(w_uniform)
      << ", \"uniform_field_first_order_sum_abs_J\": " << num(double(first_order)) << ", \"work_gradient_field_J\": " << num(w_grad)
      << ", \"input_kinetic_energy_J\": " << num(ein_total) << "}}";
    out += o.str();
    std::cout << "ledger mesh " << mesh << ": " << o.str().substr(0, 900) << " ..." << std::endl;
}

// ---------------------------------------------------------------- spatial sub-bins (relocation length)
void test_subbins(std::size_t n, std::string& out) {
    const int mesh = 45;
    const auto C = fixtures::make_set(0x1234 + std::uint64_t(mesh), n);
    const HallMesh hm(mesh);
    const double dt = 1e-14;
    const double E00[3] = {0, 0, 0}, G[3][3] = {{1e8, 5e7, 2e7}, {5e7, -2e8, 3e7}, {2e7, 3e7, 1e8}};
    const auto E = test_edge_field(hm, E00, G);
    const unsigned subs[][3] = {{1, 1, 1}, {1, 1, 4}, {1, 1, 8}, {2, 2, 8}, {4, 4, 32}};
    std::ostringstream o;
    o << "[";
    bool first = true;
    for (const auto& sb : subs) {
        MsParams P = fixtures::policy(mesh, 0x9e3779b97f4a7c15ull ^ std::uint64_t(mesh), 4242, 2 * n + 1);
        for (int d = 0; d < 3; ++d) P.sub[d] = sb[d];
        const auto t0 = Clock::now();
        const auto R = transact(P, C);
        const double s = seconds(t0, Clock::now());
        const auto Jr = relocation_current(hm, R.reloc, dt);
        const double w = relocation_work_J(E, Jr, dt);
        double maxD[3] = {0, 0, 0};
        long double rms = 0;
        for (std::size_t k = 0; k < R.reloc.size(); ++k) {
            long double r2 = 0;
            for (int d = 0; d < 3; ++d) {
                maxD[d] = std::max(maxD[d], std::abs(R.reloc.D[d][k]));
                r2 += (long double)R.reloc.D[d][k] * R.reloc.D[d][k];
            }
            rms += r2;
        }
        if (R.reloc.size()) rms = std::sqrt(rms / R.reloc.size());
        const double ein = R.ledger[0].e_in.rounded() + R.ledger[1].e_in.rounded();
        bool closure = true;
        for (const auto& L : R.ledger) {
            closure = closure && L.count_equal && (L.e_in + L.e_res) == L.e_out;
            for (int d = 0; d < 3; ++d) closure = closure && (L.p_in[d] + L.p_res[d]) == L.p_out[d];
        }
        check(closure, "sub-binned transaction ledger closure");
        o << (first ? "" : ", ") << "{\"sub\": [" << sb[0] << ", " << sb[1] << ", " << sb[2] << "], \"groups\": " << R.counters.groups
          << ", \"accepted_merge\": " << R.counters.accepted_merge << ", \"accepted_split\": " << R.counters.accepted_split
          << ", \"inputs_consumed\": " << R.counters.inputs_consumed << ", \"census_after\": " << (C.size() - R.removed.size() + R.outputs.size())
          << ", \"max_abs_D_m\": [" << num(maxD[0]) << ", " << num(maxD[1]) << ", " << num(maxD[2]) << "], \"rms_D_m\": " << num(double(rms))
          << ", \"relocation_work_gradient_field_J\": " << num(w) << ", \"merge_input_kinetic_J\": " << num(ein)
          << ", \"work_over_input_kinetic\": " << num(ein > 0 ? w / ein : 0) << ", \"cpu_s\": " << num(s) << ", \"digest\": " << q(R.digest_sha256) << "}";
        first = false;
        std::cout << "sub " << sb[0] << sb[1] << sb[2] << ": merges " << R.counters.accepted_merge << " census " << (C.size() - R.removed.size() + R.outputs.size())
                  << " rmsD " << double(rms) << " work " << w << " rel " << (ein > 0 ? w / ein : 0) << std::endl;
    }
    out += o.str() + "]";
}

// ---------------------------------------------------------------- partitions
void test_partitions(int mesh, std::size_t n, const std::vector<unsigned>& parts, std::string& out) {
    const auto C = fixtures::make_set(0x777 + std::uint64_t(mesh), n);
    const MsParams P = fixtures::policy(mesh, 0xD1CEull, 99, 2 * n + 1);
    const auto ref = transact(P, C);
    const auto ref_bytes = encode(ref);
    std::ostringstream o;
    o << "{\"mesh\": " << mesh << ", \"carriers\": " << C.size() << ", \"reference_digest\": " << q(ref.digest_sha256)
      << ", \"groups\": " << ref.counters.groups << ", \"outputs\": " << ref.counters.outputs << ", \"runs\": [";
    bool all = true;
    bool first = true;
    for (bool shuffled : {false, true})
        for (unsigned p : parts) {
            PartitionPlan plan{p, shuffled, 0xFEED0000ull + p};
            const auto t0 = Clock::now();
            const auto R = transact(P, C, plan);
            const auto t1 = Clock::now();
            const bool same = R.digest_sha256 == ref.digest_sha256 && encode(R) == ref_bytes;
            all = all && same;
            o << (first ? "" : ", ") << "{\"parts\": " << p << ", \"assignment\": " << q(shuffled ? "shuffled" : "contiguous")
              << ", \"digest\": " << q(R.digest_sha256) << ", \"bytes_identical\": " << (same ? "true" : "false") << ", \"s\": " << num(seconds(t0, t1)) << "}";
            first = false;
        }
    check(all, "partition independence on mesh " + std::to_string(mesh));
    o << "], \"all_identical\": " << (all ? "true" : "false") << "}";
    out += o.str();
    std::cout << "partitions mesh " << mesh << ": all identical " << all << std::endl;
}

std::vector<unsigned> list(const char* s) {
    std::vector<unsigned> v;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, ',')) v.push_back(unsigned(std::stoul(t)));
    return v;
}
}  // namespace

int main(int argc, char** argv) {
    std::size_t n = 1000000;
    std::vector<unsigned> meshes{45, 90, 180, 360}, pmeshes{45, 360}, parts{1, 2, 5, 8};
    std::string out_path;
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string a = argv[i];
        if (a == "--n") n = std::stoull(argv[i + 1]);
        else if (a == "--meshes") meshes = list(argv[i + 1]);
        else if (a == "--partition-meshes") pmeshes = list(argv[i + 1]);
        else if (a == "--parts") parts = list(argv[i + 1]);
        else if (a == "--out") out_path = argv[i + 1];
        else { std::cerr << "unknown option " << a << std::endl; return 2; }
    }
    const auto start = Clock::now();
    try {
        field("schema", q("merge-split-cpu-tests-v1"));
        field("carriers_per_species", std::to_string(n));
        test_unit();
        test_kinematics();
        test_split_only(std::max<std::size_t>(n / 10, 1000));
        std::string ledger = "[";
        for (std::size_t k = 0; k < meshes.size(); ++k) {
            if (k) ledger += ", ";
            test_ledger(int(meshes[k]), n, ledger);
        }
        field("ledger", ledger + "]");
        std::string sub;
        test_subbins(n, sub);
        field("spatial_subbins_mesh45", sub);
        std::string part = "[";
        for (std::size_t k = 0; k < pmeshes.size(); ++k) {
            if (k) part += ", ";
            test_partitions(int(pmeshes[k]), n, parts, part);
        }
        field("partitions", part + "]");
    } catch (const std::exception& ex) {
        ++failures;
        std::cout << "FAIL exception: " << ex.what() << std::endl;
        field("exception", q(ex.what()));
    }
    field("failures", std::to_string(failures));
    field("wall_s", num(seconds(start, Clock::now())));
    const std::string doc = "{\n" + J.str() + "\n}\n";
    if (!out_path.empty()) std::ofstream(out_path) << doc;
    std::cout << (failures ? "FAILED " : "PASSED ") << failures << " failures" << std::endl;
    return failures ? 1 : 0;
}
