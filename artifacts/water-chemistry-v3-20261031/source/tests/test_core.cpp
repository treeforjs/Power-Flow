// test_core.cpp -- CPU qualification of chemistry v3.
// usage: test_core BLOB CHECKS_JSON OUTPUT_DIR [events_target]
// Checks (each prints PASS/FAIL and is recorded in OUTPUT_DIR/core-results.json):
//  math      deterministic log/log1p/sincos against the host libm (ulp bounds), u01 range
//  rng       Philox known answers (via the arrival core), layout injectivity probes, domain refusals
//  data      blob evaluation equals every tabulated source breakpoint bit for bit; midpoint
//            linearity; majorant bounds sigma_tot v at dense probes
//  events    >= events_target real events from random carriers/cells: per-event charge and nuclei
//            conservation, per-event energy closure within 2^-40 of the event energy scale, exact
//            ledger identity in = out + q + internal + defect, statuses
//  partition record digest and exact ledger identical for 1, 2, 5, 8 partitions and reversed order
//  debit     tracked-target refusal order and truncation; quantized neutral release conservation
#include "chem_host.hpp"
#include "test_common.hpp"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <random>
#include <sstream>

using namespace pfc;
using namespace pfc_test;

namespace {
int failures = 0;
std::ostringstream results;
bool first_result = true;
void check(bool ok, const std::string& name, const std::string& detail = "") {
    std::cout << (ok ? "PASS " : "FAIL ") << name << (detail.empty() ? "" : "  " + detail) << "\n";
    if (!ok) ++failures;
    results << (first_result ? "" : ",\n") << "  {\"check\": \"" << name << "\", \"pass\": " << (ok ? "true" : "false") << ", \"detail\": \"" << detail << "\"}";
    first_result = false;
}
std::string fmt(double x) {
    char b[64];
    std::snprintf(b, sizeof b, "%.6g", x);
    return b;
}
double ulp_err(double a, double ref) {
    if (a == ref) return 0;
    const double u = std::nextafter(std::fabs(ref), 1e308) - std::fabs(ref);
    return std::fabs(a - ref) / u;
}

// ---------------------------------------------------------------- math
void test_math() {
    std::uint64_t s = 12345;
    double worst_log = 0, worst_log1p = 0, worst_sc = 0;
    for (int i = 0; i < 2000000; ++i) {
        double x;
        switch (i % 4) {
            case 0: x = u01(s); break;
            case 1: x = std::ldexp(1.0 + u01(s), int(splitmix(s) % 2000) - 1000); break;
            case 2: x = 1.0 + (u01(s) - 0.5) * 1e-6; break;
            default: x = std::bit_cast<double>(splitmix(s) & 0x000FFFFFFFFFFFFFull) + 1e-320; break;
        }
        if (!(x > 0)) continue;
        worst_log = std::max(worst_log, ulp_err(pfc::log(x), std::log(x)));
        const double y = (u01(s) - 0.5) * (i % 3 == 0 ? 1e-8 : 1.9);
        if (y > -1) worst_log1p = std::max(worst_log1p, ulp_err(pfc::log1p(y), std::log1p(y)));
        double c, sn;
        const double u = u01(s);
        pfc::sincos_2pi(u, &c, &sn);
        const long double a = 2.0L * 3.14159265358979323846264338327950288L * (long double)u;
        worst_sc = std::max(worst_sc, double(std::fabs((long double)c - std::cos(a)) + std::fabs((long double)sn - std::sin(a))));
    }
    check(worst_log <= 1.0, "math/log ulp", "max ulp " + fmt(worst_log));
    check(worst_log1p <= 4.0, "math/log1p ulp", "max ulp " + fmt(worst_log1p));
    check(worst_sc <= 4e-16, "math/sincos_2pi abs", "max abs error " + fmt(worst_sc));
    check(pfc::log(1.0) == 0.0 && pfc::log(0.0) == -INFINITY && std::isnan(pfc::log(-1.0)), "math/log special values");
}

// ---------------------------------------------------------------- rng
void test_rng() {
    // Known-answer of the shared Philox core (Random123 kat_vectors).
    pfa::Block c{{0, 0, 0, 0}};
    pfa::Key2 k{{0, 0}};
    const auto r = pfa::philox4x32_10(c, k);
    check(r.v[0] == 0x6627e8d5u && r.v[1] == 0xe169c58du && r.v[2] == 0xbc57ac4cu && r.v[3] == 0x9b00dbd8u, "rng/philox known answer");
    // Domain refusals.
    ChemKey ok{1, 5, 7, PFC_PURPOSE_MCC};
    ChemKey bad_id{1, 5, PFC_CARRIER_LIMIT, PFC_PURPOSE_MCC};
    ChemKey bad_step{1, PFA_STEP_LIMIT, 7, PFC_PURPOSE_MCC};
    ChemKey bad_purpose{1, 5, 7, 2u};
    check(chem_key_valid(ok, 0) && !chem_key_valid(ok, PFC_EVENT_LIMIT) && !chem_key_valid(bad_id, 0) && !chem_key_valid(bad_step, 0) &&
              !chem_key_valid(bad_purpose, 0),
          "rng/domain refusals");
    // Counter injectivity probes: distinct tuples give distinct first words (sampled).
    std::map<pfa_u64, int> seen;
    int collisions = 0;
    std::uint64_t st = 99;
    for (int i = 0; i < 200000; ++i) {
        ChemKey key{0xC0FFEE, splitmix(st) % 1000, splitmix(st) % (1ull << 48), (i & 1) ? PFC_PURPOSE_MCC : PFC_PURPOSE_RELEASE_COUNT};
        const pfa_u64 w = chem_word(key, pfc_u32(splitmix(st) % 4096), pfc_u32(splitmix(st) % 512));
        collisions += seen.count(w);
        seen[w] = 1;
    }
    check(collisions == 0, "rng/no word collisions in 2e5 random keys");
    // Arrival purposes never alias chemistry purposes: the c2 fields differ in the purpose byte.
    pfa::StreamKey ak{1, 5, 7 & 0xFFFF, 7, 1};
    const auto ab = pfa::stream_block(ak, 0);
    const auto cb = chem_block(ChemKey{1, 5, 7, PFC_PURPOSE_MCC}, 0, 0);
    check(!(ab.v[0] == cb.v[0] && ab.v[1] == cb.v[1] && ab.v[2] == cb.v[2] && ab.v[3] == cb.v[3]), "rng/arrival and chemistry streams differ");
    // u01 mean and range over a stream.
    double sum = 0, mn = 1, mx = 0;
    const int n = 1 << 20;
    for (int i = 0; i < n; ++i) {
        const double u = chem_u01(ok, pfc_u32(i >> 9), pfc_u32(i & 511));
        sum += u;
        mn = std::min(mn, u);
        mx = std::max(mx, u);
    }
    check(mn > 0 && mx < 1 && std::fabs(sum / n - 0.5) < 5.0 * std::sqrt(1.0 / 12 / n), "rng/u01 open range and mean", "mean " + fmt(sum / n));
}

// ---------------------------------------------------------------- data
struct CheckPoint {
    pfc_u32 process;
    double E, sigma;
};
std::vector<CheckPoint> read_checks(const std::string& path, std::string* sha) {
    // minimal parser for {"blob_sha256": "...", "checks": [{"process": n, "E": x, "sigma": y, "kind": "..."}, ...]}
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    const std::string t = ss.str();
    auto p = t.find("\"blob_sha256\"");
    p = t.find('"', t.find(':', p) + 1);
    *sha = t.substr(p + 1, 64);
    std::vector<CheckPoint> out;
    std::size_t pos = 0;
    while ((pos = t.find("\"process\":", pos)) != std::string::npos) {
        CheckPoint c{};
        c.process = pfc_u32(std::stoul(t.substr(pos + 10)));
        const auto pe = t.find("\"E\":", pos);
        c.E = std::strtod(t.c_str() + pe + 4, nullptr);
        const auto ps = t.find("\"sigma\":", pos);
        c.sigma = std::strtod(t.c_str() + ps + 8, nullptr);
        out.push_back(c);
        pos = ps;
    }
    return out;
}
double speed(double E, double m) { return PFC_C_M_S * std::sqrt(E * (E + 2 * m)) / (E + m); }

void test_data(const Blob& b, const std::string& checks_path) {
    std::string sha;
    const auto cps = read_checks(checks_path, &sha);
    check(sha == b.sha256, "data/checks bind the blob", sha.substr(0, 16));
    std::size_t bad = 0;
    for (const auto& c : cps) {
        const ProcessRec& pr = b.view.procs[c.process];
        const GroupRec& gr = b.view.groups[pr.group];
        const int gi = grid_locate(b.view.grid + gr.grid_off, gr.grid_n, c.E);
        const double v = sigma_eval(b.view, gr, c.process - gr.proc_off, gi, c.E);
        if (std::bit_cast<std::uint64_t>(v) != std::bit_cast<std::uint64_t>(c.sigma)) ++bad;
    }
    check(bad == 0 && !cps.empty(), "data/source breakpoints reproduced bit for bit", std::to_string(cps.size()) + " points, " + std::to_string(bad) + " differ");
    // Majorant: dense probes of sigma_tot v must not exceed the rmq bound of any interval containing them.
    std::uint64_t s = 777;
    std::size_t viol = 0, probes = 0;
    for (pfc_u32 g = 0; g < b.n_groups; ++g) {
        const GroupRec& gr = b.view.groups[g];
        const double* grid = b.view.grid + gr.grid_off;
        const double m = b.view.species[gr.projectile].mass_eV;
        for (int k = 0; k < 200000; ++k) {
            const double lo = std::log(grid[0]), hi = std::log(grid[gr.grid_n - 1]);
            const double E = std::exp(lo + (hi - lo) * u01(s));
            const int gi = grid_locate(grid, gr.grid_n, E);
            double F = 0;
            for (pfc_u32 j = 0; j < gr.proc_n; ++j) F += sigma_eval(b.view, gr, j, gi, E) * speed(E, m);
            const double M = majorant(b.view, gr, E, E);
            ++probes;
            if (F > M) ++viol;
            // also a random enclosing range
            const double E2 = E * (1 + u01(s));
            if (F > majorant(b.view, gr, E, E2)) ++viol;
        }
    }
    check(viol == 0, "data/majorant bounds sigma_tot v", std::to_string(probes) + " probes, " + std::to_string(viol) + " violations");
}

void test_events(const Blob& b, const std::string& outdir, std::uint64_t target) {
    StepConfig cfg{};
    cfg.seed = 0x5EED5EED12345678ull;
    cfg.step = 4242;
    cfg.dt_s = 4e-11;
    cfg.trunc_sigmas = 7.0;
    cfg.purpose = PFC_PURPOSE_MCC;
    cfg.max_candidates = 4096;
    cfg.max_events = 64;
    std::uint64_t total_events = 0, conserve_bad = 0, carriers = 0;
    Ledger total;
    total.init(b);
    std::map<pfc_u32, std::uint64_t> status;
    std::vector<std::uint64_t> per_proc(b.n_proc, 0);
    int round = 0;
    const auto t0 = std::chrono::steady_clock::now();
    while (total_events < target && round < 400) {
        const Population P = make_population(b, 20000, 97, 1000 + round);
        cfg.step = 4242 + round;
        StepResult r;
        cpu_step(b, cfg, P.carriers, P.cells, r);
        for (std::size_t i = 0; i < r.outs.size(); ++i) {
            status[r.outs[i].status]++;
            ledger_add_carrier(total, r.outs[i]);
            for (pfc_u32 k = 0; k < r.outs[i].n_events; ++k) {
                const EventRec& e = r.events[i * cfg.max_events + k];
                if (!event_conserves(b, e)) ++conserve_bad;
                ledger_add_event(b, total, e);
                per_proc[e.process]++;
            }
            total_events += r.outs[i].n_events;
        }
        carriers += P.carriers.size();
        ++round;
    }
    ledger_finish(total);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    check(total_events >= target, "events/sampled count", std::to_string(total_events) + " events from " + std::to_string(carriers) + " carriers in " + fmt(secs) + " s");
    check(conserve_bad == 0, "events/charge and nuclei conserved per event", std::to_string(conserve_bad) + " violations");
    check(total.closure_violations == 0, "events/energy closure per event within 2^-40 of scale",
          std::to_string(total.closure_violations) + " violations, max relative defect " + fmt(total.max_rel_defect));
    // exact ledger identity: recompute in - out - q - internal - defect == 0 exactly
    ExactSum z = total.ke_in_projectile;
    z.add(total.ke_in_target);
    for (auto* x : {&total.ke_out_survivor, &total.ke_out_birth, &total.ke_out_reservoir, &total.ke_out_gas, &total.q_reference, &total.defect}) z.sub(*x);
    for (auto& x : total.internal) z.sub(x);
    check(z.is_zero(), "events/exact ledger identity (in = out + q + internal + defect)", "defect total " + fmt(total.defect.approx()) + " eV of input " + fmt(total.ke_in_projectile.approx()) + " eV");
    ExactSum dq = total.charge_created;
    dq.sub(total.charge_destroyed);
    bool nuc_ok = true;
    for (int k = 0; k < 6; ++k) nuc_ok = nuc_ok && total.nuclei_created[k] == total.nuclei_destroyed[k];
    check(dq.is_zero() && nuc_ok, "events/exact charge and nuclei totals balance");
    bool all_proc = true;
    std::string missing;
    for (pfc_u32 p = 0; p < b.n_proc; ++p)
        if (!per_proc[p]) {
            all_proc = false;
            missing += b.process_names[p] + "; ";
        }
    std::cout << "INFO events/processes not sampled by the random population (covered by process/targeted): " << missing << "\n";
    (void)all_proc;
    std::ostringstream st;
    for (auto& [k, v] : status) st << k << ":" << v << " ";
    const bool status_ok = status.count(PFC_ST_MAJORANT) == 0 && status.count(PFC_ST_BAD_INPUT) == 0 && status.count(PFC_ST_KEY_DOMAIN) == 0;
    check(status_ok, "events/no majorant, input or domain refusals", st.str());
    // write per-process counts
    std::ofstream f(outdir + "/events-per-process.json");
    f << "{\n";
    for (pfc_u32 p = 0; p < b.n_proc; ++p) f << "  \"" << b.process_names[p] << "\": " << per_proc[p] << (p + 1 < b.n_proc ? ",\n" : "\n");
    f << "}\n";
    std::ofstream g(outdir + "/ledger-1e6.json");
    g << "{\"events\": " << total.events << ", \"ke_in_projectile_eV\": " << total.ke_in_projectile.approx() << ", \"ke_in_target_eV\": " << total.ke_in_target.approx()
      << ", \"ke_out_survivor_eV\": " << total.ke_out_survivor.approx() << ", \"ke_out_birth_eV\": " << total.ke_out_birth.approx()
      << ", \"ke_out_reservoir_eV\": " << total.ke_out_reservoir.approx() << ", \"ke_out_gas_eV\": " << total.ke_out_gas.approx()
      << ", \"q_reference_eV\": " << total.q_reference.approx() << ", \"internal_eV\": [";
    for (int k = 0; k < int(PFC_SINK_COUNT); ++k) g << (k ? ", " : "") << total.internal[k].approx();
    g << "], \"defect_eV\": " << total.defect.approx() << ", \"defect_exact_hex\": \"" << total.defect.hex() << "\", \"max_relative_event_defect\": " << total.max_rel_defect
      << ", \"ledger_sha256\": \"" << ledger_digest(total) << "\"}\n";
}

// Per-process kinematics: every process with a positive cross section is evaluated 2000 times at the
// energy of its largest sigma*v through the shared event_kinematics (target velocities from a 300 K
// Maxwellian), and every event must conserve charge, nuclei and energy (closure within 2^-40 of scale).
// Every process is additionally sampled through the full MCC step with its group restricted to that
// single target (one cell), which checks that real events of the process occur.
void test_processes(const Blob& b, const std::string& outdir) {
    std::size_t bad = 0, tested = 0, placeholders = 0, failed_kin = 0;
    std::ofstream f(outdir + "/process-targeted.json");
    f << "{\n";
    for (pfc_u32 p = 0; p < b.n_proc; ++p) {
        const ProcessRec& pr = b.view.procs[p];
        const GroupRec& gr = b.view.groups[pr.group];
        const pfc_u32 local = p - gr.proc_off;
        const double* g = b.view.grid + gr.grid_off;
        const double m = b.view.species[gr.projectile].mass_eV, mt = b.view.species[gr.target].mass_eV;
        int best = -1;
        double bv = 0;
        for (pfc_u32 i = 0; i < gr.grid_n; ++i) {
            const double sv = b.view.sigma[gr.sigma_off + std::size_t(local) * gr.grid_n + i] * speed(g[i], m);
            if (sv > bv) {
                bv = sv;
                best = int(i);
            }
        }
        if (best < 0) {
            ++placeholders;
            f << (p ? ",\n" : "") << "  \"" << b.process_names[p] << "\": {\"placeholder\": true}";
            continue;
        }
        ++tested;
        double E = g[best];
        if (best + 1 < int(gr.grid_n) && b.view.sigma[gr.sigma_off + std::size_t(local) * gr.grid_n + best + 1] > 0) E = 0.5 * (g[best] + g[best + 1]);
        const double st = std::sqrt(8.617333262e-5 * 300.0 / mt) * PFC_C_M_S;
        std::uint64_t s = 555 + p, ok = 0, closed = 0, viol = 0;
        double max_rel = 0;
        for (int k = 0; k < 2000; ++k) {
            const V3 u = u_from_p(m, scale(iso_dir(u01(s), u01(s)), std::sqrt(E * (E + 2 * m))));
            const V3 ut = v3(st * std::sqrt(-2 * std::log(u01(s))) * std::cos(6.283185307179586 * u01(s)), st * (u01(s) - 0.5), st * (u01(s) - 0.5));
            const V3 wt = w_from_u(ut);
            const double gt = std::sqrt(1 + norm2(wt));
            const V3 pl = p_from_u(m, u);
            const double K = ke_from_u(m, u);
            const V3 pt = boost_in(wt, gt, pl, m + K);
            const double T = kinetic_from_p2(m, norm2(pt));
            const int gi = grid_locate(g, gr.grid_n, T);
            ChemKey key{0x7E57, 9, pfa_u64(1 + k), PFC_PURPOSE_FIXTURE_TEST};
            V3 plab[PFC_MAX_OUT];
            double mk[PFC_MAX_OUT], eint;
            if (!event_kinematics(b.view, pr, gr, gi, key, pfc_u32(p), m, mt, pt, T, wt, gt, plab, mk, &eint)) {
                ++closed;
                continue;
            }
            EventRec e{};
            e.process = p;
            e.target_species = gr.target;
            e.n_out = pr.n_out;
            e.weight = 1e6;
            e.ke_in_p = K;
            e.ke_in_t = ke_from_u(mt, ut);
            e.q_eV = pr.q_eV;
            e.e_internal = eint;
            for (pfc_u32 i = 0; i < pr.n_out; ++i) {
                const V3 uo = u_from_p(mk[i], plab[i]);
                e.out_species[i] = pr.out_species[i];
                e.out_dest[i] = pr.out_dest[i];
                e.ke_out[i] = ke_from_u(b.view.species[pr.out_species[i]].mass_eV, uo);
            }
            ++ok;
            const double d = event_defect(e);
            max_rel = std::max(max_rel, std::fabs(d) / event_tolerance(e) * 0x1p-40);
            if (!event_conserves(b, e) || !(std::fabs(d) <= event_tolerance(e))) ++viol;
        }
        if (viol || ok < 1990) ++bad;
        if (closed) ++failed_kin;
        f << (p ? ",\n" : "") << "  \"" << b.process_names[p] << "\": {\"E_eV\": " << E << ", \"events\": " << ok << ", \"kinematically_closed\": " << closed
          << ", \"violations\": " << viol << ", \"max_relative_defect\": " << max_rel << "}";
    }
    f << "\n}\n";
    check(bad == 0, "process/every process with data: 2000 conserving events at its peak sigma*v energy",
          std::to_string(tested) + " processes, " + std::to_string(bad) + " failing, " + std::to_string(failed_kin) + " with kinematically closed draws, " +
              std::to_string(placeholders) + " no-data placeholders");
}

void test_partition(const Blob& b) {
    StepConfig cfg{};
    cfg.seed = 0xABCDEF0123456789ull;
    cfg.step = 77;
    cfg.dt_s = 5e-12;
    cfg.trunc_sigmas = 7.0;
    cfg.purpose = PFC_PURPOSE_MCC;
    cfg.max_candidates = 4096;
    cfg.max_events = 64;
    const Population P = make_population(b, 30000, 31, 4711);
    StepResult ref;
    cpu_step(b, cfg, P.carriers, P.cells, ref);
    const std::string dref = step_digest(ref.outs, ref.events, cfg.max_events, P.carriers);
    const std::string lref = ledger_digest(ledger_of(b, ref, cfg.max_events));
    bool all = true;
    std::string detail;
    for (int parts : {1, 2, 5, 8}) {
        // each partition is a separate contiguous slice of a permuted carrier list, run in reverse slice order
        std::vector<std::size_t> perm(P.carriers.size());
        for (std::size_t i = 0; i < perm.size(); ++i) perm[i] = i;
        std::uint64_t s = 1000 + parts;
        for (std::size_t i = perm.size(); i > 1; --i) std::swap(perm[i - 1], perm[splitmix(s) % i]);
        std::vector<CarrierOut> outs(P.carriers.size());
        std::vector<EventRec> evs(P.carriers.size() * cfg.max_events);
        for (int q = parts - 1; q >= 0; --q) {
            const std::size_t lo = perm.size() * q / parts, hi = perm.size() * (q + 1) / parts;
            std::vector<CarrierIn> slice;
            for (std::size_t i = lo; i < hi; ++i) slice.push_back(P.carriers[perm[i]]);
            StepResult r;
            cpu_step(b, cfg, slice, P.cells, r);
            for (std::size_t i = lo; i < hi; ++i) {
                outs[perm[i]] = r.outs[i - lo];
                for (pfc_u32 k = 0; k < cfg.max_events; ++k) evs[perm[i] * cfg.max_events + k] = r.events[(i - lo) * cfg.max_events + k];
            }
        }
        StepResult merged;
        merged.outs = outs;
        merged.events = evs;
        const std::string d = step_digest(outs, evs, cfg.max_events, P.carriers);
        const std::string l = ledger_digest(ledger_of(b, merged, cfg.max_events));
        all = all && d == dref && l == lref;
        detail += std::to_string(parts) + (d == dref && l == lref ? ":same " : ":DIFF ");
    }
    check(all, "partition/records and exact ledger identical for 1,2,5,8 partitions", detail + "digest " + dref.substr(0, 16));
}

void test_debit_release(const Blob& b) {
    StepConfig cfg{};
    cfg.seed = 99;
    cfg.step = 3;
    cfg.dt_s = 2e-11;
    cfg.trunc_sigmas = 7.0;
    cfg.purpose = PFC_PURPOSE_MCC;
    cfg.max_candidates = 4096;
    cfg.max_events = 64;
    Population P = make_population(b, 5000, 3, 1234);
    StepResult r;
    cpu_step(b, cfg, P.carriers, P.cells, r);
    // tracked inventory of H2O in cell 0 limited to a few events' weight
    const pfc_u32 h2o = species_index(b, "H2O");
    std::map<InventoryKey, ExactSum> inv;
    ExactSum stock;
    stock.add(3.5e7);
    inv[InventoryKey{0, h2o}] = stock;
    const DebitResult d = resolve_debits(b, P.carriers, r, cfg.max_events, inv);
    // debited <= stock, and a refusal occurred
    ExactSum left = stock;
    left.sub(d.debited.count(InventoryKey{0, h2o}) ? d.debited.at(InventoryKey{0, h2o}) : ExactSum{});
    check(!left.negative() && d.refused > 0, "debit/tracked inventory never overdrawn; refusals recorded", std::to_string(d.refused) + " refused");
    // truncation: admitted flags form a prefix per carrier
    bool prefix = true;
    for (std::size_t i = 0; i < P.carriers.size(); ++i)
        for (pfc_u32 k = 0; k < r.outs[i].n_events; ++k) prefix = prefix && (d.admitted[i * cfg.max_events + k] == (k < d.admitted_events[i] ? 1 : 0));
    check(prefix, "debit/admitted events form a per-carrier prefix");
    // reservoir: deposit admitted events, release with w0
    std::map<InventoryKey, Reservoir> res;
    ExactSum deposited;
    for (std::size_t i = 0; i < P.carriers.size(); ++i)
        for (pfc_u32 k = 0; k < d.admitted_events[i]; ++k) {
            const EventRec& e = r.events[i * cfg.max_events + k];
            deposit_reservoir(b, res, e);
            for (pfc_u32 j = 0; j < e.n_out; ++j)
                if (e.out_dest[j] == PFC_DEST_RESERVOIR) deposited.add(e.weight);
        }
    ExactSum before;
    for (auto& [k, R] : res) before.add(R.count);
    const auto rel = release_reservoirs(b, res, 2.5e6, 99, 3);
    ExactSum after, released;
    for (auto& [k, R] : res) {
        after.add(R.count);
        check(!R.count.negative(), "release/stock never negative " + std::to_string(k.cell) + ":" + b.species_names[k.species]);
    }
    for (auto& p : rel) released.add(p.weight);
    ExactSum sum = after;
    sum.add(released);
    check(before == deposited && sum == before, "release/stock + released = deposited exactly", std::to_string(rel.size()) + " packets");
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "usage: test_core BLOB CHECKS_JSON OUTPUT_DIR [events_target]\n";
        return 2;
    }
    const std::string out = argv[3];
    std::filesystem::create_directories(out);
    const std::uint64_t target = argc > 4 ? std::stoull(argv[4]) : 1000000ull;
    try {
        const Blob b = load_blob(argv[1]);
        std::cout << "blob " << b.sha256 << " species " << b.n_species << " groups " << b.n_groups << " processes " << b.n_proc << " regions " << b.n_regions << "\n";
        test_math();
        test_rng();
        test_data(b, argv[2]);
        test_events(b, out, target);
        test_processes(b, out);
        test_partition(b);
        test_debit_release(b);
        std::ofstream f(out + "/core-results.json");
        f << "{\"blob_sha256\": \"" << b.sha256 << "\", \"failures\": " << failures << ", \"results\": [\n" << results.str() << "\n]}\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR " << e.what() << "\n";
        return 1;
    }
    std::cout << (failures ? "FAILED " : "ALL PASS ") << failures << "\n";
    return failures ? 1 : 0;
}
