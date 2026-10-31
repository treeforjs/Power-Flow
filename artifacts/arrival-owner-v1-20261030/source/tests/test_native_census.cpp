// How staged births enter the native v5 census (native candidate005, source identity
// cca72822..., library build/libpic_native_v5.a, read only). Native has no birth API
// yet (plan item H2); its only carrier-creating routes are genesis (create_full_time)
// and reaction products, both of which draw IDs from one monotone counter in the
// order packets are supplied, one native ID per positive carrier. This test supplies
// one step's births, in the staged canonical order, as authored packets after the
// pre-existing packets of a genesis and checks:
//  * each birth leaf is one positive carrier, hence exactly one native ID;
//  * native assigns ID = base + rank in canonical order (base = 1 + carriers before);
//  * the census (ascending native ID) lists births in canonical order with the
//    logical packet ID, species state and weight bits of the staged birth;
//  * the receipt built from the census is accepted by Owner::commit, while a
//    census built from a shuffled order is refused (native-receipt).
// usage: test_native_census OUTPUT_DIRECTORY
#include "raikou/mfem_pic_full_time.h"
#include "pf_arrival/arrival_owner.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>

namespace m = raikou::mfem_backend;
namespace pa = pf_arrival;

namespace {
long checks = 0;
void check(bool ok, const std::string& what) {
    ++checks;
    if (!ok) throw std::runtime_error("check failed: " + what);
}
m::GapPicOptions options() {   // the native v5 storage-control fixture geometry
    m::GapPicOptions o;
    o.nx = o.ny = 45;
    o.width_m = o.gap_m = .009;
    o.depth_m = .025;
    o.finite_foils = true;
    o.finite_axial_extent = true;
    o.short_z_m = 0;
    o.feed_z_m = .025;
    o.lower_foil = {.002, .007, .0018, .002};
    o.upper_foil = {.002, .007, .007, .0072};
    return o;
}
std::vector<m::ArenaStateSpec> states() { return {{"e", 0, 0, -1, m::pic_me, 0, 0}, {"n", 2, 1, 0, 2.9915e-26, 0, 0}}; }
std::string digest(const std::string& s) { return pa::sha256_hex(s); }
std::vector<m::ArenaChannel> channels() {
    return {{"conditional-storage-control", "e", "n", {"e", "n"}, digest("data"), digest("kernel"), digest("binding"), "manufactured-storage-control-only"}};
}
m::FullTimeGenesisBinding genesis() {
    m::FullTimeGenesisBinding g;
    g.time_s = 26.31715e-9;
    g.origin_provenance = "manufactured/arrival-owner-census-check-no-simulation";
    g.field_sha256 = digest("field");
    g.mesh_sha256 = digest("mesh");
    g.material_sha256 = digest("material");
    g.frame_sha256 = digest("frame");
    g.common_clock_sha256 = digest("initial-clock");
    g.heat_sha256 = digest("heat");
    g.source_sha256 = digest("source");
    g.rng_sha256 = digest("chemistry-RNG-and-config");
    return g;
}
std::shared_ptr<const pa::Policy> policy() {
    pa::Configuration c;
    c.run_seed = 9001;
    c.streams = {{0, pa::Kind::electron, 0, "e", 1e5, 4096}, {3, pa::Kind::electron, 1, "e", 2.5e4, 4096}, {8, pa::Kind::neutral, 0, "n", 1e7, 4096},
                 {11, pa::Kind::neutral, 1, "n", 3.3e6, 4096}};
    return std::make_shared<const pa::Policy>(c);
}
pa::StepInput input(std::uint64_t step) {
    const double t = 26.31715e-9, dt = 1e-14;
    pa::StepInput in;
    in.step = step;
    in.electrons = {{0, 4.2e6, std::log(4.2e6), true, t, dt, t + dt}, {3, 7.77e5, std::log(7.77e5), true, t, dt, t + dt}};
    in.neutrals = {{8, pa::Count(3.21e8), std::log(3.21e8), pa::Count(1e20), t, t + dt}, {11, pa::Count(9.9e7), std::log(9.9e7), pa::Count(1e20), t, t + dt}};
    return in;
}
// Authored packets: two pre-existing packets (count 4 and 8: one carrier each),
// then the births in the given order at distinct interior positions.
std::vector<m::FullTimeAuthoredPacket> packets(const std::vector<pa::Birth>& births, const pa::Policy& P) {
    std::vector<m::FullTimeAuthoredPacket> out = {{"electron", "e", "manufactured/source", {4, 0}, {.004, .004, .0125}, {1e6, 0, 0}},
                                                  {"neutral", "n", "manufactured/source", {8, 0}, {.00403, .004, .0125}, {0, 0, 0}}};
    for (const auto& b : births) {
        const double dx = 1e-6 * double(b.index % 1000), dy = 1e-6 * double(b.stream_id);
        out.push_back({b.packet_id(), P.stream(b.stream_id).species_state, "arrival-v1/birth", pa::Count(b.weight), {.0041 + dx, .0042 + dy, .0125},
                       {b.kind == pa::Kind::electron ? 2e5 : 300., 0, 0}});
    }
    return out;
}
pa::NativeBirthReceipt receipt_from_census(const m::PicArena& arena) {
    auto census = arena.rate_population(arena.full_time_carrier_cap());
    std::sort(census.begin(), census.end(), [](const auto& a, const auto& b) { return a.native.id < b.native.id; });
    pa::NativeBirthReceipt r;
    std::uint64_t last = 0;
    for (const auto& e : census) {
        last = std::max(last, e.native.id);
        if (e.packet_id.rfind("arrival-v1/", 0) != 0) continue;
        if (r.entries.empty()) r.base = e.native.id;
        r.entries.push_back({e.packet_id, e.state, e.native.id, e.native.weight});
    }
    r.next_id_after = last + 1;   // genesis: next native ID follows the last carrier
    return r;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: test_native_census OUTPUT_DIRECTORY\n";
        return 2;
    }
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
    const auto t0 = std::chrono::steady_clock::now();
    std::string failure, census_digest;
    std::size_t births_n = 0;
    try {
        auto P = policy();
        const std::uint64_t step = 41;
        // Canonical order: the staged births with base = 3 (IDs 1 and 2 are the pre-existing carriers).
        pa::Owner owner(P, step);
        auto staged = owner.stage(owner.plan(input(step)), 3);
        for (const auto& d : staged.plan().decisions())
            if (d.kind == pa::Kind::neutral) {
                std::vector<double> t;
                for (const auto& b : staged.births())
                    if (b.stream_id == d.stream_id) t.push_back(d.begin_s + b.marks[5] * (d.end_s - d.begin_s));
                staged.bind_source_times(d.stream_id, t);
            }
        births_n = staged.births().size();
        check(births_n > 100, "fixture has births of all four streams");
        for (const auto& b : staged.births()) check(pa::Count(b.weight).positive_carriers().size() == 1, "one positive carrier per birth leaf");
        auto arena = m::PicArena::create_full_time(options(), states(), channels(), packets(staged.births(), *P), genesis());
        const auto r = receipt_from_census(*arena);
        check(r.base == 3 && r.entries.size() == births_n && r.next_id_after == 3 + births_n, "native assigned one contiguous ID block");
        for (std::size_t k = 0; k < births_n; ++k) {
            const auto& b = staged.births()[k];
            const auto& e = r.entries[k];
            check(e.native_id == b.native_id && e.native_id == 3 + k, "native ID = base + canonical rank");
            check(e.packet_id == b.packet_id() && e.state == P->stream(b.stream_id).species_state, "census packet ID and species");
            check(std::bit_cast<std::uint64_t>(e.weight) == std::bit_cast<std::uint64_t>(b.weight), "census carrier weight bits");
        }
        // Logical packets: exact count equals the leaf, one native ID each.
        const auto snap = arena->snapshot();
        for (const auto& b : staged.births()) {
            const auto& pk = snap->packets.at(b.packet_id());
            check(pk.count == pa::Count(b.weight) && pk.native_ids.size() == 1 && pk.native_ids[0] == b.native_id, "logical packet maps to one native ID");
        }
        census_digest = arena->state_prefix();
        check(owner.commit(std::move(staged), r).chain_head == owner.ledger_head(), "receipt from the native census accepted");
        // Shuffled supply order: native numbers the births differently; the receipt is refused.
        pa::Owner other(P, step);
        auto staged2 = other.stage(other.plan(input(step)), 3);
        for (const auto& d : staged2.plan().decisions())
            if (d.kind == pa::Kind::neutral) {
                std::vector<double> t;
                for (const auto& b : staged2.births())
                    if (b.stream_id == d.stream_id) t.push_back(d.begin_s);
                staged2.bind_source_times(d.stream_id, t);
            }
        auto shuffled = staged2.births();
        std::mt19937_64 g(5);
        std::shuffle(shuffled.begin(), shuffled.end(), g);
        auto arena2 = m::PicArena::create_full_time(options(), states(), channels(), packets(shuffled, *P), genesis());
        const auto r2 = receipt_from_census(*arena2);
        const auto before = other.state_sha256();
        bool refused = false;
        try {
            (void)other.commit(std::move(staged2), r2);
        } catch (const pa::Refusal& e) {
            refused = e.kind == pa::RefusalKind::native_receipt;
        }
        check(refused && other.state_sha256() == before, "non-canonical native numbering refused, owner unchanged");
    } catch (const std::exception& e) {
        failure = e.what();
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::ofstream r(out / "native-results.json");
    r << "{\"program\":\"test_native_census\",\"checks\":" << checks << ",\"passed\":" << (failure.empty() ? "true" : "false") << ",\"failure\":\"" << failure
      << "\",\"births\":" << births_n << ",\"native_state_prefix\":\"" << census_digest << "\",\"seconds\":" << seconds << "}\n";
    std::cout << "test_native_census checks=" << checks << " " << (failure.empty() ? "PASS" : "FAIL: " + failure) << "\n";
    return failure.empty() ? 0 : 1;
}
