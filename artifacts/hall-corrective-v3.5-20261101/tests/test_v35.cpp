// Hall corrective v3.5 qualification: ResidentOwner::restore_as_genesis (genesis at a vacuum pre-run checkpoint).
//   test_v35 genesis TRACE OUT NVRTC UUID V33_CHECKPOINT V34_CHECKPOINT
//       V33_CHECKPOINT: the sealed Hall v3.3 20 ns pre-run checkpoint (hall-v3.2-resident-checkpoint-v1);
//       V34_CHECKPOINT: the Hall v3.4 20 ns pre-run checkpoint (same physical state, divB gross ledger).
//       Word-for-word restore of the 20 ns state (field words, charge ledger, admitted Gauss ledger, slice thermal states)
//       at revision 0 with the clock kept and the provenance equal to the checkpoint's seals, into the pre-run Model and
//       into a Model with heat-only patch declarations, other carrier limits, another field-step bound and other native
//       foil planes; save/restore of the rebased owner; physical lockstep of the rebased owner against the plain restore
//       of the same checkpoint (vacuum steps) and against a fresh genesis (vacuum and coupled steps); refusals (flipped
//       field word, flipped heat word, non-vacuum checkpoint, altered or removed patch, other drive trace, other mesh,
//       transfers not admitted by GenesisTransfer); the host Owner's genesis of its own checkpoints (words, slice states,
//       provenance, save/restore, lockstep against the plain restore, transfers, refusals).
//   test_v35 genesis-run TRACE OUT NVRTC UUID CHECKPOINT / genesis-continue TRACE OUT NVRTC UUID REBASED REFERENCE
//       fresh-process restore of a rebased owner continuing bit-identically.
#include "resident_owner.hpp"
#include "gpu/gpu_backends.hpp"
#include "common/hall_gpu_common.hpp"
#include "hall-finite-foil-em-20261029/axial_particles.hpp"
namespace hc = hall_corrective_v3;
#include "motion_fixtures.hpp"
#include <fstream>
#include <iostream>
#include <set>

namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
namespace gp = hall_corrective_v3::gpu;
namespace fs = std::filesystem;
using J = hc::Json;

namespace {
struct Results {
    J cases = J::array();
    int failed = 0;
    void record(const std::string& group, const std::string& name, bool ok, J evidence = J::object()) {
        cases.push_back({{"group", group}, {"name", name}, {"passed", ok}, {"evidence", std::move(evidence)}});
        if (!ok) ++failed;
        std::cout << group << " | " << name << ": " << (ok ? "pass" : "FAIL") << std::endl;
        if (!ok) std::cout << "  evidence: " << cases.back()["evidence"].dump().substr(0, 4000) << std::endl;
    }
};
struct Boundary { std::string field_sha256, heat_sha256, source_sha256, mesh_sha256, material_sha256, common_clock_sha256; double time_s = 0; std::uint64_t accepted_step = 0; };
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }
bool same_bits(const hc::Vec& a, const hc::Vec& b) { return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * 8) == 0); }
bool same_words(const hc::WordVec& a, const hc::WordVec& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (bits(a[i].hi) != bits(b[i].hi) || bits(a[i].lo) != bits(b[i].lo)) return false;
    return true;
}
bool zero_or_empty(const hc::Vec& v) { for (double x : v) if (bits(x) != 0) return false; return true; }
hc::ChargeLedger declared(const hc::Model& m, const axial::State& s, hc::FreeCharge particles) {
    const auto& g = m.mesh(); hc::ChargeLedger q; q.particles = std::move(particles);
    q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
    q.conductor_declaration = hc::sha256("parity: declared DC conductor/PEC induced charge");
    q.exterior_declaration = hc::sha256("parity: declared DC exterior termination charge");
    const auto field = hc::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = field[k]; else if (g.outer_node[k]) q.exterior_C[k] = field[k]; }
    return q;
}
hall_driver::MotionBinding binding_of(const hc::ResidentOwner& res, const hc::Model& m, long step, double dt) {
    hall_driver::MotionBinding b; b.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = res.revision(); b.staged_arena_revision = res.revision() + 1; b.begin_s = res.time_s(); b.dt_s = dt; b.end_s = b.begin_s + dt;
    b.staged_arena_sha256 = hc::sha256("v35/staged/" + std::to_string(step)); b.mesh_sha256 = res.motion_mesh_sha256();
    b.material_sha256 = m.material_sha256(); b.species_sha256 = hc::sha256("v35/species");
    b.common_clock_sha256 = hc::sha256("v35/clock/" + std::to_string(step)); b.kinetic_policy_sha256 = res.motion_context().policy_sha256;
    return b;
}
struct Tables { std::vector<std::string> states; std::vector<std::uint64_t> keys, idx; };
Tables tables_of(const std::vector<n::ArenaRateEntry>& c) {
    Tables t; std::set<std::string> u; for (const auto& r : c) u.insert(r.state); t.states.assign(u.begin(), u.end());
    std::map<std::string, std::uint64_t> key; for (const auto& r : c) { auto [it, ins] = key.emplace(r.packet_id, r.native.id); if (!ins) it->second = std::min(it->second, r.native.id); }
    std::vector<const n::ArenaRateEntry*> rows; for (const auto& r : c) rows.push_back(&r);
    std::stable_sort(rows.begin(), rows.end(), [](auto a, auto b) { return a->native.id < b->native.id; });
    for (const auto* r : rows) { t.keys.push_back(key.at(r->packet_id)); t.idx.push_back(std::uint64_t(std::lower_bound(t.states.begin(), t.states.end(), r->state) - t.states.begin())); }
    return t;
}
gp::ResidentOptions options() { gp::ResidentOptions o; o.host_work = gp::HostWork::device; o.verify_samples = 16; o.sample_seed = 7; return o; }

// Models: the pre-run Model of test_v32 prerun (resident backend, default carrier limits, no declaration, dt_max 1e-13,
// native v7 forms) and variants.
struct Models {
    std::shared_ptr<gp::Device> dev; std::string desc; std::string trace;
    std::shared_ptr<const hc::Model> make(const std::string& trace_file, hc::MeshSpec spec, hc::CarrierLimits l, hc::ThermalDeclaration d, double dt) const {
        return std::make_shared<const hc::Model>(spec, trace_file, std::make_shared<hc::ResidentBackend>(desc), l, std::move(d), dt, hc::MotionForms{true, 1e-14});
    }
    std::shared_ptr<const hc::Model> prerun(double dt = 1e-13) const { return make(trace, {45, 4, 1, .009}, hc::CarrierLimits{}, {}, dt); }
};
// Heat-only patches on two foil slices (the coordinator's source-owner declarations: no desorption, zero inventory).
hc::ThermalDeclaration heat_only(const hc::Mesh& g, double end = 1) {
    hc::ThermalDeclaration d;
    for (int body = 0; body < 2; ++body) {
        hc::SlicePatches s; s.source_hash = hc::sha256("v35/source-owner/foil-" + std::to_string(body));
        hc::eh::Patch p; p.id = "v35-heat-only-" + std::to_string(body); p.source_hash = s.source_hash; p.face = body ? hc::eh::Face::LowY : hc::eh::Face::HighY;
        p.begin = 0; p.end = end;
        s.patches.push_back(p);
        d.slices[{body, g.kshort + 1}] = std::move(s);
    }
    return d;
}
hc::CarrierLimits coordinator_limits() { hc::CarrierLimits l; l.maximum_carriers = 2000000; l.memory_budget_bytes = l.maximum_carriers * 4096 + (512ull << 20); return l; }
// Other native foil planes (the coordinator's planes differ from the test fixture's in their binary64 values; here the
// outer face of each foil moved by one ulp): other motion mesh digest, same mesh.
n::GapPicOptions aligned(const hc::Mesh& g) {
    auto o = mf::geometry(g);
    o.lower_foil.bottom_m = std::nextafter(o.lower_foil.bottom_m, -1.0); o.upper_foil.top_m = std::nextafter(o.upper_foil.top_m, 1.0);
    return o;
}
// Checkpoint contents for the word-for-word comparison.
struct Saved { axial::State s; hc::ChargeLedger q; J j; };
Saved saved(const hc::Mesh& g, const fs::path& dir) {
    Saved c; std::ifstream in(dir / "checkpoint.json"); c.j = J::parse(in);
    c.s = hc::checkpoint_detail::field(c.j.at("field"), g); c.q = hc::checkpoint_detail::charge(c.j.at("charge"), g);
    c.q.gauss_admitted_C = hc::checkpoint_detail::ledger_words(c.j.at("admitted_gauss_C"), g.nn);
    if (c.j.contains("divb_gross_T_m")) c.q.divb_gross_T_m = hc::checkpoint_detail::divb_gross(c.j.at("divb_gross_T_m"), g.cells.size());
    return c;
}
J word_checks(const axial::State& s, const hc::ChargeLedger& q, const axial::State& r, const hc::ChargeLedger& p, bool divb) {
    J c = J::object();
    c["E"] = same_bits(s.E, r.E); c["B"] = same_bits(s.B, r.B); c["E_low"] = same_bits(s.E_low, r.E_low); c["B_low"] = same_bits(s.B_low, r.B_low);
    c["time"] = bits(s.time_s) == bits(r.time_s);
    c["free_net"] = same_words(q.particles.net_C, p.particles.net_C); c["free_gross"] = same_bits(q.particles.gross_C, p.particles.gross_C);
    c["conductor"] = same_words(q.conductor_C, p.conductor_C); c["exterior"] = same_words(q.exterior_C, p.exterior_C);
    c["declarations"] = q.conductor_declaration == p.conductor_declaration && q.exterior_declaration == p.exterior_declaration;
    c["admitted_gauss_ledger"] = same_words(q.gauss_admitted_C, p.gauss_admitted_C);
    c["residue"] = same_bits(q.residue_C, p.residue_C);
    if (divb) c["divb_gross_ledger"] = (q.divb_gross_T_m.empty() && p.divb_gross_T_m.empty()) || same_bits(q.divb_gross_T_m, p.divb_gross_T_m) ||
                                       (zero_or_empty(q.divb_gross_T_m) && zero_or_empty(p.divb_gross_T_m));
    return c;
}
bool all_true(const J& c) { for (auto& [k, v] : c.items()) if (!v.get<bool>()) return false; return true; }
// Slice thermal states: temperatures, cumulative energies, watermarks and time equal; inventory: the recorded patches'
// values followed by zeros (added heat-only patches); revisions as given.
J heat_checks(const hc::ThermalV3& a, const hc::ThermalV3& b, std::uint64_t rev_a, std::uint64_t rev_b) {
    const auto sa = a.states(), sb = b.states();
    bool T = sa.size() == sb.size(), cum = T, wm = T, time = T, rev = T, inv = T;
    for (std::size_t i = 0; T && i < sa.size(); ++i) {
        T = T && same_bits(sa[i].T, sb[i].T); wm = wm && sa[i].watermarks == sb[i].watermarks; time = time && bits(sa[i].time) == bits(sb[i].time);
        for (std::size_t k = 0; k < sa[i].cumulative.size(); ++k) cum = cum && bits(sa[i].cumulative[k]) == bits(sb[i].cumulative[k]);
        rev = rev && sa[i].revision == rev_a && sb[i].revision == rev_b;
        const auto& small = sa[i].inventory.size() <= sb[i].inventory.size() ? sa[i].inventory : sb[i].inventory;
        const auto& large = sa[i].inventory.size() <= sb[i].inventory.size() ? sb[i].inventory : sa[i].inventory;
        for (std::size_t k = 0; k < large.size(); ++k) {
            const auto x = k < small.size() ? small[k] : hc::eh::Count{};
            inv = inv && bits(large[k].hi) == bits(x.hi) && bits(large[k].lo) == bits(x.lo);
        }
    }
    return {{"temperatures", T}, {"cumulative", cum}, {"watermarks", wm}, {"time", time}, {"revisions", rev}, {"inventory", inv}};
}
std::string refusal(const std::function<void()>& f) {
    try { f(); } catch (const std::exception& e) { return e.what(); }
    return "";
}
bool refused_with(const std::string& r, const std::string& text) { return !r.empty() && r.find(text) != std::string::npos; }
void copy_dir(const fs::path& from, const fs::path& to) { fs::remove_all(to); fs::copy(from, to, fs::copy_options::recursive); }
// Vacuum step of one owner (pre-run protocol).
void vacuum_step(hc::ResidentOwner& r, double dt) {
    auto p = r.prepare_vacuum(dt);
    Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
    hc::need(r.preflight_vacuum_final(p, b), "vacuum preflight");
    r.publish_vacuum_prevalidated(std::move(p));
}
J state_of(const hc::ResidentOwner& r, axial::State& s, hc::ChargeLedger& q) { r.accepted_state(s, q); return J::object(); }


// Host Owner: vacuum steps and the genesis of its checkpoints (CPU PCG Model of the parity tests, dt_max 1e-14).
void host_vacuum_step(hc::Owner& r, double dt) {
    auto p = r.prepare_vacuum(dt);
    Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
    hc::need(r.preflight_vacuum_final(p, b), "host vacuum preflight");
    r.publish_vacuum_prevalidated(std::move(p));
}
void host_genesis(Results& R, const std::string& trace, const fs::path& out) {
    auto hm = std::make_shared<const hc::Model>(hc::MeshSpec{45, 4, 1, .009}, trace);
    const auto& g = hm->mesh(); const auto o = mf::geometry(g); const std::vector<n::ArenaRateEntry> none;
    const auto seed = hm->dc_seed();
    hc::Owner run(hm, seed, declared(*hm, seed, hc::charge_of_census(g, o, none)), o, none);
    for (int k = 0; k < 5; ++k) host_vacuum_step(run, 1e-14);
    const auto dir = out / "host-prerun-r5"; fs::remove_all(dir); run.save(dir);
    const auto ck = saved(g, dir);
    std::unique_ptr<hc::Owner> w, p;
    const auto r = refusal([&] { w = hc::Owner::restore_as_genesis(hm, o, none, dir); });
    R.record("host-genesis", "a host Owner vacuum checkpoint (revision 5) is admitted as a genesis", r.empty() && w, {{"refusal", r}});
    if (!w) return;
    {   const auto c = word_checks(w->field(), w->charge(), ck.s, ck.q, true);
        const auto h = heat_checks(w->thermal(), run.thermal(), 0, 5);
        const auto& pv = w->genesis_provenance();
        J pc{{"words", all_true(c)}, {"heat", all_true(h)}, {"revision_0", w->revision() == 0 && w->field().revision == 0},
             {"provenance", pv.field_sha256 == run.field_sha256() && pv.heat_sha256 == run.heat_sha256() && pv.revision == 5 &&
                            pv.configuration_sha256 == hm->configuration_sha256() && pv.checkpoint_schema == "hall-v3.5-host-checkpoint"},
             {"identity_differs", w->field_sha256() != run.field_sha256()}};
        R.record("host-genesis", "words, slice states, revision 0 and provenance", all_true(pc), {{"checks", pc}, {"words", c}, {"heat", h}}); }
    {   const auto d2 = out / "host-rebased"; fs::remove_all(d2); w->save(d2);
        std::unique_ptr<hc::Owner> back; const auto rr = refusal([&] { back = hc::Owner::restore(hm, o, none, d2); });
        R.record("host-genesis", "save and plain restore of the rebased host owner keep digests and provenance",
                 rr.empty() && back && back->field_sha256() == w->field_sha256() && back->heat_sha256() == w->heat_sha256() &&
                 back->genesis_provenance().sha256() == w->genesis_provenance().sha256(), {{"refusal", rr}}); }
    const auto rp = refusal([&] { p = hc::Owner::restore(hm, o, none, dir); });
    if (p) {
        bool ok = true; std::string why;
        const auto rs = refusal([&] {
            for (int k = 0; k < 20; ++k) {
                host_vacuum_step(*w, 1e-14); host_vacuum_step(*p, 1e-14);
                const auto c = word_checks(w->field(), w->charge(), p->field(), p->charge(), true);
                J h = heat_checks(w->thermal(), p->thermal(), w->revision(), p->revision()); h.erase("revisions");
                if ((!all_true(c) || !all_true(h)) && ok) { ok = false; why = c.dump() + h.dump(); }
            }
        });
        R.record("host-genesis", "20 vacuum steps: the rebased host owner equals the plain restore in every word and slice state", rs.empty() && ok,
                 {{"refusal", rs}, {"first_difference", why}});
    } else R.record("host-genesis", "plain restore for the lockstep", false, {{"refusal", rp}});
    // transfers into another host Model (heat-only patches, carrier limits, foil planes)
    {   auto cm = std::make_shared<const hc::Model>(hc::MeshSpec{45, 4, 1, .009}, trace, hc::make_cpu_pcg(), coordinator_limits(), heat_only(g), 1e-14);
        std::unique_ptr<hc::Owner> t; const auto rt = refusal([&] { t = hc::Owner::restore_as_genesis(cm, aligned(g), none, dir); });
        bool ok = rt.empty() && t;
        J ev{{"refusal", rt}};
        if (t) { const auto c = word_checks(t->field(), t->charge(), ck.s, ck.q, true); ev["words"] = c; ev["heat"] = heat_checks(t->thermal(), run.thermal(), 0, 5);
                 ok = ok && all_true(c) && all_true(ev["heat"]); ev["transfer_text"] = t->genesis_provenance().transfer_text; }
        R.record("host-genesis", "transfer into a host Model with heat-only patches, other carrier limits and foil planes", ok, ev); }
    // refusals: flipped word, non-vacuum
    {   const auto tw = out / "host-tamper"; copy_dir(dir, tw);
        {   std::ifstream in(tw / "checkpoint.json"); J j = J::parse(in); in.close();
            auto& e = j["field"]["B"][7]; e = std::bit_cast<double>(bits(e.get<double>()) ^ 1);
            std::ofstream(tw / "checkpoint.json", std::ios::trunc) << j.dump() << '\n'; }
        const auto r1 = refusal([&] { (void)hc::Owner::restore_as_genesis(hm, o, none, tw); });
        R.record("host-genesis", "a flipped field word is refused by the seal", refused_with(r1, "checkpoint field/charge seal mismatch"), {{"refusal", r1}});
        auto census = mf::pairs(40, 5, 1e5);
        hc::Owner np(hm, seed, declared(*hm, seed, hc::charge_of_census(g, o, census)), o, census);
        const auto nd = out / "host-non-vacuum"; fs::remove_all(nd); np.save(nd);
        const auto r2 = refusal([&] { (void)hc::Owner::restore_as_genesis(hm, o, census, nd); });
        R.record("host-genesis", "a non-vacuum host checkpoint is refused", refused_with(r2, "is not a vacuum state"), {{"refusal", r2}}); }
}

int genesis(const std::string& trace, const fs::path& out, const std::string& nvrtc, const std::string& uuid, const fs::path& v33, const fs::path& v34) {
    Results R;
    Models M; M.trace = trace;
    M.dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    {   auto probe = std::make_shared<const hc::Model>(hc::MeshSpec{45, 4, 1, .009}, trace); M.desc = hc::resident_descriptor_sha256(M.dev, *probe); }
    fs::create_directories(out);
    const auto pre = M.prerun();
    const auto& g = pre->mesh(); const auto o = mf::geometry(g);
    const std::vector<n::ArenaRateEntry> none;
    J info = J::object();

    // ---------------------------------------------------------------- G1: the v3.3 20 ns checkpoint into the pre-run Model
    {
        const auto ck = saved(g, v33);
        std::unique_ptr<hc::ResidentOwner> w;
        const auto r = refusal([&] { w = hc::ResidentOwner::restore_as_genesis(pre, M.dev, o, none, v33, options()); });
        R.record("genesis-v3.3", "the v3.3 20 ns checkpoint is admitted as a genesis in the pre-run Model", r.empty() && w, {{"refusal", r}});
        if (w) {
            axial::State s; hc::ChargeLedger q; w->accepted_state(s, q);
            auto c = word_checks(s, q, ck.s, ck.q, false);
            c["divb_gross_ledger_zero"] = zero_or_empty(q.divb_gross_T_m);
            R.record("genesis-v3.3", "every field and charge word, the admitted Gauss ledger and the clock equal the checkpoint", all_true(c), c);
            hc::ThermalV3 ref(g, pre->drive().start(), {}, 1e-13); ref.load(v33 / "heat");
            const auto h = heat_checks(w->thermal(), ref, 0, ck.s.revision);
            R.record("genesis-v3.3", "slice temperatures, cumulative energies, watermarks and time equal the checkpoint; slice revisions 0", all_true(h), h);
            const auto& p = w->genesis_provenance();
            J pc{{"revision_0", w->revision() == 0}, {"time_kept", bits(w->time_s()) == bits(ck.s.time_s)},
                 {"schema", p.checkpoint_schema == ck.j.at("checkpoint").get<std::string>()},
                 {"configuration", p.configuration_sha256 == ck.j.at("configuration_sha256").get<std::string>()},
                 {"field", p.field_sha256 == ck.j.at("field_sha256").get<std::string>()}, {"heat", p.heat_sha256 == ck.j.at("heat_sha256").get<std::string>()},
                 {"ledgers", p.kinetic_ledgers_sha256 == ck.j.at("kinetic_ledgers_sha256").get<std::string>()},
                 {"motion_mesh", p.motion_mesh_sha256 == ck.j.at("motion_mesh_sha256").get<std::string>()},
                 {"revision", p.revision == ck.s.revision}, {"time", bits(p.time_s) == bits(ck.s.time_s)},
                 {"identity_differs", w->field_sha256() != p.field_sha256 && w->heat_sha256() != p.heat_sha256}};
            R.record("genesis-v3.3", "revision 0, clock kept, provenance equal to the checkpoint's seals", all_true(pc),
                     {{"checks", pc}, {"transfer_text", p.transfer_text}, {"provenance_sha256", p.sha256()}});
            info["v33_genesis"] = {{"field_sha256", w->field_sha256()}, {"heat_sha256", w->heat_sha256()}, {"provenance_sha256", p.sha256()},
                                   {"transfer_text", p.transfer_text}};
            // save / restore of the rebased owner keep the provenance and the digests
            const auto dir = out / "rebased-v33"; fs::remove_all(dir); w->save(dir);
            std::unique_ptr<hc::ResidentOwner> back;
            const auto rr = refusal([&] { back = hc::ResidentOwner::restore(pre, M.dev, o, none, dir, options()); });
            J sc{{"restored", rr.empty() && back}};
            if (back) {
                axial::State s2; hc::ChargeLedger q2; back->accepted_state(s2, q2);
                sc["digests"] = back->field_sha256() == w->field_sha256() && back->heat_sha256() == w->heat_sha256() && back->revision() == 0;
                sc["provenance"] = back->genesis_provenance().sha256() == p.sha256();
                sc["words"] = all_true(word_checks(s2, q2, s, q, true));
            }
            R.record("genesis-v3.3", "save and plain restore of the rebased owner keep its digests and provenance", all_true(sc), {{"checks", sc}, {"refusal", rr}});
        }
    }
    // ---------------------------------------------------------------- G2: into a coordinator-like Model (all transfers)
    {
        const auto ck = saved(g, v33);
        const auto decl = heat_only(g);
        const auto cm = M.make(trace, {45, 4, 1, .009}, coordinator_limits(), decl, 5e-14);
        const auto oa = aligned(g);
        std::unique_ptr<hc::ResidentOwner> w;
        const auto r = refusal([&] { w = hc::ResidentOwner::restore_as_genesis(cm, M.dev, oa, none, v33, options()); });
        J ev{{"refusal", r}, {"motion_mesh_differs", cm->motion_mesh_sha256(oa) != ck.j.at("motion_mesh_sha256").get<std::string>()}};
        R.record("genesis-transfer", "admitted with heat-only patches added, other carrier limits, field-step bound and foil planes", r.empty() && w, ev);
        if (w) {
            axial::State s; hc::ChargeLedger q; w->accepted_state(s, q);
            const auto c = word_checks(s, q, ck.s, ck.q, false);
            R.record("genesis-transfer", "every field and charge word and the admitted Gauss ledger equal the checkpoint", all_true(c), c);
            hc::ThermalV3 ref(g, pre->drive().start(), {}, 1e-13); ref.load(v33 / "heat");
            const auto h = heat_checks(w->thermal(), ref, 0, ck.s.revision);
            R.record("genesis-transfer", "slice states equal the checkpoint; the added heat-only patches start with zero inventory", all_true(h), h);
            const auto& t = w->genesis_provenance().transfer_text;
            J tc{{"declaration", t.find("thermal_declaration=none->" + decl.sha256()) != std::string::npos},
                 {"limits", t.find("carrier_limits=100000/1073741824->2000000/") != std::string::npos},
                 {"mesh", t.find("motion_mesh=" + ck.j.at("motion_mesh_sha256").get<std::string>() + "->" + cm->motion_mesh_sha256(oa)) != std::string::npos},
                 {"dt", t.find("dt_max_s=1e-13->5.0000000000000002e-14") != std::string::npos},
                 {"rule_form", t.find("rule_form=v3.3") != std::string::npos}};
            R.record("genesis-transfer", "the transfer text records every operand that differs", all_true(tc), {{"checks", tc}, {"transfer_text", t}});
            // a few vacuum steps at the new bound
            const auto rs = refusal([&] { for (int k = 0; k < 20; ++k) vacuum_step(*w, 5e-14); });
            R.record("genesis-transfer", "20 vacuum steps after the genesis (revision 20, clock advanced)", rs.empty() && w->revision() == 20,
                     {{"refusal", rs}, {"revision", w->revision()}});
        }
        // transfers not admitted by GenesisTransfer
        hc::GenesisTransfer no_decl; no_decl.thermal_declaration = false;
        const auto r1 = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(M.make(trace, {45, 4, 1, .009}, hc::CarrierLimits{}, decl, 1e-13), M.dev, o, none, v33, options(), no_decl); });
        R.record("genesis-refusal", "added patches refused when the declaration transfer is not admitted", refused_with(r1, "may only add heat-only source patches"), {{"refusal", r1}});
        hc::GenesisTransfer no_lim; no_lim.carrier_limits = false;
        const auto r2 = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(M.make(trace, {45, 4, 1, .009}, coordinator_limits(), {}, 1e-13), M.dev, o, none, v33, options(), no_lim); });
        R.record("genesis-refusal", "other carrier limits refused when the limits transfer is not admitted", refused_with(r2, "genesis carrier limits differ"), {{"refusal", r2}});
        hc::GenesisTransfer no_geo; no_geo.native_geometry = false;
        const auto r3 = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(pre, M.dev, oa, none, v33, options(), no_geo); });
        R.record("genesis-refusal", "other foil planes refused when the geometry transfer is not admitted", refused_with(r3, "genesis native geometry differs"), {{"refusal", r3}});
    }
    // ---------------------------------------------------------------- G3: the v3.4 checkpoint (divB gross ledger kept)
    std::unique_ptr<hc::ResidentOwner> g34, p34, g33;
    {
        const auto ck = saved(g, v34);
        const auto r = refusal([&] { g34 = hc::ResidentOwner::restore_as_genesis(pre, M.dev, o, none, v34, options()); });
        R.record("genesis-v3.4", "the v3.4 20 ns checkpoint is admitted as a genesis", r.empty() && g34, {{"refusal", r}});
        if (g34) {
            axial::State s; hc::ChargeLedger q; g34->accepted_state(s, q);
            const auto c = word_checks(s, q, ck.s, ck.q, true);
            R.record("genesis-v3.4", "every word, the admitted Gauss ledger and the divB gross ledger equal the checkpoint", all_true(c), c);
            const auto& pv = g34->genesis_provenance();
            R.record("genesis-v3.4", "provenance: v3.4 schema, the checkpoint's seals, revision and clock",
                     pv.checkpoint_schema == "hall-v3.4-resident-checkpoint-v1" && pv.field_sha256 == ck.j.at("field_sha256").get<std::string>() &&
                     pv.revision == ck.s.revision && g34->revision() == 0 && pv.transfer_text.find("rule_form=v3.4") != std::string::npos,
                     {{"transfer_text", pv.transfer_text}});
        }
        const auto rp = refusal([&] { p34 = hc::ResidentOwner::restore(pre, M.dev, o, none, v34, options()); });
        R.record("continuation", "plain restore of the v3.4 checkpoint (revision kept)", rp.empty() && p34 && p34->revision() == ck.s.revision, {{"refusal", rp}});
        const auto r3 = refusal([&] { g33 = hc::ResidentOwner::restore_as_genesis(pre, M.dev, o, none, v33, options()); });
        R.record("continuation", "genesis of the v3.3 checkpoint for the lockstep", r3.empty() && g33, {{"refusal", r3}});
    }
    // ---------------------------------------------------------------- continuation: rebased vs plain restore, 200 vacuum steps
    if (g34 && p34 && g33) {
        bool words = true, words33 = true, heat = true, revs = true; std::string why;
        long steps = 0;
        const auto r = refusal([&] {
            for (int k = 0; k < 200; ++k) {
                vacuum_step(*g34, 1e-13); vacuum_step(*p34, 1e-13); vacuum_step(*g33, 1e-13); ++steps;
                axial::State a, b, c; hc::ChargeLedger qa, qb, qc;
                g34->accepted_state(a, qa); p34->accepted_state(b, qb); g33->accepted_state(c, qc);
                const auto wc = word_checks(a, qa, b, qb, true); if (!all_true(wc) && words) { words = false; why = wc.dump(); }
                auto wc33 = word_checks(c, qc, b, qb, false); if (!all_true(wc33) && words33) { words33 = false; why += wc33.dump(); }
                const auto h = heat_checks(g34->thermal(), p34->thermal(), std::uint64_t(k + 1), p34->revision());
                J hh = h; hh.erase("revisions"); if (!all_true(hh)) heat = false;
                const auto h33 = heat_checks(g33->thermal(), p34->thermal(), std::uint64_t(k + 1), p34->revision());
                J hh33 = h33; hh33.erase("revisions"); if (!all_true(hh33)) heat = false;
                revs = revs && g34->revision() == std::uint64_t(k + 1) && g33->revision() == std::uint64_t(k + 1);
            }
        });
        R.record("continuation", "200 vacuum steps: the rebased v3.4 owner equals the plain restore word for word (field, charge, Gauss and divB ledgers)",
                 r.empty() && words && steps == 200, {{"refusal", r}, {"steps", steps}, {"first_difference", why}});
        R.record("continuation", "200 vacuum steps: the rebased v3.3 owner equals it in every field and charge word (divB ledger from zero)",
                 r.empty() && words33 && steps == 200, {{"refusal", r}});
        R.record("continuation", "200 vacuum steps: slice thermal states equal (revisions 1..200 against N+1..N+200)", r.empty() && heat && revs, {{"refusal", r}});
        info["continuation_end"] = {{"rebased_v34_field_sha256", g34->field_sha256()}, {"plain_field_sha256", p34->field_sha256()},
                                    {"plain_revision", p34->revision()}, {"rebased_revision", g34->revision()}};
    }
    g34.reset(); p34.reset(); g33.reset();
    // ---------------------------------------------------------------- fresh genesis: vacuum and coupled lockstep
    {
        const auto seed = pre->dc_seed();
        hc::ResidentOwner fresh(pre, M.dev, seed, declared(*pre, seed, hc::charge_of_census(g, o, none)), o, none, options());
        const auto dir = out / "fresh-genesis-r0"; fs::remove_all(dir); fresh.save(dir);
        std::unique_ptr<hc::ResidentOwner> w;
        const auto r = refusal([&] { w = hc::ResidentOwner::restore_as_genesis(pre, M.dev, o, none, dir, options()); });
        R.record("fresh", "a revision-0 checkpoint of a fresh genesis is admitted as a genesis", r.empty() && w, {{"refusal", r}});
        if (w) {
            bool ok = true; std::string why;
            const auto rs = refusal([&] {
                for (int k = 0; k < 100; ++k) {
                    vacuum_step(fresh, 1e-13); vacuum_step(*w, 1e-13);
                    axial::State a, b; hc::ChargeLedger qa, qb; fresh.accepted_state(a, qa); w->accepted_state(b, qb);
                    const auto wc = word_checks(a, qa, b, qb, true); const auto h = heat_checks(fresh.thermal(), w->thermal(), fresh.revision(), w->revision());
                    if ((!all_true(wc) || !all_true(h) || fresh.revision() != w->revision()) && ok) { ok = false; why = wc.dump() + h.dump(); }
                }
            });
            R.record("fresh", "100 vacuum steps: the rebased owner equals the fresh genesis in every word and slice state (same revisions)", rs.empty() && ok,
                     {{"refusal", rs}, {"first_difference", why}});
            R.record("fresh", "the rebased owner's digests differ from the fresh genesis (provenance bound)", w->field_sha256() != fresh.field_sha256() &&
                     w->heat_sha256() != fresh.heat_sha256(), {});
        }
        // coupled: fresh genesis with co-located neutral pairs against a genesis from a vacuum revision-0 checkpoint with
        // the same census (zero net per node)
        auto census = mf::pairs(400, 11, 1e5);
        hc::ResidentOwner fresh_pairs(pre, M.dev, seed, declared(*pre, seed, hc::charge_of_census(g, o, census)), o, census, options());
        std::unique_ptr<hc::ResidentOwner> wp;
        const auto rp = refusal([&] { wp = hc::ResidentOwner::restore_as_genesis(pre, M.dev, o, census, dir, options()); });
        R.record("fresh", "a genesis from a vacuum checkpoint with a census of co-located neutral pairs is admitted", rp.empty() && wp, {{"refusal", rp}});
        if (wp) {
            gp::GpuMotionProvider provider(M.dev, fresh_pairs.motion_context(), std::max<std::size_t>(pre->limits().maximum_carriers, 1));
            bool decisions = true, motion = true, words = true, heat = true; std::string why; long accepted = 0;
            std::vector<n::ArenaRateEntry> ca = census, cb = census;
            const auto rs = refusal([&] {
                for (long k = 0; k < 6; ++k) {
                    const double dt = 1e-14;
                    const std::string src = hc::sha256("v35/source");
                    const auto ta = tables_of(ca), tb = tables_of(cb);
                    hc::StepInputs ia{{}, {}, false, ta.states, ta.keys, ta.idx}, ib{{}, {}, false, tb.states, tb.keys, tb.idx};
                    std::optional<hc::ResidentJointProposal> pa, pb;
                    const auto ea = refusal([&] { pa.emplace(hc::prepare_closed_proposal(fresh_pairs, o, ca, binding_of(fresh_pairs, *pre, k, dt), src, provider, ia)); });
                    const auto eb = refusal([&] { pb.emplace(hc::prepare_closed_proposal(*wp, o, cb, binding_of(*wp, *pre, k, dt), src, provider, ib)); });
                    if (ea != eb || bool(pa) != bool(pb)) { decisions = false; why = ea + " | " + eb; return; }
                    if (!pa) return;
                    const auto &mx = pa->motion(), &my = pb->motion();
                    motion = motion && mx.record_v7_sha256() == my.record_v7_sha256() && mx.census_sha256() == my.census_sha256() &&
                             mx.binding().trial_field_sha256 == my.binding().trial_field_sha256 && pa->iterations() == pb->iterations() &&
                             bits(pa->combined_balance_J()) == bits(pb->combined_balance_J());
                    Boundary b1{pa->candidate_field_sha256(), pa->candidate_heat_sha256(), src, mx.binding().mesh_sha256, mx.binding().material_sha256,
                                mx.binding().common_clock_sha256, pa->end_s(), pa->endpoint_revision()};
                    Boundary b2{pb->candidate_field_sha256(), pb->candidate_heat_sha256(), src, my.binding().mesh_sha256, my.binding().material_sha256,
                                my.binding().common_clock_sha256, pb->end_s(), pb->endpoint_revision()};
                    const bool fa = fresh_pairs.preflight_final(*pa, ca, b1), fb = wp->preflight_final(*pb, cb, b2);
                    if (fa != fb || !fa) { decisions = false; why = "preflight " + std::to_string(fa) + std::to_string(fb); return; }
                    auto na = hc::next_census(ca, mx); auto nb = hc::next_census(cb, my);
                    fresh_pairs.publish_prevalidated(std::move(*pa)); wp->publish_prevalidated(std::move(*pb));
                    ca = std::move(na); cb = std::move(nb); ++accepted;
                    axial::State a, b; hc::ChargeLedger qa, qb; fresh_pairs.accepted_state(a, qa); wp->accepted_state(b, qb);
                    const auto wc = word_checks(a, qa, b, qb, true); if (!all_true(wc) && words) { words = false; why = wc.dump(); }
                    const auto h = heat_checks(fresh_pairs.thermal(), wp->thermal(), fresh_pairs.revision(), wp->revision());
                    if (!all_true(h)) heat = false;
                }
            });
            R.record("fresh", "coupled steps: same decisions, iterations, census v2, record v7, trial field and balance as the fresh genesis",
                     rs.empty() && decisions && motion && accepted == 6, {{"refusal", rs}, {"accepted", accepted}, {"why", why}});
            R.record("fresh", "coupled steps: every accepted field, charge and ledger word and every slice state equal the fresh genesis",
                     rs.empty() && words && heat && accepted == 6, {{"first_difference", why}});
        }
    }
    // ---------------------------------------------------------------- refusals
    {
        // flipped field word (one bit of E[100]) and flipped heat word (one bit of a slice temperature)
        const auto tw = out / "tamper-field-word"; copy_dir(v34, tw);
        {   std::ifstream in(tw / "checkpoint.json"); J j = J::parse(in); in.close();
            auto& e = j["field"]["E"][100]; e = std::bit_cast<double>(bits(e.get<double>()) ^ 1);
            std::ofstream(tw / "checkpoint.json", std::ios::trunc) << j.dump() << '\n'; }
        const auto r1 = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(pre, M.dev, o, none, tw, options()); });
        R.record("genesis-refusal", "a flipped field word is refused by the seal", refused_with(r1, "checkpoint field/charge seal mismatch"), {{"refusal", r1}});
        const auto th = out / "tamper-heat-word"; copy_dir(v34, th);
        {   std::fstream f(th / "heat" / "0.bin", std::ios::in | std::ios::out | std::ios::binary);
            std::uint64_t n = 0; f.seekg(8); f.read(reinterpret_cast<char*>(&n), 8);   // magic, signature length
            const std::streamoff t0 = 16 + std::streamoff(n) + 8 + 8 + 8;            // signature, time, revision, count -> T[0]
            double x = 0; f.seekg(t0); f.read(reinterpret_cast<char*>(&x), 8); x = std::bit_cast<double>(bits(x) ^ 1);
            f.seekp(t0); f.write(reinterpret_cast<const char*>(&x), 8); }
        const auto r2 = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(pre, M.dev, o, none, th, options()); });
        R.record("genesis-refusal", "a flipped slice temperature is refused by the thermal seal", refused_with(r2, "checkpoint thermal seal mismatch"), {{"refusal", r2}});
        // non-vacuum checkpoint (free charge of co-located pairs) under the pre-run Model
        const auto seed = pre->dc_seed();
        auto census = mf::pairs(40, 5, 1e5);
        {   hc::ResidentOwner np(pre, M.dev, seed, declared(*pre, seed, hc::charge_of_census(g, o, census)), o, census, options());
            const auto dir = out / "non-vacuum"; fs::remove_all(dir); np.save(dir);
            const auto r3 = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(pre, M.dev, o, census, dir, options()); });
            R.record("genesis-refusal", "a non-vacuum checkpoint is refused", refused_with(r3, "is not a vacuum state"), {{"refusal", r3}}); }
        // a checkpoint written under a declaration; restoring declarations that alter or remove its patch
        const auto d1 = heat_only(g);
        const auto md = M.make(trace, {45, 4, 1, .009}, hc::CarrierLimits{}, d1, 1e-13);
        const auto dir = out / "declared-r0";
        {   hc::ResidentOwner dw(md, M.dev, seed, declared(*md, seed, hc::charge_of_census(g, o, none)), o, none, options()); fs::remove_all(dir); dw.save(dir); }
        hc::GenesisTransfer with_rec; with_rec.recorded_declaration = d1;
        const auto altered = M.make(trace, {45, 4, 1, .009}, hc::CarrierLimits{}, heat_only(g, .5), 1e-13);
        const auto r4 = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(altered, M.dev, o, none, dir, options(), with_rec); });
        R.record("genesis-refusal", "a declaration that alters a recorded patch is refused", refused_with(r4, "may only add heat-only source patches"), {{"refusal", r4}});
        const auto r5 = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(pre, M.dev, o, none, dir, options(), with_rec); });
        R.record("genesis-refusal", "a declaration that removes a recorded patch is refused", refused_with(r5, "may only add heat-only source patches"), {{"refusal", r5}});
        const auto r5b = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(altered, M.dev, o, none, dir, options()); });
        R.record("genesis-refusal", "without the recorded declaration the altered Model is not a transfer of the checkpoint",
                 refused_with(r5b, "is not a transfer of the restoring Model"), {{"refusal", r5b}});
        std::unique_ptr<hc::ResidentOwner> same;
        const auto r5c = refusal([&] { same = hc::ResidentOwner::restore_as_genesis(md, M.dev, o, none, dir, options()); });
        R.record("genesis-refusal", "control: the same declaration (v3.5 checkpoint) is admitted", r5c.empty() && same, {{"refusal", r5c}});
        same.reset();
        // other drive trace (one value of the C19 trace changed) and other mesh (axial depth)
        const auto tr = out / "trace-altered.csv";
        {   std::ifstream in(trace); std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            const auto nl = text.find('\n', text.find('\n', text.find('\n') + 1) + 1);   // third data line
            auto comma = text.rfind(',', nl); std::string tail = text.substr(comma + 1, nl - comma - 1);
            const double v = std::stod(tail); char b[64]; std::snprintf(b, sizeof b, "%.17g", v * (1 + 1e-6) + 1e-3);
            text = text.substr(0, comma + 1) + b + text.substr(nl);
            std::ofstream(tr) << text; }
        const auto r6 = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(M.make(tr.string(), {45, 4, 1, .009}, hc::CarrierLimits{}, {}, 1e-13), M.dev, o, none, v33, options()); });
        // the drive is pinned to the C19 bytes: a Model on another trace is refused at its construction (the configuration
        // binds the C19 digest as well)
        R.record("genesis-refusal", "another drive trace is refused", refused_with(r6, "C19 bytes differ") || refused_with(r6, "is not a transfer of the restoring Model"),
                 {{"refusal", r6}});
        // another mesh: the transverse resolution of the 100 um runs (other dimensions: refused at the checkpoint's word
        // dimensions) and the 200 um mesh with one more air layer (same transverse dimensions, other mesh record)
        const auto r7 = refusal([&] { (void)hc::ResidentOwner::restore_as_genesis(M.make(trace, {90, 4, 1, .009}, hc::CarrierLimits{}, {}, 1e-13), M.dev,
                                                                                   mf::geometry(M.make(trace, {90, 4, 1, .009}, hc::CarrierLimits{}, {}, 1e-13)->mesh()), none, v33, options()); });
        R.record("genesis-refusal", "another mesh (100 um) is refused", refused_with(r7, "dimensions"), {{"refusal", r7}});
        const auto r8 = refusal([&] { const auto m2 = M.make(trace, {45, 4, 2, .009}, hc::CarrierLimits{}, {}, 1e-13);
                                      (void)hc::ResidentOwner::restore_as_genesis(m2, M.dev, mf::geometry(m2->mesh()), none, v33, options()); });
        R.record("genesis-refusal", "another mesh (one more air layer) is refused", refused_with(r8, "dimensions") || refused_with(r8, "is not a transfer of the restoring Model"),
                 {{"refusal", r8}});
    }
    host_genesis(R, trace, out);
    J res{{"schema", "hall-v3.5-genesis-v1"}, {"cases", R.cases}, {"failed", R.failed}, {"info", info},
          {"v33_checkpoint", v33.string()}, {"v34_checkpoint", v34.string()}, {"device", J::parse(M.dev->device_json())}};
    std::ofstream(out / "genesis.json") << res.dump(1) << '\n';
    std::cout << "failed " << R.failed << std::endl;
    return R.failed ? 1 : 0;
}

// Fresh-process continuation of a rebased owner: genesis-run (genesis from CHECKPOINT, 50 vacuum steps, save, 50 more steps,
// reference digests) and genesis-continue (plain restore of the saved owner in a new process, 50 steps, compared with the
// reference).
struct Device1 { std::shared_ptr<gp::Device> dev; std::shared_ptr<const hc::Model> model; };
Device1 device1(const std::string& trace, const std::string& nvrtc, const std::string& uuid) {
    Device1 d; d.dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    Models M; M.trace = trace; M.dev = d.dev;
    {   auto probe = std::make_shared<const hc::Model>(hc::MeshSpec{45, 4, 1, .009}, trace); M.desc = hc::resident_descriptor_sha256(d.dev, *probe); }
    d.model = M.prerun(); return d;
}
J digests_of(const hc::ResidentOwner& r) {
    axial::State s; hc::ChargeLedger q; r.accepted_state(s, q);
    hc::Digest d("v35/physical-words"); d.f64s(s.E); d.f64s(s.B); d.f64s(s.E_low); d.f64s(s.B_low); d.f64(s.time_s);
    for (const auto& st : r.thermal().states()) { d.f64s(st.T); for (double c : st.cumulative) d.f64(c); d.f64(st.time); }
    return {{"field_sha256", r.field_sha256()}, {"heat_sha256", r.heat_sha256()}, {"kinetic_ledgers_sha256", r.kinetic_ledgers_sha256()},
            {"revision", r.revision()}, {"time_s", r.time_s()}, {"physical_sha256", d.hex()}, {"provenance_sha256", r.genesis_provenance().sha256()}};
}
int genesis_run(const std::string& trace, const fs::path& out, const std::string& nvrtc, const std::string& uuid, const fs::path& ck) {
    auto D = device1(trace, nvrtc, uuid); const auto o = mf::geometry(D.model->mesh()); const std::vector<n::ArenaRateEntry> none;
    auto w = hc::ResidentOwner::restore_as_genesis(D.model, D.dev, o, none, ck, options());
    for (int k = 0; k < 50; ++k) vacuum_step(*w, 1e-13);
    fs::create_directories(out); const auto dir = out / "rebased-50"; fs::remove_all(dir); w->save(dir);
    J at_save = digests_of(*w);
    for (int k = 0; k < 50; ++k) vacuum_step(*w, 1e-13);
    J ref{{"schema", "hall-v3.5-genesis-continuation-reference-v1"}, {"checkpoint", dir.string()}, {"at_save", at_save}, {"after_50", digests_of(*w)}};
    std::ofstream(out / "continuation-reference.json") << ref.dump(1) << '\n';
    std::cout << ref.dump(1) << std::endl;
    return 0;
}
int genesis_continue(const std::string& trace, const fs::path& out, const std::string& nvrtc, const std::string& uuid, const fs::path& dir, const fs::path& refp) {
    auto D = device1(trace, nvrtc, uuid); const auto o = mf::geometry(D.model->mesh()); const std::vector<n::ArenaRateEntry> none;
    std::ifstream in(refp); const J ref = J::parse(in);
    Results R;
    std::unique_ptr<hc::ResidentOwner> w;
    const auto r = refusal([&] { w = hc::ResidentOwner::restore(D.model, D.dev, o, none, dir, options()); });
    R.record("continuation", "fresh-process plain restore of the rebased owner's checkpoint", r.empty() && w && digests_of(*w) == ref.at("at_save"),
             {{"refusal", r}, {"restored", w ? digests_of(*w) : J()}, {"reference", ref.at("at_save")}});
    if (w) {
        for (int k = 0; k < 50; ++k) vacuum_step(*w, 1e-13);
        R.record("continuation", "50 further steps in the new process equal the uninterrupted run (digests, provenance, physical words, slice states)",
                 digests_of(*w) == ref.at("after_50"), {{"continued", digests_of(*w)}, {"reference", ref.at("after_50")}});
    }
    fs::create_directories(out);
    std::ofstream(out / "genesis-continue.json") << J{{"schema", "hall-v3.5-genesis-continue-v1"}, {"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    return R.failed ? 1 : 0;
}
} // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "genesis" && argc == 8) return genesis(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
        if (mode == "genesis-run" && argc == 7) return genesis_run(argv[2], argv[3], argv[4], argv[5], argv[6]);
        if (mode == "genesis-continue" && argc == 8) return genesis_continue(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
        std::cerr << "usage: test_v35 genesis TRACE OUT NVRTC UUID V33_CHECKPOINT V34_CHECKPOINT\n"
                     "       test_v35 genesis-run TRACE OUT NVRTC UUID CHECKPOINT\n"
                     "       test_v35 genesis-continue TRACE OUT NVRTC UUID REBASED_CHECKPOINT REFERENCE_JSON\n";
        return 2;
    } catch (const std::exception& e) { std::cerr << "error: " << e.what() << std::endl; return 1; }
}
