// Hall corrective v3.4 qualification: the device coupled step against the v3.3 host path.
//   test_v34 parity TRACE OUT NVRTC UUID [TRANSVERSE]
//       two resident owners on one device and one Model, host-work modes v33 (GpuMotionProvider, host census,
//       host support, host fold, host record, host residual) and device (DeviceStep), stepped in lockstep over the
//       fixtures; per step: decision and refusal text, candidate digests (same configuration -> equal Merkle
//       digests iff equal state), record v7 and census v2 digests, budget, work scales, sinks and births, the
//       record arrays, the independent residual values, the deposit support (device array vs the sealed
//       DepositSupportBuilder on the v33 receipt), the bank transfer (sparse vs the sealed bank_transfer), the
//       accepted states after publication.
//   test_v34 tamper TRACE OUT NVRTC UUID
//       device preflight negative controls (record words, census rows, census span, stale proposal, candidate slot).
//   test_v34 bench TRACE OUT NVRTC UUID TRANSVERSE CARRIERS STEPS MODE     MODE: owner | v33 | device
//       load fixture (49 % electrons, 49 % co-located H2O+, 2 % H2O of weight 1e6, in the gap), dt 10 fs; host
//       (thread CPU) and wall times per item; one owner per process.
//   test_v34 blocks TRACE OUT NVRTC UUID
//       block-form digests of binary64 arrays (native v7.1 tags) against the host reference (forms core).
//   test_v34 hostloop TRACE OUT NVRTC UUID CARRIERS REPS
//       the host DepositSupportBuilder (v3.3) and two reference loops (arithmetic only; memory streaming) timed in one
//       process before any CUDA context, after the device context, after the resident owner and push engine
//       buffers, and right after device steps (the 4-6x resident-process slowdown of HANDOVER-V3.4 item 2).
//   test_v34 divb TRACE OUT NVRTC UUID STEPS
//       divb-coefficient-gross-v1: host Owner (cpu-pcg) and resident owner (host PCG parity solve) in lockstep, the divB
//       gross ledger G_c compared bitwise every step; one-ulp G_c tamper (device slot, unsealed host and resident
//       checkpoints), nonzero bootstrap ledger; the seed's coefficient-form divB.
#include "resident_owner.hpp"
#include "gpu/gpu_backends.hpp"
#include "gpu/device_step_v34.hpp"
#include "common/hall_gpu_common.hpp"
#include "hall-finite-foil-em-20261029/axial_particles.hpp"
namespace hc = hall_corrective_v3;
#include "motion_fixtures.hpp"
#include "sink_fixtures.hpp"
#include <ctime>
#include <fstream>
#include <iostream>
#include <set>

namespace n = raikou::mfem_backend;
namespace mf = motion_fixtures;
namespace sf = sink_fixtures;
namespace gp = hall_corrective_v3::gpu;
using J = hc::Json;
using steady = std::chrono::steady_clock;

namespace {
double ms_since(steady::time_point t) { return std::chrono::duration<double, std::milli>(steady::now() - t).count(); }
double cpu_ms() { timespec ts{}; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); return ts.tv_sec * 1e3 + ts.tv_nsec * 1e-6; }
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
hc::ChargeLedger declared(const hc::Model& m, const axial::State& s, hc::FreeCharge particles) {
    const auto& g = m.mesh(); hc::ChargeLedger q; q.particles = std::move(particles);
    q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
    q.conductor_declaration = hc::sha256("parity: declared DC conductor/PEC induced charge");
    q.exterior_declaration = hc::sha256("parity: declared DC exterior termination charge");
    const auto field = hc::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = field[k]; else if (g.outer_node[k]) q.exterior_C[k] = field[k]; }
    return q;
}
template <class OwnerT>
hall_driver::MotionBinding binding_of(const OwnerT& res, const hc::Model& m, long step, double dt) {
    hall_driver::MotionBinding b; b.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = res.revision(); b.staged_arena_revision = res.revision() + 1; b.begin_s = res.time_s(); b.dt_s = dt; b.end_s = b.begin_s + dt;
    b.staged_arena_sha256 = hc::sha256("v34/staged/" + std::to_string(step)); b.mesh_sha256 = res.motion_mesh_sha256();
    b.material_sha256 = m.material_sha256(); b.species_sha256 = hc::sha256("v34/species");
    b.common_clock_sha256 = hc::sha256("v34/clock/" + std::to_string(step)); b.kinetic_policy_sha256 = res.motion_context().policy_sha256;
    return b;
}
struct Setup { std::shared_ptr<gp::Device> dev; std::shared_ptr<const hc::Model> model, owner_model; std::string descriptor; };
hc::CarrierLimits limits_for(std::size_t cap) { hc::CarrierLimits l; l.maximum_carriers = std::max<std::size_t>(cap, 1000); l.memory_budget_bytes = l.maximum_carriers * 4096 + (512ull << 20); return l; }
Setup setup(const std::string& trace, const std::string& nvrtc, const std::string& uuid, int transverse, std::size_t cap, bool owner_model = false) {
    Setup s;
    s.dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    const hc::MeshSpec spec{transverse, 4, 1, .009};
    const auto limits = limits_for(cap);
    if (owner_model) {   // the coordinator's Owner path: host arithmetic, gpu1 field solver, native v7 forms
        s.owner_model = std::make_shared<const hc::Model>(spec, trace, gp::make_gpu_field_solver(s.dev), limits, hc::ThermalDeclaration{}, 1e-14, hc::MotionForms{true, 1e-14});
        return s;
    }
    auto probe = std::make_shared<const hc::Model>(spec, trace);
    s.descriptor = hc::resident_descriptor_sha256(s.dev, *probe);
    s.model = std::make_shared<const hc::Model>(spec, trace, std::make_shared<hc::ResidentBackend>(s.descriptor), limits, hc::ThermalDeclaration{}, 1e-14,
                                                hc::MotionForms{true, 1e-14});
    return s;
}
gp::ResidentOptions options(gp::HostWork w) { gp::ResidentOptions o; o.host_work = w; o.verify_samples = 16; o.sample_seed = 7; return o; }
// native v7.1 device lane configuration: device path, no host record arrays, block-form field/current digests
gp::ResidentOptions lane_options() {
    auto o = options(gp::HostWork::device); o.host_record_arrays = false; o.field_digest_form = gp::FieldDigestForm::blocks_v71; return o;
}
// census v2 tables of the coordinator (native trial.census_v7()): packet key = smallest ID of the packet, state index.
struct Tables { std::vector<std::string> states; std::vector<std::uint64_t> keys, idx; };
Tables tables_of(const std::vector<n::ArenaRateEntry>& c) {
    Tables t; std::set<std::string> u; for (const auto& r : c) u.insert(r.state); t.states.assign(u.begin(), u.end());
    std::map<std::string, std::uint64_t> key; for (const auto& r : c) { auto [it, ins] = key.emplace(r.packet_id, r.native.id); if (!ins) it->second = std::min(it->second, r.native.id); }
    std::vector<const n::ArenaRateEntry*> rows; for (const auto& r : c) rows.push_back(&r);   // tables in ascending native ID (native census_v7)
    std::stable_sort(rows.begin(), rows.end(), [](auto a, auto b) { return a->native.id < b->native.id; });
    for (const auto* r : rows) { t.keys.push_back(key.at(r->packet_id)); t.idx.push_back(std::uint64_t(std::lower_bound(t.states.begin(), t.states.end(), r->state) - t.states.begin())); }
    return t;
}
J scales_json(const hc::WorkScales& w) {
    return {{"G_K", w.kinetic_gross_J}, {"C_K", w.kinetic_content_J}, {"G_W", w.field_work_gross_J}, {"Q", w.charge_abs_C}, {"E_hat", w.field_bound_V_m},
            {"tolerance", w.tolerance_J}, {"status", w.status}};
}
bool same_scales(const hc::WorkScales& a, const hc::WorkScales& b) {
    return bits(a.kinetic_gross_J) == bits(b.kinetic_gross_J) && bits(a.kinetic_content_J) == bits(b.kinetic_content_J) &&
           bits(a.field_work_gross_J) == bits(b.field_work_gross_J) && bits(a.charge_abs_C) == bits(b.charge_abs_C) && bits(a.field_bound_V_m) == bits(b.field_bound_V_m) &&
           bits(a.tolerance_J) == bits(b.tolerance_J) && bits(a.v31_tolerance_J) == bits(b.v31_tolerance_J) && a.status == b.status;
}
bool same_budget(const hall_driver::MotionBudget& a, const hall_driver::MotionBudget& b) {
    const double x[7] = {a.kinetic_change_J, a.field_work_J, a.work_defect_J, a.continuity_C, a.plasma_dt, a.maximum_gyro, a.maximum_displacement_cells};
    const double y[7] = {b.kinetic_change_J, b.field_work_J, b.work_defect_J, b.continuity_C, b.plasma_dt, b.maximum_gyro, b.maximum_displacement_cells};
    for (int k = 0; k < 7; ++k) if (bits(x[k]) != bits(y[k])) return false;
    return true;
}
bool same_sinks(std::span<const hc::SinkRecord> a, std::span<const hc::SinkRecord> b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const auto &x = a[i], &y = b[i];
        if (x.native_id != y.native_id || x.kind != y.kind || x.boundary != y.boundary || x.axis != y.axis || x.packet_id != y.packet_id || x.state != y.state ||
            x.collected != y.collected || x.conductor != y.conductor) return false;
        for (double p : {x.plane_m - y.plane_m, x.fraction - y.fraction}) (void)p;
        const double xs[] = {x.plane_m, x.fraction, x.flight_s, x.discarded_flight_s, x.charge_C, x.mass_kg, x.weight, x.kinetic_J, x.weighted_kinetic_J,
                             x.initial_weighted_kinetic_J, x.gyro, x.displacement_cells};
        const double ys[] = {y.plane_m, y.fraction, y.flight_s, y.discarded_flight_s, y.charge_C, y.mass_kg, y.weight, y.kinetic_J, y.weighted_kinetic_J,
                             y.initial_weighted_kinetic_J, y.gyro, y.displacement_cells};
        for (std::size_t k = 0; k < std::size(xs); ++k) if (bits(xs[k]) != bits(ys[k])) return false;
        for (int d = 0; d < 3; ++d)
            if (bits(x.segment_m[d]) != bits(y.segment_m[d]) || bits(x.hit_point_m[d]) != bits(y.hit_point_m[d]) ||
                bits(x.incident_proper_velocity_m_s[d]) != bits(y.incident_proper_velocity_m_s[d]) || bits(x.field_point_m[d]) != bits(y.field_point_m[d]))
                return false;
    }
    return true;
}
bool same_births(std::span<const hc::NewbornRecord> a, std::span<const hc::NewbornRecord> b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const auto &x = a[i], &y = b[i];
        if (x.native_id != y.native_id || bits(x.birth_time_s) != bits(y.birth_time_s) || bits(x.flight_s) != bits(y.flight_s) || x.origin != y.origin ||
            x.axis != y.axis || x.packet_id != y.packet_id || x.state != y.state || bits(x.charge_C) != bits(y.charge_C) || bits(x.weight) != bits(y.weight) ||
            bits(x.birth_weighted_kinetic_J) != bits(y.birth_weighted_kinetic_J)) return false;
        for (int d = 0; d < 3; ++d) if (bits(x.start_m[d]) != bits(y.start_m[d])) return false;
    }
    return true;
}
template <class P> J refusal_of(P&& f) {
    try { f(); return J{{"accepted", true}}; } catch (const std::exception& e) { return J{{"accepted", false}, {"refusal", e.what()}}; }
}

// Host reference of a binary64 block list with the native v7.1 tags (forms core primitives).
std::string host_root(const char* tag, const std::vector<double>& v) {
    const std::uint64_t n = v.size(), chunks = (n + 31) / 32, leaves = (n + 1023) / 1024;
    std::vector<std::array<pf7_u32, 8>> ch(chunks), lf(leaves);
    const std::string t(tag);
    for (std::uint64_t j = 0; j < chunks; ++j) {
        const std::uint64_t m = std::min<std::uint64_t>(32, n - 32 * j);
        pf7::Sha s; pf7::sha_init(s); const std::string ct = t + "/chunk"; pf7::sha_text(s, ct.data(), pf7_u32(ct.size()));
        pf7::sha_u64(s, j); pf7::sha_u64(s, m); for (std::uint64_t r = 0; r < m; ++r) pf7::sha_u64(s, std::bit_cast<std::uint64_t>(v[32 * j + r]));
        pf7::sha_finish(s, ch[j].data());
    }
    for (std::uint64_t k = 0; k < leaves; ++k) {
        const std::uint64_t rec = std::min<std::uint64_t>(1024, n - 1024 * k), q = (rec + 31) / 32;
        pf7::Sha s; pf7::sha_init(s); const std::string bt = t + "/block"; pf7::sha_text(s, bt.data(), pf7_u32(bt.size()));
        pf7::sha_u64(s, k); pf7::sha_u64(s, rec); pf7::sha_u64(s, q); for (std::uint64_t c = 0; c < q; ++c) pf7::sha_digest_words(s, ch[32 * k + c].data());
        pf7::sha_finish(s, lf[k].data());
    }
    std::array<pf7_u32, 8> tree{};
    if (lf.empty()) pf7::merkle_empty(tree.data());
    else {
        auto level = lf;
        while (level.size() > 1) { decltype(level) next; for (std::size_t p = 0; p + 1 < level.size(); p += 2) { std::array<pf7_u32, 8> x; pf7::merkle_node(level[p].data(), level[p + 1].data(), x.data()); next.push_back(x); }
            if (level.size() & 1) next.push_back(level.back()); level = std::move(next); }
        tree = level[0];
    }
    pf7::Sha s; pf7::sha_init(s); const std::string rt = t + "/root"; pf7::sha_text(s, rt.data(), pf7_u32(rt.size()));
    pf7::sha_u64(s, 1024); pf7::sha_u64(s, n); pf7::sha_u64(s, leaves); pf7::sha_digest_words(s, tree.data());
    std::array<pf7_u32, 8> root; pf7::sha_finish(s, root.data());
    return pf_v7::hex(pf_v7::digest_of_words(root.data()));
}

// ---------------------------------------------------------------- parity (lockstep, both modes)
// a: v33 host path (verification mode), b: device path (native v7 configuration: streaming digests, host record arrays),
// c: device path in the native v7.1 lane configuration (block-form digests, no host record arrays).
struct Pair {
    Setup& S; hc::ResidentOwner a, b, c; std::unique_ptr<gp::GpuMotionProvider> provider; const n::GapPicOptions o;
    Pair(Setup& s, const axial::State& seed, const hc::ChargeLedger& q, const n::GapPicOptions& opt, std::span<const n::ArenaRateEntry> census)
        : S(s), a(s.model, s.dev, seed, q, opt, census, options(gp::HostWork::v33)), b(s.model, s.dev, seed, q, opt, census, options(gp::HostWork::device)),
          c(s.model, s.dev, seed, q, opt, census, lane_options()), o(opt) {
        provider = std::make_unique<gp::GpuMotionProvider>(s.dev, a.motion_context(), std::max<std::size_t>(s.model->limits().maximum_carriers, 1));
    }
};
// One lockstep step; returns false (and records) on any difference. census is advanced on acceptance.
bool lockstep(Results& R, Pair& P, const std::string& name, long step, std::vector<n::ArenaRateEntry>& census, std::span<const hc::Newborn> newborns,
              double dt, J& steps_log) {
    const auto& g = P.S.model->mesh();
    const std::string src = hc::sha256("v34/source");
    const auto tables = tables_of(census);
    hc::StepInputs in{newborns, {}, !newborns.empty(), tables.states, tables.keys, tables.idx};
    auto ba = binding_of(P.a, *P.S.model, step, dt), bb = binding_of(P.b, *P.S.model, step, dt), bc = binding_of(P.c, *P.S.model, step, dt);
    std::optional<hc::ResidentJointProposal> pa, pb, pc;
    double ta = 0, tb = 0, tc = 0;
    auto t = steady::now();
    J ra = refusal_of([&] { pa.emplace(hc::prepare_closed_proposal(P.a, P.o, census, ba, src, *P.provider, in)); }); ta = ms_since(t); t = steady::now();
    J rb = refusal_of([&] { pb.emplace(hc::prepare_closed_proposal(P.b, P.o, census, bb, src, *P.provider, in)); }); tb = ms_since(t); t = steady::now();
    J rc = refusal_of([&] { pc.emplace(hc::prepare_closed_proposal(P.c, P.o, census, bc, src, *P.provider, in)); }); tc = ms_since(t);
    J ev{{"fixture", name}, {"step", step}, {"v33", ra}, {"device", rb}, {"lane", rc}, {"v33_ms", ta}, {"device_ms", tb}, {"lane_ms", tc}};
    if (ra != rb || ra != rc) { R.record("parity", name + " step " + std::to_string(step) + ": decision and refusal text", false, ev); steps_log.push_back(ev); return false; }
    if (!pa) { steps_log.push_back(ev); return true; }
    const auto &x = *pa, &y = *pb;
    const auto &mx = x.motion(), &my = y.motion();
    J c = J::object();
    c["candidate_field_sha256"] = x.candidate_field_sha256() == y.candidate_field_sha256();
    c["candidate_heat_sha256"] = x.candidate_heat_sha256() == y.candidate_heat_sha256();
    c["ledgers_sha256"] = x.field_heat().candidate_ledgers_sha256() == y.field_heat().candidate_ledgers_sha256();
    c["iterations"] = x.iterations() == y.iterations();
    c["combined_balance"] = bits(x.combined_balance_J()) == bits(y.combined_balance_J());
    c["census_v2"] = mx.census_sha256() == my.census_sha256();
    c["record_v7"] = mx.record_v7_sha256() == my.record_v7_sha256() && !mx.record_v7_sha256().empty();
    c["record_v7_host_reference"] = hc::MotionAssembly::record_v7(my) == my.record_v7_sha256();
    c["trial_field"] = mx.binding().trial_field_sha256 == my.binding().trial_field_sha256;
    c["budget"] = same_budget(mx.budget(), my.budget());
    c["work_scales"] = same_scales(mx.work_scales(), my.work_scales());
    c["sinks"] = same_sinks(mx.sinks(), my.sinks());
    c["births"] = same_births(mx.births(), my.births());
    c["endpoints"] = mx.endpoints().size() == my.endpoints().size() && (mx.endpoints().empty() || std::memcmp(mx.endpoints().data(), my.endpoints().data(), mx.endpoints().size() * 80) == 0);
    c["paths"] = mx.paths().size() == my.paths().size() && (mx.paths().empty() || std::memcmp(mx.paths().data(), my.paths().data(), mx.paths().size() * 64) == 0);
    c["current_initial_endpoint"] = same_bits(mx.current_A_m(), my.current_A_m()) && same_bits(mx.initial_charge_C(), my.initial_charge_C()) &&
                                    same_bits(mx.endpoint_charge_C(), my.endpoint_charge_C());
    const auto &sa = x.field_heat().solve_record(), &sb = y.field_heat().solve_record();
    c["independent_residual"] = bits(sa.independent_residual) == bits(sb.independent_residual) && bits(sa.rhs_norm) == bits(sb.rhs_norm) && bits(sa.target) == bits(sb.target);
    c["trial_E_B"] = same_bits(x.trial_E(), y.trial_E()) && same_bits(x.trial_B(), y.trial_B());
    // CPU vs GPU: deposit support (sealed builder on the v33 receipt vs the device array of the accepted iterate)
    {
        const auto ref = hc::DepositSupportBuilder::build(g, P.o, census, mx).absolute_C();
        const auto dev = P.b.device_step().support_download();
        long diff = 0; for (std::size_t k = 0; k < ref.size(); ++k) diff += bits(ref[k]) != bits(dev[k]);
        long nonzero = 0; for (double v : ref) nonzero += v != 0;
        c["deposit_support"] = diff == 0 && ref.size() == dev.size();
        ev["support_nodes_nonzero"] = nonzero; ev["support_nodes_differing"] = diff;
    }
    // bank transfer: sparse (8 CIC nodes) vs the sealed dense form
    {
        const auto ref = hc::bank_transfer(g, P.o, mx, P.S.model->exterior_support());
        const auto sp = gp::dense_of(g, gp::bank_transfer_sparse(g, P.o, my, P.S.model->exterior_support()));
        c["bank_transfer"] = ref.empty == sp.empty && (ref.empty || (same_words(ref.conductor_C, sp.conductor_C) && same_words(ref.exterior_C, sp.exterior_C) &&
                                                                       same_bits(ref.gross_C, sp.gross_C) && same_bits(ref.residue_C, sp.residue_C)));
    }
    // native v7.1 device lane view: block-form digests of the trial field and the current (device) vs the host reference
    {
        const auto& v = P.b.device_step().record_view();
        hc::Digest dt("total-midpoint-E-B-blocks-v1"); dt.text(host_root("native-v7.1/field-E", y.trial_E())); dt.text(host_root("native-v7.1/field-B", y.trial_B()));
        hc::Digest dj("motion-array-blocks-v1"); dj.text(host_root("native-v7.1/current-J", my.current_A_m()));
        c["v71_block_digests"] = v.trial_field_blocks_sha256 == dt.hex() && v.current_blocks_sha256 == dj.hex() && v.n_paths == mx.paths().size() &&
                                 v.n_endpoints == mx.endpoints().size() && v.n_sinks == mx.sinks().size() && v.n_births == mx.births().size();
    }
    // native v7.1 lane owner (block forms, no host arrays): same motion and state; its digests against the v7.1 CPU forms
    // on the host arrays of the device owner (same trial field and current), its record digest against a host recomputation
    const auto& z = *pc; const auto& mz = z.motion();
    {
        c["lane_state_digests"] = z.candidate_field_sha256() == x.candidate_field_sha256() && z.candidate_heat_sha256() == x.candidate_heat_sha256() &&
                                  z.field_heat().candidate_ledgers_sha256() == x.field_heat().candidate_ledgers_sha256() && z.iterations() == x.iterations() &&
                                  bits(z.combined_balance_J()) == bits(x.combined_balance_J());
        c["lane_motion"] = mz.census_sha256() == mx.census_sha256() && same_budget(mz.budget(), mx.budget()) && same_scales(mz.work_scales(), mx.work_scales()) &&
                           same_sinks(mz.sinks(), mx.sinks()) && same_births(mz.births(), mx.births());
        const auto tf = pf_v7::trial_field_blocks_sha256(y.trial_E(), y.trial_B()), cj = pf_v7::current_blocks_sha256(my.current_A_m());
        c["lane_trial_field_blocks"] = mz.binding().trial_field_sha256 == tf && z.midpoint().total_field_sha256() == tf;
        c["lane_current_blocks"] = mz.record_v7_array_digests()[0] == cj && mz.record_v7_array_digests()[1] == my.record_v7_array_digests()[1] &&
                                   mz.record_v7_array_digests()[2] == my.record_v7_array_digests()[2];
        c["lane_record_v7_host_reference"] = mz.record_v7_sha256() == hc::MotionAssembly::record_v7_forms(my, tf, cj) && mz.record_v7_sha256() != my.record_v7_sha256();
        c["lane_no_host_arrays"] = mz.endpoints().empty() && mz.paths().empty() && z.trial_E().empty() && z.trial_B().empty() && mz.current_A_m().empty();
        auto& st = P.c.device_step(); const auto& v = st.record_view();
        const auto cpuE = pf_v7::block_list(pf7::section_field_E, y.trial_E().data(), y.trial_E().size()).root;
        const auto cpuB = pf_v7::block_list(pf7::section_field_B, y.trial_B().data(), y.trial_B().size()).root;
        const auto cpuJ = pf_v7::block_list(pf7::section_current_J, my.current_A_m().data(), my.current_A_m().size()).root;
        c["lane_roots_forms_kernels_cpu"] = v.E_root == cpuE && v.B_root == cpuB && v.J_root == cpuJ && st.block_root(0, v.E, v.n_E) == cpuE &&
                                            st.block_root(1, v.B, v.n_B) == cpuB && st.block_root(2, v.J, v.n_E) == cpuJ;
        c["lane_view_counts"] = v.n_paths == mx.paths().size() && v.n_endpoints == mx.endpoints().size() && v.n_sinks == mx.sinks().size() &&
                                v.n_births == mx.births().size() && v.n_E == y.trial_E().size() && v.n_B == y.trial_B().size();
    }
    // preflight of both, then publication and the accepted states
    Boundary b1{x.candidate_field_sha256(), x.candidate_heat_sha256(), src, mx.binding().mesh_sha256, mx.binding().material_sha256, mx.binding().common_clock_sha256, x.end_s(), x.endpoint_revision()};
    Boundary b2{y.candidate_field_sha256(), y.candidate_heat_sha256(), src, my.binding().mesh_sha256, my.binding().material_sha256, my.binding().common_clock_sha256, y.end_s(), y.endpoint_revision()};
    t = steady::now(); const bool fa = P.a.preflight_final(x, census, b1); const double pfa = ms_since(t);
    t = steady::now(); const bool fb = P.b.preflight_final(y, census, b2); const double pfb = ms_since(t);
    Boundary b3{z.candidate_field_sha256(), z.candidate_heat_sha256(), src, mz.binding().mesh_sha256, mz.binding().material_sha256, mz.binding().common_clock_sha256, z.end_s(), z.endpoint_revision()};
    t = steady::now(); const bool fc = P.c.preflight_final(z, census, b3); const double pfc = ms_since(t);
    c["preflight_both"] = fa && fb;
    c["preflight_lane"] = fc;
    ev["preflight_ms"] = {pfa, pfb, pfc};
    auto next = hc::next_census(census, mx);
    auto next_b = hc::next_census(census, my);
    bool next_same = next.size() == next_b.size();
    for (std::size_t i = 0; next_same && i < next.size(); ++i) next_same = std::memcmp(&next[i].native, &next_b[i].native, 80) == 0;
    c["next_census"] = next_same;
    if (fa && fb && fc) {
        P.a.publish_prevalidated(std::move(*pa)); P.b.publish_prevalidated(std::move(*pb)); P.c.publish_prevalidated(std::move(*pc));
        axial::State s1, s2; hc::ChargeLedger q1, q2;
        P.a.accepted_state(s1, q1); P.b.accepted_state(s2, q2);
        c["accepted_state"] = same_bits(s1.E, s2.E) && same_bits(s1.B, s2.B) && same_bits(s1.E_low, s2.E_low) && same_bits(s1.B_low, s2.B_low) &&
                              same_words(q1.particles.net_C, q2.particles.net_C) && same_bits(q1.particles.gross_C, q2.particles.gross_C) &&
                              same_words(q1.conductor_C, q2.conductor_C) && same_words(q1.exterior_C, q2.exterior_C) && same_bits(q1.residue_C, q2.residue_C) &&
                              same_words(q1.gauss_admitted_C, q2.gauss_admitted_C);
        c["accepted_digests"] = P.a.field_sha256() == P.b.field_sha256() && P.a.heat_sha256() == P.b.heat_sha256() &&
                                P.a.kinetic_ledgers_sha256() == P.b.kinetic_ledgers_sha256();
        c["free_charge_mirror"] = same_words(P.a.free_charge().net_C, P.b.free_charge().net_C) && same_bits(P.a.free_charge().gross_C, P.b.free_charge().gross_C);
        c["lane_accepted_digests"] = P.c.field_sha256() == P.a.field_sha256() && P.c.heat_sha256() == P.a.heat_sha256() &&
                                     P.c.kinetic_ledgers_sha256() == P.a.kinetic_ledgers_sha256() && P.c.revision() == P.a.revision();
        census = std::move(next);
    }
    bool all = true; for (auto& [k, v] : c.items()) all = all && v.get<bool>();
    ev["checks"] = c; ev["sinks"] = mx.sinks().size(); ev["births"] = mx.births().size(); ev["carriers"] = mx.paths().size();
    ev["record_v7"] = mx.record_v7_sha256(); ev["field_sha256"] = x.candidate_field_sha256();
    ev["device_times"] = {{"census_gather", P.b.device_step().times().census_gather}, {"census_device", P.b.device_step().times().census_device},
                          {"motion_host", P.b.device_step().times().motion_host}, {"motion_wait", P.b.device_step().times().motion_wait},
                          {"record", P.b.device_step().times().record}, {"iterations", P.b.device_step().times().iterations}};
    steps_log.push_back(ev);
    if (!all) R.record("parity", name + " step " + std::to_string(step) + ": every compared item bitwise equal", false, ev);
    return all;
}

int parity(const std::string& trace, const std::filesystem::path& out, const std::string& nvrtc, const std::string& uuid, int transverse) {
    Results R;
    auto S = setup(trace, nvrtc, uuid, transverse, 20000);
    std::filesystem::create_directories(out);
    const auto& g = S.model->mesh(); const auto o = mf::geometry(g); const auto seed = S.model->dc_seed();
    J log = J::array();
    auto run = [&](const std::string& name, std::vector<n::ArenaRateEntry> census, long steps, double dt, std::vector<hc::Newborn> nb = {},
                   std::vector<n::ArenaRateEntry> bootstrap = {}) {
        if (bootstrap.empty()) bootstrap = census;
        std::unique_ptr<Pair> P;
        try { P = std::make_unique<Pair>(S, seed, declared(*S.model, seed, hc::charge_of_census(g, o, bootstrap)), o, bootstrap); }
        catch (const std::exception& e) { R.record("parity", name + ": bootstrap", false, {{"error", e.what()}}); return; }
        bool ok = true; long accepted = 0, refused = 0;
        for (long s = 1; s <= steps && ok; ++s) {
            const auto before = log.size();
            std::vector<hc::Newborn> step_nb = s == 1 ? nb : std::vector<hc::Newborn>{};
            ok = lockstep(R, *P, name, s, census, step_nb, dt, log);
            if (log.size() > before) { if (log.back()["v33"]["accepted"].get<bool>()) ++accepted; else { ++refused; break; } }
        }
        R.record("parity", name + ": decisions, refusal texts, records, digests and states bitwise equal in both host-work modes "
                 "and the native v7.1 lane configuration (block-form digests, no host arrays)", ok,
                 {{"accepted_steps", accepted}, {"refused", refused}});
    };
    run("pairs-2e3-w1e9", mf::pairs(2000, 11, 1e9), 4, 1e-14);
    run("pairs-1e4-w1", mf::pairs(10000, 17, 1.), 3, 1e-14);
    run("node-planes-500-w1", mf::node_plane_pairs(250, 3, 1.), 3, 1e-14);
    run("node-planes-500-w1e9", mf::node_plane_pairs(250, 5, 1e9), 3, 1e-14);
    {   // dense blob (coordinator v7 dense generator) with sink pairs
        std::mt19937_64 rng(20261101); std::uniform_real_distribution<double> u(0, 1); std::normal_distribution<double> nd(0, 1);
        const double se = std::sqrt(3 * n::pic_e / mf::me), si = std::sqrt(.025 * n::pic_e / 2.9915e-26);
        std::vector<n::ArenaRateEntry> blob;
        for (long k = 0; k < 400; ++k) {
            const n::PicVector x{.0026 + 8e-4 * u(rng), .0051 + 8e-4 * u(rng), -.003 + .006 * u(rng)};
            blob.push_back(mf::row(std::uint64_t(2 * k + 1), "be" + std::to_string(k), "electron", x, {se * nd(rng), se * nd(rng), se * nd(rng)}, -n::pic_e, mf::me, 3.45e10));
            blob.push_back(mf::row(std::uint64_t(2 * k + 2), "bi" + std::to_string(k), "H2O+", x, {si * nd(rng), si * nd(rng), si * nd(rng)}, n::pic_e, mf::ion_mass(), 3.45e10));
        }
        std::uint64_t id = blob.back().native.id + 1;
        for (int k = 0; k < 20; ++k) {
            const n::PicVector x{.0035 + 1e-4 * k, .002 + 2e-8 * (k + 1), -.002 + 2e-4 * k};
            auto ue = mf::electron_u(100);
            blob.push_back(mf::row(id++, "se" + std::to_string(k), "electron", x, {.1 * ue[0], -ue[0], 0}, -n::pic_e, mf::me, 3.45e10));
            blob.push_back(mf::row(id++, "si" + std::to_string(k), "H2O+", x, {0, 300, 0}, n::pic_e, mf::ion_mass(), 3.45e10));
        }
        run("dense-blob-sinks", blob, 8, 1e-14);
    }
    // sink and newborn fixtures (field 0 only; one step each)
    for (auto& f : sf::fixture_set()) {
        if (f.field != 0) continue;
        std::vector<hc::Newborn> nb; std::vector<n::ArenaRateEntry> boot;
        for (const auto& [id, frac] : f.newborns) nb.push_back({id, seed.time_s + frac * 1e-14});
        if (!nb.empty()) { std::set<std::uint64_t> ids; for (auto& x : nb) ids.insert(x.native_id); for (auto& r : f.census) if (!ids.contains(r.native.id)) boot.push_back(r); }
        run("sinks/" + f.name, f.census, f.newborns.empty() ? 2 : 1, 1e-14, nb, boot);
    }
    // refusal fixtures: the same text at the same point in both modes
    {
        const double cx = g.side / 2, cy = g.side / 2;
        auto base = mf::pairs(200, 23, 1.);
        auto dup = base; dup.push_back(dup[5]);                                            // duplicate ID
        auto empty_packet = base; empty_packet[7].packet_id.clear();                       // empty packet text
        auto divergent = base; divergent[9].packet_id = divergent[8].packet_id;            // divergent Count2 packet
        auto unsorted = base; std::swap(unsorted[3], unsorted[4]);                         // admitted after sorting
        auto in_foil = base; in_foil.push_back(mf::row(10000, "m-e", "electron", {cx, cy - .0026, 0.}, mf::electron_u(10), -n::pic_e, mf::me, 1.));
        for (auto& [name, cen] : std::vector<std::pair<std::string, std::vector<n::ArenaRateEntry>>>{
                 {"refusal-duplicate-id", dup}, {"refusal-empty-packet", empty_packet}, {"refusal-divergent-packet", divergent},
                 {"unsorted-census", unsorted}, {"refusal-carrier-in-foil", in_foil}}) {
            // bootstrap from the valid base census; the defective census is presented at the step
            std::unique_ptr<Pair> P = std::make_unique<Pair>(S, seed, declared(*S.model, seed, hc::charge_of_census(g, o, base)), o, base);
            auto c = cen;
            const auto before = log.size();
            const bool ok = lockstep(R, *P, name, 1, c, {}, 1e-14, log);
            R.record("refusals", name + ": identical decision and refusal text in both modes", ok && log.size() > before,
                     log.size() > before ? log.back() : J::object());
        }
        // newborn list defects
        auto f = sf::fixture_set();
        for (auto& x : f) if (x.name == "newborns-faces-and-interior") {
            std::vector<n::ArenaRateEntry> boot; for (auto& r : x.census) if (r.native.id <= 2) boot.push_back(r);
            std::vector<hc::Newborn> all;
            for (const auto& [id, frac] : x.newborns) all.push_back({id, seed.time_s + frac * 1e-14});
            auto twice = all; twice.push_back(all[0]);
            auto absent = all; absent.push_back({99, seed.time_s});
            auto late = all; late[1].birth_time_s = seed.time_s + 2e-14;
            for (auto& [name, nb] : std::vector<std::pair<std::string, std::vector<hc::Newborn>>>{
                     {"refusal-newborn-listed-twice", twice}, {"refusal-newborn-absent", absent}, {"refusal-newborn-time", late}}) {
                std::unique_ptr<Pair> P = std::make_unique<Pair>(S, seed, declared(*S.model, seed, hc::charge_of_census(g, o, boot)), o, boot);
                auto c = x.census;
                const auto before = log.size();
                const bool ok = lockstep(R, *P, name, 1, c, nb, 1e-14, log);
                R.record("refusals", name + ": identical decision and refusal text in both modes", ok && log.size() > before,
                         log.size() > before ? log.back() : J::object());
            }
        }
        // a binding outside the configured interval
        {
            std::unique_ptr<Pair> P = std::make_unique<Pair>(S, seed, declared(*S.model, seed, hc::charge_of_census(g, o, base)), o, base);
            auto c = base;
            const auto before = log.size();
            const bool ok = lockstep(R, *P, "refusal-interval", 1, c, {}, 2e-14, log);
            R.record("refusals", "refusal-interval: identical decision and refusal text in both modes", ok && log.size() > before,
                     log.size() > before ? log.back() : J::object());
        }
    }
    std::ofstream(out / "parity-steps.json") << log.dump(1) << '\n';
    std::ofstream(out / "parity.json") << J{{"cases", R.cases}, {"failed", R.failed}, {"transverse", transverse}}.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// ---------------------------------------------------------------- tamper (device preflight)
int tamper(const std::string& trace, const std::filesystem::path& out, const std::string& nvrtc, const std::string& uuid) {
    Results R;
    auto S = setup(trace, nvrtc, uuid, 45, 20000);
    std::filesystem::create_directories(out);
    const auto& g = S.model->mesh(); const auto o = mf::geometry(g); const auto seed = S.model->dc_seed();
    std::vector<n::ArenaRateEntry> census = mf::pairs(2000, 29, 1e9);
    {   // sink pairs so the sink section is not empty
        std::uint64_t id = census.back().native.id + 1;
        for (int k = 0; k < 10; ++k) {
            const n::PicVector x{.0035 + 1e-4 * k, .002 + 2e-8 * (k + 1), -.002 + 2e-4 * k};
            auto ue = mf::electron_u(100);
            census.push_back(mf::row(id++, "se" + std::to_string(k), "electron", x, {.1 * ue[0], -ue[0], 0}, -n::pic_e, mf::me, 1e9));
            census.push_back(mf::row(id++, "si" + std::to_string(k), "H2O+", x, {0, 300, 0}, n::pic_e, mf::ion_mass(), 1e9));
        }
    }
    hc::ResidentOwner res(S.model, S.dev, seed, declared(*S.model, seed, hc::charge_of_census(g, o, census)), o, census, options(gp::HostWork::device));
    gp::GpuMotionProvider unused(S.dev, res.motion_context(), 1);
    const std::string src = hc::sha256("v34/source");
    auto fresh = [&](long step) {
        const auto t = tables_of(census);
        hc::StepInputs in{{}, {}, false, t.states, t.keys, t.idx};
        return hc::prepare_closed_proposal(res, o, census, binding_of(res, *S.model, step, 1e-14), src, unused, in);
    };
    auto bd = [&](const hc::ResidentJointProposal& p) {
        const auto& b = p.motion().binding();
        return Boundary{p.candidate_field_sha256(), p.candidate_heat_sha256(), src, b.mesh_sha256, b.material_sha256, b.common_clock_sha256, p.end_s(), p.endpoint_revision()};
    };
    {
        auto p = fresh(1);
        R.record("tamper", "control: unmodified device proposal admitted at preflight", res.preflight_final(p, census, bd(p)),
                 {{"sinks", p.motion().sinks().size()}});
        auto& step = res.device_step();
        struct Case { const char* name; int section; std::size_t word; std::uint64_t mask; };
        for (const auto& c : {Case{"one ulp of an endpoint position word", 0, 11, 1}, Case{"one ulp of a path displacement word", 1, 8 * 7 + 4, 1},
                              Case{"one bit of a sink record hit point", 2, 9, 1}, Case{"one ulp of a census row velocity (device census)", 4, 10 * 5 + 4, 1}}) {
            step.tamper_record_word(c.section, c.word, c.mask);
            const bool pre = res.preflight_final(p, census, bd(p));
            step.tamper_record_word(c.section, c.word, c.mask);   // restore
            R.record("tamper", std::string(c.name) + " refused at preflight", !pre);
        }
        R.record("tamper", "control after restoring every word: admitted", res.preflight_final(p, census, bd(p)));
        auto copy = census;
        R.record("tamper", "a different census span at preflight is refused", !res.preflight_final(p, copy, bd(p)));
        {   // candidate slot free charge: one ulp of one net word
            axial::State st; hc::ChargeLedger q; res.candidate_state(st, q, p.field_heat());
            int node = -1; for (int k = 0; k < g.nn && node < 0; ++k) if (q.particles.net_C[k].hi != 0) node = k;
            q.particles.net_C[node].hi = std::nextafter(q.particles.net_C[node].hi, INFINITY);
            res.device_field().upload(1, st, q);
            R.record("tamper", "one ulp of a candidate free-charge word refused at preflight", !res.preflight_final(p, census, bd(p)), {{"node", node}});
        }
    }
    {   // native v7.1 lane configuration: block-form trial-field/current digests, no host record arrays
        hc::ResidentOwner lane(S.model, S.dev, seed, declared(*S.model, seed, hc::charge_of_census(g, o, census)), o, census, lane_options());
        const auto t = tables_of(census);
        hc::StepInputs in{{}, {}, false, t.states, t.keys, t.idx};
        auto p = hc::prepare_closed_proposal(lane, o, census, binding_of(lane, *S.model, 1, 1e-14), src, unused, in);
        auto& st = lane.device_step();
        R.record("tamper-lane", "control: unmodified lane proposal (block forms) admitted at preflight", lane.preflight_final(p, census, bd(p)),
                 {{"trial_field_sha256", p.motion().binding().trial_field_sha256}, {"current_sha256", p.motion().record_v7_array_digests()[0]}});
        const auto Jw = st.record_words(5);
        std::size_t k = 0; while (k < Jw.size() && Jw[k] == 0) ++k;
        struct Case { const char* name; int section; std::size_t word; };
        for (const auto& c : {Case{"one ulp of the current J (device array of the record)", 5, k}, Case{"one ulp of an endpoint position word", 0, 11},
                              Case{"one ulp of a path displacement word", 1, 8 * 7 + 4}}) {
            st.tamper_record_word(c.section, c.word, 1);
            const bool pre = lane.preflight_final(p, census, bd(p));
            st.tamper_record_word(c.section, c.word, 1);
            R.record("tamper-lane", std::string(c.name) + " refused at preflight", !pre, {{"word", c.word}});
        }
        R.record("tamper-lane", "control after restoring every word: admitted", lane.preflight_final(p, census, bd(p)));
    }
    {   // stale proposal: a later prepare overwrites the device record arrays and the candidate slot
        auto p1 = fresh(2);
        auto p2 = fresh(3);
        R.record("tamper", "a proposal superseded by a later prepare is refused", !res.preflight_final(p1, census, bd(p1)));
        R.record("tamper", "the later proposal is admitted", res.preflight_final(p2, census, bd(p2)));
    }
    std::ofstream(out / "tamper.json") << J{{"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// ---------------------------------------------------------------- bench (load fixture)
std::vector<n::ArenaRateEntry> load_fixture(std::size_t N, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> u(0, 1);
    std::normal_distribution<double> nd(0, 1);
    const double si = std::sqrt(.025 * n::pic_e / mf::mn);
    std::vector<n::ArenaRateEntry> out; out.reserve(N);
    const std::size_t neutrals = N / 50, pairs = (N - neutrals) / 2;
    std::uint64_t id = 1;
    for (std::size_t k = 0; k < pairs; ++k) {
        const n::PicVector x{.0035 + .002 * u(rng), .0025 + .004 * u(rng), -.004 + .008 * u(rng)};
        const double eV = 1 + 299 * u(rng), ct = 2 * u(rng) - 1, ph = 6.283185307179586 * u(rng), st = std::sqrt(std::max(0., 1 - ct * ct));
        const double xe = eV * n::pic_e / (mf::me * n::pic_c * n::pic_c), uu = n::pic_c * std::sqrt(xe * (xe + 2));
        const std::string tag = "L" + std::to_string(k);
        out.push_back(mf::row(id++, tag + "e", "electron", x, {uu * st * std::cos(ph), uu * st * std::sin(ph), uu * ct}, -n::pic_e, mf::me, 1.));
        out.push_back(mf::row(id++, tag + "i", "H2O+", x, {si * nd(rng), si * nd(rng), si * nd(rng)}, n::pic_e, mf::ion_mass(), 1.));
    }
    for (std::size_t k = 0; out.size() < N; ++k) {
        const n::PicVector x{.0035 + .002 * u(rng), .0025 + .004 * u(rng), -.004 + .008 * u(rng)};
        out.push_back(mf::row(id++, "W" + std::to_string(k), "H2O", x, {si * nd(rng), si * nd(rng), si * nd(rng)}, 0., mf::mn, 1e6));
    }
    return out;
}
// The next census from the accepted endpoints (test side, not timed): the receipt's host arrays, or the device record
// arrays (device lane without host record arrays).
template <class Proposal> std::vector<n::ArenaRateEntry> next_of(const std::vector<n::ArenaRateEntry>& census, const Proposal& p, gp::DeviceStep* lane) {
    if (!lane || p.motion().host_arrays()) return hc::next_census(census, p.motion());
    const auto w = lane->record_words(0);
    std::vector<n::ArenaRateEntry> out; out.reserve(w.size() / 10);
    std::size_t j = 0;
    for (const auto& row : census) {
        if (j < w.size() / 10 && w[10 * j] == row.native.id) { auto x = row; std::memcpy(&x.native, &w[10 * j], 80); out.push_back(std::move(x)); ++j; }
    }
    return out;
}
// Device census of the lane bench (the role of native v7.1's trial.census_v7_device(): rows, packet keys and state
// indices on the shared session, ascending native ID), uploaded by the test before each step (not timed).
struct DeviceCensus {
    hall_gpu::Session& s; hall_gpu::DeviceBuffer rows, keys, states; hall_gpu::HostBuffer stage;
    DeviceCensus(hall_gpu::Session& session, std::size_t N)
        : s(session), rows(s.alloc(std::max<std::size_t>(N, 1) * 80, "bench census rows")), keys(s.alloc(std::max<std::size_t>(N, 1) * 8, "bench census keys")),
          states(s.alloc(std::max<std::size_t>(N, 1) * 8, "bench census states")), stage(s.alloc_pinned(std::max<std::size_t>(N, 1) * 80, "bench census stage")) {}
    ~DeviceCensus() { s.release(stage); s.release(states); s.release(keys); s.release(rows); }
    void upload(const std::vector<n::ArenaRateEntry>& census, const Tables& t) {
        const std::size_t n = census.size();
        for (std::size_t i = 0; i < n; ++i) std::memcpy(stage.as<char>() + 80 * i, &census[i].native, 80);
        if (n) s.upload(rows, 0, stage, 0, n * 80);
        s.synchronize("bench census rows");
        if (n) { std::memcpy(stage.pointer, t.keys.data(), n * 8); s.upload(keys, 0, stage, 0, n * 8); }
        s.synchronize("bench census keys");
        if (n) { std::memcpy(stage.pointer, t.idx.data(), n * 8); s.upload(states, 0, stage, 0, n * 8); }
        s.synchronize("bench census states");
    }
};
template <class OwnerT, class Proposal>
J bench_step(OwnerT& owner, const hc::Model& model, const n::GapPicOptions& o, std::vector<n::ArenaRateEntry>& census, hc::MotionProvider& provider, long step,
             gp::DeviceStep* lane = nullptr, DeviceCensus* dc = nullptr) {
    const std::string src = hc::sha256("v34/source");
    const auto tables = tables_of(census);
    hc::StepInputs in{{}, {}, false, tables.states, tables.keys, tables.idx};
    if (dc) {   // native v7.1 configuration: the census is read in place on the device
        dc->upload(census, tables);
        in.census_packet_keys = {}; in.census_state_indices = {};
        in.device_rows = dc->rows.address; in.device_keys = dc->keys.address; in.device_states = dc->states.address; in.device_count = census.size();
        in.packet_of = [&census](std::uint64_t id) {
            const auto it = std::lower_bound(census.begin(), census.end(), id, [](const n::ArenaRateEntry& r, std::uint64_t x) { return r.native.id < x; });
            hc::need(it != census.end() && it->native.id == id, "bench packet_of: unknown carrier");
            return std::make_pair(it->packet_id, it->state);
        };
    }
    const auto b = binding_of(owner, model, step, 1e-14);
    double c0 = cpu_ms(); auto t0 = steady::now();
    Proposal p = hc::prepare_closed_proposal(owner, o, census, b, src, provider, in);
    const double prep_wall = ms_since(t0), prep_cpu = cpu_ms() - c0;
    const auto& bb = p.motion().binding();
    Boundary bd{p.candidate_field_sha256(), p.candidate_heat_sha256(), src, bb.mesh_sha256, bb.material_sha256, bb.common_clock_sha256, p.end_s(), p.endpoint_revision()};
    c0 = cpu_ms(); t0 = steady::now();
    const bool ok = owner.preflight_final(p, census, bd);
    const double pf_wall = ms_since(t0), pf_cpu = cpu_ms() - c0;
    hc::need(ok, "bench preflight");
    const auto& T = p.timings();
    J j{{"step", step}, {"carriers", census.size()}, {"iterations", p.iterations()},
        {"prepare_wall_ms", prep_wall}, {"prepare_cpu_ms", prep_cpu}, {"preflight_wall_ms", pf_wall}, {"preflight_cpu_ms", pf_cpu},
        {"admit_census_ms", T.admit_census_s * 1e3}, {"motion_ms", T.motion_s * 1e3}, {"trial_digest_ms", T.motion_trial_digest_s * 1e3},
        {"support_and_bank_ms", T.support_s * 1e3}, {"field_prepare_ms", T.field_prepare_s * 1e3}, {"solve_ms", T.solve_s * 1e3},
        {"residual_check_ms", T.residual_check_s * 1e3}, {"record_seal_ms", T.motion_v2_seal_s * 1e3}, {"field_seal_ms", T.seal_s * 1e3},
        {"heat_ms", T.heat_s * 1e3}, {"total_ms", T.total_s * 1e3}, {"provider", T.motion_provider},
        {"cpu", {{"trial_digest_ms", T.trial_digest_cpu_s * 1e3}, {"support_and_bank_ms", T.support_cpu_s * 1e3}, {"field_prepare_ms", T.field_prepare_cpu_s * 1e3},
                 {"finish_ms", T.finish_cpu_s * 1e3}}}};
    auto next = next_of(census, p, lane);
    c0 = cpu_ms(); t0 = steady::now();
    owner.publish_prevalidated(std::move(p));
    j["publish_wall_ms"] = ms_since(t0); j["publish_cpu_ms"] = cpu_ms() - c0;
    census = std::move(next);
    return j;
}
int bench(const std::string& trace, const std::filesystem::path& out, const std::string& nvrtc, const std::string& uuid, int transverse, std::size_t N,
          long steps, const std::string& mode) {
    std::filesystem::create_directories(out);
    const auto t_setup = steady::now();
    auto S = setup(trace, nvrtc, uuid, transverse, N, mode == "owner");
    const auto& model = mode == "owner" ? *S.owner_model : *S.model;
    const auto& g = model.mesh(); const auto o = mf::geometry(g); const auto seed = model.dc_seed();
    auto census = load_fixture(N, 20261031);
    const auto q = declared(model, seed, hc::charge_of_census(g, o, census, model.limits().maximum_carriers));
    J rows = J::array();
    J meta{{"mode", mode}, {"transverse", transverse}, {"carriers", N}, {"steps", steps}, {"nn", g.nn}, {"ne", g.ne}, {"nb", g.nb}};
    if (mode == "owner") {
        hc::Owner owner(S.owner_model, seed, q, o, census);
        gp::GpuMotionProvider provider(S.dev, owner.motion_context(), N);
        meta["setup_s"] = ms_since(t_setup) / 1e3;
        for (long s = 1; s <= steps; ++s) {
            auto j = bench_step<hc::Owner, hc::JointProposal>(owner, model, o, census, provider, s);
            const auto& st = provider.stats();
            j["push_upload_ms"] = st.upload_carriers_s * 1e3; j["push_admit_ms"] = st.admit_s * 1e3; j["push_device_ms"] = st.device_run_s * 1e3;
            j["push_assemble_ms"] = st.assemble_s * 1e3; j["push_gpu_total_ms"] = st.last_gpu_total_ms;
            rows.push_back(j); std::cout << j.dump() << std::endl;
        }
    } else {
        auto opt = options(mode == "v33" ? gp::HostWork::v33 : gp::HostWork::device);
        if (mode == "lane") opt = lane_options();   // native v7.1 device lane: no host record arrays, block-form field/current digests
        hc::ResidentOwner owner(S.model, S.dev, seed, q, o, census, opt);
        gp::GpuMotionProvider provider(S.dev, owner.motion_context(), mode == "v33" ? N : 1);
        std::unique_ptr<DeviceCensus> dc;
        if (mode == "lane") dc = std::make_unique<DeviceCensus>(S.dev->session(), N);
        meta["setup_s"] = ms_since(t_setup) / 1e3;
        for (long s = 1; s <= steps; ++s) {
            if (mode != "v33") owner.device_step().reset_times();
            auto j = bench_step<hc::ResidentOwner, hc::ResidentJointProposal>(owner, model, o, census, provider, s, mode == "lane" ? &owner.device_step() : nullptr,
                                                                              dc.get());
            if (mode == "v33") {
                const auto& st = provider.stats();
                j["push_upload_ms"] = st.upload_carriers_s * 1e3; j["push_admit_ms"] = st.admit_s * 1e3; j["push_device_ms"] = st.device_run_s * 1e3;
                j["push_assemble_ms"] = st.assemble_s * 1e3; j["push_gpu_total_ms"] = st.last_gpu_total_ms;
            } else {
                const auto& d = owner.device_step().times();
                j["device_step"] = {{"census_gather_ms", d.census_gather}, {"census_device_ms", d.census_device}, {"motion_host_ms", d.motion_host},
                                    {"motion_wait_ms", d.motion_wait}, {"record_ms", d.record}, {"push_gpu_total_ms", d.device_push_ms}, {"iterations", d.iterations},
                                    {"census_cpu_ms", d.census_cpu}, {"motion_cpu_ms", d.motion_cpu}, {"record_cpu_ms", d.record_cpu}, {"preflight_cpu_ms", d.preflight_cpu}};
            }
            rows.push_back(j); std::cout << j.dump() << std::endl;
        }
    }
    std::ofstream(out / "bench.json") << J{{"meta", meta}, {"steps", rows}}.dump(1) << '\n';
    return 0;
}

// ---------------------------------------------------------------- block-form digests
int blocks(const std::string& trace, const std::filesystem::path& out, const std::string& nvrtc, const std::string& uuid) {
    Results R;
    auto S = setup(trace, nvrtc, uuid, 45, 1000);
    std::filesystem::create_directories(out);
    const auto& g = S.model->mesh(); const auto o = mf::geometry(g); const auto seed = S.model->dc_seed();
    auto census = mf::pairs(200, 3, 1.);
    hc::ResidentOwner res(S.model, S.dev, seed, declared(*S.model, seed, hc::charge_of_census(g, o, census)), o, census, options(gp::HostWork::device));
    auto& step = res.device_step();

    hc::Vec E, B;
    res.device_field().accepted_rounded(E, B);
    res.device_field().begin_fixed_point();
    const auto tr = res.device_field().trial_addresses();
    const auto dE = step.block_root(0, tr[0], E.size()), dB = step.block_root(1, tr[1], B.size());
    const auto hE = host_root("native-v7.1/field-E", E), hB = host_root("native-v7.1/field-B", B);
    R.record("blocks", "field-E block root: device equals the host reference", dE == hE, {{"device", dE}, {"host", hE}, {"n", E.size()}});
    R.record("blocks", "field-B block root: device equals the host reference", dB == hB, {{"device", dB}, {"host", hB}, {"n", B.size()}});
    // native v7.1 forms library: Forms::block_lists (device) and pf_v7::block_list (CPU) give the same roots; block digests
    const auto fE = step.block_root_forms(0, tr[0], E.size()), fB = step.block_root_forms(1, tr[1], B.size());
    const auto cE = pf_v7::block_list(pf7::section_field_E, E.data(), E.size()).root, cB = pf_v7::block_list(pf7::section_field_B, B.data(), B.size()).root;
    R.record("blocks", "field-E/B roots: step kernels = Forms::block_lists (native v7.1) = pf_v7::block_list (CPU)", dE == fE && fE == cE && dB == fB && fB == cB,
             {{"forms", {fE, fB}}, {"cpu", {cE, cB}}});
    hc::Digest dt("total-midpoint-E-B-blocks-v1"); dt.text(dE); dt.text(dB);
    R.record("blocks", "trial-field block digest = pf_v7::trial_field_blocks_sha256 (native v7.1)", dt.hex() == pf_v7::trial_field_blocks_sha256(E, B));
    {   // current-J section on an arbitrary array of the edge count (the trial E), and the empty array
        const auto dJ = step.block_root(2, tr[0], E.size()), fJ = step.block_root_forms(2, tr[0], E.size());
        const auto cJ = pf_v7::block_list(pf7::section_current_J, E.data(), E.size()).root;
        hc::Digest dj("motion-array-blocks-v1"); dj.text(dJ);
        R.record("blocks", "current-J root and digest: step kernels = Forms::block_lists = CPU; = pf_v7::current_blocks_sha256",
                 dJ == fJ && fJ == cJ && dj.hex() == pf_v7::current_blocks_sha256(std::span<const double>(E)));
        const auto d0 = step.block_root(2, tr[0], 0), f0 = step.block_root_forms(2, tr[0], 0);
        R.record("blocks", "empty current-J array: step kernels = Forms::block_lists = CPU", d0 == f0 && f0 == pf_v7::block_list(pf7::section_current_J, E.data(), 0).root);
    }
    std::ofstream(out / "blocks.json") << J{{"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    return R.failed ? 1 : 0;
}

// ---------------------------------------------------------------- host-loop timing in a device process
int hostloop(const std::string& trace, const std::filesystem::path& out, const std::string& nvrtc, const std::string& uuid, std::size_t N, int reps) {
    std::filesystem::create_directories(out);
    auto cpu = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto& g = cpu->mesh(); const auto o = mf::geometry(g);
    auto census = mf::pairs(N, 41, 1.);
    std::vector<n::GapParticle> ends; ends.reserve(census.size());
    std::mt19937_64 rng(3); std::uniform_real_distribution<double> u(-1, 1);
    for (const auto& r : census) { auto p = r.native; for (int d = 0; d < 3; ++d) p.position_m[d] += 1e-6 * u(rng); ends.push_back(p); }
    std::vector<double> big(std::size_t(1) << 22, 1.0);
    const auto census0 = census;
    auto median = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size() / 2]; };
    auto measure = [&](const std::string& phase) {
        std::vector<double> sup, arith, mem, sup_cpu;
        for (int k = 0; k < reps; ++k) {
            double c0 = cpu_ms(); auto t = steady::now();
            auto s = hc::DepositSupportBuilder::build(g, o, census0, ends);
            sup.push_back(ms_since(t)); sup_cpu.push_back(cpu_ms() - c0);
            hc::need(s.absolute_C().size() == std::size_t(g.nn), "support");
            t = steady::now(); volatile double x = 1; for (int i = 0; i < 2000000; ++i) x = x * 1.0000001 + 1e-9; arith.push_back(ms_since(t));
            t = steady::now(); double acc = 0; for (double v : big) acc += v; volatile double sink = acc; (void)sink; mem.push_back(ms_since(t));
        }
        return J{{"phase", phase}, {"support_ms", median(sup)}, {"support_cpu_ms", median(sup_cpu)}, {"arith_2e6_ms", median(arith)}, {"stream_32MB_ms", median(mem)},
                 {"support_min_ms", *std::min_element(sup.begin(), sup.end())}};
    };
    J rows = J::array();
    rows.push_back(measure("no CUDA context")); std::cout << rows.back().dump() << std::endl;
    auto dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    rows.push_back(measure("device context (blocking sync)")); std::cout << rows.back().dump() << std::endl;
    const auto desc = hc::resident_descriptor_sha256(dev, *cpu);
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<hc::ResidentBackend>(desc), limits_for(N), hc::ThermalDeclaration{}, 1e-14,
                                                   hc::MotionForms{true, 1e-14});
    const auto seed = model->dc_seed();
    hc::ResidentOwner res(model, dev, seed, declared(*model, seed, hc::charge_of_census(g, o, census, N)), o, census, options(gp::HostWork::v33));
    gp::GpuMotionProvider provider(dev, res.motion_context(), N);
    rows.push_back(measure("resident owner and push engine allocated")); std::cout << rows.back().dump() << std::endl;
    for (long s = 1; s <= 3; ++s) (void)bench_step<hc::ResidentOwner, hc::ResidentJointProposal>(res, *model, o, census, provider, s);
    rows.push_back(measure("after three coupled device steps")); std::cout << rows.back().dump() << std::endl;
    std::ofstream(out / "hostloop.json") << J{{"carriers", N}, {"reps", reps}, {"rows", rows}}.dump(1) << '\n';
    return 0;
}

// ---------------------------------------------------------------- divB gross ledger
int divb(const std::string& trace, const std::filesystem::path& out, const std::string& nvrtc, const std::string& uuid, long steps) {
    Results R;
    std::filesystem::create_directories(out);
    auto dev = std::make_shared<gp::Device>(gp::DeviceOptions{uuid, nvrtc});
    auto cpu = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace);
    const auto desc = hc::resident_descriptor_sha256(dev, *cpu);
    auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, trace, std::make_shared<hc::ResidentBackend>(desc));
    const auto& g = model->mesh(); const auto o = mf::geometry(g); const auto seed = model->dc_seed();
    std::vector<n::ArenaRateEntry> none;
    const auto q0 = declared(*model, seed, hc::charge_of_census(g, o, none));
    hc::Owner host(cpu, seed, q0, o, none);
    auto ro = options(gp::HostWork::device); ro.solve = gp::ResidentOptions::Solve::host_cpu_pcg;
    hc::ResidentOwner res(model, dev, seed, q0, o, none, ro);
    {   // the seed in coefficient form
        const auto c = model->admit(seed, q0);
        R.record("divb", "DC seed admitted in coefficient form; largest |divB|/bound stated", true,
                 {{"max_divB_T_m", c.max_divB_T_m}, {"max_divB_ratio", c.max_divB_ratio}, {"max_bound_T_m", c.max_divB_bound_T_m}});
    }
    long equal_G = 0, equal_state = 0; double max_floor = 0, max_ratio = 0; J rows = J::array();
    for (long k = 1; k <= steps; ++k) {
        auto a = host.prepare_vacuum(1e-14);
        Boundary ba{a.candidate_field_sha256(), a.candidate_heat_sha256(), "", "", "", "", a.end_s(), a.endpoint_revision()};
        hc::need(host.preflight_vacuum_final(a, ba), "host vacuum preflight");
        host.publish_vacuum_prevalidated(std::move(a));
        auto b = res.prepare_vacuum(1e-14);
        Boundary bb{b.candidate_field_sha256(), b.candidate_heat_sha256(), "", "", "", "", b.end_s(), b.endpoint_revision()};
        hc::need(res.preflight_vacuum_final(b, bb), "resident vacuum preflight");
        res.publish_vacuum_prevalidated(std::move(b));
        axial::State s; hc::ChargeLedger q; res.accepted_state(s, q);
        const auto& hq = host.charge();
        equal_G += same_bits(hq.divb_gross_T_m, q.divb_gross_T_m);
        equal_state += same_bits(host.field().B, s.B) && same_bits(host.field().B_low, s.B_low) && same_bits(host.field().E, s.E);
        const auto c = model->admit(s, q);
        max_floor = std::max(max_floor, c.max_divB_floor_T_m); max_ratio = std::max(max_ratio, c.max_divB_ratio);
        if (k % std::max(1L, steps / 10) == 0) rows.push_back({{"step", k}, {"max_divB_ratio", c.max_divB_ratio}, {"max_floor_T_m", c.max_divB_floor_T_m},
                                                               {"max_G_T_m", *std::max_element(q.divb_gross_T_m.begin(), q.divb_gross_T_m.end())}});
    }
    R.record("divb", "host Model::prepare and device k_divb_gross give the same G_c bit for bit every step", equal_G == steps,
             {{"steps", steps}, {"equal_G", equal_G}, {"equal_state", equal_state}, {"max_floor_T_m", max_floor}, {"max_divB_ratio", max_ratio}, {"series", rows}});
    // one-ulp tamper of G_c
    axial::State s; hc::ChargeLedger q; res.accepted_state(s, q);
    std::size_t cell = 0; for (std::size_t c = 0; c < q.divb_gross_T_m.size(); ++c) if (q.divb_gross_T_m[c] > q.divb_gross_T_m[cell]) cell = c;
    {   // device accepted slot
        auto tq = q; tq.divb_gross_T_m[cell] = std::nextafter(tq.divb_gross_T_m[cell], INFINITY);
        res.device_field().upload(0, s, tq);
        auto p = res.prepare_vacuum(1e-14);
        Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
        const bool pre = res.preflight_vacuum_final(p, b);
        res.device_field().upload(0, s, q);
        R.record("divb-tamper", "one ulp of G_c in the accepted device slot refused at preflight (leaves recomputed)", !pre, {{"cell", cell}});
    }
    auto edit_checkpoint = [&](const std::filesystem::path& dir) {
        auto j = J::parse(std::ifstream(dir / "checkpoint.json"));
        auto& arr = j.at("divb_gross_T_m");
        double v = arr[cell].get<double>(); arr[cell] = std::nextafter(v, INFINITY);
        std::ofstream(dir / "checkpoint.json") << j.dump() << '\n';
    };
    {   // resident checkpoint without reseal
        const auto dir = out / "resident-ck"; std::filesystem::remove_all(dir); res.save(dir);
        edit_checkpoint(dir);
        std::string why;
        try { (void)hc::ResidentOwner::restore(model, dev, o, none, dir, ro); } catch (const std::exception& e) { why = e.what(); }
        R.record("divb-tamper", "one ulp of G_c in a resident checkpoint (not resealed) refused at restore", why == "checkpoint field/charge seal mismatch", {{"refusal", why}});
    }
    {   // host checkpoint without reseal
        const auto dir = out / "host-ck"; std::filesystem::remove_all(dir); host.save(dir);
        edit_checkpoint(dir);
        std::string why;
        try { (void)hc::Owner::restore(cpu, o, none, dir); } catch (const std::exception& e) { why = e.what(); }
        R.record("divb-tamper", "one ulp of G_c in a host checkpoint (not resealed) refused at restore", why == "checkpoint field/charge seal mismatch", {{"refusal", why}});
    }
    {   // unmodified checkpoints restore and continue
        const auto dir = out / "resident-ck-ok"; std::filesystem::remove_all(dir); res.save(dir);
        std::string why;
        try { auto r2 = hc::ResidentOwner::restore(model, dev, o, none, dir, ro); hc::need(r2->field_sha256() == res.field_sha256(), "restored digest"); }
        catch (const std::exception& e) { why = e.what(); }
        R.record("divb-tamper", "control: the unmodified resident checkpoint restores with the same state digest", why.empty(), {{"error", why}});
    }
    {   // bootstrap with a nonzero ledger
        auto tq = q0; tq.divb_gross_T_m.assign(g.cells.size(), 0.); tq.divb_gross_T_m[cell] = 1e-30;
        std::string a, b;
        try { hc::Owner x(cpu, seed, tq, o, none); } catch (const std::exception& e) { a = e.what(); }
        try { hc::ResidentOwner y(model, dev, seed, tq, o, none, ro); } catch (const std::exception& e) { b = e.what(); }
        R.record("divb-tamper", "bootstrap with a nonzero divB gross ledger refused (host and resident)",
                 a == "bootstrap divB gross ledger must be zero" && b == "bootstrap divB gross ledger must be zero", {{"host", a}, {"resident", b}});
    }
    std::ofstream(out / "divb.json") << J{{"cases", R.cases}, {"failed", R.failed}}.dump(1) << '\n';
    return R.failed ? 1 : 0;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "parity" && argc >= 6) return parity(argv[2], argv[3], argv[4], argv[5], argc > 6 ? std::stoi(argv[6]) : 45);
        if (mode == "tamper" && argc == 6) return tamper(argv[2], argv[3], argv[4], argv[5]);
        if (mode == "bench" && argc == 10) return bench(argv[2], argv[3], argv[4], argv[5], std::stoi(argv[6]), std::stoul(argv[7]), std::stol(argv[8]), argv[9]);
        if (mode == "blocks" && argc == 6) return blocks(argv[2], argv[3], argv[4], argv[5]);
        if (mode == "hostloop" && argc == 8) return hostloop(argv[2], argv[3], argv[4], argv[5], std::stoul(argv[6]), std::stoi(argv[7]));
        if (mode == "divb" && argc == 7) return divb(argv[2], argv[3], argv[4], argv[5], std::stol(argv[6]));
        std::cerr << "usage: see the header of test_v34.cpp\n";
        return 2;
    } catch (const std::exception& e) { std::cerr << "fatal: " << e.what() << '\n'; return 2; }
}
