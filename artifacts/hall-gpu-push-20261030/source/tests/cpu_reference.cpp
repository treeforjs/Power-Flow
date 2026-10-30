// cpu_reference.cpp -- CPU side of the push parity harness.
//
// For every fixture it writes, under <out>/<fixture>/:
//   input.bin     identical inputs for the GPU program (mesh tables, policy, frame, dt,
//                 midpoint E/B, ID-ordered carriers, binding metadata)
//   sealed.bin    the sealed reference: hall_driver::prepare_canonical_motion and
//                 hall_corrective::charge_of_endpoints at Hall commit 8dddec0, called
//                 unchanged (record arrays, or the refusal text). Not produced above the
//                 sealed 100000-carrier cap, where the sealed function refuses.
//   replica.bin   the reference loop bodies of prepare_canonical_motion /
//                 charge_of_particles transcribed statement by statement and calling
//                 the same sealed functions (axial::push, axial::Path, point_deposit,
//                 first_kinetic_hit, prepare_neutral_stream, deposit_charge, ...),
//                 without the carrier cap and with per-carrier diagnostics. Checked
//                 against sealed.bin wherever the sealed function runs.
//   hostport.bin  the GPU kernel bodies (push/push_bodies.hpp) compiled for the host
//                 and run sequentially: separates port logic from device arithmetic.
//   cpu.json      timings, seals, refusal texts.
// Usage: cpu_reference run <C19 csv> <out dir> [fixture ...]
//        cpu_reference seal <fixture dir> <encoding.bin>     (record seal from any encoding)
//        cpu_reference selftest <out dir>                     (x87 emulation vs long double)
#include "root_api.hpp"
#include "../hall-finite-foil-em-20261029/axial_particles.hpp"

#include "push/hall_adapter.hpp"
#include "push/push_bodies.hpp"
#include "push/push_fold.hpp"
#include "tests/encoding.hpp"

#include <mpi.h>
#include <openssl/evp.h>

#include <chrono>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <optional>
#include <random>
#include <set>

namespace h = hall_corrective;
namespace n = raikou::mfem_backend;
namespace hd = hall_driver;
namespace hp = hall_gpu::push;
namespace ha = hall_gpu::push::hall_adapter;
using axial::Vec;
using h::Json;

namespace {
using Clock = std::chrono::steady_clock;
double ms_since(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }

constexpr double t0 = 7.4978805549199997e-10, dt = 1e-14;
const axial::V3 audit_point{.00021, -.00017, .001};
constexpr double m_h2o = 2.99150711e-26;

std::string hex_sha256(const std::string& bytes) {
    unsigned char out[EVP_MAX_MD_SIZE];
    unsigned len = 0;
    EVP_Digest(bytes.data(), bytes.size(), out, &len, EVP_sha256(), nullptr);
    static const char* x = "0123456789abcdef";
    std::string s;
    for (unsigned i = 0; i < len; ++i) { s.push_back(x[out[i] >> 4]); s.push_back(x[out[i] & 15]); }
    return s;
}
// Same canonical binary as the reference Encoding (native_arena_motion.cpp).
struct Enc {
    std::string bytes;
    void integer(std::uint64_t v) { for (int i = 0; i < 8; ++i) bytes.push_back(char((v >> (8 * i)) & 255)); }
    void scalar(double v) { integer(std::bit_cast<std::uint64_t>(v)); }
    void text(const std::string& s) { integer(s.size()); bytes.append(s); }
    std::string finish() const { return hex_sha256(bytes); }
};
std::string array_hash(const std::vector<double>& v) {
    Enc e;
    e.text("motion-array-binary64-v1");
    e.integer(v.size());
    for (double x : v) e.scalar(x);
    return e.finish();
}
// canonical_census_sha256 without the sealed cap (equal to the sealed helper below the cap).
std::string census_sha256(const std::vector<n::ArenaRateEntry>& rows_in) {
    std::vector<const n::ArenaRateEntry*> rows;
    for (const auto& r : rows_in) rows.push_back(&r);
    std::sort(rows.begin(), rows.end(), [](auto a, auto b) { return a->native.id < b->native.id; });
    Enc e;
    e.text("canonical-motion-census-v1");
    e.integer(rows.size());
    for (const auto* r : rows) {
        e.text(r->packet_id);
        e.text(r->state);
        const auto& p = r->native;
        e.integer(p.id);
        for (double v : p.position_m) e.scalar(v);
        for (double v : p.proper_velocity_m_s) e.scalar(v);
        e.scalar(p.charge_C);
        e.scalar(p.mass_kg);
        e.scalar(p.weight);
    }
    return e.finish();
}

// ------------------------------------------------------------------ geometry and port fixture helpers
n::GapPicOptions geometry(const axial::Mesh& g) {
    n::GapPicOptions o;
    o.nx = o.ny = g.nx;
    o.width_m = o.gap_m = g.side;
    o.depth_m = g.depth;
    o.finite_foils = o.finite_axial_extent = true;
    o.short_z_m = -g.depth / 2;
    o.feed_z_m = g.depth / 2;
    o.lower_foil = {g.side / 2 - .0025, g.side / 2 + .0025, g.side / 2 - .0027, g.side / 2 - .0025};
    o.upper_foil = {g.side / 2 - .0025, g.side / 2 + .0025, g.side / 2 + .0025, g.side / 2 + .0027};
    return o;
}
bool conductor_node(const axial::Mesh& g, int k) { return g.node_body[k] >= 0 || g.short_node[k]; }
h::ChargeLedger dc_declaration(const axial::Mesh& g, const axial::State& s, h::FreeCharge supplied) {
    h::ChargeLedger q;
    q.particles = std::move(supplied);
    q.conductor_C.assign(g.nn, h::Words());
    q.exterior_C.assign(g.nn, h::Words());
    q.conductor_declaration = h::sha256("push parity: DC induced conductor/PEC charge declared from the supplied field");
    q.exterior_declaration = h::sha256("push parity: DC exterior termination charge declared from the supplied field");
    const auto Q = h::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) {
        if (conductor_node(g, k)) q.conductor_C[k] = Q[k];
        else if (g.outer_node[k]) q.exterior_C[k] = Q[k];
    }
    return q;
}
n::ArenaRateEntry carrier(std::uint64_t id, const n::GapPicOptions& o, const axial::V3& field_x, n::PicVector u, double charge,
                          double mass, double weight, const std::string& tag) {
    n::ArenaRateEntry row;
    row.packet_id = "audit-port/packet/" + tag;
    row.state = "audit-port/resident/" + tag;
    row.native.id = id;
    row.native.position_m = {field_x[0] + o.width_m / 2, field_x[1] + o.gap_m / 2, field_x[2]};
    row.native.proper_velocity_m_s = u;
    row.native.charge_C = charge;
    row.native.mass_kg = mass;
    row.native.weight = weight;
    return row;
}

hd::MotionBinding make_binding(const axial::Mesh& g, const n::GapPicOptions& o, const std::string& label, const std::vector<n::ArenaRateEntry>& rows,
                               const Vec& E, const Vec& B, double begin, std::uint64_t step) {
    hd::MotionBinding b;
    b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = step;
    b.staged_arena_revision = 1;
    b.begin_s = begin;
    b.dt_s = dt;
    b.end_s = begin + dt;
    b.staged_arena_sha256 = rows.size() <= 100000 ? hd::canonical_census_sha256(rows) : hd::text_sha256("push parity staged arena: " + label);
    b.trial_field_sha256 = hd::total_field_sha256(E, B);
    b.mesh_sha256 = hd::motion_mesh_sha256(g, o);
    b.material_sha256 = hd::material_binding_sha256();
    b.species_sha256 = hd::text_sha256("push parity fixture species: electron, proton, H2O+, H2O");
    b.common_clock_sha256 = hd::text_sha256(Json{{"begin", b.begin_s}, {"dt", dt}, {"end", b.end_s}, {"step", step}, {"fixture", label}}.dump());
    b.kinetic_policy_sha256 = hd::kinetic_policy_sha256(hd::declared_kinetic_policy(o));
    return b;
}

// ------------------------------------------------------------------ record seal from arrays (record_without_seal)
n::ArenaJson binding_json(const hd::MotionBinding& b) {
    return {{"phase", "full_time_position_and_proper_velocity"}, {"accepted_step", b.accepted_step}, {"staged_arena_revision", b.staged_arena_revision},
            {"trial_stage", b.trial_stage}, {"begin_s", b.begin_s}, {"end_s", b.end_s}, {"dt_s", b.dt_s}, {"staged_arena_sha256", b.staged_arena_sha256},
            {"trial_field_sha256", b.trial_field_sha256}, {"mesh_sha256", b.mesh_sha256}, {"material_sha256", b.material_sha256},
            {"species_sha256", b.species_sha256}, {"common_clock_sha256", b.common_clock_sha256}, {"kinetic_policy_sha256", b.kinetic_policy_sha256}};
}
std::string record_seal(const hd::MotionBinding& b, const n::GapPicOptions& o, const std::string& census, const hgp_enc::RecordArrays& r) {
    const auto k = hd::declared_kinetic_policy(o);
    const n::PicVector frame{-o.width_m / 2, -o.gap_m / 2, 0};
    n::ArenaJson endpoints = n::ArenaJson::array(), paths = n::ArenaJson::array();
    const std::size_t N = r.id.size();
    auto vec = [](const std::vector<double>& v, std::size_t i) { return n::PicVector{v[3 * i], v[3 * i + 1], v[3 * i + 2]}; };
    for (std::size_t i = 0; i < N; ++i)
        endpoints.push_back({{"native_id", r.id[i]}, {"position_m", vec(r.end_pos, i)}, {"proper_velocity_m_s", vec(r.end_u, i)},
                             {"charge_C", r.charge[i]}, {"mass_kg", r.mass[i]}, {"weight", r.weight[i]}});
    for (std::size_t i = 0; i < N; ++i)
        paths.push_back({{"native_id", r.id[i]}, {"field_start_m", vec(r.start, i)}, {"field_displacement_m", vec(r.disp, i)},
                         {"flight_s", r.flight[i]}, {"motion_owner", "ordinary_total_field"}});
    const auto& B = r.budget;
    n::ArenaJson record = {{"schema", hd::motion_api}, {"binding", binding_json(b)}, {"census_sha256", census},
        {"kinetic_policy", {{"schema", hd::kinetic_boundary_api}, {"width_m", k.width_m}, {"gap_m", k.gap_m},
            {"short_z_m", k.short_z_m}, {"feed_z_m", k.feed_z_m},
            {"lower_foil", {k.lower.left_m, k.lower.right_m, k.lower.bottom_m, k.lower.top_m}},
            {"upper_foil", {k.upper.left_m, k.upper.right_m, k.upper.bottom_m, k.upper.top_m}},
            {"short_bridge", "lower foil x span and inner y gap only"}, {"xy_policy", "absorbing escape"},
            {"feed_policy", "axial escape"}, {"tie_order", "x y lower upper axial; later exact tie wins"}}},
        {"native_to_field_m", frame}, {"endpoints", endpoints}, {"paths", paths},
        {"current_A_m_sha256", array_hash(r.current)}, {"initial_charge_C_sha256", array_hash(r.initial_charge)},
        {"endpoint_charge_C_sha256", array_hash(r.endpoint_charge)},
        {"budget", {{"kinetic_change_J", B[0]}, {"field_work_J", B[1]}, {"work_defect_J", B[2]}, {"continuity_C", B[3]},
                    {"plasma_dt", B[4]}, {"maximum_gyro", B[5]}, {"maximum_displacement_cells", B[6]}}},
        {"coverage", "all_native_carriers_no_sources_no_boundaries"}, {"neutral_motion_owner", "same_full_time_q0_drift"},
        {"native_Poisson_called", false}, {"native_DKD_called", false}};
    Enc seal;
    seal.text("motion-receipt-record-v1");
    seal.text(record.dump());
    return seal.finish();
}

// ------------------------------------------------------------------ status mapping
int code_of(const std::string& message) {
    static const int codes[] = {1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 40, 41, 42, 43, 44};
    for (int c : codes) if (message == hp::code_message(c)) return c;
    return 98;
}

struct Diag {
    std::vector<std::int32_t> code1, code2, code3;
    std::vector<std::vector<double>> real;
    std::vector<std::vector<std::int32_t>> integer;
    void init(std::size_t N) {
        code1.assign(N, 0); code2.assign(N, hgp::skipped); code3.assign(N, hgp::skipped);
        real.assign(hgp::diag_real_fields, std::vector<double>(N, 0.0));
        integer.assign(hgp::diag_int_fields, std::vector<std::int32_t>(N, 0));
        real[hgp::diag_hit_fraction].assign(N, 1.0);
        integer[hgp::diag_hit_surface].assign(N, -1);
        integer[hgp::diag_kinetic_sink].assign(N, -1);
        integer[hgp::diag_kinetic_axis].assign(N, -1);
    }
    void write(hgp_enc::Writer& w) const {
        w.i32("diag.code1", code1); w.i32("diag.code2", code2); w.i32("diag.code3", code3);
        for (int k = 0; k < hgp::diag_real_fields; ++k) w.f64(std::string("diag.real.") + hgp_enc::diag_real_names[k], real[k]);
        for (int k = 0; k < hgp::diag_int_fields; ++k) w.i32(std::string("diag.int.") + hgp_enc::diag_int_names[k], integer[k]);
    }
};

struct Result {
    hgp_enc::RecordArrays rec;
    Diag diag;
    Json timing = Json::object();
};

void fill_identity(hgp_enc::RecordArrays& r, const ha::OrderedCensus& c) {
    r.id = c.id;
    r.charge = c.charge;
    r.mass = c.mass;
    r.weight = c.weight;
    r.flight.assign(c.id.size(), dt);
}

// ------------------------------------------------------------------ sealed reference
struct Sealed {
    bool ran = false;
    hgp_enc::RecordArrays rec;
    std::string seal, census;
    Json timing = Json::object();
};
Sealed sealed_run(const axial::Mesh& g, const n::GapPicOptions& o, const std::vector<n::ArenaRateEntry>& rows, const hd::MotionBinding& b,
                  const Vec& E, const Vec& B, const hd::MotionReceipt* external) {
    Sealed s;
    std::optional<hd::MotionReceipt> own;
    const hd::MotionReceipt* r = external;
    if (!r) {
        const auto t = Clock::now();
        try {
            own.emplace(hd::prepare_canonical_motion(g, o, rows, b, {}, E, B));
            r = &*own;
        } catch (const std::exception& e) {
            s.rec.refusal = e.what();
        }
        s.timing["prepare_canonical_motion_ms"] = ms_since(t);
    }
    s.ran = rows.size() <= 100000;
    if (!r) return s;
    auto& a = s.rec;
    for (const auto& p : r->endpoints()) {
        a.id.push_back(p.id);
        for (double v : p.position_m) a.end_pos.push_back(v);
        for (double v : p.proper_velocity_m_s) a.end_u.push_back(v);
        a.charge.push_back(p.charge_C); a.mass.push_back(p.mass_kg); a.weight.push_back(p.weight);
    }
    for (const auto& p : r->paths()) {
        for (double v : p.field_start_m) a.start.push_back(v);
        for (double v : p.field_displacement_m) a.disp.push_back(v);
        a.flight.push_back(p.flight_s);
    }
    a.current = r->current_A_m(); a.initial_charge = r->initial_charge_C(); a.endpoint_charge = r->endpoint_charge_C();
    const auto& bu = r->budget();
    a.budget = {bu.kinetic_change_J, bu.field_work_J, bu.work_defect_J, bu.continuity_C, bu.plasma_dt, bu.maximum_gyro, bu.maximum_displacement_cells};
    s.seal = r->seal_sha256();
    s.census = r->census_sha256();
    const auto t = Clock::now();
    try {
        const auto fc = h::charge_of_endpoints(g, o, r->endpoints());
        for (const auto& w : fc.net_C) { a.net_hi.push_back(w.hi); a.net_lo.push_back(w.lo); }
        a.gross = fc.gross_C;
    } catch (const std::exception& e) {
        a.charge_refusal = e.what();
    }
    s.timing["charge_of_endpoints_ms"] = ms_since(t);
    return s;
}

// ------------------------------------------------------------------ replica (sealed functions, transcribed loop bodies)
Result replica(const axial::Mesh& g, const n::GapPicOptions& o, const ha::OrderedCensus& c, const Vec& E, const Vec& B, double begin,
               bool diagnostics) {
    Result out;
    auto& rec = out.rec;
    auto& dg = out.diag;
    const std::size_t N = c.id.size();
    dg.init(N);
    fill_identity(rec, c);
    const auto policy = hd::declared_kinetic_policy(o);
    const n::PicVector frame{-o.width_m / 2, -o.gap_m / 2, 0};
    const double end = begin + dt;
    auto particle = [&](std::size_t i) {
        n::GapParticle p;
        p.id = c.id[i];
        p.position_m = {c.x[i], c.y[i], c.z[i]};
        p.proper_velocity_m_s = {c.ux[i], c.uy[i], c.uz[i]};
        p.charge_C = c.charge[i]; p.mass_kg = c.mass[i]; p.weight = c.weight[i];
        return p;
    };
    long first1 = -1, first2 = -1, first3 = -1;
    // phase 1: census loop
    auto t = Clock::now();
    Vec omega(g.nn);
    for (std::size_t i = 0; i < N; ++i) {
        try {
            const auto p = particle(i);
            hall::require(std::isfinite(p.weight) && p.weight > 0 && std::isfinite(p.mass_kg) && p.mass_kg > 0 && std::isfinite(p.charge_C), "motion invalid carrier");
            for (const auto& v : {p.position_m, p.proper_velocity_m_s}) for (double x : v) hall::require(std::isfinite(x), "motion nonfinite kinematics");
            hd::require_kinetic_interior(policy, p.position_m);
            auto field_x = p.position_m;
            for (int d = 0; d < 3; ++d) field_x[d] += frame[d];
            double load = p.weight * p.charge_C * p.charge_C / (axial::eps * p.mass_kg);
            hall::require(std::isfinite(load) && (p.charge_C == 0 || load > 0), "motion unresolved plasma operand");
            axial::point_deposit(g, omega, field_x, load);
        } catch (const std::exception& e) {
            dg.code1[i] = code_of(e.what());
            if (first1 < 0) first1 = long(i);
            if (!diagnostics) break;
        }
    }
    double plasma_dt = 0;
    for (int k = 0; k < g.nn; ++k) plasma_dt = std::max(plasma_dt, dt * std::sqrt(omega[k] / g.VN[k]));
    out.timing["census_loop_ms"] = ms_since(t);
    // phase 2: motion per carrier
    t = Clock::now();
    std::vector<n::GapParticle> endpoint(N);
    std::vector<n::PicVector> start(N), D(N);
    std::vector<double> dke(N, 0.0), gyro(N, 0.0), dcells(N, 0.0);
    Vec scratch(g.nn);
    h::FreeCharge scratch_fc{h::WordVec(g.nn), Vec(g.nn)};
    const bool run2 = diagnostics || (first1 < 0 && std::isfinite(plasma_dt) && plasma_dt <= .1);
    for (std::size_t i = 0; run2 && i < N; ++i) {
        if (dg.code1[i] != 0) continue;
        dg.code2[i] = 0;
        try {
            const auto p = particle(i);
            auto field_p = p;
            for (int d = 0; d < 3; ++d) field_p.position_m[d] += frame[d];
            axial::Push moved{};
            bool moved_valid = false;
            if (p.charge_C == 0) {
                auto stream = hd::prepare_neutral_stream(policy, p, begin, end, dt);
                if (stream.hit.sink != hd::KineticSink::none) {
                    dg.integer[hgp::diag_kinetic_valid][i] = 1;
                    dg.integer[hgp::diag_kinetic_sink][i] = int(stream.hit.sink);
                    dg.integer[hgp::diag_kinetic_axis][i] = stream.hit.axis;
                    dg.real[hgp::diag_kinetic_fraction][i] = stream.hit.fraction;
                    dg.real[hgp::diag_kinetic_plane][i] = stream.hit.plane_m;
                }
                hall::require(stream.hit.sink == hd::KineticSink::none, "neutral first hit requires native sink owner");
                moved.particle = field_p;
                moved.hit = {1, -1};
                moved.flight_s = dt;
                for (int d = 0; d < 3; ++d) moved.D[d] = stream.endpoint.position_m[d] - p.position_m[d];
                moved.displacement_cells = std::sqrt(hall::dot(moved.D, moved.D)) / std::min({g.dx, g.dy, g.dz});
                moved_valid = true;
            } else {
                moved = axial::push(g, field_p, E, B, dt);
                moved_valid = true;
            }
            if (moved_valid) {
                for (int d = 0; d < 3; ++d) {
                    dg.real[hgp::diag_moved_dx + d][i] = moved.D[d];
                    dg.real[hgp::diag_moved_ux + d][i] = moved.particle.proper_velocity_m_s[d];
                }
                dg.real[hgp::diag_hit_fraction][i] = moved.hit.fraction;
                dg.real[hgp::diag_flight_s][i] = moved.flight_s;
                dg.real[hgp::diag_gyro][i] = moved.gyro;
                dg.real[hgp::diag_displacement_cells][i] = moved.displacement_cells;
                dg.integer[hgp::diag_moved_valid][i] = 1;
                dg.integer[hgp::diag_hit_surface][i] = moved.hit.surface;
                if (diagnostics) {
                    try {
                        const auto kh = hd::first_kinetic_hit(policy, p.position_m, moved.D);
                        dg.integer[hgp::diag_kinetic_valid][i] = 1;
                        dg.integer[hgp::diag_kinetic_sink][i] = int(kh.sink);
                        dg.integer[hgp::diag_kinetic_axis][i] = kh.axis;
                        dg.real[hgp::diag_kinetic_fraction][i] = kh.fraction;
                        dg.real[hgp::diag_kinetic_plane][i] = kh.plane_m;
                    } catch (const std::exception&) {
                    }
                }
            }
            if (p.charge_C == 0) hall::require(moved.displacement_cells <= .2, "neutral full-interval displacement guard");
            hall::require(moved.hit.surface < 0 && moved.flight_s == dt, "motion needs boundary owner before collection");
            hall::require(hd::first_kinetic_hit(policy, p.position_m, moved.D).sink == hd::KineticSink::none, "charged first hit requires native sink owner");
            auto ep = p;
            ep.proper_velocity_m_s = moved.particle.proper_velocity_m_s;
            for (int d = 0; d < 3; ++d) ep.position_m[d] = p.position_m[d] + moved.D[d];
            hd::require_kinetic_interior(policy, ep.position_m);
            n::PicVector Dp{};
            for (int d = 0; d < 3; ++d) Dp[d] = (ep.position_m[d] + frame[d]) - field_p.position_m[d];
            hall::require(axial::intersection(g, field_p.position_m, Dp).surface < 0, "rounded canonical endpoint crosses field boundary");
            axial::Path path(g, field_p.position_m, Dp);
            double charge = p.weight * p.charge_C;
            hall::require(std::isfinite(charge) && (p.charge_C == 0 || charge != 0), "motion unresolved macrocharge");
            auto field_end = ep.position_m;
            for (int d = 0; d < 3; ++d) field_end[d] += frame[d];
            axial::point_deposit(g, scratch, field_end, charge);   // the reference's endpoint range check, same call
            if (diagnostics) {   // per-carrier checks of charge_of_particles (phase 5 repeats them when admitted)
                dg.code3[i] = 0;
                try {
                    h::need(std::isfinite(ep.weight) && ep.weight > 0 && std::isfinite(ep.charge_C), "invalid canonical charge operand");
                    hd::require_kinetic_interior(policy, ep.position_m);
                    auto xc = ep.position_m;
                    xc[0] -= o.width_m / 2;
                    xc[1] -= o.gap_m / 2;
                    const auto qw = h::Words(ep.weight) * ep.charge_C;
                    h::need(h::canonical(qw) && (ep.charge_C == 0 || qw.hi != 0 || qw.lo != 0), "unresolved canonical charge");
                    h::deposit_charge(g, scratch_fc, xc, qw);
                } catch (const std::exception& e) {
                    dg.code3[i] = code_of(e.what());
                }
            }
            endpoint[i] = ep;
            start[i] = field_p.position_m;
            D[i] = Dp;
            dke[i] = n::pic_kinetic_energy_J(ep) - n::pic_kinetic_energy_J(p);
            gyro[i] = moved.gyro;
            dcells[i] = moved.displacement_cells;
        } catch (const std::exception& e) {
            dg.code2[i] = code_of(e.what());
            if (first2 < 0) first2 = long(i);
            if (!diagnostics) break;
        }
    }
    out.timing["motion_loop_ms"] = ms_since(t);
    // refusal order
    if (first1 >= 0) { rec.refusal = hp::code_message(dg.code1[std::size_t(first1)]); rec.refusal_rank = int(first1); }
    else if (!(std::isfinite(plasma_dt) && plasma_dt <= .1)) rec.refusal = "motion full-census plasma guard";
    else if (first2 >= 0) { rec.refusal = hp::code_message(dg.code2[std::size_t(first2)]); rec.refusal_rank = int(first2); }
    if (!rec.refusal.empty()) return out;
    // phase 3: deposits in carrier order
    t = Clock::now();
    Vec current(g.ne), initial(g.nn), endc(g.nn), gross(g.nn);
    long double kinetic = 0;
    double max_gyro = 0, max_dcells = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const double charge = c.weight[i] * c.charge[i];
        axial::Path path(g, start[i], D[i]);
        path.deposit(current, charge, D[i], dt);
        axial::point_deposit(g, initial, start[i], charge);
        auto field_end = endpoint[i].position_m;
        for (int d = 0; d < 3; ++d) field_end[d] += frame[d];
        axial::point_deposit(g, endc, field_end, charge);
        axial::point_deposit(g, gross, start[i], std::abs(charge));
        kinetic += static_cast<long double>(c.weight[i]) * dke[i];
        max_gyro = std::max(max_gyro, gyro[i]);
        max_dcells = std::max(max_dcells, dcells[i]);
    }
    out.timing["deposit_loop_ms"] = ms_since(t);
    // phase 4: budget
    t = Clock::now();
    hd::MotionBudget b;
    b.plasma_dt = plasma_dt;
    b.maximum_gyro = max_gyro;
    b.maximum_displacement_cells = max_dcells;
    b.kinetic_change_J = double(kinetic);
    b.field_work_J = dt * axial::dot(E, current);
    b.work_defect_J = b.kinetic_change_J - b.field_work_J;
    if (!(std::isfinite(b.work_defect_J) && std::abs(b.work_defect_J) <= 1e-20 + 1e-10 * (std::abs(b.kinetic_change_J) + std::abs(b.field_work_J)))) {
        rec.refusal = "canonical motion work mismatch";
        return out;
    }
    auto div = g.gradient_transpose(current);
    for (int k = 0; k < g.nn; ++k) b.continuity_C = std::max(b.continuity_C, std::abs(endc[k] - initial[k] - dt * div[k]));
    if (!(axial::finite(current) && std::isfinite(b.continuity_C) && b.continuity_C <= 1e-12 * std::max(axial::maximum(gross), 1e-30))) {
        rec.refusal = "canonical motion charge continuity mismatch";
        return out;
    }
    out.timing["budget_ms"] = ms_since(t);
    for (std::size_t i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) {
            rec.end_pos.push_back(endpoint[i].position_m[d]);
            rec.end_u.push_back(endpoint[i].proper_velocity_m_s[d]);
        }
        for (int d = 0; d < 3; ++d) rec.start.push_back(start[i][d]);
        for (int d = 0; d < 3; ++d) rec.disp.push_back(D[i][d]);
    }
    rec.current = current; rec.initial_charge = initial; rec.endpoint_charge = endc;
    rec.budget = {b.kinetic_change_J, b.field_work_J, b.work_defect_J, b.continuity_C, b.plasma_dt, b.maximum_gyro, b.maximum_displacement_cells};
    // phase 5: charge_of_particles over the endpoints (transcribed; uncapped)
    t = Clock::now();
    h::FreeCharge fc{h::WordVec(g.nn), Vec(g.nn)};
    for (std::size_t i = 0; i < N; ++i) {
        dg.code3[i] = 0;
        try {
            const auto* p = &endpoint[i];
            h::need(std::isfinite(p->weight) && p->weight > 0 && std::isfinite(p->charge_C), "invalid canonical charge operand");
            hd::require_kinetic_interior(policy, p->position_m);
            auto x = p->position_m;
            x[0] -= o.width_m / 2;
            x[1] -= o.gap_m / 2;
            const auto q = h::Words(p->weight) * p->charge_C;
            h::need(h::canonical(q) && (p->charge_C == 0 || q.hi != 0 || q.lo != 0), "unresolved canonical charge");
            h::deposit_charge(g, fc, x, q);
        } catch (const std::exception& e) {
            dg.code3[i] = code_of(e.what());
            if (rec.charge_refusal.empty()) { rec.charge_refusal = e.what(); rec.charge_refusal_rank = int(i); }
            if (!diagnostics) break;
        }
    }
    if (rec.charge_refusal.empty()) {
        try {
            h::free_charge_guard(g, fc);
            for (const auto& w : fc.net_C) { rec.net_hi.push_back(w.hi); rec.net_lo.push_back(w.lo); }
            rec.gross = fc.gross_C;
        } catch (const std::exception& e) {
            rec.charge_refusal = e.what();
        }
    }
    out.timing["free_charge_ms"] = ms_since(t);
    return out;
}

// ------------------------------------------------------------------ host port: device bodies run sequentially on the CPU
Result host_port(const hp::MeshTables& mt, const hp::Geometry& geo, const ha::OrderedCensus& c, const Vec& E, const Vec& B, bool diagnostics) {
    using hgp_u32 = std::uint32_t;
    Result out;
    const std::size_t N = c.id.size(), nn = std::size_t(mt.mesh.nn), ne = std::size_t(mt.mesh.ne);
    std::vector<int> codes(3 * N);
    std::vector<hgp_u32> fails(3 * N), counts(N), offsets(N), skey(8 * N), spay(8 * N), ekey(8 * N), epay(8 * N);
    std::vector<double> outs(15 * N), sval(3 * 8 * N), evals(3 * 8 * N), dreal(std::size_t(hgp::diag_real_fields) * N);
    std::vector<int> dint(std::size_t(hgp::diag_int_fields) * N);
    hgp::CarrierArgs a{};
    a.g = mt.mesh; a.p = geo.policy;
    for (int d = 0; d < 3; ++d) a.off[d] = geo.native_to_field[d];
    a.dt = dt; a.n = hgp_u32(N); a.node_sentinel = hgp_u32(nn); a.diagnostics = 1;
    a.x = c.x.data(); a.y = c.y.data(); a.z = c.z.data(); a.ux = c.ux.data(); a.uy = c.uy.data(); a.uz = c.uz.data();
    a.q = c.charge.data(); a.m = c.mass.data(); a.w = c.weight.data(); a.E = E.data(); a.B = B.data();
    a.code1 = codes.data(); a.code2 = codes.data() + N; a.code3 = codes.data() + 2 * N;
    a.fail1 = fails.data(); a.fail2 = fails.data() + N; a.fail3 = fails.data() + 2 * N;
    a.endpoint = outs.data(); a.start = outs.data() + 6 * N; a.disp = outs.data() + 9 * N; a.dke = outs.data() + 12 * N;
    a.gyro = outs.data() + 13 * N; a.dcells = outs.data() + 14 * N; a.edge_count = counts.data();
    a.start_key = skey.data(); a.start_pay = spay.data(); a.start_omega = sval.data(); a.start_init = sval.data() + 8 * N; a.start_gross = sval.data() + 16 * N;
    a.end_key = ekey.data(); a.end_pay = epay.data(); a.end_charge = evals.data(); a.end_hi = evals.data() + 8 * N; a.end_lo = evals.data() + 16 * N;
    a.diag_real = dreal.data(); a.diag_int = dint.data();
    for (std::size_t i = 0; i < N; ++i) hgp::carrier_body(a, hgp_u32(i));
    std::uint64_t T = 0;
    for (std::size_t i = 0; i < N; ++i) { offsets[i] = hgp_u32(T); T += counts[i]; }
    std::vector<hgp_u32> gkey(T), gpay(T);
    std::vector<double> gval(T);
    hgp::EmitArgs em{};
    em.g = mt.mesh; em.dt = dt; em.n = hgp_u32(N); em.q = a.q; em.w = a.w; em.start = a.start; em.disp = a.disp;
    em.edge_count = counts.data(); em.edge_offset = offsets.data(); em.key = gkey.data(); em.pay = gpay.data(); em.value = gval.data();
    for (std::size_t i = 0; i < N; ++i) hgp::emit_body(em, hgp_u32(i));
    auto sort_bounds = [](const std::vector<hgp_u32>& key, const std::vector<hgp_u32>& pay, std::size_t keys, std::vector<hgp_u32>& spay2,
                          std::vector<hgp_u32>& begin, std::vector<hgp_u32>& end) {
        std::vector<std::size_t> idx(key.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&](std::size_t x, std::size_t y) { return key[x] < key[y]; });
        spay2.resize(key.size());
        begin.assign(keys, 0); end.assign(keys, 0);
        for (std::size_t r = 0; r < idx.size(); ++r) {
            spay2[r] = pay[idx[r]];
            const hgp_u32 k = key[idx[r]];
            if (r == 0 || key[idx[r - 1]] != k) begin[k] = hgp_u32(r);
            if (r + 1 == idx.size() || key[idx[r + 1]] != k) end[k] = hgp_u32(r + 1);
        }
    };
    std::vector<hgp_u32> sorted, bb, be;
    std::vector<double> current(ne), node(8 * nn);
    sort_bounds(gkey, gpay, ne, sorted, bb, be);
    hgp::NodeReduceArgs re{hgp_u32(ne), 0, sorted.data(), bb.data(), be.data(), gval.data(), nullptr, nullptr, current.data(), nullptr, nullptr};
    for (std::size_t e = 0; e < ne; ++e) hgp::reduce_edges_body(re, hgp_u32(e));
    sort_bounds(skey, spay, nn + 1, sorted, bb, be);
    hgp::NodeReduceArgs rs{hgp_u32(nn), 0, sorted.data(), bb.data(), be.data(), sval.data(), sval.data() + 8 * N, sval.data() + 16 * N,
                           node.data(), node.data() + nn, node.data() + 2 * nn};
    for (std::size_t k = 0; k < nn; ++k) hgp::reduce_start_body(rs, hgp_u32(k));
    sort_bounds(ekey, epay, nn + 1, sorted, bb, be);
    hgp::EndNodeArgs rn{hgp_u32(nn), 0, sorted.data(), bb.data(), be.data(), evals.data(), evals.data() + 8 * N, evals.data() + 16 * N,
                        node.data() + 3 * nn, node.data() + 4 * nn, node.data() + 5 * nn, node.data() + 6 * nn};
    for (std::size_t k = 0; k < nn; ++k) hgp::reduce_end_body(rn, hgp_u32(k));
    hp::Budget b;
    for (std::size_t k = 0; k < nn; ++k) b.plasma_dt = std::max(b.plasma_dt, dt * std::sqrt(node[k] / mt.VN[k]));
    for (std::size_t i = 0; i < N; ++i) { b.maximum_gyro = std::max(b.maximum_gyro, outs[13 * N + i]); b.maximum_displacement_cells = std::max(b.maximum_displacement_cells, outs[14 * N + i]); }
    hp::FoldInput f;
    f.n = N; f.ids = c.id; f.weights = c.weight; f.dke = {outs.data() + 12 * N, N}; f.E = E; f.current = current;
    f.initial = {node.data() + nn, nn}; f.gross = {node.data() + 2 * nn, nn}; f.endpoint = {node.data() + 3 * nn, nn};
    f.net_hi = {node.data() + 4 * nn, nn}; f.net_lo = {node.data() + 5 * nn, nn}; f.free_gross = {node.data() + 6 * nn, nn};
    f.mt = &mt; f.dt = dt;
    for (int k = 0; k < 3; ++k) {
        std::uint32_t r = 0xffffffffu;
        for (std::size_t i = 0; i < N; ++i) r = std::min(r, fails[std::size_t(k) * N + i]);
        f.status[2 * k] = r;
        f.status[2 * k + 1] = r != 0xffffffffu ? std::uint32_t(codes[std::size_t(k) * N + r]) : 0u;
    }
    hp::Refusal refusal, charge_refusal;
    hp::fold_record(f, refusal, charge_refusal, b);
    auto& rec = out.rec;
    fill_identity(rec, c);
    rec.refusal = refusal.refused ? refusal.message : "";
    rec.refusal_rank = int(refusal.rank);
    rec.charge_refusal = charge_refusal.refused ? charge_refusal.message : "";
    rec.charge_refusal_rank = int(charge_refusal.rank);
    for (std::size_t i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) rec.end_pos.push_back(outs[d * N + i]);
        for (int d = 0; d < 3; ++d) rec.end_u.push_back(outs[(3 + d) * N + i]);
        for (int d = 0; d < 3; ++d) rec.start.push_back(outs[(6 + d) * N + i]);
        for (int d = 0; d < 3; ++d) rec.disp.push_back(outs[(9 + d) * N + i]);
    }
    rec.current = current;
    rec.initial_charge.assign(node.begin() + std::ptrdiff_t(nn), node.begin() + std::ptrdiff_t(2 * nn));
    rec.endpoint_charge.assign(node.begin() + std::ptrdiff_t(3 * nn), node.begin() + std::ptrdiff_t(4 * nn));
    rec.budget = {b.kinetic_change_J, b.field_work_J, b.work_defect_J, b.continuity_C, b.plasma_dt, b.maximum_gyro, b.maximum_displacement_cells};
    rec.net_hi.assign(node.begin() + std::ptrdiff_t(4 * nn), node.begin() + std::ptrdiff_t(5 * nn));
    rec.net_lo.assign(node.begin() + std::ptrdiff_t(5 * nn), node.begin() + std::ptrdiff_t(6 * nn));
    rec.gross.assign(node.begin() + std::ptrdiff_t(6 * nn), node.begin() + std::ptrdiff_t(7 * nn));
    if (diagnostics) {
        auto& dg = out.diag;
        dg.code1.assign(codes.begin(), codes.begin() + std::ptrdiff_t(N));
        dg.code2.assign(codes.begin() + std::ptrdiff_t(N), codes.begin() + std::ptrdiff_t(2 * N));
        dg.code3.assign(codes.begin() + std::ptrdiff_t(2 * N), codes.end());
        dg.real.assign(hgp::diag_real_fields, {});
        dg.integer.assign(hgp::diag_int_fields, {});
        for (int k = 0; k < hgp::diag_real_fields; ++k) dg.real[k].assign(dreal.begin() + std::ptrdiff_t(k * N), dreal.begin() + std::ptrdiff_t((k + 1) * N));
        for (int k = 0; k < hgp::diag_int_fields; ++k) dg.integer[k].assign(dint.begin() + std::ptrdiff_t(k * N), dint.begin() + std::ptrdiff_t((k + 1) * N));
    }
    return out;
}

// ------------------------------------------------------------------ synthetic populations
struct Rng {
    std::uint64_t s;
    std::uint64_t next() {
        std::uint64_t z = (s += 0x9e3779b97f4a7c15ull);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
        return z ^ (z >> 31);
    }
    double uniform() { return double(next() >> 11) * 0x1.0p-53; }
    double range(double a, double b) { return a + (b - a) * uniform(); }
    double log_range(double a, double b) { return std::exp(range(std::log(a), std::log(b))); }
};
struct Species { const char* name; double charge, mass, emin_eV, emax_eV; };
const Species species[4] = {{"electron", -n::pic_e, n::pic_me, .1, 1e4},
                            {"proton", n::pic_e, hall::proton_mass, .01, 1e4},
                            {"H2O+", n::pic_e, m_h2o - n::pic_me, .01, 1e4},
                            {"H2O", 0.0, m_h2o, .01, 1.0}};

struct Generator {
    const axial::Mesh& g;
    const n::GapPicOptions& o;
    hd::KineticPolicy pol;
    Rng rng;
    Generator(const axial::Mesh& mesh, const n::GapPicOptions& opt, std::uint64_t seed) : g(mesh), o(opt), pol(hd::declared_kinetic_policy(opt)), rng{seed} {}
    int pick_species() { const double r = rng.uniform(); return r < .4 ? 0 : r < .6 ? 1 : r < .8 ? 2 : 3; }
    n::PicVector velocity(const Species& s, double& speed_m_s) {
        const double K = rng.log_range(s.emin_eV, s.emax_eV) * n::pic_e, mc2 = s.mass * n::pic_c * n::pic_c;
        const double u = std::sqrt(K * (K + 2 * mc2)) / (s.mass * n::pic_c);
        const double cz = rng.range(-1, 1), ph = rng.range(0, 2 * M_PI), sz = std::sqrt(std::max(0.0, 1 - cz * cz));
        speed_m_s = u / std::sqrt(1 + u * u / (n::pic_c * n::pic_c));
        return {u * sz * std::cos(ph), u * sz * std::sin(ph), u * cz};
    }
    bool clear_of_metal(const n::PicVector& x, double margin) const {
        if (!(x[0] > margin && x[0] < pol.width_m - margin && x[1] > margin && x[1] < pol.gap_m - margin &&
              x[2] > pol.short_z_m + margin && x[2] < pol.feed_z_m - margin)) return false;
        for (const auto& f : {pol.lower, pol.upper})
            if (x[0] >= f.left_m - margin && x[0] <= f.right_m + margin && x[1] >= f.bottom_m - margin && x[1] <= f.top_m + margin) return false;
        return true;
    }
    n::PicVector interior(double margin) {
        for (;;) {
            n::PicVector x{rng.range(0, pol.width_m), rng.range(0, pol.gap_m), rng.range(pol.short_z_m, pol.feed_z_m)};
            if (clear_of_metal(x, margin)) return x;
        }
    }
    // Carrier placed so that it crosses interior cell faces (1 or 2 axes) during the step.
    void face_crosser(n::GapParticle& p, const Species& s) {
        for (;;) {
            double speed = 0;
            p.proper_velocity_m_s = velocity(s, speed);
            auto x = interior(20e-6);
            const int axes = rng.uniform() < .15 ? 2 : 1;
            const int first = rng.uniform() < .9 ? int(rng.next() % 2) : 2;
            bool ok = true;
            for (int k = 0; k < axes && ok; ++k) {
                const int a = (first + k) % 3;
                const double gam = std::sqrt(1 + hall::dot(p.proper_velocity_m_s, p.proper_velocity_m_s) / (n::pic_c * n::pic_c));
                const double step = std::abs(p.proper_velocity_m_s[a]) / gam * dt;
                if (step < 1e-15) { ok = false; break; }
                double face;
                if (a < 2) face = double(1 + rng.next() % std::uint64_t(g.nx - 1)) * (a ? g.dy : g.dx);   // native face positions
                else face = g.zmin + double(2 + rng.next() % 3) * g.dz;                                   // interior z faces
                const double sign = p.proper_velocity_m_s[a] > 0 ? 1 : -1;
                x[a] = face - sign * rng.range(.1, .9) * step;
            }
            if (!ok || !clear_of_metal(x, 5e-6)) continue;
            p.position_m = x;
            return;
        }
    }
    // Carrier placed to reach a kinetic or field boundary during the step.
    void boundary_approacher(n::GapParticle& p, const Species& s) {
        for (;;) {
            double speed = 0;
            p.proper_velocity_m_s = velocity(s, speed);
            auto& u = p.proper_velocity_m_s;
            const double gam = std::sqrt(1 + hall::dot(u, u) / (n::pic_c * n::pic_c));
            auto x = interior(20e-6);
            const int kind = int(rng.next() % 3);
            int a = 0;
            double plane = 0;
            if (kind == 0) {   // foil face
                const auto& f = rng.uniform() < .5 ? pol.lower : pol.upper;
                const int face = int(rng.next() % 4);
                a = face < 2 ? 0 : 1;
                plane = face == 0 ? f.left_m : face == 1 ? f.right_m : face == 2 ? f.bottom_m : f.top_m;
                const bool outward_low = face == 0 || face == 2;   // approach from below the low face
                u[a] = std::abs(u[a]) * (outward_low ? 1 : -1);
                const int b = 1 - a;
                x[b] = b == 0 ? rng.range(f.left_m + 1e-6, f.right_m - 1e-6) : rng.range(f.bottom_m + 1e-6, f.top_m - 1e-6);
            } else if (kind == 1) {   // outer xy wall
                a = int(rng.next() % 2);
                const bool high = rng.uniform() < .5;
                plane = high ? (a ? pol.gap_m : pol.width_m) : 0;
                u[a] = std::abs(u[a]) * (high ? 1 : -1);
            } else {   // short or feed plane; part of the short plane lies under the bridge
                a = 2;
                const bool feed = rng.uniform() < .5;
                plane = feed ? pol.feed_z_m : pol.short_z_m;
                u[2] = std::abs(u[2]) * (feed ? 1 : -1);
                if (!feed && rng.uniform() < .5) { x[0] = rng.range(pol.lower.left_m, pol.lower.right_m); x[1] = rng.range(pol.lower.top_m + 1e-6, pol.upper.bottom_m - 1e-6); }
            }
            const double step = std::abs(u[a]) / gam * dt;
            if (step < 1e-15) continue;
            x[a] = plane - (u[a] > 0 ? 1 : -1) * rng.range(.1, .9) * step;
            try { hd::require_kinetic_interior(pol, x); } catch (const std::exception&) { continue; }
            p.position_m = x;
            return;
        }
    }
    std::vector<n::ArenaRateEntry> population(std::size_t N, double crossers, double approachers) {
        std::vector<n::ArenaRateEntry> rows(N);
        std::uint64_t id = 0;
        for (std::size_t i = 0; i < N; ++i) {
            const auto& s = species[pick_species()];
            auto& r = rows[i];
            id += 1 + rng.next() % 3;
            r.native.id = id;
            r.native.charge_C = s.charge;
            r.native.mass_kg = s.mass;
            r.native.weight = rng.log_range(1e3, 1e8);
            r.packet_id = "push-parity/packet/" + std::to_string(id);
            r.state = std::string("push-parity/resident/") + s.name;
            const double c = rng.uniform();
            if (c < approachers) boundary_approacher(r.native, s);
            else if (c < approachers + crossers) face_crosser(r.native, s);
            else { double speed = 0; r.native.proper_velocity_m_s = velocity(s, speed); r.native.position_m = interior(5e-6); }
        }
        return rows;
    }
};

// ------------------------------------------------------------------ fixture export
void write_input(const std::string& path, const hp::MeshTables& mt, const hp::Geometry& geo, const ha::OrderedCensus& c, const Vec& E,
                 const Vec& B, bool diagnostics, const hd::MotionBinding& b, const std::string& census) {
    hgp_enc::Writer w;
    const auto& m = mt.mesh;
    const std::int32_t ints[10] = {m.nx, m.ny, m.nz, m.nex, m.ney, m.nbx, m.nby, m.nn, m.ne, m.nb};
    const double reals[6] = {m.side, m.depth, m.dx, m.dy, m.dz, m.zmin};
    w.i32("mesh.int", ints);
    w.f64("mesh.real", reals);
    w.f64("mesh.VN", mt.VN);
    w.i32("mesh.edge_a", mt.edge_a);
    w.i32("mesh.edge_b", mt.edge_b);
    w.f64("mesh.edge_length", mt.edge_length);
    const auto& p = geo.policy;
    const double pol[12] = {p.width, p.gap, p.short_z, p.feed_z, p.lower.left, p.lower.right, p.lower.bottom, p.lower.top,
                            p.upper.left, p.upper.right, p.upper.bottom, p.upper.top};
    w.f64("policy", pol);
    w.f64("frame", geo.native_to_field);
    const double t[1] = {dt};
    w.f64("dt", t);
    w.f64("field.E", E);
    w.f64("field.B", B);
    w.u64("carrier.id", c.id);
    w.f64("carrier.x", c.x); w.f64("carrier.y", c.y); w.f64("carrier.z", c.z);
    w.f64("carrier.ux", c.ux); w.f64("carrier.uy", c.uy); w.f64("carrier.uz", c.uz);
    w.f64("carrier.charge", c.charge); w.f64("carrier.mass", c.mass); w.f64("carrier.weight", c.weight);
    w.text("mode", diagnostics ? "diagnostics" : "record");
    const std::uint64_t bnum[3] = {b.accepted_step, b.staged_arena_revision, b.trial_stage};
    const double btime[3] = {b.begin_s, b.end_s, b.dt_s};
    w.u64("binding.steps", bnum);
    w.f64("binding.times", btime);
    for (auto [k, v] : std::initializer_list<std::pair<const char*, const std::string*>>{
             {"binding.staged_arena_sha256", &b.staged_arena_sha256}, {"binding.trial_field_sha256", &b.trial_field_sha256},
             {"binding.mesh_sha256", &b.mesh_sha256}, {"binding.material_sha256", &b.material_sha256}, {"binding.species_sha256", &b.species_sha256},
             {"binding.common_clock_sha256", &b.common_clock_sha256}, {"binding.kinetic_policy_sha256", &b.kinetic_policy_sha256}})
        w.text(k, *v);
    w.text("census_sha256", census);
    w.save(path);
}

void write_result(const std::string& path, const Result& r, bool diagnostics) {
    hgp_enc::Writer w;
    r.rec.write(w, true);
    if (diagnostics) r.diag.write(w);
    w.save(path);
}

struct Case {
    std::string name;
    bool diagnostics = false;
    std::vector<n::ArenaRateEntry> rows;
    Vec E, B;
    double begin = t0;
    std::uint64_t step = 0;
    const hd::MotionReceipt* external = nullptr;
    hd::MotionBinding external_binding;
    Json extra = Json::object();
};

void process(const std::filesystem::path& out_root, const axial::Mesh& g, const n::GapPicOptions& o, Case& k, Json& summary) {
    const auto dir = out_root / k.name;
    std::filesystem::create_directories(dir);
    Json j{{"fixture", k.name}, {"carriers", k.rows.size()}, {"mode", k.diagnostics ? "diagnostics" : "record"}};
    j.update(k.extra);
    const auto mt = ha::mesh_tables(g);
    const auto geo = ha::geometry(o);
    auto t = Clock::now();
    const auto census = ha::order_census(k.rows, false);
    j["host_order_census_ms"] = ms_since(t);
    if (census.packet_failure_rank >= 0) throw std::runtime_error("fixture has divergent Count2 packets");
    const auto binding = k.external ? k.external_binding : make_binding(g, o, k.name, k.rows, k.E, k.B, k.begin, k.step);
    t = Clock::now();
    const auto census_sha = census_sha256(k.rows);
    j["census_sha256_ms"] = ms_since(t);
    j["census_sha256"] = census_sha;
    write_input((dir / "input.bin").string(), mt, geo, census, k.E, k.B, k.diagnostics, binding, census_sha);
    // sealed reference
    Sealed sealed;
    if (k.rows.size() <= 100000) {
        sealed = sealed_run(g, o, k.rows, binding, k.E, k.B, k.external);
        hgp_enc::Writer w;
        sealed.rec.write(w, false);
        w.save((dir / "sealed.bin").string());
        j["sealed"] = {{"ran", true}, {"refusal", sealed.rec.refusal}, {"charge_refusal", sealed.rec.charge_refusal}, {"seal_sha256", sealed.seal},
                       {"census_sha256", sealed.census}, {"timing", sealed.timing}, {"external_receipt", k.external != nullptr}};
        if (!sealed.census.empty()) j["sealed"]["census_sha256_equals_uncapped_helper"] = sealed.census == census_sha;
    } else {
        std::string why;
        try { (void)hd::prepare_canonical_motion(g, o, k.rows, binding, {}, k.E, k.B); } catch (const std::exception& e) { why = e.what(); }
        j["sealed"] = {{"ran", false}, {"refusal", why}};
    }
    // replica
    t = Clock::now();
    const auto rep = replica(g, o, census, k.E, k.B, k.begin, k.diagnostics);
    j["replica"] = {{"total_ms", ms_since(t)}, {"phases_ms", rep.timing}, {"refusal", rep.rec.refusal}, {"refusal_rank", rep.rec.refusal_rank},
                    {"charge_refusal", rep.rec.charge_refusal}};
    write_result((dir / "replica.bin").string(), rep, k.diagnostics);
    if (rep.rec.refusal.empty()) {
        t = Clock::now();
        const auto seal = record_seal(binding, o, census_sha, rep.rec);
        j["replica"]["record_seal_sha256"] = seal;
        j["replica"]["record_seal_ms"] = ms_since(t);
        if (!sealed.seal.empty()) j["replica"]["record_seal_equals_sealed"] = seal == sealed.seal;
    }
    // host port
    t = Clock::now();
    const auto port = host_port(mt, geo, census, k.E, k.B, k.diagnostics);
    j["hostport"] = {{"total_ms", ms_since(t)}, {"refusal", port.rec.refusal}, {"refusal_rank", port.rec.refusal_rank}};
    write_result((dir / "hostport.bin").string(), port, k.diagnostics);
    std::ofstream((dir / "cpu.json").string()) << j.dump(2) << '\n';
    summary[k.name] = j;
    std::cout << k.name << " carriers=" << k.rows.size() << " sealed_refusal='" << (j["sealed"].value("refusal", std::string())) << "' replica_refusal='"
              << rep.rec.refusal << "' replica_ms=" << j["replica"]["total_ms"] << std::endl;
}

int run(const std::filesystem::path& trace, const std::filesystem::path& out, const std::set<std::string>& only) {
    std::filesystem::create_directories(out);
    auto want = [&](const std::string& name) { return only.empty() || only.count(name); };
    auto model = std::make_shared<const h::Model>(h::MeshSpec{}, trace);
    const auto& g = model->mesh();
    const auto o = geometry(g);
    const auto seed = model->dc_seed();
    const auto mid0 = model->at(seed);   // iteration-1 midpoint of prepare_closed_proposal
    Json summary = Json::object();
    summary["model_configuration_sha256"] = model->configuration_sha256();
    // ---------------- field-port fixture: co-located electron/proton pairs, prepare_closed_proposal
    for (double W : {1e9, 1e11, 1e12}) {
        char tag[32];
        std::snprintf(tag, sizeof tag, "%.0e", W);
        const std::string base = std::string("port_pair_w") + tag;
        if (!want(base + "_seed") && !want(base + "_final")) continue;
        const std::vector<n::ArenaRateEntry> rows{carrier(1, o, audit_point, {2e6, -1e6, 5e5}, -n::pic_e, n::pic_me, W, "pair-electron"),
                                                   carrier(2, o, audit_point, {0, 0, 0}, n::pic_e, hall::proton_mass, W, "pair-proton")};
        const auto ledger = dc_declaration(g, seed, h::charge_of_census(g, o, rows));
        h::Owner owner(model, seed, ledger, o, rows);
        hd::MotionBinding b;
        b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
        b.accepted_step = owner.revision();
        b.staged_arena_revision = 1;
        b.begin_s = owner.time_s();
        b.dt_s = dt;
        b.end_s = owner.time_s() + dt;
        b.staged_arena_sha256 = hd::canonical_census_sha256(rows);
        b.mesh_sha256 = owner.motion_mesh_sha256();
        b.material_sha256 = hd::material_binding_sha256();
        b.species_sha256 = hd::text_sha256("audit port fixture: electron and proton carriers");
        b.common_clock_sha256 = hd::text_sha256(Json{{"begin", b.begin_s}, {"dt", dt}, {"end", b.end_s}, {"step", b.accepted_step}, {"fixture", "audit port"}}.dump());
        b.kinetic_policy_sha256 = hd::kinetic_policy_sha256(hd::declared_kinetic_policy(o));
        const auto t = Clock::now();
        auto proposal = h::prepare_closed_proposal(owner, o, rows, b, h::sha256("audit port: unchanged closed source identity"));
        const double proposal_ms = ms_since(t);
        const auto& mid = proposal.midpoint();
        Json extra{{"weight", W}, {"closed_proposal_iterations", proposal.iterations()}, {"closed_proposal_ms", proposal_ms},
                   {"candidate_field_sha256", proposal.candidate_field_sha256()}, {"owner_preflight", owner.preflight(proposal, rows)}};
        if (want(base + "_final")) {
            Case k;
            k.name = base + "_final";
            k.rows = rows; k.E = mid.E(); k.B = mid.B(); k.begin = b.begin_s; k.step = b.accepted_step;
            k.external = &proposal.motion();
            k.external_binding = proposal.motion().binding();
            k.extra = extra;
            k.extra["midpoint"] = "authoritative compensated midpoint of the accepted closed proposal";
            process(out, g, o, k, summary);
            // field arrays for the GPU-in-the-loop check
            hgp_enc::Writer w;
            w.f64("seed.E", seed.E); w.f64("seed.B", seed.B);
            w.text("candidate_field_sha256", proposal.candidate_field_sha256());
            const std::int32_t it[1] = {std::int32_t(proposal.iterations())};
            w.i32("iterations", it);
            w.save((out / k.name / "closed_proposal.bin").string());
        }
        if (want(base + "_seed")) {
            Case k;
            k.name = base + "_seed";
            k.rows = rows; k.E = mid0.E(); k.B = mid0.B(); k.begin = b.begin_s; k.step = b.accepted_step;
            k.external_binding = b;
            k.external_binding.trial_field_sha256 = hd::total_field_sha256(mid0.E(), mid0.B());
            k.extra = extra;
            k.extra["midpoint"] = "Model::at(dc_seed): first-iteration midpoint of prepare_closed_proposal";
            // sealed run through prepare_canonical_motion with the first-iteration binding
            const auto t2 = Clock::now();
            auto r = hd::prepare_canonical_motion(g, o, rows, k.external_binding, {}, k.E, k.B);
            k.extra["prepare_canonical_motion_ms"] = ms_since(t2);
            k.external = &r;
            process(out, g, o, k, summary);
        }
    }
    // ---------------- synthetic populations on the DC-seed midpoint
    for (std::size_t N : {std::size_t(1000), std::size_t(10000), std::size_t(100000), std::size_t(1000000)}) {
        const std::string tag = N == 1000 ? "1e3" : N == 10000 ? "1e4" : N == 100000 ? "1e5" : "1e6";
        if (want("syn_record_" + tag)) {
            Generator gen(g, o, 0x5eed0000ull + N);
            Case k;
            k.name = "syn_record_" + tag;
            k.rows = gen.population(N, .2, 0);
            k.E = mid0.E(); k.B = mid0.B();
            k.extra = {{"population", "40% electrons 0.1 eV-10 keV, 20% protons, 20% H2O+ 0.01 eV-10 keV, 20% H2O 0.01-1 eV; weights 1e3-1e8; 20% placed to cross 1-2 interior cell faces"},
                       {"field", "Model::at(dc_seed)"}};
            process(out, g, o, k, summary);
        }
        if (want("syn_boundary_" + tag)) {
            Generator gen(g, o, 0xb0d00000ull + N);
            Case k;
            k.name = "syn_boundary_" + tag;
            k.diagnostics = true;
            k.rows = gen.population(N, .2, .02);
            k.E = mid0.E(); k.B = mid0.B();
            k.extra = {{"population", "as syn_record plus 2% placed to reach foil faces, outer xy walls, short/feed planes or the short bridge within the step"},
                       {"field", "Model::at(dc_seed)"}};
            process(out, g, o, k, summary);
        }
    }
    if (want("syn_metal_1e3")) {
        Generator gen(g, o, 0x3e7a1000ull);
        Case k;
        k.name = "syn_metal_1e3";
        k.diagnostics = true;
        k.rows = gen.population(1000, .2, .02);
        const auto pol = hd::declared_kinetic_policy(o);
        k.rows[611].native.position_m = {.5 * (pol.upper.left_m + pol.upper.right_m), .5 * (pol.upper.bottom_m + pol.upper.top_m), .001};
        k.E = mid0.E(); k.B = mid0.B();
        k.extra = {{"population", "syn_boundary population of 1000 with carrier rank 611 moved inside the upper foil (census-loop refusal)"}};
        process(out, g, o, k, summary);
    }
    std::ofstream((out / "cpu-summary.json").string()) << summary.dump(2) << '\n';
    return 0;
}

int seal_mode(const std::filesystem::path& fixture, const std::filesystem::path& encoding) {
    hgp_enc::Reader in((fixture / "input.bin").string()), r(encoding.string());
    hd::MotionBinding b;
    b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
    const auto steps = in.u64("binding.steps");
    const auto times = in.f64("binding.times");
    b.accepted_step = steps[0]; b.staged_arena_revision = steps[1]; b.trial_stage = steps[2];
    b.begin_s = times[0]; b.end_s = times[1]; b.dt_s = times[2];
    b.staged_arena_sha256 = in.text("binding.staged_arena_sha256"); b.trial_field_sha256 = in.text("binding.trial_field_sha256");
    b.mesh_sha256 = in.text("binding.mesh_sha256"); b.material_sha256 = in.text("binding.material_sha256");
    b.species_sha256 = in.text("binding.species_sha256"); b.common_clock_sha256 = in.text("binding.common_clock_sha256");
    b.kinetic_policy_sha256 = in.text("binding.kinetic_policy_sha256");
    if (!r.text("refusal.motion").empty()) { std::cout << "{\"seal_sha256\":null,\"refusal\":" << Json(r.text("refusal.motion")).dump() << "}\n"; return 0; }
    hgp_enc::RecordArrays a;
    a.id = r.u64("endpoint.native_id"); a.end_pos = r.f64("endpoint.position_m"); a.end_u = r.f64("endpoint.proper_velocity_m_s");
    a.charge = r.f64("endpoint.charge_C"); a.mass = r.f64("endpoint.mass_kg"); a.weight = r.f64("endpoint.weight");
    a.start = r.f64("path.field_start_m"); a.disp = r.f64("path.field_displacement_m"); a.flight = r.f64("path.flight_s");
    a.current = r.f64("current_A_m"); a.initial_charge = r.f64("initial_charge_C"); a.endpoint_charge = r.f64("endpoint_charge_C");
    a.budget = r.f64("budget");
    const auto mreal = in.f64("mesh.real");
    const auto mint = in.i32("mesh.int");
    const axial::Mesh g(mint[0], mint[2] - 2, mreal[0], 1);
    const auto o = geometry(g);
    std::cout << "{\"seal_sha256\":" << Json(record_seal(b, o, in.text("census_sha256"), a)).dump() << "}\n";
    return 0;
}

int selftest(const std::filesystem::path& out) {
    // gross_step (x87 emulation) against the long double statement of deposit_charge.
    Rng rng{0x80b17ull};
    std::size_t cases = 0, mismatches = 0;
    Json first = nullptr;
    auto random_double = [&]() {
        const int mode = int(rng.next() % 8);
        if (mode == 0) return std::ldexp(rng.uniform(), -1074 + int(rng.next() % 60));   // subnormal range
        if (mode == 1) return double(rng.next() % 1000) * 0x1p-20;
        const int e = int(rng.next() % 600) - 300;
        return std::ldexp(1 + rng.uniform(), e);
    };
    for (int k = 0; k < 2000000; ++k) {
        double g0 = random_double();
        h::Words v(random_double() * (rng.uniform() < .5 ? -1 : 1), 0.0);
        if (rng.uniform() < .7) { const double lo = std::ldexp(rng.uniform() - .5, std::ilogb(v.hi == 0 ? 1.0 : v.hi) - 53); v = h::Words(v.hi, lo); }
        if (!h::canonical(v)) continue;
        if (rng.uniform() < .1) g0 = std::abs(v.hi) * std::ldexp(1.0, int(rng.next() % 140) - 70);
        const long double gross = (long double)g0 + h::magnitude(v);
        double ref = double(gross);
        if ((long double)ref < gross) ref = std::nextafter(ref, INFINITY);
        const double emu = hgp::gross_step(g0, {v.hi, v.lo});
        ++cases;
        if (std::bit_cast<std::uint64_t>(ref) != std::bit_cast<std::uint64_t>(emu)) {
            ++mismatches;
            if (first.is_null()) first = {{"g", g0}, {"hi", v.hi}, {"lo", v.lo}, {"reference", ref}, {"emulated", emu}};
        }
    }
    // Words arithmetic against field_contract.hpp.
    std::size_t wcases = 0, wmismatch = 0;
    for (int k = 0; k < 1000000; ++k) {
        const double a = random_double() * (rng.uniform() < .5 ? -1 : 1), b = random_double() * (rng.uniform() < .5 ? -1 : 1);
        const double c = rng.uniform();
        const auto ref = (h::Words(a) * b) * c + h::Words(b, a);
        const auto emu = hgp::words_add(hgp::words_mul(hgp::words_mul({a, 0}, b), c), hgp::two_sum_words(b, a));
        ++wcases;
        if (std::bit_cast<std::uint64_t>(ref.hi) != std::bit_cast<std::uint64_t>(emu.hi) || std::bit_cast<std::uint64_t>(ref.lo) != std::bit_cast<std::uint64_t>(emu.lo)) ++wmismatch;
    }
    Json j{{"x87_gross_step_cases", cases}, {"x87_gross_step_mismatches", mismatches}, {"first_mismatch", first},
           {"words_cases", wcases}, {"words_mismatches", wmismatch}, {"long_double_digits", std::numeric_limits<long double>::digits}};
    std::filesystem::create_directories(out);
    std::ofstream((out / "selftest.json").string()) << j.dump(2) << '\n';
    std::cout << j.dump() << std::endl;
    return mismatches || wmismatch ? 1 : 0;
}
}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rc = 0;
    try {
        const std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "run" && argc >= 4) {
            std::set<std::string> only;
            for (int i = 4; i < argc; ++i) only.insert(argv[i]);
            rc = run(argv[2], argv[3], only);
        } else if (mode == "seal" && argc == 4) rc = seal_mode(argv[2], argv[3]);
        else if (mode == "selftest" && argc == 3) rc = selftest(argv[2]);
        else throw std::runtime_error("usage: cpu_reference run <csv> <out> [fixture...] | seal <fixture dir> <encoding> | selftest <out>");
    } catch (const std::exception& e) {
        std::cerr << "cpu_reference aborted: " << e.what() << '\n';
        rc = 1;
    }
    MPI_Finalize();
    return rc;
}
