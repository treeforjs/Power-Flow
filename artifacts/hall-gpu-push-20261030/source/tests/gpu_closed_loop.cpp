// gpu_closed_loop.cpp -- the fixed point of hall_corrective::prepare_closed_proposal
// with the motion and endpoint free charge taken from the GPU engine, compared with
// the sealed proposal on the field-port fixture (co-located electron/proton pairs).
// The loop transcribes root_api.hpp lines 48-58; the field step is the sealed
// Model::prepare on the CPU. Usage: gpu_closed_loop <C19 csv> <nvrtc lib> <uuid> <out json>
#include "root_api.hpp"
#include "push/hall_adapter.hpp"

#include <mpi.h>

#include <chrono>
#include <fstream>
#include <iostream>

namespace h = hall_corrective;
namespace n = raikou::mfem_backend;
namespace hd = hall_driver;
namespace hp = hall_gpu::push;
namespace ha = hall_gpu::push::hall_adapter;
using h::Json;

namespace {
constexpr double dt = 1e-14;
const axial::V3 audit_point{.00021, -.00017, .001};
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
n::ArenaRateEntry carrier(std::uint64_t id, const n::GapPicOptions& o, n::PicVector u, double charge, double mass, double weight, const std::string& tag) {
    n::ArenaRateEntry row;
    row.packet_id = "audit-port/packet/" + tag;
    row.state = "audit-port/resident/" + tag;
    row.native.id = id;
    row.native.position_m = {audit_point[0] + o.width_m / 2, audit_point[1] + o.gap_m / 2, audit_point[2]};
    row.native.proper_velocity_m_s = u;
    row.native.charge_C = charge;
    row.native.mass_kg = mass;
    row.native.weight = weight;
    return row;
}
h::ChargeLedger dc_declaration(const axial::Mesh& g, const axial::State& s, h::FreeCharge supplied) {
    h::ChargeLedger q;
    q.particles = std::move(supplied);
    q.conductor_C.assign(g.nn, h::Words());
    q.exterior_C.assign(g.nn, h::Words());
    q.conductor_declaration = h::sha256("push parity: DC induced conductor/PEC charge declared from the supplied field");
    q.exterior_declaration = h::sha256("push parity: DC exterior termination charge declared from the supplied field");
    const auto Q = h::field_charge(g, s);
    for (int k = 0; k < g.nn; ++k) {
        if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = Q[k];
        else if (g.outer_node[k]) q.exterior_C[k] = Q[k];
    }
    return q;
}
}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rc = 0;
    Json out = Json::object();
    try {
        h::need(argc == 5, "usage: gpu_closed_loop <C19 csv> <nvrtc lib> <uuid> <out json>");
        auto model = std::make_shared<const h::Model>(h::MeshSpec{}, argv[1]);
        const auto& g = model->mesh();
        const auto o = geometry(g);
        const auto seed = model->dc_seed();
        hall_gpu::SessionOptions so;
        so.require_uuid = argv[3];
        so.nvrtc_library = argv[2];
        hall_gpu::Session session(so);
        out["device"] = Json::parse(hall_gpu::device_json(session.device()));
        Json cases = Json::array();
        {
            hp::EngineOptions eo;
            eo.capacity = 2;
            hp::MotionEngine engine(session, ha::mesh_tables(g), ha::geometry(o), eo);
            out["compile"] = Json::parse(hall_gpu::compile_json(engine.compile_record()));
            for (double W : {1e9, 1e11, 1e12}) {
                const std::vector<n::ArenaRateEntry> rows{carrier(1, o, {2e6, -1e6, 5e5}, -n::pic_e, n::pic_me, W, "pair-electron"),
                                                           carrier(2, o, {0, 0, 0}, n::pic_e, hall::proton_mass, W, "pair-proton")};
                h::Owner owner(model, seed, dc_declaration(g, seed, h::charge_of_census(g, o, rows)), o, rows);
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
                auto sealed = h::prepare_closed_proposal(owner, o, rows, b, h::sha256("audit port: unchanged closed source identity"));
                const auto census = ha::order_census(rows, true);
                engine.set_carriers(census.view());
                const auto sigma = model->initial_conductivity();
                auto midpoint = model->at(owner.field());
                Json c{{"weight", W}, {"sealed_iterations", sealed.iterations()}, {"sealed_candidate_field_sha256", sealed.candidate_field_sha256()}};
                bool converged = false;
                const auto t = std::chrono::steady_clock::now();
                for (unsigned iteration = 0; iteration < 16; ++iteration) {
                    engine.set_field(midpoint.E(), midpoint.B());
                    const auto v = engine.run(dt, false);
                    h::need(!v.refusal.refused && !v.charge_refusal.refused, "GPU motion refused inside the closed loop");
                    axial::Vec current(v.current_A_m.begin(), v.current_A_m.end());
                    h::FreeCharge after{h::WordVec(g.nn), axial::Vec(v.free_gross_C.begin(), v.free_gross_C.end())};
                    for (int k = 0; k < g.nn; ++k) { after.net_C[k].hi = v.free_net_hi[k]; after.net_C[k].lo = v.free_net_lo[k]; }
                    // The sealed motion of the same midpoint, for an exact per-iteration comparison.
                    auto binding = b;
                    binding.trial_field_sha256 = hd::total_field_sha256(midpoint.E(), midpoint.B());
                    const auto cpu = hd::prepare_canonical_motion(g, o, rows, binding, {}, midpoint.E(), midpoint.B());
                    const auto cpu_after = h::charge_of_endpoints(g, o, cpu.endpoints());
                    bool same = cpu.current_A_m() == current && cpu_after.gross_C == after.gross_C;
                    for (int k = 0; k < g.nn && same; ++k)
                        same = std::bit_cast<std::uint64_t>(cpu_after.net_C[k].hi) == std::bit_cast<std::uint64_t>(after.net_C[k].hi) &&
                               std::bit_cast<std::uint64_t>(cpu_after.net_C[k].lo) == std::bit_cast<std::uint64_t>(after.net_C[k].lo);
                    c["iteration_" + std::to_string(iteration + 1) + "_gpu_current_and_free_charge_bitwise_equal_sealed"] = same;
                    auto field = model->prepare(owner.field(), owner.charge(), current, after, sigma, dt);
                    if (midpoint.E() == field.midpoint().E() && midpoint.B() == field.midpoint().B()) {
                        c["gpu_iterations"] = iteration + 1;
                        c["gpu_candidate_field_sha256"] = field.after_sha256();
                        c["candidate_field_sha256_equal"] = field.after_sha256() == sealed.candidate_field_sha256();
                        c["midpoint_E_bitwise_equal"] = field.midpoint().E() == sealed.midpoint().E();
                        c["midpoint_B_bitwise_equal"] = field.midpoint().B() == sealed.midpoint().B();
                        c["midpoint_seal_equal"] = field.midpoint().seal() == sealed.midpoint().seal();
                        converged = true;
                        break;
                    }
                    midpoint = field.midpoint();
                }
                c["gpu_loop_converged"] = converged;
                c["gpu_loop_wall_ms"] = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count();
                c["iterations_equal"] = converged && c["gpu_iterations"] == sealed.iterations();
                cases.push_back(c);
            }
        }
        out["cases"] = cases;
        out["teardown"] = Json::parse(hall_gpu::teardown_json(session.close()));
        bool ok = true;
        for (const auto& c : cases) ok = ok && c.value("candidate_field_sha256_equal", false) && c.value("iterations_equal", false);
        out["all_equal"] = ok;
        if (!ok) rc = 2;
    } catch (const std::exception& e) {
        out["error"] = e.what();
        std::cerr << "gpu_closed_loop failed: " << e.what() << '\n';
        rc = 1;
    }
    std::ofstream(argv[argc - 1]) << out.dump(2) << '\n';
    std::cout << out.dump(2) << std::endl;
    MPI_Finalize();
    return rc;
}
