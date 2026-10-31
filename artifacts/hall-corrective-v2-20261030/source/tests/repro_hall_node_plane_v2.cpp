// Stage 1 node-plane reproducer (source/repro/repro_hall_node_plane.cpp of the
// Stage 1 artifact) run against Hall corrective v2. Identical cases and scan;
// only the namespace changes to hall_corrective_v2 and the scan reports the v2
// tolerance (relative term plus the deposit-rounding absolute term).
// Expected with v2: every case, including z = 0, is admitted.
// For each carrier height z it reports whether prepare_closed_proposal admits
// the step and, independently, the worst per-node ratio
//   |after - before - incidence(current) dt| / (1e-12 (gross_before + gross_after + |incidence dt|))
// which is exactly the guard "actual particle endpoint charge/current mismatch".
#include "root_api.hpp"
#include <iomanip>
#include <iostream>
namespace hc = hall_corrective_v2;
namespace hd = hall_driver;
namespace n = raikou::mfem_backend;
using Json = hc::Json;

n::GapPicOptions geometry(const axial::Mesh& g) {
  n::GapPicOptions o; o.nx = o.ny = g.nx; o.width_m = o.gap_m = g.side; o.depth_m = g.depth;
  o.finite_foils = o.finite_axial_extent = true; o.short_z_m = -g.depth / 2; o.feed_z_m = g.depth / 2;
  o.lower_foil = {g.side / 2 - .0025, g.side / 2 + .0025, g.side / 2 - .0027, g.side / 2 - .0025};
  o.upper_foil = {g.side / 2 - .0025, g.side / 2 + .0025, g.side / 2 + .0025, g.side / 2 + .0027};
  return o;
}
hc::ChargeLedger declared_dc(const hc::Model& m, const axial::State& s) {
  const auto& g = m.mesh(); hc::ChargeLedger q;
  q.particles.net_C.resize(g.nn); q.particles.gross_C.resize(g.nn); q.conductor_C.resize(g.nn); q.exterior_C.resize(g.nn);
  q.conductor_declaration = hc::sha256("reproducer: declared DC conductor/PEC induced charge");
  q.exterior_declaration = hc::sha256("reproducer: declared DC exterior termination charge");
  auto field = hc::field_charge(g, s);
  for (int node = 0; node < g.nn; ++node) {
    if (g.node_body[node] >= 0 || g.short_node[node]) q.conductor_C[node] = field[node];
    else if (g.outer_node[node]) q.exterior_C[node] = field[node];
  }
  return q;
}
n::ArenaRateEntry row(std::uint64_t id, const char* packet, const char* state, n::PicVector x, n::PicVector u,
                      double charge, double mass, double weight) {
  n::ArenaRateEntry r; r.packet_id = packet; r.state = state; r.native.id = id; r.native.position_m = x;
  r.native.proper_velocity_m_s = u; r.native.charge_C = charge; r.native.mass_kg = mass; r.native.weight = weight;
  return r;
}
int main(int argc, char** argv) {
  if (argc != 2) { std::cerr << "usage: repro C19_TRACE\n"; return 2; }
  auto model = std::make_shared<const hc::Model>(hc::MeshSpec{}, argv[1]);
  const auto& g = model->mesh(); const auto o = geometry(g); const auto seed = model->dc_seed();
  const double mn = 2.9915e-26, me = 9.1093837139e-31, qe = n::pic_e, c = n::pic_c, threshold = 12.621;
  const double ion_mass = mn - me + threshold * qe / (c * c);
  const double x100 = 100 * qe / (me * c * c);
  const n::PicVector ue{c * std::sqrt(x100 * (x100 + 2)), 0, 0}, ui{0, 372, 0};
  const double dt = 1e-14;
  Json out = Json::array();
  bool expected_pattern = true;
  struct Case { double z, weight; const char* note; };
  const Case cases[] = {{0., 1., "carriers on the k=3 Hall node plane z=0 (Stage 1 plan position)"},
                        {0., 1e9, "same plane, weight 1e9"},
                        {1e-12, 1., "1e-12 m above the plane"}, {1e-9, 1., "1e-9 m above"},
                        {1e-7, 1., "1e-7 m above"}, {1e-6, 1., "1e-6 m above"},
                        {1e-5, 1., "1e-5 m above"}, {1e-4, 1., "1e-4 m above"},
                        {1e-3, 1., "1e-3 m above (generic interior)"}};
  for (const auto& cs : cases) {
    Json r{{"z_m", cs.z}, {"weight", cs.weight}, {"note", cs.note}};
    const n::PicVector x{.0045, .0045, cs.z};
    const std::vector<n::ArenaRateEntry> census{row(1, "e0", "electron", x, ue, -qe, me, cs.weight),
                                                row(2, "i0", "H2O+", x, ui, qe, ion_mass, cs.weight)};
    const double vz = (x[2] - g.zmin) / g.dz;
    r["z_in_cell_units"] = vz; r["distance_to_nearest_z_node_plane_cells"] = std::abs(vz - std::round(vz));
    auto ledger = declared_dc(*model, seed); ledger.particles = hc::charge_of_census(g, o, census);
    hc::Owner owner(model, seed, ledger, o, census);
    hd::MotionBinding b; b.phase = hd::ParticlePhase::full_time_position_and_proper_velocity;
    b.accepted_step = owner.revision(); b.staged_arena_revision = 1; b.begin_s = owner.time_s(); b.dt_s = dt;
    b.end_s = b.begin_s + dt; b.staged_arena_sha256 = hc::sha256("reproducer: staged arena");
    b.mesh_sha256 = owner.motion_mesh_sha256(); b.material_sha256 = model->material_sha256();
    b.species_sha256 = hc::sha256("reproducer: species"); b.common_clock_sha256 = hc::sha256("reproducer: clock");
    b.kinetic_policy_sha256 = owner.motion_context().policy_sha256;
    try {
      auto p = hc::prepare_closed_proposal(owner, o, census, b, hc::sha256("reproducer: unchanged closed source"));
      r["prepare_closed_proposal"] = "admitted"; r["iterations"] = p.iterations();
      r["v2_continuity_defect_over_tolerance"] = p.timings().continuity_ratio;
    } catch (const std::exception& e) { r["prepare_closed_proposal"] = std::string("refused: ") + e.what(); }
    // Independent scan of the same guard on the first fixed-point iterate.
    const auto mid = model->at(owner.field());
    auto bound = b; bound.trial_field_sha256 = hd::total_field_sha256(mid.E(), mid.B());
    const auto motion = hd::prepare_canonical_motion(g, o, census, bound, {}, mid.E(), mid.B());
    const auto after = hc::charge_of_endpoints(g, o, motion.endpoints());
    const auto divergence = hc::incidence(g, motion.current_A_m());
    const auto support = hc::DepositSupportBuilder::build(g, o, census, motion.endpoints());
    double worst_v2 = 0;
    double worst = 0; int worst_node = -1; long violating = 0;
    for (int node = 0; node < g.nn; ++node) {
      const auto defect = after.net_C[node] - ledger.particles.net_C[node] - divergence[node] * dt;
      const long double scale = ledger.particles.gross_C[node] + after.gross_C[node] + hc::magnitude(divergence[node] * dt);
      const long double d = std::abs(defect.wide());
      if (d == 0) continue;
      worst_v2 = std::max(worst_v2, double(d / (1e-12L * scale + (long double)support.absolute_C()[node])));
      const double ratio = scale > 0 ? double(d / (1e-12L * scale)) : INFINITY;
      if (ratio > 1) ++violating;
      if (ratio > worst) { worst = ratio; worst_node = node; }
    }
    r["first_iterate_worst_defect_over_tolerance"] = worst;
    r["first_iterate_worst_defect_over_v2_tolerance"] = worst_v2;
    if (r["prepare_closed_proposal"] != "admitted") expected_pattern = false; r["first_iterate_nodes_violating"] = violating;
    r["worst_node"] = worst_node;
    if (worst_node >= 0) {
      const auto defect = after.net_C[worst_node] - ledger.particles.net_C[worst_node] - divergence[worst_node] * dt;
      r["worst_node_defect_C"] = double(defect.wide());
      r["worst_node_gross_before_C"] = ledger.particles.gross_C[worst_node];
      r["worst_node_gross_after_C"] = after.gross_C[worst_node];
      r["worst_node_k"] = worst_node / ((g.nx + 1) * (g.ny + 1));
    }
    r["motion_continuity_C"] = motion.budget().continuity_C;
    r["electron_displacement_m"] = {motion.paths()[0].field_displacement_m[0], motion.paths()[0].field_displacement_m[1],
                                    motion.paths()[0].field_displacement_m[2]};
    out.push_back(r);
    std::cout << std::setprecision(6) << "z=" << cs.z << " w=" << cs.weight << " -> " << r["prepare_closed_proposal"].get<std::string>()
              << " | worst ratio (8dddec0 tolerance) " << worst << " violating nodes " << violating
              << " | worst ratio (v2 tolerance) " << worst_v2 << std::endl;
  }
  std::cout << out.dump(1) << std::endl;
  std::cout << (expected_pattern ? "all cases admitted" : "some case refused") << std::endl;
  return expected_pattern ? 0 : 1;
}
