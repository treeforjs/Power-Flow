// motion-work-gross-v1 in native v7 admission (API-007.txt 1.1 amendment).
// Dense-blob fixtures: the coordinator v7 / Hall v3.2 dense generator (seed
// 20261101, 400 co-located electron/H2O+ pairs, electrons 3 eV, ions 0.025 eV)
// at weights 3.45e10 and 8e10, moved by the genuine Hall corrective v2
// provider (the sealed 8dddec0 work-kick push) over dt = 10 fs in a fixed
// trial field on the default Hall axial mesh, entering native v7 through the
// transitional adapter. Field cases: uniform (E_y = 7.3e4 V/m) and split
// (+-7.3e4 V/m on the two x halves of the blob: small net exchange).
// Checks: native admits; native's tolerance and the five scales are bitwise
// equal to Hall's reference (hall-corrective-v3.2 source/motion_work_gross_v1.hpp,
// read only) on the same inputs, also for every tampered input set; tampers
// refused with the stated margins: one electron's kinetic change off by 1e-3
// of 1.5 kT w, the work field 1e-6 (relative) away from the push field, one
// carrier's mean gross exchange with the wrong sign, one ulp of K, one ulp of W.
// usage: test_v7_work_gross OUTPUT_DIRECTORY
#include "device_mode.hpp"
#include "v6_fixture.hpp"
#include "hall-finite-foil-em-20261029/axial_mesh.hpp"
#include "motion_receipt.hpp" // Hall corrective v2
#include "hall_receipt_v7.hpp"
#include HALL_V32_WORK_GROSS_HEADER
#include <random>

using namespace v6;
using FK = m::FullTimeFlowRefusal::Kind;
namespace {
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }
const axial::Mesh &mesh() {
  static const axial::Mesh g;
  return g;
}
constexpr double me = 9.1093837139e-31, mn = 2.9915e-26, threshold_eV = 12.621;
double ion_mass() { return mn - me + threshold_eV * m::pic_e / (m::pic_c * m::pic_c); }
std::vector<m::ArenaStateSpec> blob_states() {
  return {{"H2O+", 2, 1, 1, ion_mass(), 0, 0}, {"e", 0, 0, -1, m::pic_me, 0, 0}, {"n", 2, 1, 0, mn, 0, 0}};
}
// The coordinator v7 dense generator (Hall v3.2 test_v32 blob_census).
std::vector<m::FullTimeAuthoredPacket> blob(double weight, double eV, long pairs) {
  std::mt19937_64 rng(20261101);
  std::uniform_real_distribution<double> u(0, 1);
  std::normal_distribution<double> nd(0, 1);
  const double se = std::sqrt(eV * m::pic_e / me), si = std::sqrt(.025 * m::pic_e / mn);
  std::vector<m::FullTimeAuthoredPacket> out;
  const m::ArenaCount w(weight, 0);
  for (long k = 0; k < pairs; ++k) {
    const m::PicVector x{.0026 + 8e-4 * u(rng), .0051 + 8e-4 * u(rng), -.003 + .006 * u(rng)};
    const m::PicVector ve{se * nd(rng), se * nd(rng), se * nd(rng)}, vi{si * nd(rng), si * nd(rng), si * nd(rng)};
    out.push_back({"be" + std::to_string(k), "e", "dense-blob", w, x, ve});
    out.push_back({"bi" + std::to_string(k), "H2O+", "dense-blob", w, x, vi});
  }
  return out;
}
struct Field {
  std::vector<double> E, B;
};
Field field(bool split) {
  const auto &g = mesh();
  Field f{std::vector<double>(g.ne, 0.0), std::vector<double>(g.nb, 0.0)};
  std::mt19937_64 rng(7);
  std::normal_distribution<double> nd(0, 1);
  for (int e = 0; e < g.ne; ++e)
    f.E[e] = 1e3 * nd(rng);
  // y edges: 7.3e4 V/m; split: the sign flips at field-frame x = -1.5 mm
  // (native x = 3.0 mm, the middle of the blob).
  for (int k = 0; k <= g.nz; ++k)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i <= g.nx; ++i)
        f.E[g.ey(i, j, k)] = (split && i * g.dx - g.side / 2 >= -.0015) ? -7.3e4 : 7.3e4;
  return f;
}
// Hall's reference scales (motion_work_gross_v1.hpp) on explicit inputs.
motion_work_v1::Scales hall_scales(std::span<const pf_v7::Carrier> census, std::span<const pf_v7::Carrier> ends,
                                   std::span<const pf_v7::Sink> sinks, double dt, std::span<const double> E,
                                   std::span<const double> J, double K, double W) {
  motion_work_v1::CarrierSums c;
  std::size_t j = 0, s = 0;
  for (const auto &x : census) {
    const double *u = nullptr;
    if (s < sinks.size() && sinks[s].id == x.id)
      u = sinks[s++].incident_proper_velocity_m_s;
    else
      u = ends[j++].proper_velocity_m_s;
    const double K0 = pf7::kinetic_energy(x.mass_kg, x.proper_velocity_m_s[0], x.proper_velocity_m_s[1],
                                          x.proper_velocity_m_s[2]);
    c.add(x.weight, x.charge_C, K0, pf7::kinetic_energy(x.mass_kg, u[0], u[1], u[2]) - K0);
  }
  const auto &g = mesh();
  std::vector<int> a(g.ne), b(g.ne);
  for (int k = 0; k <= g.nz; ++k)
    for (int jj = 0; jj <= g.ny; ++jj)
      for (int i = 0; i < g.nx; ++i) {
        a[g.ex(i, jj, k)] = g.node(i, jj, k);
        b[g.ex(i, jj, k)] = g.node(i + 1, jj, k);
      }
  for (int k = 0; k <= g.nz; ++k)
    for (int jj = 0; jj < g.ny; ++jj)
      for (int i = 0; i <= g.nx; ++i) {
        a[g.ey(i, jj, k)] = g.node(i, jj, k);
        b[g.ey(i, jj, k)] = g.node(i, jj + 1, k);
      }
  for (int k = 0; k < g.nz; ++k)
    for (int jj = 0; jj <= g.ny; ++jj)
      for (int i = 0; i <= g.nx; ++i) {
        a[g.ez(i, jj, k)] = g.node(i, jj, k);
        b[g.ez(i, jj, k)] = g.node(i, jj, k + 1);
      }
  motion_work_v1::Scales s0;
  motion_work_v1::field_scales(s0, dt, E, J, a, b, std::size_t(g.nn));
  motion_work_v1::finish(s0, c, K, W);
  return s0;
}
bool same_scales(const motion_work_v1::Scales &h, const pf_v7::WorkScales &n) {
  return bits(h.kinetic_gross_J) == bits(n.kinetic_gross_J) && bits(h.kinetic_content_J) == bits(n.kinetic_content_J) &&
         bits(h.field_work_gross_J) == bits(n.field_work_gross_J) && bits(h.charge_abs_C) == bits(n.charge_abs_C) &&
         bits(h.field_bound_V_m) == bits(n.field_bound_V_m) && bits(h.tolerance_J) == bits(n.tolerance_J) &&
         h.status == n.status;
}
pf_v7::WorkScales native_scales(std::span<const pf_v7::Carrier> census, std::span<const pf_v7::Carrier> ends,
                                std::span<const pf_v7::Sink> sinks, double dt, std::span<const double> E,
                                std::span<const double> B, std::span<const double> J) {
  pf_v7::WorkCarrierSums c;
  std::size_t j = 0, s = 0;
  for (const auto &x : census) {
    const double *u = nullptr;
    if (s < sinks.size() && sinks[s].id == x.id)
      u = sinks[s++].incident_proper_velocity_m_s;
    else
      u = ends[j++].proper_velocity_m_s;
    const double K0 = pf_v7::kinetic_energy_J(x.mass_kg, x.proper_velocity_m_s);
    c.add(x.weight, x.charge_C, K0, pf_v7::kinetic_energy_J(x.mass_kg, u) - K0);
  }
  pf_v7::HallAxialMesh mesh7;
  if (!pf_v7::hall_axial_mesh_of(E.size(), B.size(), mesh7))
    throw std::runtime_error("mesh");
  pf_v7::WorkScales w;
  pf_v7::work_field_scales(w, dt, E, J, mesh7);
  pf_v7::work_finish(w, c);
  return w;
}
J scales_json(const pf_v7::WorkScales &w) {
  return {{"G_K_J", w.kinetic_gross_J}, {"C_K_J", w.kinetic_content_J}, {"G_W_J", w.field_work_gross_J},
          {"Q_C", w.charge_abs_C},      {"E_hat_V_m", w.field_bound_V_m}, {"tolerance_J", w.tolerance_J},
          {"status", w.status}};
}
} // namespace

int main(int argc, char **argv) {
  try {
    if (argc != 2)
      throw std::runtime_error("usage: test_v7_work_gross OUTPUT_DIRECTORY");
    const std::string dir = argv[1];
    v71::install_from_environment();
    struct ExecutionReport {
      std::string path;
      ~ExecutionReport() { v71::write_report(path); }
    } execution_report{dir + "/execution.json"};
    J cases = J::array();
    {
      // Hall axial mesh: native's derivation from the array lengths and its
      // edge endpoints equal the Hall mesh.
      const auto &g = mesh();
      pf_v7::HallAxialMesh h7;
      bool ok = pf_v7::hall_axial_mesh_of(std::uint64_t(g.ne), std::uint64_t(g.nb), h7) && h7.n == 45 &&
                h7.nz == std::uint64_t(g.nz) && h7.nn() == std::uint64_t(g.nn);
      for (int k = 0; ok && k <= g.nz; ++k)
        for (int j = 0; ok && j <= g.ny; ++j)
          for (int i = 0; ok && i < g.nx; ++i) {
            std::uint64_t a, b;
            h7.edge_nodes(std::uint64_t(g.ex(i, j, k)), a, b);
            ok = a == std::uint64_t(g.node(i, j, k)) && b == std::uint64_t(g.node(i + 1, j, k));
            if (j < g.ny) {
              h7.edge_nodes(std::uint64_t(g.ey(i, j, k)), a, b);
              ok = ok && a == std::uint64_t(g.node(i, j, k)) && b == std::uint64_t(g.node(i, j + 1, k));
            }
            if (k < g.nz) {
              h7.edge_nodes(std::uint64_t(g.ez(i, j, k)), a, b);
              ok = ok && a == std::uint64_t(g.node(i, j, k)) && b == std::uint64_t(g.node(i, j, k + 1));
            }
          }
      check(ok, "Hall axial mesh (n = 45, nz = 6) derived from the E/B lengths; edge endpoints equal Mesh::ex/ey/ez");
      check(pf_v7::total_field_sha256(field(false).E, field(false).B) ==
                hall_driver::total_field_sha256(field(false).E, field(false).B),
            "total-midpoint-E-B-binary64-v1 digest equals Hall total_field_sha256");
    }
    for (double weight : {3.45e10, 8e10})
      for (bool split : {false, true}) {
        const std::string name = std::string(split ? "split" : "uniform") + "-w" + (weight == 3.45e10 ? "3.45e10" : "8e10");
        const auto f = field(split);
        auto a = m::PicArena::create_full_time(options(), blob_states(), channels(), blob(weight, 3, 400), [&] {
          auto g = genesis();
          g.mesh_sha256 = hall_driver::motion_mesh_sha256(mesh(), options());
          return g;
        }());
        StepState st;
        const auto k = a->full_time_boundary()->accepted_step;
        auto b = binding_of(*a, 1e-14, sha(step_text("work-clock", k)));
        auto outer = a->begin_full_time_step(b, 1024, 32);
        m::ExposureBinding e1{b, m::FullTimeSlot::before_transport, b.begin_s, .5 * b.dt_s, 0, b.accepted_field_sha256};
        auto s1 = outer.first_exposure(e1);
        auto t1 = s1.prepare_compound_with_completion(1, [&](m::CollisionScope &) {
          return completion(e1, st.rng, sha("work-rng"));
        });
        if (!s1.try_publish(t1))
          throw std::runtime_error("first exposure");
        auto trial = outer.fork_transport_trial();
        const auto census = trial.native_census(1024);
        hall_driver::MotionBinding mb;
        mb.phase = hall_driver::ParticlePhase::full_time_position_and_proper_velocity;
        mb.accepted_step = b.accepted_step;
        mb.staged_arena_revision = trial.staged_arena_revision();
        mb.begin_s = b.begin_s;
        mb.end_s = b.end_s;
        mb.dt_s = b.dt_s;
        mb.staged_arena_sha256 = trial.state_prefix();
        mb.trial_field_sha256 = hall_driver::total_field_sha256(f.E, f.B);
        mb.mesh_sha256 = b.mesh_sha256;
        mb.material_sha256 = b.material_sha256;
        mb.species_sha256 = sha("dense-species");
        mb.common_clock_sha256 = b.common_clock_sha256;
        mb.kinetic_policy_sha256 = hall_driver::kinetic_policy_sha256(hall_driver::declared_kinetic_policy(options()));
        std::unique_ptr<hall_driver::MotionReceipt> receipt;
        std::string provider_refusal;
        try {
          receipt = std::make_unique<hall_driver::MotionReceipt>(
              hall_driver::prepare_canonical_motion(mesh(), options(), census, mb, {}, f.E, f.B));
        } catch (const std::exception &e) {
          provider_refusal = e.what();
        }
        J row = {{"case", name}, {"weight", weight}, {"split_field", split}};
        if (!receipt) {
          // The sealed provider applies the v6 net guard internally.
          row["sealed_provider_refusal"] = provider_refusal;
          check(split, "uniform-field dense blob passes the sealed provider (" + name + ")");
          cases.push_back(row);
          continue;
        }
        const auto v7 = m::v7_adapter::from_hall_receipt(*receipt, trial, f.E, f.B);
        const auto &rec = v7->record;
        const auto &c7 = trial.census_v7();
        const std::span<const pf_v7::Carrier> rows{reinterpret_cast<const pf_v7::Carrier *>(c7.carriers.data()),
                                                   c7.carriers.size()};
        const std::span<const pf_v7::Carrier> ends{reinterpret_cast<const pf_v7::Carrier *>(rec.endpoints.data()),
                                                   rec.endpoints.size()};
        const double K = rec.budget[0], W = rec.budget[1];
        const auto hall0 = hall_scales(rows, ends, {}, b.dt_s, f.E, receipt->current_A_m(), K, W);
        const double net_tol = 1e-20 + 1e-10 * (std::abs(K) + std::abs(W));
        row["K_J"] = K;
        row["W_J"] = W;
        row["defect_J"] = rec.budget[2];
        row["defect_over_gross_tolerance"] = std::abs(rec.budget[2]) / hall0.tolerance_J;
        row["defect_over_v6_net_tolerance"] = std::abs(rec.budget[2]) / net_tol;
        row["net_over_gross"] = std::abs(K) / hall0.kinetic_gross_J;
        row["receipt_x87"] = {{"K_J", v7->receipt_kinetic_J}, {"W_J", v7->receipt_work_J},
                              {"K_ulps_from_exact", std::abs(std::int64_t(bits(K)) - std::int64_t(bits(v7->receipt_kinetic_J)))}};
        auto mbn = m::FullTimeMotionBinding{b, rec.record_sha256, b.accepted_field_sha256, sha("f-after"),
                                            sha("h-after"), a->full_time_boundary()->source_sha256};
        // Tampers: a resealed record whose arrays and budget are mutually
        // consistent except the stated inconsistency.
        struct Tamper {
          std::string name;
          std::vector<m::GapParticle> ends;
          std::vector<double> E, J;
          double K, W;
          bool kinetic_override = false, work_override = false;
        };
        auto base_tamper = [&](const std::string &n) {
          Tamper t{n, {rec.endpoints.begin(), rec.endpoints.end()}, f.E, receipt->current_A_m(), K, W};
          return t;
        };
        std::vector<Tamper> tampers;
        {
          // One electron: kinetic change off by 1e-3 of 1.5 kT (x w).
          auto t = base_tamper("one electron kinetic change off by 1e-3 of 1.5 kT w");
          for (auto &e : t.ends)
            if (e.charge_C < 0) {
              const double KE = m::pic_kinetic_energy_J(e), delta = 1e-3 * 1.5 * 3 * m::pic_e;
              const double s = std::sqrt(1 + delta / KE);
              for (auto &v : e.proper_velocity_m_s)
                v *= s;
              break;
            }
          tampers.push_back(std::move(t));
        }
        {
          auto t = base_tamper("work field 1e-6 (relative) away from the push field");
          for (auto &x : t.E)
            x *= 1 + 1e-6;
          tampers.push_back(std::move(t));
        }
        {
          // One carrier's mean gross exchange with the wrong sign: 2 G_W / N on
          // the edge with the largest |E J|.
          auto t = base_tamper("one carrier's mean gross exchange with the wrong sign");
          std::size_t best = 0;
          for (std::size_t e = 0; e < t.J.size(); ++e)
            if (std::abs(t.E[e] * t.J[e]) > std::abs(t.E[best] * t.J[best]))
              best = e;
          const double dW = 2 * hall0.field_work_gross_J / double(census.size());
          t.J[best] += dW / (b.dt_s * t.E[best]);
          tampers.push_back(std::move(t));
        }
        {
          auto t = base_tamper("one ulp of K (defect adjusted)");
          t.K = std::nextafter(K, 1e300);
          t.kinetic_override = true;
          tampers.push_back(std::move(t));
        }
        {
          auto t = base_tamper("one ulp of W (defect adjusted)");
          t.W = std::nextafter(W, 1e300);
          t.work_override = true;
          tampers.push_back(std::move(t));
        }
        J rows_t = J::array();
        for (auto &t : tampers) {
          const std::span<const pf_v7::Carrier> te{reinterpret_cast<const pf_v7::Carrier *>(t.ends.data()), t.ends.size()};
          if (!t.kinetic_override)
            t.K = pf_v7::kinetic_change_J(rows, te, {}).value;
          if (!t.work_override)
            t.W = pf_v7::field_work_J(b.dt_s, t.E, t.J).value;
          const auto hs = hall_scales(rows, te, {}, b.dt_s, t.E, t.J, t.K, t.W);
          const auto ns = native_scales(rows, te, {}, b.dt_s, t.E, f.B, t.J);
          const double defect = t.K - t.W, margin = std::abs(defect) / ns.tolerance_J;
          auto x = std::make_shared<std::tuple<Tamper, m::FullTimeMotionRecordV7>>(t, rec);
          auto &tt = std::get<0>(*x);
          auto &xr = std::get<1>(*x);
          xr.endpoints = tt.ends;
          xr.trial_E = tt.E;
          xr.current_A_m = tt.J;
          xr.trial_field_sha256 = pf_v7::total_field_sha256(tt.E, f.B);
          xr.current_A_m_sha256 = pf_v7::array_sha256(tt.J);
          xr.budget[0] = tt.K;
          xr.budget[1] = tt.W;
          xr.budget[2] = defect;
          xr.record_sha256 = m::full_time_motion_record_v7_sha256(xr);
          auto mbt = mbn;
          mbt.transport_trial_sha256 = xr.record_sha256;
          const auto msg = message_of([&] { (void)trial.prepare_total_field_motion(mbt, m::bind_total_field_motion_v7(xr, x)); });
          const bool refused = !msg.empty();
          const bool exact_refusal = t.kinetic_override || t.work_override;
          check(refused && (exact_refusal || margin > 1) && same_scales(hs, ns),
                "tamper refused (" + name + "): " + t.name);
          rows_t.push_back({{"case", t.name}, {"defect_J", defect}, {"defect_over_tolerance", margin},
                            {"tolerance_bits_equal_hall", same_scales(hs, ns)}, {"refusal", msg},
                            {"native_scales", scales_json(ns)}});
        }
        row["tampers"] = rows_t;
        // Admission of the unmodified record on the same trial (control).
        J audit;
        {
          auto moved = trial.prepare_total_field_motion(mbn, m::v7_adapter::admit(v7));
          audit = moved.audit().at("work_tolerance");
        }
        const auto n0 = native_scales(rows, ends, {}, b.dt_s, f.E, f.B, receipt->current_A_m());
        const bool control_bits = same_scales(hall0, n0) && bits(audit.at("tolerance_J").get<double>()) == bits(hall0.tolerance_J) &&
                                  bits(audit.at("E_hat_V_m").get<double>()) == bits(hall0.field_bound_V_m) &&
                                  bits(audit.at("G_K_J").get<double>()) == bits(hall0.kinetic_gross_J) &&
                                  bits(audit.at("C_K_J").get<double>()) == bits(hall0.kinetic_content_J) &&
                                  bits(audit.at("G_W_J").get<double>()) == bits(hall0.field_work_gross_J) &&
                                  bits(audit.at("Q_C").get<double>()) == bits(hall0.charge_abs_C);
        check(control_bits, "native admits the dense blob (" + name + ") with tolerance and scales bitwise equal to Hall's");
        row["control"] = {{"admitted", true}, {"native_scales", audit}, {"tolerance_bits_equal_hall", control_bits}};
        cases.push_back(row);
      }
    J output = {{"status", "passed"},
                {"check_count", checks},
                {"checks", results},
                {"rule", motion_work_v1::rule},
                {"hall_reference_header", HALL_V32_WORK_GROSS_HEADER},
                {"cases", cases}};
    std::ofstream(dir + "/v7-work-gross.json") << output.dump(2) << '\n';
    std::cout << checks << " motion-work-gross-v1 controls passed\n";
    return 0;
  } catch (const std::exception &x) {
    std::cerr << "FAILED after " << checks << " checks: " << x.what() << '\n';
    return 1;
  }
}
