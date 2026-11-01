// hall_receipt_v7.hpp -- transitional host adapter: a Hall MotionReceipt that
// does not emit the native v7 forms (sealed 8dddec0 provider, Hall corrective
// v2, v3, v3.1) -> motion record v7 (API-007.txt). Header only; compiled
// against the Hall motion_receipt.hpp of the provider in use.
//
// The adapter recomputes the two budget sums in the v7 form with the CPU
// reference (pf_v7::kinetic_change_J over the trial census and the receipt's
// end states; pf_v7::field_work_J over the caller's midpoint E and the
// receipt's current_A_m), the work defect K - W and the motion-work-gross-v1
// guard (API-007.txt 1.1 amendment), the array digests, the four section roots
// and the record digest; the census digest is the trial's census v2; the
// trial field (E, B) and the current enter the record as admission inputs.
// Everything else is copied from the receipt.
// Cost: O(carriers) host work (it exists for parity runs and for a Hall that
// does not yet emit the v7 forms; HALL-LINK.txt describes the device path).
#pragma once
#include "pf_v7/forms.hpp"
#include "raikou/mfem_pic_full_time.h"
#include <memory>
#include <stdexcept>
#include <vector>

namespace raikou::mfem_backend::v7_adapter {
struct MotionV7 {
  std::vector<FullTimeMotionPath> paths;
  std::vector<FullTimeSinkRecord> sinks;
  std::vector<FullTimeBirthFlight> births;
  FullTimeMotionRecordV7 record; // spans: endpoints -> receipt, others -> above
  pf_v7::SumValue kinetic, work; // the recomputed v7 sums
  pf_v7::WorkScales scales;      // motion-work-gross-v1
  double receipt_kinetic_J = 0, receipt_work_J = 0; // the receipt's (v6-form) values
};
// Builds the record. E_mid, B_mid: the total midpoint field (edge E, face B
// values of the Hall axial mesh) the receipt was computed with (its
// trial_field_sha256); neutral_dt_max_s: the native configuration value. The
// receipt, E_mid and B_mid must outlive the admission (referenced, not copied).
template <class Receipt>
std::shared_ptr<MotionV7> from_hall_receipt(const Receipt &r, const PrivateTransportTrial &trial,
                                            std::span<const double> E_mid, std::span<const double> B_mid,
                                            double neutral_dt_max_s = 1e-14) {
  auto m = std::make_shared<MotionV7>();
  const auto paths = r.paths();
  m->paths.resize(paths.size());
  for (std::size_t i = 0; i < paths.size(); ++i)
    m->paths[i] = {paths[i].native_id, paths[i].field_start_m, paths[i].field_displacement_m, paths[i].flight_s};
  std::vector<pf_v7::Sink> canonical_sinks;
  std::vector<std::pair<std::uint64_t, const double *>> sink_u; // id -> incident u
  if constexpr (requires { r.sinks(); }) {
    for (const auto &k : r.sinks()) {
      FullTimeSinkRecord s;
      s.native_id = k.native_id;
      s.kind = FullTimeSinkKind(k.kind);
      s.surface = FullTimeSinkSurface(k.boundary);
      s.axis = std::uint32_t(k.axis);
      s.plane_m = k.plane_m;
      s.fraction = k.fraction;
      s.segment_m = k.segment_m;
      s.hit_point_m = k.hit_point_m;
      s.incident_proper_velocity_m_s = k.incident_proper_velocity_m_s;
      s.flight_s = k.flight_s;
      s.discarded_flight_s = k.discarded_flight_s;
      m->sinks.push_back(s);
      pf_v7::Sink c{};
      c.id = k.native_id;
      for (int d = 0; d < 3; ++d)
        c.incident_proper_velocity_m_s[d] = k.incident_proper_velocity_m_s[d];
      canonical_sinks.push_back(c);
    }
    for (const auto &b : r.births())
      m->births.push_back({b.native_id, b.birth_time_s, b.flight_s});
  }
  const auto &census = trial.census_v7();
  const std::span<const pf_v7::Carrier> rows{reinterpret_cast<const pf_v7::Carrier *>(census.carriers.data()),
                                             census.carriers.size()};
  const auto endpoints = r.endpoints();
  const std::span<const pf_v7::Carrier> ends{reinterpret_cast<const pf_v7::Carrier *>(endpoints.data()),
                                             endpoints.size()};
  const auto &b = r.binding();
  m->kinetic = pf_v7::kinetic_change_J(rows, ends, canonical_sinks);
  m->work = pf_v7::field_work_J(b.dt_s, E_mid, r.current_A_m());
  if (!m->kinetic.ok() || !m->work.ok())
    throw std::invalid_argument("motion v7 adapter: kinetic/work sum not representable or census mismatch");
  const auto &u = r.budget();
  m->receipt_kinetic_J = u.kinetic_change_J;
  m->receipt_work_J = u.field_work_J;
  const double K = m->kinetic.value, W = m->work.value, defect = K - W;
  {
    // motion-work-gross-v1 scales over the same (census, end state) pairs.
    pf_v7::WorkCarrierSums sums;
    std::size_t j = 0, k = 0;
    for (const auto &c : rows) {
      const double *u_end = nullptr;
      if (k < canonical_sinks.size() && canonical_sinks[k].id == c.id)
        u_end = canonical_sinks[k++].incident_proper_velocity_m_s;
      else if (j < ends.size() && ends[j].id == c.id)
        u_end = ends[j++].proper_velocity_m_s;
      else
        throw std::invalid_argument("motion v7 adapter: census carrier without end state");
      const double K0 = pf_v7::kinetic_energy_J(c.mass_kg, c.proper_velocity_m_s);
      sums.add(c.weight, c.charge_C, K0, pf_v7::kinetic_energy_J(c.mass_kg, u_end) - K0);
    }
    pf_v7::HallAxialMesh mesh;
    if (!pf_v7::hall_axial_mesh_of(E_mid.size(), B_mid.size(), mesh))
      throw std::invalid_argument("motion v7 adapter: E/B are not one Hall axial mesh");
    pf_v7::work_field_scales(m->scales, b.dt_s, E_mid, r.current_A_m(), mesh);
    pf_v7::work_finish(m->scales, sums);
    if (!pf_v7::work_admits(m->scales, defect))
      throw std::invalid_argument("canonical motion work mismatch");
  }
  auto &x = m->record;
  x.accepted_step = b.accepted_step;
  x.staged_arena_revision = b.staged_arena_revision;
  x.trial_stage = b.trial_stage;
  x.begin_s = b.begin_s;
  x.end_s = b.end_s;
  x.dt_s = b.dt_s;
  x.staged_arena_sha256 = b.staged_arena_sha256;
  x.trial_field_sha256 = b.trial_field_sha256;
  x.mesh_sha256 = b.mesh_sha256;
  x.material_sha256 = b.material_sha256;
  x.species_sha256 = b.species_sha256;
  x.common_clock_sha256 = b.common_clock_sha256;
  x.kinetic_policy_sha256 = b.kinetic_policy_sha256;
  x.budget = {K, W, defect, u.continuity_C, u.plasma_dt, u.maximum_gyro, u.maximum_displacement_cells};
  x.census_sha256 = trial.census_sha256();
  x.current_A_m_sha256 = pf_v7::array_sha256(r.current_A_m());
  x.endpoint_charge_C_sha256 = pf_v7::array_sha256(r.endpoint_charge_C());
  x.initial_charge_C_sha256 = pf_v7::array_sha256(r.initial_charge_C());
  x.endpoints = endpoints;
  const auto &k = r.kinetic_policy();
  x.kinetic_policy = {k.width_m,        k.gap_m,          k.short_z_m,     k.feed_z_m,
                      k.lower.left_m,   k.lower.right_m,  k.lower.bottom_m, k.lower.top_m,
                      k.upper.left_m,   k.upper.right_m,  k.upper.bottom_m, k.upper.top_m};
  x.native_to_field_m = r.frame().native_to_field_m;
  x.neutral_stream_dt_max_s = neutral_dt_max_s;
  x.paths = m->paths;
  x.sinks = m->sinks;
  x.births = m->births;
  x.trial_E = E_mid;
  x.trial_B = B_mid;
  x.current_A_m = r.current_A_m();
  x.record_sha256 = full_time_motion_record_v7_sha256(x);
  return m;
}
// Admission of an adapted receipt (the storage keeps the converted arrays).
inline FullTimeMotionAdmission admit(const std::shared_ptr<MotionV7> &m) {
  return bind_total_field_motion_v7(m->record, m);
}
} // namespace raikou::mfem_backend::v7_adapter
