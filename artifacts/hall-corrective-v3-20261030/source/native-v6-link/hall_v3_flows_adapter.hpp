#pragma once
// Hall corrective v3 MotionReceipt -> native v6 typed binary motion admission
// (native-full-time-candidate006, API-006.txt section 4). The native v6 Hall adapter
// (mfem_pic_hall_motion_v3_adapter.cpp) maps only the plain v3 record of a Hall
// corrective v2 receipt; this adapter maps the Hall v3 receipt to the plain record
// (no newborns, no sinks: v2-identical receipt) or to the flows record (schema
// canonical-full-time-total-field-motion-v3-flows-binary, coverage
// all_native_carriers_births_and_first_hit_sinks, sinks and births sections) and binds
// it with raikou::mfem_backend::bind_total_field_motion_v3. Compile against the native
// v6 headers and Hall v3 source/motion_receipt.hpp; no JSON value is built.
#include "raikou/mfem_pic_full_time.h"
#include "raikou/mfem_pic_canonical.h"
#include "motion_receipt.hpp"
#include <memory>
#include <string>
#include <vector>

namespace hall_v3_native_v6 {
namespace m = raikou::mfem_backend;

// Owned arrays of one record (the record's spans point into them).
struct FlowsRecord {
    std::vector<m::GapParticle> endpoints;
    std::vector<m::FullTimeMotionPath> paths;
    std::vector<m::FullTimeSinkRecord> sinks;
    std::vector<m::FullTimeBirthFlight> births;
    m::FullTimeMotionRecordV3 record;
    void bind() {
        record.endpoints = endpoints;
        record.paths = paths;
        record.sinks = sinks;
        record.births = births;
    }
};
// motion-array-binary64-v1 digest (the Hall array_hash encoding).
inline std::string array_digest(const std::vector<double>& v) {
    m::canonical::Hasher h;
    h.w.tag("motion-array-binary64-v1");
    h.w.u64(v.size());
    for (double x : v) {
        h.w.f64(x);
        h.flush_if_large();
    }
    return m::canonical::hex(h.finish());
}

inline std::shared_ptr<FlowsRecord> flows_record(const hall_driver::MotionReceipt& r) {
    auto o = std::make_shared<FlowsRecord>();
    const bool flows = !r.v2_identical();
    o->endpoints.assign(r.endpoints().begin(), r.endpoints().end());
    o->paths.reserve(r.paths().size());
    for (const auto& p : r.paths()) o->paths.push_back({p.native_id, p.field_start_m, p.field_displacement_m, p.flight_s});
    o->sinks.reserve(r.sinks().size());
    for (const auto& k : r.sinks()) {
        m::FullTimeSinkRecord s;
        s.native_id = k.native_id;
        s.kind = k.kind == 1 ? m::FullTimeSinkKind::collected_at_conductor : m::FullTimeSinkKind::escaped_exterior;
        s.surface = static_cast<m::FullTimeSinkSurface>(k.boundary);
        s.axis = std::uint32_t(k.axis);
        s.plane_m = k.plane_m;
        s.fraction = k.fraction;
        s.segment_m = k.segment_m;
        s.hit_point_m = k.hit_point_m;
        s.incident_proper_velocity_m_s = k.incident_proper_velocity_m_s;
        s.flight_s = k.flight_s;
        s.discarded_flight_s = k.discarded_flight_s;
        o->sinks.push_back(s);
    }
    o->births.reserve(r.births().size());
    for (const auto& b : r.births()) o->births.push_back({b.native_id, b.birth_time_s, b.flight_s});
    auto& x = o->record;
    const auto& b = r.binding();
    x.schema = flows ? m::motion_v3_flows_schema : m::motion_v3_schema;
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
    const auto& u = r.budget();
    x.budget = {u.kinetic_change_J, u.field_work_J, u.work_defect_J, u.continuity_C, u.plasma_dt, u.maximum_gyro, u.maximum_displacement_cells};
    x.census_sha256 = r.census_sha256();
    x.coverage = flows ? "all_native_carriers_births_and_first_hit_sinks" : "all_native_carriers_no_sources_no_boundaries";
    x.current_A_m_sha256 = array_digest(r.current_A_m());
    x.endpoint_charge_C_sha256 = array_digest(r.endpoint_charge_C());
    x.initial_charge_C_sha256 = array_digest(r.initial_charge_C());
    const auto& k = r.kinetic_policy();
    x.kinetic_policy = {k.width_m, k.gap_m, k.short_z_m, k.feed_z_m, k.lower.left_m, k.lower.right_m, k.lower.bottom_m, k.lower.top_m,
                        k.upper.left_m, k.upper.right_m, k.upper.bottom_m, k.upper.top_m};
    x.native_to_field_m = r.frame().native_to_field_m;
    x.neutral_motion_owner = "same_full_time_q0_drift";
    o->bind();
    x.record_sha256 = r.record_v3_sha256();   // the Hall producer's digest; native recomputes it
    return o;
}
// Typed binary admission of a Hall v3 receipt (the receipt need not outlive it).
inline m::FullTimeMotionAdmission bind_hall_v3_motion(const hall_driver::MotionReceipt& r) {
    auto x = flows_record(r);
    return m::bind_total_field_motion_v3(x->record, x);
}
}  // namespace hall_v3_native_v6
