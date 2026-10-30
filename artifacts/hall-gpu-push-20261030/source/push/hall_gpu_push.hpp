// hall_gpu_push.hpp -- GPU backend of the corrected Hall motion path.
//
// MotionEngine reproduces, bit for bit, the numerical content of
//   hall_driver::prepare_canonical_motion   (endpoints, canonical paths, current_A_m,
//                                            initial_charge_C, endpoint_charge_C, budget)
//   hall_corrective::charge_of_endpoints    (net words and gross of the endpoint census)
// for carriers given in strictly increasing native-ID order. Binding, census hash,
// Count2 packet identity and record sealing are host metadata handled by the caller
// (see hall_adapter.hpp); this engine owns the arithmetic, its refusal order and
// texts, and the transfers.
//
// Residency: set_carriers() and set_field() stage into pinned memory and enqueue
// asynchronous uploads on the session stream; carriers and E/B stay on the device
// until replaced. run() enqueues the kernels and the downloads of exactly the record
// arrays, then synchronizes once (plus one early wait for the edge-record count).
// Returned spans point into engine-owned pinned memory and stay valid until the next
// set_carriers()/set_field()/run() call or destruction.
#pragma once
#include "../common/hall_gpu_common.hpp"
#include "../common/hall_gpu_primitives.hpp"
#include "push_core.hpp"

#include <limits>
#include <span>
#include <string>
#include <vector>

namespace hall_gpu::push {

struct MeshTables {
    hgp::Mesh mesh{};
    std::vector<double> VN;            // nn dual node volumes (plasma guard)
    std::vector<int> edge_a, edge_b;   // ne edge endpoints (continuity)
    std::vector<double> edge_length;   // ne
};
struct Geometry {
    hgp::Policy policy{};
    double native_to_field[3] = {0, 0, 0};
};
struct CarrierView {
    std::span<const std::uint64_t> id;
    std::span<const double> x, y, z, ux, uy, uz, charge, mass, weight;
    std::size_t size() const noexcept { return id.size(); }
};
struct Budget {
    double kinetic_change_J = 0, field_work_J = 0, work_defect_J = 0;
    double continuity_C = 0, plasma_dt = 0, maximum_gyro = 0, maximum_displacement_cells = 0;
};
struct Refusal {
    bool refused = false;
    std::string message;   // identical to the reference exception text
    std::string stage;     // "field", "census", "plasma", "motion", "work", "continuity", "charge", "charge_guard"
    long rank = -1;        // carrier rank (ID order) when carrier-specific
    std::uint64_t id = 0;
    int code = 0;
};
struct Phase { std::string name; double ms = 0; };

struct MotionView {
    std::size_t n = 0;
    Refusal refusal;          // prepare_canonical_motion
    Refusal charge_refusal;   // charge_of_endpoints (evaluated only when the motion succeeded)
    std::span<const double> end_x, end_y, end_z, end_ux, end_uy, end_uz;     // native frame
    std::span<const double> start_x, start_y, start_z, disp_x, disp_y, disp_z;   // canonical field path
    std::span<const double> dke;   // per carrier KE(end) - KE(start), binary64
    std::span<const double> current_A_m, initial_charge_C, endpoint_charge_C, gross_C;
    std::span<const double> free_net_hi, free_net_lo, free_gross_C;
    Budget budget;
    std::uint32_t edge_records = 0;
    bool diagnostics = false;
    std::span<const int> code1, code2, code3;      // diagnostics only
    std::span<const double> diag_real;              // [hgp::diag_real_fields][n]
    std::span<const int> diag_int;                  // [hgp::diag_int_fields][n]
    std::vector<Phase> gpu_phases;                  // stream event timings (ms)
    double host_fold_ms = 0, run_wall_ms = 0;
};

struct EngineOptions {
    std::size_t capacity = 0;               // maximum carriers
    bool diagnostics_capable = false;       // allocate per-carrier diagnostic outputs
    double initial_edge_records_per_carrier = 30;
};

inline const char* code_message(int code) {
    switch (code) {
        case hgp::invalid_carrier: return "motion invalid carrier";
        case hgp::nonfinite_kinematics: return "motion nonfinite kinematics";
        case hgp::interior_boundary_owner: return "particle requires native boundary owner";
        case hgp::interior_metal: return "particle in/on native metal";
        case hgp::unresolved_plasma_operand: return "motion unresolved plasma operand";
        case hgp::charge_outside_mesh: return "axial charge outside mesh";
        case hgp::neutral_identity: return "neutral stream identity/species";
        case hgp::neutral_momentum_finite: return "neutral momentum finite";
        case hgp::neutral_gamma: return "neutral momentum outside resolved range";
        case hgp::kinetic_displacement_nonfinite: return "nonfinite kinetic displacement";
        case hgp::neutral_incident: return "neutral incident kinetic operand";
        case hgp::neutral_first_hit: return "neutral first hit requires native sink owner";
        case hgp::neutral_displacement: return "neutral full-interval displacement guard";
        case hgp::excessive_path: return "excessive axial path";
        case hgp::path_outside_mesh: return "axial path outside mesh";
        case hgp::boris_input: return "Boris update requires finite data and positive dt";
        case hgp::boris_overflow: return "Boris update overflow";
        case hgp::not_converged: return "axial native work path did not converge";
        case hgp::gyro_displacement_guard: return "axial full-interval gyro/displacement guard";
        case hgp::lost_surface: return "lost axial surface intersection";
        case hgp::boundary_owner: return "motion needs boundary owner before collection";
        case hgp::charged_first_hit: return "charged first hit requires native sink owner";
        case hgp::rounded_endpoint_crosses: return "rounded canonical endpoint crosses field boundary";
        case hgp::unresolved_macrocharge: return "motion unresolved macrocharge";
        case hgp::canonical_charge_operand: return "invalid canonical charge operand";
        case hgp::unresolved_canonical_charge: return "unresolved canonical charge";
        case hgp::deposited_words: return "invalid deposited charge words";
        case hgp::nonfinite_charge_position: return "nonfinite charge position";
        case hgp::canonical_charge_outside_mesh: return "canonical charge outside mesh";
        case hgp::gpu_path_capacity: return "device path cut capacity exceeded";
        case hgp::gpu_frame_mismatch: return "device endpoint frame mismatch between motion and free-charge census";
        default: return "unknown motion status";
    }
}

class MotionEngine {
public:
    MotionEngine(Session&, const MeshTables&, const Geometry&, const EngineOptions&);
    ~MotionEngine();
    MotionEngine(const MotionEngine&) = delete;
    MotionEngine& operator=(const MotionEngine&) = delete;

    const CompileRecord& compile_record() const noexcept;
    void set_carriers(const CarrierView&);
    void set_field(std::span<const double> E, std::span<const double> B);
    MotionView run(double dt, bool diagnostics);
    std::size_t device_bytes() const noexcept;

private:
    struct Impl;
    Impl* p_;
};

}  // namespace hall_gpu::push
