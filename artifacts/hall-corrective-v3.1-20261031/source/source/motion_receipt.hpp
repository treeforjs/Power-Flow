#pragma once
// Hall motion receipt, v3 (v2 plus kinetic sink and surface birth records).
// v2 description (unchanged): v2 equivalent of the sealed hall-canonical-driver-20261029/
// motion_receipt.hpp (commit 8dddec0). Every type, member, accessor and free
// function below has the sealed name, layout and meaning, so native003 and
// native v5 (whose adapters only forward-declare hall_driver::MotionReceipt and
// call hall_driver::motion_record) admit it unchanged. The one difference is
// construction: besides the CPU provider prepare_canonical_motion, the v2
// assembler hall_corrective_v3::MotionAssembly is a friend, so a device-computed
// motion (GPU push) yields the identical immutable receipt and record. The
// receipt additionally carries the binary "motion v3" record digest (same
// content, canonical binary encoding, computed for every fixed-point iterate);
// the v2 JSON-text seal required by native admission is computed only for the
// accepted motion (MotionAssembly::seal_v2). Native adapters never access the
// receipt layout (they only forward-declare it and call motion_record).
//
// A program links either this header with source/motion_provider.cpp or the
// sealed header with the sealed native_arena_motion.cpp, never both.
//
// v3: the receipt additionally carries the kinetic sink records (carriers whose
// path first hits the declared kinetic boundary within the step, ascending
// native ID; their truncated paths stay in paths(), they are absent from
// endpoints()) and the surface birth records (newborn carriers with a conductor
// surface origin, ascending native ID above every census ID). Without sinks and
// births the receipt, its v2 JSON record, its v2 seal and its v3 binary digest
// are identical to Hall v2. With sinks or births the JSON record has schema
// canonical-full-time-total-field-motion-v3-sinks (README: record format), which
// native candidate003/native v5 do not recognize (they refuse it).
#include "raikou/mfem_pic_arena.h"
#include "hall-canonical-driver-20261029/kinetic_boundary.hpp"
#include <span>
#include <string>
#include <vector>

namespace axial {struct Mesh;}
namespace hall_corrective_v3 {
class MotionAssembly;
struct VerifiedField;   // field_contract.hpp (v3.1)
namespace native=raikou::mfem_backend;
// One kinetic sink. The first fields are the native v6 FullTimeSinkRecord
// (native-full-time-candidate006 API-006.txt section 2) with the same meaning
// and values; the remaining fields are the Hall side (bank point, energies).
struct SinkRecord {
    std::uint64_t native_id=0;
    int kind=0;                          // 1 collected_at_conductor, 2 escaped_exterior
    int boundary=-1;                     // surface = hall_driver::KineticSink: 0 lower foil, 1 upper foil, 2 outer xy, 3 short bridge, 4 axial
    int axis=-1;double plane_m=0;        // hit axis and kinetic plane (native frame)
    double fraction=0;                   // first-hit fraction of the straight segment
    native::PicVector segment_m{};       // full native displacement of the flight F = dt (represented segment)
    native::PicVector hit_point_m{};     // start + fraction*segment componentwise, [axis] = plane
    native::PicVector incident_proper_velocity_m_s{};
    double flight_s=0,discarded_flight_s=0;   // F*fraction, F - flight_s
    // Hall side
    std::string packet_id,state;
    bool collected=false;int conductor=-1;    // collected: 0 lower foil, 1 upper foil, 2 short bridge
    native::PicVector field_point_m{};   // hit point in the field frame (bank CIC deposit point)
    double charge_C=0,mass_kg=0,weight=0;
    double kinetic_J=0;                  // per carrier, pic_kinetic_energy_J at the incident velocity
    double weighted_kinetic_J=0;         // fl(weight*kinetic_J): wall heat (foils), short-bridge or escaped energy
    double initial_weighted_kinetic_J=0;
    double gyro=0,displacement_cells=0;  // full-interval guard values of the represented step
};
// Newborn carrier of the step (native v6 trial.birth_flights(): native ID, birth
// time tau; flight end - tau). The newborn is in the census at its birth position,
// which may be a point on one face of one foil with the proper velocity pointing
// into vacuum (source v5 neutral packets, emitted electrons). A charged newborn
// takes its macro-charge from the conductor bank of that foil at the birth point.
struct Newborn {std::uint64_t native_id=0;double birth_time_s=0;};
struct NewbornRecord {
    std::uint64_t native_id=0;double birth_time_s=0,flight_s=0;   // flight = end - tau
    int origin=-1;int axis=-1;            // 0 lower-foil face, 1 upper-foil face, -1 interior point; face normal axis
    native::PicVector start_m{};std::string packet_id,state;
    double charge_C=0,weight=0,birth_weighted_kinetic_J=0;   // fl(w KE) at birth (energy supplied by the source)
};
// Binary record: the native v6 flows schema (API-006.txt section 4) for a step with
// newborns or first-hit sinks; the plain v3 schema otherwise.
inline constexpr char motion_v3_flows_schema[]="canonical-full-time-total-field-motion-v3-flows-binary";
inline constexpr char flows_coverage[]="all_native_carriers_births_and_first_hit_sinks";
// JSON export of a receipt with sinks/newborns (diagnostic; native v6 admits the binary record only).
inline constexpr char motion_api_sinks[]="hall-v3-total-field-motion-flows-json";
inline constexpr char coverage_sinks[]="all_native_carriers_births_and_first_hit_sinks";
inline constexpr char sink_record_api[]="native-v6-full-time-sink-record+hall-v3-bank";
}
namespace hall_driver {
namespace native=raikou::mfem_backend;
inline constexpr char motion_api[]="canonical-full-time-total-field-motion-v2";
enum class ParticlePhase {unbound,full_time_position_and_proper_velocity};
struct MotionBinding {
    ParticlePhase phase=ParticlePhase::unbound;
    std::uint64_t accepted_step=0,staged_arena_revision=0,trial_stage=0;
    double begin_s=0,end_s=0,dt_s=0;
    std::string staged_arena_sha256,trial_field_sha256,mesh_sha256;
    std::string material_sha256,species_sha256,common_clock_sha256,kinetic_policy_sha256;
};
struct TranslationFrame {native::PicVector native_to_field_m{};};
struct CarrierPath {
    std::uint64_t native_id=0;
    native::PicVector field_start_m{},field_displacement_m{};
    double flight_s=0;
};
struct MotionBudget {
    double kinetic_change_J=0,field_work_J=0,work_defect_J=0;
    double continuity_C=0,plasma_dt=0,maximum_gyro=0,maximum_displacement_cells=0;
};
class MotionReceipt;
MotionReceipt prepare_canonical_motion(const axial::Mesh&,const native::GapPicOptions&,
    std::span<const native::ArenaRateEntry>,const MotionBinding&,
    std::span<const std::uint64_t> active_layer_ids,
    const std::vector<double>& total_midpoint_E,const std::vector<double>& total_midpoint_B);

class MotionReceipt {
    friend MotionReceipt prepare_canonical_motion(const axial::Mesh&,
        const native::GapPicOptions&,std::span<const native::ArenaRateEntry>,
        const MotionBinding&,std::span<const std::uint64_t>,
        const std::vector<double>&,const std::vector<double>&);
    friend class ::hall_corrective_v3::MotionAssembly;
    MotionBinding binding_;TranslationFrame frame_;MotionBudget budget_;KineticPolicy kinetic_policy_;
    std::string census_sha256_,seal_sha256_,record_v3_sha256_;
    std::vector<native::GapParticle> endpoints_;
    std::vector<CarrierPath> paths_;
    std::vector<double> current_A_m_,initial_charge_C_,endpoint_charge_C_;
    // v3 members after every v2 member (the v2 member layout is unchanged).
    std::vector<hall_corrective_v3::SinkRecord> sinks_;std::vector<hall_corrective_v3::NewbornRecord> births_;
    bool birth_scope_=false;   // the step carries a native birth scope (possibly empty): flows record
    MotionReceipt()=default;
public:
    MotionReceipt(const MotionReceipt&)=delete;
    MotionReceipt&operator=(const MotionReceipt&)=delete;
    MotionReceipt(MotionReceipt&&)noexcept=default;
    MotionReceipt&operator=(MotionReceipt&&)noexcept=default;
    const MotionBinding& binding()const noexcept{return binding_;}
    const TranslationFrame& frame()const noexcept{return frame_;}
    const MotionBudget& budget()const noexcept{return budget_;}
    const KineticPolicy& kinetic_policy()const noexcept{return kinetic_policy_;}
    const std::string& census_sha256()const noexcept{return census_sha256_;}
    // v2 JSON-text seal (native admission); empty until MotionAssembly::seal_v2.
    const std::string& seal_sha256()const noexcept{return seal_sha256_;}
    // Binary motion v3 record digest (always present).
    const std::string& record_v3_sha256()const noexcept{return record_v3_sha256_;}
    std::span<const native::GapParticle> endpoints()const noexcept{return endpoints_;}
    std::span<const CarrierPath> paths()const noexcept{return paths_;}
    const std::vector<double>& current_A_m()const noexcept{return current_A_m_;}
    const std::vector<double>& initial_charge_C()const noexcept{return initial_charge_C_;}
    const std::vector<double>& endpoint_charge_C()const noexcept{return endpoint_charge_C_;}
    // v3 records (empty in a v2-identical receipt).
    std::span<const hall_corrective_v3::SinkRecord> sinks()const noexcept{return sinks_;}
    std::span<const hall_corrective_v3::NewbornRecord> births()const noexcept{return births_;}
    // v2-identical: no newborns, no sinks and no declared birth scope (plain v3 record).
    bool v2_identical()const noexcept{return !birth_scope_&&sinks_.empty()&&births_.empty();}
    bool birth_scope()const noexcept{return birth_scope_;}
};
std::string canonical_census_sha256(std::span<const native::ArenaRateEntry>);
std::string motion_mesh_sha256(const axial::Mesh&,const native::GapPicOptions&);
std::string total_field_sha256(const std::vector<double>&E,const std::vector<double>&B);
native::ArenaJson motion_record(const MotionReceipt&);
} // namespace hall_driver

namespace hall_corrective_v3 {
namespace native=raikou::mfem_backend;
// Immutable identities of one (mesh, native options) pair, computed once per
// owner so that no per-step call builds a JSON value. mesh_sha256 equals
// hall_driver::motion_mesh_sha256(mesh, options) bit for bit.
struct MotionContext {
    const axial::Mesh* mesh=nullptr;
    std::size_t maximum_carriers=100000;   // configured Hall carrier cap (Model CarrierLimits)
    double dt_max=1e-14;                   // v3.1: configured field-step upper bound (Model::dt_max)
    // v3.1: the binary v3 record digest is computed only for the accepted iterate
    // (MotionAssembly::seal_v3); the receipt of a discarded fixed-point iterate carries none.
    bool defer_record_v3=false;
    // v3 kinetic sink support (node-index ranges per KineticSink kind, filled by make_motion_context).
    int sink_lo[5][3]{},sink_hi[5][3]{},sink_supported[5]{};bool sink_support_aligned=false;std::string sink_support_text;
    native::GapPicOptions options{};
    std::string mesh_sha256;
    hall_driver::KineticPolicy policy{};
    std::string policy_sha256;
};
// Census identity checks of the sealed provider (carrier cap, ascending unique
// native IDs, nonempty packet and state identities), without hashing.
void require_census_identity(std::span<const native::ArenaRateEntry>,std::size_t maximum_carriers=100000);
// canonical-motion-census-v1 digest (unchanged encoding) with a configured cap.
std::string canonical_census_sha256(std::span<const native::ArenaRateEntry>,std::size_t maximum_carriers);
// mesh_json_dump must be axial::Mesh::json().dump() of the same mesh.
MotionContext make_motion_context(const axial::Mesh&,const native::GapPicOptions&,const std::string& mesh_json_dump,
                                  std::size_t maximum_carriers=100000);
std::string motion_mesh_sha256(const axial::Mesh&,const native::GapPicOptions&,const std::string& mesh_json_dump);

// Host admission shared by the CPU and the device provider (binding, layer,
// field, mesh/policy identity, census order and identity). Same checks, order
// and refusal texts as the sealed prepare_canonical_motion before its loops.
struct AdmittedCensus {
    std::vector<const native::ArenaRateEntry*> rows;   // ascending native ID
    std::string census_sha256;
};
// Device-computed numerical content of a motion record, carriers in ascending
// native ID order. Paths carry the canonical field-frame start and displacement.
struct MotionArrays {
    std::vector<native::GapParticle> endpoints;   // survivors only (sinks excluded), ascending ID
    std::vector<SinkRecord> sinks;                // ascending ID
    std::vector<NewbornRecord> births;            // ascending ID
    std::vector<hall_driver::CarrierPath> paths;
    std::vector<double> current_A_m,initial_charge_C,endpoint_charge_C;
    hall_driver::MotionBudget budget{};
    bool birth_scope=false;
};
// Process-wide accumulators of the host-side motion record work (single-threaded
// coordinator): census identity digest and the sealed JSON record seal.
struct MotionHostTimers {double census_hash_s=0,record_seal_s=0,record_v3_s=0;long census_hashes=0,record_seals=0,record_v3=0;};
MotionHostTimers& motion_host_timers()noexcept;
class MotionAssembly {
public:
    static AdmittedCensus admit(const MotionContext&,std::span<const native::ArenaRateEntry>,
        const hall_driver::MotionBinding&,std::span<const std::uint64_t> active_layer_ids,
        const std::vector<double>&E,const std::vector<double>&B,const VerifiedField*verified=nullptr);
    // Flight of every census row (ID order): dt for residents, end - tau for newborns,
    // and the foil face of a newborn face start (-1 otherwise). Refuses an unknown or
    // duplicate newborn ID and a birth time outside [begin, end].
    struct Flights {std::vector<double> flight,tau;std::vector<int> face_body,face_axis;std::vector<char> newborn;std::size_t count=0;};
    static Flights flights(const MotionContext&,const AdmittedCensus&,const hall_driver::MotionBinding&,std::span<const Newborn>);
    // The sealed CPU arithmetic (statement-for-statement copy of the 8dddec0 loops);
    // v3: first-hit sinks and newborn flights (ascending native ID).
    static hall_driver::MotionReceipt cpu(const MotionContext&,std::span<const native::ArenaRateEntry>,
        const hall_driver::MotionBinding&,std::span<const std::uint64_t> active_layer_ids,
        const std::vector<double>&E,const std::vector<double>&B,std::span<const Newborn> newborns={},bool birth_scope=false,const VerifiedField*verified=nullptr);
    // Receipt from device arrays after the provider has resolved every refusal in
    // the sealed order. Rechecks shapes, identities and the coverage invariants
    // that the record states, then seals the record with the sealed encoding.
    static hall_driver::MotionReceipt assemble(const MotionContext&,const AdmittedCensus&,
        const hall_driver::MotionBinding&,MotionArrays&&);
    // v2 JSON-text seal (schema canonical-full-time-total-field-motion-v2, seal over
    // "motion-receipt-record-v1" + record JSON), required by native003/native v5
    // admission. Computed once for the accepted motion; idempotent.
    static void seal_v2(hall_driver::MotionReceipt&);
    // v3.1: binary v3 record digest of a receipt built with MotionContext::defer_record_v3; idempotent.
    static void seal_v3(hall_driver::MotionReceipt&);
    // Recomputes the binary motion v3 record digest (README: v2 -> v3 field mapping).
    static std::string record_v3(const hall_driver::MotionReceipt&);
    // v3 JSON record extension (schema/coverage, sinks, births, ledgers) and the
    // digest of the sink/birth section appended to the binary v3 record.
    static native::ArenaJson extend_json_record(const hall_driver::MotionReceipt&,native::ArenaJson v2_record);
    static std::string sink_birth_digest(const hall_driver::MotionReceipt&);
};
} // namespace hall_corrective_v3
