#pragma once
// Hall motion receipt, v2 equivalent of the sealed hall-canonical-driver-20261029/
// motion_receipt.hpp (commit 8dddec0). Every type, member, accessor and free
// function below has the sealed name, layout and meaning, so native003 and
// native v5 (whose adapters only forward-declare hall_driver::MotionReceipt and
// call hall_driver::motion_record) admit it unchanged. The one difference is
// construction: besides the CPU provider prepare_canonical_motion, the v2
// assembler hall_corrective_v2::MotionAssembly is a friend, so a device-computed
// motion (GPU push) yields the identical immutable receipt and record. The
// receipt additionally carries the binary "motion v3" record digest (same
// content, canonical binary encoding, computed for every fixed-point iterate);
// the v2 JSON-text seal required by native admission is computed only for the
// accepted motion (MotionAssembly::seal_v2). Native adapters never access the
// receipt layout (they only forward-declare it and call motion_record).
//
// A program links either this header with source/motion_provider.cpp or the
// sealed header with the sealed native_arena_motion.cpp, never both.
#include "raikou/mfem_pic_arena.h"
#include "hall-canonical-driver-20261029/kinetic_boundary.hpp"
#include <span>
#include <string>
#include <vector>

namespace axial {struct Mesh;}
namespace hall_corrective_v2 {class MotionAssembly;}
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
    friend class ::hall_corrective_v2::MotionAssembly;
    MotionBinding binding_;TranslationFrame frame_;MotionBudget budget_;KineticPolicy kinetic_policy_;
    std::string census_sha256_,seal_sha256_,record_v3_sha256_;
    std::vector<native::GapParticle> endpoints_;
    std::vector<CarrierPath> paths_;
    std::vector<double> current_A_m_,initial_charge_C_,endpoint_charge_C_;
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
};
std::string canonical_census_sha256(std::span<const native::ArenaRateEntry>);
std::string motion_mesh_sha256(const axial::Mesh&,const native::GapPicOptions&);
std::string total_field_sha256(const std::vector<double>&E,const std::vector<double>&B);
native::ArenaJson motion_record(const MotionReceipt&);
} // namespace hall_driver

namespace hall_corrective_v2 {
namespace native=raikou::mfem_backend;
// Immutable identities of one (mesh, native options) pair, computed once per
// owner so that no per-step call builds a JSON value. mesh_sha256 equals
// hall_driver::motion_mesh_sha256(mesh, options) bit for bit.
struct MotionContext {
    const axial::Mesh* mesh=nullptr;
    std::size_t maximum_carriers=100000;   // configured Hall carrier cap (Model CarrierLimits)
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
    std::vector<native::GapParticle> endpoints;
    std::vector<hall_driver::CarrierPath> paths;
    std::vector<double> current_A_m,initial_charge_C,endpoint_charge_C;
    hall_driver::MotionBudget budget{};
};
// Process-wide accumulators of the host-side motion record work (single-threaded
// coordinator): census identity digest and the sealed JSON record seal.
struct MotionHostTimers {double census_hash_s=0,record_seal_s=0,record_v3_s=0;long census_hashes=0,record_seals=0,record_v3=0;};
MotionHostTimers& motion_host_timers()noexcept;
class MotionAssembly {
public:
    static AdmittedCensus admit(const MotionContext&,std::span<const native::ArenaRateEntry>,
        const hall_driver::MotionBinding&,std::span<const std::uint64_t> active_layer_ids,
        const std::vector<double>&E,const std::vector<double>&B);
    // The sealed CPU arithmetic (statement-for-statement copy of the 8dddec0 loops).
    static hall_driver::MotionReceipt cpu(const MotionContext&,std::span<const native::ArenaRateEntry>,
        const hall_driver::MotionBinding&,std::span<const std::uint64_t> active_layer_ids,
        const std::vector<double>&E,const std::vector<double>&B);
    // Receipt from device arrays after the provider has resolved every refusal in
    // the sealed order. Rechecks shapes, identities and the coverage invariants
    // that the record states, then seals the record with the sealed encoding.
    static hall_driver::MotionReceipt assemble(const MotionContext&,const AdmittedCensus&,
        const hall_driver::MotionBinding&,MotionArrays&&);
    // v2 JSON-text seal (schema canonical-full-time-total-field-motion-v2, seal over
    // "motion-receipt-record-v1" + record JSON), required by native003/native v5
    // admission. Computed once for the accepted motion; idempotent.
    static void seal_v2(hall_driver::MotionReceipt&);
    // Recomputes the binary motion v3 record digest (README: v2 -> v3 field mapping).
    static std::string record_v3(const hall_driver::MotionReceipt&);
};
} // namespace hall_corrective_v2
