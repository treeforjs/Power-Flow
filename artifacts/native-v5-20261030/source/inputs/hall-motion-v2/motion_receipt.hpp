#pragma once
#include "raikou/mfem_pic_arena.h"
#include "kinetic_boundary.hpp"
#include <span>

namespace axial {struct Mesh;}
namespace hall_driver {
namespace native=raikou::mfem_backend;
inline constexpr char motion_api[]="canonical-full-time-total-field-motion-v2";
enum class ParticlePhase {unbound,full_time_position_and_proper_velocity};
// Supplied by the private arena outer transaction AFTER its first exposure.
// The field component cannot turn a legacy arena phase into this phase.
struct MotionBinding {
    ParticlePhase phase=ParticlePhase::unbound;
    std::uint64_t accepted_step=0,staged_arena_revision=0,trial_stage=0;
    double begin_s=0,end_s=0,dt_s=0;
    std::string staged_arena_sha256,trial_field_sha256,mesh_sha256;
    std::string material_sha256,species_sha256,common_clock_sha256,kinetic_policy_sha256;
};
// Native (0..width,0..gap,z) -> field (-width/2..width/2,-gap/2..gap/2,z).
// This is a translation only. E, B and proper velocity retain their axes.
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

// Move-only immutable proposal. It does not own an accepted arena, allocate IDs,
// mutate Count2, solve native Poisson, invoke native advance, or publish anything.
// v2 covers every resident carrier (including q=0) exactly once and rejects any
// active layer, source population change, or boundary hit. Native owner validates
// its private baseline and consumes endpoints in its own typed motion operation.
class MotionReceipt {
    friend MotionReceipt prepare_canonical_motion(const axial::Mesh&,
        const native::GapPicOptions&,std::span<const native::ArenaRateEntry>,
        const MotionBinding&,std::span<const std::uint64_t>,
        const std::vector<double>&,const std::vector<double>&);
    MotionBinding binding_;TranslationFrame frame_;MotionBudget budget_;KineticPolicy kinetic_policy_;
    std::string census_sha256_,seal_sha256_;
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
    const std::string& seal_sha256()const noexcept{return seal_sha256_;}
    std::span<const native::GapParticle> endpoints()const noexcept{return endpoints_;}
    std::span<const CarrierPath> paths()const noexcept{return paths_;}
    const std::vector<double>& current_A_m()const noexcept{return current_A_m_;}
    const std::vector<double>& initial_charge_C()const noexcept{return initial_charge_C_;}
    const std::vector<double>& endpoint_charge_C()const noexcept{return endpoint_charge_C_;}
};
// SHA256 over a canonical binary encoding (length-prefixed UTF-8 strings,
// little-endian uint64 and IEEE binary64 bits, sorted native IDs). The native
// owner uses this exact helper to bind its before census, never rounded weights.
std::string canonical_census_sha256(std::span<const native::ArenaRateEntry>);
std::string motion_mesh_sha256(const axial::Mesh&,const native::GapPicOptions&);
std::string total_field_sha256(const std::vector<double>&E,const std::vector<double>&B);
native::ArenaJson motion_record(const MotionReceipt&);
} // namespace hall_driver
