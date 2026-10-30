#pragma once
#include "raikou/mfem_gap_pic.h"
#include <string>
namespace hall_driver {
namespace kinetic_native=raikou::mfem_backend;
inline constexpr char kinetic_boundary_api[]="native-finite-target-kinetic-policy-v1";
enum class KineticSink {none=-1,lower_foil=0,upper_foil=1,outer_xy_escape=2,short_bridge=3,axial_escape=4};
struct KineticPolicy {
    double width_m=0,gap_m=0,short_z_m=0,feed_z_m=0;
    kinetic_native::GapFoil lower{},upper{};
};
struct KineticHit {
    KineticSink sink=KineticSink::none;double fraction=1;
    int axis=-1;double plane_m=0;
    kinetic_native::PicVector outward_from_vacuum_normal{};
};
struct NeutralStreamProposal {
    kinetic_native::GapParticle endpoint;
    KineticHit hit;double begin_s=0,end_s=0,dt_s=0,flight_s=0,remaining_s=0,incident_kinetic_J=0;
    std::string kinetic_policy_sha256;
};
KineticPolicy declared_kinetic_policy(const kinetic_native::GapPicOptions&);
std::string kinetic_policy_sha256(const KineticPolicy&);
void require_kinetic_interior(const KineticPolicy&,const kinetic_native::PicVector&);
KineticHit first_kinetic_hit(const KineticPolicy&,const kinetic_native::PicVector& native_start,
    const kinetic_native::PicVector& native_displacement);
NeutralStreamProposal prepare_neutral_stream(const KineticPolicy&,
    const kinetic_native::GapParticle& full_time_neutral,double begin_s,double end_s,double dt_s);
} // namespace hall_driver
