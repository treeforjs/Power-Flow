#pragma once
// Binding of heat/source owner v6 to the Hall corrective v3 thermal owner (H5:
// one thermal owner). Include after Hall v3 root_api.hpp, the sealed Hall
// axial_surface_path_port.hpp and source owner v6 source_owner.hpp.
//
//  * One source v6 owner per Hall foil slice (body, axial cell k). Its thermal
//    geometry, material law, initial temperature and start time are those of the
//    Hall slice Coupler, and its patches are declared heat-only on that slice
//    (ThermalDeclaration, bound into the Hall material and configuration identities),
//    so the source owner reads the accepted Hall temperature through a const view.
//  * The source Trial's loads (certified emission/collection layer exchanges,
//    desorption energy debits, caller transfers) become Hall SurfaceLoads of that
//    slice; Hall applies them exactly once (watermarks) together with its Joule
//    heating and its collected-carrier wall heat.
//  * audit_single_thermal_owner refuses any second thermal store or second Joule
//    term (e.g. a source v5 private heat engine driven by a Joule coefficient).
#include <span>

namespace hall_source_binding {
namespace hc = hall_corrective_v3;
namespace s6 = electrode_heat_v6;
namespace eh = electrode_heat_v2;

// Thermal geometry of Hall slice (body, k): the geometry of ThermalV3::make_owners.
inline eh::Geometry slice_geometry(const axial::Mesh& g, int body, int k) {
    eh::Geometry geometry; geometry.nx = int(std::lround(.005 / g.dx)); geometry.ny = int(std::lround(.0002 / g.dy)); geometry.x0 = -.0025;
    geometry.y0 = body ? .0025 : -.0027; geometry.width = .005; geometry.thickness = .0002; geometry.height = g.dz;
    geometry.revision = "axial-foil-yee-heat-v2/" + std::to_string(g.nx) + "/" + std::to_string(g.nz_foil);
    geometry.electrode = "foil-" + std::to_string(body) + "-slice-" + std::to_string(k);
    return geometry;
}
inline eh::Material slice_material() {
    return eh::conditional_ss304_v2(ss304v2::Qualification::evaluated_solid_1500, ss304::DensityConvention::fixed_reference_volume_cp_as_cv_approximation,
                                    8000., eh::material_headers_sha256);
}
inline double slice_center_z(const axial::Mesh& g, int k) { return g.zmin + (k + .5) * g.dz; }
inline std::vector<eh::Patch> heat_only(std::vector<eh::Patch> p) {
    for (auto& x : p) { x.desorption_enabled = false; x.coverage_temkin = false; x.initial_particles = x.initial_low_particles = 0; }
    return p;
}
struct SliceSource {
    int body = 0, slice = 0;
    std::vector<eh::Patch> patches;
    std::vector<s6::PatchIdentity> identities;
    std::vector<s6::SourceSettings> settings;
    std::string source_hash;
};
// Heat-only declaration of every slice source, for the Hall Model.
inline hc::ThermalDeclaration declaration(std::span<const SliceSource> sources) {
    hc::ThermalDeclaration d;
    for (const auto& s : sources) {
        hc::need(!d.slices.contains({s.body, s.slice}), "one source owner per Hall foil slice");
        d.slices[{s.body, s.slice}] = {s.source_hash, heat_only(s.patches)};
    }
    return d;
}
// Source v6 owner bound to the accepted Hall slice (refuses a view of another configuration).
inline std::unique_ptr<s6::Owner> bind(const hc::Owner& hall, const SliceSource& s) {
    const auto& g = hall.mesh();
    return std::make_unique<s6::Owner>(slice_geometry(g, s.body, s.slice), slice_material(), s.patches, s.identities, s.settings, 300., s.source_hash,
                                       hall.model().drive().start(), hall.thermal().slice(s.body, s.slice));
}
// Load channel of a source transfer (ledger key in the Hall load accounting).
inline std::string channel_of(const eh::Transfer& t) {
    if (t.certified_source_exchange) return "source-layer-certified-exchange";
    if (t.producer.rfind("v6-neutral/", 0) == 0) return "desorption-neutral-release";
    if (t.incoming_K_J > 0 && t.outgoing_K_J == 0) return "particle-impact";
    return "source-transfer";
}
// Source Trial loads -> Hall surface loads of the slice (the transfers unchanged).
inline std::vector<hc::SurfaceLoad> loads(const s6::Trial& t, const SliceSource& s) {
    std::vector<hc::SurfaceLoad> out;
    for (const auto& x : t.loads()) { hc::SurfaceLoad l; l.body = s.body; l.slice = s.slice; l.channel = channel_of(x); l.transfer = x; out.push_back(std::move(l)); }
    return out;
}
// One thermal owner: Hall's trial must close Joule + wall heat + external loads, and no
// other owner may report Joule heating or a thermal change of its own.
struct OwnerReport { std::string name; double joule_J = 0, thermal_delta_J = 0; };
inline void audit_single_thermal_owner(const hc::ThermalTrialV3& hall, std::span<const OwnerReport> others) {
    const double expect = hall.admitted_joule_J + hall.wall_heat_J + hall.external_J;
    const double scale = std::abs(hall.admitted_joule_J) + std::abs(hall.wall_heat_J) + hall.external_heating_J + hall.external_cooling_J + std::abs(hall.thermal_delta_J);
    if (!(std::abs(hall.thermal_delta_J - expect) <= hall.roundoff_bound_J + 2e-12 * scale))
        throw hc::ThermalRefusal("load-accounting", "Hall thermal delta minus Joule and loads J", "hall", hall.thermal_delta_J - expect, hall.roundoff_bound_J,
                                 "the one thermal owner does not close its own Joule and loads");
    for (const auto& o : others)
        if (o.joule_J != 0 || o.thermal_delta_J != 0)
            throw hc::ThermalRefusal("double-joule", "second thermal store", o.name, o.joule_J != 0 ? o.joule_J : o.thermal_delta_J, 0,
                                     "Joule heat or electrode temperature evolved by more than one thermal owner (Hall owns sigma(T), Joule and T)");
}
}  // namespace hall_source_binding
