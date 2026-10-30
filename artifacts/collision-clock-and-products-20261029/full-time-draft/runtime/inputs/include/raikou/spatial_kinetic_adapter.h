#pragma once
#include "raikou/spatial_kinetics.h"
#include "raikou/spatial_partition.h"
#include <memory>

namespace raikou::spatial { struct Problem; }
namespace raikou::kinetics {
enum class KineticBoundaryOwner { ExteriorWall, PhaseContact };
struct KineticBoundaryRecord {
    std::size_t kinetic_face{0};
    KineticBoundaryOwner owner{KineticBoundaryOwner::ExteriorWall};
    std::size_t geometric_wall{std::numeric_limits<std::size_t>::max()};
    std::size_t geometric_face{std::numeric_limits<std::size_t>::max()};
    int side{-1};
    std::size_t gas_cell{0};
    std::size_t receiving_cell{std::numeric_limits<std::size_t>::max()};
    std::optional<std::size_t> absorbed_species;
    std::string energy_to_phase;
    std::size_t energy_cell{std::numeric_limits<std::size_t>::max()};
};
struct OwnedKineticBoundaryMoment {
    KineticBoundaryRecord boundary;
    // Signed outward from the gas cell. Negative values represent net injection.
    double particle_rate_s{0};
    double energy_rate_W{0};
    bool owns_gas_cell{false};
    bool owns_receiving_cell{false};
    bool owns_energy_cell{false};
};
struct OwnedKineticTrial {
    SpatialKineticResult kinetic;
    // Only records touching this rank's gas or receiving phase. Contact fluxes
    // cross a scalar face-moment halo when the receiving cell is remote.
    std::vector<OwnedKineticBoundaryMoment> boundary_moments;
};
struct KineticIncomingOverride {
    std::size_t kinetic_face{0};
    SphericalEedfState incoming;
};
struct KineticEmissionTemplate {
    std::size_t kinetic_face{0};
    std::size_t geometric_wall{0};
    SphericalEedfState per_unit_number_flux;
    double mean_incoming_energy_eV{0};
};
// Static geometry/phase setup and moment ownership for the fully coupled owned
// spatial solver. No global distribution, field, residual or solution is held.
class OwnedKineticAdapter {
public:
    OwnedKineticAdapter(const spatial::Problem& model,const spatial::mesh::Geometry& geometry,
        const std::map<std::string,std::vector<double>>& phase_fractions,
        spatial::partition::CellInterval owned_cells,const backend::ExecutionOptions& execution);
    [[nodiscard]] const KineticMesh& mesh() const;
    [[nodiscard]] const std::vector<KineticBoundaryRecord>& boundary_records() const;
    // Templates are allocated only for this rank's owned gas walls.
    [[nodiscard]] const std::vector<KineticEmissionTemplate>& emission_templates() const;
    [[nodiscard]] OwnedKineticTrial evaluate(std::span<const chemistry::Mechanism> owned_mechanisms,
        std::span<const chemistry::State> owned_candidates,
        std::span<const chemistry::Environment> owned_environments,std::span<const Vector3> owned_fields_V_m,
        std::span<const SphericalEedfState> previous,const SpatialKineticOptions& time_options,
        std::span<const SphericalEedfState> older={},std::span<const KineticIncomingOverride> incoming={}) const;
private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};
}
