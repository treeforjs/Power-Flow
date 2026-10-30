#pragma once
#include "raikou/spatial.h"

namespace raikou::spatial::boundary {
struct SheathFlux {
    double particles_m2_s = 0;
    double density_derivative_m_s = 0;
    double cell_potential_derivative_m2_s_V = 0;
    double wall_potential_derivative_m2_s_V = 0;
};
// The temperature and transmission are explicitly authored closure parameters.
// Maxwellian number flux includes exp(-max(q*DeltaPhi,0)/T). Bohm accepts only
// positive ions at an ion-attracting wall; invalid assumptions are rejected.
[[nodiscard]] SheathFlux evaluate_sheath(const Boundary::Sheath& closure,
    const chemistry::Species& species, double density_m3,
    double cell_potential_V, double wall_potential_V);
// The `drift` closure's thermal part: transmission times the one-way speed v/4 of a
// Maxwellian at temperature_eV, or at the local characteristic energy D/mu when
// temperature_eV is 0. The solver adds the upwind drift toward the wall.
[[nodiscard]] double drift_thermal_speed(const Boundary::Sheath& closure,
    const chemistry::Species& species, const Transport& transport);

struct IncomingTemplate {
    kinetics::SphericalEedfState per_unit_number_flux;
    double mean_incoming_energy_eV = 0;
};
// Uses the same discrete P_N half-range flux as the transport operator. Dividing
// by an analytic continuum thermal speed would mis-normalize a truncated grid.
[[nodiscard]] IncomingTemplate normalize_incoming(
    const kinetics::SphericalEedfState& authored_distribution,
    const kinetics::Vector3& outward_normal);
} // namespace raikou::spatial::boundary
