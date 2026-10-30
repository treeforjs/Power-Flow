#pragma once
#include "raikou/spherical_harmonics.h"
#include "raikou/backend.h"
#include <span>

namespace raikou::spatial { struct Mesh; }

namespace raikou::kinetics {

struct RzGeometry {
    // One value or one per spatial cell. This is the cell-average 1/r under
    // the cylindrical volume measure, A_meridional / integral(r dA), rather
    // than the reciprocal of a centroid radius. It remains finite at the axis.
    std::vector<double> inverse_radius_m1;
};

// Full-revolution physical cell volumes, V=2*pi*integral(r dA).
[[nodiscard]] RzGeometry make_rz_geometry(std::span<const double> meridional_areas_m2,
    std::span<const double> cell_volumes_m3);
[[nodiscard]] RzGeometry rz_geometry(const spatial::Mesh& mesh);

// Cylindrical velocity coordinates use x=r, y=theta, z=axial. The conservative
// spatial streaming LHS contains -1/r * d_beta(v_theta*f). Consequently this
// local RHS has angular entries -<d_beta(R_i) * v_theta * f>/r. R00 rows are
// identically zero: this geometry term changes neither number nor energy.
// Combine with radial FV fluxes using physical cylindrical areas and volumes.
// Retains every real m and the native slot*terms^2 even/odd energy staggering.
[[nodiscard]] std::vector<std::vector<double>> rz_geometry_rhs(const RzGeometry& geometry,
    const std::vector<SphericalEedfState>& states, std::size_t threads = 1);

// Per-cell linear RHS operator in the same flattened layout. Useful for an
// implicit monolithic assembly; no time discretization is included. Its even/
// odd velocity overlap is the same mass-adjoint pair used by face streaming.
[[nodiscard]] backend::CsrMatrix rz_geometry_matrix(const EedfGrid& grid,
    double inverse_radius_m1);

} // namespace raikou::kinetics
