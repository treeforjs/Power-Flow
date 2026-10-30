#pragma once
#include "raikou/spherical_harmonics.h"
#include "raikou/backend.h"
namespace raikou::kinetics {
// CPU assembly of one cell's physical continuous RHS with nonlinear recoil/
// Landau coefficients frozen at state. All real modes are rotated to/from the
// physical field direction. The resulting CSR can be placed into a row-owned
// global kinetic block and solved with resident CPU/MPI/CUDA Krylov algebra.
[[nodiscard]] backend::CsrMatrix freeze_spherical_eedf_operator(const std::vector<CollisionProcess>& processes,
    const EedfConditions& conditions,const SphericalEedfState& state,const Vector3& electric_field_V_m);
}
