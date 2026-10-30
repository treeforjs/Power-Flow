#pragma once

#include "raikou/spatial.h"
#include "raikou/backend.h"
#include <array>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

namespace raikou::spatial::mesh {

using Point = std::array<double, 3>;

struct Face {
    std::size_t a = 0, b = 0;
    double area = 0;
    Point normal{}; // unit normal directed from cell a into cell b
    Point center{};
    double normal_distance = 0; // dot(center[b] - center[a], normal)
    Point delta{}; // center[b] - center[a]
};

struct Wall {
    std::size_t cell = 0;
    int side = 0;
    double area = 0;
    Point normal{}; // unit normal directed out of the cell
    Point center{};
    double normal_distance = 0;
    Point delta{}; // wall center - cell center
};

struct Geometry {
    int dimension = 2;
    bool axisymmetric = false;
    std::vector<double> volume;
    std::vector<Point> centers;
    std::vector<Face> faces;
    std::vector<Wall> walls;
    // Minimum alignment of a center displacement and its outward face normal.
    // Positive convex cells do not require orthogonality.
    double minimum_normal_cosine = 1;
};

// An affine scalar form without a constant term. The caller may substitute any
// linear expression for a boundary slot (including a dielectric wall voltage).
// Boundary slots have exactly the same order as Geometry::walls.
struct ScalarStencil {
    std::vector<std::pair<std::size_t, double>> cells;
    std::vector<std::pair<std::size_t, double>> boundaries;
    [[nodiscard]] double evaluate(std::span<const double> cell_values,
        std::span<const double> boundary_values) const;
};

enum class BoundaryCondition { Dirichlet, Neumann }; // homogeneous Neumann

struct FluxStencils {
    std::vector<std::array<ScalarStencil, 3>> cell_gradients;
    // Scalar normal derivatives; multiply by face area and material coefficient
    // before assembly. Assemble an interior face once with opposite cell signs.
    std::vector<ScalarStencil> interior;
    std::vector<ScalarStencil> boundary;
    std::vector<double> gradient_condition_numbers;
    bool has_nonorthogonal_corrections = false;
};

// Structured Cartesian/RZ grids and supplied convex tri3/quad4/tet4/hex8 meshes
// with planar faces. All supplied exterior faces require an explicit side label.
// Uses geometric volume centroids, including radial-volume weighting in RZ.
// Throws invalid_argument on nonconvex, nonmanifold or numerically degenerate
// geometry. Hex connectivity follows the standard two four-node rings.
[[nodiscard]] Geometry geometry(const Mesh& mesh);

// Implicit inverse-distance-weighted least-squares gradients and conservative
// nonorthogonal correction. Empty conditions means all walls are Dirichlet.
// Each interior derivative is
//   (u_b-u_a)/dn + grad_face dot (normal-delta/dn).
// Boundary Dirichlet derivatives use the same correction with the cell gradient;
// homogeneous Neumann derivatives are exactly zero. The method is affine exact,
// but general skew corrections do not guarantee a monotone/M-matrix operator.
[[nodiscard]] FluxStencils normal_derivatives(const Geometry& geometry,
    std::span<const BoundaryCondition> conditions = {});

// Phase-local reconstruction. Fractions are scalar or one per cell in [0,1];
// only their zero/nonzero status is used. Active/absent interior contacts become
// homogeneous-Neumann walls of the active subdomain. Output arrays and stencil
// coefficient IDs retain the original geometry's global ordering; absent-cell,
// absent-face and barrier-face entries are zero. Geometric area/aperture and
// phase-volume scaling remain the caller's responsibility.
[[nodiscard]] FluxStencils normal_derivatives_active(const Geometry& geometry,
    std::span<const BoundaryCondition> conditions,
    std::span<const double> active_cell_fractions);

enum class EllipticBoundaryKind { Dirichlet, Neumann, Robin };

struct EllipticBoundary {
    EllipticBoundaryKind kind = EllipticBoundaryKind::Neumann;
    double value = 0; // Dirichlet value or Robin backing value
    double capacitance = 0; // nonnegative Robin coefficient C
    double source = 0; // Robin source sigma
};

struct EllipticAssembly {
    backend::CsrMatrix matrix;
    std::vector<double> rhs;
    // Cell unknowns come first. Each Robin wall then owns one unknown; other
    // wall entries are -1. A Robin wall value cannot in general be eliminated
    // by an independent two-point expression on a skew mesh.
    std::vector<int> robin_wall_dofs;
    FluxStencils stencils;
    std::vector<double> interior_conductivity;
};

// Assemble -div(k grad u) = rhs_density. Conductivity is one positive scalar or
// one value per cell. The RHS density is empty (zero) or one value per cell.
// Dirichlet values are substituted symbolically; Neumann means homogeneous.
// Each Robin wall appends the unintegrated equation
//   k_cell * du/dn + C * (u_wall - value) - source = 0.
// Interior k is the positive normal-distance-weighted harmonic average. This
// gives the exact series resistance of orthogonal layered materials. General
// skew coefficient jumps need interface-aware tangential reconstruction; the
// present continuous-gradient least-squares correction does not ensure that
// transmission property. Pure Neumann systems retain their constant nullspace.
[[nodiscard]] EllipticAssembly assemble_elliptic(const Geometry& geometry,
    std::span<const double> conductivity, std::span<const double> rhs_density,
    std::span<const EllipticBoundary> boundaries);

} // namespace raikou::spatial::mesh
