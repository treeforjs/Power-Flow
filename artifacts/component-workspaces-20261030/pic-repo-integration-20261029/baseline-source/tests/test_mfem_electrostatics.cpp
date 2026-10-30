#include "raikou/mfem_electrostatics.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace {
using raikou::mfem_backend::Electrostatics;
using raikou::mfem_backend::ElectrostaticsCoordinates;
using raikou::mfem_backend::ElectrostaticsOptions;

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
void close(double value, double expected, double tolerance, const std::string& message) {
    check(std::isfinite(value) && std::abs(value - expected) <= tolerance, message);
}
mfem::Array<int> all_boundaries(const mfem::Mesh& mesh) {
    mfem::Array<int> marker(mesh.bdr_attributes.Max());
    marker = 1;
    return marker;
}
mfem::Array<int> nonaxis_boundaries(mfem::Mesh& mesh) {
    auto marker = all_boundaries(mesh);
    mfem::Array<int> vertices;
    for (int face = 0; face < mesh.GetNBE(); ++face) {
        mesh.GetBdrElementVertices(face, vertices);
        bool axis = true;
        for (int i = 0; i < vertices.Size(); ++i)
            axis = axis && mesh.GetVertex(vertices[i])[0] == 0;
        if (axis) marker[mesh.GetBdrAttribute(face) - 1] = 0;
    }
    return marker;
}
double rz_l2_error(mfem::Mesh& mesh, const mfem::GridFunction& field,
                   mfem::Coefficient& exact, int order) {
    double square = 0;
    mfem::Vector position(2);
    for (int element = 0; element < mesh.GetNE(); ++element) {
        auto& transformation = *mesh.GetElementTransformation(element);
        const auto& rule = mfem::IntRules.Get(mesh.GetElementBaseGeometry(element), order);
        for (int q = 0; q < rule.GetNPoints(); ++q) {
            const auto& point = rule.IntPoint(q);
            transformation.SetIntPoint(&point);
            transformation.Transform(point, position);
            const double difference = field.GetValue(element, point) - exact.Eval(transformation, point);
            square += difference * difference * point.weight * transformation.Weight() *
                      2 * std::numbers::pi * position[0];
        }
    }
    return std::sqrt(square);
}
template<class Function> void invalid(Function function, const std::string& message) {
    bool rejected = false;
    try { function(); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, message);
}

void refinement_and_diagnostics() {
    double previous = 0;
    for (const int n : {4, 8, 16}) {
        auto mesh = mfem::Mesh::MakeCartesian2D(n, n, mfem::Element::QUADRILATERAL);
        ElectrostaticsOptions options;
        options.planar_depth_m = 0.25;
        Electrostatics electrostatics(mesh, options);
        mfem::ConstantCoefficient epsilon(2.5), zero(0);
        mfem::FunctionCoefficient exact([](const mfem::Vector& x) {
            return x[0] * (1 - x[0]) * x[1] * (1 - x[1]);
        });
        mfem::FunctionCoefficient charge([](const mfem::Vector& x) {
            return 5 * (x[0] * (1 - x[0]) + x[1] * (1 - x[1]));
        });
        const auto result = electrostatics.solve(epsilon, charge, all_boundaries(mesh), zero);
        check(result.converged && result.residual_norm < 1e-10, "manufactured Poisson converges");
        const double error = electrostatics.potential().ComputeL2Error(exact);
        if (previous > 0) check(previous / error > 3.8, "H1 Poisson has second order L2 refinement");
        previous = error;
        close(result.total_charge_C, 2.5 / 6, 1e-12, "total charge includes extrusion depth");
        if (n == 16)
            close(result.field_energy_J, 2.5 / 360, 0.03 * 2.5 / 360,
                  "manufactured field energy approaches analytic integral");
    }
}

void moved_geometry() {
    auto mesh = mfem::Mesh::MakeCartesian2D(4, 3, mfem::Element::QUADRILATERAL);
    ElectrostaticsOptions options;
    options.planar_depth_m = 0.3;
    Electrostatics electrostatics(mesh, options);
    mfem::ConstantCoefficient epsilon(2), zero(0);
    mfem::FunctionCoefficient exact([](const mfem::Vector& x) { return x[0]; });
    const auto marker = all_boundaries(mesh);
    const auto first = electrostatics.solve(epsilon, zero, marker, exact);
    close(electrostatics.potential().ComputeL2Error(exact), 0, 1e-11, "linear field before mesh motion");
    close(first.field_energy_J, 0.3, 1e-11, "initial field energy");
    mesh.Transform([](const mfem::Vector& x, mfem::Vector& y) {
        y = x;
        y[0] *= 2;
        y[1] *= 1.5;
    });
    const auto second = electrostatics.solve(epsilon, zero, marker, exact);
    close(electrostatics.potential().ComputeL2Error(exact), 0, 1e-11, "linear field after mesh motion");
    close(second.field_energy_J, 0.9, 1e-11, "reassembly observes new physical geometry");
    mesh.UniformRefinement();
    invalid([&] { static_cast<void>(electrostatics.solve(epsilon, zero, marker, exact)); },
            "topology changes require a new operator");
}

void dielectric_flux() {
    auto mesh = mfem::Mesh::MakeCartesian2D(8, 4, mfem::Element::QUADRILATERAL);
    Electrostatics electrostatics(mesh);
    mfem::FunctionCoefficient epsilon([](const mfem::Vector& x) { return x[0] < 0.5 ? 2.0 : 8.0; });
    mfem::FunctionCoefficient exact([](const mfem::Vector& x) {
        return 3.2 * (x[0] <= 0.5 ? x[0] / 2 : 0.25 + (x[0] - 0.5) / 8);
    });
    mfem::ConstantCoefficient zero(0);
    const auto result = electrostatics.solve(epsilon, zero, all_boundaries(mesh), exact);
    close(electrostatics.potential().ComputeL2Error(exact), 0, 1e-11, "layered dielectric potential");
    close(result.field_energy_J, 1.6, 1e-10, "layered dielectric field energy");
    mfem::Vector gradient(2);
    for (int element = 0; element < mesh.GetNE(); ++element) {
        auto& transformation = *mesh.GetElementTransformation(element);
        const auto& center = mfem::Geometries.GetCenter(mesh.GetElementBaseGeometry(element));
        transformation.SetIntPoint(&center);
        electrostatics.potential().GetGradient(transformation, gradient);
        close(epsilon.Eval(transformation, center) * gradient[0], 3.2, 1e-10,
              "normal dielectric displacement is continuous");
        close(gradient[1], 0, 1e-10, "layered dielectric transverse field");
    }
}

void three_dimensions() {
    auto mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::HEXAHEDRON);
    Electrostatics electrostatics(mesh);
    mfem::ConstantCoefficient epsilon(3), zero(0);
    mfem::FunctionCoefficient exact([](const mfem::Vector& x) { return x[0]; });
    const auto result = electrostatics.solve(epsilon, zero, all_boundaries(mesh), exact);
    close(electrostatics.potential().ComputeL2Error(exact), 0, 1e-11, "three dimensional linear field");
    close(result.field_energy_J, 1.5, 1e-11, "three dimensional energy units");
}

void axisymmetric_field_and_motion() {
    auto mesh = mfem::Mesh::MakeCartesian2D(5, 4, mfem::Element::TRIANGLE);
    ElectrostaticsOptions options;
    options.coordinates = ElectrostaticsCoordinates::AxisymmetricRZ;
    Electrostatics electrostatics(mesh, options);
    mfem::ConstantCoefficient epsilon(2.5), zero(0);
    mfem::FunctionCoefficient exact([](const mfem::Vector& x) { return x[1]; });
    const auto marker = nonaxis_boundaries(mesh);
    const auto first = electrostatics.solve(epsilon, zero, marker, exact);
    close(rz_l2_error(mesh, electrostatics.potential(), exact, 6), 0, 1e-11,
          "RZ axial field with natural axis boundary");
    close(first.total_charge_C, 0, 0, "RZ zero volume charge");
    close(first.field_energy_J, 1.25 * std::numbers::pi, 1e-10,
          "RZ axial field includes full cylinder volume");
    mesh.Transform([](const mfem::Vector& x, mfem::Vector& y) {
        y = x;
        y[0] *= 2;
        y[1] *= 1.5;
    });
    const auto moved = electrostatics.solve(epsilon, zero, marker, exact);
    close(rz_l2_error(mesh, electrostatics.potential(), exact, 6), 0, 1e-10,
          "RZ axial field after mesh motion");
    close(moved.field_energy_J, 7.5 * std::numbers::pi, 2e-10,
          "RZ moved field energy uses updated radius and volume");
}

void axisymmetric_radial_refinement() {
    for (const auto geometry : {mfem::Element::TRIANGLE, mfem::Element::QUADRILATERAL}) {
        double previous = 0;
        for (const int n : {4, 8, 16}) {
            auto mesh = mfem::Mesh::MakeCartesian2D(n, n, geometry);
            ElectrostaticsOptions options;
            options.coordinates = ElectrostaticsCoordinates::AxisymmetricRZ;
            Electrostatics electrostatics(mesh, options);
            mfem::ConstantCoefficient epsilon(2.5), charge(-10);
            mfem::FunctionCoefficient exact([](const mfem::Vector& x) { return x[0] * x[0]; });
            const auto result = electrostatics.solve(epsilon, charge, nonaxis_boundaries(mesh), exact);
            const double error = rz_l2_error(mesh, electrostatics.potential(), exact, 8);
            check(error > 0, "P1 RZ radial solution has a resolved error");
            if (previous > 0) check(previous / error > 3.6, "RZ radial potential has second order L2 refinement");
            previous = error;
            close(result.total_charge_C, -10 * std::numbers::pi, 1e-10,
                  "RZ radial charge integral uses full ring measure");
            if (n == 16)
                close(result.field_energy_J, 2.5 * std::numbers::pi,
                      0.015 * 2.5 * std::numbers::pi, "RZ radial field energy approaches analytic integral");
        }
        // Quadratic elements represent r^2 exactly. This also checks the
        // radius factor is included in higher-order triangle quadrature.
        auto mesh = mfem::Mesh::MakeCartesian2D(3, 3, geometry);
        ElectrostaticsOptions options;
        options.order = 2;
        options.coordinates = ElectrostaticsCoordinates::AxisymmetricRZ;
        Electrostatics electrostatics(mesh, options);
        mfem::ConstantCoefficient epsilon(2.5), charge(-10);
        mfem::FunctionCoefficient exact([](const mfem::Vector& x) { return x[0] * x[0]; });
        const auto result = electrostatics.solve(epsilon, charge, nonaxis_boundaries(mesh), exact);
        close(rz_l2_error(mesh, electrostatics.potential(), exact, 8), 0, 1e-10,
              "P2 RZ radial potential is exact");
        close(result.field_energy_J, 2.5 * std::numbers::pi, 1e-10,
              "P2 RZ radial field energy is exact");
    }
}

void axisymmetric_validation_and_rollback() {
    ElectrostaticsOptions options;
    options.coordinates = ElectrostaticsCoordinates::AxisymmetricRZ;
    auto mesh3 = mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::HEXAHEDRON);
    invalid([&] { Electrostatics solver(mesh3, options); }, "RZ rejects a 3D mesh");
    auto mesh = mfem::Mesh::MakeCartesian2D(2, 2, mfem::Element::TRIANGLE);
    options.planar_depth_m = 2;
    invalid([&] { Electrostatics solver(mesh, options); }, "RZ rejects extrusion depth");
    options.planar_depth_m = 1;
    Electrostatics electrostatics(mesh, options);
    mfem::ConstantCoefficient epsilon(2.5), zero(0), negative(-1);
    mfem::FunctionCoefficient exact([](const mfem::Vector& x) { return x[1]; });
    const auto marker = nonaxis_boundaries(mesh);
    static_cast<void>(electrostatics.solve(epsilon, zero, marker, exact));
    mfem::Vector accepted(electrostatics.potential());
    const auto unchanged = [&] {
        mfem::Vector difference(electrostatics.potential());
        difference -= accepted;
        close(difference.Norml2(), 0, 0, "RZ failed solve preserves accepted potential");
    };
    invalid([&] { static_cast<void>(electrostatics.solve(negative, zero, marker, exact)); },
            "RZ rejects negative permittivity");
    unchanged();
    // Only the vertex lies at negative radius; interior quadrature remains
    // positive for this small translation. A quadrature-only guard misses it.
    mesh.Transform([](const mfem::Vector& x, mfem::Vector& y) {
        y = x;
        y[0] -= 1e-8;
    });
    invalid([&] { static_cast<void>(electrostatics.solve(epsilon, zero, marker, exact)); },
            "RZ rejects negative radius at mesh vertices");
    unchanged();
    mesh.Transform([](const mfem::Vector& x, mfem::Vector& y) {
        y = x;
        y[0] += 1e-8;
        y[1] = -y[1];
    });
    invalid([&] { static_cast<void>(electrostatics.solve(epsilon, zero, marker, exact)); },
            "RZ rejects inverted geometry");
    unchanged();
}

void validation_and_rollback() {
    auto mesh = mfem::Mesh::MakeCartesian2D(6, 6, mfem::Element::QUADRILATERAL);
    ElectrostaticsOptions options;
    options.maximum_iterations = 1;
    options.relative_tolerance = 1e-15;
    Electrostatics electrostatics(mesh, options);
    mfem::ConstantCoefficient epsilon(1), zero(0), one(1), negative(-1);
    auto marker = all_boundaries(mesh);
    const auto initial = electrostatics.solve(epsilon, zero, marker, zero);
    check(initial.converged, "zero field initial solve");
    mfem::FunctionCoefficient charge([](const mfem::Vector& x) { return std::exp(x[0] + x[1]); });
    bool failed = false;
    try { static_cast<void>(electrostatics.solve(epsilon, charge, marker, one)); }
    catch (const std::runtime_error&) { failed = true; }
    check(failed, "iteration exhaustion reports failure");
    close(electrostatics.potential().Norml2(), 0, 0, "failed solve preserves accepted potential");
    invalid([&] { static_cast<void>(electrostatics.solve(negative, zero, marker, zero)); },
            "negative permittivity is rejected");
    mfem::ConstantCoefficient nonfinite(std::numeric_limits<double>::quiet_NaN());
    invalid([&] { static_cast<void>(electrostatics.solve(epsilon, nonfinite, marker, zero)); },
            "nonfinite charge is rejected");
    invalid([&] { static_cast<void>(electrostatics.solve(epsilon, zero, marker, nonfinite)); },
            "nonfinite boundary potential is rejected");
    marker = 0;
    invalid([&] { static_cast<void>(electrostatics.solve(epsilon, zero, marker, zero)); },
            "unanchored Poisson system is rejected");
    options.order = 0;
    invalid([&] { Electrostatics invalid_solver(mesh, options); }, "invalid polynomial order is rejected");
}
} // namespace

int main() {
    try {
        refinement_and_diagnostics();
        moved_geometry();
        dielectric_flux();
        three_dimensions();
        axisymmetric_field_and_motion();
        axisymmetric_radial_refinement();
        axisymmetric_validation_and_rollback();
        validation_and_rollback();
        std::cout << "MFEM electrostatics tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
