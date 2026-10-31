// relocation.hpp -- CPU reference of the relocation current of a merge transaction
// and of the Hall endpoint continuity check it must pass.
//
// The arithmetic is a port (same operations, same order) of
//  * axial::Path electric weights and Path::deposit (hall-finite-foil-em-20261029/
//    axial_particles.hpp:6-37): cut at cell planes, two-point Gauss per segment,
//    J[e] += charge*D[axis]*weight/dt in weight order;
//  * axial::Mesh indexing (axial_mesh.hpp: node, ex, ey, ez; edge a->b, length);
//  * Hall corrective v3 Words (field_contract.hpp:24-37), deposit_charge and
//    charge_of_particles (census_charge.hpp:15-45), incidence (field_contract.hpp)
//    and the endpoint continuity defect and scale of Model::prepare
//    (field_model.hpp:139-154).
// Relocation records are deposited in transaction order into a separate edge
// vector J_reloc; the field step uses plasma[e] = J_motion[e] + J_reloc[e].
#pragma once
#include "merge_split.hpp"

#include <vector>

namespace merge_split {

struct HallMesh {
    int nx, ny, nz, axial, air, nn, ne, nex, ney;
    double side, depth, dx, dy, dz, zmin;
    explicit HallMesh(int n, int axial = 4, double side = .009, int air = 1, double depth = .025);
    int node(int i, int j, int k) const { return (k * (ny + 1) + j) * (nx + 1) + i; }
    int ex(int i, int j, int k) const { return (k * (ny + 1) + j) * nx + i; }
    int ey(int i, int j, int k) const { return nex + (k * ny + j) * (nx + 1) + i; }
    int ez(int i, int j, int k) const { return nex + ney + (k * (ny + 1) + j) * (nx + 1) + i; }
};

// Port of Hall corrective v3 Words (double-word charge arithmetic).
struct HWords {
    double hi = 0, lo = 0;
    HWords() = default;
    explicit HWords(double v) : hi(v) {}
    HWords(double a, double b);
    long double wide() const { return (long double)hi + lo; }
};
HWords operator+(HWords a, HWords b);
HWords operator-(HWords a);
HWords operator-(HWords a, HWords b);
HWords operator*(HWords a, double b);
HWords operator/(HWords a, double b);

// J_reloc (A m per edge, the Hall plasma-current unit) of the relocation records.
std::vector<double> relocation_current(const HallMesh&, const Relocation&, double dt);

// Free charge of a carrier set (charge_of_particles arithmetic: ascending native ID,
// Words(w)*q, CIC words per node, gross rounded upward).
struct FreeChargeW {
    std::vector<HWords> net;
    std::vector<double> gross;
};
FreeChargeW free_charge(const HallMesh&, const MsParams&, const Carriers&);

// Endpoint continuity of the relocation alone: after (post-transaction census) -
// before (census) - incidence(J_reloc) dt, per node, with the Hall scale
// gross_before + gross_after + |incidence dt|.
struct ContinuityReport {
    double max_defect_C = 0;        // max |defect|
    double max_ratio = 0;           // max |defect| / (1e-12 scale): the Hall relative term alone must be <= 1
    double max_net_change_C = 0;    // max |after - before| (charge moved by the relocation)
    std::size_t nodes_changed = 0;
};
ContinuityReport relocation_continuity(const HallMesh&, const MsParams&, const Carriers& before, const Carriers& after,
                                       const std::vector<double>& J, double dt);

// Field work on the relocation current with the Hall plasma-work expression:
// sum over edges ascending of dt*E[e]*J[e].
double relocation_work_J(const std::vector<double>& E, const std::vector<double>& J, double dt);
// Edge field of a smooth test field: E(x) = E0 + G x (componentwise linear), sampled at edge midpoints.
std::vector<double> test_edge_field(const HallMesh&, const double E0[3], const double G[3][3]);

}  // namespace merge_split
