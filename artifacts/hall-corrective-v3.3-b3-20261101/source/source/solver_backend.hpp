#pragma once
// Explicit field-solve backend interface for the Hall implicit edge solve
//   A dE' = b,  A = I + D^-1/2 (dt/2 (S+Y) + dt^2/4 C^T H C) D^-1/2  (scaled, v = sqrt(M) dE)
// (8dddec0 field_model.hpp:79-88). A backend is selected by configuration; its
// descriptor digest is part of the Model configuration identity, so a CPU and
// a GPU run have different configuration and field-state identities.
//
// Solver contract, identical for every backend: the returned increment x must
// satisfy the INDEPENDENT FP64 residual check performed by Model::prepare with
// the sealed matrix-free apply and the sealed long-double norm,
//   ||b - A x||_2 <= max(1e-30, 2e-13 ||b||_2),
// before any other use. The model's own Gauss-increment and energy guards do
// not detect a 0.1 % increment error on near-equilibrium steps, so this check
// is never skipped, including for the sealed CPU solver.
#include "field_contract.hpp"
#include <chrono>
#include <functional>
#include <memory>

namespace hall_corrective_v3 {
using Apply=std::function<Vec(const Vec&)>;
// Operator operands handed to a backend (all computed by Model::prepare with the
// sealed expressions): mesh, dt, per-edge conduction S, sqrt(M) and the Jacobi
// diagonal. A backend that assembles its own operator must reproduce them.
struct SolveContext {
    const Mesh& mesh;double dt;const Vec& conduction;const Vec& sqrtM;const Vec& diagonal;
};
struct SolveRecord {
    std::string backend;
    int iterations=-1,corrections=0;
    double rhs_norm=0,target=0,backend_residual=-1,independent_residual=0;
    double solve_s=0,check_s=0;
};
class FieldSolveBackend {
public:
    virtual ~FieldSolveBackend()=default;
    virtual const std::string& name()const noexcept=0;
    // Digest of the canonical descriptor (algorithm, parameters, source identity).
    virtual const std::string& descriptor_sha256()const noexcept=0;
    virtual Vec solve(const Apply&apply,const Vec&rhs,int maxit,const SolveContext&,SolveRecord&)=0;
};
// The sealed foil::solve Jacobi PCG (8dddec0 foil_em.hpp:56-63), unchanged.
class CpuPcgBackend final:public FieldSolveBackend {
    std::string name_="cpu-pcg",descriptor_;
public:
    CpuPcgBackend(){Digest d("hall-v2-field-solve-backend-v1");d.text(name_);
        d.text("sealed foil::solve: Jacobi-preconditioned CG, recursive stop at 0.2*target, true residual "
               "||b-Ax||<=max(1e-30,2e-13||b||), at most three correction cycles, long-double dot products");
        d.text("hall-finite-foil-em-20261029/foil_em.hpp@8dddec0610982778bf02bdfe697b70a37129c012");descriptor_=d.hex();}
    const std::string& name()const noexcept override{return name_;}
    const std::string& descriptor_sha256()const noexcept override{return descriptor_;}
    Vec solve(const Apply&apply,const Vec&rhs,int maxit,const SolveContext&c,SolveRecord&r)override{
        const std::function<Vec(const Vec&)>&f=apply;r.backend=name_;
        return foil::solve(f,rhs,c.diagonal,maxit);
    }
};
inline std::shared_ptr<FieldSolveBackend> make_cpu_pcg(){return std::make_shared<CpuPcgBackend>();}
} // namespace hall_corrective_v3
