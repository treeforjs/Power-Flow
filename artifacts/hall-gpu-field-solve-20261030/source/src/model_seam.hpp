// Solver seam for the sealed Hall field model (commit 8dddec0), without any
// source change.
//
// Model::prepare calls foil::solve(apply, rhs, diagonal, 200) with a lambda
// (field_model.hpp:88). The sealed foil::solve takes a std::function, so the
// lambda needs a user-defined conversion. A function template declared here,
// after foil_em.hpp/axial_fields.hpp/axial_thermal_v2.hpp and before
// field_model.hpp, is an exact match and is selected for that one call. Calls
// compiled earlier (the DC seed solves in axial_fields.hpp) still resolve to
// the sealed solver. Everything Model::prepare computes before and after the
// solve, including its Gauss, divB, Gauss-increment and energy guards, runs
// unchanged on whatever increment the installed backend returns.
//
// This header must be the only route by which field_model.hpp enters a
// program, otherwise Model::prepare would have two different definitions.
#pragma once
#include "development/hall-field-correction-20261030/field_contract.hpp"
#include "development/hall-canonical-driver-20261029/axial_thermal_v2.hpp"
#include <functional>
#include <type_traits>

namespace hgf_seam {
using Vec = foil::Vec;
using Apply = std::function<Vec(const Vec&)>;
using Backend = std::function<Vec(const Apply& apply, const Vec& b, const Vec& diagonal, int maxit)>;
inline Backend& backend() { static Backend b; return b; }
inline std::size_t& calls() { static std::size_t n = 0; return n; }
}  // namespace hgf_seam

namespace foil {
template <class F>
    requires(std::is_invocable_r_v<Vec, const F&, const Vec&> &&
             !std::is_same_v<std::remove_cvref_t<F>, std::function<Vec(const Vec&)>>)
Vec solve(const F& apply, const Vec& b, const Vec& diagonal, int maxit) {
    auto& backend = hgf_seam::backend();
    if (!backend) throw std::runtime_error("field solve seam: no backend installed");
    ++hgf_seam::calls();
    return backend(hgf_seam::Apply(std::cref(apply)), b, diagonal, maxit);
}
}  // namespace foil

#include "development/hall-field-correction-20261030/field_model.hpp"

namespace hgf_seam {
// The sealed CPU solver, selected explicitly through its std::function signature.
inline Vec sealed_solve(const Apply& apply, const Vec& b, const Vec& diagonal, int maxit) {
    const std::function<Vec(const Vec&)>& f = apply;
    return foil::solve(f, b, diagonal, maxit);
}
}  // namespace hgf_seam
