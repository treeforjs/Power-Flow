#pragma once
#include <span>
#include <vector>

namespace raikou::time {
struct BdfWeights {
    double current = 1, previous = -1, older = 0;
    int order = 1;
};
// Residual: current*y_new + previous*y_old + older*y_older - dt*f_new.
// BDF1 permits previous_dt=0. BDF2 requires a positive previous timestep.
// Coefficients are valid for positive ratios; the caller must enforce a suitable
// variable-step stability bound (the classical BDF2 ratio bound is 1+sqrt(2)).
[[nodiscard]] BdfWeights coefficients(double dt, double previous_dt, int order);

struct ErrorEstimate {
    double normalized_error = 0; // max component error / (atol+rtol*scale)
    double full_error_multiplier = 0;
    std::vector<double> component_errors; // signed estimate for full-step error
};
// Controls the full-step endpoint, using full-refined as the embedded defect.
// atol is one positive scalar or one positive value per component; rtol>=0.
// Empty matching states/tolerances return zero error for an empty MPI partition.
// The timestep ratio is ignored for BDF1 startup.
// Both half solves must use the same order as the full solve. For BDF2 the first
// half solve reuses the accepted old/older history (ratio dt/(2*previous_dt));
// the second uses the new half state and old state with ratio 1. No accepted
// history may be overwritten until the outer candidate is accepted.
// The multiplier is 2 for BDF1, and for BDF2, r=dt/previous_dt,
//   12(1+r)^3 / (6r^3+21r^2+22r+8),
// which equals 32/19 at r=1. It is exact for constant leading polynomial forcing;
// stiff/nonlinear dynamics still require recomputing each trial's error.
[[nodiscard]] ErrorEstimate estimate_error(std::span<const double> full_step,
    std::span<const double> two_half_steps, std::span<const double> absolute_tolerances,
    double relative_tolerance, int order, double dt_over_previous_dt = 1);

// Uses the local order+1 exponent with bounded changes. A rejected trial always
// shrinks, by at least a factor 1/2 unless min_shrink itself is greater. Positive
// infinity denotes a failed nonlinear solve and returns min_shrink when rejected;
// NaN and invalid accepted errors are rejected as input errors.
[[nodiscard]] double recommended_factor(double error, int order, bool accepted,
    double min_shrink = .1, double max_grow = 2, double safety = .9);
} // namespace raikou::time
