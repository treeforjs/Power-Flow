#pragma once
#include <cmath>
#include <stdexcept>

namespace raikou::native_water::v8 {
// Launch-only requests are applied after replay has verified the old checkpoint.
// A checkpoint without a new request retains its original stopping condition.
struct RunTargetRequest {
  bool end_requested = false, window_requested = false, steps_requested = false;
  double end_s = 0, window_s = 0;
  long steps = 20;
  void validate() const {
    if (int(end_requested) + int(window_requested) + int(steps_requested) > 1)
      throw std::invalid_argument("choose only one of --end-time, --window, or --steps");
    if (end_requested && !(std::isfinite(end_s) && end_s > 0))
      throw std::invalid_argument("--end-time must be finite and positive");
    if (window_requested && !(std::isfinite(window_s) && window_s > 0))
      throw std::invalid_argument("--window must be finite and positive");
    if (steps_requested && steps < 0) throw std::invalid_argument("--steps must be nonnegative");
  }
  double after_restore(double saved_end_s, double accepted_time_s) const {
    validate();
    if (steps_requested) return 0;
    const double result = window_requested ? accepted_time_s + window_s : end_requested ? end_s : saved_end_s;
    if ((window_requested || end_requested) && !(std::isfinite(result) && result > accepted_time_s))
      throw std::invalid_argument("requested continuation end must be finite and later than the accepted checkpoint time");
    return result;
  }
};
}  // namespace raikou::native_water::v8
