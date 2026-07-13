#include "reduced.hpp"

#include <cmath>
#include <stdexcept>

namespace mhd_core {

namespace {
constexpr double kB = 1.380649e-23;
constexpr double pi = 3.141592653589793238462643383279502884;
}

HalfRangeMoments half_range_maxwell_moments(double temperature_k, double mass_kg) {
    if (temperature_k <= 0.0 || mass_kg <= 0.0) {
        throw std::invalid_argument("temperature and mass must be positive");
    }
    const double thermal = std::sqrt(kB * temperature_k / mass_kg);
    return {
        thermal / std::sqrt(2.0 * pi),
        std::sqrt(pi / 2.0) * thermal,
        2.0 * std::sqrt(2.0 / pi) * thermal,
    };
}

std::vector<double> laplacian2d(
    const std::vector<double>& values,
    int ny,
    int nx,
    double dx,
    double dy) {
    if (ny <= 2 || nx <= 2) {
        throw std::invalid_argument("grid must have at least 3x3 cells");
    }
    if (static_cast<int>(values.size()) != ny * nx) {
        throw std::invalid_argument("values size does not match grid");
    }
    std::vector<double> out(values.size(), 0.0);
    const double idx2 = 1.0 / (dx * dx);
    const double idy2 = 1.0 / (dy * dy);
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            const int k = j * nx + i;
            out[k] = (values[k - 1] - 2.0 * values[k] + values[k + 1]) * idx2
                   + (values[k - nx] - 2.0 * values[k] + values[k + nx]) * idy2;
        }
    }
    return out;
}

}  // namespace mhd_core
