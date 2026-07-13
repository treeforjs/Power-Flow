#pragma once

#include <vector>

namespace mhd_core {

struct HalfRangeMoments {
    double number_flux;
    double mean_normal_speed;
    double mean_speed;
};

HalfRangeMoments half_range_maxwell_moments(double temperature_k, double mass_kg);

std::vector<double> laplacian2d(
    const std::vector<double>& values,
    int ny,
    int nx,
    double dx,
    double dy);

}  // namespace mhd_core
