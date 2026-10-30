#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace raikou {

/**
 * @brief Electro-acoustic plasma spark channel expansion model (Braginskii).
 * Implements Phase 4 / SHOCK-FRACTURE Stage 1.
 */
class SparkBreakdownSolver {
public:
    struct SparkConfig {
        double peak_current{20000.0};      // 20 kA peak discharge current
        double pulse_duration_s{5e-6};     // 5 microsecond pulse
        double initial_density{1000.0};    // kg/m^3 (water/dielectric)
        double plasma_conductivity{2.5e4}; // S/m (Spitzer-like ionized channel)
        double channel_length{0.01};       // 10 mm spark gap
    };

    struct SparkState {
        double time_s{0.0};
        double current_A{0.0};
        double channel_radius_m{1e-4};
        double pressure_Pa{1e5};
        double temperature_K{300.0};
        double deposited_energy_J{0.0};
        double accumulated_action{0.0};
    };

    SparkBreakdownSolver() = default;

    /// Evaluates spark channel radius using Braginskii scaling:
    /// r(t) = alpha * (2 / (pi * rho0 * eta))^(1/6) * [integral I(tau)^2 dtau]^(1/6) * t^(1/3)
    static double compute_channel_radius(
        double integral_i2_dt,
        double rho0 = 1000.0,
        double eta = 2.5e4,
        double time_s = 5e-6) {

        if (!std::isfinite(integral_i2_dt) || integral_i2_dt < 0.0 ||
            !std::isfinite(rho0) || rho0 <= 0.0 ||
            !std::isfinite(eta) || eta <= 0.0 ||
            !std::isfinite(time_s) || time_s < 0.0) {
            throw std::invalid_argument("spark radius inputs must be finite and physical");
        }

        double factor = 2.0 / (3.141592653589793 * rho0 * std::max(0.1, eta));
        double term = factor * std::max(0.0, integral_i2_dt);
        double t_factor = (time_s > 0.0) ? std::pow(time_s, 1.0 / 3.0) : 1.0;
        return std::max(1e-4, 0.548 * std::pow(term, 1.0 / 6.0) * t_factor);
    }

    /// Advances spark discharge state by dt:
    static SparkState step_spark_discharge(
        const SparkState& current,
        double dt,
        const SparkConfig& config) {

        const auto positive = [](double value) { return std::isfinite(value) && value > 0.0; };
        const auto nonnegative = [](double value) { return std::isfinite(value) && value >= 0.0; };
        if (!positive(dt) || !nonnegative(config.peak_current) ||
            !positive(config.pulse_duration_s) || !positive(config.initial_density) ||
            !positive(config.plasma_conductivity) || !positive(config.channel_length) ||
            !nonnegative(current.time_s) || !nonnegative(current.current_A) ||
            !positive(current.channel_radius_m) || !nonnegative(current.pressure_Pa) ||
            !positive(current.temperature_K) || !nonnegative(current.deposited_energy_J) ||
            !nonnegative(current.accumulated_action)) {
            throw std::invalid_argument("spark state, configuration and time step must be finite and physical");
        }

        SparkState next = current;
        next.time_s += dt;

        // Damped sinusoidal pulse: I(t) = I_peak * sin(pi * t / T) * exp(-2t / T)
        if (next.time_s <= config.pulse_duration_s) {
            double phase = 3.141592653589793 * next.time_s / config.pulse_duration_s;
            double decay = std::exp(-2.0 * next.time_s / config.pulse_duration_s);
            next.current_A = config.peak_current * std::sin(phase) * decay;
        } else {
            next.current_A = 0.0;
        }

        // Action integral: \int I^2 dt
        next.accumulated_action += next.current_A * next.current_A * dt;
        if (next.accumulated_action > 0.0) {
            next.channel_radius_m = std::max(1e-4, compute_channel_radius(
                next.accumulated_action, config.initial_density, config.plasma_conductivity, next.time_s));
        } else {
            next.channel_radius_m = 1e-4;
        }

        // Resistance R = L / (sigma * pi * r^2)
        double area = 3.141592653589793 * next.channel_radius_m * next.channel_radius_m;
        double resistance = config.channel_length / (config.plasma_conductivity * std::max(1e-10, area));
        double power = next.current_A * next.current_A * resistance;
        next.deposited_energy_J += power * dt;

        // Plasma channel pressure P = (gamma - 1) * E_deposited / Volume
        double volume = area * config.channel_length;
        double gamma = 1.25; // Plasma effective adiabatic index
        next.pressure_Pa = std::max(1e5, (gamma - 1.0) * next.deposited_energy_J / std::max(1e-12, volume));

        // Ideal gas plasma temperature proxy
        double n_dens = 1e26; // 10^26 m^-3
        double kB = 1.380649e-23;
        next.temperature_K = std::clamp(next.pressure_Pa / (n_dens * kB), 300.0, 50000.0);

        if (!std::isfinite(next.time_s) || !std::isfinite(next.current_A) ||
            !std::isfinite(next.channel_radius_m) || !std::isfinite(next.pressure_Pa) ||
            !std::isfinite(next.temperature_K) || !std::isfinite(next.deposited_energy_J) ||
            !std::isfinite(next.accumulated_action)) {
            throw std::overflow_error("spark step produced a non-finite result");
        }

        return next;
    }
};

} // namespace raikou
