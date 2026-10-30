#pragma once
// Rotationally and vibrationally resolved electron collision sets for molecular hydrogen.
//
// H2 has no dipole: electrons change its rotation through the quadrupole, J -> J + 2 (and back),
// with the Gerjuoy-Stein Born cross section
//   sigma(e) = (8 pi / 15) a0^2 Q^2 (J + 2)(J + 1) / ((2J + 3)(2J + 1)) sqrt(1 - u / e)
// (Q in e a0^2), isotropic in angle, and its superelastic reverse from detailed balance. An
// optional table rescales every quadrupole cross section (polarisation and the 3 eV resonance
// raise the true values far above Born). By default J -> J + 2 cross sections in the caller's
// file set it, as their ratio to Born against energy above threshold, each J taking the nearest
// J below it that the file gives: every transition carries the measured energy dependence. Odd J are ortho (nuclear-spin weight 3), even J para (weight 1); electron
// collisions never change one into the other.
//
// Vibration: each vibrational excitation in the caller's set (an EXCITATION of the target with a
// threshold between 0.3 and 4.6 eV; its state v is the one whose band origin is closest to the
// threshold) is split level by level, from each populated level (0, J) to (v, J), the threshold
// moved by the change in rotational energy. Every vibrational state reached gets its own
// rotational ladder. All levels take one Boltzmann population (vibration optionally at its own
// temperature), so every superelastic collision is in detailed balance. The file's own rotational
// excitations (thresholds below 0.3 eV) are replaced by the built-in ones.
#include "raikou/water_vapour.h"
#include <string>
#include <vector>

namespace raikou::kinetics {

struct HydrogenLevel {
    int v{0}, j{0};
    double energy_ev{0};  // above (0, 0)
    bool ortho{false};
    std::string name;     // e.g. "H2(0,1)" for v = 0, J = 1
};
// Rovibrational levels of the ground electronic state: band origins of v = 0 .. 14 and
// B_v = B_e - alpha_e (v + 1/2), D_e and H_e. Names prefix(v,J).
[[nodiscard]] std::vector<HydrogenLevel> hydrogen_levels(int v_max, int j_max, const std::string& prefix = "H2");
// The vibrational state whose band origin is closest to the energy (eV).
[[nodiscard]] int hydrogen_vibrational_state(double energy_ev);

// Gerjuoy-Stein J -> J + 2 cross section (m^2) at energy e (eV), threshold u.
[[nodiscard]] double quadrupole_cross_section(double e_ev, int j, double threshold_ev, double quadrupole_ea0sq);

struct HydrogenOptions {
    double gas_temperature_k{300};
    double vibrational_temperature_k{-1};  // -1: the gas temperature
    int j_max{-1};                         // -1: up to the level with thermal population below min_population
    double min_population{1e-7};
    double quadrupole_ea0sq{0.4836};
    double ortho_para_ratio{-1};           // -1: thermal; otherwise fixed (3 for normal hydrogen)
    std::vector<double> born_scale_energy_ev, born_scale;  // rescales every quadrupole cross section
    // With no table: the file's J -> J + 2 cross sections over the Born ones set the scale.
    bool scale_to_file_rotation{true};
    bool momentum_includes_rotation{false};  // as in WaterVapourOptions
    bool rotational_ladders{true};         // off: vibrational excitations stay as given
    std::string target{"H2"};
    double mass_amu{2.01565};
};
// `other` holds the target's momentum transfer and its other processes, for example from an
// LXCat file. GasSet::levels carries J in j and v in ka.
[[nodiscard]] GasSet hydrogen_set(const std::vector<CollisionProcess>& other, const HydrogenOptions& options = {});

} // namespace raikou::kinetics
