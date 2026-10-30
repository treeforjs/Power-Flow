#pragma once
// Rotationally and vibrationally resolved electron collision sets for water vapour.
//
// Rotation: H2O is an asymmetric top (rotational constants A, B, C) with its dipole along the b
// axis. Its levels J_KaKc come from diagonalising the rigid rotor in the symmetric-top basis
// (representation Ir: z = a, x = b, y = c); ortho levels (Ka + Kc odd) carry nuclear-spin weight
// 3, para levels weight 1, and electron collisions never change one into the other. Each dipole
// (b-type) transition gets the Born point-dipole cross section
//   sigma(e) = (8 pi / 3) a0^2 (mu / e a0)^2 (Ry / e) (S / (2J + 1)) ln[(sqrt e + sqrt(e - u)) / (sqrt e - sqrt(e - u))]
// with S the line strength over mu^2 (summed over magnetic sublevels; S over the final levels of
// a given initial level sums to 2J + 1), its forward-peaked angular distribution, and the reverse
// (superelastic) collision from detailed balance. An optional table rescales every dipole cross
// section (to measured or close-coupling values, for instance J = 0 -> 1). Populations are
// thermal at the gas temperature.
//
// Vibration and the rest: the caller supplies the momentum transfer (rotationally summed, as
// usually tabulated, or its elastic part) and any vibrational, ionization, attachment, dissociation or electronic
// processes, for example from an LXCat file. When the momentum transfer is rotationally
// summed the builder takes the thermal rotational share out (rotational collisions carry their own).
// Each vibrational state the file names gets its own rotational ladder (its own rotational
// constants, levels and dipole transitions), and the vibrational excitation is split over the
// rotational levels: from each ground level to the matching level of the vibrational state, the
// rotation carried along unchanged (symmetric vibrations) or across the asymmetry doublet
// (antisymmetric ones, where nuclear spin would otherwise change), with the threshold shifted by
// the change in rotational energy. All levels take a Boltzmann population over the whole
// rovibrational set (vibration optionally at its own temperature), so every superelastic
// collision is in detailed balance. The set comes back with every target's fraction of the gas.
#include "raikou/electron_boltzmann.h"
#include <map>
#include <string>
#include <vector>

namespace raikou::kinetics {

struct RotorLevel {
    int j{0}, ka{0}, kc{0};
    double energy_ev{0};          // above the lowest level
    bool ortho{false};
    std::string name;             // e.g. "H2O(1_01)", "H2O(11_3,8)" from J = 10
    std::vector<double> vector;   // eigenvector over K = -J .. J
};
// Rigid asymmetric-top levels up to j_max, constants in cm^-1, names prefix(J_KaKc).
[[nodiscard]] std::vector<RotorLevel> asymmetric_top_levels(double a_cm, double b_cm, double c_cm, int j_max, const std::string& prefix = "H2O");

struct DipoleTransition {
    std::size_t lower{0}, upper{0};  // indices into the levels
    double strength{0};              // line strength over mu^2, summed over sublevels
};
// The b-axis dipole transitions between the levels (Delta J = 0, +-1), strongest first is not
// guaranteed; transitions weaker than min_strength are dropped.
[[nodiscard]] std::vector<DipoleTransition> b_dipole_transitions(const std::vector<RotorLevel>& levels, double min_strength = 1e-8);

// A vibrational state's rotational ladder: the state's name as the cross-section file gives it
// (the product of a vibrational excitation of the target) and its rotational constants.
struct VibrationalLadder {
    std::string state;
    double a_cm{0}, b_cm{0}, c_cm{0};
    bool antisymmetric{false};        // B2 vibration (nu3): ortho and para labels swap
};
// The bend (010), the stretches (100) and (001), and the two stretches lumped as Itikawa and
// Mason give them (100+001, the mean constants, treated as symmetric). Effective constants
// rounded from the spectroscopic literature; rigid rotor, no centrifugal distortion.
[[nodiscard]] std::vector<VibrationalLadder> water_vibrational_ladders(const std::string& target = "H2O");

struct WaterVapourOptions {
    double gas_temperature_k{300};
    int j_max{-1};                    // -1: up to the level with thermal population below min_population
    double min_population{1e-7};      // levels and transitions from levels below this are left out
    double dipole_debye{1.8546};
    double a_cm{27.8806}, b_cm{14.5216}, c_cm{9.2778};
    double ortho_para_ratio{-1};      // -1: thermal; otherwise fixed (3 for room-temperature water)
    std::vector<double> born_scale_energy_ev, born_scale;  // rescales every dipole cross section
    // Whether the supplied momentum transfer already includes the rotational collisions (beam
    // values, rotationally summed: their thermal share is taken out) or is the elastic part alone
    // (as swarm-derived sets such as Yousfi and Benabdessadok 1996 give it: used as it is).
    bool momentum_includes_rotation{false};
    std::string target{"H2O"};
    double mass_amu{18.0106};
    // Rotational ladders in the vibrational states. Off: each vibrational state is one lumped
    // level, excited from the lumped ground state.
    bool rotational_ladders{true};
    std::vector<VibrationalLadder> ladders{water_vibrational_ladders()};  // states not listed use the ground constants
    double vibrational_temperature_k{-1};  // -1: the gas temperature
};

struct GasSet {
    std::vector<CollisionProcess> processes;
    std::map<std::string, double> fractions;
    std::vector<RotorLevel> levels;    // the rotational levels included, ground state first
    double rotational_momentum_ratio_max{0};  // largest share of the given momentum transfer the
                                              // rotational collisions carry (clipped at 1)
};
// `other` holds the target's momentum transfer (ELASTIC or EFFECTIVE on options.target) and its
// other processes; vibrational EXCITATION blocks with a product and a weight ratio become
// state-to-state pairs with thermal product populations.
[[nodiscard]] GasSet water_vapour_set(const std::vector<CollisionProcess>& other, const WaterVapourOptions& options = {});

// The Born dipole cross section of one transition (m^2) at energy e (eV).
[[nodiscard]] double born_dipole_cross_section(double e_ev, double threshold_ev, double strength_over_degeneracy, double dipole_debye);

} // namespace raikou::kinetics
