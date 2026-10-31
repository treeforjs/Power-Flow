#pragma once
// Reduced vacuum (pre-breakdown) Hall port load for calibration and long circuit runs:
//   V_mid = L_H (I_{n+1} - I_n)/dt + R_H(T_n) I_mid,   I_mid = (I_n + I_{n+1})/2
// L_H: DC magnetostatic inductance of the Hall domain at the port (2 W_B / I^2 from the sealed Hall dc_initial on the
//      45 x 45 x (4+2) mesh, 200 um cells, 9 mm box): 1.1045115428330026e-08 H (tools/hall_port_load.cpp).
// R_H: two SS304 foils in series (feed upper foil -> short -> lower foil), each 25 mm long, 5 mm x 200 um, uniform
//      current density; sigma(T) and the volumetric energy U(T) from the Hall SS304 law (evaluated_solid_1500, fixed
//      reference volume, rho_ref 8000), tabulated at 1 K by tools/ss304_law_table.cpp; adiabatic Joule heating
//      (thermal diffusion length over 200 ns is ~1 um, far below the 200 um cell). The Hall DC resistance of record is
//      216.10608772358006 V / 6012.07136047 A = 0.0359454 ohm at 300 K; this formula gives 0.0359499 ohm.
// The temperature is advanced explicitly (R at the start-of-step temperature); the Joule energy of the step enters the
// foil enthalpy exactly, so the load ledger closes: 1/2 L_H (I1^2 - I0^2) + dt R I_mid^2 = dt I_mid V_mid.
// Breakdown option: after t_breakdown the foil resistance is shunted by a plasma resistance R_plasma (parallel).
#include "c19_circuit.hpp"

namespace c19_circuit_v1 {

struct LawTable {
    std::vector<double> T, sigma, U;   // K, S/m, J/m^3 (U(300 K) = 0)
    explicit LawTable(const std::filesystem::path& csv) {
        std::ifstream in(csv); require(bool(in), "missing SS304 law table");
        std::string line; std::getline(in, line);
        while (std::getline(in, line)) {
            std::istringstream r(line); std::string c; std::vector<double> v;
            while (std::getline(r, c, ',')) v.push_back(std::stod(c));
            require(v.size() == 4, "law table row");
            T.push_back(v[0]); sigma.push_back(v[1]); U.push_back(v[3]);
        }
        require(T.size() > 2, "law table empty");
    }
    static double lerp(const std::vector<double>& x, const std::vector<double>& y, double q) {
        require(q >= x.front() && q <= x.back(), "SS304 law table: outside the evaluated 300-1500 K range");
        auto hi = std::upper_bound(x.begin(), x.end(), q); if (hi == x.end()) --hi; if (hi == x.begin()) ++hi;
        auto lo = hi - 1; const std::size_t i = lo - x.begin();
        return y[i] + (q - *lo) / (*hi - *lo) * (y[i + 1] - y[i]);
    }
    double sigma_at(double TK) const { return lerp(T, sigma, TK); }
    double T_of_U(double u) const { return lerp(U, T, u); }
};

struct ReducedHallLoad {
    double L_H = 1.1045115428330026e-08;
    double length_m = .025, width_m = .005, thickness_m = .0002, foils = 2;
    double t_breakdown_s = std::numeric_limits<double>::infinity();
    double R_plasma_ohm = 0;      // shunt after breakdown (0 = ideal short of the resistive voltage)
    const LawTable* law = nullptr;
    // state
    double I = 0, U = 0, T = 300;
    Sum2 joule_J, foil_joule_J, magnetic_J, received_J;
    // Above the evaluated law range (1500 K; Hall's thermal owner refuses there) sigma is held at sigma(1500 K) and the
    // first time this happens is recorded: a stated limitation of the reduced load, not a material model.
    double t_law_exceeded_s = std::numeric_limits<double>::quiet_NaN();
    double foil_R(double TK) const { return foils * length_m / (law->sigma_at(std::min(TK, law->T.back())) * width_m * thickness_m); }
    double R_now(double t_mid) const {
        const double Rf = foil_R(T);
        if (t_mid < t_breakdown_s) return Rf;
        return R_plasma_ohm > 0 ? Rf * R_plasma_ohm / (Rf + R_plasma_ohm) : 0.;
    }
    // Affine port relation of the step: V_mid = V0 + Z I_mid.
    std::pair<double, double> relation(double dt, double t_mid) const { return {-2 * L_H * I / dt, 2 * L_H / dt + R_now(t_mid)}; }
    double voltage(double dt, double t_mid, double I_mid) const { auto [V0, Z] = relation(dt, t_mid); return V0 + Z * I_mid; }
    void accept(double dt, double t_mid, double I_mid, double V_mid) {
        const double R = R_now(t_mid), I1 = 2 * I_mid - I, Rf = foil_R(T);
        const double joule = dt * R * I_mid * I_mid;
        // The resistive voltage R I_mid drives the foil branch; after breakdown the plasma shunt takes the rest.
        const double foil_joule = dt * (R * I_mid) * (R * I_mid) / Rf;
        U += foil_joule / (foils * length_m * width_m * thickness_m);
        if (U > law->U.back()) { if (std::isnan(t_law_exceeded_s)) t_law_exceeded_s = t_mid; T = law->T.back() + 1e-9 * (U - law->U.back()); }
        else T = law->T_of_U(U);
        joule_J.add(joule); foil_joule_J.add(foil_joule);
        magnetic_J.add(.5 * L_H * (I1 - I) * (I1 + I)); received_J.add(dt * I_mid * V_mid);
        I = I1;
    }
};

}  // namespace c19_circuit_v1
