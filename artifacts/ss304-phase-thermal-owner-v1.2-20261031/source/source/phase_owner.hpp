// phase_owner.hpp -- SS304 phase-change thermal owner v1.2: the single owner of foil metal mass (Fe, Cr, Ni, Mn,
// Si, C atoms) and condensed-phase energy (thermal and melt kinetic) for all foil slices, from 298.2 K through
// melting into strong vaporization, with melt-layer flow on the two large faces.
//
// v1.2 options (Config::enable_v1_2 switches them on over v1.1; with them off every state word is bit-identical to
// the committed v1.1 owner, tests/test_regression.cpp):
//   * depth-averaged melt layer with inertia (flow_model 2, phase_kernels.h 3d-3g) along x and z (the 25 mm height
//     through the slices of one body), Rusanov fluxes with hydrostatic reconstruction, CFL-limited sub-steps (with
//     deterministic halving on a CFL refusal), exact atom/energy transfers and a kinetic-energy word per column;
//   * exact Lorentz-work and surface-force-work ledger terms (both models; the lubrication model books them as
//     dissipated heat with book_work); the Lorentz work is the exchange term the field solver books (INTEGRATION.txt);
//   * ejected-melt reservoir: open-edge outflow and liquid lifted beyond the departure limit leave the foil as
//     ejected atoms and energy (thermal + kinetic), reported per patch;
//   * thermal sub-cycling: field steps accumulate exact Joule words and (J x B) dt in the engine (host or device
//     input); one owner step consumes the window; diagnostics bound the conductivity lag of the window.
// v1.1 description:
// v1 numerics and ledgers are unchanged (README v1); with Config defaults (all v1.1 options off) every state word is
// bit-identical to v1. v1.1 options (Config::enable_v1_1 switches all of them on):
//   * law extensions (phase_law.h): pulse-heating liquid, regular-solution activities, JANAF internal energy of
//     the vapor, Knight Knudsen layer with back pressure, thermal radiation, phase-explosion temperature;
//   * Mn, Si, C inventories (Mills Table 1 composition) with JANAF vapor pressures; preferential depletion is the
//     surface cell's own composition;
//   * melt-layer flow along x on the large faces: lubrication flux from the tangential Lorentz force (Hall J x B),
//     the gradient of recoil, capillary and normal-Lorentz pressure, and Marangoni shear; exact transfers of atoms
//     and energy between surface cells (antisymmetric, so the ledgers stay exact);
//   * surface displacement per face column (mass gained through the face plus half the thermal expansion of the
//     column from the density law) with a typed refusal when it exceeds a stated fraction of the Hall EM cell;
//   * block-Merkle state digest (state_digest.h), computed on the engine (device for the GPU engine).
// Ledgers: E + K = E0 + K0 + Joule + loads + deposits + boundary + Lorentz work + surface-force work
// - evaporation(reference + kinetic + internal) - radiation - releases - ejected; N_i = N_i,0 + deposited - evaporated
// - released - ejected, exact 128-bit identities checked every step (K: melt kinetic energy, zero unless flow_model 2).
#pragma once
#include "engine.hpp"
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <span>
#include <tuple>

namespace ss304_phase_v1_2 {

struct PhaseRefusal : std::runtime_error {
    std::string kind, quantity, identity;
    double observed = 0, bound = 0;
    PhaseRefusal(std::string k, std::string q, std::string id, double obs, double b, const std::string& detail)
        : std::runtime_error("phase-owner refusal [" + k + "] " + q + " at " + id + ": " + detail),
          kind(std::move(k)), quantity(std::move(q)), identity(std::move(id)), observed(obs), bound(b) {}
};
inline void need(bool ok, const char* what) { if (!ok) throw PhaseRefusal("invalid-configuration", "configuration", "owner", 0, 0, what); }

enum Face : int { LowY = 0, HighY = 1, LowX = 2, HighX = 3 };
inline const char* face_name(int f) { return f == 0 ? "LowY" : f == 1 ? "HighY" : f == 2 ? "LowX" : "HighX"; }

struct SliceSpec { int body = 0, k = 0; double x0 = -0.0025, y0 = -0.0027; };
struct PatchSpec { int slice = 0, face = 0; double begin = 0, end = 1; };

struct Config {
    std::vector<SliceSpec> slices;
    double width = 0.005, thickness = 0.0002, height = 0.0002;
    int nx_em = 25, ny_em = 1, rx = 1;
    std::vector<double> y_edges;
    double rho_ref = 8000.0, T0 = 300.0, start_time = 0.0;
    std::array<double, PX_NEL> mass_fraction = v1_mass_fraction;
    LawOptions law;
    int joule_mode = 1, evap_mode = 1;
    std::array<int, 4> bc_kind{0, 0, 0, 0};
    std::array<double, 4> bc_T{0, 0, 0, 0};
    std::array<int, 4> evap_face{1, 1, 1, 1};
    double phi_merge = 0.25, phi_floor = 0.05, stability_safety = 0.4;
    std::int64_t max_substeps = 1000000;
    std::vector<PatchSpec> patches;
    // v1.1
    bool flow = false;                       // melt-layer flow on the large faces
    double flow_max_fraction = 0.25;         // admitted fraction of a donor surface cell moved in one sub-step
    double geometry_refusal_fraction = 0.0;  // refuse when |displacement| > fraction x EM cell thickness (0: off)
    double back_pressure_Pa = 0.0;           // default back pressure of every patch (StepInput overrides)
    double max_reynolds = 0.0;               // refuse when the melt-layer Reynolds number exceeds this (0: off; lubrication model only)
    // v1.2
    int flow_model = 1;                      // 1 lubrication (v1.1), 2 depth-averaged layer with inertia
    bool book_work = false;                  // lubrication model: book the Lorentz and surface-force work as dissipated heat
    bool eject = false;                      // depth-averaged model: ejected-melt reservoir (open edges, departure excess)
    bool open_x = true, open_z = false;      // with eject: the x edges / the z ends of a body's face are open (outflow)
    double eject_fraction = 0.25;            // departure limit: displacement / EM cell thickness above which liquid is ejected
    double cfl = 0.25;                       // admitted (|u| + a) tau / column width
    int max_cfl_retries = 8;                 // sub-step doublings after a CFL refusal inside the step
    bool diagnostics = false;                // per-step conductivity change and window heating estimate
    // Switch every v1.1 option on (recommended configuration).
    static void enable_v1_1(Config& c) {
        c.law = LawOptions::v1_1();
        c.mass_fraction = nominal_mass_fraction;
        c.flow = true;
        c.geometry_refusal_fraction = 0.5;
    }
    // Switch every v1.2 option on over v1.1 (recommended configuration).
    static void enable_v1_2(Config& c) {
        enable_v1_1(c);
        c.flow_model = 2; c.eject = true; c.diagnostics = true;
    }
};

inline std::vector<double> graded_rows(double thickness, double dy0, double ratio, double dy_max, int ny_em) {
    std::vector<double> half{0.0};
    double dy = dy0;
    while (half.back() + dy < 0.5 * thickness - 1e-15) { half.push_back(half.back() + dy); dy = std::min(dy * ratio, dy_max); }
    std::vector<double> e = half;
    e.push_back(0.5 * thickness);
    for (auto it = half.rbegin(); it != half.rend(); ++it) e.push_back(thickness - *it);
    for (int j = 1; j < ny_em; ++j) e.push_back(thickness * j / ny_em);
    std::sort(e.begin(), e.end());
    std::vector<double> out;
    for (double x : e) if (out.empty() || x - out.back() > 1e-9 * thickness) out.push_back(x);
    out.back() = thickness;
    return out;
}

struct SurfaceLoad { int slice = 0; double x = 0, y = 0; double energy_J = 0; std::string producer, channel; std::uint64_t sequence = 0; };
struct Deposit { int slice = 0; double x = 0, y = 0; int element = 0; double atoms = 0, kinetic_J = 0; std::string producer; std::uint64_t sequence = 0; };
struct Release { int patch = 0, element = 0; std::uint64_t macro_count = 0; double weight_atoms = 0, kinetic_J = 0; };
struct StepInput {
    double dt = 0;
    std::vector<double> em_power_W;           // [S ny_em nx_em]
    std::vector<double> em_energy_J;          // alternative to em_power_W (super-cycling)
    std::vector<double> em_JB;                // optional [6 n_em]: J (A/m2) x, y, z and B (T) x, y, z per EM cell
    std::vector<double> back_pressure_Pa;     // optional [n_patches]; default Config::back_pressure_Pa
    std::vector<SurfaceLoad> loads;
    std::vector<Deposit> deposits;
    std::vector<Release> releases;
    std::string quote_sha256;
    std::string field_state_sha256;
    bool use_accumulated = false;             // v1.2: Joule energy and J x B from the accumulated field steps (dt = window)
};

using Words = std::array<px_i128, PX_NEL>;
inline Words zero_words() { Words w; w.fill(px_zero()); return w; }

struct Ledger {
    px_i128 E0{}, joule{}, loads{}, deposit_E{}, boundary{}, evap_ref{}, evap_kin{}, evap_int{}, radiation{}, release_E{};
    Words N0 = zero_words(), dep_N = zero_words(), evap_N = zero_words(), rel_N = zero_words();
    // v1.2: melt kinetic energy at the start, Lorentz work (from the field), surface-force work, dissipation booked
    // as heat (internal), ejected energy (thermal + kinetic) and its kinetic part, ejected atoms
    px_i128 K0{}, lorentz{}, surface{}, dissipation{}, eject_E{}, eject_K{};
    Words eject_N = zero_words();
};
struct PatchReport {
    std::string id; int slice = 0, face = 0; double begin = 0, end = 1;
    double area_m2 = 0, mean_T = 0, max_T = 0, liquid_fraction = 0, melt_depth_m = 0, surface_density = 0;
    std::array<double, PX_NEL> flux_atoms_m2_s{}, expected_atoms{}, surface_atom_fraction{};
    std::array<double, PX_NEL> cumulative_evaporated_atoms{}, cumulative_deposited_atoms{};
    double cumulative_evaporated_kg = 0, recession_m = 0;
    double displacement_mean_m = 0, displacement_max_abs_m = 0, recoil_pressure_Pa = 0, melt_speed_m_s = 0, released_internal_J = 0, reynolds_max = 0;
    // v1.2
    double crater_depth_m = 0, layer_ux_m_s = 0, layer_uz_m_s = 0, layer_thickness_m = 0;
    std::array<double, PX_NEL> cumulative_ejected_atoms{};
    double cumulative_ejected_J = 0, cumulative_ejected_kinetic_J = 0, cumulative_lorentz_work_J = 0;
};
struct Quote { double dt = 0; std::uint64_t revision = 0; std::vector<std::array<double, PX_NEL>> expected_atoms; std::vector<double> mean_T; std::string sha256; };
struct StepRecord {
    double dt = 0; std::int64_t n_sub = 0; double tau = 0;
    px_i128 dE{}, joule{}, loads{}, deposit_E{}, boundary{}, evap_ref{}, evap_kin{}, evap_int{}, radiation{}, release_E{};
    Words dN = zero_words(), dep_N = zero_words(), evap_N = zero_words(), rel_N = zero_words();
    double joule_requested_J = 0, max_T = 0, max_displacement_fraction = 0, phase_explosion_margin = 0, max_reynolds = 0; int merges = 0;
    std::vector<Words> patch_evap_N, patch_dep_N;
    std::string inputs_sha256, field_state_sha256;
    // v1.2
    px_i128 dK{}, lorentz{}, surface{}, dissipation{}, eject_E{}, eject_K{};
    Words eject_N = zero_words();
    std::vector<Words> patch_eject_N;
    std::vector<px_i128> patch_eject_E, patch_eject_K, patch_lorentz;   // per patch: the Lorentz work is the field-energy exchange
    std::vector<px_i128> patch_dissipation;
    double max_crater_fraction = 0, max_signal_speed = 0, sigma_change_max = 0, window_heating_K = 0;
    std::int64_t window_steps = 0, cfl_retries = 0;
    double J(px_i128 w) const { return px_to_double(w, PX_E_INV); }
};

class Owner;
class Trial {
    friend class Owner;
    const Owner* owner_ = nullptr;
    std::uint64_t generation_ = 0, base_revision_ = 0;
    int slot_ = -1;
    bool merged_ = false;
    HostElements elements_;
    std::vector<int> flags_;
    std::map<std::string, std::uint64_t> watermarks_;
    StepRecord record_;
    std::vector<px_i128> sums_;
    std::vector<double> sigma_;
    double amax_ = 0;
public:
    const StepRecord& record() const { return record_; }
};

class Owner {
    Config cfg_;
    px_law law_{};
    std::string law_sha_, config_sha_;
    px_dims dims_{};
    HostGeometry geom_;
    std::vector<PatchSpec> patch_specs_;
    std::vector<std::string> patch_ids_;
    HostElements elements_;
    std::unique_ptr<Engine> engine_;
    int slot_ = 0;
    std::uint64_t revision_ = 0, generation_ = 0;
    double time_ = 0;
    Ledger ledger_;
    std::map<std::string, std::uint64_t> watermarks_;
    std::vector<Words> patch_evap_N_, patch_dep_N_;
    std::vector<double> patch_pb_;
    double tau_max_ = 0, dy_em_ = 0;
    mutable bool eval_valid_ = false;
    mutable double eval_dt_ = -1;
    mutable EngineResult eval_;
    mutable std::vector<double> sigma_em_;
    std::vector<int> accepted_flags_;
    // v1.2
    std::vector<px_i128> sums_;              // slice sums of the accepted state (E, N, K, J)
    bool sums_valid_ = false;
    std::vector<double> sigma_cache_;        // EM-cell conductivity of the accepted state (from the step)
    bool sigma_valid_ = false;
    double acc_dt_ = 0, amax_ = 0;
    std::int64_t acc_n_ = 0;
    std::vector<Words> patch_ej_N_;
    std::vector<px_i128> patch_ej_E_, patch_ej_K_, patch_WL_;
    int per_slice() const { return dims_.nx * dims_.ny; }
    std::string slice_text(int s) const { return "foil-" + std::to_string(cfg_.slices[s].body) + "-slice-" + std::to_string(cfg_.slices[s].k); }

    void build() {
        const auto& c = cfg_;
        need(!c.slices.empty() && c.slices.size() < 100000, "1..99999 slices required");
        need(c.width > 0 && c.thickness > 0 && c.height > 0 && std::isfinite(c.width * c.thickness * c.height), "positive finite geometry");
        need(c.nx_em > 0 && c.ny_em > 0 && c.rx > 0 && c.nx_em * c.rx <= 100000, "EM/thermal column counts");
        need(c.y_edges.size() >= 2 && c.y_edges.size() <= 4097 && c.y_edges.front() == 0.0 && c.y_edges.back() == c.thickness, "thermal row edges span 0..thickness");
        for (std::size_t j = 1; j < c.y_edges.size(); ++j) need(c.y_edges[j] > c.y_edges[j - 1], "thermal row edges increase");
        need(c.rho_ref > 0 && std::isfinite(c.T0) && c.phi_merge > c.phi_floor && c.phi_floor > 0 && c.phi_merge < 1, "density/merge thresholds");
        need(c.stability_safety > 0 && c.stability_safety <= 0.5 && c.max_substeps >= 1, "stability safety in (0, 0.5]");
        need(c.joule_mode == 0 || c.joule_mode == 1, "joule mode");
        need(c.evap_mode >= 0 && c.evap_mode <= 2, "evaporation mode");
        need(c.flow_max_fraction > 0 && c.flow_max_fraction < 1 && c.geometry_refusal_fraction >= 0 && c.back_pressure_Pa >= 0, "flow / geometry / back-pressure options");
        double wsum = 0; for (double w : c.mass_fraction) { need(w >= 0, "mass fractions"); wsum += w; }
        need(std::abs(wsum - 1.0) < 1e-12, "mass fractions sum to 1");
        law_ = make_law(c.law);
        law_sha_ = ss304_phase_v1_2::law_sha256(law_);
        need(c.T0 >= law_.T_min && c.T0 <= law_.T_max, "initial temperature inside the law domain");
        if (c.flow) need(c.law.kind == 0 || (c.law.c_rho > 0 && c.law.c_mu > 0 && c.law.c_sigma >= 0), "melt flow needs density, viscosity and surface tension");
        for (int f = 0; f < 4; ++f) {
            need(c.bc_kind[f] == 0 || c.bc_kind[f] == 1, "boundary kind");
            if (c.bc_kind[f] == 1) need(c.bc_T[f] >= law_.T_min && c.bc_T[f] <= law_.T_max, "prescribed temperature inside the law domain");
        }
        const int ny = int(c.y_edges.size()) - 1, nx = c.nx_em * c.rx, S = int(c.slices.size());
        need(std::int64_t(S) * nx * ny < (1LL << 30), "cell count");
        dims_ = {};
        dims_.S = S; dims_.nx = nx; dims_.ny = ny; dims_.nx_em = c.nx_em; dims_.ny_em = c.ny_em; dims_.rx = c.rx;
        dims_.n_cells = S * nx * ny; dims_.n_em = S * c.nx_em * c.ny_em;
        dims_.dx = c.width / nx; dims_.dz = c.height;
        dims_.joule_mode = c.joule_mode; dims_.evap_mode = c.evap_mode;
        for (int f = 0; f < 4; ++f) { dims_.bc_kind[f] = c.bc_kind[f]; dims_.bc_T[f] = c.bc_T[f]; dims_.evap_face[f] = c.evap_face[f] ? 1 : 0; }
        dims_.phi_floor = c.phi_floor; dims_.phi_merge = c.phi_merge;
        dims_.flow = c.flow ? 1 : 0; dims_.flow_max_fraction = c.flow_max_fraction;
        dims_.rho0 = c.law.kind == 0 ? px_density(&law_, c.T0) : c.law.c_rho;
        dims_.n_fc = S * 2 * nx; dims_.n_ic = nx > 1 ? S * 2 * (nx - 1) : 0;
        dims_.thickness = c.thickness;
        dy_em_ = c.thickness / c.ny_em;
        // v1.2
        need(c.flow_model == 1 || c.flow_model == 2, "flow model 1 (lubrication) or 2 (depth-averaged)");
        need(!(c.flow_model == 2 && c.flow) || ny >= 2, "the depth-averaged layer needs at least two thermal rows");
        need(!c.book_work || c.flow_model == 1, "book_work applies to the lubrication model (the depth-averaged model always books the work)");
        need(!c.eject || (c.flow && c.flow_model == 2), "ejection needs the depth-averaged layer");
        need(c.cfl > 0 && c.cfl <= 0.5 && c.eject_fraction > 0 && c.max_cfl_retries >= 0 && c.max_cfl_retries <= 30, "cfl in (0, 0.5], positive ejection limit, retries 0..30");
        dims_.flow_model = c.flow_model; dims_.book_work = c.book_work ? 1 : 0; dims_.eject = c.eject ? 1 : 0;
        dims_.open_x = c.open_x ? 1 : 0; dims_.open_z = c.open_z ? 1 : 0;
        dims_.eta_limit = c.eject_fraction * dy_em_; dims_.cfl = c.cfl;
        dims_.k_cap = 3.14159265358979323846 / std::min(dims_.dx, dims_.dz);
        dims_.diag = c.diagnostics ? 1 : 0;
        dims_.cp_min = law_.kind == 0 ? px_ss_cap(law_.T_min) : std::min(law_.c_cps, law_.c_cpl);
        geom_ = {};
        geom_.row_dy.resize(ny); geom_.row_em.resize(ny); geom_.gx.resize(ny); geom_.gy.resize(std::max(ny - 1, 1));
        geom_.em_row_begin.assign(c.ny_em + 1, -1);
        for (int je = 0; je <= c.ny_em; ++je) {
            const double edge = je == c.ny_em ? c.thickness : c.thickness * je / c.ny_em;
            int found = -1;
            for (int j = 0; j <= ny; ++j) if (std::abs(c.y_edges[j] - edge) <= 1e-12 * c.thickness) { found = j; break; }
            need(found >= 0, "thermal row edges contain every EM row edge");
            geom_.em_row_begin[je] = found;
        }
        for (int j = 0; j < ny; ++j) {
            geom_.row_dy[j] = c.y_edges[j + 1] - c.y_edges[j];
            int je = 0; while (je + 1 < c.ny_em && j >= geom_.em_row_begin[je + 1]) ++je;
            geom_.row_em[j] = je;
            geom_.gx[j] = geom_.row_dy[j] * dims_.dz / dims_.dx;
            if (j + 1 < ny) geom_.gy[j] = dims_.dx * dims_.dz / (0.5 * (geom_.row_dy[j] + c.y_edges[j + 2] - c.y_edges[j + 1]));
        }
        if (c.flow && c.flow_model == 2) build_topology();
        geom_.vref.resize(dims_.n_cells); geom_.mref.resize(dims_.n_cells);
        for (int cidx = 0; cidx < dims_.n_cells; ++cidx) {
            const int j = (cidx / nx) % ny;
            geom_.vref[cidx] = dims_.dx * geom_.row_dy[j] * dims_.dz;
            geom_.mref[cidx] = c.rho_ref * geom_.vref[cidx];
        }
        patch_specs_.clear();
        std::map<std::pair<int, int>, std::vector<PatchSpec>> given;
        for (const auto& p : c.patches) { need(p.slice >= 0 && p.slice < S && p.face >= 0 && p.face < 4, "patch slice/face"); given[{p.slice, p.face}].push_back(p); }
        for (int s = 0; s < S; ++s) for (int f = 0; f < 4; ++f) {
            auto it = given.find({s, f});
            if (it == given.end()) { patch_specs_.push_back({s, f, 0.0, 1.0}); continue; }
            auto v = it->second;
            std::sort(v.begin(), v.end(), [](const PatchSpec& a, const PatchSpec& b) { return a.begin < b.begin; });
            need(v.front().begin == 0.0 && v.back().end == 1.0, "patches of a face tile [0,1]");
            for (std::size_t i = 0; i < v.size(); ++i) { need(v[i].end > v[i].begin, "patch extent"); if (i) need(v[i].begin == v[i - 1].end, "patches of a face are contiguous"); patch_specs_.push_back(v[i]); }
        }
        patch_ids_.clear();
        for (const auto& p : patch_specs_) { std::ostringstream o; o << slice_text(p.slice) << "/" << face_name(p.face) << "/" << p.begin << "-" << p.end; patch_ids_.push_back(o.str()); }
        dims_.n_patches = int(patch_specs_.size());
        patch_pb_.assign(patch_specs_.size(), c.back_pressure_Pa);
        Encoder e;
        e.text("ss304-phase-owner-v1.2-config");
        e.text(law_sha_);
        e.u64(c.slices.size());
        for (const auto& s : c.slices) { e.i64(s.body); e.i64(s.k); e.f64(s.x0); e.f64(s.y0); }
        for (double v : {c.width, c.thickness, c.height, c.rho_ref, c.T0, c.start_time, c.phi_merge, c.phi_floor, c.stability_safety,
                         c.flow_max_fraction, c.geometry_refusal_fraction, c.back_pressure_Pa, c.max_reynolds}) e.f64(v);
        e.i64(c.nx_em); e.i64(c.ny_em); e.i64(c.rx); e.i64(c.joule_mode); e.i64(c.evap_mode); e.i64(c.max_substeps); e.i64(c.flow ? 1 : 0);
        e.u64(c.y_edges.size()); for (double y : c.y_edges) e.f64(y);
        for (double w : c.mass_fraction) e.f64(w);
        for (int f = 0; f < 4; ++f) { e.i64(c.bc_kind[f]); e.f64(c.bc_T[f]); e.i64(c.evap_face[f]); }
        e.u64(patch_specs_.size()); for (const auto& p : patch_specs_) { e.i64(p.slice); e.i64(p.face); e.f64(p.begin); e.f64(p.end); }
        e.i64(c.flow_model); e.i64(c.book_work); e.i64(c.eject); e.i64(c.open_x); e.i64(c.open_z); e.f64(c.eject_fraction); e.f64(c.cfl);
        e.i64(c.max_cfl_retries); e.i64(c.diagnostics);
        config_sha_ = e.digest();
        tau_max_ = stability_bound();
    }
    double stability_bound() const {
        double cp_min, k_max;
        if (law_.kind == 0) {
            cp_min = px_ss_cap(law_.T_min);
            k_max = 0;
            for (int i = 0; i <= 4000; ++i) { const double T = law_.T_min + (law_.T_max - law_.T_min) * i / 4000.0; k_max = std::max(k_max, px_k(&law_, T)); }
        } else {
            cp_min = std::min(law_.c_cps, law_.c_cpl);
            k_max = std::max(law_.c_ks, law_.c_kl);
        }
        double dq = 0;
        if (law_.vapor_enabled && cfg_.evap_mode != 0) {
            // evaporation stiffness bound over the elements present (v1: Fe, Cr, Ni), x = 1 each
            auto q = [&](double T) {
                double s = 0; int st = 0;
                for (int el = 0; el < PX_NEL; ++el) if (cfg_.mass_fraction[el] > 0)
                    s += px_flux(&law_, el, T, 1.0, &st) * (law_.eps0[el] + 2 * boltzmann * T + (law_.gas ? px_eint(&law_, el, T) : 0.0));
                return s;
            };
            dq = 2.0 * (q(law_.T_max) - q(law_.T_max - 10.0)) / 10.0;
            if (law_.radiation) dq += 2.0 * 4.0 * law_.emissivity * PX_SIGMA_SB * law_.T_max * law_.T_max * law_.T_max;
        }
        double best = INFINITY;
        const int nx = dims_.nx, ny = dims_.ny;
        for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            const int c = j * nx + i;
            double G = 0;
            for (int f = 0; f < 4; ++f) {
                const bool interior = (f == 0 && j > 0) || (f == 1 && j + 1 < ny) || (f == 2 && i > 0) || (f == 3 && i + 1 < nx);
                if (interior) { G += k_max * (f == 0 ? geom_.gy[j - 1] : f == 1 ? geom_.gy[j] : geom_.gx[j]); continue; }
                const double A = f < 2 ? dims_.dx * dims_.dz : geom_.row_dy[j] * dims_.dz;
                const double half = f < 2 ? 0.5 * geom_.row_dy[j] : 0.5 * dims_.dx;
                if (cfg_.bc_kind[f] == 1) G += k_max * A / half;
                if (cfg_.evap_face[f]) G += dq * A;
            }
            const double C = (cfg_.evap_mode != 0 ? cfg_.phi_floor : 1.0) * geom_.mref[c] * cp_min;
            if (G > 0) best = std::min(best, C / G);
        }
        // explicit capillary leveling of the melt layer: largest rate sigma delta^3 (4/dx^2)^2 / (3 mu) with the layer
        // capped at half the thickness, the largest surface tension and the smallest viscosity of the law
        if (cfg_.flow) {
            double sig_max = 0, mu_min = INFINITY;
            if (law_.kind == 0) {
                for (int i = 0; i <= 2000; ++i) {
                    const double T = 1727.0 + (law_.T_max - 1727.0) * i / 2000.0;
                    double g, dg; px_surface_tension(&law_, T, &g, &dg);
                    sig_max = std::max(sig_max, g); mu_min = std::min(mu_min, px_viscosity(T));
                }
            } else { sig_max = law_.c_sigma; mu_min = law_.c_mu; }
            const double dmax = 0.5 * cfg_.thickness, kap = 4.0 / (dims_.dx * dims_.dx);
            const double lam = sig_max * dmax * dmax * dmax * kap * kap / (3.0 * mu_min);
            if (lam > 0) best = std::min(best, 2.0 / lam);
        }
        return cfg_.stability_safety * best;
    }
    // Depth-averaged layer topology: x interfaces inside each slice face, z interfaces between the slices k, k + 1 of
    // one body, boundary interfaces at the x edges and at the z ends of each body.
    void build_topology() {
        const int S = dims_.S, nx = dims_.nx;
        std::map<std::pair<int, int>, int> at;
        for (int s = 0; s < S; ++s) need(at.emplace(std::pair{cfg_.slices[s].body, cfg_.slices[s].k}, s).second, "distinct (body, k) per slice");
        geom_.slice_zn.assign(2 * std::size_t(S), -1);
        for (int s = 0; s < S; ++s) {
            auto dn = at.find({cfg_.slices[s].body, cfg_.slices[s].k - 1}), upi = at.find({cfg_.slices[s].body, cfg_.slices[s].k + 1});
            geom_.slice_zn[2 * s] = dn == at.end() ? -1 : dn->second;
            geom_.slice_zn[2 * s + 1] = upi == at.end() ? -1 : upi->second;
        }
        geom_.if_a.clear(); geom_.if_b.clear(); geom_.if_dir.clear();
        geom_.fc_if.assign(4 * std::size_t(dims_.n_fc), -1);
        auto add = [&](int a, int b, int dir) { geom_.if_a.push_back(a); geom_.if_b.push_back(b); geom_.if_dir.push_back(dir); return int(geom_.if_a.size()) - 1; };
        auto fci = [&](int s, int f, int i) { return (s * 2 + f) * nx + i; };
        for (int s = 0; s < S; ++s) for (int f = 0; f < 2; ++f) {
            for (int i = 0; i + 1 < nx; ++i) { const int q = add(fci(s, f, i), fci(s, f, i + 1), PX_IF_X); geom_.fc_if[4 * fci(s, f, i) + 1] = q; geom_.fc_if[4 * fci(s, f, i + 1)] = q; }
            geom_.fc_if[4 * fci(s, f, 0)] = add(fci(s, f, 0), -1, PX_IF_MX);
            geom_.fc_if[4 * fci(s, f, nx - 1) + 1] = add(fci(s, f, nx - 1), -1, PX_IF_PX);
            const int up = geom_.slice_zn[2 * s + 1], dn = geom_.slice_zn[2 * s];
            for (int i = 0; i < nx; ++i) {
                if (up >= 0) { const int q = add(fci(s, f, i), fci(up, f, i), PX_IF_Z); geom_.fc_if[4 * fci(s, f, i) + 3] = q; geom_.fc_if[4 * fci(up, f, i) + 2] = q; }
                else geom_.fc_if[4 * fci(s, f, i) + 3] = add(fci(s, f, i), -1, PX_IF_PZ);
                if (dn < 0) geom_.fc_if[4 * fci(s, f, i) + 2] = add(fci(s, f, i), -1, PX_IF_MZ);
            }
        }
        dims_.n_if = int(geom_.if_a.size());
        for (int v : geom_.fc_if) need(v >= 0, "every face column has four layer interfaces");
    }
    HostElements build_elements(const std::vector<int>& flags) const {
        HostElements el;
        const int nx = dims_.nx, ny = dims_.ny, n = dims_.n_cells;
        el.el_of.assign(4 * n, -1);
        struct Item { int patch, cell, face; double area; };
        std::vector<Item> items;
        for (int c = 0; c < n; ++c) {
            if (flags[c] & 1) continue;
            const int s = c / (nx * ny), j = (c / nx) % ny, i = c % nx;
            for (int f = 0; f < 4; ++f) {
                int cc = c; bool reaches = true;
                for (;;) { const int nb = px_neighbor(&dims_, cc, f); if (nb < 0) break; if (!(flags[nb] & 1)) { reaches = false; break; } cc = nb; }
                if (!reaches) continue;
                const double u = f < 2 ? (i + 0.5) / nx : 0.5 * (cfg_.y_edges[j] + cfg_.y_edges[j + 1]) / cfg_.thickness;
                int patch = -1;
                for (std::size_t p = 0; p < patch_specs_.size(); ++p) {
                    const auto& ps = patch_specs_[p];
                    if (ps.slice == s && ps.face == f && u >= ps.begin && (u < ps.end || (ps.end == 1.0 && u <= 1.0))) { patch = int(p); break; }
                }
                if (patch < 0) throw std::logic_error("surface element without patch");
                items.push_back({patch, c, f, f < 2 ? dims_.dx * dims_.dz : geom_.row_dy[j] * dims_.dz});
            }
        }
        std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return std::tie(a.patch, a.cell, a.face) < std::tie(b.patch, b.cell, b.face); });
        el.patch_begin.assign(patch_specs_.size() + 1, 0);
        for (std::size_t q = 0; q < items.size(); ++q) {
            el.el_cell.push_back(items[q].cell); el.el_face.push_back(items[q].face); el.el_area.push_back(items[q].area); el.el_patch.push_back(items[q].patch);
            el.el_of[4 * items[q].cell + items[q].face] = int(q);
        }
        std::size_t q = 0;
        for (std::size_t p = 0; p <= patch_specs_.size(); ++p) { while (q < items.size() && std::size_t(items[q].patch) < p) ++q; el.patch_begin[p] = int(q); }
        return el;
    }
    double law_atom_mass(int el) const { return law_.atom_mass[el]; }
    HostState initial_state() const {
        HostState s;
        const int n = dims_.n_cells;
        s.E.assign(n, px_zero()); s.N.assign(PX_NEL * n, px_zero()); s.flags.assign(n, 0); s.fcm.assign(PX_NEL * std::max(dims_.n_fc, 1), px_zero());
        s.fcV.assign(3 * std::size_t(std::max(dims_.n_fc, 1)), 0.0); s.fcK.assign(std::max(dims_.n_fc, 1), px_zero());
        const double h0 = px_h(&law_, cfg_.T0);
        for (int c = 0; c < n; ++c) {
            int bad = 0;
            double m = 0;
            for (int el = 0; el < PX_NEL; ++el) {
                s.N[PX_NEL * c + el] = px_from_double(cfg_.mass_fraction[el] * geom_.mref[c] / law_atom_mass(el), PX_N_POW2, &bad);
                m += px_to_double(s.N[PX_NEL * c + el], PX_N_INV) * law_atom_mass(el);
            }
            s.E[c] = px_from_double(m * h0, PX_E_POW2, &bad);
            need(!bad, "initial state outside the fixed-point range");
        }
        return s;
    }
    void ensure_elements(const HostElements& el) { engine_->set_elements(el); dims_.n_el = int(el.el_cell.size()); }
    const EngineResult& evaluated(double dt) const {
        if (!eval_valid_ || eval_dt_ != dt) {
            px_dims d = dims_; d.dt = dt; d.n_sub = 1; d.tau = dt;
            const_cast<Owner*>(this)->engine_->evaluate(slot_, d, patch_pb_, eval_, &sigma_em_);
            eval_valid_ = true; eval_dt_ = dt;
        }
        return eval_;
    }
    static px_i128 add(px_i128 a, px_i128 b) { int ov = 0; auto r = px_add(a, b, &ov); if (ov) throw std::overflow_error("ledger word overflow"); return r; }
    static px_i128 sub(px_i128 a, px_i128 b) { int ov = 0; auto r = px_sub(a, b, &ov); if (ov) throw std::overflow_error("ledger word overflow"); return r; }
    static px_i128 total_E(const EngineResult& r, int S) { px_i128 t = px_zero(); for (int s = 0; s < S; ++s) t = add(t, r.slice_sum[PX_SLICE_SUM * s]); return t; }
    static px_i128 total_N(const EngineResult& r, int S, int el) { px_i128 t = px_zero(); for (int s = 0; s < S; ++s) t = add(t, r.slice_sum[PX_SLICE_SUM * s + 1 + el]); return t; }
    static px_i128 sum_at(const std::vector<px_i128>& v, int S, int k) { px_i128 t = px_zero(); for (int s = 0; s < S; ++s) t = add(t, v[PX_SLICE_SUM * s + k]); return t; }
    const std::vector<px_i128>& accepted_sums() {
        if (!sums_valid_) { sums_ = evaluated(eval_valid_ ? eval_dt_ : 0.0).slice_sum; sums_valid_ = true; }
        return sums_;
    }
    void state_changed() { eval_valid_ = false; sums_valid_ = false; sigma_valid_ = false; }
    int locate(int s, double x, double y, int& face, const char* what) const {
        need(s >= 0 && s < dims_.S, "slice index");
        const auto& sl = cfg_.slices[s];
        const double tol = 32 * 2.220446049250313e-16 * std::max({std::abs(sl.x0), std::abs(sl.y0), cfg_.width, cfg_.thickness});
        if (!(std::isfinite(x) && std::isfinite(y))) throw PhaseRefusal("nonfinite-input", what, slice_text(s), 0, 0, "position not finite");
        int f = -1;
        if (std::abs(y - sl.y0) <= tol) f = 0;
        else if (std::abs(y - (sl.y0 + cfg_.thickness)) <= tol) f = 1;
        else if (std::abs(x - sl.x0) <= tol) f = 2;
        else if (std::abs(x - (sl.x0 + cfg_.width)) <= tol) f = 3;
        if (f < 0 || x < sl.x0 - tol || x > sl.x0 + cfg_.width + tol || y < sl.y0 - tol || y > sl.y0 + cfg_.thickness + tol)
            throw PhaseRefusal("off-surface", what, slice_text(s), 0, 0, "point is not on the foil cross-section boundary");
        int i, j;
        if (f < 2) { i = std::clamp(int((x - sl.x0) / dims_.dx), 0, dims_.nx - 1); j = f == 0 ? 0 : dims_.ny - 1; }
        else {
            i = f == 2 ? 0 : dims_.nx - 1;
            const double yy = std::clamp(y - sl.y0, 0.0, cfg_.thickness);
            j = std::clamp(int(std::upper_bound(cfg_.y_edges.begin(), cfg_.y_edges.end(), yy) - cfg_.y_edges.begin()) - 1, 0, dims_.ny - 1);
        }
        int c = (s * dims_.ny + j) * dims_.nx + i;
        const int inward = f == 0 ? 1 : f == 1 ? 0 : f == 2 ? 3 : 2;
        while (c >= 0 && (accepted_flags_[c] & 1)) c = px_neighbor(&dims_, c, inward);
        if (c < 0) throw PhaseRefusal("column-depleted", what, slice_text(s), 0, 0, "no material along the inward normal");
        const int e = elements_.el_of[4 * c + f];
        if (e < 0) throw std::logic_error("exposed cell without element");
        face = f;
        return e;
    }

public:
    explicit Owner(Config cfg, std::unique_ptr<Engine> engine = nullptr) : cfg_(std::move(cfg)) {
        build();
        engine_ = engine ? std::move(engine) : std::make_unique<CpuEngine>();
        engine_->configure(law_, dims_, geom_);
        const HostState s0 = initial_state();
        accepted_flags_ = s0.flags;
        elements_ = build_elements(s0.flags);
        ensure_elements(elements_);
        engine_->write_state(0, s0);
        slot_ = 0; time_ = cfg_.start_time;
        const auto& r = evaluated(0.0);
        for (int s = 0; s < dims_.S; ++s) if (r.slice_status[4 * s + 1] != PX_OK)
            throw PhaseRefusal("invalid-initial-state", "initial temperature", slice_text(s), cfg_.T0, r.slice_status[4 * s + 1], "outside the law domain or not representable");
        ledger_.E0 = total_E(r, dims_.S);
        for (int el = 0; el < PX_NEL; ++el) ledger_.N0[el] = total_N(r, dims_.S, el);
        patch_evap_N_.assign(patch_specs_.size(), zero_words());
        patch_dep_N_ = patch_evap_N_;
        patch_ej_N_ = patch_evap_N_;
        patch_ej_E_.assign(patch_specs_.size(), px_zero()); patch_ej_K_ = patch_ej_E_; patch_WL_ = patch_ej_E_;
        sums_ = r.slice_sum; sums_valid_ = true;
    }
    Owner(const Owner&) = delete; Owner& operator=(const Owner&) = delete;

    const Config& config() const { return cfg_; }
    const px_law& law() const { return law_; }
    const px_dims& dims() const { return dims_; }
    const HostGeometry& geometry() const { return geom_; }
    const std::string& law_sha256() const { return law_sha_; }
    const std::string& config_sha256() const { return config_sha_; }
    std::string engine_name() const { return engine_->name(); }
    Engine& engine() { return *engine_; }
    double time_s() const { return time_; }
    std::uint64_t revision() const { return revision_; }
    double stability_tau() const { return tau_max_; }
    std::int64_t substeps_for(double dt) const { return std::max<std::int64_t>(1, std::int64_t(std::ceil(dt / tau_max_))); }
    const Ledger& ledger() const { return ledger_; }
    const std::vector<std::string>& patch_ids() const { return patch_ids_; }
    const std::vector<PatchSpec>& patch_specs() const { return patch_specs_; }
    int n_patches() const { return int(patch_specs_.size()); }
    int em_index(int slice, int je, int ie) const { return (slice * cfg_.ny_em + je) * cfg_.nx_em + ie; }
    double em_cell_thickness() const { return dy_em_; }

    px_i128 energy_word() const { return total_E(evaluated(eval_valid_ ? eval_dt_ : 0.0), dims_.S); }
    // v1.2: melt kinetic energy of the accepted state (sum of the column words)
    px_i128 kinetic_word() const { return sum_at(evaluated(eval_valid_ ? eval_dt_ : 0.0).slice_sum, dims_.S, SL_K); }
    Words atom_words() const { const auto& r = evaluated(eval_valid_ ? eval_dt_ : 0.0); Words w; for (int el = 0; el < PX_NEL; ++el) w[el] = total_N(r, dims_.S, el); return w; }
    double energy_J() const { return px_to_double(energy_word(), PX_E_INV); }
    double max_temperature_K() const { const auto& r = evaluated(eval_valid_ ? eval_dt_ : 0.0); double m = 0; for (double t : r.slice_maxT) m = std::max(m, t); return m; }
    bool ledger_closes() const {
        px_i128 e = add(ledger_.E0, ledger_.K0);
        for (auto w : {ledger_.joule, ledger_.loads, ledger_.deposit_E, ledger_.boundary, ledger_.lorentz, ledger_.surface}) e = add(e, w);
        for (auto w : {ledger_.evap_ref, ledger_.evap_kin, ledger_.evap_int, ledger_.radiation, ledger_.release_E, ledger_.eject_E}) e = sub(e, w);
        if (!px_eq(e, add(energy_word(), kinetic_word()))) return false;
        const auto N = atom_words();
        for (int el = 0; el < PX_NEL; ++el) {
            px_i128 n = add(ledger_.N0[el], ledger_.dep_N[el]);
            n = sub(sub(sub(n, ledger_.evap_N[el]), ledger_.rel_N[el]), ledger_.eject_N[el]);
            if (!px_eq(n, N[el])) return false;
        }
        return true;
    }
    // EM-cell conductivity of the accepted state (from the step that produced it, otherwise evaluated)
    std::vector<double> sigma_em() const { if (sigma_valid_) return sigma_cache_; evaluated(eval_valid_ ? eval_dt_ : 0.0); return sigma_em_; }
    std::vector<double> temperatures() const { return const_cast<Owner*>(this)->engine_->temperatures(slot_, dims_); }
    std::vector<double> face_displacements() const { return const_cast<Owner*>(this)->engine_->face_displacements(slot_, dims_); }
    HostState state() const { return const_cast<Owner*>(this)->engine_->read_state(slot_); }
    std::vector<PatchReport> patch_reports(double dt = 0.0) const {
        const auto& r = evaluated(dt);
        std::vector<PatchReport> out;
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) {
            const double* o = r.patch_out.data() + PX_PATCH_OUT * p;
            PatchReport q;
            q.id = patch_ids_[p]; q.slice = patch_specs_[p].slice; q.face = patch_specs_[p].face; q.begin = patch_specs_[p].begin; q.end = patch_specs_[p].end;
            q.area_m2 = o[XO_AREA]; q.mean_T = o[XO_MEAN_T]; q.max_T = o[XO_MAX_T]; q.liquid_fraction = o[XO_FL]; q.melt_depth_m = o[XO_MELT]; q.surface_density = o[XO_RHO];
            q.displacement_mean_m = o[XO_ETA_MEAN]; q.displacement_max_abs_m = o[XO_ETA_MAX]; q.recoil_pressure_Pa = o[XO_RECOIL];
            q.melt_speed_m_s = o[XO_SPEED]; q.released_internal_J = o[XO_EINT]; q.reynolds_max = o[XO_RE];
            q.crater_depth_m = o[XO_ETA_NEG]; q.layer_ux_m_s = o[XO_UX]; q.layer_uz_m_s = o[XO_UZ]; q.layer_thickness_m = o[XO_H];
            for (int el = 0; el < PX_NEL; ++el) q.cumulative_ejected_atoms[el] = px_to_double(patch_ej_N_[p][el], PX_N_INV);
            q.cumulative_ejected_J = px_to_double(patch_ej_E_[p], PX_E_INV); q.cumulative_ejected_kinetic_J = px_to_double(patch_ej_K_[p], PX_E_INV);
            q.cumulative_lorentz_work_J = px_to_double(patch_WL_[p], PX_E_INV);
            double kg = 0;
            for (int el = 0; el < PX_NEL; ++el) {
                q.expected_atoms[el] = o[XO_D + el]; q.flux_atoms_m2_s[el] = o[XO_J + el]; q.surface_atom_fraction[el] = o[XO_X + el];
                q.cumulative_evaporated_atoms[el] = px_to_double(patch_evap_N_[p][el], PX_N_INV);
                q.cumulative_deposited_atoms[el] = px_to_double(patch_dep_N_[p][el], PX_N_INV);
                kg += (q.cumulative_evaporated_atoms[el] - q.cumulative_deposited_atoms[el]) * law_atom_mass(el);
            }
            q.cumulative_evaporated_kg = kg;
            const double rho = law_.kind == 0 && q.surface_density > 0 ? q.surface_density : cfg_.rho_ref;
            q.recession_m = q.area_m2 > 0 ? kg / (rho * q.area_m2) : 0;
            out.push_back(q);
        }
        return out;
    }
    void set_back_pressure(const std::vector<double>& pb) {
        need(pb.size() == patch_specs_.size(), "one back pressure per patch");
        for (double x : pb) need(std::isfinite(x) && x >= 0, "back pressure finite and nonnegative");
        patch_pb_ = pb; eval_valid_ = false;
    }
    Quote quote(double dt) const {
        if (!(std::isfinite(dt) && dt > 0)) throw PhaseRefusal("invalid-step", "dt", "quote", dt, 0, "dt must be finite and positive");
        const auto& r = evaluated(dt);
        Quote q; q.dt = dt; q.revision = revision_;
        Encoder e; e.text("ss304-phase-owner-v1.2-quote"); e.text(config_sha_); e.u64(revision_); e.f64(time_); e.f64(dt);
        for (double pb : patch_pb_) e.f64(pb);
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) {
            const double* o = r.patch_out.data() + PX_PATCH_OUT * p;
            std::array<double, PX_NEL> D; for (int el = 0; el < PX_NEL; ++el) { D[el] = o[XO_D + el]; e.f64(D[el]); }
            q.expected_atoms.push_back(D); q.mean_T.push_back(o[XO_MEAN_T]);
        }
        q.sha256 = e.digest();
        return q;
    }
    // Block-Merkle roots of the accepted state words (computed by the engine; on the GPU engine on the device).
    PxStateRoots state_roots(std::vector<std::array<px_u32, 8>>* leaves = nullptr) const { return const_cast<Owner*>(this)->engine_->state_roots(slot_, dims_, leaves); }
    std::string state_sha256() const {
        const PxStateRoots r = state_roots();
        Encoder e; e.text("ss304-phase-owner-v1.2-state-merkle"); e.text(config_sha_); e.u64(revision_); e.f64(time_);
        for (auto w : r.cells) e.u64(w);
        for (auto w : r.fcm) e.u64(w);
        if (r.has_flow) for (auto w : r.flow) e.u64(w);
        for (auto w : {ledger_.E0, ledger_.joule, ledger_.loads, ledger_.deposit_E, ledger_.boundary, ledger_.evap_ref, ledger_.evap_kin, ledger_.evap_int, ledger_.radiation, ledger_.release_E}) e.w(w);
        for (int el = 0; el < PX_NEL; ++el) for (auto w : {ledger_.N0[el], ledger_.dep_N[el], ledger_.evap_N[el], ledger_.rel_N[el]}) e.w(w);
        for (auto w : {ledger_.K0, ledger_.lorentz, ledger_.surface, ledger_.dissipation, ledger_.eject_E, ledger_.eject_K}) e.w(w);
        for (int el = 0; el < PX_NEL; ++el) e.w(ledger_.eject_N[el]);
        e.f64(acc_dt_); e.i64(acc_n_);
        e.u64(watermarks_.size()); for (const auto& [k, v] : watermarks_) { e.text(k); e.u64(v); }
        return e.digest();
    }
    // Host recomputation of every cell leaf from the downloaded state words, compared with the engine's leaves.
    bool verify_state_leaves() const {
        std::vector<std::array<px_u32, 8>> leaves;
        const PxStateRoots r = state_roots(&leaves);
        const HostState s = state();
        std::vector<std::array<px_u32, 8>> host;
        const px_digest_src g{s.E.data(), s.N.data(), s.flags.data(), s.fcm.data(), s.fcV.data(), s.fcK.data()};
        const PxStateRoots h = px_host_state_roots(g, dims_.n_cells, dims_.n_fc, r.has_flow, &host);
        return h.cells == r.cells && h.fcm == r.fcm && h.flow == r.flow && host == leaves;
    }

    template <class F> void initialize_temperature(F&& profile) {
        if (revision_ != 0) throw PhaseRefusal("invalid-initialization", "revision", "owner", double(revision_), 0, "profile must be set before the first step");
        HostState s = state();
        for (int c = 0; c < dims_.n_cells; ++c) {
            const int sl = c / per_slice(), j = (c / dims_.nx) % dims_.ny, i = c % dims_.nx;
            const double x = cfg_.slices[sl].x0 + (i + 0.5) * dims_.dx, y = cfg_.slices[sl].y0 + 0.5 * (cfg_.y_edges[j] + cfg_.y_edges[j + 1]);
            const double T = profile(sl, x, y);
            if (!(T >= law_.T_min && T <= law_.T_max)) throw PhaseRefusal("invalid-initialization", "temperature", slice_text(sl), T, law_.T_max, "initial temperature outside the law domain");
            double m = 0; for (int el = 0; el < PX_NEL; ++el) m += px_to_double(s.N[PX_NEL * c + el], PX_N_INV) * law_atom_mass(el);
            int bad = 0; s.E[c] = px_from_double(m * px_h(&law_, T), PX_E_POW2, &bad);
            if (bad) throw PhaseRefusal("fixed-point-range", "initial energy", slice_text(sl), T, 0, "outside the fixed-point range");
        }
        engine_->write_state(slot_, s);
        state_changed();
        ledger_.E0 = energy_word();
        ++generation_;
    }

    // Initial surface displacement eta(slice, face, x_centre) of the large faces (benchmarks), before the first step:
    // the displacement mass eta rho0 dx dz is added to (or taken from) the surface cell with that cell's composition
    // and specific energy, and recorded as face mass so that the displacement bookkeeping stays consistent.
    template <class F> void initialize_displacement(F&& eta) {
        if (revision_ != 0) throw PhaseRefusal("invalid-initialization", "revision", "owner", double(revision_), 0, "displacement must be set before the first step");
        need(dims_.rho0 > 0, "reference density required");
        HostState s = state();
        for (int sl = 0; sl < dims_.S; ++sl) for (int f = 0; f < 2; ++f) for (int i = 0; i < dims_.nx; ++i) {
            const double x = cfg_.slices[sl].x0 + (i + 0.5) * dims_.dx;
            const double e = eta(sl, f, x);
            if (e == 0.0) continue;
            const int j = f == 0 ? 0 : dims_.ny - 1, c = (sl * dims_.ny + j) * dims_.nx + i;
            double m = 0; for (int el = 0; el < PX_NEL; ++el) m += px_to_double(s.N[PX_NEL * c + el], PX_N_INV) * law_atom_mass(el);
            const double ratio = e * dims_.rho0 * dims_.dx * dims_.dz / m;
            if (!(ratio > -0.5 && ratio < 0.5)) throw PhaseRefusal("invalid-initialization", "displacement", slice_text(sl), e, 0, "displacement mass exceeds half the surface cell");
            int bad = 0, ov = 0;
            const int fc = (sl * 2 + f) * dims_.nx + i;
            for (int el = 0; el < PX_NEL; ++el) {
                const px_i128 dn = px_from_double(px_to_double(s.N[PX_NEL * c + el], PX_N_INV) * ratio, PX_N_POW2, &bad);
                s.N[PX_NEL * c + el] = px_add(s.N[PX_NEL * c + el], dn, &ov);
                s.fcm[PX_NEL * fc + el] = px_add(s.fcm[PX_NEL * fc + el], dn, &ov);
            }
            s.E[c] = px_add(s.E[c], px_from_double(px_to_double(s.E[c], PX_E_INV) * ratio, PX_E_POW2, &bad), &ov);
            if (bad || ov) throw PhaseRefusal("fixed-point-range", "displacement", slice_text(sl), e, 0, "outside the fixed-point range");
        }
        engine_->write_state(slot_, s);
        state_changed();
        ledger_.E0 = energy_word();
        const auto N = atom_words(); for (int el = 0; el < PX_NEL; ++el) ledger_.N0[el] = N[el];
        ++generation_;
    }

    // v1.2 (benchmarks): initial layer velocity (u_x, u_z) = vel(slice, face, x) of the depth-averaged layer, before the
    // first step. The layer is found as in the step (fully liquid cells of the face's half column); K0 is booked.
    template <class F> void initialize_flow(F&& vel) {
        if (revision_ != 0) throw PhaseRefusal("invalid-initialization", "revision", "owner", double(revision_), 0, "flow must be set before the first step");
        need(cfg_.flow && cfg_.flow_model == 2, "initial flow needs the depth-averaged layer");
        HostState s = state();
        const std::vector<double> T = temperatures();
        px_i128 K = px_zero();
        for (int sl = 0; sl < dims_.S; ++sl) for (int f = 0; f < 2; ++f) for (int i = 0; i < dims_.nx; ++i) {
            const int fc = (sl * 2 + f) * dims_.nx + i;
            int c = -1;
            for (int jj = 0; jj < dims_.ny; ++jj) { const int j = f == 0 ? jj : dims_.ny - 1 - jj, cc = (sl * dims_.ny + j) * dims_.nx + i; if (!(s.flags[cc] & 1)) { c = cc; break; } }
            double m = 0;
            if (c >= 0 && px_in_half(&dims_, c, f)) for (;;) {
                if ((s.flags[c] & 1) || px_fl(&law_, T[c]) < 1.0) break;
                for (int el = 0; el < PX_NEL; ++el) m += px_to_double(s.N[PX_NEL * c + el], PX_N_INV) * law_atom_mass(el);
                const int nxt = c + px_layer_step(&dims_, f);
                if (!px_in_half(&dims_, nxt, f)) break;
                c = nxt;
            }
            const auto u = vel(sl, f, cfg_.slices[sl].x0 + (i + 0.5) * dims_.dx);
            if (!(std::isfinite(u[0]) && std::isfinite(u[1]))) throw PhaseRefusal("invalid-initialization", "velocity", slice_text(sl), 0, 0, "velocity not finite");
            if (m <= 0) { s.fcV[3 * fc] = s.fcV[3 * fc + 1] = s.fcV[3 * fc + 2] = 0; s.fcK[fc] = px_zero(); continue; }
            s.fcV[3 * fc] = m * u[0]; s.fcV[3 * fc + 1] = m * u[1]; s.fcV[3 * fc + 2] = m;
            int bad = 0;
            s.fcK[fc] = px_from_double(0.5 * m * (u[0] * u[0] + u[1] * u[1]), PX_E_POW2, &bad);
            if (bad) throw PhaseRefusal("fixed-point-range", "kinetic energy", slice_text(sl), 0, 0, "outside the fixed-point range");
            K = add(K, s.fcK[fc]);
        }
        engine_->write_state(slot_, s);
        state_changed();
        ledger_.K0 = K;
        ++generation_;
    }

    // v1.2 sub-cycling: one field step of length dt into the accumulators (Joule power W per EM cell, optional J and B
    // [6 n_em]). The owner step with StepInput::use_accumulated consumes the window (dt = sum of the field steps).
    void accumulate(const std::vector<double>& em_power_W, const std::vector<double>& em_JB, double dt) {
        if (!(std::isfinite(dt) && dt > 0)) throw PhaseRefusal("invalid-step", "dt", "accumulate", dt, 0, "field step must be finite and positive");
        if (int(em_power_W.size()) != dims_.n_em || !(em_JB.empty() || em_JB.size() == 6 * std::size_t(dims_.n_em)))
            throw PhaseRefusal("invalid-input", "accumulated vectors", "owner", double(em_power_W.size()), dims_.n_em, "one power (and optionally six J, B values) per EM cell");
        engine_->accumulate(em_power_W, em_JB, dt);
        acc_dt_ += dt; ++acc_n_;
    }
    // Same with the field data already on the device of the GPU engine's session (read in stream order).
    void accumulate_device(const double* em_power_W, const double* em_JB, double dt) {
        if (!(std::isfinite(dt) && dt > 0)) throw PhaseRefusal("invalid-step", "dt", "accumulate", dt, 0, "field step must be finite and positive");
        engine_->accumulate_device(em_power_W, em_JB, dt);
        acc_dt_ += dt; ++acc_n_;
    }
    double window_dt() const { return acc_dt_; }
    std::int64_t window_steps() const { return acc_n_; }

    // -------------------------------------------------------------------------------------------- step
    Trial prepare(const StepInput& in) {
        const bool acc = in.use_accumulated;
        if (acc && acc_n_ == 0) throw PhaseRefusal("invalid-input", "accumulated window", "owner", 0, 0, "no field step was accumulated");
        if (acc && !(in.dt == 0 || in.dt == acc_dt_)) throw PhaseRefusal("invalid-step", "dt", "owner", in.dt, acc_dt_, "an accumulated step takes the window length (dt 0 or the window)");
        const double dt = acc ? acc_dt_ : in.dt;
        if (!(std::isfinite(dt) && dt > 0 && std::isfinite(time_ + dt) && time_ + dt > time_))
            throw PhaseRefusal("invalid-step", "dt", "owner", dt, 0, "dt must be finite, positive and advance the clock");
        std::int64_t n_sub = substeps_for(dt);
        if (cfg_.flow && cfg_.flow_model == 2 && amax_ > 0) {
            // CFL estimate from the largest signal speed of the last accepted step (25 % margin)
            const double ncfl = std::ceil(1.25 * dt * amax_ / (cfg_.cfl * std::min(dims_.dx, dims_.dz)));
            if (ncfl > double(n_sub)) n_sub = ncfl < 9e18 ? std::int64_t(ncfl) : std::numeric_limits<std::int64_t>::max();
        }
        if (n_sub > cfg_.max_substeps || n_sub > 0x7fffffff)
            throw PhaseRefusal("substep-budget", "sub-steps", "owner", double(n_sub), double(cfg_.max_substeps), "the step needs more explicit sub-steps than the configured budget");
        const bool by_energy = !in.em_energy_J.empty();
        const auto& jv = by_energy ? in.em_energy_J : in.em_power_W;
        if (acc) {
            if (!in.em_energy_J.empty() || !in.em_power_W.empty() || !in.em_JB.empty())
                throw PhaseRefusal("invalid-input", "Joule vector", "owner", 0, 0, "an accumulated step takes no Joule or J, B vector");
        } else if (int(jv.size()) != dims_.n_em || (by_energy && !in.em_power_W.empty()))
            throw PhaseRefusal("invalid-input", "Joule vector", "owner", double(jv.size()), dims_.n_em, "one Joule power or one Joule energy (not both) per EM foil cell");
        StepData data;
        data.use_acc = acc;
        double joule_req = 0;
        px_i128 joule_exact = px_zero();
        if (!acc) data.em_energy.resize(dims_.n_em);
        if (!acc) for (int e = 0; e < dims_.n_em; ++e) {
            const double P = jv[e];
            if (!(std::isfinite(P) && P >= 0)) throw PhaseRefusal("invalid-input", by_energy ? "Joule energy J" : "Joule power W", "em-cell-" + std::to_string(e), P, 0, "Joule input must be finite and nonnegative");
            data.em_energy[e] = by_energy ? P : P * dt;
            joule_req += data.em_energy[e];
            int bad = 0;
            joule_exact = add(joule_exact, px_from_double(data.em_energy[e], PX_E_POW2, &bad));
            if (bad) throw PhaseRefusal("fixed-point-range", "Joule energy", "em-cell-" + std::to_string(e), data.em_energy[e], 0, "outside the 128-bit energy range");
        }
        if (acc) {}
        else if (in.em_JB.empty()) data.em_JB.assign(6 * std::size_t(dims_.n_em), 0.0);
        else {
            if (in.em_JB.size() != 6 * std::size_t(dims_.n_em)) throw PhaseRefusal("invalid-input", "J and B vector", "owner", double(in.em_JB.size()), 6.0 * dims_.n_em, "six values (J, B) per EM cell");
            for (double v : in.em_JB) if (!std::isfinite(v)) throw PhaseRefusal("nonfinite-input", "J and B", "owner", v, 0, "J and B must be finite");
            data.em_JB = in.em_JB;
        }
        if (!in.back_pressure_Pa.empty()) {
            if (in.back_pressure_Pa.size() != patch_specs_.size()) throw PhaseRefusal("invalid-input", "back pressure", "owner", double(in.back_pressure_Pa.size()), double(patch_specs_.size()), "one back pressure per patch");
            for (double v : in.back_pressure_Pa) if (!(std::isfinite(v) && v >= 0)) throw PhaseRefusal("invalid-input", "back pressure Pa", "owner", v, 0, "finite and nonnegative");
            if (in.back_pressure_Pa != patch_pb_) { patch_pb_ = in.back_pressure_Pa; eval_valid_ = false; }
        }
        data.patch_pb = patch_pb_;
        auto watermarks = watermarks_;
        std::set<std::pair<std::string, std::uint64_t>> seen;
        auto admit = [&](const std::string& producer, std::uint64_t seq, const char* what) {
            if (producer.empty()) throw PhaseRefusal("load-identity", what, "owner", 0, 0, "producer required");
            auto it = watermarks_.find(producer);
            const std::uint64_t prior = it == watermarks_.end() ? 0 : it->second;
            if (seq <= prior || !seen.emplace(producer, seq).second)
                throw PhaseRefusal("load-identity", what, producer + "#" + std::to_string(seq), double(seq), double(prior), "duplicate or stale producer sequence");
            watermarks[producer] = std::max(watermarks[producer], seq);
        };
        std::map<int, std::array<px_i128, PX_SC>> per_el;
        auto entry = [&](int e) -> std::array<px_i128, PX_SC>& { auto it = per_el.find(e); if (it == per_el.end()) { std::array<px_i128, PX_SC> z; z.fill(px_zero()); it = per_el.emplace(e, z).first; } return it->second; };
        px_i128 loads_total = px_zero(), depE_total = px_zero();
        Words depN_total = zero_words();
        std::vector<Words> patch_dep(patch_specs_.size(), zero_words());
        for (const auto& l : in.loads) {
            admit(l.producer, l.sequence, "surface load");
            if (l.channel.empty()) throw PhaseRefusal("load-identity", "surface load channel", l.producer, 0, 0, "channel required");
            if (!std::isfinite(l.energy_J)) throw PhaseRefusal("nonfinite-input", "surface load energy", l.producer, l.energy_J, 0, "load energy not finite");
            int f; const int e = locate(l.slice, l.x, l.y, f, "surface load");
            int bad = 0; const px_i128 q = px_from_double(l.energy_J, PX_E_POW2, &bad);
            if (bad) throw PhaseRefusal("fixed-point-range", "surface load energy", l.producer, l.energy_J, 0, "outside the 128-bit energy range");
            auto& a = entry(e);
            a[0] = add(a[0], q); loads_total = add(loads_total, q);
        }
        for (const auto& dpt : in.deposits) {
            admit(dpt.producer, dpt.sequence, "metal deposit");
            if (!(dpt.element >= 0 && dpt.element < PX_NEL && std::isfinite(dpt.atoms) && dpt.atoms >= 0 && std::isfinite(dpt.kinetic_J) && dpt.kinetic_J >= 0))
                throw PhaseRefusal("invalid-input", "metal deposit", dpt.producer, dpt.atoms, 0, "element, atoms and kinetic energy must be valid");
            int f; const int e = locate(dpt.slice, dpt.x, dpt.y, f, "metal deposit");
            int bad = 0;
            const px_i128 qn = px_from_double(dpt.atoms, PX_N_POW2, &bad);
            const double a_q = px_to_double(qn, PX_N_INV);
            // a returning atom brings back its reference energy and its kinetic energy (no internal excitation)
            const px_i128 qe = add(px_from_double(a_q * law_.eps0[dpt.element], PX_E_POW2, &bad), px_from_double(dpt.kinetic_J, PX_E_POW2, &bad));
            if (bad) throw PhaseRefusal("fixed-point-range", "metal deposit", dpt.producer, dpt.atoms, 0, "outside the fixed-point range");
            auto& a = entry(e);
            a[1 + dpt.element] = add(a[1 + dpt.element], qn); a[1 + PX_NEL] = add(a[1 + PX_NEL], qe);
            depN_total[dpt.element] = add(depN_total[dpt.element], qn); depE_total = add(depE_total, qe);
            const int p = elements_.el_patch[e];
            patch_dep[p][dpt.element] = add(patch_dep[p][dpt.element], qn);
        }
        for (const auto& [e, v] : per_el) { data.sc_el.push_back(e); for (auto w : v) data.sc_val.push_back(w); }
        data.rel_total_N.assign(PX_NEL * patch_specs_.size(), px_zero());
        data.rel_kin_J.assign(PX_NEL * patch_specs_.size(), 0.0);
        Words relN_total = zero_words();
        if (cfg_.evap_mode == 2) {
            const Quote q = quote(dt);
            if (!in.releases.empty() && in.quote_sha256 != q.sha256) throw PhaseRefusal("release-binding", "quote digest", "owner", 0, 0, "releases were not drawn from this owner's quote for this state and dt");
            std::set<std::pair<int, int>> once;
            for (const auto& r : in.releases) {
                if (!(r.patch >= 0 && r.patch < int(patch_specs_.size()) && r.element >= 0 && r.element < PX_NEL && once.emplace(r.patch, r.element).second))
                    throw PhaseRefusal("release-binding", "release identity", "patch-" + std::to_string(r.patch), 0, 0, "one release per patch and element");
                const double w = r.weight_atoms, n = double(r.macro_count), D = q.expected_atoms[r.patch][r.element];
                if (!(std::isfinite(w) && w >= 1 && w == std::floor(w) && n * w < 9007199254740992.0 && std::isfinite(r.kinetic_J) && r.kinetic_J >= 0))
                    throw PhaseRefusal("release-binding", "release weight", patch_ids_[r.patch], w, 0, "integer weight, exact atom count and nonnegative kinetic energy required");
                const double lo = std::floor(D / w), hi = std::ceil(D / w);
                if (!(n == lo || n == hi)) throw PhaseRefusal("release-binding", "macro count", patch_ids_[r.patch], n, D / w, "n must be floor or ceil of D/w0 (stochastic rounding of the quote)");
                if (r.macro_count == 0 && r.kinetic_J != 0) throw PhaseRefusal("release-binding", "kinetic energy", patch_ids_[r.patch], r.kinetic_J, 0, "no atoms released");
                int bad = 0;
                const px_i128 R = px_from_double(n * w, PX_N_POW2, &bad);
                data.rel_total_N[PX_NEL * r.patch + r.element] = R;
                data.rel_kin_J[PX_NEL * r.patch + r.element] = r.kinetic_J;
                relN_total[r.element] = add(relN_total[r.element], R);
            }
        } else if (!in.releases.empty()) {
            throw PhaseRefusal("release-binding", "releases", "owner", 0, 0, "releases are accepted only in quantized mode");
        }
        px_dims d = dims_;
        d.n_el = int(elements_.el_cell.size());
        d.dt = dt; d.n_sub = int(n_sub); d.tau = dt / double(n_sub); d.sub = 0;
        const std::vector<px_i128> before = accepted_sums();
        const px_i128 E_before = sum_at(before, dims_.S, 0), K_before = sum_at(before, dims_.S, SL_K);
        Words N_before; for (int el = 0; el < PX_NEL; ++el) N_before[el] = sum_at(before, dims_.S, 1 + el);
        EngineResult r;
        std::int64_t retries = 0;
        for (;;) {
            engine_->step(slot_, d, data, r);
            bool cfl = false;
            for (int s = 0; s < dims_.S; ++s) if (r.slice_status[4 * s] >= 0 && (r.slice_status[4 * s + 1] == PX_ERR_CFL || r.slice_status[4 * s + 1] == PX_ERR_FLOW)) cfl = true;
            if (!cfl || !(cfg_.flow && cfg_.flow_model == 2) || retries >= cfg_.max_cfl_retries || 2 * std::int64_t(d.n_sub) > std::min<std::int64_t>(cfg_.max_substeps, 0x7fffffff)) break;
            // a melt-layer wave crossed more than the admitted fraction of a column: halve the sub-step (deterministic)
            ++retries; d.n_sub *= 2; d.tau = dt / double(d.n_sub);
        }
        n_sub = d.n_sub;
        eval_valid_ = false;
        if (acc) for (int s = 0; s < dims_.S; ++s) if (r.slice_diag[4 * s + 2] != 0.0)
            throw PhaseRefusal("invalid-input", "accumulated Joule power or J, B", slice_text(s), 0, 0, "an accumulated field step was not finite or negative");
        for (int s = 0; s < dims_.S; ++s) {
            const int cell = r.slice_status[4 * s], code = r.slice_status[4 * s + 1];
            if (cell < 0 || code == PX_OK) continue;
            const int local = cell % per_slice();
            const std::string id = slice_text(s) + "/row-" + std::to_string(local / dims_.nx) + "/col-" + std::to_string(local % dims_.nx);
            switch (code) {
            case PX_ERR_T_HIGH: throw PhaseRefusal("temperature-above-domain", "specific energy", id, r.slice_maxT[s], law_.T_max, "the step would heat a cell above the law domain (T_max); no state change");
            case PX_ERR_T_LOW: throw PhaseRefusal("temperature-below-domain", "specific energy", id, 0, law_.T_min, "the step would cool a cell below the law domain (T_min); no state change");
            case PX_ERR_DEPLETED: throw PhaseRefusal("surface-depletion-within-step", "cell mass fraction", id, 0, cfg_.phi_floor, "a cell lost more than its admitted mass within one step; use a shorter step");
            case PX_ERR_NEGATIVE_INV: throw PhaseRefusal("release-exceeds-inventory", "element inventory", id, 0, 0, "a release or evaporation exceeds the element inventory of the surface cell");
            case PX_ERR_FIXED_RANGE: throw PhaseRefusal("fixed-point-range", "energy/atom word", id, 0, 0, "a value left the 128-bit fixed-point range");
            case PX_ERR_FLOW: throw PhaseRefusal("melt-flow-step", "donor mass fraction", id, 0, cfg_.flow_max_fraction, "melt flow would move more than the admitted fraction of a surface cell in one sub-step; use a shorter step");
            case PX_ERR_PHASE_EXPLOSION: throw PhaseRefusal("phase-explosion", "temperature", id, r.slice_maxT[s], law_.T_explosion, "a cell reached the phase-explosion temperature (0.9 T_c); no state change");
            case PX_ERR_CFL: throw PhaseRefusal("melt-flow-cfl", "melt-layer signal speed", id, double(retries), cfg_.cfl, "a melt-layer wave would cross more than the admitted column fraction even after the sub-step halvings; use a shorter step");
            default: throw PhaseRefusal("nonfinite-state", "cell state", id, 0, 0, "nonfinite intermediate value");
            }
        }
        StepRecord rec;
        rec.dt = dt; rec.n_sub = n_sub; rec.tau = d.tau; rec.joule_requested_J = joule_req; rec.field_state_sha256 = in.field_state_sha256;
        rec.joule = acc ? sum_at(r.slice_sum, dims_.S, SL_J) : joule_exact;
        if (acc) rec.joule_requested_J = px_to_double(rec.joule, PX_E_INV);
        rec.loads = loads_total; rec.deposit_E = depE_total; rec.dep_N = depN_total;
        rec.window_steps = acc ? acc_n_ : 0; rec.cfl_retries = retries;
        rec.patch_evap_N.assign(patch_specs_.size(), zero_words());
        rec.patch_dep_N = patch_dep;
        rec.patch_eject_N.assign(patch_specs_.size(), zero_words());
        rec.patch_eject_E.assign(patch_specs_.size(), px_zero()); rec.patch_eject_K = rec.patch_eject_E; rec.patch_lorentz = rec.patch_eject_E;
        rec.patch_dissipation = rec.patch_eject_E;
        double max_eta = 0, max_neg = 0;
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) {
            const px_i128* ps = r.patch_sum.data() + PX_PATCH_SUM * p;
            for (int el = 0; el < PX_NEL; ++el) { rec.evap_N[el] = add(rec.evap_N[el], ps[XS_EVN + el]); rec.patch_evap_N[p][el] = add(rec.patch_evap_N[p][el], ps[XS_EVN + el]); }
            rec.evap_ref = add(rec.evap_ref, ps[XS_EREF]); rec.evap_kin = add(rec.evap_kin, ps[XS_EKIN]); rec.evap_int = add(rec.evap_int, ps[XS_EINT]);
            rec.boundary = add(rec.boundary, ps[XS_BND]); rec.radiation = add(rec.radiation, ps[XS_RAD]);
            rec.release_E = sub(rec.release_E, ps[XS_REL]);
            for (int el = 0; el < PX_NEL; ++el) rec.patch_evap_N[p][el] = add(rec.patch_evap_N[p][el], data.rel_total_N[PX_NEL * p + el]);
            max_eta = std::max(max_eta, r.patch_out[PX_PATCH_OUT * p + XO_ETA_MAX]);
            max_neg = std::max(max_neg, r.patch_out[PX_PATCH_OUT * p + XO_ETA_NEG]);
            rec.max_reynolds = std::max(rec.max_reynolds, r.patch_out[PX_PATCH_OUT * p + XO_RE]);
            // v1.2 terms
            rec.lorentz = add(rec.lorentz, ps[XS_WL]); rec.surface = add(rec.surface, ps[XS_WS]); rec.dissipation = add(rec.dissipation, ps[XS_Q]);
            rec.eject_E = add(rec.eject_E, ps[XS_EJE]); rec.eject_K = add(rec.eject_K, ps[XS_EJK]);
            for (int el = 0; el < PX_NEL; ++el) { rec.eject_N[el] = add(rec.eject_N[el], ps[XS_EJN + el]); rec.patch_eject_N[p][el] = ps[XS_EJN + el]; }
            rec.patch_eject_E[p] = ps[XS_EJE]; rec.patch_eject_K[p] = ps[XS_EJK]; rec.patch_lorentz[p] = ps[XS_WL]; rec.patch_dissipation[p] = ps[XS_Q];
        }
        rec.rel_N = relN_total;
        const px_i128 E_after = total_E(r, dims_.S);
        rec.dE = sub(E_after, E_before);
        rec.dK = sub(sum_at(r.slice_sum, dims_.S, SL_K), K_before);
        px_i128 expect = add(add(add(rec.joule, rec.loads), rec.deposit_E), rec.boundary);
        expect = add(add(expect, rec.lorentz), rec.surface);
        for (auto w : {rec.evap_ref, rec.evap_kin, rec.evap_int, rec.radiation, rec.release_E, rec.eject_E}) expect = sub(expect, w);
        if (!px_eq(expect, add(rec.dE, rec.dK))) throw std::logic_error("exact energy ledger failed to close (implementation defect)");
        for (int el = 0; el < PX_NEL; ++el) {
            rec.dN[el] = sub(total_N(r, dims_.S, el), N_before[el]);
            const px_i128 en = sub(sub(sub(rec.dep_N[el], rec.evap_N[el]), rec.rel_N[el]), rec.eject_N[el]);
            if (!px_eq(en, rec.dN[el])) throw std::logic_error("exact atom ledger failed to close (implementation defect)");
        }
        for (double t : r.slice_maxT) rec.max_T = std::max(rec.max_T, t);
        rec.max_displacement_fraction = max_eta / dy_em_;
        rec.max_crater_fraction = max_neg / dy_em_;
        for (int s = 0; s < dims_.S; ++s) {
            rec.sigma_change_max = std::max(rec.sigma_change_max, r.slice_diag[4 * s]);
            rec.window_heating_K = std::max(rec.window_heating_K, r.slice_diag[4 * s + 1]);
            rec.max_signal_speed = std::max(rec.max_signal_speed, r.slice_diag[4 * s + 3]);
        }
        rec.phase_explosion_margin = law_.T_explosion > 0 ? rec.max_T / law_.T_explosion : 0.0;
        if (cfg_.geometry_refusal_fraction > 0 && rec.max_displacement_fraction > cfg_.geometry_refusal_fraction)
            throw PhaseRefusal("geometry-departure", "surface displacement / EM cell thickness", "owner", rec.max_displacement_fraction, cfg_.geometry_refusal_fraction,
                               "the foil surface would move further from the fixed Hall mesh than admitted; no state change");
        if (cfg_.max_reynolds > 0 && cfg_.flow_model == 1 && rec.max_reynolds > cfg_.max_reynolds)
            throw PhaseRefusal("melt-flow-reynolds", "melt-layer Reynolds number", "owner", rec.max_reynolds, cfg_.max_reynolds,
                               "the lubrication (Stokes) melt-flow model is outside its admitted Reynolds range; no state change");
        Trial trial;
        int candidates = 0; for (int s = 0; s < dims_.S; ++s) candidates += r.slice_status[4 * s + 2];
        if (candidates > 0) {
            trial.elements_ = elements_;
            HostState st = engine_->read_state(r.trial_slot);
            for (int c = 0; c < dims_.n_cells; ++c) {
                if (st.flags[c] & 1) continue;
                double m = 0; for (int el = 0; el < PX_NEL; ++el) m += px_to_double(st.N[PX_NEL * c + el], PX_N_INV) * law_atom_mass(el);
                bool exposed = false; for (int f = 0; f < 4; ++f) if (trial.elements_.el_of[4 * c + f] >= 0) exposed = true;
                if (!exposed || !(m < cfg_.phi_merge * geom_.mref[c])) continue;
                int target = -1;
                for (int f = 0; f < 4 && target < 0; ++f) {
                    if (trial.elements_.el_of[4 * c + f] < 0) continue;
                    const int nb = px_neighbor(&dims_, c, f == 0 ? 1 : f == 1 ? 0 : f == 2 ? 3 : 2);
                    if (nb >= 0 && !(st.flags[nb] & 1)) target = nb;
                }
                if (target < 0) throw PhaseRefusal("column-depleted", "surface cell mass", slice_text(c / per_slice()), m, cfg_.phi_merge * geom_.mref[c], "a depleted surface cell has no inward neighbour (burn-through)");
                int ov = 0;
                st.E[target] = px_add(st.E[target], st.E[c], &ov); st.E[c] = px_zero();
                for (int el = 0; el < PX_NEL; ++el) { st.N[PX_NEL * target + el] = px_add(st.N[PX_NEL * target + el], st.N[PX_NEL * c + el], &ov); st.N[PX_NEL * c + el] = px_zero(); }
                if (ov) throw std::overflow_error("merge overflow");
                st.flags[c] |= 1;
                trial.elements_ = build_elements(st.flags);
                ++rec.merges;
            }
            engine_->write_state(r.trial_slot, st);
            if (rec.merges > 0) { trial.merged_ = true; trial.flags_ = st.flags; }
        }
        Encoder ie; ie.text("ss304-phase-owner-v1.2-step-inputs"); ie.f64(dt); ie.text(in.field_state_sha256);
        if (acc) {
            // accumulated window: field-step count, window length and the exact Joule words per slice
            ie.text("accumulated"); ie.i64(acc_n_);
            for (int s = 0; s < dims_.S; ++s) ie.w(r.slice_sum[PX_SLICE_SUM * s + SL_J]);
        } else {
            ie.i64(by_energy ? 1 : 0);
            for (double P : jv) ie.f64(P);
            for (double v : data.em_JB) ie.f64(v);
        }
        for (double v : data.patch_pb) ie.f64(v);
        for (std::size_t q = 0; q < data.sc_el.size(); ++q) { ie.i64(data.sc_el[q]); for (int k = 0; k < PX_SC; ++k) ie.w(data.sc_val[PX_SC * q + k]); }
        for (std::size_t q = 0; q < data.rel_total_N.size(); ++q) { ie.w(data.rel_total_N[q]); ie.f64(data.rel_kin_J[q]); }
        rec.inputs_sha256 = ie.digest();
        trial.owner_ = this; trial.generation_ = generation_; trial.base_revision_ = revision_;
        trial.slot_ = r.trial_slot; trial.watermarks_ = std::move(watermarks); trial.record_ = std::move(rec);
        trial.sums_ = r.slice_sum; trial.sigma_ = r.sigma_em; trial.amax_ = trial.record_.max_signal_speed;
        return trial;
    }
    bool preflight(const Trial& t) const noexcept { return t.owner_ == this && t.generation_ == generation_ && t.base_revision_ == revision_ && t.slot_ >= 0 && t.slot_ != slot_; }
    bool commit(Trial&& t) {
        if (!preflight(t)) return false;
        const auto& r = t.record_;
        slot_ = t.slot_;
        if (t.merged_) { elements_ = std::move(t.elements_); ensure_elements(elements_); accepted_flags_ = std::move(t.flags_); }
        ledger_.joule = add(ledger_.joule, r.joule); ledger_.loads = add(ledger_.loads, r.loads); ledger_.deposit_E = add(ledger_.deposit_E, r.deposit_E);
        ledger_.boundary = add(ledger_.boundary, r.boundary); ledger_.evap_ref = add(ledger_.evap_ref, r.evap_ref); ledger_.evap_kin = add(ledger_.evap_kin, r.evap_kin);
        ledger_.evap_int = add(ledger_.evap_int, r.evap_int); ledger_.radiation = add(ledger_.radiation, r.radiation);
        ledger_.release_E = add(ledger_.release_E, r.release_E);
        for (int el = 0; el < PX_NEL; ++el) {
            ledger_.dep_N[el] = add(ledger_.dep_N[el], r.dep_N[el]); ledger_.evap_N[el] = add(ledger_.evap_N[el], r.evap_N[el]); ledger_.rel_N[el] = add(ledger_.rel_N[el], r.rel_N[el]);
        }
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) for (int el = 0; el < PX_NEL; ++el) {
            patch_evap_N_[p][el] = add(patch_evap_N_[p][el], r.patch_evap_N[p][el]);
            patch_dep_N_[p][el] = add(patch_dep_N_[p][el], r.patch_dep_N[p][el]);
        }
        ledger_.lorentz = add(ledger_.lorentz, r.lorentz); ledger_.surface = add(ledger_.surface, r.surface);
        ledger_.dissipation = add(ledger_.dissipation, r.dissipation);
        ledger_.eject_E = add(ledger_.eject_E, r.eject_E); ledger_.eject_K = add(ledger_.eject_K, r.eject_K);
        for (int el = 0; el < PX_NEL; ++el) ledger_.eject_N[el] = add(ledger_.eject_N[el], r.eject_N[el]);
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) {
            for (int el = 0; el < PX_NEL; ++el) patch_ej_N_[p][el] = add(patch_ej_N_[p][el], r.patch_eject_N[p][el]);
            patch_ej_E_[p] = add(patch_ej_E_[p], r.patch_eject_E[p]); patch_ej_K_[p] = add(patch_ej_K_[p], r.patch_eject_K[p]);
            patch_WL_[p] = add(patch_WL_[p], r.patch_lorentz[p]);
        }
        if (r.window_steps > 0) { engine_->reset_accumulator(); acc_dt_ = 0; acc_n_ = 0; }
        watermarks_ = t.watermarks_;
        time_ += r.dt; ++revision_; ++generation_;
        eval_valid_ = false;
        sums_ = std::move(t.sums_); sums_valid_ = true;   // merges move words inside a slice: the slice sums are unchanged
        sigma_cache_ = std::move(t.sigma_); sigma_valid_ = !t.merged_;
        amax_ = t.amax_;
        t.owner_ = nullptr;
        return true;
    }
    StepRecord advance(const StepInput& in) { Trial t = prepare(in); StepRecord r = t.record(); if (!commit(std::move(t))) throw std::logic_error("commit failed"); return r; }

    // -------------------------------------------------------------------------------------------- checkpoint
    void save(const std::filesystem::path& path) const {
        if (acc_n_ != 0) throw PhaseRefusal("checkpoint", "accumulated window", path.string(), double(acc_n_), 0, "save after the owner step that consumes the accumulated window");
        const HostState s = state();
        Encoder e; e.text("ss304-phase-owner-v1.2-checkpoint"); e.text(config_sha_); e.text(law_sha_); e.f64(time_); e.u64(revision_);
        e.u64(s.E.size()); for (auto w : s.E) e.w(w); for (auto w : s.N) e.w(w); for (int f : s.flags) e.i64(f);
        e.u64(s.fcm.size()); for (auto w : s.fcm) e.w(w);
        for (double v : s.fcV) e.f64(v);
        for (auto w : s.fcK) e.w(w);
        for (auto w : {ledger_.E0, ledger_.joule, ledger_.loads, ledger_.deposit_E, ledger_.boundary, ledger_.evap_ref, ledger_.evap_kin, ledger_.evap_int, ledger_.radiation, ledger_.release_E}) e.w(w);
        for (int el = 0; el < PX_NEL; ++el) for (auto w : {ledger_.N0[el], ledger_.dep_N[el], ledger_.evap_N[el], ledger_.rel_N[el]}) e.w(w);
        for (auto w : {ledger_.K0, ledger_.lorentz, ledger_.surface, ledger_.dissipation, ledger_.eject_E, ledger_.eject_K}) e.w(w);
        for (int el = 0; el < PX_NEL; ++el) e.w(ledger_.eject_N[el]);
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) for (int el = 0; el < PX_NEL; ++el) { e.w(patch_evap_N_[p][el]); e.w(patch_dep_N_[p][el]); }
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) { for (int el = 0; el < PX_NEL; ++el) e.w(patch_ej_N_[p][el]); e.w(patch_ej_E_[p]); e.w(patch_ej_K_[p]); e.w(patch_WL_[p]); }
        e.f64(amax_);
        for (double pb : patch_pb_) e.f64(pb);
        e.u64(watermarks_.size()); for (const auto& [k, v] : watermarks_) { e.text(k); e.u64(v); }
        const std::string digest = e.digest();
        if (std::filesystem::exists(path)) throw PhaseRefusal("checkpoint", "destination", path.string(), 0, 0, "checkpoint destination exists");
        std::ofstream f(path.string() + ".tmp", std::ios::binary);
        f.write(reinterpret_cast<const char*>(e.bytes.data()), std::streamsize(e.bytes.size()));
        f.write(digest.data(), 64);
        f.close();
        if (!f) throw PhaseRefusal("checkpoint", "write", path.string(), 0, 0, "write failed");
        std::filesystem::rename(path.string() + ".tmp", path);
    }
    void load(const std::filesystem::path& path) {
        if (revision_ != 0) throw PhaseRefusal("checkpoint", "restart", path.string(), 0, 0, "restart requires a fresh owner");
        std::ifstream f(path, std::ios::binary);
        std::vector<unsigned char> b((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        if (b.size() < 64) throw PhaseRefusal("checkpoint", "length", path.string(), 0, 0, "truncated");
        const std::string digest(b.end() - 64, b.end());
        b.resize(b.size() - 64);
        if (sha256(b.data(), b.size()) != digest) throw PhaseRefusal("checkpoint", "digest", path.string(), 0, 0, "checkpoint digest mismatch");
        std::size_t at = 0;
        auto u64 = [&]() { if (at + 8 > b.size()) throw PhaseRefusal("checkpoint", "length", path.string(), 0, 0, "truncated"); std::uint64_t v = 0; for (int i = 0; i < 8; ++i) v |= std::uint64_t(b[at++]) << (8 * i); return v; };
        auto text = [&]() { const std::uint64_t n = u64(); if (n > 4096 || at + n > b.size()) throw PhaseRefusal("checkpoint", "text", path.string(), 0, 0, "bad text"); std::string s(b.begin() + long(at), b.begin() + long(at + n)); at += n; return s; };
        auto f64 = [&]() { const std::uint64_t u = u64(); double v; std::memcpy(&v, &u, 8); return v; };
        auto w = [&]() { px_i128 v; v.lo = u64(); v.hi = std::int64_t(u64()); return v; };
        if (text() != "ss304-phase-owner-v1.2-checkpoint" || text() != config_sha_ || text() != law_sha_) throw PhaseRefusal("checkpoint", "identity", path.string(), 0, 0, "configuration or law differs");
        const double t = f64(); const std::uint64_t rev = u64();
        HostState s; const std::uint64_t n = u64();
        if (n != std::uint64_t(dims_.n_cells)) throw PhaseRefusal("checkpoint", "dimensions", path.string(), double(n), dims_.n_cells, "cell count differs");
        s.E.resize(n); s.N.resize(PX_NEL * n); s.flags.resize(n);
        for (auto& x : s.E) x = w();
        for (auto& x : s.N) x = w();
        for (auto& x : s.flags) x = int(std::int64_t(u64()));
        const std::uint64_t nf = u64();
        if (nf != std::uint64_t(PX_NEL * std::max(dims_.n_fc, 1))) throw PhaseRefusal("checkpoint", "dimensions", path.string(), double(nf), 0, "face-column count differs");
        s.fcm.resize(nf); for (auto& x : s.fcm) x = w();
        const std::size_t nfc = std::max(dims_.n_fc, 1);
        s.fcV.resize(3 * nfc); for (auto& x : s.fcV) x = f64();
        s.fcK.resize(nfc); for (auto& x : s.fcK) x = w();
        Ledger L;
        for (auto* x : {&L.E0, &L.joule, &L.loads, &L.deposit_E, &L.boundary, &L.evap_ref, &L.evap_kin, &L.evap_int, &L.radiation, &L.release_E}) *x = w();
        for (int el = 0; el < PX_NEL; ++el) for (auto* x : {&L.N0[el], &L.dep_N[el], &L.evap_N[el], &L.rel_N[el]}) *x = w();
        for (auto* x : {&L.K0, &L.lorentz, &L.surface, &L.dissipation, &L.eject_E, &L.eject_K}) *x = w();
        for (int el = 0; el < PX_NEL; ++el) L.eject_N[el] = w();
        auto pe = patch_evap_N_, pd = patch_dep_N_, pjn = patch_ej_N_;
        auto pje = patch_ej_E_, pjk = patch_ej_K_, pwl = patch_WL_;
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) for (int el = 0; el < PX_NEL; ++el) { pe[p][el] = w(); pd[p][el] = w(); }
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) { for (int el = 0; el < PX_NEL; ++el) pjn[p][el] = w(); pje[p] = w(); pjk[p] = w(); pwl[p] = w(); }
        const double amax = f64();
        std::vector<double> pb(patch_specs_.size()); for (auto& x : pb) x = f64();
        std::map<std::string, std::uint64_t> wm; const std::uint64_t nw = u64();
        for (std::uint64_t i = 0; i < nw; ++i) { auto k = text(); wm[k] = u64(); }
        if (at != b.size()) throw PhaseRefusal("checkpoint", "length", path.string(), 0, 0, "trailing bytes");
        const HostElements el = build_elements(s.flags);
        const int other = work_slots(slot_).first;
        engine_->write_state(other, s);
        ensure_elements(el);
        slot_ = other; elements_ = el; accepted_flags_ = s.flags;
        ledger_ = L; patch_evap_N_ = pe; patch_dep_N_ = pd; patch_pb_ = pb; watermarks_ = wm; time_ = t; revision_ = rev; ++generation_;
        patch_ej_N_ = pjn; patch_ej_E_ = pje; patch_ej_K_ = pjk; patch_WL_ = pwl; amax_ = amax;
        state_changed();
        if (!ledger_closes()) throw PhaseRefusal("checkpoint", "ledger", path.string(), 0, 0, "restored ledgers do not close");
    }
};

}  // namespace ss304_phase_v1_2
