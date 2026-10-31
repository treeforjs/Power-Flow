// phase_owner.hpp -- SS304 phase-change thermal owner v1: the single owner of foil metal mass (Fe, Cr, Ni atom
// inventories) and condensed-phase energy for all foil slices, from 298.2 K through melting into strong
// vaporization (law domain 298.2-5000 K).
//
// Numerics (per slice: 2-D cross-section, nx columns across the width x ny rows through the thickness, fixed
// reference geometry, no axial conduction, as Hall v3):
//   * state per cell: energy E (J, 128-bit fixed point, unit 2^-96 J) and atoms of Fe, Cr, Ni (unit 2^-32 atom);
//     T = h^-1(E/m) from the phase law (enthalpy method; melting is the mushy-zone enthalpy, no front tracking);
//   * explicit finite-volume conduction with harmonic-mean face conductivity; each face energy is converted to
//     the fixed-point unit once and added to one cell and subtracted from the other, so conduction conserves
//     energy exactly; sub-steps of equal length below a state-independent stability bound;
//   * sources: Hall Joule energy per EM cell split exactly over its sub-cells (by sigma*V, parallel axial
//     conductors, or by V), surface loads, metal deposits, quantized releases or expected HKL evaporation; every
//     step total is divided over the sub-steps by exact integer division;
//   * evaporation removes atoms of each element and energy eps0_i + kinetic per atom (eps0_i constant), so the
//     vapor carries exactly the energy the condensed phase loses;
//   * surface cells that fall below phi_merge of their reference mass are merged into the inward neighbour at the
//     end of the step (exact transfer); exposure is a pure function of the void flags.
// Exact ledgers: E_now = E_0 + Joule + loads + deposits + boundary - evaporation(reference + kinetic)
// - releases, and N_i,now = N_i,0 + deposited - evaporated - released, as 128-bit integer identities checked on
// every step. Refusals are typed (PhaseRefusal) and occur only at the law domain (298.2 K, 5000 K), on depletion,
// on invalid/duplicate inputs, or when a requested step exceeds the configured sub-step budget.
#pragma once
#include "engine.hpp"
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <span>
#include <tuple>

namespace ss304_phase_v1 {

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
    double width = 0.005, thickness = 0.0002, height = 0.0002;   // m (height = axial slice length dz)
    int nx_em = 25, ny_em = 1, rx = 1;                            // EM foil cells per slice, thermal columns per EM column
    std::vector<double> y_edges;                                  // thermal row edges 0..thickness (must contain the EM row edges)
    double rho_ref = 8000.0, T0 = 300.0, start_time = 0.0;
    std::array<double, 3> mass_fraction = nominal_mass_fraction;
    LawOptions law;
    int joule_mode = 1;                                           // 1 sigma-weighted (default), 0 uniform per EM cell (Hall v3)
    int evap_mode = 1;                                            // 0 none, 1 expected flux, 2 quantized releases
    std::array<int, 4> bc_kind{0, 0, 0, 0};
    std::array<double, 4> bc_T{0, 0, 0, 0};
    std::array<int, 4> evap_face{1, 1, 1, 1};
    double phi_merge = 0.25, phi_floor = 0.05, stability_safety = 0.4;
    std::int64_t max_substeps = 1000000;
    std::vector<PatchSpec> patches;                               // per (slice, face): tiling of [0,1]; default whole faces
};

// Symmetric thermal rows: geometric growth from both faces (first cell dy0, ratio r, capped at dy_max), with the
// EM row edges inserted.
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
    std::vector<double> em_power_W;           // [S * ny_em * nx_em], slice-major, then EM row, then EM column
    std::vector<double> em_energy_J;          // alternative to em_power_W: Joule energy of the step per EM cell (J), e.g.
                                              // accumulated by the caller over several field steps (super-cycling)
    std::vector<SurfaceLoad> loads;
    std::vector<Deposit> deposits;
    std::vector<Release> releases;            // quantized mode only
    std::string quote_sha256;                 // quantized mode: the quote the releases were drawn from
    std::string field_state_sha256;
};

struct Ledger {
    pt_i128 E0{}, joule{}, loads{}, deposit_E{}, boundary{}, evap_ref{}, evap_kin{}, release_E{};
    std::array<pt_i128, 3> N0{}, dep_N{}, evap_N{}, rel_N{};
};
struct PatchReport {
    std::string id; int slice = 0, face = 0; double begin = 0, end = 1;
    double area_m2 = 0, mean_T = 0, max_T = 0, liquid_fraction = 0, melt_depth_m = 0, surface_density = 0;
    std::array<double, 3> flux_atoms_m2_s{}, expected_atoms{};       // at the accepted state; expected over quote dt
    std::array<double, 3> cumulative_evaporated_atoms{}, cumulative_deposited_atoms{};
    double cumulative_evaporated_kg = 0, recession_m = 0;
};
struct Quote { double dt = 0; std::uint64_t revision = 0; std::vector<std::array<double, 3>> expected_atoms; std::vector<double> mean_T; std::string sha256; };
struct StepRecord {
    double dt = 0; std::int64_t n_sub = 0; double tau = 0;
    pt_i128 dE{}, joule{}, loads{}, deposit_E{}, boundary{}, evap_ref{}, evap_kin{}, release_E{};
    std::array<pt_i128, 3> dN{}, dep_N{}, evap_N{}, rel_N{};
    double joule_requested_J = 0, max_T = 0; int merges = 0;
    std::vector<std::array<pt_i128, 3>> patch_evap_N;   // per patch: evaporated or released atoms this step
    std::vector<std::array<pt_i128, 3>> patch_dep_N;
    std::string inputs_sha256, field_state_sha256;
    double J(pt_i128 w) const { return pt_to_double(w, PT_E_INV); }
};

class Owner;
class Trial {
    friend class Owner;
    const Owner* owner_ = nullptr;
    std::uint64_t generation_ = 0, base_revision_ = 0;
    int slot_ = -1;
    bool merged_ = false;          // end-of-step merges changed the void flags (elements_ and flags_ are set)
    HostElements elements_;
    std::vector<int> flags_;
    std::map<std::string, std::uint64_t> watermarks_;
    StepRecord record_;
public:
    const StepRecord& record() const { return record_; }
};

class Owner {
    Config cfg_;
    pt_law law_{};
    std::string law_sha_, config_sha_;
    pt_dims dims_{};
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
    std::vector<std::array<pt_i128, 3>> patch_evap_N_, patch_dep_N_;
    double tau_max_ = 0;
    // cached evaluation of the accepted state
    mutable bool eval_valid_ = false;
    mutable double eval_dt_ = -1;
    mutable EngineResult eval_;
    mutable std::vector<double> sigma_em_;
    int per_slice() const { return dims_.nx * dims_.ny; }
    std::string slice_text(int s) const { return "foil-" + std::to_string(cfg_.slices[s].body) + "-slice-" + std::to_string(cfg_.slices[s].k); }

    // -------------------------------------------------------------------------------------------- setup
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
        double wsum = 0; for (double w : c.mass_fraction) { need(w >= 0, "mass fractions"); wsum += w; }
        need(std::abs(wsum - 1.0) < 1e-12, "mass fractions sum to 1");
        law_ = make_law(c.law);
        law_sha_ = ss304_phase_v1::law_sha256(law_);
        need(c.T0 >= law_.T_min && c.T0 <= law_.T_max, "initial temperature inside the law domain");
        if (c.evap_mode != 0) need(law_.vapor_enabled || c.law.kind == 1, "vapor law disabled with evaporation requested");
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
        // rows and EM rows
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
        geom_.vref.resize(dims_.n_cells); geom_.mref.resize(dims_.n_cells);
        for (int cidx = 0; cidx < dims_.n_cells; ++cidx) {
            const int j = (cidx / nx) % ny;
            geom_.vref[cidx] = dims_.dx * geom_.row_dy[j] * dims_.dz;
            geom_.mref[cidx] = c.rho_ref * geom_.vref[cidx];
        }
        // patches
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
        for (const auto& p : patch_specs_) {
            std::ostringstream o; o << slice_text(p.slice) << "/" << face_name(p.face) << "/" << p.begin << "-" << p.end;
            patch_ids_.push_back(o.str());
        }
        dims_.n_patches = int(patch_specs_.size());
        // configuration identity
        Encoder e;
        e.text("ss304-phase-owner-v1-config");
        e.text(law_sha_);
        e.u64(c.slices.size());
        for (const auto& s : c.slices) { e.i64(s.body); e.i64(s.k); e.f64(s.x0); e.f64(s.y0); }
        for (double v : {c.width, c.thickness, c.height, c.rho_ref, c.T0, c.start_time, c.phi_merge, c.phi_floor, c.stability_safety}) e.f64(v);
        e.i64(c.nx_em); e.i64(c.ny_em); e.i64(c.rx); e.i64(c.joule_mode); e.i64(c.evap_mode); e.i64(c.max_substeps);
        e.u64(c.y_edges.size()); for (double y : c.y_edges) e.f64(y);
        for (double w : c.mass_fraction) e.f64(w);
        for (int f = 0; f < 4; ++f) { e.i64(c.bc_kind[f]); e.f64(c.bc_T[f]); e.i64(c.evap_face[f]); }
        e.u64(patch_specs_.size()); for (const auto& p : patch_specs_) { e.i64(p.slice); e.i64(p.face); e.f64(p.begin); e.f64(p.end); }
        config_sha_ = e.digest();
        tau_max_ = stability_bound();
    }
    // State-independent explicit stability bound (s): safety * min_c C_min / G_max.
    double stability_bound() const {
        double cp_min, k_max;
        if (law_.kind == 0) {
            cp_min = pt_ss_cap(law_.T_min);
            k_max = 0;
            for (int i = 0; i <= 4000; ++i) { const double T = law_.T_min + (law_.T_max - law_.T_min) * i / 4000.0; k_max = std::max(k_max, pt_k(&law_, T)); }
        } else {
            cp_min = std::min(law_.c_cps, law_.c_cpl);
            k_max = std::max(law_.c_ks, law_.c_kl);
        }
        // evaporation stiffness bound: d/dT of sum_i j_i(T, x=1) (eps0_i + 2kT) at T_max (convex increasing)
        double dq = 0;
        if (law_.vapor_enabled && cfg_.evap_mode != 0) {
            auto q = [&](double T) { double s = 0; int st = 0; for (int el = 0; el < 3; ++el) s += pt_flux(&law_, el, T, 1.0, &st) * (law_.eps0[el] + 2 * boltzmann * T); return s; };
            dq = 2.0 * (q(law_.T_max) - q(law_.T_max - 10.0)) / 10.0;
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
                // a boundary cell can also inherit up to two more exposed faces after merges; bounded by the
                // interior conductance it replaces (k_max * g >= 0), so no extra term is needed.
            }
            // cells can only lose mass through evaporation or releases (down to phi_floor within a step)
            const double C = (cfg_.evap_mode != 0 ? cfg_.phi_floor : 1.0) * geom_.mref[c] * cp_min;
            if (G > 0) best = std::min(best, C / G);
        }
        return cfg_.stability_safety * best;
    }
    // Exposure is a pure function of the void flags: face f of a non-void cell is a surface element iff every
    // cell beyond it in direction f up to the slice boundary is void.
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
                for (;;) {
                    const int nb = pt_neighbor(&dims_, cc, f);
                    if (nb < 0) break;
                    if (!(flags[nb] & 1)) { reaches = false; break; }
                    cc = nb;
                }
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
        for (std::size_t p = 0; p <= patch_specs_.size(); ++p) {
            while (q < items.size() && std::size_t(items[q].patch) < p) ++q;
            el.patch_begin[p] = int(q);
        }
        return el;
    }
    HostState initial_state() const {
        HostState s;
        const int n = dims_.n_cells;
        s.E.assign(n, pt_zero()); s.N.assign(3 * n, pt_zero()); s.flags.assign(n, 0);
        const double h0 = pt_h(&law_, cfg_.T0);
        for (int c = 0; c < n; ++c) {
            int bad = 0;
            double atoms[3], m = 0;
            for (int el = 0; el < 3; ++el) {
                s.N[3 * c + el] = pt_from_double(cfg_.mass_fraction[el] * geom_.mref[c] / law_atom_mass(el), PT_N_POW2, &bad);
                atoms[el] = pt_to_double(s.N[3 * c + el], PT_N_INV);
                m += atoms[el] * law_atom_mass(el);
            }
            s.E[c] = pt_from_double(m * h0, PT_E_POW2, &bad);
            need(!bad, "initial state outside the fixed-point range");
        }
        return s;
    }
    double law_atom_mass(int el) const { return law_.atom_mass[el]; }
    void ensure_elements(const HostElements& el) { engine_->set_elements(el); dims_.n_el = int(el.el_cell.size()); }
    const EngineResult& evaluated(double dt) const {
        if (!eval_valid_ || eval_dt_ != dt) {
            pt_dims d = dims_; d.dt = dt; d.n_sub = 1; d.tau = dt;
            const_cast<Owner*>(this)->engine_->evaluate(slot_, d, eval_, &sigma_em_);
            eval_valid_ = true; eval_dt_ = dt;
        }
        return eval_;
    }
    static pt_i128 add(pt_i128 a, pt_i128 b) { int ov = 0; auto r = pt_add(a, b, &ov); if (ov) throw std::overflow_error("ledger word overflow"); return r; }
    static pt_i128 sub(pt_i128 a, pt_i128 b) { int ov = 0; auto r = pt_sub(a, b, &ov); if (ov) throw std::overflow_error("ledger word overflow"); return r; }
    static pt_i128 total_E(const EngineResult& r, int S) { pt_i128 t = pt_zero(); for (int s = 0; s < S; ++s) t = add(t, r.slice_sum[PT_SLICE_SUM * s]); return t; }
    static pt_i128 total_N(const EngineResult& r, int S, int el) { pt_i128 t = pt_zero(); for (int s = 0; s < S; ++s) t = add(t, r.slice_sum[PT_SLICE_SUM * s + 1 + el]); return t; }
    // (slice, x, y) in the field frame -> surface element (walks inward through void cells).
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
            j = int(std::upper_bound(cfg_.y_edges.begin(), cfg_.y_edges.end(), yy) - cfg_.y_edges.begin()) - 1;
            j = std::clamp(j, 0, dims_.ny - 1);
        }
        int c = (s * dims_.ny + j) * dims_.nx + i;
        const int inward = f == 0 ? 1 : f == 1 ? 0 : f == 2 ? 3 : 2;
        const auto& flags = accepted_flags_;
        while (c >= 0 && (flags[c] & 1)) c = pt_neighbor(&dims_, c, inward);
        if (c < 0) throw PhaseRefusal("column-depleted", what, slice_text(s), 0, 0, "no material along the inward normal");
        const int e = elements_.el_of[4 * c + f];
        if (e < 0) throw std::logic_error("exposed cell without element");
        face = f;
        return e;
    }
    std::vector<int> accepted_flags_;

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
        for (int s = 0; s < dims_.S; ++s) if (r.slice_status[4 * s + 1] != PT_OK) throw PhaseRefusal("invalid-initial-state", "initial temperature", slice_text(s), cfg_.T0, r.slice_status[4 * s + 1], "outside the law domain or not representable (status " + std::to_string(r.slice_status[4 * s + 1]) + ")");
        ledger_.E0 = total_E(r, dims_.S);
        for (int el = 0; el < 3; ++el) ledger_.N0[el] = total_N(r, dims_.S, el);
        patch_evap_N_.assign(patch_specs_.size(), {pt_zero(), pt_zero(), pt_zero()});
        patch_dep_N_ = patch_evap_N_;
    }
    Owner(const Owner&) = delete; Owner& operator=(const Owner&) = delete;

    const Config& config() const { return cfg_; }
    const pt_law& law() const { return law_; }
    const pt_dims& dims() const { return dims_; }
    const HostGeometry& geometry() const { return geom_; }
    const std::string& law_sha256() const { return law_sha_; }
    const std::string& config_sha256() const { return config_sha_; }
    const std::string engine_name() const { return engine_->name(); }
    double time_s() const { return time_; }
    std::uint64_t revision() const { return revision_; }
    double stability_tau() const { return tau_max_; }
    std::int64_t substeps_for(double dt) const { return std::max<std::int64_t>(1, std::int64_t(std::ceil(dt / tau_max_))); }
    const Ledger& ledger() const { return ledger_; }
    const std::vector<std::string>& patch_ids() const { return patch_ids_; }
    const std::vector<PatchSpec>& patch_specs() const { return patch_specs_; }
    int n_patches() const { return int(patch_specs_.size()); }
    int em_index(int slice, int je, int ie) const { return (slice * cfg_.ny_em + je) * cfg_.nx_em + ie; }

    pt_i128 energy_word() const { return total_E(evaluated(eval_valid_ ? eval_dt_ : 0.0), dims_.S); }
    std::array<pt_i128, 3> atom_words() const { const auto& r = evaluated(eval_valid_ ? eval_dt_ : 0.0); return {total_N(r, dims_.S, 0), total_N(r, dims_.S, 1), total_N(r, dims_.S, 2)}; }
    double energy_J() const { return pt_to_double(energy_word(), PT_E_INV); }
    double max_temperature_K() const { const auto& r = evaluated(eval_valid_ ? eval_dt_ : 0.0); double m = 0; for (double t : r.slice_maxT) m = std::max(m, t); return m; }
    // Exact ledger identities of the accepted state.
    bool ledger_closes() const {
        pt_i128 e = ledger_.E0;
        for (auto w : {ledger_.joule, ledger_.loads, ledger_.deposit_E, ledger_.boundary}) e = add(e, w);
        for (auto w : {ledger_.evap_ref, ledger_.evap_kin, ledger_.release_E}) e = sub(e, w);
        if (!pt_eq(e, energy_word())) return false;
        const auto N = atom_words();
        for (int el = 0; el < 3; ++el) {
            pt_i128 n = add(ledger_.N0[el], ledger_.dep_N[el]);
            n = sub(sub(n, ledger_.evap_N[el]), ledger_.rel_N[el]);
            if (!pt_eq(n, N[el])) return false;
        }
        return true;
    }
    // Effective axial electrical conductivity per EM cell (S/m) at the accepted state; Hall's sigma(T) input.
    std::vector<double> sigma_em() const { evaluated(eval_valid_ ? eval_dt_ : 0.0); return sigma_em_; }
    std::vector<double> temperatures() const { return const_cast<Owner*>(this)->engine_->temperatures(slot_, dims_); }
    HostState state() const { return const_cast<Owner*>(this)->engine_->read_state(slot_); }
    std::vector<PatchReport> patch_reports(double dt = 0.0) const {
        const auto& r = evaluated(dt);
        std::vector<PatchReport> out;
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) {
            const double* o = r.patch_out.data() + PT_PATCH_OUT * p;
            PatchReport q;
            q.id = patch_ids_[p]; q.slice = patch_specs_[p].slice; q.face = patch_specs_[p].face; q.begin = patch_specs_[p].begin; q.end = patch_specs_[p].end;
            q.area_m2 = o[0]; q.mean_T = o[1]; q.max_T = o[2]; q.liquid_fraction = o[3]; q.melt_depth_m = o[4]; q.surface_density = o[11];
            double kg = 0;
            for (int el = 0; el < 3; ++el) {
                q.expected_atoms[el] = o[5 + el]; q.flux_atoms_m2_s[el] = o[8 + el];
                q.cumulative_evaporated_atoms[el] = pt_to_double(patch_evap_N_[p][el], PT_N_INV);
                q.cumulative_deposited_atoms[el] = pt_to_double(patch_dep_N_[p][el], PT_N_INV);
                kg += (q.cumulative_evaporated_atoms[el] - q.cumulative_deposited_atoms[el]) * law_atom_mass(el);
            }
            q.cumulative_evaporated_kg = kg;
            const double rho = law_.kind == 0 && q.surface_density > 0 ? q.surface_density : cfg_.rho_ref;
            q.recession_m = q.area_m2 > 0 ? kg / (rho * q.area_m2) : 0;
            out.push_back(q);
        }
        return out;
    }
    // Expected evaporated atoms per patch and element over dt at the accepted state (quantized-release quote).
    Quote quote(double dt) const {
        if (!(std::isfinite(dt) && dt > 0)) throw PhaseRefusal("invalid-step", "dt", "quote", dt, 0, "dt must be finite and positive");
        const auto& r = evaluated(dt);
        Quote q; q.dt = dt; q.revision = revision_;
        Encoder e; e.text("ss304-phase-owner-v1-quote"); e.text(config_sha_); e.u64(revision_); e.f64(time_); e.f64(dt);
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) {
            const double* o = r.patch_out.data() + PT_PATCH_OUT * p;
            q.expected_atoms.push_back({o[5], o[6], o[7]}); q.mean_T.push_back(o[1]);
            for (int el = 0; el < 3; ++el) e.f64(o[5 + el]);
        }
        q.sha256 = e.digest();
        return q;
    }
    std::string state_sha256() const {
        const HostState s = state();
        Encoder e; e.text("ss304-phase-owner-v1-state"); e.text(config_sha_); e.u64(revision_); e.f64(time_);
        e.u64(s.E.size()); for (auto w : s.E) e.w(w);
        for (auto w : s.N) e.w(w);
        for (int f : s.flags) e.i64(f);
        for (auto w : {ledger_.E0, ledger_.joule, ledger_.loads, ledger_.deposit_E, ledger_.boundary, ledger_.evap_ref, ledger_.evap_kin, ledger_.release_E}) e.w(w);
        for (int el = 0; el < 3; ++el) for (auto w : {ledger_.N0[el], ledger_.dep_N[el], ledger_.evap_N[el], ledger_.rel_N[el]}) e.w(w);
        e.u64(watermarks_.size()); for (const auto& [k, v] : watermarks_) { e.text(k); e.u64(v); }
        return e.digest();
    }

    // Initial temperature profile T(slice, x, y) at cell centres (field frame), before the first step.
    template <class F> void initialize_temperature(F&& profile) {
        if (revision_ != 0) throw PhaseRefusal("invalid-initialization", "revision", "owner", double(revision_), 0, "profile must be set before the first step");
        HostState s = state();
        for (int c = 0; c < dims_.n_cells; ++c) {
            const int sl = c / per_slice(), j = (c / dims_.nx) % dims_.ny, i = c % dims_.nx;
            const double x = cfg_.slices[sl].x0 + (i + 0.5) * dims_.dx, y = cfg_.slices[sl].y0 + 0.5 * (cfg_.y_edges[j] + cfg_.y_edges[j + 1]);
            const double T = profile(sl, x, y);
            if (!(T >= law_.T_min && T <= law_.T_max)) throw PhaseRefusal("invalid-initialization", "temperature", slice_text(sl), T, law_.T_max, "initial temperature outside the law domain");
            double m = 0; for (int el = 0; el < 3; ++el) m += pt_to_double(s.N[3 * c + el], PT_N_INV) * law_atom_mass(el);
            int bad = 0; s.E[c] = pt_from_double(m * pt_h(&law_, T), PT_E_POW2, &bad);
            if (bad) throw PhaseRefusal("fixed-point-range", "initial energy", slice_text(sl), T, 0, "outside the fixed-point range");
        }
        engine_->write_state(slot_, s);
        eval_valid_ = false;
        ledger_.E0 = energy_word();
        ++generation_;
    }

    // -------------------------------------------------------------------------------------------- step
    Trial prepare(const StepInput& in) {
        const double dt = in.dt;
        if (!(std::isfinite(dt) && dt > 0 && std::isfinite(time_ + dt) && time_ + dt > time_))
            throw PhaseRefusal("invalid-step", "dt", "owner", dt, 0, "dt must be finite, positive and advance the clock");
        const std::int64_t n_sub = substeps_for(dt);
        if (n_sub > cfg_.max_substeps || n_sub > 0x7fffffff)
            throw PhaseRefusal("substep-budget", "sub-steps", "owner", double(n_sub), double(cfg_.max_substeps), "the step needs more explicit sub-steps than the configured budget");
        const bool by_energy = !in.em_energy_J.empty();
        const auto& jv = by_energy ? in.em_energy_J : in.em_power_W;
        if (int(jv.size()) != dims_.n_em || (by_energy && !in.em_power_W.empty()))
            throw PhaseRefusal("invalid-input", "Joule vector", "owner", double(jv.size()), dims_.n_em, "one Joule power or one Joule energy (not both) per EM foil cell");
        StepData data;
        data.em_energy.resize(dims_.n_em);
        double joule_req = 0;
        pt_i128 joule_exact = pt_zero();
        for (int e = 0; e < dims_.n_em; ++e) {
            const double P = jv[e];
            if (!(std::isfinite(P) && P >= 0)) throw PhaseRefusal("invalid-input", by_energy ? "Joule energy J" : "Joule power W", "em-cell-" + std::to_string(e), P, 0, "Joule input must be finite and nonnegative");
            data.em_energy[e] = by_energy ? P : P * dt;
            joule_req += data.em_energy[e];
            int bad = 0;
            joule_exact = add(joule_exact, pt_from_double(data.em_energy[e], PT_E_POW2, &bad));
            if (bad) throw PhaseRefusal("fixed-point-range", "Joule energy", "em-cell-" + std::to_string(e), data.em_energy[e], 0, "outside the 128-bit energy range");
        }
        // exactly-once identities
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
        std::map<int, std::array<pt_i128, 5>> per_el;
        pt_i128 loads_total = pt_zero(), depE_total = pt_zero();
        std::array<pt_i128, 3> depN_total{pt_zero(), pt_zero(), pt_zero()};
        std::vector<std::array<pt_i128, 3>> patch_dep(patch_specs_.size(), {pt_zero(), pt_zero(), pt_zero()});
        for (const auto& l : in.loads) {
            admit(l.producer, l.sequence, "surface load");
            if (l.channel.empty()) throw PhaseRefusal("load-identity", "surface load channel", l.producer, 0, 0, "channel required");
            if (!std::isfinite(l.energy_J)) throw PhaseRefusal("nonfinite-input", "surface load energy", l.producer, l.energy_J, 0, "load energy not finite");
            int f; const int e = locate(l.slice, l.x, l.y, f, "surface load");
            int bad = 0; const pt_i128 q = pt_from_double(l.energy_J, PT_E_POW2, &bad);
            if (bad) throw PhaseRefusal("fixed-point-range", "surface load energy", l.producer, l.energy_J, 0, "outside the 128-bit energy range");
            auto& a = per_el.try_emplace(e, std::array<pt_i128, 5>{pt_zero(), pt_zero(), pt_zero(), pt_zero(), pt_zero()}).first->second;
            a[0] = add(a[0], q); loads_total = add(loads_total, q);
        }
        for (const auto& dpt : in.deposits) {
            admit(dpt.producer, dpt.sequence, "metal deposit");
            if (!(dpt.element >= 0 && dpt.element < 3 && std::isfinite(dpt.atoms) && dpt.atoms >= 0 && std::isfinite(dpt.kinetic_J) && dpt.kinetic_J >= 0))
                throw PhaseRefusal("invalid-input", "metal deposit", dpt.producer, dpt.atoms, 0, "element, atoms and kinetic energy must be valid");
            int f; const int e = locate(dpt.slice, dpt.x, dpt.y, f, "metal deposit");
            int bad = 0;
            const pt_i128 qn = pt_from_double(dpt.atoms, PT_N_POW2, &bad);
            const double a_q = pt_to_double(qn, PT_N_INV);
            const pt_i128 qe = add(pt_from_double(a_q * law_.eps0[dpt.element], PT_E_POW2, &bad), pt_from_double(dpt.kinetic_J, PT_E_POW2, &bad));
            if (bad) throw PhaseRefusal("fixed-point-range", "metal deposit", dpt.producer, dpt.atoms, 0, "outside the fixed-point range");
            auto& a = per_el.try_emplace(e, std::array<pt_i128, 5>{pt_zero(), pt_zero(), pt_zero(), pt_zero(), pt_zero()}).first->second;
            a[1 + dpt.element] = add(a[1 + dpt.element], qn); a[4] = add(a[4], qe);
            depN_total[dpt.element] = add(depN_total[dpt.element], qn); depE_total = add(depE_total, qe);
            const int p = elements_.el_patch[e];
            patch_dep[p][dpt.element] = add(patch_dep[p][dpt.element], qn);
        }
        for (const auto& [e, v] : per_el) { data.sc_el.push_back(e); for (auto w : v) data.sc_val.push_back(w); }
        // quantized releases
        data.rel_total_N.assign(3 * patch_specs_.size(), pt_zero());
        data.rel_kin_J.assign(3 * patch_specs_.size(), 0.0);
        std::array<pt_i128, 3> relN_total{pt_zero(), pt_zero(), pt_zero()};
        if (cfg_.evap_mode == 2) {
            const Quote q = quote(dt);
            if (!in.releases.empty() && in.quote_sha256 != q.sha256) throw PhaseRefusal("release-binding", "quote digest", "owner", 0, 0, "releases were not drawn from this owner's quote for this state and dt");
            std::set<std::pair<int, int>> once;
            for (const auto& r : in.releases) {
                if (!(r.patch >= 0 && r.patch < int(patch_specs_.size()) && r.element >= 0 && r.element < 3 && once.emplace(r.patch, r.element).second))
                    throw PhaseRefusal("release-binding", "release identity", "patch-" + std::to_string(r.patch), 0, 0, "one release per patch and element");
                const double w = r.weight_atoms, n = double(r.macro_count), D = q.expected_atoms[r.patch][r.element];
                if (!(std::isfinite(w) && w >= 1 && w == std::floor(w) && n * w < 9007199254740992.0 && std::isfinite(r.kinetic_J) && r.kinetic_J >= 0))
                    throw PhaseRefusal("release-binding", "release weight", patch_ids_[r.patch], w, 0, "integer weight, exact atom count and nonnegative kinetic energy required");
                const double lo = std::floor(D / w), hi = std::ceil(D / w);
                if (!(n == lo || n == hi)) throw PhaseRefusal("release-binding", "macro count", patch_ids_[r.patch], n, D / w, "n must be floor or ceil of D/w0 (stochastic rounding of the quote)");
                if (r.macro_count == 0 && r.kinetic_J != 0) throw PhaseRefusal("release-binding", "kinetic energy", patch_ids_[r.patch], r.kinetic_J, 0, "no atoms released");
                int bad = 0;
                const pt_i128 R = pt_from_double(n * w, PT_N_POW2, &bad);
                data.rel_total_N[3 * r.patch + r.element] = R;
                data.rel_kin_J[3 * r.patch + r.element] = r.kinetic_J;
                relN_total[r.element] = add(relN_total[r.element], R);
            }
        } else if (!in.releases.empty()) {
            throw PhaseRefusal("release-binding", "releases", "owner", 0, 0, "releases are accepted only in quantized mode");
        }
        // run
        pt_dims d = dims_;
        d.n_el = int(elements_.el_cell.size());
        d.dt = dt; d.n_sub = int(n_sub); d.tau = dt / double(n_sub); d.sub = 0;
        const EngineResult& before = evaluated(eval_valid_ ? eval_dt_ : 0.0);
        const pt_i128 E_before = total_E(before, dims_.S);
        const std::array<pt_i128, 3> N_before{total_N(before, dims_.S, 0), total_N(before, dims_.S, 1), total_N(before, dims_.S, 2)};
        EngineResult r;
        engine_->step(slot_, d, data, r);
        eval_valid_ = false;   // engine scratch now describes the trial
        // typed refusal from the first failing cell (lowest slice index first)
        for (int s = 0; s < dims_.S; ++s) {
            const int cell = r.slice_status[4 * s], code = r.slice_status[4 * s + 1];
            if (cell < 0 || code == PT_OK) continue;
            const int local = cell % per_slice();
            const std::string id = slice_text(s) + "/row-" + std::to_string(local / dims_.nx) + "/col-" + std::to_string(local % dims_.nx);
            switch (code) {
            case PT_ERR_T_HIGH: throw PhaseRefusal("temperature-above-domain", "specific energy", id, r.slice_maxT[s], law_.T_max, "the step would heat a cell above the law domain (T_max); no state change");
            case PT_ERR_T_LOW: throw PhaseRefusal("temperature-below-domain", "specific energy", id, 0, law_.T_min, "the step would cool a cell below the law domain (T_min); no state change");
            case PT_ERR_DEPLETED: throw PhaseRefusal("surface-depletion-within-step", "cell mass fraction", id, 0, cfg_.phi_floor, "a cell lost more than its admitted mass within one step; use a shorter step");
            case PT_ERR_NEGATIVE_INV: throw PhaseRefusal("release-exceeds-inventory", "element inventory", id, 0, 0, "a release or evaporation exceeds the element inventory of the surface cell");
            case PT_ERR_FIXED_RANGE: throw PhaseRefusal("fixed-point-range", "energy/atom word", id, 0, 0, "a value left the 128-bit fixed-point range");
            default: throw PhaseRefusal("nonfinite-state", "cell state", id, 0, 0, "nonfinite intermediate value");
            }
        }
        // exact step ledger
        StepRecord rec;
        rec.dt = dt; rec.n_sub = n_sub; rec.tau = d.tau; rec.joule_requested_J = joule_req; rec.field_state_sha256 = in.field_state_sha256;
        rec.joule = joule_exact; rec.loads = loads_total; rec.deposit_E = depE_total; rec.dep_N = depN_total;
        rec.patch_evap_N.assign(patch_specs_.size(), {pt_zero(), pt_zero(), pt_zero()});
        rec.patch_dep_N = patch_dep;
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) {
            const pt_i128* ps = r.patch_sum.data() + PT_PATCH_SUM * p;
            for (int el = 0; el < 3; ++el) { rec.evap_N[el] = add(rec.evap_N[el], ps[el]); rec.patch_evap_N[p][el] = add(rec.patch_evap_N[p][el], ps[el]); }
            rec.evap_ref = add(rec.evap_ref, ps[3]); rec.evap_kin = add(rec.evap_kin, ps[4]); rec.boundary = add(rec.boundary, ps[5]);
            rec.release_E = sub(rec.release_E, ps[8]);   // in_rel_E is negative (removal)
            for (int el = 0; el < 3; ++el) rec.patch_evap_N[p][el] = add(rec.patch_evap_N[p][el], data.rel_total_N[3 * p + el]);
        }
        rec.rel_N = relN_total;
        const pt_i128 E_after = total_E(r, dims_.S);
        rec.dE = sub(E_after, E_before);
        pt_i128 expect = add(add(add(rec.joule, rec.loads), rec.deposit_E), rec.boundary);
        expect = sub(sub(sub(expect, rec.evap_ref), rec.evap_kin), rec.release_E);
        if (!pt_eq(expect, rec.dE)) throw std::logic_error("exact energy ledger failed to close (implementation defect)");
        for (int el = 0; el < 3; ++el) {
            rec.dN[el] = sub(total_N(r, dims_.S, el), N_before[el]);
            const pt_i128 en = sub(sub(rec.dep_N[el], rec.evap_N[el]), rec.rel_N[el]);
            if (!pt_eq(en, rec.dN[el])) throw std::logic_error("exact atom ledger failed to close (implementation defect)");
        }
        for (double t : r.slice_maxT) rec.max_T = std::max(rec.max_T, t);
        // end-of-step merges of depleted surface cells (exact transfers)
        Trial trial;
        int candidates = 0; for (int s = 0; s < dims_.S; ++s) candidates += r.slice_status[4 * s + 2];
        if (candidates > 0) {
            trial.elements_ = elements_;
            HostState st = engine_->read_state(r.trial_slot);
            for (int c = 0; c < dims_.n_cells; ++c) {
                if (st.flags[c] & 1) continue;
                double m = 0; for (int el = 0; el < 3; ++el) m += pt_to_double(st.N[3 * c + el], PT_N_INV) * law_atom_mass(el);
                bool exposed = false; for (int f = 0; f < 4; ++f) if (trial.elements_.el_of[4 * c + f] >= 0) exposed = true;
                if (!exposed || !(m < cfg_.phi_merge * geom_.mref[c])) continue;
                int target = -1;
                for (int f = 0; f < 4 && target < 0; ++f) {
                    if (trial.elements_.el_of[4 * c + f] < 0) continue;
                    const int nb = pt_neighbor(&dims_, c, f == 0 ? 1 : f == 1 ? 0 : f == 2 ? 3 : 2);
                    if (nb >= 0 && !(st.flags[nb] & 1)) target = nb;
                }
                if (target < 0) throw PhaseRefusal("column-depleted", "surface cell mass", slice_text(c / per_slice()), m, cfg_.phi_merge * geom_.mref[c], "a depleted surface cell has no inward neighbour (burn-through)");
                int ov = 0;
                st.E[target] = pt_add(st.E[target], st.E[c], &ov); st.E[c] = pt_zero();
                for (int el = 0; el < 3; ++el) { st.N[3 * target + el] = pt_add(st.N[3 * target + el], st.N[3 * c + el], &ov); st.N[3 * c + el] = pt_zero(); }
                if (ov) throw std::overflow_error("merge overflow");
                st.flags[c] |= 1;
                trial.elements_ = build_elements(st.flags);
                ++rec.merges;
            }
            engine_->write_state(r.trial_slot, st);
            if (rec.merges > 0) { trial.merged_ = true; trial.flags_ = st.flags; }
        }
        Encoder ie; ie.text("ss304-phase-owner-v1-step-inputs"); ie.f64(dt); ie.text(in.field_state_sha256);
        ie.i64(by_energy ? 1 : 0);
        for (double P : jv) ie.f64(P);
        for (std::size_t q = 0; q < data.sc_el.size(); ++q) { ie.i64(data.sc_el[q]); for (int k = 0; k < 5; ++k) ie.w(data.sc_val[5 * q + k]); }
        for (std::size_t q = 0; q < data.rel_total_N.size(); ++q) { ie.w(data.rel_total_N[q]); ie.f64(data.rel_kin_J[q]); }
        rec.inputs_sha256 = ie.digest();
        trial.owner_ = this; trial.generation_ = generation_; trial.base_revision_ = revision_;
        trial.slot_ = r.trial_slot; trial.watermarks_ = std::move(watermarks); trial.record_ = std::move(rec);
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
        ledger_.release_E = add(ledger_.release_E, r.release_E);
        for (int el = 0; el < 3; ++el) {
            ledger_.dep_N[el] = add(ledger_.dep_N[el], r.dep_N[el]); ledger_.evap_N[el] = add(ledger_.evap_N[el], r.evap_N[el]); ledger_.rel_N[el] = add(ledger_.rel_N[el], r.rel_N[el]);
        }
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) for (int el = 0; el < 3; ++el) {
            patch_evap_N_[p][el] = add(patch_evap_N_[p][el], r.patch_evap_N[p][el]);
            patch_dep_N_[p][el] = add(patch_dep_N_[p][el], r.patch_dep_N[p][el]);
        }
        watermarks_ = t.watermarks_;
        time_ += r.dt; ++revision_; ++generation_;
        eval_valid_ = false;
        t.owner_ = nullptr;
        return true;
    }
    // Convenience: prepare + commit.
    StepRecord advance(const StepInput& in) { Trial t = prepare(in); StepRecord r = t.record(); if (!commit(std::move(t))) throw std::logic_error("commit failed"); return r; }

    // -------------------------------------------------------------------------------------------- checkpoint
    void save(const std::filesystem::path& path) const {
        const HostState s = state();
        Encoder e; e.text("ss304-phase-owner-v1-checkpoint"); e.text(config_sha_); e.text(law_sha_); e.f64(time_); e.u64(revision_);
        e.u64(s.E.size()); for (auto w : s.E) e.w(w); for (auto w : s.N) e.w(w); for (int f : s.flags) e.i64(f);
        for (auto w : {ledger_.E0, ledger_.joule, ledger_.loads, ledger_.deposit_E, ledger_.boundary, ledger_.evap_ref, ledger_.evap_kin, ledger_.release_E}) e.w(w);
        for (int el = 0; el < 3; ++el) for (auto w : {ledger_.N0[el], ledger_.dep_N[el], ledger_.evap_N[el], ledger_.rel_N[el]}) e.w(w);
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) for (int el = 0; el < 3; ++el) { e.w(patch_evap_N_[p][el]); e.w(patch_dep_N_[p][el]); }
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
        auto w = [&]() { pt_i128 v; v.lo = u64(); v.hi = std::int64_t(u64()); return v; };
        if (text() != "ss304-phase-owner-v1-checkpoint" || text() != config_sha_ || text() != law_sha_) throw PhaseRefusal("checkpoint", "identity", path.string(), 0, 0, "configuration or law differs");
        const double t = f64(); const std::uint64_t rev = u64();
        HostState s; const std::uint64_t n = u64();
        if (n != std::uint64_t(dims_.n_cells)) throw PhaseRefusal("checkpoint", "dimensions", path.string(), double(n), dims_.n_cells, "cell count differs");
        s.E.resize(n); s.N.resize(3 * n); s.flags.resize(n);
        for (auto& x : s.E) x = w();
        for (auto& x : s.N) x = w();
        for (auto& x : s.flags) x = int(std::int64_t(u64()));
        Ledger L;
        for (auto* x : {&L.E0, &L.joule, &L.loads, &L.deposit_E, &L.boundary, &L.evap_ref, &L.evap_kin, &L.release_E}) *x = w();
        for (int el = 0; el < 3; ++el) for (auto* x : {&L.N0[el], &L.dep_N[el], &L.evap_N[el], &L.rel_N[el]}) *x = w();
        auto pe = patch_evap_N_, pd = patch_dep_N_;
        for (std::size_t p = 0; p < patch_specs_.size(); ++p) for (int el = 0; el < 3; ++el) { pe[p][el] = w(); pd[p][el] = w(); }
        std::map<std::string, std::uint64_t> wm; const std::uint64_t nw = u64();
        for (std::uint64_t i = 0; i < nw; ++i) { auto k = text(); wm[k] = u64(); }
        if (at != b.size()) throw PhaseRefusal("checkpoint", "length", path.string(), 0, 0, "trailing bytes");
        const HostElements el = build_elements(s.flags);
        const int other = work_slots(slot_).first;
        engine_->write_state(other, s);
        ensure_elements(el);
        slot_ = other; elements_ = el; accepted_flags_ = s.flags;
        ledger_ = L; patch_evap_N_ = pe; patch_dep_N_ = pd; watermarks_ = wm; time_ = t; revision_ = rev; ++generation_;
        eval_valid_ = false;
        if (!ledger_closes()) throw PhaseRefusal("checkpoint", "ledger", path.string(), 0, 0, "restored ledgers do not close");
    }
};

}  // namespace ss304_phase_v1
