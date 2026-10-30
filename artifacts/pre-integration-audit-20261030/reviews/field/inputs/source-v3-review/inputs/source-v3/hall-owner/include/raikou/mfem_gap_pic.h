#pragma once

#include <mfem.hpp>
#include <array>
#include <cstdint>
#include <vector>
#include <functional>

namespace raikou::mfem_backend {

inline constexpr double pic_e = 1.602176634e-19;
inline constexpr double pic_me = 9.1093837139e-31;
inline constexpr double pic_c = 299792458.;
using PicVector = std::array<double, 3>;

void pic_parallel_workers();
void pic_parallel_finish();
void pic_parallel_record(const std::string& directory);

struct GapParticle {
    std::uint64_t id = 0; // assigned by the owning GapPic on insertion
    PicVector position_m{};
    // Proper velocity u = gamma*v; z is retained even though fields are 2D.
    PicVector proper_velocity_m_s{};
    double charge_C = -pic_e; // single physical particle, not macrocharge
    double mass_kg = pic_me;
    double weight = 1; // represented number of physical particles
};

enum class GapElectrode { Lower, Upper };
enum class GapSurface { Gap, Back, Left, Right };
struct GapFoil { double left_m=0, right_m=0, bottom_m=0, top_m=0; };

struct GapPicOptions {
    int nx = 16, ny = 32;
    double width_m = 1e-3, gap_m = 1e-3, depth_m = 1e-3;
    double lower_potential_V = 0, upper_potential_V = 0;
    // Prescribed uniform B only in this first implementation. No induction solve.
    PicVector magnetic_field_T{};
    // In this mode width_m/gap_m describe the entire vacuum box. Rectangles
    // must align to the Q1 grid. All conductor nodes impose their potential;
    // particles stop at their actual first surface intersection. The outer
    // box has homogeneous Neumann field data and absorbing particle boundaries.
    bool distributed_fields = false; // root particles, disjoint MPI field rows
    bool finite_foils = false;
    GapFoil lower_foil{}, upper_foil{};
    // Quasi-3D extrusion: Poisson remains height-averaged in x/y.
    bool finite_axial_extent = false;
    double short_z_m = -.0125, feed_z_m = .0125;

};

struct GapPicDiagnostics {
    double time_s = 0;
    double particle_charge_C = 0, deposited_charge_C = 0;
    double lower_field_charge_C = 0, upper_field_charge_C = 0;
    double gauss_defect_C = 0, particle_charge_defect_C = 0;
    double field_energy_J = 0, particle_kinetic_energy_J = 0;
    double lower_absorbed_charge_C = 0, upper_absorbed_charge_C = 0;
    double lower_absorbed_energy_J = 0, upper_absorbed_energy_J = 0;
    double lower_injected_charge_C = 0, upper_injected_charge_C = 0;
    double injected_energy_J = 0;
    double authored_charge_C = 0, authored_kinetic_energy_J = 0;
    double maximum_nodal_continuity_defect_C = 0;
    double escaped_charge_C = 0, escaped_kinetic_energy_J = 0;
    double short_absorbed_charge_C = 0, short_absorbed_energy_J = 0;
    double axial_escaped_charge_C = 0; // subset of escaped_charge_C
    double poisson_residual_C_per_m = 0;
};

// Relativistic Boris momentum update. x is not advanced by this function.
[[nodiscard]] PicVector pic_boris(const PicVector& u, double charge_mass_ratio,
                                 const PicVector& E_V_m, const PicVector& B_T,
                                 double dt_s);
[[nodiscard]] double pic_kinetic_energy_J(const GapParticle& particle);
[[nodiscard]] double pic_richardson_supply_A_m2(double temperature_K,
    double work_function_eV, double richardson_A_m2_K2);

// Serial, Cartesian rectangular Q1 2D3V research prototype, SI units.
// Both y electrodes impose Dirichlet potentials and absorb particles. x walls
// reflect particles and impose homogeneous Neumann electrostatic conditions.
// Q1 nodal charge and path-integrated edge transfers satisfy discrete continuity
// for each drift, including electrode absorption. No Maxwell induction solve.
// The side walls and uniform B remain benchmark restrictions, not an MITL model.
struct PicPoissonState;
class GapPic {
public:
    explicit GapPic(GapPicOptions options);
    ~GapPic();
    void release_field_solver();
    GapPic(const GapPic&) = delete;
    GapPic& operator=(const GapPic&) = delete;
    void add_particle(const GapParticle& particle);
    // Commit a complete externally authored population after validating all members.
    // No field solve/gather can observe a partial batch. This is NOT a collision API:
    // all kinetic energy is charged to the authored reservoir, with no reactant debit.
    [[nodiscard]] std::vector<std::uint64_t> add_authored_batch(
        const std::vector<GapParticle>& particles);
    // Diagnostic callback at the actual force midpoint. Must not throw or mutate
    // this solver. Receives the particle before the kick and its resulting u.
    using PushObserver = std::function<void(const GapParticle&, const PicVector&,
        const PicVector&, const PicVector&, double)>;
    void set_push_observer(PushObserver observer) { push_observer_=std::move(observer); }

    // Caller owns the emission law. Flux-weighted half-Maxwellian launch,
    // nonrelativistic thermal distribution (T <= 1e6 K), exact represented charge.
    // Injection is an instantaneous boundary source; caller owns its time rule.
    void inject_thermal_electrons(GapElectrode electrode, double temperature_K,
                                  double current_density_A_m2, double duration_s,
                                  int macro_particles, std::uint64_t seed,
                                  GapSurface surface = GapSurface::Gap,
                                  double start_fraction = 0, double end_fraction = 1);
    void inject_richardson_electrons(GapElectrode electrode, double temperature_K,
        double work_function_eV, double richardson_A_m2_K2, double duration_s,
        int macro_particles, std::uint64_t seed, GapSurface surface = GapSurface::Gap);
    // Uniform geometric patches; CDF sampling integrates J(T,E) before sampling.
    void inject_surface_profile(GapElectrode electrode, GapSurface surface,
        const std::vector<double>& temperature_K, const std::vector<double>& current_A_m2,
        double duration_s, int macro_particles, std::uint64_t seed);
    void set_external_fields(std::function<PicVector(double,double,double)> electric,
                             std::function<PicVector(double,double,double)> magnetic);
    void set_external_electric_3d(std::function<PicVector(double,double,double,double)> correction);
    void set_potentials(double lower_V, double upper_V);
    // Exact linear Laplace rescaling of a solved empty domain; no particle step skipped.
    void set_empty_potentials(double lower_V, double upper_V);
    void solve_fields();
    void advance(double dt_s);
    // Empty-vacuum clock advance: no particles, no source, static potentials.
    void advance_empty(double dt_s);
    // Same frozen executable/ABI restart files, atomically replaced on save.
    void save_checkpoint(const std::string& path) const;
    void load_checkpoint(const std::string& path);
    [[nodiscard]] PicVector electric_field(double x_m, double y_m);
    [[nodiscard]] GapPicDiagnostics diagnostics() const;
    [[nodiscard]] const std::vector<GapParticle>& particles() const { return particles_; }
    [[nodiscard]] const mfem::GridFunction& potential() const { return potential_; }
    [[nodiscard]] const GapPicOptions& options() const { return options_; }
    [[nodiscard]] bool inside_conductor(double x_m, double y_m) const;
    // Oriented +x/+y edge charge transfers during the last advance (C).
    // Divide by that advance's dt for edge current (A). Injection is a separate
    // instantaneous boundary source, not an interior trajectory current.
    [[nodiscard]] const std::vector<double>& horizontal_transfers_C() const { return horizontal_transfers_; }
    [[nodiscard]] const std::vector<double>& vertical_transfers_C() const { return vertical_transfers_; }
private:
    int locate(double x, double y, mfem::IntegrationPoint& point) const;
    void drift(double dt);
    std::function<PicVector(double,double,double)> external_electric_, external_magnetic_;
    std::function<PicVector(double,double,double,double)> external_electric_3d_;
    PushObserver push_observer_;
    GapPicOptions options_;
    mfem::Mesh mesh_;
    mfem::H1_FECollection collection_;
    mfem::FiniteElementSpace space_;
    mfem::GridFunction potential_;
    std::vector<GapParticle> particles_;
    std::vector<double> horizontal_transfers_, vertical_transfers_;
    GapPicDiagnostics ledger_;
    bool fields_current_ = false;
    std::uint64_t next_particle_id_ = 1;
    std::unique_ptr<PicPoissonState> poisson_;
};
} // namespace raikou::mfem_backend
