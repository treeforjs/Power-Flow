// forms.hpp -- native v7 producer forms: host C++ API (CPU reference).
// Definitions: API-007.txt section 1. The device kernels (gpu/forms_gpu.hpp)
// produce the same bits. Link: forms/src/forms_cpu.cpp, -lcrypto.
#pragma once
#include "pf_v7/forms_core.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace pf_v7 {
using Digest = std::array<unsigned char, 32>;
std::string hex(const Digest &);
Digest digest_of_words(const std::uint32_t *words8); // 8 SHA-256 state words -> bytes
void words_of_digest(const Digest &, std::uint32_t *words8);

// ---------------------------------------------------------------- records
// Canonical little-endian packed layouts (byte-identical to the native types
// noted; static_asserts in forms_cpu.cpp and in the native sources).
struct Carrier { // = raikou::mfem_backend::GapParticle, 80 bytes
  std::uint64_t id;
  double position_m[3], proper_velocity_m_s[3];
  double charge_C, mass_kg, weight;
};
struct Path { // = FullTimeMotionPath, 64 bytes
  std::uint64_t id;
  double field_start_m[3], field_displacement_m[3], flight_s;
};
struct Sink { // FullTimeSinkRecord fields widened to u64, 136 bytes
  std::uint64_t id, kind, surface, axis;
  double plane_m, fraction, segment_m[3], hit_point_m[3], incident_proper_velocity_m_s[3];
  double flight_s, discarded_flight_s;
};
struct Birth { // = FullTimeBirthFlight, 24 bytes
  std::uint64_t id;
  double birth_time_s, flight_s;
};
struct CensusRecord { // 96 bytes
  Carrier carrier;
  std::uint64_t packet_key, state_index;
};
static_assert(sizeof(Carrier) == pf7::endpoint_bytes && sizeof(Path) == pf7::path_bytes &&
              sizeof(Sink) == pf7::sink_bytes && sizeof(Birth) == pf7::birth_bytes &&
              sizeof(CensusRecord) == pf7::census_bytes);

// ---------------------------------------------------------------- (a) exact sums
// Exact sum of exact products of binary64 values; any order, any partition.
class ExactSum {
public:
  using Digits = std::array<std::int64_t, pf7::acc_digits>;
  ExactSum() { digits_.fill(0); }
  void add_product(double a, double b);   // exact a*b
  void add(double x) { add_product(x, 1.0); }
  void merge(const ExactSum &other);      // exact
  // Correctly rounded value (nearest, ties to even; exact zero -> +0.0).
  // status: pf7::sum_* bits (nonfinite input, overflow).
  double value(std::uint32_t *status = nullptr) const;
  // Normalised digits (0..132 in [0, 2^32), 133 signed): the interchange
  // form of a partial sum (equal partial sums have equal digits).
  Digits digits() const;
  static ExactSum from_digits(const Digits &, std::uint32_t status = 0, std::uint64_t terms = 0);
  std::uint32_t status() const noexcept { return status_; }
  std::uint64_t terms() const noexcept { return terms_; }

private:
  Digits digits_;
  std::uint64_t terms_ = 0, pending_ = 0;
  std::uint32_t status_ = 0;
};
struct SumValue {
  double value = 0;       // rounded exact sum (field work: dt x rounded sum)
  double sum = 0;         // rounded exact sum before the dt factor
  std::uint32_t status = 0;
  bool ok() const noexcept { return status == 0; }
};
// kinetic_change_J over the census (ascending native ID): census carrier i
// ends at endpoints[j] (same ID; survivors) or at sinks[s] (same ID; incident
// proper velocity). Term_i = w_i x RN(KE(m_i, u_end) - KE(m_i, u_start)).
// Missing or extra IDs set pf7::sum_census_missing.
SumValue kinetic_change_J(std::span<const Carrier> census, std::span<const Carrier> endpoints,
                          std::span<const Sink> sinks);
// Same with the end proper velocity of every census carrier (3 doubles each).
SumValue kinetic_change_J_aligned(std::span<const Carrier> census, std::span<const double> u_end);
// RN(sum a_i b_i) exactly; field_work_J = RN(dt x RN(sum E_e J_e)).
SumValue exact_dot(std::span<const double> a, std::span<const double> b);
SumValue field_work_J(double dt_s, std::span<const double> E, std::span<const double> J);
double kinetic_energy_J(double mass_kg, const double *u); // pic_kinetic_energy_J

// ---------------------------------------------------------------- (b) block lists
struct BlockList {
  std::string root;           // hex root
  Digest root_digest{};
  std::vector<Digest> leaves; // ascending leaf index
  std::uint64_t records = 0;
};
// Section block list over n canonical records (record size of the section).
BlockList block_list(std::uint32_t section, const void *records, std::uint64_t n);
// Root from leaf digests (verification of device leaves).
std::string block_root(std::uint32_t section, std::uint64_t records, std::span<const Digest> leaves);
// Motion record v7 (API-007.txt 1.2): every field except the per-carrier
// arrays, which enter through their section roots.
inline constexpr char motion_v7_schema[] = "canonical-full-time-total-field-motion-v7-merkle";
inline constexpr char motion_v7_coverage[] = "all_native_carriers_births_and_first_hit_sinks";
struct MotionHeaderV7 {
  std::uint64_t accepted_step = 0, staged_arena_revision = 0, trial_stage = 0;
  double begin_s = 0, end_s = 0, dt_s = 0;
  std::string staged_arena_sha256, trial_field_sha256, mesh_sha256, material_sha256,
      species_sha256, common_clock_sha256, kinetic_policy_sha256;
  // kinetic_change, field_work, work_defect, continuity, plasma_dt,
  // maximum_gyro, maximum_displacement_cells
  std::array<double, 7> budget{};
  std::string census_sha256, coverage = motion_v7_coverage;
  std::string current_A_m_sha256, endpoint_charge_C_sha256, initial_charge_C_sha256;
  std::array<double, 12> kinetic_policy{};
  std::array<double, 3> native_to_field_m{};
  std::string neutral_motion_owner = "same_full_time_q0_drift";
  bool native_DKD_called = false, native_Poisson_called = false;
  // Interval bound of the resident-neutral stream used by the producer
  // (native v7 configuration; default 1e-14). Must equal native's value.
  double neutral_stream_dt_max_s = 1e-14;
};
struct MotionSectionsV7 {
  std::uint64_t endpoints = 0, paths = 0, sinks = 0, births = 0; // record counts
  std::string endpoints_root, paths_root, sinks_root, births_root;
};
std::string motion_record_v7_sha256(const MotionHeaderV7 &, const MotionSectionsV7 &);
MotionSectionsV7 motion_sections_v7(std::span<const Carrier> endpoints, std::span<const Path> paths,
                                    std::span<const Sink> sinks, std::span<const Birth> births);
// motion-array-binary64-v1 digest (Hall array_hash; unchanged form).
std::string array_sha256(std::span<const double>);

// ---------------------------------------------------------------- (c) census
inline constexpr char census_v2_tag[] = "canonical-motion-census-v2-merkle";
// states_sha256 = SHA256(text "canonical-motion-census-v2/states", u64 S,
// S x text name), names in ascending byte order (= state index order).
std::string census_states_sha256(std::span<const std::string> states);
struct CensusCheck {
  std::uint32_t code = pf7::census_ok; // pf7::CensusCode of the first failing record
  std::uint64_t index = 0;
  bool ok() const noexcept { return code == pf7::census_ok; }
};
CensusCheck check_census(std::span<const CensusRecord>, std::uint64_t states);
// census_sha256 = SHA256(text tag, u64 N, text states_sha256, text carriers root).
std::string census_v2_sha256(std::uint64_t carriers, const std::string &states_sha256,
                             const std::string &carriers_root);
// Checks and digests one census; throws std::invalid_argument on a failed check.
std::string census_v2_sha256(std::span<const CensusRecord>, std::span<const std::string> states);
// Stable ascending-ID permutation (canonical census order).
std::vector<std::uint64_t> census_order(std::span<const CensusRecord>);
// Census records of host entries (packet ID text and state text per carrier,
// any order): packet key = smallest native ID among the carriers of the
// packet ID; state index = rank of the state text in `states` (sorted,
// unique). Output in ascending native ID.
struct HostCensusEntry {
  std::string_view packet_id, state;
  const Carrier *carrier;
};
std::vector<CensusRecord> census_records(std::span<const HostCensusEntry>,
                                         std::span<const std::string> states);
// Census carriers block list from rows + per-row keys/states (no copy of the
// records).
BlockList census_block_list(std::span<const Carrier> rows, std::span<const std::uint64_t> keys,
                            std::span<const std::uint64_t> states);

// ---------------------------------------------------------------- neutral stream
// Resident-neutral stream with a configured interval bound (API-007.txt 1.4).
// With dt_max_s = 1e-14 it reproduces the sealed Hall 8dddec0
// prepare_neutral_stream bit for bit (values, refusal texts and order).
inline constexpr double neutral_stream_dt_max_default_s = 1e-14;
struct KineticPolicy { // native-finite-target-kinetic-policy-v1 operands
  double width_m = 0, gap_m = 0, short_z_m = 0, feed_z_m = 0;
  double lower[4]{}, upper[4]{}; // left, right, bottom, top
};
struct NeutralStreamResult {
  Carrier endpoint{};
  pf7::Hit hit{-1, 1, -1, 0};
  double displacement_m[3]{};
  double begin_s = 0, end_s = 0, dt_s = 0, flight_s = 0, remaining_s = 0, incident_kinetic_J = 0;
};
// Throws std::invalid_argument with the sealed refusal texts.
NeutralStreamResult neutral_stream(const KineticPolicy &, const Carrier &, double begin_s, double end_s,
                                   double dt_s, double dt_max_s = neutral_stream_dt_max_default_s);
// Sealed kinetic policy guard (texts of kinetic_boundary.cpp policy_guard).
void require_kinetic_policy(const KineticPolicy &);

// ---------------------------------------------------------------- work tolerance
// motion-work-gross-v1 (amendment of API-007.txt 1.1; specification
// hall-corrective-v3.2-20261031/MOTION-WORK-GROSS-V1.txt):
//   tol = RN(RN(RN(1e-20 + RN(1e-10 RN(G_K + G_W))) + RN(2^-47 C_K)) + RN(RN(8e-18 E^) Q))
// over the charged census carriers (q != 0): G_K = RN(sum w |d|), C_K =
// RN(sum w K0 + w RN(K0 + d)), Q = RN(sum |q| w); field: G_W = RN(dt RN(sum
// |E_e| |J_e|)), E^ = max |E_e| over the edges with an endpoint node on an edge
// carrying J != 0 (0 when none). Exact sums, exact maximum.
// Hall axial mesh (nx = ny = n, nz): edge numbering of E and J.
struct HallAxialMesh {
  std::uint64_t n = 0, nz = 0;
  std::uint64_t ne() const noexcept { return 2 * n * (n + 1) * (nz + 1) + (n + 1) * (n + 1) * nz; }
  std::uint64_t nb() const noexcept { return 2 * n * (n + 1) * nz + n * n * (nz + 1); }
  std::uint64_t nn() const noexcept { return (n + 1) * (n + 1) * (nz + 1); }
  void edge_nodes(std::uint64_t e, std::uint64_t &a, std::uint64_t &b) const noexcept;
};
// The unique (n, nz) with ne() == edges and nb() == faces (n in [1, 4096]);
// false when none or several exist.
bool hall_axial_mesh_of(std::uint64_t edges, std::uint64_t faces, HallAxialMesh &);
// total-midpoint-E-B-binary64-v1 digest (Hall total_field_sha256; trial_field_sha256).
std::string total_field_sha256(std::span<const double> E, std::span<const double> B);
struct WorkScales {
  double kinetic_gross_J = 0, kinetic_content_J = 0, field_work_gross_J = 0, charge_abs_C = 0,
         field_bound_V_m = 0, tolerance_J = 0;
  std::uint32_t status = 0; // pf7 sum bits of the four sums; E^ nonfinite -> sum_nonfinite_input
};
class WorkCarrierSums {
public:
  // weight w > 0, charge q, K0 = KE(m, u_start), d = RN(KE(m, u_end) - K0)
  void add(double weight, double charge, double K0, double d);
  void merge(const WorkCarrierSums &);
  ExactSum gk, ck, q;
};
void work_field_scales(WorkScales &, double dt_s, std::span<const double> E, std::span<const double> J,
                       const HallAxialMesh &);
void work_finish(WorkScales &, const WorkCarrierSums &); // carrier sums + tolerance
double work_tolerance(const WorkScales &);
bool work_admits(const WorkScales &, double defect);
} // namespace pf_v7
