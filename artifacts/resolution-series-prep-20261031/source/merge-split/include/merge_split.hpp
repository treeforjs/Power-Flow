// merge_split.hpp -- CPU reference of the conservative merge/split transaction
// (schema merge-split-transaction-v1). Arithmetic is ms_core.h, shared with the
// NVRTC kernels; this file adds the orchestration (classification, canonical
// order, segments, group formation, kinematics, ID allocation, emission), the
// exact ledger, the canonical encoding and its SHA-256 digest, and the
// partitioned execution used to show rank independence.
#pragma once
#include "ms_core.h"

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace merge_split {

struct Error : std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

// Carrier set (structure of arrays, any order). IDs must be distinct and nonzero.
struct Carriers {
    std::vector<std::uint64_t> id;
    std::vector<std::uint32_t> species, flags;
    std::array<std::vector<double>, 3> x, u;   // native position (m), proper velocity (m/s)
    std::vector<double> w;                     // represented physical particles
    std::size_t size() const noexcept { return id.size(); }
    void reserve(std::size_t n);
    void push(std::uint64_t id, std::uint32_t species, std::uint32_t flags, const double x[3], const double u[3], double w);
};

// Exact value in units of 2^-1074 (normalized superaccumulator digits).
struct Exact {
    std::array<std::int64_t, MS_ACC_DIGITS> d{};
    static Exact of(MsAcc a);                   // normalizes a copy
    Exact operator-(const Exact&) const;
    Exact operator+(const Exact&) const;
    bool operator==(const Exact& o) const noexcept { return d == o.d; }
    bool zero() const noexcept;
    bool negative() const noexcept { return d[MS_ACC_DIGITS - 1] < 0; }
    double rounded() const;                     // round to nearest, ties to even
    std::vector<double> expansion() const;      // canonical: x0 = RN(v), x1 = RN(v - x0), ... (exact)
    std::string hex() const;                    // signed hexadecimal integer of units 2^-1074
};

// Per-species exact ledger of one transaction. Inputs: carriers consumed by
// accepted groups; outputs: the carriers created. Momentum is w m u (kg m/s),
// energy w K with K = pic_kinetic_energy_J (J), weight w (count). Every term is
// an exact product (TwoProduct words), so the sums are exact integers of 2^-1074.
struct SpeciesLedger {
    Exact w_in, w_out, e_in, e_out, e_res;
    std::array<Exact, 3> p_in, p_out, p_res;
    double e_res_rn = 0;                        // booked merge_energy_residual (RN of the exact value)
    std::array<double, 3> p_res_rn{};           // booked merge_momentum_residual
    std::vector<double> e_res_expansion;        // exact booked value as a canonical binary64 expansion
    std::array<std::vector<double>, 3> p_res_expansion;
    bool count_equal = false;                   // w_out == w_in exactly (ExactParticleCount equality)
};

struct Counters {
    std::uint64_t carriers = 0, segments = 0, groups = 0, merge_groups = 0, split_groups = 0;
    std::uint64_t accepted_merge = 0, accepted_split = 0, formation_refused = 0;
    std::array<std::uint64_t, MS_STATUS_KINDS> status{};
    std::uint64_t centroid_fallback = 0, isotropic = 0;
    std::uint64_t inputs_consumed = 0, outputs = 0, relocations = 0;
};

// One record per input of an accepted merge group, in transaction order (group,
// then canonical member order): straight relocation path in the Hall field frame.
struct Relocation {
    std::vector<std::uint64_t> id;          // input native ID (retired by the transaction)
    std::vector<std::uint32_t> group;       // group rank
    std::array<std::vector<double>, 3> start, D;   // start = x + shift; D = (c + shift) - (x + shift)
    std::vector<double> charge;             // w q (binary64, as the Hall path deposit operand)
    std::size_t size() const noexcept { return id.size(); }
};

struct Result {
    std::vector<MsGroup> groups;            // canonical order (cell, species, bin, first input ID)
    std::vector<std::uint64_t> group_first_id;
    std::vector<MsMergeResult> merge;       // per group (zero for splits and skipped groups)
    Carriers outputs;                       // IDs base + rank, canonical output order
    std::vector<std::uint64_t> removed;     // ascending native ID
    Relocation reloc;
    std::vector<SpeciesLedger> ledger;      // per species
    Counters counters;
    std::uint64_t next_id_after = 0;
    std::string params_sha256, digest_sha256;
};

struct PartitionPlan {
    unsigned parts = 1;
    bool shuffled = false;                  // false: contiguous ID blocks; true: hashed owner + shuffled local order
    std::uint64_t shuffle_seed = 0;
};

// Validates parameters (throws Error).
void validate(const MsParams&);
// Transaction binding (seed, step, base ID and policy) and the policy identity alone.
std::string params_sha256(const MsParams&);
std::string policy_sha256(const MsParams&);
// Mesh fields of the Hall axial mesh Mesh(n, axial, side, air) (axial_mesh.hpp arithmetic:
// dx = side/n, dz = depth/axial, zmin = -depth/2 - air dz), kinetic cells k in
// [air, air + axial), native -> field shift (-side/2, -side/2, 0); nsub = 1, no species.
MsParams hall_mesh_params(int n, int axial = 4, double side = .009, int air = 1, double depth = .025);
// Species s with mass and charge (mc2 = m c c); policy fields zero (caller sets them).
void set_species(MsParams&, unsigned s, double mass_kg, double charge_C);
// The transaction. The result depends only on the carrier set and the parameters.
Result transact(const MsParams&, const Carriers&, const PartitionPlan& = {});
// Post-transaction census: survivors (input order kept) then outputs.
Carriers apply(const Carriers& census, const Result&);
// Canonical byte encoding of a result (digest input) and its SHA-256.
std::string encode(const Result&);
std::string sha256_hex(const std::string&);

// Exact count/charge per species of an arbitrary carrier set (superaccumulator).
std::vector<Exact> exact_weights(const MsParams&, const Carriers&);

namespace detail {
// Counters derived from the groups (kinds, status, flags, outputs, relocations, next ID).
void fill_group_counters(Result&, const MsParams&);
// Ledger from per-species exact totals, layout per species:
// W_in W_out P_in[3] P_out[3] E_in E_out (the GPU ledger kernel uses the same layout).
void fill_ledger(Result&, const MsParams&, const std::vector<Exact>& totals);
}  // namespace detail

}  // namespace merge_split
