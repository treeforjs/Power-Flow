#pragma once
#include <compare>
#include <cstdint>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace water_scope {
enum class Mode { legacy_parent, legacy_diagnostic, full_water_breakdown_candidate };
enum class PairStatus { absent, supported, entirely_unbounded };
struct Pair { std::uint32_t projectile, target; auto operator<=>(const Pair&) const = default; };
struct InventoryRow {
    std::uint64_t native_id;
    std::uint64_t collision_cell;
    std::uint32_t species;
    double physical_count;
};
struct StageIdentity {
    std::string input_state_sha256;
    std::string native_registry_sha256;
    std::string configuration_sha256;
    std::string expected_census_sha256;
    std::string expected_stage_binding_sha256;
    std::uint64_t native_row_count;
    std::uint64_t step;
    std::uint32_t stage;
    double dt_s;
    double collision_cell_volume_m3;
    bool all_registered_species_before_target_filter;
};
struct Gap { std::uint64_t cell; Pair pair; PairStatus status; };
struct Admission {
    bool enforced = false;
    bool admitted = false;
    bool existing_region_admission_required = false;
    bool network_complete = false;
    bool physical_model_validated = false;
    std::string scope_sha256, census_sha256, stage_binding_sha256;
    std::vector<Gap> gaps;
};

std::string sha256(std::span<const unsigned char> bytes);
std::string sha256(const std::string& text);
std::string census_digest(std::span<const InventoryRow> rows);

// Immutable minimal catalogue read from the actual table bytes. No rates are
// invented and no cross-section interpolation or energy support is evaluated.
class Catalogue {
public:
    Catalogue(std::span<const unsigned char> bytes, const std::string& expected_blob_sha256);
    const std::string& blob_sha256() const { return blob_sha_; }
    const std::string& registry_sha256() const { return registry_sha_; }
    const std::vector<std::string>& species() const { return names_; }
    const std::vector<Pair>& pairs() const { return pairs_; }
    std::uint32_t species_id(const std::string& name) const;
    PairStatus status(Pair pair) const;
    bool has_resolved_interval(Pair pair) const;
private:
    std::string blob_sha_, registry_sha_;
    std::vector<std::string> names_;
    std::vector<Pair> pairs_;
    std::vector<PairStatus> status_;
    std::vector<bool> resolved_;
};

// This minimum scope includes every existing table pair plus e+OH/O/H2. It is
// deliberately not described as a complete water/plasma reaction network.
std::string minimum_contract(const Catalogue& table, const std::string& native_registry_sha256,
                             const std::string& configuration_sha256);
class Guard {
public:
    Guard(Catalogue table, const std::string& contract, const std::string& expected_contract_sha256);
    // The native owner computes and retains this from its immutable stage
    // snapshot, then supplies it to check. It must not re-sign mutated inputs.
    std::string bind_stage(std::span<const InventoryRow> rows, const StageIdentity& identity) const;
    Admission check(std::span<const InventoryRow> rows, const StageIdentity& identity) const;
    const Catalogue& catalogue() const { return table_; }
private:
    Catalogue table_;
    std::string contract_sha_, native_registry_sha_, configuration_sha_;
    std::vector<Pair> required_;
};

// Existing parent/diagnostic paths bypass this opt-in gate without asserting
// full-network admission or requiring a new contract. Unknown modes refuse.
Admission admit(Mode mode, const Guard* guard, std::span<const InventoryRow> rows,
                const StageIdentity& identity);
} // namespace water_scope
