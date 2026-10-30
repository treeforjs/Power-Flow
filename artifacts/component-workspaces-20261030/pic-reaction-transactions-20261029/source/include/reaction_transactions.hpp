#pragma once

#include "raikou/json.hpp"
#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// An isolated manufactured-event adapter, not a collision/rate solver. It never
// calls GapPic::add_particle, advances PIC, or evolves an electron F0.
namespace raikou::pic_reactions {
using Json = nlohmann::json;
using Vec = std::array<double, 3>;
inline constexpr double electron_charge_C = 1.602176634e-19;
inline constexpr double speed_of_light_m_s = 299792458.;
inline constexpr double ledger_roundoff_factor = 256.;

struct Species {
    std::string id;
    int hydrogen = 0, oxygen = 0, charge_e = 0;
    double mass_kg = 0;
};
struct State {
    std::string id, species;
    double reference_J = 0, internal_J = 0;
};
struct Channel {
    std::string id, target_state;
    bool scatters_primary = true;
    std::vector<std::string> product_states; // slot0 is the retained primary if scattering
};
// Construct once, return const definitions only, and bind the complete catalog
// digest into every event/restart. Mutation requires a different catalog/run.
class Catalog {
  public:
    explicit Catalog(Json definition);
    const Species &species(const std::string &id) const;
    const State &state(const std::string &id) const;
    const Channel &channel(const std::string &id) const;
    const Species &species_of(const std::string &state) const;
    const std::string &digest() const {
        return digest_;
    }
    const Json &definition() const {
        return definition_;
    }

  private:
    Json definition_;
    std::string digest_;
    std::map<std::string, Species> species_;
    std::map<std::string, State> states_;
    std::map<std::string, Channel> channels_;
};

struct Packet {
    std::string id, state, source_id, birth_event_id, last_event_id;
    double weight = 0;
    Vec position_m{}, proper_velocity_m_s{};
};
// A gas-parent pool, already released by its separate surface owner. A reaction
// debits this one lot; it must not debit the adsorbed source a second time.
struct NeutralLot {
    std::string id, state, source_id;
    double remaining_weight = 0;
    Vec position_m{}, proper_velocity_m_s{};
};
struct Product {
    std::string id, state;
    double weight = 0;
    Vec proper_velocity_m_s{};
};
struct ForcedEvent {
    std::string id, channel, parent_id, primary_id, catalog_digest;
    std::string manufactured_partition; // explicit arbitrary fixture provenance
    std::uint64_t sequence = 0;
    double time_s = 0, weight = 0;
    std::vector<Product> products;
};
struct Genesis {
    std::string run_id, coupling_prefix, electron_owner = "PIC";
    std::vector<Packet> packets;
    std::vector<NeutralLot> parents;
    // Imported external source accounting is distinct from internal reactions.
    double authored_charge_C = 0, authored_kinetic_J = 0;
};
struct Snapshot {
    std::uint64_t revision = 0, next_event_sequence = 1;
    double time_s = 0;
    std::map<std::string, Packet> packets;
    std::map<std::string, NeutralLot> parents;
    double authored_charge_C = 0, authored_kinetic_J = 0;
    std::map<std::string, std::string> accepted_events; // immutable ID -> payload SHA256
    std::set<std::string> ever_used_packet_ids;         // consumed IDs can never be recycled
    std::vector<Json> journal;
    std::string chain_digest;
};
struct ParticleDelta {
    std::string event_id, removed_primary;
    std::vector<Packet> inserted; // includes same-ID primary replacement, if any
    Json conservation;
};
enum class ManufacturedFault {
    None,
    AfterParentDebit,
    AfterPrimaryRemoval,
    AfterProductInsertion,
    BeforePublication
};
class PreparedReaction {
  public:
    const ParticleDelta &delta() const noexcept {
        return delta_;
    }
    std::shared_ptr<const Snapshot> input() const noexcept {
        return input_;
    }
    std::shared_ptr<const Snapshot> output() const noexcept {
        return output_;
    }

  private:
    friend class Transactions;
    PreparedReaction(std::shared_ptr<const Snapshot> input, std::shared_ptr<const Snapshot> output,
                     ParticleDelta delta)
        : input_(std::move(input)), output_(std::move(output)), delta_(std::move(delta)) {}
    std::shared_ptr<const Snapshot> input_, output_;
    ParticleDelta delta_;
};

class Transactions {
  public:
    Transactions(Catalog catalog, Genesis genesis);
    Transactions(const Transactions &) = delete;
    Transactions &operator=(const Transactions &) = delete;
    Transactions(Transactions &&) noexcept = default;
    Transactions &operator=(Transactions &&) noexcept = default;
    // One owner thread. All potentially throwing work occurs in a private copy;
    // the only publication is a no-throw shared_ptr swap after all checks.
    ParticleDelta apply(const ForcedEvent &event, ManufacturedFault fault = ManufacturedFault::None);
    PreparedReaction prepare(const ForcedEvent &event,
                             ManufacturedFault fault = ManufacturedFault::None) const;
    bool can_publish(const PreparedReaction &prepared) const noexcept;
    // Local owner publication only. A native multi-owner commit must preflight
    // every participant under the same exclusive accepted-interval barrier.
    bool try_publish(const PreparedReaction &prepared) noexcept;
    std::shared_ptr<const Snapshot> snapshot() const {
        return snapshot_;
    }
    const Catalog &catalog() const {
        return catalog_;
    }
    const Genesis &genesis() const {
        return genesis_;
    }
    Json checkpoint_document() const;
    void save_checkpoint(const std::string &path) const;
    static Transactions restore(const Catalog &catalog, const Json &checkpoint,
                                const std::string &expected_coupling_prefix);
    static Transactions load_checkpoint(const Catalog &catalog, const std::string &path,
                                        const std::string &expected_coupling_prefix);

  private:
    Catalog catalog_;
    Genesis genesis_;
    std::shared_ptr<const Snapshot> snapshot_;
};

double kinetic_J(const Species &species, const Vec &proper_velocity);
Json event_json(const ForcedEvent &event);
Json snapshot_json(const Snapshot &snapshot);
std::string sha256(const std::string &bytes);
} // namespace raikou::pic_reactions
