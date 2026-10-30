#include "reaction_transactions.hpp"
#include <openssl/sha.h>
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>
#include <system_error>
#include <type_traits>
#include <fcntl.h>
#include <unistd.h>

namespace raikou::pic_reactions {
namespace {
constexpr std::size_t maximum_records = 4096, maximum_checkpoint_bytes = 8 * 1024 * 1024;
static_assert(std::is_nothrow_move_constructible_v<ParticleDelta>);
void require(bool pass, const std::string &message) {
    if (!pass)
        throw std::invalid_argument(message);
}
void finite(double x, const std::string &what) {
    require(std::isfinite(x), what + " is nonfinite");
}
void identifier(const std::string &id) {
    require(!id.empty() && id.size() <= 256 &&
                std::all_of(id.begin(), id.end(), [](unsigned char c) { return c > 32 && c < 127; }),
            "invalid immutable ID");
}
void keys(const Json &j, std::initializer_list<const char *> expected) {
    require(j.is_object() && j.size() == expected.size(), "unexpected record fields");
    for (auto key : expected)
        require(j.contains(key), std::string("missing field ") + key);
}
void vector_finite(const Vec &v) {
    for (double x : v)
        finite(x, "vector operand");
}
double norm2(const Vec &u) {
    return u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
}
void close(double a, double b, double scale, const char *message) {
    finite(a, message);
    finite(b, message);
    finite(scale, message);
    const double limit =
        ledger_roundoff_factor * std::numeric_limits<double>::epsilon() * std::max(scale, 1e-300);
    require(std::abs(a - b) <= limit, message);
}
Json packet_json(const Packet &p) {
    return {{"id", p.id},
            {"state", p.state},
            {"source_id", p.source_id},
            {"birth_event_id", p.birth_event_id},
            {"last_event_id", p.last_event_id},
            {"weight", p.weight},
            {"position_m", p.position_m},
            {"proper_velocity_m_s", p.proper_velocity_m_s}};
}
Packet packet_from(const Json &j) {
    keys(j, {"id", "state", "source_id", "birth_event_id", "last_event_id", "weight", "position_m",
             "proper_velocity_m_s"});
    return {j.at("id"),
            j.at("state"),
            j.at("source_id"),
            j.at("birth_event_id"),
            j.at("last_event_id"),
            j.at("weight"),
            j.at("position_m").get<Vec>(),
            j.at("proper_velocity_m_s").get<Vec>()};
}
Json lot_json(const NeutralLot &p) {
    return {{"id", p.id},
            {"state", p.state},
            {"source_id", p.source_id},
            {"remaining_weight", p.remaining_weight},
            {"position_m", p.position_m},
            {"proper_velocity_m_s", p.proper_velocity_m_s}};
}
NeutralLot lot_from(const Json &j) {
    keys(j, {"id", "state", "source_id", "remaining_weight", "position_m", "proper_velocity_m_s"});
    return {j.at("id"),
            j.at("state"),
            j.at("source_id"),
            j.at("remaining_weight"),
            j.at("position_m").get<Vec>(),
            j.at("proper_velocity_m_s").get<Vec>()};
}
Json genesis_json(const Genesis &g) {
    Json p = Json::array(), n = Json::array();
    for (const auto &v : g.packets)
        p.push_back(packet_json(v));
    for (const auto &v : g.parents)
        n.push_back(lot_json(v));
    return {{"run_id", g.run_id},
            {"coupling_prefix", g.coupling_prefix},
            {"electron_owner", g.electron_owner},
            {"packets", p},
            {"parents", n},
            {"authored_charge_C", g.authored_charge_C},
            {"authored_kinetic_J", g.authored_kinetic_J}};
}
Genesis genesis_from(const Json &j) {
    keys(j, {"run_id", "coupling_prefix", "electron_owner", "packets", "parents", "authored_charge_C",
             "authored_kinetic_J"});
    Genesis g;
    g.run_id = j.at("run_id");
    g.coupling_prefix = j.at("coupling_prefix");
    g.electron_owner = j.at("electron_owner");
    g.authored_charge_C = j.at("authored_charge_C");
    g.authored_kinetic_J = j.at("authored_kinetic_J");
    require(j.at("packets").is_array() && j.at("parents").is_array() &&
                j.at("packets").size() + j.at("parents").size() <= maximum_records,
            "checkpoint record bound");
    for (const auto &v : j.at("packets"))
        g.packets.push_back(packet_from(v));
    for (const auto &v : j.at("parents"))
        g.parents.push_back(lot_from(v));
    return g;
}
ForcedEvent event_from(const Json &j) {
    keys(j, {"id", "channel", "parent_id", "primary_id", "catalog_digest", "manufactured_partition",
             "sequence", "time_s", "weight", "products"});
    ForcedEvent e;
    e.id = j.at("id");
    e.channel = j.at("channel");
    e.parent_id = j.at("parent_id");
    e.primary_id = j.at("primary_id");
    e.catalog_digest = j.at("catalog_digest");
    e.manufactured_partition = j.at("manufactured_partition");
    e.sequence = j.at("sequence");
    e.time_s = j.at("time_s");
    e.weight = j.at("weight");
    require(j.at("products").is_array() && j.at("products").size() <= 32, "product record bound");
    for (const auto &p : j.at("products")) {
        keys(p, {"id", "state", "weight", "proper_velocity_m_s"});
        e.products.push_back(
            {p.at("id"), p.at("state"), p.at("weight"), p.at("proper_velocity_m_s").get<Vec>()});
    }
    return e;
}
struct Moments {
    double H = 0, O = 0, Q = 0, K_e = 0, K_ion = 0, K_neutral = 0, internal = 0, reference = 0;
    Vec momentum{}, absolute_momentum{};
    double total() const {
        return K_e + K_ion + K_neutral + internal + reference;
    }
    double absolute_energy = 0;
    void add(const Catalog &c, const std::string &state, double w, const Vec &u) {
        const auto &s = c.species_of(state);
        const auto &st = c.state(state);
        H += w * s.hydrogen;
        O += w * s.oxygen;
        Q += w * s.charge_e;
        double k = w * kinetic_J(s, u);
        if (s.id == "e")
            K_e += k;
        else if (s.charge_e)
            K_ion += k;
        else
            K_neutral += k;
        internal += w * st.internal_J;
        reference += w * st.reference_J;
        absolute_energy += std::abs(k) + std::abs(w * st.internal_J) + std::abs(w * st.reference_J);
        for (int i = 0; i < 3; ++i) {
            double p = w * s.mass_kg * u[i];
            momentum[i] += p;
            absolute_momentum[i] += std::abs(p);
        }
    }
    Json json() const {
        return {{"H_nuclei", H},
                {"O_nuclei", O},
                {"charge_e", Q},
                {"momentum_kg_m_s", momentum},
                {"electron_kinetic_J", K_e},
                {"ion_kinetic_J", K_ion},
                {"neutral_kinetic_J", K_neutral},
                {"internal_J", internal},
                {"reference_J", reference},
                {"total_event_energy_J", total()}};
    }
};
void validate_packet(const Catalog &c, const Packet &p) {
    for (const auto &id : {p.id, p.state, p.source_id, p.birth_event_id, p.last_event_id})
        identifier(id);
    finite(p.weight, "packet weight");
    require(p.weight > 0, "packet weight must be positive");
    vector_finite(p.position_m);
    vector_finite(p.proper_velocity_m_s);
    finite(p.weight * kinetic_J(c.species_of(p.state), p.proper_velocity_m_s), "weighted packet energy");
}
void fault_at(ManufacturedFault requested, ManufacturedFault here) {
    if (requested == here)
        throw std::runtime_error("manufactured rollback fault");
}
} // namespace

std::string sha256(const std::string &bytes) {
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char *>(bytes.data()), bytes.size(), digest);
    static constexpr char hex[] = "0123456789abcdef";
    std::string out;
    out.reserve(64);
    for (auto b : digest) {
        out.push_back(hex[b >> 4]);
        out.push_back(hex[b & 15]);
    }
    return out;
}
double kinetic_J(const Species &s, const Vec &u) {
    vector_finite(u);
    const double u2 = norm2(u);
    finite(u2, "proper speed squared");
    const double gamma = std::sqrt(1 + u2 / (speed_of_light_m_s * speed_of_light_m_s));
    const double k = s.mass_kg * u2 / (gamma + 1);
    finite(k, "relativistic kinetic energy");
    return k;
}
Catalog::Catalog(Json d) : definition_(std::move(d)) {
    keys(definition_, {"schema", "status", "electron_owner", "species", "states", "channels"});
    require(definition_.at("schema") == "manufactured-reaction-catalog/v1" &&
                definition_.at("status") == "manufactured_forced_only" &&
                definition_.at("electron_owner") == "PIC",
            "manufactured PIC-owned catalog required");
    for (const auto &j : definition_.at("species")) {
        for (const char *field : {"H", "O", "charge_e"})
            require(j.at(field).is_number_integer(), "integer species counts required");
        require(j.at("H") >= 0 && j.at("H") <= 64 && j.at("O") >= 0 && j.at("O") <= 64 &&
                    j.at("charge_e") >= -4 && j.at("charge_e") <= 4,
                "integer species count range");
        keys(j, {"id", "H", "O", "charge_e", "mass_kg"});
        Species s{j.at("id"), j.at("H"), j.at("O"), j.at("charge_e"), j.at("mass_kg")};
        identifier(s.id);
        finite(s.mass_kg, "species mass");
        require(s.hydrogen >= 0 && s.hydrogen <= 64 && s.oxygen >= 0 && s.oxygen <= 64 &&
                    s.mass_kg > 0 && s.charge_e >= -4 && s.charge_e <= 4,
                "invalid species");
        require(species_.emplace(s.id, s).second, "duplicate species ID");
    }
    const auto &electron = species("e");
    require(electron.hydrogen == 0 && electron.oxygen == 0 && electron.charge_e == -1,
            "electron catalog mismatch");
    for (const auto &j : definition_.at("states")) {
        keys(j, {"id", "species", "reference_J", "internal_J"});
        State s{j.at("id"), j.at("species"), j.at("reference_J"), j.at("internal_J")};
        identifier(s.id);
        species(s.species);
        finite(s.reference_J, "reference energy");
        finite(s.internal_J, "internal energy");
        require(s.internal_J >= 0, "negative internal energy");
        require(states_.emplace(s.id, s).second, "duplicate state ID");
    }
    for (const auto &j : definition_.at("channels")) {
        keys(j, {"id", "target_state", "scatters_primary", "product_states"});
        Channel ch{j.at("id"), j.at("target_state"), j.at("scatters_primary"),
                   j.at("product_states").get<std::vector<std::string>>()};
        identifier(ch.id);
        const auto &target = species_of(ch.target_state);
        require(target.charge_e == 0 && target.id != "e", "one neutral target required");
        require(ch.product_states.size() >= 2 && ch.product_states.size() <= 32,
                "bounded explicit product list required");
        int H = 0, O = 0, Q = 0;
        for (const auto &st : ch.product_states) {
            const auto &s = species_of(st);
            H += s.hydrogen;
            O += s.oxygen;
            Q += s.charge_e;
        }
        require(H == target.hydrogen && O == target.oxygen && Q == target.charge_e - 1,
                "channel H/O/charge imbalance");
        require(!ch.scatters_primary || species_of(ch.product_states[0]).id == "e",
                "scattered primary slot must be electron");
        require(ch.scatters_primary ||
                    std::none_of(ch.product_states.begin(), ch.product_states.end(),
                                 [&](const auto &st) { return species_of(st).id == "e"; }),
                "attachment cannot also return primary electron");
        require(channels_.emplace(ch.id, ch).second, "duplicate channel ID");
    }
    require(species_.size() <= 128 && states_.size() <= 256 && channels_.size() <= 256,
            "catalog size bound");
    digest_ = sha256(definition_.dump());
}
const Species &Catalog::species(const std::string &id) const {
    return species_.at(id);
}
const State &Catalog::state(const std::string &id) const {
    return states_.at(id);
}
const Channel &Catalog::channel(const std::string &id) const {
    return channels_.at(id);
}
const Species &Catalog::species_of(const std::string &id) const {
    return species(state(id).species);
}

Json event_json(const ForcedEvent &e) {
    Json products = Json::array();
    for (const auto &p : e.products)
        products.push_back({{"id", p.id},
                            {"state", p.state},
                            {"weight", p.weight},
                            {"proper_velocity_m_s", p.proper_velocity_m_s}});
    return {{"id", e.id},
            {"channel", e.channel},
            {"parent_id", e.parent_id},
            {"primary_id", e.primary_id},
            {"catalog_digest", e.catalog_digest},
            {"manufactured_partition", e.manufactured_partition},
            {"sequence", e.sequence},
            {"time_s", e.time_s},
            {"weight", e.weight},
            {"products", products}};
}
Json snapshot_json(const Snapshot &s) {
    Json p = Json::array(), n = Json::array();
    for (const auto &[id, v] : s.packets)
        p.push_back(packet_json(v));
    for (const auto &[id, v] : s.parents)
        n.push_back(lot_json(v));
    return {{"revision", s.revision},
            {"next_event_sequence", s.next_event_sequence},
            {"time_s", s.time_s},
            {"packets", p},
            {"parents", n},
            {"authored_charge_C", s.authored_charge_C},
            {"authored_kinetic_J", s.authored_kinetic_J},
            {"accepted_events", s.accepted_events},
            {"ever_used_packet_ids", s.ever_used_packet_ids},
            {"chain_digest", s.chain_digest}};
}
Transactions::Transactions(Catalog c, Genesis g) : catalog_(std::move(c)), genesis_(std::move(g)) {
    identifier(genesis_.run_id);
    identifier(genesis_.coupling_prefix);
    require(genesis_.electron_owner == "PIC", "PIC is sole electron owner; F0 evolution forbidden");
    require(genesis_.packets.size() + genesis_.parents.size() <= maximum_records, "genesis size bound");
    auto s = std::make_shared<Snapshot>();
    finite(genesis_.authored_charge_C, "authored charge");
    finite(genesis_.authored_kinetic_J, "authored energy");
    require(genesis_.authored_kinetic_J >= 0, "negative authored energy");
    s->authored_charge_C = genesis_.authored_charge_C;
    s->authored_kinetic_J = genesis_.authored_kinetic_J;
    std::set<std::string> ids;
    for (const auto &p : genesis_.packets) {
        validate_packet(catalog_, p);
        require(ids.insert(p.id).second, "duplicate packet/parent ID");
        s->packets.emplace(p.id, p);
        s->ever_used_packet_ids.insert(p.id);
    }
    for (const auto &p : genesis_.parents) {
        identifier(p.id);
        identifier(p.source_id);
        finite(p.remaining_weight, "parent inventory");
        require(p.remaining_weight >= 0 && catalog_.species_of(p.state).charge_e == 0,
                "invalid neutral inventory");
        vector_finite(p.position_m);
        kinetic_J(catalog_.species_of(p.state), p.proper_velocity_m_s);
        require(ids.insert(p.id).second, "duplicate packet/parent ID");
        s->parents.emplace(p.id, p);
    }
    s->chain_digest = sha256(catalog_.digest() + genesis_json(genesis_).dump());
    snapshot_ = std::move(s);
}
PreparedReaction Transactions::prepare(const ForcedEvent &e, ManufacturedFault fault) const {
    identifier(e.id);
    const auto encoded = event_json(e);
    const auto digest = sha256(encoded.dump());
    auto replay = snapshot_->accepted_events.find(e.id);
    if (replay != snapshot_->accepted_events.end())
        throw std::invalid_argument(replay->second == digest ? "duplicate accepted event"
                                                             : "accepted event ID payload conflict");
    require(e.catalog_digest == catalog_.digest(), "catalog provenance mismatch");
    require(e.sequence == snapshot_->next_event_sequence &&
                e.id == genesis_.run_id + "/event/" + std::to_string(e.sequence),
            "event sequence/ID mismatch");
    require(e.manufactured_partition.rfind("manufactured/", 0) == 0,
            "manufactured partition provenance required");
    finite(e.time_s, "event clock");
    finite(e.weight, "represented events");
    require(e.time_s >= snapshot_->time_s && e.weight > 0, "invalid event clock/count");
    const auto &ch = catalog_.channel(e.channel);
    const auto &parent = snapshot_->parents.at(e.parent_id);
    const auto &primary = snapshot_->packets.at(e.primary_id);
    require(parent.state == ch.target_state && catalog_.species_of(primary.state).id == "e",
            "reactant state mismatch");
    require(parent.remaining_weight >= e.weight, "insufficient parent inventory");
    require(primary.weight == e.weight,
            "whole-primary matched weight required; partial splitting unsupported");
    require(primary.position_m == parent.position_m, "co-located reactants required");
    require(e.products.size() == ch.product_states.size(), "product multiplicity mismatch");
    Moments before, after;
    before.add(catalog_, primary.state, e.weight, primary.proper_velocity_m_s);
    before.add(catalog_, parent.state, e.weight, parent.proper_velocity_m_s);
    std::vector<Packet> inserted;
    std::set<std::string> new_ids;
    for (std::size_t i = 0; i < e.products.size(); ++i) {
        const auto &p = e.products[i];
        identifier(p.id);
        require(p.state == ch.product_states[i] && p.weight == e.weight,
                "product state/represented weight mismatch");
        const bool retained = ch.scatters_primary && i == 0;
        require(retained ? p.id == primary.id : p.id != primary.id,
                "primary replacement/consumption ID mismatch");
        require(new_ids.insert(p.id).second, "duplicate product ID");
        require(retained ||
                    (!snapshot_->ever_used_packet_ids.count(p.id) && !snapshot_->parents.count(p.id)),
                "product ID collision/recycling");
        Packet q{p.id,
                 p.state,
                 retained ? primary.source_id : parent.source_id,
                 retained ? primary.birth_event_id : e.id,
                 e.id,
                 p.weight,
                 parent.position_m,
                 p.proper_velocity_m_s};
        validate_packet(catalog_, q);
        inserted.push_back(q);
        after.add(catalog_, q.state, q.weight, q.proper_velocity_m_s);
    }
    require(before.total() >= after.internal + after.reference,
            "insufficient event energy for declared state products");
    close(before.H, after.H, std::abs(before.H) + std::abs(after.H), "H nuclei ledger");
    close(before.O, after.O, std::abs(before.O) + std::abs(after.O), "O nuclei ledger");
    close(before.Q, after.Q, std::abs(before.Q) + std::abs(after.Q), "charge ledger");
    for (int i = 0; i < 3; ++i)
        close(before.momentum[i], after.momentum[i],
              before.absolute_momentum[i] + after.absolute_momentum[i], "momentum ledger");
    close(before.total(), after.total(), before.absolute_energy + after.absolute_energy,
          "total event energy ledger");
    require(snapshot_->journal.size() < maximum_records &&
                snapshot_->packets.size() - 1 + inserted.size() + snapshot_->parents.size() <=
                    maximum_records,
            "transaction record bound");
    auto next = std::make_shared<Snapshot>(*snapshot_);
    auto &target = next->parents.at(e.parent_id);
    const double parent_before = target.remaining_weight;
    target.remaining_weight -= e.weight;
    require(target.remaining_weight >= 0 && target.remaining_weight < parent_before,
            "parent debit below numerical resolution");
    close(parent_before - target.remaining_weight, e.weight, std::abs(parent_before) + e.weight,
          "parent debit roundoff");
    fault_at(fault, ManufacturedFault::AfterParentDebit);
    next->packets.erase(primary.id);
    fault_at(fault, ManufacturedFault::AfterPrimaryRemoval);
    for (const auto &p : inserted) {
        next->packets.emplace(p.id, p);
        next->ever_used_packet_ids.insert(p.id);
    }
    fault_at(fault, ManufacturedFault::AfterProductInsertion);
    Json certificate = {{"before", before.json()},
                        {"after", after.json()},
                        {"parent_debit_count", 1},
                        {"represented_events", e.weight},
                        {"parent_remaining_before", parent_before},
                        {"parent_remaining_after", target.remaining_weight},
                        {"parent_debit_roundoff", parent_before - target.remaining_weight - e.weight},
                        {"energy_residual_J", after.total() - before.total()},
                        {"charge_delta_C", (after.Q - before.Q) * electron_charge_C},
                        {"authored_external_charge_delta_C", 0.},
                        {"authored_external_energy_delta_J", 0.},
                        {"external_event_work_J", 0.},
                        {"photon_port_J", 0.},
                        {"electron_owner", "PIC"},
                        {"deposition_status", "co_located_delta_only_not_native_deposition"}};
    const auto chain = sha256(next->chain_digest + encoded.dump() + certificate.dump());
    next->journal.push_back({{"event", encoded},
                             {"payload_sha256", digest},
                             {"conservation", certificate},
                             {"previous_chain", next->chain_digest},
                             {"chain_sha256", chain}});
    next->chain_digest = chain;
    next->accepted_events.emplace(e.id, digest);
    ++next->revision;
    ++next->next_event_sequence;
    next->time_s = e.time_s;
    ParticleDelta delta{e.id, primary.id, std::move(inserted), std::move(certificate)};
    fault_at(fault, ManufacturedFault::BeforePublication);
    std::shared_ptr<const Snapshot> published = std::move(next);
    return PreparedReaction(snapshot_, std::move(published), std::move(delta));
}
bool Transactions::can_publish(const PreparedReaction &prepared) const noexcept {
    return prepared.input_ == snapshot_;
}
bool Transactions::try_publish(const PreparedReaction &prepared) noexcept {
    if (!can_publish(prepared))
        return false;
    auto next = prepared.output_;
    snapshot_.swap(next);
    return true;
}
ParticleDelta Transactions::apply(const ForcedEvent &event, ManufacturedFault fault) {
    auto prepared = prepare(event, fault);
    auto delta = prepared.delta(); // All allocations finish before publication.
    require(try_publish(prepared), "stale prepared reaction");
    return delta;
}
Json Transactions::checkpoint_document() const {
    Json body = {{"schema", "manufactured-reaction-checkpoint/v1"},
                 {"catalog_sha256", catalog_.digest()},
                 {"genesis", genesis_json(genesis_)},
                 {"snapshot", snapshot_json(*snapshot_)},
                 {"journal", snapshot_->journal}};
    return {{"payload", body}, {"sha256", sha256(body.dump())}};
}
Transactions Transactions::restore(const Catalog &c, const Json &document,
                                   const std::string &expected_prefix) {
    keys(document, {"payload", "sha256"});
    const auto &body = document.at("payload");
    keys(body, {"schema", "catalog_sha256", "genesis", "snapshot", "journal"});
    require(document.at("sha256") == sha256(body.dump()), "checkpoint digest mismatch");
    require(body.at("schema") == "manufactured-reaction-checkpoint/v1" &&
                body.at("catalog_sha256") == c.digest(),
            "checkpoint schema/catalog mismatch");
    auto g = genesis_from(body.at("genesis"));
    require(g.coupling_prefix == expected_prefix, "coupling restart prefix mismatch");
    Transactions result(c, std::move(g));
    require(body.at("journal").is_array() && body.at("journal").size() <= maximum_records,
            "journal size bound");
    for (const auto &record : body.at("journal")) {
        keys(record, {"event", "payload_sha256", "conservation", "previous_chain", "chain_sha256"});
        result.apply(event_from(record.at("event")));
        require(result.snapshot_->journal.back() == record, "restart journal replay differs");
    }
    require(snapshot_json(*result.snapshot_) == body.at("snapshot"),
            "restart state differs from conservative journal replay");
    return result;
}
void Transactions::save_checkpoint(const std::string &path) const {
    const std::filesystem::path target(path),
        directory = target.parent_path().empty() ? std::filesystem::path(".") : target.parent_path();
    std::string pattern = (directory / (target.filename().string() + ".tmp.XXXXXX")).string();
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back('\0');
    int fd = mkstemp(name.data());
    if (fd < 0)
        throw std::system_error(errno, std::generic_category(), "checkpoint temporary create");
    try {
        const auto bytes = checkpoint_document().dump() + "\n";
        require(bytes.size() <= maximum_checkpoint_bytes, "checkpoint size bound");
        std::size_t written = 0;
        while (written < bytes.size()) {
            auto n = write(fd, bytes.data() + written, bytes.size() - written);
            if (n < 0 && errno == EINTR)
                continue;
            if (n <= 0)
                throw std::system_error(errno, std::generic_category(), "checkpoint write");
            written += static_cast<std::size_t>(n);
        }
        if (fsync(fd))
            throw std::system_error(errno, std::generic_category(), "checkpoint fsync");
        if (::close(fd)) {
            fd = -1;
            throw std::system_error(errno, std::generic_category(), "checkpoint close");
        }
        fd = -1;
        if (rename(name.data(), path.c_str()))
            throw std::system_error(errno, std::generic_category(), "checkpoint rename");
        int dirfd = open(directory.c_str(), O_RDONLY | O_DIRECTORY);
        if (dirfd < 0)
            throw std::system_error(errno, std::generic_category(), "checkpoint directory open");
        int status = fsync(dirfd), saved = errno;
        ::close(dirfd);
        if (status)
            throw std::system_error(saved, std::generic_category(), "checkpoint directory fsync");
    } catch (...) {
        if (fd >= 0)
            ::close(fd);
        unlink(name.data());
        throw;
    }
}
Transactions Transactions::load_checkpoint(const Catalog &c, const std::string &path,
                                           const std::string &prefix) {
    require(std::filesystem::file_size(path) <= maximum_checkpoint_bytes, "checkpoint size bound");
    std::ifstream stream(path);
    require(bool(stream), "checkpoint open");
    auto j = Json::parse(stream);
    stream >> std::ws;
    require(stream.eof(), "checkpoint trailing data");
    return restore(c, j, prefix);
}
} // namespace raikou::pic_reactions
