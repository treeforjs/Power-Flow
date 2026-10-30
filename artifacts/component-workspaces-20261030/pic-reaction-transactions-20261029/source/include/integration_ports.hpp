#pragma once
#include "reaction_transactions.hpp"
#include <cmath>
#include <stdexcept>

namespace raikou::pic_reactions {
// These are typed proposals, not native field/current/material operations.
struct AcceptedInterval {
    std::uint64_t step = 0, spatial_cell_id = 0;
    double begin_s = 0, end_s = 0;
    std::string reaction_input_sha256;
    std::string pic_input_sha256, neutral_input_sha256, em_input_sha256, material_input_sha256;
    std::string deposition_basis_id;
};
struct ChargedPointTerm {
    std::string packet_id, state_id;
    int sign = 0, charge_e = 0; // -1 removes the incoming packet, +1 inserts output
    double weight = 0, mass_kg = 0;
    Vec position_m{}, proper_velocity_m_s{};
};
struct NeutralParentDebit {
    std::string lot_id, state_id, source_id;
    double represented_count = 0;
    unsigned debit_operations = 1;
    // This is the same reserved debit already present in PreparedReaction,
    // not permission for a second neutral/surface inventory subtraction.
};
struct InternalEnergyTransfer {
    double electron_translation_J = 0, ion_translation_J = 0, neutral_translation_J = 0;
    double state_internal_J = 0, species_reference_J = 0, numerical_residual_J = 0;
    double photon_out_J = 0, external_event_work_J = 0;
    // Resolved recoil and state/reference changes are owned above. They are not
    // also wall heat or an unmodified all-heavy caloric source.
    double additional_wall_heat_J = 0;
};
struct ReactionPorts {
    AcceptedInterval interval;
    std::string event_id, channel_id, payload_sha256, catalog_sha256;
    NeutralParentDebit neutral_debit;
    std::vector<ChargedPointTerm> charged_terms;
    std::vector<Packet> neutral_products;
    InternalEnergyTransfer energy;
    double net_created_charge_C = 0;
    // A point-current moment, C*m/s. This is NOT deposited J[A/m^2] or a
    // path-current operator. The EM owner must provide all3 compatible currents.
    Vec charge_velocity_moment_change_C_m_s{};
};

inline ReactionPorts integration_ports(const Catalog &catalog, const PreparedReaction &prepared,
                                       AcceptedInterval interval) {
    auto require = [](bool condition, const char *message) {
        if (!condition)
            throw std::invalid_argument(message);
    };
    auto digest = [](const std::string &s) {
        return s.size() == 64 && s.find_first_not_of("0123456789abcdef") == std::string::npos;
    };
    require(interval.step && interval.spatial_cell_id && std::isfinite(interval.begin_s) &&
                std::isfinite(interval.end_s) && interval.begin_s >= 0 &&
                interval.end_s > interval.begin_s,
            "accepted interval/cell identity required");
    for (const auto *hash :
         {&interval.reaction_input_sha256, &interval.pic_input_sha256, &interval.neutral_input_sha256,
          &interval.em_input_sha256, &interval.material_input_sha256})
        require(digest(*hash), "all owner input-prefix SHA256 bindings required");
    require(!interval.deposition_basis_id.empty(), "explicit native deposition basis required");
    require(interval.reaction_input_sha256 == sha256(snapshot_json(*prepared.input()).dump()),
            "reaction input-prefix mismatch");
    const auto &record = prepared.output()->journal.back();
    const auto &event = record.at("event");
    double time = event.at("time_s");
    require(time >= interval.begin_s && time <= interval.end_s, "event outside accepted interval");
    require(event.at("catalog_digest") == catalog.digest(), "port catalog mismatch");
    const auto &parent = prepared.input()->parents.at(event.at("parent_id").get<std::string>());
    const auto &primary = prepared.input()->packets.at(prepared.delta().removed_primary);
    ReactionPorts ports;
    ports.interval = std::move(interval);
    ports.event_id = prepared.delta().event_id;
    ports.channel_id = event.at("channel");
    ports.payload_sha256 = record.at("payload_sha256");
    ports.catalog_sha256 = catalog.digest();
    ports.neutral_debit = {parent.id, parent.state, parent.source_id, event.at("weight"), 1};
    auto add_charged = [&](const Packet &packet, int sign) {
        const auto &species = catalog.species_of(packet.state);
        ports.charged_terms.push_back({packet.id, packet.state, sign, species.charge_e, packet.weight,
                                       species.mass_kg, packet.position_m, packet.proper_velocity_m_s});
        const double charge = sign * packet.weight * species.charge_e * electron_charge_C;
        ports.net_created_charge_C += charge;
        double u2 = 0;
        for (double component : packet.proper_velocity_m_s)
            u2 += component * component;
        double gamma = std::sqrt(1 + u2 / (speed_of_light_m_s * speed_of_light_m_s));
        for (int k = 0; k < 3; ++k)
            ports.charge_velocity_moment_change_C_m_s[k] +=
                charge * packet.proper_velocity_m_s[k] / gamma;
    };
    add_charged(primary, -1);
    for (const auto &packet : prepared.delta().inserted) {
        if (catalog.species_of(packet.state).charge_e)
            add_charged(packet, +1);
        else
            ports.neutral_products.push_back(packet);
    }
    const auto &before = prepared.delta().conservation.at("before");
    const auto &after = prepared.delta().conservation.at("after");
    auto change = [&](const char *key) {
        return after.at(key).get<double>() - before.at(key).get<double>();
    };
    ports.energy.electron_translation_J = change("electron_kinetic_J");
    ports.energy.ion_translation_J = change("ion_kinetic_J");
    ports.energy.neutral_translation_J = change("neutral_kinetic_J");
    ports.energy.state_internal_J = change("internal_J");
    ports.energy.species_reference_J = change("reference_J");
    ports.energy.numerical_residual_J = prepared.delta().conservation.at("energy_residual_J");
    return ports;
}
inline Json integration_ports_json(const ReactionPorts &p) {
    Json charged = Json::array(), neutrals = Json::array();
    for (const auto &term : p.charged_terms)
        charged.push_back({{"packet_id", term.packet_id},
                           {"state_id", term.state_id},
                           {"sign", term.sign},
                           {"charge_e", term.charge_e},
                           {"weight", term.weight},
                           {"mass_kg", term.mass_kg},
                           {"position_m", term.position_m},
                           {"proper_velocity_m_s", term.proper_velocity_m_s}});
    for (const auto &packet : p.neutral_products)
        neutrals.push_back({{"packet_id", packet.id},
                            {"state_id", packet.state},
                            {"source_id", packet.source_id},
                            {"birth_event_id", packet.birth_event_id},
                            {"weight", packet.weight},
                            {"position_m", packet.position_m},
                            {"proper_velocity_m_s", packet.proper_velocity_m_s}});
    return {{"schema", "manufactured-reaction-integration-ports/v1"},
            {"status", "staged_only_no_native_deposition_or_partner_commit"},
            {"event_id", p.event_id},
            {"channel_id", p.channel_id},
            {"payload_sha256", p.payload_sha256},
            {"catalog_sha256", p.catalog_sha256},
            {"accepted_interval",
             {{"step", p.interval.step},
              {"cell_id", p.interval.spatial_cell_id},
              {"begin_s", p.interval.begin_s},
              {"end_s", p.interval.end_s},
              {"reaction_input_sha256", p.interval.reaction_input_sha256},
              {"pic_input_sha256", p.interval.pic_input_sha256},
              {"neutral_input_sha256", p.interval.neutral_input_sha256},
              {"em_input_sha256", p.interval.em_input_sha256},
              {"material_input_sha256", p.interval.material_input_sha256},
              {"deposition_basis_id", p.interval.deposition_basis_id}}},
            {"neutral_debit",
             {{"lot_id", p.neutral_debit.lot_id},
              {"state_id", p.neutral_debit.state_id},
              {"source_id", p.neutral_debit.source_id},
              {"represented_count", p.neutral_debit.represented_count},
              {"debit_operations", p.neutral_debit.debit_operations},
              {"semantics", "same_reserved_debit_not_an_additional_subtraction"}}},
            {"charged_point_terms", charged},
            {"neutral_products", neutrals},
            {"net_created_charge_C", p.net_created_charge_C},
            {"raw_current_moment_change_C_m_s", p.charge_velocity_moment_change_C_m_s},
            {"energy_transfer_J",
             {{"electron_translation", p.energy.electron_translation_J},
              {"ion_translation", p.energy.ion_translation_J},
              {"neutral_translation", p.energy.neutral_translation_J},
              {"internal", p.energy.state_internal_J},
              {"reference", p.energy.species_reference_J},
              {"numerical_residual", p.energy.numerical_residual_J},
              {"additional_wall_heat", p.energy.additional_wall_heat_J},
              {"external_event_work", p.energy.external_event_work_J},
              {"photon_out", p.energy.photon_out_J}}}};
}
} // namespace raikou::pic_reactions
