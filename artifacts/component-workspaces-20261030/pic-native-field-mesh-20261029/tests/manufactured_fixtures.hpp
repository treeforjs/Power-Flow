#pragma once
#include "raikou/mfem_gap_pic.h"
#include "reaction_transactions.hpp"
#include <cmath>
namespace manufactured {
using namespace raikou::pic_reactions;
Json definition() {
    // Deliberately authored test masses/reference/internal partitions. These are
    // NOT water thresholds, measured branching, or a physical collision model.
    const double me = raikou::mfem_backend::pic_me, H = 1.6735575e-27, O = 16 * H, e = electron_charge_C;
    return {
        {"schema", "manufactured-reaction-catalog/v1"},
        {"status", "manufactured_forced_only"},
        {"electron_owner", "PIC"},
        {"species",
         Json::array({{{"id", "e"}, {"H", 0}, {"O", 0}, {"charge_e", -1}, {"mass_kg", me}},
                      {{"id", "H2O"}, {"H", 2}, {"O", 1}, {"charge_e", 0}, {"mass_kg", 2 * H + O}},
                      {{"id", "H2Op"}, {"H", 2}, {"O", 1}, {"charge_e", 1}, {"mass_kg", 2 * H + O - me}},
                      {{"id", "OH"}, {"H", 1}, {"O", 1}, {"charge_e", 0}, {"mass_kg", H + O}},
                      {{"id", "H"}, {"H", 1}, {"O", 0}, {"charge_e", 0}, {"mass_kg", H}},
                      {{"id", "Hm"}, {"H", 1}, {"O", 0}, {"charge_e", -1}, {"mass_kg", H + me}}})},
        {"states",
         Json::array(
             {{{"id", "e.g"}, {"species", "e"}, {"reference_J", 0.}, {"internal_J", 0.}},
              {{"id", "H2O.g"}, {"species", "H2O"}, {"reference_J", 0.}, {"internal_J", 0.}},
              {{"id", "H2Op.test"}, {"species", "H2Op"}, {"reference_J", 3 * e}, {"internal_J", .5 * e}},
              {{"id", "OH.test"}, {"species", "OH"}, {"reference_J", .75 * e}, {"internal_J", .5 * e}},
              {{"id", "H.test"}, {"species", "H"}, {"reference_J", .5 * e}, {"internal_J", 0.}},
              {{"id", "Hm.test"}, {"species", "Hm"}, {"reference_J", .25 * e}, {"internal_J", 0.}}})},
        {"channels", Json::array({{{"id", "manufactured/ionization"},
                                   {"target_state", "H2O.g"},
                                   {"scatters_primary", true},
                                   {"product_states", {"e.g", "e.g", "H2Op.test"}}},
                                  {{"id", "manufactured/dissociation"},
                                   {"target_state", "H2O.g"},
                                   {"scatters_primary", true},
                                   {"product_states", {"e.g", "OH.test", "H.test"}}},
                                  {{"id", "manufactured/attachment"},
                                   {"target_state", "H2O.g"},
                                   {"scatters_primary", false},
                                   {"product_states", {"Hm.test", "OH.test"}}}})}};
}
Genesis genesis(const Catalog &c, double parent_weight = 1., double electron_K_eV = 100.) {
    Genesis g;
    g.run_id = "manufactured-run/v1";
    g.coupling_prefix = "manufactured/accepted-interval/zero-no-native-PIC-state";
    g.parents.push_back({"gas-parent-lot/1",
                         "H2O.g",
                         "manufactured/single-parent-inventory",
                         parent_weight,
                         {.003, .004, 0},
                         {0, 0, 0}});
    const double ratio = electron_K_eV * electron_charge_C /
                         (c.species("e").mass_kg * speed_of_light_m_s * speed_of_light_m_s);
    const double u = speed_of_light_m_s * std::sqrt(ratio * (2 + ratio));
    for (int i = 0; i < 3; ++i) {
        Packet p{"primary/" + std::to_string(i),
                 "e.g",
                 "manufactured/authored-electron-source",
                 "manufactured/genesis",
                 "manufactured/genesis",
                 .25,
                 {.003, .004, 0},
                 {u, 0, 0}};
        g.packets.push_back(p);
        g.authored_charge_C -= .25 * electron_charge_C;
        g.authored_kinetic_J += .25 * kinetic_J(c.species("e"), p.proper_velocity_m_s);
    }
    return g;
}
ForcedEvent forced(const Transactions &engine, const std::string &channel,
                   const std::string &primary_id) {
    const auto state = engine.snapshot();
    const auto &c = engine.catalog();
    const auto &primary = state->packets.at(primary_id);
    const auto &parent = state->parents.begin()->second;
    const auto &ch = c.channel(channel);
    ForcedEvent event;
    event.sequence = state->next_event_sequence;
    event.id = engine.genesis().run_id + "/event/" + std::to_string(event.sequence);
    event.time_s = event.sequence * 1e-12;
    event.weight = primary.weight;
    event.channel = channel;
    event.parent_id = parent.id;
    event.primary_id = primary.id;
    event.catalog_digest = c.digest();
    event.manufactured_partition = "manufactured/arbitrary-transverse-pair-recoil/v1";
    Vec P{};
    for (int k = 0; k < 3; ++k)
        P[k] = c.species_of(primary.state).mass_kg * primary.proper_velocity_m_s[k] +
               c.species_of(parent.state).mass_kg * parent.proper_velocity_m_s[k];
    double total_mass = 0, potential = 0;
    for (const auto &st : ch.product_states) {
        total_mass += c.species_of(st).mass_kg;
        potential += c.state(st).reference_J + c.state(st).internal_J;
    }
    const auto &ps = c.state(primary.state);
    const auto &ns = c.state(parent.state);
    const double available = kinetic_J(c.species_of(primary.state), primary.proper_velocity_m_s) +
                             kinetic_J(c.species_of(parent.state), parent.proper_velocity_m_s) +
                             ps.reference_J + ps.internal_J + ns.reference_J + ns.internal_J - potential;
    if (available < 0)
        throw std::runtime_error("manufactured fixture lacks energy");
    for (std::size_t i = 0; i < ch.product_states.size(); ++i)
        event.products.push_back(
            {ch.scatters_primary && i == 0 ? primary.id : event.id + "/product/" + std::to_string(i),
             ch.product_states[i],
             event.weight,
             {}});
    auto partition = [&](double pair_momentum) {
        double K = 0;
        for (std::size_t i = 0; i < event.products.size(); ++i) {
            auto &p = event.products[i];
            const auto &species = c.species_of(p.state);
            for (int k = 0; k < 3; ++k)
                p.proper_velocity_m_s[k] = P[k] / total_mass;
            if (i < 2)
                p.proper_velocity_m_s[1] += (i == 0 ? 1. : -1.) * pair_momentum / species.mass_kg;
            K += kinetic_J(species, p.proper_velocity_m_s);
        }
        return K;
    };
    double lo = 0, hi = 1e-30;
    while (partition(hi) < available) {
        hi *= 2;
        if (!std::isfinite(hi))
            throw std::runtime_error("fixture bracket");
    }
    for (int i = 0; i < 180; ++i) {
        double mid = .5 * (lo + hi);
        if (partition(mid) < available)
            lo = mid;
        else
            hi = mid;
    }
    partition(.5 * (lo + hi));
    return event;
}

}
