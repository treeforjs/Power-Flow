#pragma once
#include "reaction_transactions.hpp"
#include <stdexcept>

namespace raikou::pic_reactions {
// Value conversion only: this does not insert into GapPic or update authored
// source ledgers. A future native accepted-interval owner must reserve uint64
// IDs, sidecars and all vector/deposition allocations before one publication.
template <class NativeParticle>
NativeParticle charged_native_view(const Catalog &catalog, const Packet &packet,
                                   std::uint64_t reserved_native_id) {
    const auto &species = catalog.species_of(packet.state);
    if (!reserved_native_id || !species.charge_e)
        throw std::invalid_argument("charged native view requires reserved ID and nonzero charge");
    NativeParticle result;
    result.id = reserved_native_id;
    result.position_m = packet.position_m;
    result.proper_velocity_m_s = packet.proper_velocity_m_s;
    result.charge_C = species.charge_e * electron_charge_C;
    result.mass_kg = species.mass_kg;
    result.weight = packet.weight;
    return result;
}
} // namespace raikou::pic_reactions
