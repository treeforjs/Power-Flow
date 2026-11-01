#pragma once
// Coordinator v8: the Hall field owner type is a template parameter of the step, selected per run
// (--hall-owner owner|resident). The step reaches the owner only through the protocol both types
// implement (prepare_closed_proposal, preflight_final, publish_prevalidated, save, the accessors
// time_s/revision/field_sha256/heat_sha256/kinetic_ledgers/thermal/model/mesh/motion_context) and through
// HallTraits for the places where the two differ:
//   accepted(h, cache)   accepted E, B (edge/face high words) and charge ledger on the host: the Owner's
//                        own arrays; the ResidentOwner's slot-0 download (accepted_state), once per accepted
//                        state (cache keyed by revision and field digest). Diagnostic and source-field input
//                        only: the same words as the Owner's arrays (Hall v3.1-v3.3 parity tests).
//   mid_E/mid_B/mid_seal  the authoritative fixed-point midpoint of a proposal and its role seal;
//   make / restore       bootstrap from an accepted host state, checkpoint restore (each owner's own format).
// Hall v3.4 implements the device side of the ResidentOwner behind the same protocol.
#include "identity_v8.hpp"
#include <cstdint>
#include <filesystem>
#include <limits>

namespace raikou::native_water::v8 {
// Device context of a resident owner (unused by the host Owner).
struct HallContext {
  std::shared_ptr<hc::gpu::Device> device;
  hc::gpu::ResidentOptions resident{};
};
struct AcceptedHost {
  const axial::State* state = nullptr;
  const hc::Vec* E = nullptr;
  const hc::Vec* B = nullptr;
  const hc::ChargeLedger* charge = nullptr;
};

template <class HallT> struct HallTraits;

template <> struct HallTraits<hc::Owner> {
  using Proposal = hc::JointProposal;
  static constexpr char name[] = "owner";
  static constexpr bool resident = false;
  struct Cache {};
  static AcceptedHost accepted(const hc::Owner& h, Cache&) { return {&h.field(), &h.field().E, &h.field().B, &h.charge()}; }
  static const hc::Vec& mid_E(const Proposal& p) { return p.midpoint().E(); }
  static const hc::Vec& mid_B(const Proposal& p) { return p.midpoint().B(); }
  static std::string mid_seal(const Proposal& p) { return p.midpoint().seal(); }
  static std::unique_ptr<hc::Owner> make(std::shared_ptr<const hc::Model> model, axial::State seed, hc::ChargeLedger q, const m::GapPicOptions& o,
                                         std::span<const m::ArenaRateEntry> census, const HallContext&) {
    return std::make_unique<hc::Owner>(std::move(model), std::move(seed), std::move(q), o, census);
  }
  static std::unique_ptr<hc::Owner> restore(std::shared_ptr<const hc::Model> model, const m::GapPicOptions& o,
                                            std::span<const m::ArenaRateEntry> census, const std::filesystem::path& dir, const HallContext&) {
    return hc::Owner::restore(std::move(model), o, census, dir);
  }
};

template <> struct HallTraits<hc::ResidentOwner> {
  using Proposal = hc::ResidentJointProposal;
  static constexpr char name[] = "resident";
  static constexpr bool resident = true;
  struct Cache {
    std::uint64_t revision = std::numeric_limits<std::uint64_t>::max();
    std::string field;
    axial::State state;
    hc::ChargeLedger charge;
  };
  static AcceptedHost accepted(const hc::ResidentOwner& h, Cache& c) {
    if (c.revision != h.revision() || c.field != h.field_sha256()) {
      h.accepted_state(c.state, c.charge);
      c.revision = h.revision();
      c.field = h.field_sha256();
    }
    return {&c.state, &c.state.E, &c.state.B, &c.charge};
  }
  static const hc::Vec& mid_E(const Proposal& p) { return p.midpoint_E(); }
  static const hc::Vec& mid_B(const Proposal& p) { return p.midpoint_B(); }
  static std::string mid_seal(const Proposal& p) { return p.field_heat().midpoint_seal(); }
  static std::unique_ptr<hc::ResidentOwner> make(std::shared_ptr<const hc::Model> model, axial::State seed, hc::ChargeLedger q, const m::GapPicOptions& o,
                                                 std::span<const m::ArenaRateEntry> census, const HallContext& c) {
    if (!c.device) throw std::runtime_error("coordinator v8: the resident Hall owner needs a GPU device");
    return std::make_unique<hc::ResidentOwner>(std::move(model), c.device, std::move(seed), std::move(q), o, census, c.resident);
  }
  static std::unique_ptr<hc::ResidentOwner> restore(std::shared_ptr<const hc::Model> model, const m::GapPicOptions& o,
                                                    std::span<const m::ArenaRateEntry> census, const std::filesystem::path& dir, const HallContext& c) {
    if (!c.device) throw std::runtime_error("coordinator v8: the resident Hall owner needs a GPU device");
    return hc::ResidentOwner::restore(std::move(model), c.device, o, census, dir, c.resident);
  }
};

// Source owner v7 bound to an accepted Hall slice: the sealed hall_source_binding_v7::bind body (owner
// geometry, material, patches, identities, settings, 300 K, source hash, the Hall drive start and the
// Hall slice thermal view) for either Hall owner type; for the host Owner the sealed function itself.
template <class HallT>
std::unique_ptr<s7::Owner> bind_source(const HallT& hall, const hb::SliceSource& s, const s7::OwnerOptions& options) {
  if constexpr (std::is_same_v<HallT, hc::Owner>) {
    return hb::bind(hall, s, options);
  } else {
    const auto& g = hall.mesh();
    return std::make_unique<s7::Owner>(hb::slice_geometry(g, s.body, s.slice), hb::slice_material(), s.patches, s.identities, s.settings, 300.,
                                       s.source_hash, hall.model().drive().start(), hall.thermal().slice(s.body, s.slice), options);
  }
}
}  // namespace raikou::native_water::v8
