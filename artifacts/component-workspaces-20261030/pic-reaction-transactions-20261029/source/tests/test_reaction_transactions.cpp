#include "reaction_transactions.hpp"
#include "native_particle_view.hpp"
#include "integration_ports.hpp"
#include "raikou/mfem_gap_pic.h" // frozen ABI type only: no GapPic object or advance
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace raikou::pic_reactions;
namespace {
Json checks = Json::array();
void check(bool pass, const std::string &name) {
    checks.push_back({{"name", name}, {"passed", pass}});
    if (!pass)
        throw std::runtime_error(name);
}
template <class Call> void rejected(Call call, const std::string &name) {
    bool failed = false;
    try {
        call();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed, name);
}
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
void reject_event(Transactions &engine, const ForcedEvent &event, const std::string &label,
                  ManufacturedFault fault = ManufacturedFault::None) {
    const auto before = engine.checkpoint_document();
    const auto pointer = engine.snapshot();
    rejected([&] { engine.apply(event, fault); }, label);
    check(engine.snapshot() == pointer && engine.checkpoint_document() == before,
          label + " leaves the complete owner unchanged");
}
void rehash(Json &document) {
    document["sha256"] = sha256(document.at("payload").dump());
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc != 2)
            throw std::runtime_error("expected evidence directory");
        const std::filesystem::path out(argv[1]);
        std::filesystem::create_directories(out);
        Catalog catalog(definition());
        auto g = genesis(catalog);
        Transactions owner(catalog, g);
        Transactions staged(catalog, g);
        const auto stage_input = staged.snapshot();
        const auto stage_event = forced(staged, "manufactured/ionization", "primary/0");
        auto proposal = staged.prepare(stage_event);
        auto competing = staged.prepare(stage_event);
        check(staged.snapshot() == stage_input, "preparation has no externally visible effects");
        AcceptedInterval interval;
        interval.step = 1;
        interval.spatial_cell_id = 42;
        interval.begin_s = 0;
        interval.end_s = 2e-12;
        interval.reaction_input_sha256 = sha256(snapshot_json(*stage_input).dump());
        interval.pic_input_sha256 = sha256("manufactured/PIC-input");
        interval.neutral_input_sha256 = sha256("manufactured/neutral-input");
        interval.em_input_sha256 = sha256("manufactured/EM-input");
        interval.material_input_sha256 = sha256("manufactured/material-input");
        interval.deposition_basis_id = "manufactured/common-point-shape";
        const auto ports = integration_ports(catalog, proposal, interval);
        std::ofstream(out / "typed-integration-example.json")
            << integration_ports_json(ports).dump(2) << '\n';
        check(ports.charged_terms.size() == 4 && ports.neutral_products.empty() &&
                  ports.net_created_charge_C == 0 && ports.neutral_debit.debit_operations == 1 &&
                  ports.neutral_debit.represented_count == .25,
              "typed ionization ports carry one debit and simultaneous charge-balanced terms");
        check(ports.energy.electron_translation_J < 0 && ports.energy.species_reference_J > 0 &&
                  ports.energy.state_internal_J > 0 && ports.energy.additional_wall_heat_J == 0 &&
                  ports.energy.external_event_work_J == 0,
              "typed energy transfers do not duplicate resolved recoil/reference as wall heat");
        auto wrong_interval = interval;
        wrong_interval.reaction_input_sha256 = sha256("different snapshot");
        rejected([&] { integration_ports(catalog, proposal, wrong_interval); },
                 "typed ports reject stale reaction prefix");
        wrong_interval = interval;
        wrong_interval.em_input_sha256.clear();
        rejected([&] { integration_ports(catalog, proposal, wrong_interval); },
                 "typed ports require every coupled owner prefix");
        wrong_interval = interval;
        wrong_interval.begin_s = 1.5e-12;
        rejected([&] { integration_ports(catalog, proposal, wrong_interval); },
                 "typed ports reject events outside accepted interval");
        check(staged.can_publish(proposal) && staged.try_publish(proposal),
              "prepared reaction publishes after partner preflight");
        const auto stage_output = staged.snapshot();
        check(!staged.can_publish(competing) && !staged.try_publish(competing) &&
                  !staged.try_publish(proposal) && staged.snapshot() == stage_output,
              "stale/duplicate prepared publications make no changes");
        Transactions foreign(catalog, g);
        check(!foreign.can_publish(proposal) && !foreign.try_publish(proposal),
              "prepared token cannot publish into a separate owner");
        auto three_dimensional = g;
        three_dimensional.packets[0].proper_velocity_m_s[2] =
            .2 * three_dimensional.packets[0].proper_velocity_m_s[0];
        Transactions z_owner(catalog, three_dimensional);
        auto z_prepared = z_owner.prepare(forced(z_owner, "manufactured/attachment", "primary/0"));
        auto z_interval = interval;
        z_interval.reaction_input_sha256 = sha256(snapshot_json(*z_owner.snapshot()).dump());
        auto z_ports = integration_ports(catalog, z_prepared, z_interval);
        std::ofstream(out / "typed-integration-Jz-example.json")
            << integration_ports_json(z_ports).dump(2) << '\n';
        check(z_ports.charge_velocity_moment_change_C_m_s[2] != 0 &&
                  z_ports.neutral_products.size() == 1 && z_ports.charged_terms.size() == 2,
              "three-component current port includes attachment Jz moment and neutral output");
        std::ofstream(out / "manufactured-catalog.json") << catalog.definition().dump(2) << '\n';
        check(owner.snapshot()->packets.size() == 3 && owner.snapshot()->parents.size() == 1,
              "one gas-parent owner and three PIC electron packets");
        auto source = definition();
        Catalog immutable(source);
        source["states"][2]["reference_J"] = 999.;
        check(immutable.digest() == catalog.digest(),
              "caller catalog mutation cannot change immutable definitions");
        for (const auto &name : {"H", "O", "charge", "fractional", "unknown-state", "duplicate-species",
                                 "physical", "F0"}) {
            auto bad = definition();
            std::string n = name;
            if (n == "H")
                bad["species"][3]["H"] = 2;
            else if (n == "O")
                bad["species"][3]["O"] = 2;
            else if (n == "charge")
                bad["species"][2]["charge_e"] = 0;
            else if (n == "fractional")
                bad["species"][3]["H"] = 1.5;
            else if (n == "unknown-state")
                bad["channels"][0]["product_states"][2] = "missing";
            else if (n == "duplicate-species")
                bad["species"].push_back(bad["species"][0]);
            else if (n == "physical")
                bad["status"] = "physical_H2O";
            else
                bad["electron_owner"] = "Raikou-F0";
            rejected([&] { Catalog invalid(bad); }, "invalid catalog: " + n);
        }
        auto wrong_owner = g;
        wrong_owner.electron_owner = "Raikou-F0";
        rejected([&] { Transactions invalid(catalog, wrong_owner); },
                 "duplicate electron/F0 owner rejected");
        auto event = forced(owner, "manufactured/ionization", "primary/0");
        for (auto fault :
             {ManufacturedFault::AfterParentDebit, ManufacturedFault::AfterPrimaryRemoval,
              ManufacturedFault::AfterProductInsertion, ManufacturedFault::BeforePublication})
            reject_event(owner, event,
                         "injected private-stage rollback/" + std::to_string(static_cast<int>(fault)),
                         fault);
        for (int kind = 0; kind < 13; ++kind) {
            auto bad = event;
            std::string label;
            switch (kind) {
            case 0:
                bad.channel = "missing";
                label = "unknown channel";
                break;
            case 1:
                bad.catalog_digest = "changed";
                label = "catalog digest";
                break;
            case 2:
                bad.sequence = 2;
                label = "sequence";
                break;
            case 3:
                bad.id = "different-event-id";
                label = "event immutable ID";
                break;
            case 4:
                bad.manufactured_partition = "physical/data";
                label = "physical partition promotion";
                break;
            case 5:
                bad.products[2].state = "H.test";
                label = "product state";
                break;
            case 6:
                bad.products[1].weight *= 2;
                label = "product matched weight";
                break;
            case 7:
                bad.products[1].id = bad.products[2].id;
                label = "duplicate product ID";
                break;
            case 8:
                bad.products[1].id = "primary/1";
                label = "existing packet ID collision";
                break;
            case 9:
                bad.products[0].id = "new-primary-id";
                label = "scattered primary replacement";
                break;
            case 10:
                bad.products[0].proper_velocity_m_s[1] *= 1.1;
                label = "momentum/energy overspend";
                break;
            case 11:
                bad.products[0].proper_velocity_m_s[0] = std::numeric_limits<double>::quiet_NaN();
                label = "nonfinite product";
                break;
            default:
                bad.time_s = -1.;
                label = "clock";
            }
            reject_event(owner, bad, label);
        }
        Transactions insufficient(catalog, genesis(catalog, .125));
        reject_event(insufficient, event, "insufficient one-parent inventory");
        Transactions tiny(catalog, genesis(catalog, 1e300));
        reject_event(tiny, event, "sub-resolution parent debit");
        Transactions low(catalog, genesis(catalog, 1., 1.));
        auto unpayable = event;
        for (auto &p : unpayable.products)
            p.proper_velocity_m_s = {};
        reject_event(low, unpayable, "insufficient declared event energy");
        auto partial = event;
        partial.weight = .125;
        for (auto &p : partial.products)
            p.weight = .125;
        reject_event(owner, partial, "unsupported partial-primary splitting");
        auto energy_overspend = event;
        for (int i = 0; i < 2; ++i)
            energy_overspend.products[i].proper_velocity_m_s[1] *= 1.1;
        reject_event(owner, energy_overspend, "momentum-balanced energy overspend");
        auto remote = g;
        remote.parents[0].position_m[2] = .001;
        Transactions separate(catalog, remote);
        reject_event(separate, event, "nonlocal reactant pair");
        const auto old_snapshot = owner.snapshot();
        auto first = owner.apply(event);
        check(old_snapshot->parents.at("gas-parent-lot/1").remaining_weight == 1. &&
                  owner.snapshot()->parents.at("gas-parent-lot/1").remaining_weight == .75,
              "snapshot and one parent debit are atomic");
        check(first.inserted.size() == 3 && owner.snapshot()->packets.size() == 5,
              "ionization replaces primary and inserts electron plus ion");
        check(first.conservation.at("parent_debit_count") == 1 &&
                  first.conservation.at("charge_delta_C") == 0.,
              "ionization charge/source bookkeeping");
        check(owner.snapshot()->authored_charge_C == g.authored_charge_C &&
                  owner.snapshot()->authored_kinetic_J == g.authored_kinetic_J,
              "internal reaction never increments authored external source ledger");
        check(owner.snapshot()->packets.at("primary/0").birth_event_id == "manufactured/genesis" &&
                  owner.snapshot()->packets.at(event.id + "/product/1").birth_event_id == event.id,
              "primary lineage retained and secondary lineage created");
        reject_event(owner, event, "duplicate accepted event");
        auto conflicting = event;
        conflicting.time_s += 1.;
        reject_event(owner, conflicting, "same ID with changed payload");
        event.products[0].proper_velocity_m_s[0] = 0;
        check(owner.snapshot()->journal.front().at("event").at("products")[0].at(
                  "proper_velocity_m_s")[0] != 0,
              "caller event mutation cannot change accepted journal");
        // Save/reload at one accepted transaction, then continue independently.
        const auto checkpoint = out / "manufactured-reactions.restart.json";
        owner.save_checkpoint(checkpoint.string());
        auto restored = Transactions::load_checkpoint(catalog, checkpoint.string(), g.coupling_prefix);
        check(restored.checkpoint_document() == owner.checkpoint_document(),
              "restart exact catalog/state/journal/ledger prefix replay");
        reject_event(restored, forced(Transactions(catalog, g), "manufactured/ionization", "primary/0"),
                     "duplicate after restarted continuation");
        auto second = forced(owner, "manufactured/dissociation", "primary/1");
        auto second_delta = owner.apply(second);
        restored.apply(second);
        check(second_delta.inserted.size() == 3 && owner.snapshot()->packets.count("primary/1"),
              "dissociation retains scattered primary with two neutral products");
        auto third = forced(owner, "manufactured/attachment", "primary/2");
        auto third_delta = owner.apply(third);
        restored.apply(third);
        check(!owner.snapshot()->packets.count("primary/2") && third_delta.inserted.size() == 2,
              "attachment consumes electron and creates anion plus neutral");
        check(owner.snapshot()->parents.at("gas-parent-lot/1").remaining_weight == .25,
              "three reactions debit one parent each exactly once");
        check(restored.checkpoint_document().dump() == owner.checkpoint_document().dump(),
              "uninterrupted and restarted continuation byte-identical");
        auto recycled = forced(owner, "manufactured/ionization", "primary/0");
        recycled.products[1].id = "primary/2";
        reject_event(owner, recycled, "consumed particle ID can never be recycled");
        auto exhausted = Transactions::restore(catalog, owner.checkpoint_document(), g.coupling_prefix);
        exhausted.apply(forced(exhausted, "manufactured/ionization", "primary/0"));
        check(exhausted.snapshot()->parents.at("gas-parent-lot/1").remaining_weight == 0.,
              "last parent slice can be consumed exactly");
        auto beyond = forced(exhausted, "manufactured/dissociation", "primary/1");
        reject_event(exhausted, beyond, "exhausted parent inventory rejects further reaction");
        for (const auto &record : owner.snapshot()->journal) {
            const auto &c = record.at("conservation");
            check(c.at("before").at("H_nuclei") == c.at("after").at("H_nuclei") &&
                      c.at("before").at("O_nuclei") == c.at("after").at("O_nuclei"),
                  "accepted event atom ledgers " + record.at("event").at("id").get<std::string>());
            check(c.at("before").at("charge_e") == c.at("after").at("charge_e") &&
                      c.at("external_event_work_J") == 0. &&
                      c.at("authored_external_energy_delta_J") == 0.,
                  "accepted event charge/internal energy ownership " +
                      record.at("event").at("id").get<std::string>());
        }
        // A shared shape at a single point preserves nodal net charge. This is an
        // algebraic witness, not a claim to have run native Q1 deposition/current.
        for (double shape : {.1, .2, .3, .4})
            for (const auto &record : owner.snapshot()->journal) {
                auto c = record.at("conservation");
                double q0 = c.at("before").at("charge_e").get<double>() * electron_charge_C * shape,
                       q1 = c.at("after").at("charge_e").get<double>() * electron_charge_C * shape;
                check(q0 == q1, "co-located common-shape charge witness");
            }
        std::uint64_t native_id = 100;
        for (const auto &p : third_delta.inserted) {
            const auto &sp = catalog.species_of(p.state);
            if (sp.charge_e) {
                auto native =
                    charged_native_view<raikou::mfem_backend::GapParticle>(catalog, p, native_id++);
                check(native.charge_C == -electron_charge_C && native.weight == p.weight &&
                          native.proper_velocity_m_s == p.proper_velocity_m_s,
                      "actual frozen GapParticle anion value view");
            } else
                rejected(
                    [&] {
                        charged_native_view<raikou::mfem_backend::GapParticle>(catalog, p, native_id);
                    },
                    "neutral product has no charged native projection");
        }
        for (const auto &p : first.inserted) {
            auto native =
                charged_native_view<raikou::mfem_backend::GapParticle>(catalog, p, native_id++);
            check(native.charge_C == catalog.species_of(p.state).charge_e * electron_charge_C &&
                      native.mass_kg == catalog.species_of(p.state).mass_kg,
                  "actual frozen GapParticle ion/electron value view");
        }
        rejected(
            [&] {
                charged_native_view<raikou::mfem_backend::GapParticle>(catalog, first.inserted.front(),
                                                                       0);
            },
            "unreserved native ID rejected");
        auto valid = owner.checkpoint_document();
        for (int kind = 0; kind < 8; ++kind) {
            auto bad = valid;
            std::string label;
            switch (kind) {
            case 0:
                bad["sha256"] = "wrong";
                label = "digest";
                break;
            case 1:
                bad["payload"]["snapshot"]["parents"][0]["remaining_weight"] = .5;
                rehash(bad);
                label = "forged snapshot with rehash";
                break;
            case 2:
                bad["payload"]["journal"][0]["conservation"]["parent_debit_count"] = 2;
                rehash(bad);
                label = "forged conservation receipt";
                break;
            case 3:
                bad["payload"]["journal"].push_back(bad["payload"]["journal"][0]);
                rehash(bad);
                label = "duplicated journal event";
                break;
            case 4:
                std::swap(bad["payload"]["journal"][0], bad["payload"]["journal"][1]);
                rehash(bad);
                label = "reordered journal";
                break;
            case 5:
                bad["payload"]["genesis"]["electron_owner"] = "Raikou-F0";
                rehash(bad);
                label = "restored double electron owner";
                break;
            case 6:
                bad["payload"]["catalog_sha256"] = "different";
                rehash(bad);
                label = "catalog provenance";
                break;
            default:
                bad["payload"]["extra"] = "unrecognized";
                rehash(bad);
                label = "unknown restart fields";
            }
            rejected([&] { Transactions::restore(catalog, bad, g.coupling_prefix); },
                     "restart rejects " + label);
        }
        rejected([&] { Transactions::restore(catalog, valid, "different-PIC-coupling-prefix"); },
                 "restart tied to same accepted coupling prefix");
        auto changed = definition();
        changed["states"][2]["reference_J"] = 4 * electron_charge_C;
        Catalog new_catalog(changed);
        rejected([&] { Transactions::restore(new_catalog, valid, g.coupling_prefix); },
                 "restart rejects state/reference-energy mutation");
        std::ofstream(out / "truncated.restart") << "{\"payload\":";
        rejected(
            [&] {
                Transactions::load_checkpoint(catalog, (out / "truncated.restart").string(),
                                              g.coupling_prefix);
            },
            "truncated restart rejected");
        const auto doc_before = owner.checkpoint_document();
        rejected([&] { owner.save_checkpoint((out / "missing-parent/path.json").string()); },
                 "checkpoint I/O failure");
        check(owner.checkpoint_document() == doc_before, "checkpoint I/O failure does not mutate owner");
        owner.save_checkpoint((out / "manufactured-continuation.restart.json").string());
        Json report = {{"status", "passed"},
                       {"checks", checks},
                       {"check_count", checks.size()},
                       {"catalog_sha256", catalog.digest()},
                       {"event_count", owner.snapshot()->journal.size()},
                       {"journal", owner.snapshot()->journal},
                       {"scope", "manufactured forced conservative reaction transactions only; no "
                                 "physical rates, PIC evolution, native deposition, GPU, or F0 solve"},
                       {"native_type_size_bytes", sizeof(raikou::mfem_backend::GapParticle)},
                       {"ledger_roundoff_factor", ledger_roundoff_factor},
                       {"electron_owner", "PIC"}};
        std::ofstream(out / "test-results.json") << report.dump(2) << '\n';
        std::cout << "Passed " << checks.size()
                  << " manufactured transaction checks;3 accepted events; no simulation advance\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        std::cerr << checks.dump(2) << '\n';
        return 1;
    }
}
