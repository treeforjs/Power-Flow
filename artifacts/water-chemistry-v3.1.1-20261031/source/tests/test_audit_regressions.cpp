#include "chem_host.hpp"
#include <iostream>

using namespace pfc;
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const Blob b = load_blob(argv[1]);
    const auto sid = [&](const std::string& name) {
        for (pfc_u32 i = 0; i < b.n_species; ++i) if (b.species_names[i] == name) return i;
        throw std::runtime_error("species not found");
    };
    int failures = 0;
    auto check = [&](bool ok, const char* name) {
        std::cout << (ok ? "PASS " : "FAIL ") << name << '\n';
        failures += !ok;
    };
    const auto electron = sid("e"), water = sid("H2O"), waterion = sid("H2O+"), hydrogen = sid("H");
    StepConfig cfg{};
    cfg.seed = 73; cfg.step = 2; cfg.dt_s = 1e-30; cfg.trunc_sigmas = 7;
    cfg.purpose = PFC_PURPOSE_MCC; cfg.max_candidates = 4096; cfg.max_events = 64;
    CellTargets cell{};
    cell.n_targets = 1; cell.species[0] = waterion; cell.n[0] = 1e20; cell.s[0] = 2e6;
    CarrierIn ci{}; ci.id = 1; ci.species = electron; ci.weight = 1;
    EventRec events[64]{}; CarrierOut out{};
    mcc_carrier(b.view, b.cover, cfg, ci, cell, events, &out);
    double exposure = 0;
    for (pfc_u32 j = 0; j < out.n_cover; ++j) exposure += out.cover_exposure[j];
    check(exposure > 0 && out.status == PFC_ST_OK, "thermal co-moving unresolved exposure is positive");
    // The mean energy is inside the measured DR range, but the target thermal
    // support crosses its upper edge; that region must not be omitted.
    const double m = b.view.species[electron].mass_eV, E = 18.0;
    ci.u[0] = PFC_C_M_S * std::sqrt(E * (E + 2 * m)) / m;
    cell.s[0] = 2e5;
    mcc_carrier(b.view, b.cover, cfg, ci, cell, events, &out);
    exposure = 0;
    for (pfc_u32 j = 0; j < out.n_cover; ++j) exposure += out.cover_exposure[j];
    check(exposure > 0 && out.status == PFC_ST_OK, "thermal threshold crossing records unresolved exposure");

    // Two actual equal-weight products must retain both velocities and positions.
    bool release_ok = true;
    for (pfc_u64 seed = 0; seed < 32; ++seed) {
        std::map<InventoryKey, Reservoir> rs;
        EventRec e{}; e.cell = 1; e.weight = 1; e.n_out = 1;
        e.out_species[0] = hydrogen; e.out_dest[0] = PFC_DEST_RESERVOIR;
        double x0[3] = {1, 2, 3}, x1[3] = {4, 5, 6};
        deposit_reservoir(b, rs, e, x0);
        e.out_u[0][0] = 2e5;
        e.ke_out[0] = ke_from_u(b.view.species[hydrogen].mass_eV, v3(2e5, 0, 0));
        deposit_reservoir(b, rs, e, x1);
        const auto packets = release_reservoirs(b, rs, 1, seed, 3);
        const auto& r = rs.at({1, hydrogen});
        bool same = packets.size() == 2 && r.count.is_zero() && r.kinetic.is_zero();
        for (const auto& p : r.momentum) same = same && p.is_zero();
        int slow = 0, fast = 0;
        for (const auto& p : packets) {
            slow += p.u[0] == 0 && p.x[0] == 1;
            fast += p.u[0] == 2e5 && p.x[0] == 4;
        }
        release_ok = release_ok && same && slow == 1 && fast == 1;
    }
    check(release_ok, "equal-weight parcels release exact velocity/position/count/energy/momentum");

    std::vector<CarrierIn> cs(1); cs[0].id = 1; cs[0].cell = 0;
    StepResult sr; sr.outs.resize(1); sr.outs[0].n_events = 2; sr.events.resize(2);
    auto& elastic = sr.events[0]; elastic.carrier = 1; elastic.cell = 0;
    elastic.weight = 1; elastic.target_species = water; elastic.n_out = 1;
    elastic.out_species[0] = water; elastic.out_dest[0] = PFC_DEST_GAS;
    auto& consume = sr.events[1]; consume.carrier = 1; consume.cell = 0;
    consume.weight = 1; consume.target_species = water; consume.n_out = 1;
    consume.out_species[0] = waterion; consume.out_dest[0] = PFC_DEST_BIRTH;
    std::map<InventoryKey, ExactSum> inventory; inventory[{0, water}].add(1);
    const auto d = resolve_debits(b, cs, sr, 2, inventory);
    check(d.refused == 0 && d.admitted_events[0] == 2 && d.debited.at({0, water}).approx() == 1,
          "GAS return leaves inventory for subsequent consuming event");
#ifdef PFC_NEUTRAL_PARCELS_VERSION
    {
        cell.species[0] = sid("OH+"); cell.s[0] = 0;
        ci.u[0] = ci.u[1] = ci.u[2] = 0;
        mcc_carrier(b.view, b.cover, cfg, ci, cell, events, &out);
        Ledger reach; reach.init(b); ledger_add_carrier(reach, out);
        std::uint64_t supports = 0;
        bool physical_zero = true;
        for (pfc_u32 region = 0; region < b.n_regions; ++region) {
            supports += reach.region_support[region];
            physical_zero = physical_zero && reach.region_exposure[region].is_zero() && reach.region_bound[region].is_zero();
        }
        check(out.status == PFC_ST_OK && supports > 0 && physical_zero,
              "cold co-moving unknown channel retains support without fabricating X or B");
        cfg.dt_s = 0;
        mcc_carrier(b.view, b.cover, cfg, ci, cell, events, &out);
        check(out.n_cover == 0, "zero-time stage has no reachable-support exposure");
        cfg.dt_s = 1e-30;
    }
    bool placeholder_bounds_ok = true;
    for (pfc_u32 g = 0; g < b.n_groups; ++g)
        for (pfc_u32 bin = 0; bin <= b.view.groups[g].grid_n; ++bin) {
            const auto region = b.cover.bin_region[b.cover.bin_off[g] + bin];
            if (region != 0xFFFFFFFFu && b.region_names[region].find("no data") != std::string::npos)
                placeholder_bounds_ok = placeholder_bounds_ok && b.cover.bin_row[b.cover.bin_off[g] + bin] == 0xFFFFFFFFu;
        }
    check(placeholder_bounds_ok, "placeholder-containing coverage bins have no finite bound row");
    {
        std::map<InventoryKey, Reservoir> rs;
        EventRec e{}; e.cell = 3; e.n_out = 1; e.out_species[0] = hydrogen; e.out_dest[0] = PFC_DEST_RESERVOIR;
        auto deposit = [&](double weight, double ux, double position) {
            e.weight = weight; e.out_u[0][0] = ux;
            e.ke_out[0] = ke_from_u(b.view.species[hydrogen].mass_eV, v3(ux, 0, 0));
            double x[3] = {position, 0, 0}; deposit_reservoir(b, rs, e, x);
        };
        deposit(0.3, 2e5, 1);
        const auto small = release_reservoirs(b, rs, 1, 0, 1);
        check(small.size() == 1 && small[0].weight == 0.3 && small[0].u[0] == 2e5,
              "sub-quantum product is released without delaying its physical source");
        deposit(2.3, -3e5, 2); deposit(0.5, 7e5, 3);
        const auto before = rs;
        bool refused = false;
        try { (void)release_reservoirs(b, rs, 1, 0, 2, 1); }
        catch (const Error&) { refused = true; }
        const auto& after_cap = rs.at({3, hydrogen});
        const auto& before_cap = before.at({3, hydrogen});
        bool untouched = after_cap.count == before_cap.count && after_cap.momentum == before_cap.momentum &&
                         after_cap.kinetic == before_cap.kinetic && after_cap.parcels.size() == before_cap.parcels.size();
        for (std::size_t j = 0; untouched && j < after_cap.parcels.size(); ++j) {
            untouched = after_cap.parcels[j].weight == before_cap.parcels[j].weight;
            for (int axis = 0; axis < 3; ++axis)
                untouched = untouched && after_cap.parcels[j].u[axis] == before_cap.parcels[j].u[axis] &&
                            after_cap.parcels[j].x[axis] == before_cap.parcels[j].x[axis];
        }
        check(refused && untouched,
              "release cap refuses atomically without dropping persistent parcels");
        auto restored = rs;
        const auto actual = release_reservoirs(b, rs, 1, 0, 2);
        const auto replay = release_reservoirs(b, restored, 1, 999, 999);
        ExactSum released_count, released_ke;
        std::array<ExactSum, 3> released_p;
        bool same = actual.size() == replay.size();
        for (std::size_t k = 0; k < actual.size(); ++k) {
            const auto& p = actual[k];
            same = same && p.weight == replay[k].weight && p.u[0] == replay[k].u[0] && p.x[0] == replay[k].x[0];
            same = same && p.weight > 0 && p.weight <= 1;
            const V3 u = v3(p.u[0], p.u[1], p.u[2]);
            const V3 mom = p_from_u(b.view.species[hydrogen].mass_eV, u);
            released_count.add(p.weight);
            released_ke.add_product(p.weight, ke_from_u(b.view.species[hydrogen].mass_eV, u));
            released_p[0].add_product(p.weight, mom.x);
            released_p[1].add_product(p.weight, mom.y);
            released_p[2].add_product(p.weight, mom.z);
        }
        const auto& stock = before.at({3, hydrogen});
        check(same && released_count == stock.count && released_ke == stock.kinetic && released_p == stock.momentum,
              "unequal parcels across steps conserve exact phase-space moments and replay deterministically");
        auto corrupt = before;
        corrupt.at({3, hydrogen}).parcels.clear();
        refused = false;
        try { (void)release_reservoirs(b, corrupt, 1, 0, 2); } catch (const Error&) { refused = true; }
        check(refused, "moments-only or legacy signed reservoir checkpoint is refused");
    }
#endif
    std::cout << "failures " << failures << '\n';
    return failures ? 1 : 0;
}
