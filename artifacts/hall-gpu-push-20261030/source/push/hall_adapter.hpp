// hall_adapter.hpp -- builds the push engine's plain tables from the sealed Hall
// objects. Include after the Hall headers (root_api.hpp or motion_receipt.hpp plus
// axial_mesh.hpp) in a translation unit compiled with the Hall flags.
#pragma once
#include "hall_gpu_push.hpp"

#include <algorithm>
#include <map>

namespace hall_gpu::push::hall_adapter {

inline MeshTables mesh_tables(const axial::Mesh& g) {
    MeshTables t;
    t.mesh = {g.nx, g.ny, g.nz, g.nex, g.ney, g.nbx, g.nby, g.nn, g.ne, g.nb, g.side, g.depth, g.dx, g.dy, g.dz, g.zmin};
    t.VN = g.VN;
    t.edge_a.resize(std::size_t(g.ne));
    t.edge_b.resize(std::size_t(g.ne));
    t.edge_length.resize(std::size_t(g.ne));
    for (int e = 0; e < g.ne; ++e) {
        t.edge_a[std::size_t(e)] = g.edges[std::size_t(e)].a;
        t.edge_b[std::size_t(e)] = g.edges[std::size_t(e)].b;
        t.edge_length[std::size_t(e)] = g.edges[std::size_t(e)].length;
    }
    return t;
}

// Same policy and translation frame as prepare_canonical_motion.
inline Geometry geometry(const raikou::mfem_backend::GapPicOptions& o) {
    const auto k = hall_driver::declared_kinetic_policy(o);
    Geometry g;
    g.policy = {k.width_m, k.gap_m, k.short_z_m, k.feed_z_m,
                {k.lower.left_m, k.lower.right_m, k.lower.bottom_m, k.lower.top_m},
                {k.upper.left_m, k.upper.right_m, k.upper.bottom_m, k.upper.top_m}};
    g.native_to_field[0] = -o.width_m / 2;
    g.native_to_field[1] = -o.gap_m / 2;
    g.native_to_field[2] = 0;
    return g;
}

// Census admission of prepare_canonical_motion done on the host: ID order and
// identity (the reference `ordered`), optional sealed 100000-carrier cap, and the
// Count2 packet consistency check. Returns SoA carriers in ID order.
struct OrderedCensus {
    std::vector<std::uint64_t> id;
    std::vector<double> x, y, z, ux, uy, uz, charge, mass, weight;
    long packet_failure_rank = -1;   // first rank whose Count2 packet differs, if any
    CarrierView view() const {
        return {id, x, y, z, ux, uy, uz, charge, mass, weight};
    }
};

inline OrderedCensus order_census(std::span<const raikou::mfem_backend::ArenaRateEntry> census, bool sealed_cap) {
    if (sealed_cap) ::hall::require(census.size() <= 100000, "motion census exceeds sealed native carrier cap");
    std::vector<const raikou::mfem_backend::ArenaRateEntry*> rows;
    rows.reserve(census.size());
    for (const auto& r : census) rows.push_back(&r);
    std::sort(rows.begin(), rows.end(), [](auto a, auto b) { return a->native.id < b->native.id; });
    std::uint64_t prior = 0;
    for (auto row : rows) {
        ::hall::require(row->native.id > prior && !row->packet_id.empty() && !row->state.empty(), "motion missing/duplicate native/packet identity");
        prior = row->native.id;
    }
    OrderedCensus c;
    const std::size_t n = rows.size();
    for (auto* v : {&c.x, &c.y, &c.z, &c.ux, &c.uy, &c.uz, &c.charge, &c.mass, &c.weight}) v->resize(n);
    c.id.resize(n);
    std::map<std::string, std::size_t> packets;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = rows[i]->native;
        c.id[i] = p.id;
        c.x[i] = p.position_m[0]; c.y[i] = p.position_m[1]; c.z[i] = p.position_m[2];
        c.ux[i] = p.proper_velocity_m_s[0]; c.uy[i] = p.proper_velocity_m_s[1]; c.uz[i] = p.proper_velocity_m_s[2];
        c.charge[i] = p.charge_C; c.mass[i] = p.mass_kg; c.weight[i] = p.weight;
        auto [where, inserted] = packets.emplace(rows[i]->packet_id, i);
        if (!inserted && c.packet_failure_rank < 0) {
            const auto& q = rows[where->second]->native;
            if (!(rows[i]->state == rows[where->second]->state && p.position_m == q.position_m &&
                  p.proper_velocity_m_s == q.proper_velocity_m_s && p.charge_C == q.charge_C && p.mass_kg == q.mass_kg))
                c.packet_failure_rank = long(i);
        }
    }
    return c;
}

}  // namespace hall_gpu::push::hall_adapter
