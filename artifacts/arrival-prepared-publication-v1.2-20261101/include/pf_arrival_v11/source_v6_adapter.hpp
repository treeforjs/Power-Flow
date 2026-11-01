// source_v6_adapter.hpp -- binds the arrival owner v1.1 to heat/source owner v6
// (pic-electrode-heat-source-v6-20261030; same entry points as source v5).
// Copied from the arrival v1 source v5 adapter with the namespaces changed.
// Include the sealed Hall axial_surface_path_port.hpp and source v6
// source_owner.hpp before this header.
//
//  * inputs: PreparedPatch expected count (scalar + finite log) and its v3
//    Exposure; ReleaseQuote exact debit + Owner::inventory exact stock.
//  * outputs: NeutralLeaf{weight, Draw} lists for Owner::prepare (the exact sum
//    is re-checked there: "v5 exact newborn leaf sum differs from exact source
//    debit"), and v3::BirthSample for Owner::begin_source (weight, marks, time
//    policy UniformHeldExposure, canonical ID = native ID, staged count-owner
//    receipt). Energy/angle/position sampling stays in the source emission law.
//  * neutral release times are read back from the source v5 Trial packets.
#pragma once
#include "pf_arrival_v11/arrival_owner.hpp"

namespace pf_arrival_v11::source_v6 {
namespace s5 = electrode_heat_v6;
namespace v3 = electrode_heat_v3;
namespace eh = electrode_heat_v2;

inline ElectronExposure electron_input(std::uint32_t stream_id, const s5::PreparedPatch& p, const v3::Exposure& x, double step_end_s) {
    ElectronExposure e;
    e.stream_id = stream_id;
    e.expected_number = p.expected_number();
    e.log_expected_number = p.log_expected_number();
    e.has_log = std::isfinite(p.log_expected_number());
    e.begin_s = x.begin_s;
    e.duration_s = x.duration_s;
    e.end_s = step_end_s;
    return e;
}

inline NeutralRelease neutral_input(const Policy& P, std::uint32_t stream_id, const s5::ReleaseQuote& q, const s5::Owner& o, const eh::Epoch& e) {
    if (P.stream(stream_id).patch_index != q.patch_index)
        throw Refusal(RefusalKind::stream_input, "neutral stream patch", "stream=" + std::to_string(stream_id), "release quote of another source patch");
    NeutralRelease r;
    r.stream_id = stream_id;
    r.debit = q.released;
    r.log_debit = q.log_released;
    r.inventory = o.inventory(q.patch_index);
    r.begin_s = e.t;
    r.end_s = e.t + e.dt;   // = TrialFrame::end_s (begin_s + dt_s)
    return r;
}

// One NeutralLeaf list per source patch, leaves in canonical (index) order.
inline std::vector<std::vector<s5::NeutralLeaf>> neutral_leaves(const Staged& s, std::size_t patches) {
    std::vector<std::vector<s5::NeutralLeaf>> out(patches);
    for (const auto& b : s.births()) {
        if (b.kind != Kind::neutral) continue;
        if (b.patch_index >= patches) throw Refusal(RefusalKind::stream_input, "neutral leaf patch", b.packet_id(), "patch outside the source owner");
        out[b.patch_index].push_back({b.weight, eh::Draw{b.marks[0], b.marks[1], b.marks[2], b.marks[3], b.marks[4], b.marks[5]}});
    }
    return out;
}

// The source v3 BirthSample of one staged electron birth.
inline v3::BirthSample birth_sample(const Staged& s, const Birth& b, const s5::PreparedPatch& p, std::uint64_t accepted_revision) {
    if (b.kind != Kind::electron) throw Refusal(RefusalKind::stream_input, "birth sample", b.packet_id(), "not an electron birth");
    v3::BirthSample x;
    x.physical_weight = b.weight;
    x.normal_uniform = b.marks[0];
    x.transverse_uniform = b.marks[1];
    x.azimuth_uniform = b.marks[2];
    x.time_uniform = b.marks[3];
    x.event_id = b.ordinal + 1;              // positive, unique over the run
    x.canonical_particle_id = b.native_id;   // the native ID of this birth
    x.count_owner_revision = accepted_revision;
    x.count_owner_base_hash = s.owner_base_sha256();
    x.count_owner_token = s.token();
    x.exposure_key = p.exposure_key();
    x.time_policy = v3::BirthTime::UniformHeldExposure;
    x.weight_provenance = v3::WeightProvenance::StagedCanonicalOwner;
    return x;
}

// Release times of one neutral stream's leaves from a source v5 Trial (packets of
// that patch, in leaf order); every packet weight must equal its leaf bit for bit.
inline std::vector<double> neutral_times(const s5::Trial& t, const Staged& s, std::uint32_t stream_id, const std::string& source_patch_id) {
    std::vector<double> leaf, packet, times;
    for (const auto& b : s.births())
        if (b.stream_id == stream_id) leaf.push_back(b.weight);
    for (const auto& p : t.packets())
        if (p.patch == source_patch_id) {
            packet.push_back(p.weight);
            times.push_back(p.time);
        }
    if (leaf.size() != packet.size() || !std::equal(leaf.begin(), leaf.end(), packet.begin(), [](double a, double b) {
            return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
        }))
        throw Refusal(RefusalKind::incomplete, "source neutral packets", "stream=" + std::to_string(stream_id), "packets differ from the staged leaves");
    return times;
}

}  // namespace pf_arrival_v11::source_v6
