#include "raikou/mfem_pic_reaction_batch.h"
#include "native_particle_view.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <stdexcept>
#include <fstream>
#include <filesystem>
#include <fcntl.h>
#include <unistd.h>
#include <system_error>
#include <cerrno>
#include <iterator>

#ifndef PIC_REACTION_NATIVE_SOURCE_SHA256
#error "native source identity must be bound by the frozen build manifest"
#endif

namespace raikou::mfem_backend {
namespace r = raikou::pic_reactions;
namespace {
void require(bool condition, const char *message) {
    if (!condition)
        throw std::invalid_argument(message);
}
bool particle_equal(const GapParticle &a, const GapParticle &b) noexcept {
    return a.id == b.id && a.position_m == b.position_m &&
           a.proper_velocity_m_s == b.proper_velocity_m_s && a.charge_C == b.charge_C &&
           a.mass_kg == b.mass_kg && a.weight == b.weight;
}
void native_particle_guard(const GapParticle &p) {
    double u2 = 0;
    for (double x : p.position_m)
        require(std::isfinite(x), "nonfinite native product position");
    for (double u : p.proper_velocity_m_s) {
        require(std::isfinite(u), "nonfinite native product velocity");
        u2 += u * u;
    }
    require(std::isfinite(u2) && std::isfinite(p.charge_C) && std::isfinite(p.mass_kg) &&
                p.mass_kg > 0 && std::isfinite(p.weight) && p.weight > 0 &&
                std::isfinite(p.charge_C * p.weight) && std::isfinite(p.charge_C / p.mass_kg) &&
                std::isfinite(pic_kinetic_energy_J(p) * p.weight),
            "native product guard");
}
bool options_equal(const GapPicOptions &a, const GapPicOptions &b) noexcept {
    const auto foil = [](const GapFoil &x, const GapFoil &y) {
        return x.left_m == y.left_m && x.right_m == y.right_m && x.bottom_m == y.bottom_m &&
               x.top_m == y.top_m;
    };
    return a.nx == b.nx && a.ny == b.ny && a.width_m == b.width_m && a.gap_m == b.gap_m &&
           a.depth_m == b.depth_m && a.lower_potential_V == b.lower_potential_V &&
           a.upper_potential_V == b.upper_potential_V && a.magnetic_field_T == b.magnetic_field_T &&
           a.distributed_fields == b.distributed_fields && a.finite_foils == b.finite_foils &&
           a.finite_axial_extent == b.finite_axial_extent && a.short_z_m == b.short_z_m &&
           a.feed_z_m == b.feed_z_m && foil(a.lower_foil, b.lower_foil) &&
           foil(a.upper_foil, b.upper_foil);
}
void fail_at(NativeReactionFault requested, NativeReactionFault point) {
    if (requested == point)
        throw std::runtime_error("manufactured native staging fault");
}
struct TemporaryNativeFile {
    std::string path;
    explicit TemporaryNativeFile(const std::string &near_path) {
        std::string pattern = near_path + ".native.XXXXXX";
        std::vector<char> name(pattern.begin(), pattern.end());
        name.push_back('\0');
        int fd = mkstemp(name.data());
        if (fd < 0)
            throw std::system_error(errno, std::generic_category(), "native temporary create");
        ::close(fd);
        path = name.data();
    }
    ~TemporaryNativeFile() {
        unlink(path.c_str());
        unlink((path + ".tmp").c_str());
    }
};
std::string hex_bytes(const std::string &bytes) {
    static constexpr char digits[] = "0123456789abcdef";
    std::string out;
    out.reserve(2 * bytes.size());
    for (unsigned char x : bytes) {
        out.push_back(digits[x >> 4]);
        out.push_back(digits[x & 15]);
    }
    return out;
}
std::string unhex(const std::string &text) {
    require(text.size() % 2 == 0 && text.size() <= 32 * 1024 * 1024, "native bundle byte bound");
    std::string bytes;
    bytes.reserve(text.size() / 2);
    auto digit = [](char c) -> unsigned {
        if (c >= '0' && c <= '9')
            return c - '0';
        if (c >= 'a' && c <= 'f')
            return c - 'a' + 10;
        throw std::invalid_argument("native hex encoding");
    };
    for (std::size_t i = 0; i < text.size(); i += 2)
        bytes.push_back(char(16 * digit(text[i]) + digit(text[i + 1])));
    return bytes;
}
void atomic_json(const r::Json &document, const std::string &path) {
    const auto bytes = document.dump() + "\n";
    require(bytes.size() <= 64 * 1024 * 1024, "native bundle size bound");
    std::string pattern = path + ".tmp.XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back('\0');
    int fd = mkstemp(name.data());
    if (fd < 0)
        throw std::system_error(errno, std::generic_category(), "bundle temporary create");
    try {
        std::size_t offset = 0;
        while (offset < bytes.size()) {
            auto n = write(fd, bytes.data() + offset, bytes.size() - offset);
            if (n < 0 && errno == EINTR)
                continue;
            if (n <= 0)
                throw std::system_error(errno, std::generic_category(), "bundle write");
            offset += std::size_t(n);
        }
        if (fsync(fd))
            throw std::system_error(errno, std::generic_category(), "bundle fsync");
        ::close(fd);
        fd = -1;
        if (rename(name.data(), path.c_str()))
            throw std::system_error(errno, std::generic_category(), "bundle rename");
        auto directory = std::filesystem::path(path).parent_path();
        if (directory.empty())
            directory = ".";
        int d = open(directory.c_str(), O_RDONLY | O_DIRECTORY);
        if (d < 0)
            throw std::system_error(errno, std::generic_category(), "bundle directory open");
        int status = fsync(d), saved = errno;
        ::close(d);
        if (status)
            throw std::system_error(saved, std::generic_category(), "bundle directory fsync");
    } catch (...) {
        if (fd >= 0)
            ::close(fd);
        unlink(name.data());
        throw;
    }
}
void bound_population(const GapPic &p, const r::Transactions &t,
                      const ReactionPacketBindings &bindings) {
    std::map<std::uint64_t, const r::Packet *> packets;
    for (const auto &[id, packet] : t.snapshot()->packets)
        if (t.catalog().species_of(packet.state).charge_e) {
            auto it = bindings.find(id);
            require(it != bindings.end() && packets.emplace(it->second, &packet).second,
                    "bundle particle binding");
        }
    require(bindings.size() == packets.size() && p.particles().size() == packets.size(),
            "bundle population count mismatch");
    for (const auto &native : p.particles()) {
        auto it = packets.find(native.id);
        require(it != packets.end(), "bundle missing particle binding");
        require(particle_equal(native,
                               r::charged_native_view<GapParticle>(t.catalog(), *it->second, native.id)),
                "bundle native/core particle mismatch");
    }
    require(p.diagnostics().time_s == t.snapshot()->time_s, "bundle native/core clock mismatch");
    const auto d = p.diagnostics();
    require(d.authored_charge_C == t.genesis().authored_charge_C &&
                d.authored_kinetic_energy_J == t.genesis().authored_kinetic_J,
            "bundle authored-source prefix mismatch");
}
} // namespace
struct PicReactionBatch::NativeState {
    GapPicOptions options;
    std::vector<GapParticle> particles;
    std::vector<double> horizontal, vertical;
    mfem::Vector potential;
    GapPicDiagnostics ledger;
    std::uint64_t next_id = 0;
    bool fields_current = false;
};
struct PreparedPicReaction::Impl {
    GapPic *pic_owner;
    r::Transactions *reaction_owner;
    ReactionPacketBindings *binding_owner;
    r::PreparedReaction reaction;
    PicReactionBatch::NativeState before, after, work;
    ReactionPacketBindings bindings_before, bindings_after, bindings_work;
    r::Json audit;
    enum class Status { Staged, Published, RolledBack };
    Status status = Status::Staged;
    Impl(GapPic &p, r::Transactions &t, ReactionPacketBindings &b, const r::PreparedReaction &e)
        : pic_owner(&p), reaction_owner(&t), binding_owner(&b), reaction(e) {}
};
PreparedPicReaction::PreparedPicReaction(std::unique_ptr<Impl> p) : impl_(std::move(p)) {}
PreparedPicReaction::PreparedPicReaction(PreparedPicReaction &&) noexcept = default;
PreparedPicReaction &PreparedPicReaction::operator=(PreparedPicReaction &&) noexcept = default;
PreparedPicReaction::~PreparedPicReaction() = default;
const ReactionPacketBindings &PreparedPicReaction::output_bindings() const {
    return impl_->bindings_after;
}
const r::Json &PreparedPicReaction::audit() const {
    return impl_->audit;
}
PicReactionBatch::NativeState PicReactionBatch::capture(const GapPic &p) {
    NativeState s;
    s.options = p.options_;
    s.particles = p.particles_;
    s.horizontal = p.horizontal_transfers_;
    s.vertical = p.vertical_transfers_;
    s.potential = p.potential_;
    s.ledger = p.ledger_;
    s.next_id = p.next_particle_id_;
    s.fields_current = p.fields_current_;
    return s;
}
void PicReactionBatch::swap_state(GapPic &p, NativeState &s) noexcept {
    p.particles_.swap(s.particles);
    p.horizontal_transfers_.swap(s.horizontal);
    p.vertical_transfers_.swap(s.vertical);
    p.potential_.Swap(s.potential);
    std::swap(p.ledger_, s.ledger);
    std::swap(p.next_particle_id_, s.next_id);
    std::swap(p.fields_current_, s.fields_current);
}
bool PicReactionBatch::matches(const GapPic &p, const NativeState &s) noexcept {
    static_assert(sizeof(GapPicDiagnostics) == 200 && sizeof(GapParticle) == 80 &&
                  sizeof(GapPicOptions) == 168);
    if (!options_equal(p.options_, s.options) || p.particles_.size() != s.particles.size() ||
        p.potential_.Size() != s.potential.Size() || p.next_particle_id_ != s.next_id ||
        p.fields_current_ != s.fields_current || p.horizontal_transfers_ != s.horizontal ||
        p.vertical_transfers_ != s.vertical || std::memcmp(&p.ledger_, &s.ledger, sizeof(s.ledger)))
        return false;
    for (std::size_t i = 0; i < s.particles.size(); ++i)
        if (!particle_equal(p.particles_[i], s.particles[i]))
            return false;
    for (int i = 0; i < s.potential.Size(); ++i)
        if (p.potential_[i] != s.potential[i])
            return false;
    return true;
}
std::string PicReactionBatch::fingerprint(const GapPic &p) {
    std::string data;
    auto add = [&](const auto &x) { data.append(reinterpret_cast<const char *>(&x), sizeof(x)); };
    add(p.options_.nx);
    add(p.options_.ny);
    for (double x :
         {p.options_.width_m, p.options_.gap_m, p.options_.depth_m, p.options_.lower_potential_V,
          p.options_.upper_potential_V, p.options_.short_z_m, p.options_.feed_z_m})
        add(x);
    for (double x : p.options_.magnetic_field_T)
        add(x);
    for (bool x :
         {p.options_.distributed_fields, p.options_.finite_foils, p.options_.finite_axial_extent})
        add(x);
    for (auto f : {p.options_.lower_foil, p.options_.upper_foil})
        for (double x : {f.left_m, f.right_m, f.bottom_m, f.top_m})
            add(x);
    add(p.ledger_);
    add(p.next_particle_id_);
    add(p.fields_current_);
    std::uint64_t count = p.particles_.size();
    add(count);
    for (const auto &particle : p.particles_)
        add(particle);
    for (int i = 0; i < p.potential_.Size(); ++i)
        add(p.potential_[i]);
    for (double x : p.horizontal_transfers_)
        add(x);
    for (double x : p.vertical_transfers_)
        add(x);
    return r::sha256(data);
}
const char *PicReactionBatch::source_identity() noexcept {
    return PIC_REACTION_NATIVE_SOURCE_SHA256;
}
void PicReactionBatch::save_bundle(const GapPic &p, const r::Transactions &t,
                                   const ReactionPacketBindings &bindings, const std::string &path) {
    require(p.fields_current_ && !p.options_.distributed_fields,
            "bundle requires solved serial native state");
    bound_population(p, t, bindings);
    TemporaryNativeFile temporary(path);
    p.save_checkpoint(temporary.path);
    std::ifstream stream(temporary.path, std::ios::binary);
    std::string bytes((std::istreambuf_iterator<char>(stream)), {});
    r::Json fields = r::Json::array();
    for (int i = 0; i < p.potential_.Size(); ++i)
        fields.push_back(p.potential_[i]);
    r::Json body = {
        {"schema", "manufactured-native-reaction-bundle/v1"},
        {"native_source_sha256", source_identity()},
        {"native_abi", {sizeof(GapParticle), sizeof(GapPicOptions), sizeof(GapPicDiagnostics)}},
        {"native_bytes_hex", hex_bytes(bytes)},
        {"native_bytes_sha256", r::sha256(bytes)},
        {"reaction_checkpoint", t.checkpoint_document()},
        {"packet_bindings", bindings},
        {"saved_potential_V", fields},
        {"scope", "native electrostatic/reaction only; no external EM/material/source closure state"}};
    atomic_json({{"payload", body}, {"sha256", r::sha256(body.dump())}}, path);
}
RestoredPicReaction PicReactionBatch::load_bundle(const r::Catalog &catalog,
                                                  const GapPicOptions &options,
                                                  const std::string &prefix, const std::string &path) {
    require(std::filesystem::file_size(path) <= 64 * 1024 * 1024, "bundle read size bound");
    std::ifstream input(path);
    auto doc = r::Json::parse(input);
    require(doc.size() == 2 && doc.contains("payload") && doc.contains("sha256"),
            "bundle wrapper fields");
    const auto &body = doc.at("payload");
    require(doc.at("sha256") == r::sha256(body.dump()), "native bundle digest mismatch");
    require(body.size() == 9 && body.at("schema") == "manufactured-native-reaction-bundle/v1" &&
                body.at("native_source_sha256") == source_identity() &&
                body.at("native_abi") == r::Json::array({sizeof(GapParticle), sizeof(GapPicOptions),
                                                         sizeof(GapPicDiagnostics)}),
            "native bundle source/schema/ABI mismatch");
    auto bytes = unhex(body.at("native_bytes_hex"));
    require(r::sha256(bytes) == body.at("native_bytes_sha256"), "native byte digest mismatch");
    auto core = r::Transactions::restore(catalog, body.at("reaction_checkpoint"), prefix);
    auto pic = std::make_unique<GapPic>(options);
    TemporaryNativeFile temporary(path);
    {
        std::ofstream output(temporary.path, std::ios::binary);
        output.write(bytes.data(), bytes.size());
        require(bool(output), "native temporary write");
    }
    pic->load_checkpoint(temporary.path);
    auto bindings = body.at("packet_bindings").get<ReactionPacketBindings>();
    bound_population(*pic, core, bindings);
    const auto &saved = body.at("saved_potential_V");
    require(saved.is_array() && saved.size() == std::size_t(pic->potential_.Size()),
            "native saved field size");
    double scale = 0, difference = 0;
    for (int i = 0; i < pic->potential_.Size(); ++i) {
        double value = saved[i];
        require(std::isfinite(value), "nonfinite saved field");
        scale = std::max(scale, std::abs(value));
        difference = std::max(difference, std::abs(value - pic->potential_[i]));
    }
    // Existing saved-state field relative budget, not a changed solver tolerance.
    require(difference <= 1e-10 * std::max(scale, 1e-300),
            "native restarted field differs from saved state budget");
    return {std::move(pic), std::make_unique<r::Transactions>(std::move(core)), std::move(bindings)};
}
PreparedPicReaction
PicReactionBatch::prepare(GapPic &p, r::Transactions &t, const r::PreparedReaction &e,
                          ReactionPacketBindings &binding, const ReactionPhase &phase,
                          const PrimaryWeightSlice &slice, NativeReactionFault fault) {
    require(p.fields_current_ && !p.options_.distributed_fields,
            "native reaction v1 requires solved CPU-serial fields");
    require(t.can_publish(e), "stale reaction preparation");
    require(phase.expected_time_s == p.ledger_.time_s, "native/reaction clock mismatch");
    const auto &event = e.output()->journal.back().at("event");
    require(event.at("time_s") == p.ledger_.time_s, "event is not at the native synchronized state");
    if (phase.kind == ReactionPhase::Kind::InitialSynchronizedFieldProof)
        require(phase.accepted_step == 0 && p.ledger_.time_s == 0 && phase.previous_dt_s == 0 &&
                    phase.next_dt_s == 0,
                "initial field-only phase requires initial synchronized state");
    else
        require(phase.kind == ReactionPhase::Kind::MidXEqualStepBoundary && phase.accepted_step > 0 &&
                    std::isfinite(phase.previous_dt_s) && phase.previous_dt_s > 0 &&
                    phase.previous_dt_s == phase.next_dt_s && p.ledger_.time_s >= phase.previous_dt_s,
                "provisional mid-x phase requires equal positive adjacent steps");
    require(slice.event_weight == event.at("weight").get<double>() &&
                slice.untouched_primary_native_id == 0,
            "v1 rejects partial-primary splitting/remainder IDs");
    const auto &primary = e.input()->packets.at(e.delta().removed_primary);
    require(slice.event_weight == primary.weight, "whole-primary weight required");
    std::map<std::uint64_t, std::string> reverse;
    std::size_t charged = 0;
    for (const auto &[key, packet] : e.input()->packets)
        if (t.catalog().species_of(packet.state).charge_e) {
            ++charged;
            auto it = binding.find(key);
            require(it != binding.end() && it->second > 0 && it->second < p.next_particle_id_,
                    "missing/stale native particle binding");
            require(reverse.emplace(it->second, key).second, "duplicate native particle binding");
        }
    require(binding.size() == charged && p.particles_.size() == charged,
            "native/core charged population mismatch");
    for (const auto &native : p.particles_) {
        auto it = reverse.find(native.id);
        require(it != reverse.end(), "unbound native particle");
        const auto &packet = e.input()->packets.at(it->second);
        auto expected = r::charged_native_view<GapParticle>(t.catalog(), packet, native.id);
        require(particle_equal(native, expected), "native/core particle values differ");
    }
    auto impl = std::make_unique<PreparedPicReaction::Impl>(p, t, binding, e);
    impl->before = capture(p);
    impl->bindings_before = binding;
    impl->bindings_after = binding;
    const auto primary_id = binding.at(primary.id);
    std::vector<GapParticle> replacement;
    replacement.reserve(p.particles_.size() + e.delta().inserted.size());
    auto next = p.next_particle_id_;
    std::size_t new_charged = 0;
    impl->bindings_after.erase(primary.id);
    for (const auto &old : p.particles_) {
        if (old.id != primary_id) {
            replacement.push_back(old);
            continue;
        }
        for (const auto &packet : e.delta().inserted)
            if (t.catalog().species_of(packet.state).charge_e) {
                std::uint64_t id;
                if (packet.id == primary.id)
                    id = old.id;
                else {
                    require(next < std::numeric_limits<std::uint64_t>::max(),
                            "native ID range exhausted");
                    id = next++;
                    ++new_charged;
                }
                auto native = r::charged_native_view<GapParticle>(t.catalog(), packet, id);
                native_particle_guard(native);
                mfem::IntegrationPoint ip;
                static_cast<void>(p.locate(native.position_m[0], native.position_m[1], ip));
                require(!p.inside_conductor(native.position_m[0], native.position_m[1]) &&
                            (!p.options_.finite_axial_extent ||
                             (native.position_m[2] >= p.options_.short_z_m &&
                              native.position_m[2] <= p.options_.feed_z_m)),
                        "reaction product outside vacuum");
                require(impl->bindings_after.emplace(packet.id, id).second, "product binding collision");
                replacement.push_back(native);
            }
    }
    require(replacement.size() <= 2000000, "native particle resource guard");
    auto work = capture(p);
    work.particles = std::move(replacement);
    work.next_id = next;
    work.fields_current = false;
    const auto before_hash = fingerprint(p);
    bool workspace_active = false;
    try {
        swap_state(p, work);
        workspace_active = true;
        fail_at(fault, NativeReactionFault::AfterPopulationStage);
        p.solve_fields();
        fail_at(fault, NativeReactionFault::AfterFieldSolve);
        impl->after = capture(p);
        swap_state(p, work);
        workspace_active = false;
        impl->work = std::move(work);
        impl->bindings_work = impl->bindings_after;
        fail_at(fault, NativeReactionFault::BeforePreparedReturn);
    } catch (...) {
        if (workspace_active)
            swap_state(p, work);
        throw;
    }
    require(matches(p, impl->before) && fingerprint(p) == before_hash,
            "native preparation failed to restore original state");
    require(
        impl->after.ledger.authored_charge_C == impl->before.ledger.authored_charge_C &&
            impl->after.ledger.authored_kinetic_energy_J ==
                impl->before.ledger.authored_kinetic_energy_J &&
            impl->after.ledger.lower_injected_charge_C == impl->before.ledger.lower_injected_charge_C &&
            impl->after.ledger.upper_injected_charge_C == impl->before.ledger.upper_injected_charge_C &&
            impl->after.ledger.injected_energy_J == impl->before.ledger.injected_energy_J,
        "reaction changed external source accounting");
    impl->audit = {{"status", "manufactured_native_field_stage"},
                   {"event_id", e.delta().event_id},
                   {"before_native_sha256", before_hash},
                   {"catalog_sha256", t.catalog().digest()},
                   {"new_charged_ids", new_charged},
                   {"particles_before", impl->before.particles.size()},
                   {"particles_after", impl->after.particles.size()},
                   {"fields_invalidated_before_actual_solve", true},
                   {"actual_CPU_Poisson_solve", true},
                   {"authored_source_unchanged", true},
                   {"phase", phase.kind == ReactionPhase::Kind::InitialSynchronizedFieldProof
                                 ? "initial_synchronized_field_proof"
                                 : "provisional_equal_dt_mid_x_boundary"},
                   {"MCC_heating_qualification", false},
                   {"Jxyz_EM_deposition_qualified", false}};
    return PreparedPicReaction(std::move(impl));
}
bool PicReactionBatch::can_publish(const GapPic &p, const r::Transactions &t,
                                   const ReactionPacketBindings &binding,
                                   const PreparedPicReaction &token) noexcept {
    if (!token.impl_)
        return false;
    const auto &i = *token.impl_;
    return i.status == PreparedPicReaction::Impl::Status::Staged && i.pic_owner == &p &&
           i.reaction_owner == &t && i.binding_owner == &binding && t.can_publish(i.reaction) &&
           binding == i.bindings_before && matches(p, i.before);
}
bool PicReactionBatch::try_publish(GapPic &p, r::Transactions &t, ReactionPacketBindings &binding,
                                   PreparedPicReaction &token) noexcept {
    if (!can_publish(p, t, binding, token))
        return false;
    auto &i = *token.impl_;
    swap_state(p, i.work);
    binding.swap(i.bindings_work);
    if (!t.try_publish(i.reaction))
        std::terminate();
    i.status = PreparedPicReaction::Impl::Status::Published;
    return true;
}
bool PicReactionBatch::try_rollback(GapPic &p, r::Transactions &t, ReactionPacketBindings &binding,
                                    PreparedPicReaction &token) noexcept {
    if (!token.impl_)
        return false;
    auto &i = *token.impl_;
    if (i.status != PreparedPicReaction::Impl::Status::Published || i.pic_owner != &p ||
        i.reaction_owner != &t || i.binding_owner != &binding || t.snapshot() != i.reaction.output() ||
        binding != i.bindings_after || !matches(p, i.after))
        return false;
    swap_state(p, i.work);
    binding.swap(i.bindings_work);
    if (!t.try_revert(i.reaction))
        std::terminate();
    i.status = PreparedPicReaction::Impl::Status::RolledBack;
    return true;
}
} // namespace raikou::mfem_backend
