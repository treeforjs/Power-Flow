#include "raikou/mfem_pic_reaction_batch.h"
#include "manufactured_fixtures.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace m = raikou::mfem_backend;
namespace r = raikou::pic_reactions;
namespace {
r::Json checks = r::Json::array();
void check(bool ok, const std::string &name) {
    checks.push_back({{"name", name}, {"passed", ok}});
    if (!ok)
        throw std::runtime_error(name);
}
template <class F> void reject(F fn, const std::string &name) {
    bool caught = false;
    try {
        fn();
    } catch (const std::exception &) {
        caught = true;
    }
    check(caught, name);
}
struct Scenario {
    r::Catalog catalog{manufactured::definition()};
    std::unique_ptr<m::GapPic> pic;
    std::unique_ptr<r::Transactions> reactions;
    m::ReactionPacketBindings binding;
    Scenario(double preceding_step = 0) {
        m::GapPicOptions o;
        o.nx = o.ny = 45;
        o.width_m = o.gap_m = .009;
        o.depth_m = .025;
        o.finite_foils = true;
        o.finite_axial_extent = true;
        o.lower_foil = {.002, .007, .0018, .002};
        o.upper_foil = {.002, .007, .007, .0072};
        pic = std::make_unique<m::GapPic>(o);
        auto g = manufactured::genesis(catalog);
        std::vector<m::GapParticle> seed;
        for (const auto &packet : g.packets) {
            m::GapParticle p;
            p.position_m = packet.position_m;
            p.proper_velocity_m_s = packet.proper_velocity_m_s;
            p.weight = packet.weight;
            seed.push_back(p);
        }
        auto ids = pic->add_authored_batch(seed);
        pic->solve_fields();
        if (preceding_step > 0)
            pic->advance(preceding_step);
        for (std::size_t i = 0; i < ids.size(); ++i) {
            binding.emplace(g.packets[i].id, ids[i]);
            g.packets[i].position_m = pic->particles()[i].position_m;
            g.packets[i].proper_velocity_m_s = pic->particles()[i].proper_velocity_m_s;
        }
        g.parents[0].position_m = g.packets[0].position_m;
        auto d = pic->diagnostics();
        g.authored_charge_C = d.authored_charge_C;
        g.authored_kinetic_J = d.authored_kinetic_energy_J;
        g.coupling_prefix = "manufactured/native-input/" + m::PicReactionBatch::fingerprint(*pic);
        reactions = std::make_unique<r::Transactions>(catalog, g);
    }
    r::ForcedEvent event(const std::string &channel) {
        auto e = manufactured::forced(*reactions, channel, "primary/0");
        e.time_s = pic->diagnostics().time_s;
        return e;
    }
};
std::vector<double> nodal(const m::GapPic &p) {
    const auto &o = p.options();
    std::vector<double> q((o.nx + 1) * (o.ny + 1), 0);
    for (const auto &part : p.particles()) {
        double sx = part.position_m[0] / o.width_m * o.nx, sy = part.position_m[1] / o.gap_m * o.ny;
        int ix = std::min(int(sx), o.nx - 1), iy = std::min(int(sy), o.ny - 1);
        double x = sx - ix, y = sy - iy, w = part.charge_C * part.weight;
        q[iy * (o.nx + 1) + ix] += w * (1 - x) * (1 - y);
        q[iy * (o.nx + 1) + ix + 1] += w * x * (1 - y);
        q[(iy + 1) * (o.nx + 1) + ix] += w * (1 - x) * y;
        q[(iy + 1) * (o.nx + 1) + ix + 1] += w * x * y;
    }
    return q;
}
bool same_sources(const m::GapPicDiagnostics &a, const m::GapPicDiagnostics &b) {
    return a.authored_charge_C == b.authored_charge_C &&
           a.authored_kinetic_energy_J == b.authored_kinetic_energy_J &&
           a.lower_injected_charge_C == b.lower_injected_charge_C &&
           a.upper_injected_charge_C == b.upper_injected_charge_C &&
           a.injected_energy_J == b.injected_energy_J;
}
} // namespace
int main(int argc, char **argv) {
    try {
        if (argc != 2)
            throw std::runtime_error("evidence directory required");
        const std::string out = argv[1];
        r::Json audits = r::Json::array();
        for (double weight : {.1, 1e-20}) {
            r::Catalog catalog(manufactured::definition());
            auto g = manufactured::genesis(catalog);
            g.authored_charge_C = 0;
            g.authored_kinetic_J = 0;
            for (auto &packet : g.packets) {
                packet.weight = weight;
                g.authored_charge_C -= weight * r::electron_charge_C;
                g.authored_kinetic_J +=
                    weight * r::kinetic_J(catalog.species("e"), packet.proper_velocity_m_s);
            }
            r::Transactions exact_debit(catalog, g);
            auto before = exact_debit.checkpoint_document();
            auto event = manufactured::forced(exact_debit, "manufactured/ionization", "primary/0");
            reject([&] { exact_debit.prepare(event); },
                   "unrepresentable one-double parent debit is explicitly rejected");
            check(exact_debit.checkpoint_document() == before,
                  "failed exact-debit admission preserves every inventory/packet/journal value");
        }
        for (const std::string channel :
             {"manufactured/ionization", "manufactured/dissociation", "manufactured/attachment"}) {
            Scenario s;
            const auto old_native = m::PicReactionBatch::fingerprint(*s.pic);
            const auto old_core = s.reactions->checkpoint_document();
            const auto old_binding = s.binding;
            const auto sources = s.pic->diagnostics();
            const auto rho0 = nodal(*s.pic);
            const double K0 = sources.particle_kinetic_energy_J;
            auto event = s.event(channel);
            auto prepared = s.reactions->prepare(event);
            m::PrimaryWeightSlice slice{event.weight, 0};
            m::ReactionPhase phase;
            for (auto fault :
                 {m::NativeReactionFault::AfterPopulationStage, m::NativeReactionFault::AfterFieldSolve,
                  m::NativeReactionFault::BeforePreparedReturn}) {
                reject(
                    [&] {
                        m::PicReactionBatch::prepare(*s.pic, *s.reactions, prepared, s.binding, phase,
                                                     slice, fault);
                    },
                    channel + " staged native fault");
                check(m::PicReactionBatch::fingerprint(*s.pic) == old_native &&
                          s.reactions->checkpoint_document() == old_core && s.binding == old_binding,
                      channel + " complete native/core/sidecar rollback on preparation failure");
            }
            auto partial = slice;
            partial.event_weight *= .5;
            partial.untouched_primary_native_id = 99;
            reject(
                [&] {
                    m::PicReactionBatch::prepare(*s.pic, *s.reactions, prepared, s.binding, phase,
                                                 partial);
                },
                channel + " partial-primary v1 path explicitly rejected");
            auto wrong = s.binding;
            wrong["primary/1"] = wrong.at("primary/0");
            reject(
                [&] {
                    m::PicReactionBatch::prepare(*s.pic, *s.reactions, prepared, wrong, phase, slice);
                },
                channel + " duplicate/missing native ID binding rejected");
            auto native =
                m::PicReactionBatch::prepare(*s.pic, *s.reactions, prepared, s.binding, phase, slice);
            auto competing =
                m::PicReactionBatch::prepare(*s.pic, *s.reactions, prepared, s.binding, phase, slice);
            check(m::PicReactionBatch::fingerprint(*s.pic) == old_native &&
                      s.reactions->checkpoint_document() == old_core,
                  channel + " native field staging has no publication");
            check(m::PicReactionBatch::can_publish(*s.pic, *s.reactions, s.binding, native) &&
                      m::PicReactionBatch::try_publish(*s.pic, *s.reactions, s.binding, native),
                  channel + " atomic native/core/binding publication");
            const auto post = s.pic->diagnostics();
            check(same_sources(sources, post),
                  channel + " external authored/thermal source ledgers unchanged");
            check(s.binding == native.output_bindings() && s.binding.size() == s.pic->particles().size(),
                  channel + " all native charged IDs have immutable state bindings");
            check(s.pic->particles().size() == (channel == "manufactured/ionization" ? 5u : 3u),
                  channel + " actual native product population");
            check(post.particle_charge_C == sources.particle_charge_C &&
                      std::abs(post.particle_charge_defect_C) < 1e-30,
                  channel + " native charge closure");
            const auto rho1 = nodal(*s.pic);
            double max_rho = 0;
            for (std::size_t i = 0; i < rho0.size(); ++i)
                max_rho = std::max(max_rho, std::abs(rho1[i] - rho0[i]));
            check(max_rho <=
                      256 * std::numeric_limits<double>::epsilon() * std::abs(sources.particle_charge_C),
                  channel + " same-point Q1 nodal charge conservation");
            check(std::isfinite(post.poisson_residual_C_per_m) && std::isfinite(post.gauss_defect_C) &&
                      std::isfinite(s.pic->electric_field(.003, .004)[0]),
                  channel + " actual native Poisson/deposition/field gather succeeds");
            const auto &certificate = prepared.delta().conservation;
            double expected_delta = certificate.at("after").at("electron_kinetic_J").get<double>() +
                                    certificate.at("after").at("ion_kinetic_J").get<double>() -
                                    certificate.at("before").at("electron_kinetic_J").get<double>();
            check(std::abs((post.particle_kinetic_energy_J - K0) - expected_delta) <
                      256 * std::numeric_limits<double>::epsilon() * K0,
                  channel + " native resident K matches internal event receipt");
            check(!m::PicReactionBatch::try_publish(*s.pic, *s.reactions, s.binding, competing) &&
                      !m::PicReactionBatch::try_publish(*s.pic, *s.reactions, s.binding, native),
                  channel + " stale/duplicate native publication rejected");
            audits.push_back(native.audit());
            const auto bundle = out + "/" + channel.substr(channel.find('/') + 1) + ".bundle.json";
            m::PicReactionBatch::save_bundle(*s.pic, *s.reactions, s.binding, bundle);
            auto restarted = m::PicReactionBatch::load_bundle(
                s.catalog, s.pic->options(), s.reactions->genesis().coupling_prefix, bundle);
            check(restarted.reactions->checkpoint_document() == s.reactions->checkpoint_document() &&
                      restarted.bindings == s.binding,
                  channel + " native/core/binding bundle restart provenance");
            bool exact_particles = restarted.pic->particles().size() == s.pic->particles().size();
            for (std::size_t i = 0; i < s.pic->particles().size(); ++i) {
                const auto &a = s.pic->particles()[i];
                const auto &b = restarted.pic->particles()[i];
                exact_particles = exact_particles && a.id == b.id && a.position_m == b.position_m &&
                                  a.proper_velocity_m_s == b.proper_velocity_m_s &&
                                  a.charge_C == b.charge_C && a.mass_kg == b.mass_kg &&
                                  a.weight == b.weight;
            }
            check(exact_particles && same_sources(post, restarted.pic->diagnostics()),
                  channel + " native binary population/source counters restored exactly");
            reject([&] { restarted.reactions->prepare(event); },
                   channel + " accepted reaction replay rejected after native bundle restart");
            std::ifstream input(bundle);
            r::Json document;
            input >> document;
            for (int mutation = 0; mutation < 9; ++mutation) {
                auto bad = document;
                auto &body = bad["payload"];
                if (mutation == 0)
                    bad["sha256"] = "wrong";
                if (mutation == 1)
                    body["native_source_sha256"] = std::string(64, '0');
                if (mutation == 2)
                    body["packet_bindings"].begin().value() = 999999;
                if (mutation == 3)
                    body["native_bytes_sha256"] = "wrong";
                if (mutation == 4)
                    body["reaction_checkpoint"]["payload"]["snapshot"]["parents"][0]
                        ["remaining_weight"] = .5;
                if (mutation == 5)
                    body["saved_potential_V"][0] = 1.;
                if (mutation == 6)
                    body["native_abi"][0] = 81;
                if (mutation == 7)
                    body["schema"] = "physical/unqualified";
                if (mutation == 8)
                    body["extra"] = "unrecognized";
                if (mutation != 0)
                    bad["sha256"] = r::sha256(body.dump());
                const auto path = out + "/corrupt.bundle.json";
                std::ofstream(path) << bad.dump();
                reject(
                    [&] {
                        m::PicReactionBatch::load_bundle(s.catalog, s.pic->options(),
                                                         s.reactions->genesis().coupling_prefix, path);
                    },
                    channel + " rejects native bundle corruption/" + std::to_string(mutation));
            }
            reject(
                [&] {
                    m::PicReactionBatch::load_bundle(s.catalog, s.pic->options(),
                                                     "different-native-prefix", bundle);
                },
                channel + " bundle requires matching native/core initial prefix");
            check(m::PicReactionBatch::try_rollback(*s.pic, *s.reactions, s.binding, native),
                  channel + " joint-barrier rollback after publication");
            check(m::PicReactionBatch::fingerprint(*s.pic) == old_native &&
                      s.reactions->checkpoint_document() == old_core && s.binding == old_binding,
                  channel + " original solved native/core/sidecar state restored exactly");
            check(!m::PicReactionBatch::try_rollback(*s.pic, *s.reactions, s.binding, native) &&
                      !m::PicReactionBatch::try_publish(*s.pic, *s.reactions, s.binding, native),
                  channel + " consumed native token cannot be reused");
        }
        Scenario phase_case(1e-13);
        auto event = phase_case.event("manufactured/ionization");
        auto prepared = phase_case.reactions->prepare(event);
        m::ReactionPhase phase{m::ReactionPhase::Kind::MidXEqualStepBoundary,
                               phase_case.pic->diagnostics().time_s, 1e-13, 1e-13, 1};
        m::PrimaryWeightSlice slice{event.weight, 0};
        auto invalid = phase;
        invalid.next_dt_s = 2e-13;
        reject(
            [&] {
                m::PicReactionBatch::prepare(*phase_case.pic, *phase_case.reactions, prepared,
                                             phase_case.binding, invalid, slice);
            },
            "unequal adjacent-step mid-x claim rejected");
        auto mid = m::PicReactionBatch::prepare(*phase_case.pic, *phase_case.reactions, prepared,
                                                phase_case.binding, phase, slice);
        check(m::PicReactionBatch::try_publish(*phase_case.pic, *phase_case.reactions,
                                               phase_case.binding, mid),
              "actual synchronized equal-dt native boundary is exposed explicitly");
        check(m::PicReactionBatch::try_rollback(*phase_case.pic, *phase_case.reactions,
                                                phase_case.binding, mid),
              "mid-x staged reaction rollback preserves existing force history");
        Scenario baseline, reactive;
        auto reaction = reactive.event("manufactured/ionization");
        auto core = reactive.reactions->prepare(reaction);
        m::ReactionPhase initial;
        auto native = m::PicReactionBatch::prepare(*reactive.pic, *reactive.reactions, core,
                                                   reactive.binding, initial, {reaction.weight, 0});
        check(m::PicReactionBatch::try_publish(*reactive.pic, *reactive.reactions, reactive.binding,
                                               native),
              "field-feedback fixture commit");
        baseline.pic->advance(1e-13);
        reactive.pic->advance(1e-13);
        double difference = 0, scale = 0;
        for (int i = 0; i < baseline.pic->potential().Size(); ++i) {
            difference = std::max(difference,
                                  std::abs(reactive.pic->potential()[i] - baseline.pic->potential()[i]));
            scale = std::max(scale, std::abs(baseline.pic->potential()[i]));
        }
        check(difference > 1e-8 * scale,
              "actual charged-product transport changes self-consistent native potential");
        check(same_sources(baseline.pic->diagnostics(), reactive.pic->diagnostics()) &&
                  std::abs(reactive.pic->diagnostics().particle_charge_defect_C) < 1e-30,
              "post-transport field response preserves external charge bookkeeping");
        check(!m::PicReactionBatch::try_rollback(*reactive.pic, *reactive.reactions, reactive.binding,
                                                 native),
              "rollback rejects native state changed by a later advance");
        reject(
            [&] {
                m::PicReactionBatch::save_bundle(*reactive.pic, *reactive.reactions, reactive.binding,
                                                 out + "/unsynchronized.bundle.json");
            },
            "bundle refuses stale core views after native transport");
        r::Json report = {{"status", "passed"},
                          {"check_count", checks.size()},
                          {"checks", checks},
                          {"audits", audits},
                          {"manufactured_only", true},
                          {"native_source_baseline", "660d1c80fb368af9ffef5244db4636a64a92210f"},
                          {"potential_response_max_abs_V", difference},
                          {"potential_reference_scale_V", scale},
                          {"physical_rates_or_breakdown_claim", false},
                          {"MCC_heating_or_EM_Jxyz_qualification", false},
                          {"ABI_sizes",
                           {{"GapPic", sizeof(m::GapPic)},
                            {"Particle", sizeof(m::GapParticle)},
                            {"Options", sizeof(m::GapPicOptions)},
                            {"Diagnostics", sizeof(m::GapPicDiagnostics)}}}};
        std::ofstream(out + "/native-tests.json") << report.dump(2) << '\n';
        std::cout << "Passed " << checks.size() << " manufactured native reaction/Poisson checks\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n' << checks.dump(2) << '\n';
        return 1;
    }
}
