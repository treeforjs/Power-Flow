// Native v7 forms against the sealed native and Hall functions:
//  * pf7::kinetic_energy = sealed native.o pic_kinetic_energy_J (bitwise);
//  * pf7::first_hit = sealed Hall first_kinetic_hit (8dddec0) on interior starts;
//  * pf_v7::neutral_stream(dt_max = 1e-14) = sealed prepare_neutral_stream
//    (values, refusal texts and order); a larger configured bound admits dt
//    above 1e-14 with the same arithmetic;
//  * native configuration: the default bound keeps the v6 configuration
//    identity, a non-default bound changes it and survives a checkpoint;
//  * native admission: a step with a resident neutral above the configured
//    bound refuses; the same step under a larger configured bound is admitted.
// usage: test_v7_forms_native OUTPUT_DIRECTORY [reduced]
#include "v6_fixture.hpp"
#include <filesystem>
#include <random>

using namespace v6;
using FK = m::FullTimeFlowRefusal::Kind;
namespace {
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }
pf_v7::KineticPolicy policy_of(const m::GapPicOptions &o) {
  return {o.width_m, o.gap_m, o.short_z_m, o.feed_z_m,
          {o.lower_foil.left_m, o.lower_foil.right_m, o.lower_foil.bottom_m, o.lower_foil.top_m},
          {o.upper_foil.left_m, o.upper_foil.right_m, o.upper_foil.bottom_m, o.upper_foil.top_m}};
}
pf7::Policy core_policy(const m::GapPicOptions &o) {
  return {o.width_m, o.gap_m, o.short_z_m, o.feed_z_m,
          {o.lower_foil.left_m, o.lower_foil.right_m, o.lower_foil.bottom_m, o.lower_foil.top_m},
          {o.upper_foil.left_m, o.upper_foil.right_m, o.upper_foil.bottom_m, o.upper_foil.top_m}};
}
} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 2)
      throw std::runtime_error("usage: test_v7_forms_native OUTPUT_DIRECTORY [reduced]");
    const std::string dir = argv[1];
    const bool reduced = argc > 2 && std::string(argv[2]) == "reduced";
    const std::size_t samples = reduced ? 20000 : 400000;
    J measurements = J::object();
    std::mt19937_64 g(20261031);
    std::uniform_real_distribution<double> U(0, 1);
    const auto o = options();

    // ---- KE parity.
    {
      std::uint64_t bad = 0;
      const double masses[] = {m::pic_me, 2.9915e-26, 1.6726e-27};
      for (std::size_t i = 0; i < samples; ++i) {
        m::GapParticle p;
        p.mass_kg = i % 7 == 6 ? std::ldexp(1 + U(g), int(g() % 40) - 110) : masses[i % 3];
        for (int d = 0; d < 3; ++d) {
          const int e = int(g() % 50) - 20; // 1e-6 .. 1e9 m/s scale
          double v = (U(g) - .5) * std::ldexp(1.0, e);
          if (g() % 17 == 0)
            v = 0;
          if (g() % 101 == 0)
            v = std::ldexp(1.0, -1074) * double(g() % 7);
          if (g() % 53 == 0)
            v = (U(g) - .5) * 4e8;
          p.proper_velocity_m_s[d] = v;
        }
        const double a = m::pic_kinetic_energy_J(p);
        const double b = pf7::kinetic_energy(p.mass_kg, p.proper_velocity_m_s[0], p.proper_velocity_m_s[1],
                                             p.proper_velocity_m_s[2]);
        bad += bits(a) != bits(b);
      }
      check(bad == 0, "pf7::kinetic_energy equals the sealed pic_kinetic_energy_J bit for bit (" +
                          std::to_string(samples) + " samples)");
      measurements["kinetic_energy_parity"] = {{"samples", samples}, {"mismatches", bad}};
    }

    // ---- First-hit parity (interior starts).
    const auto hp = hall_driver::declared_kinetic_policy(o);
    const auto cp = core_policy(o);
    {
      std::uint64_t bad = 0, n = 0, hits = 0;
      while (n < samples) {
        m::PicVector x{o.width_m * U(g), o.gap_m * U(g), o.short_z_m + o.depth_m * U(g)};
        try {
          hall_driver::require_kinetic_interior(hp, x);
        } catch (const std::exception &) {
          continue;
        }
        const double scale = std::ldexp(1.0, -int(g() % 20) - 3);
        m::PicVector D{(U(g) - .5) * scale, (U(g) - .5) * scale, (U(g) - .5) * scale * 4};
        if (g() % 9 == 0)
          D[int(g() % 3)] = 0;
        const auto h = hall_driver::first_kinetic_hit(hp, x, D);
        const auto c = pf7::first_hit(cp, x.data(), D.data());
        bad += int(h.sink) != c.sink || bits(h.fraction) != bits(c.fraction) || h.axis != c.axis ||
               bits(h.plane_m) != bits(c.plane);
        hits += c.sink >= 0;
        ++n;
      }
      check(bad == 0, "pf7::first_hit equals the sealed first_kinetic_hit (" + std::to_string(n) + " segments)");
      measurements["first_hit_parity"] = {{"segments", n}, {"hits", hits}, {"mismatches", bad}};
    }

    // ---- Neutral stream parity at the default bound; configured bound.
    {
      const auto pol = policy_of(o);
      std::uint64_t bad = 0, refused = 0, n = 0, larger_ok = 0, larger_bad = 0;
      for (std::size_t i = 0; i < samples / 4; ++i) {
        m::GapParticle p;
        p.id = 1 + g() % 1000;
        p.charge_C = g() % 50 == 0 ? -m::pic_e : 0;
        p.mass_kg = 2.9915e-26;
        p.weight = 1e9 * (1 + U(g));
        p.position_m = {o.width_m * U(g), o.gap_m * U(g), o.short_z_m + o.depth_m * U(g)};
        const double vs = std::ldexp(1.0, int(g() % 30));
        p.proper_velocity_m_s = {(U(g) - .5) * vs, (U(g) - .5) * vs, (U(g) - .5) * vs};
        if (g() % 200 == 0)
          p.proper_velocity_m_s[1] = std::numeric_limits<double>::infinity();
        const double begin = 2.6e-8 * U(g);
        const double dt = g() % 10 == 0 ? 1.0000000000000002e-14 + 1e-14 * U(g) : 1e-14 * (U(g) + 1e-3);
        const double end = begin + dt;
        std::string ha, ca;
        hall_driver::NeutralStreamProposal hs;
        pf_v7::NeutralStreamResult cs;
        try {
          hs = hall_driver::prepare_neutral_stream(hp, p, begin, end, dt);
        } catch (const std::exception &e) {
          ha = e.what();
        }
        const pf_v7::Carrier c{p.id, {p.position_m[0], p.position_m[1], p.position_m[2]},
                               {p.proper_velocity_m_s[0], p.proper_velocity_m_s[1], p.proper_velocity_m_s[2]},
                               p.charge_C, p.mass_kg, p.weight};
        try {
          cs = pf_v7::neutral_stream(pol, c, begin, end, dt);
        } catch (const std::exception &e) {
          ca = e.what();
        }
        ++n;
        if (!ha.empty() || !ca.empty()) {
          bad += ha != ca;
          ++refused;
        } else {
          bool same = int(hs.hit.sink) == cs.hit.sink && bits(hs.hit.fraction) == bits(cs.hit.fraction) &&
                      hs.hit.axis == cs.hit.axis && bits(hs.hit.plane_m) == bits(cs.hit.plane) &&
                      bits(hs.flight_s) == bits(cs.flight_s) && bits(hs.remaining_s) == bits(cs.remaining_s) &&
                      bits(hs.incident_kinetic_J) == bits(cs.incident_kinetic_J);
          for (int d = 0; d < 3; ++d)
            same = same && bits(hs.endpoint.position_m[d]) == bits(cs.endpoint.position_m[d]);
          bad += !same;
        }
        // Larger configured bound: dt up to 2e-14 admitted, arithmetic = the
        // sealed function at dt (evaluated by scaling the interval check only).
        if (ha == "neutral authoritative interval" && dt <= 2e-14 && std::isfinite(begin) && end == begin + dt) {
          try {
            const auto r = pf_v7::neutral_stream(pol, c, begin, end, dt, 2e-14);
            pf7::NeutralStream s{};
            (void)pf7::neutral_stream(cp, c.id, c.position_m, c.proper_velocity_m_s, c.charge_C, c.mass_kg,
                                      c.weight, begin, end, dt, 2e-14, s);
            double u2 = 0;
            for (double v : c.proper_velocity_m_s)
              u2 += v * v;
            const double gamma = std::sqrt(1 + u2 / (m::pic_c * m::pic_c));
            bool same = true;
            for (int d = 0; d < 3; ++d)
              same = same && bits(r.displacement_m[d]) == bits(dt * c.proper_velocity_m_s[d] / gamma);
            larger_ok += same;
            larger_bad += !same;
          } catch (const std::exception &) {
          }
        }
      }
      check(bad == 0, "pf_v7::neutral_stream (default bound) equals the sealed prepare_neutral_stream: values, "
                      "texts and order (" + std::to_string(n) + " cases, " + std::to_string(refused) + " refused)");
      check(larger_ok > 0 && larger_bad == 0, "a 2e-14 bound admits dt in (1e-14, 2e-14] with the sealed arithmetic");
      measurements["neutral_stream_parity"] = {{"cases", n}, {"refused", refused}, {"mismatches", bad},
                                               {"larger_bound_admitted", larger_ok}};
    }

    // ---- Configuration identity and checkpoint.
    {
      auto fresh_with = [&](double bound) {
        m::FullTimeConfiguration c;
        c.motion.neutral_stream_dt_max_s = bound;
        return m::PicArena::create_full_time(options(), states(), channels(), base_population(), genesis(), c);
      };
      const auto plain = m::PicArena::create_full_time(options(), states(), channels(), base_population(), genesis());
      const auto d = fresh_with(1e-14), w = fresh_with(2e-14);
      check(plain->state_prefix() == d->state_prefix() &&
                plain->full_time_boundary()->accepted_arena_sha256 == d->full_time_boundary()->accepted_arena_sha256,
            "default neutral bound keeps the v6 configuration identity");
      check(w->full_time_boundary()->accepted_arena_sha256 != d->full_time_boundary()->accepted_arena_sha256 &&
                w->full_time_motion_policy().neutral_stream_dt_max_s == 2e-14,
            "a non-default neutral bound is part of the configuration identity");
      check(refuses([&] { (void)fresh_with(0); }) && refuses([&] { (void)fresh_with(std::nan("")); }),
            "non-positive or nonfinite neutral bound refuses");
      const auto path = dir + "/bound-2e-14.checkpoint.json";
      w->save(path);
      auto restored = m::PicArena::restore_full_time(path, options(), states(), channels(), *w->full_time_boundary());
      check(restored->full_time_motion_policy().neutral_stream_dt_max_s == 2e-14 &&
                restored->full_time_boundary()->accepted_arena_sha256 == w->full_time_boundary()->accepted_arena_sha256,
            "the configured bound round-trips through a checkpoint");
      const auto path0 = dir + "/bound-default.checkpoint.json";
      d->save(path0);
      std::ifstream in(path0);
      const auto text = std::string(std::istreambuf_iterator<char>(in), {});
      check(text.find("neutral_stream_dt_max_s") == std::string::npos,
            "a default-bound checkpoint has the v6 configuration keys");
    }

    // ---- Native admission of the configured bound (reference oracle,
    //      resident neutral carrier at rest, dt = 1.5e-14).
    {
      auto run = [&](double bound, double dt) {
        m::FullTimeConfiguration c;
        c.motion.neutral_stream_dt_max_s = bound;
        auto a = m::PicArena::create_full_time(options(), states(), channels(), base_population(), genesis(), c);
        StepState st;
        StepPlan plan;
        plan.dt_s = dt;
        full_step(*a, plan, st, {nullptr, reference_hook(64, [bound](ReferenceMotion &r) {
                                   r.record.neutral_stream_dt_max_s = bound;
                                 }), nullptr});
        return a->full_time_boundary()->accepted_step;
      };
      check(flow_refuses([&] { (void)run(1e-14, 1.5e-14); }, FK::motion_record),
            "default bound: a resident neutral over dt = 1.5e-14 refuses");
      check(run(2e-14, 1.5e-14) == 1, "bound 2e-14: the same step is admitted");
      check(run(1e-14, 1e-14) == 1, "default bound: dt = 1e-14 admitted");
    }
    J output = {{"status", "passed"}, {"check_count", checks}, {"checks", results}, {"measurements", measurements}};
    std::ofstream(dir + "/v7-forms-native.json") << output.dump(2) << '\n';
    std::cout << checks << " native v7 forms controls passed\n";
    return 0;
  } catch (const std::exception &x) {
    std::cerr << "FAILED after " << checks << " checks: " << x.what() << '\n';
    return 1;
  }
}
