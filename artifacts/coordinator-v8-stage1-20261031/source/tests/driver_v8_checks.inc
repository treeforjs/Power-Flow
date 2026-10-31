// Coordinator v8 qualification checks (included by driver_v8.cpp): the coordinator v7.1 checks
// (T1-T4, T6, spot, V7, V7.1) on the host Hall owner, with source owner v7 at benchmark options.
struct Results {
  J cases = J::array();
  int failed = 0;
  void record(const std::string& test, const std::string& name, bool ok, J ev = J::object()) {
    cases.push_back({{"test", test}, {"name", name}, {"passed", ok}, {"evidence", std::move(ev)}});
    if (!ok) ++failed;
    std::cout << test << " | " << name << ": " << (ok ? "pass" : "FAIL") << std::endl;
  }
};
long scratch_counter = 0;
template <class H> std::string owner_digest(const v7::OwnersV8<H>& o, const fs::path& scratch) {
  v7::Encoding e;
  e.text("coordinator-v6-owner-snapshot-v1");
  const auto unique = scratch / ("snap-" + std::to_string(::getpid()) + "-" + std::to_string(++scratch_counter));
  fs::create_directories(unique);
  o.arena->save((unique / "native.json").string());
  e.text(nw::sha(read_bytes(unique / "native.json")));
  o.hall->save(unique / "hall");
  std::vector<fs::path> files;
  for (const auto& x : fs::recursive_directory_iterator(unique / "hall"))
    if (x.is_regular_file()) files.push_back(x.path());
  std::sort(files.begin(), files.end());
  for (const auto& f : files) { e.text(fs::relative(f, unique).string()); e.text(nw::sha(read_bytes(f))); }
  fs::remove_all(unique);
  e.text(o.arena->state_prefix());
  e.text(o.hall->field_sha256()); e.text(o.hall->heat_sha256()); e.text(o.hall->kinetic_ledgers_sha256());
  e.text(o.chemistry->checkpoint().dump());
  for (const auto& s : o.sources) e.text(s->state_sha256());
  e.text(o.arrival->state_sha256()); e.text(o.arrival->ledger_head()); e.text(o.arrival->save());
  e.text(o.accepted_source_seal);
  return e.finish();
}
int checks_mode(const fs::path& out) {
  using HO = hc::Owner;
  Results R;
  const fs::path scratch = out / "scratch";
  fs::create_directories(scratch);
  const auto in = load_inputs("qual");
  const double dt = opt.dt_s;
  auto make = [&] { return genesis<HO>(in, "qual"); };
  // T1: one positive step, zero allocations in the swap window, second step preflights; reference digest.
  std::string reference_after;
  {
    auto o = make();
    const auto d0 = owner_digest(o, scratch);
    auto t = v7::prepare_v8(o, dt);
    const auto n = t.numbers();
    const auto pub = publish_measured(o, t);
    reference_after = owner_digest(o, scratch);
    R.record("T1", "positive coupled step published (births, motion, sinks, exposures)", pub.published,
             {{"electron_births", n.electron_births}, {"h2o_births", n.neutral_births}, {"sinks", n.sinks}, {"parent_events", n.parent_events},
              {"census", n.census}});
    R.record("T1", "zero allocations between the first and last swap", pub.published && pub.swap_allocations == 0, {{"allocations", pub.swap_allocations}});
    R.record("T1", "accepted owners changed by the step", d0 != reference_after);
    auto t2 = v7::prepare_v8(o, dt);
    R.record("T1", "a second step prepares and preflights", v7::can_publish_v8(o, t2));
    R.record("T1", "no legacy native solve/advance call", forbidden_calls == 0);
  }
  // T2: mutations after sealing; each refuses with every owner unchanged, then the
  // unmutated step reproduces the T1 reference.
  {
    auto o = make();
    auto twin = make();
    // A twin whose arrival policy has another run seed (different plan and token).
    auto in2 = in;
    in2.settings.arrival_seed += 1;
    in2.setup = v7::make_source_setup(in.model->mesh(), in2.settings);
    auto twin_seed = genesis<HO>(in2, "qual");
    const auto d0 = owner_digest(o, scratch);
    using TR = v7::TrialV8<HO>;
    struct Case { std::string name; std::function<void(TR&, TR&)> mutate; };
    std::vector<Case> cases = {
        {"foreign source transition seal", [](TR& t, TR&) { TA::source_after(t) = nw::sha("foreign source seal"); }},
        {"Hall proposal of a twin owner set (dt/2)", [](TR& t, TR& x) { std::swap(TA::proposal(t), TA::proposal(x)); }},
        {"source trial of a twin source owner", [](TR& t, TR& x) { std::swap(TA::source_trials(t)[0], TA::source_trials(x)[0]); }},
        {"arrival staging of a twin arrival owner with another run seed", [](TR& t, TR& x) { std::swap(TA::staged(t), TA::staged(x)); }},
        {"native token of a twin arena", [](TR& t, TR& x) { std::swap(TA::native(t), TA::native(x)); }},
    };
    for (const auto& c : cases) {
      auto t = v7::prepare_v8(o, dt);
      auto& other = c.name.find("run seed") != std::string::npos ? twin_seed : twin;
      auto x = v7::prepare_v8(other, c.name.find("dt/2") != std::string::npos ? .5 * dt : dt);
      c.mutate(t, x);
      const bool refused = !v7::try_publish_v8(o, t);
      const int point = t.refusal_point();
      const bool unchanged = owner_digest(o, scratch) == d0;
      R.record("T2", c.name + ": refused, every owner unchanged", refused && unchanged, {{"refusal_point", point}});
    }
    {   // live source owner record differs from the accepted digest
      auto t = v7::prepare_v8(o, dt);
      const auto keep = o.accepted_source_states[0];
      o.accepted_source_states[0] = nw::sha("foreign source owner state");
      const bool refused = !v7::try_publish_v8(o, t);
      o.accepted_source_states[0] = keep;
      R.record("T2", "live source owner state differs from the accepted record: refused, owners unchanged", refused && owner_digest(o, scratch) == d0,
               {{"refusal_point", t.refusal_point()}});
    }
    auto t = v7::prepare_v8(o, dt);
    const bool pub = v7::try_publish_v8(o, t);
    R.record("T2", "after every refusal the unmutated step publishes the T1 reference state bit for bit", pub && owner_digest(o, scratch) == reference_after);
  }
  // T3: manufactured faults in preparation and in can_publish refuse with owners unchanged.
  {
    auto o = make();
    const auto d0 = owner_digest(o, scratch);
    int prepare_ok = 0, publish_ok = 0;
    for (int k = 1; k <= v7::prepare_fault_points_v7; ++k) {
      bool threw = false;
      try { (void)v7::prepare_v8(o, dt, k); } catch (const std::exception&) { threw = true; }
      prepare_ok += threw && owner_digest(o, scratch) == d0;
    }
    for (int k = 1; k <= 10; ++k) {
      auto t = v7::prepare_v8(o, dt);
      TA::set_can_publish_fault(t, k);
      const bool refused = !v7::try_publish_v8(o, t);
      publish_ok += refused && t.refusal_point() == 100 + k && owner_digest(o, scratch) == d0;
    }
    R.record("T3", "every preparation fault point refuses with every owner unchanged", prepare_ok == v7::prepare_fault_points_v7, {{"passed_points", prepare_ok}});
    R.record("T3", "every can_publish fault point refuses with every owner unchanged", publish_ok == 10, {{"passed_points", publish_ok}});
    auto t = v7::prepare_v8(o, dt);
    R.record("T3", "after the fault sweep a fresh step reproduces the T1 reference", v7::try_publish_v8(o, t) && owner_digest(o, scratch) == reference_after);
  }
  // T4: single use, sibling outer reservation, moved-from trials.
  {
    auto o = make();
    const auto d0 = owner_digest(o, scratch);
    auto a = v7::prepare_v8(o, dt);
    std::string sibling;
    try { (void)v7::prepare_v8(o, dt); } catch (const std::exception& e) { sibling = e.what(); }
    R.record("T4", "a second concurrent step preparation (sibling outer reservation) refuses, owners unchanged",
             !sibling.empty() && owner_digest(o, scratch) == d0, {{"refusal", sibling}});
    auto moved = std::move(a);
    const bool from_moved = v7::try_publish_v8(o, a);
    R.record("T4", "a moved-from trial refuses", !from_moved && owner_digest(o, scratch) == d0);
    const bool first = v7::try_publish_v8(o, moved);
    const auto d1 = owner_digest(o, scratch);
    const bool again = v7::try_publish_v8(o, moved);
    R.record("T4", "the trial publishes once; a second publication refuses with owners unchanged", first && !again && owner_digest(o, scratch) == d1);
    R.record("T4", "the published state is the T1 reference", d1 == reference_after);
  }
  // T6: genesis admission.
  {
    auto refuses = [&](const std::string& what, std::function<void(v7::GenesisRequestV8&)> edit, v7::GenesisV8<HO>::PacketHook ph = nullptr,
                       v7::GenesisV8<HO>::LedgerHook lh = nullptr) {
      auto rq = request(in, qual_packets());
      edit(rq);
      int reached = 0;
      std::string msg;
      try { (void)v7::GenesisTestAccessV7::admit<HO>(rq, ph, lh, &reached); } catch (const std::exception& e) { msg = e.what(); }
      R.record("T6", what, !msg.empty(), {{"refusal", msg.substr(0, 400)}, {"stage", reached}});
    };
    refuses("a lone electron (no compensating charge) refuses (absolute Gauss)", [](v7::GenesisRequestV8& rq) {
      rq.packets.push_back({"lone", "electron", "qual/lone", {1, 0}, {.0045, .0045, .002}, {0, 0, 0}});
    });
    refuses("genesis time other than the first C19 knot refuses", [](v7::GenesisRequestV8& rq) { rq.seed.time_s = std::nextafter(rq.seed.time_s, 1.); });
    refuses("a declared gross charge that differs from the census refuses", [](v7::GenesisRequestV8&) {}, nullptr,
            [](hc::ChargeLedger& q) { q.particles.gross_C[q.particles.gross_C.size() / 2] += 1e-30; });
    refuses("a missing foil-slice source owner refuses", [](v7::GenesisRequestV8& rq) { rq.sources.slices.pop_back(); });
    refuses("a source patch not declared on the Hall thermal owner refuses", [](v7::GenesisRequestV8& rq) {
      rq.sources.slices[0].patches[0].initial_particles *= 2;
      rq.sources.slices[0].patches[0].id += "-altered";
    });
    refuses("an authored packet in the reserved arrival namespace refuses", [](v7::GenesisRequestV8&) {},
            [](std::vector<m::FullTimeAuthoredPacket>& p) { p.push_back({"arrival-v1/s0/p1/k0", "H2O", "qual/foreign", {1, 0}, {.0045, .0045, 0}, {0, 0, 0}}); });
  }
  // Admission spot checks of the birth path and of the owners' own guards.
  {
    auto o = make();
    const auto d0 = owner_digest(o, scratch);
    auto refusal = [](auto&& f) { try { f(); } catch (const std::exception& e) { return std::string(e.what()); } return std::string(); };
    // source v6 quantized release: a neutral leaf of another weight than w0 is refused
    const auto msg = refusal([&] {
      auto& src = *o.sources[0];
      const auto e = src.make_epoch(dt, "spot-check");
      eh::Epoch ep = e;
      (void)ep;
      ::electrode_heat_v3::TrialFrame f;
      f.role = ::electrode_heat_v3::FieldRole::TotalMidpoint; f.accepted_revision = o.hall->revision(); f.trial_stage = 1; f.mesh_revision = 1;
      f.material_revision = 2; f.begin_s = o.hall->time_s(); f.dt_s = dt; f.end_s = f.begin_s + dt; f.field_stage_s = f.begin_s + .5 * dt;
      f.arena_hash = "spot"; f.field_hash = "spot"; f.mesh_hash = "spot"; f.material_hash = "spot"; f.species_hash = "spot";
      auto fld = std::make_shared<const s7::Field>(o.hall->mesh(), f, o.hall->field().E, o.hall->field().B);
      const s7::FinalFieldView view{fld->frame(), o.hall->mesh(), o.hall->field().E, o.hall->field().B};
      std::vector<std::vector<s7::NeutralLeaf>> leaves(o.slices[0].patches.size());
      leaves[0].push_back({1.5 * o.settings.w0_neutral, eh::Draw{.5, .5, .5, .5, .5, .5}});
      (void)src.prepare(e, fld, view, {}, {}, leaves);
    });
    R.record("spot", "source v6 refuses a neutral leaf whose weight is not w0", msg.find("exactly w0") != std::string::npos, {{"refusal", msg}});
    const auto plan = refusal([&] { pa::StepInput x; x.step = o.arrival->next_step() + 1; (void)o.arrival->plan(x); });
    R.record("spot", "arrival v1.1 refuses a plan for another than the next accepted step", !plan.empty(), {{"refusal", plan}});
    const auto vac = refusal([&] { (void)o.hall->prepare_vacuum(dt); });
    R.record("spot", "Hall refuses a field-only (vacuum) step while carriers are present", !vac.empty(), {{"refusal", vac}});
    const auto big = refusal([&] { (void)v7::prepare_v8(o, 2e-14); });
    R.record("spot", "the coupled step refuses dt above the unchanged 10 fs guard", !big.empty(), {{"refusal", big}});
    R.record("spot", "every spot-check refusal left every owner unchanged", owner_digest(o, scratch) == d0);
    auto t = v7::prepare_v8(o, dt);
    R.record("spot", "after the spot checks a fresh step reproduces the T1 reference", v7::try_publish_v8(o, t) && owner_digest(o, scratch) == reference_after);
  }
  // V7: variable dt, predicate v7, step-size refusals and retries, controller determinism.
  {
    auto o = make();
    const auto d0 = owner_digest(o, scratch);
    // Predicate v7 against the native v6 (V3) predicate on the same observations.
    int compared = 0, equal = 0;
    {
      auto t = v7::prepare_v8(o, dt);
      TA::compare_v6_predicate(t);
      const bool ok = v7::can_publish_v8(o, t);
      const auto [p6, p7] = TA::predicates(t);
      compared += p6 >= 0; equal += p6 >= 0 && p6 == p7;
      R.record("V7", "predicate v7 admits the 10 fs step and agrees with the V3/native-v6 predicate", ok && p6 == 1 && p7 == 1, {{"v6", p6}, {"v7", p7}});
    }
    {
      auto t = v7::prepare_v8(o, dt);
      TA::compare_v6_predicate(t);
      TA::source_after(t) = nw::sha("foreign source seal");
      const bool refused = !v7::try_publish_v8(o, t);
      const auto [p6, p7] = TA::predicates(t);
      compared += p6 >= 0; equal += p6 >= 0 && p6 == p7;
      R.record("V7", "foreign source seal: both predicates refuse, owners unchanged", refused && p6 == 0 && p7 == 0 && owner_digest(o, scratch) == d0,
               {{"refusal_point", t.refusal_point()}, {"v6", p6}, {"v7", p7}});
    }
    {
      auto t = v7::prepare_v8(o, dt);
      TA::compare_v6_predicate(t);
      auto& c = TA::clock(t);
      c.decided_dt_s = std::nextafter(c.decided_dt_s, 0.);
      const bool refused = !v7::try_publish_v8(o, t);
      const auto [p6, p7] = TA::predicates(t);
      R.record("V7", "a step whose dt differs from the controller decision by one ulp is refused by predicate v7, owners unchanged",
               refused && t.refusal_point() == 7 && p7 == 0 && p6 == 1 && owner_digest(o, scratch) == d0, {{"refusal_point", t.refusal_point()}, {"v6", p6}, {"v7", p7}});
    }
    {
      auto t = v7::prepare_v8(o, dt);
      auto& c = TA::clock(t);
      c.cap_s = 2e-14;
      const bool refused = !v7::try_publish_v8(o, t);
      R.record("V7", "a step clock with a cap above the field owner maximum is refused by predicate v7, owners unchanged",
               refused && t.refusal_point() == 7 && owner_digest(o, scratch) == d0, {{"refusal_point", t.refusal_point()}});
    }
    {
      auto t = v7::prepare_v8(o, dt);
      TA::clock(t).decision_sha256 = "not-a-digest";
      const bool refused = !v7::try_publish_v8(o, t);
      R.record("V7", "a step clock without a decision digest is refused by predicate v7, owners unchanged",
               refused && t.refusal_point() == 7 && owner_digest(o, scratch) == d0, {{"refusal_point", t.refusal_point()}});
    }
    auto refusal = [](auto&& f) { try { f(); } catch (const std::exception& e) { return std::string(e.what()); } return std::string(); };
    {
      o.dt_cap_s = 5e-15;
      const auto big = refusal([&] { (void)v7::prepare_v8(o, dt); });
      o.dt_cap_s = 2e-14;
      const auto over = refusal([&] { (void)v7::prepare_v8(o, 5e-15); });
      o.dt_cap_s = ctl.cap_s;
      R.record("V7", "prepare refuses dt above a configured cap of 5 fs and a configured cap above the field owner maximum, owners unchanged",
               !big.empty() && !over.empty() && owner_digest(o, scratch) == d0, {{"refusal_dt", big}, {"refusal_cap", over}});
    }
    // Variable steps publish: ladder rungs 3 and 8 (5 fs) and the predicates agree.
    {
      bool ok = true;
      for (int j : {3, 8}) {
        auto t = v7::prepare_v8(o, v7::ladder_dt(ctl, j));
        TA::compare_v6_predicate(t);
        ok = ok && v7::try_publish_v8(o, t);
        const auto [p6, p7] = TA::predicates(t);
        compared += p6 >= 0; equal += p6 >= 0 && p6 == p7;
        ok = ok && p6 == 1 && p7 == 1;
      }
      R.record("V7", "steps of 9.17 fs (rung 3) and 5 fs (rung 8) publish under predicate v7 (V3 predicate agrees)", ok,
               {{"time_s", o.hall->time_s()}, {"predicate_evaluations_compared", compared}, {"equal", equal}});
    }
    const bool ladder = std::bit_cast<std::uint64_t>(v7::rung_dt(1e-14, 0)) == std::bit_cast<std::uint64_t>(1e-14) &&
                        v7::rung_dt(1e-14, 8) == .5e-14 && v7::rung_dt(1e-14, 16) == .25e-14 && [] {
                          for (int j = -64; j < 64; ++j) if (!(v7::rung_dt(1e-14, j + 1) < v7::rung_dt(1e-14, j))) return false;
                          v7::ControllerConfig c;
                          c.cap_s = 1e-13;
                          const int top = v7::top_rung(c);
                          return v7::ladder_dt(c, top) <= 1e-13 && v7::ladder_dt(c, top - 1) > 1e-13 && v7::ladder_dt(c, 0) == 1e-14 &&
                                 v7::rung_dt(1e-14, -8) == 2e-14; }();
    R.record("V7", "dt ladder: rung 0 is the particle cap bit for bit, rung 8 half of it, rung -8 twice it, strictly decreasing; "
             "with a 100 fs cap the top rung is the largest rung below the cap", ladder);
  }
  // Step-size refusal inside the transaction and retry at a smaller dt.
  {
    auto a = make();
    auto b = make();
    const auto d0 = owner_digest(a, scratch);
    std::string msg;
    try { (void)v7::prepare_v8(a, dt, 108); } catch (const std::exception& e) { msg = e.what(); }
    const bool unchanged = owner_digest(a, scratch) == d0;
    const double retry = v7::ladder_dt(ctl, ctl.retry_rungs);
    auto ta = v7::prepare_v8(a, retry);
    const bool pa_ok = v7::try_publish_v8(a, ta);
    auto tb = v7::prepare_v8(b, retry);
    const bool pb_ok = v7::try_publish_v8(b, tb);
    R.record("V7", "a step-size refusal after the Hall proposal is classified as a step-size guard and leaves every owner unchanged",
             v7::step_size_refusal(msg) && unchanged, {{"refusal", msg}});
    R.record("V7", "the retried step (5 fs) equals a direct 5 fs step from the same state bit for bit", pa_ok && pb_ok && owner_digest(a, scratch) == owner_digest(b, scratch));
    R.record("V7", "a non-step-size fault is not retried", !v7::step_size_refusal("manufactured v6 preparation fault") &&
             !v7::step_size_refusal("absolute Gauss constraint disagrees with actual declared charge"));
  }
  // Controller inputs: partition/order independence, repeatability, Hall plasma operand.
  {
    const auto dense = load_inputs("dense");
    auto o = genesis<HO>(dense, "dense");
    for (int k = 0; k < 3; ++k) { auto t = v7::prepare_v8(o, dt); if (!v7::try_publish_v8(o, t)) break; }
    v7::ControllerState cst;
    v7::CellBounds cb;
    const auto base = v7::control_inputs(o, cst, ctl, cb, 1);
    auto adaptive = ctl;
    adaptive.adaptive = true;
    const auto d1 = v7::decide(adaptive, base, o.hall->model().drive(), true);
    const auto j1 = v7::decision_json(adaptive, base, d1).dump();
    bool same = true;
    J ev = J::array();
    for (int parts : {1, 2, 3, 8}) {
      const auto x = v7::control_inputs(o, cst, ctl, cb, parts);
      const auto d = v7::decide(adaptive, x, o.hall->model().drive(), true);
      const bool eq = v7::decision_json(adaptive, x, d).dump() == j1 &&
                      std::bit_cast<std::uint64_t>(x.omega_p_max) == std::bit_cast<std::uint64_t>(base.omega_p_max);
      same = same && eq;
      ev.push_back({{"parts", parts}, {"omega_p_max", hexd(x.omega_p_max)}, {"dt_s", hexd(d.dt_s)}, {"identical", eq}});
    }
    R.record("V7", "controller inputs and decision are bit-identical for the census in 1, 2, 3 and 8 parts in reversed order", same,
             {{"partitions", ev}, {"charged", base.charged}, {"criterion", v7::criterion_name(d1.criterion)}});
    auto t = v7::prepare_v8(o, d1.dt_s);
    const double measured = t.proposal().motion().budget().plasma_dt;
    const double predicted = base.omega_p_max * d1.dt_s;
    R.record("V7", "the controller's omega_p dt equals the Hall motion plasma guard value of the step (summation order only)",
             std::abs(measured - predicted) <= 1e-12 * measured && measured > 0,
             {{"predicted", predicted}, {"measured", measured}, {"dt_s", d1.dt_s}, {"criterion", v7::criterion_name(d1.criterion)},
              {"criterion_value", d1.value}});
  }
  // A real Hall plasma-guard refusal: a dense blob 8x heavier refuses 10 fs; the
  // controller's dt for that state publishes.
  {
    const double w = opt.blob_weight;
    opt.blob_weight = 8 * w;
    const auto dense = load_inputs("dense");
    auto o = genesis<HO>(dense, "dense");
    opt.blob_weight = w;
    const auto d0 = owner_digest(o, scratch);
    std::string msg;
    try { (void)v7::prepare_v8(o, dt); } catch (const std::exception& e) { msg = e.what(); }
    const bool unchanged = owner_digest(o, scratch) == d0;
    v7::ControllerState cst;
    v7::CellBounds cb;
    auto adaptive = ctl;
    adaptive.adaptive = true;
    const auto in = v7::control_inputs(o, cst, adaptive, cb);
    const auto d = v7::decide(adaptive, in, o.hall->model().drive(), true);
    auto t = v7::prepare_v8(o, d.dt_s);
    const bool pub = v7::try_publish_v8(o, t);
    R.record("V7", "Hall refuses 10 fs on a dense blob (plasma guard, step-size class), owners unchanged; the controller's dt publishes",
             !msg.empty() && v7::step_size_refusal(msg) && unchanged && pub,
             {{"refusal", msg}, {"omega_p_dt_at_10fs", in.omega_p_max * 1e-14}, {"controller_dt_s", d.dt_s}, {"criterion", v7::criterion_name(d.criterion)},
              {"criterion_value", d.value}});
  }
  // V7.1: foil planes of the generic rule equal the coordinator v6 planes on the 200 um mesh (bit for bit).
  {
    const axial::Mesh g(45, 4, .009, 1);
    const auto a = v7::hall_geometry(g), b = v7::hall_geometry_v6_200um(g);
    bool same = true;
    J ev = J::array();
    for (const auto* f : {&a.lower_foil, &a.upper_foil})
      for (double x : {f->left_m, f->right_m, f->bottom_m, f->top_m}) ev.push_back(hexd(x));
    for (int k = 0; k < 2; ++k) {
      const auto& fa = k ? a.upper_foil : a.lower_foil;
      const auto& fb = k ? b.upper_foil : b.lower_foil;
      for (auto [x, y] : {std::pair{fa.left_m, fb.left_m}, {fa.right_m, fb.right_m}, {fa.bottom_m, fb.bottom_m}, {fa.top_m, fb.top_m}})
        same = same && std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y);
    }
    R.record("V7.1", "generic foil-plane rule gives the coordinator v6 planes on the 200 um mesh bit for bit", same, {{"planes", ev}});
    R.record("V7.1", std::string("binding pair ") + v7::bind::pair_name, true, {{"hall", v7::bind::hall_binding}, {"native", v7::bind::native_binding}});
    // Exact weight sums of the ledger totals against native's exact counts.
    std::mt19937_64 rng(20261103);
    std::uniform_real_distribution<double> u01(0, 1);
    const double typical[] = {1., 1e6, 2e-23, 2., 1e15, 3.45e10, 1.6e11, 5e-324};
    bool values = true, orders = true;
    for (int t = 0; t < 2000; ++t) {
      v7::ExactWeightSum x, y;
      m::ArenaCount c;
      std::vector<double> ws;
      const int k = 1 + int(rng() % 40);
      for (int i = 0; i < k; ++i)
        ws.push_back(rng() % 2 ? typical[rng() % 8] : std::ldexp(u01(rng), int(rng() % 2000) - 1074 + 50));
      for (double w : ws) { x.add(w); c = c.credited(w); }
      for (auto it = ws.rbegin(); it != ws.rend(); ++it) y.add(*it);
      values = values && std::bit_cast<std::uint64_t>(x.value()) == std::bit_cast<std::uint64_t>(v7::exact_value(c));
      orders = orders && x == y;
    }
    R.record("V7.1", "exact weight sum equals the ArenaCount rounded value on 2000 random sums and is order independent", values && orders,
             {{"values", values}, {"orders", orders}});
  }
  J report = {{"cases", R.cases}, {"failed", R.failed}, {"forbidden_native_calls", forbidden_calls}};
  write_text(out / "checks-result.json", report.dump(1) + "\n");
  fs::remove_all(scratch);
  return R.failed == 0 ? 0 : 1;
}
