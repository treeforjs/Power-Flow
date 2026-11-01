// Coordinator v8 run loop, whole-run checkpoint and restore, Hall-only pre-run (included by driver_v8.cpp).

// ---------------------------------------------------------------- configuration checks
void finish_configuration() {
  v8::RunTargetRequest{opt.end_time_requested, opt.window_requested, opt.steps_requested,
                       opt.end_time_s, opt.window_s, opt.steps}.validate();
  ctl.fixed_dt_s = opt.dt_s;
  ctl.source_interval = opt.source_emission_error > 0 || opt.source_desorption_error > 0;
  ctl.validate(std::min({opt.hall_dt_max_s, opt.source_dt_max_s}));
  if (ctl.neutrals_need_particle_cap) {
    if (ctl.particle_cap_s > opt.neutral_dt_max_s)
      throw std::runtime_error("particle cap above the resident-neutral interval bound (--neutral-dt-max)");
  } else if (opt.neutral_dt_max_s < ctl.cap_s)
    throw std::runtime_error("--charged-only-particle-cap needs a resident-neutral interval bound >= the controller cap (--neutral-dt-max)");
  if (!(opt.source_dt_max_s > 0 && opt.source_dt_max_s <= 1e-12)) throw std::runtime_error("--source-dt-max must lie in (0, 1e-12]");
}

// Run options recorded in every checkpoint and applied at restore (the physics and controller configuration;
// operational options -- steps, end time, outputs, checkpoint cadence, test hooks -- come from the restoring launch).
J options_json() {
  J c = {{"adaptive", ctl.adaptive}, {"fixed_grid_retry", ctl.fixed_grid_retry}, {"fixed_dt_s", b64(ctl.fixed_dt_s)}, {"cap_s", b64(ctl.cap_s)},
         {"particle_cap_s", b64(ctl.particle_cap_s)}, {"neutrals_need_particle_cap", ctl.neutrals_need_particle_cap}, {"c_plasma", b64(ctl.c_plasma)},
         {"c_gyro", b64(ctl.c_gyro)}, {"c_disp", b64(ctl.c_disp)}, {"c_chem", b64(ctl.c_chem)}, {"c_drive", b64(ctl.c_drive)},
         {"drive_floor_A", b64(ctl.drive_floor_A)}, {"lte_tol", b64(ctl.lte_tol)}, {"birth_threshold", b64(ctl.birth_threshold)},
         {"growth_rungs", ctl.growth_rungs}, {"hold_after_retry", ctl.hold_after_retry}, {"retry_rungs", ctl.retry_rungs}, {"max_retries", ctl.max_retries},
         {"grow_margin", b64(ctl.grow_margin)}, {"dt_min_s", b64(ctl.dt_min_s)}, {"sync_every_s", b64(ctl.sync_every_s)},
         {"sync_origin_s", b64(ctl.sync_origin_s)}, {"source_interval", ctl.source_interval}};
  J o = {{"field", opt.field}, {"push", opt.push}, {"fixture", opt.fixture}, {"dt_s", b64(opt.dt_s)}, {"w0e", b64(opt.w0e)}, {"w0n", b64(opt.w0n)},
         {"phi", b64(opt.phi)}, {"faces", opt.faces}, {"gpu_capacity", opt.gpu_capacity}, {"seed_offset", opt.seed_offset},
         {"blob_pairs", opt.blob_pairs}, {"blob_weight", b64(opt.blob_weight)}, {"blob_eV", b64(opt.blob_eV)},
         {"manufactured_refusal_step", opt.manufactured_refusal_step}, {"hall_dt_max_s", b64(opt.hall_dt_max_s)},
         {"probe_energy_every", opt.probe_energy_every}, {"host_work", opt.host_work}, {"mesh_file", opt.mesh_file}, {"collision_grid", opt.collision_grid},
         {"load_carriers", opt.load_carriers}, {"load_neutral_fraction", b64(opt.load_neutral_fraction)}, {"hall_owner", opt.hall_owner},
         {"source_dt_max_s", b64(opt.source_dt_max_s)}, {"source_emission_error", b64(opt.source_emission_error)},
         {"source_desorption_error", b64(opt.source_desorption_error)}, {"source_safety", b64(opt.source_safety)},
         {"source_count_floor", b64(opt.source_count_floor)}, {"source_interval_preset", opt.source_interval_preset},
         {"neutral_dt_max_s", b64(opt.neutral_dt_max_s)}, {"surface", opt.surface}, {"hall_checkpoint", opt.hall_checkpoint.string()}};
  return {{"options", o}, {"controller", c}};
}
void apply_options(const J& j) {
  const auto& o = j.at("options");
  const auto& c = j.at("controller");
  opt.field = o.at("field"); opt.push = o.at("push"); opt.fixture = o.at("fixture"); opt.dt_s = from_b64(o.at("dt_s"));
  opt.w0e = from_b64(o.at("w0e")); opt.w0n = from_b64(o.at("w0n")); opt.phi = from_b64(o.at("phi")); opt.faces = o.at("faces");
  opt.gpu_capacity = o.at("gpu_capacity"); opt.seed_offset = o.at("seed_offset"); opt.blob_pairs = o.at("blob_pairs");
  opt.blob_weight = from_b64(o.at("blob_weight")); opt.blob_eV = from_b64(o.at("blob_eV"));
  opt.manufactured_refusal_step = o.at("manufactured_refusal_step"); opt.hall_dt_max_s = from_b64(o.at("hall_dt_max_s"));
  opt.probe_energy_every = o.at("probe_energy_every"); opt.host_work = o.at("host_work"); opt.mesh_file = o.at("mesh_file");
  opt.collision_grid = o.at("collision_grid"); opt.load_carriers = o.at("load_carriers"); opt.load_neutral_fraction = from_b64(o.at("load_neutral_fraction"));
  opt.hall_owner = o.at("hall_owner"); opt.source_dt_max_s = from_b64(o.at("source_dt_max_s"));
  opt.source_emission_error = from_b64(o.at("source_emission_error")); opt.source_desorption_error = from_b64(o.at("source_desorption_error"));
  opt.source_safety = from_b64(o.at("source_safety")); opt.source_count_floor = from_b64(o.at("source_count_floor"));
  opt.source_interval_preset = o.at("source_interval_preset"); opt.neutral_dt_max_s = from_b64(o.at("neutral_dt_max_s"));
  opt.surface = o.at("surface"); opt.hall_checkpoint = o.at("hall_checkpoint").template get<std::string>();
  ctl.adaptive = c.at("adaptive"); ctl.fixed_grid_retry = c.at("fixed_grid_retry"); ctl.fixed_dt_s = from_b64(c.at("fixed_dt_s"));
  ctl.cap_s = from_b64(c.at("cap_s")); ctl.particle_cap_s = from_b64(c.at("particle_cap_s")); ctl.neutrals_need_particle_cap = c.at("neutrals_need_particle_cap");
  ctl.c_plasma = from_b64(c.at("c_plasma")); ctl.c_gyro = from_b64(c.at("c_gyro")); ctl.c_disp = from_b64(c.at("c_disp"));
  ctl.c_chem = from_b64(c.at("c_chem")); ctl.c_drive = from_b64(c.at("c_drive")); ctl.drive_floor_A = from_b64(c.at("drive_floor_A"));
  ctl.lte_tol = from_b64(c.at("lte_tol")); ctl.birth_threshold = from_b64(c.at("birth_threshold")); ctl.growth_rungs = c.at("growth_rungs");
  ctl.hold_after_retry = c.at("hold_after_retry"); ctl.retry_rungs = c.at("retry_rungs"); ctl.max_retries = c.at("max_retries");
  ctl.grow_margin = from_b64(c.at("grow_margin")); ctl.dt_min_s = from_b64(c.at("dt_min_s")); ctl.sync_every_s = from_b64(c.at("sync_every_s"));
  ctl.sync_origin_s = from_b64(c.at("sync_origin_s")); ctl.source_interval = c.at("source_interval");
}

// ---------------------------------------------------------------- serialization of coordinator state
J count_json(const m::ArenaCount& c) { return b64s(c.positive_carriers()); }
m::ArenaCount count_of(const J& j) {
  m::ArenaCount c;
  for (double w : from_b64s(j)) c = c.credited(w);
  if (c.positive_carriers() != from_b64s(j)) throw std::runtime_error("checkpoint noncanonical exact count");
  return c;
}
J words_json(const hc::Words& w) { return J::array({b64(w.hi), b64(w.lo)}); }
hc::Words words_of(const J& j) { hc::Words w; w.hi = from_b64(j.at(0)); w.lo = from_b64(j.at(1)); return w; }
J totals_json(const Totals& t) {
  return {{"n_e", t.n_e}, {"n_i", t.n_i}, {"n_n", t.n_n}, {"w_e", t.w_e.hex()}, {"w_i", t.w_i.hex()}, {"w_n", t.w_n.hex()},
          {"h_nuclei", t.h_nuclei.hex()}, {"o_nuclei", t.o_nuclei.hex()}, {"stock", count_json(t.stock)}, {"kinetic_J", ld_hex(t.kinetic_J)},
          {"reference_J", ld_hex(t.reference_J)}, {"q_free", words_json(t.q_free)}, {"q_conductor", words_json(t.q_conductor)},
          {"q_exterior", words_json(t.q_exterior)}, {"q_gross", b64(t.q_gross)}};
}
Totals totals_of(const J& j) {
  Totals t;
  t.n_e = j.at("n_e"); t.n_i = j.at("n_i"); t.n_n = j.at("n_n");
  t.w_e = XS::from_hex(j.at("w_e")); t.w_i = XS::from_hex(j.at("w_i")); t.w_n = XS::from_hex(j.at("w_n"));
  t.h_nuclei = XS::from_hex(j.at("h_nuclei")); t.o_nuclei = XS::from_hex(j.at("o_nuclei")); t.stock = count_of(j.at("stock"));
  t.kinetic_J = from_ld_hex(j.at("kinetic_J")); t.reference_J = from_ld_hex(j.at("reference_J"));
  t.q_free = words_of(j.at("q_free")); t.q_conductor = words_of(j.at("q_conductor")); t.q_exterior = words_of(j.at("q_exterior"));
  t.q_gross = from_b64(j.at("q_gross"));
  return t;
}
J run_state_json(const RunState& s) {
  return {{"genesis", totals_json(s.genesis)}, {"water_initial", s.water_initial.hex()}, {"sunk_h2o", count_json(s.sunk_h2o)},
          {"sunk_ion", count_json(s.sunk_ion)}, {"q_initial", words_json(s.q_initial)}, {"reaction_defect_cum", ld_hex(s.reaction_defect_cum)},
          {"worst_step_residual", b64(s.worst_step_residual)}, {"worst_gauss_ratio", b64(s.worst_gauss_ratio)},
          {"worst_divB_ratio", b64(s.worst_divB_ratio)}, {"worst_q_change", b64(s.worst_q_change)},
          {"worst_reaction_defect", b64(s.worst_reaction_defect)}, {"water_exact", s.water_exact}, {"nuclei_exact", s.nuclei_exact},
          {"electron_births", s.electron_births}, {"neutral_births", s.neutral_births}, {"sinks", s.sinks}, {"events", s.events}};
}
RunState run_state_of(const J& j) {
  RunState s;
  s.genesis = totals_of(j.at("genesis")); s.water_initial = XS::from_hex(j.at("water_initial"));
  s.sunk_h2o = count_of(j.at("sunk_h2o")); s.sunk_ion = count_of(j.at("sunk_ion")); s.q_initial = words_of(j.at("q_initial"));
  s.reaction_defect_cum = from_ld_hex(j.at("reaction_defect_cum"));
  s.worst_step_residual = from_b64(j.at("worst_step_residual")); s.worst_gauss_ratio = from_b64(j.at("worst_gauss_ratio"));
  s.worst_divB_ratio = from_b64(j.at("worst_divB_ratio")); s.worst_q_change = from_b64(j.at("worst_q_change"));
  s.worst_reaction_defect = from_b64(j.at("worst_reaction_defect")); s.water_exact = j.at("water_exact"); s.nuclei_exact = j.at("nuclei_exact");
  s.electron_births = j.at("electron_births"); s.neutral_births = j.at("neutral_births"); s.sinks = j.at("sinks"); s.events = j.at("events");
  return s;
}
J rate_json(const std::vector<sv7::RateSample>& v) {
  J a = J::array();
  for (const auto& r : v)
    a.push_back({{"t_s", b64(r.t_s)}, {"dt_s", b64(r.dt_s)}, {"log_emission", b64s(r.log_emission)}, {"log_desorption", b64s(r.log_desorption)},
                 {"T_K", b64s(r.T_K)}});
  return a;
}
std::vector<sv7::RateSample> rates_of(const J& a) {
  std::vector<sv7::RateSample> v;
  for (const auto& r : a) {
    sv7::RateSample x;
    x.t_s = from_b64(r.at("t_s")); x.dt_s = from_b64(r.at("dt_s"));
    x.log_emission = from_b64s(r.at("log_emission")); x.log_desorption = from_b64s(r.at("log_desorption")); x.T_K = from_b64s(r.at("T_K"));
    v.push_back(std::move(x));
  }
  return v;
}
std::string ring_bytes(const v8::ControllerState& c) {
  std::string b;
  auto u64 = [&](std::uint64_t x) { for (int i = 0; i < 8; ++i) b.push_back(char((x >> (8 * i)) & 255)); };
  auto f64 = [&](double x) { u64(std::bit_cast<std::uint64_t>(x)); };
  b += "coordinator-v8-controller-ring-v1";
  u64(c.ring.size());
  for (const auto& s : c.ring) {
    f64(s.t);
    u64(s.E.size());
    for (double x : s.E) f64(x);
    u64(s.B.size());
    for (double x : s.B) f64(x);
  }
  return b;
}
void ring_of(const std::string& b, v8::ControllerState& c, std::size_t ne, std::size_t nb) {
  std::size_t at = 0;
  const std::string tag = "coordinator-v8-controller-ring-v1";
  if (b.compare(0, tag.size(), tag) != 0) throw std::runtime_error("controller ring tag");
  at = tag.size();
  auto u64 = [&] {
    if (b.size() - at < 8) throw std::runtime_error("controller ring truncated");
    std::uint64_t x = 0;
    for (int i = 0; i < 8; ++i) x |= std::uint64_t(static_cast<unsigned char>(b[at++])) << (8 * i);
    return x;
  };
  auto f64 = [&] { return std::bit_cast<double>(u64()); };
  const auto n = u64();
  if (n > 4) throw std::runtime_error("controller ring size");
  c.ring.clear();
  for (std::uint64_t k = 0; k < n; ++k) {
    v8::ControllerState::Snap s;
    s.t = f64();
    if (u64() != ne) throw std::runtime_error("controller ring E dimension");
    s.E.resize(ne);
    for (auto& x : s.E) x = f64();
    if (u64() != nb) throw std::runtime_error("controller ring B dimension");
    s.B.resize(nb);
    for (auto& x : s.B) x = f64();
    c.ring.push_back(std::move(s));
  }
  if (at != b.size()) throw std::runtime_error("controller ring trailing bytes");
}

// ---------------------------------------------------------------- the run (stepper)
struct Outputs {
  fs::path dir;
  std::ofstream ledger_bin, digests, det_bin, sub_bin, coord_bin, control_bin, receipts, source_interval;
  FrameWriter frames;
  ProbeWriter probes;
};
// Crash injection (tests): exit 86 while writing the checkpoint of opt.crash_step at opt.crash_stage.
void crash_hook(int stage) {
  if (opt.crash_stage == stage) {
    std::cout << "crash injected in checkpoint write at stage " << stage << std::endl;
    std::_Exit(86);
  }
}
template <class H> struct Run {
  const Inputs& in;
  v8::OwnersV8<H> o;
  RunState st;
  Totals before;
  std::vector<m::ArenaRateEntry> population;
  v8::ControllerState cst;
  std::string receipt_head;
  double end_time = 0;
  long step_target = 20;
  long step = 0;   // last accepted step
  bool legacy = false;
  std::map<std::string, long> criteria_histogram;
  long refused_attempts_total = 0, parity_compared = 0, parity_identical = 0, max_swap_alloc = 0, steps_this_process = 0;
  double dt_min_used = INFINITY, dt_max_used = 0;
  std::vector<double> step_ms;
  std::vector<std::vector<double>> sub_rows;
  Ledger L, C;
  CoordTimes ct;
  DeviceCounters dc_prev{};
  std::shared_ptr<hc::MotionProvider> device_motion;
  Outputs* out = nullptr;
  std::string genesis_manifest_sha256, last_manifest_sha256;
  fs::path checkpoint_root;

  Run(const Inputs& i, v8::OwnersV8<H>&& owners) : in(i), o(std::move(owners)) { device_motion = o.motion; legacy = opt.host_work == "v7"; o.legacy_host_work = legacy; }

  // Fresh start (genesis): the v7.1 order of the initial records.
  void start_fresh(J& summary) {
    population = v8::accepted_population(o);
    st.genesis = legacy ? totals_legacy(o) : totals(o, population);
    st.water_initial = st.genesis.w_n;
    st.water_initial.add(st.genesis.w_i);
    st.water_initial.add(st.genesis.stock);
    st.q_initial = st.genesis.q_free + st.genesis.q_conductor + st.genesis.q_exterior;
    summary["initial"] = {{"source_stock_molecules", v8::exact_value(st.genesis.stock)}, {"native_h2o", st.genesis.w_n.value()},
                          {"native_ion", st.genesis.w_i.value()}, {"native_electron", st.genesis.w_e.value()},
                          {"total_charge_C", st.q_initial.rounded()}};
    before = st.genesis;
    end_time = v8::RunTargetRequest{opt.end_time_requested, opt.window_requested, opt.steps_requested,
                                    opt.end_time_s, opt.window_s, opt.steps}.after_restore(0, o.hall->time_s());
    step_target = opt.steps;
    { const auto a = o.accepted(); cst.record_state(o.hall->time_s(), *a.E, *a.B); }
    receipt_head = nw::sha("coordinator-v7-step-receipts-genesis/" + ctl.sha256() + "/" + o.arena->state_prefix() + "/" + o.hall->field_sha256() + "/" +
                           o.accepted_source_seal);
    if (ctl.sync_every_s > 0 || ctl.fixed_grid_retry) ctl.sync_origin_s = o.hall->time_s();
    v8::sample_source_rates(o);
  }
  void open_outputs(Outputs& x, bool fresh) {
    out = &x;
    x.ledger_bin.open(x.dir / "ledger.bin", std::ios::binary);
    x.digests.open(x.dir / "step-digests.txt");
    x.det_bin.open(x.dir / "deterministic.bin", std::ios::binary);
    x.sub_bin.open(x.dir / "subphase-timings.bin", std::ios::binary);
    x.coord_bin.open(x.dir / "coordinator-timings.bin", std::ios::binary);
    write_text(x.dir / "coordinator-timings-schema.json", J{{"schema", "coordinator-v7.1-coordinator-timings-v1"}, {"byte_order", "little"},
                                                            {"value", "float64"}, {"columns", coordinator_columns}}.dump(1) + "\n");
    write_text(x.dir / "subphase-schema.json", J{{"schema", "coordinator-v6.1-subphase-timings-v1"}, {"byte_order", "little"}, {"value", "float64"},
                                                 {"columns", subphase_columns}}.dump(1) + "\n");
    if (opt.frame_every > 0) {
      x.frames.open(x.dir, o.hall->mesh());
      write_text(x.dir / "frames-layout.json", x.frames.layout.dump(1) + "\n");
      if (fresh) x.frames.write(o, 0, population);
    }
    x.control_bin.open(x.dir / "control-ledger.bin", std::ios::binary);
    x.receipts.open(x.dir / "step-receipts.jsonl");
    x.probes.open(x.dir, o);
    if (fresh) x.probes.write(o, 0);
    if (ctl.source_interval) x.source_interval.open(x.dir / "source-interval.jsonl");
  }
  void flush() {
    if (!out) return;
    out->ledger_bin.flush(); out->digests.flush(); out->det_bin.flush(); out->sub_bin.flush(); out->control_bin.flush(); out->receipts.flush();
    out->probes.flush(); out->coord_bin.flush();
    if (out->source_interval.is_open()) out->source_interval.flush();
  }
  bool finished(long next) const {
    if (end_time > 0) return o.hall->time_s() >= end_time - 1e-6 * std::min(ctl.cap_s, ctl.fixed_dt_s);
    return next > step_target;
  }
  // One accepted step (the v7.1 loop body). Returns "" or the stop reason.
  std::string advance() {
    const long s = step + 1;
    std::shared_ptr<ParityMotionProvider> parity;
    if (opt.push == "gpu" && opt.parity_every > 0 && (s == 1 || s % opt.parity_every == 0)) {
      parity = std::make_shared<ParityMotionProvider>(device_motion);
      o.motion = parity;
    } else o.motion = device_motion;
    const auto tc = steady::now();
    v8::CellBounds cb;
    if (legacy) cb.legacy_27_loop = true;
    const auto cin = legacy ? v8::control_inputs(o, cst, ctl, cb) : v8::control_inputs_from(o, population, cst, ctl, cb);
    const auto dec = v8::decide(ctl, cin, o.hall->model().drive(), s == 1);
    const double control_s = since(tc);
    if (dec.refused)
      return std::string("controller found no admissible dt >= dt_min at step ") + std::to_string(s) + " (" + v8::criterion_name(dec.criterion) + ")";
    auto attempt_dec = dec;
    J attempts = J::array();
    std::string decision_sha;
    std::optional<v8::TrialV8<H>> trial;
    double prepare_s = 0;
    std::string refused;
    for (int attempt = 0;; ++attempt) {
      const J record = v8::decision_json(ctl, cin, attempt_dec);
      decision_sha = nw::sha(record.dump());
      const auto t0 = steady::now();
      try {
        if (opt.manufactured_refusal_step == s && attempt == 0)
          throw std::runtime_error("manufactured step-size refusal (test hook): motion full-census plasma guard");
        trial.emplace(v8::prepare_v8(o, attempt_dec.dt_s, 0, &decision_sha));
        prepare_s = since(t0);
        attempts.push_back({{"attempt", attempt}, {"dt_s", hexd(attempt_dec.dt_s)}, {"rung", attempt_dec.rung}, {"decision_sha256", decision_sha},
                            {"outcome", "prepared"}});
        break;
      } catch (const std::exception& e) {
        attempts.push_back({{"attempt", attempt}, {"dt_s", hexd(attempt_dec.dt_s)}, {"rung", attempt_dec.rung}, {"decision_sha256", decision_sha},
                            {"outcome", "refused"}, {"refusal", e.what()}, {"step_size_guard", v8::step_size_refusal(e.what())}});
        ++refused_attempts_total;
        const int rung = attempt_dec.rung + ctl.retry_rungs;
        const double next = ctl.adaptive ? v8::ladder_dt(ctl, rung) : .5 * attempt_dec.dt_s;
        if (!(ctl.adaptive || ctl.fixed_grid_retry) || !v8::step_size_refusal(e.what()) || attempt >= ctl.max_retries || next < ctl.dt_min_s) {
          refused = std::string("prepare refused at step ") + std::to_string(s) + ": " + e.what();
          break;
        }
        attempt_dec.rung = ctl.adaptive ? rung : attempt_dec.rung;
        attempt_dec.dt_s = next;
        attempt_dec.landing = false;
        attempt_dec.criterion = v8::Criterion::retry;
        attempt_dec.value = attempt + 1;
        attempt_dec.limit = ctl.max_retries;
      }
    }
    if (!refused.empty()) {
      if (out) out->receipts << J{{"schema", "coordinator-v7-step-refusal-v1"}, {"step", s}, {"previous_sha256", receipt_head}, {"attempts", attempts}}.dump() << '\n';
      return refused;
    }
    const double t_begin = o.hall->time_s();
    const J pre = pre_numbers(*trial);
    const auto interval_records = trial->interval_records();
    const bool interval_startup = trial->interval_startup();
    const auto pub = publish_measured(o, *trial);
    if (!pub.published) {
      if (out)
        out->receipts << J{{"schema", "coordinator-v7-step-refusal-v1"}, {"step", s}, {"previous_sha256", receipt_head}, {"attempts", attempts},
                           {"publication_refusal_point", pub.refusal}}.dump() << '\n';
      return "publication refused at step " + std::to_string(s) + " (refusal point " + std::to_string(pub.refusal) + ")";
    }
    step = s;
    ++steps_this_process;
    max_swap_alloc = std::max(max_swap_alloc, long(pub.swap_allocations));
    auto tq = steady::now();
    if (!legacy) population = v8::accepted_population(o);
    else if (opt.frame_every > 0 && s % opt.frame_every == 0) population = o.arena->rate_population(o.inputs.maximum_carriers);
    const auto after = legacy ? totals_legacy(o) : totals(o, population);
    ct.totals_s = since(tq); tq = steady::now();
    fill_row(L, o, *trial, pre, pub, before, after, st, s, prepare_s);
    ct.ledger_row_s = since(tq); tq = steady::now();
    if (out && L.first)
      write_text(out->dir / "ledger-schema.json", J{{"schema", "coordinator-v6-step-ledger-v1"}, {"byte_order", "little"}, {"value", "float64"},
                                                    {"columns", L.names}}.dump(1) + "\n");
    if (out) {
      out->ledger_bin.write(reinterpret_cast<const char*>(L.row.data()), std::streamsize(L.row.size() * sizeof(double)));
      for (std::size_t k = 0; k < L.row.size(); ++k)
        if (deterministic_column(L.names[k])) out->det_bin.write(reinterpret_cast<const char*>(&L.row[k]), sizeof(double));
    }
    L.first = false;
    ct.ledger_io_s = since(tq); tq = steady::now();
    const std::string rng = o.chemistry->accepted_rng_sha256();
    if (out) out->digests << digest_line(o, s, rng);
    ct.digest_line_s = since(tq); tq = steady::now();
    {
      const int retries = int(attempts.size()) - 1;
      cst.dt_prev_s = attempt_dec.dt_s;
      cst.rung_prev = attempt_dec.rung;
      cst.hold = retries > 0 ? ctl.hold_after_retry : std::max(0, cst.hold - 1);
      cst.chem_rate = chemistry_rate(trial->chemistry_state().coverage());
      cst.births_rate = trial->numbers().expected_electrons / o.settings.w0_electron / attempt_dec.dt_s;
      { const auto a = o.accepted(); cst.record_state(o.hall->time_s(), *a.E, *a.B); }
      v8::sample_source_rates(o);
      ct.controller_state_s = since(tq); tq = steady::now();
      control_row(C, s, t_begin, o.hall->time_s(), attempt_dec, dec, cin, int(attempts.size()), pre, control_s);
      if (out && C.first)
        write_text(out->dir / "control-schema.json", J{{"schema", "coordinator-v7-control-ledger-v1"}, {"byte_order", "little"}, {"value", "float64"},
                                                       {"columns", C.names}, {"criterion_codes", criterion_codes()}}.dump(1) + "\n");
      if (out) out->control_bin.write(reinterpret_cast<const char*>(C.row.data()), std::streamsize(C.row.size() * sizeof(double)));
      C.first = false;
      ct.control_row_s = since(tq); tq = steady::now();
      if (out) out->probes.write(o, s);
      ct.probes_s = since(tq); tq = steady::now();
      J receipt = {{"schema", "coordinator-v7-step-receipt-v1"}, {"step", s}, {"previous_sha256", receipt_head},
                   {"controller_sha256", ctl.sha256()}, {"decision", v8::decision_json(ctl, cin, attempt_dec)}, {"decision_sha256", decision_sha},
                   {"attempts", attempts}, {"common_clock_sha256", trial->common_clock_sha256()}, {"dt_s", hexd(attempt_dec.dt_s)},
                   {"dt_bits", bits(attempt_dec.dt_s)}, {"criterion", v8::criterion_name(attempt_dec.criterion)},
                   {"criterion_value", attempt_dec.value}, {"criterion_limit", std::isfinite(attempt_dec.limit) ? J(attempt_dec.limit) : J("inf")},
                   {"after", {{"time_s", hexd(o.hall->time_s())}, {"native_prefix", o.arena->state_prefix()}, {"field", o.hall->field_sha256()},
                              {"heat", o.hall->heat_sha256()}, {"kinetic_ledgers", o.hall->kinetic_ledgers_sha256()},
                              {"source_seal", o.accepted_source_seal}, {"rng", legacy ? o.chemistry->accepted_rng_sha256() : rng},
                              {"arrival", o.accepted_arrival_state}}}};
      if (ctl.source_interval) receipt["source_interval"] = {{"startup", interval_startup}, {"owner_records", interval_records}};
      receipt_head = nw::sha(receipt_head + receipt.dump());
      receipt["receipt_sha256"] = receipt_head;
      if (out) out->receipts << receipt.dump() << '\n';
      if (out && ctl.source_interval)
        out->source_interval << J{{"step", s}, {"dt_s", hexd(attempt_dec.dt_s)}, {"bound_s", hexd(cin.source_bound_s)}, {"startup", cin.source_startup},
                                  {"limiting", cin.source_limiting}, {"rate_per_s", cin.source_rate_per_s}, {"owner_records", interval_records}}.dump() << '\n';
      ct.receipt_s = since(tq); tq = steady::now();
      ++criteria_histogram[v8::criterion_name(attempt_dec.criterion)];
      dt_min_used = std::min(dt_min_used, attempt_dec.dt_s);
      dt_max_used = std::max(dt_max_used, attempt_dec.dt_s);
    }
    {
      const auto row = subphase_row(*trial, s);
      if (out) out->sub_bin.write(reinterpret_cast<const char*>(row.data()), std::streamsize(row.size() * sizeof(double)));
      sub_rows.push_back(row);
    }
    if (parity) { parity_compared += parity->compared; parity_identical += parity->identical; }
    before = after;
    step_ms.push_back(1e3 * (prepare_s + pub.total_s));
    tq = steady::now();
    if (out && opt.frame_every > 0 && s % opt.frame_every == 0) out->frames.write(o, s, population);
    ct.frames_s = since(tq);
    ct.control_s = control_s;
    ct.observe_prepare_s = trial->times().observe_s;
    ct.observe_preflight_s = trial->sub_times().cp_observe_s;
    ct.step_s = (prepare_s + pub.total_s);
    {
      const auto dc = device_counters(in, device_motion);
      ct.gpu_solve_s = dc.solve - dc_prev.solve; ct.gpu_push_upload_s = dc.upload - dc_prev.upload; ct.gpu_push_admit_s = dc.admit - dc_prev.admit;
      ct.gpu_push_device_s = dc.device - dc_prev.device; ct.gpu_push_assemble_s = dc.assemble - dc_prev.assemble; ct.gpu_push_total_s = dc.total - dc_prev.total;
      dc_prev = dc;
    }
    if (out) {
      const auto row = coordinator_row(ct, s, cin.carriers);
      out->coord_bin.write(reinterpret_cast<const char*>(row.data()), std::streamsize(row.size() * sizeof(double)));
    }
    ct = {};
    if (out && opt.progress_every > 0 && s % opt.progress_every == 0) {
      flush();
      std::cout << "step " << s << " t=" << o.hall->time_s() << " ms/step(last)=" << step_ms.back() << " e=" << after.n_e << " i=" << after.n_i
                << " n=" << after.n_n << " births(e,n)=" << st.electron_births << "," << st.neutral_births << " sinks=" << st.sinks
                << " events=" << st.events << " rss=" << rss_kib() << std::endl;
    }
    return "";
  }

  // ---------------------------------------------------------------- whole-run checkpoint
  J owner_record() const {
    const auto b = o.arena->full_time_boundary();
    J nb = {{"accepted_step", b->accepted_step}, {"accepted_arena_revision", b->accepted_arena_revision}, {"time_s", b64(b->time_s)},
            {"accepted_arena_sha256", b->accepted_arena_sha256}, {"field_sha256", b->field_sha256}, {"heat_sha256", b->heat_sha256},
            {"source_sha256", b->source_sha256}, {"rng_sha256", b->rng_sha256}, {"mesh_sha256", b->mesh_sha256}, {"material_sha256", b->material_sha256},
            {"frame_sha256", b->frame_sha256}, {"common_clock_sha256", b->common_clock_sha256}, {"authored_full_time_genesis", b->authored_full_time_genesis}};
    J src_cfg = J::array(), src_state = J::array();
    for (const auto& x : o.sources) { src_cfg.push_back(x->configuration_sha256()); src_state.push_back(x->state_sha256()); }
    return {{"native", {{"version", m::PicArena::source_identity()}, {"file", "native.json"}, {"state_prefix", o.arena->state_prefix()}, {"boundary", nb}}},
            {"hall", {{"version", {{"binding", v8::bind::hall_binding}, {"owner_type", v8::HallTraits<H>::name}, {"configuration_sha256", o.hall->configuration_sha256()}}},
                      {"dir", "hall"}, {"field_sha256", o.hall->field_sha256()}, {"heat_sha256", o.hall->heat_sha256()},
                      {"kinetic_ledgers_sha256", o.hall->kinetic_ledgers_sha256()}, {"time_s", b64(o.hall->time_s())}, {"revision", o.hall->revision()}}},
            {"source", {{"version", {{"owner", "electrode_heat_v7"}, {"configurations", src_cfg}}}, {"files", "source-<i>.bin"}, {"states", src_state},
                        {"accepted_source_seal", o.accepted_source_seal}}},
            {"arrival", {{"version", o.arrival_policy->configuration_sha256()}, {"file", "arrival.bin"}, {"state_sha256", o.arrival->state_sha256()},
                         {"ledger_head", o.arrival->ledger_head()}}},
            {"chemistry", {{"version", o.chemistry->configuration_sha256()}, {"file", "chemistry.json"}, {"rng_sha256", o.chemistry->accepted_rng_sha256()}}},
            {"coordinator", {{"version", "coordinator-v8.0.1"}, {"files", {"coordinator.json", "controller-ring.bin"}}, {"receipt_head", receipt_head},
                             {"configuration_sha256", v8::observe(o).configuration}}}};
  }
  J coordinator_json() const {
    return {{"schema", "coordinator-v8-coordinator-state-v1"}, {"step", step}, {"time_s", b64(o.hall->time_s())},
            {"accepted_source_seal", o.accepted_source_seal}, {"accepted_source_states", o.accepted_source_states},
            {"accepted_arrival_state", o.accepted_arrival_state}, {"receipt_head", receipt_head}, {"end_time_s", b64(end_time)},
            {"step_target", step_target},
            {"run_state", run_state_json(st)},
            {"controller", {{"dt_prev_s", b64(cst.dt_prev_s)}, {"rung_prev", cst.rung_prev}, {"hold", cst.hold}, {"chem_rate", b64(cst.chem_rate)},
                            {"births_rate", b64(cst.births_rate)}, {"ring", cst.ring.size()}, {"text", ctl.text()}, {"sha256", ctl.sha256()}}},
            {"rates", {{"previous", rate_json(o.rate_previous)}, {"current", rate_json(o.rate_current)}}},
            {"statistics", {{"criteria_histogram", criteria_histogram}, {"refused_attempts", refused_attempts_total}, {"dt_min_used_s", b64(dt_min_used)},
                            {"dt_max_used_s", b64(dt_max_used)}}},
            {"genesis", o.genesis}};
  }
  fs::path save_checkpoint() {
    fs::create_directories(checkpoint_root);
    const auto final_dir = checkpoint_root / ck::step_name(std::uint64_t(step));
    const auto tmp = fs::path(final_dir.string() + ".tmp");
    if (fs::exists(tmp)) fs::remove_all(tmp);
    fs::create_directories(tmp);
    const bool crash = opt.crash_step == step;
    o.arena->save((tmp / "native.json").string());
    o.hall->save(tmp / "hall");
    if (crash && opt.crash_stage == 0) crash_hook(0);
    for (std::size_t i = 0; i < o.sources.size(); ++i) o.sources[i]->save((tmp / ("source-" + std::to_string(i) + ".bin")).string());
    ck::write_file(tmp / "arrival.bin", o.arrival->save());
    ck::write_file(tmp / "chemistry.json", o.chemistry->checkpoint().dump() + "\n");
    ck::write_file(tmp / "coordinator.json", coordinator_json().dump(1) + "\n");
    ck::write_file(tmp / "controller-ring.bin", ring_bytes(cst));
    ck::fsync_tree(tmp);
    J body = {{"schema", ck::manifest_schema}, {"step", step}, {"time_s", b64(o.hall->time_s())}, {"ranks_recorded", opt.ranks},
              {"mpi_ranks_executed", 1}, {"rank_metadata_is_execution_evidence", false},
              {"bindings", {{"hall", v8::bind::hall_binding}, {"native", v8::bind::native_binding}, {"pair", v8::bind::pair_name},
                            {"hall_owner", v8::HallTraits<H>::name}}},
              {"run", options_json()}, {"controller_text", ctl.text()}, {"owners", owner_record()}, {"files", ck::file_table(tmp)},
              {"genesis_checkpoint", step == 0 ? J(nullptr) : J{{"dir", ck::step_name(0)}, {"manifest_body_sha256", genesis_manifest_sha256}}},
              {"previous_manifest_body_sha256", last_manifest_sha256}};
    const auto dir = ck::commit(tmp, final_dir, body, crash ? &crash_hook : nullptr);
    last_manifest_sha256 = nw::sha(body.dump());
    if (step == 0) genesis_manifest_sha256 = last_manifest_sha256;
    return dir;
  }
};

// ---------------------------------------------------------------- restore
std::string js(const J& j) { return j.template get<std::string>(); }
// Owner versions of this binary against the manifest (refused on any difference).
template <class H> void check_versions(const J& body, const v8::OwnersV8<H>* o = nullptr) {
  const auto& ow = body.at("owners");
  for (const char* k : {"native", "hall", "source", "arrival", "chemistry", "coordinator"})
    if (!ow.contains(k)) throw std::runtime_error(std::string("coordinator v8 restore: owner missing from the manifest: ") + k);
  if (js(ow.at("native").at("version")) != m::PicArena::source_identity()) throw std::runtime_error("coordinator v8 restore: native owner version differs");
  const auto& hv = ow.at("hall").at("version");
  if (js(hv.at("binding")) != v8::bind::hall_binding || js(hv.at("owner_type")) != v8::HallTraits<H>::name)
    throw std::runtime_error("coordinator v8 restore: Hall owner binding or owner type differs");
  if (js(ow.at("coordinator").at("version")) != "coordinator-v8.0.1") throw std::runtime_error("coordinator v8 restore: coordinator version differs");
  if (js(body.at("bindings").at("pair")) != v8::bind::pair_name || js(body.at("bindings").at("native")) != v8::bind::native_binding)
    throw std::runtime_error("coordinator v8 restore: owner binding pair differs");
  if (o) {
    if (js(hv.at("configuration_sha256")) != o->hall->configuration_sha256()) throw std::runtime_error("coordinator v8 restore: Hall configuration differs");
    J cfg = J::array();
    for (const auto& x : o->sources) cfg.push_back(x->configuration_sha256());
    if (js(ow.at("source").at("version").at("owner")) != "electrode_heat_v7" || ow.at("source").at("version").at("configurations") != cfg)
      throw std::runtime_error("coordinator v8 restore: source owner version/configuration differs");
    if (js(ow.at("arrival").at("version")) != o->arrival_policy->configuration_sha256()) throw std::runtime_error("coordinator v8 restore: arrival owner version differs");
    if (js(ow.at("chemistry").at("version")) != o->chemistry->configuration_sha256()) throw std::runtime_error("coordinator v8 restore: chemistry version differs");
  }
}
m::FullTimeAcceptedBoundary boundary_of(const J& j) {
  m::FullTimeAcceptedBoundary b;
  b.accepted_step = j.at("accepted_step"); b.accepted_arena_revision = j.at("accepted_arena_revision"); b.time_s = from_b64(j.at("time_s"));
  b.accepted_arena_sha256 = j.at("accepted_arena_sha256"); b.field_sha256 = j.at("field_sha256"); b.heat_sha256 = j.at("heat_sha256");
  b.source_sha256 = j.at("source_sha256"); b.rng_sha256 = j.at("rng_sha256"); b.mesh_sha256 = j.at("mesh_sha256"); b.material_sha256 = j.at("material_sha256");
  b.frame_sha256 = j.at("frame_sha256"); b.common_clock_sha256 = j.at("common_clock_sha256"); b.authored_full_time_genesis = j.at("authored_full_time_genesis");
  return b;
}
// The owners after the step recorded in the manifest (digests) must equal the live owners.
template <class H> void check_owner_digests(const v8::OwnersV8<H>& o, const J& body, const std::string& receipt_head, const char* what) {
  const auto& ow = body.at("owners");
  const auto fail = [&](const std::string& k) { throw std::runtime_error(std::string("coordinator v8 restore (") + what + "): " + k + " differs from the manifest"); };
  if (o.arena->state_prefix() != js(ow.at("native").at("state_prefix"))) fail("native accepted state");
  const auto& hl = ow.at("hall");
  if (o.hall->field_sha256() != js(hl.at("field_sha256")) || o.hall->heat_sha256() != js(hl.at("heat_sha256")) ||
      o.hall->kinetic_ledgers_sha256() != js(hl.at("kinetic_ledgers_sha256")))
    fail("Hall field/heat/kinetic-ledger state");
  J src = J::array();
  for (const auto& x : o.sources) src.push_back(x->state_sha256());
  if (src != ow.at("source").at("states") || o.accepted_source_seal != js(ow.at("source").at("accepted_source_seal"))) fail("source owner states");
  if (o.arrival->state_sha256() != js(ow.at("arrival").at("state_sha256")) || o.arrival->ledger_head() != js(ow.at("arrival").at("ledger_head"))) fail("arrival owner state");
  if (o.chemistry->accepted_rng_sha256() != js(ow.at("chemistry").at("rng_sha256"))) fail("chemistry RNG state");
  if (receipt_head != js(ow.at("coordinator").at("receipt_head"))) fail("receipt chain head");
  if (v8::observe(o).configuration != js(ow.at("coordinator").at("configuration_sha256"))) fail("coordinator configuration");
}
// Direct restore of every owner from one checkpoint directory (native restore_full_time first). Throws
// m::FullTimeResourceRefusal (restore_deferred) when native cannot restore that state from disk.
template <class H> v8::OwnersV8<H> restore_owners(const Inputs& in, const fs::path& dir, const J& body) {
  auto rq = request(in, {});
  const auto& ow = body.at("owners");
  v8::OwnersV8<H> o;
  o.hall_context = in.hall_context;
  o.model = in.model;
  o.settings = in.settings;
  o.source_options = in.source_options;
  o.neutral_dt_max_s = opt.neutral_dt_max_s;
  o.arena = m::PicArena::restore_full_time((dir / "native.json").string(), rq.options, rq.states, rq.channels, boundary_of(ow.at("native").at("boundary")));
  const auto census = o.arena->rate_population(in.settings.maximum_carriers);
  o.hall = v8::HallTraits<H>::restore(in.model, rq.options, census, dir / "hall", in.hall_context);
  o.slices = in.setup.slices;
  o.patches = in.setup.patches;
  for (std::size_t i = 0; i < o.slices.size(); ++i) {
    o.sources.push_back(v8::bind_source(*o.hall, o.slices[i], in.source_options));
    o.sources.back()->load((dir / ("source-" + std::to_string(i) + ".bin")).string());
  }
  o.arrival_policy = in.setup.arrival;
  o.arrival = pa::Owner::restore(o.arrival_policy, ck::read_file(dir / "arrival.bin"));
  o.chemistry = std::make_unique<cv::ParentSubflow>(*o.arena, rq.grid_counts, rq.binding, rq.limits, *rq.bridge, *rq.scope, rq.clock_seed, rq.product_seed, true);
  o.chemistry->restore(J::parse(ck::read_file(dir / "chemistry.json")));
  const auto cj = J::parse(ck::read_file(dir / "coordinator.json"));
  o.accepted_source_seal = cj.at("accepted_source_seal");
  o.accepted_source_states = cj.at("accepted_source_states").template get<std::vector<std::string>>();
  o.accepted_arrival_state = cj.at("accepted_arrival_state");
  o.genesis = cj.at("genesis");
  v8::InputsV8 x;
  x.states = rq.states;
  x.channels = rq.channels;
  x.species_sha256 = v8::species_sha256(rq.states, rq.channels);
  x.frame_sha256 = v8::frame_sha256(rq.options, *in.model);
  std::vector<std::string> sig;
  for (const auto& s : o.sources) sig.push_back(s->configuration_sha256());
  x.source_configuration_sha256 = v8::source_configuration_sha256(o.slices, sig);
  x.settings_text = in.settings.json_text();
  x.maximum_native_records = in.settings.maximum_native_records;
  x.maximum_carriers = in.settings.maximum_carriers;
  o.inputs = x;
  o.rate_previous = rates_of(cj.at("rates").at("previous"));
  o.rate_current = rates_of(cj.at("rates").at("current"));
  configure_owners(o);
  check_versions(body, &o);
  const auto seen = v8::observe(o);
  v8::need(nw::common_commit::owner_values_match(seen), "restored owners disagree with the accepted native boundary");
  check_owner_digests(o, body, cj.at("receipt_head"), "direct restore");
  return o;
}
// Coordinator state of a checkpoint into a Run built on restored owners.
template <class H> void restore_coordinator(Run<H>& r, const fs::path& dir) {
  const auto cj = J::parse(ck::read_file(dir / "coordinator.json"));
  if (js(cj.at("controller").at("text")) != ctl.text()) throw std::runtime_error("coordinator v8 restore: controller configuration differs");
  r.step = cj.at("step");
  r.receipt_head = cj.at("receipt_head");
  r.end_time = from_b64(cj.at("end_time_s"));
  r.step_target = cj.at("step_target");
  if (r.step_target < 0) throw std::runtime_error("coordinator v8 restore: negative step target");
  r.st = run_state_of(cj.at("run_state"));
  const auto& c = cj.at("controller");
  r.cst.dt_prev_s = from_b64(c.at("dt_prev_s")); r.cst.rung_prev = c.at("rung_prev"); r.cst.hold = c.at("hold");
  r.cst.chem_rate = from_b64(c.at("chem_rate")); r.cst.births_rate = from_b64(c.at("births_rate"));
  ring_of(ck::read_file(dir / "controller-ring.bin"), r.cst, std::size_t(r.o.hall->mesh().ne), std::size_t(r.o.hall->mesh().nb));
  if (r.cst.ring.size() != c.at("ring").template get<std::size_t>()) throw std::runtime_error("coordinator v8 restore: controller ring size");
  const auto& s = cj.at("statistics");
  r.criteria_histogram = s.at("criteria_histogram").template get<std::map<std::string, long>>();
  r.refused_attempts_total = s.at("refused_attempts");
  r.dt_min_used = from_b64(s.at("dt_min_used_s")); r.dt_max_used = from_b64(s.at("dt_max_used_s"));
  r.population = v8::accepted_population(r.o);
  r.before = totals(r.o, r.population);
}
// The coordinator state of a Run equals the one recorded in a checkpoint (replay check).
template <class H> bool same_coordinator_state(const Run<H>& r, const fs::path& dir) {
  const auto cj = J::parse(ck::read_file(dir / "coordinator.json"));
  auto mine = r.coordinator_json();
  return mine == cj && ring_bytes(r.cst) == ck::read_file(dir / "controller-ring.bin");
}

// ---------------------------------------------------------------- run mode
template <class H> int run_mode(const fs::path& outdir) {
  std::signal(SIGTERM, on_signal);
  std::signal(SIGINT, on_signal);
  // Restore: the manifest's run options are the run's configuration.
  fs::path restore_dir;
  J restore_body, restore_info = J::object();
  if (!opt.resume.empty() || !opt.restore_checkpoint.empty()) {
    J skipped = J::array();
    restore_dir = !opt.restore_checkpoint.empty() ? opt.restore_checkpoint : ck::newest(opt.resume, &skipped);
    restore_body = ck::verify(restore_dir);
    const auto requested_owner = opt.hall_owner;
    apply_options(restore_body.at("run"));
    if (opt.hall_owner != std::string(v8::HallTraits<H>::name) || requested_owner != opt.hall_owner)
      throw std::runtime_error("coordinator v8 restore: Hall owner type of the checkpoint (" + opt.hall_owner + ") differs from the launch (" + requested_owner + ")");
    check_versions<H>(restore_body);
    restore_info = {{"checkpoint", restore_dir.string()}, {"step", restore_body.at("step")}, {"skipped", skipped},
                    {"ranks_recorded", restore_body.at("ranks_recorded")}, {"ranks_now", opt.ranks}};
  }
  finish_configuration();
  const auto in = load_inputs(opt.fixture);
  J summary = {{"schema", "coordinator-v6-run-summary-v1"}, {"fixture", opt.fixture}, {"field_backend", opt.field}, {"push", opt.push},
               {"dt_s", opt.dt_s}, {"settings", in.settings.json_text()}, {"hall_configuration_sha256", in.model->configuration_sha256()}};
  if (backends.device) summary["device"] = J::parse(backends.device->device_json());
  const auto tg = steady::now();
  std::unique_ptr<Run<H>> run;
  Outputs files;
  files.dir = outdir;
  const fs::path root = !opt.checkpoint_dir.empty() ? opt.checkpoint_dir : outdir / "checkpoints";
  if (restore_dir.empty()) {
    run = std::make_unique<Run<H>>(in, genesis<H>(in, opt.fixture));
    summary["genesis"] = run->o.genesis;
    summary["genesis_s"] = since(tg);
    summary["motion_provider"] = run->o.motion->name();
    write_text(outdir / "genesis.json", run->o.genesis.dump(1) + "\n");
    run->start_fresh(summary);
    run->open_outputs(files, true);
    run->checkpoint_root = root;
    if (opt.checkpoint_every > 0) run->save_checkpoint();
  } else {
    // Direct restore (every owner from this checkpoint); when native refuses a nonempty disk restore
    // (native v7.1: restore_deferred, A05/A11), verified replay from the run's genesis checkpoint.
    const long target = restore_body.at("step");
    std::string mode = "direct";
    try {
      run = std::make_unique<Run<H>>(in, restore_owners<H>(in, restore_dir, restore_body));
      restore_coordinator(*run, restore_dir);
    } catch (const m::FullTimeResourceRefusal& e) {
      if (e.kind() != m::FullTimeResourceRefusal::Kind::restore_deferred) throw;
      mode = "replay";
      restore_info["native_direct_restore"] = e.what();
      const auto& g = restore_body.at("genesis_checkpoint");
      if (g.is_null()) throw;
      const auto gdir = restore_dir.parent_path() / g.at("dir").template get<std::string>();
      const auto gbody = ck::verify(gdir);
      if (nw::sha(gbody.dump()) != js(g.at("manifest_body_sha256"))) throw std::runtime_error("coordinator v8 restore: genesis checkpoint differs from the one the checkpoint names");
      if (gbody.at("run") != restore_body.at("run")) throw std::runtime_error("coordinator v8 restore: genesis checkpoint run options differ");
      run = std::make_unique<Run<H>>(in, restore_owners<H>(in, gdir, gbody));
      restore_coordinator(*run, gdir);
      const auto t0 = steady::now();
      while (run->step < target) {
        const auto stop = run->advance();
        if (!stop.empty()) throw std::runtime_error("coordinator v8 restore: replay stopped at step " + std::to_string(run->step + 1) + ": " + stop);
      }
      restore_info["replayed_steps"] = run->steps_this_process;
      restore_info["replay_s"] = since(t0);
      // Stopping targets are operational metadata and do not enter advance().
      // Earlier continuation launches may have changed this target after genesis.
      const auto saved_targets = J::parse(ck::read_file(restore_dir / "coordinator.json"));
      run->end_time = from_b64(saved_targets.at("end_time_s"));
      run->step_target = saved_targets.at("step_target");
      if (run->step_target < 0) throw std::runtime_error("coordinator v8 restore: negative step target");
      // The replayed owners must be the checkpoint's: manifest digests, receipt chain, coordinator state, and
      // every owner's own checkpoint restored directly beside them (Hall, sources, arrival, chemistry, native bytes).
      check_owner_digests(run->o, restore_body, run->receipt_head, "replay");
      J cross = J::object();
      {
        const auto census = run->o.arena->rate_population(in.settings.maximum_carriers);
        auto hall = v8::HallTraits<H>::restore(in.model, run->o.arena->native_options(), census, restore_dir / "hall", in.hall_context);
        cross["hall"] = hall->field_sha256() == run->o.hall->field_sha256() && hall->heat_sha256() == run->o.hall->heat_sha256() &&
                        hall->kinetic_ledgers_sha256() == run->o.hall->kinetic_ledgers_sha256();
        bool sources = true;
        for (std::size_t i = 0; i < run->o.sources.size(); ++i) {
          auto s = v8::bind_source(*run->o.hall, run->o.slices[i], in.source_options);
          s->load((restore_dir / ("source-" + std::to_string(i) + ".bin")).string());
          sources = sources && s->state_sha256() == run->o.sources[i]->state_sha256();
        }
        cross["sources"] = sources;
        cross["arrival"] = pa::Owner::restore(run->o.arrival_policy, ck::read_file(restore_dir / "arrival.bin"))->state_sha256() == run->o.arrival->state_sha256();
        cross["chemistry"] = J::parse(ck::read_file(restore_dir / "chemistry.json")) == run->o.chemistry->checkpoint();
        const auto scratch = outdir / "replay-native-check.json";
        run->o.arena->save(scratch.string());
        cross["native_bytes"] = ck::read_file(scratch) == ck::read_file(restore_dir / "native.json");
        fs::remove(scratch);
        cross["coordinator_state"] = same_coordinator_state(*run, restore_dir);
      }
      restore_info["cross_checks"] = cross;
      for (const auto& [k, v] : cross.items())
        if (!v.template get<bool>()) throw std::runtime_error("coordinator v8 restore: replayed " + k + " differs from the checkpoint's own record");
      run->steps_this_process = 0;
      run->step_ms.clear();
      run->sub_rows.clear();
      run->parity_compared = run->parity_identical = run->max_swap_alloc = 0;
      run->L = Ledger{};
      run->C = Ledger{};
    }
    restore_info["mode"] = mode;
    restore_info["restore_s"] = since(tg);
    summary["restored_from"] = restore_info;
    summary["genesis"] = run->o.genesis;
    summary["motion_provider"] = run->o.motion->name();
    write_text(outdir / "restore.json", restore_info.dump(1) + "\n");
    write_text(outdir / "genesis.json", run->o.genesis.dump(1) + "\n");
    run->open_outputs(files, false);
    run->checkpoint_root = root;
    run->genesis_manifest_sha256 = restore_body.at("genesis_checkpoint").is_null() ? nw::sha(restore_body.dump())
                                                                                  : restore_body.at("genesis_checkpoint").at("manifest_body_sha256").template get<std::string>();
    run->last_manifest_sha256 = nw::sha(restore_body.dump());
    if (root != restore_dir.parent_path() && opt.checkpoint_every > 0) {
      // A new checkpoint root needs the genesis checkpoint of the run for later replays: link it by copy.
      const auto g = restore_dir.parent_path() / ck::step_name(0);
      if (fs::exists(g) && !fs::exists(root / ck::step_name(0))) { fs::create_directories(root); fs::copy(g, root / ck::step_name(0), fs::copy_options::recursive); }
    }
  }
  auto& r = *run;
  // Replay must compare the checkpoint under its original target. Only after
  // that verification may this launch extend the run or select a step limit.
  if (!restore_dir.empty()) {
    const double saved_end = r.end_time;
    r.end_time = v8::RunTargetRequest{opt.end_time_requested, opt.window_requested, opt.steps_requested,
                                     opt.end_time_s, opt.window_s, opt.steps}.after_restore(saved_end, r.o.hall->time_s());
    if (opt.steps_requested) r.step_target = opt.steps;
    summary["continuation_target"] = {{"saved_end_time_s", saved_end}, {"end_time_s", r.end_time},
                                      {"absolute_steps", r.step_target}, {"explicit_end_time", opt.end_time_requested},
                                      {"explicit_window", opt.window_requested}, {"explicit_steps", opt.steps_requested}};
  }
  summary["end_time_s"] = r.end_time;
  summary["step_target"] = r.step_target;
  summary["controller"] = ctl.text();
  summary["hall_binding"] = v8::bind::hall_binding;
  summary["hall_owner"] = v8::HallTraits<H>::name;
  summary["binding_pair"] = v8::bind::pair_name;
  summary["mesh"] = mesh_record;
  summary["collision_grid"] = opt.collision_grid;
  summary["host_work"] = opt.host_work;
  summary["native_binding"] = v8::bind::native_binding;
  summary["hall_configurable_cap"] = v8::bind::hall_configurable_cap;
  summary["field_owner_dt_max_s"] = v8::field_owner_dt_max(*in.model);
  summary["coupled_owner_dt_max_s"] = v8::owner_dt_max(*in.model, in.source_options);
  summary["source_owner"] = {{"owner", "electrode_heat_v7 (pic-electrode-heat-source-v7.2-prepared-state-20261101)"}, {"interval_cap_s", in.source_options.interval.cap_s},
                             {"emission_relative_error", in.source_options.interval.emission_relative_error},
                             {"desorption_relative_error", in.source_options.interval.desorption_relative_error},
                             {"safety", in.source_options.interval.safety}, {"count_floor", in.source_options.interval.count_floor},
                             {"interval_preset", opt.source_interval_preset}, {"surface", opt.surface}};
  summary["neutral_stream_dt_max_s"] = opt.neutral_dt_max_s;
  summary["controller_sha256"] = ctl.sha256();
  summary["execution"] = {{"mpi_ranks_executed", 1}, {"requested_rank_metadata", opt.ranks},
                           {"distributed_restart_qualified", false}};
  summary["physics_qualification"] = {{"fully_coupled_water_breakdown", false},
                                       {"native_configuration", "v5 compatibility; nonempty restore uses verified replay"},
                                       {"chemistry", "parent ionization only in this coordinator binding"},
                                       {"thermal", "solid-only Hall thermal binding"},
                                       {"weak_positive_emission_field", "refused; no qualified launch model for 0 < F < 1 V/m"}};
  summary["receipt_chain_genesis"] = r.receipt_head;
  r.dc_prev = device_counters(in, r.device_motion);
  std::string stop = "steps";
  const auto wall = steady::now();
  J checkpoints = J::array();
  for (;;) {
    const long next = r.step + 1;
    if (r.finished(next)) { stop = r.end_time > 0 ? "end_time" : "steps"; break; }
    if (stop_signal) { stop = "signal"; break; }
    if (!opt.stop_file.empty() && next % 100 == 0 && fs::exists(opt.stop_file)) { stop = "stop_file"; break; }
    const auto s = r.advance();
    if (!s.empty()) { stop = s; break; }
    if (opt.checkpoint_every > 0 && r.step % opt.checkpoint_every == 0) {
      const auto tc = steady::now();
      const auto dir = r.save_checkpoint();
      checkpoints.push_back({{"step", r.step}, {"dir", dir.string()}, {"write_s", since(tc)}});
    }
  }
  r.flush();
  if (r.legacy) r.population = r.o.arena->rate_population(r.o.inputs.maximum_carriers);
  if (opt.frame_every > 0 && r.step % opt.frame_every != 0) files.frames.write(r.o, r.step, r.population);
  summary["checkpoints"] = checkpoints;
  summary["receipt_chain_head"] = r.receipt_head;
  summary["criteria_histogram"] = r.criteria_histogram;
  summary["refused_attempts"] = r.refused_attempts_total;
  const double t0 = r.o.genesis.contains("time_s") ? r.o.genesis.at("time_s").template get<double>() : in.model->drive().start();
  summary["dt_used_s"] = {{"min", std::isfinite(r.dt_min_used) ? J(r.dt_min_used) : J(nullptr)}, {"max", r.dt_max_used},
                          {"mean", r.step > 0 ? (r.o.hall->time_s() - t0) / double(r.step) : 0.}};
  if (!r.sub_rows.empty()) {
    J med = J::object();
    for (std::size_t c = 1; c < subphase_columns.size(); ++c) {
      std::vector<double> v;
      for (const auto& x : r.sub_rows) v.push_back(x[c]);
      std::sort(v.begin(), v.end());
      med[subphase_columns[c]] = v[v.size() / 2];
    }
    summary["subphase_median"] = med;
  }
  summary["coordinator"] = "v8";
  summary["emission_statistics"] = emission_statistics_json();
  summary["steps_completed"] = r.step;
  summary["steps_this_process"] = r.steps_this_process;
  summary["stop"] = stop;
  summary["wall_s"] = since(wall);
  summary["final_time_s"] = r.o.hall->time_s();
  summary["simulated_s"] = r.o.hall->time_s() - t0;
  if (!r.step_ms.empty()) {
    auto s = r.step_ms;
    std::sort(s.begin(), s.end());
    double mean = 0;
    for (double x : r.step_ms) mean += x;
    summary["step_ms"] = {{"first", r.step_ms.front()}, {"median", s[s.size() / 2]}, {"mean", mean / double(s.size())},
                          {"p90", s[std::size_t(.9 * double(s.size() - 1))]}, {"max", s.back()}};
  }
  summary["max_allocations_in_swap_window"] = r.max_swap_alloc;
  const auto& st = r.st;
  summary["closure"] = {{"worst_abs_step_energy_residual_J", st.worst_step_residual}, {"worst_gauss_over_bound", st.worst_gauss_ratio},
                        {"worst_divB_over_bound", st.worst_divB_ratio}, {"worst_abs_total_charge_change_C", st.worst_q_change},
                        {"water_inventory_exact_every_step", st.water_exact}, {"H_O_nuclei_exact_every_step", st.nuclei_exact},
                        {"worst_abs_reaction_defect_J", st.worst_reaction_defect}, {"cumulative_reaction_defect_J", double(st.reaction_defect_cum)}};
  summary["totals"] = {{"electron_births", st.electron_births}, {"h2o_births", st.neutral_births}, {"sinks", st.sinks}, {"parent_events", st.events}};
  const auto fin = totals(r.o, r.population);
  summary["final"] = {{"carriers_electron", fin.n_e}, {"carriers_ion", fin.n_i}, {"carriers_h2o", fin.n_n},
                      {"source_stock_molecules", v8::exact_value(fin.stock)}, {"hall_kinetic_ledgers", r.o.hall->kinetic_ledgers().sinks.json()},
                      {"hall_birth_ledgers", r.o.hall->kinetic_ledgers().births.json()}, {"native_flow_totals_births", r.o.arena->full_time_flow_totals().births},
                      {"native_flow_totals_sunk", r.o.arena->full_time_flow_totals().sunk_carriers},
                      {"thermal_max_T_K", r.o.hall->thermal().max_temperature_K()}, {"chemistry_coverage", r.o.chemistry->accepted_coverage()}};
  {
    const auto hs = r.o.arena->full_time_history_status();
    summary["native_history"] = {{"anchors", hs.anchors}, {"anchor_bytes", hs.anchor_bytes}, {"window_rows", hs.window_rows},
                                 {"cumulative_rows", hs.cumulative_rows}, {"cumulative_bytes", hs.cumulative_bytes}};
  }
  if (opt.push == "gpu") summary["push_parity"] = {{"iterates_compared", r.parity_compared}, {"iterates_identical", r.parity_identical}};
  if (opt.field == "gpu" && opt.hall_owner == "owner") {
    const auto s = hg::gpu_field_solver_stats(in.model->solver());
    summary["gpu_field_solves"] = {{"solves", s.solves}, {"total_s", s.total_s}};
  }
  summary["forbidden_native_calls"] = forbidden_calls;
  if (const auto g = std::dynamic_pointer_cast<hg::GpuMotionProvider>(r.device_motion)) {
    const auto& x = g->stats();
    summary["gpu_push"] = {{"runs", x.runs}, {"upload_carriers_s", x.upload_carriers_s}, {"admit_s", x.admit_s}, {"device_run_s", x.device_run_s},
                           {"assemble_s", x.assemble_s}, {"total_s", x.total_s}};
  }
  summary["frames"] = files.frames.frames;
  summary["rss_kib_final"] = rss_kib();
  write_text(outdir / "run-summary.json", summary.dump(1) + "\n");
  std::cout << "run finished: steps=" << r.step << " stop=" << stop << std::endl;
  const bool ok = (stop == "steps" || stop == "end_time" || stop == "stop_file" || stop == "signal") && r.parity_compared == r.parity_identical &&
                  st.water_exact && st.nuclei_exact && forbidden_calls == 0 && r.max_swap_alloc == 0;
  return ok ? 0 : 1;
}

// ---------------------------------------------------------------- vacuum step cost (v7.1)
int vacuum_mode(const fs::path& out) {
  const auto in = load_inputs("pilot");
  const auto& model = *in.model;
  const auto seed = model.dc_seed();
  hc::ChargeLedger q;
  const auto& g = model.mesh();
  q.particles.net_C.resize(g.nn); q.particles.gross_C.resize(g.nn); q.conductor_C.resize(g.nn); q.exterior_C.resize(g.nn);
  q.conductor_declaration = hc::sha256("vacuum pre-run: DC conductor charge"); q.exterior_declaration = hc::sha256("vacuum pre-run: DC exterior charge");
  const auto field = hc::field_charge(g, seed);
  for (int n = 0; n < g.nn; ++n) {
    if (g.node_body[n] >= 0 || g.short_node[n]) q.conductor_C[n] = field[n];
    else if (g.outer_node[n]) q.exterior_C[n] = field[n];
  }
  hc::Owner hall(in.model, seed, q, v8::hall_geometry(g), {});
  struct Boundary { std::string field_sha256, heat_sha256; double time_s = 0; std::uint64_t accepted_step = 0; };
  J summary = {{"schema", "coordinator-v6-vacuum-prerun-cost-v1"}, {"field_backend", opt.field}, {"dt_s", opt.dt_s}};
  std::string refusal;
  try { (void)hall.prepare_vacuum(std::nextafter(model.dt_max(), 1.)); } catch (const std::exception& e) { refusal = e.what(); }
  summary["dt_above_cap_refusal"] = refusal;
  std::vector<double> ms;
  const auto wall = steady::now();
  for (long s = 0; s < opt.steps; ++s) {
    const auto t0 = steady::now();
    auto p = hall.prepare_vacuum(opt.dt_s);
    Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), p.end_s(), p.endpoint_revision()};
    if (!hall.preflight_vacuum_final(p, b)) throw std::runtime_error("vacuum preflight refused");
    hall.publish_vacuum_prevalidated(std::move(p));
    ms.push_back(1e3 * since(t0));
  }
  summary["steps"] = opt.steps;
  summary["wall_s"] = since(wall);
  summary["final_time_s"] = hall.time_s();
  write_text(out / "vacuum-summary.json", summary.dump(1) + "\n");
  return refusal.empty() ? 1 : 0;
}

// ---------------------------------------------------------------- Hall-only pre-run in the coordinator Model
// The Hall v3.3 pre-run (test_v32 prerun: DC seed, the same declared conductor/exterior charge, empty census, landing
// grid m * SAMPLE_S with equal sub-steps of at most DT) on the coordinator's Model (source-slice thermal declaration,
// carrier limits, step cap, motion forms of this run). Hall checkpoints every --prerun-checkpoint-every and at the end
// (the Hall owner's own format: restorable by --hall-checkpoint). With --compare-hall-checkpoint REF, REF is restored with
// its own pre-run Model (the configuration recorded in REF: ThermalDeclaration{}, CarrierLimits{}, its dt_max and forms)
// and its accepted field words, charge ledger (free, conductor, exterior, residue, admitted Gauss ledger) and every slice's
// thermal state (temperatures, cumulative energies, clock) are compared word for word with this pre-run's end state.
struct PrerunBoundary { std::string field_sha256, heat_sha256, source_sha256, mesh_sha256, material_sha256, common_clock_sha256; double time_s = 0; std::uint64_t accepted_step = 0; };
J compare_words(const char* name, const hc::Vec& a, const hc::Vec& b) {
  long diff = 0;
  for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i) diff += std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i]);
  return {{"array", name}, {"size", a.size()}, {"same_size", a.size() == b.size()}, {"differing_words", diff}, {"identical", a.size() == b.size() && diff == 0}};
}
J compare_words(const char* name, const hc::WordVec& a, const hc::WordVec& b) {
  long diff = 0;
  for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i)
    diff += std::bit_cast<std::uint64_t>(a[i].hi) != std::bit_cast<std::uint64_t>(b[i].hi) || std::bit_cast<std::uint64_t>(a[i].lo) != std::bit_cast<std::uint64_t>(b[i].lo);
  return {{"array", name}, {"size", a.size()}, {"same_size", a.size() == b.size()}, {"differing_words", diff}, {"identical", a.size() == b.size() && diff == 0}};
}
// Plain foil planes of the Hall v3.x test fixtures (tests/motion_fixtures.hpp geometry()): the native geometry of the
// Hall v3.3 pre-run, bound into its checkpoint's motion mesh digest.
m::GapPicOptions plain_geometry(const axial::Mesh& g) {
  m::GapPicOptions o;
  o.nx = o.ny = g.nx; o.width_m = o.gap_m = g.side; o.depth_m = g.depth;
  o.finite_foils = o.finite_axial_extent = true; o.short_z_m = -g.depth / 2; o.feed_z_m = g.depth / 2;
  o.lower_foil = {g.side / 2 - .0025, g.side / 2 + .0025, g.side / 2 - .0027, g.side / 2 - .0025};
  o.upper_foil = {g.side / 2 - .0025, g.side / 2 + .0025, g.side / 2 + .0025, g.side / 2 + .0027};
  return o;
}
// The reference Hall checkpoint restored with its own pre-run Model (the configuration recorded in it: no thermal
// declaration, default carrier limits, its dt_max and motion forms; the plain native geometry) and compared word for word
// with the accepted state of `ours` (field words, charge ledger incl. residue and admitted Gauss ledger, every slice's
// thermal state and clock).
template <class H> J compare_hall_checkpoint(const Inputs& in, const H& ours) {
  const auto& g = in.model->mesh();
  const auto ref_json = J::parse(read_bytes(opt.compare_hall_checkpoint / "checkpoint.json"));
  const double ref_dt_max = ref_json.at("dt_max_s").template get<double>();
  const hc::MotionForms ref_forms{ref_json.at("motion_forms_native_v7").template get<bool>(), ref_json.at("neutral_stream_dt_max_s").template get<double>()};
  const auto spec = selected_mesh();
  auto ref_model = std::make_shared<const hc::Model>(spec, paths.trace, solver(spec), hc::CarrierLimits{}, hc::ThermalDeclaration{}, ref_dt_max, ref_forms);
  auto ref = v8::HallTraits<H>::restore(ref_model, plain_geometry(g), {}, opt.compare_hall_checkpoint, in.hall_context);
  typename v8::HallTraits<H>::Cache c1, c2;
  const auto A = v8::HallTraits<H>::accepted(*ref, c1);
  const auto B = v8::HallTraits<H>::accepted(ours, c2);
  bool identical = true;
  J cmp = J::array();
  cmp.push_back(compare_words("E", A.state->E, B.state->E));
  cmp.push_back(compare_words("B", A.state->B, B.state->B));
  cmp.push_back(compare_words("E_low", A.state->E_low, B.state->E_low));
  cmp.push_back(compare_words("B_low", A.state->B_low, B.state->B_low));
  cmp.push_back(compare_words("free_net_C", A.charge->particles.net_C, B.charge->particles.net_C));
  cmp.push_back(compare_words("free_gross_C", A.charge->particles.gross_C, B.charge->particles.gross_C));
  cmp.push_back(compare_words("conductor_C", A.charge->conductor_C, B.charge->conductor_C));
  cmp.push_back(compare_words("exterior_C", A.charge->exterior_C, B.charge->exterior_C));
  cmp.push_back(compare_words("residue_C", A.charge->residue_C, B.charge->residue_C));
  cmp.push_back(compare_words("admitted_gauss_C", A.charge->gauss_admitted_C, B.charge->gauss_admitted_C));
  for (const auto& x : cmp) identical = identical && x.at("identical").template get<bool>();
  const bool clock = A.state->time_s == B.state->time_s && A.state->revision == B.state->revision;
  identical = identical && clock;
  J heat = J::array();
  bool heat_same = true;
  for (int k = g.kshort; k < g.kfeed; ++k)
    for (int body = 0; body < 2; ++body) {
      const auto& a = ref->thermal().slice(body, k).state();
      const auto& b = ours.thermal().slice(body, k).state();
      std::vector<double> ca(a.cumulative.begin(), a.cumulative.end()), cb(b.cumulative.begin(), b.cumulative.end());
      const auto t = compare_words("T", a.T, b.T), cu = compare_words("cumulative", ca, cb);
      const bool same = t.at("identical").template get<bool>() && cu.at("identical").template get<bool>() && a.time == b.time && a.revision == b.revision;
      heat_same = heat_same && same;
      heat.push_back({{"body", body}, {"slice", k}, {"T", t}, {"cumulative", cu}, {"clock_identical", a.time == b.time && a.revision == b.revision}});
    }
  identical = identical && heat_same;
  return {{"reference", opt.compare_hall_checkpoint.string()}, {"reference_configuration_sha256", ref_json.at("configuration_sha256")},
          {"reference_model_configuration_sha256", ref_model->configuration_sha256()}, {"reference_restored_field_sha256", ref->field_sha256()},
          {"reference_restored_heat_sha256", ref->heat_sha256()}, {"ours_field_sha256", ours.field_sha256()}, {"ours_heat_sha256", ours.heat_sha256()},
          {"ours_configuration_sha256", ours.configuration_sha256()}, {"field_and_charge", cmp}, {"clock_identical", clock},
          {"heat_slices", heat}, {"heat_identical", heat_same}, {"all_identical", identical}};
}
// hallcompare: restore --hall-checkpoint (a Hall checkpoint of this coordinator's Model: coordinator geometry, source-slice
// declaration) and compare it with --compare-hall-checkpoint (restored with its own pre-run Model).
template <class H> int hallcompare_mode(const fs::path& out) {
  const auto in = load_inputs("pilot");
  auto ours = v8::HallTraits<H>::restore(in.model, v8::hall_geometry(in.model->mesh()), {}, opt.hall_checkpoint, in.hall_context);
  J summary = {{"schema", "coordinator-v8-hall-checkpoint-compare-v1"}, {"ours", opt.hall_checkpoint.string()}, {"hall_owner", v8::HallTraits<H>::name},
               {"ours_model_configuration_sha256", in.model->configuration_sha256()},
               {"ours_thermal_declaration_sha256", in.model->thermal_declaration().sha256()}};
  summary["compare"] = compare_hall_checkpoint<H>(in, *ours);
  const bool same = summary["compare"].at("all_identical").template get<bool>();
  write_text(out / "hall-compare.json", summary.dump(1) + "\n");
  std::cout << summary.dump(1).substr(0, 6000) << std::endl;
  return same ? 0 : 1;
}
template <class H> int prerun_mode(const fs::path& out) {
  const auto in = load_inputs("pilot");
  const auto& model = *in.model;
  const auto& g = model.mesh();
  const auto o = v8::hall_geometry(g);
  const auto seed = model.dc_seed();
  hc::ChargeLedger q;
  q.particles = hc::charge_of_census(g, o, std::vector<m::ArenaRateEntry>{});
  q.conductor_C.assign(g.nn, hc::Words()); q.exterior_C.assign(g.nn, hc::Words());
  q.conductor_declaration = hc::sha256("parity: declared DC conductor/PEC induced charge");
  q.exterior_declaration = hc::sha256("parity: declared DC exterior termination charge");
  const auto fc = hc::field_charge(g, seed);
  for (int k = 0; k < g.nn; ++k) { if (g.node_body[k] >= 0 || g.short_node[k]) q.conductor_C[k] = fc[k]; else if (g.outer_node[k]) q.exterior_C[k] = fc[k]; }
  auto res = v8::HallTraits<H>::make(in.model, seed, q, o, {}, in.hall_context);
  const double dt = opt.dt_s, sample_s = opt.prerun_sample_s, t_end = opt.prerun_end_s;
  v8::need(dt > 0 && dt <= model.dt_max() && t_end > seed.time_s && t_end <= model.drive().end() && sample_s > 0, "pre-run arguments");
  const double per = std::round(1.0 / sample_s);
  v8::need(std::abs(per * sample_s - 1) < 1e-12, "--prerun-sample must be the reciprocal of an integer");
  const long long m_end = std::llround(t_end * per);
  v8::need(double(m_end) / per == t_end, "--prerun-end must lie on the landing grid");
  long long mm = (long long)std::floor(seed.time_s * per) + 1;
  const long long ck_every = opt.prerun_checkpoint_s > 0 ? std::llround(opt.prerun_checkpoint_s * per) : 0;
  long steps = 0;
  J checkpoints = J::array();
  const auto t0 = steady::now();
  double worst_balance = 0, worst_residual = 0;
  for (; mm <= m_end; ++mm) {
    const double target = double(mm) / per;
    while (res->time_s() < target) {
      const double remaining = target - res->time_s();
      long long nleft = (long long)std::ceil(remaining / dt);
      double hstep = remaining / double(std::max<long long>(nleft, 1));
      if (hstep > dt) hstep = remaining / double(nleft + 1);
      if (nleft <= 1 && remaining <= dt) hstep = remaining;
      auto p = res->prepare_vacuum(hstep);
      PrerunBoundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), "", "", "", "", p.end_s(), p.endpoint_revision()};
      v8::need(res->preflight_vacuum_final(p, b), "pre-run preflight");
      const auto& rc = p.field_receipt();
      const double activity = std::abs(rc.field_delta_J) + std::abs(rc.plasma_work_J) + std::abs(rc.joule_J) + std::abs(rc.outer_loss_J) +
                              std::abs(rc.outer_supply_J) + std::abs(rc.circuit_work_J);
      worst_balance = std::max(worst_balance, std::abs(rc.balance_J) / (1e-20 + 1e-10 * activity));
      worst_residual = std::max(worst_residual, p.solve_record().independent_residual / p.solve_record().target);
      res->publish_vacuum_prevalidated(std::move(p));
      ++steps;
    }
    v8::need(res->time_s() == target, "landing time not reached exactly");
    if ((ck_every > 0 && mm % ck_every == 0) || mm == m_end) {
      char name[64];
      std::snprintf(name, sizeof name, "checkpoint-%.6fns", target * 1e9);
      const auto dir = out / name;
      fs::remove_all(dir);
      res->save(dir);
      checkpoints.push_back({{"time_s", target}, {"dir", dir.string()}, {"field_sha256", res->field_sha256()}, {"heat_sha256", res->heat_sha256()},
                             {"revision", res->revision()}});
      write_text(out / "checkpoints.json", checkpoints.dump(1) + "\n");
      std::cout << "pre-run checkpoint " << name << " steps=" << steps << " wall=" << since(t0) << std::endl;
    }
  }
  J summary = {{"schema", "coordinator-v8-prerun-v1"}, {"hall_owner", v8::HallTraits<H>::name}, {"dt_s", dt}, {"start_s", seed.time_s},
               {"end_time_s", res->time_s()}, {"steps", steps}, {"wall_s", since(t0)}, {"ms_per_step", 1e3 * since(t0) / double(std::max(1L, steps))},
               {"configuration_sha256", model.configuration_sha256()}, {"final_field_sha256", res->field_sha256()}, {"final_heat_sha256", res->heat_sha256()},
               {"worst_energy_balance_over_bound", worst_balance}, {"worst_independent_residual_over_target", worst_residual}, {"checkpoints", checkpoints},
               {"thermal_declaration_sha256", model.thermal_declaration().sha256()}};
  bool identical = true;
  if (!opt.compare_hall_checkpoint.empty()) {
    summary["compare"] = compare_hall_checkpoint<H>(in, *res);
    identical = summary["compare"].at("all_identical").template get<bool>();
  }
  write_text(out / "prerun.json", summary.dump(1) + "\n");
  std::cout << summary.dump(1).substr(0, 4000) << std::endl;
  return identical ? 0 : 1;
}
