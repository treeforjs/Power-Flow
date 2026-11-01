// Coordinator v8 driver (the coordinator v7.1 driver with source owner v7, the Hall owner type chosen at
// launch, the whole-run checkpoint and restart, genesis at a Hall checkpoint and the Hall-only pre-run in the
// coordinator Model). Modes:
//   run      genesis (fixture qual|pilot|dense|load; at the first C19 knot or at --hall-checkpoint) or a restore
//            (--resume ROOT: newest committed checkpoint below ROOT; --restore-checkpoint DIR), then consecutive
//            accepted steps with the v7.1 records (ledgers, digests, frames, probes, control ledger, step
//            receipts) and, with --checkpoint-every K, whole-run checkpoints (run_checkpoint_v8.hpp)
//   checks   the coordinator v7.1 qualification checks on the host Hall owner (source v7 at benchmark options)
//   vacuum   Hall-only vacuum step cost and the dt guard (v7.1)
//   prerun   Hall-only pre-run in the coordinator Model to --prerun-end with Hall checkpoints, optionally
//            compared word for word with a reference Hall checkpoint restored with its own pre-run Model
// usage: driver_v8 MODE TRACE H2O_DATA BINDING OUTER3 KERNEL_BINDING OUTDIR NATIVE_IDENTITY [options]
#include "genesis_v8.hpp"
#include "step_controller_v8.hpp"
#include "exact_weight_sum_v8.hpp"
#include "run_checkpoint_v8.hpp"
#include "run_target_v8.hpp"
#include "gpu/gpu_backends.hpp"
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <random>

extern "C" {
long long stage1_allocation_count() noexcept;
void stage1_set_allocation_countdown(long long) noexcept;
}
namespace v8 = raikou::native_water::v8;
namespace v7 = v8;   // the v7.1 check texts name the coordinator namespace v7::
namespace ck = raikou::native_water::v8::checkpoint;
namespace emi = ::electrode_heat_v6_1::emission;
namespace v3k = ::electrode_heat_v3;
namespace se = ::surface_emission;
namespace m = raikou::mfem_backend;
namespace h = ::hall_driver;
namespace hc = ::hall_corrective_v3;
namespace hg = ::hall_corrective_v3::gpu;
namespace hb = ::hall_source_binding_v7;
namespace cv = ::raikou::native_water::coverage_v2;
namespace nw = ::raikou::native_water;
namespace pc = ::raikou::pic_collisions;
namespace w = ::raikou::pic_water;
namespace pa = ::pf_arrival_v11;
namespace s7 = ::electrode_heat_v7;
namespace sv7 = ::surface_v7;
namespace eh = ::electrode_heat_v2;
namespace fs = std::filesystem;
using J = nw::Json;
using steady = std::chrono::steady_clock;

int forbidden_calls = 0;
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic12solve_fieldsEv(m::GapPic*) { ++forbidden_calls; throw std::runtime_error("forbidden native solver trap"); }
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic7advanceEd(m::GapPic*, double) { ++forbidden_calls; throw std::runtime_error("forbidden native advance trap"); }
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic13advance_emptyEd(m::GapPic*, double) { ++forbidden_calls; throw std::runtime_error("forbidden native empty advance trap"); }
extern "C" void __wrap__ZN6raikou12mfem_backend6GapPic16advance_observedEdRKSt8functionIFvRKNS0_11GapParticleEidEE(m::GapPic*, double,
                                                                                                                const m::GapPic::SinkObserver&) {
  ++forbidden_calls;
  throw std::runtime_error("forbidden native observed advance trap");
}

namespace raikou::native_water::v8 {
struct TrialTestAccessV7 {
  template <class H> static void set_can_publish_fault(TrialV8<H>& t, int k) { t.can_publish_fault_ = k; }
  template <class H> static void set_publish_probe(TrialV8<H>& t, void (*p)(int)) { t.publish_probe_ = p; }
  template <class H> static std::string& source_after(TrialV8<H>& t) { return t.source_after_; }
  template <class H> static auto& proposal(TrialV8<H>& t) { return t.field_heat_; }
  template <class H> static std::unique_ptr<m::PreparedArena>& native(TrialV8<H>& t) { return t.native_; }
  template <class H> static std::vector<std::optional<s7::Trial>>& source_trials(TrialV8<H>& t) { return t.source_trials_; }
  template <class H> static std::optional<pa::Staged>& staged(TrialV8<H>& t) { return t.staged_; }
  template <class H> static c7::StepClock& clock(TrialV8<H>& t) { return t.clock_; }
  template <class H> static void compare_v6_predicate(TrialV8<H>& t) { t.compare_v6_predicate_ = true; }
  template <class H> static std::pair<int, int> predicates(const TrialV8<H>& t) { return {t.v6_predicate_, t.v7_predicate_}; }
};
struct GenesisTestAccessV7 {
  template <class H>
  static OwnersV8<H> admit(const GenesisRequestV8& rq, typename GenesisV8<H>::PacketHook p, typename GenesisV8<H>::LedgerHook l, int* reached) {
    return GenesisV8<H>::admit_with_hooks(rq, p, l, reached);
  }
};
}  // namespace raikou::native_water::v8
using TA = v8::TrialTestAccessV7;

namespace {
volatile std::sig_atomic_t stop_signal = 0;
void on_signal(int) { stop_signal = 1; }
double since(steady::time_point a) { return std::chrono::duration<double>(steady::now() - a).count(); }
void write_text(const fs::path& p, const std::string& s) {
  std::ofstream f(p, std::ios::binary);
  f << s;
  f.close();
  if (!f) throw std::runtime_error("write " + p.string());
}
std::string read_bytes(const fs::path& p) {
  std::ifstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error("read " + p.string());
  return {std::istreambuf_iterator<char>(f), {}};
}
long rss_kib() {
  std::ifstream f("/proc/self/status");
  std::string line;
  while (std::getline(f, line))
    if (line.rfind("VmRSS:", 0) == 0) return std::stol(line.substr(6));
  return -1;
}
std::string hexd(double x) { std::ostringstream s; s << std::hexfloat << x; return s.str(); }
std::string bits(double x) { std::ostringstream s; s << std::hex << std::bit_cast<std::uint64_t>(x); return s.str(); }
// Exact binary64 record (16 hex digits of the bit pattern) and back.
std::string b64(double x) { std::ostringstream s; s << std::hex << std::setw(16) << std::setfill('0') << std::bit_cast<std::uint64_t>(x); return s.str(); }
double from_b64(const J& j) {
  if (!j.is_string() || j.template get<std::string>().size() != 16) throw std::runtime_error("checkpoint binary64 record");
  return std::bit_cast<double>(std::stoull(j.template get<std::string>(), nullptr, 16));
}
J b64s(const std::vector<double>& v) { J a = J::array(); for (double x : v) a.push_back(b64(x)); return a; }
std::vector<double> from_b64s(const J& j) { std::vector<double> v; for (const auto& x : j) v.push_back(from_b64(x)); return v; }
// x87 extended value (the 10 significant bytes) and back.
std::string ld_hex(long double x) {
  unsigned char b[sizeof(long double)] = {};
  std::memcpy(b, &x, sizeof x);
  std::ostringstream s;
  for (int i = 0; i < 10; ++i) s << std::hex << std::setw(2) << std::setfill('0') << int(b[i]);
  return s.str();
}
long double from_ld_hex(const J& j) {
  const auto t = j.template get<std::string>();
  if (t.size() != 20) throw std::runtime_error("checkpoint extended record");
  unsigned char b[sizeof(long double)] = {};
  for (int i = 0; i < 10; ++i) b[i] = static_cast<unsigned char>(std::stoul(t.substr(std::size_t(2 * i), 2), nullptr, 16));
  long double x;
  std::memcpy(&x, b, sizeof x);
  return x;
}

// ---------------------------------------------------------------- configuration
struct Options {
  std::string field = "cpu", push = "cpu", uuid, nvrtc, fixture = "pilot";
  long steps = 20, frame_every = 0, parity_every = 0, progress_every = 1000, checkpoint_every = 0;
  double end_time_s = 0, dt_s = 1e-14;
  bool end_time_requested = false, window_requested = false, steps_requested = false;
  double w0e = -1, w0n = -1, phi = -1;
  int faces = 4;
  std::size_t gpu_capacity = 65536;
  fs::path stop_file;
  // v7
  double window_s = 0;
  long seed_offset = 0;
  long blob_pairs = 400;
  double blob_weight = 1.6e11, blob_eV = 3;
  long manufactured_refusal_step = 0;
  double hall_dt_max_s = 1e-14;
  long probe_energy_every = 1;
  std::string host_work = "v7.1";
  // v7.1
  std::string mesh_file;
  std::string collision_grid = "mesh";
  long load_carriers = 0;
  double load_neutral_fraction = .02;
  // v8
  std::string hall_owner = "owner";            // owner (host arithmetic) | resident (field state on the GPU)
  double source_dt_max_s = 1e-14;              // source v7 trial interval bound (OwnerOptions::interval.cap_s)
  double source_emission_error = 0, source_desorption_error = 0, source_safety = .5, source_count_floor = 1e-3;
  std::string source_interval_preset = "benchmark";   // benchmark: no targets (v6.1); production: 1e-3 / 1e-3
  double neutral_dt_max_s = 1e-14;             // resident-neutral interval bound (Hall motion forms and native)
  std::string surface = "benchmark";           // source v7 surface options (README section 2)
  fs::path hall_checkpoint;                    // genesis at this Hall checkpoint (pre-run then genesis)
  fs::path checkpoint_dir;                     // whole-run checkpoints (default OUTDIR/checkpoints)
  fs::path resume, restore_checkpoint;         // restore: newest committed checkpoint below ROOT / this checkpoint
  long crash_step = -1;                        // test hook: crash (exit 86) while writing the checkpoint of this step
  int crash_stage = 0;                         // 0 owners written, 1 before manifest, 2 manifest temp, 3 before commit rename
  long ranks = 1;                              // recorded in checkpoints (not required at restore)
  // prerun
  double prerun_end_s = 0, prerun_sample_s = 1e-11, prerun_checkpoint_s = 1e-9;
  fs::path compare_hall_checkpoint;
};
Options opt;
v8::ControllerConfig ctl;
struct Backends {
  std::shared_ptr<hg::Device> device;
  std::shared_ptr<hc::FieldSolveBackend> solver;
  std::shared_ptr<hg::Device> gpu() {
    if (!device) device = std::make_shared<hg::Device>(hg::DeviceOptions{opt.uuid, opt.nvrtc});
    return device;
  }
} backends;
struct Paths {
  std::string trace, data, binding, outer3, kernel_binding;
} paths;
inline constexpr char outer3_sha[] = "c436020d3f2d07902fb059cfc0a269b2e914a3cbb2789300093d3631df896551";

struct Inputs {
  std::shared_ptr<const hc::Model> model;
  std::shared_ptr<const w::Data> data;
  std::shared_ptr<const w::ParentBridge> bridge;
  std::shared_ptr<const cv::HighEnergyBands> bands;
  std::shared_ptr<const cv::CoverageScope> scope;
  nw::Binding binding;
  v8::SourceSetup setup;
  v8::PilotSettings settings;
  v8::HallContext hall_context;
  s7::OwnerOptions source_options;
};
// Source v7 surface options of a preset. Only "benchmark" (v6.1 surface) and "work-function" (oxidized SS304
// with adsorbed-water coverage; prepare_patch uses it, the coordinator's birth path is unchanged) run in stage 1;
// enhancement classes and wall responses need the integration items of source v7 INTEGRATION.txt (D4 per-class
// streams and sample_birth births; B1-B4 native wall-response births; C1-C4 Hall wall heat), and are refused here.
sv7::SurfaceOptions surface_options(const std::string& preset) {
  sv7::Registry R;
  if (preset == "benchmark") return {};
  if (preset == "work-function") { sv7::SurfaceOptions s; s.work_function = sv7::ss304_work_function(R); return s; }
  if (preset == "work-function+enhancement" || preset == "full")
    throw std::runtime_error("surface preset " + preset + " needs per-class arrival streams and sample_birth births (source v7 INTEGRATION D4)" +
                             (preset == "full" ? std::string(" and native/Hall wall-response births and wall heat (B1-B4, C1-C4)") : std::string()) +
                             "; not enabled in coordinator v8 stage 1");
  throw std::runtime_error("--surface benchmark|work-function|work-function+enhancement|full");
}
v8::PilotSettings settings_for(const std::string& fixture) {
  v8::PilotSettings S;
  if (fixture == "qual") {
    S.work_function_eV = 2.0;
    S.w0_electron = 2e-23;
    S.w0_neutral = 2;
  } else if (fixture == "dense") {
    S.work_function_eV = 2.0;
    S.w0_electron = 2e-23;
    S.w0_neutral = 2;
  } else if (fixture == "load") {
    const std::size_t n = std::size_t(opt.load_carriers);
    if (n < 2) throw std::runtime_error("load fixture needs --load-carriers >= 2");
    S.maximum_carriers = n + n / 5 + 10000;
    S.hall_memory_budget = std::uint64_t(S.maximum_carriers) * 4096 + (std::uint64_t(512) << 20);
    S.native_memory_budget = std::max<std::uint64_t>(S.native_memory_budget, std::uint64_t(S.maximum_carriers) * 2048 + (std::uint64_t(512) << 20));
    if (opt.gpu_capacity < S.maximum_carriers) opt.gpu_capacity = S.maximum_carriers;
  } else if (fixture != "pilot")
    throw std::runtime_error("fixture must be qual, dense, load or pilot");
  S.arrival_seed += std::uint64_t(opt.seed_offset);
  if (opt.phi > 0) S.work_function_eV = opt.phi;
  if (opt.w0e > 0) S.w0_electron = opt.w0e;
  if (opt.w0n > 0) S.w0_neutral = opt.w0n;
  S.faces_per_slice = opt.faces;
  S.surface_preset = opt.surface;
  S.surface = surface_options(opt.surface);
  return S;
}
J mesh_record;
hc::MeshSpec selected_mesh() {
  hc::MeshSpec spec{};
  mesh_record = {{"file", opt.mesh_file.empty() ? J(nullptr) : J(opt.mesh_file)}};
  if (opt.mesh_file.empty()) return spec;
  const auto text = read_bytes(opt.mesh_file);
  const J mj = J::parse(text);
  if (mj.at("schema") != "resolution-series-mesh-input-v1") throw std::runtime_error("mesh file schema");
  const auto& hh = mj.at("hall_mesh_spec");
  spec = {hh.at("transverse").template get<int>(), hh.at("axial").template get<int>(), hh.at("air").template get<int>(), hh.at("width_m").template get<double>()};
  mesh_record["file_sha256"] = nw::sha(text);
  mesh_record["hall_mesh_spec"] = hh;
  mesh_record["encode_mesh_sha256"] = mj.value("encode_mesh_sha256", "");
  return spec;
}
void check_mesh(const axial::Mesh& g) {
  if (opt.mesh_file.empty()) return;
  const J mj = J::parse(read_bytes(opt.mesh_file));
  const auto& c = mj.at("counts");
  const auto bx = [](double x) { std::ostringstream o; o << std::hex << std::setw(16) << std::setfill('0') << std::bit_cast<std::uint64_t>(x); return o.str(); };
  const auto& sp = mj.at("spacing");
  const bool ok = c.at("nx") == g.nx && c.at("ny") == g.ny && c.at("nz") == g.nz && c.at("nodes") == g.nn && c.at("edges") == g.ne &&
                  c.at("faces") == g.nb && c.at("cells") == long(g.cells.size()) && sp.at("dx_m").at("bits") == bx(g.dx) &&
                  sp.at("dy_m").at("bits") == bx(g.dy) && sp.at("dz_m").at("bits") == bx(g.dz) && sp.at("zmin_m").at("bits") == bx(g.zmin);
  if (!ok) throw std::runtime_error("Hall mesh differs from the mesh file (counts or spacing bits)");
  mesh_record["checked_counts_and_spacing_bits"] = true;
}
// Field-solve backend of the run's Hall owner: host Owner: sealed CPU PCG (--field cpu) or gpu1 (--field gpu);
// ResidentOwner: the resident backend record (descriptor of the resident kernels on this device); its solve is
// the device CG (--field gpu) or the sealed CPU PCG on the downloaded right-hand side (--field cpu, the parity
// mode that matches the host Owner's cpu-pcg bit for bit).
std::shared_ptr<hc::FieldSolveBackend> solver(const hc::MeshSpec& spec) {
  if (!backends.solver) {
    if (opt.hall_owner == "resident") {
      const hc::Model probe(spec, paths.trace);
      backends.solver = std::make_shared<hc::ResidentBackend>(hc::resident_descriptor_sha256(backends.gpu(), probe));
    } else if (opt.field == "gpu") backends.solver = hg::make_gpu_field_solver(backends.gpu());
    else if (opt.field == "cpu") backends.solver = hc::make_cpu_pcg();
    else throw std::runtime_error("unknown field backend");
  }
  return backends.solver;
}
s7::OwnerOptions source_owner_options() {
  s7::OwnerOptions o;
  o.interval.cap_s = opt.source_dt_max_s;
  o.interval.emission_relative_error = opt.source_emission_error;
  o.interval.desorption_relative_error = opt.source_desorption_error;
  o.interval.safety = opt.source_safety;
  o.interval.count_floor = opt.source_count_floor;
  return o;
}
Inputs load_inputs(const std::string& fixture) {
  Inputs in;
  in.settings = settings_for(fixture);
  const auto spec = selected_mesh();
  const axial::Mesh g(spec.transverse, spec.axial, spec.width_m, spec.air);
  in.setup = v8::make_source_setup(g, in.settings);
  in.model = v8::bind::make_hall_model(opt.hall_dt_max_s, opt.neutral_dt_max_s, spec, paths.trace, solver(spec),
                                       hc::CarrierLimits{in.settings.maximum_carriers, in.settings.hall_memory_budget},
                                       hb::declaration(in.setup.slices));
  check_mesh(in.model->mesh());
  if (opt.hall_owner == "resident") {
    in.hall_context.device = backends.gpu();
    in.hall_context.resident.solve = opt.field == "gpu" ? hg::ResidentOptions::Solve::device_cg : hg::ResidentOptions::Solve::host_cpu_pcg;
    in.hall_context.resident.sample_seed = 20261031;   // host re-hash sample of digest leaves (no effect on any digest)
  }
  in.source_options = source_owner_options();
  std::ifstream f(paths.binding);
  J b;
  f >> b;
  if (!f) throw std::runtime_error("binding read");
  in.binding = {b.at("dataset_sha256"), b.at("kernel_sha256"), b.at("kernel_binding_file_sha256")};
  in.data = std::make_shared<const w::Data>(w::Data::read(paths.data));
  in.bridge = std::make_shared<const w::ParentBridge>(in.data, w::parent_model, in.binding.dataset_sha256, in.binding.kernel_sha256, v8::water_mass_kg);
  in.bands = std::make_shared<const cv::HighEnergyBands>(cv::HighEnergyBands::read(paths.outer3, outer3_sha, *in.data, in.binding.dataset_sha256,
                                                                                  paths.kernel_binding, in.binding.kernel_sha256,
                                                                                  in.binding.kernel_binding_file_sha256, cv::ProvisionalOptions{}));
  in.scope = std::make_shared<const cv::CoverageScope>(in.data, in.binding.dataset_sha256, *in.bridge, in.bands, cv::CoveragePolicy{});
  return in;
}
#include "driver_v8_fixtures.inc"
v8::GenesisRequestV8 request(const Inputs& in, std::vector<m::FullTimeAuthoredPacket> packets) {
  v8::GenesisRequestV8 rq;
  rq.model = in.model;
  rq.options = v8::hall_geometry(in.model->mesh());
  rq.seed = in.model->dc_seed();
  rq.states = v8::pilot_states();
  rq.channels = v8::pilot_channels(in.scope->arena_binding(in.binding.kernel_binding_file_sha256));
  rq.packets = std::move(packets);
  rq.grid_counts = {45, 45, 5};
  if (opt.collision_grid == "mesh") rq.grid_counts = {std::uint32_t(in.model->mesh().nx), std::uint32_t(in.model->mesh().ny), 5};
  else if (opt.collision_grid != "fixed200") throw std::runtime_error("--collision-grid mesh|fixed200");
  rq.binding = in.binding;
  rq.limits = {};
  if (opt.fixture == "load") {
    rq.limits.carriers = in.settings.maximum_carriers;
    rq.limits.pairs = rq.limits.queries = std::max<std::size_t>(rq.limits.pairs, 64 * in.settings.maximum_carriers);
    rq.limits.candidates = std::max<std::size_t>(rq.limits.candidates, in.settings.maximum_carriers);
  }
  rq.bridge = in.bridge;
  rq.scope = in.scope;
  rq.clock_seed = 123 + std::uint64_t(opt.seed_offset);
  rq.product_seed = 456 + std::uint64_t(opt.seed_offset);
  rq.origin_provenance = "coordinator-v6/genesis-at-first-C19-knot";
  rq.conductor_declaration = "coordinator v6 genesis: DC conductor/PEC induced charge of the C19 first-knot seed";
  rq.exterior_declaration = "coordinator v6 genesis: DC exterior termination charge of the C19 first-knot seed";
  rq.settings = in.settings;
  rq.sources = in.setup;
  rq.hall_context = in.hall_context;
  rq.source_options = in.source_options;
  rq.neutral_dt_max_s = opt.neutral_dt_max_s;
  if (!opt.hall_checkpoint.empty()) {
    rq.hall_checkpoint = opt.hall_checkpoint;
    rq.origin_provenance = "coordinator-v8/genesis-at-hall-checkpoint";
  }
  return rq;
}
template <class H> void attach_motion(v8::OwnersV8<H>& o) {
  if (opt.push == "gpu") o.motion = std::make_shared<hg::GpuMotionProvider>(backends.gpu(), o.hall->motion_context(), opt.gpu_capacity);
  else if (opt.push != "cpu") throw std::runtime_error("unknown motion provider");
}
template <class H> void configure_owners(v8::OwnersV8<H>& o) {
  o.dt_cap_s = ctl.cap_s;
  o.rate_sample_dt_s = ctl.cap_s;
  o.interval_startup_dt_s = ctl.particle_cap_s;
  attach_motion(o);
}
template <class H> v8::OwnersV8<H> genesis(const Inputs& in, const std::string& fixture) {
  auto o = v8::GenesisV8<H>::admit(request(in, fixture_packets(fixture)));
  configure_owners(o);
  return o;
}

// ---------------------------------------------------------------- publication with probes
long long probe_alloc[2] = {-1, -1};
steady::time_point probe_time[2];
int probe_calls = 0;
void publish_probe(int k) {
  probe_alloc[k] = stage1_allocation_count();
  probe_time[k] = steady::now();
  ++probe_calls;
}
struct PublishResult {
  bool published = false;
  int refusal = 0;
  double can_publish_s = 0, swap_s = 0, post_s = 0, total_s = 0;
  long long swap_allocations = -1;
};
template <class H> PublishResult publish_measured(v8::OwnersV8<H>& o, v8::TrialV8<H>& t) {
  probe_calls = 0;
  probe_alloc[0] = probe_alloc[1] = -1;
  TA::set_publish_probe(t, &publish_probe);
  const auto t0 = steady::now();
  PublishResult r;
  r.published = v8::try_publish_v8(o, t);
  const auto t1 = steady::now();
  r.refusal = t.refusal_point();
  r.total_s = std::chrono::duration<double>(t1 - t0).count();
  if (r.published && probe_calls == 2) {
    r.can_publish_s = std::chrono::duration<double>(probe_time[0] - t0).count();
    r.swap_s = std::chrono::duration<double>(probe_time[1] - probe_time[0]).count();
    r.post_s = std::chrono::duration<double>(t1 - probe_time[1]).count();
    r.swap_allocations = probe_alloc[1] - probe_alloc[0];
  }
  return r;
}

// ---------------------------------------------------------------- accepted totals
using XS = v8::ExactWeightSum;
struct Totals {
  std::uint64_t n_e = 0, n_i = 0, n_n = 0;
  XS w_e, w_i, w_n, h_nuclei, o_nuclei;
  m::ArenaCount stock;
  long double kinetic_J = 0, reference_J = 0;
  hc::Words q_free, q_conductor, q_exterior;
  double q_gross = 0;
};
template <class H> void charge_totals(const v8::OwnersV8<H>& o, Totals& t) {
  const auto& q = *o.accepted().charge;
  for (const auto& x : q.particles.net_C) t.q_free = t.q_free + x;
  for (double x : q.particles.gross_C) t.q_gross += x;
  for (const auto& x : q.conductor_C) t.q_conductor = t.q_conductor + x;
  for (const auto& x : q.exterior_C) t.q_exterior = t.q_exterior + x;
}
template <class H> Totals totals_legacy(const v8::OwnersV8<H>& o) {
  Totals t;
  m::ArenaCount we, wi, wn, hn, on;
  const auto idx = v8::state_index(o.inputs.states);
  for (const auto& r : o.arena->rate_population(o.inputs.maximum_carriers)) {
    const auto& s = *idx.at(r.state);
    const double wgt = r.native.weight;
    if (r.state == "electron") { ++t.n_e; we = we.credited(wgt); }
    else if (r.state == "H2O+") { ++t.n_i; wi = wi.credited(wgt); }
    else { ++t.n_n; wn = wn.credited(wgt); }
    for (int k = 0; k < s.hydrogen; ++k) hn = hn.credited(wgt);
    for (int k = 0; k < s.oxygen; ++k) on = on.credited(wgt);
    t.kinetic_J += (long double)wgt * m::pic_kinetic_energy_J(r.native);
    t.reference_J += (long double)wgt * (s.internal_J + s.reference_J);
  }
  t.w_e.add(we); t.w_i.add(wi); t.w_n.add(wn); t.h_nuclei.add(hn); t.o_nuclei.add(on);
  for (std::size_t i = 0; i < o.sources.size(); ++i)
    for (std::size_t j = 0; j < o.slices[i].patches.size(); ++j) t.stock = v8::exact_add(t.stock, o.sources[i]->inventory(j));
  charge_totals(o, t);
  return t;
}
template <class H> Totals totals(const v8::OwnersV8<H>& o, std::span<const m::ArenaRateEntry> population) {
  Totals t;
  const auto idx = v8::state_index(o.inputs.states);
  for (const auto& r : population) {
    const auto& s = *idx.at(r.state);
    const double wgt = r.native.weight;
    if (r.state == "electron") { ++t.n_e; t.w_e.add(wgt); }
    else if (r.state == "H2O+") { ++t.n_i; t.w_i.add(wgt); }
    else { ++t.n_n; t.w_n.add(wgt); }
    for (int k = 0; k < s.hydrogen; ++k) t.h_nuclei.add(wgt);
    for (int k = 0; k < s.oxygen; ++k) t.o_nuclei.add(wgt);
    t.kinetic_J += (long double)wgt * m::pic_kinetic_energy_J(r.native);
    t.reference_J += (long double)wgt * (s.internal_J + s.reference_J);
  }
  for (std::size_t i = 0; i < o.sources.size(); ++i)
    for (std::size_t j = 0; j < o.slices[i].patches.size(); ++j) t.stock = v8::exact_add(t.stock, o.sources[i]->inventory(j));
  charge_totals(o, t);
  return t;
}
double jnum(const J& j, std::initializer_list<const char*> path) {
  const J* x = &j;
  for (const char* k : path) {
    if (!x->is_object() || !x->contains(k)) return std::numeric_limits<double>::quiet_NaN();
    x = &x->at(k);
  }
  if (x->is_number()) return x->get<double>();
  if (x->is_object() && x->contains("value") && x->at("value").is_number()) return x->at("value").template get<double>();
  return std::numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------------- per-step ledger
struct Ledger {
  std::vector<std::string> names;
  std::vector<double> row;
  bool first = true;
  std::size_t at = 0;
  void begin() { row.clear(); at = 0; }
  void put(const char* name, double v) {
    if (first) names.push_back(name);
    else if (at >= names.size() || names[at] != name) throw std::runtime_error(std::string("ledger schema drift at ") + name);
    row.push_back(v);
    ++at;
  }
  void put(const char* name, long double v) { put(name, double(v)); }
  void put(const char* name, std::uint64_t v) { put(name, double(v)); }
  void put(const char* name, long v) { put(name, double(v)); }
  void put(const char* name, int v) { put(name, double(v)); }
  void put(const char* name, long long v) { put(name, double(v)); }
  void put(const char* name, bool v) { put(name, v ? 1. : 0.); }
};
bool deterministic_column(const std::string& n) { return n.rfind("t_", 0) != 0 && n != "alloc_swap_window" && n != "rss_kib"; }

struct RunState {
  Totals genesis;
  XS water_initial;
  m::ArenaCount sunk_h2o, sunk_ion;
  hc::Words q_initial;
  long double reaction_defect_cum = 0;
  double worst_step_residual = 0, worst_gauss_ratio = 0, worst_divB_ratio = 0, worst_q_change = 0, worst_reaction_defect = 0;
  bool water_exact = true, nuclei_exact = true;
  long electron_births = 0, neutral_births = 0, sinks = 0, events = 0;
};
template <class H>
void fill_row(Ledger& L, const v8::OwnersV8<H>& o, const v8::TrialV8<H>& t, const J& pre, const PublishResult& pub, const Totals& before,
              const Totals& after, RunState& st, long step, double prepare_s) {
  const auto& n = t.numbers();
  const auto& tm = t.times();
  L.begin();
  L.put("step", step);
  L.put("time_s", o.hall->time_s());
  L.put("dt_s", pre.at("dt").template get<double>());
  L.put("current_A", pre.at("current_A").template get<double>());
  L.put("terminal_V", pre.at("terminal_V").template get<double>());
  L.put("carriers_electron", after.n_e);
  L.put("carriers_ion", after.n_i);
  L.put("carriers_h2o", after.n_n);
  L.put("weight_electron", after.w_e.value());
  L.put("weight_ion", after.w_i.value());
  L.put("weight_h2o", after.w_n.value());
  L.put("births_electron", n.electron_births);
  L.put("births_h2o", n.neutral_births);
  L.put("sinks", n.sinks);
  L.put("parent_events", n.parent_events);
  const auto q_total = after.q_free + after.q_conductor + after.q_exterior;
  const double q_change = (q_total - st.q_initial).rounded();
  L.put("q_free_net_C", after.q_free.rounded());
  L.put("q_free_gross_C", after.q_gross);
  L.put("q_conductor_bank_C", after.q_conductor.rounded());
  L.put("q_exterior_bank_C", after.q_exterior.rounded());
  L.put("q_total_C", q_total.rounded());
  L.put("q_total_change_since_genesis_C", q_change);
  for (const char* k : {"q_collected_step_C", "q_escaped_step_C", "gauss_max_abs_C", "gauss_bound_C", "divB_max_T_m", "divB_bound_T_m",
                        "gauss_increment_C", "motion_continuity_C"})
    L.put(k, pre.at(k).template get<double>());
  const auto& kl = o.hall->kinetic_ledgers();
  L.put("q_collected_cum_C", kl.sinks.collected_charge_C().rounded());
  L.put("q_escaped_cum_C", kl.sinks.escaped_charge_C().rounded());
  st.sunk_h2o = v8::exact_add(st.sunk_h2o, n.sunk_h2o);
  st.sunk_ion = v8::exact_add(st.sunk_ion, n.sunk_h2o_ion);
  XS water = after.w_n;
  water.add(after.w_i); water.add(after.stock); water.add(st.sunk_h2o); water.add(st.sunk_ion);
  const bool water_ok = water == st.water_initial;
  st.water_exact = st.water_exact && water_ok;
  L.put("h2o_source_stock", v8::exact_value(after.stock));
  L.put("h2o_native_neutral", after.w_n.value());
  L.put("h2o_native_ion", after.w_i.value());
  L.put("h2o_sunk_neutral_cum", v8::exact_value(st.sunk_h2o));
  L.put("h2o_sunk_ion_cum", v8::exact_value(st.sunk_ion));
  L.put("h2o_total", water.value());
  L.put("h2o_total_exactly_conserved", water_ok);
  XS Hn = after.h_nuclei, On = after.o_nuclei;
  for (int k = 0; k < 2; ++k) { Hn.add(after.stock); Hn.add(st.sunk_h2o); Hn.add(st.sunk_ion); }
  On.add(after.stock); On.add(st.sunk_h2o); On.add(st.sunk_ion);
  XS H0 = st.genesis.h_nuclei, O0 = st.genesis.o_nuclei;
  for (int k = 0; k < 2; ++k) H0.add(st.genesis.stock);
  O0.add(st.genesis.stock);
  const bool nuclei_ok = Hn == H0 && On == O0;
  st.nuclei_exact = st.nuclei_exact && nuclei_ok;
  L.put("H_nuclei", Hn.value());
  L.put("O_nuclei", On.value());
  L.put("nuclei_exactly_conserved", nuclei_ok);
  for (const char* k : {"E_field_delta_J", "E_plasma_work_J", "E_joule_J", "E_outer_loss_J", "E_outer_supply_J", "E_circuit_work_J",
                        "E_kinetic_change_J", "E_kinetic_store_change_J", "E_thermal_delta_J", "E_wall_heat_J", "E_external_loads_J",
                        "E_collected_kinetic_J", "E_short_bridge_kinetic_J", "E_escaped_kinetic_J", "E_birth_kinetic_J", "E_combined_balance_J",
                        "E_step_residual_J", "E_step_scale_J", "E_roundoff_bound_J"})
    L.put(k, pre.at(k).template get<double>());
  L.put("E_source_neutral_kinetic_J", n.neutral_kinetic_J);
  L.put("E_source_neutral_internal_J", n.neutral_internal_J);
  L.put("E_source_adsorbed_reference_J", n.neutral_adsorbed_reference_J);
  L.put("E_source_loads_J", n.source_load_J);
  L.put("E_emission_nottingham_J", n.emission_nottingham_J);
  L.put("E_emission_reference_J", n.emission_reference_J);
  L.put("E_emission_matching_work_unrepresented_J", n.emission_matching_work_J);
  L.put("E_emission_kinetic_J", n.emission_kinetic_J);
  const long double dK = after.kinetic_J - before.kinetic_J, dR = after.reference_J - before.reference_J;
  const long double chem_dK = dK - (long double)pre.at("E_kinetic_store_change_J").template get<double>();
  const long double defect = chem_dK + dR;
  st.reaction_defect_cum += defect;
  L.put("E_census_kinetic_J", after.kinetic_J);
  L.put("E_census_reference_J", after.reference_J);
  L.put("E_chemistry_kinetic_delta_J", chem_dK);
  L.put("E_reaction_reference_delta_J", dR);
  L.put("E_reaction_defect_J", defect);
  L.put("E_reaction_defect_cum_J", st.reaction_defect_cum);
  const auto& cov = t.chemistry_state().coverage();
  L.put("cov_unresolved_pair_evaluations_step", jnum(cov, {"last_step", "totals", "unresolved", "pair_evaluations"}));
  L.put("cov_unresolved_upper_events_step", jnum(cov, {"last_step", "totals", "unresolved", "upper_bound_events"}));
  L.put("cov_unresolved_upper_physical_step", jnum(cov, {"last_step", "totals", "unresolved", "upper_bound_physical_reactions"}));
  L.put("cov_applied_expected_events_step", jnum(cov, {"last_step", "totals", "applied", "expected_events"}));
  L.put("cov_applied_expected_physical_step", jnum(cov, {"last_step", "totals", "applied", "expected_physical_reactions"}));
  double cats[3][2] = {{0, 0}, {0, 0}, {0, 0}};
  if (cov.contains("last_step") && cov.at("last_step").at("totals").contains("categories"))
    for (const auto& c : cov.at("last_step").at("totals").at("categories")) {
      const auto name = c.at("category").template get<std::string>();
      const int i = name.find("song") != std::string::npos ? 0 : name.find("outer3") != std::string::npos ? 1 : 2;
      cats[i][0] += c.value("executed_events", 0.);
      cats[i][1] += jnum(c, {"applied", "expected_events"});
    }
  L.put("cov_song_executed_events", cats[0][0]);
  L.put("cov_song_expected_events", cats[0][1]);
  L.put("cov_outer3_executed_events", cats[1][0]);
  L.put("cov_outer3_expected_events", cats[1][1]);
  L.put("cov_provisional_rbeb_executed_events", cats[2][0]);
  L.put("cov_provisional_rbeb_expected_events", cats[2][1]);
  L.put("cov_unresolved_upper_events_cum", jnum(cov, {"cumulative", "unresolved", "upper_bound_events"}));
  L.put("cov_applied_expected_events_cum", jnum(cov, {"cumulative", "applied", "expected_events"}));
  L.put("arr_expected_electrons_step", n.expected_electrons);
  L.put("arr_max_log10_expected_electrons_patch", n.max_log10_expected_electrons);
  L.put("arr_unresolved_emission_patches_step", n.unresolved_emission_patches);
  L.put("arr_unresolved_emission_lower_electrons_step", n.unresolved_emission_lower_electrons);
  L.put("arr_unresolved_emission_upper_electrons_step", n.unresolved_emission_upper_electrons);
  L.put("arr_expected_h2o_molecules_step", n.expected_neutral_molecules);
  L.put("arr_released_h2o_molecules_step", v8::exact_value(n.source_debit));
  std::uint64_t ce = 0, cn = 0;
  long double ee = 0, en = 0, var = 0, se_ = 0, sn = 0;
  for (const auto& x : o.arrival->totals()) {
    const bool electron = x.stream_id < v8::neutral_stream_base;
    (electron ? ce : cn) += x.births;
    (electron ? ee : en) += v8::exact_value(x.expected_sum);
    (electron ? se_ : sn) += v8::exact_value(x.sampled_sum);
    if (electron) var += x.variance_sum;
  }
  L.put("arr_cum_electron_births", ce);
  L.put("arr_cum_h2o_births", cn);
  L.put("arr_cum_expected_electrons", ee);
  L.put("arr_cum_sampled_electrons", se_);
  L.put("arr_cum_expected_h2o_molecules", en);
  L.put("arr_cum_released_h2o_molecules", sn);
  L.put("arr_electron_standardized_deviation", var > 0 ? double((se_ - ee) / ((long double)o.settings.w0_electron * std::sqrt(var))) : 0.);
  L.put("src_patch_T_max_K", n.max_patch_T_K);
  L.put("src_patch_T_min_K", n.min_patch_T_K);
  L.put("src_extraction_field_max_V_m", n.max_extraction_field_V_m);
  L.put("src_extraction_field_min_V_m", n.min_extraction_field_V_m);
  L.put("hall_iterations", int(tm.hall_iterations));
  L.put("solve_iterations", pre.at("solve_iterations").template get<double>());
  L.put("solve_residual_over_target", pre.at("solve_residual_over_target").template get<double>());
  L.put("plasma_dt", pre.at("plasma_dt").template get<double>());
  L.put("max_gyro", pre.at("max_gyro").template get<double>());
  L.put("max_displacement_cells", pre.at("max_displacement_cells").template get<double>());
  L.put("census_rows_step", n.census);
  L.put("newborns_step", n.newborns);
  const auto& hl = tm.hall;
  L.put("t_observe_ms", 1e3 * tm.observe_s);
  L.put("t_native_begin_ms", 1e3 * tm.begin_s);
  L.put("t_first_exposure_ms", 1e3 * tm.first_exposure_s);
  L.put("t_source_field_ms", 1e3 * tm.source_field_s);
  L.put("t_source_quotes_ms", 1e3 * tm.source_quotes_s);
  L.put("t_emission_exposure_ms", 1e3 * tm.emission_exposure_s);
  L.put("t_arrival_plan_ms", 1e3 * tm.arrival_plan_s);
  L.put("t_emission_births_ms", 1e3 * tm.emission_births_s);
  L.put("t_source_trials_ms", 1e3 * tm.source_trials_s);
  L.put("t_birth_credit_ms", 1e3 * tm.birth_credit_s);
  L.put("t_census_ms", 1e3 * tm.census_s);
  L.put("t_hall_admit_census_ms", 1e3 * hl.admit_census_s);
  L.put("t_hall_motion_ms", 1e3 * hl.motion_s);
  L.put("t_hall_deposit_support_ms", 1e3 * hl.support_s);
  L.put("t_hall_field_prepare_ms", 1e3 * hl.field_prepare_s);
  L.put("t_hall_field_solve_ms", 1e3 * hl.solve_s);
  L.put("t_hall_residual_check_ms", 1e3 * hl.residual_check_s);
  L.put("t_hall_v2_seal_ms", 1e3 * hl.motion_v2_seal_s);
  L.put("t_hall_field_seal_ms", 1e3 * hl.seal_s);
  L.put("t_hall_heat_ms", 1e3 * hl.heat_s);
  L.put("t_hall_proposal_ms", 1e3 * tm.hall_proposal_s);
  L.put("t_binding_checks_ms", 1e3 * tm.binding_checks_s);
  L.put("t_native_motion_admission_ms", 1e3 * tm.motion_admission_s);
  L.put("t_second_exposure_ms", 1e3 * tm.second_exposure_s);
  L.put("t_seal_ms", 1e3 * tm.seal_s);
  L.put("t_prepare_total_ms", 1e3 * prepare_s);
  L.put("t_can_publish_ms", 1e3 * pub.can_publish_s);
  L.put("t_swap_window_ms", 1e3 * pub.swap_s);
  L.put("t_post_commit_ms", 1e3 * pub.post_s);
  L.put("t_step_total_ms", 1e3 * (prepare_s + pub.total_s));
  L.put("alloc_swap_window", pub.swap_allocations);
  L.put("rss_kib", rss_kib());
  st.worst_step_residual = std::max(st.worst_step_residual, std::abs(pre.at("E_step_residual_J").template get<double>()));
  const double gb = pre.at("gauss_bound_C").template get<double>(), db = pre.at("divB_bound_T_m").template get<double>();
  if (gb > 0) st.worst_gauss_ratio = std::max(st.worst_gauss_ratio, pre.at("gauss_max_abs_C").template get<double>() / gb);
  if (db > 0) st.worst_divB_ratio = std::max(st.worst_divB_ratio, pre.at("divB_max_T_m").template get<double>() / db);
  st.worst_q_change = std::max(st.worst_q_change, std::abs(q_change));
  st.worst_reaction_defect = std::max(st.worst_reaction_defect, double(std::abs(defect)));
  st.electron_births += long(n.electron_births);
  st.neutral_births += long(n.neutral_births);
  st.sinks += long(n.sinks);
  st.events += long(n.parent_events);
}
template <class H> J pre_numbers(const v8::TrialV8<H>& t) {
  const double dt = t.dt_s();
  const auto& p = t.proposal();
  const auto& fh = p.field_heat();
  const auto& r = fh.field_receipt();
  const auto& c = fh.candidate_constraints();
  const auto& b = p.motion().budget();
  const auto& e = p.energy();
  const auto sinks = hc::SinkLedger::of(p.motion().sinks());
  const auto& s = fh.solve_record();
  return {{"dt", dt}, {"current_A", r.current_A}, {"terminal_V", r.terminal_V},
          {"q_collected_step_C", sinks.collected_charge_C().rounded()}, {"q_escaped_step_C", sinks.escaped_charge_C().rounded()},
          {"gauss_max_abs_C", c.max_gauss_C}, {"gauss_bound_C", c.max_gauss_bound_C}, {"divB_max_T_m", c.max_divB_T_m},
          {"divB_bound_T_m", c.max_divB_bound_T_m}, {"gauss_increment_C", r.gauss_increment_C}, {"motion_continuity_C", b.continuity_C},
          {"E_field_delta_J", e.field_delta_J}, {"E_plasma_work_J", e.plasma_work_J}, {"E_joule_J", e.joule_J}, {"E_outer_loss_J", e.outer_loss_J},
          {"E_outer_supply_J", e.outer_supply_J}, {"E_circuit_work_J", e.circuit_work_J}, {"E_kinetic_change_J", e.kinetic_change_J},
          {"E_kinetic_store_change_J", e.plasma_kinetic_store_change_J}, {"E_thermal_delta_J", e.thermal_delta_J}, {"E_wall_heat_J", e.wall_heat_J},
          {"E_external_loads_J", e.external_J}, {"E_collected_kinetic_J", e.collected_kinetic_J}, {"E_short_bridge_kinetic_J", e.short_bridge_kinetic_J},
          {"E_escaped_kinetic_J", e.escaped_kinetic_J}, {"E_birth_kinetic_J", e.birth_kinetic_J}, {"E_combined_balance_J", e.combined_balance_J},
          {"E_step_residual_J", e.residual_J}, {"E_step_scale_J", e.scale_J}, {"E_roundoff_bound_J", fh.heat().roundoff_bound_J},
          {"solve_iterations", double(s.iterations)}, {"solve_residual_over_target", s.target > 0 ? s.independent_residual / s.target : 0.},
          {"plasma_dt", b.plasma_dt}, {"max_gyro", b.maximum_gyro}, {"max_displacement_cells", b.maximum_displacement_cells}};
}
const std::vector<std::string> subphase_columns = {
    "step", "t_start_field_digest_ms", "t_source_field_build_ms", "t_emission_extracting_ms", "t_emission_retarding_ms",
    "t_emission_other_ms", "extracting_patches", "retarding_patches", "unresolved_patches", "t_cp_native_ms", "t_cp_chemistry_ms",
    "t_cp_observe_ms", "t_cp_predicate_ms", "t_cp_sources_ms", "t_cp_arrival_ms", "t_cp_hall_ms", "cp_source_preflights"};
template <class H> std::vector<double> subphase_row(const v8::TrialV8<H>& t, long step) {
  const auto& x = t.sub_times();
  return {double(step), 1e3 * x.start_field_digest_s, 1e3 * x.source_field_build_s, 1e3 * x.emission_extracting_s, 1e3 * x.emission_retarding_s,
          1e3 * x.emission_other_s, double(x.extracting_patches), double(x.retarding_patches), double(x.unresolved_patches), 1e3 * x.cp_native_s,
          1e3 * x.cp_chemistry_s, 1e3 * x.cp_observe_s, 1e3 * x.cp_predicate_s, 1e3 * x.cp_sources_s, 1e3 * x.cp_arrival_s, 1e3 * x.cp_hall_s,
          double(x.cp_source_preflights)};
}
J emission_statistics_json() {
  const auto e = emi::statistics();
  return {{"calls", e.calls}, {"sealed_calls", e.sealed_calls}, {"retarding_reused", e.retarding_reused}, {"retarding_sealed", e.retarding_sealed},
          {"extracting_fast", e.extracting_fast}, {"nodes", e.nodes}, {"nodes_exact", e.nodes_exact}, {"nodes_bounded_zero", e.nodes_bounded_zero},
          {"self_checks", e.self_checks}, {"sealed_fallback_configs", e.sealed_fallback_configs}};
}
template <class H> std::string digest_line(const v8::OwnersV8<H>& o, long step, const std::string& rng) {
  std::ostringstream s;
  s << step << ' ' << std::hexfloat << o.hall->time_s() << ' ' << o.arena->state_prefix() << ' ' << o.hall->field_sha256() << ' '
    << o.hall->heat_sha256() << ' ' << o.hall->kinetic_ledgers_sha256() << ' ' << o.accepted_source_seal << ' ' << rng
    << ' ' << o.accepted_arrival_state << '\n';
  return s.str();
}

// ---------------------------------------------------------------- frames (midplane maps)
struct FrameWriter {
  std::ofstream bin, index;
  long frames = 0;
  J layout;
  void open(const fs::path& dir, const hc::Mesh& g) {
    bin.open(dir / "frames.bin", std::ios::binary);
    index.open(dir / "frames.jsonl");
    layout = {{"schema", "coordinator-v6-midplane-frames-v1"}, {"byte_order", "little"}, {"value", "float64"}, {"nx", g.nx}, {"ny", g.ny},
              {"dx_m", g.dx}, {"dy_m", g.dy}, {"slab_m", g.dz}, {"midplane_z_m", 0.},
              {"record", {"step", "time_s", "n_e_cm3[nx*ny]", "n_ion_cm3[nx*ny]", "n_h2o_cm3[nx*ny]", "mean_electron_energy_eV[nx*ny]",
                          "Ex_V_m[nx*ny]", "Ey_V_m[nx*ny]", "phi_V[(nx+1)*(ny+1)]"}},
              {"cell_order", "index i + nx*j, i along native x (cell centre (i+1/2)dx), j along native y"},
              {"density", "cloud-in-cell in x,y over the slab |z| <= dz/2 around the z=0 node plane; carrier weights / (dx dy dz), cm^-3"},
              {"field", "accepted Hall E on the z=0 node plane, edge values averaged to cell centres"},
              {"potential", "phi(i,0)=0 on the y=0 boundary; phi(i,j+1)=phi(i,j)-Ey(i,j+1/2) dy along node lines (quasi-static gap potential)"}};
  }
  template <class H> void write(const v8::OwnersV8<H>& o, long step, std::span<const m::ArenaRateEntry> population) {
    const auto& g = o.hall->mesh();
    const int nx = g.nx, ny = g.ny;
    std::vector<double> ne(nx * ny), ni(nx * ny), nn(nx * ny), ke(nx * ny), Ex(nx * ny), Ey(nx * ny), phi((nx + 1) * (ny + 1));
    std::vector<double> we(nx * ny);
    std::vector<double> w3(std::size_t(nx) * ny * g.nz), k3(std::size_t(nx) * ny * g.nz);
    const double vol = g.dx * g.dy * g.dz;
    double total_e = 0, total_ke = 0;
    for (const auto& r : population) {
      const auto& p = r.native.position_m;
      const double K_eV = m::pic_kinetic_energy_J(r.native) / pc::qe;
      if (r.state == "electron") {
        const int i = std::clamp(int(p[0] / g.dx), 0, nx - 1), j = std::clamp(int(p[1] / g.dy), 0, ny - 1),
                  k = std::clamp(int((p[2] - g.zmin) / g.dz), 0, g.nz - 1);
        const std::size_t c = std::size_t(i) + std::size_t(nx) * (std::size_t(j) + std::size_t(ny) * std::size_t(k));
        w3[c] += r.native.weight;
        k3[c] += r.native.weight * K_eV;
        total_e += r.native.weight;
        total_ke += r.native.weight * K_eV;
      }
      if (std::abs(p[2]) > .5 * g.dz) continue;
      const double fx = p[0] / g.dx - .5, fy = p[1] / g.dy - .5;
      const int i0 = int(std::floor(fx)), j0 = int(std::floor(fy));
      const double ax = fx - i0, ay = fy - j0;
      for (int di = 0; di < 2; ++di)
        for (int dj = 0; dj < 2; ++dj) {
          const int i = std::clamp(i0 + di, 0, nx - 1), j = std::clamp(j0 + dj, 0, ny - 1);
          const double f = (di ? ax : 1 - ax) * (dj ? ay : 1 - ay) * r.native.weight;
          const std::size_t c = std::size_t(i) + std::size_t(nx) * std::size_t(j);
          if (r.state == "electron") { ne[c] += f; we[c] += f; ke[c] += f * K_eV; }
          else if (r.state == "H2O+") ni[c] += f;
          else nn[c] += f;
        }
    }
    for (std::size_t c = 0; c < ne.size(); ++c) {
      ke[c] = we[c] > 0 ? ke[c] / we[c] : 0;
      ne[c] *= 1e-6 / vol; ni[c] *= 1e-6 / vol; nn[c] *= 1e-6 / vol;
    }
    const auto& E = *o.accepted().E;
    const int k0 = g.nz / 2;
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        Ex[i + nx * j] = .5 * (E[g.ex(i, j, k0)] + E[g.ex(i, j + 1, k0)]);
        Ey[i + nx * j] = .5 * (E[g.ey(i, j, k0)] + E[g.ey(i + 1, j, k0)]);
      }
    for (int i = 0; i <= nx; ++i)
      for (int j = 0; j < ny; ++j) phi[i + (nx + 1) * (j + 1)] = phi[i + (nx + 1) * j] - E[g.ey(i, j, k0)] * g.dy;
    double ne_max3 = 0, te_at = 0;
    for (std::size_t c = 0; c < w3.size(); ++c)
      if (w3[c] > 0 && w3[c] / vol > ne_max3) { ne_max3 = w3[c] / vol; te_at = (2. / 3.) * k3[c] / w3[c]; }
    double ne_max_slab = 0;
    for (double x : ne) ne_max_slab = std::max(ne_max_slab, x);
    const double lambda = ne_max3 > 0 && te_at > 0 ? std::sqrt(8.8541878128e-12 * te_at / (ne_max3 * pc::qe)) : INFINITY;
    double emax = 0;
    for (std::size_t c = 0; c < Ex.size(); ++c) emax = std::max(emax, std::hypot(Ex[c], Ey[c]));
    const double header[2] = {double(step), o.hall->time_s()};
    bin.write(reinterpret_cast<const char*>(header), sizeof header);
    for (const auto* v : {&ne, &ni, &nn, &ke, &Ex, &Ey, &phi}) bin.write(reinterpret_cast<const char*>(v->data()), std::streamsize(v->size() * sizeof(double)));
    bin.flush();
    J row = {{"frame", frames}, {"step", step}, {"time_s", o.hall->time_s()}, {"ne_max_slab_cm3", ne_max_slab},
             {"ne_max_cell_m3", ne_max3}, {"Te_at_max_eV", te_at}, {"debye_length_m", std::isfinite(lambda) ? J(lambda) : J("inf")},
             {"cell_size_m", g.dx}, {"debye_over_cell", std::isfinite(lambda) ? J(lambda / g.dx) : J("inf")},
             {"electron_weight_total", total_e}, {"electron_mean_energy_eV", total_e > 0 ? total_ke / total_e : 0.},
             {"max_midplane_E_V_m", emax}, {"gap_voltage_midline_V", phi[nx / 2 + (nx + 1) * ny]}};
    index << row.dump() << '\n';
    index.flush();
    ++frames;
  }
};

// ---------------------------------------------------------------- v7 records
double chemistry_rate(const J& cov) {
  double rate = 0;
  if (!cov.contains("last_step") || !cov.at("last_step").contains("exposures")) return 0;
  for (const auto& x : cov.at("last_step").at("exposures")) {
    const double d = x.value("duration_s", 0.);
    const auto& lim = x.at("exposure_limit").at("events_upper");
    const double up = lim.is_number() ? lim.template get<double>() : INFINITY;
    if (d > 0 && up > 0) rate = std::max(rate, up / d);
  }
  return rate;
}
J criterion_codes() {
  J c = J::object();
  // v7.1 codes 0..14; the source_interval code (15) is listed only for runs with the criterion active
  const int last = ctl.source_interval ? int(v8::Criterion::source_interval) : int(v8::Criterion::retry);
  for (int k = 0; k <= last; ++k) c[v8::criterion_name(v8::Criterion(k))] = k;
  return c;
}
void control_row(Ledger& C, long step, double t_begin, double t_end, const v8::Decision& d, const v8::Decision& first, const v8::ControlInputs& in,
                 int attempts, const J& pre, double control_s) {
  constexpr double c = 299792458.;
  C.begin();
  C.put("step", step);
  C.put("time_s", t_end);
  C.put("t_begin_s", t_begin);
  C.put("dt_s", d.dt_s);
  C.put("rung", d.rung);
  C.put("landing", d.landing);
  C.put("adaptive", ctl.adaptive);
  C.put("criterion", int(d.criterion));
  C.put("criterion_value", d.value);
  C.put("criterion_limit", d.limit);
  C.put("first_dt_s", first.dt_s);
  C.put("first_criterion", int(first.criterion));
  C.put("attempts", attempts);
  C.put("carriers", in.carriers);
  C.put("charged_carriers", in.charged);
  C.put("electron_carriers", in.electrons);
  C.put("omega_p_max_rad_s", in.omega_p_max);
  C.put("omega_p_dt", in.omega_p_max * d.dt_s);
  C.put("omega_c_max_rad_s", in.omega_c_max);
  C.put("omega_c_dt", in.omega_c_max * d.dt_s);
  C.put("v_max_m_s", in.v_max);
  C.put("a_max_m_s2", in.a_max);
  C.put("displacement_pred_cells", d.dt_s * std::min(c, in.v_max + in.a_max * d.dt_s) / in.h_min);
  C.put("chem_rate_per_s", in.chem_rate);
  C.put("chem_pred_events_exposure", in.chem_rate * .5 * d.dt_s);
  C.put("chem_limit_events_exposure", in.chem_limit);
  C.put("drive_I_A", in.drive_I);
  C.put("lte_valid", in.lte_valid);
  C.put("lte_coeff_per_s3", in.lte_coeff);
  C.put("lte_pred_rel", in.lte_coeff * d.dt_s * d.dt_s * d.dt_s / 12);
  C.put("births_rate_per_s", in.births_rate);
  C.put("births_pred", in.births_rate * d.dt_s);
  C.put("bound_plasma_s", d.bound_plasma);
  C.put("bound_gyro_s", d.bound_gyro);
  C.put("bound_displacement_s", d.bound_disp);
  C.put("bound_chemistry_s", d.bound_chem);
  C.put("bound_drive_s", d.bound_drive);
  C.put("bound_lte_s", d.bound_lte);
  C.put("measured_plasma_dt", pre.at("plasma_dt").template get<double>());
  C.put("measured_gyro", pre.at("max_gyro").template get<double>());
  C.put("measured_displacement_cells", pre.at("max_displacement_cells").template get<double>());
  C.put("t_control_ms", 1e3 * control_s);
  if (ctl.source_interval) {   // v8 columns (interval policy only; the v7.1 layout otherwise)
    C.put("bound_source_s", d.bound_source);
    C.put("source_startup", in.source_startup);
    C.put("source_rate_per_s", in.source_rate_per_s);
  }
}
struct ProbeWriter {
  std::ofstream bin;
  std::vector<std::array<int, 6>> idx;
  J layout;
  template <class H> void open(const fs::path& dir, const v8::OwnersV8<H>& o) {
    const auto& g = o.hall->mesh();
    const std::vector<std::pair<std::string, std::array<double, 3>>> points = {
        {"gap-centre", {0, 0, 0}}, {"gap-lower-half", {0, -1.25e-3, 0}}, {"gap-upper-half", {0, 1.25e-3, 0}},
        {"near-lower-foil", {0, -2.3e-3, 0}}, {"near-upper-foil", {0, 2.3e-3, 0}}, {"beyond-foil-edge", {2.7e-3, 0, 0}},
        {"toward-feed", {0, 0, 9e-3}}, {"blob-centre", {-1.5e-3, 1.0e-3, 0}}, {"blob-edge", {-.9e-3, 1.0e-3, 0}},
        {"above-blob", {-1.5e-3, 2.0e-3, 0}}};
    J pts = J::array(), cols = J::array({"step", "time_s", "W_E_J", "W_B_J"});
    for (const auto& [name, x] : points) {
      const int i = std::clamp(int(std::floor((x[0] + g.side / 2) / g.dx)), 0, g.nx - 1);
      const int j = std::clamp(int(std::floor((x[1] + g.side / 2) / g.dy)), 0, g.ny - 1);
      const int k = std::clamp(int(std::floor((x[2] - g.zmin) / g.dz)), 0, g.nz - 1);
      idx.push_back({g.ex(i, j, k), g.ey(i, j, k), g.ez(i, j, k), g.bx(i, j, k), g.by(i, j, k), g.bz(i, j, k)});
      pts.push_back({{"name", name}, {"field_point_m", x}, {"cell", {i, j, k}}});
      for (const char* c : {"Ex", "Ey", "Ez", "Bx", "By", "Bz"}) cols.push_back(name + std::string(":") + c);
    }
    layout = {{"schema", "coordinator-v7-probes-v1"}, {"byte_order", "little"}, {"value", "float64"}, {"columns", cols}, {"points", pts},
              {"note", "accepted end-of-step E (edges ex/ey/ez(i,j,k)) and B (faces bx/by/bz(i,j,k)) of the cell containing the point; "
                       "W_E = sum_e M_e E_e^2/2, W_B = sum_b H_b B_b^2/2 in index order every --probe-energy-every steps (NaN otherwise); row 0 is genesis"}};
    bin.open(dir / "probes.bin", std::ios::binary);
    write_text(dir / "probes-schema.json", layout.dump(1) + "\n");
  }
  template <class H> void write(const v8::OwnersV8<H>& o, long step) {
    const auto& g = o.hall->mesh();
    const auto a = o.accepted();
    const auto& E = *a.E;
    const auto& B = *a.B;
    double we = NAN, wb = NAN;
    if (step == 0 || (opt.probe_energy_every > 0 && step % opt.probe_energy_every == 0)) {
      we = wb = 0;
      for (int e = 0; e < g.ne; ++e) we += .5 * g.M[std::size_t(e)] * E[std::size_t(e)] * E[std::size_t(e)];
      for (int b = 0; b < g.nb; ++b) wb += .5 * g.H[std::size_t(b)] * B[std::size_t(b)] * B[std::size_t(b)];
    }
    std::vector<double> row = {double(step), o.hall->time_s(), we, wb};
    for (const auto& x : idx)
      for (int c = 0; c < 6; ++c) row.push_back(c < 3 ? E[std::size_t(x[std::size_t(c)])] : B[std::size_t(x[std::size_t(c)])]);
    bin.write(reinterpret_cast<const char*>(row.data()), std::streamsize(row.size() * sizeof(double)));
  }
  void flush() { bin.flush(); }
};
struct CoordTimes {
  double control_s = 0, totals_s = 0, ledger_row_s = 0, ledger_io_s = 0, digest_line_s = 0, controller_state_s = 0, control_row_s = 0,
         probes_s = 0, receipt_s = 0, frames_s = 0, observe_prepare_s = 0, observe_preflight_s = 0, step_s = 0;
  double gpu_solve_s = 0, gpu_push_upload_s = 0, gpu_push_admit_s = 0, gpu_push_device_s = 0, gpu_push_assemble_s = 0, gpu_push_total_s = 0;
};
struct DeviceCounters {
  double solve = 0, upload = 0, admit = 0, device = 0, assemble = 0, total = 0;
};
DeviceCounters device_counters(const Inputs& in, const std::shared_ptr<hc::MotionProvider>& motion) {
  DeviceCounters c;
  if (opt.field == "gpu" && opt.hall_owner == "owner") c.solve = hg::gpu_field_solver_stats(in.model->solver()).total_s;
  if (const auto g = std::dynamic_pointer_cast<hg::GpuMotionProvider>(motion)) {
    const auto& s = g->stats();
    c.upload = s.upload_carriers_s; c.admit = s.admit_s; c.device = s.device_run_s; c.assemble = s.assemble_s; c.total = s.total_s;
  }
  return c;
}
const std::vector<std::string> coordinator_columns = {
    "step", "carriers", "t_control_inputs_and_decision_ms", "t_totals_ms", "t_ledger_row_ms", "t_ledger_io_ms", "t_digest_line_ms",
    "t_controller_state_ms", "t_control_row_ms", "t_probes_ms", "t_receipt_ms", "t_frames_ms", "t_observe_in_prepare_ms",
    "t_observe_in_preflight_ms", "t_step_total_ms", "t_device_field_solve_ms", "t_device_push_upload_ms", "t_device_push_admit_ms",
    "t_device_push_run_ms", "t_device_push_assemble_ms", "t_device_push_total_ms"};
std::vector<double> coordinator_row(const CoordTimes& c, long step, std::uint64_t carriers) {
  return {double(step), double(carriers), 1e3 * c.control_s, 1e3 * c.totals_s, 1e3 * c.ledger_row_s, 1e3 * c.ledger_io_s, 1e3 * c.digest_line_s,
          1e3 * c.controller_state_s, 1e3 * c.control_row_s, 1e3 * c.probes_s, 1e3 * c.receipt_s, 1e3 * c.frames_s, 1e3 * c.observe_prepare_s,
          1e3 * c.observe_preflight_s, 1e3 * c.step_s, 1e3 * c.gpu_solve_s, 1e3 * c.gpu_push_upload_s, 1e3 * c.gpu_push_admit_s,
          1e3 * c.gpu_push_device_s, 1e3 * c.gpu_push_assemble_s, 1e3 * c.gpu_push_total_s};
}

#include "driver_v8_run.inc"
#include "driver_v8_checks.inc"
}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 9)
      throw std::runtime_error("usage: driver_v8 run|checks|vacuum|prerun TRACE H2O_DATA BINDING OUTER3 KERNEL_BINDING OUTDIR NATIVE_IDENTITY [options]");
    const std::string mode = argv[1];
    paths = {argv[2], argv[3], argv[4], argv[5], argv[6]};
    const fs::path out = argv[7];
    if (std::string(m::PicArena::source_identity()) != argv[8]) throw std::runtime_error("linked native source identity differs");
    for (int i = 9; i < argc; ++i) {
      const std::string a = argv[i];
      auto value = [&] { if (i + 1 >= argc) throw std::runtime_error("missing option value"); return std::string(argv[++i]); };
      if (a == "--field") opt.field = value();
      else if (a == "--push") opt.push = value();
      else if (a == "--uuid") opt.uuid = value();
      else if (a == "--nvrtc") opt.nvrtc = value();
      else if (a == "--fixture") opt.fixture = value();
      else if (a == "--steps") { opt.steps = std::stol(value()); opt.steps_requested = true; }
      else if (a == "--end-time") { opt.end_time_s = std::stod(value()); opt.end_time_requested = true; }
      else if (a == "--dt") opt.dt_s = std::stod(value());
      else if (a == "--frame-every") opt.frame_every = std::stol(value());
      else if (a == "--parity-every") opt.parity_every = std::stol(value());
      else if (a == "--progress-every") opt.progress_every = std::stol(value());
      else if (a == "--checkpoint-every") opt.checkpoint_every = std::stol(value());
      else if (a == "--w0e") opt.w0e = std::stod(value());
      else if (a == "--w0n") opt.w0n = std::stod(value());
      else if (a == "--phi") opt.phi = std::stod(value());
      else if (a == "--faces") opt.faces = std::stoi(value());
      else if (a == "--gpu-capacity") opt.gpu_capacity = std::size_t(std::stod(value()));
      else if (a == "--stop-file") opt.stop_file = value();
      else if (a == "--dt-mode") { const auto v = value(); if (v != "fixed" && v != "adaptive") throw std::runtime_error("--dt-mode fixed|adaptive"); ctl.adaptive = v == "adaptive"; }
      else if (a == "--dt-cap") ctl.cap_s = std::stod(value());
      else if (a == "--hall-dt-max") opt.hall_dt_max_s = std::stod(value());
      else if (a == "--particle-cap") ctl.particle_cap_s = std::stod(value());
      else if (a == "--charged-only-particle-cap") ctl.neutrals_need_particle_cap = false;
      else if (a == "--c-plasma") ctl.c_plasma = std::stod(value());
      else if (a == "--c-gyro") ctl.c_gyro = std::stod(value());
      else if (a == "--c-disp") ctl.c_disp = std::stod(value());
      else if (a == "--c-chem") ctl.c_chem = std::stod(value());
      else if (a == "--c-drive") ctl.c_drive = std::stod(value());
      else if (a == "--lte-tol") ctl.lte_tol = std::stod(value());
      else if (a == "--birth-threshold") ctl.birth_threshold = std::stod(value());
      else if (a == "--growth-rungs") ctl.growth_rungs = std::stoi(value());
      else if (a == "--max-retries") ctl.max_retries = std::stoi(value());
      else if (a == "--sync-every") ctl.sync_every_s = std::stod(value());
      else if (a == "--window") { opt.window_s = std::stod(value()); opt.window_requested = true; }
      else if (a == "--fixed-grid-retry") ctl.fixed_grid_retry = true;
      else if (a == "--mesh") opt.mesh_file = value();
      else if (a == "--probe-energy-every") opt.probe_energy_every = std::stol(value());
      else if (a == "--host-work") { opt.host_work = value(); if (opt.host_work != "v7" && opt.host_work != "v7.1") throw std::runtime_error("--host-work v7|v7.1"); }
      else if (a == "--collision-grid") opt.collision_grid = value();
      else if (a == "--load-carriers") opt.load_carriers = std::stol(value());
      else if (a == "--load-neutral-fraction") opt.load_neutral_fraction = std::stod(value());
      else if (a == "--seed-offset") opt.seed_offset = std::stol(value());
      else if (a == "--blob-pairs") opt.blob_pairs = std::stol(value());
      else if (a == "--blob-weight") opt.blob_weight = std::stod(value());
      else if (a == "--blob-ev") opt.blob_eV = std::stod(value());
      else if (a == "--manufactured-refusal-step") opt.manufactured_refusal_step = std::stol(value());
      // v8
      else if (a == "--hall-owner") { opt.hall_owner = value(); if (opt.hall_owner != "owner" && opt.hall_owner != "resident") throw std::runtime_error("--hall-owner owner|resident"); }
      else if (a == "--source-dt-max") opt.source_dt_max_s = std::stod(value());
      else if (a == "--source-emission-error") opt.source_emission_error = std::stod(value());
      else if (a == "--source-desorption-error") opt.source_desorption_error = std::stod(value());
      else if (a == "--source-safety") opt.source_safety = std::stod(value());
      else if (a == "--source-count-floor") opt.source_count_floor = std::stod(value());
      else if (a == "--source-interval-preset") {
        opt.source_interval_preset = value();
        if (opt.source_interval_preset == "production") { opt.source_emission_error = 1e-3; opt.source_desorption_error = 1e-3; opt.source_safety = .5; opt.source_count_floor = 1e-3; }
        else if (opt.source_interval_preset == "benchmark") { opt.source_emission_error = 0; opt.source_desorption_error = 0; }
        else throw std::runtime_error("--source-interval-preset benchmark|production");
      }
      else if (a == "--neutral-dt-max") opt.neutral_dt_max_s = std::stod(value());
      else if (a == "--surface") opt.surface = value();
      else if (a == "--hall-checkpoint") opt.hall_checkpoint = value();
      else if (a == "--checkpoint-dir") opt.checkpoint_dir = value();
      else if (a == "--resume") opt.resume = value();
      else if (a == "--restore-checkpoint") opt.restore_checkpoint = value();
      else if (a == "--crash-in-checkpoint") { const auto v = value(); const auto c = v.find(':'); opt.crash_step = std::stol(v.substr(0, c)); opt.crash_stage = c == std::string::npos ? 0 : std::stoi(v.substr(c + 1)); }
      else if (a == "--ranks") opt.ranks = std::stol(value());
      else if (a == "--prerun-end") opt.prerun_end_s = std::stod(value());
      else if (a == "--prerun-sample") opt.prerun_sample_s = std::stod(value());
      else if (a == "--prerun-checkpoint-every") opt.prerun_checkpoint_s = std::stod(value());
      else if (a == "--compare-hall-checkpoint") opt.compare_hall_checkpoint = value();
      else throw std::runtime_error("unknown option " + a);
    }
    fs::create_directories(out);
    int rc = 2;
    if (mode == "run") rc = opt.hall_owner == "resident" ? run_mode<hc::ResidentOwner>(out) : run_mode<hc::Owner>(out);
    else if (mode == "checks") { finish_configuration(); rc = checks_mode(out); }
    else if (mode == "vacuum") { finish_configuration(); rc = vacuum_mode(out); }
    else if (mode == "prerun") { finish_configuration(); rc = opt.hall_owner == "resident" ? prerun_mode<hc::ResidentOwner>(out) : prerun_mode<hc::Owner>(out); }
    else if (mode == "hallcompare") { finish_configuration(); rc = opt.hall_owner == "resident" ? hallcompare_mode<hc::ResidentOwner>(out) : hallcompare_mode<hc::Owner>(out); }
    else throw std::runtime_error("unknown mode");
    if (backends.device) write_text(out / "gpu-teardown.json", backends.device->close_json() + "\n");
    return rc;
  } catch (const std::exception& e) {
    std::cerr << "fatal: " << e.what() << '\n';
    return 2;
  }
}
