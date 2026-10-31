// Coordinator v7 driver (the coordinator v6.1 driver with the v7 step controller,
// predicate v7, the dense calibration fixture, probes, the control ledger and the
// sealed step receipts). Modes:
//   run      genesis (fixture qual|pilot) + consecutive accepted steps with per-step
//            ledgers, deterministic digests, density/field frames and timings; the
//            ledger, digest, genesis and frame files have the v6 layout and schemas;
//            v6.1 sub-phase timings go to subphase-timings.bin (subphase-schema.json)
//   checks   T1/T2/T3/T4/T6-style qualification on the qual fixture
//   vacuum   Hall-only vacuum pre-run cost per step and the dt guard
//   profile  after --steps accepted steps: sealed v6 operations against their v6.1
//            replacements on the live state (times, bit identity), profile.json
// usage: driver_v7 MODE TRACE H2O_DATA BINDING OUTER3 KERNEL_BINDING OUTDIR NATIVE_IDENTITY [options]
#include "genesis_v7_1.hpp"
#include "step_controller_v7_1.hpp"
#include "exact_weight_sum_v7_1.hpp"
#include "gpu/gpu_backends.hpp"
#include <csignal>
#include <filesystem>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <random>

extern "C" {
long long stage1_allocation_count() noexcept;
void stage1_set_allocation_countdown(long long) noexcept;
}
namespace v7 = raikou::native_water::v7_1;
namespace emi = ::electrode_heat_v6_1::emission;
namespace v3k = ::electrode_heat_v3;
namespace se = ::surface_emission;
namespace m = raikou::mfem_backend;
namespace h = ::hall_driver;
namespace hc = ::hall_corrective_v3;
namespace hg = ::hall_corrective_v3::gpu;
namespace hb = ::hall_source_binding_v6_1;
namespace cv = ::raikou::native_water::coverage_v2;
namespace nw = ::raikou::native_water;
namespace pc = ::raikou::pic_collisions;
namespace w = ::raikou::pic_water;
namespace pa = ::pf_arrival_v11;
namespace s6 = ::electrode_heat_v6_1;
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

namespace raikou::native_water::v7_1 {
struct TrialTestAccessV7 {
  static void set_can_publish_fault(TrialV7& t, int k) { t.can_publish_fault_ = k; }
  static void set_publish_probe(TrialV7& t, void (*p)(int)) { t.publish_probe_ = p; }
  static std::string& source_after(TrialV7& t) { return t.source_after_; }
  static std::unique_ptr<hc::JointProposal>& proposal(TrialV7& t) { return t.field_heat_; }
  static std::unique_ptr<m::PreparedArena>& native(TrialV7& t) { return t.native_; }
  static std::vector<std::optional<s6::Trial>>& source_trials(TrialV7& t) { return t.source_trials_; }
  static std::optional<pa::Staged>& staged(TrialV7& t) { return t.staged_; }
  static c7::StepClock& clock(TrialV7& t) { return t.clock_; }
  static void compare_v6_predicate(TrialV7& t) { t.compare_v6_predicate_ = true; }
  static std::pair<int, int> predicates(const TrialV7& t) { return {t.v6_predicate_, t.v7_predicate_}; }
};
struct GenesisTestAccessV7 {
  static OwnersV7 admit(const GenesisRequestV7& rq, GenesisV7::PacketHook p, GenesisV7::LedgerHook l, int* reached) {
    return GenesisV7::admit_with_hooks(rq, p, l, reached);
  }
};
}  // namespace raikou::native_water::v7_1
using TA = v7::TrialTestAccessV7;

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

// ---------------------------------------------------------------- configuration
struct Options {
  std::string field = "cpu", push = "cpu", uuid, nvrtc, fixture = "pilot";
  long steps = 20, frame_every = 0, parity_every = 0, progress_every = 1000, checkpoint_every = 0;
  double end_time_s = 0, dt_s = 1e-14;
  double w0e = -1, w0n = -1, phi = -1;
  int faces = 4;
  std::size_t gpu_capacity = 65536;
  fs::path stop_file;
  // v7
  double window_s = 0;                       // end time = genesis time + window (alternative to --end-time)
  long seed_offset = 0;                      // added to the arrival, collision-clock and product seeds
  long blob_pairs = 400;                     // dense fixture: co-located electron/H2O+ pairs
  double blob_weight = 1.6e11, blob_eV = 3;  // weight per carrier; electron temperature (eV)
  long manufactured_refusal_step = 0;        // test hook: the first attempt of this step is refused as a step-size guard
  double hall_dt_max_s = 1e-14;              // configured step cap of the Hall model (Hall bindings with a configurable cap)
  long probe_energy_every = 1;
  std::string host_work = "v7.1";            // v7: the coordinator v7 per-step host work (profile baseline; same outputs)               // field energies W_E, W_B in probes.bin every K steps (O(mesh)); 0: never (NaN)
  // v7.1
  std::string mesh_file;                     // resolution-series mesh (meshes/mesh-*.json); empty: the 200 um mesh
  std::string collision_grid = "mesh";       // mesh: collision cells = transverse mesh cells x 5; fixed200: 45 x 45 x 5
  long load_carriers = 0;                    // load fixture: number of authored carriers
  double load_neutral_fraction = .02;        // load fixture: share of H2O carriers (weight 1e6); the rest e/H2O+ pairs (weight 1)
};
Options opt;
v7::ControllerConfig ctl;
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
  v7::SourceSetup setup;
  v7::PilotSettings settings;
};
v7::PilotSettings settings_for(const std::string& fixture) {
  v7::PilotSettings S;
  if (fixture == "qual") {
    // Qualification-only emission and release weights (births every step); not physical.
    S.work_function_eV = 2.0;   // lower end of the kernel's 2..6 eV domain
    S.w0_electron = 2e-23;
    S.w0_neutral = 2;
  } else if (fixture == "dense") {
    // Dense calibration fixture: the qualification fixture plus a seeded plasma blob.
    S.work_function_eV = 2.0;
    S.w0_electron = 2e-23;
    S.w0_neutral = 2;
  } else if (fixture == "load") {
    // Carrier-load fixture for host profiles: the pilot sources and opt.load_carriers authored carriers; caps sized to them.
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
  return S;
}
std::shared_ptr<hc::FieldSolveBackend> solver() {
  if (!backends.solver) {
    if (opt.field == "gpu") backends.solver = hg::make_gpu_field_solver(backends.gpu());
    else if (opt.field == "cpu") backends.solver = hc::make_cpu_pcg();
    else throw std::runtime_error("unknown field backend");
  }
  return backends.solver;
}
// Mesh selected at launch (v7.1): one resolution-series mesh file sets the Hall MeshSpec, the
// source-setup mesh and (with --collision-grid mesh) the collision grid; the built Hall mesh is
// checked against the file's counts and spacing bits.
J mesh_record;
hc::MeshSpec selected_mesh() {
  hc::MeshSpec spec{};
  mesh_record = {{"file", opt.mesh_file.empty() ? J(nullptr) : J(opt.mesh_file)}};
  if (opt.mesh_file.empty()) return spec;
  const auto text = read_bytes(opt.mesh_file);
  const J m = J::parse(text);
  if (m.at("schema") != "resolution-series-mesh-input-v1") throw std::runtime_error("mesh file schema");
  const auto& h = m.at("hall_mesh_spec");
  spec = {h.at("transverse").get<int>(), h.at("axial").get<int>(), h.at("air").get<int>(), h.at("width_m").get<double>()};
  mesh_record["file_sha256"] = nw::sha(text);
  mesh_record["hall_mesh_spec"] = h;
  mesh_record["encode_mesh_sha256"] = m.value("encode_mesh_sha256", "");
  return spec;
}
void check_mesh(const axial::Mesh& g) {
  if (opt.mesh_file.empty()) return;
  const J m = J::parse(read_bytes(opt.mesh_file));
  const auto& c = m.at("counts");
  const auto bits = [](double x) { std::ostringstream o; o << std::hex << std::setw(16) << std::setfill('0') << std::bit_cast<std::uint64_t>(x); return o.str(); };
  const auto& sp = m.at("spacing");
  const bool ok = c.at("nx") == g.nx && c.at("ny") == g.ny && c.at("nz") == g.nz && c.at("nodes") == g.nn && c.at("edges") == g.ne &&
                  c.at("faces") == g.nb && c.at("cells") == long(g.cells.size()) && sp.at("dx_m").at("bits") == bits(g.dx) &&
                  sp.at("dy_m").at("bits") == bits(g.dy) && sp.at("dz_m").at("bits") == bits(g.dz) && sp.at("zmin_m").at("bits") == bits(g.zmin);
  if (!ok) throw std::runtime_error("Hall mesh differs from the mesh file (counts or spacing bits)");
  mesh_record["checked_counts_and_spacing_bits"] = true;
}
Inputs load_inputs(const std::string& fixture) {
  Inputs in;
  in.settings = settings_for(fixture);
  const auto spec = selected_mesh();
  const axial::Mesh g(spec.transverse, spec.axial, spec.width_m, spec.air);
  in.setup = v7::make_source_setup(g, in.settings);
  in.model = v7::bind::make_hall_model(opt.hall_dt_max_s, spec, paths.trace, solver(),
                                        hc::CarrierLimits{in.settings.maximum_carriers, in.settings.hall_memory_budget},
                                        hb::declaration(in.setup.slices));
  check_mesh(in.model->mesh());
  std::ifstream f(paths.binding);
  J b;
  f >> b;
  if (!f) throw std::runtime_error("binding read");
  in.binding = {b.at("dataset_sha256"), b.at("kernel_sha256"), b.at("kernel_binding_file_sha256")};
  in.data = std::make_shared<const w::Data>(w::Data::read(paths.data));
  in.bridge = std::make_shared<const w::ParentBridge>(in.data, w::parent_model, in.binding.dataset_sha256, in.binding.kernel_sha256, v7::water_mass_kg);
  in.bands = std::make_shared<const cv::HighEnergyBands>(cv::HighEnergyBands::read(paths.outer3, outer3_sha, *in.data, in.binding.dataset_sha256,
                                                                                  paths.kernel_binding, in.binding.kernel_sha256,
                                                                                  in.binding.kernel_binding_file_sha256, cv::ProvisionalOptions{}));
  in.scope = std::make_shared<const cv::CoverageScope>(in.data, in.binding.dataset_sha256, *in.bridge, in.bands, cv::CoveragePolicy{});
  return in;
}
m::PicVector electron_u(double eV) {
  const double x = eV * pc::qe / (w::me * pc::c * pc::c);
  return {pc::c * std::sqrt(x * (x + 2)), 0, 0};
}
// Qualification population: co-located electron/H2O+ pairs (every node starts
// neutral) with electrons of 5 eV .. 20 keV around four dense water packets; ten
// pairs whose electrons hit the lower foil inner face within the first step.
std::vector<m::FullTimeAuthoredPacket> qual_packets() {
  std::vector<m::FullTimeAuthoredPacket> out;
  std::mt19937_64 rng(20261031);
  std::uniform_real_distribution<double> u(0, 1);
  const m::PicVector water[4] = {{.0040, .0040, -.004}, {.0050, .0045, .001}, {.0045, .0050, .004}, {.0042, .0046, -.0005}};
  for (int k = 0; k < 4; ++k) out.push_back({"w" + std::to_string(k), "H2O", "qual/authored-water", {1e15, 0}, water[k], {0, 0, 0}});
  for (int k = 0; k < 160; ++k) {
    const auto& c = water[k % 4];
    const m::PicVector x{c[0] + 2e-4 * (u(rng) - .5), c[1] + 2e-4 * (u(rng) - .5), c[2] + 2e-4 * (u(rng) - .5)};
    const double eV = std::exp(std::log(5.) + (std::log(2e4) - std::log(5.)) * u(rng));
    const double ct = 2 * u(rng) - 1, ph = 6.283185307179586 * u(rng), st = std::sqrt(std::max(0., 1 - ct * ct));
    const double s = electron_u(eV)[0];
    out.push_back({"e" + std::to_string(k), "electron", "qual/authored-electron", {1, 0}, x, {s * st * std::cos(ph), s * st * std::sin(ph), s * ct}});
    out.push_back({"i" + std::to_string(k), "H2O+", "qual/authored-ion", {1, 0}, x, {372 * (2 * u(rng) - 1), 372 * (2 * u(rng) - 1), 372 * (2 * u(rng) - 1)}});
  }
  for (int k = 0; k < 10; ++k) {
    const m::PicVector x{.003 + 3e-4 * k, .0020 + 2e-7, -.006 + 1.2e-3 * k};
    out.push_back({"se" + std::to_string(k), "electron", "qual/authored-sink-electron", {1, 0}, x, {0, -1e8, 0}});
    out.push_back({"si" + std::to_string(k), "H2O+", "qual/authored-sink-ion", {1, 0}, x, {0, 0, 0}});
  }
  return out;
}
v7::GenesisRequestV7 request(const Inputs& in, std::vector<m::FullTimeAuthoredPacket> packets) {
  v7::GenesisRequestV7 rq;
  rq.model = in.model;
  rq.options = v7::hall_geometry(in.model->mesh());
  rq.seed = in.model->dc_seed();
  rq.states = v7::pilot_states();
  rq.channels = v7::pilot_channels(in.scope->arena_binding(in.binding.kernel_binding_file_sha256));
  rq.packets = std::move(packets);
  rq.grid_counts = {45, 45, 5};
  if (opt.collision_grid == "mesh") rq.grid_counts = {std::uint32_t(in.model->mesh().nx), std::uint32_t(in.model->mesh().ny), 5};
  else if (opt.collision_grid != "fixed200") throw std::runtime_error("--collision-grid mesh|fixed200");
  rq.binding = in.binding;
  rq.limits = {};
  if (opt.fixture == "load") {   // collision subflow sized to the load (bound into the chemistry configuration of these runs)
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
  return rq;
}
// Dense fixture: qual_packets() plus blob_pairs co-located electron/H2O+ pairs of weight
// blob_weight, uniform in the box native x 2.6-3.4 mm, y 5.1-5.9 mm, z -3..3 mm (collision
// cells without authored water), electrons Maxwellian at blob_eV, ions at 0.025 eV.
std::vector<m::FullTimeAuthoredPacket> dense_packets() {
  auto out = qual_packets();
  std::mt19937_64 rng(20261101);
  std::uniform_real_distribution<double> u(0, 1);
  std::normal_distribution<double> nd(0, 1);
  const double se = std::sqrt(opt.blob_eV * pc::qe / w::me), si = std::sqrt(.025 * pc::qe / v7::water_mass_kg);
  for (long k = 0; k < opt.blob_pairs; ++k) {
    const m::PicVector x{.0026 + 8e-4 * u(rng), .0051 + 8e-4 * u(rng), -.003 + .006 * u(rng)};
    const m::PicVector ve{se * nd(rng), se * nd(rng), se * nd(rng)}, vi{si * nd(rng), si * nd(rng), si * nd(rng)};
    out.push_back({"be" + std::to_string(k), "electron", "dense/blob-electron", {opt.blob_weight, 0}, x, ve});
    out.push_back({"bi" + std::to_string(k), "H2O+", "dense/blob-ion", {opt.blob_weight, 0}, x, vi});
  }
  return out;
}
// Load fixture (host profiles): opt.load_carriers authored carriers in the gap (native x, y 2.2-6.8 mm,
// z -12..12 mm): a share opt.load_neutral_fraction of H2O carriers (weight 1e6, 300 K), the rest co-located
// electron/H2O+ pairs (weight 1; electrons Maxwellian at 2 eV, ions 0.025 eV), fixed generator seed.
std::vector<m::FullTimeAuthoredPacket> load_packets() {
  std::vector<m::FullTimeAuthoredPacket> out;
  std::mt19937_64 rng(20261102);
  std::uniform_real_distribution<double> u(0, 1);
  std::normal_distribution<double> nd(0, 1);
  const long n = opt.load_carriers, neutrals = long(double(n) * opt.load_neutral_fraction), pairs = (n - neutrals) / 2;
  const double se = std::sqrt(2. * pc::qe / w::me), si = std::sqrt(.025 * pc::qe / v7::water_mass_kg);
  auto point = [&] { return m::PicVector{.0022 + .0046 * u(rng), .0022 + .0046 * u(rng), -.012 + .024 * u(rng)}; };
  out.reserve(std::size_t(2 * pairs + neutrals));
  for (long k = 0; k < pairs; ++k) {
    const auto x = point();
    out.push_back({"le" + std::to_string(k), "electron", "load/electron", {1, 0}, x, {se * nd(rng), se * nd(rng), se * nd(rng)}});
    out.push_back({"li" + std::to_string(k), "H2O+", "load/ion", {1, 0}, x, {si * nd(rng), si * nd(rng), si * nd(rng)}});
  }
  const double sn = std::sqrt(1.380649e-23 * 300. / v7::water_mass_kg);
  for (long k = 0; k < neutrals; ++k)
    out.push_back({"ln" + std::to_string(k), "H2O", "load/water", {1e6, 0}, point(), {sn * nd(rng), sn * nd(rng), sn * nd(rng)}});
  return out;
}
std::vector<m::FullTimeAuthoredPacket> fixture_packets(const std::string& f) {
  return f == "qual" ? qual_packets() : f == "dense" ? dense_packets() : f == "load" ? load_packets() : std::vector<m::FullTimeAuthoredPacket>{};
}

// GPU push compared with the sealed CPU arithmetic on every fixed-point iterate.
class ParityMotionProvider final : public hc::MotionProvider {
  std::shared_ptr<hc::MotionProvider> device_;
  hc::CpuMotionProvider cpu_;
  std::string name_;

 public:
  long compared = 0, identical = 0;
  explicit ParityMotionProvider(std::shared_ptr<hc::MotionProvider> d) : device_(std::move(d)), name_(device_->name() + "+cpu-parity") {}
  const std::string& name() const noexcept override { return name_; }
  void begin(const hc::MotionContext& c, std::span<const m::ArenaRateEntry> census) override {
    device_->set_newborns(newborns_, birth_scope_);
    cpu_.set_newborns(newborns_, birth_scope_);
    device_->begin(c, census);
    cpu_.begin(c, census);
  }
  h::MotionReceipt motion(const h::MotionBinding& b, const hc::Vec& E, const hc::Vec& B, hc::FreeCharge& after) override {
    auto g = device_->motion(b, E, B, after);
    hc::FreeCharge ca;
    auto c = cpu_.motion(b, E, B, ca);
    bool words = ca.net_C.size() == after.net_C.size() && ca.gross_C == after.gross_C;
    for (std::size_t k = 0; words && k < ca.net_C.size(); ++k)
      words = std::bit_cast<std::uint64_t>(ca.net_C[k].hi) == std::bit_cast<std::uint64_t>(after.net_C[k].hi) &&
              std::bit_cast<std::uint64_t>(ca.net_C[k].lo) == std::bit_cast<std::uint64_t>(after.net_C[k].lo);
    const bool same = v7::bind::iterate_record_sha256(g) == v7::bind::iterate_record_sha256(c) && words && g.current_A_m() == c.current_A_m() &&
                      g.initial_charge_C() == c.initial_charge_C() && g.endpoint_charge_C() == c.endpoint_charge_C() &&
                      g.sinks().size() == c.sinks().size() && g.births().size() == c.births().size();
    ++compared;
    if (same) ++identical;
    return g;
  }
  void end() noexcept override {
    device_->end();
    cpu_.end();
  }
};
void attach_motion(v7::OwnersV7& o) {
  if (opt.push == "gpu") o.motion = std::make_shared<hg::GpuMotionProvider>(backends.gpu(), o.hall->motion_context(), opt.gpu_capacity);
  else if (opt.push != "cpu") throw std::runtime_error("unknown motion provider");
}
v7::OwnersV7 genesis(const Inputs& in, const std::string& fixture) {
  auto o = v7::GenesisV7::admit(request(in, fixture_packets(fixture)));
  o.dt_cap_s = ctl.cap_s;
  attach_motion(o);
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
PublishResult publish_measured(v7::OwnersV7& o, v7::TrialV7& t) {
  probe_calls = 0;
  probe_alloc[0] = probe_alloc[1] = -1;
  TA::set_publish_probe(t, &publish_probe);
  const auto t0 = steady::now();
  PublishResult r;
  r.published = v7::try_publish_v7(o, t);
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
// Ledger totals of the accepted state from the step's one population traversal (v7.1): exact
// weights by species and nuclei as exact superaccumulator sums (the same exact values and rounded
// ledger values as the per-carrier ArenaCount credits of v7), energies as before.
using XS = v7::ExactWeightSum;
struct Totals {
  std::uint64_t n_e = 0, n_i = 0, n_n = 0;
  XS w_e, w_i, w_n, h_nuclei, o_nuclei;
  m::ArenaCount stock;
  long double kinetic_J = 0, reference_J = 0;
  hc::Words q_free, q_conductor, q_exterior;
  double q_gross = 0;
};
// Coordinator v7 form (--host-work v7): its own population traversal and per-carrier ArenaCount credits.
Totals totals_legacy(const v7::OwnersV7& o) {
  Totals t;
  m::ArenaCount we, wi, wn, hn, on;
  const auto idx = v7::state_index(o.inputs.states);
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
    for (std::size_t j = 0; j < o.slices[i].patches.size(); ++j) t.stock = v7::exact_add(t.stock, o.sources[i]->inventory(j));
  for (const auto& x : o.hall->charge().particles.net_C) t.q_free = t.q_free + x;
  for (double x : o.hall->charge().particles.gross_C) t.q_gross += x;
  for (const auto& x : o.hall->charge().conductor_C) t.q_conductor = t.q_conductor + x;
  for (const auto& x : o.hall->charge().exterior_C) t.q_exterior = t.q_exterior + x;
  return t;
}
Totals totals(const v7::OwnersV7& o, std::span<const m::ArenaRateEntry> population) {
  Totals t;
  const auto idx = v7::state_index(o.inputs.states);
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
    for (std::size_t j = 0; j < o.slices[i].patches.size(); ++j) t.stock = v7::exact_add(t.stock, o.sources[i]->inventory(j));
  for (const auto& x : o.hall->charge().particles.net_C) t.q_free = t.q_free + x;
  for (double x : o.hall->charge().particles.gross_C) t.q_gross += x;
  for (const auto& x : o.hall->charge().conductor_C) t.q_conductor = t.q_conductor + x;
  for (const auto& x : o.hall->charge().exterior_C) t.q_exterior = t.q_exterior + x;
  return t;
}
double jnum(const J& j, std::initializer_list<const char*> path) {
  const J* x = &j;
  for (const char* k : path) {
    if (!x->is_object() || !x->contains(k)) return std::numeric_limits<double>::quiet_NaN();
    x = &x->at(k);
  }
  if (x->is_number()) return x->get<double>();
  if (x->is_object() && x->contains("value") && x->at("value").is_number()) return x->at("value").get<double>();
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
// Deterministic columns (excluded: wall-clock timings and allocation counts).
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
void fill_row(Ledger& L, const v7::OwnersV7& o, const v7::TrialV7& t, const J& pre, const PublishResult& pub, const Totals& before,
              const Totals& after, RunState& st, long step, double prepare_s) {
  const auto& n = t.numbers();
  const auto& tm = t.times();
  L.begin();
  L.put("step", step);
  L.put("time_s", o.hall->time_s());
  L.put("dt_s", pre.at("dt").get<double>());
  L.put("current_A", pre.at("current_A").get<double>());
  L.put("terminal_V", pre.at("terminal_V").get<double>());
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
  // charge
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
    L.put(k, pre.at(k).get<double>());
  const auto& kl = o.hall->kinetic_ledgers();
  L.put("q_collected_cum_C", kl.sinks.collected_charge_C().rounded());
  L.put("q_escaped_cum_C", kl.sinks.escaped_charge_C().rounded());
  // water inventory (exact)
  st.sunk_h2o = v7::exact_add(st.sunk_h2o, n.sunk_h2o);
  st.sunk_ion = v7::exact_add(st.sunk_ion, n.sunk_h2o_ion);
  XS water = after.w_n;
  water.add(after.w_i); water.add(after.stock); water.add(st.sunk_h2o); water.add(st.sunk_ion);
  const bool water_ok = water == st.water_initial;
  st.water_exact = st.water_exact && water_ok;
  L.put("h2o_source_stock", v7::exact_value(after.stock));
  L.put("h2o_native_neutral", after.w_n.value());
  L.put("h2o_native_ion", after.w_i.value());
  L.put("h2o_sunk_neutral_cum", v7::exact_value(st.sunk_h2o));
  L.put("h2o_sunk_ion_cum", v7::exact_value(st.sunk_ion));
  L.put("h2o_total", water.value());
  L.put("h2o_total_exactly_conserved", water_ok);
  XS H = after.h_nuclei, O = after.o_nuclei;
  for (int k = 0; k < 2; ++k) { H.add(after.stock); H.add(st.sunk_h2o); H.add(st.sunk_ion); }
  O.add(after.stock); O.add(st.sunk_h2o); O.add(st.sunk_ion);
  XS H0 = st.genesis.h_nuclei, O0 = st.genesis.o_nuclei;
  for (int k = 0; k < 2; ++k) H0.add(st.genesis.stock);
  O0.add(st.genesis.stock);
  const bool nuclei_ok = H == H0 && O == O0;
  st.nuclei_exact = st.nuclei_exact && nuclei_ok;
  L.put("H_nuclei", H.value());
  L.put("O_nuclei", O.value());
  L.put("nuclei_exactly_conserved", nuclei_ok);
  // energy
  for (const char* k : {"E_field_delta_J", "E_plasma_work_J", "E_joule_J", "E_outer_loss_J", "E_outer_supply_J", "E_circuit_work_J",
                        "E_kinetic_change_J", "E_kinetic_store_change_J", "E_thermal_delta_J", "E_wall_heat_J", "E_external_loads_J",
                        "E_collected_kinetic_J", "E_short_bridge_kinetic_J", "E_escaped_kinetic_J", "E_birth_kinetic_J", "E_combined_balance_J",
                        "E_step_residual_J", "E_step_scale_J", "E_roundoff_bound_J"})
    L.put(k, pre.at(k).get<double>());
  L.put("E_source_neutral_kinetic_J", n.neutral_kinetic_J);
  L.put("E_source_neutral_internal_J", n.neutral_internal_J);
  L.put("E_source_adsorbed_reference_J", n.neutral_adsorbed_reference_J);
  L.put("E_source_loads_J", n.source_load_J);
  L.put("E_emission_nottingham_J", n.emission_nottingham_J);
  L.put("E_emission_reference_J", n.emission_reference_J);
  L.put("E_emission_matching_work_unrepresented_J", n.emission_matching_work_J);
  L.put("E_emission_kinetic_J", n.emission_kinetic_J);
  const long double dK = after.kinetic_J - before.kinetic_J, dR = after.reference_J - before.reference_J;
  const long double chem_dK = dK - (long double)pre.at("E_kinetic_store_change_J").get<double>();
  const long double defect = chem_dK + dR;
  st.reaction_defect_cum += defect;
  L.put("E_census_kinetic_J", after.kinetic_J);
  L.put("E_census_reference_J", after.reference_J);
  L.put("E_chemistry_kinetic_delta_J", chem_dK);
  L.put("E_reaction_reference_delta_J", dR);
  L.put("E_reaction_defect_J", defect);
  L.put("E_reaction_defect_cum_J", st.reaction_defect_cum);
  // coverage (parent chemistry) ledger of the step
  const auto& cov = t.chemistry_state().coverage();
  L.put("cov_unresolved_pair_evaluations_step", jnum(cov, {"last_step", "totals", "unresolved", "pair_evaluations"}));
  L.put("cov_unresolved_upper_events_step", jnum(cov, {"last_step", "totals", "unresolved", "upper_bound_events"}));
  L.put("cov_unresolved_upper_physical_step", jnum(cov, {"last_step", "totals", "unresolved", "upper_bound_physical_reactions"}));
  L.put("cov_applied_expected_events_step", jnum(cov, {"last_step", "totals", "applied", "expected_events"}));
  L.put("cov_applied_expected_physical_step", jnum(cov, {"last_step", "totals", "applied", "expected_physical_reactions"}));
  double cats[3][2] = {{0, 0}, {0, 0}, {0, 0}};
  if (cov.contains("last_step") && cov.at("last_step").at("totals").contains("categories"))
    for (const auto& c : cov.at("last_step").at("totals").at("categories")) {
      const auto name = c.at("category").get<std::string>();
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
  // arrival ledger
  L.put("arr_expected_electrons_step", n.expected_electrons);
  L.put("arr_max_log10_expected_electrons_patch", n.max_log10_expected_electrons);
  L.put("arr_unresolved_emission_patches_step", n.unresolved_emission_patches);
  L.put("arr_unresolved_emission_lower_electrons_step", n.unresolved_emission_lower_electrons);
  L.put("arr_unresolved_emission_upper_electrons_step", n.unresolved_emission_upper_electrons);
  L.put("arr_expected_h2o_molecules_step", n.expected_neutral_molecules);
  L.put("arr_released_h2o_molecules_step", v7::exact_value(n.source_debit));
  std::uint64_t ce = 0, cn = 0;
  long double ee = 0, en = 0, var = 0, se = 0, sn = 0;
  for (const auto& x : o.arrival->totals()) {
    const bool electron = x.stream_id < v7::neutral_stream_base;
    (electron ? ce : cn) += x.births;
    (electron ? ee : en) += v7::exact_value(x.expected_sum);
    (electron ? se : sn) += v7::exact_value(x.sampled_sum);
    if (electron) var += x.variance_sum;
  }
  L.put("arr_cum_electron_births", ce);
  L.put("arr_cum_h2o_births", cn);
  L.put("arr_cum_expected_electrons", ee);
  L.put("arr_cum_sampled_electrons", se);
  L.put("arr_cum_expected_h2o_molecules", en);
  L.put("arr_cum_released_h2o_molecules", sn);
  L.put("arr_electron_standardized_deviation", var > 0 ? double((se - ee) / ((long double)o.settings.w0_electron * std::sqrt(var))) : 0.);
  L.put("src_patch_T_max_K", n.max_patch_T_K);
  L.put("src_patch_T_min_K", n.min_patch_T_K);
  L.put("src_extraction_field_max_V_m", n.max_extraction_field_V_m);
  L.put("src_extraction_field_min_V_m", n.min_extraction_field_V_m);
  L.put("hall_iterations", int(tm.hall_iterations));
  L.put("solve_iterations", pre.at("solve_iterations").get<double>());
  L.put("solve_residual_over_target", pre.at("solve_residual_over_target").get<double>());
  L.put("plasma_dt", pre.at("plasma_dt").get<double>());
  L.put("max_gyro", pre.at("max_gyro").get<double>());
  L.put("max_displacement_cells", pre.at("max_displacement_cells").get<double>());
  L.put("census_rows_step", n.census);
  L.put("newborns_step", n.newborns);
  // timings (ms)
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
  // summary extremes
  st.worst_step_residual = std::max(st.worst_step_residual, std::abs(pre.at("E_step_residual_J").get<double>()));
  const double gb = pre.at("gauss_bound_C").get<double>(), db = pre.at("divB_bound_T_m").get<double>();
  if (gb > 0) st.worst_gauss_ratio = std::max(st.worst_gauss_ratio, pre.at("gauss_max_abs_C").get<double>() / gb);
  if (db > 0) st.worst_divB_ratio = std::max(st.worst_divB_ratio, pre.at("divB_max_T_m").get<double>() / db);
  st.worst_q_change = std::max(st.worst_q_change, std::abs(q_change));
  st.worst_reaction_defect = std::max(st.worst_reaction_defect, double(std::abs(defect)));
  st.electron_births += long(n.electron_births);
  st.neutral_births += long(n.neutral_births);
  st.sinks += long(n.sinks);
  st.events += long(n.parent_events);
  L.first = false;
}
// Proposal numbers extracted before publication (the proposal moves into Hall).
J pre_numbers(const v7::TrialV7& t) {
  const double dt = t.dt_s();
  const auto& p = t.proposal();
  const auto& fh = p.field_heat();
  const auto& r = fh.field_receipt();
  const auto& c = fh.candidate_constraints();
  const auto& b = p.motion().budget();
  const auto& e = p.energy();
  const auto sinks = hc::SinkLedger::of(p.sinks());
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
// v6.1 sub-phase timings (not part of the v6 ledger layout).
const std::vector<std::string> subphase_columns = {
    "step", "t_start_field_digest_ms", "t_source_field_build_ms", "t_emission_extracting_ms", "t_emission_retarding_ms",
    "t_emission_other_ms", "extracting_patches", "retarding_patches", "unresolved_patches", "t_cp_native_ms", "t_cp_chemistry_ms",
    "t_cp_observe_ms", "t_cp_predicate_ms", "t_cp_sources_ms", "t_cp_arrival_ms", "t_cp_hall_ms", "cp_source_preflights"};
std::vector<double> subphase_row(const v7::TrialV7& t, long step) {
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
std::string digest_line(const v7::OwnersV7& o, long step, const std::string& rng) {
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
  void write(const v7::OwnersV7& o, long step, std::span<const m::ArenaRateEntry> population) {
    const auto& g = o.hall->mesh();
    const int nx = g.nx, ny = g.ny;
    std::vector<double> ne(nx * ny), ni(nx * ny), nn(nx * ny), ke(nx * ny), Ex(nx * ny), Ey(nx * ny), phi((nx + 1) * (ny + 1));
    std::vector<double> we(nx * ny);
    // 3D electron cells (NGP) for the Debye length.
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
    const auto& E = o.hall->field().E;
    const int k0 = g.nz / 2;   // z = 0 node plane
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
std::string hexd(double x) { std::ostringstream s; s << std::hexfloat << x; return s.str(); }
std::string bits(double x) { std::ostringstream s; s << std::hex << std::bit_cast<std::uint64_t>(x); return s.str(); }
// Largest unresolved upper-bound parent-event rate (events per second of exposure) of
// the exposures of the last accepted step (coverage ledger of the step).
double chemistry_rate(const J& cov) {
  double rate = 0;
  if (!cov.contains("last_step") || !cov.at("last_step").contains("exposures")) return 0;
  for (const auto& x : cov.at("last_step").at("exposures")) {
    const double d = x.value("duration_s", 0.);
    const auto& lim = x.at("exposure_limit").at("events_upper");
    const double up = lim.is_number() ? lim.get<double>() : INFINITY;
    if (d > 0 && up > 0) rate = std::max(rate, up / d);
  }
  return rate;
}
J criterion_codes() {
  J c = J::object();
  for (int k = 0; k <= int(v7::Criterion::retry); ++k) c[v7::criterion_name(v7::Criterion(k))] = k;
  return c;
}
void control_row(Ledger& C, long step, double t_begin, double t_end, const v7::Decision& d, const v7::Decision& first, const v7::ControlInputs& in,
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
  C.put("measured_plasma_dt", pre.at("plasma_dt").get<double>());
  C.put("measured_gyro", pre.at("max_gyro").get<double>());
  C.put("measured_displacement_cells", pre.at("max_displacement_cells").get<double>());
  C.put("t_control_ms", 1e3 * control_s);
}
// Accepted end-of-step E and B at fixed mesh entries near probe points (no
// interpolation) and the total field energies (sums in index order).
struct ProbeWriter {
  std::ofstream bin;
  std::vector<std::array<int, 6>> idx;
  J layout;
  void open(const fs::path& dir, const v7::OwnersV7& o) {
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
  void write(const v7::OwnersV7& o, long step) {
    const auto& g = o.hall->mesh();
    const auto& E = o.hall->field().E;
    const auto& B = o.hall->field().B;
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

// v7.1: wall time of the coordinator's own per-step work (outside the owners' prepare/publish).
struct CoordTimes {
  double control_s = 0, totals_s = 0, ledger_row_s = 0, ledger_io_s = 0, digest_line_s = 0, controller_state_s = 0, control_row_s = 0,
         probes_s = 0, receipt_s = 0, frames_s = 0, observe_prepare_s = 0, observe_preflight_s = 0, step_s = 0;
  // device time of the step (deltas of the GPU backends' counters; 0 on CPU backends)
  double gpu_solve_s = 0, gpu_push_upload_s = 0, gpu_push_admit_s = 0, gpu_push_device_s = 0, gpu_push_assemble_s = 0, gpu_push_total_s = 0;
};
struct DeviceCounters {
  double solve = 0, upload = 0, admit = 0, device = 0, assemble = 0, total = 0;
};
DeviceCounters device_counters(const Inputs& in, const std::shared_ptr<hc::MotionProvider>& motion) {
  DeviceCounters c;
  if (opt.field == "gpu") c.solve = hg::gpu_field_solver_stats(in.model->solver()).total_s;
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

// ---------------------------------------------------------------- run mode
int run_mode(const fs::path& out) {
  std::signal(SIGTERM, on_signal);
  std::signal(SIGINT, on_signal);
  const auto in = load_inputs(opt.fixture);
  J summary = {{"schema", "coordinator-v6-run-summary-v1"}, {"fixture", opt.fixture}, {"field_backend", opt.field}, {"push", opt.push},
               {"dt_s", opt.dt_s}, {"settings", in.settings.json_text()}, {"hall_configuration_sha256", in.model->configuration_sha256()}};
  if (backends.device) summary["device"] = J::parse(backends.device->device_json());
  const auto tg = steady::now();
  auto o = genesis(in, opt.fixture);
  summary["genesis"] = o.genesis;
  summary["genesis_s"] = since(tg);
  summary["motion_provider"] = o.motion->name();
  write_text(out / "genesis.json", o.genesis.dump(1) + "\n");
  RunState st;
  // The one population traversal per accepted state (v7.1): ledger totals, the next decision's
  // inputs, frames and the species verification of observe() all use it.
  const bool legacy = opt.host_work == "v7";
  o.legacy_host_work = legacy;
  std::vector<m::ArenaRateEntry> population = v7::accepted_population(o);
  st.genesis = legacy ? totals_legacy(o) : totals(o, population);
  st.water_initial = st.genesis.w_n;
  st.water_initial.add(st.genesis.w_i);
  st.water_initial.add(st.genesis.stock);
  st.q_initial = st.genesis.q_free + st.genesis.q_conductor + st.genesis.q_exterior;
  summary["initial"] = {{"source_stock_molecules", v7::exact_value(st.genesis.stock)}, {"native_h2o", st.genesis.w_n.value()},
                        {"native_ion", st.genesis.w_i.value()}, {"native_electron", st.genesis.w_e.value()},
                        {"total_charge_C", st.q_initial.rounded()}};
  std::ofstream ledger_bin(out / "ledger.bin", std::ios::binary), digests(out / "step-digests.txt"), det_bin(out / "deterministic.bin", std::ios::binary);
  std::ofstream sub_bin(out / "subphase-timings.bin", std::ios::binary);
  std::ofstream coord_bin(out / "coordinator-timings.bin", std::ios::binary);
  write_text(out / "coordinator-timings-schema.json", J{{"schema", "coordinator-v7.1-coordinator-timings-v1"}, {"byte_order", "little"},
                                                        {"value", "float64"}, {"columns", coordinator_columns}}.dump(1) + "\n");
  CoordTimes ct;
  DeviceCounters dc_prev{};
  write_text(out / "subphase-schema.json", J{{"schema", "coordinator-v6.1-subphase-timings-v1"}, {"byte_order", "little"}, {"value", "float64"},
                                             {"columns", subphase_columns}}.dump(1) + "\n");
  std::vector<std::vector<double>> sub_rows;
  FrameWriter frames;
  if (opt.frame_every > 0) {
    frames.open(out, o.hall->mesh());
    write_text(out / "frames-layout.json", frames.layout.dump(1) + "\n");
    frames.write(o, 0, population);
  }
  Ledger L;
  Totals before = st.genesis;
  const auto device_motion = o.motion;
  long parity_compared = 0, parity_identical = 0, steps = 0, max_swap_alloc = 0;
  std::string stop = "steps";
  const auto wall = steady::now();
  std::vector<double> step_ms;
  const double end_time = opt.window_s > 0 ? o.hall->time_s() + opt.window_s : opt.end_time_s;
  summary["end_time_s"] = end_time;
  // v7: controller state, control ledger, probes, sealed step receipts.
  v7::ControllerState cst;
  cst.record_state(o.hall->time_s(), o.hall->field().E, o.hall->field().B);
  Ledger C;
  std::ofstream control_bin(out / "control-ledger.bin", std::ios::binary), receipts(out / "step-receipts.jsonl");
  ProbeWriter probes;
  probes.open(out, o);
  probes.write(o, 0);
  std::string receipt_head = nw::sha("coordinator-v7-step-receipts-genesis/" + ctl.sha256() + "/" + o.arena->state_prefix() + "/" +
                                     o.hall->field_sha256() + "/" + o.accepted_source_seal);
  if (ctl.sync_every_s > 0 || ctl.fixed_grid_retry) ctl.sync_origin_s = o.hall->time_s();
  summary["controller"] = ctl.text();
  summary["hall_binding"] = v7::bind::hall_binding;
  summary["binding_pair"] = v7::bind::pair_name;
  summary["mesh"] = mesh_record;
  summary["collision_grid"] = opt.collision_grid;
  summary["host_work"] = opt.host_work;
  summary["native_binding"] = v7::bind::native_binding;
  summary["hall_configurable_cap"] = v7::bind::hall_configurable_cap;
  summary["field_owner_dt_max_s"] = v7::field_owner_dt_max(*in.model);
  summary["coupled_owner_dt_max_s"] = v7::owner_dt_max(*in.model);
  summary["controller_sha256"] = ctl.sha256();
  summary["receipt_chain_genesis"] = receipt_head;
  std::map<std::string, long> criteria_histogram;
  long refused_attempts_total = 0;
  double dt_min_used = INFINITY, dt_max_used = 0;
  dc_prev = device_counters(in, device_motion);
  for (long step = 1;; ++step) {
    if (end_time > 0) {
      // a sum of steps rounds: within 1e-6 of the shortest configured step counts as the end
      if (o.hall->time_s() >= end_time - 1e-6 * std::min(ctl.cap_s, ctl.fixed_dt_s)) { stop = "end_time"; break; }
    } else if (step > opt.steps) { stop = "steps"; break; }
    if (stop_signal) { stop = "signal"; break; }
    if (!opt.stop_file.empty() && step % 100 == 0 && fs::exists(opt.stop_file)) { stop = "stop_file"; break; }
    std::shared_ptr<ParityMotionProvider> parity;
    if (opt.push == "gpu" && opt.parity_every > 0 && (step == 1 || step % opt.parity_every == 0)) {
      parity = std::make_shared<ParityMotionProvider>(device_motion);
      o.motion = parity;
    } else o.motion = device_motion;
    // Decision from the accepted state, then attempts within the step transaction.
    const auto tc = steady::now();
    v7::CellBounds cb;
    if (legacy) cb.legacy_27_loop = true;
    const auto cin = legacy ? v7::control_inputs(o, cst, ctl, cb) : v7::control_inputs_from(o, population, cst, ctl, cb);
    const auto dec = v7::decide(ctl, cin, o.hall->model().drive(), step == 1);
    const double control_s = since(tc);
    if (dec.refused) {
      stop = std::string("controller found no admissible dt >= dt_min at step ") + std::to_string(step) + " (" + v7::criterion_name(dec.criterion) + ")";
      break;
    }
    auto attempt_dec = dec;
    J attempts = J::array();
    std::string decision_sha;
    std::optional<v7::TrialV7> trial;
    double prepare_s = 0;
    std::string refused;
    for (int attempt = 0;; ++attempt) {
      const J record = v7::decision_json(ctl, cin, attempt_dec);
      decision_sha = nw::sha(record.dump());
      const auto t0 = steady::now();
      try {
        if (opt.manufactured_refusal_step == step && attempt == 0)
          throw std::runtime_error("manufactured step-size refusal (test hook): motion full-census plasma guard");
        trial.emplace(v7::prepare_v7(o, attempt_dec.dt_s, 0, &decision_sha));
        prepare_s = since(t0);
        attempts.push_back({{"attempt", attempt}, {"dt_s", hexd(attempt_dec.dt_s)}, {"rung", attempt_dec.rung}, {"decision_sha256", decision_sha},
                            {"outcome", "prepared"}});
        break;
      } catch (const std::exception& e) {
        attempts.push_back({{"attempt", attempt}, {"dt_s", hexd(attempt_dec.dt_s)}, {"rung", attempt_dec.rung}, {"decision_sha256", decision_sha},
                            {"outcome", "refused"}, {"refusal", e.what()}, {"step_size_guard", v7::step_size_refusal(e.what())}});
        ++refused_attempts_total;
        const int rung = attempt_dec.rung + ctl.retry_rungs;
        const double next = ctl.adaptive ? v7::ladder_dt(ctl, rung) : .5 * attempt_dec.dt_s;
        if (!(ctl.adaptive || ctl.fixed_grid_retry) || !v7::step_size_refusal(e.what()) || attempt >= ctl.max_retries || next < ctl.dt_min_s) {
          refused = std::string("prepare refused at step ") + std::to_string(step) + ": " + e.what();
          break;
        }
        attempt_dec.rung = ctl.adaptive ? rung : attempt_dec.rung;
        attempt_dec.dt_s = next;
        attempt_dec.landing = false;
        attempt_dec.criterion = v7::Criterion::retry;
        attempt_dec.value = attempt + 1;
        attempt_dec.limit = ctl.max_retries;
      }
    }
    if (!refused.empty()) {
      stop = refused;
      receipts << J{{"schema", "coordinator-v7-step-refusal-v1"}, {"step", step}, {"previous_sha256", receipt_head}, {"attempts", attempts}}.dump() << '\n';
      break;
    }
    const double t_begin = o.hall->time_s();
    const J pre = pre_numbers(*trial);
    const auto pub = publish_measured(o, *trial);
    if (!pub.published) {
      stop = "publication refused at step " + std::to_string(step) + " (refusal point " + std::to_string(pub.refusal) + ")";
      receipts << J{{"schema", "coordinator-v7-step-refusal-v1"}, {"step", step}, {"previous_sha256", receipt_head}, {"attempts", attempts},
                    {"publication_refusal_point", pub.refusal}}.dump() << '\n';
      break;
    }
    max_swap_alloc = std::max(max_swap_alloc, long(pub.swap_allocations));
    if (std::getenv("COORDINATOR_V6_BANK_SCAN")) {
      const auto& g = o.hall->mesh();
      const auto& q = o.hall->charge();
      for (int nn = 0; nn < g.nn; ++nn) {
        const bool cond = g.node_body[nn] >= 0 || g.short_node[nn];
        const double c = q.conductor_C[nn].rounded();
        const double res = q.residue_C.empty() ? 0. : q.residue_C[nn];
        if (!cond && (c != 0 || res != 0)) {
          const int i = nn % (g.nx + 1), j = (nn / (g.nx + 1)) % (g.ny + 1), k = nn / ((g.nx + 1) * (g.ny + 1));
          std::cout << "bank-scan step " << step << " node " << nn << " ijk " << i << "," << j << "," << k << " conductor " << c << " residue " << res << std::endl;
        }
      }
    }
    auto tq = steady::now();
    if (!legacy) population = v7::accepted_population(o);
    else if (opt.frame_every > 0 && step % opt.frame_every == 0) population = o.arena->rate_population(o.inputs.maximum_carriers);
    const auto after = legacy ? totals_legacy(o) : totals(o, population);
    ct.totals_s = since(tq); tq = steady::now();
    fill_row(L, o, *trial, pre, pub, before, after, st, step, prepare_s);
    ct.ledger_row_s = since(tq); tq = steady::now();
    if (step == 1) write_text(out / "ledger-schema.json", J{{"schema", "coordinator-v6-step-ledger-v1"}, {"byte_order", "little"},
                                                            {"value", "float64"}, {"columns", L.names}}
                                                               .dump(1) + "\n");
    ledger_bin.write(reinterpret_cast<const char*>(L.row.data()), std::streamsize(L.row.size() * sizeof(double)));
    for (std::size_t k = 0; k < L.row.size(); ++k)
      if (deterministic_column(L.names[k])) det_bin.write(reinterpret_cast<const char*>(&L.row[k]), sizeof(double));
    ct.ledger_io_s = since(tq); tq = steady::now();
    const std::string rng = o.chemistry->accepted_rng_sha256();   // one RNG digest per step (digest line and receipt)
    digests << digest_line(o, step, rng);
    ct.digest_line_s = since(tq); tq = steady::now();
    {   // v7 records: controller history, control ledger, probes, sealed step receipt
      const int retries = int(attempts.size()) - 1;
      cst.dt_prev_s = attempt_dec.dt_s;
      cst.rung_prev = attempt_dec.rung;
      cst.hold = retries > 0 ? ctl.hold_after_retry : std::max(0, cst.hold - 1);
      cst.chem_rate = chemistry_rate(trial->chemistry_state().coverage());
      cst.births_rate = trial->numbers().expected_electrons / o.settings.w0_electron / attempt_dec.dt_s;
      cst.record_state(o.hall->time_s(), o.hall->field().E, o.hall->field().B);
      ct.controller_state_s = since(tq); tq = steady::now();
      control_row(C, step, t_begin, o.hall->time_s(), attempt_dec, dec, cin, int(attempts.size()), pre, control_s);
      if (step == 1) write_text(out / "control-schema.json", J{{"schema", "coordinator-v7-control-ledger-v1"}, {"byte_order", "little"},
                                                               {"value", "float64"}, {"columns", C.names},
                                                               {"criterion_codes", criterion_codes()}}.dump(1) + "\n");
      control_bin.write(reinterpret_cast<const char*>(C.row.data()), std::streamsize(C.row.size() * sizeof(double)));
      C.first = false;
      ct.control_row_s = since(tq); tq = steady::now();
      probes.write(o, step);
      ct.probes_s = since(tq); tq = steady::now();
      J receipt = {{"schema", "coordinator-v7-step-receipt-v1"}, {"step", step}, {"previous_sha256", receipt_head},
                   {"controller_sha256", ctl.sha256()}, {"decision", v7::decision_json(ctl, cin, attempt_dec)}, {"decision_sha256", decision_sha},
                   {"attempts", attempts}, {"common_clock_sha256", trial->common_clock_sha256()}, {"dt_s", hexd(attempt_dec.dt_s)},
                   {"dt_bits", bits(attempt_dec.dt_s)}, {"criterion", v7::criterion_name(attempt_dec.criterion)},
                   {"criterion_value", attempt_dec.value}, {"criterion_limit", std::isfinite(attempt_dec.limit) ? J(attempt_dec.limit) : J("inf")},
                   {"after", {{"time_s", hexd(o.hall->time_s())}, {"native_prefix", o.arena->state_prefix()}, {"field", o.hall->field_sha256()},
                              {"heat", o.hall->heat_sha256()}, {"kinetic_ledgers", o.hall->kinetic_ledgers_sha256()},
                              {"source_seal", o.accepted_source_seal}, {"rng", legacy ? o.chemistry->accepted_rng_sha256() : rng},
                              {"arrival", o.accepted_arrival_state}}}};
      receipt_head = nw::sha(receipt_head + receipt.dump());
      receipt["receipt_sha256"] = receipt_head;
      receipts << receipt.dump() << '\n';
      ct.receipt_s = since(tq); tq = steady::now();
      ++criteria_histogram[v7::criterion_name(attempt_dec.criterion)];
      dt_min_used = std::min(dt_min_used, attempt_dec.dt_s);
      dt_max_used = std::max(dt_max_used, attempt_dec.dt_s);
    }
    {
      const auto row = subphase_row(*trial, step);
      sub_bin.write(reinterpret_cast<const char*>(row.data()), std::streamsize(row.size() * sizeof(double)));
      sub_rows.push_back(row);
    }
    if (parity) { parity_compared += parity->compared; parity_identical += parity->identical; }
    before = after;
    ++steps;
    step_ms.push_back(1e3 * (prepare_s + pub.total_s));
    tq = steady::now();
    if (opt.frame_every > 0 && step % opt.frame_every == 0) frames.write(o, step, population);
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
    {
      const auto row = coordinator_row(ct, step, cin.carriers);
      coord_bin.write(reinterpret_cast<const char*>(row.data()), std::streamsize(row.size() * sizeof(double)));
      ct = {};
    }
    if (opt.progress_every > 0 && step % opt.progress_every == 0) {
      ledger_bin.flush(); digests.flush(); det_bin.flush(); sub_bin.flush(); control_bin.flush(); receipts.flush(); probes.flush(); coord_bin.flush();
      std::cout << "step " << step << " t=" << o.hall->time_s() << " ms/step(last)=" << step_ms.back() << " e=" << after.n_e << " i=" << after.n_i
                << " n=" << after.n_n << " births(e,n)=" << st.electron_births << "," << st.neutral_births << " sinks=" << st.sinks
                << " events=" << st.events << " rss=" << rss_kib() << std::endl;
    }
    if (opt.checkpoint_every > 0 && step % opt.checkpoint_every == 0) {
      const auto dir = out / ("diagnostic-checkpoint-" + std::to_string(step));
      fs::create_directories(dir);
      try {
        o.arena->save((dir / "native.json").string());
        o.hall->save(dir / "hall");
        for (std::size_t i = 0; i < o.sources.size(); ++i) o.sources[i]->save((dir / ("source-" + std::to_string(i) + ".bin")).string());
        write_text(dir / "arrival.bin", o.arrival->save());
        write_text(dir / "chemistry.json", o.chemistry->checkpoint().dump(1) + "\n");
        write_text(dir / "coordinator.json", J{{"step", step}, {"time_s", o.hall->time_s()}, {"accepted_source_seal", o.accepted_source_seal},
                                               {"note", "diagnostic only; restart is not implemented"}}.dump(1) + "\n");
      } catch (const std::exception& e) { std::cout << "diagnostic checkpoint at step " << step << " incomplete: " << e.what() << std::endl; }
    }
  }
  ledger_bin.flush(); digests.flush(); det_bin.flush(); sub_bin.flush(); control_bin.flush(); receipts.flush(); probes.flush();
  if (legacy) population = o.arena->rate_population(o.inputs.maximum_carriers);
  if (opt.frame_every > 0 && steps % opt.frame_every != 0) frames.write(o, steps, population);
  summary["receipt_chain_head"] = receipt_head;
  summary["criteria_histogram"] = criteria_histogram;
  summary["refused_attempts"] = refused_attempts_total;
  summary["dt_used_s"] = {{"min", std::isfinite(dt_min_used) ? J(dt_min_used) : J(nullptr)}, {"max", dt_max_used},
                          {"mean", steps > 0 ? (o.hall->time_s() - in.model->drive().start()) / double(steps) : 0.}};
  if (!sub_rows.empty()) {
    J med = J::object();
    for (std::size_t c = 1; c < subphase_columns.size(); ++c) {
      std::vector<double> v;
      for (const auto& r : sub_rows) v.push_back(r[c]);
      std::sort(v.begin(), v.end());
      med[subphase_columns[c]] = v[v.size() / 2];
    }
    summary["subphase_median"] = med;
  }
  summary["coordinator"] = "v7";
  summary["source_owner"] = "electrode_heat_v6_1 (pic-electrode-heat-source-v6.1-20261031)";
  summary["emission_statistics"] = emission_statistics_json();
  summary["steps_completed"] = steps;
  summary["stop"] = stop;
  summary["wall_s"] = since(wall);
  summary["final_time_s"] = o.hall->time_s();
  summary["simulated_s"] = o.hall->time_s() - in.model->drive().start();
  if (!step_ms.empty()) {
    auto s = step_ms;
    std::sort(s.begin(), s.end());
    double mean = 0;
    for (double x : step_ms) mean += x;
    summary["step_ms"] = {{"first", step_ms.front()}, {"median", s[s.size() / 2]}, {"mean", mean / double(s.size())}, {"p90", s[std::size_t(.9 * double(s.size() - 1))]},
                          {"max", s.back()}};
  }
  summary["max_allocations_in_swap_window"] = max_swap_alloc;
  summary["closure"] = {{"worst_abs_step_energy_residual_J", st.worst_step_residual}, {"worst_gauss_over_bound", st.worst_gauss_ratio},
                        {"worst_divB_over_bound", st.worst_divB_ratio}, {"worst_abs_total_charge_change_C", st.worst_q_change},
                        {"water_inventory_exact_every_step", st.water_exact}, {"H_O_nuclei_exact_every_step", st.nuclei_exact},
                        {"worst_abs_reaction_defect_J", st.worst_reaction_defect}, {"cumulative_reaction_defect_J", double(st.reaction_defect_cum)}};
  summary["totals"] = {{"electron_births", st.electron_births}, {"h2o_births", st.neutral_births}, {"sinks", st.sinks}, {"parent_events", st.events}};
  const auto fin = totals(o, population);
  summary["final"] = {{"carriers_electron", fin.n_e}, {"carriers_ion", fin.n_i}, {"carriers_h2o", fin.n_n},
                      {"source_stock_molecules", v7::exact_value(fin.stock)}, {"hall_kinetic_ledgers", o.hall->kinetic_ledgers().sinks.json()},
                      {"hall_birth_ledgers", o.hall->kinetic_ledgers().births.json()}, {"native_flow_totals_births", o.arena->full_time_flow_totals().births},
                      {"native_flow_totals_sunk", o.arena->full_time_flow_totals().sunk_carriers}, {"thermal_max_T_K", o.hall->thermal().max_temperature_K()},
                      {"chemistry_coverage", o.chemistry->accepted_coverage()}};
  {
    const auto hs = o.arena->full_time_history_status();
    summary["native_history"] = {{"anchors", hs.anchors}, {"anchor_bytes", hs.anchor_bytes}, {"window_rows", hs.window_rows},
                                 {"cumulative_rows", hs.cumulative_rows}, {"cumulative_bytes", hs.cumulative_bytes}};
  }
  if (opt.push == "gpu") summary["push_parity"] = {{"iterates_compared", parity_compared}, {"iterates_identical", parity_identical}};
  if (opt.field == "gpu") {
    const auto s = hg::gpu_field_solver_stats(in.model->solver());
    summary["gpu_field_solves"] = {{"solves", s.solves}, {"total_s", s.total_s}};
  }
  summary["forbidden_native_calls"] = forbidden_calls;
  if (const auto g = std::dynamic_pointer_cast<hg::GpuMotionProvider>(device_motion)) {
    const auto& x = g->stats();
    summary["gpu_push"] = {{"runs", x.runs}, {"upload_carriers_s", x.upload_carriers_s}, {"admit_s", x.admit_s}, {"device_run_s", x.device_run_s},
                           {"assemble_s", x.assemble_s}, {"total_s", x.total_s}};
  }
  summary["frames"] = frames.frames;
  summary["rss_kib_final"] = rss_kib();
  write_text(out / "run-summary.json", summary.dump(1) + "\n");
  std::cout << "run finished: steps=" << steps << " stop=" << stop << std::endl;
  const bool ok = (stop == "steps" || stop == "end_time" || stop == "stop_file" || stop == "signal") && parity_compared == parity_identical &&
                  st.water_exact && st.nuclei_exact && forbidden_calls == 0 && max_swap_alloc == 0;
  return ok ? 0 : 1;
}

// ---------------------------------------------------------------- vacuum pre-run cost
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
  hc::Owner hall(in.model, seed, q, v7::hall_geometry(g), {});
  struct Boundary { std::string field_sha256, heat_sha256; double time_s = 0; std::uint64_t accepted_step = 0; };
  J summary = {{"schema", "coordinator-v6-vacuum-prerun-cost-v1"}, {"field_backend", opt.field}, {"dt_s", opt.dt_s}};
  // dt guard of the field model (not loosened): a 20 fs vacuum step must refuse.
  std::string refusal;
  try { (void)hall.prepare_vacuum(2e-14); } catch (const std::exception& e) { refusal = e.what(); }
  summary["dt_2e-14_refusal"] = refusal;
  std::vector<double> ms, prep, solve;
  const auto wall = steady::now();
  for (long s = 0; s < opt.steps; ++s) {
    const auto t0 = steady::now();
    auto p = hall.prepare_vacuum(opt.dt_s);
    const double tp = since(t0);
    Boundary b{p.candidate_field_sha256(), p.candidate_heat_sha256(), p.end_s(), p.endpoint_revision()};
    if (!hall.preflight_vacuum_final(p, b)) throw std::runtime_error("vacuum preflight refused");
    hall.publish_vacuum_prevalidated(std::move(p));
    ms.push_back(1e3 * since(t0));
    prep.push_back(1e3 * tp);
  }
  summary["steps"] = opt.steps;
  summary["wall_s"] = since(wall);
  summary["final_time_s"] = hall.time_s();
  summary["max_temperature_K"] = hall.thermal().max_temperature_K();
  auto stats = [](std::vector<double> v) {
    std::sort(v.begin(), v.end());
    double mean = 0;
    for (double x : v) mean += x;
    return J{{"median", v[v.size() / 2]}, {"mean", mean / double(v.size())}, {"min", v.front()}, {"max", v.back()}};
  };
  summary["ms_per_step"] = stats(ms);
  summary["prepare_ms"] = stats(prep);
  if (opt.field == "gpu") {
    const auto s = hg::gpu_field_solver_stats(model.solver());
    summary["gpu_field_solves"] = {{"solves", s.solves}, {"total_s", s.total_s}};
  }
  write_text(out / "vacuum-summary.json", summary.dump(1) + "\n");
  std::cout << summary.dump(1) << std::endl;
  return refusal.empty() ? 1 : 0;
}

// ---------------------------------------------------------------- qualification checks
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
std::string owner_digest(const v7::OwnersV7& o, const fs::path& scratch) {
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
  Results R;
  const fs::path scratch = out / "scratch";
  fs::create_directories(scratch);
  const auto in = load_inputs("qual");
  const double dt = opt.dt_s;
  auto make = [&] { return genesis(in, "qual"); };
  // T1: one positive step, zero allocations in the swap window, second step preflights; reference digest.
  std::string reference_after;
  {
    auto o = make();
    const auto d0 = owner_digest(o, scratch);
    auto t = v7::prepare_v7(o, dt);
    const auto n = t.numbers();
    const auto pub = publish_measured(o, t);
    reference_after = owner_digest(o, scratch);
    R.record("T1", "positive coupled step published (births, motion, sinks, exposures)", pub.published,
             {{"electron_births", n.electron_births}, {"h2o_births", n.neutral_births}, {"sinks", n.sinks}, {"parent_events", n.parent_events},
              {"census", n.census}});
    R.record("T1", "zero allocations between the first and last swap", pub.published && pub.swap_allocations == 0, {{"allocations", pub.swap_allocations}});
    R.record("T1", "accepted owners changed by the step", d0 != reference_after);
    auto t2 = v7::prepare_v7(o, dt);
    R.record("T1", "a second step prepares and preflights", v7::can_publish_v7(o, t2));
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
    auto twin_seed = genesis(in2, "qual");
    const auto d0 = owner_digest(o, scratch);
    struct Case { std::string name; std::function<void(v7::TrialV7&, v7::TrialV7&)> mutate; };
    std::vector<Case> cases = {
        {"foreign source transition seal", [](v7::TrialV7& t, v7::TrialV7&) { TA::source_after(t) = nw::sha("foreign source seal"); }},
        {"Hall proposal of a twin owner set (dt/2)", [](v7::TrialV7& t, v7::TrialV7& x) { std::swap(TA::proposal(t), TA::proposal(x)); }},
        {"source trial of a twin source owner", [](v7::TrialV7& t, v7::TrialV7& x) { std::swap(TA::source_trials(t)[0], TA::source_trials(x)[0]); }},
        {"arrival staging of a twin arrival owner with another run seed", [](v7::TrialV7& t, v7::TrialV7& x) { std::swap(TA::staged(t), TA::staged(x)); }},
        {"native token of a twin arena", [](v7::TrialV7& t, v7::TrialV7& x) { std::swap(TA::native(t), TA::native(x)); }},
    };
    for (const auto& c : cases) {
      auto t = v7::prepare_v7(o, dt);
      auto& other = c.name.find("run seed") != std::string::npos ? twin_seed : twin;
      auto x = v7::prepare_v7(other, c.name.find("dt/2") != std::string::npos ? .5 * dt : dt);
      c.mutate(t, x);
      const bool refused = !v7::try_publish_v7(o, t);
      const int point = t.refusal_point();
      const bool unchanged = owner_digest(o, scratch) == d0;
      R.record("T2", c.name + ": refused, every owner unchanged", refused && unchanged, {{"refusal_point", point}});
    }
    {   // live source owner record differs from the accepted digest
      auto t = v7::prepare_v7(o, dt);
      const auto keep = o.accepted_source_states[0];
      o.accepted_source_states[0] = nw::sha("foreign source owner state");
      const bool refused = !v7::try_publish_v7(o, t);
      o.accepted_source_states[0] = keep;
      R.record("T2", "live source owner state differs from the accepted record: refused, owners unchanged", refused && owner_digest(o, scratch) == d0,
               {{"refusal_point", t.refusal_point()}});
    }
    auto t = v7::prepare_v7(o, dt);
    const bool pub = v7::try_publish_v7(o, t);
    R.record("T2", "after every refusal the unmutated step publishes the T1 reference state bit for bit", pub && owner_digest(o, scratch) == reference_after);
  }
  // T3: manufactured faults in preparation and in can_publish refuse with owners unchanged.
  {
    auto o = make();
    const auto d0 = owner_digest(o, scratch);
    int prepare_ok = 0, publish_ok = 0;
    for (int k = 1; k <= v7::prepare_fault_points_v7; ++k) {
      bool threw = false;
      try { (void)v7::prepare_v7(o, dt, k); } catch (const std::exception&) { threw = true; }
      prepare_ok += threw && owner_digest(o, scratch) == d0;
    }
    for (int k = 1; k <= 10; ++k) {
      auto t = v7::prepare_v7(o, dt);
      TA::set_can_publish_fault(t, k);
      const bool refused = !v7::try_publish_v7(o, t);
      publish_ok += refused && t.refusal_point() == 100 + k && owner_digest(o, scratch) == d0;
    }
    R.record("T3", "every preparation fault point refuses with every owner unchanged", prepare_ok == v7::prepare_fault_points_v7, {{"passed_points", prepare_ok}});
    R.record("T3", "every can_publish fault point refuses with every owner unchanged", publish_ok == 10, {{"passed_points", publish_ok}});
    auto t = v7::prepare_v7(o, dt);
    R.record("T3", "after the fault sweep a fresh step reproduces the T1 reference", v7::try_publish_v7(o, t) && owner_digest(o, scratch) == reference_after);
  }
  // T4: single use, sibling outer reservation, moved-from trials.
  {
    auto o = make();
    const auto d0 = owner_digest(o, scratch);
    auto a = v7::prepare_v7(o, dt);
    std::string sibling;
    try { (void)v7::prepare_v7(o, dt); } catch (const std::exception& e) { sibling = e.what(); }
    R.record("T4", "a second concurrent step preparation (sibling outer reservation) refuses, owners unchanged",
             !sibling.empty() && owner_digest(o, scratch) == d0, {{"refusal", sibling}});
    auto moved = std::move(a);
    const bool from_moved = v7::try_publish_v7(o, a);
    R.record("T4", "a moved-from trial refuses", !from_moved && owner_digest(o, scratch) == d0);
    const bool first = v7::try_publish_v7(o, moved);
    const auto d1 = owner_digest(o, scratch);
    const bool again = v7::try_publish_v7(o, moved);
    R.record("T4", "the trial publishes once; a second publication refuses with owners unchanged", first && !again && owner_digest(o, scratch) == d1);
    R.record("T4", "the published state is the T1 reference", d1 == reference_after);
  }
  // T6: genesis admission.
  {
    auto refuses = [&](const std::string& what, std::function<void(v7::GenesisRequestV7&)> edit, v7::GenesisV7::PacketHook ph = nullptr,
                       v7::GenesisV7::LedgerHook lh = nullptr) {
      auto rq = request(in, qual_packets());
      edit(rq);
      int reached = 0;
      std::string msg;
      try { (void)v7::GenesisTestAccessV7::admit(rq, ph, lh, &reached); } catch (const std::exception& e) { msg = e.what(); }
      R.record("T6", what, !msg.empty(), {{"refusal", msg.substr(0, 400)}, {"stage", reached}});
    };
    refuses("a lone electron (no compensating charge) refuses (absolute Gauss)", [](v7::GenesisRequestV7& rq) {
      rq.packets.push_back({"lone", "electron", "qual/lone", {1, 0}, {.0045, .0045, .002}, {0, 0, 0}});
    });
    refuses("genesis time other than the first C19 knot refuses", [](v7::GenesisRequestV7& rq) { rq.seed.time_s = std::nextafter(rq.seed.time_s, 1.); });
    refuses("a declared gross charge that differs from the census refuses", [](v7::GenesisRequestV7&) {}, nullptr,
            [](hc::ChargeLedger& q) { q.particles.gross_C[q.particles.gross_C.size() / 2] += 1e-30; });
    refuses("a missing foil-slice source owner refuses", [](v7::GenesisRequestV7& rq) { rq.sources.slices.pop_back(); });
    refuses("a source patch not declared on the Hall thermal owner refuses", [](v7::GenesisRequestV7& rq) {
      rq.sources.slices[0].patches[0].initial_particles *= 2;
      rq.sources.slices[0].patches[0].id += "-altered";
    });
    refuses("an authored packet in the reserved arrival namespace refuses", [](v7::GenesisRequestV7&) {},
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
      auto fld = std::make_shared<const s6::Field>(o.hall->mesh(), f, o.hall->field().E, o.hall->field().B);
      const s6::FinalFieldView view{fld->frame(), o.hall->mesh(), o.hall->field().E, o.hall->field().B};
      std::vector<std::vector<s6::NeutralLeaf>> leaves(o.slices[0].patches.size());
      leaves[0].push_back({1.5 * o.settings.w0_neutral, eh::Draw{.5, .5, .5, .5, .5, .5}});
      (void)src.prepare(e, fld, view, {}, {}, leaves);
    });
    R.record("spot", "source v6 refuses a neutral leaf whose weight is not w0", msg.find("exactly w0") != std::string::npos, {{"refusal", msg}});
    const auto plan = refusal([&] { pa::StepInput x; x.step = o.arrival->next_step() + 1; (void)o.arrival->plan(x); });
    R.record("spot", "arrival v1.1 refuses a plan for another than the next accepted step", !plan.empty(), {{"refusal", plan}});
    const auto vac = refusal([&] { (void)o.hall->prepare_vacuum(dt); });
    R.record("spot", "Hall refuses a field-only (vacuum) step while carriers are present", !vac.empty(), {{"refusal", vac}});
    const auto big = refusal([&] { (void)v7::prepare_v7(o, 2e-14); });
    R.record("spot", "the coupled step refuses dt above the unchanged 10 fs guard", !big.empty(), {{"refusal", big}});
    R.record("spot", "every spot-check refusal left every owner unchanged", owner_digest(o, scratch) == d0);
    auto t = v7::prepare_v7(o, dt);
    R.record("spot", "after the spot checks a fresh step reproduces the T1 reference", v7::try_publish_v7(o, t) && owner_digest(o, scratch) == reference_after);
  }
  // V7: variable dt, predicate v7, step-size refusals and retries, controller determinism.
  {
    auto o = make();
    const auto d0 = owner_digest(o, scratch);
    // Predicate v7 against the native v6 (V3) predicate on the same observations.
    int compared = 0, equal = 0;
    {
      auto t = v7::prepare_v7(o, dt);
      TA::compare_v6_predicate(t);
      const bool ok = v7::can_publish_v7(o, t);
      const auto [p6, p7] = TA::predicates(t);
      compared += p6 >= 0; equal += p6 >= 0 && p6 == p7;
      R.record("V7", "predicate v7 admits the 10 fs step and agrees with the V3/native-v6 predicate", ok && p6 == 1 && p7 == 1, {{"v6", p6}, {"v7", p7}});
    }
    {
      auto t = v7::prepare_v7(o, dt);
      TA::compare_v6_predicate(t);
      TA::source_after(t) = nw::sha("foreign source seal");
      const bool refused = !v7::try_publish_v7(o, t);
      const auto [p6, p7] = TA::predicates(t);
      compared += p6 >= 0; equal += p6 >= 0 && p6 == p7;
      R.record("V7", "foreign source seal: both predicates refuse, owners unchanged", refused && p6 == 0 && p7 == 0 && owner_digest(o, scratch) == d0,
               {{"refusal_point", t.refusal_point()}, {"v6", p6}, {"v7", p7}});
    }
    {
      auto t = v7::prepare_v7(o, dt);
      TA::compare_v6_predicate(t);
      auto& c = TA::clock(t);
      c.decided_dt_s = std::nextafter(c.decided_dt_s, 0.);
      const bool refused = !v7::try_publish_v7(o, t);
      const auto [p6, p7] = TA::predicates(t);
      R.record("V7", "a step whose dt differs from the controller decision by one ulp is refused by predicate v7, owners unchanged",
               refused && t.refusal_point() == 7 && p7 == 0 && p6 == 1 && owner_digest(o, scratch) == d0, {{"refusal_point", t.refusal_point()}, {"v6", p6}, {"v7", p7}});
    }
    {
      auto t = v7::prepare_v7(o, dt);
      auto& c = TA::clock(t);
      c.cap_s = 2e-14;
      const bool refused = !v7::try_publish_v7(o, t);
      R.record("V7", "a step clock with a cap above the field owner maximum is refused by predicate v7, owners unchanged",
               refused && t.refusal_point() == 7 && owner_digest(o, scratch) == d0, {{"refusal_point", t.refusal_point()}});
    }
    {
      auto t = v7::prepare_v7(o, dt);
      TA::clock(t).decision_sha256 = "not-a-digest";
      const bool refused = !v7::try_publish_v7(o, t);
      R.record("V7", "a step clock without a decision digest is refused by predicate v7, owners unchanged",
               refused && t.refusal_point() == 7 && owner_digest(o, scratch) == d0, {{"refusal_point", t.refusal_point()}});
    }
    auto refusal = [](auto&& f) { try { f(); } catch (const std::exception& e) { return std::string(e.what()); } return std::string(); };
    {
      o.dt_cap_s = 5e-15;
      const auto big = refusal([&] { (void)v7::prepare_v7(o, dt); });
      o.dt_cap_s = 2e-14;
      const auto over = refusal([&] { (void)v7::prepare_v7(o, 5e-15); });
      o.dt_cap_s = ctl.cap_s;
      R.record("V7", "prepare refuses dt above a configured cap of 5 fs and a configured cap above the field owner maximum, owners unchanged",
               !big.empty() && !over.empty() && owner_digest(o, scratch) == d0, {{"refusal_dt", big}, {"refusal_cap", over}});
    }
    // Variable steps publish: ladder rungs 3 and 8 (5 fs) and the predicates agree.
    {
      bool ok = true;
      for (int j : {3, 8}) {
        auto t = v7::prepare_v7(o, v7::ladder_dt(ctl, j));
        TA::compare_v6_predicate(t);
        ok = ok && v7::try_publish_v7(o, t);
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
    try { (void)v7::prepare_v7(a, dt, 108); } catch (const std::exception& e) { msg = e.what(); }
    const bool unchanged = owner_digest(a, scratch) == d0;
    const double retry = v7::ladder_dt(ctl, ctl.retry_rungs);
    auto ta = v7::prepare_v7(a, retry);
    const bool pa_ok = v7::try_publish_v7(a, ta);
    auto tb = v7::prepare_v7(b, retry);
    const bool pb_ok = v7::try_publish_v7(b, tb);
    R.record("V7", "a step-size refusal after the Hall proposal is classified as a step-size guard and leaves every owner unchanged",
             v7::step_size_refusal(msg) && unchanged, {{"refusal", msg}});
    R.record("V7", "the retried step (5 fs) equals a direct 5 fs step from the same state bit for bit", pa_ok && pb_ok && owner_digest(a, scratch) == owner_digest(b, scratch));
    R.record("V7", "a non-step-size fault is not retried", !v7::step_size_refusal("manufactured v6 preparation fault") &&
             !v7::step_size_refusal("absolute Gauss constraint disagrees with actual declared charge"));
  }
  // Controller inputs: partition/order independence, repeatability, Hall plasma operand.
  {
    const auto dense = load_inputs("dense");
    auto o = genesis(dense, "dense");
    for (int k = 0; k < 3; ++k) { auto t = v7::prepare_v7(o, dt); if (!v7::try_publish_v7(o, t)) break; }
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
    auto t = v7::prepare_v7(o, d1.dt_s);
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
    auto o = genesis(dense, "dense");
    opt.blob_weight = w;
    const auto d0 = owner_digest(o, scratch);
    std::string msg;
    try { (void)v7::prepare_v7(o, dt); } catch (const std::exception& e) { msg = e.what(); }
    const bool unchanged = owner_digest(o, scratch) == d0;
    v7::ControllerState cst;
    v7::CellBounds cb;
    auto adaptive = ctl;
    adaptive.adaptive = true;
    const auto in = v7::control_inputs(o, cst, adaptive, cb);
    const auto d = v7::decide(adaptive, in, o.hall->model().drive(), true);
    auto t = v7::prepare_v7(o, d.dt_s);
    const bool pub = v7::try_publish_v7(o, t);
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
}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 9)
      throw std::runtime_error("usage: driver_v7 run|checks|vacuum|profile TRACE H2O_DATA BINDING OUTER3 KERNEL_BINDING OUTDIR NATIVE_IDENTITY [options]");
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
      else if (a == "--steps") opt.steps = std::stol(value());
      else if (a == "--end-time") opt.end_time_s = std::stod(value());
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
      else if (a == "--window") opt.window_s = std::stod(value());
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
      else throw std::runtime_error("unknown option " + a);
    }
    ctl.fixed_dt_s = opt.dt_s;
    ctl.validate(std::min(opt.hall_dt_max_s, v7::bind::source_owner_dt_max_s));
    fs::create_directories(out);
    int rc = 2;
    if (mode == "run") rc = run_mode(out);
    else if (mode == "checks") rc = checks_mode(out);
    else if (mode == "vacuum") rc = vacuum_mode(out);
    else throw std::runtime_error("unknown mode");
    if (backends.device) write_text(out / "gpu-teardown.json", backends.device->close_json() + "\n");
    return rc;
  } catch (const std::exception& e) {
    std::cerr << "fatal: " << e.what() << '\n';
    return 2;
  }
}
