// Coordinator v7.1 fixtures and the GPU push parity provider, unchanged (included by driver_v8.cpp).
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
std::vector<m::FullTimeAuthoredPacket> dense_packets() {
  auto out = qual_packets();
  std::mt19937_64 rng(20261101);
  std::uniform_real_distribution<double> u(0, 1);
  std::normal_distribution<double> nd(0, 1);
  const double se = std::sqrt(opt.blob_eV * pc::qe / w::me), si = std::sqrt(.025 * pc::qe / v8::water_mass_kg);
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
  const double se = std::sqrt(2. * pc::qe / w::me), si = std::sqrt(.025 * pc::qe / v8::water_mass_kg);
  auto point = [&] { return m::PicVector{.0022 + .0046 * u(rng), .0022 + .0046 * u(rng), -.012 + .024 * u(rng)}; };
  out.reserve(std::size_t(2 * pairs + neutrals));
  for (long k = 0; k < pairs; ++k) {
    const auto x = point();
    out.push_back({"le" + std::to_string(k), "electron", "load/electron", {1, 0}, x, {se * nd(rng), se * nd(rng), se * nd(rng)}});
    out.push_back({"li" + std::to_string(k), "H2O+", "load/ion", {1, 0}, x, {si * nd(rng), si * nd(rng), si * nd(rng)}});
  }
  const double sn = std::sqrt(1.380649e-23 * 300. / v8::water_mass_kg);
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
    const bool same = v8::bind::iterate_record_sha256(g) == v8::bind::iterate_record_sha256(c) && words && g.current_A_m() == c.current_A_m() &&
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
