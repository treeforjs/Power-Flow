#pragma once
#include "beb_rbeb.hpp"
#include "native_parent_subflow.hpp"
#include <fstream>
#include <iterator>
#include <optional>

// Parent-ionization rates above the Song 2021 parent table (1 keV):
//  (1000, 10000] eV  dataset-backed: H2O:conditional_outer3:parent of the frozen
//                    outer3 dataset (h2o-gas-20261029-outer3-v5), the NIST
//                    four-orbital BEB formula with the orbital_weights
//                    [2a1, 1b2, 3a1, 1b1] = [0, 0.08, 1, 1] of model
//                    outer3_BEB_Tan_parent_v1 (Tan et al. 1978 high-energy
//                    single-vacancy branching as transcribed by Murakami et al.
//                    2012, Eq. 27). Evaluated with the dataset arithmetic.
//  (10000, max] eV   PROVISIONAL, model-derived: the same orbitals and the same
//                    per-orbital parent weights applied to relativistic BEB
//                    (Kim, Santos & Parente 2000) orbital cross sections.
// No other channel, total or yield is used as a parent rate.
namespace raikou::native_water::coverage_v2 {
namespace m = ::raikou::mfem_backend;
namespace pc = ::raikou::pic_collisions;
namespace w = ::raikou::pic_water;
using Json = m::ArenaJson;
using ::raikou::native_water::sha;

inline constexpr const char* outer3_channel = "H2O:conditional_outer3:parent";
inline constexpr const char* outer3_model = "outer3_BEB_Tan_parent_v1";
inline constexpr const char* provisional_model = "provisional_RBEB_KSP2000_x_Tan_parent_v1";
inline constexpr const char* sdcs_kernel = "parent_ground_BEB_SDCS_binary_encounter_v2";
inline constexpr const char* banded_model = "parent_banded_song2021_outer3BEB_provisionalRBEB_v2";
// Electron rest energy, CODATA 2022 (the value used by the dataset's
// collision_data.h speed() and equal to w::me c^2/e to 1e-11).
inline constexpr double electron_rest_energy_eV = 510998.95069;

struct ArenaBinding {
  std::string dataset_sha256, kernel_manifest_sha256, kernel_binding_file_sha256, model;
  Json json() const {
    return {{"dataset_sha256", dataset_sha256}, {"kernel_manifest_sha256", kernel_manifest_sha256},
            {"kernel_binding_file_sha256", kernel_binding_file_sha256}, {"model", model}};
  }
};
struct ProvisionalOptions {
  bool admit = true;
  // Finite support of the provisional band. Above it a pair is unresolved
  // without a bound (refused unless admit_unbounded).
  double maximum_eV = 1e7;
};

inline std::string read_bytes(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  pc::require(bool(in), "cannot open high-energy input");
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

class HighEnergyBands {
  std::string outer3_sha_, song_sha_, kernel_manifest_, kernel_binding_file_, song_kernel_, song_binding_file_;
  Json outer3_channel_json_, outer3_model_json_;
  std::vector<beb_v2::Orbital> orbitals_;
  std::vector<double> weights_;
  beb_v2::Constants k_;
  double beb_lo_ = 0, beb_hi_ = 0;
  ProvisionalOptions provisional_;
  Json description_;
  ArenaBinding arena_;

  static double num(const Json& j) {
    pc::require(j.is_number(), "high-energy metadata number");
    const double x = j.get<double>();pc::require(std::isfinite(x), "high-energy metadata finite");return x;
  }

public:
  // outer3_bytes: exact bytes of h2o-gas-20261029-outer3.json; song_data: the
  // bound v1 dataset; kernel_binding_bytes: exact bytes of kernel-binding-v2.json.
  HighEnergyBands(const std::string& outer3_bytes, const std::string& expected_outer3_sha256,
                  const w::Data& song_data, const std::string& song_dataset_sha256,
                  const std::string& kernel_binding_bytes, const std::string& song_kernel_sha256,
                  const std::string& song_kernel_binding_file_sha256, ProvisionalOptions provisional = {})
      : provisional_(provisional) {
    outer3_sha_ = sha(outer3_bytes);
    pc::require(outer3_sha_ == expected_outer3_sha256, "outer3 dataset bytes differ from the bound digest");
    pc::require(pc::digest(song_dataset_sha256), "song dataset digest");
    song_sha_ = song_dataset_sha256;
    const Json d = Json::parse(outer3_bytes);
    pc::require(d.at("schema") == "raikou-strict-gas-collision-data-v1" && d.at("id") == "h2o-gas-20261029-outer3-v5" &&
                    d.at("phase") == "gas" && d.at("target") == "H2O",
                "outer3 dataset identity");
    const Json* ch = nullptr;const Json* song_copy = nullptr;const Json* valence = nullptr;
    for (const auto& c : d.at("channels")) {
      if (c.at("id") == outer3_channel) { pc::require(!ch, "duplicate outer3 parent channel");ch = &c; }
      if (c.at("id") == w::parent_channel) song_copy = &c;
      if (c.at("id") == "H2O:nist_BEB_valence") valence = &c;
    }
    pc::require(ch && song_copy && valence, "outer3 dataset lacks the parent, Song or valence BEB channel");
    // The outer3 freeze carries the v1 Song parent table unchanged; the 1 keV
    // seam therefore compares the same table the bound bridge uses.
    pc::require(*song_copy == song_data.channel_metadata(w::parent_channel), "outer3 dataset Song parent table differs from the bound v1 table");
    pc::require(ch->at("model") == "nist_BEB" && ch->at("unit_SI") == "m2" && ch->at("outside_support") == "error" &&
                    ch->at("zero_below_eV").is_null() && ch->at("observable") == "conditional_exclusive_event_rate" &&
                    ch->at("exclusive_products") == Json{{"e_g", 2}, {"H2Op", 1}},
                "outer3 parent channel contract");
    beb_lo_ = num(ch->at("support_eV")[0]);beb_hi_ = num(ch->at("support_eV")[1]);
    pc::require(beb_lo_ == 1000 && beb_hi_ == 10000, "outer3 parent support must be 1000..10000 eV");
    const auto& orb = ch->at("orbitals");const auto& wts = ch->at("orbital_weights");
    pc::require(orb.is_array() && orb.size() == 4 && wts.is_array() && wts.size() == 4, "outer3 orbital dimensions");
    pc::require(orb == valence->at("orbitals") && ch->at("bohr_radius_m") == valence->at("bohr_radius_m") &&
                    ch->at("rydberg_eV") == valence->at("rydberg_eV"),
                "outer3 orbital constants differ from the dataset valence BEB");
    const std::vector<std::string> labels{"2a1", "1b2", "3a1", "1b1"};
    for (std::size_t i = 0; i < 4; ++i) {
      beb_v2::Orbital o;o.label = orb[i].at("label").get<std::string>();
      o.B_eV = num(orb[i].at("B_eV"));o.U_eV = num(orb[i].at("U_eV"));o.N = num(orb[i].at("N"));o.Q = num(orb[i].at("Q"));
      pc::require(o.label == labels[i] && o.B_eV > 0 && o.U_eV > 0 && o.N > 0 && o.Q == 1, "outer3 orbital record");
      orbitals_.push_back(o);
      const double x = num(wts[i]);pc::require(x >= 0 && x <= 1, "orbital weight");weights_.push_back(x);
    }
    pc::require(weights_ == std::vector<double>{0, .08, 1, 1}, "outer3 parent weights must be Tan/Murakami [0, 0.08, 1, 1]");
    k_.bohr_radius_m = num(ch->at("bohr_radius_m"));k_.rydberg_eV = num(ch->at("rydberg_eV"));
    k_.rest_energy_eV = electron_rest_energy_eV;
    pc::require(std::abs(static_cast<long double>(w::me) * pc::c * pc::c / pc::qe - electron_rest_energy_eV) <=
                    1e-9L * electron_rest_energy_eV,
                "electron rest energy differs from the bridge mass");
    for (const auto& model : d.at("conditional_event_models"))
      if (model.at("id") == outer3_model) outer3_model_json_ = model;
    pc::require(outer3_model_json_.is_object() && outer3_model_json_.at("channel") == outer3_channel &&
                    outer3_model_json_.at("data_support_eV") == Json::array({1000., 10000.}) &&
                    outer3_model_json_.at("products") == Json{{"e_g", 2}, {"H2Op", 1}} &&
                    outer3_model_json_.at("parent_debit") == 1 && outer3_model_json_.at("electron_net_gain") == 1 &&
                    num(outer3_model_json_.at("chemical_reference_gain_eV")) == w::threshold_eV,
                "outer3 parent model contract");
    outer3_channel_json_ = *ch;
    pc::require(provisional_.maximum_eV > beb_hi_ && std::isfinite(provisional_.maximum_eV), "provisional band support");
    // Kernel binding: canonical manifest over the v2 kernel sources plus the
    // bound Song-band kernel identity.
    const Json kb = Json::parse(kernel_binding_bytes);
    pc::require(kb.at("schema") == "parent-coverage-v2-kernel-binding" && kb.at("sources").is_object() &&
                    kb.at("song_band_kernel_manifest_sha256") == song_kernel_sha256 &&
                    kb.at("song_band_kernel_binding_file_sha256") == song_kernel_binding_file_sha256 &&
                    kb.at("high_energy_kernel") == sdcs_kernel,
                "kernel binding file contract");
    const Json manifest{{"schema", kb.at("schema")}, {"sources", kb.at("sources")},
                        {"song_band_kernel_manifest_sha256", song_kernel_sha256},
                        {"song_band_kernel_binding_file_sha256", song_kernel_binding_file_sha256},
                        {"high_energy_kernel", sdcs_kernel}};
    kernel_manifest_ = sha(manifest.dump());
    pc::require(kb.at("canonical_manifest_sha256") == kernel_manifest_, "kernel binding canonical manifest differs");
    kernel_binding_file_ = sha(kernel_binding_bytes);
    song_kernel_ = song_kernel_sha256;song_binding_file_ = song_kernel_binding_file_sha256;
    Json orbitals = Json::array();
    for (std::size_t i = 0; i < 4; ++i)
      orbitals.push_back({{"label", orbitals_[i].label}, {"B_eV", orbitals_[i].B_eV}, {"U_eV", orbitals_[i].U_eV},
                          {"N", orbitals_[i].N}, {"Q", orbitals_[i].Q}, {"parent_weight", weights_[i]}});
    description_ = {
      {"schema", "parent-coverage-v2-composite-dataset"},
      {"song_band", {{"dataset_sha256", song_sha_}, {"channel", w::parent_channel}, {"support_eV", {13.5, 1000.}},
                     {"kernel", w::parent_model}, {"basis", "dataset"}}},
      {"outer3_band", {{"dataset_sha256", outer3_sha_}, {"dataset_id", d.at("id")}, {"channel", outer3_channel},
                       {"rate_model", outer3_model}, {"band_eV", {{"lower_exclusive", beb_lo_}, {"upper_inclusive", beb_hi_}}},
                       {"orbitals", orbitals}, {"bohr_radius_m", k_.bohr_radius_m}, {"rydberg_eV", k_.rydberg_eV},
                       {"rate", "dataset arithmetic: sum_i weight_i S_i/(t+u+1)[Q ln t/2 (1-1/t^2)+(2-Q)(1-1/t-ln t/(t+1))]"},
                       {"kernel", sdcs_kernel}, {"basis", "dataset"},
                       {"source", d.at("sources").at("murakami2012")}, {"orbital_source", d.at("sources").at("nist_beb")}}},
      {"provisional_band", {{"admitted", provisional_.admit}, {"model", provisional_model},
                            {"band_eV", {{"lower_exclusive", beb_hi_}, {"upper_inclusive", provisional_.maximum_eV}}},
                            {"orbitals", orbitals}, {"bohr_radius_m", k_.bohr_radius_m}, {"rydberg_eV", k_.rydberg_eV},
                            {"electron_rest_energy_eV", k_.rest_energy_eV},
                            {"formula", "Kim, Santos & Parente, Phys. Rev. A 62, 052710 (2000), RBEB with alpha^4 = (2R/mc^2)^2, "
                                        "times per-orbital parent weights"},
                            {"kernel", sdcs_kernel}, {"basis", "model_derived"},
                            {"validation_targets_not_in_dataset",
                             {"Schutten et al., J. Chem. Phys. 44, 3924 (1966), total ionization 0.1-20 keV",
                              "Rieke & Prepejchal, Phys. Rev. A 6, 1507 (1972), total ionization 0.1-2.7 MeV"}}}}};
    arena_ = {sha(description_.dump()), kernel_manifest_, kernel_binding_file_, banded_model};
  }
  static HighEnergyBands read(const std::string& outer3_path, const std::string& expected_outer3_sha256,
                              const w::Data& song_data, const std::string& song_dataset_sha256,
                              const std::string& kernel_binding_path, const std::string& song_kernel_sha256,
                              const std::string& song_kernel_binding_file_sha256, ProvisionalOptions provisional = {}) {
    return HighEnergyBands(read_bytes(outer3_path), expected_outer3_sha256, song_data, song_dataset_sha256,
                           read_bytes(kernel_binding_path), song_kernel_sha256, song_kernel_binding_file_sha256, provisional);
  }
  double beb_low_eV() const noexcept { return beb_lo_; }
  double beb_high_eV() const noexcept { return beb_hi_; }
  const ProvisionalOptions& provisional() const noexcept { return provisional_; }
  const std::vector<beb_v2::Orbital>& orbitals() const noexcept { return orbitals_; }
  const std::vector<double>& weights() const noexcept { return weights_; }
  const beb_v2::Constants& constants() const noexcept { return k_; }
  const Json& description() const noexcept { return description_; }
  const ArenaBinding& arena_binding() const noexcept { return arena_; }
  const std::string& outer3_sha256() const noexcept { return outer3_sha_; }
  const std::string& song_dataset_sha256() const noexcept { return song_sha_; }
  const std::string& song_kernel_sha256() const noexcept { return song_kernel_; }
  const std::string& song_kernel_binding_file_sha256() const noexcept { return song_binding_file_; }

  // Dataset-backed outer3 parent rate on its closed support [1000, 10000] eV.
  double outer3_parent_m2(double E) const {
    pc::require(std::isfinite(E) && E >= beb_lo_ && E <= beb_hi_, "outer3 parent rate outside 1000..10000 eV");
    return beb_v2::beb_weighted_sum_m2(orbitals_, weights_, k_, E);
  }
  // Four-orbital unweighted dataset BEB (H2O:nist_BEB_valence arithmetic).
  double valence_beb_m2(double E) const { return beb_v2::beb_weighted_sum_m2(orbitals_, {}, k_, E); }
  // Per-orbital weighted terms of the outer3 rate, in dataset arithmetic.
  std::vector<long double> outer3_terms(double E) const {
    std::vector<long double> out;
    for (std::size_t i = 0; i < 4; ++i) {
      if (weights_[i] == 0) { out.push_back(0);continue; }
      out.push_back(beb_v2::beb_weighted_sum_m2({orbitals_[i]}, {weights_[i]}, k_, E));
    }
    return out;
  }
  std::vector<long double> provisional_terms(double E) const {
    std::vector<long double> out;
    for (std::size_t i = 0; i < 4; ++i)
      out.push_back(weights_[i] == 0 ? 0.L : weights_[i] * beb_v2::rbeb_orbital_m2(orbitals_[i], k_, E));
    return out;
  }
  // Provisional parent rate (RBEB x parent weights), (10000, maximum] eV.
  double provisional_parent_m2(double E) const {
    pc::require(provisional_.admit && std::isfinite(E) && E > beb_hi_ && E <= provisional_.maximum_eV,
                "provisional parent rate outside its admitted band");
    long double sum = 0;for (auto x : provisional_terms(E)) sum += x;
    return static_cast<double>(sum);
  }
  // Diagnostics without band restriction (reports and tests only).
  long double rbeb_parent_any(long double E) const {
    long double s = 0;for (std::size_t i = 0; i < 4; ++i) if (weights_[i] > 0) s += weights_[i] * beb_v2::rbeb_orbital_m2(orbitals_[i], k_, E);return s;
  }
  long double rbeb_valence_any(long double E) const {
    long double s = 0;for (const auto& o : orbitals_) s += beb_v2::rbeb_orbital_m2(o, k_, E);return s;
  }
  long double beb_parent_any(long double E) const {
    long double s = 0;for (std::size_t i = 0; i < 4; ++i) if (weights_[i] > 0) s += weights_[i] * beb_v2::beb_orbital_ext_m2(orbitals_[i], k_, E);return s;
  }
  long double beb_valence_any(long double E) const {
    long double s = 0;for (const auto& o : orbitals_) s += beb_v2::beb_orbital_ext_m2(o, k_, E);return s;
  }
};
} // namespace raikou::native_water::coverage_v2
