// chem_v3_native.hpp -- native v8 chemistry v3 owner core (API-008.txt 3).
//
// Links the sealed chemistry v3 component (pic-water-chemistry-v3-20261031,
// headers read in place): registry of the 22 blob species, the collision grid,
// stage inputs from exact census sums, the host stage (C3 mcc_carrier), the
// canonical record digests, debit resolution (C3 resolve_debits), the
// application plan, exact ledgers (C3 Ledger plus nuclei ABI v4 and gas
// heating), the coverage-ledger-v3 policy and the quantized neutral release.
// No arena types besides GapParticle and the public configuration.
#pragma once
#include "chem_host.hpp"
#include "water_network_scope.hpp"
#include "pf_v7/forms.hpp"
#include "raikou/mfem_pic_full_time.h"
#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace raikou::mfem_backend {
// Accepted chemistry state (immutable once shared; API-008.txt 3.6).
struct FullTimeChemistryState;
namespace chem8 {
inline constexpr char owner_id[] = "native-chemistry-v3";
inline constexpr char channel_id[] = "pfc-chem-v3";
inline constexpr char model_id[] = "pfc-chem-v3-null-collision-mcc";
inline constexpr char state_prefix[] = "pfc-v3/";
inline constexpr char packet_namespace[] = "chem-v3/";
inline constexpr std::uint32_t none = 0xFFFFFFFFu;

using ESum = pfc::ExactSum;
// Nuclei ABI v4 (the shared 16-slot layout of Raikou's kinetic core,
// raikou-mfem-pic 934fb24 chem_tables.h: PFC_SPECIES_LAYOUT 4, PFC_NUCLEUS_*):
// same slot order and values; chemistry v3 records (6 slots) map slot k -> k.
inline constexpr std::uint32_t species_layout = 4, nuclei_slots = 16;
inline constexpr std::uint32_t nucleus_H = 0, nucleus_O = 1, nucleus_Fe = 2, nucleus_Cr = 3, nucleus_Ni = 4,
                               nucleus_Mn = 5, nucleus_Si = 6, nucleus_C = 7, nucleus_N = 8; // 9-15 spare
// ------------------------------------------------------------- registry
struct Registry {
  std::shared_ptr<const pfc::Blob> blob;
  std::vector<std::string> states;       // "pfc-v3/<name>", blob index order
  std::vector<double> mass_kg;
  std::vector<int> charge;
  std::vector<std::array<std::int32_t, 16>> nuclei; // ABI v4
  std::vector<std::uint8_t> projectile, target;
  std::vector<std::uint8_t> region_bounded;
  std::string sha256;
};
std::shared_ptr<const Registry> load_registry(const std::string &path, const std::string &expected_sha256);
// ------------------------------------------------------------- grid
// native-v8/collision-grid-v1: uniform Cartesian cells over [0, width] x
// [0, gap] x [short_z, feed_z]; foil cells are solid. Binary64 locate:
// i_d = min(n_d - 1, (u64)(((x_d - low_d) / (high_d - low_d)) * n_d)).
struct Grid {
  std::array<double, 3> low{}, high{};
  std::array<std::uint32_t, 3> n{};
  double volume = 0;
  std::vector<std::uint8_t> vacuum;
  std::string sha256;
  std::uint64_t cells() const noexcept { return std::uint64_t(n[0]) * n[1] * n[2]; }
  std::uint32_t locate(const double *x) const noexcept;
};
Grid make_grid(const GapPicOptions &, const std::array<std::uint32_t, 3> &counts);
// ------------------------------------------------------------- setup
struct Setup {
  std::string base_configuration_sha256;
  std::shared_ptr<const water_scope::Guard> water_scope;
  FullTimeChemistryConfiguration config;
  std::shared_ptr<const Registry> reg;
  Grid grid;
  std::vector<std::uint8_t> admitted; // per region
  std::string policy_sha256, configuration_sha256;
  ArenaChannel channel;
};
std::shared_ptr<const Setup> make_setup(const FullTimeChemistryConfiguration &, const GapPicOptions &);
// ------------------------------------------------------------- ledgers
struct Key {
  std::uint32_t cell = 0, species = 0; // atlas cell, species
  bool operator<(const Key &o) const noexcept { return cell != o.cell ? cell < o.cell : species < o.species; }
  bool operator==(const Key &o) const noexcept { return cell == o.cell && species == o.species; }
};
struct Ledger {
  pfc::Ledger c3;                                   // C3 energy/count/process/coverage
  std::array<ESum, 16> nuclei_created, nuclei_destroyed; // ABI v4
  std::map<Key, ESum> gas_heating;                  // eV, per (atlas cell, target species)
  std::uint64_t admitted_events = 0, refused_debits = 0, births = 0, consumed = 0, released_packets = 0,
                debited_packets = 0, removed_packets = 0, remainder_carriers = 0, stages = 0;
  ESum released_weight, birth_weight, consumed_weight, debited_weight;
  void init(const pfc::Blob &b) { c3.init(b); }
  void merge(const Ledger &o);
};
std::string ledger_digest(const Ledger &);
// ------------------------------------------------------------- reservoir
using Reservoir = pfc::Reservoir; // pfc-neutral-parcels-v1, persistent conservative phase-space stock
} // namespace chem8
// v8 ejected-melt reservoir (exact; immutable once shared).
struct FullTimeMeltState {
  std::uint64_t credits = 0;
  std::array<pfc::ExactSum, 6> atoms;
  pfc::ExactSum energy_J, kinetic_J;
  std::string last_record_sha256;
  mutable std::string digest;
};
// v8.1 flows state (API-008.txt 8.5; exact; immutable once shared; null until
// the first v8.1 event): internal energy carried by packets (metal vapour
// births), cumulative wall-response / born-and-collected / bridge births,
// cumulative extra motion work, merge/split cumulative counters and exact
// residuals, and the internal energies of the last step's sunk macros.
struct FullTimeFlowsV81 {
  struct Internal {
    double weight = 0, energy_J = 0; // packet weight at birth, carried energy
  };
  std::map<std::string, Internal> internal; // packet ID -> carried internal energy
  pfc::ExactSum internal_born_J, internal_sunk_J, internal_released_J;
  std::uint64_t wall_births = 0, returned_births = 0, bridges = 0;
  std::map<std::string, pfc::ExactSum> wall_weight, returned_weight; // per state
  pfc::ExactSum extra_work_J; // exact sum of the booked per-step extra works
  std::uint64_t merge_transactions = 0, merge_groups = 0, merge_inputs = 0, merge_outputs = 0,
                merge_relocations = 0, merge_formation_refused = 0, merge_centroid_fallback = 0;
  std::array<std::uint64_t, 4> merge_status{};
  // per state: momentum x, y, z (kg m/s) and energy (J) residuals (exact)
  std::map<std::string, std::array<pfc::ExactSum, 4>> merge_residual;
  bool has_last_merge = false;
  std::uint64_t last_merge_step = 0, last_merge_base = 0, last_merge_outputs = 0;
  std::string last_merge_sha256, merge_policy_sha256;
  std::uint64_t last_sink_step = 0; // accepted step + 1 of last_sink_internal (0: none)
  std::vector<std::pair<std::uint64_t, double>> last_sink_internal; // sunk native ID -> J
  mutable std::string digest;
};
struct FullTimeChemistryState {
  std::string configuration_sha256; // chemistry configuration identity
  // Latest admitted census and rolling scope history, authenticated by state/row roots.
  std::string water_scope_census_sha256, water_scope_binding_sha256, water_scope_history_sha256;
  std::string host_target_feedback_sha256, host_target_history_sha256;
  std::map<chem8::Key, chem8::Reservoir> reservoirs;
  chem8::Ledger cumulative;     // every accepted stage
  chem8::Ledger step;           // pending: stage 1 of the step in progress (empty between steps)
  chem8::Ledger last_step;      // last sealed step
  bool pending = false;         // stage 1 done, stage 2 not yet
  std::uint64_t stages = 0;
  mutable std::string digest;   // cache (immutable once shared)
};
namespace chem8 {
std::string state_digest(const FullTimeChemistryState &);
std::string melt_digest(const FullTimeMeltState &);
std::string encode_melt(const FullTimeMeltState &);
std::shared_ptr<FullTimeMeltState> decode_melt(const std::string &);
std::string flows_v81_digest(const FullTimeFlowsV81 &);
std::string encode_flows_v81(const FullTimeFlowsV81 &);
std::shared_ptr<FullTimeFlowsV81> decode_flows_v81(const std::string &);
inline std::string state_digest(const Setup &, const FullTimeChemistryState &s) { return state_digest(s); }
// Canonical binary codec of the state (checkpoint) and its inverse.
std::string encode_state(const FullTimeChemistryState &);
std::shared_ptr<FullTimeChemistryState> decode_state(const Setup &, const std::string &);
std::shared_ptr<FullTimeChemistryState> initial_state(const Setup &);
// ------------------------------------------------------------- stage
struct RowView {
  const GapParticle *rows = nullptr;
  std::size_t n = 0;
  const std::uint32_t *species = nullptr; // per row: registry index or none
  const std::uint8_t *single = nullptr;   // per row: 1 if its packet has one carrier
};
struct StageRequest {
  std::uint64_t step = 0;
  unsigned stage = 1;     // 1: before transport (0x30), 2: after (0x34)
  double dt_s = 0;        // authoritative dt / 2
  std::size_t max_release_packets = 1000000; // bounded by outer carrier cap
};
struct StageInputs {
  pfc::StepConfig cfg{};
  std::vector<pfc::CarrierIn> carriers;      // projectiles, ascending native ID
  std::vector<std::uint32_t> carrier_row;    // row index of each projectile
  std::vector<pfc::CellTargets> cells;       // compact cells (ascending atlas ID)
  std::vector<std::uint32_t> cell_atlas;     // compact -> atlas cell
  std::vector<std::uint32_t> row_cell;       // per row: atlas cell (none: no chemistry species)
  std::map<pfc::InventoryKey, ESum> inventory; // (compact cell, target species) exact counts
  std::string sha256;                        // inputs digest
};
void rehash_inputs(StageInputs &, const std::string &target_feedback_sha256 = {});
StageInputs build_inputs(const Setup &, const RowView &, const StageRequest &);
// native-v8 block roots (API-008.txt 3.3): leaf = SHA256(text tag/leaf, u64
// leaf, u64 count, record bytes); root = SHA256(text tag/root, u64 n, u64
// leaves, leaf digests). Inputs digest (native-v8/chemistry-inputs-v2) from
// the carrier (1024/leaf), cell (64/leaf) and atlas (u64, 4096/leaf) roots.
std::string block_root(const std::string &tag, const void *data, std::size_t record_bytes, std::uint64_t n,
                       std::uint64_t per_leaf);
std::string block_root_of(const std::string &tag, std::uint64_t n, const std::vector<std::array<unsigned char, 32>> &);
std::string inputs_digest(const pfc::StepConfig &, const std::string &carriers_root, std::uint64_t n_carriers,
                          const std::string &cells_root, std::uint64_t n_cells, const std::string &atlas_root);
// Records of a stage: CarrierOut per projectile and the populated events in
// canonical order (carrier-major, candidate order); event_begin[i] is the
// first event of carrier i (size carriers + 1).
struct StageRecords {
  std::vector<pfc::CarrierOut> outs;
  std::vector<pfc::EventRec> events;
  std::vector<std::uint64_t> event_begin;
};
StageRecords run_host(const Setup &, const StageInputs &, unsigned partitions);
// pfc-chem-v3-step-records-v1 (C3 step_digest) over the compact records.
std::string records_stream_sha256(const StageRecords &);
// native-v8/chemistry-records-blocks-v1: block Merkle over 256-carrier leaves.
std::string records_root(const StageRecords &);
// ------------------------------------------------------------- application plan
struct Birth {
  std::uint32_t row = 0;      // projectile row (position)
  std::uint64_t carrier = 0;  // projectile native ID
  std::uint32_t candidate = 0, product = 0, species = 0, process = 0;
  double weight = 0;
  double u[3] = {0, 0, 0};
};
struct Release {
  std::uint32_t cell = 0, species = 0, index = 0;
  double weight = 0;
  double u[3] = {0, 0, 0}, x[3] = {0, 0, 0};
};
struct TargetDebit {
  std::uint32_t row = 0;  // first row of the debited packet
  ESum amount;            // exact debit (<= the packet count)
  bool whole = false;     // the packet is removed
};
struct Plan {
  std::vector<std::pair<std::uint32_t, std::array<double, 3>>> velocity; // surviving projectiles with events
  std::vector<std::uint32_t> consumed_rows;
  std::vector<TargetDebit> debits;  // ascending row
  std::vector<Birth> births;        // canonical order
  std::vector<Release> releases;    // canonical order (stage 2)
};
// Resolution of one stage: statuses, debits (prefix rule), ledgers, coverage
// policy, plan, reservoir deposits and (stage 2) the release. `state` is the
// chemistry state of the frame before the stage; the result is the state
// after it. Throws FullTimeChemistryRefusal.
struct Resolution {
  Plan plan;
  std::shared_ptr<FullTimeChemistryState> state;
  Ledger stage;
  pfc::DebitResult debits;
  bool policy_refused = false;
  std::string policy_message;
  std::vector<FullTimeChemistryAdmission> admissions;
  std::vector<FullTimeChemistryAdmittedEvent> admitted_records;
};
// Carriers with events (ascending native ID) and the totals over every
// carrier record of the stage (host: from the records; device: reductions).
struct ActiveCarrier {
  std::uint32_t row = 0, cell = 0, species = 0, status = 0; // cell: compact index
  std::uint64_t id = 0, begin = 0, end = 0;                 // event range
  double x[3] = {0, 0, 0};
  double weight = 0;
};
struct CarrierTotals {
  std::uint64_t carriers = 0, candidates = 0;
  std::map<std::uint32_t, std::uint64_t> status_counts;
  double nu_dt_max = 0;
  std::vector<ESum> exposure, bound;  // per region (exact)
  std::vector<std::uint64_t> support; // reachable carrier/region pairs, including zero rate
  std::uint64_t bad_carrier = ~0ull;  // first carrier (ascending ID) with a refusing status
  std::uint32_t bad_status = 0;
};
struct StageActive {
  std::vector<ActiveCarrier> carriers;
  std::vector<pfc::EventRec> events;
  CarrierTotals totals;
  std::vector<std::uint32_t> cell_atlas;
};
StageActive active_of(const Setup &, const StageInputs &, const StageRecords &, const RowView &);
// inventory(compact cell, species): exact stage-start count of the target
// species' carriers in the cell; packets(atlas cell, species): (first row,
// exact count) of the species' packets in the cell, ascending native ID.
struct DebitAccess {
  std::function<ESum(std::uint32_t, std::uint32_t)> inventory;
  std::function<std::vector<std::pair<std::uint32_t, ESum>>(std::uint32_t, std::uint32_t)> packets;
};
Resolution resolve(const Setup &, const StageActive &, const StageRequest &, const FullTimeChemistryState &before,
                   const DebitAccess &);
// ------------------------------------------------------------- helpers
std::string sha_hex(const std::string &bytes);
std::string exact_digits(const ESum &);          // forms_core digits (hex text)
FullTimeExactValue exact_value(const ESum &);
pf_v7::ExactSum to_forms(const ESum &);
std::string packet_id_birth(std::uint64_t step, unsigned stage, const Birth &);
std::string packet_id_release(std::uint64_t step, const Release &);
bool reserved_packet_id(const std::string &) noexcept;
// Canonical chemistry packet ID (birth or release of the given step).
bool canonical_packet_id(const std::string &id) noexcept;
std::string status_name(std::uint32_t);
void fill_receipt(FullTimeChemistryStageReceipt &, const Setup &, const StageInputs &, const StageActive &,
                  const Resolution &, const StageRequest &);
std::string coverage_json(const Setup &, const FullTimeChemistryState &);
std::vector<std::pair<std::string, FullTimeExactValue>> energy_values(const Ledger &);
std::vector<std::pair<std::string, FullTimeExactValue>> count_values(const Setup &, const Ledger &);
} // namespace chem8
} // namespace raikou::mfem_backend
