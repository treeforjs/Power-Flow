#pragma once
// Hall v3.4 device coupled step (namespace hall_corrective_v3::gpu): the per-carrier and per-mesh host work of
// the resident owner's coupled proposal moved to the device, on the owner's session and stream:
//   census      host census rows gathered once per proposal (80-byte records; or rows already on the device),
//               identity checks, the census charge of the residents (charge_of_particles: Words CIC in carrier
//               order, x87 gross) compared with the accepted free charge (same_net_charge, x87), census v2
//               (pack, check, block-list root: native v7 forms), newborn flights and start rule;
//   motion      the v3 push (MotionEngineV3::run_device) on the device trial midpoint, the motion fold on the
//               device (exact K, W, G_K, C_K, Q, G_W; continuity; E^; free-charge guard; record operand checks),
//               sinks and births compacted on the device and downloaded (O(sinks + newborns));
//   support     the deposit support of the endpoint continuity guard (support_v34.hpp, x87 emulated, bit for
//               bit the sealed DepositSupportBuilder), written where the field step reads it;
//   record      the motion record v7 section arrays (canonical layouts) and their block-list roots on the device
//               for the accepted iterate only; host materialisation of endpoints/paths on request;
//   blocks      block-form digests of binary64 arrays (native v7.1 trial field and current).
// Refusals are resolved on the host from O(1) status words in the order of the v3.3 host path, with the same
// texts. The header exposes no CUDA types.
#include "resident_field.hpp"
#include "pf_v7/forms.hpp"
#include <functional>

namespace pf_v7::gpu { class Forms; }

namespace hall_corrective_v3::gpu {

// Census of one proposal.
struct CensusSource {
    std::span<const native::ArenaRateEntry> rows;               // host census (native v7 trial.native_census)
    std::span<const std::uint64_t> packet_keys, state_indices;  // native census v2 tables in ascending native ID (trial.census_v7());
                                                                // empty: computed from the packet/state texts (v3.3 rule)
    std::span<const std::string> states;                        // census v2 state table (empty: sorted states of the rows)
    // Rows already on the device (native v7.1 census_v7_device): 80-byte records ascending, keys, state indices.
    std::uint64_t device_rows = 0, device_keys = 0, device_states = 0;
    std::size_t device_count = 0;
    // Packet and state texts of one carrier (device rows; host rows carry them); O(sinks + newborns) calls.
    std::function<std::pair<std::string, std::string>(std::uint64_t)> packet_of;
};

// Host/device times of one proposal (ms), by item.
struct DeviceStepTimes {
    double census_gather = 0, census_device = 0, motion_host = 0, motion_wait = 0, support = 0, record = 0, field_digest = 0;
    double fold_host = 0, sinks_host = 0, materialise = 0;
    double device_push_ms = 0, device_fold_ms = 0;
    long iterations = 0;
    // thread CPU time (ms) of the host thread in each item (wall times above include device waits)
    double census_cpu = 0, motion_cpu = 0, record_cpu = 0, preflight_cpu = 0;
};

class DeviceStep {
    struct Impl;
    std::unique_ptr<Impl> p_;
public:
    DeviceStep(std::shared_ptr<Device>, const Model&, const MotionContext&, std::size_t capacity);
    ~DeviceStep();
    DeviceStep(const DeviceStep&) = delete;
    DeviceStep& operator=(const DeviceStep&) = delete;

    pf_v7::gpu::Forms& forms();                // the v7 forms instance on the session (native v7.1 device lane)
    std::string compile_json() const;
    std::size_t capacity() const noexcept;

    // ---- census admission (first statement of the coupled proposal). Throws the v3.3 texts in the v3.3 order:
    // resident identity (cap, order, IDs, texts), charge_of_particles of the residents, free_charge_guard,
    // same_net_charge with the accepted free charge (accepted slot words on the device), then the accepted
    // field admitted with the census charge (resident field, device).
    void begin(const CensusSource&, std::span<const Newborn> newborns, double begin_s, double end_s, double dt, ResidentField& field,
               const SlotMeta& accepted);
    // Census v2 digest of the proposal (device root + host header); refusal deferred to the first motion.
    const CensusV2& census_v2() const noexcept;
    std::size_t carriers() const noexcept;

    // ---- one fixed-point iterate: motion on the device trial midpoint (E, B device addresses).
    struct Iterate {
        hall_driver::MotionReceipt receipt;    // host fields O(1) + sinks/births; arrays on the device
        std::uint64_t current = 0, after_hi = 0, after_lo = 0, after_gross = 0, support = 0;   // device inputs of the field step
        double host_ms = 0, wait_ms = 0;
    };
    Iterate motion(const hall_driver::MotionBinding&, std::uint64_t E, std::uint64_t B, const VerifiedField& trial, bool birth_scope);
    // Deposit-support refusal of the last iterate (in the v3.3 position: after the bank transfer).
    void support_check() const;

    // ---- accepted iterate: record arrays and section roots on the device; seals the receipt (record v7 digest).
    // host_arrays: download endpoints and paths into the receipt (native v7 host admission).
    // form: the trial-field and current digests of the record (blocks_v71: computed on the device from the block-list
    // roots of trial_E, trial_B and the current; the binding's trial_field_sha256 is set here).
    void seal(hall_driver::MotionReceipt&, bool host_arrays, std::uint64_t trial_E = 0, std::uint64_t trial_B = 0,
              FieldDigestForm form = FieldDigestForm::streaming_v7);
    // Preflight recomputation: census v2 root and the four section roots recomputed on the device from the device
    // census and record arrays, the census charge of the residents recomputed and compared with the accepted free
    // charge, the endpoint free charge recomputed from the endpoint records and compared with the candidate slot.
    bool preflight(const hall_driver::MotionReceipt&, std::span<const native::ArenaRateEntry> census, ResidentField& field, const SlotMeta& accepted);

    // ---- block-form digests (native v7.1, API-071.txt 1.1): roots of binary64 device arrays, which = 0 field-E,
    // 1 field-B, 2 current-J (hex of the root digest).
    std::string block_root(int which, std::uint64_t values, std::size_t n);
    // The same root by the native v7.1 forms library (Forms::block_lists; qualification of the two device paths).
    std::string block_root_forms(int which, std::uint64_t values, std::size_t n);

    // ---- native v7.1 device lane (HALL-LINK-7.1.txt items 1, 3-5): the record arrays of the sealed receipt and the
    // A1 inputs on the shared session (canonical layouts, ascending native ID), valid until the next motion() of
    // this engine; the block-form digests of the trial field and the current, computed on the device at the seal:
    //   trial_field_blocks = SHA256(text "total-midpoint-E-B-blocks-v1", text E root, text B root)
    //   current_blocks     = SHA256(text "motion-array-blocks-v1", text J root)
    struct DeviceRecordView {
        std::uint64_t endpoints = 0, paths = 0, sinks = 0, births = 0;   // 80, 64, 136, 24-byte records
        std::size_t n_endpoints = 0, n_paths = 0, n_sinks = 0, n_births = 0;
        std::uint64_t E = 0, B = 0, J = 0;                                // trial midpoint E [n_E], B [n_B], current J [n_E]
        std::size_t n_E = 0, n_B = 0;
        std::string E_root, B_root, J_root, trial_field_blocks_sha256, current_blocks_sha256;
        std::uint64_t seal_generation = 0;
    };
    const DeviceRecordView& record_view() const noexcept;
    hall_gpu::Session& session();

    const DeviceStepTimes& times() const noexcept;
    void reset_times() noexcept;
    // Last iterate's device support (download; qualification).
    Vec support_download();
    // Device record arrays of the sealed receipt (download; qualification and tamper tests).
    std::vector<std::uint64_t> record_words(int section);
    void tamper_record_word(int section, std::size_t word, std::uint64_t xor_mask);   // qualification only (4: device census rows, 5: current J)
};

// Bank transfer of one step restricted to the eight CIC nodes of every deposit (census_charge.hpp bank_transfer:
// same Words/x87 operations per node in the same deposit order; O(sinks + newborns)).
struct SparseBank {
    bool empty = true;
    std::vector<int> nodes;                         // ascending
    WordVec conductor_C, exterior_C; Vec gross_C, residue_C;
    bool residue_any = false;
};
SparseBank bank_transfer_sparse(const Mesh&, const native::GapPicOptions&, const hall_driver::MotionReceipt&, const std::vector<bool>& exterior_support);
// Full-size form (comparison with bank_transfer).
BankTransfer dense_of(const Mesh&, const SparseBank&);

}  // namespace hall_corrective_v3::gpu
